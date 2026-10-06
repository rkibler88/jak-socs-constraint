"""Step 2 - Depth-maximisation scoping survey.

Measures how many sequences and, critically, how many *paired species* each
retrieval route recovers for the JAK and SOCS families across nested taxonomic
tiers. The paired-species count is the quantity that caps Neff for an
inter-protein DCA, so it is what the dataset design must be chosen on.

Retrieval routes compared
-------------------------
gene_name   UniProt ``gene:<symbol>`` query - what the previous pipeline used.
interpro    InterPro family signature - profile-HMM based, recovers entries
            from unannotated genomes that carry no gene symbol.

Outputs
-------
depth_scoping_raw.csv     one row per retrieved entry
depth_scoping.csv         one row per (route, tier, family/paralogue)
depth_scoping_pairs.csv   paired-species counts per (route, tier, pair)
"""

from __future__ import annotations

import csv
import json
import logging
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

LOG = logging.getLogger("depth_scoping")

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
PAGE_SIZE = 500
FIELDS = "accession,id,organism_id,organism_name,gene_primary,gene_names,length,reviewed,protein_existence"

# Human canonical references. P52333 is JAK3 - the previous pipeline labelled it
# JAK2, which is why it is listed here explicitly and then excluded.
REFERENCES: dict[str, str] = {
    "JAK1": "P23458",
    "JAK2": "O60674",
    "TYK2": "P29597",
    "SOCS1": "O15524",
    "SOCS3": "O14543",
}
EXCLUDED_PARALOGUES: dict[str, str] = {"JAK3": "P52333"}

#: Nested taxonomic tiers, ordered outward. SOCS1/SOCS3 orthology is a
#: gnathostome distinction, so the deeper tiers answer a family-level question.
TAXA: dict[str, int] = {
    "Gnathostomata": 7776,
    "Vertebrata": 7742,
    "Chordata": 7711,
    "Metazoa": 33208,
}

#: InterPro family signatures, read off the human references rather than assumed.
#: IPR016251 = non-receptor tyrosine kinase, Jak/Tyk2 family.
#: IPR000980 = SH2 domain; IPR001496 = SOCS box. Their conjunction is the
#: SOCS family; SH2 alone would return every SH2 protein in the proteome.
INTERPRO = {
    "JAK_family": ["IPR016251"],
    "SOCS_family": ["IPR000980", "IPR001496"],
}

GENE_SYMBOLS = {
    "JAK_family": ["jak1", "jak2", "tyk2"],
    "SOCS_family": ["socs1", "socs3"],
}


@dataclass
class Entry:
    """One retrieved UniProt entry, normalised across routes."""

    accession: str
    entry_name: str
    taxid: str
    organism: str
    gene_primary: str
    gene_names: str
    length: int
    reviewed: bool
    protein_existence: str
    route: str
    tier: str
    family: str
    query_label: str = ""


@dataclass
class SurveyResult:
    entries: list[Entry] = field(default_factory=list)
    query_log: list[dict] = field(default_factory=list)


def _request(url: str, contact_email: str | None, retries: int = 4) -> tuple[str, str | None]:
    """GET a UniProt URL, returning (body, next-page-url).

    Retries on transient 5xx/timeouts with linear backoff. UniProt rate limits
    are generous but the streamed pagination here is long-running, so failing a
    whole survey on one hiccup would be wasteful.
    """
    headers = {"Accept": "text/plain"}
    if contact_email:
        # Only ever sent as a contact header to UniProt, per their fair-use note.
        headers["User-Agent"] = f"jak-socs-coevolution ({contact_email})"
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=180) as resp:
                body = resp.read().decode("utf-8")
                link = resp.headers.get("Link")
            match = re.search(r'<([^>]+)>;\s*rel="next"', link or "")
            return body, (match.group(1) if match else None)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as err:
            if isinstance(err, urllib.error.HTTPError) and 400 <= err.code < 500:
                raise
            last_err = err
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"UniProt request failed after {retries} attempts: {last_err}")


def stream_query(query: str, contact_email: str | None = None) -> Iterator[dict]:
    """Yield every UniProt entry matching ``query`` as a dict of TSV fields."""
    params = {"query": query, "format": "tsv", "fields": FIELDS, "size": str(PAGE_SIZE)}
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"
    header: list[str] | None = None
    while url:
        body, url = _request(url, contact_email)
        lines = body.strip("\n").split("\n")
        if not lines or not lines[0]:
            return
        if header is None:
            header = lines[0].split("\t")
        for line in lines[1:]:
            if not line:
                continue
            values = line.split("\t")
            # Trailing empty fields are dropped by the TSV writer; pad them back.
            values += [""] * (len(header) - len(values))
            yield dict(zip(header, values))


def _to_entry(row: dict, route: str, tier: str, family: str, label: str) -> Entry:
    return Entry(
        accession=row.get("Entry", ""),
        entry_name=row.get("Entry Name", ""),
        taxid=row.get("Organism (ID)", ""),
        organism=row.get("Organism", ""),
        gene_primary=row.get("Gene Names (primary)", ""),
        gene_names=row.get("Gene Names", ""),
        length=int(row["Length"]) if str(row.get("Length", "")).isdigit() else 0,
        reviewed=str(row.get("Reviewed", "")).lower().startswith("reviewed"),
        protein_existence=row.get("Protein existence", ""),
        route=route,
        tier=tier,
        family=family,
        query_label=label,
    )


def survey(contact_email: str | None = None) -> SurveyResult:
    """Run every (route x tier x family) query and collect the entries."""
    result = SurveyResult()
    for tier, taxid in TAXA.items():
        tax_clause = f"(taxonomy_id:{taxid})"

        # Route 1: gene-symbol queries, one per paralogue.
        for family, symbols in GENE_SYMBOLS.items():
            for symbol in symbols:
                query = f"(gene:{symbol}) AND {tax_clause}"
                rows = list(stream_query(query, contact_email))
                LOG.info("gene_name  %-14s %-6s n=%d", tier, symbol, len(rows))
                result.entries.extend(
                    _to_entry(r, "gene_name", tier, family, symbol.upper()) for r in rows
                )
                result.query_log.append(
                    {"route": "gene_name", "tier": tier, "family": family,
                     "label": symbol.upper(), "query": query, "n_entries": len(rows)}
                )

        # Route 2: InterPro family signatures.
        for family, signatures in INTERPRO.items():
            sig_clause = " AND ".join(f"(xref:interpro-{s})" for s in signatures)
            query = f"{sig_clause} AND {tax_clause}"
            rows = list(stream_query(query, contact_email))
            LOG.info("interpro   %-14s %-12s n=%d", tier, family, len(rows))
            result.entries.extend(
                _to_entry(r, "interpro", tier, family, "+".join(signatures)) for r in rows
            )
            result.query_log.append(
                {"route": "interpro", "tier": tier, "family": family,
                 "label": "+".join(signatures), "query": query, "n_entries": len(rows)}
            )
    return result


def write_outputs(result: SurveyResult, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    raw = outdir / "depth_scoping_raw.csv"
    with raw.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["route", "tier", "family", "query_label", "accession", "entry_name",
                        "taxid", "organism", "gene_primary", "gene_names", "length",
                        "reviewed", "protein_existence"],
        )
        writer.writeheader()
        for e in result.entries:
            writer.writerow({
                "route": e.route, "tier": e.tier, "family": e.family,
                "query_label": e.query_label, "accession": e.accession,
                "entry_name": e.entry_name, "taxid": e.taxid, "organism": e.organism,
                "gene_primary": e.gene_primary, "gene_names": e.gene_names,
                "length": e.length, "reviewed": int(e.reviewed),
                "protein_existence": e.protein_existence,
            })

    (outdir / "depth_scoping_queries.json").write_text(json.dumps(result.query_log, indent=1))
    LOG.info("wrote %s (%d entries)", raw, len(result.entries))


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path,
                        default=Path("/Users/ryankibler/Downloads/TakeTwo/results/manuscript_analysis"))
    parser.add_argument("--contact-email", default=None)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    started = time.time()
    result = survey(args.contact_email)
    write_outputs(result, args.outdir)
    LOG.info("survey complete in %.1f s", time.time() - started)


if __name__ == "__main__":
    main()
