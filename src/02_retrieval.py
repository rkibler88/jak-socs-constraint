"""Step 3 - Paralogue-resolved sequence retrieval.

Retrieves JAK-family and SOCS-family sequences by InterPro signature (the route
the step-2 scoping survey selected: +66% paired species over gene-symbol
queries) and assigns each to a specific paralogue by best-match identity
against human canonical references.

Gene symbols are deliberately not trusted: 87% of the JAK-family entries that
signature retrieval recovers and symbol queries miss carry no gene symbol at
all, and the previous pipeline's paralogue mixing came from trusting annotation
over sequence. Assignment here is by mmseqs2 search against a reference panel
that includes the *non-target* paralogues (JAK3, SOCS2, SOCS4-7, CIS), so an
off-target family member is rejected on evidence rather than missed.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib

LOG = logging.getLogger("step03")

#: Human canonical references. The panel deliberately includes paralogues we do
#: NOT want, so that a SOCS2 or JAK3 sequence is assigned to its own reference
#: and excluded, instead of being mis-assigned to a target paralogue.
TARGETS = {
    "JAK1": "P23458",
    "JAK2": "O60674",
    "TYK2": "P29597",
    "SOCS1": "O15524",
    "SOCS3": "O14543",
}
DECOYS = {
    "JAK3": "P52333",     # excluded by user requirement
    "SOCS2": "O14508",    # the non-binding contrast in Zhai et al. 2026
    "SOCS4": "Q8WXH5",
    "SOCS5": "O75159",
    "SOCS6": "O14544",
    "SOCS7": "O14512",
    "CISH": "Q9NSE2",
}
REFERENCE_PANEL = {**TARGETS, **DECOYS}

JAK_LIKE = {"JAK1", "JAK2", "JAK3", "TYK2"}

#: Assignment thresholds. Orthologues of a given JAK or SOCS paralogue are far
#: closer to their own human reference than to a sibling paralogue, so a
#: bitscore margin is the discriminating statistic; coverage guards against
#: fragments that would not carry the domains the analysis needs.
MIN_TARGET_COVERAGE = 0.55
MIN_IDENTITY = 0.35
MIN_BITSCORE_MARGIN = 1.10


@dataclass
class Assignment:
    """Paralogue call for one retrieved sequence, with its audit trail."""

    accession: str
    taxid: str
    organism: str
    claimed_gene: str
    length: int
    best_reference: str
    best_identity: float
    best_target_coverage: float
    best_bitscore: float
    runner_up_reference: str
    runner_up_bitscore: float
    bitscore_margin: float
    assigned_paralogue: str
    accepted: bool
    reason: str


def build_reference_fasta(path: Path, contact_email: str | None = None) -> dict[str, str]:
    seqs = dcalib.fetch_sequences(list(REFERENCE_PANEL.values()), batch=50,
                                  contact_email=contact_email)
    missing = [k for k, a in REFERENCE_PANEL.items() if a not in seqs]
    assert not missing, f"reference sequences unavailable: {missing}"
    records = {f"{name}|{acc}": seqs[acc] for name, acc in REFERENCE_PANEL.items()}
    dcalib.write_fasta(records, path)
    LOG.info("reference panel: %d proteins", len(records))
    return records


def run_mmseqs(query_fasta: Path, ref_fasta: Path, out_tsv: Path,
               tmp: Path, threads: int = 8) -> pd.DataFrame:
    """Search every retrieved sequence against the reference panel."""
    cols = ["query", "target", "pident", "alnlen", "qlen", "tlen",
            "qcov", "tcov", "bits", "evalue"]
    cmd = ["mmseqs", "easy-search", str(query_fasta), str(ref_fasta),
           str(out_tsv), str(tmp),
           "--format-output", ",".join(cols),
           "-s", "7.5", "--max-seqs", "50", "-e", "1e-5",
           "--threads", str(threads), "-v", "1"]
    LOG.info("running mmseqs easy-search")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"mmseqs failed ({proc.returncode}): {proc.stderr[-3000:]}")
    hits = pd.read_csv(out_tsv, sep="\t", names=cols)
    # mmseqs reports pident as a fraction in recent versions and as a
    # percentage in older ones; normalise so downstream thresholds are stable.
    if hits["pident"].max() > 1.5:
        hits["pident"] = hits["pident"] / 100.0
    LOG.info("mmseqs returned %d hits for %d queries",
             len(hits), hits["query"].nunique())
    return hits


def assign_paralogues(hits: pd.DataFrame, meta: pd.DataFrame) -> list[Assignment]:
    """Call one paralogue per query from its reference-panel hits.

    A query is accepted only when its best reference hit is a target paralogue,
    identity and coverage clear the floors, and the best bitscore exceeds the
    best bitscore of any *different* paralogue by MIN_BITSCORE_MARGIN. The
    margin is what prevents the paralogue mixing that invalidated the previous
    analysis.
    """
    hits = hits.copy()
    hits["ref_name"] = hits["target"].str.split("|").str[0]
    meta_by_acc = meta.set_index("accession")[["taxid", "organism",
                                               "gene_primary", "length"]].to_dict("index")

    out: list[Assignment] = []
    for acc, grp in hits.groupby("query"):
        grp = grp.sort_values("bits", ascending=False)
        best = grp.iloc[0]
        others = grp[grp["ref_name"] != best["ref_name"]]
        runner = others.iloc[0] if len(others) else None

        margin = (best["bits"] / runner["bits"]) if runner is not None and runner["bits"] > 0 else float("inf")
        info = meta_by_acc.get(acc, {})

        accepted, reason = True, "ok"
        if best["ref_name"] not in TARGETS:
            accepted, reason = False, f"best match is non-target paralogue {best['ref_name']}"
        elif best["tcov"] < MIN_TARGET_COVERAGE:
            accepted, reason = False, f"reference coverage {best['tcov']:.2f} < {MIN_TARGET_COVERAGE}"
        elif best["pident"] < MIN_IDENTITY:
            accepted, reason = False, f"identity {best['pident']:.2f} < {MIN_IDENTITY}"
        elif margin < MIN_BITSCORE_MARGIN:
            accepted, reason = False, (f"ambiguous: {best['ref_name']} vs "
                                       f"{runner['ref_name']} margin {margin:.2f}")

        out.append(Assignment(
            accession=acc,
            taxid=str(info.get("taxid", "")),
            organism=str(info.get("organism", "")),
            claimed_gene=str(info.get("gene_primary", "") or ""),
            length=int(info.get("length", 0) or 0),
            best_reference=best["ref_name"],
            best_identity=round(float(best["pident"]), 4),
            best_target_coverage=round(float(best["tcov"]), 4),
            best_bitscore=float(best["bits"]),
            runner_up_reference=str(runner["ref_name"]) if runner is not None else "",
            runner_up_bitscore=float(runner["bits"]) if runner is not None else 0.0,
            bitscore_margin=round(margin, 3) if margin != float("inf") else -1.0,
            assigned_paralogue=best["ref_name"] if accepted else "",
            accepted=accepted,
            reason=reason,
        ))
    return out


def select_one_per_species(assignments: list[Assignment],
                           seqs: dict[str, str]) -> dict[str, dict[str, str]]:
    """Keep the best sequence per (paralogue, species).

    Ranked by reference coverage first and length second, so a complete but
    shorter orthologue beats a long fragment - the failure mode of the previous
    pipeline's 'longest protein per species' rule.
    """
    chosen: dict[str, dict[str, Assignment]] = {p: {} for p in TARGETS}
    for a in assignments:
        if not a.accepted or a.accession not in seqs:
            continue
        slot = chosen[a.assigned_paralogue]
        incumbent = slot.get(a.taxid)
        key = (a.best_target_coverage, a.length)
        if incumbent is None or key > (incumbent.best_target_coverage, incumbent.length):
            slot[a.taxid] = a
    return {p: {t: a.accession for t, a in d.items()} for p, d in chosen.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scoping-raw", type=Path,
                        default=Path("results/manuscript_analysis/depth_scoping_raw.csv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--workdir", type=Path, default=Path("work/step03"))
    parser.add_argument("--tier", default="Metazoa",
                        help="Metazoa is the superset; paralogue assignment restricts it.")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    started = time.time()
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    seqdir = args.outdir / "sequences"
    seqdir.mkdir(exist_ok=True)

    raw = pd.read_csv(args.scoping_raw, dtype={"taxid": str})
    sel = raw[(raw.route == "interpro") & (raw.tier == args.tier)].drop_duplicates("accession")
    LOG.info("InterPro route, tier=%s: %d unique entries", args.tier, len(sel))

    refs = build_reference_fasta(args.workdir / "references.fasta")

    seqs = dcalib.fetch_sequences(sel.accession.tolist(), batch=100)
    LOG.info("retrieved %d of %d sequences", len(seqs), len(sel))
    query_fasta = args.workdir / "queries.fasta"
    dcalib.write_fasta({a: s for a, s in seqs.items()}, query_fasta)

    hits = run_mmseqs(query_fasta, args.workdir / "references.fasta",
                      args.workdir / "hits.tsv", args.workdir / "tmp", args.threads)
    assignments = assign_paralogues(hits, sel)

    df = pd.DataFrame([asdict(a) for a in assignments])
    df.to_csv(args.outdir / "paralogue_assignment.csv", index=False)

    picked = select_one_per_species(assignments, seqs)
    for paralogue, mapping in picked.items():
        recs = {f"{paralogue}|{t}|{a}": seqs[a] for t, a in mapping.items()}
        dcalib.write_fasta(recs, seqdir / f"{paralogue}.fasta")
        LOG.info("%-6s %d species", paralogue, len(recs))

    summary = {
        "tier": args.tier,
        "entries_retrieved": int(len(sel)),
        "sequences_fetched": int(len(seqs)),
        "queries_with_hits": int(df.accession.nunique()),
        "accepted": int(df.accepted.sum()),
        "rejected": int((~df.accepted).sum()),
        "rejection_reasons": df.loc[~df.accepted, "reason"].str.split(":").str[0]
                               .value_counts().to_dict(),
        "best_reference_counts": df.best_reference.value_counts().to_dict(),
        "species_per_paralogue": {p: len(m) for p, m in picked.items()},
        "runtime_s": round(time.time() - started, 1),
    }
    (args.outdir / "retrieval_summary.json").write_text(json.dumps(summary, indent=2))
    LOG.info("done in %.1f s", time.time() - started)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
