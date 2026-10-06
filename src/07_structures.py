"""Step 8 - Structural interface definition from 4GL9 and 6C7Y.

Computes the experimental JAK-SOCS interfaces that the constraint analysis is
tested against. Two independent crystal structures cover two different
paralogue pairs and, as it turns out, two complementary parts of the SOCS
molecule:

  4GL9  mouse JAK2 kinase domain + mouse SOCS3 (resolved 38-128: ESS and SH2,
        KIR not modelled)                    Kershaw et al. 2013, 10.1038/nsmb.2519
  6C7Y  human JAK1 kinase domain + chicken SOCS1 (resolved 48-164: includes the
        KIR)                                 Liau et al. 2018, 10.1038/s41467-018-04013-1

Neither structure uses the human protein throughout, so residue mapping is done
by aligning each observed chain to the canonical sequence of the accession SIFTS
assigns it - never by assuming an offset.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import urllib.request
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.Data.IUPACData import protein_letters_3to1

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib

LOG = logging.getLogger("step08")

CONTACT_CUTOFF_A = 8.0

#: Structure specification. ``chains`` maps chain id -> (UniProt accession,
#: role), taken from the SIFTS uniprot_segments mapping rather than assumed.
STRUCTURES = {
    "4GL9": {
        "jak_chains": {"A": "Q62120", "B": "Q62120", "C": "Q62120", "D": "Q62120"},
        "socs_chains": {"E": "O35718", "F": "O35718", "G": "O35718", "H": "O35718"},
        "jak_paralogue": "JAK2",
        "socs_paralogue": "SOCS3",
        "citation": "Kershaw 2013 10.1038/nsmb.2519",
    },
    "6C7Y": {
        "jak_chains": {"A": "P23458"},
        "socs_chains": {"B": "B6RCQ2"},
        "jak_paralogue": "JAK1",
        "socs_paralogue": "SOCS1",
        "citation": "Liau 2018 10.1038/s41467-018-04013-1",
    },
}

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


@dataclass
class ContactPair:
    structure: str
    jak_chain: str
    socs_chain: str
    jak_paralogue: str
    socs_paralogue: str
    jak_accession: str
    socs_accession: str
    jak_canonical_residue: int
    jak_residue_letter: str
    socs_canonical_residue: int
    socs_residue_letter: str
    min_heavy_atom_distance: float


def download_cif(pdb_id: str, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{pdb_id}.cif"
    if not path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        LOG.info("downloading %s", url)
        urllib.request.urlretrieve(url, path)
    return path


def observed_residues(chain) -> list[tuple[int, str, np.ndarray]]:
    """Return [(auth_seq_id, one-letter code, heavy-atom coords)] for a chain."""
    out = []
    for res in chain:
        if res.id[0] != " ":          # skip waters and heteroatoms
            continue
        code = THREE_TO_ONE.get(res.get_resname().upper())
        if code is None:
            continue
        coords = np.array([a.coord for a in res if a.element != "H"])
        if len(coords) == 0:
            continue
        out.append((res.id[1], code, coords))
    return out


def map_to_canonical(observed: list[tuple[int, str, np.ndarray]],
                     canonical: str, label: str) -> dict[int, int]:
    """Map auth_seq_id -> canonical UniProt position by sequence alignment.

    Uses a global alignment with terminal gaps free, so a construct that covers
    an internal fragment of the canonical protein maps correctly without any
    offset assumption.
    """
    obs_seq = "".join(c for _, c, _ in observed)
    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "global"
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    aln = aligner.align(canonical, obs_seq)[0]

    mapping: dict[int, int] = {}
    identical = 0
    for (c_start, c_end), (o_start, o_end) in zip(*aln.aligned):
        for offset in range(c_end - c_start):
            canon_pos = c_start + offset + 1          # 1-based UniProt
            obs_idx = o_start + offset
            mapping[observed[obs_idx][0]] = canon_pos
            identical += canonical[canon_pos - 1] == obs_seq[obs_idx]

    coverage = len(mapping) / len(observed)
    pct_id = identical / max(len(mapping), 1)
    LOG.info("%s: mapped %d/%d observed residues to canonical (%.1f%% identical)",
             label, len(mapping), len(observed), 100 * pct_id)
    assert coverage > 0.90, f"{label}: only {coverage:.2f} of observed residues mapped"
    assert pct_id > 0.95, f"{label}: alignment only {pct_id:.2f} identical - wrong accession?"
    return mapping


def chain_contacts(jak_obs, socs_obs, jak_map, socs_map, cutoff: float
                   ) -> list[tuple[int, str, int, str, float]]:
    """Minimum heavy-atom distances below ``cutoff`` between two chains."""
    # Coarse pre-filter on residue centroids keeps this O(n*m) loop cheap.
    j_cent = np.array([c.mean(axis=0) for _, _, c in jak_obs])
    s_cent = np.array([c.mean(axis=0) for _, _, c in socs_obs])
    d_cent = np.linalg.norm(j_cent[:, None, :] - s_cent[None, :, :], axis=2)
    candidates = np.argwhere(d_cent < cutoff + 25.0)

    out = []
    for i, k in candidates:
        ji, jc, jcoord = jak_obs[i]
        si, sc, scoord = socs_obs[k]
        d = np.linalg.norm(jcoord[:, None, :] - scoord[None, :, :], axis=2).min()
        if d <= cutoff and ji in jak_map and si in socs_map:
            out.append((jak_map[ji], jc, socs_map[si], sc, float(d)))
    return out


def relative_sasa(structure_path: Path, chain_ids: list[str]) -> dict[tuple[str, int], float]:
    """Per-residue SASA for isolated chains, used to build matched controls."""
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure("s", structure_path)[0]
    sr = ShrakeRupley()
    out: dict[tuple[str, int], float] = {}
    for cid in chain_ids:
        # Copy the chain in isolation so SASA reflects the unbound monomer.
        chain = model[cid].copy()
        sr.compute(chain, level="R")
        for res in chain:
            if res.id[0] == " ":
                out[(cid, res.id[1])] = float(res.sasa)
    return out


def analyse(pdb_id: str, spec: dict, cifdir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    path = download_cif(pdb_id, cifdir)
    model = MMCIFParser(QUIET=True).get_structure(pdb_id, path)[0]

    accessions = sorted(set(spec["jak_chains"].values()) | set(spec["socs_chains"].values()))
    canon = dcalib.fetch_sequences(accessions, batch=20)
    assert all(a in canon for a in accessions), f"missing canonical sequences for {pdb_id}"

    obs, maps = {}, {}
    for cid, acc in {**spec["jak_chains"], **spec["socs_chains"]}.items():
        o = observed_residues(model[cid])
        obs[cid] = o
        maps[cid] = map_to_canonical(o, canon[acc], f"{pdb_id}:{cid}({acc})")

    pairs: list[ContactPair] = []
    per_copy = []
    for jc in spec["jak_chains"]:
        for sc in spec["socs_chains"]:
            cs = chain_contacts(obs[jc], obs[sc], maps[jc], maps[sc], CONTACT_CUTOFF_A)
            per_copy.append({"structure": pdb_id, "jak_chain": jc, "socs_chain": sc,
                             "n_contacts": len(cs)})
            for jr, jl, sr_, sl, d in cs:
                pairs.append(ContactPair(
                    structure=pdb_id, jak_chain=jc, socs_chain=sc,
                    jak_paralogue=spec["jak_paralogue"], socs_paralogue=spec["socs_paralogue"],
                    jak_accession=spec["jak_chains"][jc], socs_accession=spec["socs_chains"][sc],
                    jak_canonical_residue=jr, jak_residue_letter=jl,
                    socs_canonical_residue=sr_, socs_residue_letter=sl,
                    min_heavy_atom_distance=round(d, 2)))

    contacts = pd.DataFrame([asdict(p) for p in pairs])
    copies = pd.DataFrame(per_copy)

    # A crystal with several copies gives an internal reproducibility check:
    # keep the chain pairings that form a real interface, not lattice contacts.
    if len(copies):
        thresh = max(10, 0.25 * copies.n_contacts.max())
        real = copies[copies.n_contacts >= thresh]
        contacts = contacts.merge(real[["jak_chain", "socs_chain"]],
                                  on=["jak_chain", "socs_chain"], how="inner")
        LOG.info("%s: %d of %d chain pairings form an interface (>=%.0f contacts)",
                 pdb_id, len(real), len(copies), thresh)

    sasa = relative_sasa(path, sorted({**spec["jak_chains"], **spec["socs_chains"]}))
    sasa_rows = []
    for (cid, auth), val in sasa.items():
        if auth in maps.get(cid, {}):
            acc = {**spec["jak_chains"], **spec["socs_chains"]}[cid]
            sasa_rows.append({"structure": pdb_id, "chain": cid, "accession": acc,
                              "canonical_residue": maps[cid][auth], "sasa_unbound": round(val, 2)})
    sasa_df = pd.DataFrame(sasa_rows)

    meta = {"pdb_id": pdb_id, "citation": spec["citation"],
            "jak_paralogue": spec["jak_paralogue"], "socs_paralogue": spec["socs_paralogue"],
            "chain_pair_contacts": copies.to_dict("records"),
            "n_contact_pairs_kept": int(len(contacts)),
            "jak_interface_residues": int(contacts.jak_canonical_residue.nunique()) if len(contacts) else 0,
            "socs_interface_residues": int(contacts.socs_canonical_residue.nunique()) if len(contacts) else 0,
            "socs_resolved_range": [int(min(maps[c].values())) for c in spec["socs_chains"]][:1]
                                   + [int(max(maps[c].values())) for c in spec["socs_chains"]][:1]}
    return contacts, sasa_df, meta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--cifdir", type=Path, default=Path("work/structures"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    all_contacts, all_sasa, metas = [], [], []
    for pdb_id, spec in STRUCTURES.items():
        c, s, m = analyse(pdb_id, spec, args.cifdir)
        all_contacts.append(c); all_sasa.append(s); metas.append(m)

    contacts = pd.concat(all_contacts, ignore_index=True)
    sasa = pd.concat(all_sasa, ignore_index=True).drop_duplicates(
        ["structure", "accession", "canonical_residue"])
    contacts.to_csv(args.outdir / "structural_contacts.csv", index=False)
    sasa.to_csv(args.outdir / "structural_sasa.csv", index=False)

    # Interface residue sets, deduplicated across crystal copies.
    rows = []
    for _, g in contacts.groupby(["structure", "socs_paralogue", "socs_canonical_residue"]):
        r = g.iloc[0]
        rows.append({"structure": r.structure, "side": "SOCS", "paralogue": r.socs_paralogue,
                     "accession": r.socs_accession, "canonical_residue": int(r.socs_canonical_residue),
                     "residue_letter": r.socs_residue_letter,
                     "n_partner_residues": int(g.jak_canonical_residue.nunique()),
                     "min_distance": float(g.min_heavy_atom_distance.min())})
    for _, g in contacts.groupby(["structure", "jak_paralogue", "jak_canonical_residue"]):
        r = g.iloc[0]
        rows.append({"structure": r.structure, "side": "JAK", "paralogue": r.jak_paralogue,
                     "accession": r.jak_accession, "canonical_residue": int(r.jak_canonical_residue),
                     "residue_letter": r.jak_residue_letter,
                     "n_partner_residues": int(g.socs_canonical_residue.nunique()),
                     "min_distance": float(g.min_heavy_atom_distance.min())})
    iface = pd.DataFrame(rows).sort_values(["structure", "side", "canonical_residue"])
    iface.to_csv(args.outdir / "interface_residues.csv", index=False)

    (args.outdir / "structures_meta.json").write_text(json.dumps(metas, indent=2))
    print(json.dumps(metas, indent=2))
    print("\ninterface residue counts:")
    print(iface.groupby(["structure", "side", "paralogue"]).size().to_string())


if __name__ == "__main__":
    main()
