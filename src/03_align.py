"""Step 4 - Per-family alignment and gap-aware coordinate mapping.

One alignment per family, with all paralogues of that family in it, so that
cross-paralogue comparisons share a single column coordinate system. Concatenated
records are never aligned - doing so is what let the previous pipeline slice a
'SOCS' block that contained no SOCS residues.

Two-pass alignment: a fast full-length pass locates the domain of interest via
the human reference coordinate map, then each sequence's own residues from that
span are realigned with iterative refinement. Aligning 1,150-residue full-length
kinases to extract a 276-residue domain wastes accuracy exactly where it matters.

Sequences whose numbering the structural analysis depends on (mouse JAK2 and
SOCS3 for 4GL9, chicken SOCS1 for 6C7Y) are forced into the alignments, so the
structure -> column mapping is exact rather than inferred from a near relative.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib

LOG = logging.getLogger("step04")

#: Human canonical references, which anchor every coordinate map.
HUMAN_REFS = {
    "JAK1": "P23458", "JAK2": "O60674", "TYK2": "P29597",
    "SOCS1": "O15524", "SOCS3": "O14543", "SOCS2": "O14508",
}

#: Sequences the structural analysis maps through. Forced into the alignment.
STRUCTURAL_REFS = {
    "Q62120": ("JAK", "JAK2_MOUSE_4GL9"),
    "O35718": ("SOCS", "SOCS3_MOUSE_4GL9"),
    "B6RCQ2": ("SOCS", "SOCS1_CHICK_6C7Y"),
}

FAMILY_MEMBERS = {
    "jak": ["JAK1", "JAK2", "TYK2"],
    "socs": ["SOCS1", "SOCS3", "SOCS2"],
}

#: Domain span used to define the first-pass extraction window, in the canonical
#: numbering of the named human reference. Widened by FLANK so the second-pass
#: alignment has context beyond the annotated boundary.
DOMAIN_WINDOW = {
    "jak": ("JAK2", 849, 1124),     # Protein kinase 2 (JH1)
    "socs": ("SOCS1", 40, 215),     # KIR through SOCS box, with flanks
}
FLANK = 15


def load_paralogue_sets(seqdir: Path, assignment: pd.DataFrame,
                        socs2_from_assignment: bool = True) -> dict[str, dict[str, str]]:
    """Load the per-paralogue FASTA files written by step 3.

    SOCS2 is not among step 3's accepted targets (it is a decoy that identifies
    off-target family members), but it is the non-binding contrast for the
    specificity analysis, so it is reconstructed here from the assignment table
    using the same acceptance thresholds.
    """
    sets: dict[str, dict[str, str]] = {}
    for path in sorted(seqdir.glob("*.fasta")):
        name = path.stem
        recs = dcalib.read_fasta(path)
        sets[name] = {k.split("|")[1]: v for k, v in recs.items()}
        LOG.info("%-6s %d species from %s", name, len(sets[name]), path.name)

    if socs2_from_assignment:
        ok = assignment[(assignment.best_reference == "SOCS2")
                        & (assignment.best_target_coverage >= 0.55)
                        & (assignment.best_identity >= 0.35)
                        & (assignment.bitscore_margin >= 1.10)]
        best = (ok.sort_values(["best_target_coverage", "length"], ascending=False)
                  .drop_duplicates("taxid"))
        accs = best.accession.tolist()
        seqs = dcalib.fetch_sequences(accs, batch=100)
        sets["SOCS2"] = {str(r.taxid): seqs[r.accession]
                         for r in best.itertuples() if r.accession in seqs}
        LOG.info("SOCS2  %d species (reconstructed as specificity contrast)", len(sets["SOCS2"]))
    return sets


def assemble_family(family: str, sets: dict[str, dict[str, str]]) -> dict[str, str]:
    """Build the record set for one family, forcing in the required references."""
    recs: dict[str, str] = {}
    for paralogue in FAMILY_MEMBERS[family]:
        for taxid, seq in sets.get(paralogue, {}).items():
            recs[f"{paralogue}|{taxid}"] = seq

    # Human canonical references must be present and must be the canonical
    # accession, not whichever entry step 3 happened to pick for taxid 9606.
    needed = {p: a for p, a in HUMAN_REFS.items() if p in FAMILY_MEMBERS[family]}
    extra = dcalib.fetch_sequences(list(needed.values()), batch=20)
    for paralogue, acc in needed.items():
        recs[f"REF|{paralogue}|{acc}"] = extra[acc]

    struct = {acc: label for acc, (fam, label) in STRUCTURAL_REFS.items()
              if fam.lower() == family.lower()}
    if struct:
        sseqs = dcalib.fetch_sequences(list(struct), batch=20)
        for acc, label in struct.items():
            recs[f"STRUCT|{label}|{acc}"] = sseqs[acc]
    LOG.info("family %s: %d records assembled", family, len(recs))
    return recs


def two_pass_align(family: str, recs: dict[str, str], workdir: Path,
                   threads: int) -> tuple[dict[str, str], dict]:
    """Full-length pass to locate the domain, then realign the extracted span."""
    ref_name, start, end = DOMAIN_WINDOW[family]
    ref_key = f"REF|{ref_name}|{HUMAN_REFS[ref_name]}"

    LOG.info("pass 1: full-length alignment of %d %s records", len(recs), family)
    pass1 = dcalib.run_mafft(recs, workdir / f"{family}_pass1.fasta",
                             mode="auto", threads=threads)
    cmap1 = dcalib.build_coordinate_map(ref_name, pass1[ref_key], recs[ref_key])
    lo = max(1, start - FLANK)
    hi = min(len(recs[ref_key]), end + FLANK)
    cols = dcalib.project_range(cmap1, lo, hi)
    colset = set(cols)

    # Each sequence contributes its own residues occupying the window, ungapped.
    # For the reference rows we additionally record, per extracted residue, its
    # canonical position in that protein. Recovering it later by substring search
    # would be wrong: a paralogue's residues inside a window defined on another
    # paralogue are not contiguous in its own sequence wherever it carries an
    # insertion, so the extracted string is not a substring of its canonical form.
    extracted: dict[str, str] = {}
    canonical_positions: dict[str, list[int]] = {}
    ref_keys = [k for k in pass1 if k.startswith(("REF|", "STRUCT|"))]
    ref_cmaps = {k: dcalib.build_coordinate_map(k.split("|")[1], pass1[k], recs[k])
                 for k in ref_keys}

    for key, aligned in pass1.items():
        taken = [c for c in cols if aligned[c] not in ("-", ".")]
        sub = "".join(aligned[c] for c in taken)
        if len(sub) < 0.5 * (hi - lo):
            continue
        extracted[key] = sub
        if key in ref_cmaps:
            canonical_positions[key] = [ref_cmaps[key].msa_to_human[c] for c in taken]

    LOG.info("pass 1 window: canonical %d-%d of %s -> %d columns; "
             "%d/%d records retained at >=50%% occupancy",
             lo, hi, ref_name, len(cols), len(extracted), len(pass1))
    assert ref_key in extracted, f"{ref_key} lost during extraction"

    LOG.info("pass 2: iterative-refinement alignment of extracted domains")
    pass2 = dcalib.run_mafft(extracted, workdir / f"{family}_domain.fasta",
                             mode="nsi", threads=threads)
    meta = {"family": family, "window_reference": ref_name,
            "window_canonical": [lo, hi], "pass1_columns": len(cols),
            "records_in": len(recs), "records_retained": len(extracted),
            "pass2_columns": len(next(iter(pass2.values())))}
    return pass2, extracted, canonical_positions, meta


def build_maps(aligned: dict[str, str], recs: dict[str, str],
               extracted: dict[str, str],
               canonical_positions: dict[str, list[int]]
               ) -> dict[str, dict[str, dict[int, int]]]:
    """Canonical-residue <-> alignment-column maps for every reference row.

    Composes the pass-2 map (extracted-fragment position -> column) with the
    per-residue canonical positions recorded during extraction, giving a direct
    bijection between canonical UniProt numbering and final alignment column.
    Every mapping is validated by asserting the column holds the residue the
    canonical sequence has at that position.
    """
    maps: dict[str, dict[str, dict[int, int]]] = {}
    for key, positions in canonical_positions.items():
        label = key.split("|")[1]
        frag = extracted[key]
        assert len(frag) == len(positions), f"{key}: fragment/position length mismatch"
        cmap = dcalib.build_coordinate_map(label, aligned[key], frag)

        canon_to_col, col_to_canon = {}, {}
        full = recs[key]
        for i, canon_pos in enumerate(positions, start=1):
            col = cmap.human_to_msa[i]
            assert aligned[key][col] == full[canon_pos - 1], (
                f"{label}: canonical residue {canon_pos} ({full[canon_pos-1]}) "
                f"maps to column {col} holding {aligned[key][col]!r}"
            )
            canon_to_col[canon_pos] = col
            col_to_canon[col] = canon_pos

        maps[label] = {"canonical_to_column": canon_to_col,
                       "column_to_canonical": col_to_canon}
        LOG.info("map %-20s %d canonical residues, range %d-%d, validated",
                 label, len(canon_to_col), min(canon_to_col), max(canon_to_col))
    return maps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--workdir", type=Path, default=Path("work/step04"))
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    started = time.time()
    args.workdir.mkdir(parents=True, exist_ok=True)
    aln_dir = args.outdir / "alignments"
    aln_dir.mkdir(parents=True, exist_ok=True)

    assignment = pd.read_csv(args.outdir / "paralogue_assignment.csv", dtype={"taxid": str})
    sets = load_paralogue_sets(args.outdir / "sequences", assignment)

    all_meta, quality_rows, map_export = [], [], {}
    for family in ["jak", "socs"]:
        recs = assemble_family(family, sets)
        aligned, extracted, canon_pos, meta = two_pass_align(
            family, recs, args.workdir, args.threads)
        dcalib.write_fasta(aligned, aln_dir / f"{family}_family.fasta")
        maps = build_maps(aligned, recs, extracted, canon_pos)

        msa = dcalib.encode_msa(aligned)
        weights = dcalib.sequence_weights(msa, 0.8)
        gapfrac = dcalib.gap_fraction(msa)
        neff = float(weights.sum())

        # Mean pairwise identity on a bounded subsample - the full matrix is
        # quadratic in N and the estimate is stable well before that.
        rng = np.random.default_rng(42)
        idx = rng.choice(msa.shape[0], min(250, msa.shape[0]), replace=False)
        sub = msa[idx]
        eq = (sub[:, None, :] == sub[None, :, :]).sum(axis=2) / sub.shape[1]
        iu = np.triu_indices(len(idx), 1)

        meta.update({"n_sequences": int(msa.shape[0]), "neff_80": round(neff, 1),
                     "mean_pairwise_identity": round(float(eq[iu].mean()), 4),
                     "median_gap_fraction": round(float(np.median(gapfrac)), 4)})
        all_meta.append(meta)

        for paralogue in FAMILY_MEMBERS[family]:
            sel = [i for i, k in enumerate(aligned) if k.startswith(f"{paralogue}|")]
            if not sel:
                continue
            quality_rows.append({
                "family": family, "paralogue": paralogue, "n_sequences": len(sel),
                "neff_80": round(float(dcalib.sequence_weights(msa[sel], 0.8).sum()), 1),
                "mean_occupancy": round(float((msa[sel] != 0).mean()), 4),
            })

        map_export[family] = {"columns": meta["pass2_columns"], "references": maps}

    (args.outdir / "coordinate_maps.json").write_text(json.dumps(map_export))
    pd.DataFrame(quality_rows).to_csv(args.outdir / "alignment_quality.csv", index=False)
    (args.outdir / "alignment_meta.json").write_text(json.dumps(all_meta, indent=2))
    LOG.info("done in %.1f s", time.time() - started)
    print(json.dumps(all_meta, indent=2))
    print()
    print(pd.DataFrame(quality_rows).to_string(index=False))


if __name__ == "__main__":
    main()
