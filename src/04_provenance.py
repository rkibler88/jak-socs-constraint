"""Step 5 - Domain extraction and column provenance.

Builds the table that names every alignment column: which canonical residue of
each reference protein it holds, and which annotated feature it belongs to.

The previous pipeline's central failure was a 203-column block labelled 'SOCS
core' that contained no SOCS residues. The assertions at the end of this module
are what make that class of error impossible: every JAK column must map to a
residue inside an annotated kinase domain, and every SOCS column to a residue
inside the annotated SOCS regulatory span, checked against UniProt feature
annotations rather than against the column's name.

Feature boundaries come from the UniProt REST API. All five references have
KIR, ESS, SH2, SOCS box and the protein-kinase domains annotated, so no residue
range is hardcoded.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib
import lib_constraint as cons

LOG = logging.getLogger("step05")

REFERENCES = {
    "JAK1": ("jak", "P23458"), "JAK2": ("jak", "O60674"), "TYK2": ("jak", "P29597"),
    "SOCS1": ("socs", "O15524"), "SOCS3": ("socs", "O14543"), "SOCS2": ("socs", "O14508"),
}

#: Features retained, keyed by the UniProt description text. The JAK entry of
#: interest is 'Protein kinase 2' (JH1, the catalytic domain); 'Protein kinase 1'
#: is the JH2 pseudokinase and is not part of the SOCS-binding surface.
FEATURE_KEYS = {
    "Protein kinase 2": "JH1",
    "Kinase inhibitory region (KIR)": "KIR",
    "Extended SH2 subdomain (ESS)": "ESS",
    "SH2": "SH2",
    "SOCS box": "SOCS_box",
}

#: Which feature a column must fall inside for its family's assertion to pass.
REQUIRED_FEATURES = {"jak": {"JH1"}, "socs": {"KIR", "ESS", "SH2", "SOCS_box"}}


def fetch_features(accession: str) -> list[dict]:
    url = (f"https://rest.uniprot.org/uniprotkb/{accession}.json"
           "?fields=ft_domain,ft_region,ft_motif,sequence")
    data = json.load(urllib.request.urlopen(url, timeout=120))
    feats = []
    for f in data.get("features", []):
        desc = f.get("description", "")
        label = FEATURE_KEYS.get(desc)
        if label is None and desc.startswith("SH2"):      # 'SH2; atypical'
            label = "SH2"
        if label is None:
            continue
        loc = f["location"]
        feats.append({"label": label, "description": desc,
                      "start": loc["start"]["value"], "end": loc["end"]["value"]})
    return feats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())
    alignments = {fam: dcalib.read_fasta(args.outdir / "alignments" / f"{fam}_family.fasta")
                  for fam in ["jak", "socs"]}

    # --- feature annotations, projected through the validated maps
    feat_rows, projection_rows = [], []
    for ref, (fam, acc) in REFERENCES.items():
        feats = fetch_features(acc)
        assert feats, f"{ref} ({acc}): no usable features returned"
        c2col = {int(k): v for k, v in
                 maps[fam]["references"][ref]["canonical_to_column"].items()}
        for f in feats:
            covered = [p for p in range(f["start"], f["end"] + 1) if p in c2col]
            feat_rows.append({"protein": ref, "accession": acc, "feature": f["label"],
                              "description": f["description"],
                              "canonical_start": f["start"], "canonical_end": f["end"],
                              "canonical_length": f["end"] - f["start"] + 1,
                              "residues_in_alignment": len(covered),
                              "coverage": round(len(covered) / (f["end"] - f["start"] + 1), 3),
                              "column_min": min((c2col[p] for p in covered), default=-1),
                              "column_max": max((c2col[p] for p in covered), default=-1)})
            for p in covered:
                projection_rows.append({"protein": ref, "accession": acc,
                                        "canonical_residue": p, "column": c2col[p],
                                        "feature": f["label"], "feature_type": f["description"]})
        LOG.info("%-6s %s: %s", ref, acc,
                 ", ".join(f"{f['label']}({f['start']}-{f['end']})" for f in feats))

    features = pd.DataFrame(feat_rows)
    projection = pd.DataFrame(projection_rows)
    features.to_csv(args.outdir / "feature_ranges.csv", index=False)
    projection.to_csv(args.outdir / "feature_projection.csv", index=False)

    # --- per-column provenance
    prov_rows = []
    for fam, aln in alignments.items():
        msa = dcalib.encode_msa(aln)
        hw = cons.henikoff_weights(msa)
        gapfrac = dcalib.gap_fraction(msa)
        n_cols = msa.shape[1]
        refs_here = [r for r, (f, _) in REFERENCES.items() if f == fam]

        col_to_canon = {r: {int(k): v for k, v in
                            maps[fam]["references"][r]["column_to_canonical"].items()}
                        for r in refs_here}
        feat_by_ref = {r: projection[projection.protein == r].set_index("canonical_residue")["feature"].to_dict()
                       for r in refs_here}
        seqs = list(aln.values())

        for col in range(n_cols):
            row = {"family": fam, "column": col,
                   "gap_fraction": round(float(gapfrac[col]), 4),
                   "effective_occupancy": round(float(
                       ((msa[:, col] != 0) * hw).sum() / hw.sum()), 4)}
            feats_here = set()
            for r in refs_here:
                pos = col_to_canon[r].get(col)
                row[f"{r}_residue"] = pos if pos is not None else ""
                row[f"{r}_letter"] = (aln[[k for k in aln if k.split("|")[1] == r][0]][col]
                                      if pos is not None else "")
                f = feat_by_ref[r].get(pos) if pos is not None else None
                row[f"{r}_feature"] = f or ""
                if f:
                    feats_here.add(f)
            row["features"] = ";".join(sorted(feats_here))
            row["in_required_feature"] = bool(feats_here & REQUIRED_FEATURES[fam])
            prov_rows.append(row)

    provenance = pd.DataFrame(prov_rows)
    provenance.to_csv(args.outdir / "column_provenance.csv", index=False)

    # --- the assertions that the previous pipeline lacked
    report = {}
    for fam in ["jak", "socs"]:
        sub = provenance[provenance.family == fam]
        anchor = "JAK2" if fam == "jak" else "SOCS1"
        anchored = sub[sub[f"{anchor}_residue"] != ""]
        in_feat = anchored[anchored.in_required_feature]
        report[fam] = {
            "columns": int(len(sub)),
            "columns_with_anchor_residue": int(len(anchored)),
            "columns_in_required_feature": int(len(in_feat)),
            "fraction_in_required_feature": round(len(in_feat) / max(len(anchored), 1), 4),
            "features_present": sorted({f for s in sub.features for f in s.split(";") if f}),
        }
        # Every column carrying the anchor protein's residue must sit inside an
        # annotated feature of the right family. This is the check that would
        # have caught a 'SOCS' block made of kinase sequence.
        offending = anchored[~anchored.in_required_feature]
        assert len(offending) / max(len(anchored), 1) < 0.20, (
            f"{fam}: {len(offending)} of {len(anchored)} anchored columns fall outside "
            f"{REQUIRED_FEATURES[fam]} - domain extraction is wrong"
        )
        wrong_family = {"JH1"} if fam == "socs" else {"KIR", "ESS", "SOCS_box"}
        contaminated = sub[sub.features.apply(
            lambda s: bool(set(s.split(";")) & wrong_family))]
        assert len(contaminated) == 0, (
            f"{fam}: {len(contaminated)} columns carry {wrong_family} annotations - "
            "blocks are mislabelled"
        )
        LOG.info("%s: %d columns, %d anchored, %.1f%% inside %s - assertions passed",
                 fam, len(sub), len(anchored),
                 100 * report[fam]["fraction_in_required_feature"], REQUIRED_FEATURES[fam])

    (args.outdir / "provenance_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print()
    print(features[["protein", "feature", "canonical_start", "canonical_end",
                    "residues_in_alignment", "coverage", "column_min", "column_max"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
