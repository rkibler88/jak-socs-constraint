"""Step 7 - Specificity-determining position analysis.

Locates the variation that exists. Within a paralogue's orthologue set almost
every position is either fixed or nearly so, which is why DCA fails here; the
informative variation is *between* paralogues. This step asks, per column, how
much the residue identity tells you which paralogue you are looking at.

Two axes are analysed:

  socs_binding   SOCS1 + SOCS3 (bind and inhibit JAK) vs SOCS2 (does not),
                 the contrast that makes the KIR/ESS result in step 6
                 residue-addressable.
  socs_paralogue SOCS1 vs SOCS3, which distinguishes the two binders.
  jak_paralogue  JAK1 vs JAK2 vs TYK2.

Statistics: Miller-Madow bias-corrected mutual information between column state
and group label, against a label-permutation null that preserves both the
column's residue composition and the group sizes. Benjamini-Hochberg across
columns within each axis.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib
import lib_constraint as cons

LOG = logging.getLogger("step07")

SEED = 42

#: (family, {group label -> paralogues in it}). Group labels are the biology,
#: not the gene names, so the binding axis pools the two inhibitory paralogues.
AXES = {
    "socs_binding":   ("socs", {"binder": ["SOCS1", "SOCS3"], "nonbinder": ["SOCS2"]}),
    "socs_paralogue": ("socs", {"SOCS1": ["SOCS1"], "SOCS3": ["SOCS3"]}),
    "jak_paralogue":  ("jak",  {"JAK1": ["JAK1"], "JAK2": ["JAK2"], "TYK2": ["TYK2"]}),
}

ANCHOR = {"socs": "SOCS1", "jak": "JAK2"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--n-perm", type=int, default=10000)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())
    provenance = pd.read_csv(args.outdir / "column_provenance.csv")

    alignments, encoded = {}, {}
    for fam in ["jak", "socs"]:
        aln = dcalib.read_fasta(args.outdir / "alignments" / f"{fam}_family.fasta")
        alignments[fam] = aln
        encoded[fam] = dcalib.encode_msa(aln)

    rows, null_store = [], {}
    for axis, (family, groups) in AXES.items():
        aln, msa = alignments[family], encoded[family]
        keys = list(aln)
        idx, labels = [], []
        for label, paralogues in groups.items():
            for p in paralogues:
                sel = [i for i, k in enumerate(keys) if k.startswith(f"{p}|")]
                idx.extend(sel)
                labels.extend([label] * len(sel))
        idx = np.array(idx)
        labels = np.array(labels)
        sub = msa[idx]
        weights = cons.henikoff_weights(sub)

        LOG.info("%s: %d sequences, groups %s", axis, len(idx),
                 {g: int((labels == g).sum()) for g in groups})
        res = cons.specificity_mi(sub, weights, labels,
                                  n_perm=args.n_perm, seed=SEED)
        contrast = cons.group_contrast(sub, weights, labels, dcalib.ALPHABET)
        q = cons.benjamini_hochberg(res["p_empirical"])
        null_store[f"{axis}_null_max"] = res["null_max"]
        null_store[f"{axis}_null_mean"] = res["null_mean"]
        null_store[f"{axis}_null_sd"] = res["null_sd"]
        null_store[f"{axis}_observed"] = res["mi"]

        anchor = ANCHOR[family]
        col2canon = {int(k): v for k, v in
                     maps[family]["references"][anchor]["column_to_canonical"].items()}
        prov = provenance[provenance.family == family].set_index("column")
        group_names = list(contrast["groups"])
        consensus = contrast["consensus_per_group"]
        letters = dcalib.ALPHABET

        for col in range(sub.shape[1]):
            cons_res = {g: letters[consensus[gi, col] + 1] if consensus[gi, col] + 1 < len(letters)
                        else "?" for gi, g in enumerate(group_names)}
            rows.append({
                "axis": axis, "family": family, "column": col,
                f"{anchor}_residue": col2canon.get(col, ""),
                "anchor_letter": prov.loc[col, f"{anchor}_letter"] if col in prov.index else "",
                "feature": prov.loc[col, "features"] if col in prov.index else "",
                "mi_bits": round(float(res["mi"][col]), 5),
                "null_mean": round(float(res["null_mean"][col]), 5),
                "z": round(float(res["z"][col]), 3),
                "p_empirical": float(res["p_empirical"][col]),
                "q_bh": round(float(q[col]), 5),
                "significant": bool(q[col] < 0.05),
                "within_group_conservation": round(float(contrast["within_group_conservation"][col]), 4),
                "min_within_group_conservation": round(float(contrast["min_within_group_conservation"][col]), 4),
                "n_distinct_consensus": int(contrast["n_distinct_consensus"][col]),
                "consensus": "/".join(f"{g}:{cons_res[g]}" for g in group_names),
                "gap_fraction": round(float(prov.loc[col, "gap_fraction"]) if col in prov.index else np.nan, 4),
            })

    sdp = pd.DataFrame(rows)
    sdp.to_csv(args.outdir / "sdp_table.csv", index=False)
    np.savez_compressed(args.outdir / "sdp_nulls.npz", **null_store)

    summary = []
    for axis, g in sdp.groupby("axis"):
        summary.append({"axis": axis, "columns": len(g),
                        "significant_q05": int(g.significant.sum()),
                        "median_mi": round(float(g.mi_bits.median()), 4),
                        "max_mi": round(float(g.mi_bits.max()), 4)})
    print(pd.DataFrame(summary).to_string(index=False))

    # Where do the binding-axis SDPs concentrate?
    b = sdp[sdp.axis == "socs_binding"]
    feat = (b.assign(feature=b.feature.fillna("none").replace("", "none"))
             .groupby("feature")
             .agg(n_columns=("column", "size"), n_significant=("significant", "sum"),
                  mean_mi=("mi_bits", "mean"))
             .assign(mean_mi=lambda d: d.mean_mi.round(4)))
    feat["frac_significant"] = (feat.n_significant / feat.n_columns).round(3)
    feat.to_csv(args.outdir / "sdp_feature_summary.csv")
    print()
    print(feat.to_string())
    print("\ntop 20 binding-axis SDPs by MI:")
    cols = ["SOCS1_residue", "anchor_letter", "feature", "mi_bits", "z", "q_bh", "consensus"]
    print(b.nlargest(20, "mi_bits")[cols].to_string(index=False))


if __name__ == "__main__":
    main()
