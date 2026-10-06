"""Step 12 - Cross-paralogue and cross-structure consistency.

Replication is what converts the step-6/7/10 results from description into
evidence, so it is tested directly rather than asserted:

  (a) Do constraint profiles agree between paralogues of a family? If the
      measurement were noise-dominated they would not.
  (b) Does the interface-constraint effect replicate across two independent
      structures? The four step-9 tests are combined by inverse-variance
      meta-analysis on Cliff's delta, which is the formal version of noticing
      that three of four point the same way.
  (c) Do the declared motif results survive dropping each paralogue in turn?
  (d) Do specificity positions recur across the binding and paralogue axes above
      hypergeometric expectation?
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib
import lib_constraint as cons

LOG = logging.getLogger("step12")
SEED = 42

FAMILY_PARALOGUES = {"jak": ["JAK1", "JAK2", "TYK2"], "socs": ["SOCS1", "SOCS3", "SOCS2"]}


def profile_concordance(constraint: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for family, paralogues in FAMILY_PARALOGUES.items():
        for a, b in itertools.combinations(paralogues, 2):
            da = constraint[(constraint.family == family) & (constraint.scope == a)]
            db = constraint[(constraint.family == family) & (constraint.scope == b)]
            m = da.merge(db, on="column", suffixes=("_a", "_b"))
            rho, p = stats.spearmanr(m.jsd_a, m.jsd_b)
            pear = stats.pearsonr(m.jsd_a, m.jsd_b)
            rows.append({"family": family, "paralogue_a": a, "paralogue_b": b,
                         "n_columns": len(m), "spearman_rho": round(float(rho), 4),
                         "spearman_p": float(p), "pearson_r": round(float(pear[0]), 4)})
    return pd.DataFrame(rows)


def meta_analyse_interface(tests: pd.DataFrame) -> pd.DataFrame:
    """Inverse-variance fixed- and random-effects combination of Cliff's delta.

    SE is taken from the bootstrap CI half-width (CI/1.96), which propagates the
    resampling uncertainty already computed in step 9 rather than assuming a
    parametric form for delta.
    """
    d = tests.copy()
    d["se"] = (d.delta_ci_high - d.delta_ci_low) / (2 * 1.959964)
    d = d[d.se > 0]
    w = 1.0 / d.se ** 2
    fixed = float((w * d.cliffs_delta).sum() / w.sum())
    se_fixed = float(np.sqrt(1.0 / w.sum()))
    z = fixed / se_fixed
    p = 2 * stats.norm.sf(abs(z))

    # Cochran's Q for between-study heterogeneity
    q_stat = float((w * (d.cliffs_delta - fixed) ** 2).sum())
    df = len(d) - 1
    p_het = float(stats.chi2.sf(q_stat, df)) if df > 0 else np.nan
    i2 = max(0.0, (q_stat - df) / q_stat) if q_stat > 0 else 0.0

    rows = [{"model": "fixed_effect_all_four", "n_tests": len(d),
             "pooled_delta": round(fixed, 4),
             "ci_low": round(fixed - 1.959964 * se_fixed, 4),
             "ci_high": round(fixed + 1.959964 * se_fixed, 4),
             "z": round(float(z), 3), "p": float(p),
             "cochran_q": round(q_stat, 3), "p_heterogeneity": p_het,
             "i_squared": round(i2, 3)}]

    # Same, restricted to each structure, to show neither drives the result alone
    for pdb, g in d.groupby("structure"):
        wg = 1.0 / g.se ** 2
        f = float((wg * g.cliffs_delta).sum() / wg.sum())
        se = float(np.sqrt(1.0 / wg.sum()))
        rows.append({"model": f"fixed_effect_{pdb}_only", "n_tests": len(g),
                     "pooled_delta": round(f, 4),
                     "ci_low": round(f - 1.959964 * se, 4),
                     "ci_high": round(f + 1.959964 * se, 4),
                     "z": round(f / se, 3), "p": float(2 * stats.norm.sf(abs(f / se))),
                     "cochran_q": np.nan, "p_heterogeneity": np.nan, "i_squared": np.nan})
    return pd.DataFrame(rows)


def leave_one_paralogue_out(outdir: Path, provenance: pd.DataFrame,
                            maps: dict, n_perm: int = 2000) -> pd.DataFrame:
    """Recompute the KIR specificity result dropping each binder in turn.

    The headline specificity finding pools SOCS1 and SOCS3 as 'binders'. If it
    depended on one of them the pooling would be doing the work, so each is
    dropped in turn and the test repeated.
    """
    aln = dcalib.read_fasta(outdir / "alignments" / "socs_family.fasta")
    msa = dcalib.encode_msa(aln)
    keys = list(aln)
    prov = provenance[provenance.family == "socs"]
    kir_cols = sorted(prov.loc[prov.features.fillna("").str.split(";")
                               .apply(lambda x: "KIR" in x), "column"])
    all_cols = sorted(prov.column)
    ctrl_cols = [c for c in all_cols if c not in set(kir_cols)]

    rows = []
    for binders in [["SOCS1", "SOCS3"], ["SOCS1"], ["SOCS3"]]:
        idx, labels = [], []
        for p in binders + ["SOCS2"]:
            sel = [i for i, k in enumerate(keys) if k.startswith(f"{p}|")]
            idx += sel
            labels += [("binder" if p in binders else "nonbinder")] * len(sel)
        idx = np.array(idx)
        labels = np.array(labels)
        sub = msa[idx]
        w = cons.henikoff_weights(sub)
        mi = cons.specificity_mi(sub, w, labels, n_perm=0, seed=SEED)["mi"]

        rng = np.random.default_rng(SEED)
        obs = float(mi[kir_cols].mean())
        pool = mi[ctrl_cols]
        draws = np.array([rng.choice(pool, len(kir_cols), replace=False).mean()
                          for _ in range(n_perm)])
        p = ((np.abs(draws - pool.mean()) >= abs(obs - pool.mean())).sum() + 1) / (n_perm + 1)
        rows.append({"binders_included": "+".join(binders), "n_binder_sequences": int((labels == "binder").sum()),
                     "kir_mean_mi": round(obs, 4), "control_mean_mi": round(float(pool.mean()), 4),
                     "fold": round(obs / float(pool.mean()), 3),
                     "z": round(float((obs - draws.mean()) / draws.std()), 3),
                     "p_perm": float(p)})
        LOG.info("leave-one-out %-12s fold=%.3f p=%.4f", "+".join(binders), rows[-1]["fold"], p)
    return pd.DataFrame(rows)


def sdp_axis_overlap(sdp: pd.DataFrame, top_n: int = 30) -> pd.DataFrame:
    """Do top specificity positions recur across axes more than chance?"""
    rows = []
    socs_axes = [a for a in sdp.axis.unique() if a.startswith("socs")]
    for a, b in itertools.combinations(socs_axes, 2):
        ta = set(sdp[sdp.axis == a].nlargest(top_n, "mi_bits").column)
        tb = set(sdp[sdp.axis == b].nlargest(top_n, "mi_bits").column)
        total = sdp[sdp.axis == a].column.nunique()
        k = len(ta & tb)
        p = stats.hypergeom.sf(k - 1, total, top_n, top_n)
        rows.append({"axis_a": a, "axis_b": b, "top_n": top_n, "n_columns": total,
                     "overlap": k, "expected": round(top_n * top_n / total, 2),
                     "hypergeometric_p": float(p)})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    constraint = pd.read_csv(args.outdir / "residue_constraint.csv")
    tests = pd.read_csv(args.outdir / "interface_constraint_tests.csv")
    sdp = pd.read_csv(args.outdir / "sdp_table.csv")
    provenance = pd.read_csv(args.outdir / "column_provenance.csv")
    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())

    conc = profile_concordance(constraint)
    conc.to_csv(args.outdir / "profile_concordance.csv", index=False)
    meta = meta_analyse_interface(tests)
    meta.to_csv(args.outdir / "interface_meta_analysis.csv", index=False)
    loo = leave_one_paralogue_out(args.outdir, provenance, maps)
    loo.to_csv(args.outdir / "leave_one_paralogue_out.csv", index=False)
    ov = sdp_axis_overlap(sdp)
    ov.to_csv(args.outdir / "sdp_axis_overlap.csv", index=False)

    print("=== constraint profile concordance between paralogues ===")
    print(conc.to_string(index=False))
    print("\n=== interface effect, meta-analysed across structures ===")
    print(meta.to_string(index=False))
    print("\n=== KIR specificity, leave-one-binder-out ===")
    print(loo.to_string(index=False))
    print("\n=== specificity-axis overlap ===")
    print(ov.to_string(index=False))


if __name__ == "__main__":
    main()
