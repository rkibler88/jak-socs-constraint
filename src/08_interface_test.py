"""Step 9 - Interface constraint test against matched controls.

The primary hypothesis test: are residues that contact the partner protein in
4GL9 / 6C7Y under stronger purifying selection than comparable residues that do
not?

'Comparable' is doing the work. Interface residues are not a random sample of
the protein - they sit in particular domains, at particular solvent
accessibilities, in better-aligned columns. A naive interface-vs-rest comparison
would recover those differences rather than selection. Controls are therefore
matched on three covariates simultaneously:

  domain membership     (an interface residue in the SH2 is compared to other
                         SH2 residues, never to a disordered linker)
  solvent accessibility (from the unbound chain, binned - buried residues are
                         conserved whether or not they are in an interface)
  alignment gap content (a well-occupied column scores differently from a
                         sparsely occupied one)

Structure residues are mapped to alignment columns through the structure's own
organism (mouse JAK2/SOCS3 for 4GL9, chicken SOCS1 for 6C7Y), whose coordinate
maps were built and validated in step 4.
"""

from __future__ import annotations

import argparse
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

LOG = logging.getLogger("step09")

SEED = 42
N_PERM = 10000

#: Which coordinate map each structure chain is read through, and which
#: constraint scope (paralogue orthologue set) the test is run in.
STRUCT_MAP = {
    ("4GL9", "SOCS"): {"map": "SOCS3_MOUSE_4GL9", "family": "socs", "scope": "SOCS3"},
    ("4GL9", "JAK"):  {"map": "JAK2_MOUSE_4GL9",  "family": "jak",  "scope": "JAK2"},
    ("6C7Y", "SOCS"): {"map": "SOCS1_CHICK_6C7Y", "family": "socs", "scope": "SOCS1"},
    ("6C7Y", "JAK"):  {"map": "JAK1",             "family": "jak",  "scope": "JAK1"},
}

N_SASA_BINS = 3
N_GAP_BINS = 3


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta: P(a>b) - P(a<b). Non-parametric, no distributional assumption."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    # Rank-based computation avoids the O(n*m) pairwise matrix.
    combined = np.concatenate([a, b])
    ranks = stats.rankdata(combined)
    ra = ranks[:len(a)].sum()
    u = ra - len(a) * (len(a) + 1) / 2.0
    return float(2.0 * u / (len(a) * len(b)) - 1.0)


def bootstrap_delta(a: np.ndarray, b: np.ndarray, n_boot: int = 5000,
                    seed: int = SEED) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for i in range(n_boot):
        out[i] = cliffs_delta(rng.choice(a, len(a), replace=True),
                             rng.choice(b, len(b), replace=True))
    return float(np.nanquantile(out, 0.025)), float(np.nanquantile(out, 0.975))


def matched_permutation(values: pd.Series, is_interface: pd.Series,
                        strata: pd.Series, n_perm: int = N_PERM,
                        seed: int = SEED) -> tuple[float, float, np.ndarray]:
    """Permute interface labels *within* covariate strata.

    This is the test's null: it asks whether interface residues are more
    constrained than other residues *of the same domain, burial and gap class*.
    Shuffling globally would instead ask whether they differ from the protein
    average, which they would for reasons that have nothing to do with binding.
    """
    rng = np.random.default_rng(seed)
    observed = values[is_interface].mean() - values[~is_interface].mean()
    df = pd.DataFrame({"v": values.values, "i": is_interface.values, "s": strata.values})
    null = np.empty(n_perm)
    groups = [g for _, g in df.groupby("s")]
    for p in range(n_perm):
        lab = np.concatenate([rng.permutation(g["i"].values) for g in groups])
        vals = np.concatenate([g["v"].values for g in groups])
        null[p] = vals[lab].mean() - vals[~lab].mean()
    p_emp = ((np.abs(null) >= abs(observed)).sum() + 1) / (n_perm + 1)
    return float(observed), float(p_emp), null


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--n-perm", type=int, default=N_PERM)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())
    provenance = pd.read_csv(args.outdir / "column_provenance.csv")
    iface = pd.read_csv(args.outdir / "interface_residues.csv")
    sasa = pd.read_csv(args.outdir / "structural_sasa.csv")
    constraint = pd.read_csv(args.outdir / "residue_constraint.csv")

    rows, nulls, control_rows = [], {}, []
    for (pdb, side), spec in STRUCT_MAP.items():
        fam, scope, mapname = spec["family"], spec["scope"], spec["map"]
        c2col = {int(k): v for k, v in
                 maps[fam]["references"][mapname]["canonical_to_column"].items()}

        sel = iface[(iface.structure == pdb) & (iface.side == side)]
        if not len(sel):
            continue
        acc = sel.accession.iloc[0]
        iface_cols = {c2col[r] for r in sel.canonical_residue if r in c2col}
        mapped = len(iface_cols)

        # SASA of the unbound chain, for the burial covariate
        s = sasa[(sasa.structure == pdb) & (sasa.accession == acc)]
        sasa_by_col = {c2col[r.canonical_residue]: r.sasa_unbound
                       for r in s.itertuples() if r.canonical_residue in c2col}

        cons_scope = constraint[(constraint.family == fam) & (constraint.scope == scope)]
        cons_by_col = cons_scope.set_index("column")["jsd"].to_dict()
        prov = provenance[provenance.family == fam].set_index("column")

        recs = []
        for col, sasa_val in sasa_by_col.items():
            if col not in cons_by_col or col not in prov.index:
                continue
            feats = str(prov.loc[col, "features"])
            recs.append({"column": col, "jsd": cons_by_col[col],
                         "sasa": sasa_val,
                         "gap_fraction": float(prov.loc[col, "gap_fraction"]),
                         "domain": feats if feats and feats != "nan" else "none",
                         "is_interface": col in iface_cols})
        d = pd.DataFrame(recs)
        if d.is_interface.sum() < 5 or (~d.is_interface).sum() < 5:
            LOG.warning("%s %s: too few residues (%d interface, %d control), skipped",
                        pdb, side, int(d.is_interface.sum()), int((~d.is_interface).sum()))
            continue

        # Covariate strata: domain x SASA tertile x gap tertile
        d["sasa_bin"] = pd.qcut(d.sasa.rank(method="first"), N_SASA_BINS, labels=False)
        d["gap_bin"] = pd.qcut(d.gap_fraction.rank(method="first"), N_GAP_BINS, labels=False)
        d["stratum"] = d.domain + "|" + d.sasa_bin.astype(str) + "|" + d.gap_bin.astype(str)
        # Drop strata that cannot contribute a comparison
        usable = d.groupby("stratum").is_interface.transform(lambda x: 0 < x.sum() < len(x))
        dm = d[usable].copy()

        a = dm.loc[dm.is_interface, "jsd"].values
        b = dm.loc[~dm.is_interface, "jsd"].values
        if len(a) < 5 or len(b) < 5:
            LOG.warning("%s %s: matched comparison empty, skipped", pdb, side)
            continue

        delta = cliffs_delta(a, b)
        lo, hi = bootstrap_delta(a, b, seed=SEED)
        mwu = stats.mannwhitneyu(a, b, alternative="two-sided")
        ks = stats.ks_2samp(a, b)
        diff, p_perm, null = matched_permutation(dm.jsd, dm.is_interface, dm.stratum,
                                                 n_perm=args.n_perm, seed=SEED)
        nulls[f"{pdb}_{side}"] = null

        rows.append({
            "structure": pdb, "side": side, "paralogue": sel.paralogue.iloc[0],
            "accession": acc, "constraint_scope": scope,
            "interface_residues_structural": int(sel.canonical_residue.nunique()),
            "interface_columns_mapped": mapped,
            "n_interface_matched": len(a), "n_control_matched": len(b),
            "n_strata": int(dm.stratum.nunique()),
            "mean_jsd_interface": round(float(a.mean()), 4),
            "mean_jsd_control": round(float(b.mean()), 4),
            "mean_difference": round(diff, 4),
            "cliffs_delta": round(delta, 4),
            "delta_ci_low": round(lo, 4), "delta_ci_high": round(hi, 4),
            "mannwhitney_p": float(mwu.pvalue), "ks_p": float(ks.pvalue),
            "permutation_p": float(p_perm),
        })
        for r in dm.itertuples():
            control_rows.append({"structure": pdb, "side": side, "column": r.column,
                                 "jsd": r.jsd, "sasa": r.sasa,
                                 "gap_fraction": r.gap_fraction, "domain": r.domain,
                                 "stratum": r.stratum, "is_interface": r.is_interface})
        LOG.info("%s %-4s %-5s  iface n=%3d  ctrl n=%3d  strata=%2d  "
                 "delta=%+.3f [%+.3f,%+.3f]  p_perm=%.4f",
                 pdb, side, scope, len(a), len(b), dm.stratum.nunique(), delta, lo, hi, p_perm)

    res = pd.DataFrame(rows)
    q = cons.benjamini_hochberg(res.permutation_p.values)
    res["q_bh"] = np.round(q, 5)
    res["significant_q05"] = res.q_bh < 0.05
    res.to_csv(args.outdir / "interface_constraint_tests.csv", index=False)
    pd.DataFrame(control_rows).to_csv(args.outdir / "matched_control_sets.csv", index=False)
    np.savez_compressed(args.outdir / "interface_nulls.npz", **nulls)

    print(res[["structure", "side", "paralogue", "n_interface_matched", "n_control_matched",
               "mean_jsd_interface", "mean_jsd_control", "cliffs_delta",
               "delta_ci_low", "delta_ci_high", "permutation_p", "q_bh",
               "significant_q05"]].to_string(index=False))


if __name__ == "__main__":
    main()
