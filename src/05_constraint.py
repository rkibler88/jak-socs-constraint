"""Step 6 - Per-residue evolutionary constraint quantification.

Replaces coupling inference as the core measurement. For every alignment column
this computes sequence-weighted Shannon entropy, normalised conservation, and
Jensen-Shannon divergence against an amino-acid background, within each
paralogue group and for each family as a whole, with bootstrap confidence
intervals over the species sample.

Weighting is Henikoff position-based, not the DCA 80%-identity convention: the
latter collapses a set of orthologues of one paralogue almost completely (JAK2:
N_eff 8.8 from 625 sequences) because orthologues of a single paralogue are by
construction highly similar. The DCA-convention N_eff is still reported, in the
power analysis, where it is the relevant quantity.
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

LOG = logging.getLogger("step06")

FAMILY_PARALOGUES = {"jak": ["JAK1", "JAK2", "TYK2"], "socs": ["SOCS1", "SOCS3", "SOCS2"]}
N_BOOT = 1000
N_PERM = 10000
SEED = 42


def group_rows(aln: dict[str, str], paralogue: str) -> list[int]:
    """Row indices for one paralogue's orthologue set.

    Reference and structural rows are excluded so a protein is not counted twice
    (the human canonical is already present as its species row).
    """
    return [i for i, k in enumerate(aln) if k.startswith(f"{paralogue}|")]


def module_deficit_test(
    constraint_df: pd.DataFrame,
    family: str,
    binders: tuple[str, ...],
    contrast: str,
    modules: tuple[str, ...],
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    """Test whether the binder-minus-contrast constraint deficit is module-specific.

    The module means alone do not establish that the contrast paralogue is
    *selectively* less constrained in the inhibitory modules, because it could
    simply be less constrained everywhere. This asks a sharper question: taking
    the per-column difference between the binder mean and the contrast, is that
    difference larger inside a module than across the rest of the alignment?
    Columns are permuted, so the null preserves the observed difference
    distribution and only breaks its association with module membership.

    Three variants are reported, because the contrast paralogue may lack the
    module outright and a deficit that is purely deletion would be circular:

    ``gap_penalised``
        The published per-column score, scaled by non-gap fraction.
    ``residues_only``
        Gap scaling removed, so only divergence among residues actually
        present contributes.
    ``contrast_occupied``
        Restricted to columns the contrast paralogue occupies at >50%, i.e.
        positions it demonstrably retains.

    Args:
        constraint_df: Per-column output of this stage, all scopes.
        family: Family key to test (``"socs"`` or ``"jak"``).
        binders: Scopes pooled as the functional group.
        contrast: Scope used as the non-functional comparator.
        modules: Feature substrings to test, each against all other columns.
        n_perm: Column permutations per module.
        seed: RNG seed.

    Returns:
        One row per (variant, module) with the observed mean difference, the
        null mean, a z score, an empirical p-value and a BH-adjusted q-value
        computed within each variant.
    """
    sub = constraint_df[constraint_df.family == family]
    jsd = sub.pivot_table(index="column", columns="scope", values="jsd")
    gap = sub.pivot_table(index="column", columns="scope", values="gap_fraction")
    feature = (sub[sub.scope == contrast].set_index("column")["feature"]
               .reindex(jsd.index).fillna("").astype(str))

    missing = [s for s in (*binders, contrast) if s not in jsd.columns]
    if missing:
        raise ValueError(f"{family}: scopes absent from constraint table: {missing}")

    # Undo the non-gap scaling to recover divergence among present residues.
    residues_only = (jsd / (1.0 - gap)).replace([np.inf, -np.inf], np.nan)
    occupied = (gap[contrast] < 0.5).values

    variants = {
        "gap_penalised": (jsd[list(binders)].mean(axis=1) - jsd[contrast],
                          np.ones(len(jsd), dtype=bool)),
        "residues_only": (residues_only[list(binders)].mean(axis=1) - residues_only[contrast],
                          np.ones(len(jsd), dtype=bool)),
        "contrast_occupied": (residues_only[list(binders)].mean(axis=1) - residues_only[contrast],
                              occupied),
    }

    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for variant, (delta, eligible) in variants.items():
        values = delta.to_numpy(dtype=float)
        usable = eligible & ~np.isnan(values)
        for module in modules:
            in_module = feature.str.contains(module, regex=False).to_numpy() & usable
            n = int(in_module.sum())
            pool = values[usable]
            if n < 3 or pool.size <= n:
                LOG.warning("%s/%s/%s: %d usable columns, skipped", family, variant, module, n)
                continue
            observed = float(values[in_module].mean())
            null = np.array([rng.choice(pool, n, replace=False).mean()
                             for _ in range(n_perm)])
            rows.append({
                "family": family, "variant": variant, "module": module,
                "n_columns": n,
                "mean_delta": round(observed, 4),
                "mean_delta_other_columns": round(float(values[usable & ~in_module].mean()), 4),
                "null_mean": round(float(null.mean()), 4),
                "null_sd": round(float(null.std()), 4),
                "z": round(float((observed - null.mean()) / null.std()), 3),
                "p_perm": round(float((1 + int((null >= observed).sum())) / (n_perm + 1)), 5),
            })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # Correct within each variant: the three variants are the same hypothesis
    # asked three ways, not three independent families of tests.
    out["q_bh"] = np.nan
    for variant, grp in out.groupby("variant"):
        out.loc[grp.index, "q_bh"] = cons.benjamini_hochberg(grp.p_perm.to_numpy())
    out["q_bh"] = out.q_bh.round(5)
    out["significant_q05"] = out.q_bh < 0.05
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    parser.add_argument("--n-perm", type=int, default=N_PERM,
                        help="Column permutations for the module-deficit test.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())
    provenance = pd.read_csv(args.outdir / "column_provenance.csv")
    rng_alphabet = dcalib.ALPHABET

    rows, summaries = [], []
    for family, paralogues in FAMILY_PARALOGUES.items():
        aln = dcalib.read_fasta(args.outdir / "alignments" / f"{family}_family.fasta")
        msa = dcalib.encode_msa(aln)
        prov = provenance[provenance.family == family].set_index("column")

        scopes: dict[str, list[int]] = {p: group_rows(aln, p) for p in paralogues}
        scopes["FAMILY"] = [i for i, k in enumerate(aln)
                            if not k.startswith(("REF|", "STRUCT|"))]

        for scope, idx in scopes.items():
            if len(idx) < 20:
                LOG.warning("%s/%s: only %d sequences, skipped", family, scope, len(idx))
                continue
            sub = msa[idx]
            w_hen = cons.henikoff_weights(sub)
            w_dca = dcalib.sequence_weights(sub, 0.8)
            res = cons.column_constraint(sub, w_hen, rng_alphabet)
            lo, hi = cons.bootstrap_constraint(sub, w_hen, rng_alphabet,
                                               n_boot=args.n_boot, seed=SEED, measure="jsd")
            LOG.info("%-6s %-7s n=%4d  kish_neff=%6.1f  neff_dca80=%5.1f  "
                     "mean JSD=%.3f", family, scope, len(idx),
                     cons.kish_effective_size(w_hen), w_dca.sum(),
                     float(np.mean(res.jsd)))

            # Kish n, not sum(w): Henikoff weights are normalised to sum to N,
            # so their sum restates n_sequences and says nothing about depth.
            kish = cons.kish_effective_size(w_hen)
            summaries.append({"family": family, "scope": scope, "n_sequences": len(idx),
                              "neff_dca80": round(float(w_dca.sum()), 1),
                              "mean_jsd": round(float(np.mean(res.jsd)), 4),
                              "mean_conservation": round(float(np.mean(res.conservation_norm)), 4),
                              "kish_neff": round(kish, 1),
                              "kish_frac_of_N": round(kish / len(idx), 3)})

            # Canonical numbering for whichever reference this scope corresponds to
            ref_for_scope = scope if scope in maps[family]["references"] else None
            col_to_canon = ({int(k): v for k, v in
                             maps[family]["references"][ref_for_scope]["column_to_canonical"].items()}
                            if ref_for_scope else {})

            for col in range(msa.shape[1]):
                p = prov.loc[col] if col in prov.index else None
                rows.append({
                    "family": family, "scope": scope, "column": col,
                    "canonical_residue": col_to_canon.get(col, ""),
                    "residue_letter": (p[f"{scope}_letter"] if p is not None
                                       and f"{scope}_letter" in p.index else ""),
                    "feature": (p["features"] if p is not None else ""),
                    "entropy_bits": round(float(res.entropy_bits[col]), 4),
                    "conservation_norm": round(float(res.conservation_norm[col]), 4),
                    "jsd": round(float(res.jsd[col]), 4),
                    "jsd_ci_low": round(float(lo[col]), 4),
                    "jsd_ci_high": round(float(hi[col]), 4),
                    "gap_fraction": round(float(res.gap_fraction[col]), 4),
                })

    constraint_df = pd.DataFrame(rows)
    constraint_df.to_csv(args.outdir / "residue_constraint.csv", index=False)
    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(args.outdir / "constraint_scope_summary.csv", index=False)

    # Module-level aggregates: the level at which the effective sample supports
    # a claim, since a single column's CI is wide at this depth.
    mod_rows = []
    for (family, scope, feature), g in constraint_df.assign(
            feature=constraint_df.feature.replace("", np.nan)).dropna(subset=["feature"]) \
            .groupby(["family", "scope", "feature"]):
        boot = np.random.default_rng(SEED).choice(g.jsd.values, (2000, len(g)), replace=True)
        means = boot.mean(axis=1)
        mod_rows.append({"family": family, "scope": scope, "feature": feature,
                         "n_columns": len(g),
                         "mean_jsd": round(float(g.jsd.mean()), 4),
                         "ci_low": round(float(np.quantile(means, 0.025)), 4),
                         "ci_high": round(float(np.quantile(means, 0.975)), 4),
                         "mean_conservation": round(float(g.conservation_norm.mean()), 4),
                         "mean_gap_fraction": round(float(g.gap_fraction.mean()), 4)})
    module_df = pd.DataFrame(mod_rows).sort_values(["family", "scope", "feature"])
    module_df.to_csv(args.outdir / "module_constraint_summary.csv", index=False)

    # Is the contrast paralogue selectively less constrained in the inhibitory
    # modules, or just less constrained throughout? SOCS2 is used here only as a
    # constraint comparator; it is absent from every coupling design.
    deficit_df = module_deficit_test(
        constraint_df, family="socs",
        binders=("SOCS1", "SOCS3"), contrast="SOCS2",
        modules=("KIR", "ESS", "SH2", "SOCS_box"),
        n_perm=args.n_perm, seed=SEED,
    )
    deficit_df.to_csv(args.outdir / "module_deficit_test.csv", index=False)

    print(summary_df.to_string(index=False))
    print()
    print(module_df.to_string(index=False))
    print()
    print(deficit_df.to_string(index=False))


if __name__ == "__main__":
    main()
