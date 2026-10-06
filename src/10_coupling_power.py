"""Step 11 - DCA infeasibility: formal power analysis.

Turns the step-2 scoping finding into an empirical result. Rather than asserting
that residue-level inter-protein DCA cannot work at this depth, this runs it and
shows that what comes out is indistinguishable from what comes out of alignments
in which the JAK-SOCS pairing has been destroyed.

Method:
  - Build the best available paired alignment (the factorial pool: all six
    JAK x SOCS paralogue combinations by species, the configuration that scored
    highest on the depth gate).
  - Compute inter-block coupling with mutual information + average product
    correction, applied to the inter-block submatrix only. MI rather than
    plmDCA: at N_eff < 130 against ~40,000 column pairs a pseudolikelihood fit
    is determined by its prior, and MI is the estimator the co-primary decision
    in METHODS_REVIEW.md selects for this regime.
  - Repeat on N_REP species-scrambled replicates, which preserve both marginal
    alignments, all conservation and gap structure, and the phylogeny of each
    family, while destroying the pairing.
  - Compare the observed inter-block score distribution and its top tail against
    the scrambled ensemble.

If the observed and scrambled distributions are not separable, the couplings
carry no pairing information, which is the claim the power analysis needs.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib
import lib_constraint as cons

LOG = logging.getLogger("step11")
SEED = 42
N_REP = 30

JAK_PARALOGUES = ["JAK1", "JAK2", "TYK2"]
SOCS_PARALOGUES = ["SOCS1", "SOCS3"]


def build_pool(jak_aln: dict[str, str], socs_aln: dict[str, str],
               jak_cols: list[int], socs_cols: list[int],
               pairing: dict[str, str] | None = None) -> tuple[np.ndarray, list[str]]:
    """Concatenate JAK and SOCS blocks over all six paralogue combinations.

    ``pairing`` optionally remaps each JAK record's species to a different SOCS
    species, which is how the scrambled null is produced.
    """
    jk = {(k.split("|")[0], k.split("|")[1]): k for k in jak_aln}
    sk = {(k.split("|")[0], k.split("|")[1]): k for k in socs_aln}
    rows, labels = [], []
    for jp, sp in itertools.product(JAK_PARALOGUES, SOCS_PARALOGUES):
        jt = {t for (p, t) in jk if p == jp}
        st = {t for (p, t) in sk if p == sp}
        for t in sorted(jt & st):
            t_socs = pairing.get(t, t) if pairing else t
            if (sp, t_socs) not in sk:
                continue
            jrow = jak_aln[jk[(jp, t)]]
            srow = socs_aln[sk[(sp, t_socs)]]
            rows.append("".join(jrow[c] for c in jak_cols)
                        + "".join(srow[c] for c in socs_cols))
            labels.append(f"{jp}-{sp}|{t}")
    msa = dcalib.encode_msa(dict(zip(range(len(rows)), rows)))
    return msa, labels


def inter_block_mi_apc(msa: np.ndarray, weights: np.ndarray, n_jak: int
                       ) -> np.ndarray:
    """MI with APC for every JAK-column x SOCS-column pair.

    APC is computed on the inter-block submatrix alone. Applying it to the full
    two-block matrix would subtract a background that is partly the block
    structure itself, which systematically distorts exactly these scores.
    """
    q = cons.N_AA + 1
    n_socs = msa.shape[1] - n_jak
    w = weights / weights.sum()

    # Per-column weighted marginals and entropies
    def marginals(block: np.ndarray) -> np.ndarray:
        out = np.zeros((block.shape[1], q))
        for s in range(q):
            out[:, s] = ((block == s) * w[:, None]).sum(axis=0)
        return out

    pj = marginals(msa[:, :n_jak])
    ps = marginals(msa[:, n_jak:])

    def ent(p: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(p > 0, p * np.log2(p), 0.0)
        return -t.sum(axis=-1)

    hj, hs = ent(pj), ent(ps)

    mi = np.zeros((n_jak, n_socs))
    onehot_s = np.stack([(msa[:, n_jak:] == s) for s in range(q)], axis=-1)  # (N, Ls, q)
    for i in range(n_jak):
        col = msa[:, i]
        # joint[s_j, l, s_s] via weighted contraction
        acc = np.zeros((q, n_socs, q))
        for sj in range(q):
            m = (col == sj)
            if not m.any():
                continue
            acc[sj] = np.tensordot(w[m], onehot_s[m], axes=(0, 0))
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(acc > 0, acc * np.log2(acc), 0.0)
        h_joint = -t.sum(axis=(0, 2))
        mi[i] = np.maximum(hj[i] + hs - h_joint, 0.0)

    # Average product correction on the inter-block submatrix
    row_mean = mi.mean(axis=1, keepdims=True)
    col_mean = mi.mean(axis=0, keepdims=True)
    grand = mi.mean()
    return mi - (row_mean * col_mean) / grand if grand > 0 else mi


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--n-rep", type=int, default=N_REP)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    started = time.time()

    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())
    provenance = pd.read_csv(args.outdir / "column_provenance.csv")
    jak_aln = dcalib.read_fasta(args.outdir / "alignments" / "jak_family.fasta")
    socs_aln = dcalib.read_fasta(args.outdir / "alignments" / "socs_family.fasta")

    # Restrict to the annotated functional columns: JH1 on the JAK side, the
    # KIR-to-SOCS-box span on the SOCS side. This is the most favourable
    # configuration for the gate, not the most convenient one.
    pj = provenance[provenance.family == "jak"]
    ps = provenance[provenance.family == "socs"]
    jak_cols = sorted(pj.loc[pj.features.fillna("") == "JH1", "column"])
    socs_cols = sorted(ps.loc[ps.features.fillna("") != "", "column"])
    LOG.info("blocks: JAK %d columns, SOCS %d columns", len(jak_cols), len(socs_cols))

    msa, labels = build_pool(jak_aln, socs_aln, jak_cols, socs_cols)
    weights = dcalib.sequence_weights(msa, 0.8)
    neff = float(weights.sum())
    l1, l2 = len(jak_cols), len(socs_cols)
    gate = (l1 + l2) / 2.0
    LOG.info("pool: N=%d  N_eff=%.1f  L=%d  gate needs %.1f  ratio %.3f",
             msa.shape[0], neff, l1 + l2, gate, neff / gate)

    obs = inter_block_mi_apc(msa, weights, l1)
    LOG.info("observed inter-block scores computed (%d pairs)", obs.size)

    rng = np.random.default_rng(SEED)
    socs_species = sorted({k.split("|")[1] for k in socs_aln if "|" in k})
    null_means, null_tops, null_maxes = [], [], []
    top_k = 100
    for rep in range(args.n_rep):
        shuffled = list(socs_species)
        rng.shuffle(shuffled)
        pairing = dict(zip(socs_species, shuffled))
        m_null, _ = build_pool(jak_aln, socs_aln, jak_cols, socs_cols, pairing)
        w_null = dcalib.sequence_weights(m_null, 0.8)
        sc = inter_block_mi_apc(m_null, w_null, l1)
        null_means.append(float(sc.mean()))
        null_tops.append(float(np.sort(sc.ravel())[-top_k:].mean()))
        null_maxes.append(float(sc.max()))
        LOG.info("  scramble %2d/%d: mean %.5f  top%d %.5f  max %.5f",
                 rep + 1, args.n_rep, null_means[-1], top_k, null_tops[-1], null_maxes[-1])

    obs_mean = float(obs.mean())
    obs_top = float(np.sort(obs.ravel())[-top_k:].mean())
    obs_max = float(obs.max())
    nm, nt = np.array(null_means), np.array(null_tops)

    def z_and_p(o: float, null: np.ndarray) -> tuple[float, float]:
        z = (o - null.mean()) / null.std() if null.std() > 0 else np.nan
        p = ((np.abs(null - null.mean()) >= abs(o - null.mean())).sum() + 1) / (len(null) + 1)
        return float(z), float(p)

    z_mean, p_mean = z_and_p(obs_mean, nm)
    z_top, p_top = z_and_p(obs_top, nt)

    # Distributional comparison against the pooled scrambled scores
    rows = [
        {"statistic": "mean_inter_block_score", "observed": round(obs_mean, 6),
         "null_mean": round(float(nm.mean()), 6), "null_sd": round(float(nm.std()), 6),
         "z": round(z_mean, 3), "p_empirical": p_mean},
        {"statistic": f"mean_of_top_{top_k}", "observed": round(obs_top, 6),
         "null_mean": round(float(nt.mean()), 6), "null_sd": round(float(nt.std()), 6),
         "z": round(z_top, 3), "p_empirical": p_top},
        {"statistic": "max_inter_block_score", "observed": round(obs_max, 6),
         "null_mean": round(float(np.mean(null_maxes)), 6),
         "null_sd": round(float(np.std(null_maxes)), 6),
         "z": round(z_and_p(obs_max, np.array(null_maxes))[0], 3),
         "p_empirical": z_and_p(obs_max, np.array(null_maxes))[1]},
    ]
    res = pd.DataFrame(rows)
    res.to_csv(args.outdir / "dca_null_comparison.csv", index=False)

    depth = pd.DataFrame([{
        "dataset": "factorial_pool_JH1_x_SOCSmodules", "n_records": int(msa.shape[0]),
        "l_jak": l1, "l_socs": l2, "l_total": l1 + l2,
        "neff_theta_0.7": round(float(dcalib.sequence_weights(msa, 0.7).sum()), 1),
        "neff_theta_0.8": round(neff, 1),
        "neff_theta_0.9": round(float(dcalib.sequence_weights(msa, 0.9).sum()), 1),
        "neff_over_l": round(neff / (l1 + l2), 4),
        "ovchinnikov_threshold": gate, "gate_ratio": round(neff / gate, 4),
        "passes_gate": bool(neff > gate),
    }])
    depth.to_csv(args.outdir / "depth_power_report.csv", index=False)

    np.savez_compressed(args.outdir / "dca_null_distributions.npz",
                        observed_scores=obs.astype(np.float32),
                        null_means=nm, null_tops=nt, null_maxes=np.array(null_maxes))
    LOG.info("done in %.1f s", time.time() - started)
    print(depth.to_string(index=False))
    print()
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
