"""Per-column evolutionary constraint and specificity statistics.

Two families of measure, matched to the two kinds of signal the JAK/SOCS
alignment actually contains:

**Constraint** (within a paralogue group): how strongly purifying selection
holds a position fixed. Sequence-weighted Shannon entropy and the
Jensen-Shannon divergence against an amino-acid background, following Capra &
Singh (2007), which benchmarks as the strongest single-column conservation
measure.

**Specificity** (between paralogue groups): how much a position's residue
identity tells you which paralogue you are looking at. Mutual information
between column state and paralogue label, with a Miller-Madow bias correction
and a label-permutation null - both necessary because MI is upward-biased at
finite sample size and the bias grows with the number of occupied bins.

All measures take sequence weights, so redundancy in the alignment does not
inflate apparent conservation.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

LOG = logging.getLogger("constraint")

#: Index 0 is the gap symbol; residues occupy 1..20. Matches dcalib.ALPHABET.
N_AA = 20
GAP = 0

#: Background amino-acid frequencies (Robinson & Robinson 1991), the
#: distribution the Jensen-Shannon conservation score is measured against.
#: Ordered to match dcalib.ALPHABET[1:].
_BG = {
    "A": 0.0787, "C": 0.0157, "D": 0.0531, "E": 0.0636, "F": 0.0407,
    "G": 0.0690, "H": 0.0227, "I": 0.0591, "K": 0.0594, "L": 0.0964,
    "M": 0.0238, "N": 0.0427, "P": 0.0468, "Q": 0.0393, "R": 0.0526,
    "S": 0.0694, "T": 0.0550, "V": 0.0667, "W": 0.0118, "Y": 0.0311,
}


def background_frequencies(alphabet: str) -> np.ndarray:
    """Background distribution over the 20 residues, ordered by ``alphabet``."""
    q = np.array([_BG[c] for c in alphabet[1:]], dtype=np.float64)
    return q / q.sum()


def henikoff_weights(msa: np.ndarray, include_gaps: bool = False) -> np.ndarray:
    """Henikoff & Henikoff (1994) position-based sequence weights.

    Why not the DCA convention: identity-threshold weighting at 80% is designed
    to stop near-duplicate sequences dominating a coupling fit, and within a set
    of orthologues of a single paralogue almost every pair exceeds that
    threshold - JAK2 orthologues collapse to N_eff = 8.8 out of 625 sequences.
    Using that weighting for conservation would discard nearly all the data.

    Position-based weighting instead assigns each sequence, at each column, a
    share 1/(r_i * n_i(a)) where r_i is the number of residue types in the
    column and n_i(a) the count of that sequence's residue. Rare residues earn
    more weight, redundancy is still discounted, and the effective sample stays
    close to the number of distinct sequences. This is the convention used by
    conservation scores including Capra & Singh's JSD.
    """
    n, length = msa.shape
    weights = np.zeros(n, dtype=np.float64)
    lo = 0 if include_gaps else 1
    for j in range(length):
        col = msa[:, j]
        symbols, counts = np.unique(col, return_counts=True)
        keep = symbols >= lo
        symbols, counts = symbols[keep], counts[keep]
        if len(symbols) <= 1:
            continue                      # invariant column carries no weight info
        share = {s: 1.0 / (len(symbols) * c) for s, c in zip(symbols, counts)}
        weights += np.array([share.get(a, 0.0) for a in col])
    total = weights.sum()
    if total <= 0:                        # degenerate alignment: fall back to uniform
        return np.ones(n, dtype=np.float64)
    return weights * (n / total)          # normalised so weights sum to N


def weighted_counts(msa: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """(L, 21) weighted symbol counts, gap in column 0."""
    n, length = msa.shape
    counts = np.zeros((length, N_AA + 1), dtype=np.float64)
    for symbol in range(N_AA + 1):
        counts[:, symbol] = ((msa == symbol) * weights[:, None]).sum(axis=0)
    return counts


def residue_frequencies(counts: np.ndarray, pseudocount: float = 0.0
                        ) -> tuple[np.ndarray, np.ndarray]:
    """Return (residue frequencies over 20 AAs, per-column gap fraction).

    Frequencies are renormalised over non-gap symbols, so a position is scored
    on the residues actually present; the gap fraction is returned separately
    and used to downweight the score rather than being folded into it.
    """
    total = counts.sum(axis=1, keepdims=True)
    gap_frac = (counts[:, GAP:GAP + 1] / np.maximum(total, 1e-12)).ravel()
    res = counts[:, 1:] + pseudocount
    res_total = res.sum(axis=1, keepdims=True)
    freqs = res / np.maximum(res_total, 1e-12)
    return freqs, gap_frac


def shannon_entropy(freqs: np.ndarray) -> np.ndarray:
    """Per-column Shannon entropy in bits."""
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(freqs > 0, freqs * np.log2(freqs), 0.0)
    return -terms.sum(axis=1)


def jsd_conservation(freqs: np.ndarray, gap_frac: np.ndarray,
                     background: np.ndarray, gap_penalty: bool = True) -> np.ndarray:
    """Jensen-Shannon divergence of each column against ``background``.

    Capra & Singh (2007), Bioinformatics 23:1875. Returned in bits on [0, 1];
    1 means a perfectly conserved column whose residue is rare in the
    background. With ``gap_penalty`` the score is scaled by the non-gap
    fraction, which is the published convention.
    """
    m = 0.5 * (freqs + background[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_p = np.where(freqs > 0, freqs * np.log2(freqs / m), 0.0).sum(axis=1)
        kl_q = np.where(background[None, :] > 0,
                        background[None, :] * np.log2(background[None, :] / m), 0.0).sum(axis=1)
    jsd = 0.5 * (kl_p + kl_q)
    return jsd * (1.0 - gap_frac) if gap_penalty else jsd


@dataclass
class ConstraintResult:
    """Per-column constraint scores for one alignment (or paralogue group)."""

    entropy_bits: np.ndarray
    conservation_norm: np.ndarray     # 1 - H/log2(20), in [0, 1]
    jsd: np.ndarray
    gap_fraction: np.ndarray
    n_eff: float
    n_sequences: int


def column_constraint(msa: np.ndarray, weights: np.ndarray,
                      alphabet: str, pseudocount: float = 0.0) -> ConstraintResult:
    """Compute all constraint measures for every column of ``msa``."""
    counts = weighted_counts(msa, weights)
    freqs, gap_frac = residue_frequencies(counts, pseudocount)
    entropy = shannon_entropy(freqs)
    return ConstraintResult(
        entropy_bits=entropy,
        conservation_norm=1.0 - entropy / np.log2(N_AA),
        jsd=jsd_conservation(freqs, gap_frac, background_frequencies(alphabet)),
        gap_fraction=gap_frac,
        n_eff=float(weights.sum()),
        n_sequences=int(msa.shape[0]),
    )


def bootstrap_constraint(msa: np.ndarray, weights: np.ndarray, alphabet: str,
                         n_boot: int = 1000, seed: int = 42,
                         measure: str = "jsd") -> tuple[np.ndarray, np.ndarray]:
    """Bootstrap 95% CI for a per-column constraint measure.

    Resamples sequences with replacement, carrying their weights, so the
    interval reflects uncertainty from the finite species sample rather than
    from the alignment length.
    """
    rng = np.random.default_rng(seed)
    n = msa.shape[0]
    draws = np.empty((n_boot, msa.shape[1]), dtype=np.float64)
    bg = background_frequencies(alphabet)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        counts = weighted_counts(msa[idx], weights[idx])
        freqs, gap_frac = residue_frequencies(counts)
        draws[b] = (jsd_conservation(freqs, gap_frac, bg) if measure == "jsd"
                    else 1.0 - shannon_entropy(freqs) / np.log2(N_AA))
    return (np.quantile(draws, 0.025, axis=0), np.quantile(draws, 0.975, axis=0))


# ------------------------------------------------------- specificity (SDP)

def _entropy_mm(counts: np.ndarray, n: float) -> float:
    """Miller-Madow bias-corrected Shannon entropy of a count vector.

    The plug-in entropy estimator is downward-biased by roughly (m-1)/(2N),
    where m is the number of occupied bins. Correcting it matters here because
    columns differ widely in how many residues they use, so the raw bias would
    masquerade as a specificity signal.
    """
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts / total
    nz = p > 0
    h = -(p[nz] * np.log2(p[nz])).sum()
    return h + (nz.sum() - 1) / (2.0 * max(n, 1.0) * np.log(2))


def _entropy_mm_vec(p: np.ndarray, n_eff: float, axis: tuple[int, ...] | int) -> np.ndarray:
    """Miller-Madow corrected entropy of probability arrays, vectorised."""
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(p > 0, p * np.log2(p), 0.0)
    h = -terms.sum(axis=axis)
    occupied = (p > 0).sum(axis=axis)
    return h + (occupied - 1) / (2.0 * max(n_eff, 1.0) * np.log(2))


def _onehot_matrix(msa: np.ndarray):
    """Sparse (N, L*(Q)) indicator of 'sequence i has symbol s at column j'.

    Lets each permutation's joint contingency table be obtained with one sparse
    matrix-vector product per group - O(N*L) instead of O(N*L*Q) - which is what
    makes a 10,000-permutation null affordable.
    """
    from scipy import sparse

    n, length = msa.shape
    rows = np.repeat(np.arange(n), length)
    cols = (np.tile(np.arange(length), n) * (N_AA + 1) + msa.ravel().astype(np.int64))
    data = np.ones(n * length, dtype=np.float64)
    return sparse.csr_matrix((data, (rows, cols)),
                             shape=(n, length * (N_AA + 1))).T.tocsr()


def specificity_mi(msa: np.ndarray, weights: np.ndarray, labels: np.ndarray,
                   n_perm: int = 10000, seed: int = 42,
                   keep_null: bool = False) -> dict[str, np.ndarray]:
    """Mutual information between column residue state and group label.

    Returns Miller-Madow corrected MI per column plus a label-permutation null,
    giving a z-score and an empirical p-value. Permuting labels leaves each
    column's residue composition untouched and only breaks its association with
    the group, so the null tests exactly the association of interest.

    Gaps are retained as a symbol: a module that is deleted in one paralogue and
    present in another is a genuine specificity signal, not missing data.
    """
    n, length = msa.shape
    groups = np.unique(labels)
    q = N_AA + 1
    n_eff = float(weights.sum())
    onehot_t = _onehot_matrix(msa)          # (L*q, N)

    def mi_for(lab: np.ndarray) -> np.ndarray:
        joint = np.empty((length, len(groups), q), dtype=np.float64)
        for gi, g in enumerate(groups):
            v = weights * (lab == g)
            joint[:, gi, :] = (onehot_t @ v).reshape(length, q)
        total = joint.sum(axis=(1, 2), keepdims=True)
        total = np.maximum(total, 1e-12)
        p_joint = joint / total
        p_x = p_joint.sum(axis=1)            # (L, q)
        p_y = p_joint.sum(axis=2)            # (L, G)
        h_x = _entropy_mm_vec(p_x, n_eff, axis=1)
        h_y = _entropy_mm_vec(p_y, n_eff, axis=1)
        h_xy = _entropy_mm_vec(p_joint, n_eff, axis=(1, 2))
        return np.maximum(h_x + h_y - h_xy, 0.0)

    observed = mi_for(labels)
    rng = np.random.default_rng(seed)
    ge = np.zeros(length, dtype=np.int64)
    running_sum = np.zeros(length, dtype=np.float64)
    running_sq = np.zeros(length, dtype=np.float64)
    null_max = np.full(length, -np.inf)
    null_all = np.empty((n_perm, length), dtype=np.float32) if keep_null else None

    for p in range(n_perm):
        draw = mi_for(rng.permutation(labels))
        ge += draw >= observed
        running_sum += draw
        running_sq += draw ** 2
        null_max = np.maximum(null_max, draw)
        if keep_null:
            null_all[p] = draw
        if n_perm >= 500 and (p + 1) % max(n_perm // 5, 1) == 0:
            LOG.info("  permutations %d/%d", p + 1, n_perm)

    mu = running_sum / n_perm
    var = np.maximum(running_sq / n_perm - mu ** 2, 0.0)
    sd = np.sqrt(var)
    z = np.divide(observed - mu, sd, out=np.zeros_like(observed), where=sd > 0)
    out = {"mi": observed, "null_mean": mu, "null_sd": sd, "z": z,
           "p_empirical": (ge + 1) / (n_perm + 1), "null_max": null_max}
    if keep_null:
        out["null"] = null_all
    return out


def group_contrast(msa: np.ndarray, weights: np.ndarray, labels: np.ndarray,
                   alphabet: str) -> dict[str, np.ndarray]:
    """Within-group conserved / between-group divergent contrast per column.

    A classical specificity-determining-position signature: high conservation
    inside each paralogue group combined with different consensus residues
    between groups. Reported alongside MI because the two disagree in
    interpretable ways - MI also fires on columns that are variable in one
    group and fixed in another.
    """
    groups = np.unique(labels)
    bg = background_frequencies(alphabet)
    within, consensus = [], []
    for g in groups:
        sel = labels == g
        counts = weighted_counts(msa[sel], weights[sel])
        freqs, gap = residue_frequencies(counts)
        within.append(jsd_conservation(freqs, gap, bg))
        consensus.append(np.argmax(freqs, axis=1))
    within = np.array(within)
    consensus = np.array(consensus)
    n_distinct = np.array([len(set(consensus[:, j])) for j in range(msa.shape[1])])
    return {"within_group_conservation": within.mean(axis=0),
            "min_within_group_conservation": within.min(axis=0),
            "n_distinct_consensus": n_distinct,
            "consensus_per_group": consensus,
            "groups": groups,
            "contrast": within.min(axis=0) * (n_distinct - 1)}


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """BH-adjusted q-values for a vector of p-values."""
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.clip(ranked, 0, 1)
    return q
