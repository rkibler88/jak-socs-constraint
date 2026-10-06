"""Step 10 - Declared motif tests.

Tests the motif hypotheses fixed in advance in METHODS_REVIEW.md section 9, as a
single family under one Benjamini-Hochberg correction. The set was declared from
the literature before any constraint or specificity number was computed, so
these are confirmatory tests rather than a selection over the 235-column table.

Declared sources:
  Zhai et al. 2026 (10.1039/d6cp00092d)   KIR 'Y', 'QR', 'FF' motifs; JAK 'GQM'
  Yasukawa et al. 1999 (10.1093/emboj/18.5.1309)  SOCS1 KIR and ESS; I68, L75
  Bergamin et al. 2006 (10.1016/j.str.2006.06.011)  SH2 phosphotyrosine pocket
  Kershaw 2013 / Liau 2018                JAK substrate groove, activation loop

Each motif is tested two ways: mean constraint against size-matched non-motif
controls from the same domain, and mean binder-vs-SOCS2 specificity against the
same controls.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib
import lib_constraint as cons

LOG = logging.getLogger("step10")
SEED = 42
N_PERM = 10000


@dataclass(frozen=True)
class Motif:
    """A pre-declared motif, addressed in one reference's canonical numbering."""

    name: str
    family: str
    reference: str
    residues: tuple[int, ...]
    control_domain: str | None    # restrict controls to this feature, if any
    source: str


#: Motif positions were LOCATED IN THE SEQUENCES, not assumed. Zhai et al. place
#: their motifs in "the KIR region and BC loop"; the BC loop lies in the SH2
#: domain, so the 'QR' and 'FF' motifs are not KIR residues. Literal string
#: search of the human canonical sequences gives:
#:   'GQM'  JAK1 1097, JAK2 1071, TYK2 1117   (one occurrence each, inside JH1)
#:   'FF'   SOCS1 112, SOCS3 79               (SH2 / BC loop); absent from SOCS2
#:   'QR'   SOCS1 108 only                    (not positionally homologous in SOCS3)
#:   'FLVR' SOCS1 101 -> invariant SH2 arginine R104; absent as a literal string
#:          from SOCS3 and SOCS2, so the pocket arginine is addressed via SOCS1.
MOTIFS: tuple[Motif, ...] = (
    # --- Zhai 2026: the motifs attributed to SOCS1/SOCS3 selectivity for JAK1/TYK2
    Motif("KIR_aromatic", "socs", "SOCS1", (64,), "KIR", "Zhai 2026 ('Y' motif)"),
    Motif("BCloop_FF", "socs", "SOCS1", (112, 113), "SH2", "Zhai 2026 ('FF' motif)"),
    Motif("BCloop_QR_SOCS1only", "socs", "SOCS1", (108, 109), "SH2",
          "Zhai 2026 ('QR' motif; not homologous in SOCS3)"),
    # --- Yasukawa 1999: KIR H54-R59 implicated by mutagenesis; ESS I68 and L75
    #     conserved among JAB-related proteins, required for activation-loop binding
    Motif("KIR_H54_R59", "socs", "SOCS1", tuple(range(54, 60)), "KIR", "Yasukawa 1999"),
    Motif("ESS_I68_L75", "socs", "SOCS1", (68, 75), "ESS", "Yasukawa 1999"),
    Motif("KIR_full", "socs", "SOCS1", tuple(range(55, 67)), None, "UniProt / Yasukawa 1999"),
    Motif("ESS_full", "socs", "SOCS1", tuple(range(67, 79)), None, "UniProt / Yasukawa 1999"),
    # --- Bergamin 2006: SH2 phosphotyrosine pocket. Internal negative control -
    #     shared by every SOCS protein, so it should NOT discriminate binders.
    Motif("SH2_pTyr_arginine", "socs", "SOCS1", (104,), "SH2", "Bergamin 2006"),
    Motif("SH2_FLVR", "socs", "SOCS1", (101, 102, 103, 104), "SH2", "Bergamin 2006"),
    # --- JAK side, JAK2 numbering
    Motif("JAK_GQM_motif", "jak", "JAK2", (1071, 1072, 1073), "JH1", "Zhai 2026"),
    Motif("JAK_activation_loop", "jak", "JAK2", tuple(range(1004, 1015)), "JH1",
          "Kershaw 2013 / Liau 2018"),
)

#: Aromatic residues, for the aromaticity measure. Zhai et al.'s mechanism is
#: Met-aromatic interaction with the JAK GQM motif, so aromaticity rather than
#: exact residue identity is the property under selection.
AROMATIC = set("FYW")

#: Constraint scope each motif's test is run in.
SCOPES = {"socs": ["SOCS1", "SOCS3", "SOCS2"], "jak": ["JAK1", "JAK2", "TYK2"]}


def matched_control_columns(all_cols: np.ndarray, motif_cols: set[int],
                            domain_cols: set[int] | None) -> np.ndarray:
    pool = domain_cols if domain_cols else set(all_cols)
    return np.array(sorted(pool - motif_cols))


def permutation_test(values: pd.Series, motif_cols: list[int], control_cols: np.ndarray,
                     n_perm: int, seed: int) -> tuple[float, float, float]:
    """Mean motif value vs size-matched random draws from the control pool."""
    rng = np.random.default_rng(seed)
    v = values.to_dict()
    obs = float(np.mean([v[c] for c in motif_cols if c in v]))
    pool = np.array([v[c] for c in control_cols if c in v])
    if len(pool) < len(motif_cols) + 2:
        return obs, np.nan, np.nan
    k = len([c for c in motif_cols if c in v])
    draws = np.array([rng.choice(pool, k, replace=False).mean() for _ in range(n_perm)])
    p = ((np.abs(draws - pool.mean()) >= abs(obs - pool.mean())).sum() + 1) / (n_perm + 1)
    z = (obs - draws.mean()) / draws.std() if draws.std() > 0 else np.nan
    return obs, float(p), float(z)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    parser.add_argument("--n-perm", type=int, default=N_PERM)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())
    provenance = pd.read_csv(args.outdir / "column_provenance.csv")
    constraint = pd.read_csv(args.outdir / "residue_constraint.csv")
    sdp = pd.read_csv(args.outdir / "sdp_table.csv")
    iface = pd.read_csv(args.outdir / "interface_residues.csv")

    # interface columns per family, pooled over both structures
    STRUCT_MAP = {("4GL9", "SOCS"): ("socs", "SOCS3_MOUSE_4GL9"),
                  ("4GL9", "JAK"): ("jak", "JAK2_MOUSE_4GL9"),
                  ("6C7Y", "SOCS"): ("socs", "SOCS1_CHICK_6C7Y"),
                  ("6C7Y", "JAK"): ("jak", "JAK1")}
    iface_cols: dict[str, set[int]] = {"socs": set(), "jak": set()}
    for (pdb, side), (fam, mapname) in STRUCT_MAP.items():
        c2c = {int(k): v for k, v in maps[fam]["references"][mapname]["canonical_to_column"].items()}
        sel = iface[(iface.structure == pdb) & (iface.side == side)]
        iface_cols[fam] |= {c2c[r] for r in sel.canonical_residue if r in c2c}

    rows = []
    for m in MOTIFS:
        c2col = {int(k): v for k, v in
                 maps[m.family]["references"][m.reference]["canonical_to_column"].items()}
        motif_cols = [c2col[r] for r in m.residues if r in c2col]
        missing = [r for r in m.residues if r not in c2col]
        if not motif_cols:
            LOG.warning("%s: no residues mapped, skipped", m.name)
            continue

        prov = provenance[provenance.family == m.family]
        all_cols = prov.column.values
        domain_cols = None
        if m.control_domain:
            domain_cols = set(prov.loc[prov.features.fillna("").str.split(";")
                                       .apply(lambda x: m.control_domain in x), "column"])
        ctrl = matched_control_columns(all_cols, set(motif_cols), domain_cols)
        # A 6-residue motif inside a 12-column domain leaves too few controls for
        # a resampling test. Fall back to the whole block and record which pool
        # was used, rather than silently returning no result.
        pool_used = m.control_domain or "whole_block"
        if len(ctrl) < len(motif_cols) + 2:
            ctrl = matched_control_columns(all_cols, set(motif_cols), None)
            pool_used = f"whole_block (fallback from {m.control_domain})"

        base = {"motif": m.name, "family": m.family, "reference": m.reference,
                "residues": ",".join(map(str, m.residues)),
                "n_residues_mapped": len(motif_cols),
                "residues_unmapped": ",".join(map(str, missing)),
                "control_domain": pool_used,
                "n_control_columns": len(ctrl), "source": m.source,
                "interface_overlap": sum(c in iface_cols[m.family] for c in motif_cols)}

        # constraint, per paralogue scope
        for scope in SCOPES[m.family]:
            cs = constraint[(constraint.family == m.family) & (constraint.scope == scope)]
            vals = cs.set_index("column")["jsd"]
            obs, p, z = permutation_test(vals, motif_cols, ctrl, args.n_perm, SEED)
            ctrl_mean = float(vals.reindex(ctrl).dropna().mean())
            rows.append({**base, "measure": "constraint_jsd", "scope": scope,
                         "observed": round(obs, 4), "control_mean": round(ctrl_mean, 4),
                         "fold": round(obs / ctrl_mean, 3) if ctrl_mean else np.nan,
                         "z": round(z, 3) if z == z else np.nan, "p_perm": p})

        # Aromaticity per paralogue at the motif columns. Reported as a plain
        # weighted fraction with a Wilson interval rather than a permutation
        # test: when the contrast is 0.99 vs 0.005 a p-value adds nothing, and
        # the effect size is the result.
        aln = dcalib.read_fasta(args.outdir / "alignments" / f"{m.family}_family.fasta")
        msa = dcalib.encode_msa(aln)
        akeys = list(aln)
        arom_idx = {i for i, c in enumerate(dcalib.ALPHABET) if c in AROMATIC}
        for scope in SCOPES[m.family]:
            sel = [i for i, k in enumerate(akeys) if k.startswith(f"{scope}|")]
            if not sel:
                continue
            block = msa[np.ix_(sel, motif_cols)]
            n = block.size
            k_arom = int(np.isin(block, list(arom_idx)).sum())
            k_gap = int((block == 0).sum())
            phat = k_arom / n
            # Wilson 95% interval
            z95 = 1.959964
            denom = 1 + z95**2 / n
            centre = (phat + z95**2 / (2 * n)) / denom
            half = z95 * np.sqrt(phat * (1 - phat) / n + z95**2 / (4 * n**2)) / denom
            rows.append({**base, "measure": "aromatic_fraction", "scope": scope,
                         "observed": round(phat, 4), "control_mean": np.nan,
                         "fold": np.nan, "z": np.nan, "p_perm": np.nan,
                         "ci_low": round(max(centre - half, 0), 4),
                         "ci_high": round(min(centre + half, 1), 4),
                         "gap_fraction": round(k_gap / n, 4), "n_observations": n})

        # binder-vs-SOCS2 specificity (SOCS family only)
        if m.family == "socs":
            b = sdp[sdp.axis == "socs_binding"].set_index("column")["mi_bits"]
            obs, p, z = permutation_test(b, motif_cols, ctrl, args.n_perm, SEED)
            ctrl_mean = float(b.reindex(ctrl).dropna().mean())
            rows.append({**base, "measure": "specificity_mi", "scope": "binder_vs_SOCS2",
                         "observed": round(obs, 4), "control_mean": round(ctrl_mean, 4),
                         "fold": round(obs / ctrl_mean, 3) if ctrl_mean else np.nan,
                         "z": round(z, 3) if z == z else np.nan, "p_perm": p})

    res = pd.DataFrame(rows)
    valid = res.p_perm.notna()
    res.loc[valid, "q_bh"] = cons.benjamini_hochberg(res.loc[valid, "p_perm"].values)
    res["significant_q05"] = res.q_bh < 0.05
    res = res.sort_values(["measure", "p_perm"])
    res.to_csv(args.outdir / "motif_tests.csv", index=False)
    res[["motif", "measure", "scope", "p_perm", "q_bh", "significant_q05"]].to_csv(
        args.outdir / "multiple_testing.csv", index=False)

    LOG.info("declared family: %d tests, %d significant at q<0.05",
             int(valid.sum()), int(res.significant_q05.sum()))
    show = ["motif", "scope", "n_residues_mapped", "interface_overlap", "observed",
            "control_mean", "fold", "z", "p_perm", "q_bh", "significant_q05"]
    for measure, g in res.groupby("measure"):
        print(f"\n=== {measure} ===")
        print(g[show].to_string(index=False))


if __name__ == "__main__":
    main()
