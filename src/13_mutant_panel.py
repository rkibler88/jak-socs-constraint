"""Step 13 - Evolutionary scoring of a declared SOCS1 variant panel.

Projects an externally specified set of SOCS1 point substitutions onto the
alignment and reports, for each position, the constraint in every SOCS scope,
the binder-versus-contrast specificity with its rank among all SOCS columns,
the per-paralogue consensus residue, and the frequency of the residue the
substitution introduces among SOCS1 orthologues.

The frequency column is the one that carries interpretation. A substitution
that introduces the residue most orthologues already carry is a reversion to
the consensus state, which is a different proposition from one that introduces
a residue the family never samples, even when the two positions are equally
constrained.

The panel is declared here as a constant rather than taken from the command
line because it is a fixed experimental input to the manuscript, and recording
it in source keeps the mapping auditable. Positions are full-length UniProt
SOCS1 (O15524) numbering; the construct labels are retained so the table can be
joined against the experimental ranking.
"""

from __future__ import annotations

import argparse
import collections
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_sequence as dcalib

LOG = logging.getLogger("step13")

SOCS_SCOPES = ("SOCS1", "SOCS3", "SOCS2")
SPECIFICITY_AXIS = "socs_binding"

# SOCS chain of each experimental complex, and the coordinate map that carries
# its organism's numbering. Interface residues are in these organisms' canonical
# sequences, not human SOCS1's.
STRUCTURE_REFERENCES: dict[str, str] = {
    "6C7Y": "SOCS1_CHICK_6C7Y",   # chicken SOCS1 B6RCQ2, bound to human JAK1
    "4GL9": "SOCS3_MOUSE_4GL9",   # mouse SOCS3 O35718, bound to mouse JAK2
}


@dataclass(frozen=True)
class Variant:
    """One declared point substitution.

    Attributes:
        construct: Label in construct numbering, as used experimentally.
        wild_type: One-letter wild-type residue.
        position: Full-length SOCS1 (O15524) residue number.
        mutant: One-letter substituted residue.
    """

    construct: str
    wild_type: str
    position: int
    mutant: str

    @property
    def label(self) -> str:
        """Substitution in full-length numbering, e.g. ``"H54F"``."""
        return f"{self.wild_type}{self.position}{self.mutant}"


# AI-guided mutational screen panel. Construct numbering is offset -51 from
# full-length SOCS1; both are recorded so neither has to be inferred.
DECLARED_PANEL: tuple[Variant, ...] = (
    Variant("H3F", "H", 54, "F"),
    Variant("V35T", "V", 86, "T"),
    Variant("Q57E", "Q", 108, "E"),
    Variant("F62Y", "F", 113, "Y"),
    Variant("A69G", "A", 120, "G"),
    Variant("V77I", "V", 128, "I"),
    Variant("R90K", "R", 141, "K"),
    Variant("E91L", "E", 142, "L"),
    Variant("R142K", "R", 193, "K"),
    Variant("R150K", "R", 201, "K"),
    Variant("S154A", "S", 205, "A"),
    Variant("F156Y", "F", 207, "Y"),
)


def paralogue_rows(alignment: dict[str, str], paralogue: str) -> list[str]:
    """Alignment keys belonging to one paralogue's orthologue set."""
    return [k for k in alignment if k.startswith(f"{paralogue}|")]


def column_composition(
    alignment: dict[str, str], keys: list[str], column: int
) -> tuple[str, float, collections.Counter]:
    """Consensus residue, its frequency, and the full residue count at a column.

    Args:
        alignment: Aligned sequences keyed by record name.
        keys: Subset of keys to tabulate.
        column: Zero-based alignment column.

    Returns:
        Tuple of (consensus symbol, consensus frequency, residue counter).
    """
    counts = collections.Counter(alignment[k][column] for k in keys)
    total = sum(counts.values())
    symbol, n = counts.most_common(1)[0]
    return symbol, n / total, counts


def interface_columns(
    interface: pd.DataFrame, maps: dict, structure_references: dict[str, str]
) -> dict[str, set[int]]:
    """Alignment columns occupied by SOCS-side interface residues, per structure.

    Interface residues are numbered in the canonical sequence of the structure's
    own organism -- chicken SOCS1 (B6RCQ2) for 6C7Y, mouse SOCS3 (O35718) for
    4GL9 -- so they cannot be compared to human SOCS1 positions directly. They
    are routed through each structure's own coordinate map into alignment
    columns, which is the only coordinate system shared across organisms.

    Args:
        interface: Stage 07 interface residue table.
        maps: Parsed coordinate maps.
        structure_references: Structure label to its key in the SOCS map.

    Returns:
        Structure label to the set of alignment columns it contacts.
    """
    columns: dict[str, set[int]] = {}
    for structure, reference in structure_references.items():
        canonical_to_column = {
            int(k): v for k, v in
            maps["socs"]["references"][reference]["canonical_to_column"].items()
        }
        residues = interface[(interface.side == "SOCS")
                             & (interface.structure == structure)].canonical_residue
        mapped = {canonical_to_column[int(r)] for r in residues
                  if int(r) in canonical_to_column}
        if len(mapped) < 0.5 * len(residues):
            LOG.warning("%s: only %d of %d interface residues mapped to columns",
                        structure, len(mapped), len(residues))
        columns[structure] = mapped
        LOG.info("%s: %d SOCS interface residues -> %d alignment columns",
                 structure, len(residues), len(mapped))
    return columns


def score_panel(
    panel: tuple[Variant, ...],
    alignment: dict[str, str],
    canonical_to_column: dict[int, int],
    constraint: pd.DataFrame,
    specificity: pd.DataFrame,
    interface_cols: dict[str, set[int]],
) -> pd.DataFrame:
    """Score every declared variant against the orthologue alignment.

    Args:
        panel: Declared substitutions.
        alignment: SOCS family alignment.
        canonical_to_column: SOCS1 canonical residue to alignment column.
        constraint: Per-column constraint, all scopes (stage 05).
        specificity: Per-column specificity (stage 06).
        interface: Structural interface residues (stage 07).

    Returns:
        One row per variant, sorted by descending specificity.

    Raises:
        ValueError: If a declared wild-type residue disagrees with the
            alignment, which would mean the panel and the reference sequence
            are on different numbering.
    """
    groups = {p: paralogue_rows(alignment, p) for p in SOCS_SCOPES}
    human_key = next(k for k in alignment if k.startswith("REF|SOCS1"))
    mi = specificity[specificity.axis == SPECIFICITY_AXIS].set_index("column")

    rows: list[dict[str, object]] = []
    for variant in panel:
        column = canonical_to_column.get(variant.position)
        if column is None:
            LOG.warning("%s: SOCS1 %d is outside the alignment window",
                        variant.label, variant.position)
            continue

        observed = alignment[human_key][column]
        if observed != variant.wild_type:
            raise ValueError(
                f"{variant.label}: alignment carries {observed!r} at SOCS1 "
                f"{variant.position}, panel declares {variant.wild_type!r}; "
                "numbering disagrees"
            )

        row: dict[str, object] = {
            "construct": variant.construct,
            "mutation": variant.label,
            "SOCS1_residue": variant.position,
            "column": column,
        }
        for structure, cols in interface_cols.items():
            row[f"at_{structure}_interface_column"] = column in cols
        row["at_any_interface_column"] = any(column in c for c in interface_cols.values())
        for scope in SOCS_SCOPES:
            hit = constraint[(constraint.scope == scope) & (constraint.column == column)]
            row[f"jsd_{scope}"] = round(float(hit.jsd.iloc[0]), 4) if len(hit) else np.nan
            consensus, frequency, counts = column_composition(alignment, groups[scope], column)
            row[f"{scope}_consensus"] = consensus
            row[f"{scope}_consensus_freq"] = round(frequency, 4)
            if scope == "SOCS1":
                total = sum(counts.values())
                row["freq_of_introduced_residue"] = round(
                    counts.get(variant.mutant, 0) / total, 4)
                row["reverts_to_consensus"] = bool(consensus == variant.mutant)
        if column in mi.index:
            row["specificity_mi_bits"] = round(float(mi.loc[column, "mi_bits"]), 4)
            row["specificity_percentile"] = int(
                round(100.0 * float((mi.mi_bits < mi.loc[column, "mi_bits"]).mean())))
            row["gapfrac_contrast"] = round(float(mi.loc[column, "gapfrac_nonbinder"]), 4)
        rows.append(row)

    return pd.DataFrame(rows).sort_values("specificity_mi_bits", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("results/manuscript_analysis"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    maps = json.loads((args.outdir / "coordinate_maps.json").read_text())
    canonical_to_column = {
        int(k): v for k, v in
        maps["socs"]["references"]["SOCS1"]["canonical_to_column"].items()
    }
    alignment = dcalib.read_fasta(args.outdir / "alignments" / "socs_family.fasta")

    iface_cols = interface_columns(
        pd.read_csv(args.outdir / "interface_residues.csv"),
        maps,
        STRUCTURE_REFERENCES,
    )
    table = score_panel(
        DECLARED_PANEL,
        alignment,
        canonical_to_column,
        pd.read_csv(args.outdir / "residue_constraint.csv").query("family == 'socs'"),
        pd.read_csv(args.outdir / "sdp_table.csv"),
        iface_cols,
    )
    table.to_csv(args.outdir / "mutant_panel_scores.csv", index=False)

    reversions = table[table.reverts_to_consensus]
    LOG.info("%d of %d variants revert to the SOCS1 orthologue consensus: %s",
             len(reversions), len(table), ", ".join(reversions.mutation))
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
