"""Stage 12 - package the analysis outputs into the published tree.

The analysis stages read and write a flat working directory using short
internal names. This stage is the single place where those names are mapped to
the organised, numbered layout that ships in the repository, so the mapping is
auditable in one file rather than spread across eleven scripts.

  --pack     working dir -> data/ + results/   (default; run after stage 11)
  --unpack   published tree -> working dir     (to re-run a stage offline
                                                against the shipped outputs)

Writes `results/tables/MANIFEST.csv` describing every published file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import shutil
from pathlib import Path

LOG = logging.getLogger("package")

#: internal filename -> (published subdirectory, published filename, description)
FILEMAP: dict[str, tuple[str, str, str]] = {
    # ---------------------------------------------------------------- inputs
    "coordinate_maps.json": ("data/reference", "coordinate_maps.json",
        "Bijective canonical-residue <-> alignment-column maps for all 9 reference rows"),
    "reference_signatures.json": ("data/reference", "human_reference_signatures.json",
        "Pfam/InterPro signatures read off the six human reference proteins"),
    "feature_ranges.csv": ("data/reference", "uniprot_domain_ranges.csv",
        "UniProt domain boundaries per reference, with alignment coverage and column span"),
    "feature_projection.csv": ("data/reference", "uniprot_domain_projection.csv",
        "Every annotated residue projected to its alignment column"),
    "column_provenance.csv": ("data/reference", "alignment_column_provenance.csv",
        "Per-column identity: canonical residue and feature for each reference, gap fraction"),
    "paralogue_assignment.csv": ("data/reference", "paralogue_assignment.csv",
        "Per-entry paralogue call with identity, coverage, bitscore margin and rejection reason"),
    "references.csv": ("data/reference", "literature_references.csv",
        "32 cited sources with PMID, DOI and the role each plays in the design"),
    # -------------------------------------------------------- depth / retrieval
    "depth_scoping.csv": ("results/tables", "01_depth_scoping_designs.csv",
        "All 12 candidate designs with N, L, N_eff, depth-gate threshold and ratio"),
    "retrieval_yield.csv": ("results/tables", "01_retrieval_yield_by_route_tier.csv",
        "Entries and distinct species per retrieval route x taxonomic tier x family"),
    "per_pair_species.csv": ("results/tables", "01_paired_species_per_paralogue_pair.csv",
        "Paired-species counts for the six JAK x SOCS paralogue combinations"),
    "retrieval_summary.json": ("results/tables", "02_retrieval_summary.json",
        "Retrieval totals, rejection reasons and decoy-diversion counts"),
    "alignment_quality.csv": ("results/tables", "03_alignment_quality_per_paralogue.csv",
        "Per-paralogue sequence counts, N_eff and mean column occupancy"),
    "alignment_meta.json": ("results/tables", "03_alignment_dimensions.json",
        "Alignment dimensions, extraction window and mean pairwise identity per family"),
    "provenance_report.json": ("results/tables", "04_provenance_assertions.json",
        "Outcome of the column-provenance assertions for both families"),
    # ------------------------------------------------------------- constraint
    "residue_constraint.csv": ("results/tables", "05_residue_constraint.csv",
        "Per-residue entropy, conservation and JSD with bootstrap CIs, in UniProt numbering"),
    "module_constraint_summary.csv": ("results/tables", "05_module_constraint.csv",
        "Module-level mean constraint with bootstrap 95% CIs per paralogue"),
    "module_deficit_test.csv": ("results/tables", "05_module_deficit_test.csv",
        "Permutation test that the SOCS1/SOCS3-minus-SOCS2 constraint deficit is "
        "module-specific, in three variants: gap-penalised, residues-only, and "
        "restricted to columns SOCS2 occupies at >50%"),
    "mutant_panel_scores.csv": ("results/tables", "13_mutant_panel_scores.csv",
        "Declared SOCS1 variant panel scored for constraint, specificity rank, "
        "per-paralogue consensus and orthologue frequency of the introduced residue"),
    "module_gap_decomposition.csv": ("results/tables",
        "05_module_constraint_gap_decomposition.csv",
        "Module constraint with and without the gap penalty, separating absence from divergence"),
    "constraint_scope_summary.csv": ("results/tables",
        "05_constraint_effective_sample_sizes.csv",
        "Sequence counts, DCA-convention N_eff and Kish effective n per paralogue"),
    # ------------------------------------------------------------ specificity
    "sdp_table.csv": ("results/tables", "06_specificity_per_column.csv",
        "Per-column mutual information, permutation z/p/q, per-group gap fraction and consensus"),
    "sdp_module_summary.csv": ("results/tables", "06_specificity_by_module.csv",
        "Module-level specificity with and without deletion-driven columns"),
    "sdp_gap_decomposition.csv": ("results/tables", "06_specificity_gap_decomposition.csv",
        "Fraction of each module's specificity signal attributable to deletion"),
    "sdp_feature_summary.csv": ("results/tables", "06_specificity_by_feature.csv",
        "Specificity summary grouped by annotated feature"),
    "sdp_axis_overlap.csv": ("results/tables", "06_specificity_axis_overlap.csv",
        "Overlap of top specificity positions between contrasts vs hypergeometric expectation"),
    "sdp_nulls.npz": ("results/tables", "06_specificity_null_distributions.npz",
        "Permutation null summaries per specificity axis"),
    # ------------------------------------------------------------- structures
    "structural_contacts.csv": ("results/tables", "07_structural_contact_pairs.csv",
        "Inter-chain residue pairs within 8 A, in each structure's own canonical numbering"),
    "interface_residues.csv": ("results/tables", "07_interface_residues.csv",
        "Interface residue sets per structure and side, deduplicated across crystal copies"),
    "structural_sasa.csv": ("results/tables", "07_unbound_chain_sasa.csv",
        "Per-residue solvent-accessible surface area of the isolated chains"),
    "structures_meta.json": ("results/tables", "07_structure_metadata.json",
        "Chain pairing contact counts, resolved ranges and citations"),
    # --------------------------------------------------------- interface test
    "interface_constraint_tests.csv": ("results/tables", "08_interface_constraint_tests.csv",
        "Primary test: interface vs matched-control constraint, effect sizes and p/q values"),
    "matched_control_sets.csv": ("results/tables", "08_matched_control_residues.csv",
        "Every residue entering the matched comparison with its covariate stratum"),
    "interface_power_analysis.csv": ("results/tables", "08_interface_power_curve.csv",
        "Simulated power against constraint shift, per structure and side"),
    "interface_min_detectable_effect.csv": ("results/tables",
        "08_interface_min_detectable_effect.csv",
        "Smallest constraint shift detectable at 80% power, in control SD units"),
    "interface_nulls.npz": ("results/tables", "08_interface_null_distributions.npz",
        "Within-stratum permutation nulls for the interface tests"),
    # ------------------------------------------------------------ motif tests
    "motif_tests.csv": ("results/tables", "09_declared_motif_tests.csv",
        "The pre-declared motif family: constraint, specificity and aromaticity per motif"),
    "multiple_testing.csv": ("results/tables", "09_multiple_testing_correction.csv",
        "Benjamini-Hochberg correction across the declared motif family"),
    # --------------------------------------------------------- coupling power
    "depth_power_report.csv": ("results/tables", "10_coupling_depth_gate.csv",
        "Depth diagnostics for the best design against the published inter-protein threshold"),
    "dca_null_comparison.csv": ("results/tables", "10_coupling_vs_scrambled_null.csv",
        "Observed inter-block coupling statistics against the species-scrambled null"),
    "dca_null_distributions.npz": ("results/tables", "10_coupling_null_distributions.npz",
        "Observed score matrix and the species-scrambled replicate statistics"),
    # ------------------------------------------------------------ consistency
    "profile_concordance.csv": ("results/tables", "11_constraint_profile_concordance.csv",
        "Rank correlation of constraint profiles between paralogues within each family"),
    "interface_meta_analysis.csv": ("results/tables", "11_interface_meta_analysis.csv",
        "Inverse-variance meta-analysis of the interface effect, with heterogeneity"),
    "leave_one_paralogue_out.csv": ("results/tables", "11_kir_leave_one_binder_out.csv",
        "KIR specificity recomputed dropping each inhibitory paralogue in turn"),
    # ------------------------------------------------------------------- logs
    "depth_scoping_queries.json": ("results/logs", "uniprot_queries.json",
        "Exact UniProt query strings and result counts for every scoping query"),
}

#: directories copied wholesale (internal subdir -> published subdir)
DIRMAP = {"sequences": "data/sequences", "alignments": "data/alignments"}


def sha256(path: Path, limit: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(limit):
            h.update(chunk)
    return h.hexdigest()[:16]


def pack(workdir: Path, root: Path) -> None:
    rows = []
    for subdir, published in DIRMAP.items():
        src = workdir / subdir
        if not src.is_dir():
            LOG.warning("missing directory %s, skipped", src)
            continue
        dest = root / published
        dest.mkdir(parents=True, exist_ok=True)
        for f in sorted(src.iterdir()):
            if not f.is_file():
                continue
            out = dest / f.name.lower()
            shutil.copy2(f, out)
            rows.append({"published_path": str(out.relative_to(root)),
                         "internal_name": f"{subdir}/{f.name}",
                         "bytes": out.stat().st_size, "sha256_16": sha256(out),
                         "description": ("Paralogue orthologue set"
                                         if subdir == "sequences"
                                         else "Per-family domain alignment")})

    for internal, (subdir, published, desc) in FILEMAP.items():
        src = workdir / internal
        if not src.exists():
            LOG.warning("missing %s, skipped", internal)
            continue
        dest = root / subdir
        dest.mkdir(parents=True, exist_ok=True)
        out = dest / published
        shutil.copy2(src, out)
        rows.append({"published_path": str(out.relative_to(root)),
                     "internal_name": internal, "bytes": out.stat().st_size,
                     "sha256_16": sha256(out), "description": desc})

    manifest = root / "results/tables/MANIFEST.csv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["published_path", "internal_name",
                                           "bytes", "sha256_16", "description"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["published_path"]))
    LOG.info("packaged %d files; manifest at %s", len(rows), manifest)


def unpack(workdir: Path, root: Path) -> None:
    """Reverse the mapping, so a stage can be re-run against shipped outputs."""
    workdir.mkdir(parents=True, exist_ok=True)
    n = 0
    for subdir, published in DIRMAP.items():
        src = root / published
        if not src.is_dir():
            continue
        dest = workdir / subdir
        dest.mkdir(parents=True, exist_ok=True)
        for f in sorted(src.iterdir()):
            if f.is_file():
                # restore the capitalised paralogue filenames the scripts expect
                name = f.name.upper() if f.suffix == ".fasta" and "_" not in f.stem else f.name
                shutil.copy2(f, dest / name)
                n += 1
    for internal, (subdir, published, _) in FILEMAP.items():
        src = root / subdir / published
        if src.exists():
            shutil.copy2(src, workdir / internal)
            n += 1
    LOG.info("unpacked %d files into %s", n, workdir)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--workdir", type=Path, default=Path("work/analysis"))
    p.add_argument("--root", type=Path, default=Path("."))
    g = p.add_mutually_exclusive_group()
    g.add_argument("--pack", action="store_true", default=True)
    g.add_argument("--unpack", action="store_true")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.unpack:
        unpack(args.workdir, args.root)
    else:
        pack(args.workdir, args.root)


if __name__ == "__main__":
    main()
