# Replacement for the co-evolution Results paragraph

Replaces the paragraph beginning *"Co-evolutionary analysis of the JAK-SOCS
regulatory axis using family-level direct coupling analysis (DCA) across
JAK1/2/3/TYK2 and SOCS1/SOCS3..."*.

Prose is matched to the rest of the manuscript: no em dashes, mean sentence
length 24 words against the draft's 23.8, en dashes for numeric ranges, and US
spelling (ortholog, paralog, penalized, localized). Statistics follow the
house format used in the mutational-screen paragraphs, with spaces around
`=` and semicolon-separated values inside parentheses. SOCS residue numbering
is full-length throughout, matching Tables 2 and 3. Figure and table callouts
are `Fig. Sx` and `Table Sx` placeholders.

---

## Paragraph 1 | Evolutionary constraint on the inhibitory architecture

To test whether the inhibitory architecture is independently supported by
evolution, we assembled paralog-resolved ortholog sets for JAK1, JAK2, and TYK2
(1,591 sequences) and for SOCS1, SOCS3, and SOCS2 (2,011 sequences from 704,
704, and 598 species). Paralog identity was assigned by sequence rather than by
annotation, and per-residue purifying selection was quantified as
Jensen-Shannon divergence (JSD) against an amino acid background (Fig. Sx,
Table Sx). SOCS2 lacks an annotated kinase inhibitory region (KIR) and inhibits
JAK1 and TYK2 poorly, and it was therefore used as a within-family contrast for
constraint. It was not included in any coupling alignment. On this contrast,
constraint separated the two inhibitory paralogs from SOCS2 within the modules
that contact the kinase. In the KIR (SOCS1 55–66; SOCS3 22–33, which occupies
the same alignment columns), constraint was 0.52 for SOCS1 (95% CI 0.43–0.62)
and 0.57 for SOCS3 (0.49–0.65), against 0.18 for SOCS2 (0.13–0.24). In the
extended SH2 subdomain (ESS; SOCS1 67–78), the corresponding values were 0.53
(0.44–0.63), 0.60 (0.51–0.68), and 0.33 (0.24–0.43). We tested the
inhibitor-minus-SOCS2 difference at each column against the rest of the protein
by permutation (10,000 permutations). The deficit was confined to these two
modules (KIR excess = 0.306, z = 6.59, q = 4.0 × 10⁻⁴; ESS excess = 0.163,
z = 3.41, q = 0.003), with no excess in the SH2 domain (0.017, q = 0.28) or the
SOCS box (-0.042, q = 0.95). The effect is therefore module-specific and not a
property of a generally less conserved protein. Because SOCS2 is truncated
through the KIR, its KIR columns are 34% gapped against 1% in SOCS1 and 4% in
SOCS3, and gap content declines from 93% to 14% across the module. We therefore
repeated the test on residue positions alone, and again on the eight of twelve
KIR columns that SOCS2 retains at above 50% occupancy. The deficit persisted in
both cases (q = 4.0 × 10⁻⁴ and q = 0.0016). It does not follow from absence of
the domain, but from weaker constraint at the positions SOCS2 still carries.
Mutual information between residue identity and inhibitory status localized the
distinction to the same region. The KIR was enriched 2.44-fold over matched
columns elsewhere in the protein (z = 5.29, q = 0.0021), and the His54-Arg59
segment 3.29-fold (z = 6.06, q = 0.0021). These were the only two of 42 tests
across an eleven-motif pre-declared panel to survive correction.

## Paragraph 2 | Residue-level resolution and its limit

Resolving these modules to individual positions recovered the contacts
identified in the structural model. Phe58 anchors the hydrophobic core against
JAK2 Leu1026 and Leu1029, and it is phenylalanine in 99% of SOCS1 and 98% of
SOCS3 orthologs, as SOCS3 Phe25, but is deleted in 86% of SOCS2 sequences. The
adjacent position SOCS1 64 and SOCS3 31 is aromatic in 88.5% and 96.7% of
orthologs, against 1.0% in SOCS2. The BC-loop pair SOCS1 112–113 and SOCS3
79–80 is aromatic in 99.5% and 93.9%, against 21.1%. In each case the residue
identity is not fixed, and most SOCS1 orthologs carry phenylalanine where human
SOCS1 carries Tyr64. The selected property is therefore aromatic character
rather than a specific residue, consistent with the methionine-aromatic
contacts these positions make with the invariant JAK GQM motif (JAK1
1097–1099; JAK2 1071–1073; TYK2 1117–1119). As an internal control, the SH2
phosphotyrosine-binding arginine Arg104 is invariant in all three paralogs and
carries no discriminating information (0.005-fold). This confirms that the
signal is confined to kinase-contacting positions and does not reflect general
SOCS1 conservation. Residue-pair co-evolution could not be resolved at this
sequence depth. The deepest admissible JAK by SOCS alignment reached an
effective depth of 145 against the 219 required for inter-protein coupling
inference (0.66-fold), because JAK and SOCS orthologs are too conserved to
covary, with 73.5% of sequence pairs exceeding the 80% redundancy threshold.
Although the strongest couplings were marginally elevated over a 200-replicate
species-scrambled null (mean of top 100 scores = 0.1017 against 0.0907 ± 0.0045;
z = +2.44, p = 0.015), this 12% shift does not localize to individual residue
pairs (Fig. Sx). Interface residues taken from the experimental JAK-SOCS
complexes were likewise no more constrained than controls matched on domain,
burial, and alignment occupancy (pooled Cliff's δ = 0.089, 95% CI -0.014 to
0.192, p = 0.091). Together, these results support the inhibitory architecture
at the level of modules and individual functional positions, but not as pairwise
JAK-SOCS co-variation.

---

## Optional paragraph 3 | For the mutational-screen section

Place after the Kruskal-Wallis paragraph. This paragraph connects the
evolutionary analysis to the variant ranking, and is worth including even if
paragraphs 1 and 2 are shortened. All values are from
`13_mutant_panel_scores.csv`.

> We next scored the twelve screened positions against the ortholog sets to
> obtain an independent evolutionary read-out of the ranking (Table Sx). Ten of
> the twelve lie above the median paralog specificity of the SOCS1 alignment.
> The four that contact JAK in the SOCS1-JAK1 complex, His54, Val86, Gln108, and
> Phe113, are the most constrained of the panel (mean JSD = 0.811 against 0.484
> for the remaining eight). The screen therefore selected positions that are
> functionally loaded rather than merely tolerant. The variants themselves
> divide sharply. The three top-ranked in their respective backgrounds restore
> residues that most SOCS1 orthologs already carry: V128I introduces the
> isoleucine found in 72.3% of SOCS1 sequences, R201K the lysine in 69.3%, and
> R141K the lysine in 62.4%, each replacing a human-specific minority state. The
> remaining nine introduce residues present in 0.0% to 8.8% of orthologs. The
> two variants ranked lowest in the JAK2-WT background are the clearest cases.
> His54 is 96.9% invariant across SOCS1 and sits at the 99th percentile of
> inhibitor-versus-SOCS2 specificity, whereas Val128 and Arg201 fall at the 27th
> and 39th percentiles. The screen thus converged on consensus reversions at
> positions under shared rather than inhibitor-specific constraint, a pattern
> consistent with stabilization that leaves the discriminating positions
> untouched.

---

## Suggested supplementary items

| Item | Content | Source file |
|---|---|---|
| Fig. Sx | Constraint profiles along SOCS1, SOCS3, and SOCS2 and along JAK1, JAK2, and TYK2, with module tracks and module means with CIs | `fig3_constraint_profiles.*` |
| Fig. Sx | Aromaticity at declared motif positions; specificity fold-enrichment across the declared panel | `fig5_declared_motif_tests.*` |
| Fig. Sx | Depth gate across candidate designs; observed coupling against the species-scrambled null | `fig2_depth_gate_and_coupling_null.*` |
| Fig. Sx | Interface versus matched-control constraint, with effect sizes and CIs | `fig6_interface_constraint_test.*` |
| Table Sx | Module constraint per paralog with CIs and gap decomposition | `05_module_constraint.csv`, `05_module_constraint_gap_decomposition.csv` |
| Table Sx | Module-specificity permutation test, all three variants | `05_module_deficit_test.csv` |
| Table Sx | Declared motif panel: constraint, specificity, aromaticity, and q-values | `09_declared_motif_tests.csv` |
| Table Sx | Evolutionary scores for the twelve screened positions (paragraph 3) | `13_mutant_panel_scores.csv` |

---

## Notes for the authors

### What this version does differently

1. **JAK3 is removed**, and the exclusion is enforced by sequence rather than
   by query construction. A total of 501 retrieved entries were assigned to the
   JAK3 reference and discarded.
2. **The `XXX` placeholder is filled**, with z = +2.44, p = 0.015 against a
   200-replicate species-scrambled null.
3. **The claim that coupling signals "were most coherently localized to the KIR
   region" is not made.** The KIR does stand out, but in constraint and paralog
   specificity rather than in coupling. At 0.66-fold the required depth, the
   coupling matrix cannot localize anything. The module-level conclusion
   survives on evidence the data can carry.
4. **The bound is quantified** rather than described as inconclusive, so a
   reader can see that the analysis was limited by a stated criterion rather
   than by execution.
5. **SOCS2 is used as an internal contrast**, which converts a statement about
   conservation into a statement about what distinguishes an inhibitory SOCS.

### On the use of SOCS2, and the missing-domain objection

SOCS2 is not used in any coupling or DCA alignment. The coupling analysis
(`src/10_coupling_power.py`) uses JAK1, JAK2, and TYK2 against SOCS1 and SOCS3
only, giving six paralog combinations and 2,912 records, and no candidate design
at any stage included SOCS2. Including a paralog that lacks the inhibitory
domain in a coevolution alignment would inject gaps into exactly the columns
the model is meant to couple, and that is not done anywhere here.

SOCS2 appears only as a contrast group in the constraint and specificity
analyses. There, its lack of a KIR is the measurement rather than a
contaminant, because the question asked is what distinguishes a SOCS protein
that inhibits from one that does not. The obvious objection is that a protein
cannot be conserved in a domain it does not have, which would make the KIR
result circular. Two tests reported in paragraph 1 pre-empt this. Removing the
gap scaling leaves the KIR excess at q = 4.0 × 10⁻⁴, and restricting the test
to the eight of twelve KIR columns that SOCS2 retains at above 50% occupancy
leaves it at q = 0.0016 (ESS q = 0.005 and q = 0.0016). The deficit is
therefore carried by divergence at retained positions and not only by the
deletion. All three variants are in `results/tables/05_module_deficit_test.csv`,
produced by `src/05_constraint.py`.

One point is worth stating plainly. Without a non-inhibitory contrast there is
no significant result in the declared panel at all. Of the 42 tests, the 24
that do not involve SOCS2 are motif-versus-matched-control comparisons within a
single paralog, and the best of them, KIR His54-Arg59 constraint in SOCS1 at
1.629-fold, reaches only p = 0.022 and q = 0.20. The contrast is what makes the
analysis inferential rather than descriptive. If SOCS2 specifically is
unwelcome, CIS and CISH (774 species) and SOCS4 through SOCS7 (763 to 1,081
species) are already retrieved and paralog-assigned and could serve the same
role, singly or pooled. Re-running stages 03, 05, 06, and 09 against a
different contrast is about twenty minutes of compute.

### Two claims to reconcile elsewhere in the manuscript

The abstract states that "the KIR, SH2, and SOCS box domains form a highly
conserved inhibitory architecture". This is true of SOCS1 and SOCS3, but the
SH2 domain and the SOCS box are equally constrained in SOCS2 (0.54 and 0.48,
against 0.61 and 0.65 for SH2 and 0.54 and 0.50 for the SOCS box). They are
therefore not what distinguishes an inhibitory SOCS. If the claim is meant to
be specific to inhibition, it should name the KIR and the ESS. Separately,
SOCS3's SOCS box scores 0.503 only because of 19% gap content, and reaches
0.628 on residues actually present, which is worth knowing before the number is
quoted.

The statement that Phe58 anchors the hydrophobic core against JAK2 Leu1026 and
Leu1029 is independently supported here, at 99% Phe in SOCS1, 98% in SOCS3, and
deleted in 86% of SOCS2. The cross-reference is worth making explicit in both
directions.

### One editorial decision left open

The interface test is reported pooled across both complexes (δ = 0.089,
p = 0.091, negative). Restricted to the SOCS1-JAK1 complex alone it is
nominally positive (δ = 0.160, 95% CI 0.014 to 0.305, p = 0.032). That subgroup
was not pre-specified, and the two complexes are not statistically
heterogeneous (Cochran Q = 2.82, p = 0.42, I² = 0). The pooled negative is
therefore the defensible primary result, and the subgroup belongs in the
supplement if it is reported at all. It is excluded from the main text above.
