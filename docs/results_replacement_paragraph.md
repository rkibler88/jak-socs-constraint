# Replacement for the co-evolution Results paragraph

Replaces the paragraph beginning *"Co-evolutionary analysis of the JAK–SOCS
regulatory axis using family-level direct coupling analysis (DCA) across
JAK1/2/3/TYK2 and SOCS1/SOCS3…"*.

Two paragraphs, 716 words (395 + 321; the paragraph replaced was 198). If space
is tight, paragraph 2 can lose its last three sentences — the interface test and
the closing qualifier — and move them to the supplement, leaving 609 → ~500.
JAK3 is excluded;
the `XXX` placeholder is filled; SOCS residue numbering is full-length
throughout, matching Tables 2 and 3. Figure and table callouts are `Fig. Sx` /
`Table Sx` placeholders — suggested items are listed at the end.

A third, optional paragraph for the mutational-screen section follows. It is
the part that makes this analysis load-bearing rather than standalone, and is
worth including even if the two main paragraphs are cut further.

---

## Paragraph 1 — evolutionary constraint on the inhibitory architecture

To test whether the inhibitory architecture identified structurally is
independently supported by evolution, we assembled paralogue-resolved
orthologue sets for JAK1, JAK2 and TYK2 (1,591 sequences) and for SOCS1, SOCS3
and SOCS2 (2,011 sequences from 704, 704 and 598 species), assigning paralogue
identity by sequence rather than by annotation, and quantified per-residue
purifying selection as Jensen–Shannon divergence against an amino-acid
background (Fig. Sx, Table Sx). SOCS2, which lacks an annotated KIR and
inhibits JAK1 and TYK2 poorly, was used as a within-family contrast for
constraint only and was never included in a coupling alignment. On this
contrast, constraint separated the two inhibitory paralogues from SOCS2
specifically within the modules that contact the kinase. In the KIR (SOCS1
55–66; SOCS3 22–33, which occupies the same alignment columns), constraint was
0.52 (95% CI 0.43–0.62) for SOCS1 and 0.57 (0.49–0.65) for SOCS3 against 0.18
(0.13–0.24) for SOCS2, and in the extended SH2 subdomain (SOCS1 67–78) 0.53
(0.44–0.63), 0.60 (0.51–0.68) and 0.33 (0.24–0.43). Testing the inhibitor-minus-
SOCS2 constraint difference per column against the rest of the protein by
permutation (10,000 label permutations) confirmed that the deficit is confined
to these two modules: KIR excess 0.306 (z = 6.59, q = 4.0 × 10⁻⁴) and ESS 0.163
(z = 3.41, q = 0.003), against no excess in the SH2 domain (0.017, q = 0.28) or
SOCS box (−0.042, q = 0.95), so the effect is module-specific rather than a
property of a generally less-conserved protein. Because SOCS2 is truncated
through the KIR — its KIR columns are 34% gapped against 1% and 4% in SOCS1 and
SOCS3, with gap content declining monotonically from 93% to 14% across the
module — we repeated the test on residue positions alone and on the eight of
twelve KIR columns that SOCS2 retains at above 50% occupancy; the deficit
persists in both (q = 4.0 × 10⁻⁴ and q = 0.0016), so it does not follow from the
absence of the domain but from weaker constraint at the positions SOCS2 still
carries. Mutual information between
residue identity and inhibitory status localised the distinction to the same
region: the KIR was enriched 2.44-fold over matched columns elsewhere in the
protein (z = 5.29, q = 0.0021) and the His54–Arg59 segment 3.29-fold (z = 6.06,
q = 0.0021), the only two of 42 tests across an eleven-motif pre-declared panel
to survive correction.

## Paragraph 2 — residue-level resolution and its limit

Resolving these modules to individual positions recovered the contacts
identified in the structural model. Phe58, the hydrophobic anchor packing
against JAK2 Leu1026 and Leu1029, is phenylalanine in 99% of SOCS1 and 98% of
SOCS3 orthologues (as SOCS3 Phe25) but deleted in 86% of SOCS2 sequences. The
adjacent position SOCS1 64 / SOCS3 31 is aromatic in 88.5% and 96.7% of
orthologues against 1.0% in SOCS2, and the BC-loop pair SOCS1 112–113 / SOCS3
79–80 in 99.5% and 93.9% against 21.1%; in each case identity is not fixed —
most SOCS1 orthologues carry Phe where human SOCS1 carries Tyr64 — so the
selected property is aromatic character rather than a specific residue,
consistent with the Met–aromatic contacts these positions make with the
invariant JAK GQM motif (JAK1 1097–1099, JAK2 1071–1073, TYK2 1117–1119). As an
internal control, the SH2 phosphotyrosine-binding arginine (SOCS1 Arg104) is
invariant in all three paralogues and carries no discriminating information
(0.005-fold), confirming that the signal is confined to kinase-contacting
positions rather than reflecting general SOCS1 conservation. Residue-pair
co-evolution, however, could not be resolved: the deepest admissible JAK × SOCS
alignment reaches an effective depth of 145 against the 219 required for
inter-protein coupling inference (0.66×), because JAK and SOCS orthologues are
too conserved to covary — 73.5% of sequence pairs exceed the 80% redundancy
threshold — and although the strongest couplings are marginally elevated over a
200-replicate species-scrambled null (top-100 mean 0.1017 versus 0.0907 ±
0.0045, z = +2.44, p = 0.015), this 12% shift does not localise to individual
residue pairs (Fig. Sx). Interface residues taken from the experimental
JAK–SOCS complexes were likewise not more constrained than controls matched on
domain, burial and alignment occupancy (pooled Cliff's δ = 0.089, 95% CI −0.014
to 0.192, p = 0.091). The evolutionary evidence therefore supports the
inhibitory architecture at the level of modules and individual functional
positions, but not as pairwise JAK–SOCS co-variation.

---

## Optional paragraph 3 — for the mutational-screen section

~190 words. Place after the variant-ranking paragraph (the Kruskal–Wallis
paragraph). All values from `mutant_panel_evolutionary_scores.csv`.

> Scoring the twelve screened positions against the orthologue sets gave an
> independent evolutionary read-out of the ranking (Table Sx). Ten of the twelve
> lie above the median paralogue-specificity of the SOCS1 alignment, and the
> four that contact JAK in the SOCS1–JAK1 complex — His54, Val86, Gln108 and
> Phe113 — are the most constrained of the panel (mean 0.811 versus 0.484 for
> the remaining eight), indicating the screen selected positions that are
> functionally loaded rather than merely tolerant. The variants themselves
> divide sharply. The three top-ranked in their respective backgrounds restore
> residues that most SOCS1 orthologues already carry: V128I introduces the
> isoleucine found in 72.3% of SOCS1 sequences, R201K the lysine in 69.3% and
> R141K the lysine in 62.4%, each replacing a human-specific minority state. The
> remaining nine introduce residues present in 0.0–8.8% of orthologues. The two
> ranked lowest in the JAK2-WT background are the clearest cases: His54 is 96.9%
> invariant across SOCS1 and sits at the 99th specificity percentile, while
> Val128 and Arg201 fall at the 27th and 39th. The screen therefore converged on
> consensus reversions at positions under shared rather than
> inhibitor-specific constraint — stabilising substitutions that leave the
> discriminating positions untouched.

---

## Suggested supplementary items

| Item | Content | Source file |
|---|---|---|
| Fig. Sx | Constraint profiles along SOCS1/SOCS3/SOCS2 and JAK1/JAK2/TYK2 with module tracks; module means with CIs | `fig3_constraint_profiles.*` |
| Fig. Sx | Aromaticity at declared motif positions; specificity fold-enrichment across the declared panel | `fig5_declared_motif_tests.*` |
| Fig. Sx | Depth gate across candidate designs; observed coupling against the species-scrambled null | `fig2_depth_gate_and_coupling_null.*` |
| Fig. Sx | Interface versus matched-control constraint, effect sizes with CIs | `fig6_interface_constraint_test.*` |
| Table Sx | Module constraint per paralogue with CIs and gap decomposition | `05_module_constraint.csv`, `05_module_constraint_gap_decomposition.csv` |
| Table Sx | Declared motif panel: constraint, specificity, aromaticity, q-values | `09_declared_motif_tests.csv` |
| Table Sx | **Evolutionary scores for the twelve screened positions** (paragraph 3) | `mutant_panel_evolutionary_scores.csv` |

---

## Notes for the authors

### What this version does differently

1. **JAK3 is removed**, and the exclusion is enforced by sequence rather than by
   query construction: 501 retrieved entries were assigned to the JAK3
   reference and discarded.
2. **The `XXX` placeholder is filled** — z = +2.44, p = 0.015 against a
   200-replicate species-scrambled null.
3. **The claim that coupling signals "were most coherently localized to the KIR
   region" is not made.** The KIR does stand out, but in constraint and
   paralogue specificity, not in coupling; at 0.66× the required depth the
   coupling matrix cannot localise anything. The module-level conclusion
   survives, on evidence the data can carry.
4. **The bound is quantified** rather than described as inconclusive, so a
   reader can see the analysis was limited by a stated criterion rather than by
   execution.
5. **SOCS2 is used as an internal contrast**, which converts a statement about
   conservation into a statement about what distinguishes an inhibitory SOCS.

### On the use of SOCS2 — scope and the missing-domain objection

SOCS2 is **not** used in any coupling or DCA alignment. The coupling analysis
(`src/10_coupling_power.py`) uses JAK1/JAK2/TYK2 × SOCS1/SOCS3 only — six
paralogue combinations, 2,912 records — and no candidate design at any stage
included SOCS2. Including a paralogue that lacks the inhibitory domain in a
coevolution alignment would inject gaps into exactly the columns the model is
meant to couple, and that is not done anywhere here.

Where SOCS2 appears is as a contrast group in the constraint and specificity
analyses, and there its lack of a KIR is the measurement rather than a
contaminant: the question asked is what distinguishes a SOCS protein that
inhibits from one that does not. The obvious objection — that a protein cannot
be conserved in a domain it does not have, making the KIR result circular — is
pre-empted by two tests reported in the paragraph. Removing the gap scaling
leaves the KIR excess at q = 4.0 × 10⁻⁴, and restricting to the eight of twelve
KIR columns SOCS2 retains at above 50% occupancy leaves it at q = 0.0016 (ESS
q = 0.005 and 0.0016). The deficit is therefore carried by divergence at
retained positions, not only by the deletion. All three variants are in
`results/tables/05_module_deficit_test.csv`, produced by `src/05_constraint.py`.

Worth stating plainly: without a non-inhibitory contrast there is no
significant result in the declared panel at all. Of the 42 tests, the 24 that do
not involve SOCS2 are motif-versus-matched-control comparisons within a single
paralogue, and the best of them — KIR His54–Arg59 constraint in SOCS1, fold
1.629 — reaches only p = 0.022, q = 0.20. The contrast is what makes the
analysis inferential rather than descriptive. If SOCS2 specifically is
unwelcome, CIS/CISH (774 species) and SOCS4–7 (763–1,081 species) are already
retrieved and paralogue-assigned and could serve the same role, singly or
pooled; re-running stages 03, 05, 06 and 09 against a different contrast is
about twenty minutes of compute.

### Two claims to reconcile elsewhere in the manuscript

- The abstract states that "the KIR, SH2, and SOCS box domains form a highly
  conserved inhibitory architecture". This is true of SOCS1 and SOCS3, but the
  SH2 domain and SOCS box are **equally** constrained in SOCS2 (0.54 and 0.48
  against 0.61/0.65 and 0.54/0.50), so they are not what distinguishes an
  inhibitory SOCS. If the claim is meant to be specific to inhibition it should
  name the KIR and ESS. Separately, SOCS3's SOCS box scores 0.503 only because
  of 19% gap content and reaches 0.628 on residues actually present — worth
  knowing before the number is quoted.
- Paragraph 8's statement that Phe58 anchors the hydrophobic core against JAK2
  Leu1026/Leu1029 is independently supported here (99% Phe in SOCS1, 98% in
  SOCS3, deleted in 86% of SOCS2) and the cross-reference is worth making
  explicit in both directions.

### One editorial decision left open

The interface test is reported pooled across both complexes (δ = 0.089,
p = 0.091, negative). Restricted to the SOCS1–JAK1 complex alone it is
nominally positive (δ = 0.160, 95% CI 0.014–0.305, p = 0.032). That subgroup
was not pre-specified and the two complexes are not statistically heterogeneous
(Cochran Q = 2.82, p = 0.42, I² = 0), so the pooled negative is the defensible
primary result and the subgroup belongs in the supplement if it is reported at
all. It is excluded from the main text above.
