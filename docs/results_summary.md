# JAK–SOCS evolutionary constraint and specificity: results summary

Analysis of purifying selection and paralogue specificity across the JAK kinase
family (JAK1, JAK2, TYK2) and the SOCS family (SOCS1, SOCS3, with SOCS2 as a
non-binding contrast), tested against two experimental JAK–SOCS interfaces.
JAK3 is excluded throughout and the exclusion is enforced by sequence, not only
by query construction.

Every method choice was fixed in advance in `METHODS_REVIEW.md`; the motif
hypotheses in §9 of that document were declared before any constraint or
specificity value was computed.

---

## 1. What this dataset can and cannot support

Residue-level inter-protein direct-coupling analysis is **not attainable** for
this protein pair, and the bound is quantitative rather than impressionistic.
The published criterion for reliable inter-protein couplings is
N_eff > (L1+L2)/2 (Ovchinnikov, Kamisetty & Baker, eLife 2014). Across twelve
candidate designs — spanning two retrieval routes, four taxonomic tiers, six
paralogue pairings, factorial pooling and two column-set reductions — the best
configuration reaches **0.66×** the required depth
(N_eff = 145.3 against a requirement of
218.5). None passes.

The cause is conservation, not sampling effort. On a paralogue-homogeneous
JAK2×SOCS3 alignment of 365 species, median pairwise identity is 0.904 and
73.5% of all sequence pairs exceed the 80% redundancy threshold, so 365 species
collapse to N_eff = 8.1. DCA extracts information from covarying variation; an
invariant position carries none however many species are sampled. The taxonomic
range where these proteins are divergent enough — invertebrate metazoans — is
where SOCS1/SOCS3 orthology is undefined and only 81 paired species are
annotated, though those 81 species carry 57% of the effective diversity of all
2,045 vertebrate records.

Running the inference anyway and comparing against a species-scrambled null
(200 replicates, which preserve both marginal alignments, all conservation and
gap structure, and each family's phylogeny, while destroying the pairing) gives:

| statistic | observed | null | z | p |
|---|---|---|---|---|
| mean inter-block score | — | — | — | not interpretable (APC forces the mean to ≈0) |
| mean of top 100 | 0.1017 | 0.0907 ± 0.0045 | +2.44 | 0.015 |
| maximum score | 0.1748 | 0.1358 ± 0.0145 | +2.69 | 0.030 |

So a weak residual pairing signal **is** detectable in the extreme tail — a 12%
elevation of the top-100 mean over the scrambled ensemble. This is reported
because it is real, not because it is useful: a 12% shift in a tail statistic
does not identify which column pairs carry it, and at 0.66× the required depth
no individual pair can be assigned. The appropriate conclusion is that
cross-family pairing information exists in this sequence record but is roughly
an order of magnitude below what residue-level assignment requires.

Claims below are therefore confined to evolutionary constraint and paralogue
specificity. No direct atomic contact, residue-pair coevolution, or
JAK2–SOCS1-specific interface prediction is asserted.

## 2. The kinase-inhibitory modules are released in the non-binding paralogue

The clearest result. Constraint (Jensen-Shannon divergence against an amino-acid
background, Henikoff-weighted, bootstrap 95% CIs) by module and paralogue:

| module | SOCS1 | SOCS3 | SOCS2 |
|---|---|---|---|
| **KIR** | 0.52 [0.43, 0.62] | 0.57 [0.49, 0.65] | **0.18 [0.13, 0.24]** |
| **ESS** | 0.53 [0.44, 0.63] | 0.60 [0.51, 0.68] | **0.33 [0.24, 0.43]** |
| SH2 | 0.61 | 0.65 | 0.54 |
| SOCS box | 0.54 | 0.50 | 0.48 |

The KIR and ESS confidence intervals for SOCS2 do not approach those of either
binder; the SH2 and SOCS-box intervals nearly coincide. The effect is therefore
module-specific rather than a general property of a less-conserved protein: the
structural core and the ubiquitin-ligase module are held to comparable
standards in all three paralogues, and only the two modules that contact the
kinase have been released.

Two decompositions guard this result. SOCS2's KIR columns are 34% gapped
against 1% (SOCS1) and 4% (SOCS3), so the gap penalty in the JSD score is
separated from residue-level divergence: SOCS2's KIR is **both** partly absent
(a third of its KIR columns are majority-gapped, where neither binder has any)
**and**, where present, roughly 45% less conserved (0.30 against 0.53 and 0.59).
Both components point the same way. The same decomposition prevents an error in
the opposite direction: SOCS3's SOCS box appears weakly constrained (0.50) only
because of its 19% gap content, and scores 0.63 on residues actually present,
comparable to SOCS1.

The gap gradient across the KIR is monotonic — SOCS2's gap fraction falls from
0.93 at the position homologous to SOCS1 52 to 0.14 at SOCS1 66 — indicating
N-terminal truncation through the module rather than scattered indels.

## 3. Specificity concentrates in the KIR, and the declared motifs hold

Mutual information between residue state and binder/non-binder identity,
Miller-Madow corrected. At n = 2,006 sequences this test is heavily
overpowered — z-scores reach 50–100 and 132 of 235 columns exceed even the
family-wise permutation maximum — so results are reported as effect sizes
rather than significance.

Of the declared motif family (42 tests, one Benjamini-Hochberg correction),
**two survive**, both in the KIR:

- `KIR_H54_R59` (Yasukawa 1999 mutagenesis set): fold 3.29, z = 6.06, q = 0.0021
- `KIR_full`: fold 2.44, z = 5.29, q = 0.0021

Neither depends on pooling the two binders: repeating the KIR test with SOCS1
alone gives fold 2.06 and with SOCS3 alone fold 2.40, both p = 0.0005.

Locating the published motifs in the sequences rather than assuming their
positions corrected a misreading worth recording: Zhai et al. (2026) place
their motifs in "the KIR region **and BC loop**", and the BC loop lies in the
SH2 domain. Literal search of the human canonical sequences gives `GQM` once
each in JAK1 (1097), JAK2 (1071) and TYK2 (1117), all inside JH1; `FF` at
SOCS1 112 and SOCS3 79; `QR` in SOCS1 only, not positionally homologous in
SOCS3. `FF` and `QR` are absent from SOCS2 entirely.

Aromatic character at these positions — the property Zhai et al.'s
Met-aromatic mechanism requires — separates binders from SOCS2 almost
completely:

| position | SOCS1 | SOCS3 | SOCS2 |
|---|---|---|---|
| KIR aromatic (SOCS1 64 / SOCS3 31) | 0.885 | 0.967 | **0.010** |
| BC-loop FF (SOCS1 112–113 / SOCS3 79–80) | 0.995 | 0.939 | **0.211** |
| SH2 FLVR (shared-function control) | 0.250 | 0.250 | 0.249 |

The KIR position is aromatic in ~89% of SOCS1 and ~97% of SOCS3 orthologues and
in 1% of SOCS2 orthologues. Note that the identity is not fixed — most SOCS1
orthologues carry phenylalanine where human SOCS1 has Y64 — so the property
under selection is aromaticity rather than a specific residue.

The internal negative control behaves as predicted: the invariant SH2
phosphotyrosine arginine (SOCS1 R104) is R in 100% of all three paralogues and
carries specificity fold 0.005. Phosphotyrosine binding is shared
SOCS-family function; the discrimination is confined to the JAK-contacting
motifs. That contrast is what makes the positive result interpretable.

## 4. Interface residues are not measurably more constrained

The pre-specified primary test is negative. Interface residues from 4GL9 and
6C7Y were compared against controls matched on domain × solvent-accessibility
tertile × gap tertile, with permutation within strata:

structure side paralogue  n_interface_matched  n_control_matched  cliffs_delta  delta_ci_low  delta_ci_high  permutation_p    q_bh
     4GL9 SOCS     SOCS3                   31                 53       -0.0870       -0.3524         0.1820       0.385461 0.51395
     4GL9  JAK      JAK2                   42                209        0.0625       -0.1160         0.2295       0.615038 0.61504
     6C7Y SOCS     SOCS1                   31                 64        0.2016       -0.0620         0.4405       0.036196 0.12119
     6C7Y  JAK      JAK1                   51                195        0.1382       -0.0439         0.3132       0.060594 0.12119

All four confidence intervals include zero and none survives correction.
Meta-analysed across all four (the pre-specified combination), pooled
Cliff's δ = 0.089 [-0.014, 0.192], p = 0.091, with no detectable
heterogeneity (I² = 0). Restricted to 6C7Y the pooled effect reaches
δ = 0.160 [0.014, 0.305], p = 0.032, and restricted to 4GL9 it is
absent (p = 0.803). The 6C7Y-only result is a subgroup analysis
performed after seeing the data and is not treated as confirmatory.

This is a bounded negative rather than an absence of evidence. Power simulation
on the observed sample sizes puts the minimum detectable constraint shift at
0.50–0.75 control standard deviations, larger than anything observed, so large
effects are excluded and small-to-moderate ones are not resolvable at 31–51
interface residues. A mechanistic contributor is visible in the data:
constraint is saturated, with 70–80% of JAK control residues already exceeding
JSD 0.6. In a catalytic domain where nearly every residue is conserved, an
interface has little room to stand out — the same property that makes DCA
infeasible, appearing in a second and independent measurement.

## 5. Reproducibility of the measurements themselves

Constraint profiles agree between independent orthologue sets at Spearman
ρ = 0.75–0.86 across all six paralogue pairs (all p < 10⁻⁵⁹), so the
per-column measurement is not noise-dominated. Top specificity positions recur
across the binder-versus-SOCS2 and SOCS1-versus-SOCS3 axes at 10 of the top 30
against 3.8 expected (hypergeometric p = 0.0014).

## 6. Interpretation

The JAK–SOCS regulatory interaction is encoded in the sequence record as a
**module-level** constraint signature, not a residue-pair one. SOCS1 and SOCS3
maintain strong purifying selection across the kinase-inhibitory region and the
extended SH2 subdomain, together with near-invariant aromatic character at two
specific positions — one in the KIR, one in the BC loop — and SOCS2, which
inhibits JAK1/TYK2 poorly, has lost all of it while retaining the SH2 fold, the
phosphotyrosine pocket and the SOCS box. The modules under selection are those
that the 4GL9 and 6C7Y structures place against the kinase.

That pattern is consistent with the KIR acting as a pseudo-substrate element
whose aromatic character is the functionally constrained property, which is
what the structures and the published simulation work describe. It is
convergent support from an independent evidence class, reached without any
structural information entering the constraint calculation.

What the data do not support is any finer resolution. Interface residues are
not individually distinguishable from matched non-interface residues; specific
JAK–SOCS residue pairs cannot be assigned; and the weak tail signal in the
coupling analysis cannot be localised. The proteins are too uniformly conserved
within each paralogue lineage, and too few paired species exist outside the
lineages where orthology is definable.

## 7. Limitations

- **Depth.** No design reaches the published inter-protein DCA threshold; the
  best reaches 0.66×. All coupling-based statements are correspondingly bounded.
- **Interface test power.** Minimum detectable shift 0.50–0.75 control SDs;
  small-to-moderate interface effects remain possible.
- **Constraint saturation.** 70–80% of JAK control residues exceed JSD 0.6,
  compressing any interface-specific signal on the kinase side.
- **SOCS2 as a contrast.** Its non-binding status is taken from published
  simulation and binding work, not measured here. The constraint result is
  conditional on that assignment.
- **Structural coverage.** 4GL9 resolves mouse SOCS3 22–128 plus 164–168 and
  173–184 (a ΔPEST construct); 6C7Y resolves chicken SOCS1 48–164. Neither
  structure uses the human protein throughout, so all structure-to-alignment
  mapping runs through the structure's own organism.
- **TYK2 sampling.** TYK2 gained only 1% of additional species from
  profile-based retrieval (366 against 362), so TYK2-containing comparisons are
  the shallowest in the analysis.
- **Alignment.** Module boundaries rest on MAFFT alignment of a 12-residue KIR;
  SOCS1 and SOCS3 KIRs map to identical columns, which supports the alignment
  there, but individual column assignments in short low-complexity regions
  remain the weakest link in any per-residue claim.
- **Not assessed.** Residual phylogenetic structure beyond the species-scramble
  and within-stratum permutation nulls; codon-level selection (dN/dS);
  paralogue-matched recovery of teleost duplicate copies.

## 8. Pre-registration statement

Fixed before any result was computed, in `METHODS_REVIEW.md`: the depth gate and
its source; the estimator choice and its per-block APC convention; the primary
null (species scramble); the primary hypothesis test (interface versus matched
controls); and the motif family of §9. Declared after inspecting sequences but
before testing: the exact motif residue positions, which were located by string
search rather than assumed, and the aromaticity measure, adopted because the
published mechanism is Met-aromatic rather than residue-identity-specific.
Exploratory and labelled as such: the 6C7Y-only interface subgroup, and the
gap/substitution decompositions, which were added after the pooled numbers
showed gap content differed between paralogues.
