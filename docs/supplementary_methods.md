# Supplementary Methods: evolutionary constraint and specificity analysis

Accompanies the co-evolution paragraphs of the Results. A one-sentence pointer
in the main-text Methods is sufficient, e.g.:

> *Evolutionary constraint and specificity analyses, including sequence
> retrieval, paralogue assignment, alignment, coordinate mapping, statistical
> testing and power analysis, are described in full in Supplementary Methods
> §S1–S8. Code and data: [repository URL].*

---

## S1. Sequence retrieval

JAK-family and SOCS-family sequences were retrieved from UniProtKB (accessed
29 September 2026) by InterPro domain signature rather than by gene symbol,
using the REST queries `(xref:interpro-IPR016251) AND (taxonomy_id:33208)` for
the JAK family (IPR016251, non-receptor tyrosine kinase Jak/Tyk2) and
`(xref:interpro-IPR000980) AND (xref:interpro-IPR001496) AND
(taxonomy_id:33208)` for the SOCS family (the conjunction of the SH2 domain and
SOCS box signatures; SH2 alone returns every SH2 protein in the proteome).

Signature retrieval was selected on the basis of a scoping survey of 66,148
entries spanning two retrieval routes and four nested taxonomic tiers
(Gnathostomata 7776, Vertebrata 7742, Chordata 7711, Metazoa 33208). It
recovers 744 species carrying both families against 449 for gene-symbol
queries, and 87% of the JAK-family entries in the 260 additional species carry
no primary gene symbol — i.e. the gain comes from genomes that are sequenced
but not yet annotated. Canonical sequences for the resulting 14,338 unique
entries were fetched in batches of 100 from the `/uniprotkb/accessions`
endpoint. Exact query strings and per-query counts are in the deposited
`uniprot_queries.json`.

## S2. Paralogue assignment

Paralogue identity was assigned by sequence, not by annotation, because gene
symbols are absent or unreliable for most non-model entries. Each sequence was
searched with MMseqs2 v18.8cc5c using

```
mmseqs easy-search <queries> <panel> <out> <tmp> \
  --format-output query,target,pident,alnlen,qlen,tlen,qcov,tcov,bits,evalue \
  -s 7.5 --max-seqs 50 -e 1e-5 --threads 8
```

against a twelve-protein human reference panel: the five targets (JAK1 P23458,
JAK2 O60674, TYK2 P29597, SOCS1 O15524, SOCS3 O14543) plus seven decoys
included so that off-target family members are assigned to their own reference
and excluded on evidence rather than silently absorbed (JAK3 P52333, SOCS2
O14508, SOCS4 Q8WXH5, SOCS5 O75159, SOCS6 O14544, SOCS7 O14512, CIS Q9NSE2).

A sequence was accepted only when its highest-bitscore hit was a target
paralogue, with reference coverage ≥ 0.55, sequence identity ≥ 0.35 and a
bitscore margin ≥ 1.10 over the best hit from any different paralogue. Of
14,338 entries, 5,899 were accepted and 8,437 rejected, the rejections
dominated by decoy assignment (SOCS7 1,612; SOCS5 1,475; SOCS6 1,139; CIS
1,078; SOCS2 1,000; SOCS4 831; JAK3 501). JAK3 is therefore excluded by
sequence, not only by query construction. One sequence was retained per
(paralogue, species), ranked by reference coverage and then length — not by
length alone, which preferentially selects long fragments over complete
orthologues — yielding 617 JAK1, 648 JAK2, 366 TYK2, 704 SOCS1 and 704 SOCS3
species. SOCS2 (598 species) was reconstructed from the same table under
identical thresholds for use as the non-inhibitory contrast. Per-entry
identity, coverage, margin, runner-up reference and rejection reason are
deposited in `paralogue_assignment.csv`.

## S3. Alignment and coordinate mapping

Each family was aligned separately; concatenated records were never aligned.
Alignment used MAFFT v7.526 in two passes. A first pass over full-length
sequences (`--auto --anysymbol --thread 8`) located the region of interest
through the human reference coordinate map, using a window of JAK2 834–1132 and
SOCS1 25–211 (the annotated domain extended by 15 residues of flank). Each
sequence's own residues occupying that window were then extracted ungapped and
realigned by iterative refinement (`--retree 2 --maxiterate 1000 --anysymbol
--thread 8`). This gave a JAK alignment of 1,591 sequences × 319 columns (mean
pairwise identity 0.653, median column gap fraction 0.007) and a SOCS alignment
of 2,011 sequences × 235 columns (0.598, 0.030). Records occupying less than
50% of the first-pass window were dropped (44 of 1,635 JAK records; none of
2,011 SOCS records).

Bijective maps between canonical UniProt numbering and final alignment column
were built for nine reference rows by composing the second-pass map with the
canonical position of each extracted residue, recorded during extraction.
Recovering the offset by substring search is invalid here: a paralogue's
residues within a window defined on a different paralogue are non-contiguous
wherever it carries an insertion. Every map was validated by asserting that
each alignment column holds the residue the canonical sequence carries at that
position; all nine passed. Maps cover JAK1, JAK2, TYK2, SOCS1, SOCS2, SOCS3
and, for structural mapping, mouse JAK2 (Q62120), mouse SOCS3 (O35718) and
chicken SOCS1 (B6RCQ2), which were forced into the alignments so that
structure-to-column mapping runs through each structure's own organism.

Domain boundaries were taken from UniProt feature annotations rather than
hardcoded: JAK1 JH1 875–1153, JAK2 JH1 849–1124, TYK2 JH1 897–1176 (feature
"Protein kinase 2"; the JH2 pseudokinase was excluded); SOCS1 KIR 55–66, ESS
67–78, SH2 79–174, SOCS box 161–210; SOCS3 KIR 22–33, ESS 34–45, SH2 46–142,
SOCS box 177–224; SOCS2 SH2 48–156, SOCS box 151–197, with no KIR or ESS
annotated. Annotations were projected through the validated maps and the
resulting assignment audited: 299 of 319 JAK columns carry a JAK2 residue, of
which 281 (93.98%) fall inside JH1; 187 of 235 SOCS columns carry a SOCS1
residue, of which 156 (83.42%) fall inside an annotated KIR, ESS, SH2 or
SOCS-box span; and no column carries an annotation from the opposing family.

## S4. Constraint

Per-column constraint was quantified as the Jensen–Shannon divergence of the
weighted residue distribution against the Robinson & Robinson background
amino-acid frequencies, scaled by the non-gap fraction, with Shannon entropy
and normalised conservation reported alongside.

Sequences were weighted by Henikoff & Henikoff position-based weights rather
than by the identity-threshold scheme conventional in coevolution analysis. The
latter is inappropriate here: at an 80% identity threshold a set of orthologues
of a single paralogue collapses almost completely (JAK2: effective depth 8.8
from 625 sequences), discarding nearly all the data. Position-based weighting
gives Kish effective sample sizes of 96.3 (JAK2), 99.1 (SOCS3), 127.8 (TYK2),
135.3 (SOCS2), 191.5 (JAK1) and 195.8 (SOCS1). Per-column 95% confidence
intervals came from 1,000 sequence bootstrap resamples; module-level intervals
from 2,000 column resamples.

Because the divergence score is gap-penalised, every module result was also
computed without the gap penalty, so that absence of a region is distinguishable
from divergence within it. This matters in both directions: SOCS2's KIR deficit
is partly deletion (34% gaps, 4 of 12 columns more than half gapped) and partly
divergence (0.303 against 0.529 and 0.590 on residues actually present),
whereas SOCS3's SOCS box scores 0.503 only because of its 19% gap content and
reaches 0.628 on residues present.

## S5. Specificity

Specificity was quantified as the mutual information between column residue
state and group label, with gaps retained as an informative symbol and a
Miller–Madow correction applied to each entropy term. Three contrasts were
tested: inhibitor (SOCS1 + SOCS3) versus SOCS2, SOCS1 versus SOCS3, and JAK1
versus JAK2 versus TYK2. Significance was assessed against 10,000 permutations
of the group labels, which leave each column's residue composition and each
group's size intact and break only the association under test; each
permutation's contingency table was obtained by sparse matrix–vector product,
reducing cost from O(N·L·q) to O(N·L) per group.

At n = 2,006 sequences this test is heavily overpowered: z scores reach 50–100
and 132 of 235 columns exceed even the maximum of the permutation distribution
(MI = 0.124 bits). Specificity is therefore reported as effect size, with the
permutation maximum quoted as a family-wise reference rather than used as a
selection criterion. Columns at which either group exceeded 20% gap content
were flagged and module means reported with and without them, separating
deletion from substitution; on this basis 5 of 12 KIR columns are
substitution-driven, at which the KIR still leads all modules (0.406 bits
against 0.311, 0.259 and 0.262).

## S6. Declared motif panel and multiple testing

The motif panel was fixed from the literature before any constraint or
specificity value was computed. Motif positions were located by literal string
search of the human canonical sequences rather than assumed, which corrected
one assignment: the "QR" and "FF" motifs lie in the BC loop, within the SH2
domain, not in the KIR. Searching gives GQM once in each kinase (JAK1
1097–1099, JAK2 1071–1073, TYK2 1117–1119, all within JH1 and all occupying the
same three alignment columns), FF at SOCS1 112–113 and SOCS3 79–80, QR at SOCS1
108–109 only and not positionally homologous in SOCS3, and FLVR at SOCS1
101–104 giving the invariant SH2 arginine Arg104.

The panel comprised eleven motifs, nine in SOCS (KIR_full 55–66; KIR_H54_R59
54–59, from Yasukawa et al. 1999; KIR_aromatic 64; ESS_full 67–78; ESS_I68_L75
68 and 75; BC-loop FF 112–113; BC-loop QR 108–109, SOCS1 only; SH2_FLVR
101–104; SH2 pTyr arginine 104) and two in JAK (GQM 1071–1073; activation loop
1004–1014, JAK2 numbering), giving 42 p-value-bearing tests across the
constraint and specificity measures. Two motifs serve as internal negative
controls: the SH2 phosphotyrosine arginine and the FLVR motif containing it are
required for phosphopeptide binding by all SH2 domains and so should show
constraint without inhibitor-specific discrimination, which is what is
observed (specificity 0.005-fold and 0.260-fold, q = 0.42 and 0.24).

Each motif was tested against size-matched resamples of non-motif columns from
the same domain, falling back to the whole block where a domain-restricted pool
was smaller than the motif plus two columns; the pool used is recorded per
test. Aromaticity was scored as the weighted fraction of Phe/Tyr/Trp with
Wilson 95% intervals, because the published mechanism is a Met–aromatic
interaction and the selected property is aromatic character rather than residue
identity. Benjamini–Hochberg correction was applied across all 42
p-value-bearing tests of the declared panel as one family, and separately
across columns within each specificity contrast.

## S7. Structural interfaces, the interface test, and depth

**Interfaces.** Coordinates for PDB 4GL9 (mouse JAK2 kinase domain, mouse
SOCS3, IL-6 receptor β fragment) and 6C7Y (human JAK1 kinase domain, chicken
SOCS1) were parsed with Biopython 1.88 `MMCIFParser`. Observed residues were
mapped to canonical numbering by global pairwise alignment of each chain's
observed sequence against the canonical sequence of the SIFTS-assigned
accession (`Bio.Align.PairwiseAligner`, BLASTP scoring, terminal gap penalties
zero), requiring ≥ 90% of observed residues to map at ≥ 95% identity.
Alignment-based mapping rather than SIFTS segment arithmetic was necessary:
4GL9 chain E is a ΔPEST construct spanning mouse SOCS3 22–128, 164–168 and
173–184, whereas the single SIFTS segment reports 38–128 and so omits the
entire KIR. Inter-chain contacts were defined at 8 Å minimum heavy-atom
distance, computed for all JAK × SOCS chain pairings and retained for pairings
with at least 10 contacts or 25% of the maximum, excluding lattice contacts
among the four crystallographic copies in 4GL9. This gave 46 SOCS3 and 42 JAK2
interface residues from 4GL9 and 43 SOCS1 and 51 JAK1 from 6C7Y. Per-residue
solvent-accessible surface area of the isolated chains was computed with
`Bio.PDB.SASA.ShrakeRupley`.

**Interface test.** Interface residues were compared against non-interface
controls matched simultaneously on domain membership, solvent-accessibility
tertile and alignment gap-fraction tertile, with the interface label permuted
*within* covariate strata (10,000 permutations) so that an SH2 interface
residue is compared only to other SH2 residues of similar burial and occupancy.
Strata containing only interface or only control residues cannot contribute and
were dropped, leaving 8–12 usable strata per test. Effect size was Cliff's δ
with 95% intervals from 5,000 bootstrap resamples, reported with Mann–Whitney U
and two-sample Kolmogorov–Smirnov tests. The four tests were combined by
inverse-variance fixed-effect meta-analysis on δ, with standard errors from the
bootstrap half-width and heterogeneity assessed by Cochran's Q and I²
(Q = 2.82, p = 0.42, I² = 0). Power was established by simulation on the
observed sample sizes: control values were resampled with replacement, shifted
by 0.1–1.0 control standard deviations and re-tested 600 times per shift,
placing the smallest shift detectable at 80% power at 0.50 control SD for the
JAK-side tests and 0.75 for the SOCS side — larger than any observed shift, so
large interface-specific effects are excluded while small-to-moderate ones
remain possible. Constraint is additionally near saturation on the kinase side,
with 70–80% of JAK control residues exceeding a divergence of 0.6.

**Depth.** Depth available for inter-protein coupling inference was assessed
against N_eff > (L₁+L₂)/2, with N_eff the sum of inverse neighbour counts at
80% sequence identity across all columns, gaps included as a symbol. Twelve
candidate designs were evaluated — both retrieval routes, four taxonomic tiers,
the six paralogue-homogeneous JAK × SOCS species pairings, a fully factorial
pool of all six, an invertebrate-inclusive pool and two column-set reductions —
and all fail, at ratios 0.037–0.586 on the scoping alignments; the design
carried forward to the coupling analysis reaches 0.665 on the final alignments
and also fails. N_eff at 70% and 90% thresholds is
reported for sensitivity only; all decisions used the 80% value. The limiting
factor is conservation rather than sampling: median pairwise identity in a
JAK2 × SOCS3 alignment of 365 species is 0.904 and 73.5% of sequence pairs
exceed the 80% threshold, so those species collapse to N_eff = 8.1.

Inter-block coupling was computed on the best-scoring design (factorial pool;
2,912 records; JH1 281 columns × SOCS modules 156 columns; N_eff = 145.3
against a requirement of 218.5) as mutual information with the average product
correction applied to the inter-block submatrix alone. Mutual information
rather than pseudolikelihood maximisation was used because at this depth
against 43,836 column pairs a pseudolikelihood fit is determined by its
regularisation prior. Applying the correction to the pooled two-block matrix
would subtract a background that is partly the block structure itself; because
the correction also forces the inter-block mean to approximately zero, the mean
statistic is reported but flagged uninterpretable and inference rests on the
tail. The null was 200 species-scrambled replicates in which the JAK–SOCS
species correspondence was randomised while both marginal alignments, all
conservation and gap structure, and each family's phylogeny were preserved.

## S8. Software and reproducibility

Python 3.12.14 with NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.6, Biopython 1.88,
statsmodels 0.15.0 and Matplotlib 3.11.2; MAFFT v7.526; MMseqs2 v18.8cc5c.
Every permutation, bootstrap and resampling step was seeded at 42. Analysis
code, intermediate alignments, coordinate maps, per-column provenance, all
result tables and a run log recording software versions, alignment dimensions
and the outcome of every validation assertion are deposited at [repository
URL], together with a single-command pipeline script that reproduces every
table and figure.

---

## Note on what belongs in the main-text Methods

Nothing above needs to appear in the main text beyond the pointer sentence. If
a reviewer asks for more in-text detail, the three items most worth promoting
are: (i) that paralogue assignment is by sequence against a decoy-containing
panel, since this is what excludes JAK3 and SOCS2–7; (ii) that SOCS2 serves as
the internal non-inhibitory contrast; and (iii) that the depth gate is a
pre-specified published criterion rather than a post-hoc judgement.
