# Methodological review: inter-protein direct-coupling analysis for JAK–SOCS

Scope: what contemporary practice requires of an inter-protein DCA, and which
parameter this project adopts at each decision point. Every recommendation is
traced to a retrieved source; full bibliographic detail, including the role each
reference plays in the design, is in [`references.csv`](references.csv).

This document is written to be read before the analysis, as a pre-registration
of method choices. It is deliberately explicit about the two places where the
previous attempt went wrong for methodological rather than clerical reasons,
because both are traps that the literature warns about and that a
reasonable-looking pipeline walks into by default.

---

## 1. Depth: what "enough sequences" means, and the gate this project adopts

Inter-protein DCA is substantially more data-hungry than the monomer problem.
The relevant published criterion comes from Ovchinnikov, Kamisetty & Baker
(eLife 2014, `10.7554/eLife.02030`), who examined the 50S ribosomal subunit plus
28 bacterial complexes of known structure and found that pseudolikelihood-derived
inter-protein covarying pairs are almost always genuine contacts **provided the
number of aligned sequences exceeds the average length of the two proteins**.
Written as a gate:

```
N_eff  >  (L1 + L2) / 2
```

The monomer convention — N_eff/L ≳ 1 as a bare minimum, ≳ 5 for comfort — is the
weaker of the two constraints and is reported alongside it, but the
Ovchinnikov criterion is the one that governs an inter-protein claim and is
therefore this project's primary gate.

For calibration against the previous run: `plmc` reported N_eff = 46.6 over
L = 483 with L1 = 280, L2 = 203. That gives N_eff/L = 0.097 and requires
N_eff > 241.5, so it missed the inter-protein gate by a factor of about five and
the comfortable monomer target by a factor of about fifty. No downstream
statistical treatment could have repaired that, which is why depth is attacked
first here rather than accepted as a constraint.

Two further quantities are tracked because they diagnose *where* depth is
missing rather than just how much: per-column effective occupancy (a column
present in few sequences contributes little regardless of global N_eff), and
N_eff/N, the redundancy penalty. The previous alignment's ~10× penalty is the
signature of a mammal- and bird-dominated sample, and it is the reason
taxonomic *breadth* matters more than taxonomic *volume* (§2).

## 2. Raising effective depth: which levers are real

Because N_eff for a paired alignment counts effectively independent
*(JAK, SOCS) species-pairs*, not sequences, levers divide sharply into those
that add independent information and those that only add rows.

**Defensible, ordered by expected yield.**

1. **Profile-based retrieval instead of gene-symbol queries.** A
   `gene:jak1`-style query returns only entries that a curator or automatic
   pipeline gave that symbol, which systematically excludes proteins from
   newly sequenced non-model genomes — precisely the divergent, low-redundancy
   sequences that raise N_eff. Retrieval by InterPro/Pfam family signature, or
   by `jackhmmer`/`hmmsearch` with an HMM built from a curated orthologue set,
   recovers them. Paralogue identity is then assigned downstream by
   best-match identity rather than trusted from the annotation.
2. **Taxonomic breadth, analysed as a nested ladder.** Outward from jawed
   vertebrates: cyclostomes (lamprey, hagfish), for which IL-6/STAT-family
   signalling has recently been described (Fontenla-Iglesias et al., bioRxiv
   2026, PMID 41993459), then invertebrate metazoans, where JAK (Drosophila
   *hopscotch*) and SOCS-box/SH2 proteins are present. Each added sequence at
   these depths contributes far more effective diversity than any additional
   mammal. **The cost is a change of question, not merely of noise:** SOCS1 and
   SOCS3 are gnathostome paralogues, so orthology is undefined outside
   gnathostomes and a metazoan-inclusive alignment addresses the JAK family
   against the SOCS-box/SH2 family. Tiers are therefore analysed separately,
   each with its own declared claim scope, and the tier ceiling is set on
   measured alignment quality rather than assumed.
3. **Reducing L rather than raising N_eff.** N_eff/L and the Ovchinnikov gate
   both improve as fast from halving L as from doubling N_eff, and both
   structures identify where contacts can physically occur. The subtlety that
   matters: restricting the model to structurally selected columns destroys the
   interface-enrichment test, because only the interface was modelled. This
   project therefore runs a full-domain model for enrichment and a reduced-L
   model for residue-level resolution, the latter paired with a size-matched,
   pre-declared non-interface control region.
4. **Cross-dataset pooling with paralogue identity as a fixed effect.**
   Combining evidence across the six paralogue-pair datasets recovers up to
   ~√6 in precision for genuinely shared signal while remaining blind to
   paralogue-specific artifacts — the statistically correct way to exploit the
   dataset structure rather than choosing between pooling and splitting.
5. **Paralogue recovery in duplicated lineages.** Teleosts and salmonids carry
   additional *socs1a/socs1b* and *socs3a/socs3b* copies (Jin et al. 2007
   PMID 18029016; Chu et al. 2019 PMID 31846779; Ma et al. 2025 PMID 40402276;
   Xue et al. 2024 PMID 38849106), which a one-per-species rule discards.
   Recovering them requires solving which copy pairs with which JAK — the
   standard paralogue-matching problem (§3) — so the gain comes with
   pairing-error risk and is kept as a sensitivity analysis.

**Refused, and why.** Isoforms and splice variants are not independent
evolutionary samples. Population-level variation is not at mutation–selection–drift
equilibrium and violates the generative assumption the model rests on. Raising
the reweighting threshold from 80% to 90% inflates the printed N_eff without
adding information, so thresholds other than 80% are reported only to display
estimator sensitivity and never used for a decision. Padding with unpaired
homologs of either family adds no cross-family information. Re-admitting JAK3,
or pooling paralogues without modelling composition, is the failure mode
diagnosed in the previous attempt.

## 3. Pairing strategy: orthology, species matching, and paralogue matching

Species matching alone is legitimate only when each species contributes one
unambiguous member of each family. When it does not, the choice is between
resolving paralogue identity in advance and inferring the pairing.

- **Inferring the pairing.** Gueudré et al. (PNAS 2016,
  `10.1073/pnas.1607570113`) introduced the Iterative Pairing Algorithm, which
  simultaneously matches interacting paralogues and infers inter-protein
  contacts. Bitbol (PLoS Comput Biol 2018, `10.1371/journal.pcbi.1006401`)
  showed an approximate mutual-information maximisation performs slightly
  better for partner identification specifically. Lupo, Sgarbossa & Bitbol
  (PNAS 2024, `10.1073/pnas.2311887121`) introduced DiffPALM, which uses
  MSA Transformer's masked-token objective and outperforms coevolution-based
  pairing **on shallow alignments** — the regime this project is in — while
  being competitive with orthology-based pairing.
- **Resolving identity in advance.** Where genuine one-to-one orthology exists,
  assigning it directly is more reliable than inferring it and leaves no
  pairing error to propagate.

**Adopted.** Gnathostome JAK1/JAK2/TYK2 and SOCS1/SOCS3 are one-to-one
orthologue groups, so identity is resolved in advance from best-match identity
against human canonical references, cross-checked against orthology-resource
calls. Six paralogue-homogeneous datasets follow, in which paralogue identity is
constant *by construction*. This is stronger than controlling for composition
after the fact: the confound is absent rather than adjusted for. Paralogue
matching (IPA/MI) is used only for the teleost duplicate-recovery sensitivity
analysis.

## 4. Phylogeny: a confounder for contacts, an asset for partners

This is the single most important methodological point for this project, and it
is the reason the previous attempt produced an ESS enrichment at *q* < 0.001
that was not real.

Marmier, Weigt & Bitbol (PLoS Comput Biol 2019,
`10.1371/journal.pcbi.1007179`) generated synthetic sequences containing
**phylogeny only — no interactions and no contacts** — and found that DCA still
accurately identifies which sequences share evolutionary history. Their
conclusion is explicit in both directions: phylogenetic correlations are
*useful* for predicting interaction partners, and they *confound* the
identification of contacting residues. Bitbol (2018) makes the complementary
point that the statistical dependences enabling partner prediction are not
restricted to residue pairs in contact at the interface.

The consequences adopted here:

- An inter-block coupling signal is **not** evidence of physical contact until
  it survives a null that reproduces the phylogenetic structure. A permutation
  null that treats pairs as exchangeable cannot do this.
- Paralogue composition is a phylogenetic covariate. When JAK and SOCS
  paralogue identities are non-independently distributed across records,
  correlated identity mimics correlated residue state — the mechanism behind
  the spurious prior result.
- Added taxonomic divergence helps for a second, independent reason beyond
  N_eff: it breaks up phylogenetic correlation. This motivates the
  divergence-ladder test, in which a genuine structural constraint is expected
  to strengthen with divergence while a phylogenetic artifact dissolves.
- Language-model-based evidence is not phylogeny-free either: Lupo, Sgarbossa
  & Bitbol (Nat Commun 2022, `10.1038/s41467-022-34032-y`) showed MSA-trained
  models learn phylogenetic relationships explicitly.

## 5. Concatenation, linkers, gap handling, and column provenance

**Linkers are unnecessary, and the previous concern was misplaced.** A
pairwise maximum-entropy model over alignment columns contains no
sequence-adjacency term; sites enter the graphical model exchangeably, so the
last JAK column and first SOCS column being numerically adjacent creates no
coupling by itself. Inserting a poly-N linker and then modelling it, as the
previous pipeline did, is actively harmful — it adds columns of near-constant
state that consume regularisation budget and appear in pair counts.

**Sequence separation must not be used as a cross-block covariate.** Intra-block
coupling magnitude falls with |i−j|, but inter-block pairs have no meaningful
separation: the numerical distance between a JAK column and a SOCS column is an
artifact of concatenation order. A separation-matched comparison *between*
intra- and inter-block classes is therefore ill-defined. Separation matching is
legitimate only *within* a block, and this is a narrower role than the previous
analysis gave it.

**Gap handling.** Following EVcouplings conventions (Hopf et al.,
Bioinformatics 2019, `10.1093/bioinformatics/bty862`): drop columns above a
declared gap fraction, drop sequences with excessive gap content, and — the
step that failed before — carry the coordinate map through the filtering so no
residue index is silently invalidated by renumbering.

**Column provenance is a hard requirement, not documentation.** Every final
column must carry its canonical residue identity, and the assertion that JAK
block columns map to kinase residues and SOCS block columns to SOCS residues
must be executed, not assumed. The previous pipeline's 203-column "SOCS core"
block contained no SOCS residues at all; an executed provenance assertion would
have caught it immediately.

## 6. Inference: estimator and parameters

**Estimator.** Asymmetric pseudolikelihood maximisation (plmDCA), the method
underlying both `plmc` and CCMpred (Seemayer, Gruber & Söding, Bioinformatics
2014, `10.1093/bioinformatics/btu500`), and the method whose inter-protein
behaviour Ovchinnikov et al. characterised. Mutual information with APC is
computed alongside as a model-free baseline, motivated by Bitbol's (2018)
finding that MI is competitive for cross-family signal.

**Sequence reweighting.** Each sequence weighted by the inverse of its number
of neighbours at 80% sequence identity; N_eff is the sum of weights.

**Regularisation.** L2 on fields and couplings, λ_h = 0.01 with λ_J scaled by
(L−1) in the `plmc`/EVcouplings convention. The exact resolved values are
written to `analysis.log` at runtime rather than restated here, so the record
cannot drift from the code.

**Gauge and score.** Zero-sum (Ising) gauge, then the Frobenius norm of each
20×20 coupling block as the raw score.

**APC, applied per submatrix.** The average product correction (Dunn, Wahl &
Gloor, Bioinformatics 2007, `10.1093/bioinformatics/btm604`) subtracts the
product of row and column means over the score matrix. Applied to a pooled
two-block matrix it mixes two populations with different score scales, so the
correction removes a background that is partly the block structure itself and
systematically distorts inter-block scores. APC is therefore computed
separately within the intra-JAK, intra-SOCS and inter-block submatrices, and
both raw and corrected scores are retained.

**Inter-protein scores are expected to be smaller than intra-protein scores.**
Hopf et al. (eLife 2014, `10.7554/eLife.03430`) normalise inter-protein
couplings against the intra-protein background for exactly this reason, and
Green et al. (Nat Commun 2021, `10.1038/s41467-021-21636-z`, EVcomplex2)
replace the ad-hoc normalisation with a probabilistic score. The previous
attempt's observation that cross-block mean score was below intra-block mean
score is a property of the estimator, not a finding about JAK–SOCS biology, and
the same holds for the scarcity of cross-block pairs in a global top-200 list.
Inter-block scores are ranked against inter-block scores here, never against a
pooled distribution.

## 7. Null models

Ranked by what each can rule out:

1. **Species scramble** — re-run the full inference on alignments in which the
   JAK–SOCS species correspondence is randomised while both marginal
   alignments are preserved. This destroys pairing while retaining each
   family's composition, conservation, gap structure and phylogeny. It is the
   reference null for inter-protein DCA and it is the only one that tests the
   actual hypothesis.
2. **Phylogeny-matched resampling** — following Marmier et al. (2019),
   quantify how much inter-block signal shared evolutionary history reproduces
   on its own.
3. **Paralogue-stratified permutation** — required for any pooled dataset; the
   control that collapsed the previous ESS result.
4. **Structural precision against experimental interfaces** — the most direct
   test available, since two independent structures exist (§8), with
   random-pair and gap-content-matched baselines.

**Rejected.** Fisher exact tests over pair counts assume pairs are
exchangeable. Entries of a coupling matrix are not: they share rows, columns,
and a common fitted parameter set. Such tests are reported for comparison only
and never as a primary result.

## 8. Structural ground truth

Two independent crystal structures cover two different paralogue pairs, which
makes cross-validation possible rather than merely calibration:

- **4GL9** — SOCS3 with the JAK2 kinase domain and an IL-6 receptor β fragment
  (Kershaw et al., Nat Struct Mol Biol 2013, `10.1038/nsmb.2519`). SOCS3
  engages receptor and JAK2 on opposing surfaces; JAK2 binding is
  phospho-independent and at a non-canonical surface, and the KIR occludes the
  JAK2 substrate-binding groove. Non-competitive kinetics are characterised in
  Babon et al. (Immunity 2012, `10.1016/j.immuni.2011.12.015`).
- **6C7Y** — SOCS1 with the JAK1 kinase domain (Liau et al., Nat Commun 2018,
  `10.1038/s41467-018-04013-1`). The KIR targets the JAK substrate-binding
  groove with high specificity; SOCS1's SOCS box is compromised for Cullin5
  recruitment, so direct kinase inhibition rather than ubiquitination is the
  operative mode.

Because 4GL9 is a JAK2–SOCS3 complex and 6C7Y a JAK1–SOCS1 complex, the
JAK2×SOCS3 dataset should preferentially recover 4GL9 contacts and the
JAK1×SOCS1 dataset 6C7Y contacts. A method that recovers both equally, or
recovers the wrong one better, is detecting something other than
pair-specific interface constraint.

## 9. Pre-declared motif hypotheses

Testing named motifs declared in advance is what separates a result from a
post-hoc pick over a large coupling matrix. The declared set:

- **SOCS KIR "Y", "QR" and "FF" motifs and the JAK "GQM" motif.** Zhai et al.
  (PCCP 2026, `10.1039/d6cp00092d`) compared six JAK/SOCS complexes by
  microsecond MD with MM/GBSA and attributed the stronger binding of
  SOCS1/SOCS3 to JAK1/TYK2 — relative to SOCS2 — to these KIR and BC-loop
  motifs forming hydrogen-bond networks and Met-aromatic interactions with the
  JAK GQM motifs. This is the most specific published residue-level hypothesis
  available for this system and it is adopted as the primary declared test.
- **SOCS1 KIR and extended SH2 subdomain.** Yasukawa et al. (EMBO J 1999,
  `10.1093/emboj/18.5.1309`) defined both regions and identified Ile68 and
  Leu75 in the ESS as conserved among JAB-related proteins and required for
  binding the phosphorylated activation loop; the KIR contributes to
  high-affinity binding of the JAK2 kinase domain and is required for
  inhibition.
- **SH2 phosphotyrosine pocket** (Bergamin et al., Structure 2006,
  `10.1016/j.str.2006.06.011`) and the **JAK substrate-binding groove and
  activation-loop tyrosines**, both defined from the structures in §8.

## 10. Sequence-only alternatives at low depth, and their evidence class

At marginal depth, methods that borrow statistical strength from the PDB
outperform DCA. They are useful here, but they answer a different question and
are labelled accordingly.

| Approach | Behaviour at low depth | Admissible as independent evolutionary evidence? |
|---|---|---|
| plmDCA (this project's primary) | Degrades sharply below the §1 gate | Yes |
| MI + APC | Degrades more gracefully; competitive for cross-family signal | Yes |
| DiffPALM / MSA Transformer (Lupo 2024) | Built for shallow MSAs | For *pairing*, yes; for contacts, partly — MSA-trained models encode phylogeny (Lupo 2022) |
| FilterDCA (Muscat et al. 2020, `10.1371/journal.pcbi.1007621`) | Structure-informed denoising of low-depth coevolution | No — supervised on structure |
| Paired-sequence PLMs (Liu et al., Nat Commun 2026, `10.1038/s41467-026-70457-5`) | Strong; trained for inter-protein contacts | No — supervised on structure |
| AlphaFold 3 (Abramson et al., Nature 2024, `10.1038/s41586-024-07487-w`) | Strong | No — circular against a structural ground truth |

The project's claim is whether evolutionary covariance *independently*
supports the JAK–SOCS interface. A predictor trained on the PDB cannot supply
that evidence and would be circular against the §8 validation. The unsupervised
methods (plmDCA, MI) therefore carry the primary analysis; structure-informed
methods can appear only as clearly-labelled convergent support.

For context on the achievable ceiling: Wuyun et al. (Brief Bioinform 2018,
`10.1093/bib/bbw106`) benchmark contact-prediction methods at scale, and
Morcos & Onuchic (Curr Opin Struct Biol 2019, `10.1016/j.sbi.2019.03.024`)
review where coevolutionary signatures do and do not support complex inference.
Mazzocato et al. (ACS Cent Sci 2024, `10.1021/acscentsci.4c01428`) show what a
validated coevolutionary signal can be used for downstream, which is the
motivation for demanding residue-level output.

---

## 11. Parameters adopted

| Decision | Value | Source |
|---|---|---|
| Primary depth gate | N_eff > (L1+L2)/2 | Ovchinnikov 2014 |
| Secondary depth diagnostic | N_eff/L, reported against ≳1 and ≳5 | monomer convention |
| Retrieval | InterPro family signature (IPR016251; IPR000980 + IPR001496), paralogue assigned by identity | §2 |
| Excluded paralogue | JAK3 (P52333), logged | user requirement |
| Human references | JAK1 P23458, JAK2 O60674, TYK2 P29597, SOCS1 O15524, SOCS3 O14543 | UniProt canonical |
| Taxonomic tiers | Gnathostomata → Vertebrata → Chordata → Metazoa, analysed separately | §2 |
| Pairing | Paralogue-resolved species matching; IPA/MI only for teleost duplicates | Gueudré 2016; Bitbol 2018 |
| Alignment | MAFFT per family, never on concatenated records | §5 |
| Linker | None; no linker columns modelled | §5 |
| Sequence reweighting | 80% identity | EVcouplings convention |
| Regularisation | L2, λ_h = 0.01, λ_J scaled by (L−1) | plmc / CCMpred convention |
| Gauge / score | Zero-sum gauge, Frobenius norm | standard |
| APC | Applied separately per submatrix | Dunn 2007; §6 |
| Inter-block ranking | Against inter-block scores only | Hopf 2014; Green 2021 |
| Primary null | Species scramble, 30–50 replicates | §7 |
| Structural validation | 4GL9, 6C7Y; 8 Å heavy-atom contacts | Kershaw 2013; Liau 2018 |
| Declared motifs | KIR Y/QR/FF, JAK GQM, SOCS1 ESS I68/L75, SH2 pTyr pocket | Zhai 2026; Yasukawa 1999 |
| Multiple testing | Benjamini–Hochberg across the whole declared family | standard |
| Seed | 42 | reproducibility |

## 12. What this design can and cannot conclude

It can establish, with quantified power: whether inter-block evolutionary
covariance concentrates in the KIR/ESS/SH2 modules above a species-scramble
null; whether that concentration is reproducible across independent
paralogue-pair datasets; whether it strengthens with taxonomic divergence as a
structural constraint should; and what precision the top couplings achieve
against two experimental interfaces.

It cannot establish direct atomic contacts, nor residue-pair coevolution
specific to a fixed JAK2–SOCS1 heterodimer, from sequence alone at this
family's diversity. Claims stay at the level the depth gate supports, and a
negative result from a design that passes the gate is reported as a negative
result — which, unlike the previous attempt's negative, would be informative.
