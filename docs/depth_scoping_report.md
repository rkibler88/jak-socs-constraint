# Step 2 — Depth-maximisation scoping survey

Measured, not assumed. Every number below comes from `depth_scoping_raw.csv`
(66,148 UniProt entries retrieved across 2 routes x 4 taxonomic tiers),
`retrieval_yield.csv`, `per_pair_species.csv` and `depth_scoping.csv`.

## Headline

**No admissible design reaches the depth required for residue-level
inter-protein DCA.** The best achievable configuration reaches
**0.59x** the published requirement
(Ovchinnikov, Kamisetty & Baker, eLife 2014: `N_eff > (L1+L2)/2`). All
12 designs tested fail. The shortfall is a factor of
1.7 in the best case and up to
27 for a strictly orthologous pair.

This is not a data-collection failure. It is a property of the proteins.

## Lever 1 — retrieval route: large, real, and cheap

| Route | Paired species (Gnathostomata) | Paired species (Metazoa) |
|---|---|---|
| `gene:` symbol query (what the previous pipeline used) | 449 | 450 |
| InterPro family signature | **744** | **825** |

Profile-based retrieval recovers **+66% more paired species**. The mechanism is
confirmed rather than inferred: of the JAK-family entries in the 260 species
that InterPro finds and gene-symbol queries miss, **87% carry no primary gene
symbol at all** (84% on the SOCS side). These are proteins from
newly-sequenced, unannotated genomes — exactly the divergent sequences that
effective depth rewards.

InterPro signatures used, read off the human references rather than assumed:
`IPR016251` (non-receptor tyrosine kinase, Jak/Tyk2) and
`IPR000980 + IPR001496` (SH2 + SOCS box). SH2 alone would return every SH2
protein in the proteome; the conjunction is the SOCS family.

## Lever 2 — taxonomic breadth: small, and capped by annotation

Gnathostomata to Metazoa adds only **+81 paired species** (744 to 825), because
the JAK side does not extend: JAK-family species rise 747 to 830 while the SOCS
side rises 776 to 1,302. The limiting factor is annotation, not biology —
*Drosophila* Hopscotch (P08155) carries **no** `IPR016251` assignment and only a
single Pfam domain, so invertebrate JAKs are largely invisible to signature
retrieval.

But those 81 species matter far more than their count suggests (panel b):

| Set | Rows | N_eff | Efficiency (N_eff/N) |
|---|---|---|---|
| Gnathostome factorial pool | 2,045 | 79.6 | 3.9% |
| Invertebrate pairs (arthropods, mostly) | 81 | 45.7 | **56.4%** |

**81 arthropod species carry 57% of the effective depth of all 2,045 vertebrate
rows.** Vertebrate JAK/SOCS sequences are ~14x less informative per row. Their
functional columns remain occupied after alignment (invertebrate occupancy: JH1
0.89, SOCS KIR+ESS 0.86), so the added divergence is not purely gaps — though
whether a 12-residue KIR is *correctly* aligned across 700 Myr is a separate
question this survey cannot settle.

## Lever 3 — reducing L: helps, insufficient alone

Restricting the SOCS block from the full KIR-to-SOCS-box span (156 residues) to
KIR+ESS only (24 residues) moves the gate from 216 to 150. Combined with a
halved JAK block the gate falls to 81. Neither closes the gap.

## Lever 4 — factorial paralogue pooling

Pooling all six JAK x SOCS combinations raises N_eff from 8-20 (single
orthologous pair) to 79.6 — a ~6x gain, because JAK and SOCS paralogues are
genuinely divergent from one another. It is defensible here on biological
grounds: SOCS1 and SOCS3 both inhibit multiple JAK family members, so all six
pairings are real interactions rather than arbitrary concatenations. Paralogue
identities are independent by construction in a fully factorial design, which
removes the stratification confound that invalidated the previous ESS result.

**The N_eff of 79.6 is nonetheless optimistic**: each underlying sequence is
reused across pairings, so factorial rows are not independent observations. The
true count of independent paired species remains ~744.

## Why the depth cannot be raised further

The binding constraint is conservation. On the paralogue-homogeneous
JAK2 x SOCS3 alignment (365 species, JH1 + SOCS modules):

- median pairwise identity **0.904** (gap-inclusive), 0.944 over residue-only positions
- **73.5%** of all sequence pairs exceed the 80% redundancy threshold
- 365 species therefore collapse to N_eff = 8.1

DCA extracts information from *covarying variation*. A position that is
invariant across the alignment carries none, however many species are sampled.
JAK and SOCS are under strong purifying selection throughout jawed vertebrates,
so the sequence record that exists is the wrong shape for this estimator — and
the taxonomic range where they are divergent enough (invertebrates) is where
orthology is undefined, the KIR is unalignable with confidence, and only 81
paired species are annotated.

Note the irony in the previous attempt: its N_eff of 46.6 was *higher* than any
strictly orthologous design here (8-20) precisely because it mixed four JAK
paralogues with two SOCS paralogues. The paralogue mixing that invalidated its
statistics was also what inflated its depth.

## Design decisions taken from this survey

1. **Retrieval:** InterPro signature route, paralogue assigned downstream by
   identity. Adopted — a free 66% gain.
2. **Taxonomic ladder:** truncated. Vertebrata and Chordata add 2-8 species over
   Gnathostomata and are dropped as separate tiers. Gnathostomata is the
   orthology-resolved tier; a Metazoa-inclusive family-level tier is retained as
   the one configuration with meaningful added diversity.
3. **Paralogue-homogeneous per-pair DCA:** abandoned as the primary analysis.
   At N_eff 8-20 against a gate of 216 these datasets cannot support any
   residue-level claim, and reporting couplings from them would repeat the
   previous attempt's error with better bookkeeping.
4. **Residue-level inter-protein coupling:** not attainable. To be reported as a
   quantified negative with this power analysis as its evidence, rather than
   attempted and interpreted.

## Files

| File | Contents |
|---|---|
| `depth_scoping.csv` | all 12 designs with N, L1, L2, N_eff, N_eff/L, gate threshold and ratio |
| `retrieval_yield.csv` | entries and distinct species per route x tier x family |
| `per_pair_species.csv` | paired-species counts for the six paralogue combinations |
| `depth_scoping_raw.csv` | one row per retrieved UniProt entry |
| `depth_scoping_queries.json` | exact query strings and result counts |
| `fig_depth_scoping.{png,pdf,svg}` | panels a-c |
