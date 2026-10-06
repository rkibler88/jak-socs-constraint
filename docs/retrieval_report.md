# Step 3 — Paralogue-resolved retrieval

Route: InterPro family signature (`IPR016251`; `IPR000980`+`IPR001496`), tier Metazoa
(superset — paralogue assignment restricts it). Runtime 2749 s.

## Yield against the previous pipeline's gene-symbol route

| Paralogue | gene-symbol species | InterPro + assignment | change |
|---|---|---|---|
| JAK1 | 398 | 617 | +55% |
| JAK2 | 391 | 648 | +66% |
| TYK2 | 362 | 366 | +1% |
| SOCS1 | 440 | 704 | +60% |
| SOCS3 | 455 | 704 | +55% |

TYK2 is the exception: it gains almost nothing, so TYK2-containing datasets stay
the shallowest. SOCS1 and SOCS3 landing on the same species count (704) is
coincidence, not a bug — the two sets share 662 species and differ in 84.

## Assignment outcome

14338 entries retrieved, 14336 searched against the reference panel,
**5899 accepted / 8437 rejected**.

Rejections are dominated by the decoy panel doing its job — the family
signature returns the whole SOCS family, and off-target members are assigned to
their own reference and excluded on evidence rather than missed:

| Best-matching reference | Entries diverted |
|---|---|
| SOCS7 | 1612 |
| SOCS5 | 1475 |
| SOCS6 | 1139 |
| CIS | 1078 |
| SOCS2 | 1000 |
| SOCS4 | 831 |
| JAK3 | 501 |

JAK3 exclusion is therefore enforced by sequence, not only by query construction.
SOCS2 is retained as an assignment target for the specificity analysis, since it
is the non-binding contrast in Zhai et al. (2026).

## Thresholds

Assignment required the best panel hit to be a target paralogue, with reference
coverage ≥ 0.55, identity ≥ 0.35, and a bitscore margin ≥ 1.10 over the best hit
from any *different* paralogue. Every sequence's identity, coverage, margin,
runner-up and rejection reason is recorded in `paralogue_assignment.csv`, so the
call is auditable per record.

One sequence per (paralogue, species), ranked by reference coverage then length —
not by length alone, which is how the previous pipeline preferred long fragments
over complete orthologues.
