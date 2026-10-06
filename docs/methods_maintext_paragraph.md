# Main-text Methods subsection

Replaces the evolutionary-analysis portion of the Methods. 171 words. Full
detail moves to Supplementary Methods 1–6.

---

## Drop-in text

**Evolutionary constraint and paralogue-specificity analysis.** JAK-family and
SOCS-family sequences were retrieved from UniProtKB by InterPro domain
signature (IPR016251 for JAKs; IPR000980 ∧ IPR001496 for SOCS proteins) across
Metazoa, and paralogue identity was assigned by MMseqs2 v18.8cc5c search against
a twelve-protein human reference panel including JAK3 and SOCS2–7 as decoys, so
that off-target family members are excluded on sequence evidence rather than
annotation. Each family was aligned separately with MAFFT v7.526; concatenated
records were never aligned. Per-column constraint was the gap-scaled
Jensen–Shannon divergence against Robinson & Robinson background frequencies
under Henikoff position-based sequence weighting, with bootstrap confidence
intervals; paralogue specificity was Miller–Madow-corrected mutual information
between residue state and group label, assessed against 10,000 label
permutations. Structurally defined interface residues (PDB 4GL9, 6C7Y) were
compared with controls matched on domain, solvent accessibility and alignment
occupancy, permuting within strata. Achievable depth for inter-protein coupling
was evaluated against the N_eff > (L₁+L₂)/2 criterion. All resampling was
seeded at 42. Full parameters are given in Supplementary Methods 1–6.
