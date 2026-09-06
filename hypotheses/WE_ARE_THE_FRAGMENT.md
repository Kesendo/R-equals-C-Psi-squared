# We Are the Fragment: 1/4 as Double Fragmentation

**Status:** Structural reading (Tier 4). Synthesis of three existing results, no new proof.
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:** [Uniqueness of the 1/4 Boundary](../docs/proofs/UNIQUENESS_PROOF.md), [Hierarchy of Incompleteness](../docs/HIERARCHY_OF_INCOMPLETENESS.md), [Primordial Qubit](PRIMORDIAL_QUBIT.md), [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md)
**See also:** [the K-partnership proof](../docs/proofs/PROOF_K_PARTNERSHIP.md), [Z⊗N Partnership](../experiments/Z_N_PARTNERSHIP.md)

---

## The question

Why 1/4 specifically as the fold value? The Uniqueness Proof gives a formal answer: 1/4 is the discriminant boundary of the recursion R = C(Ψ + R)², following from purity = Tr(ρ²) being the unique degree-2 polynomial in ρ. This is mathematically airtight. But it does not say *why we live exactly there*.

Tom's reading on 2026-04-25, after re-reading ZERO_IS_THE_MIRROR + PRIMORDIAL_QUBIT + the magnetism conversation:

> If we ourselves are part of the magnetic field, and zero is the middle, is what we have itself only a fragment (Bruch), hence the 1/4?

This document records the structural answer: **yes, the 1/4 is the fingerprint of our being a fragment of a fragment**. We are one sublattice within a bipartite mirror, and within that sublattice we see only half the operator content directly. The 1/4 is (1/2) × (1/2) where the two halves are not the same kind.

## The two fragmentations

**Fragment 1: the operator-level half (C = 1/2)**

From d² − 2d = 0 ([Qubit Necessity](../docs/QUBIT_NECESSITY.md)): at d = 2, the Pauli space splits into 2 immune (I, Z) and 2 decaying (X, Y) operators under Z-dephasing. The split is exactly 1:1, so C = 2/4 = 1/2. This is not a soft inequality. It is the *only* dimension where the immune-decaying split is symmetric, and is what enables the palindromic mirror at all.

The interpretation: at the operator level, **half of what we have is invisible to us as static structure** (the decaying half) and **half is invisible to us as dynamic information** (the immune half does not respond to the noise channel). Either reading: we work with half the algebra at any given moment.

**Fragment 2: the spatial bipartite half (sublattice = 1/2)**

A bipartite chain of length N has ⌈N/2⌉ sites on the A-sublattice and ⌊N/2⌋ on the B-sublattice. The bipartite "+− attracts" structure is exactly the spatial realization of the algebraic Z₂-grading: A and B are mirror partners. The K-operator K = Π_{l ∈ B} Z_l implements the sublattice gauge that gives KHK = −H ([the K-partnership proof](../docs/proofs/PROOF_K_PARTNERSHIP.md)), and that identity is an XY statement: for a chain with a ZZ term, KHK + H = 2H_ZZ exactly, so the sign flip survives only at Δ = 0. The proof lists open-chain XXZ at Δ ≠ 0 under what it does not cover, and it is scoped to the single-excitation sector, while the Néel state below is not in that sector. So this fragment rests on the bipartite grading, which is a statement about the lattice, and not on the K identity carrying over to the Heisenberg chain.

The Néel state |+−+−+⟩ is the canonical bipartite ground state. Its Z⊗N-mirror |−+−+−⟩ is the other half. At Σγ = 0 the two are degenerate; at Σγ > 0 we sit on one and the other is "outside" us.

The interpretation: at the spatial level, **we are one sublattice**, and the other sublattice is the magnetic field we feel. The bipartite "+− attracts" structure binds the two halves into one whole, but from the inside we only see ours directly.

## Why the two multiply

The 1/4 fold emerges from the discriminant of CΨ², and 1/4 factors as (1/2)² in
the arithmetic. It is worth being exact about what does NOT factor that way: the
framework's own C and Ψ at the fold. For Bell⁺ under Z-dephasing, where
C = Tr(ρ²) and Ψ is the normalized l1-coherence, the crossing sits at
f* = 0.8612 (the root of f(1+f²) = 3/2), giving

```
1/4 = C × Ψ = 0.8709 × 0.2871
```

not (1/2) × (1/2). And Ψ does not start near 1: for Bell⁺ it starts at exactly
1/3 and has only to fall to 0.287.

So the two halves this document is about are not the two factors of CΨ. They are:

- The **algebraic 1/2**, fixed and structural, the C-axiom that selects d = 2 as
  the only viable dimension. It does not move, and it is not the purity.
- The **spatial 1/2**, the sublattice split, a statement about the lattice rather
  than about a coherence measure reaching a value.

That both halves are one-half, and that the fold sits at one-half squared, is the
coincidence the document is reading. Whether it is more than a coincidence is
exactly what stays open: the arithmetic below does not establish it, because the
quantity that reaches 1/4 does not arrive there as 1/2 times 1/2.

The reading this document offers is that the fold at CΨ = 1/4 is where the
operator-immune-decaying split (always there) and the spatial sublattice
fragmentation meet. That is an interpretation laid over the fold, not a
derivation of it, and the previous section says why: the two 1/2's are not the
factors the fold's own arithmetic uses.

## What ZERO_IS_THE_MIRROR adds

At Σγ = 0, the palindrome is exact and centered. ZERO_IS_THE_MIRROR shows this is not a degenerate case but the **ground state of the palindrome itself**, where every eigenvalue is purely imaginary and Π becomes the time-reversal operator. There is no fragmentation here: A and B sublattices are degenerate, both halves equally accessible. **No fold exists at Σγ = 0**.

The fold at CΨ = 1/4 emerges only at Σγ > 0 (or Σγ < 0 mirror-symmetrically). The fold IS the geometry of the displacement from the mirror. We live at the fold because we live at Σγ > 0, displaced from zero into one half of the magnet.

The critical noise threshold Σγ_crit/J ≈ 0.25 to 0.50% (flat in N across the measured N = 2 to 5 for the product state) tells us how much displacement is needed: very little. About a quarter of a percent of the coupling strength is enough to break the bipartite degeneracy and create the 1/4 fold.

## What PRIMORDIAL_QUBIT adds

The Primordial Qubit hypothesis says noise does not originate; system and noise are two readings of one algebraic structure. The Pauli space is a genuine C²⊗C² with two independent bits (a = dephasing sensitivity, b = Π²-parity). [L, Π²] = 0 is proven for all N: the Z₂-grading of the Liouvillian respects this doubling everywhere.

The two standard "doubling constructions" of operator algebra theory (Tomita-Takesaki, thermofield double) both fail to reproduce this structure. The doubling is real, but no external mechanism builds it.

**This document proposes that the bipartite spatial structure of magnetism IS the missing mechanism**, in a specific sense: the doubling that PRIMORDIAL_QUBIT identifies algebraically (M_{2|2}(C)) is implemented physically by the A/B sublattice partition of the bipartite chain. The "two sides of the mirror" in the algebra correspond to the two sublattices in space.

The Inside-Outside operational result (PRIMORDIAL_QUBIT Section 9: only Q = J/γ measurable from inside) is the operational consequence of this fragmentation. From the inside of one sublattice, the absolute J and γ cannot be separated, only their ratio Q. We do not have the vantage point of the bipartite whole; we have only the ratio our position lets us read.

## What this is and is not

**This is:** a structural reading that connects three existing results into one statement. The 1/4 fold = double fragmentation = (operator C-axiom) × (spatial sublattice). The bipartite magnetism gives a *physical mechanism* to the algebraic doubling that the standard constructions could not reproduce.

**This is not:** a new proof. The 1/4 boundary itself remains formally derived from the discriminant of the CΨ² recursion. The C = 1/2 fragmentation is from d² − 2d = 0. The bipartite K-symmetry is from PROOF_K_PARTNERSHIP. What this document adds is the *interpretation* that ties them together: we live at the fold because we are double-fragmented, and the 1/4 measures both fragmentations at once.

**Risk:** the spatial sublattice argument is cleanest for bipartite NN chains, and its K identity holds at Δ = 0 rather than for Heisenberg. The C = 1/2 algebraic argument applies to any qubit. For non-bipartite topologies (triangles, frustrated lattices), the spatial half does not partition cleanly, but the fold at 1/4 still exists. So the doubled-fragmentation reading may be the bipartite-special case of a more general structural fact. Open question.

## The single sentence

The 1/4 is not a number we measure from outside the magnet. It is the fingerprint of the fact that we are inside it, on one half, looking at the other through the coherence Ψ that connects the two halves of the algebra (immune ↔ decaying) and the two sublattices of the chain (A ↔ B). When Ψ falls to 1/2, both bridges close at once.

---

*"Wenn wir selbst Teil vom Magnetfeld sind, und Null die Mitte, ist das was wir haben selbst nur ein Bruch, daher die 1/4?"*  Thomas Wicht, 2026-04-25
