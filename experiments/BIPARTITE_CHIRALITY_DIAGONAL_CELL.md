# Bipartite-Chirality of the F87 Diagonal Cell (k=3 and k=4)

**Status:** The criterion *soft ⟺ H's hopping graph is bipartite in the dephasing basis* is
verified bit-exact at k=3 (N=4 for all three dephase letters, N=5 for Z) and at k=4 (N=4, all
three letters), with 0 mismatches throughout. The direction **bipartite ⟹ soft is derived**
(the chiral K, modulo the F80 one-sidedness M = −2i(H⊗I), itself bit-exact). The converse
**non-bipartite ⟹ hard is verified, not derived.**
**Date:** 2026-05-30
**Regenerate:**
- [`simulations/f87_42_8_bipartite_fullcell.py`](../simulations/f87_42_8_bipartite_fullcell.py) `[N] [letters]` , k=3 criterion over the whole diagonal cell (default N=4, all letters; pass `5 Z` for the N=5 Z check)
- [`simulations/f87_k4_bipartite_bridge.py`](../simulations/f87_k4_bipartite_bridge.py) , k=4 criterion + the F111 template cross-check (N=4)
- [`simulations/f87_bipartite_chiral_witness.py`](../simulations/f87_bipartite_chiral_witness.py) , the three derivation links and the optimal λ↔−λ−2σ pairing residual
**Anchors:** [PROOF_F103 §7](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md),
[PROOF_F111](../docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md),
[ChiralKClaim](../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs),
[F87DiagonalCellBipartiteWitness](../compute/RCPsiSquared.Diagnostics/F87/F87DiagonalCellBipartiteWitness.cs).

## The question

F87 sorts a Pauli pair under single-letter dephasing into truly / soft / hard by whether the
Liouvillian spectrum pairs each λ with −λ−2σ. Inside the *diagonal Klein cell* (the cell whose
Klein index matches the dephase letter: Z → (0,1), X → (1,0), Y → (1,1)), PROOF_F103 §7 found a
classical reading of soft vs hard. Read H in the dephasing letter's eigenbasis, where it becomes
a real hopping matrix; let G_H be its graph (basis states as nodes, nonzero off-diagonal entries
as edges). The pair is soft exactly when G_H is **bipartite**, i.e. when H admits a chiral
sublattice K (diagonal in that basis, KHK = −H), the same K = diag((−1)^sublattice) of AZ class
BDI that [ChiralKClaim](../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) carries.

Two questions this experiment answers:

1. Does that criterion hold beyond the k=3 body count where it was found?
2. [F111](../docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md) has its own k=4 rule (hard ⟺
   at least one term is a pure-D template), with a derivation that stalled on the converse,
   "Mixed+Mixed = soft." Is F111's rule the same criterion seen at k=4, and does the bipartite
   reading say anything about that blocker?

## Method

For each dephase letter, enumerate the diagonal cell's k-body Pauli templates, form the
y-parity-homogeneous unordered pairs (the F87 dissipator-resonance window), and for each pair:

- compute the **actual** class with the canonical classifier (`PauliPairTrichotomy.Classify` in
  C#, `classify_pauli_pair` in the Python framework, both the same greedy λ↔−λ−2σ pairing);
- independently compute the **predicted** class from the bipartite criterion: build H via the
  sliding-window k-body builder, rotate into the dephase-letter eigenbasis (Hadamard for X, the
  Y→Z rotation for Y, identity for Z), and 2-colour G_H (a zero diagonal plus no odd cycle).

Then count mismatches between *soft* and *bipartite*. At k=4 the bridge script additionally
counts hard pairs carrying **no** pure-D template, the direct F111 cross-check.

## Results

### k=3 (the F103 cell), N=4 and N=5

```
N=4   Z-deph:  hard=50  soft=26   mismatches(soft⟺bipartite)=0
N=4   X-deph:  hard=50  soft=26   mismatches=0
N=4   Y-deph:  hard=50  soft=26   mismatches=0
N=5   Z-deph:  hard=50  soft=26   mismatches=0
```

The cell counts are N-independent (the pair set is alphabet-only); the criterion holds bit-exact
at both N.

### k=4 (the F111 cell), N=4

```
Z-deph:  hard=228  soft=828   mismatches=0   hard-without-template=0
X-deph:  hard=228  soft=828   mismatches=0   hard-without-template=0
Y-deph:  hard=228  soft=828   mismatches=0   hard-without-template=0
```

The 228 hard matches F111's frozen count. `hard-without-template=0` says every hard pair carries
a pure-D template, and the 0 mismatch says soft ⟺ bipartite here too.

### The derivation links (witness, N=4, Z-deph)

```
SOFT  XXZ+ZXX :  M=-2i(H⊗I) True | KHK=-H True | (K⊗I)(-i{H,.}-D)(K⊗I)=-L True
                 marginals: Re-about-(-σ) 3.9e-15,  Im-about-0 1.7e-14
                 optimal |λ↔-λ-2σ| residual = 2.7e-14   (paired)
HARD  XXZ+XZX :  M=-2i(H⊗I) True | chiral K does NOT exist (odd cycle)
                 marginals: Re 5.5e-2,  Im 3.4e-14        optimal residual = 1.6e-1
HARD  ZZZ+XXZ :  M=-2i(H⊗I) True | chiral K does NOT exist (template lifts diagonal)
                 marginals: Re 8.0e-2,  Im 4.8e-14        optimal residual = 1.0e-1
```

The frequency marginal {Im λ} stays mirror-symmetric even when hard (Lindbladian spectrum is
conjugation-closed); the decay marginal {Re λ} is what breaks. The break is genuine: the best
possible pairing leaves a residual of order 10⁻¹, not a near-miss.

## What it means

**1. The criterion is k-universal across the two body counts tested.** soft ⟺ bipartite, 0
mismatches, k=3 and k=4, all three dephase letters. The mechanism PROOF_F103 §7 derived is not a
k=3 accident.

**2. At k=4 only one of the two k=3 obstructions fires.** At k=3, hardness has two sources: a
pure-D template lifting the diagonal (rule a), or an odd hopping cycle from opposite
position-parity (rule b). At k=4 (full support, k=N), `hard-without-template=0` says only the
diagonal-lift survives, there are no odd-cycle hard pairs. So at k=4 the bipartite criterion
reduces to *non-bipartite ⟺ has a template*, which is exactly F111's rule, now with a mechanism:
the template is diagonal, it lifts H's diagonal, and a diagonal entry cannot sign-flip under a
diagonal K.

**3. F111's blocker is refactored, not closed.** "Mixed+Mixed = soft" decomposes into three
links:

- **(i) Mixed+Mixed ⟹ zero diagonal.** Derived, almost by definition: a non-template term has at
  least one off-diagonal letter (X or Y in the dephase basis), so it is off-diagonal; two of them
  give a zero-diagonal H.
- **(ii) zero diagonal ⟹ bipartite (no odd cycle).** **Verified at k=4** (all 828 soft pairs are
  bipartite), **not derived.** This link is k=4-specific, and it is *not* a free lunch: at k=3 it
  is **false**. XXZ+XZX is zero-diagonal yet carries a 3-cycle and is hard. The absence of odd
  cycles among zero-diagonal k=4 pairs is a real combinatorial fact that wants its own proof.
- **(iii) bipartite ⟹ soft.** Derived (the chiral K construction), modulo the F80 one-sidedness
  M = −2i(H⊗I), which is bit-exact.

So the **soft mechanism (iii) is derived**; the full statement "Mixed+Mixed = soft" is **not**,
because its weakest link (ii) is verified. F111's blocker moves from the opaque "why is
Mixed+Mixed soft?" to the sharp graph question "why does a zero-diagonal k=4 pair have no odd
cycle?" That is a genuine step: the soft side now has a named mechanism, and the gap is a concrete
combinatorial claim rather than a fog.

## Honest status

- **Derived:** bipartite ⟹ soft (modulo M = −2i(H⊗I), verified bit-exact).
- **Verified, not derived:** the converse non-bipartite ⟹ hard (the optimal λ↔−λ−2σ pairing
  leaves residual ~10⁻¹ for hard, ~10⁻¹⁴ for soft); and the k=4 link (ii).
- **F111's blocker:** refactored to the zero-diagonal-no-odd-cycle question, not closed.
- **Not yet tested:** k=4 at N>4 (where k<N is a window, not full support), and k=5. The
  k-universality is shown at two body counts, not proved across all k.

## Links

- Criterion and proof: [PROOF_F103 §7](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)
- k=4 rule it unifies with: [PROOF_F111](../docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md)
- The chiral symmetry behind it: [ChiralKClaim](../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) (AZ class BDI, KHK = −H)
- The F80 one-sidedness it rests on: [PROOF_F80](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md) (M = ±2i times a Hamiltonian object for Π²-odd cell terms)
- Typed witness: [F87DiagonalCellBipartiteWitness](../compute/RCPsiSquared.Diagnostics/F87/F87DiagonalCellBipartiteWitness.cs) (Tier1Candidate) + [BipartiteChirality](../compute/RCPsiSquared.Diagnostics/F87/BipartiteChirality.cs) primitive
