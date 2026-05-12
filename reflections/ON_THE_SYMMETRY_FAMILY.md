# On the Symmetry Family

**Status:** Reflection. Synthesis of the discrete-symmetry inventory of the chain XY+Z-dephasing Liouvillian after F91 (γ-Z₄) and the discovery that F71 alone is insufficient for the open questions 5/6/9 (factor-4 split, N=10 push, Z_n-Twins in J/h). Written 2026-05-12 alongside the typed `SymmetryFamilyInventory` aggregator.
**Authors:** Thomas Wicht, Claude (Opus 4.7)

---

For most of this project we have treated [F71](../docs/ANALYTICAL_FORMULAS.md) (chain spatial mirror) as if it were the symmetry. It is not. It is one of a family. The family has at least seven elements that we have already touched, and probably more we have not. The Elephant has been sitting in the room since the day we first wrote [`ChainMirror.cs`](../compute/RCPsiSquared.Core/Symmetry/ChainMirror.cs); F71 is necessary but not sufficient, and what F71 cannot show lives in axes orthogonal to F71's spatial-site action.

The five axes we now know:

The first is **U(1) × U(1)**: per-side popcount conservation. Continuous, not discrete. It is the substrate: every other symmetry refines or pairs the (N+1)² blocks U(1) × U(1) makes. We typed it as [`JointPopcountSectors`](../compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs). It is the floor; everything else stands on it.

The second is **F71**: chain spatial mirror. Discrete, Z₂, operator-side, splits each U(1) × U(1) block by factor 2 when the chain is uniform-γ-palindromic. Typed as [`F71MirrorBlockRefinement`](../compute/RCPsiSquared.Core/BlockSpectrum/F71MirrorBlockRefinement.cs). The day we discovered F71 is also the day the Elephant first appeared, but we did not see it then. F71 lives on the spatial axis alone; the parameter axis was invisible to us.

The third is **[F91](../docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md)**: F71 anti-palindromic in γ. Discrete, Z₄ (with i² = palindromic, i = anti-palindromic), parameter-side on γ, leaves the diagonal-block spectrum invariant. Tier 1 derived 2026-05-12 with full algebraic proof. We discovered it because we asked: what does F71 not see in the γ-distribution? The answer was: the pair-difference content. F91 is the parameter-side twin of F71's operator-side action.

The fourth and fifth are **[F92](../docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md)** (J anti-palindromic) and **[F93](../docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md)** (h anti-palindromic). The same Z₄ structure as F91 but for J_b in bond-coupling space and h_l in per-site detuning space. Mathematically identical to F91: the same proof structure, just applied to a different parameter family. They are the obvious extensions, and they are typed in this same plan.

The sixth is **Z⊗N**: Pauli-letter parity, the global Z-string operator. Operator-side, Z₂. We have it in [`Symmetry/ZGlobalMirror.cs`](../compute/RCPsiSquared.Core/Symmetry/ZGlobalMirror.cs) since the early days but never typed it as a BlockSpectrum refinement (the typed wrapper [`ZGlobalMirrorRefinement`](../compute/RCPsiSquared.Core/SymmetryFamily/ZGlobalMirrorRefinement.cs) landed 2026-05-12). It turns out to be trivially redundant with joint-popcount-parity, a pleasant null result that the inventory makes explicit. We type it anyway, because the inventory must be complete or it stops being an inventory.

The seventh is **[X⊗N](../compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs)**: global charge-conjugation. Operator-side, Z₂. It does not split blocks; it pairs them. Sector (p_c, p_r) and sector (N−p_c, N−p_r) have identical spectra. This halves the number of distinct eigendecompositions we have to run, even when no block-size reduction is available. It is the first of the family that gives a real algorithmic gain over F71-only without further math.

These seven sit on three axes: operator-side (F71, Z⊗N, X⊗N), parameter-side (F91, F92, F93), and the substrate (U(1) × U(1)). Negative results, [F71_col × F71_row factor-4 split](../compute/RCPsiSquared.Core/BlockSpectrum/F71BilateralBlockRefinement.cs), sit alongside as documented limits. K (chiral / AZ class BDI) sits at the edge: it commutes with the Hamiltonian but not with the dissipator, so it is not a full L-symmetry; it is typed in [`Symmetry/ChiralK.cs`](../compute/RCPsiSquared.Core/Symmetry/ChiralK.cs) and lives outside this family because [BlockSpectrum](../compute/RCPsiSquared.Core/BlockSpectrum/) cannot use it.

What this means for the work: every operator symmetry of L has a parameter-side twin that the F-chain inheritance generates automatically. F91/F92/F93 are not three separate discoveries; they are one structural fact, that the F71 spatial reflection has a Z₄ rotational extension in any parameter coordinate that controls L, applied three times. Future investigations of new dissipator types or new bilinear couplings will produce more such twins by construction. The plan's task list is finite; the structural identity it instantiates is reusable.

The reason we needed the Elephant before the open-questions could close is that Item 5 (factor-4 split via F71_col × F71_row) was a category error; we were trying to extract more from F71 alone when the answer is to widen to the family. Item 6 (N=10 push) is now visible as harder: F92/F93 do not give block-reduction; the path is X⊗N pairing plus possibly a new operator symmetry we have not yet seen. Item 9 (Z_n-Twins in J_b, h_l) is now solved by construction in F92/F93.

The Elephant is not a single new symmetry. It is the discipline of refusing to treat F71 as the only one. Once the [inventory](../docs/SYMMETRY_FAMILY_INVENTORY.md) is in place, the empty rows in the table are clearly empty, and the next discovery has a typed home before it arrives.

---

**Anchors:**
- [`docs/SYMMETRY_FAMILY_INVENTORY.md`](../docs/SYMMETRY_FAMILY_INVENTORY.md): typed inventory table
- [`compute/RCPsiSquared.Core/SymmetryFamily/SymmetryFamilyInventory.cs`](../compute/RCPsiSquared.Core/SymmetryFamily/SymmetryFamilyInventory.cs): Tier1Derived aggregator
- F71 math anchor: [`docs/ANALYTICAL_FORMULAS.md` §F71](../docs/ANALYTICAL_FORMULAS.md) + proof [`PROOF_C1_MIRROR_SYMMETRY.md`](../docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md)
- F91 algebraic proof: [`docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md`](../docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md)
- F92, F93 proofs: [`docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md`](../docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [`docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md`](../docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md)
- Pi2-Z₄ operator side: [`NinetyDegreeMirrorMemoryClaim` in `Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs)
- Companion reflection: [`reflections/ON_THE_NINETY_DEGREE_GAMMA.md`](./ON_THE_NINETY_DEGREE_GAMMA.md) (the operator-side / parameter-side Pi2-Z₄ recognition that opened the family)
- Negative result: [`compute/RCPsiSquared.Core/BlockSpectrum/F71BilateralBlockRefinement.cs`](../compute/RCPsiSquared.Core/BlockSpectrum/F71BilateralBlockRefinement.cs) (Tier2Empirical, factor-4 does not hold)
