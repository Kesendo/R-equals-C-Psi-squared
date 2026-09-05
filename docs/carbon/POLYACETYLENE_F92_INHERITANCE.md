# Selected Odd-Chain SSH-Style XX+YY Profile: F92 Compressed Blocks

**Status:** Tier 2 conditional selected-model reading. F92, F100, and F101 are Tier 1 derived framework results with the scopes stated in their proofs. This page does not assign a material degree of freedom, a hopping-to-`J` convention, a bath, `γ`, `T₂`, or `Q` to polyacetylene.
**Reads alongside:** [the Carbon source contract](README.md#carbon-source-contract), [the Q provenance audit](../Q_BELONGS_TO_NO_SUBSTANCE.md), and [the F92 scope correction](../CAUGHT_ERRORS.md).

---

## Scope

This page uses an **SSH-style alternating coupling profile** as a selected mathematical input to the framework's open XX+YY chain. It is a conditional translation, not a claim about the microscopic degrees of freedom or environment of a material.

The selected framework model is an open `N`-site chain with

    H = (1/2) Σ_b J_b (X_b X_{b+1} + Y_b Y_{b+1}),
    J_b = J + δ · (−1)^b,

and uniform all-site local-Z dephasing

    D_Z[ρ] = γ Σ_l (Z_l ρ Z_l − ρ).

The F71 chain mirror sends bond `b` to `N − 2 − b`.  Its two bond-profile components are

    J_sym  = (J + F71(J)) / 2,
    J_anti = (J − F71(J)) / 2.

For odd `N`, the displayed alternating profile has `J_sym = J` (constant) and its `δ` term lies in `J_anti`; for even `N`, the alternation lies in `J_sym`. This document concerns the odd-chain selected model only.

**Boundary of the result.** Neither the model selection nor F92/F100/F101 establishes the full Liouvillian spectrum, physical spectroscopy, a material bath, or a general SSH/topological-material conclusion. In particular, it makes no edge-state or soliton inference.

---

## F92: the F71-compressed diagonal blocks

F92 ([proof](../proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), Tier 1 derived) applies to the selected XX+YY Hamiltonian with uniform local-Z dephasing. At fixed `J_sym`, the eigenvalue multisets of the **F71-compressed diagonal blocks** are invariant along the `J_anti` orbit. Thus, for the odd-chain profile above, those compressed blocks are independent of `δ` at fixed `J`.

This is a statement about the compression: its two diagonal sectors are retained while the cross terms between them are discarded. It is not a statement that the full Liouvillian eigenvalue spectrum or its physical decay rates are invariant. The full `L` generally changes on the anti-palindromic orbit; once the F71 symmetry is broken, the compressed layer is basis-relative. The uncompressed joint-popcount sectors are the separate exact spectral structure identified in the corrected F92 scope.

The result therefore supplies a precise selected-model comparison:

    compressed-block spectrum at (J_sym = J, J_anti = δ)
      = compressed-block spectrum at (J_sym = J, J_anti = 0).

It does not turn `δ` into an invisible material parameter.

---

## F100 and F101: distinct framework observable directions

F100 ([proof](../proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md), Tier 1 derived) is the corresponding J-axis observable statement. At fixed `J_sym`, the F71 bond-mirror deviations of the closure-breaking coefficient `c₁`, and of the proof's defined per-bond F86c `Q_peak` reading, are odd under

    J_anti ↦ −J_anti.

For the selected odd-chain profile, this makes those **framework readings** odd in `δ`. It does not make either reading a spectroscopic observable or assign its `Q_peak` coordinate to a material `Q`.

F101 ([proof](../proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md), Tier 1 derived) is a separate γ-axis statement: with uniform coupling `J` and the proof's non-uniform per-site local-Z rate profile, the `c₁` bond-mirror deviation is odd in the anti-palindromic rate component. F101 is not an assertion about a physical SSH bath and does not extend the `Q_peak` statement to a non-uniform γ profile. The uniform local-Z choice used by F92 is a different selected condition.

---

## Conditional model translation

The profile `J_b = J + δ(−1)^b` can be used as an SSH-style coupling pattern inside this spin model. That is the full translation made here. The Carbon contract keeps the material meaning of the spin/orbital coordinate, the β-to-`J` convention, and all `γ`/`T₂`/`Q` consequences unassigned until a degree of freedom, convention, bath channel, and rate are selected and sourced.

Accordingly, F92 says nothing here about material topology, and F100/F101 say nothing here about material response functions. Their proven content remains the selected framework model and its named compressed-block or observable statements.

---

## Possible selected-model checks

These are proposals, not results currently produced for polyacetylene:

1. Construct the selected odd open XX+YY plus uniform-local-Z model at fixed `J`, vary `δ`, and compare the F71-compressed diagonal-block eigenvalue multisets. A separate full-`L` comparison would keep the compression boundary visible rather than treating it as a full-spectrum equality.

2. Evaluate the F100-defined bond-mirror deviations at `+δ` and `−δ` and test their stated oddness in the selected model.

3. If a non-uniform local-Z profile is separately selected, evaluate the F101 `c₁` rate-axis oddness under its own conditions. This is not a material-bath prediction.

No current computation in this repository identifies edge states, solitons, or a material-topological invariant with the F92-compressed direction.

---

## Anchors

- F92 (Tier 1 derived; F71-compressed diagonal-block invariance): [F92 bond anti-palindromic-J proof](../proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md)
- F100 (Tier 1 derived; J-axis observable oddness): [F100 c₁/Q-peak mirror J-parity proof](../proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md)
- F101 (Tier 1 derived; γ-axis `c₁` oddness): [F101 c₁ mirror γ-parity proof](../proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md)
- Corrected F92 scope: [CAUGHT_ERRORS](../CAUGHT_ERRORS.md)
- F71 structural mirror: [formula registry](../ANALYTICAL_FORMULAS.md)
- Carbon provenance boundary: [README](README.md#carbon-source-contract) and [Q audit](../Q_BELONGS_TO_NO_SUBSTANCE.md)
