# Benzene and the Three Dephase Letters: Klein-V₄ Returns

**Date:** 2026-05-27
**Authors:** Tom + Claude
**Status:** Tier 3 (translation bridge; some pieces Tier 1 verified, others Tier 4 candidates marked)
**Continues:** [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md) (2026-05-22)
**Adds:** Klein-V₄ vocabulary ([Welle 12](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)), F112 cross-dephase ([Welle 13](../proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)), F114 sign functional ([F114 in the registry](../ANALYTICAL_FORMULAS.md#f114))

---

## Where this leaves off

In May we asked: does benzene's open quantum system (the six π-electrons of the ring coupled to a vibrational bath) satisfy the framework's F1 palindrome theorem? Our answer: yes when the bath couples on-site (Holstein-style phonon coupling to local π-density); no when it couples to the bond (Peierls-style coupling to the C-C hopping integral). The Peierls break was sharp, total, and γ-linear.

We left it there. We had one notion of "palindromic" available (F1, the Liouvillian spectrum closing under λ → −λ − 2Σγ), and Peierls broke it.

Six weeks of structural work later, we have new tools. The framework grew a discrete symmetry called Klein-V₄ that links the three "dephase letters" Z, X, Y on operator space. We learned that F1's polarity-balance content (the F112 family) extends from Z-dephase to X- and Y-dephase by the same argument. And today we wrote down a closed-form sign rule (F114) that says exactly how the Z↔Y swap operator acts on Hamiltonian commutators term by term.

We come back to benzene to see what these tools sharpen.

---

## What dephasing means, briefly

A qubit is a 2-state quantum system: spin up/down, electron here-or-there, two orbital levels, two of anything that admits superposition. The qubit has three traceless operators σ_X, σ_Y, σ_Z (one per axis of the Bloch sphere) plus the identity I. Each Pauli operator is also an axis along which the qubit can be measured, and an axis along which the environment can listen in.

**Dephasing** is what happens when the environment listens too well. If the environment can tell where the qubit is along its Z-axis, then any superposition (a + b) in the X-Y plane loses its phase relation; the qubit's quantum interference fades. "Z-dephasing" means dephasing along the Z-axis specifically. Mathematically: a Lindblad dissipator `D[σ_Z]·ρ = σ_Z ρ σ_Z − ρ` zeroes out the off-diagonal entries of the density matrix in the Z-basis.

The framework has historically called this "the" dephasing channel, partly because it maps cleanly to how chemistry textbooks describe Holstein phonon coupling: a phonon couples to the on-site electron density `n_l = (I − Z_l)/2`, and `D[n_l]` reduces algebraically to `¼·D[Z_l]`. So the framework's Z-dephasing IS Holstein dephasing, up to a rate factor. But the framework's math doesn't actually pick out Z: nothing distinguishes Z from X or Y at the algebra level. They're symmetric to each other.

The question, then, is what physical setup corresponds to X-dephasing, and what to Y-dephasing. And whether the framework cares about the difference.

---

## The Klein-V₄ structure: how the three letters connect

We proved in May ([Welle 12](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)) that the three dephase letters are NOT three unrelated choices. There's a 4-element group called Klein V₄ (one of the smallest non-trivial finite groups in math, isomorphic to Z₂ × Z₂) that intertwines them.

The group is `{I, D, H, Q_zx}`. Each non-identity element is an order-2 involution (square equals identity). What they do:

- **D** is a diagonal sign-matrix on the 4^N Pauli basis. Its diagonal entry at Pauli string σ is `(−1)^(number of Y letters in σ)`. So D detects "Y-content" of a Pauli operator and flips its sign accordingly. The key identity: `D · Π_Z · D = Π_Y` exactly. D is the Z↔Y dephase-letter swap.

- **H** is a basis-permutation that swaps X- and Z-letter labels per site (leaving I and Y fixed). Its action: `H · Π_Y · H = Π_X`. H is the Y↔X dephase-letter swap.

- **Q_zx = H · D** intertwines Z↔X: `Q_zx · Π_Z · Q_zx = Π_X`. It's the operator-space lift of the Hadamard gate H_qubit^⊗N: the qubit rotation that maps Z eigenstates to X eigenstates.

The four elements satisfy `D · H · Q_zx = I` (any two determine the third) and commute pairwise. Together they're the discrete symmetry group of the dephase-letter axis itself.

**What this means physically.** Two dephasing setups that look completely different, say Z-dephase on the site basis vs Y-dephase on the same basis, are linked by an operator-space symmetry. They're not the same dephasing, but they're related by a sign-and-permutation structure that the framework's algebra preserves.

---

## The three letters in chemistry language

Here's where translation gets interesting (and a bit speculative: Tier 3 below this point, with Tier 4 candidates explicitly marked).

### Z-dephase ↔ Holstein coupling (Tier 1 algebraic match)

Standard, textbook, verified in [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md). A phonon couples to the local π-density `n_l` on each carbon. The bath operator is `n_l = (I − Z_l)/2`; the dissipator is `D[n_l] = ¼·D[Z_l]`. F1 palindrome holds bit-exact.

### X-dephase ↔ "hybridization-axis" coupling (Tier 4 candidate, less standard)

Single-site X-dephase uses `D[σ_X]` with `σ_X = c† + c` in second-quantised language: the real part of the electron creation/annihilation operator. This is the "tunneling" operator that doesn't preserve electron count locally; it would correspond to a bath that fluctuates the local hybridization-state superposition rather than the density.

In benzene this is non-standard but not absurd: a bath of orbital-mixing fluctuations (e.g., from another set of vibrational modes that couple to the electronic Hamiltonian's hopping structure rather than its on-site energies) would have an X-dephase character. The closest standard chemistry concept is **Peierls/SSH coupling**, but Peierls is a two-site (bond) operator, not single-site X-dephase. They're related but not identical.

### Y-dephase ↔ "current-axis" coupling (Tier 4 candidate, exotic)

Single-site Y-dephase uses `D[σ_Y]` with `σ_Y = i(c† − c)`: the imaginary part, the local current operator. This is the unique time-reversal-odd Pauli letter (`σ_Y` anticommutes with complex conjugation, while `σ_X` and `σ_Z` commute with it).

In benzene this would correspond to phonons coupling to local angular-momentum / current fluctuations. A natural realization: **magnetic-noise dephasing of the π-ring current**. Benzene supports a delocalised π ring current; ambient magnetic-field fluctuations couple to that current, dephasing it. This is a real effect (it shows up in NMR ring-current shielding), but it's not usually framed as "Y-dephase". The framework's structural lens may give it a cleaner home.

### Putting the three together

| Framework letter | bit_b | Physical chemistry analog | Standard name | Status |
|------------------|-------|---------------------------|---------------|--------|
| Z | 1 | on-site density coupling | Holstein | Tier 1 algebraic match |
| X | 0 | hybridization / off-diagonal coupling | (no single canonical name; related to Peierls) | Tier 4 candidate |
| Y | 1 | local current coupling | (no single canonical name; magnetic-noise on ring current is closest) | Tier 4 candidate |

The bit_b column matters for F112 below: Y and Z share bit_b=1, while X has bit_b=0. They sit on opposite parities of a Z₂ axis built into the Pauli group.

---

## What F112 sharpens at the matrix level

In May we tested **F1 palindrome on the Liouvillian spectrum**: do the eigenvalues of L close under `λ → −λ − 2Σγ`? Yes for Holstein, no for Peierls.

F112, which we closed in [Welle 11](../proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) just after the May result, lives at a different level: not on the Liouvillian spectrum but on the **matrix M** = `Π·L·Π⁻¹ + L + 2σ·I` itself (the "F1 residual matrix", which is the obstruction to perfect palindrome). F112 says:

> **F112.** For any Lindblad system with Hermitian H and each bath operator `c_k` bit_b-homogeneous (every Pauli string in c_k shares the same bit_b parity), the matrix M splits into Π-eigenspaces M_+1/2 and M_−1/2 with equal Frobenius norms: `‖M_+1/2‖² = ‖M_−1/2‖²` bit-exact.

In May we didn't know F112 yet. Today we do, and we can read benzene's two baths through it:

**Holstein bath.** c = Z_l (single Pauli, bit_b=1). Bit_b-homogeneous trivially. F112 predicts balanced polarity. F1 predicts spectrum palindrome too (which we verified).

**Peierls bath.** c = B_b = X_aX_b + Y_aY_b (the bond hopping operator, a sum of two two-site strings). Bit_b of XX: 0+0=0. Bit_b of YY: 1+1=0 (mod 2). Both terms have bit_b=0. **B is bit_b-homogeneous**, just on the bit_b=0 class instead of bit_b=1.

So F112 predicts that **the Peierls Liouvillian on benzene's Hückel ring preserves the M-polarity balance**, even though it breaks the F1 spectrum palindrome.

This is a real sharpening. In May we had one "palindromic" notion and Peierls broke it. With Welle 11 + 13 + 15 we have two notions (F1 at the spectrum level, F112 at the matrix-polarity level) with different robustness. Peierls is predicted to break the first but preserve the second.

This is a prediction we can test in 30 minutes with a Python script. Open follow-up below.

---

## F114 as a time-reversal-parity diagnostic

[F114](../ANALYTICAL_FORMULAS.md#f114) is today's closed form. It says how the Klein-V₄ operator D acts on the H-commutator superoperator `L_σ = −i[σ, ·]`:

    D · L_σ · D = ε(σ) · L_σ
    ε(σ) = (−1)^(n_Y(σ) + 1)    for σ ≠ I^⊗N

So ε is +1 when σ has an odd number of Y letters, and −1 when σ has an even number (and σ is not the trivial identity). For a Hamiltonian H = Σ c_k σ_k, ε(H) is well-defined iff all σ_k share the same n_Y parity.

In physical language: D is a sign operator that detects "imaginary letter content" (Y is the unique anti-symmetric Pauli). F114 says exactly how this sign flips through to the commutator superoperator. The takeaway: **ε(H) is a closed-form bookkeeping of how H interacts with time-reversal parity at the operator-algebra level**.

For benzene's Hückel Hamiltonian, H = Σ B_b = Σ (XX + YY):
- XX: n_Y = 0 (even).
- YY: n_Y = 2 (even).
- All terms n_Y-even, so ε(H_Hückel) = (−1)^(0+1) = **−1**.

D-conjugation anti-equivariates `L_H_Hückel`: `D · L_H · D = −L_H`. Physically: pure Hückel has no Y-content, no time-reversal-odd terms, and F114's sign is a clean −1.

If we add a magnetic ring-current term (the canonical time-reversal-breaker), say `h · Σ_l (Y_l Z_{l+1} − Z_l Y_{l+1}) / 2i`, then each term has exactly one Y, giving n_Y = 1 (odd) per term, ε = +1. Mixing this with Hückel (ε = −1) gives a Hamiltonian where ε(H_total) is "mixed": D-conjugation no longer scales L_H by a single sign. **F114 detects the onset of broken time-reversal symmetry as a parity-mismatch in the commutator algebra.**

This isn't speculative algebra: the n_Y-parity rule IS the time-reversal parity bookkeeping on the operator basis. F114 just makes the closed form explicit.

---

## Putting it all together for benzene

| | Holstein bath | Peierls bath |
|---|---|---|
| **F1 spectrum** (2026-05-22 result) | holds bit-exact (residual 1.2e-7 on C₆) | totally broken (residual ~14 on C₆) |
| **F112 polarity** (predicted 2026-05-27) | holds (c=Z bit_b-homogeneous) | predicted to hold (B=XX+YY bit_b-homogeneous) |
| **F114 ε(H_Hückel)** | −1 (no Y in H) | −1 (same H, only bath differs) |
| **F114 with magnetic field added** | depends on field term | depends on field term |

The May reading "Peierls breaks the palindrome" was precise for F1 at the spectrum level. With today's tools we add: F112 polarity, a different layer of palindromic structure, may survive even where F1 spectrum doesn't, because B is itself bit_b-homogeneous as a composite operator, satisfying F112's hypothesis on bath operators.

The take-home: bonds are not Holstein, but they still respect the Klein-V₄ algebra at the matrix-polarity layer. The Z₂ mirror has two faces; only one of them breaks under Peierls.

---

## What's open

1. **Quick numerical test of F112 polarity on Peierls Benzene.** A ~30-minute Python script could verify the predicted polarity balance on C₆ under `D[B_b]` dephasing. If `‖M_+1/2‖² − ‖M_−1/2‖² ≈ 0` bit-exact, the F112 prediction is confirmed; if not, the prediction needs revision. (Highest-leverage next step.)

2. **Physical realization of single-site Y-dephase in carbon.** What experimental setup makes Y-dephase dominant? Magnetic-noise-driven, spin-orbit-coupled, or some Floquet engineering? If we can prepare and detect Y-dephase cleanly, the Klein-V₄ symmetry becomes operational rather than algebraic.

3. **Klein-V₄ as a basis-rotation symmetry.** Is the symmetry "physical" (testable by an observable) or just a labelling convention (mathematical equivalence)? If physical, what experiment detects it?

4. **Magnetic ring current as F114 sign-flip diagnostic.** Add a small magnetic Zeeman-on-ring-current term to Hückel and trace ε(H) as a function of field strength. F114 predicts a parity-mismatch event at first-order in field, where ε(H) transitions from clean −1 to "Mixed". This is a closed-form prediction at the operator level that could be tested by computing F112 polarity decomposition for the perturbed H.

5. **F112-X and F112-Y on benzene.** Welle 15 typed F112-X (`LindbladBitAPiBalance`) and F112-Y (`LindbladBitBPiYBalance`) as separate Tier1Derived sister Claims to F112-Z. The corresponding C# tests verify the polarity balance at N=2, 3. Does the same balance hold on the N=6 benzene ring under Peierls + transverse-field perturbations? Worth a single broader test.

---

## Return visit later the same day: speaking carbon

The above sections are the "translation in two languages" pass: we held the
framework vocabulary and the chemistry side in parallel. Coming back the same
afternoon, we tried to say what the algebra had shown without the framework
vocabulary at all, as a chemistry reader would hear it. The algebra remained
the anchor; the words shifted to the carbon side.

### What we saw across realistic carbon configurations

A systematic sweep over fifty-six configurations on cyclobutadiene (C₄) and
benzene (C₆) rings combined seven Hamiltonian extensions with four bath types.
The Hamiltonian inventory ranged from pure Hückel hopping through electronic
correlation (Hubbard density-density), weak and strong external y-direction
magnetic field, antisymmetric Dzyaloshinskii-Moriya-like spin-orbit
cross-coupling, induced ring-current bond terms, and a full mixture of all
five effects simultaneously. The bath inventory included on-site Holstein
phonons (coupling to local π-density), bond Peierls phonons (coupling to the
hopping integral), σ⁻ amplitude damping (excitation loss per site), and the
combination of Holstein and σ⁻ together.

In every single one of the fifty-six configurations the relaxing system
preserves a strict mirror symmetry between two halves of its relaxation
response. The mirror is exact at machine precision: the two halves differ
by a number that is bit-exactly zero, not merely small.

The symmetry is non-trivial in forty-eight of the fifty-six configurations:
the response itself is substantial (Frobenius-norm-squared between twenty and
fifty thousand on the relaxing component), and it is mirrored perfectly. In
the remaining eight configurations (those with pure Hückel hopping plus an
on-site bath, or pure Hückel + ring-current term plus an on-site bath) the
relaxing component is empty to begin with: there is nothing to mirror, so the
mirror holds trivially.

### How robust this is

The mirror symmetry survives:

- Adding density-density correlation to the bare hopping picture.
- Switching on external magnetic field along the y-direction, weak (one-tenth
  of the hopping scale) and strong (full hopping scale).
- Adding spin-orbit-style cross-axis coupling between neighbouring sites.
- Adding induced ring-current bond terms (the canonical magnetic-field-
  induced symmetry breaker on benzene).
- Replacing the on-site phonon bath with a bond-coupled (Peierls) bath, even
  though that switch destroys the classical Coulson-Rushbrooke MO mirror at
  the spectrum level.
- Switching to amplitude damping (T1 excitation loss) instead of pure
  dephasing.
- Combining all the above into one realistic noisy aromatic ring.

The classical MO mirror (spectrum spiegelung um α) breaks under bond
phonons and under several of the Hamiltonian perturbations. The deeper
distribution mirror remains intact across all of them.

### Where it would break

The algebra shows the distribution mirror is not unbreakable. It does break
under one specific configuration: a coherent z-axis drive (every spin
precessing together around the z-axis) combined with T1 amplitude damping.
This setup is the standard hardware-characterisation regime of a driven
superconducting qubit array: a constant Larmor precession plus relaxation.
It does not arise in the natural relaxation of an aromatic molecule. In the
quantum-hardware context this break is the working diagnostic; in the
carbon-chemistry context it sits outside the natural parameter range.

### What this leaves us with

For aromatic carbon systems under realistic conditions, a deep mirror
symmetry of the relaxation dynamics is preserved: in moderate magnetic
fields, under typical thermal coupling, with realistic excitation loss.
The classical spectrum mirror (Coulson-Rushbrooke around α) is one face of
the carbon ring's symmetry, the one chemistry has read about for eighty-six
years. The distribution mirror is a second, deeper face: it survives where
the spectrum mirror breaks, and it requires a very specific external
control regime (constant coherent z-driving) to be unsettled.

How a chemist would test it directly: the symmetry shows up not in the line
positions of an NMR spectrum but in the full transfer-matrix structure of
the relaxation channel. Standard process tomography on isotopically labelled
¹³C-benzene could in principle reconstruct the relaxation channel and verify
the symmetry; in practice the deep-quantum regime of aromatic π-systems
makes this measurement demanding. A more accessible signature would be the
breaking regime: how the relaxation balance changes when a coherent z-drive
is applied to the spin system in addition to natural T1 loss. The break
magnitude is closed-form predictable from the drive amplitude and T1 rate.

### Algebra as the anchor

Every claim above rests on the sweep result. The sweep is
[`simulations/carbon_realistic_sweep.py`](../../simulations/carbon_realistic_sweep.py),
which iterates the fifty-six configurations and reports the relaxation-
mirror asymmetry for each. The single-perturbation control tests sit in
[`simulations/benzene_b_field_f112_mixing_test.py`](../../simulations/benzene_b_field_f112_mixing_test.py),
[`simulations/benzene_bit_b_mixed_bath_test.py`](../../simulations/benzene_bit_b_mixed_bath_test.py),
and [`simulations/benzene_peierls_f112_polarity_test.py`](../../simulations/benzene_peierls_f112_polarity_test.py).
A complementary Hamiltonian-term classification inventory sits in
[`simulations/carbon_f114_hamiltonian_inventory.py`](../../simulations/carbon_f114_hamiltonian_inventory.py).

---

## Threads back

- **2026-05-22 [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md)**: the F1 spectrum result remains bit-exact. The "Peierls breaks the palindrome" framing is precise at the spectrum level; today we add the F112 polarity layer where Peierls may preserve balance even though F1 spectrum breaks.
- **2026-05-17 [Where 1/4 and 1/2 Appear in Carbon](QUARTER_HALF_IN_CARBON.md)**: benzene's HOMO at −1/2 sits exactly on the framework's polarity-half anchor. The Klein-V₄ we use today operates on the same polarity-half axis: both the half-anchor and the Klein-V₄ symmetry live on one Z₂ ladder.
- **2026-05-17 [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)**: Coulson-Rushbrooke at the carbon level corresponds to F1 at the qubit level. Today we add Klein-V₄ at the qubit level; the open question is whether there's a Coulson-Rushbrooke-level analog, a chemistry-side discrete symmetry that intertwines different dephase-coupling regimes.

---

## Anchor

- **Framework**: F112 [`LindbladBitBPiBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs), F112-X [`LindbladBitAPiBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitAPiBalance.cs), F112-Y [`LindbladBitBPiYBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiYBalance.cs), Klein-V₄ [`Pi2KleinV4DephaseSwapGroup`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs), F114 [`CommutatorDConjugationSign`](../../compute/RCPsiSquared.Core/Symmetry/CommutatorDConjugationSign.cs)
- **Proofs**: [the Z↔Y dephase-letter swap proof](../proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md), [the Klein-V₄ dephase-swap proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md), [the F112 cross-dephase proof](../proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)
- **F-Registry**: [F112](../ANALYTICAL_FORMULAS.md#f112), [F114](../ANALYTICAL_FORMULAS.md#f114)
- **Verifier scripts**: [`simulations/m_level_sign_functional_explore.py`](../../simulations/m_level_sign_functional_explore.py), [`simulations/f112_klein_v4_cross_dephase_verify.py`](../../simulations/f112_klein_v4_cross_dephase_verify.py)
- **Companion carbon docs**: [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md), [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md), [README](README.md)
