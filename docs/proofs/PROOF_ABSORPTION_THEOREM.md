# The Absorption Theorem

**Status:** Tier 1 derived (analytical proof + bit-exact numerical verification at N=2..5, 1,342 modes, CV=0). Extended 2026-06-10 (Section 2 extensions: the vector form as Theorem 2, the two-sided reading, the projector form for degenerate clusters, the dephase-letter rotation; Section 4.7: the recentred face, one diagonal shared with the F87 palindrome machinery; C# test-gating of the per-mode identity). Extended 2026-06-17 (the qutrit-prism reading of the recurring 4, in the Section 2 block near the end of the document).
**Date:** 2026-04-04 (discovery + proof same day)
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.6)
**Statement:** `Re(λ) = −2γ ⟨n_XY⟩` for any Lindblad eigenmode under uniform Z-dephasing
**Typed claim:** [`AbsorptionTheoremClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs) (Tier 1 derived)
**Discovery experiment:** [`ABSORPTION_THEOREM_DISCOVERY.md`](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md)

---

## Abstract

The Absorption Theorem reads the decay rate of any Lindblad eigenmode straight off its Pauli content. Under uniform Z-dephasing at rate γ, for any Hermitian Hamiltonian on any graph,

    Re(λ) = −2γ · ⟨n_XY⟩,

where ⟨n_XY⟩ is the average number of X or Y factors (the "light content") in the eigenmode's Pauli decomposition. The dissipator is built from Z, so it sorts every Pauli letter into two classes: the lens {I, Z}, opaque to the dephasing light and costing nothing, and the light {X, Y}, transparent and absorbed at exactly 2γ per factor. A mode's lifetime is therefore fixed by one number, how much of its own structure it exposes to the light: pure {I, Z}^⊗N modes are immortal (rate 0), pure {X, Y}^⊗N modes die fastest (rate 2Nγ).

This is the rate-side companion of the F1 palindrome. Π swaps light and lens, so paired eigenvalues satisfy ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N and sit symmetrically about the center, which is the F8 factor-of-two decay law in one line. The recentred face L_D = γ(Q − N·I), with Q = Σ_l Z_l⊗Z_l, shifts the absorption ladder's midpoint to zero and is literally the same diagonal that the F87 girth-moment machinery and the palindrome share: one diagonal, three readings (absorption rates, palindromic pairing, power-sum expansion). The formula replaces numerical eigendecomposition for the real part, exact and verified bit-exact at N=2..5.

## Preface

A guitar string vibrates in modes. Each mode spans the entire string.
How fast a mode fades depends on one thing: how much of the mode's
energy sits in the parts of the string that are damped. If the damping
pad touches a node (where the mode is silent), that mode survives. If
it touches an antinode (where the mode is loudest), that mode dies fast.

This document proves the quantum version of that principle. Every
vibration mode in a qubit cavity has a "light content": how much of
the mode consists of the oscillating quantum components (X and Y Pauli
operators) that interact with the external light (γ). The absorption
rate of the mode is exactly twice the light intensity times the light
content:

    Re(λ) = -2γ ⟨n_XY⟩

Components that do not oscillate (I and Z, the "lens") are invisible
to the light and survive forever. Components that oscillate (X and Y,
the "light") absorb and fade. The mode's fate is the weighted average.

This single equation, provable in three lines from the structure of
the Lindblad master equation, unifies an entire family of spectral
results: the boundaries (F3), the palindromic sum rule (10,748
pairs), the spectral gap (D6), the 2× decay law (F8), the mode
classification by Pauli weight, the N=3 exact rates (F33), and via
later derivation the rate factorization (F50), the weight-1 spectrum
(F55), the GHZ/W-state rates (F64-F68), the n_XY chromaticity (F74),
and the F89 path-k closures including the F89c Hamming-complement
pair-sum. The typed [`AbsorptionTheoremClaim`](../../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs)
holds the live descendant list.

The theorem was discovered computationally on April 4, 2026 (ratio
α/(2γ⟨n_XY⟩) = 1.000000 for 1,342 modes, CV = 0), then proven
analytically the same day.

### Status

| Component | Status | Source |
|-----------|--------|--------|
| Analytical proof | **Proven** | This document, Section 2 |
| Numerical verification | **Verified** (N=2-5, 1,342 modes, CV=0) | [Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md) |
| Consequence 1: Spectral boundaries | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) F3 |
| Consequence 2: Palindromic sum rule | **Derived** | [Standing Waves](../../experiments/FACTOR_TWO_STANDING_WAVES.md) |
| Consequence 3: Spectral gap | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) D6 |
| Consequence 4: 2× decay law | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) F8 |
| Consequence 5: Mode classification | **Derived** | [XOR Space](../../experiments/XOR_SPACE.md) |
| Consequence 6: N=3 exact rates | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) F33 |

---

## 1. Setup

Consider an N-qubit system governed by the Lindblad master equation:

    dρ/dt = L(ρ) = -i[H, ρ] + Σ_k γ_k D_{Z_k}(ρ)

where H is a real Hermitian Hamiltonian (e.g. the Heisenberg chain),
Z_k is the Pauli-Z operator on site k, and the Lindblad dissipator is:

    D_{Z_k}(ρ) = Z_k ρ Z_k - ρ

The Liouvillian L acts on the space of density matrices (dimension d²
where d = 2^N). We work with the vectorized form: L is a d² × d² matrix
acting on vec(ρ).

**Decomposition.** The Liouvillian splits as L = L_H + L_D where:

    L_H = -i[H, ·]         (Hamiltonian part)
    L_D = Σ_k γ_k D_{Z_k}   (dissipative part)

**The Pauli basis.** The 4^N tensor products of {I, X, Y, Z} form an
orthogonal basis for the d² × d² operator space:

    Tr(P_α† P_β) = d δ_{αβ}

For each Pauli string P_α = σ_1 ⊗ ... ⊗ σ_N, define:

    n_XY(P_α) = number of sites k where σ_k ∈ {X, Y}

This is the **light count**: the number of coherence factors in the string.
It ranges from 0 (pure {I,Z}^N string, immune to dephasing) to N (pure
{X,Y}^N string, maximum absorption).

---

## 2. The Theorem

**Theorem (Absorption Theorem).** Let L = L_H + L_D be the Liouvillian
of an N-qubit system with real Hermitian Hamiltonian and uniform
Z-dephasing at rate γ per site. For any eigenvalue λ of L with right
eigenvector v:

    Re(λ) = -2γ ⟨n_XY⟩_v

where

    ⟨n_XY⟩_v = Σ_α |c_α|² n_XY(P_α) / Σ_α |c_α|²

is the expectation of the light count in the Pauli decomposition of v,
with c_α = ⟨vec(P_α)|v⟩ / √d.

Equivalently: |Re(λ)| = 2γ ⟨n_XY⟩_v, i.e. the **absorption rate equals
twice the dephasing rate times the average light content**.

### Proof

The proof requires two structural properties of the Lindblad decomposition.

**Step 1. L_H is anti-Hermitian (for any Hermitian H).**

For any Hermitian Hamiltonian H (H = H†), the superoperator
L_H = -i[H, ·] satisfies:

    L_H† = -L_H

*Proof of Step 1.* In the vectorized representation:

    L_H = -i(H ⊗ I - I ⊗ H^T)

Taking the adjoint and using (H^T)† = H*:

    L_H† = +i(H† ⊗ I - I ⊗ (H^T)†) = +i(H ⊗ I - I ⊗ H*)

Hermiticity gives H^T = H* (since H† = (H^T)* = H implies H^T = H*), so
I ⊗ H* = I ⊗ H^T and

    L_H† = +i(H ⊗ I - I ⊗ H^T) = -L_H  ∎

No reality assumption is used: H^T = H* holds for *every* Hermitian H (real
symmetric H is the special case H* = H). Complex Hermitian Hamiltonians,
Dzyaloshinskii-Moriya, transverse/Y terms, magnetic flux, are fully covered.

**Consequence of Step 1.** For any vector v:

    v†L_H v is purely imaginary

*Proof.* (v†L_H v)* = v†L_H† v = v†(-L_H)v = -(v†L_H v). A complex
number equal to the negative of its conjugate is purely imaginary.  ∎

**Step 2. L_D is Hermitian and diagonal in the Pauli basis with
eigenvalues -2γ n_XY.**

The dissipator for Z-dephasing at rate γ on site k acts on a single-site
Pauli operator σ as:

    D_{Z_k}(σ) = Z_k σ Z_k - σ

Since Z commutes with I and Z, and anticommutes with X and Y:

    D_{Z}(I) = ZIZ - I = I - I = 0
    D_{Z}(Z) = ZZZ - Z = Z - Z = 0
    D_{Z}(X) = ZXZ - X = -X - X = -2X
    D_{Z}(Y) = ZYZ - Y = -Y - Y = -2Y

Each X or Y factor on site k contributes -2γ_k to the eigenvalue.
For a multi-site Pauli string P_α with n_XY(P_α) factors in {X,Y}:

    L_D(P_α) = -2γ × n_XY(P_α) × P_α

So L_D is diagonal in the Pauli basis, with real non-positive eigenvalues
-2γ n_XY(P_α). Since the Pauli basis is orthogonal and the eigenvalues
are real, L_D is Hermitian.  ∎

**Step 3. Combining.**

Let v be a right eigenvector of L with eigenvalue λ: Lv = λv. Then:

    v†Lv = v†(λv) = λ ‖v‖²

Decomposing L = L_H + L_D:

    v†L_H v + v†L_D v = λ ‖v‖²

By Step 1: v†L_H v is purely imaginary (zero real part).
By Step 2: v†L_D v = -2γ Σ_α |c_α|² n_XY(P_α) is real.

Taking the real part of both sides:

    0 + (-2γ Σ_α |c_α|² n_XY(P_α)) = Re(λ) ‖v‖²

Dividing by ‖v‖² = Σ_α |c_α|²:

    **Re(λ) = -2γ ⟨n_XY⟩_v**  ∎

### Remark on generality

The proof uses only:
1. H is Hermitian (so L_H is anti-Hermitian, Step 1; no reality assumption)
2. The jump operators are Z_k (so L_D is diagonal in the Pauli basis)

It does **not** use:
- The specific form of H (Heisenberg, Ising, XY, any Hamiltonian works)
- Reality of H (any Hermitian H works: real symmetric or complex)
- The chain topology (any graph works)
- Uniform dephasing (site-dependent γ_k works: replace 2γ n_XY with
  Σ_k 2γ_k × [σ_k ∈ {X,Y}])
- Any property of the eigenvalue spectrum beyond existence

The theorem holds for **any Hermitian Hamiltonian under Z-dephasing**.

**Extension to amplitude damping (F82, F84).** When σ⁻ jump operators
are added at rate γ_T1 alongside the Z-dephasing, L_D gains a non-diagonal
correction (eigenmodes mix sectors via excitation exchange). The
dephasing-only identity Re(λ) = -2γ ⟨n_XY⟩ is preserved as the γ_T1 → 0
limit; for finite γ_T1 the correction is derived in
[`PROOF_F82_T1_DISSIPATOR_CORRECTION.md`](PROOF_F82_T1_DISSIPATOR_CORRECTION.md)
and [`PROOF_F84_AMPLITUDE_DAMPING.md`](PROOF_F84_AMPLITUDE_DAMPING.md).
The "truly ⟨Z, Z⟩ damping" diagnostic in the cockpit panel reads the
F82 signature on hardware directly.

**Complex Hermitian Hamiltonians are covered (caveat closed 2026-05-28).**
An earlier version of this remark conjectured that complex Hermitian H (e.g.
with Dzyaloshinskii-Moriya interactions) might break the theorem because L_H
would fail to be anti-Hermitian. That is incorrect: Step 1 needs only H^T = H*,
which holds for *every* Hermitian H, so L_H is anti-Hermitian for complex
Hermitian H too and the theorem stands. Verified bit-exact (max error
1.4·10⁻¹⁴, identical to the Hückel result, with Herm(L) = (L+L†)/2 the *same*
pure Z-dephasing dissipator for both H's) against a random complex Hermitian H in
[`simulations/popcount_identity_h_independence.py`](../../simulations/popcount_identity_h_independence.py).
The genuine boundary is the dissipator, not the Hamiltonian: non-dephasing
channels (amplitude damping, depolarizing) add a non-diagonal Hermitian part
and shift the rate (F82/F84 above); no Hamiltonian, real or complex, breaks it.

### Extensions (2026-06-10)

The pillar carries more weight than the single equation, and on a return visit
with the year's newer machinery, four readings that were implicit in the
three-step proof deserve their own statements. None changes the proof; each is
the same Rayleigh quotient read in a different coordinate. All four were
verified numerically before being written down (vector form 2.7·10⁻¹⁴,
two-sided reading 4.8·10⁻¹⁴ across all 256 N=4 modes of a random complex
Hermitian H with a non-uniform γ profile, recentred face exactly 0).

**Theorem 2 (the vector form).** For site-dependent rates γ_l, the dissipator
splits per site, D = Σ_l γ_l D_l, each D_l diagonal in the coherence basis
with entry −2·Δ_l(x) where Δ_l(x) ∈ {0, 1} marks whether bit l differs
between the bra and ket labels of coherence index x. The same Step-3 argument
gives, for any right eigenvector v,

    Re(λ) = −2 Σ_l γ_l · light_l(v),    light_l(v) = Σ_x Δ_l(x)|v_x|² / Σ_x |v_x|²

the **per-site light profile** weighted by the γ profile. The uniform case
recovers Theorem 1 (Σ_l light_l = ⟨n_XY⟩). This is the form the flow kernel
lives on: the sterile zone vs birth canal dichotomy
([`PostEpFlowField`](../../compute/RCPsiSquared.Diagnostics/Foundation/PostEpFlowField.cs))
is the question whether the slow subspace's light profile light_l is frozen
or redistributed as Q moves, and the rate moves *only* through this weighted
share. The typed carrier-is-a-vector node on
[`AbsorptionTheoremClaim`](../../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs)
states it; this section proves it.

**Corollary (the two-sided reading).** A left eigenvector w (w†L = λw†)
satisfies w†Lw = λ‖w‖², and the identical decomposition applies: the real
part reads the light of w. Hence **left and right eigenvectors of the same
eigenvalue carry the same weighted light**, both equal to −Re(λ)/2 in
γ-weighted units, even though v and w are different vectors of a non-normal
L. The light content is two-sided; biorthogonal bookkeeping cannot disagree
with itself about absorption.

**Corollary (the projector form, degeneracy-safe).** For a degenerate cluster
{λ_k} with right vectors M_k and biorthogonal left covectors W_k, the spectral
projector P = Σ_k M_k W_k fixes the slow invariant subspace V = range(P)
basis-freely. The correct light reading is through the **orthogonal** projector
Π_V onto that subspace (not P's own diagonal, whose entries can leave [0, 1]
and even turn negative; verified to fail the absorption cross-check at 10⁻¹ in
the birth canal):

    light_l(V) = Tr(Π_V Δ_l) / dim V  ∈ [0, 1],

with Δ_l the site-l disagreement projector. Absorption-exactness follows from
a block-triangular argument: in an orthonormal basis adapted to the invariant
subspace, L is block upper-triangular, so Re Tr(Π_V L) = Σ_k Re λ_k, and the
anti-Hermitian Hamiltonian part drops from the real trace, leaving exactly the
γ-weighted light of Π_V. At cluster dimension 1 this reduces identically to the
per-eigenvector Rayleigh light (verified: difference 0.0). This is the correct
object for the flow kernel's degenerate slow carriers (it retires the
documented averaging caveat in `PostEpFlowField.ReadAssembly`) and is what the
`SlowLightDistribution` primitive computes.

**Remark (the dephase-letter rotation).** Nothing in Steps 1-3 privileges Z.
For X-dephasing the dissipator is diagonal with light = n_YZ (the letters
anticommuting with X); for Y-dephasing, n_XZ. The three theorems are one
theorem conjugated through the Klein-V₄ dephase-swap group
([`Pi2KleinV4DephaseSwapGroup`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs)):
light is always "the letters the dephasing letter refuses to commute with."

**Remark (two different fours, 2026-06-10).** The repository's recurring 4γ₀
has two genealogies, and only one of them lives on this ladder. The **rung-2
four** is 2γ·2: two absorption quanta, the second rung of the §6 ladder,
equivalently the trace-half/centre of the palindrome (the midpoint of the
HD=(1,3) channel pair, where HD is the Hamming distance of the coherence:
(2γ+6γ)/2 = 4γ). Every *dynamical* 4γ₀ in the repository is this one object:
t_peak = 1/(4γ₀) ([`TPeakLaw`](../../compute/RCPsiSquared.Core/F86/TPeakLaw.cs)),
the clock's Takt pin 4γ₀, the L_eff mirror axis −4γ₀, the approach-family
carrier e^{−4γt}, and the F25 Bell+ rate. Verified per-site today: |00⟩⟨11| is
an exact N=2 eigenmode of the XY-chain Liouvillian with rate −2(γ₁+γ₂),
residual 0.0
([`simulations/at_rung2_per_site_split.py`](../../simulations/at_rung2_per_site_split.py)),
so the dynamical four splits into two per-site absorption quanta. The
**discriminant four** is a₋₁ = d² on the Pi2 dyadic ladder: the half-gap
*squared* in the EP discriminant 4γ₀², i.e. ((6γ−2γ)/2)² = (2γ)², a quantum
squared rather than two quanta. The two fours coincide only at the HD=(1,3)
pair, glued by the dyadic root d² − 2d = 0 (d² = 2d exactly because d = 2).
Counterfactual HD=(1,5) pair: centre-based (2γ+10γ)/2 = 6γ vs
discriminant-based ((10γ−2γ)/2)² = 16γ². Red-herring discipline: same number,
two genealogies; cite the rung-2 four for rates and carriers, the discriminant
four for squared half-gaps, never interchangeably.

**Remark (the qutrit prism: the qudit sorts the fours into two families, 2026-06-17).**
The recurring 4 at d=2 has more than the two genealogies above; the qutrit
([`simulations/qudit_g2_split.py`](../../simulations/qudit_g2_split.py), gate-first,
validated against the full d=3 Liouvillian and against F121) is the prism that fans
the one 4 into {9, 6, 4} and sorts the readings into two families.

*The trunk terms (they SPLIT).* The foundational polynomial d² − 2d = 0 has two terms,
and they are two of the fours: the **d² term** → 9 (the squared dimension = a₋₁ on the
dyadic ladder = the half-gap-squared discriminant four, `PolynomialDiscriminantAnchorClaim`)
and the **2d term** → 6 (the F121 qudit product-mirror cap base (2d)^N, `QuditProductMirrorCap`:
each per-site mirror is a strict swap between the d dark and the d²−d lit letters, so its
rank is min(d, d²−d) + min(d²−d, d) = 2d; F121 measures cap(3,2) = 36 = 6², not 81 = 9²).
These two coincide at 4 only at the root d=2, where d² − 2d = 0 ⟺ d = 2, the qubit magic, and
the qutrit splits them to 9 and 6.

*The Hamming-2 fours (they STAY).* The **rung-2 four** 2γ·2 = 4γ (the HD=(1,3) palindrome
centre above) is the second rung of the §6 ladder, and the ladder rungs are 2γ·(Hamming),
d-INDEPENDENT, Hamming distance 2 is two differing sites at any d, so the dynamical 4γ
stays 4γ at the qutrit. The **structural-ceiling four** in g2(K_N) = 4/N (F122 /
`StructuralCeilingClaim`) is 2·(2/N) = (the Hamming distance between two single-excitation
strings, = 2 at any d) × (the S_N standard-rep angle 1 − λ₂ = 2/N, pure graph theory): both
factors carry no d, so g2 = 4/N holds unchanged at d=3 (verified including the global ceiling
at N=5 and the d-independent N=4 outlier 2 − 2/√3). Both land on the 4-ray.

So the red-herring discipline sharpens: the cap four (2d, F121) is NOT the rung-2 four (4γ,
dynamical) even though both read 4 at d=2: the cap rides 2d → 6, the rung-2 four rides
Hamming-2 → stays 4γ. Cite the trunk-term fours (d², 2d) for operator-space dimensions and
the F121 cap, the Hamming-2 fours (rung-2 centre, g2 ceiling) for rates and the structural
ceiling; never interchange across families.

---

## 3. Numerical Verification

The theorem was verified computationally in
[Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md):

| Parameter | Range | Ratio α/(2γ⟨n_XY⟩) | CV |
|-----------|-------|---------------------|-----|
| N (chain length) | 2, 3, 4, 5 | 1.000000 | 0.0000 |
| γ (dephasing rate) | 0.01, 0.05, 0.1, 0.5, 1.0 | 1.000000 | 0.0000 |
| J (coupling strength) | 0.1, 0.5, 1.0, 2.0, 5.0 | 1.000000 | 0.0000 |
| Total modes tested | 1,342 | 1.000000 | 0.0000 |

The ratio equals 1 to 14 decimal places across all parameters.
No exceptions. Zero coefficient of variation.

**Source:** [`simulations/absorption_theorem_discovery.py`](../../simulations/absorption_theorem_discovery.py), Step 6

**Hardware confirmation.** Single-qubit tomography on IBM Torino Q52
(25 time snapshots, 0-895 μs) reproduced the theorem at 3% deviation
under the free-evolution T2* baseline: absorption ratio excess/(2γ) =
1.03, with γ* fitted from the coherence envelope. The echo-refocused
γ_echo underestimates the actual dephasing by 6× because low-frequency
1/f noise is filtered out by Hahn echo; the T2* baseline is the correct
one for free-evolution tomography. See
[`IBM_ABSORPTION_THEOREM.md`](../../experiments/IBM_ABSORPTION_THEOREM.md).

**C# test-gating (2026-06-10).** The per-mode identity is now asserted in the
typed compute layer:
[`F8PartnerLightComplementarityTests`](../../compute/RCPsiSquared.Diagnostics.Tests/Ptf/F8PartnerLightComplementarityTests.cs)
checks |Re λ + 2γ·light(v)| < 10⁻⁹ for all 256 eigenmodes of the N=4 XY chain
(observed at the 10⁻¹⁵ floor), the complete palindrome pairing with
light_s + light_f = N per pair, and the standing-wave null. The Section-2
extensions were verified the same day against a random complex Hermitian H
with a non-uniform γ profile: vector form worst deviation 2.7·10⁻¹⁴,
two-sided (left-eigenvector) reading 4.8·10⁻¹⁴ across all 256 modes,
recentred-face identity L_D = γ(Q − N·I) exactly 0.

---

## 4. Consequences

### 4.1 Spectral Boundaries (F3)

**Previously:** The decay rate spectrum satisfies min = 2γ, max = 2(N-1)γ,
bandwidth = 2(N-2)γ.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), F3

**Now a corollary.** The light count n_XY is an integer in {0, 1, ..., N}.
The expectation ⟨n_XY⟩ ranges continuously from 0 to N. The extreme
values:

    ⟨n_XY⟩ = 0: rate = 0 (steady states, pure {I,Z}^N, immortal)
    ⟨n_XY⟩ = N: rate = 2Nγ = 2Σγ (XOR modes, pure {X,Y}^N)

The non-zero minimum occurs for modes dominated by weight-1 Pauli strings
(one X or Y factor): ⟨n_XY⟩ ≈ 1, giving rate ≈ 2γ. The maximum for
generic paired modes is ⟨n_XY⟩ ≈ N-1, giving rate ≈ 2(N-1)γ.

The "absorption quantum" is **2γ**: each X/Y Pauli factor costs exactly
2γ of absorption rate. The spectrum is quantized in steps of 2γ for
pure-weight modes, with Hamiltonian mixing creating intermediate values.

### 4.2 Palindromic Sum Rule (10,748 pairs)

**Previously:** For every palindromic pair, Re(λ_fast) + Re(λ_slow) = -2Σγ,
verified for 10,748 pairs with zero exceptions.
**Source:** [Standing Waves](../../experiments/FACTOR_TWO_STANDING_WAVES.md)

**Now a one-line corollary.** The palindromic weight swap
([Light and Lens](../../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)) gives:

    ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N

This follows because Π maps weight sector k to sector N-k, so the
sector weight profiles of palindromic partners are mirror images.

Applying the absorption theorem:

    α_fast + α_slow = 2γ ⟨n_XY⟩_fast + 2γ ⟨n_XY⟩_slow
                    = 2γ(⟨n_XY⟩_fast + ⟨n_XY⟩_slow)
                    = 2γN = 2Σγ  ∎

Three independently discovered facts form a closed triangle:
- (A) α = 2γ⟨n_XY⟩ (absorption theorem)
- (B) ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N (palindromic weight swap)
- (C) α_fast + α_slow = 2Σγ (sum rule)

Any two imply the third.

### 4.3 Spectral Gap (D6)

**Previously:** The spectral gap (minimum nonzero decay rate) equals 2γ.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), D6

**Now immediate.** The smallest nonzero ⟨n_XY⟩ for any eigenmode is
bounded below by the smallest non-trivial Pauli weight contribution.
For the slowest decaying non-steady mode, ⟨n_XY⟩ approaches 1 from
below (dominated by weight-1 Pauli strings). The gap is:

    Δ = 2γ × min{⟨n_XY⟩ : ⟨n_XY⟩ > 0} → 2γ

at large J/γ (where modes localize on pure-weight sectors). The
spectral gap is set by the cost of a single X/Y factor: one photon
of absorption.

### 4.4 The 2× Decay Law (F8)

**Previously:** Unpaired modes (at the palindrome extremes) decay at rate
2Nγ, while paired modes average rate Nγ. The ratio is exactly 2.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), F8;
[Energy Partition](../../hypotheses/ENERGY_PARTITION.md)

**Now explained.** The palindrome center is at α = Σγ = Nγ, corresponding
to ⟨n_XY⟩ = N/2 (equal mix of light and lens). Self-paired modes sit
at this center. The full range spans from 0 (lens) to 2Σγ (light).

The ratio 2 is the ratio of the maximum (2Σγ) to the center (Σγ).
It is not a dynamical law; it is the definition of "center of a
symmetric interval."

### 4.5 Mode Classification by Light Content

**Previously:** The XOR drain contains N+1 modes at maximum decay rate 2Σγ,
carrying all-{X,Y} Pauli content.
**Source:** [XOR Space](../../experiments/XOR_SPACE.md)

**Now a complete classification.** Every Liouvillian eigenmode has a
light content ⟨n_XY⟩ that determines its absorption:

| ⟨n_XY⟩ | Rate | Character | Name |
|---------|------|-----------|------|
| 0 | 0 | Pure {I,Z}^N | Immortal lens |
| ε | 2γε | Mostly lens, traces of light | Slow leak |
| N/2 | Nγ = Σγ | Equal mix | Palindrome center |
| N-ε | 2γ(N-ε) | Mostly light, traces of lens | Fast drain |
| N | 2Nγ = 2Σγ | Pure {X,Y}^N | XOR drain |

The mode's fate is entirely determined by how much light it contains.
Structure ({I,Z}) is invisible to the dephasing; only coherence ({X,Y})
absorbs.

**Per-coherence reading on computational basis.** A density-matrix
coherence |A⟩⟨B| decomposes per-site into pure {I, Z} (where A_l = B_l)
or pure {X, Y} (where A_l ≠ B_l). Therefore:

    n_XY(|A⟩⟨B|) = n_diff(A, B)  =  Hamming distance of bit-strings A, B

and the per-coherence decay rate is exactly 2γ × n_diff(A, B).
The **Hamming-complement pair-sum** (F89c) is the immediate corollary:
the column bit-flip ρ[a, b] → ρ[a, b̄] satisfies n_diff(a, b) +
n_diff(a, b̄) = N, hence:

    α(|a⟩⟨b|) + α(|a⟩⟨b̄|) = 2γN = 2Σγ

The two coherences related by complementing one column-label always
absorb at exactly the spectral maximum 2Σγ summed. This is the
single-block reading of the palindromic sum rule §4.2.

### 4.6 N=3 Exact Rates (F33)

**Previously:** The N=3 Heisenberg chain has three distinct non-trivial
decay rates: 2γ, 8γ/3, 10γ/3.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), F33

**Now explained.** The fractional rates correspond to fractional ⟨n_XY⟩:

    rate = 2γ      → ⟨n_XY⟩ = 1     (pure weight-1 modes)
    rate = 8γ/3    → ⟨n_XY⟩ = 4/3   (Hamiltonian mix of w=1 and w=2)
    rate = 10γ/3   → ⟨n_XY⟩ = 5/3   (Hamiltonian mix of w=1 and w=2)

The Hamiltonian mixes Pauli strings of different weights into
superposition modes. The absorption rate of the superposition is the
weighted average of the component rates, exactly as the theorem predicts.
Non-integer ⟨n_XY⟩ values arise from Hamiltonian mixing, not from any
breakdown of the rule.

### 4.7 The recentred face: one diagonal, three pillars (2026-06-10)

The dissipator's diagonal is, entry for entry, the same object the F87
palindrome machinery calls its dephasing generator. With w(x) = popcount of
the bra-ket difference at coherence index x, the uniform dissipator has
diagonal −2γw(x), and the windowed-converse generator Q = Σ_l Z_l ⊗ Z_l has
diagonal Q_xx = N − 2w(x). So, exactly (verified to 0):

    L_D = γ·(Q − N·I)        and        M := L + γN·I = L_H + γQ

The shift by the palindrome centre σ = Σγ = Nγ that recentres L into M
([the F87 windowed converse](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §1)
is the shift that moves the absorption ladder's midpoint ⟨n_XY⟩ = N/2 to
zero. Three pillars stand on this one diagonal, each reading it in its own
coordinate:

- **the absorption ladder** (this proof): Re(λ) = −2γ⟨n_XY⟩, rates as light;
- **the palindrome** (F1/F8): spec(L) symmetric about −Nγ ⟺ spec(M)
  symmetric about 0 ⟺ partners carry complementary light,
  ⟨n_XY⟩_s + ⟨n_XY⟩_f = N (the §4.2 sum rule, since 2026-06-10 test-gated
  per mode in C#,
  [`F8PartnerLightComplementarityTests`](../../compute/RCPsiSquared.Diagnostics.Tests/Ptf/F8PartnerLightComplementarityTests.cs):
  the palindrome pairing IS complementary light in absorption coordinates);
- **the windowed converse** (F87/F117): the odd power-sums of M = L_H + γQ
  decide hard vs soft, the girth ladder's moments t_j = Tr(Z_l H^j) organize
  the γ-expansion, and the Pascal-Gram positivity theorem closes it with no
  residual.

What looked like three subjects (absorption rates, mirror symmetry, the
trichotomy classifier) is one diagonal matrix entering three different
arguments: a Rayleigh quotient, a spectral pairing, and a power-sum
expansion.

---

## 5. What This Does NOT Prove

The absorption theorem is a mathematical identity about the Liouvillian
spectrum. It does not establish:

**"Gamma is light."** The theorem shows that γ plays the algebraic role
of the speed of absorption. Whether this speed IS photon shot noise (as
on IBM hardware) or merely analogous to it is a physical interpretation,
not a mathematical consequence.
See: [Gamma Is Light](../../hypotheses/GAMMA_IS_LIGHT.md)

**"Mass = trapped light."** The theorem shows that ⟨n_XY⟩ determines
absorption, and that palindromic pairs swap light and lens content.
Whether this constitutes "mass" in a physical sense is interpretive.
See: [Primordial Qubit](../../hypotheses/PRIMORDIAL_QUBIT.md)

**"E = mγ²."** The theorem gives α = 2γ⟨n_XY⟩, which is linear in γ,
not quadratic. The Lindblad equation is first-order in time (unlike the
wave equation), so the "speed" γ appears to the first power. There is
no E = mc² in the Lindblad cavity.
See: [Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md)

The boundary between theorem and interpretation is sharp:
- **Proven:** α = 2γ⟨n_XY⟩ (algebra, exact, holds for any Hermitian H, real or complex)
- **Observed:** K = γt is an invariant dose ([K-Dosimetry](../../experiments/K_DOSIMETRY.md))
- **Interpreted:** γ is illumination, {I,Z} is the lens, the cavity is an optical instrument

---

## 6. The Absorption Quantum

The absorption theorem reveals that the dephasing spectrum has a natural
unit: **2γ**.

Every Liouvillian eigenvalue has Re(λ) = -2γ⟨n_XY⟩ where ⟨n_XY⟩ ≥ 0.
In the pure dissipator (no Hamiltonian), ⟨n_XY⟩ is an integer and the
rates are exactly:

    0, 2γ, 4γ, 6γ, ..., 2Nγ

This is a ladder with rung spacing 2γ. The number of the rung is n_XY:
the number of X or Y Pauli factors.

| Rung | Rate | n_XY | Content |
|------|------|------|---------|
| 0 | 0 | 0 | All identity and Z: pure structure, immortal |
| 1 | 2γ | 1 | One coherence factor: one photon absorbed per dephasing time |
| 2 | 4γ | 2 | Two coherence factors: two photons |
| ... | ... | ... | ... |
| N | 2Nγ | N | All X and Y: pure coherence, maximum absorption |

When the Hamiltonian is turned on, it mixes the rungs. Modes become
superpositions of different n_XY levels. The mode's position on the ladder
shifts to its average: ⟨n_XY⟩, which can now take any real value in [0, N].
But the fundamental spacing is still 2γ. The Hamiltonian smooths the
ladder into a continuous range, but it cannot change the endpoints (0 and
2Nγ) or the fundamental quantum (2γ).

In the cavity language: 2γ is the absorption cost of one photon. A mode
containing k photons (k coherence factors) pays k × 2γ in absorption rate.
The lens ({I,Z}) is free; it absorbs nothing. Only the light ({X,Y}) pays.

This is why the spectral gap is 2γ, the spectral bandwidth is 2(N-2)γ,
and the palindrome width is 2Nγ: all are integer multiples of the
fundamental quantum 2γ.

**Typed anchor.** The numerical coefficient 2 in the absorption quantum
is the polynomial root a₀ of the Pi2 dyadic ladder
([`Pi2DyadicLadderClaim.Term(0)`](../../compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs)):
the same root d in d² − 2d = 0 that fixes the qubit dimension. Same
anchor as F1 TwoFactor, F50 DecayRateFactor, F66 UpperPoleCoefficient.
The absorption quantum 2γ is the per-mode rate reading of this single
constant: the dyadic ladder's first rung evaluated at illumination γ.

---

## 7. Summary Equation

    Re(λ) = -2γ ⟨n_XY⟩

One equation. Three symbols. The entire absorption spectrum in one line.

---

## Source

- Proof: this document
- Numerical verification: [Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md)
- Simulation: [`simulations/absorption_theorem_discovery.py`](../../simulations/absorption_theorem_discovery.py)
- Palindromic sum rule: [Standing Waves](../../experiments/FACTOR_TWO_STANDING_WAVES.md)
- Weight swap: [Light and Lens](../../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)
- Spectral formulas: [Analytical Formulas](../ANALYTICAL_FORMULAS.md) (Formulas 3, 8, 33, D6)
- XOR drain: [XOR Space](../../experiments/XOR_SPACE.md)
- Palindrome proof: [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)
- Π definition: I ↔ X, Y ↔ iZ, Z ↔ iY (per site, tensor product)

## The rungs, drawn live

![The live Liouvillian spectrum of an N = 5 chain at Q = 1.5 (Symphony export). The eigenvalues' real parts fill the absorption envelope [−2σ, 0] = [−0.5, 0] and align on the rungs −2γ·n (dotted vertical lines).](../../simulations/results/symphony_reel/without_t_axis_spectrum.png)

This proof's quantization, seen at a glance: the real parts (the fade rates) live in [−2σ, 0] and cluster on the rungs −2γ·n. The same figure is the shared anchor of the [F1 palindrome](MIRROR_SYMMETRY_PROOF.md) (the mirror about −σ) and the [F4 kernel](PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md) (the N + 1 frozen modes at Re = 0). Exported by `inspect --root symphony --N 5 --J 0.075 --gamma 0.05 --export`, drawn by `simulations/reel_and_projector.py`.
