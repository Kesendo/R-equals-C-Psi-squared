# The Absorption Theorem

**Status:** Tier 1 derived (analytical proof + per-mode numerical verification at N=2..5, 1,342 modes, ratio 1.000000, CV=0).
**Date:** 2026-04-04 (discovery + proof same day)
**Last refreshed:** 2026-07-19 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.6)
**Statement:** `Re(λ) = −2γ ⟨n_XY⟩` for any Lindblad eigenmode under uniform Z-dephasing
**Typed claim:** [`AbsorptionTheoremClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs) (Tier 1 derived)
**Discovery experiment:** [`ABSORPTION_THEOREM_DISCOVERY.md`](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md)

---

## Abstract

The Absorption Theorem reads the decay rate of any Lindblad eigenmode straight off its Pauli content. Under uniform Z-dephasing at rate γ, for any Hermitian Hamiltonian on any graph,

    Re(λ) = −2γ · ⟨n_XY⟩,

where ⟨n_XY⟩ is the average number of X or Y factors (the "light content") in the eigenmode's Pauli decomposition. The dissipator is built from Z, so it sorts every Pauli letter into two classes: the lens {I, Z}, which the dephasing passes through without cost, and the light {X, Y}, absorbed at exactly 2γ per factor. A mode's lifetime is therefore fixed by one number, how much of its own structure it exposes to the light: eigenmodes carrying no light are immortal (rate 0), eigenmodes of pure {X, Y}^⊗N content die fastest (rate 2Nγ). The number is a property of the eigenmode, not of Pauli content as such: H rotates {I, Z} strings into the light, and for the Heisenberg chain only two {I, Z} strings are frozen outright (I^⊗N and Z^⊗N), the rest of the (N+1)-dimensional F4 kernel being superpositions.

This is the rate-side companion of the F1 palindrome. Π swaps light and lens, so paired eigenvalues satisfy ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N and sit symmetrically about the center, which is the F8 factor-of-two decay law in one line. The recentred face L_D = γ(Q − N·I), with Q = Σ_l Z_l⊗Z_l, shifts the absorption ladder's midpoint to zero and is literally the same diagonal that the F87 girth-moment machinery and the palindrome share: one diagonal, three readings (absorption rates, palindromic pairing, power-sum expansion).

The identity is exact and verified per mode at N=2..5. It is a reading, not a shortcut: ⟨n_XY⟩ is defined from the eigenvector, so evaluating it means already having diagonalized. What it replaces is not the computation but the mystery, by saying what the real part *is*. Where the light content is known in advance, as it is for pure-weight modes, it does predict the rate outright.

## Preface

A guitar string vibrates in modes. Each mode spans the entire string.
How fast a mode fades depends on how much of that mode the damping can
actually grip. Damping at one point of the string grips by position:
a mode with a node under the pad survives, a mode with an antinode
there dies fast. The damping in this document is of the other kind. It
sits on every site equally, and it grips not by where a mode lives but
by what a mode is made of, so what decides a mode's fate is its
composition rather than its shape. (The position-sensitive picture does
return, exactly, once the rates differ from site to site: that is
Theorem 2.)

This document proves the quantum version of that principle. Every
vibration mode in a qubit cavity has a "light content": how much of
the mode consists of the oscillating quantum components (X and Y Pauli
operators) that interact with the external light (γ). The absorption
rate of the mode is exactly twice the light intensity times the light
content:

    Re(λ) = -2γ ⟨n_XY⟩

Components that commute with the dephasing letter (I and Z, the
"lens") are invisible to it and pay nothing. Components that
anticommute (X and Y, the "light") fade. Commuting is the criterion,
not standing still: a Z string moves plenty under H, it simply costs
nothing. Under X-dephasing the free pair is {I, X} instead. The mode's fate is the weighted average,
so an eigenmode carrying no light at all lives forever, and any light
in it costs. Note the word eigenmode: an {I, Z} string that the
Hamiltonian rotates into the light is not one, and it decays.

This single equation, provable in three lines from the structure of
the Lindblad master equation, unifies an entire family of spectral
results: the boundaries (F3), the palindromic sum rule (10,748
pairs), the 2× decay law (F8), the mode classification by Pauli
weight, the N=3 rate ladder (F33), and via later derivation the
weight-1 degeneracy count (F50), the weight-1 dispersion relation
(F2), the GHZ XOR-drain (F22), the n_XY chromaticity (F74), and the
F89 path-k closures including the F89c Hamming-complement pair-sum.
Two relatives sit just outside that list: the spectral gap (D6) and
the dose that rests on it (F55) are *relocated* by the theorem rather
than derived from it, and both carry a coupling-regime condition
(§4.3). The typed
[`AbsorptionTheoremClaim`](../../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs)
holds its own list, which is maintained by hand and does not match this one in
either direction. The live typed edges run the other way: fourteen or more claims,
`StructuralCeilingClaim` and `ClockHandLadderClaim` among them, inject the
absorption claim as a parent.

### Status

| Component | Status | Source |
|-----------|--------|--------|
| Analytical proof | **Proven** | This document, §2, Steps 1-3 |
| Numerical verification | **Verified** (N=2-5, 1,342 modes, CV=0) | [Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md) |
| Consequence 1: Spectral boundaries | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) F3 |
| Consequence 2: Palindromic sum rule | **Derived** | [Standing Waves](../../experiments/FACTOR_TWO_STANDING_WAVES.md) |
| Consequence 3: Spectral gap | **Not derived** (§4.3: 2γ only above an N-dependent Q*_gap) | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) D6 |
| Consequence 4: 2× decay law | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) F8 |
| Consequence 5: Mode classification | **Derived** | [XOR Space](../../experiments/XOR_SPACE.md) |
| Consequence 6: N=3 rate ladder | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) F33 |
| Consequence 7: The recentred face | **Derived** | This document, Section 4.7 |

---

## 1. Setup

Consider an N-qubit system governed by the Lindblad master equation:

    dρ/dt = L(ρ) = -i[H, ρ] + Σ_k γ_k D_{Z_k}(ρ)

where H is any Hermitian Hamiltonian (e.g. the Heisenberg chain),
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

**The coherence basis.** The same count reads off the computational basis, and
the two readings are used interchangeably below. A coherence |A⟩⟨B| decomposes
per site into {I, Z} where A_l = B_l and into {X, Y} where A_l ≠ B_l, so

    n_XY(|A⟩⟨B|) = n_diff(A, B) = Hamming distance of the bit-strings A, B

Write x for a coherence index (a bra-ket label pair) and Δ_l(x) ∈ {0, 1} for the
indicator that bit l differs between its bra and ket labels. Then
n_XY(x) = Σ_l Δ_l(x), and the dissipator is diagonal in this basis too, with
entry −2γ·n_XY(x).

**Notation, including four symbols this repository overloads.**

- **α** indexes a Pauli string (P_α, c_α) in Sections 1 and 2, and denotes a
  decay rate |Re λ| elsewhere, including the Preface, §3's table and the formula
  registry.
  Both usages are established elsewhere in the repo, so both are kept here.
- **Q** is the coupling ratio J/γ₀ everywhere except §4.7, where it is the F87
  dephasing generator Σ_l Z_l ⊗ Z_l. §4.7 says so again at the point of use.
- **Q\*_gap(N)** is not a running variable but the N-dependent threshold value
  of that ratio above which the spectral gap equals 2γ (§4.3). It is distinct
  from the coherence horizon Q\*(N) (see the [Glossary](../GLOSSARY.md), the
  factor-2 convention section).
- **σ** is a single-site Pauli operator in Sections 1 and 2, the lowering
  operator in the compound σ⁻, and the palindrome centre Σγ in §4.7 and the
  figure caption.
- **γ₀** is the same rate as **γ**. The repository writes γ₀ where a single
  uniform rate is meant and γ where the per-site rates γ_k may differ; both
  appear below, and nothing here depends on the distinction.
- **Π** is the palindrome conjugation of F1, acting per site as
  I ↔ X, Y ↔ iZ, Z ↔ iY, extended by tensor product. Each per-site leg swaps a
  lens letter for a light letter, which is the only property this document uses:
  n_XY → N − n_XY. Its role as the eigenvector map is F1's, proved in
  [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md).

---

## 2. The Theorem

**Theorem 1 (Absorption Theorem).** Let L = L_H + L_D be the Liouvillian
of an N-qubit system with any Hermitian Hamiltonian and uniform
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
Checked against a random complex Hermitian H in
[`simulations/popcount_identity_h_independence.py`](../../simulations/popcount_identity_h_independence.py):
the decisive number is that Herm(L) = (L+L†)/2 is the *same* pure Z-dephasing
dissipator for it and for the Hückel H, difference exactly 0.

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
-2γ n_XY(P_α). The Pauli strings all have the same norm, ‖vec(P_α)‖² = d, so
dividing by √d turns the basis orthonormal; L_D is then diagonal with real
entries in an orthonormal basis, hence Hermitian.  ∎

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

**Where the boundary actually runs.** Steps 1 and 3 are not where the physics
is. Together they say Re(λ) = v†Herm(L)v/‖v‖² for any right eigenvector, which
is true of *every* matrix (Bendixson), Lindblad or not, Hermitian H or not:
verified to 1.9·10⁻¹⁵ even under amplitude damping, where the theorem's rate
formula misses by 0.06. All the content sits in Step 2, which says that for
Z-dephasing Herm(L) is *diagonal in the Pauli basis with entries −2γ·n_XY*.

So the boundary is not about which channels the identity covers, it is about
which channels give that diagonal a name. It keeps one whenever the jump
operators are proportional to Pauli strings, because then every conjugation
L_k ρ L_k† sends a Pauli string to ± itself:

- **Z-dephasing**: Re(λ) = −2γ⟨n_XY⟩, this document.
- **X- or Y-dephasing**: the same statement with n_YZ resp. n_XZ (the
  dephase-letter rotation below).
- **Depolarizing** (X, Y, Z jumps at rate γ each per site): still exactly
  diagonal, and Re(λ) = −4γ⟨n_nonI⟩, where n_nonI counts the non-identity
  letters. The light/lens split collapses to identity versus everything.

Two channels are where the naming stops, and only one of them is obvious.

**Amplitude damping**: σ⁻ is not a multiple of a Pauli string, so Herm(L) gains
a non-diagonal part in this basis and eigenmodes mix excitation sectors.  That
correction is F82/F84 above.

**Any jump operator that is a SUM of Pauli strings**, collective dephasing
L = Σ_k Z_k being the case worth naming. It passes every informal description of
"the channel this theorem covers": pure dephasing, built from Z alone, moving no
population. But the cross terms Z_j ρ Z_k survive, so Herm(L) is not diagonal in
the Pauli basis (largest off-diagonal 0.2 at N=3, γ=0.05) and the reading fails.

The fence is therefore *one Pauli string per jump operator*, and it is not the
same as locality. Two measurements at N=3, γ=0.05, J=0.075 make the point from
both sides. A strictly local, uncorrelated channel with jumps (X_k + Z_k)/√2
**breaks** it (off-diagonal 0.05): each jump is a sum on its own site. A
genuinely two-site channel with jumps Z_k Z_{k+1} **keeps** Herm(L) exactly
diagonal (off-diagonal 0.0) and obeys a theorem of the same form, with the count
running over bonds instead of sites. Correlated is fine; a sum is not.

In both cases the Rayleigh identity is untouched; what is lost is the
closed-form reading of its value.

### Extensions

Four readings are implicit in the three-step proof and deserve their own
statements: the vector form, the two-sided reading, the projector form, and the
dephase-letter rotation. None changes the proof; each is the same Rayleigh
quotient read in a different coordinate, so each is re-derivable from Steps 1-3
alone. Where a reading is machine-checked, the test is named at the reading.

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
or redistributed as the coupling ratio Q = J/γ₀ moves, and the rate moves
*only* through this weighted share. The typed carrier-is-a-vector node on
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
and even turn negative; it fails the absorption cross-check by 1.7·10⁻¹ on the
flat-bulk-edge γ-profile of
[`BirthCanalSurfaceWitness`](../../compute/RCPsiSquared.Diagnostics/Foundation/BirthCanalSurfaceWitness.cs),
recorded in
[`SlowLightDistribution.cs`](../../compute/RCPsiSquared.Diagnostics/Ptf/SlowLightDistribution.cs)):

    light_l(V) = Tr(Π_V Δ_l) / dim V  ∈ [0, 1],

with Δ_l the site-l disagreement projector. Absorption-exactness follows from
a block-triangular argument: in an orthonormal basis adapted to the invariant
subspace, L is block upper-triangular, so Re Tr(Π_V L) = Σ_k Re λ_k, and the
anti-Hermitian Hamiltonian part drops from the real trace, leaving exactly the
γ-weighted light of Π_V. At cluster dimension 1 this reduces identically to the
per-eigenvector Rayleigh light (gated at 10⁻¹⁰ by
`SlowLightDistributionTests.Projector_MatchesEigenvectorAverage_WhereDegeneracyIsOne`).
This is the correct object for the flow kernel's degenerate slow carriers, and
is what the `SlowLightDistribution` primitive computes.

**Remark (the dephase-letter rotation).** Nothing in Steps 1-3 privileges Z.
For X-dephasing the dissipator is diagonal with light = n_YZ (the letters
anticommuting with X); for Y-dephasing, n_XZ. Light is always "the letters the
dephasing letter refuses to commute with." The three readings are permuted by
the same Klein V₄ that permutes the palindrome family
([`Pi2KleinV4DephaseSwapGroup`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs),
which types the Π-family swap; the light-count rotation is stated here, not
there).

**Remark (the recurring four).** Two different objects in this repository both
read "4" at d=2, and they must not be interchanged. The **rung-2 four** is
2γ·2, the second rung of the §6 ladder and the centre of the HD=(1,3) palindrome
pair, (2γ+6γ)/2 = 4γ. It is Hamming-indexed, and it splits per site: |00⟩⟨11| is an exact N=2 eigenmode of the
XY-chain Liouvillian at rate −2(γ₁+γ₂), residual 0.0
([`simulations/at_rung2_per_site_split.py`](../../simulations/at_rung2_per_site_split.py)).
Every *dynamical* 4γ₀ in the repository is this one object, t_peak = 1/(4γ₀)
([`TPeakLaw`](../../compute/RCPsiSquared.Core/F86/TPeakLaw.cs)) among them. The
**discriminant four** is dimensionless, d² on the Pi2 dyadic ladder, and rides d
to 9 at the qutrit. Cite the rung-2 four for rates and carriers, the discriminant
four for squared half-gaps, never interchangeably. The full qudit sorting of the
fours, with the F121 cap and the F122 ceiling, belongs to those claims rather
than to this ladder; see
[`simulations/qudit_g2_split.py`](../../simulations/qudit_g2_split.py).

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

The 1,342 are the modes on which the ratio is defined: 4²+4³+4⁴+4⁵ = 1,360
minus the 18 = Σ(N+1) kernel modes, where α/(2γ⟨n_XY⟩) is 0/0. The kernel is
not an exception to the theorem, it is the ⟨n_XY⟩ = 0 end of it.

**Source:** [`simulations/absorption_theorem_discovery.py`](../../simulations/absorption_theorem_discovery.py).
The N and γ rows are per-mode (Step 6); the J row is the median of the
per-*pair* ratio across the five J values (Step 4).

**Hardware confirmation, and what it does not reach.** Single-qubit tomography
on IBM Torino Q52 (25 time snapshots, 0-895 μs) reproduced the theorem at 3%
deviation under the free-evolution T2* baseline: absorption ratio excess/(2γ) =
1.03, with γ* fitted from the coherence envelope. The echo-refocused γ_echo
underestimates the dephasing this qubit actually saw by 6.2× because
low-frequency 1/f noise is filtered out by Hahn echo; the T2* baseline is the
correct one for free-evolution tomography. See
[`IBM_ABSORPTION_THEOREM.md`](../../experiments/IBM_ABSORPTION_THEOREM.md).

The scope of that number is narrow, and worth stating so nobody reads it as
more. At N=1 with n_XY = 1 the theorem says the coherence decays at 2γ, while γ*
is itself extracted from that same coherence envelope, so the 1.03 measures the
consistency of two fits to one decay, not the ladder. The content of the theorem
is the *spacing*: a Hamming-2 coherence must decay at exactly twice the rate of a
Hamming-1 coherence on the same device. That measurement needs N ≥ 2 and has not
been made.

**The regime map.** Where each corollary of §4 holds, and where it stops, is a
separate question from the identity itself, and every number §4 quotes about it
is asserted in
[`simulations/absorption_ladder_regimes.py`](../../simulations/absorption_ladder_regimes.py):
the N=3 ladder with multiplicities, the J-dependence of the two fractional
rates, the two gap regimes, the palindromic band erosion at N=5, the kernel
dimension against the {I,Z} sector size, the per-coherence residuals, and the
depolarizing channel. Run it to reproduce §4 from scratch.

**C# test-gating.** The per-mode identity is asserted in the typed compute
layer:
[`F8PartnerLightComplementarityTests`](../../compute/RCPsiSquared.Diagnostics.Tests/Ptf/F8PartnerLightComplementarityTests.cs)
checks |Re λ + 2γ·light(v)| < 10⁻⁹ for all 256 eigenmodes of the N=4 XY chain,
the complete palindrome pairing with light_s + light_f = N per pair (128 pairs),
and the standing-wave null. That test covers the right-eigenvector reading; the
two-sided reading of Section 2 is not yet gated anywhere.

---

## 4. Consequences

**Scope of this section.** Section 2's theorem holds for any Hermitian H. The
consequences below do not: they describe the number-conserving XY/Heisenberg
family the rest of this repository works in, and they inherit that family's
structure through F4 (kernel dimension) and F1 (the pairing). The theorem
converts each of those structural facts into a statement about rates; it does
not supply them. For a generic real symmetric H at N=3, γ=0.05, the kernel is
one-dimensional, no mode reaches 2Nγ (4.20γ on the draw the regime script fixes),
the rates form a continuum with no 2γ spacing, and the palindrome is absent. By
contrast a pure Ising ZZ Hamiltonian freezes all 2^N of the {I,Z} strings, so the
kernel is 8-dimensional at N=3 rather than 4. Nothing there
contradicts Section 2, and every rate still equals −2γ⟨n_XY⟩ exactly. Read
Section 4 as the ladder's shape *in this family*.

### 4.1 Spectral Boundaries (F3)

**Previously:** The decay rate spectrum satisfies min = 2γ, max = 2(N-1)γ,
bandwidth = 2(N-2)γ, *for the generic band*: the N+1 kernel modes at 0 and the
N+1 XOR-drain modes at 2Nγ sit outside it at both ends.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), F3

**Now a corollary.** The light count n_XY is an integer in {0, 1, ..., N}, and
⟨n_XY⟩ takes values in [0, N]. The two ends of the full spectrum:

    ⟨n_XY⟩ = 0: rate = 0 (the F4 kernel, N+1 modes)
    ⟨n_XY⟩ = N: rate = 2Nγ = 2Σγ (the XOR drain, N+1 modes)

so the full spectral width is 2Nγ, and F3's band is what remains after both
extremes are removed: from 2γ up to 2(N-1)γ, width 2(N-2)γ. Keeping the two
conventions apart matters, because 2(N-1)γ and 2Nγ are both correct answers to
"the maximum" for different sets, and at N=2 the band degenerates to the single
rate 2γ (width 0).

F3's band edges are exact where they hold: at N=3 and N=4 in the canonical
regime the 2γ and 2(N-1)γ levels are attained to machine precision (each with
multiplicity 14 at N=3), with ⟨n_XY⟩ exactly 1 and exactly N-1, no approximation.
But they are not a universal bound, and the theorem is what explains why. Nothing
in Re(λ) = −2γ⟨n_XY⟩ forbids ⟨n_XY⟩ from falling between 0 and 1: strong
Hamiltonian mixing spreads a mode's light over sectors, and modes then leave the
band. On the **Heisenberg** chain at N=5, γ=0.05, J=0.075 (Q = 1.5, below that
chain's threshold; the XY chain at the same point shows no erosion at all, §4.3)
the extreme band rates are 1.2754γ and 8.7246γ,
outside [2γ, 8γ] at *both* ends at once, and their sum is exactly 10γ = 2Nγ.
That is the signature: modes leave the band in palindromic pairs (§4.2), so the
band erodes symmetrically and the full width 2Nγ is untouched. The band edges
sharpen back to 2γ and 2(N-1)γ as J/γ grows and modes relocalize on pure weight
sectors.

### 4.2 Palindromic Sum Rule (10,748 pairs)

**Previously:** For every palindromic pair, Re(λ_fast) + Re(λ_slow) = -2Σγ,
verified for 10,748 pairs across N=2..7 with zero exceptions (a wider range
than this document's own N=2..5 per-mode verification).
**Source:** [Standing Waves](../../experiments/FACTOR_TWO_STANDING_WAVES.md)

**Now a one-line corollary.** The palindromic weight swap
([Light and Lens](../../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)) gives:

    ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N

This needs two inputs, not one. Π maps weight sector k to sector N-k, which is
the per-site statement checked at the end of this document; and Πv is the
partner eigenvector, which is F1's content, proved in
[Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md). The weight map alone does not
suffice: it is F1 that puts the two modes in the same pair.

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

**Not derived.** The theorem turns the gap into a question about light and then
stops:

    Δ = 2γ × min{⟨n_XY⟩ : ⟨n_XY⟩ > 0}

It places no lower bound on that minimum, so it does not by itself yield 2γ.
What the minimum is depends on H, and it has two regimes.

Above a threshold in Q = J/γ the slowest non-kernel modes sit on pure weight-1
strings, ⟨n_XY⟩ is exactly 1, and the gap is exactly 2γ: the cost of a single X
or Y factor. The threshold is not universal and grows with N. Located by
bisection on the **Heisenberg** chain, Q*_gap = 0.5000 at N=2, 0.8002 at N=3,
1.3422 at N=4, 1.8194 at N=5. No closed form is claimed; the successive gaps
are 0.300, 0.542 and 0.477, so it is not linear in N over this range.

The threshold belongs to the Hamiltonian, not to N alone. On the **XY** chain
the same bisection gives 0.7071 at N=3 (that is 1/√2 to five places), 0.9393 at
N=4 and 1.1861 at N=5, all below the canonical Q = 1.5, and the XY band at N=5,
Q = 1.5 is exactly [2γ, 8γ] with no erosion at all. So §4.1's erosion is a fact
about the ZZ term, not about reaching N=5: on Heisenberg the canonical regime
sits above threshold at N=3 and N=4 and below it at N=5, while on XY it stays
above throughout. Note also that §3's C# gate runs the XY chain while §4's
regime numbers are Heisenberg.

As J/γ → 0 the gap follows a different law. The {I,Z} strings that H lifts out
of the J=0 kernel acquire light only perturbatively, ⟨n_XY⟩ ~ (J/γ)², and

    Δ → 2(1 − cos(π/N)) · 2J²/γ        (asymptotic, not a second regime)

the second factor being the Zeno-suppressed rate and the first the spectral gap
of the path-graph Laplacian on N sites, the classical hopping generator that
survives strong dephasing, not a coherent band edge. Measured against 2J²/γ at
γ=0.3, J=0.001: 2.000022 at N=2, 1.000003 at N=3, 0.585784 at N=4, 0.381964 at
N=5, against 2(1−cos(π/N)) = 2, 1, 0.585786, 0.381966. The bare form Δ = 2J²/γ
is the N=3 case only. At N=3, γ=0.3, J=0.001 the gap is 1.1·10⁻⁵ of 2γ: five
orders below the quoted value.

This is an asymptote, not the other half of a partition of the Q axis. All
error percentages here are measured against the true gap, (asymptote − gap)/gap.
The asymptote is accurate to about 1% only for Q ≲ 0.1, overshoots by +8% at
Q = 0.5 and +35% at Q = 1.5 (N=5), and crosses the 10% line near Q ≈ 0.57.
From there up to a narrow sliver just below the threshold (the 2γ form itself
comes within 10% only above Q ≈ 1.75 at N=5, against Q*_gap = 1.82) neither
closed form is good to 10%, and the gap is simply a number one computes: at
N=5, Q = 1.5 the asymptotic form would give 1.7188γ against the true 1.2754γ
that §4.1 quotes. The overshoot is not monotone in Q: it peaks near Q ≈ 1.55
(+35%) and falls back to +29% at the threshold itself, where the true gap is
climbing to meet 2γ, so the two closed forms approach each other from a
worst point in between rather than trading off cleanly.

This is Zeno suppression, read from the light side. Strong dephasing freezes the
transport that would give a mode its light, and a mode with almost no light
almost cannot decay. The D6 value 2γ is the strong-coupling face of the gap, not
its universal value.

**Which modes fall, and which never do.** The erosion is confined to one sector,
and this bounds how far the correction reaches. Sorting the spectrum by
|Δpopcount|, the difference in excitation number between a coherence's bra and
ket labels, everything below 2γ sits at |Δpopcount| = 0: the population modes,
which is exactly where dephasing freezes transport. Every coherence sector keeps
its minimum at the ladder rung, at every coupling, and this one is derived, not
sampled:

For a number-conserving H the Liouvillian is exactly block-diagonal in
|Δpopcount| (verified entry-wise: the largest matrix element between different
|Δpopcount| sectors is 0.0 at N = 3 and N = 4), so each sector's eigenmodes are
superpositions of basis coherences |A⟩⟨B| drawn from that sector alone. Every
such coherence has n_XY = Hamming(A, B) ≥ |popcount(A) − popcount(B)| = k, since
flipping k excitations takes at least k bit flips. By the theorem the rate of a
mode is 2γ⟨n_XY⟩, and ⟨n_XY⟩ is a convex combination of those Hamming distances,
hence ≥ k. So a |Δpopcount| = k mode decays at least as fast as 2kγ, at every
coupling and for every number-conserving H. ∎

The measurement agrees: at Q = 0.2, 0.5 and 1.5 for N = 3, 4, 5, the
|Δpopcount| = 1 minimum is 2.000000γ and the |Δpopcount| = 2 minimum is
4.000000γ, without exception, while the population minimum runs 0.0301γ to
2.4607γ over the same range. Note the bound is saturated in these sectors, which
the argument alone does not force.

So two statements that look like one must be kept apart. *The spectral gap is
2γ* is regime-bounded. *A coherence decays no slower than 2γ* is not: it holds
at every coupling. Anything reasoning about coherent windows, transfer times,
echo cycles or T2 is using the second statement and is untouched by the
threshold; only claims about the slowest mode of the full Liouvillian, and about
mixing to the steady manifold, need the Q*_gap condition.

### 4.4 The 2× Decay Law (F8)

**Previously:** Unpaired modes (at the palindrome extremes) decay at rate
2Nγ, while paired modes average rate Nγ. The ratio is exactly 2.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), F8;
[Energy Partition](../../hypotheses/ENERGY_PARTITION.md)

**Now explained.** Two words have to be kept apart first, because F8's phrasing
compresses them. Every mode is palindromically paired; the drain modes at 2Nγ
are the partners of the kernel modes at 0 (0 + 2Nγ = 2Σγ). What F8 calls
"unpaired" means *not self-paired*: a mode whose Π-partner is a different mode,
which is why it can sit at an extreme. A **self-paired** mode is its own partner
and therefore has nowhere to sit but the fixed point of the reflection.

That fixed point is the palindrome center α = Σγ = Nγ, corresponding to
⟨n_XY⟩ = N/2, an equal mix of light and lens. The full range spans 0 (lens) to
2Σγ (light).

The ratio 2 is the ratio of the maximum (2Σγ) to the center (Σγ), and it is not
a dynamical law. Two different averages meet at that centre and it is worth
keeping them apart. The mean over the *whole spectrum* is Σγ for any Hermitian
H, with no symmetry needed: Tr(L_H) = 0 and Tr(L_D) = −γN·d², since summing
n_XY over all 4^N strings gives N·d²/2. Verified at 3.0000000000·γ and
4.0000000000·γ for a random complex Hermitian H at N=3 and N=4. F8's statement
is the stronger one, that *each palindromic pair* averages Σγ, and that does need
F1. The two coincide exactly when the spectrum is palindromic. The numerator,
max = 2Nγ, belongs to the number-conserving family either way, so for a generic
H the ratio is not 2 at all.

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

The first and last rows name eigenmodes, not Pauli sectors. The {I,Z} sector
holds 2^N strings, and on the Heisenberg chain exactly two of them are frozen
outright, I^⊗N and Z^⊗N; the F4 kernel is an (N+1)-dimensional subspace whose
remaining directions are superpositions, not single strings. H rotates the other
strings into the light and they decay. The pure string Z₁ at N=3 is the plain
example, retaining only its 1/N share of the conserved total Z:
⟨Z₁(t)|Z₁(0)⟩ → 1/3 exactly. Immortality is a property of ⟨n_XY⟩ = 0, which H
must permit, not of carrying {I,Z} letters.

**Per-coherence reading on the computational basis.** By the coherence-basis
identity of §1, n_XY(|A⟩⟨B|) = n_diff(A, B), so the dissipator assigns the
coherence |A⟩⟨B| the diagonal entry −2γ·n_diff(A, B).

This is a statement about the diagonal of L_D, and it becomes a *decay rate*
only when |A⟩⟨B| is an eigenmode of the full L, which requires |A⟩ and |B⟩ to be
eigenstates of H. That holds for |00⟩⟨11| at N=2 (residual exactly 0.0, XY and
Heisenberg) and fails generically: |001⟩⟨010| at N=3 has residual 0.30 at
γ=0.05, J=0.075 and 4.0 at J=1. Where H does not fix the two labels, the
coherence is a superposition of eigenmodes and 2γ·n_diff is its initial slope,
not its rate.
The **Hamming-complement pair-sum** (F89c) is the immediate corollary:
the column bit-flip ρ[a, b] → ρ[a, b̄] satisfies n_diff(a, b) +
n_diff(a, b̄) = N, hence:

    α(|a⟩⟨b|) + α(|a⟩⟨b̄|) = 2γN = 2Σγ

The two coherences related by complementing one column-label always sit at
diagonal entries summing to exactly the spectral maximum 2Σγ. This is the
single-block reading of the palindromic sum rule §4.2, and it carries the same
caveat: it is an identity on the dissipator's diagonal, exact at any H, and a
statement about decay rates only where those coherences are eigenmodes.

### 4.6 The N=3 Rate Ladder (F33)

**Previously:** The N=3 Heisenberg chain has three distinct non-trivial
decay rates: 2γ, 8γ/3, 10γ/3.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), F33

**Now explained, and split into two kinds.** The theorem reads every rate as
light, and doing so separates the N=3 spectrum into rates that are exact at
every coupling and rates that are only a limit.

*Exact at every J*, the four pure-weight rungs:

    rate = 0      → ⟨n_XY⟩ = 0     multiplicity 4   (the F4 kernel)
    rate = 2γ     → ⟨n_XY⟩ = 1     multiplicity 14  (pure weight-1)
    rate = 4γ     → ⟨n_XY⟩ = 2     multiplicity 14  (pure weight-2)
    rate = 6γ     → ⟨n_XY⟩ = 3     multiplicity 4   (the XOR drain)

*A strong-coupling limit only:* 8γ/3 and 10γ/3. Each is not one level but a
triple, and the triples merge only as J/γ → ∞:

    J/γ = 1.5 (canonical):  2.4607γ  2.6040γ  2.6980γ
    J/γ = 20:               2.6655γ  2.6663γ  2.6668γ
    J/γ → ∞:                8γ/3 = 2.6667γ

Each band carries multiplicity 14, so the ladder closes: 36 modes on the four
pure rungs plus 28 in the two mixed bands is the full 64. The band's lowest
level approaches its limit as 0.46·(γ/J)², and the band's own spread closes as
0.53·(γ/J)², both flat over three decades. So 8γ/3 is a J/γ → ∞
face of the N=3 spectrum, not an exact rational rate, and in the regime this
repository actually works in (Q = 1.5) it is off by up to 7.7%.

Both kinds obey the theorem exactly; that is the point. Non-integer ⟨n_XY⟩
arises from Hamiltonian mixing across weight sectors, and the mixing weights
themselves depend on J, which is precisely why the resulting rate does too. What
is J-independent is the ladder of pure-weight levels 2γ·n, because a mode that
sits on one weight sector has nothing left to mix.

The mixing is parity-preserving, which is worth stating because the obvious
guess is wrong. XX+YY+ZZ moves XY-weight in steps of two: the parity operator
(−1)^{n_XY} commutes with L at every J (residual 1·10⁻¹⁶), so a Δw = 1 mixture
is structurally forbidden. The 14-dimensional cluster splits into 10
even-parity and 4 odd-parity eigenvectors, a w=0/w=2 family and a w=1/w=3
family, never w=1 with w=2.

The weight shares themselves follow the same limit discipline as the rate.
Read basis-freely through the orthogonal projector onto the cluster, ×21:

    shares over w = 0, 1, 2, 3, scaled by 21:

    Q = 1.5 (canonical):  5.309, 5.309, 9.691, 0.691   ⟨n_XY⟩ = 1.2745
    Q → ∞:                5,     5,     10,    1       ⟨n_XY⟩ = 4/3

This is forced, not a separate fact: ⟨n_XY⟩ = rate/2γ, so if 8γ/3 is a limit
then 4/3 is the same limit. To check it against the rates above, weight the
three canonical levels by their multiplicities within the 14, which are 8, 2 and
4: (8·2.4607 + 2·2.6040 + 4·2.6980)/14 = 2.5490, and half of that is 1.2745.
The unweighted mean of the three levels is 1.2938 and is not the light content. The 10γ/3 cluster is the mirror, (1, 10, 5, 5)/21
and 5/3, in the same limit.

### 4.7 The recentred face: one diagonal, three pillars

Per the notation block in §1: in this section **Q** is the F87 dephasing
generator Σ_l Z_l ⊗ Z_l, not the coupling ratio, and **σ** is the palindrome
centre Σγ, not a single-site Pauli operator.

The dissipator's diagonal is, entry for entry, the same object the F87
palindrome machinery calls its dephasing generator. With w(x) = popcount of
the bra-ket difference at coherence index x, the uniform dissipator has
diagonal −2γw(x), and the windowed-converse generator Q = Σ_l Z_l ⊗ Z_l has
diagonal Q_xx = N − 2w(x), since Q_xx = Σ_l (−1)^{Δ_l(x)}. So, entry for entry:

    L_D = γ·(Q − N·I)        and        M := L + γN·I = L_H + γQ

The shift by the palindrome centre σ = Σγ = Nγ that recentres L into M
([the F87 windowed converse](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §1)
is the shift that moves the absorption ladder's midpoint ⟨n_XY⟩ = N/2 to
zero. Three pillars stand on this one diagonal, each reading it in its own
coordinate:

- **the absorption ladder** (this proof): Re(λ) = −2γ⟨n_XY⟩, rates as light;
- **the palindrome** (F1/F8), which stands on this diagonal only in part: the
  reflection λ → −λ* − 2σ of spec(L) about −Nγ is equivalent to spec(M)
  symmetric about 0, and partners then carry complementary light. That the
  partner is Πv is F1's, a statement about L_H (§4.2), not about this diagonal;
  ⟨n_XY⟩_s + ⟨n_XY⟩_f = N (the §4.2 sum rule, test-gated per mode in C#,
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

**"Absorption" is a name, not a process.** The word is in the title and runs
through this document, and nothing is absorbed. The dissipator annihilates every
diagonal ρ exactly, so it moves no population: what decays is coherence, and
n_XY counts coherence factors, not excitations. 2γ is a decoherence rate. The
name records that the rate is *read off* a content, the way an absorption
coefficient is; it does not claim a photon.

Note what this argument does *not* say. Dephasing conserves populations, but it
does not conserve energy unless H commutes with every Z_k: Tr(H·D[ρ]) vanishes identically on Ising ZZ, at every state, because D
preserves the Z-diagonal; on the Heisenberg chain it does not, ranging over
about ±0.04 for Haar-random pure states at N=3, γ=0.05, J=0.075, and ±0.5 at
J=1. The chain
heats toward infinite temperature. So the case against the photon reading rests
on population preservation alone, not on energy conservation.

Two further limits are structural rather than interpretive.

**Re λ is the asymptote of an envelope, not the envelope.** L is non-normal, so
‖e^{Lt}‖ can grow substantially while every Re λ is negative, and at a defective
point the decay carries a polynomial prefactor t^k·e^{Re λ·t}. The theorem
covers the whole spectrum, since every eigenvalue has at least one genuine
eigenvector, but a mode's light fixes where it ends up, not the path it takes.

**It constrains only Re λ.** Nothing here says anything about Im λ, about how
the rates are distributed, or about the structure of the spectrum beyond its
real parts.

The boundary between theorem and interpretation is sharp:
- **Proven:** α = 2γ⟨n_XY⟩ for Z-dephasing, any Hermitian H, real or complex.
  Other Pauli-string jump operators keep the form but change the count
  (n_YZ, n_XZ, or 4γ·n_nonI for depolarizing), per §2.
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

When the Hamiltonian is turned on, it mixes the rungs. Modes become
superpositions of different n_XY levels. The mode's position on the ladder
shifts to its average: ⟨n_XY⟩, which can now take any real value in [0, N].
But the rungs themselves stay where they are: 2γ is a property of the
dissipator, and no Hamiltonian moves it.

What a Hamiltonian *can* do is decide which rungs get occupied. Within the
number-conserving family the endpoints 0 and 2Nγ are always reached, so the
ladder is spanned end to end; outside it they need not be (§4's scope note).
The quantum 2γ is universal; the occupancy of the ladder is the Hamiltonian's
business.

Each X or Y factor costs exactly 2γ, each I or Z factor costs nothing; a mode
carrying k coherence factors pays k × 2γ. Keep the accounting in those words and
not in photons: what is counted is coherence, not excitation (§5). The channel
that would earn the photon language is amplitude damping, and that is exactly
the channel this theorem does not cover.

The three widths follow, each on its own set: the full palindrome width is 2Nγ,
F3's generic band is 2(N-2)γ wide once both extremes are removed, and the gap
above the kernel is 2γ wherever J/γ is above threshold (§4.3). All are integer
multiples of the fundamental quantum 2γ.

**Where the 2 comes from.** Directly, from Step 2: Z anticommutes with X and Y,
so ZXZ − X = −2X. The 2 is the spread of the spectrum of the conjugation map,
|1 − (−1)|, and nothing about the dimension enters.

**Typed anchor, and an open question about it.** The repository also carries the
2 as the polynomial root a₀ of the Pi2 dyadic ladder
([`Pi2DyadicLadderClaim.Term(0)`](../../compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs)),
the root d of d² − 2d = 0 that fixes the qubit dimension, shared with F1
TwoFactor, F50 DecayRateFactor and F66 UpperPoleCoefficient. Those four
coefficients do all equal 2, and the anchor records that.

Whether they equal 2 for the *same reason* is not settled here. The direct
derivation above needs no dimension, which is reason enough to doubt the
identification, but the doubt cannot be settled by asking "does the rung stay 2
at d=3", because that question has no unique answer. A Hamming-indexed ladder
d·γ·(differing sites) is still Hamming-indexed with a d-dependent spacing, so
"indexed by Hamming distance" and "spacing independent of d" are not the same
claim, and three natural qutrit dephasers give three different rungs: a single
clock generator Z = diag(1, ω, ω²) gives 1 − cos(2πm/d), i.e. 1.5; all powers
Z^k at rate γ each give d = 3; the equidistant full-Cartan construction used in
`qudit_g2_split.py` gives 2 by construction. At d=4 the single-clock rates are
not even monotone in m, so "the rung spacing" is not well posed there.

So the open question is not a measurement but a modeling choice: *which qudit
channel is the right analogue of single-qubit Z-dephasing*. Until that is fixed,
the shared anchor should be read as recording that four coefficients equal 2,
not as a claim that one causes the others. Moving it would touch F1, F50 and
F66, so it stays as it is and this paragraph states the reservation.

---

## Source

- Proof: this document
- Numerical verification: [Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md)
- Simulation: [`simulations/absorption_theorem_discovery.py`](../../simulations/absorption_theorem_discovery.py)
- Regime map for Section 4: [`simulations/absorption_ladder_regimes.py`](../../simulations/absorption_ladder_regimes.py)
- Palindromic sum rule: [Standing Waves](../../experiments/FACTOR_TWO_STANDING_WAVES.md)
- Weight swap: [Light and Lens](../../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)
- Spectral formulas: [Analytical Formulas](../ANALYTICAL_FORMULAS.md) (Formulas 3, 8, 33, D6)
- XOR drain: [XOR Space](../../experiments/XOR_SPACE.md)
- Palindrome proof: [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)
- Π definition: see the notation block in Section 1

## The rungs, drawn live

![The live Liouvillian spectrum of an N = 5 chain at Q = 1.5 (Symphony export). The eigenvalues' real parts fill the absorption envelope [−2σ, 0] = [−0.5, 0] and align on the rungs −2γ·n (dotted vertical lines).](../../simulations/results/symphony_reel/without_t_axis_spectrum.png)

This proof's quantization, seen at a glance: the real parts (the fade rates) live in [−2σ, 0] and cluster on the rungs −2γ·n. The same figure is the shared anchor of the [F1 palindrome](MIRROR_SYMMETRY_PROOF.md) (the mirror about −σ) and the [F4 kernel](PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md) (the N + 1 frozen modes at Re = 0). Exported by `inspect --root symphony --N 5 --J 0.075 --gamma 0.05 --export`, drawn by `simulations/reel_and_projector.py`.
