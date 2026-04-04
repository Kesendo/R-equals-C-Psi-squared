# The Absorption Theorem

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
the Lindblad master equation, unifies six previously separate results:
the spectral boundaries (Formula 3), the palindromic sum rule (10,748
pairs), the spectral gap (D6), the 2× decay law (Formula 8), the mode
classification by Pauli weight, and the N=3 exact rates (Formula 33).

The theorem was discovered computationally on April 4, 2026 (ratio
α/(2γ⟨n_XY⟩) = 1.000000 for 1,343 modes, CV = 0), then proven
analytically the same day.

### Status

| Component | Status | Source |
|-----------|--------|--------|
| Analytical proof | **Proven** | This document, Section 2 |
| Numerical verification | **Verified** (N=2-5, 1,343 modes, CV=0) | [Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md) |
| Consequence 1: Spectral boundaries | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) Formula 3 |
| Consequence 2: Palindromic sum rule | **Derived** | [Standing Waves](../../experiments/FACTOR_TWO_STANDING_WAVES.md) |
| Consequence 3: Spectral gap | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) D6 |
| Consequence 4: 2× decay law | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) Formula 8 |
| Consequence 5: Mode classification | **Derived** | [XOR Space](../../experiments/XOR_SPACE.md) |
| Consequence 6: N=3 exact rates | **Derived** | [Analytical Formulas](../ANALYTICAL_FORMULAS.md) Formula 33 |

---

## 1. Setup

Consider an N-qubit system governed by the Lindblad master equation:

    dρ/dt = L(ρ) = -i[H, ρ] + Σ_k γ_k D[Z_k](ρ)

where H is a real Hermitian Hamiltonian (e.g. the Heisenberg chain),
Z_k is the Pauli-Z operator on site k, and the Lindblad dissipator is:

    D[Z_k](ρ) = Z_k ρ Z_k - ρ

The Liouvillian L acts on the space of density matrices (dimension d²
where d = 2^N). We work with the vectorized form: L is a d² × d² matrix
acting on vec(ρ).

**Decomposition.** The Liouvillian splits as L = L_H + L_D where:

    L_H = -i[H, ·]         (Hamiltonian part)
    L_D = Σ_k γ_k D[Z_k]   (dissipative part)

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

**Step 1. L_H is anti-Hermitian.**

For a real Hermitian Hamiltonian H (H = H† = H̄ = H^T), the superoperator
L_H = -i[H, ·] satisfies:

    L_H† = -L_H

*Proof of Step 1.* In the vectorized representation:

    L_H = -i(H ⊗ I - I ⊗ H^T) = -i(H ⊗ I - I ⊗ H)

where the second equality uses H^T = H (real). Taking the adjoint:

    L_H† = +i(H† ⊗ I - I ⊗ H†) = +i(H ⊗ I - I ⊗ H) = -L_H  ∎

**Consequence of Step 1.** For any vector v:

    v†L_H v is purely imaginary

*Proof.* (v†L_H v)* = v†L_H† v = v†(-L_H)v = -(v†L_H v). A complex
number equal to the negative of its conjugate is purely imaginary.  ∎

**Step 2. L_D is Hermitian and diagonal in the Pauli basis with
eigenvalues -2γ n_XY.**

The dissipator for Z-dephasing at rate γ on site k acts on a single-site
Pauli operator σ as:

    D[Z_k](σ) = Z_k σ Z_k - σ

Since Z commutes with I and Z, and anticommutes with X and Y:

    D[Z](I) = ZIZ - I = I - I = 0
    D[Z](Z) = ZZZ - Z = Z - Z = 0
    D[Z](X) = ZXZ - X = -X - X = -2X
    D[Z](Y) = ZYZ - Y = -Y - Y = -2Y

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
1. H is real Hermitian (so L_H is anti-Hermitian)
2. The jump operators are Z_k (so L_D is diagonal in the Pauli basis)

It does **not** use:
- The specific form of H (Heisenberg, Ising, XY, any Hamiltonian works)
- The chain topology (any graph works)
- Uniform dephasing (site-dependent γ_k works: replace 2γ n_XY with
  Σ_k 2γ_k × [σ_k ∈ {X,Y}])
- Any property of the eigenvalue spectrum beyond existence

The theorem holds for **any real Hermitian Hamiltonian under Z-dephasing**.

For complex Hermitian Hamiltonians (e.g. with Dzyaloshinskii-Moriya
interactions), L_H is not anti-Hermitian, and the theorem may not hold.
This boundary is precise and testable.

---

## 3. Numerical Verification

The theorem was verified computationally in
[Absorption Theorem Discovery](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md):

| Parameter | Range | Ratio α/(2γ⟨n_XY⟩) | CV |
|-----------|-------|---------------------|-----|
| N (chain length) | 2, 3, 4, 5 | 1.000000 | 0.0000 |
| γ (dephasing rate) | 0.01, 0.05, 0.1, 0.5, 1.0 | 1.000000 | 0.0000 |
| J (coupling strength) | 0.1, 0.5, 1.0, 2.0, 5.0 | 1.000000 | 0.0000 |
| Total modes tested | 1,343 | 1.000000 | 0.0000 |

The ratio equals 1 to 14 decimal places across all parameters.
No exceptions. Zero coefficient of variation.

**Source:** [`simulations/absorption_theorem_discovery.py`](../../simulations/absorption_theorem_discovery.py), Step 6

---

## 4. Consequences

### 4.1 Spectral Boundaries (Formula 3)

**Previously:** The decay rate spectrum satisfies min = 2γ, max = 2(N-1)γ,
bandwidth = 2(N-2)γ.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), Formula 3

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

### 4.4 The 2× Decay Law (Formula 8)

**Previously:** Unpaired modes (at the palindrome extremes) decay at rate
2Nγ, while paired modes average rate Nγ. The ratio is exactly 2.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), Formula 8;
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

### 4.6 N=3 Exact Rates (Formula 33)

**Previously:** The N=3 Heisenberg chain has three distinct non-trivial
decay rates: 2γ, 8γ/3, 10γ/3.
**Source:** [Analytical Formulas](../ANALYTICAL_FORMULAS.md), Formula 33

**Now explained.** The fractional rates correspond to fractional ⟨n_XY⟩:

    rate = 2γ      → ⟨n_XY⟩ = 1     (pure weight-1 modes)
    rate = 8γ/3    → ⟨n_XY⟩ = 4/3   (Hamiltonian mix of w=1 and w=2)
    rate = 10γ/3   → ⟨n_XY⟩ = 5/3   (Hamiltonian mix of w=1 and w=2)

The Hamiltonian mixes Pauli strings of different weights into
superposition modes. The absorption rate of the superposition is the
weighted average of the component rates, exactly as the theorem predicts.
Non-integer ⟨n_XY⟩ values arise from Hamiltonian mixing, not from any
breakdown of the rule.

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
- **Proven:** α = 2γ⟨n_XY⟩ (algebra, exact, holds for any real H)
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
