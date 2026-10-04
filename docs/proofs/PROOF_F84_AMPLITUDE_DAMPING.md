# Proof of F84: F82 Generalized to Thermal Amplitude Damping

**Tier:** 1 (closed-form algebraic proof + numerical verification at machine precision).
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F82_T1_DISSIPATOR_CORRECTION.md](PROOF_F82_T1_DISSIPATOR_CORRECTION.md) (single-channel σ⁻ analysis)
- [PROOF_F81_PI_CONJUGATION_OF_M.md](PROOF_F81_PI_CONJUGATION_OF_M.md) (Π² action on Pauli basis)
- [`framework/lindblad.py`](../../simulations/framework/lindblad.py) (`lindbladian_general`)

## Abstract

F82 closed the F81 correction for T1 amplitude damping at zero temperature (pure σ⁻ cooling, no thermal heating). On a real qubit running in a finite-temperature bath, the noise process is amplitude damping in both directions: σ⁻ cooling at one rate and σ⁺ heating at another, with the rates set by the bath occupation. The natural question after F82 is what happens to the F81-violation closed form when both directions are turned on.

The answer is the cleanest possible generalization. The correction term in F81's identity remains the antisymmetric part of the amplitude-damping dissipator. The closed form scales with the NET cooling rate per site: cooling rate minus heating rate, squared and summed. At zero temperature (heating = 0) it recovers F82 exactly. At detailed balance (cooling = heating) the correction vanishes entirely; the F81 identity holds exactly even in the presence of amplitude damping, because the cooling and heating contributions cancel in the antisymmetric piece. Between the two limits, the correction interpolates smoothly with the net rate.

A second observation falls out for free. Pure Pauli-channel dissipators (Z, X, Y dephasing on a single Pauli letter) are Π²-symmetric in the Pauli basis and contribute zero to f81_violation. Π²-antisymmetric content needs a jump whose Pauli parts mix bit_b parity: σ⁻ (cooling) and σ⁺ (heating) pair X with Y, register one by one and cancel at equal rates. So within the amplitude-damping and Pauli channels, F81 violations point to net cooling and never to phase-only or bit-flip-only noise along a Pauli axis; which other jumps register, [the F82 proof's section on other dissipators](PROOF_F82_T1_DISSIPATOR_CORRECTION.md#other-dissipators) reads off the jump.

The diagnostic upshot extends F82 to realistic thermal hardware. A measured `f81_violation` reads off the magnitude of the net population-pumping rate, not just the cooling rate. At detailed balance the diagnostic reports zero even though both cooling and heating are present; this is the right answer because both contribute equally and cancel in the antisymmetric piece. For hardware at very low effective temperature (T₁ ≫ T₁_pump), F82's closed form is recovered. The full F84 generalization is what you need when the system is not zero-temperature, and the structural reading is the same: F81-violation isolates the energy-flow-asymmetric noise content.

**Statement (Theorem F84):** For any 2-bilinear Hamiltonian H = H_even + H_odd under Z-dephasing plus thermal amplitude damping with per-site cooling rate γ_↓_l and heating rate γ_↑_l:

    Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{AmplDamp, odd}

with closed form:

    ‖D_{AmplDamp, odd}‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)
                          = |Δγ|_RMS · √N · 2^(N−1)         (uniform Δγ = γ_↓ − γ_↑)

where Δγ_l = γ_↓_l − γ_↑_l is the *net* cooling rate at site l. F82 is recovered as the special case γ_↑_l = 0 (vacuum bath / zero temperature).

**Corollary (single-Pauli channel cancellation):** Pure Pauli-channel dissipators D_{Z}, D_{X}, D_{Y} are Π²-symmetric and contribute zero to f81_violation. Π²-anti-symmetric content needs a jump whose Pauli components mix bit_b parity; among the channels treated here those are σ⁻ (cooling) and σ⁺ (heating). Consequently, within this family F81 violations on hardware are diagnostic of **population-inverting channels** (energy-emitting/absorbing transitions), not of phase-only or bit-flip-only noise along a Pauli axis.

---

## Numerical verification (N=3 chain, Z-dephasing γ_z=0.1, all residuals at machine precision)

| Configuration (γ_↓, γ_↑) | \|Δγ\| | Predicted ‖D_odd‖ | Measured f81_violation |
|--------------------------|--------|-------------------|------------------------|
| (0.10, 0.00) cooling only | 0.10 | 0.10·√3·4 = 0.6928 | 0.6928 ✓ |
| (0.00, 0.10) heating only | 0.10 | 0.10·√3·4 = 0.6928 | 0.6928 ✓ |
| (0.10, 0.10) detailed balance | 0.00 | 0.0000 | 0.0000 ✓ |
| (0.10, 0.05) net cooling | 0.05 | 0.05·√3·4 = 0.3464 | 0.3464 ✓ |
| (0.05, 0.10) net heating | 0.05 | 0.05·√3·4 = 0.3464 | 0.3464 ✓ |
| (0.20, 0.05) strong cooling | 0.15 | 0.15·√3·4 = 1.0392 | 1.0392 ✓ |

The violation is symmetric in γ_↓ ↔ γ_↑: only the net |γ_↓ − γ_↑| matters. Heating-only and cooling-only at equal rates produce identical violations.

---

## Pauli-Channel Cancellation Lemma (used in Step 3)

**Lemma:** Each single-qubit Pauli-channel dissipator D[c] with c ∈ {Z, X, Y} is fully Π²-symmetric in the Pauli basis: Π² · D[c] · Π⁻² = D[c]. Hence ‖D[c]_odd‖_F = 0 for c ∈ {Z, X, Y}.

**Proof.** For c a single-qubit Pauli operator (Hermitian with c² = I, and, the operative property, conjugation by c maps every Pauli to ± itself; a tilted axis like (X+Z)/√2 has the first two properties and fails the third):

    D[c](ρ) = γ · (c ρ c − ρ).

For each single-qubit Pauli σ_β, c σ_β c is again a single-qubit Pauli (with sign ±1 depending on whether [c, σ_β] = 0 or {c, σ_β} = 0). Hence D[c] is diagonal in the single-qubit Pauli basis:

  - D_{Z}(σ_β) = γ · (Z σ_β Z − σ_β) = 0 if σ_β commutes with Z (β ∈ {I, Z}); = −2γ σ_β if σ_β anticommutes with Z (β ∈ {X, Y}).
  - D_{X}(σ_β) = 0 for β ∈ {I, X}; = −2γ σ_β for β ∈ {Y, Z}.
  - D_{Y}(σ_β) = 0 for β ∈ {I, Y}; = −2γ σ_β for β ∈ {X, Z}.

In each case the matrix is diagonal in Pauli basis. The Π² conjugation factor (-1)^{bit_b(γ)+bit_b(β)} on a diagonal entry (γ = β) is (-1)^{2·bit_b(β)} = +1 always. Hence Π²-symmetric, ‖D[c]_odd‖ = 0. ∎

This Lemma confirms empirically what was verified by direct computation: D_{Z}, D_{X}, D_{Y} all give ‖D_odd‖_F = 0 exactly.

**Scope:** the same argument covers any Hermitian N-qubit Pauli STRING c = P (e.g. Z⊗Z, X⊗Y): P σ_β P = ±σ_β for every Pauli string σ_β, so D[P] is diagonal in the string basis, and the diagonal Π²-factor is again (+1) always. Hence every Pauli-string dissipator, and every mixture of them (Pauli-twirled noise), is fully Π²-symmetric and contributes zero to the violation. Verified numerically at N = 2 for Z⊗Z and X⊗Y (violation 0.0 exactly; pinned in `simulations/framework/tests/diagnostics/test_f84_amplitude_damping.py`, the Pauli-string-boundary tests). The boundary of the lemma is neither unitality nor the Pauli axis; [the F82 proof's section on other dissipators](PROOF_F82_T1_DISSIPATOR_CORRECTION.md#other-dissipators) reads the odd part off the jump. The tilted-axis dephasing operator c = (X+Z)/√2 is unital, not diagonal in the string basis, and carries a nonzero Π²-odd part (√2·γ at N = 1), while (Y+Z)/√2, no Pauli axis either, carries none. Correlated decay σ⁻⊗σ⁻ does contribute odd content; its closed form is in the first open-generalization bullet below.

---

## Proof of F84

### Step 1: D_AmplDamp_odd structure for combined cooling + heating

The combined amplitude-damping dissipator on site l with rates γ_↓_l and γ_↑_l:

    D_{AmplDamp, l}(ρ) = γ_↓_l · D[σ⁻_l](ρ) + γ_↑_l · D[σ⁺_l](ρ)

where σ⁻ = (X+iY)/2 (lowering, framework convention) and σ⁺ = (σ⁻)†. From PROOF_F82 Step 3, D[σ⁻]'s only Π²-anti-symmetric matrix element is the (Z, I) entry of value +γ_↓_l. By analogous computation for D[σ⁺]:

    D[σ⁺_local](I) = γ · [σ⁺ I σ⁻ − ½ {σ⁻σ⁺, I}]
                   = γ · [(I−Z)/2 − (I+Z)/2]
                   = −γ · Z.

(Using σ⁺σ⁻ = (I−Z)/2 and σ⁻σ⁺ = (I+Z)/2 in the framework's lowering convention.) The X, Y, Z rows mirror PROOF_F82 Step 3 with sign convention reversed from cooling. The single-site D[σ⁺] Pauli-basis matrix has only one Π²-anti-symmetric entry: (Z, I) = −γ_↑_l (opposite sign to cooling).

Combining cooling + heating at site l:

    D_{AmplDamp, l, odd} : (Z_l, I_l) entry = +γ_↓_l − γ_↑_l = +Δγ_l.

The Π²-anti-symmetric part depends only on the *net* cooling rate Δγ_l. ∎

### Step 2: Frobenius norm closed form

Following PROOF_F82 Step 4 with the substitution γ_l → Δγ_l:

    ‖D_{AmplDamp, l, odd}‖²_F = (Δγ_l)² · 4^(N−1)

(Per-site, 4^(N−1) "rest of qubits unchanged" entries in the framework's normalized Pauli basis.) Since D_{AmplDamp, l, odd} for different l have disjoint matrix-element supports, summing over sites:

    ‖D_{AmplDamp, odd}‖²_F = Σ_l (Δγ_l)² · 4^(N−1).

Taking the square root:

    ‖D_{AmplDamp, odd}‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1). ∎

### Step 3: F84 identity

The full Lindbladian under H + Z-dephasing + amplitude damping decomposes as L = L_H + L_Z + L_{AmplDamp}. The Π²-conjugation analysis from PROOF_F82 Step 4 generalizes:

  - L_Z (Z-dephasing) commutes with Π² (Pauli-Channel Cancellation Lemma applied to D_{Z}). Π²·L_Z·Π⁻² = L_Z, no contribution to D_odd.
  - L_{AmplDamp} = γ_↓ · L_{σ⁻} + γ_↑ · L_{σ⁺}. Both σ⁻ and σ⁺ contribute Π²-anti-symmetric parts; their combined D_odd has the (Z, I) site-l entry of value Δγ_l = γ_↓_l − γ_↑_l (Step 1).

Substituting into the palindrome equation as in PROOF_F82 Step 5:

    Π² · L · Π⁻² = L − 2 · L_{H_odd} − 2 · D_{AmplDamp, odd}.

(L_{H_odd} and D_{AmplDamp, odd} are mutually orthogonal in operator space: the Hamiltonian commutator part has different matrix-element support from the amplitude-damping channel part.) Then:

    Π · M · Π⁻¹ = Π² · L · Π⁻² + Π · L · Π⁻¹ + 2Σγ·I
                = (L − 2·L_{H_odd} − 2·D_{AmplDamp, odd}) + Π·L·Π⁻¹ + 2Σγ·I
                = M − 2·L_{H_odd} − 2·D_{AmplDamp, odd}. ∎

F82 is recovered as the special case γ_↑_l = 0 (vacuum bath, T=0): Δγ_l = γ_↓_l, ‖D‖_F = √(Σγ²_↓_l)·2^(N−1).

### Step 4: Diagnostic interpretation (thermodynamic reading)

The f81_violation primitive (in `pi_decompose_M`) computes ‖M_anti − L_{H_odd}‖_F. By F84 this equals ‖D_{AmplDamp, odd}‖_F:

    f81_violation = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1).

**Three thermodynamic regimes:**

| Regime | γ_↓ vs γ_↑ | f81_violation |
|--------|------------|---------------|
| Vacuum bath (T = 0) | γ_↑ = 0 | full F82: √(Σγ²_↓)·2^(N−1) |
| Detailed balance (T → ∞) | γ_↓ = γ_↑ | 0 (no signature) |
| Finite temperature | γ_↓ > γ_↑ > 0 | proportional to vacuum-fluctuation contribution |

For a thermal photon bath at frequency ω and temperature T:
- Mean occupation n_th = 1/(exp(ℏω/k_B T) − 1)
- γ_↓ = γ_0·(n_th + 1) (spontaneous + stimulated emission)
- γ_↑ = γ_0·n_th (stimulated absorption)
- Δγ = γ_↓ − γ_↑ = γ_0 (vacuum contribution; temperature-independent)

So at any temperature, f81_violation = γ_0·√N·2^(N−1) (uniform vacuum γ_0). The temperature only adds *symmetric* contributions (γ_↑ = γ_↓ component) which cancel under the F81 anti-symmetric projection. **f81_violation directly measures the vacuum (spontaneous emission) component of amplitude damping**, independent of temperature.

This is the deep thermodynamic content of F84: among the amplitude-damping and Pauli channels, only the *vacuum fluctuation* component (which exists even at T=0) reaches M_anti. Thermal photons (which couple γ_↓ and γ_↑ symmetrically) do not contribute. The Π²-anti-symmetric residue is a quantum-statistical fingerprint of zero-point fluctuations.

---

## Reading

F84 closes the dissipator side of the Π-decomposition picture:

| Channel | Π²-action | Contribution to f81_violation |
|---------|-----------|-------------------------------|
| Z-dephasing D_{Z} | symmetric | 0 |
| X-noise D_{X} | symmetric | 0 |
| Y-noise D_{Y} | symmetric | 0 |
| T1 cooling σ⁻ | anti-symmetric | +γ_↓_l per site |
| T1 heating σ⁺ | anti-symmetric | −γ_↑_l per site |
| Combined amplitude damping | anti-symmetric (net) | Δγ_l = γ_↓_l − γ_↑_l per site |

**Why no Pauli channel breaks F81 (but σ⁻, σ⁺ do):** Pauli-channel dissipators D[c] with c a Pauli operator (c² = I alone does not suffice; the scope note above) act diagonally in Pauli basis (each Pauli string is mapped to a scalar multiple of itself, with the multiplier being 0 or −2γ). Diagonal operators commute with Π² (which is also diagonal). Hence D_odd = 0. In contrast, σ⁻ = (X+iY)/2 and σ⁺ have Pauli components of both bit_b parities (X has bit_b 0, Y has 1), and the cross terms between those components are exactly the off-diagonal entries that mix Pauli strings of different Π²-parity ([the F82 proof, other dissipators](PROOF_F82_T1_DISSIPATOR_CORRECTION.md#other-dissipators)). Hermiticity is not what decides it: the Hermitian tilted axis (X+Z)/√2 of the scope note mixes in the same way.

**Hardware implication:** among the channels treated here, F84 says hardware F81 violations are quantitatively constrained to the *vacuum amplitude damping* component. Phase noise, bit-flip noise, and thermal photon equilibrium all give zero violation. A jump outside this family can register too, and [the F82 proof's section on other dissipators](PROOF_F82_T1_DISSIPATOR_CORRECTION.md#other-dissipators) says which (an axis with an X and a Y or Z component does, decay in the x basis does not), so the readout assumes the noise is amplitude damping plus Pauli channels. A measured nonzero f81_violation on hardware quantifies the spontaneous emission rate, which is the relevant rate for circuit decoherence at typical T ≪ ℏω/k_B (qubit transitions are typically several GHz, much larger than k_B T at mK temperatures).

**Verified:** N=3 chain at six (γ_↓, γ_↑) configurations, machine-precision residual.
**Open generalizations:**
- Two-qubit dissipators beyond single correlated jumps. ZZ-dephasing is a Pauli string and adds nothing. One correlated decay jump σ⁻⊗σ⁻ on a pair of sites (a pair-decay jump; collective emission σ⁻₁ + σ⁻₂ is a different one) has, by [the F82 proof's rule](PROOF_F82_T1_DISSIPATOR_CORRECTION.md#other-dissipators), an odd part whose entries are one-site I ↔ Z moves; counting them (at N = 2 four of size γ/2 and eight of size γ/4, the identity on the other sites) gives ‖D_odd‖²_F = (3/2)·γ²·4^(N−2), the γ²·6·4^(N−3) that [`2qubit_dissipator_exploration.py`](../../simulations/2qubit_dissipator_exploration.py) measured for it and its three siblings. That script also finds that such jumps on overlapping pairs do not add, their odd entries landing on the same transitions; a closed form for them, and the hardware calibration that would turn the number into a rate, are open.
- Higher-body dissipators: same conceptual framework; explicit closed forms case-by-case.
- Hardware extraction from process tomography: combining F84 with process-tomography-derived L gives a temperature-independent vacuum-rate readout.

---

## F1 T1-residual Pythagorean closure (2026-05-18)

The F82/F84 amplitude-damping content is exactly the Π²-antisymmetric half of the full F1 T1-residual `M = Π·L·Π⁻¹ + L + 2Σγ·I`. Under Π²-orthogonality the residual norm splits Pythagorean,

    ‖M(T1)‖²_F      = 4^(N−1) · [ 3·Σγ²_T1 + 4·(Σγ_T1)² ]
    ‖M_anti(T1)‖²_F = 4^(N−1) · Σγ²_T1                             (this proof; = ‖D_{T1, odd}‖²_F)
    ‖M_sym(T1)‖²_F  = 4^(N−1) · [ 2·Σγ²_T1 + 4·(Σγ_T1)² ]           (Π²-even complement)

with `‖M_anti‖² + ‖M_sym‖² = ‖M‖²` bit-exact. The anti side is purely local (coefficients (1, 0)); the cooperative cross-site `4·(Σγ)²` term lives entirely in the Π²-symmetric complement and is invisible to F81 / F82 / F84 diagnostics. Derivation: [`PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md`](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md) Step 7 (uses the F81 identity Π·M·Π⁻¹ = M − 2·D_{T1, odd} to collapse M_anti onto D_{T1, odd}). Typed claims: [`F1T1ResidualClosedForm.cs`](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs) (parent total) and [`F1T1ResidualPi2Decomposition.cs`](../../compute/RCPsiSquared.Core/F1/F1T1ResidualPi2Decomposition.cs) (this split).
