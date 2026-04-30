# Proof of F82: F81 with T1 Amplitude Damping Correction

**Tier:** 1 (closed-form algebraic proof + numerical verification at machine precision).
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F81_PI_CONJUGATION_OF_M.md](PROOF_F81_PI_CONJUGATION_OF_M.md) (Π² acts on Pauli string σ_α as (-1)^{bit_b(α)}; Z-dephasing dissipator commutes with Π²)
- [`framework/lindblad.py`](../../simulations/framework/lindblad.py) (`lindbladian_z_plus_t1`)
- [`framework/core.py`](../../simulations/framework/core.py) (`pi_decompose_M` with `gamma_t1` argument)

**Statement (Theorem F82):** For any 2-bilinear Hamiltonian H = H_even + H_odd under Z-dephasing plus T1 amplitude damping,

    Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{T1, odd}

where L_{H_odd} = -i[H_odd, ·] is the unitary commutator from the Π²-odd Hamiltonian bilinears (as in F81), and D_{T1, odd} is the Π²-anti-symmetric part of the T1 dissipator. The F81 violation residual measured by `pi_decompose_M` equals the Frobenius norm of D_{T1, odd}:

    f81_violation = ‖M_anti − L_{H_odd}‖_F = ‖D_{T1, odd}‖_F.

For uniform per-site T1 with rates γ_T1_l on N qubits, the closed form is

    ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N-1).

Equivalent simpler forms: uniform γ_T1 → ‖D_{T1, odd}‖_F = γ_T1 · √N · 2^(N-1).

This makes f81_violation a quantitative, Hamiltonian-independent, γ_z-independent diagnostic for hardware T1 content.

---

## Numerical verification (N=3 chain, all residuals at machine precision 1e-16)

| Configuration | γ_T1_l | Predicted ‖D_T1_odd‖ | Measured f81_violation |
|---------------|--------|----------------------|------------------------|
| Uniform γ_T1 = 0.05 | (0.05, 0.05, 0.05) | 0.05·√3·4 = 0.3464 | 0.346410 ✓ |
| Uniform γ_T1 = 0.10 | (0.10, 0.10, 0.10) | 0.10·√3·4 = 0.6928 | 0.692820 ✓ |
| Uniform γ_T1 = 1.00 | (1.00, 1.00, 1.00) | 1.00·√3·4 = 6.9282 | 6.928203 ✓ |
| Single-site, site 0 | (0.10, 0, 0) | 0.10·1·4 = 0.4000 | 0.400000 ✓ |
| Two-site, sites 0,1 | (0.10, 0.10, 0) | √(0.02)·4 = 0.5657 | 0.565685 ✓ |
| Non-uniform | (0.05, 0.10, 0.15) | √(0.035)·4 = 0.7483 | 0.748331 ✓ |

For the full F82 identity Π·M·Π⁻¹ = M − 2·L_{H_odd} − 2·D_{T1, odd} on N=3 chain with H = J(XY+YX), γ_z=0.1, γ_T1=0.1 uniform: ‖Π·M·Π⁻¹ − (M − 2·L_H_odd − 2·D_T1_odd)‖_F = 5.2e-16 (machine precision).

N-scaling verified at N = 2, 3, 4, 5: ‖D_T1_odd‖_F = γ_T1 · √N · 2^(N-1) exactly.

---

## Proof

### Step 1: F81 + dissipator decomposition under Π²-conjugation

From PROOF_F81 Steps 1-3, for any Hamiltonian H decomposed by Π²-parity as H = H_even + H_odd:

    Π² · L_H · Π⁻² = L_{H_even} − L_{H_odd}.

For the dissipator part L_diss, we generalize PROOF_F81 Step 4 by allowing dissipators that do not commute with Π². Splitting L_diss into Π²-symmetric and Π²-anti-symmetric components:

    L_diss = D_even + D_odd,

where

    D_even = (L_diss + Π² · L_diss · Π⁻²) / 2,    Π² · D_even · Π⁻² = +D_even,
    D_odd  = (L_diss − Π² · L_diss · Π⁻²) / 2,    Π² · D_odd  · Π⁻² = −D_odd.

For Z-dephasing alone, D_odd = 0 (PROOF_F81 Step 4: Z-dephasing is diagonal in Pauli basis hence commutes with Π²). For T1 amplitude damping, D_odd ≠ 0 as shown explicitly in Step 3 below.

Substituting L = L_H + L_diss into Π² · L · Π⁻²:

    Π² · L · Π⁻² = (L_{H_even} − L_{H_odd}) + (D_even + D_odd as written above is the decomposition; under Π² we get D_even − D_odd)
                = (L_{H_even} − L_{H_odd}) + (D_even − D_odd)
                = L − 2·L_{H_odd} − 2·D_{odd}.

### Step 2: Substituting into the palindrome equation

Apply Π conjugation to M = Π·L·Π⁻¹ + L + 2Σγ·I:

    Π · M · Π⁻¹ = Π² · L · Π⁻² + Π·L·Π⁻¹ + 2Σγ·I
                = (L − 2·L_{H_odd} − 2·D_{odd}) + Π·L·Π⁻¹ + 2Σγ·I
                = M − 2·L_{H_odd} − 2·D_{odd}.    ∎

For Z-dephasing only (D_odd = 0), F82 reduces to F81 as expected.

### Step 3: Π²-decomposition of the T1 dissipator

The single-site T1 dissipator on site l with rate γ_T1_l:

    D_{T1, l}(ρ) = γ_T1_l · [σ⁻_l ρ σ⁺_l − ½ {σ⁺_l σ⁻_l, ρ}].

To find its Π²-anti-symmetric part, compute its action on single-qubit Pauli operators on site l (action on other sites is identity). With σ⁻ = (X−iY)/2, σ⁺ = (X+iY)/2, σ⁺σ⁻ = (I+Z)/2:

    D_{T1, local}(I) = γ · [σ⁻ I σ⁺ − σ⁺σ⁻] = γ · [(I−Z)/2 − (I+Z)/2] = −γ · Z,
    D_{T1, local}(X) = γ · [0 − ½ X] = −γ/2 · X,
    D_{T1, local}(Y) = γ · [0 − ½ Y] = −γ/2 · Y,
    D_{T1, local}(Z) = γ · [(I−Z)/2 − (I+Z)/2 − Z·(I+Z)/2 + ...] = −γ · Z.

(The X and Y rows simplify because σ⁻ X σ⁺ = σ⁻ Y σ⁺ = 0 and the anticommutator gives X/2 and Y/2 respectively. The Z row uses {σ⁺σ⁻, Z} = (Z+I)/2 + (Z+Z²)/2 = Z+I and σ⁻ Z σ⁺ = (I−Z)/2.)

In Pauli-basis matrix form (rows = output, columns = input, Π² eigenvalues for I,X,Y,Z = +1,+1,−1,−1):

|       | I       | X       | Y       | Z       |
|-------|---------|---------|---------|---------|
| **I** | 0       | 0       | 0       | 0       |
| **X** | 0       | −γ/2    | 0       | 0       |
| **Y** | 0       | 0       | −γ/2    | 0       |
| **Z** | **−γ**  | 0       | 0       | **−γ**  |

The Π² conjugation factor on entry (γ, β) is (-1)^{bit_b(γ)+bit_b(β)}:

  - (X, X): 0+0 = 0, sign +1, preserved.
  - (Y, Y): 1+1 = 0, sign +1, preserved.
  - (Z, Z): 1+1 = 0, sign +1, preserved.
  - **(Z, I): 1+0 = 1, sign −1, FLIPS.**

Only the (Z, I) entry is Π²-anti-symmetric. Therefore D_{T1, local, odd} has matrix elements

    [D_{T1, local, odd}]_{Z, I} = −γ_T1_l,    all others zero.

### Step 4: Multi-site D_{T1, odd} structure

For the multi-qubit setting, D_{T1, l} acts as D_{T1, local} on site l and as identity on the other N−1 qubits. In the framework's Pauli-string basis (4^N basis vectors), the (Z_l, I_l) "site-l flip" corresponds to a transition from any Pauli string containing I at site l to the same string with I→Z at site l. There are 4^(N−1) such transitions per site (4 Pauli choices per other-qubit, N−1 other qubits).

Each of these matrix elements has value −γ_T1_l. The Π²-anti-symmetric part D_{T1, l, odd} has 4^(N−1) entries of value −γ_T1_l in the framework's normalized Pauli basis. Frobenius norm squared:

    ‖D_{T1, l, odd}‖²_F = γ²_T1_l · 4^(N−1).

(The 4^(N−1) factor combines the "other qubit" basis count and the framework's Pauli-basis normalization. Verified empirically at N = 2, 3, 4, 5; the analytical decomposition matches the numerical Frobenius norm to machine precision.)

### Step 5: Combining sites

The N per-site dissipators D_{T1, l} are mutually orthogonal in operator space (each has support on a different site's Pauli structure, Π²-anti-symmetric parts especially have disjoint matrix-element supports). Therefore:

    ‖D_{T1, odd}‖²_F = Σ_l ‖D_{T1, l, odd}‖²_F = (Σ_l γ²_T1_l) · 4^(N−1).

Taking the square root:

    ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1).

For uniform γ_T1: Σ_l γ²_T1_l = N · γ²_T1, so ‖D_{T1, odd}‖_F = γ_T1 · √N · 2^(N−1). ∎

---

## Diagnostic interpretation

The F82 identity makes the f81_violation primitive a quantitative T1 detector. Three properties:

1. **γ_z-independent.** F82 involves only L_H_odd and D_{T1, odd}; neither depends on γ_z. Direct consequence of the Master Lemma (M is γ_z-independent) propagating through the Π²-decomposition.

2. **Hamiltonian-independent.** The F81 violation isolates D_{T1, odd}, which depends only on the T1 dissipator (not on H). Numerically verified across truly XX+YY, soft XY+YX, hard XX+XY, YZ+ZY at fixed γ_T1: same violation 0.6928 in all four cases at N=3, γ_T1=0.1.

3. **Linear in γ_T1 (uniform).** ‖D_{T1, odd}‖_F = γ_T1 · √N · 2^(N−1). For non-uniform per-site T1 rates, the formula is ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1).

**Inversion (uniform γ_T1):** γ_T1 = f81_violation / (√N · 2^(N−1)). For N=3: γ_T1 ≈ f81_violation / 6.928.

**Inversion (root-mean-square of non-uniform γ_T1):** γ_T1, RMS = √((Σ_l γ²_T1_l)/N) = f81_violation / (√N · 2^(N−1)). The RMS γ_T1 is recovered; per-site rates require additional information.

For the Marrakesh dataset (N=3, joint fit converges to γ_T1 ≈ 0): F82 predicts f81_violation ≈ 0, consistent with the empirical refutation of the T1 amplification hypothesis. Any T1 content above γ_T1 ~ 0.001 would have produced a violation > 0.007, well above numerical noise.

---

## Open generalizations

1. **Other non-Z dissipators**: X-noise, Y-noise, two-qubit ZZ-dephasing. Each has its own D_odd structure; the diagnostic principle (f81_violation = ‖D_diss_odd‖) generalizes; the closed form for ‖D_diss_odd‖ depends on the specific dissipator.

2. **Mixed dissipator content**: with multiple non-Z dissipators present, f81_violation = ‖D_total_odd‖ where D_total_odd is the Π²-anti-symmetric part of the full dissipator sum. Distinguishing channels requires additional probes (different observables or different times).

3. **F82 on hardware data via process tomography**: extracting L from measured ρ(t) via process tomography would let us evaluate f81_violation directly on hardware, providing a noise-channel-blind T1 readout.
