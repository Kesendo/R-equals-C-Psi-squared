# Proof of F82: F81 with T1 Amplitude Damping Correction

**Tier:** 1 (closed-form algebraic proof + numerical verification at machine precision).
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F81_PI_CONJUGATION_OF_M.md](PROOF_F81_PI_CONJUGATION_OF_M.md) (Π² acts on Pauli string σ_α as (-1)^{bit_b(α)}; Z-dephasing dissipator commutes with Π²)
- [`framework/lindblad.py`](../../simulations/framework/lindblad.py) (`lindbladian_z_plus_t1`)
- [`framework/diagnostics/f81_pi_decomposition.py`](../../simulations/framework/diagnostics/f81_pi_decomposition.py) (`pi_decompose_M` with `gamma_t1` argument)

## Abstract

F81 says that under pure Z-dephasing, the antisymmetric part of M (the Π-anti piece) equals exactly the Π²-odd Hamiltonian commutator. When you measure the F81 identity on a real system, any deviation from that equality is a fingerprint of something beyond pure Z-dephasing. The natural first thing to add is T1 amplitude damping, the canonical energy-loss process. The question is: what shape does the F81 violation take when T1 is present?

This proof writes the closed form. The F81 identity gets corrected by a single additional term: minus twice the Π²-antisymmetric part of the T1 dissipator. The violation residual that the framework's `pi_decompose_M` primitive reports as `f81_violation` equals exactly the Frobenius norm of that antisymmetric T1 piece. For uniform T1 across N qubits, the residual scales as a clean closed form (gamma_T1 times a known prefactor depending on N).

The mechanism is structural. T1 is implemented via the lowering operator σ⁻, which has one foot in the bit_b = 0 sector (the X component) and one foot in the bit_b = 1 sector (the Y component). The dissipator built from σ⁻ therefore carries content in both Π² sectors, and the Π²-antisymmetric piece is what survives the F81 decomposition as a non-Hamiltonian source. By contrast, pure Pauli-channel dissipators (Z, X, or Y dephasing, depolarizing) sit cleanly in the Π²-symmetric sector and contribute zero to f81_violation.

The diagnostic upshot is sharp. F81's identity is exact for pure Z-dephasing, exact for pure Pauli-channel dissipators (Z, X, Y, depolarizing per F84's generalization), and gets ONE specific correction term when T1 amplitude damping is added. A measured `f81_violation` therefore reads off the T1 amplitude-damping component of a noise process made of amplitude damping and Pauli channels (other jumps can add odd content of their own; see Other dissipators), independent of the Hamiltonian and independent of the Z-dephasing rate. F82 turns F81 into a working hardware-noise diagnostic for population-inverting channels.

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

To find its Π²-anti-symmetric part, compute its action on single-qubit Pauli operators on site l (action on other sites is identity). The framework's [`lindbladian_z_plus_t1`](../../simulations/framework/lindblad.py) uses the lowering convention σ⁻ = (X+iY)/2 = [[0, 1], [0, 0]] (taking |1⟩ → |0⟩). Then σ⁺ = (X−iY)/2 = [[0, 0], [1, 0]], σ⁻σ⁺ = (I+Z)/2 (= |0⟩⟨0|), and σ⁺σ⁻ = (I−Z)/2 (= |1⟩⟨1|):

    D_{T1, local}(I) = γ · [σ⁻ I σ⁺ − σ⁺σ⁻] = γ · [(I+Z)/2 − (I−Z)/2] = +γ · Z,
    D_{T1, local}(X) = γ · [0 − ½ X] = −γ/2 · X,
    D_{T1, local}(Y) = γ · [0 − ½ Y] = −γ/2 · Y,
    D_{T1, local}(Z) = γ · [σ⁻ Z σ⁺ − ½ {σ⁺σ⁻, Z}] = γ · [−(I+Z)/2 − (Z−I)/2] = −γ · Z.

(The X and Y rows simplify because σ⁻ X σ⁺ = σ⁻ Y σ⁺ = 0 and the anticommutator with X or Y reduces to X/2, Y/2 respectively since {Z, X} = {Z, Y} = 0. The Z row uses σ⁻ Z σ⁺ = −(I+Z)/2 (Z anticommutes with σ⁺ in the lowering decomposition) and {σ⁺σ⁻, Z} = (I−Z)/2 · Z + Z · (I−Z)/2 = Z − I (using Z² = I).)

In Pauli-basis matrix form (rows = output, columns = input, Π² eigenvalues for I,X,Y,Z = +1,+1,−1,−1):

|       | I       | X       | Y       | Z       |
|-------|---------|---------|---------|---------|
| **I** | 0       | 0       | 0       | 0       |
| **X** | 0       | −γ/2    | 0       | 0       |
| **Y** | 0       | 0       | −γ/2    | 0       |
| **Z** | **+γ**  | 0       | 0       | −γ      |

The Π² conjugation factor on entry (γ, β) is (-1)^{bit_b(γ)+bit_b(β)}:

  - (X, X): 0+0 = 0, sign +1, preserved.
  - (Y, Y): 1+1 = 0, sign +1, preserved.
  - (Z, Z): 1+1 = 0, sign +1, preserved.
  - **(Z, I): 1+0 = 1, sign −1, FLIPS.**

Only the (Z, I) entry is Π²-anti-symmetric. Therefore D_{T1, local, odd} has matrix element

    [D_{T1, local, odd}]_{Z, I} = +γ_T1_l,    all others zero.

(The sign of this entry depends on the σ⁻ convention. With σ⁻ = (X−iY)/2 (raising-into-|0⟩ convention used in some physics texts), the sign would be −γ. The closed form for ‖D_{T1, odd}‖_F derived below is convention-independent since it depends only on the magnitude γ.)

### Step 4: Multi-site D_{T1, odd} structure

For the multi-qubit setting, D_{T1, l} acts as D_{T1, local} on site l and as identity on the other N−1 qubits. In the framework's Pauli-string basis (4^N basis vectors), the (Z_l, I_l) "site-l flip" corresponds to a transition from any Pauli string containing I at site l to the same string with I→Z at site l. There are 4^(N−1) such transitions per site (4 Pauli choices per other-qubit, N−1 other qubits).

Each of these matrix elements has value ±γ_T1_l (sign per Step 3's convention note). The Π²-anti-symmetric part D_{T1, l, odd} has 4^(N−1) entries of magnitude γ_T1_l in the framework's normalized Pauli basis. Frobenius norm squared:

    ‖D_{T1, l, odd}‖²_F = γ²_T1_l · 4^(N−1).

The 4^(N−1) factor follows from the framework's Pauli-basis transform conventions: the transform M (from `_vec_to_pauli_basis_transform`) satisfies M†M = 2^N · I, so each Pauli string σ_α of Hilbert-Schmidt norm ‖σ_α‖²_HS = 2^N corresponds to a unit vector in the Pauli basis. The transform definition `palindrome_residual` line 107 reads L_pauli = M† L_vec M / 2^N, which divides matrix elements by 2^N, exactly canceling the 2^N factor coming from each Pauli string's HS norm. The net effect: a single (γ, β) Hilbert-Schmidt entry of value γ_T1_l in the operator-space representation maps to a single (γ, β) matrix element of value γ_T1_l in the framework's normalized Pauli basis. With 4^(N−1) such "rest of qubits unchanged" entries per site, ‖D_{T1, l, odd}‖²_F = γ²_T1_l · 4^(N−1). Verified empirically at N = 2, 3, 4, 5 to machine precision.

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

## Other dissipators

Step 1's identity holds for any dissipator, so what remains is to say which part of a given dissipator is Π²-odd. It can be read off the jump operator. Write a jump as a sum of Pauli strings, J = Σ_i c_i P_i. Then

    D[J](ρ) = Σ_{i,j} c_i c̄_j · (P_i ρ P_j† − ½{P_j† P_i, ρ}),

and the (i, j) term sends a string σ_α to multiples of P_i σ_α P_j†, P_j†P_i σ_α and σ_α P_j†P_i. Each of these carries bit_b(α) + bit_b(P_i) + bit_b(P_j) mod 2, and Π² is the sign (−1)^bit_b, so

    D_odd[J] = the sum of the (i, j) terms with bit_b(P_i) ≠ bit_b(P_j).

The sandwich and the anticommutator of a pair belong together: for T1, each of them alone moves I to Z and Z to I, and only their sum is Step 3's single (Z, I) entry. Four readings follow.

1. **A single Pauli string has no odd part.** Then only i = j occurs. X-, Y- and Z-noise, two-qubit ZZ-dephasing, any Pauli string, and every mixture of them (depolarizing, Pauli-twirled noise) have D_odd = 0 and leave f81_violation at exactly zero. This is what this proof's abstract says and what [the F84 proof's Pauli-string scope](PROOF_F84_AMPLITUDE_DAMPING.md) proves; the rule above is the other side of it.
2. **The odd part is where a jump mixes parities, Hermitian or not.** σ⁻ = (X + iY)/2 pairs X (bit_b = 0) with Y (bit_b = 1), and those cross terms are T1's whole odd part. A dephasing axis mixes the same way when it has an X component and a Y or Z component: c = (X + Z)/√2 is Hermitian and unital, and its X/Z cross terms are the √2·γ that the F84 proof pins at N = 1 (√2·γ·2^(N−1) on one site of N), while (Y + Z)/√2 pairs two letters of bit_b 1 and has no odd part although it is no Pauli axis either. The pair-decay jump σ⁻⊗σ⁻ = ¼(XX + iXY + iYX − YY), which the F84 proof calls correlated decay (collective emission σ⁻₁ + σ⁻₂ is a different jump), pairs {XX, YY} with {XY, YX}; every odd entry it has is a one-site I ↔ Z move, and counting them (at N = 2 four of size γ/2 and eight of size γ/4) gives ‖D_odd‖²_F = (3/2)·γ²·4^(N−2) for one such jump on one pair of sites. No list of channel names draws this line: decay in the x basis, (Z − iY)/2, has no odd part, and relaxation toward a y eigenstate, (Z + iX)/2, has one. For a single jump without an identity component, odd pairs present means odd part present: the sandwich terms of different pairs are independent superoperators, and the anticommutator terms, which carry an identity leg, cannot cancel them either.
3. **Odd parts belong to jumps and can cancel between them.** σ⁻ and σ⁺ at equal rates, or the pair {X + Y, X − Y}, carry odd parts of opposite sign that sum to zero; the F84 proof's |γ↓ − γ↑| is this cancellation.
4. **An identity component adds a commutator, not decay.** An (I, P) pair contributes i·Im(c_P c̄_I)·[P, ·], which vanishes for real coefficients.

The same rule holds for the other two turns of [the three turns](PROOF_BIT_B_PARITY_SYMMETRY.md#the-three-turns), with the count of letters that anticommute with Y or with Z in place of bit_b; that section's statement that jumps built from strings of one parity keep a turn is the case in which no odd pair exists. In [F112's](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) terms, D_odd is the dissipator's ±i content under Π-conjugation, and F112's bit_b-homogeneous collapse operators are again the jumps without an odd pair.

What the repo already held, from a sweep of the registry, the proofs, the experiments, the outbound documents, the OpenArcs registry, the caught-errors ledger, the glossary, fw.Confirmations, the Core claims, the Diagnostics witnesses and the simulation scripts. The Pauli-channel half stood in this proof's abstract, the F84 scope note, the three turns and F112's Step 2. The criterion for a single traceless jump stood in full in the F61 claim's letter rule (it breaks a parity iff its Pauli components are inhomogeneous in that bit). The tilted axis's √2·γ stood pinned in the F84 tests and in the bridge experiment, the outbound noise-asymmetry note and the OpenArcs arc on outbound label adapters, the last two also with σ⁺⊗σ⁺ writing the same odd cells as σ⁻⊗σ⁻. The pair-decay number stood in [`2qubit_dissipator_exploration.py`](../../simulations/2qubit_dissipator_exploration.py) as γ²·6·4^(N−3), for σ⁻⊗σ⁻ and its three siblings, beside detailed balance and the non-additivity over overlapping pairs; this section's count gives the same number. fw.Confirmations holds nothing on it. New here: the cross-term form of the odd part, the two pieces of T1, and the axes of the Y–Z plane as the non-Pauli axes with no odd part. Adjacent and not citing this proof: [the F113 coefficient derivation](PROOF_F113_COEFFICIENT_DERIVATION.md) rebuilds Step 3's (Z, I) entry on its own.

All of this is checked exactly, with integer and Gaussian-integer data and no tolerance, by [`f82_dissipator_odd_part_gate.py`](../../simulations/f82_dissipator_odd_part_gate.py): every single-string jump at N = 3 under each of the three turns, 60 random string sums per turn, T1's closed form at N = 2 to 4, the sandwich and anticommutator pieces of T1 separately, single jumps on both sides of the line (X + Z, X + Y and (Z + iX)/2 register; Y + Z, I + Z, (I − Z)/2 and (Z − iY)/2 do not; σ⁻ and σ⁺ register one by one and cancel as a pair), the tilted axis on one site of N = 2 and 3, pair decay σ⁻⊗σ⁻ and σ⁺⊗σ⁺ at N = 2 to 4, and F81's identity under X-noise, Y-noise, ZZ-dephasing and depolarizing.

## Open generalizations

1. **Mixed dissipator content**: with several dissipators present, f81_violation = ‖D_total_odd‖, the sum of the per-jump odd parts above. Pauli channels add nothing to it, so it sees only parity-mixing jumps; telling those apart (decay against a tilted axis, say) requires additional probes (different observables or different times).

2. **F82 on hardware data via process tomography**: extracting L from measured ρ(t) via process tomography would let us evaluate f81_violation directly on hardware, providing a T1 readout blind to every Pauli channel.
