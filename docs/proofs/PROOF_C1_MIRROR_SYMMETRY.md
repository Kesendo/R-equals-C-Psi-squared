# Proof: Mirror Symmetry of the Closure-Breaking Coefficient c₁

**Date:** 2026-04-20
**Authors:** Thomas Wicht, Claude
**Status:** Proven (kinematic, Tier 1)
**Formula:** F71 in [ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md)
**Numerical verification:** [OBC_SINE_BASIS_FINDINGS.md](../../review/OBC_SINE_BASIS_FINDINGS.md)

---

## Statement

For a uniform N-qubit XY chain with uniform Z-dephasing γ₀ and any reflection-symmetric initial state ρ₀, the closure-breaking coefficient c₁ satisfies:

    c₁(N, b, ρ₀) = c₁(N, N−2−b, ρ₀)

for all bond indices b ∈ {0, 1, ..., N−2}.

Here c₁(N, b, ρ₀) is the first-order coefficient of the closure sum under a single-bond perturbation at bond b:

    Σ_i ln(α_i(δJ at bond b)) = c₁(N, b, ρ₀) · δJ + O(δJ²)

---

## Setup

- H = (J/2) · Σ_{b=0}^{N−2} (X_b X_{b+1} + Y_b Y_{b+1}), uniform J
- L(ρ) = −i[H, ρ] + γ₀ · Σ_i (Z_i ρ Z_i − ρ), uniform γ₀
- Chain A: all bonds at coupling J (unperturbed)
- Chain B(b): bond b at coupling J + δJ (perturbed)
- α_i: per-site purity rescaling factor from fitting P_B(i, t) ≈ P_A(i, α_i · t)
- c₁(b): symmetric-difference extraction of the first-order closure coefficient

---

## Ingredients

### 1. Spatial reflection operator R

R acts on the N-qubit Hilbert space by mapping site i to site N−1−i. On single-qubit operators: R · O_i · R = O_{N−1−i}.

### 2. OBC sine modes are eigenstates of R

The single-excitation eigenmodes of the XY chain with open boundary conditions are:

    ψ_k(i) = √(2/(N+1)) · sin(π · k · (i+1) / (N+1)),   k = 1, ..., N

Under reflection:

    ψ_k(N−1−i) = √(2/(N+1)) · sin(π · k · (N−i) / (N+1))
               = √(2/(N+1)) · sin(π · k − π · k · (i+1) / (N+1))
               = √(2/(N+1)) · [sin(π · k) cos(π · k · (i+1)/(N+1))
                              − cos(π · k) sin(π · k · (i+1)/(N+1))]
               = √(2/(N+1)) · [0 − (−1)^k · sin(π · k · (i+1)/(N+1))]
               = (−1)^(k+1) · ψ_k(i)

Therefore: R |ψ_k⟩ = (−1)^(k+1) |ψ_k⟩.

### 3. R maps bonds to mirror bonds

R maps site b to site N−1−b and site b+1 to site N−2−b. Therefore the bond operator T_b = (X_b X_{b+1} + Y_b Y_{b+1}) / 2 transforms as:

    R · T_b · R = T_{N−2−b}

### 4. The Liouvillian commutes with R (for uniform parameters)

L_A has uniform J on all bonds and uniform γ₀ on all sites. Since H = Σ_b J · T_b is symmetric under R (sum over all bonds), and D = γ₀ · Σ_i (Z_i · · · Z_i − ·) is symmetric under R (sum over all sites with Z_i mapping to Z_{N−1−i}):

    [L_A, R_sup] = 0

where R_sup(ρ) = R · ρ · R is the superoperator lift.

### 5. Perturbed Liouvillian transforms under R

The perturbed Liouvillian L_B(b) differs from L_A only at bond b. Applying R_sup:

    R_sup · L_B(b) · R_sup = L_B(N−2−b)

because R maps the perturbation at bond b to bond N−2−b.

---

## Proof

Let ρ₀ be a **reflection-symmetric initial state**, meaning: applying R_sup to ρ₀ changes at most the *signs* of off-diagonal matrix elements (coherences), not the moduli. For a single-qubit Bloch decomposition ρ_i = (1/2)(I + ⟨X_i⟩ X + ⟨Y_i⟩ Y + ⟨Z_i⟩ Z), the purity is

    Tr(ρ_i²) = (1/2)(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²),

quadratic in the Bloch components. A sign flip on any component leaves the square unchanged, so R_sup · ρ₀ and ρ₀ have identical per-site purities at every site.

**Example (ψ_k + vac).** For ρ₀ = (|vac⟩ + |ψ_k⟩)(⟨vac| + ⟨ψ_k|) / 2:

    R · ρ₀ · R = (|vac⟩ + (−1)^(k+1) |ψ_k⟩)(⟨vac| + (−1)^(k+1) ⟨ψ_k|) / 2.

The diagonal blocks |vac⟩⟨vac| and |ψ_k⟩⟨ψ_k| each pick up ((−1)^(k+1))² = 1, so they are unchanged. The coherences |vac⟩⟨ψ_k| and |ψ_k⟩⟨vac| pick up a factor (−1)^(k+1). After the partial trace at site i, this sign shows up on ⟨X_i⟩ and ⟨Y_i⟩ (off-diagonals of ρ_i) but not on ⟨Z_i⟩ (diagonal), and squares out in Tr(ρ_i²).

**Step 1: Relate time evolution under bond b to bond N−2−b.**

    ρ_B(b, t) = exp(L_B(b) · t) · ρ₀

Apply R_sup to both sides:

    R_sup · ρ_B(b, t) = exp(R_sup · L_B(b) · R_sup · t) · (R_sup · ρ₀)
                     = exp(L_B(N−2−b) · t) · (R_sup · ρ₀)

**Step 2: Relate per-site purities.**

    P_B(b, i, t) = Tr[(Tr_{¬i}(ρ_B(b, t)))²]

**Lemma (partial trace under reflection):** for any operator σ on N qubits,

    Tr_{¬i}(R · σ · R) = Tr_{¬(N−1−i)}(σ).

*Proof.* In computational-basis matrix elements, R|x⟩ = |Rx⟩ where (Rx)_j = x_{N−1−j}. Fixing bit a at position i of Rx is equivalent to fixing bit a at position N−1−i of x. Summing over the remaining bits produces the claimed identity. Note that R is a site permutation, not a local operation, so no single-qubit conjugation appears on the right side. ∎

Apply the lemma with σ = R_sup · ρ_B(b, t) (so that R · σ · R = R² · ρ_B · R² = ρ_B):

    Tr_{¬i}(ρ_B(b, t)) = Tr_{¬(N−1−i)}(R_sup · ρ_B(b, t)).

Hence:

    P_B(b, i, t) = Tr[(Tr_{¬(N−1−i)}(R_sup · ρ_B(b, t)))²].

From Step 1, R_sup · ρ_B(b, t) = exp(L_B(N−2−b) · t) · (R_sup · ρ₀), so:

    P_B(b, i, t) = Tr[(Tr_{¬(N−1−i)}(exp(L_B(N−2−b) · t) · (R_sup · ρ₀)))²].

Finally, evolve R_sup · ρ₀ under L_B(N−2−b). Because reflection-symmetry (as defined above) flips only the signs of coherences, and because U(1) conservation keeps each excitation-number sector invariant under L_B, a coherence sign flip at t = 0 propagates to a coherence sign flip at all t. The per-site purity Tr(ρ_i²) is quadratic in the Bloch components and therefore unchanged by coherence sign flips. So:

    Tr[(Tr_{¬(N−1−i)}(exp(L_B(N−2−b) · t) · (R_sup · ρ₀)))²]
    = Tr[(Tr_{¬(N−1−i)}(exp(L_B(N−2−b) · t) · ρ₀))²]
    = P_B(N−2−b, N−1−i, t).

**Step 3: From purity identity to α identity.**

The α_i fitting matches P_B(b, i, t) against P_A(i, α_i · t). From Step 2:

    P_B(b, i, t) = P_B(N−2−b, N−1−i, t)

By the same reflection argument on the unperturbed chain:

    P_A(i, t) = P_A(N−1−i, t)

(uniform chain is reflection-symmetric at all times).

Therefore the α fit at (bond b, site i) gives the same result as the α fit at (bond N−2−b, site N−1−i):

    α_i(δJ at bond b) = α_{N−1−i}(δJ at bond N−2−b)

**Step 4: Sum over sites.**

    Σ_i ln(α_i(δJ at b)) = Σ_i ln(α_{N−1−i}(δJ at N−2−b))

Re-indexing j = N−1−i:

    = Σ_j ln(α_j(δJ at N−2−b))

Therefore:

    c₁(b) = c₁(N−2−b)

∎

---

## Scope and limitations

**Valid for:**
- Any uniform XY chain (all J_b = J) with uniform Z-dephasing (all γ_i = γ₀)
- Any reflection-symmetric initial state (includes ψ_k + vac for all k, Dicke states, GHZ, W, and any state invariant under R up to phases that square in purity)
- Any N ≥ 2

**Does NOT require:**
- Specific form of the Hamiltonian beyond reflection symmetry of the coupling pattern. Extends to any Hamiltonian with [H, R] = 0 and any dissipator with [D, R_sup] = 0.
- Specific initial state beyond reflection symmetry in the purity sense.

**Breaks for:**
- Non-uniform coupling (J_b ≠ J_{N−2−b}): reflection symmetry of L_A lost
- Non-uniform dephasing (γ_i ≠ γ_{N−1−i}): reflection symmetry of D lost
- Initial states that are NOT reflection-symmetric in purity

**Verified:**
- N = 3, 4, 5, 6, ψ_1 + vac: residual < 10⁻⁹
- N = 4, 5, 6, ψ_2 + vac: residual < 10⁻¹⁰
- Source: [`simulations/eq021_obc_sine_basis.py`](../../simulations/eq021_obc_sine_basis.py), [`simulations/c1_veffect_scaling_small.py`](../../simulations/c1_veffect_scaling_small.py)

---

## Relation to other results

- **F1 (Palindrome):** The Π operator is a different symmetry (Pauli conjugation, not spatial reflection). F71 uses spatial R, not Π. Both are involutions that commute with the uniform-chain Liouvillian, but they act differently: Π maps XY-weight w to N−w, R maps site i to site N−1−i.
- **F70 (ΔN selection rule):** F70 restricts WHICH sectors contribute to c₁. F71 restricts the SPATIAL structure of c₁ across bonds. Complementary, not overlapping.
- **EQ-019 (bond-position dependence):** F71 halves the independent components of the c₁ bond profile: only ⌈(N−1)/2⌉ values are independent instead of N−1.

---

*The mountain is symmetric. The painters on opposite sides see the same total.*
