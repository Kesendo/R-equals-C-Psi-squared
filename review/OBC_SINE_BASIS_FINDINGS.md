# OBC Sine-Basis Structure of the c_1 Bond Profile

**Date:** 2026-04-20
**Authors:** Thomas Wicht, Claude
**Status:** Phase 1-4 complete. Mirror symmetry proven kinematically.
Spectral constant verdict: XY uses E_k = 2J·cos(πk/(N+1)), F2 does not match
Heisenberg SE either.
**Relates to:** [TASK_OBC_SINE_BASIS_STRUCTURE.md](../ClaudeTasks/TASK_OBC_SINE_BASIS_STRUCTURE.md),
[EQ021_FINDINGS.md](EQ021_FINDINGS.md), [RESULT_TASK_C1_VEFFECT.md](../ClaudeTasks/RESULT_TASK_C1_VEFFECT.md),
[PROOF_DELTA_N_SELECTION_RULE.md](../docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md)

## TL;DR

Four phases executed, all numerical verifications to machine precision:

1. **Phase 1 (spectral constants):** XY chain SE eigenvalues follow
   E_k = 2J·cos(πk/(N+1)), k=1..N exactly (residual 10⁻¹⁵). The w=1
   Liouvillian oscillatory modes sit at Im(λ) = |E_k|, confirming the
   mode-decomposition finding from EQ-021 Phase 1. **F2 does NOT match
   Heisenberg SE eigenvalues either**: F2 predicts N−1 values, the
   Heisenberg SE Hamiltonian has N eigenvalues, and the values
   differ. F2 must describe a different quantity (w=1 Liouvillian
   spacing or SFF peak; needs separate verification). Recommendation
   below.

2. **Phase 2 (matrix elements):** The closed form
   `⟨ψ_k|T_b|ψ_m⟩ = ψ_k(b)·ψ_m(b+1) + ψ_k(b+1)·ψ_m(b)` is exact
   (residual 10⁻¹⁷). Key geometric fact: `⟨ψ_N|T_0|ψ_1⟩ = 0` at every
   N. The k=N mode decouples entirely from the endpoint bond for
   ψ_1 input. This is Dirichlet-boundary geometry.

3. **Phase 3 (mirror symmetry):** R|ψ_k⟩ = (−1)^(k+1)·|ψ_k⟩ verified
   to 10⁻¹⁶. R·T_b·R = T_{N−2−b} exact. Combined with [L_A, R] = 0
   and the fact that per-site purity is quadratic in ρ (so phases on
   coherences square away), this proves c_1(b) = c_1(N−2−b) for ANY
   reflection-symmetric initial state, not only ψ_1+vac. This is a
   kinematic Tier-1 result, same class as F70.

4. **Phase 4 (sign pattern and ψ_2):** The bond profile is strongly
   state-dependent. For ψ_2+vac:
   - N=4: (+0.898, −0.311, +0.898) identical to ψ_1+vac
   - N=5: (+0.281, +0.357, +0.357, +0.281) all positive, interior
     exceeds endpoint
   - N=6: (+0.213, +0.740, +0.224, +0.740, +0.213) next-to-endpoint
     dominates
   
   ψ_1 has endpoints + / centre −; ψ_2 has different patterns depending
   on N. Mirror symmetry holds in both cases (the (−1)^(k+1) phase in
   R|ψ_k⟩ squares out in the purity).

## What was computed

Script: `simulations/eq021_obc_sine_basis.py`
Results: `simulations/results/eq021_obc_sine_basis/`

## Phase 1: spectral constants

### XY chain, open boundaries

The single-excitation Hamiltonian for H = (J/2)·Σ (X_b X_{b+1} + Y_b Y_{b+1})
is an N×N tridiagonal matrix with off-diagonals J and zero diagonal:

    H_SE[i, i+1] = J,  H_SE[i+1, i] = J,  H_SE[i, i] = 0

Eigenvalues (verified numerically to 10⁻¹⁵ at N=3..6):

    E_k = 2J · cos(πk / (N+1)),  k = 1, ..., N

| N | E_k (numeric)                          |
|---|----------------------------------------|
| 3 | ±1.414, 0                              |
| 4 | ±1.618, ±0.618                         |
| 5 | ±1.732, ±1.000, 0                      |
| 6 | ±1.802, ±1.247, ±0.445                 |

### w=1 Liouvillian oscillatory frequencies (XY)

Modes with Re(λ) = −2γ₀ (one-excitation coherences |vac⟩⟨ψ_k|) have:

    Im(λ) = ±|E_k| = ±2J · cos(πk / (N+1))

Exactly matches the XY SE spectrum. Each |E_k| value appears with
multiplicity consistent with the number of coherence pairs in that
sector.

### Heisenberg chain (separate check)

H_SE for Heisenberg H = J · Σ (XX+YY+ZZ) has hopping 2J plus a
diagonal ZZ shift that is NOT translation invariant (boundary sites
and bulk sites differ). Numerical eigenvalues do not follow any simple
cos formula:

| N | Heisenberg SE eigenvalues             |
|---|---------------------------------------|
| 3 | −4, 0, 2                              |
| 4 | −3.828, −1, 1.828, 3                  |
| 5 | −3.236, −1.236, 1.236, 3.236, 4       |
| 6 | −2.464, −1, 1, 3, 4.464, 5            |

### F2 comparison

F2 claims ω_k = 4J(1 − cos(πk/N)), k=1..N−1 for Heisenberg. Evaluated
numerically:

| N | F2 values                             |
|---|---------------------------------------|
| 3 | 2, 6                                  |
| 4 | 1.172, 4, 6.828                       |
| 5 | 0.764, 2.764, 5.236, 7.236            |
| 6 | 0.536, 2, 4, 6, 7.464                 |

F2 has N−1 values while Heis SE has N eigenvalues, and the values
themselves do not coincide. F2 therefore does not describe Heisenberg
SE spectrum directly.

**Verdict on F2:** F2 is valid (per the proof D10 and the SFF match
cited in ANALYTICAL_FORMULAS.md), but it describes a w=1 Liouvillian
quantity that is NOT the SE Hamiltonian eigenvalues. A plausible
interpretation is that F2 describes energy DIFFERENCES between
adjacent SE modes, or the w=1 sector of the Lindblad spectrum under
a specific Pauli-string convention. This task does not resolve
that interpretation; it only shows that the SE eigenvalues are the
quantity actually appearing in the XY mode decomposition, and they
are 2J·cos(πk/(N+1)), not F2.

## Phase 2: sine-basis matrix elements

### Closed form

In the single-excitation sector, T_b = (X_b X_{b+1} + Y_b Y_{b+1})/2
acts as the nearest-neighbour hop: T_b |i⟩ = δ_{i,b} |b+1⟩ + δ_{i,b+1} |b⟩.

Therefore in the OBC sine basis ψ_k(i) = √(2/(N+1)) · sin(πk(i+1)/(N+1)):

    ⟨ψ_k | T_b | ψ_m⟩ = ψ_k(b) · ψ_m(b+1) + ψ_k(b+1) · ψ_m(b)

verified to 10⁻¹⁷ against the direct Hamiltonian matrix elements.

### Endpoint coupling structure

`⟨ψ_k | T_0 | ψ_1⟩` for k=1..N at bond 0:

| N | k=1    | k=2    | k=3    | k=4    | k=5    | k=6    |
|---|--------|--------|--------|--------|--------|--------|
| 4 | 0.447  | 0.500  | 0.224  | 0      |        |        |
| 5 | 0.289  | 0.394  | 0.289  | 0.106  | 0      |        |
| 6 | 0.194  | 0.296  | 0.272  | 0.164  | 0.054  | 0      |

**Key observation:** ⟨ψ_N | T_0 | ψ_1⟩ = 0 exactly at every N. The
highest-frequency mode k=N never couples to the endpoint bond via
T_0 when ψ_1 is the input. This is a Dirichlet-boundary fact:
ψ_N(0)·ψ_1(1) + ψ_N(1)·ψ_1(0) vanishes because ψ_N is the staggered
version of ψ_1 (up to sign).

### Bond density ⟨ψ_1 | T_b | ψ_1⟩

| N | b=0    | b=1    | b=2    | b=3    | b=4    |
|---|--------|--------|--------|--------|--------|
| 4 | 0.447  | 0.724  | 0.447  |        |        |
| 5 | 0.289  | 0.577  | 0.577  | 0.289  |        |
| 6 | 0.194  | 0.436  | 0.543  | 0.436  | 0.194  |

Peaked in the chain centre, small at endpoints. This is the sin²
envelope of ψ_1.

### Attempt to reconstruct c_1(b) from sine-basis overlaps

Not pursued in this phase. A first-order perturbation-theory expansion
in δJ involves:

- δE_k = δJ · ⟨ψ_k | T_b | ψ_k⟩ (energy shifts)
- δ|ψ_k⟩ via Σ_{m≠k} (⟨ψ_m|T_b|ψ_k⟩ / (E_k − E_m)) · |ψ_m⟩ (eigenvector mixing)
- Feed both into the Liouvillian eigenmode expansion and the per-site
  purity kernel to get δα_i

The individual ingredients are all closed-form, but combining them
into a single analytic expression for c_1(b) requires threading
through the biorthogonal Liouvillian eigenvector expansion and the
non-linear α fit. Left as future work; the mode-decomposition data
from EQ-021 Phase 1 shows the three contributions (energy shift,
eigenvector rotation, steady-state rebalance) are comparable and
none alone reproduces c_1.

## Phase 3: mirror symmetry proof

### Ingredients (all verified)

1. **R|ψ_k⟩ = (−1)^(k+1) · |ψ_k⟩** (residual 10⁻¹⁶)
   
   Proof: sin(πk·(N−i)/(N+1)) = sin(π·k − πk(i+1)/(N+1))
   = sin(πk)cos(πk(i+1)/(N+1)) − cos(πk)sin(πk(i+1)/(N+1))
   = 0 − (−1)^k · sin(πk(i+1)/(N+1))
   = (−1)^(k+1) · sin(πk(i+1)/(N+1)).
   
   So R maps ψ_k(i) → ψ_k(N−1−i) = (−1)^(k+1) · ψ_k(i). ✓

2. **R·T_b·R = T_{N−2−b}** (residual 0, exact)
   
   R maps site b to site N−1−b, so bond (b, b+1) maps to bond
   (N−1−b, N−2−b) = (N−2−b, N−1−b) after reordering, which is bond
   N−2−b. ✓

3. **R·H_XY·R = H_XY** for uniform J (sum over bonds is symmetric).
   ✓

4. **R·D·R = D** for uniform γ₀ Z-dephasing (sum over sites is
   symmetric, and R·Z_i·R = Z_{N−1−i}). ✓

5. Combined: **[L_A, R_sup] = 0** where R_sup ρ = R·ρ·R is the
   superoperator lift. Under bond-b perturbation,
   R_sup · L_B(b) · R_sup = L_B(N−2−b). ✓

### The proof for c_1(b) = c_1(N−2−b)

Let ρ_0 be any reflection-symmetric initial state, meaning
R · ρ_0 · R differs from ρ_0 only by phases on coherences that
square out in per-site purities:
Tr[(ρ_i)²] = Tr[(R · ρ · R)_i²] where the RHS marginal is taken
at site N−1−i.

For ψ_k + vac initial states, R·|ψ_k+vac⟩ = |vac⟩ + (−1)^(k+1) · |ψ_k⟩.
For k=1 this equals the original state (R-invariant). For k=2, 4, ...
the sign flips on the coherence but squares away in ρ_i(t) purity.

Claim: under bond-b perturbation, the per-site purity at site i in
the B chain equals the per-site purity at site N−1−i in the (N−2−b)
chain:

    P_B(b, i, t) = P_B(N−2−b, N−1−i, t)

Proof: L_B(b) propagation under ρ_0 gives ρ_B(b, t) =
exp(L_B(b) t) ρ_0. Apply R_sup:
R_sup · ρ_B(b, t) = exp(R_sup · L_B(b) · R_sup · t) · R_sup · ρ_0
                  = exp(L_B(N−2−b) · t) · R_sup · ρ_0.

For the purity of the marginal at site i:
P_B(b, i, t) = Tr[(Tr_{¬i}(ρ_B(b, t)))²]
            = Tr[(Tr_{¬i}(R_sup·R_sup·ρ_B(b, t)))²] (insert R² = I)
            = Tr[(Tr_{¬(N−1−i)}(R_sup·ρ_B(b, t)))²] (reflection maps
              site i → site N−1−i in the marginal)
            = Tr[(Tr_{¬(N−1−i)}(exp(L_B(N−2−b) · t) · R_sup · ρ_0))²]
            = Tr[(Tr_{¬(N−1−i)}(ρ_B(N−2−b, t)))²] (R_sup · ρ_0 has
              same purities as ρ_0 because coherence phases square)
            = P_B(N−2−b, N−1−i, t). ✓

Therefore α_i(δJ at bond b) = α_{N−1−i}(δJ at bond N−2−b), and

    Σ_i ln(α_i(δJ at b)) = Σ_i ln(α_{N−1−i}(δJ at N−2−b))
                         = Σ_i ln(α_i(δJ at N−2−b))

so c_1(b) = c_1(N−2−b). QED.

This holds for ANY reflection-symmetric initial state
(in particular ψ_k+vac for all k) and any uniform J chain.

**Status:** kinematic Tier-1 result. Same class as F70 (ΔN selection
rule). Should be added to ANALYTICAL_FORMULAS.md as F71 or similar
(via chat).

## Phase 4: sign patterns and state dependence

### ψ_2+vac bond profile (new data)

| N | bond 0 | bond 1 | bond 2 | bond 3 | bond 4 | pattern       |
|---|--------|--------|--------|--------|--------|---------------|
| 4 | +0.898 | −0.311 | +0.898 |        |        | same as ψ_1   |
| 5 | +0.281 | +0.357 | +0.357 | +0.281 |        | all +, interior dominates |
| 6 | +0.213 | +0.740 | +0.224 | +0.740 | +0.213 | next-to-endpoint dominant |

Mirror symmetry holds (residuals 10⁻¹⁰ or better). But the
sign/magnitude pattern differs drastically from ψ_1+vac.

### Comparison ψ_1 vs ψ_2

| N | ψ_1 profile                          | ψ_2 profile                          |
|---|--------------------------------------|--------------------------------------|
| 3 | (+0.264, +0.264)                     | (+0.570, +0.570)                     |
| 4 | (+0.898, −0.311, +0.898)             | (+0.898, −0.311, +0.898)             |
| 5 | (+0.922, −0.212, −0.212, +0.922)     | (+0.281, +0.357, +0.357, +0.281)     |
| 6 | (+1.019, +0.660, −1.243, +0.660, +1.019) | (+0.213, +0.740, +0.224, +0.740, +0.213) |

At N=4 ψ_1 and ψ_2 give identical profiles, a non-trivial coincidence
that needs structural explanation (two matrix elements align). At N=5,
6 the profiles differ qualitatively: ψ_1 has sign alternation, ψ_2
has purely positive profile with interior dominance.

### Why sign patterns differ (qualitative)

The bond profile c_1(b) depends on:

- ⟨ψ_k | T_b | ψ_k⟩ (energy shift for the populated mode): this is
  ψ_k(b)·ψ_k(b+1) + h.c., largest where both sine amplitudes are
  large (bulk for ψ_1, specific offsets for ψ_2).
- ⟨ψ_m | T_b | ψ_k⟩ for m ≠ k (eigenvector mixing): these determine
  the direction of the rotation.

For ψ_1 (no nodes inside the chain), the density ψ_1(b)·ψ_1(b+1) is
strictly positive but peaked in the centre. For ψ_2 (one node at
i = (N−1)/2 for odd N, between two sites for even N), ψ_2(b) changes
sign at the node, so ψ_2(b)·ψ_2(b+1) is negative near the node and
positive away from it. This already predicts different bond profiles
at first order.

A quantitative reconstruction of c_1(b) from these overlaps was not
attempted in this task (Phase 2, last paragraph) but is the natural
next step.

### The N=4 coincidence

c_1(ψ_1+vac, b) = c_1(ψ_2+vac, b) exactly at N=4 for every b. At
N=5, 6 this fails. At N=4 the ψ_1 and ψ_2 modes have specific
amplitudes that make ⟨ψ_1|T_b|ψ_1⟩² + cross-terms equal to
⟨ψ_2|T_b|ψ_2⟩² + cross-terms for each b. This is an algebraic
identity at N=4 but not at other N. Not investigated further in this
task.

## Recommendations for ANALYTICAL_FORMULAS.md

**Do not edit directly.** Report through chat:

1. **Add a new formula entry for the XY OBC single-excitation
   spectrum:**

       E_k = 2J · cos(πk / (N+1)),   k = 1, ..., N
       (XY chain with OBC, hopping J from H = (J/2)(XX+YY))

   This is the quantity appearing in the w=1 Liouvillian
   oscillatory modes for XY, distinct from F2 (which is for
   Heisenberg and describes a different quantity).

2. **Clarify F2:** add a note that F2 refers to the Heisenberg w=1
   sector with a specific boundary/Pauli-string convention, NOT the
   single-excitation Hamiltonian eigenvalues. The relation
   ω_k = 4J(1−cos(πk/N)), k=1..N−1 gives N−1 values while the
   Heis SE Hamiltonian has N eigenvalues; the actual quantity F2
   describes should be resolved (via proof D10 review) and stated
   explicitly.

3. **Add the mirror-symmetry identity for c_1:**

       c_1(N, b, ρ_0) = c_1(N, N−2−b, ρ_0)
       for any reflection-symmetric ρ_0 on the uniform XY+γ₀ chain.

   This is a kinematic Tier-1 result, proof in Phase 3 above. Same
   rigour as F70.

## Files

- `simulations/eq021_obc_sine_basis.py` (Phases 1-4)
- `simulations/results/eq021_obc_sine_basis/phase1_spectral_comparison.json`
- `simulations/results/eq021_obc_sine_basis/phase2_matrix_elements.json`
- `simulations/results/eq021_obc_sine_basis/phase3_mirror_verification.json`
- `simulations/results/eq021_obc_sine_basis/phase4_psi2_bond_profile.json`
- `simulations/results/eq021_obc_sine_basis/summary.json`
- `simulations/results/eq021_obc_sine_basis/run.log`
