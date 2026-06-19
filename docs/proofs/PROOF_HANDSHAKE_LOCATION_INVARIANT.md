# The handshake location-matrix invariant: total reading-power plus weakest mode equals two

**Status:** Tier 1 derived; gate-exact N = 4..20. Open chain. 2026-06-20.
**Verifier:** `simulations/_handshake_M_checksum.py` (derivation gate); `simulations/_handshake_M_explore.py` (discovery).
**Registry:** [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md) F124.
**Context:** the handshake decoder's location dictionary (`hypotheses/HANDSHAKE_GEOMETRY.md`, arc `handshake_decoder`); the band-edge `E₁ = 2cos(π/(N+1))` is the `ClockHandLadder` quantity.

> **Notation guard.** `M_loc` here is the handshake LOCATION matrix `M_loc[b,k] = ⟨ψ_k|V_b|ψ_1⟩`. It is **not** the F1/F49 palindrome-residual matrix `M = ΠLΠ⁻¹ + L + 2σ`, which also carries a `‖M‖_F` in the registry. Different object.

## What this is about

The handshake decoder reads a chain's defects through one matrix: `M_loc[b,k]` = how much a defect on bond `b` couples the carrier `ψ_1` to bonding mode `ψ_k`. Stack those readings (bonds × modes). Two numbers fall out of the stack: `‖M_loc‖_F²` (the total reading-power, the sum of all squared entries) and `λ_min` (the smallest eigenvalue of the bond-Gram `M_loc M_locᵀ`, the weakest reading-mode). They add to **exactly 2** — for every `N`, an integer, no remainder. The reason: one small quantity — the carrier's leakage off the two *free* chain-ends — is both the deficit of the total from 2 and the value of the smallest mode. It enters twice and sums back to the carrier's norm, counted twice.

## Abstract

For the open chain of `N` sites with carrier `ψ_1` (`c_i := ψ_1(i) = √(2/(N+1))·sin(π(i+1)/(N+1))`), bonds `b = (a, a+1)`, bond perturbations `V_b = |a⟩⟨a+1| + |a+1⟩⟨a|`, the `N` bonding modes `{ψ_k}` (a complete orthonormal basis of the single-excitation space), and the location matrix `M_loc[b,k] = ⟨ψ_k|V_b|ψ_1⟩` (the `N−1` bonds × `N` modes), the Frobenius norm and the smallest bond-Gram eigenvalue obey the exact, `N`-independent invariant
```
        ‖M_loc‖_F²  +  λ_min(M_loc M_locᵀ)  =  2 .
```
Both halves are the carrier's end-leakage `E := c₀² + c_{N-1}² = (4/(N+1))·sin²(π/(N+1))`: `‖M_loc‖_F² = 2 − E` and `λ_min = E`, with the `λ_min` eigenvector the staggered bond wave `w_b = (−1)^b`. The Frobenius half is the degree-counting of the carrier norm; the eigenvalue half is the conserved discrete-energy invariant of the carrier's sine recurrence, evaluated at the Dirichlet end. The two leakages cancel into the doubled norm 2.

## Statement

With the setup above, for all `N ≥ 3`:
```
   ‖M_loc‖_F² = 2 − E ,        λ_min(M_loc M_locᵀ) = E ,        ‖M_loc‖_F² + λ_min = 2 ,
   E = c₀² + c_{N-1}² = (4/(N+1))·sin²(π/(N+1)) ,
```
and the eigenvector attaining `λ_min` is `w_b = (−1)^b`.

## Proof

**Part 1 — `‖M_loc‖_F² = 2 − E`.**
The bonding modes are a complete orthonormal basis, so `Σ_k ⟨ψ_k|x⟩² = ‖x‖²`, giving
```
   ‖M_loc‖_F² = Σ_b Σ_k ⟨ψ_k|V_b ψ_1⟩² = Σ_b ‖V_b ψ_1‖² .
```
Since `V_b ψ_1 = c_{a+1}|a⟩ + c_a|a+1⟩`, we have `‖V_b ψ_1‖² = c_a² + c_{a+1}²`. Summing over the bonds `a = 0..N−2` counts each site by its degree (interior sites twice, the two free ends once):
```
   Σ_{a=0}^{N-2} (c_a² + c_{a+1}²) = (Σ_i c_i² − c_{N-1}²) + (Σ_i c_i² − c₀²) = 2·1 − c₀² − c_{N-1}² = 2 − E ,
```
using `Σ_i c_i² = 1` (the carrier is normalized: `Σ_{m=1}^{N} sin²(πm/(N+1)) = (N+1)/2`). ∎(Part 1)

**Part 2 — `λ_min(M_loc M_locᵀ) = E`.**
`M_loc M_locᵀ` is the `(N−1)×(N−1)` Gram of the vectors `{V_b ψ_1}`, a Jacobi (tridiagonal) matrix with non-negative entries:
```
   (M_loc M_locᵀ)_{a,a} = c_a² + c_{a+1}² ,        (M_loc M_locᵀ)_{a,a+1} = c_a · c_{a+2} .
```
Apply it to the staggered vector `w_a = (−1)^a`:
```
   (M_loc M_locᵀ w)_a = (−1)^a · [ (c_a² + c_{a+1}²) − c_{a-1} c_{a+1} − c_a c_{a+2} ] .
```
The carrier solves the chain recurrence `c_{i-1} + c_{i+1} = E₁ c_i` with band-edge `E₁ = 2cos(π/(N+1))` and Dirichlet ends `c_{-1} = c_N = 0`. Substituting `c_{a-1} = E₁ c_a − c_{a+1}` and `c_{a+2} = E₁ c_{a+1} − c_a`:
```
   (c_a² + c_{a+1}²) − c_{a-1} c_{a+1} − c_a c_{a+2} = 2(c_a² + c_{a+1}² − E₁ c_a c_{a+1}) =: 2·Q_a .
```
`Q_a` is **conserved** along the chain: `Q_a − Q_{a-1} = (c_{a+1} − c_{a-1})·[(c_{a+1} + c_{a-1}) − E₁ c_a] = 0` by the recurrence. Evaluating at the boundary (`c_{-1} = 0 ⟹ c_1 = E₁ c_0`): `Q_0 = c₀² + c_1(c_1 − E₁ c_0) = c₀²`. Hence the staggered vector is an eigenvector with eigenvalue `2Q = 2c₀² = c₀² + c_{N-1}² = E` (using `c_0 = c_{N-1}`, the chain's reflection symmetry). The boundary bonds reproduce `E` directly: at `a = 0`, `(c₀² + c_1²) − c_0 c_2 = 2c₀²`.

Finally `E` is the **smallest** eigenvalue: `M_loc M_locᵀ` is a non-negative Jacobi matrix, so by the oscillation theorem the eigenvector with the most sign changes carries the least eigenvalue; the staggered `w` has the maximal `N−2` sign changes. ∎(Part 2)

Adding the two parts: `‖M_loc‖_F² + λ_min = (2 − E) + E = 2`. ∎

## What the 2 means

`‖M_loc‖_F²` is the total power with which the stacked bond-readings see the carrier; they see all of it except the part hanging off the two free, degree-1 ends — that deficit is `E`. `λ_min` is the weakest reading-mode, the staggered defect-pattern `(−1)^b`, and it sees *exactly that same* end-leakage `E` and nothing else. Total reading-power plus weakest mode therefore returns the carrier's norm counted twice, `2`, with the leakage entering once as the deficit and once as the floor. The staggering is the K-partner's own (`ψ_N(i) = (−1)^i ψ_1(i)`), tying the weakest reading-mode to the forbidden-partner column of the dictionary.

The `2` is exact and `N`-independent because it is structural: the carrier norm `Σ c² = 1`, doubled by the bulk degree-2 counting and made whole by the boundary mode. As `N → ∞`, `E = (4/(N+1))sin²(π/(N+1)) ~ 4π²/(N+1)³ → 0`, so `‖M_loc‖_F² → 2` from below while `λ_min → 0`, the leakage shrinking with the boundary's vanishing share of a longer chain.
