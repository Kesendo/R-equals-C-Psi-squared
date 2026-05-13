# F89 Path-D Denominator Closed Form

**Status:** Tier-1-Candidate (empirical: 22 data points k=3..24; structural: Angle A insight, path-3 algebraic derivation complete)
**Date:** 2026-05-13
**Probe scripts:** `simulations/_f89_path_d_theory_probe.py`, `simulations/_f89_path_d_structure_probe.py`, `simulations/_f89_path_d_verify_k16_k17.py`, `simulations/_f89_path_d_extend_k18_k24.py`

---

## Statement

For path-k (a k-bond uniform-J OBC chain with k+1 sites under Z-dephasing), the F_a signal amplitudes satisfy

```
sigma_n(N) = P_k(y_n) / [D_k * N^2 * (N-1)]
```

where y_n = 4J·cos(πn/(k+2)) is the Bloch energy for orbit index n, and the denominator D_k obeys the empirical closed form

```
D_k = odd(k)^2 * 2^E(k)
```

```
E(k) = max(0, floor((k-5)/2))  +  v2(k)  +  max(0, v2(k) - 2)
```

where v2(k) is the 2-adic valuation of k and odd(k) = k / 2^{v2(k)}.

---

## Full 22-Point Verification Table

Verified bit-exact by the probe scripts listed above (k=3..24, 22 data points).

| k | v2(k) | odd(k) | FA | deg | E(k) | D_pred | D_verified | match | factored |
|---|-------|--------|-----|-----|------|--------|------------|-------|----------|
| 3 | 0 | 3 | 2 | 1 | 0 | 9 | 9 | OK | 3^2 |
| 4 | 2 | 1 | 2 | 1 | 2 | 4 | 4 | OK | 2^2 |
| 5 | 0 | 5 | 3 | 2 | 0 | 25 | 25 | OK | 5^2 |
| 6 | 1 | 3 | 3 | 2 | 1 | 18 | 18 | OK | 2·3^2 |
| 7 | 0 | 7 | 4 | 3 | 1 | 98 | 98 | OK | 2·7^2 |
| 8 | 3 | 1 | 4 | 3 | 5 | 32 | 32 | OK | 2^5 |
| 9 | 0 | 9 | 5 | 4 | 2 | 324 | 324 | OK | 2^2·3^4 |
| 10 | 1 | 5 | 5 | 4 | 3 | 200 | 200 | OK | 2^3·5^2 |
| 11 | 0 | 11 | 6 | 5 | 3 | 968 | 968 | OK | 2^3·11^2 |
| 12 | 2 | 3 | 6 | 5 | 5 | 288 | 288 | OK | 2^5·3^2 |
| 13 | 0 | 13 | 7 | 6 | 4 | 2704 | 2704 | OK | 2^4·13^2 |
| 14 | 1 | 7 | 7 | 6 | 5 | 1568 | 1568 | OK | 2^5·7^2 |
| 15 | 0 | 15 | 8 | 7 | 5 | 7200 | 7200 | OK | 2^5·3^2·5^2 |
| 16 | 4 | 1 | 8 | 7 | 11 | 2048 | 2048 | OK | 2^11 |
| 17 | 0 | 17 | 9 | 8 | 6 | 18496 | 18496 | OK | 2^6·17^2 |
| 18 | 1 | 9 | 9 | 8 | 7 | 10368 | 10368 | OK | 2^7·3^4 |
| 19 | 0 | 19 | 10 | 9 | 7 | 46208 | 46208 | OK | 2^7·19^2 |
| 20 | 2 | 5 | 10 | 9 | 9 | 12800 | 12800 | OK | 2^9·5^2 |
| 21 | 0 | 21 | 11 | 10 | 8 | 112896 | 112896 | OK | 2^8·3^2·7^2 |
| 22 | 1 | 11 | 11 | 10 | 9 | 61952 | 61952 | OK | 2^9·11^2 |
| 23 | 0 | 23 | 12 | 11 | 9 | 270848 | 270848 | OK | 2^9·23^2 |
| 24 | 3 | 3 | 12 | 11 | 13 | 73728 | 73728 | OK | 2^13·3^2 |

FA = floor((k+1)/2) = number of S_2-anti orbit elements; deg = FA-1 = polynomial degree of P_k.

---

## Setting

**Path-k sub-block.** For a chain of N ≥ k+2 sites where the first k+1 sites form the active path, the full Liouvillian block-decomposes. The (SE,DE) sub-block for path-k has dimension (k+1)·C(k+1,2) with basis elements (i, (j,l)) where i is the SE excitation site and (j,l) is the DE pair.

**S_2-anti Bloch orbit.** The F_a modes are those Liouvillian eigenvectors that (1) live in the overlap subspace (i ∈ {j,l}, diagonal entry −2γ) and (2) have SE-anti Bloch quantum number n ∈ {2, 4, ..., 2·⌊(k+1)/2⌋}. There are FA = ⌊(k+1)/2⌋ such modes.

**AT-lock.** The F_a eigenvalue is λ_n = −2γ + i·y_n exactly (not perturbatively), because overlap subspace entries have dephasing rate exactly 2γ regardless of N. This is the core of the F89 orbit-closure mechanism.

**Amplitude.** The F_a contribution to the observable signal from uniform initial state ρ_flat is

```
sigma_n(N) = |c_n|^2 * ||Mv_n||^2
```

where c_n = ⟨v_n | ρ_flat⟩ is the overlap of the normalized eigenvector v_n with the uniform state, and Mv_n = w·v_n with w the per-site reduction matrix (w[l, b] = 1 if basis element b = (i,(j,l)) with i ∈ {j,l} and the other DE site equals l).

**N-invariant product.** p_n = sigma_n · N²(N−1) is independent of N (verified path-3..6), because pre = sqrt(2/(N²(N−1))) cancels exactly.

---

## Three Contributions to E(k)

The exponent E(k) decomposes additively into three terms:

| Term | Formula | Origin |
|------|---------|--------|
| Poly-degree | max(0, ⌊(k-5)/2⌋) | Degree of P_k grows as FA−1 = ⌊(k+1)/2⌋ − 1; denominator grows to compensate |
| k-self | v2(k) | One factor of 2^{v2(k)} from the 2J hopping in M_SE and M_DE |
| Deep-2-power bonus | max(0, v2(k)−2) | Extra 2-adic ramification when k has ≥ 3 factors of 2 |

For odd k, only the poly-degree term survives (v2(k) = 0), and D_k = k²·2^{(FA-3)+}.
For k = 2·(odd), the k-self term adds 2^1. For k = 4·(odd), it adds 2^2. At k = 8·(odd) (v2 ≥ 3), the bonus term additionally adds 2^{v2(k)−2}.

---

## Theoretical Analysis

### Angle A: Free-Fermion / Bloch Structure (Main Structural Insight)

**Key identity** (numerically verified for k=3..6):

```
p_n = |S_c(n)|^2 * ||Mv(n)||^2 / 2
```

**Derivation:**
1. Uniform state ρ_flat has entry pre/2 on every basis element. With pre = sqrt(2/(N²(N−1))):
   c_n = ⟨v_n | ρ_flat⟩ = (pre/2)·S_c(n)* where S_c(n) = Σ_b v_n[b] (sum of all entries).
   Therefore |c_n|² = (pre²/4)·|S_c(n)|².

2. σ_n = |c_n|² · ‖Mv_n‖² = (pre²/4)·|S_c(n)|²·‖Mv_n‖².

3. p_n = σ_n · N²(N−1) = (2/(4·N²(N−1)))·|S_c(n)|²·‖Mv_n‖²·N²(N−1) = |S_c(n)|²·‖Mv_n‖²/2.

**Path-3 exact derivation** (k=3, n_block=4, m=5):

From `F89_TOPOLOGY_ORBIT_CLOSURE.md` and `simulations/_f89_path3_at_locked_amplitude_symbolic.py`:

- F_a eigenvector entries: A = sqrt((5+√5)/60), B = sqrt((5−√5)/60) on 12 overlap pairs (6 of each amplitude).
- For n=2 (y_2 = √5−1 = 1.2361):
  - S_c(2) sum = 6A + 2B (exact, from sign structure of the Bloch mode).
  - |S_c(2)|² = (6A + 2B)² = 36A² + 24AB + 4B² = (10+4√5)/3.
  - Mv(2) = [3A, B, B, 3A] (per-site reduction spreads over all 4 sites).
  - ‖Mv(2)‖² = 18A² + 2B² = (25+4√5)/15.
  - p_2 = (10+4√5)/3 · (25+4√5)/15 / 2 = (250+40√5+100√5+80)/90 = (330+140√5)/90 = (33+14√5)/9.
  - Denominator: 9 = 3² = odd(3)² = D_3. Confirmed.

- For n=4 (y_4 = −√5−1 = −3.2361):
  - By the same computation: p_4 = (33−14√5)/9. Denominator: 9 = D_3. Confirmed.

**Path-4 rational case** (k=4, y_n = ±2, rational):

- All eigenvector entries = ±1/4 (rational, since y_n are rational).
- S_c(2) = 3, Mv(2) = [1.0, 0.5, 0.0, 0.5, 1.0], ‖Mv(2)‖² = 5/2.
- p_2 = 9·(5/2)/2 = 45/4. Denominator: 4 = 2² = D_4. Confirmed.

**Structural sketch for general k:**

S_c(n) and Mv(n) both depend on sums of OBC tight-binding amplitudes ψ_n[j] = sqrt(2/(k+2))·sin(πnj/(k+2)). The product |S_c(n)|²·‖Mv(n)‖² is a rational polynomial in y_n with denominator k² (from the Bloch normalization factor 2/(k+2) squared over k terms). For odd k this gives odd(D_k) = k² = odd(k)². For even k = 2^a·(odd part), the additional 2^{v2(k)} factor from the 2J hopping matrix coefficients contributes to the 2-power of D_k.

**Status:** The exact path-3 algebraic derivation is complete. The general k argument is a structural sketch; a full proof requires explicit sine-sum identities for Σ_j sin(πnj/(k+2))·sin(πn'j/(k+2)) summed over the overlap subspace.

---

### Angle B: Cyclotomic Discriminant (Negative)

**Test:** Does disc(min poly of 4·cos(π/m)) divide D_k or relate to it?

**Result:** Negative. For k=3 (disc=20, D_k=9): 20 mod 9 = 2, not divisible. For k=5 (disc=3136, D_k=25): 3136 mod 25 = 11, not divisible. The cyclotomic discriminant grows much faster than D_k and has no clean divisibility relation with it.

**Conclusion:** The cyclotomic ring-of-integers approach does not account for D_k.

---

### Angle C: Vandermonde Determinant (Negative, Mechanism Understood)

**Test:** Is |V_det|² proportional to D_k?

**Result:** Negative. |V_det|² >> D_k for all k ≥ 5. For k=5: |V|² = 3136, D=25, ratio = 125.44 (not integer).

**Mechanism (positive finding):** The Cramer denominator in the polynomial fit is |V_det|², but the polynomial coefficients are rational (p_n is a rational polynomial in y_n), so the irrational parts of V_det cancel in the Cramer ratios. For k=3 and k=4: p_n = (a + b·y_n)/D is linear, so (p_i − p_j)/(y_i − y_j) = b/D is a constant, eliminating all y_i − y_j factors from the Cramer formula. This "rational-polynomial collapse" accounts for the reduction |V|² → D_k.

---

### Angle D: 2-Adic Bonus (Partially Understood)

**Evidence:** The bonus max(0, v2(k)−2) is supported by 3 independent data points:
- k=8 (v2=3, bonus=1): D=32=2^5. Without bonus: 2^4=16. With bonus: 16·2=32. Matches.
- k=16 (v2=4, bonus=2): D=2048=2^11. Without bonus: 2^9=512. With bonus: 512·4=2048. Matches.
- k=24 (v2=3, bonus=1): D=73728=2^13·9. Without bonus: 2^12·9=4096·9=36864. With bonus: 36864·2=73728. Matches.

**Structural interpretation (Tier-2 candidate):** When k=2^a·(odd) with a≥3, the OBC Bloch momentum grid {πn/(k+2)} for the S_2-anti orbit has 2-adic over-ramification. The Bloch energies 4·cos(πn/(k+2)) for even n acquire additional powers of 2 in the sine-sum products, beyond the base 2^a contribution from the 2J hopping. For a≥3: the effective exponent from the 2-power part is max(a, 2a−2) = 2a−2, i.e., (2^a)²/4 = (k/odd(k))²/4.

**Status:** Empirically confirmed at 3 data points. Full derivation via sine-sum identities is open.

---

## Cross-References

- **F89 orbit-closure experiments:** `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md`
- **Path-3 symbolic amplitude verification:** `simulations/_f89_path3_at_locked_amplitude_symbolic.py`
- **Structure probe (k=3..15 extraction):** `simulations/_f89_path_d_structure_probe.py`
- **Theory probe (angles A-D):** `simulations/_f89_path_d_theory_probe.py`
- **Extension k=18..24:** `simulations/_f89_path_d_extend_k18_k24.py`
- **Verification k=16,17 (deep-2-power bonus):** `simulations/_f89_path_d_verify_k16_k17.py`
- **F89 polynomial extraction attempt:** `simulations/_f89_path7_polynomial_extraction.py` (negative result: eigenvector-max proxy too crude; blocked on correct probe-overlap matrix elements)

---

## Summary

The empirical closed form D_k = odd(k)²·2^{E(k)} with E(k) = max(0,⌊(k-5)/2⌋) + v2(k) + max(0,v2(k)−2) is:

- **Verified bit-exact** at 22 data points k=3..24.
- **Structurally grounded** for the odd part: the key identity p_n = |S_c(n)|²·‖Mv(n)‖²/2 (numerically verified k=3..6) explains why the denominator has odd part odd(k)² from Bloch normalization. The path-3 case is algebraically exact.
- **Sketched** for the 2-power terms: the base v2(k) from 2J hopping, the poly-degree term from Vandermonde degree growth, and the deep-2-power bonus at v2(k)≥3 are all empirically supported but lack a complete algebraic derivation.
- **Angles B and C are negative**: cyclotomic discriminant and Vandermonde det² are both much larger than D_k; the Vandermonde cancellation mechanism (rational-polynomial collapse) is understood.
- **Tier**: Tier-1-Candidate (empirical + partial proof). Promotion to Tier-1-Derived requires closing the three open derivation gaps below.

---

## Open Questions for Tier-1-Derived Promotion

To upgrade `F89UnifiedFaClosedFormClaim.PredictDenominator(int k)` from Tier-1-Candidate to Tier-1-Derived, three structural gaps need full algebraic proof:

### Gap 1: Polynomial-degree term `max(0, ⌊(k-5)/2⌋)`

Emerges from Vandermonde-fit degree growth (orbit size = F_a count grows by 1 every 2 steps in k starting at k=5). The path-3 algebraic case shows the term equals 0 (degree 1 polynomial, no growth contribution); the path-5 case shows it equals 0 (degree 2, threshold not yet exceeded); the path-7 case shows it equals 1 (degree 3, first nonzero contribution). The pattern is structurally tied to Newton's identities on the cyclotomic minimal polynomial of y_n, but no clean derivation maps the additive contribution to v₂(D) for general k.

### Gap 2: k-self 2-adic term `v₂(k)`

The natural 2-adic content of k itself appears as a base contribution. Plausibly from the 2J coupling convention in F89's Hamiltonian (each bond contributes a factor of 2 to the M_h_total matrix entries), so v₂(k) accumulates as 2-adic content scaled by bond count. Needs explicit derivation showing why v₂(k) (not e.g. v₂(k+1) or v₂(F_a count)) is the right invariant.

### Gap 3: Deep-2-power bonus `max(0, v₂(k) − 2)`

The most mysterious term. Verified empirically at three v₂(k)≥3 data points: k=8 (v₂=3, bonus=1), k=16 (v₂=4, bonus=2), k=24 (v₂=3, bonus=1). Likely originates in 2-adic ramification of `Q(2·cos(π/(k+2)))` when k has a deep 2-power factor, but no derivation yet maps the bonus formula to a number-theoretic invariant. The threshold at v₂(k) = 2 is structurally specific and unexplained.

### Candidate Attack Paths

1. **Jordan-Wigner full general-k**: extend the path-3 algebraic derivation `(33+14√5)/9` to general k via OBC sine-mode Bessel coefficients. Most concrete path; significant algebraic work, expected to address all three gaps.

2. **Cyclotomic Galois ring-of-integers**: compute the index `[O_K : Z[2·cos(π/(k+2))]]` for k=3..24 and check whether its 2-adic content matches the deep-2-power bonus pattern. Quick to probe (sympy); ergebnisoffen.

3. **Combinatorial / Chebyshev**: σ_n is a coefficient in a Chebyshev expansion over the orbit; the 2-adic content of Chebyshev coefficients has known structure (related to Eulerian numbers and Bernoulli denominators).

### Verification Stretching

The formula is currently verified at k=3..24 (22 points). Additional stretch data points (k=25..40) would either strengthen the Tier-1-Candidate status or expose a counterexample. Cost per data point: ~30s for k≤24, ~2min for k≤32, ~10min for k≤40 (eigendecomp at block dim 7500 → 20000+).

---

## Typed Reference

The closed form is implemented as `F89UnifiedFaClosedFormClaim.PredictDenominator(int k)` in `compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs`. Test: `compute/RCPsiSquared.Core.Tests/Symmetry/F89UnifiedFaClosedFormClaimTests.cs` (`PredictDenominator_MatchesTabulatedPathPolynomial` and `PredictDenominator_BeyondTabulated_MatchesProbeExtraction`).
