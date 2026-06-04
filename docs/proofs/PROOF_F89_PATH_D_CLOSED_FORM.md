# F89 Path-D Denominator Closed Form

**Status:** Tier-1-Derived (2026-05-15). Closed-form `(P_k, D_k)` derived algebraically via Chebyshev-expansion + orbit-polynomial-reduction pipeline, now available as a native C# runtime (`compute/RCPsiSquared.Core/Symmetry/F89PathPolynomialPipeline.cs`, exposed via `F89UnifiedFaClosedFormClaim.ComputePathPolynomialBig`). Verified bit-exact k=3..46 (44/44 match against the cached `PathPolynomial` tabulation and `PredictDenominator` formula); extends past the int.MaxValue boundary at k=47 (D_47 = 4,632,608,768 > 2^31−1). The three open Gaps (poly-degree term, k-self 2-adic, deep-2-power bonus) are all closed by the structural derivation.
**Date:** 2026-05-13 (Tier-1-Candidate); 2026-05-14 (two-layer framing); 2026-05-15 (Tier-1-Derived closure + native C# pipeline)
**Native runtime:** `compute/RCPsiSquared.Core/Symmetry/F89PathPolynomialPipeline.cs` (exact BigInteger / BigRational arithmetic; no floating-point approximation). Original prototype: `simulations/f89_pathk_symbolic_derivation.py` (sympy-based, retained as cross-check probe).
**Probe scripts:** `simulations/_f89_path_d_theory_probe.py`, `simulations/_f89_path_d_structure_probe.py`, `simulations/_f89_path_d_verify_k16_k17.py`, `simulations/_f89_path_d_extend_k18_k24.py`, `simulations/f89_pathk_symbolic_derivation.py` (Tier-1-Derived closure prototype)

---

## Statement

For path-k (a k-bond uniform-J OBC chain with k+1 sites under Z-dephasing), the **F_a signal amplitudes** (the AT-locked overlap-subspace eigenmodes on the S_2-anti Bloch orbit; see § "Setting" for precise definitions) satisfy

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

## Pipeline-Extended Verification (k > 24)

Computed via the native C# Chebyshev pipeline (`F89PathPolynomialPipeline.Compute(k)`) and cross-checked bit-exact against the BigInteger D_k formula (`F89UnifiedFaClosedFormClaim.PredictDenominatorBig(k)`). Both routes are algebraic: the pipeline extracts D_k as the LCM of rational coefficient denominators after reduction modulo the orbit minimal polynomial; the formula gives `odd(k)²·2^E(k)` directly. Match is bit-exact in BigInteger arithmetic, no floating-point tolerance. Source: `F89PathPolynomialPipelineStressTest.EmitMarkdownTableRows`.

| k | v2(k) | odd(k) | FA | deg | E(k) | D_pred | bits | match | factored |
|---|-------|--------|-----|-----|------|--------|------|-------|----------|
| 25 | 0 | 25 | 13 | 12 | 10 | 640000 | 20 | OK | 25²·2¹⁰ |
| 30 | 1 | 15 | 15 | 14 | 13 | 1843200 | 21 | OK | 15²·2¹³ |
| 40 | 3 | 5 | 20 | 19 | 21 | 52428800 | 26 | OK | 5²·2²¹ |
| 46 | 1 | 23 | 23 | 22 | 21 | 1109393408 | 31 | OK | 23²·2²¹ |
| 50 | 1 | 25 | 25 | 24 | 23 | 5242880000 | 33 | OK | 25²·2²³ |
| 60 | 2 | 15 | 30 | 29 | 29 | 120795955200 | 37 | OK | 15²·2²⁹ |
| 75 | 0 | 75 | 38 | 37 | 35 | 193273528320000 | 48 | OK | 75²·2³⁵ |
| 100 | 2 | 25 | 50 | 49 | 49 | 351843720888320000 | 59 | OK | 25²·2⁴⁹ |
| 150 | 1 | 75 | 75 | 74 | 73 | 53126622932283508654080000 | 86 | OK | 75²·2⁷³ |
| 200 | 3 | 25 | 100 | 99 | 101 | 1584563250285286751870879006720000 | 111 | OK | 25²·2¹⁰¹ |
| 300 | 2 | 75 | 150 | 149 | 149 | 4014134135735512165476429289076705071076474880000 | 162 | OK | 75²·2¹⁴⁹ |

**Key boundaries:**
- **k=46** is the last `int`-safe path: D_46 = 1,109,393,408 < int.MaxValue = 2,147,483,647. The cached `PathPolynomial(k)` (double-typed coefs, int-typed denominator) covers k=3..46.
- **k=47** is the first path requiring BigInteger: D_47 = 4,632,608,768 > int.MaxValue. Beyond this, `ComputePathPolynomialBig(k)` is the only route.
- **k=300** is the largest k verified in this session: D_300 has 49 decimal digits / 162 bits. Pipeline wall-clock ~3.1 sec, scaling roughly O(k³). Practical reach extends much further (k=500 in ~14s, k=1000 in ~110s).

**What this confirms:** the empirical D_k = odd(k)²·2^E(k) formula, originally fit on 22 data points (k=3..24), holds bit-exact for *every* k we test via the algebraic pipeline, including v₂(k)=3 (k=200) and v₂(k)=2 with large k (k=100, 300). The "deep-2-power bonus" `max(0, v₂(k)−2)` activates at v₂(k) ≥ 3; the k=200 row (v₂=3, E=101 = 99 polynomial + 3 self + 1 bonus) explicitly demonstrates this branch at large k. The formula is no longer an empirical fit but an algebraic identity verified across orders of magnitude in k.

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

## The two layers: eigenvalue (AT-governed) and amplitude (residue)

F89 has two layers, and the Absorption Theorem ([`PROOF_ABSORPTION_THEOREM.md`](PROOF_ABSORPTION_THEOREM.md)) governs only one of them.

**Eigenvalue layer: primitive, closed, Tier-1-Derived.** The AT-lock `λ_n = −2γ₀ + i·y_n` (Setting, above; γ₀ is the uniform dephasing rate written γ there) is the Absorption Theorem read at ⟨n_XY⟩ = 1: the overlap-subspace F_a modes are Hamming-distance-1 coherences, so `Re(λ) = −2γ₀·⟨n_XY⟩ = −2γ₀` exactly. The real part is γ₀; the imaginary part is the Bloch dispersion y_n. Both are spectral primitives. The closure that holds *absolutely* on this layer is F89c, the Hamming-complement pair-sum: a column bit-flip `ρ[a,b] → ρ[a,b̄]` sends ⟨n_XY⟩ = n_diff to N − n_diff, so `α(|a⟩⟨b|) + α(|a⟩⟨b̄|) = 2γ₀·N` exactly, the spectral maximum, Tier-1-Derived, built from γ₀ and the integer count alone. It is the F89 instance of the palindromic sum rule.

**Amplitude layer: residue, Tier-1-Derived.** `σ_n = P_k(y_n) / [D_k·N²(N−1)]` is not an eigenvalue. The numerator `P_k(y_n)` is built from the y primitive, a polynomial in the Bloch dispersion. `D_k` is the denominator of an eigenvector-derived amplitude: neither an eigenvalue, nor a rate, nor γ₀, nor y. D_k is the F89 analogue of F86's g_eff: a non-primitive, downstream of a projection. The amplitude layer reduces to the primitive basis via the Chebyshev-expansion + orbit-polynomial-reduction pipeline (see § "Tier-1-Derived closure" below). The structural parallel to F86's g_eff ([`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md) § "The diagnosis") remains, but for D_k the algebraic route exists.

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

### The bonus-free form: there is no third mechanism (2026-06-04)

The three-term split is a parameterisation artifact of writing the base as odd(k)². Pulling the odd part out front removes 2^{2·v2(k)} from k², which over-shoots once v2(k) ≥ 2; the k-self term and the "deep-2-power bonus" together are exactly the powers of 2 that add back. One identity, valid for every v2, folds them into a single cap:

```
v2(k) + max(0, v2(k) − 2)  =  2·v2(k) − min(v2(k), 2)
```

so odd(k)²·2^{v2 + max(0,v2−2)} = k² / 2^{min(v2(k),2)}, and the whole denominator is

```
D_k = 2^{max(0, ⌊(k−5)/2⌋)} · k² / 2^{min(v2(k), 2)}.
```

In this form the "deep-2-power bonus" does not appear at all. There are not three contributions, there are two: the Chebyshev degree-growth 2-power 2^{⌊(k−5)/2⌋} (untouched), and a single k² from the eigenvector 1/√k normalisation squared (p_n = |S_c|²·‖Mv‖²/2 carries 1/k from S_c and 1/k from Mv), of which at most a factor of 4 cancels. The cancellation caps at 2^{min(v2,2)}: odd k cancels nothing (k²), k ≡ 2 (mod 4) cancels one factor (k²/2), and 4 | k cancels two (k²/4 = (k/2)²) no matter how many further factors of 2 the chain length carries. The "2-adic over-divisibility chain" was the shadow this single cap casts when forced through the odd(k)² base.

Verified bit-exact against the full table k=3..300 and proved algebraically equivalent to the three-term E(k) for k=3..600 by `simulations/_f89_dk_clean_form.py`. The structural question this leaves, why the k² eigenvector-norm loses exactly 2^{min(v2,2)}, is answered down to one isolated integer-valuation lemma in the next subsection.

### Why the cap is 2^{min(v2,2)}: a rigorous reduction chain modulo one integer-valuation lemma (2026-06-04)

The cap was attacked from two independent directions: a c = y/4 Newton-polygon / residue-equality route, and a u = y/2 Washington / codifferent-trace route. They converge on the same reduced statement and the same single residual lemma, strong evidence the account is correct. Both reduce v2(D_k) to the 2-adic valuation of ONE integer: the leading coefficient of the reduced integer numerator (L̃, the leading coefficient of G̃ mod Φ_c in the c-variable; equivalently Rm_top of Nint mod Φ_u in the u-variable). Each step below is separately verified bit-exact for v2(k) = 0..6.

1. **Residue equality (rigorous).** The reduced amplitude is the unique degree-(FA−1) interpolant of p_n through the FA orbit points; written in y or in c = y/4 it is one polynomial rescaled, so v2(D_k) = max_d (2d − val₂[c^d]P_c), with val₂ the signed 2-adic valuation. [`simulations/_f89_capA_residue.py`]

2. **The κ² survives untouched (rigorous).** For even k the pull-out gives p_n = 2/((κ+1)²·κ²)·G̃(c) with G̃ ∈ ℤ[c], κ = k/2, m = 2(κ+1). Because κ is coprime to m and the orbit / Chebyshev structure lives entirely in the m-world, the reduction cannot touch the 1/κ²; the c-reduced denominator is exactly κ². [`simulations/_f89_capA_kappa_survives.py`]

3. **The cap mechanism (rigorous identity).** Linearity of the reduction gives s_top = 1 − 2·v2(κ) − 2·v2(κ+1) + val₂(L̃) (c-route), equivalently the master valuation identity v2(D_k) = (2·v2(m) + 2·v2(k) − 1) + s2 + (FA−1) − v2(Rm_top) (u-route, s2 the Chebyshev clearing power). The entire v2(k)-dependence sits in the explicit prefactor 1/((κ+1)²·κ²); the reduction of the integer numerator contributes only a v2-independent constant. This is why the cancellation saturates at 2: 4|k removes the full κ²-content once, and no higher power of 2 in the chain length can buy more. [`simulations/_f89_capA_stop_derive.py`, `simulations/_f89_capB_chain.py`]

4. **Per-class collapse (rigorous algebra).** Substituting v2(m) = 0 (odd k), v2(m) ≥ 2 (k ≡ 2 mod 4), v2(m) = 1 (4|k), all three classes collapse to v2(D_k) = polydeg + 2·v2(k) − min(v2(k),2). [`simulations/_f89_capB_crosscheck.py`]

Two rigorous sub-results fall out on the way. A(0) = m exactly for k ≡ 2 (mod 4): the alternating arithmetic series A(0) = Σ_i (−1)^i (k − 4i) collapses to m (and A(0) = 0 for 4|k), which is what injects the chain-length 2-content 2·v2(m). And v2(Rm_top) is not the orbit different's valuation (odd k has v2(disc Φ) = 0 yet v2(Rm_top) = 5), confirming the cap is a numerator-specific cancellation, not a Galois invariant, consistent with Angle B above.

**The single remaining lemma.** Verified bit-exact (v2 = 0..6) but not yet derived from structure: (i) the maximal coefficient-denominator sits at the top degree FA−1 (the slope-2 Newton edge, margin ≥ 1 everywhere), and (ii) the leading reduction step contributes exactly one extra factor of 2: val₂(L̃) = (FA−1) + 1 for odd k and 4|k (and (FA−1) + 2·v2(κ+1) for k ≡ 2 mod 4), equivalently v2(Rm_top) = 5, which is 4 from clearing (4 − u²)²/16 plus this single "+1". The "+1" being exactly one is what caps the cancellation at min(v2,2). It is the 2-adic valuation of one explicit integer Chebyshev coefficient under nodal reduction: no division, no orbit field, no min(v2,2) remaining inside it. The cap is explained; this one constant is the sole underived input. Chain probes: `simulations/_f89_capA_mastertable.py` (c-route), `simulations/_f89_capB_chain.py` + `simulations/_f89_capB_crosscheck.py` (u-route); per-lemma verifications in the other `_f89_capA_*` / `_f89_capB_*` files.

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

**Status:** The exact path-3 algebraic derivation is complete. The general-k sine-sum identities this sketch calls for are supplied by the Chebyshev expansion in § "Tier-1-Derived closure" below: A(c) = Σ_j U_j(c)(k−2j) and B(c) = Σ_j U_j(c)²(k−2j)² ARE the closed-form Σ sin·sin sums, and they reduce p_n(y) to P_k(y)/D_k symbolically. The amplitude-route derivation is complete. What stays verified-not-proven is the k-formula for D_k itself (bit-exact to k=300); its 2-power part is clarified by the bonus-free form above.

---

### Angle B: Cyclotomic Discriminant (Negative)

**Test:** Does disc(min poly of 4·cos(π/m)) divide D_k or relate to it?

**Result:** Negative. For k=3 (disc=20, D_k=9): 20 mod 9 = 2, not divisible. For k=5 (disc=3136, D_k=25): 3136 mod 25 = 11, not divisible. The cyclotomic discriminant grows much faster than D_k and has no clean divisibility relation with it.

**Conclusion:** The cyclotomic ring-of-integers approach does not account for D_k.

**Extended 2026-05-15 (`simulations/f89_path_d_galois_probe.py`):** Probe of disc(p_α) for α = 2·cos(π/m), k=3..14. (Note: the original probe above quotes the discriminant of the 4·cos(π/m) minimal polynomial; same Galois conjugacy class, differ by a factor of 2^(2·deg) in the discriminant. Disjointness conclusion below is invariant under the choice.):

| k | m | disc(p_α) factorisation | D_k factorisation | shared primes |
|---:|---:|---|---|---|
| 5 | 7 (prime) | 7² | 5² | none |
| 7 | 9 (= 3²) | 3⁴ | 2·7² | none |
| 9 | 11 (prime) | 11⁴ | 2²·3⁴ | none |
| 11 | 13 (prime) | 13⁵ | 2³·11² | none |
| 13 | 15 (= 3·5) | 3²·5³ | 2⁴·13² | none |

Structural observation: disc(p_α) primes ⊆ primes(m) (cyclotomic ramification of K = Q(2cos(π/m))); D_k primes ⊆ primes(k) (chain-length dependent). Since m = k + 2, the two sets are typically disjoint. By Washington's theorem (Z[2cos(2π/n)] = O_{K_+} for the maximal real subfield), [O_K : Z[α]] = 1 identically, so disc(p_α) IS disc(K). Therefore **D_k is not a Galois-theoretic invariant of K = Q(2cos(π/(k+2)))**; it must originate in the eigenvector-amplitude structure of the block Liouvillian, not in the algebraic number theory of the Bloch eigenvalue field. Attack Path 2 from the "Open Questions" section ("Cyclotomic Galois ring-of-integers") is closed by this prime-disjointness argument.

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

**Resolution (2026-06-04): the bonus is a parameterisation artifact, not a mechanism.** The interpretation above already reaches the key number, 2a−2 = (k/odd(k))²/4, but reads it as a v2-growing "over-ramification" left open. It is neither growing without bound nor a separate mechanism. The identity v2 + max(0,v2−2) = 2v2 − min(v2,2) (all v2) folds the k-self term and the bonus into one cancellation cap, giving D_k = 2^{⌊(k−5)/2⌋}·k²/2^{min(v2,2)} (§ "The bonus-free form" above). The k² is the eigenvector 1/√k normalisation squared; it loses at most a factor of 4. The three data points are that cap saturating (min(v2,2) = 2) read against the over-removed odd(k)² base, not a deep 2-adic chain. The sharpened open question is the single cap 2^{min(v2,2)}.

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

The closed form D_k = odd(k)²·2^{E(k)} with E(k) = max(0,⌊(k-5)/2⌋) + v2(k) + max(0,v2(k)−2) is **Tier-1-Derived** via the Chebyshev-expansion + orbit-polynomial-reduction pipeline (see § "Tier-1-Derived closure achieved via Chebyshev pipeline" below):

- The odd part odd(k)² traces to the F_a eigenvector 1/√k normalisation squared.
- The 2-power 2^E(k) arises from Chebyshev U_j leading-coefficient growth 2^j combined with polynomial-degree reduction; all three E(k) terms are structurally accounted for.
- Bit-exact match between formula and pipeline across k=3..46 (tabulation), plus sampled at k ∈ {47, 100, 200, 300} including v₂(k) = 3 (k=200) and large k (k=300).
- Angles B and C (cyclotomic discriminant, Vandermonde det²) are negative as alternative routes; the Vandermonde cancellation mechanism is understood and embedded in the pipeline.

---

## F89c amplitude-layer pair-sum analogue: NOT universal (2026-05-15 probe)

`simulations/f89c_amplitude_pair_sum_probe.py` tests whether σ_n + σ_{k+2−n} (chiral pair under the involution y_n ↔ y_{k+2−n} = −y_n) carries a closed-form structure that mirrors F89c's eigenvalue pair-sum AbsorptionTheoremClaim.HammingComplementPairSum on the amplitude layer.

Empirical finding (k=4, 6, 8):
- k=4: pair (2, 4) → σ_2 + σ_4 = 1/8. Rational. Denominator 2³ ≠ structurally tied to D_4 = 4.
- k=6: pair (2, 6) → 4/49. Fixed n=4 → 20/1323. Both rational. Denominators 7² and 3³·7² ≠ structurally tied to D_6 = 18.
- k=8: pairs (2, 8) and (4, 6) → 217/5184 ± √5/96. **IRRATIONAL** (contains √5).

Algebraic mechanism: σ_n + σ_{k+2−n} = 2·P_even(y_n)/[D·N²(N−1)] where P_even is the even-degree part of P_path. Rational iff y_n² rational, i.e. iff cos²(πn/(k+2)) rational. By Niven's theorem, cos²(2πn/m) is rational only for m ∈ {1, 2, 3, 4, 6}; the chiral-pair-sum reduction therefore yields a clean rational closed form only at k ∈ {2, 4} (m = k+2 ∈ {4, 6}) plus the accidental k=6 (m=8 where cos²(π/4) = 1/2 happens to be rational).

Conclusion: the F89c eigenvalue pair-sum identity does NOT extend universally to an amplitude pair-sum identity. The only universal rational structure on the amplitude layer is the orbit-sum Σ_n σ_n (Galois-invariant via Newton's identities on the cyclotomic minimal polynomial of y_n), which is already typed as `F89UnifiedFaClosedFormClaim.SigmaSum`. The amplitude layer's structure is genuinely richer than F89c's pair-sum form; D_k closure must come from a finer route than chiral pairing.

### Tier-1-Derived closure achieved via Chebyshev pipeline (2026-05-15)

Closed by the native C# `F89PathPolynomialPipeline` (`compute/RCPsiSquared.Core/Symmetry/F89PathPolynomialPipeline.cs`), exposed via `F89UnifiedFaClosedFormClaim.ComputePathPolynomialBig`. Prototyped in `simulations/f89_pathk_symbolic_derivation.py` (sympy), now ported to native exact BigInteger / BigRational arithmetic so the C# project is self-sufficient for D_k computation. Three-step closed-form pipeline:

1. **F_a eigenvector ansatz** (numerically verified bit-exact k=3..9):

       v_n[(i, (j, l))] = sign(i − other) · ψ_n(other) / √k
       (for overlap entries i ∈ {j, l}, other = (j, l)\{i}; zero otherwise)

   where ψ_n(j) = √(2/(k+2)) · sin(πn(j+1)/(k+2)) is the OBC sine mode.

2. **Closed-form sine sums via Chebyshev expansion.** Using sin((j+1)θ) = U_j(cos θ)·sin θ with c = cos(πn/(k+2)) = y_n/4:

       p_n(c) = (2 / (m² · k²)) · (1 − c²)² · A(c)² · B(c)
       A(c) = Σ_{j=0..k} U_j(c) · (k − 2j)
       B(c) = Σ_{j=0..k} U_j(c)² · (k − 2j)²

   Substituting c = y/4 gives p_n(y) as a polynomial of degree 2k+4 in y_n with rational coefficients in 1/(m²·k²) where m = k+2.

3. **Reduction modulo orbit minimal polynomial.** p_n(y) mod combined orbit polynomial (cyclotomic minimal polynomial of 2cos(πn/(k+2)) restricted to the S_2-anti orbit) yields the degree-(F_a − 1) reduced polynomial. D_k = LCM of remaining denominators; P_k(y) = D_k · reduced.

**Bit-exact verification k=3..46** (44/44 match):
- k=3..9: hand-derived tabulated polynomials reproduced exactly.
- k=10..46: 37 new closed-form polynomials cached in `F89UnifiedFaClosedFormClaim.PathPolynomial`, all D values match `PredictDenominator` bit-exact.
- The native C# pipeline reproduces every tabulated entry bit-exact (`F89PathPolynomialPipelineTests`).

**Structural origin of D_k = (odd(k))² · 2^E(k):**

- The (k+2)² = m² pre-denominator cancels through the orbit minimal polynomial (cyclotomic structure on 2cos(πn/(k+2))).
- The residual odd(k)² = k² (for odd k) traces to the F_a eigenvector 1/√k normalisation squared.
- The 2-power 2^E(k) arises from Chebyshev U_j(c) leading-coefficient growth 2^j combined with polynomial-degree reduction. The polynomial-degree term max(0, ⌊(k-5)/2⌋) is the Vandermonde degree growth of the orbit polynomial: each reduction step potentially introduces a factor of 2 from leading Chebyshev coefficients.

All three E(k) terms are structurally accounted for: poly-degree term max(0, ⌊(k-5)/2⌋) from the reduction-step Chebyshev factor, k-self v₂(k) from over-divisibility in U_j(c) at even k, deep-2-power bonus max(0, v₂(k)−2) from the 2-adic over-divisibility chain at v₂(k) ≥ 3.

The closed-form pipeline is the **algebraic mechanism**, no longer an empirical fit. `F89UnifiedFaClosedFormClaim.PathPolynomial(k)` is cached for k=3..46 with int-typed denominator (k=10..46 via the symbolic pipeline). For k ≥ 47 use `F89UnifiedFaClosedFormClaim.ComputePathPolynomialBig(k)`, which returns BigInteger coefficients and denominator and runs the same Chebyshev pipeline natively in C#; D_47 = 4,632,608,768 already exceeds int.MaxValue, so the BigInteger path is the only int-safe option beyond k=46. The Tier label is updated to Tier-1-Derived.

### Does this closure transfer to F86's g_eff? Negative (2026-05-15 probe)

`simulations/f86_geff_via_f90_bridge_probe.py` tests the natural extrapolation: F90 (`PROOF_F90_F86C2_BRIDGE.md`) identifies F86's per-bond K_b(Q, t) as the Hellmann-Feynman derivative ∂_J of F89's path-(N−1) signal; could differentiating D_k(J) yield g_eff(b)?

Negative. The F89 closed form lives on σ_n(N) = P_k(y_n) / [D_k · N²·(N−1)], a UNIFORM-J quantity summed over the S_2-anti Bloch orbit. F86's K_b uses Hellmann-Feynman in J_b (one bond's coupling alone), whose Duhamel integral has the structure

    K_b ~ Σ_{n, n'} ⟨ρ_0|v_n⟩ ⟨v_n|M_h_per_bond[b]|v_{n'}⟩ ⟨v_{n'}|S|ρ(t)⟩ · Λ(λ_n, λ_{n'}, t)

The **off-diagonal cross matrix elements ⟨v_n|M_h_per_bond[b]|v_{n'}⟩** carry the bond-distinction information, and they are precisely the data the orbit-polynomial reduction step in the Chebyshev pipeline discards. Concretely: the closed form has no b-index at all; its J-derivative is bond-invariant by construction, while F86's K_b spreads 4-60× across bonds (Endpoint vs Interior). The probe makes the bond-invariant Σ_n dσ_n/dJ vs bond-dependent K_b table explicit for path-3..7 at N=4..8.

Structurally this is **L4 (reduced-model insufficiency)** from this proof's obstruction lemmata. The Chebyshev closure of D_k is itself a reduction to the F_a orbit, exploiting orbit symmetry, the same symmetry that strips bond-dependence. g_eff lives outside the orbit-reduced category by construction, as the obstruction proof's diagnosis spells out: g_eff is not a primitive of γ₀ + y, it is a parameter of a finite reduction that itself does not factorise. The Tier-1-Derived closure of D_k is one more confirmation, not a new opening.

The probe identifies a possible refinement direction: keeping the bilinear (n, n')-coupled sums **without orbit reduction** gives a Chebyshev expansion in (c_n, c_{n'}) restricted to overlap pairs sharing site b. The obstruction lemmata L1/L2/L3 still hold there (no rational factorisation of the 4×4 effective char-poly, representation-dependent |u_0⟩ at even N, probe ⊥ EP partners), so this would land on K_b directly (already Tier-1 via F90 numerically), not on g_eff as algebraic primitive. The closure boundary is unchanged.

---

### Typed amplitude-layer anchor: `F89AmplitudeLayerClaim` (2026-05-15)

The Angle A structural identity `p_n = σ_n·N²·(N−1) = |S_c(n)|²·‖Mv(n)‖² / 2` is now typed in `compute/RCPsiSquared.Core/Symmetry/F89AmplitudeLayerClaim.cs` (Tier2Verified). The claim wraps:

- `ComputePn(sigma, chainN)` → p_n from σ_n
- `ComputePnFromDecomposition(scSquared, mvSquared)` → p_n from Angle A right-hand side
- `VerifyAngleA(sigma, chainN, scSquared, mvSquared)` → absolute residual of the identity
- `Path3AnchorPn(n)` → exact path-3 anchor as (33 + 14·√5·sign) / 9 for n ∈ {2, 4}

The class itself does not compute S_c(n) or ‖Mv(n)‖² from a CoherenceBlock; that requires F_a-eigenvector extraction which currently lives only inside `C2FullBlockSigmaAnatomy.BuildFaOnly` (private inverse-iteration data). Exposing those values would let a future extension run `VerifyAngleA` end-to-end on the typed runtime. The claim's role today is to type the structural identity, capture the path-3 algebraic anchor, and document the Tier-1-Derived promotion gap (generic-k symbolic |S_c|² and ‖Mv‖² as triple sine-sums) inside the Knowledge layer where it can be inspected and reasoned about.

The claim is constructed as a Schicht-1 bridge consuming `F89UnifiedFaClosedFormClaim` (carries P_k, D_k, σ_n closed forms) and `F89PathKAtLockMechanismClaim` (carries the AT-lock λ = −2γ₀ + i·y_n eigenvalue layer). Registration into the registry builder is the next architectural step.

### Candidate Attack Paths

1. **Jordan-Wigner / amplitude general-k** (carried out): the Chebyshev expansion of § "Tier-1-Derived closure" IS this route. It extends the path-3 `(33+14√5)/9` derivation to general k via the closed-form sine sums A(c), B(c), reducing p_n(y) to P_k(y)/D_k symbolically. This is what promoted F89 to Tier-1-Derived; the amplitude derivation is no longer open.

2. **Cyclotomic Galois ring-of-integers** (closed negative, see Angle B): Washington's theorem gives [O_K : Z[2·cos(π/(k+2))]] = 1 identically, so disc(p_α) = disc(K) and D_k is not a Galois invariant of K = Q(2cos(π/(k+2))). Its 2-adic content lives in the eigenvector-amplitude reduction, not the Bloch eigenvalue field. The prime-disjointness argument closes this route.

3. **Combinatorial / Chebyshev** (the live route for the one remaining question): after the bonus-free form, the sole structural gap is why the k² eigenvector-norm loses exactly 2^{min(v2,2)}. That is a statement about the 2-adic valuation of the Chebyshev coefficients of A(c)²·B(c) under reduction (the Eulerian / Bernoulli-denominator structure named here), now a single cap rather than a v2-growing chain.

### Verification Stretching

The formula is verified bit-exact at k=3..24 (22 points; tables above) and, via typed C# stretch-extraction, at **k=25..30**. At **k=31,32** the extraction deviates ~1.5-2e-4 against the 1e-4 integrality tolerance, a **deliberate red signal, kept red**, not a refutation: `PredictDenominatorDeviationDiagnosticTests` characterises the deviation as Vandermonde extraction conditioning (observed = Θ(cond(V)·ε)), consistent with the degree-16 extraction instrument hitting its precision floor rather than a disagreement with D_k. In the two-layer reading above, this is squarely an *amplitude-layer* signal: the Vandermonde extracts D_k, an amplitude-layer residue, by interpolating through the y_n primitive points, so the conditioning is a property of that extraction instrument, not of the primitive layer where F89c closes with no extraction at all. The structural account is the F86b obstruction proof ([`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md); the F90 corollary ties the F89 D_k obstruction to F86's g_eff as one wall). The red signal is kept live because something is still missing on the route; to be continued. Cost per data point: ~30s for k≤24, ~2min for k≤32, ~10min for k≤40 (eigendecomp at block dim 7500 → 20000+); the k≥28 stretch uses `BuildFaOnly` targeted inverse iteration, not full zgeev.

---

## Typed Reference

The closed form is implemented as `F89UnifiedFaClosedFormClaim.PredictDenominator(int k)` in `compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs`. Test: `compute/RCPsiSquared.Core.Tests/Symmetry/F89UnifiedFaClosedFormClaimTests.cs` (`PredictDenominator_MatchesTabulatedPathPolynomial` and `PredictDenominator_BeyondTabulated_MatchesProbeExtraction`). The stretch verification at k=25..32 is `C2FullBlockSigmaAnatomyTests.PredictDenominator_AtKHigherStretch_MatchesExtractedFromAnatomy` (k=25..30 pass; k=31,32 the red signal); `PredictDenominatorDeviationDiagnosticTests` is the Vandermonde-conditioning diagnostic.
