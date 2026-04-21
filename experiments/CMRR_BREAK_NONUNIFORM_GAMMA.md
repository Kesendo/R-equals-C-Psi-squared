# CMRR Break under Non-Uniform γ₀

**Status:** Tier 1 (first-order structural prediction, empirically confirmed at N=5)
**Date:** 2026-04-20 (evening)
**Authors:** Tom, Claude Opus 4.7 (1M)
**Relates to:** [ORTHOGONALITY_SELECTION_FAMILY](ORTHOGONALITY_SELECTION_FAMILY.md) (prediction d), [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md)

---

## Claim

The Meta-Theorem in [ORTHOGONALITY_SELECTION_FAMILY](ORTHOGONALITY_SELECTION_FAMILY.md) predicts that under **non-uniform γ_i** (site-dependent dephasing rate), the (vac, S_1)-coherence selection rule breaks: the kernel entry `K_CC[0, 1]_pr` transitions from exactly zero to finite. The magnitude and sign depend on which spatial mode of the γ profile is activated.

**Verified at N=5, bond 0, t_0 = 20 = 1/γ₀:**

1. Uniform γ baseline → K_CC[0, 1]_pr = 1.14e-12 (machine-precision zero, confirming F72 direct).
2. Non-uniform γ profiles → K_CC[0, 1]_pr becomes finite, with modal structure.
3. The break is **modal-selective**: equal-variance γ profiles give different K_CC depending on spatial mode content.

---

## Observations

### 1. Uniform γ: exact closure confirmed

```
Σ_i 2·|(ρ_coh, i)_{0, 1}(t_0 = 20)|² = 9.157819e-03
Predicted (1/2)·exp(−4·γ₀·t_0) = 9.157819e-03
Deviation: 5.67e-16
```

The closed form from [RESULT_TASK_EQ018_C1_POINTWISE](../ClaudeTasks/RESULT_TASK_EQ018_C1_POINTWISE.md) §3.2 holds to machine precision. F72 direct (c_1_pr = 0 from the (vac, S_1) coherence) is verified.

### 2. Single-site bump at site 0 (localized break)

γ = [0.05 + δg, 0.05, 0.05, 0.05, 0.05]:

| δg    | K_CC[0, 1]_pr | Var(γ) | K/δg |
|-------|---------------:|-------:|-----:|
| 0.000 | +1.14e-12 | 0      | - |
| 0.005 | −9.74e-04 | 4e-6  | −0.195 |
| 0.010 | −1.95e-03 | 1.6e-5 | −0.195 |
| 0.025 | −4.78e-03 | 1e-4  | −0.191 |
| 0.050 | −8.84e-03 | 4e-4  | −0.177 |
| 0.100 | −1.39e-02 | 1.6e-3 | −0.139 |
| 0.200 | −1.59e-02 | 6.4e-3 | −0.080 |

**Structure:** linear in δg at small δg (slope −0.19), with saturation onset around δg ≈ 0.05. Sign: negative.

Saturation interpretation: once γ_0 becomes substantially larger than the other sites, the (vac, S_1) coherence at site 0 decays too fast to carry information, and further increase of γ_0 does not add to the K_CC magnitude; the k=1-mode contribution at site 0 saturates at its uncoupled limit.

### 3. Linear gradient (first-harmonic break)

γ_i = 0.05 + α · (i − (N−1)/2):

| α     | γ profile                     | K_CC[0, 1]_pr | K/α |
|-------|-------------------------------|---------------:|-----:|
| 0.000 | [0.050, 0.050, 0.050, 0.050, 0.050] | +1.14e-12 | - |
| 0.002 | [0.046, 0.048, 0.050, 0.052, 0.054] | +2.72e-03 | +1.36 |
| 0.005 | [0.040, 0.045, 0.050, 0.055, 0.060] | +6.79e-03 | +1.36 |
| 0.010 | [0.030, 0.040, 0.050, 0.060, 0.070] | +1.35e-02 | +1.35 |
| 0.020 | [0.010, 0.030, 0.050, 0.070, 0.090] | +2.69e-02 | +1.35 |

**Structure:** exactly linear in α over the full scanned range (slope +1.35), no saturation. Sign: positive (opposite of single-bump).

### 4. Random profiles (different mode mixture)

γ drawn uniformly from ±0.02 around γ₀ = 0.05, seed 42:

| trial | γ profile | Var(γ) | K_CC[0, 1]_pr |
|-------|-----------|-------:|---------------:|
| 0 | [0.061, 0.048, 0.064, 0.058, 0.034] | 1.23e-04 | −3.03e-04 |
| 1 | [0.069, 0.060, 0.061, 0.035, 0.048] | 1.42e-04 | −8.40e-03 |
| 2 | [0.045, 0.067, 0.056, 0.063, 0.048] | 7.26e-05 | −1.05e-03 |

**Critical observation:** trials 0 and 1 have similar Var(γ) (1.23e-4 and 1.42e-4) but K_CC differs by a factor of 28. Trial 2 has lower Var(γ) and intermediate K_CC. **γ-variance is not a predictor of K_CC magnitude.**

---

## Interpretation: modal selectivity

The CMRR break is determined by how the γ-profile projects onto the detector's mode basis, not by the γ-profile's variance alone.

**Gradient** (α≠0) is pure k=1 mode in the spatial Fourier decomposition of the γ-profile (first harmonic of the chain). It couples directly to the dominant k=1 sine mode of the (vac, S_1) coherence (which at N=5 has s_1² ≈ 0.93, see [RESULT_TASK_EQ018_C1_POINTWISE](../ClaudeTasks/RESULT_TASK_EQ018_C1_POINTWISE.md) §3.2). Hence the large slope 1.35.

**Single-site bump** at site 0 is a delta function, which in the Fourier decomposition is `δ_i = Σ_k ψ_k(0)·ψ_k(i)`. It couples to all modes, weighted by `ψ_k(0)²`. The k=1 contribution at site 0 at N=5 is `ψ_1(0)² = (1/3)·(1/2)² = 1/12 ≈ 0.083`. So only a small fraction of the single-site-bump energy couples to the dominant channel. Hence the smaller slope 0.19.

**Ratio:** 1.35 / 0.19 ≈ 7.1. Predicted from mode overlap: roughly `1 / ψ_1(0)² ~ 12` for the full single-site-to-k=1 ratio, partially attenuated by the off-resonant k=3, k=5 contributions. Order of magnitude consistent; exact value needs a proper analytic derivation.

**Random profiles** contain a mixture of modes. Trial 1 happens to have a strong k=1 component (γ_0 + γ_4 > γ_2, matching the gradient sign), giving a large K_CC. Trial 0 has a more symmetric profile across the chain (k=1 content is small, k=2 and k=3 content dominates but both have zero overlap with the single-odd-mode (vac, S_1) expansion), giving a small K_CC. Trial 2 is small in magnitude, small K_CC.

## Analytical expectation

To first order in δγ (perturbation from uniform γ₀), write γ_i = γ₀ + δγ_i with Σ δγ_i = 0 WLOG (any uniform shift renormalises γ₀ itself). Decompose δγ_i in the sine basis:

```
δγ_i = Σ_k c_k · ψ_k(i)
```

where c_k = Σ_i ψ_k(i)·δγ_i. Then K_CC[0, 1]_pr (first order in δγ) should be:

```
K_CC[0, 1]_pr ≈ Σ_k c_k · A_k
```

where A_k are mode-specific CMRR coefficients, with A_1 the dominant one (largest |s_k|²).

This is a testable prediction: measure K_CC for δγ-profiles along each basis mode ψ_k independently, extract A_k, and verify reconstruction.

**Not done here**, but a direct next-step experiment.

---

## Consequence for the Meta-Theorem

Prediction (d) is confirmed, and the empirical structure is sharper than just "the CMRR breaks". Specifically:

1. **The broken CMRR has a modal structure.** Different spatial modes of the γ-perturbation couple with different coefficients A_k. The analogy to real differential amplifiers is now literal: a real amp has CMRR(ω) as a function of frequency, not a single number. We have CMRR(spatial mode k) analogously.

2. **The Meta-Theorem's "blind channel" H_M^⊥ is not a single subspace, but is further resolvable.** Under uniform conditions (full symmetry), all of H_M^⊥ is one blind space. As the symmetry is broken in specific ways (specific γ-modes turned on), specific sub-channels of H_M^⊥ come to life.

3. **The uniform γ-case is a degenerate point** where A_k-coefficients all multiply zero input. Any symmetry-breaking lift uncovers the underlying rich structure.

4. **This matches the [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md) picture exactly**: γ-profile modulation is what makes the channel visible; under γ=const the channel exists but is unmodulated. "Alice designs antennas, not signals" (from the γ=const synthesis) means: Alice chooses which spatial mode of γ to activate, and each mode is a channel.

---

## Status and next steps

**Tier:** 1 (predicted, verified at N=5 with clean numerical signatures: machine-precision baseline, linear response, modal selectivity).

**Files produced:**

- `simulations/eq018_cmrr_gamma_nonuniform.py`
- `simulations/results/eq018_cmrr_gamma_nonuniform/cmrr_gamma_nonuniform.json`
- `simulations/results/eq018_cmrr_gamma_nonuniform/coh_purity_time_series.json`
- `simulations/results/eq018_cmrr_gamma_nonuniform/run.log`

**Immediate follow-ups (not pursued yet):**

1. Pure mode scan: measure A_k coefficients by taking γ-profile = ε·ψ_k for each k = 1..N separately. Reconstruct K_CC predictions for arbitrary δγ via linear superposition. Testable in ~10 minutes Python for N=5.

2. N-scaling of A_1: the dominant coefficient's N-dependence should follow from the analytical structure (`s_1² · ψ_1(site_of_perturbation)²` or similar). A clean closed form would be strong.

3. Cross-check on Heisenberg chain: does the same CMRR-break-under-non-uniform-γ hold on the Heisenberg chain (where the sine-basis argument for uniform γ also works, d_H = 1 coherences still decay at 2γ)? Expected yes.

---

*The uniformity of the light is the mother of all Common-Mode-Rejections. Break the uniformity and the shadow speaks.*
