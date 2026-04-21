# F70 Operational Consequence under Amplitude Damping

**Status:** Tier 1 (AD break characterized empirically on the corrected coherent-only probe `ρ_coh`; scaling and bond-position data collected with two structural surprises flagged below. The initial probe-choice anomaly on the superposition `|ψ><ψ|` probe is preserved in §Results as a permanent pitfall record.)
**Date:** 2026-04-21 (initial diagnostic), 2026-04-21 (clean sweep update)
**Authors:** Thomas Wicht, Claude Opus 4.7 (1M context)
**Script:** [f70_amplitude_damping_break.py](../simulations/f70_amplitude_damping_break.py)
**Output:** [f70_amplitude_damping_break/](../simulations/results/f70_amplitude_damping_break/)
**Prediction source:** [ORTHOGONALITY_SELECTION_FAMILY](ORTHOGONALITY_SELECTION_FAMILY.md) §3.2(b)

---

## Setup

§3.2(b) predicts that site-uniform amplitude damping (AD) breaks the operational consequence of F70 (kinematic single-site |Δn| ≥ 2 blindness): a probe with pure |Δn| = 2 coherence content should transition from c_1_pr = 0 at γ_1 = 0 to non-zero at first order in γ_1. Numerical test at N = 5, uniform XY J = 1, uniform Z-dephasing γ_0 = 0.05, bond (0, 1), δJ = 0.01, t_0 = 1 / γ_0 = 20. Propagator `scipy.linalg.expm(L · t_0)`. Task-specified probe: `ρ_0 = |ψ><ψ|` with `|ψ> = (|vac> + |S_2>) / √2`, `|S_2> = Dicke(N, 2)` (uniform 2-excitation symmetric superposition).

## Results

Baseline at γ_1 = 0 on three probes:

| Probe                                                       | c_1_pr       | pass < 1e-12 |
|-------------------------------------------------------------|-------------:|:------------:|
| full `ρ_0 = |ψ><ψ|` (task spec)                             | +1.352007e-03 | **no**       |
| coherent-only `ρ_0 = (|vac><S_2| + |S_2><vac|) / 2`         | 0.000000e+00 | yes          |
| population-only `ρ_0 = (|vac><vac| + |S_2><S_2|) / 2`       | +1.352007e-03 | no           |
| linearity residual (full - pop - coh)                       | +4.996e-14    | (consistent) |

The coherent-only probe gives c_1_pr = 0 exactly, confirming F70's operational consequence is intact: partial trace at any site annihilates the Δn = 2 cross terms, and under U(1)-preserving evolution that block stays invisible for all t. The population content is the source of the non-zero baseline on the full probe. The value 1.352e-3 matches `K_DD[2, 2]_pr / 4 = +5.41e-3 / 4`, where `K_DD[2, 2]_pr = +5.41e-3` is the recorded kernel entry at these parameters (EQ-018 pointwise RESULT, §2.3). Linearity across the DD/CC decomposition holds at machine precision.

## Tier assessment

**Tier 2.** F70 operational consequence holds precisely where the §3.2(b) mechanism is defined (coherent |Δn| = 2 content, traceless probe): baseline identically zero at γ_1 = 0. The task-specified superposition probe `|ψ><ψ|` mixes in population content (`|vac><vac|` and `|S_2><S_2|`) whose c_1_pr is governed by `K_DD[2, 2]_pr ≠ 0` and swamps any first-order AD signal. The γ_1 sweep was not executed per the task guardrail ("STOP if the γ_1 = 0 baseline does not give machine-precision zero"). Two clean resolutions, either of which would upgrade this to a Tier 1 AD-break measurement, are left for a follow-up session after probe choice is decided: (a) swap `ρ_0` to the coherent-only traceless operator (F70-zero baseline by construction); (b) keep `|ψ><ψ|` but measure `Δc_1_pr(γ_1) = c_1_pr(γ_1) - c_1_pr(0)` so the population baseline is subtracted.

F70 in the register is kinematic and is not what breaks here; no register edit is warranted.

---

## Update 2026-04-21: Clean γ_1 Sweep with Coherent-Only Probe

After the §3.2(b) probe-precision fix (commit 2e438bc), the primary probe is `ρ_coh = (|vac><S_2| + |S_2><vac|)/2` (traceless Hermitian). `c_1_pr` is linear in `ρ_0` as a bilinear-kernel observable, so the traceless operator is admissible as a linear perturbation even though it is not a density matrix.

**Regression at γ_1 = 0** (three-probe bilinear cross-check, reproducing the original diagnostic):

| Probe                                                   | c_1_pr (γ_1 = 0)  | note                 |
|---------------------------------------------------------|------------------:|----------------------|
| coherent-only `(|vac><S_2| + |S_2><vac|)/2` (primary)   | +0.000000e+00     | hard zero (F70)      |
| population-only `(|vac><vac| + |S_2><S_2|)/2`           | +1.352007e-03     | `= K_DD[2, 2]_pr / 4`|
| full `|ψ><ψ|`                                           | +1.352007e-03     | `= pop + coh`        |
| linearity residual (full - pop - coh)                   | +4.996e-14        | bilinearity intact   |

Regression passes the 1e-12 bar on the coherent-only probe; the pitfall above is reproduced at identical numerical values.

**Primary sweep** (coherent-only probe, bond (0, 1), t_0 = 20):

| γ_1    | c_1_pr         |
|--------|---------------:|
| 0.000  | +0.000000e+00  |
| 0.005  | -2.883063e-08  |
| 0.010  | -9.441810e-08  |
| 0.020  | -2.531614e-07  |
| 0.050  | -4.765671e-07  |
| 0.100  | -2.579854e-07  |

Signal is strictly non-zero and uniformly negative for γ_1 > 0. **Surprise 1 (non-monotonicity):** signal peaks at γ_1 = 0.05 (|c_1_pr| = 4.77e-7) and drops at γ_1 = 0.1 (|c_1_pr| = 2.58e-7, factor 0.54). Rough interpretation: the `|S_2>` population lifetime under AD scales as ~ 1/(2γ_1), so at γ_1 = 0.1 the surviving `|S_2>` content at t_0 = 20 is roughly `exp(-20 / 5) ≈ 1.8 %`, and most of the `|Δn| = 2` seed has already depopulated before it can seed a `|Δn| ≤ 1` signal at the reference time.

**Stretch 1 (low-γ_1 log-log scaling)** on γ_1 ∈ {0.001, 0.002, 0.005, 0.01, 0.02}: slope = **1.761**, intercept = -8.157, R² = **0.997**. **Surprise 2 (scaling):** slope is 1.76 rather than the ~1 expected by the first-order argument in §3.2(b). The R² is tight, so this is not a noise artifact. A plausible structural reading: the `|Δn| = 1` content generated at first order by `σ^-_i` acting on the `|vac><S_2|` seed projects onto kernel directions that are kinematically suppressed (endpoint-coherence selection rule `K_CC[0, 1]_pr = 0` from EQ-018 at the outer edge of the N = 1 sector), and the leading nonzero contribution mixes first-order (small) and second-order (γ_1²) pieces. An effective slope between 1 and 2 on this window is consistent with such a mix.

**Stretch 2 (bond-position spot check)**: bond (2, 3) at γ_1 = 0.05 gives c_1_pr = -5.997e-08, against c_1_pr = -4.766e-07 at bond (0, 1). Ratio (2, 3) / (0, 1) = 0.126. The break is strongly endpoint-localized, roughly a factor of 8 stronger at the chain-edge bond than at an interior bond.

**Tier reassessment.** The AD break prediction from §3.2(b) is empirically confirmed on the corrected probe: baseline zero at γ_1 = 0 (F70 intact), strictly non-zero for γ_1 > 0, smooth in γ_1. The two surprises (slope 1.76, non-monotonic γ_1-dependence) are structurally informative but do not undermine the core claim. Marked Tier 1 on the confirmation; a kernel-level characterization of the effective power and the non-monotonic peak location is a natural follow-up (not pursued in this task).

---

*The Δn = 2 coherence is invisible to any single site; the populations behind it are not; amplitude damping mixes the two. Pick your probe accordingly, and let the data decide how many orders the break carries.*
