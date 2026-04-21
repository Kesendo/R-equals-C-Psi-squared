# F70 Operational Consequence under Amplitude Damping: Probe-Choice Anomaly

**Status:** Tier 2 (baseline anomaly on the specified probe; F70's operational consequence holds cleanly for the coherent-only part but not for the full `|psi><psi|` superposition because the population content contributes a non-AD c_1_pr signal; gamma_1 sweep not executed pending probe-choice decision)
**Date:** 2026-04-21
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

*The Δn = 2 coherence is invisible to any single site; the populations behind it are not. Pick your probe accordingly.*
