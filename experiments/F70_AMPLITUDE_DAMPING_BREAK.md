# F70 Operational Consequence under Amplitude Damping: Clean Null for Pure Coherence Probes

**Status:** Tier 1 (clean null result: under the correct AD Lindblad, the coherent-only probe `ρ_coh = (|vac><S_2| + |S_2><vac|)/2` gives `c_1_pr = +0.000000e+00` exactly at every `γ_1 ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.1}`. Perfectly matches the analytical derivation `D_AD ρ_coh = -γ_1 ρ_coh`. This **refutes** §3.2(b)'s prediction for pure coherence probes; the prediction survives for population probes via `K_DD[2, 2]_pr`. §3.2(b) text rewrite is Chat territory. The original probe-choice pitfall from `|ψ><ψ|` is preserved in §Results below as a permanent record.)
**Date:** 2026-04-21 (initial diagnostic), 2026-04-21 (fix-and-rerun)
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

## Update 2026-04-21 (second pass): AD does NOT break operational F70 for pure coherence probes

**Context.** A first AD-break sweep on the corrected coherent-only probe was performed (reverted commits `d21eaf4` + `e226e0a`) and reported a non-zero signal (slope 1.76, bond ratio 0.126). A post-hoc bug check by Chat found that `build_liouvillian` was using `np.kron(SMi.T, SMi)`, which computes `σ^+ ⊗ σ^-` (the vec form of `σ^- ρ σ^-`, not `σ^- ρ σ^+`) in column-stacking convention. That operator is not trace-preserving (`Tr(ρ(t)) = 0.368 ≠ 1` on the single-qubit `|1><1|` test at `γ_1 = 0.5`, `t = 2`) and the entire reverted sweep was a deterministic artifact of the broken Lindblad. The bug is now fixed (commit `9c5c371`), a `lindblad_sanity_check` regression is run at every `main()` entry (trace and Hermiticity preserved to 1e-10 on vac / S_1 / S_2 test states), and the sweep has been rerun against the corrected Liouvillian.

**Result on the fixed sweep.**

| γ_1    | c_1_pr (coherent-only, bond (0, 1)) |
|--------|------------------------------------:|
| 0.000  | +0.000000e+00                       |
| 0.005  | +0.000000e+00                       |
| 0.010  | +0.000000e+00                       |
| 0.020  | +0.000000e+00                       |
| 0.050  | +0.000000e+00                       |
| 0.100  | +0.000000e+00                       |

Every entry is a **hard** zero, not float noise at the 1e-17 floor. Low-γ_1 log-log fit therefore has no positive data and no power law to report (slope reported as NaN in the result JSON). Bond (2, 3) at γ_1 = 0.05 is also `+0.000000e+00`; no ratio to bond (0, 1). Lindblad sanity at this setup: `max |Tr(L ρ)| = 3.47e-18`, `max Hermiticity residual = 2.78e-17`, both 8 orders below the 1e-10 bar.

**Analytical derivation that matches the numerics.** Starting from `ρ_coh = (|vac><S_2| + |S_2><vac|)/2`, the site-uniform AD dissipator acts as:

- Term 1, `σ^-_i ρ_coh σ^+_i`: `σ^-_i |vac> = 0` (no excitation at `i` to lower) annihilates the `|vac><S_2|` half on the ket side; `<vac| σ^+_i = 0` (Hermitian conjugate) annihilates the `|S_2><vac|` half on the bra side. Term 1 = 0.
- Term 2 + Term 3, `-½ (n_i ρ_coh + ρ_coh n_i)`: `n_i |vac> = 0` kills the bra-side contribution from the `|vac><S_2|` half; on the ket-side contribution `|S_2><vac| n_i = 0` similarly. The surviving pieces `|vac><S_2| n_i` and `n_i |S_2><vac|` involve `<S_2| n_i` and `n_i |S_2>`. Summing over `i` with `Σ_i n_i |S_2> = 2|S_2>` (each basis state in the `|S_2>` superposition has exactly 2 excited sites), the total gives `-γ_1 ρ_coh`.

So `D_AD ρ_coh = -γ_1 ρ_coh` exactly: pure exponential decay of the whole (vac, 2-exc) block, no sector-changing transfer. Under the XY `H` (U(1)-preserving) and uniform Z-dephasing, `ρ_coh(t)` mixes within the (vac, 2-exc) block but never leaves it. F70: single-site partial trace of any (vac, 2-exc) block is exactly zero, so `ρ_i(t) = 0` for every `i, t, δJ`, giving `c_1_pr = 0` identically. The hard zero in the numerics is not cancellation; it is the kinematic zero from F70 acting on a block that truly stays `|Δn| = 2`.

**What §3.2(b) got right.** F70 itself (kinematic, single-site `|Δn| ≥ 2` blindness) is untouched. Population probes do show a break: `c_1_pr(|S_2><S_2|) = K_DD[2, 2]_pr = +5.41e-3 ≠ 0` at U(1), so the mixed `ρ_pop = (|vac><vac| + |S_2><S_2|)/2` and the full `|ψ><ψ|` give the +1.352e-3 baseline documented in §Results above. Under AD the population content is where sector transfer actually happens (`σ^-_i |S_2><S_2| σ^+_i` is non-zero and drops content into `(|S_1><S_1|)` via `(N - 1) = 4` paths per site), which would produce a real AD-break signal on population probes.

**What §3.2(b) got wrong.** The heuristic "AD leaks `|Δn| ≤ 1` content from `|Δn| = 2` seed, which becomes site-local visible" is misleading for the pure coherence probe. The asymmetry between bra and ket in `|vac><S_2|` (zero on the vacuum side) combined with the Lindblad structure collapses AD to pure in-block decay. No leakage out of the `|Δn| = 2` block happens on this specific operator, so F70 is not operationally broken. The §3.2(b) text needs a Chat-side refinement distinguishing population probes (break present) from coherence probes (break absent under AD).

**Tier reassessment.** Marked Tier 1 on the clean null result: the numerical answer is unambiguous (hard zero at machine precision across the full sweep and stretches), the analytical derivation is exact and matches, and the physics story separating population-content from coherence-content AD responses is sharp. `K_DD[2, 2]_pr`-channel break for population probes is a natural follow-up (separate task) if Chat wants a positive AD-break experiment.

---

*The Δn = 2 coherence is invisible to any single site under F70; under the correct AD Lindblad, it also stays Δn = 2 throughout its evolution. No leak, no break. Sector transfer under AD is a property of the population content, not the coherence content.*
