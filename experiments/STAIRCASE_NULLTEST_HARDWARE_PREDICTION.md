# Staircase Null-Test: Hardware Prediction (pre-registration)

**Status: PREDICTION, v3, 2026-07-25. Registered before any hardware data exists: no measured results from any flight appear in this document; the planning constants below are inputs to the frozen σ bands only. v1 → v2 → v3 fold two empty design-review rounds of three lenses each (physics, spec, statistics); the operative changes are recorded in the gate script header.**

*Thomas Wicht and Claude. The rate side of the conditional Ramsey: the dephasing-exact-null and its two readings.*

---

## The claim under test

Under local Z-dephasing, a coherence cell |A⟩⟨B| pays only for the sites where A and B differ ([the Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), site-resolved form). Adding one excitation to BOTH patterns at the same site (one step of the η-pairing staircase, [the XY frozen band](XY_FROZEN_BAND.md), the map Φ) leaves the difference set untouched, so the Z-dephasing contribution to the decay-rate difference between the stepped and unstepped cell is exactly zero: for arbitrary site-dependent rates, and also for collective/correlated Z-noise (verified in `simulations/staircase_nulltest_feasibility.py`, check V3, cancellation at 1e-17). What remains in the difference is exactly the flipped site's net amplitude relaxation:

    Δ(m) = Γ↓,m − Γ↑,m = Γ_tot,m · (1 − 2 p∞,m)

with Γ_tot,m and the excited equilibrium population p∞,m readable in-situ (F84-consistent net cooling; the identity was re-derived independently in both design-review rounds, and it is exactly the t = 0 slope of the committed log-ratio). ZZ crosstalk is diagonal, so it enters the FREQUENCY difference (the conditional shift, for a bonded spectator), never the rate difference.

Predecessors: [PRICE_PAIR_HARDWARE_PREDICTION](PRICE_PAIR_HARDWARE_PREDICTION.md) P3 established "pattern, not carrier" and its run 4 measured the ZZ frequency side by conditional Ramsey; [IBM_CONCENTRATOR_RELOADED](IBM_CONCENTRATOR_RELOADED.md) established the matched-pair slope-difference estimator this design commits to. This flight adds the rate side: the exact null of all dephasing in the difference, the in-situ T1-reference identity, and additivity.

## Protocol (frozen)

**Backend and qubits.** ibm_kingston (Heron r2). Three adjacent qubits on one line, target q0 at the chain END: q0 (target), q1 (near spectator, bonded to q0), q2 (far spectator, bonded only to q1). Day-of selection from the fresh calibration among rule-passing lines (floors below); the selection rule is a HARD abort, no override flag.

**Arms.** Four spectator patterns s = (s1, s2) ∈ {00, 10, 01, 11}. Preparation: H on q0; X on each spectator with s = 1. Free delay t; virtual-Z detuning 30 kHz on q0 (rz(2π · 0.030 MHz · t) after the delay); measure q0 in X and in Y (two settings), spectators always in Z. **Dynamical decoupling and any scheduler-inserted echo are explicitly disabled; the delays are bare idle. The runner asserts this and it is a day-of guard.**

**Delays (µs).** Main grid {0, 1, 2, 4, 7, 11, 16, 24, 36, 54, 80} (rates). Dense-early grid {0, 0.5, 1, 1.5, 2, 3, 4} (frequency/ZZ fits ONLY; the main grid aliases a 30 kHz fringe, found and frozen at the gate stage).

**In-situ relaxation block (B-block), bracketed.** Prepare q1 = q2 = |1⟩ (q0 = |0⟩), delays {0, 40, 80, 160, 240, 400, 600} µs, measure all-Z. The full block runs TWICE, as the first and the last circuits of the batch (T1-telegraphing guard). Committed fit, per spectator, per bracket and pooled: z(t) = z∞ + (z0 − z∞)e^(−Γ_tot t) by grid search over Γ_tot ∈ linspace(1/700, 1/50, 400) µs⁻¹ with linear amplitude solve, p∞ = clip((1 − z∞)/2, 0, 0.5); the prediction uses the POOLED fit. The prediction uses these in-situ values, never calibration T1.

**Calibration circuits.** CAL0 = |000⟩, CAL1 = |111⟩, measured in Z: per-qubit asymmetric readout confusion, inverted in analysis.

**Budget.** 14 delays × 4 arms × 2 settings + 2 × 7 B-block + 2 CAL = 128 circuits × 8192 shots = 1,048,576 shots, one Batch. Billing anchor MEASURED, not estimated: the concentrator A-vs-B jobs (908k shots, 227 PUBs, same backend class and Batch mode) billed 257 s = 4.28 min each; scaling by shots and adding the delay content projects ≈ 5.6 QPU min. Re-projected before submission; HARD abort if the projection exceeds 8 QPU min.

## Committed estimator (frozen)

Per arm s, delay t, basis b ∈ {X, Y}: readout-mitigated outcome distribution p̂; the cell estimate is RAW:

    Ĉ_s(t) = Σ_out p̂(out) · (−1)^{x0(out)} · 1(spectator bits = s)     for b = X, plus i · (same for b = Y)

**The estimator divides by TOTAL shots. Post-selection lives in the numerator only.** Renormalizing by the kept count divides out exactly e^(−Γ_m t), the spectator's own survival, and makes the conditional coherence Γ-blind (this is why standard conditional Ramsey reads only frequency; found at the gate stage when the normalized estimator returned Δ ≈ 0).

**Rate differences are MATCHED-PAIR: for each arm s, ordinary least squares on ln(|Ĉ_s(t)| / |Ĉ_00(t)|) over the main-grid points where both moduli exceed 0.02 (frozen); Δ̂(s) = −slope.** Under the claim this log-ratio is exactly free of the common dephasing envelope, whatever its shape (Markovian or quasi-static; per-arm separate exponential fits are NOT used, design review quantified their Gaussian-envelope bias at −1.5σ); the residual spectator-T1 curvature is O(p∞), modeled by the gate and absorbed into the centers m_i. Frequency differences: OLS on unwrap(arg(Ĉ_s/Ĉ_00)) over the dense-early grid. Verdict statistics:

    r1 = Δ̂(10) − Γ̂_tot,1 (1 − 2 p̂∞,1)
    r2 = Δ̂(11) − Δ̂(10) − Δ̂(01)
    r3 = Δ̂(01) − Γ̂_tot,2 (1 − 2 p̂∞,2)
    δf_near = slope of arg(Ĉ_10/Ĉ_00),   δf_far = slope of arg(Ĉ_01/Ĉ_00)

## Predictions and frozen verdict constants

Constants from the committed gate `simulations/staircase_nulltest_gate.py` (v3, seed 31415, N_MC = 2000, MC error ≈ 1.6% per σ; quoted to 2 significant figures). Full planning input vector, inlined: T1 = [322, 142, 237] µs; T2* = [136, 90, 110] µs; p∞ = [1.5, 2.0, 1.5] %; readout (ε01, ε10) = [(0.6, 1.2), (1.0, 1.8), (0.8, 1.5)] %; ζ_shift = −3.9 kHz (the CONDITIONAL-RAMSEY SHIFT, the price-pair run-4 convention; the Hamiltonian coefficient is ζ_c = ζ_shift/4); quasi-static fraction of q0 pure dephasing = 0.5 matched at t_ref = 30 µs (a fixed modeling assumption; the estimator MEAN is insensitive to this split by construction, it shapes only late-point SNR and hence the σ's). Each verdict statistic is CENTERED on its gate mean; the m_i are a small, day-of-rate-dependent curvature systematic (the finite-grid T1-curvature of the log-ratio against the clean-exponential B-block rate, order 1e-5), removed by centering and recomputed by the day-of re-gate:

    m1 = −1.3e-5   s1 = 4.4e-4   /µs
    m2 = +1.6e-5   s2 = 7.0e-4   /µs
    m3 = −1.1e-5   s3 = 4.3e-4   /µs
    sf_near = 4.0e-3   sf_far = 4.1e-3   rad/µs
    s_bracket = [8.5e-5, 5.0e-5] /µs (near, far)

- **P1 (near null).** z1 = |r1 − m1|/s1: HOLDS if z1 ≤ 2; VIOLATED if z1 > 3; else MARGINAL.
- **P2 (additivity).** z2 likewise with (m2, s2).
- **P3 (far null).** z3 likewise with (m3, s3).
- **P4 (ZZ sorting; reported context, not a verdict).** δf_near expected at ζ̂_shift, i.e. −2π·0.0039 rad/µs ≈ −0.025 rad/µs under planning (nonzero at > 3 sf_near; ONLY THE MAGNITUDE is physically predicted, the sign is bra/ket-quadrature convention and is read against the price-pair sign convention after validating the Y-quadrature sign on the CAL/early points); δf_far expected small (residual next-nearest ZZ is a measured context number, not asserted-zero).

**Verdict truth table (frozen).** With the B-block CLEAN (definition below): the claim is **CONFIRMED** iff P1, P2 and P3 all HOLD; **FALSIFIED** iff any of P1/P2/P3 is VIOLATED; otherwise **INCONCLUSIVE**. If the B-block is NOT clean, P1/P3 are void (no T1 reference) and only P2 is read; the experiment is then reported as B-block-invalid regardless of P2. Under the gate's H0-true MC, P(≥1 MARGINAL among P1-P3) = 0.10 (a lone MARGINAL is expected in ~1 of 10 true-null experiments and is not evidence against the law); P(≥1 VIOLATED) = 0.007; P(B-block not clean) = 0.006.

**B-block CLEAN (frozen).** All four per-spectator conditions: (i) pooled-fit Γ_tot not at a grid rail; (ii) p̂∞ ∈ [0, 0.10]; (iii) bracket consistency |Γ̂↓−Γ̂↑ (bracket 1) − (bracket 2)| ≤ 3 s_bracket for that spectator; (iv) pooled-fit RMS residual ≤ 0.02 in z.

**Detection context (gate, H0-true).** Δ̂(10) detected at 15.6σ, Δ̂(01) at 9.5σ. Certification level: s1 corresponds to 5.7% of the near spectator's pure-dephasing rate Γφ,1 = 1/T2*₁ − 1/(2T1₁) = 7.6e-3 /µs (formula pinned: s1/Γφ,1 with the unrounded band, 4.36e-4/7.59e-3 = 5.7%). **Falsification means VIOLATED (3σ): dephasing leakage into the rate difference at ≥ 17% of the spectator's pure-dephasing rate would be flagged FALSIFIED; leakage between 11% (2σ) and 17% lands MARGINAL → INCONCLUSIVE, honestly.**

## Day-of guards (hard aborts)

1. **Fresh calibration pull.** Floors: every selected qubit T2* ≥ 60 µs, T1 ≥ 100 µs, readout error ≤ 2%.
2. **Line rule.** Three adjacent qubits, q0 at an end of the selected 3-line; q0-q2 not directly coupled on the live coupling map.
3. **Mandatory day-of re-gate (committed addendum).** Before the Batch opens, `staircase_nulltest_gate.py` is re-run with the day-of measured inputs (T1, T2*, readout from the fresh calibration; p∞ from CALIBRATION-class data only, never from the batch's own B-block, which would be circular; ζ_shift if available); the resulting m_i, s_i, sf, s_bracket REPLACE the planning constants above and are committed as an addendum to this document. The committed day-of constants govern the verdict. The re-gate also prints the min-kept-points floor (HARD abort if any arm keeps < 5 main-grid points) and the B-block self-void rate. (Statistics review: the σ bands scale ~1.25× when T2* drops 30% within the allowed floor window; planning bands do not preserve 2σ/3σ coverage across that window; a re-gate does, and it consumes only calibration inputs, so it is not data-peeking.)
4. **Band-validity window (multi-dimensional).** If any selected qubit's day-of T1 ∉ [100, 500] µs, or T2* < 60 µs, or readout > 2%, the design itself is invalid: abort, redesign, re-register (grids and G-grid were sized for this window; the T1 floor coincides with guard 1's).
5. **DD/echo off.** The runner asserts no dynamical decoupling and no scheduling-inserted echo on the delays; bare idle only. Backend max-delay must admit 600 µs.
6. **Billing projection** re-measured against comparable-run billing (anchor above); abort above 8 QPU min.
7. **Commit before flight.** This document (with the day-of addendum) AND both cited scripts (`simulations/staircase_nulltest_gate.py`, `simulations/staircase_nulltest_feasibility.py`) committed, real hash, BEFORE the Batch is opened; the committed constants govern the verdict even against any runner printout.

## Honesty notes

- Hardware-wise the circuits are standard conditional Ramsey plus relaxation and calibration blocks; the novelty is entirely in the claim tested (exact dephasing null in the rate difference, the F84 net-cooling identity for the splitting, additivity, near/far ZZ sorting) and in the RAW matched-pair estimator that keeps the rate signal alive.
- The dephasing null is proven for arbitrary local and collective Z-noise. It is NOT shielded against correlated amplitude damping (cross-relaxation between target and spectator); on Heron this is expected negligible, and it would enter r1/r3 as real physics, not as an analysis artifact.
- The gate models uneven Markovian Z-dephasing, a quasi-static (Gaussian-envelope) dephasing component on the target, T1 with heating, ZZ, and asymmetric readout. It does not model leakage to |2⟩ or mid-batch T1 telegraphing; the first would appear as structure in the fit residuals, the second is guarded by the bracketed B-block (clean-condition iii). Neither is a licence to move a band after the fact.
- Readout-confusion inversion commutes with the numerator-only post-selection (linear estimator; verified in design review), and readout error is t-independent, so it enters the contrast, not the slopes.
- The B-block prepares BOTH spectators excited (|011⟩) while the arms excite one at a time; if excitation-dependent T1 or spectator-spectator cross-relaxation exists, the reference Γ_tot,m can differ slightly from the single-excited in-situ value. Expected negligible on Heron, partially caught by clean-condition (iv); the B-block fit residual structure is the day-of eye on it.

*Gate lineage (committed together with this document, guard 7): `simulations/staircase_nulltest_feasibility.py` (v0 scout, 10/10 GREEN: exactness, uneven and collective-Z blindness at 1e-17, ZZ sorting, XY-on robustness, and the discovery of the post-selection normalization trap) and `simulations/staircase_nulltest_gate.py` (v3: flown construction, counts level, matched-pair estimator, quasi-static component, bracketed B-block; the constants above). Theory: [XY_FROZEN_BAND](XY_FROZEN_BAND.md), [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), F84.*
