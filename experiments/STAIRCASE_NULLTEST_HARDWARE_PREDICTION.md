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
4. **Band-validity window (multi-dimensional).** If any selected qubit's day-of T1 ∉ [100, 500] µs, or the TARGET qubit's T2* < 60 µs (scope per Amendment 1: spectator T2* enters no observable and does not gate), or any readout > 2%, the design itself is invalid: abort, redesign, re-register (grids and G-grid were sized for this window; the T1 floor coincides with guard 1's).
5. **DD/echo off.** The runner asserts no dynamical decoupling and no scheduling-inserted echo on the delays; bare idle only. Backend max-delay must admit 600 µs.
6. **Billing projection** re-measured against comparable-run billing (anchor above); abort above 8 QPU min.
7. **Commit before flight.** This document (with the day-of addendum) AND both cited scripts (`simulations/staircase_nulltest_gate.py`, `simulations/staircase_nulltest_feasibility.py`) committed, real hash, BEFORE the Batch is opened; the committed constants govern the verdict even against any runner printout.

## Honesty notes

- Hardware-wise the circuits are standard conditional Ramsey plus relaxation and calibration blocks; the novelty is entirely in the claim tested (exact dephasing null in the rate difference, the F84 net-cooling identity for the splitting, additivity, near/far ZZ sorting) and in the RAW matched-pair estimator that keeps the rate signal alive.
- The dephasing null is proven for arbitrary local and collective Z-noise. It is NOT shielded against correlated amplitude damping (cross-relaxation between target and spectator); on Heron this is expected negligible, and it would enter r1/r3 as real physics, not as an analysis artifact.
- The gate models uneven Markovian Z-dephasing, a quasi-static (Gaussian-envelope) dephasing component on the target, T1 with heating, ZZ, and asymmetric readout. It does not model leakage to |2⟩ or mid-batch T1 telegraphing; the first would appear as structure in the fit residuals, the second is guarded by the bracketed B-block (clean-condition iii). Neither is a licence to move a band after the fact.
- Readout-confusion inversion commutes with the numerator-only post-selection (linear estimator; verified in design review), and readout error is t-independent, so it enters the contrast, not the slopes.
- The B-block prepares BOTH spectators excited (|011⟩) while the arms excite one at a time; if excitation-dependent T1 or spectator-spectator cross-relaxation exists, the reference Γ_tot,m can differ slightly from the single-excited in-situ value. Expected negligible on Heron, partially caught by clean-condition (iv); the B-block fit residual structure is the day-of eye on it.

---

## Amendment 1 (2026-07-26, pre-data; backend switch and floor scope)

Recorded BEFORE any science data exists. The 2026-07-25 Kingston submission (job d9iidarhdfks73ch30ag) sat 8+ hours in a 167-job queue and was CANCELLED at 0 billed QPU seconds; no counts were ever produced. The flight moves to the twin Heron r2 backend, **ibm_marrakesh** (the repo's precedent for a pinned twin is the record-parity pre-registration; Marrakesh is also where the ζ = −3.9 kHz context value was measured, so P4 is read on its home chip).

**Floor scope correction, with the physics reason:** guard 1's "every selected qubit T2* ≥ 60 µs" is hereby scoped to the TARGET qubit only. Spectator dephasing enters no observable of this experiment: the arm cells carry the spectators diagonally (a Z-channel acts trivially there, the structural fact all three review rounds verified), and the B-block reads populations, which are dephasing-blind. A spectator T2* floor is therefore physically vacuous, and under the frozen 0.4 T2*-derating it excluded every line on both twin backends for no physical reason. Spectator T2* values are still recorded in the day-of snapshot and in the re-gate inputs. All other floors (T1 window, readout, target T2*) are unchanged; the runner's presubmit check is updated to match (comment "Amendment 1" at the check).

The Kingston day-of addendum below is retained as the record of the aborted attempt; the GOVERNING constants for the actual flight are in the Marrakesh addendum that follows it.

## Day-of addendum (2026-07-25, KINGSTON, superseded by the cancelled queue; retained as record)

**Backend and line.** ibm_kingston, fresh calibration `last_update 2026-07-25 20:33:34+02:00` (pulled 23:06 local, read-only). Selected line: **q0 = 104 (target, chain end), q1 = 105 (near spectator), q2 = 106 (far spectator)**; live coupling map confirms bonds (104,105) and (105,106) and NO (104,106). All floors pass (see inputs).

**Day-of re-gate inputs** (the DAY-OF INPUT block of the committed gate, seed 31415, N_MC = 2000, unchanged construction):

- T1 = [242, 190, 118] µs (backend properties, within [100, 500] ✓)
- T2* = [115, 70, 66] µs = **0.4 × reported T2 [288, 174, 165]**. Derating note, recorded before any science data: `props.t2` on this stack is not guaranteed to be Ramsey T2*, and this pipeline's own price-pair flight measured T2* = 45/88/61 µs against calibration T2 ≈ 136-220 (ratio ≈ 0.3-0.5). 0.4 is the frozen midpoint of that measured anchor. All derated values ≥ the 60 µs floor ✓. The same values are passed to the runner as `--t2star`.
- p∞ = [1.5, 2.0, 2.0] % (calibration-class assumption; no direct day-of source; NOT taken from the batch B-block)
- readout (ε01, ε10) = [(0.15, 1.37), (0.68, 1.42), (0.73, 1.07)] % (prob_meas1_prep0 / prob_meas0_prep1, day-of properties)
- ζ_shift = −3.9 kHz (price-pair Marrakesh measured value; context only)

**Day-of frozen constants (GOVERNING; gate output verbatim, quoted to 2 significant figures):**

    m1 = +4.7e-6   s1 = 5.1e-4   /µs
    m2 = −2.3e-5   s2 = 8.7e-4   /µs
    m3 = −1.6e-5   s3 = 5.6e-4   /µs
    sf_near = 4.1e-3   sf_far = 3.9e-3   rad/µs
    s_bracket = [6.3e-5, 1.0e-4] /µs (near, far)

Gate verdict on these inputs: PASS. Detection context: Δ̂(10) at 9.9σ, Δ̂(01) at 14.7σ (the far spectator's shorter T1 gives the larger splitting). Family rates: P(≥1 MARGINAL) = 0.10, P(≥1 VIOLATED) = 0.009, P(B-block not clean) = 0.005. Min kept points 11. Certification: s1/Γφ,1 with Γφ,1 = 1/70 − 1/380 = 1.17e-2 /µs gives 1σ = 4.4% of the near spectator's pure-dephasing rate (2σ = 8.8%, VIOLATED threshold 3σ = 13%).

## Day-of addendum 2 (2026-07-26, MARRAKESH; these constants GOVERN the verdict)

**Backend and line.** ibm_marrakesh (Heron r2), fresh calibration `last_update 2026-07-26 05:26:42+02:00`. Selected line under Amendment 1: **q0 = 35 (target, chain end), q1 = 34 (near spectator), q2 = 33 (far spectator)**; live coupling map confirms bonds (35,34), (34,33), NO (35,33). Queue at selection: 4 pending.

**Day-of re-gate inputs** (DAY-OF INPUT block of the committed gate, seed 31415, N_MC = 2000, construction unchanged):

- T1 = [238, 168, 129] µs (within [100, 500] ✓)
- T2* = [124, 32, 62] µs = 0.4 × reported T2 [310, 80, 154] (same frozen derating as addendum 1; target 124 ≥ 60 ✓; spectator values recorded, not gating, per Amendment 1)
- p∞ = [1.5, 2.0, 2.0] % (calibration-class assumption; not from the batch B-block)
- readout (ε01, ε10) = [(0.29, 0.34), (0.24, 0.39), (0.10, 2.05)] %
- ζ_shift = −3.9 kHz (measured on THIS backend, price-pair run 4)

**Day-of frozen constants (GOVERNING; gate output, 2 significant figures):**

    m1 = −1.0e-5   s1 = 4.5e-4   /µs
    m2 = +1.2e-5   s2 = 7.5e-4   /µs
    m3 = −3.0e-5   s3 = 4.7e-4   /µs
    sf_near = 4.0e-3   sf_far = 4.1e-3   rad/µs
    s_bracket = [7.1e-5, 9.1e-5] /µs (near, far)

Gate verdict: PASS. Detection context: Δ̂(10) at 12.6σ, Δ̂(01) at 15.8σ. Family rates: P(≥1 MARGINAL) = 0.11, P(≥1 VIOLATED) = 0.007, P(B-block not clean) = 0.008. Min kept points 11. Certification: Γφ,1 = 1/32 − 1/336 = 2.8e-2 /µs, so 1σ = 1.6% of the near spectator's pure-dephasing rate (3σ = 4.8%); the derated spectator T2* makes this the tightest certification yardstick of the three gates.

---

## HARDWARE RECORD (2026-07-26, ibm_marrakesh, flown once)

**Flight.** Job `d9ip713hdfks73chbfg0`, one Batch, line (35, 34, 33), 128 PUBs × 8192 shots, 0 PUB failures, billed 324 s = 5.4 QPU min (projection 5.6, cap 8). Counts persisted before analysis; archive complete (all 128 PUBs, day-of snapshot, constants source). An earlier same-design Kingston submission (`d9iidarhdfks73ch30ag`) was cancelled unstarted at 0 billed seconds after 8+ hours in a 167-job queue; see Amendment 1.

**Committed verdict: B-BLOCK-INVALID.** Independently re-derived from the raw counts by an empty-context audit (machine-precision agreement with the persisted analysis, ≤ 8e-15 on every quantity). The B-block CLEAN conditions fail on both spectators: near (q34) fails bracket consistency, |bc1| = 3.87e-4 > 3 s_bracket = 2.12e-4 /µs (its effective T1 ran 160.7 → 172.3 µs between batch start and end); far (q33) fails bracket consistency, |bc2| = 7.90e-4 > 2.72e-4 /µs (T1 ran 132.9 → 115.1 µs, relaxation speeding up mid-batch) AND the pooled-fit residual, rms = 0.028 > 0.02, with a structured, non-single-exponential bow. Per the frozen truth table, P1/P3 are void (no trustworthy in-situ T1 reference) and the experiment is reported B-block-invalid regardless of P2. This is the mid-batch T1-telegraphing failure mode the bracketed B-block was added to catch (design round 1); the guard held, including against the temptation of the numbers below.

**The voided numbers, recorded as observations, not verdicts:**

    Δ̂(10) = +5.817e-3 /µs   in-situ pred = +5.836e-3   r1 = −1.9e-5   z1 = 0.02  (would HOLD)
    Δ̂(11) − Δ̂(10) − Δ̂(01):  r2 = +1.1e-4   z2 = 0.13                            (HOLDS; P2 is B-block-free)
    Δ̂(01) = +8.757e-3 /µs   in-situ pred = +7.558e-3   r3 = +1.2e-3   z3 = 2.60  (would be MARGINAL)
    δf_near = −0.02373 rad/µs = −3.78 kHz   (P4 context: pred −3.90 kHz, sign and magnitude reproduced, 0.19 sf under)
    δf_far  = +0.00575 rad/µs = +0.91 kHz   (P4 context: 1.4 sf_far from zero, within band)

The near-arm dephasing-null lands on its reference to 0.3% and the additivity statistic holds at 0.13σ; the far-arm excess (+16%) sits exactly on the spectator whose T1 reference telegraphed hardest and carries the structured residual, consistent with the invalidation rather than with dephasing leakage. None of this is promoted to a verdict: the committed rule governs.

**Standing after this flight.** The claim is neither confirmed nor falsified; the instrument worked end-to-end and the guard did its job. A re-flight draws a new T1-weather card at ~5.4 QPU min; if one is flown, candidate hardening (requires a re-registered amendment BEFORE that flight): interleave B-block circuits through the batch instead of two end brackets, and/or a single-excited B-block variant per spectator (the |011⟩-vs-single-excited caveat in the honesty notes, now with a matching observed residual signature on q33).

---

## Amendment 2 (2026-07-26, pre-data for flight 2; B-block hardening)

Registered BEFORE any flight-2 data exists, in response to flight 1's B-BLOCK-INVALID (the record above): the T1 reference, not the claim, failed. Two changes, both confined to the B-block; arms, estimator, grids, verdict machinery, clean-condition FORM, and truth table are unchanged.

1. **Single-excited B-block, one preparation per spectator.** Prep 1 = |010⟩ (near excited), prep 2 = |001⟩ (far excited), each at the same 7 delays, all-Z readout; spectator q reads ONLY its own prep. This removes the |011⟩ both-excited caveat of the honesty notes, whose predicted signature (excitation-dependent structure in the fit residual) flight 1 observed on q33 (rms 0.028).
2. **Interleaved placement, frozen pattern.** Each prep flies as two half-sets; the B sequence half0 = [prep1 t0..t6, prep2 t0..t6], half1 = same again. With the 14 arm delay-groups sorted ascending, the PUB order is: for k = 0..13, the 8 arm circuits of group k, then B-sequence circuits 2k and 2k+1; CAL0, CAL1 last. Half0 thus spans the first seven groups, half1 the last seven: the reference samples the same T1 weather as the arms, and the half-split keeps drift detection in FORM (bracket-consistency statistic unchanged, s_bracket refrozen). Two honest tradeoffs, disclosed: the interleaved half-split has a shorter time lever than flight 1's end brackets, so condition (iii) is less sensitive per σ to slow monotonic drift (condition (iv), the whole-batch rms, partly compensates, and the pooled reference is now batch-matched to the arms, which is the point); and within each half, prep 1 samples earlier groups than prep 2, so under global drift the two spectators' references are weighted toward slightly different parts of the batch.

Budget becomes 112 + 28 + 2 = 142 circuits × 8192 = 1,163,264 shots, projected ≈ 6.6 QPU min (the runner's two-term projection, the governing number; pure anchor-scaling of the flight-1 billing gives 6.0), still under the 8-min hard abort. Gate v4 (committed) implements this construction; under the planning inputs it passes with detection 15.1σ / 9.8σ, self-void rate 0.004, and visibly smaller centering systematics (m3: +8.7e-7 vs v3's −1.1e-5, the single-excited reference removing curvature). The governing constants for flight 2 will be a day-of addendum 3 from the v4 gate, committed before that Batch opens; all flight-1 guards, Amendment 1 included, remain in force.

## Day-of addendum 3 (2026-07-26, MARRAKESH, flight 2 under Amendment 2; these constants GOVERN)

**Backend and line.** ibm_marrakesh, calibration `last_update 2026-07-26 05:26:42+02:00` (unchanged since addendum 2; re-verified at selection). Line (35, 34, 33) re-selected as the top rule-passing candidate; bonds re-confirmed. **Day-of re-gate inputs identical to addendum 2** (T1 = [238, 168, 129]; T2* = [124, 32, 62] = 0.4 × reported, target-only gating per Amendment 1; p∞ = [1.5, 2.0, 2.0] %; readout as addendum 2; ζ_shift = −3.9 kHz); the construction is gate v4 (Amendment 2).

**Day-of frozen constants (GOVERNING, gate v4 output, 2 significant figures):**

    m1 = −1.1e-5   s1 = 4.6e-4   /µs
    m2 = +3.8e-6   s2 = 7.9e-4   /µs
    m3 = −1.1e-5   s3 = 4.9e-4   /µs
    sf_near = 4.0e-3   sf_far = 4.0e-3   rad/µs
    s_bracket = [6.8e-5, 9.1e-5] /µs (near, far)

Gate verdict: PASS. Detection context: Δ̂(10) at 12.4σ, Δ̂(01) at 15.2σ. Family rates: P(≥1 MARGINAL) = 0.11, P(≥1 VIOLATED) = 0.007, P(B-block not clean) = 0.006. Min kept points 11.

---

## HARDWARE RECORD 2 (2026-07-26, ibm_marrakesh, flight 2 under Amendment 2)

**Flight.** Job `d9iptdoii2cc73edtm20`, one Batch, line (35, 34, 33), 142 PUBs × 8192 shots, 0 PUB failures, billed 380 s = 6.3 QPU min (projection 6.6, cap 8). Counts persisted before analysis; archive complete; an empty-context audit re-derived the full committed pipeline from raw counts at ≤ 1e-15 on every quantity.

**Committed verdict: B-BLOCK-INVALID** (again, by a different failure than flight 1). The hardened B-block PASSED everything Amendment 2 fixed: bracket consistency now clean on BOTH spectators (bc = −2.8e-5 / −1.2e-4, was +3.9e-4 / −7.9e-4), far rms cut 0.028 → 0.0209. But 0.0209 > 0.02: the far spectator (q33) fails clean-condition (iv) by 0.0009, a razor-thin margin on a still-structured residual, so P1/P3 are void per the frozen truth table. P2 (additivity, B-block-free) HOLDS at z2 = 1.19.

**The voided numbers, recorded as observations, not verdicts:**

    Δ̂(10) = +1.3346e-2 /µs   in-situ pred = +1.0962e-2   r1 = +2.38e-3   z1 = 5.20  (would be VIOLATED)
    Δ̂(01) = +1.0247e-2 /µs   in-situ pred = +0.9574e-2   r3 = +6.7e-4    z3 = 1.39  (would HOLD)
    δf_near = −0.02617 rad/µs = −4.17 kHz  (P4: pred −3.9 kHz, sign and magnitude again reproduced, 0.42 sf)
    δf_far  = +0.01073 rad/µs = +1.71 kHz  (P4 context: 2.7 sf_far from zero, largest of both flights)

**Diagnosis (audit, quantitative; no bands moved).** A genuine device excursion: q34's in-situ effective T1 collapsed to 89.3 µs against the 05:26 calibration's 167.7 (and against flight 1's in-situ ~171 two hours earlier); q33 fell 129.5 → 102. The device thereby left the design's validity window (T1 floor 100 µs) MID-BATCH, invisibly to guard 4 (which reads calibration) and visibly to the B-block gate, which is the intended backstop and fired. The voided z1 = 5.2 is consistent with this weather, not with dephasing leakage: the arm grid (0-80 µs) and the reference grid (40-600 µs) weight a non-exponential, collapsing T1 differently by construction (arm-window effective T1 ≈ 75 µs vs reference 89 µs → r1 > 0), and the sign ordering r1 > r3 > 0 matches Amendment 2's disclosed prep staggering under drift toward faster relaxation. The audit names the residual design tension honestly: with the reference window barely overlapping the arm window, a rate leak and a non-exponential reference are CONFOUNDED by construction; the separation test would be reference points inside the arm window (short B delays), a candidate Amendment 3 if a flight 3 is ever flown, together with an immediately-pre-batch calibration re-pull.

**Standing after two flights.** The claim is neither confirmed nor falsified. Twice the guard refused to read a verdict off an invalid T1 reference, once against a 0.02σ beauty (flight 1) and once against a 5.2σ scare (flight 2); the additivity statistic held both times; the ZZ frequency sorting reproduced sign and magnitude both times. What the experiment certifiably needs is not a better claim but calmer T1 weather, and a reference that samples the arms' own time window.

*Gate lineage (committed together with this document, guard 7): `simulations/staircase_nulltest_feasibility.py` (v0 scout, 10/10 GREEN: exactness, uneven and collective-Z blindness at 1e-17, ZZ sorting, XY-on robustness, and the discovery of the post-selection normalization trap) and `simulations/staircase_nulltest_gate.py` (v3: flown construction, counts level, matched-pair estimator, quasi-static component, bracketed B-block; the constants above). Theory: [XY_FROZEN_BAND](XY_FROZEN_BAND.md), [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), F84.*
