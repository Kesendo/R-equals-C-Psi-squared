# The Concentrator A-vs-B Mechanism Test (pre-registration draft, v2 after physics review)

<!-- Keywords: selective dynamical decoupling mechanism test, gate-error avoidance vs
noise contrast, uniform T2 chain control, engineered sink injected dephasing random RZ,
concentrator Sum-MI ibm hardware, noise-seeded observable Heisenberg fixed point,
A vs B interpretation April test, R=CPsi2 -->

**Status:** DESIGN + PRE-REGISTRATION DRAFT, v2 (2026-07-05, not yet flown). v1 went
through a physics-first review the same day (verdict FIX-THEN-FLY, two blockers, all
findings incorporated below; the review's exact-simulation numbers are quoted where
they set design choices). This is the open item named in
[IBM Concentrator](IBM_CONCENTRATOR.md) ("A vs B test", planned April 9, never run)
and handed to the DD community as the falsifiable prediction of
[Selective Decoupling Selection Rule](../docs/outbound/SELECTIVE_DECOUPLING_SELECTION_RULE.md)
§5. The simulator gate (Section "Simulator validation") must be recorded here and this
document committed BEFORE anything flies.
**Date:** 2026-07-05
**Authors:** Thomas Wicht, Claude (Fable 5)
**Budget context:** ≈ 100 QPU min available (Tom, 2026-07-05); this design targets ≈ 30-40 min.

## The question (from the March run's own honesty section)

The March 2026 ibm_torino result (selective DD beats uniform DD 2.0-3.2×) has two
readings. **A (gate-cost):** leaving the terrible qubit (T2 = 5 µs) undecoupled wins
because DD gates on a lost cause only add errors. **B (concentrator):** the win is the
noise CONTRAST between the noisy edge and the protected interior, the hardware face of
the γ_edge = N·γ_base − (N−1)·ε profile. The discriminating experiment, stated in both
source docs: a chain where NO qubit is a natural sink, plus an ENGINEERED sink.

## The structural fact that frames the observable (review finding, verified exactly)

|+⟩^⊗5 is an exact fixed point of every Heisenberg bond gate (the bond unitary is
e^{iθ/2}(cos θ·I − i sin θ·SWAP), and identical product states are SWAP-invariant), so
under UNIFORM noise the state stays an exact product at all times and Sum-MI ≡ 0
(verified ≤ 2·10⁻¹⁵ through 10 Trotter steps in exact 32×32 simulation). **Sum-MI in
this experiment is not a signal that noise degrades; it is a signal that only
NON-uniform noise creates.** Two consequences: (i) the paired sink-vs-no-sink
comparison, not any ratio against uniform DD, carries the physics; (ii) the March
uniform-DD row (Sum-MI 0.013-0.048) is plausibly dominated by the MI estimator's
shot-noise bias floor (+0.014-0.028 at 4000 shots, measured through the March
reconstruction pipeline in the review), a retro-insight recorded here and to be noted
in the March doc after this experiment reports.

## Design (v2)

**Chain:** one 5-qubit connected line, uniformly good by same-day API calibration
(`backend.properties()`): min T2echo ≥ 150 µs, max/min T2echo ≤ 2, readout error
≤ 2%; ranked by uniformity = min T2 / max T2, tie-break higher min T2. Because the
loaded variable is T2* (the undecoupled qubit sits at T2*, not T2echo), position 0
(the sink end) is chosen as the edge with the HIGHER estimated T2*, and the residual
natural contrast (min interior T2echo / edge T2*) is recorded as a design parameter.
The mirrored-end repeat is pre-registered as the first follow-up.

**Five configurations, same Trotterized Heisenberg circuit as March (|+⟩^5, RXX+RYY+RZZ
per bond per step, Trotter steps 2..10, 9-basis pair tomography, Sum-MI over the 4
neighbour pairs):**

1. `no_dd` (March baseline), 4000 shots/circuit.
2. `uniform_dd`: X-X DD on all 5 (March baseline), 4000 shots.
3. `selective_dd`: DD on the interior 4, edge undecoupled, no sink, 4000 shots.
4. `selective_dd_sink`: as 3 plus the engineered sink on the edge.
5. `uniform_dd_sink`: as 2 plus the same sink (the same-layout partner of config 2;
   it retires the objection that the selective LAYOUT, not the contrast, carries any
   effect; X-X DD conjugates RZ(φ) → RZ(−φ) and the phase set below is
   negation-symmetric, so the sink composes with DD unchanged).

**The engineered sink (v2; the v1 construction was a review blocker):** per Trotter
step, a virtual RZ(φ) on the edge qubit. K = 16 STRUCTURED frozen phase paths per
time point: 4 random offsets × a per-step-permuted equispaced quartet
{0, π/2, π, 3π/2}, negation-symmetric; the SAME 16 paths shared across all 9
tomography bases of that time point (seeds pre-registered in the runner output);
250 shots per instance, counts pooled per (t, basis) to 4000 effective. Review
numbers that force this: K = 4 with fresh draws per basis inflates Sum-MI by
+35-118% (the nine bases would measure nine different ensembles, and the pooled
object at small K is shared classical randomness, which even Reading A boosts);
the structured K = 16 performs like iid K ≈ 32 (residual bias +0.004-0.016,
below the effect and inside the counts-level bands). Implementation: ONE
parameterized circuit per (t, basis), transpiled once with fixed seed_transpiler;
config 3 = the φ⃗ = 0 binding, configs 4/5 = the 16 frozen vectors as SamplerV2
parameter bindings; pre-flight asserts identical op counts and durations across
bindings. The RZ sits adjacent to the last entangling gate of its step, never
splitting an idle window.

**Sink strength, stated:** uniform-quartet scramble is the strobe CEILING (per-step
phase retention r = |E e^{iφ}| = 0, continuous-equivalent γ_sim = −ln r → ∞ per
step). The review's exact Trotter sweep shows Sum-MI rises MONOTONICALLY in per-step
phase variance and saturates at this ceiling (no Zeno turnover at θ = 1 rad of
coherent exchange per step), so the ceiling is the maximum-effect operating point;
it is a different object from the concentrator formula's finite-budget optimum
γ_edge = N·γ_base, and the physical Trotter-step duration (from the scheduled
circuits) is recorded so the continuous-γ mapping can be stated. An optional
partial-strength dose point (σ ≈ 0.9 rad) is pre-registered as a follow-up
fingerprint (B predicts a monotone RISE of the paired effect with injection
variance; a noise-is-harm frame predicts a fall), not flown in run 1.

**Job structure:** ALL configs' circuits round-robin-interleaved in ONE SamplerV2
job (drift on the minutes scale must hit all configs alike; if a PUB limit forces a
split, split by time point, never by config). `backend.properties()` snapshotted
before and after. Two readout-cal circuits recorded (not used by the primary March
pipeline; for a later mitigated re-analysis).

**Circuit/execution count:** 45 + 45 + 45 + 45×16 + 45×16 + 2 = 1577 executions
in 227 PUBs (908k shots), sink instances at 250 shots, the rest at 4000. QPU
estimate (printed by the runner pre-flight): Model B (per-PUB anchor) 13.2 min,
the planning number; Model C (no parameter-set amortization) 92 min, the flagged
worst case. Flight rule, pinned: fly on the Model-B planning number against the
available budget; the Model-C tail risk is accepted explicitly at the go
decision (recorded in the RECORD section). The earlier K = 8 fallback is struck:
the runner's structured construction is K = 16 by design and no fallback is
implemented.

## Pre-registered predictions and verdict logic (v2)

**Primary (the paired statistics; gate-, schedule-, transpiler-identical pairs, so
hardware junk and estimator bias enter both sides alike):**

    Δ_sel(t) = SumMI(selective_dd_sink) − SumMI(selective_dd)
    Δ_uni(t) = SumMI(uniform_dd_sink)  − SumMI(uniform_dd)

- **Reading B (concentrator):** Δ_sel > 0 beyond its counts-level null band on ≥ 3
  of 5 time points (the sink CREATES interior correlations; the review's exact
  simulation of the flying construction predicts the effect at multiples of the
  null spread); Δ_uni > 0 likewise (the sink needs only the contrast, not the
  selective layout).
- **Reading A (gate-cost):** both Δ ≤ 0 within bands (injected noise cannot create
  signal; the March advantage was gate-economics on a lost cause).
- **Scope, stated before the shot:** B-CONFIRMED here means "the concentrator
  mechanism is live and engineerable on this hardware, and B is SUFFICIENT to
  explain March"; it does not retroactively prove March WAS B (Q85's gate-cost
  channel also existed). A-CONFIRMED corrects the outbound adapter §5 honestly.
- **Exhaustiveness clause:** ANY outcome pattern not matching the B rule
  (both Δ legs beyond band, ≥ 3 of the returned time points each) or the A rule
  (both legs within bands at all returned time points) is INCONCLUSIVE, named
  as such, with the follow-up stated; this includes significantly NEGATIVE Δ,
  mixed legs (one fires, the other does not), and partial job returns (if a
  time point drops for all configs, the "of 5" rescales to the returned count,
  minimum 3 points for any verdict).
  A dephasing-only exact simulation already produces the boost in the ideal limit
  (review, Sum-MI 0.002 → 0.22 at t = 2), so the hardware question is whether the
  mechanism SURVIVES real junk at the predicted size, which is exactly the
  adapter's open question.

**Context lines (reported, not verdict-carrying):** R_nosink(t) =
SumMI(sel)/SumMI(uni) and R_sink(t) = SumMI(sel+sink)/SumMI(uni). Both are
normalized by a junk-plus-bias floor with no theory value (the structural fact
above); they connect to the March table and nothing more. The v1 "shared null
R_nosink ≈ 1" logic is retired: on a uniform line the undecoupled edge is a weak
natural sink at T2* (so B predicts R_nosink slightly > 1, A slightly < 1, both
within the floor's noise), and no verdict rests on it.

**Bands:** every band is computed counts-level through the SAME estimator and the
SAME frozen sink construction as the flight (the MI bias floor, +0.014-0.028 at
4000 shots, is state-dependent and does not cancel in ratios; it enters bands and
data identically only if the bands go through the pipeline). The simulated bias
floor is drawn on every Sum-MI plot. Bootstrap (≥ 200 multinomial resamples) gives
the CIs; the null spread of Δ under a no-sink world sets the 2σ thresholds.

## Simulator validation (the gate; recorded 2026-07-05, runner `run_ab_test.py`)

The runner (external tomography pipeline, built to this spec the same day, smoke +
full modes) records, counts-level with the LITERAL frozen sink construction:

- **Predicted effect on the FLIGHT chain (see chain paragraph below; seed
  20260705):** Δ_sel = +0.30 / +0.31 / +0.13 / +0.13 / +0.05 at t = 1..5, all
  five bootstrap CIs excluding 0; Δ_u similar (+0.39 to +0.05); R_boost ≈
  3.5-11.8. On the earlier synthetic-uniform gate run (T1 = 250 / T2 = 180 µs)
  the per-pair table shows the largest sink-created MI at the pair FARTHEST
  from the sink (0.103 at (3,4) vs 0.012 at (0,1), t = 3): the created
  correlation propagates through the chain, the concentrator transport
  signature.
- **Null bands (no-sink world, through the same estimator), BINDING at N = 100
  on the flight chain [109, 108, 107, 106, 105] (recorded before the shot,
  `ab_test_null_20260705_172115.json`; the pooled p99 of |Δ|, not a 2σ; the
  pooled band is anti-conservative at early t and conservative at late t,
  harmless under the ≥ 3-of-5 rule at 0/100 null false-fires):** Δ significant
  outside **± 0.0486**, Δ_u outside **± 0.0462**; ratio nulls R_boost ± 4.2,
  R_sink ± 8.6, R_nosink ± 12.9 (context only). Pre-shot check PASSED: the
  predicted Δ on this chain clears the band at t = 1-4 (6.2× / 6.3× / 2.7× /
  2.6×), marginal at t = 5 (1.1×), so B-CONFIRMED is reachable under the
  ≥ 3-of-5 rule with margin. Measured MI bias floor on this chain: mean 0.034
  (median 0.034), drawn on all plots. (An earlier N = 100 pass, ± 0.047 /
  ± 0.053, was computed on the Q66-70 chain that later failed the day's rule;
  superseded together with that chain.)
- **Aer parity (the actual circuits, noiseless): PASS.** All four φ⃗ = 0-equivalent
  arms give Sum-MI exactly 0.0000 (the fixed-point theorem, visible in the real
  circuits); the bound sink deviates as it must (+0.14-0.50 bits noiseless).
- **Chain selection (same-day API calibration, 2026-07-05), with the
  stale-calibration trap demonstrated LIVE:** ibm_marrakesh FAILS the
  uniformity rule (all 428 five-qubit chains rejected on min-T2echo ≥ 150 µs;
  the runner refuses to relax, as designed). ibm_kingston at the 12:58
  calibration passed with 13 chains (best Q66-70, uniformity 0.84); by the
  ~16:45 recalibration that chain READ T2echo 118/239/176/121/117 and FAILED
  the rule (30% intraday drift within three hours), so it and its bands are
  superseded. **The FLIGHT chain is the fresh 16:45 selection:**
  [109, 108, 107, 106, 105], T2echo = 195/172/304/194/174 µs, uniformity 0.565
  (max/min 1.77, inside the ≤ 2 rule), readout within rule, sink edge Q109
  (the higher-T2* end, T2* ≈ 78 µs), residual natural contrast 2.21×
  recorded; one alternative candidate passed. The runner now HARD-ABORTS at
  submission if the loaded chain violates the rule on day-of calibration (no
  override flag, by design). Device for the shot: ibm_kingston.
- **QPU estimate (printed pre-flight by the runner, no submission without
  --yes):** 1577 executions in 227 PUBs, 908k shots. Model A (pure shots) 1.5
  min; Model B (per-PUB anchor) 13.2 min, the planning number; Model C (worst
  case, no parameter-set amortization) 92 min, flagged. The flight decision is
  taken against the printed estimate and the remaining budget.

## HARDWARE RECORD (2026-07-05, ibm_kingston, job d957dfcql68s73caa4q0)

Flown ≈ 17:22-17:40 local, chain [109, 108, 107, 106, 105], all pre-flight
asserts PASS (day-of rule, frozen phases regenerated, binding identity,
added-X ratio 0.833 ≈ 4/5); the go was taken on the Model-B planning number
with the Model-C tail risk accepted (Tom, "Du kannst starten"). Data:
[`data/ibm_ab_test_july2026/`](../data/ibm_ab_test_july2026/) (hardware master
JSON + the binding null, chain, and aer-parity artifacts). This RECORD is a
transcription of the runner's output; its empty-session pass is the next
session's first move and is named here as pending.

| t | Δ_sel (band ± 0.0486) | Δ_uni (band ± 0.0462) | R_boost |
|---|---:|---:|---:|
| 1.0 | **+0.228** [+0.155, +0.280], 4.7× band | **+0.181** [+0.113, +0.240], 3.9× | 3.58 |
| 2.0 | **+0.091** [+0.042, +0.124], 1.9× | **+0.137** [+0.082, +0.168], 3.0× | 1.85 |
| 3.0 | +0.019, within | +0.007, within | 1.18 |
| 4.0 | +0.032, within | +0.034, within | 1.39 |
| 5.0 | −0.047, within (marginal) | −0.046, at the band edge | **0.52** [0.40, 0.70] |

Context: R_nosink = 0.73-1.03 at all t (no selective advantage without a sink
on the uniform line, as both readings expect there); the no-sink arms sit at
Sum-MI ≈ 0.08-0.12, well above the simulated junk floor 0.034 (real-device
junk seeds more MI than the thermal noise model); per-pair at t = 3 the
largest sink-created MI is again at the pair FARTHEST from the sink (0.085 at
(3,4) vs 0.012-0.021 elsewhere), the transport signature, now on hardware.

**Verdict, by the pre-registered rules: INCONCLUSIVE** (the exhaustiveness
clause). B requires both Δ legs beyond band at ≥ 3 of 5 points: both fire at
exactly 2 of 5 (t = 1, 2). A requires both legs within bands everywhere:
t = 1, 2 violate it by 2-5×. The pattern is structured, not noisy:

- **The mechanism IS live at early depth:** the engineered sink CREATES
  interior correlations in both layouts at 2-5× the null band, at close to
  the simulated size (+0.228 measured vs +0.30 simulated at t = 1), with the
  transport signature to the far pair. Reading A ("injected noise only
  removes signal") is REFUTED as a universal statement at these depths.
- **The effect inverts at depth:** by t = 5 the sink configs sit at HALF the
  no-sink Sum-MI (R_boost 0.52, CI excluding 1). The gate-error-free
  simulation predicted monotone-positive Δ through t = 5; the crossover is
  real-hardware physics the sim did not carry (cumulative injected-phase
  damage compounding with gate/crosstalk noise at depth).
- **Post-hoc reading, labeled as such:** creation wins while the sink's
  injected disorder is still localized; at depth its cumulative damage
  overtakes the creation. This is a dose/depth crossover, exactly the shape
  the pre-registered (unflown) partial-strength dose point exists to map.

**Named follow-up (per the INCONCLUSIVE rule):** the dose/depth-resolved
version: the σ ≈ 0.9 rad partial-strength point plus one or two finer time
points around the crossover (t ∈ [2, 5]), same instrument, same bands
recipe. B's sharp signature there: the early-t boost grows with dose while
the crossover point moves earlier; a noise-is-harm frame has no regime with
a positive leg at all, and that regime is now measured.

## Traps carried in (campaign ledger + this review)

Stale calibration (same-day `backend.properties()` selection); no error bars in
March (bootstrap + counts-level bands mandatory here); MI bias floor ≈ the March
uniform-DD row (state it, draw it); drift between configs (one interleaved job);
the sink must be a channel, not shared classical randomness (K ≥ 16 structured,
frozen, basis-shared); the transpiler cannot cancel a virtual RZ blocked by
RXX/RYY (verified in review) but configs must share one transpilation to make the
pair identity checkable; Trotter content is the OBJECT here (DD-during-circuits,
not idle storage); single line, single day, fixed sink end (mirrored-end repeat
pre-registered as follow-up).
