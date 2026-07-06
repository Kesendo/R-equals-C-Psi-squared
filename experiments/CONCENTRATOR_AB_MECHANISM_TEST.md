# The Concentrator A-vs-B Mechanism Test (design v2, three flown runs, and the reckoning)

<!-- Keywords: selective dynamical decoupling mechanism test, gate-error avoidance vs
noise contrast, uniform T2 chain control, engineered sink injected dephasing random RZ,
concentrator Sum-MI ibm hardware, noise-seeded observable Heisenberg fixed point,
A vs B interpretation April test, R=CPsi2 -->

**Status: NEGATIVE RESULT + methods lesson (2026-07-05); arc PARKED.** This was
built to settle A-vs-B (is the concentrator win the noise-contrast mechanism,
Reading B, or gate-cost avoidance, Reading A?). It did NOT: after three runs on
ibm_kingston (ceiling, N·γ₀/2, N·γ₀ injected) and the same evening's
empty-session review cascade, A-vs-B is OPEN again. The lasting value is
negative and methodological: (1) the created-MI observable is 56-96%
classical-mixing artifact of the frozen K=16 sink at the two partial doses (the
trap: a frozen phase-path construction at partial dose is mostly shared
classical randomness, not a dephasing channel; the ceiling construction is
channel-dominated); (2) created MI measures TRANSPORT, not protection, and the
arc's own 139-360× headline is a peak created-MI figure whose "protection" label
was this arc's error (Downgrade 2); (3) the sibling March run's audit caught
that its "natural concentrator" Q85 was ~93% amplitude damping (T₁), outside the
Z-dephasing theorem's scope, so that run does not validate the mechanism either
([IBM Concentrator](IBM_CONCENTRATOR.md)). The one positive survivor ("injected
noise does not strictly remove signal") is close to textbook ENAQT and is
guaranteed by the |+⟩^⊗5 fixed-point theorem. The factor-2 dose correction is
sound. The outbound adapter built on this arc is now PARKED (not outreach-ready,
its banner). See the reckoning section at the end. The design below is v2; v1 went
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

The March 2026 ibm_torino result (selective DD beats uniform DD at all five
time points, 2.0× on average, up to 3.2×) has two
readings. **A (gate-cost):** leaving the terrible qubit (T2 = 5 µs) undecoupled wins
because DD gates on a lost cause only add errors. **B (concentrator):** the win is the
noise CONTRAST between the noisy edge and the protected interior, the hardware face of
the γ_edge = N·γ_base − (N−1)·ε profile. The discriminating experiment, stated in both
source docs: a chain where NO qubit is a natural sink, plus an ENGINEERED sink.

## The structural fact that frames the observable (review finding, verified exactly)

|+⟩^⊗5 is an exact fixed point of every Heisenberg bond gate (the bond unitary is
e^{iθ/2}(cos θ·I − i sin θ·SWAP), and identical product states are SWAP-invariant), so
under UNIFORM noise the state stays an exact product at all times and Sum-MI ≡ 0
(verified numerically zero through 10 Trotter steps in exact 32×32 simulation;
the residual, 10⁻¹⁵-10⁻¹⁴, depends on the eigenvalue-cutoff convention). **Sum-MI in
this experiment is not a signal that noise degrades; it is a signal that only
NON-uniform noise creates.** Two consequences: (i) the paired sink-vs-no-sink
comparison, not any ratio against uniform DD, carries the physics; (ii) the March
uniform-DD row (Sum-MI 0.013-0.048) is plausibly dominated by the MI estimator's
shot-noise bias floor (+0.014-0.028 at 4000 shots, measured through the March
reconstruction pipeline in the review), a retro-insight recorded here and, now that
this experiment has reported, noted in the March doc (2026-07-05, its Results
section).

## Design (v2)

**Chain:** one 5-qubit connected line, uniformly good by same-day API calibration
(`backend.properties()`): min T2echo ≥ 150 µs, max/min T2echo ≤ 2, readout error
≤ 2%; ranked by uniformity = min T2 / max T2, tie-break higher min T2. Because the
loaded variable is T2* (the undecoupled qubit sits at T2*, not T2echo), position 0
(the sink end) is chosen as the edge with the HIGHER estimated T2*, and the residual
natural contrast (min interior T2echo / edge T2*) is recorded as a design parameter.
(Disclosure, post-flight: T2* is ESTIMATED throughout as T2echo/2.5, a fixed
assumed Ramsey ratio, not a device Ramsey measurement; the sink-end choice, the
recorded contrasts, and the sim noise model all inherit that assumption.)
The mirrored-end repeat is pre-registered as the first follow-up.

**Five configurations, same Trotterized Heisenberg circuit as March (|+⟩^5, RXX+RYY+RZZ
per bond per step, Trotter steps at the even counts {2, 4, 6, 8, 10} = five time
points t = 1..5 in J-time, 9-basis pair tomography, Sum-MI over the 4
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
time point: 4 offsets (2 random plus their negations, per the flown
`sink_construction` record) × a per-step-permuted equispaced quartet
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
step; post-flight note: −ln r is the coherence-rate reading of γ, the
Lindblad-carrier mapping is r = e^{−2γ·Δt}, the factor 2 the RUN 3
pre-registration corrects, immaterial at the ceiling where r = 0). The review's exact Trotter sweep shows Sum-MI rises MONOTONICALLY in per-step
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
hardware junk enters both sides alike, and the estimator bias only approximately:
the sink changes the state and the MI bias is state-dependent, so it does not cancel
exactly in the paired difference, leaving a small residual):**

    Δ_sel(t)  = SumMI(selective_dd_sink) − SumMI(selective_dd)
    Δ_uni(t)  = SumMI(uniform_dd_sink)  − SumMI(uniform_dd)
    R_boost(t) = SumMI(selective_dd_sink) / SumMI(selective_dd)

(R_boost is the ratio form of Δ_sel, same paired legs; reported beside the
Δs in every table.)

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
the CIs; the null spread of Δ under a no-sink world sets the thresholds (written
"2σ" here in v2; implemented as the pooled p99 of |Δ|, disclosed in the
simulator-validation bullet, this sentence left unamended until the post-flight
pass).

## Simulator validation (the gate; recorded 2026-07-05, runner `run_ab_test.py`)

The runner lives in the EXTERNAL tomography pipeline (not this repo), beside
`run_price_pair.py` and `run_sacrifice_zone.py`:
`D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\run_ab_test.py`
(path anchor precedent: `experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md`). Built
to this spec the same day, smoke + full modes; it records, counts-level with the
LITERAL frozen sink construction:

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
  pooled band is anti-conservative at early t for three of the four legs
  (Δ_sel at t = 1 is conservative) and conservative at late t,
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
added-X ratio 0.833 ≈ 5/6); the go was taken on the Model-B planning number
with the Model-C tail risk accepted (Tom, "Du kannst starten"). Data:
[`data/ibm_ab_test_july2026/ab_test_hardware_20260705_174007.json`](../data/ibm_ab_test_july2026/ab_test_hardware_20260705_174007.json)
(the run-1 master; the directory also holds the binding null, chain, and
aer-parity artifacts). This RECORD is a
transcription of the runner's output; its empty-session pass ran the same
evening (it caught the factor-2 in the later runs' dose and fed the
reckoning section).

| t | Δ_sel (band ± 0.0486) | Δ_uni (band ± 0.0462) | R_boost |
|---|---:|---:|---:|
| 1.0 | **+0.228** [+0.155, +0.280], 4.7× band | **+0.181** [+0.113, +0.240], 3.9× | 3.58 |
| 2.0 | **+0.091** [+0.042, +0.124], 1.9× | **+0.137** [+0.082, +0.168], 3.0× | 1.85 |
| 3.0 | +0.019, within | +0.007, within | 1.18 |
| 4.0 | +0.032, within | +0.034, within | 1.39 |
| 5.0 | −0.047, within (0.96× band) | −0.046, marginally beyond (1.0003× band) | **0.52** [0.40, 0.70] |

Context: R_nosink = 0.73-1.03 at all t (no selective advantage without a sink
on the uniform line, as both readings expect there); the no-sink arms sit at
Sum-MI ≈ 0.08-0.12, well above the simulated junk floor 0.034 (real-device
junk seeds more MI than the thermal noise model); per-pair at t = 3 the
largest sink-arm pair MI is at the pair FARTHEST from the sink (0.085 at
(3,4) vs 0.005-0.021 elsewhere; per-pair Δ_sel +0.049 there), suggestive of
transport but at a within-band depth only (the firing depths t = 1-2 are
per-pair mixed and create nearest the sink, so this is not a clean transport
signature; see the verdict bullet).

**Verdict, by the pre-registered rules: INCONCLUSIVE** (the exhaustiveness
clause). B requires both Δ legs beyond band at ≥ 3 of 5 points: both fire at
exactly 2 of 5 (t = 1, 2). A requires both legs within bands everywhere:
t = 1, 2 violate it by 1.9-4.7×, and at t = 5 the Δ_uni leg sits marginally
beyond on the negative side (1.0003× band). The pattern is structured, not
noisy:

- **The sink CREATES interior correlations at early depth:** in both layouts
  at 1.9-4.7× the null band, at close to the simulated size at t = 1 (+0.228
  measured vs +0.30 simulated; at t = 2 the measured +0.091 sits 3.4× below
  the simulated +0.31). Per pair the creation pattern at the firing depths
  is mixed and supports no near/far transport story: at t = 1 the pair (1,2)
  +0.096 leads, the far pair (3,4) +0.053 is second, and the sink pair (0,1)
  is smallest (+0.037); at t = 2 (0,1) +0.054 leads; the far-pair maximum
  appears at t = 3-4, where Δ is within band. Reading A ("injected noise
  only removes signal") is REFUTED as a universal statement at these depths.
  (What the created MI is made of is the reckoning's finding; since the
  CEILING construction turns out channel-dominated there, this early-depth
  creation is the campaign's cleanest device fact, and the device-surviving
  fact stands; the mechanism reading does not.)
- **The effect inverts at depth (suggestive only):** by t = 5 the sink
  configs sit at HALF the no-sink Sum-MI (R_boost 0.52, bootstrap CI
  [0.40, 0.70]). Against the pre-registered calibrations this is not a
  band-level finding: Δ_sel(5) is within the binding band (0.96×), Δ_uni(5)
  marginally beyond (1.0003×), and the ratio null (p99 |R_boost − 1| = 4.2,
  p90 = 1.0) is nowhere near cleared by |0.52 − 1| = 0.48; only the
  non-binding per-t nulls (± 0.025/0.028 at t = 5) support it. The
  gate-error-free simulation predicted monotone-positive Δ through t = 5;
  the crossover reading (cumulative injected-phase damage compounding with
  device noise at depth) is post-hoc and carries the same suggestive-only
  weight as the exposure curve it later fed.
- **Post-hoc reading, labeled as such:** creation wins while the sink's
  injected disorder is still localized; at depth its cumulative damage
  overtakes the creation. This is a dose/depth crossover, exactly the shape
  the pre-registered (unflown) partial-strength dose point exists to map.

**Named follow-up (per the INCONCLUSIVE rule), SHARPENED by Tom's catch
(2026-07-05, same evening: "es braucht nur Gamma 0, nicht das
Maschinen-Gamma"):** the flown sink was the scramble CEILING (per-step
retention r = 0, γ_sink → ∞), but the concentrator formula never asks for
the ceiling; its optimum is FINITE, γ_edge = N·γ_base (post-flight note: the
reckoning's Downgrade 3 finds this a category slip; N·γ_base is the formula's
ε → 0 corner value under a conserved budget, and the additive-sink setting has
no formula optimum, so it serves only as a reference dose). Under Reading B an
overdosed sink must do exactly what the hardware showed: create early,
destroy at depth (post-flight note: this sentence reversed the pre-flight
prediction of a monotone rise with no turnover; see the reckoning's
signature post-mortem). The correct B-test flies AT the formula dose, and Tom's
second catch makes the dose MACHINE-FREE: the repo already measured γ₀ on
IBM (the April chain experiment, `data/ibm_chain_gamma0_april2026/`, and the
Confirmation `gamma0_off_the_lever_kingston_may2026`: critical damping
exactly at J = γ₀, the typed `UniversalCarrierClaim.DefaultGammaZero`
= 0.05). In the circuit's own J-units (the Trotter step carries
J·dt = θ/2 = 0.5), the target is a PURE NUMBER: γ_edge = N·γ₀ = 0.25,
γ_step = γ_edge·(J·dt) = 0.125, per-step retention
**r* = e^{−0.125} = 0.8825** (post-flight: this mapping dropped the Lindblad 2
and delivers HALF the stated rate; caught and corrected in the RUN 3
pre-registration, r* = e^{−0.25}); the dose scale S is solved from the frozen
construction's exact retention r(S) = r*. No T2, no τ_step, no machine
gamma anywhere: the earlier T2-matched recipe in this paragraph's first
draft was still the machine-gamma reading and is demoted to a secondary
candidate (useful as the device-relative cross-point on the dose curve). The dosed run's pre-registered
B-signature: the early-t boost PERSISTS to larger t (no crossover inside the
window, or the crossover moves out as the dose drops toward the optimum),
while the ceiling point (now measured) anchors the overdose end of the
curve. Runner change needed: a σ parameter for the phase paths (the
structured-quartet construction generalizes: quartet × scale). Same
instrument otherwise, same bands recipe, ~13 QPU min.

## RUN 2 pre-registration (2026-07-05, before the shot): the γ₀ dose

Flies the sharpened follow-up immediately (Tom: "es ist Dein Run"). Same
instrument, same chain file (day-of rule hard-abort protects), same binding
null bands (the null binds φ⃗ = 0 and is dose-independent). The one change:
**--dose-scale S = 0.270756**, solved so the frozen construction's exact
per-step retention meets the γ₀ anchor: r(S) = 0.882497 = e^{−0.125},
γ_step = N·γ₀·(J·dt) = 5 × 0.05 × 0.5 = 0.125000 exactly, machine-free
(post-flight: this mapping dropped the Lindblad 2, so this run flew HALF the
intended edge rate; caught in the RUN 3 pre-registration)
(anchor: `UniversalCarrierClaim.DefaultGammaZero`, hardware-confirmed in
`gamma0_off_the_lever_kingston_may2026` and `data/ibm_chain_gamma0_april2026/`).
The runner's dose knob was validated: S = 1 reproduces the ceiling paths
bit-exactly; the scaled retention table is recorded per (t, step); the mild
per-t retention variation and the S≠1 wrap-induced loss of exact negation
closure are recorded as informational (the pooled damping magnitude stays
X-parity-invariant, which is what DD conjugation acts on).

**Sim gate at the γ₀ dose (flight chain, seed 20260705):** Δ_sel =
+0.055 / +0.121 / +0.161 / +0.186 / **+0.191** at t = 1..5, monotonically
GROWING through the window (the ceiling run crushed early and died; the dose
accumulates), CIs exclude 0 at all five, clearing the binding band at
t = 2-5 by 2.5-3.9× (t = 1 marginal at 1.1×); Δ_u similar (+0.047 to
+0.243).

**Verdict rule (pre-registered):** B-dose-CONFIRMED (the persistence
signature) if both Δ legs sit beyond the binding band at ≥ 3 of the returned
time points AND neither leg is significantly negative anywhere in the
window; the ceiling run's measured crossover anchors the overdose end of
the dose curve, so this pair of runs maps the curve's two ends. A
noise-is-harm frame has no positive regime at any dose. Anything else:
INCONCLUSIVE (the exhaustiveness clause carries over verbatim).

## RUN 2 RECORD (2026-07-05, ibm_kingston, job d9581isql68s73caav1g): rule fired; verdict downgraded post-flight (see the reckoning)

Same chain [109, 108, 107, 106, 105] (day-of rule re-passed, hard-abort armed),
dose S = 0.270756 (γ_step = 0.125000 = N·γ₀·J·dt exactly), all pre-flight
asserts PASS. Data:
[`data/ibm_ab_test_july2026/ab_test_hardware_20260705_181107.json`](../data/ibm_ab_test_july2026/ab_test_hardware_20260705_181107.json).
This RECORD is a transcription of the runner output; its empty-session pass
(with run 1's) ran the same evening (genre: minimally framed fresh reviewer
agents recomputing from the raw JSONs, per the skill
`becoming-your-own-outside`) and its findings are the RUN 3 pre-registration
(the factor-2) and the reckoning section (the observable).

| t | Δ_sel (band ± 0.0486) | Δ_uni (band ± 0.0462) | R_boost |
|---|---:|---:|---:|
| 1.0 | +0.014, within | +0.035, within | 1.16 |
| 2.0 | +0.028, within | **+0.066** ✓ | 1.29 |
| 3.0 | **+0.059** ✓ | **+0.120** ✓ | 1.61 |
| 4.0 | **+0.095** ✓ | **+0.110** ✓ | 2.16 |
| 5.0 | **+0.071** ✓ | **+0.092** ✓ | 1.86 |

**By the pre-registered rule the pattern read B-dose-CONFIRMED.** Δ_sel beyond
the binding band at 3 of 5 points (t = 3, 4, 5), Δ_uni at 4 of 5 (t = 2-5), and
neither leg significantly negative anywhere in the window: the persistence
signature, exactly the shape the γ₀-dose simulation predicted (measured
maximum +0.095 vs simulated +0.186, ≈ half the size, same monotone-growth
form; the estimator-bias floor and real-device junk plausibly absorb the
difference, stated not analyzed). (Two later downgrades apply: this run flew
HALF the formula dose, the factor-2 caught in the RUN 3 pre-registration; and
the observable itself is largely a construction artifact, the reckoning
section. The rule verdict is withdrawn as a mechanism confirmation.) R_nosink ≈ 0.92-1.31 (no selective
advantage without a sink, third consistent reading); the far-pair maximum
repeats (largest sink-ARM pair MI at (3,4), t = 3: 0.050/0.080 raw; the
created, differenced values there are +0.020/+0.054).

**The two runs together map the dose curve's two ends, as pre-registered,
and the curve's right name is an EXPOSURE curve (Tom's third catch, same
evening: "Gamma zerstört doch nicht; wenn es eines macht, dann vielleicht
überbelichten"):** γ is the light, and it destroys nothing. No sink means
no illumination contrast, and a contrast-free scene records no image (the
fixed-point theorem, photographically). The γ₀ dose is the correct
exposure: the image accumulates monotonically and persists through the
window (run 2, on both layouts). The ceiling dose is
OVERexposure: r = 0 per step is a full re-exposure of the edge every step;
the image forms early and then washes out at depth (run 1's crossover),
saturation, not destruction. "Injected noise only removes signal" is
refuted at both doses; that part survives. The original conclusion here
("the concentrator mechanism is live, engineerable, and exposure-controlled
at the formula's own dose") did NOT survive the same evening's empty passes:
this run flew HALF the formula dose (the factor-2, caught and corrected in
RUN 3), the run-1-vs-run-2 comparison shares a chain but not a job or hour
(and run 3 not even the chain, so the cross-run curve is confounded), and
the created-MI observable is largely a classical-mixing artifact of the
frozen sink and not the protection metric (the reckoning section). The
exposure picture is kept as a suggestive post-hoc reading of that confounded,
artifact-dominated observable, nothing more. The Confirmations-registry
entry planned here was NOT added, for exactly that reason. Scope carried
over from the run-1 rules either way: nothing here retroactively proves
March WAS B.

## RUN 3 pre-registration (2026-07-05, before the shot): the factor-2-corrected dose

An empty-session review of both RECORDs caught a factor of 2 in the γ₀ dose.
Runs 1-2 anchored the sink to γ₀ but mapped it to a per-step retention as
`γ_step = N·γ₀·(J·dt)`, dropping the Lindblad 2: γ₀ is the framework carrier
with Re(λ) = −2γ₀·k (a k=1 edge coherence decays e^{−2γ₀t}; `UniversalCarrierClaim`:
γ₀ = 1/(2·T₂)), so the per-step edge retention is e^{−2·γ_edge·(J·dt)}, not
e^{−γ_edge·(J·dt)}. Run 2 therefore flew γ_edge = N·γ₀/2 = 0.125, **half** the
formula dose. The runner's own `recommend_dose` uses gamma_base = 1/T2 (a
coherence rate = 2γ) and is correct; only the manual γ₀ anchor dropped the 2.
This is the Lindblad coherence-rate factor of 2 (2γ₀ = 1/T₂, the −2γ of the
Absorption Theorem): a coherence rate put in a slot that wants γ₀. The same
confusion sits in the March Formula-Validation retro-note
([`IBM_CONCENTRATOR.md`](IBM_CONCENTRATOR.md)), a genuine twin. It is a
*different* factor of 2 from the q = Q/2 relabeling (a Hamiltonian J
normalization, GLOSSARY, not a rate convention) and from the ≈2 in
`gamma_0_marrakesh_calibration` (a continuous-vs-Trotter discretization bias
absorbed into a fitted γ_Z, only coincidentally near 2 for that dataset). Same
numeral, distinct mechanisms; do not lump them.

**The corrected dose:** γ_edge = N·γ₀ = 0.25 (edge Lindblad rate, the actual
formula dose), per-step retention r* = e^{−2·N·γ₀·(J·dt)} = e^{−0.25} = 0.778801.
Solved through the runner's exact frozen-path construction at the run seed (which
reproduces the Run-2 dose r = 0.882497 bit-for-bit as a control):
**`--dose-scale 0.378105`** (target r* = 0.778801; achieved mean per-step
r = 0.778800, γ_step = 0.250001).
Same instrument, same day-of chain rule (hard-abort armed), same binding null
bands (φ⃗ = 0, dose-independent: Δ_sel ± 0.0486, Δ_uni ± 0.0462).

**Sim gate at the corrected dose (flight chain, seed 20260705, --quick):** Δ_sel =
+0.150 / +0.238 / +0.302 / +0.266 / +0.217 at t = 1..5, all five CIs excluding 0
and all far beyond the binding band; Δ_uni +0.191 / +0.252 / +0.286 / +0.315 /
+0.240. Δ_sel PEAKS at t = 3 and declines (Δ_uni at t = 4), where the half-dose
Run-2 sim grew
monotonically to t = 5: doubling the dose moves the peak forward, so this
dose already sits near the exposure turnover.

**Caveat, pre-shot (the sim omits the device's own rate):** the counts-level sim
carries only the injected sink, not the machine's natural dephasing. On hardware
the effective edge rate is the injected N·γ₀ PLUS the device γ, so the effective
dose is HIGHER than the sim and the peak/turnover is expected EARLIER (and the
late-t washout stronger) than the sim's t = 3. Run 2 already turned over earlier
on hardware than in its own sim, consistent with this. (Post-flight correction:
this caveat's premise was WRONG as stated; the flown JSONs' own `noise_model`
field records "T2echo_vs_T2star (Ramsey ratios)" in the sim gate, i.e. the sim
DID carry a device-noise model, and the gate's Δ ≈ half the noiseless struct-K16
values is consistent with that. What survives of the caveat is only the weaker
under-modeled-junk argument: the device no-sink arms sit at 0.08-0.14 vs the
simulated floor 0.034.)

**Verdict rule (pre-registered, same as Run 2):** B-dose-CONFIRMED if both Δ legs
sit beyond the binding band at ≥ 3 of the returned time points AND neither leg is
significantly negative in the window; else INCONCLUSIVE by the exhaustiveness
clause. Run 3 flies the actual formula-scaled edge rate (N·γ₀), so it becomes the
"formula dose" point of the exposure curve; Run 2 stays on record as the half-dose
point, Run 1 as the ceiling.

## RUN 3 RECORD (2026-07-05, ibm_kingston, job d95a33sql68s73cad7fg): rule fires 5/5; verdict downgraded post-flight (see the reckoning)

Flown ~20:25-20:32 local at the corrected dose S = 0.378105 (per-step r =
0.767-0.788, γ_step ≈ 0.24-0.27, mean 0.250 = 2·N·γ₀·(J·dt), the ×2 fix), all pre-flight
asserts PASS (day-of rule, binding identity, added-X ratio 0.833 ≈ 5/6, 270
parameterized RZ survive transpile). Fresh day-of chain [50, 51, 58, 71, 72]
(uniformity 0.617, min T2echo 171 µs, sink edge Q50, residual natural contrast
1.54×), the runner's best pick today; Runs 1-2 flew [109..105], so Run 3 is a
different (better-uniformity) line. Data:
[`ab_test_hardware_20260705_203200.json`](../data/ibm_ab_test_july2026/ab_test_hardware_20260705_203200.json),
chain [`ab_test_chain_run3_20260705.json`](../data/ibm_ab_test_july2026/ab_test_chain_run3_20260705.json).
Bands: the Run-2 binding null (φ⃗ = 0, dose-independent) reused, ± 0.0486 /
± 0.0462; it was calibrated on [109..105] and is an estimator-noise spread
(chain-independent to first order), applied here as pre-registered.

| t | Δ_sel (band ± 0.0486) | Δ_uni (band ± 0.0462) | R_boost |
|---|---:|---:|---:|
| 1.0 | **+0.094** [+0.031, +0.151] ✓ | **+0.077** [+0.007, +0.123] ✓ | 1.94 |
| 2.0 | **+0.070** [+0.018, +0.108] ✓ | **+0.128** [+0.077, +0.161] ✓ | 1.60 |
| 3.0 | **+0.088** [+0.039, +0.122] ✓ | **+0.094** [+0.051, +0.120] ✓ | 1.79 |
| 4.0 | **+0.112** [+0.073, +0.145] ✓ | **+0.089** [+0.048, +0.122] ✓ | 2.27 |
| 5.0 | **+0.073** [+0.035, +0.111] ✓ | **+0.080** [+0.043, +0.109] ✓ | 1.75 |

**By the pre-registered rule the pattern reads B-dose-CONFIRMED** (both Δ legs
beyond the binding band at 5 of 5 time points, no significantly-negative leg),
cleaner than the half-dose Run 2 (3/5 and 4/5). **The same evening's
empty-session passes then showed that the rule's OBSERVABLE does not measure the
mechanism the rule was written to confirm, so that verdict is withdrawn; the
full accounting is the reckoning section below.** What the flight itself
establishes:

- **The ×2 dose fix is sound.** Run 3 injected the γ₀-anchored rate
  γ_edge = N·γ₀ as mapped (the fix verified three independent ways plus the
  runner code; the frozen-path construction reproduces the Run-2 dose
  bit-for-bit as a control). Whether N·γ₀ deserves the name "the formula
  dose" is a separate question the reckoning answers no to (Downgrade 3).
- **The positive Δ is beyond shot-sampling noise, thinly at the closest
  points.** Both legs clear the pre-registered sim band; a stricter 2.33×
  bootstrap-SE threshold on the SAME bootstrap distribution (a
  re-parameterization, not an independent hurdle) is also cleared at all 5
  points, but by as little as Δ_uni 1.11× at t = 1 and Δ_sel 1.31× at t = 2,
  and that null is shot-noise-only, blind to drift and config-level
  systematics. What that Δ is made of is the reckoning's finding.
- **Persistent, no crossover to negative.** Both legs stay beyond band across
  the whole window; Δ_sel peaks at t = 4 and declines ~35% to t = 5 (the
  CIs at t = 3-5 overlap, so the peak position is indistinguishable within
  that range), unlike Run 1's ceiling washout to negative. Noted honestly:
  the pre-shot caveat expected the peak EARLIER than the sim's t = 3; within
  resolution it came at t = 4 or later, not earlier, unexplained (and the
  caveat's premise was itself wrong, see the post-flight correction above).
- **Amplitude modest.** Measured Δ_sel 0.07-0.11 (Δ_uni up to 0.13) vs the
  sim gate's 0.15-0.32; the residual compression is attributable only to the
  under-modeled junk floor (the sim gate DID carry a noise model, see the
  corrected caveat); the no-sink arms sit at Sum-MI
  0.08-0.14 (device junk floor) with the sink arms ~1.7-2.4× above them
  (R_sink).
- **Far-pair count, stated honestly:** the far pair (3,4) carries the largest
  sink-created pair-MI (per-pair Δ_sel) at 3 of 5 time points (t = 1, 3, 4:
  +0.042 / +0.041 / +0.053; at t = 2 it is (0,1), at t = 5 it is (2,3)), sink
  at Q50. Suggestive of transport, not universal; the earlier "a fourth time"
  phrasing over-counted. R_nosink ≈ 0.86-1.32 (no selective advantage without
  a sink, again).

The three runs sit at three doses (Run 1 = ceiling γ→∞, Run 2 = N·γ₀/2,
Run 3 = N·γ₀), but they differ in chain, job, and time of day, so the
three-point "exposure curve" is CONFOUNDED and only suggestive. And the
reckoning below shows the observable itself is dominated at the partial doses
by a construction artifact, so no dose-curve reading of created MI carries
mechanism weight.

## The reckoning (2026-07-05, late evening): what the created-MI observable measures

The empty-session cascade on the three RECORDs (an external-referee pass, a
record audit recomputing from the raw JSONs, and a zero-QPU reproduction of the
flown circuit) ended with one sound correction, three downgrades, and a
signature post-mortem. Recorded with
the same rigor as the runs; this section is the arc's honest conclusion.

**Sound, kept: the factor-2 dose correction.** The RUN 3 pre-registration's
analysis stands in full: runs 1-2 put the Lindblad carrier γ₀ = 1/(2·T₂) in a
coherence-rate slot, Run 2 flew half the formula dose, and Run 3 corrected it.
The math was checked three independent ways plus the runner code; the runner's
own `recommend_dose` was always correct. This is the exact Lindblad rate factor
(the −2γ eigenvalue, 2γ₀ = 1/T₂ coherence rate, twice the carrier γ₀); the same confusion recurs in the
March Formula-Validation retro-note. It is NOT the q = Q/2 J-relabeling, and NOT
the continuous-vs-Trotter factor of `gamma_0_marrakesh_calibration` (a different
≈2): the repo's factor-of-2 traps are several distinct objects, not one costume.

**Downgrade 1: at the flown partial doses the created MI is 56-96%
classical-mixing artifact of the frozen sink construction; the artifact
fraction is dose-dependent.** Zero-QPU check
([`simulations/ab_classical_mixing_check.py`](../simulations/ab_classical_mixing_check.py),
a noiseless exact density-matrix simulation of the flown circuit at all three
doses; its φ⃗ = 0
gate reproduces Sum-MI 0, the fixed-point theorem): at the Run-3 dose the flown
struct-K16 sink gives Sum-MI 0.32-0.59 (t = 1..5), while a TRUE dephasing
channel at the same per-step retention gives 0.02-0.12; the classical-mixing
floor (struct-K16 minus channel) is 62% / 82% / 95% / 95% / 95% of the signal.
At the Run-2 half dose the floor is 56-96%. At the Run-1 CEILING the same
decomposition gives 3.5% / 27% / 20% / 24% / −9%: the ceiling construction is
channel-DOMINATED (not channel-equivalent: only the per-step marginal is exact
scramble at r = 0, the frozen multi-step correlation survives, hence the
20-27% mid-t floor and the unexplained −9% at t = 5; the v2 design's ceiling
validation checked a different, counts-level quantity, residual bias
+0.004-0.016 vs iid K ≈ 32), so Run 1's early-t creation is the
campaign's cleanest device fact, and the artifact owns the partial doses only.
"Floor = struct − channel" is an additive split of a nonlinear functional
(MI), so the fractions are heuristic, which scope statement (i) below also
concedes. One further consequence, recorded: X-X DD conjugation preserves the
sink's per-step CHANNEL magnitude but reshuffles the frozen intra-path
correlations, so at the partial doses the selective and uniform layouts share
the channel but NOT the artifact, and cross-layout comparisons (Δ_sel vs
Δ_uni) are artifact-confounded there.
The iid control (independent-per-step phases at the same marginal
retention) converges to the true channel as K grows, so the partial-dose
artifact is the frozen intra-path correlation of the 16 shared
phase paths, NOT finite K: more shots or more paths do not fix it. The design's
own trap line ("the sink must be a channel, not shared classical randomness")
fired after all, at the doses its validation never covered. Three scope
statements, for honesty: (i) these fractions are noiseless decompositions of
the flown ensemble; on the device the sink arms measure Sum-MI 0.16-0.23, well
below the noiseless 0.32-0.59, so the split does not transfer linearly and the
device-level fraction is not directly measured; (ii) the measured Δ is above
shot-sampling noise either way (Run 3 clears a hardware-native null at all 5
points, though that null is shot-noise-only and blind to drift/systematics);
(iii) the noiseless channel-only signal clears the binding band only at
t = 1-2 (0.122, 0.073 vs ± 0.0486) and sits below it at t = 3-5
(0.024-0.030), so the late-t persistence that fired the Run-2/3 rule is not
reachable by the intended channel at these doses: the rule-firing shape is
the artifact's.

**Downgrade 2: created interior MI measures transport, and "protection" was
this arc's own label error on it.** Tom's question, verbatim: "how do we know
large MI is better?" We do not, and the question reaches further than these
runs. The 139-360× headline of [Resonant Return](RESONANT_RETURN.md) is a
PEAK CREATED-MI ratio (its "Formula vs Optimizers" table's Peak SumMI column:
0.230 vs 0.000639 at N = 5; |+⟩^⊗N initial state, where the uniform profile
gives Sum-MI ≡ 0); the formula was derived by maximizing the Mode-2 TRANSPORT
projection, and no interior-coherence-lifetime figure was ever computed
anywhere in this arc (the "rings far longer" reading entered in
[Inside/Outside the Sacrifice Zone](../docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md)
and propagated). So created MI is not a proxy for the concentrator's metric;
it IS the concentrator's historical metric, and the "protection / lifetime"
name this experiment (and the outbound adapter) attached to the 139-360× was
a label error, corrected today. What protection content exists is mode-level
and exact (the Absorption Theorem: interior modes with low ⟨n_XY⟩ decay
slowly; the profile's ε-input "protects" by construction), but a perfectly
protecting concentrator KEEPS the interior a product state and gives
Sum-MI = 0, so the transport metric and the protection objective PULL
OPPOSITE WAYS at the extremes, and no computed number connects them.
Consequences: (i) these runs measured the right observable CLASS for the
transport claim, and at the partial doses measured mostly Downgrade 1's
construction artifact of it; (ii) they measured nothing about protection;
(iii) whether channel-contrast-created transport MI is a resource (the ENAQT
reading) or merely detectable structure is untested; (iv) a protection
(lifetime) objective for the concentrator profile has never been computed,
in simulation or on hardware, and is now the named open instrument.

**Downgrade 3: "the formula dose" was a reference value, not the formula's
optimum.** The concentrator formula γ_edge = N·γ_base − (N−1)·ε is the corner
of a BUDGET-CONSTRAINED optimization ([Resonant Return](RESONANT_RETURN.md):
"the total noise across the chain stays the same, only its distribution
changes"), and it prescribes the TOTAL edge rate of a chain whose interior is
protected at ε. The flown runs conserve no budget (the sink ADDS noise to an
interior that stays at the device rate), and the injected N·γ₀ sits on top of
the edge's own device rate, so the physical edge rate was ≈ (N+1)·γ₀; while
the formula, at the natural hardware identification ε ≈ γ_base ≈ γ₀, would
prescribe γ_edge = γ₀, i.e. no injection at all. The dropped −(N−1)·ε term
(≈ 0.20 at N = 5, γ₀ = 0.05) is larger than the factor-2 that was caught and
corrected. So N·γ₀ is the formula's ε → 0 corner value, used here as a
machine-free REFERENCE dose on an axis the formula does not govern; "flew the
formula's own dose" over-read it, and the run-1 follow-up's "its optimum is
FINITE" over-read it first. The factor-2 correction itself is untouched by
this: it fixes the mapping from a stated γ_edge to a per-step retention and
is correct whatever the target.

**The signature was never a discriminator (post-mortem).** Pre-flight, the
design pinned Reading B to a monotone RISE of the paired effect with
injection variance, and the review sweep found no Zeno turnover; after
Run 1's measured crossover, the named follow-up asserted an overdosed sink
"must do exactly" that, and Runs 2-3 were then scored against the revised
persistence signature. A signature revised after the data it is scored on is
not a discriminator; and nothing in the concentrator theory was shown to
predict either the crossover or its absence, so persistence never separated
Reading B from generic injected-dephasing-spreading-through-coupling
dynamics. The A rule had the twin weakness: "both Δ ≤ 0" operationalizes a
strawman (the gate-cost hypothesis is about DD pulse economics on a bad
qubit, not injected virtual phases), and since the ideal sim already produces
the boost, an A verdict could only ever have measured detectability. The v1
review had in fact already named the outcome: "the pooled object at small K
is shared classical randomness, which even Reading A boosts"; the
partial-dose positive Δ confirms that warning, not physics beyond it. (Band
methodology, moot post-withdrawal but recorded: the sim null's junk floor
0.034 sat 3× below the device's 0.08-0.14; the pooled band is
anti-conservative at exactly the early t where Run 1's headline firings sit,
mitigated by the bootstrap CIs; Run 3 reused bands calibrated
on [109..105], asserted chain-independent to first order, cross-checked there
by the hardware-native null, which is itself shot-noise-only and blind to
drift and config-level systematics. Provenance, also recorded: the pre-flight
sim-gate numbers and same-day calibration narratives quoted in the RECORDs
survive only as transcriptions in this document; the external pipeline
archived no simulate artifact for the flight chains, so they are permanently
unverifiable. Everything flown, the bands, and the dose constructions
recompute bit-level from the committed JSONs. One caveat on that: the
classical-mixing check (`ab_classical_mixing_check.py`, Downgrade 1) re-derives
the frozen phase paths through the external `run_ab_test` runner rather than
reading the `sink_phases` block committed in the flight JSONs, so that specific
script is not runnable from a repo clone as written (its inputs are in the
committed data; the reproduction is not yet vendored into the repo).)

**The narrow claim that survives:** an engineered per-step phase sink on one
edge qubit produces a positive, device-surviving increase of nearest-neighbour
pair-MI (summed over the 4 adjacent pairs, so not purely interior, one pair
contains the sink), on both DD layouts: beyond the pre-registered bands at
early depth in Run 1 (1.9-4.7×, on the channel-dominated ceiling construction;
its t = 5 Δ_uni leg turns marginally beyond-band negative, Δ_sel stays within at
0.96×) and at 3-5 of 5 points in Runs 2-3 (with the
hardware-native-null cross-check on record for Run 3, both legs, in the
committed audit script). That refutes "injected noise strictly removes signal"
as a universal statement. The strongest single piece is Run 1's early-depth
creation: a genuine Z-rotation (correct-channel) edge sink seeds interior MI
through the Heisenberg coupling at 4.7× the band, matching the noiseless
prediction in shape, the campaign's one clean non-artifact device signal,
though at the off-regime γ→∞ ceiling and testing neither A/B nor protection.
Beyond that, little: at the partial doses the increase is largely the
construction's classical mixing, and none of it is evidence that the
concentrator mechanism protects anything.

**What stands, restated with the corrected label:** the simulation result
stands AS A TRANSPORT RESULT: under true Lindblad profiles at a fixed budget,
the concentrated profile creates more peak summed nearest-neighbour MI (360× at
N=5, declining with N to ~63× by N=15) than the
V-profile ([Resonant Return](RESONANT_RETURN.md)); Downgrade 1 does not touch
it (the sim used true channels, not the frozen construction). The March
natural-sink run's 5-of-5 Sum-MI ordering stands as recorded, with its
bias-floor and attribution caveats ([IBM Concentrator](IBM_CONCENTRATOR.md)).
The A-vs-B question this experiment was built to answer is OPEN again, and
the arc now owes TWO instruments, named: (1) a protection-metric computation,
interior coherence lifetime under the concentrator profile, first in
simulation from below, then, only if it survives there, on hardware; this
would test the protection reading for the first time anywhere in the arc;
and (2) for any created-MI re-flight, an unfrozen, per-shot-randomized sink
so the injected noise is a channel rather than a shared classical phase.
Neither is pre-registered here.

## Traps carried in (campaign ledger + this review)

Stale calibration (same-day `backend.properties()` selection); no error bars in
March (bootstrap + counts-level bands mandatory here); MI bias floor ≈ the March
uniform-DD row (state it, draw it); drift between configs (one interleaved job);
the sink must be a channel, not shared classical randomness (K ≥ 16 structured,
frozen, basis-shared; post-flight addition: this trap FIRED anyway at the
partial doses, where the frozen construction departs from the channel by
56-96%, see the
reckoning); the transpiler cannot cancel a virtual RZ blocked by
RXX/RYY (verified in review) but configs must share one transpilation to make the
pair identity checkable; Trotter content is the OBJECT here (DD-during-circuits,
not idle storage); single line, single day, fixed sink end (mirrored-end repeat
pre-registered as follow-up).
