# IBM Concentrator Reloaded: the Site-Resolved Rate Flies (Pre-Registration Draft v3)

<!-- Keywords: concentrator protection metric hardware, site-resolved Absorption Theorem test,
interior coherence lifetime, per-binding randomized Z-sink channel, dose-scaled phases,
weak-coupling regime, ibm_kingston uniform chain, edge vs payload dose contrast -->

**Status:** PRE-REGISTRATION (committed 2026-07-11, before any shot).
Nothing has flown. The design went through five empty review rounds (v1:
three lenses; v2: two; v3: two; v3.1 convergence; all with from-below
Lindblad checks; the revision notes record every finding and its fix).
Stage 7a (from-below predictions + the observable search) and stage 7b
(counts-level gate + certificates) are recorded below; the runner
`run_concentrator_reloaded.py` (external tomography pipeline) was built,
went through its own pre-commit reviews (which found and fixed a missing
hardware counts-persistence path, among others; the 7b block records
them), and its gates are green; the billing question is resolved by
measurement (section 7); an independent C# cross-check (MirrorWorld) is
recorded below. The one remaining step before the flight is the day-of
calibration pull + chain selection with the hard-abort armed.
**Date:** 2026-07-11 (design day)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Hardware target:** ibm_kingston (Heron r2), one uniform-good-T2 5-qubit line
**Relation:** the repeat the March run owes, as its own document.
[IBM concentrator March run](IBM_CONCENTRATOR.md) (ordering-only, mechanism
open) and [the A-vs-B mechanism test](CONCENTRATOR_AB_MECHANISM_TEST.md)
(reckoning: transport observable artifact-dominated, protection metric never
computed) are the two parents; this run is the S2-specific un-park
prerequisite named in the `outbound_label_adapters` ledger entry (necessary,
not sufficient: the general un-park additionally needs a real outreach
trigger, a clean empty review, and citation checks).

---

## 1. The question

The Absorption Theorem's vector form ([the
proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), Theorem 2) prices a
coherence site by site: Re λ = −2 Σ_l γ_l · light_l(v). This run tests, on
hardware, the sharpest flyable consequence: **an engineered dephasing dose
on the qubit that HOLDS a coherence prices it at the theorem's rate, while
the identical dose on a far site prices it at only the small,
coupling-mediated leakage rate.** That is the site-resolved pricing,
measured as a lifetime, for the first time in the arc.

Honest framing (v2 physics review): the distance-GRADED protection reading
(a near sink costs more than a far sink) turned out to be physically
ill-posed for the payload-coherence observable (non-monotone in sink
distance: the excitation sloshes both ways). It is not flown blind; stage
7a runs an observable search for it (zero QPU), and it returns as an arm
only if a resolvable, sign-stable observable exists. This flight's
confirmatory content is the on-payload vs far-edge contrast.

Reconciliation with the un-park banner (v3 audit): the adapter banner
demands "an engineered Z-sink at the formula dose ... plus a
protection/lifetime metric". The reckoning's own Downgrade 3 retired "the
formula dose" (no optimum on the additive axis) and Downgrade 2 retired the
transport-as-protection reading; this run flies the exact site-resolved
pricing at a resolvability-chosen dose instead, with the protection GRADING
carried by the 7a observable search. That is the banner's demand as
corrected by its own later source; the "necessary, not sufficient" hedge
above covers the remainder.

## 2. What the theory predicts (from below), and the regime choice

**The exact anchor (Tier 1, cited as theorem, not flown).** At J = 0 a
coherence held by site 2 collects exactly zero rate from a sink on any
other site (light_l = 0 there; Theorem 2). With coupling, a remote sink
costs a small coupling-mediated leakage instead of zero.

**The regime choice (v2 physics BLOCKER, folded in):** the site-resolved
contrast exists only while the coherence stays localized relative to the
sink. The v2 pin (J·dt = 0.5) delocalizes the payload in ONE Trotter step
(coh 1.0 → 0.29; from-below, exact 32×32 sim of the flown gates), collapsing
the contrast to ≤ 1.4σ and putting the deepest grid point on a transport
node where the pre-registered verdict would have MISFIRED with the wrong
sign. v3 therefore flies the weak-coupling regime: **J·dt = 0.05 per Trotter
step, sink dose UNCHANGED per step**. From-below scouting (v2 + v3 reviews,
superseded by 7a): MP−E ≈ −0.073 per step; projected significance
binding-only ≈ 8.4σ, binding+shot ≈ 6.9σ (counts-level MC at M = 256 × 32);
minimum coh_0 across the grid 0.456, no node.

**What the measured number IS (v3 physics review, folded in):** the payload
is a transported wavepacket, not a Liouvillian eigenmode, so slope(MP) is a
transport-DRESSED effective rate: the static theorem rate −2γ = −0.05 per
step plus a transport-dressing excess (+47% at these settings; scouting
−0.0735/step total, superseded by 7a; the excess makes the effective rate
FASTER, and it leans on the two deepest grid points: without depths 6 and 8
the contrast is −0.054 and the excess only +8%, so 7a must confirm the deep
points survive the device background; A-sign is robust to their loss, the
power and the +47% are not). The theorem content this run tests
is the site-resolved CONTRAST (the on-payload dose prices the coherence,
the far dose prices ≈ 0); the magnitude is recorded as the effective rate
against the 7a prediction, never as the bare −2γ.

**Novelty, stated honestly:** weak coupling sits between the proven J = 0
anchor and the strong-transport regime. What hardware adds: the site
contrast tested with a REAL engineered channel on a REAL device, against
device background. The far-sink leakage is 4th-order tiny here
(scouting slope(E) ≈ −0.0003/step, far below every noise term), so its
verdict is a null-consistency readout, not an expected detection. The kindred repo anchor stays correctly typed: the
concentrator/uniform Q_max ratio 4.9× at N = 5 ([concentrator
optics](CONCENTRATOR_OPTICS.md), Result 1) is a matched-budget slow-mode
quality factor, a different quantity; it motivates, it does not predict.

## 3. Instruments: three arms, one payload, one interleaved job

All arms share one transpiled circuit skeleton, one chain, one job,
round-robin interleaved (nothing split by config; the strongest arm must
not run first).

**Payload and observable:** |+⟩ prepared on position 2 (center of the
5-line), all other qubits |0⟩; 2-basis Ramsey tomography of qubit 2 (X and
Y). The pooled coherence per depth is |mean(⟨X₂⟩ + i⟨Y₂⟩)| (the correct
pooled channel-coherence estimator; its Rician bias is negligible away from
nodes, and the grid has none).

**Estimator (pinned):** the matched-depth paired ratio R_a(t) =
coh_a(t)/coh_0(t); the primary number per sink arm is the fitted slope of
ln R_a(t) over the depth grid. The verdict statistic is the slope
DIFFERENCE between the two sink arms.

- **Arm 0 (baseline):** Trotter chain at J·dt = 0.05, no injected sink.
- **Arm E (far sink):** same circuits + the randomized Z-sink on position 0
  (distance 2 from the payload) at the pinned per-step dose.
- **Arm MP (on-payload sink):** the same per-step dose on position 2. Its
  slope is the transport-dressed effective rate (section 2); it is CHECKED
  against the 7a prediction, never against the bare −2γ (a −0.05 check
  would flag a false +47% mismatch). The injected dose itself is pinned
  independently by the retention table; no circularity. The runner asserts
  exactly ONE sink RZ per sink site per step (a double insertion doubles
  the effect; from-below verified).

**Depth grid:** {1, 2, 3, 4, 6, 8} Trotter steps at J·dt = 0.05. Floor
gate (v2 physics BLOCKER 2 fix): coh_0(t) ≥ 0.3 · coh(t = 0) at EVERY grid
point, referenced to t = 0, checked in 7a and again at 7b; any failing
point is dropped from the grid BEFORE pre-registration, never after data.

**Cut from this flight:** the distance-graded arm (v2's M1; ill-posed for
this observable) and the created-MI transport block (cannot carry a
verdict). Both are named follow-ups; M1 returns only through the 7a
observable search + a design amendment with its own review round.

## 4. The sink: per-binding randomized, dose-scaled, certified

SamplerV2 cannot randomize per shot; the honest channel is **M = 256 fresh
i.i.d. phase vectors per (arm, depth)** (M = 256 FINAL: billing measured
per-PUB, section 7; no re-gate needed), shots split 8192/M = 32 per
binding, counts pooled. Phases are independent ACROSS TROTTER STEPS within each vector and
across bindings (the reckoning's trap was frozen intra-path correlation,
not finite K).

**The dose pin:** per step, φ ~ N(0, σ²) with σ = √(4·d·Δt_step) where the
per-step Lindblad retention is e^{−2·d·Δt_step} (the factor-2 that broke
A-vs-B runs 1-2, verified from below in both v2 reviews). The pinned
per-step retention is e^{−0.05} ≈ 0.9512 (σ ≈ 0.3162 rad), UNCHANGED from
v2. Honest bookkeeping (v2 physics MAJOR): with J·dt lowered to 0.05, this
per-step dose corresponds to ~10× the γ₀ carrier rate per unit J-time; the
dose axis is per-STEP, chosen for resolvability of the site contrast, and
is not claimed as the γ₀-per-J-time carrier operating point.

**Channel certificates (both hard pre-flight aborts):** the
`retention_table` (r = |mean_k e^{iφ_k}|, the realized channel strength;
this is what pins the dose, independently of MP) and the
`ab_classical_mixing_check` artifact-fraction decomposition reproduced on
THIS construction with target < 10% on the payload observable.

**Vacuity bookkeeping:** the sink is a virtual RZ, physical only through
entangling gates; arms are duration-identical by construction; a final-step
sink RZ is a Z-basis no-op and X/Y-physical, which is consistent with a
dephasing channel and expected by the runner.

## 5. Statistics (pinned before the runner exists)

- **Hierarchical bootstrap:** resample the M bindings with replacement,
  then shots within bindings, AND the Arm-0 shots (the denominator leg of
  every R_a); ≥ 500 resamples. This is the only error bar.
- **Null model (v2 audit fix):** H0 = slope(MP) − slope(E) = 0, realized as
  zero-injection synthetic counts for BOTH sink arms through the identical
  pooled pipeline (same M, same shots/binding, same estimator, same
  readout-confusion model, representative device background); ≥ 200 seeds;
  bands frozen as constants with percentile provenance.
- **Roles (v2 audit fix):** the matched-pool NULL band is the significance
  gate; the hierarchical BOOTSTRAP is the measurement error bar; verdict A
  requires both (below).
- **Band provenance (v2 audit MAJOR 1):** bands are frozen at
  pre-registration on a REPRESENTATIVE Kingston background (the freshest
  calibration CSV at freeze time), with a recorded first-order argument for
  chain-independence; the day-of calibration then only feeds the chain-rule
  HARD abort, never re-tunes bands (the price-pair precedent).
- **Confirmatory bar (v2 audit MAJOR 2c):** one-sided 99.87% (3σ-level)
  bootstrap percentile, AND outside the frozen null band; the 7a band then
  gates only the MAGNITUDE verdict A-mag, not the sign verdict A-sign
  (v3.1 review: a correct-sign contrast off the predicted magnitude reads
  "confirmed, magnitude off-prediction", never INCONCLUSIVE).
  Drift caveat carried verbatim: the bootstrap covers shot and binding
  noise only; between-arm drift inside the one job is bounded by the
  interleaving, not by the error bars.
- **Readout mitigation:** CAL0/CAL1 tensor confusion inversion applied to
  the POOLED counts; 7b checks the alternative ordering shifts the slope by
  < 10% of its SE, else escalate to a declared instrument condition.
- **Verdicts:**
  - **A-sign (the site-resolved contrast, the primary claim):** slope(MP)
    − slope(E) < 0 at the 99.87% level and outside the null band → the
    site-resolved pricing CONTRAST is confirmed on hardware (the
    QUALITATIVE site-resolution is the sign; the theorem's QUANTITATIVE
    rate is A-mag; the sign verdict must not be hostage to the magnitude
    landing in a ±2SE window). The null band is pool-structure-matched,
    not signal-level-matched (zero-injection synthetics sit at Arm-0
    coherence), so the low-coherence MP noise is carried by the bootstrap
    leg of the AND, deliberately.
  - **A-mag (magnitude agreement, recorded separately; evaluated only when
    A-sign holds):** the measured
    slope(MP) − slope(E) additionally inside the OPERATIVE band, which is
    the 7b-reconciled band [−0.09167, −0.05482] (it supersedes the 7a
    band; the 7b record explains why) → the
    transport-dressed magnitude matches the from-below prediction; A-sign
    + A-mag together are the Confirmations-candidate (the 7a number as the
    prediction). A-sign without A-mag is recorded as "contrast confirmed,
    magnitude off-prediction" with the deviation quantified, not as
    INCONCLUSIVE.
  - **Anti-A:** slope(MP) − slope(E) significantly POSITIVE at the same
    level and outside the null band → the pricing contrast is falsified
    for this instrument; recorded as such.
  - **L (leakage, null-consistency readout):** slope(E) is expected ≈ 0
    here (scouting −0.0003/step, far below shot, binding, and drift
    noise); slope(E) inside the null band is the expected outcome. A
    significant slope(E) at this regime is treated FIRST as an instrument
    suspicion (drift, crosstalk), investigated before any physics reading;
    it cannot confirm anything on its own.
  - **Instrument conditions:** any failed fit, NaN, guard trip (retention,
    artifact fraction, day-of chain rule, parity floor, mitigation-order
    check) → declared instrument failure, never a physics verdict.
  - **Exhaustiveness:** any other pattern is INCONCLUSIVE, named as such,
    with the follow-up stated.

## 6. The two gates, 7a and 7b (stage names kept from v2; nothing flies before both are recorded here)

**7a. From-below stage (the un-park's "sim first", its own stage):**
exact density-matrix simulation of the FLOWN circuits (Trotter gates at
J·dt = 0.05 + per-step sink channel at the pinned dose + representative
device background), producing predicted paired-ratio slopes and gaps with
bands for all arms. **The 7a band's construction is pinned here (v3 audit):
central slope-difference prediction ± 2× the projected combined SE
(binding + shot, counts-level MC at the FLOWN M and allocation), frozen as
constants with provenance; this band is the A-mag boundary (magnitude
agreement), separate from the A-sign contrast verdict. The OPERATIVE
A-mag band is the 7b-reconciled one, [−0.09167, −0.05482].** Hardware happens only if slope(MP) − slope(E) survives at ≥ 3×
the projected combined (binding + shot) error and every grid point passes
the coh_0 ≥ 0.3·coh(0) floor. PLUS the observable search: is there a
payload/readout combination monotone and resolvable in sink DISTANCE (the
protection grading)? Outcome recorded either way; found → design amendment
+ its own review round; not found → A flies alone and the search result is
part of the record.

**7a RECORDED (2026-07-11; [the 7a
script](../simulations/concentrator_reloaded_7a.py), seeds
20260711/770711/110711; output reproduced exactly on an independent
rerun; promoted from its WIP name at this commit):** representative background = the median of the 56 rule-passing
Kingston qubits from the 2026-07-11T10:02Z calibration (T1 278.1 µs, T2
246.8 µs, readout 0.68%; two rule-passing connected 5-lines exist, best
47-57-67-66-65, whose parameters bracket the medians); per-step wall time
1.2 µs parallel-layered (2.4 µs sequential bound). Central
slope(MP) − slope(E) = **−0.07325/step**; slope(E) = −0.00029/step
(null-consistent, as pre-declared); floor gate PASS at all six grid points
(min coh_0 = 0.438); deep-point sensitivity: without depths 6, 8 the
contrast is −0.05387, same sign (A-sign robust; magnitude and the +47%
dressing lean on the deep points, as pre-flagged). The paired-ratio
contrast is background-independent to 5 digits and wall-time-independent:
the common-mode cancellation (section 8) verified from below. Counts-level
MC (M = 256 × 32, per-qubit asymmetric readout confusion, 500 hierarchical
bootstrap replicates): combined SE = 0.01111, **7a band =
[−0.09546, −0.05104]**, **power margin 6.60 ≥ 3 → hardware-permissible**.
**Observable search: NOT FOUND.** The end payload gives |slope| monotone in
sink distance but wrong-signed for a cost reading (a near sink slightly
RAISES the payload coherence, +0.0056 at distance 1: a Zeno-like leakage
suppression) and unresolvable (gap 0.0057 vs 3·SE = 0.0107); the center
payload confirms the sign flip (+0.0090 at d = 1, −0.0003 at d = 2). Root
tension recorded: the weak coupling that makes the on-payload contrast
clean makes every distance-graded signal 4th-order tiny at this budget.
The distance-graded arm stays cut; Arm A flies alone.

**7b. Counts-level gate:** counts-level simulation of the actual circuits
through the actual estimator (≥ 200 seeds, per-qubit asymmetric confusion,
hierarchical bootstrap), Aer noiseless parity of the flown circuits, the
retention table, the artifact-fraction certificate (< 10%, hard abort),
null bands frozen with provenance, the mitigation-order check, and the
Rician-bias check specifically covering the deepest MP point (the smallest
coherence with the most fit leverage). **Both channel certificates and the
7a band are M-dependent (v3 audit): if billing later forces M down, 7a and
7b re-run at the FLOWN M before commit; certificates computed at a
provisional M never carry.**

**7b RECORDED (2026-07-11, post-fix state after the pre-commit reviews;
runner `run_concentrator_reloaded.py` in the external tomography pipeline,
gate outputs `concentrator_reloaded_gate_{simulate,null,aer,certify}.txt`
beside it):** the two pre-commit reviews (runner code; records audit)
found and fixed: the hardware path originally persisted NO counts and had
no counts→verdict reduction (a paid flight would have been un-analyzable
locally; built now: raw counts per (arm, depth, basis, binding) + CAL +
T1/T2* written into the hardware JSON with fault-tolerant partials, and
`--analyze` runs the full pinned reduction, CAL inversion on pooled counts
→ pooled coherence → paired ratio → hierarchical bootstrap with the Arm-0
leg → section-5 verdicts, instrument conditions first; proven end-to-end
on a hardware-shaped synthetic JSON, which returns real verdicts through
the real path). Also fixed: Aer seeding (now byte-identical across runs;
the earlier record cited a superseded unseeded JSON), the simulate
sampling model reconciled to the hardware's per-depth-independent
bindings, a NaN guard, a post-transpile sink-parameter assertion, the
day-of abort tripping on missing readout, T2*-based sink-end orientation.
Gate results at the fixed state: `--simulate` central −0.07325 (|Δ| vs 7a
< 1e-6), reconciled SE 0.00921, **power margin 7.95**, **A-mag band
[−0.09167, −0.05482]** (this supersedes the 7a band per the pinned
construction rule, ± 2 × projected combined SE at the flown M under the
hardware-faithful sampling model; the shift is printed as a reconciliation
table, nothing silent); 500/500 bootstrap replicates finite. `--null 200`:
200/200 seeds, H0 3σ band [−0.01050, +0.01310]; the predicted contrast
lies far outside → detectable. `--aer` PASS and deterministic (two runs
byte-identical; max |arm0 − clean transport| = 0.01855, max |φ = 0 sink −
arm0| = 0.00000, threshold 0.0442). `--certify` PASS: realized retention
0.9530 vs target 0.9512 (every step within 0.02); classical-mixing
artifact on the payload observable: STRUCTURAL (M → ∞) component exactly
0.000%, systematic Rician 0.58% < 10%; the certificate now PRINTS the full
decomposition itself: single-realization slope gap 0.0208 (28.3% of the
channel slope), seed scatter 12.6%, SE 15.1%, combined
√(scatter² + SE²) = 19.7%, the "double-counting" figure, artifact-backed,
recorded as an interpretive decomposition, not a gate number; the
single-realization gap shrinks with M (0.0208 → 0.0032 → 0.0010 at
M = 256/1024/4096). Runner invariants asserted in code: exactly one sink
RZ per site per step (logical AND post-transpile); one shared skeleton
carrying virtual RZ slots on BOTH candidate sink sites with the inactive
site bound to 0 (identity); identical arm durations. `--hardware` built
and dry-stopped behind `--yes`; `--calibrate` built; neither run.

**Independent C# cross-check (2026-07-11, MirrorWorld, committed with this
pre-registration):** the repo's own continuous-Lindblad engine
(`compute/MirrorWorld`, new additive `siteGammas`/`zz` parameters +
`concentrator N` run mode + 6 from-below tests, suite 124/124 green)
independently gives slope(MP) − slope(E) = **−0.073618/step** at N = 5
(pinned by `ConcentratorTests`, the live witness); the gap to the flown
Trotter value (−0.073249) is −0.000369 = 0.5%, inside the O((J·dt)²)
bound: the discretization error between the flown circuit and its
continuous Lindblad limit, quantified. ZZ moves the
contrast by only 0.15% (the common-mode cancellation, measured). Past the
N = 8 wall: at N = 9 the contrast persists at −0.073941/step and the
far-sink leakage vanishes (distance 4): the site contrast is N-robust
and the far arm's null-consistency reading strengthens with N, as the
Absorption Theorem predicts. The J-convention trap (C# hopping = 2× the
unit-Pauli XY coefficient) was established from the code and pinned by
test, not assumed.

## 7. Cost, budget, and the gate ORDER (one sequence)

Job size: 3 arms × 6 depths × 2 bases = 36 science circuits + CAL0/CAL1 +
in-situ T1/T2* blocks ≈ 48 PUBs (A-vs-B flew 227). Sink PUBs carry M = 256
bindings × 32 shots; total shots ≈ 0.30-0.35M including CAL and T1/T2*
blocks. **BILLING RESOLVED (2026-07-11, measured, zero QPU;
`_query_billing.py` in the pipeline, job.usage() on the July-5 jobs):**
billing is per-PUB / shots-proportional (each A-vs-B job billed exactly
257 s for 908k shots = 0.283 ms/shot, matching the runtime's own 271 s
estimate within 5.5%; per-binding models ruled out by 21×). The reload
projects to **≈ 85-99 s ≈ 1.4-1.7 QPU min**. M = 256 flies unchanged; no
re-gate needed. Side-correction to the budget ledger: the July-5 A-vs-B
campaign billed ~12.8 min total (3 × 4.28 min), well under its ~40-min
estimate, so the 2026 reserve holds more headroom than the earlier
~60-100 min estimate.

**The order (v2 audit MAJORs 4+5 folded in):** design v3 empty review →
fixes → 7a (from-below + observable search) → runner built → 7b gate +
certificates → **empty review of runner + gate records** → **billing
resolved** (retrieve the ACTUAL billed minutes of the July-5 A-vs-B jobs
via the service API, zero QPU; if still ambiguous, one tiny calibration
PUB; per-binding billing → M drops and the statistics re-gate; if that
breaks the 3× power margin, the flight does not happen) → this document
committed as the pre-registration → fresh calibration pull + chain
selection (day-of hard-abort armed) → Tom's go stands (given 2026-07-11) →
ONE job → analysis → RECORD → post-flight empty rounds → the March doc
gets its link.

## 8. Chain rule and guards (carried in verbatim)

The A-vs-B uniform-line rule, unchanged: every qubit T2echo ≥ 150 µs,
max/min T2echo ≤ 2, readout ≤ 2%; HARD abort at submission if the day-of
calibration violates the rule (no override flag). Sink-end orientation:
position 0 = the higher-T2* end. Guards ported from the price-pair runner:
failures[] → one gate; in-situ T1/T2* blocks in the same job; CAL0/CAL1
mitigation; backend.properties() snapshotted before and after; interleaved
everything; 8192 shots per science circuit. **ZZ robustness, correctly based (v3.1 review):** coherent ZZ crosstalk
cancels COMMON-MODE in the paired contrast slope(MP) − slope(E): the sink
RZ commutes with diagonal ZZ, ZZ is identical in both sink arms, and its
residual effect on the transported payload is common to all arms (the
eigenmode argument "Hamiltonian → Im λ only" does NOT apply to a
transported wavepacket, so the contrast, not the theorem, carries the
immunity). The device family's ZZ is a measured known (≈ −3.9/−3.6 kHz
nearest-neighbour, null next-nearest; the price-pair campaign, runs 3-4). Traps carried in from the AB ledger: stale calibration, bias floor
drawn, drift between configs, channel-not-frozen, one shared transpilation,
single line single day, fixed sink end (mirrored-end repeat as follow-up).

## 9. Scope fence

This flight: N = 5, one uniform Kingston line, one day, three arms, one
payload, site-resolved rate confirmatory (A) + leakage secondary (L).
OUT of scope, named follow-ups: the distance-graded protection arm
(pending the 7a observable search), the created-MI transport re-flight
with the randomized sink, N = 7, multiple chains, other-day
reproducibility, the budget-conserving F9-profile comparison, the
mirrored-end repeat.

## RECORD: the flight (2026-07-11, ibm_kingston)

**Flight.** Job `d99a970tcv6s73dn2atg`, chain [109, 108, 107, 106, 105]
(sink end Q109, payload Q107), submitted 22:16 local after a fresh
calibration pull (day-of snapshot in the JSON: T1 = 322/142/237/140/169 µs,
T2 = 265/181/278/176/175 µs, readout 0.55-1.06%; day-of chain rule PASS;
transpiled sink-RZ scheduled duration 0.0 = virtual, PASS). 46 PUBs
(36 science round-robin-interleaved by arm + 2 CAL + 8 in-situ T1/T2*),
376,832 shots, **billed 119 s ≈ 2.0 QPU min** (job.usage() API query,
2026-07-11 22:21 local; the flight log itself carries no usage field, a
runner-revision item). Counts persisted (pooled per cell, see the
instrument deviation below):
`results/reloaded_hardware_20260711_221840.json` (external pipeline).

**Reduction (the real counts path).** CAL confusion (payload) p01 = 0.0088,
p10 = 0.0092; pooled coherences per depth:
coh_0 = [0.930, 0.886, 0.853, 0.788, 0.695, 0.599],
coh_E = [0.920, 0.898, 0.821, 0.772, 0.675, 0.572],
coh_MP = [0.890, 0.796, 0.718, 0.663, 0.503, 0.374].
Bootstrap 500/500 finite, SE = 0.00436. Analysis JSON:
`results/reloaded_analyze_20260711_221936.json`.

**Instrument deviation, disclosed (post-flight audit, its own catch):** the
hardware JSON persisted POOLED per-cell histograms; SamplerV2's
get_counts() aggregated the M = 256 bindings at persistence, so the
pre-registered binding-resample leg was a NO-OP on this artifact and the
realized SE 0.00436 is SHOT-ONLY. The channel itself WAS randomized
correctly (the parameter sweep was built and bound at (256, n_params); the
pooling happened only at save time), so every point estimate is
unaffected. Honest significance against the pre-registered binding+shot
projection (SE 0.00921, the 7b gate): A-sign at ≈ 5.8σ, still far beyond
the 99.87% bar; the shot-only reading would say ~12σ and is not used. In
projected-SE units the A-mag marginality is 0.16 SE above the operative
band edge (even more marginal than the shot-only 0.33). The
mitigation-order "0.00% shift PASS" is VACUOUS with one binding row
(pinned and alternative orders identical by construction) and is withdrawn
as a check; readout-order sensitivity is unquantified on this artifact.
Follow-up pinned for any future flight: persist per-binding counts (bypass
the get_counts() aggregation). Related: the persisted analyze JSON stores
the runner's own stale-band verdict ("A-sign + A-mag CONFIRMED",
a_mag: true); this RECORD, applying the committed operative band and the
projected SE, governs.

**Measured:** slope(MP) − slope(E) = **−0.05337/step**; one-sided 99.87%
bootstrap CI [−0.06505, −0.04107] (entirely negative); null band
[−0.01050, +0.01310], measured value far outside. slope(E) = −0.00585,
inside its null band [−0.01235, +0.00829].

**Verdicts, per the committed rules:**
- **A-sign: CONFIRMED.** The site-resolved pricing contrast is real on
  hardware at the 99.87% level (≈ 5.8σ against the pre-registered
  projected SE, see the instrument deviation) and outside the null band.
  The theorem's qualitative content, the first lifetime-rate contrast of
  the arc; foregrounded honestly: qualitative site-resolution is also what
  any local-dephasing model predicts, the quantitative discriminator was
  A-mag.
- **A-mag: off-prediction, marginal.** Measured −0.05337 vs the OPERATIVE
  7b band [−0.09167, −0.05482]: 0.00145 above the upper edge (0.33 of the
  measured SE). Per the pre-registered rule this reads "contrast
  confirmed, magnitude off-prediction", NOT inconclusive. Instrument note,
  recorded honestly: the runner's `--analyze` printout compared against
  the superseded 7a band ([−0.09546, −0.05104], inside which the value
  falls) and printed "A-mag CONFIRMED"; the committed pre-registration's
  operative band governs this record. The stale band constant in the
  analyze printout is a labeling bug, not a data issue.
- **L: null-consistent**, as pre-declared. No leakage detection claimed.

**The physics reading (a reading, labeled; tightened by the post-flight
referee):** the measured contrast is the dressed prediction compressed to
≈ 0.73× (−0.0534 vs −0.0733; 4.6 shot-only SE = 2.2 projected-SE below
it), and the compression is
UNIFORM across depth slices (matched-grid: measured no-deep slope −0.0369
vs ideal no-deep −0.0539, the same ≈ 0.68-0.73×; the deep points still
steepen the slope on hardware, +45%, as in the ideal sim, +36%; the
earlier draft's "sits on the no-deep prediction" was a cross-grid
coincidence and is withdrawn). At this SNR, "the bare theorem rate
survives" and "generic gate-error compression of the dressed prediction"
are numerically DEGENERATE (bare −2γ lies 0.8 shot-only SE = 0.37
projected-SE from the measured value precisely because 0.68× of the
dressed prediction happens to land there);
the data cannot separate them, and the pre-registration's own §2/§3 forbid
the bare-−2γ yardstick. No theorem claim is made from the magnitude.
Honest foregrounding: the CONFIRMED content (A-sign) is the qualitative
site-resolution, which any local-dephasing model also predicts; the
theorem's quantitative discriminator was A-mag, and it missed, marginally.

**Systematics caveat, carried where it bites (§5, verbatim rule):** on
THIS artifact the bootstrap covers shot noise only (the instrument
deviation above); binding noise enters via the pre-registered projection,
and between-arm drift inside the one job is bounded by the interleaving,
not by the error bars (the before/after properties snapshots returned
identical calibration VALUES, timestamps aside, so they bound nothing
within-job). A-sign at ≈ 5.8σ is safe against it; the A-mag marginality
(0.16-0.33 SE depending on the SE used) and any magnitude interpretation
sit INSIDE that unquantified scale. A structural strength, stated: arms E and MP measure
the SAME payload qubit Q107, so Q107's readout confusion and T1/T2 are
common-mode and cancel in the verdict statistic; a static per-qubit
asymmetry cannot fake the contrast. Scope note: this payload-self-coherence
contrast is a different object from the still-owed "interior coherence
lifetime under the full concentrator PROFILE" instrument (the AB reckoning's
Downgrade 2); that one remains open. Billed 119 s vs the 85-99 s
projection: the in-situ T1/T2* delay blocks run slower per shot, same
pattern as the heating-leg jobs; transparent, immaterial to budget.
(Label fix from the audit: the CAL figures 0.0088/0.0092 are the exact
magnitudes with p01/p10 read in the runner's own convention.)

**Confirmations question (open for the post-flight round):** the
pre-registration made A-sign + A-mag together the Confirmations-candidate;
with A-mag marginally off, the entry decision (register the contrast with
the honest band story, or hold) goes to the post-flight empty rounds and
Tom.

## Revision notes

**v1 → v2 (2026-07-11).** Three empty design reviews on v1 (physics with
from-below Lindblad checks, requirements audit, measurement statistics;
all FIX-THEN-BUILD). Folded in: sink distribution pinned dose-scaled (v1
had U[0,2π) = the ceiling, ~20× the stated dose); estimator pinned to the
matched-depth paired ratio (per-arm envelope fits sign-unstable under
coherent transport); hierarchical binding+shot bootstrap + matched-pool
null bands (shots-only bootstrap understates SE ~3×); mid sink moved off
the payload and the on-payload arm named as the theorem's direct rate
(near-tautology + arm-asymmetric bias); Bell pair → single |+⟩ with
2-basis readout (concurrence unfittable, 9-basis wasteful); M 64 → 256;
transport block dropped; the from-below stage made its own gate and the
unsourced 4.87× anchor replaced by the correctly-typed OPTICS Q_max 4.9×;
channel certificate added; null model, tail rule, mitigation order, time
axis, gate order pinned. Plus a self-caught factor-2 in the v2 σ-pin
(σ = √(4·d·Δt), not √(2·d·Δt)).

**v2 → v3 (2026-07-11).** Two empty design reviews on v2 (physics with
exact from-below sim of the flown construction; requirements+statistics
audit; both FIX-THEN-BUILD, no blockers in the audit, two physics
blockers). Folded in: the REGIME (J·dt 0.5 → 0.05 with the per-step dose
kept; at 0.5 the payload delocalizes in one step, the contrast collapses
to ≤ 1.4σ, and the pinned grid's deepest point sat on a transport node
that would have MISFIRED the pre-registered Anti-A with the wrong sign);
the floor gate re-referenced to coh(t = 0) at every grid point; the
distance-graded arm CUT as physically ill-posed for this observable
(non-monotone in sink distance) and moved to the 7a observable search; the
honest ~10×-carrier dose bookkeeping; H0 restated against the sink-arm
difference; null-band vs bootstrap roles separated; confirmatory bar
raised to 99.87% one-sided + drift caveat; MP-vs-retention circularity
resolved (retention pins the dose, MP is checked against it); band
provenance pinned to a representative background + day-of hard-abort only;
the pre-commit empty review of runner + gate records and the
billing-before-commit step inserted into the order; Arm-0 bootstrap leg,
per-step phase independence, hard-abort certificates, ZZ-immunity
strength, the named leakage verdict L, M-provisional and Rician-bias notes
added.

**v3 fix round (2026-07-11, same day).** Two empty design reviews on v3
(physics with exact from-below + counts-level MC; requirements+statistics
audit; both FIX-THEN-BUILD, no blockers). Folded in: the headline rate
correctly typed as transport-DRESSED (−0.05 static + 47% transport excess;
claiming it as "the theorem's rate" would have been the arc's label error
again, and a runner check against bare −2γ would have flagged a false
mismatch); projected significance restated with shot noise (≈ 6.9σ
binding+shot, 8.4σ binding-only at v3 settings; the v2 scouting figure had
been 8.9σ binding-only); verdict L demoted to a null-consistency readout
(slope(E) ≈ −0.0003/step, unreachable; a firing L is instrument-suspect
first); the 7a band construction pinned (± 2× projected combined SE at the
flown M); certificates and bands re-run at the flown M if billing moves it;
cost stated billing-conditional; the banner reconciliation added; exactly
one sink RZ per site per step asserted; the deep-MP Rician check added to
7b; Anti-A wording tightened.

**v3.1 convergence round (2026-07-11).** One combined empty review on the
post-fix state (all from-below numbers reproduced; FIX-THEN-BUILD, no
blockers). Folded in: the ZZ-robustness argument rebased on common-mode
cancellation in the paired contrast (the eigenmode "Im λ only" argument
does not apply to a transported wavepacket); verdict A split into A-sign
(the contrast, the primary claim) and A-mag (magnitude vs the 7a band), so
a correct-sign contrast off the predicted magnitude reads "confirmed,
magnitude off-prediction" instead of INCONCLUSIVE; the σ-figures
reconciled; shots total tightened; these revision notes ordered
chronologically. The state gated by the NEXT review round is this one
(post-v3.1); per the house rule, post-fix states count as new work.

## Related

- [IBM concentrator March run](IBM_CONCENTRATOR.md), the parent whose What
  Remains this run partially closes
- [The A-vs-B mechanism test](CONCENTRATOR_AB_MECHANISM_TEST.md), the
  reckoning that specified this run
- [The absorption theorem proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md),
  the exact anchor
- [Concentrator optics](CONCENTRATOR_OPTICS.md), the kindred matched-budget
  Q_max anchor
- [The parked DD adapter](../docs/outbound/SELECTIVE_DECOUPLING_SELECTION_RULE.md),
  whose S2-specific un-park prerequisite this run is
- [The price-pair pre-registration](PRICE_PAIR_HARDWARE_PREDICTION.md), the
  genre template
