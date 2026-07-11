# IBM Hardware: Selective DD Yields 1.4-3.2x Higher Sum-MI than Uniform DD at All Five Time Points (Concentrator Mechanism Open)

<!-- Keywords: IBM Torino hardware validation concentrator profile, selective dynamic
decoupling beats uniform DD, first spatial dephasing profile hardware test,
Q85 natural concentrator qubit T2=5us, Sum-MI selective vs uniform 2x-3.2x,
palindrome-derived noise engineering quantum hardware, R=CPsi2 IBM experiment -->

**Status:** Hardware, Tier 2: the selective > uniform Sum-MI ORDERING holds at
all 5 time points, now with error bars (2026-07-11 re-analysis in Results:
joint bootstrap ordering confidence ≈ 0.93 floor-bias-adjusted, ≈ 0.95 raw).
Single run, one chain. The ratio SIZE (1.4-3.2× per point, avg 2.0×) is
floor-sensitive (floor-correcting raises it; late points soft) and Sum-MI
here is a noise-seeded transport proxy (2026-07-05 retro-note in Results);
the A-vs-B attribution is open.
**Date:** March 24, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Hardware:** ibm_torino (Heron r2, 133 qubits)
**Chain:** Q85-Q86-Q87-Q88-Q94 (concentrator: Q85, T2echo=5.2 us)
**QPU time used:** ~163s of 210s budget (47s remaining)
**Data:** [data/ibm_sacrifice_zone_march2026/](../data/ibm_sacrifice_zone_march2026/)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

**Naming note (2026-07-05):** renamed from "IBM Hardware: Sacrifice-Zone
Selective DD...". The edge qubit (Q85) sacrifices nothing; it concentrates the
noise and lets the interior sit at the fold. "Concentrator" is the corrected
term (the "sacrifice zone" misnomer was resolved 2026-03-28, see [Inside and
Outside the Sacrifice Zone](../docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md), where
"sacrifice" survives as the inside-the-boundary reading). The frozen data
directory keeps its original name `data/ibm_sacrifice_zone_march2026/` as the
provenance of the actual run.

---

## Abstract

First hardware test of the concentrator formula's spatial noise profile,
approximated by qubit-selective Dynamic Decoupling. (Selective DD itself is
established practice, and the doc's own Interpretation A notes an
experimentalist would skip DD on a bad qubit anyway; the formula-derived
edge-concentrator profile is the new part.) On a 5-qubit Heisenberg chain (ibm_torino,
Q85-86-87-88-94), removing DD from one edge qubit (Q85, natural T2=5.2 us)
while protecting the other four produces 2.0x higher Sum-MI (summed
nearest-neighbour mutual information over the four adjacent pairs) on average
and up to 3.2x at t=4.0, compared to standard uniform DD on all qubits.
The selective > uniform Sum-MI ordering holds at all 5 measured time points
("higher Sum-MI" is a transport-signal reading, not a protection win; see
the caveats).

**Caveats:** Single run on one chain. Error bars were added by the 2026-07-11
re-analysis (Results): the 5-of-5 ordering is robust to shot noise, joint
bootstrap ordering confidence ≈ 0.93 floor-bias-adjusted (≈ 0.95 raw). The
advantage may partly reflect that DD on Q85 (T2=5.2 us)
adds gate errors without benefit, rather than a pure concentrator
effect. Reproducibility across chains and days is untested. The connection
to the palindromic eigenstructure is indirect: the formula predicted the
optimal noise profile, but the hardware implementation approximates it
via DD rather than controlling dephasing rates directly. Added 2026-07-05
(retro-note in Results): Sum-MI on |+⟩^5 is a noise-seeded transport proxy,
not a protection metric. The 07-05 note's floor worry is now quantified
(2026-07-11): floor-correcting the ratios RAISES them; the late-time
corrected values are soft.


## Setup

### Chain Selection

ibm_torino calibration (March 24, 2026) identified Q85 as an extreme
outlier: T2(echo) = 5.2 us while neighboring qubits have T2 = 123-244 us.
This provides a natural concentrator qubit with 51.9x effective contrast
(mean T2echo of the four protected qubits / T2*_concentrator).

**Channel caveat (2026-07-05 review, important; numbers recomputed and the
un-echoed slice added 2026-07-11 from the stored chain_info, the 07-05 note
had used rounded inputs and compared only the echo slice):** Q85's short
T2echo is ~92% amplitude damping, not Z-dephasing. With T1 = 2.836 us,
1/(2·T1) = 0.176 us⁻¹ is 92% of 1/T2echo = 0.191 us⁻¹, so Q85's echoed
PURE-dephasing time is Tφ ≈ 66 us, only ~2.6-6.7x shorter than the interior
(Tφ ≈ 170-440 us), not 52x. In the winning SELECTIVE config, though, Q85
runs UN-echoed at the T2* level: Tφ* ≈ 10.9 us (amplitude-damping fraction
~66% there), so the Z-dephasing rate contrast actually present against the
echoed interior is ~16-40x.
The concentrator selection rule and its [Absorption
Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) are proven for
Z-dephasing only (amplitude damping is the separate F82/F84 object, and DD
cannot refocus T1). So this natural-sink run does not cleanly test the
Z-dephasing rule: the T1 sink is the larger channel at both slices (~92%
echoed, ~66% un-echoed), which keeps the gate-cost reading (Interpretation
A, DD-on-a-T1-limited-qubit hurts) nearly forced for the uniform-vs-selective
comparison; but a real ~16-40x Z-dephasing contrast coexists with the T1
sink, so the noise-contrast reading (B) is weakened on THIS chain, not
inapplicable. The
engineered-sink A-vs-B test used a genuine Z-rotation sink for exactly this
reason (see [the concentrator A-vs-B mechanism test](CONCENTRATOR_AB_MECHANISM_TEST.md)).

| Qubit | Role | T1 (us) | T2echo (us) | T2* est (us) |
|-------|------|---------|-------------|-------------|
| Q85 | **Concentrator** | 2.8 | **5.2** | ~4 |
| Q86 | Protected | 211.2 | 122.7 | ~61 |
| Q87 | Protected | 272.9 | 243.9 | ~98 |
| Q88 | Protected | 155.5 | 170.0 | ~68 |
| Q94 | Protected | 295.9 | 237.6 | ~95 |

### Configurations

Three configurations on the same chain, same Trotterized Heisenberg circuit:

1. **No DD** - All qubits at natural T2* (no refocusing)
2. **Uniform DD** - DD (X-X echo) on all 5 qubits (standard practice)
3. **Selective DD** - DD on Q86,87,88,94 only. No DD on Q85 (concentrator)

### Circuit

- Initial state: |+>^5 (Hadamard on each qubit)
- Heisenberg coupling: XX + YY + ZZ between neighbors (Trotterized)
- Trotter steps: 2, 4, 6, 8, 10 (times t = 1.0 to 5.0)
- Measurement: 9-circuit tomography for all 4 neighbor pairs
- Shots: 4000 per circuit
- Total: 135 circuits (5 time points x 9 tomo x 3 configs)

### DD Verification

Gate count (from the identical-circuit rebuild, provenance below) is
consistent with selective DD acting as intended:
- Selective DD: 2700 X-gates across 45 circuits (DD on 4 qubits)
- Uniform DD: 3240 X-gates across 45 circuits (DD on 5 qubits)
- No DD: 0 X-gates (only Heisenberg + tomography gates)
- Difference: 540 X-gates = DD on Q85 in uniform but not selective

(The counts are not in this run's frozen JSONs; the A-vs-B campaign's
identical-circuit build reproduces them, `added_x_selective = 2700` and
`added_x_uniform = 3240` in [data/ibm_ab_test_july2026/](../data/ibm_ab_test_july2026/).)


## Results

### Sum-MI Across Configurations

| t | No DD | Uniform DD | Selective DD | Sel/Uni | Sel/NoDD |
|---|-------|-----------|-------------|---------|---------|
| 1.0 | 0.0743 | 0.0475 | **0.0759** | **1.60x** | 1.02x |
| 2.0 | 0.0431 | 0.0352 | **0.0504** | **1.43x** | 1.17x |
| 3.0 | 0.0366 | 0.0212 | **0.0542** | **2.56x** | 1.48x |
| 4.0 | 0.0352 | 0.0157 | **0.0501** | **3.19x** | 1.42x |
| 5.0 | 0.0373 | 0.0131 | **0.0377** | **2.87x** | 1.01x |
| **avg** | **0.0453** | **0.0265** | **0.0537** | **2.02x** | 1.18x |

**The selective > uniform Sum-MI ordering holds at ALL 5 time points.** Average: 2.0x. Peak: 3.2x at t=4.0.

**Retro-insight (2026-07-05, from the A-vs-B campaign's review; the note
promised for this table):** the MI estimator carries a shot-noise bias floor of
+0.014-0.028 at 4000 shots (measured through this experiment's own
reconstruction pipeline), and the uniform-DD row above (0.013-0.048) is
comparable to that floor. The ratios' SIZE is therefore floor-contaminated
(at t=5 the denominator 0.0131 is at the floor; superseded 2026-07-11: the
re-analysis below, point 2, puts the t=5 floor at 0.0067, so 0.0131 is ~2×
the floor, not at it). The 5-of-5 ORDERING is more
stable than the ratio size, but not floor-independent either: the bias floor is
state-dependent and does not cancel between the selective and uniform configs
(different states), so with no error bars this is a suggestive direction, not a
significance claim. Also recorded there: |+⟩^5 is an exact fixed point of the
Heisenberg gates (exact also for the flown rxx·ryy·rzz per-bond decomposition,
since XX, YY, ZZ commute on a bond), so Sum-MI in this design is a noise-seeded observable (only
NON-uniform noise creates it), and created MI is a detectability proxy, not
the concentrator's protection metric. And the ordering's SIGN is expected under
EITHER reading, at least in the mild-contrast regime here: since Sum-MI is
noise-seeded and selective DD is by construction the MORE non-uniform
configuration, "selective creates more than uniform" is largely a design
artifact of which config is more non-uniform (the "more non-uniform → more
created MI" step is monotone only locally; it turns over at strong dose, as
Run 1's ceiling shows). So the ordering does not discriminate the concentrator
mechanism (B) from gate-cost (A) even in principle. See
[Concentrator A-vs-B Mechanism Test](CONCENTRATOR_AB_MECHANISM_TEST.md), the
structural-fact section and the reckoning.

### Re-analysis with error bars (2026-07-11, zero QPU)

The retro-note above left three quantitative questions open: the floor's size
per time point, the ordering's robustness, and what the ratio reads after
floor correction. All three are now answered from the frozen raw counts alone
([concentrator re-analysis script](../simulations/concentrator_reanalysis.py)
on [data/ibm_sacrifice_zone_march2026/](../data/ibm_sacrifice_zone_march2026/)).
Gate first: the script reproduces every published sum_mi and pair_mi cell from
raw_counts with deviation 0.0 before computing any statistics. Two seeds,
1000 bootstrap replicates and 1000 null draws each; all verdicts seed-stable.

**1. The ordering now carries error bars.** Bootstrapping the full 5-qubit
counts (multinomial per tomography setting) gives the ordering confidence:
the fraction of shot-noise resamples preserving selective > uniform. (This
is a bootstrap stability probability, not a p-value against an equal-MI
null; no such null is constructible from two distinct-state configs.) Raw:
0.99 / 0.96 / 1.00 / 1.00 / 1.00 at t = 1..5, jointly ≈ 0.95. One
refinement matters: the bias floor (point 2) is config-DEPENDENT and larger
for selective at every t (selective preserves more single-qubit coherence,
hence more linear-inversion bias), so the raw comparison is slightly
anti-conservative. Shifting by the floor-mean difference gives the honest
headline: jointly ≈ 0.93 (0.938 / 0.932 across the two seeds; the
adjustment bites only at t = 1, 2). The joint figure treats the five points
as independent, which they are under shot-noise resampling; they share one
calibration and one batch, so drift-like systematics sit outside it (the
closing paragraph). The ordering stands either way. What Remains item 1 is
closed.

**2. The uniform row is NOT floor-dominated; the 07-05 floor range was too
coarse.** The zero-true-MI null built from THIS run's measured per-setting
marginals (product of marginals, same shot numbers, same estimator) puts the
Sum-MI floor at 0.013-0.016 at t = 1 (config-dependent), falling to
0.0067-0.0082 at t = 5. The earlier +0.014-0.028 range matches only the purest
case (t = 1) and overestimates the later floors by about 2×. (The same
upward bias shows directly in the bootstrap: the resample means sit
≈ 0.007-0.011 above the published point estimates; the published values are
the point estimates, the bootstrap spread the error bar.) The measured
uniform-DD row sits ABOVE its own floor at all five points
(P(floor ≥ measured) ≤ 0.003 everywhere; at t = 5: 0.0131 vs floor
0.0067 ± 0.0016, about 4σ in null-spread units, ≈ 1.7σ once the measured
value's own shot noise is folded in; the weight sits on the early points,
t = 1, 2 at ≈ 3.8σ and 4.4σ combined). These σ are measured against the
marginals-conditioned floor, the heuristic baseline that point 3 qualifies. A from-below cross-check: the floor falls only
gently (factor ~2, essentially all of it between t = 1 and t = 2, flat
afterwards) and readout crosstalk would be flat in t, but the uniform row
falls steadily by factor 3.6; the steep t-dependent part is dynamically
created correlation, not estimator artifact. This does not contradict the fixed-point
argument: uniform DD is not uniform NOISE (Q85 stays at T2echo = 5.2 µs under
DD while the interior reaches 123-244 µs), so |+⟩^5 staying an exact product
is not forced in any of the three configs.

**3. Floor-correcting RAISES the ratios.** The floor subtraction is a
heuristic bias correction (the floor conditions on the measured marginals
only, while the true estimator bias also depends on the state's correlation
structure), so corrected values are approximate excesses, not unbiased MI
estimates. With that stated: subtracting each config's own floor mean,
Sel/Uni goes 1.60 / 1.43 / 2.56 / 3.19 / 2.87 (raw) to
1.8 / 1.5 / 3.2 / 4.8 / 4.6 (corrected point estimates, default seed). The
07-05 worry ("the denominator is at the floor") had the direction
backwards: the floor eats proportionally more of the smaller denominator.
Honest limit: from t = 4 the corrected denominator is small (uniform excess
+0.0087 and +0.0064 with shot-noise SE 0.0037 and 0.0035, i.e. 43% and 55%
relative error), so the late corrected values are factor-2 readings, not
two-digit numbers.

**4. The quiet fourth finding: selective vs no-DD is NOT resolved at the
ends** (ordering confidence 0.57 at t = 1, 0.54 at t = 5; only mid-window
t = 3, 4 reach ≈ 0.94-0.96; default seed, seed-stable). The real outlier of
the three-config table is the uniform row being LOW: uniform DD suppresses
the noise-seeded signal for everyone.

**What this changes and what it does not:** the record sharpens, the
attribution does not move. Sum-MI remains a noise-seeded transport proxy,
the ordering's sign remains expected under either reading (the retro-note
above), and Q85 remains a mostly-T1 sink outside the Z-dephasing rule's
scope. The March data carry more real correlation signal than the 07-05
caveat granted; what that signal measures (A vs B) is still open. One
systematic the bootstrap cannot see, and for a single-run ordering claim the
dominant one: the three configs ran as three separate jobs inside one batch,
submitted and collected in the order selective → uniform → no-DD (file stamps
50-60 s apart), so the config with the highest Sum-MI also ran first, and
slow calibration or thermal drift across those minutes is neither in the
error bars nor bounded by anything in this run. Only a repeat (What Remains
items 5, 6) bounds it.

### Surprise: Uniform DD Comes Last - Two Interpretations

Uniform DD performs worst of all three configurations. This is
counterintuitive (more protection should help). Two interpretations:

**Interpretation A (conservative): DD on bad qubits is harmful.**
DD on Q85 (T2=5.2 us) adds X-gates to a qubit that cannot be saved.
These gates introduce errors (gate infidelity, crosstalk to neighbors)
without extending Q85's coherence meaningfully. Any experienced
experimentalist would skip DD on Q85. The selective DD advantage
may simply reflect "don't waste gates on lost causes."

**Interpretation B (concentrator): Spatial noise contrast helps.**
The concentrator formula predicts that concentrating noise on one
edge maximizes information transfer through the chain. Removing DD
from Q85 increases the noise contrast between protected interior and
noisy edge. This contrast, not just the absence of wasted gates, is
what improves Sum-MI. The per-pair table below shows even the concentrator
pair (0,1) at 2.2x improvement under selective DD; note honestly that the
same table's own explanation for that number (Q86 freed of crosstalk from
Q85's DD gates) is an Interpretation-A mechanism, so the 2.2x does not
discriminate between the readings.

**What would distinguish A from B:** Test selective DD on a chain where
ALL qubits have good T2 (>150 us). If selective DD still wins by
removing DD from one good edge qubit, it's the contrast (B). If
selective DD only wins when the concentrator qubit is naturally bad, it's
the gate-error effect (A). Originally planned for April 9; run 2026-07-05
with an engineered sink and NOT settled there (see What Remains item 2).

### Per-Pair MI at t=3.0

| Pair | No DD | Uniform DD | Selective DD | Sel/Uni |
|------|-------|-----------|-------------|---------|
| (0,1) Q85-Q86 | 0.0062 | 0.0054 | **0.0117** | 2.2x |
| (1,2) Q86-Q87 | 0.0097 | 0.0042 | **0.0171** | 4.1x |
| (2,3) Q87-Q88 | 0.0149 | 0.0079 | **0.0149** | 1.9x |
| (3,4) Q88-Q94 | 0.0058 | 0.0037 | **0.0106** | 2.9x |

Protected pair (1,2) shows 4.1x improvement under selective DD.
The concentrator pair (0,1) still shows 2.2x improvement because the
neighboring Q86 benefits from not having crosstalk from Q85's DD gates.


## Analysis

### Why 2-3x Instead of 360x?

The simulation (N=5, formula with epsilon->0) predicted 360x improvement
over V-shape. Hardware delivers 2-3.2x over uniform DD. The gap comes from:

1. **DD is not perfect protection.** Simulation assumes epsilon->0 (zero
   noise on protected qubits). Real DD reduces noise but doesn't eliminate it.
2. **Gate errors dominate.** Each DD gate (X-gate) has ~0.1% error; at 60
   DD gates per circuit that is a ~6% chance of a DD-gate error per executed
   circuit (the earlier "2700 gates ≈ 2.7 errors" summed across 45
   independent circuits, an aggregation that is not commensurable with a
   per-circuit Sum-MI; corrected 2026-07-05).
3. **Crosstalk.** Gates on one qubit affect neighbors. DD on Q86 affects Q85.
4. **Trotter error.** The Trotterized Heisenberg circuit is an approximation.
5. **Measurement error.** Readout fidelity ~98-99% adds noise to all configs.

The simulation tests the *principle* (spatial dephasing profiles help).
The hardware tests the *implementation* (can DD approximate it?).
Both show the same direction on the same observable class (created Sum-MI);
whether the hardware advantage is the concentrator mechanism or gate-cost
avoidance is the open A-vs-B attribution (see "Does not yet show" and What
Remains item 2).

### The T2/T2* Effect

The hardware effect DD exploits: DD refocuses slow noise, pushing T2* toward
T2(echo). Without DD, qubits remain at their natural T2*. The contrast
between DD-protected and unprotected qubits comes from this gap:

- Protected (with DD): T2_eff ~ T2(echo) ~ 123-244 us
- Concentrator (no DD): T2_eff ~ T2* ≈ 3.7 us (Q85 is naturally terrible)

This ~33-65x contrast in effective coherence times is the setup the concentrator
profile prescribes (a noise contrast between protected interior and unprotected
edge). Two caveats stack on it, both above: (1) Q85's side of this contrast is mostly
amplitude damping (T1; ~66% at the un-echoed T2* level, ~92% echoed), NOT the
Z-dephasing the rule is proven for (Chain Selection channel caveat); the
pure-Z-dephasing part of the same contrast is ~16-40x (echo-to-echo it would
be only ~2.6-6.7x). (2)
Whether the contrast, rather than the gate-cost of decoupling a bad qubit,
carries the measured advantage is the open A-vs-B question (see "Does not yet
show" #1). (The T2* values here are estimated, not device Ramsey measurements on this
chain: T2echo divided by per-qubit Ramsey ratios 1.4 / 2.0 / 2.5 / 2.5 / 2.5,
anchored to the March 18 Ramsey runs and stored in chain_info. The A-vs-B
doc's flat /2.5 disclosure describes its own Kingston estimates, not these.)

### Comparison with Literature

| Method | System | Improvement | What's optimized |
|--------|--------|------------|-----------------|
| ENAQT (Plenio & Huelga 2008) | N=3 theory | ~2x | scalar uniform gamma |
| IBM Bayesian PST (2025) | N=4 hardware | +8% | J-couplings |
| **This work (simulation)** | **N=5 theory** | **360x** | **spatial gamma profile** |
| **This work (hardware)** | **N=5 ibm_torino** | **2-3.2x** | **selective DD** |

The rows are different quantities against different baselines (selective-vs-
uniform DD here; optimal-γ-vs-no-γ in ENAQT; coupling tuning at IBM), so the
table places this work in the landscape rather than ranking it; no direct
"exceeds" comparison is claimed (wording corrected 2026-07-05; the ratio-size
caveat in the Results section applies here too).

### Connection to the Palindromic Formula

The discovery path was: palindromic eigenstructure -> SVD response matrix
-> numerical optimizer -> analytical formula -> hardware test. The formula
(gamma_edge = N*gamma_base - (N-1)*epsilon) was derived from the palindromic
sensitivity structure and predicts 360x improvement at N=5 in simulation.

The hardware experiment does NOT test the formula directly, and (2026-07-05
review) does not confirm its mechanism either. It sets up the qualitative
picture "concentrate noise on one edge, protect the rest", but the measured
2-3.2x selective-vs-uniform ordering does NOT discriminate the concentrator
mechanism (Reading B) from gate-cost avoidance (Reading A) even in principle
(see the Results retro-note), and the natural sink Q85 was dominated by
amplitude damping (the T1 channel is the larger one at both slices, Chain
Selection caveat), so it does not cleanly test the Z-dephasing rule. So the
run is consistent with the direction; it does not confirm it. The quantitative
prediction (360x) additionally assumes perfect noise control (epsilon->0),
which DD cannot achieve.

A direct test of the formula would require controlling individual qubit
dephasing rates, which IBM hardware does not currently support. The DD
approximation is the closest available implementation.


## What This Shows (and What It Doesn't)

**Shows:**
1. **Selective DD > Uniform DD on real hardware.** 5/5 time points, average
   2.0x; the ordering is robust to shot noise (2026-07-11 re-analysis: joint
   bootstrap ordering confidence ≈ 0.93 floor-bias-adjusted, ≈ 0.95 raw; the
   ratio SIZE stays floor-sensitive, and floor-correcting raises it).
2. **DD on bad qubits appears harmful here.** The likely reading is that gate
   errors on Q85 outweigh DD's benefit (itself a gate-cost, Interpretation-A
   mechanism); single run, no error bars, so read as directional not established.
3. **Natural T2 variation is exploitable.** Q85's T2=5us is a feature, not a bug.

**Does not yet show:**
1. **That the concentrator formula is the reason.** Could be gate-error
   avoidance (Interpretation A) rather than noise-contrast benefit (B).
2. **Significance beyond shot noise.** The 2026-07-11 bootstrap covers shot
   noise only; between-job drift and readout systematics are not in the
   error bars. Single run. One chain.
3. **Reproducibility.** Untested on other chains or other days.
4. **Scaling.** Only N=5. The formula predicts similar advantage at N=7, 9.

## What Remains

1. **Error bars:** Bootstrap or jackknife on Sum-MI from shot counts.
   **Done 2026-07-11** ([concentrator re-analysis
   script](../simulations/concentrator_reanalysis.py), Results re-analysis
   subsection): ordering robust (joint ordering confidence ≈ 0.93
   floor-bias-adjusted), per-time-point floors quantified, floor-corrected
   ratios read higher, not lower.
2. **A vs B test (April 9):** Selective DD on a UNIFORM-T2 chain. If it still
   wins, it's the contrast (concentrator). If not, it's gate-error avoidance.
   **RUN 2026-07-05, still open** (three pre-registered runs, engineered sink;
   the post-flight audit found the created-MI observable 56-96% construction
   artifact at the two partial doses (channel-dominated at the ceiling) and
   a transport quantity rather than a protection metric, so A vs B is not
   settled; "injected noise strictly removes signal" is refuted, nothing more):
   [Concentrator A-vs-B Mechanism Test](CONCENTRATOR_AB_MECHANISM_TEST.md),
   the reckoning section.
3. **N=7 on hardware:** Longer chain with 10:00 April budget
4. **Noise injection:** Intentional Z-rotations on concentrator qubit for more
   contrast. **Done 2026-07-05** in the A-vs-B campaign (virtual-RZ engineered
   sink); before re-proposing, read its instrument-level pitfall: a FROZEN
   phase-path construction at partial dose is mostly shared classical
   randomness, not a dephasing channel (the AB doc's reckoning, Downgrade 1).
5. **Multiple chains:** Same experiment on different chip regions
6. **Reproducibility:** Repeat on a different day (calibration fluctuates)

## QPU Budget

IBM Quantum offers 10 minutes of free QPU time per 28-day rolling window
on ibm_torino (133 qubits). Beyond that: $96/minute. This entire project
runs on the free tier. Previous experiments (CΨ crossing validation,
Ramsey T2*, shadow hunt) consumed the March allocation. This experiment
ran inside the remaining 3:30 of the current cycle (~163 s of it, ~47 s left).

| Run | QPU time | Budget | Remaining |
|-----|----------|--------|-----------|
| Selective DD | ~50s | 210s | 160s |
| Uniform DD | ~50s | 160s | 110s |
| No DD | ~50s | 110s | ~47s |
| **Total used** | **~163s** | | **~47s left** |

(Row times are rounded estimates; the 163 s total and 47 s remainder are the
billed figures, which the ~50 s rows do not exactly sum to.)

Next allocation: 10:00 on April 9, 2026. Planned: A/B test on uniform-T2
chain and N=7 experiment. (March plan as written; the A/B test ran
2026-07-05 with an engineered sink, see What Remains item 2; N=7 remains
unrun.)

---

## Formula Validation (March 28, 2026)

**Retro-note (2026-07-05, from the A-vs-B reckoning's review; read before the
table):** this section's rates put COHERENCE rates in Lindblad slots: the
script builds c = √γ·Z (coherence decays e^{−2γt}) and feeds it γ = 1/T2*,
double the physical rate; the repo's canonical convention is γ = 1/(2·T₂),
and the script's own comment concedes the simplification. It also subtracts
no T1 contribution (for the T1-limited Q85 the pure-dephasing rate would be
≈ 0.046, not 0.268 MHz) and contains no amplitude damping. Further, J was
selected by fitting to the same hardware ratios the table reports agreement
against (the fit objective is the t = 1 and t = 5 hardware ratios, so both
those rows are circular; t = 1 is the only row that both entered the fit and
agrees), and the "t (μs)" axis rests
on an unstated 1 J-time = 1 µs identification. So read the 6-13% figures as
the residual of a one-parameter fit of a doubled-rate dephasing-only model,
not as an independent match of the formula. This is the same rate-convention
twin the A-vs-B reckoning documents (its factor-2 finding), in an earlier
costume.

The concentrator formula (Lindblad simulation with per-qubit gamma)
was compared quantitatively against the hardware measurements.

Echo-aware gamma profiles used:
- Selective DD: Q85 at gamma = 1/T2* = 0.268 MHz, others at 1/T2 ~ 0.004-0.008 MHz
- Uniform DD: all at 1/T2 (Q85 at 0.192 MHz, others same as selective)

| t (μs) | Simulation ratio | Hardware ratio | Agreement |
|--------|-----------------|----------------|-----------|
| 1.0 | 1.39x | 1.60x | 13% |
| 2.0 | 1.35x | 1.43x | 6% |
| 3.0 | 1.31x | 2.56x | Diverges |
| 4.0 | 1.29x | 3.19x | Diverges |
| 5.0 | 1.33x | 2.87x | Diverges |

(The t = 4.0 row was missing from the March table; added 2026-07-11 from the
same script, so the table shows all five measured points: 2 of 5 within the
fit residual, 3 of 5 diverging.)

At t=1-2 μs the formula's ratio sits within 6-13% of hardware, but per the
retro-note above this is the residual of a circular one-parameter (J) fit to
these same ratios, not an independent match. At later times the
hardware advantage grows beyond the formula prediction because DD
gates on Q85 (T1 = 2.84 μs, gate time ~ 0.5 μs) introduce T1
losses that accumulate with each DD cycle. The formula models DD
as perfect gamma reduction; real DD on a T1-limited qubit adds
errors that worsen the uniform configuration.

The hardware imperfections amplify the measured advantage of selective DD:
they make the concentrator qubit WORSE under uniform DD. Note honestly that
this amplification channel (DD-gate T1 losses on Q85) is itself a gate-cost
(Interpretation-A) mechanism, so it is not evidence for the concentrator
reading (clarified 2026-07-05).

Script: [`simulations/ibm_formula_test.py`](../simulations/ibm_formula_test.py).

---

## References

- [Resonant Return (formula + optimization)](RESONANT_RETURN.md)
- [IBM Run 3 (1.9% CΨ validation)](IBM_RUN3_PALINDROME.md)
- [gamma as Signal](GAMMA_AS_SIGNAL.md): 15.5 bits channel capacity
- [gamma Control](GAMMA_CONTROL.md): V-shape 21.5x baseline
- [Inside and Outside the Sacrifice Zone](../docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md): this hardware run read from both sides of the γ-boundary (sacrifice vs concentration)
- Raw data: [data/ibm_sacrifice_zone_march2026/](../data/ibm_sacrifice_zone_march2026/)
