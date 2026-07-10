# IBM Hardware: Selective DD Beats Uniform DD at All Five Time Points (Concentrator Mechanism Open)

<!-- Keywords: IBM Torino hardware validation concentrator profile, selective dynamic
decoupling beats uniform DD, first spatial dephasing profile hardware test,
Q85 natural concentrator qubit T2=5us, Sum-MI selective vs uniform 2x-3.2x,
palindrome-derived noise engineering quantum hardware, R=CPsi2 IBM experiment -->

**Status:** Hardware, Tier 2: the selective > uniform Sum-MI ORDERING holds at
all 5 time points. Single run, one chain, no error bars. The ratio SIZE
(2-3.2×) is floor-caveated and Sum-MI here is a noise-seeded transport proxy
(2026-07-05 retro-note in Results); the A-vs-B attribution is open.
**Date:** March 24, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Hardware:** ibm_torino (Heron r2, 133 qubits)
**Chain:** Q85-Q86-Q87-Q88-Q94 (concentrator: Q85, T2echo=5 us)
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

First hardware test of spatially selective Dynamic Decoupling on a
quantum computer. On a 5-qubit Heisenberg chain (ibm_torino,
Q85-86-87-88-94), removing DD from one edge qubit (Q85, natural T2=5 us)
while protecting the other four produces 2.0x higher Sum-MI on average
and up to 3.2x at t=4.0, compared to standard uniform DD on all qubits.
Selective DD beats uniform DD at all 5 measured time points.

**Caveats:** Single run on one chain. No error bars yet (4000 shots per
circuit). The advantage may partly reflect that DD on Q85 (T2=5 us)
adds gate errors without benefit, rather than a pure concentrator
effect. Reproducibility across chains and days is untested. The connection
to the palindromic eigenstructure is indirect: the formula predicted the
optimal noise profile, but the hardware implementation approximates it
via DD rather than controlling dephasing rates directly. Added 2026-07-05
(retro-note in Results): the ratio SIZE is floor-contaminated, and Sum-MI on
|+⟩^5 is a noise-seeded transport proxy, not a protection metric.


## Setup

### Chain Selection

ibm_torino calibration (March 24, 2026) identified Q85 as an extreme
outlier: T2(echo) = 5.2 us while neighboring qubits have T2 = 123-244 us.
This provides a natural concentrator qubit with 51.9x effective contrast
(T2echo_protected / T2*_concentrator).

**Channel caveat (2026-07-05 review, important):** Q85's short T2echo is
~93% amplitude damping, not Z-dephasing. With T1 = 2.8 us, 1/(2·T1) = 0.179
us⁻¹ is 92% of 1/T2echo = 0.192 us⁻¹, so Q85's PURE-dephasing time is
Tφ ≈ 73 us, only ~2-6x shorter than the interior (Tφ ≈ 170-440 us), not 52x.
The concentrator selection rule and its Absorption Theorem are proven for
Z-dephasing only (amplitude damping is the separate F82/F84 object, and DD
cannot refocus T1). So this natural-sink run does not cleanly test the
Z-dephasing rule: Q85 is mostly a T1 sink, which makes the gate-cost reading
(Interpretation A, DD-on-a-T1-limited-qubit hurts) nearly forced and leaves
the noise-contrast reading (B) largely inapplicable to THIS chain. The
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

Gate count confirms selective DD works as intended:
- Selective DD: 2700 X-gates across 45 circuits (DD on 4 qubits)
- Uniform DD: 3240 X-gates across 45 circuits (DD on 5 qubits)
- No DD: 0 X-gates (only Heisenberg + tomography gates)
- Difference: 540 X-gates = DD on Q85 in uniform but not selective


## Results

### Sum-MI Across Configurations

| t | No DD | Uniform DD | Selective DD | Sel/Uni | Sel/NoDD |
|---|-------|-----------|-------------|---------|---------|
| 1.0 | 0.0743 | 0.0475 | **0.0759** | **1.60x** | 1.02x |
| 2.0 | 0.0431 | 0.0352 | **0.0504** | **1.43x** | 1.17x |
| 3.0 | 0.0366 | 0.0212 | **0.0542** | **2.56x** | 1.48x |
| 4.0 | 0.0352 | 0.0157 | **0.0501** | **3.19x** | 1.42x |
| 5.0 | 0.0373 | 0.0131 | **0.0377** | **2.87x** | 1.01x |
| **avg** | **0.0453** | **0.0265** | **0.0537** | **2.02x** | 1.19x |

**Selective DD beats Uniform DD at ALL 5 time points.** Average: 2.0x. Peak: 3.2x at t=4.0.

**Retro-insight (2026-07-05, from the A-vs-B campaign's review; the note this
table was promised):** the MI estimator carries a shot-noise bias floor of
+0.014-0.028 at 4000 shots (measured through this experiment's own
reconstruction pipeline), and the uniform-DD row above (0.013-0.048) is
comparable to that floor. The ratios' SIZE is therefore floor-contaminated
(at t=5 the denominator 0.0131 is at the floor). The 5-of-5 ORDERING is more
stable than the ratio size, but not floor-independent either: the bias floor is
state-dependent and does not cancel between the selective and uniform configs
(different states), so with no error bars this is a suggestive direction, not a
significance claim. Also recorded there: |+⟩^5 is an exact fixed point of the
Heisenberg gates, so Sum-MI in this design is a noise-seeded observable (only
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

### Surprise: No DD Also Beats Uniform DD - Two Interpretations

Uniform DD performs worst of all three configurations. This is
counterintuitive (more protection should help). Two interpretations:

**Interpretation A (conservative): DD on bad qubits is harmful.**
DD on Q85 (T2=5 us) adds X-gates to a qubit that cannot be saved.
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
over V-shape. Hardware delivers 2-3x over uniform DD. The gap comes from:

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

- Protected (with DD): T2_eff ~ T2(echo) ~ 170-244 us
- Concentrator (no DD): T2_eff ~ T2* ~ 4 us (Q85 is naturally terrible)

This 40-60x contrast in effective coherence times is the setup the concentrator
profile prescribes (a noise contrast between protected interior and unprotected
edge). Two caveats stack on it, both above: (1) this is a T2echo contrast, ~93%
of which is amplitude damping (T1), NOT the Z-dephasing the rule is proven for
(Chain Selection channel caveat); the pure-dephasing contrast is only ~2-6x. (2)
Whether the contrast, rather than the gate-cost of decoupling a bad qubit,
carries the measured advantage is the open A-vs-B question (see "Does not yet
show" #1). (The T2* values here are estimated, not device Ramsey measurements;
the A-vs-B doc discloses the T2echo/2.5 assumption they inherit.)

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
2-3x selective-vs-uniform ordering does NOT discriminate the concentrator
mechanism (Reading B) from gate-cost avoidance (Reading A) even in principle
(see the Results retro-note), and the natural sink Q85 was ~93% amplitude
damping, outside the Z-dephasing rule's scope (Chain Selection caveat). So the
run is consistent with the direction; it does not confirm it. The quantitative
prediction (360x) additionally assumes perfect noise control (epsilon->0),
which DD cannot achieve.

A direct test of the formula would require controlling individual qubit
dephasing rates, which IBM hardware does not currently support. The DD
approximation is the closest available implementation.


## What This Shows (and What It Doesn't)

**Shows:**
1. **Selective DD > Uniform DD on real hardware.** 5/5 time points, average
   2.0x (a point-estimate direction; the bias floor is state-dependent and
   does not cleanly cancel between configs, and there are no error bars, so no
   formal significance is claimed).
2. **DD on bad qubits appears harmful here.** The likely reading is that gate
   errors on Q85 outweigh DD's benefit (itself a gate-cost, Interpretation-A
   mechanism); single run, no error bars, so read as directional not established.
3. **Natural T2 variation is exploitable.** Q85's T2=5us is a feature, not a bug.

**Does not yet show:**
1. **That the concentrator formula is the reason.** Could be gate-error
   avoidance (Interpretation A) rather than noise-contrast benefit (B).
2. **Statistical significance.** No error bars. Single run. One chain.
3. **Reproducibility.** Untested on other chains or other days.
4. **Scaling.** Only N=5. The formula predicts similar advantage at N=7, 9.

## What Remains

1. **Error bars:** Bootstrap or jackknife on Sum-MI from shot counts
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
used the remaining 3:30 from the current cycle.

| Run | QPU time | Budget | Remaining |
|-----|----------|--------|-----------|
| Selective DD | ~50s | 210s | 160s |
| Uniform DD | ~50s | 160s | 110s |
| No DD | ~50s | 110s | ~47s |
| **Total used** | **~163s** | | **~47s left** |

(Row times are rounded estimates; the 163 s total and 47 s remainder are the
billed figures, which the ~50 s rows do not exactly sum to.)

Next allocation: 10:00 on April 9, 2026. Planned: A/B test on uniform-T2
chain and N=7 experiment.

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
against (circular for at least the t = 1 row), and the "t (μs)" axis rests
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
| 5.0 | 1.33x | 2.87x | Diverges |

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
