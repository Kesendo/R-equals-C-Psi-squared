# IBM F129: The Standing Fringe

**Status:** FLOWN AND CONFIRMED (2026-07-15, ibm_kingston, job
d9br4vmg26ic73dgbgk0; see RECORD). Pre-registered at commit 5bf3b79 (v2.5, 7a +
7b recorded) BEFORE the shot; the flight followed the pinned order: fresh
calibration + §8 hard aborts → Tom's explicit go → ONE job → the RECORD below.
The runner: `run_f129_ramsey_fringe.py` in the external IBM pipeline
(AIEvolution.UI/experiments/ibm_quantum_tomography), importing this repo's
committed gates as the single source of circuits, estimator and verdict.
**Date:** 2026-07-15
**Authors:** Thomas Wicht & Claude
**Design gate:** [`simulations/f129_ramsey_fringe_design.py`](../simulations/f129_ramsey_fringe_design.py)
(every verdict-determining number in §5 is printed by it; the few doc-side
arithmetic steps on gate primitives are marked as such)
**Hardware target:** one Heron r2 device (ibm_kingston class), one job, ~6.6-7.7 QPU minutes
**Relation:** the law on display is [F129](../docs/ANALYTICAL_FORMULAS.md) (proof:
[PROOF_F129_LEVEL_COLLISION_LAW.md](../docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md));
the protocol this replaces is recorded dead in
[F130_HW_INFEASIBILITY.md](F130_HW_INFEASIBILITY.md); genre parent:
[IBM_CONCENTRATOR_RELOADED.md](IBM_CONCENTRATOR_RELOADED.md).

## What this document is about

F129 says the comb's level map is injective away from 3|n and 10|n: two distinct clean
mode triples can share an energy level only there, and at n = 9 the collision
(1,5,7) ~ (2,4,8) is exact. On hardware that collision is audible as a fringe that
stands still. Prepare the superposition of the two triples as 3-magnon states, let the
chain step forward M times, and read the relative phase: for the colliding pair the
phase slope is pinned to a small computed third-order drift, while every detuned pair
winds at first order. The concert hall's law, played on a real device: one note that
does not beat. Scoping, said plainly: one flight at one n tests the EXISTENCE half of
F129 at its smallest firing n (an exact clean-clean coincidence standing against the
winding census); the injectivity half stays with the proof.

## 1. The question

Does the physical device show the F129 collision as a standing Ramsey fringe: the
collision arm's phase slope consistent with the computed third-order Trotter drift
(and nothing else), while the verdict control winds at its predicted first-order
rate, at ≥ 5σ separation?

## 2. What the theory predicts (from below), and the regime choice

The comb at n = the single-particle spectrum of the open XX chain with N = n−1 sites:
modes k = 1..N at energies cos(kπ/n). A triple's level S(τ) = Σ cos(kᵢπ/n). On the
stepped device the energies deform into the Floquet phases of the Trotter step
(odd/even Givens layering); the deformation is smooth and small (worst mode overlap
with the sine modes 0.9820 at N = 8, θ = 0.5), and the collision survives it to third
order: the primary pair's drift has fitted θ-exponent 3.0054 over θ ∈ [0.02, 0.5]
(design gate scan), with coefficient 1/6 in the θ → 0 limit rising to the effective
0.170 at the flown θ = 0.5.

**The lab, n = 9 (N = 8):** the collision (1,5,7) ~ (2,4,8). This is the F129 law's
clean exemplar (both triples clean, 3|9, the smallest firing n), it is disjoint (six
distinct modes, so the two-branch superposition is a 6-qubit cat in the mode register),
and it is the MIRROR pair (σ = n−τ), which buys the flight its best systematic. The
active symmetry is the CHIRAL (sublattice) map behind the mirror: the Floquet modes
obey V_{n−k} = C·conj(V_k) with C = diag((−1)^j), so the two branches have identical
Z-diagonal observables PER SITE and for EVERY pair (i, j), any range, not just
nearest-neighbour bonds (design gate certificate: machine-zero at every θ; the
all-range statement covers real-device ZZ crosstalk to spectator qubits). The
consequence is stronger than a sum rule: the collision arm's first-order drift
vanishes for the always-on ZZ coupling AND for single-qubit Z-detuning, even when
ζ_ij and δ_i are SITE-DEPENDENT (each pair and site cancels term by term), at any
θ. A degeneracy has one further first-order channel the diagonal argument does not
touch: coherent BRANCH MIXING, ⟨σ|z_i z_j|τ⟩ or ⟨σ|n_i|τ⟩ Rabi-coupling the two
exactly degenerate branches. It vanishes identically by Slater-Condon: the branches
differ in three occupied orbitals while the operators are at most 2-body (design
gate certificate: ≤ 1.8·10⁻¹⁵ over all i, j; disjointness is load-bearing here, not
just cat bookkeeping). The dangerous coherent per-step systematics of a real chain
are retired at FIRST order by the pair choice itself. The price the mirror charges
(a 7a find): the chiral map negates the Floquet spectrum, so the SECOND-order ZZ
shifts of the two branches are OPPOSITE and the fringe sees twice one of them (the
symmetry fixes the factor 2, not the sign); the 7a statevector scan pinned the law
under the flown estimator, bias = 0.00257·(ζ/3.8 kHz)²·(τ_step/1.2 µs)² rad/step
(1.2 µs = the law's normalization point). The bias direction is a property of the
same-sign always-on regime: over 40 all-positive ζ draws the gate finds the bias
strictly positive, while mixed-sign ζ patterns push at most ~1.3·10⁻³ downward
at the 1.2 µs normalization point (gate-printed, 40 mixed-sign draws at
|ζ| ≤ 7.6 kHz), i.e. ~5·10⁻⁴ at the flown τ_step = 0.7 µs, which the two-sided
b_qs budget covers. Convention, pinned: ζ is the Ramsey-shift ZZ, H_ZZ = πζ·Z_iZ_j,
the convention of the price_pair conditional-Ramsey measurement that anchors
the 3.6-3.9 kHz transfer. At the flown τ_step = 0.7 µs the nominal bias is
9·10⁻⁴; it is budgeted in clause (a), not assumed away. The pair is resonant (S = 0); what the
fringe measures is the level EQUALITY, and resonance is simply this pair's value of it.

**Regime:** hop angle θ = 0.5 per Trotter step. The XX chain is free-fermionic, so the
matchgate step is EXACT at any θ (no many-body Trotter error; the one-particle Floquet
spectrum is the whole truth), and the larger angle nearly doubles every first-order
winding rate while the collision drift stays third-order: predicted +0.0213 rad/step,
itself a from-below θ³ signature the flight will measure. Pure hopping, NO engineered
dephasing: F129 is a spectral law of the closed chain, the fringe is a coherent phase
measurement, and the right regime is the longest coherence. Shots 16384 per circuit.

The n = 12 equal-level pair (1,2,10) ~ (3,5,6) is NOT part of this pre-registration
(see §9): (1,2,10) is not clean (2+10 = n), the pair tests F130's equal-level
hypothesis rather than F129's law, and the design gate's appendix shows its nearest
impostor (dS = 0.018) is inseparable at this budget.

## 3. Instruments: three arms, one network, one job

**State preparation (the enabler):** the two branch states are prepared as exact
eigenstates of the stepped evolution itself. The Givens network is compiled from the
eigenvectors of the actual Trotter step (classical preprocessing; NOT the sine modes,
so there is no prep leakage by construction), 28 Givens rotations at N = 8. The seed
is a 6-qubit GHZ-type superposition of the two occupation bitstrings (τ-modes vs
σ-modes) in the mode register, walked through the network into the site basis.

**Evolution:** M Trotter steps of the pure hopping chain, M ∈ {0, 1, 2, ..., 8}.

**Readout:** the inverse network, un-GHZ, then the seed qubit measured in X and Y
(two quadratures per (arm, M) point), readout-corrected by the in-job calibration
PUBs. The fringe phase is exactly M·ΔΦ(pair), ΔΦ the difference of the two 3-mode
Floquet phase sums.

**Arms, all in one interleaved job (round-robin, the strongest arm must not run
first; concentrator rule carried in verbatim):**

| Arm | Pair | dS | Role |
|-----|------|-----|------|
| A0 | (1,5,7) ~ (2,4,8) | 0 (exact) | the collision |
| A1 | (1,5,7) ~ (2,5,6) | 0.092 | secondary dial (near detuning, opposite winding sign) |
| A2 | (1,5,7) ~ (1,6,7) | 0.326 | VERDICT control |

3 arms × 9 M-points × 2 quadratures = 54 science circuits, plus 2 readout-calibration
PUBs and 8 in-situ T1/T2* PUBs interleaved (concentrator pattern), ≈ 64 PUBs total.
(A second near-control at dS = 0.092, (3,4,7), was carried through v2 and dropped in
v2.2 for the billing cushion; it duplicated A1's detuning magnitude and decided
nothing.)

**Instruments (v2.5 status):** the Givens network compiler EXISTS and is certified
machine-zero ([`f129_givens_compiler.py`](../simulations/f129_givens_compiler.py):
per-arm column permutations put each arm's differing modes contiguous, so the cat
costs 5/3/1 CX for A0/A1/A2; branch signs tracked from the column-permutation
parity; circuit-vs-direct Slater < 6·10⁻¹⁵, full-chain fringe = M·ΔΦ exact); the
seed + estimator run end-to-end in the 7a gate
([`f129_ramsey_7a_gate.py`](../simulations/f129_ramsey_7a_gate.py)); and the
HARDWARE runner exists (`run_f129_ramsey_fringe.py`, external pipeline), with
modes --calibrate / --certify / --simulate (the 7b gate) / --hardware (gated
behind --yes) / --analyze, counts persisted BEFORE any reduction (the
concentrator lesson), and the repo gates imported as the single source of
circuits, estimator and verdict (no copied physics). Per-arm CX counts:
112 + 2(cat−1) + 14M, so the shared A0-formula V-model over-counts the smaller
control cats, conservative as pinned.

## 4. Dose

There is no engineered dephasing dose in this experiment. The knobs are pinned:
θ = 0.5 per step, M-grid {0, 1, ..., 8}, 16384 shots per circuit, one Heron r2
device, all arms transpiled from one skeleton. Timing (v2.4): steps run at gate
speed with NO idle padding, τ_step ≈ 0.7 µs wall at Heron-era durations (150 ns
2q, 50 ns 1q; the 1.2 µs padded step of v1-v2.3 was concentrator legacy, and
padding only feeds the ζ² systematic and the dephasing); durations reconciled at
7b against the actual backend.

## 5. Statistics (pinned before the runner exists; verdict-determining numbers gate-printed)

**Estimator:** per (arm, M) the readout-corrected fringe phase φ = atan2(⟨Y⟩, ⟨X⟩)
and amplitude V_meas = √(⟨X⟩² + ⟨Y⟩²). Per arm, the phase is unwrapped over the
M-grid (worst predicted increment 0.32 rad/step × ΔM = 1 < π/2, ample margin) and
fitted by weighted least squares, φ(M) = slope·M + b. **Weights are computed from the
PREDICTED visibility model, not from V_meas** (removes weight-estimation bias);
V_meas is used only by the floor rule below. The primary numbers are the fitted
slopes; σ_slope is the weighted-LSQ error, and every two-arm comparison uses the
combined error √2·σ_slope.

**Error model (from below, re-frozen once at 7a):** σ_φ = 1/(V·√shots) per point;
V(M) = 0.75 · (1−p₂)^CX(M) · e^(−3t/T1) · e^(−m_deph·t/T2*), with m_deph the exactly
computed FIRST MOMENT of the branch-pair Hamming distance (the two Slater branches'
site-configuration distributions; 3.68 for A0, each control arm uses its own pair's
value, ~3.8) and 0.75 the one-qubit-gate attenuation factor the 7a gate re-froze
(transpiled sx count 336 + 42·M at 50 ns each; sim/old-model ratio 0.79-0.82,
pinned conservative, and the 7a mixture check verifies V_sim ≥ model at every M).
The 6-qubit branch cat dephases ~3.7× faster than one qubit; modeling it as one
collective DOF would understate σ, and the single-exponential model built on the
first moment is conservative (Jensen: the true multi-exponential envelope lies
above it). Two-arm comparisons combine √(σ_A0² + σ_arm²). Projected at p₂ = 0.3 %:
V(M=8) = 0.215, σ_slope(A0) = 0.0033 rad/step. Sensitivity, pinned: a uniform
V-model miss scales σ_slope linearly (V/2 → σ ×2); it erodes the margins below,
not the verdicts, which survive a full 2× miss.

**ZZ and detuning systematic policy (one treatment, used everywhere):** the
collision arm carries ZERO budget for ZZ and single-qubit detuning, by the chiral
all-range/per-site cancellation of §2 plus the Slater-Condon branch-mixing zero
(exact for site-dependent ζ_ij, δ_i; both design-gate certificates at machine
precision). Control arms carry exact first-order Wick values
(one-body density matrices of the 3-magnon branches), floored at 5× their magnitude
as the device-transfer safety factor (ζ measured 3.6-3.9 kHz on ibm_marrakesh; 5×
covers ζ up to ~20 kHz, far beyond a flyable chain). The one per-step coherent
systematic the symmetry does NOT retire on A0 is hop-angle miscalibration:
d(drift)/dθ = 0.1294, so a 1 % angle error moves A0 by 0.20σ_slope (doc-side arithmetic on gate primitives); it is closed by
the θ̂ rule in clause (a).

**Predictions (θ = 0.5, from the design gate):**

| Arm | Floquet slope (rad/step) | ZZ (Wick, ζ = 4 kHz, τ = 0.7 µs) | projected separation from A0 |
|-----|--------------------------|-----------------------------------|------------------------------|
| A0 | +0.0213 (= 0.170·θ³) | 0 at first order (chiral); ζ² budgeted | |
| A1 | −0.0879 | −0.0024 | 20.7σ (secondary) |
| A2 | +0.3204 | +0.0010 | 62.4σ (verdict; margin 12.5× vs 5σ) |

Separations are the true slope gaps |slope_arm − slope_A0| (sign-safe; A1 winds
opposite to the A0 drift) in combined error (per-arm m_deph), floored by 5× the
Wick ZZ.

**Clause-(a) systematic budgets (pinned, gate-printed by
[`f129_ramsey_7a_gate.py`](../simulations/f129_ramsey_7a_gate.py)):**
b_zz2 = +0.0043 one-sided (the ζ² anti-protection at 2× the transferred ζ, i.e.
ζ ≤ 7.6 kHz, at τ_step = 0.7 µs, with a 1.22 estimator-pattern uplift; exact
statevector law 0.00257·(ζ/3.8 kHz)²·(τ/1.2 µs)², gate mode zz2_scan; one-sided
because the always-on regime is same-sign, checked over 40 positive draws), and
b_qs = ±0.008 two-sided (quasi-static site-disorder at the T2*-consistent scale
σ_δ ≈ 3 kHz pulls the A0 slope through prep-basis mismatch; 7a measured
−0.0075 ± 0.0037 over ten disorder draws; on one device it is a single frozen draw
of unknown sign, hence two-sided; it also covers the ≤ 10⁻³ mixed-sign ζ
excursions and the in-sim decoherence pull ~−0.002).

**Verdict rules (v2; the verdict is ONE conjunction, the secondaries carry no verdict
weight, so there is no multiple-comparison inflation):**

- **(a) The standing fringe:** A0's fitted slope lies within
  center + [−(3σ_a + b_qs), +(3σ_a + b_zz2 + b_qs)] (the budgets above; the ζ²
  term pushes only upward). The center slope_A0(θ̂) is a pinned FUNCTION, not a
  constant: in the linearized form emitted by the design gate,
  center = 0.0213 + 0.1294·(θ̂ − 0.5) rad/step (= +0.0213 at θ̂ = 0.5;
  linearization error ≤ 4.6·10⁻⁵ over the admissible band). The θ̂ rule
  (pinned, deterministic, closes the hop-angle exposure):
  θ̂ = 0.5 + (slope_A2_measured − 0.3204)/0.6164; A2 is 4.8× more θ-sensitive than
  A0, so it reads the effective angle out of the same job, and clause (c) bounds
  |θ̂ − 0.5| ≤ 0.024 whenever it holds. σ_a is the FROZEN design projection
  0.0033 rad/step inflated by the propagated θ̂ variance,
  σ_a = σ_slope·√(1 + (0.1294·σ_θ̂/σ_slope)²) = 1.022·σ_slope (gate-printed; A0 and
  A2 shots are independent, so there is no estimator correlation, only this
  variance term); re-frozen once at 7a (done, see 7a RECORDED) and finally at 7b
  per the null-band rule below, never after data. Discrimination power: the
  nearest impostor by this clause's own metric (the minimum of
  |slope_c − slope_A0|/σ over all 54 distinct triples, design gate) is (3,4,7) at
  dS = 0.092, predicted slope +0.1092, landing 26.5σ from the window center
  (margin 5.3× vs 5σ) and ~19σ beyond even the budget-widened window edge, so
  clause (a) alone rejects every "accidental near-degeneracy" alternative F129
  permits.
- **(b) The winding control:** A2's slope differs from A0's by ≥ 5σ combined
  (realized post-floor errors, per-arm m_deph), with the predicted sign.
  Projected 62.4σ, margin 12.5×.
- **(c) Validity:** A2's slope agrees with its NOMINAL θ = 0.5 prediction
  (+0.3204 + ZZ) within ±(3σ_slope + 5×|ZZ_Wick|), realized errors. Nominal, not
  θ̂-evaluated: comparing A2 to a prediction at the angle inferred from A2 itself
  would be vacuous; against the nominal angle, (c) is exactly what bounds
  |θ̂ − 0.5| ≤ 0.024. If (c) fails, the fringe clock is broken (a > 2.4 % angle
  miss or worse): the run is VOID (instrument), never a physics verdict.
- **STANDING FRINGE CONFIRMED** iff (a) ∧ (b) ∧ (c).
- **VIOLATED (collision winding)** iff A0's slope is outside
  center + [−(5σ_a + b_qs), +(5σ_a + b_zz2 + b_qs)] while (c) holds. The gap
  between (a) and VIOLATED is a deliberate dead band: it reads inconclusive.
- Anything else: **instrument failure or inconclusive**, never a physics verdict.
- **Secondary dial A1:** prediction pinned above, result reported with the same
  estimator; it sharpens the story (the clock resolves dS = 0.092, winding the
  other way) but decides nothing.

**Point-floor rule (pinned):** a science point with V_meas < 0.08 is dropped
(predicted minimum 0.215, margin ~2.7×). An arm is valid iff ≥ 6 of its 9 points
survive, including ≥ 2 with M ≤ 2 and ≥ 2 with M ≥ 6. A0 or A2 invalid → run VOID;
A1 invalid → recorded, no verdict impact. Clause (a)'s σ_a stays the frozen
projection regardless of dropped points (the drop inflates only the realized errors
of (b)/(c), conservatively); a drop pattern that passes the quorum inflates the
realized σ_slope by at most ×1.54 (gate-printed worst case over all
quorum-passing survivor sets), which the (b) margin 12.5× absorbs.

**False-positive note:** CONFIRMED is a conjunction whose clause (a) window center
sits 26.5σ from the nearest alternative (≈ 19σ beyond even the budget-widened
edge); the projected false-CONFIRM probability under the design's own noise model
is ≈ 0. The design's exposure is false-negative (power), which the margins above
cover, including at the p₂ stress point 0.5 % and under a uniform 2× V-model miss
(margins halve, clauses survive).

**Null band and bootstrap (DISCHARGED at 7a/7b, see §6):** the 7a gate produced
(i) the frozen null band for clause (a) from a zero-drift synthetic arm through
this exact estimator, and (ii) hierarchical bootstrap SEs; the re-freeze trigger
(bootstrap SE > 1.3× analytic) did not fire (ratios 0.90/0.63/0.50), so the
analytic bands stand, confirmed unchanged at the 7b gate.

## 6. The two gates (RECORDED below, before any shot)

- **7a (from-below sim):** full-circuit Aer simulation of the flown construction
  under the anchored noise model (p₂, T1/T2*, readout, ZZ = 4 kHz), through the §5
  estimator. Required: the CONFIRMED verdict fires in sim with ≥ 3× power margin on
  clauses (a) and (b); the null band and bootstrap of §5 recorded. Named 7a checks
  beyond the verdict: branch-dependent decoherence (the two Slater branches see
  different site weights), quasi-static vs Markovian dephasing (T2* model), the
  visibility mixture model (m_deph is a first moment of a POSITIVE mixture; the
  physical fringe sums complex amplitudes, which 7a checks end-to-end), the
  branch-mixing zero end-to-end (a simulated site-dependent ZZ + detuning layer must
  leave the A0 fringe slope at its θ³ value), readout correction efficacy, and
  magnon-number leakage feeding spurious phase. RECORDED below.
- **7b (counts level):** the actual runner's circuits, the actual transpilation,
  counts through the actual estimator code; band reconciliation recorded. Hard
  pre-flight abort if 7b contradicts 7a. RECORDED below.

### 7a RECORDED (2026-07-15, gate [`f129_ramsey_7a_gate.py`](../simulations/f129_ramsey_7a_gate.py), seeds pinned in the gate)

Full Aer simulation of the flown construction (per-arm compiled Givens networks,
certified machine-zero noiselessly in [`f129_givens_compiler.py`](../simulations/f129_givens_compiler.py));
anchored noise model p₂ = 0.3 % + T1/T2* = 200/70 µs at Heron-era durations +
1 % readout (corrected in-estimator) + coherent site-dependent NN ZZ
ζᵢ = 3.8 kHz ± 30 % (pinned seed); transpiled CX count 234 at M = 8 = the budget
model exactly.

- **Verdict: CONFIRMED fires on every seed** (six independent noise seeds;
  dev_a ∈ [−0.41, +1.29] σ_a, sep_b ≈ 63σ). Power margins: clause (a) impostor
  5.2×, clause (b) 12.6×, both ≥ 3× as required. The verdict logic itself is
  exercised in-gate: synthetic A0 slopes at the window edges route to
  CONFIRMED / inconclusive (both sides) / VIOLATED (both sides) as pinned, and
  the budgets b_zz2, b_qs are gate-printed constants of the same run.
- **Null band (frozen):** H0 parametric bootstrap through this estimator,
  [−0.0088, +0.0089] rad/step around the θ̂ center (3σ quantiles, 400 draws);
  the analytic 3σ_a = 0.0102 encloses it (conservative).
- **Hierarchical bootstrap vs analytic σ:** ratios A0 0.90, A1 0.63, A2 0.50,
  all < 1.3, so per the pinned rule the analytic bands stand.
- **The six named checks:** (1) branch decoherence and (2) mixture model: V_sim ≥
  the re-frozen V-model at every M (the re-freeze, factor 0.75, is recorded in §5;
  the pre-freeze model was optimistic, caught by this check); (3) quasi-static vs
  Markovian T2*: quasi-static site disorder at σ_δ ≈ 3 kHz pulls the A0 slope by
  −0.0075 ± 0.0037 → the two-sided b_qs budget of §5; (4) branch-mixing/protection
  end-to-end: site-dependent ZZ + detuning shift the A0 slope by ≤ 1.6σ, consistent
  with the budgeted ζ² law (the exact statevector scan pinned
  bias = 0.00242·(ζ/3.8 kHz)²·(τ/1.2 µs)², quadratic to 4 digits, and its
  mechanism: the chiral map negates the Floquet gaps, so second-order shifts are
  opposite between branches; the law under the flown estimator is
  0.00257·(ζ/3.8 kHz)²·(τ/1.2 µs)², gate mode zz2_scan, with the same-sign
  positivity check); (5) readout correction: slope shift < 0.1σ, V recovered;
  (6) leakage/estimator bias noiseless: −0.0004 ≈ 0.1σ (negligible,
  gate-printed).
- **Instrument decisions this gate forced (all folded into §2/§4/§5):** idle
  padding removed (τ_step 1.2 → 0.7 µs; the padding fed the ζ² systematic
  quadratically and bought nothing), the V-model re-freeze, and the two
  clause-(a) budgets b_zz2, b_qs.

### 7b RECORDED (2026-07-15, runner `run_f129_ramsey_fringe.py --simulate`, external pipeline)

The runner's OWN circuit list (54 science, round-robin by M-point, + 2 cal
PUBs), its own transpilation, counts through the estimator and verdict IMPORTED
from the 7a gate (single code path, nothing re-implemented): **CONFIRMED on all
three 7b seeds** (dev_a +0.70 to +1.67, sep_b ≈ 62σ), band reconciliation
σ_a = 0.00339 = the 7a value exactly, budgets b_zz2 = 0.00427 / b_qs = 0.008
unchanged; the analytic bands stand. Runner `--certify` PASS: per-arm CX counts
match 112 + 2(cat−1) + 14M exactly, noiseless slopes within sampling noise of
the Floquet ΔΦ, and the conservatism ordering (A0 ≥ A1 ≥ A2 CX) holds. Billing:
64 PUBs × 16384 shots ≈ 1.05M shots → 5.5 min shots-ratio, 6.6-7.7 min with the
depth uplift (unchanged; re-checked against the account meter on flight day).
Zero-QPU dress rehearsal of the calibration path ran against the live backend
(2026-07-15, ibm_kingston): 41 chains pass the §8 v2.5 rule, chain file
written.

## 7. Cost, budget, and the gate order

54 science circuits × 16384 shots + ~10 cal PUBs ≈ 1.05M shots. Billing by
measurement, not estimate: the flown comparable (IBM_CONCENTRATOR_RELOADED: 376,832
shots billed 119 s) gives 5.5 QPU minutes by shots-ratio (4.7 science-only, the
design gate's printed number, plus the cal PUBs). Shots-ratio is known-optimistic
for deep circuits (the comparable itself billed 119 s against its own 85-99 s
projection), so the operative projection carries a ×1.2-1.4 depth uplift:
**6.6-7.7 QPU minutes**, leaving ≥ 2.3 minutes of cushion under the free
10 minutes/month; re-checked against the account meter before submission. The order,
verbatim discipline from the concentrator flight: design review rounds → fixes → 7a
recorded → runner built → 7b + certificates recorded → empty review of runner +
records → billing re-checked → THIS DOCUMENT COMMITTED as the pre-registration →
fresh calibration + chain selection (hard abort armed) → Tom's explicit go → ONE job
→ analysis → RECORD appended.

## 8. Chain rule and guards (carried in verbatim, extended)

Day-of hard aborts, no override flag (v2.5 rule; every threshold is an anchor of
a pinned model, not an aesthetic): an 8-qubit line with every qubit
T2echo ≥ 100 µs (the pinned 2× V-sensitivity absorbs a shortfall), the SEED-end
readout error ≤ 2 % (only the seed qubit is measured; the chain is oriented so
the better-readout end is the seed), median two-qubit (CZ, the Kingston native) error of the chain's bonds
≤ 0.5 % (the p₂ stress point where the §5 margins still hold 3.0× / 7.0×), and
the chain's median T2* estimate ≥ 70 µs (the visibility-model anchor, which
T2echo does not guarantee; where the calibration data carries no T2*, the
free-window M = 0 network Ramsey of the de-risking rule below stands in). The
v1-v2.4 rule inherited the concentrator's all-qubit RO ≤ 2 % and max/min
T2echo ≤ 2 verbatim; on an 8-qubit line that uniformity aesthetic excluded
every chain on a healthy device (0 of 950 on the 2026-07-15 ibm_kingston
calibration) while protecting nothing this design measures, and it was replaced
by the anchors above (41 of 950 pass the same day, best chain median 2q 0.15 %,
median T2* estimate 136 µs, seed RO 0.5 %).
Any failed fit, NaN, or guard trip is declared an instrument failure, never a
physics verdict. Between-arm drift is bounded by the round-robin
interleaving, not by the error bars (concentrator caveat carried in). The
good-window clock is not the readiness clock: a good calibration before the gates
are recorded is spent on free de-risking (characterization shots on the network
alone, M = 0), never on the flight.

## 9. Scope fence

- **Not a test of F130's q⁴ law.** That protocol is dead on amplitude
  ([F130_HW_INFEASIBILITY.md](F130_HW_INFEASIBILITY.md)); this flight measures a
  first-order phase, not a second-order coupling. (That document's "safe angle
  θ = 0.32" was the dead beat protocol's linewidth constraint; this flight has no
  engineered linewidth to respect, which is why θ = 0.5 is available here.)
- **The n = 12 lab does not fly under this pre-registration.** The equal-level pair
  (1,2,10) ~ (3,5,6) has a non-clean member and tests F130's hypothesis, not F129's
  law; its nearest impostor (dS = 0.018) is inseparable at this budget (design-gate
  appendix: 1.3σ). Any future n = 12 flight requires its own pre-registration with
  its own bands.
- **The n = 9 statement is F129's:** a clean-clean collision standing still against
  a winding clock, with the nearest permitted alternative rejected at 26.5σ.
  Precisely: the device directly measures the STEPPED chain's collision (the
  Floquet near-degeneracy with its computed θ³ drift); the continuum law's
  collision is inferred through the verified θ³ scaling, and cleanness is a
  construction property of the loaded triples, not a measured quantity.
- One device, one job, one lab; no cross-device claims.

## RECORD: the flight

**Flown 2026-07-15, ibm_kingston, job `d9br4vmg26ic73dgbgk0`, ONE Batch, 64 PUBs
× 16384 shots, billed 297 s ≈ 5.0 QPU min** (projection 6.6-7.7; under). Chain
[11, 12, 13, 14, 15, 19, 35, 34] (fresh §8 rule: min T2e 101 µs, med 2q 0.15 %,
med T2* est 136 µs, seed RO 0.5 %), submitted 18:17 after the pinned order
completed (pre-registration commit 5bf3b79, fresh calibration, dry-run aborts
green, Tom's explicit go). A local power failure killed the waiting process
mid-queue; the job was unaffected (0 QPU s at reconnect) and the results were
fetched by the submit-record job ID through the runner's --recover path, counts
persisted before any reduction
([`data/ibm_f129_ramsey_fringe_july2026/`](../data/ibm_f129_ramsey_fringe_july2026/),
the persisted counts + submit record; written by the external pipeline as
`results_f129/f129_flight_20260715_222348.json`).

**Measured (the committed §5 estimator on the persisted counts):**

| Arm | fitted slope (rad/step) | prediction | V(M=8) |
|-----|--------------------------|------------|--------|
| A0 | +0.0326 ± 0.0033 | center(θ̂) +0.0212 + budgets | 0.328 |
| A1 | −0.0928 ± 0.0034 | −0.0903 (0.7σ) | 0.459 |
| A2 | +0.3199 ± 0.0033 | +0.3214 (0.4σ) | 0.444 |

θ̂ = 0.4992 (the device hit the pinned angle to 0.16 %; clause (c) bound 0.024).

**Verdicts (the committed rules, applied by hand against the persisted counts;
the runner printout agreed):**
- **(a) PASS:** A0 deviation from the θ̂ center +0.0114, inside the budgeted
  window [−0.0182, +0.0225]. The excess is POSITIVE and sized within
  b_zz2 + b_qs, exactly the direction the ζ² anti-protection law predicts for
  same-sign always-on ZZ; purely statistically it is +3.36σ_a, which is why the
  budgets exist and were pre-registered.
- **(b) PASS:** A2 separates from A0 by 61.0σ combined, predicted sign.
- **(c) PASS:** A2 within 0.4σ of its nominal prediction.
- **STANDING FRINGE CONFIRMED.** The nearest impostor hypothesis (the smallest
  detuned triple in the census, predicted slope +0.109) lies ~22σ above the
  measured A0. F129's existence half has a hardware sighting: the clean mirror
  collision of the n = 9 comb stands still (to its computed θ³ drift and
  budgeted systematics) while the winding clock turns at 61σ and the
  opposite-winding dial A1 tracks its prediction at 0.7σ.
- **Instrument notes:** V(M=8) measured 0.33-0.46 across arms vs the
  conservative model 0.215 (the in-situ Ramsey PUBs put the chain's T2* well
  above the 70 µs anchor, consistent); σ_slope is therefore the frozen
  projection, conservative as pinned. No point fell below the 0.08 floor
  (min V = 0.328, all 27 points, 9/9 per arm); no PUB returned empty; billed
  usage under projection. Readout is asymmetric (P(err|0) = 0.02 %,
  P(err|1) = 0.84 %), handled by the 2×2 correction as pinned. Because the
  results were fetched through --recover after the power failure, the
  persisted snapshot_after matches snapshot_before (same calibration epoch);
  no independent post-flight drift snapshot exists, and between-arm drift is
  bounded by the round-robin interleaving (§8), not by the snapshots. The
  post-flight empty round recomputed every verdict number independently from
  the raw counts (clause (c) exact: 0.45σ) and the RECORD stands.

## Revision notes

- **v1 (2026-07-15):** initial design draft from the feasibility scout
  (θ = 0.32, even M-grid, 8192 shots, all arms verdict-bearing, ±0.030 ZZ band
  language mixed with exact-Wick numbers, single-arm σ).
- **v2 (2026-07-15, review round 1 folded; three empty lenses: physics from-below,
  spec audit, measurement statistics):** the physics lens returned CLEAN (mirror-ZZ
  zero, θ³ order, protocol exactness, bound robustness all recomputed
  independently). The statistics and spec lenses converged on one real defect: the
  v1 separations used single-arm σ and exact-Wick ZZ while the prose promised
  combined errors and a ±0.030 band, and the visibility model treated the 6-qubit
  branch cat as one collective DOF. Folded from below: the exact dephasing
  multiplier m_deph = 3.68 was computed and adopted (under it the v1 design fails
  its own 5σ bar: near controls ~4.3σ/3.7σ under the v2-era arm set and error
  model, superseded numbers); the design was
  re-powered (θ 0.32 → 0.5, integer M-grid, 16384 shots) and the verdict
  restructured so the discriminating clauses are A0's own impostor rejection
  (21.7σ, margin 4.3×) and the far control A3 (50.6σ, margin 10.1×), with A1/A2
  demoted to secondary dials. One ZZ treatment pinned (exact Wick, 5× transfer
  safety). The n = 12 lab was removed from flyable scope (unpinnable free arm).
  Added: predicted-V weights, per-point floor rule with arm-validity quorum, void
  mapping, cal/T1T2* PUBs, drift caveat, ECR day-of abort, null-band + bootstrap
  obligations on 7a/7b, billing by shots-ratio (~7 QPU min).
- **v2.1 (2026-07-15, review round 2 folded; three fresh lenses):** physics CLEAN
  and stronger than written: the mirror pair's protection was recomputed as the
  CHIRAL per-site/per-bond identity V_{n−k} = C·conj(V_k), which retires ZZ AND
  single-qubit Z-detuning even site-dependently; the spec lens's worry that
  site-dependent ζ breaks the A0 zero was resolved in the strong direction by that
  recomputation (per-bond machine-zero certified in the design gate), so clause (a)
  keeps its zero ZZ budget with a certificate instead of gaining a fudge band. The
  hop-angle exposure the chirality does not cover (d drift/dθ = 0.129) was closed
  by the pinned θ̂-from-A3 rule. Folded from the statistics lens: per-arm m_deph in
  the combined errors (A3 50.6σ → 50.2σ), the V-model linear-sensitivity statement
  (margins erode, verdicts survive 2×), and the first-moment/Jensen wording for
  m_deph. Folded from the spec lens: the n = 9 θ-exponent scan added to the gate
  (fitted 3.0054 over [0.02, 0.5]), the T2* day-of guard (T2echo does not imply
  T2* = 70 µs), the existence-half scoping sentence in the banner, the cal-PUB
  reconciliation of 6.2 vs 7.1 QPU minutes, and the θ = 0.32-vs-0.5 cross-note.
- **v2.2 (2026-07-15, review round 3 folded; three fresh lenses, physics CLEAN,
  spec and statistics MINOR/NIT-only):** pinning-precision pass. Clause (a)'s
  window center pinned as the FUNCTION slope_A0(θ̂) in the gate-emitted linear
  form (spec: the constant-vs-function ambiguity moved the center by up to
  1.07σ); the θ̂ inversion frozen as a gate artifact with its linearization error
  (9.0·10⁻⁵ ≈ 0.02σ); frozen-vs-realized σ stated per clause, with clause (c)
  explicitly against the NOMINAL θ = 0.5 prediction (the vacuous reading closed);
  the θ̂ variance propagated into σ_a (×1.023). From the statistics lens: the
  separation formula corrected to the true slope gap |slope_arm − slope_A0|
  (the old |·|−|·| form under-stated opposite-sign arms: the (2,5,6) dial
  8.0σ → 15.4σ), and
  the impostor search now minimizes the clause-(a) metric itself over all 54
  triples (same winner, robust method). From the spec lens: billing's shots-ratio
  carries the ×1.2-1.4 depth uplift the comparable itself showed, and the
  redundant second near-control (3,4,7) was dropped for the cushion (arms
  relabeled: the verdict control, formerly A3, is now A2; 72 → 54 science
  circuits, operative projection 6.6-7.7 QPU min, cushion ≥ 2.3 min). Scope
  fence: the stepped-device nuance and cleanness-as-construction made explicit.
- **v2.5 (2026-07-15, the runner + 7b session):** the hardware runner built in
  the external pipeline (modes --calibrate/--certify/--simulate/--hardware
  --yes/--analyze; counts persisted before reduction; repo gates imported as
  the single source). 7b RECORDED (§6): CONFIRMED on the runner's own
  transpiled circuits, bands identical to 7a. Two instrument findings folded:
  (i) the per-arm CX counts are 112 + 2(cat−1) + 14M (the smaller control cats
  spend fewer CX; the shared A0-formula V-model is thereby conservative, and
  the runner's --certify pins the exact per-arm counts); (ii) the §8 chain rule
  as inherited from the concentrator was structurally infeasible for an
  8-qubit line (0 of 950 chains on the live 2026-07-15 Kingston calibration,
  killed by the all-qubit RO and max/min T2echo uniformity clauses that
  protect nothing this design measures) and was replaced by the model-anchor
  rule of §8 (min T2echo 100 µs under the pinned 2× V-sensitivity, seed-end
  readout, median 2q ≤ 0.5 %, median T2* ≥ 70 µs), under which 41 chains pass
  the same day. Status: PRE-REGISTRATION; remaining before the flight are only
  the day-of calibration + hard abort and Tom's explicit go.
- **v2.4.1 (2026-07-15, 7a-landing review round 1 folded):** the spec lens's
  MAJOR was real: the budgeted clause-(a) window existed only in prose; the 7a
  gate's verdict() now implements the asymmetric budgeted window, the VIOLATED
  path and the dead band (exercised in-gate on synthetic slopes, both sides),
  and prints b_zz2, b_qs, zz_A2 (Wick, computed not hardcoded) and sigma_a from
  the emitted inflation. The zz2 statevector scan moved INTO the gate
  (mode zz2_scan); under the flown estimator weights the law reads 0.00257
  (estimator-defined; the 0.00242 of the padding-era weights is superseded),
  and b_zz2 rose 0.0033 → 0.0043 (1.22 estimator-pattern uplift, strictly
  conservative per the physics lens's recompute). The physics lens's substantive
  point folded: one-sidedness of b_zz2 is a property of the same-sign always-on
  regime (gate-checked over 40 positive draws), NOT of the chiral symmetry
  (which fixes the factor 2, not the sign); mixed-sign excursions ≤ 10⁻³ fall
  under b_qs. Staleness swept: dev_a envelope [−0.41, +1.29], 0.20σ per 1 %
  angle error, noiseless estimator bias gate-printed (−0.0004), sx-count
  provenance printed, compiler residual < 6·10⁻¹⁵.
- **v2.4 (2026-07-15, the 7a build session):** the Givens compiler exists and is
  certified (per-arm column permutations put each arm's cat block contiguous:
  GHZ-6/4/2 for A0/A1/A2; branch signs from the column-permutation parity tracked;
  circuit-vs-direct Slater machine-zero, full chain fringe = M·ΔΦ exact). 7a ran
  and is RECORDED in §6. What 7a changed: (i) the beat-era idle padding was
  dropped, τ_step = 0.7 µs (the ζ² systematic scales with τ², and the padding was
  concentrator legacy); (ii) the V-model gained the 0.75 one-qubit-gate factor
  (re-frozen, conservative, mixture-check enforced); (iii) clause (a) gained two
  pinned systematic budgets, b_zz2 = +0.0033 one-sided (the chiral mirror
  ANTI-protects at second order in ζ: mirrored Floquet gaps make the two branches'
  second-order shifts opposite; exact law 0.00242·(ζ/3.8 kHz)²·(τ/1.2 µs)², a
  physics find of this gate) and b_qs = ±0.008 two-sided (quasi-static site
  disorder at the T2*-consistent scale pulls the A0 slope via prep-basis
  mismatch); (iv) all §5 numbers re-pinned at the new timing (σ_slope 0.0033,
  impostor 26.5σ/5.3×, A1 20.7σ, A2 62.4σ/12.5×, V(M=8) = 0.215, θ̂ band 0.024,
  worst drop ×1.54). Remaining before the pre-registration commit: the 7b counts
  gate on the actual runner (external pipeline), billing re-check, and the final
  empty rounds.
- **v2.3 (2026-07-15, review round 4 folded; combined spec+statistics lens
  MINOR/NIT-only, physics lens with one F130 MAJOR):** traceability pass on this
  doc (every previously doc-side number is now gate-printed: the θ̂ band 0.0338
  derived instead of hardcoded, σ_a ×1.023, the worst quorum-drop inflation
  corrected ×1.45 → ×1.55, billing 5.5 incl. cal → 6.6-7.7 total). The physics
  lens named the one first-order channel at an exact degeneracy the diagonal
  argument does not touch, coherent branch mixing; it is identically zero by
  Slater-Condon (3-orbital difference vs ≤ 2-body operators), now stated in §2,
  certified in the gate (≤ 1.8·10⁻¹⁵), and added to the 7a checks; the chiral
  certificate was extended from NN bonds to all-range pairs (spectator ZZ
  crosstalk). The F130 companion doc's MAJOR: its "best working point q = 2" was
  the edge of a scan box, not a maximum (the optimistic bound grows monotonically
  in q); the argument was re-scoped to the coupling law's validity domain q ≤ 1
  (second-order perturbation theory), where the verdict STRENGTHENS to 360× short
  of 5σ, with q = 2 shown as a generous over-extrapolation that still fails 69×.
