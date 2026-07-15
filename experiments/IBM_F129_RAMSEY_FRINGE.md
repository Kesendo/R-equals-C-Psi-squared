# IBM F129: The Standing Fringe

**Status:** DESIGN DRAFT v2, pre-registration candidate. NOT flight-ready: the 7a/7b
gates are not recorded, the runner is not built. Per the one-shot discipline the
flight question exists only after: design reviews converge → 7a sim gate → runner +
7b counts gate → billing re-checked → this document committed as the pre-registration
→ fresh calibration + chain rule → Tom's explicit go.
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
are retired by the pair choice itself. The pair is resonant (S = 0); what the
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

**Instruments that do not exist yet (build list for the runner stage):** the Givens
network compiler (step eigenvectors → 2-qubit rotation schedule + inverse), the
GHZ-in-mode-register seed, the readout correction, and the two-quadrature fringe
estimator. The external pipeline (bonding cascade, ab_test Aer harness, calibration
pull) is reusable; a 2/3-magnon preparation has never been flown there.

## 4. Dose

There is no engineered dephasing dose in this experiment. The knobs are pinned:
θ = 0.5 per step, M-grid {0, 1, ..., 8}, 16384 shots per circuit, one Heron r2
device, all arms transpiled from one skeleton.

## 5. Statistics (pinned before the runner exists; verdict-determining numbers gate-printed)

**Estimator:** per (arm, M) the readout-corrected fringe phase φ = atan2(⟨Y⟩, ⟨X⟩)
and amplitude V_meas = √(⟨X⟩² + ⟨Y⟩²). Per arm, the phase is unwrapped over the
M-grid (worst predicted increment 0.32 rad/step × ΔM = 1 < π/2, ample margin) and
fitted by weighted least squares, φ(M) = slope·M + b. **Weights are computed from the
PREDICTED visibility model, not from V_meas** (removes weight-estimation bias);
V_meas is used only by the floor rule below. The primary numbers are the fitted
slopes; σ_slope is the weighted-LSQ error, and every two-arm comparison uses the
combined error √2·σ_slope.

**Error model (from below):** σ_φ = 1/(V·√shots) per point;
V(M) = (1−p₂)^CX(M) · e^(−3t/T1) · e^(−m_deph·t/T2*), with m_deph the exactly
computed FIRST MOMENT of the branch-pair Hamming distance (the two Slater branches'
site-configuration distributions; 3.68 for A0, each control arm uses its own pair's
value, ~3.8). The 6-qubit branch cat dephases ~3.7× faster than one qubit; modeling
it as one collective DOF would understate σ, and the single-exponential model built
on the first moment is conservative (Jensen: the true multi-exponential envelope
lies above it). Two-arm comparisons combine √(σ_A0² + σ_arm²). Projected at
p₂ = 0.3 %: V(M=8) = 0.155, σ_slope(A0) = 0.0041 rad/step. Sensitivity, pinned: a
uniform V-model miss scales σ_slope linearly (V/2 → σ ×2); it erodes the margins
below, not the verdicts, which survive a full 2× miss.

**ZZ and detuning systematic policy (one treatment, used everywhere):** the
collision arm carries ZERO budget for ZZ and single-qubit detuning, by the chiral
all-range/per-site cancellation of §2 plus the Slater-Condon branch-mixing zero
(exact for site-dependent ζ_ij, δ_i; both design-gate certificates at machine
precision). Control arms carry exact first-order Wick values
(one-body density matrices of the 3-magnon branches), floored at 5× their magnitude
as the device-transfer safety factor (ζ measured 3.6-3.9 kHz on ibm_marrakesh; 5×
covers ζ up to ~20 kHz, far beyond a flyable chain). The one per-step coherent
systematic the symmetry does NOT retire on A0 is hop-angle miscalibration:
d(drift)/dθ = 0.1294, so a 1 % angle error moves A0 by 0.16σ_slope (doc-side arithmetic on gate primitives); it is closed by
the θ̂ rule in clause (a).

**Predictions (θ = 0.5, from the design gate):**

| Arm | Floquet slope (rad/step) | ZZ (Wick, ζ = 4 kHz) | projected separation from A0 |
|-----|--------------------------|----------------------|------------------------------|
| A0 | +0.0213 (= 0.170·θ³) | 0 (exact, chiral) | |
| A1 | −0.0879 | −0.0041 | 15.4σ (secondary) |
| A2 | +0.3204 | +0.0017 | 50.2σ (verdict; margin 10.0× vs 5σ) |

Separations are the true slope gaps |slope_arm − slope_A0| (sign-safe; A1 winds
opposite to the A0 drift) in combined error (per-arm m_deph), floored by 5× the
Wick ZZ.

**Verdict rules (v2; the verdict is ONE conjunction, the secondaries carry no verdict
weight, so there is no multiple-comparison inflation):**

- **(a) The standing fringe:** A0's fitted slope lies within ±3σ_a of the window
  center slope_A0(θ̂). The center is a pinned FUNCTION, not a constant: in the
  linearized form emitted by the design gate,
  center = 0.0213 + 0.1294·(θ̂ − 0.5) rad/step (= +0.0213 at θ̂ = 0.5;
  linearization error ≤ 9.0·10⁻⁵ ≈ 0.02σ over the admissible band). The θ̂ rule
  (pinned, deterministic, closes the hop-angle exposure):
  θ̂ = 0.5 + (slope_A2_measured − 0.3204)/0.6164; A2 is 4.8× more θ-sensitive than
  A0, so it reads the effective angle out of the same job, and clause (c) bounds
  |θ̂ − 0.5| ≤ 0.034 whenever it holds. σ_a is the FROZEN design projection
  0.0041 rad/step inflated by the propagated θ̂ variance,
  σ_a = σ_slope·√(1 + (0.1294·σ_θ̂/σ_slope)²) = 1.023·σ_slope (gate-printed; A0 and A2 shots are
  independent, so there is no estimator correlation, only this variance term);
  re-frozen once at 7a/7b per the null-band rule below, never after data.
  Discrimination power: the nearest impostor by this clause's own metric (the
  minimum of |slope_c − slope_A0|/σ over all 54 distinct triples, design gate) is
  (3,4,7) at dS = 0.092, predicted slope +0.1092, landing 21.7σ from the window
  center (margin 4.3× vs 5σ), so clause (a) alone rejects every "accidental
  near-degeneracy" alternative F129 permits.
- **(b) The winding control:** A2's slope differs from A0's by ≥ 5σ combined
  (realized post-floor errors, per-arm m_deph), with the predicted sign.
  Projected 50.2σ, margin 10.0×.
- **(c) Validity:** A2's slope agrees with its NOMINAL θ = 0.5 prediction
  (+0.3204 + ZZ) within ±(3σ_slope + 5×|ZZ_Wick|), realized errors. Nominal, not
  θ̂-evaluated: comparing A2 to a prediction at the angle inferred from A2 itself
  would be vacuous; against the nominal angle, (c) is exactly what bounds
  |θ̂ − 0.5| ≤ 0.034. If (c) fails, the fringe clock is broken (a > 3.4 % angle
  miss or worse): the run is VOID (instrument), never a physics verdict.
- **STANDING FRINGE CONFIRMED** iff (a) ∧ (b) ∧ (c).
- **VIOLATED (collision winding)** iff A0's slope is outside ±5σ_a of the same
  θ̂-evaluated center while (c) holds. The 3-5σ gap between (a) and VIOLATED is a
  deliberate dead band: it reads inconclusive.
- Anything else: **instrument failure or inconclusive**, never a physics verdict.
- **Secondary dial A1:** prediction pinned above, result reported with the same
  estimator; it sharpens the story (the clock resolves dS = 0.092, winding the
  other way) but decides nothing.

**Point-floor rule (pinned):** a science point with V_meas < 0.08 is dropped
(predicted minimum 0.155, margin ~2×). An arm is valid iff ≥ 6 of its 9 points
survive, including ≥ 2 with M ≤ 2 and ≥ 2 with M ≥ 6. A0 or A2 invalid → run VOID;
A1 invalid → recorded, no verdict impact. Clause (a)'s σ_a stays the frozen
projection regardless of dropped points (the drop inflates only the realized errors
of (b)/(c), conservatively); a drop pattern that passes the quorum inflates the
realized σ_slope by at most ×1.55 (gate-printed worst case over all
quorum-passing survivor sets), which the (b) margin
10.0× absorbs.

**False-positive note:** CONFIRMED is a conjunction whose clause (a) window center
sits 21.7σ from the nearest alternative; the projected false-CONFIRM probability
under the design's own noise model is ≈ 0. The design's exposure is false-negative
(power), which the margins above cover, including at the p₂ stress point 0.5 %
(margins 3.0× / 7.0×, design gate) and under a uniform 2× V-model miss (margins
halve, clauses survive).

**Null band and bootstrap (owed to 7a/7b):** the 7a sim gate must produce (i) a
frozen null band for clause (a) from a zero-drift synthetic arm through this exact
estimator, and (ii) hierarchical bootstrap SEs; if the bootstrap SE exceeds 1.3×
the analytic σ_slope, the bands are re-frozen from the bootstrap at the 7b gate,
recorded in the Revision notes BEFORE the pre-registration commit, never after data.

## 6. The two gates (to be recorded here before any shot)

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
  magnon-number leakage feeding spurious phase. Not recorded yet.
- **7b (counts level):** the actual runner's circuits, the actual transpilation,
  counts through the actual estimator code; band reconciliation recorded. Hard
  pre-flight abort if 7b contradicts 7a. Not recorded yet.

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

Day-of hard aborts, no override flag: an 8-qubit line with every qubit
T2echo ≥ 150 µs, max/min T2echo ≤ 2, readout error ≤ 2 %, and median two-qubit (ECR)
error of the chain ≤ 0.5 % (the p₂ stress point where the §5 margins still hold
3.0× / 7.0×). The visibility model additionally rests on T2* = 70 µs, which T2echo
does not guarantee: the chain's calibration T2* (or, where the calibration data
carries none, the free-window M = 0 network Ramsey of the de-risking rule below)
must read ≥ 70 µs median, else abort; the §5 sensitivity line covers a residual
2× miss.
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
  a winding clock, with the nearest permitted alternative rejected at 21.7σ.
  Precisely: the device directly measures the STEPPED chain's collision (the
  Floquet near-degeneracy with its computed θ³ drift); the continuum law's
  collision is inferred through the verified θ³ scaling, and cleanness is a
  construction property of the loaded triples, not a measured quantity.
- One device, one job, one lab; no cross-device claims.

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
