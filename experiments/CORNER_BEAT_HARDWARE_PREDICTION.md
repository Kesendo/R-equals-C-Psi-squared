# The corner beat: a hardware prediction (pre-registration, DRAFT v7.2)

**Status: DRAFT v7.2, gate-build stage (Tom's go 2026-08-16). Design
converged over seven rounds; the gate went through its own empty review
(two rounds; all findings folded) and its v2.1 runs stand in §8a.
Verdict: provisionally flyable, CONDITIONAL ON FRACTIONAL-RZZ, at
the frozen 21-point grid with 16384 shots, ~22 QPU min (worst-end
power 2.9, P(detect) 0.92). The committed estimator is the
GAUGE-PINNED eigenchannel fit (§6); pinning the gauge dropped the
8192-shot detection probability from the pre-pin 0.985 to the honest
0.815, and Tom re-decided on the honest numbers (2026-08-16): **the
frozen configuration is the 21 × 3 grid at 16384 shots, ~22 QPU min:
verdict-grade worst-end P(detect) = 0.920 at H0 false rate 0.000
held-out, other bracket end 1.000 at 0.020.** Any flight remains
Tom's separate go.**

## What this is about

Take six spins in a row, all coupled equally, with exactly one spin flipped
up. That one flip can ride along the row in six different wave patterns,
each with its own pitch. Prepare two patterns at once and the row hums: the
site occupations beat at the difference of the two pitches. Three such
pattern pairs share the same difference, so the row offers three hums at
one frequency. Shine engineered noise on the row and the hums fade. This
experiment moves WHERE the noise shines without changing how much of it
there is in total. Spread evenly over all six sites, all three hums fade at
one single rate. Pushed to a corner, all of it on one site of each mirror
pair and none on the partner, the fading splits, and WHICH corner decides
by how much: half the corners give a wide split whose fastest rate sits at
the edge of what the proven containment allows, the other half a split
exactly 1/√3 as wide that stays clear of the edge. The machine is asked to
fly one corner of each kind against the even spread and show all three
numbers: the wide split, the narrow split, and the unchanged average.

## What the repo already holds

Swept 2026-08-16 (two store-sweep agents + the session's own reads), by
store:

- `experiments/THE_TWO_SPIN_ZEROS.md` (45e2f29) and
  `experiments/THE_MIRROR_TRANSVERSAL_CERTIFICATE.md` (e3dbab0): the theory
  this flight rides. The (1,1) ω = −2 room at N = 6, Δ = 1 is dim 3,
  breaking (no C_l protection on (1,1)), with gated w = 1/3, δ = 0,
  TIGHT at the UPPER edge of the centre interval
  (`simulations/results/locus_containment_gate.txt`, the
  `(1,1) om=-2 ... TIGHT(upper)` line). The gate prints bound = 1/3 on the
  interval [0, 2]; bound = min(s_max − m̄, m̄ − s_min) admits m̄ = 1/3 and
  m̄ = 5/3, and TIGHT(upper) selects m̄ = 5/3 (recomputed from below).
  w = max_σ ‖Σ_l σ_l C_l‖ (l ranging over MIRROR PAIRS in this sum) is
  an ADVERSARIAL MAXIMUM over the transversal sign choices; §2's two transversal classes are that definition made
  concrete, not a new result.
- `docs/ANALYTICAL_FORMULAS.md`: F154 (the locus saturation law; its own
  text still says the breaking interiors are "measured with slack", stale
  against e3dbab0's certificate, so this doc cites the certificate, not
  F154, for the containment); F140 (the frozen divisor: λ = −4γ̄ on the
  (0,1) corner block, a NEIGHBOUR of this experiment's 4γ̄, different
  block, do not conflate); F91 (the R₉₀ locus definition and its "not the
  decay rates of anything" fence, which the certificate answers); F122
  (the ceiling; the compression machinery's home, "high-Q degenerate PT");
  F67 (spectral encoding is the protecting mechanism; receiver territory).
- `experiments/` prior flights: `IBM_CONCENTRATOR_RELOADED.md` (the genre
  parent: randomized-RZ dephasing injection, M-binding channel, billing law
  0.283 ms/shot measured, sign/magnitude verdict split, the stale-band and
  pooled-binding traps, the transport-dressed trap, the shared-skeleton
  invariant, and the channel-realization ladder quoted verbatim in §5);
  `IBM_F129_RAMSEY_FRINGE.md` (compile the preparation from the ACTUAL
  Trotter step's modes; its C1/C2 certificates are free-fermion- and
  Slater-specific and do NOT transfer here, §5);
  `STAIRCASE_NULLTEST_HARDWARE_PREDICTION.md` (two B-BLOCK-INVALID
  flights: mid-batch T1 telegraphing invisible to calibration-reading
  guards is the recorded killer; frozen CLEAN conditions, family error
  rates, VOID precedence, the calibration-only day-of re-gate with its
  anti-circularity fence, the DD-disabled guard, the ζ convention pin
  ζ_c = ζ_shift/4, and per-shot billing 0.309-0.327 ms on delay-bearing
  jobs); `RECORD_PARITY_HARDWARE_PREDICTION.md` (Class 1 pre-submit vs
  Class 2 in-job aborts; the guard bank = pooled statistics PLUS
  individual floors, the retired "pooled bank" shorthand; the
  factor-collapse floors under ratio statistics; per-device band sets for
  twin backends; the no-runtime-mitigation requirement with the
  resilience flag secondary);
  `PRICE_PAIR_HARDWARE_PREDICTION.md` (the conditional-Ramsey ζ
  measurement this pipeline's ZZ numbers come from).
- `fw.Confirmations` / C# `ConfirmationsRegistry` (24 entries): entry 23
  `concentrator_site_contrast_kingston_july2026` (engineered per-site
  dephasing, rate observable, differential arms, ibm_kingston); entry 24
  `f129_standing_fringe_kingston_july2026` (Givens-prepared two-mode
  magnon superposition); `ibm_ep_onset_may2026` (per-site ⟨n_l⟩ from
  Z-basis counts only, no tomography). Nothing measuring a rate SPLIT
  under a moved profile has flown: this experiment is new.
- OpenArcs `compressed_density_laws` (Open): NextStep (3) closed both
  halves; open residue named there and inherited here: M0 diagonality and
  the 2160/2304 denominators measured, underived. `sideways_spin_ladder`:
  arithmetic kin (su(2) Clebsch-Gordan), different algebra, no overlap;
  its strip-retention precondition binds ITS device run, not this one.
- `docs/GLOSSARY.md`, `docs/CAUGHT_ERRORS.md`: no prior "corner beat"; the
  recorded dose-factor and convention error shapes feed §5 and §10.
- Local session scouts (gitignored, re-runnable, deliberately not named:
  a tracked document may not cite a gitignored path; every load-bearing
  number moves into the committed gate §8 before any flight). Numbers
  corrected across the rounds by recompute: v1's split fraction 96% and
  secular error 4% at Q = 10 (fit-window artifacts); v2's ω = −4 null
  (unpopulated dyads) and secular sign; v3's envelope-immunity claim,
  (1,2)-null form, M-ladder indexing, and C-only corridor (round 3).

## 1. The claim, one line

On a six-site uniform-J open XXZ chain at Δ = 1, sliding engineered
per-site dephasing from the uniform profile to a mirror-transversal corner
profile at the SAME total dose splits the one beat-decay rate of the
(1,1) ω = −2 room, (10/3)γ̄, into three, with a width set by the
transversal's CLASS: γ̄·{8/3, 10/3, 4} for the maximizing class (split
statistic s = (2/3)γ̄, top rate at the certified centre-interval edge) and
γ̄·{10/3 − 2/√27, 10/3, 10/3 + 2/√27} for the non-maximizing class
(s = (2/√27)γ̄, exactly 1/√3 of the maximizing width; the middle rate is
exactly 10/3 in both classes), the mean over the room exactly unchanged in
every case. What the arms TEST is the width law (s², both classes: the
C width against its band AND against the non-maximizing alternative
level, the ideal s²-ratio 3 entering through the dressed alternative
f_{C′}(B̂), §7 as restated in v6), and the mean anchor; the full rate
triples and
the edge value are reported context (§7), since the second-moment
statistic is blind to the triple's shape; the shape (one rate exactly at
the mean, the decoupled dyad) is reported through the CENTRED third
invariant (§6.4) with a frozen band, not as a CONFIRMED conjunct.

Wording fences:
- "Edge" means the edge of the size-class CENTRE interval
  [s_min, s_max]/2 · (−4γ̄) = [0, −4γ̄], NOT the edge of
  PROOF_CODIM1_BY_ADDITIVITY §6's cell-rate window, which is wider under a
  profile.
- The certificate's TIGHT(upper) attainment is an exact statement about
  the COMPRESSION; the finite-Q Liouvillian eigenvalue sits within ~0.1%
  of the edge on EITHER side (measured: +0.06% to +0.09% above at the
  continuous candidate points, ~0.1% below under the Strang step), so a
  dressed top rate slightly above 4γ̄ is expected finite-Q behaviour, not
  a violation. No verdict arm reads absolute attainment.

## 2. The objects

- Chain: N = 6, open, uniform J, Δ = 1 (Pauli book H = Σ_bonds XX+YY+ZZ).
- One-magnon modes k = 0..5 ascending, E = {−2.4641, −1, 1, 3, 4.4641, 5}.
- The room: ω = −2 dyads (1,2), (2,3), (3,5) spanning the (1,1) block's
  dim-3 mixed eigenspace (modes 1,2,3,5 at E = −1, 1, 3, 5).
- Notation: mirror SITE pairs are written {a,b} ({0,5}, {1,4}, {2,3});
  mode DYADS are written (i,j) ((1,2), (2,3), (3,5)); the collision on
  "(2,3)"-shaped strings is real, so the brackets carry the type.
- Profiles on the R₉₀ locus (γ_l + γ_{5−l} = 2γ̄, total dose 6γ̄ in every
  engineered arm):
  - UNIFORM: γ_l = γ̄ on every site. Compressed matrix EXACTLY
    −(10/3)γ̄·Id (scalar; rests on Σ_l ψ_i(l)²ψ_j(l)² = 1/6 exactly FOR
    THE THREE ROOM DYADS (general mode pairs differ, e.g. 1/12), an
    exact route per the house no-rounding rule, compared exactly in G1).
  - CORNER: a mirror transversal, one site of each pair
    {0,5}, {1,4}, {2,3} at 2γ̄, partner dark. **The eight transversals
    split into two classes** (recomputed from below in all three rounds):
    - MAXIMIZING (σ_outer ≠ σ_middle: the lit sites of pairs {0,5} and
      {2,3} on OPPOSITE ends): {0,1,3}, {0,3,4}, {1,2,5}, {2,4,5} →
      rates γ̄·{8/3, 10/3, 4}, s = (2/3)γ̄, attains w = 1/3. The split
      lives in the (2,3)↔(3,5) coupling; dyad (1,2) is EXACTLY decoupled
      and keeps the uniform rate.
    - NON-MAXIMIZING (σ_outer = σ_middle): {0,1,2}, {0,2,4}, {1,3,5},
      {3,4,5} → the rate triple γ̄·{10/3 − 2/√27, 10/3, 10/3 + 2/√27}
      (the MIDDLE rate is exactly 10/3: here the DECOUPLED dyad is
      (2,3)), s = (2/√27)γ̄ ≈ 0.385γ̄ (rate
      spread = 4γ̄ × the NORM of the transversal's odd part, in the
      pair-odd convention C_p = comp(N_low − N_high) per mirror pair p
      with comp the certificate's NUMBER-OPERATOR compression (the
      dissipator's dose factor lives entirely in the 4γ̄ prefactor;
      building C_p from the rate compression double-counts it),
      under which C_1 ≡ 0 on this room, the maximizing combination is
      the DIFFERENCE with norm 1/3 = w and the non-maximizing the SUM
      with norm 1/√27 (round-5 repair of a swapped sign); the numeric
      kinship with the census's resC = 1/√27 for (1,1) is an observation,
      not a provenance: the census norm is a single C_l, this one a
      two-pair combination), top rate 7% inside the edge. Here
      the coupling sits on (1,2)↔(3,5): the (1,2) dyad mixes, and WHERE
      the coupling sits is itself class physics.
    The pair {1,4} is free in both classes (C_1 ≡ 0 on this room).
  - Pinned profiles: **C = A = {0,3,4}** (maximizing) and
    **C′ = A′ = {0,1,2}** (non-maximizing). Class membership of the
    REALIZED profiles on the PHYSICAL qubits is a post-transpile hard
    assertion (§9).
- The conserved mean is the compressed trace: Tr = −10γ̄ for the uniform
  profile and all eight transversals. Off the locus at equal dose the
  trace ranges over [−10.67, −8.67]γ̄ (single-site vertices; −6.7%/+13.3%
  around −10γ̄), so the mean anchors the locus but does not confirm it
  alone.
- The two rooms are anti-correlated across classes: Σσ_l C_l vanishes on
  the ω = −4 room for exactly the four ω = −2-maximizing transversals.
  The ω = −4 room is NOT read by this experiment (no pinned preparation
  populates its dyads (1,3)/(2,5); the v2 "free null" was refuted).

## 3. Arms

Four shot-bearing arms plus one in-circuit null. All four arms share ONE
transpiled skeleton, differing only in the bound phase tables (φ = 0 for
N0; the concentrator's shared-skeleton invariant, asserted in code), so
they are gate- and duration-identical by construction. Interleave order is
DETERMINISTIC and pinned: within each depth, N0 → U → C′ → C, depths
ascending. The honest rationale (round-4 repair of a leftover pre-s²
argument): under the verdict statistic s², UNIFORM monotone drift has
exactly zero effect (shift invariance), and SCATTERED drift is unsigned
(G4 budgets it symmetrically), so no order can be claimed
drift-conservative; the order is pinned for determinism (the
record-parity requirement), and the W pair C′/C is placed adjacent
because W carries the two-class law, the design's sharpest comparison.

- **U**: uniform engineered profile γ̄ per site.
- **C**: corner profile A = {0,3,4} (maximizing class).
- **C′**: corner profile A′ = {0,1,2} (non-maximizing class), same total
  dose, same depth grid: the live alternative hypothesis that makes the
  W-arm a discrimination between two from-below predictions.
- **N0**: no engineered dephasing (φ = 0 bindings). N0's fitted FULL 3×3
  room generator B̂ (not a profile, not a scalar: the site-profile
  inverse is rank-deficient, 5 of 6, and the real background need not be
  of dephasing form) enters the dressed predictions (§7). Because s² is
  shift-invariant, **s²_dressed(U) = s²(B̂) up to the §4 secular term
  (exact if the U engineered part were exactly scalar): the U arm's
  predicted split IS the background's own split**, so the U band is set
  by N0's fit precision, not by U's shots; var(B̂) propagates into every
  band (G1), and N0 carries its own CLEAN conditions (§7).
- **In-circuit null**: the decoupled (1,2) dyad, read as a
  **C-minus-U generator difference**: [Ĝ(C) − Ĝ(U)] must have (1,2) row
  and column consistent with zero (exact in the model for ANY common
  background B̂, which cancels in the difference; stated as a difference
  because under generic B̂ neither arm's own (1,2) row is zero and the
  (1,2) direction is not an eigenvector of either single-arm generator;
  caught in round 3). Design asset: sites 1 and 4 read dyad (1,2) nearly free of the
  other two dyads (rows (±0.2357, 0, 0) EXACTLY in the continuous-H
  basis; in the flown Strang basis the zeros become ±0.0023, 1%
  contamination, which the G2 margin carries). The DIFFERENCE reading is a wiring/class-membership check and
  may VOID (threshold G2); the (1,2) RATE itself is class physics and is
  REPORTED, never a VOID trigger (a genuine falsification must not be
  recorded as instrument failure). In C′ the (1,2) dyad mixes by class
  physics, but C′ has its OWN exactly decoupled dyad, (2,3): the
  analogous C′-minus-U difference null on the (2,3) row/column applies
  there (round-5 addition: the only in-job wiring check behind C′'s
  Class-1 class assertion), with the caveat that no site reads (2,3)
  contamination-free (rows 1 and 4 are zero on it), so it needs the
  full three-dyad inversion, unlike the (1,2) null.

## 4. The working point and the corridor (two factors, opposite directions, per arm)

The corridor is the PRODUCT of two effects (measured estimator-free in
rounds 1-3, agreeing on mechanism and magnitude; the committed gate
re-establishes the numbers through the committed estimator):

- **Secular factor, per arm and opposite in sign:** the finite-Q
  compression error moves the two arms in OPPOSITE directions: it
  INCREASES s(C) (+0.60% at Q = 20, +1.72% at Q = 12, +2.51% at Q = 10;
  re-signed in round 2) and DECREASES s(C′) by roughly twice as much
  (−1.3% at Q = 20, −3.5% at Q = 12, −5.1% at Q = 10; caught in round 4,
  where the v3 bullet was still C-only). The per-rate error is smaller
  (≤ 1.1% down to Q = 10), sits mostly on the lowest rate, and is NOT
  the band-setting quantity; the split's error is.
- **Trotter factor, per arm:** the Strang step detunes dyad (3,5) by
  δ_T = 2.0·(J·dt)², an ANGULAR FREQUENCY in units of J, commensurate
  with the rates (0.0201 at J·dt = 0.10, 0.0453 at 0.15; renamed from
  Δ in round 6: Δ is the XXZ anisotropy alone, and a reader taking the
  detuning as a per-step phase concludes falsely that the candidate is
  past the EP). For C the
  detuning competes with the coupling c = (2/3)γ̄ in its 2×2 block:
  half-split √(c² − (Δ/2)²), a factor ≤ 1 that reproduces C's rows
  within the method spread. For C′ (c′ = c/√3, and the coupling on a DIFFERENT dyad pair)
  the 2×2 model holds only at Q ≥ 15; at the candidate (J·dt = 0.15,
  Q = 10) and the retired point the C′ Trotter factor is measured
  ABOVE 1 (+8.5% and +6.5% respectively). The table is the authority
  for the ROOM SPECTRUM; the committed estimator adds its own
  measured dressing on top (pinned-gauge, noiseless, the committed
  gate's frozen mode: **s² = 1.21× bare ideal for C, 0.54× for C′,
  ~0 for U**; from the class-basis channel projection and a
  prep-dependent secular population leak, deterministic, carried by
  f_C/f_{C′}; §6),
  so the model is C's mechanism, not a per-arm law, and the TABLE below
  (direct Floquet-Lindblad) is the authority. **Each arm has its own
  corridor and its own exceptional point**, C′'s EP at markedly lower
  Q; beyond an EP that arm's split VANISHES. Higher Q is NOT safer
  here.

Full Floquet-Lindblad values, in s (NOT s²; the ratio column is
context only since v6 retired the ratio arm):

| J·dt | Q | s(C)/ideal | s(C′)/ideal | (s(C)/s(C′))/√3 |
|------|----|-----------|-------------|------------------|
| 0.10 | 12 (retired) | 1.008 | 1.027 | 0.982 |
| 0.15 | 10 (CANDIDATE) | 0.984 | 1.030 | 0.955 |
| 0.15 | 15 | 0.893 | 0.816 | 1.094 |
| 0.15 | 20 | 0.766 | 0.274 | 2.80 (C′ near EP) |

EP locations at J·dt = 0.15 (the candidate's row): C at Q ≈ 30 exact
(2×2 model 29.5, conservative by ~1 unit in Q, not more: Q = 32 is
already collapsed); **C′ at Q ≈ 21, the BINDING margin at the flown
Q = 10 is C′'s, a factor 2**. (At the retired J·dt = 0.10 point the
margins were Q ≈ 67 / ≈ 38.) All bands come from G1's counts-level sim,
never from this table alone (§7: the post-selection envelope biases
absolute s² and no bare-corridor band survives it).

Method spread: compressed-map and exact-eigenvalue readings differ by
2-4%; the committed gate's counts-level estimator is the arbiter.

Candidate working point (chosen by G1 v2's two-stage selection on the
VERDICT metric P(d > θ_D), held-out seed, both f_leak ends; §8a): **Q =
10, J·dt = 0.15, grid 21 points × 3 steps (deep end 60), 16384 total
shots per (arm, depth, prep) (Tom's second budget decision, on the
post-pin numbers), FRACTIONAL-RZZ gates (2 two-qubit gates
per XXZ block; the machinery exists in the pipeline's record-parity
scripts, use_fractional_gates=True with the CZ = NO-FLIGHT
assertion)**: verdict-grade at ~22 QPU min under the GAUGE-PINNED committed
estimator: **worst-end P(detect) = 0.920 (power 2.9, H0 false 0.000
held-out), other bracket end P(detect) = 1.000** (the 8192-shot
variant: P(detect) 0.815 at ~11 min; the pre-pin 0.985 was
gauge-favored and is withdrawn). **Grid AND shots
FROZEN by Tom's two budget decisions (2026-08-16), the second taken
on the honest post-pin numbers.** The v4 fear that the
deep-end exponent kills the flight was WRONG in mechanism: the weighted
estimator devalues the dying deep points instead of dying with them.
Standard gates (3 per block) reach worst-end power only ~1.7 at 11
min and ~1.4 at 22 min under the pinned estimator: fractional-RZZ is
LOAD-BEARING and its availability on the
flown backend is a Class-1 guard, as is the calibrated 2q error on the
used edges (p2 ≤ 0.6%; the measured p2 ladder is §8a's). The former candidate Q = 12, J·dt = 0.10
underperforms at every tested budget and is retired.

**J-calibration sensitivity** (G5 v2.1, estimator level, nominal-J
preparation basis, pinned gauge, the flight's situation): the response
is EVEN, ∂s²/s² = +4.5% for both ±5% J (a one-sided inflation bias;
§8a); bond-J scatter σ_J = 1% costs ≤ 10% of s²(C) (95%, frozen-draw
systematic), σ_J = 2% ≤ 25%. The band carries these terms (G5). J is VERIFIED IN-ANALYSIS from
the science arms' own fitted dyad frequencies (the Floquet frequencies
are J-sensitive and fitted anyway; estimator frozen by G5; no extra
PUBs), a precondition on D-mag.

## 5. Protocol

**Backend and line.** ibm_kingston or its pinned twin ibm_marrakesh
(Heron r2). Six-qubit line by the concentrator uniform-line rule (§9,
Class 1). **Per-device band sets:** G1 freezes a band set per twin; the
flown device's set governs (the record-parity twin rule). BOTH twins
must expose fractional rzz for the twin rule to hold; the pre-gate
properties query checks it on both, and a twin without it drops out
of the pinned pair (no silent standard-gate fallback). **Runtime
mitigation OFF:** the pinned requirement is NO runtime mitigation and
no scheduler dynamical decoupling (DD would rewrite the very dephasing
profile under test; offline CAL inversion would double-correct);
per the record-parity wording the resilience flag is evidence, not the
assertion (it is not guaranteed a Sampler knob), so the runner asserts
the requirement on the transpiled circuits, and the flag is set where
exposed.

**Preparations.** Three, one per room dyad: (ψ_i + ψ_j)/√2 is itself a
SINGLE-PARTICLE orbital, prepared by N − 1 = 5 Givens rotations from
|100000⟩ (no cat block; v1's cat was F129's three-particle structure
mis-imported). Orbitals from the eigenvectors of the ACTUAL Strang step
(the f129 lesson that does transfer). The f129 C1/C2 certificates do NOT
transfer (C2 is a Slater-Condon zero needing 3 differing orbitals; C1 is
the pure-hopping chiral relation, absent at Δ = 1). Replacement
certificates, committed in the gate: prep fidelity of the compiled
circuit against the Strang-step eigenvector, and a leakage bound out of
the one-magnon sector under the interacting step.

**Time evolution.** Second-order (Strang) odd/even Trotter at the frozen
J·dt. First-order is excluded by measurement (spurious uniform-arm
half-split 13-25% of ideal across the grid vs Strang's ~1-2%). The final
injected-RZ layer is a Z-basis no-op (realized dose at depth n is n−1
layers), absorbed by the estimator's free amplitude. The depth grid is
UNIFORM in Trotter steps: 21 points at spacing k = 3 steps, starting at
step 0 (the realized dose at depth n is n−1 injection layers) and
ending at step 60 (§4/§8a; frozen by Tom's budget decision). The
committed estimator (§6) fits on that grid; the per-grid-interval room
phase is |ω|·k·dt ≈ 0.9 rad ≪ π (frequency identification safe). No
grid point is dropped after data exists (the concentrator rule: failing
points are dropped before pre-registration, never after).

**Dephasing injection.** The concentrator recipe verbatim, convention
pinned: RZ(φ) = exp(−iφZ/2); per step, per site, φ ~ N(0, σ_l²),
σ_l = √(4·γ_l·Δt_step), giving per-step retention e^{−2γ_l·Δt} per site
and e^{−2(γ_l+γ_m)·Δt} on a two-site coherence (recomputed, rounds 1-3).
**Channel-realization ladder, quoted verbatim from the parent flight**
(`IBM_CONCENTRATOR_RELOADED.md`: "single-realization slope gap 0.0208
(28.3% of the channel slope) ... the single-realization gap shrinks with
M (0.0208 → 0.0032 → 0.0010 at M = 256/1024/4096)", recorded there as an
interpretive decomposition, not a gate number): so ~28% at M = 256,
~4.4% at M = 1024, ~1.4% at M = 4096, on THAT flight's slope observable.
**M = 1024 candidate** (realization error ~4.4% against a ±20% verdict
quantity, a real but affordable term); G3 freezes M from the realization
error measured on s² itself, and the §10 cost band re-freezes with it
(certificates computed at a provisional M never carry; note the
concentrator's flown FINAL was M = 256, and its billing law showed no
per-binding term, but §10 prices the per-binding overhead question to
G3 anyway). Shots are TOTAL per (arm, depth, prep), split across the M
bindings (16384/1024 = 16 per binding at the frozen numbers; the
concentrator flew 32/binding); per-binding counts persisted (bypassing
get_counts() pooling, the concentrator's recorded instrument deviation).
Exactly one injected RZ per site per step, asserted post-transpile.

**Readout.** Z-basis counts only, per-site ⟨n_l(t)⟩. Per-site resolution
is MANDATORY: the three dyads' site-amplitude columns sum to exactly
zero (total magnetization carries no beat). Sites 1 and 4 are exact
nodes of mode 3 and blind to dyads (2,3) and (3,5): four of six sites
carry the split (and exactly those two sites carry the clean (1,2)-null
readout, §3). CAL0/CAL1 confusion PUBs fly in the same job; the
confusion inversion is applied to per-binding counts before
post-selection (linear, so algebraically order-free until the fit; G2
owns the point where the chain becomes nonlinear), and G2 must show the
post-selection identity survives asymmetric per-qubit confusion.
Post-selection on total excitation = 1 is pinned in: it removes the
LEAKAGE component of gate error exactly. The within-sector component
survives (7 of 15 two-qubit Pauli errors preserve the excitation count
ON a differing pair: IZ/ZI/ZZ outright, XX/XY/YX/YY on the one-magnon
pair subspace only) and is a G1 line item. The residual post-selection
envelope (1−p)/(1−0.90625p) is common to all dyads and arms but
NON-EXPONENTIAL: it biases ABSOLUTE s² by up to ~−20% at deep grids
(measured, round 3) while cancelling in the s²-ratio (to < 0.1%) and in
r̄ differences; which is why D-mag's band comes from G1's counts-level
sim and never from §4's envelope-free corridor (§7).

## 6. Committed estimator (frozen in structure; constants from the gate)

1. Counts → confusion-inverted per-binding counts → post-selected
   one-magnon site occupations n_l(t).
2. Site→dyad inversion: the 6×3 map A[l, d] = ψ_i(l)ψ_j(l) has exactly
   orthogonal columns (Gram = (1/6)·Id, cond = 1.000); fitted WITH an
   intercept, which does not degrade the conditioning (the intercept
   column is orthogonal to all three dyad columns: Σ_l ψ_i(l)ψ_j(l) = 0
   for i ≠ j).
3. THE COMMITTED FIT (amended in v7; the formerly pinned complex
   one-step propagator fit was withdrawn for measured cause, recorded
   in §14: blind pole retrieval cannot separate three near-degenerate
   frequencies over ~2 periods, and quadrature demodulation is biased
   by the decaying envelope at the 8%-per-sample level, both measured
   in the gate's estimator trials): the EIGENCHANNEL damped-cosine
   fit. Per arm, the class-predicted eigenchannels (C: y₁ decoupled,
   v± = y₂ ± y₃; C′: y₂ decoupled, u± = y₁ ± y₃; U: y₁, y₂, y₃) are
   each fit as y(t) = b + a·e^{−rt}·cos((ω+δω)t + φ), linear in
   (b, a·cosφ, a·sinφ) at fixed (r, δω), a bounded 2-D grid search
   with refinement over (r, δω), traces weighted by kept counts
   (√w in the design matrix), joint over the preparations that carry
   the channel. The estimator ENCODES the arm's class hypothesis
   (legitimate: the profile is chosen, not measured, and a wrong basis
   MIXES rates and shrinks the fitted split, the conservative
   direction); its deterministic channel bias (~±0.2γ̄ on the ±
   channels) is part of the dressed centres per the corridor-dressing
   rule.
4. **s² ≡ ((r₁−r₂)² + (r₁−r₃)² + (r₂−r₃)²)/6** on the three fitted
   channel rates: sort-free by construction (symmetric), and
   NON-NEGATIVE by construction (the v5 negative-s² machinery is
   withdrawn with the estimator that produced it; under H0 the
   statistic has a POSITIVE noise floor, measured by G1's H0 runs and
   absorbed into θ_D, which is why θ_D is frozen on the H0
   distribution of the DIFFERENCE and not at zero). Mean rate
   r̄ = (r₁+r₂+r₃)/3. The shape statistic Π(r_i − r̄) (shift-
   invariant, predicted exactly 0 for both classes' symmetric
   triples) is REPORTED with a frozen band, not a CONFIRMED conjunct.
   The fitted δω per channel is the detuning readout feeding the G5
   J-estimator.
5. Rates are invariant under any fixed LINEAR mis-calibration of the
   site→dyad map (similarity transform); the affine part is handled by
   the intercept in step 2.
6. Resampling unit for all error bars: the BINDING, resampled jointly
   across depths and preparations within an arm. s² resamples cleanly
   (no folding); the sorted rates do not and are presentation only.

## 7. Verdict rule (structure frozen; bands, thresholds, margins, and family rates from the gate)

Arm types are heterogeneous and each carries its own decision rule
(round-3 repair: the staircase's HOLDS/MARGINAL/VIOLATED grading is a
NULL-statistic shape and applies only where a null is tested):

- **D-sign, one-sided detection.** s²(C) − s²(U) > θ_D (θ_D frozen
  by G1 ON the H0 distribution of this difference, both arms carry a
  positive noise floor and the s²(B̂) offset, which cancel in
  expectation but not in variance; the detection is on the difference,
  never against an assumed-zero baseline; the frozen-config
  measurement: θ_D = 0.0025, held-out H0 false rate 0.000 at the
  worst f_leak end, 0.020 at the other; §8a).
- **W, the class law as a LEVEL discrimination (restated in v6: the
  ratio form failed its own G1 reachability check; re-measured in
  v7.1 with the committed gauge-pinned estimator).** The gate measures
  s²(C′) BELOW the U arm at budget (s²(C′) − s²(U) = −0.0015: the
  hop-noise floor of the degenerate U arm exceeds the whole C′
  prediction), so neither a ratio nor a C′-detection arm is freezable.
  The class law lives in the exceedance: (i) s²(C) − s²(C′) > θ_W
  (one-sided, θ_W frozen by G1; measured at the frozen config
  d_W = +0.0063 ± 0.0013, power 4.7 at the worst f_leak end); The C′
  containment itself lives ONLY in D-mag (v7 repair: a
  duplicated band in two arms broke the branch partition and the
  registration rationale); W's verdict content is the one-sided
  exceedance alone, which keeps W a non-magnitude arm and the
  registration scope (D-sign + W) coherent. Both centres are frozen
  functions of B̂ (f_C, f_{C′}, committed as code). The s²(C) floor
  and its frozen adjudication stay (a floor trip with dose
  certificates passing and N0-CLEAN holding is ¬W class anomalous,
  else VOID); the former s²(C′) floor is retired (a low s²(C′) is
  expected, not pathological).
- **D-mag, two-sided agreement, BOTH corner magnitudes.** s²(C)
  inside its frozen band around f_C(B̂) AND s²(C′) inside its band
  around f_{C′}(B̂). s² is NOT additive in the generator
  (s²(Ĝ + B̂) carries cross terms unless B̂ is scalar), so both
  CENTRES are frozen FUNCTIONS committed as code with the same status
  as §9's width scaling functions, evaluated on the measured B̂
  (round-4 repair: a frozen NUMBER cannot exist pre-flight). The
  functions include the post-selection envelope from G1's counts-level
  sim (which biases absolute s² up to ~−20% at deep grids); a
  bare-corridor band would fire VIOLATED on a perfect device. G1
  verifies the reachability of the "D-sign ∧ W ∧ ¬D-mag" branch (the
  W exceedance and the D-mag band both read s²(C); the branch exists
  only if the band widths admit it). PRECONDITION: the G5 in-analysis J
  estimator inside its frozen pass band; a J-band failure is a VOID
  trigger (G5), not a silent widening.
- **A, equivalence test (the anchor; renamed from M in round 7
  to end the collision with the binding count M, G3's object).** The
  CIs on r̄(C) − r̄(U) and
  r̄(C′) − r̄(U) are CONTAINED in frozen equivalence margins (a TOST
  shape: "consistent with zero" alone would be confirmation by
  imprecision). Exact-in-the-compressed-model under any common additive
  background; the model's own secular error is §4's 1-2% and the margin
  carries it.
- **In-circuit nulls, equivalence tests**, one per corner arm: the
  (1,2) row/column of Ĝ(C) − Ĝ(U), and the (2,3) row/column of
  Ĝ(C′) − Ĝ(U) (§3; C′'s wiring check, needing the full three-dyad
  inversion, so G2 freezes its margin SEPARATELY from the
  contamination-free (1,2) case). Both VOID-capable (wiring), while
  the decoupled-dyad RATES are report-only (physics).

**N0-CLEAN (frozen by G2; failing it is VOID):** B̂ physical (spectrum
in the admissible half-plane, magnitude envelope), fit residual bound,
and first-half vs second-half batch consistency of B̂ (the staircase
bracket shape). N0 is verdict-input data consumed ONLY through the
frozen committed functions (the dressings and f_C); no analyst-side
adjustment ever; the day-of re-gate (§9) is calibration-only and N0 is
not among its inputs.

**VOID (Class 2, §9).** Triggers, each with its freezing owner: N0-CLEAN
(G2), T1 CLEAN conditions (G4), dose certificates (G3), the (1,2) and (2,3) difference-nulls moved (G2,
separate margins), fit health: channel-fit residual or frequency-offset excursion (G2),
the s²(C) floor under W (G1), band-validity window left (G1). The bank
combines POOLED statistics with INDIVIDUAL floors (the record-parity
correction of its own earlier shorthand); with the G5 J pass band added the bank holds NINE triggers and the
union-bound false-VOID accounting is measured JOINTLY by G1 (the
verdict arms are built from the same FOUR fits, N0's B̂ entering every
dressed prediction, and are strongly correlated: family rates are
measured, never multiplied).

**Verdicts.** ALL arms are always computed (round-4 repair: a
short-circuit would make every A-dependent branch below unreachable);
the order VOID → D-sign → W → D-mag → A is REPORTING precedence only,
naming which failure leads the verdict when several co-occur:

- CONFIRMED = D-sign ∧ W ∧ D-mag ∧ A.
- "Split confirmed, width off-prediction" = D-sign ∧ W ∧ ¬D-mag ∧ A.
- "Split confirmed, class anomalous" = D-sign ∧ ¬W ∧ A (follow-up: the
  realized-profile assertion records which classes flew).
- "Split observed, anchor failed" = D-sign ∧ ¬A (any W/D-mag);
  the symmetric ¬D-sign ∧ ¬A is named "anchor failed, no split
  observed" and likewise reports the failing anchor:
  instrument-suspect INCONCLUSIVE; follow-up: the dose or the locus
  condition did not realize; the dose certificates and N0 record decide
  which.
- ¬D-sign partitions (disjoint by construction, all three under the
  same preconditions: guards clean, A holds, and the power condition
  SE[s²(C) − s²(U)] ≤ κ × G1-projected, κ frozen, on that NAMED
  statistic): **FALSIFIED** = the CI on s²(C) − s²(U) is contained in
  [−θ_F, θ_F], where θ_F is a separate PHYSICS-SIZED equivalence
  margin frozen by G1 as a fraction of the predicted difference
  (round-5 repair: containment in [−θ_D, θ_D] is scale-invariantly
  impossible when θ_D is the detection critical value ~1.6-2 SE
  against a ~3.9-SE-wide CI; G1 additionally verifies θ_F's
  reachability at the frozen shot budget); **Anti-D** = the CI lies
  entirely below −θ_F, recorded
  as falsified-with-inverted-sign; **INCONCLUSIVE-underpowered / indeterminate** = the
  power condition fails, or ANY other ¬D-sign CI position (straddling
  ±θ_F, or lying entirely above θ_F, reachable when θ_F < θ_D) (and any ¬D-sign case whose preconditions fail reports
  the failing guard/anchor instead).
- Any unlisted pattern: INCONCLUSIVE, named, with follow-up stated.

**On CONFIRMED:** Confirmations registration scoped to D-sign + W only
(the concentrator precedent: magnitude arms are not registered); the
entry records that D-mag and A held as CONFIRMED's other conjuncts and
that W's predicted centre was the B̂-dressed function, since the
registry stores predicted-vs-measured.

## 8. Gate items before bands freeze (the committed gate: `simulations/corner_beat_gate.py`)

- G1: counts-level end-to-end sim of the flown construction (shots,
  readout confusion, gate-error as depolarizing per 2q gate including
  its within-sector component, the post-selection envelope, M-finite
  channel realization, var(B̂) propagation) through the §6 estimator;
  joint choice of (Q, J·dt, depth grid endpoint, shots, and the depth
  question, answered by v2: which grid survives the deep-end pricing);
  freezes: all bands, θ_D, θ_F (with its reachability verified at the
  frozen shot budget, §7), the "W ∧ ¬D-mag" branch reachability, the
  s²(C) floor, the A-arm equivalence margins (both layers: the statistical
  part and the §4 secular physics term the margin carries; the null
  margins are G2's), the projected SE of s²(C) − s²(U) at the frozen
  budget (the κ condition's reference), the frozen
  centre-functions f_C(B̂) AND f_{C′}(B̂) and the dressing
  compositions, κ, the family error rates (P(≥1 VOID),
  P(false CONFIRMED) under H0, measured jointly), the frozen
  depth-point list, the band-validity window (dimensions and its
  relation to the Class-1 line rule), the day-of band-WIDTH scaling
  functions (§9's re-gate has no other source), the billing cap, and
  the centred-third-invariant report band (§6.4); and the H0
  definition: H0 = the class physics absent, both corner arms
  responding as U does, all arms dressed by the same measured B̂ (the
  s²(B̂) contribution then cancels in the D-sign difference in
  expectation, and G1 measures its variance contribution).
- G2: estimator well-posedness at counts level: the N0-dressing
  residual calibration and N0-CLEAN thresholds; post-selection under
  asymmetric readout confusion; the U-arm demodulation treatment
  (§6.3, pre-decided); the (1,2)- and (2,3)-difference-null margins
  (separately: the (2,3) one rides the full inversion); fit-health
  thresholds; invariants-vs-sorted bias demonstration; where the
  linear chain becomes nonlinear.
- G3: dose certificates: per-site retention table AND the two-site
  realized retention |E[e^{−i(φ_l−φ_m)}]| over the M bindings for a
  lit/lit, a lit/dark, and a dark/dark pair (the two-site face has
  never been flown; a correlated-phase bug is invisible per-site and
  reads as a rate error); pass criteria with provenance (concentrator:
  every step within 0.02, artifact < 10%); M frozen from the
  realization error measured on s² itself, with the per-binding billing
  overhead question answered by measurement (query a past multi-binding
  job); the exactly-one-RZ post-transpile assertion.
- G4: the scattered-background budget, BOTH arms: a site-scattered
  profile (T1 as Γ_l/4-equivalent dephasing, gate-error scatter)
  splits the degenerate U arm from zero AND moves s²(C) at first order
  (the v1-gate per-arm figures, superseded by G4 v2's
  difference-quantity bands in §8a, kept only as provenance: U
  spurious ~12% of s(C), s(C) 1σ ±7%, both at the 0.25γ̄ scatter
  row). G4 additionally quantifies the same scatter on s²(C′) and on the
  W exceedance (v2 measures the D-sign difference; the C′-level
  quantities remain for the committed gate). G4 sets the T1 CLEAN conditions (staircase shape:
  named, individually stated, in-situ, sampling the arms' own time
  window) and the scattered-background band widening.
- G5: J-disorder and J-calibration: detuning-aware treatment under σ_J
  (priced; one frozen draw = a systematic widening, the f129 b_qs
  shape); the ∂s²/∂(J·dt) band term (§4); the in-analysis J estimator
  from the fitted dyad frequencies (§4) with its frozen pass band as a
  precondition on D-mag.
- G7: the runner gate (the concentrator 7b lesson): runner built,
  per-binding raw counts persisted before any reduction, `--analyze`
  proven end-to-end on a hardware-shaped synthetic artifact, the
  shared-skeleton and one-RZ-per-site-per-step and DD-off assertions in
  code, then an empty review of runner + gate records as its own stage.

(G-numbering note: G6 existed in v2 as the ω = −4 second demodulation
and is deleted, refuted in round 2; the number is retired, not reused.)

## 8a. Gate runs (2026-08-16, `simulations/corner_beat_gate.py` v2.1 after two promotion review rounds; outputs in `simulations/results/corner_beat/`)

The gate's own review round found and the v2 re-run fixed: a
winner's-curse selection (now two-stage on the VERDICT metric
P(d > θ_D), held-out seed, both f_leak ends), a gate-error bracket with
no bias channel (now a stochastic XX/YY magnon-hop channel inside the
evolution plus the Z-like half at the gate-count profile), a G3 that
never measured what §8 freezes M from (now the s²-realization spread
itself), a G4 on the wrong quantity (now the difference under a common
background), a G5 whose preparation basis followed the disorder (now
nominal-J basis, the flight's situation), missing readout confusion
(now symmetric 1.5% per qubit in the counts layer), and a
√w-weighting bug (~3% on the centre).

- **G1 v2.1 (two-stage, GAUGE-PINNED estimator; the pre-pin v2
  numbers rode a J·dt-unstable mode-sign gauge and are superseded):**
  stage-1, worst f_leak end: fractional beats standard in every cell;
  no standard-gate config reaches worst-end power 2 at 11 min.
  Stage-2 (200 reps, held-out seed), the frozen 21 × 3 grid, 8192
  shots, fractional: **d = +0.0049 ± 0.0017, power 2.9, θ_D = 0.0035,
  P(d > θ_D) = 0.815, H0 false rate 0.000 in-sample AND held-out,
  ~11 QPU min** for the 8192-shot variant. The C′ lines land where
  the promotion round predicted: s²(C′) = +0.0002 ± 0.0002, BELOW the
  hop-noise U floor (s²(C′) − s²(U) = −0.0015), while the W exceedance
  is robust. **THE FROZEN CONFIG (16384 shots, Tom's second decision),
  verdict-grade at 200 reps, both f_leak ends, held-out H0:
  worst end (f_leak = 8/15) d = +0.0048 ± 0.0017, power 2.9,
  θ_D = 0.0025, P(d > θ_D) = 0.920, H0 false 0.000; other end
  (f_leak = 0.9) power 4.2, P(detect) = 1.000, H0 false 0.020;
  W exceedance power 4.7 / 3.8.** p2 sensitivity at the frozen grid:
  power 2.7 at the recorded 0.5%/2q, 1.1 at 0.75%, dead at 1% (the
  Class-1 guard input; measured at 8192, conservative for 16384).
- **G3 v2:** retention pass probability over 40 table draws: M = 256
  passes 20%, M = 512 75%, M = 1024 100% (criterion 0.02, transplanted
  from the concentrator at a different dose and noted as such); the
  quantity §8 freezes M from, the s²-realization spread across frozen
  channels (re-run under the pinned estimator): **31% of centre at
  M = 256, 15.3% at M = 1024, 8.1% at M = 4096**. M = 1024 holds, and
  the 15% enters the D-mag band as a systematic (M = 4096 would halve
  it at 4 shots/binding; a G3 freeze-time trade).
- **G4 v2 (the verdict quantity, common scattered background):** the
  95% band on Δ(s²C − s²U) is 14% of the base difference at scatter
  scale 0.125γ̄, 28% at 0.25γ̄, 58% at 0.5γ̄, 111% at 1.0γ̄. **Proposed
  T1-CLEAN bound: scattered Γ_l/4-equivalent profile scale ≤ 0.25γ̄**
  (band under half the base).
- **G5 v2.1 (estimator level, nominal basis, pinned gauge):**
  σ_J = 1% costs ≤ 10% of s²(C) (95%), σ_J = 2% ≤ 25%; and the
  J-calibration response is **EVEN under the pinned estimator:
  ∂s²/s² = +4.5% for BOTH ±5% J** (any J mis-calibration INFLATES s²,
  a one-sided bias toward false D-sign, which is exactly why the G5
  J pass band is a precondition and enters θ_D's H0 model). The v6
  draft's compressed-level −27% was an artifact of letting the basis
  follow the disorder; the pre-pin ∓4% antisymmetric form was a
  gauge artifact.

Still outstanding for the committed gate: the N0 arm and var(B̂)
propagation (G2), the dressed centre functions f_C, f_{C′} as code
(f_{C′} carries the measured estimator dressing, §6 note below),
family rates at ≥ 500 H0 reps, the per-M retention criterion
replacing the transplanted 0.02, the lit/dark and dark/dark retention
pairs, the exactly-one-RZ post-transpile assertion, asymmetric
readout confusion (the v2 model is symmetric), and the within-sector
error split between the Z-like and hop faces (the gate pins 50/50;
s²(U)'s floor moves 0.0014 → 0.0003 across hop fractions 0.5 → 0, so
the C′-below-U statement is model-set in NUMBER, though C′ stays
below U across the whole bracket).

## 9. Guards

**Class 1 (pre-submit, prevents the spend):** line rule (every qubit
T2echo ≥ 150 µs, max/min T2echo ≤ 2, readout ≤ 2%); the realized
physical-qubit profiles of C and C′ in their asserted classes
(post-transpile assertion); FRACTIONAL-RZZ exposed and taken
(use_fractional_gates honored, rzz in the target's operation names, no
CZ decomposition in the transpiled circuits, CZ = NO-FLIGHT, the
record-parity assertion; load-bearing per §4/§8a); DD/resilience OFF
asserted; pending-queue
depth check; billing projection under the G1-frozen cap; the
pre-registration, gate, and runner committed at a real hash BEFORE the
Batch opens.

**Class 2 (in-job, protects the verdict, can only VOID):** the §7
bank: N0-CLEAN, T1 CLEAN on interleaved in-situ T1/T2* PUBs, dose
certificates on the flown phase tables, the (1,2) and (2,3)
difference-nulls,
fit health, the s²(C) floor, the G5 J pass band, band-validity window (if the device leaves the
window the design is invalid: abort, redesign, re-register).

**Day-of re-gate (pre-committed rule, anti-circular):** consumes
PRE-SUBMISSION calibration-class inputs ONLY (backend properties
snapshot; at most a separate pre-job calibration PUB, never the science
job's own CAL or T1/T2* PUBs, which are Class-2 guard inputs and would
be circular; the staircase's own fence); replaces band WIDTHS per the
G1-frozen scaling functions, never centers, never thresholds; its
output is committed as the Day-of addendum BEFORE the Batch opens. The
committed document governs, even against the runner's own printout
(the concentrator stale-band lesson).

**The order:** gate green → empty review of runner + gate records
(G7) → document + runner + constants committed at a hash → day-of
calibration pull → re-gate → Day-of addendum committed → Class 1
guards → submission → in-job Class 2 → analyze from persisted
counts → RECORD (with its own empty round, the house rule).

## 10. Costs (at the candidate point; billing anchors by measurement)

- Billing anchors: 0.283 ms/shot (concentrator, shallow), 0.309-0.327
  ms/shot (the two staircase flights, delay-bearing like this design):
  the band uses the delay-bearing anchors.
- Per-PUB budget table (frozen, §4/§8a): 4 shot-bearing arms ×
  21 depths × 3 preps = 252 science PUBs at 16384 total shots each
  (M = 1024 bindings × 16 shots) ≈ 4.13M shots, plus 2 CAL PUBs and
  ~8 in-situ T1/T2* PUBs ≈ 0.15M shots → ~4.3M shots ≈ **22-23 QPU
  min** at the delay-bearing anchors. Band quoted 21-25 QPU min,
  accepted by Tom 2026-08-16 (second decision, post-pin numbers); M
  is a G3 freeze output and the band re-freezes with it.
- Circuit depth at the candidate point: the RZ injection layer sits
  between Strang steps and does not commute with XX+YY, so half-layers
  do NOT merge across steps: 8 two-qubit XXZ blocks per step × 60 steps
  ≈ 480 blocks; at the pinned FRACTIONAL gates (2 per block) ≈ **960
  two-qubit gates**, deep-end error exponent ≈ 4.8 (standard gates:
  ~1440 and ≈ 7.2, the retired configuration) against a signal
  exponent γ̄T·(10/3) = 3.0 at the frozen point. The bridge from exponent to survival runs through
  the out-of-sector fraction f_leak of the error weight (post-selection
  removes only what leaves the sector): §5's 7-of-15 in-sector count on
  a differing pair gives f_leak = 8/15 → survival e^{−4.8·8/15} ≈ 7.7%,
  while X-type errors on non-differing pairs also leave the sector,
  pushing f_leak toward ~0.9 → ~1.3%; the honest range at the 60-step
  point is **~1.3-7.7% post-selected survival** (fractional; the
  retired standard-gate range was 0.2-2%). The v4-v5 reading that
  this kills the flight was WRONG in mechanism (§4, §8a): the deepest
  points do not need to survive at full SNR, because the estimator's
  kept-counts weights devalue them smoothly; G1's counts-level MC
  carries the whole depth profile and still resolves the split at the
  candidate budget WITH fractional-RZZ (2 gates per block ≈ 960 2q at
  the deep end). Standard gates (~1440 2q) remain unflyable in power
  terms at any tested budget. If fractional gates are unavailable on
  the flown backend (Class-1 guard, §9), the candidate is PARKED
  (recorded, not deleted); shallower-grid redesigns reopen it.

## 11. Honesty notes

- w = 1/3 and m̄ = 5/3 are gated-numeric at the eigensolver floor, not
  derived (M0 diagonality and the 2160/2304 denominators are e3dbab0's
  open residue). The rates inherit that status. A derivation landing
  before the flight upgrades the claim, not the protocol.
- The room is a uniform-chain resonance; bond disorder detunes it.
  The device requirement, remeasured at the estimator level under the
  pinned gauge (G5 v2.1): σ_J = 1% costs ≤ 10% of s²(C), 2% ≤ 25%,
  and the J response is even (any error inflates); the requirement
  stays σ_J ≲ 1-2%, owned by G5, cost carried in the D-mag band.
- No repo document states the Q-validity fence for this room's
  compression; this doc carries both corridor factors itself (§4),
  measured across three rounds with a 2-4% method spread, gate to
  arbitrate.
- The engineered γ̄ = J/10 = 0.1·J is a resolvability choice, exactly
  2× the canonical hardware-anchored γ₀ = 0.05 regime.
- v1 cited a ζ factor-2 from the price-pair doc; the pinned house form
  is ζ_c = ζ_shift/4 in the staircase doc. No ζ observable flies here;
  kept as a convention-trap reminder.

## 12. Scope fence (named follow-ups, out of scope here)

The mirrored maximizing transversal {1,2,5} (the l → 5−l image of the
pinned C) as a repeat arm; a second
non-maximizing corner as a deliberate positive control; other-N mixed
rooms (at Δ = 1 there are no mixed spaces at ω ≠ 0 at any N in 3..8
except 6, so the follow-up is ω = 0, Δ ≠ 1, or a topology variant, all
unexplored); a second backend / second day repeat; the derivation
upgrading w = 1/3 from gated to derived; the standard-gate variant if a target without fractional RZZ ever
must fly (retired by G1: worst-end power ≤ 1.7 at any tested budget
under the pinned estimator).

## 13. Amendment protocol

After freeze (§9's commit), changes ONLY as numbered pre-data
amendments appended here, each recording what changed, why, its
committed hash, and that no science data existed when it landed (the
staircase Amendments 1-2 shape). Post-data, nothing changes; the
RECORD applies the committed rules and names deviations as instrument
deviations.

## Day-of addendum (empty until flight day)

Filled by the §9 re-gate BEFORE the Batch opens: day-of governing band
widths, calibration snapshot, line selection with scores, Class-1 guard
outcomes.

## 14. Revision notes

- v1 (2026-08-16): first draft after the two-agent store sweep and the
  depth scout.
- v2 (2026-08-16): round 1 (physics recompute, spec audit, measurement
  statistics), ~35 findings folded; the two transversal classes, the
  Trotter EP, sort-free invariants, the prep 3× correction, the U-arm
  rebuilt.
- v3 (2026-08-16): round 2 (physics-math, spec-consistency), ~30
  findings folded; the ω = −4 null refuted and replaced, the secular
  factor re-signed, the split statistic pinned, D split into sign/mag,
  the N0 dressing moved to the full 3×3 generator, genre organs added.
- v4 (2026-08-16): round 3 (physics-math, spec-consistency), ~35
  findings folded, the load-bearing ones re-verified by the writer at
  the sources: **s² replaces s as the verdict statistic** (signed,
  real, resampleable; negative = conjugate pair, a reading not a
  failure); the envelope-immunity claim corrected (absolute s² biased
  up to −20% by the non-exponential post-selection envelope; ratios
  and r̄-differences immune; D-mag banded by G1 only); the (1,2) null
  restated as the C-minus-U generator difference (B̂ cancels) and
  split into a VOID-capable wiring check vs a report-only physics
  reading; the C′ arm given its own corridor and EP (Q ≈ 21 at
  J·dt = 0.15) and the W centre corridor-corrected (0.982·√3 ratio →
  2.9 on s²); the M-ladder re-anchored at the source verbatim (28.3%
  IS the M = 256 value; v3's "was M = 1" withdrawn); the candidate
  point declared priced-as-failing with the live alternatives named;
  arm-typed verdict rules replacing the mis-imported ternary grading,
  with equivalence tests for M and the null, floors under the W ratio,
  a partitioned ¬D-sign (FALSIFIED / Anti-D / underpowered), and named
  D∧¬M outcomes; N0-CLEAN conditions; the day-of re-gate made
  anti-circular (pre-submission inputs only) and its addendum moved
  before the Batch opens; DD/resilience-off, shared-skeleton, and
  deterministic arm order pinned; per-device twin band sets; billing
  re-anchored on delay-bearing flights.
- v4.1 (2026-08-16): round 4 (convergence check; every prior-round
  number re-verified to the digit, the seams caught): the secular
  factor stated per arm with C′'s OPPOSITE sign; the 2×2 Trotter model
  scoped to C (C′'s factor is ABOVE 1 at both live points; the table
  is the authority); the verdict evaluation made total with the order
  demoted to reporting precedence (FALSIFIED was unreachable);
  the ¬D-sign partition made disjoint and executable; D-mag's centre
  made a frozen function f_C(B̂); the arm-order rationale rewritten
  for s² (shift invariance kills the drift argument; W pair adjacent);
  s²'s sign statement made one-directional; the shape statistic moved
  to the centred third invariant; the s-vs-s² units pinned at the §4
  table and the J-term; the survival bridge (f_leak ∈ [8/15, ~0.9],
  0.2-2%) put on the page; the resilience pin softened to the house
  wording.
- v5 (2026-08-16): round 5 (final convergence check; physics verified a
  fifth time): the W arm's centres made frozen functions f_C(B̂),
  f_{C′}(B̂) (the last one-corner-arm number); s² pinned as
  Re[(e₁²−3e₂)/3] of the complex fitted R (decidable under the U-arm
  complex fit; Im parts = the detuning readout); FALSIFIED given its
  own physics-sized equivalence margin θ_F (containment in the
  detection band was scale-invariantly impossible); C′'s EPs at both
  live points on the page (Q ≈ 38 candidate, Q ≈ 21 vs flown 10 at
  fallback); the pair-odd C_p convention pinned and the swapped
  SUM/DIFFERENCE aside fixed (maximizing = difference, norm 1/3 = w;
  non-maximizing = sum, norm 1/√27); the non-maximizing triple written
  as three rates with (2,3) named as its decoupled dyad, and the (2,3)
  C′-minus-U difference null added as C′'s wiring check; the s²-floor
  VOID marked physics-or-instrument ambiguous with the RECORD naming
  both; H0 defined; the exact 1/6 route scoped to the room dyads; G4
  extended to s²(C′) and the ratio; survival 0.8% replaced by the
  f_leak bracket; s-vs-s² factor-2 notes at G4.
- v5.1 (2026-08-16): round 6 (propagation check; physics verified a
  sixth time, all repairs confirmed correct where present): the (2,3)
  C′ null propagated into §7 (its own equivalence test), G2 (separate
  margin), and §9 Class 2, trigger count 8; f_{C′}(B̂) and θ_F and the
  branch-reachability checks added to G1's freeze list; the W arm
  narrowed to the RATIO alone (the registered class arm) with both
  corner magnitudes moved to D-mag, resolving the collision with the
  registration rule; the Trotter detuning renamed δ_T with units (Δ is
  the anisotropy alone); site-pair {a,b} vs mode-dyad (i,j) notation
  pinned; the §12 mirror image corrected to {1,2,5}; the middle-rate
  ordinal, the number-operator-compression clause, and the M
  single-owner (G3) fixed.
- v5.2 (2026-08-16): round 7 (verdict-logic diff check, ZERO blockers):
  the floor-trip adjudication frozen (dose+N0 clean → ¬W class
  anomalous, else VOID; §3's partition restored); the G5 J pass band
  registered as VOID trigger nine; the anchor arm renamed M → A (the
  binding count keeps M, G3's object); the ¬D-sign partition made
  exhaustive (the above-θ_F case named); the named-verdict table
  declared governing over the reporting order; the null bullet made
  plural with separate margins; G1's equivalence margins scoped to the
  A arm and the projected SE added as its deliverable; the mirror-pair
  index clause; the ideal tag on the §1 ratio; the registration
  sentence extended to record the other conjuncts. Design stage
  CONVERGED: the next work is the committed gate (G1-G5, G7),
  Tom-gated.

- v6 (2026-08-16, gate-build stage): the gate's first G1/G3/G4/G5 runs
  folded (§8a): working point moved to Q = 10, J·dt = 0.15, 13 × 5
  grid, fractional-RZZ REQUIRED (now a Class-1 guard) and the
  priced-as-failing status REVERSED (the deep points are weight-
  devalued, not fatal); the W arm restated from the ratio to the level
  discrimination after failing its own reachability check (s²(C′)
  detection power 0.6 at budget); the s²(C′) floor retired; M = 1024
  confirmed from below (M = 256 FAILS the two-site retention
  criterion); budget re-tabled at 156 science PUBs ≈ 7 QPU min.
- v7 (2026-08-16, gate v2): the gate's own empty review round (four
  blockers: non-committed estimator, bias-free error bracket,
  winner's-curse freeze, G3 measuring the wrong M quantity) folded and
  the gate re-run; §6 AMENDED to pin the eigenchannel damped-cosine
  fit as the committed estimator (the complex one-step fit withdrawn
  for measured cause) and the negative-s² machinery withdrawn with it;
  the freeze moved to the 21 × 3 grid on the verdict metric
  (P(detect) 0.985 vs 0.84, ~11 min vs 6.8); §8a rewritten from the
  v2 run; p2 ≤ 0.6% and both-twins fractional checks added to
  Class 1; the T1-CLEAN scale ≤ 0.25γ̄ proposed from G4 v2; σ_J
  re-measured at the estimator level (8%/21% at 1%/2%); the v6
  consistency round's fix list (labels, γ̄ = 2×γ₀, grid remnant,
  exponent-7 remnants, W(ii)/D-mag duplication, twin rule) applied.
- v7.1 (2026-08-16): Tom's budget decision recorded (21 × 3 at ~11
  QPU min); the promotion round ran and BLOCKED promotion: the §7 C′
  numbers were from the withdrawn estimator (sign reverses under the
  committed one: s²(C′) reads below the hop-noise U floor, while the
  W exceedance s²(C) − s²(C′) survives at power ~4.7), and the
  committed estimator's dressing is large (s²-level ~1.5×/0.3× of
  ideal for C/C′) and belongs to f_C/f_{C′}; the writer's diagnostic
  resolved the dressing's mechanism (a J·dt-dependent mode-sign gauge,
  now PINNED to the continuous-H overlap in the gate, plus a
  prep-dependent secular population leak, real and deterministic);
  gate v2.1 adds the gauge pin, stage-2 C′/W lines, a held-out H0
  split, and the p2 row at the frozen grid; §7's W bullet is restated
  from the v2.1 run; em dashes swept, M→A completed, textual fixes
  from the promotion round applied.

- v7.2 (2026-08-16): Tom's SECOND budget decision on the honest
  post-pin numbers (option b: 21 × 3 at 16384 shots, ~22 QPU min);
  promotion round 2 blocked again on the same shape one level deeper
  (the gauge pin carried through G1 only), so G3/G5/dressing were
  re-run under the pinned gate and re-quoted (M-realization 15.3% at
  M = 1024; σ_J 10%/25%; dressing 1.21×/0.54×; the J response found
  EVEN, +4.5% both directions, a one-sided inflation bias entering
  the G5 pass band's rationale); the frozen-configuration numbers
  moved to the committed `--frozen` mode (outputs archived under
  `simulations/results/corner_beat/`); the pre-pin leftovers swept
  (status, §5 16/binding, §7 θ_D and W lines, §10 exponent 3.0,
  standard-gate ≤ 1.7, G4 provenance row relabeled); the hop-split
  bracket added to the outstanding list; §8/§8a retitled to the
  committed state.

## Open questions for Tom (spend decisions only)

1. ~~The budget fork~~ **DECIDED TWICE, finally (Tom, 2026-08-16
   ~15:53, on the honest post-pin numbers): option (b), the 21 × 3
   grid at 16384 shots, ~22 QPU min.** The frozen-config
   verdict-grade confirmation (200 reps, both f_leak ends, held-out
   H0) is IN §8a: worst-end P(detect) 0.920 at H0 false 0.000.
2. The C′ arm stays per your earlier go (its cost is inside both
   numbers above); the parking rule now reads: park if fractional-RZZ
   is unavailable on the flown backend, if day-of p2 > 0.6%/2q on the
   used edges, or if the committed gate fails to reproduce §8a.
   Confirm.
3. Committing the gate + this document (the promotion rename and the
   first commit of the arc) happens after the next empty round on the
   re-run gate; no QPU spend before a separate go from you.
