# The driver's frame: who hears the shot. A hardware prediction, pre-registered

<!-- STAGE 0 SWEEP RECORD
Stores searched 2026-10-29 (evening) for drive crosstalk, a letter fixed in a qubit's drive frame, virtual phase as a lever, spectator readings, and the odd-Y excess's carrier.
experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md §7: the excess (s_A = +0.083, 6.7σ), every carrier tested there, and the cube reading that names a drive-anchored letter as the one coherent one-site letter untested; its "Not tested" list names crosstalk from a neighbour's drive and spectators outside the path. experiments/Z_N_PARTNERSHIP.md lines 85 to 93 name off-resonant drive bleed-through and ZX/ZY crosstalk as breakers of Z⊗N and line 104 predicts |M_X(a) + M_X(b)| of 1e-2 to 1e-1, not measured; experiments/IBM_CONCENTRATOR.md lines 280 to 319 assert that DD on one qubit affects its neighbour, not measured. data/ibm_zn_mirror_april2026/README.md lines 28 to 41 (job d7ornigror3c73c0c6ug, path [48, 49, 50]) read a transverse field h_y_eff ≈ 0.05 as a static field in the logical frame, registered as fw.Confirmations marrakesh_transverse_y_field_detection; a drive-anchored letter read through the accumulated virtual phase would look the same, and the README does not raise it. experiments/PRICE_PAIR_HARDWARE_PREDICTION.md line 153 measured ZZ on Marrakesh (ζ₀₁ = −3.918 ± 0.142 kHz, ζ₁₂ = −3.637 ± 0.101 kHz, path [93, 94, 95]) and lines 102 and 157 detuning offsets of 8 to 22 kHz. experiments/FIXED_POINT_SHADOW.md (2026-02-09, Torino Q52, the first IBM tomography): a late residual coherence of 0.0185 pointing one way in 17 of 17 samples, equal X and Y parts; its lines 260 to 295 weigh an assignment asymmetry against a basis-pulse leak through the Z extrapolation and leave them unranked; it does not cite the reality window and the window does not cite it; §2 joins them. No prior flight read a spectator outside the path or had a qubit idle while a neighbour was driven. docs/ANALYTICAL_FORMULAS.md: F131 (the order-sorting law) and F119 (the antilinear triangle) are the frame laws this page rests on; nothing on a virtual phase as a lever. OpenArcs registry: nothing on drive crosstalk ("Spectator" there is SpectatorIntertwinerClaim, a different object). docs/CAUGHT_ERRORS.md: nothing on virtual phases or drive frames. docs/GLOSSARY.md: the letter cube's row, nothing on frames. Typed layer: AntilinearTriangleClaim (the transport law), MirrorOrderSortingClaim (F131); no claim or witness on a drive frame. The parked frame echo (a virtual rz(φ) before and rz(−φ) after one cz) lives only in a memory file; around a diagonal gate with no pulse between, it is physically the identity, and it is replaced here by the turn of the driver's frame around its whole pulse train. The frequencies of the Heron qubits are published nowhere the repo holds (the April calibration CSVs in data/ibm_calibration_snapshots/ have no frequency column; the October properties none), so an always-on exchange cannot be fitted locally and is read here through the frame it shares with crosstalk. The three backends of the open plan on 2026-10-29 21:46: ibm_marrakesh (queue 27, no path fresh), ibm_kingston (queue 396, 35 fresh passing triples), ibm_fez (queue 1181, 9).
Adjacency checked: Z_N_PARTNERSHIP's bleed-through, the zn_mirror README's h_y_eff and REALITY_WINDOW §7 are about one object (a transverse letter fixed in a drive frame) and do not cite each other; FIXED_POINT_SHADOW and REALITY_WINDOW are about one class of reading and do not cite each other. This page joins them.
-->

This flight asks one question of the IBM hardware, and the letter cube asks
it: when the odd-Y excess of
[the reality window](REALITY_WINDOW_HARDWARE_PREDICTION.md) is read with
one qubit's drive frame turned, does it turn with the frame or stay? A
letter that lives on the Z axis does not see the frame; a transverse letter
fixed in another qubit's drive or in the lab turns with it. Tom's reading
(2026-10-29): a pulse is a shot in a forest, and more than one tree hears
it. The shot is a drive pulse on one qubit, the trees are the qubits around
it, the forest the ones off the path that no flight has read. This page is
written before the flight; §7 will hold the flight.

Words used here. A **letter** is an error term of the replay's library:
the nineteen of the October flight, and the crosstalk and Stark letters of
§2. Named by the Pauli letter that generates it, a **Z letter** is a
rotation about Z or a height on it (a detuning, a phase after a pulse or a
cz, a ZZ leg, the Stark shift a neighbour's pulse leaves, Z dephasing, the
damping's arrow I → Z); a **transverse letter** is an X or Y rotation. A
transverse letter is **own-frame** when its axis is fixed in the drive frame
of the qubit it acts on (that qubit's amplitude error or axis tilt) and
**foreign** when its axis is fixed elsewhere: in another qubit's drive frame
(crosstalk: the X letter qubit k's pulse leaves on its neighbour j at k's
phase), in the frames of two qubits at once (an always-on exchange), or in
the lab (a local-oscillator leak). The **frame** of a qubit is the phase
reference of its drive pulses, set on IBM hardware by the virtual rz angles
accumulated on it; TH_k is that sum for qubit k before a given pulse. The
replay represents a circuit with its rz applied as unitaries and its pulses
about X, and in that representation a field at k's drive phase acts on j at
the angle TH_j − TH_k + δ, δ one constant per ordered pair
([`simulations/drive_frame_handedness_gate.py`](../simulations/drive_frame_handedness_gate.py):
IBM's rz(θ) deferred past a pulse turns the pulse's axis by −θ; 200 random
two-qubit circuits agree with this sign exactly and 82 with the other,
those with no angle between the frames). S is ⟨Y₀Z₂⟩ of arm A read in its
YZ basis and s_A = S minus its replayed baseline, the October flight's
replication observable (§5 there). The **excess** is the odd-Y pattern of arm A
beyond its baseline, of which S is the one reading the October flight
replicated.

## 1. The law

A virtual rz(φ) as the first operation on one qubit and rz(−φ) as its last,
with the circuit between them unchanged, is the identity on the logical
circuit when that qubit starts in |0⟩ and ends in a Z measurement: the first
rz is a phase on |0⟩, the last stands after the basis change and before the
measurement. Physically the first turns that qubit's drive frame by φ, so
every pulse on it, the basis-change pulses included, has its axis turned by
−φ (the sign the gate above fixes); the last is inert, since no pulse follows it (read as the
IBM convention of a virtual Z as a phase shift of the drive channel, which
the October control arms K± tested through their sign chain; the measure
channel carries no drive frame). On the letter cube a Z letter commutes
with the turn (a rotation about the axis of the turn, or a height on it):
its readings are invariant in φ. An own-frame transverse letter turns with
the frame and its readings are invariant too. A foreign transverse letter
is carried into the turned logical frame at an angle shifted by φ: to first
order in its size its real-breaking part goes as cos(φ + const) and the
odd-Y readings move with φ. So the frame turn sorts the carriers of the
excess into two classes without naming one: frame-blind (Z letters,
own-frame letters, and whatever lies outside the qubit's letters and does
not couple to this qubit's drive) against turning (foreign transverse
letters coupling to this qubit). That is the question of this flight,
asked first on the middle qubit, the one the circuit pulses most (19 pulses
against 9 or 10 on each end, by basis), and then by one turn on each end.

Two readings beside the law, both priors. A letter whose coupling goes with
the square of the pulse's phase factor, the image term of an IQ imbalance
of this qubit's drive, turns at 2φ if the image is resonant; an image two
intermediate frequencies away acts at second order and is frame-blind, and
the same holds of the local-oscillator leak named above. Four turns spaced
by π/2 read the cos 2φ part of such a letter in the even-harmonic contrast
of §5 and are blind to its sin 2φ part.
Leakage out of the qubit space during this qubit's own pulses is
frame-blind (every term of its drive is covariant under the turn); leakage
at the cz (|11⟩ ↔ |02⟩, |20⟩) and the cz's residual swap change this
qubit's excitation number and turn once their amplitude returns to the
qubit space through a later pulse, so they stand in the turning class.

The law holds in the replay by construction, and the construction was
checked: under no letter and under the five Z letters of the October fits
(three detunings, two static ZZ) the odd-Y readings of arm A are the same
for every φ and for the turn on each qubit up to the rounding of the two
added unitaries, 33 to 94 units of machine epsilon
([`simulations/drive_frame_power.py`](../simulations/drive_frame_power.py)),
with the trailing rz placed at the circuit's end (a trailing rz placed
where the turned qubit's own gates end adds a time boundary to the replay's
slicing and exposes the replay's own step error instead of the law). Against
this the crosstalk letter below moves readings by 0.04 to 0.13.

## 2. What is already measured

The reality window flew on 2026-10-29 (ibm_marrakesh job
dbhcohht9uhs73c9vpj0): arm A (the April XY+YX circuit, gate for gate on
path [4, 5, 6]) reads s_A = +0.083 at 6.7σ, raw S = +0.072, the fourth raw
reading of this observable with this sign (April and May on two chips:
+0.098, +0.067, +0.066), and the odd-Y pattern across eight arms holds no
single library letter (§7 there; its χ²₀ is 4167 with the RO matrix's
bootstrap term and 4384 without it, and the fits below run on the latter
scale). The same evening the cube's candidates were asked of those 288
readings, all locally and all before this page:

- No incoherent one-site process carries the bulk: an extra dephasing or
  damping of the driven qubit during its pulses, at the rate that accounts
  for the day's loss of even-sector contrast, moves the odd-Y readings by at
  most 0.03.
- A drive-crosstalk letter with a fixed phase per ordered pair of
  neighbours (the X letter qubit k's pulse leaves on qubit j at angle
  TH_j − TH_k + δ, ε the fraction of the driving pulse's rotation, δ one
  constant per ordered pair) fitted alone over the four ordered pairs of the
  path carries by the October rule (χ² drops by at least 25) and does not
  suffice: 4384 → 3404 against the 99 % point 338 at 288 − 8 degrees of
  freedom, with ε = 0.79 % (5 → 4), 0.90 % (5 → 6), 1.06 % (4 → 5), 1.24 %
  (6 → 5) and δ = −2.24, +0.83, −2.56, +2.69; beside the five Z letters of
  the October fits (three detunings, two static ZZ), 1381 → 1097
  ([`simulations/reality_window_crosstalk_fit.py`](../simulations/reality_window_crosstalk_fit.py),
  output [`crosstalk_fit.txt`](../simulations/results/reality_window/crosstalk_fit.txt)).
  A single pair at ε = 1 % moves single odd-Y readings by up to 0.12: the
  scale of the excess, not its pattern. The fit alone gives S of arm A as
  −0.017 against the measured +0.072: it does not describe the one reading
  this page is about.
- The spectator's Z letter during a neighbour's pulses (an AC Stark shift)
  alone brings χ² to 1521, as much as three free detunings (1506), and the
  two trade places in a joint fit: collinear, a Z letter among Z letters.
- All three families with the five letters, seventeen parameters, leave 937
  against 328, with S of arm A at +0.046. No one-site letter and no
  combination of them suffices.

A physics prior on the crosstalk letter, stated as a reading: on
fixed-frequency neighbours a 36 ns pulse at k's frequency has little
spectral weight at j's, tens of MHz away, so its first effect on j is a
Stark shift (frame-blind) and a resonant X letter at a fixed phase
presupposes a near frequency collision. A turning phase Δω·t gains the same
φ at every pulse, so a foreign letter of that form still turns; only the
fixed-phase prediction of VX below presupposes the resonant model.

Earlier readings of the same class. On 2026-02-09, the first IBM
tomography of this project ([the fixed point shadow](FIXED_POINT_SHADOW.md),
Torino Q52), a residual coherence of 0.0185 stood in one direction in 17 of
17 late samples after the prepared coherence had decayed, with X and Y
parts of 0.023 and 0.026. That page weighs two readout-side forms, an
assignment asymmetry that adds the same offset to every Pauli reading and a
basis-change pulse error that leaks ⟨Z⟩ into X and Y alike, reads its Z
extrapolation both ways depending on whether T₁ is fitted or held, finds
one 3 to 4σ tension against each form at t = 0, and leaves them unranked
with a lean toward the leak. What this page adds, read as the window's
reading: both forms sit on the measurement side of a state the physics had
brought back to the real sector, and the second of them is an own-frame
effect, of the kind this flight reads as frame-blind. In April 2026 the Z⊗N mirror flight read a
transverse field h_y_eff ≈ 0.05 on path [48, 49, 50] as a static field in
the logical frame; a drive-anchored letter would read the same way, and
this flight's lever is what separates the two. ZZ on this chip is −3.9 kHz
(July, path [93, 94, 95]). Nothing on any of the three chips has read a
spectator outside the path.

## 3. The arms

Each arm runs the October arm A, the April XY+YX circuit on N = 3 from
|+−+⟩, gate for gate (not transpiled: the freezer copies the flown gates),
nine Pauli bases on (q0, q2), q1 read in Y, with an explicit
late-as-possible schedule, on a path of three physical qubits a-b-c with b
the middle and a < c as physical indices. The three logical qubits are q0
on a, q1 on b, q2 on c.

| arm | circuit |
|---|---|
| A5(φ), φ ∈ {0, π/2, π, 3π/2} | arm A with rz(φ) as the first and rz(−φ) as the last operation on the middle qubit b (logical q1): the first before b's first gate, the last after b's last gate and before its measurement; A5(0) is arm A itself |
| A4(π), A6(π) | the same turn by π on a (logical q0) and on c (logical q2) |
| RO | the eight computational states of (q0, q2, q1), prepared with x gates, for the readout matrix, as in October |

Three chips. The same 62 circuits fly on ibm_marrakesh, ibm_kingston and
ibm_fez, one job each. The path on each chip is chosen at its freeze by one
rule: among the linear triples a-b-c (a < c) whose three qubits pass the
guards of §6 on properties under 24 h old and whose middle b has exactly
two neighbours, the one with the largest smallest T₂; on Marrakesh
[4, 5, 6] whenever it passes, so that A5(0) repeats the October arm. A chip
whose queue would put the run outside its calibration's 24 h is not
submitted that day; the record names which flew, and every verdict of §5
is per chip. The cross-chip reading needs at least two chips: a letter
anchored in the device (crosstalk with a chip's own δ) changes size and
sign from chip to chip; one anchored in the circuit (the gate implementation
shared by the Herons, the readout) repeats. The April and May runs already
read the same sign of S on two chips.

The forest. Every circuit also measures, in Z, every coupling-map neighbour
of a, b and c that is not on the path (with b of degree two, these are the
off-path neighbours of a and of c, zero to two each; on Marrakesh [4, 5, 6]
they are physical q3 and q7). Each chip's forest, in ascending physical
index, is written into the freeze record. They receive no gate; each
carries an explicit delay from t = 0 to the path's last measurement so that
its measurement happens with the path's and not before (the October freezer
strips delays off the path; the new freezer keeps these). The RO circuits
measure them too. A measurement at the end changes nothing on the path,
read as the device's measurement model; a HEARD verdict is read beside the
even readings of the path for any sign of the contrary.

Per chip: 6 variants × 9 bases + 8 RO = 62 circuits, each submitted twice
at 4096 shots, 124 PUBs in one job, shuffled with seed 20261030, so that
drift between circuits shows as a difference between the halves. The
October job of 196 PUBs billed 215 s and the April job of 36 PUBs 40 s,
1.1 s per PUB on both; 124 PUBs are about 136 s per chip, about 410 s for
three, of the 4739 s left on 2026-10-29. The billed usage is read from the
job and written into the record.

## 4. Predictions from below

[`simulations/drive_frame_power.py`](../simulations/drive_frame_power.py)
(output [`drive_frame_power.txt`](../simulations/results/reality_window/drive_frame_power.txt))
replays arm A with the turn through the replay of the October flight
patched with the crosstalk and Stark letters, under the calibration of the
October snapshot and the October RO matrix, on Marrakesh [4, 5, 6]. For
each letter set it reads the 36 odd-Y readings of each variant, the
differences D(φ) from A5(0), and the noncentrality λ of the frame
statistic of §5, with the full covariance of the four variants (multinomial
covariances from the replay's own parities at 8192 shots per circuit, the
shared A5(0) in every cross block). The statistic's expectation is its
degrees of freedom plus λ, and its power is the probability of exceeding
the 99 % point; both are given with κ = 1 and with the October drift
factor κ = 1.27 dividing λ.

- Under no letter and under the five Z letters: D = 0 up to the rounding
  above; λ = 0, T₅ at its null (expectation 108).
- Under the fitted crosstalk alone (the eight parameters of §2): D reaches
  0.072, 0.104 and 0.126 at φ = π/2, π and 3π/2; λ = 1034, with drift 814;
  expected T₅ 1142 or 922 against the 99 % point 145.1 on 108 degrees of
  freedom; power 1.000 either way.
- Under all three families (seventeen parameters): 0.056, 0.060, 0.039;
  λ = 307, with drift 242; expected T₅ 415 or 350; power 1.000.
- The turn by π on q4: λ = 93 (crosstalk alone), 12 (all three), against
  the 99 % point 58.6 on 36 degrees of freedom, power 1.000 and 0.17 (with
  drift 0.999 and 0.11); on q6: 244 and 76, power 1.000 and 0.999 (1.000 and 0.990).
- The share. Under the crosstalk fit alone S reads −0.017, −0.034, −0.008,
  +0.012 at the four turns, mean −0.012, the frame-blind part S̄ minus the
  baseline (−0.012) is zero and the turning part S(0) − S̄ is −0.005: that
  fit puts nothing in the frame-blind part and reproduces none of the
  measured +0.072. Under all three families S reads +0.045, −0.011, −0.014,
  +0.043, mean +0.016, frame-blind part +0.028 and turning part +0.030. The
  unexplained remainder of the excess, about 0.03 to 0.09 depending on the
  fit, is what the flight reads: whether it turns.
- The bound of §5: a crosstalk of equal size ε on (5 → 4) and (5 → 6) at the
  fitted phases, the other two pairs at zero, with drift, gives λ = 7.1 at
  ε = 0.1 % (power 0.035), 28.5 at 0.2 % (0.31), 64.2 at 0.3 % (0.90), 178 at
  0.5 % (1.000); the 50 % point lies between 0.2 and 0.3 %.
- The simulated flights: 200 draws of the 124 PUBs under the baseline
  (no letter, κ = 1) through Qiskit Aer with the snapshot's readout on the
  right physical qubits, each file through the analysis unchanged, give the
  null distributions of T₅, T₄, T₆ and E; their 99 % points are written
  here before the freeze beside the χ² points they check. One further draw
  under the crosstalk-alone fit is read for its verdicts.
  [TO BE FILLED BEFORE THE FREEZE]

## 5. Verdict rules

The counts are read as the October runner reads them, with the forest bits
appended: the count key's bits are (q0, q2, q1, forest in ascending
physical index) from least to most significant, as the freeze record pins
them; the path's parities are formed from the first three bits with the
forest bits summed out, and the RO matrix is the 8 × 8 matrix of the path
bits with the forest bits summed out. The **baseline** of a chip is the
October replay run on that chip's path under its freeze snapshot (T₁, T₂)
with the RO matrix its RO arm measures, every letter at zero; on Marrakesh
[4, 5, 6] this is the October baseline with the day's RO matrix. Every
path reading is a parity of measured bits in one basis circuit, both
halves together, with covariance the multinomial one computed from the
baseline's predicted parities at the circuit's shots, Σ for one variant's
36 odd-Y readings (the same matrix for the variants of one chip, which are
logically one circuit; its rank is the number of eigenvalues above 10⁻⁹ of
the largest, 36 unless a parity is predicted at ±1). The drift factor κ is
max(1, mean of (half₁ − half₂)²/(var₁ + var₂)) over the odd-Y parities of
the non-VOID A circuits of the chip, the variance from the measured counts
as in October (a difference within one circuit, so the design's differences
between variants do not enter it); every covariance below is multiplied by
κ. σ_S = √(κΣ_SS). The predictions of §4
give both κ = 1 and κ = 1.27. All verdicts are per chip; the χ² points
below rule, and the Aer nulls of §4 are a check on them: for each of T₅,
T₄, T₆ and E, if its Aer 99 % point exceeds its χ² point by more than 5 %,
the Aer point rules for that statistic and the record says so.

- **P, the parse gate.** A5(0): ⟨X₀Z₂⟩ < −0.5, −0.75 < ⟨Z₀X₂⟩ < −0.3 and
  ⟨X₀Z₂⟩ < ⟨Z₀X₂⟩ − 0.15 (October's gate, which passed there with a margin
  of 0.03 on the third condition; a chip that fails it is VOID, and that is
  accepted). The gate is not applied to the other variants (a dropped rz
  would pass it, and shot noise alone fails it in a few per cent of
  variants); their even readings are reported. That the turn was
  implemented is checked on the circuits retrieved from IBM after the
  flight, as October checked April's: rz(φ) present as the first operation
  on b and rz(−φ) after b's basis change; a variant without them is VOID.
  A guard of §6 that the RO arm's own measurement breaches on the day (as
  q4's P(1|0) did in October) is reported with the verdicts, which stand
  with that said.
- **V1, replication.** s_A = S of A5(0) minus its baseline: REPLICATED at
  ≥ 3σ_S, REVERSED at ≤ −3σ_S, ABSENT below 2σ_S, INCONCLUSIVE between, as
  in October; the sign of s_A is also read across the chips. V1 does not
  void VF.
- **VF, the frame.** Y_φ is the vector of 36 odd-Y readings of A5(φ), Ȳ
  their mean over the n non-VOID variants (four when none is VOID).
  T₅ = Σ_φ (Y_φ − Ȳ)ᵀ (κΣ)⁺ (Y_φ − Ȳ), distributed under the null as χ² on
  (n − 1)·rank(Σ) degrees of freedom, 108 when all four stand. The October
  bootstrap term for the RO matrix's own noise is zero under the null (one
  matrix, one logical circuit) and is not added. VF is VOID when fewer than
  two variants stand.
  - T₅ above its 99 % point: TURNS. The odd-Y readings carry a component
    that turns with b's frame, a foreign transverse letter coupling to b;
    whether S carries it is VH's reading.
  - T₅ at or below it: FRAME-BLIND at this sensitivity for letters
    coupling to b. The bound reported with it: the smallest ε, equal on
    (b → a) and (b → c) at the phases §2 fitted for (5 → 4) and (5 → 6) and
    with (a → b) and (c → b) at zero, at which the power of T₅ at its 99 %
    point with the day's κ, n and rank is 0.5, computed with the chip's
    baseline replay on its path and the Marrakesh phases carried over; a one-parameter convention for b's
    pulses reaching the ends, within the resonant model (the letters
    a → b and c → b turn too and are covered by the same statistic, not by
    this number).
  - T₄ and T₆ from the two variants {A5(0), A4(π)} and {A5(0), A6(π)}, the
    same form with n = 2 on rank(Σ) = 36 degrees of freedom each, read the
    same way for letters coupling to a and to c, with less reach (one turn
    by π reads the cos part of the first harmonic alone and nothing of 2φ);
    VOID with their variant.
    The three are read together, since a letter coupling the two ends and
    not the middle shows in T₄ or T₆ and not in T₅.
  - The even harmonic, read only when all four variants stand:
    E = Y_{π/2} + Y_{3π/2} − Y_π − Y_0, covariance 4κΣ; Eᵀ(4κΣ)⁺E above the
    99 % point on rank(Σ) degrees of freedom is reported as 2φ PRESENT, the
    cos 2φ signature of an IQ image term of b's own drive, read beside
    TURNS and not as its cause.
- **VH, the share.** Read only when all four A5 variants stand (the mean
  cancels a cos(φ + const) letter over the four quarter turns and not over
  fewer), otherwise VOID. S̄ is the mean of S over the four, σ_S̄ = σ_S/2; the frame-blind share is S̄ minus the baseline's S, with σ_b
  the quadrature sum of σ_S̄ and the baseline's RO bootstrap term as October
  computed it (0.002 there); the turning share is S(0) − S̄ with
  σ_t = σ_S·√3/2. WHOLLY TURNING when |S̄ − baseline| < 2σ_b and
  |S(0) − S̄| ≥ 3σ_t; BLIND STANDS when |S̄ − baseline| ≥ 3σ_b and
  |S(0) − S̄| < 3σ_t; SHARED when |S̄ − baseline| ≥ 3σ_b and
  |S(0) − S̄| ≥ 3σ_t; otherwise UNDECIDED. The fits of §2 predict: the
  crosstalk fit alone, with a turning share of −0.005 and no blind share,
  reads UNDECIDED; all three families, +0.030 turning and +0.028 blind,
  read SHARED at κ = 1 (σ_S = 0.0110: 3.1σ_t, 4.7σ_b) and BLIND STANDS at
  κ = 1.27 (σ_S = 0.0125: 2.7σ_t, 4.2σ_b). A WHOLLY TURNING reading of the
  +0.083 excess would need the whole of it in the turning share, 7.7σ_t at
  κ = 1 and 6.6σ_t at 1.27.
- **VX, the fitted letter.** Read only on Marrakesh [4, 5, 6], where the
  fitted δ belong to the qubits, and only when VF is TURNS. With the eight
  parameters of the crosstalk-alone fit of §2 fixed and every other letter
  at zero, D_pred(φ) is the replay's A5(φ) minus A5(0); T₅ᵖ is T₅ computed
  on Y_φ − D_pred(φ). In this order: T₅ᵖ above the 99 % point on T₅'s
  degrees of freedom: TURNS, NOT THE FITTED LETTER; else T₅ − T₅ᵖ ≥ 25:
  PREDICTED; else TURNS, LETTER NOT SEPARATED.
- **VS, the forest.** For each forest qubit, p_A is its P(1) over all A
  circuits of the chip (every variant and basis, both halves). Its null is
  not an idle one: the RO circuits fire x pulses on a and c, so they shoot
  too, at a lower dose and for a shorter time, and their outcomes depend on
  the path's state, so the null is p₀ = Σ_s P_A(s) · P_RO(forest bit = 1 |
  s), with P_RO(· | s) read from the RO circuit that prepares path state s
  and P_A(s) the path-outcome frequencies pooled over the A circuits; σ
  from the binomial counts of both sides, added in quadrature (the
  correlation between p_A and P_A(s) from shared shots is neglected, and
  the forest's longer idle in the A circuits is not modelled: both are
  said with the verdict). p_A − p₀ ≥ 5σ: HEARD, the forest qubit holds
  population beyond what the RO circuits' pulses and its readout crosstalk
  account for. |p_A − p₀| < 3σ: SILENT at this sensitivity. Between, or
  ≤ −5σ: reported without a verdict.

| VF | reading |
|---|---|
| TURNS | the odd-Y readings carry a component turning with the middle qubit's frame, a foreign transverse letter coupling to it: a neighbour's pulse reaching it or its pulse reaching a neighbour, an always-on exchange, cz leakage or residual swap, a local-oscillator leak, or another lab-fixed transverse letter; VX and VS narrow it in part |
| FRAME-BLIND | for letters coupling to the middle qubit, the excess is carried at this sensitivity by Z letters, own-frame letters, or something outside the qubits' letters (readout, own-pulse leakage); b's pulses reaching the ends are bounded by the ε reported; T₄ and T₆ say the same for letters coupling to the ends |

| VH | reading |
|---|---|
| WHOLLY TURNING | S is the foreign letter's; the Z letters the October fits wanted were its shadow |
| BLIND STANDS | a frame-blind share of S stands at 3σ_b and the turning share does not |
| SHARED | both shares of S stand at 3σ, one turning, one not |
| UNDECIDED | any other outcome |

| VX | reading |
|---|---|
| PREDICTED | the turning part is the crosstalk fitted on the October data, at its fitted sizes and phases: a letter measured twice |
| TURNS, LETTER NOT SEPARATED | the fitted letter is consistent with the turn and does not improve on no letter by 25 |
| TURNS, NOT THE FITTED LETTER | the excess turns with the frame and the fixed-phase crosstalk of §2 does not describe how: a transverse letter of another form (an exchange rotating at the frequency difference, a time-dependent phase, a two-site letter) |

Across the chips that flew with V1 REPLICATED or REVERSED and a VH class
other than UNDECIDED, the sign of s_A and the VH class are tabled side by
side: the same sign and class on every such chip is read as a carrier
anchored in the circuit, different signs or classes as one anchored in the
device; fewer than two such chips give no cross-chip reading. Everything else is read without a verdict; a single reading beyond
3σ is not a finding, and beyond 4σ it is reported as a sighting.

## 6. Flight conditions

The freezer is `simulations/flight/drive_frame_freeze.py`, the runner
`simulations/flight/run_drive_frame.py`; both are new files beside the
October pair and reuse its backend, readout-circuit and scheduling
functions. The runner refuses to submit on a chip unless all of these hold,
with no override flag: the chip's frozen circuit file's sha256 matches the hash
written here for that chip before its flight [TO BE FILLED AT EACH FREEZE:
Marrakesh, Kingston, Fez]; a path that ages out of its 24 h before
submission is re-frozen, and the new hash is a new commit of this page
before that chip flies; the frozen
file holds, for every turned variant, rz(φ) as the first operation on the
turned qubit and rz(−φ) after its basis change (asserted at the freeze and
again at submission); the chip's properties for the three path qubits are
under 24 h old; T₁ ≥ 100 µs and T₂ ≥ 80 µs on each; readout error ≤ 3 % and
P(1|0) ≤ 1 % on each; sx error ≤ 0.1 %; cz error ≤ 0.8 % on both bonds; the
forest qubits are read and not gated; at least 600 s of allocation left
before each submission; one job of 124 PUBs at 4096 shots in job mode (no
session), no dynamical decoupling, no twirling, init_qubits on, classified
measurement, max_execution_time 600 s; the calibration snapshot, the
shuffled order, the job id, and after the run the counts per circuit as two
halves and the billed usage from the job are written into the record under
`data/ibm_drive_frame_october2026/` (one file per chip). The
pre-registration (this page through §6, with the placeholder of §4 and the
hashes of §6 filled) is committed together with the freezer, the runner
and the frozen files before the first submission, and its commit hash is
written into §7.

## 7. The flight

Not flown.
