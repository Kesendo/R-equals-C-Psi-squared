# The driver's frame: who hears the shot. A hardware prediction, pre-registered

<!-- STAGE 0 SWEEP RECORD
Stores searched 2026-10-29 (evening) for drive crosstalk, a letter fixed in a qubit's drive frame, virtual phase as a lever, spectator readings, and the odd-Y excess's carrier.
experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md §7: the excess (s_A = +0.083, 6.7σ), every carrier tested there, and the cube reading that names a drive-anchored letter as the one coherent one-site letter untested; its "Not tested" list names crosstalk from a neighbour's drive and spectators outside the path. experiments/Z_N_PARTNERSHIP.md lines 85 to 93 name off-resonant drive bleed-through and ZX/ZY crosstalk as breakers of Z⊗N and line 104 predicts |M_X(a) + M_X(b)| of 1e-2 to 1e-1, not measured; experiments/IBM_CONCENTRATOR.md lines 280 to 319 assert that DD on one qubit affects its neighbour, not measured. data/ibm_zn_mirror_april2026/README.md lines 28 to 41 (job d7ornigror3c73c0c6ug, path [48, 49, 50]) read a transverse field h_y_eff ≈ 0.05 as a static field in the logical frame, registered as fw.Confirmations marrakesh_transverse_y_field_detection; a drive-anchored letter read through the accumulated virtual phase would look the same, and the README does not raise it. experiments/PRICE_PAIR_HARDWARE_PREDICTION.md line 153 measured ZZ on Marrakesh (ζ₀₁ = −3.918 ± 0.142 kHz, ζ₁₂ = −3.637 ± 0.101 kHz, path [93, 94, 95]) and lines 102 and 157 detuning offsets of 8 to 22 kHz. experiments/FIXED_POINT_SHADOW.md (2026-02-09, Torino Q52, the first IBM tomography): a late residual coherence of 0.0185 pointing one way in 17 of 17 samples, equal X and Y parts; its lines 260 to 295 weigh an assignment asymmetry against a basis-pulse leak through the Z extrapolation and leave them unranked; it does not cite the reality window and the window does not cite it; §2 joins them. No prior flight read a spectator outside the path or had a qubit idle while a neighbour was driven. docs/ANALYTICAL_FORMULAS.md: F131 (the order-sorting law) and F119 (the antilinear triangle) are the frame laws this page rests on; nothing on a virtual phase as a lever. OpenArcs registry: nothing on drive crosstalk ("Spectator" there is SpectatorIntertwinerClaim, a different object). docs/CAUGHT_ERRORS.md: nothing on virtual phases or drive frames. docs/GLOSSARY.md: the letter cube's row, nothing on frames. Typed layer: AntilinearTriangleClaim (the transport law), MirrorOrderSortingClaim (F131); no claim or witness on a drive frame. The parked frame echo (a virtual rz(φ) before and rz(−φ) after one cz) lives only in a memory file; around a diagonal gate with no pulse between, it is physically the identity, and it is replaced here by the turn of the driver's frame around its whole pulse train. The frequencies of the Heron qubits are published nowhere the repo holds (the April calibration CSVs in data/ibm_calibration_snapshots/ have no frequency column; the October properties none), so an always-on exchange cannot be fitted locally and is read here through the frame it shares with crosstalk. The three backends of the open plan on 2026-10-29 21:46: ibm_marrakesh (queue 27, no path fresh), ibm_kingston (queue 396, 35 fresh passing triples), ibm_fez (queue 1181, 9).
Adjacency checked: Z_N_PARTNERSHIP's bleed-through, the zn_mirror README's h_y_eff and REALITY_WINDOW §7 are about one object (a transverse letter fixed in a drive frame) and do not cite each other; FIXED_POINT_SHADOW and REALITY_WINDOW are about one class of reading and do not cite each other. This page joins them.
For the arms riding in these jobs, 2026-11-13: data/ibm_f83_signature_april2026/README.md (the X^⊗3 window on YZ+ZY, on Marrakesh [4, 5, 6] and on Kingston [43, 56, 63], data/ibm_soft_break_april2026/) and experiments/X3_WINDOW_HARDWARE_PREDICTION.md, which takes the frame turn to the second Hamiltonian; the Kingston path rule of §3 follows May's path.
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
fixed-phase reach of a FRAME-BLIND (below) presupposes the resonant model.

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
ibm_fez, one job each (on Marrakesh and Kingston beside the arms of the
X^⊗3 window page). The path on each chip is chosen at its freeze by one
rule: among the linear triples a-b-c (a < c) whose three qubits pass the
guards of §6 and whose middle b has exactly
two neighbours, the one with the largest smallest T₂ (a tie going to the
lower physical indices); on Marrakesh
[4, 5, 6] whenever it passes, so that A5(0) repeats the October arm, and on
Kingston [43, 56, 63] whenever it passes, the May path of
[the X^⊗3 window page](X3_WINDOW_HARDWARE_PREDICTION.md), whose arms ride
in the Marrakesh and Kingston jobs; should [43, 56, 63] fail the guards,
the rule above picks Kingston's path, the riding arms fly there and the
record says that the May path was not repeated. The ages of the
vendor's properties are recorded and not gated: the frame verdicts compare
the turns among themselves and need no calibration, the RO arm measures
the day's readout, and on Marrakesh and Kingston the colouring of the
X^⊗3 window page and arm Q read the day's coherence (on Fez no arm
does); the baseline's T₁ and T₂ enter V1 and W weakly (halving both on
the paths of 2026-11-13 moves S by 0.04σ_S on Marrakesh, 0.01σ_S on
Kingston and 0.003σ_S on Fez, and W's noncentrality by 1.3 at most), and
the simulated flights, the reaches and the drift check, drawn from the
same replay, are read as predictions; the record states how old the
properties were and when the job ran, and the path rule and the guards
read them whatever their age. The record names
which chips flew, and every verdict of §5 is per chip. The chips are set side by side without a verdict: a letter
anchored in the device (crosstalk with a chip's own δ) can turn on every
chip with its own size and sign, so a shared class does not tell device
from circuit. The April and May runs read the same sign of S on two chips.

The forest. Every circuit also measures, in Z, every coupling-map neighbour
of a, b and c that is not on the path (with b of degree two, these are the
off-path neighbours of a and of c, zero to two each; on Marrakesh [4, 5, 6]
they are physical q3 and q7). Each chip's forest, in ascending physical
index, is written into the freeze record. They receive no gate; each
carries an explicit delay from t = 0 to the path's last measurement so that
its measurement happens with the path's and not before (the October
freezer kept the schedule on the path and dropped the padding on every
other qubit; the new freezer keeps it on the forest as well). The RO circuits
measure them too. A measurement at the end changes nothing on the path,
read as the device's measurement model; a HEARD verdict is read beside the
even readings of the path for any sign of the contrary.

Per chip: 6 variants × 9 bases + 8 RO = 62 circuits, each submitted twice
at 4096 shots, 124 PUBs in one job, shuffled with seed 20261030, so that
drift between circuits shows as a difference between the halves. On Fez
that is the whole job. On Marrakesh and Kingston the job also carries the 78
circuits of the X^⊗3 window page's arms B and Q, and the 280 PUBs of the
two pages are shuffled as one with seed 20261113; those PUBs enter no
verdict of this page, and every circuit of them reads the forest as these
do. The October job of 196 PUBs billed 215 s and the April job of 36 PUBs
40 s, 1.1 s per PUB on both; that is about 136 s on Fez and 308 s on each
of the other two, about 750 s, under the 15 minutes set for the flight, of
the 4739 s the allocation reported on 2026-11-13 (read from IBM's usage record that day; the runner writes the day's remainder into the record before each submission; the remainder after the
October flight, no job having run since). The billed usage is read from the job and written into the
record.

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
the 99 % point and, where §5 calls TURNS, the 99.9 % point; both are given
with κ = 1 and with October's drift factors dividing λ, 1.27 over all its
fit arms and 1.96 over arm A's readings alone (both from the October
halves, [`x3_window_power.py`](../simulations/x3_window_power.py) part
7). The letters are yardsticks
for the statistics' reach; none is a prediction of what the chip carries.

- Under no letter and under the five Z letters: D = 0 up to the rounding
  above (on Marrakesh [45, 46, 47] up to 1 × 10⁻⁵, where a float sum of
  durations moves one of the replay's time boundaries by a rounding and
  exposes its step error); λ = 0, T₅ at its null (expectation 108).
- Under the fitted crosstalk alone (the eight parameters of §2): D reaches
  0.072, 0.104 and 0.126 at φ = π/2, π and 3π/2; λ = 1034, with drift 814;
  expected T₅ 1142 or 922 against the 99 % point 145.1 on 108 degrees of
  freedom; power 1.000 either way, and at the 99.9 % point 159.2 with
  either drift factor.
- Under all three families (seventeen parameters): 0.056, 0.060, 0.039;
  λ = 307, with drift 242; expected T₅ 415 or 350; power 1.000.
- The turn by π on q4: λ = 93 (crosstalk alone), 12 (all three); at the
  99.9 % point 68.0 on 36 degrees of freedom, power 1.000 and 0.046, with
  the drift factor 1.27 0.992 and 0.026, with 1.96 0.833 and 0.011; on q6:
  244 and 76, power 1.000 and 0.994, with 1.27 1.000 and 0.952, with 1.96
  1.000 and 0.651.
- The reach of §5: a crosstalk of equal size ε on (5 → 4) and (5 → 6) at the
  fitted phases, the other two pairs at zero, gives T₅ power at its 99.9 %
  point of 0.11, 0.72 and 0.996 at ε = 0.2, 0.3 and 0.4 % with the drift
  factor 1.27, and 0.03, 0.30 and 0.84 with 1.96; the 50 % point lies near
  0.27 % and 0.34 %. These use one drift factor; the reach §5 reports is
  computed with the drift matrix of §5, which divides even a drift-free
  chip's readings by more than one. Drawn with each simulated flight's own
  Σ and drift matrix on the paths frozen on 2026-11-13, its 50 % point lies
  at 0.285 % on Marrakesh [45, 46, 47], 0.281 % on Kingston [13, 14, 15]
  and 0.277 % on Fez [123, 136, 143].
- The simulated flights: 10 000 draws of each chip's job from the baseline
  replay's own distributions (multinomial at each circuit's shots, the RO
  circuits' counts drawn from the snapshot's RO matrix, which is then
  re-estimated in each draw from those counts), so that their mean is the
  baseline by construction; each file goes through the analysis unchanged
  and gives the null distributions of s_A/σ_S, T₅, T₄, T₆ and E (and those
  of the X^⊗3 window page); their 95 %, 99 % and 99.9 % points are written
  here beside the χ² points they check. They carry no drift, so they
  check the pipeline and the χ² points and not the drift matrix; the
  flights with drift are §5's. One further draw
  under the crosstalk-alone fit is run through the analysis as a check of
  the pipeline, not as a prediction. The flights are drawn at each chip's
  freeze by [`run_drive_frame.py`](../simulations/flight/run_drive_frame.py)
  `--nulls` (output in [`simulations/results/drive_frame/`](../simulations/results/drive_frame/)).
  At the freezes of 2026-11-13, 19:13 UTC, neither preferred path passed the guards
  of §6 (on Marrakesh [4, 5, 6] the P(1|0) of q4 read 2.15 %; on Kingston
  [43, 56, 63] the T₂ of q43 and q56 read 76 and 45 µs and the P(1|0) of
  q43 1.17 %), so the rule of §3 chose Marrakesh [45, 46, 47] (forest q37,
  q44, q48, q57), the one passing triple on that chip, and Kingston
  [13, 14, 15] (forest q12, q19), and on Fez [123, 136, 143] (forest
  q122, q124, q142, q144), whose properties were last measured on
  2026-11-09. The simulated points, 95 / 99 / 99.9 %:

  | statistic (degrees of freedom) | Marrakesh | Kingston | Fez | χ² |
  |---|---|---|---|---|
  | T₅ (108) | 111.5 / 122.1 / 134.4 | 112.5 / 122.1 / 138.2 | 111.9 / 122.6 / 134.4 | 133.3 / 145.1 / 159.2 |
  | T₄ (36) | 42.5 / 49.4 / 56.1 | 42.8 / 49.7 / 59.1 | 43.1 / 49.6 / 57.3 | 51.0 / 58.6 / 68.0 |
  | T₆ (36) | 42.6 / 49.8 / 56.9 | 42.6 / 49.2 / 57.0 | 42.8 / 49.5 / 57.6 | 51.0 / 58.6 / 68.0 |
  | E (36) | 42.7 / 49.1 / 57.5 | 42.7 / 50.3 / 62.0 | 42.9 / 49.5 / 58.5 | 51.0 / 58.6 / 68.0 |

  and the 99.7 % point of |s_A|/σ_S 2.70, 2.73 and 2.81 against 3. Every
  simulated point lies below its χ² point, so the χ² points rule: the drift
  matrix, whose eigenvalues are at least κ ≥ 1, makes every statistic
  conservative on a drift-free chip (T₅ averages 88 against 108). The
  forest's z, with the RO circuits' forest bits drawn as the path circuits'
  are, has mean 0.0 and spread 0.94 to 0.98 on every forest qubit with a nonzero
  P(1|0) in the snapshot; q57's is zero there, and with no flip on either
  side it gives no verdict. The draw under the crosstalk-alone fit carried
  by role reads T₅ 902 on Marrakesh, 1151 on Kingston and 964 on Fez,
  TURNS on T₅, T₄ and T₆, and V1 ABSENT on the first two and INCONCLUSIVE
  on Fez.

## 5. Verdict rules

The counts are read as the October runner reads them, with the forest bits
appended: the count key's bits are (q0, q2, q1, forest in ascending
physical index) from least to most significant, as the freeze record pins
them; the path's parities are formed from the first three bits with the
forest bits summed out, and the RO matrix is the 8 × 8 matrix of the path
bits with the forest bits summed out. The **baseline** of a chip is the
October replay run on that chip's path under its freeze snapshot (T₁, T₂)
with the RO matrix its RO arm measures, every letter at zero. Every path
reading is a parity of measured bits in one basis circuit, both halves
together. Σ is the multinomial covariance of one variant's 36 odd-Y
readings at the circuit's shots, computed from the measured parities
pooled over the standing A variants, which are logically one circuit and
share their parities under the null (the replay's parities would carry
its contrast, and October's even readings stood at 0.72 to 0.94 of it);
its rank is the number of eigenvalues above 10⁻⁹ of the largest, 36 unless
a parity is ±1. Σ has no covariance between two basis circuits, and every
PUB runs at its own shuffled position in the job, so a drift moves the
readings of one circuit half together and those of two circuits
independently. The drift is weighed block by block: in each basis
circuit's block Σ is whitened and S_b is the mean outer product of the
standing variants' whitened half differences (half₁ − half₂)/2, which have
covariance Σ when nothing drifts; κ is max(1, the sum of the traces of the
S_b over the rank); K_b is S_b with every eigenvalue below κ raised to κ;
and every covariance below written κΣ is Σ_b^½ K_b Σ_b^½ on each block and
zero between blocks. With S_b diagonal this is max(κ, κ_i) reading by
reading; as a matrix it takes a drift that moves the readings of one
circuit together, as a wandering frequency does, on A's correlated
readings (Σ's eigenvalues spread over a factor of 21 to 46), where
reading-by-reading weights fall behind (below). σ_S = √(κΣ_SS + Σ_RO,SS), Σ_RO the noise of the day's RO
matrix carried into the baseline (column by column, multinomial, J·C·Jᵀ,
as [the X^⊗3 window page](X3_WINDOW_HARDWARE_PREDICTION.md)'s §5 carries
it), times the RO arm's own drift factor κ_RO, max(1, the χ² of its eight
circuits' half differences against their multinomial variance over its
degrees of freedom); it cancels in every difference between variants. The predictions of
§4 give κ = 1, 1.27 and 1.96. The flights of §4 carry no drift;
[`drive_frame_drift_check.py`](../simulations/drive_frame_drift_check.py) (output for
[Marrakesh](../simulations/results/drive_frame/ibm_marrakesh_drift_check.txt),
[Kingston](../simulations/results/drive_frame/ibm_kingston_drift_check.txt) and
[Fez](../simulations/results/drive_frame/ibm_fez_drift_check.txt)) draws 1000 flights per
case on each chip's frozen path through this analysis with what it was not
built from: the October library's letters (detunings, amplitudes, the
phases after the pulses and after the cz gates, the cz's conditional phase,
the spectator's phase at a cz, static ZZ) wandering together from PUB to
PUB at sizes that give κ up to 9 on A and 47 on B, and one at a time (a
detuning of 3 to 16 kHz on a or on b, the phase on b after the b–c cz);
the readout's P(1|0) and P(0|1) wandering from PUB to PUB with a spread of
up to 0.02; the state depolarized by 20 to 34 % against the baseline, with
and without drift; and the forest at flip rates from 2 × 10⁻⁵ to 10⁻⁴ or
at 5 % wandering by 20 %. Each rate below is a count among 1000 flights
(0.9 % is 9 of them), the largest over the cases drawn on the three chips. On arm A, T₅, T₄, T₆ and E pass their 99.9 %
points by chance in 0.3 % of the flights or fewer; on arm B, the X^⊗3
window page's T₅ᴮ, T₄ᴮ, T₆ᴮ and Eᴮ, whose blocks hold one to three
readings, in 0.9 % or fewer; every frame statistic passes its 99 % point
in 1.9 % or fewer, which reads INCONCLUSIVE. V1 reaches 3σ_S in 1.5 % or
fewer, and in up to 2.9 % (REPLICATED alone, 1.5 %) when a single letter's
drift falls on S, as the phase on b after the b–c cz does, wandering by
0.12 rad from PUB to PUB: S's block weighs its drift from six half
differences, and z carries their tail. W passes its 99.9 % point in 0.7 %
or fewer and its 95 % point, below which alone SHUT is called, in up to
5 % with the readout wandering by 0.02. VQ calls DETUNED by chance on one
of the three qubits in 1.3 % of the flights or fewer, a single qubit's
wandering detuning included. VS calls HEARD in none. Reading-by-reading
weights on the same flights let T₅ pass its 99.9 % point in up to 0.8 %
and its 99 % point in up to 2.6 % with the letters at κ_A ≈ 3, where the
block form gives 0.1 % or less at both. The frame
statistics call TURNS only above their 99.9 % points. All verdicts
are per chip; the χ² points below rule, and the simulated flights of §4
are a check on them: for each of T₅, T₄, T₆ and E, if a simulated 95 %,
99 % or 99.9 % point (numpy's default quantile) exceeds its χ² point by
more than 5 %, the simulated point rules for that statistic and the record says so; for s_A the same with
the simulated 99.7 % point of |s_A|/σ_S against 3.

- **P, the parse gate.** A5(0): ⟨X₀Z₂⟩ < −0.5, −0.75 < ⟨Z₀X₂⟩ < −0.3 and
  ⟨X₀Z₂⟩ < ⟨Z₀X₂⟩ − 0.05. That is October's gate with the margin of its
  third condition narrowed from 0.15: on the paths of 2026-11-13 the
  baseline's gap is 0.23 and 0.25, and with 0.15 the drift check above
  voids a chip at a contrast of 0.78, 0.70 and 0.66 in 1, 22 and 48 % of
  flights on Marrakesh (0, 6 and 17 % on Kingston, 0.1, 4 and 21 % on
  Fez), with 0.05 in none,
  while an exchange of q0 and q2, which the gate is there to catch,
  reverses the difference and moves ⟨Z₀X₂⟩ out of its range. A chip that
  fails it is VOID, and that is accepted. The gate is not applied to the
  other variants (a dropped rz would pass it); their even readings are
  reported. The circuits retrieved from the job after the flight are
  compared with the frozen ones gate for gate on the path, and a variant
  with any difference, a missing turn included, is VOID. They are the
  circuits the client submitted as IBM stores them, so a turn dropped
  below them is not seen there (FRAME-BLIND covers it, below); if they
  cannot be read, the record says so and the variants stand on the
  structure asserted at submission.
  A guard of §6 that the RO arm's own measurement breaches on the day (as
  q4's P(1|0) did in October; P(1|0) and the readout error, the mean of
  P(1|0) and P(0|1), per path qubit) is reported with the verdicts, which stand
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
  - T₅ above its 99.9 % point: TURNS. The odd-Y readings carry a
    component that turns with b's frame, a foreign transverse letter
    coupling to b; whether S carries it is not decided (the share, below).
  - T₅ at or below its 99 % point: FRAME-BLIND at this sensitivity for
    letters coupling to b; between the two points, INCONCLUSIVE. FRAME-BLIND
    also covers a turn the control stack dropped below the submitted
    circuit, which no logical reading can show. The reach reported with it:
    the ε, equal on (b → a) and (b → c) at the phases §2 fitted
    for (5 → 4) and (5 → 6), a fitted letter's phases used as a yardstick
    and not as a prediction, and with (a → b) and (c → b) at zero, at which
    T₅, formed with the day's κΣ, passes its 99.9 %
    point in half of 2000 simulated flights drawn under that letter from the
    baseline replay (a 50 % point, not a bound: at that ε a share of the
    flights still reads FRAME-BLIND), computed with the chip's
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
    99.9 % point on rank(Σ) degrees of freedom is reported as 2φ PRESENT, the
    cos 2φ signature of an IQ image term of b's own drive, read beside
    TURNS and not as its cause.
- **The share, a reading.** When all four A5 variants stand, S̄ (the mean
  of S over the four) minus the baseline's S, with σ = √(κΣ_SS/4 +
  Σ_RO,SS), and S(0) − S̄, with σ_t = √(κΣ_SS)·√3/2, are reported without a
  verdict: a transverse letter puts part of itself into S̄ at second order,
  so S̄ is not a frame-blind share.
- **VS, the forest.** For each forest qubit, p_A is its P(1) over all A
  circuits of the chip (every variant and basis, both halves). Its null is
  not an idle one: the RO circuits fire x pulses on a and c, so they shoot
  too, at a lower dose and for a shorter time, and their outcomes depend on
  the path's state, so the null is p₀ = Σ_s P_A(s) · P_RO(forest bit = 1 |
  s), with P_RO(· | s) read from the RO circuit that prepares path state s
  and P_A(s) the path-outcome frequencies pooled over the A circuits. σ²
  is the sum of the binomial variances of both sides taken at the null's
  rates, with each RO circuit's rate taken as (ones + ½)/(n + 1) (and p₀
  formed from those on the A side), so that an RO circuit in which the
  qubit never flips is not given a zero variance, times the forest qubit's drift factor, max(1, the mean over
  every A and RO circuit of (half₁ − half₂)²/(var₁ + var₂) at that
  circuit's pooled rate). The correlation between p_A and P_A(s) from
  shared shots is neglected, and the forest's longer idle in the A circuits
  is not modelled: both are said with the verdict. With no forest flip on
  either side no verdict is given. p_A − p₀ ≥ 5σ: HEARD, the forest qubit holds
  population beyond what the RO circuits' pulses and its readout crosstalk
  account for. |p_A − p₀| < 3σ: SILENT at this sensitivity. Between, or
  ≤ −5σ: reported without a verdict.

| VF | reading |
|---|---|
| TURNS | the odd-Y readings carry a component turning with the middle qubit's frame, a foreign transverse letter coupling to it: a neighbour's pulse reaching it or its pulse reaching a neighbour, an always-on exchange, cz leakage or residual swap, a local-oscillator leak, or another lab-fixed transverse letter; VS narrows it in part |
| FRAME-BLIND | for letters coupling to the middle qubit, the excess is carried at this sensitivity by Z letters, own-frame letters, or something outside the qubits' letters (readout, own-pulse leakage); b's pulses reaching the ends are below the reported reach in the sense of its 50 % point; T₄ and T₆ say the same for letters coupling to the ends, for the cosine part alone |
| INCONCLUSIVE or VOID | no class |

Across the chips that flew with V1 REPLICATED or REVERSED, the sign of s_A
and the VF classes are tabled side by side without a verdict: a carrier
anchored in the device can read TURNS on several chips with different
sizes and signs. Everything else is read without a verdict; a single reading beyond
3σ is not a finding, and beyond 4σ it is reported as a sighting.

## 6. Flight conditions

The freezer is `simulations/flight/drive_frame_freeze.py`, the runner
`simulations/flight/run_drive_frame.py`; both are new files beside the
October pair, whose backend, readout-circuit and scheduling functions they
carry over, generalised to a chosen path and its forest; on Marrakesh and Kingston the freezer also builds the X^⊗3
window page's arms B and Q into the same frozen file and asserts their
turns as these and the idles of Q. The runner refuses to submit on a chip unless all of these hold,
with no override flag: the chip's path, its frozen circuit file's sha256
and its freeze record's sha256 (of the record's canonical JSON, so that no
line ending moves it) are the ones in the table below; no submission from
that frozen file is on record already (one is recovered with `--recover`,
not resubmitted); the target still carries every operation of the frozen
circuits (a forest qubit gone faulty aborts the chip) at its dt; this page, the X^⊗3 window page, the freezer, the
runner, the replay, the frozen file, its freeze record and its file of
simulated flights (`drive_frame_<chip>_nulls.json`, the ruling points of
10 000 flights drawn for that circuit file and that freeze record) are
committed and unchanged, and the drift check of §5 names that freeze
record; the target's gate and measurement durations on the
path and the forest are the freeze's; the frozen
file holds, for every turned variant, rz(φ) as the first operation on the
turned qubit and rz(−φ) after its basis change (asserted at the freeze and
again at submission); on the properties read at submission, whatever
their age, which is recorded: T₁ ≥ 100 µs and T₂ ≥ 80 µs on each; readout error ≤ 3 % and
P(1|0) ≤ 1 % on each; sx error ≤ 0.1 %; cz error ≤ 0.8 % on both bonds; the
forest qubits are read and not gated; at least 600 s of allocation left
before each submission; one job of 124 PUBs (280 on Marrakesh and
Kingston) at 4096 shots in job mode (no session), no dynamical decoupling, no twirling, init_qubits on, classified
measurement, max_execution_time 600 s; the calibration snapshot, the
shuffled order and the job id are written before the run; after it the
counts per circuit as two halves are written first, then the circuits
retrieved from the job and their comparison, the result's metadata (the
execution spans), the job's timestamps and its billed usage, into the
record under `data/ibm_drive_frame_<month><year>/`, the month the flight
flies in, which is the month of its freeze, written at the freeze (one
file per chip). The frozen files, 2026-11-13, 19:13 UTC:

| chip | path | sha256 of the frozen circuit file | sha256 of the freeze record (canonical JSON) |
|---|---|---|---|
| ibm_marrakesh | [45, 46, 47] | `fcee7bd194e4fcf63c366654792492291fa53435b20cf97faf721d248e92a30c` | `a65fc571102e98ca59a8b0e9b5a315c112b117c52c2cb41593e9d177c66a25a5` |
| ibm_kingston | [13, 14, 15] | `956cc064521f46f808f0cc461c9d8421e564b042facb5927f99ef51972bd2093` | `04727b56cc92d738d703d4f410e304eb46c63a87a57901352924571a53610a73` |
| ibm_fez | [123, 136, 143] | `df61e1b82a9a8e18fdbae1f778ec8a89540487046ee742751611b72120e22de1` | `90e1852a4b13b8143ee0b04d31998ca3061e3eb1bfc4240dcd08606e0aabd9d1` |

A flight from another freeze runs `--nulls`, `--freeze-numbers` and the
drift check again on it and rewrites this table, §4's numbers, §5's rates
and the parse gate's figures, and the X^⊗3 window page's §4 and its count
of verdicts per chip from them; the runner refuses a chip whose drift
check does not name its freeze record. The
pre-registration (this page through §6, with §4's numbers and §6's
hashes those of the freezes that fly) is committed together with the
freezer, the runner and the frozen files before the first submission, and
its commit hash is written into §7.

## 7. The flight

Not flown.
