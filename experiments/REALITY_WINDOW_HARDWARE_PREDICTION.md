# The reality window: a hardware prediction (pre-registration)

<!-- STAGE 0 SWEEP RECORD
Stores searched 2026-10-28 for Y-parity, n_Y parity, reality, complex conjugation as a symmetry of the generator, and any prior reading of odd-Y observables on hardware.
docs/ANALYTICAL_FORMULAS.md: F119 (the antilinear triangle, docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md lines 43 to 45: conj(A) carries (−1)^{n_Y}, θ is the mirror D = diag((−1)^{n_Y})); line 6879 notes that σ⁻ and σ⁺ are real. Neither is read as a noise detector.
experiments/: MARRAKESH_THREE_LAYERS.md line 92 states the rule for the 2026-04-26 run (XY+YX keeps a real ρ₀ real, every odd-Y string vanishes, the Y-parity fw.y_parity_panel reads); line 95 reports the predicted zeros reading up to 0.098 for soft on hardware; line 198 lists soft ⟨Y₀Z₂⟩ = +0.098 and ⟨Y₀⟩ = −0.079 among residuals "what carries them we do not know".
data/: ibm_f83_signature_april2026/README.md lists pi2_odd ⟨Y₀Z₂⟩ = +0.067 and ⟨Z₀Y₂⟩ = +0.053 as "what carries those we do not know" and reads pi2_even's ⟨Y₀⟩, ⟨Y₂⟩ only through X^⊗3; ibm_soft_break_april2026/README.md line 78 lists soft ⟨Y₀Z₂⟩ = +0.098 without comment, and its lines 66 to 67 say continuous Lindblad, Trotter and hardware agree in sign on ⟨Y₀Z₂⟩, ⟨I₀Y₂⟩, ⟨X₀Y₂⟩, which for soft cannot hold (the models give exactly zero there): a sentence to repair, not repaired here. The Kingston 2026-05-05 JSON (job d7sqjpiudops73976960) has no README reading its odd-Y values.
simulations/: f83_hy_field_check.py tested an h_y field, which keeps every imaginary Hamiltonian imaginary and so cannot light an odd-Y reading in three of the four classes; no Z-detuning, ZZ or per-gate model was fitted to these runs before this document; framework/symmetry.py y_parity_panel builds L from Z-dephasing and T1 only and reports Y-parity preserved by construction.
fw.Confirmations: f83_pi2_class_signature_marrakesh (the 2026-04-30 run and its Kingston reproduction) and regime_uniformity_kingston_uniform_quantum (the same Kingston run through the F88b lens) register these runs; marrakesh_transverse_y_field_detection (h_y_eff ≈ 0.05, Z⊗N mirror test, path [48,49,50]); none reads Y-parity.
OpenArcs registry: no arc on the F83 X^⊗3 window or on Y-parity on hardware; line 2591 names Y-parity as a remaining calibration instrument.
docs/CAUGHT_ERRORS.md: the 2026-10-28 Trotter-order entry records the YZ+ZY window through X^⊗3; line 4003 (a detuning does not freeze a phase) is the nearest precedent on detuning readings.
Typed layer: the eight YParity-axis claims in compute/RCPsiSquared.Core/Symmetry (F102, F103, F105, F106, F107, F109, F110, F111; HardCellYInversionPattern.cs line 145 lists them) classify Hamiltonian terms by y_par, MotherSoftYParityOnePurity being the typed home of the premise that XY+YX carries odd y_par; none has dynamics or a hardware reading; no Diagnostics witness reads odd-Y observables.
docs/GLOSSARY.md: Y-parity appears only as a term label.
Adjacency checked: the F83 README, MARRAKESH_THREE_LAYERS and the soft_break README describe the same zero and do not cite each other; the h_y rejection and the reality rule are about the same object and do not meet. This document joins them.
-->

This flight asks one question of the IBM hardware, and the letter cube asks
it: where does the reading come from that the circuit forbids?

## 1. The zero

Three of the four Hamiltonians flown on 2026-04-30 (XY+YX, YZ+ZY, XY+YZ)
carry an odd number of Y letters in every term. Each such H is purely
imaginary, every Trotter factor exp(−iδt·P) is a real matrix, the initial
state |+−+⟩ is real, and every noise channel in the model (Z dephasing, σ⁻,
σ⁺, any Pauli channel) maps real density matrices to real ones. A real ρ
has zero expectation on every Pauli string with an odd number of Y letters,
because those strings are imaginary antisymmetric matrices. On the letter
cube this is a weight parity: a Hamiltonian letter P with n_Y(P) odd keeps
the real (even-Y) sector closed, a letter with n_Y(P) even (a Z field, a ZZ
coupling, an X field) carries a state out of it, and a real jump (Z
dephasing, σ⁻, σ⁺, any Pauli channel) never does.

So for XY+YX an odd-Y reading measures a real-breaking letter and nothing
else, once the letters are read in the logical frame. Physical noise acting
inside the compiled circuit's intermediate frames is conjugated by the gates
around it, and a cooling jump there is no longer the logical σ⁻; with the
day's T1 and T2 this puts a floor of about 0.002 under the zero (§4).
Readout is not real-breaking but it is not invisible either: an asymmetric
readout error shifts a measured single-bit parity by P(0|1) − P(1|0) and
leaks even readings into odd ones, so the readout is measured in the job
and modelled (§3, §5).

## 2. What is already measured

The same observable reads nonzero three times, on two chips, with one sign
(⟨Y₀Z₂⟩ and ⟨Z₀Y₂⟩ from their own basis; ⟨Y₀⟩ averaged over the three
bases that carry it):

| run | job | path | ⟨Y₀Z₂⟩ | ⟨Z₀Y₂⟩ | ⟨Y₀⟩ |
|---|---|---|---|---|---|
| Marrakesh 2026-04-26 | d7mjnjjaq2pc73a1pk4g | [48, 49, 50] | +0.098 | +0.017 | −0.076 |
| Marrakesh 2026-04-30 | d7pol1e7g7gs73cf7j90 | [4, 5, 6] | +0.067 | +0.053 | −0.030 |
| Kingston 2026-05-05 | d7sqjpiudops73976960 | [43, 56, 63] | +0.066 | +0.032 | −0.032 |

The circuits actually flown were retrieved from IBM on 2026-10-28 and are
stored beside their data (`data/ibm_f83_signature_april2026/flown_circuits_d7pol1e7g7gs73cf7j90.qpy`
and the two files `flown_circuits_d7sqjpiudops73976960.qpy`,
`flown_circuits_d7mjnjjaq2pc73a1pk4g.qpy` in `data/ibm_soft_break_april2026/`).
All three runs flew the same XY+YX gate sequence, op for op on their paths,
and on 2026-04-30 the 36 flown circuits implement the runner's first-order
Trotter circuits exactly (outcome probabilities equal to 1.4e−14). They
carry no delays: where IBM placed each gate in time was the backend's
choice and is not recorded.

[`simulations/reality_window_april_letters.py`](../simulations/reality_window_april_letters.py)
replays the flown circuits gate by gate with the day's T1, T2 and readout
and reads every class through its forbidden sectors, once with the gates
as late as possible and once as soon as possible (q1, which April did
not read, loses its gates after its last cz, single-qubit gates on q1
alone after its last interaction). With ideal gates, ⟨Y₀Z₂⟩ comes out
−0.002 under either schedule against the measured +0.067, and the
forbidden readings give χ² = 24, 50, 70 and 12 under the late schedule
(24, 52, 71, 12 under the early one; truly, XY+YX, YZ+ZY, XY+YZ, over 8,
6, 10 and 6 readings). All nineteen physical letters of §5, fitted jointly
to those 30 readings, bring the total to 20 (late) or 19 (early) over 11 remaining
degrees of freedom, but only with a detuning of about −2.2 rad/µs on the
middle qubit (some 350 kHz, twenty times its bound in §5) and further
letters beyond their bounds; even so they leave ⟨Y₀Z₂⟩ at +0.035 under
either schedule, 2.1σ short of the measured value (σ = 0.0156 at
April's 4096 shots). Static ZZ alone, fitted
to XY+YX alone, gives χ² 39 and ⟨Y₀Z₂⟩ −0.015 under the late schedule,
and χ² 18 and +0.024 at −0.72 rad/µs on pair 4-5 under the early one.
What carries the XY+YX reading we do not know.

## 3. The arms

All arms run XY+YX on N = 3 from |+−+⟩, J = 1, δt = 0.8/3, nine Pauli
bases on (q0, q2), on ibm_marrakesh path [4, 5, 6]. In every arm q1 is also
read, in the Y basis (sdg then h), so that its own odd-Y readings exist.

| arm | circuit |
|---|---|
| A | the 2026-04-30 flown circuits, gate for gate on every qubit, then q1's readout |
| R | the same gates with physical 4 and 6 exchanged |
| D0, D1, D2, D3, D4, D6 | the logical circuit at n_trotter = 0, 1, 2, 3, 4, 6, compiled 2026-10-28 (level 3, seed 11, layout [4, 5, 6]) |
| K+, K− | arm A with a virtual rz(±0.2) on q1 right after its sixth cz |
| RO | the eight computational states of (q0, q2, q1), prepared with x gates, for the readout matrix |

Every circuit is scheduled as late as possible against the ibm_marrakesh
target (sx and x 36 ns, cz 68 ns) and carries the schedule as explicit
delays on the three path qubits, so the backend has no timing left to
choose and the replay follows the same timeline. The 98 circuits are
frozen by [`simulations/flight/reality_window_freeze.py`](../simulations/flight/reality_window_freeze.py).
Each is submitted twice at 4096 shots, all 196 in one job, in an order
shuffled with seed 20261028, so that drift between circuits shows as a
difference between the halves.

A virtual rz is a frame change in the control electronics, exact on IBM
hardware, so K± place a Z letter of known sign and size at a known point.
It is a control of the sign chain and the pipeline, not of the mechanism:
a frame change is exactly what a physical error is not.

Arm A is not April's circuit in every respect: its gates carry an explicit
schedule where April's had the backend's, and reading q1 changes its final
state (April left it with ⟨Z₁⟩ ≈ −0.6) and so any readout crosstalk that
depends on q1. The single-site readings are therefore not comparable with
April; S, a two-body reading of q0 and q2, would need a large crosstalk to
move. D3 and A are the same
logical circuit compiled twice, and the two compiles differ in several
single-qubit gates on every qubit and in where x pulses sit, by basis.

The arms are not separate verdicts. Each physical letter leaves its own
pattern across all of them: R moves letters bound to an end qubit or a CZ
pair and leaves the middle qubit's letters in place, D changes where each
gate sits in time, D0 is blind to detuning and to every Z-phase letter but
not to amplitude errors. The pattern is the measurement.

## 4. Predictions from below

[`simulations/flight/run_reality_window.py`](../simulations/flight/run_reality_window.py) `--gate`
replays the frozen circuits with ideal gates under the calibration of
2026-10-28 (T1 = 184, 202, 226 µs and T2 = 175, 146, 225 µs on q4, q5, q6;
readout P(1|0) = 0.0083, 0.0005, 0.0000 and P(0|1) = 0.0308, 0.0073, 0.0083).
The replay is time ordered: the schedule fixes every gate's start; sx and x
are finite 36 ns pulses during which detuning and static ZZ keep acting, so
a large letter turns the pulse as it would on the device; rz is a frame
change at its instant, and a cz acts at its start (the coherent terms on
its two qubits during it are diagonal and commute with it; decay does not,
and moving every cz to its end changes an odd-Y reading by less than
1e−3); decay acts throughout. The gate checks it against a second
evaluation, a pure-state one written apart from it (finite pulses, cz at
its midpoint), in 36 cases up to 1.3 rad/µs of detuning or ZZ: they agree
to 5e−12. That route shares the schedule extraction and the elementary
matrices and carries neither decay nor the letters other than detuning and
ZZ; it catches errors in the time slicing and the pulses, not in those.

- Without readout, every odd-Y reading of A, R and D0 to D6 is at most
  0.0023, a fifth of the shot noise (the gate holds it below σ/4): this is
  the floor T1 leaves under the zero of §1. With the readout forward model the asymmetry leaks into odd-Y
  readings, up to 0.034 (a three-body reading of D4) and 0.022 on q4's
  single bit, and S = ⟨Y₀Z₂⟩ picks up at most 0.006 from ⟨Z₂⟩; this is the
  baseline subtracted below.
- K±: S = +0.131 and −0.142; the control C = +0.137, band (0.087, 0.160).
- Arm A's even readings ⟨X₀Z₂⟩ = −0.836 and ⟨Z₀X₂⟩ = −0.595 (April
  measured −0.849 and −0.566).
- A simulated flight: Qiskit Aer runs all 98 frozen circuits with ideal
  gates and the snapshot's readout errors on the right physical qubits,
  writes them in the hardware's two-half shuffled format, and the file goes
  through the analysis unchanged: parse gate PASS, control PASS
  (C = +0.141), V2 NO SIGNAL, V1 INCONCLUSIVE (s_A = +0.029, 2.4σ, and D0
  at 2.5σ, the largest arm, in one draw of null data). This checks
  the bit order from Qiskit's count keys to the verdicts by a second route,
  and the whole analysis on data of the shape the flight returns.
- Letter power: each of the nineteen letters placed alone through the full
  replay at the smallest size, searched on a grid of both signs, that
  brings s_A to +0.067 and again to +0.033, with the readout matrix
  re-estimated from simulated calibration counts and every letter refitted
  through the replay, one draw per scenario. At +0.033 all nineteen come
  out LETTER FOUND with the true letter best and none unseparated (for
  detuning on q6, which reaches only +0.015, at that size). Eighteen
  letters reach +0.067 (the sx phase on q4 stops at +0.033); of these,
  fifteen come out LETTER FOUND with the true letter best and none
  unseparated, static ZZ on pair 4-5 (at −1.37 rad/µs) and the
  conditional phase on pair 5-6 (at 1.74 rad) miss the sufficing line
  narrowly (χ² 344 and 342 against 339, the true letter still best:
  SIGNAL, NO SINGLE LETTER), and detuning on q6 gets there only at
  5.8 rad/µs, where its fit fails.
- The sizes these single letters need are themselves a reading: to reach
  +0.067 alone, every letter but one needs a size beyond its bound in §5.
  The exception is the local Z phase on q5 after each cz of pair 4-5
  (0.034 rad); at +0.033 the local phase on q5 from pair 5-6 and both
  spectator phases also come within their bounds.
- Errors outside the library, each scaled toward s_A = +0.067: a tilt of the
  sx pulse axis on q4 or q5 moves no odd-Y reading at all (it is a choice
  of frame, which the reality window does not see); detuning on q5 with
  static ZZ on 4-5 in equal parts, the local phase on q5 from both pairs,
  and an equal amplitude error on all three qubits (the first reaching
  +0.060, the others +0.067) each come out SIGNAL, NO SINGLE LETTER, the
  composite read as composite and not as one of its parts.
- What the power runs leave out: no drift (κ = 1) and no bootstrap term
  for the readout matrix's own noise, although the matrix is re-estimated
  in each run. Carrying is then slightly easier and sufficing slightly
  harder than under the rule of §5, so the runs bound the rule only
  roughly, in neither direction exactly.

Shot noise: a two-body reading from one basis has σ = 1/√8192 = 0.0110
before the drift factor.

## 5. Verdict rules

The counts are read as the frozen runner reads them: count key c2c1c0 with
bit 0 = q0, bit 1 = q2, bit 2 = q1. Every reading is a parity of measured
bits in one basis circuit, from both halves together. Its covariance with
the other parities of the same circuit is the multinomial one computed from
the baseline's predicted parities (not from the counts, where a nearly
unanimous bit would carry a variance estimated from a handful of flips).

- **Drift.** The factor κ = max(1, mean of (half₁ − half₂)² / (var₁ +
  var₂)), over the 288 odd-Y parities of the fit arms with nonzero
  variance. Every σ below is σ·√κ, and V2's covariance is multiplied by κ.
  A drift that is the same in both halves does not show here. κ is an
  estimate: without any drift it still comes out above 1 by chance (1.14
  in the simulated null flight), which makes the rule more conservative,
  not less.
- **The baseline.** The replay of §4 under the calibration snapshot the
  runner records before submitting, with the readout matrix measured by
  arm RO (the full 8 × 8 matrix, so correlated errors and crosstalk between
  the three readouts are in it; in arm R bits 0 and 1 are exchanged, since
  there q0 sits on q6). The matrix's own shot noise is carried into V2:
  the RO counts are resampled 200 times, the baseline rebuilt from each,
  and the covariance of those baselines is added to the readings'. V0 and
  V1 do not carry this term: its size there is about 0.0005 against
  σ = 0.011.
- **P, the parse gate.** Arm A: ⟨X₀Z₂⟩ < −0.5, −0.75 < ⟨Z₀X₂⟩ < −0.3 and
  ⟨X₀Z₂⟩ < ⟨Z₀X₂⟩ − 0.15. A swap of q0 and q2, or a parse that lands on q1,
  fails it. FAIL makes V1 and V2 VOID.
- **V0, the control.** C = (S_K+ − S_K−)/2 against its baseline prediction
  C_pred, σ_C = σ/√2. C ≤ −3σ_C: SIGN CHAIN BROKEN. |C| < 3σ_C: CONTROL
  DEAD. Both make V1 and V2 VOID. C in [0.75·C_pred − 2σ_C,
  C_pred + 3σ_C]: PASS. Otherwise MAGNITUDE OFF, and V1 and V2 stand. The
  lower edge is the worst contrast either April Marrakesh run kept against
  the noiseless Trotter model on a large reading (0.75, ⟨Z₀X₂⟩ of XY+YX on
  2026-04-26); the readout is already inside C_pred. The half-sum of K±
  minus its baseline is a second estimate of s_A, reported.
- **V1, replication.** s = S_A minus its baseline. s ≥ 3σ: REPLICATED.
  s ≤ −3σ: REVERSED. |s| < 2σ: ABSENT. Otherwise INCONCLUSIVE. The April
  letter fit of §2 reached +0.035, and only with
  letters beyond their bounds; a REPLICATED s is read beside that.
- **V2, the letters.** The readings are every odd-Y parity of every circuit
  of A, R, D0 to D6: 288 readings; the dof is the rank of their covariance
  on the day (282 in the simulation, where three D0 bits are unanimous;
  288 when they are not). The library, with the size each letter can
  plausibly have on a calibrated Heron (an order-of-magnitude bound, fixed
  here):

  | letter | where | bound |
  |---|---|---|
  | Z detuning | q4, q5, q6; during every gate and delay | 0.1 rad/µs (16 kHz) |
  | sx/x amplitude error | q4, q5, q6 | 1 % of the rotation |
  | Z phase right after each sx/x pulse | q4, q5, q6 | 0.02 rad |
  | local Z phase after each cz | each qubit of each pair | 0.05 rad |
  | CZ conditional phase (an rzz angle κ, the phase error is 2κ) | pairs 4-5, 5-6 | κ ≤ 0.05 rad |
  | static ZZ (rate of the rzz angle; in the ζ convention, twice this) | pairs 4-5, 5-6 | 0.1 rad/µs |
  | Z phase on the spectator at each cz | q6 at 4-5, q4 at 5-6 | 0.05 rad |

  For each letter one amplitude is fitted through the full replay
  (weighted least squares from four starts: the linear estimate, half and
  twice it, and its negative; the best converged fit is kept, and a letter
  none of whose starts converged is kept at its best value and flagged),
  giving χ²_j; χ²₀ is the residual with no letter, and b the letter of
  smallest χ². The fit is nonlinear because the sizes a single letter
  needs to reach the April reading lie outside the linear range.
  - The letter b carries when χ²₀ − χ²_b ≥ 25.
  - χ²₀ ≤ the 99 % point of χ² at the rank and b does not carry: NO SIGNAL.
  - b carries and χ²_b is at most the 99 % point at rank − 1: LETTER FOUND.
    b carries and χ²_b is above it: SIGNAL, NO SINGLE LETTER. χ²₀ above
    its point and b does not carry: SIGNAL OUTSIDE THE LIBRARY.
  - Every letter j with χ²_j − χ²_b < 9 is listed as not separated from b.
  - b's fitted size is reported against its bound as plausible or
    implausible; the verdict does not depend on it, the reading does.
  - When the linear estimate of a letter is exactly zero, its fit starts
    from ± five times the letter's finite-difference step instead.
  - Static ZZ runs through the cz's own 68 ns as well; a real device folds
    that part into the cz's calibration, so between them the zz and czk
    letters share it: which of the two the fit names is not read.

The K arms are excluded from V2: their large controlled signal would weigh
the fit by the contrast the model does not carry.

What each outcome means. The instrument comes first: a parse gate FAIL,
or a control SIGN CHAIN BROKEN or CONTROL DEAD, voids V1 and V2. It says
the instrument or the day failed, not what the device does to the reality
window: the parse gate also fails if arm A's contrast falls below about
0.6, a device far off its calibration. A control
MAGNITUDE OFF says the sign chain holds and the contrast is not the
model's; V1 and V2 stand, with that said. Then V1 and V2 are read
separately and their readings put side by side:

| V1 | reading |
|---|---|
| REPLICATED | April's ⟨Y₀Z₂⟩ is present today, with April's sign |
| REVERSED | present with the opposite sign: a finding in itself, since April's sign held on two chips |
| ABSENT | not present at today's sensitivity: drift of the device since April, April's reading not what it seemed, or the explicit schedule here placing the gates where April's unrecorded one did not; this flight cannot tell which |
| INCONCLUSIVE | undecided at this sensitivity |

| V2 | reading |
|---|---|
| NO SIGNAL | the odd-Y pattern across all arms holds no real-breaking letter at this sensitivity; the best letter, its size and the unseparated list are printed but not read |
| LETTER FOUND, plausible size | consistent with letter b, and with any letter listed as not separated from it; with none listed, with no other single library letter |
| LETTER FOUND, implausible size | the pattern is b's (or an unseparated letter's), the size is not a calibrated device's: something outside the library that looks like b |
| SIGNAL, NO SINGLE LETTER | a real-breaking signal that no single library letter carries whole: composite, or a letter outside the library close to one inside |
| SIGNAL OUTSIDE THE LIBRARY | a real-breaking signal that is none of these nineteen physical letters |

Three pairings need a sentence of their own. REPLICATED with NO SIGNAL
means the anomaly is there but too localized to stand out across 288
readings: the letter test cannot attribute it, and that is the reading.
A SIGNAL with V1 ABSENT is read with the readout first in mind: the
baseline carries the measured readout matrix, but a readout effect that
changes between the RO circuits and the others would land in V2 and not in
V1. So would a residual excitation at the start, which the RO arm cannot
tell from a readout flip: counted as readout, it leaks into odd-Y readings
that a mixed real start would leave at zero. At q4's 0.83 % today this adds
about 14 to χ²₀; §6 caps P(1|0) at 1 % so that it stays below the carrying
line of 25. A best letter whose fit is flagged unconverged keeps its verdict, and
the flag is reported beside it.

Everything else is read without a verdict. There are 288 odd-Y readings
plus the even ones; at 3σ about one sighting is expected by chance, so a
single reading beyond 3σ is not reported as a finding, and beyond 4σ only
as a sighting.

## 6. Flight conditions

`run_reality_window.py --hardware` refuses to submit unless all of these
hold: the frozen circuit file `simulations/flight/reality_window_circuits.qpy`
has sha256 `e37a378f97070f6077b00a520993c3364f33847df77d9955645ff0b0479d3f62`;
every T1, T2, readout and CZ property of the path is less than 24 hours old
by its own date; on q4, q5, q6 T1 ≥ 100 µs, T2 ≥ 80 µs, readout error
≤ 3 %, P(1|0) ≤ 1 %, both readout asymmetry fields present, and sx
error ≤ 0.1 %; both
CZ errors ≤ 0.8 % (the April baseline had 0.633 % on 4-5); the allocation
reports at least 10 minutes remaining. There is no fallback path or
backend; if a condition fails, the flight waits.

Execution options are set explicitly, where the April runner left them to
the service defaults: job mode, no dynamical decoupling, no gate or
measurement twirling, qubit initialization on, classified measurement,
max_execution_time 600 s.

Cost: the 2026-04-30 job (36 circuits at 4096 shots) billed 40 s of usage
(`job.usage()`, read 2026-10-28). This job carries 196 circuits at 4096
shots, about five and a half times the shot volume, so about 220 s; the
allocation holds 4954 s and runs until 2027-04-20.

The runner writes a stub with the submission order by circuit name, the
snapshot and the hash before it submits, adds the job id as soon as the job
is accepted, then writes the raw counts of both halves to `simulations/flight/results/` (gitignored)
before any analysis; `--recover STUB` rebuilds the counts file from the job
id if the process dies in between. The counts file also carries the job's
own timestamps; a recalibration between the snapshot and the run would
show there, and no rule reads it. After the flight the files are copied
to `data/ibm_reality_window_october2026/` and committed with the record.
