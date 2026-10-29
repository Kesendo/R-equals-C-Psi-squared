# The reality window: a hardware prediction, pre-registered and flown

<!-- STAGE 0 SWEEP RECORD
Stores searched 2026-10-28 for Y-parity, n_Y parity, reality, complex conjugation as a symmetry of the generator, and any prior reading of odd-Y observables on hardware.
docs/ANALYTICAL_FORMULAS.md: F119 (the antilinear triangle, docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md lines 43 to 45: conj(A) carries (−1)^{n_Y}, θ is the mirror D = diag((−1)^{n_Y})); line 6879 notes that σ⁻ and σ⁺ are real. Neither is read as a noise detector.
experiments/: MARRAKESH_THREE_LAYERS.md line 92 states the rule for the 2026-04-26 run (XY+YX keeps a real ρ₀ real, every odd-Y string vanishes, the Y-parity fw.y_parity_panel reads); line 95 reports the predicted zeros reading up to 0.098 for soft on hardware; line 198 lists soft ⟨Y₀Z₂⟩ = +0.098 and ⟨Y₀⟩ = −0.079 among residuals "what carries them we do not know".
data/: ibm_f83_signature_april2026/README.md lists pi2_odd ⟨Y₀Z₂⟩ = +0.067 and ⟨Z₀Y₂⟩ = +0.053 as "what carries those we do not know" and reads pi2_even's ⟨Y₀⟩, ⟨Y₂⟩ only through X^⊗3; ibm_soft_break_april2026/README.md line 78 lists soft ⟨Y₀Z₂⟩ = +0.098 without comment, and its lines 66 to 67 say continuous Lindblad, Trotter and hardware agree in sign on ⟨Y₀Z₂⟩, ⟨I₀Y₂⟩, ⟨X₀Y₂⟩, which for soft cannot hold (the models give exactly zero there); repaired with the flight record. The Kingston 2026-05-05 JSON (job d7sqjpiudops73976960) has no README reading its odd-Y values.
simulations/: f83_hy_field_check.py tested an h_y field, which keeps every imaginary Hamiltonian imaginary and so cannot light an odd-Y reading in three of the four classes; no Z-detuning, ZZ or per-gate model was fitted to these runs before this document; framework/symmetry.py y_parity_panel builds L from Z-dephasing and T1 only and reports Y-parity preserved by construction.
fw.Confirmations: f83_pi2_class_signature_marrakesh (the 2026-04-30 run and its Kingston reproduction) and regime_uniformity_kingston_uniform_quantum (the same Kingston run through the F88b lens) register these runs; marrakesh_transverse_y_field_detection (h_y_eff ≈ 0.05, Z⊗N mirror test, path [48,49,50]); none reads Y-parity.
OpenArcs registry: no arc on the F83 X^⊗3 window or on Y-parity on hardware; line 2591 names Y-parity as a remaining calibration instrument.
docs/CAUGHT_ERRORS.md: the 2026-10-28 Trotter-order entry records the YZ+ZY window through X^⊗3; line 4003 (a detuning does not freeze a phase) is the nearest precedent on detuning readings.
Typed layer: the eight YParity-axis claims in compute/RCPsiSquared.Core/Symmetry (F102, F103, F105, F106, F107, F109, F110, F111; HardCellYInversionPattern.cs line 145 lists them) classify Hamiltonian terms by y_par, MotherSoftYParityOnePurity being the typed home of the premise that XY+YX carries odd y_par; none has dynamics or a hardware reading; no Diagnostics witness reads odd-Y observables.
docs/GLOSSARY.md: Y-parity appears only as a term label.
Post-flight, 2026-10-29, for the reading of section 7: experiments/PRICE_PAIR_HARDWARE_PREDICTION.md (lines 100, 102, 157) measured residual detuning offsets of 8 to 22 kHz on ibm_marrakesh lines [2, 3, 4] and [93, 94, 95] in July, +21.5 kHz of them on logical q2 of [2, 3, 4], which its runner (initial_layout = path) puts on physical q4, and coherent nearest-neighbour ZZ near 4 kHz. This is the adjacency the pre-flight sweep missed: the detuning bound of section 5 lies below an offset the repo had measured on this chip. fw.Confirmations (f83_pi2_class_signature_marrakesh, named above) and its C# mirror record YZ+ZY's ⟨Y₀⟩, ⟨Y₂⟩ near 3σ with the cause unknown and nothing on XY+YX's odd-Y readings; ANALYTICAL_FORMULAS holds nothing on either; the OpenArcs registry has no arc for them.
Post-flight, 2026-10-29, for the cube reading of section 7: docs/THE_ONE_SQUARE.md §9 (the damping row and gate M8: the arrow I → Z under the three turns, the quarter-turn about Z, the shift, the swap, D and R; no quarter-turn about X appears there; compute/MirrorWorld/LetterCube.cs TableMoves names the quarter-turn about X as a uniform move and nothing applies it to the arrow); compute/MirrorWorld.Tests/LetterCubeTests.cs lines 420 to 421 (an arrow leaves I along [J, J†], exactly when J is not normal); compute/MirrorWorld/LetterCube.cs TransferTimesFour, HeightsTimesFour, Arrows (the PTM in Gaussian integers; no method carries a jump by a move or flags an arrow into an odd-n_Y letter), Theodolite.Setting.Rotated (two Clifford axes, X + Z and 1 − iX, no sx†, no angle), TiltedAxes (Hermitian integer axes only), Sextant (copy and reflection verdicts, no letters); simulations/framework/cube_tables.py transfer_times_four (the Python twin; simulations/framework/diagnostics/cube_letters.py holds no PTM); framework/symmetry.py y_parity_panel (Z dephasing and T1 in the logical frame only); docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md line 183 (σ⁻ in a turned frame graded by bit_b: the adjacency joined in section 7, neither document had cited the other); docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md line 91 (conj fixes a real rotation, so a Hadamard keeps the window and an sx does not) and PROOF_F84_AMPLITUDE_DAMPING.md (tilted Hermitian axes only); experiments/THE_CUBE_SITE_BY_SITE.md line 154 lists weights in the Hadamard and quarter-turn frames as open, a different question (the Theodolite's weight search) left as it stands; reflections/, hypotheses/, data/*/README.md: nothing on a jump carried by a gate; fw.Confirmations: nothing on frames or rotations; the typed layer: AntilinearTriangleClaim carries the transport law, MirrorOrderSortingClaim F131, no claim or Diagnostics witness reads a frame-carried jump.
Adjacency checked: the F83 README, MARRAKESH_THREE_LAYERS and the soft_break README describe the same zero and do not cite each other; the h_y rejection and the reality rule are about the same object and do not meet. This document joins them.
Post-flight, 2026-10-29, for the three carriers tested in section 7 and the ones left untested, the same stores searched again. A residual exchange: experiments/F81_VIOLATION_HARDWARE_BRIDGE.md line 125 bounds the idle exchange to a spectator neighbour by (J/Δ)² ≈ 10⁻⁸, and its line 80 notes that an exchange term breaks the per-site Bloch form; the sweep found no model of an exchange during the CZ. Drive crosstalk: the sweep found no model on hardware. CZ leakage: experiments/CORNER_BEAT_HARDWARE_PREDICTION.md lines 2399 to 2404 bound it at 1e-4 to 1e-3 per gate, unbudgeted, in a flight never flown; the sweep found no measurement of it. Quasi-static dephasing: experiments/PRICE_PAIR_HARDWARE_PREDICTION.md lines 54 and 77 read local quasi-static noise on this chip as a contrast effect, rates adding in quadrature, and experiments/IBM_F129_RAMSEY_FRINGE.md lines 220 and 307 carry it as quasi-static site disorder in a T2* model; docs/CAUGHT_ERRORS.md line 3286 keeps quasi-static K as a planning model; none reads an odd-Y observable. That an odd-Y reading is odd in the real-breaking letters is F131 (docs/ANALYTICAL_FORMULAS.md line 6335, the order-sorting law, typed as MirrorOrderSortingClaim in compute/RCPsiSquared.Core/Symmetry): its ODD-response cell, with conjugation as the antiunitary mirror (χ_M = −1, from F119: docs/ANALYTICAL_FORMULAS.md line 5560, docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md line 91, conj∘Ad_U∘conj = Ad_Ū), a real detuning as the direction it fixes, and an odd-Y readout. The sweep found neither applied to a hardware reading or a noise average before; section 7 applies them, and this pair (F131 and the reality window of section 1) is a second adjacency the pre-flight sweep missed. experiments/ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md lines 44 to 46 find a residual population unidentifiable from calibration arms alone, and docs/CAUGHT_ERRORS.md lines 3255 to 3266 record the same: the same object as the residual excitation counted as readout here, and the reason section 7 asks the odd-Y readings rather than the RO arm. The Marrakesh snapshots in data/ibm_calibration_snapshots/ put q4's P(1|0) at 0.10, 0.12 and 0.05 % on 2026-04-25, 04-30 and 05-05 and its T1 at 177, 206 and 102 µs; experiments/F81_VIOLATION_HARDWARE_BRIDGE.md line 58 read physical q4 in July (run 1 on [2, 3, 4]) without reading a defect, and the sweep found no store that reads one. The OpenArcs registry, the glossary and the Diagnostics witnesses returned nothing on these carriers.
-->

This flight asks one question of the IBM hardware, and the letter cube asks
it: where does the reading come from that the circuit forbids? Sections 1
to 6 were committed before the flight (commit 043fa18f); it flew on
2026-10-29, and §7 holds the verdicts and a reading beside them.

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

## 7. The flight

Flown 2026-10-29 on ibm_marrakesh, job `dbhcohht9uhs73c9vpj0`: queued
04:32 UTC, ran 05:12 to 05:16 UTC, billed 215 s. The conditions of §6
held at submission on the backend's numbers: last update 03:02 (+01:00),
T1 and T2 measured 2026-10-28 10:14 (+01:00) and the readout 09:37, under
24 h old at the run. The properties the backend reports for 05:14 UTC are
the snapshot's: no recalibration the backend reports fell between the
snapshot and the run.
The record is in `data/ibm_reality_window_october2026/` (the counts, the
submission stub, the properties at the run), the runner's printout in
[`simulations/results/reality_window/flight.txt`](../simulations/results/reality_window/flight.txt).

**The verdicts**, by the rules of §5 (the runner's printout agrees):

- **P: PASS.** Arm A reads ⟨X₀Z₂⟩ = −0.664 and ⟨Z₀X₂⟩ = −0.485, apart by
  0.179. The day's baseline, with the readout the RO arm measured, predicts
  −0.747 and −0.516; §4, with the snapshot's readout, −0.836 and −0.595.
  April read −0.849 and −0.566 on the same path, close to §4's values.
- **V0: PASS.** C = +0.105 against C_pred = +0.121 on the day's baseline,
  band (0.073, 0.147); §4's +0.137 assumed the snapshot's readout. The sign
  chain from the gates to the verdicts holds.
- **V1: REPLICATED.** s_A = +0.083 at σ = 0.0125 (κ = 1.27), 6.7σ, with
  April's sign. The raw reading is +0.072; the baseline it is taken
  against, −0.012, is the day's readout leaking even readings into S,
  twice §4's 0.006 because the readout was worse (below). For the same reason the readout matrix's own
  noise, which §5 leaves out of V1 at about 0.0005, is 0.002 on the day;
  with it s_A stands at 6.6σ. The half-sum of K± gives +0.100. The earlier
  three raw readings of this observable were +0.098, +0.067 and +0.066.
  Read beside §2: today's s is more than twice the
  +0.035 the April letter fit reached, and that only with letters beyond
  their bounds.
- **V2: SIGNAL, NO SINGLE LETTER.** χ²₀ = 4167 against its 99 % point
  346.8 at rank 288. The best letter, detuning on q6 at +0.175 rad/µs
  (beyond its bound of 0.1), carries (χ²_b = 2015, a drop of 2152) and does
  not suffice (2015 against the 99 % point 345.7 at rank − 1); the next
  letter sits 58 above it, so none is unseparated, and no fit is flagged
  unconverged. A real-breaking signal is present that no single library
  letter accounts for whole. One guard of this reading was not met on the
  day: §6 caps P(1|0) at 1 % so that a residual excitation, counted as
  readout, stays below the carrying line, and the RO arm put q4 at 3.5 %
  (below). §5 estimated such a term at 14 for 0.83 %; computed through the
  replay it stays at a few units at most (below, under *Three carriers
  tested*).

**The readout of the day.** The RO arm measured q4 far off the snapshot:
P(1|0) = 3.5 % and P(0|1) = 10.1 %, against 0.83 % and 3.1 %; q6 read
0.7 % and 2.0 % against 0 and 0.83 %, q5 near its snapshot. The flip rates
of q4 and q6 also differ between the job's two interleaved halves
(Σz² = 16.5 and 16.8 on 8 dof each; pooled over the three qubits 35.1
on 24 dof, p ≈ 0.07): within what chance allows, and not read. The baseline
carries the measured matrix, so the verdicts read against the day's
readout, but by the RO arm's own number q4's P(1|0) stood above the 1 %
that §6 caps on the backend's number. That cap rested on §5's estimate
that a residual excitation at the start, counted as readout, adds about 14
to χ²₀ at 0.83 %; at four times the flip rate the backend's numbers do not
secure the guard's purpose on the day (its own init_error for q4 reads
0.84 %, dated 2026-10-27). Computed through the replay, such a term stays
far below the carrying line, so the purpose holds although §5's estimate
does not (below, under *Three carriers tested*). How much of q4's
departure is such an excitation is not known: the odd-Y readings cannot
tell, and the even readings were not asked. The readout does
not carry the bulk of the signal: A and D3
share their readout and S baseline (−0.012) and read S raw +0.072 and
−0.017, and the largest asymmetry the RO arm measured, q4's 0.066 (q6's
is 0.013), moves a single-bit reading by at most that, where D6 carries
residuals of ±0.2.

**Read beside the verdicts** (not pre-registered; the script is
[`simulations/reality_window_post_flight.py`](../simulations/reality_window_post_flight.py),
its output
[`simulations/results/reality_window/post_flight.txt`](../simulations/results/reality_window/post_flight.txt);
its fits carry no bootstrap term, so its total is 4384 rather than the
verdict's 4167).

The excess is not in ⟨Y₀Z₂⟩ alone. Each arm's 36 odd-Y readings, taken by
themselves against the baseline (columns in order of χ²; s is ⟨Y₀Z₂⟩
minus its baseline, per arm):

| arm | D0 | D1 | D2 | D3 | A | D4 | R | D6 |
|---|---|---|---|---|---|---|---|---|
| χ² over 36 | 68 | 73 | 163 | 205 | 272 | 677 | 934 | 1992 |
| s | −0.023 | −0.017 | +0.003 | −0.005 | +0.083 | −0.091 | +0.107 | −0.080 |

Every arm is above its 99 % point of 58.6. In D0 the window cannot tell a
real-breaking letter from a population error: a D0 circuit carries at most
one sx per qubit, and the circuit that reads q2's odd ⟨Y₂⟩ is physically
the one that reads its even ⟨Z₂⟩ (the rz around the pulse act on |0⟩ and
before a Z measurement). Both read the same offset, −0.03 to −0.05 against
the baseline: the population after one pulse on q6 is off, or its readout
differs from the RO arm's, and the same physical error reads as real in
one circuit's logical frame and as real-breaking in the other's. Fitted
to D0 alone, sx amplitude errors on the three qubits bring it from 67.8 to
35.0, below its 99 % point of 54.8, at +2.3 % on q6, beyond its bound of 1 %; fitted over all arms beside the
detunings and static ZZ, sx amplitude on q6 comes out at +0.1 %.

Among the fresh compiles the excess grows with the Trotter depth, up to
residuals of ±0.2 at six steps (⟨Y₂⟩ +0.22 and ⟨X₀Y₁⟩ −0.20 in the XY
basis of D6). That growth does not single out a carrier. The real sector
loses contrast with depth on the same day (the even readings against the
baseline: slope 1.00 at D0, 0.98, 0.92, 0.83, 0.79, 0.78 at D6, A and R
0.85 and 0.84), any error repeated per step grows the same way, and depth
is not the only order: R, A's own gates with physical 4 and 6 exchanged,
is above D4, so where the gates run sets the reading too. So does the
compile: A and D3 are one logical circuit at the same depth (1.50 and 1.46 µs long); their 36 odd-Y residuals correlate
at 0.29, and s reads +0.083 in one and −0.005 in the other.

No set of library letters fitted here suffices. On the post-flight scale
the three detunings leave 1506, with static ZZ on both pairs 1381 (against
341), with sx amplitude added 938; all nineteen, fitted nonlinearly from
the five-letter optimum, leave 707 against 326. These
are local fits within this replay, whose even-sector contrast is itself
above the day's; they do not exclude a better fit elsewhere. Fitted arm by
arm (D0 left out: §3 reads it as blind to them), the five letters
leave χ² 70 to 182 against 52, at sizes that change from arm to arm. The
sizes come from fits that do not suffice and are not read as the device's.
Apart from those fits, §5's detuning bound of 0.1 rad/µs (16 kHz) is
itself below an offset the repo has measured on this chip: +21.5 kHz on
physical q4 in July
([PRICE_PAIR_HARDWARE_PREDICTION](PRICE_PAIR_HARDWARE_PREDICTION.md)), on
the qubit where today's joint fits put no detuning. The sentence is about
the bound, not support for the fitted sizes.

The job's two halves, interleaved by the shuffle, carry the same excess:
their odd-Y z-scores correlate at 0.886, and each half alone gives a
diagonal χ² of about 3200. The bulk of the carrier is not a drift
during the job.

**Three carriers tested**, through a copy of the replay with the start
state and an exchange letter patched in, on the post-flight scale
(χ²₀ = 4384). With no letter the copy equals the frozen replay exactly,
a gate on the copy's baseline path and not on the new terms; the exchange
generator's convention is gated apart, and the frame bookkeeping and the
excitation's start and readout steps were checked by reading and a test
run, not gated. The script is
[`simulations/reality_window_carriers.py`](../simulations/reality_window_carriers.py),
its output
[`simulations/results/reality_window/carriers.txt`](../simulations/results/reality_window/carriers.txt).

- *A residual excitation at the start.* The RO arm's x gates act on the
  same start, so the matrix it measures is the true readout after the
  start's flips, A_RO = A_true·M(p), and the model starts q4 in
  diag(1 − p, p) and reads out with A_true = A_RO·M(p)⁻¹. The odd-Y
  expectations are blind to it, and the measured parities move only
  through the readout: by at most 0.0034 at p = 3.4 % (q4's P(1|0) in the
  RO arm's all-zero preparation; averaged over the preparations it is the
  3.5 % above). The reason, read by hand: a start in |1⟩ instead of |0⟩
  is a real state, so the window stays closed; and a symmetric flip leaves
  each bit's readout asymmetry P(0|1) − P(1|0) exactly unchanged and only
  scales its contrast by 1 − 2p (the same division by 1 − 2p₀ that
  [ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION](ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md)
  meets), so counting it as readout leaks into the odd readings only a
  term of order p times that asymmetry times a change of even readings.
  On readings that held nothing else, the excitation would add
  its noncentrality d·W·d to χ²₀: 0.04 at p = 0.83 % and 0.70 at 3.4 %. On
  the measured readings χ² moves by +0.6 and +2.5. §5's estimate of 14 at
  0.83 %, the premise of the cap in §6, is not what this replay gives, and
  how it was obtained is not on record. Fitted on all three qubits, each p
  bounded by the flip rate on its qubit, p goes to zero on q4 and to its
  bound on q5 and q6, and χ² moves from 4384 to 4382. The bound is
  approximate: the RO matrix has cells that read zero, so A_true holds
  small negative entries at every p > 0, down to −0.006 at q4's bound.
  Whether q4's departure was an excitation the odd-Y readings cannot tell;
  the even readings, which an excitation shifts at first order, were not
  asked.
- *Quasi-static detuning.* In the logical frame an odd-Y reading is odd
  in the real-breaking letters taken together: complex conjugation keeps
  the generator −i[H,·] of the imaginary H and every real jump, and turns
  the evolution under a detuning δ into the one under −δ (F119:
  conj∘Ad_U∘conj = Ad_Ū), while the odd-Y Pauli string changes sign. This
  is the ODD-response cell of
  [F131](../docs/ANALYTICAL_FORMULAS.md#f131), the mirror's order-sorting
  law: conjugation the antiunitary mirror, a real detuning the direction it
  fixes, an odd-Y readout of parity −1. A spread of detunings symmetric
  about zero cancels there. In the compiled circuit F131's premise fails
  and the spread need not cancel exactly, since conjugation also turns
  each sx into its inverse and each rz(θ) into rz(−θ). Averaged over ±σ on each of
  the three qubits through the replay, the largest shift of an odd-Y
  reading is 0.0002 at σ = 0.03 rad/µs, 0.0023 at 0.1, 0.009 at 0.2 and
  0.019 at 0.3, and χ² rises at every σ (4386, 4407, 4483, 4628). About a
  nonzero mean the spread mixes with the mean at second order, so that
  case was computed on its own: about the three-detuning optimum
  (χ² 1505.5) the spread moves χ² by less than 4 up to σ = 0.2 and worsens
  it beyond (1562 at 0.3). Within this replay, a spread independent on
  each qubit, two-point, up to σ = 0.3 rad/µs, about zero or about those
  detunings, does not carry the excess; a spread common to the qubits was
  not tried.
- *A residual exchange at the CZ*, U = exp(−iεG) at every cz with
  G = e^{iφ}σ⁺_aσ⁻_b + h.c., a the pair's lower physical qubit: real at
  φ = 0, (XX+YY)/2, and imaginary at φ = π/2, (XY−YX)/2. In the logical
  frame the first breaks the window and the second keeps it; inside the
  compiled circuit the gates around each cz decide which part breaks it.
  With the phase fixed in the hardware frame (the virtual rz accumulated
  before each cz subtracted), the exchange alone leaves χ² = 4115, and
  beside the three detunings and static ZZ on both pairs 1331, against the
  1381 of those five alone: about fifty for four parameters, something and
  not enough. With the phase fixed in the circuit's own frame instead the
  two fits leave 2743 and 1102. Such a phase would follow every virtual
  rz, as an error tied to the drive's frame would; a static coupling
  carries neither fixed phase (below), so the better fit is a reading, not
  support for that frame. Each fit is the best of five starts (the
  perturbed ones in the ε cos φ₀ components only) and, beside the five
  letters, a sixth at their optimum; the best converged (least-squares
  status 2), and they are local fits. The fitted sizes, ε from 0.007 to 0.025,
  keep ε² below the day's CZ errors of 0.22 %; no fit comes near its 99 %
  point of 342 or 337. Left out are two phases: the (ω_a − ω_b)t that a
  static coupling carries in the rotating frame, which moves with the
  absolute time of each cz, and any virtual phase the calibrated cz carries
  inside itself, which would move the hardware frame at every cz. The
  backend no longer publishes its qubit frequencies (none in its
  properties, target or configuration on 2026-10-29), so an exchange whose
  phase moves that way is not tested.

**Read on the cube** (one site, by hand, after the fits; the second route is
[`simulations/reality_window_cube_frames.py`](../simulations/reality_window_cube_frames.py),
exact in sympy, 22 rows with two anti-vacuity rows, which must print "cube
frame gate: ALL GREEN"; the house instruments carry the same moduli in
Gaussian integers, the Theodolite's `Setting.Rotated("Y")`, which turns a
site by 1 − iX, through the LetterCube's `TransferTimesFour`, pinned in
`TheodoliteTests.In_The_Quarter_Turn_Frame_The_Damping_Arrow_Lands_On_Y`,
with the arrow's sign reversed there because `Rotated` carries in the
opposite direction). The zero of §1 is a statement about letters in the
logical frame, and a physical process acts in the frame the gates before it
have set, so the circuit carries every physical letter into the logical
frame by the Bloch rotation of those gates, Ad_U† ∘ P ∘ Ad_U for the gate
product U: an rz(θ) turns about Z by θ, an sx is the quarter-turn about X,
a cz keeps every one-site letter. On one site the real letters are I, X
and Z and Y is the odd one. A one-site process is a height, an arrow or a
field ([THE_ONE_SQUARE](../docs/THE_ONE_SQUARE.md) §9), and the three fare
differently when carried into the logical frame.

- A Pauli dephasing is a height along its axis, and a height moves nothing.
  The carried axis stays a Pauli axis as long as no rz of a generic angle
  precedes an sx in the frame (Z under rz(θ) stays Z, under sx it becomes
  Y), so physical Z dephasing keeps the window closed in those frames at
  any strength. After an rz(θ) followed by an sx the axis lies tilted in
  the XY plane, sin θ·X + cos θ·Y, and a dephasing along it carries the
  arrow X ↔ Y with coefficient γ·sin 2θ, γ the rate of D[Z], whose heights
  on X and Y are 2γ cos²θ and 2γ sin²θ: the arrow is the geometric mean of
  the two heights (a dephasing along a general axis n carries X ↔ Y at
  2γn_x·n_y and Z ↔ Y at 2γn_y·n_z, and never moves I). So dephasing does
  break the window in a generic frame, and it pays for it with the same
  rate that takes contrast from X and Y: its odd-Y weight is at most the
  time spent in such frames times |sin 2θ| over 2T_φ, T_φ the pure
  dephasing time (1/T_φ = 1/T₂ − 1/2T₁).
- The damping σ⁻ is two lights of γ/2 and one arrow, I → Z, of the full
  rate γ, the arrow pointing along [J, J†] (`LetterCubeTests`). Carried by
  sx the arrow is I → Y: the trace, weight one in every state, feeds the
  odd letter at the rate 1/T₁ for as long as the frame holds Z on the Y
  axis. Carried by a Hadamard (rz·sx·rz) it is I → X, real. This is §4's
  floor of 0.0023 read as a why: the source is the trace, the rate is
  1/T₁, and the floor is bounded by the frame time over T₁ (with the tilted
  dephasing's share beside it).
- A field is an arrow between the two letters perpendicular to its carried
  axis, at twice its coefficient on the Pauli letter (−i[hZ, X] = 2hY), and
  pays no rate. A Z letter (a detuning, a phase, one leg of a ZZ) is the
  arrow X ↔ Y in the logical and rz frames, real-breaking; inside an sx
  frame it has become a Y field, Z ↔ X, real-preserving; after rz(θ) then
  sx it lies in the XY plane and its X component sin θ is the
  real-breaking part (Z ↔ Y). An X letter (a drive: an amplitude error, a
  neighbour's pulse reaching this qubit) in a frame rz(φ) is real-breaking
  through cos φ (Z ↔ Y) and real-preserving through sin φ.

Read together: every incoherent one-site process, in every frame, pays for
its odd-Y weight with the rate that takes contrast from the real sector.
At the snapshot's T₁ and T₂ that budget is the circuit's duration over
them: 1.50 µs for A, R and D3, 2.82 µs for D6, against 146 to 226 µs, so
0.01 to 0.02 at most, and the replay's 0.0023 is the value for these
circuits. The measured +0.083 in A and the ±0.2 residuals of D6 stand
eight and ten times above those budgets. The rates need not be the
snapshot's: against the baseline the real sector kept a contrast of 0.78
to 0.85 in A, R, D3, D4 and D6 (the even-reading slopes under *Among the
fresh compiles*; 0.92 to 0.98 in D1 and D2, 1.00 in D0, which has one
pulse and nothing to lose), which is what an incoherent process thirty to
forty-five times the idle dephasing rate, acting while the qubit is
driven, would take, and the drive is where the frame is tilted. So the
cube's two incoherent candidates were asked of the data
([`simulations/reality_window_pulse_dephasing.py`](../simulations/reality_window_pulse_dephasing.py),
its output
[`pulse_dephasing.txt`](../simulations/results/reality_window/pulse_dephasing.txt);
the frozen replay patched with an extra process on the driven qubit during
its own pulses, gated to equal the frozen replay exactly with the extra
rate at zero, which checks the baseline path and not the new terms). An
extra Z dephasing, coherence rate G, does not carry: at G = 0.2 /µs the
even slopes of every arm but D0 come to 1.00 ± 0.05, so the contrast loss
is accounted for, while no odd-Y model value moves by more than 0.012
against the baseline (the median shift is 0.002) and χ² rises from 4384 to
4550; at G = 0.05, 0.5 and 1.0 it rises as well (4424, 4789, 5084); the
global fit runs to −0.32 /µs and the per-qubit fit puts two of three rates
below zero, lowering χ² to 4226 and 3715 only through rates that are no
channel, and the same extra dephasing at all times, idle included, does
the same. An extra amplitude damping of the driven qubit during its
pulses, rate D, is the stronger of the two, as the cube says (the trace feeds its arrow): at D = 0.2 /µs,
where the even slopes again come to 1.00 ± 0.05, the largest odd-Y shift is
0.034 (median 0.006), three times the dephasing's, and χ² falls from 4384
to 4176, a drop that carries by §5's rule and does not suffice (the 99 %
point is 346); the global fit sits at 0.18 /µs with χ² 4173, the per-qubit
fit at 0.12, 0.06 and 0.67 /µs leaves 3539 with the even slopes at
1.00 ± 0.03, and in arm A the model's ⟨Y₀Z₂⟩ moves away from the measured
+0.072, from −0.011 to −0.018 at D = 0.2. So the damping leaves a part and
not the bulk, and its sizes are not read. The incoherent bound is loose,
then, by the geometry of the frames: the odd-Y weight an incoherent process leaves is a
small fraction of the contrast it takes, because the frames with sin 2θ
near one are a small part of a qubit's pulse time (0.32 to 0.68 µs of A's
1.50 µs is pulse) and the arrows of different pulses carry different
signs; which of the two is not separated here. Within this replay the
carrier is coherent (a field, a two-site phase, the phase a leaked level
brings back) or lies outside the qubit's letters (readout, leakage). What
the cube says a coherent carrier would leave: a letter anchored in the
drive rather than in Z reads in the logical frame through the accumulated
virtual phase of its qubit, which differs between two compiles of one
logical circuit and between the basis circuits of one arm, so one physical
size shows as different logical sizes from arm to arm, the shape of the
arm-by-arm sizes of the five-letter fits above; crosstalk from a
neighbour's drive, untested above, is such a letter, and the exchange fit
anchored in the circuit's frame beating the one anchored in the hardware's
(2743 against 4115) is a reading in the same direction, not support for
either carrier. Adjacency: the T1 proof's second reading
([PROOF_F82](../docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md))
holds σ⁻ in a turned frame, "relaxation toward a y eigenstate, (Z + iX)/2,
has one" (an odd part), graded by Π² = (−1)^bit_b, where bit_b is the
cube's k_X mod 2; the window grades by conjugation, n_Y = (k_Z + k_X −
k_Y)/2, the half-sum no coordinate parity sees. The two gradings read one
object two ways: Π² sees the arrow I → Z of the logical frame as odd, the
window sees it as real and the carried arrow I → Y of the quarter-turn
frame as odd. THE_ONE_SQUARE §9 applies its gated moves to the arrow, the
three turns, the quarter-turn about Z, the shift, the swap, D and R, and
not the quarter-turn about X; that turn is what the compiled circuit
applies.

What carries it we do not know. Not sufficient, within this replay and
these fits: an extra dephasing of the driven qubit during its pulses, or
at all times (the even sector fixes it near 0.2 /µs during the drive, where
it moves no odd-Y reading by more than 0.012); an extra amplitude damping of
the driven qubit during its pulses (it carries, a drop of 208 at the rate
the even sector fixes, and leaves 4176, per qubit 3539); each library
letter alone, the detunings with static ZZ,
those with sx amplitude, the nineteen jointly from one start, an exchange
at the CZ with a fixed phase, alone or beside five letters, a spread of
detunings independent on each qubit, about zero or about the
three-detuning optimum; the odd-Y expectations are blind to a residual
excitation at the start. Ruled out: a
recalibration the backend reports between snapshot and run, a drift
during the job as the bulk, a readout of the size the RO arm measured (or one shared
by A and D3) as the bulk.
Not tested: a readout that differs between the RO circuits and the others as
a contributor,
leakage at the CZ, an exchange whose phase moves with the time of each cz
or with a virtual phase inside the calibrated cz, other coherent two-qubit
errors that are not diagonal, a spread of detunings common to the qubits,
excited or leaked spectators outside the
path, crosstalk from the drive of one qubit onto its neighbour, a defect
on q4 other than an excitation at the start (its readout was far off its
calibration), and the baseline's own
readout leak, computed from modelled even readings that are stronger than
the measured ones. Whether the April
readings and today's share a carrier is not known. April's nineteen
letters together sufficed only at sizes far beyond their bounds (§2), on
30 readings; today no fit tried here suffices, on 288. Neither day's
fitted sizes are a measurement to compare.
