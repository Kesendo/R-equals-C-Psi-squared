# The X^⊗3 window under a turned frame: a hardware prediction, pre-registered

<!-- STAGE 0 SWEEP RECORD
Stores searched 2026-11-13 for the X^⊗3 window of YZ+ZY, sign twins, an echo or delay before the last pulse, twirling, a lit string decaying at the price, a conserved Z parity on hardware, and the dephasing rates prior flights inferred.
data/ibm_f83_signature_april2026/README.md, section "The X^⊗3 window": the break measured on Marrakesh [4, 5, 6] 2026-04-30 and on Kingston [43, 56, 63] 2026-05-05 (the Kingston JSON in data/ibm_soft_break_april2026/, whose README documents only the 2026-04-26 Marrakesh run), what it is not (the Trotter circuit, stochastic Pauli noise between the factors), the letter fingerprint; simulations/f83_xn_window.py its post-hoc fit over 18 readings of two classes (94 with ideal gates, 70 with an end detuning, 3.6 with all nineteen letters, seven at their bounds; that 70 is a different sum from the ten-reading 70 of §2 below).
experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md and data/ibm_reality_window_october2026/: §2 fitted the nineteen letters to the forbidden readings of all four April classes at once, XY+YX and YZ+ZY from the same job among them (30 readings), and reached a small sum only with a middle-qubit detuning twenty times its bound; flown 2026-10-29 on the same Marrakesh path, XY+YX only; §5 the letter bounds; §6 the execution options set explicitly; §7 the odd-Y excess no library letter carries, the three carriers tested, the price-pair detuning of +21.5 kHz on physical q4 above §5's bound.
experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md: the frame-turn law (§1), the arms, the forest, the flight conditions; not flown; its freezer and runner (simulations/flight/drive_frame_freeze.py, run_drive_frame.py) build both pages' arms. This page rides in its Marrakesh and Kingston jobs, and that page carries the riding arms in its §3, §4 and §6 and shares its §5 rules (the drift weighed block by block, the simulated flights drawn from the replay, the 99.9 % points for TURNS, verdicts that need no fitted letter) with this one.
experiments/Z_N_PARTNERSHIP.md lines 85 to 104: drive bleed-through and ZX or ZY crosstalk named as breakers of the Z⊗N mirror, the turning class of this page; Z⊗N carries |+−+⟩ to |−+−⟩ and YZ+ZY to −(YZ+ZY), the algebra of the sign twin of §3. experiments/MARRAKESH_THREE_LAYERS.md line 91: truly's X^⊗3 zero set on [48, 49, 50].
data/ibm_zn_mirror_april2026/README.md and fw.Confirmations marrakesh_transverse_y_field_detection: the |+−+⟩ and |−+−⟩ twin flown 2026-04-29 under Heisenberg XX+YY+ZZ on [48, 49, 50], read as h_y_eff ≈ 0.05.
experiments/PRICE_PAIR_HARDWARE_PREDICTION.md and data/ibm_price_pair_july2026/README.md: block C read ⟨XXX⟩ of a GHZ state on Marrakesh in July, its coherence outliving the product of its singles (ln E_GHZ − Σ ln E_j = +1.6 to +1.9 at 65 µs), carried by coherent nearest-neighbour ZZ that XXX commutes with; a colouring reading that experiments/THE_PALINDROME_AS_A_COLOURING.md does not cite.
fw.Confirmations: f83_pi2_class_signature_marrakesh records the window with the cause unknown; the October flight is in neither registry. The OpenArcs registry: no arc on the window, twins, echoes or colouring hardware; its line 3128 notes that a full Pauli twirl kills a flux. docs/CAUGHT_ERRORS.md: the Trotter-order entry (γ_Z_eff 0.03 and 0.10), a Y-odd sign flip from vectorisation, a channel-both-ways design never flown because its discriminating arm was untested. docs/ANALYTICAL_FORMULAS.md: nothing on the window or on hardware for F173. The glossary, the Claim graph in compute/RCPsiSquared.Core and the Diagnostics witnesses: nothing on these readings.
Adjacency checked: the price-pair ⟨XXX⟩ and the colouring page are one object and do not cite each other (left for their own change); Z_N_PARTNERSHIP's mirror and the sign twin are one map (recorded here); the drive-frame page and the X^⊗3 window ask one question of two Hamiltonians (joined here).
-->

The chip that breaks the reality window on XY+YX also breaks a second
selection rule on YZ+ZY. Fits of error letters have not said what opens
either: the October library failed letter after letter, and every carrier
number since has been a guess. This page asks a question that needs no
letter. The April YZ+ZY circuits, gate for gate, are read with one qubit's
drive frame turned, which is the identity on the logical circuit
([the driver's frame](DRIVE_FRAME_HARDWARE_PREDICTION.md) §1), in the same
jobs as that page's XY+YX arms, beside a direct measurement of each
qubit's idle detuning. It is written before the flight; §7 will hold the
flight.

## 1. The zero

YZ+ZY on the open chain of three commutes with X^⊗3 term by term (YZ and
ZY each carry two letters that anticommute with X), and so does every
Trotter factor. Conjugation by X^⊗3 commutes with Z dephasing, which flips
X^⊗3 itself (it is lit) but carries the generator into itself, and
|+, −, +⟩ is invariant under it. So every reading whose Pauli string
anticommutes with X^⊗3 vanishes, in the continuous evolution and in the
Trotter circuit
([the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md),
"Against each term"): on the letter cube X^⊗3 is a colouring of YZ+ZY. Each
term of YZ+ZY carries one Y, so H is imaginary, the start is real, and the
reality window of [the October flight](REALITY_WINDOW_HARDWARE_PREDICTION.md)
§1 closes the odd-Y readings as well. With q0 and q2 read in any of nine
bases and q1 read in X, the forbidden readings are the parities whose
string carries an odd number of Y and Z letters on q0 and q2, or an odd
number of Y: 36 of the 63 parities over the nine bases. They fall into
three kinds, read on the letters of q0 and q2 (q1's X keeps both rules):
one Z and no Y (⟨Z₀⟩, ⟨Z₂⟩, ⟨Z₀X₂⟩, ⟨X₀Z₂⟩, with or without X₁) breaks
X^⊗3 alone; one Y and no Z (⟨Y₀⟩, ⟨Y₂⟩, ⟨Y₀X₂⟩, ⟨X₀Y₂⟩) breaks both; one Y
and one Z (⟨Y₀Z₂⟩, ⟨Z₀Y₂⟩) breaks reality alone.

More is true of this start: YZ+ZY annihilates it. On one bond
(Y⊗Z + Z⊗Y)|+, −⟩ = −i|−, +⟩ + i|−, +⟩ = 0, and so on the other; the
April circuits take the bonds one after the other, and on one bond YZ and
ZY commute, so every Trotter factor leaves |+, −, +⟩ where it is: in the logical circuit the
state is the start again at every boundary between factors. The compiled
circuit merges the factors and carries the state through frames the
single-qubit gates turn, away from the start between its cz gates, and
that is where errors act. Everything a B circuit shows beyond its
preparation and readout is the hardware's. It also means that ⟨X₀⟩,
⟨X₂⟩ and ⟨X₁⟩ stay near ±1, so a forbidden reading with an X partner nearly repeats the one
without it (⟨Z₀X₂⟩ ≈ ⟨Z₀⟩), and the difference of the two lives on rare
outcomes, a direction of the covariance whose true size any incoherent
error the replay leaves out would move. The statistics of §5 therefore
read the 14 forbidden readings that carry no X letter: ⟨Z₀⟩, ⟨Z₂⟩, ⟨Y₀⟩,
⟨Y₂⟩ in each basis that reads them alone, and ⟨Y₀Z₂⟩, ⟨Z₀Y₂⟩; all three
kinds stand among them.

What can open it, in the words of the driver's frame. On the logical
circuit: a coherent Z rotation on any of the three qubits, a turn about Y,
damping (no Pauli channel); not a Pauli channel between the Trotter
factors, nor a preparation error that leaves a real, incoherent mixture
of eigenstates of X^⊗3. Inside the
compiled circuit more: a Pauli error between native gates passes the
non-Clifford rz of a factor and leaves as something else (the replay puts
the October snapshot's T₂ dephasing there at 3.4 × 10⁻⁴), a Z letter of the gates (a
phase after a cz, the replay's detuning acting through every gate) opens
it, and so does an own-frame amplitude error, by an amount that depends on
the qubit (§4). A tilt common to every sx and x pulse on one qubit does
not: it is a turn of that qubit's drive frame and replays to machine zero.
A tilt of the sx pulses against the x pulses does open it, and is
own-frame too (§4).

## 2. What is already measured

Pooled over the three bases that read each letter (the README's estimator,
binomial σ about 0.009 at 4096 shots per basis):

| run | path | ⟨Y₀⟩ | ⟨Z₀⟩ | ⟨Y₂⟩ | ⟨Z₂⟩ |
|---|---|---|---|---|---|
| Marrakesh 2026-04-30, job d7pol1e7g7gs73cf7j90 | [4, 5, 6] | +0.039 | +0.023 | +0.042 | +0.057 |
| Kingston 2026-05-05, job d7sqjpiudops73976960 | [43, 56, 63] | +0.022 | +0.061 | +0.048 | +0.009 |

truly (XX+YY), flown beside it, reads nothing beyond 2σ on Kingston and
⟨Z₀⟩ at 3.8σ, ⟨Z₂⟩ at 3.1σ on Marrakesh (binomial σ; the README's
scatter between bases takes ⟨Z₀⟩ to near 3σ and ⟨Z₂⟩ below it)
([`f83_xn_window.py`](../simulations/f83_xn_window.py) →
[`f83_xn_window.txt`](../simulations/results/f83_xn_window.txt)). The
large ⟨Z⟩ sits on q2 on one chip and on q0 on the other.

The April and the May job each flew XY+YX and YZ+ZY together. The October page's §2
fitted its nineteen letters to the forbidden readings of all four April
classes at once (30 readings) and reached a small sum only with a
middle-qubit detuning twenty times its bound; fitted to YZ+ZY and truly
alone, the same letters reach 3.6 with seven of them at their bounds and
more letters than readings. The letter sets fitted on
the October XY+YX flight do not carry April's YZ+ZY window either:
replayed on April's flown YZ+ZY circuits with April's T₁, T₂ and readout,
over the ten readings the class forbids on (q0, q2) at April's binomial σ,
the sum of squared residuals reads 70 with no letter, 736 under October's
five Z letters, 87 under its crosstalk alone and 406 under all three
families (as late as possible; 71, 885, 99, 479 as soon as possible;
[`x3_window_power.py`](../simulations/x3_window_power.py) part 0). A
further fit would be a further guess. The frame turn sorts the carriers
into two classes without naming one, which no fit has done.

## 3. The arms

| arm | circuit |
|---|---|
| B5(φ), φ ∈ {0, π/2, π, 3π/2} | the 2026-04-30 flown YZ+ZY circuits (class pi2_even_nontruly), nine bases on (q0, q2), gate for gate, then q1 read in X; rz(φ) as the first and rz(−φ) as the last operation on the middle qubit b (logical q1), the last after b's basis change |
| B4(π), B6(π) | the same turn by π on a (logical q0) and on c (logical q2) |
| Q | one path qubit at a time, the other two left in \|0⟩: sx, an idle of τ ∈ {0, 1, 2, 4} µs, read in X or in Y: 24 circuits, built at the gate level (sx, the delay, then the basis change in rz and sx gates), not transpiled, so that the τ = 0 circuits carry the same pulses as the others |

B is built as the October freezer builds arm A
([`reality_window_freeze.py`](../simulations/flight/reality_window_freeze.py),
`from_april`), with the X basis change on q1 in place of the Y one, and
scheduled as late as possible with explicit delays. The turn is the law of
the driver's frame §1, which holds here unchanged: the leading rz acts on
|0⟩, the trailing one stands after b's basis change and before its
measurement. A Z letter and an own-frame letter are blind to it; a foreign
transverse letter coupling to the turned qubit, a neighbour's pulse
reaching it or its pulse reaching a neighbour, moves with it, to first
order as a·cos φ + b·sin φ. Two more carriers turn: an image term of the
turned qubit's own drive, at 2φ, and leakage or a residual swap at a cz.
The turn by π on an end qubit reads only the cosine part, so a foreign
letter in quadrature there reads as blind. Q measures each qubit's idle
precession against its own frame with its neighbours in |0⟩, where the
device calibrates each frequency; an always-on coupling ζ|11⟩⟨11| does not
shift it there, while with a neighbour on the equator it would, by ζ/2.

A sign twin from |−, +, −⟩ (the Z⊗N mirror of
[Z_N_PARTNERSHIP](Z_N_PARTNERSHIP.md), which carries YZ+ZY to its
negative) and a twirl of the cz gates are not flown. Replayed through the
letter library on a fresh compile (level 3, seed 11, against
FakeMarrakesh) under April's T₁ and T₂, the twin changes the sign of some
single-site readings and keeps others: under a phase of 0.05 after every
cz ⟨Z₂⟩
goes from −0.341 to +0.109 while ⟨Z₀⟩ goes from −0.155 to −0.136, so a
sign
comparison does not sort coherent letters from damping and readout. A
twirl suppresses single-qubit letters as well as cz letters
([`x3_window_arms_scout.py`](../simulations/x3_window_arms_scout.py) →
[`x3_window_arms_scout.txt`](../simulations/results/x3_window/x3_window_arms_scout.txt)).
An idle before the last pulse is asked directly by Q.

The B arms and Q ride in the drive-frame jobs on ibm_marrakesh and
ibm_kingston, each circuit twice at 4096 shots, shuffled with the
drive-frame circuits as one permutation by seed 20261113, the RO arm
shared and the forest read in every circuit with the drive-frame page's
count-key layout. On Kingston the path is [43, 56, 63] whenever it passes
the guards (the drive-frame page's rule), May's path; if it fails,
that page's rule picks the path, B and Q fly there, and the record says
the May path was not repeated. On both chips the circuits are April's
Marrakesh gates, not May's Kingston compile, as the drive-frame arms are.
Per chip 6 × 9 + 24 = 78 circuits and 156 PUBs; with the drive-frame's 124
that is 280 PUBs in one job, about 308 s at the 1.1 s per PUB both earlier
jobs billed. Fez flies the drive-frame job alone. The flight is about
750 s, under the 15 minutes Tom set for it, of the 4739 s the allocation
reported on 2026-11-13 (usage_remaining_seconds; the remainder after the
October flight, no job having run since).

## 4. Predictions from below

[`x3_window_power.py`](../simulations/x3_window_power.py)
(output [`x3_window_power.txt`](../simulations/results/x3_window/x3_window_power.txt))
replays B and Q on Marrakesh [4, 5, 6] under the October snapshot and RO
matrix, through the patched replay of
[`reality_window_crosstalk_fit.py`](../simulations/reality_window_crosstalk_fit.py),
with multinomial covariance at 8192 shots per circuit. The letters used
below are yardsticks for sensitivity; none is a prediction of what the
chip carries.

- The baseline (no letter). B5(0)'s XX basis reads ⟨X₀⟩ +0.918, ⟨X₂⟩
  +0.974, ⟨X₁⟩ −0.980, ⟨X₀X₁X₂⟩ −0.881. In the Z bases ⟨Z₀⟩ reads +0.0715,
  mostly from the readout (the October RO matrix's asymmetry on q4,
  marginal over the other bits, is 0.066, as the October page's §7 gives)
  and the rest damping. A single-site reading is read against the day's RO matrix,
  never against zero.
- The RO matrix's own noise. Measured from eight circuits of 8192 shots,
  it moves the baseline's 14 readings by σ 0.0027 (median; largest 0.0040)
  against a shot σ of 0.011. Left out of the covariance, W of §5 would
  exceed its 99 % point under noise alone at a rate of 0.016 and the
  turn-averaged mean against the baseline at 0.057 (on all 36 forbidden
  readings, 0.048 and 0.36); with it in, the
  rate is 1 % by construction in these Gaussian draws, and the simulated
  flights below are the check.
- A frame-blind null the replay leaves out, for the covariance alone. A
  global depolarization of 1, 2, 4 and 10 %, applied to the parities at the
  end and drawn for its fluctuation only, keeps the false-positive rate of
  W and T₅ᴮ on the 14 readings at 1 % (0.006 to 0.013). All 36 forbidden
  readings would not: at 1 % depolarization T₅ᴮ would fire in 67 % of
  draws, at 4 % in all of them, through the repeated readings of §1. What
  incoherent gate errors inside the compiled circuit do to W's mean is not
  drawn here; by §1 they can open the window a little, and they are among
  the frame-blind carriers of §5.
- Frame-blind yardsticks. Under October's five Z letters, a detuning of
  0.135 rad/µs (21.5 kHz) on q4 or on q6, 0.1 rad/µs on q5, an sx
  and x amplitude error of 1 % on q4 or on q6, a phase of 0.03 after every
  cz, and the sx pulses on q4 and q6 tilted by 0.06 rad against the x
  pulses, B's forbidden readings are invariant under every turn (largest
  change 2.5 × 10⁻¹⁵) and W's noncentrality is 654, 142, 245, 33, 130, 6,
  1119 and 184 against the 99.9 % point 36.1 on 14 degrees of freedom. A
  tilt of 0.06 rad common to every sx and x pulse on q4 and q6 moves them
  by 1.1 × 10⁻¹⁵.
- A turning yardstick. Under October's crosstalk fit alone, the turn on q5
  moves B's forbidden readings by up to 0.164. On the 14 readings T₅ᴮ's
  noncentrality is 775 against the 99.9 % point 76.1 on 42 degrees of
  freedom, power 1.000, and so it stays with a drift factor of 1.27
  (October's over all its fit arms) or 1.96 (October's over arm A's
  readings alone; both from the October halves, part 7); the turn by π on q4: 51.0 against 36.1 on 14, power
  0.98 (0.92 and 0.61 with the two drift factors); on q6: 373, power
  1.000. W's noncentrality
  is 259. ⟨Z₂⟩ pooled reads +0.096, +0.071, −0.067, −0.049 at the four
  turns against a baseline of +0.013. The turn-averaged mean is not the
  baseline: a transverse letter carries a part independent of φ at second
  order, here up to 0.0014 per reading. Under all three October families
  T₄ᴮ is 8.8, power 0.06 (0.04 and 0.02).
- The reach. A crosstalk ε on (5 → 4) and (5 → 6) at the October phases,
  the other pairs at zero, gives T₅ᴮ power at its 99.9 % point of 0.31,
  0.81 and 0.99 at ε = 0.3, 0.4 and 0.5 % with the drift factor 1.27, and
  0.10, 0.42 and 0.82 with 1.96; the 50 % point lies near 0.34 % and
  0.42 %. These use one drift factor; the reach §5 reports is computed
  with §5's κΣ, which divides even a drift-free chip's readings by more
  than one, and on the paths frozen on 2026-11-13 it lies at 0.34 %
  (below). The phases are the October fit's, a yardstick and not
  a prediction.
- The drift factor estimated from the data. When the drift is real, a χ²
  test divided by one global factor is anticonservative, most of all when
  the drift sits on a few readings, as one qubit's wandering frequency
  would put it: with three of the 14 readings drifting at ten times their
  variance T₅ᴮ fires by chance at 13 % at its 99 % point and 7 % at its
  99.9 % point (part 7). §5 therefore weighs the drift block by block, and
  the drive-frame page's §5 gives the false-positive rates of the verdicts
  under drift drawn through the analysis itself. The rules of
  §5 call TURNS and OPEN only above the 99.9 % point.
- Q. The replay of Q under detunings of up to 0.2 rad/µs, RO-corrected and
  fitted per qubit by (δ, θ₀, contrast, decay rate), returns the detuning
  set, sign included; over 200 noisy redraws at 8192 shots the fitted δ
  scatter by 0.0032 to 0.0045 rad/µs about the set values. That fit is
  unweighted with the exact RO matrix; the weighted fit of §5 is checked by
  the simulated flights.
- At each chip's freeze. The numbers above are replayed under the October
  snapshot on Marrakesh [4, 5, 6]. At each freeze they are recomputed on
  the chip's path from its freeze snapshot and RO matrix with the gate
  durations of its target, every letter carried over by role
  ([`run_drive_frame.py`](../simulations/flight/run_drive_frame.py)
  `--freeze-numbers`, output in
  [`simulations/results/drive_frame/`](../simulations/results/drive_frame/)).
  On 2026-11-13 the rule chose Marrakesh [45, 46, 47] and Kingston
  [13, 14, 15] (the drive-frame page's §4 says why):

  | | Marrakesh [45, 46, 47] | Kingston [13, 14, 15] |
  |---|---|---|
  | B5(0), XX basis: ⟨X₀⟩, ⟨X₂⟩, ⟨X₀X₁X₂⟩ | +0.975, +0.987, −0.933 | +0.986, +0.990, −0.932 |
  | pooled ⟨Y₀⟩, ⟨Y₂⟩, ⟨Z₀⟩, ⟨Z₂⟩ | +0.011, +0.029, +0.006, +0.025 | +0.005, +0.007, +0.003, +0.003 |
  | W's noncentrality: five Z letters; detuning 0.135 rad/µs on a, on c; 0.1 rad/µs on b | 660; 207, 249; 42 | 649; 194, 243; 41 |
  | the same: amplitude +1 % on a, on c; cz phase 0.03; sx tilted 0.06 rad against x on a and c | 146, 7.0; 1188; 258 | 158, 7.5; 1275; 270 |
  | reach of T₅ᴮ at its 99.9 % point | ε = 0.343 % | ε = 0.334 % |
  | reach of a SHUT, detuning on c | 0.0334 rad/µs (5.3 kHz) | 0.0341 rad/µs (5.4 kHz) |

  The reaches are drawn with each simulated flight's own Σ and drift
  matrix, as a drift-free chip would give them, from 2000 flights per size.
- The simulated flights of the drive-frame page's §4, 10 000 per chip,
  drawn from the baseline replay's own distributions so that their mean is
  y₀ by construction, the RO matrix re-estimated in each draw, every draw
  through the analysis unchanged. On the paths of 2026-11-13, 95 / 99 /
  99.9 %:

  | statistic (degrees of freedom) | Marrakesh | Kingston | χ² |
  |---|---|---|---|
  | W (14) | 20.1 / 25.2 / 30.8 | 20.5 / 25.4 / 31.0 | 23.7 / 29.1 / 36.1 |
  | T₅ᴮ (42) | 49.8 / 57.6 / 65.5 | 50.4 / 57.5 / 66.9 | 58.1 / 66.2 / 76.1 |
  | T₄ᴮ (14) | 20.1 / 25.3 / 32.8 | 20.3 / 25.5 / 31.6 | 23.7 / 29.1 / 36.1 |
  | T₆ᴮ (14) | 20.2 / 24.7 / 31.4 | 20.1 / 24.9 / 32.1 | 23.7 / 29.1 / 36.1 |
  | Eᴮ (14) | 20.2 / 25.7 / 31.4 | 20.3 / 25.1 / 31.0 | 23.7 / 29.1 / 36.1 |

  Every simulated point lies below its χ² point, so the χ² points rule.
  The fitted δ scatter by 0.0039 to 0.0040 rad/µs on the path qubits
  against a median σ_δ of 0.0044 to 0.0046, a ratio of 0.86 to 0.87: σ_δ
  is not scaled.

## 5. Verdict rules

Every reading and RO matrix is formed as the drive-frame page's §5 forms
them, per chip. The baseline y₀ is the replay of a circuit on the chip's
path under its freeze snapshot with the day's RO matrix and no letter. Σ
is the multinomial covariance at the circuit's shots over the 14 readings
of §1, from the measured parities pooled over B's standing variants; Σ_RO
is the covariance
that the RO arm's counts give y₀ through the RO matrix (column by column,
multinomial, J·C·Jᵀ, times the RO arm's own drift factor κ_RO of the
drive-frame page's §5), one matrix for all the B variants, which are
logically one circuit. κΣ is the drive-frame page's, the drift weighed
block by block from the half differences of B's standing variants.
Every comparison against y₀ uses κΣ + Σ_RO; a difference between variants
of B uses κΣ alone, Σ_RO cancelling there. Pseudo-inverses keep
eigenvalues above 10⁻⁹ of the largest. The χ² points rule; where a 95 %,
99 % or 99.9 % point of the 10 000 simulated flights of §4 (numpy's
default quantile) exceeds its χ² point by more than 5 %, the simulated
point rules for that statistic and the record says so. Where the spread
of the simulated δ exceeds the median of the simulated σ_δ by more than
5 %, the day's σ_δ is multiplied by their ratio.

- **PB, the parse gate.** In B5(0)'s XX basis: ⟨X₀⟩ > 0.7, ⟨X₂⟩ > 0.7 and
  ⟨X₀X₁X₂⟩ < −0.6 (each freeze's baseline is in §4). PB asks for more contrast than the
  drive-frame page's P: in its drift check PB fails at a contrast of 0.70
  in 94 to 98 % of flights, where P passes. A chip that fails
  it has no B verdicts. A variant whose turn is missing from the circuits
  retrieved after the flight is VOID.
- **WB, the window.** W = rᵀ(κΣ + Σ_RO)⁺r over B5(0)'s 14 readings of §1,
  r = measured minus y₀, on the rank: OPEN above its 99.9 % point, SHUT
  below its 95 % point, INCONCLUSIVE between. SHUT is reported with its
  reach: the detuning on c (logical q2) at which W, with the day's weights,
  passes its 95 % point in half of 2000 simulated flights drawn under it
  from the baseline replay. W restricted to each of the three
  kinds of §1, the other 22 forbidden readings and the four pooled
  single-site readings are reported with their σ and no verdict.
- **VFB, the frame.** T₅ᴮ, T₄ᴮ, T₆ᴮ exactly as the drive-frame page's T₅,
  T₄, T₆, on B's 14 readings: for each turned qubit, TURNS above the
  99.9 % point, FRAME-BLIND at or below the 99 % point (on a and c, for
  the cosine part alone), INCONCLUSIVE between, VOID with fewer than two
  standing variants. These verdicts stand whatever WB reads: a turning
  letter in quadrature at φ = 0 leaves B5(0) at its baseline and still
  turns. FRAME-BLIND also covers a turn the control stack dropped below the
  submitted circuit, which no logical reading can show. The reach of a
  FRAME-BLIND on b is reported as the ε of §4's yardstick at which T₅ᴮ,
  formed with the day's κΣ, passes its 99.9 % point
  in half of 2000 simulated flights drawn under that letter from the
  baseline replay (with the October phases carried over by role), a 50 %
  point and not a bound; it concerns that letter form alone, not a → b, c → b, an exchange or a leak, and its
  phases are a fitted letter's, a yardstick. On a and c no reach is
  reported: one turn by π reads the cosine part alone. The even harmonic
  Eᴮ, formed as the drive-frame page's E on the 14 readings, is reported
  as 2φ PRESENT above its 99.9 % point, beside TURNS
  and not as its cause.
- **VQ, the idle detuning.** Read whether or not PB passes. Every Q circuit
  measures all three path qubits and the forest, with the drive-frame
  page's count-key layout. Q's τ delays must stand in the frozen circuits
  and in the circuits retrieved after the flight; a Q circuit without them
  is VOID. Each Q circuit's outcome distribution, from its measured counts,
  is corrected by the inverse of the day's RO matrix. For each
  path qubit its eight readings, from its own eight circuits, are fitted by
  ⟨X⟩ = C·e^(−gτ)·cos(θ₀ + δτ), ⟨Y⟩ = C·e^(−gτ)·sin(θ₀ + δτ), starting at
  δ = 0, θ₀ from τ = 0, C = 1, g = 0, bounded to |δ| ≤ 1.5 rad/µs (the 1 µs
  steps alias at 2π rad/µs), 0 ≤ C ≤ 1.2, 0 ≤ g ≤ 1 per µs, |θ₀| ≤ 10,
  scipy's least_squares at its default tolerances, weighted by the
  covariance of the corrected readings: the multinomial covariance of the
  measured counts carried through the inverse RO matrix A⁻¹ (as A⁻¹ C A⁻ᵀ
  on the outcome probabilities, then onto the parity), scaled by κ_Q
  (max(1, the mean over all 24 readings of (half₁ − half₂)²/(var₁ +
  var₂))); plus a phase that wanders from PUB to PUB, which moves a reading
  by its angle times τ perpendicular to the Bloch vector (d⟨X⟩ = −⟨Y⟩ dθ,
  d⟨Y⟩ = ⟨X⟩ dθ): per qubit one variance s² per unit τ², the excess of
  its six readings' (half₁ − half₂)² with τ > 0 over var₁ + var₂, summed,
  over the sum of 2τ²g² (g the perpendicular reading, ⟨Y⟩ for an X reading
  and ⟨X⟩ for a Y reading; s² at least zero), entering each reading as
  s²τ²g²/2; plus the RO arm's own noise, times κ_RO, carried through A⁻¹
  by the derivative of A⁻¹p in each column of A. δ is reported
  with σ_δ from the weighted Jacobian, inflated by √(χ²_fit/4) where that
  exceeds one, and with the fit's χ². A qubit whose fit ends on a bound of
  δ, does not converge, or has a VOID circuit is VOID. The verdict, per
  qubit: DETUNED when |δ| ≥ 4σ_δ, NOT DETUNED when |δ| < 2σ_δ,
  INCONCLUSIVE between. The wander rests on six half differences, which
  gives z a tail beyond the normal's; the rates at which DETUNED is called
  by chance are in the drive-frame page's §5. The
  replay of B5(0) with the three δ as its only letters is reported beside
  W, as a reading: the replay's detuning acts through every gate as well as
  every idle, and Q sees only the idle.

| WB | VFB on a qubit | reading for that qubit |
|---|---|---|
| OPEN | TURNS | something that turns with its frame stands in the window: a foreign transverse letter coupling to it, an image term of its own drive, or leakage or a residual swap at a cz (2φ PRESENT is read beside it, not as its cause); whether it opens the window whole is not decided |
| OPEN | FRAME-BLIND | nothing that turns with its frame is seen (on b at the reported reach for b's pulses reaching a and c; on a and c for the cosine part alone); Z letters, a static ZZ, own-frame letters, incoherent gate errors, damping and the readout are not separated from one another |
| OPEN | INCONCLUSIVE or VOID | nothing for that qubit |
| SHUT or INCONCLUSIVE | TURNS | something turns with its frame without an excess over the baseline in B5(0), as a letter in quadrature at φ = 0 would |
| SHUT or INCONCLUSIVE | FRAME-BLIND, INCONCLUSIVE or VOID | no reading beyond the classes themselves |

The mean of the 14 readings over the four B5 turns is reported
against y₀ with its χ² under κΣ/4 + Σ_RO and no verdict: a transverse letter puts part of itself there
at second order (§4), so that mean is not a frame-blind share. The two
chips' classes are tabled side by side with no verdict: a device-anchored
carrier can read TURNS on both chips with different sizes and signs, and
two different Z letters can read FRAME-BLIND on both. Over both pages a
Marrakesh or Kingston chip carries thirteen verdicts and one more per
forest qubit (V1, T₅, T₄, T₆, E and VS on each forest qubit; WB, T₅ᴮ,
T₄ᴮ, T₆ᴮ, Eᴮ and VQ on each path qubit), seventeen on Marrakesh
[45, 46, 47] and fifteen on Kingston [13, 14, 15], the paths of
2026-11-13, each firing by chance at about 1 % or less under the rules
above (under drift the drive-frame page's §5 gives up to 0.9 % for B's
frame statistics, 0.7 % for W, 1.3 % for VQ and 2.9 % for V1 when a drift
falls on S), so one of them fires on a
chip by chance at up to about 15 %, and a lone firing is read that way.
The forest bits of B and Q are summed
out of the path's parities; the forest's P(1) in B and in Q is reported
without a verdict, VS being the drive-frame page's, on its A circuits.

⟨X₀X₁X₂⟩ of B5(0) in its XX basis, the colouring of §1, is reported beside
its baseline replay, without a verdict. In the logical evolution it
commutes with H and pays the dephasing of its three sites; here the
noiseless state is a product and ⟨X₀X₁X₂⟩ the product of its single-site
readings, unlike block C of
[PRICE_PAIR_HARDWARE_PREDICTION](PRICE_PAIR_HARDWARE_PREDICTION.md), whose
GHZ ⟨XXX⟩ outlived that product. Everything else is read without a verdict; a single
reading beyond 3σ is not a finding, and beyond 4σ it is reported as a
sighting.

## 6. Flight conditions

The drive-frame page's §6 holds for every job that carries B and Q. Its
freezer builds them from the April QPY on Marrakesh and Kingston into the
same frozen file, asserts their turns as for its own arms and asserts Q's
idles. The record goes with the drive-frame page's, in the same counts
files. The pre-registration (this page through §6, with §4's
numbers those of the freezes that fly, and the drive-frame page with its hashes) is
committed together with the freezer, the runner and the frozen files
before the first submission, and its commit hash is written into §7.

## 7. The flight

Not flown.
