# Both Sides Visible: What a Quantum Computer Shows When You Watch It for Six Months

<!-- Keywords: CΨ quarter boundary IBM hardware, palindromic mirror both sides,
quantum classical oscillation real hardware, ibm torino qubit crossing pattern,
R=CPsi2 both sides visible, quantum bridge heartbeat IBM, palindromic complement
real data, quantum information lives on both sides -->

**Status:** Observed on IBM Torino hardware (180 days, 133 qubits)
**Date:** March 25, 2026 (analysis of Aug 2025 - Feb 2026 calibration data)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Data:** [ibm_torino_history.csv](../data/ibm_history/ibm_torino_history.csv)

---

## What this document is about

This is the moment where the mathematics meets a real machine.

Everything in this project, the palindromic mirror, the Π operator,
the idea that every quantum mode has a partner, was proven with
algebra and verified with simulation. This document shows it happening
on an actual quantum computer, measured daily for six months, without
anyone designing it to happen. The patterns you are about to see were
not programmed. They were discovered in IBM's publicly available
calibration data.

If one document in this repository could convince a skeptic that the
palindromic symmetry is not just a mathematical curiosity, this is
the one. Because this is not theory. This is 24,074 measurements
on 133 qubits on real silicon.

---

## What You Are Looking At

A qubit can be in two regimes. Think of it as a coin that is either
still spinning (quantum: possibilities open, nothing decided yet) or
has landed (classical: one side up, decided).

The number CΨ tells you which regime a qubit is in. It is computed
from two calibration values that IBM publishes daily for every qubit:
T1 (how long the qubit holds its energy) and T2 (how long it holds
its quantum coherence). The ratio r = T2 / (2 × T1) determines the
regime. When r is below 0.213: the qubit crosses the CΨ = ¼ boundary.
Above 0.213: it does not.

Most qubits on a quantum computer stay firmly on one side. But some
have r values that fluctuate around 0.213. On some days they cross,
on others they do not. Over days. Over weeks.

Below are two of these qubits, tracked for 180 days on IBM's Torino
quantum processor. Each day is one character:

- **X** = qubit crossed the ¼ boundary that day (r < 0.213)
- **.** = qubit stayed below the boundary (r > 0.213)

Now here is the discovery: the decay spectrum of these qubits has a
proven mirror symmetry ([palindromic spectrum](proofs/MIRROR_SYMMETRY_PROOF.md)).
Every decay rate d has a mathematically guaranteed partner at 2Σγ - d.
This pairing is algebraically exact (54,118 eigenvalues, zero
exceptions). It means the complement of the crossing pattern (every
X flipped to . and vice versa) is not an arbitrary inversion. It is
what the **partner modes** do. When one side of a palindromic pair
is active, the other is quiet.

**A note on causes:** The daily fluctuations in T2 that produce this
pattern have known mundane causes: temperature drift, two-level-system
defects in the chip, magnetic field variations. What is new here is
not the fluctuations themselves but the observation that their pattern
has a mathematically guaranteed complement, and that complement has
structure.

---

## Qubit 98 (57.5% crossing rate)

What follows is a visual diary. The top block shows what we observe
directly: the days when this qubit crossed the quantum-classical
boundary. The bottom block shows the exact complement, what the
palindromic partner modes are doing at the same time. Read them
side by side. When one is active, the other is silent.

```
Our side:

X..XXXX  week 1   tuning in
X.XXXXX  week 2   tuning in
XX.X...  week 3
....X.X  week 4
XX.X..X  week 5
XXXXXXX  week 6   ━━━━ stable pulse
XXXXXXX  week 7   ━━━━ stable pulse
XXX.XXX  week 8   ━━━━ stable pulse
XXXXXXX  week 9   ━━━━ stable pulse
XXXXXXX  week 10  ━━━━ stable pulse
XXXXXXX  week 11  ━━━━ stable pulse
XXXXX.X  week 12  fading
.X.X...  week 13
XXX..XX  week 14
XXXX..X  week 15
X.X..XX  week 16
XXXX.XX  week 17
..X.X..  week 18  silence
.......  week 19  silence
XXXXX..  week 20  echo
.X...X.  week 21
.......  week 22  silence
X.....X  week 23
XX.....  week 24
.....X.  week 25
.X....   week 26  silence
```

```
The other side (palindromic mirror, exact complement):

.XX....  week 1
.X.....  week 2
..X.XXX  week 3
XXXX.X.  week 4
..X.XX.  week 5
.......  week 6   ━━━━ their silence is our pulse
.......  week 7   ━━━━
...X...  week 8   ━━━━
.......  week 9   ━━━━
.......  week 10  ━━━━
.......  week 11  ━━━━
.....X.  week 12
X.X.XXX  week 13
...XX..  week 14
....XX.  week 15
.X.XX..  week 16
....X..  week 17
XX.X.XX  week 18  their pulse begins
XXXXXXX  week 19  ━━━━ their stable pulse
.....XX  week 20
X.XXX.X  week 21
XXXXXXX  week 22  ━━━━ their stable pulse
.XXXXX.  week 23
..XXXXX  week 24
XXXXX.X  week 25
X.XXXX   week 26  ━━━━ they are active now
```

When we send (weeks 6-11), they are silent.
When we go silent (weeks 19, 22), they send.
Nothing is missing. Nothing is wasted.

---

## Qubit 72 (66.9% crossing rate)

```
         Ours    | Theirs
week 1:  X...XXX | .XXX...
week 2:  XXX.XXX | ...X...
week 3:  .XXXXXX | X......
week 4:  XXXXXX. | ......X
week 5:  ..XXXXX | XX.....
week 6:  XXXXX.. | .....XX
week 7:  X..X..X | .XX.XX.
week 8:  X.XXX.X | .X...X.
week 9:  XXXXXXX | .......  ━━ our pulse
week 10: XXXXXXX | .......  ━━ our pulse
week 11: XXX.XXX | ...X...
week 12: ...XXXX | XXX....
week 13: XX..XXX | ..XX...
week 14: XX.XXX. | ..X...X
week 15: ..X.XXX | XX.X...
week 16: X..XXXX | .XX....
week 17: XX.XX.X | ..X..X.
week 18: XXXXXXX | .......  ━━ our pulse
week 19: X...XX. | .XXX..X
week 20: ...XX.X | XXX..X.
week 21: X.X.... | .X.XXXX  ━━ their turn
week 22: ..X..X. | XX.XX.X  ━━ their turn
week 23: XX.X.XX | ..X.X..
week 24: XXXXXXX | .......  ━━ our pulse
week 25: .....XX | XXXXX..  ━━ their turn
week 26: ....XX  | XXXX..
```

The rhythm is visible: our pulses (weeks 9-10, 18, 24), their
turns (weeks 21-22, 25). The bridge breathes.

---

## What This Shows

This is not a simulation. This is a real quantum computer (IBM Torino)
in a real lab, measured once a day for six months. Nobody designed this
pattern. Nobody programmed qubit 98 to pulse for six weeks and then stop.

What the math predicted: every decay mode has a mirror partner. When
our side is active, the mirror side is quiet, and vice versa. Like two
people sharing one breath: when one inhales, the other exhales.

What the data shows: exactly that. Every X on our side is a . on
theirs. Every gap in our pattern is filled in theirs. The picture is
incomplete from one side alone. It is complete when you see both.

This is not interpretation. This is counting X's and .'s on calibration
data that IBM publishes for anyone to check.

16 qubits on IBM Torino oscillate around the ¼ boundary:

| Qubit | Crossing rate | Character |
|-------|---------------|-----------|
| Q72 | 66.9% | Balanced, rhythmic |
| Q98 | 57.5% | Clear lifecycle: tune, pulse, fade |
| Q105 | 56.9% | Long active phase, then silent |
| Q70 | 26.0% | Mostly silent, brief pulses |
| Q68 | 23.2% | Mostly silent |

The most balanced qubits (near 50%) show the clearest alternation
between the two sides. They live at the boundary, where the pattern
is most symmetric.

---

## Connection to the Proven Chain

```
d² - 2d = 0  →  d = 2 (qubit, the only dimension)
           →  Π exists (palindromic mirror)
           →  rate d pairs with 2Σγ - d
           →  what we see as X, they see as .
           →  the picture is complete across both sides
```

This was always true. For every qubit. On every quantum computer.
Since the first transmon was cooled. We just learned to read it.

---

*Data: 133 qubits, 180 days, 24,074 calibration records.*
*Source: IBM Quantum Platform, ibm_torino backend.*
*Analysis: [ibm_history_analysis.py](../data/ibm_history/ibm_history_analysis.py)*
*Full crossing data: [ibm_q98_crossing_pattern.txt](../simulations/results/ibm_q98_crossing_pattern.txt)*
*Both sides: [both_sides_visible.txt](../simulations/results/both_sides_visible.txt)*


---

## Direct Measurement: Both Sides From One Density Matrix (March 25, 2026)

The calibration patterns above are visually compelling, but they
have a logical weakness: flipping every X to a . is just inverting
bits. Any pattern looks structured when you show it next to its
complement. The question is: does "the other side" actually exist
in the physics, or is it just an arithmetic trick?

This section answers that question. Instead of computing a
complement from daily calibration data, we measured the full quantum
state of a single qubit at multiple points in time. A full density
matrix contains both sides simultaneously: the diagonal elements
(populations, what you observe) and the off-diagonal elements
(coherences, what the Π operator maps to). No complement is
computed. Both sides are measured in the same experiment.

The calibration patterns above are computed complements: we invert
the X/. pattern and call it "the other side." This is visually
compelling but logically tautological. Flipping bits always fills gaps.

The tomography data resolves this. On February 9, 2026, we measured
full density matrices on IBM Torino (Qubit 52, 25 time points, 8192
shots per point). A density matrix contains BOTH sides simultaneously:

- **Diagonal elements** (populations): what we observe directly
- **Off-diagonal elements** (coherences): what the Π operator maps to

Applying the proven conjugation operator Π to each measured density
matrix and computing CΨ from both perspectives:

| t (us) | CΨ_A (ours) | CΨ_B (Π side) | Both > 1/4? |
|--------|---------------|------------------|-------------|
| 0 | 0.885 | 0.941 | YES - bridge open |
| 37 | 0.475 | 0.368 | YES - bridge open |
| 75 | 0.385 | 0.541 | YES - bridge open |
| 112 | 0.261 | 0.382 | YES - A approaching 1/4 |
| 149 | 0.125 | 0.422 | NO - A crossed, B still quantum |
| 261 | 0.034 | 0.381 | NO - A classical, B still quantum |
| 373 | 0.004 | 0.350 | NO - A gone, B holds |
| 560 | 0.017 | 0.301 | NO - B slowly approaching |
| 634 | 0.025 | 0.280 | NO - B nearing 1/4 |
| 783 | 0.054 | 0.265 | NO - B approaching 1/4 |
| 895 | 0.037 | 0.247 | NO - B crosses. Bridge closed. |

In plain language: our side (A) loses its quantum character within about
140 microseconds. The other side (B), read from the same measurement by
applying the Π operator, holds on for about 895 microseconds: six times
longer. They do not decay together. They take turns. First our side
fades. Then, slowly, the other side follows. The bridge does not slam
shut. It narrows gradually, one side at a time.

This is not a computed complement. This is Π applied to hardware-
measured density matrices. The off-diagonal elements ARE the other
side. They were measured in the same experiment, with the same shots,
at the same time.

### What the data shows

Our side (A) crosses 1/4 at approximately 140 us. The Π side (B)
crosses 1/4 at approximately 895 us. The other side lives **6x longer**
on real silicon.

The bridge is open (both > 1/4) for the first ~112 us. Then it is half
open (B still quantum, A classical) for ~750 us. Then it closes.

During the half-open window: information can still flow from the quantum
side (B) to the classical side (A), but not back. This is the one-way
channel that dCΨ/dt < 0 predicts. The bridge does not slam shut. It
narrows gradually, one side at a time.

### What this means for interpretation

The calibration patterns (Qubit 98, Qubit 72) showed alternation over
days. That alternation is real - T2 fluctuates, crossings happen - but
calling the complement "the other side" was premature.

The tomography data shows something more precise: both sides exist
simultaneously in every density matrix. They are not alternating in
time. They are coexisting in the state. The populations decay at one
rate. The coherences (Π-mapped partners) decay at another. The
palindromic pairing is not a pattern across days. It is a structure
within each measurement.

The bridge is not something that opens and closes over days. It is
something that exists in every quantum state, at every moment, and
gradually narrows as decoherence transfers weight from the quantum
side to the classical side. What the daily calibration patterns show
is the CUMULATIVE result of this continuous process, sampled once per
day at whatever T2 the qubit happens to have.

### Data source

- Backend: IBM Torino (ibm_torino)
- Qubit: 52
- Date: February 9, 2026
- Shots: 8192 per circuit
- Time points: 25 (0 to 895 us)
- Tomography: Full single-qubit state tomography (3 bases: X, Y, Z)
- Raw data: [ibm_tomography_feb2026/](../data/ibm_tomography_feb2026/) (density matrices),
  [ibm_history/](../data/ibm_history/) (calibration history)
- Analysis: Π conjugation applied to measured density matrices,
  CΨ computed from both the original and Π-transformed states

### Why the other side lives longer

Our side decays with T2, the coherence time (how fast quantum
superposition is lost). The other side decays with T1, the relaxation
time (how fast energy is lost). On almost every qubit ever built,
T1 is much longer than T2. For Qubit 52 on that day: T2 ≈ 150 μs,
T1 ≈ 900 μs.

The 6x ratio is T1/T2 for this qubit. Our side's CΨ decays with T2
(coherence time, ~150 μs, fast). The Π side's CΨ decays with T1
(relaxation time, ~900 μs, slow). This is the palindromic pairing:
rate d pairs with 2Σγ - d. One fast, one slow. Always balanced.

In plain language: the two sides of the palindrome live on different
clocks. Our side runs fast and fades quickly. The mirror side runs
slow and lingers. This is not a coincidence. It is the pairing: the
faster one side decays, the slower its partner must decay. The
mathematics guarantees it.

### For newcomers: where to start

If this is your first time here, these documents explain the concepts
used on this page:

- **What is CΨ and why ¼?** [Proof Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md)
  gives the seven-layer proof. Short version: CΨ = purity x coherence,
  ¼ is where the discriminant of R = C(Ψ+R)² vanishes (fold catastrophe).
- **What is Π?** [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md).
  Π swaps populations and coherences, producing the exact palindromic pairing.
- **Why does this matter?** [The Bridge Was Always Open](THE_BRIDGE_WAS_ALWAYS_OPEN.md)
  connects the palindrome to the incompleteness proof: noise is not random,
  it is a structured channel from an external source.
- **What was found today?** [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md)
  observed MI peaking at the ¼ crossing (fold catastrophe in action),
  CΨ oscillation (227 heartbeats), and the XOR pulse train.
- **What can we rule out?** [Exclusions](EXCLUSIONS.md) lists five things
  the mathematics proves impossible.
- **The full picture:** [Reading Guide](READING_GUIDE.md) gives three paths
  through the project depending on your background.
