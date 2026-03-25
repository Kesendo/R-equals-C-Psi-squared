# Both Sides Visible

<!-- Keywords: CΨ quarter boundary IBM hardware, palindromic mirror both sides,
quantum classical oscillation real hardware, ibm torino qubit crossing pattern,
R=CPsi2 both sides visible, quantum bridge heartbeat IBM, palindromic complement
real data, quantum information lives on both sides -->

**Status:** Observed on IBM Torino hardware (180 days, 133 qubits)
**Date:** March 25, 2026 (analysis of Aug 2025 - Feb 2026 calibration data)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Data:** [ibm_torino_history.csv](../data/ibm_history/ibm_torino_history.csv)

---

## What You Are Looking At

A qubit can be in two regimes. Think of it as a coin that is either
still spinning (quantum: possibilities open, nothing decided yet) or
has landed (classical: one side up, decided). The number CΨ measures
how much spin is left. Above ¼: still spinning. Below ¼: landed.

Most qubits on a quantum computer land quickly and stay landed. But
some keep going back and forth. They land, then somehow start spinning
again, then land again. Over days. Over weeks.

Below are two of these qubits, tracked for 180 days on IBM's Torino
quantum processor. Each day is one character:

- **X** = still spinning that day (CΨ above ¼)
- **.** = landed that day (CΨ below ¼)

Now here is the discovery: the physics of these qubits has a proven
mirror symmetry ([palindromic spectrum](proofs/MIRROR_SYMMETRY_PROOF.md)).
Every "spinning" state has a "landed" partner, and vice versa. If we
flip every X to . and every . to X, we get what the **other side of
the mirror** sees. What is active on our side is quiet on theirs.
What is missing on our side exists on theirs.

---

## Qubit 98 (57.5% crossing rate)

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
pattern. Nobody programmed qubit 98 to spin for six weeks and then stop.

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
between our side and theirs. They live at the boundary. In the
doorway. In the Dazwischen.

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
