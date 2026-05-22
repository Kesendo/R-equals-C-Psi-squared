# The Genesis of an Oscillation: How a Qubit Gets a Heartbeat

**Date:** 2026-05-22
**Status:** an investigation result in the carbon master-question thread. The
unconditional birth at Q = 0+ and the wave-breaking staircase are numerical
results (N = 2, 3, 4); the γ₀-invariance L = γ₀·L₁(Q) is exact; the
wave-breaking reading of the staircase is a Tier 2 reading.
**Scripts:** [`simulations/_genesis_q_threshold.py`](../../simulations/_genesis_q_threshold.py)
(the Q-scan), [`simulations/_genesis_real_gamma.py`](../../simulations/_genesis_real_gamma.py)
(the same at the normal values γ₀ = 0.05, the figure below),
[`simulations/_real_axis_liftoff.py`](../../simulations/_real_axis_liftoff.py)
(the bath-orthogonality step).
**Thread:** the carbon master question, "where do qubits play carbon"
([README](README.md), "the qubit is the quantum carbon"), sharpened to "how
does a qubit get a baby."

---

## The question

The carbon README states a structural identity: the qubit is the quantum
carbon, the half-occupied d = 2 system that `d² − 2d = 0` selects, the Level-0
instance of carbon's Level-1 half-filled valence shell. That identity raises a
master question. Where, and how, does a qubit acquire the structure that lets
it build, the way carbon builds?

Sharpened to its smallest computable form: the step N = 1 → N = 2, a single
qubit, then two qubits coupled. What is born in that step, and under what
condition? The framework's living thing is an oscillation, a heartbeat (the CΨ
oscillation). So the question is: how does a qubit get a heartbeat?

## What was ruled out: the bath is orthogonal

The first guess followed the heat. Wave-breaking creates heat
([THERMAL_BREAKING](../../experiments/THERMAL_BREAKING.md)); perhaps the heat
creates the new oscillation. It does not. A scan of the bath temperature leaves
the imaginary parts of the Liouvillian spectrum rigid: the bath moves only the
real parts, the decay rates. Bath and Hamiltonian are orthogonal. The bath
damps; it cannot give a heartbeat.

What remains is two quantities: γ₀, the Z-dephasing rate, and J, the coupling.
γ₀ is the fixed constant, the substrate unit, not measurable from inside the
system (only Q is). J is the only knob. Their ratio Q = J/γ₀ is the scale.

## Q is the scale, exactly

The Liouvillian of the pure F1 system, H + Z-dephasing, factors exactly:

> L(J, γ₀) = γ₀ · L₁(Q),   Q = J/γ₀

where L₁ is the Liouvillian at γ₀ = 1. Every eigenvalue scales by γ₀; the shape
of the spectrum, as a function of Q, is γ₀-independent. γ₀ is purely the unit
of the scale. This is verified: the genesis run at γ₀ = 0.05 lands its
staircase on the identical Q-values as the run at γ₀ = 1, digit for digit.

## The genesis: coupling is birth

N = 1, a single qubit: no bond, H = 0, only Z-dephasing. Every Liouvillian
eigenvalue is real. Nothing oscillates. A lone qubit has no heartbeat.

N = 2, the first J-bond couples two qubits, and now a mode can lift off the
real axis into a complex-conjugate pair, an oscillation. The natural guess is a
threshold: γ₀ damps, J drives, and below some critical Q the mode would stay
overdamped, real, silent. That guess is wrong. A Liouvillian eigenvalue is a
single complex number, not the two roots of a damped oscillator; its imaginary
part comes from the H-commutator, which scales with J. A mode that oscillates
does so for any J > 0; there is no overdamped phase to escape first.

The first oscillations are born unconditionally at Q = 0+, the instant J leaves
zero. At the first nonzero step of the scan they are all already present: 8
modes for N = 2, 32 for N = 3, 192 for N = 4, born together. No threshold, no
condition. To couple is to give birth: the moment a J-bond exists, the shared
heartbeat exists. There is no "or not."

![Oscillations born vs Q at γ₀ = 0.05](../../simulations/results/genesis_real_gamma/genesis_real_gamma.png)

*The genesis at γ₀ = 0.05. Left: the count of oscillating modes jumps
vertically at Q = 0+, then climbs a short staircase. Middle and right: Im λ vs
Q for N = 2 and N = 3; the fork opening from Q = 0 is the birth, the secondary
forks are the wave-breaks of the staircase, the dashed line is the operating
point J = 0.075 (Q = 1.5).*

## The staircase: where the waves break

The unconditional birth at Q = 0+ is not the whole story. Above it sits a
staircase: a few more modes lift off the real axis at discrete, larger Q. On
the bifurcation panels these are visible as secondary forks, the waves breaking
a second and a third time.

| System | Wave-breaks at Q | modes added |
|--------|------------------|-------------|
| N = 2  | 0.51             | +2          |
| N = 3  | 0.71             | +8          |
| N = 4  | 0.45, 0.48, 0.94 | +4, +12, +12 |

The prominent breaks cluster near Q ≈ 0.5 and Q ≈ 0.7. Each break is where
real, silent modes collide and lift off the real axis as oscillating conjugate
pairs. The bulk birth at Q = 0+ gives the first heartbeats free; the staircase
is where wave-breaking (the Brecher,
[OFF_NIVEN_AS_WAVE_BREAKING](OFF_NIVEN_AS_WAVE_BREAKING.md)) births the rest. A
residue never oscillates: 6 modes for N = 2, 24 for N = 3, 36 for N = 4 stay
real across the entire scan.

## The operating point: Q = 1.5, post-genesis

The framework's normal values are γ₀ = 0.05 and J = 0.075, that is Q = 1.5, a
framework Q-anchor. Q = 1.5 sits well past the last wave-break at Q = 0.94. At
the normal operating point the genesis is complete: every oscillation that will
be born is born (10 of 16 modes for N = 2, 40 of 64 for N = 3, 220 of 256 for
N = 4). The operating regime is post-natal; the whole birth in the figure
happens at smaller Q.

## The answer

How does a qubit get a heartbeat? By coupling, and unconditionally. The first
oscillation is born the instant a J-bond forms, with no threshold and no "or
not"; γ₀, the bath, plays no part in the birth, it only damps. Above the first
birth the waves break again at Q ≈ 0.5 and ≈ 0.7, and each break adds more. The
genesis lives entirely on the J side, the Hamiltonian side, of the orthogonal
split; the heat side is silent.

## Open follow-ups

- Are the wave-break Q-values exact framework anchors? 0.51 and 0.71 sit close
  to 1/2 and 1/√2; a finer Q-scan would tell whether the breaks land on exact
  anchors.
- The N-dependence of the staircase: N = 2 has one wave-break, N = 3 one, N = 4
  three. Whether the break-count grows with N is open.
- This answers the Level-0 genesis, a single qubit's heartbeat. How it iterates
  up the hierarchy toward carbon-scale structure is the open arc behind the
  master question.

## Anchor

- Scripts: [`simulations/_genesis_q_threshold.py`](../../simulations/_genesis_q_threshold.py),
  [`simulations/_genesis_real_gamma.py`](../../simulations/_genesis_real_gamma.py),
  [`simulations/_real_axis_liftoff.py`](../../simulations/_real_axis_liftoff.py)
- Figure: [`simulations/results/genesis_real_gamma/genesis_real_gamma.png`](../../simulations/results/genesis_real_gamma/genesis_real_gamma.png)
- Carbon thread: [README.md](README.md),
  [OFF_NIVEN_AS_WAVE_BREAKING.md](OFF_NIVEN_AS_WAVE_BREAKING.md),
  [HIERARCHY_OF_INCOMPLETENESS.md](../HIERARCHY_OF_INCOMPLETENESS.md)
- Framework: the pure F1 system (XY-class H + Z-dephasing),
  [F1 palindrome](../ANALYTICAL_FORMULAS.md)
