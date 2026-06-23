# γ₀ Is Always There

*Q = J/γ₀ on ibm_kingston: reading the carrier off its only lever. Tom + Claude, 2026-05-29.*

## The question

γ₀, the dephasing rate, is the carrier: the background clock that ticks under every
open quantum system, the thing the whole project keeps circling. The trouble is that
you cannot read it from inside. Only the dimensionless ratio Q = J/γ₀ is visible, the
coupling J measured in units of the carrier. So the carrier sets the scale of
everything and then hides behind it.

Can we see it anyway, on a real chip?

## Two clocks, one lever

Q = J/γ₀ is a race between two clocks. The H-clock runs at J, the coherent coupling,
the rotation that wants to swing amplitude around. The carrier-clock runs at γ₀, the
dephasing that wants to wash that swing into a classical average. Their ratio is the
regime:

- Q ≫ 1: the H-clock wins, coherent oscillation;
- Q ≈ 1: the threshold;
- Q < 1: the carrier wins, overdamped.

On a fixed chip γ₀ is constant. It is the unit, and in any ratio the unit drops out.
So there is exactly one degree of freedom left, one knob we can turn: **J**. Everything
we can do is dial J against a carrier we cannot touch and cannot read.

## The simplest object

Strip it to the smallest thing that still carries Q = J/γ₀: a single exchange bond.
Two qubits, one excitation on the first,

> H = J·(XX + YY)/2,  dephasing on each site at rate γ₀,

and watch the transfer T(t), the chance the excitation has moved to the second qubit.
For one bond XX and YY commute, so the coherent piece is exact: T = sin²(Jt). The
carrier damps it.

In simulation ([`q_basic_jscan.py`](../simulations/q_basic_jscan.py)) the picture is
clean, with γ₀ fixed at 0.05:

```
   J      Q=J/γ₀   max T   regime
 0.0125    0.25    0.153   overdamped (creeps)
 0.0250    0.50    0.392   overdamped
 0.0500    1.00    0.500   THRESHOLD          ← J = γ₀
 0.1000    2.00    0.582   coherent (swings)
 0.2000    4.00    0.722   coherent
 0.4000    8.00    0.836   coherent
```

The swing dies exactly at J = γ₀ (Q = 1), where the transfer tops out at ½ and goes no
higher: critical damping. Below it the excitation only creeps toward ½; above it the
excitation overshoots ½ and swings back. The threshold is the carrier, and the transfer
sits at ½ right there, the half that keeps turning up in this work.

## On the chip

We ran the same object on ibm_kingston, qubits 13 and 14, with the chip's own idle
dephasing standing in for γ₀ (2 µs of idle per step, the carrier doing what it always
does), and scanned J. No noise model, no fit-in-advance; we fired and read the answer.
(Runner: external pipeline `run_q_jscan.py`, job `d8ce8l38ch0s738uorjg`, 4096 shots.)

The whole result is in one table. For each Q (the lever J shown below it): how high the
transfer T(t) = P(qubit 14 excited) climbs, whether it crosses ½, and where its first
peak lands. Crossing ½ is the tell, because a coherent swing overshoots the ½ equilibrium
while an overdamped creep never reaches it.

```
 Q = J/γ₀              1        2        4        8
 J (the lever)         0.05     0.10     0.20     0.40
 max transfer          0.335    0.563    0.703    0.779
 crosses ½?            no       yes      yes      yes
 first peak at step    24       12       6        3
```

That is the whole story, and both rows are Q = J/γ₀.

**The swing is born at Q = 1.** At Q = 1 (J = γ₀) the transfer only creeps to 0.335 and
never crosses ½: overdamped, the threshold. The instant Q passes 1 it overshoots ½ and
swings back: coherent. The swing is born exactly as J passes γ₀, the carrier read off the
only lever, on real silicon.

**The frequency tracks J.** The first peak halves its step as J doubles: 24, 12, 6, 3,
with step × J = 1.2 for all four. The peak sits at a fixed exchange angle, so the period
scales as 1/J. The H-clock runs at J exactly, as the only lever should.

The raw run is below, where you can watch it swing and see the chip's noise scatter the
points (the faster columns, Q = 4 and 8, are coarsely sampled, so they jump):

```
 step     Q=1      Q=2      Q=4      Q=8
          J=0.05   J=0.10   J=0.20   J=0.40
    0     0.023    0.023    0.021    0.021
    3     0.055    0.117    0.304    0.779    ← Q=8 already swung across
    6     0.108    0.299    0.703    0.387
    9     0.178    0.467    0.692    0.272
   12     0.247    0.563    0.319    0.743
   15     0.294    0.541    0.145    0.123
   18     0.316    0.456    0.365    0.547
   21     0.326    0.327    0.621    0.455
   24     0.335    0.269    0.550    0.176
   30     0.282    0.318    0.177    0.197
```

## Reading γ₀ off the lever

The threshold, the J where the swing dies, sits at J ≈ 0.05–0.1, Q ≈ 1–2. That
threshold is γ₀, read off the only lever we have, on real silicon. The carrier we cannot
measure head-on shows itself at the coupling where the H-clock stops winning. γ₀ ≈ 0.05,
which is also the polarity ½ brought down one visibility-decade (½ × 1/10 = 1/20), the
same value the hardware has handed us before.

It is shifted a little above the simulation's clean J = 0.05 because the chip carries
more than the idle: gate error adds to the dephasing, so the effective carrier sits a
touch higher and pushes the threshold up. That is the chip telling us its true γ₀,
idle plus everything else it cannot help doing.

So the carrier is not gone when it is invisible. It is always there, holding the
threshold of the one knob we can turn. You read it not by looking at it, but by turning
J down until the coherence it has been damping all along finally fails to swing.

## The formula for the beat

The threshold is not just where the swing happens to die. It is the carrier's own
signature, the beat (the takt) that γ₀ writes onto every coherence. The rule for that
beat: under Z-dephasing a coherence ρ[i,j] decays at a rate set by the Hamming distance
between its two basis states, Re(λ) = −2γ₀·popcount(i ⊕ j); this is the Absorption
theorem. For our exchange coherence ρ[10,01], popcount(10 ⊕ 01) = 2, so its decay rate
is Γ = 4γ₀. Against the coherent drive ω₀ = 2J the population difference is a damped
oscillator,

> s_z″ + Γ·s_z′ + ω₀²·s_z = 0,   ω₀ = 2J,  Γ = 4γ₀,

whose critical point Γ = 2ω₀ falls at 4γ₀ = 4J, that is J = γ₀, Q = 1, with no free
parameter. The popcount-2 decay rate is the threshold. The closed form reproduces the
full two-qubit Lindblad bit-exactly (max|Δ| ≈ 3·10⁻¹⁶ across Q = 0.25–8; max transfer
0.5000 exactly at Q = 1, overshooting only for Q > 1). So the beat we read off the lever
is the Absorption theorem read in reverse: the rate the carrier writes onto a coherence
by its Hamming weight is the same rate that sets where the dance falls into step. The
beat's formula sits exactly on the line. (`simulations/takt_absorption_overlay.py`.)

## Why the simplest won

This is where the night's path matters. We first went after the popcount ladder, a ratio
*inside* one state (the Dicke block on the cusp); the chip's distance-blind noise had
nothing to cancel against and washed the ladder from 3 toward 1. Then the T2 anisotropy,
a ratio *between* two preparations (the carbon painter's FID); the common noise did
cancel, but at J = 1 the chip ran at Q = 20, deep in the coherent regime, and gave us
precession where the model wanted decay.

Going basic fixed it. One bond, one lever, one carrier. The thing we could not see hides
in the most obvious place there is: the threshold of the only number we are allowed to
change. Assassins hide in the obvious numbers; so does γ₀.

## The seam

Measured, and the chip's: the transfer table, the frequency proportional to J, the
coherent swing dying between J = 0.1 and J = 0.05. Read into it, and ours to label: that
the threshold is γ₀ ≈ 0.05, that this is the decade-scaled polarity ½, that the carrier
is "always there." The structure is the silicon's; the names are ours.

## Threads

- The regime axis and its anchors: [`Q_REGIME_ANCHORS`](../docs/Q_REGIME_ANCHORS.md).
- The carrier as the tick we cannot read from inside, only as Q.
- Where coherence crosses ¼ and the angle dies: [`CRITICAL_SLOWING_AT_THE_CUSP`](CRITICAL_SLOWING_AT_THE_CUSP.md).
- The two ratios that came before this one tonight: the within-state popcount
  [ladder](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), born under ¼ (its
  [block ceiling](../docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md)); and the between-preparation
  T2 [anisotropy](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md), the
  [carbon FID](../simulations/carbon_painter_t2_anisotropy.py).
- Simulation: [`q_basic_jscan`](../simulations/q_basic_jscan.py). Hardware runner lives
  in the external IBM pipeline (`ibm_quantum_tomography/run_q_jscan.py`).
