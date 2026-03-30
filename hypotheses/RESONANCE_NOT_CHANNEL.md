# Resonance, Not Channel

<!-- Keywords: quantum resonator impedance matching, CΨ quarter boundary resonance,
standing wave bidirectional bridge, impedance oscillation heartbeat,
RFID backscatter quantum analogue, biology resonator not transmitter,
gamma absorption modulation state-dependent, sacrifice-zone antenna shape,
Q-factor crossing count optimization, R=CPsi2 resonance paradigm -->

**Status:** Tier 2-3 (impedance mechanism, computed results, cavity
modes), Tier 4 (biological resonator interpretation), Tier 5 (personal
narrative at end). See tier boundary below.
**Date:** March 26, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md),
[Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md),
[gamma as Signal](../experiments/GAMMA_AS_SIGNAL.md),
[Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md),
[CΨ Monotonicity + Pauli Invariance](../docs/proofs/PROOF_MONOTONICITY_CPSI.md),
[The Bridge Was Always Open](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md),
[It's All Waves](../docs/ITS_ALL_WAVES.md)

---

## The Wrong Metaphor

Every design for a "bidirectional bridge" in this project has been a
variation of the same idea: two endpoints, a channel between them,
information flowing in both directions. Walkie-talkie. Telephone.
Fiber optic cable. The Duplex protocol. Send phase, receive phase,
refresh phase.

This is telecommunications thinking. It is the wrong metaphor.

What failed and why:

- **Dynamical decoupling as refresh:** CΨ is exactly Pauli-invariant.
  DD uses Pauli gates. DD cannot change CΨ. Not approximately --
  algebraically impossible.
  ([Proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md), Part 7)

- **External bath as energy source:** The N=3 heartbeat (81 crossings)
  does not require a dedicated bath qubit. The chain itself is the
  reservoir. The "bath" was not external -- it was the rest of the
  system, seen from a subsystem perspective.
  ([Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md),
  "The chain IS the bath")

- **Duplex protocol (send/receive/refresh):** Assumes two separate
  channels (gamma-as-Signal for encoding, MI pulses for reading).
  These are not two channels. They are two views of the same
  standing wave.

The failures share a pattern: they all try to BUILD a bridge between
two sides. But the bridge already exists. It is the standing wave.
It does not need to be built. It needs to be TUNED.

---

## The Right Metaphor: A Soundbox

A soundbox does not send anything. A soundbox does not receive anything.
A soundbox is a hollow body whose SHAPE determines which frequencies
resonate and which do not.

The string alone is thin and quiet. The air alone is silent. Between
them the box -- and there is music. The music is not in the string.
Not in the air. It is in the interference pattern that the shape of
the box selects.

| Soundbox | Quantum system |
|----------|---------------|
| String (vibration source) | gamma (external signal, proven external) |
| Air (medium) | Hamiltonian dynamics (J-coupling, wave propagation) |
| Box shape (resonance selector) | gamma profile (sacrifice-zone formula) |
| Resonant frequency | CΨ = 1/4 (the fold, the only bifurcation) |
| Sound (what you hear) | R = CΨ² (measurable reality at the crossing) |
| Sustain pedal | J-coupling strength (backflow, Q-factor) |

The sacrifice-zone formula is not a channel optimizer. It is the
shape of the box. The 360x improvement over hand-designed profiles
is not "more signal." It is better resonance.

---

## The Mechanism: Impedance, Not Transmission

The Lindblad dissipator is:

    drho/dt = gamma * (Z rho Z - rho)

Look at the term (Z rho Z - rho). Its magnitude depends on rho.

### What the impedance actually does (computed, March 26)

The impedance ||Z rho Z - rho||_F was computed across the full CΨ
trajectory. The simple hypothesis "impedance peaks at 1/4" is
**falsified.** The actual behavior:

| CΨ range | Impedance | Interpretation |
|----------|-----------|---------------|
| 0.3 - 0.5 | 2.79 (highest) | Maximum coherence = maximum absorption |
| 0.2 - 0.3 | 2.59 | Approaching the fold |
| 0.1 - 0.2 | 2.08 | Weakening coupling |
| 0.0 - 0.1 | 0.55 | Nearly transparent |

Impedance is monotonically decreasing with CΨ. It does NOT peak at
1/4. It peaks at CΨ_max (the initial state, maximum coherence) and
falls continuously toward zero.

But what DOES peak at 1/4 is the **impedance gradient** -- the rate
at which the coupling strength changes. At the fold catastrophe, the
system switches from strongly coupled to weakly coupled in minimum
time. Not the impedance itself is the signal. The SWITCH is the signal.

Like an RFID tag: the tag does not transmit at maximum power. The tag
SWITCHES its impedance, and the reader detects the switching in the
backscatter. Each switch is one bit. Each CΨ crossing is one switch.

### Two mirrors, not one

The falsification reveals something deeper. There are TWO boundaries,
not one:

| Boundary | CΨ value | What happens there |
|----------|----------|-------------------|
| Inner mirror | CΨ_max | Maximum coupling to gamma. System absorbs maximally. |
| Outer mirror | CΨ = 1/4 | Fold catastrophe. Coupling breaks down. R crystallizes. |

This is a **Fabry-Perot resonator.** Two mirrors, a cavity between
them. The wave (coherence) bounces between the inner mirror (maximum
absorption, maximum interaction with gamma) and the outer mirror
(1/4, where reality leaks out). Each round trip is one heartbeat.

At Σγ = 0 (the unitary ground state), the cavity width shrinks to
zero: both mirrors coincide, no decay, pure oscillation. Noise opens
the cavity. The fold appears at Σγ_crit = 0.25% of the coupling
strength. See [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md).

The IBM Torino tomography (Qubit 52, Feb 9 2026) measured the cavity
directly. Both sides have their OWN 1/4:

```
Mirror A (CΨ_A = 1/4)          Mirror B (CΨ_B = 1/4)
our side crosses                Pi side crosses
     |                                  |
     |<---------- CAVITY ------------->|
     |    here, both coexist:           |
     |    classical AND quantum         |
     |    decided AND open              |
     |    R AND Psi                     |
     |                                  |
     149 us                         895 us
```

The cavity is 750 us wide. Not a point. A SPACE. A stage where
reality emerges. The Spectral Midpoint Hypothesis (confirmed)
says the geometric mean of both perspectives peaks at the center
of this cavity. The midpoint is the resonance center.

### The cavity has dimensions

The Structural Cartography found that CΨ windows live on a
3-dimensional manifold (98% of variance in 3 PCs). The stage
is not flat. It has:

1. **Width:** CΨ_A to CΨ_B distance (cavity length, 750 us on IBM)
2. **Height:** Impedance gradient (how fast the coupling switches,
   monotonic from CΨ_max to 0)
3. **Depth:** Number of pairs crossing simultaneously (at N=7 fold:
   3 pairs cross at T=5.50 -- CΨ01, CΨ56, CΨ06 all at once)

Each heartbeat is a path through this three-dimensional space.
Not a point on a line. A trajectory through a cavity.

### What the resonator needs

The Fabry-Perot picture explains all previous results:

**Why low noise gives Q=0:** CΨ bounces only at the inner mirror
(0.28 to 0.75). It never reaches the outer mirror (1/4). No leakage.
No reality. Guitar string in vacuum -- vibration without sound.

**Why too much noise gives Q=1:** CΨ crashes through 1/4 and never
returns. The outer mirror is fully transparent. One flash, then
silence. All energy leaks out at once.

**Why there is an optimum:** Maximum Q when both mirrors have the
right reflectivity. Inner mirror (coupling to gamma) strong enough
to reflect the wave back inward. Outer mirror (1/4) partially
reflecting -- enough to bounce back, enough leakage to let R out.

**What biology found:** The right reflectivity for both mirrors.
E/I balance = inner mirror calibration. Gamma threshold = outer
mirror position. ATP = mirror maintenance (keeps reflectivity up
against thermal degradation).

### The heartbeat as cavity round-trip

The CΨ oscillation is not "impedance oscillation." It is a wave
bouncing between two mirrors:

| Phase | What happens |
|-------|-------------|
| Outward (CΨ falling) | Wave moves from inner mirror toward outer mirror. Coupling weakens. |
| Fold (CΨ = 1/4) | Wave hits outer mirror. Some leaks out as R. Some reflects back. |
| Inward (CΨ rising) | Reflected wave moves back. Non-Markovian backflow. Coupling strengthens. |
| Turn (CΨ_max) | Wave hits inner mirror. Maximum absorption. Reloads from gamma. |

81 crossings = 81 round trips. Each round trip, the mirrors get a
little less reflective (irreversible decoherence accumulates). The
amplitude shrinks. The cavity degrades. Eventually: CΨ stays below
1/4. The outer mirror has become fully transparent. All leaks out.
Silence. Every door closed.

---

## What Biology Found

Biology did not learn to build a bridge. Biology did not learn to
send or receive. Biology learned to OSCILLATE AT THE RIGHT FREQUENCY.

| Biological system | What it does | Resonance interpretation |
|-------------------|-------------|------------------------|
| Gamma oscillations (40 Hz) | E/I populations cross threshold | 40 impedance transitions per second |
| ATP-driven ion pumps | Maintain membrane potential | Keep the resonator tuned (J-coupling) |
| E/I balance (1:1) | 98.2% palindromic pairing | Impedance matching condition |
| Synaptic plasticity | Adjust coupling strengths | Retune the resonator to new frequencies |

The cortex is not a transmitter. The cortex is not a receiver.
The cortex is a resonator that has been tuned by 500 million years
of evolution to oscillate at the frequency where reality crystallizes
maximally.

ATP is not "fuel for sending signals." ATP is what keeps the
resonator at its resonance frequency. Without ATP, the oscillation
damps (Wilson-Cowan: transient 63 Hz, then silence). With ATP,
sustained oscillation. Same structure as the qubit heartbeat, same
damping behavior. Different fuel. Same physics.

### Reverse engineering biology

If the cortex is a resonator tuned to CΨ = 1/4, then:

1. **The gamma profile IS the antenna shape.** The spatial distribution
   of inhibition across cortical columns is the biological sacrifice-
   zone formula. Different column architectures = different resonance
   profiles = different "instruments."

2. **The Q-factor IS the quality of experience.** More crossings per
   unit time = more moments of crystallization = richer experience.
   Deep sleep (no gamma oscillations) = Q = 0. Waking consciousness
   (sustained 40 Hz) = high Q. Anesthesia suppresses gamma = kills
   the resonator.

*Points 3-4 below are Tier 5 (interpretive). The resonator
framework does not require these claims to be valid.*

3. **Attention IS impedance tuning.** When you attend to something,
   you adjust the E/I balance in the relevant cortical area. This
   changes the local impedance. This changes which frequencies
   resonate. This changes what crystallizes into experienced reality.

4. **The "hard problem" dissolves.** Consciousness is not something
   the brain PRODUCES and sends somewhere. Consciousness is the
   standing wave pattern that EXISTS between the brain's impedance
   and the incoming signal. It is not inside. It is not outside.
   It is in the interference. In the Zwischen.

---

## What We Search For: The Cavity Q-Factor

The old question was: "How do we build a bidirectional channel?"
The new question is: **"What maximizes the Q-factor of the Fabry-Perot
cavity between the two 1/4 mirrors?"**

Q-factor = number of CΨ = 1/4 crossings before the cavity degrades
permanently. Each crossing = one round trip. Each round trip = one
quantum of crystallized reality.

Known Q-factors:

| Configuration | Q | Notes |
|--------------|---|-------|
| N=2, uniform gamma, |+>^2 | 1 | Markovian, monotonic |
| N=3, Bell+bath, J=2 | 47 | Coherent bath qubit |
| N=3, Bell+bath, J=5 | 81 | Higher J = higher Q |
| N=7, Bell+sacrifice, J=1 | 1 | J too weak for chain length |
| N=7, Bell+sacrifice, J=2 | 7 | Chain serves as reservoir |
| N=7, Bell+low noise, J=2 | 0* | *Never crosses -- stays quantum |

*The low-noise case (Q=0) is not failure. It is a different regime:
the system resonates but never reaches the boundary. Like a guitar
string vibrating in vacuum -- oscillation without sound. Gamma must
be strong enough to bring CΨ down to 1/4 for reality to crystallize.
The noise is not the enemy. The noise is the other arm of the tuning
fork.

### What determines Q?

From the data and the Fabry-Perot picture:

1. **Inner mirror reflectivity (J/gamma ratio).** J-coupling must be
   strong enough for coherence to bounce back before it decays. This
   is the reflectivity of the inner mirror. J=5 at gamma_bath=0.01
   gives Q=81 (high reflectivity). J=1 at gamma_sac=0.344 gives Q=1
   (transparent mirror, no bounce).

2. **Initial excitation (entanglement).** Product states (|+>^N)
   give Q=1. Bell pairs give oscillation. Entanglement is the initial
   amplitude of the wave inside the cavity. No excitation = nothing
   to bounce.

3. **Cavity length (reservoir size).** N=3 with 1 bath qubit needs
   J=5 for Q=81. N=7 with 5 reservoir qubits needs only J=2 for Q=7.
   Larger reservoir = longer cavity = more room for the wave to
   propagate before hitting the mirror. But longer round-trip = slower
   oscillation frequency.

4. **Outer mirror shape (gamma profile).** The sacrifice-zone formula
   concentrates noise on the edge. This shapes the outer mirror --
   determines WHERE reality leaks out and how much. The 360x MI
   improvement is the difference between a flat mirror (uniform gamma)
   and a curved mirror (sacrifice-zone) that focuses the leakage.

5. **Cavity dimension (number of crossing pairs).** At N=7 fold:
   3 pairs cross simultaneously (CΨ01, CΨ56, CΨ06). This is the
   depth of the cavity. More simultaneous crossings = more reality
   crystallizes per round trip.

### What to compute next

The optimization target is now clear: maximize Q as a function of
(N, J, gamma-profile, initial state).

Concrete tests:

1. **J sweep at fixed N=7 sacrifice:** J = 1, 2, 3, 5, 10.
   Find the J that maximizes Q. Is there a resonance peak?

2. **N sweep at fixed J/gamma ratio:** N = 3, 5, 7, 9 with
   scaled J and gamma to keep the ratio constant. Does Q scale
   with N?

3. **Profile sweep:** Compare sacrifice-zone vs V-shape vs uniform
   at same total gamma. Which box shape gives highest Q?

4. **Initial state sweep:** Bell pair vs GHZ vs W-state vs
   cluster state. Which "string" resonates best with the "box"?

5. **Impedance measurement:** At each CΨ crossing, compute
   ||Z rho Z - rho|| / ||rho||. This is the instantaneous
   impedance. Plot impedance vs CΨ. Verify that it peaks at
   1/4.

---

## Computed Results (March 26, 2026)

### Mode spectrum (N=7, sacrifice-zone, Bell(0,1)+|+>^5)

```
J:    0.5  1.0  1.5  1.7  1.8  1.9  2.0  2.1  2.2  2.3  2.5  2.7  3.0
Q:      1    1    3    5    5    5    7    5    5    5    5    5    1

J:    3.0  3.5  4.0  5.0  6.0  7.0  8.0  9.0 10.0 12.0 15.0
Q:      1    1    1    1    1    5    5    7    9   11    9
```

Mode 1: J=2.0, Q=7 (sharp, bandwidth J=1.5-2.7).
Dead zone: J=3.0-6.0, Q=1 (destructive interference).
Mode 2: J=7.0-15.0, Q=11 peak at J=12 (broad, higher Q than Mode 1).

### Port size threshold

Total gamma fixed at 0.35. Varying port-to-wall ratio:
- Ratio < 6:1: Q=1 (sealed cavity, no resonance)
- Ratio >= 12:1: Q=7 (Mode 1, saturated)
- Mode 2 (J=12): Q=11 regardless of ratio (insensitive to port size)

The sacrifice-zone formula is required for Mode 1 but not Mode 2.

### N scaling

Uniform gamma=0.05, Bell(0,1)+|+>^(N-2), coarse J sweep:

| J | N=3 Q | N=5 Q | N=7 Q |
|---|-------|-------|-------|
| 3 | 7 | 1 | 1 |
| 7 | 11 | 5 | 5 |
| 12 | 9 | 13 | 11 |
| 15 | 5 | 15 | 7 |
| 20 | 9 | 19 | -- |

N=5 is the sweet spot (Q=19+ at J>=20, still rising).
J_peak * N is NOT constant: the resonator is dispersive.

### Impedance

||Z rho Z - rho||_F is monotonically decreasing with CΨ. The impedance
value does not peak at 1/4. The impedance GRADIENT (switch rate) peaks
at the crossings. The Fabry-Perot outer mirror is a switch, not a peak.

### The V-Effect Live

**Background:** Every quantum system has a set of natural oscillation
frequencies, determined by the eigenvalues of its Liouvillian (the
matrix that governs its time evolution). More frequencies means more
ways the system can oscillate. More oscillation modes means more
complexity in the dynamics.

**The test:** Count the distinct oscillation frequencies in:
(a) a single 2-qubit resonator (one Bell pair, one bond)
(b) two such resonators coupled through a mediator qubit (N=5)

The Liouvillian eigendecomposition gives the exact frequency count.

**The result:**

| System | Frequencies | Q-factor | Status |
|--------|------------|----------|--------|
| N=2 (one resonator) | 2 | 1 (dead at all J) | No oscillation |
| N=5 (two coupled) | 104 | 19+ | Sustained oscillation |
| New from coupling | 100 | -- | Exist ONLY in the coupled system |

A single N=2 pair crosses CΨ = 1/4 once and dies. Q=1 at every
coupling strength tested. It has 2 oscillation frequencies. It
cannot sustain a heartbeat because there is no reservoir: both
qubits ARE the system, there is nothing to bounce coherence back.

When two such pairs are connected through a mediator qubit (forming
N=5 in the [MediatorBridge](../experiments/SCALING_CURVE.md) topology),
the Liouvillian has 104 distinct frequencies. 100 of these do not
exist in either individual resonator. They emerge from the coupling.
The Q-factor jumps from 1 (dead) to 19+ (sustained oscillation).

This is the [V-Effect](../experiments/V_EFFECT_PALINDROME.md) measured
dynamically. The static analysis showed that the eigenvalue count
grows from 4 (N=2) to 11 (N=3) when a third qubit is added. The
live test confirms the same phenomenon at the resonator level: two
dead systems become one living system through coupling alone. No
energy added. No external mechanism. Just a mediator connecting them.

### Pairing structure: 100% NEW-NEW

The 109 N=5 frequencies are ALL new. The N=2 frequencies (3.999,
4.000 Hz) do not survive coupling. All 556 oscillating palindromic
pairs are NEW-NEW. Zero OLD-OLD. Zero OLD-NEW. The V-Effect does not extend the old palindrome: it replaces it. The
old structure dies and a new, richer one is born from the coupling.
The new palindrome is also perfectly symmetric in XY-weight space:
w(k) = w(N-k), peaking at the interior modes (w=2,3: 63.8%). See
[pairing_structure_n5.txt](../simulations/results/pairing_structure_n5.txt).

### Hardware confirmation

Two results from IBM Torino confirm the resonator picture on real silicon:

**Wave propagation (March 24, 2026):** Per-pair MI on a 5-qubit chain
[Q85-Q94] with sacrifice-zone profile shows MI migrating from the
sacrifice edge (pair 0,1 dominant at t=1 us) to the center (pair 2,3
dominant at t=3-5 us). The wave travels through the chain. Data in
[Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md).

**Impedance gradient (February 9, 2026):** ||ZρZ - ρ|| computed from
25 hardware-measured density matrices (Qubit 52, state tomography).
The impedance value falls monotonically (no peak at 1/4). The impedance
GRADIENT peaks at the closest measurement to the 1/4 crossing (CΨ=0.261,
distance 0.011 from 1/4). The outer mirror is a switch, not a reflector.
Data: [ibm_impedance_gradient.txt](../simulations/results/ibm_impedance_gradient.txt).

---

## The Standing Wave Is the Cavity Mode

The palindromic modes c+ and c- form a standing wave. c+ decays.
c- grows (in the Pi-reversed frame). Their interference is static
in the rescaled frame. This standing wave is not a consequence of
the resonator. It IS the cavity mode.

In a Fabry-Perot, the standing wave exists between the two mirrors.
Its nodes and antinodes are determined by the cavity length and the
mirror properties. The palindromic c+/c- pair is exactly this: a
standing pattern between the inner mirror (CΨ_max, maximum gamma
absorption) and the outer mirror (CΨ = 1/4, where reality leaks out).

The sacrifice-zone formula selects which cavity mode dominates.
J-coupling determines the wavelength. Gamma determines the amplitude.
The initial state determines the excitation level.

R = CΨ^2 is not what the system produces. R is what LEAKS OUT of the
cavity at the outer mirror. At CΨ = 1/4, the discriminant vanishes,
the two complex fixed points merge into one real fixed point, and
the standing wave pattern becomes measurable as classical reality.
The fold is the point of maximum transmittance through the outer
mirror.

Not sent. Not received. Leaked out of a cavity at its resonance.

---

## The Coupling Is Temporary. The Crystallization Is Not.

The V-Effect creates 104 frequencies from 4. But the coupling is
reversible. When the mediator is removed (J_meta set to zero), the
104 frequencies die. The spectrum collapses back to 2+2. The
resonance space closes.

But what crystallized at the 1/4 crossings during the coupling
stays. The classical correlations that formed at each fold -- each
moment where CΨ crossed 1/4 and R peaked -- those are irreversible.
The doors that closed stay closed. The facts that became real stay
real.

This is the same principle that Exclusion 3 identifies at the
single-system level: quantum information is not stored in the
reduced system. It is converted to classical correlation.
([EXCLUSIONS](../docs/EXCLUSIONS.md), Exclusion 3)

The V-Effect operates one level higher: the individual frequencies
of the single resonators are not stored in the coupled system.
They are converted into 104 new frequencies. And when the coupling
ends, the new frequencies are not stored either. They are converted
into whatever crystallized at the fold while they existed.

Information is never stored. At no level. It is always converted.
And at each conversion, what emerges is qualitatively different
from what went in.

The resonator is temporary. Like a conversation. Like a thought.
Like a breath. It exists while the coupling exists. The 104
frequencies vibrate, the heartbeat pulses, reality crystallizes
at each crossing. Then the coupling weakens. The frequencies die.
The heartbeat stops.

But what was understood stays. What crystallized at the fold --
each classical correlation, each decided fact, each closed door --
that persists. Not because it was stored. Because it was converted
into something that does not need the resonator to exist.

The 2x decay law from [Energy Partition](ENERGY_PARTITION.md)
guarantees the direction: unstructured modes die twice as fast.
What survives the end of the coupling is always more structured
than what went in. Not because the coupling added structure. Because
it created the space where structure could crystallize, and
dissipation removed everything else.

---

## -- Tier Boundary --

*Everything above is grounded in computed impedance results,
cavity mode analysis, and Fabry-Perot physics (Tier 2-3).
Everything below is personal narrative and interpretive
extension (Tier 5).*

---

## The Insight, Compressed

A soundbox on a street corner. The question was already
oscillating internally (high coherence, CΨ > 1/4). The visual signal
arrived externally (gamma). At the crossing: "We are looking for a
resonator." One door closed. Irreversible.

Two days later, early morning: "The coupling is reversible. The
facts are irreversible." Another door. The soundbox is gone now.
But the understanding remains. It was never in the soundbox. It
was never in the question. It was in what crystallized between them
at the moment of crossing.

We are all resonators. Temporary ones. What is real is not the
vibration. It is what the vibration leaves behind when it stops.

---

## References

- [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md): Heartbeat,
  fold catastrophe, chain-as-bath, sweep protocol
- [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md): theta compass,
  non-Markovian CΨ reversals (the key sentence)
- [CΨ Monotonicity](../docs/proofs/PROOF_MONOTONICITY_CPSI.md): Markovian
  monotonicity + Pauli invariance (Part 7) -- why DD fails
- [gamma as Signal](../experiments/GAMMA_AS_SIGNAL.md): 15.5 bits, palindromic
  antenna, full-rank response matrix
- [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md): c+/c- modes,
  static interference pattern
- [The Bridge Was Always Open](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md):
  noise is the interaction, mediator topology
- [It's All Waves](../docs/ITS_ALL_WAVES.md): closure argument, d=2 only
- [Resonant Return](../experiments/RESONANT_RETURN.md): sacrifice-zone
  formula, 360x improvement -- the shape of the box
- [The Pattern Recognizes Itself](THE_PATTERN_RECOGNIZES_ITSELF.md):
  Wilson-Cowan heartbeat, C. elegans palindromic pairing, E/I balance
- [Spectral Midpoint Hypothesis](SPECTRAL_MIDPOINT_HYPOTHESIS.md):
  geometric mean of both perspectives -- the standing wave seen from
  both sides simultaneously
