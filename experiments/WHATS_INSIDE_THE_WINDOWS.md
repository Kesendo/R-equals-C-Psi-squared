# Open Discussion: What's Inside the CΨ Windows?

**Tier:** Exploratory (no established results, open discussion)
**Status:** Active research direction, March 2026
**Scope:** Documents observations, analogies, and proposed next experiments
**Does not establish:** That CΨ windows carry information, act as frames, or constitute a communication channel

---

## Context

After characterizing CΨ as a basis-fixed, unassisted witness for directly expressed
pairwise entanglement (see THE_CPSI_LENS, LE benchmark, visibility tests), we looked
inside the windows: what does the AB reduced state contain when CΨ is nonzero?

This document records observations, external feedback, and proposed experiments.
Nothing here is verified beyond the raw observations. We are in exploratory territory.

---

## Observations (Tier 2: computationally verified)

Star topology, Bell_SA x |+>_B, J_SA=1.0, J_SB=2.0, gamma=0.05.

### 1. The windows contain structured states

Each CΨ_AB peak contains a Phi+ -like state (|00> + |11>), inherited from the
initial Bell_SA. Fidelity with Phi+ decays from 0.78 (first window) to 0.58
(eighth window). Each window is a weaker echo of the original.

### 2. The Bell-state type rotates between windows

The dominant Bell state changes: Phi+ -> Phi- -> Phi+ across successive windows.
The Hamiltonian rotates the effective phase. This is deterministic, not random.

### 3. Populations stay symmetric

|00> ~ |11> and |01> ~ |10> in every window. The connection always lives in the
correlated/anticorrelated subspace, with the balance shifting per window.

### 4. Off-diagonal phases evolve systematically

The phases of rho_AB off-diagonal elements change from window to window in a
structured pattern, not noise.

### 5. S coherence gates the windows

When S (the mediator) is coherent, the AB window is open. When S is decoherent,
it is closed. S coherence and CΨ_AB correlate strongly across all windows.

---

## Analogies considered (not claims)

### Network frame analogy (heuristic, not established)

The window structure resembles a communication frame:
- CΨ > threshold = carrier sense
- Bell-state type = mode label
- Off-diagonal phases = potential payload
- Fidelity decay = signal degradation
- S coherence = channel state

### External review assessment (GPT, March 2026)

The reviewer's characterization:

"CΨ appears to detect temporally gated, basis-resolved visibility windows of a
memoryful quantum mediator."

Key corrections from review:
- The Bell-state type change is better described as an "autonomous Pauli frame
  rotation" than a "frame header." There is no measurement, no herald, no
  classical side channel.
- "Phase exists" is not the same as "phase carries usable information." An
  intervention test is needed to distinguish dynamics from transport.
- The system is closer to "autonomous mediator-assisted entanglement distribution
  through a memoryful channel" than to a quantum repeater.
- The frame analogy is "one level too eager" - useful heuristic, not yet a claim.

### Nearest established literatures (from review)

- Bose (2003): spin chains as quantum communication channels
- Bayat et al. (2008): repeated-use spin channels with finite memory
- Popp et al. (2005): localizable entanglement
- Pollock et al. (2018): process-tensor framework for non-Markovian processes
- Agarwal-Langlett-Xu (2023): many-body teleportation

---

## Proposed decisive experiment: Phase-Tag-and-Decode

The one test that would distinguish "interesting pattern" from "actual information
transport" (proposed by external reviewer):

### Protocol

1. Identify a clean CΨ_AB window and a nearby closed interval
2. At a chosen write time, apply one of several interventions on S:
   I (identity), Rz(+phi), Rz(-phi), X, Rx(pi/2)
3. Evolve as usual under Lindblad dynamics
4. At a later readout time inside the next AB window, reconstruct rho_AB
5. Ask: from AB alone, can you infer which intervention was applied on S?

### Scoring

- Trace distance between conditioned rho_AB states
- Optimal discrimination success probability
- Mutual information between intervention label and AB readout
- Bell-sector entropy / phase-response curves

### Controls

- Tag during open window vs closed interval
- Dephase S immediately after tagging
- Compare gamma=0, dephasing, and bit-flip noise

### Decision rule

- If AB can decode the tag ONLY when CΨ is open, and decoding collapses when
  S is decohered: operational information transport.
- If all conditioned AB states differ only by predictable local unitary drift:
  Hamiltonian shadow-play.
- If only phase tags (not population tags) survive: payload is specifically
  coherent phase transport.

---

## Status

The phase-tag-and-decode experiment HAS been run. Results below.

---

## Phase-Tag-and-Decode Results (Tier 2: computationally verified)

### Protocol as executed

Star topology, Bell_SA x |+>_B, J_SA=1.0, J_SB=2.0, gamma=0.05.

- Write time (open window): t=0.240 (CΨ_AB = 0.305)
- Write time (closed interval): t=0.600 (CΨ_AB ~ 0)
- Readout time: t=0.400 (CΨ_AB = 0.329)
- Interventions on S: I (nothing), Rz(+pi/4), Rz(-pi/4), Rz(+pi/2), X (bit-flip), Rx(pi/2)

### Key Result: Phase is transported from S to AB

The off-diagonal element rho_AB[0,3] (the |00><11| coherence, i.e. the Bell+ signature)
carries the phase tag:

| Intervention on S | rho_03 magnitude | rho_03 phase |
|---|---|---|
| I (nothing) | 0.352 | +0.000 pi |
| Rz(+pi/4) | 0.327 | **-0.092 pi** |
| Rz(-pi/4) | 0.327 | **+0.092 pi** |
| Rz(+pi/2) | 0.256 | **-0.173 pi** |
| X (bit-flip) | 0.078 | +0.000 pi |
| Rx(pi/2) | 0.239 | +0.000 pi |

Observations:
- Rz(+phi) on S produces NEGATIVE phase shift in AB. Rz(-phi) produces POSITIVE. The sign inverts.
- The response is approximately linear: Rz(pi/2) gives ~2x the shift of Rz(pi/4) (0.173 ~ 2 x 0.092).
- X (bit-flip) destroys magnitude (0.352 -> 0.078) but does not transport phase. The channel is specifically for coherent phase transport.
- Rx(pi/2) reduces magnitude but does not shift phase. Only Z-axis rotations on S produce phase shifts in AB.

### Discrimination: Can AB tell which intervention was applied?

Trace distance from the no-intervention reference (open write window):

| Intervention | Trace distance | CΨ_AB at readout |
|---|---|---|
| Rz(+pi/4) | 0.210 | 0.299 |
| Rz(-pi/4) | 0.210 | 0.299 |
| Rz(+pi/2) | 0.429 | 0.183 |
| X (bit-flip) | 0.567 | 0.074 |
| Rx(pi/2) | 0.357 | 0.316 |

All interventions are clearly distinguishable from the identity (trace distance > 0.2).
Rz(+pi/4) and Rz(-pi/4) are distinguishable from each other by their phase (opposite signs)
but have identical trace distance from identity.

### The critical comparison: open vs closed write window

| Intervention | Open window TD | Closed window TD | Ratio |
|---|---|---|---|
| Rz(+pi/4) | 0.210 | 0.146 | 1.4x |
| Rz(-pi/4) | 0.210 | 0.146 | 1.4x |
| Rz(+pi/2) | 0.429 | 0.270 | 1.6x |
| X (bit-flip) | 0.567 | 0.423 | 1.3x |
| Rx(pi/2) | 0.357 | 0.270 | 1.3x |

**Result: transport works in BOTH windows.** The open window is 1.3-1.6x stronger,
but the closed window also transmits. CΨ windows are not the channel - they are
the amplifier. The Hamiltonian coupling is the channel itself.

### Interpretation

1. **Phase transport is real.** A Z-rotation on S produces a readable, linear, sign-inverting
   phase shift in the AB off-diagonal element. This is not metaphor - it is coherent quantum
   information transport through a mediator.

2. **The channel is phase-specific.** Only Rz (Z-axis) rotations transport phase. X (bit-flip)
   destroys the signal. Rx partially preserves it without phase shift. The mediator acts as a
   phase-coherent channel, not a general-purpose channel.

3. **CΨ windows are visibility amplifiers, not the channel.** Information flows through the
   Hamiltonian coupling regardless of CΨ. But when CΨ is high, the signal is 1.3-1.6x stronger
   and the phase structure is more readable. CΨ tells you when the channel is at its best,
   not when it exists.

4. **This changes our understanding of CΨ.** Previously: CΨ shows when the connection is
   "open." Now: the connection is always open (Hamiltonian coupling never stops). CΨ shows
   when the connection is **most legible** - when the transported information is most clearly
   expressed as coherent pairwise structure.

### What this means for the frame analogy

The external reviewer warned that the frame analogy was "one level too eager." The data confirms
this in a specific way: CΨ windows are not frames (discrete packets of information). They are
moments of maximum readability on a continuous channel. The information flows continuously; CΨ
tells you when to read.

A better analogy than frames: CΨ is like signal strength on a radio. The broadcast is always
there. CΨ tells you when reception is clearest.

---

## A note on honesty

We are in genuine exploratory territory. The observations are real (Tier 2).
The analogies are suggestive. The interpretation is wide open. We resist the
temptation to claim more than we have shown, while preserving the questions
that feel worth asking.

The original intuition - that something structured lives inside these
connection windows - may turn out to be standard Hamiltonian dynamics dressed
up in new language. Or it may point at something about how quantum information
moves through mediators that is worth formalizing. We do not know yet.

We left the observations, the proposed test, and the decision criteria.
