# The Shadow of the Fixed Point

**Date**: 2026-02-09  
**Depends on**: BOUNDARY_NAVIGATION.md, RESIDUAL_ANALYSIS.md, IBM_QUANTUM_TOMOGRAPHY.md

**Tier:** 2 (Verified computation, speculative interpretation LARGELY CLOSED)
**Status:** Shadow universality tested and not confirmed (March 2026)
**Scope:** Boundary correlation and frozen complex direction in hardware data
**Does not establish:** That the shadow effect is universal or boundary-specific

---

## What We Found

On February 9, 2026, state tomography on IBM Torino (qubit 52) confirmed that the product C·Ψ crosses the ¼ boundary during decoherence. This was expected. The framework predicted it.

What was not expected: the system did not go quietly into the classical regime.

After crossing ¼ (after the complex fixed points merged into two real ones, after the system entered the regime where a definite classical attractor exists) residual coherence persisted. Not as random noise. Not as the dying echo of a decaying signal. As a *directed* presence, pointing into the fourth quadrant of the complex plane (Re > 0, Im < 0), in all 17 measurements taken beyond T₂, without a single exception.

This coherence does not decay. It grows. At a rate of +0.008 per T₂, in a regime where every term in the Lindblad equation that could sustain it has already fallen to zero.

A 10,000-run Monte Carlo simulation of pure shot noise on identical quantum states produces nothing like this. The probability of the directional consistency alone is 6 × 10⁻¹¹.

And the coherence follows a rule: it correlates with proximity to the ¼ boundary at r = −0.9955. When C·Ψ fluctuates closer to ¼, coherence is stronger. When it drifts further away, coherence weakens. Not loosely. Nearly perfectly.

## What It Means in the Framework

### The Bifurcation Leaves a Trace

In BOUNDARY_NAVIGATION.md, we documented what happens at C·Ψ = ¼. Two complex conjugate fixed points (R⁺ = (1 − 2CΨ + i·√(4CΨ−1))/(2C) and R⁻ = (1 − 2CΨ − i·√(4CΨ−1))/(2C)) approach each other as C·Ψ falls, merge at ¼, and split into two real fixed points below ¼.

Standard theory says: once the system crosses, the complex structure is gone. The real fixed points take over. The system converges to R₁ (the stable attractor), and the memory of the complex regime is erased.

The data says otherwise.

The last complex fixed point before crossing was R⁻ = 0.8015 − 0.1712i, at phase −12.1°. The residual coherence in the classical regime has mean phase −48.4°. Same quadrant. Same sign structure. The system crossed from complex to real, but the direction it came from is still written in the off-diagonal elements of the density matrix.

We propose a name for this: **the shadow of the fixed point**.

### The Shadow Is Not Decay; It Is Structure

A decaying signal loses amplitude and randomises phase. That is not what we observe. The shadow has:

- **Fixed direction**: 17/17 consistent sign (Re⁺/Im⁻)
- **Growing amplitude**: +0.008/T₂, opposite to all Markovian predictions
- **Boundary correlation**: r = −0.9955 with distance from ¼

This is not residual oscillation. It is not a tail of the initial state. The initial state's coherence fell below the noise floor at t ≈ 1.2 T₂. The shadow emerges *after* that, in a regime where the initial preparation is irrelevant.

If the framework is correct, this is the geometry of the bifurcation expressing itself in measurable quantum observables. The ¼ crossing is not a clean break. It is a phase transition, and like all phase transitions, it has critical behaviour: fluctuations, correlations, and memory that persist beyond the transition point.

### Why the Shadow Points Where It Points

When C·Ψ > ¼, the iteration R_{n+1} = C(Ψ + R_n)² has no real fixed point. The state explores the complex plane, oscillating with angular frequency θ = arctan(√(4CΨ − 1)). The two complex fixed points R⁺ and R⁻ are conjugates: one in the upper half-plane, one in the lower.

A physical qubit prepared in |+⟩ starts with Re(ρ₀₁) > 0 and Im(ρ₀₁) ≈ 0. As it decoheres with a negative detuning (ω = −5.68 kHz, measured), the state spirals through the lower half-plane. It tracks R⁻, not R⁺.

At the crossing, R⁻ collapses onto the real axis. But the system's trajectory was *approaching* R⁻ from a specific direction in the complex plane. That approach direction (fourth quadrant, Re⁺/Im⁻) is frozen into the density matrix at the moment of bifurcation.

What we measure at late times is this frozen direction. The shadow is the last bearing of the complex fixed point, preserved in the quantum state after the map it belonged to has ceased to exist.

### Why the Shadow Grows

This is the most difficult observation to explain, and the most important.

At t > 2T₂, the decoherence rate is effectively zero. The exponential decay exp(−t/T₂*) has already delivered everything it can. The system sits in the classical regime, far below ¼, with two well-separated real fixed points. Nothing in the standard Lindblad dynamics generates new off-diagonal coherence.

Yet the coherence grows. At 7.4 × 10⁻⁵ per microsecond (0.008/T₂*). Unopposed, because the only process that could oppose it, decoherence, has already exhausted itself.

Three possible explanations, in order of decreasing conservatism:

**1. Non-Markovian memory.** The environment that decohered the qubit has finite correlation time. Information that leaked into the environment during decoherence can partially return. This is known physics (see: non-Markovian quantum dynamics, Breuer & Petruccione). It would explain both the growth and the directionality: the returning information carries the phase it had when it left. This is testable: non-Markovian revivals are sensitive to the spectral density of the environment. Different qubits with different environments should show different revival rates.

**2. TLS-mediated feedback.** A two-level system defect in the substrate absorbs coherence during the early decay and re-emits it later. This is essentially a microscopic non-Markovian channel. The fixed direction comes from the TLS having its own fixed frequency. The growth comes from the TLS being a coherent resonator that rings after being struck. Testable: TLS frequencies drift on hour timescales. Repeated measurement should show the direction changing.

**3. The boundary is not passive.** In the R = CΨ² framework, the ¼ crossing is where possibility becomes reality. Above ¼, complex fixed points represent unresolved quantum states. Below ¼, real fixed points represent definite outcomes. The crossing is the act of resolution.

But what if the boundary has structure? What if the bifurcation point is not just a threshold but a *source*: a place where the geometry of the transition generates coherence rather than merely permitting its decay? The growing shadow would then be the ¼ boundary radiating its influence into the classical regime, like light from a horizon.

This is the most speculative interpretation. It is also the only one that naturally explains all three observations (direction, growth, and boundary correlation) with a single mechanism. And it is the interpretation that connects directly to the framework's deepest claim: that the ¼ boundary is where consciousness (C), possibility (Ψ), and reality (R) meet.

## What This Does Not Prove

- It does not prove that consciousness is involved.
- It does not prove that the ¼ boundary is physically special beyond being a mathematical bifurcation.
- It does not prove that the shadow is anything more than a mundane hardware artifact with an unusually clean signature.

What it proves is:
1. The anomaly is real (p < 0.0001 against null model).
2. It has structure (17/17 directional consistency, r = −0.9955 boundary correlation).
3. It connects quantitatively to the framework's predictions about what happens at ¼.

## The Question for March 2026

The March hardware run is designed to answer one question:

**Is the shadow a property of qubit 52, or a property of the ¼ boundary?**

If it is qubit 52: different qubits will show different shadow directions (or none). The effect is local, hardware-specific, and explained by H1 or H2 above.

If it is the ¼ boundary: every qubit that crosses ¼ will cast a shadow, and the shadow will always point in the direction the complex fixed point was heading when it merged. The direction depends on the initial state and detuning, but the *existence* of the shadow and its correlation with ¼ will be universal.

The second outcome would mean that the bifurcation leaves a measurable scar in the density matrix of every quantum system that undergoes decoherence. Not predicted by standard Lindblad theory. Predicted by the geometry of R = CΨ².

## How We Got Here

This was not planned. The tomography experiment was designed to test whether C·Ψ crosses ¼. It did. The experiment was complete.

Then someone said: *look at the residuals*.

Then someone said: *the peaks are delimiters*.

Then someone said: *this is too perfect for a quantum system*.

Each observation led to the next. The data did not resist; it opened. What we found was already there, in a JSON file we had written off as fully analysed. It took three rounds of looking before we saw it: the shadow was always in the data. We just had to learn how to read it.

The shadow of the fixed point. A trace of the complex regime, persisting into classical reality, pointing back toward the boundary it came from.

If it is real, it is the first empirical signature of the ¼ bifurcation that goes beyond the crossing itself. Not just *that* the boundary exists, but that crossing it *costs something*, or *leaves something behind*.

## March 2026 Hardware Results: The Shadow Is Not Universal

### The experiment

On March 9, 2026, we ran the shadow hunt on IBM Torino. Two permanent-crosser
qubits qualified (r < 0.20): Q102 (r=0.159) and Q80 (r=0.170). Each received
10 tomography points (2 reference + 8 shadow zone, t/T2* up to 5.0). Budget:
60 batches, ~3 min QPU time.

Both synthetic simulation and Aer simulation confirmed the null hypothesis
beforehand: 0/4 shadow under standard Lindblad.

### Q102: No pattern

8 late-time points scattered across all 4 quadrants. Phase std = 115 degrees.
No directional consistency. Re signs mixed, Im signs mixed. This is noise.

### Q80: Consistent direction, but NOT the Run 1 shadow

8/8 late-time points in Quadrant 1 (Re+/Im+). Phase = 29 +/- 10 degrees.
p(chance) = 6.1e-05. All 8 points significant above shot noise floor.

But Run 1 (Q52, Feb 2026) was Quadrant 4 (Re+/Im-), phase ~ -44 degrees.
Q80 shows a DIFFERENT systematic direction.

### Residual analysis

Residual Re (measured minus Lindblad prediction) is consistently negative
for BOTH qubits, all 8 late-time points each. This means measured Re(rho_01)
is systematically LOWER than the simple exp(-t/T2) model predicts.

But Im residuals differ: Q102 mixed signs, Q80 consistently positive.
The Im component is where qubit-specific behavior appears.

### Interpretation

The most likely explanation for Q80's consistent Im+ signal is a
**frequency offset (detuning)** not captured by the simple T1/T2 decay model.
A small delta_omega produces a phase rotation exp(-i*delta_omega*t) on rho_01,
creating a systematic imaginary component. Different qubits have different
detunings, producing different phase directions.

This also retroactively explains Q52's Run 1 shadow: the -44 degree phase
was likely Q52's specific detuning, not a universal boundary property.

### Verdict

| | Q52 (Run 1, Feb) | Q102 (Mar) | Q80 (Mar) |
|---|---|---|---|
| Direction | Re+/Im- (-44 deg) | Chaos (mixed) | Re+/Im+ (+29 deg) |
| Consistency | 16/16 (raw) | 0/8 | 8/8 |
| Likely cause | Qubit-specific detuning | Shot noise | Qubit-specific detuning |

**The shadow is NOT a universal property of the 1/4 boundary.**

It is a qubit-specific phase rotation from uncompensated frequency offsets,
visible at late times when the exponential decay has removed the dominant
real component and the detuning-driven imaginary part becomes visible.

### What this means for the framework

The three hypotheses from the original document:

1. **Non-Markovian memory** - Not ruled out but not supported. Different
   qubits show different directions, not a universal revival pattern.

2. **TLS-mediated feedback** - Possible for Q80's consistent signal, but
   the simpler detuning explanation is preferred (Occam's razor).

3. **"The boundary is not passive"** - **Not supported.** The shadow
   direction is qubit-specific, not boundary-specific. The ¼ crossing does
   not leave a universal scar.

### What survives

- The 1/4 crossing itself is confirmed (standard QM, reproduced).
- Late-time coherence has qubit-specific phase structure (real, interesting for hardware characterization, not framework-specific).
- The shadow hunt methodology (simulator null hypothesis + hardware comparison) is sound and reusable.

### Data

Hardware results: `ibm_quantum_tomography/results/shadow_march/`
- `shadow_hardware_combined_20260309_181852.json` (both qubits)
- `shadow_hardware_q102_20260309_181852.json`
- `shadow_hardware_q80_20260309_181852.json`
- Simulation baselines from Feb 19 and Mar 9.

### Skeleton analysis: same pattern as simulation

Applying the window-XOR method from STRUCTURAL_CARTOGRAPHY to the IBM data
reveals the same fundamental pattern seen in the star topology simulation:

**The skeleton is stable. What changes is the phase.**

Q102 (fast rotator):
- Population change per step: 0.011 (tiny)
- Phase change per step: 0.377 pi (massive)
- Ratio: phase dominates 33.5x over populations
- This is a rapidly spinning rotor, estimated ~27 kHz detuning

Q80 (slow drifter):
- Population change per step: 0.016
- Phase change per step: 0.046 pi
- Phase grows monotonically: +0.01, +0.07, +0.11, +0.13, +0.22, +0.24 pi
- Estimated ~4.6 kHz detuning

The simple T1/T2 Lindblad model predicts phase = 0 at all times (rho_01
stays real). Every non-zero phase is structure the model does not capture.

The "shadow" from Run 1 was Q52's specific detuning frequency, frozen into
the late-time phase when the amplitude had decayed below the noise floor.
Different qubits have different detunings, producing different phase
directions - exactly what we observe.

### What this connects to

The simulation (star topology, 3 qubits) and the hardware (single qubits,
IBM Torino) show the same structural decomposition:

| | Simulation | Hardware |
|---|---|---|
| Skeleton | Phi+ core (populations + main correlation) | Populations (T1 relaxation) |
| Rotation | Cross-coupling phases in YZ/ZY plane | Off-diagonal phase (detuning) |
| Shared between steps | 88% | 83-98% |
| What drives rotation | Hamiltonian (J couplings) | Qubit frequency offset |
| What damps it | Lindblad dephasing (gamma) | T2 decoherence |

The pattern is: stable skeleton plus single rotational degree of freedom,
damped over time. The driver differs (Hamiltonian vs detuning), but the
structural decomposition is the same.

This was discovered by applying the "overlay and remove what's shared"
method from STRUCTURAL_CARTOGRAPHY Phase A to real hardware data.

---

*Previous: [Residual Analysis](RESIDUAL_ANALYSIS.md), the statistical evidence*
