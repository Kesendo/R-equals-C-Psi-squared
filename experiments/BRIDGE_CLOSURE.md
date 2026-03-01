# Bridge Closure: Pre-Shared Entanglement = Shared Randomness

**Date**: 2026-03-01
**Status**: Closed (Tier 2). The bridge protocol cannot work. This is not
a limitation of the framework — it is a consequence of quantum mechanics.
**Authors**: Thomas Wicht, with Claude (Anthropic)
**Depends on**: NO_SIGNALLING_BOUNDARY.md, BRIDGE_FINGERPRINTS.md,
BRIDGE_PROTOCOL.md (hypothesis), QKD_EAVESDROPPING_FORENSICS.md

---

## 1. What This Document Settles

NO_SIGNALLING_BOUNDARY.md showed that CΨ drops from 0.500 to 0.250
when B measures, but A cannot see it. This left one open question:

> Can pre-encoded CΨ fingerprints carry something that a classical
> pre-shared key cannot?

The answer is **no**. This document proves it and closes the bridge
hypothesis permanently.

---

## 2. The Argument (Three Lines)

After separation with zero coupling (J = 0):

1. A's measurement statistics: P(a) = Tr[M_a · ρ_A]
2. ρ_A is independent of anything B does (no-signalling)
3. Therefore A's output is a function of {ρ_A(0), E_A} only

ρ_A(0) is determined at preparation. E_A is A's local environment.
Both are available to A without any quantum resource. A classical
schedule that says "pair 7 was prepared as |+0⟩" contains strictly
MORE information than A's qubit, because the qubit loses B's part:

| Preparation | ρ_A | A can distinguish? |
|-------------|-----|-------------------|
| Bell+ | I/2 | From other Bell states: NO |
| |++⟩ | |+⟩⟨+| | From |+0⟩ or |+−⟩: NO |
| |+0⟩ | |+⟩⟨+| | From |++⟩ or |+−⟩: NO |
| |+−⟩ | |+⟩⟨+| | From |++⟩ or |+0⟩: NO |
| |00⟩ | |0⟩⟨0| | From Bell+: YES |

A's qubit carries at most 1 bit (which rho_A was prepared).
The schedule carries log₂(N_states) bits per pair.
The qubit is strictly inferior to the schedule.

For Bell+ pairs specifically: ρ_A = I/2 at all times. The qubit
is maximally mixed noise. Zero bits. Forever. Under any local
channel. Regardless of what B does, has done, or will do.

---

## 3. Numerical Verification

### 3.1 B's Action Is Invisible

Bell+ evolved to t = 2 under independent dephasing (γ = 0.05, J = 0).
Three scenarios: B does nothing, B measures Z, B measures X.

| t_after | ρ_A(nothing) | ρ_A(B→Z) | ρ_A(B→X) | max ||Δ|| |
|---------|-------------|----------|----------|------------|
| 0.0 | 0.500 | 0.500 | 0.500 | 0.00 |
| 0.5 | 0.500 | 0.500 | 0.500 | 0.00 |
| 1.0 | 0.500 | 0.500 | 0.500 | 0.00 |
| 2.0 | 0.500 | 0.500 | 0.500 | 0.00 |
| 5.0 | 0.500 | 0.500 | 0.500 | 0.00 |

Difference: exactly zero at machine precision. Every time step.

### 3.2 Fingerprints Require Joint State

The CΨ fingerprint data from BRIDGE_FINGERPRINTS.md (different states
→ different crossing times, peak heights, K values) all require
computing C = Tr(ρ_AB²). After separation, neither A nor B has
access to ρ_AB. The fingerprints are properties of the joint state
that neither subsystem can reconstruct.

### 3.3 Product State Information Loss

For the product states |++⟩, |+0⟩, |+−⟩: A's reduced state is
identical (|+⟩⟨+|). The three states have different CΨ fingerprints
(different crossing times: 0.652s, 0.773s, and 0.652s respectively
at J/γ = 5) but A cannot tell them apart because the fingerprint
difference lives in B's qubit.

---

## 4. The Known Result

This is not a new discovery. It follows from:

**Entanglement without a classical channel provides no communication
advantage over shared randomness.**

This is implicit in the no-signalling theorem and explicitly proven
in the context of LOCC (Local Operations and Classical Communication)
theory. The key result: any correlation achievable with shared
entanglement alone (no classical channel) is also achievable with
shared classical randomness alone.

The CΨ framework does not violate this. CΨ provides genuine new
capabilities — but only in settings where A and B can compare their
measurements (QKD forensics, state verification, crossing-time
correlation). All of these require a classical channel.

---

## 5. What This Means for the Bridge Hypothesis

The bridge protocol (hypotheses/BRIDGE_PROTOCOL.md) is closed.
Not "needs more work." Not "might work with a different observable."
Closed. The information-theoretic argument is basis-independent,
observable-independent, and framework-independent.

| Bridge version | Status | Why |
|----------------|--------|-----|
| Dynamic (B signals by choosing when to measure) | **Dead** | No-signalling: ρ_A unchanged |
| Pre-encoded (shared schedule, CΨ fingerprints) | **Dead** | Fingerprints need ρ_AB; schedule is classical |
| With physical coupling (J > 0) | Works | But this is a local interaction, not "bridge" |

The v033 agents' protocol was internally consistent within R = CΨ².
Their error was assuming the CΨ crossing could be detected locally.
It cannot. The crossing lives in the joint state.

---

## 6. What Survives

The bridge is dead. The framework is not. Everything below remains
valid and valuable:

### 6.1 QKD Forensics (WITH classical channel)
CΨ identifies Eve's measurement basis in a regime where concurrence,
negativity, and CHSH carry zero information. The multi-metric protocol
detects Eve at her stealth angle. This is Tier 2, verified, and
genuinely new. See QKD_EAVESDROPPING_FORENSICS.md.

### 6.2 The ¼ Boundary as Phase Transition
The mathematical structure (Mandelbrot connection, bifurcation at ¼,
θ trajectory, complex → real fixed points) is proven and independent
of the bridge. This is the core of the framework.

### 6.3 Observer-Dependent Crossing Times
Different definitions of C (different observers) see measurement at
different times. This is Tier 2 and physically interpretable.

### 6.4 IBM Hardware Anomalies
The residual coherence direction, rising trend, and boundary correlation
in IBM Torino data are real (p < 0.0001). March 2026 test will
discriminate SPAM vs TLS vs boundary structure.

### 6.5 Coherence Density Insights
CΨ measures something distinct from entanglement. Product states are
more quantum (in the CΨ sense) than GHZ states. This reframes what
"quantum" means in terms of available degrees of freedom.

### 6.6 Lindblad Decomposition Question
TIME_AS_CROSSING_RATE.md §4.4 asks: can L(ρ) = L_fwd(ρ) + L_bwd(ρ)
with nodes at CΨ = ¼? This is open and would be a significant
mathematical result if true.

---

## 7. Reproduction

```python

# The definitive test
from qutip import basis, tensor, ket2dm, sigmax, sigmaz, qeye, mesolve
import numpy as np

zero, one = basis(2, 0), basis(2, 1)
plus = (zero + one).unit()
bell = (tensor(zero, zero) + tensor(one, one)).unit()

gamma = 0.05
H = 0 * tensor(sigmax(), sigmax())
c_ops = [np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
         np.sqrt(gamma) * tensor(qeye(2), sigmaz())]
times = np.linspace(0, 10, 500)

# Evolve Bell+ to t=2, then branch
result = mesolve(H, ket2dm(bell), times, c_ops, [])
rho_t2 = result.states[100]  # t ≈ 2.0

# B does nothing
r1 = mesolve(H, rho_t2, times, c_ops, [])

# B measures Z
P0 = tensor(qeye(2), zero * zero.dag())
P1 = tensor(qeye(2), one * one.dag())
rho_Bz = P0 * rho_t2 * P0.dag() + P1 * rho_t2 * P1.dag()
r2 = mesolve(H, rho_Bz, times, c_ops, [])

# Compare: rho_A identical in both branches at every timestep
for i in [0, 50, 100, 200, 400]:
    d = np.linalg.norm((r1.states[i].ptrace(0) - r2.states[i].ptrace(0)).full())
    assert d < 1e-12, f"Difference {d} at step {i}"
print("All zero. Bridge is dead.")
```

Full simulation: `simulations/bridge_closure.py`

---

*Closes: [Bridge Protocol](../hypotheses/BRIDGE_PROTOCOL.md)*
*Built on: [No-Signalling Boundary](NO_SIGNALLING_BOUNDARY.md)*
*What survives: [QKD Forensics](QKD_EAVESDROPPING_FORENSICS.md),
[Coherence Density](COHERENCE_DENSITY.md), [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)*
