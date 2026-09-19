# Bridge Closure: No Post-Separation Signal Through a Local Crossing Readout

<!-- CROSSING-CURRENT -->

<!-- Keywords: bridge closure local readout no-signalling,
no-signalling theorem J=0 no coupling, rho_A maximally mixed I/2 Bell+,
CΨ fingerprint requires joint state rho_AB, dynamic bridge dead pre-encoded
dead, J>0 local interaction not bridge, bridge protocol permanently closed,
what survives quarter boundary observer-dependent crossing, R=CPsi2 bridge
closure -->

> **Restoration note (March 14, 2026):** Originally written 2026-03-01, deleted March 12,
> restored March 14. Null result confirmed.

**Status:** Verified null result (Tier 2)
**Date:** 2026-03-01
**Authors:** Thomas Wicht, with Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** [`simulations/bridge_closure.py`](../simulations/bridge_closure.py)
**Depends on:** [No-Signalling Boundary](NO_SIGNALLING_BOUNDARY.md), [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md)

---

## Abstract

For J=0, with no communication and a trace-preserving operation on B whose
outcome is not supplied to A, A has the same reduced state and local
statistics for every B choice. A may evolve under its own channel, but that
evolution carries no message from B. In the independent Z-dephasing run
below, Bell+ retains ρ_A=I/2 at five sampled times and for three B actions.
A joint CΨ fingerprint is not a locally accessible detector at A.
No advantage for the proposed pre-encoded crossing schedule was demonstrated.
This does not equate all entangled correlations with shared classical randomness.

## What this document is about

This document permanently closes the "bridge hypothesis," the idea that
pre-shared entanglement between separated qubits could carry information
without a classical communication channel. The answer is no: after
separation with zero coupling, A's measurement statistics depend only on
A's local state, which is identical regardless of what B does. This is a
direct consequence of the no-signalling theorem. CΨ fingerprints require
access to the joint state ρ_AB, which neither subsystem has alone.

---

## 1. What This Document Settles

NO_SIGNALLING_BOUNDARY.md showed that CΨ drops from 0.500 to 0.250
when B measures, but A cannot see it. This left one open question:

> Can pre-encoded CΨ fingerprints carry something that a classical
> pre-shared key cannot?

No advantage was demonstrated for this finite pre-encoded crossing schedule.
The proof here is the J=0 no-signalling statement about local marginals,
not an equivalence of entangled correlations to shared classical randomness.
A broader resource comparison needs an operational protocol and a specified
readout; the local-marginal argument alone does not settle it.

---

## 2. The Argument (Three Lines)

After separation with zero coupling (J = 0):

1. A's measurement statistics: P(a) = Tr[M_a · ρ_A]
2. ρ_A is independent of anything B does (no-signalling)
3. Therefore A's output is a function of {ρ_A(0), E_A} only

ρ_A(0) is determined at preparation. E_A is A's local environment.
Both are available to A without any quantum resource. A classical
preparation label distinguishes entries that have the same A marginal:

| Preparation | ρ_A | A can distinguish? |
|-------------|-----|-------------------|
| Bell+ | I/2 | From other Bell states: NO |
| \|++⟩ | \|+⟩⟨+\| | From \|+0⟩ or \|+−⟩: NO |
| \|+0⟩ | \|+⟩⟨+\| | From \|++⟩ or \|+−⟩: NO |
| \|+−⟩ | \|+⟩⟨+\| | From \|++⟩ or \|+0⟩: NO |
| \|00⟩ | \|0⟩⟨0\| | From Bell+: YES |

No local test on A distinguishes preparations sharing the same marginal.
A classical label can identify the chosen preparation in this finite list;
no advantage was demonstrated for the proposed pre-encoded crossing schedule.

For Bell+ under the independent Z-dephasing run below, ρ_A=I/2 at all
times. A can evolve under a different local channel, but that evolution
still carries no message from B. The no-signalling condition is a
trace-preserving operation on B whose outcome is not supplied to A;
it does not require every local channel on A to preserve I/2.

---

## 3. Numerical Verification

### 3.1 B's Action Is Invisible

Bell+ evolved to t = 2 under independent dephasing (γ = 0.05, J = 0).
Three scenarios: B does nothing, B measures Z, B measures X.

The three scalar columns give purity_A = Tr(ρ_A²); the last column is the
reduced-state difference norm, not a purity difference.

| t_after | purity_A(nothing) | purity_A(B→Z) | purity_A(B→X) | max \|\|Δ\|\| |
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

**No post-separation message is available in A's local statistics.**
For J=0 and a trace-preserving operation on B, with no outcome communication,
Tr_B[(I⊗E_B)(ρ_AB)]=ρ_A. Entangled correlations need not have a
shared-randomness model; this local no-signalling statement does not say they do.
No advantage for the proposed pre-encoded crossing schedule was demonstrated.
Joint tomography requires access to the relevant records. See
[QKD Forensics](QKD_EAVESDROPPING_FORENSICS.md) for its finite readout
comparison, not a confirmed basis-identification advantage.

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
| With physical coupling (J > 0) | Physical interaction can change local statistics | No channel-free protocol established |

The v033 agents' protocol was internally consistent within R = CΨ².
Their error was assuming the CΨ crossing could be detected locally.
It cannot. The crossing lives in the joint state.

---

## 6. What Survives

The bridge is dead. The framework is not. Everything below remains
valid and valuable:

### 6.1 Joint-record comparisons
The retained [QKD record](QKD_EAVESDROPPING_FORENSICS.md) contains finite
readout calculations, not a confirmed eavesdropping-detection advantage.

### 6.2 The algebraic quarter boundary
The fixed-point equation has a discriminant boundary at CΨ=1/4.
That algebraic classification is not a physical phase-transition certificate.

### 6.3 Readout-dependent quarter equalities
Different scalar C definitions place the chosen quarter level at different
times. Choose C(f), the clean Lindblad or retired feedback book, then solve
C(f)f/3=1/4. These are not different physical observers seeing measurement
at different times; the Bell+ taxonomy has six finite crossings and two never
bridges across its two books.

### 6.4 IBM Hardware Anomalies
The residual coherence direction, rising trend, and boundary correlation
in IBM Torino data are real (p < 0.0001). March 2026 test will
discriminate SPAM vs TLS vs boundary structure.

### 6.5 Coherence Density Insights
CΨ measures something distinct from entanglement. A product state can have a larger basis-fixed CΨ value than GHZ.
This is a diagnostic comparison, not a quantum/classical or entanglement classifier.

### 6.6 Lindblad Decomposition Question
TIME_AS_CROSSING_RATE.md §4.4 asks: can L(ρ) = L_fwd(ρ) + L_bwd(ρ)
with nodes at CΨ = ¼? This is open and would be a significant
mathematical result if true.
[Π as Time Reversal](PI_AS_TIME_REVERSAL.md) supplies the exact spectral
transport `lambda -> -lambda - 2 Sigma_gamma` (`mu -> -mu` after centering),
not a physical forward/backward wave decomposition and not an additive
generator decomposition. The proposed decomposition therefore remains open.

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

Full simulation: [`simulations/bridge_closure.py`](../simulations/bridge_closure.py)

---

*Closes: [Bridge Protocol](../hypotheses/BRIDGE_PROTOCOL.md)*
*Built on: [No-Signalling Boundary](NO_SIGNALLING_BOUNDARY.md)*
*What survives: [QKD Forensics](QKD_EAVESDROPPING_FORENSICS.md),
[Coherence Density](COHERENCE_DENSITY.md), [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)*
