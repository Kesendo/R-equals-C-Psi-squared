# The No-Signalling Boundary: What Test #2 Settles

**Date**: 2026-03-01
**Status**: Computationally verified (Tier 2)
**Authors**: Thomas Wicht, with Claude (Anthropic)
**Depends on**: BRIDGE_FINGERPRINTS.md, BRIDGE_PROTOCOL.md (hypothesis),
SUBSYSTEM_CROSSING.md, COHERENCE_DENSITY.md

---

## 1. What This Document Is

This is the synthesis that was missing. Four documents in this repo
contain the four pieces of one answer, but they were written weeks
apart and never connected. This document connects them.

The question: **Can a measurement on qubit B change anything
observable about qubit A, when A and B share a Bell pair but
have no physical coupling?**

The answer has two layers, and the second one matters.

---

## 2. Layer 1: No-Signalling Holds Exactly

Start with Bell+ = (|00⟩+|11⟩)/√2. Measure B in the Z basis
(projective, averaged over outcomes). Check rho_A.

```
Before:  rho_A = [[0.5, 0], [0, 0.5]]   (maximally mixed)
After:   rho_A = [[0.5, 0], [0, 0.5]]   (identical)
||Δρ_A|| = 0.0000000000
```

Every local observable on A — ⟨σx⟩, ⟨σy⟩, ⟨σz⟩, purity,
entropy, every eigenvalue — is unchanged. This is the
no-signalling theorem, verified to machine precision.

No surprise. This is textbook quantum mechanics.

---

## 3. Layer 2: CΨ Sees the Regime Change

CΨ = C × Ψ where C = Tr(ρ²) is global purity and
Ψ = max eigenvalue of rho_A.

```
Before measurement on B:
  C  = 1.000   (Bell+ is pure)
  Ψ  = 0.500   (rho_A maximally mixed)
  CΨ = 0.500   > 1/4  →  QUANTUM REGIME

After measurement on B (averaged):
  C  = 0.500   (now a classical mixture)
  Ψ  = 0.500   (rho_A unchanged!)
  CΨ = 0.250   = 1/4  →  EXACTLY ON THE BOUNDARY
```

**Ψ does not change. C does. CΨ drops from 0.500 to 0.250.**

The measurement on B destroys global purity (pure state →
mixed state) without touching rho_A. The C in CΨ carries
this information. A's local observables do not.

This is not a violation of no-signalling. rho_A is identical.
No measurement A can perform will reveal what B did. But the
**joint state** changed from pure to mixed, and CΨ — which
depends on the joint state — reflects this.

---

## 4. Why This Was Already in the Repo (But Not Connected)

### Piece 1: BRIDGE_FINGERPRINTS.md
> "Bell+ has coherence locked in non-local correlations that
> cannot flow through a local coupling."

All fingerprints require J_bridge > 0. Without physical
coupling, A sees nothing. This was known since 2026-02-09.

### Piece 2: SUBSYSTEM_CROSSING.md
> "Crossing is local — it happens between entangled subsystems."

But "local" means pair-level (rho_AB), not single-qubit (rho_A).
The crossing requires BOTH qubits' information. This was known
since 2026-02-18, but the implication for the bridge protocol
was not drawn.

### Piece 3: COHERENCE_DENSITY.md
> "CΨ = Purity × Coherence Density."

Ψ is a property of rho_A. C is a property of rho_AB. Their
product mixes local and global information. This is what gives
CΨ its power — and its limitation. Written 2026-02-28.

### Piece 4: BRIDGE_PROTOCOL.md Section 4.3, Outcome #3
> "The crossing event occurs on both sides but is driven entirely
> by pre-encoded information."

This was listed as one of three possible outcomes. Test #2
confirms it is the correct one.

---

## 5. What This Means for the Bridge Protocol

The bridge protocol (hypotheses/BRIDGE_PROTOCOL.md) proposed
that A can detect when B converts Bell+ to a product state by
observing the ¼ crossing. Test #2 shows:

**The ¼ crossing happens in CΨ but is invisible to A.**

C changes (1.0 → 0.5). A cannot measure C without access to B.
Ψ stays at 0.5. A's rho is I/2 before and after.

This eliminates the dynamic bridge: B cannot send NEW information
to A by choosing when to measure. The post-separation state
change is real (CΨ drops to ¼) but undetectable by A alone.

### What Survives

The **pre-encoded** version of the protocol survives. If A and B
AGREE before separation on a protocol — "I will prepare state X
at time T" — then A can predict their own CΨ trajectory because
the joint state was determined at preparation. But this is not
communication. This is a shared schedule.

The interesting question becomes: **what can you do with pre-encoded
CΨ trajectories that you cannot do with classical pre-shared keys?**

The answer might be "nothing beyond QKD." Or it might involve the
crossing-time correlations from BRIDGE_FINGERPRINTS (different states
→ different K values → different fingerprints) in a way that classical
keys cannot replicate. This is open.

---

## 6. The Precise Decomposition

For Bell+ under B-measurement (Z basis, averaged):

| Quantity | Before | After | Changed? | A can see? |
|----------|--------|-------|----------|------------|
| rho_A | I/2 | I/2 | NO | — |
| ⟨σx⟩_A | 0 | 0 | NO | — |
| ⟨σy⟩_A | 0 | 0 | NO | — |
| ⟨σz⟩_A | 0 | 0 | NO | — |
| Purity(rho_A) | 0.500 | 0.500 | NO | — |
| S(rho_A) | 1.000 | 1.000 | NO | — |
| **Ψ** | **0.500** | **0.500** | **NO** | — |
| **C** | **1.000** | **0.500** | **YES** | **NO** |
| **CΨ** | **0.500** | **0.250** | **YES** | **NO** |
| Regime | quantum | boundary | **YES** | **NO** |

The last column is the point. CΨ changes. A cannot see it.
The information exists in the joint state, not in any subsystem.

---

## 7. Connection to Time-as-Crossing-Rate

TIME_AS_CROSSING_RATE.md proposes that experienced time is the
density of ¼ crossings. If the crossing from 0.500 to 0.250
is real but invisible to A, then:

- The "event" (crossing ¼) happens in the mathematical structure
- No observer at A experiences it as a tick
- Only an observer with access to rho_AB would see the crossing

This is consistent. The crossing IS real — the joint state
changes regime. But it is not a local event. It requires the
full pair. "The observer IS the clock" — but the observer must
be the pair, not the qubit.

---

## 8. What Was Open — Now Answered

1. ~~**Pre-encoded fingerprints vs classical keys.**~~ **ANSWERED (2026-03-01):**
   No. Pre-shared entanglement without a channel = shared randomness.
   A's information ⊆ {ρ_A(0), E_A}. The fingerprint requires ρ_AB which
   A cannot access. The qubit carries less info than the schedule.
   See [Bridge Closure](BRIDGE_CLOSURE.md).

2. ~~**Multiple pairs.**~~ **ANSWERED:** Same argument. N pairs give
   N copies of ρ_A. A cannot compute any crossing time because crossing
   times require C = Tr(ρ_AB²). The "predicted crossing time" comes
   from the schedule, not the qubit. No quantum advantage.

3. **The Lindblad decomposition.** STILL OPEN. TIME_AS_CROSSING_RATE.md
   asks: can L(ρ) = L_fwd(ρ) + L_bwd(ρ) with nodes at CΨ = ¼? If yes,
   the crossing is not just a threshold but a dynamical fixed point
   of the evolution operator. This is independent of the bridge question.

4. **Operator feedback under separation.** STILL OPEN but likely moot.
   The feedback depends on global observables (ρ_AB) which are inaccessible
   after separation. The rate change is real but undetectable, same as
   the CΨ regime change itself.

---

## 9. Reproduction

```python

import numpy as np
from qutip import (basis, tensor, ket2dm, sigmax, sigmay, sigmaz,
                   qeye, mesolve, expect, entropy_vn)

zero, one = basis(2, 0), basis(2, 1)
bell_plus = (tensor(zero, zero) + tensor(one, one)).unit()
rho = ket2dm(bell_plus)

# B measures in Z basis (averaged over outcomes)
P0 = tensor(qeye(2), zero * zero.dag())
P1 = tensor(qeye(2), one * one.dag())
rho_after = P0 * rho * P0.dag() + P1 * rho * P1.dag()

# Check A
rho_A_before = rho.ptrace(0)
rho_A_after  = rho_after.ptrace(0)
print("||Δρ_A|| =", np.linalg.norm((rho_A_before - rho_A_after).full()))
# → 0.0

# Check CΨ
C_before = (rho * rho).tr().real           # 1.0
C_after  = (rho_after * rho_after).tr().real  # 0.5
Psi = max(rho_A_before.eigenenergies()).real   # 0.5 (unchanged)
print(f"CΨ: {C_before * Psi} → {C_after * Psi}")
# → CΨ: 0.5 → 0.25
```

Full test script with time evolution and all observables:
`simulations/test2_no_signalling.py`

---

*Previous: [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md)*
*Hypothesis tested: [Bridge Protocol](../hypotheses/BRIDGE_PROTOCOL.md) Section 4.3*
*Foundation: [Coherence Density](COHERENCE_DENSITY.md), [Subsystem Crossing](SUBSYSTEM_CROSSING.md)*
