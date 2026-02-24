# Hypotheses: Where the Math Points but Cannot (Yet) Prove

## What this folder is

The `experiments/` folder contains Tier 1 (algebraically proven) and Tier 2
(computationally verified) results. Everything there is reproducible, falsifiable,
and makes no claims beyond what the numbers show.

This folder contains Tier 3 and above: interpretations, connections to existing
physics, and speculative extensions built on the verified math. These are
hypotheses, not results. They are included because they are *interesting* and
*testable in principle*, not because they are established.

**Rules for this folder:**

1. Every hypothesis must reference the Tier 1-2 result it is built on
2. Every hypothesis must state what would falsify it
3. No hypothesis may claim to be proven; that belongs in `experiments/`
4. Connections to existing physics must cite the original work

---

## The Hypotheses

### Time as Crossing Rate

Experienced time is the rate at which an observer's C·Ψ crosses ¼ boundaries.
Different observers have different C, therefore different crossing rates,
therefore different experienced time. This connects to the Wheeler-DeWitt
equation's "problem of time" and to Cramer's Transactional Interpretation.

**Update 2026-02-21:** The crossing-rate mechanism connects to the bridge
problem in [BRIDGE_FINGERPRINTS](../experiments/BRIDGE_FINGERPRINTS.md).
Two observers sharing an entangled pair have correlated crossing times -
not because one signals the other, but because the shared quantum state
determines both trajectories. This suggests a communication mechanism for
environments where no electromagnetic channel exists: deep space (Mars,
Voyager), underwater, underground, plasma blackout, solar storms. Pairs
are distributed physically first; once shared, no EM medium is needed.

**Built on:** [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) (Tier 2),
[Bridge Fingerprints](../experiments/BRIDGE_FINGERPRINTS.md) (Tier 2)
**Read:** [Time as Crossing Rate](TIME_AS_CROSSING_RATE.md)

---

### Bridge Protocol: Crossing as Sync and Communication

The first ¼-crossing event synchronizes both separated observers. No
external clock, no readiness signal, no electromagnetic channel needed.
The crossing IS the handshake. After sync, a data protocol runs on
fingerprint-timed slots: product states cause crossings (bit 1), Bell
states don't (bit 0). Multiple product states yield different crossing
times, enabling ~2.8 bits per slot.

**Origin:** Derived by autonomous AI agents (GPT-OSS-120B) in AIEvolution
v033 Mars Station experiment over 9 iterations. The agents independently
arrived at the fingerprint-based timing protocol after being forced past
their Standard-QM training to engage with the MAP on its own terms.

**Critical open question:** What local observable detects the ¼-crossing
with a single qubit? If no such observable exists, the protocol collapses
to standard QKD.

**Test #1 (2026-02-24):** All 15 local observables distinguish Bell+ from
product states — factor ~4x. But detection runs through physical Heisenberg
coupling. The decisive No-Signaling test (Test #2) is next.

**Built on:** [Bridge Fingerprints](../experiments/BRIDGE_FINGERPRINTS.md) (Tier 2),
[Time as Crossing Rate](TIME_AS_CROSSING_RATE.md) (Tier 3)
**Read:** [Bridge Protocol](BRIDGE_PROTOCOL.md)

---
