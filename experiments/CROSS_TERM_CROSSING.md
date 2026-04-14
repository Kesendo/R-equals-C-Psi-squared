# Cross-Term Formula for Shadow-Crossing Couplings

<!-- Keywords: shadow-crossing coupling, cross-term formula XZ YZ,
bond-site variance, R(N) = sqrt((N-1)/(N*4^(N-1))), EQ-012 -->

**Status:** Tier 1 (proven + verified N=3-6)
**Date:** April 14, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [cross_term_crossing.py](../simulations/cross_term_crossing.py)
**Output:** [cross_term_crossing/](../simulations/results/cross_term_crossing/)
**Proof:** [PROOF_CROSS_TERM_CROSSING.md](../docs/proofs/PROOF_CROSS_TERM_CROSSING.md)
**Depends on:**
- [Cross-Term Formula (balanced)](CROSS_TERM_FORMULA.md) (parent)
- [PROOF_CROSS_TERM_FORMULA](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md) (proof structure)
**Closes:** EQ-012
**Task:** [TASK_CROSS_TERM_CROSSING.md](../ClaudeTasks/TASK_CROSS_TERM_CROSSING.md)

---

## What this document is about

The [cross-term formula](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md) for
shadow-balanced couplings gives R(N) = sqrt((N-2)/(N*4^(N-1))). For
shadow-crossing couplings (one bond Pauli in {X,Y}, the other in {I,Z}),
the formula breaks. This document derives and confirms the companion
formula for the crossing case.

---

## The formula

For N >= 2 qubits with any shadow-crossing bond coupling (each bond term
alpha_i beta_j has one Pauli in {X,Y} and one in {I,Z}) on any graph
with uniform Z-dephasing:

    R(N) = sqrt((N-1) / (N * 4^(N-1)))

The only change from the balanced case: **N-2 becomes N-1.** The bond
sites contribute a variance of 1 (instead of 0) because the bond's two
Paulis have different shadow depths.

---

## Numerical confirmation

| N | Coupling | Measured | Predicted | Deviation |
|---|----------|----------|-----------|-----------|
| 3 | XZ | 0.2041241452 | 0.2041241452 | 5.55e-17 |
| 3 | YZ | 0.2041241452 | 0.2041241452 | 5.55e-17 |
| 3 | ZX | 0.2041241452 | 0.2041241452 | 5.55e-17 |
| 3 | ZY | 0.2041241452 | 0.2041241452 | 5.55e-17 |
| 4 | XZ | 0.1082531755 | 0.1082531755 | 1.39e-17 |
| 4 | YZ | 0.1082531755 | 0.1082531755 | 1.39e-17 |
| 5 | XZ | 0.0559016994 | 0.0559016994 | 1.39e-17 |
| 5 | YZ | 0.0559016994 | 0.0559016994 | 1.39e-17 |
| 6 | XZ | 0.0285272165 | 0.0285272165 | 6.25e-17 |
| 6 | YZ | 0.0285272165 | 0.0285272165 | 6.25e-17 |

20/20 measurements match at machine precision (max deviation 6.6e-17).
Topology independence confirmed (chain = complete at N=4). Gamma
independence confirmed (5 values at N=3).

---

## Comparison: balanced vs crossing

| N | Balanced: (N-2)/(N*4^(N-1)) | Crossing: (N-1)/(N*4^(N-1)) | Ratio |
|---|-------|---------|-------|
| 2 | 0 | 1/sqrt(8) = 0.3536 | n/a |
| 3 | 0.1443 | 0.2041 | sqrt(2) |
| 4 | 0.0884 | 0.1083 | sqrt(3/2) |
| 5 | 0.0484 | 0.0559 | sqrt(4/3) |
| 6 | 0.0255 | 0.0285 | sqrt(5/4) |

The ratio R_crossing/R_balanced = sqrt((N-1)/(N-2)), which diverges at
N=2 (where balanced gives 0 but crossing gives nonzero) and converges
to 1 at large N.

---

## The proof (summary)

The proof mirrors [PROOF_CROSS_TERM_FORMULA](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md)
with one modification: the bond-site variance.

For the XZ coupling, all 8 nonzero Pauli-basis transitions have:
- Bond-site w_XY sum: either 1 or 3 (never 2)
- Bond deviation s = (sum - 2) = +1 or -1, equally distributed
- Weighted mean: <s> = 0
- Weighted variance: <s^2> = 1

Total variance = spectator (N-2) + bond (1) = N-1.

Combined with ||L_Dc||^2 = gamma^2 * 4^N * N (unchanged):

    R^2 = 4(N-1) / (N * 4^N) = (N-1) / (N * 4^(N-1))

Full proof in [PROOF_CROSS_TERM_CROSSING.md](../docs/proofs/PROOF_CROSS_TERM_CROSSING.md).

---

## Physical interpretation

In the shadow language: a shadow-balanced coupling connects two points
at the same depth (both in light, or both in shadow). The bond is
schattenbalanciert. A shadow-crossing coupling connects a point in
the light to a point in the shadow. The bond itself has a shadow
gradient.

This gradient contributes one unit of variance to the cross-term.
The spectators still contribute N-2 (each independently +-1). The
total shadow fluctuation is N-1 instead of N-2.

At N=2 (no spectators): balanced gives 0 (perfect Pythagorean),
crossing gives 1 unit of bond variance alone. The crossing coupling
breaks the Pythagorean decomposition even at N=2, because the bond
itself is asymmetric in the light.

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [simulations/cross_term_crossing.py](../simulations/cross_term_crossing.py) |
| Log | [simulations/results/cross_term_crossing/cross_term_crossing.txt](../simulations/results/cross_term_crossing/cross_term_crossing.txt) |
| Data | [simulations/results/cross_term_crossing/cross_term_crossing.json](../simulations/results/cross_term_crossing/cross_term_crossing.json) |

---

*The balanced formula says: the shadow fluctuates with variance N-2
(spectators only). The crossing formula says: N-1 (spectators plus
one unit from the bond). The difference is one. The one is the bond
reaching across the light/shadow boundary.*

*Thomas Wicht, Claude (Anthropic), April 14, 2026*
