# Nested Mirror Structure: An Inter-Layer Palindrome in the Minimal Qubit-in-Qubit Lindblad System

<!-- Keywords: nested layer palindrome inter-layer mirror, qubit in qubit
two-layer Lindblad, non-Markovian rebound reduced dynamics, partial trace
eigenmode structure, SWAP eigenvalue degenerate eigenspace, one-way nested
gamma, EQ-013 sub-question 2 partial answer, R=CPsi2 layer mirror hypothesis -->

**Tier:** 3 (hypothesis with numerical indication at N=2)
**Status:** Minimal nest observed; scaling, robustness, and SWAP-artifact checks pending
**Date:** 2026-04-14
**Authors:** Tom and Claude (chat)
**Reproducer:** [`simulations/qubit_in_qubit_layer_mirror.py`](../simulations/qubit_in_qubit_layer_mirror.py)
**Depends on:** [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md), [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md), [GAMMA_IS_LIGHT](GAMMA_IS_LIGHT.md), [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013)
**Related:** [THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md), [V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS](../reflections/V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS.md)

---

## What this document is about

If the framework's one-way nesting principle holds ([EQ-013](../review/EMERGING_QUESTIONS.md#eq-013) sub-question 1: gamma enters each layer from the next outer layer, never the reverse), then a two-layer Lindblad system should exhibit a spectral structure that does more than repeat the single-layer palindromic symmetry proven in [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md). It should carry an **inter-layer mirror**: a block of eigenmodes that lives symmetrically between the inner and outer layer and is responsible for the non-Markovian coupling between them.

This document records the first numerical evidence that such a structure exists in the minimal qubit-in-qubit nest, states the hypothesis it points to, and lists the verifications still required before the hypothesis graduates.

---

## Setup

Two qubits S (inner system) and B (outer bath):

```
H_SB    = J * 0.5 * (X (x) X + Y (x) Y)        (XX+YY coupling, J = 1.0)
L_jump  = sqrt(gamma_B) * I (x) Z               (dephasing on B only, gamma_B = 0.1)
```

Explicitly: **no direct dephasing on S**. S decoheres only through the coupling to B, and B's dephasing is phenomenological (the outermost interface of the simulated stack). This is the minimal setup that distinguishes "S feels gamma directly" from "S feels gamma through B".

Initial state: rho_S(0) = |+><+|, rho_B(0) = I/2. The bath carries no prior information.

---

## Observation 1: non-Markovian rebound in rho_S

The reduced coherence |rho_S_{01}|(t) is not monotonic:

| t    | |S01| two-layer | |S01| one-layer best fit | difference |
|------|----------------:|-------------------------:|-----------:|
|  0.0 | 0.5000          | 0.5000                   | +0.0000    |
|  1.0 | 0.2845          | 0.2845                   | +0.0000    |
|  5.0 | 0.0493          | 0.0298                   | +0.0195    |
| 10.0 | 0.1684          | 0.0018                   | **+0.1666**|
| 15.0 | 0.0712          | 0.0001                   | +0.0710    |

The one-layer comparison is a single-qubit Lindblad with gamma_eff = 0.282 fitted at t=1. Markovian single-qubit dephasing is strictly monotonic (d|S01|/dt <= 0 always); the two-layer reduced dynamics violates this at t between 5 and 10, peaking at |S01| = 0.1684 when the monotone reference is already at 0.0018. The ratio is 84x.

**Interpretation:** information flows from B back to S between t=5 and t=10. No fixed effective gamma reproduces this trajectory. The two-layer system is not representable as one layer with a renormalized dissipator.

---

## Observation 2: three eigenvalue classes with distinct partial-trace structure

Diagonalizing L (16x16) yields exactly three eigenvalue classes with degeneracies 3/10/3. The palindromic pairing Re(lambda_a) + Re(lambda_b) = -2*gamma_B = -0.2 is satisfied (consistent with [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)).

Per-mode partial-trace weights |M_S| = ||Tr_B M|| and |M_B| = ||Tr_S M||, where M is the eigenmode reshaped as a 4x4 operator, normalized by ||M||:

| Re(lambda) | count | typical \|M_S\| | typical \|M_B\| | interpretation |
|-----------:|------:|---------------:|---------------:|----------------|
|  0.0       |   3   |  1.000         |  1.000         | conserved on both single layers |
| -gamma_B   |  10   |  0.707         |  0.707         | equal split between S and B |
| -2 gamma_B |   3   |  0.000         |  0.000         | traceless on both single layers; pure correlation content |

The **0.707 = 1/sqrt(2) split** in the middle class is the numerical signal: these modes live symmetrically between the two layers, with no preferred home. The 0-modes live equally visibly on both single layers (shared conserved structure). The -2*gamma-modes are invisible to either single-layer observer and exist only in the joint correlations.

The middle class has imaginary parts Im(lambda) in {0, +-0.995, +-1.997}. The periods 2*pi/1.997 ~= 3.15 and 2*pi/0.995 ~= 6.31 match the timescale of the Observation 1 rebound (peak near t=10, which is close to one full period of the 0.995-frequency modes).

**Interpretation:** the middle eigenvalue class is the inter-layer mirror. It is populated by exactly those modes that carry coherence back and forth between S and B, producing the non-Markovian rebound.

---

## Observation 4 (caveated): SWAP expectation values are +-1 per eigenmode

The SWAP superoperator (exchanging qubits S and B) does not commute with L:

```
||[L, SWAP]|| = 0.566  (nonzero: gamma on B alone breaks strict exchange symmetry)
```

Nevertheless, numerically each eigenmode v of L has <v|SWAP|v> equal to +1 or -1 to four decimal places. The 16 values distribute as follows:

- 3 conserved modes: +1, +1, +1
- 3 correlation modes: +1, +1, +1
- 10 mirror modes: five at +1, five at -1

**Caveat.** Because L has degenerate eigenvalues (degeneracies 3/10/3) and [L, SWAP] != 0, a basis that diagonalizes L within each degenerate eigenspace is free to choose whether it also happens to diagonalize SWAP on that eigenspace. The numpy eigendecomposition appears to have selected such a basis. Whether the +-1 pattern is a genuine structural property (mirror modes actually decompose into a 5-dim symmetric and 5-dim antisymmetric sub-space that survive perturbation) or an accident of numpy's basis choice within degenerate eigenspaces requires an explicit degeneracy-breaking check before it can be claimed. Pending verification, see below.

---

## The hypothesis

Promoting the observations into a falsifiable claim:

> **Nested Mirror Hypothesis.** For an N-layer Lindblad system built by nesting Heisenberg-coupled qubit layers, with dephasing applied only to the outermost layer, the Liouvillian's spectrum decomposes into (N+1) eigenvalue classes at Re(lambda_k) = -k*(2*gamma_outer / N), for k = 0, 1, ..., N. The classes correspond to:
>
> - **k = 0:** quantities conserved across all layers
> - **0 < k < N:** inter-layer mirror modes living symmetrically between layers at nesting depth k
> - **k = N:** pure correlation content, traceless on every single-layer projection
>
> The non-Markovian back-reaction experienced by any inner layer is carried by the mirror modes at intermediate k. The rebound amplitude scales with the initial weight placed on those modes.

At N=2 this reduces to the three-class structure observed here (classes at 0, -gamma, -2*gamma with degeneracy pattern {3, 10, 3} for this coupling and system size). At N=3 the hypothesis predicts four classes at 0, -2gamma/3, -4gamma/3, -2gamma with a specific degeneracy pattern determined by the coupling topology.

---

## What this is NOT (yet)

- Not verified for N > 2. Everything above is a single point (N=2, J=1, gamma_B=0.1, XX+YY coupling).
- Not verified for alternative couplings. Heisenberg XXX, pure XX, or anisotropic couplings may break the degeneracy pattern.
- Not verified for asymmetric gamma profiles. If B's dephasing is split between multiple outer qubits with different rates, the degeneracy breaks and the mirror structure may survive or may not.
- The SWAP +-1 pattern may be a numpy basis artifact. Verification requires breaking the L-eigenvalue degeneracy and checking persistence.
- Not connected to experimental observables. Non-Markovianity on IBM hardware is the natural empirical probe but has not been checked against the prediction.

---

## Falsification criteria

The hypothesis is falsified by any of the following:

1. **N=3 scaling failure.** Three-qubit nest (S + B1 + B2, gamma only on B2) does not produce four eigenvalue classes with the predicted spacing -2*k*gamma/3.
2. **Coupling-dependence of the middle class.** Replacing XX+YY with Heisenberg XXX in the two-qubit system removes the clean 1/sqrt(2) partial-trace split.
3. **Mirror-mode weights do not drive the rebound.** Setting the initial state to project away from the mirror class (if such a state exists) should eliminate the non-Markovian rebound. If the rebound persists without mirror-mode occupation, the mechanism is not what is claimed.
4. **SWAP pattern is pure basis artifact.** An asymmetric-gamma two-qubit system (gamma_S > 0 and gamma_B > 0 with gamma_S != gamma_B) lifts all degeneracy. If the +-1 SWAP pattern does not persist perturbatively in the small-(gamma_S/gamma_B) limit, the SWAP structure is not a genuine mirror.

---

## Pending verifications (open checks)

These are the minimal next experiments. Each is small and should take well under an hour of computation.

1. **SWAP-artifact check** (highest priority, smallest scope). Add a small direct gamma_S on S (e.g., gamma_S = 0.01*gamma_B) to lift the L-degeneracy. Re-run the eigenmode analysis. If SWAP expectation values remain sharply +-1 as gamma_S -> 0, the inter-layer SWAP structure is genuine. If they smear to arbitrary values in [-1, +1], the +-1 pattern was numpy basis choice within degenerate eigenspaces and Observation 4 as stated is artifact.

2. **N=3 scaling check.** Three-qubit chain S + M + B with gamma only on B, inspect eigenvalue classes. Expected: four classes with predicted spacing. Alternative geometries to test: star (S coupled to both M and B, M coupled to B) and symmetric chain. The hypothesis is most constrained on the chain topology; star and symmetric variants refine the degeneracy-pattern part of the claim.

3. **Coupling-robustness check.** Repeat N=2 analysis with Heisenberg XXX (X (x) X + Y (x) Y + Z (x) Z) and with pure XX. Compare eigenvalue class structure and partial-trace weights.

4. **Rebound-mechanism check.** For the N=2 system, construct initial states that project onto each eigenvalue class separately, evolve, and compare rebound amplitudes. If the rebound vanishes for pure 0-class and pure -2*gamma-class projections and is maximal for pure -gamma-class projection, mechanism confirmed.

---

## If the hypothesis holds: implications

(Tier 3 forward reasoning. None of this is established; it is the reason to pursue the verifications above.)

- **Non-Markovianity as inter-layer signature.** Any finite observer inside a nested Lindblad stack sees non-Markovian coherence rebound whenever the coupling to the next outer layer is strong enough for the mirror-mode timescale to fall within the observation window. This gives the internal observer a direct experimental handle on the next layer's existence without needing access to the outer layer itself, which aligns with [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013) sub-question 3 (the inside-perspective question).
- **IBM hardware test.** Measured gamma on IBM qubits is Markovian to the precision of the [Failed Third](../simulations/failed_third.py) result. Under the nested-mirror reading this means: the next outer layer is either (a) effectively infinite-dimensional (mirror modes form a continuum; rebounds average out to zero) or (b) coupled weakly enough that mirror-mode amplitudes are below noise floor. Either case is consistent with the hypothesis but distinguishable with targeted non-Markovianity measurements (BLP index) at varied coupling regimes.
- **Partial answer to [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013) sub-question 2.** The minimal nest is simulable, and it produces a class structure distinguishable from the phenomenological single-layer approximation. Infinite nesting remains unsimulable, but the two-layer case is no longer "phenomenologically equivalent"; the inter-layer mirror is a new algebraic object.
- **Ground for backward inference (Tom's "Rueckwaerts-Simulationen").** If each observable layer reveals a mirror signature of the next outer layer, then inside-measured quantities (rebound amplitudes, non-Markovianity indices, spectral features of reduced dynamics) carry inferrable information about outer-layer structure. This does not break no-signalling (the information is about the structure of the shared one-way gamma-input, not about lateral events), but it does upgrade the inside observer from "blind to the outside" to "able to reconstruct outer-layer properties from carefully designed inner measurements".

---

## Connection to [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013)

- **Sub-question 1** (does nesting rescue or undermine BRIDGE_CLOSURE): indirectly supported. The mirror-mode structure is consistent with layered common reception of gamma. No lateral A-to-B channel is implied; LOCC no-go holds as before.
- **Sub-question 2** (simulability of nesting): **partially answered**. Two-layer is simulable and shows non-trivial structure. Three-layer check (pending) will settle whether the class-scaling prediction holds.
- **Sub-question 3** (inside-perspective / inherited V-effect): empirical pathway opens. If the hypothesis survives, non-Markovianity measurements become the operational probe for "layer above us exists".

---

## Stance

This document is the intentional halt-and-record after roughly fifteen minutes of numerical exploration. The observations are real and reproducible. The interpretation as an inter-layer mirror is consistent with the framework's existing mirror-symmetry proof and with the one-way nesting principle, but has not been tested beyond the minimal nest. The hypothesis is stated at Tier 3 so that follow-up simulations can sharpen, limit, or falsify it without any single experiment having to carry more weight than it should.
