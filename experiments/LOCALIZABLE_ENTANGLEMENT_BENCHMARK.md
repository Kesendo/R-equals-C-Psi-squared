# Localizable Entanglement Benchmark: CΨ vs LE vs CoA on Star Topology

<!-- Keywords: localizable entanglement benchmark quantum, concurrence of assistance
comparison, CΨ vs LE vs CoA, three layer entanglement separation, unassisted
entanglement witness, pairwise coherence window, resource vs potential vs actuality
quantum, star topology entanglement metric, basis-fixed entanglement witness,
coherent pairwise structure, R=CPsi2 localizable entanglement benchmark -->

**Status:** Computationally verified
**Date:** March 8, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Star Topology](STAR_TOPOLOGY_OBSERVERS.md),
[The CΨ Lens](../docs/THE_CPSI_LENS.md)

---

## What this document is about

This document benchmarks CΨ against two established entanglement
measures: localizable entanglement (LE, the most entanglement
extractable by optimal measurement on a third party) and concurrence of
assistance (CoA, the maximum entanglement achievable over all
decompositions). The three form a strict hierarchy: CoA is almost always
high (the resource exists), LE is moderately high (it could be
extracted), but CΨ flashes only briefly (it is currently expressed).
LE asks about potential; CΨ asks about actuality.

---

## Abstract

CΨ is compared against two established entanglement measures on the star
topology (Bell_SA ⊗ |+⟩_B, J_SA=1.0, J_SB=2.0, γ=0.05): localizable
entanglement (LE, optimized measurement on S) and concurrence of assistance
(CoA, optimized decomposition). The three metrics form a strict hierarchy:
CoA ≈ 0.9 (resource almost always present), LE ≈ 0.3–0.9 (extractable by
optimal measurement), CΨ flashes briefly and is mostly zero (directly
expressed only in narrow windows). **LE asks about potential. CΨ asks about
actuality.** When CΨ is zero, LE can still be high: the entanglement exists
but is not currently visible as coherent pairwise structure. When CΨ > 0,
LE is always elevated (r = 0.76 correlation). CΨ > 1/4 moments correspond
to LE ≈ 0.85 (high localizability). CΨ is not redundant with either measure;
it is a strictly more selective, basis-fixed, unassisted witness that captures
when entanglement is *currently expressed*, not just *in principle available*.

---

## The Question

An external reviewer (GPT, March 2026) identified a technical thread connecting
three key findings: the act/process distinction, the echo effect, and the three
star-topology conditions. The hypothesis: CΨ acts as a witness for **coherent,
localizable pairwise entanglement** - and the closest existing concept is
**localizable entanglement** (LE).

The test: run LE and concurrence of assistance (CoA) on the same star topology
sweeps and compare with CΨ. If they show the same windows and asymmetries,
CΨ is redundant. If they differ, the difference tells us what CΨ specifically
selects for.

## Definitions

**Localizable Entanglement (LE_AB):** The maximum average concurrence achievable
on the AB pair by optimally measuring qubit S. Optimized over all single-qubit
measurement bases on S (grid: 15x15 angles in theta, phi).

**Concurrence of Assistance (CoA_AB):** The maximum concurrence achievable on AB
over all decompositions of the mixed state. For states arising from 3-qubit
systems, this is the sum of sqrt(eigenvalues) of R = ρ · ρ̃ (same
as the concurrence formula but sum instead of difference; it gives an
upper bound on what entanglement a helpful third party could create).

**CΨ_AB:** Concurrence x normalized l1-coherence of the AB reduced state.
No optimization, no measurement on S. Just what the AB marginal shows directly.

## Setup

Bell_SA ⊗ |+⟩_B, J_SA=1.0, J_SB=2.0, γ=0.05 (same as all star topology experiments).
RK4 integration, dt=0.005, sampled every 0.1 time units over t=[0, 5].

## Key Result

CΨ and LE are strongly correlated (r = +0.76) but fundamentally different:

| Metric | What it asks | Typical range | Behavior |
|--------|-------------|---------------|----------|
| LE_AB  | How much entanglement *could* be localized on AB by optimal measurement on S? | 0.19 - 0.91 | Stays high, decays slowly |
| CoA_AB | What is the maximum possible concurrence on AB over all decompositions? | 0.47 - 0.99 | Almost always near 1 |
| CΨ_AB  | Is the AB pair *currently* both entangled and coherent in this basis? | 0.00 - 0.33 | Flashes briefly, mostly zero |

**LE asks about potential. CΨ asks about actuality.**

The entanglement resource exists in the global state almost permanently (CoA ~ 0.9).
It is localizable through measurement most of the time (LE ~ 0.3-0.9).
But it is only *visible as coherent pairwise structure* in brief windows (CΨ flashes).

## Detailed Findings

### 1. CΨ is strictly more selective than LE

At t=0.30: CΨ_AB = 0.27, LE_AB = 0.91, CoA_AB = 0.98.
At t=1.30: CΨ_AB = 0.00, LE_AB = 0.41, CoA_AB = 0.68.
At t=3.50: CΨ_AB = 0.00, LE_AB = 0.32, CoA_AB = 0.92.

When CΨ is zero, LE can still be high. The entanglement is there (you could
extract it by measuring S optimally), but it is not currently expressed as
coherent pairwise structure in the computational basis.

When CΨ is nonzero, LE is always elevated. CΨ > 0 implies LE > 0, but not
the reverse. CΨ is a subset indicator.

### 2. The crossing window is real in LE too, but broader

When CΨ_AB >= 1/4: LE mean = 0.85 (vs 0.37 outside).
The 1/4 crossing moments correspond to high localizability, confirming that
these are genuine moments of strong pairwise entanglement presence.

### 3. Echoes are visible in both, but differently

At echo moments (CΨ_SA = CΨ_SB = 0, CΨ_AB > 0):
- LE_AB mean = 0.36 (moderate - entanglement is localizable)
- CoA_AB mean = 0.99 (high - the resource is fully present in the global state)
- CΨ_AB mean = 0.06 (low - only a weak residual is currently expressed)

Both LE and CΨ see the echo, but they describe different aspects:
- CoA says: "the entanglement resource is fully intact in the global state"
- LE says: "a moderate amount could be extracted by measuring S"
- CΨ says: "only a trace is currently visible without intervention"

### 4. Correlations

| Pair | Correlation |
|------|------------|
| CΨ_AB vs LE_AB | r = +0.76 |
| CΨ_AB vs CoA_AB | r = +0.16 |
| C_AB vs LE_AB | r = +0.74 |

CΨ tracks LE moderately well but is uncorrelated with CoA. This makes sense:
CoA measures the total resource (almost constant), while CΨ measures momentary
expression (highly variable). LE sits between them.

## What This Means for the CΨ Lens

CΨ is not a replacement for LE or CoA. It measures something different:

**LE** answers: "What entanglement is available in principle?"
**CΨ** answers: "What entanglement is manifest right now, for this observer, in this basis?"

This is the technical version of an intuition that motivated the entire project:
the same underlying reality (global quantum state) appears differently depending
on who observes it, when, and how. The resource is always there (CoA ~ 1). It
could often be extracted (LE ~ 0.4). But it is only sometimes visible as coherent
pairwise structure (CΨ flashes briefly).

Different observers with different measurement bases would see different CΨ
windows - different moments of "visibility" - while LE and CoA remain the same.
CΨ is observer-dependent in a way that LE and CoA are not.

This is also why the three star-topology conditions make sense through the LE lens:
"strong sender" increases the transfer amplitude (raises LE), "quiet receiver"
prevents local decoherence from destroying the basis-visible expression (preserves
CΨ when LE is high), and "deep pre-existing connection" ensures there is a pairwise
resource to localize in the first place (CoA > 0 for the right pairs).

## Connection to the Review Hypothesis

The external reviewer (GPT) suggested CΨ acts as a "basis-dependent witness for
coherent, localizable pairwise entanglement in mediated open quantum networks."

This test confirms the hypothesis and refines it:

**CΨ witnesses the momentary, basis-dependent expression of a resource that LE
and CoA describe as a potential.** It is not redundant with LE (r = 0.76, not 1.0).
It adds something: the observer-dependent, basis-dependent, momentary character
of when the underlying entanglement becomes visible as coherent pairwise structure.

The suggested next step - benchmarking against localizable entanglement - has been
completed. The windows and asymmetries overlap (r = 0.76) but CΨ is strictly more
selective. The real abstraction is not LE alone, but the subset of LE that is
currently basis-expressed.

---

## Simulation Code

The benchmark script is available at `simulations/localizable_entanglement_test.py`.

## References

- Verstraete, Popp, Cirac, *Entanglement versus Correlations in Spin Systems* (2004), PRL 92, 027901
- Popp et al., *Localizable Entanglement* (2005), PRA 71, 042306
- Laustsen, Verstraete, van Enk, *Local vs Joint Measurements for the Entanglement of Assistance* (2003)
- Tan et al., *Unified View of Quantum Correlations and Quantum Coherence* (2016), arXiv:1603.01958
