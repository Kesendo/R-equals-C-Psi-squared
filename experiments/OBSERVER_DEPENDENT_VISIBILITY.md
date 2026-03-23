# Observer-Dependent Visibility: How Local Noise Profiles Shape CΨ Windows

<!-- Keywords: observer dependent visibility quantum, local noise profile CΨ window,
dephasing rate visibility threshold, proper time non-universality quantum, multi-scale
open system dynamics, unassisted entanglement witness, basis-fixed pairwise coherence,
CoA LE CΨ three layers, quiet receiver visibility window, noise environment quantum
crossing, R=CPsi2 observer dependent visibility -->

**Status:** Computationally verified
**Date:** March 8, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Star Topology](STAR_TOPOLOGY_OBSERVERS.md),
[Localizable Entanglement Benchmark](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md)

---

## Abstract

The same star topology (Bell_SA ⊗ |+⟩_B, J_SA=1.0, J_SB=2.0) is evolved
under four different local noise profiles (γ_A = 0.03, 0.05, 0.10, 0.20,
other γ fixed at 0.05). At the same coordinate time, observers with different
γ_A see different CΨ_AB values: peak CΨ ranges from 0.340 (slow observer)
to 0.261 (fast observer). A slightly noisier observer would see the peak
drop below 1/4 entirely. Plotting CΨ_AB against proper time τ_A = γ_A·t
does **not** collapse the curves: CΨ is not a function of local proper time
alone. It depends on the absolute noise rate because the dynamics have
multiple competing timescales (J, γ_A, γ_B, γ_S). **Important correction:**
different γ_A values produce different Lindblad trajectories, not different
views of the same evolved state. The correct statement: different local noise
profiles open or close windows in which a persistent entanglement resource
(CoA ≈ 1) is directly visible in a given pair and basis. CΨ is a basis-fixed,
unassisted witness of directly expressed pairwise entanglement.

> **Update March 14, 2026:** The mirror symmetry discussed here has been
> proven analytically. Visibility windows relate to palindromic mode pairs
> under Π. See [MIRROR_SYMMETRY_PROOF.md](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).

---

## The Question

The localizable entanglement benchmark showed that CoA (the total entanglement resource)
is nearly constant (~1) while CΨ flashes briefly. LE sits between them. The resource is
always there, but CΨ only sees it in certain moments.

The next question: if two observers have different local time rates (different γ), do they
see *different* CΨ windows? Does the formula "rotate" per observer?

## Setup

Star topology (S + A + B), Bell_SA x |+⟩_B, J_SA=1.0, J_SB=2.0. Four runs with
different γ_A values (γ_S = γ_B = 0.05 held constant):

| Config | γ_A | γ_B | γ_S |
|--------|-----|-----|-----|
| slow A | 0.03 | 0.05 | 0.05 |
| equal | 0.05 | 0.05 | 0.05 |
| fast A | 0.10 | 0.05 | 0.05 |
| very fast A | 0.20 | 0.05 | 0.05 |

## Key Result

At the same coordinate time t, observers with different γ_A see different CΨ_AB:

| t | slow A (γ=0.03) | equal (γ=0.05) | fast A (γ=0.10) | very fast A (γ=0.20) |
|---|---|---|---|---|
| 0.2 | 0.298 | 0.293 | 0.281 | 0.259 |
| 0.4 | **0.340** | **0.329** | **0.303** | **0.257** |
| 1.0 | 0.157 | 0.141 | 0.108 | 0.060 |
| 2.0 | 0.037 | 0.021 | 0.000 | 0.000 |

The faster the observer's time runs (higher γ_A), the **less** they see of the connection.

## The Critical Threshold

Peak CΨ_AB values per observer:

| Config | Peak CΨ_AB | Above 1/4? |
|--------|-----------|------------|
| slow A (γ=0.03) | 0.340 | Yes, clearly |
| equal (γ=0.05) | 0.329 | Yes |
| fast A (γ=0.10) | 0.303 | Yes, barely |
| very fast A (γ=0.20) | 0.261 | Barely |

A slightly noisier observer would see the peak drop *below* 1/4. The connection
would never cross the threshold. The same entanglement, the same global state,
the same CoA ~ 1 resource - but one observer sees the crossing, the other does not.

## Non-Universality in Proper Time

If CΨ_AB were a function of proper time τ_A = γ_A * t alone, then plotting
CΨ_AB(τ_A) should give the same curve for all observers. It does not:

| τ_A | slow A | equal | fast A | very fast A |
|-----|--------|-------|--------|-------------|
| 0.005 | 0.237 | 0.147 | 0.059 | --- |
| 0.010 | 0.308 | 0.293 | 0.143 | 0.058 |
| 0.020 | 0.095 | 0.329 | 0.281 | 0.137 |
| 0.050 | 0.131 | 0.141 | 0.169 | 0.261 |

At the same proper time, different observers see completely different values.
CΨ_AB is not a function of τ_A alone. It depends on the *absolute* noise rate,
not just the elapsed proper time.

This means: the visibility of the connection depends on **who you are** (your γ),
not just **where you are** in your own timeline.

## What This Means

Three layers of the same initial entanglement setup:

| Metric | What it shows | Depends on noise profile? |
|--------|--------------|--------------------------|
| CoA ~ 1 | The resource exists in the global state | Weakly |
| LE ~ 0.4 | The resource could be extracted by optimal measurement | Moderately |
| CΨ flashes | The resource is currently expressed as coherent pairwise structure | **Strongly** |

CoA and LE are relatively robust to changes in local noise.
CΨ is highly sensitive: different noise profiles open or close the visibility window.

Note: this is not "same state, different observers" in the strict sense. Different
γ_A values produce different Lindblad trajectories. The correct framing is:
different noise environments produce qualitatively different direct-visibility
windows for an entanglement resource that persists under assistance.

## Important Correction (from external review)

When γ_A is changed, the Lindblad generator changes. At fixed coordinate time t,
the four runs are **not** comparing the same evolved global state seen by different
observers. They are comparing **different open-system trajectories** from the same
initial state and topology.

This distinction matters: "same state, different observers" would be a much stronger
claim than "different noise profiles, different trajectories." The correct statement is:

**Different local noise profiles produce different direct-visibility windows for a
pairwise resource that remains available under assistance.**

This is narrower than "different observers see the same reality differently" but
technically stronger and defensible.

## Why No Proper-Time Universality

The non-collapse onto a universal curve in τ_A = γ_A * t is not mysterious. In
the rescaled equation, the dimensionless ratios J/γ_A, γ_B/γ_A, and γ_S/γ_A all
change when γ_A changes. The dynamics have **multiple competing timescales**
(Hamiltonian coupling, three separate noise rates), not a single local decoherence
clock. Universal collapse only occurs when the entire generator is controlled by
one timescale, which is not the case here.

The cleanest way to state this:

**"The absence of τ_A-collapse shows that CΨ windows are controlled by multi-scale
open-system dynamics, not by a single local decoherence clock."**

## What CΨ Witnesses (Precise Statement)

Based on the LE benchmark and this visibility test, the cleanest characterization:

**CΨ is a basis-fixed, unassisted witness of directly expressed pairwise entanglement.**

The three-layer version:

- **CoA** quantifies assisted entanglement capacity (observer-independent)
- **LE** quantifies localizable entanglement potential (observer-independent)
- **CΨ** quantifies unassisted, basis-fixed direct expression in the observed marginal

One sentence that survives hostile reading:

**"The same initial network can retain high assisted entanglement capacity while
local noise profiles open or close windows in which that resource is directly
visible in a chosen pair and basis."**



This result connects the three key findings identified in the external review:

1. **Act vs process (99% vs 69%):** Different *types* of intervention produce
   different visibility. Now we see that different *rates of noise* also produce
   different visibility. Both are about the channel structure, not just the resource.

2. **Echoes:** When SA and SB are zero but AB is nonzero, the resource is still
   there (CoA ~ 1). CΨ sees a faint residual. A slower observer (lower γ) would
   see a *stronger* residual of the same echo.

3. **Three conditions:** "Quiet receiver" (low γ_A) is now understood more precisely:
   it means the receiver's time runs slowly enough that the CΨ window stays open
   long enough to see the connection. A fast receiver's window closes before the
   connection becomes visible.

## The Technical Thread

The external reviewer (GPT) identified CΨ as a "basis-dependent witness for
coherent, localizable pairwise entanglement." After this test and the correction,
the refined characterization:

**CΨ is a basis-fixed, unassisted witness of directly expressed pairwise
correlated coherence/entanglement.**

It sits at the intersection of two known concepts:
- **Localizable entanglement / entanglement of assistance** (resource available in principle)
- **Quantum coherence / correlated coherence** (basis-dependent structure already present)

CΨ is not redundant with either. It is more selective than both: it requires
the entanglement to be *already directly expressed* in the reduced pair and basis,
without optimization or intervention.

No standard quantum-information term exists for exactly this quantity. The closest
neighborhoods are localizable entanglement, correlated coherence, and (distantly)
hidden entanglement - but CΨ = concurrence x normalized l1-coherence is its own
diagnostic at the intersection of these concepts.

**Warning:** The phrase "observer-dependent entanglement" exists in relativistic
quantum information (dependence from motion, gravity, mode accessibility). That is
a different phenomenon from varying a local dephasing rate. This document avoids
that terminology.

---

## Simulation Code

The proper time intersection test is at `simulations/proper_time_intersection_test.py`.

## References

- [Localizable Entanglement Benchmark](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md) - CΨ vs LE vs CoA
- [Star Topology](STAR_TOPOLOGY_OBSERVERS.md) - Three conditions, receiver noise
- [THE_CPSI_LENS](../docs/THE_CPSI_LENS.md) - Canonical description
- Verstraete, Popp, Cirac, *Entanglement versus Correlations in Spin Systems* (2004), arXiv:quant-ph/0411123
- Streltsov, Adesso, Plenio, *Colloquium: Quantum coherence as a resource* (2017), Rev. Mod. Phys. 89, 041003
- Alsina, Razavi, *Observer dependent entanglement* (2012), arXiv:1210.2223 (note: different phenomenon - relativistic, not noise-rate-dependent)
