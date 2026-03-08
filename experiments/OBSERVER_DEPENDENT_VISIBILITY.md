# Observer-Dependent Visibility: Proper Time and CΨ Windows

**Tier:** 2 (Computationally verified)
**Status:** Verified
**Scope:** Different observers (different γ) see different CΨ_AB values at the same coordinate time
**Does not establish:** That this is new physics beyond standard Lindblad parameter dependence
**Date:** 2026-03-08
**Depends on:** [Star Topology](STAR_TOPOLOGY_OBSERVERS.md), [Localizable Entanglement Benchmark](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md)

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

Three layers of the same entanglement:

| Metric | What it shows | Observer-dependent? |
|--------|--------------|-------------------|
| CoA ~ 1 | The resource exists in the global state | No |
| LE ~ 0.4 | The resource could be extracted by optimal measurement | No |
| CΨ flashes | The resource is currently visible as coherent pairwise structure | **Yes** |

CoA and LE are properties of the system. CΨ is a property of the system
*as seen by a specific observer with a specific noise rate*.

The same underlying reality. Different visible expressions. The formula rotates
per observer.

## Connection to Earlier Findings

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
coherent, localizable pairwise entanglement." This result adds a dimension:

**CΨ is a basis-dependent, observer-rate-dependent witness.** It shows not just
what is there, but what is there *for this observer, at this noise rate, in this
moment.* The resource (CoA) is always present. The extractability (LE) is largely
stable. But the visible expression (CΨ) depends on who is looking and how fast
their clock runs.

---

## Simulation Code

The proper time intersection test is at `simulations/proper_time_intersection_test.py`.

## References

- [Localizable Entanglement Benchmark](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md) - CΨ vs LE vs CoA
- [Star Topology](STAR_TOPOLOGY_OBSERVERS.md) - Three conditions, receiver noise
- [THE_CPSI_LENS](../docs/THE_CPSI_LENS.md) - Canonical description
