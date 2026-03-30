# Combined Optimization: What We Know and What We Don't

<!-- Keywords: sacrifice zone selective dynamic decoupling combined,
chain selection mode protection IBM Torino, cavity mode localization
intra-chain inter-chain comparison, April 9 hardware test experiment plan,
R=CPsi2 combined optimization -->

**Status:** Tier 2-3 (eigenvalue analysis proven, time evolution
simulated, hardware interpretation open). Hardware test planned April 9.
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)
**Scripts:** [combined_optimization.py](../simulations/combined_optimization.py),
[time_evolution_6scenarios.py](../simulations/time_evolution_6scenarios.py),
[time_evolution_neel.py](../simulations/time_evolution_neel.py)

---

## Three levels of analysis (do not confuse them)

This document separates three questions that are easy to conflate:

### Level 1: Eigenvalue efficiency (proven)

*Within a fixed total noise budget, how efficiently does the noise
distribution protect the slowest cavity modes?*

| Configuration | Σγ | vs own uniform |
|--------------|-----|---------------|
| Sacrifice, Selective DD | 0.290 | **3.72x** |
| Sacrifice, Uniform DD | 0.214 | 3.35x |
| Sacrifice, No DD | 0.320 | 2.81x |
| Mean-T2, Selective DD | 0.033 | 1.25x |
| Mean-T2, Uniform DD | 0.024 | 1.06x |
| Mean-T2, No DD | 0.048 | 1.06x |

**vs own uniform** compares each configuration against the same total
noise spread equally. This is a fair, apples-to-apples comparison.
The sacrifice chain's noise gradient protects modes 3.72x better
than equal spreading. The mean-T2 chain's flat profile provides
almost no spatial advantage (1.06x).

This is mathematically sound. The cavity mode localization (r = 0.994)
explains the mechanism: center-localized modes see less edge noise.

### Level 2: Intra-chain DD comparison (hardware-confirmed)

*On the SAME chain, does selective DD beat uniform DD?*

This is what the March 24 hardware experiment tested on Q85-Q94:

| DD strategy | Avg SumMI | vs Uniform DD |
|------------|----------|--------------|
| Selective DD | 0.054 | **2.02x** |
| No DD | 0.045 | 1.71x |
| Uniform DD | 0.027 | 1.00x |

**Hardware-confirmed at all 5 time points.** Selective DD outperforms
uniform DD by 2.0x average, up to 3.2x at t=4.0. This is on a
single chain with |+>^5 initial state.

Open question: Is this because DD on Q85 wastes gates (Interpretation A)
or because noise contrast helps modes (Interpretation B)? See
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md) for the A/B discussion.

### Level 3: Inter-chain comparison (simulation only, depends on initial state)

*Does a sacrifice chain with more total noise outperform a quiet chain
with less total noise in absolute SumMI?*

This depends entirely on the initial state:

**With |+>^5 (Heisenberg eigenstate):**
The quiet chain produces near-zero SumMI (eigenstate frozen, no dynamics).
The sacrifice chain produces SumMI ~ 0.05-0.20 (noise breaks symmetry).
The sacrifice chain "wins" by 200x, but this is an artifact: |+>^5 does
not evolve without noise.

**With |01010> (Neel state, fair comparison):**
Both chains show rich dynamics. The quiet chain wins in absolute SumMI:

| t=2.5 | SumMI | Purity |
|-------|-------|--------|
| Mean-T2, Uniform DD | **1.650** | 0.907 |
| Mean-T2, No DD | 1.510 | 0.824 |
| Sacrifice, Selective DD | 0.970 | 0.386 |
| Sacrifice, No DD | 0.908 | 0.351 |

**Less total noise wins.** The sacrifice chain's 3.72x eigenvalue
efficiency does not compensate for 6x more total noise. This is
not surprising: a chain with Σγ = 0.03 decays slower than one with
Σγ = 0.29, regardless of how cleverly the noise is distributed.

---

## What is established

1. **Eigenvalue efficiency of spatial noise gradients** (Tier 1-2):
   The sacrifice-zone formula produces 2.81-3.72x mode protection
   within a fixed noise budget. Verified N=2-7.

2. **Intra-chain selective DD advantage** (Tier 2, hardware):
   Selective DD beats uniform DD by 2.0x on the sacrifice chain.
   Single run, caveats apply.

3. **Center-localized modes are protected** (Tier 2, computed):
   r = 0.994 correlation between edge weight and decay rate.
   Palindromic partners are spatially complementary.

4. **The formula is correct** (Tier 1, algebraic):
   Stationary(N) = Sum_J m(J,N)*(2J+1)^2. Verified N=2-7.

## What is NOT established

1. **That sacrifice chains beat quiet chains in absolute performance.**
   They don't, with the Neel state. The sacrifice-zone is about
   EFFICIENCY, not absolute advantage.

2. **Whether the hardware advantage is mode protection (B) or gate-error
   avoidance (A).** The A/B test has not been performed.

3. **Whether the results hold for initial states other than |+>^5.**
   The March 24 hardware data used |+>^5. With Neel state, the
   dynamics are fundamentally different.

4. **Whether the eigenvalue efficiency translates to observable advantage
   at fixed total noise.** This requires comparing two chains with
   SIMILAR Σγ but different spatial profiles.

---

## IBM experiment plan (April 9, 2026)

### The right experiment

The most informative test is NOT "sacrifice vs mean-T2." Those chains
have 6x different total noise. The informative tests are:

**Test 1: A/B test (Interpretation A vs B)**
On a chain with GOOD qubits everywhere (all T2 > 100 us), apply
selective DD (remove DD from one edge qubit). If selective DD still
beats uniform DD, it is the noise CONTRAST that helps (B), not
gate-error avoidance (A). If selective DD loses, A was correct.

**Test 2: Intra-chain replication**
Reproduce the March 24 result on chain Q85-Q94. Same chain, same DD
strategies, but with Neel initial state |01010>. Does the 2.0x
advantage persist with a different initial state?

**Test 3: Equal-noise comparison**
Find two chains with SIMILAR Σγ but different spatial profiles. One
with an edge sacrifice (high gradient), one with uniform noise. Compare
SumMI(t). This isolates the mode-protection effect from the total-noise
effect.

### Practical setup

| Parameter | Value |
|-----------|-------|
| Chain A | Q85-Q86-Q87-Q88-Q94 (sacrifice, Σγ ~ 0.32) |
| Chain B | To be selected: similar Σγ, uniform profile |
| Chain C | Best mean-T2 chain (for A/B test with selective DD) |
| Initial states | |+>^5 (for March 24 replication) AND |01010> |
| DD configs | None, Uniform, Selective (per chain) |
| Trotter | dt = 0.5 us, steps = [2, 4, 6, 8, 10] |
| Shots | 4000 per config per time point |

### Success criteria

| Test | Success | Failure |
|------|---------|---------|
| A/B test | Selective DD wins on good chain (B confirmed) | Selective DD loses on good chain (A confirmed) |
| Replication | Selective > Uniform at 3+ time points | No consistent advantage |
| Equal-noise | Sacrifice profile > uniform profile at same Σγ | No advantage at equal noise |

**Both A and B outcomes are honest results.** The theory predicts B
(mode localization). If A is correct (gate-error avoidance), the
cavity-mode explanation is incomplete and needs revision.

---

*See also:*
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md) (1.97x measured March 24),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994),
[Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md) (chain selection),
[Resonant Return](RESONANT_RETURN.md) (the sacrifice-zone formula)
