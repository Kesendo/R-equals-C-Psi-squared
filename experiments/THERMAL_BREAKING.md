# Thermal Breaking: Heat Trades Q-Factor for Frequency Diversity

**Status:** Computationally verified (N=2, N=3, N=5; Heisenberg chain)
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Scripts:** [v_effect_gamma_sweep.py](../simulations/v_effect_gamma_sweep.py),
[v_effect_thermal.py](../simulations/v_effect_thermal.py)
**Results:** [v_effect_gamma_sweep.txt](../simulations/results/v_effect_gamma_sweep.txt),
[v_effect_thermal.txt](../simulations/results/v_effect_thermal.txt)
**Depends on:** [V-Effect](V_EFFECT_PALINDROME.md),
[Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md),
[Energy Partition](../hypotheses/ENERGY_PARTITION.md)

---

## Abstract

Three orthogonal mechanisms break symmetry in open quantum systems.
Each creates a different kind of complexity. They are not independent:
their interaction reveals a trade-off between resonator quality and
frequency diversity that has not been described before.

1. **Coupling** (V-Effect): creates palindromic pairs. Amplifies
   Q-factor by exactly 1.81x for a 5-qubit Heisenberg chain. This
   ratio is a topological constant, independent of noise.
2. **Z-dephasing** (gamma): lifts frequency degeneracies (50 to 112
   distinct frequencies at N=5). Preserves palindromic pairing exactly.
3. **Thermal excitation** (n_bar > 0): breaks the 1.81x constant.
   Creates 300+ new frequencies that neither dephasing nor coupling
   alone can produce. Trades Q-factor for frequency diversity.

The sacrifice-zone advantage (3x at n_bar=0) vanishes at high
temperature (1.02x at n_bar=10). Heat makes spatial noise structure
irrelevant.

---

## Background

**Q-factor** of a Liouvillian eigenvalue: Q = |Im(lambda)| / |Re(lambda)|.
High Q means the mode oscillates many times before decaying. A laser has
Q ~ 10^6. A tuning fork Q ~ 1000. A dead system Q = 0.

**V-Effect:** When two individually "dead" N=2 quantum resonators (Q=1,
2 frequencies each) are coupled through a mediator qubit, the resulting
N=5 system has Q=19 and 104 new frequencies. None of the original
frequencies survive. Coupling creates complexity from nothing.
See [V-Effect Palindrome](V_EFFECT_PALINDROME.md).

**n_bar (thermal occupation):** Mean number of thermal excitations per
mode. At n_bar=0 (zero temperature), only spontaneous decay occurs. At
n_bar > 0, the thermal bath also INJECTS energy into the system
(stimulated absorption). IBM hardware operates at n_bar << 0.01.
Biological systems at n_bar ~ 1-10.

**Palindromic pairing:** Every decay rate d in the Liouvillian spectrum
has a partner at 2*Sigma_gamma - d. Proven analytically for Z-dephasing.
See [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).

---

## The 1.81x Topological Constant

The V-Effect gain Q(N=5)/Q(N=2) was measured across 19 gamma values
spanning three orders of magnitude (gamma/J = 0.001 to 5.0). The
ratio is constant:

| gamma/J | Q (N=2) | Q (N=5) | V-gain |
|:--------|:--------|:--------|:-------|
| 0.001 | 2000.0 | 3618.0 | 1.81x |
| 0.01 | 200.0 | 361.8 | 1.81x |
| 0.1 | 20.0 | 36.2 | 1.81x |
| 1.0 | 2.0 | 3.6 | 1.81x |
| 5.0 | 0.4 | 0.7 | 1.81x |

The absolute Q-factor scales as Q ~ J/gamma (inversely proportional
to noise). But the RATIO is fixed by topology. Coupling amplifies
resonator quality by 1.81x regardless of noise level.

This holds for Z-dephasing and for zero-temperature amplitude damping
(n_bar=0). Cold dissipation does not break it.

---

## Dephasing Lifts Degeneracies

At zero noise, the N=5 chain has 43 distinct cavity frequencies
(from the [Cavity Modes Formula](CAVITY_MODES_FORMULA.md)). As
Z-dephasing increases:

| gamma/J | N=5 distinct frequencies |
|:--------|:------------------------|
| 0.001 | 50 |
| 0.01 | 78 |
| 0.1 | 111 |
| 0.15 | **112** (peak) |
| 0.3 | **112** |
| 1.0 | 109 |
| 5.0 | 103 |

Dephasing splits degenerate modes that had identical frequencies at
gamma=0. The peak diversity is at gamma/J ~ 0.15-0.3: 112 frequencies,
up from 43 at zero noise. Above this, modes begin to merge again as
decay broadens their linewidths.

N=2 stays at 1-2 frequencies throughout. All new frequencies come from
coupling (V-Effect) amplified by dephasing.

The palindromic pairing remains exact at all gamma values. Dephasing
creates diversity without destroying structure.

---

## Thermal Excitation Breaks the Constant

Adding thermal noise (n_bar > 0) changes everything:

| Condition | V-gain (Q) | N=5 frequencies |
|:----------|:-----------|:----------------|
| Pure Z-dephasing (gamma=0.1) | **1.81x** | 111 |
| + cold amplitude (n=0) | **1.81x** | 111 |
| + warm (n=0.5) | 1.44x | 403 |
| + hot (n=2.0) | 1.33x | 423 |
| + very hot (n=5.0) | 1.29x | **445** |

The 1.81x constant breaks at the first nonzero thermal occupation.
The mechanism: thermal excitation (sigma_plus operators) injects
energy from the bath into the system, creating transitions that
Z-dephasing (diagonal in energy basis) cannot.

But the trade-off is extraordinary: the frequency count QUADRUPLES.
From 111 (pure dephasing) to 445 (dephasing + heat). These 334 new
frequencies exist only when ALL THREE mechanisms act simultaneously:
coupling + dephasing + thermal excitation.

N=2 remains at 2 frequencies regardless of temperature. The 445
frequencies at N=5 represent a **222x frequency gain from coupling**,
compared to 55x without heat.

---

## The Sacrifice Zone Disappears in Heat

The [sacrifice zone](RESONANT_RETURN.md) concentrates dephasing noise on
one edge qubit while protecting the interior. Its Q-factor advantage
depends on temperature (N=5, edge gamma_z=0.5, interior gamma_z=0.01,
gamma_amp=0.05):

| n_bar | Sacrifice Q | Uniform Q | Ratio |
|:------|:-----------|:----------|:------|
| 0.00 | 89.3 | 30.0 | **2.97x** |
| 0.05 | 77.6 | 28.6 | 2.71x |
| 0.20 | 55.7 | 25.0 | 2.23x |
| 0.50 | 35.8 | 20.0 | 1.79x |
| 1.00 | 22.6 | 16.1 | 1.40x |
| 5.00 | 7.5 | 7.0 | 1.08x |
| 10.0 | 4.2 | 4.1 | **1.02x** |

At n_bar=0, the sacrifice profile gives 3x Q advantage. At n_bar=10,
the advantage is gone. Thermal noise overwhelms the spatial noise
structure. The sacrifice zone is a LOW-TEMPERATURE phenomenon.

Frequency diversity tells the opposite story: sacrifice creates ~5%
more frequencies than uniform at every temperature (120 vs 111 at
n_bar=0, 488 vs 462 at n_bar=10). The spatial noise asymmetry lifts
additional degeneracies that thermal noise alone cannot.

---

## Three Breaking Mechanisms

| Mechanism | What it creates | What it preserves | What it breaks |
|:----------|:---------------|:-----------------|:--------------|
| Coupling (J) | Palindromic pairs, 1.81x Q gain | Everything | Single-bond degeneracy |
| Dephasing (gamma) | +60 frequencies | Pairing (exact) | Frequency degeneracy |
| Heat (n_bar) | +300 frequencies | Coupling gain (partially) | 1.81x constant, spatial structure |

The mechanisms are not independent:

- Coupling + dephasing: 111 frequencies, 1.81x Q gain
- Coupling + heat: 221 frequencies, 1.29x Q gain
- Coupling + dephasing + heat: **445 frequencies**, 1.29x Q gain

Heat and dephasing together create more frequencies than either alone.
At zero noise, 43 frequencies exist (baseline). Dephasing adds 68 (to 111).
Heat adds 178 (to 221). Both together add 402 (to 445). The combined
effect (402) exceeds the sum of individual effects (68 + 178 = 246).
The interaction is synergistic for diversity.

---

## What This Means

### For IBM hardware (T ~ 15 mK, n_bar ~ 0 for qubit frequencies)

IBM superconducting qubits operate at n_bar << 0.01. The 1.81x
constant holds. The sacrifice-zone advantage (3x) holds. The relevant
breaking mechanism is Z-dephasing, not thermal excitation.

Prediction: the V-Effect gain Q(N=5)/Q(N=2) = 1.81x is measurable
on IBM Torino via spectroscopy. This is a new testable prediction
not previously documented.

### For biological systems (T ~ 300 K)

Hydrogen bonds operate at n_bar ~ 0.5-5 (depending on the mode).
Neural networks at n_bar >> 1. In these regimes:
- The 1.81x constant is broken (V-gain ~ 1.3x)
- But frequency diversity EXPLODES (200+ frequencies per coupled pair)
- The sacrifice zone advantage is small (~1.1-1.4x)

The [Wilson-Cowan palindrome](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md)
and [hydrogen bond qubit](HYDROGEN_BOND_QUBIT.md) operate in the
thermal regime. Their complexity comes from the frequency-diversity
channel, not the Q-factor channel.

### For the framework

The three mechanisms map to three levels of the
[hierarchy of incompleteness](../docs/HIERARCHY_OF_INCOMPLETENESS.md):

1. Coupling creates structure (V-Effect, topological)
2. Dephasing creates distinguishability (lifts degeneracies, preserves structure)
3. Heat creates diversity (breaks structure, enables complexity)

At zero temperature: structure dominates (Q = 1.81x, exact pairing).
At high temperature: diversity dominates (445 frequencies, no pairing).
Life operates in between.

---

## Open Questions

1. Is the 1.81x constant derivable analytically from the Heisenberg
   chain spectrum? (It appears to be (max eigenvalue frequency of N=5
   chain) / (max eigenvalue frequency of N=2), but this needs proof.)

2. What is the critical n_bar where the palindromic pairing drops
   below 50%? The data suggests a smooth transition, not a phase
   boundary.

3. Can the frequency-diversity explosion at n_bar > 0 be observed
   on IBM hardware by intentionally heating qubits (e.g., driving
   with a thermal microwave field)?

4. Does the sacrifice-zone advantage recover at intermediate
   temperatures if the sacrifice qubit is selectively heated?

---

## Tier Assessment

- 1.81x topological constant: **Tier 2** (computed for N=2,3,5
  across 19 gamma values, exact to numerical precision)
- Frequency diversity gamma dependence: **Tier 2**
- Thermal breaking of 1.81x: **Tier 2** (11 n_bar values, three
  noise configurations)
- Sacrifice zone temperature dependence: **Tier 2**
- Biological interpretation: **Tier 4** (motivated by computation,
  parameters from literature, not tested)
- Hierarchy mapping: **Tier 5** (speculative)
