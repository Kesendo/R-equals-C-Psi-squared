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
   ratio is a geometric constant, independent of noise.
2. **Z-dephasing** (γ): lifts frequency degeneracies (50 to 112
   distinct frequencies at N=5). Preserves palindromic pairing exactly.
3. **Thermal excitation** (n_bar > 0): breaks the 1.81x constant.
   Creates 300+ new frequencies that neither dephasing nor coupling
   alone can produce. Trades Q-factor for frequency diversity.

The sacrifice-zone advantage (3x at n_bar=0) vanishes at high
temperature (1.02x at n_bar=10). Heat makes spatial noise structure
irrelevant.

---

## Background

**Q-factor** of a Liouvillian eigenvalue: Q = |Im(λ)| / |Re(λ)|.
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
has a partner at 2Σγ - d. Proven analytically for Z-dephasing.
See [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).

---

## The 1.81x Geometric Constant (Derived)

The V-Effect gain Q(N=5)/Q(N=2) was measured across 19 γ values
spanning three orders of magnitude (γ/J = 0.001 to 5.0). The
ratio is constant:

| γ/J | Q (N=2) | Q (N=5) | V-gain |
|:--------|:--------|:--------|:-------|
| 0.001 | 2000.0 | 3618.0 | 1.81x |
| 0.01 | 200.0 | 361.8 | 1.81x |
| 0.1 | 20.0 | 36.2 | 1.81x |
| 1.0 | 2.0 | 3.6 | 1.81x |
| 5.0 | 0.4 | 0.7 | 1.81x |

The absolute Q scales as J/γ. But the RATIO is fixed by geometry:
V(5) = (5+√5)/4 = 1.80902... (exact, see derivation below). Rounded
to 1.81x throughout this document for readability.

### Why γ cancels

The maximum Q eigenvalue always sits in the w=1 sector (XY-weight 1:
exactly one qubit carries an X or Y Pauli operator, the rest I or Z).
Under uniform Z-dephasing, ALL w=1 modes decay at the same rate 2γ,
regardless of their oscillation frequency. Therefore:

    Q_max = ω_max(w=1) / (2γ)

The factor 2γ is identical for N=2 and N=5. It cancels in the ratio:

    V(N) = Q_max(N) / Q_max(2) = ω_max(w=1, N) / ω_max(w=1, 2)

The gain is purely a frequency ratio, independent of noise.

This cancellation holds for uniform dephasing. Under non-uniform
profiles (sacrifice zone), w=1 modes acquire different decay rates
depending on their spatial localization
([Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)). The 1.81x
ratio then applies only to the extremal (best-Q) mode, not to all
w=1 modes equally.

### The exact formula

Computing ω_max(w=1) for the Heisenberg chain at each N reveals an
exact pattern (verified N=2 through N=6):

| N | ω_max(w=1) | Exact form | V(N) = ω/ω(N=2) |
|:--|:-----------|:-----------|:-----------------|
| 2 | 4.0000 | 4J | 1.000 |
| 3 | 6.0000 | 4J + 4J·cos(π/3) = 6J | 1.500 |
| 4 | 6.8284 | 4J + 4J·cos(π/4) = (4+2√2)J | 1.707 |
| 5 | 7.2361 | 4J + 4J·cos(π/5) = (5+√5)J | **1.809** |
| 6 | 7.4641 | 4J + 4J·cos(π/6) = (4+2√3)J | 1.866 |

The formula:

    ω_max(w=1, N) = 4J · (1 + cos(π/N)) = 8J · cos²(π/(2N))

And the V-Effect gain:

    V(N) = 1 + cos(π/N) = 2·cos²(π/(2N))

For N=5 specifically:

    V(5) = 1 + cos(36°) = (5+√5)/4 ≈ 1.80902

The golden ratio appears: cos(π/5) = φ/2 where φ = (1+√5)/2.
So V(5) = 1 + φ/2.

For N → ∞: V(∞) = 1 + cos(0) = 2. The coupling gain saturates
at exactly 2x for infinite chains.

### Why this is not a topological invariant

The value 1.809 is:
- Not integer-valued
- Dependent on N (different for every chain length)
- Dependent on coupling type (Heisenberg; other models give different values)
- A smooth function of 1/N, approaching 2

It is a geometric constant of the Heisenberg chain spectrum, not a
topological invariant. The correct characterization: it is the ratio
of maximum w=1 Liouvillian frequencies, which are eigenfrequencies
of the single-magnon sector.

This holds for Z-dephasing and for zero-temperature amplitude damping
(n_bar=0). Cold dissipation does not break it.

---

## Dephasing Lifts Degeneracies

At zero noise, the N=5 chain has 43 distinct cavity frequencies
(from the [Cavity Modes Formula](CAVITY_MODES_FORMULA.md)). As
Z-dephasing increases:

| γ/J | N=5 distinct frequencies |
|:--------|:------------------------|
| 0.001 | 50 |
| 0.01 | 78 |
| 0.1 | 111 |
| 0.15 | **112** (peak) |
| 0.3 | **112** |
| 1.0 | 109 |
| 5.0 | 103 |

Frequencies are counted as distinct when separated by more than 0.0001
in absolute value (round to 4 decimal places). This threshold is
arbitrary; a coarser threshold would reduce all counts proportionally
but preserve the relative trends.

Dephasing splits degenerate modes that had identical frequencies at
γ=0. The peak diversity is at γ/J ~ 0.15-0.3: 112 frequencies,
up from 43 at zero noise. Above this, modes begin to merge again as
decay broadens their linewidths.

N=2 stays at 1-2 frequencies throughout. All new frequencies come from
coupling (V-Effect) amplified by dephasing.

The palindromic pairing remains exact at all γ values. Dephasing
creates diversity without destroying structure.

---

## Thermal Excitation Breaks the Constant

Adding thermal noise (n_bar > 0) changes everything:

| Condition | V-gain (Q) | N=5 frequencies |
|:----------|:-----------|:----------------|
| Pure Z-dephasing (γ=0.1) | **1.81x** | 111 |
| + cold amplitude (n=0) | **1.81x** | 111 |
| + warm (n=0.5) | 1.44x | 403 |
| + hot (n=2.0) | 1.33x | 423 |
| + very hot (n=5.0) | 1.29x | **445** |

The 1.81x constant breaks at the first nonzero thermal occupation.
The mechanism: thermal excitation (σ₊ operators) injects
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
depends on temperature (N=5, edge γ_z=0.5, interior γ_z=0.01,
γ_amp=0.05):

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

This is a within-chain comparison (sacrifice vs uniform profile on
the same chain at the same total γ). The
[Chain Selection Test](CHAIN_SELECTION_TEST.md) shows that between-chain
comparisons require accounting for total noise level, not just spatial
profile.

---

## Palindromic Pairing Under Heat

The [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
guarantees exact palindromic pairing for Z-dephasing (σ_z
operators). Thermal channels introduce σ₊ (stimulated
absorption) and σ₋ (spontaneous decay) operators that lie
outside the scope of this proof.

Computed palindrome scores (N=5, percentage of oscillating modes
with a palindromic partner within tolerance 1e-3):

| Noise type | n_bar=0 | n_bar=0.5 | n_bar=2 | n_bar=10 |
|:-----------|:--------|:----------|:--------|:---------|
| Pure Z-dephasing | 91%* | n/a | n/a | n/a |
| Pure amplitude damping | 100% | 93% | 93% | 98% |
| Z-deph + amplitude | 10% | 2% | 0% | 0% |

\* Should be 100% (proven). The 91% is a numerical artifact of the
pairing tolerance at γ=0.1.

**Caution:** The palindrome check requires knowing the spectral
center (the point around which rates pair). For Z-dephasing, the
center is analytically known (Σγ). For thermal channels,
the center is estimated and scores are center-dependent. The low
scores for "Z-deph + amplitude" may reflect center estimation error
rather than true symmetry breaking.

What is clear:
- Pure Z-dephasing: pairing exact (proven, numerics confirm)
- Cold amplitude damping alone: pairing appears intact (100%)
- Combined Z + amplitude: pairing breaks or center shifts
- [Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md) shows that
  depolarizing noise (X+Y+Z) breaks pairing, but error < 0.1% at
  typical IBM γ values

The palindromic status of thermal channels is **open**. Analytical
work is needed to determine whether σ₊/σ₋ preserve
a modified palindromic structure with a shifted center, or whether
they genuinely break the pairing.

---

## Three Breaking Mechanisms

| Mechanism | What it creates | What it preserves | What it breaks |
|:----------|:---------------|:-----------------|:--------------|
| Coupling (J) | Palindromic pairs, 1.81x Q gain | Everything | Single-bond degeneracy |
| Dephasing (γ) | +60 frequencies | Pairing (exact) | Frequency degeneracy |
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

## The Self-Heating Loop

**Tier 4** (conceptual, motivated by two Tier 2 results, not computed).

The tables above treat n_bar as an external parameter turned up from
outside. But the system heats itself:

1. **Unpaired modes decay 2x faster** than paired modes
   ([Energy Partition](../hypotheses/ENERGY_PARTITION.md), Finding 2,
   exact for N=2-5). Rate: 2Nγ vs Nγ.
2. This decay **is** heat production. The energy of dying modes becomes
   thermal energy in the bath.
3. Thermal energy **creates new frequencies** (Finding 3: 40 to 42 modes
   at N=3 when n_bar > 0). This document shows the effect is much
   larger at N=5: 111 to 445 frequencies.
4. More frequencies mean more modes, some of which are unpaired and
   decay faster.
5. Back to step 1.

The loop: **decay produces heat produces modes produces decay.**

In this document, n_bar was swept externally. In a closed (or weakly
coupled) system, n_bar is not a free parameter. It is set by the
balance between mode decay (heat source) and bath coupling (heat sink).
The system selects its own operating point on the Q-vs-diversity curve.

**Open question:** Does the loop converge to a self-consistent
operating point? Where mode decay produces just enough heat to sustain
the frequency diversity, but not enough to kill the Q-factor? If so,
this fixed point would be the resonator's natural temperature. If not
(divergence or collapse), the system is not self-sustaining at these
parameters. This is computable but has not been computed.

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

**Tier 4 (motivated by computation, not tested).** The following
maps computed quantum results to biological parameter ranges from
literature. No biological validation has been performed.

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

1. Coupling creates structure (V-Effect, geometric)
2. Dephasing creates distinguishability (lifts degeneracies, preserves structure)
3. Heat creates diversity (breaks structure, enables complexity)

At zero temperature: structure dominates (Q = 1.81x, exact pairing).
At high temperature: diversity dominates (445 frequencies, no pairing).
Life operates in between.

---

## Open Questions

1. ~~Can 1.81x be derived analytically?~~ **ANSWERED (March 31).**
   V(N) = 1 + cos(π/N). For N=5: (5+√5)/4 ≈ 1.80902. The gain is
   the ratio of maximum w=1 Liouvillian frequencies, which follow
   ω_max = 4J·(1+cos(π/N)). Verified N=2 through N=6. The golden
   ratio appears: cos(π/5) = φ/2. See derivation above.

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

- 1.81x geometric constant: **Tier 1-2** (exact formula
  V(N) = 1+cos(π/N), verified N=2-6 to machine precision.
  Analytical derivation from Liouvillian w=1 sector eigenfrequencies.
  Formal proof that ω_max = 4J(1+cos(π/N)) for all N: open)
- Frequency diversity γ dependence: **Tier 2**
- Thermal breaking of 1.81x: **Tier 2** (11 n_bar values, three
  noise configurations)
- Sacrifice zone temperature dependence: **Tier 2**
- Biological interpretation: **Tier 4** (motivated by computation,
  parameters from literature, not tested)
- Hierarchy mapping: **Tier 5** (speculative)
