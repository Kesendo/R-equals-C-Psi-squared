# Thermal Breaking: Heat Trades Q-Factor for Frequency Diversity

**Status:** Computationally verified (N=2 to N=5; Heisenberg chain). The
palindrome section is stronger: exact at N=2 and N=3, by rational characteristic
polynomial rather than by an eigensolver.
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Scripts:** [v_effect_gamma_sweep.py](../simulations/v_effect_gamma_sweep.py),
[v_effect_thermal.py](../simulations/v_effect_thermal.py),
[self_heating_fixpoint.py](../simulations/self_heating_fixpoint.py),
[thermal_palindrome_centre.py](../simulations/thermal_palindrome_centre.py)
**Results:** [v_effect_gamma_sweep.txt](../simulations/results/v_effect_gamma_sweep.txt),
[v_effect_thermal.txt](../simulations/results/v_effect_thermal.txt),
[self_heating_fixpoint.txt](../simulations/results/self_heating_fixpoint.txt),
[thermal_palindrome_centre.txt](../simulations/results/thermal_palindrome_centre.txt)
**Depends on:** [V-Effect](V_EFFECT_PALINDROME.md),
[Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md),
[Energy Partition](../hypotheses/ENERGY_PARTITION.md)

---

## What this document is about

A quantum resonator has two qualities you might want: high Q-factor (how
many times it oscillates before dying) and frequency diversity (how many
different frequencies it can sustain). This document shows that you cannot
maximize both at once. Cold systems have high Q but few frequencies. Hot
systems have many frequencies but low Q. The trade-off is mediated by
three independent mechanisms: coupling creates structured oscillation,
dephasing lifts degeneracies, and thermal excitation explodes the
frequency count but degrades quality. The sacrifice-zone advantage (3×
at zero temperature) vanishes in heat. Biology operates in the middle:
enough heat for diversity, enough structure for quality.

## Abstract

Three orthogonal mechanisms break symmetry in open quantum systems.
Each creates a different kind of complexity. They are not independent:
their interaction reveals a trade-off between resonator quality and
frequency diversity that has not been described before.

1. **Coupling** (V-Effect): creates palindromic pairs. Amplifies
   Q-factor by exactly 1.81x for a 5-qubit Heisenberg chain. Under
   uniform dephasing this ratio is a geometric constant and the rate γ
   cancels out of it; under a γ profile it does not survive.
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
N=5 system has Q=19 and 109 new frequencies. None of the original
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

The Q_max in the V(N) argument is the best Q on the chain's **(0,1)
coherence block**: the N-dimensional span of the |0⟩⟨j| between the
ferromagnet and the single excitations, which is the object F2's
dispersion describes
([D10](../docs/proofs/derivations/D10_W1_DISPERSION.md)). Under uniform
Z-dephasing every mode of that block decays at the same rate 2γ,
regardless of its oscillation frequency, because a |0⟩⟨j| coherence
disagrees at exactly one site, so the dissipator is the scalar −2γ·I
there (D10 Step 1). Therefore:

    Q_max = ω_max / (2γ)

The factor 2γ is identical for N=2 and N=5. It cancels in the ratio:

    V(N) = Q_max(N) / Q_max(2) = ω_max(N) / ω_max(2)

The gain is purely a frequency ratio, independent of the dephasing
RATE. It is not independent of the dephasing PROFILE; see the
non-uniform paragraph below.

The block is not the XY-weight-1 Pauli sector, which is far larger (160
dimensions against 5 at N=5), is not L-invariant, and therefore has no
spectrum of its own; D10 Step 6 carries that scope. Whether some other
joint-popcount block reaches a higher Q than this one at some N is a
separate question and is not settled here (Open Question 7 below).

This cancellation holds for uniform dephasing and for that alone. Under
a non-uniform profile (sacrifice zone) diag(γ) stops commuting with the
block's Laplacian, and both halves of the argument go at once: the decay
rates spread ([Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md))
AND the frequencies themselves move, measured level by level and growing
with the profile's unevenness ([Concentrator Optics](CONCENTRATOR_OPTICS.md)
Result 3). So the 1.81x ratio does not survive by retreating to the
extremal mode: ω_max moves too, and V(N) = 1 + cos(π/N) is a uniform-γ
statement.

What replaces it under real conditions is computed rather than derived, and
computed at assumed inputs rather than measured ones: the two numbers below are
Q_max = |Im λ|/|Re λ| at the stipulated substrate parameters, whose denominators
are unsourced ([Q Belongs to No Substance](../docs/Q_BELONGS_TO_NO_SUBSTANCE.md)).
At biological temperature (310 K, n̄ ~ 1.5-2):
G-C DNA base pair has Q_max = 0.57 (down from 1.95 cold), with
frequency diversity increasing from 15 to 26
([DNA Base Pairing](DNA_BASE_PAIRING.md)). N=3 water chain: Q_max
= 0.43 at 300 K ([Proton Water Chain](PROTON_WATER_CHAIN.md)).

### The exact formula

Computing the block's ω_max for the Heisenberg chain at each N reveals
an exact pattern (verified N=2 through N=6):

| N | ω_max | Exact form | V(N) = ω/ω(N=2) |
|:--|:-----------|:-----------|:-----------------|
| 2 | 4.0000 | 4J | 1.000 |
| 3 | 6.0000 | 4J + 4J·cos(π/3) = 6J | 1.500 |
| 4 | 6.8284 | 4J + 4J·cos(π/4) = (4+2√2)J | 1.707 |
| 5 | 7.2361 | 4J + 4J·cos(π/5) = (5+√5)J | **1.809** |
| 6 | 7.4641 | 4J + 4J·cos(π/6) = (4+2√3)J | 1.866 |

The formula:

    ω_max(N) = 4J · (1 + cos(π/N)) = 8J · cos²(π/(2N))

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
topological invariant. The correct characterization: it is the ratio of
the (0,1) coherence block's maximum Liouvillian frequencies, and those
are the single-magnon energies measured from the ferromagnetic vacuum
(states where exactly one spin is flipped, the simplest excitations of
the chain), because the block's generator is 2J times the chain's graph
Laplacian and the Laplacian is that one-magnon Hamiltonian.

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
decay broadens their linewidths (the frequency range over which each mode responds, like the width of a bell curve).

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

The reading below is the canonical F1 pairing distance about the spectrum's own
centre, in units of ε · spectral radius (N=5). At this N the eigensolver floor
is O(10-100), so a value in that range means "palindromic, to the accuracy the
solver can offer", and 1e14 means broken:

| Noise type | n̄=0 | n̄=0.5 | n̄=2 | n̄=10 |
|:-----------|:----|:------|:----|:-----|
| Pure Z-dephasing | 55 | n/a | n/a | n/a |
| Pure amplitude damping | 98 | 57 | 72 | 69 |
| Z-dephasing **and** amplitude | 4.1e14 | 1.1e15 | 1.1e15 | 2.7e15 |

**The column used to report a greedy first-fit percentage inside an absolute
tolerance of 1e-3, and that percentage measured its own matcher.** The spectrum
is far denser than the tolerance, so a rate rarely grabbed its true mirror, and
tightening the tolerance *saturates* the printed score long before it fixes the
pairing, which inverts what a high score means: measured on the sibling scan in
[Concentrator Qubit Mapping](CONCENTRATOR_MAPPING.md), at 1e-6 the scorer prints
100% while 56 of 480 accepted pairs are still wrong by more than 1e-12. The
amplitude-damping row therefore read as a symmetry weakening with temperature
when the pairing is exact at every n̄.

Worth naming, because it is the more interesting half: that row was scored
against the **right** centre all along. The script computed
`N·γ_amp·(1 + 2n̄)/2`, which is exactly Σ(γ↓+γ↑)/2, before anyone had written
that expression down as a centre. The centre was in the code and the claim was
in no document. Two of the other call sites did pass the Z-dephasing centre for
rows with an amplitude channel running; all of them read the centre off the
trace now, which removes the choice.

**Answered 2026-08-05, and the answer is a shifted centre.** The paragraph that
stood here called the thermal status open and cautioned that the centre was
only estimated, so the low scores might be centre error rather than broken
symmetry. Both halves are settled now, and the caution was the more important
of the two.

**The centre is not estimated; it is read off the trace.** If a multiset is
closed under λ ↦ 2c − λ then every pair sums to 2c, so Σλ = n·c and
c = mean(λ) exactly. There is one candidate centre, not a fitted one, which
turns "is it palindromic?" from a search into a check: a large distance at that
centre means **no** centre works.

**σ₊/σ₋ do preserve a modified palindromic structure with a shifted centre**,
which is the first of the two possibilities this section named:

| Channel | Palindrome centre | Status |
|:--------|:------------------|:-------|
| Z-dephasing | −Σγ | F1, proven |
| Amplitude damping, T = 0 | −Σγ/2 | [F137](../docs/ANALYTICAL_FORMULAS.md#f137) |
| Thermal bath, σ₋ and σ₊ | −Σ(γ↓ + γ↑)/2 | F137 extended, 2026-08-05 |
| Z-dephasing **and** amplitude damping | none exists | breaks |

Dephasing pays the full shift, amplitude damping half of it, and the thermal
bath half of the **total** per-site rate γ↓ + γ↑. Measured with the Heisenberg
H at N = 2 through 5 with independently drawn per-site rates: predicted centre
matching the trace to 1e-15, pairing at the eigensolver floor
([`thermal_palindrome_centre.py`](../simulations/thermal_palindrome_centre.py)).
The H = 0 case needs nothing new; it follows from the thermal per-site rates
[0, r/2, r/2, r], r = γ↓ + γ↑, which
[KMS_DETAILED_BALANCE](../docs/KMS_DETAILED_BALANCE.md) had already computed.

Better than measured, where the rates are rational: the characteristic
polynomial is then exact and the palindrome is the identity p(2c − x) ≡ p(x),
which holds at N=2 and N=3 with the Heisenberg H and no eigensolver anywhere.
The same test returns "neither" for Z + amplitude, so that break is proven, not
merely large.

**The "Z-deph + amplitude" row is the hardware case, and it genuinely breaks:**
no centre pairs it.
[MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) has the reason,
and it is sharper than "two channels are worse than one". *Transverse* dephasing
composes with T1 exactly, 64 of 64 configurations there; only **co-axial**
Z-dephasing breaks it, 8 of 64. Two channels break the mirror when they share an
axis, not because there are two of them.

**A trap, since a centre is now easy to compute.** T1 with co-axial Z and T1
with transverse X give the *same* centre to every digit, and only the first
breaks the pairing: the centre is the trace, and all three additions carry the
same total rate. The centre must be computed **and** the distance checked;
neither alone decides.

Where this leaves temperature: the palindrome's *existence* does not depend on
n̄, only on the total rate γ↓ + γ↑, and the centre moves with it since
γ↓ + γ↑ = γ(2n̄ + 1). There is no critical temperature to look for. And the
hardware operating point is not the hot one: a flown two-leg protocol
(`fw.Confirmations.lookup('f84_heating_leg_attribution_kingston_july2026')`)
found γ↑ consistent with zero and thermal populations of 0.23 to 0.83%. The
regime that matters on real machines is T1 beside co-axial Z, the one
combination that fails.

Still open: a proof of the H ≠ 0 pairing at general N. The *centre* at H ≠ 0 is
not open; it is a trace identity, and the commutator part of the Liouvillian is
traceless, so the centre never depended on H at all.

For contrast, [Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md) shows that
depolarizing noise (X+Y+Z) breaks pairing, with error < 0.1% at typical IBM γ
values.

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

**Tier 2** (loop concept from two Tier 2 results, fixed-point
computation performed March 31).

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

### Does the loop converge?

**Computed** ([self_heating_fixpoint.py](../simulations/self_heating_fixpoint.py),
[results](../simulations/results/self_heating_fixpoint.txt)).

Method: fixed-point iteration. At each n_bar, compute the Liouvillian
steady state ρ_ss and compare its energy E_ss = Tr(H·ρ_ss) with the
thermal energy E_th(n_bar). Adjust n_bar until E_ss = E_th.

Result for all 6 configurations tested (N=3 and N=5, pure amplitude
damping, Z+amplitude, sacrifice profile):

**The loop diverges.** n_bar runs away to infinity. The steady state
is always hotter than the bath at any finite n_bar.

| Config | E_ss (n_bar→0) | E_ground | Gap | Outcome |
|:-------|:---------------|:---------|:----|:--------|
| N=3, pure amp | +1.99 | -4.00 | 5.99 | diverges |
| N=3, Z+amp | +1.99 | -4.00 | 5.99 | diverges |
| N=5, pure amp | +3.98 | -7.71 | 11.70 | diverges |
| N=5, Z+amp | +3.98 | -7.71 | 11.70 | diverges |
| N=5, Z+amp (weak) | +3.98 | -7.71 | 11.70 | diverges |
| N=5, sacrifice+amp | +3.98 | -7.71 | 11.70 | diverges |

The reason: the Lindblad steady state under dephasing is the maximally
mixed state (E ≈ Tr(H)/d, near the spectral center). The thermal
ground state is at E_min. This gap never closes at finite n_bar.

### What divergence means

Without external cooling, the resonator thermalizes to maximum
entropy (n_bar → ∞, Q → 0, frequency diversity → maximum, structure
→ zero). The system is not self-sustaining.

To maintain structure, an external mechanism must hold n_bar at a
finite value. In quantum hardware: the cryostat (15 mK). In biology:
metabolism (ATP hydrolysis pumps heat out while coupling pumps
structure in).

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
At high temperature: diversity dominates (445 frequencies). The pairing is
not what heat costs: under the thermal bath alone it survives at every n_bar,
at a centre that moves with the total rate.
Life operates in between.

---

## Answered Questions

1. ~~Can 1.81x be derived analytically?~~ **ANSWERED (March 31).**
   V(N) = 1 + cos(π/N). For N=5: (5+√5)/4 ≈ 1.80902. The gain is
   the ratio of the (0,1) coherence block's maximum Liouvillian
   frequencies, which follow
   ω_max = 4J·(1+cos(π/N)). Verified N=2 through N=6. The golden
   ratio appears: cos(π/5) = φ/2. See derivation above.

2. ~~Does the self-heating loop converge?~~ **ANSWERED (March 31).**
   No. The loop diverges in all 6 configurations tested. The system
   thermalizes to maximum entropy without external cooling. See
   "The Self-Heating Loop" section above.

## Open Questions

### Computable (no hardware needed)

3. ~~What is the critical n_bar where the palindromic pairing drops
   below 50%?~~ Answered 2026-08-05: there is none. The pairing under a
   thermal bath is exact at every n_bar, and the smooth transition the old
   data suggested was the greedy scorer, not the physics. See the palindrome
   section above.

4. What external cooling rate stabilizes the system at a given
   n_bar? The fixed-point computation shows divergence without
   cooling. The inverse question: how much cooling for a target
   operating point? (Model: add a cold bath channel that competes
   with self-heating.)

### Requires controlled thermal injection on hardware

5. Can the frequency-diversity explosion at n_bar > 0 be observed
   on superconducting qubit hardware by intentionally heating
   qubits (e.g., driving with a thermal microwave field)?
   Standard cryogenic operation holds n_bar ≈ 0; this would
   require deliberate thermal injection outside normal operating
   conditions.

6. Does the sacrifice-zone advantage recover at intermediate
   temperatures if the sacrifice qubit is selectively heated?
   Same requirement: controlled per-qubit thermal injection.

### Analytical (proof needed)

7. Does any other joint-popcount block of the Liouvillian reach a
   higher Q than the (0,1) coherence block? The V(N) argument reads
   Q_max off that block alone, and the step from "best on the block"
   to "best in the Liouvillian" has never been taken. It is the one
   the V-Effect reading needs. (The thermal and sacrifice-zone Q_max
   values quoted above come from other documents' runs under
   amplitude damping or a gamma profile, where neither the block's
   closed form nor this question applies as stated.)

   What used to stand here, a formal proof that
   ω_max = 4J·(1+cos(π/N)) for all N, was settled in April and is no
   longer open: ω_max is the k = N−1 member of F2's dispersion
   ω_k = 4J·(1 − cos(πk/N)), derived from the block's tight-binding
   reduction in
   [D10](../docs/proofs/derivations/D10_W1_DISPERSION.md).

---

## Tier Assessment

- 1.81x geometric constant: **Tier 1-2** (exact formula
  V(N) = 1+cos(π/N), verified N=2-6 to machine precision.
  Analytical derivation from the (0,1) coherence block's Liouvillian
  eigenfrequencies; ω_max = 4J(1+cos(π/N)) is proven for all N in D10.
  What keeps it off a clean Tier 1 is the unproven step from that
  block's best Q to the Liouvillian's, Open Question 7)
- Frequency diversity γ dependence: **Tier 2**
- Thermal breaking of 1.81x: **Tier 2** (11 n_bar values, three
  noise configurations)
- Sacrifice zone temperature dependence: **Tier 2**
- Self-heating divergence: **Tier 2** (6 configs, N=3 and N=5,
  fixed-point iteration, all diverge to n_bar → ∞)
- Biological interpretation: **Tier 4** (motivated by computation,
  parameters from literature, not tested)
- Hierarchy mapping: **Tier 5** (speculative)
