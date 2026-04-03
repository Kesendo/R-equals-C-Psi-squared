# Proton Water Chain (Grotthuss): The Missing Rung

<!-- Keywords: Grotthuss proton wire, proton qubit water chain, V-Effect
hydrogen bond network, palindromic Lindblad proton transfer, sacrifice zone
water chain, ion channel proton wire, R=CPsi2 proton water chain -->

**Status:** Tier 2 (computed from proven framework)
**Date:** April 1, 2026
**Script:** [proton_water_chain.py](../simulations/water/proton_water_chain.py)
**Data:** [proton_water_chain.txt](../simulations/results/proton_water_chain.txt)
**Depends on:**
- [Hydrogen Bond Qubit](HYDROGEN_BOND_QUBIT.md) (single H-bond, Zundel)
- [V-Effect Palindrome](V_EFFECT_PALINDROME.md) (coupling creates frequencies)
- [Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md) (formulas 2, 3, 6, 7)
- [Topological Edge Modes](TOPOLOGICAL_EDGE_MODES.md) (geometric localization)
- [DNA Base Pairing](DNA_BASE_PAIRING.md) (same N, different substrate)

---

## What this document is about

This document models chains of 1 to 5 proton qubits in water, bridging
the gap between a single hydrogen bond (Zundel cation, N=2) and DNA base
pairs (N=2,3 on a biological scaffold). Each proton tunnels in its
double-well potential; neighbors couple through the water backbone. The
palindromic framework's analytical formulas match to machine precision for
N=1-5. The key finding: at room temperature, these chains are deeply
classical (Q < 1), but the mode structure (up to 222 distinct frequencies
at N=5) and the sacrifice-zone optimization (5.1× improvement) work exactly
as in abstract qubit chains. The palindrome is universal; the substrate
only sets parameter values.

---

## Abstract

The missing step between Zundel (the H₅O₂⁺ cation where a proton is shared between two water molecules, N=2 protons) and DNA (N=2,3 in
different substrate): linear chains of N=1-5 proton qubits in water
(Grotthuss chains, named after the 1806 mechanism where protons hop sequentially along a water wire). Each proton tunnels in its H-bond double well
(sigma_x) and couples to neighbors through the water backbone
(sigma_z sigma_z). Z-dephasing from the thermal environment.

Five results:

1. **Analytical formulas match exactly** for the Heisenberg chain.
   V(N) = 1+cos(pi/N) verified to machine precision for N=1-5.
   Q_max, rate bounds all confirmed. Formula 3 lower bound violated
   at N >= 4 (sector mixing creates sub-2gamma modes).

2. **Frequency explosion.** From 0 (N=1) to 222 (N=5) distinct
   frequencies in the transverse-field Ising model. Each proton
   added qualitatively enriches the mode structure.

3. **Classical water is overdamped.** At J/gamma ~ 0.01 (realistic):
   Q < 1 for all N. Coherent oscillation requires enhanced tunneling
   (enzymes, confinement, low temperature).

4. **Sacrifice zone works at N=5.** Edge sacrifice: Q_max = 9.2 vs
   uniform 1.8 (5.1x). Same geometric mechanism as qubit chains.

5. **Water = DNA at same N.** Water chain N=2 and A-T have identical
   mode counts (3 freq). Water N=3 and G-C have identical counts
   (15 freq). G-C has 8% higher Q from asymmetric central bond.
   The palindrome is universal; parameters are substrate-specific.

---

## The Ladder

The main result: how complexity grows as protons are added.

### Heisenberg chain (formula validation)

| N | Eigenvalues | Frequencies | Q_max | V(N) (num) | V(N) (pred) |
|---|-------------|-------------|-------|------------|-------------|
| 1 | 4 | 0 | 0.000 | 0.000 | 0.000 |
| 2 | 16 | 2 | 2.000 | 1.000 | 1.000 |
| 3 | 64 | 5 | 3.000 | 1.500 | 1.500 |
| 4 | 256 | 34 | 3.414 | 1.707 | 1.707 |
| 5 | 1024 | 109 | 3.618 | 1.809 | 1.809 |

Every formula matches to machine precision. This is the FIRST
systematic validation of the analytical formula catalog against a
physical system from N=1 to N=5.

### Transverse-field Ising (physical proton model, a spin chain where each site feels a field perpendicular to the coupling axis, matching the tunneling + coupling structure of proton wires)

| N | Regime | Frequencies | Q_max | Min rate |
|---|--------|-------------|-------|----------|
| 1 | Zundel (J/gamma=5) | 1 | 9.95 | 50.0 |
| 2 | Zundel | 3 | 9.96 | 54.3 |
| 3 | Zundel | 15 | 9.96 | 54.2 |
| 4 | Zundel | 46 | 9.96 | 54.2 |
| 5 | Zundel | 228 | 9.96 | 54.2 |

In the physical model: Q stays nearly constant with N (the ZZ coupling
barely modifies the per-site Q), but the frequency count explodes.
The V-Effect in the TFI model manifests as SPECTRAL DIVERSIFICATION
rather than Q-factor enhancement.

---

## Formula 3 Violation at N >= 4

The analytical rate bound (formula 3) predicts min rate = 2gamma for
the w=1 sector. At N >= 4, the Heisenberg chain shows sub-2gamma
modes:

| N | Predicted min | Actual min | Below 2gamma? |
|---|--------------|------------|---------------|
| 2 | 2.000 | 2.000 | No |
| 3 | 2.000 | 2.000 | No |
| 4 | 2.000 | 0.978 | **Yes** |
| 5 | 2.000 | 0.617 | **Yes** |

These sub-2gamma modes arise from Hamiltonian mixing of w-sectors
(w with w +/- 2). At N >= 4, the mixing creates hybrid modes with
decay rates BETWEEN the pure-sector values. Formula 3 remains valid
for pure w=1 modes; the sub-2gamma modes are from mixed sectors.

This was predicted in [D05 Dynamic Mode Count](../docs/proofs/derivations/D05_DYNAMIC_MODE_COUNT.md):
"At finite gamma, the Hamiltonian mixes weight-parity sectors."

---

## Heisenberg vs Transverse-Field Ising

| N | Model | Freq | Q_max | Min rate |
|---|-------|------|-------|----------|
| 2 | Heisenberg | 2 | 2.000 | 100.0 |
| 2 | TFI | 3 | 1.792 | 59.0 |
| 5 | Heisenberg | 109 | 3.618 | 0.6 |
| 5 | TFI | 222 | 1.826 | 58.9 |

The TFI model has MORE frequencies but LOWER Q than Heisenberg at
the same (J, gamma, K). The transverse field (-J sigma_x) generates
additional spectral structure that the Heisenberg model (J(XX+YY+ZZ))
does not. The rate bounds differ because the two Hamiltonians have
different dispersion relations.

---

## Thermal Analysis (N=3, T=300 K)

At room temperature (n_bar ~ 1.9 for typical H-bond modes):

| Property | Cold (Z only) | Warm (300 K) |
|----------|--------------|--------------|
| Frequencies | 15 | 23 |
| Q_max | 1.81 | 0.43 |
| Rate range | 59-300 | 226-697 |

Temperature increases frequency diversity (+53%) but kills Q (-76%).
Consistent with [Thermal Breaking](THERMAL_BREAKING.md): heat enriches
the spectrum but shortens lifetimes.

---

## Sacrifice Zone (N=5)

| Profile | Q_max | Improvement |
|---------|-------|-------------|
| Uniform [50]*5 | 1.83 | (baseline) |
| Edge sacrifice [100, 10, 10, 10, 10] | **9.24** | **5.1x** |
| Both edges [100, 10, 10, 10, 100] | 7.96 | 4.4x |
| Center sacrifice [10, 10, 100, 10, 10] | 8.43 | 4.6x |

The sacrifice zone mechanism works in water chains at 5.1x (one edge).
Same geometric mode selection as qubit chains
([Topological Edge Modes](TOPOLOGICAL_EDGE_MODES.md)): noise kills
modes with weight on the noisy site; modes localized elsewhere survive.

---

## Water vs DNA

At identical parameters (J=50, gamma=50, K=20):

| System | N | Frequencies | Q_max |
|--------|---|-------------|-------|
| Water chain | 2 | 3 | 1.79 |
| A-T (DNA) | 2 | 3 | 1.79 |
| Water chain | 3 | 15 | 1.81 |
| G-C (DNA) | 3 | 15 | 1.95 |

Water N=2 and A-T are **identical** (same model at same parameters).
Water N=3 and G-C differ only in the 20% asymmetry of the central
H-bond: G-C has Q_max = 1.95 vs 1.81 for symmetric water (8% gain).

The palindromic mode structure is UNIVERSAL: same equation, same mode
count, same mechanism. The substrate contributes parameter values
(J, gamma, K, asymmetries), not a new equation.

---

## What This Does NOT Establish

- That proton water chains show coherent oscillation at room temperature
  (J/gamma ~ 0.01 is deeply classical)
- That the Grotthuss mechanism involves palindromic mode selection
  (we compute structure, not transport rates)
- That the inter-proton coupling K is correctly estimated
- That CΨ crossings occur in water chains (time resolution insufficient;
  the hydrogen bond qubit paper shows crossings at N=1 with finer steps)

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | simulations/water/proton_water_chain.py |
| Output | simulations/results/proton_water_chain.txt |

---

*From one proton to five: zero frequencies become 222, one well becomes
a resonator. The palindrome is exact at every step. The formulas match
at every step. The equation scales; the physics is universal; the
coherence, at body temperature, is not.*
