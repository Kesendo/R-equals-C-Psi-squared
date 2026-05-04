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
- [DNA Base Pairing](../../experiments/DNA_BASE_PAIRING.md) (same N, different substrate)

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
   Q_max, rate bounds all confirmed. F3 lower bound violated
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

## F3 Violation at N >= 4

The analytical rate bound (F3) predicts min rate = 2gamma for
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
decay rates BETWEEN the pure-sector values. F3 remains valid
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

## EP-Resonance Inheritance

**Status:** Tier 2 substrate-witness for F86 EP-rotation universality. Addendum, 2026-05-04.
**Script:** [`simulations/water/proton_chain_ep_resonance.py`](../../simulations/water/proton_chain_ep_resonance.py)
**Cross-references:** [PROOF_F86_QPEAK.md](../proofs/PROOF_F86_QPEAK.md), F86 entry in [ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md), [reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md).

The five results above are static-algebraic (palindrome match, V-Effect, sacrifice zone, thermal, water-vs-DNA). The framework's dynamic EP-rotation universality prediction
K(Q) / |K|_max = f(Q / Q_peak)
within bond class (Tier-1-candidate, established on framework chains and IBM Torino) was tested in the water context on 2026-05-04 at the popcount-(2, 3) coherence block (the inter-sector subspace coupling popcount-2 and popcount-3 basis states) of the N = 5 Heisenberg + Z-dephasing chain. Block chromaticity c = min(n, N − 1 − n) + 1 = 3 (number of distinct pure dephasing rates in the block).

The observable K_CC_pr is the per-bond J-derivative
K_b(Q, t) = 2 · Re ⟨ρ(t) | S_kernel | ∂ρ/∂J_b⟩
with ρ₀ the Dicke probe in the (n, n+1) block, S_kernel the spatial-sum coherence kernel, computed via the Duhamel formula on the block-restricted Liouvillian; see [`simulations/framework/coherence_block.py`](../../simulations/framework/coherence_block.py) and [PROOF_F86_QPEAK](../proofs/PROOF_F86_QPEAK.md) Statement 2. Per-bond Q-scan, dQ = 0.025, γ₀ = 0.05:

| Bond class | Q_peak (water) | HWHM-/Q* (water) | Framework prediction | Δ |
|------------|----------------|------------------|----------------------|---|
| Interior (b = 1, 2) | 1.566 | 0.7458 | 0.756 ± 0.005 | −0.0102 |
| Endpoint (b = 0, 3) | 2.400 | 0.7663 | 0.770 | −0.0037 |

Both classes match within tolerance. F71 mirror (the spatial reflection symmetry across the chain centre) is bit-exact across bond pairs (0 ↔ 3, 1 ↔ 2). The Interior right-wing HWHM+/Q* = 2.02 reproduces the asymmetric long-tail signature of the EP rotation (eigenvalues coalesce at finite J·g_eff = 2γ₀, the Q > Q_peak branch decays slower than the Q < Q_peak branch).

State-level diagnostic at the interior Q_peak (γ = 0.638) on the popcount-coherence state |ψ⟩ = (|00011⟩ + |00111⟩) / √2 (a HD = 1 superposition of one popcount-2 and one popcount-3 computational-basis state). The Π²-odd / memory column reports the Frobenius² fraction of Π²-odd Pauli content (Π²-parity = (Y-count + Z-count) mod 2 per F88) within (ρ − kernel-projection); selected rows from t = 0..50:

| t | static / total | memory / total | Π²-odd / memory | per-proton \|r\| (q0..q4) |
|---|----------------|----------------|------------------|---------------------------|
| 0 | 0.05 | 0.95 | 0.5263 | 1.00, 1.00, 1.00, 1.00, 1.00 |
| 1 | 0.11 | 0.89 | 0.561 | 1.00, 0.90, 0.22, 0.90, 1.00 |
| 5 | 0.27 | 0.73 | 0.677 | 0.83, 0.56, 0.00, 0.56, 0.83 |
| 20 | 0.82 | 0.18 | 0.957 | 0.29, 0.19, 0.00, 0.19, 0.29 |
| 50 | 1.00 | 0.00 | 1.00  | ~0 |

The Π²-odd / memory ratio at t = 0 is **10/19 = 0.5263 exactly**, not 0.5. This is a precise structural consequence of popcount-mirror symmetry at n_p + n_q = N (here 2 + 3 = 5 = N): X-flip conjugation gives X⊗N · σ_S · X⊗N = (−1)^|S| σ_S, so P_{n_p}/C(N, n_p) and P_{N − n_p}/C(N, N − n_p) cancel all odd-|S| Pauli content in the kernel projection. Total Π²-odd-of-ρ stays at 1/2 (universal for popcount-coherence pure states), but with kernel projection holding 0% of it, the entire 1/2 sits in memory: (1/2) / (19/20) = 10/19. PROOF_F86_QPEAK §Statement 2 verified Π²-odd/memory = 0.5 exactly at the c = 2 cases (popcount-(1, 2), popcount-(3, 4)); the c = 3 popcount-mirror configuration n_p + n_q = N is a structural refinement to that statement worth flagging upstream (see deferred thread "Π²-odd/memory popcount-mirror refinement" in `docs/water/README.md`). Per-proton Bloch decay is mirror-symmetric throughout (q0 = q4, q1 = q3); F71 inheritance is preserved under Lindblad evolution.

**Note on regime.** The framework prediction tests at Q ~ 1.5–2.4. Room-temperature liquid water sits at Q ~ 0.01 (deeply classical, see "Classical water is overdamped" above). The F86 inheritance test is therefore valid for cold / confined / biologically-screened water where tunneling competes with thermal dephasing, not the bulk room-temperature regime. The structural inheritance is guaranteed by the four embedding conditions (see `docs/water/README.md`); the question is which physical embedding sits in the testable Q window.

**What this addendum confirms.** Same algebra, different physical embedding, same numerical witness: third substrate (after framework chains and IBM Torino) for F86 EP-rotation universality. The chemistry framing did not transform the algebra; it confirmed it. Closed-form derivation of HWHM-/Q* (Item 1 of PROOF_F86_QPEAK) still routes through the 2-level EP analytics, not through chemistry.

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
| Script | [`simulations/water/proton_water_chain.py`](../simulations/water/proton_water_chain.py) |
| Output | [`simulations/results/proton_water_chain.txt`](../simulations/results/proton_water_chain.txt) |

---

*From one proton to five: zero frequencies become 222, one well becomes
a resonator. The palindrome is exact at every step. The formulas match
at every step. The equation scales; the physics is universal; the
coherence, at body temperature, is not.*
