# Gamma as Binding Parameter: Per-Sector Rate Sensitivity

> **Superseded by [LIGHT_DOSE_RESPONSE](LIGHT_DOSE_RESPONSE.md).** This document (V1) correctly measured per-sector rates and the 134% nonlinearity, but its mechanistic explanation (mode-crossing hypothesis, Zeno language) was disproved by V2. The numbers in the tables below are correct; the "Why" section is not. See V2 for the correct mechanism (eigenvector rotation).

**Status:** Superseded. EQ-008 closed via V2.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/gamma_sector_sensitivity.py`
**Data:** `simulations/results/gamma_binding/`
**Related:** [Symmetry Census](SYMMETRY_CENSUS.md) (sector enumeration), [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (rate formula)

---

## Motivation (EQ-008)

The V2 boundary straddling analysis showed that the Bell-pair coherence (cross-sector block) and the slow-mode coupling (SE block) are orthogonal: they do not see each other. But both are shaped by the same physical parameter: the dephasing profile γ. EQ-008 asks: how does γ act on each sector, and is this "binding" a simple linear scale or something richer?

---

## Frage 1: Per-sector slowest rates at N=5

Three γ profiles with matched sum Σγ = 2.608 (so differences come from spatial distribution, not total noise):

### Diagonal sectors (w = w', contain the N+1 steady-state exits)

| w | Uniform | Sacrifice (IBM) | Moderate gradient |
|---|---------|-----------------|-------------------|
| 0 | 0.000 | 0.000 | 0.000 |
| 1 | 1.502 | 0.318 | 1.564 |
| 2 | 1.154 | 0.375 | 1.131 |
| 3 | 1.154 | 0.375 | 1.131 |
| 4 | 1.502 | 0.318 | 1.564 |
| 5 | 0.000 | 0.000 | 0.000 |

The sacrifice profile has dramatically slower rates in the SE sector (0.318 vs 1.502 for uniform). This is the lens effect: the sacrifice geometry creates one slow mode in the SE sector that is absent under uniform dephasing. All three profiles have w=0 and w=5 as trivially stationary (one-dimensional sectors, nothing to decay).

The palindromic symmetry rate(w) = rate(N−w) holds exactly for all three profiles. This follows from spin-flip symmetry (X⊗N conjugation commutes with L for any γ profile, mapping sector (w,w) to (N−w,N−w) with identical spectra).

### Off-diagonal sectors (w ≠ w', coherence sectors)

The slowest rate in each coherence sector measures how fast cross-sector coherences decay. Key pattern: rates increase with |w − w'| (more distant sectors decohere faster). The (0,5) and (5,0) sectors always have rate = 2Σγ (the maximum), because they connect maximally distant sectors (the single basis element |00...0⟩⟨11...1| acquires a phase factor from every Z_k).

The sacrifice profile compresses all coherence rates toward lower values compared to uniform (at matched Σγ). This is a global consequence of the sacrifice geometry: concentrating noise on one site slows down the entire Liouvillian.

---

## Frage 2: Scaling with α (γ_k → α·γ_k)

If rates scaled linearly with α, the ratio rate(α) / (α × rate(1)) would be 1.0 everywhere.

### Linearity deviations

| Profile | Max deviation from linearity |
|---------|------------------------------|
| Uniform | 0.94 (94%) |
| Sacrifice | 1.34 (134%) |
| Moderate | 0.94 (94%) |

**Scaling is massively nonlinear.** The deviations exceed 100%, meaning rates can be off by more than a factor of 2 from the linear prediction.

### Why: eigenvector rotation (corrected in V2)

> **V1 hypothesis (below) was wrong.** V1 proposed mode-crossing as the mechanism. [LIGHT_DOSE_RESPONSE](LIGHT_DOSE_RESPONSE.md) (V2) tracked individual eigenvalue curves and found zero level crossings. The actual mechanism is eigenvector rotation: as α changes, each mode's Pauli content shifts continuously, making Re(λ) = −2α Σ γ_k ⟨1_XY(k)⟩ nonlinear because ⟨1_XY(k)⟩ itself depends on α.

~~The absorption theorem gives Re(λ) = −2 Σ γ_k ⟨1_XY(k)⟩ for each eigenmode. If γ scales by α, each eigenvalue scales by α. This is exact, per mode.~~

~~But the SLOWEST MODE PER SECTOR can change identity as α varies.~~ The Liouvillian is L(α) = L_H + α·L_D. At small α, the Hamiltonian dominates; at large α, the dephasing dominates. The eigenvector structure changes with α, rotating the Pauli content of each mode and making individual rates nonlinear.

### Per-sector scaling pattern (sacrifice profile)

| Sector | rate(0.5)/expected | rate(4.0)/expected | Behavior |
|--------|--------------------|--------------------|----------|
| (1,1) SE | 1.72 | 0.37 | Hamiltonian-protected at high gamma |
| (2,2) interior | 1.60 | 0.26 | Same pattern, stronger |
| (0,1) edge | 1.67 | 0.74 | Moderate deviation |
| (0,5) extremal | 1.00 | 1.00 | Exactly linear (trivial sector) |

At weak dephasing (α = 0.5): rates are HIGHER than linear prediction (the Hamiltonian shapes the eigenvectors toward higher XY-weight content).

At strong dephasing (α = 4.0): rates are LOWER than linear prediction (the eigenvectors rotate toward the Pauli basis, reducing their effective XY-weight).

The only perfectly linear sector is (0,5): this sector has dimension 1 and its rate is 2Σγ, which scales trivially.

---

## Conclusion: gamma as binding parameter (EQ-008 answer)

Gamma connects the sectors not by coupling them (sector conservation is exact), but by **differentially modulating** their dynamics:

1. The spatial profile of γ (which sites get more noise) determines the rate hierarchy across sectors. This is the sacrifice effect.

2. The global scale of γ determines the balance between Hamiltonian (coherent) and dissipative (dephasing) dynamics within each sector. This is nonlinear: doubling γ does NOT double all sector rates.

3. Different sectors respond with different sensitivity to γ scaling. The interior sector (2,2) shows the strongest nonlinearity (eigenvector rotation, see [V2](LIGHT_DOSE_RESPONSE.md)). The edge sectors (0,1) show moderate nonlinearity. The extremal sector (0,5) is trivially linear.

γ is therefore a **nonlinear common modulator**: it shapes all sectors simultaneously but with sector-dependent sensitivity. It does not couple sectors (no population transfer), but the fact that different sectors respond differently to the same γ change means there is a structural relationship between sectors mediated by γ.

---

## Files

- `simulations/gamma_sector_sensitivity.py` (computation)
- `simulations/results/gamma_binding/gamma_binding_results.json` (raw data)
- `simulations/results/gamma_binding/sector_rate_heatmaps.png` (heatmaps)
- `simulations/results/gamma_binding/scaling_loglog.png` (scaling plot)
- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (per-mode linear scaling)
- [Symmetry Census](SYMMETRY_CENSUS.md) (sector enumeration)

---

*April 12, 2026. EQ-008 closed: γ is a nonlinear common modulator, not a linear scale parameter. Mechanism identified in [V2](LIGHT_DOSE_RESPONSE.md): eigenvector rotation.*
