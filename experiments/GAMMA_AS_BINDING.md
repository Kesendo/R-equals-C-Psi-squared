# Gamma as Binding Parameter: Per-Sector Rate Sensitivity

**Status:** Complete. Closes EQ-008.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/gamma_sector_sensitivity.py`
**Data:** `simulations/results/gamma_binding/`
**Related:** [Symmetry Census](SYMMETRY_CENSUS.md) (sector enumeration), [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (rate formula)

---

## Motivation (EQ-008)

The V2 boundary straddling analysis showed that the Bell-pair coherence (cross-sector block) and the slow-mode coupling (SE block) are orthogonal: they do not see each other. But both are shaped by the same physical parameter: the dephasing profile gamma. EQ-008 asks: how does gamma act on each sector, and is this "binding" a simple linear scale or something richer?

---

## Frage 1: Per-sector slowest rates at N=5

Three gamma profiles with matched sum Sigma_gamma = 2.608 (so differences come from spatial distribution, not total noise):

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

The palindromic symmetry rate(w) = rate(N-w) holds for uniform (by spin-flip symmetry) and sacrifice profile. It holds for moderate gradient only approximately (broken by the asymmetric profile).

### Off-diagonal sectors (w != w', coherence sectors)

The slowest rate in each coherence sector measures how fast cross-sector coherences decay. Key pattern: rates increase with |w - w'| (more distant sectors decohere faster). The (0,5) and (5,0) sectors always have rate = Sigma_gamma (the maximum), because they connect maximally distant sectors and require N bit flips.

The sacrifice profile compresses all coherence rates toward lower values compared to uniform (at matched Sigma_gamma). This is a global consequence of the sacrifice geometry: concentrating noise on one site slows down the entire Liouvillian.

---

## Frage 2: Scaling with alpha (gamma_k -> alpha * gamma_k)

If rates scaled linearly with alpha, the ratio rate(alpha) / (alpha * rate(1)) would be 1.0 everywhere.

### Linearity deviations

| Profile | Max deviation from linearity |
|---------|------------------------------|
| Uniform | 0.94 (94%) |
| Sacrifice | 1.34 (134%) |
| Moderate | 0.94 (94%) |

**Scaling is massively nonlinear.** The deviations exceed 100%, meaning rates can be off by more than a factor of 2 from the linear prediction.

### Why: the absorption theorem is per-mode, not per-sector

The absorption theorem gives Re(lambda) = -2 sum gamma_k <1_XY(k)> for each eigenmode. If gamma scales by alpha, each eigenvalue scales by alpha. This is exact, per mode.

But the SLOWEST MODE PER SECTOR can change identity as alpha varies. The Liouvillian is L(alpha) = L_H + alpha * L_D. At small alpha, the Hamiltonian dominates; at large alpha, the dephasing dominates. The eigenvector structure changes with alpha, so different modes become the "slowest" in a sector at different alpha values.

### Per-sector scaling pattern (sacrifice profile)

| Sector | rate(0.5)/expected | rate(4.0)/expected | Behavior |
|--------|--------------------|--------------------|----------|
| (1,1) SE | 1.72 | 0.37 | Hamiltonian-protected at high gamma |
| (2,2) interior | 1.60 | 0.26 | Same pattern, stronger |
| (0,1) edge | 1.67 | 0.74 | Moderate deviation |
| (0,5) extremal | 1.00 | 1.00 | Exactly linear (trivial sector) |

At weak dephasing (alpha = 0.5): rates are HIGHER than linear prediction (the Hamiltonian redistributes energy and makes modes decay faster than pure dephasing would predict).

At strong dephasing (alpha = 4.0): rates are LOWER than linear prediction (a quantum-Zeno-like effect where strong measurement slows down the dynamics).

The only perfectly linear sector is (0,5): this sector has dimension 1 and its rate is simply Sigma_gamma, which scales trivially.

---

## Conclusion: gamma as binding parameter (EQ-008 answer)

Gamma connects the sectors not by coupling them (sector conservation is exact), but by **differentially modulating** their dynamics:

1. The spatial profile of gamma (which sites get more noise) determines the rate hierarchy across sectors. This is the sacrifice effect.

2. The global scale of gamma determines the balance between Hamiltonian (coherent) and dissipative (dephasing) dynamics within each sector. This is nonlinear: doubling gamma does NOT double all sector rates.

3. Different sectors respond with different sensitivity to gamma scaling. The interior sectors (2,2) and (3,3) show the strongest nonlinearity (Zeno effect). The edge sectors (0,1) show moderate nonlinearity. The extremal sector (0,5) is trivially linear.

Gamma is therefore a **nonlinear common modulator**: it shapes all sectors simultaneously but with sector-dependent sensitivity. It does not couple sectors (no population transfer), but the fact that different sectors respond differently to the same gamma change means there is a structural relationship between sectors mediated by gamma.

---

## Files

- `simulations/gamma_sector_sensitivity.py` (computation)
- `simulations/results/gamma_binding/gamma_binding_results.json` (raw data)
- `simulations/results/gamma_binding/sector_rate_heatmaps.png` (heatmaps)
- `simulations/results/gamma_binding/scaling_loglog.png` (scaling plot)
- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (per-mode linear scaling)
- [Symmetry Census](SYMMETRY_CENSUS.md) (sector enumeration)

---

*April 12, 2026. EQ-008 closed: gamma is a nonlinear common modulator, not a linear scale parameter.*
