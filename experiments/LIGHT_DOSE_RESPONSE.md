# Light Dose Response Per Sector

**Status:** Complete. Mechanism identified: eigenvector rotation, not mode-crossing.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/light_dose_response.py`
**Data:** `simulations/results/light_dose_response/`
**Supersedes:** [GAMMA_AS_BINDING](GAMMA_AS_BINDING.md) (V1, correct numbers, wrong framing)
**Depends on:**
- [EXCLUSIONS](../docs/EXCLUSIONS.md) (Exclusion 2: gamma is external, five internal candidates eliminated)
- [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md) (chain is a passive optical cavity, 4/5 tests)
- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (per-mode rate formula)
- [Symmetry Census](SYMMETRY_CENSUS.md) (sector enumeration)

**Interpreted by:** [On the Light and What Casts Shadows in It](../reflections/ON_THE_LIGHT_AND_WHAT_CASTS_SHADOWS_IN_IT.md) (names the mechanism: shadows cast by self-structure, by sector profile, and by eigenvector rotation)

---

## Context

γ is external ([EXCLUSIONS](../docs/EXCLUSIONS.md) Exclusion 2). On IBM hardware, γ is literally light: microwave photons in a physical resonator cause dephasing through photon shot noise (Sears et al., Phys. Rev. B 86, 180504, 2012). The qubit chain is a passive optical cavity (4/5 cavity tests, R² = 0.998 for beam profile).

V1 ([GAMMA_AS_BINDING](GAMMA_AS_BINDING.md)) measured per-sector rates and found 134% deviation from linear scaling. V1's numbers were correct but its framing (γ as neutral "binding parameter") ignored the established status of γ as light. V2 reframes: γ is light, the table is a light-dose response curve, and the question is what mechanism produces the nonlinearity.

---

## Frage 1: Light dose per sector (sacrifice profile, alpha = 1.0)

The sacrifice profile concentrates light on site 0 (γ₀ = 2.34, quiet sites γ_k ≈ 0.05-0.10). The per-sector response:

| Sector type | Example | Slowest rate | Interpretation |
|-------------|---------|-------------|----------------|
| SE diagonal (1,1) | 0.318 | Lens mode: slowest in the system |
| Interior diagonal (2,2) | 0.375 | Slightly faster than SE |
| Adjacent coherence (0,1) | 0.167 | Cross-sector coherence, slowest off-diagonal |
| Distant coherence (0,5) | 5.216 = 2Σγ | Maximum rate: full dephasing on all N sites |
| Extremal diagonal (0,0) | 0.000 | Stationary (vacuum, nothing to decay) |

The [sacrifice geometry](SACRIFICE_GEOMETRY.md) protects the SE sector (rate 0.318) by concentrating light away from the cavity's slow standing-wave mode. This is the lens effect quantified per sector.

**Read as shadows** ([reflection](../reflections/ON_THE_LIGHT_AND_WHAT_CASTS_SHADOWS_IN_IT.md)): the rate in the Interpretation column is the depth of the shadow that sector's slowest mode casts on itself. The extremal diagonal (0,0) is pure self-shadow; no XY content, no light reaches it, rate zero. The distant coherence (0,5) is fully illuminated; every site contributes XY weight, rate 2Σγ. The SE and interior sectors sit in between, with partial shadows shaped by their mode space and the incident γ profile. The sacrifice geometry works by aligning the light with the sites where the SE slow mode has already placed its own shadow.

---

## Frage 2: Linear vs nonlinear light response

### V1 finding: 134% deviation from linearity

V1 scaled all γ_k by α and found per-sector rates deviate up to 134% from the linear prediction rate(α) = α × rate(1).

### V2 finding: mechanism is eigenvector rotation, not mode-crossing

V2 tracked individual eigenvalue curves across α ∈ [0.1, 5.0] with 50 steps. Two hypotheses tested:

**Hypothesis A (mode-crossing):** Individual modes are linear, but different modes become the "slowest" at different α, making the minimum nonlinear.

**Hypothesis B (eigenvector rotation):** Individual modes are themselves nonlinear because the Hamiltonian-dissipator balance changes their Pauli content.

### Evidence

| Sector | Per-mode R² (mean) | Per-mode R² (min) | Sector-min R² | Level crossings |
|--------|---------------------|---------------------|----------------|-----------------|
| (1,1) SE | 0.889 | 0.063 | 0.404 | 0 |
| (2,2) interior | 0.801 | 0.001 | 0.144 | 0 |
| (0,1) edge | 0.935 | 0.807 | 0.996 | 0 |
| (1,2) cross | 0.876 | 0.001 | 0.967 | 0 |
| (2,3) cross | 0.798 | 0.006 | 0.785 | 0 |

**Zero level crossings in all sectors.** The slowest mode maintains its identity across the full α range. Hypothesis A is ruled out.

**Individual modes are themselves nonlinear** (R² as low as 0.001). Hypothesis B is confirmed.

### Mechanism: eigenvector rotation

The Liouvillian is L(α) = L_H + α·L_D. As α increases:

1. At small α (Hamiltonian-dominated): eigenvectors resemble the Hamiltonian's normal modes. These have mixed Pauli content, so the [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) rate Re(λ) = −2α Σ γ_k ⟨1_XY(k)⟩ is sensitive to the eigenvector's XY-weight.

2. At large α (dephasing-dominated): eigenvectors rotate toward the Pauli basis (each mode approaches a definite Pauli string). The XY-weight ⟨1_XY(k)⟩ stabilizes.

3. The transition between these regimes is continuous (no crossings), but the eigenvector rotation changes ⟨1_XY(k)⟩ nonlinearly in α. This makes Re(λ) nonlinear even for a fixed mode.

The interior sector (2,2) shows the strongest nonlinearity (R² = 0.14) because it has the largest mode space (dim = 100) and the most room for eigenvector rotation. The edge sector (0,1) is nearly linear (R² = 0.996) because its small dimension (dim = 5) constrains rotation.

---

## Conclusion

The qubit chain responds to light (γ) with a **nonlinear dose-response curve** that varies by sector. The nonlinearity is not a mode-crossing artifact; it is genuine eigenvector rotation driven by Hamiltonian-dissipator competition. Each sector has a characteristic dose-response shape determined by its mode space dimension and the competition between coherent (Hamiltonian) and dissipative (dephasing) dynamics within that sector.

In cavity language: the cavity's internal mode structure rotates as the incident light intensity changes. At low light, the cavity's standing-wave pattern is set by the Hamiltonian. At high light, the pattern is set by the dephasing geometry. The transition between these regimes is smooth but nonlinear.

V1's measurement (134% deviation from linearity) was quantitatively correct. V2 identifies the mechanism: eigenvector rotation within each sector, not mode-crossing between sectors.

---

## Files

- `simulations/light_dose_response.py` (V2 analysis with mode tracking)
- `simulations/results/light_dose_response/light_dose_results.json` (raw data)
- `simulations/results/light_dose_response/light_dose_response.png` (plots)
- [GAMMA_AS_BINDING](GAMMA_AS_BINDING.md) (V1: first measurement, wrong framing)
- [EXCLUSIONS](../docs/EXCLUSIONS.md) (gamma is external, Exclusion 2)
- [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md) (chain is cavity, 4/5 tests)

---

*April 12, 2026. The cavity responds nonlinearly to light because its internal mode structure rotates with light intensity.*
