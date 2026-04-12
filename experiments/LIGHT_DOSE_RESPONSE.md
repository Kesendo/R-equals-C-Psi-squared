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

---

## Context

Gamma is external (EXCLUSIONS Exclusion 2). On IBM hardware, gamma is literally light: microwave photons in a physical resonator cause dephasing through photon shot noise (Sears et al., Phys. Rev. B 86, 180504, 2012). The qubit chain is a passive optical cavity (4/5 cavity tests, R^2 = 0.998 for beam profile).

V1 (GAMMA_AS_BINDING.md) measured per-sector rates and found 134% deviation from linear scaling. V1's numbers were correct but its framing (gamma as neutral "binding parameter") ignored the established status of gamma as light. V2 reframes: gamma is light, the table is a light-dose response curve, and the question is what mechanism produces the nonlinearity.

---

## Frage 1: Light dose per sector (sacrifice profile, alpha = 1.0)

The sacrifice profile concentrates light on site 0 (gamma_0 = 2.34, quiet sites gamma_k ~ 0.05-0.10). The per-sector response:

| Sector type | Example | Slowest rate | Interpretation |
|-------------|---------|-------------|----------------|
| SE diagonal (1,1) | 0.318 | Lens mode: slowest in the system |
| Interior diagonal (2,2) | 0.375 | Slightly faster than SE |
| Adjacent coherence (0,1) | 0.167 | Cross-sector coherence, slowest off-diagonal |
| Distant coherence (0,5) | 5.216 = Sigma_gamma | Maximum rate: total light dose |
| Extremal diagonal (0,0) | 0.000 | Stationary (vacuum, nothing to decay) |

The sacrifice geometry protects the SE sector (rate 0.318) by concentrating light away from the cavity's slow standing-wave mode. This is the lens effect quantified per sector.

---

## Frage 2: Linear vs nonlinear light response

### V1 finding: 134% deviation from linearity

V1 scaled all gamma_k by alpha and found per-sector rates deviate up to 134% from the linear prediction rate(alpha) = alpha * rate(1).

### V2 finding: mechanism is eigenvector rotation, not mode-crossing

V2 tracked individual eigenvalue curves across alpha in [0.1, 5.0] with 50 steps. Two hypotheses tested:

**Hypothesis A (mode-crossing):** Individual modes are linear, but different modes become the "slowest" at different alpha, making the minimum nonlinear.

**Hypothesis B (eigenvector rotation):** Individual modes are themselves nonlinear because the Hamiltonian-dissipator balance changes their Pauli content.

### Evidence

| Sector | Per-mode R^2 (mean) | Per-mode R^2 (min) | Sector-min R^2 | Level crossings |
|--------|---------------------|---------------------|----------------|-----------------|
| (1,1) SE | 0.889 | 0.063 | 0.404 | 0 |
| (2,2) interior | 0.801 | 0.001 | 0.144 | 0 |
| (0,1) edge | 0.935 | 0.807 | 0.996 | 0 |
| (1,2) cross | 0.876 | 0.001 | 0.967 | 0 |
| (2,3) cross | 0.798 | 0.006 | 0.785 | 0 |

**Zero level crossings in all sectors.** The slowest mode maintains its identity across the full alpha range. Hypothesis A is ruled out.

**Individual modes are themselves nonlinear** (R^2 as low as 0.001). Hypothesis B is confirmed.

### Mechanism: eigenvector rotation

The Liouvillian is L(alpha) = L_H + alpha * L_D. As alpha increases:

1. At small alpha (Hamiltonian-dominated): eigenvectors resemble the Hamiltonian's normal modes. These have mixed Pauli content, so the absorption theorem rate Re(lambda) = -2 alpha sum gamma_k <1_XY(k)> is sensitive to the eigenvector's XY-weight.

2. At large alpha (dephasing-dominated): eigenvectors rotate toward the Pauli basis (each mode approaches a definite Pauli string). The XY-weight <1_XY(k)> stabilizes.

3. The transition between these regimes is continuous (no crossings), but the eigenvector rotation changes <1_XY(k)> nonlinearly in alpha. This makes Re(lambda) nonlinear even for a fixed mode.

The interior sectors (2,2) and (3,3) show the strongest nonlinearity (R^2 = 0.14) because they have the largest mode spaces and the most room for eigenvector rotation. The edge sectors (0,1) are nearly linear (R^2 = 0.996) because their small dimension constrains rotation.

---

## Conclusion

The qubit chain responds to light (gamma) with a **nonlinear dose-response curve** that varies by sector. The nonlinearity is not a mode-crossing artifact; it is genuine eigenvector rotation driven by Hamiltonian-dissipator competition. Each sector has a characteristic dose-response shape determined by its mode space dimension and the competition between coherent (Hamiltonian) and dissipative (dephasing) dynamics within that sector.

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
