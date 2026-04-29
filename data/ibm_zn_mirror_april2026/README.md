# IBM Marrakesh Z⊗N-Mirror Hardware Test, April 29, 2026

Live hardware verification of the framework's Z⊗N-Mirror diagnostic for
transverse-field detection on IBM Heron r2.

## Files

| File | Description |
|------|-------------|
| `zn_mirror_ibm_marrakesh_20260429_102824.json` | Raw counts and 16 Pauli expectations on (q0=48, q2=50) for both states ρ_a=\|+−+⟩ and ρ_b=\|−+−⟩ |

## Experiment summary

- **Backend:** ibm_marrakesh (Heron r2)
- **Date:** 2026-04-29, 10:28 UTC
- **Job ID:** `d7ornigror3c73c0c6ug`
- **Path:** [48, 49, 50] (same as soft_break April 26)
- **States:** ρ_a = \|+−+⟩ (X-basis Néel), ρ_b = Z⊗N · \|+−+⟩ = \|−+−⟩
- **Hamiltonian:** Heisenberg XX+YY+ZZ, J = 1.0
- **Evolution:** t = 0.8, n_trotter = 3
- **Tomography:** 9 Pauli bases on (q0, q2), 4096 shots/basis
- **Total:** 18 circuits, **23 s billed QPU** (Batch-mode efficient, faster than the per-job session overhead estimate)

## Results

| Quantity | Value |
|----------|-------|
| Max violation | **0.182** (worst Pauli: Z,Z) |
| RMS violation | 0.087 |

### Predictions vs Hardware

| Scenario | Predicted max_violation | Hardware match? |
|----------|------------------------|-----------------|
| Clean Heisenberg + T1 + Tφ + ZZ-crosstalk | < 1e-3 | NO (0.182) |
| Effective h_x = 0.05 (transverse X-field) | ≈ 4e-3 | NO (too small) |
| **Effective h_y = 0.05 (transverse Y-field)** | **≈ 0.18** | **YES (matches 0.182)** |

**Conclusion:** Marrakesh has an effective transverse Y-field of magnitude
h_y_eff ≈ 0.05 at Hamiltonian level, NOT a transverse X-field. This matches
the framework's predicted 40× asymmetry between Y and X (Y is bit_b-odd
like Z, the dephasing axis, so it mixes more strongly).

### Top-5 violating Pauli strings

| Pauli | ⟨P⟩_a | ⟨P⟩_b | expected_b | violation |
|-------|-------|-------|------------|-----------|
| Z,Z | 0.085 | 0.267 | 0.085 | **0.182** |
| Y,Z | 0.252 | -0.096 | -0.252 | 0.156 |
| Y,Y | 0.118 | 0.259 | 0.118 | 0.141 |
| Z,Y | -0.237 | 0.124 | 0.237 | 0.113 |
| X,Y | -0.047 | 0.043 | -0.047 | 0.090 |

The Z,Y / Y,Z asymmetry confirms a Y-axis rotation between the two
states' evolutions — consistent with a single-site h_y field.

## How to reproduce

Hardware-run script: `AIEvolution.UI/experiments/ibm_quantum_tomography/run_zn_mirror.py`
```
python run_zn_mirror.py --hardware --backend ibm_marrakesh --path 48,49,50
```

Analysis: `simulations/_zn_mirror_hardware_analysis.py`
```
python simulations/_zn_mirror_hardware_analysis.py \
       data/ibm_zn_mirror_april2026/zn_mirror_ibm_marrakesh_20260429_102824.json
```

Reference: Confirmation entry `marrakesh_transverse_y_field_detection` in
`framework.Confirmations`. Memory: `project_zn_mirror_diagnostic.md`.
