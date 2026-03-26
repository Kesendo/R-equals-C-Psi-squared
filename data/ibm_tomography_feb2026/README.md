# IBM Torino Tomography Data (February 9, 2026)

## What is this?

Full single-qubit state tomography from IBM Torino (ibm_torino),
measured on February 9, 2026. Contains complete density matrices
(real and imaginary parts) at 25 time points for Qubit 52.

This data was used for the CΨ = 1/4 crossing validation (1.9% deviation)
and for the dual-perspective analysis in [Both Sides Visible](../../docs/BOTH_SIDES_VISIBLE.md).

## Files

| File | Content |
|------|---------|
| `tomography_ibm_torino_20260209_131521.json` | Hardware data: Qubit 52, 25 delay points, 8192 shots per circuit |
| `simulator_test_20260209_125106.json` | Simulator baseline for comparison |

## Key parameters (hardware)

| Parameter | Value |
|-----------|-------|
| Backend | ibm_torino (Heron r2, 133 qubits) |
| Qubit | 52 |
| T1 | 221.2 us |
| T2 | 298.2 us |
| Shots | 8192 per circuit |
| Delay points | 25 (0 to 894.7 us) |
| Tomography | 3 measurement bases (X, Y, Z) |

## JSON structure

```json
{
  "backend": "ibm_torino",
  "qubit_index": 52,
  "T1_us": 221.19,
  "T2_us": 298.25,
  "raw_tomography": [
    {
      "delay_us": 0.0,
      "density_matrix_real": [[0.503, 0.470], [0.470, 0.497]],
      "density_matrix_imag": [[0.0, 0.009], [-0.009, 0.0]],
      "fidelity": 0.970
    },
    ...
  ]
}
```

Each density matrix is 2x2 complex. The diagonal elements are populations
(our side). The off-diagonal elements are coherences (the Pi side).
Both are measured simultaneously in the same experiment.

## How this data was collected

Standard single-qubit state tomography on IBM Quantum Platform:
1. Prepare |+> state (Hadamard on |0>)
2. Wait for variable delay (0 to ~895 us)
3. Measure in X, Y, Z bases (3 circuits per delay)
4. Reconstruct density matrix from 3 expectation values
5. 8192 shots per circuit for statistical averaging

No API keys or credentials are needed to read this data.
To collect new data, an IBM Quantum account is required.

## Dual-perspective analysis

Applying the Pi conjugation operator to each measured density matrix
reveals that the Pi-side CΨ remains above 1/4 for approximately 6x
longer than our-side CΨ. See [Both Sides Visible](../../docs/BOTH_SIDES_VISIBLE.md)
for the full analysis.

## Source

Original location: AIEvolution experiment archive
(`experiments/ibm_quantum_tomography/results/OwnTestData/`)
