# IBM Calibration Snapshots

Single-day snapshots of IBM Quantum backend calibration data, used as test
fixtures for the C# `RCPsiSquared.Core.Calibration` module and as static
anchors for analysis scripts.

## Contents

| File | Backend | Date | Qubits |
|------|---------|------|--------|
| `ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv` | ibm_marrakesh (Heron r2) | 2026-04-25 | 156 |
| `ibm_marrakesh_calibrations_2026-04-30T16_25_19Z.csv` | ibm_marrakesh (Heron r2) | 2026-04-30 | 156 |
| `ibm_marrakesh_calibrations_2026-05-05T04_18_17Z.csv` | ibm_marrakesh (Heron r2) | 2026-05-05 | 156 |

These are point-in-time snapshots, distinct from the multi-day time series
in `data/ibm_history/results/ibm_<backend>_history.csv`. The April 25 file
brackets the 2026-04-26 framework_snapshots and soft_break experiment runs
on the high side (1 day before); the April 30 file brackets them on the low
side (4 days after). The pair documents the 5-day drift window discussed in
the Marrakesh path-biography review.

## Source

Originally pulled from the IBM Quantum platform via `target_history()` and
saved as CSV by an external pipeline. These copies are committed to the repo
so the test suite is reproducible without depending on the gitignored
`ClaudeTasks/` workspace.

## Fields

Standard IBM calibration CSV: Qubit ID, T1 (μs), T2 (μs), readout error,
single-qubit gate errors (sx, X), Pauli-X error, CZ neighbour map, RZZ
neighbour map, operational flag. Parsed by `IbmCalibration.Load` in
`compute/RCPsiSquared.Core/Calibration/IbmCalibration.cs`.

## Anchored test values (Apr-25)

- 156 operational qubits
- Best 3-chain: [4, 3, 2], score ≈ 867.07
- Best 5-chain: [1, 2, 3, 4, 5], score ≈ 1246.58
- Documented soft_break path [48, 49, 50] score ≈ 682.50
- Original framework_snapshots path [0, 1, 2] score ≈ 597.27
- Q0: pulse-stable quantum-side (r ≈ 0.074)
- Q1, Q48, Q49, Q50: classical-side (r > 0.4)
