# IBM Run 3 Data: Palindrome Validation (March 18, 2026)

## Files

| File | Description |
|------|-------------|
| `palindrome_predictions.json` | Locked predictions BEFORE hardware run. CΨ trajectory, crossing time, delay points. |
| `palindrome_ibm_torino_20260318_191348.json` | Hardware results: 8-point tomography on Q80, density matrices, CΨ analysis. |
| `palindrome_simulator.json` | Simulator comparison (Aer with same T1/T2*). |
| `palindrome_validation_ibm_torino.png` | Comparison plot: measured vs predicted CΨ(t). |
| `ramsey_sameday_20260318.json` | Same-day Ramsey T2* measurement (Q80: 17.36 μs, Q102: 26.4 μs). |
| `ramsey_march12_20260312.json` | March 12 Ramsey T2* for drift comparison (Q80: 11.0 μs, Q102: 15.4 μs). |

## Key Result

CΨ=1/4 crossing validated at **1.9% deviation** with same-day T2*.

| Parameter | Value |
|-----------|-------|
| Backend | ibm_torino (Heron r2) |
| Qubit | 80 (permanent crosser) |
| T1 (day-of) | 143.1 μs |
| T2* (same-day Ramsey) | 17.36 μs |
| Predicted crossing | t* = 15.01 μs |
| Measured crossing | t* = 15.29 μs |
| Deviation | 0.28 μs = 1.9% |

## T2* Drift

T2* increased 58% in 6 days (11.0 → 17.36 μs), explaining the initial
61% deviation when using the stale March 12 value. Same-day Ramsey
is essential for accurate predictions.
