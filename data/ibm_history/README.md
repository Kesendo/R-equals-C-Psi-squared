# IBM Torino Calibration History

## What is this?

Daily calibration data for all 133 qubits on IBM's Torino quantum
processor (Eagle r3), from August 14, 2025 to February 10, 2026.

IBM publishes these values for every backend on the IBM Quantum Platform.
They are measured automatically by the system, not by us.

## Files

| File | Content |
|------|---------|
| `ibm_torino_history.csv` | 24,074 records: daily T1, T2, frequency, r\_param, CΨ metrics for 133 qubits |
| `ibm_history_analysis.py` | Python script to fetch and analyze the data (requires IBM Quantum account) |

## CSV Columns

| Column | Meaning |
|--------|---------|
| date | Calibration date (YYYY-MM-DD) |
| qubit | Qubit index (0-132) |
| T1\_us | Energy relaxation time in microseconds |
| T2\_us | Coherence time in microseconds |
| frequency\_GHz | Qubit transition frequency |
| r\_param | T2 / (2 × T1), the ratio that determines CΨ crossing |
| crosses\_quarter | True if r < 0.213 (qubit crosses CΨ = ¼ boundary) |

## How to read it

The key column is **r\_param**. It determines whether a qubit is in the
quantum regime (r < 0.213, CΨ crosses ¼) or the classical regime (r > 0.213).

Most qubits have r between 0.3 and 0.8 and never cross. But 16 qubits
oscillate: their r fluctuates around 0.213, crossing some days and not
others. These are the qubits visible in
[Both Sides Visible](../../docs/BOTH_SIDES_VISIBLE.md).

## How to reproduce

```bash
# Requires: pip install qiskit-ibm-runtime pandas
# Requires: IBM Quantum account token
export IBM_QUANTUM_TOKEN=your_token_here
python ibm_history_analysis.py
```

The script fetches the latest calibration data from IBM's API and
computes the derived columns (r\_param, crosses\_quarter).

## Key findings

- **Q98** (57.5% crossing rate): clear lifecycle pattern (tune, pulse, fade)
- **Q72** (66.9%): balanced rhythmic alternation
- **Q105** (56.9%): long active phase, then silent
- **Q80** (validated in [IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md)):
  consistently crosses (r ≈ 0.07-0.11), 1.9% deviation from theory

See [Both Sides Visible](../../docs/BOTH_SIDES_VISIBLE.md) for the
palindromic complement analysis.
