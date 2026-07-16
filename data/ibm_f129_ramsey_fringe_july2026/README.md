# The F129 standing-fringe campaign, ibm_kingston, 2026-07-15

The persisted counts of the pre-registered F129 flight: the n = 9 clean-clean level
collision measured as a standing Ramsey fringe (slope ≈ 0) against a winding clock,
job `d9br4vmg26ic73dgbgk0`, ONE Batch, 64 PUBs × 16384 shots, billed ≈ 5.0 QPU min.
Pre-registration (v2.5, frozen at commit 5bf3b79), the compiler and 7a/7b gate
records, the flight RECORD, and the three-clause verdict:
[experiments/IBM_F129_RAMSEY_FRINGE.md](../../experiments/IBM_F129_RAMSEY_FRINGE.md).
Runner (external pipeline): `AIEvolution.UI/experiments/ibm_quantum_tomography/run_f129_ramsey_fringe.py`
(modes `--calibrate/--certify/--simulate/--hardware/--recover/--analyze`).

| File | What |
|------|------|
| `f129_flight_20260715_222348.json` | The flight: raw counts for all 64 PUBs (arms A0/A1/A2 × 9 steps in X/Y, CAL 0/1, T1 + RAMSEY at τ = 20/60/120/200 µs), plus before/after calibration snapshots and the `recovered: true` flag (fetched via `--recover` by the submit-record job ID after a local power failure). Chain [11,12,13,14,15,19,35,34], 16384 shots. |
| `f129_submit_20260715_181658.json` | The submit record written at 18:16:58 before the shot: job ID, backend, chain, the 64 PUB keys in submission order, shots, and the pre-flight calibration snapshot. This is the file the `--recover` path keyed on when the waiting process died mid-queue. |

**Campaign verdict, one line:** STANDING FRINGE CONFIRMED. The clean n = 9 mirror
collision (A0) stands still to its computed θ³ drift and budgeted systematics
(measured slope +0.0326 ± 0.0033, inside the budgeted window; the excess is positive
and sized within b_zz2 + b_qs, the direction the ζ² anti-protection law predicts),
while the winding clock A2 turns at 61.0σ and the opposite-winding dial A1 tracks its
prediction at 0.7σ; the nearest impostor triple (predicted +0.109) lies ~22σ above the
measured A0. θ̂ = 0.4992 (the device hit the pinned angle to 0.16 %). Confirmation
entry 24 in both registries (Python `Confirmations`, C# `ConfirmationsRegistry`).
