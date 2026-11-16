# The drive-frame flight, ibm_fez, ibm_marrakesh, ibm_kingston, 2026-11-13 to 2026-11-16

The raw record of the flight pre-registered and read in
[`experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md`](../../experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md)
(arm A, the October XY+YX circuit with one qubit's drive frame turned, a
readout arm and the forest qubits read beside the path) and, on Marrakesh
and Kingston, in
[`experiments/X3_WINDOW_HARDWARE_PREDICTION.md`](../../experiments/X3_WINDOW_HARDWARE_PREDICTION.md)
(arm B, the April YZ+ZY circuit under the same turns, and arm Q, each path
qubit's idle precession). One job per chip, each circuit twice at 4096
shots, submitted from commit `39defc7c`; §7 of both pages holds the
verdicts and the readings.

| chip, path | job | PUBs | billed |
|---|---|---|---|
| ibm_fez [123, 136, 143] | `dbrmdepabdes73fq8jm0` | 124 | 134 s |
| ibm_marrakesh [45, 46, 47] | `dbrmdf1abdes73fq8jmg` | 280 | 307 s |
| ibm_kingston [13, 14, 15] | `dbrmdlgk09pc73av8gog` | 280 | 304 s |

Per chip, with `<stamp>` the submission time in the file name:

| File | What it is |
|---|---|
| `drive_frame_<chip>_<stamp>.submitted.json` | the stub written before submission (the shuffled order, the day-of snapshot, the hashes of the frozen circuits and of the freeze record) with the job id added on acceptance |
| `drive_frame_<chip>_<stamp>.json` | the counts of both halves per circuit, the submission order, the freeze snapshot and the day-of snapshot, the result's metadata, the job's own timestamps and metrics, the billed usage, and the comparison of the retrieved circuits with the frozen ones (`turn_ok`) |
| `drive_frame_<chip>_<stamp>.analysis.json` | the analysis of §5 of both pages on those counts |
| `drive_frame_<chip>_<stamp>.analysis.txt` | the runner's printout of the same |

Beside them, `backend_versions_read_afterwards.json`: each backend's
version and the date of the calibration its properties carry, queried
from IBM's property history on 2026-11-16 for 2026-11-01, 2026-11-10,
2026-11-13 19:00 UTC, the middle of each run and the time of the query.
The record itself does not carry the version.

The frozen circuits and freeze records are
`simulations/flight/drive_frame_<chip>_circuits.qpy` and
`drive_frame_<chip>_freeze.json`; the analysis is
`simulations/flight/run_drive_frame.py --analyze <counts file> --reach`
(the Kingston job, which waited 60.5 h in the queue, was collected with
`--recover` on its stub). The record keeps the gate-for-gate comparison of
the retrieved circuits with the frozen ones, not the retrieved circuits.
Two post-flight readings: the colouring X^⊗3 behind the window page's
§5, `simulations/flight/drive_frame_colouring_reading.py` (output in
`simulations/results/drive_frame/<chip>_colouring_reading.txt`), and the
gate errors the replay leaves out with the baseline at the day's
coherence, `simulations/flight/drive_frame_gate_error_reading.py` (output
in `simulations/results/drive_frame/gate_error_reading.txt`).
