# The reality-window flight, ibm_marrakesh, 2026-10-29

The raw record of the flight pre-registered and read in
[`experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md`](../../experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md):
XY+YX on path [4, 5, 6] from |+−+⟩, where every reading with an odd number
of Y letters is zero for any real channel (up to a small floor from noise
inside the compiled frames and the readout's leak): ten arms of nine bases
each and a readout arm of eight circuits, 98 frozen circuits flown twice at
4096 shots each.

| File | What it is |
|---|---|
| `reality_window_ibm_marrakesh_20261029_053203.json` | job `dbhcohht9uhs73c9vpj0`: the counts of both halves per circuit, the submission order, the day-of snapshot, the job's own timestamps, the sha256 of the frozen circuits |
| `reality_window_20261029_053203.submitted.json` | the stub written before submission (order, snapshot, hash) with the job id added on acceptance |
| `properties_at_run.json` | the backend's properties for q4, q5, q6 and their gates as reported for 05:14 UTC (during the run) and at a later pull, with each property's own date; `identical: false` reflects only date stamps of a 06:25 (+01:00) update after the run |

The circuits are `simulations/flight/reality_window_circuits.qpy`; the
analysis is `simulations/flight/run_reality_window.py --analyze` on the
counts file, and the reading beside the verdicts is
`simulations/reality_window_post_flight.py`.
