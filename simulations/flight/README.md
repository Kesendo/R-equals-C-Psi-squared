# The corner-beat flight pipeline, credential-free

This directory is a **copy for the record**. The live pipeline that talks to
IBM lives outside this repository, deliberately, because that working tree also
holds an API key. The copy exists so the procedure behind a paid hardware run is
readable, runnable and reviewable by anyone who has the repository, including a
future session that no longer has the other machine.

Nothing here contains a credential. There is no key file, no `.env`, and no
token literal: the runner reads `IBM_QUANTUM_TOKEN` from the environment and
stops with a message if it is unset. That is the placeholder, and it is the only
one needed.

## What is here, and what runs without hardware

| File | What it is |
|---|---|
| `run_corner_beat.py` | the flight runner: circuit construction, the local certificates, the Aer path, the day-of gates, submission, analysis |
| `corner_beat_refreeze.py` | the threshold harness: re-measures the registered thresholds from simulation banks |
| `corner_beat_verdict.py` | the verdict composer: applies the committed rules to measured numbers |
| `corner_beat_constants.json` | the frozen constants the runner refuses to start without |
| `test_corner_beat_*.py` | 135 tests over the above |

Three things run with no network and no account, and they are most of what makes
the design checkable:

```bash
cd simulations/flight
python -m pytest test_corner_beat_*.py -q        # 135 tests
python run_corner_beat.py --certify              # local certificates, no network
python run_corner_beat.py --aer                  # the simulator path
```

`--certify` is the interesting one for a reader: it re-derives the Strang step's
sector parity, the prepared dyads, the transpiled circuit shape, the dose
certificates over all 180 swept tables, and a bit-level parity of the fit
against the committed gate at `simulations/corner_beat_gate.py`. It was run from
this copy on 2026-08-18 and passed every check.

`--hardware` needs `IBM_QUANTUM_TOKEN` in the environment (or a `.env` in this
directory, which is gitignored). It costs money and is not repeatable; the
pre-registration governs whether it may be run at all.

## Provenance, and the drift this creates

Copied 2026-08-18 from
`AIEvolution/AIEvolution.UI/experiments/ibm_quantum_tomography`, at these source
hashes (sha256, first 16):

| file | sha256[:16] |
|---|---|
| `run_corner_beat.py` | `ebac22020a5bc398` |
| `corner_beat_refreeze.py` | `71a3e2e08d21a1ba` |
| `corner_beat_verdict.py` | `1454bf028ce11bb0` |
| `corner_beat_constants.json` | `ee603481799eb41d` |
| `test_corner_beat_constants.py` | `a3a825ec888bf4e6` |
| `test_corner_beat_floor.py` | `229b2250191bef7a` |
| `test_corner_beat_verdict.py` | `972015f00b06d6b2` |
| `test_corner_beat_t1_clean.py` | `6fedd36ef59e5cf7` |

**One edit was made to the copy**, so `run_corner_beat.py` will not match its
source hash: `RCPSI_REPO` is derived from `__file__` here (the repository root is
two levels up from this directory) instead of being a hardcoded absolute path to
one machine. Everything else is byte-identical at the time of copying.

**Two copies of a 3473-line runner will drift, and pretending otherwise is how a
record becomes a lie.** The rule for this directory: the copy is the RECORD, the
outside pipeline is the WORKING tree. When the working tree changes, re-copy and
re-state the hashes above in the same commit. A copy whose hashes are stale is
worse than no copy, because it looks authoritative.

## What is deliberately not here

- `results/` (gitignored): flight artifacts, produced by a paid run and carrying
  raw counts and job metadata. The pre-registration cites them by filename and
  sha, never by pasting content.
- the API key file and the billing query script that reads it.
- the other runners in that tree, which belong to other experiments.

## Where the rest of the flight lives

- The pre-registration, which governs everything here:
  [`experiments/CORNER_BEAT_HARDWARE_PREDICTION.md`](../../experiments/CORNER_BEAT_HARDWARE_PREDICTION.md)
- The committed gate the runner checks its fit against:
  [`simulations/corner_beat_gate.py`](../corner_beat_gate.py)
- Open questions about the conventions this runner uses: the arc
  `gamma_book_enforced_nowhere` in
  [`OpenArcsRegistry.cs`](../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs).
  The flight's own arc, `corner_beat`, does not exist yet; it is due with the
  freeze commit.

## A note on the split, since it caused a real problem

The runner reaches back into this repository for the committed gate and for the
pre-registration commit check. That direction is fine. The direction that bit us
is the other one: the runner's own DEFINITIONS, in particular which convention a
rate is written in, live only in the pipeline tree, so the repository's
"look here first" convention could not reach them. One factor-of-two seam was
found this way on 2026-08-18 (arc `gamma_book_enforced_nowhere`). Having the
runner in the repository is part of closing that hole.
