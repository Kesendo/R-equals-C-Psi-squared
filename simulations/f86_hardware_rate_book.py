"""Read-only current-truth loader for the 2026-05-31 Kingston population records."""
from __future__ import annotations

import json
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "ibm_ep_onset_may2026"
HARDWARE_SCAN = DATA_DIR / "ep_onset_hardware_ep_ibm_kingston_20260531_064022.json"
TWIRL_SCAN = DATA_DIR / "ep_onset_simulate_twirl_20260531_063048.json"
PART_A = DATA_DIR / "ep_onset_hardware_ibm_kingston_20260531_060943.json"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def population_scan(path: Path = HARDWARE_SCAN, stored_key: str = "revival") -> tuple[list[float], list[float], list[float]]:
    """Return Q_label, canonical Q_Lindblad=2Q_label, and recomputed max n0(t>=2 us)."""
    rows = _read(path)["scan"]
    q_label = [float(row["Q"]) for row in rows]
    revival = []
    for row in rows:
        measured = max(float(pops[0]) for t, pops in row["pops"].items() if float(t) >= 2.0)
        if stored_key in row and abs(measured - float(row[stored_key])) > 1e-12:
            raise ValueError(f"stored {stored_key} disagrees with populations at Q_label={row['Q']}")
        revival.append(measured)
    return q_label, [2.0 * q for q in q_label], revival


def part_a_endpoint() -> tuple[float, list[float]]:
    record = _read(PART_A)
    return float(record["ts_us"][-1]), [float(x) for x in record["pops"][-1]]
