"""IBM Heron r2 calibration loader, qubit scoring, and chain finder.

Loads a calibration CSV (e.g. from `ClaudeTasks/IBM_R2_calibration_ibm_marrakesh/`)
and provides:

- `load_calibration(csv_path)`: parse rows into QubitData records
- `score_qubit(q)`: composite quality score (higher = better) based on T1, T2,
  readout-error, single-qubit-gate-error
- `best_qubits(qubits, k)`: top-k qubits by score
- `parse_couplings(qubits)`: directed-edge coupling graph (qubit → neighbour →
  CZ-error, RZZ-error)
- `best_chain(qubits, length)`: best contiguous path of given length, scored as
  qubit-quality sum + bond-quality sum

Motivation: 2026-04-26 F87 trichotomy run on path [48, 49, 50] gave a 23×
cleaner truly Π²-odd baseline than the same-day run on path [0, 1, 2]
(0.0013 vs 0.0297). The truly baseline is qubit-quality-dependent; chain
selection materially improves signal-to-noise. Calibration-driven chain
choice replaces ad-hoc qubit picking.
"""
from __future__ import annotations

import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class QubitData:
    qubit: int
    t1_us: float
    t2_us: float
    readout_error: float
    sx_error: float
    pauli_x_error: float
    operational: bool
    cz_neighbours: dict[int, float] = field(default_factory=dict)
    rzz_neighbours: dict[int, float] = field(default_factory=dict)


def _parse_neighbour_field(field_str: str) -> dict[int, float]:
    """Parse "1:0.0023;0:0.0024" → {1: 0.0023, 0: 0.0024}."""
    out: dict[int, float] = {}
    if not field_str:
        return out
    for entry in field_str.split(";"):
        if ":" not in entry:
            continue
        nbr, val = entry.split(":")
        try:
            out[int(nbr)] = float(val)
        except ValueError:
            continue
    return out


def _safe_float(value: str, default: float = 0.0) -> float:
    """IBM calibration CSV occasionally has empty fields for non-applicable
    metrics (e.g. SX error on a permanently-down qubit). Parse leniently."""
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def load_calibration(csv_path: Path | str) -> list[QubitData]:
    """Parse IBM calibration CSV into QubitData records. Handles empty fields
    (non-operational qubits sometimes have missing gate-error entries) by
    defaulting to 0."""
    qubits: list[QubitData] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            q = QubitData(
                qubit=int(row["Qubit"]),
                t1_us=_safe_float(row["T1 (us)"]),
                t2_us=_safe_float(row["T2 (us)"]),
                readout_error=_safe_float(row["Readout assignment error"]),
                sx_error=_safe_float(row["√x (sx) error"]),
                pauli_x_error=_safe_float(row["Pauli-X error"]),
                operational=row["Operational"].strip().lower() == "yes",
                cz_neighbours=_parse_neighbour_field(row["CZ error"]),
                rzz_neighbours=_parse_neighbour_field(row["RZZ error"]),
            )
            qubits.append(q)
    return qubits


def score_qubit(q: QubitData) -> float:
    """Composite quality score: rewards high T1/T2, penalises gate + readout
    errors. The dominant term is T2 (decoherence drives Π²-odd noise floor at
    truly Hamiltonians); T1 secondary; gate/readout errors as multiplicative
    penalties. Non-operational qubits return -inf.
    """
    if not q.operational:
        return float("-inf")
    coherence = min(q.t2_us, 2.0 * q.t1_us)
    gate_quality = (1.0 - q.sx_error) ** 4 * (1.0 - q.pauli_x_error) ** 4
    readout_quality = 1.0 - q.readout_error
    return coherence * gate_quality * readout_quality


def best_qubits(qubits: list[QubitData], k: int) -> list[QubitData]:
    """Top-k qubits by score."""
    return sorted(qubits, key=score_qubit, reverse=True)[:k]


def best_chain(qubits: list[QubitData], length: int) -> tuple[list[int], float]:
    """Find best contiguous (CZ-coupled) path of `length` qubits.

    Score = sum of per-qubit scores + sum of -log(CZ-error) along bonds (so
    weaker errors contribute more). Searches all paths via DFS without revisits.
    Returns (path, score). Path is a list of qubit indices in chain order.
    """
    if length < 1:
        raise ValueError("length must be ≥ 1")
    by_id = {q.qubit: q for q in qubits}
    best_path: list[int] = []
    best_score = float("-inf")

    def dfs(path: list[int], visited: set[int], partial: float) -> None:
        nonlocal best_path, best_score
        if len(path) == length:
            if partial > best_score:
                best_score = partial
                best_path = list(path)
            return
        last = by_id[path[-1]]
        for nbr_id, cz_err in last.cz_neighbours.items():
            if nbr_id in visited or nbr_id not in by_id:
                continue
            nbr = by_id[nbr_id]
            if not nbr.operational:
                continue
            bond_score = -1e3 * cz_err  # weaker CZ error → less negative → higher score
            qubit_score = score_qubit(nbr)
            visited.add(nbr_id)
            path.append(nbr_id)
            dfs(path, visited, partial + qubit_score + bond_score)
            path.pop()
            visited.remove(nbr_id)

    for q in qubits:
        if not q.operational:
            continue
        dfs([q.qubit], {q.qubit}, score_qubit(q))
    return best_path, best_score


def chain_score(qubits: list[QubitData], path: list[int]) -> float:
    """Score a specific chain (for comparison against `best_chain`)."""
    by_id = {q.qubit: q for q in qubits}
    score = sum(score_qubit(by_id[p]) for p in path)
    for a, b in zip(path[:-1], path[1:]):
        cz_err = by_id[a].cz_neighbours.get(b, by_id[b].cz_neighbours.get(a, 1.0))
        score += -1e3 * cz_err
    return score


# ---------- demo / self-test ----------

def main():
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    repo_root = Path(__file__).resolve().parents[1]
    calib_path = (
        repo_root
        / "ClaudeTasks"
        / "IBM_R2_calibration_ibm_marrakesh"
        / "ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv"
    )
    print(f"Loading {calib_path.name}")
    qubits = load_calibration(calib_path)
    print(f"  {len(qubits)} qubits, {sum(1 for q in qubits if q.operational)} operational")
    print()

    print("Top-10 qubits by composite score (T2-dominated):")
    print(f"  {'qubit':>6} {'T1 (us)':>9} {'T2 (us)':>9} {'sx err':>10} {'readout':>9} {'score':>10}")
    for q in best_qubits(qubits, 10):
        print(f"  {q.qubit:>6} {q.t1_us:>9.1f} {q.t2_us:>9.1f} {q.sx_error:>10.5f} {q.readout_error:>9.4f} {score_qubit(q):>10.2f}")
    print()

    print("Best 3-qubit contiguous chain:")
    path, score = best_chain(qubits, 3)
    print(f"  path = {path}, score = {score:.2f}")
    print()

    print("Validation against runs on this calibration (2026-04-26 datasets):")
    for label, p in [("good run [48, 49, 50]", [48, 49, 50]), ("first run [0, 1, 2]", [0, 1, 2])]:
        s = chain_score(qubits, p)
        print(f"  {label:<25} score = {s:>8.2f}")
    print()

    print("Best 5-qubit chain (for K_CC_pr / F86 EP-resonance experiments):")
    path5, score5 = best_chain(qubits, 5)
    print(f"  path = {path5}, score = {score5:.2f}")


if __name__ == "__main__":
    main()
