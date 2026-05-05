"""Cross-backend Heron-r2 comparison: Marrakesh vs Kingston vs Fez on the
2026-05-05 snapshots.

Three Heron-r2 chips, 156 qubits each, same CZ-coupling pattern (heavy-hex
honeycomb). What changes between them is the per-qubit T1/T2 distribution
and the spatial pattern of which qubits are stably quantum-side. This
script asks four questions:

1. **Score distribution**: how do the chips compare on aggregate
   per-qubit quality? Median + top-10 + bottom-10 score percentiles.

2. **Stable-quantum population**: how many qubits sit at r < R* on each
   chip (instantaneous snapshot proxy for the 91-day PulseStable archetype)?

3. **Topology-addressable stable-quantum chains**: for each chip, find
   CZ-coupled pairs and triples among the stable-quantum qubits. This
   is the addressable-uniform-quantum question we ran on Marrakesh:
   does Kingston or Fez open chains we couldn't build there?

4. **Best 3- and 5-chains by score**: today's pick on each chip,
   addressable + stable. Candidates for the next hardware run.

Hardware access via AIEvolution (external pipeline). This script informs
the chip-and-path decision; the actual job submission lives elsewhere.
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ibm_calibration import load_calibration, score_qubit, best_chain, chain_score
from _qubit_biography import R_STAR

REPO_ROOT = Path(__file__).resolve().parents[1]
SNAP_DIR = REPO_ROOT / "data" / "ibm_calibration_snapshots"
SNAPS = {
    "marrakesh": SNAP_DIR / "ibm_marrakesh_calibrations_2026-05-05T04_18_17Z.csv",
    "kingston":  SNAP_DIR / "ibm_kingston_calibrations_2026-05-05T07_48_33Z.csv",
    "fez":       SNAP_DIR / "ibm_fez_calibrations_2026-05-05T08_00_06Z.csv",
}


def stable_quantum(qubits, r_threshold=R_STAR):
    """Qubits with r = T2/(2*T1) < R_STAR (instantaneous quantum-side)."""
    return [q for q in qubits if q.t1_us > 0 and (q.t2_us / (2 * q.t1_us)) < r_threshold]


def cz_pairs_among(quantum_qubits, full_qubits):
    """CZ-coupled pairs where both qubits are in the stable-quantum set."""
    quantum_ids = {q.qubit for q in quantum_qubits}
    by_id = {q.qubit: q for q in full_qubits}
    pairs = set()
    for q in quantum_qubits:
        for nbr_id in by_id[q.qubit].cz_neighbours:
            if nbr_id in quantum_ids and nbr_id != q.qubit:
                pair = tuple(sorted([q.qubit, nbr_id]))
                pairs.add(pair)
    return sorted(pairs)


def cz_triples_among(quantum_qubits, full_qubits):
    """CZ-coupled paths of three qubits where all are stable-quantum."""
    quantum_ids = {q.qubit for q in quantum_qubits}
    by_id = {q.qubit: q for q in full_qubits}
    triples = []
    seen = set()
    for q in quantum_qubits:
        for n1 in by_id[q.qubit].cz_neighbours:
            if n1 not in quantum_ids:
                continue
            for n2 in by_id[n1].cz_neighbours:
                if n2 == q.qubit or n2 not in quantum_ids:
                    continue
                key = (q.qubit, n1, n2)
                rev = (n2, n1, q.qubit)
                if key in seen or rev in seen:
                    continue
                seen.add(key)
                triples.append([q.qubit, n1, n2])
    return triples


def main():
    print("Cross-backend Heron r2 comparison (2026-05-05 snapshots)")
    print("=" * 78)

    loaded = {name: load_calibration(p) for name, p in SNAPS.items()}

    # --- Section 1: aggregate score distribution ---
    print()
    print("1. AGGREGATE PER-QUBIT SCORE DISTRIBUTION")
    print()
    print(f"  {'backend':<12} {'qubits':>7} {'op':>4} {'median':>9} {'top-10':>9} {'bot-10':>9} {'range':>14}")
    print("  " + "-" * 70)
    for name, qubits in loaded.items():
        scores = sorted([score_qubit(q) for q in qubits if q.operational], reverse=True)
        op = sum(1 for q in qubits if q.operational)
        median = scores[len(scores) // 2]
        top10 = np.mean(scores[:10])
        bot10 = np.mean(scores[-10:])
        rng = f"[{scores[-1]:.0f} .. {scores[0]:.0f}]"
        print(f"  {name:<12} {len(qubits):>7} {op:>4} {median:>9.2f} {top10:>9.2f} {bot10:>9.2f} {rng:>14}")

    # --- Section 2: stable-quantum population ---
    print()
    print(f"2. STABLE-QUANTUM POPULATION (r < R* = {R_STAR:.4f})")
    print()
    print(f"  {'backend':<12} {'sq count':>9} {'sq fraction':>12}")
    print("  " + "-" * 40)
    sq_per_backend = {}
    for name, qubits in loaded.items():
        sq = stable_quantum(qubits)
        sq_per_backend[name] = sq
        frac = len(sq) / len(qubits) * 100
        print(f"  {name:<12} {len(sq):>9} {frac:>11.1f}%")

    # --- Section 3: topology-addressable stable-quantum CZ structures ---
    print()
    print("3. STABLE-QUANTUM CZ-COUPLED STRUCTURES PER BACKEND")
    print()
    print(f"  {'backend':<12} {'sq pairs':>9} {'sq triples':>11}")
    print("  " + "-" * 40)
    for name, qubits in loaded.items():
        sq = sq_per_backend[name]
        pairs = cz_pairs_among(sq, qubits)
        triples = cz_triples_among(sq, qubits)
        print(f"  {name:<12} {len(pairs):>9} {len(triples):>11}")
        for p in pairs[:5]:
            scores = [score_qubit(next(q for q in qubits if q.qubit == qid)) for qid in p]
            print(f"    pair  {list(p)}    scores [{scores[0]:.0f}, {scores[1]:.0f}]")
        for t in triples[:5]:
            scores = [score_qubit(next(q for q in qubits if q.qubit == qid)) for qid in t]
            print(f"    triple {t}  scores [{scores[0]:.0f}, {scores[1]:.0f}, {scores[2]:.0f}]")
        if len(pairs) > 5 or len(triples) > 5:
            print(f"    ... (showing first 5 of {len(pairs)} pairs, {len(triples)} triples)")

    # --- Section 4: best 3- and 5-chains per backend ---
    print()
    print("4. BEST 3- AND 5-CHAIN BY SCORE PER BACKEND")
    print()
    print(f"  {'backend':<12} {'best 3-chain':<22} {'score':>9}  {'best 5-chain':<22} {'score':>9}")
    print("  " + "-" * 90)
    for name, qubits in loaded.items():
        p3, s3 = best_chain(qubits, 3)
        p5, s5 = best_chain(qubits, 5)
        p3s = str(list(p3))
        p5s = str(list(p5))
        print(f"  {name:<12} {p3s:<22} {s3:>9.2f}  {p5s:<22} {s5:>9.2f}")

    # --- Section 5: known experimental paths re-evaluated on each backend ---
    print()
    print("5. NAMED EXPERIMENTAL PATHS RE-SCORED ON EACH BACKEND")
    print()
    print(f"  {'path':<28} {'backend':<12} {'score':>9} {'note'}")
    print("  " + "-" * 78)
    named = [
        ("[0, 1, 2]", [0, 1, 2], "framework_snapshots Marrakesh"),
        ("[48, 49, 50]", [48, 49, 50], "soft_break / zn_mirror Marrakesh"),
        ("[14, 15]", [14, 15], "f25_cusp_trajectory Kingston"),
        ("[4, 5, 6]", [4, 5, 6], "f83_pi2_class_signature Marrakesh"),
    ]
    for label, path, note in named:
        for name, qubits in loaded.items():
            try:
                s = chain_score(qubits, path)
            except (KeyError, AttributeError):
                s = float("nan")
            print(f"  {label:<28} {name:<12} {s:>9.2f} {note if name == 'marrakesh' else ''}")

    # --- Section 6: reading ---
    print()
    print("=" * 78)
    print("READING")
    print()
    sq_counts = {name: len(sq) for name, sq in sq_per_backend.items()}
    print(f"  Stable-quantum populations: marrakesh {sq_counts['marrakesh']}, "
          f"kingston {sq_counts['kingston']}, fez {sq_counts['fez']}")
    print(f"  All three are 156-qubit Heron r2 chips with the same heavy-hex CZ graph,")
    print(f"  so structural differences come from per-qubit T1/T2 variation, not topology.")
    print()
    print(f"  Pick logic for the next hardware run:")
    print(f"    A. continuity (same path as prior runs)  → marrakesh [48, 49, 50]")
    print(f"    B. best current score on any chip        → see Section 4")
    print(f"    C. uniform-quantum experiment            → see Section 3 pairs/triples")


if __name__ == "__main__":
    main()
