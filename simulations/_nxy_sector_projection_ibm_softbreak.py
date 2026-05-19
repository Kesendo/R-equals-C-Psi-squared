"""n_XY-graded sector projection of the in-repo IBM soft-break tomography.

The admixture lens (`experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md`,
`reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md`) says the chain slow mode is
95% n_XY=0 + 5% n_XY=2 in the Pauli-basis weight distribution. This script
asks whether the F87 / Klein-Vierergruppe Hamiltonian categories that
Marrakesh + Kingston measured leave a corresponding n_XY signature in the
2-qubit reduced (q0, q2) Pauli tomography.

Per Pauli string `P = (P_0, P_2)` on 2 qubits, n_XY counts the number of
X or Y letters: I → 0, Z → 0, X → 1, Y → 1. Three bins per category:

  bin n_XY=0: {I,I; I,Z; Z,I; Z,Z}                       (4 strings, "dark")
  bin n_XY=1: {I,X; I,Y; Z,X; Z,Y; X,I; X,Z; Y,I; Y,Z}   (8 strings, "mixed")
  bin n_XY=2: {X,X; X,Y; Y,X; Y,Y}                       (4 strings, "light")

For each Hamiltonian category in each JSON, this script reports:
  - Σ |⟨P⟩|² per bin (the operator-content weight in that bin)
  - RMS ⟨P⟩ per bin (mean magnitude per string, normalised by bin size)
  - Bin 2 / bin 0 ratio (relative "light" vs "dark" weight)
  - The single largest contributor per bin

The admixture-lens prediction:
  - truly_unbroken (XX+YY): H preserves popcount; weight concentrates in
    bin 0 (the "memory" content) + bin 2 (the admixture-channel content)
  - hard_broken (XX+XY): H mixes parities; bin 1 should pick up weight
  - pi2_odd_pure: maximal admixture coupling; bin 2 should peak
  - pi2_even_nontruly: less admixture coupling than truly

If the n_XY-graded bins show clear separation by F87 category, the
admixture lens has a hardware shadow on the existing tomography data
(no new QPU minutes needed). If not, the connection is more subtle than
"bin-grade weight tracks F87 class".
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def n_xy_of(letter: str) -> int:
    return 1 if letter in ("X", "Y") else 0


def n_xy_pair(key: str) -> int:
    a, b = key.split(",")
    return n_xy_of(a) + n_xy_of(b)


def analyse(filepath: str) -> None:
    with open(filepath, "r", encoding="utf-8-sig") as f:
        d = json.load(f)
    print(f"File: {filepath}")
    print(f"  backend: {d['backend']}, path: {d['path']}, t={d['parameters']['t']}, "
          f"shots: {d['parameters']['shots']}")
    print()

    categories = list(d["expectations"].keys())
    # Header
    print(f"  {'Category':>22}  | "
          f"{'Σ|P|² n=0':>10} {'Σ|P|² n=1':>10} {'Σ|P|² n=2':>10} | "
          f"{'RMS n=0':>8} {'RMS n=1':>8} {'RMS n=2':>8} | "
          f"{'n=2/n=0':>8}")
    print("  " + "-" * 110)
    for cat in categories:
        bins_w = {0: 0.0, 1: 0.0, 2: 0.0}
        bins_count = {0: 0, 1: 0, 2: 0}
        for key, val in d["expectations"].items() if False else d["expectations"][cat].items():
            n = n_xy_pair(key)
            bins_w[n] += val * val
            bins_count[n] += 1
        # Normalise RMS per bin (divide weight by bin size, then sqrt)
        rms = {n: (bins_w[n] / bins_count[n]) ** 0.5 for n in (0, 1, 2)}
        # Exclude II from bin 0 weight for the n=2/n=0 ratio (since II=1 always)
        bin0_no_ii = bins_w[0] - 1.0
        ratio = bins_w[2] / bin0_no_ii if bin0_no_ii > 1e-9 else float("nan")
        print(f"  {cat:>22}  | "
              f"{bins_w[0]:10.4f} {bins_w[1]:10.4f} {bins_w[2]:10.4f} | "
              f"{rms[0]:8.4f} {rms[1]:8.4f} {rms[2]:8.4f} | "
              f"{ratio:8.4f}")

    # Largest contributor per bin per category
    print()
    print("  Top-magnitude Pauli pair per bin per category:")
    print(f"  {'Category':>22}  | {'bin n=0 top':>22} {'bin n=1 top':>22} {'bin n=2 top':>22}")
    print("  " + "-" * 100)
    for cat in categories:
        exps = d["expectations"][cat]
        per_bin = {0: [], 1: [], 2: []}
        for key, val in exps.items():
            n = n_xy_pair(key)
            per_bin[n].append((abs(val), key, val))
        tops = {}
        for n in (0, 1, 2):
            per_bin[n].sort(reverse=True)
            # Skip II for bin 0
            for absval, key, val in per_bin[n]:
                if key == "I,I":
                    continue
                tops[n] = f"{key}={val:+.3f}"
                break
            else:
                tops[n] = "—"
        print(f"  {cat:>22}  | {tops[0]:>22} {tops[1]:>22} {tops[2]:>22}")
    print()


def main() -> None:
    print("=" * 110)
    print("n_XY-graded sector projection of IBM soft-break + Klein-Vierergruppe tomography")
    print("=" * 110)
    print()
    print("Bin definitions:")
    print("  bin n=0 (dark):  {I,I; I,Z; Z,I; Z,Z}          — pure population content")
    print("  bin n=1 (mixed): one X/Y letter, one I/Z letter — off-diagonal popcount sector")
    print("  bin n=2 (light): {X,X; X,Y; Y,X; Y,Y}          — admixture-channel content")
    print()

    files = [
        "data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json",
        "data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json",
    ]
    for fp in files:
        if Path(fp).exists():
            analyse(fp)
            print()
        else:
            print(f"  [SKIP] {fp} not found")

    print("=" * 110)
    print("Admixture-lens predictions to check against the table:")
    print("  truly_unbroken (XX+YY):     n=0 and n=2 weights both substantial, n=1 small")
    print("  soft_broken / pi2_odd_pure: n=2 weight maximised (cleanest admixture coupling)")
    print("  hard_broken / pi2_even:     n=1 picks up weight (parity mixing visible)")
    print("  mixed_anti_one_sixth:       similar to hard or pi2_even pattern expected")


if __name__ == "__main__":
    main()
