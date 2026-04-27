#!/usr/bin/env python3
"""Y-parity protection verification on existing 7 IBM hardware snapshots.

The 4-panel cockpit (commit 8496344) identified Y-parity as a universal
algebraic protection: when L preserves Y-parity AND ρ_0 has no Y-odd
content, the 28 Y-odd Pauli observables stay at zero forever.

For Snapshot D's three Hamiltonians at N=3, |+−+⟩, the y_parity_panel
predicts:
  truly XX+YY  → Y-parity NOT preserved (XX has 0 Y, YY has 2 Y, mixed)
                 → some Y-odd cells should be NONZERO on hardware
  soft  XY+YX  → Y-parity PRESERVED ('1 Y per term') → all 6 Y-odd cells
                 should be in noise band
  hard  XX+XY  → Y-parity NOT preserved (XX has 0 Y, XY has 1 Y, mixed)
                 → some Y-odd cells should be NONZERO

The 6 Y-odd cells on (q0, q2) are: IY, XY, ZY, YI, YX, YZ
The 10 Y-even cells are: II, IX, IZ, XI, XX, XZ, YY, ZI, ZX, ZZ

This script:
  1. Loads all 7 framework_snapshots_ibm_*.json
  2. For each (backend, category), extract the 6 Y-odd cells
  3. Compute max |⟨P⟩| per category, compared to noise threshold ~0.2
  4. Verify: soft has small max, truly + hard have large max somewhere
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

RESULTS_DIR = Path(
    r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI"
    r"\experiments\ibm_quantum_tomography\results"
)

Y_ODD_KEYS = ['I,Y', 'X,Y', 'Z,Y', 'Y,I', 'Y,X', 'Y,Z']
Y_EVEN_KEYS = ['I,X', 'I,Z', 'X,I', 'X,X', 'X,Z', 'Y,Y', 'Z,I', 'Z,X', 'Z,Z']

NOISE_THRESHOLD = 0.20  # generous noise band for Heron r2 tomography


def load_snapshot_d(json_path):
    with open(json_path, encoding='utf-8') as f:
        d = json.load(f)
    if 'snapshot_d_softbreak_trichotomy' not in d:
        return None
    backend = d.get('backend', json_path.stem)
    ts = json_path.stem.split('_')[-1]
    return {
        'backend': backend,
        'timestamp': ts,
        'key': f"{backend}_{ts}",
        'expectations': d['snapshot_d_softbreak_trichotomy']['expectations_per_category'],
    }


def main():
    files = sorted(RESULTS_DIR.glob("framework_snapshots_ibm_*.json"))
    snapshots = []
    for f in files:
        snap = load_snapshot_d(f)
        if snap is not None:
            snapshots.append(snap)

    print(f"Loaded {len(snapshots)} hardware Snapshot-D runs")
    print()
    print(f"Y-parity prediction:")
    print(f"  truly XX+YY  →  Y-parity NOT preserved → some Y-odd cells nonzero")
    print(f"  soft  XY+YX  →  Y-parity PRESERVED     → all 6 Y-odd cells small")
    print(f"  hard  XX+XY  →  Y-parity NOT preserved → some Y-odd cells nonzero")
    print()
    print(f"Hardware noise threshold: |⟨P⟩| < {NOISE_THRESHOLD} considered 'in band'")
    print()

    print(f"{'backend':<22s}  {'cat':<6s}  ", end='')
    for k in Y_ODD_KEYS:
        print(f"{k:>7s}  ", end='')
    print(f"{'max|Y-odd|':>10s}  {'verdict':<24s}")
    print('-' * 145)

    summary = {'truly': [], 'soft': [], 'hard': []}
    for snap in snapshots:
        for cat in ['truly', 'soft', 'hard']:
            exps = snap['expectations'].get(cat, {})
            if not exps:
                continue
            y_odd_vals = [exps.get(k, float('nan')) for k in Y_ODD_KEYS]
            max_y_odd = max(abs(v) for v in y_odd_vals if not math.isnan(v))

            # Y-parity verdict: small max → preserved (consistent with prediction)
            if cat == 'soft':
                # Predicted preserved
                verdict = "✓ Y-parity preserved" if max_y_odd < NOISE_THRESHOLD else "✗ exceeds noise"
            else:
                # Predicted NOT preserved → expect some large
                verdict = "✓ Y-odd nonzero (expected)" if max_y_odd > NOISE_THRESHOLD else "  small (within noise)"

            print(f"{snap['key']:<22s}  {cat:<6s}  ", end='')
            for v in y_odd_vals:
                print(f"{v:>+7.3f}  ", end='')
            print(f"{max_y_odd:>10.4f}  {verdict:<24s}")

            summary[cat].append({
                'backend': snap['key'],
                'max_y_odd': max_y_odd,
                'y_odd_vals': y_odd_vals,
            })
        print()

    print()
    print("=" * 100)
    print("Per-category summary:")
    print("=" * 100)
    print()

    for cat in ['truly', 'soft', 'hard']:
        rows = summary[cat]
        if not rows:
            continue
        max_vals = [r['max_y_odd'] for r in rows]
        mean_max = float(np.mean(max_vals))
        std_max = float(np.std(max_vals))
        print(f"{cat:<6s}: n_runs={len(rows)}, mean max|Y-odd|={mean_max:.4f} ± {std_max:.4f}")
        n_in_band = sum(1 for v in max_vals if v < NOISE_THRESHOLD)
        n_out = len(rows) - n_in_band
        print(f"        runs with max|Y-odd| < {NOISE_THRESHOLD}: {n_in_band}/{len(rows)}")
        print(f"        runs with max|Y-odd| ≥ {NOISE_THRESHOLD}: {n_out}/{len(rows)}")
    print()

    print("Y-parity verification verdict:")
    print()
    print("  soft XY+YX (predicted Y-parity preserved):")
    soft_max_all = [r['max_y_odd'] for r in summary['soft']]
    soft_in_band = all(v < NOISE_THRESHOLD for v in soft_max_all) if soft_max_all else False
    if soft_in_band:
        print(f"    ✓ ALL {len(soft_max_all)} backends show max|Y-odd| < {NOISE_THRESHOLD}")
        print(f"    ✓ Y-parity preservation hardware-verified across 3 Heron r2 backends")
    else:
        print(f"    Some runs exceed noise — review individually")
    print()
    print("  truly + hard (predicted Y-parity NOT preserved):")
    other_max = ([r['max_y_odd'] for r in summary['truly']]
                  + [r['max_y_odd'] for r in summary['hard']])
    n_large = sum(1 for v in other_max if v > NOISE_THRESHOLD)
    print(f"    {n_large}/{len(other_max)} runs show max|Y-odd| > {NOISE_THRESHOLD}")
    print(f"    (expected: Y-parity broken → some cells significantly above noise)")


if __name__ == "__main__":
    main()
