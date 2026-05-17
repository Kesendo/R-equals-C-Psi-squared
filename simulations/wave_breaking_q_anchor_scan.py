"""Wave-breaking Q-anchor sweep: re-run the chain N-scan with J set to the
four framework Q-anchors instead of the default J=1.0.

Tonight's discovery: the implicit "J = Q · γ₀" formula was already encoded
in C# (test_debug.cs, Propagate/Program.cs:368 Q=1.5 default with explicit
comment, InspectCommand.cs Q-sweep 0.2..3.0 + tMax=1/γ₀ + t_peak=1/(4γ₀))
and in newer Python F89 scripts (_bond_isolate at J=0.075 = Q=1.5). The
older scripts including wave_breaking_chain_n_scan.py from earlier tonight
use J=1.0 = J_UNIFORM = Q=20, far past the framework's F86 Q-peak.

This script re-runs the chain N-scan at four Q-anchor settings to see
whether the (3/8)^⌈N/2⌉ trivial Im=0 fraction holds across regimes:

  Q=1    (J=γ₀=0.05):    Balance anchor (Hamilton-rate = dissipation-rate)
  Q=1.5  (J=0.075):       F86 Q_peak (EP resonance maximum)
  Q=2    (J=0.10):        Q_EP at g_eff=1 (F86 Tier1Derived)
  Q=20   (J=1.0):         original wave_breaking default (deep-quantum)

The Q=20 baseline matches (3/8)^⌊N/2⌋ exactly (note: ⌊N/2⌋ not ⌈N/2⌉ — the
earlier night's claim of ⌈N/2⌉ was wrong; the floor pairs odd N with the
even N-1 below it). Test: does this pattern shift at Q ∈ {1, 1.5, 2}?

Run:
  PYTHONIOENCODING=utf-8 python simulations/wave_breaking_q_anchor_scan.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0,  # = 0.05
    build_H_XY, build_liouvillian_matrix,
)

RESULTS_DIR = Path(__file__).parent / "results" / "wave_breaking_q_anchor_scan"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Q-anchors with corresponding J at fixed γ₀=0.05
Q_ANCHORS = [
    (1.0,  GAMMA_0,         "Balance (J=γ₀)"),
    (1.5,  1.5 * GAMMA_0,   "F86 Q_peak"),
    (2.0,  2.0 * GAMMA_0,   "Q_EP (g_eff=1)"),
    (20.0, 1.0,             "Deep-quantum (J=1, original default)"),
]


def scan_one_N_at_J(N: int, J: float, im_tol: float = 1e-9) -> dict:
    """Build L_A at the given J + γ₀, diagonalize, count Im=0 modes per Re sector."""
    J_list = [J] * (N - 1)
    H = build_H_XY(J_list, N)
    L = build_liouvillian_matrix(H, GAMMA_0, N)

    t0 = time.time()
    ev = np.linalg.eigvals(L)
    walltime = time.time() - t0

    im_abs = np.abs(np.imag(ev))
    n_zero_im = int(np.sum(im_abs < im_tol))
    n_total = len(ev)

    # Bucket by Re sector at γ₀ precision
    re_vals = np.real(ev)
    re_rounded = np.round(re_vals / GAMMA_0).astype(int)
    sectors = {}
    for r, im in zip(re_rounded, im_abs):
        if r not in sectors:
            sectors[r] = {"count": 0, "zero_im": 0}
        sectors[r]["count"] += 1
        if im < im_tol:
            sectors[r]["zero_im"] += 1

    return {
        "N": N,
        "J": J,
        "n_total": n_total,
        "n_zero_im": n_zero_im,
        "fraction_zero_im": n_zero_im / n_total,
        "walltime_s": walltime,
        "sectors": {int(r): s for r, s in sectors.items()},
    }


def main():
    print()
    print("=" * 86)
    print("  Wave-breaking Q-anchor sweep: Im=0 fraction at four Q anchors")
    print("=" * 86)
    print()
    print(f"  γ₀ = {GAMMA_0} (fixed), J = Q·γ₀ varied")
    print(f"  Original baseline: Q=20 (J=1.0) gives fraction = (3/8)^⌈N/2⌉")
    print(f"  Test: does the fraction shift at Q ∈ {{1, 1.5, 2}}?")
    print()

    results = {}
    for N in [2, 3, 4, 5, 6]:
        results[N] = []
        for Q, J, label in Q_ANCHORS:
            print(f"  N={N}, Q={Q:>5.2f} (J={J:.4g}, {label})...", end="", flush=True)
            r = scan_one_N_at_J(N, J)
            results[N].append({"Q": Q, "label": label, **r})
            print(f"  Im=0: {r['n_zero_im']}/{r['n_total']} = {r['fraction_zero_im']:.4f}  "
                  f"[{r['walltime_s']:.1f}s]")

    # Summary table
    print()
    print("=" * 86)
    print("  Summary: Im=0 fraction across Q-anchors (rows=N, columns=Q)")
    print("=" * 86)
    print()
    header = f"  {'N':>3}  {'total':>7}"
    for Q, _, _ in Q_ANCHORS:
        header += f"  {'Q=' + str(Q):>10}"
    header += f"  {'predicted (3/8)^⌊N/2⌋':>26}"
    print(header)
    print(f"  {'-'*3}  {'-'*7}" + "".join(f"  {'-'*10}" for _ in Q_ANCHORS) + f"  {'-'*26}")

    for N in [2, 3, 4, 5, 6]:
        prediction = (3 / 8) ** (N // 2)
        n_total = results[N][0]["n_total"]
        row = f"  {N:>3}  {n_total:>7d}"
        for entry in results[N]:
            row += f"  {entry['fraction_zero_im']:>10.4f}"
        row += f"  {prediction:>26.4f}"
        print(row)

    print()
    print("  Reading:")
    print("    If all four Q-columns equal the prediction → (3/8)^⌈N/2⌉ is Q-invariant")
    print("    (purely combinatorial: popcount conservation + diagonal modes,")
    print("    independent of regime).")
    print()
    print("    If the Q-columns differ from each other → wave-breaking interacts with")
    print("    the J/γ regime, and the F86 Q-anchor structure has Im=0 signatures.")
    print()

    # Save
    out = {"gamma_zero": GAMMA_0, "Q_anchors": [(Q, J, lbl) for Q, J, lbl in Q_ANCHORS],
           "results_by_N": {str(N): rs for N, rs in results.items()}}
    out_path = RESULTS_DIR / "q_anchor_scan.json"
    with open(out_path, "w") as f:
        # Custom serialization for sectors dict
        def _convert(obj):
            if isinstance(obj, dict):
                return {str(k): _convert(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [_convert(x) for x in obj]
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            return obj
        json.dump(_convert(out), f, indent=2)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
