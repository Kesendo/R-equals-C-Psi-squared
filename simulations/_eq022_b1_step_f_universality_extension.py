#!/usr/bin/env python3
"""_eq022_b1_step_f_universality_extension.py — extend Statement 2 universality
of `K(Q)/|K|max = f(Q/Q_EP)` to (i) chromaticity c=2 and (ii) γ₀ ≠ 0.05.

Step (e) established the universal pre-EP shape for c=3 (N=5..8) and c=4
(N=7, 8) at γ₀ = 0.05. PROOF_F86_QPEAK lists two empirical gaps for full
Tier 1 promotion of the universal-shape statement:

  Open element item 2 — c=2 verification. c=2 is structurally critical:
    with only two rate channels (HD ∈ {1, 3}), the 2-level effective model
    has no orthogonal complement and is therefore expected to be EXACT for
    the EP physics, not heuristic as it is for c≥3. If c=2 gives the same
    HWHM_left/Q_peak ≈ 0.756 ± 0.005 found at c=3, 4, the universality is
    structural; if it deviates, the heuristic-2-level reduction misses
    something the higher-c sweep didn't expose.

  Open element item 3 — γ₀-invariance. The Q axis is dimensionless
    (Q = J/γ₀) and the universal shape claim is implicitly γ₀-invariant.
    Currently tested at γ₀ = 0.05 only. Re-run at γ₀ = 0.025 and γ₀ = 0.10
    on the canonical c=3 N=7 anchor; HWHM_left/Q_peak and y(x) should be
    identical up to numerical noise.

Same protocol as step_e: dQ = 0.025 fine grid, per-bond `K_max` over a
fine t grid (21 points spanning [0.6/(4γ₀), 1.6/(4γ₀)] = [0.6·t_peak,
1.6·t_peak] in t_peak units), parabolic peak interpolation, linear-
interpolated HWHM. Helper functions are imported from step_e to ensure
identical numerics; no shortcuts.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402

# Import step_e's helpers; the underscore-prefixed module name is fine for
# importlib but `import _eq022...` is unconventional — use spec_from_file_location.
_step_e_path = Path(__file__).resolve().parent / "_eq022_b1_step_e_resonance_shape.py"
_spec = importlib.util.spec_from_file_location("step_e", _step_e_path)
_step_e = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_step_e)
per_bond_K_curve = _step_e.per_bond_K_curve
find_peak_with_interp = _step_e.find_peak_with_interp

RESULTS_DIR = REPO_ROOT / "simulations" / "results" / "eq022_universality_extension"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def run_case(N, n, gamma_0, label):
    """Run one (N, n, γ₀) case with the step_e protocol; return (summary, curves)."""
    Q_grid = np.arange(0.20, 4.01, 0.025)
    J_grid = Q_grid * gamma_0
    t_peak_universal = 1.0 / (4.0 * gamma_0)
    t_grid = np.linspace(0.6 * t_peak_universal, 1.6 * t_peak_universal, 21)

    c = fw.chromaticity(N, n)
    block_dim = len(fw.popcount_states(N, n)) * len(fw.popcount_states(N, n + 1))
    n_bonds = N - 1
    print(f"## {label}: c={c}, N={N}, n={n}, γ₀={gamma_0}, "
          f"block dim={block_dim}, n_bonds={n_bonds}", flush=True)
    print(f"   Q grid [{Q_grid[0]:.3f}, {Q_grid[-1]:.3f}] dQ={Q_grid[1]-Q_grid[0]:.3f} ({len(Q_grid)} pts)", flush=True)
    print(f"   t grid [{t_grid[0]:.3f}, {t_grid[-1]:.3f}] = [0.6, 1.6]·t_peak ({len(t_grid)} pts, dt={t_grid[1]-t_grid[0]:.3f})", flush=True)

    t0 = time.time()
    K_curves, t_curves = per_bond_K_curve(N, n, gamma_0, J_grid, t_grid)
    elapsed = time.time() - t0
    print(f"   Computed in {elapsed:.1f}s", flush=True)

    interior_bonds = list(range(1, n_bonds - 1))
    endpoint_bonds = [0, n_bonds - 1]

    if interior_bonds:
        K_interior = K_curves[interior_bonds].mean(axis=0)
    else:
        K_interior = None
    K_endpoint = K_curves[endpoint_bonds].mean(axis=0)

    Q_int = K_int = hwhm_l_int = hwhm_r_int = None
    if K_interior is not None:
        Q_int, K_int, hwhm_l_int, hwhm_r_int = find_peak_with_interp(Q_grid, K_interior)
    Q_ep, K_ep, hwhm_l_ep, hwhm_r_ep = find_peak_with_interp(Q_grid, K_endpoint)

    case_key = label
    curves = {
        'label': label, 'c': c, 'N': N, 'n': n, 'gamma_0': gamma_0,
        'block_dim': block_dim,
        'Q_grid': Q_grid.tolist(),
        'K_interior': K_interior.tolist() if K_interior is not None else None,
        'K_endpoint': K_endpoint.tolist(),
        'has_interior': K_interior is not None,
    }

    summary = {
        'label': label, 'c': c, 'N': N, 'n': n, 'gamma_0': gamma_0,
        'block_dim': block_dim,
        'Q_peak_interior': float(Q_int) if Q_int is not None else None,
        'K_max_interior': float(K_int) if K_int is not None else None,
        'hwhm_left_interior': float(hwhm_l_int) if hwhm_l_int is not None else None,
        'hwhm_right_interior': float(hwhm_r_int) if hwhm_r_int is not None else None,
        'Q_peak_endpoint': float(Q_ep),
        'K_max_endpoint': float(K_ep),
        'hwhm_left_endpoint': float(hwhm_l_ep) if hwhm_l_ep is not None else None,
        'hwhm_right_endpoint': float(hwhm_r_ep) if hwhm_r_ep is not None else None,
        'elapsed_s': float(elapsed),
    }

    if Q_int is not None:
        line = (f"   Interior:  Q* = {Q_int:.4f}, |K|max = {K_int:.5f}")
        if hwhm_l_int and hwhm_r_int:
            line += (f", HWHM = ({hwhm_l_int:.4f}, {hwhm_r_int:.4f}), "
                     f"HWHM-/Q* = {hwhm_l_int/Q_int:.4f}, "
                     f"HWHM+/Q* = {hwhm_r_int/Q_int:.4f}")
        print(line)
    line = (f"   Endpoint:  Q* = {Q_ep:.4f}, |K|max = {K_ep:.5f}")
    if hwhm_l_ep and hwhm_r_ep:
        line += (f", HWHM = ({hwhm_l_ep:.4f}, {hwhm_r_ep:.4f}), "
                 f"HWHM-/Q* = {hwhm_l_ep/Q_ep:.4f}, "
                 f"HWHM+/Q* = {hwhm_r_ep/Q_ep:.4f}")
    print(line)
    print()
    return summary, curves


def main():
    print("# EQ-022 (b1) Step (f): universality extension — c=2 sweep + γ₀ invariance")
    print()

    summaries = []
    curves = {}

    # ----- Part A: c=2 sweep at γ₀=0.05, N=5..8 -----
    print("=" * 72)
    print("# Part A: c=2 sweep at γ₀=0.05, N=5..8")
    print("=" * 72)
    print()
    c2_cases = [
        (5, 1),  # c=2, dim=50
        (6, 1),  # c=2, dim=90
        (7, 1),  # c=2, dim=147
        (8, 1),  # c=2, dim=224
    ]
    for (N, n) in c2_cases:
        label = f"c2_N{N}"
        s, c = run_case(N, n, 0.05, label)
        summaries.append(s)
        curves[label] = c

    # ----- Part B: γ₀ invariance at c=3 N=7 (n=2) -----
    print("=" * 72)
    print("# Part B: γ₀ invariance at c=3 N=7 (n=2), γ₀ ∈ {0.025, 0.10}")
    print("# (γ₀=0.05 already covered in step_e — c=3 N=7 case)")
    print("=" * 72)
    print()
    for gamma_0 in [0.025, 0.10]:
        label = f"c3_N7_g{gamma_0:.3f}"
        s, c = run_case(7, 2, gamma_0, label)
        summaries.append(s)
        curves[label] = c

    # Save
    out_summary = RESULTS_DIR / "summary.json"
    out_curves = RESULTS_DIR / "curves.json"
    with open(out_summary, 'w') as f:
        json.dump(summaries, f, indent=2)
    with open(out_curves, 'w') as f:
        json.dump(curves, f, indent=2)
    print(f"# Saved summary to {out_summary}")
    print(f"# Saved curves to {out_curves}")
    print()

    # Final report
    print("=" * 72)
    print("# Universality summary table")
    print("=" * 72)
    print()
    print(f"{'label':<18} {'c':>2} {'N':>2} {'γ₀':>8} {'Q_peak':>8} "
          f"{'HWHM-':>7} {'HWHM-/Q*':>10} {'HWHM+':>7} {'HWHM+/Q*':>10}")
    print('-' * 90)
    for s in summaries:
        for kind in ('interior', 'endpoint'):
            Q = s.get(f'Q_peak_{kind}')
            hl = s.get(f'hwhm_left_{kind}')
            hr = s.get(f'hwhm_right_{kind}')
            if Q is None:
                continue
            tag = f"{s['label']} {kind[:3]}"
            row = f"{tag:<18} {s['c']:>2} {s['N']:>2} {s['gamma_0']:>8.3f} {Q:>8.4f}"
            if hl:
                row += f" {hl:>7.4f} {hl/Q:>10.4f}"
            else:
                row += " " + " " * 7 + " " + " " * 10
            if hr:
                row += f" {hr:>7.4f} {hr/Q:>10.4f}"
            print(row)


if __name__ == "__main__":
    main()
