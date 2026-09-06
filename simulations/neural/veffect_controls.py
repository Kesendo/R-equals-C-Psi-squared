#!/usr/bin/env python3
"""Executable resolution and transpose controls for the neural V census.

The transposed matrices have the same exact characteristic polynomial.  Any
bin-count change therefore reads numerical eigensolver/resolution sensitivity,
not a change in physical modes.
"""

from pathlib import Path
import sys

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from neural_palindrome import (
    build_exact_weights,
    build_linear_jacobian,
    make_balanced_dale_network,
)
from veffect_and_heat import build_jacobian_with_sigmoid, frequency_counts


TAU_E, TAU_I, ALPHA = 5.0, 10.0, 0.5
COUPLINGS = (0.0, 0.05)
RESOLUTIONS = (1e-6, 2.5e-7, 6.25e-8)
DRIVE_GRID = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0)


def coupled_jacobian(n, coupling):
    """Build the same odd-seat bridge used by ``veffect_exact.py``."""
    weights_a, signs_a, _ = build_exact_weights(
        n, n // 2, TAU_E, TAU_I, density=0.3, seed=42
    )
    weights_b, signs_b, _ = build_exact_weights(
        n, n // 2, TAU_E, TAU_I, density=0.3, seed=99
    )
    size = 2 * n + 1
    weights = np.zeros((size, size))
    signs = np.zeros(size)
    weights[:n, :n] = weights_a
    weights[n:2 * n, n:2 * n] = weights_b
    signs[:n] = signs_a
    signs[n:2 * n] = signs_b
    signs[2 * n] = 1.0
    for offset in (0, n):
        for endpoint in (offset, offset + n - 1):
            weights[2 * n, endpoint] = coupling
            weights[endpoint, 2 * n] = coupling
    return build_linear_jacobian(weights, signs, TAU_E, TAU_I, ALPHA)


def coupling_controls():
    """Return keyed (K_act, K_corr) reads without mutating module bindings."""
    rows = {}
    for n in (10, 20):
        for coupling in COUPLINGS:
            jacobian = coupled_jacobian(n, coupling)
            direct_values = np.linalg.eigvals(jacobian)
            for resolution in RESOLUTIONS:
                rows[(n, coupling, resolution, False)] = frequency_counts(
                    direct_values, resolution
                )
            transposed_values = np.linalg.eigvals(jacobian.T)
            rows[(n, coupling, RESOLUTIONS[0], True)] = frequency_counts(
                transposed_values, RESOLUTIONS[0]
            )
    return rows


def drive_solver_controls():
    """Compare the recorded drive-grid counts at solver tolerances 1e-12/1e-14."""
    weights, signs = make_balanced_dale_network(50, 25, density=0.3, seed=42)
    rows = {}
    for drive in DRIVE_GRID:
        baseline, _, baseline_residual = build_jacobian_with_sigmoid(
            weights, signs, TAU_E, TAU_I, 0.3, drive, solver_tol=1e-12
        )
        refined, _, refined_residual = build_jacobian_with_sigmoid(
            weights, signs, TAU_E, TAU_I, 0.3, drive, solver_tol=1e-14
        )
        rows[drive] = {
            "baseline_counts": frequency_counts(np.linalg.eigvals(baseline), 1e-4),
            "refined_counts": frequency_counts(np.linalg.eigvals(refined), 1e-4),
            "baseline_residual": baseline_residual,
            "refined_residual": refined_residual,
        }
    return rows


def render_report():
    rows = coupling_controls()
    drive_rows = drive_solver_controls()
    lines = [
        "NEURAL V-CENSUS NUMERICAL CONTROLS",
        "Transposition preserves the exact spectrum; count changes expose numerical sensitivity.",
    ]
    for key in sorted(rows):
        n, coupling, resolution, transposed = key
        k_act, k_corr = rows[key]
        route = "transpose" if transposed else "direct"
        lines.append(
            f"N={n} c={coupling:g} eps={resolution:g} {route}: "
            f"K_act={k_act} K_corr={k_corr}"
        )
    stable = all(
        row["baseline_counts"] == row["refined_counts"]
        for row in drive_rows.values()
    )
    if not stable:
        raise RuntimeError("drive-grid counts changed under solver refinement")
    max_residual = max(row["refined_residual"] for row in drive_rows.values())
    lines.append(
        "solver tol 1e-14: all 13 drive-grid counts stable; "
        f"max fresh equation residual={max_residual:.2e}"
    )
    return "\n".join(lines) + "\n"


def main():
    print(render_report(), end="")


if __name__ == "__main__":
    main()
