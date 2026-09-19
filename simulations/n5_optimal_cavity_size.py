"""N=5 and the Goldilocks question: a bounded finite cavity census.

The dense inventory and sacrifice scan cover exactly N=(3,4,5).  N=5 is the
largest sampled row, so a metric winner here is not a global optimum.  No row
outside that domain is used to select a winner.

Output: simulations/results/n5_optimal_cavity_size.txt
"""

from fractions import Fraction
from pathlib import Path
import os
import sys

import numpy as np

try:
    from simulations.veffect_cavity_modes import (
        J, TOL_FREQ, build_hamiltonian, chain_bonds, distinct_frequencies,
    )
except ModuleNotFoundError as exc:  # direct ``python simulations/...py`` execution
    if exc.name != "simulations":
        raise
    from veffect_cavity_modes import (
        J, TOL_FREQ, build_hamiltonian, chain_bonds, distinct_frequencies,
    )


PHI = (1 + np.sqrt(5)) / 2
EPS = 0.001
DENSE_CENSUS_N = (3, 4, 5)
SACRIFICE_SCAN_N = (3, 4, 5)
SACRIFICE_GAMMAS = (0.01, 0.02, 0.05, 0.1, 0.2)


def f6_q_edge_gain(N):
    """F6 within-N Q_max/Q_mean for the declared single-excitation model."""
    return 1 + np.cos(np.pi / N)


def cold_frequencies(N):
    m = np.arange(1, N, dtype=float)
    return 4 * J * (1 - np.cos(np.pi * m / N))


def cold_frequency_ratios(N):
    frequencies = cold_frequencies(N)
    return frequencies[1:] / frequencies[:-1]


def fixed_index_ratio_limit(numerator_index, denominator_index):
    """lim N->infinity omega_k/omega_j = k^2/j^2 at fixed indices."""
    return Fraction(numerator_index**2, denominator_index**2)


_PROFILE_BLOCK_SPECTRUM_CACHE = {}


def profile_block_spectrum(N, gammas):
    """Joint-popcount block spectrum for a site-dependent dephasing profile."""
    gammas = tuple(float(value) for value in gammas)
    key = (N, gammas)
    if key in _PROFILE_BLOCK_SPECTRUM_CACHE:
        return _PROFILE_BLOCK_SPECTRUM_CACHE[key]
    hamiltonian = build_hamiltonian(N, chain_bonds(N))
    basis = tuple(
        tuple(i for i in range(2**N) if i.bit_count() == population)
        for population in range(N + 1)
    )
    spectra = []
    for rows in basis:
        hp = hamiltonian[np.ix_(rows, rows)]
        for columns in basis:
            hq = hamiltonian[np.ix_(columns, columns)]
            rp, cq = len(rows), len(columns)
            block = -1j * (
                np.kron(hp, np.eye(cq)) - np.kron(np.eye(rp), hq.T)
            )
            costs = np.array([
                -2.0 * sum(
                    gamma for site, gamma in enumerate(gammas)
                    if ((row ^ column) >> site) & 1
                )
                for row in rows for column in columns
            ])
            block.flat[:: block.shape[0] + 1] += costs
            spectra.append(np.linalg.eigvals(block))
    result = np.concatenate(spectra)
    _PROFILE_BLOCK_SPECTRUM_CACHE[key] = result
    return result


def spectrum_metrics(N, gammas):
    eigenvalues = profile_block_spectrum(N, tuple(gammas))
    frequencies = distinct_frequencies(eigenvalues.imag)
    oscillatory = eigenvalues[np.abs(eigenvalues.imag) > TOL_FREQ]
    q_values = np.abs(oscillatory.imag) / np.abs(oscillatory.real)
    gaps = np.diff(frequencies)
    return {
        "frequency_bins": len(frequencies),
        "q_max": float(np.max(q_values)),
        "q_times_bins": float(np.max(q_values)) * len(frequencies),
        "q_times_sqrt_bins": float(np.max(q_values)) * np.sqrt(len(frequencies)),
        "mean_gap": float(np.mean(gaps)) if len(gaps) else 0.0,
        "min_gap": float(np.min(gaps)) if len(gaps) else 0.0,
        "min_mean_ratio": float(np.min(gaps) / np.mean(gaps)) if len(gaps) else 0.0,
    }


def _winner(rows, metric, *, maximize=True):
    chooser = max if maximize else min
    value = chooser(row[metric] for row in rows.values())
    return tuple(N for N, row in rows.items() if row[metric] == value)


def render_report():
    lines = [
        "=" * 78,
        "N=5 AND THE GOLDILOCKS QUESTION: A FINITE CAVITY CENSUS",
        "=" * 78,
        "",
        "Dense mode/Q census and sacrifice domain: sampled N=(3,4,5).",
        "N=5 is the largest tested row; this scan does not select N=5 as an optimum.",
        "",
        "EXACT PENTAGON IDENTITIES",
        f"  phi = (1+sqrt(5))/2 = {PHI:.10f}",
        f"  F6 Q-edge gain V(5) = 1+phi/2 = (5+sqrt(5))/4 = {f6_q_edge_gain(5):.10f}",
        f"  omega_2/omega_1 at N=5 = 2 + phi = {cold_frequency_ratios(5)[0]:.10f}",
        "  Both contain phi; they are different ratios.",
        "",
        "FIXED-INDEX LIMITS (indices held fixed)",
        "  omega_4/omega_3 -> 16/9",
        "  omega_3/omega_2 -> 9/4",
        "",
        "FINITE FREQUENCY RATIOS",
    ]
    for N in DENSE_CENSUS_N:
        ratios = cold_frequency_ratios(N)
        lines.append(
            f"  N={N}: " + ", ".join(f"{value:.6f}" for value in ratios)
        )
    lines.extend((
        "",
        "DENSE CENSUS (gamma=0.05)",
        f"  {'N':>2} {'bins':>6} {'Q_max':>10} {'Q*bins':>12} {'Q*sqrt(bins)':>14} {'min/mean gap':>13}",
    ))
    dense = {
        N: spectrum_metrics(N, (0.05,) * N) for N in DENSE_CENSUS_N
    }
    for N, row in dense.items():
        lines.append(
            f"  {N:2d} {row['frequency_bins']:6d} {row['q_max']:10.3f} "
            f"{row['q_times_bins']:12.3f} {row['q_times_sqrt_bins']:14.3f} "
            f"{row['min_mean_ratio']:13.6f}"
        )
    lines.extend((
        "",
        "FINITE METRIC WINNERS WITHIN N=(3,4,5)",
    ))
    for label, key, maximize in (
        ("frequency bins", "frequency_bins", True),
        ("Q_max", "q_max", True),
        ("Q*bins", "q_times_bins", True),
        ("Q*sqrt(bins)", "q_times_sqrt_bins", True),
        ("mean gap", "mean_gap", True),
        ("minimum gap", "min_gap", True),
        ("min/mean gap ratio", "min_mean_ratio", True),
    ):
        lines.append(f"  {label}: N={','.join(map(str, _winner(dense, key, maximize=maximize)))}")

    lines.extend((
        "",
        "SACRIFICE-PROFILE Q_max IMPROVEMENT (sampled N=(3,4,5))",
        f"  {'gamma':>8} {'N=3':>10} {'N=4':>10} {'N=5':>10} {'largest sampled':>17}",
    ))
    for gamma in SACRIFICE_GAMMAS:
        improvements = {}
        for N in SACRIFICE_SCAN_N:
            uniform = spectrum_metrics(N, (gamma,) * N)["q_max"]
            edge = N * gamma - (N - 1) * EPS
            sacrifice = spectrum_metrics(N, (edge,) + (EPS,) * (N - 1))["q_max"]
            improvements[N] = sacrifice / uniform
        winner = max(improvements, key=improvements.get)
        lines.append(
            f"  {gamma:8.3f} {improvements[3]:10.3f} {improvements[4]:10.3f} "
            f"{improvements[5]:10.3f} {('N='+str(winner)):>17}"
        )
    lines.extend((
        "",
        "BOUNDED VERDICT",
        "The exact phi identities survive. Every finite metric winner is reported,",
        "but the producer samples only N=(3,4,5), so it cannot locate a global N-optimum or",
        "establish a universal Goldilocks mechanism. The N-optimum question remains OPEN.",
    ))
    return "\n".join(lines) + "\n"


def main(output_path=None):
    if sys.platform == "win32":
        os.environ.setdefault("PYTHONIOENCODING", "utf-8")
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(encoding="utf-8")
    output_path = (Path(output_path) if output_path is not None else
                   Path(__file__).parent / "results" / "n5_optimal_cavity_size.txt")
    text = render_report()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")
    print(f">>> Results saved to: {output_path}")


if __name__ == "__main__":
    main()
