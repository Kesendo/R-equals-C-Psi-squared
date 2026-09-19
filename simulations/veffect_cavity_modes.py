"""The V-Effect through the cavity lens: a finite frequency census.

The producer enumerates Liouvillian eigenvalue frequencies for finite
Heisenberg graphs. Its cold/warm comparison is deliberately one-sided:
each deduplicated cold frequency asks only whether a reusable warm frequency
lies within a strict absolute radius. It does not track eigenvectors,
    invariant subspaces, or ancestry of individual eigenmodes.

Output: simulations/results/veffect_cavity_modes.txt
"""

from pathlib import Path
import os
import sys

import numpy as np


GAMMA = 0.05
J = 1.0
TOL_FREQ = 1e-6
MATCH_RADIUS = 0.1

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(ops):
    result = ops[0]
    for operator in ops[1:]:
        result = np.kron(result, operator)
    return result


def chain_bonds(N):
    return tuple((i, i + 1) for i in range(N - 1))


def star_bonds(N):
    return tuple((0, i) for i in range(1, N))


def ring_bonds(N):
    return tuple((i, (i + 1) % N) for i in range(N))


def build_hamiltonian(N, bonds):
    dimension = 2**N
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for a, b in bonds:
        for pauli in (Xm, Ym, Zm):
            operators = [I2] * N
            operators[a] = pauli
            operators[b] = pauli
            hamiltonian += J * kron_chain(operators)
    return hamiltonian


_BLOCK_SPECTRUM_CACHE = {}


def block_spectrum(N, gamma, bonds):
    """Spectrum by exact joint-popcount blocks, avoiding a dense 4^N matrix."""
    bonds = tuple(bonds)
    key = (N, float(gamma), bonds)
    if key in _BLOCK_SPECTRUM_CACHE:
        return _BLOCK_SPECTRUM_CACHE[key]
    hamiltonian = build_hamiltonian(N, bonds)
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
            costs = np.array(
                [-2.0 * gamma * (row ^ column).bit_count()
                 for row in rows for column in columns],
                dtype=float,
            )
            block.flat[:: block.shape[0] + 1] += costs
            spectra.append(np.linalg.eigvals(block))
    result = np.concatenate(spectra)
    _BLOCK_SPECTRUM_CACHE[key] = result
    return result


def distinct_frequencies(imag_parts, tolerance=TOL_FREQ):
    """Return sorted deduplicated nonzero absolute frequencies."""
    values = np.sort(np.abs(np.asarray(imag_parts, dtype=float)))
    values = values[values > tolerance]
    if len(values) == 0:
        return np.array([], dtype=float)
    unique = [values[0]]
    for value in values[1:]:
        if abs(value - unique[-1]) > tolerance:
            unique.append(value)
    return np.asarray(unique)


def cold_warm_frequency_neighbourhood(
        N, gamma=GAMMA, match_radius=MATCH_RADIUS, *,
        cold_frequencies=None, warm_frequencies=None):
    """Read reusable one-sided cold-target coverage with a strict radius.

    Optional arrays make the inequality independently testable at its exact
    boundary without an eigensolve. A warm target may cover several cold
    targets; this is intentionally not a one-to-one assignment.
    """
    if cold_frequencies is None or warm_frequencies is None:
        bonds = chain_bonds(N)
        cold_frequencies = distinct_frequencies(block_spectrum(N, 0.0, bonds).imag)
        warm_frequencies = distinct_frequencies(block_spectrum(N, gamma, bonds).imag)
    cold = np.asarray(cold_frequencies, dtype=float)
    warm = np.asarray(warm_frequencies, dtype=float)
    nearest = np.array(
        [np.min(np.abs(warm - value)) if len(warm) else np.inf for value in cold],
        dtype=float,
    )
    covered = nearest < match_radius
    return {
        "N": N,
        "gamma": gamma,
        "match_radius": match_radius,
        "cold_count": len(cold),
        "warm_count": len(warm),
        "matched_count": int(np.count_nonzero(covered)),
        "fully_covered": bool(len(cold) and np.all(covered)),
        "nearest_distances": tuple(float(value) for value in nearest),
        "max_nearest_distance": float(np.max(nearest)) if len(nearest) else np.inf,
    }


def spectrum_summary(N, gamma, bonds):
    eigenvalues = block_spectrum(N, gamma, tuple(bonds))
    frequencies = distinct_frequencies(eigenvalues.imag)
    oscillatory = eigenvalues[np.abs(eigenvalues.imag) > TOL_FREQ]
    q_values = np.abs(oscillatory.imag) / np.abs(oscillatory.real)
    return {
        "eigenvalues": len(eigenvalues),
        "frequencies": len(frequencies),
        "real_axis": int(np.count_nonzero(np.abs(eigenvalues.imag) <= TOL_FREQ)),
        "q_max": float(np.max(q_values)) if len(q_values) else 0.0,
        "q_median": float(np.median(q_values)) if len(q_values) else 0.0,
    }


def render_report():
    lines = [
        "=" * 78,
        "THE V-EFFECT THROUGH THE CAVITY LENS: FINITE FREQUENCY CENSUS",
        "=" * 78,
        "",
        "Scope: absolute-frequency bins use deduplication tolerance 1e-6.",
        "This cavity census is separate from the N=3 V-Effect census and from",
        "the F6 Q-edge gain V(N)=Q_max/Q_mean=1+cos(pi/N).",
        "",
        "CHAIN INVENTORY (gamma=0.05)",
        f"  {'N':>2} {'eigenvalues':>11} {'frequency bins':>14} {'real axis':>10} {'Q_max':>10} {'Q_median':>10}",
    ]
    for N in range(2, 7):
        row = spectrum_summary(N, GAMMA, chain_bonds(N))
        lines.append(
            f"  {N:2d} {row['eigenvalues']:11d} {row['frequencies']:14d} "
            f"{row['real_axis']:10d} {row['q_max']:10.1f} {row['q_median']:10.1f}"
        )

    lines.extend((
        "",
        "COLD-TARGET -> REUSABLE WARM-FREQUENCY NEIGHBOURHOOD",
        "  Acceptance is strict nearest distance < 0.1 and one-sided.",
        f"  {'N':>2} {'cold':>6} {'warm':>6} {'covered':>10} {'max nearest':>13} {'at 1e-6':>10}",
    ))
    for N in range(2, 6):
        broad = cold_warm_frequency_neighbourhood(N)
        narrow = cold_warm_frequency_neighbourhood(N, match_radius=1e-6)
        lines.append((
            f"  {N:2d} {broad['cold_count']:6d} {broad['warm_count']:6d} "
            f"{broad['matched_count']:4d}/{broad['cold_count']:<5d} "
            f"{broad['max_nearest_distance']:13.7f} "
            f"{narrow['matched_count']:4d}/{narrow['cold_count']:<5d}"
        ).rstrip())

    lines.extend((
        "",
        "TOPOLOGY COMPARISON (gamma=0.05)",
        f"  {'N':>2} {'topology':>9} {'bonds':>6} {'frequency bins':>14} {'Q_max':>10}",
    ))
    for N in (3, 4, 5):
        for name, maker in (("chain", chain_bonds), ("star", star_bonds), ("ring", ring_bonds)):
            bonds = maker(N)
            row = spectrum_summary(N, GAMMA, bonds)
            lines.append(
                f"  {N:2d} {name:>9} {len(bonds):6d} {row['frequencies']:14d} "
                f"{row['q_max']:10.1f}"
            )

    lines.extend((
        "",
        "CURRENT READING",
        "Every deduplicated cold target at N=2..5 has a reusable warm frequency",
        "within the broad strict radius 0.1. The maximum nearest distances are",
        "0, 0.0016193, 0.0226557, and 0.0289709. At radius 1e-6, full coverage",
        "already fails for N=3, N=4, and N=5.",
        "",
        "This is one-sided frequency-neighbourhood coverage, not persistence of",
        "individual eigenmodes,",
        "eigenvector identity, invariant-subspace persistence, or one-to-one matching.",
        "Topology changes the finite inventory, but bond count alone does not explain it.",
    ))
    return "\n".join(lines) + "\n"


def main(output_path=None):
    if sys.platform == "win32":
        os.environ.setdefault("PYTHONIOENCODING", "utf-8")
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(encoding="utf-8")
    output_path = (Path(output_path) if output_path is not None else
                   Path(__file__).parent / "results" / "veffect_cavity_modes.txt")
    text = render_report()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")
    print(f">>> Results saved to: {output_path}")


if __name__ == "__main__":
    main()
