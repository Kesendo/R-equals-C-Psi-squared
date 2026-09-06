#!/usr/bin/env python3
"""Multiplicity-aware census of two different spectral involutions.

The linear F1 map is ``lambda -> -lambda - 2*Sigma``. Combining F1 with
independent conjugation closure gives ``lambda -> -conj(lambda) - 2*Sigma``.
Their fixed sets differ: the former fixes only ``lambda=-Sigma``; the latter
fixes the full centre line ``Re(lambda)=-Sigma``.

Output: simulations/results/factor_two_standing_waves.txt

The N=2..7 census consumes committed CSV spectra.  Those artifacts contain
only ``Re``/``Im`` columns, so generator parameters and the regeneration
command are reported from the current C# source rather than inferred as
embedded CSV metadata.
"""

from dataclasses import dataclass
from pathlib import Path
import os
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial import cKDTree


if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

RESULTS_DIR = Path(__file__).parent / "results"
J = 1.0
GAMMA = 0.05
TOL = 1e-8
TOLERANCES = (1e-6, 1e-8, 1e-10)

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)


@dataclass(frozen=True)
class OrbitCensus:
    pairs: int
    fixed: int
    assigned: int
    max_assignment_residual: float


def partner(values, sigma, kind):
    if kind == "linear":
        return -values - 2 * sigma
    if kind == "composite":
        return -values.conjugate() - 2 * sigma
    raise ValueError(f"unknown orbit kind: {kind}")


def orbit_census(values, sigma, kind, tol=TOL):
    """Assign every spectral occurrence to one partner occurrence.

    A sparse bipartite maximum matching is used, so repeated eigenvalues are
    consumed with their multiplicity rather than merely found by membership.
    """
    values = np.asarray(values, dtype=complex)
    targets = partner(values, sigma, kind)
    points = np.column_stack((values.real, values.imag))
    target_points = np.column_stack((targets.real, targets.imag))
    neighbours = cKDTree(points).query_ball_point(target_points, tol)
    rows = []
    cols = []
    for row, candidates in enumerate(neighbours):
        rows.extend([row] * len(candidates))
        cols.extend(candidates)
    graph = csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)),
        shape=(len(values), len(values)),
    )
    assignment = maximum_bipartite_matching(graph, perm_type="column")
    assigned = int(np.count_nonzero(assignment >= 0))
    residual = (
        float(np.max(np.abs(targets - values[assignment])))
        if assigned == len(values)
        else float("inf")
    )
    fixed = int(np.count_nonzero(np.abs(targets - values) <= tol))
    remainder = len(values) - fixed
    if assigned != len(values) or remainder % 2:
        raise RuntimeError(
            f"{kind} map failed multiplicity assignment: assigned={assigned}/{len(values)}, "
            f"fixed={fixed}, residual={residual:.3e}"
        )
    return OrbitCensus(remainder // 2, fixed, assigned, residual)


def kron_chain(ops):
    result = ops[0]
    for operator in ops[1:]:
        result = np.kron(result, operator)
    return result


def build_liouvillian(n, gammas, bonds):
    dimension = 2**n
    identity = np.eye(dimension, dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for a, b in bonds:
        for pauli in (Xm, Ym, Zm):
            ops = [I2] * n
            ops[a] = pauli
            ops[b] = pauli
            hamiltonian += J * kron_chain(ops)
    generator = -1j * (
        np.kron(hamiltonian, identity) - np.kron(identity, hamiltonian.T)
    )
    for site in range(n):
        ops = [I2] * n
        ops[site] = Zm
        jump = np.sqrt(gammas[site]) * kron_chain(ops)
        jump_norm = jump.conj().T @ jump
        generator += np.kron(jump, jump.conj()) - 0.5 * (
            np.kron(jump_norm, identity) + np.kron(identity, jump_norm.T)
        )
    return generator


def chain_bonds(n):
    return [(site, site + 1) for site in range(n - 1)]


def star_bonds(n):
    return [(0, site) for site in range(1, n)]


def ring_bonds(n):
    return [(site, (site + 1) % n) for site in range(n)]


def load_eigenvalues(n):
    path = RESULTS_DIR / f"rmt_eigenvalues_N{n}.csv"
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    return data[:, 0] + 1j * data[:, 1]


def main():
    lines = []

    def log(message=""):
        print(message)
        lines.append(message)

    log("SPECTRAL ORBIT CENSUS: LINEAR F1 AND CONJUGATE-COMPOSITE MAPS")
    log("=" * 72)
    log("linear F1:            lambda -> -lambda - 2*Sigma")
    log("conjugate-composite:  lambda -> -conj(lambda) - 2*Sigma")
    log("All assignments consume spectral occurrences with multiplicity.")
    log(f"Primary numerical tolerance TOL={TOL:.0e} (not exact arithmetic).")
    log("Committed CSV inputs: simulations/results/rmt_eigenvalues_N{2..7}.csv")
    log("CSV metadata: columns Re, Im only.")
    log("Current regeneration route: dotnet run -c Release --project compute/RCPsiSquared.Compute -- rmt chain")
    log("Current generator parameters in current source: chain, J=1.0, uniform gamma=0.05.")
    log("Provenance limit: the CSV artifacts do not encode backend, source revision, command, or timestamp.")
    log()
    linear_total_pairs = linear_total_fixed = 0
    composite_total_pairs = composite_total_fixed = 0
    for n in range(2, 8):
        values = load_eigenvalues(n)
        sigma = n * GAMMA
        linear = orbit_census(values, sigma, "linear")
        composite = orbit_census(values, sigma, "composite")
        linear_total_pairs += linear.pairs
        linear_total_fixed += linear.fixed
        composite_total_pairs += composite.pairs
        composite_total_fixed += composite.fixed
        log(
            f"N={n}: linear pairs={linear.pairs} fixed={linear.fixed}; "
            f"composite pairs={composite.pairs} fixed={composite.fixed}; "
            f"max residuals={linear.max_assignment_residual:.2e},"
            f"{composite.max_assignment_residual:.2e}"
        )
        log(f"     mean decay={-float(np.mean(values.real)):.6f}; Sigma={sigma:.6f}")
    log()
    log(f"linear F1 total:       pairs={linear_total_pairs} fixed={linear_total_fixed}")
    log(
        f"conjugate-composite:   pairs={composite_total_pairs} "
        f"fixed={composite_total_fixed}"
    )
    log()
    log("TOLERANCE STABILITY (numerical multiplicities)")
    reference = None
    for tol in TOLERANCES:
        row = []
        for n in range(2, 8):
            values = load_eigenvalues(n)
            sigma = n * GAMMA
            linear = orbit_census(values, sigma, "linear", tol=tol)
            composite = orbit_census(values, sigma, "composite", tol=tol)
            row.append((linear.pairs, linear.fixed, composite.pairs, composite.fixed))
        stable = reference is None or row == reference
        if reference is None:
            reference = row
        log(f"tol={tol:.0e}: {'BASELINE' if tol == TOLERANCES[0] else ('PASS stable' if stable else 'FAIL changed')}")
        if not stable:
            raise RuntimeError(f"orbit census changes at tol={tol:.0e}")
    log()
    log("N=4 topology negative control (the fixed counts are topology-sensitive):")
    n = 4
    sigma = n * GAMMA
    for name, bond_builder in (("chain", chain_bonds), ("star", star_bonds), ("ring", ring_bonds)):
        values = np.linalg.eigvals(build_liouvillian(n, [GAMMA] * n, bond_builder(n)))
        linear = orbit_census(values, sigma, "linear")
        composite = orbit_census(values, sigma, "composite")
        log(f"N=4 {name:<9} linear: pairs={linear.pairs} fixed={linear.fixed}")
        log(f"{'':13}{'composite:':<11} pairs={composite.pairs} fixed={composite.fixed}")
    log()
    log("Scope: this is a spectral orbit census. It assigns no mode mechanism")
    log("and no transport, damping, or biological interpretation.")

    output_path = RESULTS_DIR / "factor_two_standing_waves.txt"
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
