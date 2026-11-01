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
import itertools
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
CLUSTER = 1e-4

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


def build_liouvillian(n, gammas, bonds, couplings=None):
    dimension = 2**n
    identity = np.eye(dimension, dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for index, (a, b) in enumerate(bonds):
        coupling = J if couplings is None else couplings[index]
        for pauli in (Xm, Ym, Zm):
            ops = [I2] * n
            ops[a] = pauli
            ops[b] = pauli
            hamiltonian += coupling * kron_chain(ops)
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


def render_report():
    lines = []

    def log(message=""):
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
    log("CONDITION-1 CANDIDATES: the centre line minus the F1 point")
    log("Eigenvalues on Re(lambda) = -Sigma with nonzero frequency are the composite-")
    log("fixed set minus the linear-fixed set. They are candidates for the flat-envelope")
    log("condition of docs/STANDING_WAVE_THEORY.md, not standing waves.")
    expected = {2: 6, 3: 0, 4: 138, 5: 0, 6: 1820, 7: 0}
    for n in range(2, 8):
        values = load_eigenvalues(n)
        sigma = n * GAMMA
        linear = orbit_census(values, sigma, "linear")
        composite = orbit_census(values, sigma, "composite")
        candidates = composite.fixed - linear.fixed
        distance = np.abs(values.real + sigma)
        on_line = distance <= TOL
        nearest_off = float(distance[~on_line].min())
        if on_line.any():
            farthest_on = float(distance[on_line].max())
            # An eigensolver count. What decides it is the tolerance stability of
            # both censuses (gated above at 1e-6, 1e-8, 1e-10); the spread and the
            # nearest off-line distance are printed as the reading, not gated.
            log(f"N={n}: candidates={candidates}; on-line spread {farthest_on:.1e}, "
                f"nearest off-line {nearest_off:.1e}")
        else:
            log(f"N={n}: candidates=0; nothing on the line, nearest {nearest_off:.6f}")
        assert candidates == expected[n], (n, candidates)
    for name, bond_builder in (("star", star_bonds), ("ring", ring_bonds)):
        values = np.linalg.eigvals(build_liouvillian(4, [GAMMA] * 4, bond_builder(4)))
        linear = orbit_census(values, 4 * GAMMA, "linear")
        composite = orbit_census(values, 4 * GAMMA, "composite")
        log(f"N=4 {name}: candidates={composite.fixed - linear.fixed}")
        assert composite.fixed - linear.fixed == {"star": 110, "ring": 132}[name]
    log()
    log("Semisimplicity on the line (chain): each on-line cluster of multiplicity m")
    log("has m singular values of L - lambda at the eigensolver floor and the next far")
    log("above it, so geometric = algebraic multiplicity (guards Jordan blocks up to")
    log("size 3; a larger one would split wider than the cluster radius 1e-4).")
    for n in (2, 4):
        generator = build_liouvillian(n, [GAMMA] * n, chain_bonds(n))
        values = np.linalg.eigvals(generator)
        sigma = n * GAMMA
        line = values[np.abs(values.real + sigma) <= TOL]
        used = np.zeros(len(line), dtype=bool)
        worst = np.inf
        clusters = 0
        representatives = []
        identity = np.eye(len(generator))
        for i, value in enumerate(line):
            if used[i]:
                continue
            # The cluster radius must exceed the split of a hidden Jordan block of
            # size m, about (eps*||L||)^(1/m): 3e-5 for m = 3, 5e-4 for m = 4. The
            # multiplicity is counted over the whole spectrum, so a partner split
            # off the line is still inside. Distinct levels on the line come as
            # close as 2.0e-4 at N=4 (printed and asserted above the radius), so
            # 1e-4 guards Jordan blocks up to size 3 and no further.
            used |= np.abs(line - value) <= CLUSTER
            representatives.append(value)
            m = int(np.count_nonzero(np.abs(values - value) <= CLUSTER))
            sv = np.linalg.svd(generator - value * identity, compute_uv=False)
            worst = min(worst, sv[-m - 1] / max(sv[-m], 1e-300))
            clusters += 1
        representatives.sort(key=lambda z: (z.imag, z.real))
        spacing = min(abs(a - c) for k, a in enumerate(representatives) for c in representatives[k + 1:])
        log(f"N={n}: {clusters} clusters, smallest ratio (next / largest null) {worst:.1e}, "
            f"closest two levels {spacing:.1e}")
        assert spacing > CLUSTER, spacing
        assert worst > 1e6, (n, worst)
    log()
    log("The cube's reading at one bond: in the orthonormal Pauli basis Re(lambda) =")
    log("-2*gamma*<k_Z>, so the line is the mid-height plane <k_Z> = N/2. At N=2 the")
    log("eight strings of XY weight 1 form a sector L_H does not leave:")
    labels = [a + b for a in "IXYZ" for b in "IXYZ"]
    single = {"I": I2, "X": Xm, "Y": Ym, "Z": Zm}
    basis = np.column_stack([
        np.kron(single[l[0]], single[l[1]]).reshape(-1) / 2.0 for l in labels
    ])
    pauli_generator = basis.conj().T @ build_liouvillian(2, [GAMMA] * 2, chain_bonds(2)) @ basis
    weight = np.array([sum(c in "XY" for c in l) for l in labels])
    inside = np.where(weight == 1)[0]
    outside = np.where(weight != 1)[0]
    leak = max(np.abs(pauli_generator[np.ix_(inside, outside)]).max(),
               np.abs(pauli_generator[np.ix_(outside, inside)]).max())
    sector = np.linalg.eigvals(pauli_generator[np.ix_(inside, inside)])
    log(f"coupling out of the weight-1 sector: {leak:.1e}; its eigenvalues: "
        + ", ".join(f"{v.real:+.3f}{v.imag:+.3f}i" for v in sorted(sector, key=lambda z: z.imag)))
    # Exact, not gated: at one bond every L_H transition joins strings whose XY
    # weights sum to 2 (experiments/PRIMORDIAL_QUBIT_ALGEBRA.md, Step 10) and D is
    # diagonal, so weight 1 couples only to itself and the entry is zero by structure.
    assert leak == 0.0, leak
    log()
    log("Why odd N is empty: L_H keeps the parity of k_Z, so L splits into an even-k_Z")
    log("and an odd-k_Z sector. Pi sends k_Z to N - k_Z: at odd N it carries one sector")
    log("onto the other, at even N each onto itself. Random bond and rate profiles:")
    rng = np.random.default_rng(20261101)
    for n in (3, 4, 5):
        labels = ["".join(t) for t in itertools.product("IXYZ", repeat=n)]
        basis = np.column_stack([
            kron_chain([single[c] for c in l]).reshape(-1) / np.sqrt(2**n) for l in labels
        ])
        weight = np.array([sum(c in "XY" for c in l) for l in labels])
        even = np.where(weight % 2 == 0)[0]
        odd = np.where(weight % 2 == 1)[0]
        for trial in range(5):
            gammas = rng.uniform(0.02, 0.10, n)
            sigma = float(gammas.sum())
            couplings = rng.uniform(0.5, 1.5, n - 1)
            generator = build_liouvillian(n, gammas, chain_bonds(n), couplings)
            pauli = basis.conj().T @ generator @ basis
            leak = max(np.abs(pauli[np.ix_(even, odd)]).max(), np.abs(pauli[np.ix_(odd, even)]).max())
            # Exact: L preserves the parity of i XOR j for |i><j| (the bond terms
            # flip two bits or none, the dephasing is diagonal), and strings of
            # different k_Z parity have different X-mask parity, so their supports
            # in the vectorised basis are disjoint and every product is an exact zero.
            assert leak == 0.0, leak
            spec_even = np.linalg.eigvals(pauli[np.ix_(even, even)])
            spec_odd = np.linalg.eigvals(pauli[np.ix_(odd, odd)])
            on = [int(np.count_nonzero(np.abs(s.real + sigma) <= TOL)) for s in (spec_even, spec_odd)]
            nearest = min(float(np.abs(s.real + sigma)[np.abs(s.real + sigma) > TOL].min())
                          for s in (spec_even, spec_odd))
            line = f"N={n} trial {trial}: on the line even/odd sector {on[0]}/{on[1]}, nearest off {nearest:.1e}"
            if n % 2:
                # Multiplicity-aware, both directions: every odd-sector eigenvalue is
                # assigned to one F1 image of an even-sector eigenvalue.
                images = -spec_even - 2 * sigma
                pairs_found = cKDTree(np.column_stack((images.real, images.imag))).query_ball_point(
                    np.column_stack((spec_odd.real, spec_odd.imag)), TOL)
                rows = [i for i, c in enumerate(pairs_found) for _ in c]
                cols = [j for c in pairs_found for j in c]
                graph = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                                   shape=(len(spec_odd), len(images)))
                assignment = maximum_bipartite_matching(graph, perm_type="column")
                assert np.all(assignment >= 0), "odd sector is not the F1 image of the even one"
                mismatch = float(np.max(np.abs(spec_odd - images[assignment])))
                line += f"; odd sector = F1 image of even, matched with multiplicity, worst {mismatch:.1e}"
                assert on == [0, 0], (n, on)
            else:
                # The even-N half of the reason: each sector is F1-symmetric by itself
                # (orbit_census raises if its spectrum does not pair with multiplicity).
                for spectrum in (spec_even, spec_odd):
                    orbit_census(spectrum, sigma, "linear")
                assert min(on) > 0, (n, on)
            log(line)
    log()
    log("Where the argument needs its premise: a second mirror inside a sector. A bonded")
    log("pair beside an isolated unwatched site (N=3) carries the pair's own mirror,")
    log("and eigenvalues sit on the line at odd N:")
    labels = ["".join(t) for t in itertools.product("IXYZ", repeat=3)]
    gammas = [GAMMA, GAMMA, 0.0]
    values = np.linalg.eigvals(build_liouvillian(3, gammas, [(0, 1)]))
    sigma = float(sum(gammas))
    on = np.abs(values.real + sigma) <= TOL
    oscillating = int(np.count_nonzero(on & (np.abs(values.imag) > TOL)))
    log(f"N=3, bond (0,1), rates (0.05, 0.05, 0): {int(on.sum())} on the line, {oscillating} oscillating")
    # Exact: the spectrum is the pair's N=2 spectrum (10 on its line, 6 oscillating)
    # repeated over the four zero-rate letters of the isolated site.
    assert (int(on.sum()), oscillating) == (40, 24)
    log()
    log("Scope: this is a spectral orbit census. It assigns no mode mechanism")
    log("and no transport, damping, or biological interpretation.")

    return "\n".join(lines) + "\n"


def main():
    output = render_report()
    print(output, end="")
    output_path = RESULTS_DIR / "factor_two_standing_waves.txt"
    output_path.write_text(output, encoding="utf-8")
    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
