"""Bond-to-population response with a fixed positive dephasing rate.

The earlier coherent producer is the gamma=0 *calibration limit*. This one
keeps gamma0>0 fixed and varies hopping J, hence Q=J/gamma0. It runs the
single-excitation (1,1) density block, not a Hamiltonian mode sum. Two
explicit illumination profiles use the same gamma0: every site, or the centre
site alone. The reported ranks belong to the declared finite time sample.

Run: python simulations/handshake_bond_seat_fixed_gamma.py
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm_frechet

from framework.sine_slater import sine_mode_matrix
from handshake_bond_seat_readout import predicted_slopes, seat_inventory


def liouvillian(n: int, hopping: float, gamma0: float,
                watched_sites: tuple[int, ...]) -> np.ndarray:
    """(1,1) block in row-major |a><b| order.

    Every watched site carries the same fixed rate gamma0. For a!=b the
    coherence rate is -2(gamma_a+gamma_b); populations pay zero.
    """
    if n < 3 or not np.isfinite(hopping) or hopping <= 0:
        raise ValueError("require N>=3 and finite positive hopping J")
    if not np.isfinite(gamma0) or gamma0 <= 0:
        raise ValueError("the operating producer requires fixed gamma0>0")
    if not watched_sites or len(set(watched_sites)) != len(watched_sites) or any(
        not 0 <= site < n for site in watched_sites
    ):
        raise ValueError("watched sites must be distinct chain sites")

    h = hopping * (np.diag(np.ones(n - 1), 1) +
                   np.diag(np.ones(n - 1), -1))
    eye = np.eye(n)
    site_rates = np.zeros(n)
    site_rates[list(watched_sites)] = gamma0
    damping = -2.0 * (site_rates[:, None] + site_rates[None, :])
    np.fill_diagonal(damping, 0.0)
    return (-1j * (np.kron(h, eye) - np.kron(eye, h.T)) +
            np.diag(damping.ravel()))


def bond_superoperator(n: int, bond: int) -> np.ndarray:
    """Derivative of -i[H,rho] with respect to additive bond hopping."""
    if not 0 <= bond < n - 1:
        raise ValueError("bond outside chain")
    v = np.zeros((n, n))
    v[bond, bond + 1] = v[bond + 1, bond] = 1.0
    eye = np.eye(n)
    return -1j * (np.kron(v, eye) - np.kron(eye, v.T))


def response_tensor(n: int, times, hopping: float, gamma0: float,
                    watched_sites: tuple[int, ...]) -> np.ndarray:
    """d population_site(t)/d bond_b: [time, bond, readout site]."""
    generator = liouvillian(n, hopping, gamma0, watched_sites)
    grid = np.atleast_1d(np.asarray(times, dtype=float))
    if grid.ndim != 1 or not np.all(np.isfinite(grid)):
        raise ValueError("times must be a finite one-dimensional array")
    psi = sine_mode_matrix(n)[0]
    initial = np.outer(psi, psi).ravel()
    bonds = [bond_superoperator(n, bond) for bond in range(n - 1)]
    readings = np.zeros((len(grid), n - 1, n))
    for ti, time in enumerate(grid):
        for bond, direction in enumerate(bonds):
            tangent = expm_frechet(time * generator, time * direction,
                                    compute_expm=False) @ initial
            populations = np.diag(tangent.reshape(n, n))
            if np.max(np.abs(populations.imag)) > 1e-10:
                raise ArithmeticError("population derivative lost Hermiticity")
            readings[ti, bond, :] = populations.real
    return readings


def independent_tangent_slope(n: int, readout: int, bond: int, time: float,
                               hopping: float, gamma0: float,
                               watched_sites: tuple[int, ...]) -> float:
    """Separate matrix-shaped ODE for (rho, d rho/d bond), with no kron/Frechet."""
    if n < 3 or not np.isfinite(hopping) or hopping <= 0:
        raise ValueError("require N>=3 and finite positive hopping J")
    if not 0 <= readout < n or not 0 <= bond < n - 1:
        raise ValueError("readout or bond outside chain")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(gamma0) or gamma0 <= 0:
        raise ValueError("require gamma0>0")
    if not watched_sites or len(set(watched_sites)) != len(watched_sites) or any(
        not 0 <= site < n for site in watched_sites
    ):
        raise ValueError("watched sites must be distinct chain sites")
    h = np.zeros((n, n))
    for edge in range(n - 1):
        h[edge, edge + 1] = h[edge + 1, edge] = hopping
    v = np.zeros((n, n))
    v[bond, bond + 1] = v[bond + 1, bond] = 1.0
    light = set(watched_sites)
    rates = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            if a != b:
                rates[a, b] = -2.0 * gamma0 * (int(a in light) + int(b in light))
    psi = sine_mode_matrix(n)[0]
    initial = np.concatenate((np.outer(psi, psi).ravel(),
                              np.zeros(n * n, dtype=complex)))

    def rhs(_, state):
        rho = state[:n * n].reshape(n, n)
        tangent = state[n * n:].reshape(n, n)
        rho_dot = -1j * (h @ rho - rho @ h) + rates * rho
        tangent_dot = (-1j * (h @ tangent - tangent @ h +
                               v @ rho - rho @ v) + rates * tangent)
        return np.concatenate((rho_dot.ravel(), tangent_dot.ravel()))

    solution = solve_ivp(rhs, (0.0, time), initial, method="DOP853",
                         rtol=1e-13, atol=1e-15)
    if not solution.success:
        raise ArithmeticError(solution.message)
    tangent = solution.y[n * n:, -1].reshape(n, n)
    return float(tangent[readout, readout].real)


def _sample_ranks(tensor: np.ndarray) -> tuple[list[int], float]:
    """Numerical rank with an explicit roundoff scale and gap to retained SVs."""
    ranks = []
    margins = []
    for site in range(tensor.shape[2]):
        matrix = tensor[:, :, site]
        singular = np.linalg.svd(matrix, compute_uv=False)
        tol = 64.0 * np.finfo(float).eps * max(matrix.shape) * singular[0]
        retained = singular[singular > tol]
        ranks.append(len(retained))
        margins.append(float(retained[-1] / tol))
    return ranks, min(margins)


def fixed_gamma_gate() -> dict:
    """Fixed-gamma Q=10 fixture with Q=1,2 rank controls from changing J."""
    n, gamma0, hopping = 7, 1.0, 10.0
    times = np.array([.02, .04, .07, .10, .15, .23, .35, .50, .70, 1.0])
    profiles = {
        "uniform": tuple(range(n)),
        "single_centre": (3,),
    }
    tensors = {name: response_tensor(n, times, hopping, gamma0, sites)
               for name, sites in profiles.items()}
    uniform_ranks, uniform_margin = _sample_ranks(tensors["uniform"])
    centre_ranks, centre_margin = _sample_ranks(tensors["single_centre"])
    q_rank_profiles = {10.0: {"uniform": uniform_ranks,
                              "single_centre": centre_ranks}}
    q_time_grids = {10.0: times.tolist()}
    margins = [uniform_margin, centre_margin]
    for other_hopping in (1.0, 2.0):
        # Keep J*t fixed; only J/gamma0 changes because gamma0 stays at one.
        other_times = times * (hopping / other_hopping)
        q_time_grids[other_hopping] = other_times.tolist()
        q_rank_profiles[other_hopping] = {}
        for name, sites in profiles.items():
            other_tensor = response_tensor(n, other_times, other_hopping,
                                           gamma0, sites)
            ranks, margin = _sample_ranks(other_tensor)
            q_rank_profiles[other_hopping][name] = ranks
            margins.append(margin)

    checks = []
    for name, readout, bond in (("uniform", 2, 0), ("uniform", 3, 2),
                                ("single_centre", 2, 5),
                                ("single_centre", 3, 0)):
        time = float(times[3])
        predicted = float(tensors[name][3, bond, readout])
        direct = independent_tangent_slope(n, readout, bond, time,
                                            hopping, gamma0, profiles[name])
        error_scale = (np.finfo(float).eps * n *
                       (abs(time) * (1 + abs(hopping * time)) +
                        abs(predicted) + abs(direct)))
        ratio = (abs(predicted - direct) / error_scale
                 if np.isfinite(predicted) and np.isfinite(direct)
                 else float("inf"))
        checks.append({"profile": name, "readout": readout, "bond": bond,
                       "time": time, "frechet": predicted, "direct": direct,
                       "error_ratio": ratio})
    max_ratio = max(row["error_ratio"] for row in checks)

    mirror_residual = max(
        float(np.max(np.abs(tensor[:, b, 3] - tensor[:, n - 2 - b, 3])))
        for tensor in tensors.values() for b in range((n - 1) // 2)
    )
    reflection_budget = (64.0 * np.finfo(float).eps * n *
                         max(float(np.max(np.abs(tensor)))
                             for tensor in tensors.values()))
    permutation = np.eye(n)[::-1]
    lift = np.kron(permutation, permutation)
    exact_reflection = all(
        np.array_equal(lift @ liouvillian(n, hopping, gamma0, sites) @ lift.T,
                       liouvillian(n, hopping, gamma0, sites))
        for sites in profiles.values()
    )
    uniform_sum = float(np.sum(tensors["uniform"][3, :, 2]))
    centre_sum = float(np.sum(tensors["single_centre"][3, :, 2]))
    offcentre_difference = float(abs(tensors["uniform"][3, 0, 2] -
                                     tensors["uniform"][3, 5, 2]))
    coherent_ranks = [seat_inventory(n, site)["rank"] for site in range(n)]
    coherent_sampled_ranks = []
    for site in range(n):
        sampled = predicted_slopes(n, site, times, hopping).T
        singular = np.linalg.svd(sampled, compute_uv=False)
        tol = 64.0 * np.finfo(float).eps * max(sampled.shape) * singular[0]
        coherent_sampled_ranks.append(int(np.count_nonzero(singular > tol)))
    expected_ranks = [6, 6, 6, 3, 6, 6, 6]
    passes = (all(ranks == expected_ranks
                  for rows in q_rank_profiles.values()
                  for ranks in rows.values()) and
              min(margins) > 1e8 and
              coherent_sampled_ranks == coherent_ranks and
              max_ratio <= 16.0 and
              max(abs(row["direct"]) for row in checks) > 0.005 and
              exact_reflection and mirror_residual <= reflection_budget and
              offcentre_difference > 0.001 and
              abs(uniform_sum) > 1e-4 and abs(centre_sum) > 1e-4)
    return {
        "gamma0": gamma0, "hopping": hopping, "times": times.tolist(),
        "uniform_ranks": uniform_ranks,
        "single_centre_ranks": centre_ranks,
        "q_rank_profiles": q_rank_profiles,
        "q_time_grids": q_time_grids,
        "coherent_ranks": coherent_ranks,
        "coherent_sampled_ranks": coherent_sampled_ranks,
        "rank_tolerance_margin": min(margins),
        "tangent_checks": checks, "max_tangent_error_ratio": max_ratio,
        "centre_mirror_residual": mirror_residual,
        "reflection_budget": reflection_budget,
        "offcentre_mirror_difference": offcentre_difference,
        "exact_generator_reflection": exact_reflection,
        "uniform_bond_sum": uniform_sum,
        "single_centre_bond_sum": centre_sum,
        "passes": passes,
    }


def main() -> int:
    result = fixed_gamma_gate()
    def times_text(q):
        return "[" + ", ".join(f"{time:.2f}" for time in result["q_time_grids"][q]) + "]"

    print("=== Fixed positive gamma: bond-to-seat population response ===")
    print("gamma0=1 fixed (time unit 1/gamma0); J=10, Q=J/gamma0=10")
    print("Profiles: uniform gamma0 at every site; single-centre F157 control")
    print(f"Q=10 gamma0*t: {times_text(10.0)}")
    print(f"uniform ranks: {result['uniform_ranks']}")
    print(f"single-centre ranks: {result['single_centre_ranks']}")
    for q in (1.0, 2.0):
        rows = result["q_rank_profiles"][q]
        print(f"Q={int(q)} gamma0*t: {times_text(q)}")
        print(f"Q={int(q)} uniform ranks: {rows['uniform']}")
        print(f"Q={int(q)} single-centre ranks: {rows['single_centre']}")
    print("Across Q=1,2,10, gamma0 stays fixed; J and the listed observation "
          "times change to hold J*t fixed.")
    print(f"coherent-limit ranks: {result['coherent_ranks']} (mathematical control)")
    print(f"coherent-limit ranks on the same Q=10 time grid: "
          f"{result['coherent_sampled_ranks']}")
    print(f"Smallest retained singular-value / rank tolerance: "
          f"{result['rank_tolerance_margin']:.2e}")
    print(f"Independent tangent ODE: max error / (eps*N*scale) "
          f"{result['max_tangent_error_ratio']:.3f} (budget 16)")
    print(f"At gamma0*t=0.1, site 2: uniform-bond sum "
          f"{result['uniform_bond_sum']:+.9f} (uniform illumination), "
          f"{result['single_centre_bond_sum']:+.9f} (centre illumination)")
    print(f"Centre mirror residual {result['centre_mirror_residual']:.2e}; "
          f"off-centre end-bond difference "
          f"{result['offcentre_mirror_difference']:.6f}; "
          f"exact generator reflection {result['exact_generator_reflection']}")
    print("Scope: finite-sample numerical rank at Q=1,2,10; no all-Q rank law, "
          "hardware measurement, or finite-noise recovery claim.")
    print(f"VERDICT: {'PASS' if result['passes'] else 'A CHECK FIRED'}")
    return 0 if result["passes"] else 1


if __name__ == "__main__":
    sys.exit(main())
