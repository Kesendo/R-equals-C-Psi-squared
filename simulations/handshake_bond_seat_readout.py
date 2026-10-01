"""Coherent bond -> mode -> seat readout for an open uniform XY chain.

This producer connects F124's full bond-transition matrix with the sine-mode
nodes counted by F157. The observable is local single-excitation population
from a prepared psi_1 state at gamma=0, not the dissipative PTF alpha painter.
J is the single-excitation hopping: H_SE = J*A_path, equivalent to the
full-chain convention H = (J/2)*sum_b(XX+YY). The defect delta J is additive.

Run: python simulations/handshake_bond_seat_readout.py
"""

from __future__ import annotations

from functools import reduce
from itertools import combinations
from math import gcd, lcm, log10
import sys

import mpmath as mp
import numpy as np
from scipy.linalg import expm, expm_frechet

from framework.sine_slater import sine_dispersion, sine_mode_matrix


def transition_matrix(n: int) -> np.ndarray:
    """F124 M[b,k-1] = <psi_k|V_b|psi_1>, including its strength column."""
    u = sine_mode_matrix(n)
    return (u[:, :-1] * u[0, 1:] + u[:, 1:] * u[0, :-1]).T


def seat_inventory(n: int, seat: int) -> dict:
    """F157 sine nodes and the k=2..N-1 modes in a psi_1 population trace."""
    if n < 3 or not 0 <= seat < n:
        raise ValueError("require N >= 3 and 0 <= seat < N")
    nodes = tuple(k for k in range(1, n + 1) if k * (seat + 1) % (n + 1) == 0)
    visible = tuple(k for k in range(2, n) if k not in nodes)
    h = gcd(seat + 1, n + 1)
    return {
        "h": h,
        "node_modes": nodes,
        "visible_modes": visible,
        "rank": n - h - 1,
    }


def _high_precision_modal_slope(n: int, seat: int, bond: int, time: float,
                                j_coupling: float, modes: tuple[int, ...]) -> float:
    """Rebuild the sine-mode sum when cancellation defeats double precision.

    Reusing double-precision M or u here would preserve a spurious t^2 term,
    even if the subsequent sum were done with arbitrary precision.
    """
    if time == 0.0 or not modes:
        return 0.0
    small_time_digits = max(0.0, -log10(abs(j_coupling)) - log10(abs(time)))
    precision = max(80, int(40 + n * small_time_digits))
    with mp.workdps(precision):
        scale = mp.sqrt(mp.mpf(2) / (n + 1))
        step = mp.pi / (n + 1)
        coupling = mp.mpf(j_coupling)
        duration = mp.mpf(time)

        def amplitude(k: int, site: int):
            return scale * mp.sin(k * (site + 1) * step)

        carrier_at_seat = amplitude(1, seat)
        carrier_at_bond = amplitude(1, bond)
        carrier_at_next = amplitude(1, bond + 1)
        carrier_energy = 2 * coupling * mp.cos(step)
        total = mp.mpf(0)
        for k in modes:
            matrix_element = (amplitude(k, bond) * carrier_at_next +
                              amplitude(k, bond + 1) * carrier_at_bond)
            gap = carrier_energy - 2 * coupling * mp.cos(k * step)
            factor = 2 * mp.sin(gap * duration / 2) ** 2 / gap
            total += amplitude(k, seat) * matrix_element * factor
        return float(2 * carrier_at_seat * total)


def predicted_slopes(n: int, seat: int, times, j_coupling: float = 1.0) -> np.ndarray:
    """First-order d population_j(t)/d bond_b, one row per bond.

    The prepared state is the unperturbed psi_1. A bond perturbation has
    V_b[b,b+1] = V_b[b+1,b] = 1, so the derivative is with respect to an
    additive hopping change, not a fractional change in J.
    """
    inventory = seat_inventory(n, seat)
    if not np.isfinite(j_coupling) or j_coupling <= 0:
        raise ValueError("require a finite positive uniform hopping J")
    time_grid = np.atleast_1d(np.asarray(times, dtype=float))
    if time_grid.ndim != 1 or not np.all(np.isfinite(time_grid)):
        raise ValueError("times must be a finite one-dimensional array")
    u = sine_mode_matrix(n)
    energy = sine_dispersion(n, j_coupling)
    m = transition_matrix(n)
    mode_indices = np.array([k - 1 for k in inventory["visible_modes"]], dtype=int)
    gaps = energy[0] - energy[mode_indices]
    # 1-cos(x) cancels to zero for small nonzero x; 2*sin(x/2)^2 keeps the
    # quadratic-time response that the direct Hamiltonian derivative sees.
    time_factors = 2.0 * np.sin(np.outer(gaps, time_grid) / 2.0) ** 2 / gaps[:, None]
    terms = (2.0 * u[0, seat] *
             (m[:, mode_indices] * u[mode_indices, seat]))[:, :, None] * time_factors[None, :, :]
    slopes = np.sum(terms, axis=1)
    absolute_sum = np.sum(np.abs(terms), axis=1)
    # Terms at distant seats can cancel through t^2 while the true response
    # starts at t^4. Detect loss of relative precision and rebuild *all*
    # modal ingredients at higher precision, not merely the final sum.
    unstable = ((absolute_sum > 0.0) &
                (np.abs(slopes) < np.sqrt(np.finfo(float).eps) * absolute_sum))
    for bond, time_index in zip(*np.nonzero(unstable)):
        slopes[bond, time_index] = _high_precision_modal_slope(
            n, seat, int(bond), float(time_grid[time_index]), j_coupling,
            inventory["visible_modes"]
        )
    return slopes


def modal_rank_report(n: int, seat: int) -> dict:
    """Read numerical ranks of the F124 matrix after the F157 seat mask.

    The integer predictions come from gcd arithmetic; SVD is a separate
    finite-N check, with its default error scale reported rather than used as
    a proof of the all-N statement.
    """
    inventory = seat_inventory(n, seat)
    u = sine_mode_matrix(n)
    m = transition_matrix(n)
    b = m * u[:, seat][None, :]
    full = np.linalg.svd(b, compute_uv=False)
    location = np.linalg.svd(b[:, 1:], compute_uv=False)
    eps = np.finfo(float).eps
    # A nodal seat can make the *entire* location matrix rounding-sized. Its
    # own largest singular value is then no error scale. Use the unmasked
    # transition norm and the seat's largest mode amplitude instead.
    tol = (4.0 * eps * n * np.linalg.svd(m, compute_uv=False)[0] *
           np.max(np.abs(u[:, seat])))
    full_tol = location_tol = tol
    full_rank = int(np.count_nonzero(full > full_tol))
    location_rank = int(np.count_nonzero(location > location_tol))
    expected_full = n - inventory["h"]
    expected_location = inventory["rank"]
    return {
        "expected_full_rank": expected_full,
        "expected_location_rank": expected_location,
        "observed_full_rank": full_rank,
        "observed_location_rank": location_rank,
        "full_tolerance": float(full_tol),
        "location_tolerance": float(location_tol),
        "passes": (full_rank == expected_full and
                   location_rank == expected_location and
                   len(inventory["visible_modes"]) == expected_location),
    }


def direct_population_slope(n: int, seat: int, bond: int, time: float,
                            j_coupling: float = 1.0,
                            prepared_mode: int = 1) -> float:
    """Independent Hamiltonian/Fréchet derivative of local population.

    This route diagonalizes the real tridiagonal H only to prepare its mode;
    the time derivative uses a matrix exponential, never the F124 dictionary,
    F157 mask, or F2b dispersion.
    """
    seat_inventory(n, seat)
    if not 0 <= bond < n - 1 or not 1 <= prepared_mode <= n:
        raise ValueError("bond or prepared mode outside the chain")
    if not np.isfinite(j_coupling) or j_coupling <= 0 or not np.isfinite(time):
        raise ValueError("require finite time and positive uniform hopping J")
    h = j_coupling * (
        np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
    )
    v = np.zeros((n, n))
    v[bond, bond + 1] = v[bond + 1, bond] = 1.0
    _, eigenvectors = np.linalg.eigh(h)
    initial = eigenvectors[:, n - prepared_mode]
    propagator, derivative = expm_frechet(
        -1j * time * h, -1j * time * v, compute_expm=True
    )
    amplitude = (propagator @ initial)[seat]
    response = (derivative @ initial)[seat]
    return float(2.0 * np.real(np.conj(amplitude) * response))


def response_gate() -> dict:
    """Compare the modal prediction to independent Hamiltonian derivatives.

    Cases include a truly zero N=3 centre response, quadratic and fourth-order
    short-time germs, all N=7 bonds at two seats and times, a nodal off-centre
    seat, and fixed-Jt runs over two decades of J.
    The Fréchet input is t*V_b, so the residual model is eps_machine * N
    times |t|*(1+|Jt|) plus the response magnitudes. The run prints the
    measured maximum ratio next to the 16-fold numerical budget.
    """
    cases = [
        (3, 0, 0, 0.7, 1.0),
        (3, 1, 0, 1.0, 1.0),
        (7, 2, 0, 1e-9, 1.0),
        (7, 3, 0, 1e-8, 1.0),
        (5, 1, 0, 1.37, 1.0),
        (9, 1, 6, 0.73, 0.7),
        (7, 2, 0, 13.0, 0.1),
        (7, 2, 0, 0.13, 10.0),
        (3, 1, 0, 1000.0, 0.001),
        (3, 1, 0, 10000.0, 0.0001),
    ]
    cases.extend((7, seat, bond, time, 1.0)
                 for seat in (2, 3) for bond in range(6)
                 for time in (1.0, 1.3))
    readings = []
    for n, seat, bond, time, j_coupling in cases:
        predicted = float(predicted_slopes(n, seat, [time], j_coupling)[bond, 0])
        direct = direct_population_slope(n, seat, bond, time, j_coupling)
        if np.isfinite(predicted) and np.isfinite(direct):
            error = abs(predicted - direct)
            error_model = (np.finfo(float).eps * n *
                           (abs(time) * (1.0 + abs(j_coupling * time)) +
                            abs(predicted) + abs(direct)))
            error_ratio = error / error_model
        else:
            error = error_ratio = float("inf")
        readings.append({
            "n": n, "seat": seat, "bond": bond, "time": time,
            "J": j_coupling, "predicted": predicted, "direct": direct,
            "error": error, "error_ratio": error_ratio,
        })
    max_ratio = max(row["error_ratio"] for row in readings)
    finite_rows = [row for row in readings
                   if np.isfinite(row["predicted"]) and np.isfinite(row["direct"])]
    max_direct = max((abs(row["direct"]) for row in finite_rows), default=0.0)
    all_finite = len(finite_rows) == len(readings)
    fourth_order = next(row for row in readings
                        if (row["n"], row["seat"], row["bond"], row["time"])
                        == (7, 3, 0, 1e-8))
    leading = -np.sqrt(2.0) * (1e-8) ** 4 / 48.0
    fourth_order_relative_error = max(
        abs(fourth_order["predicted"] - leading),
        abs(fourth_order["direct"] - leading)
    ) / abs(leading)
    return {
        "readings": readings,
        "max_error_ratio": max_ratio,
        "max_direct_slope": max_direct,
        "fourth_order_leading": leading,
        "fourth_order_reading": fourth_order,
        "fourth_order_relative_error": fourth_order_relative_error,
        "all_finite": all_finite,
        "passes": (all_finite and max_ratio <= 16.0 and max_direct > 0.01
                   and fourth_order_relative_error <= 1e-6),
    }


def reflection_report(paired_bond: int = 5) -> dict:
    """N=7 finite-defect centre mirror identity and off-centre control.

    The intentionally exposed paired-bond input lets a wrong mirror pair pass
    through the same comparison door; bond 4 must fail the bond-0/5 claim.
    """
    n, time, delta = 7, 1.3, 0.2
    if not 0 <= paired_bond < n - 1:
        raise ValueError("paired bond outside the N=7 chain")
    h = np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
    _, eigenvectors = np.linalg.eigh(h)
    initial = eigenvectors[:, -1]

    def moved(bond):
        result = h.copy()
        result[bond, bond + 1] += delta
        result[bond + 1, bond] += delta
        return result

    left_h, right_h = moved(0), moved(paired_bond)
    left = np.abs(expm(-1j * time * left_h) @ initial) ** 2
    right = np.abs(expm(-1j * time * right_h) @ initial) ** 2
    centre_difference = float(abs(left[3] - right[3]))
    offcentre_difference = float(abs(left[2] - right[2]))
    reflected_offseat_error = float(abs(left[2] - right[4]))
    budget = 16.0 * np.finfo(float).eps * n * (1.0 + time * (1.0 + delta))
    exact_h_reflection = bool(np.array_equal(left_h[::-1, ::-1], right_h))
    return {
        "paired_bond": paired_bond,
        "centre_difference": centre_difference,
        "offcentre_difference": offcentre_difference,
        "reflected_offseat_error": reflected_offseat_error,
        "roundoff_budget": budget,
        "exact_h_reflection": exact_h_reflection,
        "passes": (exact_h_reflection and centre_difference <= budget and
                   reflected_offseat_error <= budget and
                   offcentre_difference > 0.01),
    }


def set_inventory(n: int, seats) -> dict:
    """The same count for a SET of read sites, populations read at all of them.

    Mode k reaches the stacked trace unless it has a node at EVERY site of S;
    the modes dark at all of S are the multiples of (N+1)/G with
    G = gcd(N+1, {j+1 : j in S}), G-1 of them, none equal to 1 or N.
    """
    seats = tuple(sorted(set(seats)))
    if not seats:
        raise ValueError("a set of read sites needs at least one site")
    for seat in seats:
        seat_inventory(n, seat)
    g = reduce(gcd, (seat + 1 for seat in seats), n + 1)
    dark = tuple(k for k in range(1, n + 1)
                 if all(k * (seat + 1) % (n + 1) == 0 for seat in seats))
    return {"G": g, "dark_modes": dark, "rank": n - 1 - g,
            "dark_count_matches": len(dark) == g - 1 and 1 not in dark and n not in dark}


def set_rank_report(n: int, seats) -> dict:
    """Modal route for a set: the stacked masked F124 matrices, rows (j, k).

    Mode k contributes u_1(j) u_k(j) M[.,k] f_k(t) at site j, and the f_k are
    independent, so the trace rank is the rank of the matrix with one row per
    (site, mode k = 2..N). Each visible mode adds the one bond direction
    M[.,k] however many sites of S see it, so the rank counts modes. Error
    scale: modal_rank_report's bound for one site's block, times sqrt(|S|)
    for stacking |S| such blocks (their errors add in the Frobenius norm).
    The margins on both sides are returned so the run can print them.
    """
    seats = tuple(sorted(set(seats)))
    inventory = set_inventory(n, seats)
    u = sine_mode_matrix(n)
    m = transition_matrix(n)
    stacked = np.vstack([(m * u[:, seat][None, :])[:, 1:].T for seat in seats])
    sv = np.linalg.svd(stacked, compute_uv=False)
    tol = (4.0 * np.finfo(float).eps * n * np.linalg.svd(m, compute_uv=False)[0] *
           max(np.max(np.abs(u[:, seat])) for seat in seats) * np.sqrt(len(seats)))
    observed = int(np.count_nonzero(sv > tol))
    expected = inventory["rank"]
    return {"seats": seats, "expected": expected, "observed": observed,
            "kept_over_tol": float(sv[expected - 1] / tol) if expected > 0 else float("inf"),
            "dropped_over_tol": float(sv[expected] / tol) if expected < len(sv) else 0.0,
            "passes": observed == expected}


def _all_site_slopes(n: int, times, j_coupling: float = 1.0) -> np.ndarray:
    """Fréchet slopes at every site: D[t, site, bond], from the propagator.

    The same independent route as direct_population_slope (no F124 matrix, no
    F157 mask, no sine dispersion), read at all sites from one derivative.
    """
    h = j_coupling * (np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1))
    initial = np.linalg.eigh(h)[1][:, n - 1]
    out = np.zeros((len(times), n, n - 1))
    for ti, time in enumerate(times):
        for bond in range(n - 1):
            v = np.zeros((n, n))
            v[bond, bond + 1] = v[bond + 1, bond] = 1.0
            propagator, derivative = expm_frechet(
                -1j * time * h, -1j * time * v, compute_expm=True)
            amplitude = propagator @ initial
            response = derivative @ initial
            out[ti, :, bond] = 2.0 * np.real(np.conj(amplitude) * response)
    return out


def set_trace_gate(n_range=range(4, 12), max_size: int = 3) -> dict:
    """Rank of the sampled time trace at a SET of read sites, two routes.

    The Fréchet route stacks the slopes of every site in S at 2N times into a
    (|S|*2N) x (N-1) matrix. Its entry error model is the response gate's,
    eps*N*(t(1+Jt) + |modal| + |direct|) at budget 16, and it is checked here
    entry by entry against the modal slopes on the whole grid it is used on,
    every site, bond and time. A singular value is then resolved above that
    bound carried through the Frobenius norm (Weyl), 16*eps*N*scale*sqrt(size).
    The run prints the separation on both sides, and the dropped values in
    units of the SVD's own floor eps*sigma_max*max(shape), where they sit.
    Two wrong rules are run against the OBSERVED Frechet rank beside the gcd:
    the lcm of the site divisors, and the largest single-site rank (no gain
    from reading several sites). Both agree on every single site and must miss
    on some sets. The dark-mode set found by a direct node check must have
    G-1 members, none of them k=1 or k=N.
    """
    eps = np.finfo(float).eps
    rows, modal_fail, frechet_fail, dark_fail = 0, 0, 0, 0
    lcm_misses, max_single_misses = 0, 0
    min_kept_ratio, max_dropped_ratio, max_dropped_floor = float("inf"), 0.0, 0.0
    modal_min_kept, modal_max_dropped = float("inf"), 0.0
    max_entry_ratio = 0.0
    for n in n_range:
        times = [0.37 * (i + 1) for i in range(2 * n)]
        slopes = _all_site_slopes(n, times)
        for seat in range(n):
            modal = predicted_slopes(n, seat, times)
            direct = slopes[:, seat, :].T
            model = eps * n * (np.array(times)[None, :] * (1.0 + np.array(times)[None, :]) +
                               np.abs(modal) + np.abs(direct))
            max_entry_ratio = max(max_entry_ratio, float(np.max(np.abs(modal - direct) / model)))
        scale = max(t * (1.0 + t) for t in times) + 2.0 * float(np.max(np.abs(slopes)))
        single = [seat_inventory(n, seat)["rank"] for seat in range(n)]
        for size in range(1, max_size + 1):
            for seats in combinations(range(n), size):
                inventory = set_inventory(n, seats)
                expected = inventory["rank"]
                dark_fail += not inventory["dark_count_matches"]
                modal = set_rank_report(n, seats)
                modal_fail += not modal["passes"]
                modal_min_kept = min(modal_min_kept, modal["kept_over_tol"])
                modal_max_dropped = max(modal_max_dropped, modal["dropped_over_tol"])
                block = slopes[:, list(seats), :].reshape(-1, n - 1)
                sv = np.linalg.svd(block, compute_uv=False)
                tol = 16.0 * eps * n * scale * np.sqrt(block.size)
                observed = int(np.count_nonzero(sv > tol))
                frechet_fail += observed != expected
                g_lcm = gcd(n + 1, lcm(*(seat + 1 for seat in seats)))
                lcm_misses += observed != n - 1 - g_lcm
                max_single_misses += observed != max(single[seat] for seat in seats)
                if expected > 0:
                    min_kept_ratio = min(min_kept_ratio, float(sv[expected - 1] / tol))
                if expected < len(sv):
                    max_dropped_ratio = max(max_dropped_ratio, float(sv[expected] / tol))
                    max_dropped_floor = max(max_dropped_floor, float(
                        sv[expected] / (eps * sv[0] * max(block.shape))))
                rows += 1
    return {"sets": rows, "modal_failures": modal_fail, "dark_count_failures": dark_fail,
            "frechet_failures": frechet_fail, "lcm_misses": lcm_misses,
            "max_single_misses": max_single_misses,
            "max_entry_error_ratio": max_entry_ratio,
            "min_kept_over_bound": min_kept_ratio,
            "max_dropped_over_bound": max_dropped_ratio,
            "max_dropped_over_svd_floor": max_dropped_floor,
            "modal_min_kept_over_tol": modal_min_kept,
            "modal_max_dropped_over_tol": modal_max_dropped,
            "passes": (modal_fail == 0 and frechet_fail == 0 and dark_fail == 0
                       and max_entry_ratio <= 16.0
                       and lcm_misses > 0 and max_single_misses > 0)}


def main() -> int:
    """Run the scoped workflow and independent controls; return 1 on a firing."""
    print("=== F124 -> F157 coherent bond/seat producer (gamma=0) ===")
    print("Convention: H_SE=J_hop*A_path = (J_hop/2)*sum_b(XX+YY); "
          "delta J is additive hopping.")
    all_ranks = [modal_rank_report(n, seat)
                 for n in range(3, 12) for seat in range(n)]
    ranks_ok = all(row["passes"] for row in all_ranks)
    print(f"SVD rank checks N=3..11: {sum(row['passes'] for row in all_ranks)}"
          f"/{len(all_ranks)}")
    profile = [seat_inventory(7, seat)["rank"] for seat in range(7)]
    centre = seat_inventory(7, 3)
    outside = seat_inventory(7, 2)
    print(f"N=7 location ranks: {profile}")
    print(f"centre modes: {centre['visible_modes']}  nodes: {centre['node_modes']}")
    print(f"seat 2 modes: {outside['visible_modes']}  nodes: {outside['node_modes']}")

    sets = set_trace_gate()
    print(f"Read-site sets, sizes 1..3, N=4..11: {sets['sets']} sets; rank "
          f"N-1-gcd(N+1, {{j+1}}) misses: modal {sets['modal_failures']}, "
          f"Frechet trace {sets['frechet_failures']}; dark-mode count "
          f"!= G-1: {sets['dark_count_failures']}")
    print(f"  controls against the observed rank: lcm of the divisors misses "
          f"{sets['lcm_misses']}, "
          f"largest single-site rank misses {sets['max_single_misses']}")
    print(f"  Frechet entries vs modal on the whole grid: max error / model "
          f"{sets['max_entry_error_ratio']:.3f} (budget 16)")
    print(f"  Frechet margins: smallest kept singular value "
          f"{sets['min_kept_over_bound']:.3g} x the bound, largest dropped "
          f"{sets['max_dropped_over_bound']:.3g} x "
          f"({sets['max_dropped_over_svd_floor']:.3g} x eps*sigma_max*max(shape))")
    print(f"  modal margins: smallest kept {sets['modal_min_kept_over_tol']:.3g} x tol, "
          f"largest dropped {sets['modal_max_dropped_over_tol']:.3g} x tol")
    pair = set_inventory(11, (2, 7))
    print(f"N=11 sites {{2,7}}: G={pair['G']}, dark modes {pair['dark_modes']}, "
          f"rank {pair['rank']} = N-2 (each alone: "
          f"{seat_inventory(11, 2)['rank']}, {seat_inventory(11, 7)['rank']})")

    response = response_gate()
    print(f"Independent expm-Frechet response: max error / "
          f"(eps*N*scale) = {response['max_error_ratio']:.3f} "
          f"(budget 16; {'PASS' if response['passes'] else 'FIRED'})")
    all_bond_rows = [row for row in response["readings"]
                     if row["n"] == 7 and row["seat"] in (2, 3)
                     and row["time"] in (1.0, 1.3)]
    print(f"{len(all_bond_rows)} N=7 all-bond comparisons "
          f"(seats 2,3; t=1,1.3; J=1)")
    germ = response["fourth_order_reading"]
    print(f"N=7 centre bond=0 t=1e-8, J=1: "
          f"modal={germ['predicted']:+.6e} "
          f"direct={germ['direct']:+.6e} "
          f"leading={response['fourth_order_leading']:+.6e}; "
          f"relative error {response['fourth_order_relative_error']:.2e} "
          f"(budget 1e-6)")
    for row in response["readings"]:
        if row["n"] == 7 and row["time"] == 1.0 and row["bond"] in (0, 5):
            print(f"  N=7 J=1, t=1 seat={row['seat']} bond={row['bond']}: "
                  f"modal={row['predicted']:+.9f} "
                  f"direct={row['direct']:+.9f}")

    reflection = reflection_report()
    print(f"Finite defect delta=0.2, t=1.3: centre mirror diff "
          f"{reflection['centre_difference']:.2e}; "
          f"seat-2 mirror diff {reflection['offcentre_difference']:.6f}; "
          f"exact H reflection {reflection['exact_h_reflection']}")

    uniform = predicted_slopes(7, 2, [0.7, 1.3, 2.1]).sum(axis=0)
    uniform_residual = float(np.max(np.abs(uniform)))
    uniform_budget = 16.0 * np.finfo(float).eps * 7
    print(f"Uniform bond rescaling, seat 2: max linear response "
          f"{uniform_residual:.2e}")

    psi1_slope = float(predicted_slopes(7, 2, [1.0])[0, 0])
    psi2_slope = direct_population_slope(7, 2, 0, 1.0, prepared_mode=2)
    preparation_mismatch = abs(psi1_slope - psi2_slope)
    print(f"Changed preparation psi_2: mismatch with psi_1 formula "
          f"{preparation_mismatch:.6f}")
    print("Scope: ideal unitary single-excitation population; PTF alpha and "
          "finite-gamma purity are outside this result.")

    passed = (ranks_ok and sets["passes"] and response["passes"] and reflection["passes"] and
              uniform_residual <= uniform_budget and
              preparation_mismatch > 0.02)
    print(f"VERDICT: {'PASS' if passed else 'A CHECK FIRED'}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
