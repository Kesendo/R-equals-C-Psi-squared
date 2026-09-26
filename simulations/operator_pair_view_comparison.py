"""Finite comparison of pair-flow readout ranks and Hilbert reflection motion."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from numbers import Integral
from pathlib import Path
from typing import Sequence

import numpy as np
import sympy as sp
from scipy.linalg import expm

from operator_pair_flow_atlas import exact_heisenberg_se_h, exact_readout_rank


PRIMES = (1_000_003, 1_000_033)
CASES = {
    (3, (0, 1, 0)): "centre-only",
    (3, (1, 1, 1)): "uniform",
    (3, (1, 1, 0)): "broken",
    (3, (0, 3, 0)): "matched-total-control",
    (4, (0, 1, 1, 0)): "symmetric",
    (4, (1, 1, 1, 0)): "broken",
    (5, (0, 0, 1, 0, 0)): "centre-only",
    (5, (1, 0, 1, 0, 0)): "broken",
}


@dataclass(frozen=True)
class RankMapCertificate:
    n: int
    gamma: tuple[int, ...]
    ranks: tuple[tuple[int, ...], ...]
    modular_ranks: tuple[tuple[tuple[int, ...], ...], ...]
    upper_bounds: tuple[tuple[int, ...], ...]
    method: str
    generator_rank: int | None


def _integer_profile(n: int, gamma_per_site: Sequence[int]) -> tuple[int, ...]:
    if n < 2 or len(gamma_per_site) != n:
        raise ValueError("n >= 2 and one site rate per site required")
    if any(not isinstance(rate, Integral) or rate < 0 for rate in gamma_per_site):
        raise ValueError("nonnegative integer site rates required")
    return tuple(int(rate) for rate in gamma_per_site)


def _single_excitation_h(n: int) -> np.ndarray:
    """Independent integer H for J sum(XX+YY+ZZ) on an open chain."""
    h = np.zeros((n, n), dtype=np.int64)
    for site in range(n):
        degree = int(site > 0) + int(site < n - 1)
        h[site, site] = (n - 1) - 2 * degree
    for site in range(n - 1):
        h[site, site + 1] = h[site + 1, site] = 2
    return h


def integer_hermitian_generator(
    n: int, gamma_per_site: Sequence[int]
) -> np.ndarray:
    """Exact integer L in coordinates (populations, Re upper, Im upper)."""
    gamma = _integer_profile(n, gamma_per_site)
    h = _single_excitation_h(n)
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    labels = [("population", a, a) for a in range(n)]
    labels += [("real", a, b) for a, b in pairs]
    labels += [("imag", a, b) for a, b in pairs]
    d = n * n
    generator = np.zeros((d, d), dtype=np.int64)

    for source, (kind, a, b) in enumerate(labels):
        real = np.zeros((n, n), dtype=np.int64)
        imag = np.zeros((n, n), dtype=np.int64)
        if kind == "population":
            real[a, a] = 1
        elif kind == "real":
            real[a, b] = real[b, a] = 1
        else:
            imag[a, b], imag[b, a] = 1, -1

        # -i[H, real+i*imag] = [H, imag] - i[H, real].
        image_real = h @ imag - imag @ h
        image_imag = -(h @ real - real @ h)
        for left in range(n):
            for right in range(n):
                if left != right:
                    cost = 2 * (gamma[left] + gamma[right])
                    image_real[left, right] -= cost * real[left, right]
                    image_imag[left, right] -= cost * imag[left, right]
        for target, (target_kind, left, right) in enumerate(labels):
            generator[target, source] = (
                image_real[left, right]
                if target_kind != "imag"
                else image_imag[left, right]
            )
    return generator


def _rank_mod_prime(matrix: np.ndarray, prime: int) -> int:
    work = np.array(matrix, dtype=np.int64, copy=True) % prime
    rows, cols = work.shape
    rank = 0
    for col in range(cols):
        pivot = next((row for row in range(rank, rows) if work[row, col]), None)
        if pivot is None:
            continue
        if pivot != rank:
            work[[rank, pivot]] = work[[pivot, rank]]
        inverse = pow(int(work[rank, col]), -1, prime)
        work[rank, col:] = (work[rank, col:] * inverse) % prime
        for row in range(rank + 1, rows):
            factor = work[row, col]
            if factor:
                work[row, col:] = (
                    work[row, col:] - factor * work[rank, col:]
                ) % prime
        rank += 1
        if rank == rows:
            break
    return rank


def _modular_rank_map(generator: np.ndarray, n: int, prime: int) -> tuple[tuple[int, ...], ...]:
    """Field ranks are lower bounds on rational moment-Hankel ranks."""
    d = n * n
    reduced = generator % prime
    result = []
    hankel_indices = np.add.outer(np.arange(d), np.arange(d))
    for preparation in range(n):
        state = np.zeros(d, dtype=np.int64)
        state[preparation] = 1
        moments = np.zeros((2 * d - 1, n), dtype=np.int64)
        for order in range(2 * d - 1):
            moments[order, :] = state[:n]
            state = (reduced @ state) % prime
        result.append(tuple(
            _rank_mod_prime(moments[hankel_indices, readout], prime)
            for readout in range(n)
        ))
    return tuple(result)


def _centre_symmetry_certified(n: int, gamma: tuple[int, ...]) -> bool:
    centre = n // 2
    reversal = sp.zeros(n)
    for site in range(n):
        reversal[site, n - 1 - site] = 1
    h = exact_heisenberg_se_h(n, 1)
    jump = sp.eye(n)
    jump[centre, centre] = -1
    return (
        n % 2 == 1
        and gamma == tuple(int(site == centre) for site in range(n))
        and h * reversal == reversal * h
        and jump * reversal == reversal * jump
    )


def certified_rank_map(n: int, gamma_per_site: Sequence[int]) -> RankMapCertificate:
    """Exact ranks for seven predeclared cases and one matched-rate control."""
    gamma = _integer_profile(n, gamma_per_site)
    if (n, gamma) not in CASES:
        raise ValueError("only the predeclared cases and matched-rate control have exact-rank certificates")
    generator = integer_hermitian_generator(n, gamma)
    modular_ranks = tuple(_modular_rank_map(generator, n, p) for p in PRIMES)
    d = n * n

    if n <= 4:
        h = exact_heisenberg_se_h(n, 1)
        ranks = tuple(tuple(
            exact_readout_rank(h, gamma, preparation=a, readout=b).rank
            for b in range(n)
        ) for a in range(n))
        if any(residue != ranks for residue in modular_ranks):
            raise ArithmeticError("modular rank disagrees with exact small-N rank")
        upper = tuple((d,) * n for _ in range(n))
        return RankMapCertificate(
            n, gamma, ranks, modular_ranks, upper, "exact-sympy-and-modular", None
        )

    if gamma == (0, 0, 1, 0, 0):
        if not _centre_symmetry_certified(n, gamma):
            raise ArithmeticError("centre-only Hilbert parity certificate failed")
        generator_rank = int(sp.Matrix(generator.tolist()).rank())
        if generator_rank != 22:
            raise ArithmeticError("expected three independent stationary directions")
        even_operator_dim = ((n + 1) // 2) ** 2
        upper = tuple(tuple(
            even_operator_dim if a == n // 2 or b == n // 2
            else generator_rank + 1
            for b in range(n)
        ) for a in range(n))
    else:
        generator_rank = None
        upper = tuple((d,) * n for _ in range(n))

    if any(residue != upper for residue in modular_ranks):
        raise ArithmeticError("modular lower bound does not reach exact upper bound")
    return RankMapCertificate(
        n, gamma, upper, modular_ranks, upper,
        "modular-lower-plus-exact-upper", generator_rank,
    )


def _twice_odd_readout() -> np.ndarray:
    """2*Tr(P_odd rho) for N=3 in the integer Hermitian coordinates."""
    readout = np.zeros(9, dtype=np.int64)
    readout[0] = readout[2] = 1
    readout[4] = -2  # Re rho_02, after populations and the (0,1) coordinate.
    return readout


def odd_weight_derivatives(gamma_per_site: Sequence[int], count: int) -> tuple[sp.Rational, ...]:
    """Exact derivatives at t=0 from a centre population preparation, N=3."""
    gamma = _integer_profile(3, gamma_per_site)
    if count < 1:
        raise ValueError("positive derivative count required")
    generator = sp.Matrix(integer_hermitian_generator(3, gamma).tolist())
    state = sp.zeros(9, 1)
    state[1] = 1
    readout = sp.Matrix([int(x) for x in _twice_odd_readout()]).T
    values = []
    for _ in range(count):
        values.append(sp.Rational((readout * state)[0], 2))
        state = generator * state
    return tuple(values)


def odd_weight_curve(gamma_per_site: Sequence[int], times: Sequence[float]) -> tuple[float, ...]:
    """Hilbert reflection-odd population after preparing site 1 at N=3."""
    generator = integer_hermitian_generator(3, gamma_per_site).astype(float)
    initial = np.zeros(9)
    initial[1] = 1
    readout = _twice_odd_readout().astype(float) / 2
    return tuple(float(readout @ (expm(float(t) * generator) @ initial)) for t in times)


def centre_return_signature(gamma_per_site: Sequence[int]) -> dict[str, float | int]:
    """N=3 centre population: an exact third derivative and one finite-time read."""
    generator = integer_hermitian_generator(3, gamma_per_site)
    state = np.zeros(9, dtype=np.int64)
    state[1] = 1
    initial = state.copy()
    for _ in range(3):
        state = generator @ state
    return {
        "third_derivative_at_zero": int(state[1]),
        "population_at_t_0_5": float((expm(0.5 * generator.astype(float)) @ initial)[1]),
    }


def _draw_rank_map(ax, ranks: tuple[tuple[int, ...], ...], *, vmin: int, vmax: int, title: str):
    image = ax.imshow(ranks, cmap="viridis", vmin=vmin, vmax=vmax)
    n = len(ranks)
    ax.set(
        title=title,
        xlabel="Mess-Site b",
        ylabel="Präparation a",
        xticks=range(n),
        yticks=range(n),
    )
    threshold = (vmin + vmax) / 2
    for a in range(n):
        for b in range(n):
            value = ranks[a][b]
            ax.text(b, a, str(value), ha="center", va="center",
                    color="black" if value > threshold else "white",
                    fontsize=11, fontweight="bold")
    return image


def _render_n3_comparison(cases: list[RankMapCertificate], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    palette = {"centre-only": "#2864a6", "uniform": "#6b3c91", "broken": "#c45c3c"}
    labels = {
        "centre-only": "nur Mitte γ=(0,1,0)",
        "uniform": "uniform γ=(1,1,1)",
        "broken": "einseitig γ=(1,1,0)",
    }
    fig = plt.figure(figsize=(13.5, 8.1))
    grid = fig.add_gridspec(2, 4, width_ratios=(1, 1, 1, 0.055),
                            height_ratios=(1, 1.05), hspace=0.4, wspace=0.32)
    fig.suptitle("N=3: gleicher Signalrang kann anderen Paritätsfluss verbergen", fontsize=14)
    image = None
    for col, case in enumerate(cases):
        label = CASES[(case.n, case.gamma)]
        image = _draw_rank_map(
            fig.add_subplot(grid[0, col]), case.ranks,
            vmin=4, vmax=9, title=labels[label],
        )
    fig.colorbar(image, cax=fig.add_subplot(grid[0, 3]),
                 label="exakter Signalrang")

    ax = fig.add_subplot(grid[1, :3])
    times = np.linspace(0, 1.4, 141)
    for case in cases:
        label = CASES[(case.n, case.gamma)]
        values = odd_weight_curve(case.gamma, times)
        ax.plot(times, values, color=palette[label], lw=2.4, label=labels[label])
    ax.axvline(0.5, color="#77838a", lw=1, ls="--")
    ax.set(
        title="Mitte präpariert: Gewicht im ungeraden Hilbertraum",
        xlabel="Zeit t",
        ylabel="Tr(P_ungerade ρ(t))",
        xlim=(0, 1.4),
        ylim=(-0.015, 0.52),
    )
    ax.legend(loc="upper left", frameon=False)
    ax.grid(alpha=0.22)
    fig.text(0.5, 0.025,
             "Die beiden linken Rangkarten sind identisch; "
             "unter den drei gezeichneten Profilen bleibt das ungerade Gewicht "
             "nur bei γ=(0,1,0) exakt null.",
             ha="center", fontsize=9, color="#34424b")
    fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.12)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _render_n5_comparison(cases: list[RankMapCertificate], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11.8, 5.6))
    grid = fig.add_gridspec(1, 3, width_ratios=(1, 1, 0.055), wspace=0.33)
    axes = [fig.add_subplot(grid[0, col]) for col in (0, 1)]
    fig.suptitle("N=5: der feste Mittelsitz zeichnet ein Rangkreuz", fontsize=14)
    titles = (
        "nur Mitte γ=(0,0,1,0,0)",
        "einseitig gebrochen γ=(1,0,1,0,0)",
    )
    for ax, case, title in zip(axes, cases, titles):
        image = _draw_rank_map(ax, case.ranks, vmin=9, vmax=25, title=title)
    fig.colorbar(image, cax=fig.add_subplot(grid[0, 2]),
                 label="exakter Signalrang")
    fig.text(0.5, 0.04,
             "Rang 9: Präparation oder Messung am Reflexionsfixpunkt; "
             "Rang 23: andere Paare. Einseitige Rate: Rang 25 überall.",
             ha="center", fontsize=9, color="#34424b")
    fig.subplots_adjust(left=0.09, right=0.94, top=0.84, bottom=0.18)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def run_view_comparison(output_dir: Path) -> dict:
    """Write the finite certificate and two figures; no all-N law is inferred."""
    if not all(sp.isprime(prime) for prime in PRIMES):
        raise ArithmeticError("the modular lower-bound route requires primes")
    results = [certified_rank_map(n, gamma) for n, gamma in CASES]
    by_key = {(item.n, item.gamma): item for item in results}
    n3_cases = [by_key[(3, gamma)] for gamma in ((0, 1, 0), (1, 1, 1), (1, 1, 0), (0, 3, 0))]
    n5_cases = [by_key[(5, gamma)] for gamma in ((0, 0, 1, 0, 0), (1, 0, 1, 0, 0))]
    output_dir.mkdir(parents=True, exist_ok=True)
    n3_figure = output_dir / "n3_rank_and_parity.png"
    n5_figure = output_dir / "n5_rank_maps.png"
    _render_n3_comparison(n3_cases[:3], n3_figure)
    _render_n5_comparison(n5_cases, n5_figure)

    odd_observation = {}
    for case in n3_cases:
        label = CASES[(case.n, case.gamma)]
        odd_observation[label] = {
            "gamma": case.gamma,
            "derivatives_at_zero": [int(value) for value in odd_weight_derivatives(case.gamma, 4)],
            "weight_at_t_0_5": odd_weight_curve(case.gamma, (0.5,))[0],
        }
    summary = {
        "scope": "finite open Heisenberg chains N=3,4,5, J=1, integer local Z rates",
        "prime_moduli": PRIMES,
        "rank_cases": [
            {"label": CASES[(item.n, item.gamma)], **asdict(item)}
            for item in results
        ],
        "odd_weight_from_centre_n3": odd_observation,
        "centre_return_from_centre_n3": {
            CASES[(case.n, case.gamma)]: centre_return_signature(case.gamma)
            for case in n3_cases
        },
        "figures": [n3_figure.name, n5_figure.name],
    }
    (output_dir / "operator_pair_view_comparison.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).parent / "results" / "operator_pair_view_comparison",
    )
    args = parser.parse_args()
    summary = run_view_comparison(args.output_dir)
    print(f"Certified {len(summary['rank_cases'])} finite rank maps and wrote 2 figures to {args.output_dir}")


if __name__ == "__main__":
    main()
