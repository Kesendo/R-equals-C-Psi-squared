"""Finite operator-pair flow and preparation/readout atlas pilot."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import sympy as sp
from scipy.linalg import expm


@dataclass(frozen=True)
class PairFlow:
    states: tuple[int, ...]
    flat_indices: tuple[int, ...]
    generator: np.ndarray


@dataclass(frozen=True)
class ReadoutRank:
    dimension: int
    rank: int
    moments: tuple[int, ...]


def one_excitation_basis(n: int) -> tuple[int, ...]:
    """Computational states with one occupied site, site 0 at the leftmost bit."""
    if n < 1:
        raise ValueError("n must be positive")
    return tuple(1 << (n - 1 - site) for site in range(n))


def single_excitation_flow(h: np.ndarray, gamma_per_site: Sequence[float]) -> PairFlow:
    """Act on E_ab=|1_a><1_b| by signed H edges and local Z costs.

    This is a closed sector only if H preserves the one-excitation subspace.
    The pair basis is row-major: (a,b) has index a*n+b.
    """
    gamma = tuple(float(rate) for rate in gamma_per_site)
    n = len(gamma)
    states = one_excitation_basis(n)
    d = 1 << n
    h = np.asarray(h, dtype=complex)
    if h.shape != (d, d):
        raise ValueError(f"H must have shape {(d, d)}")
    if not np.all(np.isfinite(h)) or not np.all(np.isfinite(gamma)):
        raise ValueError("H and gamma must be finite")
    if not np.allclose(h, h.conj().T, atol=1e-12, rtol=0):
        raise ValueError("H must be Hermitian")
    outside = [state for state in range(d) if state not in states]
    leakage = np.max(np.abs(h[np.ix_(outside, states)])) if outside else 0.0
    if leakage > 1e-12:
        raise ValueError(f"one-excitation sector is not invariant under H (leakage {leakage:g})")

    generator = np.zeros((n * n, n * n), dtype=complex)
    for a, state_a in enumerate(states):
        for b, state_b in enumerate(states):
            source = a * n + b
            difference = state_a ^ state_b
            for site, rate in enumerate(gamma):
                if difference & (1 << (n - 1 - site)):
                    generator[source, source] -= 2 * rate
            for c, state_c in enumerate(states):
                # -i H|a><b| acts on the first index; +i |a><b|H on the second.
                generator[c * n + b, source] -= 1j * h[state_c, state_a]
                generator[a * n + c, source] += 1j * h[state_b, state_c]

    flat = tuple(a * d + b for a in states for b in states)
    return PairFlow(states=states, flat_indices=flat, generator=generator)


def exact_heisenberg_se_h(n: int, coupling: int) -> sp.Matrix:
    """Integer single-excitation H for J sum(XX+YY+ZZ) on an open chain."""
    if n < 2 or not isinstance(coupling, int):
        raise ValueError("n >= 2 and integer coupling required")
    h = sp.zeros(n)
    bonds = n - 1
    for site in range(n):
        degree = int(site > 0) + int(site < n - 1)
        h[site, site] = coupling * (bonds - 2 * degree)
    for site in range(n - 1):
        h[site, site + 1] = 2 * coupling
        h[site + 1, site] = 2 * coupling
    return h


def exact_readout_rank(
    h: sp.Matrix,
    gamma_per_site: Sequence[int],
    *,
    preparation: int,
    readout: int | str,
) -> ReadoutRank:
    """Exact moment-Hankel rank of one population signal in an integer SE model.

    The real Hermitian basis has N populations, then Re(rho_ab), then
    Im(rho_ab), a<b. It is N²-dimensional over the reals. Only integer
    H and gamma are accepted, so the reported rank has no float threshold.
    """
    n = h.rows
    if h.cols != n or len(gamma_per_site) != n:
        raise ValueError("H must be square and gamma must have one entry per site")
    if any(not isinstance(rate, (int, sp.Integer)) for rate in gamma_per_site):
        raise ValueError("exact readout rank requires integer gamma")
    if any(not value.is_Integer for value in h):
        raise ValueError("exact readout rank requires integer H")
    if h != h.T:
        raise ValueError("exact readout rank requires real symmetric H")
    if not 0 <= preparation < n:
        raise ValueError("preparation site out of range")
    if readout != "trace" and not (isinstance(readout, int) and 0 <= readout < n):
        raise ValueError("readout must be a site or 'trace'")

    gamma = tuple(sp.Integer(rate) for rate in gamma_per_site)
    labels = [("population", a, a) for a in range(n)]
    labels += [("real", a, b) for a in range(n) for b in range(a + 1, n)]
    labels += [("imag", a, b) for a in range(n) for b in range(a + 1, n)]
    basis = []
    for kind, a, b in labels:
        cell = sp.zeros(n)
        if kind == "population":
            cell[a, a] = 1
        elif kind == "real":
            cell[a, b] = cell[b, a] = 1
        else:
            cell[a, b] = sp.I
            cell[b, a] = -sp.I
        basis.append(cell)

    dimension = n * n
    real_generator = sp.zeros(dimension)
    for source, cell in enumerate(basis):
        image = -sp.I * (h * cell - cell * h)
        for a in range(n):
            for b in range(n):
                if a != b:
                    image[a, b] -= 2 * (gamma[a] + gamma[b]) * cell[a, b]
        for target, (kind, a, b) in enumerate(labels):
            value = image[a, b]
            real_generator[target, source] = (
                value if kind == "population" else
                sp.re(value) if kind == "real" else sp.im(value)
            )

    initial = sp.zeros(dimension, 1)
    initial[preparation] = 1
    measurement = sp.zeros(1, dimension)
    if readout == "trace":
        for site in range(n):
            measurement[0, site] = 1
    else:
        measurement[0, readout] = 1
    moments = []
    state = initial
    for _ in range(2 * dimension - 1):
        moments.append(int((measurement * state)[0]))
        state = real_generator * state
    hankel = sp.Matrix(dimension, dimension, lambda i, j: moments[i + j])
    return ReadoutRank(dimension, hankel.rank(), tuple(moments))


def exact_readout_ode(result: ReadoutRank) -> tuple[sp.Rational, ...]:
    """Return a_j in y^(r) + sum_{j<r} a_j y^(j) = 0 for this signal.

    The leading r-by-r Hankel minor must be nonsingular. This is true for
    the finite pilot; other moment sequences may require a shifted minor.
    Every available later moment is checked against the recurrence.
    """
    r = result.rank
    if r < 1:
        raise ValueError("zero signal has no positive-order ODE")
    moments = result.moments
    leading = sp.Matrix(r, r, lambda i, j: moments[i + j])
    if leading.det() == 0:
        raise ValueError("leading Hankel minor is singular")
    rhs = sp.Matrix([-moments[r + i] for i in range(r)])
    coefficients = tuple(leading.LUsolve(rhs))
    for start in range(len(moments) - r):
        residual = moments[start + r] + sum(
            coefficients[j] * moments[start + j] for j in range(r)
        )
        if residual != 0:
            raise ArithmeticError("moment recurrence fails exact validation")
    return coefficients


def _render_n3_atlas(flow: PairFlow, end_rank: int, centre_rank: int, path: Path) -> None:
    """Draw signed pair moves and two input-output traces for the fixed N=3 pilot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle

    n = 3
    fig, (ax, ax_signal) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    fig.suptitle("Operatorpaar-Fluss: Heisenberg-Kette, N=3, J=1, γ=(0,1,0)", fontsize=13)

    # Every source-to-target H move changes exactly one pair index. Draw each
    # two-way edge once, retaining its distinct left/right complex sign.
    for a in range(n):
        for b in range(n):
            for c in range(a + 1, n):
                if abs(flow.generator[(c * n + b), (a * n + b)]) > 1e-12:
                    ax.plot([b, b], [-a, -c], color="#2864a6", lw=2.0, zorder=1)
            for c in range(b + 1, n):
                if abs(flow.generator[(a * n + c), (a * n + b)]) > 1e-12:
                    ax.plot([b, c], [-a, -a], color="#cb632f", lw=2.0, zorder=1)

    for a in range(n):
        for b in range(n):
            diagonal = flow.generator[a * n + b, a * n + b]
            cost = -diagonal.real
            phase = diagonal.imag
            face = "#f2b778" if cost > 0 else "#e9f2e9"
            ax.add_patch(Circle((b, -a), 0.34, facecolor=face,
                                edgecolor="#34424b", linewidth=1.25, zorder=3))
            ax.text(b, -a + 0.14, f"({a},{b})", ha="center", va="center",
                    fontsize=10, fontweight="bold", zorder=4)
            ax.text(b, -a - 0.015, f"D={diagonal.real:g}", ha="center", va="center",
                    fontsize=8.2, zorder=4)
            phase_label = f"{phase:+g}" if abs(phase) > 1e-12 else "0"
            ax.text(b, -a - 0.17, f"ω={phase_label}", ha="center", va="center",
                    fontsize=8.2, zorder=4)
    ax.add_patch(Circle((0, 0), 0.40, fill=False, edgecolor="#6b3c91", lw=3, zorder=5))
    ax.add_patch(Circle((1, -1), 0.40, fill=False, edgecolor="#158069", lw=3, zorder=5))
    ax.set(xlim=(-0.65, 2.65), ylim=(-2.65, 0.65), xlabel="zweiter Index b",
           ylabel="erster Index a", title="(1,1)-Block: 9 Paarknoten")
    ax.set_xticks(range(n))
    ax.set_yticks([-a for a in range(n)], labels=[str(a) for a in range(n)])
    ax.set_aspect("equal")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.legend(handles=[
        Line2D([0], [0], color="#2864a6", lw=2, label="erster Index: −iH"),
        Line2D([0], [0], color="#cb632f", lw=2, label="zweiter Index: +iH"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#f2b778",
               markeredgecolor="#34424b", markersize=11, label="Dephasierung D<0"),
    ], loc="lower center", bbox_to_anchor=(0.5, 0.085), ncol=3,
       frameon=False, fontsize=9)

    times = np.linspace(0, 2.5, 151)
    for site, rank, color, name in (
        (0, end_rank, "#6b3c91", "Rückkehr an Site 0"),
        (1, centre_rank, "#158069", "Rückkehr an Site 1"),
    ):
        initial = np.zeros(n * n, dtype=complex)
        initial[site * n + site] = 1
        signal = [(expm(t * flow.generator) @ initial)[site * n + site].real
                  for t in times]
        ax_signal.plot(times, signal, color=color, lw=2.3,
                       label=f"{name}: exakter Hankel-Rang {rank}")
    ax_signal.set(xlabel="Zeit t", ylabel="gemessene Population",
                  title="Ein Generator, verschiedene Lesefenster", ylim=(-0.03, 1.08))
    ax_signal.legend(loc="upper right", frameon=False, fontsize=9)
    ax_signal.grid(alpha=0.2)
    fig.text(0.5, 0.035,
             "Knotendiagonale D+iω mit ω=−(Eₐ−Eᵦ) aus ZZ; "
             "der Rang beschreibt nur das jeweils gezeigte skalare Signal.",
             ha="center", fontsize=9, color="#35434a")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.23, wspace=0.35)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def run_pilot(output_dir: Path) -> dict:
    """Reproduce finite checks and save the figure/JSON, without claiming a law."""
    from framework.chain_system import ChainSystem
    from framework.lindblad import lindbladian_z_dephasing

    comparisons = []
    for n in (3, 4, 5):
        h = ChainSystem(n, J=1.0, H_type="heisenberg").H
        gamma = [round(0.13 + 0.07 * site, 2) for site in range(n)]
        pair = single_excitation_flow(h, gamma)
        full = lindbladian_z_dephasing(h, gamma)
        flat = list(pair.flat_indices)
        comparisons.append({
            "N": n,
            "site_gamma": gamma,
            "block_dimension": n * n,
            "max_entry_difference": float(np.max(np.abs(
                pair.generator - full[np.ix_(flat, flat)]
            ))),
        })

    n = 3
    h = ChainSystem(n, J=1.0, H_type="heisenberg").H
    gamma = [0, 1, 0]
    pair = single_excitation_flow(h, gamma)
    exact_h = exact_heisenberg_se_h(n, 1)
    end = exact_readout_rank(exact_h, gamma, preparation=0, readout=0)
    centre = exact_readout_rank(exact_h, gamma, preparation=1, readout=1)
    trace = exact_readout_rank(exact_h, gamma, preparation=0, readout="trace")
    end_ode = exact_readout_ode(end)
    centre_ode = exact_readout_ode(centre)

    # A non-real preparation and an off-diagonal Hermitian observable make a
    # wrong row/column vec convention visible in the finite-time comparison.
    d = 1 << n
    psi = np.zeros(d, dtype=complex)
    psi[list(pair.states)] = np.array([1, 1j, 0.5 + 0.25j])
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())
    observable = np.zeros((d, d), dtype=complex)
    a, b, c = pair.states
    observable[a, b], observable[b, a], observable[c, c] = 1j, -1j, 0.5
    t = 0.37
    full = lindbladian_z_dephasing(h, gamma)
    full_rho = (expm(t * full) @ rho.flatten()).reshape(d, d)
    pair_rho = (expm(t * pair.generator) @ rho[np.ix_(pair.states, pair.states)].flatten()).reshape(n, n)
    full_signal = np.trace(observable @ full_rho)
    pair_signal = np.trace(observable[np.ix_(pair.states, pair.states)] @ pair_rho)

    output_dir.mkdir(parents=True, exist_ok=True)
    figure_path = output_dir / "operator_pair_flow_atlas_n3.png"
    _render_n3_atlas(pair, end.rank, centre.rank, figure_path)
    summary = {
        "scope": "finite N=3..5 Heisenberg open chain, local Z dephasing; rank case N=3 J=1 gamma=(0,1,0)",
        "generator_comparisons": comparisons,
        "readout_ranks": {"site0_return": end.rank, "site1_return": centre.rank,
                          "trace_from_site0": trace.rank,
                          "block_dimension": end.dimension},
        "ode_coefficients_ascending": {
            "site0_return": [int(value) for value in end_ode],
            "site1_return": [int(value) for value in centre_ode],
        },
        "complex_signal_at_t_0_37": {
            "full_real": float(full_signal.real),
            "pair_real": float(pair_signal.real),
            "absolute_difference": float(abs(full_signal - pair_signal)),
        },
        "figure": figure_path.name,
    }
    (output_dir / "operator_pair_flow_atlas.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).parent / "results" / "operator_pair_flow_atlas")
    args = parser.parse_args()
    summary = run_pilot(args.output_dir)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
