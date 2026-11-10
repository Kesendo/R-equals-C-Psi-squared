"""Draw the fifth centred trace under loss and gain as flux times the area of the rate triangle (F137;
docs/proofs/PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN.md §4).

Three sites on a triangle of hopping bonds, t01 = t12 = 1 and t20 = e^{iΦ}, so that the loop
0 -> 1 -> 2 -> 0 has Im(t01 t12 t20) = sin Φ. Every site has a rate point (γ↓, γ↑). The closed form is

    Tr((L + c/2)^5) = 60 · 4^N · sin Φ · det[(1, 1, 1), γ↓, γ↑],

the determinant being twice the signed area of the triangle the three rate points span.

What is calculated: the trace, directly as the fifth power of the 64 × 64 generator (imported from
simulations/fifth_moment_flux_area.py), and the closed form beside it. What is a drawing choice: the
two fixed rate points, the circle the third one runs on, the flux values, the colours.

Checked exactly before anything is drawn (integer rates, Φ = π/2, where every entry is a Gaussian
integer): the direct trace equals the closed form at a rate point on either side of the line through
the other two and at one on it, where both are zero. Printed, not gated: the largest deviation
between the two along the float sweeps.

  python simulations/flux_area_picture.py            # PNG and GIF
  python simulations/flux_area_picture.py --no-gif

Writes visualizations/flux_area_fifth_trace_neon.png and .gif.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyArrowPatch, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fifth_moment_flux_area import generator, hamiltonian  # noqa: E402

OUT = os.path.join(HERE, "..", "visualizations")

BG, PANEL, TEXT, FAINT, TITLE = "#080b12", "#0c1018", "#7ef9ff", "#1d3a46", "#39ff14"
CYAN, PINK, AMBER, VIOLET = "#08f7fe", "#fe53bb", "#ffe700", "#b14aed"
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": TEXT, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT, "ytick.color": TEXT, "axes.edgecolor": TEXT,
                     "figure.facecolor": BG, "savefig.facecolor": BG, "axes.facecolor": PANEL})

N = 3
P0, P1 = (4.0, 1.0), (8.0, 3.0)            # the two fixed rate points (γ↓, γ↑)
CENTRE, RADIUS = (6.0, 3.2), 1.8           # the circle the third rate point runs on
PHI0 = np.pi / 2                           # the flux of the sweep over the rate point
ANGLE0 = np.pi / 2                         # where the third point stands in the sweep over the flux


def point(angle):
    return CENTRE[0] + RADIUS * np.cos(angle), CENTRE[1] + RADIUS * np.sin(angle)


def det(p2):
    """det[(1, 1, 1), γ↓, γ↑]: twice the signed area of the rate triangle 0, 1, 2."""
    (x0, y0), (x1, y1), (x2, y2) = P0, P1, p2
    return (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)


def closed(p2, phi):
    return 60 * 4 ** N * np.sin(phi) * det(p2)


def direct(p2, phase):
    """Tr((L + c/2)^5) from the generator M = 2L + c; `phase` is t20, the amplitude of the hop 2 -> 0."""
    down = np.array([P0[0], P1[0], p2[0]]); up = np.array([P0[1], P1[1], p2[1]])
    H = hamiltonian(N, {(0, 1): 1, (1, 2): 1, (2, 0): phase}, {})
    M = generator(N, H, down + up, down - up)
    return np.trace(np.linalg.matrix_power(M, 5)) / 32


def checks():
    for p2 in ((6, 5), (6, 2), (6, 1)):                      # above the line, on it, below it
        down = [int(P0[0]), int(P1[0]), p2[0]]; up = [int(P0[1]), int(P1[1]), p2[1]]
        H = hamiltonian(N, {(0, 1): 1, (1, 2): 1, (2, 0): 1j}, {})
        M = generator(N, H, [d + u for d, u in zip(down, up)], [d - u for d, u in zip(down, up)])
        tr = np.trace(np.linalg.matrix_power(M, 5))
        want = 32 * 60 * 4 ** N * int(det(p2))
        assert tr.imag == 0.0 and tr.real == want, (p2, tr, want)
        print(f"rate point 2 at {p2}: Tr((L + c/2)^5) = {want // 32}, exact")
    assert det((6, 2)) == 0


def sweeps():
    angles = np.linspace(0, 2 * np.pi, 181)
    over_point = np.array([direct(point(a), np.exp(1j * PHI0)) for a in angles])
    fluxes = np.linspace(0, 2 * np.pi, 121)
    over_flux = np.array([direct(point(ANGLE0), np.exp(1j * f)) for f in fluxes])
    dev = max(np.abs(over_point - [closed(point(a), PHI0) for a in angles]).max(),
              np.abs(over_flux - closed(point(ANGLE0), fluxes)).max())
    print(f"largest |direct − closed form| along the two float sweeps: {dev:.3e} (values up to {np.abs(over_point).max():.0f})")
    return angles, over_point.real, fluxes, over_flux.real


def rate_plane(ax, p2, title=None, small=False):
    d = det(p2)
    col = CYAN if d > 1e-9 else PINK if d < -1e-9 else AMBER
    for m in (0.15, 0.3, 0.5, 0.8, 1.3):                      # lines of one ratio γ↑/γ↓, one temperature
        ax.plot([0, 12], [0, 12 * m], color=FAINT, lw=0.6, zorder=1)
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(CENTRE[0] + RADIUS * np.cos(t), CENTRE[1] + RADIUS * np.sin(t), color=VIOLET, lw=0.8, ls=(0, (2, 3)), zorder=2)
    x = np.array([0, 12.0]); slope = (P1[1] - P0[1]) / (P1[0] - P0[0])
    ax.plot(x, P0[1] + slope * (x - P0[0]), color=AMBER, lw=0.8, ls=(0, (5, 4)), zorder=2)
    ax.add_patch(Polygon([P0, P1, p2], closed=True, facecolor=col, alpha=0.28, edgecolor=col, lw=1.6, zorder=3))
    for i, p in enumerate((P0, P1, p2)):
        ax.scatter(*p, s=46 if small else 60, color=TEXT if i < 2 else AMBER, edgecolors=BG, zorder=5)
        ax.annotate(str(i), p, textcoords="offset points", xytext=(7, -12 if i < 2 else 7), fontsize=10, color=TEXT)
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.set_aspect("equal")
    ax.set_xlabel(r"loss $\gamma_\downarrow$"); ax.set_ylabel(r"gain $\gamma_\uparrow$")
    if title:
        ax.set_title(title, color=col, fontsize=10.5, loc="left")
    return col


def bonds(ax, phi_label=r"$\Phi$"):
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.87]])
    for a, b, lab in ((0, 1, r"$t_{01} = 1$"), (1, 2, r"$t_{12} = 1$"), (2, 0, r"$t_{20} = e^{i\Phi}$")):
        ax.add_patch(FancyArrowPatch(pts[a], pts[b], arrowstyle="-|>", mutation_scale=14, color=CYAN, lw=1.8,
                                     shrinkA=12, shrinkB=12, zorder=3))
        mid = (pts[a] + pts[b]) / 2; out = mid - pts.mean(axis=0)
        ax.text(*(mid + 0.75 * out), lab, ha="center", va="center", fontsize=10.5, color=TEXT)
    ax.scatter(pts[:, 0], pts[:, 1], s=330, color=PANEL, edgecolors=TEXT, linewidths=1.4, zorder=4)
    for i, p in enumerate(pts):
        ax.text(*p, str(i), ha="center", va="center", fontsize=11, color=TEXT, zorder=5)
    c = pts.mean(axis=0)
    ax.add_patch(FancyArrowPatch(c + (0.11, -0.02), c + (0.02, 0.11), connectionstyle="arc3,rad=1.6", arrowstyle="-|>",
                                 mutation_scale=10, color=AMBER, lw=1.2, zorder=3))
    ax.text(c[0], c[1] - 0.005, phi_label, ha="center", va="center", fontsize=12, color=AMBER)
    ax.set_xlim(-0.35, 1.35); ax.set_ylim(-0.25, 1.1); ax.set_aspect("equal"); ax.set_axis_off()


def figure(angles, over_point, fluxes, over_flux):
    fig = plt.figure(figsize=(14, 8.2))
    gs = fig.add_gridspec(2, 3, left=0.05, right=0.98, top=0.9, bottom=0.08, wspace=0.3, hspace=0.46)
    fig.text(0.05, 0.955, r"The fifth centred trace under loss and gain:  flux $\times$ the area of the rate triangle", color=TITLE, fontsize=14)
    for i, (a, word) in enumerate(((np.pi / 2, "one way round"), (None, "on the line of the other two: area 0"), (-np.pi / 2, "the other way round"))):
        if a is None:                                         # where the circle meets the line through 0 and 1
            p2 = point(brentq(lambda b: det(point(b)), -np.pi / 2, -np.pi))
        else:
            p2 = point(a)
        ax = fig.add_subplot(gs[0, i])
        rate_plane(ax, p2, rf"{word}:  $\mathrm{{Tr}}\,(L + c/2)^5$ = {closed(p2, PHI0) + 0.0:.0f}")
    ax = fig.add_subplot(gs[1, 0]); bonds(ax)
    ax.set_title("three sites, three hopping bonds, a flux through the loop", color=TITLE, fontsize=10.5, loc="left")
    fig.text(0.05, 0.475, r"faint lines: one ratio $\gamma_\uparrow/\gamma_\downarrow$, one temperature;   dashed: the line through rate points 0 and 1;   dotted: the circle rate point 2 runs on", fontsize=10)

    ax = fig.add_subplot(gs[1, 1])
    ax.axhline(0, color=FAINT, lw=0.8)
    ax.plot(np.degrees(angles), [closed(point(a), PHI0) for a in angles], color=AMBER, lw=1.4, label="closed form")
    ax.scatter(np.degrees(angles[::6]), over_point[::6], s=16, color=CYAN, zorder=4, label="fifth power of the generator")
    ax.set_xlabel("where rate point 2 stands on its circle (degrees)"); ax.set_ylabel(r"$\mathrm{Tr}\,(L + c/2)^5$")
    ax.set_title(r"moving one rate point, $\Phi = \pi/2$", color=TITLE, fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=9, labelcolor=TEXT); ax.grid(color=FAINT, lw=0.5)

    ax = fig.add_subplot(gs[1, 2])
    ax.axhline(0, color=FAINT, lw=0.8)
    ax.plot(fluxes / np.pi, closed(point(ANGLE0), fluxes), color=AMBER, lw=1.4, label="closed form")
    ax.scatter(fluxes[::5] / np.pi, over_flux[::5], s=16, color=CYAN, zorder=4, label="fifth power of the generator")
    ax.set_xlabel(r"flux $\Phi\,/\,\pi$"); ax.set_ylabel(r"$\mathrm{Tr}\,(L + c/2)^5$")
    ax.set_title("turning the flux, rate point 2 at the top of its circle", color=TITLE, fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=9, labelcolor=TEXT); ax.grid(color=FAINT, lw=0.5)
    path = os.path.join(OUT, "flux_area_fifth_trace_neon.png")
    fig.savefig(path, dpi=150); plt.close(fig)
    print("wrote", os.path.normpath(path))


def gif(angles, over_point, frames=90):
    fig = plt.figure(figsize=(9.2, 4.3))
    gs = fig.add_gridspec(1, 2, left=0.07, right=0.98, top=0.88, bottom=0.14, wspace=0.28)
    a1, a2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    curve = [closed(point(a), PHI0) for a in angles]
    lim = 1.12 * max(abs(min(curve)), abs(max(curve)))

    def frame(i):
        a = np.pi / 2 + 2 * np.pi * i / frames
        p2 = point(a)
        a1.clear(); a2.clear()
        col = rate_plane(a1, p2, small=True)
        a1.set_title(f"twice the signed area: {det(p2):+.2f}", color=col, fontsize=10.5, loc="left")
        a2.axhline(0, color=FAINT, lw=0.8)
        turn = np.linspace(np.pi / 2, a, max(2, int(180 * i / frames) + 2))
        a2.plot(np.degrees(turn - np.pi / 2), [closed(point(b), PHI0) for b in turn], color=AMBER, lw=1.4)
        a2.scatter([np.degrees(a - np.pi / 2)], [direct(p2, np.exp(1j * PHI0)).real], s=46, color=col, edgecolors=BG, zorder=5)
        a2.set_xlim(0, 360); a2.set_ylim(-lim, lim)
        a2.set_xlabel("how far rate point 2 has run round its circle (degrees)")
        a2.set_title(r"$\mathrm{Tr}\,(L + c/2)^5$ at flux $\pi/2$", color=TITLE, fontsize=10.5, loc="left")
        a2.grid(color=FAINT, lw=0.5)
        return []

    path = os.path.join(OUT, "flux_area_fifth_trace_neon.gif")
    FuncAnimation(fig, frame, frames=frames, blit=False).save(path, writer=PillowWriter(fps=12), dpi=90)
    plt.close(fig)
    print("wrote", os.path.normpath(path), f"({os.path.getsize(path) / 1e6:.2f} MB)")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    checks()
    angles, over_point, fluxes, over_flux = sweeps()
    figure(angles, over_point, fluxes, over_flux)
    if "--no-gif" not in sys.argv:
        gif(angles, over_point)
