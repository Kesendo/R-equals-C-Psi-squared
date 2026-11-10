"""Draw the letter cube's tetrahedron with the walk of Z₁ on it (docs/THE_ONE_SQUARE.md §7, §8).

What is calculated: the lattice points of the Pauli strings, (k_Z, k_X, k_Y), and how many strings
sit at each (framework.letter_cube); the operator Z₁ propagated under the Heisenberg chain with
Z-dephasing at one rate γ on every site, its centroid ⟨(k_Z, k_X, k_Y)⟩ weighted by the squared
string coefficients, and its Hilbert-Schmidt norm. What is a drawing choice: the camera (k_Z upward),
the colours, the dot sizes, N, J, γ and the time span.

Checked exactly before anything is drawn: the two edges are the strings with k_Z = 0 and with
k_Z = N, and R's letter flip carries one onto the other.
Printed, not gated: the largest |⟨k_X⟩ − ⟨k_Y⟩| along the walk of Z₁, the largest difference of its
shares at (k_Z, k_X, k_Y) and (k_Z, k_Y, k_X), and the largest deviation of ln‖ρ(t)‖² from
−4γ∫⟨k_Z⟩dt (trapezoid on the one grid dt).

  python simulations/one_square_tetrahedron.py                  # neon, PNG and GIF
  python simulations/one_square_tetrahedron.py --no-gif
  python simulations/one_square_tetrahedron.py --all-palettes   # also the okabe variant

Writes visualizations/one_square_tetrahedron_neon.png and .gif.
"""
from __future__ import annotations

import os
import sys
from itertools import product

import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from framework.letter_cube import cube_coords, cube_points  # noqa: E402

OUT = os.path.join(HERE, "..", "visualizations")

N, J, GAMMA = 4, 1.0, 0.05
T_MAX, DT = 40.0, 0.01
OFF_T = 8.0            # how long the off-plane walk of X₁Y₂ is drawn

PALETTES = {
    "neon": dict(bg="#080b12", panel="#0c1018", text="#7ef9ff", faint="#1d3a46", title="#39ff14",
                 quiet="#08f7fe", far="#fe53bb", walk="#ffe700", turn="#b14aed", off="#39ff14",
                 cmap=LinearSegmentedColormap.from_list("h", ["#08f7fe", "#4a6cf7", "#b14aed", "#fe53bb"])),
    "okabe": dict(bg="#f4f6f8", panel="#ffffff", text="#1b2530", faint="#c3ccd5", title="#1b2530",
                  quiet="#0072B2", far="#D55E00", walk="#000000", turn="#009E73", off="#CC79A7",
                  cmap=LinearSegmentedColormap.from_list("h", ["#0072B2", "#56B4E9", "#E69F00", "#D55E00"])),
}

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def string(label):
    m = np.eye(1, dtype=complex)
    for c in label:
        m = np.kron(m, PAULI[c])
    return m


def site(letter, l, n):
    return string("".join(letter if k == l else "I" for k in range(n)))


def walk(start="Z" + "I" * (N - 1)):
    """Propagate the string `start` (Z₁ by default) and return times, centroids (k_Z, k_X, k_Y), weights w, ln‖ρ‖² and the
    share of Σ|c_σ|² at each lattice point (the points in cube_points order)."""
    d = 2 ** N
    H = sum(J * site(P, l, N) @ site(P, l + 1, N) for l in range(N - 1) for P in "XYZ")
    eye = np.eye(d)
    # row-major vec: vec(AρB) = (A ⊗ Bᵀ) vec(ρ)
    L = -1j * (np.kron(H, eye) - np.kron(eye, H.T))
    for l in range(N):
        Z = site("Z", l, N)
        L += GAMMA * (np.kron(Z, Z.T) - np.eye(d * d))
    U = expm(L * DT)

    labels = ["".join(p) for p in product("IXYZ", repeat=N)]
    coords = np.array([cube_coords(s) for s in labels], dtype=float)
    weight = coords.sum(axis=1) / 2
    read = np.array([string(s).conj().reshape(-1) for s in labels]) / d     # c_σ = Tr(σρ)/2^N

    v = string(start).reshape(-1)
    steps = int(round(T_MAX / DT))
    t = np.arange(steps + 1) * DT
    cen = np.empty((steps + 1, 3)); w = np.empty(steps + 1); lognorm = np.empty(steps + 1)
    index = {k: i for i, k in enumerate(cube_points(N))}
    member = np.zeros((len(labels), len(index)))
    for j, s in enumerate(labels):
        member[j, index[cube_coords(s)]] = 1.0
    occ = np.empty((steps + 1, len(index)))
    for i in range(steps + 1):
        p = np.abs(read @ v) ** 2
        cen[i] = p @ coords / p.sum()
        w[i] = p @ weight / p.sum()
        occ[i] = p @ member / p.sum()
        lognorm[i] = np.log(np.vdot(v, v).real)
        v = U @ v
    return t, cen, w, lognorm, occ


def checks(t, cen, lognorm):
    pts = cube_points(N)
    labels = ["".join(p) for p in product("IXYZ", repeat=N)]
    quiet = {s for s in labels if set(s) <= set("IZ")}
    far = {s for s in labels if set(s) <= set("XY")}
    assert quiet == {s for s in labels if cube_coords(s)[0] == 0}
    assert far == {s for s in labels if cube_coords(s)[0] == N}
    flip = str.maketrans("IXYZ", "XIZY")
    assert {s.translate(flip) for s in quiet} == far
    for s in labels:
        kz, kx, ky = cube_coords(s)
        assert cube_coords(s.translate(flip)) == (N - kz, kx, N - ky)
    assert sum(pts.values()) == 4 ** N
    print(f"lattice points {len(pts)}, strings {sum(pts.values())}; edges and R's flip: exact")
    print(f"max |⟨k_X⟩ − ⟨k_Y⟩| along the walk: {np.abs(cen[:, 1] - cen[:, 2]).max():.3e}")
    budget = -4 * GAMMA * np.concatenate([[0], np.cumsum((cen[1:, 0] + cen[:-1, 0]) / 2) * DT])
    dev = np.abs((lognorm - lognorm[0]) - budget).max()
    print(f"max |ln‖ρ‖² − ln‖ρ₀‖² + 4γ∫⟨k_Z⟩dt| = {dev:.3e} at dt = {DT} (trapezoid)")
    return budget


def xyz(k):
    """Plot coordinates: k_X, k_Y on the floor, k_Z upward."""
    k = np.asarray(k, dtype=float)
    return k[..., 1], k[..., 2], k[..., 0]


def draw_tetra(ax, pal, cen, head=None, labels=True, off=None, occ=None):
    pts = cube_points(N)
    V = {"I": (0, 0, 0), "X": (N, 0, N), "Y": (N, N, 0), "Z": (0, N, N)}   # (k_Z, k_X, k_Y)
    for a, b in [("I", "X"), ("I", "Y"), ("Z", "X"), ("Z", "Y")]:
        ax.plot(*xyz([V[a], V[b]]), color=pal["faint"], lw=0.9, zorder=1)
    ax.plot(*xyz([V["I"], V["Z"]]), color=pal["quiet"], lw=2.6, zorder=2)
    ax.plot(*xyz([V["X"], V["Y"]]), color=pal["far"], lw=2.6, zorder=2)
    mid = [(N / 2, N / 2, N / 2 - 0.0)]
    ax.plot(*xyz([(N / 2, -0.35, N / 2), (N / 2, N + 0.35, N / 2)]), color=pal["turn"], lw=1.0, ls=(0, (4, 3)), zorder=1)

    k = np.array(list(pts.keys()), dtype=float)
    cnt = np.array(list(pts.values()), dtype=float)
    if occ is None:
        ax.scatter(*xyz(k), s=14 + 16 * np.sqrt(cnt), c=k[:, 0] / N, cmap=pal["cmap"], vmin=0, vmax=1,
                   edgecolors=pal["bg"], linewidths=0.4, alpha=0.95, depthshade=False, zorder=3)
    else:                                                  # disc area: the share of the operator at the point
        ax.scatter(*xyz(k), s=5, color=pal["faint"], depthshade=False, zorder=2)
        ax.scatter(*xyz(k), s=1500 * occ, c=k[:, 0] / N, cmap=pal["cmap"], vmin=0, vmax=1,
                   edgecolors=pal["bg"], linewidths=0.4, alpha=0.9, depthshade=False, zorder=3)
    ax.scatter(*xyz(np.array(mid)), s=40, facecolors="none", edgecolors=pal["turn"], linewidths=1.0, zorder=3)

    plane = np.array([V["I"], V["Z"], (N, N / 2, N / 2)], dtype=float)      # the plane k_X = k_Y inside the tetrahedron
    px, py, pz = xyz(plane)
    ax.add_collection3d(Poly3DCollection([list(zip(px, py, pz))], facecolor=pal["walk"], alpha=0.07, edgecolor=pal["walk"], linewidths=0.5))
    ax.plot(*xyz(cen), color=pal["walk"], lw=1.6, zorder=5)
    if off is not None:                                    # a start the quarter-turn does not fix
        ax.plot(*xyz(off), color=pal["off"], lw=1.4, zorder=5)
        ax.scatter(*xyz(off[:1]), s=26, color=pal["off"], edgecolors=pal["bg"], linewidths=0.5, depthshade=False, zorder=6)
    if head is not None:
        ax.scatter(*xyz(cen[head:head + 1]), s=70, facecolors="none", edgecolors=pal["walk"], linewidths=1.6,
                   depthshade=False, zorder=6)
    if labels:
        off = {"I": (-0.1, -0.35, -0.35), "Z": (-0.1, N + 0.3, N + 0.3), "X": (N + 0.35, -0.3, N + 0.3), "Y": (N + 0.35, N + 0.3, -0.3)}
        for name, o in off.items():
            x, y, z = xyz(np.array(o))
            ax.text(x, y, z, f"{name}$^{{{N}}}$", color=pal["text"], fontsize=11, ha="center", va="center")
    ax.set_xlim(0, N); ax.set_ylim(0, N); ax.set_zlim(0, N)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()


def figure(name, pal, t, cen, w, lognorm, budget, off):
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": pal["text"], "axes.labelcolor": pal["text"],
                         "xtick.color": pal["text"], "ytick.color": pal["text"], "axes.edgecolor": pal["text"]})
    fig = plt.figure(figsize=(13.5, 7.2), facecolor=pal["bg"])
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], left=0.0, right=0.96, top=0.9, bottom=0.12, wspace=0.12, hspace=0.34)

    ax = fig.add_subplot(gs[:, 0], projection="3d", facecolor=pal["bg"])
    draw_tetra(ax, pal, cen, off=off)
    ax.view_init(elev=14, azim=-62)
    fig.text(0.03, 0.94, f"The Pauli strings of N = {N} sites on the tetrahedron", color=pal["title"], fontsize=13)
    fig.text(0.03, 0.905, r"dot size: how many strings sit at the point;  $k_Z$ upward, the coordinate the rate $2\gamma k_Z$ reads", color=pal["text"], fontsize=10)
    fig.text(0.03, 0.145, f"green line:  the centroid of X₁Y₂ up to t = {OFF_T:g}, which starts on the plane and leaves it sideways", color=pal["off"], fontsize=10)
    fig.text(0.03, 0.115, r"quiet edge  $k_Z = 0$:  the populations $\{I, Z\}^N$, rate 0", color=pal["quiet"], fontsize=10)
    fig.text(0.03, 0.085, r"far edge  $k_Z = N$:  the strings of X and Y, rate $2N\gamma$", color=pal["far"], fontsize=10)
    fig.text(0.03, 0.055, "dashed:  the axis of R's half-turn, which carries one edge onto the other", color=pal["turn"], fontsize=10)
    fig.text(0.03, 0.025, "line:  the centroid of Z₁ under the Heisenberg chain, inside the plane $k_X = k_Y$ (shaded)", color=pal["walk"], fontsize=10)

    a1 = fig.add_subplot(gs[0, 1], facecolor=pal["panel"])
    a1.plot(w, cen[:, 0], color=pal["walk"], lw=1.3)
    a1.scatter([w[0]], [cen[0, 0]], color=pal["quiet"], s=40, zorder=5, label="start: Z₁")
    a1.scatter([w[-1]], [cen[-1, 0]], facecolors="none", edgecolors=pal["far"], linewidths=1.6, s=150, zorder=5, label=f"t = {T_MAX:g}")
    a1.axhline(0, color=pal["quiet"], lw=1.0, alpha=0.7)
    a1.set_xlabel(r"size $\langle w\rangle$, the number of letters other than I")
    a1.set_ylabel(r"height $\langle k_Z\rangle$")
    a1.set_title(r"The walk in its own plane $k_X = k_Y$", color=pal["title"], fontsize=11, loc="left")
    a1.legend(frameon=False, fontsize=9, labelcolor=pal["text"])
    a1.grid(color=pal["faint"], lw=0.5)

    a2 = fig.add_subplot(gs[1, 1], facecolor=pal["panel"])
    a2.plot(t, lognorm - lognorm[0], color=pal["far"], lw=3.2, alpha=0.45, label=r"$\ln\|\rho(t)\|^2 - \ln\|\rho(0)\|^2$")
    a2.plot(t, budget, color=pal["walk"], lw=1.1, label=r"$-4\gamma\int\langle k_Z\rangle\,dt$")
    a2.set_xlabel(rf"t   (J = {J:g}, $\gamma$ = {GAMMA:g})")
    a2.set_ylabel("what the operator has lost")
    a2.set_title("The budget: every loss is paid above the quiet edge", color=pal["title"], fontsize=11, loc="left")
    a2.legend(frameon=False, fontsize=9, labelcolor=pal["text"])
    a2.grid(color=pal["faint"], lw=0.5)

    path = os.path.join(OUT, f"one_square_tetrahedron_{name}.png")
    fig.savefig(path, dpi=150, facecolor=pal["bg"])
    plt.close(fig)
    print("wrote", os.path.normpath(path))


def gif(name, pal, cen, occ, t_end=20.0, frames=125):
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": pal["text"]})
    fig = plt.figure(figsize=(5.2, 5.2), facecolor=pal["bg"])
    ax = fig.add_axes([0, 0, 1, 1], projection="3d", facecolor=pal["bg"])
    stride = int(round(t_end / DT)) // frames

    def frame(i):
        ax.clear()
        draw_tetra(ax, pal, cen[: i * stride + 1], head=i * stride, occ=occ[i * stride])
        ax.view_init(elev=14, azim=-62 + 360 * i / frames)
        ax.text2D(0.04, 0.95, f"t = {i * stride * DT:4.1f}", transform=ax.transAxes, color=pal["text"], fontsize=11)
        ax.text2D(0.04, 0.04, "disc area: the share of Z₁ at the point;  ring: its centroid", transform=ax.transAxes, color=pal["text"], fontsize=8.5)
        return []

    path = os.path.join(OUT, f"one_square_tetrahedron_{name}.gif")
    FuncAnimation(fig, frame, frames=frames + 1, blit=False).save(path, writer=PillowWriter(fps=12), dpi=100)
    plt.close(fig)
    print("wrote", os.path.normpath(path), f"({os.path.getsize(path) / 1e6:.2f} MB)")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    t, cen, w, lognorm, occ = walk()
    budget = checks(t, cen, lognorm)
    keys = list(cube_points(N))
    swap = [keys.index((kz, ky, kx)) for kz, kx, ky in keys]
    print(f"max |share at (k_Z, k_X, k_Y) − share at (k_Z, k_Y, k_X)| for Z₁: {np.abs(occ - occ[:, swap]).max():.3e}")
    off = walk("XY" + "I" * (N - 2))[1][: int(round(OFF_T / DT)) + 1]
    assert off[0, 1] == off[0, 2]
    print(f"X₁Y₂ starts on the plane; max |⟨k_X⟩ − ⟨k_Y⟩| along its walk: {np.abs(off[:, 1] - off[:, 2]).max():.4f}")
    names = list(PALETTES) if "--all-palettes" in sys.argv else ["neon"]
    for name in names:
        figure(name, PALETTES[name], t, cen, w, lognorm, budget, off)
        if "--no-gif" not in sys.argv:
            gif(name, PALETTES[name], cen, occ)
