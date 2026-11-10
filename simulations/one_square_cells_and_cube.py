"""Draw the cell square with its eight moves (docs/THE_ONE_SQUARE.md §2) and the one-site letter
cube with its three half-turns (§7).

What is calculated: where each of the eight moves sends a cell |a⟩⟨b| and the disagreement count
popcount(a ⊕ b) it lands on; the corners (k_Z, k_X, k_Y) of the four letters and the pair sums of
the three pairings. What is a drawing choice: N = 3 for the square, the camera, the colours.

Checked before anything is drawn, all exact (integers):
  §2's table: the first four moves keep a ⊕ b and the other four complement it; the fixed cells are
  all, none, the diagonal, the anti-diagonal, and none for each paying move; the eight are distinct;
  §7: the letters fill the even tetrahedron, the half-turn about the k_P axis pairs I with P, and
  its two pair sums differ by 4γ_P for every rate triple in {0, …, 3}³.

  python simulations/one_square_cells_and_cube.py

Writes visualizations/one_square_cell_square_neon.png and visualizations/one_square_letter_cube_neon.png.
"""
from __future__ import annotations

import os
import sys
from itertools import product

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from framework.letter_cube import cube_coords  # noqa: E402

OUT = os.path.join(HERE, "..", "visualizations")

BG, PANEL, TEXT, FAINT, TITLE = "#080b12", "#0c1018", "#7ef9ff", "#1d3a46", "#39ff14"
CYAN, PINK, AMBER, VIOLET = "#08f7fe", "#fe53bb", "#ffe700", "#b14aed"
HEIGHT = LinearSegmentedColormap.from_list("h", ["#0c1018", "#2b2f7a", "#8a4fd8", "#eadcff"])
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": TEXT, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT, "ytick.color": TEXT, "axes.edgecolor": TEXT,
                     "figure.facecolor": BG, "savefig.facecolor": BG})

N = 3
D = 2 ** N
FULL = D - 1

# each move as (a, b) ↦ (a', b'); the bar is the complement of every bit
MOVES = [
    ("1", lambda a, b: (a, b)),
    (r"$\mathcal{F}$", lambda a, b: (a ^ FULL, b ^ FULL)),
    ("D", lambda a, b: (b, a)),
    (r"$\mathcal{F}$D", lambda a, b: (b ^ FULL, a ^ FULL)),
    ("R", lambda a, b: (a, b ^ FULL)),
    (r"$\mathcal{F}$R", lambda a, b: (a ^ FULL, b)),
    (r"$\Pi_Z = R\cdot D$", lambda a, b: (b, a ^ FULL)),
    (r"$\Pi_Y$", lambda a, b: (b ^ FULL, a)),
]
FIXED = ["all", "none", "diagonal", "anti-diagonal", "none", "none", "none", "none"]
PROBE = (0b001, 0b011)          # one cell followed through the eight moves


def pop(x):
    return bin(x).count("1")


def check_square():
    cells = list(product(range(D), repeat=2))
    images = []
    for i, (name, f) in enumerate(MOVES):
        img = {c: f(*c) for c in cells}
        assert sorted(img.values()) == cells, name                       # a permutation
        for (a, b), (a2, b2) in img.items():
            assert (a2 ^ b2) == ((a ^ b) if i < 4 else (a ^ b) ^ FULL), name
        fixed = {c for c in cells if img[c] == c}
        want = {"all": set(cells), "none": set(), "diagonal": {(a, a) for a in range(D)},
                "anti-diagonal": {(a, a ^ FULL) for a in range(D)}}[FIXED[i]]
        assert fixed == want, name
        if i >= 4:
            assert {img[(a, a)] for a in range(D)} == {(a, a ^ FULL) for a in range(D)}, name
        images.append(tuple(img[c] for c in cells))
    assert len(set(images)) == 8
    print("cell square: eight distinct moves, disagreement kept / complemented, fixed cells: exact")


def check_cube():
    corners = {P: cube_coords(P) for P in "IXYZ"}
    assert set(corners.values()) == {c for c in product((0, 1), repeat=3) if sum(c) % 2 == 0}
    axis = {"Z": 0, "X": 1, "Y": 2}
    for P, ax in axis.items():
        turn = lambda k: tuple(v if i == ax else 1 - v for i, v in enumerate(k))
        assert turn(corners["I"]) == corners[P]
        others = [Q for Q in "XYZ" if Q != P]
        assert turn(corners[others[0]]) == corners[others[1]]
        for gz, gx, gy in product(range(4), repeat=3):
            rate = {Q: 2 * (gz * k[0] + gx * k[1] + gy * k[2]) for Q, k in corners.items()}
            g = {"Z": gz, "X": gx, "Y": gy}[P]
            assert rate["I"] + rate[P] == 2 * (gz + gx + gy - g)
            assert rate[others[0]] + rate[others[1]] == 2 * (gz + gx + gy + g)
    print("letter cube: even tetrahedron, three half-turns, pair sums differ by 4γ_P: exact")
    return corners


def tile(ax, f, title, big=False):
    grid = np.empty((D, D))
    for a, b in product(range(D), repeat=2):
        a2, b2 = f(a, b)
        grid[a2, b2] = pop(a ^ b)                 # the cell carries its disagreement count with it
    ax.imshow(grid, cmap=HEIGHT, vmin=0, vmax=N, interpolation="nearest")
    for a in range(D):                            # where the diagonal and the anti-diagonal land
        for (src, col) in (((a, a), CYAN), ((a, a ^ FULL), PINK)):
            r, c = f(*src)
            ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False, ec=col, lw=2.2 if big else 1.5, zorder=4 if col == PINK else 3))
    r, c = f(*PROBE)
    ax.scatter([c], [r], marker="o", s=70 if big else 26, color=AMBER, edgecolors=BG, linewidths=0.6, zorder=5)
    ax.set_title(title, color=TITLE, fontsize=12 if big else 10)
    if big:
        lab = [format(i, f"0{N}b") for i in range(D)]
        ax.set_xticks(range(D), lab, fontsize=8); ax.set_yticks(range(D), lab, fontsize=8)
        ax.set_xlabel("bra b"); ax.set_ylabel("ket a")
        for a, b in product(range(D), repeat=2):
            k = pop(a ^ b)
            ax.text(b, a, str(k), ha="center", va="center", fontsize=8, color=TEXT if k < 2 else BG)
    else:
        ax.set_xticks([]); ax.set_yticks([])


def draw_square():
    fig = plt.figure(figsize=(14, 6.6))
    gs = fig.add_gridspec(2, 6, width_ratios=[2.25, 0.12, 1, 1, 1, 1], left=0.06, right=0.985, top=0.84, bottom=0.2, wspace=0.16, hspace=0.42)
    tile(fig.add_subplot(gs[:, 0]), MOVES[0][1], f"The cells |a⟩⟨b| of N = {N} sites", big=True)
    heads = [("the four that keep the two lines: the dephasing is unchanged", CYAN),
             (r"the four that exchange them: the dephasing becomes $-L_D - 2\Sigma\gamma$", PINK)]
    for i, (name, f) in enumerate(MOVES):
        a = fig.add_subplot(gs[i // 4, 2 + i % 4])
        tile(a, f, name)
        if i % 4 == 0:
            a.text(0, 1.2, heads[i // 4][0], color=heads[i // 4][1], fontsize=11, transform=a.transAxes)
    fig.text(0.06, 0.105, "number in a cell: at how many sites ket and bra disagree; the cell decays at 2γ times it (one rate γ on every site)", fontsize=10)
    fig.text(0.06, 0.07, "cyan frame: the diagonal a = b, the populations, and where a move carries it", color=CYAN, fontsize=10)
    fig.text(0.06, 0.035, "pink frame: the anti-diagonal b = ā, every site disagreeing, and where a move carries it", color=PINK, fontsize=10)
    fig.text(0.56, 0.035, "dot: the one cell |001⟩⟨011| followed through the eight moves", color=AMBER, fontsize=10)
    path = os.path.join(OUT, "one_square_cell_square_neon.png")
    fig.savefig(path, dpi=150); plt.close(fig)
    print("wrote", os.path.normpath(path))


def xyz(k):
    k = np.asarray(k, dtype=float)
    return k[..., 1], k[..., 2], k[..., 0]          # k_X, k_Y on the floor, k_Z upward


def cube(ax, corners, turn=None, big=False):
    for p, q in product(product((0, 1), repeat=3), repeat=2):
        if sum(abs(i - j) for i, j in zip(p, q)) == 1 and p < q:
            ax.plot(*xyz([p, q]), color=FAINT, lw=0.9)
    odd = np.array([c for c in product((0, 1), repeat=3) if sum(c) % 2])
    ax.scatter(*xyz(odd), s=70 if big else 34, facecolors="none", edgecolors=FAINT, linewidths=1.2, depthshade=False)
    letters = "IXYZ"
    if turn is None:
        for i, P in enumerate(letters):
            for Q in letters[i + 1:]:
                ax.plot(*xyz([corners[P], corners[Q]]), color=VIOLET, lw=1.3, alpha=0.8)
    else:
        others = [Q for Q in "XYZ" if Q != turn]
        ax.plot(*xyz([corners["I"], corners[turn]]), color=CYAN, lw=2.4)
        ax.plot(*xyz([corners[others[0]], corners[others[1]]]), color=PINK, lw=2.4)
        i = {"Z": 0, "X": 1, "Y": 2}[turn]
        lo, hi = [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]
        lo[i], hi[i] = -0.3, 1.45
        ax.plot(*xyz([lo, hi]), color=AMBER, lw=1.1, ls=(0, (4, 3)))
    k = np.array([corners[P] for P in letters])
    ax.scatter(*xyz(k), s=150 if big else 70, color=TEXT, edgecolors=BG, depthshade=False, zorder=5)
    for P in letters:
        x, y, z = xyz(np.array(corners[P]) * 1.24 - 0.12 if big else np.array(corners[P]) + (0.17, 0.0, 0.0))
        ax.text(x, y, z, P, fontsize=15 if big else 11, ha="center", va="center", color=TEXT, clip_on=False)
    if big:
        for name, pos in ((r"$k_X$", (0, 0.5, -0.14)), (r"$k_Y$", (0, -0.14, 0.5)), (r"$k_Z$", (0.5, -0.1, -0.1))):
            x, y, z = xyz(np.array(pos))
            ax.text(x, y, z, name, fontsize=10, color=FAINT if False else "#4f8794", ha="center", va="center")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
    ax.view_init(elev=18, azim=-58)


def draw_cube(corners):
    fig = plt.figure(figsize=(14, 6.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.7, 1, 1, 1], left=0.0, right=0.99, top=0.86, bottom=0.3, wspace=0.0)
    ax = fig.add_subplot(gs[0], projection="3d", facecolor=BG)
    cube(ax, corners, big=True)
    fig.text(0.03, 0.92, "One site: the four letters on the even corners of the cube", color=TITLE, fontsize=13)
    fig.text(0.03, 0.105, r"filled: I (0,0,0), X (1,0,1), Y (1,1,0), Z (0,1,1) in $(k_Z, k_X, k_Y)$; $k_P = 1$ where the letter anticommutes with P", fontsize=10)
    fig.text(0.03, 0.07, "hollow: the four odd corners, which no letter reaches; a move flipping one coordinate alone would land there", color="#4f8794", fontsize=10)
    fig.text(0.03, 0.035, r"a letter decays at $2(\gamma_Z k_Z + \gamma_X k_X + \gamma_Y k_Y)$ under dephasing along all three", fontsize=10)
    fig.text(0.47, 0.92, "The three pairings of the letters are the three half-turns", color=TITLE, fontsize=13)
    for i, P in enumerate("ZXY"):
        a = fig.add_subplot(gs[1 + i], projection="3d", facecolor=BG)
        cube(a, corners, turn=P)
        others = [Q for Q in "XYZ" if Q != P]
        for row, (txt, col) in enumerate([(rf"about the $k_{P}$ axis", AMBER),
                                          (rf"I + {P}:  $2(\gamma_{others[0]} + \gamma_{others[1]})$", CYAN),
                                          (rf"{others[0]} + {others[1]}:  $2(\gamma_{others[0]} + \gamma_{others[1]}) + 4\gamma_{P}$", PINK),
                                          (rf"a pairing exactly when $\gamma_{P} = 0$", TEXT)]):
            a.text2D(0.2, 0.0 - 0.075 * row, txt, color=col, fontsize=10.5, transform=a.transAxes)
    path = os.path.join(OUT, "one_square_letter_cube_neon.png")
    fig.savefig(path, dpi=150); plt.close(fig)
    print("wrote", os.path.normpath(path))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    check_square()
    corners = check_cube()
    draw_square()
    draw_cube(corners)
