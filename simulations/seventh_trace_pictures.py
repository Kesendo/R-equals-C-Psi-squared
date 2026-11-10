"""Draw three pictures for the seventh centred trace at second order in the rates
(docs/proofs/PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN.md, Theorem B and §3).

  1. the weights: how a hopping bond counts beside a triangle (l, k, j), by how it touches the triangle,
     with the product of the fifth and second traces in and taken out, and the borrowed triangle;
  2. the ring of five: the five rate points as a pentagon and as a pentagram, the part of degree two
     being 2520 · 4^N · Im(loop) · (pentagram − pentagon);
  3. the time circle of §3 with its two marks, for one two-step path in two placements of its hops.

What is calculated: the weights, from the constants of the gate simulations/seventh_moment_transport.py;
for the ring of five the two sums and, exactly, the part of degree two of Tr((L + c/2)^7) from the
1024 × 1024 generator (the gate's own routine), set against the closed form. The time circle is a
drawing of the proof's rule (a site is excited from its σ⁺ forward to its next σ⁻) and calculates nothing.
What is a drawing choice: the graphs, the rate points, the placements of the hops, the colours.

  python simulations/seventh_trace_pictures.py

Writes visualizations/seventh_trace_weights_neon.png, seventh_trace_ring_of_five_neon.png and
seventh_trace_time_circle_neon.png.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Arc, FancyArrowPatch, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fifth_moment_flux_area as g  # noqa: E402
import seventh_moment_transport as gate  # noqa: E402

OUT = os.path.join(HERE, "..", "visualizations")

BG, PANEL, TEXT, FAINT, TITLE = "#080b12", "#0c1018", "#7ef9ff", "#1d3a46", "#39ff14"
CYAN, PINK, AMBER, VIOLET, GREY = "#08f7fe", "#fe53bb", "#ffe700", "#b14aed", "#4f8794"
WEIGHT = LinearSegmentedColormap.from_list("w", [PINK, VIOLET, "#4a6cf7", CYAN])
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": TEXT, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT, "ytick.color": TEXT, "axes.edgecolor": TEXT,
                     "figure.facecolor": BG, "savefig.facecolor": BG, "axes.facecolor": PANEL})


# ---- 1. the weights ---------------------------------------------------------------------------------------

def weights():
    """The weight of |t|² for a hopping bond, by how it touches the triangle (l, k, j): Theorem B's W."""
    C = gate.CONSTANTS
    w = {"off": C["bond"],
         "at l or k": C["bond"] - C["touch"],
         "at j": C["bond"] - C["third"],
         "edge through j": C["bond"] - C["touch"] - C["third"],
         "edge l–k": C["bond"] - 2 * C["touch"] - C["own"]}
    assert w == {"off": 30, "at l or k": 20, "at j": 25, "edge through j": 15, "edge l–k": 8}
    return w


NODES = {"l": (0.0, 0.0), "k": (1.0, 0.0), "j": (0.5, 0.87), "a": (-0.62, -0.42), "b": (1.62, -0.42),
         "c": (0.5, 1.62), "d": (1.55, 1.25), "e": (2.15, 0.75)}
BONDS = [("l", "k", "edge l–k"), ("l", "j", "edge through j"), ("j", "k", "edge through j"),
         ("l", "a", "at l or k"), ("k", "b", "at l or k"), ("j", "c", "at j"), ("d", "e", "off")]


def node(ax, p, name, strong=True):
    ax.scatter(*p, s=300, color=PANEL, edgecolors=TEXT if strong else GREY, linewidths=1.4, zorder=4)
    if name:
        ax.text(*p, name, ha="center", va="center", fontsize=11, color=TEXT if strong else GREY, zorder=5)


def weight_panel(ax, w, shift, title):
    lo, hi = min(w.values()) - shift, max(w.values()) - shift
    for a, b, kind in BONDS:
        v = w[kind] - shift
        pa, pb = np.array(NODES[a]), np.array(NODES[b])
        ax.plot(*zip(pa, pb), color=WEIGHT((v - lo) / (hi - lo)), lw=3.2, zorder=2, solid_capstyle="round")
        mid = (pa + pb) / 2
        nrm = np.array([-(pb - pa)[1], (pb - pa)[0]]); nrm = nrm / np.linalg.norm(nrm)
        if np.dot(nrm, mid - np.array([0.5, 0.29])) < 0:
            nrm = -nrm
        ax.text(*(mid + 0.13 * nrm), f"{v:+d}" if shift and v else str(v), ha="center", va="center", fontsize=12, color=TEXT)
    for name, p in NODES.items():
        node(ax, p, name if name in "lkj" else "", strong=name in "lkj")
    ax.set_title(title, color=TITLE, fontsize=10.5, loc="left")
    ax.set_xlim(-0.95, 2.45); ax.set_ylim(-0.75, 1.9); ax.set_aspect("equal"); ax.set_axis_off()


def borrowed_panel(ax):
    P = {"i": (0.0, 0.0), "k": (1.0, 0.0), "j": (0.5, 0.87), "l": (-0.95, 0.35)}
    for a, b in (("i", "k"), ("k", "j"), ("j", "i")):
        ax.add_patch(FancyArrowPatch(P[a], P[b], arrowstyle="-|>", mutation_scale=13, color=CYAN, lw=2.0, shrinkA=11, shrinkB=11, zorder=3))
    ax.plot(*zip(P["l"], P["i"]), color=PINK, lw=3.2, zorder=2)
    ax.text(-0.5, 0.33, r"$-15\,|t_{li}|^2$", ha="center", fontsize=11.5, color=PINK)
    ax.add_patch(FancyArrowPatch(P["l"], P["k"], connectionstyle="arc3,rad=0.5", arrowstyle="-", color=AMBER, lw=1.1, ls=(0, (4, 3)), zorder=1))
    ax.text(0.0, -0.62, "the pair (l, k)", ha="center", fontsize=10, color=AMBER)
    ax.text(0.5, 0.3, r"$\mathrm{Im}(t_{ik}\,t_{kj}\,t_{ji})$", ha="center", va="center", fontsize=10.5, color=CYAN)
    for name, p in P.items():
        node(ax, p, name)
    ax.set_title("the borrowed triangle: l joined to the corner i", color=TITLE, fontsize=10.5, loc="left")
    ax.set_xlim(-1.4, 1.5); ax.set_ylim(-0.85, 1.35); ax.set_aspect("equal"); ax.set_axis_off()


def draw_weights():
    w = weights()
    fig = plt.figure(figsize=(14, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1], left=0.02, right=0.99, top=0.84, bottom=0.2, wspace=0.05)
    weight_panel(fig.add_subplot(gs[0]), w, 0, r"the weight of $|t|^2$ for every hopping bond of the graph")
    weight_panel(fig.add_subplot(gs[1]), w, w["off"], "with the fifth trace times the second taken out")
    borrowed_panel(fig.add_subplot(gs[2]))
    fig.text(0.03, 0.93, r"The seventh centred trace at second order in the rates: what a bond adds beside a triangle (l, k, j)", color=TITLE, fontsize=13)
    fig.text(0.03, 0.115, r"the pair (l, k) reads the loop $\mathrm{Im}(t_{lk}\,t_{kj}\,t_{jl})$ times the sum of these weights:  30 for every bond, 10 less for each touch of l or k, 5 less for a touch of j, the edge l–k two less", fontsize=10)
    fig.text(0.03, 0.075, "middle: 30 on every bond is the fifth trace times the second; what is left, the cumulant, is zero for a bond that does not touch the triangle", fontsize=10)
    fig.text(0.03, 0.035, "right: a site l joined to one corner i of a triangle pairs with each of the two corners it is not joined to, reading that triangle's loop in i's place", color=PINK, fontsize=10)
    path = os.path.join(OUT, "seventh_trace_weights_neon.png")
    fig.savefig(path, dpi=150); plt.close(fig)
    print("wrote", os.path.normpath(path))


# ---- 2. the ring of five -----------------------------------------------------------------------------------

RING = {(0, 1): 1 + 0j, (1, 2): 1 + 0j, (2, 3): 1 + 0j, (3, 4): 1 + 0j, (4, 0): 1j}     # flux π/2 through the ring
DOWN, UP = [2, 6, 8, 5, 1], [1, 1, 4, 7, 5]                                             # the five rate points


def ring_numbers():
    cross = lambda x, y: DOWN[x] * UP[y] - DOWN[y] * UP[x]
    pentagon = sum(cross(x, (x + 1) % 5) for x in range(5))
    pentagram = sum(cross(x, (x + 2) % 5) for x in range(5))
    loop = 1
    for x in range(5):
        loop *= g.t_of(RING, x, (x + 1) % 5)
    want = 2520 * 4 ** 5 * int(loop.imag) * (pentagram - pentagon)
    s = [d + u for d, u in zip(DOWN, UP)]; a = [d - u for d, u in zip(DOWN, UP)]
    parts = gate.rate_degrees(5, g.hamiltonian(5, RING, {}), s, a)
    assert parts[2] == 128 * want and want != 0, (parts, want)
    assert gate.closed_form(5, RING, {}, s, a) == 128 * want
    print(f"ring of five: pentagon {pentagon}, pentagram {pentagram}, Im(loop) {int(loop.imag)}; "
          f"degree-two part of Tr((L + c/2)^7) = {want}, exact against the generator and the gate's closed form")
    return pentagon, pentagram, int(loop.imag), want, parts


def cut(p1, p2, p3, p4):
    """Where the lines p1 p2 and p3 p4 meet."""
    (x1, y1), (x2, y2), (x3, y3), (x4, y4) = p1, p2, p3, p4
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    a, b = x1 * y2 - y1 * x2, x3 * y4 - y3 * x4
    return ((a * (x3 - x4) - (x1 - x2) * b) / den, (a * (y3 - y4) - (y1 - y2) * b) / den)


def rate_panel(ax, order_step, col, title):
    pts = list(zip(DOWN, UP))
    for m in (0.15, 0.3, 0.5, 0.8, 1.3, 2.5, 6):
        ax.plot([0, 12], [0, 12 * m], color=FAINT, lw=0.6, zorder=1)
    poly = [pts[(i * order_step) % 5] for i in range(5)]
    ax.add_patch(Polygon(poly, closed=True, facecolor=col, alpha=0.22, edgecolor=col, lw=1.8, zorder=3))
    if order_step == 2:                                      # the middle of the star is wound round twice
        inner = [cut(pts[i], pts[(i + 2) % 5], pts[(i + 1) % 5], pts[(i + 3) % 5]) for i in range(5)]
        ax.add_patch(Polygon(inner, closed=True, facecolor=col, alpha=0.22, edgecolor="none", zorder=3))
    for i, p in enumerate(pts):
        ax.scatter(*p, s=60, color=TEXT, edgecolors=BG, zorder=5)
        ax.annotate(str(i), p, textcoords="offset points", xytext=(7, 6), fontsize=10, color=TEXT)
    ax.set_xlim(0, 9.5); ax.set_ylim(0, 8); ax.set_aspect("equal")
    ax.set_xlabel(r"loss $\gamma_\downarrow$"); ax.set_ylabel(r"gain $\gamma_\uparrow$")
    ax.set_title(title, color=col, fontsize=10.5, loc="left")


def ring_panel(ax):
    ang = np.pi / 2 - 2 * np.pi * np.arange(5) / 5
    pts = np.c_[np.cos(ang), np.sin(ang)]
    for x in range(5):
        y = (x + 1) % 5
        ax.add_patch(FancyArrowPatch(pts[x], pts[y], arrowstyle="-|>", mutation_scale=13, color=CYAN, lw=1.9, shrinkA=12, shrinkB=12, zorder=3))
        mid = (pts[x] + pts[y]) / 2
        ax.text(*(1.3 * mid), "i" if x == 4 else "1", ha="center", va="center", fontsize=11, color=TEXT)
    for i, p in enumerate(pts):
        node(ax, p, str(i))
    ax.text(0, 0, r"$\Phi = \pi/2$", ha="center", va="center", fontsize=12, color=AMBER)
    ax.set_title("a ring of five hopping bonds with a flux", color=TITLE, fontsize=10.5, loc="left")
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.35, 1.4); ax.set_aspect("equal"); ax.set_axis_off()


def draw_ring():
    pentagon, pentagram, loop, want, _ = ring_numbers()
    fig = plt.figure(figsize=(14, 5.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1, 1], left=0.02, right=0.98, top=0.84, bottom=0.26, wspace=0.2)
    ring_panel(fig.add_subplot(gs[0]))
    rate_panel(fig.add_subplot(gs[1]), 1, CYAN, f"the pentagon, neighbours joined: {pentagon}")
    rate_panel(fig.add_subplot(gs[2]), 2, PINK, f"the pentagram, every second joined: {pentagram}")
    fig.text(0.03, 0.93, "The ring of five in the seventh centred trace: the pentagram of the rate points minus their pentagon", color=TITLE, fontsize=13)
    fig.text(0.03, 0.13, rf"part of degree two in the rates of $\mathrm{{Tr}}\,(L + c/2)^7$  =  $2520 \cdot 4^5 \cdot \mathrm{{Im(loop)}} \cdot$ (pentagram $-$ pentagon)  =  "
                         rf"$2520 \cdot 1024 \cdot {loop} \cdot ({pentagram} - {pentagon})$  =  {want}", fontsize=11.5, color=AMBER)
    fig.text(0.03, 0.085, "the five sites' rate points, numbered the way the loop runs; each number is twice a signed area, and the middle of the star counts twice (shaded twice)", fontsize=10)
    fig.text(0.03, 0.045, "two sites of the ring count +15 as neighbours and −15 two steps apart: that is the pentagon against the pentagram", fontsize=10)
    path = os.path.join(OUT, "seventh_trace_ring_of_five_neon.png")
    fig.savefig(path, dpi=150); plt.close(fig)
    print("wrote", os.path.normpath(path))


# ---- 3. the time circle ------------------------------------------------------------------------------------

def on_circle(frac, r=1.0):
    """The point at the fraction `frac` of the circumference, clockwise from the top."""
    a = np.pi / 2 - 2 * np.pi * frac
    return np.array([r * np.cos(a), r * np.sin(a)])


def arc(ax, start, end, r, col, lw=5):
    """The arc from the fraction `start` forward (clockwise) to `end`, wrapping past the top if end < start."""
    if end < start:
        end += 1
    f = np.linspace(start, end, 120)
    a = np.pi / 2 - 2 * np.pi * f
    ax.plot(r * np.cos(a), r * np.sin(a), color=col, lw=lw, solid_capstyle="butt", zorder=3)


def circle_panel(ax, y, p1, p2, title):
    t = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(t), np.sin(t), color=TEXT, lw=1.2, zorder=2)
    for frac, lab, col in ((0.0, r"$\sigma^+_l$ at 0", AMBER), (y, r"$\sigma^-_k$ at y", AMBER)):
        p = on_circle(frac)
        ax.scatter(*p, s=110, marker="D", color=col, edgecolors=BG, zorder=6)
        ax.text(*on_circle(frac, 0.74), lab, ha="center", va="center", fontsize=10.5, color=col)
    for frac, lab in ((p1, r"hop $l \to j$"), (p2, r"hop $j \to k$")):
        ax.scatter(*on_circle(frac), s=60, color=TEXT, edgecolors=BG, zorder=6)
        ax.text(*on_circle(frac, 0.72), lab, ha="center", va="center", fontsize=9.5, color=TEXT)
    arc(ax, 0.0, p1, 1.13, CYAN); arc(ax, p1, p2, 1.26, VIOLET); arc(ax, p2, y, 1.39, PINK)
    ax.text(*on_circle(y / 2, 0.22), "arc y", ha="center", va="center", fontsize=10, color=GREY)
    ax.text(*on_circle((1 + y) / 2, 0.22), "arc x", ha="center", va="center", fontsize=10, color=GREY)
    ax.add_patch(FancyArrowPatch(on_circle(0.02, 1.58), on_circle(0.1, 1.58), connectionstyle="arc3,rad=-0.15", arrowstyle="-|>",
                                 mutation_scale=11, color=GREY, lw=1.0))
    ax.set_title(title, color=TITLE, fontsize=10.5, loc="left")
    ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.6, 1.75); ax.set_aspect("equal"); ax.set_axis_off()


def draw_circle():
    fig = plt.figure(figsize=(12.5, 6.6))
    gs = fig.add_gridspec(1, 2, left=0.02, right=0.98, top=0.86, bottom=0.2, wspace=0.05)
    circle_panel(fig.add_subplot(gs[0]), 0.62, 0.2, 0.42, "both hops on the arc y")
    circle_panel(fig.add_subplot(gs[1]), 0.45, 0.8, 0.25, "the first hop on the arc x: the excitations wrap round")
    fig.text(0.03, 0.94, r"The time circle: one circumference $\tau = x + y$, two marks, and the two-step path $l \to j \to k$", color=TITLE, fontsize=13)
    fig.text(0.03, 0.13, "outer arcs: where a site is excited, from its raising forward to its lowering:", fontsize=10)
    fig.text(0.50, 0.13, "l", color=CYAN, fontsize=11); fig.text(0.52, 0.13, "j", color=VIOLET, fontsize=11); fig.text(0.54, 0.13, "k", color=PINK, fontsize=11)
    fig.text(0.03, 0.095, "every placement of the two hops round the circle is allowed, each site's raising and lowering alternating;", fontsize=10)
    fig.text(0.03, 0.06, "so a simple path of h hops weighs the full volume, a power of the circumference alone", fontsize=10)
    fig.text(0.03, 0.025, "a site visited twice carries four marks, and only the placements in which they alternate count: the weights of the bonds beside a path come from there", fontsize=10)
    path = os.path.join(OUT, "seventh_trace_time_circle_neon.png")
    fig.savefig(path, dpi=150); plt.close(fig)
    print("wrote", os.path.normpath(path))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    draw_weights()
    draw_ring()
    draw_circle()
