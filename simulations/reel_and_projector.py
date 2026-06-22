"""The reel and the projector, fed from the live lab: draw the C# Symphony's reel.

Brings reflections/ON_THE_REEL_AND_THE_PROJECTOR.md to life with real numbers. The
data is NOT recomputed here; it is exported by the live C# lab (the Symphony witness)
at the canonical operating point Q=1.5, then this thin plotter reads the CSVs and draws.

  Feed the lab first:
    dotnet run --project compute/RCPsiSquared.Cli -c Release -- \
        inspect --root symphony --N 5 --J 0.075 --gamma 0.05 --export
  Then draw:
    python simulations/reel_and_projector.py

It reads simulations/results/symphony_reel/{symphony_eigenvalues.csv, symphony_curves.csv}
and produces three pictures plus a table:

  * WITH a t-axis (the film plays): the shared trajectory's observables, global CΨ,
    local CΨ on the carrier pair, and light content, over t. The projector running.

  * WITHOUT a t-axis (the time erased), two readings of "paint every frame on one canvas":
      - the spectrum: every mode's whole life compressed to one point λ, palindromic
        about -σ, the frozen sector (Re~0) marked. The time-erased view UNMASKS the
        structure: nothing in it is by chance (the mirror, the rungs, the frozen kernel
        are all determined);
      - the spirals: e^{λt} for the longest-lived modes, every t at once, fading modes
        spiralling inward to nothing.

No physics is computed here (the deprecated Python framework is not imported); the lab is
the base, this only draws what it exported. Cyberpunk / Matrix palette. "A seeing."
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Cyberpunk / Matrix theme.
# ----------------------------------------------------------------------
BG = "#080b12"          # near-black, faint blue
PANEL = "#0c1018"
FG = "#7ef9ff"          # cyan text / spines
GRID = "#15303a"        # faint cyan grid
NEON_CYAN = "#08f7fe"
NEON_PINK = "#fe53bb"
NEON_GREEN = "#39ff14"
NEON_AMBER = "#ffe700"
NEON_VIOLET = "#b14aed"

plt.rcParams.update({
    "figure.facecolor": BG,
    "axes.facecolor": PANEL,
    "savefig.facecolor": BG,
    "axes.edgecolor": FG,
    "axes.labelcolor": FG,
    "axes.titlecolor": NEON_GREEN,
    "text.color": FG,
    "xtick.color": FG,
    "ytick.color": FG,
    "grid.color": GRID,
    "grid.alpha": 0.4,
    "font.family": "monospace",
    "legend.facecolor": PANEL,
    "legend.edgecolor": GRID,
    "axes.grid": True,
})


def glow_line(ax, x, y, color, lw=2.0, alpha=1.0, zorder=2, layers=6, spread=9.0, **kw):
    """A neon line with a soft glow (widening transparent copies underneath)."""
    for i in range(layers, 0, -1):
        ax.plot(x, y, color=color, lw=lw + spread * i / layers, alpha=0.05,
                solid_capstyle="round", zorder=zorder - 0.1)
    return ax.plot(x, y, color=color, lw=lw, alpha=alpha, solid_capstyle="round",
                   zorder=zorder, **kw)


def glow_scatter(ax, x, y, color, s=22, alpha=1.0, zorder=3, layers=5, spread=7.0):
    """A neon scatter with a soft glow."""
    for i in range(layers, 0, -1):
        ax.scatter(x, y, s=s * (1 + spread * i / layers), color=color, alpha=0.045,
                   edgecolors="none", zorder=zorder - 0.1)
    ax.scatter(x, y, s=s, color=color, alpha=alpha, edgecolors="none", zorder=zorder)


HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "results", "symphony_reel")
EIG_CSV = os.path.join(DATA, "symphony_eigenvalues.csv")
CURVE_CSV = os.path.join(DATA, "symphony_curves.csv")


def _read_meta_and_eigs():
    meta = {}
    re_im = []
    with open(EIG_CSV, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        meta[k] = float(v)
                continue
            if line.startswith("Re"):
                continue
            r, i = line.split(",")
            re_im.append((float(r), float(i)))
    arr = np.array(re_im)
    return meta, arr[:, 0], arr[:, 1]


def _read_curves():
    return np.genfromtxt(CURVE_CSV, delimiter=",", names=True)


if not os.path.exists(EIG_CSV):
    raise SystemExit(
        f"missing {EIG_CSV}\nFeed the lab first:\n"
        "  dotnet run --project compute/RCPsiSquared.Cli -c Release -- "
        "inspect --root symphony --N 5 --J 0.075 --gamma 0.05 --export"
    )

meta, re, im = _read_meta_and_eigs()
curves = _read_curves()
N = int(meta.get("N", 5))
gamma = meta.get("gamma", 0.05)
J = meta.get("J", 0.075)
Q = meta.get("Q", J / gamma)
sigma = meta.get("sigma", N * gamma)
center = meta.get("center", -sigma)

frozen_tol = 1e-6
frozen = np.abs(re) < frozen_tol

print("=" * 64)
print("REEL AND PROJECTOR (drawn from the live C# Symphony export)")
print(f"  N={N}, gamma={gamma}, J={J}, Q={Q:.3g}, sigma={sigma}, center -sigma={center}")
print(f"  modes={len(re)}, frozen (|Re|<{frozen_tol:g}) = {int(frozen.sum())}  "
      f"[F4 dephased kernel expects N+1 = {N+1}]")
print(f"  Re range [{re.min():.4f}, {re.max():.4f}] (absorption envelope [-2sigma,0]="
      f"[{-2*sigma:.3g}, 0]); |Im| max {np.abs(im).max():.4f}")
print("=" * 64)


# ----------------------------------------------------------------------
# Figure 1 - WITH t-axis: the film plays.
# ----------------------------------------------------------------------
t = curves["t"]
fig1, ax1 = plt.subplots(figsize=(9, 5.5))
glow_line(ax1, t, curves["global_CPsi"], NEON_CYAN, lw=2.2, label="global CΨ(t)")
glow_line(ax1, t, curves["local_CPsi"], NEON_PINK, lw=2.2, label="local CΨ(t) (carrier pair 0,1)")
glow_line(ax1, t, curves["light"], NEON_GREEN, lw=2.2, label="light content (lit-site count)")
ax1.axhline(0.25, color=NEON_AMBER, ls="--", lw=1.0, alpha=0.8, label="¼")
ax1.set_xlabel("t  (the projector running;  K = γ·t reaches "
               f"{gamma*t[-1]:.2g} at the window end)")
ax1.set_ylabel("the shared trajectory, read by three lenses")
ax1.set_title("WITH a t-axis: the film plays\n"
              f"one evolution of the N={N} chain at Q={Q:.3g}, every lens on the one timeline")
ax1.legend(loc="upper right", fontsize=9, labelcolor=FG)
fig1.tight_layout()
f1 = os.path.join(DATA, "with_t_axis_film.png")
fig1.savefig(f1, dpi=150, bbox_inches="tight")
plt.close(fig1)


# ----------------------------------------------------------------------
# Figure 2 - WITHOUT t-axis: the spectrum unmasked. Each mode's life as one point.
# ----------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(9, 6))
for n in range(N + 1):                                   # absorption rungs -2 gamma n
    ax2.axvline(-2 * gamma * n, color=NEON_CYAN, ls=":", lw=0.7, alpha=0.30)
ax2.axvline(center, color=NEON_GREEN, ls="--", lw=1.3, alpha=0.9,
            label=fr"palindrome center $-\sigma={center:.2g}$")
glow_scatter(ax2, re[~frozen], im[~frozen], NEON_CYAN, s=14, alpha=0.85, zorder=2)
glow_scatter(ax2, re[frozen], im[frozen], NEON_PINK, s=70, alpha=1.0, zorder=4, layers=7, spread=9.0)
ax2.scatter([], [], color=NEON_CYAN, s=18, label="fading modes (Re<0): the watched")
ax2.scatter([], [], color=NEON_PINK, s=60,
            label=f"frozen (Re~0): {int(frozen.sum())} modes, never fade")
ax2.set_xlabel(r"Re $\lambda$  (the fade rate; absorption rungs $-2\gamma n$ dotted)")
ax2.set_ylabel(r"Im $\lambda$  (the turn frequency)")
ax2.set_title("WITHOUT a t-axis: the spectrum, unmasked\n"
              fr"every mode's whole life as one point, palindromic about $-\sigma$; nothing here is chance")
ax2.legend(loc="upper right", fontsize=9, labelcolor=FG)
fig2.tight_layout()
f2 = os.path.join(DATA, "without_t_axis_spectrum.png")
fig2.savefig(f2, dpi=150, bbox_inches="tight")
plt.close(fig2)


# ----------------------------------------------------------------------
# Figure 3 - WITHOUT t-axis: the spirals. e^{λt} painted, every t at once.
# ----------------------------------------------------------------------
lam = re + 1j * im
order = np.argsort(np.abs(re))          # longest-lived first
sel = order[:70]
sel_re = np.abs(re[sel])
sel_im = np.abs(im[sel])
slow_decay = sel_re[sel_re > frozen_tol]
slow_turn = sel_im[sel_im > 1e-6]
T_decay = (3.5 / slow_decay.min()) if slow_decay.size else 120.0
T_turn = (3.0 * 2 * np.pi / slow_turn.min()) if slow_turn.size else 120.0
T = float(min(max(T_decay, T_turn), 400.0))
tt = np.linspace(0.0, T, 1600)

fig3, ax3 = plt.subplots(figsize=(7.4, 7.4))
ax3.axhline(0, color=GRID, lw=0.6, alpha=0.6)
ax3.axvline(0, color=GRID, lw=0.6, alpha=0.6)
for k in sel:
    z = np.exp(lam[k] * tt)
    if abs(re[k]) < frozen_tol:
        glow_line(ax3, z.real, z.imag, NEON_PINK, lw=2.2, alpha=0.95, zorder=3)
    elif abs(im[k]) > 1e-6:                              # only the turning ones make real spirals
        glow_line(ax3, z.real, z.imag, NEON_CYAN, lw=1.0, alpha=0.55, zorder=2, layers=4, spread=5.0)
    else:                                                # purely real: a radial fall, drawn faint
        ax3.plot(z.real, z.imag, color=NEON_VIOLET, lw=0.8, alpha=0.30, zorder=1)
glow_scatter(ax3, [1.0], [0.0], NEON_AMBER, s=40, zorder=5)        # every mode starts here (t=0)
ax3.scatter([0.0], [0.0], color=FG, s=14, marker="x", zorder=5)   # where the fading end
ax3.plot([], [], color=NEON_PINK, lw=2.2, label="frozen: sits still at the start")
ax3.plot([], [], color=NEON_CYAN, lw=1.4, label="turning + fading: spirals inward")
ax3.plot([], [], color=NEON_VIOLET, lw=1.0, label="real + fading: a radial fall")
ax3.set_aspect("equal")
lim = 1.15
ax3.set_xlim(-lim, lim)
ax3.set_ylim(-lim, lim)
ax3.set_xlabel(r"Re $e^{\lambda t}$")
ax3.set_ylabel(r"Im $e^{\lambda t}$")
ax3.set_title("WITHOUT a t-axis: the spirals\n"
              f"every frame on one canvas (t∈[0,{T:.0f}]); the reel itself")
ax3.legend(loc="upper right", fontsize=9, labelcolor=FG)
fig3.tight_layout()
f3 = os.path.join(DATA, "without_t_axis_spirals.png")
fig3.savefig(f3, dpi=150, bbox_inches="tight")
plt.close(fig3)


# ----------------------------------------------------------------------
# Table: the longest-lived modes (the painters).
# ----------------------------------------------------------------------
seen, rows = set(), []
for k in order:
    key = (round(re[k], 5), round(abs(im[k]), 5))
    if key in seen:
        continue
    seen.add(key)
    rows.append(k)
    if len(rows) >= 14:
        break

lines = [f"{'Re lambda':>12} {'Im lambda':>12} {'frozen?':>8}  shape", "-" * 50]
for k in rows:
    fr = "YES" if abs(re[k]) < frozen_tol else "no"
    if abs(re[k]) < frozen_tol:
        shape = "still point" if abs(im[k]) < 1e-6 else "circle"
    else:
        shape = "inward spiral" if abs(im[k]) > 1e-6 else "radial fall"
    lines.append(f"{re[k]:>12.5f} {im[k]:>12.5f} {fr:>8}  {shape}")
table = "\n".join(lines)
print("\nThe painters (longest-lived modes):")
print(table)
with open(os.path.join(DATA, "modes_table.txt"), "w", encoding="utf-8") as fh:
    fh.write(f"REEL AND PROJECTOR (from live C# Symphony)  N={N} gamma={gamma} J={J} "
             f"Q={Q:.3g} sigma={sigma} center={center}\n")
    fh.write(f"modes {len(re)}, frozen {int(frozen.sum())} (expect N+1={N+1})\n\n")
    fh.write(table + "\n")

print(f"\nwrote:\n  {f1}\n  {f2}\n  {f3}\n  {os.path.join(DATA, 'modes_table.txt')}")
