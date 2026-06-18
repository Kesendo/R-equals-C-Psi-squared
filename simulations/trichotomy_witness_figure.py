"""The chain/ring/star trichotomy in one figure: survivor darkness <n_XY>(Q) at N=6.

Data is verbatim from the live witness:
    dotnet run --project compute/RCPsiSquared.Cli -c Release -- inspect --root trichotomy --N 6
(the carbon un-freeze RouteSweep). Re-run that command to regenerate the numbers.

The trichotomy in one reading: chain and ring REACH the -2g Absorption floor (<n_XY> -> 1)
and un-freeze there (the survivor switches to the oscillating (0,1) band edge); the star
SATURATES below the floor, on the structural ceiling g2 = 4/(N-1) = 0.8 at N=6, and stays
frozen at every Q. Three topologies, three freeze-routes, one picture.

Output: docs/figures/trichotomy_nxy_vs_q.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Q = np.array([1, 1.5, 2, 3, 6, 12, 25, 50], dtype=float)

# <n_XY>(Q) and the frozen(0)/oscillating-band-edge(1) flag, read off the witness render (N=6):
TOPO = {
    "chain": dict(nxy=[0.069, 0.161, 0.307, 1.0, 1.0, 1.0, 1.0, 1.0],
                  osc=[0, 0, 0, 1, 1, 1, 1, 1], color="#1f77b4"),
    "ring":  dict(nxy=[0.266, 0.632, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                  osc=[0, 0, 0, 1, 1, 1, 1, 1], color="#2ca02c"),
    "star":  dict(nxy=[0.228, 0.425, 0.567, 0.695, 0.774, 0.794, 0.799, 0.8],
                  osc=[0, 0, 0, 0, 0, 0, 0, 0], color="#d62728"),
}
CEILING = 4.0 / 5.0  # g2 = 4/(N-1) at N=6

fig, ax = plt.subplots(figsize=(8.4, 5.6))

ax.axhline(1.0, ls="--", lw=1.0, color="0.55")
ax.text(0.97, 1.008, "the -2γ floor  ⟨n_XY⟩ = 1  (the band edge)",
        color="0.4", fontsize=9, va="bottom")
ax.axhline(CEILING, ls=":", lw=1.3, color="#d62728", alpha=0.75)
ax.text(52, CEILING - 0.006, "star ceiling  g₂ = 4/(N−1) = 0.8",
        color="#d62728", fontsize=9, ha="right", va="top")

for name, d in TOPO.items():
    nxy = np.array(d["nxy"]); osc = np.array(d["osc"]); c = d["color"]
    ax.plot(Q, nxy, "-", color=c, lw=1.8, alpha=0.45, zorder=1, label=name)
    fz, ob = osc == 0, osc == 1
    ax.scatter(Q[fz], nxy[fz], s=58, color=c, zorder=3)                       # frozen interior
    ax.scatter(Q[ob], nxy[ob], s=58, facecolors="none", edgecolors=c,        # oscillating band edge
               marker="s", linewidths=1.6, zorder=3)

ax.annotate("chain un-freezes\nat Q*(6) ≈ 2.88", xy=(3, 1.0), xytext=(3.6, 0.52),
            fontsize=8.5, color="#1f77b4", ha="left",
            arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1.1))
ax.annotate("ring (2,2) reaches the\nfloor at Q_h = 2", xy=(2, 1.0), xytext=(1.0, 0.80),
            fontsize=8.5, color="#2ca02c", ha="left",
            arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1.1))
ax.annotate("star saturates on the\nceiling — never un-freezes", xy=(50, 0.8), xytext=(5.5, 0.27),
            fontsize=8.5, color="#d62728", ha="left",
            arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.1))

# marker key
ax.scatter([], [], s=58, color="0.3", label="●  frozen interior (Δn=0)")
ax.scatter([], [], s=58, facecolors="none", edgecolors="0.3", marker="s",
           linewidths=1.6, label="□  oscillating (0,1) band edge (Δn=1)")

ax.set_xscale("log")
ax.set_xticks([1, 2, 3, 6, 12, 25, 50])
ax.set_xticklabels(["1", "2", "3", "6", "12", "25", "50"])
ax.set_xlim(0.93, 58)
ax.set_ylim(0, 1.09)
ax.set_xlabel("Q = J/γ   (coupling-to-dephasing ratio, log scale)")
ax.set_ylabel("survivor darkness   ⟨n_XY⟩ = Re(λ) / (−2γ)")
ax.set_title("The chain / ring / star trichotomy in one view  (N = 6)\n"
             "does the longest-lived survivor reach the floor (un-freeze) or the ceiling (stay frozen)?",
             fontsize=11)
ax.legend(loc="center right", fontsize=8.5, framealpha=0.92)
ax.grid(True, which="both", alpha=0.2)
fig.tight_layout()

os.makedirs("docs/figures", exist_ok=True)
out = os.path.join("docs", "figures", "trichotomy_nxy_vs_q.png")
fig.savefig(out, dpi=150)
print("wrote", out)
