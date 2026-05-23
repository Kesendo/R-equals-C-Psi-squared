"""Two perspectives on the same Q-spiral: from above (F50 world) and from below (Cusp world).

The same trajectory data is rendered with two camera positions:

  Left:  elev=+22°  — observer sits above the F50 ceiling, looks down at the cusp.
         This is the default 3D plot view. Positions us in the dissipative world.

  Right: elev=−22°  — observer sits below the cusp, looks up at the F50 ceiling.
         This is what we *would* see if we were classical observers in the
         conserved (kernel) substrate looking out at the dissipative dynamics.

Tom 2026-05-19: "wenn es so wäre wie sie sagen, und wir klassisch sind, müsste ich
1/4 von unten sehen und nicht aus der 0.75 Welt von oberhalb"

Output: simulations/results/q_scaling_visual/q_sweep_two_perspectives.png
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def main() -> None:
    N = 5
    gamma = 0.05
    out_dir = Path("simulations/results/q_scaling_visual")
    cache_file = out_dir / f"q_sweep_data_N{N}.npz"
    if not cache_file.exists():
        raise FileNotFoundError(f"Cache missing: {cache_file}")

    data = np.load(cache_file)
    J_values = data["J_values"]
    slow_re = np.maximum(data["slow_re"], 1e-16)
    slow_nxy = np.maximum(data["slow_nxy"], 1e-16)
    Q_values = J_values / gamma

    log_J = np.log10(J_values)
    theta = (log_J - log_J.min()) / (log_J.max() - log_J.min()) * 2 * math.pi
    r = slow_nxy
    z = np.log10(slow_re)
    x_polar = r * np.cos(theta)
    y_polar = r * np.sin(theta)

    fig = plt.figure(figsize=(20, 9))
    perspectives = [
        (22, -55, "Sicht VON OBEN  (elev=+22°)\n"
                  "wir sitzen in der F50-Welt (dissipative Schicht)\n"
                  "Cusp liegt unten als Tiefe vor uns"),
        (-22, -55, "Sicht VON UNTEN  (elev=−22°)\n"
                   "wir sitzen im Cusp (Kernel-Schicht, klassisch-erhalten)\n"
                   "F50-Decke hängt über uns, Spirale steigt hinauf"),
    ]

    for idx, (elev, azim, title_suffix) in enumerate(perspectives):
        ax = fig.add_subplot(1, 2, idx + 1, projection="3d")

        theta_ref = np.linspace(0, 2 * math.pi, 120)
        for z_ref, lbl in [(-1, "F50 floor"), (-6, "10⁻⁶"),
                             (-12, "10⁻¹²"), (-15, "10⁻¹⁵")]:
            ax.plot(np.cos(theta_ref), np.sin(theta_ref),
                    np.full_like(theta_ref, z_ref),
                    color="black", alpha=0.15, lw=0.5)

        ax.scatter([1], [0], [-1], color="red", s=300, marker="X",
                    edgecolors="black", linewidths=1.5,
                    label="F50 vertex", zorder=10)

        sc = ax.scatter(x_polar, y_polar, z,
                         c=log_J, cmap="viridis", s=70,
                         edgecolors="black", linewidths=0.4, zorder=5)
        ax.plot(x_polar, y_polar, z, color="gray", alpha=0.5, lw=1.2)

        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_zlim(-16.5, 0.2)
        ax.set_xlabel("r · cos(θ)")
        ax.set_ylabel("r · sin(θ)")
        ax.set_zlabel("log₁₀|Re(λ)|")
        ax.set_title(title_suffix, fontsize=10)
        ax.view_init(elev=elev, azim=azim)
        if idx == 0:
            ax.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        f"Zwei Perspektiven auf dieselbe Q-Spirale  "
        f"(N={N} chain, γ₀={gamma} konstant)\n"
        f"Links: Beobachter in der dissipativen Welt blickt auf den Cusp herab.\n"
        f"Rechts: Beobachter im klassisch-erhaltenen Cusp blickt zur F50-Decke hinauf.",
        fontsize=12)
    plt.tight_layout()

    out_file = out_dir / "q_sweep_two_perspectives.png"
    plt.savefig(out_file, dpi=120, bbox_inches="tight")
    print(f"Saved: {out_file}")
    plt.close(fig)


if __name__ == "__main__":
    main()
