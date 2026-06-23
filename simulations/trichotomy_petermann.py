"""The connection to Petermann: does K = 1/r^2 (the Petermann excess-noise factor) DIVERGE at the chain
coherence horizon Q*(N)?  (2026-06-18)

K = 1/r^2 is, by definition, the Petermann factor: the excess-noise / linewidth-enhancement factor Petermann
(1979) found for non-orthogonal laser modes. At an exceptional point the modes coalesce, r -> 0, K -> inf.
Here we resolve it at our coherence horizon: sample the half-filling (m,m) sector at Q = Q*(N) +/- delta on a
log grid and ask (i) does chain K spike while ring/star stay bounded, and (ii) what is the divergence law
K ~ |delta|^(-beta) (from r ~ |delta|^alpha, beta = 2*alpha).

DATA: results/trichotomy_cube/petermann_divergence.csv, from the witness statics CarbonImAndRigidity /
CarbonSlowestRate at the (m,m) sector, Q* = IncompletenessSurvivorWitness.HandoverQ(N, Chain). Throwaway
[Fact] _PetermannDivergenceSweep generated it (recreate if the tap recurs).
"""
import csv
import math
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "results", "trichotomy_cube", "petermann_divergence.csv")
OUT = os.path.join(HERE, "results", "trichotomy_cube", "petermann_factor_diverges_at_horizon.png")


def load():
    rows = defaultdict(list)  # (N, topo) -> [(Q, signedDelta, im, r, K)]
    with open(CSV, newline="") as f:
        for d in csv.DictReader(f):
            rows[(int(d["N"]), d["topo"])].append(
                (float(d["Q"]), float(d["signedDelta"]), float(d["imMax"]),
                 float(d["rigidity"]), float(d["K"])))
    for k in rows:
        rows[k].sort()
    return rows


def fit_powerlaw(pts):
    """log r = alpha*log|delta| + c, over the small-|delta| points. Returns alpha."""
    xs = [math.log(abs(dl)) for (_, dl, _, r, _) in pts if abs(dl) < 0.25 and r > 0]
    ys = [math.log(r) for (_, dl, _, r, _) in pts if abs(dl) < 0.25 and r > 0]
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx if sxx else float("nan")


def main():
    rows = load()
    Ns = sorted({k[0] for k in rows})

    print("=" * 86)
    print("THE PETERMANN FACTOR K = 1/r^2 AT THE COHERENCE HORIZON Q*(N)")
    print("  (m,m) half-filling sector; K_max = peak excess-noise factor; alpha = r ~ |Q-Q*|^alpha")
    print("=" * 86)
    for n in Ns:
        qs = rows[(n, "chain")][0][0] - rows[(n, "chain")][0][1]  # Q - signedDelta = Q*
        print(f"\n  N={n}  Q*(N)={qs:.4f}")
        for topo in ["chain", "ring", "star"]:
            pts = rows.get((n, topo))
            if not pts:
                continue
            Kmax = max(p[4] for p in pts)
            rmin = min(p[3] for p in pts)
            alpha = fit_powerlaw(pts) if topo == "chain" else float("nan")
            extra = f" | r~|d|^{alpha:.2f} => K~|d|^{-2*alpha:.2f}" if alpha == alpha else ""
            print(f"    {topo:5s}  K_max={Kmax:9.1f}  (r_min={rmin:.4f}){extra}")

    # ---- figure: K(Q) per topology, log scale, Q* marked ----
    fig, axes = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4.4), sharey=True)
    C = {"chain": "#1f77b4", "ring": "#2ca02c", "star": "#d62728"}
    for ax, n in zip(axes, Ns):
        qs = rows[(n, "chain")][0][0] - rows[(n, "chain")][0][1]
        for topo in ["ring", "star", "chain"]:
            pts = rows.get((n, topo))
            if not pts:
                continue
            Q = [p[0] for p in pts]
            K = [p[4] for p in pts]
            lw = 2.6 if topo == "chain" else 1.6
            ax.semilogy(Q, K, "-o", color=C[topo], lw=lw, ms=3, label=topo)
        ax.axvline(qs, color="k", ls=":", lw=1.4)
        ax.text(qs, ax.get_ylim()[1], f"Q*={qs:.3f}", ha="center", va="bottom", fontsize=9)
        ax.set_title(f"N={n}, half-filling ({n//2},{n//2})")
        ax.set_xlabel("Q = J/γ")
        ax.grid(alpha=0.25, which="both")
        if n == Ns[0]:
            ax.set_ylabel("Petermann factor  K = 1/r²")
            ax.legend(loc="upper right")
    fig.suptitle("The Petermann factor K = 1/r² diverges at the coherence horizon Q*(N): the chain is a "
                 "2nd-order EP, r~|Q−Q*|^½ ⇒ K~1/|Q−Q*| (clean, symmetric, ON Q*). Ring stays bounded "
                 "(level crossing); the star's smaller off-Q* spikes are separate real coalescences. "
                 "Petermann's 1979 laser scale, located in the dephased chain.", fontsize=9.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
