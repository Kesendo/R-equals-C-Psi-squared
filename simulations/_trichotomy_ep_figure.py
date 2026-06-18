"""Figure: the Petermann rigidity makes the chain coherence-horizon Q*(N) a MEASURED exceptional point.

DATA: simulations/results/trichotomy_cube/rigidity_sweep.csv, produced by driving the SHIPPED witness'
public statics over a Q-grid (no re-derivation):
    TrichotomyWitness.CarbonImAndRigidity(topo, N, Q, pc, pr) -> (|Im|, rigidity)   [carbon Qh=0.5, γ=1/Q]
    TrichotomyWitness.CarbonSlowestRate  (topo, N, Q, pc, pr) -> rate
for topo in {chain,ring,star}, N in 4..6, sector in {(0,1), (m,m)=half-filling}, Q in linspace(0.6,6,60).
To regenerate the CSV: a ~25-line xUnit [Fact] looping those calls with a row-by-row StreamWriter
(the throwaway _TrichotomyRigiditySweep harness; recreate if the tap recurs).

Per topology, in the half-filling (m,m) coherence sector, plot rigidity r(Q) (solid, left axis) and the
chain's |Im| (dashed, right axis). The chain r collapses to ~0 coincident with the |Im| birth at Q*(N)
(an exceptional point); the ring keeps r bounded (a level crossing); the star dips but never oscillates.
"""
import csv
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = os.path.join(HERE, "results", "trichotomy_cube", "rigidity_sweep.csv")
THRESH = os.path.join(HERE, "results", "trichotomy_cube", "thresholds.csv")
OUT = os.path.join(HERE, "results", "trichotomy_cube", "coherence_horizon_is_an_EP.png")

QSTAR = {}
with open(THRESH, newline="") as f:
    for r in csv.DictReader(f):
        QSTAR[int(r["N"])] = float(r["chain_Qstar"])

series = defaultdict(list)
with open(SWEEP, newline="") as f:
    for r in csv.DictReader(f):
        series[(int(r["N"]), r["topo"], r["sector"])].append(
            (float(r["Q"]), float(r["imMax"]), float(r["rigidity"])))
for k in series:
    series[k].sort()

Ns = [4, 5, 6]
C = {"chain": "#1f77b4", "ring": "#2ca02c", "star": "#d62728"}
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)

for ax, n in zip(axes, Ns):
    half = n // 2
    sector = f"{half}_{half}"
    ax2 = ax.twinx()
    # rigidity: chain bold, ring/star thinner
    for topo, lw, a in [("ring", 1.6, 0.85), ("star", 1.6, 0.85), ("chain", 3.0, 1.0)]:
        seq = series.get((n, topo, sector))
        if not seq:
            continue
        Q = [s[0] for s in seq]
        rg = [s[2] for s in seq]
        ax.plot(Q, rg, "-", color=C[topo], lw=lw, alpha=a, label=topo)
        if topo == "chain":
            qrm, rm = min(((s[0], s[2]) for s in seq), key=lambda t: t[1])
            ax.plot([qrm], [rm], "o", color=C["chain"], ms=8, zorder=5)
            ax.annotate(f"r_min={rm:.3f}", (qrm, rm), textcoords="offset points",
                        xytext=(6, -2), fontsize=8, color=C["chain"])
    # chain |Im| on the right axis, to show the EP coincidence (r->0 AT |Im| birth)
    seq = series.get((n, "chain", sector))
    if seq:
        ax2.plot([s[0] for s in seq], [s[1] for s in seq], "--", color=C["chain"], lw=1.3, alpha=0.6)
    qs = QSTAR.get(n)
    if qs == qs:
        ax.axvline(qs, color="k", ls=":", lw=1.4)
        ax.text(qs, 1.03, f"Q*={qs:.3f}", ha="center", va="bottom", fontsize=9)
    ax.set_title(f"N={n}, half-filling sector ({half},{half})")
    ax.set_xlabel("Q = J/γ")
    ax.set_ylim(-0.03, 1.06)
    ax.grid(alpha=0.25)
    ax2.set_ylim(bottom=0)
    if n == Ns[0]:
        ax.set_ylabel("phase rigidity r  (solid; r→0 ⇒ EP)")
        ax.legend(loc="center right", framealpha=0.9)
    if n == Ns[-1]:
        ax2.set_ylabel("chain |Im λ_slow|  (dashed; oscillation onset)")

fig.suptitle("The coherence horizon Q*(N) is a measured exceptional point — chain r→0 at the |Im| birth; "
             "ring stays bounded (level crossing); star never oscillates (frozen commutant)",
             y=1.0, fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print(f"wrote {OUT}")
