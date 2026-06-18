"""Derive-or-debunk: does <n_XY> ~ c*Q^2/N^2 hold, and is c the dispersion constant (no incompleteness 1/2)?

c_eff = <n_XY> * N^2 / Q^2. If the law holds: c_eff is Q-flat (at low Q) and N-flat; ring/chain -> 4;
star breaks (flat band, no dispersion). The dispersion+perturbation prediction is c_chain = pi^2/2 * (a
diffusion-constant factor), so check the MEASURED c against pi^2/2 and whether a clean 1/2 is needed.
DATA: results/survivor_scaling/survivor_scaling.csv (IncompletenessSurvivorWitness.Survivor, low Q, N=4..7).
"""
import csv, math, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "results", "survivor_scaling", "survivor_scaling.csv")

rows = defaultdict(list)  # (topo, N) -> [(Q, nXy, c_eff, pc, pr)]
with open(CSV, newline="") as f:
    for r in csv.DictReader(f):
        rows[(r["topo"], int(r["N"]))].append(
            (float(r["Q"]), float(r["nXy"]), float(r["c_eff"]), int(r["pc"]), int(r["pr"])))
for k in rows:
    rows[k].sort()

Ns = sorted({k[1] for k in rows})
print("pi^2/2 =", round(math.pi**2 / 2, 4), "   pi^2 =", round(math.pi**2, 4))
print("=" * 78)
print("c_eff = <n_XY>*N^2/Q^2  across Q  (should plateau at low Q if <n_XY> ~ c*Q^2/N^2)")
print("=" * 78)
for topo in ["chain", "ring", "star"]:
    print(f"\n  {topo}")
    print("    N |   " + "   ".join(f"Q={q}" for q in [0.05, 0.08, 0.12, 0.18, 0.25, 0.35]))
    for n in Ns:
        seq = rows.get((topo, n))
        if not seq:
            continue
        cells = "  ".join(f"{c:6.2f}" for (_, _, c, _, _) in seq)
        sectors = "/".join(sorted({f"{pc}{pr}" for (_, _, _, pc, pr) in seq}))
        print(f"   {n}  |  {cells}   sec:{sectors}")

# low-Q plateau estimate: mean c_eff over the two lowest Q
print("\n" + "=" * 78)
print("c (low-Q plateau, mean of Q=0.05,0.08) vs N  -- N-flat? ring/chain -> 4?")
print("=" * 78)
cval = {}
for topo in ["chain", "ring", "star"]:
    print(f"\n  {topo}")
    for n in Ns:
        seq = rows.get((topo, n))
        if not seq:
            continue
        lowq = [c for (q, _, c, _, _) in seq if q <= 0.085]
        c = sum(lowq) / len(lowq) if lowq else float("nan")
        cval[(topo, n)] = c
        print(f"    N={n}: c={c:7.3f}   (c/(pi^2/2)={c/(math.pi**2/2):5.2f},  c*N^2-form)")

print("\n  ring/chain ratio of c (predict 4):")
for n in Ns:
    rc = cval.get(("ring", n))
    cc = cval.get(("chain", n))
    if rc and cc and cc == cc and rc == rc:
        print(f"    N={n}: ring/chain = {rc/cc:5.2f}")

print("\n  star: is c N-flat (law holds) or drifting (flat band, law breaks)?")
for n in Ns:
    print(f"    N={n}: c_star={cval.get(('star', n), float('nan')):7.3f}")
