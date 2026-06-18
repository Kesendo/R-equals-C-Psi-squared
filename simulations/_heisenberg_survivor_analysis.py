"""Dir fehlt Z: does ZZ lift the XY filling-degeneracy and pin the survivor to half-filling (incompleteness)?

Per (model, N, topo, Q): print the darkness <n_XY> for each interior (p,p) sector, mark the survivor
(min darkness). XY should be ~p-flat (degenerate); Heisenberg should dip to a minimum at p=N/2 (half-filling).
"""
import csv, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "results", "survivor_scaling", "heisenberg_survivor.csv")

d = defaultdict(dict)  # (model,N,topo,Q) -> {p: nXy}
ce = defaultdict(dict)
with open(CSV, newline="") as f:
    for r in csv.DictReader(f):
        key = (r["model"], int(r["N"]), r["topo"], float(r["Q"]))
        d[key][int(r["p"])] = float(r["nXy"])
        ce[key][int(r["p"])] = float(r["cEff"])

print("=" * 84)
print("PER-SECTOR survivor darkness <n_XY>  (survivor = MIN, marked *).  half-filling p=N/2")
print("=" * 84)
for topo in ["chain", "ring"]:
    for n in [4, 5, 6]:
        half = n / 2
        print(f"\n  {topo} N={n}   (half-filling p={half:g})")
        for model in ["XY", "Heis"]:
            for Q in [0.05, 0.1, 0.2]:
                row = d.get((model, n, topo, Q))
                if not row:
                    continue
                pmin = min(row, key=row.get)
                cells = []
                for p in sorted(row):
                    mark = "*" if p == pmin else " "
                    cells.append(f"p{p}={row[p]:.5f}{mark}")
                spread = max(row.values()) - min(row.values())
                tag = "DEGENERATE" if spread < 1e-9 else f"survivor p={pmin}"
                print(f"    {model:4s} Q={Q:<4}  " + "  ".join(cells) + f"   | {tag}")

# scaling of the survivor (min-sector) darkness, Heisenberg vs XY
print("\n" + "=" * 84)
print("SURVIVOR c_eff = <n_XY>*N^2/Q^2 at the min sector (XY pi^2/4~2.47; Heisenberg retunes)")
print("=" * 84)
for topo in ["chain", "ring"]:
    print(f"\n  {topo}")
    for model in ["XY", "Heis"]:
        for n in [4, 5, 6]:
            vals = []
            for Q in [0.05, 0.1, 0.2]:
                row = ce.get((model, n, topo, Q))
                drow = d.get((model, n, topo, Q))
                if not row:
                    continue
                pmin = min(drow, key=drow.get)
                vals.append((Q, pmin, row[pmin]))
            s = "  ".join(f"Q={q}:c={c:.3f}(p{p})" for q, p, c in vals)
            print(f"    {model:4s} N={n}:  {s}")
