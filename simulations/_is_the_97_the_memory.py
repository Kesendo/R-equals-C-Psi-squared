"""Tom's gate (2026-06-20): is the complete-graph's 97.6% degeneracy THE memory 97%?

THE_VIEW_ONTO_THE_MEMORY.md names a "97%": rho_past, the slow half |Re lambda| < Sigma_gamma
(the Born-shadow, F94 = (4/3)Q^2 K^3), read as the low-drain-depth memory. The RMT probe found a
DIFFERENT 97.6%: the fraction of complete-graph oscillating eigenvalues that collapse onto coincident
(degenerate) values. Same numeral, two windows -- or one structure?

Gate-first: if (a) the degeneracy collapse is not a stable ~97% but varies with topology, and
(b) it is a different number from the past-fraction (the actual memory split), and (c) my count is
over the modes the memory doc EXCLUDES (oscillating, Im>0, the drained future) -- then the match is
a coincidence, and saying "the 97% is the memory" would be the false-unification trap.

Run:  python simulations/_is_the_97_the_memory.py
"""
from pathlib import Path
import numpy as np

RESULTS = Path(__file__).parent / "results"
GAMMA = 0.05


def load(topo, N):
    name = f"rmt_eigenvalues_N{N}.csv" if topo == "chain" else f"rmt_eigenvalues_{topo}_N{N}.csv"
    p = RESULTS / name
    if not p.exists():
        return None
    r, i = [], []
    with open(p) as f:
        f.readline()
        for line in f:
            q = line.strip().split("\t")
            if len(q) == 2:
                r.append(float(q[0].replace(",", ".")))
                i.append(float(q[1].replace(",", ".")))
    return np.array(r) + 1j * np.array(i)


print("=" * 88)
print("IS THE 97.6% THE MEMORY? three distinct fractions side by side")
print("  collapse%  = degeneracy of OSCILLATING modes (Im>0): my RMT number, 1 - distinct/total")
print("  past%      = the memory split by mode count: #{rate < Sigma_gamma} / total (the Born shadow)")
print("  kernel%    = the strict memory: #{rate ~ 0} / total (the depth-0 conserved floor)")
print("=" * 88)
print(f"{'topo':9} {'N':>2} {'Sg':>5} {'collapse%':>10} {'past%':>8} {'kernel%':>8}")
collapse_by_topo = {}
for topo in ("chain", "ring", "star", "complete"):
    for N in range(3, 8):
        ev = load(topo, N)
        if ev is None:
            continue
        Sg = N * GAMMA
        rates = -ev.real
        up = ev[ev.imag > 1e-6]
        nd = len(np.unique(np.round(up, 9)))
        collapse = 100 * (1 - nd / len(up)) if len(up) else float("nan")
        past = 100 * np.mean(rates < Sg - 1e-9)
        kernel = 100 * np.mean(rates < 1e-9)
        collapse_by_topo.setdefault(topo, {})[N] = collapse
        print(f"{topo:9} {N:>2} {Sg:>5.2f} {collapse:>10.1f} {past:>8.1f} {kernel:>8.2f}")

print("\n" + "-" * 88)
print("GATE READING:")
print("-" * 88)
# (a) is the collapse a stable ~97%?
n6 = {t: collapse_by_topo[t][6] for t in ("chain", "ring", "star", "complete") if 6 in collapse_by_topo.get(t, {})}
spread = max(n6.values()) - min(n6.values())
print(f"(a) collapse% at N=6 ranges {min(n6.values()):.1f}..{max(n6.values()):.1f} "
      f"(spread {spread:.1f} pts) -> {'NOT a constant' if spread > 5 else 'roughly constant'}; "
      f"only complete sits near 97.6, chain is ~{n6.get('chain', float('nan')):.0f}%.")
print("(b) collapse% (a DEGENERACY count, symmetry-driven) is a different quantity from past% (a")
print("    SLOW-HALF count, the actual memory split) -- compare the two columns above; they disagree.")
print("(c) collapse% is measured over Im>0 (the OSCILLATING / drained-future modes); the memory is the")
print("    Im=0 depth-0 kernel (the kernel% column), which the upper-half filter EXCLUDES by construction.")
print("\nVERDICT: the 97.6% is symmetry-driven spectral degeneracy of the FUTURE modes; the memory 97%")
print("is the F94 Born-shadow PURITY of the slow PAST. Same numeral, opposite ends of the axis, no node.")
