"""felt_time D, extend the sin^2 shape-miss to high N (Tom 2026-06-19): does the density-mode profile
keep converging to the continuum cos/sin^2 harmonic as N grows? Proof table has 0.17/0.12/0.08/0.06 at
N=4/5/6/7; push to N=8, N=9 (the dense ceiling on 128GB). PYTHON/numpy => no .NET managed-array/int32
limit (the half-filling block is C(N,N//2)^2: N=8->4900, N=9->15876; 15876^2=2.5e8 elements < int.MaxValue,
~4GB dense, eig ~16GB -- feasible on 128GB, and numpy allocates natively anyway).

Matches the proof's method exactly: the survivor = slowest mode of the half-filling (p,p) block (p=N//2);
its density profile n(j); grad^2(b)=(n(j)-n(j+1))^2; sin-pred(b)=|sin(pi*(b0+1)/N)|; shape-miss =
max | grad^2/max(grad^2) - sinpred^2/max(sinpred^2) | (the exact _felt_time_amplitude_law.py formula)."""
import sys
import time
import importlib.util
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


st = _load("stone", "simulations/stone_survivor_alpha_closure.py")


def density_profile(N, p, J, g, bnds):
    """Per-site occupation n(j) of the slowest (p,p)-block mode (dense eig inside st.mode_embed)."""
    M, reMode = st.mode_embed(N, p, p, J, g, bnds)
    diag = np.real(np.diag(M))
    n = np.zeros(N)
    for s in range(1 << N):
        if abs(diag[s]) > 0.0:
            for j in range(N):
                if (s >> j) & 1:
                    n[j] += diag[s]
    return n, float(reMode)


def nz(a):
    a = np.asarray(a, float)
    return a / max(np.max(np.abs(a)), 1e-15)


def main():
    Q, J = 1.5, 1.0
    g = 1.0 / Q
    print("=== felt_time D: sin^2 shape-miss vs N (half-filling survivor, Q=1.5 chain) ===")
    print(f"{'N':>3} {'sector':>7} {'C(N,p)^2':>9} {'<n_XY>':>7} {'shape-miss':>11} {'sec':>7}", flush=True)
    for N in (4, 5, 6, 7, 8, 9):
        p = N // 2
        bnds = bonds(N, "chain")
        from math import comb
        D = comb(N, p) ** 2
        t0 = time.time()
        n, reMode = density_profile(N, p, J, g, bnds)
        grad2 = np.array([(n[a] - n[c]) ** 2 for a, c in bnds])
        sinpred = np.array([abs(np.sin(np.pi * (b[0] + 1) / N)) for b in bnds])
        miss = float(np.max(np.abs(nz(grad2) - nz(sinpred ** 2))))
        nxy = -reMode / (2 * g)
        dt = time.time() - t0
        print(f"{N:>3} ({p},{p}) {D:>9} {nxy:>7.3f} {miss:>11.3f} {dt:>7.1f}", flush=True)
    print("\nREAD: proof table is 0.17/0.12/0.08/0.06 for N=4..7. If N=8, N=9 continue downward, the",
          flush=True)
    print("continuum cos/sin^2 harmonic is confirmed as the large-N limit of the density profile.", flush=True)


if __name__ == "__main__":
    main()
