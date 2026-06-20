"""handshake M2c read-cost DIAGNOSTIC (2026-06-20): the gate in _handshake_read_cost.py FIRED
(slope -0.67, non-monotone jump at Q=8, K_peak*Q drifts 1.65->3.57). Verdict flagged two
suspects: the ~2/Q LAW is wrong, OR the operational definition of "cost of one read" is wrong.

This script diagnoses the operational definition WITHOUT loosening the gate (gate-flow: diagnose,
do not relax the pass criterion). The original took K_peak = argmax_K FI via find_peaks with a
height>=5%*max, prominence>=2%*max filter. At high Q the FI(K) curve develops large LATE coherent
revivals; the genuine FIRST dark->bright rotation is a SMALL early bump that falls below the 5%
threshold, so find_peaks's "first peak" silently becomes a later big revival (note: in the original
run K_peak == K_global in every row -- the tell). "Cost of ONE read" should be the dose to the
FIRST rotation regardless of how small its FI signal is.

DIAGNOSTIC: re-extract the cost three ways, no pass-criterion change, and SHOW which definition (if
any) the hypothesis K_cost*Q ~ const ~2 actually picks out:
  A  argmax FI            (the original; expected to drift -- the artifact)
  B  first LOCAL max of FI, no height threshold   (the honest "first rotation")
  C  first dose where dFI/dK first turns over (inflection->plateau onset), no threshold
Also print the theory anchor: the closed-system first-rotation dwell t ~ 1/J -> K_dwell = g/J = 1/Q,
so K_cost*Q ~ O(1) if the mechanism is right.

A firing of THIS diagnostic (all three drift) => the ~2/Q law is genuinely refuted on the FI
apparatus, report it as a negative result. A clean B or C (K_cost*Q ~ const) => the original
operational definition was the bug; fix _handshake_read_cost.py to use it and re-run the real gate.
"""
import importlib.util
import sys

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


blk = _load("blk", "simulations/_handshake_rk_block.py")
psi, hopping, L11, populations = blk.psi, blk.hopping, blk.L11, blk.populations


def fi_over_k(N, J, g, dJ, bnds, defect, kgrid):
    """Z-population Fisher info of a delta-J bond defect, per dose K, in the (1,1) block."""
    from scipy.linalg import eig
    psi1 = np.array([psi(1, a, N) for a in range(N)])
    v0 = np.outer(psi1, psi1).astype(complex).flatten()
    tg = kgrid / g
    HA = hopping(N, J, bnds)
    wA, RA = eig(L11(HA, g))
    nA = populations(wA, RA, np.linalg.inv(RA) @ v0, tg, N)
    HB = HA.copy()
    HB[defect[0], defect[1]] += dJ
    HB[defect[1], defect[0]] += dJ
    wB, RB = eig(L11(HB, g))
    nB = populations(wB, RB, np.linalg.inv(RB) @ v0, tg, N)
    dp = (nB - nA) / dJ
    fi = np.sum(dp ** 2 / np.maximum(nA, 1e-12), axis=0)
    return fi


def first_local_max(k, y):
    """First interior local maximum, NO height/prominence threshold. None if monotone."""
    for i in range(1, len(y) - 1):
        if y[i] >= y[i - 1] and y[i] > y[i + 1]:
            return k[i]
    return k[int(np.argmax(y))]


def first_turnover(k, y):
    """First dose where the rising FI slope first decreases (onset of the first plateau/peak),
    i.e. first i with dy[i] < dy[i-1] while dy>0. Catches the first rotation even pre-peak."""
    dy = np.diff(y)
    for i in range(1, len(dy)):
        if dy[i - 1] > 0 and dy[i] < dy[i - 1]:
            return k[i]
    return k[int(np.argmax(y))]


def main():
    N, J, dJ = 5, 1.0, 0.02
    bnds = bonds(N, "chain")
    defect = bnds[len(bnds) // 2]
    Qs = [1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 13.0, 20.0, 35.0]
    print(f"=== M2c read-cost DIAGNOSTIC: three cost definitions (N={N}, defect {defect}, Z) ===\n",
          flush=True)
    print(f"{'Q':>6} {'A:argmax':>9} {'A*Q':>7} {'B:1stMax':>9} {'B*Q':>7} "
          f"{'C:turnov':>9} {'C*Q':>7} {'1/Q':>7}", flush=True)
    A, B, C = [], [], []
    for Q in Qs:
        g = 1.0 / Q
        kgrid = np.linspace(0.001, 6.0, 6000)
        fi = fi_over_k(N, J, g, dJ, bnds, defect, kgrid)
        kA = kgrid[int(np.argmax(fi))]
        kB = first_local_max(kgrid, fi)
        kC = first_turnover(kgrid, fi)
        A.append(kA); B.append(kB); C.append(kC)
        print(f"{Q:>6.1f} {kA:>9.4f} {kA*Q:>7.2f} {kB:>9.4f} {kB*Q:>7.2f} "
              f"{kC:>9.4f} {kC*Q:>7.2f} {1.0/Q:>7.4f}", flush=True)
    Q = np.array(Qs)
    for name, K in (("A argmax", np.array(A)), ("B firstMax", np.array(B)), ("C turnover", np.array(C))):
        K = np.array(K)
        slope = float(np.polyfit(np.log(Q), np.log(K), 1)[0])
        prod = K * Q
        cv = float(np.std(prod) / max(abs(np.mean(prod)), 1e-12))
        mono = bool(np.all(np.diff(K) < 1e-9))
        clean = mono and abs(slope + 1.0) < 0.35 and cv < 0.35
        print(f"\n{name:>11}: slope={slope:+.2f}  K*Q med={np.median(prod):.2f}  CV={cv:.2f}  "
              f"mono={mono}  [{'CLEAN ~c/Q' if clean else 'drifts'}]", flush=True)
    print("\nREAD: if B or C is CLEAN, the operational definition was the bug -> fix the gate and re-run.",
          flush=True)
    print("      if all three drift, the ~2/Q law is genuinely refuted on the FI apparatus.", flush=True)


if __name__ == "__main__":
    main()
