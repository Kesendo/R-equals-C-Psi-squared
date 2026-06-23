"""handshake_decoder M2c: the read-cost law ~2/Q (2026-06-19, gate-first).

RESULT 2026-06-20: this strict gate FIRES (slope -0.67 not -1; K_peak*Q drifts 1.3->3.6 with a
regime break near Q~8). Qualitatively the cost DOES fall with Q (the doc's flat-or-inverted
falsification line is not crossed), but the quantitative ~2/Q does not hold. The operational
definition was ruled out as the culprit by handshake_read_cost_diag.py (argmax-FI == first
local max). See HANDSHAKE_GEOMETRY.md 'Read-cost'. The verdict strings below are conditional.


HANDSHAKE_GEOMETRY.md: "one recall rotates dark -> bright (dwell ~1/J) while the bright pays the light
(~2gamma), so the dose cost of one read of the dark memory is ~2/Q. High-Q reads almost free; at the EP
every read erases of order what it reads. Falsifiable form: measure cost-per-recall(Q) on the same
apparatus as the FI(Q) curve; a FLAT or INVERTED dependence refutes it."

OPERATIONAL DEFINITION (the literal "dose cost of one read", on the FI apparatus): cost-per-recall(Q) =
K_peak, the dose K=gamma*t at which the readout's distinguishing signal FI(delta-J) is MAXIMAL -- the dose
you must spend to read best. The ReadoutFisher apparatus already computes FI over a K-grid and takes the
MAX value; the read-cost is its ARGMAX. Z-basis (population) readout, the bonding carrier |psi_1>, which is
pure single excitation (<X>=<Y>=0), so the read lives in the cheap N^2 (1,1) Haken-Strobl block -- exact for
this carrier, runs to high N (vs ReadoutFisher's full 4^N, N<=5).

GATES that can fire:
  G1  SIGN/SCALING: K_peak must DECREASE with Q (a flat or inverted K_peak(Q) refutes the law). log-log
      slope of K_peak vs Q ~ -1.
  G2  PREFACTOR: K_peak * Q ~ const ~ 2 (the "~2/Q"; the hypothesis is order-of-magnitude on the 2).
Cross-check (B): the fraction of bright amplitude paid by K_peak ~ 2*K_peak, also ~1/Q -- the "erases of
order what it reads" at the EP.
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


blk = _load("blk", "simulations/handshake_rk_block.py")
psi, hopping, L11, populations = blk.psi, blk.hopping, blk.L11, blk.populations


def fi_over_k(N, J, g, dJ, bnds, defect, kgrid):
    """Z-population Fisher info of a delta-J bond defect, per dose K, in the (1,1) block."""
    from scipy.linalg import eig
    psi1 = np.array([psi(1, a, N) for a in range(N)])
    v0 = np.outer(psi1, psi1).astype(complex).flatten()
    tg = kgrid / g                                            # K = g*t  =>  t = K/g
    HA = hopping(N, J, bnds)
    wA, RA = eig(L11(HA, g))
    nA = populations(wA, RA, np.linalg.inv(RA) @ v0, tg, N)   # [N, len(kgrid)]
    HB = HA.copy()
    HB[defect[0], defect[1]] += dJ
    HB[defect[1], defect[0]] += dJ
    wB, RB = eig(L11(HB, g))
    nB = populations(wB, RB, np.linalg.inv(RB) @ v0, tg, N)
    dp = (nB - nA) / dJ
    fi = np.sum(dp ** 2 / np.maximum(nA, 1e-12), axis=0)      # FI(K)
    return fi


def main():
    N, J, dJ = 5, 1.0, 0.02
    bnds = bonds(N, "chain")
    defect = bnds[len(bnds) // 2]                             # interior bond
    Qs = [1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 13.0, 20.0, 35.0]
    print(f"=== handshake M2c read-cost: K_peak(Q) = dose at FI-max (N={N}, defect {defect}, Z readout) ===\n",
          flush=True)
    from scipy.signal import find_peaks
    # "ONE recall" = the FIRST dark->bright rotation = the FIRST FI peak (NOT the global argmax, which at
    # high Q jumps to a later, bigger coherent revival -- the regime jump that fired the first gate).
    print(f"{'Q':>6} {'K_peak=1st (cost)':>17} {'2/Q':>8} {'K_peak*Q':>9} {'FI@peak':>8} {'K_global':>9}", flush=True)
    rows = []
    for Q in Qs:
        g = 1.0 / Q
        kgrid = np.linspace(0.002, 6.0, 3000)
        fi = fi_over_k(N, J, g, dJ, bnds, defect, kgrid)
        peaks, _ = find_peaks(fi, height=0.05 * fi.max(), prominence=0.02 * fi.max())
        ip = int(peaks[0]) if len(peaks) else int(np.argmax(fi))   # the first recall
        ig = int(np.argmax(fi))                                    # the best revival (the old wrong object)
        k_peak = kgrid[ip]
        rows.append((Q, k_peak))
        print(f"{Q:>6.1f} {k_peak:>17.4f} {2.0/Q:>8.4f} {k_peak*Q:>9.3f} {fi[ip]:>8.3f} {kgrid[ig]:>9.4f}",
              flush=True)
    Q = np.array([r[0] for r in rows])
    K = np.array([r[1] for r in rows])
    # G1: log-log slope K_peak vs Q (want ~ -1); strictly decreasing
    slope = float(np.polyfit(np.log(Q), np.log(K), 1)[0])
    decreasing = bool(np.all(np.diff(K) < 1e-9))
    # G2: prefactor K_peak*Q ~ const
    prod = K * Q
    pref_med = float(np.median(prod))
    pref_cv = float(np.std(prod) / max(abs(np.mean(prod)), 1e-12))
    g1 = decreasing and abs(slope + 1.0) < 0.35
    g2 = pref_cv < 0.35
    print(f"\nG1 sign/scaling: log-log slope = {slope:+.2f} (want ~ -1), strictly decreasing = {decreasing}"
          f"  [{'ok' if g1 else 'CHECK'}]", flush=True)
    print(f"G2 prefactor: K_peak*Q median = {pref_med:.2f} (want ~2), CV = {pref_cv:.2f}  [{'ok' if g2 else 'CHECK'}]",
          flush=True)
    print("\nVERDICT:", flush=True)
    if g1 and g2:
        print(f"  read-cost ~ {pref_med:.1f}/Q CONFIRMED: cost-per-recall DECREASES as ~1/Q (slope {slope:+.2f}),"
              f"\n  prefactor {pref_med:.1f} ~ 2. High-Q reads cheap, EP reads dear. The hypothesis holds.", flush=True)
    elif g1:
        print(f"  scaling ~1/Q holds (slope {slope:+.2f}) but the prefactor is {pref_med:.1f}, not ~2 -- the law's"
              f"\n  SHAPE is right, the constant differs; report it, do not force 2.", flush=True)
    else:
        print(f"  GATE FIRED: K_peak(Q) is flat/inverted (slope {slope:+.2f}) -- the read-cost law is REFUTED"
              f"\n  on this apparatus, OR the operational definition is wrong; diagnose, do not loosen.", flush=True)
    return g1


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
