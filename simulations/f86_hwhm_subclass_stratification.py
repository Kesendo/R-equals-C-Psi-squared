"""F86 HWHM_ratio + g_eff stratified by F71-orbit bond sub-class (Phase 1, 2026-05-13).

Reads existing F90 bridge empirical anchors (per-bond HWHM_ratio + Q_peak at
N=5..8 from PROOF_F90_F86C2_BRIDGE.md numerical verification table) and
classifies bonds into sub-classes (Endpoint, Flanking, Mid, CentralSelfPaired,
Orbit2Escape, CentralEscapeOrbit3). Computes g_eff(N, b) via empirical inversion
g_eff = BareDoubledPtfXPeak / (Q_peak + 2) per PROOF_F86_QPEAK.md:90.

Output: per-bond stratified table + per-sub-class summary (count, mean, min,
max of HWHM_ratio + g_eff). Phase 1 of F86 HWHM closed-form attack.
"""
import numpy as np
from collections import defaultdict

# Empirical anchors from PROOF_F90_F86C2_BRIDGE.md "Numerical verification" tables
EMPIRICAL_ANCHORS = [
    # (N, bond_index, Q_peak, HWHM_ratio)
    (5, 0, 2.40, 0.7700), (5, 1, 1.50, 0.7454), (5, 2, 1.50, 0.7454), (5, 3, 2.40, 0.7700),
    (6, 0, 2.52, 0.7737), (6, 1, 1.65, 0.7503), (6, 2, 1.44, 0.7449),
    (6, 3, 1.65, 0.7503), (6, 4, 2.52, 0.7737),
    (7, 0, 2.53, 0.7738), (7, 1, 7.27, 0.9162), (7, 2, 1.54, 0.7469),
    (7, 3, 1.54, 0.7469), (7, 4, 7.27, 0.9162), (7, 5, 2.53, 0.7738),
    (8, 0, 2.53, 0.7734), (8, 1, 8.07, 0.8899), (8, 2, 1.51, 0.7475),
    (8, 3, 16.79, 0.5778), (8, 4, 1.51, 0.7475), (8, 5, 8.07, 0.8899), (8, 6, 2.53, 0.7734),
]

BARE_DOUBLED_PTF_XPEAK = 4.39382  # from C2EffectiveSpectrum / PROOF_F86_QPEAK


def classify_bond(N: int, b: int, Q_peak: float) -> str:
    n_bonds = N - 1
    if b == 0 or b == n_bonds - 1:
        return "Endpoint"
    if Q_peak > 12.0:
        return "CentralEscapeOrbit3"
    if Q_peak > 4.0:
        return "Orbit2Escape"
    if N % 2 == 0 and b == n_bonds // 2:
        return "CentralSelfPaired"
    if b == 1 or b == n_bonds - 2:
        return "Flanking"
    return "Mid"


def compute_g_eff(Q_peak: float) -> float:
    return BARE_DOUBLED_PTF_XPEAK / (Q_peak + 2.0)


def main():
    by_subclass = defaultdict(list)
    for N, b, Q, hwhm in EMPIRICAL_ANCHORS:
        sub = classify_bond(N, b, Q)
        g = compute_g_eff(Q)
        by_subclass[sub].append((N, b, Q, hwhm, g))

    print("=" * 100)
    print(f"Per-bond stratified anchors (N=5..8 from F90 bridge)")
    print("=" * 100)
    print(f"{'Sub-class':<22} {'N':>3} {'b':>3} {'Q_peak':>10} {'HWHM_ratio':>11} {'g_eff':>8}")
    print("-" * 100)
    for sub in ["Endpoint", "Flanking", "Mid", "CentralSelfPaired",
                "Orbit2Escape", "CentralEscapeOrbit3"]:
        for N, b, Q, hwhm, g in by_subclass.get(sub, []):
            print(f"{sub:<22} {N:>3} {b:>3} {Q:>10.4f} {hwhm:>11.6f} {g:>8.4f}")

    print()
    print("=" * 100)
    print("Per-sub-class summary")
    print("=" * 100)
    print(f"{'Sub-class':<22} {'count':>5} {'HWHM_mean':>11} {'HWHM_min':>10} {'HWHM_max':>10} {'g_eff_mean':>11}")
    print("-" * 100)
    for sub in ["Endpoint", "Flanking", "Mid", "CentralSelfPaired",
                "Orbit2Escape", "CentralEscapeOrbit3"]:
        rows = by_subclass.get(sub, [])
        if not rows:
            continue
        hwhms = [r[3] for r in rows]
        gs = [r[4] for r in rows]
        print(f"{sub:<22} {len(rows):>5} {np.mean(hwhms):>11.6f} {np.min(hwhms):>10.6f} "
              f"{np.max(hwhms):>10.6f} {np.mean(gs):>11.6f}")


if __name__ == "__main__":
    main()
