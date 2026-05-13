"""HWHM_ratio(g_eff) closed-form mapping via 2-level EP-rotation eigenvector
algebra + bare 4-mode floor + H_B-mixed octic lift (Phase 3, 2026-05-13).

The 2-level EP-rotation per F86 c=2 gives:
  lambda_+/-(k) = -4*gamma_0*k +/- sqrt(4*gamma_0^2 - (J*g_eff)^2)
The eigenvector rotation tan(theta) = J*g_eff / (2*gamma_0) = Q*g_eff/2,
so the K(Q) observable depends on Q only via x = Q*g_eff/2.

Bare 4-mode contribution gives HWHM_left/Q_peak = 0.671535 floor (the
BareDoubledPtfHwhmRatio constant). H_B-mixed octic residual lifts this
to the empirical class anchor.

Empirical residual per sub-class (lift above floor):
  Endpoint: 0.7728 - 0.671535 = 0.101265 (~10% lift)
  Flanking/Mid Interior: 0.7506 - 0.671535 = 0.079065 (~8% lift)
  Central self-paired: 0.7449 - 0.671535 = 0.073365 (~7% lift)
  Orbit2 escape: 0.9031 - 0.671535 = 0.231565 (~23% lift, at large Q_peak)
  Central escape orbit3: 0.5778 - 0.671535 = -0.093735 (BELOW floor, dramatic)

Hypothesis: lift scales monotonically with g_eff. Test by linear regression
(or power-law fit) on (g_eff, lift) per sub-class.
"""
import numpy as np
from collections import defaultdict

BARE_FLOOR = 0.671535
BARE_DOUBLED_PTF_XPEAK = 4.39382

# Same anchors as Phase 1
EMPIRICAL_ANCHORS = [
    (5, 0, 2.40, 0.7700), (5, 1, 1.50, 0.7454), (5, 2, 1.50, 0.7454), (5, 3, 2.40, 0.7700),
    (6, 0, 2.52, 0.7737), (6, 1, 1.65, 0.7503), (6, 2, 1.44, 0.7449),
    (6, 3, 1.65, 0.7503), (6, 4, 2.52, 0.7737),
    (7, 0, 2.53, 0.7738), (7, 1, 7.27, 0.9162), (7, 2, 1.54, 0.7469),
    (7, 3, 1.54, 0.7469), (7, 4, 7.27, 0.9162), (7, 5, 2.53, 0.7738),
    (8, 0, 2.53, 0.7734), (8, 1, 8.07, 0.8899), (8, 2, 1.51, 0.7475),
    (8, 3, 16.79, 0.5778), (8, 4, 1.51, 0.7475), (8, 5, 8.07, 0.8899), (8, 6, 2.53, 0.7734),
]


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


def main():
    by_sub = defaultdict(list)
    for N, b, Q, hwhm in EMPIRICAL_ANCHORS:
        sub = classify_bond(N, b, Q)
        g = BARE_DOUBLED_PTF_XPEAK / (Q + 2.0)
        lift = hwhm - BARE_FLOOR
        by_sub[sub].append((N, b, Q, hwhm, g, lift))

    print("=" * 110)
    print("HWHM_ratio = bare_floor + lift, with lift correlated to g_eff per sub-class")
    print("=" * 110)
    print(f"{'Sub-class':<22} {'<g_eff>':>9} {'<HWHM>':>9} {'<lift>':>9} "
          f"{'lift/g_eff':>11} {'lift correlation type':<30}")
    print("-" * 110)
    for sub in ["Endpoint", "Flanking", "Mid", "CentralSelfPaired",
                "Orbit2Escape", "CentralEscapeOrbit3"]:
        rows = by_sub.get(sub, [])
        if not rows:
            continue
        gs = np.array([r[4] for r in rows])
        lifts = np.array([r[5] for r in rows])
        hwhms = np.array([r[3] for r in rows])
        # check lift / g_eff ratio constant?
        ratios = lifts / gs if (gs > 0).all() else np.array([np.nan])
        ratio_mean = np.mean(ratios)
        ratio_std = np.std(ratios)
        relation = f"~constant lift/g_eff = {ratio_mean:.4f} +/- {ratio_std:.4f}"
        print(f"{sub:<22} {np.mean(gs):>9.4f} {np.mean(hwhms):>9.6f} {np.mean(lifts):>9.6f} "
              f"{ratio_mean:>11.4f} {relation:<30}")

    print()
    print("=" * 110)
    print("Hypothesis test: HWHM_ratio = bare_floor + alpha * g_eff (per sub-class)")
    print("=" * 110)
    print(f"{'Sub-class':<22} {'alpha (slope)':>13} {'beta (offset)':>13} {'R^2':>8} {'predict_at_g=1':>16}")
    print("-" * 110)
    for sub in ["Endpoint", "Flanking", "Mid"]:
        rows = by_sub.get(sub, [])
        if len(rows) < 2:
            continue
        gs = np.array([r[4] for r in rows])
        lifts = np.array([r[5] for r in rows])
        # linear fit: lift = alpha * g_eff + beta
        coefs = np.polyfit(gs, lifts, deg=1)
        alpha, beta = coefs[0], coefs[1]
        predicted = alpha * gs + beta
        ss_res = np.sum((lifts - predicted) ** 2)
        ss_tot = np.sum((lifts - np.mean(lifts)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        predict_at_1 = BARE_FLOOR + alpha * 1.0 + beta
        print(f"{sub:<22} {alpha:>13.6f} {beta:>13.6f} {r2:>8.4f} {predict_at_1:>16.6f}")


if __name__ == "__main__":
    main()
