"""
Primordial Gamma: Collapse Axis Reanalysis
============================================

Chat-side reanalysis of the V1 probe (commit fe35f46, script
primordial_gamma_layer_3qubit.py, data results/primordial_gamma/probe1_results.txt).

V1 plotted gamma_eff/gamma_B vs Q_MB = J_MB/gamma_B and concluded "no collapse,
40x variation at fixed Q_MB". This script tests the alternative reading: the
correct dimensionless axis is the structural ratio (J_SM/J_MB), not Q_MB.

Two parts:

Part A: Re-projection of V1 data (J_SM=1.0 fixed, sweep J_MB and gamma_B).
  At fixed J_MB across gamma_B variation, does gamma_eff/gamma_B stay constant?
  (Claim: yes, in the high-R^2 regime.)
  In the strong-J_MB limit, does gamma_eff/gamma_B approach (J_SM/J_MB)^2?
  (Claim: yes, asymptotically.)

Part B: New small sweep to verify the axis: vary J_SM at fixed (J_MB, gamma_B).
  If the collapse axis is really (J_SM/J_MB), then varying J_SM should produce
  gamma_eff that scales as J_SM^2 (in the strong-J_MB regime).

Date: 2026-04-15
Authors: Tom and Claude (chat)
Context: Reaction to RESULT_PRIMORDIAL_GAMMA_OPERATIONAL.md, which interpreted
the V1 sweep as falsification. The V1 task was over-specified by chat to
test a single dimensionless axis (Q_MB). The data itself is clean; the
projection in V1 was not.
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import curve_fit
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def liouvillian(H, jumps):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def partial_trace_MB(rho_8x8):
    """Trace out qubits 1 (M) and 2 (B), keep qubit 0 (S)."""
    r = rho_8x8.reshape(2, 2, 2, 2, 2, 2)
    return np.einsum('ibjcjc->ib', r)


def extract_gamma_eff(L_super, rho0_3q, t_max, n_points=200):
    """Same extraction as V1: |rho_S_{01}(t)| ~ A exp(-2 gamma_eff t),
    fit to first 1/3 of trajectory."""
    ts = np.linspace(0, t_max, n_points + 1)
    coh = []
    for t in ts:
        v = rho0_3q.reshape(-1, order='F')
        vt = expm(L_super * t) @ v
        rho_t = vt.reshape((8, 8), order='F')
        rho_S = partial_trace_MB(rho_t)
        coh.append(abs(rho_S[0, 1]))
    coh = np.array(coh)

    # Detect rebound
    has_rebound = False
    rebound_amp = 0.0
    min_coh = coh[0]
    for c in coh[1:]:
        if c < min_coh:
            min_coh = c
        elif c > min_coh + 0.005:
            has_rebound = True
            rebound_amp = max(rebound_amp, c - min_coh)

    n_fit = max(10, n_points // 3)
    ts_fit = ts[:n_fit]
    coh_fit = coh[:n_fit]
    valid = coh_fit > 1e-15
    if np.sum(valid) < 3:
        return 0.0, 0.0, has_rebound, rebound_amp

    try:
        def decay_model(t, A, gamma_eff):
            return A * np.exp(-2 * gamma_eff * t)
        popt, _ = curve_fit(decay_model, ts_fit[valid], coh_fit[valid],
                            p0=[coh[0], 0.01], bounds=([0, 0], [1, 100]))
        gamma_eff = popt[1]
        predicted = decay_model(ts_fit[valid], *popt)
        ss_res = np.sum((coh_fit[valid] - predicted) ** 2)
        ss_tot = np.sum((coh_fit[valid] - np.mean(coh_fit[valid])) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    except Exception:
        log_coh = np.log(coh_fit[valid])
        slope = np.polyfit(ts_fit[valid], log_coh, 1)[0]
        gamma_eff = -slope / 2
        r2 = 0.5
    return gamma_eff, r2, has_rebound, rebound_amp


def run_3qubit(J_SM, J_MB, gamma_B, t_max=None):
    """Build N=3 S-M-B system, return gamma_eff, R^2, rebound info."""
    H = (J_SM * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2))
         + J_MB * 0.5 * (kron(I2, X, X) + kron(I2, Y, Y)))
    jumps = [np.sqrt(gamma_B) * kron(I2, I2, Z)]
    L_super = liouvillian(H, jumps)
    rho_plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
    rho_mix = 0.5 * I2
    rho0 = np.kron(rho_plus, np.kron(rho_mix, rho_mix))
    if t_max is None:
        # Heuristic: t_max scales inversely with the slower of two timescales
        t_max = max(50.0, 5.0 / max(gamma_B, 1e-6))
    return extract_gamma_eff(L_super, rho0, t_max)


# ============================================================
# Part A: Re-projection of V1 data
# ============================================================
print("=" * 78)
print("Part A: Re-projection of V1 data on (J_SM/J_MB) axis")
print("=" * 78)
print()
print("V1 data (J_SM = 1.0 fixed, sweep J_MB and gamma_B).")
print("Hypothesis to test: at fixed J_MB, gamma_eff/gamma_B is constant")
print("across gamma_B variation. Asymptotic form: -> (J_SM/J_MB)^2.")
print()

# V1 data hardcoded (from probe1_results.txt)
v1_data = [
    # (J_MB, gamma_B, gamma_eff/gamma_B, R^2)
    (0.01, 0.01, 0.000058, -0.0000),
    (0.01, 0.03, 0.000085,  0.0004),
    (0.01, 0.10, 0.000069,  0.0049),
    (0.01, 0.30, 0.000084,  0.0222),
    (0.01, 1.00, 0.000096,  0.1126),
    (0.03, 0.01, 0.000515,  0.0001),
    (0.03, 0.03, 0.000574,  0.0024),
    (0.03, 0.10, 0.000504,  0.0271),
    (0.03, 0.30, 0.000478,  0.0725),
    (0.03, 1.00, 0.000349,  0.1557),
    (0.10, 0.01, 0.008557,  0.0028),
    (0.10, 0.03, 0.006191,  0.0208),
    (0.10, 0.10, 0.005209,  0.1864),
    (0.10, 0.30, 0.004204,  0.3225),
    (0.10, 1.00, 0.001897,  0.3117),
    (0.30, 0.01, 0.079234,  0.0163),
    (0.30, 0.03, 0.051228,  0.0886),
    (0.30, 0.10, 0.042780,  0.4991),
    (0.30, 0.30, 0.035103,  0.7046),
    (0.30, 1.00, 0.012096,  0.6193),
    (1.00, 0.01, 0.435337,  0.7682),
    (1.00, 0.03, 0.414734,  0.9261),
    (1.00, 0.10, 0.388501,  0.9597),
    (1.00, 0.30, 0.397457,  0.9248),
    (1.00, 1.00, 0.135722,  0.8814),
    (3.00, 0.01, 0.080710,  0.0706),
    (3.00, 0.03, 0.096635,  0.5635),
    (3.00, 0.10, 0.098643,  0.9455),
    (3.00, 0.30, 0.099554,  0.9820),
    (3.00, 1.00, 0.100524,  0.9721),
    (10.0, 0.01, 0.009641,  0.0069),
    (10.0, 0.03, 0.009473,  0.1127),
    (10.0, 0.10, 0.009641,  0.7907),
    (10.0, 0.30, 0.009673,  0.9184),
    (10.0, 1.00, 0.009672,  0.9170),
]

J_SM_fixed = 1.0

print(f"  {'J_MB':>5}  {'(J_SM/J_MB)^2':>14}  {'rows w/ R^2>0.9':>18}")
print(f"  {' ':>5}  {'(prediction)':>14}  {'g_eff/g_B values':>18}")
print("  " + "-" * 60)
J_MB_values = sorted(set(d[0] for d in v1_data))
for J_MB in J_MB_values:
    pred = (J_SM_fixed / J_MB) ** 2
    rows_high_r2 = [(gB, ratio, r2) for (J, gB, ratio, r2) in v1_data
                    if J == J_MB and r2 > 0.9]
    if rows_high_r2:
        ratios = [r[1] for r in rows_high_r2]
        n = len(rows_high_r2)
        mean = np.mean(ratios)
        spread = (max(ratios) - min(ratios)) / mean if mean > 0 else 0
        ratios_str = ", ".join(f"{r:.4f}" for r in ratios)
        match = mean / pred
        print(f"  {J_MB:>5.2f}  {pred:>14.4f}  n={n}, mean={mean:.4f}, spread={spread:.1%}, mean/pred={match:.3f}")
        print(f"  {' ':>5}  {' ':>14}    [{ratios_str}]")
    else:
        print(f"  {J_MB:>5.2f}  {pred:>14.4f}  (no rows with R^2 > 0.9)")
print()
print("Interpretation:")
print("  - At J_MB=10, mean(g_eff/g_B) = 0.0097, prediction (1/10)^2 = 0.01.")
print("    Match within 3%, spread across 100x gamma_B variation < 3%.")
print("  - At J_MB=3, mean ~ 0.099, prediction ~ 0.111. Match ~ 90%.")
print("  - At J_MB=1, mean ~ 0.40, prediction = 1.0. Match ~ 40%.")
print("    We are at J_MB ~ J_SM here, no longer in strong-J_MB asymptotic regime.")


# ============================================================
# Part B: New sweep, vary J_SM at fixed (J_MB, gamma_B)
# ============================================================
print()
print("=" * 78)
print("Part B: New sweep, vary J_SM at fixed J_MB, gamma_B")
print("=" * 78)
print()
print("If the collapse axis is (J_SM/J_MB), then varying J_SM at fixed J_MB")
print("should give gamma_eff propto J_SM^2 (in strong-J_MB regime).")
print()

# Use J_MB=10, gamma_B=0.1: this is deep in the asymptotic regime in V1.
J_MB = 10.0
gamma_B = 0.1

J_SM_sweep = [0.1, 0.3, 1.0, 3.0]
print(f"  Fixed: J_MB = {J_MB}, gamma_B = {gamma_B}")
print(f"  Sweep: J_SM in {J_SM_sweep}")
print()
print(f"  {'J_SM':>5}  {'gamma_eff':>11}  {'g_eff/g_B':>11}  {'(J_SM/J_MB)^2':>14}  {'ratio':>7}  {'R^2':>6}")
print("  " + "-" * 68)

new_data = []
for J_SM in J_SM_sweep:
    gamma_eff, r2, has_rebound, reb_amp = run_3qubit(J_SM, J_MB, gamma_B)
    ratio = gamma_eff / gamma_B
    pred = (J_SM / J_MB) ** 2
    match = ratio / pred if pred > 0 else 0
    new_data.append((J_SM, gamma_eff, ratio, pred, match, r2))
    print(f"  {J_SM:>5.2f}  {gamma_eff:>11.6f}  {ratio:>11.6f}  {pred:>14.6f}  {match:>7.3f}  {r2:>6.3f}")

print()
print("Interpretation:")
print("  If gamma_eff/gamma_B ~ (J_SM/J_MB)^2 holds, the ratio column should")
print("  be ~ 1.0 across the J_SM sweep (within fit-quality tolerance).")
print()

# Cross-check: gamma_eff should scale as J_SM^2
print("Cross-check (gamma_eff scaling with J_SM):")
print(f"  {'J_SM_pair':>10}  {'g_eff ratio':>12}  {'J_SM^2 ratio':>13}  {'match':>7}")
print("  " + "-" * 50)
for i in range(1, len(new_data)):
    J_SM_a, g_a = new_data[i-1][0], new_data[i-1][1]
    J_SM_b, g_b = new_data[i][0], new_data[i][1]
    g_ratio = g_b / g_a
    J_ratio_sq = (J_SM_b / J_SM_a) ** 2
    match = g_ratio / J_ratio_sq
    print(f"  {J_SM_a:.2f}->{J_SM_b:.2f}  {g_ratio:>12.4f}  {J_ratio_sq:>13.4f}  {match:>7.3f}")


# ============================================================
# Part C: Two-axis dimensional analysis check
# ============================================================
print()
print("=" * 78)
print("Part C: Test g_eff/g_B = h(J_SM/g_B, J_MB/g_B) collapse")
print("=" * 78)
print()
print("Buckingham Pi: with units 1/time for J_SM, J_MB, gamma_B,")
print("the dimensionless function must be h(J_SM/g_B, J_MB/g_B).")
print("Test: pick two configurations with the same (J_SM/g_B, J_MB/g_B) ratios")
print("and check that g_eff/g_B agrees.")
print()

# Two configurations with the same dimensionless ratios:
# Config 1: gamma_B = 0.1, J_SM = 1.0, J_MB = 5.0  -> (10, 50)
# Config 2: gamma_B = 0.5, J_SM = 5.0, J_MB = 25.0 -> (10, 50)
configs = [
    (1.0, 5.0, 0.1),
    (5.0, 25.0, 0.5),
    (0.5, 2.5, 0.05),
]
print(f"  All configurations have J_SM/g_B = 10, J_MB/g_B = 50.")
print(f"  Predicted g_eff/g_B should be identical (Buckingham Pi).")
print()
print(f"  {'J_SM':>5}  {'J_MB':>5}  {'gamma_B':>8}  {'gamma_eff':>11}  {'g_eff/g_B':>11}  {'R^2':>6}")
print("  " + "-" * 60)
ratios_pi = []
for J_SM, J_MB, gB in configs:
    g_eff, r2, _, _ = run_3qubit(J_SM, J_MB, gB)
    ratio = g_eff / gB
    ratios_pi.append(ratio)
    print(f"  {J_SM:>5.2f}  {J_MB:>5.2f}  {gB:>8.4f}  {g_eff:>11.6f}  {ratio:>11.6f}  {r2:>6.3f}")

print()
spread = (max(ratios_pi) - min(ratios_pi)) / np.mean(ratios_pi)
print(f"  Spread across the three: {spread:.2%}")
print(f"  Mean g_eff/g_B = {np.mean(ratios_pi):.6f}")
print(f"  Predicted from (J_SM/J_MB)^2 = {(1.0/5.0)**2:.6f}")
print()

# ============================================================
# Conclusion (honest)
# ============================================================
print()
print("=" * 78)
print("CONCLUSION (honest read of the data)")
print("=" * 78)
print()
print("Part A (re-projection of V1 data):")
print("  STRONG support for collapse. At fixed J_MB the ratio g_eff/g_B is")
print("  approximately constant across gamma_B variation in the high-R^2 rows.")
print("  At J_MB=10: spread 0%, match to (J_SM/J_MB)^2 within 3%.")
print("  At J_MB=3:  spread 2%, match within 10%.")
print("  At J_MB=1:  spread 7%, match 40% (we are at J_MB ~ J_SM, outside")
print("              the strong-J_MB asymptotic regime where (J_SM/J_MB)^2 holds).")
print()
print("Part B (vary J_SM at fixed J_MB=10):")
print("  PARTIAL support. The J_SM=1->3 cross-check matches within 1.3%.")
print("  Smaller J_SM values give R^2 < 0.1: the Markovian fit fails because")
print("  the decay timescale (~ J_SM^2/J_MB^2 * gamma_B) becomes comparable to")
print("  or longer than the chosen fit window. This is a method limitation,")
print("  not necessarily a hypothesis failure. A wider time window or a")
print("  different extraction method would be needed to test these regimes.")
print()
print("Part C (Buckingham Pi check):")
print("  PARTIAL support. Two of three identical-dimensionless-ratio configs")
print("  give g_eff/g_B = 0.0388 to four digits; the third (with R^2=0.97)")
print("  gives 0.0329, ~16% off. Either the extraction method has")
print("  scale-dependent artifacts (different fit windows relative to decay")
print("  scale at different absolute scales), or the dimensional analysis")
print("  predicts only an approximate scale invariance and there are O(1)")
print("  scale-dependent corrections at finite parameter values.")
print()
print("Net verdict:")
print("  The refractive composition g_eff = g_B * f(J_SM/J_MB) is supported")
print("  in the strong-J_MB asymptotic regime (where the Markovian fit is")
print("  reliable). Asymptotic form: f -> (J_SM/J_MB)^2. The V1 verdict")
print("  ('no operational support') overlooked this collapse axis. Outside")
print("  the asymptotic regime, the function f either has corrections that")
print("  break clean dimensional scaling, or the Markovian extraction method")
print("  is inadequate. Both possibilities are open and worth investigating.")
print()
print("The hypothesis is NOT cleanly confirmed and NOT cleanly falsified.")
print("It survives in a sharpened, scope-limited form: refractive composition")
print("with (J_SM/J_MB) as the structural axis, asymptotically (J_SM/J_MB)^2,")
print("with non-trivial corrections at finite J_MB/J_SM that need further work.")
