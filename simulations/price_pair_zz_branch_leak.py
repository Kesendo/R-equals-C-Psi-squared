"""How much of a coherent-ZZ mask leaks into the price-pair c-tomograph's c02.

The price-pair runner (external pipeline, AIEvolution.UI/experiments/ibm_quantum_tomography/
run_price_pair.py, `_fit_branch`) fits each pair branch with a FREE Gaussian term,
y = A exp(-G t - (K t)^2) cos(2 pi f t + phi), and reads c_ij = (G+ - G-)/4 from the exponential
parts only. In |+++> the outer pair (0, 2) shares its neighbour q1, so under H = sum (pi zeta_b / 2) Z_i Z_j
the aligned (sum) branch carries the mask cos(pi (zeta01 + zeta12) t) and the anti-aligned (diff)
branch cos(pi (zeta01 - zeta12) t). This script fits noiseless branches with that same model on
run 3's delay grid, with the runner's per-point sigma at 8192 shots, and prints the apparent c02
the mask alone produces (no covariance in the model).

Result: the mask goes almost entirely into K+; the leak into G is small and negative, about
-0.0005/us at the run-4 shifts, against run 3's c02 sigma of about 0.0004-0.0005. A ZZ mask
therefore cannot produce run 2's c02 = +0.0025 (wrong sign, five times too large).
"""
import numpy as np
from scipy.optimize import curve_fit

TS = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0, 16.0, 22.0, 30.0, 45.0, 65.0])
SHOTS = 8192
F0, F2 = 0.0331, 0.0817          # run 3's measured single fringes on q0, q2 (MHz)
GAMMA = 0.0097                   # gamma_0 + gamma_2 (1/us), run 3's singles


def fit_branch(ts, ys, sig, f, phi, seed):
    """The runner's _fit_branch: frozen frequency and phase, free A, G, K >= 0."""
    def model(t, a, g, k):
        return a * np.exp(-g * t - (k * t) ** 2) * np.cos(2 * np.pi * f * t + phi)
    p0 = [max(min(abs(ys[0]) * 2, 1.0), 0.1), max(seed, 1e-4), 1e-4]
    popt, _ = curve_fit(model, ts, ys, p0=p0, sigma=sig, absolute_sigma=True,
                        bounds=([0.02, 0.0, 0.0], [1.5, 1.0, 0.5]), maxfev=20000)
    return popt


def apparent_c02(zeta01_khz, zeta12_khz):
    zs, zd = (zeta01_khz + zeta12_khz) * 1e-3, (zeta01_khz - zeta12_khz) * 1e-3
    plus = 0.5 * np.exp(-GAMMA * TS) * np.cos(np.pi * zs * TS) * np.cos(2 * np.pi * (F0 + F2) * TS)
    minus = 0.5 * np.exp(-GAMMA * TS) * np.cos(np.pi * zd * TS) * np.cos(2 * np.pi * (F0 - F2) * TS)
    sig = np.sqrt(2) * np.sqrt(np.clip(1 - (4 * plus) ** 2, 0.05, 1.0) / SHOTS) / 4
    gp, gm = fit_branch(TS, plus, sig, F0 + F2, 0.0, GAMMA), fit_branch(TS, minus, sig, F0 - F2, 0.0, GAMMA)
    return gp, gm, (gp[1] - gm[1]) / 4


if __name__ == "__main__":
    for z01, z12 in [(-3.918, -3.637), (-4.2, -4.2), (-3.0, -3.0), (0.0, 0.0)]:
        gp, gm, c = apparent_c02(z01, z12)
        print(f"zeta01 = {z01:+.3f}, zeta12 = {z12:+.3f} kHz: G+ = {gp[1]:.5f}, K+ = {gp[2]:.5f}, "
              f"G- = {gm[1]:.5f}, K- = {gm[2]:.5f}, apparent c02 = {c:+.5f} /us")
    _, _, c_null = apparent_c02(0.0, 0.0)
    assert abs(c_null) < 1e-6, "no mask must give no apparent covariance"
    _, _, c_run4 = apparent_c02(-3.918, -3.637)
    assert -0.001 < c_run4 < 0.0, "the run-4 mask leaks a small negative c02"
    print("checks: zero mask gives zero; the run-4 mask leaks a small negative c02")
