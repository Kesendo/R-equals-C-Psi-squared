#!/usr/bin/env python3
"""Diagnostic for the c_1 alpha-fit outlier (systematic-debugging Phase 1).

The c_1-via-alpha-fit pipeline produces sporadic wild outliers (reference
c1_bond_scan_multi_N at N=4 uniform: c_1(|S_2>, bond1) ~ +213, bond2 ~ -403,
own F71 mirror check fails with diff 4e+02). This script instruments ONE
outlier cell and one working cell to localise WHY fit_alpha returns a bad
alpha -- flat purity (C1), multimodal mse landscape (C3), or shape mismatch
(C2).

Repro: N=4, uniform J=1.0, bond 1 perturbed by +/-dJ.
  |S_2> -> c_1 ~ +213  (outlier)
  |S_1> -> c_1 ~ -0.13 (working)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.interpolate import interp1d

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (  # noqa: E402
    GAMMA_0, T_FINAL, N_STEPS, T_FIT_MAX,
    build_H_XY, build_liouvillian_matrix, density_matrix,
    per_site_purity, fit_alpha,
)
from c1_bilinearity_test import dicke_state  # noqa: E402

DJ = 0.01
N = 4
BOND = 1
FLAT_TOL = 1e-9   # per-site purity span below this -> site carries no closure signal


def eig_and_inv(L):
    ev, VR = eig(L)
    return ev, VR, np.linalg.inv(VR)


def propagate(dec, rho0, times):
    ev, VR, VLi = dec
    c0 = VLi @ rho0.flatten(order='F')
    d = rho0.shape[0]
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        out[k] = (VR @ (np.exp(ev * t) * c0)).reshape(d, d, order='F')
    return out


def mse_landscape(t, pA_i, pB_i, alphas, t_max=T_FIT_MAX):
    """mse(alpha) exactly as fit_alpha computes it, on an alpha grid."""
    interp = interp1d(t, pA_i, bounds_error=False,
                      fill_value=(float(pA_i[0]), float(pA_i[-1])), kind='cubic')
    mask = t <= t_max
    te, b = t[mask], pB_i[mask]
    return np.array([float(np.mean((interp(a * te) - b) ** 2)) for a in alphas])


def main():
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    mask = times <= T_FIT_MAX
    J = [1.0] * (N - 1)
    print(f"N={N} uniform J=1.0 bond={BOND} dJ={DJ} gamma_0={GAMMA_0}")
    print(f"fit window t<={T_FIT_MAX}: {int(mask.sum())} of {len(times)} steps\n")

    decA = eig_and_inv(build_liouvillian_matrix(build_H_XY(J, N), GAMMA_0, N))
    Jp = list(J); Jp[BOND] += DJ
    Jm = list(J); Jm[BOND] -= DJ
    decBp = eig_and_inv(build_liouvillian_matrix(build_H_XY(Jp, N), GAMMA_0, N))
    decBm = eig_and_inv(build_liouvillian_matrix(build_H_XY(Jm, N), GAMMA_0, N))

    for n in (1, 2):
        rho0 = density_matrix(dicke_state(N, n))
        PA = per_site_purity(propagate(decA, rho0, times), N)
        PBp = per_site_purity(propagate(decBp, rho0, times), N)
        PBm = per_site_purity(propagate(decBm, rho0, times), N)
        tag = "OUTLIER expected" if n == 2 else "working"
        print(f"=== Dicke |S_{n}>  ({tag}) ===")
        ln_p = ln_m = 0.0
        alphas = []
        for i in range(N):
            spanA = float(np.ptp(PA[mask, i]))
            spanBp = float(np.ptp(PBp[mask, i]))
            spanBm = float(np.ptp(PBm[mask, i]))
            ap, rp = fit_alpha(times, PA[:, i], PBp[:, i])
            am, rm = fit_alpha(times, PA[:, i], PBm[:, i])
            ln_p += np.log(ap); ln_m += np.log(am)
            alphas.append((ap, am))
            dd = np.diff(PA[mask, i])
            dd = dd[np.abs(dd) > 1e-12]
            n_extrema = int(np.sum(np.diff(np.sign(dd)) != 0)) if len(dd) > 1 else 0
            at_bound = " BOUND!" if (min(ap, am) < 0.105 or max(ap, am) > 9.9) else ""
            print(f"  site {i}: spanA={spanA:.2e} spanB+={spanBp:.2e} "
                  f"spanB-={spanBm:.2e}  a+={ap:.5f}(rmse {rp:.1e}) "
                  f"a-={am:.5f}(rmse {rm:.1e})  P_A extrema={n_extrema}  "
                  f"dln_a={np.log(ap) - np.log(am):+.4f}{at_bound}")
        c1 = (ln_p - ln_m) / (2 * DJ)
        print(f"  closure+={ln_p:+.6f}  closure-={ln_m:+.6f}  c_1={c1:+.4f}  "
              f"(UNGUARDED)")

        # --- Fix 1 test: flat-site guard (span_A < FLAT_TOL -> alpha := 1) ---
        ln_p_g = ln_m_g = 0.0
        n_guarded = 0
        for i in range(N):
            if float(np.ptp(PA[mask, i])) < FLAT_TOL:
                n_guarded += 1
                continue                       # alpha = 1, ln(1) = 0 contribution
            ap, _ = fit_alpha(times, PA[:, i], PBp[:, i])
            am, _ = fit_alpha(times, PA[:, i], PBm[:, i])
            ln_p_g += np.log(ap); ln_m_g += np.log(am)
        c1_g = (ln_p_g - ln_m_g) / (2 * DJ)
        print(f"  GUARDED: {n_guarded}/{N} flat sites set to alpha=1  ->  "
              f"c_1={c1_g:+.4f}")

        culprit = max(range(N),
                      key=lambda i: abs(np.log(alphas[i][0]) - np.log(alphas[i][1])))
        ag = np.linspace(0.1, 10.0, 100)
        msp = mse_landscape(times, PA[:, culprit], PBp[:, culprit], ag)
        msm = mse_landscape(times, PA[:, culprit], PBm[:, culprit], ag)
        print(f"  culprit site = {culprit}")
        sample = np.linspace(0, mask.sum() - 1, 9).astype(int)
        print(f"    P_A[site {culprit}](t) over window: "
              + " ".join(f"{PA[mask, culprit][k]:.4f}" for k in sample))
        print(f"    mse(alpha) landscape, 12 samples (argmin = best-fit alpha):")
        idx = np.linspace(0, len(ag) - 1, 12).astype(int)
        amin_p, amin_m = int(np.argmin(msp)), int(np.argmin(msm))
        for j in idx:
            mp = " <-min+" if j == amin_p else ""
            mm = " <-min-" if j == amin_m else ""
            print(f"      a={ag[j]:5.2f}  mse+={msp[j]:.3e}{mp:8s}  "
                  f"mse-={msm[j]:.3e}{mm}")
        print(f"    grid-argmin: alpha+={ag[amin_p]:.3f}  alpha-={ag[amin_m]:.3f}"
              f"   (fit_alpha returned a+={alphas[culprit][0]:.3f} "
              f"a-={alphas[culprit][1]:.3f})")
        # count local minima of mse+ -> multimodality test
        loc_min = int(np.sum((msp[1:-1] < msp[:-2]) & (msp[1:-1] < msp[2:])))
        print(f"    mse+ has {loc_min} interior local minima "
              f"(>1 => multimodal landscape)\n")


if __name__ == "__main__":
    main()
