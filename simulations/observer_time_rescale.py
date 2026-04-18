#!/usr/bin/env python3
"""observer_time_rescale.py

Observer-time rescaling analysis for the N=7 coupling-defect overlay data.

Hypothesis under test: under a local J-perturbation (experiment B) relative
to the uniform chain (experiment A), site-resolved single-qubit purity
P_B(i, t) is well approximated by a time rescaling of the control:
    P_B(i, t) = P_A(i, alpha_i * t)
for some site-specific alpha_i. If true, alpha_i < 1 means site i
"experiences time slower" in B, alpha_i > 1 means "faster". The
dimensionless dose K = gamma * t is invariant under this rescaling
(F14-extended per site).

For each site i and each defect variant, fit:
  (1) one-parameter alpha_i via bounded scalar minimisation of
      MSE(alpha) = <[P_A(i, alpha * t) - P_B(i, t)]^2> on t in [0, T_FIT].
  (2) control two-parameter (alpha, beta) with t_i = alpha * t + beta,
      via L-BFGS-B on the same MSE.

Then report:
  * alpha_i per site and J_mod
  * residual RMS per (site, J_mod) for one- and two-parameter fits
  * log-log slope and correlation r of ln(alpha_i) vs ln(J_mod) per site
  * multiplicative conservation Sum_i ln(alpha_i) per J_mod
  * boundary-hit flag (alpha hit the bounds, meaning fit is ill-posed)

Usage:
    python observer_time_rescale.py [scan_dir]

Default scan_dir: simulations/results/n7_coupling_defect_overlay/

Expected file structure in scan_dir:
    times.npy                         (T,) time vector
    experiment_A.npz                  keys: purity (T, N), ...
    experiment_B_<label>.npz          keys: purity (T, N), ...

<label> is whatever string follows the `experiment_B_` prefix; e.g.
"0.5", "1.5", "J0.5_bond23". Each B file is fitted against the single A
control.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize, minimize_scalar

T_FIT_DEFAULT = 20.0        # fit over t in [0, T_FIT]
ALPHA_BOUNDS = (0.1, 10.0)  # scalar-fit bounds
ALPHA_BOUNDS_2P = (0.05, 20.0)  # looser for the 2-param control
BETA_BOUNDS = (-10.0, 10.0)

_BOUND_EPS = 1e-4  # how close to the bound counts as a "boundary hit"


# ---------------------------------------------------------------------------
# Fitters
# ---------------------------------------------------------------------------
def _interp_A(t, p_A):
    return interp1d(t, p_A, bounds_error=False,
                    fill_value=(float(p_A[0]), float(p_A[-1])),
                    kind='cubic')


def alpha_fit(t, p_A_i, p_B_i, t_max=T_FIT_DEFAULT,
              alpha_bounds=ALPHA_BOUNDS):
    """One-parameter time-rescale fit. Returns (alpha, rmse, boundary_hit)."""
    interp = _interp_A(t, p_A_i)
    mask = t <= t_max
    t_eval = t[mask]
    b = p_B_i[mask]

    def mse(alpha):
        d = interp(alpha * t_eval) - b
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=alpha_bounds, method='bounded',
                          options={'xatol': 1e-6})
    alpha = float(res.x)
    rmse = float(np.sqrt(res.fun))
    lo, hi = alpha_bounds
    boundary = bool(abs(alpha - lo) < _BOUND_EPS or abs(alpha - hi) < _BOUND_EPS)
    return alpha, rmse, boundary


def alpha_beta_fit(t, p_A_i, p_B_i, t_max=T_FIT_DEFAULT):
    """Two-parameter control: t_i = alpha * t + beta. Returns
    (alpha, beta, rmse, boundary_hit)."""
    interp = _interp_A(t, p_A_i)
    mask = t <= t_max
    t_eval = t[mask]
    b = p_B_i[mask]

    def mse(params):
        a, bt = params
        d = interp(a * t_eval + bt) - b
        return float(np.mean(d * d))

    # Start from the 1-parameter optimum as a warm start.
    a0, _, _ = alpha_fit(t, p_A_i, p_B_i, t_max, ALPHA_BOUNDS_2P)
    res = minimize(mse, x0=[a0, 0.0],
                   bounds=[ALPHA_BOUNDS_2P, BETA_BOUNDS],
                   method='L-BFGS-B',
                   options={'ftol': 1e-12, 'gtol': 1e-10, 'maxiter': 500})
    alpha, beta = float(res.x[0]), float(res.x[1])
    rmse = float(np.sqrt(res.fun))
    a_lo, a_hi = ALPHA_BOUNDS_2P
    b_lo, b_hi = BETA_BOUNDS
    boundary = bool(abs(alpha - a_lo) < _BOUND_EPS or
                    abs(alpha - a_hi) < _BOUND_EPS or
                    abs(beta - b_lo) < _BOUND_EPS or
                    abs(beta - b_hi) < _BOUND_EPS)
    return alpha, beta, rmse, boundary


def log_log_slope_and_r(x_vals, y_vals):
    """Slope and Pearson r of ln(y) vs ln(x)."""
    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)
    m_ok = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if m_ok.sum() < 2:
        return float('nan'), float('nan')
    lx = np.log(x[m_ok])
    ly = np.log(y[m_ok])
    dx = lx - lx.mean()
    dy = ly - ly.mean()
    var_x = float(np.sum(dx * dx))
    slope = float(np.sum(dx * dy) / var_x) if var_x > 0 else float('nan')
    denom = float(np.sqrt(var_x * float(np.sum(dy * dy))))
    r = float(np.sum(dx * dy) / denom) if denom > 0 else float('nan')
    return slope, r


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_scan(scan_dir):
    """Load times.npy, experiment_A.npz, experiment_B_*.npz. Returns
    (times, purity_A, dict label -> purity_B)."""
    d = Path(scan_dir)
    t = np.load(d / "times.npy")
    A = np.load(d / "experiment_A.npz")
    p_A = A['purity']
    Bs = {}
    for f in sorted(d.glob("experiment_B_*.npz")):
        label = f.stem.removeprefix("experiment_B_")
        Bs[label] = np.load(f)['purity']
    return t, p_A, Bs


# ---------------------------------------------------------------------------
# Label parsing (for sorting the J_mod-only scan)
# ---------------------------------------------------------------------------
def try_float(label):
    try:
        return float(label)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def analyze_scan(scan_dir, t_fit=T_FIT_DEFAULT, out_dir=None, log=None):
    if log is None:
        log = print

    t, p_A, Bs = load_scan(scan_dir)
    N_sites = p_A.shape[1]
    labels = sorted(Bs.keys(), key=lambda l: (try_float(l) is None, try_float(l) or 0, l))

    log(f"=== Scan: {scan_dir} ===")
    log(f"  {len(labels)} B-variants, N_sites = {N_sites}, T = {t[-1]:.1f}, "
        f"fit window [0, {t_fit}], n_points in window = {int((t <= t_fit).sum())}")
    log()

    rows_single = []
    rows_ab = []
    alphas_mat = np.full((len(labels), N_sites), np.nan)
    rmse_mat = np.full((len(labels), N_sites), np.nan)
    boundary_mat = np.zeros((len(labels), N_sites), dtype=bool)
    rmse2p_mat = np.full((len(labels), N_sites), np.nan)

    for row_i, label in enumerate(labels):
        p_B = Bs[label]
        for site in range(N_sites):
            a, rmse, bd = alpha_fit(t, p_A[:, site], p_B[:, site], t_max=t_fit)
            a2, b2, rmse2, bd2 = alpha_beta_fit(t, p_A[:, site], p_B[:, site],
                                                t_max=t_fit)
            alphas_mat[row_i, site] = a
            rmse_mat[row_i, site] = rmse
            boundary_mat[row_i, site] = bd
            rmse2p_mat[row_i, site] = rmse2
            rows_single.append(dict(label=label, site=site, alpha=a, rmse=rmse,
                                    boundary=int(bd)))
            rows_ab.append(dict(label=label, site=site, alpha=a2, beta=b2,
                                rmse=rmse2, boundary=int(bd2)))

    # --- Table 1: alpha per (label, site) ----------------------------
    log("--- Table 1: one-parameter alpha_i fit ---")
    header = (f"  {'J_mod / label':<16} " +
              "  ".join(f"site {s}".rjust(6) for s in range(N_sites)))
    log(header)
    for row_i, label in enumerate(labels):
        cells = []
        for s in range(N_sites):
            a = alphas_mat[row_i, s]
            flag = "*" if boundary_mat[row_i, s] else " "
            cells.append(f"{a:6.3f}{flag}")
        log(f"  {label:<16} " + "  ".join(cells))
    log("  (* = fit hit the alpha bounds)")
    log()

    # --- Table 2: rmse per (label, site) -----------------------------
    log("--- Table 2: one-parameter fit RMSE (purity units) ---")
    log(header)
    for row_i, label in enumerate(labels):
        cells = [f"{rmse_mat[row_i, s]:6.4f}" for s in range(N_sites)]
        log(f"  {label:<16} " + "  ".join(cells))
    log()

    # --- Table 3: two-parameter RMSE as control ----------------------
    log("--- Table 3: two-parameter (alpha, beta) fit RMSE ---")
    log(header)
    for row_i, label in enumerate(labels):
        cells = [f"{rmse2p_mat[row_i, s]:6.4f}" for s in range(N_sites)]
        log(f"  {label:<16} " + "  ".join(cells))
    log()

    # --- Table 4: log-log slope + r vs J_mod (numeric labels only) ---
    numeric = [(l, try_float(l)) for l in labels]
    numeric = [(l, v) for l, v in numeric if v is not None]
    if len(numeric) >= 2:
        log("--- Table 4: log-log slope / r of alpha_i vs numeric label ---")
        log(f"  (numeric labels parsed: {[l for l, v in numeric]})")
        log(f"  {'site':<6} {'slope (ln a / ln J)':>22} {'r':>8}")
        numeric_alphas = np.array([alphas_mat[labels.index(l)]
                                   for l, v in numeric])
        xs = [v for l, v in numeric]
        slopes_rs = []
        for s in range(N_sites):
            slope, r = log_log_slope_and_r(xs, numeric_alphas[:, s])
            slopes_rs.append((slope, r))
            log(f"  {s:<6} {slope:>22.4f} {r:>8.3f}")
        log()
    else:
        slopes_rs = None
        log("--- Table 4: skipped (fewer than 2 numeric labels) ---\n")

    # --- Table 5: multiplicative conservation Sum_i ln(alpha_i) ------
    log("--- Table 5: multiplicative conservation ---")
    log(f"  {'label':<16} {'Sum_i ln(alpha_i)':>22}  interpretation")
    sums = []
    for row_i, label in enumerate(labels):
        s = float(np.sum(np.log(np.clip(alphas_mat[row_i], 1e-9, None))))
        sums.append((label, s))
        interp = ("conserved" if abs(s) < 0.2 else
                  "mildly broken" if abs(s) < 1.0 else "broken")
        log(f"  {label:<16} {s:>22.3f}  {interp}")
    log()

    # --- Save outputs ------------------------------------------------
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        # Per-scan CSVs
        tag = Path(scan_dir).name
        fn1 = out_dir / f"{tag}_alpha_fit.csv"
        with open(fn1, 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=['label', 'site', 'alpha',
                                               'rmse', 'boundary'])
            w.writeheader()
            for r in rows_single:
                w.writerow(r)
        fn2 = out_dir / f"{tag}_alpha_beta_fit.csv"
        with open(fn2, 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=['label', 'site', 'alpha',
                                               'beta', 'rmse', 'boundary'])
            w.writeheader()
            for r in rows_ab:
                w.writerow(r)
        fn3 = out_dir / f"{tag}_summary.json"
        summary = dict(
            scan_dir=str(scan_dir), t_fit=t_fit, labels=labels,
            alphas=alphas_mat.tolist(), rmse=rmse_mat.tolist(),
            rmse_2p=rmse2p_mat.tolist(), boundary=boundary_mat.tolist(),
            log_log_slopes_rs=slopes_rs,
            sum_ln_alpha=[s for _, s in sums],
        )
        with open(fn3, 'w', encoding='utf-8') as fh:
            json.dump(summary, fh, indent=2)
        log(f"  -> {fn1}")
        log(f"  -> {fn2}")
        log(f"  -> {fn3}")
    log()

    return dict(
        labels=labels,
        alphas=alphas_mat,
        rmse=rmse_mat,
        rmse_2p=rmse2p_mat,
        boundary=boundary_mat,
        log_log=slopes_rs,
        sums=sums,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        default = Path(__file__).parent / "results" / "n7_coupling_defect_overlay"
        scan_dirs = [default]
    else:
        scan_dirs = [Path(p) for p in sys.argv[1:]]

    out_root = Path(__file__).parent / "results" / "observer_time_rescale"
    out_root.mkdir(parents=True, exist_ok=True)

    log_lines = []
    def log(msg=""):
        print(msg, flush=True)
        log_lines.append(msg)

    log("OBSERVER-TIME RESCALE ANALYSIS")
    log("=" * 72)
    for sd in scan_dirs:
        analyze_scan(sd, t_fit=T_FIT_DEFAULT, out_dir=out_root, log=log)
    (out_root / "run_log.txt").write_text("\n".join(log_lines) + "\n",
                                         encoding="utf-8")
    print(f"\nLog: {out_root / 'run_log.txt'}")


if __name__ == "__main__":
    main()
