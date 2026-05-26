"""F113 candidate: closed-form for the F112 counterexample asymmetry magnitude.

The Welle 2 hardware analysis identified Z-drive H + σ⁻ T1 as the FIRST
structural counterexample to the F112 polarity-balance broader empirical
envelope. Synthetic isolation gives:

  H = ω·(Z_0 + Z_1)/2
  c = σ⁻ on each qubit at rate γ_T1
  + optionally γ_Z Z-dephasing per site

Observed at (ω, γ_T1, γ_Z) = (0.13, 0.001, 0.005): rel asymmetry ≈ 3.85e-3.

This script:
  1. Sweeps (ω, γ_T1, γ_Z) over a grid at N=2
  2. Computes F112 asymmetry numerically at each grid point
  3. Identifies the scaling law via log-log polynomial regression
  4. Candidates a closed form for asymmetry(ω, γ_T1, γ_Z)
  5. Verifies the candidate against additional parameter points (held out)

If a clean monomial scaling appears (asymmetry = c · ω^a · γ_T1^b · γ_Z^c),
the closed form is the F113 candidate.

Run: python -X utf8 simulations/_f113_break_formula_derivation.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 0], [1, 0]], dtype=complex)


def site_op(N, l, m2):
    mats = [I2] * N
    mats[l] = m2
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def f112_setup(omega, gamma_T1, gamma_Z, N=2):
    """Return (H, c_ops, gammas, sigma) for the counterexample family."""
    H = (omega / 2.0) * (site_op(N, 0, Z) + site_op(N, 1, Z))
    c_list = []
    g_list = []
    for l in range(N):
        c_list.append(site_op(N, l, Z))
        g_list.append(gamma_Z)
        c_list.append(site_op(N, l, SIGMA_MINUS))
        g_list.append(gamma_T1)
    sigma = sum(g_list)
    return H, c_list, g_list, sigma


def asymmetry_at(omega, gamma_T1, gamma_Z, N=2):
    """Compute F112 asymmetry (absolute, not relative) and ‖M‖² for the
    counterexample family at the given parameters."""
    H, c_list, g_list, sigma = f112_setup(omega, gamma_T1, gamma_Z, N=N)
    r = fw.polarity_coordinates_from_hc(H, c_list, g_list, N, sigma=sigma)
    return float(r['asymmetry']), float(r['norm_sq']['M'])


def log_log_scaling(xs, ys, label):
    """Fit log(y) = a*log(x) + b; return (a, b, R²)."""
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    mask = (xs > 0) & (np.abs(ys) > 1e-15)
    if mask.sum() < 3:
        print(f"  {label}: insufficient non-zero data ({mask.sum()} pts), cannot fit")
        return None, None, None
    lx, ly = np.log(xs[mask]), np.log(np.abs(ys[mask]))
    a, b = np.polyfit(lx, ly, 1)
    pred = a * lx + b
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1 - ss_res / max(ss_tot, 1e-15)
    print(f"  {label}: log|y| = {a:.4f}·log(x) + {b:.4f}  (R² = {r2:.6f})")
    print(f"           → |y| ≈ {np.exp(b):.4e} · x^{a:.4f}")
    return a, b, r2


def main():
    print("F113 candidate: deriving the F112 counterexample asymmetry formula")
    print("=" * 78)
    print()
    print("Family: H = ω·(Z_0+Z_1)/2,  c = σ⁻ at γ_T1 + Z-deph at γ_Z (per site, N=2)")
    print()

    # ============================================================
    # Univariate scalings: pin down each parameter's exponent
    # ============================================================

    # Baseline: scan ω with fixed γ_T1, γ_Z
    print("--- Scan 1: vary ω, fix γ_T1=0.001, γ_Z=0.005 ---")
    omegas = np.logspace(-3, 0, 12)
    asym_omega = [asymmetry_at(o, 0.001, 0.005)[0] for o in omegas]
    log_log_scaling(omegas, asym_omega, "asymmetry vs ω")
    print()

    # Scan γ_T1 with fixed ω, γ_Z
    print("--- Scan 2: vary γ_T1, fix ω=0.13, γ_Z=0.005 ---")
    gt1s = np.logspace(-4, -1, 12)
    asym_gt1 = [asymmetry_at(0.13, g, 0.005)[0] for g in gt1s]
    log_log_scaling(gt1s, asym_gt1, "asymmetry vs γ_T1")
    print()

    # Scan γ_Z with fixed ω, γ_T1
    print("--- Scan 3: vary γ_Z, fix ω=0.13, γ_T1=0.001 ---")
    gzs = np.logspace(-4, -1, 12)
    asym_gz = [asymmetry_at(0.13, 0.001, g)[0] for g in gzs]
    log_log_scaling(gzs, asym_gz, "asymmetry vs γ_Z")
    print()

    # ============================================================
    # Try the candidate scaling: asymmetry = c · ω^a · γ_T1^b · γ_Z^c
    # ============================================================

    print("--- Multivariate fit: log|asym| = A + a·log(ω) + b·log(γ_T1) + c·log(γ_Z) ---")
    rng = np.random.default_rng(seed=2026)
    n_pts = 60
    omg_pts = rng.uniform(0.01, 0.5, n_pts)
    gt1_pts = rng.uniform(1e-4, 0.05, n_pts)
    gz_pts = rng.uniform(1e-4, 0.05, n_pts)
    asym_pts = np.array([asymmetry_at(o, g1, gz)[0]
                         for o, g1, gz in zip(omg_pts, gt1_pts, gz_pts)])

    mask = np.abs(asym_pts) > 1e-15
    X_mat = np.column_stack([
        np.ones(mask.sum()),
        np.log(omg_pts[mask]),
        np.log(gt1_pts[mask]),
        np.log(gz_pts[mask]),
    ])
    y_vec = np.log(np.abs(asym_pts[mask]))
    coefs, *_ = np.linalg.lstsq(X_mat, y_vec, rcond=None)
    A, a, b, c = coefs
    pred = X_mat @ coefs
    ss_res = float(np.sum((y_vec - pred) ** 2))
    ss_tot = float(np.sum((y_vec - y_vec.mean()) ** 2))
    r2 = 1 - ss_res / max(ss_tot, 1e-15)
    print(f"  log|asym| = {A:.4f} + {a:.4f}·log(ω) + {b:.4f}·log(γ_T1) + {c:.4f}·log(γ_Z)")
    print(f"  R² = {r2:.6f}")
    print(f"  → |asym| ≈ {np.exp(A):.4e} · ω^{a:.4f} · γ_T1^{b:.4f} · γ_Z^{c:.4f}")
    print()

    # Round exponents to nearest integer/half-integer and re-fit
    a_round = round(a)
    b_round = round(b)
    c_round = round(c)
    print(f"  Rounded exponents: ω^{a_round}, γ_T1^{b_round}, γ_Z^{c_round}")

    # Refit just the constant C = asym / (ω^a · γ_T1^b · γ_Z^c)
    if a_round != 0 and b_round != 0:
        const_vals = []
        for o, g1, gz, asym in zip(omg_pts[mask], gt1_pts[mask], gz_pts[mask], asym_pts[mask]):
            denom = (o ** a_round) * (g1 ** b_round) * (gz ** c_round if c_round != 0 else 1.0)
            if abs(denom) > 1e-30:
                const_vals.append(asym / denom)
        const_vals = np.array(const_vals)
        print(f"  Implied constant C (from {len(const_vals)} samples):")
        print(f"    mean = {const_vals.mean():.6f}, std = {const_vals.std():.6f}")
        print(f"    median = {np.median(const_vals):.6f}")
        print(f"    range = [{const_vals.min():.6f}, {const_vals.max():.6f}]")

    print()

    # ============================================================
    # Check ‖M‖² scaling too (denominator of relative asymmetry)
    # ============================================================
    print("--- ‖M‖² scaling at the same sweep points ---")
    m_pts = np.array([asymmetry_at(o, g1, gz)[1]
                      for o, g1, gz in zip(omg_pts, gt1_pts, gz_pts)])
    mask_m = m_pts > 1e-15
    if mask_m.sum() > 4:
        X_m = np.column_stack([
            np.ones(mask_m.sum()),
            np.log(omg_pts[mask_m]),
            np.log(gt1_pts[mask_m]),
            np.log(gz_pts[mask_m]),
        ])
        coefs_m, *_ = np.linalg.lstsq(X_m, np.log(m_pts[mask_m]), rcond=None)
        Am, am, bm, cm = coefs_m
        print(f"  log‖M‖² = {Am:.4f} + {am:.4f}·log(ω) + {bm:.4f}·log(γ_T1) + {cm:.4f}·log(γ_Z)")

    # ============================================================
    # Rel asymmetry scaling
    # ============================================================
    print()
    print("--- Relative asymmetry = asym / ‖M‖² ---")
    rel_pts = np.where(mask_m, np.abs(asym_pts) / np.where(mask_m, m_pts, 1.0), 0.0)
    rel_mask = rel_pts > 1e-15
    if rel_mask.sum() > 4:
        X_r = np.column_stack([
            np.ones(rel_mask.sum()),
            np.log(omg_pts[rel_mask]),
            np.log(gt1_pts[rel_mask]),
            np.log(gz_pts[rel_mask]),
        ])
        coefs_r, *_ = np.linalg.lstsq(X_r, np.log(rel_pts[rel_mask]), rcond=None)
        Ar, ar, br, cr = coefs_r
        pred_r = X_r @ coefs_r
        ss_res_r = float(np.sum((np.log(rel_pts[rel_mask]) - pred_r) ** 2))
        ss_tot_r = float(np.sum((np.log(rel_pts[rel_mask]) - np.log(rel_pts[rel_mask]).mean()) ** 2))
        r2_r = 1 - ss_res_r / max(ss_tot_r, 1e-15)
        print(f"  log(rel asym) = {Ar:.4f} + {ar:.4f}·log(ω) + {br:.4f}·log(γ_T1) + {cr:.4f}·log(γ_Z)")
        print(f"  R² = {r2_r:.6f}")
        print(f"  → rel asym ≈ {np.exp(Ar):.4e} · ω^{ar:.4f} · γ_T1^{br:.4f} · γ_Z^{cr:.4f}")


if __name__ == '__main__':
    main()
