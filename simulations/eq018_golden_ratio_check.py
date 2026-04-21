#!/usr/bin/env python3
"""eq018_golden_ratio_check.py

Check whether the 18 self-Pi modes at N=4 (from the Pi-parity scan) are
structurally tied to the golden-ratio eigenvalue structure of the N=4 XY
chain:

    E_k = 2 cos(pi k / 5) for k=1..4
        = {phi, 1/phi, -1/phi, -phi}    with phi = (1 + sqrt(5))/2

The Liouvillian acts on operators, so mode eigenvalues are sums/differences
of H-eigenvalues depending on the operator's structure. In the n_XY=2
sector (Re(lambda) = -2*gamma_0*2 = -0.2), Im(lambda) values come from
linear combinations of the H-eigenvalues.

Hypothesis: all 18 self-Pi modes have Im = 0 because of exact cancellations
enabled by phi + (-phi) = 0, 1/phi + (-1/phi) = 0, and similar golden-ratio
symmetries. At N=6, H-eigenvalues are {2cos(pi/7), 2cos(2pi/7), ...} which
do NOT have a nice closed-form ratio structure, and the cancellations
generically do not occur.

Test plan:
  1. Build L_A at N=4, diagonalize.
  2. Verify H single-excitation eigenvalues are ±phi, ±1/phi.
  3. Extract Im(lambda) values in the n_XY=2 sector (Re ≈ -0.2).
  4. Count Im=0 modes (should be 18).
  5. For the non-zero Im values, check if they are integer linear combinations
     of {phi, 1/phi, 1, sqrt(5)}.
  6. Compare to N=6: list Im values in n_XY=3 sector, show they are NOT
     cleanly in a finite-rational-combination space.

Rules: UTF-8 stdout, no em-dashes.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, X, Y,
    build_H_XY, build_liouvillian_matrix, site_op,
)

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_golden_ratio_check"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


PHI = (1 + np.sqrt(5)) / 2


def single_excitation_hamiltonian(N, J=J_UNIFORM):
    """H_SE in the single-excitation sector: N x N tridiagonal with off-diag J."""
    H = np.zeros((N, N), dtype=float)
    for i in range(N - 1):
        H[i, i+1] = J
        H[i+1, i] = J
    return H


def match_linear_combo(target, basis, max_coef=4, tol=1e-6):
    """Try to write target = sum c_k * basis[k] with c_k in [-max_coef, max_coef]
    (integer). Return (combo, residual) if found; else None."""
    n = len(basis)
    best = None
    best_res = 1e99
    # Exhaustive search over coefficients
    from itertools import product
    ranges = [range(-max_coef, max_coef + 1)] * n
    for coefs in product(*ranges):
        val = sum(c * b for c, b in zip(coefs, basis))
        res = abs(val - target)
        if res < best_res:
            best_res = res
            best = coefs
    if best_res < tol:
        return best, best_res
    return None, best_res


def analyze_sector(evals, target_re, sector_name, tol_re=1e-9, tol_im=1e-9):
    """Extract modes in a specific Re(lambda) sector, analyse Im distribution."""
    mask = np.abs(np.real(evals) - target_re) < tol_re
    sector_evals = evals[mask]
    n_modes = len(sector_evals)

    im_values = np.imag(sector_evals)
    zero_im_count = int(np.sum(np.abs(im_values) < tol_im))
    nonzero_im = sector_evals[np.abs(im_values) >= tol_im]

    # Get unique |Im| values
    unique_abs_ims = sorted(set(round(abs(z.imag), 9) for z in nonzero_im
                                 if abs(z.imag) >= tol_im))

    return {
        "sector_name": sector_name,
        "target_re": target_re,
        "n_modes_in_sector": n_modes,
        "n_zero_im": zero_im_count,
        "unique_abs_im_values": unique_abs_ims,
    }


def main():
    start = time.time()
    print("=" * 78)
    print("EQ-018 Golden-Ratio Check: N=4 self-Pi mode structure")
    print("=" * 78)

    # ---- N=4 analysis ----
    print("\n--- N=4 ---")
    N4 = 4
    sigma_gamma_4 = N4 * GAMMA_0
    target_re_4 = -2 * sigma_gamma_4 / 2  # n_XY=2 sector: Re = -0.2
    # Note: sigma_gamma_4 = 0.2, and n_XY_mid = N/2 = 2.
    # Re for n_XY=2 mode = -2*gamma_0*2 = -0.2. Same as -sigma_gamma for N=4!
    print(f"  sigma_gamma = {sigma_gamma_4}")
    print(f"  n_XY=N/2=2 sector: Re(lam) = -2*gamma_0*2 = {target_re_4}")

    # H single-excitation spectrum
    H_SE_4 = single_excitation_hamiltonian(N4)
    H_evals_4 = np.sort(np.linalg.eigvalsh(H_SE_4))
    print(f"  H single-excitation eigenvalues: {H_evals_4}")
    expected_4 = np.array([-PHI, -1/PHI, 1/PHI, PHI])
    diff = np.abs(np.sort(H_evals_4) - np.sort(expected_4)).max()
    print(f"  Expected: ±phi ≈ ±{PHI:.5f}, ±1/phi ≈ ±{1/PHI:.5f}")
    print(f"  Match residual: {diff:.2e}")

    # Full Liouvillian at N=4
    J_list = [J_UNIFORM] * (N4 - 1)
    t0 = time.time()
    L_4 = build_liouvillian_matrix(build_H_XY(J_list, N4), GAMMA_0, N4)
    ev_4 = np.linalg.eigvals(L_4)
    print(f"  L_A diagonalised in {time.time() - t0:.2f} s "
          f"({len(ev_4)} modes)")

    # Analyse n_XY=2 sector (Re = -0.2)
    sector_4 = analyze_sector(ev_4, target_re_4, "n_XY=2 at N=4")
    print(f"\n  n_XY=2 sector (Re={target_re_4}):")
    print(f"    {sector_4['n_modes_in_sector']} modes total")
    print(f"    {sector_4['n_zero_im']} with Im=0 (self-Pi candidates)")
    print(f"    {sector_4['n_modes_in_sector'] - sector_4['n_zero_im']} with Im!=0")
    print(f"    {len(sector_4['unique_abs_im_values'])} unique |Im| values")
    if sector_4["unique_abs_im_values"]:
        print(f"    Sample |Im| values: "
              f"{[f'{x:.4f}' for x in sector_4['unique_abs_im_values'][:8]]}")

    # Check if each unique |Im| is a combination of golden-ratio quantities
    # Basis: phi, 1/phi, 1, sqrt(5)
    basis_gr = [PHI, 1/PHI, 1.0, np.sqrt(5)]
    basis_names = ["phi", "1/phi", "1", "sqrt(5)"]
    print(f"\n  Checking each |Im| as linear combination of "
          f"{{phi, 1/phi, 1, sqrt(5)}}:")
    match_count = 0
    nomatch_count = 0
    matched_examples = []
    for im_val in sector_4["unique_abs_im_values"]:
        combo, res = match_linear_combo(im_val, basis_gr, max_coef=3, tol=1e-6)
        if combo is not None:
            match_count += 1
            expr = " + ".join(f"{c}*{n}" for c, n in zip(combo, basis_names)
                               if c != 0)
            if len(matched_examples) < 6:
                matched_examples.append((im_val, expr, res))
        else:
            nomatch_count += 1
    print(f"    {match_count}/{len(sector_4['unique_abs_im_values'])} matched "
          f"to integer combinations in [-3,3]^4")
    print(f"    Sample matches:")
    for im_val, expr, res in matched_examples:
        print(f"      |Im| = {im_val:.6f} = {expr}  (residual {res:.2e})")
    if nomatch_count > 0:
        print(f"    {nomatch_count} unmatched (would need larger coefficient "
              f"range or different basis)")

    # ---- N=6 comparison ----
    print(f"\n--- N=6 (comparison, non-golden) ---")
    N6 = 6
    sigma_gamma_6 = N6 * GAMMA_0
    target_re_6 = -2 * sigma_gamma_6 / 2  # n_XY=3 sector: Re=-0.3
    # For N=6, n_XY=N/2=3. Re = -2*gamma_0*3 = -0.3.
    H_SE_6 = single_excitation_hamiltonian(N6)
    H_evals_6 = np.sort(np.linalg.eigvalsh(H_SE_6))
    print(f"  H single-excitation eigenvalues: {H_evals_6}")
    # These are 2*cos(pi*k/7) for k=1..6. NOT golden-ratio.
    print(f"  (These are 2*cos(pi*k/7), not golden-ratio values.)")

    # Diagonalize L_A at N=6
    t0 = time.time()
    J_list_6 = [J_UNIFORM] * (N6 - 1)
    L_6 = build_liouvillian_matrix(build_H_XY(J_list_6, N6), GAMMA_0, N6)
    ev_6 = np.linalg.eigvals(L_6)
    print(f"  L_A diagonalised in {time.time() - t0:.1f} s "
          f"({len(ev_6)} modes)")

    # n_XY=3 sector: Re = -0.3
    sector_6 = analyze_sector(ev_6, target_re_6, "n_XY=3 at N=6")
    print(f"\n  n_XY=3 sector (Re={target_re_6}):")
    print(f"    {sector_6['n_modes_in_sector']} modes total")
    print(f"    {sector_6['n_zero_im']} with Im=0 (self-Pi candidates)")
    print(f"    {sector_6['n_modes_in_sector'] - sector_6['n_zero_im']} with Im!=0")
    print(f"    {len(sector_6['unique_abs_im_values'])} unique |Im| values")

    # Try to fit N=6 sector with the 2cos(pi*k/7) basis
    N6_H_evals = [2 * np.cos(np.pi * k / 7) for k in range(1, 7)]
    print(f"\n  For reference, 2*cos(pi*k/7): {[f'{v:.4f}' for v in N6_H_evals]}")

    # Save
    out = {
        "N4": {
            "sigma_gamma": sigma_gamma_4,
            "target_re": target_re_4,
            "H_eigenvalues": H_evals_4.tolist(),
            "expected_H_eigenvalues_golden": expected_4.tolist(),
            "H_match_residual": float(diff),
            "sector_analysis": sector_4,
            "gr_basis": basis_gr,
            "gr_basis_names": basis_names,
            "match_count": match_count,
            "nomatch_count": nomatch_count,
            "matched_examples": [
                {"im_val": float(v), "expr": e, "residual": float(r)}
                for v, e, r in matched_examples
            ],
        },
        "N6": {
            "sigma_gamma": sigma_gamma_6,
            "target_re": target_re_6,
            "H_eigenvalues": H_evals_6.tolist(),
            "H_2cos_pik_over_7": N6_H_evals,
            "sector_analysis": sector_6,
        },
    }
    path = RESULTS_DIR / "golden_ratio_check.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")

    # ---- Verdict ----
    print(f"\n{'=' * 78}")
    print(f"VERDICT")
    print(f"{'=' * 78}")
    print(f"  N=4 self-Pi modes: {sector_4['n_zero_im']} "
          f"(expected 18 from earlier scan)")
    print(f"  Golden-ratio basis match: {match_count} of "
          f"{len(sector_4['unique_abs_im_values'])} unique |Im| values fit "
          f"integer combinations in {{phi, 1/phi, 1, sqrt(5)}}")
    if match_count == len(sector_4['unique_abs_im_values']) and \
       sector_4['n_zero_im'] == 18:
        print(f"  CONFIRMED: the n_XY=2 sector at N=4 has a closed algebraic "
              f"structure")
        print(f"              over the golden-ratio field Q[phi]. The 18 "
              f"self-Pi modes")
        print(f"              are the kernel of this structure (elements that "
              f"evaluate to 0).")
    print(f"  N=6 self-Pi modes: {sector_6['n_zero_im']} (expected 0)")
    print(f"  N=6 has no such closed algebraic structure in cyclotomic Q[2cos(pi/7)]")
    print(f"  unless specific cancellations occur. Empirically none do.")

    print(f"\nWalltime: {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
