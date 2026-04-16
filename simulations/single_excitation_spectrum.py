#!/usr/bin/env python3
"""
Single-Excitation Spectrum for Endpoint-B Chains
=================================================

Verify the closed-form formula for single-excitation dissipation rates
of the uniform open XX chain with B at endpoint:

    alpha_k / gamma_0 = (4 / (N+1)) * sin^2(k*pi / (N+1)),   k = 1..N

Three verifications:
1. Formula vs tridiagonal eigenvector computation (N=3..30)
2. Cross-check against full Liouvillian spectrum (N=3..7)
3. Niven rationality classification

Date: 2026-04-16
"""

import numpy as np
from fractions import Fraction
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results"
OUT_PATH = RESULTS_DIR / "single_excitation_spectrum.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


def formula_alphas(N):
    """Closed-form: alpha_k/gamma_0 = (4/(N+1)) * sin^2(k*pi/(N+1))."""
    return np.array([(4.0 / (N + 1)) * np.sin(k * np.pi / (N + 1))**2
                     for k in range(1, N + 1)])


def eigenvector_alphas(N):
    """
    Diagonalize N x N tridiagonal single-excitation Hamiltonian.
    H_ij = J/2 * (delta_{i,j+1} + delta_{i,j-1}) for open chain.
    alpha_k = 2 * gamma_0 * |psi_k(N-1)|^2 (Absorption Theorem).
    Returns alpha_k / gamma_0.
    """
    # tridiagonal: diagonal = 0, off-diagonal = J/2
    diag = np.zeros(N)
    offdiag = np.full(N - 1, 0.5)  # J/2 = 0.5 for J=1

    from scipy.linalg import eigh_tridiagonal
    eigenvalues, eigenvectors = eigh_tridiagonal(diag, offdiag)

    # |a_B|^2 at endpoint B = site N-1
    aB_sq = np.abs(eigenvectors[N - 1, :]) ** 2

    # alpha_k / gamma_0 = 2 * |a_B|^2 (Absorption Theorem: alpha = 2*gamma*<n_XY>)
    return 2.0 * aB_sq


def load_full_spectrum_rates(results_path, N_target):
    """
    Load distinct alpha/gamma_0 values from the full Liouvillian scan
    (structure_points_large_n.txt) for a given N.
    """
    rates = []
    in_block = False
    past_header = False
    with open(results_path) as f:
        for line in f:
            if f"N={N_target}, topology=chain, B={N_target-1}" in line:
                in_block = True
                past_header = False
                continue
            if in_block and not past_header:
                if line.startswith("---"):
                    past_header = True
                continue
            if in_block and past_header:
                if line.startswith("===") or line.startswith("###") or line.startswith("CHAIN"):
                    break
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        val = float(parts[0])
                        rates.append(val)
                    except ValueError:
                        pass
    return np.array(rates)


# =====================================================================
if __name__ == "__main__":
    log("SINGLE-EXCITATION SPECTRUM VERIFICATION")
    log("Formula: alpha_k/gamma_0 = (4/(N+1)) * sin^2(k*pi/(N+1))")
    log()

    # === Part 1: Formula vs eigenvector computation ===
    log("=" * 70)
    log("PART 1: Formula vs eigenvector (N=3..30)")
    log("=" * 70)
    log()

    max_errs = []
    for N in range(3, 31):
        alpha_formula = np.sort(formula_alphas(N))
        alpha_eigvec = np.sort(eigenvector_alphas(N))
        max_err = np.max(np.abs(alpha_formula - alpha_eigvec))
        max_errs.append(max_err)
        status = "OK" if max_err < 1e-12 else "FAIL"
        if N <= 12 or N == 30 or max_err > 1e-12:
            log(f"  N={N:2d}: {N} modes, max error = {max_err:.2e}  [{status}]")

    all_ok = all(e < 1e-12 for e in max_errs)
    log(f"\n  Formula matches eigenvector computation: "
        f"{'YES for all N=3..30' if all_ok else 'FAILS somewhere'}")
    log(f"  Worst error across all N: {max(max_errs):.2e}")

    # === Part 2: Cross-check against full Liouvillian ===
    log("\n" + "=" * 70)
    log("PART 2: Single-excitation rates inside full Liouvillian (N=3..7)")
    log("=" * 70)
    log()

    full_spectrum_path = RESULTS_DIR / "structure_points_large_n.txt"
    if not full_spectrum_path.exists():
        log("  WARNING: structure_points_large_n.txt not found, skipping")
    else:
        for N in range(3, 8):
            alpha_se = np.sort(np.unique(np.round(formula_alphas(N), 10)))
            full_rates = load_full_spectrum_rates(full_spectrum_path, N)

            if len(full_rates) == 0:
                log(f"  N={N}: no full-spectrum data found")
                continue

            # for each single-excitation rate, find closest in full spectrum
            all_found = True
            worst_match = 0
            for a in alpha_se:
                dists = np.abs(full_rates - a)
                best = np.min(dists)
                worst_match = max(worst_match, best)
                if best > 1e-3:
                    all_found = False

            status = "ALL FOUND" if all_found else "MISSING"
            log(f"  N={N}: {len(alpha_se)} SE rates, "
                f"{len(full_rates)} full rates, "
                f"worst match = {worst_match:.2e}  [{status}]")

            if N <= 5:
                log(f"        SE rates (alpha/gamma_0): "
                    f"{' '.join(f'{a:.6f}' for a in alpha_se)}")

    # === Part 3: Niven rationality classification ===
    log("\n" + "=" * 70)
    log("PART 3: Niven rationality classification")
    log("=" * 70)
    log()
    log("sin^2(k*pi/(N+1)) is rational iff (N+1) divides into {1,2,3,4,6}.")
    log("Niven's theorem: sin(r*pi) is rational only for")
    log("  sin = 0 (r integer), sin = +/-1/2 (r = 1/6, 5/6, ...), sin = +/-1.")
    log("So sin^2 rational => sin^2 in {0, 1/4, 1/2, 3/4, 1}.")
    log()

    log(f"{'N':>4} {'N+1':>4} {'all rational?':>14} {'alpha_k/gamma_0 values':>40}")
    log("-" * 70)

    for N in range(2, 21):
        alphas = formula_alphas(N)
        # check rationality: try to express as fraction with small denominator
        all_rational = True
        fracs = []
        for a in sorted(set(np.round(alphas, 12))):
            f = Fraction(a).limit_denominator(1000)
            err = abs(float(f) - a)
            if err > 1e-10:
                all_rational = False
                fracs.append(f"{a:.6f}")
            else:
                fracs.append(str(f))

        marker = "YES" if all_rational else "no"
        vals = ", ".join(fracs[:6])
        if len(fracs) > 6:
            vals += ", ..."
        log(f"{N:4d} {N+1:4d} {marker:>14} {vals:>40}")

    log()
    log("Fully rational N (all alpha_k/gamma_0 rational):")
    rational_Ns = []
    for N in range(2, 100):
        alphas = formula_alphas(N)
        all_rat = True
        for a in alphas:
            f = Fraction(a).limit_denominator(1000)
            if abs(float(f) - a) > 1e-10:
                all_rat = False
                break
        if all_rat:
            rational_Ns.append(N)
    log(f"  N in {{2..99}}: {rational_Ns}")
    log(f"  Corresponding N+1: {[n+1 for n in rational_Ns]}")
    log(f"  These are exactly N+1 in {{1, 2, 3, 4, 6}} => N in {{0, 1, 2, 3, 5}}")
    log(f"  (N=0,1 trivial; N=2,3,5 are the non-trivial Niven cases)")

    # === Part 4: Print verified values for ANALYTICAL_FORMULAS.md ===
    log("\n" + "=" * 70)
    log("PART 4: Verified values for documentation")
    log("=" * 70)
    log()

    for N in [3, 4, 5, 6]:
        alphas = sorted(set(np.round(formula_alphas(N), 10)))
        log(f"N={N} (N+1={N+1}):")
        for k in range(1, N + 1):
            a = (4.0 / (N + 1)) * np.sin(k * np.pi / (N + 1))**2
            f = Fraction(a).limit_denominator(1000)
            err = abs(float(f) - a)
            rational = f"= {f}" if err < 1e-10 else f"~ {a:.6f}"
            log(f"  k={k}: alpha/gamma_0 = (4/{N+1})*sin^2({k}*pi/{N+1}) {rational}")
        log()

    _outf.close()
    print(f"\nResults written to {OUT_PATH}")
