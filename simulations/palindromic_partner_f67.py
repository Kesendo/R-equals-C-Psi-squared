#!/usr/bin/env python3
"""
Palindromic partner of the F67 bonding-mode encoding.

The bonding mode k=1 decays at alpha_1. By F1 (palindromic pairing),
there must exist a partner mode at 2*gamma_0 - alpha_1. This script
finds it, analyzes its Pauli structure, and tests whether it admits
an operational Bell-pair encoding.

Date: 2026-04-16
"""

import numpy as np
from itertools import product as iprod
from scipy.linalg import expm
from scipy.optimize import curve_fit
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results"
OUT_PATH = RESULTS_DIR / "palindromic_partner_f67.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, X, Y, Z]
PAULI_NAMES = ['I', 'X', 'Y', 'Z']


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron_chain(*factors)


def build_xy_chain(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        H += J * 0.5 * (site_op(X, i, N) @ site_op(X, i + 1, N)
                        + site_op(Y, i, N) @ site_op(Y, i + 1, N))
    return H


def liouvillian(H, jump_ops):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Idd, LdL)
              - 0.5 * np.kron(LdL.T, Idd))
    return L


def formula_alpha_1(N, gamma_0):
    return (4.0 * gamma_0 / (N + 1)) * np.sin(np.pi / (N + 1))**2


def pauli_label(idx):
    return ''.join(PAULI_NAMES[i] for i in idx)


def xy_weight_at_site(idx, site):
    return 1 if idx[site] in (1, 2) else 0


def total_xy_weight(idx):
    return sum(1 for i in idx if i in (1, 2))


def pauli_decompose(rho_mat, N):
    """Decompose a d x d matrix in the Pauli basis. Returns dict idx -> coeff."""
    d = 2 ** N
    all_idx = list(iprod(range(4), repeat=N))
    coeffs = {}
    for idx in all_idx:
        P = PAULIS[idx[0]]
        for k in range(1, N):
            P = np.kron(P, PAULIS[idx[k]])
        c = np.trace(P @ rho_mat) / d
        if abs(c) > 1e-12:
            coeffs[idx] = c
    return coeffs


def xy_weight_distribution(coeffs, N):
    """Sum |c|^2 by total XY-weight."""
    dist = {}
    for idx, c in coeffs.items():
        w = total_xy_weight(idx)
        dist[w] = dist.get(w, 0) + abs(c)**2
    return dist


# =====================================================================
if __name__ == "__main__":
    log("PALINDROMIC PARTNER OF THE F67 BONDING-MODE ENCODING")
    log("=" * 70)

    gamma_0 = 0.05
    J = 1.0

    for N in [3, 4, 5]:
        log(f"\n{'=' * 70}")
        log(f"N = {N}, dim(L) = {4**N}")
        log(f"{'=' * 70}")

        d = 2 ** N
        alpha_1 = formula_alpha_1(N, gamma_0)
        alpha_partner = 2 * gamma_0 - alpha_1

        log(f"gamma_0 = {gamma_0}, J = {J}")
        log(f"alpha_1 (F65) = {alpha_1:.10f}")
        log(f"predicted partner = 2*gamma_0 - alpha_1 = {alpha_partner:.10f}")

        H = build_xy_chain(N, J)
        Lop = np.sqrt(gamma_0) * site_op(Z, N - 1, N)
        L = liouvillian(H, [Lop])

        eigenvalues, eigenvectors = np.linalg.eig(L)

        # === Part 1: Spectral check (H1) ===
        log(f"\n--- Part 1: Spectral pairing (H1) ---")

        # find eigenvalues closest to -alpha_1
        # use wider tolerance: gamma_0=0.05 is not perturbative, eigenvalues
        # shift from the gamma_0->0 formula by O(gamma_0^2/J)
        target_b = -alpha_1
        dists_b = np.abs(eigenvalues.real - target_b)
        # first find THE closest eigenvalue, then use its neighborhood
        best_b = np.argmin(dists_b)
        actual_alpha_b = -eigenvalues[best_b].real
        log(f"Closest to -alpha_1: Re = {eigenvalues[best_b].real:.10f}, "
            f"actual alpha = {actual_alpha_b:.10f}, "
            f"shift from formula = {actual_alpha_b - alpha_1:.2e}")

        # group all eigenvalues at the same real part (within 1e-8)
        tol = max(1e-8, abs(actual_alpha_b) * 1e-4)
        bonding_idx = np.where(np.abs(eigenvalues.real - eigenvalues[best_b].real) < tol)[0]
        mult_b = len(bonding_idx)

        # partner: use actual bonding alpha to compute expected partner
        actual_partner = 2 * gamma_0 - actual_alpha_b
        target_p = -actual_partner
        dists_p = np.abs(eigenvalues.real - target_p)
        best_p = np.argmin(dists_p)
        actual_alpha_p = -eigenvalues[best_p].real
        log(f"Closest to partner: Re = {eigenvalues[best_p].real:.10f}, "
            f"actual alpha = {actual_alpha_p:.10f}, "
            f"shift from predicted = {actual_alpha_p - actual_partner:.2e}")

        partner_idx = np.where(np.abs(eigenvalues.real - eigenvalues[best_p].real) < tol)[0]
        mult_p = len(partner_idx)

        log(f"Bonding mode at Re = {target_b:.10f}: multiplicity {mult_b}")
        log(f"Partner mode at Re = {target_p:.10f}: multiplicity {mult_p}")

        log(f"Bonding multiplicity: {mult_b}")
        log(f"Partner multiplicity: {mult_p}")

        if mult_b > 0 and mult_p > 0:
            pair_sum = actual_alpha_b + actual_alpha_p
            err = abs(pair_sum - 2 * gamma_0)
            log(f"Pairing check: alpha_b + alpha_p = {pair_sum:.12f}")
            log(f"  expected 2*gamma_0 = {2*gamma_0:.12f}")
            log(f"  error = {err:.2e}")
            log(f"  H1: {'PASS' if err < 1e-8 else 'FAIL'}")

        # === Part 2: Structural analysis (H2) ===
        log(f"\n--- Part 2: Structural analysis (H2) ---")

        for label, idx_set in [("bonding", bonding_idx), ("partner", partner_idx)]:
            if len(idx_set) == 0:
                log(f"  {label}: no eigenvector found")
                continue

            ev_idx = idx_set[0]
            evec = eigenvectors[:, ev_idx]
            rho_mat = evec.reshape(d, d, order='F')

            # <n_XY>_B via Absorption Theorem
            alpha_val = -eigenvalues[ev_idx].real
            nxy_B = alpha_val / (2 * gamma_0)
            log(f"\n  [{label}] eigenvalue Re = {eigenvalues[ev_idx].real:.10f}")
            log(f"    alpha = {alpha_val:.10f}")
            log(f"    <n_XY>_B = alpha / (2*gamma_0) = {nxy_B:.6f}")

            # Pauli decomposition
            coeffs = pauli_decompose(rho_mat, N)
            xy_dist = xy_weight_distribution(coeffs, N)

            log(f"    XY-weight distribution:")
            total_norm = sum(xy_dist.values())
            for w in sorted(xy_dist.keys()):
                frac = xy_dist[w] / total_norm * 100
                log(f"      w={w}: {frac:6.2f}%")

            # top Pauli strings
            sorted_coeffs = sorted(coeffs.items(), key=lambda x: -abs(x[1]))
            log(f"    Top 8 Pauli strings:")
            for idx, c in sorted_coeffs[:8]:
                label_str = pauli_label(idx)
                w = total_xy_weight(idx)
                wB = xy_weight_at_site(idx, N - 1)
                log(f"      {label_str}: |c|={abs(c):.4f}, w={w}, n_XY(B)={wB}")

        # === Part 3: Operational check (H3) ===
        log(f"\n--- Part 3: Rank-1 check (H3) ---")

        if len(partner_idx) > 0:
            ev_idx = partner_idx[0]
            evec = eigenvectors[:, ev_idx]
            V_p = evec.reshape(d, d, order='F')

            # SVD of V_p
            U, sigmas, Vt = np.linalg.svd(V_p)
            log(f"  Partner mode SVD singular values (top 5):")
            for i in range(min(5, len(sigmas))):
                log(f"    sigma_{i} = {sigmas[i]:.6f}")

            ratio_12 = sigmas[1] / sigmas[0] if sigmas[0] > 1e-15 else np.inf
            log(f"  sigma_1/sigma_0 = {ratio_12:.6f}")
            is_rank1 = ratio_12 < 0.01
            log(f"  Rank-1? {'YES' if is_rank1 else 'NO'}")

            # also check bonding mode
            b_evec = eigenvectors[:, bonding_idx[0]]
            V_b = b_evec.reshape(d, d, order='F')
            U_b, sig_b, _ = np.linalg.svd(V_b)
            ratio_b = sig_b[1] / sig_b[0] if sig_b[0] > 1e-15 else np.inf
            log(f"  Bonding mode sigma_1/sigma_0 = {ratio_b:.6f} "
                f"({'rank-1' if ratio_b < 0.01 else 'not rank-1'})")

        log()

    # === Summary ===
    log("\n" + "=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()
    log("H1 (spectral pairing): check results above per N.")
    log("H2 (structural): compare XY-weight distributions.")
    log("H3 (operational): check rank-1 assessment.")

    _outf.close()
    print(f"\nResults written to {OUT_PATH}")
