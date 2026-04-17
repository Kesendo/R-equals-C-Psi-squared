#!/usr/bin/env python3
"""
Bell-pair-like encoding of the palindromic partner mode (F68 H3).

For N in {3, 4, 5}, extract the partner eigenvector V_p of the chain
Liouvillian (eigenvalue -alpha_p with alpha_p = 2 gamma_0 - alpha_b).
SVD V_p to get V_p = sigma_0 |u><v| (rank-1 for N >= 4 per F68 H2).
Build an R-C Bell-pair-like pure state (|0>_R|u> + |1>_R|v>)/sqrt(2),
propagate on R + chain (R isolated, dephasing at chain site N-1), and fit
the off-diagonal |0><1|_R block decay rate.

In the same script, propagate the F67 Variant B bonding-mode Bell pair
(|0>_R|vac> + |1>_R|psi_1>)/sqrt(2) and check the dynamical palindrome

    alpha_fit(bonding) + alpha_fit(partner) =? 2 gamma_0.

N = 3 is a negative control (rank-2 partner): the rank-1 approximations
do not give a clean single-exponential decay.

Because R is fully decoupled from the chain (no Hamiltonian, no jump),
each R-block of rho_ext evolves independently under the chain Liouvillian.
We exploit this to avoid building the 4^(N+1)-dim extended superoperator:
blk_00 = (1/2)|u><u|, blk_01 = (1/2)|u><v|, blk_11 = (1/2)|v><v| each
evolve under L_chain. Reassembly gives the full rho_ext(t) for negativity.

Date: 2026-04-17
"""

import numpy as np
from pathlib import Path
import sys
import time

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)
OUT_PATH = RESULTS_DIR / "bell_pair_partner_mode.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ---- Paulis ----
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron_chain(*factors)


def xy_hamiltonian(N, J=1.0):
    """XY Hamiltonian on N chain qubits (site 0 most significant)."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        H += J * 0.5 * (site_op(X, i, N) @ site_op(X, i + 1, N) +
                        site_op(Y, i, N) @ site_op(Y, i + 1, N))
    return H


def liouvillian_superop(H, jump_ops):
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
    """F65 bonding-mode rate (perturbative in gamma_0)."""
    return (4.0 * gamma_0 / (N + 1)) * np.sin(np.pi / (N + 1))**2


def single_excitation_mode(N, k=1):
    """|psi_k>_chain = sqrt(2/(N+1)) sum_i sin(pi k (i+1)/(N+1)) |i excited>.
    Site 0 most significant; chain qubit i excited = index 2^(N-1-i)."""
    psi = np.zeros(2**N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
        psi[2**(N - 1 - i)] = amp
    return psi


def vacuum_state(N):
    vac = np.zeros(2**N, dtype=complex)
    vac[0] = 1.0
    return vac


def build_chain_liouvillian(N, gamma_0, J):
    H = xy_hamiltonian(N, J)
    L_jump = np.sqrt(gamma_0) * site_op(Z, N - 1, N)
    return liouvillian_superop(H, [L_jump])


def find_bonding_partner(L_chain, N, gamma_0, alpha_1_formula):
    """Find bonding and partner eigenvectors of L_chain.

    Bonding: eigenvalue group closest to -alpha_1 (the F65 perturbative rate).
    Partner: eigenvalue group closest to -(2*gamma_0 - alpha_b_actual).

    To match the F67 convention (and get a rank-1 basis for degenerate
    eigenspaces), we pick the LOWEST-INDEX eigenvector within the
    degenerate group (partner_idx[0]), not argmin of the distance array.
    For multiplicity-16/20 subspaces, numerical noise at 1e-15 level can
    otherwise reroute argmin to a non-rank-1 basis element.

    Returns (alpha_b, alpha_p, mult_b, mult_p, V_b, V_p)."""
    d = 2**N
    eigenvalues, eigenvectors = np.linalg.eig(L_chain)

    # Bonding
    dists_b = np.abs(eigenvalues.real + alpha_1_formula)
    best_b = int(np.argmin(dists_b))
    alpha_b_center = -float(eigenvalues[best_b].real)
    tol_b = max(1e-8, abs(alpha_b_center) * 1e-4)
    bonding_idx = np.where(np.abs(eigenvalues.real - eigenvalues[best_b].real) < tol_b)[0]
    mult_b = int(len(bonding_idx))
    ev_b = int(bonding_idx[0])
    alpha_b = -float(eigenvalues[ev_b].real)
    V_b = eigenvectors[:, ev_b].reshape(d, d, order='F')

    # Partner (use actual alpha_b, not the F65 formula)
    alpha_p_target = 2.0 * gamma_0 - alpha_b
    dists_p = np.abs(eigenvalues.real + alpha_p_target)
    best_p = int(np.argmin(dists_p))
    alpha_p_center = -float(eigenvalues[best_p].real)
    tol_p = max(1e-8, abs(alpha_p_center) * 1e-4)
    partner_idx = np.where(np.abs(eigenvalues.real - eigenvalues[best_p].real) < tol_p)[0]
    mult_p = int(len(partner_idx))
    ev_p = int(partner_idx[0])
    alpha_p = -float(eigenvalues[ev_p].real)
    V_p = eigenvectors[:, ev_p].reshape(d, d, order='F')

    return alpha_b, alpha_p, mult_b, mult_p, V_b, V_p


def svd_rank1_factors(V):
    """V = sum_i sigma_i |u_i><v_i|. Returns (sigmas, u_0, v_0) with u_0, v_0
    unit vectors (columns) such that V approx sigma_0 |u_0><v_0|."""
    U, sigmas, Vh = np.linalg.svd(V)
    u = U[:, 0]
    v = Vh[0, :].conj()
    return sigmas, u, v


def propagate_mat_via_eig(mat_0, lam, Vmat, Vmat_inv, times, N):
    """Given eigendecomposition (lam, Vmat, Vmat_inv) of L_chain, propagate
    a 2^N x 2^N chain matrix mat_0. Returns an array of shape (nt, d, d)."""
    d = 2**N
    vec0 = mat_0.flatten(order='F')
    coefs = Vmat_inv @ vec0
    out = np.empty((len(times), d, d), dtype=complex)
    for i, t in enumerate(times):
        vec_t = Vmat @ (np.exp(lam * t) * coefs)
        out[i] = vec_t.reshape(d, d, order='F')
    return out


def assemble_rho_ext(blk_00, blk_01, blk_11, N):
    """Assemble 2^(N+1) x 2^(N+1) density matrix from its R-blocks.
    R at extended site 0 (MSB): rho[0..d, 0..d] = blk_00, etc."""
    d = 2**N
    D = 2 * d
    rho = np.empty((D, D), dtype=complex)
    rho[:d, :d] = blk_00
    rho[:d, d:] = blk_01
    rho[d:, :d] = blk_01.conj().T
    rho[d:, d:] = blk_11
    return rho


def partial_transpose_R(rho_ext, N):
    """Partial transpose w.r.t. R (R at extended site 0, MSB)."""
    dimR = 2
    dimC = 2**N
    rho_reshape = rho_ext.reshape(dimR, dimC, dimR, dimC)
    rho_PT = rho_reshape.transpose(2, 1, 0, 3)
    return rho_PT.reshape(dimR * dimC, dimR * dimC)


def rc_negativity(rho_ext, N):
    """Negativity N(R:C) = sum |lambda_-| over negative PT eigenvalues."""
    rho_PT = partial_transpose_R(rho_ext, N)
    rho_PT_H = (rho_PT + rho_PT.conj().T) / 2.0
    eigs = np.linalg.eigvalsh(rho_PT_H)
    return float(np.sum(np.abs(eigs[eigs < 0])))


def propagate_bellpair(u, v, L_chain_eig, times, N):
    """Propagate rho_0 = (|0>_R|u> + |1>_R|v>)(h.c.)/2 under the extended
    Liouvillian (R decoupled + chain with dephasing). Returns arrays of
    (off-diag R-block Frobenius norm, R:C negativity, Hermiticity residual)."""
    lam, Vmat, Vmat_inv = L_chain_eig
    blk_00_0 = 0.5 * np.outer(u, u.conj())
    blk_01_0 = 0.5 * np.outer(u, v.conj())
    blk_11_0 = 0.5 * np.outer(v, v.conj())

    arr_00 = propagate_mat_via_eig(blk_00_0, lam, Vmat, Vmat_inv, times, N)
    arr_01 = propagate_mat_via_eig(blk_01_0, lam, Vmat, Vmat_inv, times, N)
    arr_11 = propagate_mat_via_eig(blk_11_0, lam, Vmat, Vmat_inv, times, N)

    nt = len(times)
    norms = np.empty(nt)
    negs = np.empty(nt)
    herm_errs = np.empty(nt)
    for i in range(nt):
        norms[i] = float(np.linalg.norm(arr_01[i], ord='fro'))
        rho_ext = assemble_rho_ext(arr_00[i], arr_01[i], arr_11[i], N)
        herm_errs[i] = float(np.linalg.norm(rho_ext - rho_ext.conj().T, ord='fro'))
        negs[i] = rc_negativity(rho_ext, N)
    return norms, negs, herm_errs


def fit_exponential(times, norms, floor=1e-12):
    """Fit log(max(norm, floor)) = intercept - alpha * t linearly.
    Returns (alpha_fit, residual_rms, log_norms)."""
    norms_clipped = np.maximum(norms, floor)
    log_norms = np.log(norms_clipped)
    slope, intercept = np.polyfit(times, log_norms, 1)
    alpha_fit = float(-slope)
    predicted = intercept + slope * times
    residuals = log_norms - predicted
    residual_rms = float(np.sqrt(np.mean(residuals**2)))
    return alpha_fit, residual_rms, log_norms


# =========================================================================
if __name__ == "__main__":
    log("BELL-PAIR-LIKE ENCODING OF THE PALINDROMIC PARTNER MODE (F68 H3)")
    log("=" * 72)
    log("Setup: R + N-site XY chain, Z-dephasing gamma_0 at chain site N-1.")
    log("Metric: Frobenius norm of the |0><1|_R block of rho_ext(t),")
    log("        fit log||block||(t) = a - alpha_fit * t.")
    log("Target: alpha_fit(partner) vs spectral alpha_p = 2*gamma_0 - alpha_b.")
    log()

    gamma_0 = 0.05
    J = 1.0
    n_points = 80
    n_efolds = 5.0

    log(f"Parameters: gamma_0 = {gamma_0}, J = {J}")
    log(f"            2*gamma_0 = {2*gamma_0:.10f}")
    log(f"            {n_points} time points over {n_efolds} e-folds")
    log()

    summary = {}

    for N in [3, 4, 5]:
        log("=" * 72)
        log(f"N = {N}   (L_chain dim = {4**N}, rho_ext dim = {2**(N+1)})")
        log("=" * 72)

        t0 = time.time()
        L_chain = build_chain_liouvillian(N, gamma_0, J)
        alpha_1 = formula_alpha_1(N, gamma_0)
        alpha_b, alpha_p, mult_b, mult_p, V_b, V_p = find_bonding_partner(
            L_chain, N, gamma_0, alpha_1)
        log(f"  spectral extraction + eig: {time.time() - t0:.2f} s")

        log(f"  alpha_1 (F65 perturbative)  = {alpha_1:.10f}")
        log(f"  alpha_b (full L_chain)      = {alpha_b:.10f}   (mult = {mult_b})")
        log(f"  alpha_p (full L_chain)      = {alpha_p:.10f}   (mult = {mult_p})")
        log(f"  alpha_b + alpha_p           = {alpha_b + alpha_p:.12f}")
        log(f"  |sum - 2*gamma_0|           = {abs(alpha_b + alpha_p - 2*gamma_0):.2e}")

        t0 = time.time()
        lam_c, Vmat = np.linalg.eig(L_chain)
        Vmat_inv = np.linalg.inv(Vmat)
        L_chain_eig = (lam_c, Vmat, Vmat_inv)
        log(f"  L_chain diagonalization: {time.time() - t0:.2f} s")

        sigmas_p, u_p, v_p = svd_rank1_factors(V_p)
        sigmas_b, u_b_svd, v_b_svd = svd_rank1_factors(V_b)
        log(f"\n  Part 1: SVD spectrum of V_p")
        for i in range(min(3, len(sigmas_p))):
            log(f"    sigma_{i}(V_p) = {sigmas_p[i]:.6e}")
        rank_ratio_p = float(sigmas_p[1] / sigmas_p[0]) if sigmas_p[0] > 1e-15 else np.inf
        log(f"    sigma_1/sigma_0 = {rank_ratio_p:.2e}")
        # Task spec threshold is 1e-10 (machine precision). In practice LAPACK
        # zgeev leaves residuals up to ~1e-7 on highly degenerate subspaces
        # even when the true rank is 1. We treat ratio < 1e-6 as "numerically
        # rank-1" and report the exact ratio against the strict 1e-10 mark.
        is_rank1_strict = rank_ratio_p < 1e-10
        is_rank1_numerical = rank_ratio_p < 1e-6
        log(f"    Partner rank-1 strict (< 1e-10)?    {'YES' if is_rank1_strict else 'NO'}")
        log(f"    Partner rank-1 numerical (< 1e-6)?  {'YES' if is_rank1_numerical else 'NO'}")
        log(f"  SVD of V_b (bonding):")
        log(f"    sigma_0 = {sigmas_b[0]:.6e}, sigma_1 = {sigmas_b[1]:.6e}, "
            f"ratio = {sigmas_b[1]/sigmas_b[0]:.2e}")

        # Time grids (5 e-folds of the respective rate)
        times_p = np.linspace(0.0, n_efolds / alpha_p, n_points)
        times_b = np.linspace(0.0, n_efolds / alpha_b, n_points)

        # ---- Part 2: partner R-C encoding ---------------------------------
        alpha_fit_p = rel_err_p = resid_p = None
        min_neg_p = max_neg_p = None
        # Propagate partner encoding when numerically rank-1 (N >= 4).
        # For N = 3 the partner is genuinely rank-2 (sigma_1/sigma_0 ~ 1),
        # handled separately in Part 5 as a negative control.
        if is_rank1_numerical:
            log(f"\n  Part 2: Partner R-C encoding -> propagation")
            log(f"    rho_0 = (|0>_R|u_0> + |1>_R|v_0>)(h.c.) / 2")
            t0 = time.time()
            norms_p, negs_p, herm_p = propagate_bellpair(
                u_p, v_p, L_chain_eig, times_p, N)
            log(f"    propagation time: {time.time() - t0:.2f} s")
            alpha_fit_p, resid_p, _ = fit_exponential(times_p, norms_p)
            rel_err_p = abs(alpha_fit_p - alpha_p) / alpha_p
            min_neg_p = float(np.min(negs_p))
            max_neg_p = float(np.max(negs_p))
            log(f"    alpha_fit(partner) = {alpha_fit_p:.10f}")
            log(f"    alpha_p (spectral) = {alpha_p:.10f}")
            log(f"    relative error     = {rel_err_p:.2e}")
            log(f"    log-fit resid RMS  = {resid_p:.2e}")
            log(f"    ||rho_ext - rho_ext^dag||_F range = "
                f"[{float(herm_p.min()):.2e}, {float(herm_p.max()):.2e}]")
            log(f"    negativity N(R:C) at t=0       = {negs_p[0]:.6f}")
            log(f"    negativity N(R:C) at t_max     = {negs_p[-1]:.6f}")
            log(f"    negativity N(R:C) min throughout = {min_neg_p:.6e}")
            log(f"    negativity > 0 over whole window? "
                f"{'YES' if min_neg_p > 1e-8 else 'NO'}")
        else:
            log(f"\n  Part 2: Partner R-C encoding SKIPPED (partner not rank-1)")

        # ---- Part 3: bonding R-C encoding (F67 Variant B) -----------------
        log(f"\n  Part 3: Bonding-mode R-C encoding (F67 Variant B)")
        log(f"    rho_0 = (|0>_R|vac> + |1>_R|psi_1>)(h.c.) / 2")
        u_b_enc = vacuum_state(N)
        v_b_enc = single_excitation_mode(N, k=1)
        t0 = time.time()
        norms_b, negs_b, herm_b = propagate_bellpair(
            u_b_enc, v_b_enc, L_chain_eig, times_b, N)
        log(f"    propagation time: {time.time() - t0:.2f} s")
        alpha_fit_b, resid_b, _ = fit_exponential(times_b, norms_b)
        rel_err_b = abs(alpha_fit_b - alpha_b) / alpha_b
        min_neg_b = float(np.min(negs_b))
        max_neg_b = float(np.max(negs_b))
        log(f"    alpha_fit(bonding) = {alpha_fit_b:.10f}")
        log(f"    alpha_b (spectral) = {alpha_b:.10f}")
        log(f"    relative error     = {rel_err_b:.2e}")
        log(f"    log-fit resid RMS  = {resid_b:.2e}")
        log(f"    negativity N(R:C) at t=0        = {negs_b[0]:.6f}")
        log(f"    negativity N(R:C) at t_max      = {negs_b[-1]:.6f}")
        log(f"    negativity N(R:C) min           = {min_neg_b:.6e}")

        # ---- Part 4: dynamical palindromic identity -----------------------
        if alpha_fit_p is not None:
            log(f"\n  Part 4: Dynamical palindromic identity")
            dyn_sum = alpha_fit_b + alpha_fit_p
            dyn_err = abs(dyn_sum - 2.0 * gamma_0)
            rel_dyn = dyn_err / (2.0 * gamma_0)
            log(f"    alpha_fit(bonding) + alpha_fit(partner) = {dyn_sum:.10f}")
            log(f"    2*gamma_0                               = {2*gamma_0:.10f}")
            log(f"    absolute error = {dyn_err:.2e}")
            log(f"    relative error = {rel_dyn:.2e}")

            summary[N] = dict(
                alpha_1=alpha_1, alpha_b=alpha_b, alpha_p=alpha_p,
                alpha_fit_b=alpha_fit_b, alpha_fit_p=alpha_fit_p,
                rel_err_b=rel_err_b, rel_err_p=rel_err_p,
                dyn_sum=dyn_sum, rel_dyn=rel_dyn,
                resid_b=resid_b, resid_p=resid_p,
                rank_ratio=rank_ratio_p,
                min_neg_p=min_neg_p, min_neg_b=min_neg_b,
                is_rank1=True,
            )
        else:
            summary[N] = dict(
                alpha_1=alpha_1, alpha_b=alpha_b, alpha_p=alpha_p,
                alpha_fit_b=alpha_fit_b, alpha_fit_p=None,
                rel_err_b=rel_err_b, rel_err_p=None,
                dyn_sum=None, rel_dyn=None,
                resid_b=resid_b, resid_p=None,
                rank_ratio=rank_ratio_p,
                min_neg_p=None, min_neg_b=min_neg_b,
                is_rank1=False,
            )

        # ---- Part 5: N=3 rank-2 negative control --------------------------
        if N == 3:
            log(f"\n  Part 5: N=3 rank-2 negative control")
            U_p_full, sig_p_full, Vh_p_full = np.linalg.svd(V_p)
            log(f"    V_p singular spectrum (top 4):")
            for i in range(min(4, len(sig_p_full))):
                log(f"      sigma_{i} = {sig_p_full[i]:.6e}")
            log(f"    sigma_1/sigma_0 = {sig_p_full[1]/sig_p_full[0]:.4f}  "
                f"-> genuinely rank-2, not close to rank-1")
            log()
            log(f"    Try each rank-1 approximation |u_i><v_i| as R-C off-diag:")
            n3_negcontrol = []
            for i in [0, 1]:
                u_i = U_p_full[:, i]
                v_i = Vh_p_full[i, :].conj()
                norms_i, _, _ = propagate_bellpair(
                    u_i, v_i, L_chain_eig, times_p, N)
                alpha_fit_i, resid_i, _ = fit_exponential(times_p, norms_i)
                rel_err_i = abs(alpha_fit_i - alpha_p) / alpha_p
                n3_negcontrol.append(dict(
                    i=i, alpha_fit=alpha_fit_i, resid=resid_i,
                    rel_err=rel_err_i,
                    norm_series=[float(norms_i[0]),
                                 float(norms_i[n_points // 4]),
                                 float(norms_i[n_points // 2]),
                                 float(norms_i[3 * n_points // 4]),
                                 float(norms_i[-1])],
                ))
                log(f"    - Using |u_{i}><v_{i}|:")
                log(f"        alpha_fit          = {alpha_fit_i:.10f}")
                log(f"        alpha_p (target)   = {alpha_p:.10f}")
                log(f"        relative error     = {rel_err_i:.2e}")
                log(f"        log-fit resid RMS  = {resid_i:.2e}   "
                    f"({'exponential' if resid_i < 0.02 else 'non-exponential: beats / multi-mode'})")
                log(f"        norm samples: t=0  -> {norms_i[0]:.4e}")
                log(f"                      t/4  -> {norms_i[n_points//4]:.4e}")
                log(f"                      t/2  -> {norms_i[n_points//2]:.4e}")
                log(f"                      3t/4 -> {norms_i[3*n_points//4]:.4e}")
                log(f"                      t_f  -> {norms_i[-1]:.4e}")
            summary[N]['n3_negcontrol'] = n3_negcontrol
            summary[N]['sig_spectrum'] = [float(x) for x in sig_p_full[:4]]

        log()

    # -----------------------------------------------------------------------
    log("=" * 72)
    log("FINAL SUMMARY")
    log("=" * 72)
    log()
    log("Spectral + dynamical rates")
    log(f"  {'N':>3} {'alpha_b':>12} {'alpha_p':>12} "
        f"{'fit(bond)':>12} {'fit(part)':>12} {'dyn sum':>12} {'rel err':>10}")
    log(f"  {'-'*3} {'-'*12} {'-'*12} {'-'*12} {'-'*12} {'-'*12} {'-'*10}")
    for N in [3, 4, 5]:
        s = summary[N]
        fit_p_str = f"{s['alpha_fit_p']:12.6f}" if s['alpha_fit_p'] is not None else f"{'rank-2':>12}"
        sum_str = f"{s['dyn_sum']:12.6f}" if s['dyn_sum'] is not None else f"{'---':>12}"
        rel_str = f"{s['rel_dyn']:10.2e}" if s['rel_dyn'] is not None else f"{'---':>10}"
        log(f"  {N:3d} {s['alpha_b']:12.6f} {s['alpha_p']:12.6f} "
            f"{s['alpha_fit_b']:12.6f} {fit_p_str} {sum_str} {rel_str}")

    log()
    log("Acceptance criteria (task spec):")
    for N in [4, 5]:
        s = summary[N]
        r1_strict = s['rank_ratio'] < 1e-10
        r1_num = s['rank_ratio'] < 1e-6
        if s['rel_err_p'] is not None:
            r2 = s['rel_err_p'] < 1e-3
            r4 = s['rel_dyn'] < 1e-3
        else:
            r2 = r4 = False
        log(f"  N={N}:")
        log(f"    rank-1 strict (< 1e-10):            "
            f"{'PASS' if r1_strict else 'FAIL'}   (ratio = {s['rank_ratio']:.2e})")
        log(f"    rank-1 numerical (< 1e-6):          "
            f"{'PASS' if r1_num else 'FAIL'}   (ratio = {s['rank_ratio']:.2e})")
        if s['rel_err_p'] is not None:
            log(f"    alpha_fit(part)/alpha_p < 1e-3 rel: "
                f"{'PASS' if r2 else 'FAIL'}   (rel err = {s['rel_err_p']:.2e})")
            log(f"    palindrome sum within 1e-3 rel:   "
                f"{'PASS' if r4 else 'FAIL'}   (rel err = {s['rel_dyn']:.2e})")
        else:
            log(f"    (Part 2/4 skipped: partner not numerically rank-1)")

    log()
    log("Bonus criteria:")
    for N in [4, 5]:
        s = summary[N]
        if s['resid_p'] is not None:
            good_resid = (s['resid_p'] < 0.01) and (s['resid_b'] < 0.01)
            pos_neg = s['min_neg_p'] > 1e-8
            log(f"  N={N}:")
            log(f"    log-fit residuals < 0.01 for both:  "
                f"{'YES' if good_resid else 'NO'}  "
                f"(bond resid {s['resid_b']:.2e}, part resid {s['resid_p']:.2e})")
            log(f"    partner N(R:C) > 0 throughout:      "
                f"{'YES' if pos_neg else 'NO'}  (min N = {s['min_neg_p']:.2e})")
        else:
            log(f"  N={N}: (skipped, partner not rank-1)")

    log()
    log("N=3 (rank-2 degeneracy, negative control):")
    s = summary[3]
    log(f"    V_p SVD top 4:           {s['sig_spectrum']}")
    log(f"    bonding fit: alpha_fit(b) = {s['alpha_fit_b']:.6f} "
        f"vs alpha_b = {s['alpha_b']:.6f}  (rel err = {s['rel_err_b']:.2e}, "
        f"resid {s['resid_b']:.2e})")
    log(f"    rank-1 approximations of V_p:")
    for entry in s['n3_negcontrol']:
        log(f"      i={entry['i']}: alpha_fit = {entry['alpha_fit']:.6f}  "
            f"(rel err {entry['rel_err']:.2e}, resid {entry['resid']:.2e})")

    log()
    log("=" * 72)
    _outf.close()
    print(f"\nResults written to {OUT_PATH}")
