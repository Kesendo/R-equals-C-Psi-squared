"""
Ask the Lens: extracting the optimal state from the slow mode eigenvector.
==========================================================================
Steps 1-7 as specified in TASK_SLOW_MODE_LENS.md.

Convention: ibm_april_predictions Liouvillian (np.kron(H, Id) for commutator).
Cached eigendecomposition from optimal_state_n5_sacrifice.PROFILES.

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 9, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import time as _time
import numpy as np
from scipy import linalg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import (
    heisenberg_H, build_liouvillian, pauli_basis, n_xy_string,
    partial_trace_to_pair, wootters_concurrence,
)
from optimal_state_n5_sacrifice import (
    PROFILES, N, gamma_sacrifice, gamma_uniform, Sg_sac,
    max_adjacent_concurrence,
)

D = 2**N  # 32
D2 = D**2  # 1024

# Single-excitation basis indices in the d=32 computational basis
SE_IDX = [1 << (N - 1 - k) for k in range(N)]  # [16, 8, 4, 2, 1]


def safe_wootters(rho):
    c = wootters_concurrence(rho)
    return 0.0 if not np.isfinite(c) else max(0.0, float(c))


def concurrence_trajectory(psi0, profile, times):
    rho0 = np.outer(psi0, psi0.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    traj = np.zeros(len(times))
    for i, t in enumerate(times):
        rho_vec = profile['R'] @ (c0 * np.exp(profile['eigvals'] * t))
        rho_t = rho_vec.reshape(D, D)
        rho_t = (rho_t + rho_t.conj().T) / 2
        concs = [safe_wootters(partial_trace_to_pair(rho_t, N, k, k+1))
                 for k in range(N - 1)]
        traj[i] = np.mean(concs)
    return traj


def auc_windows(t, traj, windows=(2.0, 10.0, 30.0)):
    return {f"auc_{int(T)}": float(np.trapezoid(traj[t <= T], t[t <= T]))
            for T in windows}


# ===================================================================
if __name__ == "__main__":
    t_start = _time.time()
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "results", "slow_mode_lens")
    os.makedirs(out_dir, exist_ok=True)

    prof = PROFILES['sacrifice']
    eigvals = prof['eigvals']
    R = prof['R']
    R_inv = prof['R_inv']
    rates = prof['rates']

    print("=" * 80)
    print("ASK THE LENS: Slow Mode Eigenvector Analysis")
    print("=" * 80)

    # ==================================================================
    # Step 1: Locate slow modes (both the absolute slowest AND the W5-relevant one)
    # ==================================================================
    print("\n--- Step 1: Locate slow modes ---")
    nonstat = rates > 1e-10
    nonstat_idx = np.where(nonstat)[0]

    # Absolute slowest non-stationary mode
    slowest_k = nonstat_idx[np.argmin(rates[nonstat_idx])]
    print(f"  Absolute slowest: index {slowest_k}, rate {rates[slowest_k]:.6f}, "
          f"ev = {eigvals[slowest_k].real:+.6f} {eigvals[slowest_k].imag:+.6f}j")

    # The W5-relevant mode at rate ~0.318 (from Phase 1 diagnosis)
    target_rate = 0.318
    diffs = np.abs(rates[nonstat_idx] - target_rate)
    w5_mode_k = nonstat_idx[np.argmin(diffs)]
    ev_w5mode = eigvals[w5_mode_k]
    print(f"  W5-relevant mode: index {w5_mode_k}, rate {rates[w5_mode_k]:.6f}, "
          f"ev = {ev_w5mode.real:+.6f} {ev_w5mode.imag:+.6f}j")

    # Use the W5-relevant mode for Steps 2-7 (the operationally relevant one)
    slow_k = w5_mode_k
    ev_slow = ev_w5mode
    print(f"\n  Using W5-relevant mode (rate {rates[slow_k]:.6f}) for analysis."
          f"\n  (The absolute slowest at rate {rates[slowest_k]:.6f} has zero SE content,"
          f" see below.)")

    # Quick check: SE block of the absolute slowest mode
    B_slowest = np.zeros((N, N), dtype=complex)
    for m in range(N):
        for n in range(N):
            B_slowest[m, n] = R_inv[slowest_k, SE_IDX[m] * D + SE_IDX[n]]
    print(f"\n  SE block Frobenius norm of absolute slowest mode: "
          f"{np.linalg.norm(B_slowest, 'fro'):.2e}"
          f" ({'ZERO - inaccessible from W states' if np.linalg.norm(B_slowest, 'fro') < 1e-10 else 'nonzero'})")

    # ==================================================================
    # Step 2: Left eigenvector via R_inv
    # ==================================================================
    print("\n--- Step 2: Left eigenvector ---")
    # R_inv[k, :] is the k-th left co-vector.
    # c_slow = R_inv[slow_k, :] @ rho0.flatten()
    # V_left = R_inv[slow_k, :] reshaped to d x d
    V_left = R_inv[slow_k, :].reshape(D, D)
    V_right = R[:, slow_k].reshape(D, D)

    # Bi-orthogonality check: R_inv[k,:] @ R[:,k] should = 1
    biorth = R_inv[slow_k, :] @ R[:, slow_k]
    print(f"  Bi-orthogonality <L_k|R_k> = {biorth:.6f} (should be 1.0)")

    # Eigenvalue check: R_inv[k,:] @ L @ R[:,k] should = lambda_k
    H = heisenberg_H(N, 1.0)
    L = build_liouvillian(H, gamma_sacrifice)
    ev_check = R_inv[slow_k, :] @ L @ R[:, slow_k]
    ev_err = abs(ev_check - ev_slow)
    print(f"  Eigenvalue check: <L_k|L|R_k> = {ev_check.real:+.6f} {ev_check.imag:+.6f}j")
    print(f"  Error vs cached eigenvalue: {ev_err:.2e}")

    # ==================================================================
    # Step 3: Pauli decomposition of right and left eigenvectors
    # ==================================================================
    print("\n--- Step 3: Pauli decomposition ---")
    strings, P_mats = pauli_basis(N)
    P_flat = np.array([P.flatten() for P in P_mats])

    for tag, V_vec in [("RIGHT", R[:, slow_k]), ("LEFT", R_inv[slow_k, :])]:
        V = V_vec.reshape(D, D)
        # Pauli coefficients: c_P = Tr(P† V) / d
        # In flattened form: c_P = P_flat.conj() @ V_vec / d
        coeffs = P_flat.conj() @ V_vec / D
        probs = np.abs(coeffs)**2
        total = probs.sum()

        # Group by Hamming weight (number of non-I sites)
        hw_groups = {}
        for s, p in zip(strings, probs):
            hw = sum(1 for c in s if c != 'I')
            hw_groups[hw] = hw_groups.get(hw, 0) + p

        # Group by n_XY
        nxy_groups = {}
        for s, p in zip(strings, probs):
            nxy = n_xy_string(s)
            nxy_groups[nxy] = nxy_groups.get(nxy, 0) + p

        # Per-qubit weight
        qubit_wt = np.zeros(N)
        for s, p in zip(strings, probs):
            for k, c in enumerate(s):
                if c != 'I':
                    qubit_wt[k] += p
        qubit_wt_norm = qubit_wt / qubit_wt.sum() if qubit_wt.sum() > 0 else qubit_wt

        print(f"\n  {tag} eigenvector (total Pauli weight = {total:.6f}):")
        print(f"    Hamming weight groups:")
        for hw in sorted(hw_groups):
            print(f"      weight {hw}: {hw_groups[hw]/total*100:.2f}%")
        print(f"    n_XY groups:")
        for nxy in sorted(nxy_groups):
            print(f"      n_XY={nxy}: {nxy_groups[nxy]/total*100:.2f}%")
        print(f"    Per-qubit weight (normalized):")
        print(f"      [{', '.join(f'{w:.4f}' for w in qubit_wt_norm)}]")

    # ==================================================================
    # Step 4: Single-excitation block
    # ==================================================================
    print("\n--- Step 4: Single-excitation sector ---")

    # Build 5x5 block: B[m,n] = R_inv[slow_k, se_idx[m]*D + se_idx[n]]
    # c_slow = a† B^T a for single-excitation state psi = sum_k a_k |e_k>
    B = np.zeros((N, N), dtype=complex)
    for m in range(N):
        for n in range(N):
            B[m, n] = R_inv[slow_k, SE_IDX[m] * D + SE_IDX[n]]

    # Also the "natural" block from V_left reshaped
    V_left_block = np.zeros((N, N), dtype=complex)
    for m in range(N):
        for n in range(N):
            V_left_block[m, n] = V_left[SE_IDX[m], SE_IDX[n]]

    # They should be the same
    assert np.allclose(B, V_left_block), "Block extraction mismatch"
    print(f"  Block B = V_left_block check: OK")

    # Frobenius norms
    frob_block = np.linalg.norm(V_left_block, 'fro')
    frob_full = np.linalg.norm(V_left, 'fro')
    frob_ratio = frob_block / frob_full
    print(f"  ||V_left_block||_F = {frob_block:.6f}")
    print(f"  ||V_left_full||_F  = {frob_full:.6f}")
    print(f"  Frobenius ratio: {frob_ratio:.4f}")
    print(f"  Verdict: {'justified' if frob_ratio > 0.3 else 'marginal' if frob_ratio > 0.1 else 'wrong'}"
          f" (threshold: 0.3 justified, 0.1 marginal)")

    # ==================================================================
    # Step 5: Optimal single-excitation amplitudes
    # ==================================================================
    print("\n--- Step 5: Optimal amplitudes ---")

    # c_slow = a† B^T a  (derived from c_slow = R_inv[k,:] @ rho0.flatten())
    # where B[m,n] = V_left_block[m,n]
    Q = B.T  # The matrix in the quadratic form a† Q a

    # Method 1: Hermitian part
    M_H = (Q + Q.conj().T) / 2
    evals_H, evecs_H = np.linalg.eigh(M_H)
    # Take eigenvector with largest absolute eigenvalue
    idx_max = np.argmax(np.abs(evals_H))
    a_herm = evecs_H[:, idx_max]
    c_slow_herm = a_herm.conj() @ Q @ a_herm
    print(f"  Method 1 (Hermitian part of Q = B^T):")
    print(f"    Eigenvalues of M_H: {np.round(evals_H, 6)}")
    print(f"    Best eigenvalue: {evals_H[idx_max]:.6f} (index {idx_max})")
    print(f"    c_slow = {c_slow_herm:.6f},  |c_slow| = {abs(c_slow_herm):.6f}")

    # Method 2: SVD of Q
    U, s, Vh = np.linalg.svd(Q)
    a_svd = Vh[0].conj()  # right singular vector (maximizes ||Q a||)
    c_slow_svd = a_svd.conj() @ Q @ a_svd
    print(f"\n  Method 2 (SVD of Q):")
    print(f"    Singular values: {np.round(s, 6)}")
    print(f"    c_slow = {c_slow_svd:.6f},  |c_slow| = {abs(c_slow_svd):.6f}")

    # Pick the method with larger |c_slow|
    if abs(c_slow_herm) >= abs(c_slow_svd):
        a_opt = a_herm
        method_used = "Hermitian"
        print(f"\n  Selected: Hermitian method (|c_slow| = {abs(c_slow_herm):.6f})")
    else:
        a_opt = a_svd
        method_used = "SVD"
        print(f"\n  Selected: SVD method (|c_slow| = {abs(c_slow_svd):.6f})")

    # Ensure positive convention for comparison (make largest component real positive)
    phase = np.exp(-1j * np.angle(a_opt[np.argmax(np.abs(a_opt))]))
    a_opt = a_opt * phase

    print(f"\n  Optimal amplitudes a_opt (site 0..4):")
    for k in range(N):
        print(f"    site {k}: {a_opt[k].real:+.6f} {a_opt[k].imag:+.6f}j  "
              f"|a| = {abs(a_opt[k]):.6f}")

    # Construct psi_opt in the full 32-dim space
    psi_opt = np.zeros(D, dtype=complex)
    for k in range(N):
        psi_opt[SE_IDX[k]] = a_opt[k]
    psi_opt /= np.linalg.norm(psi_opt)

    # ==================================================================
    # Step 6: Compare to heuristic and Phase 2 winners
    # ==================================================================
    print("\n--- Step 6: Compare to heuristic ---")

    # Heuristic: a_k = sqrt(gamma_min / gamma_k)
    a_heur = np.sqrt(gamma_sacrifice.min() / gamma_sacrifice)
    a_heur = a_heur / np.linalg.norm(a_heur)

    # W5: uniform
    a_w5 = np.ones(N, dtype=complex) / np.sqrt(N)

    # W4 on {1,2,3,4}: zero on site 0
    a_w4 = np.zeros(N, dtype=complex)
    a_w4[1:] = 1.0 / 2.0
    a_w4 /= np.linalg.norm(a_w4)

    # W2 on {3,4}
    a_w2 = np.zeros(N, dtype=complex)
    a_w2[3] = a_w2[4] = 1.0 / np.sqrt(2)

    print(f"\n  Amplitude comparison (magnitudes):")
    print(f"  {'site':>4}  {'a_opt':>10}  {'a_heur':>10}  {'a_W5':>10}  {'ratio opt/heur':>14}")
    for k in range(N):
        ratio = abs(a_opt[k]) / a_heur[k] if a_heur[k] > 1e-10 else float('inf')
        print(f"  {k:>4}  {abs(a_opt[k]):>10.6f}  {a_heur[k]:>10.6f}  "
              f"{abs(a_w5[k]):>10.6f}  {ratio:>14.4f}")

    # Overlaps
    overlap_heur = abs(a_opt.conj() @ a_heur)
    overlap_w5 = abs(a_opt.conj() @ a_w5)
    overlap_w4 = abs(a_opt.conj() @ a_w4)
    overlap_w2 = abs(a_opt.conj() @ a_w2)

    print(f"\n  Overlaps |<a_opt | a_X>|:")
    print(f"    a_heuristic (sqrt(gmin/gk)): {overlap_heur:.6f}")
    print(f"    a_W5 (uniform):              {overlap_w5:.6f}")
    print(f"    a_W4_1234:                   {overlap_w4:.6f}")
    print(f"    a_W2_34:                     {overlap_w2:.6f}")

    # Slow-mode projection for each
    print(f"\n  Slow-mode projection |c_slow| for each state:")
    for name, a in [("a_opt", a_opt), ("a_heur", a_heur), ("a_W5", a_w5),
                     ("a_W4_1234", a_w4), ("a_W2_34", a_w2)]:
        c = a.conj() @ Q @ a
        print(f"    {name:<12}: |c_slow| = {abs(c):.6f}")

    # ==================================================================
    # Step 7: Time evolution verification
    # ==================================================================
    print("\n--- Step 7: Time evolution verification ---")
    times = np.linspace(0, 30, 301)

    # psi_opt
    traj_opt = concurrence_trajectory(psi_opt, prof, times)
    auc_opt = auc_windows(times, traj_opt)
    rho_opt = np.outer(psi_opt, psi_opt.conj())
    c_init_opt = float(max_adjacent_concurrence(rho_opt, N))

    # Sacrifice-tuned W5 (for comparison)
    psi_heur = np.zeros(D, dtype=complex)
    for k in range(N):
        psi_heur[SE_IDX[k]] = a_heur[k]
    psi_heur /= np.linalg.norm(psi_heur)
    traj_heur = concurrence_trajectory(psi_heur, prof, times)
    auc_heur = auc_windows(times, traj_heur)
    c_init_heur = float(max_adjacent_concurrence(np.outer(psi_heur, psi_heur.conj()), N))

    # W5
    psi_w5 = np.zeros(D, dtype=complex)
    for k in range(N):
        psi_w5[SE_IDX[k]] = a_w5[k]
    psi_w5 /= np.linalg.norm(psi_w5)
    traj_w5 = concurrence_trajectory(psi_w5, prof, times)
    auc_w5 = auc_windows(times, traj_w5)
    c_init_w5 = float(max_adjacent_concurrence(np.outer(psi_w5, psi_w5.conj()), N))

    print(f"\n  {'State':<22}  {'C_init':>6}  {'AUC(2)':>8}  {'AUC(10)':>8}  {'AUC(30)':>8}")
    print(f"  {'-' * 58}")
    for name, auc, ci in [("psi_opt (lens)", auc_opt, c_init_opt),
                           ("sacrifice_tuned_W5", auc_heur, c_init_heur),
                           ("W5_full", auc_w5, c_init_w5)]:
        print(f"  {name:<22}  {ci:>6.3f}  {auc['auc_2']:>8.4f}  "
              f"{auc['auc_10']:>8.4f}  {auc['auc_30']:>8.4f}")

    beats_heur = auc_opt['auc_10'] > auc_heur['auc_10'] * 1.01
    ties_heur = abs(auc_opt['auc_10'] - auc_heur['auc_10']) / auc_heur['auc_10'] < 0.01

    if beats_heur:
        print(f"\n  RESULT: psi_opt BEATS sacrifice_tuned_W5 by "
              f"{(auc_opt['auc_10']/auc_heur['auc_10'] - 1)*100:.1f}%")
    elif ties_heur:
        print(f"\n  RESULT: psi_opt TIES sacrifice_tuned_W5 within 1%")
    else:
        print(f"\n  RESULT: psi_opt LOSES to sacrifice_tuned_W5 by "
              f"{(1 - auc_opt['auc_10']/auc_heur['auc_10'])*100:.1f}%")

    # ==================================================================
    # Save results
    # ==================================================================
    out_path = os.path.join(out_dir, "slow_mode_lens_results.txt")
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("Slow Mode Lens Analysis Results\n")
        f.write(f"Computed: {np.datetime64('today')}\n\n")
        f.write(f"Slow mode: index {slow_k}, eigenvalue {ev_slow.real:+.8f} {ev_slow.imag:+.8f}j\n")
        f.write(f"Bi-orthogonality check: {biorth:.8f}\n")
        f.write(f"Eigenvalue check error: {ev_err:.2e}\n\n")
        f.write(f"Frobenius ratio (SE block / full): {frob_ratio:.6f}\n")
        f.write(f"Method used: {method_used}\n\n")
        f.write("Optimal amplitudes:\n")
        for k in range(N):
            f.write(f"  site {k}: {a_opt[k].real:+.8f} {a_opt[k].imag:+.8f}j  "
                    f"|a| = {abs(a_opt[k]):.8f}\n")
        f.write(f"\nOverlap |<a_opt | a_heuristic>| = {overlap_heur:.8f}\n")
        f.write(f"Overlap |<a_opt | a_W5>|          = {overlap_w5:.8f}\n\n")
        f.write("AUC comparison:\n")
        for name, auc, ci in [("psi_opt (lens)", auc_opt, c_init_opt),
                               ("sacrifice_tuned_W5", auc_heur, c_init_heur),
                               ("W5_full", auc_w5, c_init_w5)]:
            f.write(f"  {name:<22}  C_init={ci:.4f}  AUC(2)={auc['auc_2']:.4f}"
                    f"  AUC(10)={auc['auc_10']:.4f}  AUC(30)={auc['auc_30']:.4f}\n")

    print(f"\n  Saved {out_path}")
    print(f"\nTotal runtime: {_time.time() - t_start:.1f}s")
