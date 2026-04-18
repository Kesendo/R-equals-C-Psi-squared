#!/usr/bin/env python3
"""observer_time_slow_mode_analysis.py

Slow-mode Liouvillian eigendecomposition + first-order f_i prediction
for the N=7 observer-time investigation.

Pipeline:
  1. Build sparse Liouvillian L_A (uniform J) and the perturbation V_L
     such that L_B(dJ) ≈ L_A + dJ · V_L (dJ = J_mod - 1 on bond (0, 1)).
  2. Find slow modes of L_A: scipy.sparse.linalg.eigs with sigma = 0,
     k = 50 (complex eigenvalues closest to zero).
  3. For each slow mode, compute site profile (single-site marginal,
     Pauli decomposition per site), initial-state overlap with rho_0
     (for psi_1 and psi_2 initial states).
  4. First-order eigenvalue shift: delta_lambda_s = < W_s | V_L | M_s >
     where W_s is the left eigenvector, M_s the right eigenvector.
  5. Analytical f_i prediction based on the dominant slow mode per site.
  6. Extract empirical f_i from existing .npz data (J_mod = 0.9, 1.1) and
     compare.

Outputs:
  simulations/results/observer_time_slow_mode/
    slow_modes_A.json                 - eigenvalues + site profile + overlaps
    f_i_prediction.json               - analytical prediction vs empirical
    analysis_log.txt                  - human-readable log

Date: 2026-04-18
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).parent))
from n7_coupling_defect_overlay import (
    N, GAMMA_0, single_excitation_mode, vacuum,
)

RESULTS_DIR = Path(__file__).parent / "results" / "observer_time_slow_mode"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

D = 2**N
DD = D * D
SLOW_CUTOFF = 2.0 * GAMMA_0    # slow = |Re(lambda)| < cutoff
K_EIGS = 80                    # number of slow eigenvalues to compute


# ---------------------------------------------------------------------------
# Sparse Hamiltonian + Liouvillian
# ---------------------------------------------------------------------------
def build_sparse_H_XY(J_list, N):
    """Sparse H = sum_i (J_i / 2) (X_i X_{i+1} + Y_i Y_{i+1}) in the
    computational basis. Non-zero entries: for each pair (i, i+1), add a
    hopping term between basis states differing only by swap of bits i, i+1.

    Big-endian convention: site 0 is the most significant bit
    (state index has bit N-1-site)."""
    rows, cols, vals = [], [], []
    for i in range(N - 1):
        J_i = J_list[i]
        bit_i = 1 << (N - 1 - i)
        bit_j = 1 << (N - 1 - (i + 1))
        mask = bit_i | bit_j
        for a in range(D):
            two_bits = a & mask
            if two_bits == bit_i:      # bit i = 1, bit j = 0 -> |..10..>
                b = a ^ mask            # flip both -> |..01..>
                rows.append(a); cols.append(b); vals.append(J_i)
                rows.append(b); cols.append(a); vals.append(J_i)
    return sps.csr_matrix((vals, (rows, cols)), shape=(D, D), dtype=complex)


def build_sparse_liouvillian(J_list, N, gamma_0):
    """L in vec-row-major convention: vec(rho)[a*d + b] = rho[a, b].
    L = -i (kron(H, I) - kron(I, H^T)) + diag(-2*gamma_0 * popcount(a XOR b)).
    """
    H = build_sparse_H_XY(J_list, N)
    Id = sps.eye(D, dtype=complex, format='csr')
    L_h = -1j * (sps.kron(H, Id, format='csr')
                 - sps.kron(Id, H.T, format='csr'))
    a_idx = np.arange(DD) // D
    b_idx = np.arange(DD) % D
    h = np.array([bin(int(av ^ bv)).count('1')
                  for av, bv in zip(a_idx, b_idx)])
    deph_diag = -2.0 * gamma_0 * h
    L_d = sps.diags(deph_diag.astype(complex), format='csr')
    return (L_h + L_d).tocsr()


def build_perturbation_V_L(N):
    """V_L such that L(J_mod) = L(1.0) + (J_mod - 1) * V_L.
    Only the bond (0, 1) term in H depends on J_mod, so V_L = -i (kron(V, I)
    - kron(I, V^T)) where V is the H contribution at J_{0,1} with J = 1."""
    V = build_sparse_H_XY([1.0] + [0.0] * (N - 2), N)
    Id = sps.eye(D, dtype=complex, format='csr')
    return -1j * (sps.kron(V, Id, format='csr')
                  - sps.kron(Id, V.T, format='csr')).tocsr()


# ---------------------------------------------------------------------------
# Slow-mode extraction (shift-invert around sigma = 0)
# ---------------------------------------------------------------------------
def slow_modes(L, k=K_EIGS, sigma=-1e-3):
    """sigma slightly negative to avoid the exact-zero stationary manifold
    that makes (L - 0*I) singular."""
    vals, vecs = spla.eigs(L, k=k, sigma=sigma, which='LM',
                           tol=1e-10, maxiter=5000)
    order = np.argsort(-vals.real)
    return vals[order], vecs[:, order]


def left_eigenvectors(L, vals, right_vecs, tol=1e-6):
    """Compute left eigenvectors W_s for each eigenvalue v_s by sparse
    shift-invert on L^H with sigma = conj(v_s). Then normalize biorthogonal
    so that W_s^T * V_s = delta_{s s'}. """
    DD_local = L.shape[0]
    W = np.zeros_like(right_vecs)
    # For each val, solve L^H x = conj(val) x near conj(val).
    # Simpler: use dense LU on (L^H - sigma I) for a moderate k.
    # We'll use spla.eigs on L.getH() around conj(v_s).
    # Batch by groups of unique vals; for first pass, just loop.
    LH = L.conj().T.tocsr()
    for s, v in enumerate(vals):
        try:
            w_vals, w_vecs = spla.eigs(LH, k=1,
                                        sigma=np.conj(v), which='LM',
                                        tol=1e-10, maxiter=3000)
            W[:, s] = w_vecs[:, 0]
        except Exception:
            W[:, s] = 0.0
    # Biorthogonalize: scale each W[:, s] so W[:, s] @ right[:, s] = 1.
    for s in range(len(vals)):
        z = W[:, s] @ right_vecs[:, s]
        if abs(z) > 1e-12:
            W[:, s] /= z
    return W


# ---------------------------------------------------------------------------
# Site profiles and initial-state overlaps
# ---------------------------------------------------------------------------
def mode_to_rho(v):
    """Reshape a vec-eigenvector (length DD) to a dxd operator."""
    return v.reshape(D, D)


_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

def reduced_single(rho):
    """Trace all but qubit 0 out of an arbitrary rho (not necessarily
    density matrix). Returns a list of 2x2 marginals, one per qubit."""
    out = []
    t = rho.reshape([2] * (2 * N))
    row = [_LETTERS[q] for q in range(N)]
    col = [_LETTERS[q + N] for q in range(N)]
    for q in range(N):
        r = list(row); c = list(col)
        for p in range(N):
            if p != q:
                c[p] = r[p]
        spec = ''.join(r) + ''.join(c) + '->' + r[q] + c[q]
        out.append(np.einsum(spec, t))
    return out


def pauli_decomposition(rho_i):
    """Decompose a 2x2 operator rho_i = 0.5 * (a_I I + a_X X + a_Y Y + a_Z Z).
    Returns (a_I, a_X, a_Y, a_Z)."""
    a_I = float(np.real(np.trace(rho_i)))
    a_X = float(np.real(np.trace(rho_i @ np.array([[0, 1], [1, 0]],
                                                  dtype=complex))))
    a_Y = float(np.real(np.trace(rho_i @ np.array([[0, -1j], [1j, 0]],
                                                  dtype=complex))))
    a_Z = float(np.real(np.trace(rho_i @ np.array([[1, 0], [0, -1]],
                                                  dtype=complex))))
    return (a_I, a_X, a_Y, a_Z)


def site_profile_norm(rho):
    """Single number per site quantifying how much of the operator rho is
    localised at qubit q. Use the L2 norm of the Pauli vector (X, Y, Z
    components only; the Identity component is trivial)."""
    out = np.zeros(N)
    marginals = reduced_single(rho)
    for q in range(N):
        aI, aX, aY, aZ = pauli_decomposition(marginals[q])
        out[q] = float(np.sqrt(aX * aX + aY * aY + aZ * aZ))
    return out


# ---------------------------------------------------------------------------
# Initial state construction (as vectorized rho_0)
# ---------------------------------------------------------------------------
def initial_state_vec(k):
    """rho_0 = |phi><phi| with phi = (|vac> + |psi_k>)/sqrt(2)."""
    phi = vacuum(N) + single_excitation_mode(N, k)
    phi /= np.linalg.norm(phi)
    rho = np.outer(phi, phi.conj())
    return rho.reshape(DD)


# ---------------------------------------------------------------------------
# First-order prediction of f_i via dominant-mode picture
# ---------------------------------------------------------------------------
def predict_fi_from_slow_modes(vals, right_vecs, left_vecs, V_L, rho_0_vec,
                               site_profiles, stationary_eps=1e-5):
    """Analytical first-order f_i prediction from the dominant slowly-
    decaying mode per site. Stationary modes (|Re(lambda)| < eps) are
    excluded: they do not drive the purity transient and therefore do
    not contribute to the alpha rescale."""
    overlaps = np.array([complex(left_vecs[:, s] @ rho_0_vec)
                         for s in range(len(vals))])
    delta_lambda = np.array([complex(left_vecs[:, s] @ (V_L @ right_vecs[:, s]))
                             for s in range(len(vals))])
    fi_pred = np.full(N, np.nan)
    dominant = [None] * N
    decaying = np.abs(vals.real) > stationary_eps
    for i in range(N):
        # Score: |c_s| * |profile[i]| for decaying modes only.
        score = (np.abs(overlaps) * np.abs(site_profiles[:, i])
                 * decaying.astype(float))
        if score.max() < 1e-12:
            continue
        best = int(np.argmax(score))
        dominant[i] = (best, float(vals[best].real), float(vals[best].imag),
                       float(np.abs(overlaps[best])),
                       float(site_profiles[best, i]),
                       float(score[best]))
        lam = vals[best]
        dlam = delta_lambda[best]
        if abs(lam.real) > stationary_eps:
            fi_pred[i] = float(dlam.real / lam.real)
    return fi_pred, dominant, overlaps, delta_lambda


# ---------------------------------------------------------------------------
# Empirical f_i from existing .npz data
# ---------------------------------------------------------------------------
def empirical_fi_from_dir(data_dir):
    d = Path(data_dir)
    t = np.load(d / "times.npy")
    p_A = np.load(d / "experiment_A.npz")['purity']
    # Load 0.9 and 1.1
    try:
        p_09 = np.load(d / "experiment_B_0.9.npz")['purity']
        p_11 = np.load(d / "experiment_B_1.1.npz")['purity']
    except FileNotFoundError:
        return None
    # alpha-fit at 0.9 and 1.1 via simple bounded search
    from scipy.interpolate import interp1d
    from scipy.optimize import minimize_scalar

    fi = np.empty(N)
    alphas_09 = np.empty(N)
    alphas_11 = np.empty(N)
    t_max = 20.0
    mask = t <= t_max
    for i in range(N):
        interp_A = interp1d(t, p_A[:, i], bounds_error=False,
                            fill_value=(p_A[0, i], p_A[-1, i]), kind='cubic')
        def mse_09(a):
            d = interp_A(a * t[mask]) - p_09[mask, i]
            return float(np.mean(d * d))
        def mse_11(a):
            d = interp_A(a * t[mask]) - p_11[mask, i]
            return float(np.mean(d * d))
        r09 = minimize_scalar(mse_09, bounds=(0.1, 10.0), method='bounded',
                              options={'xatol': 1e-6})
        r11 = minimize_scalar(mse_11, bounds=(0.1, 10.0), method='bounded',
                              options={'xatol': 1e-6})
        alphas_09[i] = float(r09.x)
        alphas_11[i] = float(r11.x)
        fi[i] = (r11.x - r09.x) / 0.2
    return dict(alphas_09=alphas_09, alphas_11=alphas_11, fi=fi)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log_lines = []
    def log(msg=""):
        print(msg, flush=True)
        log_lines.append(msg)

    log("=" * 72)
    log("SLOW-MODE ANALYSIS (N = 7, gamma_0 = 0.05, defect bond (0, 1))")
    log("=" * 72)
    log(f"  Hilbert dim D = {D}, Liouvillian dim = {DD}")
    log(f"  Slow cutoff |Re(lambda)| < {SLOW_CUTOFF}, k = {K_EIGS}")
    log()

    # Build sparse L_A and V_L
    t0 = time.time()
    L_A = build_sparse_liouvillian([1.0] * (N - 1), N, GAMMA_0)
    V_L = build_perturbation_V_L(N)
    log(f"  Built sparse L_A  : {L_A.nnz} nonzeros, {time.time() - t0:.2f} s")

    # Find slow modes of L_A (sigma slightly negative to bypass the exact-0
    # stationary manifold)
    t0 = time.time()
    vals_A, right_A = slow_modes(L_A, k=K_EIGS, sigma=-1e-3)
    log(f"  eigs(L_A) for k = {K_EIGS} slow modes: {time.time() - t0:.2f} s")
    log(f"    Re(lambda) range: [{vals_A.real.min():.4f}, {vals_A.real.max():.4f}]")
    log(f"    Number with |Re(lambda)| < {SLOW_CUTOFF}: "
        f"{int(np.sum(np.abs(vals_A.real) < SLOW_CUTOFF))}")

    # Count stationary modes (|Re(lambda)| < 1e-8)
    n_stationary = int(np.sum(np.abs(vals_A.real) < 1e-8))
    log(f"    Strict stationary (|Re| < 1e-8): {n_stationary}")

    # Left eigenvectors (biorthogonal)
    t0 = time.time()
    left_A = left_eigenvectors(L_A, vals_A, right_A)
    log(f"  Left eigenvectors + biorthogonalisation: {time.time() - t0:.2f} s")

    # Site profiles for each slow mode
    t0 = time.time()
    site_profiles = np.zeros((len(vals_A), N))
    for s in range(len(vals_A)):
        rho_s = mode_to_rho(right_A[:, s])
        site_profiles[s] = site_profile_norm(rho_s)
    log(f"  Site profiles for {len(vals_A)} slow modes: "
        f"{time.time() - t0:.2f} s")

    # First-order prediction for both psi_1 and psi_2 initial states
    log()
    log("=" * 72)
    log("f_i PREDICTION (first-order dominant-slow-mode picture)")
    log("=" * 72)

    preds = {}
    for k in (1, 2):
        log()
        log(f"--- Initial state: phi = (|vac> + |psi_{k}>)/sqrt(2) ---")
        rho_0_vec = initial_state_vec(k)
        fi_pred, dominant, overlaps, dlam = predict_fi_from_slow_modes(
            vals_A, right_A, left_A, V_L, rho_0_vec, site_profiles)
        log(f"  dominant mode per site (site, mode_idx, Re(lambda), Im(lambda), "
            f"|overlap|, profile[i], score):")
        for i in range(N):
            m, rl, il, ov, pf, sc = dominant[i]
            log(f"    site {i}: mode {m:3d}  lambda = {rl:+.4f} "
                f"{il:+.4f}i   |<W|rho_0>| = {ov:.3e}   "
                f"profile = {pf:.3f}   score = {sc:.3e}")
        log(f"  predicted f_i: {fi_pred}")
        preds[k] = dict(fi_pred=fi_pred.tolist(),
                        dominant=[list(x) for x in dominant])

    # Compare against empirical
    log()
    log("=" * 72)
    log("EMPIRICAL vs PREDICTED f_i")
    log("=" * 72)
    # ψ_1 empirical: from finer_Jmod scan
    log()
    log("--- Initial state psi_1 ---")
    emp1 = empirical_fi_from_dir(
        Path(__file__).parent / "results"
        / "n7_coupling_defect_overlay_extended" / "finer_Jmod")
    if emp1 is not None:
        log(f"  empirical f_i (from J_mod = 0.9, 1.1 alpha-fit, "
            f"finer_Jmod scan):")
        log(f"    {emp1['fi']}")
        log(f"  predicted f_i (dominant-mode, psi_1 init):")
        log(f"    {np.array(preds[1]['fi_pred'])}")
    else:
        log(f"  (finer_Jmod data not found)")

    # ψ_2 empirical: from psi2_init scan
    log()
    log("--- Initial state psi_2 ---")
    emp2 = empirical_fi_from_dir(
        Path(__file__).parent / "results"
        / "n7_coupling_defect_overlay_extended" / "psi2_init")
    if emp2 is not None:
        log(f"  empirical f_i (from J_mod = 0.9, 1.1, psi2 scan):")
        log(f"    {emp2['fi']}")
        log(f"  predicted f_i (dominant-mode, psi_2 init):")
        log(f"    {np.array(preds[2]['fi_pred'])}")
    else:
        log(f"  (psi2_init data not found)")

    # Save summary
    summary = {
        'N': N, 'gamma_0': GAMMA_0,
        'lambda_range': [float(vals_A.real.min()),
                         float(vals_A.real.max())],
        'n_slow': int(np.sum(np.abs(vals_A.real) < SLOW_CUTOFF)),
        'n_stationary': n_stationary,
        'eigenvalues': [[float(v.real), float(v.imag)] for v in vals_A],
        'predictions': preds,
        'empirical_psi1': {k: v.tolist() for k, v in (emp1 or {}).items()},
        'empirical_psi2': {k: v.tolist() for k, v in (emp2 or {}).items()},
    }
    with open(RESULTS_DIR / 'slow_modes_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    (RESULTS_DIR / 'analysis_log.txt').write_text("\n".join(log_lines) + "\n",
                                                 encoding='utf-8')
    log()
    log(f"Results: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
