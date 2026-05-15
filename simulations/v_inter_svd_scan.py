"""Quick V_inter SVD R-parity scan across N (no Q-sweep, just SVD).

Reports σ_0⁺, σ_0⁻ for the top singular values of V_inter restricted to each
R-parity block, at multiple N. Q_EP per parity = 2/σ_0.

Run: python simulations/_v_inter_svd_scan.py
"""
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

PAULI = np.array([
    [[1, 0], [0, 1]],
    [[0, 1], [1, 0]],
    [[0, -1j], [1j, 0]],
    [[1, 0], [0, -1]],
], dtype=complex)


def pauli_string(alpha):
    result = PAULI[alpha[0]]
    for a in alpha[1:]:
        result = np.kron(result, PAULI[a])
    return result


def reverse_bits(b, n):
    out = 0
    for i in range(n):
        if (b >> i) & 1:
            out |= 1 << (n - 1 - i)
    return out


def analyse(N, J_ref=1.0):
    dim = 2 ** N
    op_dim = dim * dim

    H = np.zeros((dim, dim), dtype=complex)
    for l in range(N - 1):
        for axis_idx in (1, 2):
            op_alpha = [0] * N
            op_alpha[l] = axis_idx
            op_alpha[l + 1] = axis_idx
            H += (J_ref / 2) * pauli_string(op_alpha)
    I_dim = np.eye(dim, dtype=complex)
    M_H = -1j * (np.kron(I_dim, H) - np.kron(H.T, I_dim))

    def hd_basis(hd):
        indices = []
        for a in range(dim):
            for b in range(dim):
                if bin(a ^ b).count('1') == hd:
                    indices.append(b * dim + a)
        P = np.zeros((op_dim, len(indices)), dtype=float)
        for i, idx in enumerate(indices):
            P[idx, i] = 1.0
        return P

    P_hd1 = hd_basis(1)
    P_hd3 = hd_basis(3)

    R_h = np.zeros((dim, dim), dtype=float)
    for b in range(dim):
        R_h[reverse_bits(b, N), b] = 1.0
    R_op = np.kron(R_h, R_h)

    R_hd1 = P_hd1.T @ R_op @ P_hd1
    R_hd3 = P_hd3.T @ R_op @ P_hd3
    r_vals_1, r_vecs_1 = np.linalg.eigh(R_hd1)
    r_vals_3, r_vecs_3 = np.linalg.eigh(R_hd3)
    even_1 = r_vecs_1[:, np.abs(r_vals_1 - 1) < 1e-6]
    odd_1 = r_vecs_1[:, np.abs(r_vals_1 + 1) < 1e-6]
    even_3 = r_vecs_3[:, np.abs(r_vals_3 - 1) < 1e-6]
    odd_3 = r_vecs_3[:, np.abs(r_vals_3 + 1) < 1e-6]

    V_inter = P_hd1.T @ M_H @ P_hd3
    V_pp = even_1.T @ V_inter @ even_3
    V_mm = odd_1.T @ V_inter @ odd_3
    V_pm = even_1.T @ V_inter @ odd_3
    V_mp = odd_1.T @ V_inter @ even_3

    s_pp = np.linalg.svd(V_pp, compute_uv=False)
    s_mm = np.linalg.svd(V_mm, compute_uv=False)

    return {
        "N": N,
        "hd1_dim": P_hd1.shape[1],
        "hd3_dim": P_hd3.shape[1],
        "even_hd1": even_1.shape[1],
        "odd_hd1": odd_1.shape[1],
        "even_hd3": even_3.shape[1],
        "odd_hd3": odd_3.shape[1],
        "V_pp_norm": np.linalg.norm(V_pp),
        "V_mm_norm": np.linalg.norm(V_mm),
        "V_cross_norm": np.linalg.norm(V_pm) + np.linalg.norm(V_mp),
        "s_pp": s_pp,
        "s_mm": s_mm,
    }


for N in [3, 4, 5, 6]:
    r = analyse(N)
    print(f"\n# N = {N}")
    print(f"#   HD=1 dim: {r['hd1_dim']}, HD=3 dim: {r['hd3_dim']}")
    print(f"#   R-even HD=1: {r['even_hd1']}, R-odd HD=1: {r['odd_hd1']}")
    print(f"#   R-even HD=3: {r['even_hd3']}, R-odd HD=3: {r['odd_hd3']}")
    print(f"#   ‖V_++‖_F = {r['V_pp_norm']:.4f}, ‖V_--‖_F = {r['V_mm_norm']:.4f}")
    print(f"#   ‖V_cross‖_F = {r['V_cross_norm']:.2e}   (should ~ 0)")
    n_show = min(8, len(r['s_pp']))
    print(f"#   Top {n_show} σ of V_++:  {[f'{s:.4f}' for s in r['s_pp'][:n_show]]}")
    n_show_m = min(8, len(r['s_mm']))
    print(f"#   Top {n_show_m} σ of V_--:  {[f'{s:.4f}' for s in r['s_mm'][:n_show_m]]}")
    if r['s_pp'][0] > 1e-9:
        print(f"#   Q_EP (R-even) = 2/{r['s_pp'][0]:.4f} = {2/r['s_pp'][0]:.4f}")
    if r['s_mm'][0] > 1e-9:
        print(f"#   Q_EP (R-odd)  = 2/{r['s_mm'][0]:.4f} = {2/r['s_mm'][0]:.4f}")
    delta = abs(r['s_pp'][0] - r['s_mm'][0])
    print(f"#   |σ_0⁺ − σ_0⁻| = {delta:.6e}")
