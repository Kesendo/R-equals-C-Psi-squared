"""V_inter SVD with R-parity decomposition: the correct F86-style L_eff EP test.

V_inter = P_{HD=1}^† M_H_total P_{HD=3} is the HD=1 × HD=3 cross-block of L_H.
R-symmetric (since [L_H, R] = 0), so V_inter block-diagonalizes by R:

    V_inter = V_inter^{++} ⊕ V_inter^{--}

SVD each R-block independently:
- top σ_0⁺ with vectors u_0⁺ ∈ R-even HD=1, v_0⁺ ∈ R-even HD=3
- top σ_0⁻ with vectors u_0⁻ ∈ R-odd HD=1, v_0⁻ ∈ R-odd HD=3

For each R-parity, the 2×2 L_eff on span{u_0, v_0} has structure:

    L_eff = [[ −2γ, iJ·σ_0 ],
             [ iJ·σ_0, −6γ ]]

EP at J·σ_0 = 2γ, i.e. Q_EP = 2/σ_0.

If σ_0⁺ ≠ σ_0⁻, R-parity-resolved EP positions exist: a genuinely new diagnostic.

Run: python simulations/v_inter_svd_r_parity.py [N]
"""
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
gamma = 0.05
dim = 2 ** N
op_dim = dim * dim

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


def build_H(N, J):
    H = np.zeros((dim, dim), dtype=complex)
    for l in range(N - 1):
        for axis_idx in (1, 2):
            op_alpha = [0] * N
            op_alpha[l] = axis_idx
            op_alpha[l + 1] = axis_idx
            H += (J / 2) * pauli_string(op_alpha)
    return H


def build_L(N, gamma, J):
    H = build_H(N, J)
    I_dim = np.eye(dim, dtype=complex)
    L_H = -1j * (np.kron(I_dim, H) - np.kron(H.T, I_dim))
    L_D = np.zeros((op_dim, op_dim), dtype=complex)
    for l in range(N):
        op_alpha = [0] * N
        op_alpha[l] = 3
        sig_z_l = pauli_string(op_alpha)
        L_D += gamma * (np.kron(sig_z_l, sig_z_l) - np.eye(op_dim, dtype=complex))
    return L_H + L_D


def build_hd_basis(N, hd):
    indices = []
    for a in range(dim):
        for b in range(dim):
            if bin(a ^ b).count('1') == hd:
                indices.append(b * dim + a)
    P = np.zeros((op_dim, len(indices)), dtype=float)
    for i, idx in enumerate(indices):
        P[idx, i] = 1.0
    return P


print(f"# N = {N}, γ₀ = {gamma}")
print(f"# Building V_inter at reference J = 1.0")

J_ref = 1.0
H_ref = build_H(N, J_ref)
I_dim = np.eye(dim, dtype=complex)
M_H_ref = -1j * (np.kron(I_dim, H_ref) - np.kron(H_ref.T, I_dim))

P_hd1 = build_hd_basis(N, 1)
P_hd3 = build_hd_basis(N, 3)
print(f"#   HD=1 subspace dim: {P_hd1.shape[1]}")
print(f"#   HD=3 subspace dim: {P_hd3.shape[1]}")

V_inter = P_hd1.T @ M_H_ref @ P_hd3
print(f"#   V_inter shape: {V_inter.shape}")

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

print(f"#   R-decomp HD=1: even {even_1.shape[1]}, odd {odd_1.shape[1]}")
print(f"#   R-decomp HD=3: even {even_3.shape[1]}, odd {odd_3.shape[1]}")

V_pp = even_1.T @ V_inter @ even_3
V_mm = odd_1.T @ V_inter @ odd_3
V_pm = even_1.T @ V_inter @ odd_3
V_mp = odd_1.T @ V_inter @ even_3

print(f"\n# V_inter R-block Frobenius norms (verify R-symmetry of V_inter):")
print(f"#   ‖V_++‖_F = {np.linalg.norm(V_pp):.6f}")
print(f"#   ‖V_--‖_F = {np.linalg.norm(V_mm):.6f}")
print(f"#   ‖V_+-‖_F = {np.linalg.norm(V_pm):.6e}   (should ~ 0)")
print(f"#   ‖V_-+‖_F = {np.linalg.norm(V_mp):.6e}   (should ~ 0)")

U_pp, s_pp, Vh_pp = np.linalg.svd(V_pp)
U_mm, s_mm, Vh_mm = np.linalg.svd(V_mm)

print(f"\n# Top 8 singular values of V_++: {[f'{s:.4f}' for s in s_pp[:8]]}")
print(f"# Top 8 singular values of V_--: {[f'{s:.4f}' for s in s_mm[:8]]}")

if s_pp[0] > 1e-9:
    Q_EP_p = 2.0 / s_pp[0]
    print(f"\n# Predicted Q_EP (R-even) = 2/σ_0⁺ = 2/{s_pp[0]:.4f} = {Q_EP_p:.4f}")
else:
    print(f"\n# R-even top σ_0 ≈ 0; no EP from this pair")

if s_mm[0] > 1e-9:
    Q_EP_m = 2.0 / s_mm[0]
    print(f"# Predicted Q_EP (R-odd)  = 2/σ_0⁻ = 2/{s_mm[0]:.4f} = {Q_EP_m:.4f}")
else:
    print(f"# R-odd top σ_0 ≈ 0; no EP from this pair")

# Lift top singular vectors to op-space and normalize
u_top_p = P_hd1 @ even_1 @ U_pp[:, 0]
v_top_p = P_hd3 @ even_3 @ Vh_pp[0, :].conj()
u_top_p = u_top_p / np.linalg.norm(u_top_p)
v_top_p = v_top_p / np.linalg.norm(v_top_p)

u_top_m = P_hd1 @ odd_1 @ U_mm[:, 0]
v_top_m = P_hd3 @ odd_3 @ Vh_mm[0, :].conj()
u_top_m = u_top_m / np.linalg.norm(u_top_m)
v_top_m = v_top_m / np.linalg.norm(v_top_m)

# Verify R-parity of lifted vectors
print(f"\n# R-eigenvalue verification (should be ±1):")
print(f"#   ⟨u⁺|R|u⁺⟩ = {(u_top_p.conj() @ R_op @ u_top_p).real:+.4f}")
print(f"#   ⟨v⁺|R|v⁺⟩ = {(v_top_p.conj() @ R_op @ v_top_p).real:+.4f}")
print(f"#   ⟨u⁻|R|u⁻⟩ = {(u_top_m.conj() @ R_op @ u_top_m).real:+.4f}")
print(f"#   ⟨v⁻|R|v⁻⟩ = {(v_top_m.conj() @ R_op @ v_top_m).real:+.4f}")

# Q-sweep with 2×2 L_eff per parity
print(f"\n# Q-sweep with 2×2 L_eff per R-parity")
print(f"#   Q    | R-even λ_+ (Re,Im)   λ_- (Re,Im)     |ΔRe| | R-odd λ_+ (Re,Im)   λ_- (Re,Im)     |ΔRe|")
print(f"#   -----|----------------------------------------------|----------------------------------------------")

Q_grid = list(np.linspace(0.3, 1.5, 13)) + list(np.linspace(1.55, 3.0, 30))
for Q in Q_grid:
    J = Q * gamma
    L = build_L(N, gamma, J)
    L_eff_p = np.array([
        [u_top_p.conj() @ L @ u_top_p, u_top_p.conj() @ L @ v_top_p],
        [v_top_p.conj() @ L @ u_top_p, v_top_p.conj() @ L @ v_top_p],
    ])
    L_eff_m = np.array([
        [u_top_m.conj() @ L @ u_top_m, u_top_m.conj() @ L @ v_top_m],
        [v_top_m.conj() @ L @ u_top_m, v_top_m.conj() @ L @ v_top_m],
    ])
    eigs_p = np.linalg.eigvals(L_eff_p)
    eigs_m = np.linalg.eigvals(L_eff_m)
    re_gap_p = abs(eigs_p[0].real - eigs_p[1].real)
    re_gap_m = abs(eigs_m[0].real - eigs_m[1].real)
    e_p = sorted(eigs_p, key=lambda e: -e.real)
    e_m = sorted(eigs_m, key=lambda e: -e.real)
    print(f"   {Q:.2f} | ({e_p[0].real:+.4f},{e_p[0].imag:+.4f}) ({e_p[1].real:+.4f},{e_p[1].imag:+.4f}) "
          f"{re_gap_p:.4f} | ({e_m[0].real:+.4f},{e_m[0].imag:+.4f}) ({e_m[1].real:+.4f},{e_m[1].imag:+.4f}) {re_gap_m:.4f}")
