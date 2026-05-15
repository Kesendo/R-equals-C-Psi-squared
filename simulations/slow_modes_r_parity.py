"""Slow-mode R-parity decomposition for the F86 diagnostic question.

R = site-reflection on the spin chain = momentum-reversal k ↔ N+1-k on JW
Bogoliubov modes. R commutes with the full Liouvillian L for uniform-γ
Z-dephasing on the XY chain, so L block-diagonalizes by R-parity.

This script:
1. Builds L = L_H + L_D for N-site XY chain.
2. Block-diagonalizes by R-parity (R-even / R-odd subspaces).
3. Reports the slowest modes (smallest |Re(λ)|) per R-parity block.
4. Tabulates Re(λ) distribution per R-parity.

The diagnostic question: does the F86 slow-mode landscape (Q_peak, EP) live
entirely in one R-parity, or do both contain slow modes? Channel-uniform
vectors are R-even by construction, so F86's L_eff lives in R-even; if R-odd
also contains slow modes, they are a separate spectral feature not visible
to the channel-uniform probe.

Run: python simulations/_slow_modes_r_parity.py [N] [K]
  N = chain length (default 4; practical: 4, 5, 6)
  K = how many slowest modes to display per parity (default 16)
"""
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
K = int(sys.argv[2]) if len(sys.argv) > 2 else 16
J = 1.0
GAMMA = 0.05

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


def build_liouvillian(N, gamma, J):
    dim = 2 ** N
    H = np.zeros((dim, dim), dtype=complex)
    for l in range(N - 1):
        for axis_idx in (1, 2):
            op_alpha = [0] * N
            op_alpha[l] = axis_idx
            op_alpha[l + 1] = axis_idx
            H += (J / 2) * pauli_string(op_alpha)
    I_dim = np.eye(dim, dtype=complex)
    L_H = -1j * (np.kron(I_dim, H) - np.kron(H.T, I_dim))
    L_D = np.zeros((dim * dim, dim * dim), dtype=complex)
    for l in range(N):
        op_alpha = [0] * N
        op_alpha[l] = 3
        sig_z_l = pauli_string(op_alpha)
        L_D += gamma * (np.kron(sig_z_l, sig_z_l) - np.eye(dim * dim, dtype=complex))
    return L_H + L_D


def reverse_bits(b, n):
    out = 0
    for i in range(n):
        if (b >> i) & 1:
            out |= 1 << (n - 1 - i)
    return out


print(f"# N = {N}, γ₀ = {GAMMA}, J = {J}, XY chain")
print(f"# Slow-mode R-parity decomposition")
print(f"#   L dim = {4**N} x {4**N}")

dim = 2 ** N
op_dim = dim * dim
print(f"# Building L...")
L = build_liouvillian(N, GAMMA, J)

R_h = np.zeros((dim, dim), dtype=float)
for b in range(dim):
    R_h[reverse_bits(b, N), b] = 1.0
R_op = np.kron(R_h, R_h)

print(f"# Computing R-parity eigenbasis (op_dim = {op_dim})...")
r_vals, r_vecs = np.linalg.eigh(R_op)
even_mask = np.abs(r_vals - 1) < 1e-6
odd_mask = np.abs(r_vals + 1) < 1e-6
even_basis = r_vecs[:, even_mask]
odd_basis = r_vecs[:, odd_mask]
print(f"#   R-even subspace dim: {even_basis.shape[1]}")
print(f"#   R-odd  subspace dim: {odd_basis.shape[1]}")

print(f"# Projecting and diagonalizing L within R-even...")
L_even_proj = even_basis.T @ L @ even_basis
eigvals_even = np.linalg.eigvals(L_even_proj)

print(f"# Projecting and diagonalizing L within R-odd...")
L_odd_proj = odd_basis.T @ L @ odd_basis
eigvals_odd = np.linalg.eigvals(L_odd_proj)

order_even = np.argsort(np.abs(eigvals_even.real))
order_odd = np.argsort(np.abs(eigvals_odd.real))

print(f"\n# Slowest {K} R-even modes (smallest |Re(λ)|):")
for i in range(min(K, len(eigvals_even))):
    e = eigvals_even[order_even[i]]
    print(f"  #{i:>3}: Re = {e.real:+.6f}  Im = {e.imag:+.6f}")

print(f"\n# Slowest {K} R-odd modes (smallest |Re(λ)|):")
for i in range(min(K, len(eigvals_odd))):
    e = eigvals_odd[order_odd[i]]
    print(f"  #{i:>3}: Re = {e.real:+.6f}  Im = {e.imag:+.6f}")

stat_thresh = 1e-6
n_stat_even = int(np.sum(np.abs(eigvals_even.real) < stat_thresh))
n_stat_odd = int(np.sum(np.abs(eigvals_odd.real) < stat_thresh))
print(f"\n# Stationary subspace (|Re(λ)| < {stat_thresh}):")
print(f"  R-even: {n_stat_even} modes")
print(f"  R-odd:  {n_stat_odd} modes")

print(f"\n# Re(λ) distribution per R-parity:")
re_even = np.round(eigvals_even.real, 4)
re_odd = np.round(eigvals_odd.real, 4)
all_re = sorted(set(re_even.tolist()) | set(re_odd.tolist()), reverse=True)
print(f"#   Re(λ)         R-even   R-odd")
for r in all_re:
    n_e = int(np.sum(re_even == r))
    n_o = int(np.sum(re_odd == r))
    if n_e + n_o > 0:
        marker = ""
        if abs(r) < 1e-6:
            marker = "   <-- stationary"
        elif abs(r - (-N * GAMMA)) < 1e-6:
            marker = "   <-- axis (n_XY = N/2)"
        print(f"  Re = {r:+.4f}    {n_e:>5}   {n_o:>5}{marker}")
