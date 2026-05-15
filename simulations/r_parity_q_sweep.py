"""R-parity-resolved Q-sweep: does R-odd have its own Q_peak / EP?

For each Q = J/γ₀ in a sweep, block-diagonalize the full Liouvillian L by
R-parity and track the two slowest non-stationary modes in each block.
EP coalescence appears as the Re-gap |Re(λ_1) − Re(λ_2)| reaching a local
minimum while the Im-gap opens (complex-conjugate pair regime).

Standard F86 Endpoint Q_peak at small N ≈ 1.15 (= 2/g_eff with g_eff ≈ 1.74).
The R-even block contains F86's channel-uniform L_eff so its EP should
correspond to the standard F86 prediction. The R-odd block is independent;
if its slowest pair coalesces at a different Q, that is a new diagnostic.

Run: python simulations/_r_parity_q_sweep.py [N]
  N = chain length (default 4; practical 4, 5)
"""
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
gamma = 0.05

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


def build_L(N, gamma, J):
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


dim = 2 ** N
R_h = np.zeros((dim, dim), dtype=float)
for b in range(dim):
    R_h[reverse_bits(b, N), b] = 1.0
R_op = np.kron(R_h, R_h)
r_vals, r_vecs = np.linalg.eigh(R_op)
even_basis = r_vecs[:, np.abs(r_vals - 1) < 1e-6]
odd_basis = r_vecs[:, np.abs(r_vals + 1) < 1e-6]

print(f"# N = {N}, γ₀ = {gamma}")
print(f"# R-even dim: {even_basis.shape[1]}, R-odd dim: {odd_basis.shape[1]}")
print(f"# Q-sweep: tracking slowest 2 non-stationary modes per R-block")
print()
print(f"  Q    | R-even #1            | R-even #2            | |ΔRe|_even |"
      f" R-odd #1             | R-odd #2             | |ΔRe|_odd")
print(f"  -----|----------------------|----------------------|------------|"
      f"----------------------|----------------------|----------")

Q_grid = np.linspace(0.3, 3.0, 28)
for Q in Q_grid:
    J = Q * gamma
    L = build_L(N, gamma, J)
    L_even = even_basis.T @ L @ even_basis
    L_odd = odd_basis.T @ L @ odd_basis
    eig_even = np.linalg.eigvals(L_even)
    eig_odd = np.linalg.eigvals(L_odd)
    eig_even_nonst = eig_even[np.abs(eig_even.real) > 1e-9]
    eig_odd_nonst = eig_odd[np.abs(eig_odd.real) > 1e-9]
    slow_even = eig_even_nonst[np.argsort(np.abs(eig_even_nonst.real))]
    slow_odd = eig_odd_nonst[np.argsort(np.abs(eig_odd_nonst.real))]
    e1, e2 = slow_even[0], slow_even[1]
    o1, o2 = slow_odd[0], slow_odd[1]
    re_gap_e = abs(e1.real - e2.real)
    re_gap_o = abs(o1.real - o2.real)
    print(f"  {Q:.2f} | ({e1.real:+.5f},{e1.imag:+.5f}) | "
          f"({e2.real:+.5f},{e2.imag:+.5f}) | {re_gap_e:.6f}   | "
          f"({o1.real:+.5f},{o1.imag:+.5f}) | "
          f"({o2.real:+.5f},{o2.imag:+.5f}) | {re_gap_o:.6f}")
