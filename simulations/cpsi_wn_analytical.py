"""
CΨ(0) for W_N: analytical formula vs numerical verification.
Formula: CΨ(0) = 2(N^2 - 4N + 8) / (3N^3)

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 10, 2026
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy import linalg


def cpsi_wn_formula(N):
    """Closed-form CΨ(0) for W_N on any adjacent pair."""
    return 2.0 * (N**2 - 4*N + 8) / (3.0 * N**3)


def cpsi_wn_numerical(N):
    """Numerical CΨ(0) for W_N via explicit density matrix construction."""
    d = 2**N
    # W_N state
    psi = np.zeros(d, dtype=complex)
    for k in range(N):
        psi[1 << k] = 1.0 / np.sqrt(N)
    rho = np.outer(psi, psi.conj())

    # Partial trace to pair (0, 1) -- representative adjacent pair
    # Reshape to N-qubit tensor, trace out qubits 2..N-1
    rho_t = rho.reshape([2]*N + [2]*N)
    keep = [0, 1]
    in_idx = list(range(2*N))
    out_idx = [0, 1, N, N+1]
    for k in range(2, N):
        in_idx[k + N] = in_idx[k]  # trace over qubit k
    rho_01 = np.einsum(rho_t, in_idx, out_idx).reshape(4, 4)

    # CΨ = Tr(rho^2) * L1 / 3
    purity = float(np.real(np.trace(rho_01 @ rho_01)))
    l1 = float(np.sum(np.abs(rho_01)) - np.sum(np.abs(np.diag(rho_01))))
    return purity * l1 / 3.0


print("CΨ(0) for W_N: analytical vs numerical")
print(f"Formula: CΨ(0) = 2(N^2 - 4N + 8) / (3N^3)")
print()
print(f"{'N':>3}  {'analytical':>12}  {'numerical':>12}  {'diff':>12}  {'< 1/4?':>6}")
print(f"{'-' * 52}")

all_pass = True
for N in range(2, 11):
    ana = cpsi_wn_formula(N)
    num = cpsi_wn_numerical(N)
    diff = abs(ana - num)
    below = "YES" if ana < 0.25 else "no"
    print(f"{N:>3}  {ana:>12.8f}  {num:>12.8f}  {diff:>12.2e}  {below:>6}")
    if diff > 1e-10:
        print(f"  ** MISMATCH at N={N} **")
        all_pass = False

print()
if all_pass:
    print("All N=2-10 match within 1e-10. Formula verified.")
else:
    print("VERIFICATION FAILED.")

print()
print("Corollary: CΨ(0) < 1/4 for all N >= 3.")
print(f"  N=2: CΨ = {cpsi_wn_formula(2):.6f} = 1/3 > 1/4 (Bell+ crosses)")
print(f"  N=3: CΨ = {cpsi_wn_formula(3):.6f} = 10/81 < 1/4 (below fold)")
print(f"  lim N->inf: CΨ -> 2/(3N) -> 0")
