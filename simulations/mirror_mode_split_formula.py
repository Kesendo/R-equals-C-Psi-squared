"""
Mirror-Mode Split Formula
=========================

Probe 1 (commit f708923) measured mode counts per (Re-class, Pi^2-parity)
sector at N=2-5. The asymmetry between even and odd sectors (4:6 at N=2,
122:124 at N=4) appears for even N and disappears for odd N. This script
derives and verifies the closed-form formula.

Formula:
  conserved per sector:  even = floor(N/2) + 1,  odd = ceil(N/2)
  correlation per sector: same as conserved (palindrome symmetry)
  mirror per sector:     2^(2N-1) - 2 * (conserved per sector)

Mechanism: the conserved modes (Re=0) are exactly the (N+1) elementary
symmetric polynomials e_d(Z_1, ..., Z_N) for d = 0, 1, ..., N. Their w_YZ
parities are (d mod 2), giving the observed sector split. The asymmetry
for even N comes from e_N (product of all Z's) having even parity when N
is even, odd parity when N is odd.

Date: 2026-04-15
"""

import numpy as np
from itertools import combinations, product as iproduct
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
paulis = [I2, X, Y, Z]
w_yz = [0, 0, 1, 1]


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def liouvillian(H, jumps):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def chain_xy_hamiltonian(N, J=1.0):
    H = np.zeros((2**N, 2**N), dtype=complex)
    for k in range(N - 1):
        ops_x = [I2] * N
        ops_y = [I2] * N
        ops_x[k] = X; ops_x[k+1] = X
        ops_y[k] = Y; ops_y[k+1] = Y
        H += 0.5 * J * (kron(*ops_x) + kron(*ops_y))
    return H


def elementary_sym_poly_in_Z(N, d):
    """e_d(Z_1, ..., Z_N) = sum over subsets S of size d of prod_{k in S} Z_k."""
    M = np.zeros((2**N, 2**N), dtype=complex)
    for subset in combinations(range(N), d):
        ops = [I2] * N
        for k in subset:
            ops[k] = Z
        M += kron(*ops)
    return M


def count_modes_per_sector(N, gamma_B=0.1):
    """Count modes per (Re-class, w_YZ-parity) for XY chain, gamma on last site."""
    H = chain_xy_hamiltonian(N)
    ops = [I2] * N
    ops[N - 1] = Z
    jumps = [np.sqrt(gamma_B) * kron(*ops)]
    L = liouvillian(H, jumps)

    # Pauli-basis change-of-basis and Pi^2-parity diagonal
    P_to_comp = np.zeros((4**N, 4**N), dtype=complex)
    parity_diag = np.zeros(4**N, dtype=complex)
    for k, tup in enumerate(iproduct(range(4), repeat=N)):
        op = paulis[tup[0]]
        for j in tup[1:]:
            op = np.kron(op, paulis[j])
        P_to_comp[:, k] = op.reshape(-1, order='F')
        parity_diag[k] = (-1) ** sum(w_yz[t] for t in tup)
    L_pauli = np.linalg.inv(P_to_comp) @ L @ P_to_comp

    even_idx = [k for k in range(4**N) if parity_diag[k].real > 0]
    odd_idx = [k for k in range(4**N) if parity_diag[k].real < 0]

    def classify(eigs):
        cons = mirror = corr = 0
        for e in eigs:
            re = e.real
            if abs(re) < 1e-6: cons += 1
            elif abs(re + 2 * gamma_B) < 1e-6: corr += 1
            else: mirror += 1
        return cons, mirror, corr

    return {
        'even': classify(np.linalg.eigvals(L_pauli[np.ix_(even_idx, even_idx)])),
        'odd': classify(np.linalg.eigvals(L_pauli[np.ix_(odd_idx, odd_idx)])),
    }


def predicted_split(N):
    sector = 2 ** (2 * N - 1)
    cons_e = N // 2 + 1
    cons_o = (N + 1) // 2
    return {
        'even': (cons_e, sector - 2 * cons_e, cons_e),
        'odd': (cons_o, sector - 2 * cons_o, cons_o),
    }


# ============================================================
# Test 1: formula vs direct computation, N=2..5
# ============================================================
print("=" * 64)
print("Test 1: Mirror-mode split formula vs direct computation")
print("=" * 64)
print()
print(f"  {'N':>2}  {'sector':>6} | {'even (c+m+r)':>14} | {'odd (c+m+r)':>14} | match")
print("  " + "-" * 56)
all_match = True
for N in range(2, 6):
    measured = count_modes_per_sector(N)
    predicted = predicted_split(N)
    sector = 2 ** (2 * N - 1)
    me = measured['even']; mo = measured['odd']
    pe = predicted['even']; po = predicted['odd']
    e_str = f"{me[0]}+{me[1]}+{me[2]}"
    o_str = f"{mo[0]}+{mo[1]}+{mo[2]}"
    match = (me == pe and mo == po)
    if not match: all_match = False
    print(f"  {N:>2}  {sector:>6} | {e_str:>14} | {o_str:>14} | {'OK' if match else 'FAIL'}")

print()
print("All N=2..5 match." if all_match else "MISMATCH.")


# ============================================================
# Test 2: conserved modes are e_d(Z_1, ..., Z_N) for d=0..N
# ============================================================
print()
print("=" * 64)
print("Test 2: Conserved modes are elementary symmetric polynomials e_d(Z)")
print("=" * 64)
print()
print(f"  {'N':>2} | {'cons modes':>10} | {'all e_d in cons subspace?':>26} | {'parities':>15}")
print("  " + "-" * 60)
for N in range(2, 5):
    H = chain_xy_hamiltonian(N)
    ops = [I2] * N
    ops[N - 1] = Z
    jumps = [np.sqrt(0.1) * kron(*ops)]
    L = liouvillian(H, jumps)
    eigvals, eigvecs = np.linalg.eig(L)
    cons_modes = [eigvecs[:, k].reshape((2**N, 2**N), order='F')
                  for k in range(4**N)
                  if abs(eigvals[k].real) < 1e-6 and abs(eigvals[k].imag) < 1e-6]
    cons_basis = np.array([M.flatten(order='F') for M in cons_modes]).T
    Q, _ = np.linalg.qr(cons_basis)
    max_residual = 0.0
    for d in range(N + 1):
        e_d = elementary_sym_poly_in_Z(N, d).flatten(order='F')
        proj = Q @ (Q.conj().T @ e_d)
        residual = np.linalg.norm(e_d - proj) / max(np.linalg.norm(e_d), 1e-14)
        max_residual = max(max_residual, residual)
    parities = ", ".join(str(d % 2) for d in range(N + 1))
    print(f"  {N:>2} | {len(cons_modes):>10} | "
          f"{'yes' if max_residual < 1e-10 else 'no':>26} | {parities:>15}")


# ============================================================
# Predictions, N=2 to 7
# ============================================================
print()
print("=" * 64)
print("Predictions, N = 2 to 7")
print("=" * 64)
print()
print(f"  {'N':>2} | {'cons (e,o)':>11} | {'mirror (e,o)':>16} | {'sector':>8}")
print("  " + "-" * 50)
for N in range(2, 8):
    cons_e = N // 2 + 1
    cons_o = (N + 1) // 2
    sector = 2 ** (2 * N - 1)
    mirror_e = sector - 2 * cons_e
    mirror_o = sector - 2 * cons_o
    print(f"  {N:>2} | ({cons_e:>2}, {cons_o:>2})    | "
          f"({mirror_e:>5}, {mirror_o:>5})  | {sector:>8}")

print()
print("Asymmetry pattern: even N has 1 extra even-parity conserved mode")
print("(from e_N having even parity when N is even). Odd N is balanced.")
