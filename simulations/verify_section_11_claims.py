"""Quick verification of the three numerical claims in COCKPIT_SCALING.md Section 11:
  1. exactly 50 Liouvillian modes at Re(lambda) = -0.2
  2. these 50 carry ~92% of Bell pair n_XY=2 sector projection
  3. 22 of them carry 90%
Reuses the setup from path_d_bell_pair_absorption.py.
"""
import numpy as np
from itertools import product

N = 5
J = 1.0
gamma = 0.05
BELL_QUBITS = (2, 3)

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': Xm, 'Y': Ym, 'Z': Zm}
LABELS = ['I', 'X', 'Y', 'Z']

def kron_chain(ops):
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out

def pauli_string_matrix(s):
    return kron_chain([PAULIS[c] for c in s])

def n_xy(s):
    return sum(1 for c in s if c in ('X', 'Y'))

def build_H(N, J):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in (Xm, Ym, Zm):
            ops = [I2] * N
            ops[i] = P
            ops[i + 1] = P
            H += J * kron_chain(ops)
    return H

def build_L(N, H, gamma):
    d = 2**N
    Id = np.eye(d, dtype=complex)
    L_H = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    L_D = np.zeros((d * d, d * d), dtype=complex)
    for k in range(N):
        ops = [I2] * N
        ops[k] = Zm
        Zk = kron_chain(ops)
        L_D += gamma * (np.kron(Zk, Zk.T) - np.kron(Id, Id))
    return L_H + L_D

def build_pauli_basis(N):
    d = 2**N
    strings = [''.join(p) for p in product(LABELS, repeat=N)]
    basis = np.zeros((d * d, len(strings)), dtype=complex)
    for j, s in enumerate(strings):
        M = pauli_string_matrix(s)
        basis[:, j] = M.flatten(order='F') / np.sqrt(d)
    return strings, basis

def mode_n_xy_value(coeffs, strings):
    p = np.abs(coeffs) ** 2
    p_sum = p.sum()
    if p_sum < 1e-15:
        return 0.0
    weights = p / p_sum
    return sum(w * n_xy(s) for w, s in zip(weights, strings))

# ----- Setup -----
d = 2**N
H = build_H(N, J)
L = build_L(N, H, gamma)
print(f"L shape: {L.shape}")

eigvals, R = np.linalg.eig(L)
R_inv = np.linalg.inv(R)
strings, P_basis = build_pauli_basis(N)
print(f"Eigenvalues computed: {len(eigvals)}")

# n_XY per mode
mode_n = np.zeros(len(eigvals))
for k in range(len(eigvals)):
    coeffs = P_basis.conj().T @ R[:, k]
    mode_n[k] = mode_n_xy_value(coeffs, strings)

# Initial state
def build_initial_state(N, bell_qubits):
    Phi_plus = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    rho_bell = np.outer(Phi_plus, Phi_plus.conj())
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    rho_plus = np.outer(plus, plus.conj())
    b0, b1 = bell_qubits
    factors = []
    i = 0
    while i < N:
        if i == b0:
            factors.append(rho_bell)
            i += 2
        else:
            factors.append(rho_plus)
            i += 1
    return kron_chain(factors)

rho0 = build_initial_state(N, BELL_QUBITS)
vec_rho0 = rho0.flatten(order='F')
c = R_inv @ vec_rho0

# ----- CLAIMS TO VERIFY -----
print()
print("=" * 60)
print("VERIFICATION OF SECTION 11 CLAIMS")
print("=" * 60)

# CLAIM 1: exactly 50 modes at Re(lambda) = -0.2
# Be precise about "at Re=-0.2" — use floating point tolerance
tol = 1e-9
re_eigs = eigvals.real
mask_at_minus_0p2 = np.abs(re_eigs - (-0.2)) < tol
n_at_plateau = mask_at_minus_0p2.sum()
print(f"\nClaim 1: exactly 50 modes at Re(lambda) = -0.2")
print(f"  Found (tolerance {tol}): {n_at_plateau}")
print(f"  Claim correct? {n_at_plateau == 50}")

# Sanity: are these all also n_XY = 2 modes?
n_xy_at_plateau = mode_n[mask_at_minus_0p2]
print(f"  All n_XY values for these modes: min={n_xy_at_plateau.min():.6f}, "
      f"max={n_xy_at_plateau.max():.6f}, mean={n_xy_at_plateau.mean():.6f}")

# CLAIM 2: these 50 modes carry ~92% of Bell pair's projection onto n_XY=2 sector
# First: what IS the Bell pair's projection onto the n_XY=2 sector?
overlap = np.abs(c) ** 2
total_overlap = overlap.sum()  # should be approximately something well-defined

# n_XY=2 sector means: modes with n_xy ~ 2 (use a wider tolerance for plateau identification)
n_xy_2_mask = (mode_n > 1.5) & (mode_n < 2.5)
overlap_in_n_xy_2_sector = overlap[n_xy_2_mask].sum()
print(f"\nClaim 2: 50 plateau modes carry ~92% of Bell pair n_XY=2 projection")
print(f"  Total |c|^2 across all modes: {total_overlap:.6f}")
print(f"  |c|^2 in n_XY=2 sector (1.5 < n_XY < 2.5): {overlap_in_n_xy_2_sector:.6f}")
print(f"  |c|^2 in 50 plateau modes (Re=-0.2 exactly): {overlap[mask_at_minus_0p2].sum():.6f}")
fraction_of_sector = overlap[mask_at_minus_0p2].sum() / overlap_in_n_xy_2_sector
print(f"  Fraction = {fraction_of_sector*100:.2f}%")
print(f"  Claim of ~92% correct? {abs(fraction_of_sector - 0.92) < 0.05}")

# CLAIM 3: 22 of them carry 90%
plateau_overlaps = overlap[mask_at_minus_0p2]
plateau_overlaps_sorted = np.sort(plateau_overlaps)[::-1]  # descending
cumulative = np.cumsum(plateau_overlaps_sorted) / plateau_overlaps_sorted.sum()
n_for_90pct = np.searchsorted(cumulative, 0.90) + 1
print(f"\nClaim 3: 22 of 50 modes already carry 90%")
print(f"  Sorted plateau overlaps (top 5): {plateau_overlaps_sorted[:5]}")
print(f"  Number needed for 90% of plateau-internal mass: {n_for_90pct}")
print(f"  Claim of 22 correct? {n_for_90pct == 22}")

# Also: 90% of WHAT? Let me also compute 90% of the total n_XY=2 sector mass:
cumulative_global = np.cumsum(plateau_overlaps_sorted) / overlap_in_n_xy_2_sector
n_for_90pct_global = np.searchsorted(cumulative_global, 0.90) + 1
print(f"  Number needed for 90% of n_XY=2 sector mass: {n_for_90pct_global}")

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  Claim 1 (50 modes at Re=-0.2):     {'PASS' if n_at_plateau == 50 else f'FAIL ({n_at_plateau})'}")
pct_str = f"{fraction_of_sector*100:.1f}%"
print(f"  Claim 2 (~92% of n_XY=2 sector):  {pct_str}")
print(f"  Claim 3a (22 for 90% plateau-internal):  {n_for_90pct}")
print(f"  Claim 3b (22 for 90% sector-global):     {n_for_90pct_global}")
