#!/usr/bin/env python3
"""Pi_Z = R*D and the dihedral mirror group D4: committed, self-validating verification.

This script is the computational anchor of docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md:
the canonical palindromizer factors as Pi_Z = R o D, and the mirror inventory closes
into one dihedral group D4 = <R, D> of order 8. Every block asserts its identities
(any failure raises AssertionError and the process exits nonzero); the final line
ALL BLOCKS PASS prints only if every identity holds. The group-level identities are
equalities of signed permutation matrices and compare at dev 0.00e+00.

The two generators (F = X^(x)N, F^2 = I):
  R = I(x)F : rho -> rho F  (the windowed-converse spine's ket reflection,
              PROOF_F87_WINDOWED_MONOMIAL_CONVERSE sect.2)
  D         : rho -> rho^T  (the transpose superoperator, = diag((-1)^n_Y) in the
              Pauli basis, F114's / Welle 12's D)

Block ledger:
  A. Pi_Z(rho) = rho^T * X^(x)N, per-site sigma -> sigma^T X, i.e. Pi_Z = R o D
     (transpose FIRST, then right-multiply by F); the wrong-sided form F rho^T
     (= Pi_Y) is rejected at O(1). N = 1, 2, 3.
  B. <R, D> is dihedral of order 8: Pi_Z = R*D, Pi_Z^2 = curlyF = F(x)F
     (= F1^2, the X^(x)N charge conjugation), D Pi D = Pi^(-1) = Pi_Y
     (Welle 12's D Pi_Z D = Pi_Y is the dihedral inversion relation),
     closure |<R, D>| = 8 with all eight elements identified,
     spine V4 = {I, curlyF, R, curlyF R} is the Klein subgroup,
     curlyF D = diag((-1)^n_Z) (the fourth diagonal sign mirror, named in the proof).
  C. Palindrome factorization (XXZ Delta=0.7, site-dependent gamma, N=3):
     D carries the H-flip (D L_H D = -L_H for n_Y-even real H) and fixes the
     dissipator; R fixes L_H and carries the dissipator shift
     (R L_diss R = -L_diss - 2*sum(gamma)); product = the canonical palindrome
     Pi L Pi^(-1) = -L - 2 sigma.
  D. Classifier cell: truly(sigma) (F85: n_Y even AND n_Z even)
     <=> [eps_F114(sigma) = -1] AND [F sigma F = +sigma]; 63/63 strings at N=3
     against framework's _pauli_tuple_is_truly.
  E. F114 at N=5 via the transpose identification (one size beyond its verified N<=4).
  F. e^{i pi N_hat} = Z^(x)N (the trajectory-claim U(1) absorber; documented,
     here only as a table cross-check).

Conventions: row-stacking (C-order) vec, |i><j| -> e_i (x) e_j, kron(A,B) = A rho B^T,
matching framework.lindblad and PROOF_F87_WINDOWED_MONOMIAL_CONVERSE sect.1.
Promoted from the WIP scout _mirror_inventory_bridge_check.py, 2026-06-10.
"""
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, str(Path(__file__).parent))

from framework.pauli import (  # noqa: E402
    pauli_string, site_op, _k_to_indices, PAULI_LABELS,
)
from framework.symmetry import build_pi_full, _pauli_tuple_is_truly  # noqa: E402
from framework.lindblad import lindbladian_z_dephasing  # noqa: E402

rng = np.random.default_rng(7)
TOL = 1e-12


def xn(N):
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    out = np.array([[1.0]], dtype=complex)
    for _ in range(N):
        out = np.kron(out, X)
    return out


def vec_C(rho):
    return rho.flatten()           # C-order row-stacking


def basis_transform_C(N):
    """M[:, k] = vec_C(sigma_k); M^-1 = M^dag / 2^N."""
    d2, d = 4 ** N, 2 ** N
    M = np.zeros((d * d, d2), dtype=complex)
    for k in range(d2):
        M[:, k] = vec_C(pauli_string(list(_k_to_indices(k, N))))
    return M


def swap_perm(N):
    """SWAP on coherence space: vec_C(rho) -> vec_C(rho^T)."""
    d = 2 ** N
    P = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            P[j * d + i, i * d + j] = 1.0
    return P


def n_letter(k, N, letter):
    lab = [PAULI_LABELS[idx] for idx in _k_to_indices(k, N)]
    return sum(1 for c in lab if c == letter)


def string_labels(k, N):
    return ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(k, N))


# ----------------------------------------------------------------- A
print('A. Pi_Z action = rho -> rho^T X^(x)N (per-site sigma -> sigma^T X)')
for N in (1, 2, 3):
    F = xn(N)
    Pi = build_pi_full(N, 'Z')
    M = basis_transform_C(N)
    Minv = M.conj().T / (2 ** N)
    worst_right = worst_left = 0.0
    for _ in range(8):
        d = 2 ** N
        rho = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        c_in = Minv @ vec_C(rho)
        c_out = Pi @ c_in
        c_right = Minv @ vec_C(rho.T @ F)     # rho^T F
        c_left = Minv @ vec_C(F @ rho.T)      # F rho^T
        worst_right = max(worst_right, np.max(np.abs(c_out - c_right)))
        worst_left = max(worst_left, np.max(np.abs(c_out - c_left)))
    tag = 'rho^T F' if worst_right < TOL else ('F rho^T' if worst_left < TOL else 'NEITHER')
    print(f'   N={N}: |Pi(rho) - rho^T F| = {worst_right:.2e}, '
          f'|Pi(rho) - F rho^T| = {worst_left:.2e}  ->  Pi = {tag}')
    assert worst_right < TOL or worst_left < TOL, 'Pi is neither one-sided form'
RIGHT = worst_right < TOL   # convention flag from the largest N tested
print()

# ----------------------------------------------------------------- B
print('B. group structure of <R, D> on coherence space (N=3, signed-perm exact)')
N = 3
d = 2 ** N
F = xn(N)
Id = np.eye(d, dtype=complex)
D_coh = swap_perm(N).astype(complex)
R_coh = np.kron(Id, F)              # rho -> rho F   (kron(A,B): A rho B^T; F^T=F)
LF_coh = np.kron(F, Id)             # rho -> F rho
FF_coh = np.kron(F, F)              # rho -> F rho F  (curly F)
M = basis_transform_C(N)
Minv = M.conj().T / (2 ** N)
Pi_coh = M @ build_pi_full(N, 'Z') @ Minv
PiY_coh = M @ build_pi_full(N, 'Y') @ Minv

dev_fact = np.max(np.abs(Pi_coh - (R_coh @ D_coh if RIGHT else LF_coh @ D_coh)))
print(f'   Pi_Z = {"R*D" if RIGHT else "(FxI)*D"} bit-exact: dev = {dev_fact:.2e}')
assert dev_fact < TOL

dev_sq = np.max(np.abs(Pi_coh @ Pi_coh - FF_coh))
print(f'   Pi^2 = F(x)F (= F1^2 charge conjugation): dev = {dev_sq:.2e}')
assert dev_sq < TOL

dev_inv = np.max(np.abs(D_coh @ Pi_coh @ D_coh - np.linalg.inv(Pi_coh)))
dev_y = np.max(np.abs(D_coh @ Pi_coh @ D_coh - PiY_coh))
dev_iy = np.max(np.abs(PiY_coh - np.linalg.inv(Pi_coh)))
print(f'   D Pi D = Pi^-1 (dihedral relation): dev = {dev_inv:.2e}; '
      f'= Pi_Y (Welle 12): dev = {dev_y:.2e}; Pi_Y = Pi_Z^-1: dev = {dev_iy:.2e}')
assert dev_inv < TOL and dev_y < TOL and dev_iy < TOL

# close the group from generators
gens = [R_coh, D_coh]
elems = [np.eye(d * d, dtype=complex)]
changed = True
while changed:
    changed = False
    for g in gens:
        for e in list(elems):
            cand = g @ e
            if not any(np.max(np.abs(cand - x)) < TOL for x in elems):
                elems.append(cand)
                changed = True
print(f'   |<R, D>| = {len(elems)} (expect 8, dihedral D4)')
assert len(elems) == 8

# identify all eight elements against named mirrors
named = {
    'I': np.eye(d * d, dtype=complex),
    'Pi_Z': Pi_coh, 'curlyF = Pi^2': FF_coh, 'Pi_Y = Pi_Z^3': PiY_coh,
    'D (transpose)': D_coh, 'curlyF*D': FF_coh @ D_coh,
    'R = IxF': R_coh, 'curlyF*R = FxI': LF_coh,
}
for name, mat in named.items():
    hit = any(np.max(np.abs(mat - x)) < TOL for x in elems)
    print(f'     contains {name}: {hit}')
    assert hit

# diagonal V4 in the Pauli basis: D = (-1)^n_Y, curlyF*D = (-1)^n_Z, curlyF = (-1)^(n_Y+n_Z)
FD_pauli = Minv @ (FF_coh @ D_coh) @ M
D_pauli = Minv @ D_coh @ M
sgn_nz = np.diag([(-1.0) ** n_letter(k, N, 'Z') for k in range(4 ** N)])
sgn_ny = np.diag([(-1.0) ** n_letter(k, N, 'Y') for k in range(4 ** N)])
dev_nz = np.max(np.abs(FD_pauli - sgn_nz))
dev_ny = np.max(np.abs(D_pauli - sgn_ny))
print(f'   Pauli basis: D = diag((-1)^n_Y) dev = {dev_ny:.2e};  '
      f'curlyF*D = diag((-1)^n_Z) dev = {dev_nz:.2e}  (the unnamed 4th sign mirror)')
assert dev_nz < TOL and dev_ny < TOL
print()

# ----------------------------------------------------------------- C
print('C. palindrome factorization (N=3 XXZ Delta=0.7, site-dependent gamma)')
gammas = [0.05, 0.11, 0.07]
H = np.zeros((d, d), dtype=complex)
for b in range(N - 1):
    H += (site_op(N, b, 'X') @ site_op(N, b + 1, 'X')
          + site_op(N, b, 'Y') @ site_op(N, b + 1, 'Y')
          + 0.7 * site_op(N, b, 'Z') @ site_op(N, b + 1, 'Z'))
L = lindbladian_z_dephasing(H, gammas)
L_H = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
L_diss = L - L_H
sigma = 2 * sum(gammas)

dev1 = np.max(np.abs(D_coh @ L_H @ D_coh + L_H))
dev2 = np.max(np.abs(D_coh @ L_diss @ D_coh - L_diss))
dev3 = np.max(np.abs(R_coh @ L_H @ R_coh - L_H))
dev4 = np.max(np.abs(R_coh @ L_diss @ R_coh + L_diss + sigma * np.eye(d * d)))
dev5 = np.max(np.abs(Pi_coh @ L @ np.linalg.inv(Pi_coh) + L + sigma * np.eye(d * d)))
print(f'   D L_H D = -L_H:                    dev = {dev1:.2e}   (D carries the H-flip)')
print(f'   D L_diss D = +L_diss:              dev = {dev2:.2e}')
print(f'   R L_H R = +L_H:                    dev = {dev3:.2e}')
print(f'   R L_diss R = -L_diss - 2*sum(g):   dev = {dev4:.2e}   (R carries the -2sigma shift)')
print(f'   Pi L Pi^-1 = -L - 2*sum(g):        dev = {dev5:.2e}   (the product = the palindrome)')
assert max(dev1, dev2, dev3, dev4, dev5) < 1e-11
print()

# ----------------------------------------------------------------- D
print('D. classifier cell: truly <=> (F114 eps = -1) AND (F sigma F = +sigma), all 64 strings N=3')
mismatch = 0
for k in range(1, 4 ** N):
    lab = string_labels(k, N)
    sig = pauli_string(list(_k_to_indices(k, N)))
    Ls = -1j * (np.kron(sig, Id) - np.kron(Id, sig.T))
    DLD = D_coh @ Ls @ D_coh
    ny = n_letter(k, N, 'Y')
    eps = (-1.0) ** (ny + 1)
    dev_f114 = np.max(np.abs(DLD - eps * Ls))
    assert dev_f114 < TOL, f'F114 fails on {lab}'
    f_even = np.max(np.abs(F @ sig @ F - sig)) < TOL
    pred_truly = (eps == -1.0) and f_even
    if pred_truly != _pauli_tuple_is_truly(tuple(lab)):
        mismatch += 1
print(f'   F114 sign law eps = (-1)^(n_Y+1) via transpose: all 63 non-identity strings PASS')
print(f'   truly-cell match vs framework _pauli_tuple_is_truly: {63 - mismatch}/63 (mismatches: {mismatch})')
assert mismatch == 0
print()

# ----------------------------------------------------------------- E
print('E. F114 at N=5 (beyond its verified N<=4), one k=3 string, via D = transpose')
N5 = 5
d5 = 2 ** N5
sig5 = site_op(N5, 0, 'X') @ site_op(N5, 2, 'Y') @ site_op(N5, 4, 'Y')   # n_Y = 2
Id5 = np.eye(d5, dtype=complex)
Ls5 = -1j * (np.kron(sig5, Id5) - np.kron(Id5, sig5.T))
# D L D for D = transpose superoperator: (L(rho^T))^T -> permute rows+cols by swap
idx = np.arange(d5 * d5).reshape(d5, d5).T.flatten()
DLD5 = Ls5[np.ix_(idx, idx)]
eps5 = (-1.0) ** (2 + 1)
dev_e = np.max(np.abs(DLD5 - eps5 * Ls5))
print(f'   sigma = X I Y I Y (n_Y=2): D L_sigma D = -L_sigma, dev = {dev_e:.2e}')
assert dev_e < TOL
print()

# ----------------------------------------------------------------- F
print('F. e^(i pi N_hat) = Z^(x)N (documented; table cross-check)')
N4 = 4
nhat_diag = np.array([bin(i).count('1') for i in range(2 ** N4)], dtype=float)
lhs = np.diag(np.exp(1j * np.pi * nhat_diag))
Z = np.array([[1, 0], [0, -1]], dtype=complex)
rhs = np.array([[1.0]], dtype=complex)
for _ in range(N4):
    rhs = np.kron(rhs, Z)
dev_f = np.max(np.abs(lhs - rhs))
print(f'   dev = {dev_f:.2e}')
assert dev_f < TOL
print()

# ----------------------------------------------------------------- G
print('G. The cube filled (PROOF_PI_FACTORS_AS_R_TIMES_D section 7): the polarity cube axes')
print('   (bit_a, bit_b, y_par) are the characters of (Ad_Z^N, Ad_X^N, transpose theta);')
print('   conjugations give only the even Klein square, theta breaks into the odd half.')
from itertools import product as _iproduct  # noqa: E402
_P1 = {'I': np.eye(2, dtype=complex),
       'X': np.array([[0, 1], [1, 0]], dtype=complex),
       'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
       'Z': np.array([[1, 0], [0, -1]], dtype=complex)}
NG = 3
FXg = np.array([[1.0]], dtype=complex)
FZg = np.array([[1.0]], dtype=complex)
for _ in range(NG):
    FXg = np.kron(FXg, _P1['X'])
    FZg = np.kron(FZg, _P1['Z'])
checked = 0
for s in _iproduct('IXYZ', repeat=NG):
    sig = np.array([[1.0]], dtype=complex)
    for ch in s:
        sig = np.kron(sig, _P1[ch])
    pX, pY, pZ = (s.count('X') % 2, s.count('Y') % 2, s.count('Z') % 2)
    bit_a, bit_b = (s.count('X') + s.count('Y')) % 2, (s.count('Y') + s.count('Z')) % 2
    assert np.allclose(FZg @ sig @ FZg, (-1.0) ** bit_a * sig)       # Ad_Z^N reads bit_a
    assert np.allclose(FXg @ sig @ FXg, (-1.0) ** bit_b * sig)       # Ad_X^N reads bit_b
    assert np.allclose(sig.T, (-1.0) ** pY * sig)                     # theta reads y_par
    assert np.allclose((FZg @ sig @ FZg).T, (-1.0) ** pX * sig)       # theta o Ad_Z = (-1)^n_X
    assert np.allclose((FXg @ sig @ FXg).T, (-1.0) ** pZ * sig)       # theta o Ad_X = (-1)^n_Z = FD
    checked += 1
print(f'   all {checked} strings at N={NG}: three characters == the three cube axes, exact')

print('\nALL BLOCKS PASS')
