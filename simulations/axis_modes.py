"""Axis modes at general even N: the n_XY = N/2 layer of the Liouvillian.

Parameterized version of `axis_modes_n4.py` (which is the frozen N=4 reference).
Same physics; N taken from CLI:

    python simulations/axis_modes.py [N]

Default N = 4. Practical range: N = 4, 6. N = 8 would need a 65536² L matrix
(out of dense-eigendecomp range on workstation memory). Odd N gives an empty
axis subspace at Re(λ) = −Nγ₀ (the integer-γ₀ absorption grid skips it); the
script reports this and proceeds with whatever modes do hit the layer.

Tests the parity finding from `project_minus_is_the_mirror` and applies the
Majorana lens (see experiments/MAJORANA_AXIS_MODES.md) at the chosen N:
- Site-reflection R sorts axis modes into R-even / R-odd
- Im(λ) clusters decompose into integer combinations of the single-particle
  Bloch dispersion ε(k) = 2J·cos(πk/(N+1))
- Silent modes (Im = 0) are operator-space Majorana self-conjugate σ_(a,b) with
  λ_a = λ_b in the many-body H spectrum

XY chain, J = 1, uniform Z-dephasing γ₀ = 0.05.
"""
import itertools
import sys
from math import comb

import numpy as np

# Windows console defaults to cp1252; force UTF-8 so γ₀, λ, σ, ⟨⟩ render.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
J = 1.0
GAMMA = 0.05

PAULI = np.array([
    [[1, 0], [0, 1]],     # I
    [[0, 1], [1, 0]],     # X
    [[0, -1j], [1j, 0]],  # Y
    [[1, 0], [0, -1]],    # Z
], dtype=complex)
PAULI_NAMES = "IXYZ"


def pauli_string(alpha):
    result = PAULI[alpha[0]]
    for a in alpha[1:]:
        result = np.kron(result, PAULI[a])
    return result


def n_xy(alpha):
    return sum(1 for a in alpha if a in (1, 2))


def build_liouvillian(N, gamma, J):
    """L = L_H + L_D for an N-site XY chain under uniform Z-dephasing.
    Vec convention: column-stacking, so vec(H ρ − ρ H) = (I⊗H − H^T⊗I) vec(ρ).
    """
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


print(f"# γ₀ = {GAMMA}, J = {J}, XY chain")
print(f"# N = {N} ({'even' if N % 2 == 0 else 'odd'})")
print(f"# Axis: Re(λ) = −Σγ = {-N * GAMMA:+.4f}")
print(f"# Predicted axis layer dim (n_XY = N/2): "
      f"{comb(N, N // 2) * 2**N if N % 2 == 0 else 0}")

dim = 2 ** N
op_dim = dim * dim
print(f"# L matrix dim: {op_dim} x {op_dim} "
      f"({op_dim * op_dim * 16 / 1e9:.2f} GB at complex128)")

L = build_liouvillian(N, GAMMA, J)
eigvals, eigvecs = np.linalg.eig(L)
axis = -N * GAMMA
re = eigvals.real
im = eigvals.imag
mask = np.abs(re - axis) < 1e-9
axis_indices = np.where(mask)[0]

print(f"\n# Re distribution across all {len(eigvals)} modes:")
re_rounded = np.round(re, 9)
for r, c in zip(*np.unique(re_rounded, return_counts=True)):
    flag = " <-- axis" if abs(r - axis) < 1e-9 else ""
    print(f"  Re = {r:+.4f}  count {c}{flag}")

print(f"\n# Im distribution among the {len(axis_indices)} axis modes:")
if len(axis_indices) == 0:
    print("  (empty; no axis layer at this N)")
    sys.exit(0)

axis_ims = np.round(im[mask], 6)
unique_im, im_counts = np.unique(axis_ims, return_counts=True)
print(f"  {len(unique_im)} unique Im values:")
for v, c in sorted(zip(unique_im, im_counts)):
    print(f"  Im = {v:+.4f}  count {c}")

# === Pauli-basis decomposition ===
all_alphas = list(itertools.product(range(4), repeat=N))
pauli_basis = np.zeros((op_dim, op_dim), dtype=complex)
for k, alpha in enumerate(all_alphas):
    pauli_basis[:, k] = pauli_string(alpha).flatten('F')
nxy_per_alpha = np.array([n_xy(a) for a in all_alphas])

# === Im = 0 axis modes ===
im0_mask = mask & (np.abs(im) < 1e-9)
im0_indices = np.where(im0_mask)[0]

print(f"\n# Im=0 axis modes (silent): {len(im0_indices)} of {len(axis_indices)} axis modes")

# === Site-reflection R-parity ===
def reverse_bits(b, n):
    out = 0
    for i in range(n):
        if (b >> i) & 1:
            out |= 1 << (n - 1 - i)
    return out


R_h = np.zeros((dim, dim), dtype=complex)
for b in range(dim):
    R_h[reverse_bits(b, N), b] = 1.0

# Sanity: σ_{rev(α)} == R_h σ_α R_h on a representative
test_alpha = tuple([1, 2, 3, 0] + [0] * (N - 4)) if N >= 4 else (1, 2, 3)[:N]
assert np.allclose(pauli_string(test_alpha[::-1]),
                   R_h @ pauli_string(test_alpha) @ R_h), \
    "site-reflection inconsistent"

R_op = np.kron(R_h, R_h)

# Theoretical R-parity split of the n_XY = N/2 layer
if N % 2 == 0:
    nxy_palindromic = 0
    nxy_pairs = 0
    seen = set()
    for k, alpha in enumerate(all_alphas):
        if nxy_per_alpha[k] != N // 2 or k in seen:
            continue
        rev_k = all_alphas.index(alpha[::-1])
        if rev_k == k:
            nxy_palindromic += 1
            seen.add(k)
        else:
            nxy_pairs += 1
            seen.add(k)
            seen.add(rev_k)
    even_dim_pred = nxy_palindromic + nxy_pairs
    odd_dim_pred = nxy_pairs
    print(f"\n# Theoretical R-parity of n_XY = {N//2} layer "
          f"({nxy_palindromic} palindromic + {nxy_pairs} non-palindromic pairs):")
    print(f"  R-even subspace dim: {even_dim_pred}")
    print(f"  R-odd  subspace dim: {odd_dim_pred}")
    print(f"  total: {even_dim_pred + odd_dim_pred} "
          f"(= C(N, N/2)·2^N = {comb(N, N // 2) * 2**N})")

# R on Im=0 kernel
if len(im0_indices) > 0:
    V_K = eigvecs[:, im0_indices]
    U_K, _, _ = np.linalg.svd(V_K, full_matrices=False)
    R_sub_sym = (U_K.conj().T @ R_op @ U_K + (U_K.conj().T @ R_op @ U_K).conj().T) / 2
    r_eigvals = np.linalg.eigvalsh(R_sub_sym)
    n_even_im0 = int(np.sum(np.abs(r_eigvals - 1) < 1e-6))
    n_odd_im0 = int(np.sum(np.abs(r_eigvals + 1) < 1e-6))
    print(f"\n# Site-reflection in {len(im0_indices)}-dim Im=0 kernel:")
    print(f"  R-even (+1): {n_even_im0}")
    print(f"  R-odd  (−1): {n_odd_im0}")
    print(f"  ambig (|λ_R| ≠ 1): {len(r_eigvals) - n_even_im0 - n_odd_im0}")

# Full axis-subspace R decomposition with Im(λ) per parity
V_axis = eigvecs[:, axis_indices]
U_axis, _, _ = np.linalg.svd(V_axis, full_matrices=False)
R_axis_sym = (U_axis.conj().T @ R_op @ U_axis + (U_axis.conj().T @ R_op @ U_axis).conj().T) / 2
r_eigvals_axis, r_eigvecs_axis = np.linalg.eigh(R_axis_sym)

n_even_axis = int(np.sum(np.abs(r_eigvals_axis - 1) < 1e-6))
n_odd_axis = int(np.sum(np.abs(r_eigvals_axis + 1) < 1e-6))
ambig = len(r_eigvals_axis) - n_even_axis - n_odd_axis

print(f"\n# {'=' * 64}")
print(f"# R-parity decomposition of the full {len(axis_indices)}-dim axis subspace")
print(f"# {'=' * 64}")
print(f"  R-even axis modes: {n_even_axis}")
print(f"  R-odd  axis modes: {n_odd_axis}")
print(f"  ambig: {ambig}")

if n_even_axis > 0:
    even_mask = np.abs(r_eigvals_axis - 1) < 1e-6
    U_even = U_axis @ r_eigvecs_axis[:, even_mask]
    im_even = np.linalg.eigvals(U_even.conj().T @ L @ U_even).imag
    print(f"\n# Im(λ) within R-even axis subspace ({n_even_axis} modes):")
    for v, c in zip(*np.unique(np.round(im_even, 6), return_counts=True)):
        flag = " <-- silent" if abs(v) < 1e-6 else ""
        print(f"  Im = {v:+.4f}  count {c}{flag}")

if n_odd_axis > 0:
    odd_mask = np.abs(r_eigvals_axis + 1) < 1e-6
    U_odd = U_axis @ r_eigvecs_axis[:, odd_mask]
    im_odd = np.linalg.eigvals(U_odd.conj().T @ L @ U_odd).imag
    print(f"\n# Im(λ) within R-odd axis subspace ({n_odd_axis} modes):")
    for v, c in zip(*np.unique(np.round(im_odd, 6), return_counts=True)):
        flag = " <-- silent" if abs(v) < 1e-6 else ""
        print(f"  Im = {v:+.4f}  count {c}{flag}")

# === Single-particle Bloch dispersion (for reference) ===
print(f"\n# Single-particle Bloch dispersion ε(k) = 2J·cos(πk/(N+1)) at N = {N}:")
eps = [2 * J * np.cos(np.pi * k / (N + 1)) for k in range(1, N + 1)]
for k, e in enumerate(eps, 1):
    print(f"  k = {k}: ε = {e:+.4f}")
print(f"  Distinct positive: {[f'{e:+.4f}' for e in eps if e > 1e-9]}")
print(f"  Expected Im(λ) clusters: integer combinations of these single-particle energies.")
