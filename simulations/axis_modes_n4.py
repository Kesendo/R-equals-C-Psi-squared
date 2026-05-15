"""Axis modes at small N: the n_XY = N/2 layer of the Liouvillian.

Tests the parity finding from `project_minus_is_the_mirror`:
- For even N, Re(λ) = −Σγ = −Nγ₀ has a self-mirror subspace of `C(N, N/2)·2^N` modes
  (the n_XY = N/2 absorption-grid layer).
- For odd N there is no such layer (the axis falls between integer-γ₀ steps).

At N=4 we expect 96 axis modes of 256 total; we then look at their Im distribution
and Pauli-basis decomposition to see whether L_H preserves the n_XY=2 layer or
mixes it with n_XY=0,4 averaging to 2.

XY chain, J = 1, uniform Z-dephasing γ₀ = 0.05.

Run: python simulations/axis_modes_n4.py

Frozen N=4 reference run. For parameterized exploration at other even N, see
`axis_modes.py` (same physics, N from CLI). See also experiments/MAJORANA_AXIS_MODES.md
for the documented findings (R-parity sorting + Majorana operator-space lens).
"""
import itertools
import sys
from math import comb

import numpy as np

# Windows console defaults to cp1252; force UTF-8 so γ₀, λ, σ, ⟨⟩ render.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

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
        for axis_idx in (1, 2):  # X, Y
            op_alpha = [0] * N
            op_alpha[l] = axis_idx
            op_alpha[l + 1] = axis_idx
            H += (J / 2) * pauli_string(op_alpha)
    I_dim = np.eye(dim, dtype=complex)
    L_H = -1j * (np.kron(I_dim, H) - np.kron(H.T, I_dim))
    L_D = np.zeros((dim * dim, dim * dim), dtype=complex)
    for l in range(N):
        op_alpha = [0] * N
        op_alpha[l] = 3  # Z
        sig_z_l = pauli_string(op_alpha)
        L_D += gamma * (np.kron(sig_z_l, sig_z_l) - np.eye(dim * dim, dtype=complex))
    return L_H + L_D


# === Parity check: N = 3 (odd, predicted 0) .. N = 5 (odd, predicted 0) ===
print(f"# γ₀ = {GAMMA}, J = {J}, XY chain")
print(f"# Parity check: axis modes at Re(λ) = −Σγ = −Nγ₀")
print()
for N in (3, 4, 5):
    L = build_liouvillian(N, GAMMA, J)
    eigvals = np.linalg.eigvals(L)
    axis = -N * GAMMA
    mask = np.abs(eigvals.real - axis) < 1e-9
    count = int(mask.sum())
    predicted = comb(N, N // 2) * (2 ** N) if N % 2 == 0 else 0
    parity = "even" if N % 2 == 0 else "odd"
    total = 4 ** N
    print(f"  N={N} ({parity:>4}): axis Re = {axis:+.4f}, "
          f"found {count:>3} / {total} at axis, predicted {predicted}")


# === Detailed analysis at N = 4 ===
print()
print("# " + "=" * 64)
print("# Detailed analysis at N = 4")
print("# " + "=" * 64)

N = 4
dim = 2 ** N
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
axis_ims = np.round(im[mask], 6)
unique_im, im_counts = np.unique(axis_ims, return_counts=True)
print(f"  {len(unique_im)} unique Im values:")
for v, c in sorted(zip(unique_im, im_counts)):
    print(f"  Im = {v:+.4f}  count {c}")

# === Pauli-basis decomposition: do axis modes live purely in n_XY=2? ===
all_alphas = list(itertools.product(range(4), repeat=N))
pauli_basis = np.zeros((dim * dim, dim * dim), dtype=complex)
for k, alpha in enumerate(all_alphas):
    pauli_basis[:, k] = pauli_string(alpha).flatten('F')
nxy_per_alpha = np.array([n_xy(a) for a in all_alphas])

nxy_supports = np.zeros((len(axis_indices), N + 1))
for i, idx in enumerate(axis_indices):
    v = eigvecs[:, idx]
    c = pauli_basis.conj().T @ v / (2 ** N)
    abs_c2 = np.abs(c) ** 2
    if abs_c2.sum() > 0:
        abs_c2 /= abs_c2.sum()
    for k in range(N + 1):
        nxy_supports[i, k] = abs_c2[nxy_per_alpha == k].sum()

pure_count = int(np.sum(nxy_supports[:, N // 2] > 0.9999))
print(f"\n# n_XY support of axis modes:")
print(f"  Axis modes with > 99.99% support on n_XY = {N//2} layer: "
      f"{pure_count} of {len(axis_indices)}")
print(f"  Mean Pauli-support distribution across axis modes:")
for k in range(N + 1):
    print(f"    n_XY = {k}: {nxy_supports.mean(axis=0)[k]:.4f}")

# === Im=0 axis modes: the full self-mirror kernel ===
im0_mask = mask & (np.abs(im) < 1e-9)
im0_indices = np.where(im0_mask)[0]
y_per_alpha = np.array([sum(1 for a in alpha if a == 2) for alpha in all_alphas])

print(f"\n# {'=' * 64}")
print(f"# Im=0 axis modes: λ = −Σγ + 0i exactly ({len(im0_indices)} of {len(axis_indices)} axis modes)")
print(f"# {'=' * 64}")

for i, idx in enumerate(im0_indices):
    v = eigvecs[:, idx]
    c = pauli_basis.conj().T @ v / (2 ** N)
    abs_c2 = np.abs(c) ** 2
    if abs_c2.sum() > 0:
        abs_c2 /= abs_c2.sum()
    nxy_dist = [abs_c2[nxy_per_alpha == k].sum() for k in range(N + 1)]
    y_dist = [abs_c2[y_per_alpha == k].sum() for k in range(N + 1)]

    abs_c = np.abs(c)
    abs_c_norm = abs_c / abs_c.max() if abs_c.max() > 0 else abs_c
    top = np.argsort(abs_c_norm)[::-1][:4]
    pauli_str = ""
    for k in top:
        if abs_c_norm[k] < 0.2:
            break
        name = ''.join(PAULI_NAMES[a] for a in all_alphas[k])
        pauli_str += f"σ_{name}({abs_c_norm[k]:.2f}) "

    nxy_str = "[" + " ".join(f"{x:.2f}" for x in nxy_dist) + "]"
    y_str = "[" + " ".join(f"{x:.2f}" for x in y_dist) + "]"
    print(f"  #{i:2d}: n_XY {nxy_str}  Y-count {y_str}")
    print(f"       top Pauli: {pauli_str}")

# === Aggregate over the 18 Im=0 modes ===
print(f"\n# Aggregate over the {len(im0_indices)} Im=0 modes:")
purity_count = 0
mean_nxy_dist = np.zeros(N + 1)
mean_y_dist = np.zeros(N + 1)
for idx in im0_indices:
    v = eigvecs[:, idx]
    c = pauli_basis.conj().T @ v / (2 ** N)
    abs_c2 = np.abs(c) ** 2
    if abs_c2.sum() > 0:
        abs_c2 /= abs_c2.sum()
    if abs_c2[nxy_per_alpha == N // 2].sum() > 0.9999:
        purity_count += 1
    for k in range(N + 1):
        mean_nxy_dist[k] += abs_c2[nxy_per_alpha == k].sum()
        mean_y_dist[k] += abs_c2[y_per_alpha == k].sum()
mean_nxy_dist /= len(im0_indices)
mean_y_dist /= len(im0_indices)
print(f"  Pure n_XY={N//2}: {purity_count} of {len(im0_indices)}")
print(f"  Mean n_XY support: {[f'{x:.3f}' for x in mean_nxy_dist]}")
print(f"  Mean Y-count support: {[f'{x:.3f}' for x in mean_y_dist]}")

# === Which Pauli strings carry the Im=0 subspace? ===
pauli_im0_content = np.zeros(len(all_alphas))
for idx in im0_indices:
    v = eigvecs[:, idx]
    c = pauli_basis.conj().T @ v / (2 ** N)
    pauli_im0_content += np.abs(c) ** 2

non_trivial = [(k, pauli_im0_content[k]) for k in range(len(all_alphas))
               if pauli_im0_content[k] > 1e-6]
non_trivial.sort(key=lambda x: -x[1])
print(f"\n# Pauli strings with non-trivial weight in Im=0 subspace: "
      f"{len(non_trivial)} of {4**N}")
print(f"  (Σ_modes |c_α|² across the {len(im0_indices)} Im=0 modes, top 30)")
for k, content in non_trivial[:30]:
    alpha = all_alphas[k]
    name = ''.join(PAULI_NAMES[a] for a in alpha)
    print(f"    σ_{name}  n_XY={nxy_per_alpha[k]}  n_Y={y_per_alpha[k]}  "
          f"weight = {content:.4f}")

# === Site-reflection (F71) parity of the 18 Im=0 modes ===
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
test_alpha = (1, 2, 3, 0)
assert np.allclose(pauli_string(test_alpha[::-1]),
                   R_h @ pauli_string(test_alpha) @ R_h), \
    "site-reflection inconsistent"

R_op = np.kron(R_h, R_h)  # vec(R ρ R) = (R ⊗ R) vec(ρ) since R is real involutive

# Theoretical R-parity split of the n_XY=2 layer:
#   palindromic strings (α == rev(α)): R-eigenvalue +1, contribute to even subspace.
#   non-palindromic pairs: each contributes one +1 and one −1 eigenvector.
nxy2_palindromic = 0
nxy2_pairs = 0
seen = set()
for k, alpha in enumerate(all_alphas):
    if nxy_per_alpha[k] != N // 2:
        continue
    if k in seen:
        continue
    rev_k = all_alphas.index(alpha[::-1])
    if rev_k == k:
        nxy2_palindromic += 1
        seen.add(k)
    else:
        nxy2_pairs += 1
        seen.add(k)
        seen.add(rev_k)
even_dim_nxy2 = nxy2_palindromic + nxy2_pairs
odd_dim_nxy2 = nxy2_pairs
print(f"\n# Theoretical R-parity decomposition of the n_XY={N//2} layer "
      f"({nxy2_palindromic} palindromic + {nxy2_pairs} non-palindromic pairs):")
print(f"  R-even subspace dim: {even_dim_nxy2}")
print(f"  R-odd  subspace dim: {odd_dim_nxy2}")
print(f"  total: {even_dim_nxy2 + odd_dim_nxy2} (= C(N,N/2)·2^N = {comb(N, N//2) * 2**N})")

# Project R_op onto the Im=0 kernel and diagonalize
V_K = eigvecs[:, im0_indices]
U_K, _, _ = np.linalg.svd(V_K, full_matrices=False)  # orthonormal basis of the kernel
R_sub = U_K.conj().T @ R_op @ U_K
R_sub_sym = (R_sub + R_sub.conj().T) / 2  # symmetrise numerical noise
r_eigvals = np.linalg.eigvalsh(R_sub_sym)

print(f"\n# Site-reflection eigenvalues within the {len(im0_indices)}-dim Im=0 kernel:")
print(f"  {np.round(r_eigvals, 6)}")
n_even = int(np.sum(np.abs(r_eigvals - 1) < 1e-6))
n_odd = int(np.sum(np.abs(r_eigvals + 1) < 1e-6))
print(f"  R-even (+1): {n_even}")
print(f"  R-odd  (−1): {n_odd}")
print(f"  ambiguous (|λ_R| ≠ 1): {len(r_eigvals) - n_even - n_odd}")

# === R-parity decomposition of the FULL 94-dim axis subspace ===
print(f"\n# {'=' * 64}")
print(f"# R-parity decomposition of the full {len(axis_indices)}-dim axis subspace")
print(f"# {'=' * 64}")

V_axis = eigvecs[:, axis_indices]
U_axis, _, _ = np.linalg.svd(V_axis, full_matrices=False)
R_axis = U_axis.conj().T @ R_op @ U_axis
R_axis_sym = (R_axis + R_axis.conj().T) / 2
r_eigvals_axis, r_eigvecs_axis = np.linalg.eigh(R_axis_sym)

n_even_axis = int(np.sum(np.abs(r_eigvals_axis - 1) < 1e-6))
n_odd_axis = int(np.sum(np.abs(r_eigvals_axis + 1) < 1e-6))
ambig_axis = len(r_eigvals_axis) - n_even_axis - n_odd_axis
print(f"  R-even axis modes: {n_even_axis}")
print(f"  R-odd  axis modes: {n_odd_axis}")
print(f"  ambig (|λ_R| ≠ 1): {ambig_axis}")
print(f"  total: {n_even_axis + n_odd_axis} of {len(axis_indices)} axis modes; "
      f"layer-predicted 52 R-even + 44 R-odd = 96, "
      f"{96 - len(axis_indices)} mode(s) leaked off-axis")

# Im(λ) within R-even axis subspace
even_mask = np.abs(r_eigvals_axis - 1) < 1e-6
U_even = U_axis @ r_eigvecs_axis[:, even_mask]
L_even_sub = U_even.conj().T @ L @ U_even
im_even = np.linalg.eigvals(L_even_sub).imag

print(f"\n# Im(λ) within R-even axis subspace ({n_even_axis} modes):")
for v, c in zip(*np.unique(np.round(im_even, 6), return_counts=True)):
    flag = " <-- silent" if abs(v) < 1e-6 else ""
    print(f"  Im = {v:+.4f}  count {c}{flag}")

# Im(λ) within R-odd axis subspace
odd_mask = np.abs(r_eigvals_axis + 1) < 1e-6
U_odd = U_axis @ r_eigvecs_axis[:, odd_mask]
L_odd_sub = U_odd.conj().T @ L @ U_odd
im_odd = np.linalg.eigvals(L_odd_sub).imag

print(f"\n# Im(λ) within R-odd axis subspace ({n_odd_axis} modes):")
for v, c in zip(*np.unique(np.round(im_odd, 6), return_counts=True)):
    flag = " <-- silent" if abs(v) < 1e-6 else ""
    print(f"  Im = {v:+.4f}  count {c}{flag}")
