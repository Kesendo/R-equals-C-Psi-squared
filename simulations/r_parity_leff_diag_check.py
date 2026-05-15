"""F86 standard L_eff Q-sweep vs R-parity block Q-sweep at N=4.

Two parallel calculations:
1. Standard F86 L_eff: project L onto channel-uniform basis {|c_1⟩, |c_3⟩},
   find 2x2 eigenvalues vs Q. EP at Q_peak = 2/g_eff (≈ 1.15 for Endpoint N=4).
2. R-odd L_eff candidate: project L onto a chosen R-antisymmetric basis at
   HD=1 and HD=3, find 2x2 eigenvalues vs Q. Check for an R-odd EP at
   the same or different Q value.

The R-antisymmetric basis vectors at HD=k are constructed as the slowest
R-odd L-eigenmode within the HD=k subspace at γ=0 perturbed to small Q
(i.e., the analogue of channel-uniform but R-antisymmetric).

Run: python simulations/_r_parity_q_sweep_leff.py
"""
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

N = 4
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


def build_L(N, gamma, J):
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


def channel_uniform_vec(N, hd):
    """vec(|c_k⟩) = uniform sum over (a, b) with HD(a, b) = hd.
    Column-stacking convention: vec(|a⟩⟨b|)[b*dim + a] = 1.
    """
    v = np.zeros(op_dim, dtype=complex)
    count = 0
    for a in range(dim):
        for b in range(dim):
            if bin(a ^ b).count('1') == hd:
                v[b * dim + a] = 1.0
                count += 1
    return v / np.sqrt(count)


R_h = np.zeros((dim, dim), dtype=float)
for b in range(dim):
    R_h[reverse_bits(b, N), b] = 1.0
R_op = np.kron(R_h, R_h)

# Channel-uniform basis
c1 = channel_uniform_vec(N, 1)
c3 = channel_uniform_vec(N, 3)

# Verify R-evenness
r1 = R_op @ c1
r3 = R_op @ c3
print(f"# Channel-uniform R-eigenvalues:")
print(f"#   ⟨c_1|R|c_1⟩ = {c1.conj() @ r1:.4f}")
print(f"#   ⟨c_3|R|c_3⟩ = {c3.conj() @ r3:.4f}")

# === Part 1: Standard F86 L_eff (channel-uniform projection) ===
print(f"\n# Part 1: Standard F86 L_eff (channel-uniform projection) at N={N}, γ={gamma}")
print(f"#")
print(f"#   Q    | λ_+ (Re, Im)        | λ_- (Re, Im)        | |Re gap| | |Im gap|")
print(f"#   -----|---------------------|---------------------|----------|---------")

Q_grid = np.linspace(0.3, 3.0, 28)
for Q in Q_grid:
    J = Q * gamma
    L = build_L(N, gamma, J)
    L11 = c1.conj() @ L @ c1
    L13 = c1.conj() @ L @ c3
    L31 = c3.conj() @ L @ c1
    L33 = c3.conj() @ L @ c3
    L_eff = np.array([[L11, L13], [L31, L33]])
    eigs = np.linalg.eigvals(L_eff)
    eigs = sorted(eigs, key=lambda e: -e.real)
    re_gap = abs(eigs[0].real - eigs[1].real)
    im_gap = abs(eigs[0].imag - eigs[1].imag)
    print(f"   {Q:.2f} | ({eigs[0].real:+.5f}, {eigs[0].imag:+.5f}) | "
          f"({eigs[1].real:+.5f}, {eigs[1].imag:+.5f}) | {re_gap:.5f}  | {im_gap:.5f}")

# === Part 2: R-odd analogue (find R-odd vectors with slowest L-decay at HD=k) ===
print(f"\n# Part 2: R-odd analogue L_eff candidate")
print(f"#")

# Build R-projection onto HD=k subspace, find slowest R-odd L-eigenmode there
def hd_mask(N, hd):
    """Return projector onto operators with HD(a,b) = hd (column-stacking)."""
    mask = np.zeros(op_dim, dtype=bool)
    for a in range(dim):
        for b in range(dim):
            if bin(a ^ b).count('1') == hd:
                mask[b * dim + a] = True
    return mask


# Compute R-eigenbasis on full op space
r_vals, r_vecs = np.linalg.eigh(R_op)
odd_basis = r_vecs[:, np.abs(r_vals + 1) < 1e-6]

# Find slowest R-odd L-eigenmode at HD=1, HD=3 at small Q
print(f"#   Finding slowest R-odd L-eigenmode per HD layer at Q = 0.5...")
J_probe = 0.5 * gamma
L_probe = build_L(N, gamma, J_probe)
L_odd_probe = odd_basis.T @ L_probe @ odd_basis
eig_odd_vals, eig_odd_vecs = np.linalg.eig(L_odd_probe)
# Lift R-odd eigenvectors to full op space
eig_odd_full = odd_basis @ eig_odd_vecs

# For each R-odd eigenmode, compute HD-content
hd1_proj = hd_mask(N, 1)
hd3_proj = hd_mask(N, 3)
hd_content = np.zeros((len(eig_odd_vals), 5))  # HD=0..4
for k in range(5):
    mask = hd_mask(N, k)
    for i in range(len(eig_odd_vals)):
        v = eig_odd_full[:, i]
        hd_content[i, k] = np.sum(np.abs(v[mask]) ** 2)

# Find slowest R-odd mode with dominant HD=1 content (>0.5)
slowest_idx = np.argsort(np.abs(eig_odd_vals.real))
hd1_candidate = None
hd3_candidate = None
for idx in slowest_idx:
    if hd_content[idx, 1] > 0.5 and hd1_candidate is None:
        hd1_candidate = idx
    if hd_content[idx, 3] > 0.5 and hd3_candidate is None:
        hd3_candidate = idx
    if hd1_candidate is not None and hd3_candidate is not None:
        break

if hd1_candidate is None or hd3_candidate is None:
    print(f"#   Could not find dominant-HD R-odd modes; using slowest two.")
    hd1_candidate = slowest_idx[0]
    hd3_candidate = slowest_idx[1]

v1_odd = eig_odd_full[:, hd1_candidate]
v3_odd = eig_odd_full[:, hd3_candidate]
print(f"#   R-odd HD=1 candidate: eigval = {eig_odd_vals[hd1_candidate]:.5f}, "
      f"HD-content = {[f'{x:.2f}' for x in hd_content[hd1_candidate]]}")
print(f"#   R-odd HD=3 candidate: eigval = {eig_odd_vals[hd3_candidate]:.5f}, "
      f"HD-content = {[f'{x:.2f}' for x in hd_content[hd3_candidate]]}")

# Now sweep Q with this fixed R-odd basis
print(f"\n#   Q    | R-odd λ_+ (Re, Im)  | R-odd λ_- (Re, Im)  | |Re gap| | |Im gap|")
print(f"#   -----|---------------------|---------------------|----------|---------")
for Q in Q_grid:
    J = Q * gamma
    L = build_L(N, gamma, J)
    L11_o = v1_odd.conj() @ L @ v1_odd
    L13_o = v1_odd.conj() @ L @ v3_odd
    L31_o = v3_odd.conj() @ L @ v1_odd
    L33_o = v3_odd.conj() @ L @ v3_odd
    Leff_o = np.array([[L11_o, L13_o], [L31_o, L33_o]])
    eigs_o = np.linalg.eigvals(Leff_o)
    eigs_o = sorted(eigs_o, key=lambda e: -e.real)
    re_gap_o = abs(eigs_o[0].real - eigs_o[1].real)
    im_gap_o = abs(eigs_o[0].imag - eigs_o[1].imag)
    print(f"   {Q:.2f} | ({eigs_o[0].real:+.5f}, {eigs_o[0].imag:+.5f}) | "
          f"({eigs_o[1].real:+.5f}, {eigs_o[1].imag:+.5f}) | {re_gap_o:.5f}  | {im_gap_o:.5f}")
