"""
Nested Mirror Structure: Absorption Theorem Test
=================================================

Tests whether the three eigenvalue classes of the minimal qubit-in-qubit
nest inherit from the Absorption Theorem:

    Re(lambda) = -2 * Sum_k gamma_k * <n_XY>_k

where <n_XY>_k is the XY-Pauli-weight of the eigenmode at site k.

In our setup: gamma_S = 0, gamma_B = 0.1. So the prediction is:
    Re(lambda) = -2 * gamma_B * <n_XY>_B = -0.2 * <n_XY>_B

Per-class prediction:
    Class Re = 0     -> <n_XY>_B = 0    (pure I/Z content at B)
    Class Re = -0.1  -> <n_XY>_B = 0.5  (half XY content at B)
    Class Re = -0.2  -> <n_XY>_B = 1.0  (full XY content at B)

If this holds, the three-class structure is not a new algebraic object but
an inherited consequence of the Absorption Theorem applied per-site.

Date: 2026-04-14
Authors: Tom and Claude (chat)
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True)

# Pauli basis
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

paulis = [I2, X, Y, Z]
pauli_names = ['I', 'X', 'Y', 'Z']
# n_XY per site: 1 if pauli in {X,Y}, 0 otherwise
n_xy_site = [0, 1, 1, 0]  # I, X, Y, Z

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


# Build the 16 Pauli basis operators on 2 qubits
# Labels ordered (S-site, B-site): II, IX, IY, IZ, XI, XX, ...
basis_ops = []
basis_labels = []
basis_nxy_S = []
basis_nxy_B = []
for i in range(4):  # S-site pauli
    for j in range(4):  # B-site pauli
        basis_ops.append(np.kron(paulis[i], paulis[j]))
        basis_labels.append(pauli_names[i] + pauli_names[j])
        basis_nxy_S.append(n_xy_site[i])
        basis_nxy_B.append(n_xy_site[j])


def pauli_expand(M):
    """Expand 4x4 operator M in the 16 Pauli-string basis.
    c_ij = (1/4) * Tr[(sigma_i (x) sigma_j)^dagger @ M]
    Returns array of 16 complex coefficients.
    """
    coeffs = np.zeros(16, dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / 4.0
    return coeffs


def nxy_weights(M):
    """Compute average n_XY at S and B sites for operator M."""
    c = pauli_expand(M)
    weights = np.abs(c) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0, 0.0
    nxy_S = sum(w * ns for w, ns in zip(weights, basis_nxy_S)) / total
    nxy_B = sum(w * nb for w, nb in zip(weights, basis_nxy_B)) / total
    return nxy_S, nxy_B


# System setup (identical to qubit_in_qubit_layer_mirror.py)
J = 1.0
gamma_B = 0.1
H = J * 0.5 * (kron(X, X) + kron(Y, Y))
L_jump = np.sqrt(gamma_B) * kron(I2, Z)
L_super = liouvillian(H, [L_jump])

# Diagonalize
eigvals, eigvecs = np.linalg.eig(L_super)
order = sorted(range(16), key=lambda i: (eigvals[i].real, eigvals[i].imag))

print("Absorption Theorem test: per-mode <n_XY>_S and <n_XY>_B")
print("=" * 75)
print(f"Prediction: Re(lambda) = -2 * gamma_B * <n_XY>_B = -{2*gamma_B} * <n_XY>_B")
print()
print(f"  {'Re(lam)':>9} {'Im(lam)':>10}  {'<n_XY>_S':>10} {'<n_XY>_B':>10}  "
      f"{'predicted Re':>14}  {'match':>6}")
print("  " + "-" * 75)

results = []
for idx in order:
    lam = eigvals[idx]
    v = eigvecs[:, idx]
    M = v.reshape((4, 4), order='F')
    nxy_S, nxy_B = nxy_weights(M)
    predicted_re = -2 * gamma_B * nxy_B
    match = abs(lam.real - predicted_re) < 1e-6
    match_str = "EXACT" if match else f"d={abs(lam.real-predicted_re):.4f}"
    results.append((lam.real, lam.imag, nxy_S, nxy_B, predicted_re, match))
    print(f"  {lam.real:9.4f} {lam.imag:10.4f}  {nxy_S:10.4f} {nxy_B:10.4f}  "
          f"{predicted_re:14.4f}  {match_str:>6}")

# Per-class summary
print()
print("Per-class averages:")
print("  " + "-" * 60)
classes = {'Re=0 (conserved)': [], 'Re=-0.1 (mirror)': [], 'Re=-0.2 (correlation)': []}
for lam_r, lam_i, nxy_S, nxy_B, pred, match in results:
    if abs(lam_r) < 1e-6:
        classes['Re=0 (conserved)'].append((nxy_S, nxy_B))
    elif abs(lam_r + 0.1) < 1e-6:
        classes['Re=-0.1 (mirror)'].append((nxy_S, nxy_B))
    elif abs(lam_r + 0.2) < 1e-6:
        classes['Re=-0.2 (correlation)'].append((nxy_S, nxy_B))

for cls, items in classes.items():
    if items:
        avg_S = np.mean([s for s, b in items])
        avg_B = np.mean([b for s, b in items])
        count = len(items)
        print(f"  {cls}: count={count}  avg <n_XY>_S={avg_S:.4f}  "
              f"avg <n_XY>_B={avg_B:.4f}")

# Universal check
all_match = all(r[5] for r in results)
print()
print(f"Absorption theorem exact match for all 16 modes: {all_match}")
if all_match:
    print("The three-class structure IS the three-level quantization of <n_XY>_B.")
    print("Nested mirror structure = Absorption Theorem applied per site.")
else:
    print("Partial match. The three-class structure is related to but not")
    print("identical with per-site <n_XY> quantization.")
