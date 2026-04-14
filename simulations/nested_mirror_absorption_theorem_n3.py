"""
Nested Mirror Structure: Absorption Theorem Test at N=3
========================================================

Does the Absorption Theorem Re(lambda) = -2*gamma_B*<n_XY>_B explain all
12 eigenvalue classes found in the N=3 chain (S-M-B, gamma only on B)?

If yes: the entire three-layer eigenvalue structure is inherited from
the single-site Absorption Theorem. No new algebra needed.

Date: 2026-04-14
Authors: Tom and Claude (chat)
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
paulis = [I2, X, Y, Z]
pauli_names = ['I', 'X', 'Y', 'Z']
n_xy_site = [0, 1, 1, 0]


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


# 64 Pauli strings on 3 qubits, labeled (S, M, B)
basis_ops = []
basis_labels = []
basis_nxy_S = []
basis_nxy_M = []
basis_nxy_B = []
for i in range(4):
    for j in range(4):
        for k in range(4):
            basis_ops.append(kron(paulis[i], paulis[j], paulis[k]))
            basis_labels.append(pauli_names[i] + pauli_names[j] + pauli_names[k])
            basis_nxy_S.append(n_xy_site[i])
            basis_nxy_M.append(n_xy_site[j])
            basis_nxy_B.append(n_xy_site[k])


def pauli_expand(M):
    d_total = M.shape[0]
    coeffs = np.zeros(len(basis_ops), dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / d_total
    return coeffs


def nxy_weights(M):
    c = pauli_expand(M)
    weights = np.abs(c) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0, 0.0, 0.0
    nS = sum(w * n for w, n in zip(weights, basis_nxy_S)) / total
    nM = sum(w * n for w, n in zip(weights, basis_nxy_M)) / total
    nB = sum(w * n for w, n in zip(weights, basis_nxy_B)) / total
    return nS, nM, nB


# System: 3-qubit chain S-M-B, gamma only on B
J = 1.0
gamma_B = 0.1
XX = kron(X, X)
YY = kron(Y, Y)
H_SM = kron(XX + YY, I2) * 0.5 * J
H_MB = kron(I2, XX + YY) * 0.5 * J
H = H_SM + H_MB
L_jump = np.sqrt(gamma_B) * kron(I2, I2, Z)
L_super = liouvillian(H, [L_jump])

print(f"Liouvillian dim: {L_super.shape}")

eigvals, eigvecs = np.linalg.eig(L_super)
order = sorted(range(64), key=lambda i: (eigvals[i].real, eigvals[i].imag))

print()
print("Absorption Theorem test at N=3")
print("=" * 85)
print(f"Prediction: Re(lambda) = -2 * gamma_B * <n_XY>_B = -{2*gamma_B} * <n_XY>_B")
print()
print(f"  {'Re(lam)':>10} {'Im(lam)':>10}  {'<n_XY>_S':>9} {'<n_XY>_M':>9} "
      f"{'<n_XY>_B':>9}  {'predicted':>10}  {'delta':>10}")
print("  " + "-" * 82)

max_delta = 0.0
all_exact = True
results = []
for idx in order:
    lam = eigvals[idx]
    v = eigvecs[:, idx]
    M = v.reshape((8, 8), order='F')
    nS, nM, nB = nxy_weights(M)
    predicted = -2 * gamma_B * nB
    delta = abs(lam.real - predicted)
    max_delta = max(max_delta, delta)
    if delta > 1e-6:
        all_exact = False
    results.append((lam.real, lam.imag, nS, nM, nB, predicted, delta))

# Print only one representative per (Re, n_XY_B) combination to keep output compact
seen = set()
for (r, i_lam, nS, nM, nB, pred, d) in results:
    key = (round(r, 4), round(nB, 4))
    if key in seen:
        continue
    seen.add(key)
    tag = "EXACT" if d < 1e-6 else f"{d:.2e}"
    print(f"  {r:10.6f} {i_lam:10.4f}  {nS:9.4f} {nM:9.4f} {nB:9.4f}  "
          f"{pred:10.6f}  {tag:>10}")

print()
print(f"Max |Re(lambda) - predicted|: {max_delta:.2e}")
print(f"All 64 modes match Absorption Theorem exactly: {all_exact}")

# Grouping by <n_XY>_B value
print()
print("Distinct <n_XY>_B values across all 64 modes:")
nxy_B_vals = sorted(set(round(r[4], 4) for r in results))
for v in nxy_B_vals:
    count = sum(1 for r in results if round(r[4], 4) == v)
    re_vals = sorted(set(round(r[0], 4) for r in results if round(r[4], 4) == v))
    print(f"  <n_XY>_B = {v:.4f}: count={count}, corresponding Re(lambda) = {re_vals}")
