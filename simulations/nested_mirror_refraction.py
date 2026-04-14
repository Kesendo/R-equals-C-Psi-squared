"""
Nested Mirror Structure: Refraction Test
=========================================

Test Tom's intuition: the 12 <n_XY>_B values at N=3 are refraction of
the basic {0, 1} quantization by the J-coupling. At J=0 the chain has no
bonds, no hopping, no refraction - we should see only 2 distinct
<n_XY>_B values. As J increases, the values split into finer structure.

Prediction:
    J=0:    2 distinct <n_XY>_B values: {0, 1}
    J small: 2 values slightly broadened
    J=1.0:  12 distinct values (our standard case)
    J large: further splitting or saturation

Date: 2026-04-14
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
paulis = [I2, X, Y, Z]
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


# Pauli basis on 3 qubits
basis_ops = []
basis_nxy_B = []
for i in range(4):
    for j in range(4):
        for k in range(4):
            basis_ops.append(kron(paulis[i], paulis[j], paulis[k]))
            basis_nxy_B.append(n_xy_site[k])


def nxy_B(M):
    coeffs = np.zeros(64, dtype=complex)
    for m, P in enumerate(basis_ops):
        coeffs[m] = np.trace(P.conj().T @ M) / 8.0
    w = np.abs(coeffs) ** 2
    total = w.sum()
    if total < 1e-14:
        return 0.0
    return sum(w[m] * basis_nxy_B[m] for m in range(64)) / total



gamma_B = 0.1
XX = kron(X, X)
YY = kron(Y, Y)

print("Refraction test: how many distinct <n_XY>_B values at each J?")
print("=" * 70)
print(f"Fixed: gamma_B = {gamma_B}, 3-qubit chain S-M-B")
print()

J_values = [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

for J in J_values:
    H_SM = kron(XX + YY, I2) * 0.5 * J
    H_MB = kron(I2, XX + YY) * 0.5 * J
    H = H_SM + H_MB
    L_jump = np.sqrt(gamma_B) * kron(I2, I2, Z)
    L = liouvillian(H, [L_jump])
    eigvals, eigvecs = np.linalg.eig(L)

    nxy_vals = []
    for idx in range(64):
        v = eigvecs[:, idx]
        M = v.reshape((8, 8), order='F')
        nxy_vals.append(nxy_B(M))

    # Round to 4 decimal places for classification
    unique = sorted(set(round(v, 4) for v in nxy_vals))
    n_unique = len(unique)
    # Unique Re(lambda) values
    re_unique = sorted(set(round(e.real, 4) for e in eigvals))
    n_re = len(re_unique)

    if n_unique <= 6:
        vals_str = str(unique)
    else:
        vals_str = f"{unique[0]}, {unique[1]}, ..., {unique[-2]}, {unique[-1]}"

    print(f"  J = {J:6.3f}:  {n_unique:2d} distinct <n_XY>_B, "
          f"{n_re:2d} distinct Re(lambda)")
    print(f"              values: {vals_str}")
    print()
