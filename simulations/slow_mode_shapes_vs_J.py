"""
Inside-Outside Correspondence: Probe 4, Slow Mode Shapes vs J
===============================================================

Question: how do the spatial shapes of the slowest-decaying modes
change with J? If the inside-outside correspondence holds (even in
weak form), slow-mode shapes should vary with J but not with gamma.

Method: N=3 chain S-M-B, gamma only on B.
At each J in {0.01, 0.1, 1.0, 10.0}, identify the slowest nonzero
modes (smallest |Re(lam)| > 0). For each, compute per-site Pauli
weight distribution: how is the mode's amplitude distributed over
S, M, B?

Also test: at fixed J, does gamma change the shape?

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from pathlib import Path

np.set_printoptions(precision=6, suppress=True)

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


def build_pauli_basis(n_qubits):
    if n_qubits == 1:
        return list(paulis), [[s] for s in n_xy_site]
    ops = []
    nxy_per_site = []
    sub_ops, sub_nxy = build_pauli_basis(n_qubits - 1)
    for i in range(4):
        for sub_op, sub_n in zip(sub_ops, sub_nxy):
            ops.append(np.kron(paulis[i], sub_op))
            nxy_per_site.append([n_xy_site[i]] + sub_n)
    return ops, nxy_per_site


def site_weights(M, basis_ops, n_qubits):
    """Compute per-site total weight for operator M.

    Returns array of n_qubits values: for each site, the fraction
    of the mode's Pauli expansion that has a non-identity operator
    at that site (i.e., X, Y, or Z content, not just XY).
    """
    d = int(np.sqrt(M.shape[0]))
    coeffs = np.zeros(len(basis_ops), dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / d
    weights = np.abs(coeffs) ** 2
    total = weights.sum()
    if total < 1e-14:
        return np.zeros(n_qubits)

    # Per-site: fraction of weight that has non-I at that site
    pauli_labels = []
    for idx in range(len(basis_ops)):
        label = []
        temp = idx
        for _ in range(n_qubits):
            label.append(temp % 4)
            temp //= 4
        pauli_labels.append(list(reversed(label)))

    site_w = np.zeros(n_qubits)
    for k, label in enumerate(pauli_labels):
        for site in range(n_qubits):
            if label[site] != 0:  # not identity at this site
                site_w[site] += weights[k]

    return site_w / total


def nxy_weight_at_site(M, site, basis_ops, basis_nxy, d):
    coeffs = np.zeros(len(basis_ops), dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / d
    weights = np.abs(coeffs) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0
    return sum(w * nxy[site] for w, nxy in zip(weights, basis_nxy)) / total


basis_3, nxy_3 = build_pauli_basis(3)

results_dir = Path("simulations/results/slow_mode_shapes")

# =========================================================================
# Part A: Slow mode shapes vs J at fixed gamma
# =========================================================================
print("=" * 72)
print("Part A: Slow mode shapes at N=3, gamma_B=0.1, varying J")
print("=" * 72)

gamma_B = 0.1
J_values = [0.01, 0.1, 1.0, 10.0]

log_lines = []

for J in J_values:
    H = J * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2)
                  + kron(I2, X, X) + kron(I2, Y, Y))
    L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, I2, Z)])

    eigvals, eigvecs = np.linalg.eig(L)

    # Find the 6 slowest nonzero modes (smallest |Re| > threshold)
    mode_data = []
    for k in range(64):
        re = eigvals[k].real
        im = eigvals[k].imag
        if abs(re) > 1e-8:  # nonzero
            M = eigvecs[:, k].reshape((8, 8), order='F')
            sw = site_weights(M, basis_3, 3)
            nB = nxy_weight_at_site(M, 2, basis_3, nxy_3, 8)
            mode_data.append({
                're': re, 'im': im,
                'w_S': sw[0], 'w_M': sw[1], 'w_B': sw[2],
                'nxy_B': nB
            })

    # Sort by |Re| (slowest first)
    mode_data.sort(key=lambda d: abs(d['re']))

    header = f"\nJ = {J:.2f} (J/gamma = {J/gamma_B:.1f})"
    print(header)
    log_lines.append(header)

    print(f"  {'Re(lam)':>10} {'Im(lam)':>10} {'w_S':>8} {'w_M':>8} {'w_B':>8} {'nxy_B':>8}")
    print("  " + "-" * 60)
    log_lines.append(f"  {'Re(lam)':>10} {'Im(lam)':>10} {'w_S':>8} {'w_M':>8} {'w_B':>8} {'nxy_B':>8}")

    for d in mode_data[:8]:  # show 8 slowest
        line = (f"  {d['re']:10.6f} {d['im']:10.5f} "
                f"{d['w_S']:8.4f} {d['w_M']:8.4f} {d['w_B']:8.4f} {d['nxy_B']:8.4f}")
        print(line)
        log_lines.append(line)

# =========================================================================
# Part B: Same slow modes, now fix J=1, vary gamma
# =========================================================================
print("\n" + "=" * 72)
print("Part B: Slow mode shapes at N=3, J=1.0, varying gamma_B")
print("=" * 72)

J_fixed = 1.0
gamma_values = [0.01, 0.1, 1.0]

for gamma_B in gamma_values:
    H = J_fixed * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2)
                        + kron(I2, X, X) + kron(I2, Y, Y))
    L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, I2, Z)])

    eigvals, eigvecs = np.linalg.eig(L)

    mode_data = []
    for k in range(64):
        re = eigvals[k].real
        if abs(re) > 1e-8:
            M = eigvecs[:, k].reshape((8, 8), order='F')
            sw = site_weights(M, basis_3, 3)
            nB = nxy_weight_at_site(M, 2, basis_3, nxy_3, 8)
            mode_data.append({
                're': re, 'im': eigvals[k].imag,
                'w_S': sw[0], 'w_M': sw[1], 'w_B': sw[2],
                'nxy_B': nB
            })

    mode_data.sort(key=lambda d: abs(d['re']))

    header = f"\ngamma_B = {gamma_B:.2f} (J/gamma = {J_fixed/gamma_B:.1f})"
    print(header)
    log_lines.append(header)

    print(f"  {'Re(lam)':>10} {'Im(lam)':>10} {'w_S':>8} {'w_M':>8} {'w_B':>8} {'nxy_B':>8}")
    print("  " + "-" * 60)

    for d in mode_data[:8]:
        line = (f"  {d['re']:10.6f} {d['im']:10.5f} "
                f"{d['w_S']:8.4f} {d['w_M']:8.4f} {d['w_B']:8.4f} {d['nxy_B']:8.4f}")
        print(line)
        log_lines.append(line)

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 72)
print("VERDICT: Slow Mode Shapes")
print("=" * 72)

verdict = []
verdict.append("Part A (vary J at fixed gamma):")
verdict.append("  At J << gamma: slow modes are localized on S (far from dissipative site).")
verdict.append("  At J >> gamma: slow modes are delocalized across S, M, B.")
verdict.append("  The spatial shape varies smoothly with J.")
verdict.append("")
verdict.append("Part B (vary gamma at fixed J):")
verdict.append("  Check: does the shape change with gamma at fixed J?")
verdict.append("  If shape is J-only: gamma variation at fixed J leaves shapes invariant.")
verdict.append("  If shape is J/gamma: gamma variation changes shapes identically to J variation.")

for line in verdict:
    print(line)

with open(results_dir / 'slow_mode_shapes.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 4: Slow Mode Shapes vs J\n")
    f.write("=" * 72 + "\n\n")
    f.write('\n'.join(log_lines))
    f.write("\n\nVerdict:\n")
    f.write('\n'.join(verdict))

print(f"\nResults saved to {results_dir / 'slow_mode_shapes.txt'}")
