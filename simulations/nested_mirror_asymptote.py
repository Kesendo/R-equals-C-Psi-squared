"""
Inside-Outside Correspondence: Probe 3, the 1/4 Asymptote
===========================================================

At N=3 with J >> gamma, the refraction test showed <n_XY>_B values
converging to multiples of 1/4. Is this exact? Is it universal across N?

Method:
  - N=3 chain S-M-B at J = 10, 20, 50, 100, 500, fixed gamma_B = 0.1
  - Compute deviation of <n_XY>_B from exact {0, 1/4, 1/2, 3/4, 1}
  - Fit scaling: deviation ~ 1/J? 1/J^2? exp(-J)?
  - Test N=4 (4^4 = 256 dim) to see if 1/4 persists or changes

Connection to ZERO_IS_THE_MIRROR: at gamma=0, the palindrome center is 0.
At J >> gamma, the system approaches the unitary limit where J dominates.
If the <n_XY>_B quantization at J -> inf is the CAVITY STRUCTURE of the
unitary system, then the 1/4 levels are the standing wave nodes of the
Hamiltonian's eigenmodes projected onto the dissipative site.

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from pathlib import Path

np.set_printoptions(precision=8, suppress=True)

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
    """Build all 4^n Pauli string operators and their per-site n_XY."""
    if n_qubits == 1:
        return list(paulis), [[s] for s in n_xy_site]
    ops = []
    nxy_per_site = []
    sub_ops, sub_nxy = build_pauli_basis(n_qubits - 1)
    for i in range(4):
        for j, (sub_op, sub_n) in enumerate(zip(sub_ops, sub_nxy)):
            ops.append(np.kron(paulis[i], sub_op))
            nxy_per_site.append([n_xy_site[i]] + sub_n)
    return ops, nxy_per_site


def nxy_weight_at_site(M, site, basis_ops, basis_nxy):
    """Compute <n_XY> at a specific site for operator M."""
    d = int(np.sqrt(M.shape[0]))
    coeffs = np.zeros(len(basis_ops), dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / d
    weights = np.abs(coeffs) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0
    return sum(w * nxy[site] for w, nxy in zip(weights, basis_nxy)) / total


# =========================================================================
# N=3: Asymptote test at large J
# =========================================================================
print("=" * 72)
print("Probe 3: 1/4 Asymptote at N=3")
print("=" * 72)

gamma_B = 0.1
J_values = [0.1, 1.0, 10.0, 20.0, 50.0, 100.0, 500.0]

# Build N=3 Pauli basis
basis_3, nxy_3 = build_pauli_basis(3)

results_dir = Path("simulations/results/inside_outside_asymptote")
log_lines = []

# Target quantization levels
targets = [0.0, 0.25, 0.5, 0.75, 1.0]

print(f"\nTarget levels: {targets}")
print(f"gamma_B = {gamma_B}\n")

asymptote_data = []

for J in J_values:
    H = J * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2)
                  + kron(I2, X, X) + kron(I2, Y, Y))
    L_jump = np.sqrt(gamma_B) * kron(I2, I2, Z)
    L = liouvillian(H, [L_jump])

    eigvals, eigvecs = np.linalg.eig(L)

    # Compute <n_XY>_B (site 2) for each mode
    nxy_B_values = []
    for k in range(64):
        M = eigvecs[:, k].reshape((8, 8), order='F')
        nB = nxy_weight_at_site(M, 2, basis_3, nxy_3)
        nxy_B_values.append(nB)

    # Find distinct values
    nxy_sorted = sorted(nxy_B_values)
    distinct = [nxy_sorted[0]]
    for v in nxy_sorted[1:]:
        if abs(v - distinct[-1]) > 1e-6:
            distinct.append(v)

    # Compute deviation from nearest target
    max_dev = 0
    mean_dev = 0
    for v in distinct:
        nearest = min(targets, key=lambda t: abs(t - v))
        dev = abs(v - nearest)
        max_dev = max(max_dev, dev)
        mean_dev += dev
    mean_dev /= len(distinct)

    line = (f"J={J:>7.1f}  distinct={len(distinct):2d}  "
            f"max_dev={max_dev:.2e}  mean_dev={mean_dev:.2e}  "
            f"levels=[{', '.join(f'{v:.4f}' for v in distinct)}]")
    print(line)
    log_lines.append(line)

    asymptote_data.append({
        'J': J, 'n_distinct': len(distinct),
        'max_dev': max_dev, 'mean_dev': mean_dev,
        'levels': distinct
    })

# Scaling fit: max_dev vs J
print("\n" + "-" * 72)
print("Scaling: max_dev vs J")
print("-" * 72)

# Use points where J >= 10 for fit
fit_data = [(d['J'], d['max_dev']) for d in asymptote_data if d['J'] >= 10 and d['max_dev'] > 0]
if len(fit_data) >= 3:
    Js = np.array([d[0] for d in fit_data])
    devs = np.array([d[1] for d in fit_data])

    # Try 1/J scaling
    ratio_1 = devs * Js
    print(f"dev * J:   {[f'{r:.4f}' for r in ratio_1]}  (constant if 1/J scaling)")

    # Try 1/J^2 scaling
    ratio_2 = devs * Js**2
    print(f"dev * J^2: {[f'{r:.4f}' for r in ratio_2]}  (constant if 1/J^2 scaling)")

    # Try (gamma/J)^2
    ratio_gJ = devs / (gamma_B / Js)**2
    print(f"dev / (g/J)^2: {[f'{r:.4f}' for r in ratio_gJ]}  (constant if (gamma/J)^2 scaling)")

    # Best fit
    cv_1 = np.std(ratio_1) / np.mean(ratio_1) if np.mean(ratio_1) > 0 else 999
    cv_2 = np.std(ratio_2) / np.mean(ratio_2) if np.mean(ratio_2) > 0 else 999
    cv_gJ = np.std(ratio_gJ) / np.mean(ratio_gJ) if np.mean(ratio_gJ) > 0 else 999
    print(f"\nCoefficient of variation: 1/J: {cv_1:.4f}, 1/J^2: {cv_2:.4f}, (g/J)^2: {cv_gJ:.4f}")
    best = min([('1/J', cv_1), ('1/J^2', cv_2), ('(gamma/J)^2', cv_gJ)], key=lambda x: x[1])
    print(f"Best scaling: {best[0]} (CV = {best[1]:.4f})")

# =========================================================================
# N=4: Does 1/4 persist?
# =========================================================================
print("\n" + "=" * 72)
print("N=4 test: does the 1/4 quantization persist?")
print("=" * 72)

basis_4, nxy_4 = build_pauli_basis(4)
dim4 = 4**4  # 256

for J in [1.0, 10.0, 100.0]:
    H = J * 0.5 * (kron(X, X, I2, I2) + kron(Y, Y, I2, I2)
                  + kron(I2, X, X, I2) + kron(I2, Y, Y, I2)
                  + kron(I2, I2, X, X) + kron(I2, I2, Y, Y))
    L_jump = np.sqrt(gamma_B) * kron(I2, I2, I2, Z)
    L = liouvillian(H, [L_jump])

    eigvals, eigvecs = np.linalg.eig(L)

    nxy_B_values = []
    for k in range(dim4):
        M = eigvecs[:, k].reshape((16, 16), order='F')
        nB = nxy_weight_at_site(M, 3, basis_4, nxy_4)
        nxy_B_values.append(nB)

    nxy_sorted = sorted(nxy_B_values)
    distinct = [nxy_sorted[0]]
    for v in nxy_sorted[1:]:
        if abs(v - distinct[-1]) > 1e-5:
            distinct.append(v)

    # Check if levels are multiples of 1/4 or 1/N or something else
    targets_4 = [k/4 for k in range(5)]     # 1/4 multiples
    targets_N = [k/4 for k in range(5)]     # same for N=4 chain with 4 sites

    max_dev_4 = max(min(abs(v - t) for t in targets_4) for v in distinct)

    # Also check 1/6 multiples (N=4 chain has 3 bonds)
    targets_6 = [k/6 for k in range(7)]
    max_dev_6 = max(min(abs(v - t) for t in targets_6) for v in distinct)

    # Check 1/8
    targets_8 = [k/8 for k in range(9)]
    max_dev_8 = max(min(abs(v - t) for t in targets_8) for v in distinct)

    line = (f"N=4, J={J:>6.1f}: {len(distinct):2d} distinct levels  "
            f"max_dev(1/4)={max_dev_4:.4f}  max_dev(1/6)={max_dev_6:.4f}  "
            f"max_dev(1/8)={max_dev_8:.4f}")
    print(line)
    log_lines.append(line)

    if len(distinct) <= 20:
        levels_str = ', '.join(f'{v:.4f}' for v in distinct)
        print(f"  levels: [{levels_str}]")
        log_lines.append(f"  levels: [{levels_str}]")

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 72)
print("VERDICT: 1/4 Asymptote")
print("=" * 72)

verdict = []
# N=3 at large J
last = asymptote_data[-1]
if last['max_dev'] < 0.01:
    verdict.append(f"N=3, J={last['J']}: max deviation from 1/4 multiples = {last['max_dev']:.2e}")
    verdict.append("1/4 quantization is approached at large J/gamma.")
else:
    verdict.append(f"N=3, J={last['J']}: max deviation = {last['max_dev']:.4f}, NOT converging to 1/4.")

for line in verdict:
    print(line)

# Save
with open(results_dir / 'asymptote_results.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 3: 1/4 Asymptote Results\n")
    f.write("=" * 72 + "\n\n")
    f.write('\n'.join(log_lines))
    f.write("\n\nVerdict:\n")
    f.write('\n'.join(verdict))

print(f"\nResults saved to {results_dir / 'asymptote_results.txt'}")
