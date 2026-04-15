"""
Inside-Outside Correspondence: Probe 1, Lifetime-Structure Separation
======================================================================

Central question: does the Liouvillian spectrum decompose into two
independently parameterized families?
  - Lifetimes (Re(lam)) parameterized by gamma, invariant under J
  - Spatial structure (Pauli weights) parameterized by J, invariant under gamma

If yes: inside-outside correspondence is mathematically supported.
If no: both parameters act on both quantities, correspondence is metaphor.

Method: sweep gamma_B and J independently over a decade each at N=2.
For each (gamma_B, J) pair, extract per eigenmode:
  - Re(lam) (lifetime)
  - <n_XY>_S, <n_XY>_B (Pauli XY-weight per site)
  - Per-site Pauli expansion coefficients

Test the Absorption Theorem prediction: Re(lam) = -2gamma_B · <n_XY>_B.
If <n_XY>_B is J-dependent, then Re(lam) depends on BOTH gamma and J,
and the separation fails.

Date: 2026-04-14
Authors: Tom and Claude (Code)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

np.set_printoptions(precision=6, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

paulis = [I2, X, Y, Z]
pauli_names = ['I', 'X', 'Y', 'Z']
n_xy_site = [0, 1, 1, 0]

basis_ops = []
basis_nxy_S = []
basis_nxy_B = []
for i in range(4):
    for j in range(4):
        basis_ops.append(np.kron(paulis[i], paulis[j]))
        basis_nxy_S.append(n_xy_site[i])
        basis_nxy_B.append(n_xy_site[j])


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


def pauli_expand(M):
    coeffs = np.zeros(16, dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / 4.0
    return coeffs


def nxy_weights(M):
    c = pauli_expand(M)
    weights = np.abs(c) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0, 0.0
    nxy_S = sum(w * ns for w, ns in zip(weights, basis_nxy_S)) / total
    nxy_B = sum(w * nb for w, nb in zip(weights, basis_nxy_B)) / total
    return nxy_S, nxy_B


results_dir = Path("simulations/results/inside_outside_n2")

# =========================================================================
# Sweep 1: Fix J=1.0, vary gamma_B over a decade
# =========================================================================
print("=" * 72)
print("Sweep 1: Fix J=1.0, vary gamma_B from 0.01 to 1.0")
print("Question: do <n_XY> values stay invariant under gamma?")
print("=" * 72)

J_fixed = 1.0
gamma_values = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

sweep1_data = []
print(f"\n{'gamma_B':>8} | {'class':>12} | {'<n_XY>_S':>10} | {'<n_XY>_B':>10} | {'Re(lam)':>10} | {'Re/(-2gamma)':>10}")
print("-" * 75)

for gamma_B in gamma_values:
    H = J_fixed * 0.5 * (kron(X, X) + kron(Y, Y))
    L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, Z)])
    eigvals, eigvecs = np.linalg.eig(L)

    # Cluster by Re(lam)
    mode_data = []
    for k in range(16):
        M = eigvecs[:, k].reshape((4, 4), order='F')
        nS, nB = nxy_weights(M)
        mode_data.append((eigvals[k].real, nS, nB))

    # Sort into classes
    re_vals = sorted(set(round(d[0], 6) for d in mode_data))
    for re_class in re_vals:
        members = [d for d in mode_data if abs(d[0] - re_class) < 1e-5]
        avg_nS = np.mean([d[1] for d in members])
        avg_nB = np.mean([d[2] for d in members])
        ratio = re_class / (-2 * gamma_B) if gamma_B > 0 else 0
        count = len(members)
        if abs(re_class) < 1e-8:
            cls = "conserved"
        elif abs(re_class + 2 * gamma_B) < 1e-5:
            cls = "correlation"
        else:
            cls = "mirror"
        print(f"{gamma_B:8.3f} | {cls:>12s} | {avg_nS:10.6f} | {avg_nB:10.6f} | {re_class:10.6f} | {ratio:10.6f}")
        sweep1_data.append({
            'gamma': gamma_B, 'J': J_fixed, 'class': cls,
            'nxy_S': avg_nS, 'nxy_B': avg_nB, 'Re': re_class,
            'ratio': ratio, 'count': count
        })

# Check: does <n_XY>_B change with gamma?
mirror_nxy_B = [d['nxy_B'] for d in sweep1_data if d['class'] == 'mirror']
print(f"\nMirror class <n_XY>_B across gamma sweep: {[f'{v:.6f}' for v in mirror_nxy_B]}")
print(f"  Range: {max(mirror_nxy_B) - min(mirror_nxy_B):.2e}")
if max(mirror_nxy_B) - min(mirror_nxy_B) < 1e-6:
    print("  RESULT: <n_XY>_B is INVARIANT under gamma variation. Spatial structure independent of gamma.")
else:
    print("  RESULT: <n_XY>_B VARIES with gamma. Spatial structure depends on gamma.")

# =========================================================================
# Sweep 2: Fix gamma_B=0.1, vary J over two decades
# =========================================================================
print("\n" + "=" * 72)
print("Sweep 2: Fix gamma_B=0.1, vary J from 0.01 to 10.0")
print("Question: do <n_XY> values change with J?")
print("=" * 72)

gamma_fixed = 0.1
J_values = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

sweep2_data = []
print(f"\n{'J':>8} | {'class':>12} | {'<n_XY>_S':>10} | {'<n_XY>_B':>10} | {'Re(lam)':>10} | {'<n_XY>_B':>10}")
print("-" * 75)

for J in J_values:
    H = J * 0.5 * (kron(X, X) + kron(Y, Y))
    L = liouvillian(H, [np.sqrt(gamma_fixed) * kron(I2, Z)])
    eigvals, eigvecs = np.linalg.eig(L)

    mode_data = []
    for k in range(16):
        M = eigvecs[:, k].reshape((4, 4), order='F')
        nS, nB = nxy_weights(M)
        mode_data.append((eigvals[k].real, nS, nB))

    re_vals = sorted(set(round(d[0], 6) for d in mode_data))
    for re_class in re_vals:
        members = [d for d in mode_data if abs(d[0] - re_class) < 1e-5]
        avg_nS = np.mean([d[1] for d in members])
        avg_nB = np.mean([d[2] for d in members])
        count = len(members)
        if abs(re_class) < 1e-8:
            cls = "conserved"
        elif abs(re_class + 2 * gamma_fixed) < 1e-5:
            cls = "correlation"
        else:
            cls = "mirror"
        print(f"{J:8.3f} | {cls:>12s} | {avg_nS:10.6f} | {avg_nB:10.6f} | {re_class:10.6f} | {avg_nB:10.6f}")
        sweep2_data.append({
            'gamma': gamma_fixed, 'J': J, 'class': cls,
            'nxy_S': avg_nS, 'nxy_B': avg_nB, 'Re': re_class,
            'count': count
        })

mirror_nxy_B_J = [d['nxy_B'] for d in sweep2_data if d['class'] == 'mirror']
print(f"\nMirror class <n_XY>_B across J sweep: {[f'{v:.6f}' for v in mirror_nxy_B_J]}")
print(f"  Range: {max(mirror_nxy_B_J) - min(mirror_nxy_B_J):.2e}")

# =========================================================================
# Sweep 3: The decisive 2D grid (gamma, J)
# =========================================================================
print("\n" + "=" * 72)
print("Sweep 3: 2D grid (gamma_B, J), tracking <n_XY>_B of mirror class")
print("Question: Re(lam) = -2gamma · <n_XY>_B. Is <n_XY>_B = f(J) only?")
print("=" * 72)

gamma_grid = [0.01, 0.05, 0.1, 0.5, 1.0]
J_grid = [0.01, 0.1, 1.0, 10.0]

print(f"\n{'':>8}", end="")
for J in J_grid:
    print(f" | J={J:<6}", end="")
print()
print("-" * (10 + 11 * len(J_grid)))

grid_data = {}
for gamma_B in gamma_grid:
    print(f"gamma={gamma_B:<6.2f}", end="")
    for J in J_grid:
        H = J * 0.5 * (kron(X, X) + kron(Y, Y))
        L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, Z)])
        eigvals, eigvecs = np.linalg.eig(L)

        # Find mirror modes (Re between 0 and -2gamma, exclusive)
        mirror_nxy = []
        for k in range(16):
            re = eigvals[k].real
            if abs(re) > 1e-6 and abs(re + 2 * gamma_B) > 1e-6:
                M = eigvecs[:, k].reshape((4, 4), order='F')
                _, nB = nxy_weights(M)
                mirror_nxy.append(nB)

        avg = np.mean(mirror_nxy) if mirror_nxy else -1
        grid_data[(gamma_B, J)] = avg
        print(f" | {avg:7.4f} ", end="")
    print()

# Check: for fixed J, does <n_XY>_B vary with gamma?
print("\nColumn-wise variation (fixed J, varying gamma):")
for J in J_grid:
    col = [grid_data[(g, J)] for g in gamma_grid]
    spread = max(col) - min(col)
    print(f"  J={J:<6}: <n_XY>_B range = {spread:.2e}  {'INVARIANT' if spread < 1e-4 else 'VARIES'}")

# Check: for fixed gamma, does <n_XY>_B vary with J?
print("\nRow-wise variation (fixed gamma, varying J):")
for gamma_B in gamma_grid:
    row = [grid_data[(gamma_B, J)] for J in J_grid]
    spread = max(row) - min(row)
    print(f"  gamma={gamma_B:<6.2f}: <n_XY>_B range = {spread:.2e}  {'INVARIANT' if spread < 1e-4 else 'VARIES'}")

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 72)
print("VERDICT: Lifetime-Structure Separation at N=2")
print("=" * 72)

# Column invariance = spatial structure independent of gamma
col_invariant = all(
    max(grid_data[(g, J)] for g in gamma_grid) - min(grid_data[(g, J)] for g in gamma_grid) < 1e-4
    for J in J_grid
)

# Row invariance = spatial structure independent of J
row_invariant = all(
    max(grid_data[(g, J)] for J in J_grid) - min(grid_data[(g, J)] for J in J_grid) < 1e-4
    for g in gamma_grid
)

verdict_lines = []
if col_invariant and not row_invariant:
    verdict_lines.append("CLEAN SEPARATION at N=2.")
    verdict_lines.append("<n_XY>_B depends on J, NOT on gamma.")
    verdict_lines.append("Re(lam) = -2gamma · f(J). Lifetime = gamma × spatial structure.")
    verdict_lines.append("gamma controls the timescale. J controls the mode shape.")
    verdict_lines.append("The inside-outside correspondence is mathematically supported:")
    verdict_lines.append("  An observer can measure Re(lam) and spatial structure independently.")
    verdict_lines.append("  Lifetime reveals gamma (the outside). Shape reveals J (the inside).")
elif col_invariant and row_invariant:
    verdict_lines.append("<n_XY>_B is invariant under BOTH gamma AND J.")
    verdict_lines.append("The spatial structure is fixed by the Pauli algebra alone.")
    verdict_lines.append("Re(lam) = -2gamma · const. Pure gamma-scaling. J enters only through Im(lam).")
elif not col_invariant:
    verdict_lines.append("NO SEPARATION. <n_XY>_B depends on gamma.")
    verdict_lines.append("The spatial structure is shaped by gamma, not just J.")
    verdict_lines.append("Inside-outside correspondence fails: both parameters entangled.")

for line in verdict_lines:
    print(line)

# Save
with open(results_dir / 'probe1_separation.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 1: Lifetime-Structure Separation at N=2\n")
    f.write("=" * 72 + "\n\n")
    f.write("2D grid: <n_XY>_B of mirror class\n")
    f.write(f"{'':>8}")
    for J in J_grid:
        f.write(f" | J={J:<6}")
    f.write("\n")
    for gamma_B in gamma_grid:
        f.write(f"gamma={gamma_B:<6.2f}")
        for J in J_grid:
            f.write(f" | {grid_data[(gamma_B, J)]:7.4f} ")
        f.write("\n")
    f.write("\nVerdict:\n")
    for line in verdict_lines:
        f.write(line + "\n")

print(f"\nResults saved to {results_dir / 'probe1_separation.txt'}")
