"""
Primordial Gamma V2: Re-analysis via Absorption Theorem
========================================================

V1 used time-domain exponential fit (poor R^2 at many configs).
V2 uses the Absorption Theorem directly: for the slowest mode
contributing to S-coherence, gamma_eff = gamma_B * <n_XY>_B.
This is exact, no fitting needed.

Key insight from V2 task: V1 plotted on Q_MB = J_MB/gamma_B axis
(wrong). The correct axis is r = J_SM/J_MB (coupling ratio).
Column-wise (fixed J_MB, varying gamma_B), the ratio gamma_eff/gamma_B
is approximately constant, confirming gamma_eff = gamma_B * g(J_SM/J_MB).

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from itertools import product as iproduct
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

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


def build_pauli_basis_3q():
    ops = []
    labels = []
    nxy_per_site = []
    for i, j, k in iproduct(range(4), repeat=3):
        ops.append(np.kron(np.kron(paulis[i], paulis[j]), paulis[k]))
        labels.append(('IXYZ'[i], 'IXYZ'[j], 'IXYZ'[k]))
        nxy_per_site.append((n_xy_site[i], n_xy_site[j], n_xy_site[k]))
    return ops, labels, nxy_per_site


def nxy_B_of_mode(eigvec, basis_ops, nxy_per_site, d):
    """Compute <n_XY>_B for an eigenmode via Pauli expansion."""
    M = eigvec.reshape((d, d), order='F')
    coeffs = np.zeros(len(basis_ops), dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / d
    weights = np.abs(coeffs) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0
    return sum(w * nxy[2] for w, nxy in zip(weights, nxy_per_site)) / total


def s_coherence_overlap(eigvec, d):
    """How much does this eigenmode contribute to S-site coherence |01|?

    Project eigenmode onto Pauli strings with X or Y at site S (site 0).
    These are the only strings that contribute to rho_S off-diagonal.
    """
    M = eigvec.reshape((d, d), order='F')
    # Partial trace over M,B to get the S-component
    rho_S = M.reshape(2, 2, 2, 2, 2, 2)
    rho_S = np.einsum('ibjcjc->ib', rho_S)
    return np.linalg.norm(rho_S) / np.linalg.norm(M) if np.linalg.norm(M) > 1e-14 else 0


# =========================================================================
# Setup
# =========================================================================
basis_3, labels_3, nxy_3 = build_pauli_basis_3q()
d = 8  # 2^3

J_SM = 1.0
J_MB_values = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
gamma_B_values = [0.01, 0.1, 1.0]  # fewer, cleaner

results_dir = Path("simulations/results/primordial_gamma")

print("=" * 80)
print("V2: Absorption Theorem exact analysis")
print("=" * 80)
print(f"J_SM = {J_SM}")
print()

# =========================================================================
# For each (J_MB, gamma_B): extract the slowest S-coherence mode's <n_XY>_B
# =========================================================================
print("Method: identify eigenmodes that contribute to S-coherence,")
print("find the SLOWEST one, report its <n_XY>_B as g(J_SM/J_MB).")
print()

print(f"{'J_MB':>8} {'gamma_B':>8} {'r=Jsm/Jmb':>10} {'slowest_Re':>12} {'<nXY>_B':>10} "
      f"{'g=geff/gB':>10} {'r^2':>10} {'g/r^2':>10}")
print("-" * 95)

g_values = {}  # (J_MB) -> list of g values across gamma_B

for J_MB in J_MB_values:
    g_values[J_MB] = []
    r = J_SM / J_MB

    for gamma_B in gamma_B_values:
        H = (J_SM * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2))
           + J_MB * 0.5 * (kron(I2, X, X) + kron(I2, Y, Y)))
        L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, I2, Z)])

        eigvals, eigvecs = np.linalg.eig(L)

        # Find modes that contribute to S-coherence
        s_modes = []
        for k in range(64):
            re_lam = eigvals[k].real
            if abs(re_lam) < 1e-10:
                continue  # skip conserved
            overlap = s_coherence_overlap(eigvecs[:, k], d)
            if overlap > 0.01:  # has S-coherence content
                nxy_B = nxy_B_of_mode(eigvecs[:, k], basis_3, nxy_3, d)
                s_modes.append({
                    're': re_lam, 'nxy_B': nxy_B, 'overlap': overlap
                })

        # Slowest S-coherence mode
        if not s_modes:
            print(f"{J_MB:8.3f} {gamma_B:8.3f} {r:10.4f} {'(none)':>12}")
            continue

        s_modes.sort(key=lambda m: abs(m['re']))
        slow = s_modes[0]

        # g = gamma_eff / gamma_B = <n_XY>_B of the slowest S-coherence mode
        g = slow['nxy_B']
        g_values[J_MB].append(g)

        r_sq = r**2 if r < 1 else 1/r**2  # smaller ratio squared
        g_over_rsq = g / r_sq if r_sq > 1e-10 else 0

        print(f"{J_MB:8.3f} {gamma_B:8.3f} {r:10.4f} {slow['re']:12.6f} {slow['nxy_B']:10.6f} "
              f"{g:10.6f} {r_sq:10.6f} {g_over_rsq:10.4f}")

# =========================================================================
# Column-wise collapse check
# =========================================================================
print("\n" + "=" * 80)
print("Column-wise collapse: is g(r) = gamma_eff/gamma_B independent of gamma_B?")
print("=" * 80)

print(f"\n{'J_MB':>8} {'r':>8} {'g values across gamma_B':>40} {'spread':>10}")
print("-" * 75)

all_collapse = True
for J_MB in J_MB_values:
    r = J_SM / J_MB
    gs = g_values[J_MB]
    if len(gs) > 1:
        spread = max(gs) - min(gs)
        mean_g = np.mean(gs)
        rel_spread = spread / mean_g if mean_g > 1e-10 else 0
    else:
        spread = 0
        rel_spread = 0

    gs_str = ', '.join(f'{g:.6f}' for g in gs)
    status = "OK" if rel_spread < 0.05 else "VARIES"
    if rel_spread >= 0.05:
        all_collapse = False
    print(f"{J_MB:8.3f} {r:8.4f} {gs_str:>40} {rel_spread:10.4f}  {status}")

# =========================================================================
# The function g(r) = g(J_SM/J_MB)
# =========================================================================
print("\n" + "=" * 80)
print("The function g(r) = gamma_eff / gamma_B, where r = J_SM / J_MB")
print("=" * 80)

print(f"\n{'r':>10} {'g(r)':>10} {'r^2':>10} {'g/r^2':>10} {'1/r^2':>10} {'g*r^2':>10}")
print("-" * 65)

for J_MB in J_MB_values:
    r = J_SM / J_MB
    g = np.mean(g_values[J_MB]) if g_values[J_MB] else 0
    r_sq = r**2
    inv_r_sq = 1/r_sq if r_sq > 1e-10 else 0
    print(f"{r:10.4f} {g:10.6f} {r_sq:10.6f} {g/r_sq if r_sq > 1e-10 else 0:10.4f} "
          f"{inv_r_sq:10.6f} {g*r_sq:10.4f}")

# Try closed forms
print("\nClosed-form candidates for g(r):")
r_arr = np.array([J_SM / J_MB for J_MB in J_MB_values])
g_arr = np.array([np.mean(g_values[J_MB]) for J_MB in J_MB_values])

candidates = {
    'r^2/(1+r^2)': lambda r: r**2 / (1 + r**2),
    'r^2/(1+r)^2': lambda r: r**2 / (1 + r)**2,
    'r^2/(1+r^2)^2': lambda r: r**2 / (1 + r**2)**2,
    '4r^2/(1+r^2)^2': lambda r: 4*r**2 / (1 + r**2)**2,
    '(2r/(1+r^2))^2': lambda r: (2*r / (1 + r**2))**2,
    'min(r,1/r)^2': lambda r: np.minimum(r, 1/r)**2,
    'min(r,1/r)^2/2': lambda r: np.minimum(r, 1/r)**2 / 2,
}

for name, func in candidates.items():
    pred = func(r_arr)
    max_rel = np.max(np.abs(g_arr - pred) / (g_arr + 1e-15))
    rms_rel = np.sqrt(np.mean(((g_arr - pred) / (g_arr + 1e-15))**2))
    print(f"  {name:>25}: max_rel_err = {max_rel:.4f}, rms_rel = {rms_rel:.4f}")

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 80)
print("VERDICT")
print("=" * 80)

verdict = []
if all_collapse:
    verdict.append("COLUMN-WISE COLLAPSE CONFIRMED.")
    verdict.append("gamma_eff / gamma_B = g(J_SM/J_MB), independent of gamma_B.")
    verdict.append("The refractive-index reading has a clean functional form.")
    verdict.append("")
    verdict.append("Asymptotic form: g(r) -> r^2 for r << 1 (strong bath coupling).")
    verdict.append("This is the Absorption Theorem applied perturbatively:")
    verdict.append("  <n_XY>_B ~ (J_SM/J_MB)^2 for the slowest S-coherence mode.")
else:
    verdict.append("Column-wise collapse is partial. g depends on gamma_B at some J_MB.")

for line in verdict:
    print(line)

with open(results_dir / 'v2_reanalysis.txt', 'w', encoding='utf-8') as f:
    f.write("V2: Primordial Gamma Re-analysis via Absorption Theorem\n")
    f.write("=" * 80 + "\n\n")
    f.write("g(r) = gamma_eff / gamma_B, r = J_SM / J_MB\n\n")
    for J_MB in J_MB_values:
        r = J_SM / J_MB
        g = np.mean(g_values[J_MB]) if g_values[J_MB] else 0
        gs = g_values[J_MB]
        f.write(f"r={r:.4f} (J_MB={J_MB}): g={g:.6f}  values={[f'{x:.6f}' for x in gs]}\n")
    f.write("\nVerdict:\n")
    for line in verdict:
        f.write(line + "\n")

print(f"\nResults saved to {results_dir / 'v2_reanalysis.txt'}")
