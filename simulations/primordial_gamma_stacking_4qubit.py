"""
Primordial Gamma: Probe 3, Does Refraction Stack at N=4?
=========================================================

4-qubit chain S - M1 - M2 - B, XX+YY coupling, gamma only on B.

Question: can we predict gamma_eff_S from layered composition of
the 3-qubit g(r) formula?

Two approaches computed and compared:

(a) DIRECT: diagonalize the full 256x256 Liouvillian, extract
    <n_XY>_B of the slowest S-coherence mode. This is exact.

(b) STACKED: treat the chain as two nested 3-qubit problems:
    Step 1: M2 sees gamma_eff_M2 = gamma_B * g(J_M1M2 / J_M2B)
    Step 2: S sees gamma_eff_S = gamma_eff_M2 * g(J_SM1 / J_M1M2)
    So: gamma_eff_S = gamma_B * g(J_M1M2/J_M2B) * g(J_SM1/J_M1M2)

If (a) = (b): refraction stacks multiplicatively.
If (a) != (b): the layered picture breaks at 4 qubits.

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from scipy.linalg import expm
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


def build_pauli_basis(n_qubits):
    ops = []
    nxy_per_site = []
    for idx in iproduct(range(4), repeat=n_qubits):
        op = paulis[idx[0]]
        for k in idx[1:]:
            op = np.kron(op, paulis[k])
        ops.append(op)
        nxy_per_site.append(tuple(n_xy_site[i] for i in idx))
    return ops, nxy_per_site


def nxy_at_site(eigvec, site, basis_ops, nxy_per_site, d):
    M = eigvec.reshape((d, d), order='F')
    coeffs = np.zeros(len(basis_ops), dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / d
    weights = np.abs(coeffs) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0
    return sum(w * nxy[site] for w, nxy in zip(weights, nxy_per_site)) / total


def s_overlap(eigvec, d, n_qubits):
    """How much does this mode contribute to S-site (site 0) coherence?"""
    M = eigvec.reshape((d, d), order='F')
    shape = [2] * (2 * n_qubits)
    Mr = M.reshape(shape)
    # Trace over all qubits except site 0
    indices = list(range(2 * n_qubits))
    for q in range(n_qubits - 1, 0, -1):
        Mr = np.trace(Mr, axis1=q, axis2=q + n_qubits)
        n_qubits -= 1
    rho_S = Mr
    return np.linalg.norm(rho_S) / np.linalg.norm(M) if np.linalg.norm(M) > 1e-14 else 0


def g_formula(r):
    """Closed-form g(r) from the analytical derivation."""
    r_cross = 1.0 / np.sqrt(2)
    if r < r_cross:
        return r**2 / (r**2 + 1)
    else:
        return 1.0 / (2 * (r**2 + 1))


# =========================================================================
# Also derive the 4-qubit analytical formula for comparison
# =========================================================================
def g_4qubit_analytical(J_SM1, J_M1M2, J_M2B):
    """Analytical <n_XY>_B for slowest S-coherence mode of a 4-site chain.

    Single-excitation Hamiltonian:
    H_1 = [[0,      J_SM1,  0,      0     ],
            [J_SM1,  0,      J_M1M2, 0     ],
            [0,      J_M1M2, 0,      J_M2B ],
            [0,      0,      J_M2B,  0     ]]

    Eigenvalues from det(H_1 - E*I) = 0:
    E^4 - (J_SM1^2 + J_M1M2^2 + J_M2B^2)*E^2 + J_SM1^2*J_M2B^2 = 0

    Quadratic in E^2:
    E^2 = (S +- sqrt(S^2 - 4*P)) / 2
    where S = J_SM1^2 + J_M1M2^2 + J_M2B^2, P = J_SM1^2 * J_M2B^2
    """
    S = J_SM1**2 + J_M1M2**2 + J_M2B**2
    P = J_SM1**2 * J_M2B**2
    disc = S**2 - 4*P
    if disc < 0:
        disc = 0

    E2_plus = (S + np.sqrt(disc)) / 2
    E2_minus = (S - np.sqrt(disc)) / 2

    # Four eigenvalues: +-sqrt(E2_plus), +-sqrt(E2_minus)
    # For each, compute eigenvector and extract |a_B|^2

    results = []
    for E2 in [E2_plus, E2_minus]:
        if E2 < 1e-20:
            continue
        E = np.sqrt(E2)
        # Eigenvector from (H-EI)v = 0:
        # -E*v0 + J_SM1*v1 = 0 -> v1 = E*v0/J_SM1
        # J_SM1*v0 - E*v1 + J_M1M2*v2 = 0
        # -> J_SM1*v0 - E^2*v0/J_SM1 + J_M1M2*v2 = 0
        # -> v2 = (E^2/J_SM1 - J_SM1)*v0 / J_M1M2 = (E^2 - J_SM1^2)/(J_SM1*J_M1M2) * v0
        # J_M1M2*v2 - E*v3 + 0 = 0 -> not quite right
        # Let me just build and diagonalize numerically
        pass

    # Numerical approach (cleaner for 4x4)
    H1 = np.array([
        [0,      J_SM1,  0,      0     ],
        [J_SM1,  0,      J_M1M2, 0     ],
        [0,      J_M1M2, 0,      J_M2B ],
        [0,      0,      J_M2B,  0     ]
    ])
    eigvals, eigvecs = np.linalg.eigh(H1)

    # For each eigenmode: |a_B|^2 = |v[3]|^2, |a_S|^2 = |v[0]|^2
    modes = []
    for k in range(4):
        a_S_sq = abs(eigvecs[0, k])**2
        a_B_sq = abs(eigvecs[3, k])**2
        modes.append({'E': eigvals[k], 'a_S_sq': a_S_sq, 'a_B_sq': a_B_sq})

    # Slowest S-coherence mode: smallest a_B_sq among modes with a_S_sq > threshold
    s_modes = [m for m in modes if m['a_S_sq'] > 0.001]
    if not s_modes:
        s_modes = modes
    s_modes.sort(key=lambda m: m['a_B_sq'])
    return s_modes[0]['a_B_sq']


# =========================================================================
# Main test
# =========================================================================
basis_4, nxy_4 = build_pauli_basis(4)
d4 = 16  # 2^4
dim_L = 256  # 4^4

gamma_B = 0.01  # good-cavity regime

print("=" * 80)
print("Probe 3: Does refraction stack at N=4?")
print("=" * 80)
print(f"gamma_B = {gamma_B} (good-cavity regime)")
print()

configs = [
    # (J_SM1, J_M1M2, J_M2B, description)
    (1.0, 1.0, 1.0, "uniform coupling"),
    (1.0, 1.0, 3.0, "strong bath coupling"),
    (1.0, 3.0, 1.0, "strong middle coupling"),
    (0.3, 1.0, 3.0, "weak S, strong B"),
    (1.0, 0.3, 1.0, "weak middle link"),
    (3.0, 1.0, 0.3, "strong S, weak B"),
    (0.1, 1.0, 10.0, "very asymmetric"),
    (1.0, 1.0, 10.0, "strong bath, uniform inner"),
    (10.0, 1.0, 1.0, "strong inner, uniform bath"),
]

print(f"{'config':>25} {'g_direct':>10} {'g_stacked':>10} {'g_4q_anal':>10} "
      f"{'dir/stack':>10} {'dir/anal':>10}")
print("-" * 80)

results = []
for J_SM1, J_M1M2, J_M2B, desc in configs:
    # (a) DIRECT: full 256x256 Liouvillian
    H = (J_SM1 * 0.5 * (kron(X, X, I2, I2) + kron(Y, Y, I2, I2))
       + J_M1M2 * 0.5 * (kron(I2, X, X, I2) + kron(I2, Y, Y, I2))
       + J_M2B * 0.5 * (kron(I2, I2, X, X) + kron(I2, I2, Y, Y)))

    L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, I2, I2, Z)])
    eigvals, eigvecs = np.linalg.eig(L)

    # Find slowest S-coherence mode
    s_modes = []
    for k in range(dim_L):
        re = eigvals[k].real
        if abs(re) < 1e-10:
            continue
        overlap = s_overlap(eigvecs[:, k], d4, 4)
        if overlap > 0.01:
            nB = nxy_at_site(eigvecs[:, k], 3, basis_4, nxy_4, d4)
            s_modes.append({'re': re, 'nxy_B': nB, 'overlap': overlap})

    s_modes.sort(key=lambda m: abs(m['re']))
    g_direct = s_modes[0]['nxy_B'] if s_modes else 0

    # (b) STACKED: two layers of g(r)
    # Layer 1: M2 sees gamma_B through J_M2B coupling
    # This is a 2-qubit problem M2-B, so g = 1/(2(r^2+1)) with r = J_M1M2/J_M2B?
    # No: the 3-qubit formula uses the ratio of INNER coupling to BATH coupling
    # For the M2-B interface: "inner" = M1-M2, "bath" = M2-B
    r1 = J_M1M2 / J_M2B
    g1 = g_formula(r1)
    gamma_eff_M2 = gamma_B * g1

    # Layer 2: S sees gamma_eff_M2 through J_SM1 coupling
    # "inner" = S-M1, "bath" = M1-M2 (now M2 carries gamma_eff_M2)
    # But wait: gamma_eff_M2 is on M2, not on M1. The stacking isn't straightforward.
    # The 3-qubit formula assumes gamma on the END qubit.
    # For the S-M1-M2 sub-chain with gamma_eff on M2: this IS a 3-qubit problem.
    r2 = J_SM1 / J_M1M2
    g2 = g_formula(r2)
    g_stacked = g1 * g2

    # (c) 4-qubit analytical (single-excitation eigenvectors)
    g_anal = g_4qubit_analytical(J_SM1, J_M1M2, J_M2B)

    ratio_ds = g_direct / g_stacked if g_stacked > 1e-15 else 0
    ratio_da = g_direct / g_anal if g_anal > 1e-15 else 0

    results.append({
        'desc': desc, 'J_SM1': J_SM1, 'J_M1M2': J_M1M2, 'J_M2B': J_M2B,
        'g_direct': g_direct, 'g_stacked': g_stacked, 'g_anal': g_anal,
        'ratio_ds': ratio_ds, 'ratio_da': ratio_da
    })

    print(f"{desc:>25} {g_direct:10.6f} {g_stacked:10.6f} {g_anal:10.6f} "
          f"{ratio_ds:10.4f} {ratio_da:10.4f}")

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 80)
print("VERDICT")
print("=" * 80)

ratios_ds = [r['ratio_ds'] for r in results if r['ratio_ds'] > 0]
ratios_da = [r['ratio_da'] for r in results if r['ratio_da'] > 0]

spread_ds = max(ratios_ds) - min(ratios_ds) if ratios_ds else 0
mean_ds = np.mean(ratios_ds) if ratios_ds else 0

spread_da = max(ratios_da) - min(ratios_da) if ratios_da else 0
mean_da = np.mean(ratios_da) if ratios_da else 0

print(f"\nDirect / Stacked ratio: mean={mean_ds:.4f}, spread={spread_ds:.4f}")
print(f"Direct / 4q-analytical ratio: mean={mean_da:.4f}, spread={spread_da:.4f}")

verdict = []
if abs(mean_ds - 1.0) < 0.1 and spread_ds < 0.2:
    verdict.append("STACKING WORKS: direct ~ stacked within 10-20%.")
    verdict.append("The refractive-index reading composes multiplicatively.")
elif abs(mean_da - 1.0) < 0.05:
    verdict.append("STACKING FAILS but 4-qubit analytical formula works.")
    verdict.append("The eigenvector approach generalizes; the layered composition does not.")
else:
    verdict.append("BOTH stacking and 4q-analytical deviate from direct.")
    verdict.append("The slowest-S-coherence-mode identification may differ at N=4.")

for line in verdict:
    print(line)

# Save
results_dir = Path("simulations/results/primordial_gamma")
with open(results_dir / 'stacking_4qubit.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 3: Refraction stacking at N=4\n")
    f.write("=" * 80 + "\n\n")
    for r in results:
        f.write(f"{r['desc']:>25}: direct={r['g_direct']:.6f}, "
                f"stacked={r['g_stacked']:.6f}, anal={r['g_anal']:.6f}, "
                f"d/s={r['ratio_ds']:.4f}, d/a={r['ratio_da']:.4f}\n")
    f.write("\nVerdict:\n")
    for line in verdict:
        f.write(line + "\n")

print(f"\nResults saved to {results_dir / 'stacking_4qubit.txt'}")
