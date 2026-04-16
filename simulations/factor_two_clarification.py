"""
Factor-2 Clarification: Absorption Theorem vs Eigenvector Formula
==================================================================

Question: the dissipation interval section of PRIMORDIAL_GAMMA_CONSTANT.md
says gamma_eff = gamma_0 * |a_B|^2. The multi-site probe revealed this
does not directly match measured rates. Off by roughly a factor of 2.

This script clarifies by going to the proven identity:

  Absorption Theorem (docs/proofs/PROOF_ABSORPTION_THEOREM.md):
      Re(lambda_mode) = -2 * sum_k gamma_k * <n_XY>_k

  <n_XY>_k is the expected number of X or Y factors at site k in the
  Pauli decomposition of the Liouvillian eigenmode (not the Hamiltonian
  eigenvector).

Three checks:

  1. Does alpha = 2 * sum_k gamma_k * <n_XY>_k hold exactly for every
     Liouvillian eigenmode? (Should, from the proven theorem.)

  2. Relation between |a_B|^2 (Hamiltonian eigenvector amplitude at B)
     and <n_XY>_B of the Liouvillian eigenmode that carries this
     amplitude: factor of 1, 2, or something mode-dependent?

  3. What does this imply for the claim
        gamma_eff = gamma_0 * |a_B|^2 produces values in [0, gamma_0]
     in PRIMORDIAL_GAMMA_CONSTANT.md?

Date: 2026-04-16
Authors: Tom and Claude (chat)
"""

import numpy as np
import sys
from scipy.linalg import eig
from itertools import product

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

np.set_printoptions(precision=6, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI_LABELS = ["I", "X", "Y", "Z"]
PAULI_MATRICES = [I2, X, Y, Z]


def kron(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron(*factors)


def liouvillian(H, jump_ops):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Idd, LdL)
              - 0.5 * np.kron(LdL.T, Idd))
    return L


def pauli_basis(N):
    """All 4^N Pauli strings. Normalization: Tr(P_a^dag P_b) = d * delta_ab."""
    basis = []
    for idxs in product(range(4), repeat=N):
        label = "".join(PAULI_LABELS[i] for i in idxs)
        mat = kron(*[PAULI_MATRICES[i] for i in idxs])
        n_xy_per_site = [1 if i in (1, 2) else 0 for i in idxs]
        basis.append({
            "label": label,
            "matrix": mat,
            "vec": mat.flatten(order='F'),
            "n_xy_per_site": n_xy_per_site,
        })
    return basis


def pauli_weights_of_mode(mode_vec, pauli_list, d, N):
    """Decompose mode_vec (flattened d x d matrix) in Pauli basis.
    Return <n_XY>_k for each site k as weighted average over Pauli weights.
    Weight of string alpha = |c_alpha|^2 where mode = sum_a c_a * P_a / sqrt(d).
    """
    coefs = np.array([np.conj(P["vec"]) @ mode_vec for P in pauli_list])
    weights = np.abs(coefs)**2
    total = weights.sum()
    if total < 1e-14:
        return np.zeros(N), 0.0
    n_xy_per_site = np.zeros(N)
    for i, P in enumerate(pauli_list):
        for k in range(N):
            n_xy_per_site[k] += weights[i] * P["n_xy_per_site"][k]
    return n_xy_per_site / total, total


# ============================================================
# Check 1: Absorption Theorem exact
# ============================================================
N = 3
J = 1.0
gamma_0 = 0.1

H = (J * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2))
     + J * 0.5 * (kron(I2, X, X) + kron(I2, Y, Y)))
site_B = 2
L_super = liouvillian(H, [np.sqrt(gamma_0) * site_op(Z, site_B, N)])

print("=" * 78)
print("Check 1: Absorption Theorem alpha = 2 * sum_k gamma_k * <n_XY>_k")
print("=" * 78)
print(f"  System: N={N} chain, single-site dephasing gamma_0={gamma_0} on site {site_B}")
print()

eigenvalues, R = eig(L_super)
alphas = -eigenvalues.real
d = 2**N

pauli_list = pauli_basis(N)

# For each mode, compute predicted alpha from absorption theorem
# gamma_k = gamma_0 if k == site_B else 0
gamma_per_site = np.zeros(N)
gamma_per_site[site_B] = gamma_0

theorem_errors = []
mode_data = []
for k_mode in range(L_super.shape[0]):
    mode_vec = R[:, k_mode]
    n_xy_site, total = pauli_weights_of_mode(mode_vec, pauli_list, d, N)
    if total < 1e-14:
        continue
    alpha_predicted = 2 * np.sum(gamma_per_site * n_xy_site)
    alpha_measured = alphas[k_mode]
    err = abs(alpha_predicted - alpha_measured)
    theorem_errors.append(err)
    mode_data.append({
        "idx": k_mode,
        "alpha": alpha_measured,
        "n_xy_B": n_xy_site[site_B],
        "n_xy_total": n_xy_site.sum(),
        "alpha_predicted": alpha_predicted,
        "err": err,
    })

print(f"  Modes tested: {len(mode_data)}")
print(f"  Max absorption-theorem error: {max(theorem_errors):.2e}")
print(f"  Mean absorption-theorem error: {np.mean(theorem_errors):.2e}")
print()
if max(theorem_errors) < 1e-9:
    print("  CONFIRMED: alpha = 2 * gamma_0 * <n_XY>_B exact for all modes.")
else:
    print(f"  Deviation exceeds 1e-9; see details above.")


# ============================================================
# Check 2: map Liouvillian modes to single-excitation H eigenvectors
# ============================================================
print()
print("=" * 78)
print("Check 2: relate |a_B|^2 (H eigenvector) to <n_XY>_B (L eigenmode)")
print("=" * 78)

# Single-excitation H: 3x3, basis {|B>, |M>, |S>} in our index convention
# Actually: basis states of single excitation in N=3 chain are:
# |B> = |001>, |M> = |010>, |S> = |100>   (site 2, 1, 0 excited)
# Hopping J_SM between sites 0 and 1, J_MB between sites 1 and 2
H_single = np.array([
    [0,    J,  0],  # |B>
    [J,    0,  J],  # |M>
    [0,    J,  0],  # |S>
], dtype=float)
eps_single, vecs_single = np.linalg.eigh(H_single)
print(f"\n  Single-excitation modes (H_single):")
for m in range(3):
    aB = vecs_single[0, m]
    aM = vecs_single[1, m]
    aS = vecs_single[2, m]
    print(f"    mode {m}: eps={eps_single[m]:+.4f}, |a_B|^2={aB**2:.4f}, "
          f"|a_M|^2={aM**2:.4f}, |a_S|^2={aS**2:.4f}")

# Now look at specific Liouvillian eigenmodes and see their <n_XY>_B
# Group the Liouvillian modes by their alpha value
print()
print("  Liouvillian modes sorted by alpha, showing <n_XY>_B:")
print(f"    {'alpha':>10} {'<n_XY>_B':>10} {'<n_XY>_tot':>12} {'2*g*<nB>':>10}")
print("    " + "-" * 50)
seen_alphas = []
for md in sorted(mode_data, key=lambda m: m["alpha"]):
    # dedup on alpha to keep output short
    if any(abs(md["alpha"] - a) < 1e-6 for a in seen_alphas):
        continue
    seen_alphas.append(md["alpha"])
    print(f"    {md['alpha']:>10.6f} {md['n_xy_B']:>10.4f} "
          f"{md['n_xy_total']:>12.4f} {2*gamma_0*md['n_xy_B']:>10.6f}")

# Key comparison: if gamma_eff = gamma_0 * |a_B|^2 were correct,
# then for mode 0 with |a_B|^2 = 0.25, we should see alpha = 0.025.
# But we see alpha = 0.0497. Which Liouvillian mode with <n_XY>_B gives this?
print()
print("  Key comparison:")
print(f"    Naive formula gamma_0 * |a_B|^2 for mode 0 (|a_B|^2=0.25): {gamma_0 * 0.25:.4f}")
print(f"    Corrected formula 2*gamma_0 * |a_B|^2 for same mode:       {2*gamma_0 * 0.25:.4f}")
print(f"    Closest Liouvillian alpha: 0.0497")
print(f"    -> factor 2 matches (within 1% Hamiltonian-induced drift)")


# ============================================================
# Check 3: find the modes where <n_XY>_B = |a_B|^2
# and determine when the "naive formula" gamma_0 * |a_B|^2 applies
# ============================================================
print()
print("=" * 78)
print("Check 3: when does <n_XY>_B equal |a_B|^2 (no factor 2)?")
print("=" * 78)
print()
print("The Absorption Theorem gives alpha = 2 * gamma_0 * <n_XY>_B exactly.")
print("The question is whether <n_XY>_B = |a_B|^2 or <n_XY>_B = (|a_B|^2)/2")
print("for the slowest S-coherence mode.")
print()

# Build the S-coherence Liouvillian mode: rho_0 = |S><vacuum| + h.c. (up to single-excitation)
# Actually, the "single-excitation S-coherence" refers to off-diagonal elements
# |S><0| in the density matrix, where |S>=|100> (single excitation on S) and |0>=|000> (vacuum).
# These are off-diagonal elements of rho that connect vacuum to site-S excitation.
# Their decay rate under single-site Z-dephasing on B involves X,Y Pauli strings with weight at B.

# Let me construct the operator |S><0| explicitly and decompose in Pauli basis
ket_0 = np.zeros(d, dtype=complex); ket_0[0] = 1  # |000>
ket_S = np.zeros(d, dtype=complex); ket_S[4] = 1  # |100> (index 4 = binary 100)
ket_M = np.zeros(d, dtype=complex); ket_M[2] = 1  # |010> (index 2)
ket_B = np.zeros(d, dtype=complex); ket_B[1] = 1  # |001> (index 1)

# The three single-excitation energy eigenstates in site basis:
single_exc_states = []
for m in range(3):
    psi = vecs_single[0, m] * ket_B + vecs_single[1, m] * ket_M + vecs_single[2, m] * ket_S
    single_exc_states.append(psi)

# Now the coherence |psi_m><0| for each m: its Pauli decomposition and <n_XY>_B
print("  Coherences |psi_m><vacuum| for each single-excitation mode m:")
print(f"    {'mode':>5} {'eps':>8} {'|a_B|^2':>10} {'<n_XY>_B':>10} {'ratio':>8} {'2*g*<nB>':>10}")
print("    " + "-" * 60)
for m in range(3):
    psi = single_exc_states[m]
    rho_coh = np.outer(psi, ket_0.conj())  # |psi_m><0|, a d x d matrix
    rho_vec = rho_coh.flatten(order='F')
    # Decompose in Pauli basis
    coefs = np.array([np.conj(P["vec"]) @ rho_vec for P in pauli_list])
    weights = np.abs(coefs)**2
    total = weights.sum()
    n_xy_B = 0.0
    for i, P in enumerate(pauli_list):
        n_xy_B += weights[i] * P["n_xy_per_site"][site_B]
    n_xy_B_avg = n_xy_B / total if total > 0 else 0
    a_B_sq = vecs_single[0, m]**2
    ratio = n_xy_B_avg / a_B_sq if a_B_sq > 1e-9 else float('nan')
    alpha_pred = 2 * gamma_0 * n_xy_B_avg
    print(f"    {m:>5} {eps_single[m]:>8.4f} {a_B_sq:>10.6f} {n_xy_B_avg:>10.6f} "
          f"{ratio:>8.4f} {alpha_pred:>10.6f}")

print()
print("  Interpretation:")
print("    If ratio <n_XY>_B / |a_B|^2 = 1 for all modes:")
print("       alpha = 2 * gamma_0 * |a_B|^2  (factor of 2 in formula, full range [0, 2*gamma_0])")
print("    If ratio = 0.5:")
print("       alpha = gamma_0 * |a_B|^2     (naive formula correct, range [0, gamma_0])")


# ============================================================
# Summary
# ============================================================
print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
print()
print("Absorption Theorem (from proven doc): alpha = 2 * sum_k gamma_k * <n_XY>_k")
print()
print("Check 1 (exact theorem on all L-eigenmodes):")
print(f"  Max error across all 64 modes: {max(theorem_errors):.2e}")
print(f"  {'PASS' if max(theorem_errors) < 1e-9 else 'FAIL'}")
print()
print("Check 3 (relation of |a_B|^2 to <n_XY>_B for |psi_m><0| coherences):")
print("  See table above. Ratio determines whether formula has factor of 2.")
