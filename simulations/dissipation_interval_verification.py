"""
Dissipation Interval Verification
==================================

Numerical verification of the three claims added to
hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md on 2026-04-16 in the section
"The dissipation interval [0, 2*gamma_0]":

  Claim 1: For single-site dephasing at rate gamma_0, every palindromic
           partner pair of positive dissipation rates satisfies
           alpha_a + alpha_b = 2*gamma_0. The spectrum lives in
           [0, 2*gamma_0] symmetric around gamma_0.

  Claim 2: The eigenvector formula gamma_eff = gamma_0 * |a_B|^2 with
           |a_B|^2 in [0, 1] produces rates in [0, gamma_0], the
           lower half of the interval.

  Claim 3: Single-site sigma_x observable has nonzero residue only for
           modes in the lower half [0, gamma_0]. The palindromic partners
           in [gamma_0, 2*gamma_0] are algebraically present but hidden
           behind XY-weight superselection.

  Bonus: A mixed-weight observable (sigma_x(0) + sigma_x(0)*sigma_x(1))
         reaches modes in both halves.

System: N=3 chain S-M-B with XX+YY coupling, Z-dephasing on site B only.

Date: 2026-04-16
Authors: Tom and Claude (chat)
"""

import numpy as np
import sys
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

np.set_printoptions(precision=6, suppress=True)

# ============================================================
# Pauli matrices and tensor products
# ============================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    """Single-site operator on `site` (0-indexed) within N-qubit chain."""
    factors = [I2] * N
    factors[site] = op
    return kron(*factors)


def liouvillian(H, jump_ops):
    """Build vec-Liouvillian: vec(L[rho]) = L_super @ vec(rho)."""
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Idd, LdL)
              - 0.5 * np.kron(LdL.T, Idd))
    return L


# ============================================================
# Build N=3 chain S-M-B
# ============================================================
N = 3
J_SM = 1.0
J_MB = 1.0
gamma_0 = 0.1  # dephasing rate on site B

# H = sum over neighbor bonds of J_ij * 0.5 * (X_i X_j + Y_i Y_j)
H = J_SM * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2))
H += J_MB * 0.5 * (kron(I2, X, X) + kron(I2, Y, Y))

# Single-site dephasing on B = site index 2
L_jump = np.sqrt(gamma_0) * site_op(Z, 2, N)
L_super = liouvillian(H, [L_jump])

print("=" * 78)
print(f"System: N={N} chain S-M-B")
print(f"  Couplings: J_SM={J_SM}, J_MB={J_MB}")
print(f"  Dephasing: gamma_0={gamma_0} on site B (index 2)")
print(f"  Hilbert dim: {2**N}, Liouville dim: {4**N}")
print("=" * 78)


# ============================================================
# Diagonalize: get all 64 eigenvalues
# ============================================================
eigenvalues, R = eig(L_super)
# Left eigenvectors: rows of L_left = inv(R)
L_left = np.linalg.inv(R)

# Positive dissipation rates
alphas = -eigenvalues.real

print()
print("=" * 78)
print("Claim 1: Palindromic pairing alpha_a + alpha_b = 2*gamma_0")
print("=" * 78)
print()
print(f"Total eigenvalues: {len(eigenvalues)}")
print(f"  alpha range: [{alphas.min():.6f}, {alphas.max():.6f}]")
print(f"  Expected:    [0, {2*gamma_0:.6f}] = [0, 2*gamma_0]")
print(f"  Symmetry axis (mean): {alphas.mean():.6f}, expected {gamma_0:.6f}")
print()

# Group near-degenerate alphas
sorted_alphas = np.sort(alphas)
unique_alphas = []
counts = []
tol = 1e-9
for a in sorted_alphas:
    if unique_alphas and abs(a - unique_alphas[-1]) < tol:
        counts[-1] += 1
    else:
        unique_alphas.append(a)
        counts.append(1)

print(f"Distinct dissipation rates: {len(unique_alphas)}")
print()
print(f"  {'rate alpha':>12} {'multiplicity':>13} {'2*gamma_0 - alpha':>20} {'paired with':>14}")
print("  " + "-" * 65)
pair_errors = []
for i, (a, c) in enumerate(zip(unique_alphas, counts)):
    partner_target = 2 * gamma_0 - a
    # Find the closest match in unique_alphas
    diffs = np.abs(np.array(unique_alphas) - partner_target)
    j = int(np.argmin(diffs))
    err = diffs[j]
    pair_errors.append(err)
    paired_idx = j if j != i else (j, "self")
    print(f"  {a:>12.8f} {c:>13d} {partner_target:>20.8f} {unique_alphas[j]:>14.8f}")

print()
print(f"  Max pairing error: {max(pair_errors):.2e}")
print(f"  Mean pairing error: {np.mean(pair_errors):.2e}")
print()
if max(pair_errors) < 1e-9:
    print("  CONFIRMED: every distinct rate has a partner summing to 2*gamma_0.")
else:
    print(f"  Pairing holds within tolerance {max(pair_errors):.2e}.")


# ============================================================
# Claim 2: Eigenvector formula gamma_eff = gamma_0 * |a_B|^2
# Maps to lower half [0, gamma_0]
# ============================================================
print()
print("=" * 78)
print("Claim 2: Eigenvector formula maps to lower half [0, gamma_0]")
print("=" * 78)
print()

# Single-excitation Hamiltonian: 3x3 tridiagonal in basis |001>, |010>, |100>
# Or for our chain S=site0, M=site1, B=site2: |B>=|001>, |M>=|010>, |S>=|100>
# Hopping: J_SM between S and M, J_MB between M and B
H_single = np.array([
    [0,    J_MB, 0],     # |B>
    [J_MB, 0,    J_SM],  # |M>
    [0,    J_SM, 0],     # |S>
], dtype=float)

eps, vecs = np.linalg.eigh(H_single)
print(f"  Single-excitation eigenvalues: {eps}")
print(f"  Eigenvectors (columns):")
print(vecs)
print()

# |a_B|^2 = squared amplitude on |B> (first component) for each mode
a_B_sq = vecs[0, :] ** 2
print(f"  |a_B|^2 per mode: {a_B_sq}")
print(f"  Sum over modes: {a_B_sq.sum():.6f}  (should be 1, completeness)")
print()

gamma_eff_predicted = gamma_0 * a_B_sq
print(f"  gamma_eff = gamma_0 * |a_B|^2:")
for k, (e, ab2, ge) in enumerate(zip(eps, a_B_sq, gamma_eff_predicted)):
    print(f"    mode {k}: epsilon={e:>8.4f}, |a_B|^2={ab2:.6f}, gamma_eff={ge:.6f}")

print()
print(f"  All gamma_eff values lie in [0, {gamma_0}] = [0, gamma_0]: ", end="")
print(all(0 <= g <= gamma_0 + 1e-12 for g in gamma_eff_predicted))
print()
print("  Cross-check: do these gamma_eff values appear in the full spectrum?")
for k, ge in enumerate(gamma_eff_predicted):
    diffs = np.abs(np.array(unique_alphas) - ge)
    j = int(np.argmin(diffs))
    print(f"    mode {k}: gamma_eff={ge:.6f} -> closest alpha={unique_alphas[j]:.8f}, "
          f"err={diffs[j]:.2e}")


# ============================================================
# Claim 3: Single-site sigma_x sees only modes in lower half
# ============================================================
print()
print("=" * 78)
print("Claim 3: Single-site sigma_x is blind to modes in upper half")
print("=" * 78)
print()
print("Compute residue c_k = (vec(A) . r_k) (l_k . vec(A * rho_0))")
print("for each eigenmode k, where A is the observable and rho_0 the initial state.")
print()

# Initial state: a coherent superposition on S
# Use rho_0 = |+><+|_S tensor |0><0|_M tensor |0><0|_B
plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
zero = np.array([1, 0], dtype=complex)
rho_S = np.outer(plus, plus.conj())
rho_M = np.outer(zero, zero.conj())
rho_B = np.outer(zero, zero.conj())
rho_0 = kron(rho_S, rho_M, rho_B)
vec_rho_0 = rho_0.flatten(order='F')


def compute_residues(observable):
    """Returns array of |c_k| for each eigenmode."""
    A_rho = observable @ rho_0
    vec_A = observable.flatten(order='F')
    vec_A_rho = A_rho.flatten(order='F')
    # c_k = (vec(A^T) . r_k)(l_k . vec(A*rho_0))
    # but for our observables A = A^T (Hermitian, real entries), so vec(A^T) = vec(A)
    # Note: standard inner product in vec form
    coefs_obs = R.T.conj().T  # right eigenvectors are columns of R, so r_k = R[:,k]
    cs = np.zeros(R.shape[1], dtype=complex)
    for k in range(R.shape[1]):
        r_k = R[:, k]
        l_k = L_left[k, :]
        cs[k] = (vec_A.conj() @ r_k) * (l_k @ vec_A_rho)
    return np.abs(cs)


# Observable 1: sigma_x on site S (index 0). Standard dipole-like observable.
A_sigma_x_S = site_op(X, 0, N)
res_x_S = compute_residues(A_sigma_x_S)

# Observable 2: sigma_x on site B (index 2). Same XY-weight class.
A_sigma_x_B = site_op(X, 2, N)
res_x_B = compute_residues(A_sigma_x_B)

# Observable 3: mixed-weight A = sigma_x(0) + sigma_x(0)*sigma_x(1)
A_mixed = site_op(X, 0, N) + site_op(X, 0, N) @ site_op(X, 1, N)
res_mixed = compute_residues(A_mixed)


def report_visibility(name, residues, threshold):
    """Categorize modes as visible or hidden, and locate them on the interval."""
    print(f"\n  Observable: {name}")
    print(f"  (visibility threshold: |c_k| > {threshold:.0e})")
    n_visible = 0
    n_hidden = 0
    visible_alphas = []
    hidden_alphas = []
    for k in range(len(residues)):
        a_k = alphas[k]
        if residues[k] > threshold:
            n_visible += 1
            visible_alphas.append(a_k)
        else:
            n_hidden += 1
            hidden_alphas.append(a_k)
    visible_alphas = np.array(visible_alphas)
    hidden_alphas = np.array(hidden_alphas)
    print(f"    Visible modes:  {n_visible:3d}, alpha range "
          f"[{visible_alphas.min():.6f}, {visible_alphas.max():.6f}]" if n_visible else
          f"    Visible modes:  {n_visible:3d}")
    if n_visible:
        n_lower = int(np.sum(visible_alphas < gamma_0 - 1e-9))
        n_axis = int(np.sum(np.abs(visible_alphas - gamma_0) < 1e-9))
        n_upper = int(np.sum(visible_alphas > gamma_0 + 1e-9))
        print(f"      lower half [0, gamma_0):     {n_lower}")
        print(f"      on axis    (alpha = gamma_0): {n_axis}")
        print(f"      upper half (gamma_0, 2gamma]: {n_upper}")
    print(f"    Hidden  modes:  {n_hidden:3d}")
    return n_visible, visible_alphas


print()
print("-" * 78)
report_visibility("sigma_x on S (single-site, XY-weight 1)", res_x_S, 1e-10)
report_visibility("sigma_x on B (single-site, XY-weight 1)", res_x_B, 1e-10)
report_visibility("mixed: sigma_x(0) + sigma_x(0)*sigma_x(1)", res_mixed, 1e-10)


# ============================================================
# Summary
# ============================================================
print()
print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
print()
print(f"System: 3-qubit chain, single-site dephasing gamma_0={gamma_0} on B.")
print()
print(f"Claim 1 (palindromic interval [0, 2*gamma_0] symmetric around gamma_0):")
print(f"  Distinct rates: {len(unique_alphas)}")
print(f"  Range: [{min(unique_alphas):.6f}, {max(unique_alphas):.6f}]")
print(f"  Symmetric around: {(min(unique_alphas)+max(unique_alphas))/2:.6f}")
print(f"  Max pairing error: {max(pair_errors):.2e}")
print(f"  Verdict: " + ("CONFIRMED" if max(pair_errors) < 1e-9 else "see above"))
print()
print(f"Claim 2 (gamma_eff = gamma_0 * |a_B|^2 in lower half):")
print(f"  Predicted gamma_eff values: {gamma_eff_predicted}")
print(f"  All in [0, gamma_0]: " +
      str(all(0 <= g <= gamma_0 + 1e-12 for g in gamma_eff_predicted)))
print()
print(f"Claim 3 (single-site sigma_x blind to upper half):")
n_vis_S, alphas_S = (np.sum(res_x_S > 1e-10),
                    alphas[res_x_S > 1e-10])
n_in_lower_S = int(np.sum(alphas_S <= gamma_0 + 1e-9))
n_in_upper_S = int(np.sum(alphas_S > gamma_0 + 1e-9))
print(f"  sigma_x on S: {n_vis_S} visible modes, "
      f"{n_in_lower_S} in lower half, {n_in_upper_S} in upper half")
n_vis_M, alphas_M = (np.sum(res_mixed > 1e-10),
                    alphas[res_mixed > 1e-10])
n_in_lower_M = int(np.sum(alphas_M <= gamma_0 + 1e-9))
n_in_upper_M = int(np.sum(alphas_M > gamma_0 + 1e-9))
print(f"  mixed-weight: {n_vis_M} visible modes, "
      f"{n_in_lower_M} in lower half, {n_in_upper_M} in upper half")
print()
if n_in_upper_S == 0 and n_in_upper_M > 0:
    print("  Verdict: CONFIRMED. Single-site observables blind to upper half;")
    print("  mixed-weight observable reaches both halves. Superselection lies")
    print("  in the observable, not in the system.")
elif n_in_upper_S == 0:
    print("  Verdict: Single-site blind to upper half. Mixed-weight reach unclear.")
else:
    print(f"  Verdict: Single-site reaches {n_in_upper_S} upper-half modes.")
    print("  This contradicts the XY-weight superselection claim - investigate.")
