"""Polarity probe #3: Mixed dephase letters (Z + X simultaneously).

Tests whether the +/-1/2 polarity balance survives mixed-letter dephasing
(gamma_z * D[Z] + gamma_x * D[X] applied at every site simultaneously).
pi_decompose_M is single-letter (one dephase_letter at a time); this case
requires manual L construction.

The Pi used for the polarity decomposition is canonical Z-deph Pi (= build_pi_full
with default dephase_letter='Z'). The X-deph dissipator does NOT commute with
Pi_Z's full conjugation structure (only with parts of it). The asymmetry probes
whether mixing dephase letters breaks the +/-1/2 balance with respect to Pi_Z.

Construction: H = Heisenberg chain (truly), dissipator =
gamma_z * sum_l (Z_l rho Z_l - rho) + gamma_x * sum_l (X_l rho X_l - rho).

Sweep ratio r = gamma_x / gamma_z in {0, 0.1, 0.5, 1.0, 2.0}, gamma_z = 0.05 fixed.

Hermitian-conjugate prediction: balance preserved (both Z- and X-dephasing
are hermiticity-preserving so L is still Hermitian as a superoperator).
"""

import sys
sys.path.insert(0, 'simulations')
import framework as fw

import numpy as np

N = 3
J = 1.0
gamma_z = 0.05

d = 2 ** N
Id_d = np.eye(d, dtype=complex)

M_basis = fw._vec_to_pauli_basis_transform(N)
Pi = fw.build_pi_full(N, dephase_letter='Z')

# Heisenberg chain Hamiltonian (truly, so kappa=0 baseline gives M ~ 0)
H_heis = np.zeros((d, d), dtype=complex)
for l in range(N - 1):
    H_heis = H_heis + J * (
        fw.site_op(N, l, 'X') @ fw.site_op(N, l + 1, 'X')
        + fw.site_op(N, l, 'Y') @ fw.site_op(N, l + 1, 'Y')
        + fw.site_op(N, l, 'Z') @ fw.site_op(N, l + 1, 'Z')
    )

L_H = -1j * (np.kron(H_heis, Id_d) - np.kron(Id_d, H_heis.T))

# Per-site Z and X projectors for the dephase channels
Z_per_site = [fw.site_op(N, l, 'Z') for l in range(N)]
X_per_site = [fw.site_op(N, l, 'X') for l in range(N)]


def build_L(gamma_z_val, gamma_x_val):
    L = L_H.copy()
    for Zl in Z_per_site:
        L = L + gamma_z_val * (np.kron(Zl, Zl.conj()) - np.kron(Id_d, Id_d))
    for Xl in X_per_site:
        L = L + gamma_x_val * (np.kron(Xl, Xl.conj()) - np.kron(Id_d, Id_d))
    return L


CASE_W = 16
TABLE_W = 88
print("=" * TABLE_W)
print(f"Polarity probe #3: mixed Z+X dephase. H = Heisenberg, gz={gamma_z}, sweep r = gx/gz")
print(f"Pi used for decomposition: canonical Z-deph Pi at N={N}")
print("=" * TABLE_W)
print(f"{'r=gx/gz':<{CASE_W}}  {'||M||^2':>12}  {'0%':>7}  {'+1/2%':>7}  {'-1/2%':>7}  {'asym':>12}")
print("-" * TABLE_W)

for ratio in [0.0, 0.1, 0.5, 1.0, 2.0]:
    gamma_x = ratio * gamma_z
    L = build_L(gamma_z, gamma_x)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    # Sigma_gamma for the canonical Z-deph Pi is the total Z-dephasing rate
    # only; X-dephasing has its own Pi_X but is being measured under Pi_Z here.
    # The mixed setting is not Pi-Z-palindromic in general, so we use the
    # Z-only Sigma_gamma as the "expected" shift point of the F1 palindrome.
    sigma = N * gamma_z
    result = fw.polarity_coordinates_from_L(L_pauli, N, sigma=sigma, Pi=Pi)
    ns_M = result['norm_sq']['M']
    if ns_M > 1e-12:
        p0 = 100.0 * result['norm_sq']['M_zero'] / ns_M
        pp = 100.0 * result['norm_sq']['M_plus_half'] / ns_M
        pm = 100.0 * result['norm_sq']['M_minus_half'] / ns_M
        print(f"{ratio:<{CASE_W}.3f}  {ns_M:>12.4f}  {p0:>6.2f}%  {pp:>6.2f}%  {pm:>6.2f}%  {result['asymmetry']:+.3e}")
    else:
        print(f"{ratio:<{CASE_W}.3f}  {ns_M:>12.2e}    ---      ---      ---  {result['asymmetry']:+.3e}")

print()
print("=" * TABLE_W)
print("Reading:")
print("=" * TABLE_W)
print("  r = 0   -> pure Z-dephasing, truly Heisenberg: M ~ 0 (palindrome preserved)")
print("  r > 0   -> adds X-dephasing on top, breaks the Pi_Z palindrome alignment.")
print("             The X-dissipator is Pi_X-palindromic but NOT Pi_Z-palindromic,")
print("             so M (measured w.r.t. Pi_Z, sigma_z) acquires the X-deph contribution")
print("             as a full residual.")
print("  Hermitian-conjugate prediction: balance preserved (both Z and X dephasing")
print("             preserve hermiticity of rho; total L is Hermitian as a superoperator).")
print("  asym != 0 -> mixed-letter dephase DOES break +/-1/2 balance against Pi_Z; new effect.")