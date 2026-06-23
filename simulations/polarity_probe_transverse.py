"""Polarity probe #2: Single-site transverse field (h_y * sum_l Y_l).

Tests whether the +/-1/2 polarity balance survives single-site Pauli terms.
pi_decompose_M only accepts bilinear / k-body Pauli tuples; single-site
contributions silently drop, so this case has never been measured.

Construction: H = J*(XX+YY+ZZ chain) + h_y * sum_l Y_l at N=3. The transverse
Y field breaks Z(x)N-mirror (zn_mirror_diagnostic detects this) and is the
canonical hardware-relevant "leakage" signal (Marrakesh h_y_eff=0.05).

Predicted by the Hermitian-conjugate argument in reflections/POLARITY_COORDINATES.md:
single-site Pauli terms still give Hermitian H, so balance should hold.
This probe tests that prediction.

Sweep h_y in {0, 0.05, 0.1, 0.5}. h_y=0 reproduces pure Heisenberg.
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
sigma = N * gamma_z

M_basis = fw._vec_to_pauli_basis_transform(N)
Pi = fw.build_pi_full(N)

# Heisenberg chain (XX+YY+ZZ chain), Hermitian
H_heis = np.zeros((d, d), dtype=complex)
for l in range(N - 1):
    H_heis = H_heis + (
        fw.site_op(N, l, 'X') @ fw.site_op(N, l + 1, 'X')
        + fw.site_op(N, l, 'Y') @ fw.site_op(N, l + 1, 'Y')
        + fw.site_op(N, l, 'Z') @ fw.site_op(N, l + 1, 'Z')
    )

# Transverse-field carrier: sum_l Y_l
H_yfield = np.zeros((d, d), dtype=complex)
for l in range(N):
    H_yfield = H_yfield + fw.site_op(N, l, 'Y')

# Z-dephasing dissipator part of L (independent of h_y)
D_part = np.zeros((d * d, d * d), dtype=complex)
for l in range(N):
    Zl = fw.site_op(N, l, 'Z')
    D_part = D_part + gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id_d, Id_d))


def build_L(h_y):
    H_total = J * H_heis + h_y * H_yfield
    L_H = -1j * (np.kron(H_total, Id_d) - np.kron(Id_d, H_total.T))
    return L_H + D_part


CASE_W = 16
TABLE_W = 84
print("=" * TABLE_W)
print(f"Polarity probe #2: H = J*Heisenberg + h_y * sum_l Y_l at N={N}, gz={gamma_z}")
print("=" * TABLE_W)
print(f"{'h_y':<{CASE_W}}  {'||M||^2':>12}  {'0%':>7}  {'+1/2%':>7}  {'-1/2%':>7}  {'asym':>12}")
print("-" * TABLE_W)

for h_y in [0.0, 0.05, 0.1, 0.5]:
    L = build_L(h_y)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    result = fw.polarity_coordinates_from_L(L_pauli, N, sigma=sigma, Pi=Pi)
    ns_M = result['norm_sq']['M']
    if ns_M > 1e-12:
        p0 = 100.0 * result['norm_sq']['M_zero'] / ns_M
        pp = 100.0 * result['norm_sq']['M_plus_half'] / ns_M
        pm = 100.0 * result['norm_sq']['M_minus_half'] / ns_M
        print(f"{h_y:<{CASE_W}.3f}  {ns_M:>12.4f}  {p0:>6.2f}%  {pp:>6.2f}%  {pm:>6.2f}%  {result['asymmetry']:+.3e}")
    else:
        print(f"{h_y:<{CASE_W}.3f}  {ns_M:>12.2e}    ---      ---      ---  {result['asymmetry']:+.3e}")

print()
print("=" * TABLE_W)
print("Reading:")
print("=" * TABLE_W)
print("  h_y = 0   -> pure Heisenberg (truly): M ~ 0 by F77/F85, no polarity content")
print("  h_y > 0   -> transverse field breaks zn_mirror; H still Hermitian.")
print("  Hermitian-conjugate hypothesis: asymmetry stays 0 (single-site is still Hermitian H).")
print("  asym != 0 -> single-site Y carries an unexpected +/-1/2 asymmetry beyond Hermitian-")
print("               conjugate symmetry; would be a hardware-relevant discovery (Marrakesh).")