"""Polarity probe #1: Non-Hermitian H.

Tests whether the +/-1/2 polarity balance survives a non-Hermitian Hamiltonian.
The Hermitian-conjugate symmetry argument in reflections/POLARITY_COORDINATES.md
predicts: when H is non-Hermitian, the superoperator M is no longer Hermitian
under the Frobenius adjoint, and the +i/-i Pi-eigenspace balance should break.

Construction: two sweeps to cover both M sectors.
  (A) H = J*XX-chain + i*kappa*YZ-chain : XX is Pi2Z-even truly, YZ is Pi2Z-
      even non-truly. i*kappa makes H non-Hermitian. Tests if non-Hermiticity
      breaks balance on a Pi2Z-even substrate.
  (B) H = J*XY-chain + i*kappa*XZ-chain : XY is Pi2Z-odd, XZ is Pi2Z-odd. Both
      produce non-trivial M_anti and M_plus/minus_half content; i*kappa makes
      H non-Hermitian. This is the sharp test for the +/-1/2 balance under
      non-Hermiticity.

We build L = -i(H (x) I - I (x) H^T) + Z-dephasing manually since
lindbladian_pauli_dephasing enforces Hermitian H.

Sweep kappa in {0, 0.05, 0.1, 0.5, 1.0}; kappa=0 is the Hermitian baseline
(should give asymmetry = 0 bit-exact).
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
sigma = N * gamma_z  # uniform per-site dephase, Sigma_gamma = N * gamma_z

# Pre-build the (Pauli) basis transform M and Pi once
M_basis = fw._vec_to_pauli_basis_transform(N)
Pi = fw.build_pi_full(N)


def build_chain_bilinear(N, letters):
    """Build Sum_l P_a^l P_b^{l+1} as a 2^N x 2^N matrix."""
    a, b = letters
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for l in range(N - 1):
        H = H + fw.site_op(N, l, a) @ fw.site_op(N, l + 1, b)
    return H


# Pre-build the Z-dephasing dissipator part of L (independent of kappa)
D_part = np.zeros((d * d, d * d), dtype=complex)
for l in range(N):
    Zl = fw.site_op(N, l, 'Z')
    D_part = D_part + gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id_d, Id_d))


def build_L(H_hermitian, H_anti_hermitian_carrier, kappa):
    """L for H = H_hermitian + i*kappa*H_anti_hermitian_carrier + Z-deph.

    H_anti_hermitian_carrier is itself Hermitian; multiplying by i*kappa makes
    the added term anti-Hermitian, so the total H is non-Hermitian.
    """
    H_total = H_hermitian + 1j * kappa * H_anti_hermitian_carrier
    L_H = -1j * (np.kron(H_total, Id_d) - np.kron(Id_d, H_total.T))
    return L_H + D_part


def run_sweep(label, H_hermitian, H_carrier, kappas):
    print("=" * TABLE_W)
    print(label)
    print("=" * TABLE_W)
    print(f"{'kappa':<{CASE_W}}  {'||M||^2':>12}  {'0%':>7}  {'+1/2%':>7}  {'-1/2%':>7}  {'asym':>12}")
    print("-" * TABLE_W)
    for kappa in kappas:
        L = build_L(H_hermitian, H_carrier, kappa)
        L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
        result = fw.polarity_coordinates_from_L(L_pauli, N, sigma=sigma, Pi=Pi)
        ns_M = result['norm_sq']['M']
        if ns_M > 1e-12:
            p0 = 100.0 * result['norm_sq']['M_zero'] / ns_M
            pp = 100.0 * result['norm_sq']['M_plus_half'] / ns_M
            pm = 100.0 * result['norm_sq']['M_minus_half'] / ns_M
            print(f"{kappa:<{CASE_W}.3f}  {ns_M:>12.4f}  {p0:>6.2f}%  {pp:>6.2f}%  {pm:>6.2f}%  {result['asymmetry']:+.3e}")
        else:
            print(f"{kappa:<{CASE_W}.3f}  {ns_M:>12.2e}    ---      ---      ---  {result['asymmetry']:+.3e}")
    print()


CASE_W = 16
TABLE_W = 84

kappa_sweep = [0.0, 0.05, 0.1, 0.5, 1.0]

# Sweep A: Pi2Z-even substrate (XX truly + i*kappa*YZ non-truly Pi2Z-even)
run_sweep(
    f"Probe #1A: H = J*XX-chain + i*kappa*YZ-chain at N={N}, gz={gamma_z}",
    J * build_chain_bilinear(N, ('X', 'X')),
    build_chain_bilinear(N, ('Y', 'Z')),
    kappa_sweep,
)

# Sweep B: Pi2Z-odd substrate (XY Pi2Z-odd + i*kappa*XZ Pi2Z-odd)
# This is where the +/-1/2 balance has actual non-trivial content to be broken.
run_sweep(
    f"Probe #1B: H = J*XY-chain + i*kappa*XZ-chain at N={N}, gz={gamma_z}",
    J * build_chain_bilinear(N, ('X', 'Y')),
    build_chain_bilinear(N, ('X', 'Z')),
    kappa_sweep,
)

print("=" * TABLE_W)
print("Reading:")
print("=" * TABLE_W)
print("  kappa = 0 -> Hermitian baseline; asymmetry = 0 bit-exact (confirmed in earlier wave)")
print("  Sweep A (XX + i*kappa*YZ): Pi2Z-even substrate. Non-trivial M_zero content,")
print("           tests if non-Hermiticity creates an unbalanced Pi2Z-odd projection.")
print("  Sweep B (XY + i*kappa*XZ): Pi2Z-odd substrate. M_plus/minus_half non-trivial,")
print("           SHARP test for +/-1/2 balance under non-Hermitian H.")
print("  asym != 0 in B -> Hermitian-conjugate explanation IS the mechanism (real discovery).")
print("  asym = 0 in B  -> balance has a DEEPER protection (further discovery).")
