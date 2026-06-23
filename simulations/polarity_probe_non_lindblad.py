"""Polarity probe #5: Non-Lindblad dissipator.

Tests whether the +/-1/2 polarity balance survives a non-Lindblad-form
dissipator. Standard Lindblad form: D[L]rho = L rho L^dagger - (1/2){L^dagger L, rho}.
A non-Lindblad dissipator drops the anti-commutator and/or uses non-(L, L^dagger)
pairs, breaking trace preservation. The Hermitian-conjugate explanation in
reflections/POLARITY_COORDINATES.md predicts non-Lindblad forms should break
the +/-1/2 balance.

Construction: H = Heisenberg chain at N=3 (truly so baseline is M ~ 0).
Two non-Lindblad dissipator forms (vec-convention: kron(L, R.conj()) represents
the action L * rho * R^dagger):
  (A) D_A[rho] = sigma+_l rho sigma-_l - rho  via kron(sp, sp.conj()).
      Standard sigma+ jump-action MINUS the {sigma- sigma+, rho}/2 anti-commutator.
      Trace-non-preserving (-rho instead of -(1/2){sigma- sigma+, rho}) and
      non-CPTP. Hermiticity-preserving (sp and sp.conj() = sm are a proper pair).
  (B) D_B[rho] = sigma+_l rho sigma+_l - rho  via kron(sp, sm.conj()).
      Both left- and right-multiplication by raising operator sigma+. Truly
      non-Lindblad: not a (L, L^dagger) pair on either side, hermiticity-
      NON-preserving (sigma+ rho sigma+ doesn't yield a Hermitian matrix).

Both forms construct on a per-site basis with rate gamma. Sweep gamma in
{0, 0.05, 0.1, 0.5}; gamma=0 is the pure Heisenberg baseline.

Reading:
  asym != 0 in either A or B -> Hermitian-conjugate / Lindblad-CPTP structure IS
              the mechanism behind the balance; non-Lindblad breaks it.
  asym = 0 in both -> balance has a DEEPER protection than Lindblad-CPTP form.
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
sigma_shift = N * gamma_z

M_basis = fw._vec_to_pauli_basis_transform(N)
Pi = fw.build_pi_full(N)

# Heisenberg chain (truly)
H_heis = np.zeros((d, d), dtype=complex)
for l in range(N - 1):
    H_heis = H_heis + J * (
        fw.site_op(N, l, 'X') @ fw.site_op(N, l + 1, 'X')
        + fw.site_op(N, l, 'Y') @ fw.site_op(N, l + 1, 'Y')
        + fw.site_op(N, l, 'Z') @ fw.site_op(N, l + 1, 'Z')
    )
L_H = -1j * (np.kron(H_heis, Id_d) - np.kron(Id_d, H_heis.T))

# Pure-Z-deph dissipator (baseline; gives Sigma_gamma = N*gamma_z F1 shift)
D_Z = np.zeros((d * d, d * d), dtype=complex)
for l in range(N):
    Zl = fw.site_op(N, l, 'Z')
    D_Z = D_Z + gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id_d, Id_d))

# Build per-site sigma+ = (X + iY)/2 in N-qubit space
sigma_plus_2 = np.array([[0, 1], [0, 0]], dtype=complex)  # raising op
sigma_minus_2 = np.array([[0, 0], [1, 0]], dtype=complex)  # lowering op


def site_op_2x2(N, l, op2x2):
    factors = [np.eye(2, dtype=complex)] * N
    factors[l] = op2x2
    out = factors[0]
    for f in factors[1:]:
        out = np.kron(out, f)
    return out


sigma_plus_l = [site_op_2x2(N, l, sigma_plus_2) for l in range(N)]
sigma_minus_l = [site_op_2x2(N, l, sigma_minus_2) for l in range(N)]


def build_L_nonlindblad_A(gamma):
    """Form A: sigma+_l rho sigma-_l - rho (jump-only, missing anti-commutator)."""
    L = L_H + D_Z
    for sp in sigma_plus_l:
        # vec(L rho R) = (R^T (x) L) vec(rho), so kron(sp, sp.conj()) acts as
        # sp * rho * sp^dagger = sigma+ rho sigma-. This is the standard jump
        # action for the c_op = sigma+; the missing anti-commutator term
        # -(1/2){sigma- sigma+, rho} makes it non-Lindblad (trace-non-preserving).
        L = L + gamma * (np.kron(sp, sp.conj()) - np.kron(Id_d, Id_d))
    return L


def build_L_nonlindblad_B(gamma):
    """Form B: sigma+_l rho sigma+_l - rho (raising on BOTH sides; not (L, L^dagger))."""
    L = L_H + D_Z
    for sp, sm in zip(sigma_plus_l, sigma_minus_l):
        # kron(sp, sm.conj()) acts as sp * rho * sm^dagger = sigma+ rho sigma+.
        # Truly non-Lindblad: not a (L, L^dagger) jump pair, hermiticity-non-
        # preserving (sigma+ rho sigma+ is not Hermitian even for Hermitian rho).
        L = L + gamma * (np.kron(sp, sm.conj()) - np.kron(Id_d, Id_d))
    return L


def run_sweep(label, build_fn, gammas):
    print("=" * TABLE_W)
    print(label)
    print("=" * TABLE_W)
    print(f"{'gamma_nonlb':<{CASE_W}}  {'||M||^2':>12}  {'0%':>7}  {'+1/2%':>7}  {'-1/2%':>7}  {'asym':>12}")
    print("-" * TABLE_W)
    for gamma in gammas:
        L = build_fn(gamma)
        L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
        result = fw.polarity_coordinates_from_L(L_pauli, N, sigma=sigma_shift, Pi=Pi)
        ns_M = result['norm_sq']['M']
        if ns_M > 1e-12:
            p0 = 100.0 * result['norm_sq']['M_zero'] / ns_M
            pp = 100.0 * result['norm_sq']['M_plus_half'] / ns_M
            pm = 100.0 * result['norm_sq']['M_minus_half'] / ns_M
            print(f"{gamma:<{CASE_W}.3f}  {ns_M:>12.4f}  {p0:>6.2f}%  {pp:>6.2f}%  {pm:>6.2f}%  {result['asymmetry']:+.3e}")
        else:
            print(f"{gamma:<{CASE_W}.3f}  {ns_M:>12.2e}    ---      ---      ---  {result['asymmetry']:+.3e}")
    print()


CASE_W = 16
TABLE_W = 84

gamma_sweep = [0.0, 0.05, 0.1, 0.5]

run_sweep(
    f"Probe #5A: H = Heisenberg + Z-deph + gamma*(sigma+_l rho sigma-_l - rho) at N={N}, gz={gamma_z}",
    build_L_nonlindblad_A,
    gamma_sweep,
)

run_sweep(
    f"Probe #5B: H = Heisenberg + Z-deph + gamma*(sigma+_l rho sigma+_l - rho) at N={N}, gz={gamma_z}",
    build_L_nonlindblad_B,
    gamma_sweep,
)

print("=" * TABLE_W)
print("Reading:")
print("=" * TABLE_W)
print("  gamma=0 -> truly Heisenberg + pure Z-deph, M ~ 0 baseline.")
print("  Form A (sigma+ rho sigma- - rho): jump-only T1 cooling missing the")
print("           anti-commutator term; trace-non-preserving but hermiticity-preserving.")
print("  Form B (sigma+ rho sigma+ - rho): raising-on-both-sides; hermiticity-")
print("           NON-preserving (the strongest break: rho doesn't stay Hermitian).")
print("  asym != 0 -> Hermitian-conjugate / Lindblad-CPTP structure IS the mechanism (real).")
print("  asym = 0 (both forms) -> balance protected even outside Lindblad CPTP (deeper).")