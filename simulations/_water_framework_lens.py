#!/usr/bin/env python3
"""Water/proton chain through the framework's bottom-up lens.

In `simulations/water/`, the proton in O–H...O is modeled as a true
two-state qubit: |L⟩ (proton on donor side), |R⟩ (proton on acceptor
side). Tunneling through the H-bond barrier acts as a transverse field.
Dephasing acts in the |L⟩/|R⟩ basis. This is a real qubit substrate
at d = 2, so framework.py's primitives apply DIRECTLY without
substrate translation. The framework's bit_a/bit_b decomposition,
Π conjugation, and palindrome residual all act on the proton
operator space the same way they act on a generic qubit operator
space.

What the existing water scripts test:
- `hydrogen_bond_qubit.py`: builds the proton Liouvillian with
  H = −J·(X_0 + X_1) + K·Z_0·Z_1 (transverse-field Ising) for two
  protons, dephasing γ·σ_z on each. Reports "Palindrome: EXACT" via
  eigenvalue pairing.
- `proton_water_chain.py`: builds two model Hamiltonians for N = 1..5
  protons: (a) Heisenberg J·(XX+YY+ZZ), (b) transverse-field Ising
  −J·X + K·ZZ. Compares them.

What today's framework adds:
The eigenvalue-pairing test passes for many Hamiltonians that fail
the OPERATOR equation Π·L·Π⁻¹ + L + 2Σγ·I = 0. The "EXACT" verdict
in `hydrogen_bond_qubit.py` is a SPECTRAL claim. With framework.py
we can ask the operator-level question: are the proton Hamiltonians
truly palindromic, soft-broken (spectrum pairs but operator doesn't),
or hard-broken?

This script tests both Hamiltonians via framework's `palindrome_residual`.

Run: `python simulations/_water_framework_lens.py`
"""
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


# ════════════════════════════════════════════════════════════════════
#  Build the two water Hamiltonians using framework's primitives
# ════════════════════════════════════════════════════════════════════

def H_heisenberg_chain(N, J=1.0):
    """H = J · Σ_bond (XX + YY + ZZ)  — Stage 6 admissible."""
    return fw.ur_heisenberg(N, J=J)


def H_transverse_field_ising(N, J_tunnel=1.0, K_coupling=0.2):
    """H = −J · Σ_l X_l  +  K · Σ_bond Z_l Z_{l+1}.

    This is the actual model in `hydrogen_bond_qubit.py` and the (b)
    branch of `proton_water_chain.py`. The X-term is one-body; the
    ZZ-term is two-body both-parity-even (one of the four Stage-6
    admissible bilinears modulo identity). The framework's selection
    rule is for two-body bilinears only — one-body terms get tested
    separately.
    """
    bonds = [(i, i + 1) for i in range(N - 1)]
    H = fw._build_bilinear(N, bonds, [('Z', 'Z', K_coupling)])
    for l in range(N):
        H = H + (-J_tunnel) * fw.site_op(N, l, 'X')
    return H


# ════════════════════════════════════════════════════════════════════
#  Spectral and operator tests via framework primitives
# ════════════════════════════════════════════════════════════════════

def spectral_pair_max_error(L, Sigma_gamma):
    """V-Effect-style: max |λ_i + λ_j + 2Σγ| over best pairing."""
    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * Sigma_gamma
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        best_j = int(np.argmin(dists))
        if best_j != i:
            used[i] = True
            used[best_j] = True
        else:
            used[i] = True
        max_err = max(max_err, float(dists[best_j]))
    return max_err


def classify(L, Sigma_gamma, N):
    """Apply the trichotomy: truly / soft / hard."""
    M_pauli = fw.palindrome_residual(L, Sigma_gamma, N)
    op_norm = float(np.linalg.norm(M_pauli))
    spec_err = spectral_pair_max_error(L, Sigma_gamma)
    spec_ok = spec_err < 1e-6
    op_ok = op_norm < 1e-10
    if op_ok:
        verdict = "truly_unbroken"
    elif spec_ok:
        verdict = "soft_broken"
    else:
        verdict = "hard_broken"
    return verdict, op_norm, spec_err


# ════════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════════

def main():
    print()
    print("█" * 78)
    print("█" + "  Water / proton chain through the framework's bottom-up lens  ".center(76) + "█")
    print("█" * 78)
    print()
    print("Stages 1-6 (algebra → selection) apply directly to the proton qubit:")
    print("  Stage 1  d=2          one proton, two-state |L⟩, |R⟩")
    print("  Stage 2  Pauli σ_a    proton displacement, momentum, parity")
    print("  Stage 3  (σ·σ)² = 3I − 2(σ·σ)    holds on adjacent proton pairs")
    print("  Stage 4  +1 / −3      triplet / singlet of two coupled protons")
    print("  Stage 5  bit_a × bit_b   I↔X, Z↔Y per proton")
    print("  Stage 6  {II, XX, YY, ZZ}    Heisenberg form is admissible;")
    print("           the transverse field −J·X is bit_a-odd one-body —")
    print("           outside the Stage 6 selection.")
    print()
    print("Now Stages 7-9 (Lindbladian → palindrome residual) on three Hamiltonians.")
    print()

    gamma = 0.1
    print(f"  {'N':>3s}  {'Hamiltonian':>34s}  {'‖M‖_op':>12s}  {'spec err':>12s}  {'verdict':>16s}")
    print(f"  {'-' * 3}  {'-' * 34}  {'-' * 12}  {'-' * 12}  {'-' * 16}")

    for N in [2, 3, 4]:
        Sigma_gamma = N * gamma
        # (1) Pure Heisenberg — Stage 6 admissible
        H = H_heisenberg_chain(N, J=1.0)
        L = fw.lindbladian_z_dephasing(H, [gamma] * N)
        v, op_n, sp_e = classify(L, Sigma_gamma, N)
        print(f"  {N:>3d}  {'J(XX+YY+ZZ)':>34s}  {op_n:>12.4e}  {sp_e:>12.4e}  {v:>16s}")

        # (2) Transverse-field Ising at K = 0 (pure transverse field)
        H = H_transverse_field_ising(N, J_tunnel=1.0, K_coupling=0.0)
        L = fw.lindbladian_z_dephasing(H, [gamma] * N)
        v, op_n, sp_e = classify(L, Sigma_gamma, N)
        print(f"  {N:>3d}  {'-J·X (pure transverse)':>34s}  {op_n:>12.4e}  {sp_e:>12.4e}  {v:>16s}")

        # (3) TFI at K=0.2 (the water-script default)
        H = H_transverse_field_ising(N, J_tunnel=1.0, K_coupling=0.2)
        L = fw.lindbladian_z_dephasing(H, [gamma] * N)
        v, op_n, sp_e = classify(L, Sigma_gamma, N)
        print(f"  {N:>3d}  {'-J·X + 0.2·Z·Z (water TFI)':>34s}  {op_n:>12.4e}  {sp_e:>12.4e}  {v:>16s}")

        # (4) Pure ZZ coupling (Stage 6 admissible bilinear, no transverse)
        bonds = [(i, i + 1) for i in range(N - 1)]
        H = fw._build_bilinear(N, bonds, [('Z', 'Z', 0.2)])
        L = fw.lindbladian_z_dephasing(H, [gamma] * N)
        v, op_n, sp_e = classify(L, Sigma_gamma, N)
        print(f"  {N:>3d}  {'0.2·Z·Z (pure ZZ Ising)':>34s}  {op_n:>12.4e}  {sp_e:>12.4e}  {v:>16s}")
        print()

    print("Reading the table:")
    print()
    print("  truly_unbroken  ‖M‖ ≈ 0  AND  spec err ≈ 0       fully palindromic")
    print("  soft_broken     ‖M‖ ≠ 0  BUT  spec err ≈ 0       spectrum lies, vectors scrambled")
    print("  hard_broken     ‖M‖ ≠ 0  AND  spec err ≠ 0       both pairing tests fail")
    print()
    print("  Note: spec err values ~10⁻¹⁵ above are floating-point noise from")
    print("  np.linalg.eigvals; the verdict threshold is 10⁻⁶, so they classify")
    print("  as 'palindromic spectrum'. Pure ZZ gives spec err = 0 exactly because")
    print("  its eigenvalues are integer multiples of γ, no diagonalization needed.")
    print()
    print("The proton-qubit substrate maps directly: framework.py's")
    print("`palindrome_residual` runs on the water Liouvillian without translation.")
    print("The 1925/1976 spectral test (V-Effect's eigenvalue pairing) and the")
    print("framework's operator test give DIFFERENT verdicts on the same H,")
    print("just like on the soft-break trichotomy verified today on Marrakesh.")


if __name__ == "__main__":
    main()
