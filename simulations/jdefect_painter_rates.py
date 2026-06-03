"""The painter-rates of the J-defect: the phenomenological twin of the mixing matrix.

The C# telescope (JDefectField, `inspect --root between --axis jdefect`) reads the J-defect's
in-between as the first-order eigenvector mixing <W_s|V_L|M_s'>: a kernel protected by U(1) and a
live off-diagonal. This script reads the SAME defect from the other side, the phenomenological side
PTF named: the per-site painter-rates alpha_i (how far the defect rescales each site's purity
trajectory, P_B(i, t) ~ P_A(i, alpha_i * t)) and their closure Sum ln(alpha_i) ~ 0, the contract
among the painters.

The two are one structure. The matrix elements <W_s|V_L|M_s'> printed here (kernel shifts ~0,
off-diagonal alive) are the same object the C# SlowModeMixing computes; the painter-rates are what
that mixing does to the site-resolved purity over finite time. "Between N painters around a mountain,
the mountain happens" (hypotheses/PERSPECTIVAL_TIME_FIELD.md).

Run: python simulations/jdefect_painter_rates.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw


def bonding_mode_state(N):
    """psi = (|vac> + |psi_1>) / sqrt(2), the canonical PTF initial state (the bonding mode)."""
    psi = np.zeros(2 ** N, dtype=complex)
    psi[0] = 1.0
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * (i + 1) / (N + 1))
        psi[2 ** (N - 1 - i)] += amp
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def painter_rates(N, gamma, defect_bond, delta_J, t_max=20.0, n_t=400):
    """The painter panel for one J-defect, plus the mixing-matrix connection."""
    chain = fw.ChainSystem(N=N, gamma_0=gamma, J=1.0, H_type='xy')
    rho_0 = bonding_mode_state(N)

    panel = fw.perspectives_panel(chain, rho_0, defect_bond=defect_bond,
                                  delta_J=delta_J, t_max=t_max, n_t=n_t)

    clock = panel['clock']
    alphas = panel['alphas']
    f = panel['f']
    reliable = panel['reliable']

    print(f"=== J-defect painter-rates (N={N}, gamma={gamma}, bond {defect_bond}=(0,1), dJ={delta_J}) ===")
    print(f"clock: Takt gap={clock['gap']:.4f} (=2*gamma), tau={clock['tau']:.3f}; "
          f"Rotation omega_mem={clock['omega_mem']:.4f}, theta_mem={clock['theta_mem_deg']:.1f} deg")
    print()
    print("painter   alpha_i    f_i=(a-1)/dJ   reliable")
    for i in range(N):
        tag = "faster" if alphas[i] > 1.0 else ("slower" if alphas[i] < 1.0 else "same")
        print(f"  site {i}   {alphas[i]:.4f}    {f[i]:+7.3f}      {str(bool(reliable[i])):5s}  ({tag})")
    print()
    print(f"closure  Sum ln(alpha) = {panel['sigma_log_alpha_all']:+.4f} (all sites); "
          f"{panel['sigma_log_alpha_reliable']:+.4f} (reliable, {int(reliable.sum())}/{N})")
    print("         the contract among the painters: ~0 means no painter invents or omits")
    print()

    # The connection: the matrix elements <W_s|V_L|M_s'> are the same object the C# SlowModeMixing
    # computes (same V_L = dL/dJ on the defect bond, same slow modes WITH the kernel). The kernel is
    # protected; the off-diagonal is the mixing that drives the painter-rates above.
    V_L = fw.bond_perturbation(N, chain.bonds[defect_bond], kind='XY')
    sm = fw.slow_modes(chain, n_slow=16, exclude_stationary=False)
    ME = fw.pt_matrix_elements(sm, V_L)
    evals = sm['eigenvalues']
    shifts = np.diag(ME)
    kernel = np.abs(evals) < 1e-6
    off = ME - np.diag(np.diag(ME))
    off_mass = float(np.linalg.norm(off))
    max_kernel_shift = float(np.max(np.abs(shifts[kernel]))) if kernel.any() else float('nan')
    print("connection to the C# mixing matrix (the same <W_s|V_L|M_s'>):")
    print(f"  kernel modes {int(kernel.sum())}: max |first-order shift| = {max_kernel_shift:.2e} (protected, U(1))")
    print(f"  off-diagonal mixing mass = {off_mass:.3f} (alive: the eigenvectors mix)")
    print("  -> the painter-rates ARE this mixing, read over finite time on the site purities")
    return panel


if __name__ == '__main__':
    # gamma = 0.5 matches the C# telescope's J-defect axis.
    painter_rates(N=5, gamma=0.5, defect_bond=0, delta_J=0.05)
    print()
    # gamma = 0.05 is PTF's canonical reference, where the closure law was validated.
    painter_rates(N=5, gamma=0.05, defect_bond=0, delta_J=0.05)
