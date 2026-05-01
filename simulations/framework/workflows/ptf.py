"""PTF (Perspectival Time Field) workflow: per-site α_i fit + closure.

For an N-site chain under L_A (uniform J=1, Z-dephasing) and L_B (J on a single
defect bond perturbed to J_mod, others unchanged), each site's purity trajectory
P_B(i, t) is approximately a one-parameter time rescaling of P_A(i, t):

    P_B(i, t) ≈ P_A(i, α_i · t)

The α_i depend on the initial state. The closure law Σ_i ln(α_i) ≈ 0 holds
within ±0.05 in the perturbative window |δJ| ≤ 0.1 across multiple initial
states tested at N=7 (Tier 2; downgraded by EQ-014 from "first-order theorem"
to "perturbative empirical bound"). See hypotheses/PERSPECTIVAL_TIME_FIELD.md.

Public API:
  ptf_alpha_fit(chain, rho_0, defect_bond, J_mod, t_max=20.0, n_t=400)
  ptf_painter_panel(chain, rho_0, defect_bond, J_mod, t_max=20.0, n_t=400, ...)
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

from ..lindblad import lindbladian_z_dephasing, bond_perturbation
from ..pauli import _build_bilinear, site_op


# Per-site Pauli operators are reused across multiple time samples; cache by N.
_SITE_PAULI_CACHE = {}


def _site_paulis(N):
    """Return list of (X_i, Y_i, Z_i) operators for i = 0..N-1, cached by N."""
    if N not in _SITE_PAULI_CACHE:
        _SITE_PAULI_CACHE[N] = [
            (site_op(N, i, 'X'), site_op(N, i, 'Y'), site_op(N, i, 'Z'))
            for i in range(N)
        ]
    return _SITE_PAULI_CACHE[N]


def _build_xy_chain_H(chain, J_per_bond):
    """Build XY-chain Hamiltonian H = ½ Σ_b J_b (X_b X_{b+1} + Y_b Y_{b+1})."""
    bonds = list(chain.bonds)
    if len(J_per_bond) != len(bonds):
        raise ValueError(
            f"J_per_bond has length {len(J_per_bond)}, expected {len(bonds)}"
        )
    terms = []
    for b, (i, j) in enumerate(bonds):
        terms.append(('X', 'X', 0.5 * J_per_bond[b]))
        terms.append(('Y', 'Y', 0.5 * J_per_bond[b]))
    # _build_bilinear expects same J across all terms per bond, so build per bond
    H = np.zeros((chain.d, chain.d), dtype=complex)
    for b, (i, j) in enumerate(bonds):
        bond_terms = [('X', 'X', 0.5 * J_per_bond[b]),
                      ('Y', 'Y', 0.5 * J_per_bond[b])]
        H = H + _build_bilinear(chain.N, [(i, j)], bond_terms)
    return H


def _purity_trajectory(L, rho_0, t_grid, site_paulis):
    """Per-site purity trajectory P_i(t) = ½(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²).

    Uses spectral decomposition of L for closed-form propagation:
      ρ(t) = R · diag(exp(λ_s · t)) · R⁻¹ · ρ_0

    Args:
        L: 4^N × 4^N Liouvillian.
        rho_0: 2^N × 2^N initial density matrix.
        t_grid: time samples.
        site_paulis: list of (X_i, Y_i, Z_i) per site.

    Returns:
        N × len(t_grid) array of per-site purities.
    """
    d = rho_0.shape[0]
    N = len(site_paulis)
    rho_vec = rho_0.flatten('F').astype(complex)
    evals, R = np.linalg.eig(L)
    c0 = np.linalg.solve(R, rho_vec)
    purities = np.zeros((N, len(t_grid)))
    for ti, t in enumerate(t_grid):
        rho_t = (R @ (np.exp(evals * t) * c0)).reshape(d, d, order='F')
        rho_t = 0.5 * (rho_t + rho_t.conj().T)  # numerical Hermiticity
        for i, (Xi, Yi, Zi) in enumerate(site_paulis):
            x = float(np.real(np.trace(Xi @ rho_t)))
            y = float(np.real(np.trace(Yi @ rho_t)))
            z = float(np.real(np.trace(Zi @ rho_t)))
            purities[i, ti] = 0.5 * (1.0 + x * x + y * y + z * z)
    return purities


def _alpha_fit_one_site(t_grid, P_A, P_B, alpha_bounds=(0.1, 10.0)):
    """Fit α minimizing Σ_t (P_A(α·t) − P_B(t))².

    Returns (alpha, rmse, on_boundary).
    """
    interp = interp1d(t_grid, P_A, bounds_error=False,
                      fill_value=(float(P_A[0]), float(P_A[-1])),
                      kind='cubic')

    def mse(alpha):
        d = interp(alpha * t_grid) - P_B
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=alpha_bounds, method='bounded',
                          options={'xatol': 1e-7})
    alpha = float(res.x)
    rmse = float(np.sqrt(res.fun))
    boundary = (abs(alpha - alpha_bounds[0]) < 1e-3
                or abs(alpha - alpha_bounds[1]) < 1e-3)
    return alpha, rmse, boundary


def ptf_alpha_fit(chain, rho_0, defect_bond, J_mod, t_max=20.0, n_t=400):
    """Per-site α_i fit + closure law for a single-bond J-defect.

    Builds L_A (uniform J=1) and L_B (J=J_mod at defect_bond, J=1 elsewhere)
    on the chain's XY Hamiltonian + Z-dephasing. Propagates ρ_0 under both,
    samples per-site purity at n_t points up to t_max, and fits α_i for each
    site such that P_B(i, t) ≈ P_A(i, α_i·t).

    Args:
        chain: ChainSystem (chain.H_type='xy' is assumed for the PTF context;
            other H_types still work but the canonical PTF reference is XY).
        rho_0: 2^N × 2^N initial density matrix.
        defect_bond: bond index b ∈ [0, N-2) into chain.bonds (NOT a (i,j) tuple).
        J_mod: J value at the defect bond. Uniform J=1 elsewhere.
        t_max: trajectory horizon.
        n_t: number of time samples (uniform grid).

    Returns:
        dict with:
          'alphas': array of per-site α_i (length N).
          'rmses': per-site fit RMSE.
          'on_boundary': per-site boolean — True if α hit the [0.1, 10] bound.
          'sigma_log_alpha': float, Σ_i ln(α_i) — the closure-law residual.
          'defect_bond_index': b.
          'J_mod': J value applied.
          't_grid': time samples used.
          'P_A': N × n_t per-site purity trajectories under L_A.
          'P_B': N × n_t per-site purity trajectories under L_B.
    """
    if defect_bond < 0 or defect_bond >= len(chain.bonds):
        raise ValueError(
            f"defect_bond {defect_bond} outside [0, {len(chain.bonds)})"
        )
    N = chain.N
    n_bonds = len(chain.bonds)

    # Hamiltonians: A uniform, B with one bond modulated
    J_A = [1.0] * n_bonds
    J_B = list(J_A)
    J_B[defect_bond] = J_mod
    H_A = _build_xy_chain_H(chain, J_A)
    H_B = _build_xy_chain_H(chain, J_B)

    gamma_l = [chain.gamma_0] * N
    L_A = lindbladian_z_dephasing(H_A, gamma_l)
    L_B = lindbladian_z_dephasing(H_B, gamma_l)

    t_grid = np.linspace(0, t_max, n_t)
    site_paulis = _site_paulis(N)
    P_A = _purity_trajectory(L_A, rho_0, t_grid, site_paulis)
    P_B = _purity_trajectory(L_B, rho_0, t_grid, site_paulis)

    alphas = np.zeros(N)
    rmses = np.zeros(N)
    on_boundary = np.zeros(N, dtype=bool)
    for i in range(N):
        a, r, bd = _alpha_fit_one_site(t_grid, P_A[i], P_B[i])
        alphas[i] = a
        rmses[i] = r
        on_boundary[i] = bd

    safe_alphas = np.clip(alphas, 1e-30, None)
    sigma_log_alpha = float(np.sum(np.log(safe_alphas)))

    return {
        'alphas': alphas,
        'rmses': rmses,
        'on_boundary': on_boundary,
        'sigma_log_alpha': sigma_log_alpha,
        'defect_bond_index': defect_bond,
        'J_mod': float(J_mod),
        't_grid': t_grid,
        'P_A': P_A,
        'P_B': P_B,
    }


def ptf_painter_panel(chain, rho_0, defect_bond, J_mod, t_max=20.0, n_t=400,
                     include_matrix_elements=False, n_slow=None):
    """Full PTF reading: α_i fit + closure + (optional) V_L matrix elements.

    Combines `ptf_alpha_fit` with the diagnostics-layer perturbation matrix
    elements ⟨W_s | V_L | M_{s'}⟩ on the slow modes, giving access to both
    the empirical α_i pattern and the structural drivers (V_L, slow modes).

    Args:
        chain, rho_0, defect_bond, J_mod, t_max, n_t: see `ptf_alpha_fit`.
        include_matrix_elements: if True, also compute the slow-mode matrix
            elements of the bond-perturbation V_L (the "dynamics-of-dynamics"
            superoperator). Adds an eigendecomposition of L_A.
        n_slow: how many slow modes to include in the matrix-element panel.
            Default min(d²/4, 50).

    Returns:
        dict combining ptf_alpha_fit output with (if requested):
          'V_L': 4^N × 4^N variation-Liouvillian for the defect bond.
          'slow_modes': output of `framework.slow_modes(chain, ...)` on L_A.
          'matrix_elements': n_slow × n_slow matrix ⟨W_s | V_L | M_{s'}⟩.
          'eigenvalue_shifts': diag of matrix_elements — first-order δλ_s.
    """
    out = ptf_alpha_fit(chain, rho_0, defect_bond, J_mod, t_max=t_max, n_t=n_t)

    if include_matrix_elements:
        from .lens import slow_modes
        from ..diagnostics.ptf import pt_matrix_elements

        i, j = chain.bonds[defect_bond]
        V_L = bond_perturbation(chain.N, (i, j), kind='XY')

        # Slow modes of L_A (uniform-J Lindbladian)
        gamma_l = [chain.gamma_0] * chain.N
        H_A = _build_xy_chain_H(chain, [1.0] * len(chain.bonds))
        L_A = lindbladian_z_dephasing(H_A, gamma_l)
        sm = slow_modes(chain, n_slow=n_slow, L=L_A)
        ME = pt_matrix_elements(sm, V_L)

        out['V_L'] = V_L
        out['slow_modes'] = sm
        out['matrix_elements'] = ME
        out['eigenvalue_shifts'] = np.diag(ME)

    return out
