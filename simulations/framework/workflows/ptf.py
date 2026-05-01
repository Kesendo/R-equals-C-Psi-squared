"""PTF (Perspectival Time Field) workflow: per-site α_i fit + closure.

For an N-site chain under L_A (uniform J, Z-dephasing) and L_B (one bond
perturbed to J_mod), each site's purity trajectory P_B(i, t) is approximately
a one-parameter time rescaling of P_A(i, t):

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

from ..lindblad import bond_perturbation
from ..pauli import site_paulis as _site_paulis
from ._propagation import propagation_setup as _propagation_setup, per_site_bloch_trajectory


_ALPHA_BOUNDS = (0.1, 10.0)
_BOUNDARY_TOL = 1e-3


def _purity_trajectory(evals, R, c0, t_grid, site_paulis):
    """Per-site purity P_i(t) = ½(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²).

    Returns N × len(t_grid) array. Reduces the per-site Bloch trajectory
    via P_i = ½(1 + ‖Bloch_i‖²).
    """
    blochs = per_site_bloch_trajectory(evals, R, c0, t_grid, site_paulis)
    return 0.5 * (1.0 + np.sum(blochs ** 2, axis=2))


def _alpha_fit_one_site(t_grid, P_A, P_B):
    """Fit α minimizing Σ_t (P_A(α·t) − P_B(t))²."""
    interp = interp1d(t_grid, P_A, bounds_error=False,
                      fill_value=(float(P_A[0]), float(P_A[-1])),
                      kind='cubic')

    def mse(alpha):
        d = interp(alpha * t_grid) - P_B
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=_ALPHA_BOUNDS, method='bounded',
                          options={'xatol': 1e-7})
    alpha = float(res.x)
    rmse = float(np.sqrt(res.fun))
    boundary = (abs(alpha - _ALPHA_BOUNDS[0]) < _BOUNDARY_TOL
                or abs(alpha - _ALPHA_BOUNDS[1]) < _BOUNDARY_TOL)
    return alpha, rmse, boundary


def _validate_ptf_inputs(chain, defect_bond):
    if chain.H_type != 'xy':
        raise ValueError(
            f"PTF requires chain.H_type='xy'; got {chain.H_type!r}. "
            "The canonical PTF reference is the XY chain ½(XX+YY) per bond."
        )
    if defect_bond < 0 or defect_bond >= len(chain.bonds):
        raise ValueError(
            f"defect_bond {defect_bond} outside [0, {len(chain.bonds)})"
        )


def ptf_alpha_fit(chain, rho_0, defect_bond, J_mod, t_max=20.0, n_t=400):
    """Per-site α_i fit + closure law for a single-bond J-defect.

    Uses chain.L as L_A (uniform-J XY Heisenberg + Z-dephasing) and
    L_B = L_A + (J_mod − chain.J) · V_L^b where V_L^b is the variation-
    Liouvillian at the defect bond. Propagates ρ_0 under both, samples
    per-site purity at n_t points up to t_max, and fits α_i for each
    site such that P_B(i, t) ≈ P_A(i, α_i·t).

    Args:
        chain: ChainSystem with H_type='xy'.
        rho_0: 2^N × 2^N initial density matrix.
        defect_bond: bond index b ∈ [0, len(chain.bonds)) into chain.bonds.
        J_mod: J value at the defect bond (absolute, not relative).
        t_max: trajectory horizon.
        n_t: number of time samples (uniform grid).

    Returns:
        dict with:
          'alphas': array of per-site α_i (length N).
          'rmses': per-site fit RMSE.
          'on_boundary': per-site boolean — True if α hit the bound.
          'sigma_log_alpha': float, Σ_i ln(α_i) — closure-law residual.
          'defect_bond_index': b.
          'J_mod': J value applied.
          't_grid': time samples used.
          'P_A', 'P_B': N × n_t per-site purity trajectories.
          '_decomp_A': (evals, R, R_inv, c0) of L_A for downstream reuse.
    """
    _validate_ptf_inputs(chain, defect_bond)
    N = chain.N
    bond_pair = chain.bonds[defect_bond]

    L_A = chain.L
    V_L = bond_perturbation(N, bond_pair, kind='XY')
    L_B = L_A + (J_mod - chain.J) * V_L

    t_grid = np.linspace(0, t_max, n_t)
    site_paulis = _site_paulis(N)
    decomp_A = _propagation_setup(L_A, rho_0)
    decomp_B = _propagation_setup(L_B, rho_0)
    P_A = _purity_trajectory(decomp_A[0], decomp_A[1], decomp_A[3], t_grid, site_paulis)
    P_B = _purity_trajectory(decomp_B[0], decomp_B[1], decomp_B[3], t_grid, site_paulis)

    alphas = np.zeros(N)
    rmses = np.zeros(N)
    on_boundary = np.zeros(N, dtype=bool)
    for i in range(N):
        a, r, bd = _alpha_fit_one_site(t_grid, P_A[i], P_B[i])
        alphas[i] = a
        rmses[i] = r
        on_boundary[i] = bd

    sigma_log_alpha = float(np.sum(np.log(np.clip(alphas, 1e-30, None))))

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
        '_decomp_A': decomp_A,
    }


def _slow_modes_from_decomp(decomp_A, n_slow, exclude_stationary=True,
                             stationary_tol=1e-10):
    """Build a slow_modes-shaped dict from a precomputed L_A eigendecomposition,
    avoiding a second eig() call when ptf_alpha_fit already computed it.
    """
    evals, R, R_inv, _c0 = decomp_A
    re_parts = -np.real(evals)
    candidates = np.arange(len(evals))
    if exclude_stationary:
        candidates = candidates[np.abs(re_parts[candidates]) > stationary_tol]
    order = candidates[np.argsort(re_parts[candidates])]
    sel = order[:n_slow]
    return {
        'eigenvalues': evals[sel],
        'right_eigvecs': R[:, sel],
        'left_covecs': R_inv[sel, :],
        'indices': sel,
        'rates': re_parts[sel],
    }


def ptf_painter_panel(chain, rho_0, defect_bond, J_mod, t_max=20.0, n_t=400,
                     include_matrix_elements=False, n_slow=None):
    """Full PTF reading: α_i fit + closure + (optional) V_L matrix elements.

    Combines `ptf_alpha_fit` with the diagnostics-layer perturbation matrix
    elements ⟨W_s | V_L | M_{s'}⟩ on the slow modes, reusing the L_A
    eigendecomposition from `ptf_alpha_fit` to avoid a second eig() call.

    Args:
        chain, rho_0, defect_bond, J_mod, t_max, n_t: see `ptf_alpha_fit`.
        include_matrix_elements: if True, also compute slow-mode matrix
            elements of the bond-perturbation V_L (the "dynamics-of-dynamics").
        n_slow: how many slow modes to include. Default min(d²/4, 50).

    Returns:
        dict combining ptf_alpha_fit output with (if requested):
          'V_L': 4^N × 4^N variation-Liouvillian for the defect bond.
          'slow_modes': slow_modes-shaped dict from L_A.
          'matrix_elements': n_slow × n_slow matrix ⟨W_s | V_L | M_{s'}⟩.
          'eigenvalue_shifts': diag of matrix_elements — first-order δλ_s.
    """
    out = ptf_alpha_fit(chain, rho_0, defect_bond, J_mod, t_max=t_max, n_t=n_t)

    if include_matrix_elements:
        from .lens import _default_n_slow
        from ..diagnostics.ptf import pt_matrix_elements

        bond_pair = chain.bonds[defect_bond]
        V_L = bond_perturbation(chain.N, bond_pair, kind='XY')
        if n_slow is None:
            n_slow = _default_n_slow(chain.N)
        sm = _slow_modes_from_decomp(out['_decomp_A'], n_slow)
        ME = pt_matrix_elements(sm, V_L)

        out['V_L'] = V_L
        out['slow_modes'] = sm
        out['matrix_elements'] = ME
        out['eigenvalue_shifts'] = np.diag(ME)

    return out
