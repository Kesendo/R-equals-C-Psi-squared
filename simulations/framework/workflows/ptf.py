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


def _global_clock_from_evals(evals, tol=1e-9, gap_tol=1e-6):
    """The global clock (Takt + Rotation) from an L_A spectrum, the same reading the
    C# MirrorSystem voices give: gap = slowest nonzero decay rate (the Takt floor,
    = 2γ for a dephasing chain), omega_mem = max |Im λ| among the modes at the gap
    (the Rotation hand of the memory mode), theta_mem = arctan(omega_mem/gap).
    """
    rate = -np.real(evals)
    nonzero = rate[rate > tol]
    if nonzero.size == 0:                      # γ = 0: the clock is stopped
        return {'stopped': True, 'gap': 0.0, 'tau': float('inf'),
                'omega_mem': 0.0, 'theta_mem_deg': 0.0}
    gap = float(nonzero.min())
    at_gap = np.abs(rate - gap) <= gap_tol
    omega = float(np.max(np.abs(np.imag(evals[at_gap]))))
    return {'stopped': False, 'gap': gap, 'tau': 1.0 / gap, 'omega_mem': omega,
            'theta_mem_deg': float(np.degrees(np.arctan2(omega, gap)))}


def perspectives_panel(chain, rho_0, defect_bond, delta_J=0.02, guard_delta_J=0.01,
                       f_max=10.0, consistency_tol=0.5, f_floor=0.5,
                       t_max=20.0, n_t=400):
    """The conductor's stand for one chain and one event: the global clock beside the
    per-perspective rate-of-painting field, with an identifiability guard.

    Two readings on one stand:

    1. The global clock (Takt + Rotation), the mountain's one proper-time hand, read
       from the unperturbed L_A spectrum (gap = 2γ floor, theta_mem = arctan(omega/gap)),
       reusing the eigendecomposition `ptf_alpha_fit` already computed.
    2. The per-perspective field alpha_i, how each site feels the event (a J-defect on
       `defect_bond`), from `ptf_alpha_fit`. A genuine perturbative rescaling has
       (alpha_i - 1) proportional to delta_J, so f_i = (alpha_i - 1)/delta_J should be
       (a) of sane magnitude and (b) stable across two small delta_J. Sites that fail
       either test (a featureless / plateaued trajectory, e.g. one far from the defect,
       fits a confident but meaningless huge alpha at low RMSE; see the (a)-lesson in
       _the_dial_at_many_body and PTF) are flagged unreliable and excluded from the
       closure Sigma ln(alpha).

    HONEST SCOPE: alpha_i is per-(site, observable = purity, state, event); it is NOT an
    intrinsic per-site clock (the "site-local time" reading PERSPECTIVAL_TIME_FIELD
    falsified). The closure Sigma ln(alpha) approx 0 is an empirical regularity (EQ-014),
    reported here over the reliable painters only. Inherits ptf_alpha_fit's dense eig, so
    this is practical for N <= 6; N >= 7 needs the sparse / RK4 path.

    Args:
        chain: ChainSystem with H_type='xy'.
        rho_0: 2^N x 2^N initial density matrix (the state the event is witnessed from).
        defect_bond: bond index of the J-defect (the event).
        delta_J: reported defect strength (relative to chain.J); the alpha field is read here.
        guard_delta_J: the second, smaller delta_J used only for the identifiability guard.
        f_max: a site is "sane" iff |f_i| = |alpha_i - 1|/delta_J <= f_max.
        consistency_tol: a site is "linear" iff |f_i - f_guard_i| <= consistency_tol*(|f_guard_i| + f_floor).
        f_floor: floor in the linearity test so alpha approx 1 (no rescaling) is not mis-flagged.
        t_max, n_t: trajectory horizon and sampling (passed to ptf_alpha_fit).

    Returns dict:
        'clock': {stopped, gap, tau, omega_mem, theta_mem_deg} - the global clock.
        'alphas': per-site alpha_i at delta_J.
        'f', 'f_guard': per-site f_i at delta_J and guard_delta_J.
        'sane', 'linear', 'reliable': per-site booleans (reliable = sane & linear).
        'sigma_log_alpha_reliable': Sigma ln(alpha) over reliable sites (nan if none).
        'sigma_log_alpha_all': Sigma ln(alpha) over all sites.
        'n_unreliable': count of filtered sites.
        'delta_J', 'guard_delta_J', 'defect_bond_index'.
        '_fit': the full ptf_alpha_fit output at delta_J (trajectories, rmses, ...).
    """
    _validate_ptf_inputs(chain, defect_bond)

    out = ptf_alpha_fit(chain, rho_0, defect_bond, J_mod=chain.J + delta_J,
                        t_max=t_max, n_t=n_t)
    out_guard = ptf_alpha_fit(chain, rho_0, defect_bond, J_mod=chain.J + guard_delta_J,
                              t_max=t_max, n_t=n_t)

    clock = _global_clock_from_evals(out['_decomp_A'][0])

    alphas = np.asarray(out['alphas'])
    f = (alphas - 1.0) / delta_J
    f_guard = (np.asarray(out_guard['alphas']) - 1.0) / guard_delta_J
    sane = np.abs(f) <= f_max
    linear = np.abs(f - f_guard) <= consistency_tol * (np.abs(f_guard) + f_floor)
    reliable = sane & linear

    log_a = np.log(np.clip(alphas, 1e-30, None))
    sigma_reliable = float(np.sum(log_a[reliable])) if reliable.any() else float('nan')

    return {
        'clock': clock,
        'alphas': alphas,
        'f': f,
        'f_guard': f_guard,
        'sane': sane,
        'linear': linear,
        'reliable': reliable,
        'sigma_log_alpha_reliable': sigma_reliable,
        'sigma_log_alpha_all': float(np.sum(log_a)),
        'n_unreliable': int((~reliable).sum()),
        'delta_J': float(delta_J),
        'guard_delta_J': float(guard_delta_J),
        'defect_bond_index': defect_bond,
        '_fit': out,
    }
