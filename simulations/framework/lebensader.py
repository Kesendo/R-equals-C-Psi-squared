"""Lebensader cockpit_panel — composes skeleton + trace + cusp + chiral + Y-parity.

The Lebensader is the bridge concept that holds operator-level skeleton
(Π-protected observable count) and state-level trace (θ-trajectory geometry)
together via Π·L·Π⁻¹ + L + 2Σγ·I = 0. cockpit_panel computes both views
in one pass plus the cusp-pattern classifier, chiral panel, and Y-parity panel.

Hardware-confirmed: EQ-030 closure on Marrakesh April 26 (drop=28 for YZ+ZY,
Pearson(drop, Δ∫θ) = +0.85, Bures velocity null as third axis).

Public API:
  cockpit_panel(H, gamma_l, rho_0, N, gamma_t1_l=None, t_max, dt, ...)
"""
from __future__ import annotations

import math

import numpy as np

from .lindblad import lindbladian_z_dephasing, lindbladian_z_plus_t1
from .observables import pi_protected_observables
from .pauli import _vec_to_pauli_basis_transform, pauli_basis_vector
from .symmetry import chiral_panel, y_parity_panel


def _purity_psi_norm_cpsi(rho):
    """CΨ glossary: Purity·Ψ-norm. Ψ-norm = L1(ρ)/(d-1)."""
    p = float(np.real(np.trace(rho @ rho)))
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi_n = l1 / (d - 1)
    return p, psi_n, p * psi_n


def _theta_from_cpsi(c):
    if c <= 0.25:
        return 0.0
    return math.degrees(math.atan(math.sqrt(4 * c - 1)))


def _cluster_eigenvalues(evals, tol=1e-8):
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    clusters = []
    for i in range(n):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n):
            if not used[j] and abs(evals[j] - evals[i]) < tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)
    return clusters


def _trajectory_via_eigendecomp(L_pauli, V, Vinv, rho_0, N, times):
    """ρ(t) at sample times, via precomputed eigendecomposition."""
    evals = np.diag(Vinv @ L_pauli @ V).copy()
    M_basis = _vec_to_pauli_basis_transform(N)
    rho_pauli_0 = pauli_basis_vector(rho_0, N)
    c0 = Vinv @ rho_pauli_0
    out = []
    d = 2 ** N
    for t in times:
        rho_pauli_t = V @ (c0 * np.exp(evals * t))
        rho_vec = M_basis @ rho_pauli_t
        out.append(rho_vec.reshape(d, d).T)
    return out


def _find_cpsi_crossings(times, cpsi_t, threshold=0.25):
    """Linear-interpolate CΨ - threshold zeros."""
    out = []
    for i in range(len(times) - 1):
        a, b = cpsi_t[i] - threshold, cpsi_t[i + 1] - threshold
        if a * b < 0:
            frac = a / (a - b)
            t_cross = times[i] + frac * (times[i + 1] - times[i])
            out.append((float(t_cross), -1 if a > 0 else +1, i))
    return out


def cockpit_panel(H, gamma_l, rho_0, N,
                  gamma_t1_l=None,
                  t_max=10.0, dt=0.005,
                  threshold=1e-9, cluster_tol=1e-8):
    """Lebensader cockpit panel: skeleton + trace + cusp + chiral + Y-parity.

    Returns dict with:
      'lebensader':  {skeleton, trace, rating, skeleton_status, trace_status}
      'cusp':        {n_crossings, dominant_eigenvalue, pattern, mode_type}
      'chiral':      chiral_panel output
      'y_parity':    y_parity_panel output
      '_trajectory_for_inspection': {times, cpsi, theta}

    Skeleton: Π-protected counts under pure-Z and (optionally) +T1; drop is
    the difference. Trace: θ-trajectory metrics (max, tail duration, α
    descent exponent). Rating combines both with status labels.
    """
    if gamma_t1_l is None:
        gamma_t1_l = [0.0] * N
    times = np.linspace(0.0, t_max, int(t_max / dt) + 1)

    # Skeleton: Π-protected counts pure-Z + T1
    pi_pure = pi_protected_observables(
        H, gamma_l, rho_0, N, threshold=threshold, cluster_tol=cluster_tol,
    )
    n_pure = len(pi_pure['protected'])

    if any(g != 0 for g in gamma_t1_l):
        L_t1 = lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
        M_basis = _vec_to_pauli_basis_transform(N)
        L_t1_pauli = (M_basis.conj().T @ L_t1 @ M_basis) / (2 ** N)
        evals_t1, V_t1 = np.linalg.eig(L_t1_pauli)
        Vinv_t1 = np.linalg.inv(V_t1)
        c_t1 = Vinv_t1 @ pauli_basis_vector(rho_0, N)
        clusters_t1 = _cluster_eigenvalues(evals_t1, tol=cluster_tol)
        n_t1 = 0
        for alpha in range(1, 4 ** N):
            max_S = 0.0
            for cl in clusters_t1:
                S = sum(V_t1[alpha, k] * c_t1[k] for k in cl)
                max_S = max(max_S, abs(S))
            if max_S < threshold:
                n_t1 += 1
    else:
        n_t1 = n_pure

    skeleton = {'n_protected_pure': n_pure, 'n_protected_t1': n_t1,
                'drop': n_pure - n_t1}

    # Trace: θ-trajectory under L_t1 (if T1) else L_pure
    if any(g != 0 for g in gamma_t1_l):
        L_active = L_t1
        M_basis = _vec_to_pauli_basis_transform(N)
        L_active_pauli = (M_basis.conj().T @ L_active @ M_basis) / (2 ** N)
        evals_active, V_active = np.linalg.eig(L_active_pauli)
        Vinv_active = np.linalg.inv(V_active)
    else:
        L_active = lindbladian_z_dephasing(H, gamma_l)
        M_basis = _vec_to_pauli_basis_transform(N)
        L_active_pauli = (M_basis.conj().T @ L_active @ M_basis) / (2 ** N)
        evals_active, V_active = np.linalg.eig(L_active_pauli)
        Vinv_active = np.linalg.inv(V_active)

    traj = _trajectory_via_eigendecomp(L_active_pauli, V_active, Vinv_active,
                                        rho_0, N, times)
    cpsi_t = np.array([_purity_psi_norm_cpsi(r)[2] for r in traj])
    theta_t = np.array([_theta_from_cpsi(c) for c in cpsi_t])

    crossings = _find_cpsi_crossings(times, cpsi_t)
    n_crossings = len(crossings)
    theta_max = float(theta_t.max())

    above = theta_t > 0
    if above.any():
        last_above_idx = len(times) - 1 - int(np.argmax(above[::-1]))
        t_last_above = float(times[last_above_idx])
    else:
        t_last_above = 0.0

    tail_mask = (theta_t > 0) & (theta_t < 5.0)
    if tail_mask.any():
        tail_idx = np.where(tail_mask)[0]
        tail_duration = float(times[tail_idx[-1]] - times[tail_idx[0]])
    else:
        tail_duration = 0.0

    alpha_descent = None
    if crossings:
        t_cross_last = crossings[-1][0]
        win_mask = ((times >= t_cross_last - 0.3) & (times < t_cross_last)
                    & (theta_t > 0.01) & (theta_t < 10.0))
        if win_mask.sum() >= 4:
            x = np.log(np.maximum(t_cross_last - times[win_mask], 1e-9))
            y = np.log(theta_t[win_mask])
            coef = np.polyfit(x, y, 1)
            alpha_descent = float(coef[0])

    # Cusp pattern + dominant mode
    if not crossings:
        pattern = 'never crosses'
        dom_eigval = None
        dom_eigval_nonzero = None
        mode_type = 'n/a'
    else:
        t_first, _, idx_first = crossings[0]
        rho_at_cross = traj[idx_first]
        rho_pauli_at = pauli_basis_vector(rho_at_cross, N)
        c_at = Vinv_active @ rho_pauli_at
        clusters = _cluster_eigenvalues(evals_active, tol=cluster_tol)
        cluster_norms = []
        for cl in clusters:
            norm_cl = float(np.linalg.norm([c_at[k] for k in cl]))
            cluster_norms.append((norm_cl, cl))
        cluster_norms.sort(key=lambda x: -x[0])
        dom_eigval = complex(evals_active[cluster_norms[0][1][0]])
        nonzero_norms = [(n, cl) for n, cl in cluster_norms
                          if abs(complex(evals_active[cl[0]])) > 1e-6]
        if nonzero_norms:
            dom_eigval_nonzero = complex(evals_active[nonzero_norms[0][1][0]])
        else:
            dom_eigval_nonzero = None

        pattern = 'monotonic' if n_crossings == 1 else 'heartbeat'

        ref_eigval = (dom_eigval_nonzero if abs(dom_eigval) < 1e-6
                       else dom_eigval)
        if ref_eigval is None:
            mode_type = 'pure steady state'
        elif abs(dom_eigval) < 1e-6:
            if abs(ref_eigval.imag) > 1e-3:
                mode_type = 'steady-state + oscillatory sub-mode'
            else:
                mode_type = 'steady-state + real-decay sub-mode'
        elif abs(ref_eigval.imag) > 1e-3:
            mode_type = 'oscillatory'
        else:
            mode_type = 'real-decay'

    classification = (f"{pattern} | {mode_type}"
                      if pattern != 'never crosses' else 'never crosses')

    trace = {'theta_max': theta_max,
             'n_crossings': n_crossings,
             't_last_theta_positive': t_last_above,
             'tail_duration_sub5deg': tail_duration,
             'alpha_descent': alpha_descent}

    cusp = {'n_crossings': n_crossings,
            'dominant_eigenvalue_at_first_crossing': dom_eigval,
            'dominant_nonzero_eigenvalue': dom_eigval_nonzero,
            'pattern': pattern,
            'mode_type': mode_type,
            'classification': classification}

    skel_status = 'intact' if skeleton['drop'] <= 1 else (
        'partial' if skeleton['drop'] <= 5 else 'broken'
    )
    tail_threshold = max(0.005, 0.86 / (2 ** max(0, N - 3)))
    trace_status = 'long' if tail_duration > tail_threshold else (
        'short' if tail_duration > tail_threshold / 10 else 'absent'
    )
    if skel_status == 'intact' and trace_status == 'long':
        leb_rating = 'intact (skeleton holds, trace lives)'
    elif skel_status == 'intact' and trace_status != 'long':
        leb_rating = 'partial (skeleton holds, trace short)'
    elif skel_status != 'intact' and trace_status == 'long':
        leb_rating = 'partial (skeleton bleeds, trace persists)'
    else:
        leb_rating = 'collapsed (both skeleton and trace gone)'

    chiral = chiral_panel(H, rho_0, N)
    y_parity = y_parity_panel(H, gamma_l, rho_0, N, gamma_t1_l=gamma_t1_l)

    return {
        'lebensader': {'skeleton': skeleton, 'trace': trace,
                       'rating': leb_rating,
                       'skeleton_status': skel_status,
                       'trace_status': trace_status,
                       'tail_threshold_for_N': tail_threshold},
        'cusp': cusp,
        'chiral': chiral,
        'y_parity': y_parity,
        '_trajectory_for_inspection': {'times': times, 'cpsi': cpsi_t,
                                        'theta': theta_t},
    }
