#!/usr/bin/env python3
"""eq021_obc_sine_basis.py

TASK_OBC_SINE_BASIS_STRUCTURE: investigate the OBC sine-basis origin of the
c_1 bond profile. Four phases:

Phase 1: Verify spectral constants for XY vs Heisenberg, reconcile with F2.
Phase 2: Closed-form sine-basis matrix elements <psi_k | T_b | psi_m> and
         attempt analytical reconstruction of c_1(b) at first order.
Phase 3: Prove mirror symmetry c_1(b) = c_1(N-2-b) from spatial reflection R.
Phase 4: Explain sign pattern; bond profile for psi_2+vac.

No changes to existing repo files. Results go to
simulations/results/eq021_obc_sine_basis/.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig, eigh

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS,
    X, Y, Z, I2, site_op,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode,
    density_matrix, per_site_purity, fit_alpha,
)

RESULTS_DIR = Path(__file__).parent / "results" / "eq021_obc_sine_basis"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------
# Phase 1: spectral constants
# ----------------------------------------------------------------------

def build_H_XY_SE(N, J=1.0):
    """Single-excitation restriction of XY = (J/2) Sum (XX+YY).

    In SE sector, (X_b X_{b+1} + Y_b Y_{b+1})/2 = sigma^+_b sigma^-_{b+1} + h.c.,
    which hops excitation with amplitude 1. Multiplied by J gives hopping J.
    """
    H = np.zeros((N, N), dtype=float)
    for b in range(N - 1):
        H[b, b + 1] = J
        H[b + 1, b] = J
    return H


def build_H_Heis_SE(N, J=1.0):
    """Single-excitation restriction of Heisenberg = J Sum (XX+YY+ZZ).

    (XX+YY) contributes hopping amplitude 2 per bond, times J gives 2J.
    (ZZ) in SE sector: at site i the state |i> has Z_i = -1, others +1.
    Diagonal energy = J * sum_{b=0}^{N-2} <Z_b Z_{b+1}>.
    For state |i>, bonds touching i: b = i-1 (if exists) and b = i (if exists)
    give ZZ = -1. All other N-3 or N-2 bonds give ZZ = +1.
    """
    H = np.zeros((N, N), dtype=float)
    # Hopping
    for b in range(N - 1):
        H[b, b + 1] = 2.0 * J
        H[b + 1, b] = 2.0 * J
    # Diagonal ZZ
    for i in range(N):
        n_bonds_touching = (1 if i > 0 else 0) + (1 if i < N - 1 else 0)
        n_bonds_other = (N - 1) - n_bonds_touching
        H[i, i] = J * (n_bonds_other - n_bonds_touching)
    return H


def build_H_Heis_full(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        Xb = site_op(X, b, N)
        Xb1 = site_op(X, b + 1, N)
        Yb = site_op(Y, b, N)
        Yb1 = site_op(Y, b + 1, N)
        Zb = site_op(Z, b, N)
        Zb1 = site_op(Z, b + 1, N)
        H += J * (Xb @ Xb1 + Yb @ Yb1 + Zb @ Zb1)
    return H


def se_eigvals_closed_form_XY(N, J=1.0):
    """OBC tight-binding with hopping J: E_k = 2J cos(pi k / (N+1)), k=1..N."""
    return np.array([2 * J * np.cos(np.pi * k / (N + 1)) for k in range(1, N + 1)])


def se_eigvals_closed_form_Heis(N, J=1.0):
    """OBC tight-binding with hopping 2J plus diagonal ZZ shift.

    Diagonal is H[i,i] = J * (n_other - n_touch). Not translation invariant,
    so the standard tight-binding formula does NOT apply directly. We just
    diagonalize numerically and report.
    """
    return np.linalg.eigvalsh(build_H_Heis_SE(N, J))


def f2_formula(N, J=1.0):
    """F2 claim: omega_k = 4J (1 - cos(pi k / N)), k = 1, ..., N-1."""
    return np.array([4 * J * (1 - np.cos(np.pi * k / N)) for k in range(1, N)])


def phase1(N_list=(3, 4, 5, 6)):
    print(f"\n{'=' * 70}\nPhase 1: spectral constants\n{'=' * 70}")
    results = {}
    for N in N_list:
        print(f"\nN={N}, J=1:")

        # XY SE
        H_xy_se = build_H_XY_SE(N, J=1.0)
        ev_xy_se = np.linalg.eigvalsh(H_xy_se)
        ev_xy_cf = se_eigvals_closed_form_XY(N, J=1.0)
        print(f"  XY SE eigenvalues (numeric):    {np.sort(ev_xy_se).round(4)}")
        print(f"  XY SE E_k = 2J cos(pi k/(N+1)): {np.sort(ev_xy_cf).round(4)}")
        diff_xy = np.max(np.abs(np.sort(ev_xy_se) - np.sort(ev_xy_cf)))
        print(f"  XY SE |numeric - closed form|:  max = {diff_xy:.2e}")

        # Heisenberg SE
        H_heis_se = build_H_Heis_SE(N, J=1.0)
        ev_heis_se = np.linalg.eigvalsh(H_heis_se)
        print(f"  Heis SE eigenvalues (numeric):  {np.sort(ev_heis_se).round(4)}")

        # F2 formula
        f2 = f2_formula(N, J=1.0)
        print(f"  F2: omega_k = 4J(1-cos(pi k/N)) k=1..N-1: {np.sort(f2).round(4)}")

        # w=1 Liouvillian oscillatory frequencies (XY, take the ones with
        # smallest Re(lambda) > ... to avoid stationary modes)
        H_full = build_H_XY([1.0] * (N - 1), N)
        L = build_liouvillian_matrix(H_full, GAMMA_0, N)
        eigvals = np.linalg.eigvals(L)
        # Sort by Im then by |Re|
        # The w=1 oscillatory modes have Re(lambda) = -2 gamma and Im(lambda) in the SE spectrum
        mask = np.abs(np.real(eigvals) + 2 * GAMMA_0) < 1e-8
        im_vals = np.sort(np.imag(eigvals[mask]))
        im_positive = im_vals[im_vals > 1e-8]
        print(f"  L w=1 osc modes (Re(lam)=-2gamma_0): Im(lam) = {im_positive.round(4)}")

        results[N] = {
            "xy_se_numeric": ev_xy_se.tolist(),
            "xy_se_closed_form": ev_xy_cf.tolist(),
            "xy_se_match_max_diff": float(diff_xy),
            "heis_se_numeric": ev_heis_se.tolist(),
            "f2_formula": f2.tolist(),
            "L_w1_osc_freqs_positive": im_positive.tolist(),
        }

    (RESULTS_DIR / "phase1_spectral_comparison.json").write_text(
        json.dumps(results, indent=2, default=str)
    )
    print(f"\nSaved: {RESULTS_DIR / 'phase1_spectral_comparison.json'}")
    return results


# ----------------------------------------------------------------------
# Phase 2: closed-form matrix elements
# ----------------------------------------------------------------------

def psi_k_amp(k, i, N):
    """Amplitude of |psi_k> at site i (0-indexed)."""
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (i + 1) / (N + 1))


def T_b_matrix_element(k, m, b, N):
    """<psi_k | T_b | psi_m> where T_b = (X_b X_{b+1} + Y_b Y_{b+1}) / 2.

    In SE basis, T_b |i> = delta_{i,b} |b+1> + delta_{i,b+1} |b>.
    So <psi_k | T_b | psi_m> = psi_k(b) psi_m(b+1) + psi_k(b+1) psi_m(b).
    """
    return (psi_k_amp(k, b, N) * psi_k_amp(m, b + 1, N)
            + psi_k_amp(k, b + 1, N) * psi_k_amp(m, b, N))


def phase2_matrix_elements(N_list=(4, 5, 6)):
    print(f"\n{'=' * 70}\nPhase 2: sine-basis matrix elements\n{'=' * 70}")
    results = {}
    for N in N_list:
        print(f"\nN={N}:")
        # Full matrix <psi_k | T_b | psi_m> for k, m = 1..N, b = 0..N-2
        mat_els = np.zeros((N, N, N - 1))
        for k in range(1, N + 1):
            for m in range(1, N + 1):
                for b in range(N - 1):
                    mat_els[k - 1, m - 1, b] = T_b_matrix_element(k, m, b, N)

        # Verify against direct Hamiltonian SE matrix elements
        max_diff = 0.0
        for b in range(N - 1):
            # Build the T_b operator in SE basis directly
            T_b_SE = np.zeros((N, N), dtype=float)
            T_b_SE[b, b + 1] = 1.0
            T_b_SE[b + 1, b] = 1.0
            # Transform: <psi_k | T_b | psi_m> = U^T @ T_b_SE @ U where U[i, k-1] = psi_k(i)
            U = np.zeros((N, N), dtype=float)
            for k in range(1, N + 1):
                for i in range(N):
                    U[i, k - 1] = psi_k_amp(k, i, N)
            mat_els_direct = U.T @ T_b_SE @ U
            diff = np.max(np.abs(mat_els_direct - mat_els[:, :, b]))
            max_diff = max(max_diff, diff)
        print(f"  closed-form vs direct: max |diff| = {max_diff:.2e}")

        # Report a few key matrix elements: <psi_1 | T_b | psi_1> (bond-density)
        diag_11 = mat_els[0, 0, :]
        print(f"  <psi_1 | T_b | psi_1> over bonds: {diag_11.round(4)}")

        # <psi_k | T_0 | psi_1> for k=1..N (endpoint couplings)
        endpoint_col = mat_els[:, 0, 0]
        print(f"  <psi_k | T_0 | psi_1> for k=1..N: {endpoint_col.round(4)}")

        results[N] = {
            "bond_density_psi1": diag_11.tolist(),
            "endpoint_k_couplings_to_psi1": endpoint_col.tolist(),
            "max_verification_diff": float(max_diff),
        }

    (RESULTS_DIR / "phase2_matrix_elements.json").write_text(
        json.dumps(results, indent=2, default=str)
    )
    print(f"\nSaved: {RESULTS_DIR / 'phase2_matrix_elements.json'}")
    return results


# ----------------------------------------------------------------------
# Phase 3: mirror symmetry proof ingredients
# ----------------------------------------------------------------------

def phase3_mirror_check(N_list=(4, 5, 6)):
    print(f"\n{'=' * 70}\nPhase 3: mirror symmetry verification\n{'=' * 70}")
    results = {}
    for N in N_list:
        print(f"\nN={N}:")
        # Build reflection R on SE basis: R|i> = |N-1-i>
        R = np.zeros((N, N), dtype=float)
        for i in range(N):
            R[N - 1 - i, i] = 1.0

        # Check R|psi_k> = (-1)^{k+1} |psi_k>
        psi_phase = []
        for k in range(1, N + 1):
            psi = np.array([psi_k_amp(k, i, N) for i in range(N)])
            Rpsi = R @ psi
            # Compare to (-1)^{k+1} * psi
            sign_pred = (-1) ** (k + 1)
            diff = np.max(np.abs(Rpsi - sign_pred * psi))
            psi_phase.append({"k": k, "predicted_sign": sign_pred,
                              "residual": float(diff)})
            print(f"  R|psi_{k}> vs (-1)^{k+1} |psi_{k}>: residual = {diff:.2e}")

        # Check R T_b R = T_{N-2-b}
        bond_reflection_check = []
        for b in range(N - 1):
            T_b_SE = np.zeros((N, N))
            T_b_SE[b, b + 1] = 1.0
            T_b_SE[b + 1, b] = 1.0
            RTbR = R @ T_b_SE @ R
            # This should equal T_{N-2-b}_SE
            b_ref = N - 2 - b
            T_ref_SE = np.zeros((N, N))
            T_ref_SE[b_ref, b_ref + 1] = 1.0
            T_ref_SE[b_ref + 1, b_ref] = 1.0
            diff = np.max(np.abs(RTbR - T_ref_SE))
            bond_reflection_check.append({"b": b, "b_ref": b_ref,
                                          "residual": float(diff)})
        print(f"  R T_b R = T_{{N-2-b}} residuals: "
              f"max = {max(c['residual'] for c in bond_reflection_check):.2e}")

        results[N] = {
            "psi_phase_check": psi_phase,
            "bond_reflection_check": bond_reflection_check,
        }

    (RESULTS_DIR / "phase3_mirror_verification.json").write_text(
        json.dumps(results, indent=2, default=str)
    )
    print(f"\nSaved: {RESULTS_DIR / 'phase3_mirror_verification.json'}")
    return results


# ----------------------------------------------------------------------
# Phase 4: psi_2+vac bond profile
# ----------------------------------------------------------------------

def propagate_via_eig(eigvals, V_R, V_L_inv, rho_0, times):
    d = rho_0.shape[0]
    rho0_vec = rho_0.flatten(order='F')
    c0 = V_L_inv @ rho0_vec
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order='F')
    return out


def compute_c1_for_bond_state(N, bond_idx, state, times, dJ_list):
    J_A = [J_UNIFORM] * (N - 1)
    H_A = build_H_XY(J_A, N)
    L_A = build_liouvillian_matrix(H_A, GAMMA_0, N)
    ev_A, VR_A = eig(L_A)
    VLinv_A = np.linalg.inv(VR_A)

    rho_0 = density_matrix(state)
    rho_A = propagate_via_eig(ev_A, VR_A, VLinv_A, rho_0, times)
    P_A = per_site_purity(rho_A, N)

    closures = []
    for dJ in dJ_list:
        J_B = list(J_A); J_B[bond_idx] += dJ
        H_B = build_H_XY(J_B, N)
        L_B = build_liouvillian_matrix(H_B, GAMMA_0, N)
        ev_B, VR_B = eig(L_B)
        VLinv_B = np.linalg.inv(VR_B)
        rho_B = propagate_via_eig(ev_B, VR_B, VLinv_B, rho_0, times)
        P_B = per_site_purity(rho_B, N)
        alpha = np.zeros(N)
        for i in range(N):
            a, _ = fit_alpha(times, P_A[:, i], P_B[:, i])
            alpha[i] = a
        closures.append((dJ, float(np.sum(np.log(alpha)))))

    plus = next(c for dJ, c in closures if dJ > 0)
    minus = next(c for dJ, c in closures if dJ < 0)
    dJ_pos = next(dJ for dJ, _ in closures if dJ > 0)
    dJ_neg = next(dJ for dJ, _ in closures if dJ < 0)
    c1 = (plus - minus) / (dJ_pos - dJ_neg)
    return c1


def phase4_psi2_profile(N_list=(4, 5, 6), dJ=0.01):
    print(f"\n{'=' * 70}\nPhase 4: full bond profile for psi_2+vac\n{'=' * 70}")
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    results = {}
    for N in N_list:
        print(f"\nN={N}:")
        state = (vacuum_ket(N) + single_excitation_mode(N, k=2)) / np.sqrt(2.0)
        c1_vec = np.zeros(N - 1)
        for b in range(N - 1):
            t0 = time.time()
            c1 = compute_c1_for_bond_state(N, b, state, times, [-dJ, +dJ])
            c1_vec[b] = c1
            dt = time.time() - t0
            print(f"  bond ({b},{b+1}): c_1 = {c1:+.5f}  ({dt:.1f}s)")
        mirror = float(np.linalg.norm(c1_vec - c1_vec[::-1]))
        print(f"  mirror residual: {mirror:.2e}")
        results[N] = {
            "c_1_vector_psi2": c1_vec.tolist(),
            "mirror_residual": mirror,
        }
        # incremental save
        (RESULTS_DIR / "phase4_psi2_bond_profile.json").write_text(
            json.dumps(results, indent=2, default=str)
        )
    print(f"\nSaved: {RESULTS_DIR / 'phase4_psi2_bond_profile.json'}")
    return results


# ----------------------------------------------------------------------

def main():
    print("=" * 70)
    print("EQ-021 OBC Sine-Basis Structure")
    print("=" * 70)
    p1 = phase1(N_list=(3, 4, 5, 6))
    p2 = phase2_matrix_elements(N_list=(4, 5, 6))
    p3 = phase3_mirror_check(N_list=(4, 5, 6))
    p4 = phase4_psi2_profile(N_list=(4, 5, 6))

    summary = {
        "phase1": p1,
        "phase2": p2,
        "phase3": p3,
        "phase4": p4,
    }
    (RESULTS_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    print(f"\nFinal summary: {RESULTS_DIR / 'summary.json'}")


if __name__ == "__main__":
    main()
