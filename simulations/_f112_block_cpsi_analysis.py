"""F112 polarity-asymmetry lens on ibm_kingston block-CPsi-saturation trajectory.

Dataset: `data/ibm_block_cpsi_saturation_may2026/block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json`
- Initial state: (|D_0⟩ + |D_1⟩)/√2 = (|00⟩ + (|01⟩+|10⟩)/√2)/√2 on Kingston qubits 13, 14
- 5 t-points: 0, 120, 240, 360, 480 μs
- 16 Pauli expectations per t-point (full 2-qubit tomography)
- T2_min calibration: 480 μs; γ_eff = 1/T2 ≈ 0.00208 per μs
- Documented anomaly: hardware C_block decays ~1.72× faster than pure-T2 predicts

F112 hypothesis to test:
  Does the hardware-effective noise model that explains the trajectory sit
  inside F112's typed Tier1Derived scope (Hermitian H + bit_b-homogeneous c)?

Method:
  For each of several candidate L models (pure Z-deph, +ZZ crosstalk, +T1 σ⁻
  amplitude damping, +transverse h_y field), fit parameters to the 5-point
  ρ(t) trajectory via least-squares, compute:
    1. Trajectory fit residual: ‖ρ_predicted(t) − ρ_observed(t)‖
    2. F112 polarity asymmetry on the fitted L
    3. bit_b-homogeneity classification of each c_k in the model

  Models with all-bit_b-homogeneous c (pure-Z + ZZ crosstalk) are in F112's
  Tier1Derived scope; F112 says asymmetry = 0 bit-exact. Models with σ⁻ T1
  fall outside the typed scope but observed empirically balanced (probe 5).
  Transverse h_y is Hermitian H + bit_b-mixed H → also outside Hermitian
  bit_b-homogeneous H requirement of F112's exact statement.

  The diagnostic: which model fits the data, and where does it sit relative
  to F112's typed scope?

Output:
  Table per model: fit_RMS, F112 asymmetry, in_F112_scope flag.
  Interpretation: which model best explains the 1.72× anomaly and what
  F112 says about it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw  # noqa: E402

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 0], [1, 0]], dtype=complex)  # |0><1|
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

DATA_PATH = Path('data/ibm_block_cpsi_saturation_may2026/'
                 'block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json')


def load_trajectory():
    """Load 5-snapshot trajectory; return (t_us_array, [rho_4x4 per t])."""
    with open(DATA_PATH, encoding='utf-8') as f:
        d = json.load(f)
    t_grid = np.array(d['t_grid_us'], dtype=float)
    rhos = []
    for snap in d['t_snapshots']:
        exps = snap['expectations']
        rho = np.zeros((4, 4), dtype=complex)
        for key, val in exps.items():
            a, b = key.split(',')
            P = np.kron(PAULI[a], PAULI[b])
            rho = rho + float(val) * P
        rho = rho / 4.0  # 2 qubits: factor 1/2^N
        rhos.append(rho)
    return t_grid, rhos, d


def site_op(N, site, mat2):
    """2-site Pauli operator placed at `site`, identity elsewhere."""
    mats = [I2] * N
    mats[site] = mat2
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def vec(rho):
    return rho.flatten('F')


def devec(v, d=4):
    return v.reshape((d, d), order='F')


def build_L_model(params, model):
    """Return (L_vec, c_list, gamma_list, H_pauli_terms_list).

    H Hilbert form, c_ops as Hilbert matrices, gammas as scalars.
    For F112-scope check we also return the Pauli-string letter-tuples
    that constitute c (where applicable; some c are not single-Pauli).
    """
    N = 2
    d = 4
    Id = np.eye(d, dtype=complex)
    c_list = []
    g_list = []
    c_kind = []
    H = np.zeros((d, d), dtype=complex)

    if model == 'pure_Z':
        gz0, gz1 = params  # per-qubit Z-deph rates
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z)]
        g_list = [gz0, gz1]
        c_kind = ['Z@0', 'Z@1']

    elif model == 'Z_plus_T1':
        gz0, gz1, gt0, gt1 = params  # +per-qubit σ⁻ T1 rates
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z),
                  site_op(N, 0, SIGMA_MINUS), site_op(N, 1, SIGMA_MINUS)]
        g_list = [gz0, gz1, gt0, gt1]
        c_kind = ['Z@0', 'Z@1', 'sigma_minus@0', 'sigma_minus@1']

    elif model == 'Z_plus_ZZ':
        gz0, gz1, j_zz = params  # + ZZ-crosstalk in H
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z)]
        g_list = [gz0, gz1]
        c_kind = ['Z@0', 'Z@1']
        H = j_zz * np.kron(Z, Z)

    elif model == 'Z_plus_T1_plus_ZZ':
        gz0, gz1, gt0, gt1, j_zz = params
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z),
                  site_op(N, 0, SIGMA_MINUS), site_op(N, 1, SIGMA_MINUS)]
        g_list = [gz0, gz1, gt0, gt1]
        c_kind = ['Z@0', 'Z@1', 'sigma_minus@0', 'sigma_minus@1']
        H = j_zz * np.kron(Z, Z)

    elif model == 'Z_plus_hy':
        gz0, gz1, hy = params  # + single-site Y field
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z)]
        g_list = [gz0, gz1]
        c_kind = ['Z@0', 'Z@1']
        H = hy * (site_op(N, 0, Y) + site_op(N, 1, Y))

    else:
        raise ValueError(f"Unknown model: {model}")

    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c, g in zip(c_list, g_list):
        c_dag_c = c.conj().T @ c
        anti = 0.5 * (np.kron(c_dag_c, Id) + np.kron(Id, c_dag_c.T))
        L_vec = L_vec + g * (np.kron(c, c.conj()) - anti)
    return L_vec, c_list, g_list, c_kind, H


def evolve(L_vec, rho0, t):
    """exp(L*t) acting on vec(rho0); returns ρ(t)."""
    v0 = vec(rho0)
    vt = expm(L_vec * t) @ v0
    return devec(vt)


def fit_residual(params, model, t_us, rhos):
    """Sum of Frobenius² distances between predicted and observed ρ(t)."""
    # Clip params to non-negative for rates; allow free sign for H couplings
    if model in ('pure_Z', 'Z_plus_T1'):
        if np.any(np.array(params) < 0):
            return 1e6
    elif model == 'Z_plus_ZZ':
        if np.any(np.array(params[:2]) < 0):
            return 1e6
    elif model == 'Z_plus_T1_plus_ZZ':
        if np.any(np.array(params[:4]) < 0):
            return 1e6
    elif model == 'Z_plus_hy':
        if np.any(np.array(params[:2]) < 0):
            return 1e6
    try:
        L_vec, _, _, _, _ = build_L_model(params, model)
    except Exception:
        return 1e6
    total = 0.0
    rho0 = rhos[0]
    for i, t in enumerate(t_us):
        if i == 0:
            continue
        rho_pred = evolve(L_vec, rho0, t)
        diff = rho_pred - rhos[i]
        total += float(np.sum(np.abs(diff) ** 2))
    return total


def fit_model(model, t_us, rhos, x0):
    """Local minimization from initial guess x0."""
    result = minimize(
        fit_residual, x0, args=(model, t_us, rhos),
        method='Nelder-Mead',
        options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 5000},
    )
    return result.x, float(result.fun)


def is_bit_b_homogeneous_pauli_label(label):
    """For c labelled as 'X@k', 'Y@k', 'Z@k', 'sigma_minus@k', etc., return
    True iff the operator is bit_b-homogeneous as a Pauli sum.

    bit_b = (#Y + #Z) mod 2 per Pauli string. Single-Pauli letters:
    X bit_b = 0; Y, Z bit_b = 1; I bit_b = 0.
    sigma_minus = (X − iY)/2: mixed bit_b (X is 0, Y is 1) → False.
    sigma_plus = (X + iY)/2: same.
    """
    if label.startswith('X@') or label.startswith('Y@') or label.startswith('Z@'):
        return True
    if label.startswith('sigma_'):
        return False
    return False


def run_polarity_on_L(L_vec, N=2, sigma=None):
    """Transform L_vec to Pauli basis, call polarity_coordinates_from_L."""
    T = fw.pauli._vec_to_pauli_basis_transform(N)
    L_pauli = (T.conj().T @ L_vec @ T) / (2 ** N)
    if sigma is None:
        sigma = 0.0
    return fw.polarity_coordinates_from_L(L_pauli, N, sigma)


def main():
    t_us, rhos, raw = load_trajectory()
    print(f"Loaded {len(rhos)} ρ snapshots at t_us = {list(t_us)}")
    print(f"backend = {raw['backend']}, path = {raw['path']}, job_id = {raw['job_id']}")
    print(f"T2_min_cal = {raw['t2_min_us_calibration']} μs, "
          f"γ_eff_cal = {raw['gamma_eff_per_us_calibration']:.6f} /μs")
    print()

    # Initial guesses (from T2 calibration)
    g_z_cal = raw['gamma_eff_per_us_calibration']  # ~0.00208 /μs
    g_t1_cal = g_z_cal * 0.5  # weaker initial guess
    j_zz_cal = 0.0005
    hy_cal = 0.001

    models_and_x0 = [
        ('pure_Z', [g_z_cal, g_z_cal]),
        ('Z_plus_T1', [g_z_cal, g_z_cal, g_t1_cal, g_t1_cal]),
        ('Z_plus_ZZ', [g_z_cal, g_z_cal, j_zz_cal]),
        ('Z_plus_T1_plus_ZZ', [g_z_cal, g_z_cal, g_t1_cal, g_t1_cal, j_zz_cal]),
        ('Z_plus_hy', [g_z_cal, g_z_cal, hy_cal]),
    ]

    print(f"{'Model':<22} {'fit RMS':>12} {'in F112 scope':>15} {'F112 asym':>15} {'F112 rel asym':>15}  fitted params")
    print('-' * 130)

    results = {}
    for model, x0 in models_and_x0:
        x_fit, fit_loss = fit_model(model, t_us, rhos, x0)
        rms = float(np.sqrt(fit_loss / max(len(t_us) - 1, 1)))

        L_vec, c_list, g_list, c_kind, H = build_L_model(x_fit, model)
        sigma = float(np.real(sum(g_list)))
        pol = run_polarity_on_L(L_vec, N=2, sigma=sigma)
        m_sq = pol['norm_sq']['M']
        asym = pol['asymmetry']
        rel = abs(asym) / max(m_sq, 1e-15)

        # F112 scope: Hermitian H AND every c bit_b-homogeneous
        h_is_hermitian = np.allclose(H, H.conj().T)
        all_c_bit_b_homog = all(is_bit_b_homogeneous_pauli_label(k) for k in c_kind)
        in_scope = h_is_hermitian and all_c_bit_b_homog

        results[model] = {
            'fit_rms': rms,
            'fit_loss': fit_loss,
            'params': x_fit.tolist(),
            'param_kinds': c_kind + (['ZZ'] if 'ZZ' in model else []) + (['hy'] if 'hy' in model else []),
            'in_F112_scope': in_scope,
            'asymmetry': asym,
            'rel_asymmetry': rel,
            'M_norm_sq': m_sq,
        }

        scope_str = 'YES' if in_scope else 'no'
        print(f"{model:<22} {rms:>12.6f} {scope_str:>15} {asym:>+15.6e} {rel:>15.4e}  {[f'{p:.5f}' for p in x_fit]}")

    print()
    print("=" * 130)
    print("Interpretation")
    print("=" * 130)

    # Find best fit
    best_model = min(results, key=lambda m: results[m]['fit_rms'])
    print(f"\nBest-fit model: {best_model}  (RMS = {results[best_model]['fit_rms']:.6f})")
    print(f"  In F112 typed scope (Hermitian H + bit_b-homogeneous c): "
          f"{'YES' if results[best_model]['in_F112_scope'] else 'NO'}")
    print(f"  F112 polarity asymmetry on fitted L: {results[best_model]['asymmetry']:+.6e} "
          f"(rel {results[best_model]['rel_asymmetry']:.4e})")

    # Compare pure-Z baseline to best
    pure_rms = results['pure_Z']['fit_rms']
    best_rms = results[best_model]['fit_rms']
    improvement = (pure_rms - best_rms) / max(pure_rms, 1e-15)
    print(f"\nFit improvement over pure-Z baseline: {improvement * 100:.2f}% "
          f"({pure_rms:.6f} → {best_rms:.6f})")

    print()
    print("F112 reading per model:")
    for model, r in results.items():
        verdict_scope = 'in scope (Tier1Derived asymmetry = 0)' if r['in_F112_scope'] else 'outside typed scope'
        asym_verdict = 'BALANCED bit-exact' if r['rel_asymmetry'] < 1e-10 else (
            f"asymmetry rel = {r['rel_asymmetry']:.2e} (BROKEN)"
        )
        print(f"  {model:<22}: {verdict_scope}; observed F112 {asym_verdict}")

    print()
    print("Note: F112 says HERMITIAN H + EACH c_k bit_b-homogeneous → asymmetry = 0 bit-exact.")
    print("      σ⁻ amplitude damping (T1) c = σ⁻ = (X − iY)/2 has bit_b ∈ {0, 1} (mixed).")
    print("      All Z-only and Z+ZZ models are bit_b-homogeneous on c side.")
    print("      Single-site h_y · Y_l Hamiltonian is Hermitian; Y has bit_b=1 → H is bit_b-homogeneous.")


if __name__ == '__main__':
    main()
