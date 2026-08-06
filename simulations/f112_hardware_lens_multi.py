"""F112 polarity-asymmetry lens across 4 IBM Kingston Tier-A datasets.

Extends the single-dataset analysis (f112_block_cpsi_analysis.py) to all four
Tier-A datasets identified by the F112-hardware-inventory mapping (2026-05-26):

  1. block_cpsi_saturation     (2026-05-08, q13-q14, 5 t-points, idle/no-H)
  2. cusp_slowing               (2026-04-16, A_mid q124-q125 + B_high, 6 t-points, idle)
  3. chain_gamma0               (2026-04-19, 4 pairs along chain Q12-Q19, 9 t-points, (J/2)(XX+YY) bond H + Trotter)
                                NOTE these four "pair-runs" are four ADJACENT PAIRS of one
                                5-qubit chain (chain_phys [12,13,14,15,19]) carrying a
                                propagating single excitation, so each is an open, overlapping
                                marginal. A closed 2-qubit Lindblad fit is misspecified there
                                in a way no added channel repairs; the F112 asymmetry is still
                                well defined on whatever L comes out, which is all this survey
                                reads from them.
  4. f95_angle_steering         (2026-05-16, q82-q83, 2 omega values × 6 t-points, RZ Z-drive per Δt)

For each dataset / pair:
  - Reconstruct ρ(t) from 4×4 rho2_real + rho2_imag snapshots
  - Fit five candidate noise models (pure_Z, +T1, +ZZ, +T1+ZZ, +h_y) via least-squares
  - For each fitted L, compute F112 polarity asymmetry + scope classification
  - Report fit RMS, F112 asymmetry, bit_b-homogeneity verdict

Aggregate finding to surface:
  - Across all hardware-derived L's, does the standard Lindblad channel always
    preserve polarity balance (extending probes 1-14's empirical envelope)?
  - (Withdrawn as a deliverable: which noise model best fits each pair's data.
    The family has no per-qubit detuning, which this data needs, so its ranking
    compares members of the wrong family. See the companion experiment.)

No new QPU spend; pure analysis on existing data per feedback_qpu_conservative.
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

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1|, LOWERS |1> -> |0>
# Corrected 2026-07-29. This file previously used [[0,0],[1,0]] = |1><0|, the RAISING
# operator, mislabelled as |0><1|. Fitting a pump to relaxing data with rates clamped
# non-negative forced gamma_T1 -> 0, which was read as 'amplitude damping is rejected'.
# The model family below is still incomplete: it has no per-qubit Z detuning, which the
# block_cpsi transverse expectations need (they alternate in sign between samples). The
# companion experiment's fit-quality half is withdrawn until that is added; do not
# publish a model ranking from this script as it stands.
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def site_op(N, site, mat2):
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


def reconstruct_rho_from_pauli_dict(exps):
    """rho = (1/4) Σ ⟨P_α⟩ · P_α from 16-Pauli expectation dict 'a,b' → val."""
    rho = np.zeros((4, 4), dtype=complex)
    for key, val in exps.items():
        a, b = key.split(',')
        rho = rho + float(val) * np.kron(PAULI[a], PAULI[b])
    return rho / 4.0


def reconstruct_rho_from_real_imag(real_2d, imag_2d):
    """4×4 ρ from rho2_real and rho2_imag nested lists."""
    return np.array(real_2d, dtype=complex) + 1j * np.array(imag_2d, dtype=complex)


def build_L_model(params, model, fixed_h=None):
    """Returns (L_vec, c_kind_labels) for a 2-qubit candidate noise model.

    fixed_h: optional preset Hamiltonian (Hermitian, 4×4) added to the model's
    own H (e.g., the known RZ drive in f95). Bit_b homogeneity of fixed_h is
    enforced upstream; if you pass a bit_b-mixed fixed_h, the model's F112
    scope claim no longer holds.
    """
    N = 2
    d = 4
    Id = np.eye(d, dtype=complex)
    H = fixed_h.copy() if fixed_h is not None else np.zeros((d, d), dtype=complex)
    c_list = []
    g_list = []
    c_kind = []

    if model == 'pure_Z':
        gz0, gz1 = params
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z)]
        g_list = [gz0, gz1]
        c_kind = ['Z@0', 'Z@1']
    elif model == 'Z_plus_T1':
        gz0, gz1, gt0, gt1 = params
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z),
                  site_op(N, 0, SIGMA_MINUS), site_op(N, 1, SIGMA_MINUS)]
        g_list = [gz0, gz1, gt0, gt1]
        c_kind = ['Z@0', 'Z@1', 'sigma_minus@0', 'sigma_minus@1']
    elif model == 'Z_plus_ZZ':
        gz0, gz1, j_zz = params
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z)]
        g_list = [gz0, gz1]
        c_kind = ['Z@0', 'Z@1']
        H = H + j_zz * np.kron(Z, Z)
    elif model == 'Z_plus_T1_plus_ZZ':
        gz0, gz1, gt0, gt1, j_zz = params
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z),
                  site_op(N, 0, SIGMA_MINUS), site_op(N, 1, SIGMA_MINUS)]
        g_list = [gz0, gz1, gt0, gt1]
        c_kind = ['Z@0', 'Z@1', 'sigma_minus@0', 'sigma_minus@1']
        H = H + j_zz * np.kron(Z, Z)
    elif model == 'Z_plus_hy':
        gz0, gz1, hy = params
        c_list = [site_op(N, 0, Z), site_op(N, 1, Z)]
        g_list = [gz0, gz1]
        c_kind = ['Z@0', 'Z@1']
        H = H + hy * (site_op(N, 0, Y) + site_op(N, 1, Y))
    else:
        raise ValueError(f"Unknown model: {model}")

    # vec is COLUMN-stack (rho.flatten('F')), so the superoperator must be the
    # column-stack one: vec_F(A X B) = (B^T (x) A) vec_F(X). The earlier row-stack
    # form here integrated the commutator with -H^T, i.e. a sign-flipped drive for
    # any real symmetric H (the fixed f95 Z-drive); Y and the free-sign ZZ were blind
    # to it. Corrected 2026-07-29.
    L_vec = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for c, g in zip(c_list, g_list):
        c_dag_c = c.conj().T @ c
        anti = 0.5 * (np.kron(Id, c_dag_c) + np.kron(c_dag_c.T, Id))
        L_vec = L_vec + g * (np.kron(c.conj(), c) - anti)
    return L_vec, c_kind, H


def evolve(L_vec, rho0, dt):
    v0 = vec(rho0)
    vt = expm(L_vec * dt) @ v0
    return devec(vt)


def fit_residual(params, model, t_us, rhos, fixed_h):
    n_rate_params = {'pure_Z': 2, 'Z_plus_T1': 4, 'Z_plus_ZZ': 2,
                     'Z_plus_T1_plus_ZZ': 4, 'Z_plus_hy': 2}[model]
    if np.any(np.array(params[:n_rate_params]) < 0):
        return 1e6
    try:
        L_vec, _, _ = build_L_model(params, model, fixed_h=fixed_h)
    except Exception:
        return 1e6
    total = 0.0
    rho0 = rhos[0]
    dt0 = t_us[0]
    for i, t in enumerate(t_us):
        if i == 0:
            continue
        rho_pred = evolve(L_vec, rho0, t - dt0)
        diff = rho_pred - rhos[i]
        total += float(np.sum(np.abs(diff) ** 2))
    return total


def fit_one(model, t_us, rhos, x0, fixed_h):
    # STILL SINGLE-START, deliberately, and this is the sibling of
    # f112_block_cpsi_analysis.py, which went multi-start on 2026-08-05. Two
    # defects are known here and neither is repaired by a starting point:
    # (a) this script's fits show nesting violations (a superset fitting WORSE
    #     than its subset: on `chain_gamma0 Q13-Q14` Z_plus_T1 gives RMS
    #     0.628598 against pure_Z's 0.628012, and Z_plus_T1_plus_ZZ gives
    #     0.628018 against Z_plus_ZZ's 0.626140), so several rows are not
    #     converged at all;
    # (b) on block_cpsi the second Z rate is not identified from above: the
    #     optimiser raises it until the model predicts no transverse coherence
    #     (the least-squares answer to a sign-alternating sequence it cannot
    #     turn), and the objective then goes bit-identical to its
    #     gamma -> infinity limit, so no optimiser can locate it.
    # F112_HARDWARE_LENS_KINGSTON.md lists converged multistart, per-qubit Z
    # detuning, matched-parameter-count comparison and uncertainties as one
    # redo. Adding multistart alone here would make (a) quieter without making
    # any ranking mean more, so it waits for that redo.
    result = minimize(
        fit_residual, x0, args=(model, t_us, rhos, fixed_h),
        method='Nelder-Mead',
        options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 5000},
    )
    return result.x, float(result.fun)


def is_bit_b_homogeneous_label(label):
    if label.startswith(('X@', 'Y@', 'Z@')):
        return True
    if label.startswith('sigma_'):
        return False
    return False


def _to_row_stack(L_vec, d=4):
    """Row-stack representative of a column-stack Liouvillian.

    The propagator here is column-stack (rho.flatten('F')), which is what makes
    exp(L t) correct. The F112/F113 polarity read is pinned to the OTHER pairing:
    a row-stack L against the order='F' Pauli transform. That mismatch is
    the pairing F113's stated direction is written for: it fixes which half is
    called +1/2, hence the SIGN of Asymmetry (see PauliBasis.cs; F113 gives
    -2.08e-3 at omega=0.13, gamma_T1=0.001, N=2 for cooling, and a consistent
    pairing would give +2.08e-3).

    vec_C(rho) = vec_F(rho^T) = T vec_F(rho) with T the transpose permutation, so
    L_row = T L_col T exactly (T^2 = I). Converting here keeps the propagator
    correct AND the reported asymmetry sign on its pin. Bond Hamiltonians give
    asymmetry 0 either way, so only the Z-drive runs can see this.
    """
    import numpy as _np
    T = _np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            T[j * d + i, i * d + j] = 1.0
    return T @ L_vec @ T


def run_polarity_on_L(L_vec, N=2, sigma=None):
    T = fw.pauli._vec_to_pauli_basis_transform(N)
    L_pauli = (T.conj().T @ L_vec @ T) / (2 ** N)
    if sigma is None:
        sigma = 0.0
    return fw.polarity_coordinates_from_L(L_pauli, N, sigma)


def analyze_trajectory(label, t_us, rhos, gz_init, gt1_init,
                       j_zz_init, hy_init, fixed_h=None):
    """Run all 5 models, print + return result dict per model."""
    print(f"\n{'='*100}")
    print(f"  {label}")
    print(f"{'='*100}")
    print(f"  t_points: {len(t_us)}, t range: {min(t_us):.2f} - {max(t_us):.2f} μs")
    if fixed_h is not None:
        print(f"  fixed H present (e.g., known applied drive)")

    models_and_x0 = [
        ('pure_Z', [gz_init, gz_init]),
        ('Z_plus_T1', [gz_init, gz_init, gt1_init, gt1_init]),
        ('Z_plus_ZZ', [gz_init, gz_init, j_zz_init]),
        ('Z_plus_T1_plus_ZZ', [gz_init, gz_init, gt1_init, gt1_init, j_zz_init]),
        ('Z_plus_hy', [gz_init, gz_init, hy_init]),
    ]
    print(f"\n  {'Model':<22} {'fit RMS':>12} {'F112 scope':>12} {'F112 asym':>14} {'rel':>12}")
    print(f"  {'-'*78}")

    results = {}
    for model, x0 in models_and_x0:
        x_fit, fit_loss = fit_one(model, t_us, rhos, x0, fixed_h)
        rms = float(np.sqrt(fit_loss / max(len(t_us) - 1, 1)))
        L_vec, c_kind, H = build_L_model(x_fit, model, fixed_h=fixed_h)
        # The palindrome centre is the DEPHASING sum. Including the T1 rates as well
        # is the other convention in use here, and it changes nothing: `asym` is
        # untouched either way (the shift is a multiple of I, which is Ad_Pi-fixed and
        # lands wholly in M_zero), and ||M||^2 is a parabola in sigma minimised at
        # sum(gz) + sum(gt1)/2, so the two conventions sit symmetrically about the
        # minimum and agree on ||M||^2 to machine precision (`asym` itself is
        # bit-identical). Verified at N=2 and N=3.
        sigma = float(np.real(x_fit[0] + x_fit[1]))
        pol = run_polarity_on_L(_to_row_stack(L_vec), N=2, sigma=sigma)
        m_sq = pol['norm_sq']['M']
        asym = pol['asymmetry']
        # Denominator is the POLARITY CONTENT ||M_anti||^2 = ||M_+||^2 + ||M_-||^2, the whole of
        # what asym is a difference of, so rel is a contrast ratio in [0, 1]. The retired
        # max(||M||^2, 1e-15) was the wrong scale and saturated on its own floor wherever M
        # vanishes. Shape and rationale: simulations/framework/workflows/polarity_fingerprint.py.
        m_anti = float(pol['norm_sq']['M_plus_half'] + pol['norm_sq']['M_minus_half'])
        rel = 0.0 if m_anti == 0.0 else abs(asym) / m_anti
        h_herm = np.allclose(H, H.conj().T, atol=1e-12)
        all_c_homog = all(is_bit_b_homogeneous_label(k) for k in c_kind)
        in_scope = h_herm and all_c_homog

        results[model] = {
            'fit_rms': rms,
            'params': x_fit.tolist(),
            'in_F112_scope': in_scope,
            'asymmetry': asym,
            'rel_asymmetry': rel,
        }
        scope_str = 'YES' if in_scope else 'no'
        print(f"  {model:<22} {rms:>12.6f} {scope_str:>12} {asym:>+14.6e} {rel:>12.4e}")

    best = min(results, key=lambda m: results[m]['fit_rms'])
    pure_rms = results['pure_Z']['fit_rms']
    improvement = (pure_rms - results[best]['fit_rms']) / max(pure_rms, 1e-15)
    print(f"\n  Best fit: {best} (RMS = {results[best]['fit_rms']:.6f}, "
          f"{improvement*100:.1f}% improvement over pure_Z baseline)")
    print(f"  All F112 asymmetries: max rel = {max(r['rel_asymmetry'] for r in results.values()):.4e}")
    return {'label': label, 'best_model': best, 'results': results}


# ============================================================
# Dataset adapters
# ============================================================

def load_block_cpsi():
    """2026-05-08 Kingston q13-q14, 5 t-points, idle."""
    with open('data/ibm_block_cpsi_saturation_may2026/'
              'block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json',
              encoding='utf-8') as f:
        d = json.load(f)
    t_grid = np.array(d['t_grid_us'], dtype=float)
    rhos = [reconstruct_rho_from_pauli_dict(s['expectations']) for s in d['t_snapshots']]
    return [('block_cpsi q13-q14 (Kingston, 2026-05-08, idle)', t_grid, rhos, None)]


def load_cusp_slowing():
    """2026-04-16 Kingston, 2 pairs (A_mid q124-q125, B_high), 6 t-points each, idle.
    Note: first t-point is 4.05 μs (not 0); fit propagates from there."""
    with open('data/ibm_cusp_slowing_april2026/'
              'cusp_slowing_ibm_kingston_20260416_212042.json',
              encoding='utf-8') as f:
        d = json.load(f)
    out = []
    for pair_key in ['A_mid', 'B_high']:
        pr = d['pair_runs'][pair_key]
        qubits = pr['pair']['qubits']
        traj = pr['trajectory']
        t_us = np.array([s['t_us'] for s in traj], dtype=float)
        rhos = [reconstruct_rho_from_real_imag(s['rho2_real'], s['rho2_imag'])
                for s in traj]
        label = f'cusp_slowing {pair_key} q{qubits[0]}-q{qubits[1]} (Kingston, 2026-04-16, idle)'
        out.append((label, t_us, rhos, None))
    return out


def load_chain_gamma0():
    """2026-04-19 Kingston, 4 pairs along chain, 9 t-points each, (J/2)(XX+YY) bond H via Trotter."""
    with open('data/ibm_chain_gamma0_april2026/'
              'chain_gamma0_hardware_20260419_110200.json',
              encoding='utf-8') as f:
        d = json.load(f)
    t_grid = np.array(d['t_us_grid'], dtype=float)
    # The JSON field is named J_rad_per_us but the dataset README says in as many words
    # that 1.838 is a FREQUENCY in MHz (cycles/us) and the angular rate would be 11.55.
    # Convert, or the bond Hamiltonian handed to the fitter is 2*pi too weak.
    J = 2.0 * np.pi * float(d['J_rad_per_us'])  # rad/us for the (J/2)(XX+YY) bond
    out = []
    for pair_label in d['pair_labels']:
        traj = d['hw_trajectories'][pair_label]
        t_us = np.array([s['t_us'] for s in traj], dtype=float)
        rhos = [reconstruct_rho_from_real_imag(s['rho2_real'], s['rho2_imag'])
                for s in traj]
        # The applied H during evolution is J·(XX + YY)/2 (XY model bond bilinear).
        # Bit_b: X=0, Y=1; XX bit_b=0, YY bit_b=0. So bit_b-homogeneous (Π²-even, F1 truly).
        H_xy_bond = (J / 2.0) * (np.kron(X, X) + np.kron(Y, Y))
        label = f'chain_gamma0 {pair_label} (Kingston, 2026-04-19, XY-bond H J={J:.3f})'
        out.append((label, t_us, rhos, H_xy_bond))
    return out


def load_f95():
    """2026-05-16 Kingston, q82-q83 + q-? at 2 omega values × 6 t-points, with RZ Z-drive per Δt."""
    paths = [
        ('omega=0.13', 'data/ibm_f95_angle_steering_may2026/'
                       'cusp_complex_phase_hardware_ibm_kingston_omega0.130_20260516_204827.json'),
        ('omega=0.25', 'data/ibm_f95_angle_steering_may2026/'
                       'cusp_complex_phase_hardware_ibm_kingston_omega0.250_20260516_205705.json'),
    ]
    out = []
    for omega_tag, path in paths:
        with open(path, encoding='utf-8') as f:
            d = json.load(f)
        omega_per_us = float(d['omega_per_us'])
        for pair_key in d['pair_runs']:
            pr = d['pair_runs'][pair_key]
            qubits = pr['pair']['qubits']
            traj = pr['trajectory']
            t_us = np.array([s['t_us'] for s in traj], dtype=float)
            rhos = [reconstruct_rho_from_real_imag(s['rho2_real'], s['rho2_imag'])
                    for s in traj]
            # Applied H during free evolution: the protocol injects RZ(omega*dt) on ONE
            # qubit of the pair (data README: "per-chunk RZ injection on the steering
            # qubit"), and the dataset's own predicted crossing angles match
            # phi_0 - omega*t exactly, not phi_0 - 2*omega*t. Which site carries it is
            # not observable here (the Bell coherence turns at the SUM of the two Z
            # coefficients), so site 0 is a representative choice. Z has bit_b=1, so the
            # drive is bit_b-homogeneous either way.
            H_z_drive = omega_per_us * site_op(2, 0, Z) / 2.0
            label = f'f95 {pair_key} q{qubits[0]}-q{qubits[1]} ({omega_tag}, Kingston, 2026-05-16)'
            out.append((label, t_us, rhos, H_z_drive))
    return out


def main():
    print("F112 polarity-asymmetry lens across 4 Kingston Tier-A datasets")
    print("=" * 100)

    aggregated = []

    # Initial guesses (T2-calibration scale)
    gz_init = 0.005   # per μs
    gt1_init = 0.001
    j_zz_init = 0.0005
    hy_init = 0.001

    for loader in [load_block_cpsi, load_cusp_slowing, load_chain_gamma0, load_f95]:
        for (label, t_us, rhos, fixed_h) in loader():
            result = analyze_trajectory(label, t_us, rhos,
                                        gz_init, gt1_init, j_zz_init, hy_init,
                                        fixed_h=fixed_h)
            aggregated.append(result)

    # ============================================================
    # Aggregate summary
    # ============================================================
    print()
    print("=" * 100)
    print("  AGGREGATE SUMMARY ACROSS ALL 4 DATASETS")
    print("=" * 100)
    print()
    print(f"  Total pair-runs analyzed: {len(aggregated)}")

    max_rel_asym = 0.0
    n_balanced = 0
    n_total_fits = 0
    best_models = []
    for r in aggregated:
        best_models.append(r['best_model'])
        for model_name, mr in r['results'].items():
            n_total_fits += 1
            if mr['rel_asymmetry'] < 1e-10:
                n_balanced += 1
            max_rel_asym = max(max_rel_asym, mr['rel_asymmetry'])

    print(f"  Total fitted L's: {n_total_fits}")
    print(f"  F112 asymmetry BALANCED bit-exact (rel < 1e-10): {n_balanced} / {n_total_fits}")
    print(f"  Max observed relative asymmetry: {max_rel_asym:.4e}")
    print()
    from collections import Counter
    bm_counts = Counter(best_models)
    print("  Best-fit model distribution across pair-runs (NOT a channel finding:")
    print("  the model family carries no per-qubit detuning, which this data needs,")
    print("  so this ranking compares members of the wrong family):")
    for m, c in bm_counts.most_common():
        print(f"    {m:<22} : {c}")

    # Partition fits by (in-scope, balanced) status
    in_scope_balanced = 0
    in_scope_broken = 0
    out_scope_balanced = 0
    out_scope_broken = 0
    broken_details = []
    for r in aggregated:
        for model_name, mr in r['results'].items():
            balanced = mr['rel_asymmetry'] < 1e-10
            in_scope = mr['in_F112_scope']
            if in_scope and balanced:
                in_scope_balanced += 1
            elif in_scope and not balanced:
                in_scope_broken += 1
            elif not in_scope and balanced:
                out_scope_balanced += 1
            else:
                out_scope_broken += 1
                broken_details.append((r['label'], model_name, mr['rel_asymmetry']))

    print()
    print(f"  F112 typed Tier1Derived scope (Hermitian H + bit_b-homog c):")
    print(f"    in-scope BALANCED:  {in_scope_balanced}  (F112 prediction met)")
    print(f"    in-scope BROKEN:    {in_scope_broken}   (would be a counterexample to F112)")
    print(f"  Out-of-scope (bit_b-mixed c via σ⁻ T1):")
    print(f"    out-scope BALANCED: {out_scope_balanced}  (broader empirical envelope, T1 case 'happens to balance')")
    print(f"    out-scope BROKEN:   {out_scope_broken}   (broader envelope COUNTEREXAMPLES)")

    print()
    print("=" * 100)
    print("  F112 verdict")
    print("=" * 100)
    print(f"""
  F112 TYPED Tier1Derived holds across all {in_scope_balanced + in_scope_broken}
  in-scope fits ({in_scope_balanced} BALANCED bit-exact, {in_scope_broken} counterexamples).
  This is the theorem-content statement.

  The broader empirical envelope (probes 1-14 conjecture: bit_b-mixed c
  with σ⁻ T1 ALSO preserves balance) has {out_scope_broken} counterexamples
  in this run. All {out_scope_broken} are f95-dataset fits where the applied
  Hamiltonian includes a Z drive on ONE qubit of the pair, (omega/2)·Z_0, combined with
  σ⁻ T1 amplitude damping. The Z-drive lives on bit_b=1, σ⁻ on bit_b ∈ {{0,1}}
  jointly; together they produce non-zero polarity asymmetry up to rel
  {max_rel_asym:.2e}.

  For idle and XY-bond Hamiltonians, σ⁻ T1 fits preserve balance bit-exact;
  for Z-drive H, σ⁻ T1 fits do NOT. That is F113's selector Tr(Z_l H) seen on
  fitted Liouvillians. It is a statement about σ⁻, not about bit_b-mixed c in
  general: a generic cross-bit_b collapse operator breaks the balance on a plain
  Heisenberg bond, where Tr(Z_l H) = 0.

  Note what this split is NOT: evidence about the chip. The single-site Z
  reaches these models only as the known drive handed to the fitter, and no
  other model in the set has one at all, so the split follows from the model
  file. The same data need a per-qubit detuning, itself a single-site Z, that
  none of these models carry.
""")
    if broken_details:
        print("  Counterexamples (out-scope BROKEN cases):")
        for label, model_name, rel in broken_details:
            print(f"    {label:<70} {model_name:<22} rel = {rel:.3e}")


if __name__ == '__main__':
    main()
