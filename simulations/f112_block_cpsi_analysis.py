"""F112 polarity-asymmetry lens on ibm_kingston block-CPsi-saturation trajectory.

Dataset: `data/ibm_block_cpsi_saturation_may2026/block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json`
- Initial state: (|D_0⟩ + |D_1⟩)/√2 = (|00⟩ + (|01⟩+|10⟩)/√2)/√2 on Kingston qubits 13, 14
- 5 t-points: 0, 120, 240, 360, 480 μs
- 16 Pauli expectations per t-point (full 2-qubit tomography)
- T2_min calibration: 480 μs; the D[Z] rate reproducing it is
  γ_eff = 1/(2·T2) ≈ 0.00104 per μs. (The flown JSON stores 0.00208, which is
  1/T2, the superseded convention. It is a record and is not edited; see the
  note on g_z_cal in main().)
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

  The diagnostic: where does a fitted model sit relative to F112's typed scope?
  NOT which model fits best: the family has no per-qubit detuning, which this
  dataset needs, so its ranking compares members of the wrong family.

Fitting (repaired 2026-08-05):
  The fit is MULTI-START (`seed_set`), because single-start was seed-dependent:
  `Z_plus_T1` reached RMS 0.274109 from the 1/T2 seed and 0.307329 from the
  1/(2*T2) seed, a different local minimum and a worse fit. Both conventions are
  now seeds among many, so which one a caller passes cannot decide the answer.
  Multi-start changes no published RMS: all five models land on the same minima.

  It does NOT fix, and cannot fix, the second Z rate, which this data does not
  identify from above at all. `identifiability_profile` reports that beside the
  table rather than hiding it in a converged-looking number.

Output:
  Table per model: fit_RMS, F112 asymmetry, in_F112_scope flag.
  Identifiability profile per parameter: relative objective change under a
  factor-of-ten rescaling, up and down.
  Interpretation: what F112 says about where the fitted models sit. The 1.72×
  anomaly is NOT explained here; see the companion experiment for why.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# This script prints rho, gamma and mu; the Windows console default is cp1252
# and raises on the first one. Same guard as block_cpsi_run_planner_2026_05_08.
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

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
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1|, LOWERS |1> -> |0>
# Corrected 2026-07-29. This file previously used [[0,0],[1,0]] = |1><0|, the RAISING
# operator, mislabelled as |0><1|. Fitting a pump to relaxing data with rates clamped
# non-negative forced gamma_T1 -> 0, which was read as 'amplitude damping is rejected'.
# The model family below is still incomplete: it has no per-qubit Z detuning, which this
# dataset needs (qubit 14's transverse expectations alternate in sign between samples).
# The companion experiment's fit-quality half is withdrawn until that is added; do not
# publish a model ranking from this script as it stands.
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

DATA_PATH = Path('data/ibm_block_cpsi_saturation_may2026/'
                 'block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json')

# The flown calibration. Asserted against the JSON in main(), so this constant
# cannot silently drift away from the record it describes. The house D[Z] rate
# is 1/(2*T2); the JSON stores 1/T2 under the superseded convention. Both are
# used as multi-start seeds, which is why neither is privileged here.
T2_MIN_US = 480.0


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


def _fit_once(model, t_us, rhos, x0):
    """Local minimization from a single initial guess x0."""
    result = minimize(
        fit_residual, x0, args=(model, t_us, rhos),
        method='Nelder-Mead',
        options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 5000},
    )
    return result.x, float(result.fun)


def seed_set(model, n_random=40, rng_seed=20260805):
    """Deterministic multi-start seeds for `model`.

    Both T2 -> gamma conventions are included as seeds (1/T2 and 1/(2*T2)), so
    the fit no longer depends on which one a caller passes: they are two members
    of one set. The rest are log-uniform over 1e-4 .. 0.5 per us.

    That range does NOT bracket everything, and this says so rather than claiming
    a coverage it does not have: `Z_plus_T1_plus_ZZ` wins at a second Z rate of
    13.59, 27x above the range's top. It gets there because the simplex is free
    to leave its starting box, and because that parameter sits on a flat plateau
    (see `identifiability_profile`) where any large value is as good as any
    other. Basin counts REACHABLE FROM THIS SEED SET, measured by running it
    (losses distinct at 6 dp): 2, 6, 3, 21 and 3 for the five models in the order
    used in main(). That is a count of what these seeds find, not a claim about
    how many basins exist.
    """
    n_rate = {'pure_Z': 2, 'Z_plus_T1': 4, 'Z_plus_ZZ': 2,
              'Z_plus_T1_plus_ZZ': 4, 'Z_plus_hy': 2}[model]
    n_par = {'pure_Z': 2, 'Z_plus_T1': 4, 'Z_plus_ZZ': 3,
             'Z_plus_T1_plus_ZZ': 5, 'Z_plus_hy': 3}[model]

    seeds = []
    # Written as 2/T2 and 1/T2 when first landed, which is neither pair: the two
    # conventions are 1/T2 (old book) and 1/(2*T2) (house). Corrected 2026-08-05.
    # The factor of two slipped into the repair that exists because of it.
    for g in (1.0 / T2_MIN_US, 1.0 / (2.0 * T2_MIN_US)):   # old book, house book
        rates = [g, g] + ([0.5 * g, 0.5 * g] if 'T1' in model else [])
        extra = [0.0005] * ('ZZ' in model) + [0.001] * ('hy' in model)
        seeds.append(np.array(rates + extra, dtype=float))

    rng = np.random.default_rng(rng_seed)
    for _ in range(n_random):
        rates = 10 ** rng.uniform(-4, -0.3, size=n_rate)
        extra = rng.uniform(-0.01, 0.01, size=n_par - n_rate)
        seeds.append(np.concatenate([rates, extra]))
    return seeds


def fit_model(model, t_us, rhos, x0=None):
    """Multi-start minimization; returns the best (x, loss) over `seed_set`.

    Single-start was seed-dependent and quietly so: `Z_plus_T1` reaches RMS
    0.274109 from one convention's seed and 0.307329 from the other's, a
    different local minimum and a worse fit. Multi-start removes the dependence
    on WHICH seed; it does not remove a flat direction, which is what
    `identifiability_profile` below reports separately.
    """
    seeds = seed_set(model)
    if x0 is not None:
        seeds = [np.asarray(x0, dtype=float)] + seeds
    best_x, best_loss = None, np.inf
    for s in seeds:
        x, loss = _fit_once(model, t_us, rhos, s)
        if loss < best_loss:
            best_x, best_loss = x, loss
    return best_x, best_loss


def identifiability_profile(x, model, t_us, rhos, factor=10.0):
    """Per-parameter relative loss change when the parameter is scaled.

    Returned per index: (up, down), the RELATIVE objective CHANGE on scaling that
    parameter by `factor` and by 1/`factor`, the other parameters held. The change
    can be negative when the fit sits just below a flat region (`Z_plus_T1`'s
    second Z rate prints -1.8e-16), so this is a change, not an increase.

    This is READ, not gated, and deliberately so. In FOUR of the five models here
    (`Z_plus_hy` is the exception, below) the second Z rate is not identified from
    above. Against a sign-alternating sequence the least-squares optimum for a
    purely decaying model is zero, so the optimiser raises that rate until the
    model predicts no transverse coherence at all: it is switching a prediction
    off, not fitting a rate. The data itself says the coherence is alive there
    (57% of its t=0 magnitude at the first 120 us delay, decaying log-linearly
    with T2 about 170 us). `Z_plus_hy` is the exception: its second Z rate sits at
    an interior minimum near 0.0152, with gamma -> infinity about 6% WORSE. Do not
    read that as the Y-field representing the alternation. It does not turn the
    transverse coherence at all (it tips X toward Z), and the interior minimum's
    whole advantage sits in the slot-1 LONGITUDINAL Paulis: +0.0277 of the +0.0275
    total, while on the transverse block the fit is 0.0002 worse than the
    plateau's exact zero. Two tidy causal stories for this exception were written
    and withdrawn under review on 2026-08-05; it is left as a decomposition.

    Where the objective goes exactly flat is ARITHMETIC, not physics: it is where
    the model's own coherence factor exp(-2*g*120) falls under the ULP of the
    objective, near 1e-16. No value is quoted for that onset because it does not
    carry digits: it is model-dependent, and at fixed model it moves by 8% under a
    relative 1e-8 jitter of the OTHER rate. In higher precision it moves again. So
    a binary "is it bit-identical" test would report a property of float64, and a
    threshold would only move the arbitrariness elsewhere. The precision-free
    statement is the flattening itself: at 0.05 per us the objective is already
    within 2e-6 relative of its limit.

    The numbers below let a reader see that directly: `up` at or near zero means
    the fitted digits are where the simplex stopped, and the value has stopped
    denoting a rate. It is NOT a bound on that qubit's T2.

    One blind spot, since the probe is multiplicative: a parameter fitted at
    exactly 0 returns (0.0, 0.0) and so wears the plateau's signature without
    being on a plateau. No fitted parameter here is 0; check before reading the
    legend if that ever changes.
    """
    base = fit_residual(x, model, t_us, rhos)
    prof = []
    for i in range(len(x)):
        vals = []
        for f in (factor, 1.0 / factor):
            probe = np.array(x, dtype=float)
            probe[i] = probe[i] * f
            vals.append((fit_residual(probe, model, t_us, rhos) - base) / base)
        prof.append(tuple(vals))
    return prof


def is_bit_b_homogeneous_pauli_label(label):
    """For c labelled as 'X@k', 'Y@k', 'Z@k', 'sigma_minus@k', etc., return
    True iff the operator is bit_b-homogeneous as a Pauli sum.

    bit_b = (#Y + #Z) mod 2 per Pauli string. Single-Pauli letters:
    X bit_b = 0; Y, Z bit_b = 1; I bit_b = 0.
    sigma_minus = (X + iY)/2 = [[0,1],[0,0]]: mixed bit_b (X is 0, Y is 1) → False.
    sigma_plus = (X - iY)/2 = [[0,0],[1,0]]: same mixture, same verdict.
    """
    if label.startswith('X@') or label.startswith('Y@') or label.startswith('Z@'):
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

    # Starting points: NONE is privileged any more. `fit_model` runs a
    # deterministic multi-start (see `seed_set`), because single-start was
    # seed-dependent: `Z_plus_T1` reached RMS 0.274109 from the 1/T2 seed and
    # 0.307329 from the 1/(2*T2) seed, a different local minimum and a worse
    # fit. Which T2 -> gamma convention seeds the fit is therefore no longer a
    # question this script can get wrong; both conventions are seeds.
    #
    # What multi-start does NOT fix, and what `identifiability_profile` reports
    # beside the table: in four of the five models the second Z rate is not
    # identified from above. The optimiser raises it until the model predicts no
    # transverse coherence, which is the least-squares answer to a
    # sign-alternating sequence it cannot turn; the objective then goes
    # bit-identical to its gamma -> infinity limit, at a threshold set by
    # double-precision epsilon rather than by the data. A value printed inside
    # that plateau is a stopping point, not a measurement, and it is NOT a bound
    # on that qubit's T2: the same data gives that qubit T2 about 170 us from the
    # decay of its coherence MAGNITUDE. `Z_plus_hy` has no plateau at all, and the
    # docstring of `identifiability_profile` says what that does and does NOT
    # mean. See docs/GLOSSARY.md, "The T2 -> gamma conversion", for the convention
    # this whole thread hangs on.
    assert raw['t2_min_us_calibration'] == T2_MIN_US, (
        f"calibration record moved: JSON says {raw['t2_min_us_calibration']}, "
        f"this script's constant says {T2_MIN_US}")

    models = ['pure_Z', 'Z_plus_T1', 'Z_plus_ZZ', 'Z_plus_T1_plus_ZZ', 'Z_plus_hy']

    print(f"{'Model':<22} {'fit RMS':>12} {'in F112 scope':>15} {'F112 asym':>15} {'F112 rel asym':>15}  fitted params")
    print('-' * 130)

    results = {}
    for model in models:
        x_fit, fit_loss = fit_model(model, t_us, rhos)
        rms = float(np.sqrt(fit_loss / max(len(t_us) - 1, 1)))

        L_vec, c_list, g_list, c_kind, H = build_L_model(x_fit, model)
        # The palindrome centre is the DEPHASING sum, so take the Z entries only.
        # Including the T1 rates changes nothing measurable: `asym` is unchanged (the
        # shift is a multiple of I, Ad_Pi-fixed, landing wholly in M_zero), and
        # ||M||^2 is a parabola in sigma minimised at sum(gz) + sum(gt1)/2, so the two
        # conventions sit symmetrically about the minimum, agreeing on ||M||^2 to
        # machine precision (`asym` itself is bit-identical).
        sigma = float(np.real(sum(g for g, k in zip(g_list, c_kind)
                                  if k.startswith('Z@'))))
        pol = run_polarity_on_L(_to_row_stack(L_vec), N=2, sigma=sigma)
        m_sq = pol['norm_sq']['M']
        asym = pol['asymmetry']
        rel = abs(asym) / max(m_sq, 1e-15)

        # F112 scope: Hermitian H AND every c bit_b-homogeneous
        h_is_hermitian = np.allclose(H, H.conj().T)
        all_c_bit_b_homog = all(is_bit_b_homogeneous_pauli_label(k) for k in c_kind)
        in_scope = h_is_hermitian and all_c_bit_b_homog

        prof = identifiability_profile(x_fit, model, t_us, rhos)

        results[model] = {
            'identifiability': prof,
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
    print("Identifiability: relative change of the objective when one parameter is scaled")
    print("=" * 130)
    print("A parameter whose x10 column is ~0 is NOT identified from above: the fitted digits")
    print("are where the simplex stopped on a flat plateau, and the data carries a lower bound")
    print("with no upper one. Read the numbers; there is no threshold here.")
    print()
    for model in models:
        kinds = results[model]['param_kinds']
        print(f"  {model}")
        for i, (up, down) in enumerate(results[model]['identifiability']):
            name = kinds[i] if i < len(kinds) else f'p{i}'
            print(f"    {name:<18} value {results[model]['params'][i]:>12.5f}   "
                  f"x10: {up:>12.3e}   /10: {down:>12.3e}")

    print()
    print("=" * 130)
    print("Interpretation")
    print("=" * 130)

    # Find best fit
    best_model = min(results, key=lambda m: results[m]['fit_rms'])
    print("\nNOTE: the ranking below is NOT a channel finding. This model family has no")
    print("      per-qubit Z detuning, which this dataset needs (qubit 14's transverse")
    print("      expectations alternate in sign between samples), so it compares members")
    print("      of the wrong family.")
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
    print("      σ⁻ amplitude damping (T1) c = σ⁻ = (X + iY)/2 = [[0,1],[0,0]] has bit_b ∈ {0, 1} (mixed).")
    print("      All Z-only and Z+ZZ models are bit_b-homogeneous on c side.")
    print("      Single-site h_y · Y_l Hamiltonian is Hermitian; Y has bit_b=1 → H is bit_b-homogeneous.")


if __name__ == '__main__':
    main()
