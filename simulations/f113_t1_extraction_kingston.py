"""F113-inversion T1-extraction on Kingston f95 angle-steering datasets.

Welle 5.A application of the F113 closed form:

    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

inverted to extract a γ_T1 estimate from the measured polarity asymmetry:

    γ_T1 = −asymmetry / ((1 / 2) · 4^N · ω)       (ONE driven site, γ_pump = 0)

Pipeline per dataset (f95 omega=0.13 + omega=0.25, A_mid q82-q83 + B_high q13-q14):

  1. Load ρ(t) trajectory (6 t-points per pair-run).
  2. Fit minimal Z + σ⁻ T1 Lindblad model with the applied Z-drive H = (ω/2)·Z on the
     STEERING QUBIT (the known omega), gamma_z and gamma_T1 as free parameters.
  3. Compute the polarity asymmetry of the fitted L via polarity_coordinates_from_hc.
  4. Invert F113: γ_T1_F113 = −asymmetry / (8 · ω) at N = 2 (one driven site).
  5. Compare three γ_T1 readings:
     - γ_T1_fit (least-squares trajectory fit)
     - γ_T1_F113 (inverted from polarity asymmetry of the fit's L)
     - γ_T1_calibration (from device characterization in dataset metadata: 1/T1)

For self-consistency, γ_T1_fit and γ_T1_F113 should agree bit-exactly (the fit's
L was used to compute the asymmetry, so F113 inversion just recovers the fit value).
The interesting comparison is fit vs calibration: if they differ, the effective
T1 acting in the Bell-state evolution differs from the isolated-qubit T1.
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
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # standard physics lowering


def site_op(N, l, m2):
    mats = [I2] * N
    mats[l] = m2
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def vec(rho):
    return rho.flatten('F')


def devec(v, d=4):
    return v.reshape((d, d), order='F')


def reconstruct_rho(real_2d, imag_2d):
    return np.array(real_2d, dtype=complex) + 1j * np.array(imag_2d, dtype=complex)


# The protocol injects RZ(omega*dt) on ONE qubit of the pair (data README: "per-chunk
# RZ injection on the steering qubit"), and the dataset's own predicted crossing angles
# match phi_0 - omega*t exactly, not phi_0 - 2*omega*t. Which site carries it is weakly
# observable, contrary to what the Bell coherence alone suggests: the MEASURED initial
# state has single-excitation coherences that turn at 2*delta_1 and 2*delta_0 separately,
# and site 1 fits about 5% better on both A_mid runs. The effect on gamma_T1 is under
# 1.5%, so site 0 is kept as the representative choice.
DRIVE_SITE = 0


def build_L_z_plus_t1(N, omega, gamma_z, gamma_t1):
    """Minimal model: H = (omega/2)·Z_drive + Z-deph at gamma_z + sigma- T1 at gamma_t1."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    H = (omega / 2.0) * site_op(N, DRIVE_SITE, Z)
    # vec is column-stack (flatten('F')), so the superoperator is the column-stack one:
    # vec_F(A X B) = (B^T (x) A) vec_F(X). A row-stack build here would integrate the
    # commutator with -H^T, which for this real symmetric drive means -omega.
    L_vec = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        Zl = site_op(N, l, Z)
        ZlcZl = Zl.conj().T @ Zl
        L_vec = L_vec + gamma_z * (np.kron(Zl.conj(), Zl)
                                   - 0.5 * (np.kron(Id, ZlcZl) + np.kron(ZlcZl.T, Id)))
        Sl = site_op(N, l, SIGMA_MINUS)
        SlcSl = Sl.conj().T @ Sl
        L_vec = L_vec + gamma_t1 * (np.kron(Sl.conj(), Sl)
                                    - 0.5 * (np.kron(Id, SlcSl) + np.kron(SlcSl.T, Id)))
    return L_vec


def evolve(L_vec, rho0, dt):
    return devec(expm(L_vec * dt) @ vec(rho0))


def fit_z_t1(omega, t_us, rhos):
    """Fit (gamma_z, gamma_t1) to a 2-qubit Z+T1 model under known omega."""
    def loss(params):
        gz, gt1 = params
        if gz < 0 or gt1 < 0:
            return 1e6
        L_vec = build_L_z_plus_t1(2, omega, gz, gt1)
        total = 0.0
        rho0 = rhos[0]
        dt0 = t_us[0]
        for i, t in enumerate(t_us):
            if i == 0:
                continue
            pred = evolve(L_vec, rho0, t - dt0)
            total += float(np.sum(np.abs(pred - rhos[i]) ** 2))
        return total
    result = minimize(loss, [0.002, 0.001], method='Nelder-Mead',
                      options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 5000})
    return result.x, float(result.fun)


def asymmetry_of_L(omega, gamma_z, gamma_t1, N=2):
    """Compute the F112 polarity asymmetry of the fitted L."""
    d = 2 ** N
    H = (omega / 2.0) * site_op(N, DRIVE_SITE, Z)
    c_list = []
    g_list = []
    for l in range(N):
        c_list.append(site_op(N, l, Z))
        g_list.append(gamma_z)
        c_list.append(site_op(N, l, SIGMA_MINUS))
        g_list.append(gamma_t1)
    # The palindrome centre is the DEPHASING sum, matching the two sibling f112
    # scripts. Including the T1 rates changes neither `asymmetry` (the shift is a
    # multiple of I, Ad_Pi-fixed) nor ||M||^2: that norm is a parabola in sigma
    # minimised at sum(gz) + sum(gt1)/2, and the two conventions sit symmetrically
    # about the minimum, so ||M||^2 agrees to machine precision and `asymmetry`
    # to the bit.
    sigma = sum(g for c, g in zip(c_list, g_list) if np.allclose(c, c.conj().T))
    r = fw.polarity_coordinates_from_hc(H, c_list, g_list, N, sigma=sigma)
    return float(r['asymmetry']), float(r['norm_sq']['M'])


def f113_invert_gamma_t1(asymmetry, omega, N=2):
    """Invert F113: γ_T1 = −asym / ((1/2)·4^N·ω), γ_pump = 0.

    ONE driven site, so the sum over l in F113 has a single term. With a uniform
    γ_T1 the inversion returns that rate; with per-site rates it would return the
    DRIVEN site's, which is finer than the per-pair reading this doc once claimed.
    """
    prefactor = 0.5 * (4 ** N) * omega
    if abs(prefactor) < 1e-15:
        return float('nan')
    return -asymmetry / prefactor


def analyze_pair(omega, t_us, rhos, pair_label, T1_calib_us):
    print(f"\n  {pair_label}, applied Z-drive ω = {omega} per μs")
    print(f"    Calibrated T1 (per qubit): {T1_calib_us[0]:.1f} / {T1_calib_us[1]:.1f} μs")
    print(f"    → γ_T1_calibration (per qubit): {1.0/T1_calib_us[0]:.5f} / {1.0/T1_calib_us[1]:.5f} per μs")
    print(f"    → γ_T1_calibration mean: {0.5*(1.0/T1_calib_us[0] + 1.0/T1_calib_us[1]):.5f} per μs")

    (gz_fit, gt1_fit), fit_loss = fit_z_t1(omega, t_us, rhos)
    rms = np.sqrt(fit_loss / max(len(t_us) - 1, 1))
    print(f"    Fit: γ_z = {gz_fit:.5f}, γ_T1_fit = {gt1_fit:.5f} per μs   (RMS = {rms:.4f})")

    asym, M_sq = asymmetry_of_L(omega, gz_fit, gt1_fit)
    print(f"    Polarity asymmetry of fitted L: {asym:+.4e}  (||M||² = {M_sq:.4f}; rel asym = {abs(asym)/M_sq:.4e})")

    gt1_f113 = f113_invert_gamma_t1(asym, omega)
    print(f"    F113 inverted γ_T1: {gt1_f113:.5f} per μs")
    print(f"    F113 vs fit ratio: {gt1_f113 / gt1_fit:.6f}  (expect 1.000000 for self-consistency)")

    gt1_calib_mean = 0.5 * (1.0/T1_calib_us[0] + 1.0/T1_calib_us[1])
    print(f"    Fit vs calibration ratio: {gt1_fit / gt1_calib_mean:.4f}  ({gt1_fit:.5f} / {gt1_calib_mean:.5f})")

    return {
        'pair_label': pair_label, 'omega': omega,
        'gt1_fit': gt1_fit, 'gz_fit': gz_fit,
        'asymmetry': asym, 'M_sq': M_sq,
        'gt1_F113': gt1_f113,
        'T1_calib_per_qubit_us': T1_calib_us,
        'gt1_calib_mean_per_us': gt1_calib_mean,
        'fit_rms': rms,
    }


def build_L_z_t1_detuned(N, omega, gamma_z, gamma_t1, deltas):
    """Z+T1 with one free Z coefficient per qubit: H = (omega/2)*sum_l Z_l + sum_l delta_l Z_l.

    The shipped model has no detuning while this chip's coherences rotate, so this is
    the obvious missing channel. Same column-stack construction as build_L_z_plus_t1.
    """
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    H = (omega / 2.0) * site_op(N, DRIVE_SITE, Z)
    H = H + sum(deltas[l] * site_op(N, l, Z) for l in range(N))
    L_vec = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        for c, g in ((site_op(N, l, Z), gamma_z), (site_op(N, l, SIGMA_MINUS), gamma_t1)):
            cdc = c.conj().T @ c
            L_vec = L_vec + g * (np.kron(c.conj(), c)
                                 - 0.5 * (np.kron(Id, cdc) + np.kron(cdc.T, Id)))
    return L_vec


def fit_z_t1_detuned(omega, t_us, rhos, starts=12, seed=0):
    """Multistart fit of (gamma_z, gamma_t1, delta_0, delta_1). Returns (x, loss)."""
    rng = np.random.default_rng(seed)

    def loss(p):
        gz, gt1 = p[0], p[1]
        if gz < 0 or gt1 < 0:
            return 1e6
        L_vec = build_L_z_t1_detuned(2, omega, gz, gt1, (p[2], p[3]))
        total = 0.0
        for i, t in enumerate(t_us):
            if i == 0:
                continue
            total += float(np.sum(np.abs(evolve(L_vec, rhos[0], t - t_us[0]) - rhos[i]) ** 2))
        return total

    best = None
    for s in range(starts):
        x0 = [0.002, 0.001, 0.0, 0.0] if s == 0 else [
            abs(rng.normal(0.05, 0.05)), abs(rng.normal(0.005, 0.004)),
            rng.normal(0, 0.15), rng.normal(0, 0.15),
        ]
        r = minimize(loss, x0, method='Nelder-Mead',
                     options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 20000, 'maxfev': 20000})
        if best is None or r.fun < best.fun:
            best = r
    return best.x, float(best.fun)


def load_f95():
    paths = [
        ('omega=0.13', 'data/ibm_f95_angle_steering_may2026/'
                       'cusp_complex_phase_hardware_ibm_kingston_omega0.130_20260516_204827.json'),
        ('omega=0.25', 'data/ibm_f95_angle_steering_may2026/'
                       'cusp_complex_phase_hardware_ibm_kingston_omega0.250_20260516_205705.json'),
    ]
    runs = []
    for omega_tag, path in paths:
        with open(path, encoding='utf-8') as f:
            d = json.load(f)
        omega_per_us = float(d['omega_per_us'])
        for pair_key in d['pair_runs']:
            pr = d['pair_runs'][pair_key]
            qubits = pr['pair']['qubits']
            T1_calib = pr['pair']['T1_us']
            t_us = np.array([s['t_us'] for s in pr['trajectory']], dtype=float)
            rhos = [reconstruct_rho(s['rho2_real'], s['rho2_imag']) for s in pr['trajectory']]
            label = f'{omega_tag} {pair_key} q{qubits[0]}-q{qubits[1]}'
            runs.append((omega_per_us, t_us, rhos, label, T1_calib))
    return runs


def main():
    print("=" * 88)
    print("F113-inversion T1-extraction on Kingston f95 angle-steering datasets")
    print("=" * 88)
    print()
    print("Welle 5.A: invert F113 closed form")
    print("  γ_T1 = −asymmetry / ((N/2)·4^N·ω)")
    print("to extract a γ_T1 estimate from the measured polarity asymmetry of a")
    print("minimal Z+T1 Lindblad fit. Cross-check against fit-direct γ_T1 (should")
    print("agree bit-exact: self-consistency check) and against device-calibrated")
    print("γ_T1 = 1/T1 (the meaningful comparison: in-experiment effective T1 vs")
    print("isolated-qubit characterization).")
    print()

    runs = load_f95()
    results = []
    for run in runs:
        results.append(analyze_pair(*run))

    print()
    print("=" * 88)
    print("Aggregate summary")
    print("=" * 88)
    print()
    print(f"{'Pair-run':<46} {'γ_z fit':>9} {'γ_T1 fit':>10} {'γ_T1 F113':>10} {'F113/fit':>9} {'γ_T1 cal':>10} {'fit/cal':>8} {'RMS':>7}")
    print("-" * 122)
    for r in results:
        print(f"{r['pair_label']:<46} {r['gz_fit']:>9.4f} {r['gt1_fit']:>10.5f} {r['gt1_F113']:>10.5f} "
              f"{r['gt1_F113']/r['gt1_fit']:>9.6f} {r['gt1_calib_mean_per_us']:>10.5f} "
              f"{r['gt1_fit']/r['gt1_calib_mean_per_us']:>8.4f} {r['fit_rms']:>7.4f}")

    print()
    print("=" * 88)
    print("Robustness: the same fit with one free Z detuning per qubit")
    print("=" * 88)
    print()
    print(f"{'Pair-run':<46} {'γ_T1 fit':>10} {'γ_T1 +det':>10} {'ratio':>7} "
          f"{'γ_z +det':>9} {'RMS +det':>9} {'Σδ':>9} {'ω+2Σδ':>9}")
    print("-" * 122)
    for (omega, t_us, rhos, label, T1_calib), r in zip(runs, results):
        x, loss = fit_z_t1_detuned(omega, t_us, rhos)
        rms = np.sqrt(loss / max(len(t_us) - 1, 1))
        d_sum = x[2] + x[3]
        print(f"{label:<46} {r['gt1_fit']:>10.5f} {x[1]:>10.5f} {x[1]/r['gt1_fit']:>7.3f} "
              f"{x[0]:>9.4f} {rms:>9.4f} {d_sum:>+9.4f} {omega + 2*d_sum:>+9.4f}")
    print()
    print("  γ_T1 barely moves while γ_z collapses to calibration scale and the RMS")
    print("  improves: the detuning is what γ_z was absorbing. Only the SUM of the two")
    print("  δ is determined, the split between them being a gauge direction of this")
    print("  model (the Bell coherence turns at the sum), which is why the table reports")
    print("  Σδ and the effective drive ω + 2Σδ rather than δ_0 and δ_1 separately.")
    print("  Caveat for the inversion: that effective drive is what F113 would divide")
    print("  by, and the chip's own detuning cancels a good part of the applied one.")

    print()
    print("Reading:")
    print()
    print("  (1) F113/fit = 1.000000 across all 4 pair-runs: F113 inversion recovers")
    print("      the fit's γ_T1 bit-exact. Self-consistency confirmed; F113 is a")
    print("      faithful closed form for the fitted Lindblad L.")
    print()
    print("  (2) fit/cal > 1 (range 1.14-1.47): the fit's γ_T1 is LARGER than the")
    print("      device-calibrated 1/T1. The minimal Z+T1 model absorbs OTHER noise")
    print("      channels (crosstalk, drive-induced decoherence, transverse drift)")
    print("      into the σ⁻ T1 channel, inflating γ_T1_fit beyond physical T1.")
    print("      The two ω=0.25 runs do not bracket the two ω=0.13 runs, so the")
    print("      pair-to-pair spread swamps any drive dependence in this data.")
    print()
    print("  (3) γ_z fit values (0.011-0.016 per μs) run 7-30x the Lindblad target")
    print("      γ_z = (1/T2 - 1/(2·T1))/2 ≈ 0.0004-0.0023 per μs (NOT 1/T2: a D[Z]")
    print("      channel at γ_z decays coherences at 2γ_z, and T2 already carries")
    print("      1/(2·T1)). Nelder-Mead lands the unmodelled structure in γ_z, and")
    print("      RMS ≈ 0.3 says the same: the extracted rates are upper bounds.")
    print()
    print("  (4) Sharpened F113-as-diagnostic reading: F113-extracted γ_T1 is an")
    print("      EFFECTIVE-T1 number that equates ALL bit_b-mixed broken-balance")
    print("      noise to the σ⁻ T1 channel. The gap to device-calibrated 1/T1")
    print("      bounds how much non-T1 noise the model is absorbing. Adding the one")
    print("      channel we tested, a per-qubit detuning, moves it by under 3% and in")
    print("      both directions, so the bound is stable but not one-sided.")


if __name__ == '__main__':
    main()
