"""F83 hardware drift: does Marrakesh's known h_y ≈ 0.05 unify the X₀Z₂ and Y-basis surprises?

Background:
  zn_mirror Hardware-Diagnostic (2026-04-29, Job d7ornigror3c73c0c6ug) found
  Marrakesh has a transverse Y-field at the Hamiltonian level:
    h_y_eff ≈ 0.05  (single-site Σ Y_l per qubit)
  See project memory project_zn_mirror_diagnostic and Confirmations
  entry 'marrakesh_transverse_y_field_detection'.

  The April 30 F83 4-Hamiltonian run on path [4,5,6] showed two surprises:
    1. pi2_odd_pure ⟨X₀Z₂⟩ = -0.849, expected -0.723 (Trotter+γ_Z=0.1)
    2. Y-basis amplification: truly ⟨Y₀Z₂⟩ = +0.670, expected 0.381
  Y-amplification was also visible on April 26 [48,49,50]: ⟨Y₀Z₂⟩ = 0.583
  vs same 0.381 prediction. So Y-bias is backend-wide, not path-specific.

Hypothesis: the 0.05 h_y field, baked into the Hamiltonian during the
~5 μs circuit, is what shifts both the Y-basis observables AND the X,Z
amplitude. If yes, we have a coherent cross-confirmation linking the F83
result to the zn_mirror finding.

This script:
  - Adds H_drift = h_y · Σ_l Y_l to each Trotter step (small-step approx).
  - Re-computes Trotter+γ_Z=0.1 predictions for all 4 F83 categories.
  - Compares to hardware on [4,5,6].
  - Sweeps h_y ∈ [0, 0.1] to find the best-fit value vs hardware.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw
from framework.lindblad import lindbladian_z_dephasing
from framework.pauli import _build_bilinear, ur_pauli, site_op


N = 3
J = 1.0
T_EVAL = 0.8
N_TROTTER = 3
GAMMA_Z = 0.1

CATEGORIES = [
    ('truly_unbroken',       [('X', 'X'), ('Y', 'Y')]),
    ('pi2_odd_pure',         [('X', 'Y'), ('Y', 'X')]),
    ('pi2_even_nontruly',    [('Y', 'Z'), ('Z', 'Y')]),
    ('mixed_anti_one_sixth', [('X', 'Y'), ('Y', 'Z')]),
]

# Hardware values from data/ibm_f83_signature_april2026/...
HW_ON_456 = {
    'truly_unbroken':       {('X','Z'): -0.001, ('Z','X'): -0.006, ('Y','Z'): +0.670, ('Z','Y'): +0.185, ('X','X'): +0.010, ('Y','Y'): +0.535, ('Z','Z'): +0.215},
    'pi2_odd_pure':         {('X','Z'): -0.849, ('Z','X'): -0.566, ('Y','Z'): +0.067, ('Z','Y'): +0.053, ('X','X'): +0.009, ('Y','Y'): +0.537, ('Z','Z'): +0.210},
    'pi2_even_nontruly':    {('X','Z'): +0.030, ('Z','X'): +0.014, ('Y','Z'): +0.001, ('Z','Y'): -0.006, ('X','X'): +0.919, ('Y','Y'): +0.007, ('Z','Z'): -0.004},
    'mixed_anti_one_sixth': {('X','Z'): +0.154, ('Z','X'): -0.721, ('Y','Z'): -0.016, ('Z','Y'): +0.007, ('X','X'): +0.023, ('Y','Y'): -0.066, ('Z','Z'): -0.229},
}


def vec_F(M):
    return M.flatten('F')


def unvec_F(v, d):
    return v.reshape((d, d), order='F')


def build_h_drift_y(N, h_y):
    """H_drift = h_y · Σ_l Y_l (single-site transverse Y-field on each qubit)."""
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for l in range(N):
        H += h_y * site_op(N, l, 'Y')
    return H


def trotter_with_drift(N, bonds, terms, h_y, n_trot, t, gamma_z, rho_0):
    """Trotter circuit with bilinear J·H_bonds + drift h_y·Σ Y_l + per-step Z-dephasing.

    The drift acts as a small extra term inside each Trotter step:
        U_step = exp(-i·H_drift·δt) · ∏_bonds exp(-i·c·H_bond·δt)
    """
    delta_t = t / n_trot
    Id = np.eye(2 ** N, dtype=complex)
    H_drift = build_h_drift_y(N, h_y)
    U_drift_step = expm(-1j * H_drift * delta_t)

    U_step = np.eye(2 ** N, dtype=complex)
    for (P, Q, c) in terms:
        for (l, m) in bonds:
            ops = [ur_pauli('I')] * N
            ops[l] = ur_pauli(P); ops[m] = ur_pauli(Q)
            op_full = ops[0]
            for op in ops[1:]:
                op_full = np.kron(op_full, op)
            U_step = expm(-1j * c * delta_t * op_full) @ U_step
    U_step = U_drift_step @ U_step

    L_deph = np.zeros((4 ** N, 4 ** N), dtype=complex)
    for l in range(N):
        Zl = site_op(N, l, 'Z')
        L_deph += gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    M_deph = expm(L_deph * delta_t)

    rho = rho_0
    for _ in range(n_trot):
        rho = U_step @ rho @ U_step.conj().T
        rho = unvec_F(M_deph @ vec_F(rho), 2 ** N)
    return rho


def predict_category(terms, h_y):
    chain = fw.ChainSystem(N=N)
    bonds = chain.bonds
    bilinear = [(a, b, J) for (a, b) in terms]
    ket_p = np.array([1, 1], dtype=complex) / np.sqrt(2)
    ket_m = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi0 = np.kron(np.kron(ket_p, ket_m), ket_p)
    rho_0 = np.outer(psi0, psi0.conj())
    rho_t = trotter_with_drift(N, bonds, bilinear, h_y, N_TROTTER, T_EVAL, GAMMA_Z, rho_0)
    out = {}
    for p0 in 'IXYZ':
        for p2 in 'IXYZ':
            obs = np.kron(np.kron(ur_pauli(p0), np.eye(2)), ur_pauli(p2))
            out[(p0, p2)] = float(np.real(np.trace(rho_t @ obs)))
    return out


def total_residual(predictions, hw, observables):
    res2 = 0.0
    for cat, hwexp in hw.items():
        for pp in observables:
            if pp in hwexp and pp in predictions[cat]:
                res2 += (predictions[cat][pp] - hwexp[pp]) ** 2
    return np.sqrt(res2)


def main():
    print("=" * 78)
    print("F83 drift unification test: does h_y ≈ 0.05 (zn_mirror) explain the surprises?")
    print("=" * 78)
    print()

    OBSERVABLES = [('X','Z'), ('Z','X'), ('Y','Z'), ('Z','Y'), ('X','X'), ('Y','Y'), ('Z','Z')]

    # Reference: h_y = 0
    print("Step 1. Trotter+γ_Z=0.1 with h_y = 0 (current model):")
    print()
    pred_h0 = {cat: predict_category(terms, h_y=0.0) for cat, terms in CATEGORIES}
    res0 = total_residual(pred_h0, HW_ON_456, OBSERVABLES)
    print(f"  RMS residual vs hardware (28 datapoints): {res0:.4f}")
    print()

    # With zn_mirror's measured h_y_eff = 0.05
    print("Step 2. Trotter+γ_Z=0.1 with h_y = 0.05 (zn_mirror April 29 finding):")
    print()
    pred_h05 = {cat: predict_category(terms, h_y=0.05) for cat, terms in CATEGORIES}
    res05 = total_residual(pred_h05, HW_ON_456, OBSERVABLES)
    print(f"  RMS residual vs hardware: {res05:.4f}")
    print(f"  Improvement: {(res0 - res05)/res0*100:+.1f}%")
    print()

    # Per-observable diagnostic at h_y=0.05
    print("Step 3. Per-observable comparison at h_y=0.05:")
    print()
    print(f"  {'Cat':<22} {'Pauli':<6} {'h_y=0':>10} {'h_y=0.05':>10} {'HW':>10} {'h0 err':>8} {'h05 err':>8}")
    print('  ' + '-' * 84)
    for cat, _ in CATEGORIES:
        for pp in OBSERVABLES:
            v0 = pred_h0[cat][pp]
            v05 = pred_h05[cat][pp]
            hw = HW_ON_456[cat].get(pp, np.nan)
            err0 = abs(v0 - hw) if not np.isnan(hw) else np.nan
            err05 = abs(v05 - hw) if not np.isnan(hw) else np.nan
            print(f"  {cat[:20]:<22} {pp[0]},{pp[1]:<4} {v0:>+10.4f} {v05:>+10.4f} {hw:>+10.4f} {err0:>8.4f} {err05:>8.4f}")
        print()

    # Sweep h_y
    print("Step 4. Sweep h_y ∈ [0, 0.10] to find best fit:")
    print()
    h_grid = np.linspace(0, 0.10, 21)
    residuals = []
    for h in h_grid:
        preds = {cat: predict_category(terms, h_y=h) for cat, terms in CATEGORIES}
        residuals.append(total_residual(preds, HW_ON_456, OBSERVABLES))
    best_idx = int(np.argmin(residuals))
    print(f"  {'h_y':>6} {'RMS residual':>14}")
    print('  ' + '-' * 22)
    for i, h in enumerate(h_grid):
        marker = ' ←' if i == best_idx else ''
        print(f"  {h:>6.3f} {residuals[i]:>14.4f}{marker}")
    print()
    print(f"  Best h_y = {h_grid[best_idx]:.3f}, RMS = {residuals[best_idx]:.4f}")
    print(f"  zn_mirror April 29 measured: h_y_eff ≈ 0.05")
    print()

    # Anchor observables
    print("Step 5. Headline: how do the two surprises move?")
    print()
    print(f"  Surprise 1 (pi2_odd_pure ⟨X₀Z₂⟩):")
    print(f"    h_y=0:    {pred_h0['pi2_odd_pure'][('X','Z')]:>+.4f}")
    print(f"    h_y=0.05: {pred_h05['pi2_odd_pure'][('X','Z')]:>+.4f}")
    print(f"    HW:       {HW_ON_456['pi2_odd_pure'][('X','Z')]:>+.4f}")
    print()
    print(f"  Surprise 2 (truly ⟨Y₀Z₂⟩):")
    print(f"    h_y=0:    {pred_h0['truly_unbroken'][('Y','Z')]:>+.4f}")
    print(f"    h_y=0.05: {pred_h05['truly_unbroken'][('Y','Z')]:>+.4f}")
    print(f"    HW:       {HW_ON_456['truly_unbroken'][('Y','Z')]:>+.4f}")
    print()
    print(f"  Surprise 3 (truly ⟨Y₀Y₂⟩):")
    print(f"    h_y=0:    {pred_h0['truly_unbroken'][('Y','Y')]:>+.4f}")
    print(f"    h_y=0.05: {pred_h05['truly_unbroken'][('Y','Y')]:>+.4f}")
    print(f"    HW:       {HW_ON_456['truly_unbroken'][('Y','Y')]:>+.4f}")
    print()


if __name__ == "__main__":
    main()
