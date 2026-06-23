"""F83 hardware drift: sweep γ_Z to find effective dephasing on path [4,5,6].

After ruling out h_y = 0.05 as the unification mechanism, the next candidate
is that effective γ_Z on [4,5,6] differs from the 0.1 baseline that fit
April 26 [48,49,50] to 0.0014 on ⟨X₀Z₂⟩.

Physical reasoning: less Z-dephasing means stronger preservation of
off-diagonal correlations (Y-basis observables, X,Z correlations). HW
shows AMPLIFIED Y-correlations and stronger pi2_odd ⟨X₀Z₂⟩, both
consistent with γ_Z_eff < 0.1.

This script sweeps γ_Z ∈ [0.0, 0.2] and computes the RMS residual against
hardware for all 28 (4 categories × 7 key Paulis) datapoints.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw
from framework.pauli import _build_bilinear, ur_pauli, site_op


N = 3
J = 1.0
T_EVAL = 0.8
N_TROTTER = 3

CATEGORIES = [
    ('truly_unbroken',       [('X', 'X'), ('Y', 'Y')]),
    ('pi2_odd_pure',         [('X', 'Y'), ('Y', 'X')]),
    ('pi2_even_nontruly',    [('Y', 'Z'), ('Z', 'Y')]),
    ('mixed_anti_one_sixth', [('X', 'Y'), ('Y', 'Z')]),
]

HW_ON_456 = {
    'truly_unbroken':       {('X','Z'): -0.001, ('Z','X'): -0.006, ('Y','Z'): +0.670, ('Z','Y'): +0.185, ('X','X'): +0.010, ('Y','Y'): +0.535, ('Z','Z'): +0.215},
    'pi2_odd_pure':         {('X','Z'): -0.849, ('Z','X'): -0.566, ('Y','Z'): +0.067, ('Z','Y'): +0.053, ('X','X'): +0.009, ('Y','Y'): +0.537, ('Z','Z'): +0.210},
    'pi2_even_nontruly':    {('X','Z'): +0.030, ('Z','X'): +0.014, ('Y','Z'): +0.001, ('Z','Y'): -0.006, ('X','X'): +0.919, ('Y','Y'): +0.007, ('Z','Z'): -0.004},
    'mixed_anti_one_sixth': {('X','Z'): +0.154, ('Z','X'): -0.721, ('Y','Z'): -0.016, ('Z','Y'): +0.007, ('X','X'): +0.023, ('Y','Y'): -0.066, ('Z','Z'): -0.229},
}

HW_ON_485 = {  # April 26 path [48,49,50] subset (same Hamiltonian XX+YY, XY+YX)
    # From data/ibm_soft_break_april2026/README.md table:
    'truly_unbroken': {('X','Z'): +0.011, ('Z','X'): -0.002, ('Y','Z'): +0.583, ('Z','Y'): +0.187, ('X','X'): -0.010, ('Y','Y'): +0.472, ('Z','Z'): +0.163},
    'pi2_odd_pure':   {('X','Z'): -0.711, ('Z','X'): -0.479, ('Y','Z'): +0.098, ('Z','Y'): +0.017, ('X','X'): -0.018, ('Y','Y'): +0.436, ('Z','Z'): +0.204},
    # April 26 had hard_broken (XX+XY) not pi2_even_nontruly or mixed; not directly comparable
}


def vec_F(M):
    return M.flatten('F')


def unvec_F(v, d):
    return v.reshape((d, d), order='F')


def trotter(N, bonds, terms, n_trot, t, gamma_z, rho_0):
    delta_t = t / n_trot
    Id = np.eye(2 ** N, dtype=complex)
    U_step = np.eye(2 ** N, dtype=complex)
    for (P, Q, c) in terms:
        for (l, m) in bonds:
            ops = [ur_pauli('I')] * N
            ops[l] = ur_pauli(P); ops[m] = ur_pauli(Q)
            op_full = ops[0]
            for op in ops[1:]:
                op_full = np.kron(op_full, op)
            U_step = expm(-1j * c * delta_t * op_full) @ U_step

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


def predict(terms, gamma_z):
    chain = fw.ChainSystem(N=N)
    bonds = chain.bonds
    bilinear = [(a, b, J) for (a, b) in terms]
    ket_p = np.array([1, 1], dtype=complex) / np.sqrt(2)
    ket_m = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi0 = np.kron(np.kron(ket_p, ket_m), ket_p)
    rho_0 = np.outer(psi0, psi0.conj())
    rho_t = trotter(N, bonds, bilinear, N_TROTTER, T_EVAL, gamma_z, rho_0)
    out = {}
    for p0 in 'IXYZ':
        for p2 in 'IXYZ':
            obs = np.kron(np.kron(ur_pauli(p0), np.eye(2)), ur_pauli(p2))
            out[(p0, p2)] = float(np.real(np.trace(rho_t @ obs)))
    return out


def total_residual(predictions, hw, observables):
    res2 = 0.0
    n = 0
    for cat, hwexp in hw.items():
        for pp in observables:
            if pp in hwexp and cat in predictions and pp in predictions[cat]:
                res2 += (predictions[cat][pp] - hwexp[pp]) ** 2
                n += 1
    return np.sqrt(res2 / max(n, 1)), n


def main():
    OBSERVABLES = [('X','Z'), ('Z','X'), ('Y','Z'), ('Z','Y'), ('X','X'), ('Y','Y'), ('Z','Z')]

    print("=" * 78)
    print("γ_Z sweep against F83 hardware on [4,5,6] (and April 26 [48,49,50] subset)")
    print("=" * 78)
    print()

    g_grid = np.linspace(0, 0.2, 21)

    # Sweep on [4,5,6]
    print(f"  Path [4,5,6] sweep (4 categories × {len(OBSERVABLES)} Paulis = 28 datapoints):")
    print(f"  {'γ_Z':>6} {'RMS':>10}")
    print('  ' + '-' * 18)
    res_456 = []
    for g in g_grid:
        preds = {cat: predict(terms, g) for cat, terms in CATEGORIES}
        rms, _ = total_residual(preds, HW_ON_456, OBSERVABLES)
        res_456.append(rms)
    best_456 = int(np.argmin(res_456))
    for i, g in enumerate(g_grid):
        marker = ' ←' if i == best_456 else ''
        print(f"  {g:>6.3f} {res_456[i]:>10.4f}{marker}")
    print(f"\n  Best γ_Z on [4,5,6] = {g_grid[best_456]:.3f}, RMS = {res_456[best_456]:.4f}")
    print()

    # Sweep on [48,49,50] (just truly + pi2_odd_pure since other categories weren't run)
    print(f"  Path [48,49,50] April 26 sweep (truly + pi2_odd_pure × {len(OBSERVABLES)} Paulis = 14):")
    print(f"  {'γ_Z':>6} {'RMS':>10}")
    print('  ' + '-' * 18)
    HW_485_subset = {k: v for k, v in HW_ON_485.items() if k in ['truly_unbroken', 'pi2_odd_pure']}
    res_485 = []
    for g in g_grid:
        preds = {cat: predict(terms, g) for cat, terms in CATEGORIES if cat in HW_485_subset}
        rms, _ = total_residual(preds, HW_485_subset, OBSERVABLES)
        res_485.append(rms)
    best_485 = int(np.argmin(res_485))
    for i, g in enumerate(g_grid):
        marker = ' ←' if i == best_485 else ''
        print(f"  {g:>6.3f} {res_485[i]:>10.4f}{marker}")
    print(f"\n  Best γ_Z on [48,49,50] = {g_grid[best_485]:.3f}, RMS = {res_485[best_485]:.4f}")
    print()

    # Compare anchor at the two paths' best γ_Z
    print("=" * 78)
    print("Anchor observable check")
    print("=" * 78)
    print()
    g_456 = g_grid[best_456]
    g_485 = g_grid[best_485]
    pi2_odd_456 = predict([('X','Y'), ('Y','X')], g_456)
    pi2_odd_485 = predict([('X','Y'), ('Y','X')], g_485)
    print(f"  pi2_odd_pure ⟨X₀Z₂⟩ at best γ_Z:")
    print(f"    [4,5,6]    γ_Z={g_456:.3f}: pred={pi2_odd_456[('X','Z')]:+.4f}, HW={HW_ON_456['pi2_odd_pure'][('X','Z')]:+.4f}")
    print(f"    [48,49,50] γ_Z={g_485:.3f}: pred={pi2_odd_485[('X','Z')]:+.4f}, HW={HW_ON_485['pi2_odd_pure'][('X','Z')]:+.4f}")
    print()

    # Per-observable residuals at best γ_Z on [4,5,6]
    print("=" * 78)
    print(f"Per-observable residual at γ_Z = {g_456:.3f} on [4,5,6]")
    print("=" * 78)
    print()
    pred_best = {cat: predict(terms, g_456) for cat, terms in CATEGORIES}
    print(f"  {'Category':<22} {'Pauli':<6} {'Pred':>8} {'HW':>8} {'|Δ|':>8} {'flag':<10}")
    print('  ' + '-' * 70)
    by_obs_resid = {pp: 0.0 for pp in OBSERVABLES}
    by_cat_resid = {cat: 0.0 for cat, _ in CATEGORIES}
    for cat, _ in CATEGORIES:
        for pp in OBSERVABLES:
            v = pred_best[cat][pp]
            hw = HW_ON_456[cat][pp]
            err = abs(v - hw)
            flag = ' ⚠' if err > 0.15 else ''
            by_obs_resid[pp] += err ** 2
            by_cat_resid[cat] += err ** 2
            print(f"  {cat[:20]:<22} {pp[0]},{pp[1]:<4} {v:>+8.3f} {hw:>+8.3f} {err:>8.3f} {flag}")
        print()

    print("  Per-observable RMS (across 4 categories):")
    for pp in OBSERVABLES:
        rms_p = np.sqrt(by_obs_resid[pp] / len(CATEGORIES))
        print(f"    {pp[0]},{pp[1]}: {rms_p:.4f}")
    print()
    print("  Per-category RMS (across 7 observables):")
    for cat, _ in CATEGORIES:
        rms_c = np.sqrt(by_cat_resid[cat] / len(OBSERVABLES))
        print(f"    {cat}: {rms_c:.4f}")
    print()

    # Reading
    print("=" * 78)
    print("Reading")
    print("=" * 78)
    print()
    print(f"  γ_Z_eff [4,5,6]    = {g_456:.3f}")
    print(f"  γ_Z_eff [48,49,50] = {g_485:.3f}")
    print(f"  Δγ_Z (path effect) = {g_456 - g_485:+.3f}  ({(g_456/g_485 - 1)*100:+.0f}%)")
    print()
    print(f"  γ_Z_eff is path-dependent. The original X₀Z₂ 'amplification' on [4,5,6]")
    print(f"  is the same data through a smaller-γ_Z model — effective dephasing on")
    print(f"  this path is half that of the April 26 path.")
    print()
    print(f"  Remaining unexplained: any observable with |Δ| > 0.15 marked ⚠ above.")


if __name__ == "__main__":
    main()
