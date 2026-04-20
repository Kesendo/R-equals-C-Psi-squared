#!/usr/bin/env python3
"""c1_bilinearity_test.py

Approach D from the "how to solve the magnitudes puzzle" analysis.

The prior investigation showed that c_1 is bilinear in rho_0 (per-site
purity Tr(rho_i^2) is bilinear in rho, and c_1 = Sum_i f_i comes from
that bilinear). So:

    c_1(rho_0) = Sum_{mu, nu} K^{mu, nu} (rho_0)_mu (rho_0)_nu

This script tests the bilinearity by using Dicke states as a probe basis
and measuring c_1 for several pure-sector and cross-sector initial states.

Dicke state S_n at N=5:
    |S_n> = (1/sqrt(C(N,n))) * Sum_{|x|=n} |x>
where the sum is over all basis states with n excitations. At N=5 the
sector sizes are C(5,n) = {1, 5, 10, 10, 5, 1} for n=0..5.

Initial states tested:
  A. (|vac> + |psi_1>)/sqrt(2)              PTF standard bonding, known c_1
  B. (|vac> + |S_1>)/sqrt(2)                W_5 + vac superposition
  C. (|vac> + |S_2>)/sqrt(2)                Dicke-2 + vac
  D. (|vac> + |S_3>)/sqrt(2)                Dicke-3 + vac
  E. (|S_1> + |S_2>)/sqrt(2)                single cross-sector coherence
  F. (|S_2> + |S_3>)/sqrt(2)                another cross-sector
  G. |+>^5                                  full uniform superposition
  H. |S_1> (pure Dicke state, no superpos)  sector-diagonal, (1,1) only
  I. |S_2>                                  (2,2) only
  J. (|vac> + |S_1> + |S_2>)/sqrt(3)        three-component superposition

Bilinearity test: if c_1 is bilinear, then the cross-sector state
coherences should give predictable combinations of the diagonal c_1
values. Concrete test: does c_1(E) relate to c_1(H) + c_1(I) + (cross-term)?

N=5 with d^2 = 1024, three eigendecompositions take ~15 seconds total.
Propagation per state is fast. Whole run ~1-2 minutes.
"""
from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS, T_FIT_MAX,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, density_matrix,
    per_site_purity, fit_alpha,
)

RESULTS_DIR = Path(__file__).parent / "results" / "c1_bilinearity_test"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

N = 5
DJ_EXTRACT = 0.01


def dicke_state(N, n):
    """Return the |S_n> Dicke state on N qubits with n excitations.
    Big-endian: qubit 0 is the MSB."""
    d = 2**N
    v = np.zeros(d, dtype=complex)
    count = comb(N, n)
    for bits in range(d):
        if bin(bits).count('1') == n:
            v[bits] = 1.0
    v /= np.sqrt(count)
    return v


def plus_state(N):
    """|+>^N = uniform superposition of all 2^N basis states."""
    d = 2**N
    return np.ones(d, dtype=complex) / np.sqrt(d)


def bonding_plus_vacuum(N, k=1):
    v = vacuum_ket(N)
    psi = single_excitation_mode(N, k=k)
    return (v + psi) / np.sqrt(2.0)


def eig_and_inv(L):
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    return eigvals, V_R, V_Linv


def propagate(eigvals, V_R, V_Linv, rho_0, times):
    d = rho_0.shape[0]
    rho0_vec = rho_0.flatten(order='F')
    c0 = V_Linv @ rho0_vec
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order='F')
    return out


def measure_c1(label, ket, evA, VRA, VLA, evBp, VRBp, VLBp, evBm, VRBm, VLBm, N, times):
    rho_0 = density_matrix(ket)
    rho_A = propagate(evA, VRA, VLA, rho_0, times)
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    P_A = per_site_purity(rho_A, N)
    P_Bp = per_site_purity(rho_Bp, N)
    P_Bm = per_site_purity(rho_Bm, N)
    alpha_p = np.zeros(N); alpha_m = np.zeros(N)
    rmse_p = np.zeros(N); rmse_m = np.zeros(N)
    for i in range(N):
        a, r = fit_alpha(times, P_A[:, i], P_Bp[:, i])
        alpha_p[i] = a; rmse_p[i] = r
        a, r = fit_alpha(times, P_A[:, i], P_Bm[:, i])
        alpha_m[i] = a; rmse_m[i] = r
    closure_p = float(np.sum(np.log(alpha_p)))
    closure_m = float(np.sum(np.log(alpha_m)))
    c_1 = (closure_p - closure_m) / (2 * DJ_EXTRACT)
    return {
        "label": label,
        "c_1": c_1,
        "closure_plus": closure_p,
        "closure_minus": closure_m,
        "alpha_plus": alpha_p.tolist(),
        "alpha_minus": alpha_m.tolist(),
        "rmse_plus_max": float(rmse_p.max()),
        "rmse_minus_max": float(rmse_m.max()),
    }


def main():
    print("=" * 70)
    print(f"c_1 bilinearity test at N = {N}")
    print("=" * 70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}")
    print(f"  d^2 = {4**N} = {4**N}")
    print()

    # Build Liouvillians
    J_A = [J_UNIFORM] * (N - 1)
    J_B_plus = list(J_A); J_B_plus[0] = J_UNIFORM + DJ_EXTRACT
    J_B_minus = list(J_A); J_B_minus[0] = J_UNIFORM - DJ_EXTRACT

    print("  Building Liouvillians and eigendecomposing...")
    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    evA, VRA, VLA = eig_and_inv(L_A)
    del L_A
    L_Bp = build_liouvillian_matrix(build_H_XY(J_B_plus, N), GAMMA_0, N)
    evBp, VRBp, VLBp = eig_and_inv(L_Bp)
    del L_Bp
    L_Bm = build_liouvillian_matrix(build_H_XY(J_B_minus, N), GAMMA_0, N)
    evBm, VRBm, VLBm = eig_and_inv(L_Bm)
    del L_Bm
    print(f"  Built in {time.time()-t0:.1f} s")

    times_arr = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    def m(label, ket):
        return measure_c1(label, ket, evA, VRA, VLA,
                          evBp, VRBp, VLBp, evBm, VRBm, VLBm,
                          N, times_arr)

    # ========== Define test states ==========
    S = {n: dicke_state(N, n) for n in range(N + 1)}
    vac = vacuum_ket(N)
    psi_1 = single_excitation_mode(N, k=1)

    def norm(v):
        return v / np.linalg.norm(v)

    tests = [
        ("A. (|vac>+|psi_1>)/sqrt(2) [PTF standard]", norm(vac + psi_1)),
        ("B. (|vac>+|S_1>)/sqrt(2) [W+vac]",          norm(vac + S[1])),
        ("C. (|vac>+|S_2>)/sqrt(2)",                  norm(vac + S[2])),
        ("D. (|vac>+|S_3>)/sqrt(2)",                  norm(vac + S[3])),
        ("E. (|S_1>+|S_2>)/sqrt(2)",                  norm(S[1] + S[2])),
        ("F. (|S_2>+|S_3>)/sqrt(2)",                  norm(S[2] + S[3])),
        ("G. |+>^N",                                  plus_state(N)),
        ("H. |S_1> pure (W_5)",                       S[1]),
        ("I. |S_2> pure",                             S[2]),
        ("J. |S_3> pure",                             S[3]),
        ("K. (|vac>+|S_1>+|S_2>)/sqrt(3)",            norm(vac + S[1] + S[2])),
        ("L. (|vac>+|S_4>)/sqrt(2) [Pi partner of C?]", norm(vac + S[4])),
        ("M. |S_4> pure",                             S[4]),
        ("N. |S_5> pure [all excited]",               S[5]),
        ("O. |vac> pure [ground]",                    vac),
    ]

    print(f"\n  Per-state c_1 measurements:")
    print(f"  {'label':>50} {'c_1':>10} {'RMSE+':>10}")
    results = {}
    for label, ket in tests:
        r = m(label, ket)
        print(f"  {label:>50} {r['c_1']:>+10.4f} {r['rmse_plus_max']:>10.1e}")
        results[label] = r

    # ========== Bilinearity analysis ==========
    def get(label_prefix):
        for k, v in results.items():
            if k.startswith(label_prefix):
                return v["c_1"]
        return None

    c1_A = get("A."); c1_B = get("B."); c1_C = get("C."); c1_D = get("D.")
    c1_E = get("E."); c1_F = get("F."); c1_G = get("G."); c1_H = get("H.")
    c1_I = get("I."); c1_J = get("J."); c1_K = get("K."); c1_L = get("L.")
    c1_M = get("M."); c1_N = get("N."); c1_O = get("O.")

    print(f"\n  Bilinearity-structure tests:")

    # Test 1: pure Dicke states (|S_n><S_n|), diagonal sectors only.
    # If c_1 lives mainly in (n, n) diagonal, these give diagonal K_nn contributions.
    print(f"    Pure Dicke state c_1 values (diagonal sectors only):")
    for n, label in [(0, "O"), (1, "H"), (2, "I"), (3, "J"), (4, "M"), (5, "N")]:
        c = get(label + ".")
        print(f"      |S_{n}> (sector {n}, dim {comb(N, n)}): c_1 = {c:+.5f}")

    # Test 2: bilinear decomposition prediction
    # For rho = (|a> + |b>)/sqrt(2), rho = (1/2)[|a><a| + |b><b| + |a><b| + |b><a|]
    # c_1(rho) = (1/4)[c_1(|a><a|) + c_1(|b><b|) + (cross-term)]
    # Specifically: if c_1 is bilinear, c_1(|a>+|b|)^2 contains:
    #   |a|^4 term: c_1(|a><a| normalized unit) = c_1(pure |a>)
    #   |b|^4 term: c_1(pure |b>)
    #   2|a|^2|b|^2 + cross: from off-diagonal coherence
    # At equal amplitudes a=b=1/sqrt(2), the full state is 1/sqrt(2)(|a>+|b>) with
    # rho = (1/2)(|a><a|+|a><b|+|b><a|+|b><b|)
    # Each diagonal contributes c_1(pure)/4 via bilinear (rho has norm |amp|^2 = 1/2
    # for each diagonal), so the prediction is:
    #   c_1(mix) = (1/4)*c_1(|a><a|) + (1/4)*c_1(|b><b|) + (cross-terms from off-diagonals)
    #
    # Measurable: compute predicted diagonal sum and the residual is the cross-term.
    print(f"\n    Sum vs superposition test (if bilinear, diagonal contributions scale):")
    # (|vac>+|S_1>)/sqrt(2) = B. Pure components: |vac> (=O) and |S_1> (=H).
    # rho_B = (1/2)[|vac><vac| + |S_1><S_1| + |vac><S_1| + |S_1><vac|]
    # For bilinear c_1, the (|vac><vac|, |vac><vac|) self-bilinear gives c_1(O)/4 contribution
    # scaled... actually the scaling is subtle because c_1 is over Tr(rho_i^2) which is bilinear
    # in rho but with normalization. Pure |vac> has rho_vac = |vac><vac|, Tr(rho_vac^2) = 1.
    # Superposition state B has rho_B with Tr(rho_B^2) = 1 (still pure). So the normalization
    # is the same. The bilinear prediction: c_1(B) = (1/4)[c_1(O) + c_1(H)] + cross.
    # But c_1(O) might be 0 (|vac> is stationary).
    pred_B_diag = 0.25 * (c1_O + c1_H)
    cross_B = c1_B - pred_B_diag
    print(f"      B = (|vac>+|S_1>)/sqrt(2): c_1 = {c1_B:+.5f}")
    print(f"        pred from diag: (c_O + c_H)/4 = {pred_B_diag:+.5f}")
    print(f"        cross-term residual: {cross_B:+.5f}")

    pred_C_diag = 0.25 * (c1_O + c1_I)
    cross_C = c1_C - pred_C_diag
    print(f"      C = (|vac>+|S_2>)/sqrt(2): c_1 = {c1_C:+.5f}")
    print(f"        pred from diag: (c_O + c_I)/4 = {pred_C_diag:+.5f}")
    print(f"        cross-term residual: {cross_C:+.5f}")

    pred_E_diag = 0.25 * (c1_H + c1_I)
    cross_E = c1_E - pred_E_diag
    print(f"      E = (|S_1>+|S_2>)/sqrt(2): c_1 = {c1_E:+.5f}")
    print(f"        pred from diag: (c_H + c_I)/4 = {pred_E_diag:+.5f}")
    print(f"        cross-term residual: {cross_E:+.5f}")

    # Test 3: Pi-pair check across sectors. S_n <-> S_{N-n} under excitation inversion
    # (simultaneous bit flip) which is an important symmetry.
    print(f"\n    Sector-inversion symmetry (S_n <-> S_{{N-n}}):")
    print(f"      c_1(|S_1>) = {c1_H:+.5f},  c_1(|S_4>) = {c1_M:+.5f},  diff = {c1_H - c1_M:+.1e}")
    print(f"      c_1(|S_2>) = {c1_I:+.5f},  c_1(|S_3>) = {c1_J:+.5f},  diff = {c1_I - c1_J:+.1e}")

    # Test 4: |+>^N decomposition. |+>^N = Sum_n sqrt(C(N,n)/2^N) |S_n>.
    # Coefficients at N=5: sqrt(1/32), sqrt(5/32), sqrt(10/32), sqrt(10/32), sqrt(5/32), sqrt(1/32).
    # rho_|+> = Sum_{n,m} c_n c_m* |S_n><S_m|
    # Diagonal part: Sum_n |c_n|^2 |S_n><S_n| = Sum_n (C(N,n)/2^N) |S_n><S_n|
    # c_1(diagonal-only of |+>) = Sum_n (C(N,n)/2^N)^2 * c_1(|S_n>) [if c_1 pure bilinear diagonal]
    # Actually c_1 is bilinear in rho, so:
    # c_1(diagonal rho) = Sum_n (coef_n)^2 c_1(|S_n>) where coef_n = C(N,n)/2^N is weight
    # Hmm wait this doesn't quite work. Let me think. If rho = Sum_n w_n |S_n><S_n|, diagonal,
    # then c_1(rho) involves Tr(rho_i^2) = Sum_{n,m} w_n w_m Tr(Tr_{not_i}(|S_n><S_n|) Tr_{not_i}(|S_m><S_m|)).
    # That's bilinear in w_n. Not just Sum w_n^2.
    weights = [comb(N, n) / 2**N for n in range(N + 1)]
    # Bilinear prediction for |+>^N c_1 from diagonal Dicke-only terms
    # would require the per-pair (n, m) cross-term Tr(Tr_i(|S_n><S_n|) Tr_i(|S_m><S_m|))
    # which we do not directly measure from pure |S_n> states alone. But an upper-bound
    # "if all cross-terms were same as diagonal" is:
    c1_pure = [c1_O, c1_H, c1_I, c1_J, c1_M, c1_N]
    pred_plus_diag_only = sum(w * w * c for w, c in zip(weights, c1_pure))
    print(f"\n    |+>^N = Sum_n sqrt(C(N,n)/2^N) |S_n>. Weights w_n = C(N,n)/2^N:")
    print(f"      weights (n=0..N): {[f'{w:.4f}' for w in weights]}")
    print(f"      c_1 pure |S_n>:   {[f'{c:+.4f}' for c in c1_pure]}")
    print(f"      Sum w_n^2 c_1(|S_n>) (diagonal-only prediction): {pred_plus_diag_only:+.5f}")
    print(f"      c_1(|+>^N) measured:                               {c1_G:+.5f}")
    print(f"      cross-sector contribution (|+>^N − diag):         {c1_G - pred_plus_diag_only:+.5f}")

    # Save
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "defect_bond": [0, 1],
        "results": results,
        "bilinearity_tests": {
            "B_decomp": {"c_1": c1_B, "pred_diag": pred_B_diag, "cross": cross_B},
            "C_decomp": {"c_1": c1_C, "pred_diag": pred_C_diag, "cross": cross_C},
            "E_decomp": {"c_1": c1_E, "pred_diag": pred_E_diag, "cross": cross_E},
            "plus_diag_prediction": {"c_1_plus": c1_G, "diag_pred": pred_plus_diag_only,
                                      "cross_contribution": c1_G - pred_plus_diag_only,
                                      "weights": weights,
                                      "c_1_pure_Sn": c1_pure},
            "sector_inversion": {
                "S1_vs_S4_diff": c1_H - c1_M,
                "S2_vs_S3_diff": c1_I - c1_J,
            },
        },
    }
    path = RESULTS_DIR / "bilinearity_test.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
