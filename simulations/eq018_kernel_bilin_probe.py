#!/usr/bin/env python3
"""eq018_kernel_bilin_probe.py

Tests whether c_1 is strictly bilinear in rho_0 or only approximately so.

The claim in EQ-018 is:
    c_1[rho_0] = sum K_{A,B} rho_0^A rho_0^B      (strict bilinearity)

If strict, then K^CC[0, 1] extracted two ways must agree:
  way A (subtraction): K^CC = 2 * (c_1[coherent] - c_1[mixed])
      using rho_coh = (|vac> + |S_1>)(<vac|+<S_1|)/2 and
            rho_mix = (|vac><vac| + |S_1><S_1|)/2
  way B (direct probe, only works if strictly bilinear because
         the "probe" rho_coh_only is NOT a valid density matrix):
      Propagate rho_coh_only = (|vac><S_1| + |S_1><vac|)/2 directly
      through L_A and L_B, fit alpha_i on the resulting Tr(rho_i^2)
      trajectories, compute c_1. Then K^CC = 2 * c_1^{pure coherence}.

A second test: scale the coherent state weights away from 1/2.
  rho(w) = |psi(w)><psi(w)| where
      |psi(w)> = (cos(w) |vac> + sin(w) |S_1>)
  The coherence-block contribution scales as sin(w) cos(w), and the
  diagonal contributions scale as cos^2(w) and sin^2(w). Under strict
  bilinearity c_1 should be a specific quadratic function of w. Fit
  a quadratic and see if it matches prediction.

Runs on N=5, bond=0 (canonical scientific target). Expected runtime
< 30 s.

Rules from task: em-dashes forbidden. Hyphens only. UTF-8 stdout.
"""
from __future__ import annotations

import json
import sys
import time
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
    vacuum_ket, density_matrix,
    per_site_purity, fit_alpha,
)
from c1_bilinearity_test import dicke_state
from eq018_kernel_extract import build_decomps, propagate, measure_c1

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_kernel_bilin_probe"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


DJ = 0.01


def measure_c1_raw_op(op, decomps, times, N):
    """Measure c_1 for an arbitrary Hermitian operator (not necessarily
    a valid density matrix). Uses Tr(rho_i^2) fit as in measure_c1."""
    return measure_c1(op, decomps, times, N)


def test_bilinearity_on_scaled_coherent_state(N, bond, decomps, times):
    """Scan w in [0, pi/2] for |psi(w)> = cos(w)|vac> + sin(w)|S_1>
    and check bilinearity prediction."""
    S0 = dicke_state(N, 0)
    S1 = dicke_state(N, 1)

    # Reference measurements
    c_1_vac = measure_c1_raw_op(density_matrix(S0), decomps, times, N)["c_1"]
    c_1_S1 = measure_c1_raw_op(density_matrix(S1), decomps, times, N)["c_1"]
    c_1_mixed = measure_c1_raw_op(0.5 * (density_matrix(S0) + density_matrix(S1)),
                                   decomps, times, N)["c_1"]
    c_1_coherent = measure_c1_raw_op(density_matrix((S0 + S1)/np.sqrt(2)),
                                      decomps, times, N)["c_1"]
    K_DD_diag_0 = c_1_vac
    K_DD_diag_1 = c_1_S1
    K_DD_cross_01 = 2.0 * (c_1_mixed - 0.25 * (c_1_vac + c_1_S1))
    K_CC_01 = 2.0 * (c_1_coherent - c_1_mixed)

    print(f"  Reference measurements (at N={N}, bond={bond}):")
    print(f"    c_1[vac]      = {c_1_vac:+.5f}")
    print(f"    c_1[S_1]      = {c_1_S1:+.5f}")
    print(f"    c_1[mixed]    = {c_1_mixed:+.5f}")
    print(f"    c_1[coherent] = {c_1_coherent:+.5f}")
    print(f"    K_DD[0,0] = {K_DD_diag_0:+.5f}")
    print(f"    K_DD[1,1] = {K_DD_diag_1:+.5f}")
    print(f"    K_DD[0,1] = {K_DD_cross_01:+.5f}")
    print(f"    K_CC[0,1] = {K_CC_01:+.5f}")

    # Scaled state scan
    print(f"\n  Scaled state scan: rho(w) = |psi(w)><psi(w)|,"
          f" |psi(w)> = cos(w)|vac> + sin(w)|S_1>")
    print(f"  If strictly bilinear:")
    print(f"    c_1(w) = cos^4(w) K_DD[0,0] + sin^4(w) K_DD[1,1]")
    print(f"            + 2 cos^2(w) sin^2(w) K_DD[0,1]")
    print(f"            + 2 cos^2(w) sin^2(w) K_CC[0,1]")

    w_values = np.linspace(0, np.pi/2, 11)
    results = []
    print(f"  {'w/pi':>8} {'c_1 meas':>12} {'c_1 pred':>12} {'residual':>12}")
    for w in w_values:
        c = np.cos(w); s = np.sin(w)
        psi = c * S0 + s * S1
        rho_w = density_matrix(psi)
        c_1_w = measure_c1_raw_op(rho_w, decomps, times, N)["c_1"]
        # Bilinear prediction
        # rho_w components: rho^(0,0) = c^2, rho^(1,1) = s^2, rho^(0,1) = c*s (coh)
        # c_1 = (c^2)^2 K[0,0] + (s^2)^2 K[1,1] + 2*(c^2)*(s^2)*K[0,1] (DD cross)
        #     + 2 * (c*s)^2 * K_CC[0,1] (coherence self-pair, factor of 2 from hermitian)
        c_1_pred = (c**4) * K_DD_diag_0 + (s**4) * K_DD_diag_1 + \
                   2 * (c**2) * (s**2) * K_DD_cross_01 / 2 + \
                   2 * (c**2) * (s**2) * K_CC_01 / 2
        # Hmm wait, need to think about bilinear normalization.
        # For rho = rho_A |A><A| + rho_B |B><B| + rho_AB |A><B| + rho_BA |B><A|,
        # where the rho_XX are amplitudes (complex) and rho_AB = (rho_BA)^*,
        # bilinear form c_1 = sum K_{(XY),(UV)} rho_XY rho_UV with K indexed by pairs
        # of sector-blocks. For our states: rho_A = c^2, rho_B = s^2, rho_AB = c*s.
        # K_DD_diag_0 is c_1 measured on rho = |vac><vac| => rho_A = 1 => K[(0,0)(0,0)]=c_1_vac
        # K_DD_diag_1 is c_1 measured on rho = |S_1><S_1| => rho_B = 1 => K[(1,1)(1,1)]=c_1_S1
        # K_DD_cross_01 is 2*K_{(0,0)(1,1)} (defined so bilinear on diag gives cross = (2*ρ_A*ρ_B*K_{(0,0)(1,1)})/something)
        #
        # Actually: mixed state has rho_A = rho_B = 1/2. Bilinear:
        #   c_1[mix] = K[(0,0)(0,0)] * (1/2)^2 + K[(1,1)(1,1)] * (1/2)^2 + 2*K[(0,0)(1,1)]*(1/2)*(1/2)
        #           = 1/4 * K_DD_diag_0 + 1/4 * K_DD_diag_1 + 1/2 * K[(0,0)(1,1)]
        #   So K[(0,0)(1,1)] = 2*(c_1[mix] - 1/4*K_DD_diag_0 - 1/4*K_DD_diag_1)
        #   = K_DD_cross_01 as defined above. OK.
        # So K[(0,0)(1,1)] = K_DD_cross_01 in our notation.
        #
        # For rho(w), rho_{(0,0)} = c^2, rho_{(1,1)} = s^2, rho_{(0,1)} = c*s = rho_{(1,0)}.
        # Bilinear full:
        #   c_1 = (c^2)^2 K_{(0,0)(0,0)} + (s^2)^2 K_{(1,1)(1,1)}
        #       + 2*(c^2)(s^2) K_{(0,0)(1,1)}          [diagonal-cross]
        #       + 2*(c*s)(c*s) K_{(0,1)(1,0)}          [coherence self-pair, factor 2 from (0,1)(1,0)+(1,0)(0,1)]
        #       (assuming K_{(0,1)(1,0)} = K_{(1,0)(0,1)} by hermiticity of K)
        # K_CC[0,1] in our notation = K_{(0,1)(1,0)} * (some factor).
        # Let me recheck: K_CC[0,1] was defined as 2*(c_1[coh]-c_1[mix]).
        # c_1[coh] - c_1[mix] = contributions from coherence-coherence ONLY:
        #   = 2*(rho_{(0,1)})(rho_{(1,0)}) * K_{(0,1)(1,0)}
        #   = 2*(1/2)*(1/2) * K_{(0,1)(1,0)} = (1/2)*K_{(0,1)(1,0)}
        # So K_{(0,1)(1,0)} = 2*(c_1[coh]-c_1[mix]) = K_CC[0,1] as defined. Good.
        #
        # For rho(w), coh contribution:
        #   2*(c*s)(c*s) * K_{(0,1)(1,0)} = 2*c^2*s^2 * K_CC[0,1]
        # Full:
        #   c_1(w) = c^4*K_DD_diag_0 + s^4*K_DD_diag_1 + 2*c^2*s^2*K_DD_cross_01 + 2*c^2*s^2*K_CC[0,1]
        c_1_pred = (c**4) * K_DD_diag_0 + (s**4) * K_DD_diag_1 + \
                   2 * (c**2) * (s**2) * K_DD_cross_01 + \
                   2 * (c**2) * (s**2) * K_CC_01
        residual = c_1_w - c_1_pred
        results.append({"w": float(w), "c_1_measured": c_1_w,
                        "c_1_predicted": c_1_pred, "residual": residual})
        print(f"  {w/np.pi:>8.4f} {c_1_w:>+12.5f} {c_1_pred:>+12.5f} "
              f"{residual:>+12.2e}")

    max_residual = max(abs(r["residual"]) for r in results)
    print(f"\n  Max |residual|: {max_residual:.2e}")
    return {
        "references": {
            "K_DD_diag_0": K_DD_diag_0, "K_DD_diag_1": K_DD_diag_1,
            "K_DD_cross_01": K_DD_cross_01, "K_CC_01": K_CC_01
        },
        "scan": results,
        "max_residual": max_residual,
    }


def test_direct_coherence_probe(N, bond, decomps, times):
    """Direct probe of the 'pure coherence operator'
    rho_coh_only = (|vac><S_1| + |S_1><vac|)/2.

    This is NOT a valid density matrix (trace = 0, has negative
    eigenvalues). But Tr(rho_i^2) is still computable as a scalar
    function of t, and fit_alpha can be applied to it. The resulting
    c_1 is the 'direct' coherence-only c_1, which by bilinearity
    should equal K_CC[0,1] / 2 (the factor of 2 from the hermitian
    combination)."""
    S0 = dicke_state(N, 0)
    S1 = dicke_state(N, 1)
    # rho_coh_only = (|vac><S_1| + |S_1><vac|)/2 (Hermitian, trace 0)
    rho_coh_only = 0.5 * (np.outer(S0, S1.conj()) + np.outer(S1, S0.conj()))
    print(f"  rho_coh_only: trace = {np.trace(rho_coh_only):.2e} (should be ~0)")
    print(f"  rho_coh_only: hermitian residual = "
          f"{np.max(np.abs(rho_coh_only - rho_coh_only.conj().T)):.2e}")
    # Propagate and measure c_1
    r = measure_c1_raw_op(rho_coh_only, decomps, times, N)
    print(f"  c_1[rho_coh_only] = {r['c_1']:+.5f}")
    print(f"  alpha_plus: {['%.4f' % a for a in r['alpha_plus']]}")
    print(f"  alpha_minus: {['%.4f' % a for a in r['alpha_minus']]}")
    print(f"  RMSE max: {r['rmse_max']:.2e}")
    print(f"  alpha_bound_hit: {r['alpha_bound_hit']}")
    return r


def main():
    N = 5
    bond = 0
    print("=" * 78)
    print(f"EQ-018 bilinearity probe at N={N}, bond={bond}")
    print("=" * 78)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ}")

    t0 = time.time()
    decomps = build_decomps(N, bond)
    print(f"  L_A, L_B+/- built in {time.time() - t0:.1f} s")
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    print(f"\n--- Test 1: scaled-coherent-state bilinearity check ---")
    bilin_test = test_bilinearity_on_scaled_coherent_state(N, bond, decomps, times)

    print(f"\n--- Test 2: direct coherence operator probe ---")
    print(f"  If bilinearity holds STRICTLY, c_1[rho_coh_only] should equal")
    print(f"  K_CC[0,1] / 2 = {bilin_test['references']['K_CC_01']/2:+.5f}")
    direct_probe = test_direct_coherence_probe(N, bond, decomps, times)
    ratio_direct_to_half_KCC = direct_probe["c_1"] / (bilin_test["references"]["K_CC_01"] / 2) \
        if abs(bilin_test["references"]["K_CC_01"]) > 1e-10 else None
    if ratio_direct_to_half_KCC is not None:
        print(f"\n  Ratio c_1[direct] / (K_CC[0,1]/2) = {ratio_direct_to_half_KCC:.6f}")
        print(f"  (1.0 = strict bilinearity; deviation = cross-block alpha-fit mixing)")

    # Save
    out = {
        "N": N, "bond": bond, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ,
        "bilinearity_test": bilin_test,
        "direct_probe": direct_probe,
        "ratio_direct_to_half_KCC": ratio_direct_to_half_KCC,
    }
    path = RESULTS_DIR / "bilin_probe.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
