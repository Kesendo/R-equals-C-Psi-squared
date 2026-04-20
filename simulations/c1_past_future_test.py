#!/usr/bin/env python3
"""c1_past_future_test.py

Test of the hypothesis: single-excitation modes psi_k split into
past-like (E<0, k < (N+1)/2) and future-like (E>0, k > (N+1)/2) under
the Pi = time reversal operator, with the zero mode psi_{(N+1)/2} (odd N
only) as the Standing Wave middle.

For N=7, H eigenvalues in single-excitation sector are:
  E_k = 2 * cos(pi * k / 8) for k = 1..7
  = {+1.848, +1.414, +0.765,  0,  -0.765, -1.414, -1.848}  (for k=1..7)

Wait: cos(pi/8) > 0 so E_1 = +1.848 is the HIGHEST energy, E_7 = -1.848
the LOWEST (using the F65 convention k=1..N with sin-amplitudes).

Actually F65 gives amplitudes sqrt(2/(N+1)) * sin(pi k (i+1)/(N+1))
with H_1 tridiagonal having eigenvalues 2*J*cos(pi k / (N+1)).
For N=7 and J=1: E_k = 2*cos(pi*k/8), so
  k=1: +2cos(pi/8) = +1.848  (bonding: lowest energy in AFM convention,
                             but with +J convention this is HIGHEST)
  ...
We test: for antiferromagnetic XY (+J (XX+YY)/2), the spectrum is
"inverted" compared to the usual bonding/antibonding naming.

The PROJECT's "bonding" per PTF and PRIMORDIAL_QUBIT is psi_1 (smallest
k), which has E = +1.848 in our J=+1 convention. PTF calls this the
"bonding mode" which is a QUANTUM-CHEMISTRY term; in condensed-matter
language with +J it would be "anti-bonding" (highest antisymmetric
combination).

So our past/future interpretation inverts:
  k=1..3:  positive E -> "future-like" (retarded) in the J>0 convention
  k=4:     E = 0 -> Standing Wave MIDDLE
  k=5..7:  negative E -> "past-like" (advanced)

But this only matters for direction. The symmetry around k=(N+1)/2 is what
we really test.

Predictions:
  1. c_1(psi_4+vac) is structurally special (exactly 0 or maximum).
  2. c_1(psi_k+vac) and c_1(psi_{N+1-k}+vac) are related:
     - identical if Pi-symmetric initial state
     - opposite if Pi-antisymmetric
     - or something else entirely
  3. ψ_4 has zero amplitude at sites 1, 3, 5, 7 (nodes), so dephasing
     on these sites does not couple. Under bond (0,1) perturbation,
     which touches sites 0 and 1, psi_4 has amplitude at site 0
     (nonzero) but NOT site 1 (node). So the perturbation "sees" only
     half the usual bonding-mode amplitude -> smaller c_1 expected.
"""
from __future__ import annotations

import json
import sys
import time
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

RESULTS_DIR = Path(__file__).parent / "results" / "c1_past_future_test"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

N = 7
DJ_EXTRACT = 0.01
DJ_PROBE = [-DJ_EXTRACT, +DJ_EXTRACT]


def v_effect(N):
    return 1.0 + np.cos(np.pi / N)


def bonding_plus_vacuum(N, k):
    v = vacuum_ket(N)
    psi = single_excitation_mode(N, k=k)
    return (v + psi) / np.sqrt(2.0)


def single_excitation_H_spectrum(N, J):
    """Energies of the N single-excitation modes of the open XY chain."""
    return np.array([2.0 * J * np.cos(np.pi * k / (N + 1)) for k in range(1, N + 1)])


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


def main():
    print("="*70)
    print(f"c_1 past/future test at N = {N}")
    print("="*70)
    d = 2**N
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}")
    print(f"  V({N}) = {v_effect(N):.4f}")
    print(f"  Prediction for psi_1+vac: c_1 ~ 0.5 * V(N) = {0.5 * v_effect(N):.4f}")
    print(f"  EQ-014 reports c_1(psi_1) = 0.97 at N=7 (for comparison)")

    Es = single_excitation_H_spectrum(N, J_UNIFORM)
    print(f"\n  Single-excitation H-spectrum at N={N}, J=1:")
    for k in range(1, N + 1):
        label = ""
        if k < (N + 1) / 2:
            label = "future-like (positive E in J>0)"
        elif k > (N + 1) / 2:
            label = "past-like (negative E in J>0)"
        else:
            label = "zero mode / Standing Wave"
        print(f"    psi_{k}: E = {Es[k-1]:+.4f}  [{label}]")

    # ---- Heavy lifting: eigendecompose L_A, L_B+ and L_B- at N=7 ----
    J_A = [J_UNIFORM] * (N - 1)
    J_B_plus = list(J_A); J_B_plus[0] = J_UNIFORM + DJ_EXTRACT
    J_B_minus = list(J_A); J_B_minus[0] = J_UNIFORM - DJ_EXTRACT

    print(f"\n  Building Liouvillians (d^2 = {4**N} = 16384)...")
    t0 = time.time()
    H_A = build_H_XY(J_A, N)
    L_A = build_liouvillian_matrix(H_A, GAMMA_0, N)
    print(f"  L_A built in {time.time()-t0:.1f} s")

    t0 = time.time()
    ev_A, VR_A, VLinv_A = eig_and_inv(L_A)
    print(f"  L_A eigendecomp: {time.time()-t0:.1f} s")
    # Free H_A memory
    del H_A

    t0 = time.time()
    H_B_plus = build_H_XY(J_B_plus, N)
    L_B_plus = build_liouvillian_matrix(H_B_plus, GAMMA_0, N)
    ev_Bp, VR_Bp, VLinv_Bp = eig_and_inv(L_B_plus)
    del H_B_plus, L_B_plus
    print(f"  L_B+ eigendecomp: {time.time()-t0:.1f} s")

    t0 = time.time()
    H_B_minus = build_H_XY(J_B_minus, N)
    L_B_minus = build_liouvillian_matrix(H_B_minus, GAMMA_0, N)
    ev_Bm, VR_Bm, VLinv_Bm = eig_and_inv(L_B_minus)
    del H_B_minus, L_B_minus
    print(f"  L_B- eigendecomp: {time.time()-t0:.1f} s")

    # Also free L_A now (we keep the decomp)
    del L_A

    # ---- Per-state loop ----
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    results = {}
    print(f"\n  Running through states (k = 1..{N}):")
    print(f"  {'k':>2} {'E_k':>10} {'role':>25} {'c_1':>10} {'alpha_i per site':>60}")
    for k in range(1, N + 1):
        t0 = time.time()
        phi = bonding_plus_vacuum(N, k)
        rho_0 = density_matrix(phi)
        rho_A = propagate(ev_A, VR_A, VLinv_A, rho_0, times)
        rho_Bp = propagate(ev_Bp, VR_Bp, VLinv_Bp, rho_0, times)
        rho_Bm = propagate(ev_Bm, VR_Bm, VLinv_Bm, rho_0, times)
        P_A = per_site_purity(rho_A, N)
        P_Bp = per_site_purity(rho_Bp, N)
        P_Bm = per_site_purity(rho_Bm, N)

        alpha_p = np.zeros(N); alpha_m = np.zeros(N)
        for i in range(N):
            a, _ = fit_alpha(times, P_A[:, i], P_Bp[:, i])
            alpha_p[i] = a
            a, _ = fit_alpha(times, P_A[:, i], P_Bm[:, i])
            alpha_m[i] = a
        closure_p = float(np.sum(np.log(alpha_p)))
        closure_m = float(np.sum(np.log(alpha_m)))
        c_1 = (closure_p - closure_m) / (2 * DJ_EXTRACT)

        Ek = Es[k - 1]
        if k < (N + 1) / 2:
            role = "future (E>0)"
        elif k > (N + 1) / 2:
            role = "past (E<0)"
        else:
            role = "STANDING WAVE (E=0)"
        # Show alpha at dJ=+0.01 only for readability
        alpha_str = "  ".join(f"{a:.3f}" for a in alpha_p)
        print(f"  {k:>2d} {Ek:>+10.4f} {role:>25} {c_1:>+10.4f}  {alpha_str}")

        results[f"psi_{k}"] = {
            "k": k, "E_k": float(Ek), "role": role, "c_1": float(c_1),
            "closure_plus": closure_p, "closure_minus": closure_m,
            "alpha_plus": alpha_p.tolist(), "alpha_minus": alpha_m.tolist(),
            "time_per_state_s": time.time() - t0,
        }

    # ---- Pi-pair analysis ----
    print(f"\n  Pi-pair analysis: c_1(psi_k) vs c_1(psi_{{N+1-k}})")
    print(f"  {'k':>2} {'c_1(k)':>10} {'c_1(N+1-k)':>12} {'sum':>10} {'diff':>10} {'ratio':>10}")
    pi_pair_data = {}
    for k in range(1, N + 1):
        k_mirror = N + 1 - k
        if k > k_mirror:
            continue
        c_k = results[f"psi_{k}"]["c_1"]
        c_m = results[f"psi_{k_mirror}"]["c_1"]
        s = c_k + c_m
        d = c_k - c_m
        r = c_k / c_m if abs(c_m) > 1e-10 else float('inf')
        note = ""
        if k == k_mirror:
            note = "  (self-mirror: STANDING WAVE)"
        print(f"  {k:>2d}<->{k_mirror:<2d} {c_k:>+10.4f} {c_m:>+12.4f} "
              f"{s:>+10.4f} {d:>+10.4f} {r:>10.4f}{note}")
        pi_pair_data[f"pair_{k}_{k_mirror}"] = {
            "c_1_k": c_k, "c_1_mirror": c_m, "sum": s, "diff": d, "ratio": r,
        }

    # ---- Verdict ----
    print(f"\n  Verdict checks:")
    c_standing = results.get(f"psi_{(N+1)//2}", {}).get("c_1")
    if c_standing is not None:
        print(f"    psi_{(N+1)//2} (standing wave, E=0):  c_1 = {c_standing:+.4f}")
        if abs(c_standing) < 0.01:
            print(f"    VERDICT: Standing wave is closure-PROTECTED (c_1 ~ 0)")
        elif abs(c_standing) > 1.0:
            print(f"    VERDICT: Standing wave is MAXIMALLY sensitive (c_1 large)")
        else:
            print(f"    VERDICT: Standing wave is intermediate, not obviously special")

    # Save
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "V_N": v_effect(N),
        "H_spectrum": Es.tolist(),
        "results_by_state": results,
        "pi_pair_analysis": pi_pair_data,
        "prediction_psi_1": 0.5 * v_effect(N),
        "eq014_psi_1_match": 0.97,
    }
    path = RESULTS_DIR / "past_future_test.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
