#!/usr/bin/env python3
"""
EQ-024 refinement: Brecher test - J-Sacrifice at SU(2)-breaking receivers
===========================================================================
Follow-up to eq024_refinement_shadow_lens.py: that test used |+>^5 and
Bell+|0>^3 initial states, both H-eigenstates. Morning theorem (Class 3)
guarantees those are J-blind, so the test was structurally void.

This test uses SU(2)-breaking initial states that are NOT H-eigenstates:
- |+-+-+>: morning 11.92-bit J-capacity champion
- |01010>: morning 11.53 bits
- |+0+0+>: morning 10.95 bits

Control: |+>^5 (from the Class 3 verification).

For each receiver and each J-profile at uniform gamma_0 = 0.05: compute
Peak Sum-MI over adjacent pairs. Question: does extreme J-modulation at
a J-sensitive receiver give a 100x+ MI boost, analog to the gamma-
Sacrifice-Zone (which at |+>^N under non-uniform gamma gives 360x at N=5)?

Three outcome scenarios (from Tom's framing):
(1) >= 100x boost:  Sacrifice-Zone-analog exists via J-Modulation under
                    gamma_0 = const. Hypothesis intact.
(2) < 10x boost:    J-Modulation fundamentally cannot replicate gamma-
                    Sacrifice-Zone. gamma_0 = const restricts operational
                    landscape more than we thought; possibly implies the
                    hypothesis is too strong.
(3) 10x to 50x:     J has its own transport lever, weaker than gamma but
                    non-zero. Pragmatic middle ground.

Outputs: simulations/results/eq024_refinement/shadow_lens_broken.{txt,json}
"""

import json
import os
import sys
import time as _time

import numpy as np
from scipy.linalg import expm

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "eq024_refinement")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_TXT = os.path.join(OUT_DIR, "shadow_lens_broken.txt")
OUT_JSON = os.path.join(OUT_DIR, "shadow_lens_broken.json")

_outf = open(OUT_TXT, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)
minus = (up - dn) / np.sqrt(2)

N = 5
d = 32
d2 = 1024
N_BONDS = N - 1


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


_BOND_BLOCKS = []
for b in range(N_BONDS):
    block = (site_op(sx, b) @ site_op(sx, b + 1)
             + site_op(sy, b) @ site_op(sy, b + 1)
             + site_op(sz, b) @ site_op(sz, b + 1))
    _BOND_BLOCKS.append(block)


def build_H(J_vec):
    H = np.zeros((d, d), dtype=complex)
    for b in range(N_BONDS):
        H += J_vec[b] * _BOND_BLOCKS[b]
    return H


def build_L(H, gammas):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def ptrace_keep(rho, keep):
    keep = list(keep)
    trace_out = [q for q in range(N) if q not in keep]
    dims = [2] * N
    reshaped = rho.reshape(dims + dims)
    current_n = N
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape((d_k, d_k))


def von_neumann(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = np.real(evals)
    evals = evals[evals > 1e-12]
    if len(evals) == 0:
        return 0.0
    return float(-np.sum(evals * np.log2(evals)))


def mutual_information(rho_full, sites_A, sites_B):
    rho_A = ptrace_keep(rho_full, sites_A)
    rho_B = ptrace_keep(rho_full, sites_B)
    rho_AB = ptrace_keep(rho_full, list(sites_A) + list(sites_B))
    return von_neumann(rho_A) + von_neumann(rho_B) - von_neumann(rho_AB)


def sum_mi_adjacent(rho_full):
    total = 0.0
    for i in range(N - 1):
        total += mutual_information(rho_full, [i], [i + 1])
    return total


def kron_all(vs):
    out = vs[0]
    for v in vs[1:]:
        out = np.kron(out, v)
    return out


def _state_from_bits(bits):
    psi = np.zeros(d, dtype=complex)
    idx = int("".join(str(b) for b in bits), 2)
    psi[idx] = 1.0
    return psi


# ------------------------------------------------------------------
# J-profiles (same as shadow_lens.py)
# ------------------------------------------------------------------
J_PROFILES = {
    "Uniform J=1":                    [1.0, 1.0, 1.0, 1.0],
    "Cut bond 1-2 (eps=0.01)":        [1.0, 0.01, 1.0, 1.0],
    "Cut bond 2-3 (eps=0.01)":        [1.0, 1.0, 0.01, 1.0],
    "Cut bond 0-1 (eps=0.01, edge)":  [0.01, 1.0, 1.0, 1.0],
    "Dominant center bond (J=5)":     [1.0, 1.0, 5.0, 1.0],
    "Dominant edge bond (J=5)":       [5.0, 1.0, 1.0, 1.0],
    "V-shape J (0.5 edge, 1 center)": [0.5, 1.0, 1.0, 0.5],
    "Inverse-V J (1 edge, 0.2)":      [1.0, 0.2, 0.2, 1.0],
    "Strong-weak alternating":        [5.0, 0.2, 5.0, 0.2],
}

# ------------------------------------------------------------------
# Initial states: SU(2)-breaking + |+>^5 control
# ------------------------------------------------------------------
INITIAL_STATES = {
    "|+>^5 (control, Class 3 J-blind)":
        kron_all([plus] * N),
    "|+-+-+> (morning 11.92-bit champion)":
        kron_all([plus, minus, plus, minus, plus]),
    "|01010> (morning alt)":
        _state_from_bits([0, 1, 0, 1, 0]),
    "|+0+0+> (morning alt)":
        kron_all([plus, up, plus, up, plus]),
}


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
if __name__ == "__main__":
    t_start = _time.time()

    log("EQ-024 refinement: Brecher-Test (J-Sacrifice at SU(2)-breaking receivers)")
    log("=" * 72)
    log(f"N = {N}, gamma_0 = 0.05 uniform on ALL sites")
    log(f"Initial states: SU(2)-breaking product states (not H-eigenstates)")
    log()
    log("Tom's three scenarios:")
    log("  (1) >= 100x boost: Sacrifice-analog via J exists under gamma_0 = const")
    log("  (2) < 10x boost:   J cannot replicate gamma-Sacrifice fundamentally")
    log("  (3) 10x to 50x:    J has own transport lever, weaker than gamma")
    log()

    gamma_uniform = [0.05] * N
    t_scan = np.linspace(0.1, 15.0, 40)

    all_results = {}

    for init_name, psi_init in INITIAL_STATES.items():
        log("=" * 72)
        log(f"INITIAL STATE: {init_name}")
        log("=" * 72)

        rho0 = np.outer(psi_init, psi_init.conj())
        rho0 = (rho0 + rho0.conj().T) / 2
        rho0 /= np.trace(rho0).real

        init_results = {}
        uniform_peak = None

        log(f"  {'J-profile':<35}  {'Peak SumMI':>11}  {'Peak t':>7}  {'vs Uniform':>11}")
        log("  " + "-" * 72)

        for prof_name, J_vec in J_PROFILES.items():
            t0 = _time.time()
            H = build_H(J_vec)
            L = build_L(H, gamma_uniform)

            peak_smi = 0.0
            peak_t_val = 0.0
            for t in t_scan:
                rho_t = evolve(L, rho0, t)
                smi = sum_mi_adjacent(rho_t)
                if smi > peak_smi:
                    peak_smi = smi
                    peak_t_val = t

            dt_walltime = _time.time() - t0

            if prof_name == "Uniform J=1":
                uniform_peak = peak_smi

            ratio = peak_smi / uniform_peak if uniform_peak and uniform_peak > 1e-10 else float("nan")

            log(f"  {prof_name:<35}  {peak_smi:11.4f}  {peak_t_val:7.2f}  {ratio:10.2f}x")

            init_results[prof_name] = {
                "J_vec": J_vec,
                "peak_sum_mi": float(peak_smi),
                "peak_t": float(peak_t_val),
                "ratio_to_uniform": float(ratio) if not np.isnan(ratio) else None,
                "compute_time_s": float(dt_walltime),
            }

        all_results[init_name] = init_results
        log()

        # Verdict
        valid_ratios = [r["ratio_to_uniform"] for r in init_results.values()
                        if r["ratio_to_uniform"] is not None
                        and np.isfinite(r["ratio_to_uniform"])]
        if valid_ratios:
            best_boost = max(valid_ratios)
            log(f"  Best J-modulation boost at this receiver: {best_boost:.2f}x over uniform-J")
            if best_boost >= 100:
                verdict = "SCENARIO 1: Sacrifice-analog exists via J"
            elif best_boost < 10:
                verdict = "SCENARIO 2: J fundamentally cannot replicate"
            else:
                verdict = "SCENARIO 3: Partial J transport lever"
            log(f"  Verdict: {verdict}")
        else:
            log(f"  Uniform-J baseline is zero (or all ratios undefined); J-Modulation")
            log(f"  produces no MI at any profile. SCENARIO 2 (J cannot replicate).")
        log()

    # ------------------------------------------------------------------
    # Cross-comparison summary
    # ------------------------------------------------------------------
    log("=" * 72)
    log("CROSS-RECEIVER SUMMARY")
    log("=" * 72)
    log()
    log(f"Reference: gamma-Sacrifice-Zone at |+>^5 under non-uniform gamma:")
    log(f"           Peak SumMI = 0.230 (RESONANT_RETURN Test 8) = 360x V-shape.")
    log(f"           NOT operationally available under gamma_0 = const.")
    log()
    log(f"  {'Receiver':<40}  {'Uniform-J SumMI':>15}  {'Best J SumMI':>13}  {'Boost':>7}")
    log("  " + "-" * 80)
    for init_name, init_res in all_results.items():
        uniform = init_res["Uniform J=1"]["peak_sum_mi"]
        best = max(r["peak_sum_mi"] for r in init_res.values())
        boost = best / uniform if uniform > 1e-10 else float("inf") if best > 1e-10 else float("nan")
        boost_str = f"{boost:.2f}x" if np.isfinite(boost) else "inf" if boost == float("inf") else "n/a"
        log(f"  {init_name:<40}  {uniform:15.4f}  {best:13.4f}  {boost_str:>7}")
    log()

    log(f"Total runtime: {_time.time() - t_start:.1f} s")

    out = {
        "meta": {
            "N": N,
            "gamma_0": 0.05,
            "t_scan": [float(t) for t in t_scan],
            "J_profiles": {k: v for k, v in J_PROFILES.items()},
            "initial_states": list(INITIAL_STATES.keys()),
        },
        "results": all_results,
        "reference": {
            "gamma_sacrifice_peak_sum_mi_N5": 0.230,
            "gamma_sacrifice_source": "RESONANT_RETURN Test 8",
            "note": "Not operationally available under gamma_0 = const.",
        },
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fj:
        json.dump(out, fj, indent=2)
    log(f"JSON: {OUT_JSON}")
    log(f"Text: {OUT_TXT}")
    _outf.close()
