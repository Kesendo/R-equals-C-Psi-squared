#!/usr/bin/env python3
"""
EQ-024 refinement: Brecher test at N=6 (N-scaling of receiver-engineering)
===========================================================================
Follow-up to eq024_refinement_shadow_lens_broken.py (N=5). Tests whether the
receiver-engineering advantage over gamma-Sacrifice-Zone scales.

At N=5:
- |+-+-+> + uniform gamma + uniform J:       Peak SumMI = 1.32
- |+-+-+> + uniform gamma + extreme J:       Peak SumMI = 3.30
- gamma-Sacrifice at |+>^5 (RESONANT_RETURN): Peak SumMI = 0.230

Receiver engineering beats gamma-Sacrifice absolute by 5.7x to 14x at N=5.

At N=6:
- F71 mirrors: sites (0,5), (1,4), (2,3)
- 5 bonds
- RESONANT_RETURN scales: N=5 formula Peak SumMI = 0.230, N=7 = 0.408, N=9 = 0.619
  (interpolating N=6 is roughly 0.3 by their scaling pattern)
- We test: |+-++-+> (SU(2)-breaking F71-symmetric analog of |+-+-+>)
           |010010> (SU(2)-breaking F71-symmetric analog of |01010>)

Reduced scope for tractability at d^2 = 4096:
- 4 J-profiles
- 3 initial states (including |+>^6 control)
- 20 t-points in [0.1, 15]
- Expected runtime ~30-60 min

Outputs: simulations/results/eq024_refinement/brecher_n6.{txt,json}
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
OUT_TXT = os.path.join(OUT_DIR, "brecher_n6.txt")
OUT_JSON = os.path.join(OUT_DIR, "brecher_n6.json")

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

N = 6
d = 2 ** N   # 64
d2 = d * d   # 4096
N_BONDS = N - 1  # 5


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


log(f"Pre-computing {N_BONDS} Heisenberg bond blocks at N={N} (d={d}) ...")
_t_pre = _time.time()
_BOND_BLOCKS = []
for b in range(N_BONDS):
    block = (site_op(sx, b) @ site_op(sx, b + 1)
             + site_op(sy, b) @ site_op(sy, b + 1)
             + site_op(sz, b) @ site_op(sz, b + 1))
    _BOND_BLOCKS.append(block)
log(f"Bond blocks ready in {_time.time() - _t_pre:.1f} s")


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


def evolve(L, rho0_flat, t):
    v = expm(L * t) @ rho0_flat
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
# J-profiles (5 bonds at N=6), reduced set
# ------------------------------------------------------------------
J_PROFILES = {
    "Uniform J=1":                          [1.0, 1.0, 1.0, 1.0, 1.0],
    "Cut bond 2-3 (eps=0.01, center)":      [1.0, 1.0, 0.01, 1.0, 1.0],
    "Cut bond 0-1 (eps=0.01, edge)":        [0.01, 1.0, 1.0, 1.0, 1.0],
    "Strong-weak alternating":              [5.0, 0.2, 5.0, 0.2, 5.0],
}

# ------------------------------------------------------------------
# Initial states: |+>^6 control + 2 SU(2)-breaking F71-symmetric
# At N=6, F71 mirror pairs are (0,5), (1,4), (2,3).
# |+-++-+>: (+, -, +, +, -, +). sites 0=5=+, 1=4=-, 2=3=+. F71-sym.
# |010010>: sites 0=5=|0>, 1=4=|1>, 2=3=|0>. F71-sym.
# ------------------------------------------------------------------
INITIAL_STATES = {
    "|+>^6 (control, Class 3 J-blind)":
        kron_all([plus] * N),
    "|+-++-+> (N=6 F71-sym analog of |+-+-+>)":
        kron_all([plus, minus, plus, plus, minus, plus]),
    "|010010> (N=6 F71-sym analog of |01010>)":
        _state_from_bits([0, 1, 0, 0, 1, 0]),
}


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
if __name__ == "__main__":
    t_start = _time.time()

    log()
    log("EQ-024 refinement: Brecher-Test at N=6 (N-scaling)")
    log("=" * 72)
    log(f"N = {N}, d^2 = {d2}, gamma_0 = 0.05 uniform, N_BONDS = {N_BONDS}")
    log()
    log("Reduced scope: 4 J-profiles x 3 initial states x 20 t-points.")
    log("Expected runtime: ~30-60 min.")
    log()

    gamma_uniform = [0.05] * N
    t_scan = np.linspace(0.1, 15.0, 20)

    all_results = {}

    for init_name, psi_init in INITIAL_STATES.items():
        log("=" * 72)
        log(f"INITIAL STATE: {init_name}")
        log("=" * 72)

        rho0 = np.outer(psi_init, psi_init.conj())
        rho0 = (rho0 + rho0.conj().T) / 2
        rho0 /= np.trace(rho0).real
        rho0_flat = rho0.flatten()

        init_results = {}
        uniform_peak = None

        log(f"  {'J-profile':<35}  {'Peak SumMI':>11}  {'Peak t':>7}  {'vs Uniform':>11}  {'time':>8}")
        log("  " + "-" * 82)

        for prof_name, J_vec in J_PROFILES.items():
            t0 = _time.time()
            H = build_H(J_vec)
            L = build_L(H, gamma_uniform)

            peak_smi = 0.0
            peak_t_val = 0.0
            for t in t_scan:
                rho_t = evolve(L, rho0_flat, t)
                smi = sum_mi_adjacent(rho_t)
                if smi > peak_smi:
                    peak_smi = smi
                    peak_t_val = t

            dt_walltime = _time.time() - t0

            if prof_name == "Uniform J=1":
                uniform_peak = peak_smi

            ratio = peak_smi / uniform_peak if uniform_peak and uniform_peak > 1e-10 else float("nan")

            log(f"  {prof_name:<35}  {peak_smi:11.4f}  {peak_t_val:7.2f}  "
                f"{ratio:10.2f}x  {dt_walltime:7.1f}s")

            init_results[prof_name] = {
                "J_vec": J_vec,
                "peak_sum_mi": float(peak_smi),
                "peak_t": float(peak_t_val),
                "ratio_to_uniform": float(ratio) if not np.isnan(ratio) else None,
                "compute_time_s": float(dt_walltime),
            }

        all_results[init_name] = init_results
        log()

    # ------------------------------------------------------------------
    # Cross-receiver summary
    # ------------------------------------------------------------------
    log("=" * 72)
    log("CROSS-RECEIVER SUMMARY AT N=6")
    log("=" * 72)
    log()
    log(f"Reference: gamma-Sacrifice-Zone at |+>^N under non-uniform gamma:")
    log(f"  N=5: Peak SumMI = 0.230 (RESONANT_RETURN Test 8)")
    log(f"  N=7: Peak SumMI = 0.408 (RESONANT_RETURN Test 8)")
    log(f"  (N=6 not directly tabulated; interpolate ~0.3)")
    log()
    log(f"  {'Receiver':<50}  {'Uniform-J SumMI':>15}  {'Best J SumMI':>13}  {'Boost':>7}")
    log("  " + "-" * 94)
    for init_name, init_res in all_results.items():
        uniform_peak_val = init_res["Uniform J=1"]["peak_sum_mi"]
        best_peak = max(r["peak_sum_mi"] for r in init_res.values())
        if uniform_peak_val > 1e-10:
            boost = best_peak / uniform_peak_val
            boost_str = f"{boost:.2f}x"
        elif best_peak > 1e-10:
            boost_str = "inf"
        else:
            boost_str = "n/a"
        log(f"  {init_name:<50}  {uniform_peak_val:15.4f}  {best_peak:13.4f}  {boost_str:>7}")
    log()

    log("N=5 comparison (from shadow_lens_broken.py, commit bf080a3):")
    log(f"  |+>^5       : uniform J = 0.00,   best J = 0.00,   boost = n/a")
    log(f"  |+-+-+>     : uniform J = 1.32,   best J = 3.30,   boost = 2.50x")
    log(f"  |01010>     : uniform J = 1.38,   best J = 3.38,   boost = 2.46x")
    log()

    total_t = _time.time() - t_start
    log(f"Total runtime: {total_t:.1f} s ({total_t/60:.1f} min)")

    out = {
        "meta": {
            "N": N,
            "gamma_0": 0.05,
            "t_scan": [float(t) for t in t_scan],
            "J_profiles": {k: v for k, v in J_PROFILES.items()},
            "initial_states": list(INITIAL_STATES.keys()),
            "total_runtime_s": total_t,
        },
        "results": all_results,
        "reference_n5": {
            "|+>^5": {"uniform_J_peak_smi": 0.0, "best_J_peak_smi": 0.0, "boost": None},
            "|+-+-+>": {"uniform_J_peak_smi": 1.32, "best_J_peak_smi": 3.30, "boost": 2.50},
            "|01010>": {"uniform_J_peak_smi": 1.38, "best_J_peak_smi": 3.38, "boost": 2.46},
        },
        "reference_resonant_return": {
            "N5_sacrifice_peak_smi": 0.230,
            "N7_sacrifice_peak_smi": 0.408,
            "source": "RESONANT_RETURN Test 8",
        },
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fj:
        json.dump(out, fj, indent=2)
    log(f"JSON: {OUT_JSON}")
    log(f"Text: {OUT_TXT}")
    _outf.close()
