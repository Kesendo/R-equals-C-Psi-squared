#!/usr/bin/env python3
"""
EQ-024 refinement: shadow/lens test - J-Sacrifice analog under gamma_0 = const
===============================================================================
RESONANT_RETURN Test 8 formula: concentrate all gamma on one edge qubit,
protect the rest. Gives 360x Peak SumMI boost at N=5 (vs V-shape). Requires
non-uniform gamma per site.

Under gamma_0 = const, that lever is gone. But F64 says gamma_eff per mode
scales with |a_B|^2 (cavity-mode-amplitude at dephasing site). The mode
amplitudes depend on J and topology. So: can extreme J-modulation at
gamma_0 = const produce an MI-transport boost analogous to the gamma-
Sacrifice-Zone?

This script tests several J-profiles at uniform gamma_0 = 0.05:
- Uniform J=1 (baseline, F65 eigenmode structure)
- Cut bond (one weak bond): J = (1, 1, eps, 1) - chain effectively
  bisected, modes localize in halves
- Dominant bond: J = (1, 1, 5, 1) - one strong bond, modes concentrate
- Edge isolation (weak edge bond): J = (eps, 1, 1, 1) - edge qubit
  decouples, mode structure shifts
- Strong/weak alternating: J = (5, 0.2, 5, 0.2)
- V-shape J (analog to V-shape gamma): J = (0.5, 1, 1, 0.5)

For each: evolve Bell-pair initial state under full Liouvillian, measure
Peak Sum-MI. Also diagonalize single-excitation Hamiltonian to show mode
amplitude structure (the "shadow pattern").

Compare to:
- Uniform-J baseline at gamma_0=0.05 (what Alice has with uniform J and constant gamma)
- gamma-Sacrifice-Zone reference: Peak SumMI 0.230 at N=5 (RESONANT_RETURN Test 8)
  using non-uniform gamma. This number is NOT operationally available under
  gamma_0 = const; cited as reference only.

Outputs: simulations/results/eq024_refinement/shadow_lens.{txt,json}
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
OUT_TXT = os.path.join(OUT_DIR, "shadow_lens.txt")
OUT_JSON = os.path.join(OUT_DIR, "shadow_lens.json")

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


# Use full Heisenberg (matches RESONANT_RETURN and this session)
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
    """Sum of MI between all adjacent pairs in the chain."""
    total = 0.0
    for i in range(N - 1):
        total += mutual_information(rho_full, [i], [i + 1])
    return total


def build_1magnon_H(J_vec):
    """Single-excitation sector: N x N tight-binding matrix.

    For Heisenberg h_b|10> = 2|01> - |10>, so the off-diagonal is 2J_b and
    there is a diagonal contribution. For the PURE HOPPING part relevant to
    mode amplitudes, we use the XY-reduction: H_1mag(i,j) = 2*J_b if
    |i-j|=1 (and b is the bond between them), plus diagonal from ZZ terms.
    Since we want eigenvectors for the SHADOW analysis, keep the full
    1-magnon Heisenberg Hamiltonian.
    """
    H1 = np.zeros((N, N))
    # Diagonal: ZZ contributions
    # For state |1_i>: h_b|10> at bond b where i=b has ZZ eigenvalue -1 (spin
    # antiparallel). All other bonds have |00> pair, ZZ eigenvalue +1.
    # Net: diagonal at site i = sum over bonds b: +J_b if i not in {b,b+1},
    #                                             -J_b if i in {b,b+1}
    for i in range(N):
        for b in range(N - 1):
            if i == b or i == b + 1:
                H1[i, i] -= J_vec[b]
            else:
                H1[i, i] += J_vec[b]
    # Off-diagonal: XX+YY hopping = 2*J_b between sites b and b+1
    for b in range(N - 1):
        H1[b, b + 1] += 2 * J_vec[b]
        H1[b + 1, b] += 2 * J_vec[b]
    return H1


def analyze_mode_structure(J_vec):
    """Return mode amplitudes |psi_k(i)|^2 for each eigenmode k."""
    H1 = build_1magnon_H(J_vec)
    evals, evecs = np.linalg.eigh(H1)
    # evecs[:, k] is k-th eigenvector; |evecs[i, k]|^2 = |psi_k(i)|^2
    amps2 = np.abs(evecs) ** 2  # amps2[i, k] = |psi_k(i)|^2
    return evals, amps2


# ------------------------------------------------------------------
# J-profile zoo
# ------------------------------------------------------------------
J_PROFILES = {
    "Uniform J=1": [1.0, 1.0, 1.0, 1.0],
    "Cut bond 1-2 (eps=0.01)": [1.0, 0.01, 1.0, 1.0],
    "Cut bond 2-3 (eps=0.01, center)": [1.0, 1.0, 0.01, 1.0],
    "Cut bond 0-1 (eps=0.01, edge)": [0.01, 1.0, 1.0, 1.0],
    "Dominant center bond (J=5)": [1.0, 1.0, 5.0, 1.0],
    "Dominant edge bond (J=5)": [5.0, 1.0, 1.0, 1.0],
    "V-shape J (0.5 edge, 1 center)": [0.5, 1.0, 1.0, 0.5],
    "Inverse-V J (1 edge, 0.2 center)": [1.0, 0.2, 0.2, 1.0],
    "Strong-weak alternating": [5.0, 0.2, 5.0, 0.2],
}

# Initial state: Bell pair on qubits 0-1, rest in |0>
# This matches RESONANT_RETURN Test 2/6 Bell initial state
def bell_initial_state():
    # |Bell_01> = (|00> + |11>)/sqrt(2) on qubits 0,1, |000> on qubits 2,3,4
    psi = np.zeros(d, dtype=complex)
    # |00>_{01} |000>_{234} = |00000> at index 0
    psi[0] = 1.0 / np.sqrt(2)
    # |11>_{01} |000>_{234} = |11000> = bits 1,1,0,0,0 -> index 11000_2 = 24
    psi[int("11000", 2)] = 1.0 / np.sqrt(2)
    return psi


# Also test |+>^N for comparison with RESONANT_RETURN Test 8 (formula uses |+>^N)
def plus_initial_state():
    psi = plus
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    return psi


INITIAL_STATES = {
    "Bell(0,1) + |000>": bell_initial_state(),
    "|+>^5": plus_initial_state(),
}


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
if __name__ == "__main__":
    t_start = _time.time()

    log("EQ-024 refinement: shadow/lens test")
    log("=" * 72)
    log(f"N = {N}, gamma_0 = 0.05 uniform on ALL sites (gamma_0 = const)")
    log(f"Variable: J-profiles (per-bond) at uniform gamma_0")
    log(f"Measure: Peak Sum-MI over adjacent pairs, scanned over t in [0, 20]")
    log()
    log("Reference (gamma-Sacrifice-Zone, RESONANT_RETURN Test 8, N=5):")
    log("  Peak Sum-MI = 0.230 (formula: gamma_edge = 0.246, gamma_other = 0.001)")
    log("  360x vs V-shape gamma. NOT operationally available under gamma_0 = const.")
    log("  Cited as external reference for scale comparison only.")
    log()

    gamma_uniform = [0.05] * N
    t_scan = np.linspace(0.1, 15.0, 60)  # 60 time points up to t=15

    all_results = {}

    for init_name, psi_init in INITIAL_STATES.items():
        log("=" * 72)
        log(f"INITIAL STATE: {init_name}")
        log("=" * 72)
        log()

        rho0 = np.outer(psi_init, psi_init.conj())
        rho0 = (rho0 + rho0.conj().T) / 2
        rho0 /= np.trace(rho0).real

        init_results = {}

        # For reference: Peak SumMI at uniform-J
        uniform_peak = None

        for prof_name, J_vec in J_PROFILES.items():
            t0 = _time.time()

            # Mode structure analysis (cheap, N x N diagonalization)
            evals_1mag, amps2 = analyze_mode_structure(J_vec)

            # Build L and propagate
            H = build_H(J_vec)
            L = build_L(H, gamma_uniform)

            sum_mi_traj = []
            for t in t_scan:
                rho_t = evolve(L, rho0, t)
                smi = sum_mi_adjacent(rho_t)
                sum_mi_traj.append(smi)
            sum_mi_traj = np.array(sum_mi_traj)

            peak_smi = float(np.max(sum_mi_traj))
            peak_t = float(t_scan[np.argmax(sum_mi_traj)])
            smi_at_t5 = float(sum_mi_traj[np.argmin(np.abs(t_scan - 5.0))])

            dt_walltime = _time.time() - t0

            if prof_name == "Uniform J=1":
                uniform_peak = peak_smi

            ratio = peak_smi / uniform_peak if uniform_peak and uniform_peak > 1e-10 else float("nan")

            log(f"  {prof_name}:")
            log(f"    J = {J_vec}")
            log(f"    1-magnon energies: [{', '.join(f'{e:+.3f}' for e in evals_1mag)}]")
            log(f"    Mode amplitudes |psi_k(i)|^2 (rows=sites, cols=modes k=1..{N}):")
            for i in range(N):
                row = f"      site {i}: "
                row += "  ".join(f"{amps2[i, k]:.3f}" for k in range(N))
                log(row)
            log(f"    Peak Sum-MI: {peak_smi:.4f} at t={peak_t:.2f}  "
                f"(SumMI@t=5: {smi_at_t5:.4f})")
            log(f"    vs Uniform-J baseline: {ratio:.2f}x")
            log(f"    compute time: {dt_walltime:.1f} s")
            log()

            init_results[prof_name] = {
                "J_vec": J_vec,
                "mode_energies": evals_1mag.tolist(),
                "mode_amplitudes_sq": amps2.tolist(),
                "peak_sum_mi": peak_smi,
                "peak_t": peak_t,
                "sum_mi_at_t5": smi_at_t5,
                "ratio_to_uniform": ratio,
                "sum_mi_trajectory": sum_mi_traj.tolist(),
                "t_scan": t_scan.tolist(),
            }

        all_results[init_name] = init_results

        # ----------------------------------------------------------
        # Summary per initial state
        # ----------------------------------------------------------
        log(f"  Summary for {init_name}:")
        log(f"    {'profile':<40}  {'Peak SumMI':>11}  {'vs Uniform':>12}")
        log("  " + "-" * 72)
        sorted_profiles = sorted(init_results.items(),
                                 key=lambda kv: -kv[1]["peak_sum_mi"])
        for name, res in sorted_profiles:
            log(f"    {name:<40}  {res['peak_sum_mi']:11.4f}  "
                f"{res['ratio_to_uniform']:11.2f}x")
        log()

    # ------------------------------------------------------------------
    # Cross-comparison
    # ------------------------------------------------------------------
    log("=" * 72)
    log("CROSS-COMPARISON")
    log("=" * 72)
    log()
    log("gamma-Sacrifice-Zone reference (non-uniform gamma, uniform J):")
    log("  Peak Sum-MI = 0.230 at N=5 (RESONANT_RETURN Test 8)")
    log("  This requires gamma_edge = 0.246, gamma_other = 0.001.")
    log("  NOT available under gamma_0 = const.")
    log()
    log("J-modulation at gamma_0 = const uniform: best Peak Sum-MI per initial state:")
    for init_name, init_res in all_results.items():
        best = max(init_res.values(), key=lambda r: r["peak_sum_mi"])
        best_name = [k for k, v in init_res.items() if v["peak_sum_mi"] == best["peak_sum_mi"]][0]
        uniform_peak = init_res["Uniform J=1"]["peak_sum_mi"]
        best_peak = best["peak_sum_mi"]
        log(f"  {init_name}:")
        log(f"    uniform-J baseline Peak SumMI = {uniform_peak:.4f}")
        log(f"    best J-profile Peak SumMI    = {best_peak:.4f} ({best_name})")
        log(f"    J-Sacrifice gain over uniform J: {best_peak/uniform_peak:.2f}x "
            f"(vs gamma-Sacrifice 360x)")
    log()

    log(f"Total runtime: {_time.time() - t_start:.1f} s")

    out = {
        "meta": {
            "N": N,
            "gamma_0": 0.05,
            "gamma_profile": "uniform gamma_0 on all sites",
            "J_profiles_tested": {k: v for k, v in J_PROFILES.items()},
            "initial_states": list(INITIAL_STATES.keys()),
            "t_scan_range": [float(t_scan[0]), float(t_scan[-1])],
            "t_scan_n": len(t_scan),
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
