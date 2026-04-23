#!/usr/bin/env python3
"""
EQ-024: J-modulation channel capacity at N=5
=============================================
Direct analog of simulations/channel_capacity.py (F30, GAMMA_AS_SIGNAL).
F30 Jacobian: A_gamma = d(observables)/d(gamma_i), shape 150 x 5, SVD -> 15.45 bits at sigma=0.01.
This script: A_J = d(observables)/d(J_b), shape 150 x 4, across several initial states.

Surprise found on first run at |+>^5: A_J = 0 exactly.
Reason: |+>^5 is an H-eigenstate of the Heisenberg (XX+YY+ZZ) chain because
  h_{b,b+1} |++> = |++>  (eigenvalue +1 on every parallel-x pair).
Therefore [H, rho_0] = 0 at every uniform J, and perturbing J_b gives
  [h_b, rho_0] = 0 too (h_b still has |+>^5 as +1-eigenstate).
Duhamel expansion needs [h_b, L^n rho_0]; but L^n rho_0 = L_D^n rho_0 decomposes
into sums of projectors on p-flip-from-|+>^5 subspaces (|psi_{i1..ip}> = prod_k sz^{ik} |+>^5),
each of which is H-invariant (XXX is SU(2)-symmetric, conserves total Sx, hence
"number of flips from |+>^5" is conserved by H). So [H, L_D^n rho_0] = 0 at all n,
and the J-derivative of rho(t) is identically zero.

This makes |+>^5 the J-BLINDEST possible receiver -- inverse of its F30 status
as the OPTIMAL gamma receiver. To get a non-trivial J-capacity number, we need
an initial state that breaks SU(2) (not an H-eigenstate).

This script computes A_J for several initial states:
 - |+>^5 (F30 benchmark, expected 0 by the argument above)
 - |0>^5 (another H-eigenstate: h|00> = |00>, expected 0)
 - GHZ = (|0>^5 + |1>^5)/sqrt(2) (H-eigenstate, expected 0)
 - |01010> (staggered, NOT H-eigenstate, expected J-sensitive)
 - |+0+0+> (z-modulated, breaks SU(2), expected J-sensitive)
 - Dicke |S_1> (lowest-popcount symmetric, NOT generally H-eigenstate for open XXX)

For every state: Jacobian shape 150 x 4, SVD, waterfilling at spread=0.02, sigma=0.01.

Script:  simulations/eq024_j_channel_capacity.py
Outputs: simulations/results/eq024_j_channel_capacity/
"""

import json
import os
import sys
import time as _time

import numpy as np
from scipy.linalg import expm

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "eq024_j_channel_capacity")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_TXT = os.path.join(OUT_DIR, "eq024_j_channel_capacity.txt")
OUT_JSON = os.path.join(OUT_DIR, "eq024_jacobian_n5_multistate.json")

_outf = open(OUT_TXT, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# Pauli + basis
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


# Pre-compute per-bond Heisenberg blocks
_BOND_BLOCKS = []
for b in range(N_BONDS):
    block = (site_op(sx, b) @ site_op(sx, b + 1)
             + site_op(sy, b) @ site_op(sy, b + 1)
             + site_op(sz, b) @ site_op(sz, b + 1))
    _BOND_BLOCKS.append(block)


def build_H_perbond(J_vec):
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


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d_r = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d_r - 1) if d_r > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


# Initial state constructors --------------------------------------------------

def kron_all(vs):
    out = vs[0]
    for v in vs[1:]:
        out = np.kron(out, v)
    return out


def state_plus_N():
    return kron_all([plus] * N)


def state_zero_N():
    return kron_all([up] * N)


def state_ghz():
    psi = (kron_all([up] * N) + kron_all([dn] * N)) / np.sqrt(2)
    return psi


def state_staggered_01010():
    # |01010>
    seq = [up, dn, up, dn, up]  # site 0 = |0>, 1 = |1>, 2 = |0>, 3 = |1>, 4 = |0>
    return kron_all(seq)


def state_plus_zero_plus_zero_plus():
    # |+0+0+>
    seq = [plus, up, plus, up, plus]
    return kron_all(seq)


def state_plus_minus_plus_minus_plus():
    # |+-+-+>  (product of x-eigenstates with alternating sign)
    minus = (up - dn) / np.sqrt(2)
    seq = [plus, minus, plus, minus, plus]
    return kron_all(seq)


def state_dicke_S1():
    # |S_1> = 1/sqrt(N) sum_k sigma_+^k |0>^N  (one-excitation symmetric state)
    psi = np.zeros(d, dtype=complex)
    for k in range(N):
        bits = np.zeros(N, dtype=int)
        bits[k] = 1
        idx = int("".join(str(b) for b in bits), 2)
        psi[idx] += 1.0 / np.sqrt(N)
    return psi


STATE_LIBRARY = {
    "plus^N": state_plus_N,
    "zero^N": state_zero_N,
    "GHZ": state_ghz,
    "|01010>": state_staggered_01010,
    "|+0+0+>": state_plus_zero_plus_zero_plus,
    "|+-+-+>": state_plus_minus_plus_minus_plus,
    "Dicke |S_1>": state_dicke_S1,
}


def verify_H_eigenstate(psi, J_vec):
    """Return (is_eigenstate, eigenvalue_if_yes) for H on psi."""
    H = build_H_perbond(J_vec)
    Hpsi = H @ psi
    # Project Hpsi back onto psi to get the would-be eigenvalue
    norm = np.linalg.norm(psi)
    if norm < 1e-12:
        return False, None
    proj = np.vdot(psi, Hpsi) / np.vdot(psi, psi)
    residual = Hpsi - proj * psi
    is_eig = np.linalg.norm(residual) < 1e-10 * np.linalg.norm(Hpsi + 1e-300)
    return bool(is_eig), complex(proj)


def extract_features(rho):
    f = []
    for q in range(N):
        f.append(purity(ptrace_keep(rho, [q])))
    for i in range(N):
        for j in range(i + 1, N):
            f.append(cpsi(ptrace_keep(rho, [i, j])))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sz @ rho_q))))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sx @ rho_q))))
    return np.array(f)


def extract_timeseries(rho0, J_vec, gammas, t_points):
    H = build_H_perbond(J_vec)
    L = build_L(H, gammas)
    f = []
    for t in t_points:
        rho = evolve(L, rho0, t)
        f.extend(extract_features(rho).tolist())
    return np.array(f)


def compute_jacobian_J(rho0, J_ref, gammas, t_points, dJ=1e-4, method="central"):
    n_bonds = len(J_ref)
    y_ref = extract_timeseries(rho0, J_ref, gammas, t_points)
    n_feat = len(y_ref)
    A = np.zeros((n_feat, n_bonds))

    for b in range(n_bonds):
        J_plus = list(J_ref)
        J_plus[b] += dJ
        y_plus = extract_timeseries(rho0, J_plus, gammas, t_points)
        if method == "forward":
            A[:, b] = (y_plus - y_ref) / dJ
        elif method == "central":
            J_minus = list(J_ref)
            J_minus[b] -= dJ
            y_minus = extract_timeseries(rho0, J_minus, gammas, t_points)
            A[:, b] = (y_plus - y_minus) / (2.0 * dJ)
        else:
            raise ValueError(f"unknown FD method: {method}")
    return A, y_ref


def waterfilling_capacity(singular_values, P_total, sigma_noise):
    n = len(singular_values)
    gains = singular_values ** 2 / sigma_noise ** 2
    idx = np.argsort(gains)[::-1]
    gains_sorted = gains[idx]
    inv_gains = 1.0 / gains_sorted

    k_used = 1
    powers = None
    for k in range(n, 0, -1):
        mu = (P_total + np.sum(inv_gains[:k])) / k
        powers_try = mu - inv_gains[:k]
        if np.all(powers_try >= 0):
            k_used = k
            powers = powers_try
            break

    C = 0.0
    channel_info = []
    for i in range(k_used):
        snr_i = powers[i] * gains_sorted[i]
        c_i = 0.5 * np.log2(1 + snr_i)
        C += c_i
        channel_info.append((float(singular_values[idx[i]]), float(powers[i]),
                             float(snr_i), float(c_i)))
    return float(C), channel_info


def bond_f71_decomposition(v):
    assert len(v) == 4, "bond_f71_decomposition assumes N=5 (4 bonds)"
    e_s1 = np.array([1, 0, 0, 1]) / np.sqrt(2)
    e_s2 = np.array([0, 1, 1, 0]) / np.sqrt(2)
    e_a1 = np.array([1, 0, 0, -1]) / np.sqrt(2)
    e_a2 = np.array([0, 1, -1, 0]) / np.sqrt(2)
    sym_w = np.hypot(np.dot(v, e_s1), np.dot(v, e_s2))
    anti_w = np.hypot(np.dot(v, e_a1), np.dot(v, e_a2))
    total = np.hypot(sym_w, anti_w)
    if total < 1e-10:
        return 0.0, 0.0, "zero"
    sym_frac = sym_w / total
    if sym_frac > 0.99:
        tag = "F71-sym"
    elif sym_frac < 0.01:
        tag = "F71-anti"
    else:
        tag = f"mixed(sym={sym_frac:.3f})"
    return float(sym_w), float(anti_w), tag


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("EQ-024: J-modulation channel capacity at N=5 (multi-state sweep)")
    log("=" * 70)
    log(f"gamma_0 = 0.05 uniform")
    log(f"J_ref = 1.0 uniform (full Heisenberg XX+YY+ZZ)")
    log(f"t_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]")
    log(f"Features per time: 25  ->  total 150")
    log(f"Input: {N_BONDS} per-bond J perturbations")
    log()

    J_ref = [1.0] * N_BONDS
    gamma_ref = [0.05] * N
    t_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    dJ = 1e-4
    spread_ref = 0.02
    sigma_ref = 0.01
    P_ref = N_BONDS * spread_ref ** 2

    F30_SV = [21.3935, 4.5342, 3.2165, 2.8291, 1.4430]
    F30_CAP_0p01 = 15.45

    # -----------------------------------------------------------------
    # H-eigenstate sanity check for each candidate state
    # -----------------------------------------------------------------
    log("STEP 0: H-EIGENSTATE CHECK (at uniform J=1)")
    log("-" * 50)
    for name, ctor in STATE_LIBRARY.items():
        psi = ctor()
        is_eig, eig_val = verify_H_eigenstate(psi, J_ref)
        tag = f"YES  E={eig_val.real:+.4f}" if is_eig else "NO"
        log(f"  {name:>15}  H-eigenstate: {tag}")
    log()
    log("  Prediction: every H-eigenstate has A_J = 0 identically (see module docstring).")
    log()

    # -----------------------------------------------------------------
    # Per-state Jacobian + SVD + capacity
    # -----------------------------------------------------------------
    all_results = {}

    for name, ctor in STATE_LIBRARY.items():
        log("=" * 70)
        log(f"STATE: {name}")
        log("=" * 70)

        psi = ctor()
        rho0_state = np.outer(psi, psi.conj())
        # Ensure rho0 is hermitian + trace-1 (should be by construction)
        rho0_state = (rho0_state + rho0_state.conj().T) / 2
        rho0_state /= np.trace(rho0_state).real

        is_eig, eig_val = verify_H_eigenstate(psi, J_ref)

        t0_state = _time.time()
        A, y_ref = compute_jacobian_J(rho0_state, J_ref, gamma_ref,
                                      t_points, dJ=dJ, method="central")
        dt = _time.time() - t0_state

        U_mat, sv, Vt = np.linalg.svd(A, full_matrices=False)
        frob = np.linalg.norm(A)
        sv_max = float(sv[0]) if len(sv) > 0 else 0.0
        eff_rank = int(np.sum(sv > 0.01 * max(sv_max, 1e-300)))

        log(f"  Jacobian shape: {A.shape}")
        log(f"  Frobenius norm: {frob:.6e}")
        log(f"  Singular values: {', '.join(f'{s:.6e}' for s in sv)}")
        log(f"  Effective rank (sv > 1% of max): {eff_rank}")
        log(f"  Computation time: {dt:.1f} s")

        # If Jacobian is essentially zero, capacity is zero
        if sv_max < 1e-8:
            C_ref = 0.0
            channels = []
            log(f"  Capacity (sigma=0.01, spread=0.02): 0.00 bits  "
                f"[structural zero, max |sv| = {sv_max:.2e}]")
        else:
            C_ref, channels = waterfilling_capacity(sv, P_ref, sigma_ref)
            log(f"  Capacity (sigma=0.01, spread=0.02): {C_ref:.2f} bits  "
                f"(F30 gamma: {F30_CAP_0p01:.2f} bits)")

        # F71 analysis per singular vector (bond space)
        if sv_max > 1e-8:
            log("  SVD bond modes (right singular vectors):")
            for i, s in enumerate(sv):
                v_i = Vt[i]
                sym_w, anti_w, tag = bond_f71_decomposition(v_i)
                log(f"    ch{i+1}  sv={s:.4e}  "
                    f"V=[{', '.join(f'{x:+.3f}' for x in v_i)}]  [{tag}]")

        all_results[name] = {
            "is_H_eigenstate": is_eig,
            "eigenvalue_real": float(eig_val.real) if is_eig else None,
            "jacobian_shape": list(A.shape),
            "frobenius_norm": float(frob),
            "singular_values": sv.tolist(),
            "Vt": Vt.tolist(),
            "effective_rank_1pct": eff_rank,
            "capacity_ref_bits": float(C_ref),
            "channels_ref": [
                {"gain": g, "power": p, "snr": snr, "bits": b}
                for (g, p, snr, b) in channels
            ],
        }
        log()

    # -----------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------
    log("=" * 70)
    log("SUMMARY TABLE  (spread=0.02, sigma=0.01)")
    log("=" * 70)
    log()
    log(f"  {'state':>18}  {'H-eig?':>7}  "
        f"{'|A|_F':>11}  {'sv_max':>11}  {'rank':>5}  {'C (bits)':>9}  {'vs F30 15.45':>14}")
    log("  " + "-" * 92)
    for name, res in all_results.items():
        eig_tag = "yes" if res["is_H_eigenstate"] else "no"
        fro = res["frobenius_norm"]
        sv_max = res["singular_values"][0] if res["singular_values"] else 0.0
        rk = res["effective_rank_1pct"]
        C = res["capacity_ref_bits"]
        delta = C - F30_CAP_0p01
        log(f"  {name:>18}  {eig_tag:>7}  "
            f"{fro:11.4e}  {sv_max:11.4e}  {rk:>5d}  {C:9.2f}  {delta:+14.2f}")
    log()

    log(f"Reference (F30 gamma, |+>^5): gamma-Jacobian 150 x 5,"
        f" SVs = {F30_SV}, C = {F30_CAP_0p01} bits")
    log()
    log(f"Total runtime: {_time.time() - t_start:.1f} s")

    # -----------------------------------------------------------------
    # JSON dump
    # -----------------------------------------------------------------
    out = {
        "meta": {
            "N": N,
            "N_bonds": N_BONDS,
            "J_ref": J_ref,
            "gamma_ref": gamma_ref,
            "t_points": t_points,
            "n_features": 150,
            "dJ": dJ,
            "fd_method": "central",
            "spread_ref": spread_ref,
            "sigma_ref": sigma_ref,
        },
        "states": all_results,
        "f30_reference": {
            "initial_state": "|+>^5",
            "singular_values_gamma": F30_SV,
            "total_capacity_bits_sigma0p01_spread0p02": F30_CAP_0p01,
        },
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fj:
        json.dump(out, fj, indent=2)
    log(f"JSON dump: {OUT_JSON}")
    log(f"Text log:  {OUT_TXT}")

    _outf.close()
