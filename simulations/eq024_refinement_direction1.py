#!/usr/bin/env python3
"""
EQ-024 refinement, Direction 1: H-scope of M_x-polynomial blindness
====================================================================
Morning session (`eq024_j_channel_capacity.py`) proved that for full
Heisenberg H = Sum (XX+YY+ZZ), every H-eigenstate has zero J-Jacobian
at all orders. Proof Step 5 invokes [H, M_x] = 0 via SU(2) symmetry,
which fails for XY-only H = Sum (XX+YY).

This script tests the H-scope: swap XX+YY+ZZ with XX+YY, rerun the four
morning-Heisenberg-eigenstate receivers, report J-Jacobian per state.

Analytical predictions for XY-only:
- |0>^N  :  H_XY|0>^N = 0 (eigenvalue 0, because h_XY|00> = |11> - |11> = 0).
            sigma_z|0> = |0>, so |0>^N is in the DFS of L_D:
            L_D * rho_0 = 0 literally. rho(t) = rho_0 exact, J-blind strictly.
- GHZ    :  (|0>^N + |1>^N)/sqrt(2). Both |0>^N and |1>^N are H_XY-eigenstates
            with eigenvalue 0. GHZ is H_XY-eigenstate with eigenvalue 0.
            L_D acts on GHZ coherence as -2*gamma*N rate (exponential decay in
            the 2-dim {|GHZ><GHZ|, |GHZ_-><GHZ_-|} subspace). Both endpoints
            are H_XY-eigenoperators, so [H_XY, rho(t)] = 0 for all t. J-blind.
- |+>^N  :  h_XY|++> = XX|++> + YY|++> = |++> - |-->.
            NOT H_XY-eigenstate. Should be J-sensitive under XY-only.
- Dicke |S_1> : 1-magnon sector under XY is pure tight-binding (no diagonal,
                since h_XY|01> = 2|10> and h_XY|00> = 0). Eigenstates are
                sine waves sin(k*pi*j/(N+1)), not the uniform |S_1>. So
                |S_1> is NOT H_XY-eigenstate. Should be J-sensitive.

Outputs: simulations/results/eq024_refinement/
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
OUT_TXT = os.path.join(OUT_DIR, "direction1_xy_only.txt")
OUT_JSON = os.path.join(OUT_DIR, "direction1_xy_only.json")

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


# Build bond blocks for both H choices; selected at build time by the flag.
_BOND_BLOCKS_HEIS = []
_BOND_BLOCKS_XY = []
for b in range(N_BONDS):
    xx = site_op(sx, b) @ site_op(sx, b + 1)
    yy = site_op(sy, b) @ site_op(sy, b + 1)
    zz = site_op(sz, b) @ site_op(sz, b + 1)
    _BOND_BLOCKS_HEIS.append(xx + yy + zz)
    _BOND_BLOCKS_XY.append(xx + yy)


def build_H_perbond(J_vec, H_type="XY"):
    blocks = _BOND_BLOCKS_XY if H_type == "XY" else _BOND_BLOCKS_HEIS
    H = np.zeros((d, d), dtype=complex)
    for b in range(N_BONDS):
        H += J_vec[b] * blocks[b]
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


def kron_all(vs):
    out = vs[0]
    for v in vs[1:]:
        out = np.kron(out, v)
    return out


STATE_LIBRARY = {
    "plus^N":       lambda: kron_all([plus] * N),
    "zero^N":       lambda: kron_all([up] * N),
    "GHZ":          lambda: (kron_all([up] * N) + kron_all([dn] * N)) / np.sqrt(2),
    "Dicke |S_1>":  lambda: _build_dicke_S1(),
}


def _build_dicke_S1():
    psi = np.zeros(d, dtype=complex)
    for k in range(N):
        bits = np.zeros(N, dtype=int)
        bits[k] = 1
        idx = int("".join(str(b) for b in bits), 2)
        psi[idx] += 1.0 / np.sqrt(N)
    return psi


def verify_H_eigenstate(psi, J_vec, H_type):
    H = build_H_perbond(J_vec, H_type=H_type)
    Hpsi = H @ psi
    proj = np.vdot(psi, Hpsi) / np.vdot(psi, psi)
    residual = Hpsi - proj * psi
    res_norm = np.linalg.norm(residual)
    total = np.linalg.norm(Hpsi) + 1e-300
    is_eig = res_norm < 1e-10 * total
    return bool(is_eig), complex(proj), float(res_norm / total)


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


def extract_timeseries(rho0, J_vec, gammas, t_points, H_type):
    H = build_H_perbond(J_vec, H_type=H_type)
    L = build_L(H, gammas)
    f = []
    for t in t_points:
        rho = evolve(L, rho0, t)
        f.extend(extract_features(rho).tolist())
    return np.array(f)


def compute_jacobian_J(rho0, J_ref, gammas, t_points, H_type, dJ=1e-4):
    n_bonds = len(J_ref)
    y_ref = extract_timeseries(rho0, J_ref, gammas, t_points, H_type)
    n_feat = len(y_ref)
    A = np.zeros((n_feat, n_bonds))
    for b in range(n_bonds):
        J_plus = list(J_ref)
        J_plus[b] += dJ
        y_plus = extract_timeseries(rho0, J_plus, gammas, t_points, H_type)
        J_minus = list(J_ref)
        J_minus[b] -= dJ
        y_minus = extract_timeseries(rho0, J_minus, gammas, t_points, H_type)
        A[:, b] = (y_plus - y_minus) / (2.0 * dJ)
    return A, y_ref


def waterfilling_capacity(sv, P_total, sigma_noise):
    n = len(sv)
    if np.max(sv) < 1e-8:
        return 0.0
    gains = sv ** 2 / sigma_noise ** 2
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
    for i in range(k_used):
        snr_i = powers[i] * gains_sorted[i]
        C += 0.5 * np.log2(1 + snr_i)
    return float(C)


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("EQ-024 Direction 1: H-scope of M_x-polynomial blindness")
    log("=" * 70)
    log(f"Swap H = XX+YY+ZZ (Heisenberg, morning) vs H = XX+YY (XY-only).")
    log(f"Four receivers that were Heisenberg-eigenstates in the morning.")
    log()
    log(f"N = {N}, gamma_0 = 0.05 uniform, J = 1.0 uniform, "
        f"t_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]")
    log(f"Features per time: 25  ->  total 150")
    log(f"dJ = 1e-4 (central differences)")
    log(f"Waterfilling ref: spread=0.02, sigma=0.01")
    log()

    J_ref = [1.0] * N_BONDS
    gamma_ref = [0.05] * N
    t_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    dJ = 1e-4
    spread_ref = 0.02
    sigma_ref = 0.01
    P_ref = N_BONDS * spread_ref ** 2

    all_results = {}

    for H_type, H_tag in [("Heisenberg", "XX+YY+ZZ"), ("XY", "XX+YY")]:
        log("=" * 70)
        log(f"H = {H_tag}  (H_type='{H_type}')")
        log("=" * 70)

        # Per-state H-eigenstate check
        log()
        log(f"  {'state':>14}  {'H-eig?':>7}  {'E':>10}  {'residual':>12}")
        log("  " + "-" * 50)
        for name, ctor in STATE_LIBRARY.items():
            psi = ctor()
            is_eig, eig_val, residual = verify_H_eigenstate(psi, J_ref, H_type)
            tag = f"yes" if is_eig else "no"
            log(f"  {name:>14}  {tag:>7}  {eig_val.real:+10.4f}  {residual:12.4e}")
        log()

        # Per-state Jacobian + SVD + capacity
        log(f"  {'state':>14}  {'|A|_F':>12}  {'sv_max':>12}  {'rank_1pct':>10}  "
            f"{'C (bits)':>9}  {'verdict':>14}")
        log("  " + "-" * 88)
        h_results = {}
        for name, ctor in STATE_LIBRARY.items():
            psi = ctor()
            rho0_state = np.outer(psi, psi.conj())
            rho0_state = (rho0_state + rho0_state.conj().T) / 2
            rho0_state /= np.trace(rho0_state).real

            t0 = _time.time()
            A, _ = compute_jacobian_J(rho0_state, J_ref, gamma_ref,
                                      t_points, H_type=H_type, dJ=dJ)
            dt = _time.time() - t0
            sv = np.linalg.svd(A, compute_uv=False)
            frob = np.linalg.norm(A)
            sv_max = float(sv[0]) if len(sv) > 0 else 0.0
            rank = int(np.sum(sv > 0.01 * max(sv_max, 1e-300)))
            C = waterfilling_capacity(sv, P_ref, sigma_ref)

            if sv_max < 1e-8:
                verdict = "J-blind"
            elif sv_max > 1.0:
                verdict = "J-sensitive"
            else:
                verdict = "weak"

            log(f"  {name:>14}  {frob:12.4e}  {sv_max:12.4e}  {rank:>10d}  "
                f"{C:9.2f}  {verdict:>14}")

            h_results[name] = {
                "frobenius_norm": float(frob),
                "sv_max": float(sv_max),
                "singular_values": sv.tolist(),
                "rank_1pct": rank,
                "capacity_bits": float(C),
                "verdict": verdict,
                "compute_time_s": float(dt),
            }
        all_results[H_type] = h_results
        log()

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log("=" * 70)
    log("SUMMARY: morning predictions vs measured")
    log("=" * 70)
    log()
    log("Predictions for XY-only:")
    log("  |0>^5       : J-blind (DFS of L_D, H-independent)")
    log("  GHZ         : J-blind (2-dim H-degenerate subspace)")
    log("  |+>^5       : J-sensitive (not an H_XY-eigenstate)")
    log("  Dicke |S_1> : J-sensitive (not an H_XY-eigenstate)")
    log()
    log(f"  {'state':>14}  {'Heisenberg':>14}  {'XY-only':>14}  {'match?':>12}")
    log("  " + "-" * 62)
    xy_predictions = {
        "plus^N": "J-sensitive",
        "zero^N": "J-blind",
        "GHZ": "J-blind",
        "Dicke |S_1>": "J-sensitive",
    }
    for name in STATE_LIBRARY:
        h_v = all_results["Heisenberg"][name]["verdict"]
        xy_v = all_results["XY"][name]["verdict"]
        match = "yes" if xy_v == xy_predictions[name] else "NO"
        log(f"  {name:>14}  {h_v:>14}  {xy_v:>14}  {match:>12}")
    log()

    log(f"Total runtime: {_time.time() - t_start:.1f} s")

    out = {
        "meta": {
            "N": N,
            "N_bonds": N_BONDS,
            "J_ref": J_ref,
            "gamma_ref": gamma_ref,
            "t_points": t_points,
            "dJ": dJ,
            "fd_method": "central",
            "spread_ref": spread_ref,
            "sigma_ref": sigma_ref,
        },
        "results": all_results,
        "xy_predictions": xy_predictions,
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fj:
        json.dump(out, fj, indent=2)
    log(f"JSON: {OUT_JSON}")
    log(f"Text: {OUT_TXT}")

    _outf.close()
