#!/usr/bin/env python3
"""
EQ-024 refinement, Direction 3: J-receiver optimization at N=5 (Heisenberg)
============================================================================
Under gamma_0 = const, the J-channel capacity is the OPERATIONAL capacity
(F30's gamma-side 15.45 is a kinematic reference, see Framing correction in
ClaudeTasks/RESULT_TASK_EQ024_J_MOD_CHANNEL_CAPACITY.md). This script sweeps
receiver states to find the maximum achievable J-capacity at N=5, gamma_0=0.05,
J_ref=1.0 (uniform Heisenberg XX+YY+ZZ).

Constraints:
- Receiver must be F71-symmetric (bond-mirror invariant at N=5, pairs (0,4) and
  (1,3), site 2 self-mirror). Non-F71 probes could unlock larger rank but blow
  up the search space; tackle separately if needed.
- Receiver must NOT be an H-eigenstate (morning theorem: those are J-blind).
- Receiver search: product states parameterized by (theta_a, phi_a) for pair
  (0,4), (theta_b, phi_b) for pair (1,3), (theta_c, phi_c) for site 2. Plus
  explicit non-product F71-symmetric states and morning anchors.

Bloch state: |theta, phi> = cos(theta/2)|0> + exp(i*phi) sin(theta/2)|1>.

Phase 1: structured sweep (~35 receivers, ~18 min).
Phase 2 (conditional): local Nelder-Mead from top 3 if max < 13 bits.

Outputs: simulations/results/eq024_refinement/direction3_*.{txt,json}
"""

import json
import os
import sys
import time as _time

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "eq024_refinement")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_TXT = os.path.join(OUT_DIR, "direction3_receiver_opt.txt")
OUT_JSON = os.path.join(OUT_DIR, "direction3_receiver_opt.json")

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


# Heisenberg bond blocks (full XX+YY+ZZ)
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


def kron_all(vs):
    out = vs[0]
    for v in vs[1:]:
        out = np.kron(out, v)
    return out


def bloch_state(theta, phi):
    return np.array([np.cos(theta / 2),
                     np.exp(1j * phi) * np.sin(theta / 2)], dtype=complex)


def f71_product(theta_a, phi_a, theta_b, phi_b, theta_c, phi_c):
    """Build F71-symmetric product state with (theta, phi) for pair (0,4),
    pair (1,3), and site 2."""
    a = bloch_state(theta_a, phi_a)
    b = bloch_state(theta_b, phi_b)
    c = bloch_state(theta_c, phi_c)
    return kron_all([a, b, c, b, a])


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


def compute_jacobian_J(rho0, J_ref, gammas, t_points, dJ=1e-4):
    n_bonds = len(J_ref)
    y_ref = extract_timeseries(rho0, J_ref, gammas, t_points)
    A = np.zeros((len(y_ref), n_bonds))
    for b in range(n_bonds):
        J_plus = list(J_ref); J_plus[b] += dJ
        y_plus = extract_timeseries(rho0, J_plus, gammas, t_points)
        J_minus = list(J_ref); J_minus[b] -= dJ
        y_minus = extract_timeseries(rho0, J_minus, gammas, t_points)
        A[:, b] = (y_plus - y_minus) / (2.0 * dJ)
    return A


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


# Configuration
J_ref = [1.0] * N_BONDS
gamma_ref = [0.05] * N
t_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
dJ = 1e-4
spread_ref = 0.02
sigma_ref = 0.01
P_ref = N_BONDS * spread_ref ** 2


def capacity_from_psi(psi):
    """Full capacity evaluation for a state vector."""
    rho0 = np.outer(psi, psi.conj())
    rho0 = (rho0 + rho0.conj().T) / 2
    rho0 /= np.trace(rho0).real
    A = compute_jacobian_J(rho0, J_ref, gamma_ref, t_points, dJ=dJ)
    sv = np.linalg.svd(A, compute_uv=False)
    C = waterfilling_capacity(sv, P_ref, sigma_ref)
    return C, sv, A


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("EQ-024 Direction 3: J-receiver optimization at N=5 (Heisenberg)")
    log("=" * 70)
    log(f"Under gamma_0=0.05 const, uniform J=1, search for max J-capacity")
    log(f"over F71-symmetric receivers. Operational Shannon capacity at")
    log(f"spread={spread_ref}, sigma={sigma_ref}, P_total={P_ref:.4e}.")
    log()

    # ----------------------------------------------------------------
    # Phase 1: structured sweep
    # ----------------------------------------------------------------
    log("PHASE 1: structured sweep")
    log("-" * 70)
    log()

    receivers = []  # list of (name, psi)

    # Morning anchors
    pi = np.pi
    receivers.append(("|01010> (morning)",
                      f71_product(0, 0, pi, 0, 0, 0)))
    receivers.append(("|+0+0+> (morning)",
                      f71_product(pi/2, 0, 0, 0, pi/2, 0)))
    receivers.append(("|+-+-+> (morning)",
                      f71_product(pi/2, 0, pi/2, pi, pi/2, 0)))

    # Grid: (theta_a, theta_b, theta_c) in 3-value grid, phi=0 everywhere
    thetas_grid = [pi/4, pi/2, 3*pi/4]
    for ta in thetas_grid:
        for tb in thetas_grid:
            for tc in thetas_grid:
                if ta == pi/2 and tb == pi/2 and tc == pi/2:
                    # |+>^5, H-eigenstate, J-blind; skip
                    continue
                name = f"prod(a={ta:.3f},b={tb:.3f},c={tc:.3f},phi=0)"
                psi = f71_product(ta, 0, tb, 0, tc, 0)
                receivers.append((name, psi))

    # Phase variations at a non-uniform theta: (theta_a, theta_b, theta_c) = (pi/2, 0, pi/2)
    # This is the |+0+0+> morning anchor structure; vary phi_a and phi_c
    for pha, phc in [(pi/2, 0), (0, pi/2), (pi/2, pi/2), (pi, 0), (0, pi)]:
        name = f"phase(a={pi/2:.3f},b={0:.0f},c={pi/2:.3f},phi_a={pha:.2f},phi_c={phc:.2f})"
        psi = f71_product(pi/2, pha, 0, 0, pi/2, phc)
        receivers.append((name, psi))

    # Non-product F71-symmetric superpositions
    def basis(bits):
        idx = int("".join(str(b) for b in bits), 2)
        v = np.zeros(d, dtype=complex)
        v[idx] = 1.0
        return v

    # Antiferromagnetic-like: (|01010> + |10101>)/sqrt(2)
    psi_afm = (basis([0,1,0,1,0]) + basis([1,0,1,0,1])) / np.sqrt(2)
    receivers.append(("(|01010>+|10101>)/sqrt2", psi_afm))

    # Center-only flipped: (|00100> + |11011>)/sqrt(2)
    psi_center = (basis([0,0,1,0,0]) + basis([1,1,0,1,1])) / np.sqrt(2)
    receivers.append(("(|00100>+|11011>)/sqrt2", psi_center))

    # 2-excitation antiferromagnetic: |01010> alone (= morning |01010>)
    # Already included

    # Random F71-symmetric product with complex phases
    rng = np.random.default_rng(42)
    for trial in range(3):
        ta = rng.uniform(0, pi)
        tb = rng.uniform(0, pi)
        tc = rng.uniform(0, pi)
        pha = rng.uniform(0, 2*pi)
        phb = rng.uniform(0, 2*pi)
        phc = rng.uniform(0, 2*pi)
        name = (f"rand#{trial}(a={ta:.2f}/{pha:.2f},"
                f"b={tb:.2f}/{phb:.2f},c={tc:.2f}/{phc:.2f})")
        psi = f71_product(ta, pha, tb, phb, tc, phc)
        receivers.append((name, psi))

    log(f"  Scheduled {len(receivers)} receivers; ~30 s each.")
    log()
    log(f"  {'#':>3}  {'name':<60}  {'|A|_F':>10}  {'sv_max':>10}  {'C (bits)':>9}  {'t (s)':>6}")
    log("  " + "-" * 108)

    results = []
    best_C = 0.0
    best_receiver = None

    for i, (name, psi) in enumerate(receivers):
        t0 = _time.time()
        C, sv, A = capacity_from_psi(psi)
        dt = _time.time() - t0
        frob = np.linalg.norm(A)
        sv_max = float(sv[0]) if len(sv) > 0 else 0.0
        results.append({
            "name": name,
            "frobenius_norm": float(frob),
            "sv_max": sv_max,
            "singular_values": sv.tolist(),
            "capacity_bits": float(C),
            "compute_time_s": float(dt),
        })
        if C > best_C:
            best_C = C
            best_receiver = name
        log(f"  {i+1:>3}  {name:<60}  {frob:10.3e}  {sv_max:10.3e}  {C:9.2f}  {dt:6.1f}")

    log()
    log(f"  Phase 1 best: {best_C:.2f} bits at '{best_receiver}'.")
    log()

    # Sort and show top 10
    results_sorted = sorted(results, key=lambda r: -r["capacity_bits"])
    log("  Top 10 Phase 1 receivers:")
    for j, r in enumerate(results_sorted[:10]):
        log(f"    {j+1:>2}. {r['capacity_bits']:6.2f} bits    {r['name']}")
    log()

    # ----------------------------------------------------------------
    # Phase 2: local optimization (if max < 13 bits)
    # ----------------------------------------------------------------
    phase2_results = []
    if best_C < 13.0:
        log("PHASE 2: local Nelder-Mead from top 3 Phase 1 product-state seeds")
        log("-" * 70)
        log()

        def neg_capacity(params):
            ta, pha, tb, phb, tc, phc = params
            psi = f71_product(ta, pha, tb, phb, tc, phc)
            C, _, _ = capacity_from_psi(psi)
            return -C

        # Pick top 3 Phase 1 results that are product states (have parseable theta/phi)
        # Easier: just pick top 3 from the grid sweep and use their thetas/phis
        # (anchors and superpositions don't parse cleanly; skip them)
        # Better approach: map each Phase-1 result back to its params if known.
        # Simplest: re-seed Nelder-Mead with a few hand-picked starts near top grid.

        # Build starts from top grid points
        grid_tops = [r for r in results_sorted
                     if r["name"].startswith("prod(a=")][:3]
        starts = []
        for r in grid_tops:
            name = r["name"]
            ta = float(name.split("a=")[1].split(",")[0])
            tb = float(name.split("b=")[1].split(",")[0])
            tc = float(name.split("c=")[1].split(",")[0])
            starts.append([ta, 0.0, tb, 0.0, tc, 0.0])

        # Fallback start if grid seeds unavailable
        if not starts:
            starts = [[pi/2, 0, 0, 0, pi/2, 0]]

        # Add one random start
        starts.append([rng.uniform(0, pi), rng.uniform(0, 2*pi),
                       rng.uniform(0, pi), rng.uniform(0, 2*pi),
                       rng.uniform(0, pi), rng.uniform(0, 2*pi)])

        for k, x0 in enumerate(starts):
            log(f"  Start {k+1}: x0 = [{', '.join(f'{x:.3f}' for x in x0)}]")
            log(f"    seed C = {-neg_capacity(x0):.2f} bits")
            t_opt0 = _time.time()
            # Limit iterations to keep budget in check
            result = minimize(neg_capacity, x0, method="Nelder-Mead",
                              options={"maxiter": 40, "xatol": 1e-3, "fatol": 1e-2,
                                       "disp": False})
            dt_opt = _time.time() - t_opt0
            log(f"    optimized C = {-result.fun:.2f} bits "
                f"after {result.nfev} evals ({dt_opt:.1f} s)")
            log(f"    final params = [{', '.join(f'{x:.3f}' for x in result.x)}]")
            phase2_results.append({
                "seed_C": float(-neg_capacity(x0)),
                "x0": list(x0),
                "x_opt": list(result.x),
                "C_opt": float(-result.fun),
                "nfev": int(result.nfev),
                "opt_time_s": float(dt_opt),
            })
            if -result.fun > best_C:
                best_C = -result.fun
                best_receiver = f"Phase2 start {k+1}"
        log()

    log("=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()
    log(f"Phase 1 receivers tested: {len(receivers)}")
    log(f"Phase 1 best C: {results_sorted[0]['capacity_bits']:.2f} bits at "
        f"'{results_sorted[0]['name']}'")
    if phase2_results:
        best_phase2 = max(phase2_results, key=lambda r: r["C_opt"])
        log(f"Phase 2 best C: {best_phase2['C_opt']:.2f} bits")
    log(f"Overall best C: {best_C:.2f} bits ({best_receiver})")
    log(f"F30 gamma-side kinematic reference: 15.45 bits (orphaned operationally)")
    log()
    log(f"Morning anchors: |01010> 11.53, |+0+0+> 10.95, |+-+-+> 11.92")
    log()
    log(f"Total runtime: {_time.time() - t_start:.1f} s")

    out = {
        "meta": {
            "N": N, "N_bonds": N_BONDS, "J_ref": J_ref, "gamma_ref": gamma_ref,
            "t_points": t_points, "dJ": dJ, "fd_method": "central",
            "spread_ref": spread_ref, "sigma_ref": sigma_ref, "P_ref": P_ref,
        },
        "phase1_results": results,
        "phase1_results_sorted": results_sorted,
        "phase2_results": phase2_results,
        "best_C_bits": float(best_C),
        "best_receiver_name": best_receiver,
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fj:
        json.dump(out, fj, indent=2)
    log(f"JSON: {OUT_JSON}")
    log(f"Text: {OUT_TXT}")
    _outf.close()
