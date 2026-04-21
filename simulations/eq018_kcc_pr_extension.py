#!/usr/bin/env python3
"""eq018_kcc_pr_extension.py

K_CC[n, n+1]_pr extension to n >= 1 sector blocks at N = 5.

F73 closes the spatial-sum coherence purity functional
    S(t) = Sum_i 2 * |(rho_coh,i)_{0,1}(t)|^2 = (1/2) * exp(-4 gamma_0 t)
for the (vac, S_1) probe under any U(1)-preserving Hermitian H with
uniform gamma_0. Consequence: K_CC[0,1]_pr = dS/d(delta J_b) = 0 exactly.

This script tests whether the closure extends to interior sector blocks
(n, n+1) with n >= 1. Structural sketch (to test, not trust):

- rho_mix = (|S_n><S_n| + |S_{n+1}><S_{n+1}|)/2 stays block-diagonal in
  popcount under L = H-commutator + Z-dephasing, so S[rho_mix(t)] = 0
  for all t.
- Consequence: M_CC_pr(J, t) := S[rho_coh(t)] - S[rho_mix(t)]
                              = S[rho_cc(t)]
  where rho_cc_init = (|S_n><S_{n+1}| + h.c.) / 2 is the coherence-only
  perturbation (traceless, Hermitian, not a density matrix). L is linear,
  so propagation is well-defined.
- If amplitude stays in the 1-site-differing slice under H: S(t) decays
  cleanly as S(0) * exp(-4 gamma_0 t), K_CC_pr = 0.
- If H leaks amplitude into the 3-site-differing slice (invisible under
  single-site partial trace): S(t) decays faster than exp(-4 gamma_0 t),
  and the decay depends on J, giving K_CC_pr != 0.

For the vac-SE block (n=0), H cannot produce 3-site-differing at any
order (SE block is closed at popcount 1, and the only 1-site-differing
coherences from vac are (|vac><w_i| + h.c.)). That is F73's structural
privilege.

For (n=1, n=2) at N=5: by explicit calculation, H at first order keeps
amplitude in the 1-site-differing slice. At second order H couples to
3-site-differing: e.g. H^2 |w_0><S_{01}| produces |w_2><S_{01}| (3-site
different). At gamma_0 t = 1 we are non-perturbative in H t, so the
compound effect is fully developed.

Prediction: S(t) for (n >= 1, n+1) probe deviates from S(0) exp(-4 g t)
at late time, and K_CC[n, n+1]_pr(t) != 0 in general.

Refutation, a second structural privilege, or an orthogonal structure
are all valid outcomes.

Closed-form S(0) for any (n, n+1) probe, any N:
    S(0) = (N - n) * (n + 1) / (2 * N)
(derived in prompt; numerical sanity check included below.)

Rules: em-dashes forbidden; hyphens only. UTF-8 stdout.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Constants (aligned with f73_u1_generalization_sweep.py)
# ---------------------------------------------------------------------------
N_DEFAULT = 5
GAMMA_0 = 0.05
J_UNIFORM = 1.0
DJ = 0.01
T_MAX = 40.0
N_TIMES = 81
T_POINTWISE = 20.0  # anchor from RESULT_TASK_EQ018_C1_POINTWISE

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_kcc_pr_extension"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Pauli operators and chain helpers (big-endian: site 0 is MSB)
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, n):
    factors = [I2] * n
    factors[site] = op
    return kron_chain(*factors)


def build_H_XY(J_list, n):
    """H = Sum_b (J_b / 2) (X_b X_{b+1} + Y_b Y_{b+1}) with OBC."""
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        Jb = float(J_list[b])
        H += (Jb / 2.0) * (
            site_op(X, b, n) @ site_op(X, b + 1, n)
            + site_op(Y, b, n) @ site_op(Y, b + 1, n)
        )
    return H


def build_liouvillian(H, gamma_0, n):
    """L such that dvec(rho)/dt = L @ vec(rho), column-stacking convention."""
    d = 2 ** n
    I_d = np.eye(d, dtype=complex)
    L = -1j * (np.kron(I_d, H) - np.kron(H.T, I_d))
    for i in range(n):
        Zi = site_op(Z, i, n)
        L += gamma_0 * (np.kron(Zi.T, Zi) - np.kron(I_d, I_d))
    return L


# ---------------------------------------------------------------------------
# Dicke states (big-endian)
# ---------------------------------------------------------------------------
def dicke_state(n, weight):
    """|S_w> = (1/sqrt(C(n,w))) Sum_{popcount(x) = w} |x>, big-endian."""
    d = 2 ** n
    v = np.zeros(d, dtype=complex)
    norm = 1.0 / np.sqrt(comb(n, weight))
    for bits in range(d):
        if bin(bits).count("1") == weight:
            v[bits] = norm
    return v


# ---------------------------------------------------------------------------
# Propagator (step-wise expm(L dt))
# ---------------------------------------------------------------------------
def propagate(L, rho0, times):
    """Propagate rho(t) for each t in (uniform) times via U_dt = expm(L dt).
    rho0 is a Hermitian operator (not necessarily positive).
    """
    d = rho0.shape[0]
    if len(times) < 2:
        raise ValueError("Need at least 2 time points.")
    dt = float(times[1] - times[0])
    assert np.allclose(np.diff(times), dt), "times must be uniform"
    U_dt = expm(L * dt)
    out = np.empty((len(times), d, d), dtype=complex)
    rho_vec = rho0.flatten(order="F").astype(complex)
    for k in range(len(times)):
        if k > 0:
            rho_vec = U_dt @ rho_vec
        out[k] = rho_vec.reshape(d, d, order="F")
    return out


# ---------------------------------------------------------------------------
# Partial trace (keep one site)
# ---------------------------------------------------------------------------
def partial_trace_keep_site(rho, i, n):
    """Tr_{~i}(rho), returning a 2x2 matrix. Big-endian."""
    shape_2N = [2] * (2 * n)
    out = rho.reshape(shape_2N)
    ket_axes = list(range(n))
    bra_axes = list(range(n, 2 * n))
    for j in range(n - 1, -1, -1):
        if j == i:
            continue
        a_k = ket_axes[j]
        a_b = bra_axes[j]
        out = np.trace(out, axis1=a_k, axis2=a_b)
        lo, hi = sorted((a_k, a_b))
        for k in range(n):
            if k == j:
                continue
            if ket_axes[k] > hi:
                ket_axes[k] -= 2
            elif ket_axes[k] > lo:
                ket_axes[k] -= 1
            if bra_axes[k] > hi:
                bra_axes[k] -= 2
            elif bra_axes[k] > lo:
                bra_axes[k] -= 1
    a_k = ket_axes[i]
    a_b = bra_axes[i]
    if a_k == 1 and a_b == 0:
        out = out.T
    return out


def spatial_sum_coh_purity(rho, n):
    """S = Sum_i 2 * |(rho_i)_{0,1}|^2."""
    acc = 0.0
    for i in range(n):
        ri = partial_trace_keep_site(rho, i, n)
        acc += 2.0 * abs(ri[0, 1]) ** 2
    return acc


# ---------------------------------------------------------------------------
# Probe builders
# ---------------------------------------------------------------------------
def rho_cc_dicke(n, n_qubits):
    """rho_cc = (|S_n><S_{n+1}| + |S_{n+1}><S_n|) / 2 at N = n_qubits."""
    S_n = dicke_state(n_qubits, n)
    S_np1 = dicke_state(n_qubits, n + 1)
    block = 0.5 * (np.outer(S_n, S_np1.conj()) + np.outer(S_np1, S_n.conj()))
    return block


def rho_cc_site_local(site, spectator_bits, n_qubits):
    """Site-local spectator probe:
        rho_cc_site = (|0_site, s><1_site, s| + h.c.) / 2
    where s is a computational-basis spectator pattern of fixed popcount on
    the other N-1 sites. Pure 1-site-differing-at-'site' content, no
    3-site-differing contamination.

    spectator_bits is an (N-1)-long list of 0/1 values for sites
    {0..N-1}\{site}, in ascending site-index order.
    """
    d = 2 ** n_qubits
    assert len(spectator_bits) == n_qubits - 1
    # Build ket_0 and ket_1 computational basis states
    idx_0 = 0
    idx_1 = 0
    spec_iter = iter(spectator_bits)
    for s in range(n_qubits):
        if s == site:
            bit_0 = 0
            bit_1 = 1
        else:
            spec_bit = next(spec_iter)
            bit_0 = spec_bit
            bit_1 = spec_bit
        idx_0 |= bit_0 << (n_qubits - 1 - s)
        idx_1 |= bit_1 << (n_qubits - 1 - s)
    ket_0 = np.zeros(d, dtype=complex); ket_0[idx_0] = 1.0
    ket_1 = np.zeros(d, dtype=complex); ket_1[idx_1] = 1.0
    return 0.5 * (np.outer(ket_0, ket_1.conj()) + np.outer(ket_1, ket_0.conj()))


# ---------------------------------------------------------------------------
# K_CC_pr extraction
# ---------------------------------------------------------------------------
def extract_S_trajectory(rho_cc_init, L, times, n_qubits):
    """Propagate rho_cc_init through L and return S(t_k) for each t_k."""
    rho_traj = propagate(L, rho_cc_init, times)
    S = np.array([spatial_sum_coh_purity(rho_traj[k], n_qubits) for k in range(len(times))])
    return S


def compute_KCC_pr(rho_cc_init, L_A, L_Bp, L_Bm, times, n_qubits, dJ):
    """K_CC_pr(t) = [S(J+dJ, t) - S(J-dJ, t)] / (2 dJ), time series."""
    S_A = extract_S_trajectory(rho_cc_init, L_A, times, n_qubits)
    S_Bp = extract_S_trajectory(rho_cc_init, L_Bp, times, n_qubits)
    S_Bm = extract_S_trajectory(rho_cc_init, L_Bm, times, n_qubits)
    K = (S_Bp - S_Bm) / (2.0 * dJ)
    return {"S_A": S_A, "S_Bp": S_Bp, "S_Bm": S_Bm, "K": K}


# ---------------------------------------------------------------------------
# Analytical reference: S_ref(t) = S(0) * exp(-4 gamma_0 t)
# ---------------------------------------------------------------------------
def S_zero_closed_form(n, n_qubits):
    """S(0) = (N - n) * (n + 1) / (2 N) for (n, n+1) Dicke probe."""
    return (n_qubits - n) * (n + 1) / (2.0 * n_qubits)


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------
def run(n_qubits, bonds, n_values, gamma_0, J, dJ, t_max, n_times):
    times = np.linspace(0.0, t_max, n_times)
    out_bondwise = {}

    print("=" * 78)
    print(f"K_CC[n, n+1]_pr extension to n >= 1")
    print(f"N = {n_qubits}, gamma_0 = {gamma_0}, J = {J}, dJ = {dJ}")
    print(f"t in [0, {t_max}], {n_times} points, dt = {t_max / (n_times - 1):.4f}")
    print(f"Bonds: {bonds}")
    print(f"n values (coherence block (n, n+1)): {n_values}")
    print("=" * 78)

    # --- Analytical reference S(0) sanity check ---
    print("\n[sanity] Closed-form S(0) vs numerical for each (n, n+1) probe:")
    for n_coh in n_values:
        rho_cc_init = rho_cc_dicke(n_coh, n_qubits)
        S0_num = spatial_sum_coh_purity(rho_cc_init, n_qubits)
        S0_cf = S_zero_closed_form(n_coh, n_qubits)
        print(f"  n={n_coh}, (n,n+1)=({n_coh},{n_coh+1}): "
              f"S(0) numerical = {S0_num:.6f}  closed-form = {S0_cf:.6f}  "
              f"residual = {abs(S0_num - S0_cf):.2e}")

    # --- Build L_A once (shared across bonds) ---
    t0 = time.time()
    J_A = [J] * (n_qubits - 1)
    H_A = build_H_XY(J_A, n_qubits)
    L_A = build_liouvillian(H_A, gamma_0, n_qubits)
    print(f"\n[build] L_A built in {time.time() - t0:.1f} s")

    for bond in bonds:
        print(f"\n{'-' * 78}")
        print(f"Bond b = {bond}  (perturbation on sites {bond}, {bond + 1})")
        print(f"{'-' * 78}")
        t0 = time.time()
        J_Bp = list(J_A); J_Bp[bond] = J + dJ
        J_Bm = list(J_A); J_Bm[bond] = J - dJ
        L_Bp = build_liouvillian(build_H_XY(J_Bp, n_qubits), gamma_0, n_qubits)
        L_Bm = build_liouvillian(build_H_XY(J_Bm, n_qubits), gamma_0, n_qubits)
        print(f"  L_B+/- built in {time.time() - t0:.1f} s")

        bond_results = {}
        for n_coh in n_values:
            rho_cc_init = rho_cc_dicke(n_coh, n_qubits)
            res = compute_KCC_pr(rho_cc_init, L_A, L_Bp, L_Bm, times,
                                  n_qubits, dJ)
            S_A = res["S_A"]
            K = res["K"]
            # S_ref(t) = S(0) * exp(-4 gamma_0 t)
            S0 = S_zero_closed_form(n_coh, n_qubits)
            S_ref = S0 * np.exp(-4.0 * gamma_0 * times)
            S_deviation = S_A - S_ref

            # Pointwise values at t = T_POINTWISE
            idx_pointwise = int(np.argmin(np.abs(times - T_POINTWISE)))
            t_pw = times[idx_pointwise]

            # Peak (in absolute value) of K across t
            idx_peak = int(np.argmax(np.abs(K)))
            t_peak = times[idx_peak]

            bond_results[f"n_{n_coh}"] = {
                "n_coh": n_coh,
                "S_0_closed_form": float(S0),
                "S_A_pointwise_t20": float(S_A[idx_pointwise]),
                "S_ref_pointwise_t20": float(S_ref[idx_pointwise]),
                "S_deviation_pointwise_t20": float(S_deviation[idx_pointwise]),
                "K_pointwise_t20": float(K[idx_pointwise]),
                "K_peak_abs": float(K[idx_peak]),
                "K_peak_t": float(t_peak),
                "times": times.tolist(),
                "S_A": S_A.tolist(),
                "S_Bp": res["S_Bp"].tolist(),
                "S_Bm": res["S_Bm"].tolist(),
                "S_ref": S_ref.tolist(),
                "S_deviation": S_deviation.tolist(),
                "K": K.tolist(),
            }

            print(f"\n  (n, n+1) = ({n_coh}, {n_coh+1}):")
            print(f"    S(0) closed-form = {S0:.6f}")
            print(f"    S(t=20) numerical = {S_A[idx_pointwise]:.6e}")
            print(f"    S_ref(t=20) = S(0) exp(-4 gamma_0 t=20) = {S_ref[idx_pointwise]:.6e}")
            print(f"    S deviation at t=20: {S_deviation[idx_pointwise]:+.6e}")
            print(f"    K_CC[{n_coh},{n_coh+1}]_pr at t=20: {K[idx_pointwise]:+.6e}")
            print(f"    Peak |K| across t: {abs(K[idx_peak]):.6e} at t = {t_peak:.2f}")

        out_bondwise[f"bond_{bond}"] = bond_results

    return {
        "config": {
            "N": n_qubits,
            "gamma_0": gamma_0,
            "J": J,
            "dJ": dJ,
            "t_max": t_max,
            "n_times": n_times,
            "T_pointwise": T_POINTWISE,
            "bonds": list(bonds),
            "n_values": list(n_values),
        },
        "bonds": out_bondwise,
    }


def main():
    parser = argparse.ArgumentParser(
        description="K_CC[n, n+1]_pr extension for n >= 1 at configurable N.")
    parser.add_argument("--N", type=int, default=N_DEFAULT,
                        help=f"Number of qubits (default {N_DEFAULT}).")
    parser.add_argument("--bonds", type=str, default=None,
                        help="Comma-separated bond indices. Default: all N-1 bonds.")
    parser.add_argument("--n-values", type=str, default=None,
                        help="Comma-separated n values. Default: 0..N-1 (n+1 must be <= N).")
    parser.add_argument("--gamma-0", type=float, default=GAMMA_0)
    parser.add_argument("--J", type=float, default=J_UNIFORM)
    parser.add_argument("--dJ", type=float, default=DJ)
    parser.add_argument("--t-max", type=float, default=T_MAX)
    parser.add_argument("--n-times", type=int, default=N_TIMES)
    parser.add_argument("--tag", type=str, default="default",
                        help="Tag for output filename: kcc_pr_{tag}.json")
    args = parser.parse_args()

    N = args.N
    if args.bonds is None:
        bonds = list(range(N - 1))
    else:
        bonds = [int(x) for x in args.bonds.split(",")]
    if args.n_values is None:
        n_values = list(range(N))  # 0..N-1 so (n, n+1) makes sense
    else:
        n_values = [int(x) for x in args.n_values.split(",")]
    # Sanity: n+1 must be <= N
    n_values = [n for n in n_values if n + 1 <= N]

    t_start = time.time()
    result = run(
        n_qubits=N,
        bonds=bonds,
        n_values=n_values,
        gamma_0=args.gamma_0,
        J=args.J,
        dJ=args.dJ,
        t_max=args.t_max,
        n_times=args.n_times,
    )

    out_path = RESULTS_DIR / f"kcc_pr_{args.tag}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved: {out_path}")
    print(f"Total walltime: {time.time() - t_start:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
