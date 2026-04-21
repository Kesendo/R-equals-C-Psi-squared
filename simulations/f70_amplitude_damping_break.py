#!/usr/bin/env python3
"""F70 operational consequence under amplitude damping.

ORTHOGONALITY_SELECTION_FAMILY.md §3.2(b) (revised 2026-04-21) prediction:
amplitude damping (AD) breaks the *operational* consequence of F70 for the
coherent-only traceless Hermitian probe rho_coh = (|vac><S_2| + |S_2><vac|)/2.
F70 itself (kinematic: single-site partial trace annihilates |Delta N| >= 2
off-diagonal blocks at any fixed t) does not break. Under U(1)-preserving
dynamics rho_coh evolves within the |Delta n| = 2 block and stays site-local
invisible for all t, so c_1_pr(rho_coh) = 0 exactly. AD leaks |Delta n| <= 1
content out of the rho_coh seed, which is site-local visible, so c_1_pr
becomes non-zero at gamma_1 > 0.

Setup at N = 5, uniform XY J = 1, uniform Z-dephasing gamma_0 = 0.05,
site-uniform amplitude damping gamma_1 (swept), open boundary.

**Primary probe:** rho_coh = (|vac><S_2| + |S_2><vac|) / 2 (traceless).
c_1_pr is linear in rho_0 as an observable, so the traceless operator is
admissible as a linear perturbation even though it is not a density matrix.

**Regression check at gamma_1 = 0** (retained from the previous task's
baseline diagnostic): all three probes are measured and the bilinear
decomposition full = population-only + coherent-only is cross-checked.
  - coherent-only  rho_coh = (|vac><S_2| + |S_2><vac|)/2  (primary; expect c_1_pr = 0 exactly)
  - population-only rho_pop = (|vac><vac| + |S_2><S_2|)/2 (expect ~1.352e-3 = K_DD[2,2]_pr/4)
  - full           rho_full = |psi><psi|, |psi>=(|vac>+|S_2>)/sqrt(2)  (expect pop + coh = 1.352e-3)

Observable: c_1_pr at bond (0, 1), t_0 = 1 / gamma_0 = 20,
    c_1_pr = Sum_i dP_B(i, t_0) / dJ  via symmetric difference delta_J = 0.01.

Propagator: scipy.linalg.expm(L * t_0) (per the F73 lesson on eig-path noise
floor with AD-induced imaginary structure).

If coherent-only regression at gamma_1 = 0 does NOT give machine-precision
zero, STOP and document; pipeline state changed since commit d9cf462. The
previous run (commit 54dd994, baseline_diagnostics.json in this results dir)
established that regression at zero.

Rules: UTF-8 stdout, ASCII hyphens only (no em-dashes).
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
# Constants
# ---------------------------------------------------------------------------
N = 5
D = 2 ** N
GAMMA_0 = 0.05
J_UNIFORM = 1.0
T_0 = 1.0 / GAMMA_0   # = 20.0
DELTA_J = 0.01

GAMMA_1_PRIMARY = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1]
GAMMA_1_LOGSWEEP = [0.001, 0.002, 0.005, 0.01, 0.02]

REGRESSION_BAR = 1e-12   # coherent-only c_1_pr at gamma_1 = 0 must be below this
NONZERO_BAR = 1e-6       # coherent-only c_1_pr at gamma_1 > 0 must exceed this

RESULTS_DIR = Path(__file__).parent / "results" / "f70_amplitude_damping_break"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Pauli operators (big-endian: site 0 is MSB in kron ordering)
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)   # sigma^- (lowering)
SP = np.array([[0, 0], [1, 0]], dtype=complex)   # sigma^+ (raising)
NUM = SP @ SM                                     # n = [[0, 0], [0, 1]]


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, n):
    factors = [I2] * n
    factors[site] = op
    return kron_chain(*factors)


# ---------------------------------------------------------------------------
# Hamiltonian and Liouvillian builders
# ---------------------------------------------------------------------------
def build_H_XY(J_list, n):
    """H = Sum_b (J_b / 2) (X_b X_{b+1} + Y_b Y_{b+1}), open boundary."""
    assert len(J_list) == n - 1
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        Jb = float(J_list[b])
        Xb, Xb1 = site_op(X, b, n), site_op(X, b + 1, n)
        Yb, Yb1 = site_op(Y, b, n), site_op(Y, b + 1, n)
        H += (Jb / 2.0) * (Xb @ Xb1 + Yb @ Yb1)
    return H


def build_liouvillian(H, gamma_0, gamma_1, n):
    """L in column-stacking vec convention:
        dvec(rho)/dt = L @ vec(rho)
    with
        L[rho] = -i [H, rho]
                 + gamma_0 * Sum_i (Z_i rho Z_i - rho)
                 + gamma_1 * Sum_i (sigma^-_i rho sigma^+_i
                                    - (1/2) {n_i, rho}).
    """
    d = 2 ** n
    I_d = np.eye(d, dtype=complex)

    # -i [H, rho] -> -i (I kron H - H^T kron I)
    L = -1j * (np.kron(I_d, H) - np.kron(H.T, I_d))

    # Z-dephasing dissipator, uniform gamma_0
    for i in range(n):
        Zi = site_op(Z, i, n)
        L += gamma_0 * (np.kron(Zi.T, Zi) - np.kron(I_d, I_d))

    # Amplitude damping dissipator, uniform gamma_1 (if non-zero)
    if gamma_1 != 0.0:
        for i in range(n):
            SMi = site_op(SM, i, n)
            Ni = site_op(NUM, i, n)
            # vec(sigma^-_i rho sigma^+_i) = ((sigma^+_i)^T kron sigma^-_i) vec(rho)
            # (sigma^+_i)^T = sigma^-_i by construction, so use SMi.T
            L += gamma_1 * np.kron(SMi.T, SMi)
            # -(1/2) {n_i, rho}: vec form = -(1/2)(I kron n_i + n_i^T kron I)
            L -= 0.5 * gamma_1 * (np.kron(I_d, Ni) + np.kron(Ni.T, I_d))
    return L


# ---------------------------------------------------------------------------
# State builders
# ---------------------------------------------------------------------------
def vacuum_ket(n):
    v = np.zeros(2 ** n, dtype=complex)
    v[0] = 1.0
    return v


def dicke_state(n, k):
    """|D(n, k)>: uniform symmetric k-excitation superposition, real
    positive amplitudes. Big-endian basis."""
    v = np.zeros(2 ** n, dtype=complex)
    if k < 0 or k > n:
        return v
    count = comb(n, k)
    norm = 1.0 / np.sqrt(count)
    for bits in range(2 ** n):
        if bin(bits).count("1") == k:
            v[bits] = norm
    return v


def density_matrix(ket):
    return np.outer(ket, ket.conj())


# ---------------------------------------------------------------------------
# Partial trace (big-endian) -> 2x2 reduced state at site i
# ---------------------------------------------------------------------------
def partial_trace_keep_site(rho, i, n):
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


def per_site_purity(rho, n):
    """Return length-n vector of Tr(rho_i^2)."""
    P = np.zeros(n)
    for i in range(n):
        ri = partial_trace_keep_site(rho, i, n)
        P[i] = float(np.trace(ri @ ri).real)
    return P


# ---------------------------------------------------------------------------
# c_1_pr measurement via expm propagation
# ---------------------------------------------------------------------------
def propagate_to(L, rho_0, t):
    """rho(t) = expm(L * t) @ vec(rho_0) reshaped back to (d, d)."""
    d = rho_0.shape[0]
    rho_vec = expm(L * t) @ rho_0.flatten(order="F")
    return rho_vec.reshape(d, d, order="F")


def measure_c1_pr(rho_0, J_list, bond, gamma_1, delta_J, t_0, n, gamma_0=GAMMA_0):
    """c_1_pr(t_0) = Sum_i [P_B(i, t_0, J + dJ) - P_B(i, t_0, J - dJ)] / (2 dJ)
    at bond `bond`. Returns (c_1_pr_sum, per_site_contributions)."""
    J_Bp = list(J_list); J_Bp[bond] = J_list[bond] + delta_J
    J_Bm = list(J_list); J_Bm[bond] = J_list[bond] - delta_J
    H_Bp = build_H_XY(J_Bp, n)
    H_Bm = build_H_XY(J_Bm, n)
    L_Bp = build_liouvillian(H_Bp, gamma_0, gamma_1, n)
    L_Bm = build_liouvillian(H_Bm, gamma_0, gamma_1, n)
    rho_Bp = propagate_to(L_Bp, rho_0, t_0)
    rho_Bm = propagate_to(L_Bm, rho_0, t_0)
    P_Bp = per_site_purity(rho_Bp, n)
    P_Bm = per_site_purity(rho_Bm, n)
    dP_B = (P_Bp - P_Bm) / (2.0 * delta_J)
    return float(np.sum(dP_B)), dP_B.tolist()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def build_probes(n):
    """Return the three relevant probes as density matrices."""
    vac = vacuum_ket(n)
    S2 = dicke_state(n, 2)
    psi = (vac + S2) / np.sqrt(2.0)
    rho_full = density_matrix(psi)                         # task spec
    rho_coh = 0.5 * (np.outer(vac, S2.conj())
                     + np.outer(S2, vac.conj()))           # diagnostic (traceless)
    rho_pop = 0.5 * (density_matrix(vac) + density_matrix(S2))  # diagnostic
    # Cross-check: rho_full == rho_pop + rho_coh
    assert np.allclose(rho_full, rho_pop + rho_coh), "probe decomposition mismatch"
    return {"full": rho_full, "coherent_only": rho_coh, "population_only": rho_pop}


def run_baseline_and_diagnostics(n, J_list, bond, delta_J, t_0):
    """Run c_1_pr at gamma_1 = 0 on all three probes. Returns a dict."""
    probes = build_probes(n)
    out = {}
    for name, rho in probes.items():
        c1, per_site = measure_c1_pr(rho, J_list, bond, 0.0, delta_J, t_0, n)
        out[name] = {
            "c_1_pr": c1,
            "abs_c_1_pr": abs(c1),
            "per_site_contributions": per_site,
        }
    # Linearity sanity check (bilinear kernel is additive across DD/CC sectors
    # under U(1), so c_1(full) should equal c_1(pop) + c_1(coh) up to float noise)
    out["linearity_residual"] = (out["full"]["c_1_pr"]
                                 - out["population_only"]["c_1_pr"]
                                 - out["coherent_only"]["c_1_pr"])
    return out


def run_sweep(n, J_list, bond, delta_J, t_0, gamma_1_list, rho_0,
              baseline_c1=0.0):
    """Sweep over gamma_1 values. For rows with gamma_1 > 0, pass the nonzero
    bar if |c_1_pr - baseline_c1| > NONZERO_BAR; the subtraction is trivial
    (zero) for the coherent-only probe whose baseline is already 0."""
    rows = []
    for g1 in gamma_1_list:
        c1, per_site = measure_c1_pr(rho_0, J_list, bond, g1, delta_J, t_0, n)
        signal = c1 - baseline_c1
        rows.append({
            "gamma_1": g1,
            "c_1_pr": c1,
            "c_1_pr_minus_baseline": signal,
            "abs_signal": abs(signal),
            "per_site_contributions": per_site,
            "pass_nonzero_bar": (g1 == 0.0) or (abs(signal) > NONZERO_BAR),
        })
    return rows


def fit_log_log(rows):
    """Fit log(|signal|) vs log(gamma_1). Returns slope, intercept, R^2."""
    g = np.array([r["gamma_1"] for r in rows])
    c = np.array([r["abs_signal"] for r in rows])
    mask = (g > 0) & (c > 0)
    log_g = np.log(g[mask])
    log_c = np.log(c[mask])
    slope, intercept = np.polyfit(log_g, log_c, 1)
    pred = slope * log_g + intercept
    ss_res = float(np.sum((log_c - pred) ** 2))
    ss_tot = float(np.sum((log_c - log_c.mean()) ** 2))
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), float(r_sq)


def main():
    t_start = time.time()

    print("=" * 78)
    print(f"F70 AD break: c_1_pr under amplitude damping at N = {N}")
    print(f"Primary probe: rho_coh = (|vac><S_2| + |S_2><vac|)/2  (traceless Hermitian)")
    print(f"Bond: (0, 1), gamma_0 = {GAMMA_0}, t_0 = {T_0}, delta_J = {DELTA_J}")
    print(f"Regression bar: |c_1_pr(coherent-only, gamma_1=0)| < {REGRESSION_BAR:.0e}")
    print(f"Nonzero bar:    |signal(gamma_1>0)| > {NONZERO_BAR:.0e}")
    print(f"Propagator: scipy.linalg.expm(L * t_0)")
    print("=" * 78)

    J_uniform = [J_UNIFORM] * (N - 1)

    # ----- Regression check: all three probes at gamma_1 = 0 -----
    print("\n--- Regression check at gamma_1 = 0 (three-probe bilinear cross-check) ---")
    diag = run_baseline_and_diagnostics(N, J_uniform, 0, DELTA_J, T_0)

    c1_full = diag["full"]["c_1_pr"]
    c1_coh = diag["coherent_only"]["c_1_pr"]
    c1_pop = diag["population_only"]["c_1_pr"]
    lin_res = diag["linearity_residual"]
    regression_pass = abs(c1_coh) < REGRESSION_BAR

    print(f"  coherent-only (primary):  c_1_pr = {c1_coh:+.6e}  |c_1_pr| = {abs(c1_coh):.3e}")
    print(f"  population-only (cross):  c_1_pr = {c1_pop:+.6e}  |c_1_pr| = {abs(c1_pop):.3e}")
    print(f"  full |psi><psi| (cross):  c_1_pr = {c1_full:+.6e}  |c_1_pr| = {abs(c1_full):.3e}")
    print(f"  linearity residual (full - pop - coh) = {lin_res:+.3e}")
    print(f"  regression pass (< {REGRESSION_BAR:.0e}): {regression_pass}")

    regression_payload = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "delta_J": DELTA_J,
        "t_0": T_0, "bond": [0, 1],
        "primary_probe": "coherent_only",
        "diagnostics": diag,
        "regression_bar": REGRESSION_BAR,
        "regression_pass": regression_pass,
    }

    if not regression_pass:
        print()
        print("!" * 78)
        print("REGRESSION FAIL: coherent-only baseline at gamma_1 = 0 not machine zero.")
        print("Pipeline state may have changed since commit d9cf462. STOP.")
        print("!" * 78)
        regression_path = RESULTS_DIR / "regression_fail.json"
        with open(regression_path, "w", encoding="utf-8") as f:
            json.dump(regression_payload, f, indent=2)
        print(f"\nsaved {regression_path}")
        print(f"\nTotal walltime: {time.time() - t_start:.1f} s")
        return 2

    probes = build_probes(N)
    rho_coh = probes["coherent_only"]
    # baseline for coherent-only is exactly zero; no subtraction needed
    baseline_c1 = c1_coh

    # ----- Primary sweep -----
    print(f"\n--- Primary sweep: coherent-only probe, bond (0, 1), gamma_1 in {GAMMA_1_PRIMARY} ---")
    primary = run_sweep(N, J_uniform, 0, DELTA_J, T_0, GAMMA_1_PRIMARY, rho_coh,
                        baseline_c1=baseline_c1)
    for r in primary:
        print(f"  gamma_1 = {r['gamma_1']:>6.4f}: c_1_pr = {r['c_1_pr']:+.6e}  "
              f"|signal| = {r['abs_signal']:.3e}  nonzero = {r['pass_nonzero_bar']}")

    # ----- Stretch 1: low-gamma_1 log-log fit -----
    print(f"\n--- Stretch 1: low-gamma_1 log-log scaling ---")
    low = run_sweep(N, J_uniform, 0, DELTA_J, T_0, GAMMA_1_LOGSWEEP, rho_coh,
                    baseline_c1=baseline_c1)
    for r in low:
        print(f"  gamma_1 = {r['gamma_1']:>6.4f}: c_1_pr = {r['c_1_pr']:+.6e}  |signal| = {r['abs_signal']:.3e}")
    slope, intercept, r_sq = fit_log_log(low)
    print(f"  log-log slope = {slope:.4f}  intercept = {intercept:.4f}  R^2 = {r_sq:.6f}")

    # ----- Stretch 2: bond (2, 3) spot check -----
    print(f"\n--- Stretch 2: bond (2, 3), gamma_1 = 0.05 ---")
    c1_23, ps_23 = measure_c1_pr(rho_coh, J_uniform, 2, 0.05, DELTA_J, T_0, N)
    c1_01 = next(r["c_1_pr"] for r in primary if r["gamma_1"] == 0.05)
    print(f"  c_1_pr(bond (2, 3), gamma_1 = 0.05) = {c1_23:+.6e}")
    print(f"  c_1_pr(bond (0, 1), gamma_1 = 0.05) = {c1_01:+.6e}  (reference)")
    ratio = c1_23 / c1_01 if c1_01 != 0.0 else float("nan")
    print(f"  ratio (2,3) / (0,1) = {ratio:+.4f}")

    # ----- Save JSONs -----
    with open(RESULTS_DIR / "sweep_coherent.json", "w", encoding="utf-8") as f:
        json.dump({
            "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "delta_J": DELTA_J,
            "t_0": T_0, "bond": [0, 1],
            "probe": "coherent_only",
            "probe_definition": "rho_coh = (|vac><S_2| + |S_2><vac|) / 2",
            "regression_at_gamma_1_zero": regression_payload,
            "gamma_1_values": GAMMA_1_PRIMARY,
            "sweep": primary,
        }, f, indent=2)
    with open(RESULTS_DIR / "low_gamma_scaling.json", "w", encoding="utf-8") as f:
        json.dump({
            "N": N, "gamma_0": GAMMA_0, "bond": [0, 1],
            "probe": "coherent_only",
            "gamma_1_values": GAMMA_1_LOGSWEEP,
            "sweep_points": low,
            "log_log_slope": slope,
            "log_log_intercept": intercept,
            "R_squared": r_sq,
        }, f, indent=2)
    with open(RESULTS_DIR / "bond_2_3_spot.json", "w", encoding="utf-8") as f:
        json.dump({
            "N": N, "gamma_0": GAMMA_0,
            "probe": "coherent_only",
            "bond": [2, 3], "gamma_1": 0.05,
            "c_1_pr": c1_23, "per_site_contributions": ps_23,
            "bond_0_1_reference_at_same_gamma_1": c1_01,
            "ratio_2_3_over_0_1": ratio,
        }, f, indent=2)
    print(f"\nsaved sweep_coherent.json, low_gamma_scaling.json, bond_2_3_spot.json")

    # ----- Optional overlay plot -----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        g1_p = [r["gamma_1"] for r in primary]
        c1_p = [r["c_1_pr"] for r in primary]
        s_p = [r["abs_signal"] for r in primary]
        ax1.plot(g1_p, c1_p, "o-", linewidth=1.4, markersize=8, label="c_1_pr(gamma_1)")
        ax1.axhline(0.0, color="gray", linewidth=0.7)
        ax1.set_xlabel("gamma_1  (amplitude damping rate, site-uniform)")
        ax1.set_ylabel("c_1_pr (coherent-only probe, bond 0,1)")
        ax1.set_title(f"Primary sweep, N = {N}, t_0 = {T_0}")
        ax1.legend(fontsize=9, loc="best")
        ax1.grid(True, alpha=0.3)

        g1_s = np.array([r["gamma_1"] for r in low if r["gamma_1"] > 0])
        s_s = np.array([r["abs_signal"] for r in low if r["gamma_1"] > 0])
        ax2.loglog(g1_s, s_s, "o", markersize=8, label="data")
        fit_line = np.exp(intercept) * g1_s ** slope
        ax2.loglog(g1_s, fit_line, "r--", linewidth=1,
                   label=f"fit: slope = {slope:.3f}, R^2 = {r_sq:.4f}")
        ax2.set_xlabel("gamma_1")
        ax2.set_ylabel("|c_1_pr|  (coherent-only probe)")
        ax2.set_title("Stretch 1: low-gamma_1 log-log scaling")
        ax2.legend(fontsize=9, loc="best")
        ax2.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        plot_path = RESULTS_DIR / "f70_amplitude_damping_break_overlay.png"
        fig.savefig(plot_path, dpi=120)
        plt.close(fig)
        print(f"saved plot {plot_path}")
    except ImportError:
        print("matplotlib not available, skipping plot")

    all_nonzero = all(r["pass_nonzero_bar"] for r in primary if r["gamma_1"] > 0)
    print(f"\nAll gamma_1 > 0 pass nonzero bar: {all_nonzero}")
    print(f"Total walltime: {time.time() - t_start:.1f} s")
    return 0 if all_nonzero else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="F70 operational consequence under amplitude damping at N = 5.")
    args = parser.parse_args()
    sys.exit(main())
