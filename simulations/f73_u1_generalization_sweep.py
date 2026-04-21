#!/usr/bin/env python3
"""F73 U(1) generalization sweep at N = 5.

F73 (docs/ANALYTICAL_FORMULAS.md line ~1570) closes the spatial-sum
coherence purity for the (vac, S_1) probe under uniform XY with uniform
Z-dephasing:

    Sum_i 2 * |(rho_coh,i)_{0,1}(t)|^2 = (1/2) * exp(-4 * gamma_0 * t)

The structural derivation only uses (i) [H, N_total] = 0 so dynamics
stay in the single-excitation (SE) sector, (ii) H Hermitian so U_SE
unitary, (iii) uniform Z-dephasing gamma_0. It does not invoke XY
structure. This script is the numerical sanity check that F73 extends
to any U(1)-preserving Hermitian H.

Six residual checks at N = 5, gamma_0 = 0.05, t in [0, 40], 81 points:
  - XXZ family, Delta in {0, 0.5, 1.0, 2.0}.
  - Stretch 1: random Haar SE superposition at Delta = 1.0.
  - Stretch 2: inhomogeneous XY (Delta = 0, J_b uniform in [0.5, 1.5]).

Probe: rho_0 = |psi><psi| with |psi> = (|vac> + |S_1>) / sqrt(2),
|S_1> = Sum_i |S_1,i>/sqrt(N) (except Stretch 1, which replaces |S_1>
by a random Haar SE superposition).

Pass bar: max_t |S(t) - S_pred(t)| < 1e-14 for every setup.

Rules: UTF-8 stdout, ASCII hyphens only (no em-dashes).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig, expm

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N = 5
D = 2 ** N
GAMMA_0 = 0.05
J_UNIFORM = 1.0
T_MAX = 40.0
N_TIMES = 81
RESIDUAL_BAR = 1e-14
SURPRISE_BAR = 1e-12
DEFAULT_SEED = 20260421

RESULTS_DIR = Path(__file__).parent / "results" / "f73_u1_generalization"
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


# ---------------------------------------------------------------------------
# State builders
# ---------------------------------------------------------------------------
def vacuum_ket(n):
    v = np.zeros(2 ** n, dtype=complex)
    v[0] = 1.0
    return v


def uniform_S1_ket(n):
    """|S_1> = Sum_i |S_1,i> / sqrt(n), where |S_1,i> has qubit i excited
    (big-endian)."""
    v = np.zeros(2 ** n, dtype=complex)
    norm = 1.0 / np.sqrt(n)
    for bits in range(2 ** n):
        if bin(bits).count("1") == 1:
            v[bits] = norm
    return v


def se_site_kets(n):
    """Return kets [|S_1,0>, ..., |S_1,n-1>] where |S_1,i> has qubit i
    excited, in big-endian convention."""
    kets = []
    for i in range(n):
        v = np.zeros(2 ** n, dtype=complex)
        v[1 << (n - 1 - i)] = 1.0
        kets.append(v)
    return kets


def custom_SE_ket(coeffs, n):
    """|psi_SE> = Sum_i coeffs[i] |S_1,i>."""
    assert len(coeffs) == n
    kets = se_site_kets(n)
    v = np.zeros(2 ** n, dtype=complex)
    for i, c in enumerate(coeffs):
        v += c * kets[i]
    return v


def density_matrix(ket):
    return np.outer(ket, ket.conj())


# ---------------------------------------------------------------------------
# Hamiltonian builders
# ---------------------------------------------------------------------------
def build_H_XXZ(J_list, delta, n):
    """H = Sum_b J_b * [ (X_b X_{b+1} + Y_b Y_{b+1}) / 2 + delta * Z_b Z_{b+1} ]
    with open boundary conditions. J_list has length n - 1.
    """
    assert len(J_list) == n - 1
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        Jb = float(J_list[b])
        Xb, Xb1 = site_op(X, b, n), site_op(X, b + 1, n)
        Yb, Yb1 = site_op(Y, b, n), site_op(Y, b + 1, n)
        Zb, Zb1 = site_op(Z, b, n), site_op(Z, b + 1, n)
        H += Jb * (0.5 * (Xb @ Xb1 + Yb @ Yb1) + delta * (Zb @ Zb1))
    # Enforce Hermiticity up to round-off.
    H = 0.5 * (H + H.conj().T)
    return H


# ---------------------------------------------------------------------------
# Liouvillian and propagation (column-stacking vec convention)
# ---------------------------------------------------------------------------
def build_liouvillian(H, gamma_0, n):
    """Vectorised L such that dvec(rho)/dt = L @ vec(rho), column-major.
    Lindblad: -i [H, rho] + gamma_0 * Sum_i (Z_i rho Z_i - rho)."""
    d = 2 ** n
    I_d = np.eye(d, dtype=complex)
    L = -1j * (np.kron(I_d, H) - np.kron(H.T, I_d))
    for i in range(n):
        Zi = site_op(Z, i, n)
        L += gamma_0 * (np.kron(Zi.T, Zi) - np.kron(I_d, I_d))
    return L


def propagate(L, rho_0, times):
    """Propagate rho(t) = expm(L t) vec(rho_0) via step-wise application of
    U_dt = expm(L * dt). Uses expm (scaling-and-squaring Pade) rather than
    eigendecomposition: for Liouvillians with large spectral extent (XXZ at
    Delta >= 2 pushes |Im(lambda)| into the tens), np.linalg.eig's right
    eigenvector matrix becomes poorly conditioned and the eigendecomposition
    propagation floor rises to ~1e-14. expm hits ~1e-15 or better uniformly
    across all Delta tested here. Returns an (n_times, d, d) complex array.

    For eigendecomposition-propagation comparison, see propagate_eig below.
    """
    d = rho_0.shape[0]
    if len(times) < 2:
        raise ValueError("Need at least 2 time points for propagate.")
    dt = float(times[1] - times[0])
    assert np.allclose(np.diff(times), dt), "times must be uniform for step propagator"
    U_dt = expm(L * dt)
    out = np.empty((len(times), d, d), dtype=complex)
    rho_vec = rho_0.flatten(order="F").astype(complex)
    for k in range(len(times)):
        if k > 0:
            rho_vec = U_dt @ rho_vec
        out[k] = rho_vec.reshape(d, d, order="F")
    return out


def propagate_eig(L, rho_0, times):
    """Alternative propagator via eigendecomposition of L. Retained as an
    independent numerical cross-check. Gives identical results at machine
    precision for small Liouvillian spectral extent; noise floor rises to
    ~1e-14 for XXZ Delta >= 2 at N = 5."""
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    d = rho_0.shape[0]
    rho0_vec = rho_0.flatten(order="F")
    c0 = V_Linv @ rho0_vec
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order="F")
    return out


# ---------------------------------------------------------------------------
# Partial trace: 2x2 reduced density matrix at site i (big-endian)
# ---------------------------------------------------------------------------
def partial_trace_keep_site(rho, i, n):
    """Trace out all qubits except site i. rho is (2^n, 2^n) big-endian."""
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


# ---------------------------------------------------------------------------
# F73 closure measurement
# ---------------------------------------------------------------------------
def sum_coh_purity(rho, n):
    """S = Sum_i 2 * |(rho_coh,i)_{0,1}|^2."""
    acc = 0.0
    for i in range(n):
        ri = partial_trace_keep_site(rho, i, n)
        acc += 2.0 * abs(ri[0, 1]) ** 2
    return acc


# ---------------------------------------------------------------------------
# Per-setup runner
# ---------------------------------------------------------------------------
def run_setup(label, H, rho_0, times, gamma_0, config_extra=None):
    L = build_liouvillian(H, gamma_0, N)
    rho_traj = propagate(L, rho_0, times)
    S_t = np.array([sum_coh_purity(rho_traj[k], N) for k in range(len(times))])
    S_pred = 0.5 * np.exp(-4.0 * gamma_0 * times)
    residual = np.abs(S_t - S_pred)
    max_res = float(np.max(residual))
    argmax_idx = int(np.argmax(residual))
    pass_flag = bool(max_res < RESIDUAL_BAR)
    surprise_flag = bool(max_res > SURPRISE_BAR)

    config = {
        "N": N,
        "gamma_0": gamma_0,
        "t_max": T_MAX,
        "n_times": N_TIMES,
        "residual_bar": RESIDUAL_BAR,
        "surprise_bar": SURPRISE_BAR,
    }
    if config_extra:
        config.update(config_extra)

    return {
        "label": label,
        "config": config,
        "times": times.tolist(),
        "S_t": S_t.tolist(),
        "S_pred": S_pred.tolist(),
        "residual": residual.tolist(),
        "max_residual": max_res,
        "max_residual_t": float(times[argmax_idx]),
        "pass_flag_1e14": pass_flag,
        "surprise_flag_1e12": surprise_flag,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(seed: int | None) -> int:
    t_start = time.time()
    seed_effective = int(seed) if seed is not None else DEFAULT_SEED

    print("=" * 76)
    print(f"F73 U(1) generalization sweep at N = {N}, gamma_0 = {GAMMA_0}")
    print(f"Probe: (|vac> + |S_1>)/sqrt(2)  (Stretch 1 replaces |S_1>)")
    print(f"t in [0, {T_MAX}], {N_TIMES} points; dt = {T_MAX / (N_TIMES - 1)}")
    print(f"Residual bar: {RESIDUAL_BAR:.0e}   Surprise bar: {SURPRISE_BAR:.0e}")
    print(f"Seed (random SE + inhomogeneous J): {seed_effective}")
    print("=" * 76)

    times = np.linspace(0.0, T_MAX, N_TIMES)
    vac = vacuum_ket(N)
    S1 = uniform_S1_ket(N)
    psi_default = (vac + S1) / np.sqrt(2.0)
    rho_default = density_matrix(psi_default)
    J_uniform = [J_UNIFORM] * (N - 1)

    results = {}

    # --- Primary: XXZ family ---
    for delta in (0.0, 0.5, 1.0, 2.0):
        label = f"xxz_delta_{delta:g}".replace(".", "p")
        H = build_H_XXZ(J_uniform, delta, N)
        extra = {
            "family": "xxz",
            "delta": float(delta),
            "J_list": J_uniform,
            "probe": "uniform_S1_plus_vac",
        }
        print(f"\n[setup] {label}: XXZ delta = {delta}, uniform J, uniform probe")
        r = run_setup(label, H, rho_default, times, GAMMA_0, extra)
        print(f"  max residual = {r['max_residual']:.3e} at t = {r['max_residual_t']:.2f}  "
              f"pass(1e-14) = {r['pass_flag_1e14']}"
              f"{'  [SURPRISE]' if r['surprise_flag_1e12'] else ''}")
        results[label] = r

    # --- Stretch 1: random Haar SE superposition at Delta = 1.0 ---
    rng = np.random.default_rng(seed_effective)
    c = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    c /= np.linalg.norm(c)
    psi_rand_SE = custom_SE_ket(c.tolist(), N)
    psi_rand = (vac + psi_rand_SE) / np.sqrt(2.0)
    rho_rand = density_matrix(psi_rand)
    label = "xxz_delta_1p0_random_SE"
    H_heis = build_H_XXZ(J_uniform, 1.0, N)
    extra = {
        "family": "xxz",
        "delta": 1.0,
        "J_list": J_uniform,
        "probe": "random_Haar_SE_plus_vac",
        "SE_coefficients_real": c.real.tolist(),
        "SE_coefficients_imag": c.imag.tolist(),
        "seed": seed_effective,
    }
    print(f"\n[setup] {label}: Heisenberg delta = 1.0, random Haar SE probe, seed = {seed_effective}")
    r = run_setup(label, H_heis, rho_rand, times, GAMMA_0, extra)
    print(f"  max residual = {r['max_residual']:.3e} at t = {r['max_residual_t']:.2f}  "
          f"pass(1e-14) = {r['pass_flag_1e14']}"
          f"{'  [SURPRISE]' if r['surprise_flag_1e12'] else ''}")
    results[label] = r

    # --- Stretch 2: inhomogeneous XY (Delta = 0, J_b ~ U[0.5, 1.5]) ---
    rng2 = np.random.default_rng(seed_effective + 1)
    J_inhom = rng2.uniform(0.5, 1.5, size=N - 1).tolist()
    label = "xy_delta_0_inhomogeneous_J"
    H_inh = build_H_XXZ(J_inhom, 0.0, N)
    extra = {
        "family": "xy_inhomogeneous",
        "delta": 0.0,
        "J_list": J_inhom,
        "probe": "uniform_S1_plus_vac",
        "seed": seed_effective + 1,
    }
    J_str = ", ".join(f"{j:.4f}" for j in J_inhom)
    print(f"\n[setup] {label}: XY delta = 0, J = [{J_str}], seed = {seed_effective + 1}")
    r = run_setup(label, H_inh, rho_default, times, GAMMA_0, extra)
    print(f"  max residual = {r['max_residual']:.3e} at t = {r['max_residual_t']:.2f}  "
          f"pass(1e-14) = {r['pass_flag_1e14']}"
          f"{'  [SURPRISE]' if r['surprise_flag_1e12'] else ''}")
    results[label] = r

    # --- Save per-setup JSONs ---
    for label, r in results.items():
        path = RESULTS_DIR / f"{label}.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(r, f, indent=2)
        print(f"saved {path}")

    # --- Summary table ---
    print("\n" + "=" * 76)
    print("SUMMARY")
    print("=" * 76)
    header = f"{'setup':<38} {'max residual':>14} {'pass 1e-14':>11}  {'flag':>9}"
    print(header)
    print("-" * len(header))
    for label, r in results.items():
        flag = "SURPRISE" if r["surprise_flag_1e12"] else ""
        print(f"{label:<38} {r['max_residual']:>14.3e} "
              f"{str(r['pass_flag_1e14']):>11}  {flag:>9}")

    all_pass = all(r["pass_flag_1e14"] for r in results.values())
    any_surprise = any(r["surprise_flag_1e12"] for r in results.values())
    print("-" * len(header))
    print(f"All setups pass 1e-14: {all_pass}")
    print(f"Any surprise (> 1e-12): {any_surprise}")

    # --- Optional overlay plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
        for label, r in results.items():
            ax1.plot(r["times"], r["S_t"], label=label, alpha=0.75, linewidth=1.2)
        ax1.plot(times, 0.5 * np.exp(-4.0 * GAMMA_0 * times),
                 "k--", linewidth=1.5, label="S_pred = (1/2) exp(-4 gamma_0 t)")
        ax1.set_yscale("log")
        ax1.set_ylabel("S(t)")
        ax1.set_title(f"F73 U(1) generalization at N = {N}, gamma_0 = {GAMMA_0}")
        ax1.legend(fontsize=8, loc="lower left")
        ax1.grid(True, which="both", alpha=0.3)

        for label, r in results.items():
            ax2.semilogy(r["times"], np.maximum(r["residual"], 1e-20),
                         label=label, alpha=0.75, linewidth=1.2)
        ax2.axhline(RESIDUAL_BAR, color="red", linestyle=":",
                    label=f"pass bar {RESIDUAL_BAR:.0e}")
        ax2.axhline(SURPRISE_BAR, color="orange", linestyle=":",
                    label=f"surprise bar {SURPRISE_BAR:.0e}")
        ax2.set_xlabel("t")
        ax2.set_ylabel("|S(t) - S_pred(t)|")
        ax2.legend(fontsize=8, loc="lower right")
        ax2.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        plot_path = RESULTS_DIR / "f73_u1_generalization_overlay.png"
        fig.savefig(plot_path, dpi=120)
        plt.close(fig)
        print(f"saved plot {plot_path}")
    except ImportError:
        print("matplotlib not available, skipping plot")

    print(f"\nTotal walltime: {time.time() - t_start:.1f} s")
    return 0 if all_pass else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="F73 U(1) generalization sweep (N = 5, 6 setups).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Seed for the random Haar SE vector and the "
                             "inhomogeneous J realization. Default: deterministic.")
    args = parser.parse_args()
    sys.exit(main(seed=args.seed))
