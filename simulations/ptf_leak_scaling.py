#!/usr/bin/env python3
"""Edge 2 of the PTF fresh-eyes chain (2026-06-10): leak-channel scaling probe.

QUESTION. The Π-break experiment (experiments/PTF_PALINDROME_BREAKING_PERTURBATIONS.md)
found the Z-field arm keeps the PTF closure Σ_i ln α_i ≈ 0, but with a small honest
residual that wanders O(0.5), peaking near ε ≈ 1 before receding. The leak channel
got a closed form the same day (simulations/f87_deg1_face_cell_free.py): the palindrome
asymmetry enters the odd power-sums first at m = 3 with p₃(γ) = 6·4^N·ε²·γ. Does the
small-ε closure residual scale like ε²·γ (the leak channel's leading coefficient)?

RESULT (the answer is NO, a clean refutation, banked in the experiment doc's 2026-06-10
update): the residual is FIRST order in ε with a nearly γ-independent coefficient,
S ≈ k·ε with k ≈ −0.58…−0.61 across γ ∈ [0.025, 0.1] (ε-exponent 0.88-0.99 at
R² ≥ 0.997 down to ε = 0.005; the ε²γ collapse table varies ~50× with a sign change).
The wander is not the spectral-palindrome break leaking in; it is EQ-014's first-order
non-closure (review/EQ014_FINDINGS.md), seen here for a Z-field perturbation instead of
a J-defect, with a site-dependent first-order profile whose imperfect cancellation IS
the residual. γ-blindness is the fingerprint of the Hamiltonian-side mechanism.

PROTOCOL (reused verbatim from simulations/ptf_transverse_field_pi_break.py, the
committed experiment):
  - uniform open XY chain H = Σ_b (J/2)(X_b X_{b+1} + Y_b Y_{b+1}), J = 1.0;
  - Z-dephasing via the Hamming mask, dρ −= 2γ·h⊙ρ; RK4 with dt = 0.05, T_max = 20,
    401-point uniform time grid;
  - bonding-mode initial state φ = (|vac⟩ + |ψ₁⟩)/√2 (fw.bonding_mode_pair_state(N, 1));
  - per-site purity P_i(t) = ½(1 + ⟨X⟩² + ⟨Y⟩² + ⟨Z⟩²), α-fit P_A(α·t) ≈ P_B(t) with
    cubic interpolation, bounded minimize_scalar on [0.1, 10], xatol = 1e-7;
    S = Σ_i ln α_i (ptf.fit_all_sites, imported, not re-implemented);
  - Z-field on SITE 0 (the committed Phase-4/5 site, the site of the hard-drive scan
    that produced the O(0.5) wandering note).

IMPORTANT PROTOCOL FACT (determined from the committed script): the Z-field arm has NO
J-defect. The closure comparison is baseline (uniform XY, no field) vs the SAME chain
+ ε·Z₀; the Z-field is present only in run B. The J-defect (δJ on bond (0,1)) is the
separate CONTROL arm of the experiment. Hence S(0, γ) = 0 identically (A-vs-A null
control) and the leak ΔS(ε, γ) = S(ε, γ) − S(0, γ) = S(ε, γ) up to fit noise.
A secondary panel below ALSO runs the other reading (Z-field present in BOTH runs of
the δJ = 0.1 J-defect closure) so both designs are on the table.

Panels:
  1  faithfulness check: local γ-parameterized propagator == ptf.purity_trajectory
     at γ = 0.05 (bit-near-exact), and the A-vs-A null control.
  2  MAIN GRID: N = 5, ε ∈ {0, 0.02, 0.05, 0.1, 0.2, 0.4} × γ ∈ {0.025, 0.05, 0.1},
     Z on site 0. Table (S, ΔS, ΔS/(ε²γ), sign, noise flags), exponent regressions
     ln|ΔS| ~ a·ln ε + b (fixed γ; hypothesis a ≈ 2) and ln|ΔS| ~ c·ln γ + d (fixed ε;
     hypothesis c ≈ 1), collapse table ΔS/(ε²γ).
     Noise handling: (i) integration floor = |S| of baseline(dt=0.05) vs
     baseline(dt=0.025) fitted on the same grid; (ii) fit-window systematic = spread
     of S across two window variants (t ≤ 15 truncation; every-2nd-point subsample);
     (iii) all regressions repeated per window variant - a stable exponent across
     variants is the honest signal, not any single S value.
  3  N = 6 (only if N = 5 ran in seconds): same grid.
  4  robustness: Z on the CENTER site (N = 5, γ = 0.05) - same ε exponent?
  5  other design: Z-field in BOTH runs of the δJ = 0.1 J-defect closure (N = 5,
     γ = 0.05): ΔS_B(ε) = S_J(ε) − S_J(0).
  6  small-ε refinement (N = 5, site 0): ε down to 0.005; dS/ε converges to a
     constant k(γ) ⟹ the leak's leading order is ε¹, with k nearly γ-independent;
     Richardson residual (dS/ε − k)/ε probes any subleading ε² term and its γ-trend;
     γ = 0 degeneracy note (the α-ansatz needs decaying trajectories); per-site
     ln α_i / ε profile showing the first-order response site-by-site.

Run: python simulations/_ptf_leak_scaling.py
WIP probe - not committed, no files written.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ptf_transverse_field_pi_break as ptf  # noqa: E402  the committed machinery
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

J_UNIFORM = ptf.J_UNIFORM            # 1.0
T_MAX = ptf.T_MAX                    # 20.0
DT = ptf.DT                          # 0.05
N_STEPS = ptf.N_STEPS                # 400
T_GRID = np.linspace(0.0, T_MAX, N_STEPS + 1)

EPS_GRID = [0.0, 0.02, 0.05, 0.10, 0.20, 0.40]
GAMMA_GRID = [0.025, 0.05, 0.10]
DEFECT_BOND = ptf.DEFECT_BOND        # (0, 1) - control arm only
DELTA_J = 0.10                       # the experiment's headline control δJ

VARIANTS = ("std", "t<=15", "sub2")  # fit-window variants


# ---------------------------------------------------------------------------
# Propagation: ptf.purity_trajectory with γ and dt as explicit parameters
# ---------------------------------------------------------------------------
def purity_traj(N, H, rho_0, bloch_ops, hamming, gamma, dt=DT, n_steps=N_STEPS):
    """Replicates ptf.purity_trajectory line-for-line, with γ/dt as arguments
    instead of module globals. Verified against the committed function below."""
    rho = rho_0.copy()
    P = np.zeros((n_steps + 1, N))

    def rhs(r):
        return -1j * (H @ r - r @ H) - 2.0 * gamma * hamming * r

    def record(r, row):
        for i, (Xi, Yi, Zi) in enumerate(bloch_ops):
            x = np.real(np.trace(Xi @ r))
            y = np.real(np.trace(Yi @ r))
            z = np.real(np.trace(Zi @ r))
            P[row, i] = 0.5 * (1.0 + x * x + y * y + z * z)

    record(rho, 0)
    for step in range(1, n_steps + 1):
        k1 = rhs(rho)
        k2 = rhs(rho + 0.5 * dt * k1)
        k3 = rhs(rho + 0.5 * dt * k2)
        k4 = rhs(rho + dt * k3)
        rho = rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        rho = 0.5 * (rho + rho.conj().T)
        record(rho, step)
    return P


def fit_S_variants(P_A, P_B):
    """S = Σ ln α under the standard window and two fit-window variants.

    Returns dict variant -> S, plus max per-site RMSE of the standard fit."""
    _, r_std, S_std = ptf.fit_all_sites(T_GRID, P_A, P_B)
    _, _, S_w15 = ptf.fit_all_sites(T_GRID[:301], P_A[:301], P_B[:301])
    _, _, S_sub = ptf.fit_all_sites(T_GRID[::2], P_A[::2], P_B[::2])
    return {"std": S_std, "t<=15": S_w15, "sub2": S_sub}, float(r_std.max())


def loglog_fit(xs, ys):
    """Least-squares fit ln|y| = slope·ln x + icept. Returns slope, icept, R²."""
    xs = np.log(np.asarray(xs, dtype=float))
    ys = np.log(np.abs(np.asarray(ys, dtype=float)))
    A = np.vstack([xs, np.ones_like(xs)]).T
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((ys - pred) ** 2))
    ss_tot = float(np.sum((ys - ys.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(coef[0]), float(coef[1]), r2


# ---------------------------------------------------------------------------
# Panel 1: faithfulness + null control
# ---------------------------------------------------------------------------
def panel1():
    print("=" * 78)
    print("PANEL 1  Faithfulness: local propagator == committed ptf.purity_trajectory")
    print("=" * 78)
    N = 5
    bloch = ptf.site_bloch_ops(N)
    ham = ptf.hamming_matrix(N)
    psi0 = fw.bonding_mode_pair_state(N, 1)
    rho0 = np.outer(psi0, psi0.conj())
    H_A = ptf.build_H_xy(N, [J_UNIFORM] * (N - 1))

    P_ref = ptf.purity_trajectory(N, H_A, rho0, bloch, ham)   # GAMMA_0 = 0.05 baked in
    P_new = purity_traj(N, H_A, rho0, bloch, ham, gamma=0.05)
    err = float(np.max(np.abs(P_ref - P_new)))
    assert err < 1e-12, f"propagator mismatch: {err:.2e}"
    print(f"  max|P_committed − P_local| at γ=0.05: {err:.2e}  (OK)")

    _, _, S0 = ptf.fit_all_sites(T_GRID, P_new, P_new)
    assert abs(S0) < 1e-4, f"null control S(A,A) = {S0}"
    print(f"  A-vs-A null control: Σ ln α = {S0:+.2e}  (OK, optimizer floor)")
    print()


# ---------------------------------------------------------------------------
# Panel 2/3: the main grid at one N
# ---------------------------------------------------------------------------
def run_grid(N, field_site=0):
    bloch = ptf.site_bloch_ops(N)
    ham = ptf.hamming_matrix(N)
    psi0 = fw.bonding_mode_pair_state(N, 1)
    rho0 = np.outer(psi0, psi0.conj())
    H_A = ptf.build_H_xy(N, [J_UNIFORM] * (N - 1))
    Z_m = fw.site_op(N, field_site, 'Z')

    cells = {}        # (eps, gamma) -> {'S': {variant: S}, 'rmse': float}
    floors = {}       # gamma -> {variant: |S| integration floor}
    nulls = {}        # gamma -> S(A,A) standard
    baselines = {}    # gamma -> P_A (reused by panels 4/5)

    for gamma in GAMMA_GRID:
        P_A = purity_traj(N, H_A, rho0, bloch, ham, gamma)
        baselines[gamma] = P_A

        _, _, S0 = ptf.fit_all_sites(T_GRID, P_A, P_A)
        nulls[gamma] = S0
        assert abs(S0) < 1e-4, f"null control broke at N={N}, γ={gamma}: {S0}"

        # integration-noise floor: same physics, halved RK4 step
        P_A_fine = purity_traj(N, H_A, rho0, bloch, ham, gamma,
                               dt=DT / 2, n_steps=2 * N_STEPS)[::2]
        fl, _ = fit_S_variants(P_A, P_A_fine)
        floors[gamma] = {v: abs(fl[v]) for v in VARIANTS}

        for eps in EPS_GRID[1:]:
            H_B = H_A + eps * Z_m
            P_B = purity_traj(N, H_B, rho0, bloch, ham, gamma)
            Sv, rmse = fit_S_variants(P_A, P_B)
            cells[(eps, gamma)] = {'S': Sv, 'rmse': rmse}

    return {'cells': cells, 'floors': floors, 'nulls': nulls,
            'baselines': baselines, 'H_A': H_A, 'rho0': rho0,
            'bloch': bloch, 'ham': ham}


def report_grid(N, grid, field_site=0):
    cells, floors, nulls = grid['cells'], grid['floors'], grid['nulls']

    print("=" * 78)
    print(f"GRID  N={N}, Z-field on site {field_site}: S(ε,γ) = Σ ln α "
          f"(baseline vs +ε·Z_{field_site})")
    print("=" * 78)
    print("  ΔS = S(ε,γ) − S(0,γ); S(0,γ) = A-vs-A null ≈ 0 by protocol "
          "(Z-field only in run B).")
    for gamma in GAMMA_GRID:
        print(f"  γ={gamma}: null S(0,γ)={nulls[gamma]:+.2e}, "
              f"integration floor |S|={floors[gamma]['std']:.2e} (std window)")
    print()

    header = (f"  {'ε':>6s} {'γ':>7s} {'S (std)':>12s} {'ΔS':>12s} "
              f"{'ΔS/(ε²γ)':>11s} {'maxRMSE':>9s} {'win-spread':>11s}  flag")
    print(header)
    for gamma in GAMMA_GRID:
        for eps in EPS_GRID[1:]:
            c = cells[(eps, gamma)]
            S = c['S']['std']
            dS = S - nulls[gamma]
            coll = dS / (eps * eps * gamma)
            spread = max(abs(c['S'][v] - S) for v in VARIANTS[1:])
            floor = floors[gamma]['std']
            below = abs(dS) < max(5.0 * floor, 1e-6)
            shaky = spread > 0.5 * abs(dS)
            flag = ("BELOW-NOISE" if below
                    else "window-sensitive" if shaky else "")
            print(f"  {eps:>6.2f} {gamma:>7.3f} {S:>+12.3e} {dS:>+12.3e} "
                  f"{coll:>11.3f} {c['rmse']:>9.1e} {spread:>11.1e}  {flag}")
        print()

    # --- sign structure ---
    print("  SIGN of ΔS over the grid (rows ε, cols γ):")
    print(f"    {'ε \\ γ':>8s}" + "".join(f"{g:>8.3f}" for g in GAMMA_GRID))
    for eps in EPS_GRID[1:]:
        row = ""
        for gamma in GAMMA_GRID:
            dS = cells[(eps, gamma)]['S']['std'] - nulls[gamma]
            row += f"{'+' if dS > 0 else '−':>8s}"
        print(f"    {eps:>8.2f}" + row)
    print()

    # --- (a) ε exponent at fixed γ, per window variant ---
    print("  (a) ε-exponent: ln|ΔS| ~ a·ln ε at fixed γ  (hypothesis a ≈ 2)")
    print(f"      {'γ':>7s} {'window':>8s} {'ε-range':>16s} {'slope a':>9s} {'R²':>9s}")
    for gamma in GAMMA_GRID:
        for variant in VARIANTS:
            for label, eps_set in (("full", [e for e in EPS_GRID[1:]]),
                                   ("ε≤0.1", [e for e in EPS_GRID[1:] if e <= 0.10])):
                eps_use, dS_use = [], []
                for eps in eps_set:
                    dS = cells[(eps, gamma)]['S'][variant] - nulls[gamma]
                    if abs(dS) > max(5.0 * floors[gamma][variant], 1e-6):
                        eps_use.append(eps)
                        dS_use.append(dS)
                if len(eps_use) >= 3:
                    a, _, r2 = loglog_fit(eps_use, dS_use)
                    rng = f"[{min(eps_use)},{max(eps_use)}]"
                    print(f"      {gamma:>7.3f} {variant:>8s} {rng:>16s} "
                          f"{a:>9.3f} {r2:>9.5f}   ({label}, n={len(eps_use)})")
                else:
                    print(f"      {gamma:>7.3f} {variant:>8s} "
                          f"{'<3 pts above noise':>16s}   ({label})")
    print()

    # --- (b) γ exponent at fixed ε (standard window) ---
    print("  (b) γ-exponent: ln|ΔS| ~ c·ln γ at fixed ε  (hypothesis c ≈ 1)")
    print(f"      {'ε':>6s} {'window':>8s} {'slope c':>9s} {'R²':>9s}")
    for eps in EPS_GRID[1:]:
        for variant in VARIANTS:
            g_use, dS_use = [], []
            for gamma in GAMMA_GRID:
                dS = cells[(eps, gamma)]['S'][variant] - nulls[gamma]
                if abs(dS) > max(5.0 * floors[gamma][variant], 1e-6):
                    g_use.append(gamma)
                    dS_use.append(dS)
            if len(g_use) >= 3:
                c_exp, _, r2 = loglog_fit(g_use, dS_use)
                print(f"      {eps:>6.2f} {variant:>8s} {c_exp:>9.3f} {r2:>9.5f}")
            else:
                print(f"      {eps:>6.2f} {variant:>8s}   <3 γ above noise")
    print()

    # --- (c) collapse table ---
    print("  (c) COLLAPSE ΔS/(ε²γ) (std window; constant across the grid ⟺ the leak")
    print("      rides the m=3 channel):")
    print(f"    {'ε \\ γ':>8s}" + "".join(f"{g:>10.3f}" for g in GAMMA_GRID))
    for eps in EPS_GRID[1:]:
        row = ""
        for gamma in GAMMA_GRID:
            dS = cells[(eps, gamma)]['S']['std'] - nulls[gamma]
            row += f"{dS / (eps * eps * gamma):>10.3f}"
        print(f"    {eps:>8.2f}" + row)
    print()


# ---------------------------------------------------------------------------
# Panel 4: center-site robustness (does the ε exponent care which site?)
# ---------------------------------------------------------------------------
def panel4(N=5, gamma=0.05):
    center = N // 2
    print("=" * 78)
    print(f"PANEL 4  Robustness: Z on the CENTER site {center} (N={N}, γ={gamma})")
    print("=" * 78)
    grid = run_grid(N, field_site=center)
    cells, floors, nulls = grid['cells'], grid['floors'], grid['nulls']
    print(f"  {'ε':>6s} {'ΔS':>12s} {'ΔS/(ε²γ)':>11s}")
    eps_use, dS_use = [], []
    for eps in EPS_GRID[1:]:
        dS = cells[(eps, gamma)]['S']['std'] - nulls[gamma]
        print(f"  {eps:>6.2f} {dS:>+12.3e} {dS / (eps * eps * gamma):>11.3f}")
        if abs(dS) > max(5.0 * floors[gamma]['std'], 1e-6):
            eps_use.append(eps)
            dS_use.append(dS)
    if len(eps_use) >= 3:
        a, _, r2 = loglog_fit(eps_use, dS_use)
        print(f"  ε-exponent (above noise, n={len(eps_use)}): a={a:.3f}, R²={r2:.5f}")
    # also report the γ-exponent across the full grid this run produced
    for eps in (0.10, 0.20):
        g_use = [g for g in GAMMA_GRID]
        dS_g = [cells[(eps, g)]['S']['std'] - nulls[g] for g in g_use]
        c_exp, _, r2 = loglog_fit(g_use, dS_g)
        print(f"  γ-exponent at ε={eps}: c={c_exp:.3f}, R²={r2:.5f}")
    print()
    return grid


# ---------------------------------------------------------------------------
# Panel 5: the other design - Z-field in BOTH runs of the J-defect closure
# ---------------------------------------------------------------------------
def panel5(N=5, gamma=0.05):
    print("=" * 78)
    print(f"PANEL 5  Other design: Z-field ε·Z₀ in BOTH runs of the δJ={DELTA_J} "
          f"J-closure (N={N}, γ={gamma})")
    print("=" * 78)
    print("  S_J(ε) = Σ ln α between A' = uniform+ε·Z₀ and B' = A' + δJ on bond "
          f"{DEFECT_BOND}; ΔS_B = S_J(ε) − S_J(0).")
    bloch = ptf.site_bloch_ops(N)
    ham = ptf.hamming_matrix(N)
    psi0 = fw.bonding_mode_pair_state(N, 1)
    rho0 = np.outer(psi0, psi0.conj())
    Z0 = fw.site_op(N, 0, 'Z')
    J_list_B = [J_UNIFORM] * (N - 1)
    J_list_B[DEFECT_BOND[0]] = J_UNIFORM + DELTA_J
    H_A = ptf.build_H_xy(N, [J_UNIFORM] * (N - 1))
    H_B = ptf.build_H_xy(N, J_list_B)

    S_of_eps = {}
    print(f"  {'ε':>6s} {'S_J(ε)':>12s} {'ΔS_B':>12s} {'ΔS_B/(ε²γ)':>12s} {'maxRMSE':>9s}")
    for eps in EPS_GRID:
        P_Ap = purity_traj(N, H_A + eps * Z0, rho0, bloch, ham, gamma)
        P_Bp = purity_traj(N, H_B + eps * Z0, rho0, bloch, ham, gamma)
        _, r, S = ptf.fit_all_sites(T_GRID, P_Ap, P_Bp)
        S_of_eps[eps] = S
        if eps == 0.0:
            print(f"  {eps:>6.2f} {S:>+12.3e} {'—':>12s} {'—':>12s} {r.max():>9.1e}"
                  f"   (the committed control-arm closure)")
        else:
            dS = S - S_of_eps[0.0]
            print(f"  {eps:>6.2f} {S:>+12.3e} {dS:>+12.3e} "
                  f"{dS / (eps * eps * gamma):>12.3f} {r.max():>9.1e}")
    eps_use = [e for e in EPS_GRID[1:] if abs(S_of_eps[e] - S_of_eps[0.0]) > 1e-6]
    if len(eps_use) >= 3:
        a, _, r2 = loglog_fit(eps_use, [S_of_eps[e] - S_of_eps[0.0] for e in eps_use])
        print(f"  ε-exponent of ΔS_B (n={len(eps_use)}): a={a:.3f}, R²={r2:.5f}")
    print()


# ---------------------------------------------------------------------------
# Panel 6: small-ε refinement - leading order, linear coefficient, residual
# ---------------------------------------------------------------------------
def panel6(N=5):
    print("=" * 78)
    print(f"PANEL 6  Small-ε refinement (N={N}, Z on site 0): the leak's true order")
    print("=" * 78)
    bloch = ptf.site_bloch_ops(N)
    ham = ptf.hamming_matrix(N)
    psi0 = fw.bonding_mode_pair_state(N, 1)
    rho0 = np.outer(psi0, psi0.conj())
    H_A = ptf.build_H_xy(N, [J_UNIFORM] * (N - 1))
    Z0 = fw.site_op(N, 0, 'Z')

    eps_small = [0.005, 0.01, 0.02, 0.05, 0.10]
    for gamma in GAMMA_GRID:
        P_A = purity_traj(N, H_A, rho0, bloch, ham, gamma)
        ratio = []
        for eps in eps_small:
            P_B = purity_traj(N, H_A + eps * Z0, rho0, bloch, ham, gamma)
            _, _, S = ptf.fit_all_sites(T_GRID, P_A, P_B)
            ratio.append(S / eps)
        a, _, r2 = loglog_fit(eps_small, [r * e for r, e in zip(ratio, eps_small)])
        k = 2.0 * ratio[0] - ratio[1]   # Richardson: removes the O(ε) part of dS/ε
        resid = ", ".join(f"{e}:{(r - k) / e:+.3f}"
                          for e, r in zip(eps_small, ratio))
        print(f"  γ={gamma}:")
        print("    dS/ε at ε=" + ", ".join(f"{e}:{r:+.4f}"
                                           for e, r in zip(eps_small, ratio)))
        print(f"    ε-exponent: a={a:.3f} (R²={r2:.5f});  linear coeff "
              f"k=lim dS/ε = {k:+.5f}")
        print(f"    Richardson residual (dS/ε − k)/ε: {resid}")
    print()

    # γ = 0: the fit is DEGENERATE (no decay ⟹ the α-ansatz has nothing to fit;
    # dS locks to an ε-independent constant). Reported, not used as evidence.
    P_A0 = purity_traj(N, H_A, rho0, bloch, ham, 0.0)
    vals = []
    for eps in (0.005, 0.05):
        P_B = purity_traj(N, H_A + eps * Z0, rho0, bloch, ham, 0.0)
        _, _, S = ptf.fit_all_sites(T_GRID, P_A0, P_B)
        vals.append(S)
    print(f"  γ=0 degeneracy check: S(ε=0.005)={vals[0]:+.4f}, "
          f"S(ε=0.05)={vals[1]:+.4f}  (ε-independent lock ⟹ fit degenerate, "
          "excluded from scaling)")
    print()

    # per-site first-order profile
    gamma = 0.05
    P_A = purity_traj(N, H_A, rho0, bloch, ham, gamma)
    print(f"  Per-site first-order profile ln α_i / ε at γ={gamma}:")
    for site in (0, N // 2):
        Zm = fw.site_op(N, site, 'Z')
        for eps in (0.01, 0.02):
            P_B = purity_traj(N, H_A + eps * Zm, rho0, bloch, ham, gamma)
            al, _, S = ptf.fit_all_sites(T_GRID, P_A, P_B)
            prof = np.log(al) / eps
            print(f"    field@site{site}, ε={eps}: "
                  f"{np.array2string(prof, precision=3)}  (Σ·ε = S = {S:+.4e})")
    print("  Stable profiles across ε ⟹ a genuine per-site FIRST-order α response;")
    print("  S is its imperfect cancellation. Sign and magnitude depend on the")
    print("  field site (site 0: mixed signs; center: all positive, symmetric).")
    print()


# ---------------------------------------------------------------------------
def main():
    panel1()

    t0 = time.time()
    grid5 = run_grid(5, field_site=0)
    t5 = time.time() - t0
    report_grid(5, grid5, field_site=0)
    print(f"  [N=5 grid wall time: {t5:.1f} s]")
    print()

    if t5 < 120.0:
        t0 = time.time()
        grid6 = run_grid(6, field_site=0)
        t6 = time.time() - t0
        report_grid(6, grid6, field_site=0)
        print(f"  [N=6 grid wall time: {t6:.1f} s]")
        print()
    else:
        print("  [N=6 skipped: N=5 grid exceeded 120 s]")
        print()

    panel4(N=5, gamma=0.05)
    panel5(N=5, gamma=0.05)
    panel6(N=5)

    print("=" * 78)
    print("Done. Verdict to be read off the exponent fits + collapse tables above.")
    print("=" * 78)


if __name__ == "__main__":
    main()
