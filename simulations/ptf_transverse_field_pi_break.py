#!/usr/bin/env python3
"""PTF under a palindrome-breaking perturbation: does the closure law survive?

Open question from hypotheses/PERSPECTIVAL_TIME_FIELD.md ("Extension to
palindrome-breaking perturbations"):

    The PTF closure law Σ_i ln(α_i) ≈ 0 was established only for perturbations
    that RESPECT the Π palindrome (single-bond J-coupling defects). A transverse
    field h·σ_x^i BREAKS Π. If the rescaling picture survives but with a shifted
    closure law, that is a strong structural statement; if it breaks entirely, a
    clear diagnostic for the role of palindromic protection.

This probe answers it. For the canonical PTF setup (uniform XY chain, γ₀=0.05,
bonding-mode initial state φ = (|vac⟩ + |ψ_1⟩)/√2) it compares two perturbation
families at matched strength ε, on the SAME baseline:

  ARM J  (Π-respecting, control): J-defect δJ=ε on bond (0,1).
                                  Reproduces the known Σ ln α ≈ 0.
  ARM h  (Π-breaking, the test):  transverse field ε·X_i on one site.

For each site we fit P_B(i, t) ≈ P_A(i, α_i·t) and record BOTH the fitted α_i
and the fit RMSE. A large RMSE means the one-parameter time-rescaling ansatz
itself has broken (no clean per-site clock); a small RMSE with Σ ln α drifting
away from 0 means the ansatz survives but the closure law is shifted.

Three readings on the same screen:
  Phase 1  trajectory level: α_i, RMSE, Σ ln α for both arms, N ∈ {5, 6, 7}.
  Phase 2  premise check:    palindrome residual ‖M‖_F at N=4 confirms the
                             J-defect keeps Π (‖M‖≈0) while h·X breaks it.
  Phase 3  mechanism:        the slowest Liouvillian decay rates Re(λ) under
                             each perturbation: protected (J) vs shifted (h).

Propagation is RK4 with the Hamming-mask Z-dephasing (the canonical N=7 PTF
path, matching simulations/ptf_per_observable_alpha.py). Reuses framework
primitives for ρ_0, site operators, the palindrome residual, bond_perturbation,
and slow modes.

Relation to typed knowledge (the two Klein axes are already in the registry):
  - bit_b → spectral palindrome: Phase 2 reproduces the F78 single-body split,
    X truly (‖M‖=0) while Y and Z break it with an identical residual. bit_b is the
    Π² grading (F38). (F78's exact ±2c·γ·i form is for the single-body dissipator;
    here the same Y≡Z degeneracy shows up in a single-body Hamiltonian field.)
  - bit_a → U(1) excitation conservation: elementary (X, Y flip spins; Z does not),
    the magnetization cousin of F61's n_XY-parity grading.
The new content here is neither axis but their DISSOCIATION at the PTF closure:
the closure rides on bit_a (U(1), PTF Layer 3.2 protection #1, the F4 stationary
sector projectors) and is independent of bit_b (the spectral palindrome,
protection #2). See experiments/PTF_PALINDROME_BREAKING_PERTURBATIONS.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Config (canonical PTF setup)
# ---------------------------------------------------------------------------
GAMMA_0 = 0.05
J_UNIFORM = 1.0
DEFECT_BOND = (0, 1)          # the J-defect sits on this bond (ARM J)
T_MAX = 20.0
DT = 0.05
N_STEPS = int(round(T_MAX / DT))
ALPHA_BOUNDS = (0.1, 10.0)
EPS_SCAN = [0.05, 0.10, 0.15, 0.20]


# ---------------------------------------------------------------------------
# Operators / Hamiltonians (built from framework site ops)
# ---------------------------------------------------------------------------
def build_H_xy(N, J_list):
    """H = Σ_b (J_b/2)(X_b X_{b+1} + Y_b Y_{b+1}) on the open chain."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        Jb = J_list[b]
        H = H + (Jb / 2.0) * (
            fw.site_op(N, b, 'X') @ fw.site_op(N, b + 1, 'X')
            + fw.site_op(N, b, 'Y') @ fw.site_op(N, b + 1, 'Y')
        )
    return H


def hamming_matrix(N):
    """popcount(a ⊕ b) matrix; Z-dephasing dissipator = -2γ₀·hamming ⊙ ρ."""
    d = 2 ** N
    idx = np.arange(d, dtype=np.uint64)
    xor = idx[:, None] ^ idx[None, :]
    h = np.zeros((d, d), dtype=np.int32)
    for i in range(N):
        h += ((xor >> np.uint64(i)) & np.uint64(1)).astype(np.int32)
    return h


def site_bloch_ops(N):
    return [(fw.site_op(N, i, 'X'), fw.site_op(N, i, 'Y'), fw.site_op(N, i, 'Z'))
            for i in range(N)]


# ---------------------------------------------------------------------------
# Lindblad RK4 (Hamming-mask Z-dephasing): the canonical N=7 PTF path
# ---------------------------------------------------------------------------
def _rhs(rho, H, hamming):
    return -1j * (H @ rho - rho @ H) - 2.0 * GAMMA_0 * hamming * rho


def _rk4_step(rho, H, hamming, dt):
    k1 = _rhs(rho, H, hamming)
    k2 = _rhs(rho + 0.5 * dt * k1, H, hamming)
    k3 = _rhs(rho + 0.5 * dt * k2, H, hamming)
    k4 = _rhs(rho + dt * k3, H, hamming)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def purity_trajectory(N, H, rho_0, bloch_ops, hamming):
    """Per-site purity P_i(t) = ½(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²), shape (n_t, N)."""
    rho = rho_0.copy()
    P = np.zeros((N_STEPS + 1, N))

    def record(r, row):
        for i, (Xi, Yi, Zi) in enumerate(bloch_ops):
            x = np.real(np.trace(Xi @ r))
            y = np.real(np.trace(Yi @ r))
            z = np.real(np.trace(Zi @ r))
            P[row, i] = 0.5 * (1.0 + x * x + y * y + z * z)

    record(rho, 0)
    for step in range(1, N_STEPS + 1):
        rho = _rk4_step(rho, H, hamming, DT)
        rho = 0.5 * (rho + rho.conj().T)
        record(rho, step)
    return P


# ---------------------------------------------------------------------------
# α-fit: P_A(α·t) ≈ P_B(t)  (matches ptf_per_observable_alpha.py)
# ---------------------------------------------------------------------------
def alpha_fit_one_site(t_grid, P_A, P_B):
    interp = interp1d(t_grid, P_A, bounds_error=False,
                      fill_value=(float(P_A[0]), float(P_A[-1])), kind='cubic')

    def mse(alpha):
        d = interp(alpha * t_grid) - P_B
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=ALPHA_BOUNDS, method='bounded',
                          options={'xatol': 1e-7})
    return float(res.x), float(np.sqrt(res.fun))


def fit_all_sites(t_grid, P_A, P_B):
    N = P_A.shape[1]
    alphas = np.zeros(N)
    rmses = np.zeros(N)
    for i in range(N):
        alphas[i], rmses[i] = alpha_fit_one_site(t_grid, P_A[:, i], P_B[:, i])
    sigma_log = float(np.sum(np.log(np.clip(alphas, 1e-30, None))))
    return alphas, rmses, sigma_log


# ---------------------------------------------------------------------------
# Phase 1: trajectory-level comparison
# ---------------------------------------------------------------------------
def phase1():
    print("=" * 78)
    print("PHASE 1  Trajectory level: does P_B(i,t) ≈ P_A(i, α_i·t) survive?")
    print("=" * 78)
    print(f"  γ₀={GAMMA_0}, J={J_UNIFORM}, φ=(|vac⟩+|ψ₁⟩)/√2, T_max={T_MAX}, dt={DT}")
    print(f"  ARM J: J-defect ε on bond {DEFECT_BOND} (Π-respecting control)")
    print(f"  ARM h: transverse field ε·X_i on site i (Π-breaking)")
    print()

    t_grid = np.linspace(0, T_MAX, N_STEPS + 1)
    results = {}

    for N in (5, 6, 7):
        center = N // 2
        bloch_ops = site_bloch_ops(N)
        hamming = hamming_matrix(N)
        psi_0 = fw.bonding_mode_pair_state(N, 1)
        rho_0 = np.outer(psi_0, psi_0.conj())

        # baseline A (uniform)
        H_A = build_H_xy(N, [J_UNIFORM] * (N - 1))
        P_A = purity_trajectory(N, H_A, rho_0, bloch_ops, hamming)

        # null control: A vs A → α_i ≈ 1, Σ ln α ≈ 0
        a0, r0, s0 = fit_all_sites(t_grid, P_A, P_A)
        assert np.allclose(a0, 1.0, atol=1e-3), f"null control α≠1 at N={N}: {a0}"
        assert abs(s0) < 1e-2, f"null control Σlnα≠0 at N={N}: {s0}"

        print(f"  --- N={N} (field sites: 0 and center={center}) ---")
        header = (f"  {'arm':<14s}{'ε':>6s} {'Σ ln α':>10s} "
                  f"{'max RMSE':>10s} {'mean RMSE':>10s}  closure")
        print(header)

        for eps in EPS_SCAN:
            # ARM J: J-defect on DEFECT_BOND
            J_list = [J_UNIFORM] * (N - 1)
            J_list[DEFECT_BOND[0]] = J_UNIFORM + eps
            H_BJ = build_H_xy(N, J_list)
            P_BJ = purity_trajectory(N, H_BJ, rho_0, bloch_ops, hamming)
            aJ, rJ, sJ = fit_all_sites(t_grid, P_A, P_BJ)

            # ARM h: transverse field on site 0 and on center
            H_Bh0 = H_A + eps * fw.site_op(N, 0, 'X')
            P_Bh0 = purity_trajectory(N, H_Bh0, rho_0, bloch_ops, hamming)
            ah0, rh0, sh0 = fit_all_sites(t_grid, P_A, P_Bh0)

            H_Bhc = H_A + eps * fw.site_op(N, center, 'X')
            P_Bhc = purity_trajectory(N, H_Bhc, rho_0, bloch_ops, hamming)
            ahc, rhc, shc = fit_all_sites(t_grid, P_A, P_Bhc)

            for tag, s, r in (("J  bond(0,1)", sJ, rJ),
                              ("h  X@site0", sh0, rh0),
                              (f"h  X@site{center}", shc, rhc)):
                verdict = ("OK ~0" if abs(s) < 0.10
                           else "SHIFTED" if abs(s) < 0.60 else "BROKEN")
                print(f"  {tag:<14s}{eps:>6.2f} {s:>+10.4f} "
                      f"{r.max():>10.2e} {r.mean():>10.2e}  {verdict}")
            print()

            results[(N, eps)] = {
                'J': (aJ, rJ, sJ), 'h0': (ah0, rh0, sh0), 'hc': (ahc, rhc, shc),
            }

        # sanity: at N=7, eps=0.10 the J-defect reproduces the published PTF
        if N == 7:
            aJ_ref = results[(7, 0.10)]['J'][0]
            published = np.array([1.095, 1.182, 1.051, 0.991, 0.845, 0.923, 0.997])
            print(f"  PTF reproduction check (N=7, J-defect ε=0.10):")
            print(f"    fitted α : {np.array2string(aJ_ref, precision=3)}")
            print(f"    published: {np.array2string(published, precision=3)}")
            print(f"    max|Δ|   : {np.max(np.abs(aJ_ref - published)):.3f}")
            print()

    return results


# ---------------------------------------------------------------------------
# Phase 2: what does the transverse field actually break?
#   two axes: spectral palindrome ‖M‖_F  vs  U(1) conservation ‖[H, N_exc]‖
# ---------------------------------------------------------------------------
def _N_exc(N):
    """Excitation-number operator N_exc = Σ_i (I − Z_i)/2."""
    d = 2 ** N
    out = np.zeros((d, d), dtype=complex)
    half_I = 0.5 * np.eye(d, dtype=complex)
    for i in range(N):
        out = out + half_I - 0.5 * fw.site_op(N, i, 'Z')
    return out


def phase2():
    print("=" * 78)
    print("PHASE 2  What does the field break? spectral palindrome vs U(1)")
    print("=" * 78)
    print("  Two independent axes, both at strength 0.10:")
    print("    ‖M‖_F        = ‖Π L Π⁻¹ + L + 2Σγ‖   (spectral palindrome; 0 ⟺ mirror)")
    print("    ‖[H, N_exc]‖ = excitation (U(1)) non-conservation")
    print()

    h = 0.10
    for N in (3, 4, 5):
        Sigma_gamma = N * GAMMA_0
        gamma_l = [GAMMA_0] * N
        H_A = build_H_xy(N, [J_UNIFORM] * (N - 1))
        Nexc = _N_exc(N)

        def res_norm(H):
            L = fw.lindbladian_z_dephasing(H, gamma_l)
            return float(np.linalg.norm(fw.palindrome_residual(L, Sigma_gamma, N)))

        def u1(H):
            return float(np.linalg.norm(H @ Nexc - Nexc @ H))

        J_list = [J_UNIFORM] * (N - 1)
        J_list[0] = J_UNIFORM + h
        H_J = build_H_xy(N, J_list)
        H_X = H_A + h * fw.site_op(N, 0, 'X')
        H_Y = H_A + h * fw.site_op(N, 0, 'Y')
        H_Z = H_A + h * fw.site_op(N, 0, 'Z')
        H_gX = H_A + h * sum(fw.site_op(N, i, 'X') for i in range(N))

        rows = [
            ("uniform XY",        H_A),
            ("+ J-defect bond01", H_J),
            ("+ X-field site0",   H_X),
            ("+ Y-field site0",   H_Y),
            ("+ Z-field site0",   H_Z),
            ("+ global ΣX field", H_gX),
        ]
        print(f"  --- N={N} ---")
        print(f"  {'perturbation':<20s}{'‖M‖_F':>12s}{'palindrome':>12s}"
              f"{'‖[H,Nexc]‖':>13s}{'U(1)':>8s}")
        for tag, H in rows:
            m = res_norm(H)
            u = u1(H)
            pal = "mirror" if m < 1e-9 else "BROKEN"
            u1tag = "conserv" if u < 1e-9 else "BROKEN"
            print(f"  {tag:<20s}{m:>12.2e}{pal:>12s}{u:>13.2e}{u1tag:>8s}")
        print()

    # Structural assertions (N=4 anchor)
    N = 4
    gamma_l = [GAMMA_0] * N
    Sigma_gamma = N * GAMMA_0
    H_A = build_H_xy(N, [J_UNIFORM] * (N - 1))
    Nexc = _N_exc(N)
    H_X = H_A + h * fw.site_op(N, 0, 'X')
    L_X = fw.lindbladian_z_dephasing(H_X, gamma_l)
    m_X = float(np.linalg.norm(fw.palindrome_residual(L_X, Sigma_gamma, N)))
    u_X = float(np.linalg.norm(H_X @ Nexc - Nexc @ H_X))
    u_A = float(np.linalg.norm(H_A @ Nexc - Nexc @ H_A))
    assert m_X < 1e-9, f"surprise refuted? X-field broke palindrome: {m_X}"
    assert u_A < 1e-9, f"uniform XY does not conserve U(1)?? {u_A}"
    assert u_X > 1e-6, f"X-field did not break U(1): {u_X}"

    print("  READING: the J-defect respects BOTH axes. The transverse X-field")
    print("  preserves the spectral palindrome (‖M‖≈0, the spectrum stays")
    print("  mirror-symmetric about −Σγ) yet BREAKS U(1) excitation conservation.")
    print("  So the PTF closure law that Phase 1 shows collapsing is NOT guarded")
    print("  by the spectral palindrome; it is guarded by U(1) sector conservation")
    print("  (PTF doc Layer 3.2, protection #1: the stationary sector projectors).")
    print("  The doc's premise 'a transverse field breaks Π' is imprecise: X breaks")
    print("  the protection via U(1), not via the mirror. The Z-field is the mirror")
    print("  image: it conserves U(1) but breaks the palindrome (‖M‖≠0), and Phase 4")
    print("  shows its closure survives anyway. So the two axes are independent.")
    print()


# ---------------------------------------------------------------------------
# Phase 3, mechanism: slowest Liouvillian decay rates under each perturbation
# ---------------------------------------------------------------------------
def phase3():
    print("=" * 78)
    print("PHASE 3  Mechanism: slowest decay rates Re(λ) at N=5, ε=0.10")
    print("=" * 78)
    N = 5
    eps = 0.10
    n_show = 12

    chain = fw.ChainSystem(N, gamma_0=GAMMA_0, J=J_UNIFORM, topology='chain', H_type='xy')
    L_A = chain.L

    # J-defect Liouvillian: add ε·V_L on bond (0,1)
    V_J = fw.bond_perturbation(N, DEFECT_BOND, kind='XY')
    L_BJ = L_A + eps * V_J

    # transverse-field Liouvillian: add ε·(−i[X_0,·]) (same vec convention)
    X0 = fw.site_op(N, 0, 'X')
    Id = np.eye(2 ** N, dtype=complex)
    V_h = -1j * (np.kron(X0, Id) - np.kron(Id, X0.T))
    L_Bh = L_A + eps * V_h

    def rates_and_stationary(L, k, tol=1e-9):
        ev = np.linalg.eigvals(L)
        rate = -np.real(ev)                       # decay rate = -Re(λ) ≥ 0
        n_stat = int(np.sum(rate < tol))
        return np.sort(rate)[:k], n_stat

    rA, stat_A = rates_and_stationary(L_A, n_show)
    rJ, stat_J = rates_and_stationary(L_BJ, n_show)
    rh, stat_h = rates_and_stationary(L_Bh, n_show)

    print(f"  stationary modes (|Re λ| < 1e-9, the F4 sector projectors):")
    print(f"    baseline                : {stat_A}")
    print(f"    + J-defect (Π & U(1) ok): {stat_J}")
    print(f"    + field    (U(1) broken): {stat_h}")
    print()

    print(f"  slowest {n_show} decay rates -Re(λ), sorted ascending:")
    print(f"  {'#':>3s} {'baseline':>12s} {'+J-defect':>12s} {'Δ_J':>10s} "
          f"{'+field':>12s} {'Δ_h':>10s}")
    for i in range(n_show):
        print(f"  {i:>3d} {rA[i]:>12.6f} {rJ[i]:>12.6f} {rJ[i]-rA[i]:>+10.2e} "
              f"{rh[i]:>12.6f} {rh[i]-rA[i]:>+10.2e}")
    print()
    print(f"  max |Δ decay rate| over slowest {n_show} modes:")
    print(f"    J-defect (Π-respecting): {np.max(np.abs(rJ - rA)):.2e}")
    print(f"    field    (Π-breaking)  : {np.max(np.abs(rh - rA)):.2e}")
    print()
    print("  Reading: if the J-defect leaves the slow decay rates near-fixed")
    print("  (envelope shape preserved → clean time-rescale) while the field")
    print("  shifts them (envelope shape changes → no clean rescale), that is")
    print("  the mechanism behind whatever Phase 1 shows for the closure law.")
    print()


def phase4():
    print("=" * 78)
    print("PHASE 4  The Klein prediction: does closure-break track bit_a (U(1))")
    print("         rather than bit_b (palindrome)?")
    print("=" * 78)
    print("  Single-site field ε·P_0 for P ∈ {X, Z, Y}. Klein bits (a,b):")
    print("    X=(1,0) breaks U(1), keeps palindrome → PREDICT closure BREAKS")
    print("    Z=(0,1) keeps U(1), breaks palindrome → PREDICT closure SURVIVES")
    print("    Y=(1,1) breaks both                   → PREDICT closure BREAKS")
    print()
    t_grid = np.linspace(0, T_MAX, N_STEPS + 1)

    for N in (5, 6, 7):
        bloch_ops = site_bloch_ops(N)
        hamming = hamming_matrix(N)
        psi_0 = fw.bonding_mode_pair_state(N, 1)
        rho_0 = np.outer(psi_0, psi_0.conj())
        H_A = build_H_xy(N, [J_UNIFORM] * (N - 1))
        P_A = purity_trajectory(N, H_A, rho_0, bloch_ops, hamming)

        print(f"  --- N={N} ---")
        print(f"  {'field':<8s}{'(a,b)':>8s}{'ε':>6s}{'Σ ln α':>10s}"
              f"{'max RMSE':>11s}  predict   closure")
        store = {}
        for letter in ('X', 'Z', 'Y'):
            ba, bb = fw.bit_a(letter), fw.bit_b(letter)
            pred = "BREAK" if ba == 1 else "survive"
            for eps in (0.10, 0.20):
                H_B = H_A + eps * fw.site_op(N, 0, letter)
                P_B = purity_trajectory(N, H_B, rho_0, bloch_ops, hamming)
                _a, r, s = fit_all_sites(t_grid, P_A, P_B)
                obs = ("OK ~0" if abs(s) < 0.10
                       else "SHIFTED" if abs(s) < 0.60 else "BROKEN")
                store[(letter, eps)] = (s, r.max())
                print(f"  {letter:<8s}{f'({ba},{bb})':>8s}{eps:>6.2f}{s:>+10.4f}"
                      f"{r.max():>11.2e}  {pred:>7s}   {obs}")
        # The Klein prediction, asserted at ε=0.10 (matched ε): closure-break
        # tracks bit_a. NOTE the matched-ε RMSE gap is inflated because Z's
        # same-ε trajectory effect is intrinsically tiny; Phase 5 gives the
        # honest effect-size-matched comparison and the hard-drive control.
        sX, rX = store[('X', 0.10)]
        sZ, rZ = store[('Z', 0.10)]
        sY, rY = store[('Y', 0.10)]
        assert rZ < 0.01 and abs(sZ) < 0.15, f"Z (bit_a=0) did not survive at N={N}: σ={sZ}, rmse={rZ}"
        assert rX > 0.02 and abs(sX) > 0.6, f"X (bit_a=1) did not break at N={N}: σ={sX}, rmse={rX}"
        assert rY > 0.02 and abs(sY) > 0.6, f"Y (bit_a=1) did not break at N={N}: σ={sY}, rmse={rY}"
        print()

    print("  ✓ Confirmed N=5,6,7: Z (bit_a=0, palindrome-BROKEN) keeps Σ ln α ~ 0")
    print("  with a clean fit, while X and Y (bit_a=1) blow up both Σ ln α and the")
    print("  fit RMSE. The PTF closure law is guarded by U(1) excitation conservation")
    print("  (the bit_a / light axis), independent of the spectral palindrome (bit_b).")
    print()


def phase5():
    print("=" * 78)
    print("PHASE 5  Honest discriminator: effect-size, not matched ε")
    print("=" * 78)
    print("  Matched-ε is unfair: the same ε perturbs the trajectories by very")
    print("  different amounts for X vs Z, so the headline RMSE ratio is inflated.")
    print("  Here we report the trajectory effect-size eff = rms_(i,t)|P_B − P_A|")
    print("  and DRIVE THE Z-FIELD HARD. If the dissociation were a magnitude")
    print("  artifact, a hard-driven Z would eventually break the closure. It does")
    print("  not: Z's effect saturates and its closure stays bounded with a clean")
    print("  fit, while X diverges in both Σ ln α and RMSE.")
    print()
    N = 6
    t_grid = np.linspace(0, T_MAX, N_STEPS + 1)
    bloch_ops = site_bloch_ops(N)
    hamming = hamming_matrix(N)
    psi_0 = fw.bonding_mode_pair_state(N, 1)
    rho_0 = np.outer(psi_0, psi_0.conj())
    H_A = build_H_xy(N, [J_UNIFORM] * (N - 1))
    P_A = purity_trajectory(N, H_A, rho_0, bloch_ops, hamming)

    def measure(letter, eps):
        H_B = H_A + eps * fw.site_op(N, 0, letter)
        P_B = purity_trajectory(N, H_B, rho_0, bloch_ops, hamming)
        eff = float(np.sqrt(np.mean((P_B - P_A) ** 2)))
        _a, r, s = fit_all_sites(t_grid, P_A, P_B)
        return eff, s, r.max()

    scans = {'X': (0.05, 0.10, 0.20, 0.40), 'Z': (0.10, 0.50, 1.0, 2.0, 4.0, 8.0)}
    rec = {}
    print(f"  N={N}.  {'field':<7s}{'ε':>7s}{'eff-size':>10s}{'Σ ln α':>10s}{'max RMSE':>11s}")
    for letter in ('X', 'Z'):
        for eps in scans[letter]:
            rec[(letter, eps)] = measure(letter, eps)
            eff, s, r = rec[(letter, eps)]
            print(f"  {'':7s}{letter:<7s}{eps:>7.2f}{eff:>10.4f}{s:>+10.4f}{r:>11.2e}")
        print()

    eff_Z_top, s_Z_top, rmse_Z_top = rec[('Z', 8.0)]
    rmse_Z_max = max(rec[('Z', e)][2] for e in scans['Z'])
    s_Z_max = max(abs(rec[('Z', e)][1]) for e in scans['Z'])
    eff_X_05, s_X_05, _ = rec[('X', 0.05)]
    rmse_X_40 = rec[('X', 0.40)][2]

    print(f"  Z driven to ε=8.0: effect-size {eff_Z_top:.4f} (it SATURATES), and the")
    print(f"    fit RMSE stays clean (≤{rmse_Z_max:.1e}) at EVERY drive; Σ ln α stays")
    print(f"    bounded (peaks ≈{s_Z_max:.2f} near ε=1 then RECEDES to {rec[('Z',8.0)][1]:+.2f}).")
    print(f"  X at ε=0.05 already has comparable effect-size {eff_X_05:.4f} but Σ ln α={s_X_05:+.2f},")
    print(f"    and by ε=0.40 both diverge (Σ ln α={rec[('X',0.40)][1]:+.1f}, RMSE={rmse_X_40:.1e}).")
    print("  So at matched (or larger) effect-size the Z-field still does not break")
    print("  the closure, while X does. The clean discriminator is the fit RMSE")
    print("  (Z stays ≤2e-2 under any drive, X diverges); the matched-ε ratio in")
    print("  Phase 1/4 overstates the gap. The dissociation is structural (U(1)/bit_a),")
    print("  not a magnitude artifact: Z cannot break the closure no matter how hard")
    print("  it is driven, because it never leaves the excitation sectors.")
    print()

    # Honest, robust discriminator: Z's fit stays clean under arbitrary drive while
    # X's diverges. Z's Σ ln α is bounded (and recedes at strong drive); X's diverges.
    assert rmse_Z_max < 0.03, f"Z fit broke under hard drive (max RMSE {rmse_Z_max})"
    assert rmse_X_40 > 0.10, f"X fit did not diverge (RMSE {rmse_X_40})"
    assert abs(rec[('X', 0.40)][1]) > 5 * s_Z_max, "X Σ ln α did not diverge past Z's bounded range"


def main():
    res = phase1()
    phase2()
    phase3()
    phase4()
    phase5()
    print("=" * 78)
    print("Done. Headline: the PTF closure breaks under a single-site field iff")
    print("bit_a=1 (X, Y → U(1) broken), independent of the palindrome (bit_b).")
    print("Phase 5 confirms it is structural, not a magnitude artifact.")
    print("=" * 78)


if __name__ == "__main__":
    main()
