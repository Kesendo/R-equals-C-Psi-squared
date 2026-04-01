#!/usr/bin/env python3
"""
Phase 3: Does Palindrome Breaking = Chiral Phase Breaking?
============================================================
Add depolarizing noise to the fragile bridge gain-loss system.
Depolarizing breaks the palindrome (known: err ∝ ε).
Question: does γ_crit shift proportionally to the palindrome error?

Script: simulations/pt_palindrome_breaking.py
Output: simulations/results/pt_palindrome_breaking.txt
"""

import numpy as np
from scipy.linalg import eigvals
from itertools import product as iproduct
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "pt_palindrome_breaking.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ========================================================================
# Pauli infrastructure (self-contained)
# ========================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]

PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def op_at(op, qubit, n_qubits):
    result = np.array([[1]], dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == qubit else I2)
    return result


def build_coupled_liouvillian_depol(n_per_chain, gamma, eps_depol,
                                     J=1.0, J_bridge=1.0):
    """Fragile bridge with Z-dephasing (±γ) plus depolarizing (ε) on all."""
    n_total = 2 * n_per_chain
    d = 2**n_total
    d2 = d * d

    # Hamiltonian
    H = np.zeros((d, d), dtype=complex)
    for chain_start in [0, n_per_chain]:
        for i in range(chain_start, chain_start + n_per_chain - 1):
            for P in [sx, sy, sz]:
                H += J * op_at(P, i, n_total) @ op_at(P, i + 1, n_total)
    for P in [sx, sy, sz]:
        H += J_bridge * op_at(P, n_per_chain - 1, n_total) @ \
             op_at(P, n_per_chain, n_total)

    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))

    # Z-dephasing: +γ on chain A, −γ on chain B
    for k in range(n_per_chain):
        Zk = op_at(sz, k, n_total)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    for k in range(n_per_chain, 2 * n_per_chain):
        Zk = op_at(sz, k, n_total)
        L += (-gamma) * (np.kron(Zk, Zk.conj()) - np.eye(d2))

    # Depolarizing noise on ALL qubits at rate ε
    # Each axis (X, Y, Z) gets rate ε/3
    if eps_depol > 0:
        for k in range(n_total):
            for P in [sx, sy, sz]:
                Pk = op_at(P, k, n_total)
                L += (eps_depol / 3.0) * (
                    np.kron(Pk, Pk.conj()) - np.eye(d2))

    return L


def find_gamma_crit_depol(n_per_chain, J_bridge, eps_depol,
                           J=1.0, tol=1e-7):
    """Find γ_crit for fragile bridge with depolarizing perturbation."""
    def max_re(gamma):
        L = build_coupled_liouvillian_depol(
            n_per_chain, gamma, eps_depol, J, J_bridge)
        return float(np.max(eigvals(L).real))

    # Check if system is ever unstable
    if max_re(10.0) < 1e-10:
        return None

    g_lo, g_hi = 0.0, 10.0
    while max_re(g_hi) < 1e-10:
        g_hi *= 2

    while (g_hi - g_lo) > tol:
        g_mid = (g_lo + g_hi) / 2
        if max_re(g_mid) > 1e-10:
            g_hi = g_mid
        else:
            g_lo = g_mid

    return (g_lo + g_hi) / 2


def measure_palindrome_quality(n_per_chain, gamma, eps_depol,
                                J=1.0, J_bridge=1.0):
    """Measure palindrome quality: eigenvalue pairing error and max|Re|."""
    L = build_coupled_liouvillian_depol(
        n_per_chain, gamma, eps_depol, J, J_bridge)
    ev = eigvals(L)

    # For ε=0: center is 0 (Σγ_Z=0). For ε>0: center shifts.
    # Find optimal center by minimizing pairing error
    # Try c=0 (pure Z) and c based on total dephasing
    best_err = np.inf
    best_c = 0.0
    n_total = 2 * n_per_chain
    # Candidate centers: 0, and (4/3)·n_total·ε (from Z+Y anti-commuting parts)
    candidates = [0.0]
    if eps_depol > 0:
        # Z-depol: ε/3 per qubit → Σγ_Z_depol = n_total·ε/3
        # Y-depol: ε/3 per qubit → Σγ_Y_depol = n_total·ε/3
        # Total anti-commuting: Σγ_Z + Σγ_Y = 0 + 2·n_total·ε/3
        c_anti = 2 * (0 + 2 * n_total * eps_depol / 3)
        candidates.append(c_anti)
        # Also try: just the depol Z contribution
        c_z = 2 * n_total * eps_depol / 3
        candidates.append(c_z)

    for c_try in candidates:
        center = c_try / 2
        pair_errs = []
        for lam in ev:
            target = -(lam + c_try)
            pair_errs.append(np.min(np.abs(ev - target)))
        err = np.max(pair_errs)
        if err < best_err:
            best_err = err
            best_c = c_try

    # Also measure: max|Re| for eigenvalues (stability)
    max_re = np.max(ev.real)

    # Fraction of eigenvalues with |Re| < threshold
    nonzero = ev[np.abs(ev) > 1e-10]
    on_axis_frac = np.mean(np.abs(nonzero.real) < 1e-4) if len(nonzero) > 0 else 1.0

    return best_err, best_c, max_re, on_axis_frac


# ========================================================================
log("=" * 72)
log("PHASE 3: PALINDROME BREAKING AND CHIRAL PHASE STABILITY")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)

N_chain = 2
J = 1.0
J_br = 1.0

# ----------------------------------------------------------------
# 3a. Reference: γ_crit without depolarizing
# ----------------------------------------------------------------
log()
log("─" * 72)
log("3a. Reference system (no depolarizing)")
log("─" * 72)
log()

t0 = clock.time()
gc_ref = find_gamma_crit_depol(N_chain, J_br, 0.0, J=J)
log(f"  γ_crit(ε=0) = {gc_ref:.7f}  ({clock.time()-t0:.1f}s)")

# ----------------------------------------------------------------
# 3b. Sweep ε: palindrome error + γ_crit
# ----------------------------------------------------------------
log()
log("─" * 72)
log("3b. Depolarizing perturbation sweep")
log("─" * 72)
log()
log("  ε = depolarizing rate per qubit (X+Y+Z, each at ε/3)")
log("  Palindrome error measured at γ = 0.5·γ_crit(ε=0) (stable regime)")
log()

gamma_test = 0.5 * gc_ref  # measure palindrome in stable regime

eps_values = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]

log(f"  {'ε':>8}  {'pal. err':>10}  {'opt. c':>8}  {'γ_crit':>10}"
    f"  {'Δγ_c':>10}  {'Δγ_c%':>7}  {'on-axis%':>8}")
log(f"  {'─'*70}")

gc_arr = []
pe_arr = []

for eps in eps_values:
    t0 = clock.time()

    # Palindrome quality at fixed γ in stable regime
    pal_err, opt_c, _, on_axis = measure_palindrome_quality(
        N_chain, gamma_test, eps, J=J, J_bridge=J_br)

    # Find γ_crit
    gc = find_gamma_crit_depol(N_chain, J_br, eps, J=J)

    gc_arr.append(gc)
    pe_arr.append(pal_err)

    if gc is not None:
        delta_gc = gc - gc_ref
        delta_pct = 100 * delta_gc / gc_ref
        log(f"  {eps:>8.4f}  {pal_err:>10.2e}  {opt_c:>8.4f}"
            f"  {gc:>10.7f}  {delta_gc:>+10.7f}  {delta_pct:>+7.2f}%"
            f"  {100*on_axis:>7.1f}%")
    else:
        log(f"  {eps:>8.4f}  {pal_err:>10.2e}  {opt_c:>8.4f}"
            f"  {'never':>10}  {'N/A':>10}  {'N/A':>7}  {100*on_axis:>7.1f}%")

# ----------------------------------------------------------------
# 3c. Correlation analysis
# ----------------------------------------------------------------
log()
log("─" * 72)
log("3c. Correlation: palindrome error vs γ_crit shift")
log("─" * 72)
log()

valid = [(pe, gc) for pe, gc in zip(pe_arr, gc_arr)
         if gc is not None and pe > 1e-12]
if len(valid) >= 3:
    pe_v = np.array([v[0] for v in valid])
    gc_v = np.array([v[1] for v in valid])
    delta_gc_v = gc_v - gc_ref

    # Pearson correlation
    if np.std(pe_v) > 0 and np.std(delta_gc_v) > 0:
        r_corr = np.corrcoef(pe_v, delta_gc_v)[0, 1]
        log(f"  Pearson r(palindrome_error, Δγ_crit) = {r_corr:+.6f}")
        log(f"  N = {len(valid)} data points (excluding ε=0)")
        log()

        if abs(r_corr) > 0.9:
            log("  STRONG correlation: palindrome breaking and γ_crit shift")
            log("  are closely linked. Breaking the palindrome destabilizes")
            log("  the chiral phase.")
        elif abs(r_corr) > 0.5:
            log("  MODERATE correlation: partial link between palindrome")
            log("  quality and chiral phase stability.")
        else:
            log("  WEAK correlation: palindrome quality and γ_crit shift")
            log("  are largely independent. The chiral phase is protected")
            log("  by something other than (or in addition to) the palindrome.")

        # Direction
        if r_corr > 0:
            log()
            log("  Direction: POSITIVE. Larger palindrome error = larger γ_crit.")
            log("  Breaking the palindrome makes the system MORE stable (!)")
            log("  Interpretation: depolarizing noise adds dissipation to the")
            log("  gain side, reducing effective gain. The stabilization effect")
            log("  of added damping dominates over the palindrome breaking.")
        else:
            log()
            log("  Direction: NEGATIVE. Larger palindrome error = smaller γ_crit.")
            log("  Breaking the palindrome makes the system LESS stable.")
            log("  The palindrome IS the protection mechanism.")

    # Linear fit
    if len(valid) >= 3:
        slope, intercept = np.polyfit(pe_v, delta_gc_v, 1)
        log()
        log(f"  Linear fit: Δγ_crit = {slope:.4f} × pal_err + {intercept:.6f}")

# ----------------------------------------------------------------
# 3d. Direct test: does chiral phase (Im axis) break with palindrome?
# ----------------------------------------------------------------
log()
log("─" * 72)
log("3d. Chiral phase quality at fixed γ = 0.5·γ_crit(ε=0)")
log("─" * 72)
log()
log("  In the pure system: ALL eigenvalues on the imaginary axis.")
log("  Does depolarizing noise push eigenvalues off the axis?")
log()

log(f"  {'ε':>8}  {'max|Re(λ)|':>12}  {'pal. err':>10}  {'still stable?':>13}")
log(f"  {'─'*50}")

for eps in eps_values:
    L = build_coupled_liouvillian_depol(
        N_chain, gamma_test, eps, J=J, J_bridge=J_br)
    ev = eigvals(L)
    max_re = np.max(ev.real)
    pal_err = pe_arr[eps_values.index(eps)]
    stable = "YES" if max_re < 1e-6 else f"NO (Re={max_re:.2e})"
    log(f"  {eps:>8.4f}  {np.max(np.abs(ev.real)):>12.2e}"
        f"  {pal_err:>10.2e}  {stable:>13}")

# ----------------------------------------------------------------
# 3e. Instability mechanism: still Hopf?
# ----------------------------------------------------------------
log()
log("─" * 72)
log("3e. Instability mechanism check (ε = 0.05)")
log("─" * 72)
log()

eps_check = 0.05
gc_check = find_gamma_crit_depol(N_chain, J_br, eps_check, J=J)
if gc_check is not None:
    L_at_crit = build_coupled_liouvillian_depol(
        N_chain, gc_check * 1.01, eps_check, J=J, J_bridge=J_br)
    ev_crit = eigvals(L_at_crit)
    idx_max = np.argmax(ev_crit.real)
    lam_unstable = ev_crit[idx_max]

    log(f"  ε = {eps_check}, γ_crit = {gc_check:.6f}")
    log(f"  Most unstable eigenvalue at 1.01×γ_crit:")
    log(f"    λ = {lam_unstable.real:.6f} + {lam_unstable.imag:.6f}i")
    log(f"    |Im(λ)| = {abs(lam_unstable.imag):.4f}")
    log()
    if abs(lam_unstable.imag) > 0.01:
        log("  Im(λ) ≠ 0 at onset: still a HOPF bifurcation.")
        log("  Depolarizing noise does not change the instability type.")
    else:
        log("  Im(λ) ≈ 0 at onset: changed to SADDLE-NODE.")
        log("  Depolarizing noise changes the instability mechanism.")

    # Check pairing: is the partner at −λ?
    target = -lam_unstable
    dists = np.abs(ev_crit - target)
    dists[idx_max] = np.inf
    partner = ev_crit[np.argmin(dists)]
    pair_err = np.min(dists)

    # At ε>0, Σγ≠0, so true partner is at -(λ+c)
    log()
    log(f"  Pairing check (λ↔−λ): |λ + λ'| = {abs(lam_unstable + partner):.2e}")
    log(f"  (ε=0 predicts λ+λ'=0; ε>0 shifts the center)")


# ========================================================================
# SUMMARY
# ========================================================================
log()
log()
log("=" * 72)
log("SUMMARY")
log("=" * 72)
log()
log("Phase 3: Depolarizing noise simultaneously:")
log("  1. Breaks the palindrome (error ∝ ε, from X-dephasing component)")
log("  2. Shifts the eigenvalue center (Σγ no longer 0)")
log("  3. Modifies γ_crit")
log()

valid_gc = [(e, g) for e, g in zip(eps_values, gc_arr) if g is not None]
if len(valid_gc) >= 2:
    min_gc = min(g for _, g in valid_gc)
    max_gc = max(g for _, g in valid_gc)
    log(f"  γ_crit range: {min_gc:.6f} to {max_gc:.6f}")
    log(f"  γ_crit(ε=0) = {gc_ref:.6f}")

log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
