#!/usr/bin/env python3
"""
PT-Symmetry Analysis of the Palindromic Liouvillian
=====================================================
Phase 1: Shifted-generator structure of Π (Π², linearity, det, Π·L†·Π⁻¹)
Phase 2: Palindrome-axis departure in the fragile-bridge system

Script: simulations/pt_symmetry_analysis.py
Output: simulations/results/pt_symmetry_analysis.txt
"""

import numpy as np
from scipy.linalg import eigvals
from itertools import product as iproduct
import os, sys, time as clock
from pt_multiset_matching import multiset_reflection_error

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "pt_symmetry_analysis.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ========================================================================
# Pauli infrastructure
# ========================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]

# Π per-site map: I↔X (+1), Y↔Z (×i)
PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def build_pi(N):
    """Build Π as matrix in 4^N Pauli basis."""
    num = 4**N
    all_idx = list(iproduct(range(4), repeat=N))
    Pi = np.zeros((num, num), dtype=complex)
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        phase = 1
        for i in idx_b:
            phase *= PI_SIGN[i]
        a = all_idx.index(mapped)
        Pi[a, b] = phase
    return Pi, all_idx


def yz_weight(indices):
    """Count of Y(2) or Z(3) in a Pauli index tuple."""
    return sum(1 for i in indices if i in (2, 3))


def build_hamiltonian(N, bonds, J=1.0, delta=1.0):
    """Heisenberg/XXZ Hamiltonian, 2^N × 2^N."""
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    ops = [sx, sy, sz]
    for (i, j) in bonds:
        for pidx, P in enumerate(ops):
            coeff = J * (delta if pidx == 2 else 1.0)
            opi = np.eye(1, dtype=complex)
            for k in range(N):
                opi = np.kron(opi, P if k == i else I2)
            opj = np.eye(1, dtype=complex)
            for k in range(N):
                opj = np.kron(opj, P if k == j else I2)
            H += coeff * (opi @ opj)
    return H


def build_liouvillian_pauli(N, H, gamma_per_site):
    """Liouvillian in Pauli basis. Returns L_H, L_D, L."""
    dim = 2**N
    num = 4**N
    all_idx = list(iproduct(range(4), repeat=N))

    # Pre-build all Pauli matrices and stack
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in idx[1:]:
            m = np.kron(m, PAULIS[k])
        pmats.append(m)
    pstack = np.array(pmats)  # (num, dim, dim)

    # L_H via vectorized trace: L_H[a,b] = Tr(σ_a · (-i[H, σ_b])) / dim
    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        # Tr(σ_a · comm) = sum over ij of σ_a[i,j] * comm[j,i]
        L_H[:, b] = np.einsum('aij,ji->a', pstack, comm) / dim

    # L_D: Z-dephasing (diagonal in Pauli basis)
    L_D = np.zeros((num, num), dtype=complex)
    for a, idx in enumerate(all_idx):
        rate = 0.0
        for site in range(N):
            if idx[site] in (1, 2):  # X or Y at this site
                rate += 2 * gamma_per_site[site]
        L_D[a, a] = -rate

    return L_H, L_D, L_H + L_D


def op_at(op, qubit, n_qubits):
    """Place single-qubit operator on qubit k in n-qubit system."""
    result = np.array([[1]], dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == qubit else I2)
    return result


def build_coupled_liouvillian(n_per_chain, gamma, J=1.0, J_bridge=1.0):
    """Fragile bridge: two N-qubit chains, +γ (decay) and −γ (gain)."""
    n_total = 2 * n_per_chain
    d = 2**n_total
    d2 = d * d

    H = np.zeros((d, d), dtype=complex)
    # Internal bonds
    for chain_start in [0, n_per_chain]:
        for i in range(chain_start, chain_start + n_per_chain - 1):
            for P in [sx, sy, sz]:
                H += J * op_at(P, i, n_total) @ op_at(P, i + 1, n_total)
    # Bridge bond
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

    return L


def find_gamma_crit(n_per_chain, J_bridge, J=1.0, tol=1e-7):
    """Bisect for γ where max Re(λ) first exceeds 0."""
    def max_re(gamma):
        L = build_coupled_liouvillian(n_per_chain, gamma, J, J_bridge)
        return float(np.max(eigvals(L).real))

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


# ========================================================================
log("=" * 72)
log("PT-SYMMETRY ANALYSIS OF THE PALINDROMIC LIOUVILLIAN")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)


# ========================================================================
# PHASE 1: SHIFTED-GENERATOR STRUCTURE OF Π
# ========================================================================
log()
log("=" * 72)
log("PHASE 1: SHIFTED-GENERATOR STRUCTURE OF Π")
log("=" * 72)

# ----------------------------------------------------------------
# 1a. Π²
# ----------------------------------------------------------------
log()
log("─" * 72)
log("1a. Π² computation")
log("─" * 72)
log()
log("Analytical prediction (per-site):")
log("  Π(I)=X, Π(X)=I  → Π²(I)=I (+1),  Π²(X)=X (+1)")
log("  Π(Y)=iZ, Π(Z)=iY → Π²(Y)=−Y (i²), Π²(Z)=−Z (i²)")
log("  N-qubit: Π² = (−1)^{w_YZ}  where w_YZ = count of Y,Z entries")
log()

for N in [2, 3]:
    Pi, all_idx = build_pi(N)
    Pi2 = Pi @ Pi
    Pi4 = Pi2 @ Pi2
    num = 4**N

    diag_Pi2 = np.diag(Pi2)
    off_diag_max = np.max(np.abs(Pi2 - np.diag(diag_Pi2)))

    predicted = np.array([(-1.0)**yz_weight(idx) for idx in all_idx],
                         dtype=complex)
    pred_err = np.max(np.abs(diag_Pi2 - predicted))
    err_Pi4 = np.max(np.abs(Pi4 - np.eye(num)))

    n_plus = int(np.sum(np.abs(diag_Pi2 - 1) < 1e-10))
    n_minus = int(np.sum(np.abs(diag_Pi2 + 1) < 1e-10))

    log(f"  N={N} ({num}×{num}):")
    log(f"    Π² diagonal?      off-diag max = {off_diag_max:.2e}")
    log(f"    Π² = (−1)^w_YZ?   prediction error = {pred_err:.2e}")
    log(f"    Eigenvalues:       +1 ({n_plus})  −1 ({n_minus})")
    log(f"    Π⁴ = I?            error = {err_Pi4:.2e}")
    log()

log("  Result: Π² = (−1)^{w_YZ}. NOT ±I. Π has order 4.")
log("  Resolve Π² first: phase normalization gives an involutive P on each character sector.")

# ----------------------------------------------------------------
# 1b. Linearity
# ----------------------------------------------------------------
log()
log("─" * 72)
log("1b. Linear vs anti-linear")
log("─" * 72)
log()
log("  Π is a MATRIX acting on coefficient vectors in Pauli basis.")
log("  Π(αv) = α·Π(v) for all complex α. No conjugation involved.")
log("  Therefore: Π is LINEAR.")
log()
log("  Consequence: {Π, L_c} = 0 is a linear shifted-generator anticommutation.")
log("  Sectorwise P follows after resolving Π²; the full irreducible class stays open.")
log("  This is not PT-symmetry, which requires an anti-linear PT operator.")

# ----------------------------------------------------------------
# 1c. det(Π)
# ----------------------------------------------------------------
log()
log("─" * 72)
log("1c. det(Π)")
log("─" * 72)
log()
log("  Prediction: det(Π₁) = −1 for single site (4×4).")
log("  det(Π₁^{⊗N}) = det(Π₁)^{N·4^{N−1}} = (−1)^{N·4^{N−1}}")
log("  Since 4^{N−1} is even for N≥2: det = +1 for N≥2.")
log()

for N in [1, 2, 3, 4]:
    Pi, _ = build_pi(N)
    det_val = np.linalg.det(Pi)
    exp = N * 4**(N - 1)
    predicted = (-1)**exp
    log(f"  N={N}: det(Π) = {det_val.real:+.4f}{det_val.imag:+.4f}i"
        f"   predicted: {predicted:+d}   exponent: {exp}")

# ----------------------------------------------------------------
# 1d. Π·L†·Π⁻¹
# ----------------------------------------------------------------
log()
log("─" * 72)
log("1d. Π·L†·Π⁻¹ for N=3 Heisenberg chain")
log("─" * 72)
log()
log("  Analytical derivation:")
log("    L† = −L_H + L_D  (L_H anti-Hermitian, L_D real diagonal)")
log("    Π·L†·Π⁻¹ = Π(−L_H+L_D)Π⁻¹ = L_H + (−L_D − 2Σγ·I)")
log("             = −(−L_H+L_D) − 2Σγ·I = −L† − 2Σγ·I")
log()

N = 3
bonds = [(i, i + 1) for i in range(N - 1)]
gammas = [0.05] * N
H = build_hamiltonian(N, bonds)
t0 = clock.time()
L_H, L_D, L = build_liouvillian_pauli(N, H, gammas)
Pi, _ = build_pi(N)
Pi_inv = np.linalg.inv(Pi)
c = 2 * sum(gammas)
num = 4**N

L_dag = L.conj().T
err_L = np.max(np.abs(Pi @ L @ Pi_inv + L + c * np.eye(num)))
err_Ld = np.max(np.abs(Pi @ L_dag @ Pi_inv + L_dag + c * np.eye(num)))

log(f"  N=3, γ=[0.05]*3, Σγ={sum(gammas):.2f}, c=2Σγ={c:.2f}")
log(f"  ‖Π·L·Π⁻¹ + L + c·I‖   = {err_L:.2e}   (palindrome)")
log(f"  ‖Π·L†·Π⁻¹ + L† + c·I‖ = {err_Ld:.2e}   (adjoint)")
log(f"  ({clock.time()-t0:.1f}s)")
log()
log("  Both hold to machine precision.")
log("  For L_c = L+Σγ·I: Π·L_c·Π⁻¹ = −L_c AND Π·L_c†·Π⁻¹ = −L_c†.")

# ----------------------------------------------------------------
# 1e. Fragile bridge at Σγ = 0
# ----------------------------------------------------------------
log()
log("─" * 72)
log("1e. Fragile bridge palindrome (4 qubits, Σγ = 0)")
log("─" * 72)
log()

N_fb = 4
bonds_fb = [(0, 1), (1, 2), (2, 3)]
gammas_fb = [0.1, 0.1, -0.1, -0.1]
H_fb = build_hamiltonian(N_fb, bonds_fb)
t0 = clock.time()
_, _, L_fb = build_liouvillian_pauli(N_fb, H_fb, gammas_fb)
Pi_fb, _ = build_pi(N_fb)
Pi_fb_inv = np.linalg.inv(Pi_fb)
c_fb = 2 * sum(gammas_fb)
num_fb = 4**N_fb

err_fb = np.max(np.abs(Pi_fb @ L_fb @ Pi_fb_inv + L_fb
                       + c_fb * np.eye(num_fb)))

log(f"  γ = [+0.1, +0.1, −0.1, −0.1],  Σγ = {sum(gammas_fb):.1f},  c = {c_fb:.1f}")
log(f"  ‖Π·L·Π⁻¹ + L‖ = {err_fb:.2e}")
log(f"  ({clock.time()-t0:.1f}s)")
log()
log("  Exact at Σγ=0: every eigenvalue λ pairs with −λ.")
log("  This is the shifted-generator anticommutation carried by Π.")
log("  After resolving a Π² character p_x, P_px=sqrt(p_x)·Π is an involution: sectorwise P.")

# ----------------------------------------------------------------
# Phase 1 summary table
# ----------------------------------------------------------------
log()
log("─" * 72)
log("PHASE 1 SUMMARY TABLE")
log("─" * 72)
log()
log("  | Property       | Value                          | Consequence                  |")
log("  |----------------|--------------------------------|------------------------------|")
log("  | Π²             | (−1)^{w_YZ}, diagonal parity   | Not involution; order 4      |")
log("  | Linear/anti    | Linear (matrix, no conjugation) | Sectorwise P, NOT PT          |")
log("  | det(Π)         | −1 (N=1), +1 (N≥2)             | Even parity for multi-qubit  |")
log("  | Π·L†·Π⁻¹      | −L† − 2Σγ·I (same as for L)    | Shifted anticommutation      |")
log()
log("  Π is linear and order 4; it anti-commutes with L_c = L + Σγ·I.")
log("  On each resolved Π²-character sector, P_px=sqrt(p_x)·Π is involutive")
log("  and anti-commutes with L_c: sectorwise P is secured.")
log("  The full irreducible SRP class remains OPEN; this is not a global class assignment.")
log("  It is NOT PT-symmetry, which requires anti-linearity.")


# ========================================================================
# PHASE 2: PALINDROME AXIS DEPARTURE IN THE FRAGILE BRIDGE
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 2: PALINDROME AXIS DEPARTURE IN THE FRAGILE BRIDGE")
log("=" * 72)

N_chain = 2
J = 1.0
J_br = 1.0

# ----------------------------------------------------------------
# 2a. Find γ_crit
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2a. Finding γ_crit")
log("─" * 72)
log()
log(f"  N=2 per chain (4 qubits), J={J}, J_bridge={J_br}")

t0 = clock.time()
gamma_crit = find_gamma_crit(N_chain, J_br, J=J)
log(f"  γ_crit = {gamma_crit:.7f}  ({clock.time()-t0:.1f}s)")

# ----------------------------------------------------------------
# 2b. Eigenvalue pairing λ ↔ −λ
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2b. Eigenvalue pairing λ ↔ −λ")
log("─" * 72)
log()

test_gammas = [0.01, 0.05, 0.10, gamma_crit * 0.5,
               gamma_crit * 0.9, gamma_crit * 1.1, 0.5]
log(f"  {'γ':>8}  {'γ/γ_c':>6}  {'max|Re|':>12}  {'pair err':>12}  {'Im axis?':>8}")
log(f"  {'─'*52}")

for g in test_gammas:
    L_g = build_coupled_liouvillian(N_chain, g, J, J_br)
    ev = eigvals(L_g)
    max_re = np.max(np.abs(ev.real))

    # Bijective bottleneck matching preserves algebraic multiplicity.  A
    # nearest-neighbour loop can reuse one target and miss a multiplicity error.
    max_pe = multiset_reflection_error(ev, midpoint=0.0)

    nonzero = ev[np.abs(ev) > 1e-10]
    on_axis = "YES" if np.all(np.abs(nonzero.real) < 1e-6) else "NO"

    log(f"  {g:>8.4f}  {g/gamma_crit:>6.3f}  {max_re:>12.2e}"
        f"  {max_pe:>12.2e}  {on_axis:>8}")

log()
log("  λ ↔ −λ pairing is exact (machine precision) at ALL γ.")
log("  Below γ_crit: eigenvalues remain within palindrome axis.")
log("  Above γ_crit: eigenvalues move off palindrome axis; operator relation remains exact.")

# ----------------------------------------------------------------
# 2c. Spectral-abscissa scan near γ_crit
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2c. Spectral abscissa (50 independent spectra, 0.5γ_c to 1.5γ_c)")
log("─" * 72)
log()

n_pts = 50
gammas_fine = np.linspace(0.5 * gamma_crit, 1.5 * gamma_crit, n_pts)
max_re_arr = np.zeros(n_pts)
log(f"  {'γ':>9}  {'γ/γ_c':>6}  {'max Re(λ)':>12}")
log(f"  {'─'*32}")

for i, g in enumerate(gammas_fine):
    ev = eigvals(build_coupled_liouvillian(N_chain, g, J, J_br))
    idx_max = np.argmax(ev.real)
    lam_c = ev[idx_max]
    max_re_arr[i] = lam_c.real
    if i % 5 == 0 or i == n_pts - 1:
        log(f"  {g:>9.6f}  {g/gamma_crit:>6.3f}  {lam_c.real:>12.4e}")

# Transition point
cross = np.where(max_re_arr > 1e-8)[0]
if len(cross) > 0:
    g_cross = gammas_fine[cross[0]]
    log()
    log(f"  Eigenvalues leave imaginary axis at γ ≈ {g_cross:.6f}")
    log(f"  γ_crit (bisection):                     {gamma_crit:.6f}")

# ----------------------------------------------------------------
# 2d. Basis-invariant evidence boundary
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2d. BASIS-INVARIANT EVIDENCE BOUNDARY")
log("─" * 72)
log()
log("  This producer does not report single-eigenvector Petermann factors or angles.")
log("  Degenerate eigenspaces make those values basis-dependent; a subspace-level")
log("  condition measure or a threshold Jordan-rank test would be required.")

# ----------------------------------------------------------------
# 2e. Reinterpretation
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2e. PALINDROME-AXIS DEPARTURE")
log("─" * 72)
log()
log("  The measured transition is an axis departure while the palindrome remains exact:")
log()
log("  1. Σγ = 0 forces exact λ ↔ −λ palindrome pairing")
log("  2. Below γ_crit: ALL eigenvalues on the imaginary axis")
log("     (the axis-confined phase: Re(λ) = 0 for all λ)")
log("  3. Above γ_crit: off-axis eigenvalue quartets appear")
log("     (Π gives λ↔−λ; Hermiticity preservation gives λ↔λ*; together {λ,λ*,−λ,−λ*})")
log("  4. The measured event is the existence of an off-axis pair with max Re>0")
log()
log("  In Hamiltonian PT: real eigenvalues → complex (at EP).")
log("  In this Liouvillian palindrome: imaginary eigenvalues → symmetric off-axis pairs.")
log("  The spectrum changes axis occupancy without breaking the palindrome relation.")
log()
log("  Key: Π is LINEAR, not anti-linear. Resolving Π² and phase-normalizing")
log("  gives a sectorwise P generator, not a global irreducible class assignment.")
log("  The spectral consequence is the ± pairing; the full SRP class stays open.")
log()
log("  This scan independently selects max Re(λ) at each γ; it does not track a branch.")
log("  It reports the spectral-abscissa axis departure only.")


# ========================================================================
# SUMMARY
# ========================================================================
log()
log()
log("=" * 72)
log("SUMMARY")
log("=" * 72)
log()
log("Phase 1: Π is a LINEAR, ORDER-4 palindrome operator with det(Π) = +1")
log("  (N≥2). It anti-commutes with L_c = L+Σγ·I and with L_c†.")
log("  Resolving Π² gives sectorwise P; full irreducible SRP class remains OPEN.")
log("  This does not assign a global class and is not PT-symmetry.")
log()
log("Phase 2: The fragile bridge (Σγ=0) has:")
log(f"  - Exact λ ↔ −λ palindrome pairing at all γ")
log(f"  - Eigenvalues on imaginary axis for γ < γ_crit = {gamma_crit:.6f}")
log(f"  - Palindrome-axis departure (eigenvalues leave Im axis) at γ_crit")
log(f"  - The spectral-abscissa axis departure occurs while the operator relation remains exact")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
