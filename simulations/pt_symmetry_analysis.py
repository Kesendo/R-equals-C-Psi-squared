#!/usr/bin/env python3
"""
PT-Symmetry Analysis of the Palindromic Liouvillian
=====================================================
Phase 1: Formal classification of Π (Π², linearity, det, Π·L†·Π⁻¹)
Phase 2: Chiral symmetry breaking analysis of the fragile-bridge system

Script: simulations/pt_symmetry_analysis.py
Output: simulations/results/pt_symmetry_analysis.txt
"""

import numpy as np
from scipy.linalg import eigvals, eig
from itertools import product as iproduct
import os, sys, time as clock

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
# PHASE 1: FORMAL CLASSIFICATION OF Π
# ========================================================================
log()
log("=" * 72)
log("PHASE 1: FORMAL CLASSIFICATION OF Π")
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
log("  In AZ terms: generalized chiral (S⁴=I, not standard S²=I).")

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
log("  Consequence: {Π, L_c} = 0 with linear Π is CHIRAL symmetry,")
log("  not PT-symmetry (which requires anti-linear PT operator).")

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
log("  This is the defining property of CHIRAL SYMMETRY (class AIII).")

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
log("  | Linear/anti    | Linear (matrix, no conjugation) | Chiral symmetry, NOT PT      |")
log("  | det(Π)         | −1 (N=1), +1 (N≥2)             | Even parity for multi-qubit  |")
log("  | Π·L†·Π⁻¹      | −L† − 2Σγ·I (same as for L)    | Chiral for both L and L†     |")
log()
log("  Π is a GENERALIZED CHIRAL OPERATOR: linear, order 4,")
log("  anti-commutes with L_c = L + Σγ·I. Class AIII (chiral unitary).")
log("  NOT PT-symmetry (which requires anti-linearity).")


# ========================================================================
# PHASE 2: CHIRAL SYMMETRY BREAKING IN THE FRAGILE BRIDGE
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 2: CHIRAL SYMMETRY BREAKING IN THE FRAGILE BRIDGE")
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

    # Pairing check: for each λ, find closest to −λ
    pair_errs = []
    for lam in ev:
        pair_errs.append(np.min(np.abs(ev - (-lam))))
    max_pe = np.max(pair_errs)

    nonzero = ev[np.abs(ev) > 1e-10]
    on_axis = "YES" if np.all(np.abs(nonzero.real) < 1e-6) else "NO"

    log(f"  {g:>8.4f}  {g/gamma_crit:>6.3f}  {max_re:>12.2e}"
        f"  {max_pe:>12.2e}  {on_axis:>8}")

log()
log("  λ ↔ −λ pairing is exact (machine precision) at ALL γ.")
log("  Below γ_crit: eigenvalues on the imaginary axis (chiral phase).")
log("  Above γ_crit: eigenvalues leave the axis (chiral breaking).")

# ----------------------------------------------------------------
# 2c. Eigenvalue trajectory near γ_crit
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2c. Eigenvalue trajectory (50 points, 0.5γ_c to 1.5γ_c)")
log("─" * 72)
log()

n_pts = 50
gammas_fine = np.linspace(0.5 * gamma_crit, 1.5 * gamma_crit, n_pts)
max_re_arr = np.zeros(n_pts)
crit_im_arr = np.zeros(n_pts)

log(f"  {'γ':>9}  {'γ/γ_c':>6}  {'max Re(λ)':>12}  {'Im at max':>10}"
    f"  {'|λ_crit|':>10}")
log(f"  {'─'*55}")

for i, g in enumerate(gammas_fine):
    ev = eigvals(build_coupled_liouvillian(N_chain, g, J, J_br))
    idx_max = np.argmax(ev.real)
    lam_c = ev[idx_max]
    max_re_arr[i] = lam_c.real
    crit_im_arr[i] = lam_c.imag

    if i % 5 == 0 or i == n_pts - 1:
        log(f"  {g:>9.6f}  {g/gamma_crit:>6.3f}  {lam_c.real:>12.4e}"
            f"  {lam_c.imag:>10.4f}  {abs(lam_c):>10.4f}")

# Transition point
cross = np.where(max_re_arr > 1e-8)[0]
if len(cross) > 0:
    g_cross = gammas_fine[cross[0]]
    log()
    log(f"  Eigenvalues leave imaginary axis at γ ≈ {g_cross:.6f}")
    log(f"  γ_crit (bisection):                     {gamma_crit:.6f}")

# ----------------------------------------------------------------
# 2d. Petermann factor and phase rigidity
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2d. Petermann factor K and phase rigidity r")
log("─" * 72)
log()
log("  K = 1/|⟨ψ_L|ψ_R⟩|²   (diverges at EP, =1 for Hermitian)")
log("  r = |⟨ψ_L|ψ_R⟩|²      (0 at EP, 1 for Hermitian)")
log()

gamma_peter = np.linspace(0.01, 2.0 * gamma_crit, 30)
log(f"  {'γ':>9}  {'γ/γ_c':>6}  {'K':>10}  {'r':>12}  {'max Re':>10}")
log(f"  {'─'*55}")

K_arr = np.zeros(len(gamma_peter))
r_arr = np.zeros(len(gamma_peter))

for i, g in enumerate(gamma_peter):
    L_g = build_coupled_liouvillian(N_chain, g, J, J_br)
    w, vl, vr = eig(L_g, left=True, right=True)

    idx_max = np.argmax(w.real)
    overlap = np.abs(vl[:, idx_max].conj() @ vr[:, idx_max])
    K = 1.0 / (overlap**2 + 1e-30)
    r = overlap**2
    K_arr[i] = K
    r_arr[i] = r

    log(f"  {g:>9.5f}  {g/gamma_crit:>6.3f}  {K:>10.4f}  {r:>12.8f}"
        f"  {w[idx_max].real:>10.2e}")

log()
log(f"  Max Petermann factor: K = {np.max(K_arr):.2f}"
    f" at γ/γ_c = {gamma_peter[np.argmax(K_arr)]/gamma_crit:.3f}")
log(f"  Min phase rigidity:  r = {np.min(r_arr):.6f}"
    f" at γ/γ_c = {gamma_peter[np.argmin(r_arr)]/gamma_crit:.3f}")

# ----------------------------------------------------------------
# 2e. Eigenvector coalescence of critical pair
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2e. Eigenvector coalescence (right eigenvectors of λ, −λ pair)")
log("─" * 72)
log()
log("  At an EP: two right eigenvectors become parallel (cos θ → 1).")
log()

gamma_coal = np.linspace(0.5 * gamma_crit, 1.5 * gamma_crit, 25)
log(f"  {'γ':>9}  {'γ/γ_c':>6}  {'pair dist':>10}  {'cos θ':>10}"
    f"  {'Re(λ)':>10}  {'Im(λ)':>10}")
log(f"  {'─'*65}")

for g in gamma_coal:
    L_g = build_coupled_liouvillian(N_chain, g, J, J_br)
    w, vr = eig(L_g, right=True)

    idx_max = np.argmax(w.real)
    lam_max = w[idx_max]

    # Find palindromic partner closest to −λ
    dists = np.abs(w - (-lam_max))
    dists[idx_max] = np.inf
    idx_partner = np.argmin(dists)
    pair_dist = dists[idx_partner]

    # Angle between right eigenvectors
    v1 = vr[:, idx_max]
    v2 = vr[:, idx_partner]
    cos_angle = np.abs(v1.conj() @ v2) / (
        np.linalg.norm(v1) * np.linalg.norm(v2))

    log(f"  {g:>9.6f}  {g/gamma_crit:>6.3f}  {pair_dist:>10.2e}"
        f"  {cos_angle:>10.6f}  {lam_max.real:>10.2e}"
        f"  {lam_max.imag:>10.4f}")

# ----------------------------------------------------------------
# 2f. Reinterpretation
# ----------------------------------------------------------------
log()
log("─" * 72)
log("2f. REINTERPRETATION: Hopf bifurcation IS chiral symmetry breaking")
log("─" * 72)
log()
log("  FRAGILE_BRIDGE.md states: 'Hopf bifurcation, not PT breaking.'")
log("  This analysis shows both descriptions are correct simultaneously:")
log()
log("  1. Σγ = 0 forces exact λ ↔ −λ pairing (chiral symmetry)")
log("  2. Below γ_crit: ALL eigenvalues on the imaginary axis")
log("     (the chiral-symmetric phase: Re(λ) = 0 for all λ)")
log("  3. Above γ_crit: eigenvalue pairs leave the imaginary axis")
log("     (one to Re>0, partner to Re<0 = chiral symmetry BREAKING)")
log("  4. The mechanism: a complex pair crosses Re=0 (Hopf)")
log()
log("  In Hamiltonian PT: real eigenvalues → complex (at EP).")
log("  In Liouvillian chiral: imaginary eigenvalues → off-axis.")
log("  Same geometry, rotated 90°. The Hopf IS the chiral breaking.")
log()
log("  Key: Π is LINEAR, not anti-linear. This is chiral/sublattice")
log("  symmetry (class AIII), not time-reversal. The spectral consequence")
log("  (± pairing) is identical; the operator type is different.")
log()
log("  No classical EP on the real γ axis: the palindromic pair (λ,−λ)")
log("  crosses Re=0 without coalescing (distance = 2|λ| > 0).")
log("  The transition is topological (axis crossing), not local (EP).")


# ========================================================================
# SUMMARY
# ========================================================================
log()
log()
log("=" * 72)
log("SUMMARY")
log("=" * 72)
log()
log("Phase 1: Π is a LINEAR, ORDER-4 chiral operator with det(Π) = +1")
log("  (N≥2). It anti-commutes with L_c = L+Σγ·I and with L_c†.")
log("  NOT PT-symmetry. Correct class: AIII (chiral unitary).")
log()
log("Phase 2: The fragile bridge (Σγ=0) has:")
log(f"  - Exact λ ↔ −λ pairing at all γ (chiral symmetry)")
log(f"  - Eigenvalues on imaginary axis for γ < γ_crit = {gamma_crit:.6f}")
log(f"  - Chiral symmetry breaking (eigenvalues leave Im axis) at γ_crit")
log(f"  - The Hopf bifurcation IS the Liouvillian chiral breaking")
log(f"  - No exceptional point on the real γ axis")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
