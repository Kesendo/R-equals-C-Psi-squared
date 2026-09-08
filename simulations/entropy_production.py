#!/usr/bin/env python3
"""
Algebraic and Dynamical Diagnostics of the Palindromic Liouvillian
===================================================================
Phase 1: Decay-rate pairing and state-dependent entropy trajectories
Phase 2: Descriptive exponential transforms of paired rate differences
Phase 3: Descriptive rate-ratio regression
Phase 4: Formal balanced gain-loss generator diagnostics
Phase 5: CΨ and sampled occupation-number variance

Script: simulations/entropy_production.py
Output: simulations/results/entropy_production.txt
"""

import numpy as np
from scipy.linalg import eigvals, eig, expm, logm
from scipy.optimize import brentq
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "entropy_production.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ========================================================================
# Infrastructure
# ========================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_list(mats):
    r = mats[0]
    for m in mats[1:]:
        r = np.kron(r, m)
    return r


def op_n(pauli, site, N):
    ops = [I2] * N
    ops[site] = pauli
    return kron_list(ops)


def build_heisenberg_liouvillian(N, J=1.0, gamma=0.05):
    """Heisenberg chain with Z-dephasing."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * op_n(P, i, N) @ op_n(P, i + 1, N)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = op_n(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d**2))
    return L


def build_coupled_liouvillian(n_per_chain, gamma, J=1.0, J_bridge=1.0):
    """Fragile bridge: two chains, +γ decay and -γ gain."""
    n_total = 2 * n_per_chain
    d = 2**n_total
    d2 = d * d
    H = np.zeros((d, d), dtype=complex)
    for cs in [0, n_per_chain]:
        for i in range(cs, cs + n_per_chain - 1):
            for P in [sx, sy, sz]:
                H += J * op_n(P, i, n_total) @ op_n(P, i + 1, n_total)
    for P in [sx, sy, sz]:
        H += J_bridge * op_n(P, n_per_chain - 1, n_total) @ \
             op_n(P, n_per_chain, n_total)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(n_per_chain):
        Zk = op_n(sz, k, n_total)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    for k in range(n_per_chain, 2 * n_per_chain):
        Zk = op_n(sz, k, n_total)
        L += (-gamma) * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve_rho(L, rho0, t):
    """Evolve density matrix to time t."""
    d2 = L.shape[0]
    d = int(np.sqrt(d2))
    rho_vec = expm(L * t) @ rho0.flatten()
    return rho_vec.reshape(d, d)


def von_neumann_entropy(rho):
    """S(ρ) = -Tr(ρ ln ρ)."""
    eigv = np.linalg.eigvalsh(rho)
    eigv = eigv[eigv > 1e-30]
    return -np.sum(eigv * np.log(eigv))


def relative_entropy(rho, sigma):
    """D(ρ || σ) = Tr(ρ(ln ρ - ln σ))."""
    eigr, Ur = np.linalg.eigh(rho)
    eigs, Us = np.linalg.eigh(sigma)
    eigr = np.maximum(eigr, 1e-30)
    eigs = np.maximum(eigs, 1e-30)
    ln_rho = Ur @ np.diag(np.log(eigr)) @ Ur.conj().T
    ln_sigma = Us @ np.diag(np.log(eigs)) @ Us.conj().T
    return np.real(np.trace(rho @ (ln_rho - ln_sigma)))


def von_neumann_entropy_rate(L, rho0, t):
    """State-dependent dS/dt = -Tr(L[rho(t)] ln rho(t))."""
    d2 = L.shape[0]
    d = int(np.sqrt(d2))
    rho = evolve_rho(L, rho0, t)
    # L[ρ]
    drho_vec = L @ rho.flatten()
    drho = drho_vec.reshape(d, d)
    # ln(ρ)
    eigv, U = np.linalg.eigh(rho)
    eigv = np.maximum(eigv, 1e-30)
    ln_rho = U @ np.diag(np.log(eigv)) @ U.conj().T
    return -np.real(np.trace(drho @ ln_rho))


def compute_cpsi(rho):
    d = rho.shape[0]
    purity = np.real(np.trace(rho @ rho))
    L1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    return purity * L1 / (d - 1)


def palindromic_pairs(eigenvalues, sigma_gamma):
    """Find palindromic pairs: λ + λ' ≈ -2Σγ."""
    pairs = []
    used = set()
    target_sum = -2 * sigma_gamma
    for i in range(len(eigenvalues)):
        if i in used:
            continue
        best_j, best_d = -1, np.inf
        for j in range(len(eigenvalues)):
            if j != i and j not in used:
                d = abs((eigenvalues[i] + eigenvalues[j]) - target_sum)
                if d < best_d:
                    best_j = j
                    best_d = d
        if best_j >= 0 and best_d < 0.01:
            pairs.append((i, best_j))
            used.update([i, best_j])
    return pairs


# ========================================================================
log("=" * 72)
log("ALGEBRAIC AND DYNAMICAL DIAGNOSTICS OF THE PALINDROMIC LIOUVILLIAN")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)


# ========================================================================
# PHASE 1: DECAY-RATE PAIRING AND STATE-DEPENDENT ENTROPY TRAJECTORIES
# ========================================================================
log()
log("=" * 72)
log("PHASE 1: DECAY-RATE PAIRING AND STATE-DEPENDENT ENTROPY TRAJECTORIES")
log("=" * 72)
log()

gamma = 0.05
J = 1.0

for N in [2, 3, 4]:
    t0 = clock.time()
    L = build_heisenberg_liouvillian(N, J, gamma)
    ev = eigvals(L)
    sigma_gamma = N * gamma
    d = 2**N

    pairs = palindromic_pairs(ev, sigma_gamma)

    # Decay rate pairing check
    rate_sums = []
    for i, j in pairs:
        d_i = -ev[i].real
        d_j = -ev[j].real
        rate_sums.append(d_i + d_j)

    log(f"  N={N} ({len(ev)} eigenvalues, {len(pairs)} pairs, "
        f"Σγ={sigma_gamma:.3f}):")
    log(f"    Rate pair sums: mean={np.mean(rate_sums):.6f},"
        f" std={np.std(rate_sums):.2e} (should be 2Σγ={2*sigma_gamma:.3f})")

    # State-dependent trajectory: Bell+ at N=2, |0...0> pointer/H eigenstate otherwise.
    rho0 = np.zeros((d, d), dtype=complex)
    if N == 2:
        psi = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
        rho0 = np.outer(psi, psi.conj())
    else:
        rho0[0, 0] = 1.0  # |0...0⟩

    rho_ref = np.eye(d, dtype=complex) / d

    n_t = 50
    times = np.linspace(0.01, 30, n_t)
    sigma_arr = np.zeros(n_t)
    S_arr = np.zeros(n_t)
    D_arr = np.zeros(n_t)
    P_arr = np.zeros(n_t)

    for ti, t in enumerate(times):
        rho_t = evolve_rho(L, rho0, t)
        sigma_arr[ti] = von_neumann_entropy_rate(L, rho0, t)
        S_arr[ti] = von_neumann_entropy(rho_t)
        D_arr[ti] = relative_entropy(rho_t, rho_ref)
        P_arr[ti] = np.real(np.trace(rho_t @ rho_t))

    # dD/dt should be -σ
    dD_dt = np.gradient(D_arr, times)

    log(f"    dS/dt at t=0.01 = {sigma_arr[0]:.6f}")
    log(f"    largest sampled dS/dt = {np.max(sigma_arr):.6f} at t = {times[np.argmax(sigma_arr)]:.2f} (50-point grid)")
    log(f"    dS/dt at t=30 = {sigma_arr[-1]:.6f}")
    log(f"    D(ρ||I/d), chosen I/d reference: {D_arr[0]:.4f} → {D_arr[-1]:.4f}")
    log(f"    Purity: {P_arr[0]:.4f} → {P_arr[-1]:.4f}"
        f" (1/d = {1/d:.4f})")
    log(f"    ({clock.time()-t0:.1f}s)")
    log()


# ========================================================================
# PHASE 2: EXPONENTIAL TRANSFORMS OF PAIRED DECAY-RATE DIFFERENCES
# ========================================================================
log()
log("=" * 72)
log("PHASE 2: EXPONENTIAL TRANSFORMS OF PAIRED DECAY-RATE DIFFERENCES")
log("=" * 72)
log()
log("  Descriptive averages over palindromic decay-rate pairs.")
log("  The algebraic decay-rate pair sum d_fast+d_slow=2Σγ is not a thermodynamic entropy-production scale.")
log()

for N in [2, 3, 4, 5]:
    L = build_heisenberg_liouvillian(N, J, gamma)
    ev = eigvals(L)
    sigma_gamma = N * gamma
    pairs = palindromic_pairs(ev, sigma_gamma)

    if len(pairs) == 0:
        continue

    delta_d = []
    for i, j in pairs:
        d_i = -ev[i].real
        d_j = -ev[j].real
        d_fast = max(d_i, d_j)
        d_slow = min(d_i, d_j)
        # Phase 2 intentionally retains every palindrome pair, including the
        # stationary/max-rate endpoint pairs with d_slow = 0.  Phase 3 below is
        # the separate positive-rate population used by the logarithmic identity.
        delta_d.append(d_fast - d_slow)

    if len(delta_d) == 0:
        continue

    delta_d = np.array(delta_d)

    # Descriptive exponential transforms; no fluctuation-theorem interpretation.
    exp_neg = np.mean(np.exp(-delta_d))
    exp_neg_norm = np.mean(np.exp(-delta_d / (2 * sigma_gamma)))

    log(f"  N={N} ({len(delta_d)} all palindrome pairs, Σγ={sigma_gamma:.3f}):")
    log(f"    Δd range: [{np.min(delta_d):.4f}, {np.max(delta_d):.4f}]")
    log(f"    ⟨Δd⟩ = {np.mean(delta_d):.4f}")
    log(f"    ⟨exp(-Δd)⟩ = {exp_neg:.6f}  (descriptive)")
    log(f"    ⟨exp(-Δd / 2Σγ)⟩ = {exp_neg_norm:.6f}  (normalized)")
    log()


# ========================================================================
# PHASE 3: DESCRIPTIVE RATE-RATIO REGRESSION
# ========================================================================
log()
log("=" * 72)
log("PHASE 3: DESCRIPTIVE RATE-RATIO REGRESSION")
log("=" * 72)
log()
log("  For each pair: ln(d_fast / d_slow) vs (d_fast - d_slow)")
log("  This is a numerical regression of decay-rate pairs, not a temperature fit.")
log()

for N in [3, 4, 5]:
    L = build_heisenberg_liouvillian(N, J, gamma)
    ev = eigvals(L)
    sigma_gamma = N * gamma
    pairs = palindromic_pairs(ev, sigma_gamma)

    log_ratios = []
    deltas = []

    for i, j in pairs:
        d_i = -ev[i].real
        d_j = -ev[j].real
        d_fast = max(d_i, d_j)
        d_slow = min(d_i, d_j)
        if d_fast > 0.01 and d_slow > 0.01:
            log_ratios.append(np.log(d_fast / d_slow))
            deltas.append(d_fast - d_slow)

    if len(log_ratios) < 3:
        continue

    log_ratios = np.array(log_ratios)
    deltas = np.array(deltas)

    # Linear fit: ln(ratio) = a * delta + b
    if np.std(deltas) > 1e-10:
        a, b = np.polyfit(deltas, log_ratios, 1)
        residual = np.std(log_ratios - (a * deltas + b))
    else:
        a, b, residual = 0, 0, 0

    log(f"  N={N} ({len(log_ratios)} pairs):")
    log(f"    Linear fit: ln(d_fast/d_slow) = {a:.4f}·Δd + {b:.4f}")
    log(f"    Residual: {residual:.4f}")
    if residual < 0.1 * np.std(log_ratios):
        log(f"    residual is small relative to this sample's log-ratio spread")
    else:
        log(f"    residual is not small relative to this sample's log-ratio spread")
    log()


# ========================================================================
# PHASE 4: FORMAL BALANCED GAIN-LOSS GENERATOR
# ========================================================================
log()
log("=" * 72)
log("PHASE 4: FORMAL BALANCED GAIN-LOSS GENERATOR")
log("=" * 72)
log()

N_chain = 2
J_br = 1.0
gamma_fb = 0.10  # well below γ_crit ≈ 0.187

L_fb = build_coupled_liouvillian(N_chain, gamma_fb, J=1.0, J_bridge=J_br)
ev_fb = eigvals(L_fb)
d_fb = 2**(2 * N_chain)

# At Σγ = 0 the algebraic spectrum pairs as λ ↔ -λ. Negative dephasing rates
# make this a formal generator, not a completely positive physical channel.

# Initial state: four-qubit GHZ/cat coherence spanning both chains.
psi_fb = np.zeros(d_fb, dtype=complex)
psi_fb[0] = 1.0 / np.sqrt(2)  # |0000⟩
psi_fb[d_fb - 1] = 1.0 / np.sqrt(2)  # |1111⟩
rho0_fb = np.outer(psi_fb, psi_fb.conj())

times_fb = np.linspace(0.01, 20, 40)
S_fb = np.zeros(len(times_fb))
for ti, t in enumerate(times_fb):
    rho_t = evolve_rho(L_fb, rho0_fb, t)
    S_fb[ti] = von_neumann_entropy(rho_t)

log("  Initial state: four-qubit GHZ/cat coherence (|0000⟩+|1111⟩)/sqrt(2)")
log(f"  Fragile bridge: N=2/chain, γ=±{gamma_fb}, J_bridge={J_br}")
log(f"  Σγ = 0 (algebraic rate-profile sum)")
log("  Scope: formal gain-loss generator is not a physical Lindblad channel.")
log(f"  positive-eigenvalue entropy diagnostic at t=0.01: {S_fb[0]:.4f}")
log(f"  positive-eigenvalue entropy diagnostic at t≈10: {S_fb[len(times_fb)//2]:.4f}")
log(f"  positive-eigenvalue entropy diagnostic at t=20: {S_fb[-1]:.4f}")
log(f"  ln(d) reference = {np.log(d_fb):.4f}")
log()

max_re = np.max(ev_fb.real)
log(f"  Max Re(λ) = {max_re:.2e} ({'stable' if max_re < 1e-6 else 'UNSTABLE'})")
log(f"  Eigenvalue real-part range: [{np.min(ev_fb.real):.2e}, {np.max(ev_fb.real):.2e}]")
log("  No bath temperature, steady-state uniqueness, efficiency, or thermodynamic entropy-production conclusion is inferred.")


# ========================================================================
# PHASE 5: CΨ = 1/4 AND OCCUPATION VARIANCE
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 5: CΨ AND SAMPLED OCCUPATION VARIANCE")
log("=" * 72)
log()
log("  Compare the sampled CΨ=1/4 crossing time with the sampled variance maximum.")
log()

for N in [2, 3]:
    d = 2**N
    L = build_heisenberg_liouvillian(N, J, gamma)

    # Bell+ for N=2, product |+⟩^N for N>2
    if N == 2:
        psi = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    else:
        psi = np.ones(d, dtype=complex) / np.sqrt(d)
    rho0 = np.outer(psi, psi.conj())

    times_c = np.linspace(0.01, 50, 200)
    cpsi_arr = np.zeros(len(times_c))
    var_arr = np.zeros(len(times_c))
    purity_arr = np.zeros(len(times_c))

    for ti, t in enumerate(times_c):
        rho_t = evolve_rho(L, rho0, t)
        cpsi_arr[ti] = compute_cpsi(rho_t)
        purity_arr[ti] = np.real(np.trace(rho_t @ rho_t))

        # Occupation number variance: ⟨(n_k - ⟨n_k⟩)²⟩ summed over qubits
        total_var = 0
        for k in range(N):
            Zk = op_n(sz, k, N)
            nk = (np.eye(d) - Zk) / 2  # projector |1⟩⟨1| on qubit k
            exp_n = np.real(np.trace(nk @ rho_t))
            exp_n2 = np.real(np.trace(nk @ rho_t @ nk))
            total_var += exp_n2 - exp_n**2
        var_arr[ti] = total_var / N  # per-qubit average

    # Bracket and refine the CΨ = 1/4 crossing.  Reporting the left grid
    # endpoint as the crossing biases this coarse 200-point scan by one bin.
    cross_bracket = None
    for k in range(len(cpsi_arr) - 1):
        if cpsi_arr[k] > 0.25 and cpsi_arr[k + 1] <= 0.25:
            cross_bracket = (times_c[k], times_c[k + 1])
            break

    log(f"  N={N}:")
    if cross_bracket is not None:
        crossing = brentq(
            lambda t: compute_cpsi(evolve_rho(L, rho0, t)) - 0.25,
            *cross_bracket)
        rho_cross = evolve_rho(L, rho0, crossing)
        crossing_var = 0.0
        for k in range(N):
            Zk = op_n(sz, k, N)
            nk = (np.eye(d) - Zk) / 2
            exp_n = np.real(np.trace(nk @ rho_cross))
            exp_n2 = np.real(np.trace(nk @ rho_cross @ nk))
            crossing_var += exp_n2 - exp_n**2
        crossing_var /= N
        log(f"    sampled crossing bracket: [{cross_bracket[0]:.4f}, {cross_bracket[1]:.4f}]")
        log(f"    refined CΨ=1/4 crossing: t = {crossing:.6f}")
        log(f"    variance at refined crossing: {crossing_var:.6f}")
    else:
        log(f"    CΨ does not cross 1/4")

    variance_span = np.ptp(var_arr)
    log(f"    sampled variance range: [{np.min(var_arr):.6f}, {np.max(var_arr):.6f}]")
    if variance_span <= 1e-10:
        log("    sampled variance is constant; no peak time is defined")
    else:
        max_var_idx = np.argmax(var_arr)
        log(f"    largest sampled variance: {var_arr[max_var_idx]:.6f}"
            f" at t = {times_c[max_var_idx]:.2f}")
    log(f"    Theoretical max (f(1-f)): 0.250000")
    log()


# ========================================================================
# SUMMARY
# ========================================================================
log()
log("=" * 72)
log("SUMMARY")
log("=" * 72)
log()
log("Phase 1: The algebraic decay-rate pair sum d_k+d_k'=2Σγ is checked")
log("  separately from the state-dependent, grid-sampled dS/dt trajectory.")
log("  I/d is a chosen I/d reference, not a uniquely selected steady state.")
log()
log("Phase 2: Exponential transforms of paired rate differences are reported")
log("  descriptively; 2Σγ is not a thermodynamic entropy-production scale.")
log()
log("Phase 3: ln(d_fast/d_slow) vs Δd is a descriptive linear regression.")
log()
log("Phase 4: The formal gain-loss generator is not a physical Lindblad channel;")
log("  only its computed spectrum and sampled matrix-evolution diagnostics are reported.")
log()
log("Phase 5: CΨ = 1/4 and occupation variance.")
log("  Whether they coincide: the data speaks.")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
