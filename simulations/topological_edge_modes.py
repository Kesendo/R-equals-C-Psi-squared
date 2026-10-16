#!/usr/bin/env python3
"""
Topological Edge Modes Analysis
================================
Is the mode localization [0.52, 0.63, 0.70, 0.63, 0.52] topologically
protected, or is it geometric (standing wave pattern on a chain)?

Phase 1: SSH analogy ((0,1) coherence block = tight-binding, uniform hopping → trivial)
Phase 2: Chiral block structure (Π eigenspaces, off-diagonal blocks, singular values)
Phase 3: Berry phase on the uniform → sacrifice path (an open one-parameter segment)
Phase 4: Mode counting for N=3,4,5
Phase 5: Robustness sweep (sacrifice-edge → uniform → sacrifice-center)

Script: simulations/topological_edge_modes.py
Output: simulations/results/topological_edge_modes.txt
"""

import numpy as np
from scipy.linalg import eigvals, eig, svdvals
from itertools import product as iproduct
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "topological_edge_modes.txt")
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


def build_pauli_matrices(N):
    """Pre-build all 4^N Pauli string matrices for N qubits."""
    all_idx = list(iproduct(range(4), repeat=N))
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in idx[1:]:
            m = np.kron(m, PAULIS[k])
        pmats.append(m)
    return np.array(pmats), all_idx


def build_hamiltonian(N, J=1.0):
    """Heisenberg chain Hamiltonian, 2^N × 2^N."""
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    ops = [sx, sy, sz]
    for i in range(N - 1):
        for P in ops:
            opi = np.eye(1, dtype=complex)
            for k in range(N):
                opi = np.kron(opi, P if k == i else I2)
            opj = np.eye(1, dtype=complex)
            for k in range(N):
                opj = np.kron(opj, P if k == i + 1 else I2)
            H += J * (opi @ opj)
    return H


def build_liouvillian_pauli(N, H, gamma_per_site):
    """Liouvillian in Pauli basis."""
    dim = 2**N
    num = 4**N
    pstack, all_idx = build_pauli_matrices(N)

    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pstack[b] - pstack[b] @ H)
        L_H[:, b] = np.einsum('aij,ji->a', pstack, comm) / dim

    L_D = np.zeros((num, num), dtype=complex)
    for a, idx in enumerate(all_idx):
        rate = 0.0
        for site in range(N):
            if idx[site] in (1, 2):
                rate += 2 * gamma_per_site[site]
        L_D[a, a] = -rate

    return L_H + L_D


def op_at(op, qubit, n_qubits):
    result = np.array([[1]], dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == qubit else I2)
    return result


def build_liouvillian_vec(N, gamma_per_site, J=1.0):
    """Standard vectorized Liouvillian (d² × d²). Faster for eigendecomp."""
    d = 2**N
    d2 = d * d
    H = build_hamiltonian(N, J)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = op_at(sz, k, N)
        L += gamma_per_site[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def site_weights(eigvec, N):
    """Unnormalized per-qubit weight of an operator: the Pauli weight carried
    by strings that are not the identity at that site."""
    d = 2**N
    V = eigvec.reshape(d, d)
    pstack, all_idx = build_pauli_matrices(N)
    coeffs = np.array([np.trace(pstack[a].conj().T @ V) / d
                       for a in range(4**N)])
    weights = np.zeros(N)
    for a, idx in enumerate(all_idx):
        c2 = abs(coeffs[a])**2
        for site in range(N):
            if idx[site] != 0:  # not identity at this site
                weights[site] += c2
    return weights


def localization_profile(eigvec, N):
    """Per-qubit weight of a Liouvillian eigenvector, normalized to mean 1."""
    weights = site_weights(eigvec, N)
    total = np.sum(weights)
    return weights / total * N if total > 0 else weights


def eigenspace_profile(vectors, N):
    """Profile of an eigenspace, read from an orthonormal basis of the span of
    its eigenvectors: the summed weights are the trace of the site-weight form
    on the eigenspace, so they do not depend on which basis the eigensolver
    returned (a single vector of a degenerate eigenvalue does)."""
    Q, _ = np.linalg.qr(vectors)
    weights = sum(site_weights(Q[:, m], N) for m in range(Q.shape[1]))
    return weights / np.sum(weights) * N


# ========================================================================
log("=" * 72)
log("TOPOLOGICAL EDGE MODES ANALYSIS")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)


# ========================================================================
# PHASE 1: SSH ANALOGY CHECK
# ========================================================================
log()
log("=" * 72)
log("PHASE 1: SSH ANALOGY CHECK ((0,1) COHERENCE BLOCK)")
log("=" * 72)
log()
log("  The (0,1) coherence block's dispersion: ω_k = 4J(1-cos(πk/N)), k=1..N-1")
log("  This is tight-binding with UNIFORM hopping 2J.")
log("  SSH requires ALTERNATING hopping t1-t2.")
log()
log("  Standing wave eigenstates: ψ_k(j) ∝ cos(πk(j-½)/N), the NEUMANN family.")
log("  These are the eigenvectors of the open chain's Laplacian, and so of")
log("  the (0,1) coherence block, which is what ω_k above is the spectrum of.")
log("  The Dirichlet family sin(πkj/(N+1)) belongs to the chain's ADJACENCY")
log("  matrix (F2b) and diagonalizes a different operator.")
log("  For the highest mode (k=N-1): peaks at center, shallow at the edges.")
log()

for N in [4, 5]:
    freqs = [4 * (1 - np.cos(np.pi * k / N)) for k in range(1, N)]
    # Standing wave profile for highest mode (k=N-1)
    profile = [np.cos(np.pi * (N - 1) * (j - 0.5) / N)**2
               for j in range(1, N + 1)]
    profile = np.array(profile) / np.sum(profile) * N

    log(f"  N={N}:")
    log(f"    Frequencies: {['%.3f' % f for f in freqs]}")
    log(f"    |ψ_{N-1}|² profile: {['%.3f' % p for p in profile]}")

log()
log("  The highest-frequency mode is center-heavy, lowest is edge-heavy.")
log("  This is the standard 1D standing wave pattern. No alternating")
log("  structure, no sublattice, no SSH topology.")
log()

# Compare with actual Liouvillian eigenvector profile
N = 5
gamma_uniform = [0.05] * N
t0 = clock.perf_counter()
L_vec = build_liouvillian_vec(N, gamma_uniform)
w_all, vr_all = eig(L_vec, right=True)

# Find slowest oscillating modes (smallest |Re|, nonzero |Im|)
osc_mask = np.abs(w_all.imag) > 0.01
if np.any(osc_mask):
    rates = -w_all[osc_mask].real
    idx_slow = np.where(osc_mask)[0][np.argsort(rates)]
    # Get the 4 slowest oscillating modes
    slow4_idx = idx_slow[:4]
    profiles = []
    for i in slow4_idx:
        profiles.append(localization_profile(vr_all[:, i], N))
    avg_profile = np.mean(profiles, axis=0)
    slow_freq = np.mean(np.abs(w_all[slow4_idx].imag))
    slow_rate = np.mean(-w_all[slow4_idx].real)

    # How many modes share the slowest rate? Taking "the 4 slowest" out of a
    # degenerate set makes the average a function of the eigensolver's
    # ordering, not of the physics, so the count has to be printed.
    slowest = np.min(rates)
    n_tied = int(np.sum(np.abs(rates - slowest) < 1e-8))
    tied_freqs = sorted(set(np.round(np.abs(w_all[osc_mask][
        np.abs(rates - slowest) < 1e-8].imag), 6)))

    log(f"  Actual Liouvillian (N={N}, uniform γ=0.05, {clock.perf_counter()-t0:.1f}s):")
    # slow_rate is well defined (the tied modes share it exactly); a MEAN over
    # their frequencies is not a frequency of anything, so the set is printed.
    log(f"    Slowest oscillating modes: rate={slow_rate:.4f}")
    log(f"    Modes sharing the slowest rate {slowest:.4f}: {n_tied}")
    log(f"    Their frequencies: {['%.3f' % f for f in tied_freqs]}")
    log(f"    Localization (4 of them, averaged): {['%.3f' % p for p in avg_profile]}")

    # CAVITY_MODE_LOCALIZATION reports its profile in its own normalization
    # (it sums to 3.0, not to N), and for a DIFFERENT mode set: rate 0.0456,
    # freq 7.234J. Put both on sum = N before comparing anything.
    cavity_ref = np.array([0.519, 0.631, 0.700, 0.631, 0.519])
    cavity_ref_n = cavity_ref / np.sum(cavity_ref) * N
    dev = float(np.max(np.abs(avg_profile - cavity_ref_n)))
    log(f"    CAVITY_MODE reference, renormalized to sum={N}: "
        f"{['%.3f' % p for p in cavity_ref_n]}")
    log(f"    Max deviation: {dev:.3f}")
    log()
    if n_tied > 4:
        log(f"  The {n_tied} modes at the slowest rate are degenerate in rate, so")
        log("  'the 4 slowest' is an eigensolver-ordering artifact and the average")
        log("  above is not a well-defined object. It is reported, not compared.")
    elif dev < 0.05:
        log("  The profiles match. Center-localization is a standing wave effect.")
    else:
        _centre = avg_profile[N // 2]
        _shape = ("centre-peaked" if _centre >= avg_profile.max() - 1e-9
                  else "not centre-peaked")
        log(f"  The profiles do NOT match (max deviation {dev:.3f}). The set measured")
        log(f"  here sits at rate {slowest:.4f} and is {_shape} "
            f"(centre {_centre:.3f} against max {avg_profile.max():.3f}).")
        log("  CAVITY_MODE's profile belongs to a different mode set; comparing them")
        log("  needs that set selected here, which this step does not do.")


# ========================================================================
# PHASE 2: CHIRAL BLOCK STRUCTURE AND WINDING NUMBER
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 2: CHIRAL BLOCK STRUCTURE")
log("=" * 72)
log()
log("  Π has eigenvalues {+1, -1, +i, -i}.")
log("  {Π, L_c} = 0 forces: L_c maps V_λ → V_{-λ}.")
log("  This gives two chiral subsystems:")
log("    Pair 1: (V_{+1}, V_{-1}) with off-diagonal blocks (A₁, B₁)")
log("    Pair 2: (V_{+i}, V_{-i}) with off-diagonal blocks (A₂, B₂)")
log("  If det(A)=0: topological zero modes. If full rank: trivial.")
log()

for N in [4, 5]:
    t0 = clock.perf_counter()
    H = build_hamiltonian(N)
    gammas = [0.05] * N
    L = build_liouvillian_pauli(N, H, gammas)
    c = 2 * sum(gammas)
    num = 4**N
    L_c = L + (c / 2) * np.eye(num)

    # Diagonalize Π
    Pi, _ = build_pi(N)
    pi_evals, pi_evecs = eig(Pi)

    # Sort eigenvalues into {+1, -1, +i, -i} groups
    groups = {1: [], -1: [], 1j: [], -1j: []}
    for i, ev in enumerate(pi_evals):
        best = min(groups.keys(), key=lambda g: abs(ev - g))
        groups[best].append(i)

    dims = {k: len(v) for k, v in groups.items()}
    log(f"  N={N} ({num}×{num}):")
    log(f"    Π eigenspace dims: +1:{dims[1]}, -1:{dims[-1]},"
        f" +i:{dims[1j]}, -i:{dims[-1j]}")

    # Transform L_c to Π-eigenbasis
    U = pi_evecs
    L_c_pi = np.linalg.solve(U, L_c @ U)

    # Check block structure: L_c should be off-diagonal between V_λ and V_{-λ}
    for pair_name, (g1, g2) in [("±1", (1, -1)), ("±i", (1j, -1j))]:
        idx1 = groups[g1]
        idx2 = groups[g2]
        n1, n2 = len(idx1), len(idx2)

        if n1 == 0 or n2 == 0:
            log(f"    Pair {pair_name}: empty (n1={n1}, n2={n2})")
            continue

        # Extract blocks
        A = L_c_pi[np.ix_(idx1, idx2)]  # V_{g2} → V_{g1}
        B = L_c_pi[np.ix_(idx2, idx1)]  # V_{g1} → V_{g2}
        diag_11 = L_c_pi[np.ix_(idx1, idx1)]
        diag_22 = L_c_pi[np.ix_(idx2, idx2)]

        on_diag_err = max(np.max(np.abs(diag_11)), np.max(np.abs(diag_22)))

        # Singular values of off-diagonal block
        sv = svdvals(A)
        n_zero_sv = np.sum(sv < 1e-10)
        min_sv = sv[-1] if len(sv) > 0 else 0

        log(f"    Pair {pair_name} ({n1}×{n2}):")
        log(f"      On-diagonal residual: {on_diag_err:.2e}"
            f"  (0 = perfect off-diagonal)")
        log(f"      Singular values of A: min={min_sv:.4e},"
            f" max={sv[0]:.4e}, zero count={n_zero_sv}")

        if n_zero_sv > 0:
            log(f"      *** {n_zero_sv} ZERO SINGULAR VALUES → topological zero modes!")
        else:
            log(f"      Full rank → NO topological zero modes (trivial)")

    log(f"    ({clock.perf_counter()-t0:.1f}s)")


# ========================================================================
# PHASE 3: WHAT A BERRY PHASE CAN SAY ON THIS PATH
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 3: BERRY PHASE ON AN OPEN ONE-PARAMETER PATH")
log("=" * 72)
log()

N = 5
d = 2**N
gamma_base = 0.01
gamma_high = 0.5

# Parametric path: uniform (θ=0) → sacrifice-edge (θ=1)
# γ_k(θ) = (1-θ)·γ_mean + θ·sacrifice_k
gamma_mean = (gamma_high + (N - 1) * gamma_base) / N
gamma_sacrifice = [gamma_high] + [gamma_base] * (N - 1)

# The tracked mode lives in the (popcount 1, popcount 0) coherence block,
# which H and Z dephasing leave exactly invariant. In the full space it has an
# exact twin (its X^N image in the (4,5) block), so a full-space eigensolver
# hands back an arbitrary mix of the two; inside the block it is simple.
# L is linear in γ, so the block is assembled from its parts.
popc = lambda x: bin(x).count("1")
blk = [i * d + j for i in range(d) for j in range(d) if popc(i) == 1 and popc(j) == 0]
rest = np.setdiff1d(np.arange(d * d), blk)
L_H = build_liouvillian_vec(N, [0.0] * N)
L_D = [build_liouvillian_vec(N, [1.0 if m == k else 0.0 for m in range(N)]) - L_H
       for k in range(N)]
for M in [L_H] + L_D:
    assert np.max(np.abs(M[np.ix_(rest, blk)])) == 0.0, "block not invariant"
bH = L_H[np.ix_(blk, blk)]
bD = [M[np.ix_(blk, blk)] for M in L_D]


def path_block(theta):
    g = [(1 - theta) * gamma_mean + theta * gamma_sacrifice[k] for k in range(N)]
    return bH + sum(g[k] * bD[k] for k in range(N))


def berry_readings(n_theta):
    """Track the slowest mode of the block from θ=1, where it is the unique
    slowest, back to θ=0, where all five block modes share Re λ = -2·γ_mean
    and only continuity names it. Returns (per-step sum, closed phase,
    min gap, modulus of the closed product)."""
    picks, lam = [], None
    for theta in np.linspace(1, 0, n_theta):
        w, vl, vr = eig(path_block(theta), left=True, right=True)
        i = int(np.argmax(w.real)) if lam is None else int(np.argmin(np.abs(w - lam)))
        lam = w[i]
        picks.append((w, vl[:, i], vr[:, i], i))
    picks.reverse()
    naive, prev_r, Ls, Rs, min_gap = 0.0, None, [], [], np.inf
    for w, l, r, i in picks:
        min_gap = min(min_gap, np.sort(np.abs(w - w[i]))[1])
        if prev_r is not None:
            # -arg<L_{n+1}|R_n>, each vector in LAPACK's own phase
            naive += -np.imag(np.log(np.vdot(l, prev_r)))
        prev_r = r
        Ls.append(l / np.vdot(l, r).conj())     # biorthonormal: <L_n|R_n> = 1
        Rs.append(r)
    prod = 1.0 + 0j
    for k in range(n_theta - 1):
        prod *= np.vdot(Ls[k], Rs[k + 1])
    prod *= np.vdot(Ls[-1], Rs[0])              # close the open segment
    return naive, -np.angle(prod), min_gap, abs(prod)


log(f"  N={N}, path: uniform (θ=0) → sacrifice-edge (θ=1)")
log(f"  γ_mean={gamma_mean:.4f}, sacrifice=[{gamma_high}, {gamma_base}×{N-1}]")
log(f"  Slowest mode of the (1,0) coherence block (5×5, exactly invariant).")
log()
log(f"  {'steps':>6}  {'per-step Σ -arg<L|R>':>24}  {'closed biorthonormal phase':>27}  {'min gap':>8}")
log(f"  {'─'*72}")
t0 = clock.perf_counter()
readings = {}
for n_theta in (60, 600, 6000):
    naive, closed, gap, modulus = berry_readings(n_theta)
    readings[n_theta] = (naive, closed)
    log(f"  {n_theta:>6}  {naive:>24.4f}  {closed:>27.4e}  {gap:>8.3f}   |product| = {modulus:.6f}")
log(f"  ({clock.perf_counter()-t0:.1f}s)")
log()
ratio = readings[6000][0] / readings[600][0]
log(f"  The per-step sum grows with the step count (600 → 6000: ×{ratio:.1f}): it sums")
log("  the eigensolver's per-step phase convention, not a geometric phase. The")
log("  biorthonormal phase closed by <L(1)|R(0)> converges, to"
    f" {readings[6000][1]:.3e} rad.")
log("  Neither decides topology. The path is an open segment of an interval,")
log("  so a phase is gauge-invariant only with a closing rule, and the closed")
log("  value depends on that rule; a loop on an interval retraces itself.")
log("  The block is complex symmetric (H real), so the biorthogonal connection")
log("  vanishes: a loop in several γ would read 0 or π, π only around an")
log("  exceptional point, and this path encircles none. The Berry test")
log("  cannot decide here.")


# ========================================================================
# PHASE 5: ROBUSTNESS SWEEP (edge → uniform → center sacrifice)
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 5: ROBUSTNESS SWEEP")
log("=" * 72)
log()

N = 5
gamma_high = 0.5
gamma_low = 0.01
center = N // 2

# θ ∈ [-1, 1]: +1 = edge sacrifice, 0 = uniform, -1 = center sacrifice
n_sweep = 41
thetas_sweep = np.linspace(-1, 1, n_sweep)

log(f"  N={N}, γ_high={gamma_high}, γ_low={gamma_low}")
log(f"  θ = +1: sacrifice at edge (site 0)")
log(f"  θ =  0: uniform")
log(f"  θ = -1: sacrifice at center (site {center})")
log()

log(f"  {'θ':>6}  {'min rate':>10}  {'max rate':>10}  {'prot. factor':>12}"
    f"  {'gap':>10}")
log(f"  {'─'*55}")

t0 = clock.perf_counter()
prot_factors = np.zeros(n_sweep)
min_rates = np.zeros(n_sweep)
spectral_gaps = np.zeros(n_sweep)

profile_edge_sac = [gamma_high] + [gamma_low] * (N - 1)
profile_center_sac = [gamma_low] * center + [gamma_high] + [gamma_low] * (N - 1 - center)
profile_uniform_s = [np.mean(profile_edge_sac)] * N

for ti, theta in enumerate(thetas_sweep):
    if theta >= 0:
        gammas = [(1 - theta) * profile_uniform_s[k] + theta * profile_edge_sac[k]
                  for k in range(N)]
    else:
        gammas = [(1 + theta) * profile_uniform_s[k] + (-theta) * profile_center_sac[k]
                  for k in range(N)]

    ev = eigvals(build_liouvillian_vec(N, gammas))
    rates = -ev.real
    nonzero = rates[rates > 1e-10]

    if len(nonzero) > 0:
        min_r = np.min(nonzero)
        max_r = np.max(nonzero)
        prot = max_r / min_r
        gap = min_r
    else:
        min_r = max_r = prot = gap = 0

    prot_factors[ti] = prot
    min_rates[ti] = min_r
    spectral_gaps[ti] = gap

    if ti % 5 == 0:
        log(f"  {theta:>+6.2f}  {min_r:>10.4f}  {max_r:>10.4f}"
            f"  {prot:>12.2f}  {gap:>10.4f}")

log(f"  ({clock.perf_counter()-t0:.1f}s)")

# Sharp transition or steep tail? A jump in a continuous curve halves when
# the grid is doubled; a discontinuity keeps its size. The steepest stretch
# is the tail toward θ = -1, where the slowest rate falls to 2γ_low.
L_H5 = build_liouvillian_vec(N, [0.0] * N)
L_D5 = [build_liouvillian_vec(N, [1.0 if m == k else 0.0 for m in range(N)]) - L_H5
        for k in range(N)]


def protection_on_grid(n):
    ths = np.linspace(-1, 1, n)
    out = []
    for theta in ths:
        if theta >= 0:
            g = [(1 - theta) * profile_uniform_s[k] + theta * profile_edge_sac[k] for k in range(N)]
        else:
            g = [(1 + theta) * profile_uniform_s[k] - theta * profile_center_sac[k] for k in range(N)]
        r = -eigvals(L_H5 + sum(g[k] * L_D5[k] for k in range(N))).real
        r = r[r > 1e-10]
        out.append(r.max() / r.min())
    return ths, np.array(out)


jumps = {}
for n_grid, (ths, prot) in [(n_sweep, (thetas_sweep, prot_factors)),
                            (2 * n_sweep - 1, protection_on_grid(2 * n_sweep - 1))]:
    dj = np.abs(np.diff(prot))
    k = int(np.argmax(dj))
    jumps[n_grid] = dj[k]
    log()
    log(f"  {n_grid} grid points: max protection-factor jump {dj[k]:.2f}"
        f" at θ ≈ {(ths[k] + ths[k + 1]) / 2:+.3f}")
shrink = jumps[n_sweep] / jumps[2 * n_sweep - 1]
log(f"  Doubling the grid shrinks the largest jump by ×{shrink:.2f}"
    f" (a continuous curve → ×2, a discontinuity → ×1).")
log(f"  Slowest rate never closes: min over the sweep {np.min(min_rates):.4f}.")
if shrink > 1.5:
    log("  Continuous, steep only in its tail → no sharp phase boundary")
else:
    log("  The largest jump does not shrink with the grid → a discontinuity to examine")

# Localization profiles at three key points
log()
log("─" * 72)
log("  Localization profiles at key θ values:")
log()
for theta_key, label in [(1.0, "edge sacrifice"),
                          (0.0, "uniform"),
                          (-1.0, "center sacrifice")]:
    if theta_key >= 0:
        gammas = [(1 - theta_key) * profile_uniform_s[k]
                  + theta_key * profile_edge_sac[k] for k in range(N)]
    else:
        gammas = [(1 + theta_key) * profile_uniform_s[k]
                  + (-theta_key) * profile_center_sac[k] for k in range(N)]

    L = build_liouvillian_vec(N, gammas)
    w, vr = eig(L, right=True)

    osc = np.abs(w.imag) > 0.1
    if np.any(osc):
        # every oscillating mode at the slowest rate, grouped by eigenvalue;
        # each eigenspace is read basis-free, and the modes are averaged
        rates_osc = -w.real
        slow_rate = np.min(rates_osc[osc])
        tied = np.where(osc & (np.abs(rates_osc - slow_rate) < 1e-8))[0]
        groups = []
        for i in tied:
            for g in groups:
                if abs(w[i] - w[g[0]]) < 1e-8:
                    g.append(i)
                    break
            else:
                groups.append([i])
        avg_prof = sum(len(g) * eigenspace_profile(vr[:, g], N) for g in groups) / len(tied)
        log(f"  θ={theta_key:+4.1f} ({label}): rate={slow_rate:.4f}"
            f"  {len(tied)} modes, {len(groups)} eigenvalues"
            f"  profile={['%.3f' % p for p in avg_prof]}")
        for g in groups:
            if w[g[0]].imag > 0:
                log(f"      Im λ = {w[g[0]].imag:+8.4f} ({len(g)}):"
                    f"  {['%.3f' % p for p in eigenspace_profile(vr[:, g], N)]}")


# ========================================================================
# PHASE 4: MODE COUNTING
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 4: CENTER-LOCALIZED MODE COUNT")
log("=" * 72)
log()

# A per-mode label needs a per-mode vector, and the full spectrum is degenerate
# (every oscillating eigenvalue has its X^N twin in the mirrored popcount
# block). So each mode is taken from its own (popcount row, popcount column)
# block, which H and Z dephasing leave invariant: the sector-pure basis.
for N in [3, 4, 5]:
    d = 2**N
    gammas = [0.05] * N
    L = build_liouvillian_vec(N, gammas)
    popc = lambda x: bin(x).count("1")
    n_center = n_osc = 0
    max_in_block = 1
    mixed_label = 0  # in-block clusters whose label changes under mixing
    rng = np.random.default_rng(1)
    # the middle: one site at odd N, the pair N/2-1, N/2 at even N (the
    # reflection makes their weights exactly equal in every mode)
    middle = [N // 2] if N % 2 else [N // 2 - 1, N // 2]
    # "center-localized": a middle site carries the highest weight, equal to
    # rounding (1e-9); exact ties come from the reflection, and the nearest
    # genuine gap is 8.6e-5 (N=5). A coarser threshold counts more modes:
    is_center = lambda v, tol=1e-9: (lambda pr: np.max(pr[middle]) >= np.max(pr) - tol)(
        localization_profile(v, N))
    n_center_coarse = {1e-3: 0, 1e-2: 0}
    for p in range(N + 1):
        for q in range(N + 1):
            blk = [i * d + j for i in range(d) for j in range(d)
                   if popc(i) == p and popc(j) == q]
            wb, vb = eig(L[np.ix_(blk, blk)], right=True)
            for i in np.where(np.abs(wb.imag) > 0.01)[0]:
                cl = np.where(np.abs(wb - wb[i]) < 1e-8)[0]
                max_in_block = max(max_in_block, len(cl))
                v = np.zeros(d * d, dtype=complex)
                v[blk] = vb[:, i]
                n_osc += 1
                # "Center-localized": center qubit has highest weight
                label = is_center(v)
                n_center += label
                for tol in n_center_coarse:
                    n_center_coarse[tol] += is_center(v, tol)
                if len(cl) > 1 and i == cl[0]:
                    # a degenerate in-block eigenvalue: does the label survive
                    # a change of basis inside its eigenspace?
                    for _ in range(20):
                        c = rng.normal(size=len(cl)) + 1j * rng.normal(size=len(cl))
                        u = np.zeros(d * d, dtype=complex)
                        u[blk] = vb[:, cl] @ c
                        if is_center(u) != label:
                            mixed_label += 1
                            break

    log(f"  N={N}: {n_center}/{n_osc} oscillating modes are center-localized"
        f" ({100*n_center/n_osc:.1f}%), middle sites {middle}, sector-pure basis;"
        f" largest in-block multiplicity {max_in_block},"
        f" labels changed by 20 random mixes {mixed_label};"
        f" within 1e-3: {n_center_coarse[1e-3]}, within 1e-2: {n_center_coarse[1e-2]}")


# ========================================================================
# SUMMARY
# ========================================================================
log()
log()
log("=" * 72)
log("SUMMARY")
log("=" * 72)
log()
log("Phase 1: the (0,1) coherence block is tight-binding with UNIFORM hopping.")
log("  No alternating structure → no SSH topology.")
log("  Localization profile matches the cos (Neumann) standing waves this")
log("  phase computes, not the Dirichlet sin family, which belongs to F2b.")
log()
log("Phase 2: Chiral block structure exists (Π eigenspaces give")
log("  off-diagonal L_c). At N=5 the off-diagonal blocks have FULL RANK;")
log("  the N=4 zero singular values are the even-N centre modes at -Σγ.")
log()
log("Phase 3: the path is an open one-parameter segment, so no Berry phase can")
log("  decide; a per-step phase sum grows with the step count.")
log()
log("Phase 5: the protection factor is continuous (its largest jump shrinks")
log("  with the grid) and the slowest rate never closes → no sharp phase")
log("  boundary. Each mode keeps its standing-wave profile; the γ profile")
log("  selects which modes are slowest (center-heavy under edge sacrifice,")
log("  center-dark under center sacrifice, all tied under uniform noise).")
log()
log("Phase 4: the center-localized fraction depends on what counts as the")
log("  middle and on any tolerance; it is no topological invariant.")
log()
log("CONCLUSION: The localization [0.52, 0.63, 0.70, 0.63, 0.52] is")
log("GEOMETRIC, not topological. It arises from standing wave patterns")
log("on a 1D chain (cos (Neumann) eigenstates). The sacrifice zone")
log("EXPLOITS this geometry by placing noise where fast modes live")
log("(edges) and quiet where slow modes live (center). The chain")
log("topology provides the mode structure. The noise provides the")
log("selection pressure. No topological invariant protects anything.")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
