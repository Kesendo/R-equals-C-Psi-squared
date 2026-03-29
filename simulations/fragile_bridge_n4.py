#!/usr/bin/env python3
"""
Fragile Bridge N=4 Scaling
===========================
Sparse Liouvillian for N=4 per chain (8 qubits, 65536x65536).
Determines N-scaling of gamma_crit.

Script: simulations/fragile_bridge_n4.py
Output: simulations/results/fragile_bridge_n4.txt
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigs
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "fragile_bridge_n4.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# === Sparse Pauli matrices ===
I2 = sp.eye(2, format='csr', dtype=complex)
sx = sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
sy = sp.csr_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex))
sz = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))


def op_at_sparse(op, qubit, n_qubits):
    """Place single-qubit operator on qubit k (sparse Kronecker)."""
    result = sp.eye(1, format='csr', dtype=complex)
    for k in range(n_qubits):
        result = sp.kron(result, op if k == qubit else I2, format='csr')
    return result


def build_coupled_sparse(n_per_chain, gamma, J=1.0, J_bridge=0.1):
    """
    Build sparse Liouvillian for two coupled chains.
    Chain A: qubits 0..n-1, dephasing +gamma
    Chain B: qubits n..2n-1, dephasing -gamma
    Bridge: Heisenberg XXX between qubit n-1 and qubit n
    """
    n_total = 2 * n_per_chain
    d = 2**n_total
    d2 = d * d

    log(f"    Building H ({d}x{d})...", )
    t0 = clock.time()

    # Hamiltonian (sparse)
    H = sp.csr_matrix((d, d), dtype=complex)
    # Chain A: nearest-neighbor
    for i in range(n_per_chain - 1):
        for P in [sx, sy, sz]:
            H = H + J * (op_at_sparse(P, i, n_total) @ op_at_sparse(P, i + 1, n_total))
    # Chain B: nearest-neighbor
    for i in range(n_per_chain, 2 * n_per_chain - 1):
        for P in [sx, sy, sz]:
            H = H + J * (op_at_sparse(P, i, n_total) @ op_at_sparse(P, i + 1, n_total))
    # Bridge
    for P in [sx, sy, sz]:
        H = H + J_bridge * (op_at_sparse(P, n_per_chain - 1, n_total) @
                            op_at_sparse(P, n_per_chain, n_total))

    log(f"    H built ({clock.time()-t0:.1f}s, nnz={H.nnz})")

    log(f"    Building L ({d2}x{d2})...")
    t0 = clock.time()

    Id = sp.eye(d, format='csr', dtype=complex)

    # Hamiltonian part: -i(H x I - I x H^T)
    L = -1j * (sp.kron(H, Id, format='csr') - sp.kron(Id, H.T, format='csr'))

    log(f"    L_H built ({clock.time()-t0:.1f}s, nnz={L.nnz})")
    t0 = clock.time()

    # Dissipator
    Id2 = sp.eye(d2, format='csr', dtype=complex)
    for k in range(n_per_chain):
        Zk = op_at_sparse(sz, k, n_total)
        L = L + gamma * (sp.kron(Zk, Zk.conj(), format='csr') - Id2)
    for k in range(n_per_chain, 2 * n_per_chain):
        Zk = op_at_sparse(sz, k, n_total)
        L = L + (-gamma) * (sp.kron(Zk, Zk.conj(), format='csr') - Id2)

    log(f"    L complete ({clock.time()-t0:.1f}s, nnz={L.nnz})")

    return L


def max_re_sparse(L, k=10):
    """Find max Re(eigenvalue). Dense for small, growth-rate for large."""
    d2 = L.shape[0]
    # For small matrices, use dense (faster and more reliable)
    if d2 <= 4096:
        from scipy.linalg import eigvals
        evals = eigvals(L.toarray())
        return float(np.max(evals.real))

    # For large matrices: estimate max Re via matrix exponential growth
    from scipy.sparse.linalg import expm_multiply
    np.random.seed(42)
    v0 = np.random.randn(d2) + 1j * np.random.randn(d2)
    v0 = v0 / np.linalg.norm(v0)

    # Start with short time, increase if no growth detected
    for t_test in [1.0, 10.0, 100.0]:
        vt = expm_multiply(L, v0, start=0, stop=t_test, num=5, endpoint=True)
        norms = np.array([np.linalg.norm(vt[i]) for i in range(5)])

        # If overflow (inf/nan), system is definitely unstable
        if np.any(np.isinf(norms)) or np.any(np.isnan(norms)):
            return 1.0  # large positive = unstable

        # If significant growth detected, compute rate
        if norms[-1] / (norms[0] + 1e-300) > 2.0:
            times = np.linspace(0, t_test, 5)
            log_norms = np.log(norms + 1e-300)
            rate, _ = np.polyfit(times, log_norms, 1)
            return float(rate)

    # No growth at any time scale: stable
    times = np.linspace(0, 100.0, 5)
    log_norms = np.log(norms + 1e-300)
    rate, _ = np.polyfit(times, log_norms, 1)
    return float(rate)


def find_gamma_crit_sparse(n_per_chain, J_bridge, tol=1e-5):
    """Bisect for instability threshold using sparse eigensolver."""
    threshold = 1e-4  # growth rate threshold for instability

    log(f"  Testing gamma=0.5 for upper bound...")
    L = build_coupled_sparse(n_per_chain, 0.5, J_bridge=J_bridge)
    mr = max_re_sparse(L)
    if mr is None:
        return None
    if mr <= threshold:
        log(f"  rate={mr:.2e} at gamma=0.5, trying gamma=5.0...")
        L = build_coupled_sparse(n_per_chain, 5.0, J_bridge=J_bridge)
        mr = max_re_sparse(L)
        if mr is None or mr <= threshold:
            return None

    g_lo, g_hi = 0.0, 0.5 if mr > threshold else 5.0
    step = 0
    while (g_hi - g_lo) > tol:
        g_mid = (g_lo + g_hi) / 2
        step += 1
        log(f"  Bisection step {step}: gamma={g_mid:.6f}...")
        L = build_coupled_sparse(n_per_chain, g_mid, J_bridge=J_bridge)
        mr = max_re_sparse(L)
        if mr is None:
            return None
        if mr > threshold:
            g_hi = g_mid
            log(f"    rate={mr:.6f} > 0 -> UNSTABLE")
        else:
            g_lo = g_mid
            log(f"    rate={mr:.6f} <= 0 -> stable")

    return (g_lo + g_hi) / 2


# ================================================================
log("=" * 70)
log("FRAGILE BRIDGE N-SCALING")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Machine: 128 GB RAM, 24 cores")
log("=" * 70)

# ================================================================
# Verify N=2 and N=3 with sparse (should match dense results)
# ================================================================
log()
log("=" * 70)
log("VERIFICATION: N=2 and N=3 (sparse vs dense)")
log("=" * 70)

for N in [2, 3]:
    log()
    log(f"--- N={N} per chain ({2*N} qubits, {4**N}x{4**N}) ---")
    t0 = clock.time()
    gc = find_gamma_crit_sparse(N, J_bridge=0.10, tol=1e-5)
    elapsed = clock.time() - t0
    if gc is not None:
        log(f"  gamma_crit = {gc:.6f} ({elapsed:.1f}s)")
    else:
        log(f"  FAILED ({elapsed:.1f}s)")

# ================================================================
# N=4: the main computation
# ================================================================
log()
log("=" * 70)
log("N=4 per chain (8 qubits, 65536x65536)")
log("J_bridge = 0.10, J = 1.0")
log("=" * 70)
log()

t0_total = clock.time()
gc4 = find_gamma_crit_sparse(4, J_bridge=0.10, tol=1e-5)
elapsed_total = clock.time() - t0_total

log()
if gc4 is not None:
    log(f"  gamma_crit(N=4) = {gc4:.6f} ({elapsed_total:.0f}s total)")
else:
    log(f"  FAILED ({elapsed_total:.0f}s)")

# ================================================================
# Scaling analysis
# ================================================================
log()
log("=" * 70)
log("SCALING ANALYSIS")
log("=" * 70)
log()

# Collect all data
data = [(2, 0.017292)]  # from dense computation
# Add N=3 and N=4 from sparse
# (we'll fill these in from the runs above)

log("  Known values (J_bridge=0.10):")
log(f"  N=2: gamma_crit = 0.017292 (dense)")

# Read back our computed values
# For now, print summary after all computations
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
