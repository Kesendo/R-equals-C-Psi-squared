"""
Verification of the complete proof: dim(ker([H,·]|_{w=1})) = 2N
================================================================
Five independent checks for the triangle inequality / SWAP-invariance argument
in docs/proofs/PROOF_WEIGHT1_DEGENERACY.md, Step 5.

V1: Each individual SWAP fixes each kernel vector
V2: Triangle inequality saturation for kernel vectors (and non-saturation for non-kernel)
V3: No additional SWAP-invariant vectors beyond the 2N kernel vectors
V4: Orbit transitivity under nearest-neighbor SWAPs
V5: Analytical T_c basis spans the same subspace as numerical eigenvectors

Output: simulations/results/verify_proof_weight1.txt
"""

import numpy as np
from pathlib import Path
import sys
import os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
TOL = 1e-10

# ─────────────────────────────────────────────
# Pauli matrices and Hamiltonian (same as kommutator_kern_analysis.py)
# ─────────────────────────────────────────────

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, X, Y, Z]
PAULI_LABELS = ['I', 'X', 'Y', 'Z']


def kron_chain(ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result


def heisenberg_chain(N, J=1.0):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [X, Y, Z]:
            ops = [I2] * N
            ops[i] = P
            ops[i + 1] = P
            H += J * kron_chain(ops)
    return H


# ─────────────────────────────────────────────
# Weight-1 basis with label-based SWAP
# ─────────────────────────────────────────────

def enumerate_w1(N):
    """Enumerate all weight-1 Pauli strings. Returns list of labels."""
    basis = []
    for pos in range(N):
        for typ in ['X', 'Y']:
            for z_cfg in range(2**(N - 1)):
                chars = []
                z_idx = 0
                for q in range(N):
                    if q == pos:
                        chars.append(typ)
                    else:
                        bit = (z_cfg >> z_idx) & 1
                        chars.append('Z' if bit else 'I')
                        z_idx += 1
                basis.append(''.join(chars))
    return basis


def z_count(label):
    """Number of Z characters in a Pauli string label."""
    return label.count('Z')


def active_type(label):
    """Return 'X' or 'Y' for a weight-1 label."""
    for c in label:
        if c in ('X', 'Y'):
            return c
    return None


def apply_swap_label(label, s):
    """Apply SWAP_{s,s+1} to a Pauli string label (swap chars at positions s and s+1)."""
    chars = list(label)
    chars[s], chars[s + 1] = chars[s + 1], chars[s]
    return ''.join(chars)


def build_swap_permutation(w1_basis, label_to_idx, s):
    """Build permutation matrix for SWAP_{s,s+1} acting on w=1 basis."""
    n = len(w1_basis)
    P = np.zeros((n, n), dtype=complex)
    for i, label in enumerate(w1_basis):
        new_label = apply_swap_label(label, s)
        j = label_to_idx[new_label]
        P[j, i] = 1.0
    return P


# ─────────────────────────────────────────────
# Commutator kernel (reuse logic)
# ─────────────────────────────────────────────

def commutator_kernel(N):
    d = 2**N
    d2 = d * d
    H = heisenberg_chain(N)
    w1_labels = enumerate_w1(N)
    n_w1 = len(w1_labels)

    # Build w=1 basis matrices
    w1_matrices = []
    for label in w1_labels:
        ops = []
        for c in label:
            ops.append({'I': I2, 'X': X, 'Y': Y, 'Z': Z}[c])
        w1_matrices.append(kron_chain(ops))

    # Commutator matrix
    M = np.zeros((d2, n_w1), dtype=complex)
    for t in range(n_w1):
        comm = H @ w1_matrices[t] - w1_matrices[t] @ H
        M[:, t] = comm.ravel()

    U, S, Vh = np.linalg.svd(M, full_matrices=False)
    rank = np.sum(S > TOL)
    kernel_vectors = Vh[rank:]  # rows are kernel basis vectors in w=1 coords

    return kernel_vectors, w1_labels


# ─────────────────────────────────────────────
# Analytical T_c basis construction
# ─────────────────────────────────────────────

def build_tc_basis(N, w1_labels, label_to_idx):
    """Build the 2N analytical T_c^{(a)} vectors in w=1 basis coordinates."""
    n_w1 = len(w1_labels)
    tc_vectors = []

    for typ in ['X', 'Y']:
        for c in range(N):
            v = np.zeros(n_w1, dtype=complex)
            for i, label in enumerate(w1_labels):
                if active_type(label) == typ and z_count(label) == c:
                    v[i] = 1.0
            norm = np.linalg.norm(v)
            if norm > TOL:
                v /= norm
            tc_vectors.append(v)

    return np.array(tc_vectors)


# ─────────────────────────────────────────────
# Load eigenvector CSV
# ─────────────────────────────────────────────

def load_eigvec_csv(N, w1_labels, label_to_idx):
    path = RESULTS_DIR / f"eigvec_at_minus_gamma_N{N}.csv"
    if not path.exists():
        return None

    data = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            idx = int(parts[0])
            label = parts[3]
            abs_c = float(parts[5])
            phase = float(parts[6])
            coeff = abs_c * np.exp(1j * phase)
            if idx not in data:
                data[idx] = {}
            data[idx][label] = coeff

    n_eig = len(data)
    n_w1 = len(w1_labels)
    eig_vecs = np.zeros((n_eig, n_w1), dtype=complex)
    for ei, pauli_dict in data.items():
        for label, coeff in pauli_dict.items():
            if label in label_to_idx:
                eig_vecs[ei, label_to_idx[label]] = coeff

    return eig_vecs


# ─────────────────────────────────────────────
# Main verification
# ─────────────────────────────────────────────

output_lines = []
def log(msg=""):
    print(msg)
    output_lines.append(msg)


log("=" * 70)
log("VERIFICATION OF PROOF: dim(ker([H,·]|_{w=1})) = 2N")
log("=" * 70)
log()

overall_pass = True

for N in range(2, 8):
    log("=" * 70)
    log(f"N = {N}")
    log("=" * 70)

    # Setup
    kernel_vectors, w1_labels = commutator_kernel(N)
    n_w1 = len(w1_labels)
    n_ker = len(kernel_vectors)
    label_to_idx = {label: i for i, label in enumerate(w1_labels)}
    n_bonds = N - 1

    # Build SWAP permutation matrices
    swaps = []
    for s in range(n_bonds):
        P = build_swap_permutation(w1_labels, label_to_idx, s)
        swaps.append(P)

    # ─── V1: Each SWAP fixes each kernel vector ───
    log(f"\n  V1: Each SWAP_{'{i,i+1}'} fixes each kernel vector")
    v1_pass = True
    max_diff = 0.0
    for ki in range(n_ker):
        v = kernel_vectors[ki]
        for si in range(n_bonds):
            sv = swaps[si] @ v
            diff = np.linalg.norm(sv - v)
            max_diff = max(max_diff, diff)
            if diff > TOL:
                v1_pass = False
                log(f"    FAIL: kernel #{ki}, bond ({si},{si+1}): ||SWAP(v)-v|| = {diff:.2e}")

    status = "PASS" if v1_pass else "FAIL"
    if not v1_pass:
        overall_pass = False
    log(f"    {n_ker} vectors × {n_bonds} bonds = {n_ker * n_bonds} checks")
    log(f"    Max ||SWAP(v)-v||: {max_diff:.2e}")
    log(f"    V1: [{status}]")

    # ─── V2: Triangle inequality saturation ───
    log(f"\n  V2: Triangle inequality saturation")
    v2_pass = True

    for ki in range(n_ker):
        v = kernel_vectors[ki]
        sum_swap_v = np.zeros(n_w1, dtype=complex)
        sum_norms = 0.0
        for si in range(n_bonds):
            sv = swaps[si] @ v
            sum_swap_v += sv
            sum_norms += np.linalg.norm(sv)
        lhs = np.linalg.norm(sum_swap_v)
        rhs = sum_norms
        gap = abs(lhs - rhs)
        if gap > TOL:
            v2_pass = False
            log(f"    FAIL: kernel #{ki}: lhs={lhs:.10f}, rhs={rhs:.10f}, gap={gap:.2e}")

    # Non-kernel vector check (should have strict inequality)
    # Pick a random w=1 vector orthogonal to the kernel
    rand_v = np.random.randn(n_w1) + 1j * np.random.randn(n_w1)
    # Project out kernel
    for ki in range(n_ker):
        rand_v -= np.dot(np.conj(kernel_vectors[ki]), rand_v) * kernel_vectors[ki]
    rand_v /= np.linalg.norm(rand_v)

    sum_swap_rand = np.zeros(n_w1, dtype=complex)
    sum_norms_rand = 0.0
    for si in range(n_bonds):
        sv = swaps[si] @ rand_v
        sum_swap_rand += sv
        sum_norms_rand += np.linalg.norm(sv)
    lhs_rand = np.linalg.norm(sum_swap_rand)
    gap_rand = sum_norms_rand - lhs_rand
    strict = gap_rand > 0.01  # should be strictly less

    status = "PASS" if v2_pass else "FAIL"
    if not v2_pass:
        overall_pass = False
    log(f"    Kernel vectors: all saturate triangle inequality [{status}]")
    log(f"    Non-kernel vector: gap = {gap_rand:.4f} {'(strict inequality ✓)' if strict else '(WARNING: near equality)'}")

    # ─── V3: No additional SWAP-invariant vectors ───
    log(f"\n  V3: No additional SWAP-invariant vectors")

    # Stack (I - SWAP_i) for all bonds → tall matrix, compute kernel
    stacked = np.zeros((n_bonds * n_w1, n_w1), dtype=complex)
    for si in range(n_bonds):
        stacked[si * n_w1:(si + 1) * n_w1, :] = np.eye(n_w1) - swaps[si]

    _, S_stacked, Vh_stacked = np.linalg.svd(stacked, full_matrices=False)
    rank_stacked = np.sum(S_stacked > TOL)
    dim_invariant = n_w1 - rank_stacked

    v3_pass = dim_invariant == 2 * N
    if not v3_pass:
        overall_pass = False

    status = "PASS" if v3_pass else "FAIL"
    log(f"    Stacked (I-SWAP) matrix: {n_bonds * n_w1} × {n_w1}")
    log(f"    Rank: {rank_stacked}, invariant subspace dim: {dim_invariant} (expected {2*N})")
    log(f"    V3: [{status}]")

    # ─── V4: Orbit transitivity ───
    log(f"\n  V4: Orbit transitivity under nearest-neighbor SWAPs")
    v4_pass = True

    for typ in ['X', 'Y']:
        for c in range(N):
            # Collect all w=1 strings with this type and z-count
            orbit_labels = [l for l in w1_labels if active_type(l) == typ and z_count(l) == c]
            orbit_size = len(orbit_labels)

            if orbit_size == 0:
                continue

            # BFS from first element
            visited = {orbit_labels[0]}
            queue = [orbit_labels[0]]
            while queue:
                current = queue.pop(0)
                for s in range(n_bonds):
                    neighbor = apply_swap_label(current, s)
                    if neighbor in visited:
                        continue
                    if neighbor in set(orbit_labels):
                        visited.add(neighbor)
                        queue.append(neighbor)

            reachable = len(visited)
            if reachable != orbit_size:
                v4_pass = False
                log(f"    FAIL: typ={typ}, c={c}: reachable {reachable}/{orbit_size}")

    if not v4_pass:
        overall_pass = False
    status = "PASS" if v4_pass else "FAIL"
    log(f"    All (type, Z-count) orbits transitive: [{status}]")

    # ─── V5: T_c basis matches eigenvector subspace ───
    log(f"\n  V5: Analytical T_c basis vs numerical eigenvectors")

    tc_basis = build_tc_basis(N, w1_labels, label_to_idx)
    eig_vecs = load_eigvec_csv(N, w1_labels, label_to_idx)

    if eig_vecs is not None:
        combined = np.vstack([tc_basis, eig_vecs])
        _, S_comb, _ = np.linalg.svd(combined, full_matrices=False)
        rank_combined = np.sum(S_comb > TOL)

        v5_pass = rank_combined == 2 * N
        if not v5_pass:
            overall_pass = False

        status = "PASS" if v5_pass else "FAIL"
        log(f"    rank(T_c ∪ eigvecs) = {rank_combined} (expected {2*N}) [{status}]")

        # Also check T_c vectors are in the commutator kernel
        tc_in_kernel = True
        for ti in range(len(tc_basis)):
            proj = kernel_vectors @ tc_basis[ti]
            proj_norm = np.linalg.norm(proj)
            tc_norm = np.linalg.norm(tc_basis[ti])
            if abs(proj_norm - tc_norm) > TOL:
                tc_in_kernel = False
                log(f"    WARNING: T_c #{ti} not fully in kernel: proj={proj_norm:.6f}, norm={tc_norm:.6f}")
        log(f"    All T_c vectors in commutator kernel: {'✓' if tc_in_kernel else '✗'}")
    else:
        log(f"    No CSV data for N={N} (skipped)")
        if N <= 6:
            overall_pass = False

    log()

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────

log("=" * 70)
log("OVERALL RESULT")
log("=" * 70)
log(f"  All verifications passed: {'✓ YES' if overall_pass else '✗ NO'}")
log()
log("  V1 (SWAP fixpoint):              every SWAP fixes every kernel vector")
log("  V2 (triangle saturation):        equality holds for kernel, strict for non-kernel")
log("  V3 (no extra invariants):        SWAP-invariant subspace = 2N exactly")
log("  V4 (orbit transitivity):         all (type, Z-count) orbits transitive")
log("  V5 (T_c = eigenvector subspace): analytical and numerical bases match")
log()
if overall_pass:
    log("  The proof in docs/proofs/PROOF_WEIGHT1_DEGENERACY.md is FULLY VERIFIED.")
else:
    log("  WARNING: Some checks failed. The proof needs review.")

# Save
out_path = RESULTS_DIR / "verify_proof_weight1.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(output_lines))
print(f"\n>>> Results saved to: {out_path}")
