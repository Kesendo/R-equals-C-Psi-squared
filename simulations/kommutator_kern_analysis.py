"""
Commutator kernel analysis: dim(ker([H, ·]|_{w=1})) = 2N
==========================================================
For the Heisenberg chain with Z-dephasing, the 2N purely-real eigenvalues
at Re = -2γ correspond to weight-1 operators that commute with H.

This script computes the kernel of the commutator map [H, ·] restricted
to the weight-1 Pauli sector, verifies dim(ker) = 2N, and identifies
the structure of the 2N kernel vectors.

Output: simulations/results/kommutator_kern.txt + console
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
# Pauli matrices
# ─────────────────────────────────────────────

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, X, Y, Z]
PAULI_LABELS = ['I', 'X', 'Y', 'Z']


def kron_chain(ops):
    """Kronecker product of a list of 2x2 matrices."""
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result


def heisenberg_chain(N, J=1.0):
    """Build Heisenberg XXX chain Hamiltonian, open boundary conditions."""
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
# Weight-1 Pauli string enumeration
# ─────────────────────────────────────────────

def enumerate_w1(N):
    """Enumerate all weight-1 Pauli strings for N qubits.
    Returns list of (label, d×d matrix) tuples.
    Weight-1: exactly one X or Y, rest I or Z.
    """
    basis = []
    for pos in range(N):
        for typ in [1, 2]:  # 1=X, 2=Y
            for z_cfg in range(2**(N - 1)):
                ops = []
                z_idx = 0
                label_chars = []
                for q in range(N):
                    if q == pos:
                        ops.append(PAULIS[typ])
                        label_chars.append(PAULI_LABELS[typ])
                    else:
                        bit = (z_cfg >> z_idx) & 1
                        ops.append(Z if bit else I2)
                        label_chars.append('Z' if bit else 'I')
                        z_idx += 1
                label = ''.join(label_chars)
                mat = kron_chain(ops)
                basis.append((label, mat))
    return basis


def xy_weight(label):
    """Count X and Y characters in a Pauli string label."""
    return sum(1 for c in label if c in ('X', 'Y'))


# ─────────────────────────────────────────────
# Commutator matrix and kernel
# ─────────────────────────────────────────────

def commutator_kernel(N, J=1.0):
    """Compute kernel of [H, ·] restricted to weight-1 sector.
    Returns (kernel_vectors, w1_basis, rank, singular_values).
    kernel_vectors: array of shape (dim_ker, len(w1_basis)), coefficients in w1 basis.
    """
    d = 2**N
    d2 = d * d
    H = heisenberg_chain(N, J)
    w1_basis = enumerate_w1(N)
    n_w1 = len(w1_basis)

    # Build commutator matrix: M[:, t] = vec([H, P_t])
    M = np.zeros((d2, n_w1), dtype=complex)
    for t, (label, P) in enumerate(w1_basis):
        comm = H @ P - P @ H
        M[:, t] = comm.ravel()  # row-major vectorization

    # SVD
    U, S, Vh = np.linalg.svd(M, full_matrices=False)

    # Kernel = right singular vectors with singular value < TOL
    rank = np.sum(S > TOL)
    dim_ker = n_w1 - rank
    kernel_vectors = Vh[rank:]  # rows of Vh corresponding to zero singular values

    return kernel_vectors, w1_basis, rank, S


# ─────────────────────────────────────────────
# Analysis
# ─────────────────────────────────────────────

def analyze_kernel_vector(coeffs, w1_basis, N, top_k=5):
    """Analyze a single kernel vector: dominant Pauli strings, structure."""
    indices = np.argsort(-np.abs(coeffs))
    result = []
    for i in indices[:top_k]:
        if np.abs(coeffs[i]) > TOL:
            label = w1_basis[i][0]
            result.append((label, coeffs[i]))
    return result


def verify_commutator_zero(H, v_op, label=""):
    """Verify [H, v_op] = 0 by computing Frobenius norm."""
    comm = H @ v_op - v_op @ H
    frob = np.linalg.norm(comm, 'fro')
    return frob


def reconstruct_operator(coeffs, w1_basis):
    """Reconstruct d×d operator from kernel vector coefficients."""
    d = w1_basis[0][1].shape[0]
    op = np.zeros((d, d), dtype=complex)
    for c, (label, mat) in zip(coeffs, w1_basis):
        if np.abs(c) > TOL:
            op += c * mat
    return op


def load_eigvec_csv(N):
    """Load eigenvector Pauli projections from C# export."""
    path = RESULTS_DIR / f"eigvec_at_minus_gamma_N{N}.csv"
    if not path.exists():
        return None

    data = {}  # eigvec_idx -> dict of pauli_label -> complex_coeff
    with open(path) as f:
        header = f.readline()
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
    return data


def compare_subspaces(kernel_vectors, w1_basis, csv_data, N):
    """Check if kernel vectors span the same subspace as eigenvectors from CSV."""
    if csv_data is None:
        return None, "no CSV data"

    n_eig = len(csv_data)
    n_ker = len(kernel_vectors)
    n_w1 = len(w1_basis)

    # Build label-to-index map for w1 basis
    label_to_idx = {w1_basis[i][0]: i for i in range(n_w1)}

    # Convert CSV eigenvectors to w1-basis coefficient vectors
    eig_vecs = np.zeros((n_eig, n_w1), dtype=complex)
    for ei, pauli_dict in csv_data.items():
        for label, coeff in pauli_dict.items():
            if label in label_to_idx:
                eig_vecs[ei, label_to_idx[label]] = coeff

    # Compute overlap matrix: kernel_vecs (n_ker × n_w1) vs eig_vecs (n_eig × n_w1)
    # Both should span the same n_ker-dimensional subspace
    # Check via rank of combined matrix
    combined = np.vstack([kernel_vectors, eig_vecs])
    _, S_comb, _ = np.linalg.svd(combined, full_matrices=False)
    rank_combined = np.sum(S_comb > TOL)

    return rank_combined, f"rank(kernel ∪ eigvec) = {rank_combined} (expect {n_ker})"


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

output_lines = []

def log(msg=""):
    print(msg)
    output_lines.append(msg)


log("=" * 70)
log("COMMUTATOR KERNEL ANALYSIS: dim(ker([H, ·]|_{w=1})) = 2N")
log("=" * 70)
log()

all_pass = True
all_kernels = {}

for N in range(2, 8):
    d = 2**N
    n_w1_expected = N * 2**N
    ker_dim_expected = 2 * N

    log("=" * 70)
    log(f"N={N}: d={d}, w=1 sector dim={n_w1_expected}, expected kernel dim={ker_dim_expected}")
    log("=" * 70)

    kernel_vectors, w1_basis, rank, sing_vals = commutator_kernel(N)
    dim_ker = len(kernel_vectors)
    all_kernels[N] = (kernel_vectors, w1_basis)

    status = "✓" if dim_ker == ker_dim_expected else "✗ MISMATCH"
    if dim_ker != ker_dim_expected:
        all_pass = False

    log(f"  w=1 basis size: {len(w1_basis)} (expected {n_w1_expected})")
    log(f"  Commutator matrix rank: {rank}")
    log(f"  Kernel dimension: {dim_ker} (expected {ker_dim_expected}) {status}")
    log(f"  Smallest nonzero SV: {sing_vals[rank-1]:.6e}")
    if dim_ker > 0:
        log(f"  Largest 'zero' SV:    {sing_vals[rank]:.6e}")
    log()

    # Verify [H, v] = 0 for each kernel vector
    H = heisenberg_chain(N)
    log(f"  Verification [H, v] = 0:")
    for ki in range(dim_ker):
        v_op = reconstruct_operator(kernel_vectors[ki], w1_basis)
        frob = verify_commutator_zero(H, v_op)
        status_v = "PASS" if frob < TOL else f"FAIL ({frob:.2e})"
        log(f"    Kernel vector #{ki}: ||[H,v]||_F = {frob:.2e} [{status_v}]")
    log()

    # Show kernel vectors in Pauli-string notation
    log(f"  Kernel vectors (top components):")
    for ki in range(dim_ker):
        components = analyze_kernel_vector(kernel_vectors[ki], w1_basis, N, top_k=6)
        parts = [f"{label}({c.real:+.4f}{c.imag:+.4f}j)" for label, c in components]
        log(f"    #{ki}: {', '.join(parts)}")
    log()

    # Analyze structure: which positions are active? Z-dressing pattern?
    log(f"  Structure analysis:")
    for ki in range(dim_ker):
        coeffs = kernel_vectors[ki]
        active_positions = set()
        active_types = set()
        for ci, (label, _) in zip(coeffs, w1_basis):
            if np.abs(ci) > TOL:
                for q, ch in enumerate(label):
                    if ch in ('X', 'Y'):
                        active_positions.add(q)
                        active_types.add(ch)

        n_nonzero = np.sum(np.abs(coeffs) > TOL)
        log(f"    #{ki}: {n_nonzero} nonzero Pauli strings, "
            f"active positions={sorted(active_positions)}, "
            f"types={sorted(active_types)}")
    log()

    # Compare with eigenvector CSVs
    csv_data = load_eigvec_csv(N)
    if csv_data is not None:
        rank_combined, msg = compare_subspaces(kernel_vectors, w1_basis, csv_data, N)
        match = rank_combined == dim_ker if rank_combined is not None else False
        status_csv = "✓ same subspace" if match else "✗ DIFFERENT"
        log(f"  CSV comparison: {msg} {status_csv}")
    else:
        log(f"  CSV comparison: no CSV file for N={N}")
    log()

# ─────────────────────────────────────────────
# Cross-N pattern analysis
# ─────────────────────────────────────────────

log("=" * 70)
log("CROSS-N PATTERN ANALYSIS")
log("=" * 70)
log()

# Check: are S_x and S_y always in the kernel?
log("Check: global spin operators S_x, S_y in kernel?")
for N in range(2, 8):
    kernel_vectors, w1_basis = all_kernels[N]
    n_w1 = len(w1_basis)

    # S_x = (1/2) Σ_j X_j ⊗ I_rest (only I on other positions)
    sx_coeffs = np.zeros(n_w1, dtype=complex)
    sy_coeffs = np.zeros(n_w1, dtype=complex)
    for i, (label, _) in enumerate(w1_basis):
        # Check if all non-active positions are I (not Z)
        n_x = label.count('X')
        n_y = label.count('Y')
        n_z = label.count('Z')
        if n_z == 0 and n_x == 1 and n_y == 0:
            sx_coeffs[i] = 0.5
        if n_z == 0 and n_y == 1 and n_x == 0:
            sy_coeffs[i] = 0.5

    # Project S_x onto kernel
    sx_proj = kernel_vectors @ sx_coeffs
    sy_proj = kernel_vectors @ sy_coeffs

    sx_in_ker = np.linalg.norm(sx_proj) > TOL
    sy_in_ker = np.linalg.norm(sy_proj) > TOL

    log(f"  N={N}: S_x in kernel: {'✓' if sx_in_ker else '✗'} "
        f"(||proj||={np.linalg.norm(sx_proj):.6f})  "
        f"S_y in kernel: {'✓' if sy_in_ker else '✗'} "
        f"(||proj||={np.linalg.norm(sy_proj):.6f})")
log()

# Check: how many kernel vectors involve ALL-I dressing vs Z-dressing?
log("Kernel vector Z-dressing analysis:")
for N in range(2, 8):
    kernel_vectors, w1_basis = all_kernels[N]
    dim_ker = len(kernel_vectors)

    for ki in range(dim_ker):
        coeffs = kernel_vectors[ki]
        # Count strings with all-I dressing (no Z)
        n_pure = 0  # strings like XII, IYI (no Z)
        n_dressed = 0  # strings like XZI, ZYZ (has Z)
        for ci, (label, _) in zip(coeffs, w1_basis):
            if np.abs(ci) > TOL:
                if 'Z' in label:
                    n_dressed += 1
                else:
                    n_pure += 1
        if N <= 4:
            log(f"  N={N} #{ki}: {n_pure} pure (no Z) + {n_dressed} Z-dressed = {n_pure + n_dressed} total")
log()

# Rank formula verification
log("Rank formula: rank([H,·]|_{w=1}) = N·2^N - 2N")
for N in range(2, 8):
    kernel_vectors, w1_basis = all_kernels[N]
    n_w1 = len(w1_basis)
    dim_ker = len(kernel_vectors)
    rank = n_w1 - dim_ker
    expected_rank = N * 2**N - 2 * N
    match = "✓" if rank == expected_rank else "✗"
    log(f"  N={N}: rank={rank}, N·2^N - 2N = {expected_rank} {match}")
log()

# Nonzero strings per kernel vector as function of N
log("Nonzero Pauli strings per kernel vector:")
for N in range(2, 8):
    kernel_vectors, w1_basis = all_kernels[N]
    counts = []
    for ki in range(len(kernel_vectors)):
        n_nz = np.sum(np.abs(kernel_vectors[ki]) > TOL)
        counts.append(n_nz)
    log(f"  N={N}: {counts}  (w=1 sector size: {len(w1_basis)})")
log()

# ─────────────────────────────────────────────
# Deep structure: decompose into X-type and Y-type
# ─────────────────────────────────────────────

log("=" * 70)
log("X/Y DECOMPOSITION OF KERNEL VECTORS")
log("=" * 70)
log()

for N in range(2, 6):  # detailed for small N
    kernel_vectors, w1_basis = all_kernels[N]
    log(f"N={N}: {len(kernel_vectors)} kernel vectors")

    for ki in range(len(kernel_vectors)):
        coeffs = kernel_vectors[ki]
        x_parts = []
        y_parts = []
        for ci, (label, _) in zip(coeffs, w1_basis):
            if np.abs(ci) > TOL:
                if 'X' in label:
                    x_parts.append((label, ci))
                else:
                    y_parts.append((label, ci))

        log(f"  #{ki}:")
        if x_parts:
            x_str = ", ".join(f"{l}({c.real:+.5f})" for l, c in x_parts[:8])
            log(f"    X-type ({len(x_parts)}): {x_str}")
        if y_parts:
            y_str = ", ".join(f"{l}({c.real:+.5f})" for l, c in y_parts[:8])
            log(f"    Y-type ({len(y_parts)}): {y_str}")

        # Check: are X and Y parts proportional? (same Z-dressing, different active type)
        if x_parts and y_parts:
            x_labels_stripped = set(l.replace('X', '*') for l, _ in x_parts)
            y_labels_stripped = set(l.replace('Y', '*') for l, _ in y_parts)
            overlap = x_labels_stripped & y_labels_stripped
            log(f"    X/Y pattern overlap: {len(overlap)}/{max(len(x_labels_stripped), len(y_labels_stripped))}")
    log()

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────

log("=" * 70)
log("SUMMARY")
log("=" * 70)
log(f"  dim(ker([H,·]|_{{w=1}})) = 2N confirmed: {'✓ ALL N' if all_pass else '✗ FAILED'}")
log()
log("  Results by N:")
for N in range(2, 8):
    kernel_vectors, _ = all_kernels[N]
    log(f"    N={N}: kernel dim = {len(kernel_vectors)}, expected = {2*N}")
log()

# Save output
out_path = RESULTS_DIR / "kommutator_kern.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(output_lines))
print(f"\n>>> Results saved to: {out_path}")
