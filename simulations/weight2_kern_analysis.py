"""
Weight-2 commutator kernel analysis
====================================
Computes ker([H, ·]|_{w=2}) for Chain, Star, Complete topologies.
Decomposes kernel vectors into type classes (XX, YY, XY+YX, XY-YX),
tests SWAP eigenvalues, and counts orbits.

Output: simulations/results/weight2_kern_analysis.txt
"""

import numpy as np
from itertools import combinations
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
TOL = 1e-10

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_hamiltonian(N, bonds):
    """bonds: list of (i, j) pairs."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        for P in [X, Y, Z]:
            ops = [I2] * N
            ops[a] = P
            ops[b] = P
            H += kron_chain(ops)
    return H


def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]

def star_bonds(N):
    return [(0, i) for i in range(1, N)]

def ring_bonds(N):
    return [(i, (i + 1) % N) for i in range(N)]

def complete_bonds(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


# ─────────────────────────────────────────────
# Weight-2 basis
# ─────────────────────────────────────────────

def enumerate_w2(N):
    """All weight-2 Pauli strings: exactly 2 positions with X or Y."""
    basis = []
    for pos_pair in combinations(range(N), 2):
        p0, p1 = pos_pair
        for t0 in ['X', 'Y']:
            for t1 in ['X', 'Y']:
                passive = [q for q in range(N) if q not in pos_pair]
                for z_cfg in range(2**len(passive)):
                    chars = ['I'] * N
                    chars[p0] = t0
                    chars[p1] = t1
                    for zi, pq in enumerate(passive):
                        if (z_cfg >> zi) & 1:
                            chars[pq] = 'Z'
                    basis.append(''.join(chars))
    return basis


def type_class(label):
    """Extract the active type pair from a w=2 label."""
    active = [c for c in label if c in ('X', 'Y')]
    return ''.join(active)  # XX, XY, YX, YY


def z_count(label):
    return label.count('Z')


def active_positions(label):
    return tuple(i for i, c in enumerate(label) if c in ('X', 'Y'))


def apply_swap_label(label, s):
    chars = list(label)
    chars[s], chars[s + 1] = chars[s + 1], chars[s]
    return ''.join(chars)


def build_matrix(label):
    return kron_chain([PAULI[c] for c in label])


# ─────────────────────────────────────────────
# Commutator kernel
# ─────────────────────────────────────────────

def commutator_kernel_w2(N, bonds):
    d = 2**N
    d2 = d * d
    H = build_hamiltonian(N, bonds)
    w2_basis = enumerate_w2(N)
    n_w2 = len(w2_basis)

    M = np.zeros((d2, n_w2), dtype=complex)
    for t, label in enumerate(w2_basis):
        P = build_matrix(label)
        comm = H @ P - P @ H
        M[:, t] = comm.ravel()

    U, S, Vh = np.linalg.svd(M, full_matrices=False)
    rank = np.sum(S > TOL)
    kernel = Vh[rank:]
    return kernel, w2_basis, rank, S


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

out = []
def log(msg=""):
    print(msg)
    out.append(msg)


log("=" * 75)
log("WEIGHT-2 COMMUTATOR KERNEL ANALYSIS")
log("=" * 75)
log()

# Expected d_real(2) for Chain from DEGENERACY_PALINDROME
expected_chain = {3: 6, 4: 14, 5: 14, 6: 19}
# From topology comparison
expected_star = {4: 16, 5: 30}
expected_complete = {4: 36, 5: 54}

all_kernels = {}

for topo_name, bond_fn, expected, n_range in [
    ("Chain", chain_bonds, expected_chain, range(3, 7)),
    ("Star", star_bonds, expected_star, [4, 5]),
    ("Complete", complete_bonds, expected_complete, [4, 5]),
]:
    log("=" * 75)
    log(f"TOPOLOGY: {topo_name}")
    log("=" * 75)

    for N in n_range:
        bonds = bond_fn(N)
        n_w2 = len(enumerate_w2(N))

        log(f"\n  N={N}: w=2 sector dim = {n_w2}, bonds = {len(bonds)}")

        kernel, w2_basis, rank, sv = commutator_kernel_w2(N, bonds)
        dim_ker = len(kernel)
        label_to_idx = {l: i for i, l in enumerate(w2_basis)}

        exp = expected.get(N)
        match = f"✓ (expected {exp})" if exp and dim_ker == exp else (
            f"✗ MISMATCH (expected {exp})" if exp else "")
        log(f"  Kernel dim: {dim_ker} {match}")
        log(f"  Rank: {rank}, gap: {sv[rank-1]:.2e} / {sv[rank]:.2e}" if dim_ker > 0
            else f"  Rank: {rank}")

        # Store for Chain detailed analysis
        if topo_name == "Chain":
            all_kernels[(topo_name, N)] = (kernel, w2_basis, bonds)

        # ─── Verify [H, v] = 0 ───
        H = build_hamiltonian(N, bonds)
        max_res = 0
        for ki in range(dim_ker):
            op = np.zeros((2**N, 2**N), dtype=complex)
            for ci, label in zip(kernel[ki], w2_basis):
                if abs(ci) > TOL:
                    op += ci * build_matrix(label)
            res = np.linalg.norm(H @ op - op @ H, 'fro')
            max_res = max(max_res, res)
        log(f"  Verification ||[H,v]||: max = {max_res:.2e} {'PASS' if max_res < 1e-8 else 'FAIL'}")

        # ─── SWAP eigenvalue analysis ───
        if topo_name == "Chain" and N <= 6:
            log(f"\n  SWAP eigenvalue analysis:")
            n_bonds = len(bonds)
            swap_perms = []
            for s, _ in bonds:
                perm = {}
                for i, label in enumerate(w2_basis):
                    new_label = apply_swap_label(label, s)
                    perm[i] = label_to_idx[new_label]
                swap_perms.append(perm)

            # Build permutation matrices
            swap_mats = []
            for perm in swap_perms:
                P = np.zeros((n_w2, n_w2))
                for old, new in perm.items():
                    P[new, old] = 1
                swap_mats.append(P)

            for ki in range(dim_ker):
                v = kernel[ki]
                eigenvals = []
                for si, P in enumerate(swap_mats):
                    sv_test = P @ v
                    # Check if sv = +v or sv = -v
                    if np.linalg.norm(sv_test - v) < TOL:
                        eigenvals.append("+1")
                    elif np.linalg.norm(sv_test + v) < TOL:
                        eigenvals.append("-1")
                    else:
                        ratio = sv_test[np.argmax(np.abs(v))] / v[np.argmax(np.abs(v))]
                        eigenvals.append(f"~{ratio.real:+.3f}")

                log(f"    #{ki}: SWAP eigenvalues = [{', '.join(eigenvals)}]")

            # Count by type
            n_trivial = 0  # all +1
            n_alternating = 0  # all -1
            n_mixed = 0
            for ki in range(dim_ker):
                v = kernel[ki]
                all_plus = all(np.linalg.norm(P @ v - v) < TOL for P in swap_mats)
                all_minus = all(np.linalg.norm(P @ v + v) < TOL for P in swap_mats)
                if all_plus:
                    n_trivial += 1
                elif all_minus:
                    n_alternating += 1
                else:
                    n_mixed += 1

            log(f"\n    Summary: {n_trivial} trivial (+1), {n_alternating} alternating (-1), {n_mixed} mixed")

        # ─── Type class decomposition ───
        if topo_name == "Chain" and N <= 6:
            log(f"\n  Type class decomposition:")
            for ki in range(min(dim_ker, 8)):  # show up to 8
                v = kernel[ki]
                counts = {'XX': 0, 'XY': 0, 'YX': 0, 'YY': 0}
                for ci, label in zip(v, w2_basis):
                    if abs(ci) > TOL:
                        tc = type_class(label)
                        counts[tc] += 1

                # Check symmetry pattern
                sym_count = 0
                anti_count = 0
                for ci, label in zip(v, w2_basis):
                    if abs(ci) > TOL:
                        tc = type_class(label)
                        # Find the partner with swapped types
                        ap = active_positions(label)
                        partner_chars = list(label)
                        partner_chars[ap[0]], partner_chars[ap[1]] = partner_chars[ap[1]], partner_chars[ap[0]]
                        # Swap the active types
                        if tc == 'XY':
                            partner_label = label  # need to rebuild
                            pc = list(label)
                            pc[ap[0]] = 'Y'
                            pc[ap[1]] = 'X'
                            partner_label = ''.join(pc)
                        elif tc == 'YX':
                            pc = list(label)
                            pc[ap[0]] = 'X'
                            pc[ap[1]] = 'Y'
                            partner_label = ''.join(pc)
                        else:
                            partner_label = None

                        if partner_label and partner_label in label_to_idx:
                            pi = label_to_idx[partner_label]
                            cj = v[pi]
                            if abs(ci - cj) < TOL:
                                sym_count += 1
                            elif abs(ci + cj) < TOL:
                                anti_count += 1

                log(f"    #{ki}: XX={counts['XX']} XY={counts['XY']} YX={counts['YX']} YY={counts['YY']}  "
                    f"sym_pairs={sym_count//2} anti_pairs={anti_count//2}")

    log()

# ─────────────────────────────────────────────
# Cross-topology comparison at N=4
# ─────────────────────────────────────────────

log("=" * 75)
log("CROSS-TOPOLOGY SUMMARY")
log("=" * 75)
log()

for N in [4, 5]:
    log(f"N={N}:")
    for topo_name, bond_fn, expected in [
        ("Chain", chain_bonds, expected_chain),
        ("Star", star_bonds, expected_star),
        ("Complete", complete_bonds, expected_complete),
    ]:
        if N not in expected:
            continue
        bonds = bond_fn(N)
        kernel, w2_basis, rank, _ = commutator_kernel_w2(N, bonds)
        log(f"  {topo_name:10s}: dim(ker) = {len(kernel):3d}  (expected {expected[N]})  "
            f"bonds = {len(bonds)}")
    log()

# ─────────────────────────────────────────────
# Orbit counting for Chain
# ─────────────────────────────────────────────

log("=" * 75)
log("ORBIT STRUCTURE (Chain)")
log("=" * 75)
log()

for N in range(3, 7):
    kernel, w2_basis, bonds = all_kernels[("Chain", N)]
    label_to_idx = {l: i for i, l in enumerate(w2_basis)}
    n_bonds = len(bonds)

    # Build SWAP permutation matrices for chain bonds
    swap_mats = []
    for s in range(n_bonds):
        n_w2 = len(w2_basis)
        P = np.zeros((n_w2, n_w2))
        for i, label in enumerate(w2_basis):
            new_label = apply_swap_label(label, s)
            P[label_to_idx[new_label], i] = 1
        swap_mats.append(P)

    # For each kernel vector: classify as trivial/alternating/mixed
    n_w2 = len(w2_basis)
    trivial_vecs = []
    alt_vecs = []
    mixed_vecs = []

    for ki in range(len(kernel)):
        v = kernel[ki]
        is_trivial = all(np.linalg.norm(P @ v - v) < TOL for P in swap_mats)
        is_alt = all(np.linalg.norm(P @ v + v) < TOL for P in swap_mats)
        if is_trivial:
            trivial_vecs.append(ki)
        elif is_alt:
            alt_vecs.append(ki)
        else:
            mixed_vecs.append(ki)

    # Count Z-count range for w=2
    z_range = N - 2 + 1  # 0 to N-2

    # Symmetric types: XX, YY, XY+YX → 3 types × z_range
    sym_expected = 3 * z_range
    # Antisymmetric type: XY-YX → need to count orbits

    log(f"N={N}: kernel dim = {len(kernel)}")
    log(f"  Trivial (all SWAPs = +1): {len(trivial_vecs)}")
    log(f"  Alternating (all SWAPs = -1): {len(alt_vecs)}")
    log(f"  Mixed: {len(mixed_vecs)}")
    log(f"  Naive symmetric count (3 × {z_range}): {sym_expected}")
    log(f"  Excess over symmetric: {len(kernel) - sym_expected}")
    log(f"  Alternating count: {len(alt_vecs)}")
    log(f"  Check: trivial + alternating = kernel? {len(trivial_vecs) + len(alt_vecs) == len(kernel)}")
    log()

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────

log("=" * 75)
log("FINAL SUMMARY")
log("=" * 75)
log()
log("Chain d_real(2) decomposition:")
log(f"{'N':>3s} {'dim_ker':>8s} {'trivial':>8s} {'altern.':>8s} {'mixed':>8s} {'3(N-2)+1':>9s}")
for N in range(3, 7):
    kernel, w2_basis, bonds = all_kernels[("Chain", N)]
    label_to_idx = {l: i for i, l in enumerate(w2_basis)}
    swap_mats_local = []
    for s in range(len(bonds)):
        n_w2 = len(w2_basis)
        P = np.zeros((n_w2, n_w2))
        for i, label in enumerate(w2_basis):
            new_label = apply_swap_label(label, s)
            P[label_to_idx[new_label], i] = 1
        swap_mats_local.append(P)

    nt = sum(1 for ki in range(len(kernel))
             if all(np.linalg.norm(P @ kernel[ki] - kernel[ki]) < TOL for P in swap_mats_local))
    na = sum(1 for ki in range(len(kernel))
             if all(np.linalg.norm(P @ kernel[ki] + kernel[ki]) < TOL for P in swap_mats_local))
    nm = len(kernel) - nt - na
    log(f"{N:3d} {len(kernel):8d} {nt:8d} {na:8d} {nm:8d} {3*(N-2)+1:9d}")

# Save
out_path = RESULTS_DIR / "weight2_kern_analysis.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
print(f"\n>>> Results saved to: {out_path}")
