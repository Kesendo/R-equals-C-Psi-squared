#!/usr/bin/env python3
"""F50 Q4: efficient n_XY=1 commutant sweep — bypasses dense Liouvillian.

The F50 formula `d_real(-2γ) = 2N + δ_G(N)` counts Pauli operators A satisfying
  L A = -2γ A  ⟺  [H_G, A] = 0  AND  n_XY(A) = 1
where n_XY counts X- and Y-letters in A. The Z-dephasing dissipator acts diagonally
on Pauli strings: D A = -2γ · n_XY(A) · A, so the eigenvalue λ = -2γ pins n_XY = 1.

The n_XY=1 subspace has dimension 2N · 2^(N-1):
  - choose one of 2N (site, letter) positions for the X or Y
  - choose Z-or-I for each of the remaining N-1 sites
The commutator [H_Heisenberg, ·] maps this subspace into the n_XY ∈ {1, 3} sub-
space (Heisenberg ZZ bond can convert Z·I or I·Z into Y·X or X·Y, raising n_XY
by 2). We build the matrix of [H_G, ·] restricted to n_XY=1 (output basis is the
union n_XY ∈ {1, 3} subspace), then rank → kernel dim → δ_G(N).

  N | n_XY=1 dim | output basis bound
  --|------------|--------------------
  3 |         12 |                  ~50
  5 |         80 |                ~700
  7 |        448 |               ~5000
  8 |       2048 |              ~17000
 10 |      10240 |              ~110000

For Petersen N=10 the matrix is ~10240 × 110000 ≈ 1 G floats ≈ 8 GB — manageable
on a 32 GB machine via SciPy sparse; this script uses dense numpy for clarity at
the cost of memory. Fall back to N ≤ 8 if Petersen exceeds available RAM.

Open per PROOF_WEIGHT1_DEGENERACY.md Q4 (2026-05-17 evening):
"Is K_3 N=3 the ONLY weight-1 anomaly, or are there other small-graph cases?"
"Petersen graph at N=10, hypercube Q_3 at N=8 UNTESTED per the doc."
"""
from __future__ import annotations

import sys
import numpy as np
from itertools import combinations, product
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import svds

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Pauli string encoding + commutator algebra
# ---------------------------------------------------------------------------

# We encode a Pauli string as a tuple of letters ('I', 'X', 'Y', 'Z') of length N.
# Multiplication σ_a · σ_b = c_ab · σ_{ab} where c_ab is a phase factor and σ_{ab}
# is the resulting Pauli letter.

PAULI_TABLE = {
    ("I", "I"): ( 1, "I"), ("I", "X"): ( 1, "X"), ("I", "Y"): ( 1, "Y"), ("I", "Z"): ( 1, "Z"),
    ("X", "I"): ( 1, "X"), ("X", "X"): ( 1, "I"), ("X", "Y"): ( 1j, "Z"), ("X", "Z"): (-1j, "Y"),
    ("Y", "I"): ( 1, "Y"), ("Y", "X"): (-1j, "Z"), ("Y", "Y"): ( 1, "I"), ("Y", "Z"): ( 1j, "X"),
    ("Z", "I"): ( 1, "Z"), ("Z", "X"): ( 1j, "Y"), ("Z", "Y"): (-1j, "X"), ("Z", "Z"): ( 1, "I"),
}


def pauli_string_mul(s1, s2):
    """Multiply two Pauli strings. Returns (phase, result_string)."""
    if len(s1) != len(s2):
        raise ValueError(f"length mismatch: {len(s1)} vs {len(s2)}")
    phase = 1 + 0j
    out = []
    for a, b in zip(s1, s2):
        p, c = PAULI_TABLE[(a, b)]
        phase *= p
        out.append(c)
    return phase, tuple(out)


def commutator_strings(s1, s2):
    """[s1, s2] = s1 s2 - s2 s1. Returns (coeff, result_string) or None if zero.
    For two Pauli strings, the commutator is either 0 or 2·(phase) · result."""
    p1, r1 = pauli_string_mul(s1, s2)
    p2, r2 = pauli_string_mul(s2, s1)
    if r1 != r2:
        # Shouldn't happen for Pauli products on same N: results are unique up to phase.
        raise RuntimeError(f"unexpected: {r1} vs {r2}")
    coeff = p1 - p2
    if abs(coeff) < 1e-14:
        return None
    return coeff, r1


def n_xy(string):
    return sum(1 for c in string if c in ("X", "Y"))


# ---------------------------------------------------------------------------
# Heisenberg bond operators + Liouvillian commutator action
# ---------------------------------------------------------------------------

def bond_string(N, i, j, a, b):
    """Single Pauli string for σ_a^i σ_b^j (Heisenberg bond term: (XX + YY + ZZ) / 4)."""
    s = ["I"] * N
    s[i] = a
    s[j] = b
    return tuple(s)


def heisenberg_terms(N, bonds, J=1.0):
    """List of (coeff, Pauli string) terms in H_G = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1})."""
    terms = []
    for (i, j) in bonds:
        for letter in ("X", "Y", "Z"):
            terms.append((J / 4.0, bond_string(N, i, j, letter, letter)))
    return terms


def commutator_with_H(A, H_terms):
    """[H, A] = sum over H terms of [h, A]. Returns dict {string: complex coeff}."""
    out = {}
    for (h_coeff, h_string) in H_terms:
        c = commutator_strings(h_string, A)
        if c is None:
            continue
        coeff, result = c
        full = h_coeff * coeff
        if result in out:
            out[result] += full
        else:
            out[result] = full
    return {k: v for k, v in out.items() if abs(v) > 1e-12}


# ---------------------------------------------------------------------------
# n_XY = 1 subspace enumeration + rank computation
# ---------------------------------------------------------------------------

def enumerate_nxy_eq_1(N):
    """Yield all Pauli strings of length N with exactly 1 X or Y letter (and
    remaining letters in {I, Z}). Total count = 2N · 2^(N-1)."""
    for l in range(N):  # site of the X/Y
        for letter in ("X", "Y"):
            other_sites = [k for k in range(N) if k != l]
            for zi_mask in range(1 << len(other_sites)):
                s = ["I"] * N
                s[l] = letter
                for idx, other in enumerate(other_sites):
                    if (zi_mask >> idx) & 1:
                        s[other] = "Z"
                yield tuple(s)


def commutant_dim_nxy1(bonds, N, verbose=False):
    """dim ker[H_G, ·] restricted to the n_XY = 1 subspace.

    Builds a sparse matrix mapping n_XY=1 basis ops to their commutator-image
    expressed in the union n_XY ∈ {1, 3} basis."""
    H_terms = heisenberg_terms(N, bonds)
    input_basis = list(enumerate_nxy_eq_1(N))
    input_dim = len(input_basis)
    input_idx = {s: i for i, s in enumerate(input_basis)}

    # Build sparse [H, ·] matrix; output rows accumulate dynamically.
    output_idx = dict(input_idx)  # n_XY=1 ops first; n_XY=3 appended on the fly
    # We use a transposed (rows = input, cols = output) layout, then take rank
    # of the transposed matrix. Use lil_matrix for incremental build.

    rows, cols, data = [], [], []
    for i, A in enumerate(input_basis):
        commut = commutator_with_H(A, H_terms)
        for B, coeff in commut.items():
            if B not in output_idx:
                output_idx[B] = len(output_idx)
            j = output_idx[B]
            rows.append(i)
            cols.append(j)
            data.append(coeff)

    output_dim = len(output_idx)
    if verbose:
        print(f"    [H, ·] matrix: input dim {input_dim}, output dim {output_dim}, nnz {len(data)}")

    M = csr_matrix((data, (rows, cols)), shape=(input_dim, output_dim), dtype=complex)

    # For rank: SVD-based via numpy on dense. At small N this is fine; at N=10
    # input_dim = 10240, output_dim could be ~50k → 4 GB dense. Too big.
    # Better: use scipy sparse SVD via svds (gives only top-k singular values, not
    # full rank). For rank computation, we need the count of σ_i > tol. The
    # robust approach: dense numpy on the rectangular matrix (input × output)
    # or its Gram matrix M·M† (input × input, smaller).
    #
    # M·M† has the same non-zero singular values squared. dim input × input
    # = 10240² × 16 bytes = 1.7 GB. Doable.
    #
    # For very large N, replace with iterative rank estimate (LSQR-based).

    if input_dim <= 16384:
        gram = (M @ M.conjugate().T).toarray()
        eigs = np.linalg.eigvalsh(gram)
        rank = int(np.sum(eigs > 1e-12 * eigs.max() if eigs.max() > 0 else 1e-12))
    else:
        raise NotImplementedError(f"input_dim {input_dim} too large; need iterative rank")

    kernel = input_dim - rank
    return kernel, input_dim, rank


# ---------------------------------------------------------------------------
# Graph constructors
# ---------------------------------------------------------------------------

def chain(N):
    return [(l, l + 1) for l in range(N - 1)]


def ring(N):
    return [(l, (l + 1) % N) for l in range(N)]


def star(N):
    return [(0, l) for l in range(1, N)]


def complete(N):
    return list(combinations(range(N), 2))


def petersen():
    """Petersen graph: 10 vertices, 15 edges."""
    outer = [(i, (i + 1) % 5) for i in range(5)]
    inner = [(5 + i, 5 + ((i + 2) % 5)) for i in range(5)]
    spokes = [(i, 5 + i) for i in range(5)]
    return outer + inner + spokes


def q3_cube():
    """3-cube Q_3: 8 vertices, 12 edges (one-bit-flip neighbors)."""
    edges = []
    for a in range(8):
        for b in range(a + 1, 8):
            if bin(a ^ b).count("1") == 1:
                edges.append((a, b))
    return edges


def k33():
    """K_{3,3} bipartite: 6 vertices, 9 edges."""
    return [(i, j) for i in range(3) for j in range(3, 6)]


def mobius_kantor_M8():
    """Cayley graph of Z_8 with generators ±1, ±3: edge-transitive, Aut size 96."""
    edges = set()
    for i in range(8):
        for d in [1, 3]:
            edges.add(tuple(sorted([i, (i + d) % 8])))
    return sorted(edges)


GRAPH_SUITE = [
    # (name, N, expected 2N, bonds factory)
    ("chain N=3", 3, 6, lambda: chain(3)),
    ("K_3 N=3 (known anomaly: 2N+2 = 8)", 3, 6, lambda: complete(3)),
    ("chain N=4", 4, 8, lambda: chain(4)),
    ("K_4", 4, 8, lambda: complete(4)),
    ("chain N=5", 5, 10, lambda: chain(5)),
    ("K_5", 5, 10, lambda: complete(5)),
    ("chain N=6", 6, 12, lambda: chain(6)),
    ("K_{3,3} (bipartite, Aut 72)", 6, 12, lambda: k33()),
    ("Q_3 cube (N=8, edge-transitive, Aut 48)", 8, 16, lambda: q3_cube()),
    ("Möbius-Kantor M_8 (N=8, edge-transitive, Aut 96)", 8, 16, lambda: mobius_kantor_M8()),
    ("Petersen (N=10, vertex+edge transitive, Aut 120)", 10, 20, lambda: petersen()),
]


def main():
    print("=" * 80)
    print("F50 Q4: n_XY=1 commutant sweep — answers 'is K_3 N=3 the only anomaly?'")
    print("=" * 80)
    print()
    print(f"  {'graph':<54} {'edges':>5} {'2N':>4} {'dim':>5} {'δ':>5}")
    print(f"  {'-' * 54} {'-' * 5} {'-' * 4} {'-' * 5} {'-' * 5}")

    for name, N, expected, factory in GRAPH_SUITE:
        bonds = factory()
        try:
            kernel, input_dim, rank = commutant_dim_nxy1(bonds, N)
        except NotImplementedError as e:
            print(f"  {name:<54}   SKIP (too large: {e})")
            continue
        delta = kernel - expected
        mark = "OK" if delta == 0 else (f"+{delta}" if delta > 0 else f"{delta}")
        print(f"  {name:<54} {len(bonds):>5} {expected:>4} {kernel:>5} {mark:>5}")

    print()
    print("=" * 80)
    print("Interpretation")
    print("=" * 80)
    print("  δ = (n_XY=1 commutant dim) − 2N")
    print("  δ = 0: F50 universal formula d_real(-2γ) = 2N holds")
    print("  δ > 0: NEW anomaly beyond K_3 N=3")
    print()
    print("  Per PROOF_WEIGHT1_DEGENERACY.md Q4: Petersen + Q_3 cube + K_{3,3} are")
    print("  the high-symmetry untested cases. If all show δ = 0, K_3 N=3 remains")
    print("  the unique anomaly: the small-N face of the universal central-weight")
    print("  excess pattern at the N=3-only coincidence (central weight = 1).")


if __name__ == "__main__":
    main()
