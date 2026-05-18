"""Verify the F1 H-block residual closed form on arbitrary graphs.

Closes the last F1 OpenQuestion ("general topology beyond chain/ring/star/K_N")
by extending verification from the chain/ring/star/K_N suite at N=4, 5 to:

  (1) anchor c_H at N=2 single-bond non-truly Hamiltonian (XX+YZ, Π²-mixed)
  (2) named graphs (path, cycle, star, K_N, K_{2,N−2}) at N=5, 6
  (3) random connected Erdős-Rényi graphs at N=5, 6 (30 per N across p ∈ {0.3, 0.5, 0.7})
  (4) disconnected components (two disjoint 3-site chains at N=6)
  (5) weighted edges (chain at N=4 with non-uniform per-bond couplings)
  (6) single-body class spot check at N=5 chain (F = (D2/2)·4^(N−2))

The math is already settled analytically: PROOF_CROSS_TERM_FORMULA Lemma 3 Corollary
shows that different bonds have disjoint Pauli-string transition supports for any
graph topology. The (B, D2) parameterisation is therefore universal: connected or
disconnected, weighted or unweighted, any graph. This script is the verification-
record extension that closes the OpenQuestion: it confirms the closed form

    ‖M(N, G)‖²_F = c_H · B(G) · 4^(N − 2)            (main class)
    ‖M(N, G)‖²_F = c_H · (D2(G) / 2) · 4^(N − 2)     (single-body class)

bit-exactly for arbitrary graphs at the verified-N range.

Hamiltonian choice: Heisenberg (XX+YY+ZZ) is "truly" (Π²-even), so the F1 residual
vanishes identically and c_H = 0 — useless for verifying the scaling. We use the
canonical non-truly bilinear XX+YZ per bond, matching PalindromeResidualScalingClaim's
anchor in compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs (XX+YZ at
N=2 gives c_H). For weighted edges (section 5), the bond count B is interpreted as
Σ_b J²_b / J²_ref because each bond's Hamiltonian H_b contributes ‖H_b‖²_F ∝ J²_b.
For the single-body class (section 6), we use the bond-bilinear I·Y on each chain
bond (IY + YI), matching the convention from OPERATOR_RIGIDITY_ACROSS_CUSP.md
section "Algebraic origin" (single-body bilinears Iσ, σI).
"""
from __future__ import annotations

import sys

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

import framework as fw  # noqa: E402
from framework.lindblad import (  # noqa: E402
    lindbladian_z_dephasing,
    palindrome_residual,
)


# --------------------------------------------------------------------------- #
# Hamiltonian builders                                                        #
# --------------------------------------------------------------------------- #


def xx_plus_yz_graph_h(N: int, bonds: list[tuple[int, int]],
                       weights: list[float] | None = None) -> np.ndarray:
    """Build H = Σ_b J_b · (X_i X_j + Y_i Z_j) for the given bond list.

    The XX+YZ per-bond Hamiltonian is the canonical non-truly choice used by the
    PalindromeResidualScalingClaim.Verify anchor in
    compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs; it gives a
    non-zero c_H so the scaling test is meaningful (Heisenberg XX+YY+ZZ is truly
    and would give c_H = 0).
    """
    if weights is None:
        weights = [1.0] * len(bonds)
    if len(weights) != len(bonds):
        raise ValueError(
            f"weights length {len(weights)} must match bonds length {len(bonds)}"
        )
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for J_b, (i, j) in zip(weights, bonds):
        H = H + J_b * fw.site_op(N, i, "X") @ fw.site_op(N, j, "X")
        H = H + J_b * fw.site_op(N, i, "Y") @ fw.site_op(N, j, "Z")
    return H


def single_body_iy_plus_yi_graph_h(N: int, bonds: list[tuple[int, int]]) -> np.ndarray:
    """Build single-body bond-bilinear H = Σ_b (I_i Y_j + Y_i I_j) on the bond set.

    Single-body bilinears (Iσ, σI) generate the operator J·(σ_0 + Σ_{interior} 2σ
    + σ_{N−1}) on a chain: middle sites doubled because they appear in two bonds.
    On a general graph G this becomes Σ_i deg(i)·Y_i, so the per-site effective
    coefficient is deg(i) and the single-body D2 invariant captures the squared
    Frobenius norm via Σ_i deg(i)². Matches the convention from
    OPERATOR_RIGIDITY_ACROSS_CUSP.md "Algebraic origin" section.
    """
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        H = H + fw.site_op(N, i, "I") @ fw.site_op(N, j, "Y")  # = Y_j alone
        H = H + fw.site_op(N, i, "Y") @ fw.site_op(N, j, "I")  # = Y_i alone
    return H


def m_norm_squared(N: int, H: np.ndarray, gamma: list[float]) -> float:
    """‖M‖²_F for M = Π·L·Π⁻¹ + L + 2Σγ·I under uniform Z-dephasing."""
    L = lindbladian_z_dephasing(H, gamma)
    M = palindrome_residual(L, sum(gamma), N)
    return float(np.real(np.trace(M.conj().T @ M)))


def bond_count_and_d2(N: int, bonds: list[tuple[int, int]]) -> tuple[int, int]:
    """Compute (B, D2) for an unweighted bond list: B = len(bonds),
    D2 = Σ_i deg(i)² counting both endpoints."""
    deg = [0] * N
    for (i, j) in bonds:
        deg[i] += 1
        deg[j] += 1
    return len(bonds), sum(d * d for d in deg)


# --------------------------------------------------------------------------- #
# Anchor                                                                      #
# --------------------------------------------------------------------------- #


def anchor_c_h_main() -> float:
    """Extract c_H (main class) from the smallest non-truly bilinear at N=2.

    XX+YZ per bond at N=2 (one bond), uniform γ=0.1. F(2, chain, main) = B·4^(N−2)
    = 1·1 = 1, so c_H = ‖M(N=2)‖² directly. Matches the c_H = 128 anchor for
    XX+YZ documented in OPERATOR_RIGIDITY_ACROSS_CUSP.md.
    """
    N = 2
    bonds = [(0, 1)]
    H = xx_plus_yz_graph_h(N, bonds)
    cH = m_norm_squared(N, H, [0.1] * N)
    print(f"  anchor c_H (main, XX+YZ, γ=0.1, J=1) from N=2 single bond: c_H = {cH:.10f}")
    return cH


def anchor_c_h_single_body() -> float:
    """Extract c_H (single-body class) from N=2 single-bond IY+YI.

    F(2, chain, single_body) = (D2/2)·4^(N−2) = (2/2)·1 = 1, so
    c_H_single = ‖M(N=2)‖² directly.
    """
    N = 2
    bonds = [(0, 1)]
    H = single_body_iy_plus_yi_graph_h(N, bonds)
    cH = m_norm_squared(N, H, [0.1] * N)
    print(f"  anchor c_H_single (IY+YI, γ=0.1) from N=2 single bond: c_H_single = {cH:.10f}")
    return cH


# --------------------------------------------------------------------------- #
# Graph generators                                                            #
# --------------------------------------------------------------------------- #


def path_bonds(N: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(N - 1)]


def cycle_bonds(N: int) -> list[tuple[int, int]]:
    return [(i, (i + 1) % N) for i in range(N)]


def star_bonds(N: int) -> list[tuple[int, int]]:
    return [(0, i) for i in range(1, N)]


def complete_bonds(N: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


def complete_bipartite_bonds(N: int, parts: tuple[int, int]) -> list[tuple[int, int]]:
    """K_{a,b} with a + b = N; left part = sites 0..a−1, right part = sites a..N−1."""
    a, b = parts
    if a + b != N:
        raise ValueError(f"K_{{a,b}} needs a+b = N; got a+b = {a+b}, N = {N}")
    return [(i, j) for i in range(a) for j in range(a, N)]


def is_connected(N: int, bonds: list[tuple[int, int]]) -> bool:
    """Trivial union-find connectedness for small N."""
    parent = list(range(N))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (i, j) in bonds:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    return len({find(i) for i in range(N)}) == 1


def random_connected_erdos_renyi(N: int, p: float, rng: np.random.RandomState,
                                  max_tries: int = 200) -> list[tuple[int, int]]:
    """Sample G(N, p) and reject until connected. p is per-edge inclusion probability."""
    for _ in range(max_tries):
        bonds = [(i, j) for i in range(N) for j in range(i + 1, N)
                 if rng.uniform() < p]
        if len(bonds) >= N - 1 and is_connected(N, bonds):
            return bonds
    raise RuntimeError(
        f"Failed to sample connected G({N}, {p}) after {max_tries} tries"
    )


# --------------------------------------------------------------------------- #
# Sections                                                                    #
# --------------------------------------------------------------------------- #


def section_2_named_graphs(cH: float) -> None:
    """Named graphs (path, cycle, star, K_N, K_{2,N−2}) at N=5, 6: bit-exact match
    against c_H · B · 4^(N−2)."""
    print("\n(2) Named graphs at N=5, 6: path, cycle, star, K_N, K_{2,N−2}.")
    print("    Closed form: ‖M‖² = c_H · B · 4^(N−2). Bit-exact assertion.")
    cases: list[tuple[int, str, list[tuple[int, int]]]] = []
    for N in (5, 6):
        cases.extend([
            (N, "path",      path_bonds(N)),
            (N, "cycle",     cycle_bonds(N)),
            (N, "star",      star_bonds(N)),
            (N, "K_N",       complete_bonds(N)),
            (N, "K_{2,N−2}", complete_bipartite_bonds(N, (2, N - 2))),
        ])
    for N, label, bonds in cases:
        B, _ = bond_count_and_d2(N, bonds)
        predicted = cH * B * 4 ** (N - 2)
        H = xx_plus_yz_graph_h(N, bonds)
        observed = m_norm_squared(N, H, [0.1] * N)
        diff = abs(observed - predicted)
        rel = diff / predicted if predicted else diff
        print(f"  N={N} {label:12s} B={B:2d}: observed = {observed:.6f}  "
              f"predicted = {predicted:.6f}  |Δ| = {diff:.3e}  (rel = {rel:.2e})")
        # Bit-exact on small dim; allow ~4^N · ε rounding noise via 1e-9 relative.
        assert rel < 1e-9, (
            f"N={N} {label}: rel = {rel:.3e} exceeds 1e-9; "
            f"observed = {observed}, predicted = {predicted}"
        )


def section_3_random_graphs(cH: float) -> None:
    """30 random connected Erdős-Rényi graphs at N=5, 6, 10 per p ∈ {0.3, 0.5, 0.7}."""
    print("\n(3) Random connected Erdős-Rényi graphs at N=5, 6 (30 per N).")
    print("    p ∈ {0.3, 0.5, 0.7}, 10 samples each; reject disconnected and resample.")
    for N in (5, 6):
        rng = np.random.RandomState(42 + N)
        n_pass = 0
        worst_rel = 0.0
        for p in (0.3, 0.5, 0.7):
            for k in range(10):
                bonds = random_connected_erdos_renyi(N, p, rng)
                B, _ = bond_count_and_d2(N, bonds)
                predicted = cH * B * 4 ** (N - 2)
                H = xx_plus_yz_graph_h(N, bonds)
                observed = m_norm_squared(N, H, [0.1] * N)
                diff = abs(observed - predicted)
                rel = diff / predicted if predicted else diff
                worst_rel = max(worst_rel, rel)
                assert rel < 1e-9, (
                    f"N={N} p={p} sample {k}: rel = {rel:.3e} exceeds 1e-9; "
                    f"bonds = {bonds}, observed = {observed}, predicted = {predicted}"
                )
                n_pass += 1
        print(f"  N={N}: {n_pass} random graphs verified; worst |Δ|/predicted = {worst_rel:.3e}")


def section_4_disconnected_components(cH: float) -> None:
    """N=6 with two disjoint 3-site chains; the global B counts bonds across both
    components and the closed form remains ‖M‖² = c_H · B · 4^(N−2)."""
    print("\n(4) Disconnected components: N=6 = two disjoint 3-site chains.")
    print("    Global B = 2 + 2 = 4 (each chain has 2 bonds); closed form unchanged.")
    N = 6
    bonds = [(0, 1), (1, 2),    # left 3-chain on {0, 1, 2}
             (3, 4), (4, 5)]    # right 3-chain on {3, 4, 5}
    assert not is_connected(N, bonds), "expected the two-component graph to be disconnected"
    B, _ = bond_count_and_d2(N, bonds)
    predicted = cH * B * 4 ** (N - 2)
    H = xx_plus_yz_graph_h(N, bonds)
    observed = m_norm_squared(N, H, [0.1] * N)
    diff = abs(observed - predicted)
    rel = diff / predicted if predicted else diff
    print(f"  N=6 two-chain B={B}: observed = {observed:.6f}  predicted = {predicted:.6f}  "
          f"|Δ| = {diff:.3e}  (rel = {rel:.2e})")
    assert rel < 1e-9, (
        f"disconnected components: rel = {rel:.3e} exceeds 1e-9; "
        f"observed = {observed}, predicted = {predicted}"
    )


def section_5_weighted_edges(cH: float) -> None:
    """N=4 chain with per-bond weights J = [1.0, 2.0, 3.0]; predicted norm
    re-weights B → Σ_b J²_b at fixed c_H (anchored at J=1)."""
    print("\n(5) Weighted edges: N=4 chain with J = [1.0, 2.0, 3.0].")
    print("    Bilinear bond Hamiltonians contribute ‖H_b‖²_F ∝ J²_b; so the natural")
    print("    extension is B → Σ_b J²_b with c_H still anchored at J=1.")
    N = 4
    bonds = path_bonds(N)
    weights = [1.0, 2.0, 3.0]
    H = xx_plus_yz_graph_h(N, bonds, weights=weights)
    observed = m_norm_squared(N, H, [0.1] * N)
    B_weighted = sum(J * J for J in weights)
    predicted = cH * B_weighted * 4 ** (N - 2)
    diff = abs(observed - predicted)
    rel = diff / predicted if predicted else diff
    print(f"  N=4 weighted: Σ J²_b = {B_weighted:.3f}  "
          f"observed = {observed:.6f}  predicted = {predicted:.6f}  "
          f"|Δ| = {diff:.3e}  (rel = {rel:.2e})")
    assert rel < 1e-9, (
        f"weighted edges: rel = {rel:.3e} exceeds 1e-9; "
        f"observed = {observed}, predicted = {predicted}"
    )


def section_6_single_body_spot_check() -> None:
    """N=5 chain with bond-bilinear single-body H = Σ_b (I_i Y_j + Y_i I_j) on the
    chain bonds: F = (D2/2)·4^(N−2) for chain default D2 = 4N − 6."""
    print("\n(6) Single-body class spot check: H = Σ_b (I_i Y_j + Y_i I_j) on N=5 chain bonds.")
    print("    Anchor c_H_single from N=2 single-body bond; predicted F(N=5) = (D2/2)·4^(N−2)")
    print("    with chain D2 = 4N − 6 = 14.")
    cH_single = anchor_c_h_single_body()

    N = 5
    bonds = path_bonds(N)
    # Chain single-body default: D2 = 4N − 6 = 14, so F = (D2/2)·4^(N−2) = 7·64 = 448.
    _, D2 = bond_count_and_d2(N, bonds)
    # Chain bond list has D2 = Σ deg² = (degree-1 endpoints count 2 each, interior 2) = 4N − 6.
    assert D2 == 4 * N - 6, f"chain D2 mismatch: got {D2}, expected {4*N-6}"
    F_predicted = (D2 / 2) * 4 ** (N - 2)
    predicted = cH_single * F_predicted
    H = single_body_iy_plus_yi_graph_h(N, bonds)
    observed = m_norm_squared(N, H, [0.1] * N)
    diff = abs(observed - predicted)
    rel = diff / predicted if predicted else diff
    print(f"  N=5 chain single-body IY+YI: D2 = {D2}, F = (D2/2)·4^(N−2) = {F_predicted}")
    print(f"           observed = {observed:.6f}  predicted = {predicted:.6f}  "
          f"|Δ| = {diff:.3e}  (rel = {rel:.2e})")
    assert rel < 1e-9, (
        f"single-body N=5 chain: rel = {rel:.3e} exceeds 1e-9; "
        f"observed = {observed}, predicted = {predicted}"
    )


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #


def main() -> None:
    print("F1 general-topology verification")
    print("=" * 78)
    print("\n(1) Anchor c_H at N=2 single-bond XX+YZ (Π²-mixed, non-truly).")
    cH = anchor_c_h_main()
    section_2_named_graphs(cH)
    section_3_random_graphs(cH)
    section_4_disconnected_components(cH)
    section_5_weighted_edges(cH)
    section_6_single_body_spot_check()
    print("\nAll sections complete. The (B, D2) parameterisation of")
    print("  ‖M(N, G)‖²_F = c_H · F(N, G)")
    print("is bit-exact verified for connected and disconnected graphs, weighted")
    print("and unweighted edges, and the single-body class. F1 OpenQuestion closed.")


if __name__ == "__main__":
    main()
