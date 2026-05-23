"""F86e chain link: V_inter in the F89 (SE,DE) Bloch / OBC-sine basis.

Careful step (Tom 2026-05-20): build V_inter in F89's Bloch basis and see
whether it becomes a clean sine-sum matrix, the structure F89's Chebyshev
machinery handles.

σ_0 is uniform-J (uses M_H_total, not per-bond) → it sits in F89's natural
orbit-reducible category, NOT in F86b₂'s per-bond obstruction territory.

Transform:
  SE side (popcount-1):  |i⟩ = Σ_k ψ_k(i)|k⟩,  ψ = N×N OBC sine matrix.
  DE side (popcount-2):  |j={s1,s2}⟩ = Σ_{k1<k2} D[j,(k1,k2)] |k1 k2⟩
     with D[j,(k1,k2)] = ψ_{k1}(s1)ψ_{k2}(s2) − ψ_{k1}(s2)ψ_{k2}(s1).
     NO 1/√2: both bases are already 2-particle Slater states; the overlap
     is the full 2×2 determinant. (This was the Phase-B factor-2 bug.)

For each N: transform U_top, V_top (the singular vectors carrying σ_0) into
the Bloch operator basis |k⟩⟨k1 k2|, verify norm preservation (= transform
is unitary = bug fixed), and read the structure.
"""
from __future__ import annotations

import math
import sys
from itertools import combinations

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def popcount(x: int) -> int:
    return bin(x).count("1")


def build_h_block(N: int, p: int):
    states = sorted(sum(1 << b for b in c) for c in combinations(range(N), p))
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s in states:
        for b in range(N - 1):
            if ((s >> b) & 1) != ((s >> (b + 1)) & 1):
                s2 = s ^ ((1 << b) | (1 << (b + 1)))
                H[idx[s2], idx[s]] += 1.0
    return H, states, idx


def obc_sine(N: int) -> np.ndarray:
    """psi[k-1, site] = √(2/(N+1)) sin(πk(site+1)/(N+1))."""
    psi = np.zeros((N, N))
    for k in range(1, N + 1):
        for i in range(N):
            psi[k - 1, i] = math.sqrt(2.0 / (N + 1)) * math.sin(
                math.pi * k * (i + 1) / (N + 1))
    return psi


def sites_of(integer: int, N: int) -> list[int]:
    return [b for b in range(N) if (integer >> b) & 1]


def build_v_inter(N: int):
    """V_inter (HD-1 rows, HD-3 cols) in computational basis + labels."""
    H_p1, st_p1, ix_p1 = build_h_block(N, 1)
    H_p2, st_p2, ix_p2 = build_h_block(N, 2)
    hd1 = [(ip, jp) for ip in st_p1 for jp in st_p2 if (ip & jp) == ip]
    hd3 = [(ii, jj) for ii in st_p1 for jj in st_p2 if (ii & jj) == 0]
    i1 = np.array([p[0] for p in hd1]); j1 = np.array([p[1] for p in hd1])
    i3 = np.array([p[0] for p in hd3]); j3 = np.array([p[1] for p in hd3])
    i1x = np.array([ix_p1[x] for x in i1]); i3x = np.array([ix_p1[x] for x in i3])
    j1x = np.array([ix_p2[x] for x in j1]); j3x = np.array([ix_p2[x] for x in j3])
    term1 = H_p1[i1x[:, None], i3x[None, :]] * (j1[:, None] == j3[None, :])
    term2 = -H_p2[j1x[:, None], j3x[None, :]] * (i1[:, None] == i3[None, :])
    V = term1 + term2
    return V, hd1, hd3, st_p1, st_p2


def slater_matrix(N: int, st_p2: list[int], psi: np.ndarray) -> np.ndarray:
    """D[j_idx, pair_idx] = 2x2 Slater determinant, NO 1/√2."""
    pairs = [(k1, k2) for k1 in range(1, N + 1) for k2 in range(k1 + 1, N + 1)]
    D = np.zeros((len(st_p2), len(pairs)))
    for jx, j in enumerate(st_p2):
        s1, s2 = sites_of(j, N)
        for px, (k1, k2) in enumerate(pairs):
            D[jx, px] = (psi[k1 - 1, s1] * psi[k2 - 1, s2]
                         - psi[k1 - 1, s2] * psi[k2 - 1, s1])
    return D, pairs


def to_bloch_operator(vec, labels, N, st_p1, st_p2, psi, D):
    """Transform an HD-subspace operator-vector into the Bloch operator basis
    |k⟩⟨k1 k2|. Returns array shape (N, n_pairs)."""
    ix_p1 = {s: i for i, s in enumerate(st_p1)}
    ix_p2 = {s: i for i, s in enumerate(st_p2)}
    M_p2 = len(st_p2)
    # operator as N × M_p2 matrix
    op = np.zeros((N, M_p2), dtype=complex)
    for (i_int, j_int), amp in zip(labels, vec):
        # i_int is popcount-1 → site index; map to row by SITE
        site = sites_of(i_int, N)[0]
        op[site, ix_p2[j_int]] = amp
    # α[k, pair] = Σ_site Σ_j ψ[k,site]·op[site,j]·D[j,pair]
    alpha = psi @ op @ D
    return alpha


def main() -> None:
    print("=" * 88)
    print("F86e: V_inter in the F89 (SE,DE) Bloch basis, careful chain link")
    print("=" * 88)

    for N in (5, 6, 7):
        V, hd1, hd3, st_p1, st_p2 = build_v_inter(N)
        psi = obc_sine(N)
        D, pairs = slater_matrix(N, st_p2, psi)

        # unitarity check on D
        DDt = D @ D.T
        d_unit = float(np.linalg.norm(DDt - np.eye(len(pairs))))

        U, S, Vh = np.linalg.svd(V, full_matrices=False)
        sigma_0 = float(S[0])
        U_top = U[:, 0]
        V_top = Vh[0, :].conj()

        # transform U_top (HD-1) and V_top (HD-3) to Bloch operator basis
        U_bloch = to_bloch_operator(U_top, hd1, N, st_p1, st_p2, psi, D)
        V_bloch = to_bloch_operator(V_top, hd3, N, st_p1, st_p2, psi, D)

        print()
        print(f"--- N={N}  σ_0={sigma_0:.8f}  (D unitarity residual {d_unit:.1e}) ---")
        print(f"  |U_top|² computational = {np.sum(np.abs(U_top)**2):.6f}")
        print(f"  |U_top|² in Bloch basis = {np.sum(np.abs(U_bloch)**2):.6f}  "
              f"(must equal 1.0 → transform unitary, Phase-B bug fixed)")
        print(f"  |V_top|² in Bloch basis = {np.sum(np.abs(V_bloch)**2):.6f}")

        # single-particle marginal of U_bloch (sum over Slater pair index)
        u_k = np.sum(np.abs(U_bloch) ** 2, axis=1)
        v_k = np.sum(np.abs(V_bloch) ** 2, axis=1)
        print(f"  U_top single-particle marginal |c_k|² (k=1..N):")
        print("    " + "  ".join(f"k{k+1}:{u_k[k]:.4f}" for k in range(N)))
        print(f"  V_top single-particle marginal |c_k|² (k=1..N):")
        print("    " + "  ".join(f"k{k+1}:{v_k[k]:.4f}" for k in range(N)))

        # Slater-pair marginal of U_bloch
        u_pair = np.sum(np.abs(U_bloch) ** 2, axis=0)
        top_pairs = np.argsort(u_pair)[::-1][:5]
        print(f"  U_top top-5 Slater pairs:")
        for px in top_pairs:
            print(f"    {pairs[px]}: |c|²={u_pair[px]:.4f}")


if __name__ == "__main__":
    main()
