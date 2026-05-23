"""F86e chain link 2: the full V_inter matrix in the F89 Bloch basis.

V_inter = Π_HD1 · M_H · Π_HD3  (on the c=2 operator block).
In the Bloch (OBC-sine operator) basis |k⟩⟨k₁k₂|:
  - M_H is DIAGONAL: M_H|k⟩⟨k₁k₂| = Δ·|...|, Δ = E_k − E_{k₁} − E_{k₂}
  - Π_HD1, Π_HD3 (computational-basis projectors) become Slater-overlap
    matrices Π̃_HD1, Π̃_HD3; these carry the structure.
  - V_inter_bloch = Π̃_HD1 · diag(Δ) · Π̃_HD3.

Question: are Π̃_HD1 / Π̃_HD3 clean, banded, sine-sum-structured matrices,
the kind whose product's operator norm F89's Chebyshev pipeline can close?

Operator-basis transform: an operator A (N × M_p2) → A_bloch = ψ · A · D,
so the vec-space super-operator transform is T = kron(ψ, Dᵀ).

Output: simulations/results/f86e_vinter_bloch/
"""
from __future__ import annotations

import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

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
    psi = np.zeros((N, N))
    for k in range(1, N + 1):
        for i in range(N):
            psi[k - 1, i] = math.sqrt(2.0 / (N + 1)) * math.sin(
                math.pi * k * (i + 1) / (N + 1))
    return psi


def sites_of(integer: int, N: int) -> list[int]:
    return [b for b in range(N) if (integer >> b) & 1]


def main() -> None:
    out_dir = Path("simulations/results/f86e_vinter_bloch")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 86)
    print("F86e: V_inter as a full matrix in the F89 Bloch basis")
    print("=" * 86)

    for N in (5, 6, 7):
        H_p1, st_p1, _ = build_h_block(N, 1)   # popcount-1, dim N
        H_p2, st_p2, _ = build_h_block(N, 2)   # popcount-2, dim C(N,2)
        Mp1, Mp2 = len(st_p1), len(st_p2)
        psi = obc_sine(N)

        # Slater matrix D[j_idx, pair_idx], NO 1/√2
        pairs = [(k1, k2) for k1 in range(1, N + 1) for k2 in range(k1 + 1, N + 1)]
        D = np.zeros((Mp2, len(pairs)))
        for jx, j in enumerate(st_p2):
            s1, s2 = sites_of(j, N)
            for px, (k1, k2) in enumerate(pairs):
                D[jx, px] = (psi[k1 - 1, s1] * psi[k2 - 1, s2]
                             - psi[k1 - 1, s2] * psi[k2 - 1, s1])

        # M_H = commutator [H,·] on the c=2 block (real symmetric)
        I_p1 = np.eye(Mp1); I_p2 = np.eye(Mp2)
        M_H = np.kron(H_p1, I_p2) - np.kron(I_p1, H_p2)

        # HD projectors (diagonal in computational op basis)
        block_dim = Mp1 * Mp2
        hd = np.zeros(block_dim, dtype=int)
        for ix, i in enumerate(st_p1):
            for jx, j in enumerate(st_p2):
                hd[ix * Mp2 + jx] = popcount(i ^ j)
        P_HD1 = np.diag((hd == 1).astype(float))
        P_HD3 = np.diag((hd == 3).astype(float))

        # vec-space Bloch transform T = kron(ψ, Dᵀ)
        T = np.kron(psi, D.T)

        M_H_bloch = T @ M_H @ T.T
        P_HD1_bloch = T @ P_HD1 @ T.T
        P_HD3_bloch = T @ P_HD3 @ T.T
        V_bloch = P_HD1_bloch @ M_H_bloch @ P_HD3_bloch

        # diagnostics
        offdiag = M_H_bloch - np.diag(np.diag(M_H_bloch))
        mh_offdiag_norm = float(np.linalg.norm(offdiag))
        # analytic Δ
        E = np.array([2 * math.cos(math.pi * k / (N + 1)) for k in range(1, N + 1)])
        Delta = np.array([E[k - 1] - E[k1 - 1] - E[k2 - 1]
                          for k in range(1, N + 1) for (k1, k2) in pairs])
        delta_match = float(np.linalg.norm(np.diag(M_H_bloch) - Delta))

        sigma_0 = float(np.linalg.svd(V_bloch, compute_uv=False)[0])
        # projector sparsity
        p1_nnz = int(np.sum(np.abs(P_HD1_bloch) > 1e-9))
        p3_nnz = int(np.sum(np.abs(P_HD3_bloch) > 1e-9))
        v_nnz = int(np.sum(np.abs(V_bloch) > 1e-9))

        print()
        print(f"--- N={N}  block_dim={block_dim} ---")
        print(f"  M_H_bloch off-diagonal norm: {mh_offdiag_norm:.2e}  "
              f"(0 → M_H diagonal in Bloch basis ✓)")
        print(f"  diag(M_H_bloch) vs analytic Δ=E_k−E_k1−E_k2: residual {delta_match:.2e}")
        print(f"  σ_0 from V_bloch = {sigma_0:.8f}")
        print(f"  Π̃_HD1 nonzeros: {p1_nnz}/{block_dim**2}  "
              f"({100*p1_nnz/block_dim**2:.1f}% dense)")
        print(f"  Π̃_HD3 nonzeros: {p3_nnz}/{block_dim**2}  "
              f"({100*p3_nnz/block_dim**2:.1f}% dense)")
        print(f"  V_bloch nonzeros: {v_nnz}/{block_dim**2}  "
              f"({100*v_nnz/block_dim**2:.1f}% dense)")
        # is Π̃_HD1 close to a projector?
        p1_proj_resid = float(np.linalg.norm(
            P_HD1_bloch @ P_HD1_bloch - P_HD1_bloch))
        print(f"  Π̃_HD1 idempotency ‖Π²−Π‖ = {p1_proj_resid:.2e} "
              f"(projector preserved under unitary transform)")

        # heatmaps
        fig, axes = plt.subplots(1, 4, figsize=(22, 5.2))
        for ax, M, title in [
            (axes[0], M_H_bloch, "M_H (Bloch): diagonal Δ"),
            (axes[1], P_HD1_bloch, "Π̃_HD1 (Bloch)"),
            (axes[2], P_HD3_bloch, "Π̃_HD3 (Bloch)"),
            (axes[3], V_bloch, "V_inter (Bloch)"),
        ]:
            im = ax.imshow(np.abs(M), cmap="magma", aspect="auto")
            ax.set_title(f"{title}\n|entry| magnitude")
            plt.colorbar(im, ax=ax, shrink=0.7)
        fig.suptitle(f"V_inter in the F89 Bloch basis, N={N} "
                     f"(block_dim {block_dim})", fontsize=12)
        plt.tight_layout()
        f = out_dir / f"vinter_bloch_N{N}.png"
        plt.savefig(f, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"  heatmap saved: {f}")

    print()
    print("Reading: M_H diagonal in Bloch basis (expected). The structure of")
    print("σ_0 lives in Π̃_HD1 and Π̃_HD3; if they are banded / low-rank /")
    print("sine-sum-patterned, σ_0 = ‖Π̃_HD1 · diagΔ · Π̃_HD3‖ is Chebyshev-closable.")


if __name__ == "__main__":
    main()
