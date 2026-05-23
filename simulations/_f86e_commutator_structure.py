"""F86e chain link 3: structure of the commutator [Π_HD1, M_H].

σ_0 = ‖[Π_HD1, M_H]‖ (verified, typed in master as SigmaZeroCommutatorNormClaim).

Next link: in the Bloch basis M_H = diag(Δ), so
  [Π_HD1, M_H]_{ab} = Π̃_HD1[a,b]·Δ_b − Δ_a·Π̃_HD1[a,b] = Π̃_HD1[a,b]·(Δ_b − Δ_a).
The commutator is the Hadamard product  Π̃_HD1 ⊙ ΔDiff,  ΔDiff[a,b] = Δ_b − Δ_a.

This script:
  1. builds [Π_HD1, M_H] in the Bloch basis,
  2. verifies the Hadamard identity bit-exact,
  3. orders the basis by Δ and examines the matrix structure:
     Toeplitz (constant diagonals) / Hankel (constant anti-diagonals) /
     banded, the structures whose operator-norm asymptotics are closed-form.

Output: simulations/results/f86e_commutator/
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
    return H, states


def obc_sine(N: int) -> np.ndarray:
    psi = np.zeros((N, N))
    for k in range(1, N + 1):
        for i in range(N):
            psi[k - 1, i] = math.sqrt(2.0 / (N + 1)) * math.sin(
                math.pi * k * (i + 1) / (N + 1))
    return psi


def sites_of(integer: int, N: int) -> list[int]:
    return [b for b in range(N) if (integer >> b) & 1]


def diagonal_band_uniformity(M: np.ndarray) -> float:
    """Mean coefficient of variation along diagonals (Toeplitz → 0)."""
    n = M.shape[0]
    covs = []
    for d in range(-(n - 1), n):
        diag = np.diag(M, d)
        big = np.abs(diag) > 1e-6
        if np.sum(big) >= 3:
            vals = np.abs(diag[big])
            covs.append(np.std(vals) / (np.mean(vals) + 1e-15))
    return float(np.mean(covs)) if covs else float("nan")


def antidiagonal_band_uniformity(M: np.ndarray) -> float:
    """Mean CV along anti-diagonals (Hankel → 0)."""
    return diagonal_band_uniformity(np.fliplr(M))


def main() -> None:
    out_dir = Path("simulations/results/f86e_commutator")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 84)
    print("F86e: structure of the commutator [Π_HD1, M_H]")
    print("=" * 84)

    for N in (5, 6, 7):
        H_p1, st_p1 = build_h_block(N, 1)
        H_p2, st_p2 = build_h_block(N, 2)
        Mp1, Mp2 = len(st_p1), len(st_p2)
        psi = obc_sine(N)

        pairs = [(k1, k2) for k1 in range(1, N + 1) for k2 in range(k1 + 1, N + 1)]
        D = np.zeros((Mp2, len(pairs)))
        for jx, j in enumerate(st_p2):
            s1, s2 = sites_of(j, N)
            for px, (k1, k2) in enumerate(pairs):
                D[jx, px] = (psi[k1 - 1, s1] * psi[k2 - 1, s2]
                             - psi[k1 - 1, s2] * psi[k2 - 1, s1])

        M_H = np.kron(H_p1, np.eye(Mp2)) - np.kron(np.eye(Mp1), H_p2)
        block_dim = Mp1 * Mp2
        hd = np.array([popcount(i ^ j) for i in st_p1 for j in st_p2])
        P_HD1 = np.diag((hd == 1).astype(float))

        # commutator in computational basis
        comm_comp = P_HD1 @ M_H - M_H @ P_HD1
        sigma_0 = float(np.linalg.svd(comm_comp, compute_uv=False)[0])

        # transform to Bloch basis
        T = np.kron(psi, D.T)
        comm_bloch = T @ comm_comp @ T.T
        P_HD1_bloch = T @ P_HD1 @ T.T

        # Δ vector (diagonal of M_H in Bloch basis)
        E = np.array([2 * math.cos(math.pi * k / (N + 1)) for k in range(1, N + 1)])
        Delta = np.array([E[k - 1] - E[k1 - 1] - E[k2 - 1]
                          for k in range(1, N + 1) for (k1, k2) in pairs])

        # Hadamard identity: comm_bloch[a,b] should equal P_HD1_bloch[a,b]·(Δ_b − Δ_a)
        delta_diff = Delta[None, :] - Delta[:, None]
        hadamard = P_HD1_bloch * delta_diff
        hadamard_resid = float(np.linalg.norm(comm_bloch - hadamard))

        # order basis by Δ
        order = np.argsort(Delta)
        comm_ord = comm_bloch[np.ix_(order, order)]

        toeplitz_cv = diagonal_band_uniformity(comm_ord)
        hankel_cv = antidiagonal_band_uniformity(comm_ord)

        print()
        print(f"--- N={N}  block_dim={block_dim} ---")
        print(f"  σ_0 = ‖[Π_HD1, M_H]‖ = {sigma_0:.8f}")
        print(f"  Hadamard identity  [Π,M]_bloch = Π̃_HD1 ⊙ ΔDiff:  "
              f"residual {hadamard_resid:.2e}  "
              f"{'✓ exact' if hadamard_resid < 1e-9 else '✗'}")
        print(f"  Δ-ordered commutator structure:")
        print(f"    Toeplitz CV (diagonals, →0 if Toeplitz):     {toeplitz_cv:.4f}")
        print(f"    Hankel   CV (anti-diagonals, →0 if Hankel):  {hankel_cv:.4f}")

        # heatmaps
        fig, axes = plt.subplots(1, 3, figsize=(17, 5))
        for ax, M, ttl in [
            (axes[0], comm_bloch, "[Π_HD1, M_H] (Bloch basis)"),
            (axes[1], comm_ord, "[Π_HD1, M_H] (Δ-ordered)"),
            (axes[2], delta_diff[np.ix_(order, order)], "ΔDiff = Δ_b−Δ_a (Δ-ordered)"),
        ]:
            im = ax.imshow(np.abs(M) if "Diff" not in ttl else M,
                           cmap="magma" if "Diff" not in ttl else "RdBu",
                           aspect="auto")
            ax.set_title(ttl)
            plt.colorbar(im, ax=ax, shrink=0.7)
        fig.suptitle(f"Commutator [Π_HD1, M_H] structure, N={N}", fontsize=12)
        plt.tight_layout()
        f = out_dir / f"commutator_N{N}.png"
        plt.savefig(f, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"  heatmap: {f}")

    print()
    print("Reading: Hadamard identity confirms [Π,M] = Π̃_HD1 ⊙ ΔDiff.")
    print("Δ-ordered: if Toeplitz CV → 0, the commutator is a convolution operator")
    print("and σ_0(∞) = sup of its Fourier symbol (closed form). If Hankel CV → 0,")
    print("Nehari gives σ_0(∞) as a symbol distance. If neither, the structure is")
    print("a genuine Schur-multiplier norm.")


if __name__ == "__main__":
    main()
