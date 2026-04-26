#!/usr/bin/env python3
"""Generate the super-operator block-norm heatmap as a real image file.

Produces `visualizations/super_operator_blocks.png`: three panels showing
the (w_row, w_col) Frobenius-norm matrix of M = Π·L·Π⁻¹ + L + 2Σγ·I in
Pauli-string basis for the three Hamiltonian categories at N=3, γ=0.1:

  - truly_unbroken  (XX + YY)   →  M ≡ 0 at machine precision
  - soft_broken     (XY + YX)   →  ‖M‖ = 32, structured break
  - hard_broken     (XX + XY)   →  ‖M‖ = 22.6, eigenvalue spectrum also broken

Same matrices that `_derive_from_below.py` Stage 8b prints as ASCII.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def pauli_string_xy_weight(k, N):
    weight = 0
    kk = k
    for _ in range(N):
        if kk % 4 & 1:
            weight += 1
        kk //= 4
    return weight


def block_norms_for(label, terms, N=3, gamma=0.1):
    bonds = [(i, i + 1) for i in range(N - 1)]
    bilinear_terms = [(t[0], t[1], 1.0) for t in terms]
    H = fw._build_bilinear(N, bonds, bilinear_terms)
    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    M = fw.palindrome_residual(L, N * gamma, N)
    weights = np.array([pauli_string_xy_weight(k, N) for k in range(4 ** N)])
    blocks = np.zeros((N + 1, N + 1))
    for w_row in range(N + 1):
        rows = np.where(weights == w_row)[0]
        for w_col in range(N + 1):
            cols = np.where(weights == w_col)[0]
            blocks[w_row, w_col] = float(np.linalg.norm(M[np.ix_(rows, cols)]))
    return blocks, float(np.linalg.norm(M))


def main():
    N = 3
    gamma = 0.1
    cases = [
        ('truly_unbroken\nH = J(XX + YY)', [('X', 'X'), ('Y', 'Y')]),
        ('soft_broken\nH = J(XY + YX)', [('X', 'Y'), ('Y', 'X')]),
        ('hard_broken\nH = J(XX + XY)', [('X', 'X'), ('X', 'Y')]),
    ]
    block_data = []
    for label, terms in cases:
        blocks, norm = block_norms_for(label, terms, N, gamma)
        block_data.append((label, blocks, norm))

    # Common color scale across the three panels for fair comparison.
    # Use log scale so machine-zero ('truly') still renders distinctly.
    finite_max = max(b.max() for _, b, _ in block_data if b.max() > 0)
    finite_min_nonzero = min(
        b[b > 0].min() if (b > 0).any() else finite_max
        for _, b, _ in block_data
    )
    vmin = max(finite_min_nonzero, 1e-15)
    vmax = finite_max

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    for ax, (label, blocks, total_norm) in zip(axes, block_data):
        # Replace zeros with a very small value so log norm doesn't fail
        plot_blocks = np.where(blocks > 0, blocks, vmin / 10)
        im = ax.imshow(
            plot_blocks,
            cmap='viridis',
            norm=LogNorm(vmin=vmin, vmax=vmax),
            origin='upper',
            interpolation='nearest',
        )
        ax.set_title(f"{label}\n‖M‖ = {total_norm:.2e}", fontsize=10)
        ax.set_xticks(range(N + 1))
        ax.set_yticks(range(N + 1))
        ax.set_xlabel("XY-weight w_col")
        ax.set_ylabel("XY-weight w_row")
        for i in range(N + 1):
            for j in range(N + 1):
                v = blocks[i, j]
                txt = "0" if v < 1e-12 else f"{v:.1f}"
                color = "white" if v > vmax * 0.3 or v < 1e-12 else "black"
                ax.text(j, i, txt, ha='center', va='center', fontsize=8, color=color)

    fig.colorbar(im, ax=axes, label="‖M_(w_row, w_col)‖ (Pauli-block Frobenius norm)",
                 shrink=0.8, location='right')
    fig.suptitle(
        f"Super-operator residual M = Π·L·Π⁻¹ + L + 2Σγ·I  by Pauli weight blocks  (N={N}, γ={gamma})",
        fontsize=12,
    )

    out_path = SCRIPT_DIR.parent / "visualizations" / "super_operator_blocks.png"
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
