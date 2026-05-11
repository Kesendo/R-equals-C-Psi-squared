"""F89 path-4, path-5: AT-lock pattern + closed-form-frequency identification.

Generalisation of path-3 results to larger blocks. For each path-k (block of
k+1 sites with k bonds + Z-deph):

  1. Build the (SE, DE) sub-block of L_super (dim = (k+1)·C(k+1,2))
  2. Project to S_2-sym subspace (mirror j ↔ k-j)
  3. Diagonalise; identify AT-locked modes at rate exactly 2γ (overlap)
     and rate exactly 6γ (no-overlap)
  4. Verify their frequencies match SE-anti single-particle Bloch eigenvalues
     E_n = 4J·cos(πn/(N_block+1)) for n in S_2-anti orbit (n even when
     N_block-1 is odd, etc.)
  5. Verify eigvec lives entirely in overlap-only / no-overlap-only subspace
  6. Count: total = 2·|S_2-anti SE Bloch orbit| (matches path-3's 4 = 2·2)

Expectations:
  - Path-4 (N_block=5): SE-anti Bloch n=2,4 with E=±2J (rational, no √).
    F_a/F_b each give 2 roots → 4 AT-locked total.
  - Path-5 (N_block=6): SE-anti Bloch n=2,4 (5 was self-mirror at k=3 for
    odd N_block; for N_block=6 even, anti modes are n=2, 4, 6, all distinct).
    cos(πn/7) for n=2,4,6 are Cardano-cubic roots.
    F_a/F_b each give 3 roots → 6 AT-locked total.
"""

from __future__ import annotations

import sys
from itertools import combinations

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_pathk_se_de_sym(J: float, gamma: float, n_block: int):
    """Build (SE, DE) S_2-sym sub-block of L_super for path-k where k = n_block-1.

    Returns (L_sym, basis_24/etc, P_projector).
    """
    de_pairs = list(combinations(range(n_block), 2))
    basis = [(i, jk) for i in range(n_block) for jk in de_pairs]
    n_basis = len(basis)

    # H_B^SE: tridiagonal n_block × n_block, off-diag 2J
    M_SE = np.zeros((n_block, n_block))
    for a in range(n_block):
        for b in range(n_block):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J

    # H_B^DE: hopping in 2-particle space
    n_de = len(de_pairs)
    M_DE = np.zeros((n_de, n_de))
    for idx, (j, k) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j < n_block and new_j != k:
                new_pair = tuple(sorted([new_j, k]))
                if abs(new_j - j) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2 * J
        for new_k in [k - 1, k + 1]:
            if 0 <= new_k < n_block and new_k != j:
                new_pair = tuple(sorted([j, new_k]))
                if abs(new_k - k) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2 * J

    L = np.zeros((n_basis, n_basis), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(n_block):
            if M_SE[i2, i] != 0:
                idx2 = basis.index((i2, jk))
                L[idx2, idx] += -1j * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(n_de):
            if M_DE[jk_idx, jk2_idx] != 0:
                jk2 = de_pairs[jk2_idx]
                idx2 = basis.index((i, jk2))
                L[idx2, idx] += 1j * M_DE[jk_idx, jk2_idx]
        L[idx, idx] += -2 * gamma if i in jk else -6 * gamma

    # S_2 mirror: site l ↔ n_block - 1 - l
    perm_se = lambda i: n_block - 1 - i
    perm_de = lambda jk: tuple(sorted([perm_se(jk[0]), perm_se(jk[1])]))

    cols_sym = []
    cols_anti = []
    handled = set()
    for idx, (i, jk) in enumerate(basis):
        if idx in handled:
            continue
        idx2 = basis.index((perm_se(i), perm_de(jk)))
        v_sym = np.zeros(n_basis, dtype=complex)
        if idx == idx2:
            v_sym[idx] = 1.0
            cols_sym.append(v_sym)
        else:
            v_sym[idx] = 1.0 / np.sqrt(2)
            v_sym[idx2] = 1.0 / np.sqrt(2)
            cols_sym.append(v_sym)
            v_anti = np.zeros(n_basis, dtype=complex)
            v_anti[idx] = 1.0 / np.sqrt(2)
            v_anti[idx2] = -1.0 / np.sqrt(2)
            cols_anti.append(v_anti)
            handled.add(idx2)
        handled.add(idx)
    P_sym = np.column_stack(cols_sym)
    L_sym = P_sym.conj().T @ L @ P_sym

    return L_sym, basis, P_sym, de_pairs


def analyse_pathk(k: int, J: float, gamma: float):
    """Analyse path-k (n_block = k+1) AT-lock structure."""
    n_block = k + 1
    L_sym, basis, P_sym, de_pairs = build_pathk_se_de_sym(J, gamma, n_block)
    n_sub = L_sym.shape[0]
    print(f"\n## Path-{k} (N_block={n_block}): (SE, DE) sub-block dim {len(basis)}, S_2-sym dim {n_sub}")

    # Predicted SE single-particle Bloch eigenvalues
    E_se = [4 * J * np.cos(np.pi * n / (n_block + 1)) for n in range(1, n_block + 1)]
    print(f"# SE single-particle Bloch eigenvalues E_n = 4J·cos(πn/{n_block+1}):")
    for n, e in enumerate(E_se, 1):
        sym_class = "sym" if (n % 2 == 1) else "anti"  # ψ_n S_2-sym for n odd at N_block-1 even
        # Actually: ψ_k(N_block-1-j) = (-1)^(k+1) ψ_k(j), so sym for k odd (when k+1 even).
        # Mirror sites 0 ↔ n_block-1: for ψ_k(j) = sin(πk(j+1)/(n_block+1)), mirror gives sin(πk(n_block-j)/(n_block+1))
        # = sin(πk - πk(j+1)/(n_block+1)) = (-1)^(k+1) sin(πk(j+1)/(n_block+1)).
        # So sym for k odd (sign = +1 when (-1)^(k+1) = +1, i.e., k+1 even, k odd).
        print(f"#   E_{n} = {e:+.4f}  (S_2-{sym_class})")

    # Diagonalise L_sym, find AT-locked eigenvalues
    eigvals, eigvecs_sub = np.linalg.eig(L_sym)
    rates = -eigvals.real / gamma
    freqs = eigvals.imag / J

    at_2g_idx = np.where(np.abs(rates - 2.0) < 1e-6)[0]
    at_6g_idx = np.where(np.abs(rates - 6.0) < 1e-6)[0]
    n_at = len(at_2g_idx) + len(at_6g_idx)
    print(f"\n# AT-locked count: rate-2γ: {len(at_2g_idx)}, rate-6γ: {len(at_6g_idx)}, total: {n_at}")

    # SE-anti Bloch eigenvalues (predicted F_a/F_b frequencies)
    se_anti_E = [E_se[n - 1] for n in range(1, n_block + 1) if n % 2 == 0]
    print(f"# SE-anti Bloch eigenvalues (predicted F_a/F_b ω): {[f'{e:+.4f}' for e in se_anti_E]}")
    print(f"#   Predicted AT-locked count: 2 (F_a + F_b) × {len(se_anti_E)} = {2*len(se_anti_E)}")

    # Match each AT-locked eigenvalue to a SE-anti Bloch eigenvalue
    print("\n# Verifying AT-locked frequencies match SE-anti Bloch eigenvalues:")
    print("| Eigenvalue λ | rate Γ/γ | freq ω/J | Closest SE-anti E_n? |")
    print("|---|---|---|---|")
    for k_idx in sorted(np.concatenate([at_2g_idx, at_6g_idx])):
        ω = eigvals[k_idx].imag
        # Find closest SE-anti Bloch eigenvalue
        diffs = [abs(ω - e) for e in se_anti_E]
        closest = se_anti_E[np.argmin(diffs)]
        print(f"| {eigvals[k_idx]:.4f} | {rates[k_idx]:.4f} | {freqs[k_idx]:+.4f} | E_? = {closest:+.4f}, diff = {abs(ω - closest):.2e} |")

    # Verify eigvec overlap-only / no-overlap-only support
    overlap_mask = np.array([i in jk for i, jk in basis])
    eigvecs_full = P_sym @ eigvecs_sub  # lift to full 2D basis

    print("\n# Eigvec overlap support for AT-locked modes:")
    print("| Eigenvalue | overlap fraction | no-overlap fraction |")
    print("|---|---|---|")
    for k_idx in sorted(np.concatenate([at_2g_idx, at_6g_idx])):
        v = eigvecs_full[:, k_idx]
        v = v / np.linalg.norm(v)
        ov = np.sum(np.abs(v[overlap_mask]) ** 2)
        no = np.sum(np.abs(v[~overlap_mask]) ** 2)
        print(f"| {eigvals[k_idx]:.4f} | {ov:.4f} | {no:.4f} |")

    print(f"\n# H_B-mixed (non-AT-locked) sub-factor: degree {n_sub - n_at} (likely irreducible / Galois-non-solvable for ≥5)")


def main() -> None:
    J, gamma = 1.5, 1.0  # q = 1.5
    print(f"# F89 path-4 / path-5 AT-lock structure scan (J={J}, γ={gamma}, q={J/gamma})\n")

    for k in [3, 4, 5, 6]:  # path-3 (sanity check), path-4, path-5, path-6
        analyse_pathk(k, J, gamma)


if __name__ == "__main__":
    main()
