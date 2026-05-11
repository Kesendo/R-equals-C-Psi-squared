"""F89 path-6: closed-form for the F_a (AT-locked rate 2γ) mode amplitudes.

Path-6 (N_block=7) has SE Bloch eigenvalues E_n = 4J·cos(πn/8) for n=1..7.
S_2-anti modes are n=2, 4, 6 with:
  E_2 = 4J·cos(π/4) = 2√2·J
  E_4 = 4J·cos(π/2) = 0
  E_6 = 4J·cos(3π/4) = -2√2·J

So three F_a modes at ω/J ∈ {+2√2, 0, -2√2} — algebraic extension Q[√2].

Expected: closed-form amplitudes in Q[√2] (analog of path-3's Q[√5]).
"""

from __future__ import annotations

import sys
from itertools import combinations

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

N_BLOCK = 7


def build_se_de(J: float, gamma: float):
    de_pairs = list(combinations(range(N_BLOCK), 2))
    basis = [(i, jk) for i in range(N_BLOCK) for jk in de_pairs]
    n_basis = len(basis)

    M_SE = np.zeros((N_BLOCK, N_BLOCK))
    for a in range(N_BLOCK):
        for b in range(N_BLOCK):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J

    n_de = len(de_pairs)
    M_DE = np.zeros((n_de, n_de))
    for idx, (j, k) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j < N_BLOCK and new_j != k:
                new_pair = tuple(sorted([new_j, k]))
                if abs(new_j - j) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2 * J
        for new_k in [k - 1, k + 1]:
            if 0 <= new_k < N_BLOCK and new_k != j:
                new_pair = tuple(sorted([j, new_k]))
                if abs(new_k - k) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2 * J

    L = np.zeros((n_basis, n_basis), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(N_BLOCK):
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

    return L, basis, de_pairs


def per_site_reduction(basis, n_block=N_BLOCK):
    w = np.zeros((n_block, len(basis)), dtype=float)
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue
            other = jk[1] if jk[0] == l else jk[0]
            if other == i:
                w[l, idx] = 1.0
    return w


def main() -> None:
    Ns = [8, 9, 10, 11, 12, 14, 16, 20]  # N ≥ N_block + 1 = 8
    print(f"# F89 path-6: AT-locked amplitude closed-form derivation")
    print(f"# N_block = {N_BLOCK} (path-6); Ns = {Ns}\n")

    gamma = 1.0
    q = 1.5
    J = q * gamma
    L, basis, _ = build_se_de(J, gamma)
    eigvals, eigvecs = np.linalg.eig(L)
    rates = -eigvals.real / gamma
    w = per_site_reduction(basis)

    at_idx = np.where((np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6))[0]
    print(f"## {len(at_idx)} AT-locked modes at q=1.5\n")

    sigs_data = {}
    for k in at_idx:
        v = eigvecs[:, k]
        v = v / np.linalg.norm(v)
        Mv = w @ v
        sigs_per_N = []
        for N in Ns:
            pre = np.sqrt(2 / (N**2 * (N - 1)))
            rho = np.full(len(basis), pre / 2, dtype=complex)
            c = np.vdot(v, rho)
            sigs = abs(c) ** 2 * np.sum(np.abs(Mv) ** 2)
            scaled = sigs * (N * N * (N - 1))
            sigs_per_N.append(scaled)
        sigs_data[k] = (eigvals[k], sigs_per_N)

    nontrivial = [(k, lam, sigs) for k, (lam, sigs) in sigs_data.items()
                  if abs(np.mean(sigs)) > 1e-15 and abs(rates[k] - 2.0) < 1e-6]

    print(f"# {len(nontrivial)} F_a modes with non-trivial amplitude:\n")
    print("| Eigenvalue | freq ω/J | sigs · N²(N-1) |")
    print("|---|---|---|")
    for k, lam, sigs in nontrivial:
        const = all(abs(sigs[i] - sigs[0]) < 1e-9 for i in range(len(sigs)))
        val = np.mean(sigs)
        print(f"| {lam:+.6f} | {lam.imag/J:+.6f} | {val:.10f} ({'const' if const else 'varies'}) |")

    fa_sum = sum(np.mean(sigs) for k, lam, sigs in nontrivial)
    print(f"\n# Σ F_a sigs · N²(N-1) = {fa_sum:.10f}")

    print("\n## Recognising in algebraic extensions:")
    extensions = [
        (None, "rational"),
        ([sp.sqrt(2)], "Q[√2]"),
        ([sp.sqrt(2), sp.sqrt(3)], "Q[√2, √3]"),
        ([sp.sqrt(2), sp.sqrt(5)], "Q[√2, √5]"),
    ]
    for ext, label in extensions:
        try:
            rat = sp.nsimplify(fa_sum, ext if ext else [], rational=False if ext else True, tolerance=1e-9)
            if not rat.has(sp.Float):
                print(f"#   Σ F_a in {label}: {rat} ≈ {float(rat):.10f}")
                if ext is not None and label != "rational":
                    break
        except Exception:
            pass

    print("\n## Per-mode symbolic recognition:")
    for k, lam, sigs in nontrivial:
        val = np.mean(sigs)
        print(f"# λ={lam:+.4f}, ω/J={lam.imag/J:+.4f}: sigs · N²(N-1) = {val:.10f}")
        for ext, label in extensions:
            try:
                rat = sp.nsimplify(val, ext if ext else [], rational=False if ext else True, tolerance=1e-9)
                if not rat.has(sp.Float):
                    print(f"#   {label}: {rat}")
                    break
            except Exception:
                pass


if __name__ == "__main__":
    main()
