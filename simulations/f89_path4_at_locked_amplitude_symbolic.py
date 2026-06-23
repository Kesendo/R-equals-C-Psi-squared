"""F89 path-4: closed-form for the F_a (AT-locked rate 2γ) mode amplitudes.

Path-4 (N_block=5): F_a has 2 AT-locked modes at λ = -2γ ± i·3J (rational
half-integer Bloch eigenvalues E_2 = 2J, E_4 = -2J → ω = ±2J... wait
actually checking: SE-anti Bloch at N_block=5 are n=2, 4 with
E_n = 4J·cos(πn/6). E_2 = 4J·cos(π/3) = 2J, E_4 = 4J·cos(2π/3) = -2J.

But the AT-lock scan reported F_a eigenvalues at -2γ ± 3iJ (not ±2iJ).
Recompute: at q=1.5, J=1.5, the F_a freq was +3 (in absolute units).
That's 3 = 2J (since J=1.5, 2J=3.0). ✓ So F_a freq = 2J = E_2.

Strategy (analog of path-3 amplitude script):
  1. Numerical eigvec extraction at multiple q (verify q-independence).
  2. Inner product with ρ_block(0)|_(SE,DE) (uniform vector pre/2 on 50 entries).
  3. Per-site reduction sigs for F_a modes.
  4. Fit sigs(N) · N²(N−1) for closed form (expected: rational with maybe √3
     if F_b multi-particle contributes).
"""

from __future__ import annotations

import sys
from itertools import combinations

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

N_BLOCK = 5


def build_se_de(J: float, gamma: float):
    """Build (SE, DE) sub-block of L_super for path-4 (N_block=5)."""
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
    """Per-site reduction matrix for (SE, DE) basis pairs."""
    w = np.zeros((n_block, len(basis)), dtype=float)
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue
            other = jk[1] if jk[0] == l else jk[0]
            if other == i:
                w[l, idx] = 1.0
    return w


def find_at_locked(L, gamma, J):
    """Find F_a, F_b eigenvalues + eigenvectors. Returns dict label → (lam, v)."""
    eigvals, eigvecs = np.linalg.eig(L)
    rates = -eigvals.real / gamma

    matched = {}
    for k in range(len(eigvals)):
        if abs(rates[k] - 2.0) < 1e-6:  # F_a
            label = f"F_a:λ={eigvals[k]:+.4f}"
            matched[label] = (eigvals[k], eigvecs[:, k])
        elif abs(rates[k] - 6.0) < 1e-6:  # F_b
            label = f"F_b:λ={eigvals[k]:+.4f}"
            matched[label] = (eigvals[k], eigvecs[:, k])
    return matched


def main() -> None:
    qs = [0.5, 1.0, 1.5, 2.0, 3.0]
    Ns = [6, 7, 8, 9, 10, 11, 12, 14, 16, 20]  # N ≥ N_block + 1 = 6
    print(f"# F89 path-4: AT-locked amplitude closed-form derivation")
    print(f"# N_block = {N_BLOCK} (path-4); qs = {qs}; Ns = {Ns}\n")

    gamma = 1.0

    # 1. Verify F_a/F_b eigvec q-independence
    print("## Step 1: q-independence of F_a/F_b eigenvectors")
    L0, basis0, _ = build_se_de(qs[0], gamma)
    matched0 = find_at_locked(L0, gamma, qs[0])
    print(f"# At q={qs[0]}: {len(matched0)} AT-locked modes found")

    # Check at other q values that the rate-2γ subspace is preserved
    for q in qs[1:]:
        L, basis, _ = build_se_de(q, gamma)
        matched = find_at_locked(L, gamma, q)
        same_count = (len(matched) == len(matched0))
        print(f"# At q={q}: {len(matched)} AT-locked  ({'count consistent' if same_count else 'COUNT CHANGED'})")
    print()

    # 2. Compute sigs[F_a/F_b] per mode at multiple (q, N)
    print("## Step 2: sigs · N²(N-1) per AT-locked mode at q=1.5")
    q = 1.5
    J = q * gamma
    L, basis, _ = build_se_de(J, gamma)
    matched = find_at_locked(L, gamma, J)
    w = per_site_reduction(basis)

    print(f"\n# {len(matched)} AT-locked modes. Per-mode sigs · N²(N-1) at q=1.5:\n")
    print("| Eigenvalue λ | rate | freq | N=6 | N=7 | N=11 | N=20 |")
    print("|---|---|---|---|---|---|---|")
    sigs_data = {}
    for label, (lam, v) in matched.items():
        v = v / np.linalg.norm(v)
        Mv = w @ v
        row = [str(label), f"{-lam.real:.2f}γ", f"{lam.imag/J:+.4f}J"]
        sigs_per_N = []
        for N in Ns:
            pre = np.sqrt(2 / (N**2 * (N - 1)))
            rho_24 = np.full(len(basis), pre / 2, dtype=complex)
            c = np.vdot(v, rho_24)
            sigs = abs(c) ** 2 * np.sum(np.abs(Mv) ** 2)
            scaled = sigs * (N * N * (N - 1))
            sigs_per_N.append(scaled)
        sigs_data[label] = (lam, sigs_per_N)
        # Show subset
        row_subset = [row[0], row[1], row[2]] + [f"{sigs_per_N[i]:.4e}" for i in [0, 1, 5, 9]]
        print("| " + " | ".join(row_subset) + " |")

    # 3. Sum F_a, sum F_b — looks for clean rational
    print("\n## Step 3: Sum sigs · N²(N-1) over all F_a / all F_b")
    fa_sum_per_N = [0.0] * len(Ns)
    fb_sum_per_N = [0.0] * len(Ns)
    for label, (lam, sigs_per_N) in sigs_data.items():
        rate = -lam.real
        for i, s in enumerate(sigs_per_N):
            if abs(rate - 2.0) < 1e-6:
                fa_sum_per_N[i] += s
            elif abs(rate - 6.0) < 1e-6:
                fb_sum_per_N[i] += s
    print(f"# Σ F_a sigs·N²(N-1): {[f'{s:.6f}' for s in fa_sum_per_N]}")
    print(f"# Σ F_b sigs·N²(N-1): {[f'{s:.6f}' for s in fb_sum_per_N]}")

    # Check constancy across N
    fa_const = all(abs(fa_sum_per_N[i] - fa_sum_per_N[0]) < 1e-10 for i in range(len(Ns)))
    fb_const = all(abs(fb_sum_per_N[i] - fb_sum_per_N[0]) < 1e-10 for i in range(len(Ns)))
    print(f"# F_a sum N-independent: {fa_const} (mean: {np.mean(fa_sum_per_N):.10f})")
    print(f"# F_b sum N-independent: {fb_const} (mean: {np.mean(fb_sum_per_N):.10f})")

    # 4. Symbolic recognition
    print("\n## Step 4: Symbolic recognition of constants")
    if fa_const:
        val = np.mean(fa_sum_per_N)
        for ext in [None, [sp.sqrt(2)], [sp.sqrt(3)], [sp.sqrt(5)], [sp.sqrt(2), sp.sqrt(3)]]:
            label = "rational" if ext is None else f"with {ext}"
            try:
                rat = sp.nsimplify(val, ext if ext else [], rational=False if ext else True, tolerance=1e-9)
                if not rat.has(sp.Float):
                    print(f"#   F_a sum {label}: {rat} ≈ {float(rat):.10f}  (input {val:.10f})")
            except Exception:
                pass
    if fb_const:
        val = np.mean(fb_sum_per_N)
        for ext in [None, [sp.sqrt(2)], [sp.sqrt(3)], [sp.sqrt(5)], [sp.sqrt(2), sp.sqrt(3)]]:
            label = "rational" if ext is None else f"with {ext}"
            try:
                rat = sp.nsimplify(val, ext if ext else [], rational=False if ext else True, tolerance=1e-9)
                if not rat.has(sp.Float):
                    print(f"#   F_b sum {label}: {rat} ≈ {float(rat):.10f}  (input {val:.10f})")
            except Exception:
                pass

    # 5. Per-mode symbolic recognition
    print("\n## Step 5: Per-mode symbolic recognition")
    for label, (lam, sigs_per_N) in sigs_data.items():
        if all(abs(sigs_per_N[i] - sigs_per_N[0]) < 1e-10 for i in range(len(Ns))):
            val = np.mean(sigs_per_N)
            print(f"# {label}: sigs · N²(N-1) = {val:.10f} (constant in N)")
            for ext in [None, [sp.sqrt(2)], [sp.sqrt(3)], [sp.sqrt(5)], [sp.sqrt(2), sp.sqrt(3)]]:
                ext_label = "rational" if ext is None else f"with {ext}"
                try:
                    rat = sp.nsimplify(val, ext if ext else [], rational=False if ext else True, tolerance=1e-9)
                    if not rat.has(sp.Float):
                        print(f"#     {ext_label}: {rat}")
                        break
                except Exception:
                    pass


if __name__ == "__main__":
    main()
