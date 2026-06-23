"""F89 path-3: octic-mode amplitude q-dependence symbolic probe.

For each of the 8 octic modes, the per-mode sigs(N) follows const(q)/[N²(N−1)]
empirically (verified at q=1.5). Question: what is const(q)?

Strategy:
  1. Sample sigs[mode] at multiple q values (q ∈ {0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 5.0}).
  2. For each q: identify the 8 octic modes by exclusion (all (SE, DE) S_2-sym
     eigenvalues except the 4 closed-form F_a/F_b roots).
  3. Track each mode across q (modes shift in rate/freq with q).
  4. For each "tracked" mode, fit sigs(q) against rational functions of (q, √5).
  5. Hamming-complement pair sums (Γ_a + Γ_b = 8γ) — check if pair sigs sum is
     a clean function (might be cleaner than individual modes).
"""

from __future__ import annotations

import sys

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_path3_se_de(J: float, gamma: float):
    """24-dim (SE, DE) L_super on path-3."""
    de_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    basis = [(i, jk) for i in range(4) for jk in de_pairs]
    n = len(basis)

    M_SE = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J

    M_DE = np.zeros((6, 6))
    for idx, (j, k) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j <= 3 and new_j != k:
                new_pair = tuple(sorted([new_j, k]))
                if abs(new_j - j) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2 * J
        for new_k in [k - 1, k + 1]:
            if 0 <= new_k <= 3 and new_k != j:
                new_pair = tuple(sorted([j, new_k]))
                if abs(new_k - k) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2 * J

    L = np.zeros((n, n), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(4):
            if M_SE[i2, i] != 0:
                idx2 = basis.index((i2, jk))
                L[idx2, idx] += -1j * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(6):
            if M_DE[jk_idx, jk2_idx] != 0:
                jk2 = de_pairs[jk2_idx]
                idx2 = basis.index((i, jk2))
                L[idx2, idx] += 1j * M_DE[jk_idx, jk2_idx]
        L[idx, idx] += -2 * gamma if i in jk else -6 * gamma

    return L, basis, de_pairs


def build_s2_sym_projector(basis):
    perm_se = {0: 3, 1: 2, 2: 1, 3: 0}
    perm_de = lambda jk: tuple(sorted([perm_se[jk[0]], perm_se[jk[1]]]))
    cols = []
    handled = set()
    for idx, (i, jk) in enumerate(basis):
        if idx in handled:
            continue
        idx2 = basis.index((perm_se[i], perm_de(jk)))
        v = np.zeros(len(basis), dtype=complex)
        if idx == idx2:
            v[idx] = 1.0
        else:
            v[idx] = 1.0 / np.sqrt(2)
            v[idx2] = 1.0 / np.sqrt(2)
            handled.add(idx2)
        cols.append(v)
        handled.add(idx)
    return np.column_stack(cols)


def per_site_reduction(basis, n_block=4):
    w = np.zeros((n_block, len(basis)), dtype=float)
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue
            other = jk[1] if jk[0] == l else jk[0]
            if other == i:
                w[l, idx] = 1.0
    return w


def octic_modes_sigs(q: float, gamma: float = 1.0):
    """Return list of (rate Γ/γ, |freq|/J, sigs · N²(N-1)) for the 8 octic modes
    (S_2-sym 12-dim sub-block minus 4 AT-locked F_a/F_b)."""
    J = q * gamma
    L, basis, de_pairs = build_path3_se_de(J, gamma)
    P = build_s2_sym_projector(basis)
    L_sym = P.conj().T @ L @ P  # 12×12

    eigvals, eigvecs_12 = np.linalg.eig(L_sym)
    rates = -eigvals.real / gamma
    freqs = np.abs(eigvals.imag) / J

    is_at_locked = (np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6)
    octic_idx = np.where(~is_at_locked)[0]
    assert len(octic_idx) == 8, f"Expected 8 octic modes, got {len(octic_idx)} (q={q})"

    # Lift S_2-sym eigvecs to 24-dim (SE, DE) basis
    eigvecs_24 = P @ eigvecs_12

    N = 11
    pre = np.sqrt(2 / (N**2 * (N - 1)))
    rho_24 = np.full(len(basis), pre / 2, dtype=complex)
    w = per_site_reduction(basis)

    results = []
    for k in octic_idx:
        v = eigvecs_24[:, k]
        v = v / np.linalg.norm(v)
        c = np.vdot(v, rho_24)
        Mv = w @ v
        sigs = abs(c) ** 2 * np.sum(np.abs(Mv) ** 2)
        sigs_scaled = sigs * (N * N * (N - 1))
        results.append((rates[k], freqs[k], sigs_scaled))
    return results


def track_modes_across_q(qs):
    """Track 8 octic modes across q values by closest (rate, freq) matching."""
    print("## Octic-mode sigs · N²(N-1) tracked across q (8 modes)\n")

    # Get q=1.5 modes as reference
    ref = octic_modes_sigs(1.5)
    ref_sorted = sorted(ref, key=lambda x: x[0])
    print("# Reference (q=1.5):")
    for i, (r, f, s) in enumerate(ref_sorted):
        print(f"#   mode_{i}: rate={r:.4f}, |freq|/J={f:.4f}, sigs·N²(N-1)={s:.6e}")
    print()

    print("## sigs · N²(N-1) for each mode at varying q:")
    print("| q | mode_0 | mode_1 | mode_2 | mode_3 | mode_4 | mode_5 | mode_6 | mode_7 | sum |")
    print("|---|" + "|".join(["---"] * 9) + "|")

    all_data = {}
    for q in qs:
        modes = octic_modes_sigs(q)
        # Sort by rate to match indexing across q
        # Note: mode-tracking is fragile when rates cross; this works only for "monotonically separated" modes
        modes_sorted = sorted(modes, key=lambda x: x[0])
        row = [f"{q}"]
        sigs_sum = 0
        for r, f, s in modes_sorted:
            row.append(f"{s:.4e}")
            sigs_sum += s
        row.append(f"{sigs_sum:.4e}")
        print("| " + " | ".join(row) + " |")
        all_data[q] = modes_sorted
    return all_data


def hamming_complement_pair_sums(qs):
    """For each q, pair octic modes by Γ_a + Γ_b = 8γ at fixed |freq|, sum sigs."""
    print("\n## Hamming-complement pair sigs sum (paired by Γ_a + Γ_b = 8γ at same |freq|/J):")
    print("# Pair structure: (3.349, 4.651), (3.599, 4.401), (3.777, 4.223), (4.0, 4.0) — these are q-shifted.\n")

    print("| q | pair_1 sum | pair_2 sum | pair_3 sum | pair_4 sum | total |")
    print("|---|---|---|---|---|---|")
    for q in qs:
        modes = octic_modes_sigs(q)
        # Group by |freq|/J (round to 4 decimals to handle q-shift mismatch)
        from collections import defaultdict
        groups = defaultdict(list)
        for r, f, s in modes:
            groups[round(f, 3)].append((r, s))
        # Pair within each |freq| group
        pair_sums = []
        for f, items in groups.items():
            if len(items) == 2:
                pair_sums.append(sum(s for _, s in items))
            elif len(items) == 1:
                # Singleton — could be a mode that doesn't pair (rate 4γ self-pair)
                pair_sums.append(items[0][1])
            else:
                pair_sums.append(sum(s for _, s in items))  # multi-paired
        # Sort and print first 4 pair sums
        pair_sums_sorted = sorted(pair_sums, reverse=True)
        row = [f"{q}"] + [f"{p:.4e}" for p in pair_sums_sorted[:4]] + [f"{sum(pair_sums):.4e}"]
        print("| " + " | ".join(row) + " |")


def fit_total_octic_sigs_in_q(qs):
    """Try to fit Σ_8octic sigs · N²(N-1) as polynomial in q."""
    print("\n## Sum of all 8 octic mode sigs · N²(N-1) vs q — polynomial fit?")
    qs_arr = np.array(qs, dtype=float)
    sigs_sums = np.array([sum(s for _, _, s in octic_modes_sigs(q)) for q in qs])
    print("| q | Σ sigs · N²(N-1) |")
    print("|---|---|")
    for q, s in zip(qs, sigs_sums):
        print(f"| {q} | {s:.6e} |")

    print("\n## Polynomial fit of Σ vs q:")
    for degree in range(0, 6):
        coefs = np.polynomial.polynomial.polyfit(qs_arr, sigs_sums, degree)
        residual = sigs_sums - np.polynomial.polynomial.polyval(qs_arr, coefs)
        rel_err = np.max(np.abs(residual)) / max(1e-15, np.max(np.abs(sigs_sums)))
        if rel_err < 1e-10:
            print(f"# Degree {degree} fit: rel_err = {rel_err:.2e}")
            for i, c in enumerate(coefs):
                rat = sp.nsimplify(c, [sp.sqrt(5)], rational=False, tolerance=1e-8)
                print(f"#   coef of q^{i}: numeric {c:.10f}  →  symbolic {rat}")
            break
    else:
        print(f"# No polynomial fit ≤ degree 5 (max rel_err: {rel_err:.2e})")
        # Try rational function fits
        for power in [-2, -1, 1, 2]:
            transformed = sigs_sums * qs_arr**(-power) if power < 0 else sigs_sums * qs_arr**power
            print(f"# Σ · q^({-power}): {transformed}")


def main() -> None:
    qs = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 5.0]
    print("# F89 path-3 octic-mode amplitude q-dependence probe")
    print(f"# qs = {qs}\n")

    track_modes_across_q(qs)
    hamming_complement_pair_sums(qs)
    fit_total_octic_sigs_in_q(qs)


if __name__ == "__main__":
    main()
