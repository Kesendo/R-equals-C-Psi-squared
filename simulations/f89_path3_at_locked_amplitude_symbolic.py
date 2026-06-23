"""F89 path-3: closed-form for the F_a (AT-locked rate 2γ) mode amplitudes.

Strategy: numerical eigvec extraction at multiple q values + nsimplify with √5.

Key claims tested:
  1. F_a, F_b eigenvectors are q-INDEPENDENT (only frequency depends on q;
     overlap-only support is geometric/topological, q-free).
  2. Inner product ⟨v_F_a | ρ_block(0)|_(SE,DE)⟩ is q-independent rational
     in N (with √5 from F_a structure).
  3. Per-site reduced amplitude sums to a clean closed form.
"""

from __future__ import annotations

import sys

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_path3_se_de(J: float, gamma: float):
    """24-dim (SE, DE) L_super on path-3 (4-qubit block)."""
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


def find_at_locked_eigvecs(L, gamma, J, sqrt5_target_sign):
    """Find F_a eigenvectors with eigenvalue -2γ + iJ·sqrt5_target_sign·(... ±√5).

    Returns dict {label: (lambda, eigvec_24)}.
    """
    eigvals, eigvecs = np.linalg.eig(L)
    sqrt5 = np.sqrt(5)
    targets = {
        "F_a:λ=−2γ+iJ(−1+√5)": -2 * gamma + 1j * J * (-1 + sqrt5),
        "F_a:λ=−2γ+iJ(−1−√5)": -2 * gamma + 1j * J * (-1 - sqrt5),
        "F_b:λ=−6γ+iJ(−1+√5)": -6 * gamma + 1j * J * (-1 + sqrt5),
        "F_b:λ=−6γ+iJ(−1−√5)": -6 * gamma + 1j * J * (-1 - sqrt5),
    }
    matched = {}
    for label, target in targets.items():
        k = np.argmin(np.abs(eigvals - target))
        matched[label] = (eigvals[k], eigvecs[:, k])
    return matched


def compare_eigvecs_across_q(qs):
    """Verify F_a, F_b eigenvectors are q-independent (modulo phase + degenerate-subspace rotation)."""
    print("## Are F_a, F_b eigenvectors q-independent?")
    print("# Compare overlap-pattern of eigvec at q=q1 vs q=q2 (gauge-fix by aligning largest-magnitude entry phase)\n")
    gamma = 1.0
    J0 = qs[0] * gamma
    L0, basis0, _ = build_path3_se_de(J0, gamma)
    matched0 = find_at_locked_eigvecs(L0, gamma, J0, +1)

    overlap_mask = np.array([i in jk for i, jk in basis0])

    for label, (lam0, v0) in matched0.items():
        # Gauge-fix: rotate so largest-magnitude entry is real positive
        idx_max = np.argmax(np.abs(v0))
        v0 = v0 * np.conj(v0[idx_max] / np.abs(v0[idx_max]))
        v0 = v0 / np.linalg.norm(v0)

        same_pattern = True
        deviations = []
        for q in qs[1:]:
            J = q * gamma
            L, basis, _ = build_path3_se_de(J, gamma)
            matched = find_at_locked_eigvecs(L, gamma, J, +1)
            v = matched[label][1]
            v = v * np.conj(v[idx_max] / np.abs(v[idx_max]))
            v = v / np.linalg.norm(v)
            # Compare overlap-pattern entries (real part since phase-aligned)
            diff = np.max(np.abs(np.abs(v) - np.abs(v0)))
            deviations.append((q, diff))
            if diff > 1e-6:
                same_pattern = False
        max_dev = max(d for _, d in deviations)
        print(f"#   {label}: overlap-pattern stable across q={qs[1:]} (max |Δv| = {max_dev:.2e})  "
              f"{'✓ q-independent' if same_pattern else '✗ q-dependent'}")
        # Show the overlap-pattern
        ov_pattern = sorted(np.abs(v0[overlap_mask]).round(4).tolist(), reverse=True)
        no_pattern = sorted(np.abs(v0[~overlap_mask]).round(4).tolist(), reverse=True)
        print(f"#     overlap |v| pattern (sorted): {ov_pattern}")
        print(f"#     no-overlap |v| pattern (sorted): {no_pattern}")


def extract_inner_products(qs, Ns):
    """Compute |⟨v_F | ρ_block(0)⟩|² for AT-locked F_a, F_b at multiple (q, N) values."""
    print("\n## Inner products |⟨v_F_a | ρ_block(0)|_(SE,DE)⟩|² across q, N")
    gamma = 1.0
    sqrt5 = np.sqrt(5)

    for q in qs:
        J = q * gamma
        L, basis, _ = build_path3_se_de(J, gamma)
        matched = find_at_locked_eigvecs(L, gamma, J, +1)
        print(f"\n### q = {q}")
        print("| N | |⟨F_a:E_2|ρ⟩|² | |⟨F_a:E_4|ρ⟩|² | |⟨F_b:E_2|ρ⟩|² | |⟨F_b:E_4|ρ⟩|² | sum F_a (E_2+E_4) |")
        print("|---|---|---|---|---|---|")
        for N in Ns:
            # ρ_block(0) (SE, DE) part: each entry = pre/2, where pre = √(2/(N²(N-1)))
            pre = np.sqrt(2 / (N**2 * (N - 1)))
            rho_24 = np.full(len(basis), pre / 2, dtype=complex)
            row = [N]
            sum_fa = 0.0
            for label, (lam, v) in matched.items():
                v = v / np.linalg.norm(v)
                ip = np.vdot(v, rho_24)
                ip_sq = abs(ip) ** 2
                row.append(f"{ip_sq:.4e}")
                if "F_a" in label:
                    sum_fa += ip_sq
            row.append(f"{sum_fa:.4e}")
            print("| " + " | ".join(map(str, row)) + " |")


def fit_amplitude_in_n(Ns, vals, label):
    """Fit vals[i] · N²(N−1) = polynomial(N), report rationalized coefficients."""
    n_arr = np.array(Ns, dtype=float)
    v_arr = np.array(vals, dtype=float)
    target = v_arr * (n_arr * n_arr * (n_arr - 1))
    for degree in range(0, 5):
        coefs = np.polynomial.polynomial.polyfit(n_arr, target, degree)
        residual = target - np.polynomial.polynomial.polyval(n_arr, coefs)
        rel_err = np.max(np.abs(residual)) / max(1e-15, np.max(np.abs(target)))
        if rel_err < 1e-9:
            print(f"# {label}: amplitude · N²(N−1) = polynomial degree {degree}, rel_err {rel_err:.2e}")
            for i, c in enumerate(coefs):
                # Try to recognize as a + b·√5 over integers
                rat = sp.nsimplify(c, [sp.sqrt(5)], rational=False, tolerance=1e-8)
                print(f"#   coef of N^{i}: numeric {c:.10f}  →  symbolic {rat}")
            return coefs, degree
    print(f"# {label}: no polynomial fit ≤ degree 4")
    return None, None


def fit_amplitude_with_sqrt5(qs, Ns):
    """Fit |⟨v_F_a|ρ⟩|² · N²(N−1) = polynomial in N with coefficients in Q[√5]."""
    print("\n## Symbolic fit: A_F_a(N, q) · N²(N−1) = polynomial in (N, q, √5)?\n")
    gamma = 1.0

    for q in qs:
        J = q * gamma
        L, basis, _ = build_path3_se_de(J, gamma)
        matched = find_at_locked_eigvecs(L, gamma, J, +1)
        pre = np.array([np.sqrt(2 / (N**2 * (N - 1))) for N in Ns])

        print(f"\n### q = {q}")
        for label, (lam, v) in matched.items():
            v = v / np.linalg.norm(v)
            vals = []
            for N, pre_N in zip(Ns, pre):
                rho_24 = np.full(len(basis), pre_N / 2, dtype=complex)
                ip = np.vdot(v, rho_24)
                vals.append(abs(ip) ** 2)
            fit_amplitude_in_n(Ns, vals, label)


def per_site_reduction_within_block_se_de(basis, n_block=4):
    """Per-site reduction matrix w[l, idx] for (SE, DE) basis pairs.

    Acts on a length-24 (SE, DE) flat-vector v: (w @ v)[l] = (single-site reduced
    block ρ at site l)_{0,1}, computed by tracing out the other (n_block-1) sites.

    For (SE, DE) basis pair (i, jk): w[l] is non-zero only when l ∈ {j, k} AND
    SE site i = the OTHER element of {j, k} (i.e., overlap pair with l ≠ i).
    """
    w = np.zeros((n_block, len(basis)), dtype=float)
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue  # would not contribute (need bit pattern: row 0_l, col 1_l, others match)
            other_in_jk = jk[1] if jk[0] == l else jk[0]
            if other_in_jk == i:
                w[l, idx] = 1.0
    return w


def compute_sigs_per_F(qs, Ns):
    """Compute sigs[F_a/b] = |c[F]|² · Σ_l |w[l] · v_F|² for each AT-locked mode.

    This is the actual contribution to S(t) from each (SE, DE) F_a/b mode.
    """
    print("\n## Per-site-reduced amplitude sigs[F] = |c|²·Σ_l|w·v|² (S(t) contribution)")
    print("# Should match the (SE, DE) component of the path-3 numerical multi-exp at "
          "(rate 2γ, |freq|/J ∈ {√5−1, √5+1}) modes.")
    gamma = 1.0
    sqrt5 = np.sqrt(5)

    for q in qs:
        J = q * gamma
        L, basis, _ = build_path3_se_de(J, gamma)
        matched = find_at_locked_eigvecs(L, gamma, J, +1)
        w = per_site_reduction_within_block_se_de(basis, n_block=4)

        print(f"\n### q = {q}")
        print("| N | sigs[F_a:E_2] | sigs[F_a:E_4] | sigs[F_b:E_2] | sigs[F_b:E_4] |")
        print("|---|---|---|---|---|")
        for N in Ns:
            pre = np.sqrt(2 / (N**2 * (N - 1)))
            rho_24 = np.full(len(basis), pre / 2, dtype=complex)
            row = [N]
            for label, (lam, v) in matched.items():
                v = v / np.linalg.norm(v)
                # c[F] = ⟨v | ρ_0⟩ for L_oo restricted (anti-Hermitian within overlap)
                c = np.vdot(v, rho_24)
                # M[l, F] = (w · v)[l] = per-site reduction of v
                Mv = w @ v
                # sigs[F] = |c|² · Σ_l |M[l, F]|²
                sig = abs(c) ** 2 * np.sum(np.abs(Mv) ** 2)
                row.append(f"{sig:.4e}")
            print("| " + " | ".join(map(str, row)) + " |")


def verify_closed_form_sigs(qs, Ns):
    """Verify closed-form sigs[F_a:E_n] = (33 ± 14√5) / [9·N²(N−1)]."""
    print("\n## Closed-form verification:")
    print("#   sigs[F_a:E_2](N) = (33 + 14√5) / [9·N²(N−1)]")
    print("#   sigs[F_a:E_4](N) = (33 − 14√5) / [9·N²(N−1)]")
    print("#   Sum F_a = 22 / [3·N²(N−1)]   (rational, √5 cancels)")
    sqrt5 = np.sqrt(5)
    gamma = 1.0
    J = 1.5 * gamma  # arbitrary q (verified q-independent)
    L, basis, _ = build_path3_se_de(J, gamma)
    matched = find_at_locked_eigvecs(L, gamma, J, +1)
    w = per_site_reduction_within_block_se_de(basis, n_block=4)

    print("\n| N | sigs_E2 (numerical) | (33+14√5)/[9N²(N-1)] | diff | sigs_E4 (numerical) | (33-14√5)/[9N²(N-1)] | diff | sum | 22/[3N²(N-1)] |")
    print("|---|---|---|---|---|---|---|---|---|")
    for N in Ns:
        pre = np.sqrt(2 / (N**2 * (N - 1)))
        rho_24 = np.full(len(basis), pre / 2, dtype=complex)
        denom = N * N * (N - 1)
        pred_e2 = (33 + 14 * sqrt5) / (9 * denom)
        pred_e4 = (33 - 14 * sqrt5) / (9 * denom)
        sum_pred = 22 / (3 * denom)

        v_e2 = matched["F_a:λ=−2γ+iJ(−1+√5)"][1]
        v_e2 = v_e2 / np.linalg.norm(v_e2)
        sigs_e2 = abs(np.vdot(v_e2, rho_24)) ** 2 * np.sum(np.abs(w @ v_e2) ** 2)

        v_e4 = matched["F_a:λ=−2γ+iJ(−1−√5)"][1]
        v_e4 = v_e4 / np.linalg.norm(v_e4)
        sigs_e4 = abs(np.vdot(v_e4, rho_24)) ** 2 * np.sum(np.abs(w @ v_e4) ** 2)

        sum_num = sigs_e2 + sigs_e4
        print(f"| {N} | {sigs_e2:.5e} | {pred_e2:.5e} | {sigs_e2 - pred_e2:+.2e} | "
              f"{sigs_e4:.5e} | {pred_e4:.5e} | {sigs_e4 - pred_e4:+.2e} | "
              f"{sum_num:.5e} | {sum_pred:.5e} |")


def main() -> None:
    qs = [0.5, 1.0, 1.5, 2.0, 3.0]
    Ns = [5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20]

    print("# F89 path-3: AT-locked amplitude closed-form derivation")
    print(f"# qs = {qs}")
    print(f"# Ns = {Ns}\n")

    compare_eigvecs_across_q(qs)
    extract_inner_products(qs, Ns)
    fit_amplitude_with_sqrt5(qs, Ns)
    compute_sigs_per_F(qs, Ns)
    verify_closed_form_sigs(qs, Ns)


if __name__ == "__main__":
    main()
