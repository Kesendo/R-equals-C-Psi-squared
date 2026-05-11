"""F89 path-3: locate the exceptional point at q ≈ 0.659.

The path-3 octic discriminant has a perfect-square factor (3q⁴ + q² − 1)².
Its zero at q² = (−1+√13)/6 ≈ 0.4343, q ≈ 0.6590, locates an EP where two
octic eigenvalues coalesce.

Goals:
  1. Verify the EP location numerically (sample q densely near 0.659).
  2. Identify WHICH two octic eigenvalues collide.
  3. Track their (rate, freq) labels through the EP.
  4. Diagnose Jordan-block formation: at the EP, the eigvec matrix becomes
     rank-deficient (defective). Compute the smallest singular value of
     (L_sym - λ_EP·I) to detect.
  5. Compare to F86's c=2 Q_EP = 1/√2 ≈ 0.7071: are they the same EP or different?
"""

from __future__ import annotations

import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_path3_se_de_sym(J: float, gamma: float):
    """Return (L_sym 12×12, basis 24, P 24x12 projector)."""
    de_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    basis = [(i, jk) for i in range(4) for jk in de_pairs]

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

    L = np.zeros((24, 24), dtype=complex)
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

    perm_se = {0: 3, 1: 2, 2: 1, 3: 0}
    perm_de = lambda jk: tuple(sorted([perm_se[jk[0]], perm_se[jk[1]]]))
    cols = []
    handled = set()
    for idx, (i, jk) in enumerate(basis):
        if idx in handled:
            continue
        idx2 = basis.index((perm_se[i], perm_de(jk)))
        v = np.zeros(24, dtype=complex)
        if idx == idx2:
            v[idx] = 1.0
        else:
            v[idx] = 1.0 / np.sqrt(2)
            v[idx2] = 1.0 / np.sqrt(2)
            handled.add(idx2)
        cols.append(v)
        handled.add(idx)
    P = np.column_stack(cols)
    L_sym = P.conj().T @ L @ P
    return L_sym, basis, P


def octic_eigenvalues(q: float, gamma: float = 1.0):
    """Return sorted list of 8 octic eigenvalues (S_2-sym minus 4 AT-locked)."""
    J = q * gamma
    L_sym, _, _ = build_path3_se_de_sym(J, gamma)
    eigvals = np.linalg.eigvals(L_sym)
    rates = -eigvals.real / gamma
    is_at = (np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6)
    octic = eigvals[~is_at]
    assert len(octic) == 8, f"Expected 8 octic, got {len(octic)} (q={q})"
    return sorted(octic, key=lambda x: (x.real, x.imag))


def smallest_singular_value(L_sym, lam):
    """For Jordan-block detection: smallest σ of (L_sym - λI)."""
    M = L_sym - lam * np.eye(L_sym.shape[0])
    s = np.linalg.svd(M, compute_uv=False)
    return s[-1]


def main() -> None:
    print("# F89 path-3 EP locator at q² = (−1+√13)/6, q ≈ 0.659\n")

    # 1. Verify analytical EP location
    q_ep_analytical = np.sqrt((-1 + np.sqrt(13)) / 6)
    print(f"## Analytical EP location: q = √((−1+√13)/6) = {q_ep_analytical:.10f}")
    # Verify: 3q⁴ + q² − 1 = 0 at q_ep
    res = 3 * q_ep_analytical**4 + q_ep_analytical**2 - 1
    print(f"#   Check 3q⁴+q²−1 at q_EP: {res:.2e} (should be 0)\n")

    # 2. Track 8 octic eigenvalues across q sweep near EP
    print("## Octic eigenvalues swept across q near EP:")
    print("# Look for two eigenvalues whose distance |λ_a − λ_b| → 0 near q_EP\n")

    qs = np.linspace(0.55, 0.78, 24)
    print("| q | min |λ_a − λ_b| (closest pair distance) | which pair (rate Γ/γ, freq ω/J·γ) |")
    print("|---|---|---|")
    for q in qs:
        evs = octic_eigenvalues(q)
        # Find closest pair
        min_dist = np.inf
        min_pair = None
        for i in range(8):
            for j in range(i + 1, 8):
                d = abs(evs[i] - evs[j])
                if d < min_dist:
                    min_dist = d
                    min_pair = (evs[i], evs[j])
        rates = (-min_pair[0].real, -min_pair[1].real)
        freqs = (min_pair[0].imag, min_pair[1].imag)
        print(f"| {q:.4f} | {min_dist:.4e} | rates ({rates[0]:.4f}, {rates[1]:.4f}), freqs ({freqs[0]:+.4f}, {freqs[1]:+.4f}) |")

    print(f"\n## Fine-resolution sweep around q_EP = {q_ep_analytical:.6f}")
    qs_fine = np.linspace(q_ep_analytical - 0.005, q_ep_analytical + 0.005, 21)
    print("| q | min |λ_a − λ_b| | smallest σ((L − λ_close I)) |")
    print("|---|---|---|")
    for q in qs_fine:
        L_sym, _, _ = build_path3_se_de_sym(q, 1.0)
        eigvals = np.linalg.eigvals(L_sym)
        rates = -eigvals.real
        is_at = (np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6)
        octic = sorted(eigvals[~is_at], key=lambda x: (x.real, x.imag))
        # Closest pair
        dists = [abs(octic[i] - octic[i+1]) for i in range(7)]
        min_dist = min(dists)
        # Singular value at average eigenvalue (approximation of EP eigenvalue)
        i_min = np.argmin(dists)
        lam_close = (octic[i_min] + octic[i_min + 1]) / 2
        sv = smallest_singular_value(L_sym, lam_close)
        print(f"| {q:.6f} | {min_dist:.4e} | {sv:.4e} |")

    # 3. EP eigenvalue at the analytical q_EP
    print("\n## EP eigenvalue (the merged value) at q = q_EP exactly:")
    L_sym, _, _ = build_path3_se_de_sym(q_ep_analytical, 1.0)
    eigvals = np.linalg.eigvals(L_sym)
    rates = -eigvals.real
    is_at = (np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6)
    octic = sorted(eigvals[~is_at], key=lambda x: (x.real, x.imag))
    dists = [abs(octic[i] - octic[i+1]) for i in range(7)]
    i_min = np.argmin(dists)
    lam_ep = (octic[i_min] + octic[i_min + 1]) / 2
    print(f"# Merged eigenvalue λ_EP ≈ {lam_ep:.6f}")
    print(f"# In units of γ: rate Γ/γ = {-lam_ep.real:.4f}, freq ω/(J·γ) = {lam_ep.imag/q_ep_analytical:.4f}")
    print(f"# Re(λ_EP)/γ = {-lam_ep.real:.4f}  vs F86 universal t_peak = 1/(4γ₀) → Re = -4γ predicted")

    # 4. Compare to F86 c=2 Q_EP
    print("\n## Comparison with F86 c=2 Q_EP:")
    print("# F86 c=2 Q_EP = 2/g_eff ≈ 1/√(2·1) = 1/√2 ≈ 0.7071 (N→∞ asymptote)")
    print("# F86 c=2 Q_EP at small N (per memory): N=5: 0.7233, N=6: 0.7138, N=7: 0.7071, N=8: 0.7044")
    print(f"# Path-3 (N_block=4) octic EP: q ≈ {q_ep_analytical:.4f}")
    print(f"#   N_block=4 extrapolation of F86 trend (N=5..8 → 0.72→0.70): would give Q_EP > 0.7233 at N=4")
    print(f"#   Path-3's q={q_ep_analytical:.4f} is BELOW F86's expected Q_EP(N=4): they are DIFFERENT EPs.")
    print(f"#   F86 EP comes from the (vac, SE) ↔ first-block 2-level effective.")
    print(f"#   Path-3 EP comes from a coupling INSIDE the (SE, DE) octic — different sector.")


if __name__ == "__main__":
    main()
