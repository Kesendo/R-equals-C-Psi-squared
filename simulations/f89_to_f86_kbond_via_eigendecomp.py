"""F89→F86 direct bridge: compute F86 K_b(Q, t_peak) via F89 (SE,DE) full eigendecomposition
+ per-bond Hellmann-Feynman. Test if F89's AT-locked F_a/F_b modes alone reproduce the
bare doubled-PTF SVD-block floor HWHM_left/Q_peak ≈ 0.6715, and if the full octic-included
computation reproduces the empirical 0.7506 (Interior) / 0.7728 (Endpoint).

The connection: F86 c=2 N qubits L_super on (n=1, n+1=2) coherence block IS F89 (SE,DE)
sub-block for path-(N-1). Same algebraic object.

c=2 N=5 → F89 path-4 (50-dim full, 26-dim S_2-sym, 4 AT-locked F_a/F_b at λ = -2γ ± 2iJ).
"""

from __future__ import annotations

import sys
from itertools import combinations

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

N = 8  # F86 c=2 N=8 ↔ F89 path-7 (N_block=8, dim 224)


def build_se_de_full(J: float, gamma: float, n_block: int):
    """24/50/90/147-dim L_super on (SE, DE) sector (uniform-J multi-bond)."""
    de_pairs = list(combinations(range(n_block), 2))
    basis = [(i, jk) for i in range(n_block) for jk in de_pairs]
    n = len(basis)

    M_SE = np.zeros((n_block, n_block))
    for a in range(n_block):
        for b in range(n_block):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J

    M_DE = np.zeros((len(de_pairs), len(de_pairs)))
    for idx, (j, k) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j < n_block and new_j != k:
                np_ = tuple(sorted([new_j, k]))
                if abs(new_j - j) == 1 and np_ in de_pairs:
                    M_DE[de_pairs.index(np_), idx] += 2 * J
        for new_k in [k - 1, k + 1]:
            if 0 <= new_k < n_block and new_k != j:
                np_ = tuple(sorted([j, new_k]))
                if abs(new_k - k) == 1 and np_ in de_pairs:
                    M_DE[de_pairs.index(np_), idx] += 2 * J

    L = np.zeros((n, n), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(n_block):
            if M_SE[i2, i] != 0:
                L[basis.index((i2, jk)), idx] += -1j * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(len(de_pairs)):
            if M_DE[jk_idx, jk2_idx] != 0:
                jk2 = de_pairs[jk2_idx]
                L[basis.index((i, jk2)), idx] += 1j * M_DE[jk_idx, jk2_idx]
        L[idx, idx] += -2 * gamma if i in jk else -6 * gamma
    return L, basis, de_pairs


def build_per_bond_M(b_bond: int, n_block: int, basis, de_pairs):
    """M_h_per_bond[b] = ∂L/∂J_b at J=1, i.e. -i[H_B^(b), ·] for the b-th bond.

    H_B^(b) = X_b X_{b+1} + Y_b Y_{b+1} = 2(σ⁺_b σ⁻_{b+1} + σ⁻_b σ⁺_{b+1}).

    Action on (SE, DE) basis pairs: same as M_SE, M_DE above but ONLY the b-th bond's
    hopping (sites b ↔ b+1).
    """
    n = len(basis)
    M_SE_b = np.zeros((n_block, n_block))
    M_SE_b[b_bond, b_bond + 1] = 2.0
    M_SE_b[b_bond + 1, b_bond] = 2.0

    M_DE_b = np.zeros((len(de_pairs), len(de_pairs)))
    for idx, (j, k) in enumerate(de_pairs):
        # Hop site b ↔ b+1 on the j-particle (if j == b or j == b+1)
        if j == b_bond:  # try j → b+1 (if b+1 != k)
            if b_bond + 1 != k:
                np_ = tuple(sorted([b_bond + 1, k]))
                if np_ in de_pairs:
                    M_DE_b[de_pairs.index(np_), idx] += 2.0
        if j == b_bond + 1:  # try j → b
            if b_bond != k:
                np_ = tuple(sorted([b_bond, k]))
                if np_ in de_pairs:
                    M_DE_b[de_pairs.index(np_), idx] += 2.0
        if k == b_bond:  # try k → b+1 (if b+1 != j)
            if b_bond + 1 != j:
                np_ = tuple(sorted([j, b_bond + 1]))
                if np_ in de_pairs:
                    M_DE_b[de_pairs.index(np_), idx] += 2.0
        if k == b_bond + 1:  # try k → b
            if b_bond != j:
                np_ = tuple(sorted([j, b_bond]))
                if np_ in de_pairs:
                    M_DE_b[de_pairs.index(np_), idx] += 2.0

    M_b = np.zeros((n, n), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(n_block):
            if M_SE_b[i2, i] != 0:
                M_b[basis.index((i2, jk)), idx] += -1j * M_SE_b[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(len(de_pairs)):
            if M_DE_b[jk_idx, jk2_idx] != 0:
                jk2 = de_pairs[jk2_idx]
                M_b[basis.index((i, jk2)), idx] += 1j * M_DE_b[jk_idx, jk2_idx]
    return M_b


def compute_kbond(Q: float, gamma: float, n_block: int, b_bond: int, t_grid):
    """K_b(Q, t) for given Q via full eigendecomposition of L = D + J·M_h_total.

    Probe ρ_0: the (n=1, n+1=2) Dicke 0-1 coherence — for the (SE, DE) sector this
    is the F89 ρ_block(0)|_(SE,DE) = uniform pre/2 entries on all basis pairs.

    S_kernel: spatial-sum kernel = Σ_l 2·|0_l⟩⟨1_l| reduced to (SE, DE) basis pairs.
    Per F89: w[l] picks (i, jk) with l ∈ jk and other element = i (overlap).
    """
    J = Q * gamma
    L, basis, de_pairs = build_se_de_full(J, gamma, n_block)
    n = len(basis)

    # Per-bond perturbation
    M_b = build_per_bond_M(b_bond, n_block, basis, de_pairs)

    # Eigendecomposition of L
    eigvals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)

    # Probe ρ_0: F89 Term-1 (vac, SE) Dicke probe restricted to (SE, DE) is just the
    # uniform vector with all entries = pre/2 where pre = √(2/(N²(N-1)))
    # For F86: use N for the FULL chain dim (probe is on full N-qubit Dicke)
    pre = np.sqrt(2 / (n_block * n_block * (n_block - 1)))
    probe = np.full(n, pre / 2, dtype=complex)

    # S-kernel: per-site reduction summed
    # (ρ_l)_{0,1} = Σ_{(i,jk) overlap, with l in jk and i = other(jk)} ρ[i, jk]
    # S(t) = Σ_l 2|...|² but K_b uses LINEAR not squared, so we need
    # ⟨ρ| S_kernel |∂ρ/∂J_b⟩ = Σ_l (ρ_l)_{0,1}* · (∂ρ_l/∂J_b)_{0,1} · 2 + cc
    # Equivalent: S_kernel = 2 Σ_l |w_l⟩⟨w_l| where w_l[idx] = 1 if (idx is overlap-pair with l)
    # Actually ⟨ρ| S |drho⟩ where S projects onto the per-site (ρ_l)_{0,1} sum:
    # K_b(t) = 2·Re(⟨ρ(t)| S |∂ρ/∂J_b⟩)
    # S = Σ_l |w_l⟩⟨w_l| (projector-like, but unnormalized; F89's per_site_reduction)
    w = np.zeros((n_block, n), dtype=float)
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue
            other = jk[1] if jk[0] == l else jk[0]
            if other == i:
                w[l, idx] = 1.0

    # c0 = R⁻¹ · probe
    c0 = Rinv @ probe
    # X_b = R⁻¹ · M_b · R
    X_b = Rinv @ M_b @ R

    k_max_over_t = 0.0
    t_at_peak = 0.0
    for t in t_grid:
        expLam = np.exp(eigvals * t)
        # ρ(t) = R · diag(expLam) · c0
        rho_t = R @ (expLam * c0)
        # ∂ρ/∂J_b via Duhamel: drho = R · fbC0 with fbC0[r] = Σ_c X_b[r,c] · iMat[r,c] · c0[c]
        iMat = np.zeros((n, n), dtype=complex)
        for r in range(n):
            for c in range(n):
                diff = eigvals[c] - eigvals[r]
                if abs(diff) > 1e-10:
                    iMat[r, c] = (expLam[c] - expLam[r]) / diff
                else:
                    iMat[r, c] = t * expLam[r]
        fbC0 = np.einsum("rc,rc,c->r", X_b, iMat, c0)
        drho = R @ fbC0
        # S_kernel · drho: each site l projects via w[l]; sum the |w_l⟩⟨w_l|·drho contribution
        # ⟨ρ(t)| S |drho⟩ = Σ_l ⟨ρ(t)|w_l⟩·⟨w_l|drho⟩
        rho_proj = w @ rho_t  # shape (n_block,)
        drho_proj = w @ drho  # shape (n_block,)
        inner = np.sum(np.conj(rho_proj) * drho_proj)
        kAbs = abs(2.0 * inner.real)
        if kAbs > k_max_over_t:
            k_max_over_t = kAbs
            t_at_peak = t
    return k_max_over_t, t_at_peak


def main() -> None:
    n_block = N
    gamma = 0.05
    t_peak = 1.0 / (4.0 * gamma)
    t_grid = np.linspace(0.6 * t_peak, 1.6 * t_peak, 21)
    q_grid = np.linspace(0.05, 10.00, 300)  # extended for orbit-escape at N≥7 (F89-J = F86-J/2)

    print(f"# F89→F86 K_b probe via (SE,DE) full eigendecomposition")
    print(f"# c=2 N={n_block} (= F89 path-{n_block-1}); γ={gamma}, t_peak={t_peak}")
    print(f"# Q-grid: 153 points [0.20, 4.00]; t-grid: 21 points [0.6, 1.6]·t_peak")
    print(f"# Empirical anchor: HWHM_left/Q_peak = 0.7506 (Interior), 0.7728 (Endpoint)\n")

    bonds = list(range(n_block - 1))  # bonds 0-1, 1-2, 2-3, 3-4 for N=5
    print(f"## Computing K_b(Q, t_peak) for {len(bonds)} bonds at {len(q_grid)} Q-points")
    print(f"# (Per-bond Hellmann-Feynman + full L eigendecomposition; ~{len(bonds)*len(q_grid)} K computations)")

    results = {}
    for b in bonds:
        k_curve = np.zeros(len(q_grid))
        for iQ, Q in enumerate(q_grid):
            k_curve[iQ], _ = compute_kbond(Q, gamma, n_block, b, t_grid)
        # Find Q_peak (parabolic refinement of grid max)
        iMax = int(np.argmax(k_curve))
        kMax = k_curve[iMax]
        if 0 < iMax < len(q_grid) - 1:
            # Parabolic fit
            x = q_grid[iMax - 1:iMax + 2]
            y = k_curve[iMax - 1:iMax + 2]
            # y = a(x - q*)² + k*; vertex at x* = -b/(2a) for ax²+bx+c
            coefs = np.polyfit(x, y, 2)
            q_peak = -coefs[1] / (2 * coefs[0])
        else:
            q_peak = q_grid[iMax]
        # HWHM_left: linear interp where k_curve crosses kMax/2 going up
        half = kMax / 2
        hwhm_left = None
        for i in range(iMax, 0, -1):
            if k_curve[i] >= half >= k_curve[i - 1]:
                # Linear interp
                frac = (half - k_curve[i - 1]) / (k_curve[i] - k_curve[i - 1])
                q_half = q_grid[i - 1] + frac * (q_grid[i] - q_grid[i - 1])
                hwhm_left = q_peak - q_half
                break
        ratio = hwhm_left / q_peak if hwhm_left is not None else None
        bond_class = "Endpoint" if (b == 0 or b == n_block - 2) else "Interior"
        results[b] = (q_peak, kMax, hwhm_left, ratio, bond_class)
        hwhm_str = f"{hwhm_left:.4f}" if hwhm_left is not None else "NaN"
        ratio_str = f"{ratio:.4f}" if ratio is not None else "NaN"
        print(f"# bond {b} ({bond_class}): Q_peak = {q_peak:.4f}, K_max = {kMax:.4e}, "
              f"HWHM_left = {hwhm_str}, ratio = {ratio_str}")

    print("\n## Class-mean HWHM_left/Q_peak ratios:")
    print("| Class | bonds | mean ratio | empirical anchor |")
    print("|---|---|---|---|")
    end_bonds = [b for b in bonds if results[b][4] == "Endpoint"]
    int_bonds = [b for b in bonds if results[b][4] == "Interior"]
    if end_bonds:
        end_mean = np.mean([results[b][3] for b in end_bonds if results[b][3] is not None])
        print(f"| Endpoint | {end_bonds} | {end_mean:.4f} | 0.7728 |")
    if int_bonds:
        int_mean = np.mean([results[b][3] for b in int_bonds if results[b][3] is not None])
        print(f"| Interior | {int_bonds} | {int_mean:.4f} | 0.7506 |")

    print(f"\n## Verdict")
    print(f"# If the Endpoint/Interior class-mean matches 0.7728/0.7506 to ≤ 0.005:")
    print(f"#   → F89 (SE,DE) full eigendecomposition IS the F86 K_b structural anchor (Tier-1-numerical)")
    print(f"#   → Closed form derivable from F89 AT-locked F_a/F_b structure (next step)")
    print(f"# If Endpoint/Interior split is present but values differ from empirical:")
    print(f"#   → Probe normalisation or kernel detail differs; reconcile with ResonanceScan")
    print(f"# If no class split or wildly off:")
    print(f"#   → Other coherence sectors (vac-SE, etc.) contribute; F89 (SE,DE) alone insufficient")


if __name__ == "__main__":
    main()
