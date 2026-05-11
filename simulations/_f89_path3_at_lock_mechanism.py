"""F89 path-3: WHY are F_a, F_b eigenvectors AT-rate-locked?

Hypothesis: the 4 closed-form eigenvalues from F_a (rate 2γ) and F_b
(rate 6γ) correspond to eigenvectors that are supported PURELY on the
overlap-only (resp no-overlap-only) basis pairs of (SE, DE). H_B mixes
within overlap (rate 2γ uniformly) and within no-overlap (rate 6γ
uniformly), but the cross-coupling overlap↔no-overlap is what generally
shifts rates. AT-locking happens iff the eigenvector is in the kernel of
the cross-coupling.

This script:
  1. Builds the 24-dim (SE, DE) basis with explicit overlap/no-overlap labels
  2. Constructs the 12-dim S_2-sym L_super sub-block at q=1.5 numerically
  3. Diagonalises and identifies the 4 AT-locked eigenvalues
  4. Reports the overlap vs no-overlap support of each eigenvector

If the hypothesis holds: F_a-eigenvectors have 100% support on overlap,
F_b-eigenvectors have 100% support on no-overlap.
"""

from __future__ import annotations

import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_path3_se_de_numerical(J: float, gamma: float):
    """Build 24-dim L_super on (SE, DE) for path-3 numerically.

    Returns (L, basis, overlap_mask) where:
      basis: list of (i, (j, k))
      overlap_mask: bool array, True if i in {j, k}
    """
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

    overlap_mask = np.array([i in jk for i, jk in basis])
    return L, basis, overlap_mask


def build_s2_sym_projector(basis):
    """Return P (24×12) such that columns are S_2-sym basis vectors."""
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


def main() -> None:
    J, gamma = 0.075, 0.05
    print(f"# F89 path-3 AT-lock mechanism, J={J}, γ={gamma}, q={J/gamma}\n")

    L, basis, overlap_mask = build_path3_se_de_numerical(J, gamma)
    n_overlap = overlap_mask.sum()
    n_no_overlap = len(basis) - n_overlap
    print(f"## (SE, DE) basis: {len(basis)} total = {n_overlap} overlap + {n_no_overlap} no-overlap")

    P = build_s2_sym_projector(basis)
    print(f"## S_2-sym projector P: {P.shape}")
    L_sym = P.conj().T @ L @ P
    print(f"## S_2-sym L sub-block: {L_sym.shape}")

    eigvals, eigvecs = np.linalg.eig(L_sym)
    rates = -eigvals.real / gamma
    freqs = eigvals.imag / J

    # Identify F_a, F_b roots (closed-form predictions)
    sqrt5 = np.sqrt(5)
    expected = {
        "F_a:λ=−2γ+iJ(−1+√5)": (-2 * gamma, J * (-1 + sqrt5)),
        "F_a:λ=−2γ+iJ(−1−√5)": (-2 * gamma, J * (-1 - sqrt5)),
        "F_b:λ=−6γ+iJ(−1+√5)": (-6 * gamma, J * (-1 + sqrt5)),
        "F_b:λ=−6γ+iJ(−1−√5)": (-6 * gamma, J * (-1 - sqrt5)),
    }

    print("\n## Matching closed-form eigenvalues to numerical L_sym diagonalisation")
    print("| Closed-form prediction | Matched λ (numerical) | Eigvec overlap-support | Eigvec no-overlap-support |")
    print("|---|---|---|---|")

    # Lift S_2-sym eigenvectors back to 24-dim basis to compute overlap/no-overlap support
    eigvecs_24 = P @ eigvecs  # (24, 12)

    for label, (re_pred, im_pred) in expected.items():
        # Find numerical eigenvalue closest to prediction
        target = re_pred + 1j * im_pred
        dists = np.abs(eigvals - target)
        k = np.argmin(dists)
        v = eigvecs_24[:, k]
        v_norm = np.linalg.norm(v)
        v = v / v_norm
        overlap_support = np.sum(np.abs(v[overlap_mask]) ** 2)
        no_overlap_support = np.sum(np.abs(v[~overlap_mask]) ** 2)
        print(f"| {label} | {eigvals[k]:.4f} | {overlap_support:.6f} | {no_overlap_support:.6f} |")

    print("\n## Hypothesis check:")
    print("# If F_a eigenvectors live PURELY in overlap subspace and F_b PURELY in no-overlap,")
    print("# then AT-locking is explained: H_B does not couple overlap↔no-overlap WITHIN those modes.")

    # Sanity: also report the other 8 (octic-derived) eigenvectors' overlap/no-overlap support
    print("\n## All 12 S_2-sym eigenvalues with overlap/no-overlap support split:")
    print("| λ | rate Γ/γ | freq ω/J | overlap support | no-overlap support |")
    print("|---|---|---|---|---|")
    sort_idx = np.argsort(rates)
    for k in sort_idx:
        v = eigvecs_24[:, k]
        v_norm = np.linalg.norm(v)
        v = v / v_norm
        ov = np.sum(np.abs(v[overlap_mask]) ** 2)
        no = np.sum(np.abs(v[~overlap_mask]) ** 2)
        print(f"| {eigvals[k]:+.4f} | {rates[k]:.4f} | {freqs[k]:+.4f} | {ov:.4f} | {no:.4f} |")

    # ===== 2nd verification: F_a, F_b live in SE-anti × DE-anti sub-block =====
    print("\n## Bloch decomposition: do F_a, F_b live in SE-anti × DE-anti?")
    print("# SE Bloch modes: ψ_n(j) = √(2/5)·sin(π·n·(j+1)/5), n=1..4; n odd → S_2-sym, n even → S_2-anti")
    print("# DE: 6 2-particle states; under S_2 mirror split into 4-sym + 2-anti")

    # SE Bloch matrix: U_SE[j, n] = √(2/5)·sin(πn(j+1)/5), j=0..3, n=1..4
    N_block = 4
    U_SE = np.array([
        [np.sqrt(2 / 5) * np.sin(np.pi * n * (j + 1) / 5) for n in range(1, 5)]
        for j in range(N_block)
    ])
    # SE single-particle E_n = 4J·cos(πn/5)
    E_SE = np.array([4 * J * np.cos(np.pi * n / 5) for n in range(1, 5)])
    print(f"# SE single-particle eigenvalues E_n: {E_SE}")
    print(f"#   E_1 = {E_SE[0]:.4f}J·units (S_2-sym, n=1)")
    print(f"#   E_2 = {E_SE[1]:.4f}J·units (S_2-anti, n=2)  ← matches F_a freq J(√5-1) ≈ {J*(np.sqrt(5)-1):.4f}")
    print(f"#   E_3 = {E_SE[2]:.4f}J·units (S_2-sym, n=3)")
    print(f"#   E_4 = {E_SE[3]:.4f}J·units (S_2-anti, n=4)  ← matches F_a freq -J(1+√5) ≈ {-J*(np.sqrt(5)+1):.4f}")

    # DE diagonalisation
    de_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    M_DE_phys = np.zeros((6, 6))
    for idx, (j, k) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j <= 3 and new_j != k:
                new_pair = tuple(sorted([new_j, k]))
                if abs(new_j - j) == 1 and new_pair in de_pairs:
                    M_DE_phys[de_pairs.index(new_pair), idx] += 2 * J
        for new_k in [k - 1, k + 1]:
            if 0 <= new_k <= 3 and new_k != j:
                new_pair = tuple(sorted([j, new_k]))
                if abs(new_k - k) == 1 and new_pair in de_pairs:
                    M_DE_phys[de_pairs.index(new_pair), idx] += 2 * J
    E_DE, U_DE = np.linalg.eigh(M_DE_phys)
    print(f"\n# DE 2-particle eigenvalues: {E_DE}")

    # Check S_2-symmetry of each SE and DE Bloch mode
    perm_se = [3, 2, 1, 0]
    perm_de = [de_pairs.index(tuple(sorted([3 - jk[0], 3 - jk[1]]))) for jk in de_pairs]

    se_sym_class = []
    for n in range(4):
        v = U_SE[:, n]
        v_mirrored = v[perm_se]
        if np.allclose(v_mirrored, v):
            se_sym_class.append("sym")
        elif np.allclose(v_mirrored, -v):
            se_sym_class.append("anti")
        else:
            se_sym_class.append("mixed")
    print(f"# SE Bloch S_2 class: n=1:{se_sym_class[0]}, n=2:{se_sym_class[1]}, n=3:{se_sym_class[2]}, n=4:{se_sym_class[3]}")

    de_sym_class = []
    for m in range(6):
        v = U_DE[:, m]
        v_mirrored = v[perm_de]
        if np.allclose(v_mirrored, v):
            de_sym_class.append("sym")
        elif np.allclose(v_mirrored, -v):
            de_sym_class.append("anti")
        else:
            de_sym_class.append("mixed")
    print(f"# DE Bloch S_2 class (E sorted): {de_sym_class} (energies {[round(e, 4) for e in E_DE]})")

    # Project F_a, F_b eigenvectors onto SE-anti × DE-anti sub-block
    se_anti_idx = [n for n in range(4) if se_sym_class[n] == "anti"]
    de_anti_idx = [m for m in range(6) if de_sym_class[m] == "anti"]
    print(f"\n# SE-anti modes: n = {[i+1 for i in se_anti_idx]}, total {len(se_anti_idx)}")
    print(f"# DE-anti modes: m at energies {[round(E_DE[i], 4) for i in de_anti_idx]}, total {len(de_anti_idx)}")

    # Build SE-anti × DE-anti basis as outer products in (SE, DE) basis
    # |ψ_n^SE⟩⟨ψ_m^DE| in basis-pair representation: M[i, jk] = U_SE[i, n] · U_DE[jk, m]
    # vec(M) in flat 24-dim basis
    print("\n## Support of F_a, F_b eigenvectors on SE-anti × DE-anti sub-block:")
    print("| Eigenvalue | SE-sym×DE-sym | SE-anti×DE-anti | SE-sym×DE-anti | SE-anti×DE-sym |")
    print("|---|---|---|---|---|")

    # Build the 4 sub-block projectors
    def build_subspace_projector(se_idx_list, de_idx_list):
        cols = []
        for n in se_idx_list:
            for m in de_idx_list:
                M = np.outer(U_SE[:, n], U_DE[:, m])  # (4, 6) matrix
                v = np.zeros(len(basis), dtype=complex)
                for idx, (i, jk) in enumerate(basis):
                    v[idx] = M[i, de_pairs.index(jk)]
                cols.append(v)
        return np.column_stack(cols) if cols else None

    se_sym_idx = [n for n in range(4) if se_sym_class[n] == "sym"]
    de_sym_idx = [m for m in range(6) if de_sym_class[m] == "sym"]

    P_ss = build_subspace_projector(se_sym_idx, de_sym_idx)
    P_aa = build_subspace_projector(se_anti_idx, de_anti_idx)
    P_sa = build_subspace_projector(se_sym_idx, de_anti_idx)
    P_as = build_subspace_projector(se_anti_idx, de_sym_idx)

    target_labels = ["F_a:E_2(=J(√5−1))", "F_a:E_4(=−J(1+√5))", "F_b:E_2", "F_b:E_4"]
    targets = [
        -2 * gamma + 1j * J * (np.sqrt(5) - 1),
        -2 * gamma - 1j * J * (1 + np.sqrt(5)),
        -6 * gamma + 1j * J * (np.sqrt(5) - 1),
        -6 * gamma - 1j * J * (1 + np.sqrt(5)),
    ]
    for label, target in zip(target_labels, targets):
        k = np.argmin(np.abs(eigvals - target))
        v = eigvecs_24[:, k]
        v = v / np.linalg.norm(v)
        # Project onto each sub-block (using least-squares to handle non-orthonormality)
        def support(P):
            if P is None or P.shape[1] == 0:
                return 0.0
            coefs, _, _, _ = np.linalg.lstsq(P, v, rcond=None)
            v_proj = P @ coefs
            return float(np.real(np.vdot(v_proj, v_proj)))
        s_ss = support(P_ss)
        s_aa = support(P_aa)
        s_sa = support(P_sa)
        s_as = support(P_as)
        print(f"| {label}: λ={eigvals[k]:.4f} | {s_ss:.4f} | {s_aa:.4f} | {s_sa:.4f} | {s_as:.4f} |")


if __name__ == "__main__":
    main()
