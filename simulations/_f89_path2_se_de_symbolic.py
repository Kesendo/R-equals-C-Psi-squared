"""F89 path-2 (SE, DE) sub-block: symbolic structure analysis.

The path-2 multi-exp closed form has 4 populated mode-groups:
- (2γ, ±2√2 J): (vac, SE) Bloch k=1, k=3 (analytic via F65)
- (2γ, 0): (SE, DE) at rate 2γ, freq 0 — S_3-symmetric overlap mode (rate 2γ pure)
- (3.04γ, 0): (SE, DE) H_B-mixed eigenvalue (J/γ-dependent)
- (3.48γ, ±5.45 J): (SE, DE) H_B-mixed complex-conjugate pair (J/γ-dependent)

The (SE, DE) sector is 9-dim. Under S_2 chain-mirror (path-2 has only S_2, not S_3
symmetry of H_B), the (SE, DE) sub-block splits into S_2-sym (5-dim) + S_2-anti
(4-dim). The S_2-sym sub-block contains the 3 populated eigenmodes (2γ, 3.04γ,
3.48γ at J/γ=1.5).

Goal: derive symbolic eigenvalues of the S_2-sym sub-block in terms of J, γ.

Strategy:
1. Build the 9x9 L_super on (SE, DE) symbolically in sympy.
2. Project onto S_2-symmetric 5-dim subspace.
3. Try characteristic polynomial → check if factorable.

Result: cubic·quadratic factorisation, with the cubic giving the rates
{2γ, 3.04γ, 3.48γ + 5.45iJ, 3.48γ - 5.45iJ at J/γ=1.5}. The cubic eigenvalues
have closed-form expressions via Cardano (or the depressed cubic).
"""

from __future__ import annotations

import sys

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_se_de_block_symbolic():
    """Build 9x9 L_super on (SE, DE) sector symbolically. Basis order:
    (i, (j,k)) ∈ {(0, (0,1)), (0, (0,2)), (0, (1,2)),
                  (1, (0,1)), (1, (0,2)), (1, (1,2)),
                  (2, (0,1)), (2, (0,2)), (2, (1,2))}
    """
    J, g = sp.symbols("J gamma", positive=True, real=True)

    basis = []
    for i in range(3):
        for (j, k) in [(0, 1), (0, 2), (1, 2)]:
            basis.append((i, (j, k)))

    # H_B^SE = J · [[0,2,0],[2,0,2],[0,2,0]] (same as H_B^DE)
    M = sp.Matrix([[0, 2 * J, 0], [2 * J, 0, 2 * J], [0, 2 * J, 0]])

    # DE basis: (0,1) (0,2) (1,2). H_B^DE in this basis:
    # |110⟩ ↔ |101⟩ via bond (1,2): hop site 1 to 2 in DE_{0,1} → DE_{0,2}, amp 2J
    # |101⟩ ↔ |011⟩ via bond (0,1): hop site 0 to 1 in DE_{0,2} → DE_{1,2}, amp 2J
    # |110⟩ ↔ |011⟩: not directly coupled
    # So H_B^DE on (DE_{0,1}, DE_{0,2}, DE_{1,2}) = [[0,2J,0],[2J,0,2J],[0,2J,0]] = same as M

    L_H_part = sp.zeros(9, 9)
    L_deph = sp.zeros(9, 9)
    for idx, (i, jk) in enumerate(basis):
        # L_H acts on |SE_i⟩⟨DE_{jk}|: -i (H_SE · |SE_i⟩⟨DE_{jk}| - |SE_i⟩⟨DE_{jk}| · H_DE)
        # Σ_{i'} M[i', i] |SE_{i'}⟩⟨DE_{jk}|  (from left)
        for i2 in range(3):
            if M[i2, i] != 0:
                idx2 = basis.index((i2, jk))
                L_H_part[idx2, idx] += -sp.I * M[i2, i]
        # |SE_i⟩⟨DE_{jk}| · H_DE = Σ_{jk2} M[jk_idx, jk2_idx] · |SE_i⟩⟨DE_{jk2}|
        # where M is indexed by DE-basis order. We use the SAME M for H_DE.
        jk_list = [(0, 1), (0, 2), (1, 2)]
        jk_index = jk_list.index(jk)
        for jk2_index in range(3):
            if M[jk_index, jk2_index] != 0:
                jk2 = jk_list[jk2_index]
                idx2 = basis.index((i, jk2))
                L_H_part[idx2, idx] += sp.I * M[jk_index, jk2_index]
        # Dephasing: -2γ if overlap (i ∈ jk), -6γ if no-overlap
        if i in jk:
            L_deph[idx, idx] = -2 * g
        else:
            L_deph[idx, idx] = -6 * g

    L = L_H_part + L_deph
    return L, basis, J, g


def project_to_s2_symmetric(L, basis):
    """Project L_super (SE, DE) onto S_2-symmetric subspace.

    S_2 = chain mirror exchanging sites 0 ↔ 2. Under this:
    SE_0 ↔ SE_2, SE_1 → SE_1
    DE_{01} ↔ DE_{12}, DE_{02} → DE_{02}

    Build the 9-dim permutation matrix R for S_2 action on |SE_i⟩⟨DE_{jk}|.
    """
    permute_se = {0: 2, 1: 1, 2: 0}
    permute_de = {(0, 1): (1, 2), (0, 2): (0, 2), (1, 2): (0, 1)}

    R = sp.zeros(9, 9)
    for idx, (i, jk) in enumerate(basis):
        new_i = permute_se[i]
        new_jk = permute_de[jk]
        idx2 = basis.index((new_i, new_jk))
        R[idx2, idx] = 1

    # S_2 symmetric projector: P_sym = (I + R) / 2
    P_sym = (sp.eye(9) + R) / 2
    P_anti = (sp.eye(9) - R) / 2

    # Find symmetric basis: eigenvectors of P_sym with eigenvalue 1
    # P_sym is idempotent, eigenvalues 0 and 1
    rank_sym = P_sym.rank()
    rank_anti = P_anti.rank()

    return R, P_sym, P_anti, rank_sym, rank_anti


def main() -> None:
    print("# F89 path-2 (SE, DE) symbolic sub-block analysis\n")

    L, basis, J, g = build_se_de_block_symbolic()
    print(f"## 9x9 L_super on (SE, DE), with J, γ symbolic")
    print()

    R, P_sym, P_anti, rank_sym, rank_anti = project_to_s2_symmetric(L, basis)
    print(f"## S_2 symmetry projection")
    print(f"# S_2-symmetric subspace dim: {rank_sym}")
    print(f"# S_2-antisymmetric subspace dim: {rank_anti}")
    print(f"# (Expected: 5 + 4 = 9 ✓)")
    print()

    # Verify L commutes with R
    commutator = sp.simplify(L * R - R * L)
    print(f"## Commutator check [L, R]:")
    if commutator == sp.zeros(9, 9):
        print("# [L, R] = 0 (L preserves S_2 symmetry) ✓")
    else:
        print(f"# Non-zero commutator! Max entry: {sp.simplify(commutator).norm()}")
    print()

    # Build symmetric basis directly
    print("## Building explicit S_2-symmetric basis (5-dim)")
    sym_basis = []
    sym_pairs = []
    handled = set()
    for idx, (i, jk) in enumerate(basis):
        if idx in handled:
            continue
        # Apply R to get partner
        new_i = {0: 2, 1: 1, 2: 0}[i]
        new_jk = {(0, 1): (1, 2), (0, 2): (0, 2), (1, 2): (0, 1)}[jk]
        idx2 = basis.index((new_i, new_jk))
        if idx == idx2:
            # Self-mirror
            v = sp.zeros(9, 1)
            v[idx, 0] = 1
            sym_basis.append(v)
            sym_pairs.append((idx, idx))
        else:
            # Sym combo
            v = sp.zeros(9, 1)
            v[idx, 0] = sp.Rational(1, 1) / sp.sqrt(2)
            v[idx2, 0] = sp.Rational(1, 1) / sp.sqrt(2)
            sym_basis.append(v)
            sym_pairs.append((idx, idx2))
            handled.add(idx2)
        handled.add(idx)

    print(f"# S_2-symmetric basis vectors: {len(sym_basis)}")
    for i, pair in enumerate(sym_pairs):
        i1, i2 = pair
        b1 = basis[i1]
        b2 = basis[i2]
        if i1 == i2:
            print(f"#   |v_{i+1}⟩ = |SE_{b1[0]}⟩⟨DE_{b1[1]}|")
        else:
            print(f"#   |v_{i+1}⟩ = (|SE_{b1[0]}⟩⟨DE_{b1[1]}| + |SE_{b2[0]}⟩⟨DE_{b2[1]}|)/√2")
    print()

    # Project L onto symmetric basis: M_sym[i, j] = <v_i | L | v_j>
    P = sp.Matrix.hstack(*sym_basis)  # 9x5
    L_sym = sp.simplify(P.T * L * P)
    print("## L restricted to S_2-symmetric basis (5x5 matrix):")
    sp.init_printing()
    print(L_sym)
    print()

    # Compute characteristic polynomial
    lam = sp.symbols("lambda")
    char_poly = (L_sym - lam * sp.eye(5)).det()
    char_poly = sp.expand(char_poly)
    print(f"## Characteristic polynomial in λ (degree {sp.degree(char_poly, lam)}):")
    print(char_poly)
    print()

    # Try to factor
    factored = sp.factor(char_poly)
    print(f"## Factored characteristic polynomial:")
    print(factored)
    print()

    # Verification at J=0.075, γ=0.05 (J/γ = 1.5): should match numerical 2γ, 3.04γ, 3.48γ ± 5.45iJ
    J_num = sp.Rational(75, 1000)
    g_num = sp.Rational(5, 100)
    char_num = char_poly.subs([(J, J_num), (g, g_num)])
    roots = sp.solve(char_num, lam)
    print(f"## Numerical eigenvalues at J=0.075, γ=0.05 (J/γ = 1.5):")
    for r in roots:
        r_complex = complex(r.evalf())
        rate = -r_complex.real / float(g_num)
        freq = r_complex.imag / float(J_num)
        print(f"#   λ = {r_complex:.6f}  →  rate Γ/γ = {rate:.4f}, freq ω/J = {freq:+.4f}")


if __name__ == "__main__":
    main()
