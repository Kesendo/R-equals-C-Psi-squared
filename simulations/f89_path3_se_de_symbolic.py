"""F89 path-3 (SE, DE) sub-block: symbolic factorisation attempt.

Path-3 (N_block=4) has 4 SE states × 6 DE states = 24-dim (SE, DE) sub-block.
Under S_2 chain-mirror (0↔3, 1↔2), splits into 12-dim S_2-sym + 12-dim S_2-anti.

Goal: build the 12×12 S_2-sym L_super sub-block symbolically (γ=1, q=J/γ as
the only symbol), compute char(λ) via Faddeev-Leverrier, and check whether
it factors into lower-degree pieces (analog of path-2's linear · linear · cubic).

Per F89c at γ=J=1: path-3 has rates {0, 2, 4, 6, 8}γ corresponding to
n_diff ∈ {0, 1, 2, 3, 4}. For (SE, DE) coherences, n_diff ∈ {1, 3}
(overlap, no-overlap). Expected pure-AT linear factors at λ = -2γ and
λ = -6γ. The H_B-mixed factors carry q-dependence.
"""

from __future__ import annotations

import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_se_de_path3_symbolic():
    """Build 24-dim L_super on (SE, DE) sector for path-3 symbolically.

    γ=1; q := J/γ as the only symbolic variable.
    """
    q = sp.symbols("q", positive=True, real=True)
    g = sp.Integer(1)

    de_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    basis = []
    for i in range(4):
        for jk in de_pairs:
            basis.append((i, jk))
    n_basis = len(basis)  # 24

    # H_B^SE for path-3 chain (4 sites): tridiagonal 4×4, off-diag 2q
    M_SE = sp.Matrix(4, 4, lambda a, b: 2 * q if abs(a - b) == 1 else 0)

    # H_B^DE for path-3: 6×6 in DE basis order
    M_DE = sp.zeros(6, 6)
    for idx, (j, k) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if new_j < 0 or new_j > 3 or new_j == k:
                continue
            new_pair = tuple(sorted([new_j, k]))
            if abs(new_j - j) == 1 and new_pair in de_pairs:
                M_DE[de_pairs.index(new_pair), idx] += 2 * q
        for new_k in [k - 1, k + 1]:
            if new_k < 0 or new_k > 3 or new_k == j:
                continue
            new_pair = tuple(sorted([j, new_k]))
            if abs(new_k - k) == 1 and new_pair in de_pairs:
                M_DE[de_pairs.index(new_pair), idx] += 2 * q

    L = sp.zeros(n_basis, n_basis)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(4):
            if M_SE[i2, i] != 0:
                idx2 = basis.index((i2, jk))
                L[idx2, idx] += -sp.I * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(6):
            if M_DE[jk_idx, jk2_idx] != 0:
                jk2 = de_pairs[jk2_idx]
                idx2 = basis.index((i, jk2))
                L[idx2, idx] += sp.I * M_DE[jk_idx, jk2_idx]
        L[idx, idx] += -2 * g if i in jk else -6 * g

    return L, basis, q, de_pairs


def build_s2_sym_basis(basis, de_pairs):
    """12-dim S_2-sym + 12-dim S_2-anti basis for path-3 (mirror sites 0↔3, 1↔2).

    Returns sym_basis only; for path-3 there are no fixed (i, jk) pairs since
    R fixes no SE state, so all 24 basis pairs split into 12 size-2 orbits.
    """
    perm_se = {0: 3, 1: 2, 2: 1, 3: 0}
    perm_de = {jk: tuple(sorted([perm_se[jk[0]], perm_se[jk[1]]])) for jk in de_pairs}

    sym = []
    handled = set()
    for idx, (i, jk) in enumerate(basis):
        if idx in handled:
            continue
        idx2 = basis.index((perm_se[i], perm_de[jk]))
        v = sp.zeros(len(basis), 1)
        if idx == idx2:
            v[idx, 0] = 1
        else:
            v[idx, 0] = sp.Rational(1, 1) / sp.sqrt(2)
            v[idx2, 0] = sp.Rational(1, 1) / sp.sqrt(2)
            handled.add(idx2)
        sym.append(v)
        handled.add(idx)
    return sym


def main() -> None:
    print("# F89 path-3 (SE, DE) symbolic factorisation attempt", flush=True)
    print("# γ=1, q := J/γ as only symbol\n", flush=True)

    L, basis, q, de_pairs = build_se_de_path3_symbolic()
    print(f"## (SE, DE) sub-block dim: {L.shape[0]}", flush=True)

    sym = build_s2_sym_basis(basis, de_pairs)
    print(f"## S_2-symmetric subspace dim: {len(sym)} (expected 12)", flush=True)

    P = sp.Matrix.hstack(*sym)  # 24 × 12
    print("\n## Projecting L_super onto 12-dim S_2-sym basis...", flush=True)
    L_sym = P.T * L * P  # NO simplify — kills perf; sympy handles rationals fine
    print(f"# L_sym shape: {L_sym.shape}", flush=True)

    print("\n## Computing characteristic polynomial via Matrix.charpoly...", flush=True)
    lam = sp.symbols("lambda")
    cp = L_sym.charpoly(lam)
    char_poly = sp.expand(cp.as_expr())
    print(f"# Polynomial degree in λ: {sp.degree(char_poly, lam)}", flush=True)

    print("\n## Factoring characteristic polynomial...", flush=True)
    factored = sp.factor(char_poly, lam)
    print()
    print(factored)
    print()

    print("## Factor structure:")
    if factored.func == sp.Mul:
        for f in factored.args:
            if f.is_number:
                print(f"  const: {f}")
            else:
                p = sp.Poly(f, lam) if not f.is_Pow else sp.Poly(f.base, lam)
                deg = p.degree()
                if f.is_Pow:
                    print(f"  ({f.base})^{f.exp}  →  degree-{deg} factor with multiplicity {f.exp}")
                else:
                    print(f"  degree-{deg}: {f}")
    else:
        print(f"  Single irreducible factor of degree {sp.degree(factored, lam)}: {factored}")

    print("\n## Numerical eigenvalues at q = 1.5 (in units of γ):")
    char_num = char_poly.subs(q, sp.Rational(3, 2))
    roots = sp.nroots(sp.Poly(char_num, lam), n=10)
    print(f"# {len(roots)} eigenvalues:")
    for r in sorted(roots, key=lambda x: (float(sp.re(x)), float(sp.im(x)))):
        rc = complex(r)
        rate = -rc.real
        freq_over_J = rc.imag / 1.5
        print(f"#   λ/γ = {rc.real:+.6f}{rc.imag:+.6f}j  →  rate Γ/γ = {rate:.4f}, freq ω/J = {freq_over_J:+.4f}")


if __name__ == "__main__":
    main()
