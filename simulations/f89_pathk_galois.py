"""F89 path-4/5/6 H_B-mixed Galois groups: F_18 / F_32 / F_53 = S_18 / S_32 / S_53.

Closes the Tier-2 "conjecturally Galois-non-solvable for degree ≥ 5" flag: each
path-k H_B-mixed Liouvillian factor is IRREDUCIBLE over Q(i)(q) with the FULL
symmetric group S_n, hence non-solvable (the eigenvalues λ_k(q) admit no radical
expression in q). Generalises the path-3 octic = S_8 result (f89_path3_octic_galois.py).

Method (gate-first, same chain as path-3, generalised Jordan):
  - build the path-k (SE,DE) S_2-sym block at integer q0 EXACTLY over Q(i) using the
    integer mirror basis e_n+e_m, so M = (BᵀB)⁻¹ BᵀLB has /2 denominators, no √2
    (its charpoly is a true factor of L's charpoly — verified: it divides exactly);
  - charpoly, factor over Q(i); bucket factors by decay rate: AT-locked roots have
    Re(λ) ∈ {−2,−6} (rate 2γ/6γ); the single remaining factor is the H_B-mixed F_d;
  - Frobenius cycle types of F_d over split primes 𝔭|p (p≡1 mod4, i↦√−1; Dedekind);
  - transitive (F_d irreducible) + a cycle type with a part that is a PRIME p with
    n/2 < p ≤ n−3  (⟹ a p-cycle ⟹ primitive [p>n/2] ⟹ ⊇A_n [Jordan, p≤n−3])
    + an ODD cycle type (n−#parts odd ⟹ ⊄A_n)  ⟹  S_n.
  Specialization can only SHRINK the group, so G_{q0}=S_n forces the generic G=S_n.

RESULTS (q0=2, corroborated at q0=3 with different witness primes; adversarially
reviewed, PROVEN — no hole):
  path-4 (N_block=5): sym dim 26 = 8 AT-locked + F_18, witness prime 13 ⟹ Gal = S_18
  path-5 (N_block=6): sym dim 45 = 13 AT-locked + F_32, witness prime 23 ⟹ Gal = S_32
  path-6 (N_block=7): sym dim 75 = 22 AT-locked + F_53, witness prime 37 ⟹ Gal = S_53

S_n is the GENERIC group of an irreducible degree-n polynomial (van der Waerden;
Bhargava, Annals 2025): the content is negative — integrability spends itself on the
AT-locked F_a/F_b factorisation (single-particle frequencies in radicals), leaving the
H_B-mixed residue maximally generic; the closed-form program ends at the AT-protected half.
"""

from __future__ import annotations

import os
import sys
import time
from itertools import combinations

import numpy as np
import sympy as sp
from sympy.ntheory import isprime

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from f89_path3_octic_galois import frobenius_cycle_types, is_irreducible_over_qi  # noqa: E402

LAM = sp.symbols("lambda")


def build_pathk_sym_over_qi(q0: int, n_block: int):
    """L|_sym (S_2-symmetric (SE,DE) sub-block) at integer q0, γ=1, as a sympy Matrix
    over Q(i). Integer mirror basis (e_n + e_m): entries are Gaussian rationals with
    denominator 1 or 2 — no √2, so the charpoly stays over Q(i). Eigenvalue/charpoly
    cross-checked against the orthonormal-basis builder in
    f89_path4_path5_at_lock_scan.build_pathk_se_de_sym."""
    q = sp.Integer(q0)
    de_pairs = list(combinations(range(n_block), 2))
    basis = [(i, jk) for i in range(n_block) for jk in de_pairs]
    nb = len(basis)
    idx_of = {b: t for t, b in enumerate(basis)}

    L = sp.zeros(nb, nb)
    for col, (i, jk) in enumerate(basis):
        for i2 in (i - 1, i + 1):                                  # SE hop 2q (ket): −i·M_SE
            if 0 <= i2 < n_block:
                L[idx_of[(i2, jk)], col] += -sp.I * 2 * q
        j, k = jk                                                  # DE hop 2q (bra): +i·M_DE
        for new_j in (j - 1, j + 1):
            if 0 <= new_j < n_block and new_j != k:
                npair = tuple(sorted((new_j, k)))
                if npair in de_pairs:
                    L[idx_of[(i, npair)], col] += sp.I * 2 * q
        for new_k in (k - 1, k + 1):
            if 0 <= new_k < n_block and new_k != j:
                npair = tuple(sorted((j, new_k)))
                if npair in de_pairs:
                    L[idx_of[(i, npair)], col] += sp.I * 2 * q
        L[col, col] += (-2 if i in jk else -6)                     # γ=1 diagonal

    def perm_se(i): return n_block - 1 - i
    def perm_de(p): return tuple(sorted((perm_se(p[0]), perm_se(p[1]))))
    cols, handled = [], set()
    for t, (i, jk) in enumerate(basis):
        if t in handled:
            continue
        t2 = idx_of[(perm_se(i), perm_de(jk))]
        v = [0] * nb
        v[t] = 1
        if t2 != t:
            v[t2] = 1
            handled.add(t2)
        handled.add(t)
        cols.append(v)
    B = sp.Matrix.hstack(*[sp.Matrix(nb, 1, c) for c in cols])      # nb × dim_sym
    BtB = B.T * B                                                  # diagonal (1 or 2)
    return BtB.inv() * (B.T * L * B)                               # L|_sym over Q(i)


def rate_bucket_factors(charpoly_expr):
    """Factor charpoly over Q(i); split into (AT-locked factors, H_B-mixed factors).
    AT-locked = every root has Re(λ) ∈ {−2,−6} (rate 2γ/6γ). Roots via numpy.roots
    (double precision distinguishes Re=−2/−6 with ~5 orders of margin; nroots fails to
    converge at high degree)."""
    factored = sp.factor(charpoly_expr, LAM, extension=[sp.I])
    args = factored.args if factored.func == sp.Mul else [factored]
    at, hb = [], []
    for a in args:
        base = a.base if a.func == sp.Pow else a
        if not base.has(LAM):
            continue
        mult = int(a.exp) if a.func == sp.Pow else 1
        roots = np.roots([complex(c) for c in sp.Poly(base, LAM).all_coeffs()])
        is_at = all(abs(r.real + 2) < 1e-7 or abs(r.real + 6) < 1e-7 for r in roots)
        (at if is_at else hb).extend([base] * mult)
    return at, hb


def jordan_verdict(cycle_types, n):
    """Transitivity assumed (F_d irreducible). A prime p∈(n/2,n−3] cycle ⟹ ⊇A_n
    (primitive + Jordan); an odd cycle type ⟹ ⊄A_n. Both ⟹ S_n."""
    cts = list(cycle_types)
    big_prime = next((p for ct in cts for p in ct if isprime(p) and n / 2 < p <= n - 3), None)
    odd = any((n - len(ct)) % 2 == 1 for ct in cts)
    if big_prime and odd:
        return f"S_{n}  (non-solvable)", big_prime
    if big_prime:
        return f"⊇A_{n} (A_{n} or S_{n}; no odd element seen yet)", big_prime
    return (f"INDETERMINATE — no prime p∈({int(n/2)+1}..{n-3}) cycle found "
            f"(more primes, or a proper/solvable group)"), None


def analyse_path(k: int, q0: int, prime_hi: int):
    n_block = k + 1
    print("=" * 78)
    print(f"PATH-{k}  (N_block={n_block})  at q0={q0}")
    print("=" * 78)
    t0 = time.time()
    M = build_pathk_sym_over_qi(q0, n_block)
    cp = sp.expand(M.charpoly(LAM).as_expr())
    at, hb = rate_bucket_factors(cp)
    print(f"  L|_sym {M.shape[0]}×{M.shape[1]} over Q(i); charpoly degree {sp.degree(cp, LAM)} "
          f"(built+factored in {time.time()-t0:.1f}s)")
    print(f"  {len(at)} AT-locked factor(s) deg {sorted((sp.Poly(f,LAM).degree() for f in at), reverse=True)} "
          f"(sum {sum(sp.Poly(f,LAM).degree() for f in at)}), "
          f"{len(hb)} H_B-mixed factor(s) deg {sorted((sp.Poly(f,LAM).degree() for f in hb), reverse=True)}")

    primes = list(sp.primerange(5, prime_hi))
    for i, fd in enumerate(sorted(hb, key=lambda f: -sp.Poly(f, LAM).degree())):
        n = sp.Poly(fd, LAM).degree()
        irr = is_irreducible_over_qi(fd, LAM)
        types = frobenius_cycle_types(fd, LAM, primes)
        verdict, bp = jordan_verdict(types.keys(), n)
        print(f"\n  H_B-mixed factor #{i+1}: degree {n}, irreducible/Q(i)={irr}, "
              f"{len(types)} distinct cycle types over {len([p for p in primes if p%4==1])} split primes")
        if bp:
            ex = next(ct for ct in types if bp in ct)
            print(f"     prime-{bp} cycle present (e.g. {ex}); odd element present: "
                  f"{any((n-len(ct))%2==1 for ct in types)}")
        print(f"     ===> Gal(F_{n}) over Q(i) = {verdict}")


GATES = {
    "x^18 - x - 1 (Selmer, S_18)": (sp.Symbol("x")**18 - sp.Symbol("x") - 1, 18, "S_18"),
    "x^18 - 2     (solvable)":     (sp.Symbol("x")**18 - 2, 18, "PROPER"),
}


def run_gates():
    print("=" * 78)
    print("GATE  generalised Jordan verdict on known degree-18 answers")
    print("=" * 78)
    x = sp.Symbol("x")
    primes = list(sp.primerange(5, 600))
    ok = True
    for name, (poly, n, expect) in GATES.items():
        types = frobenius_cycle_types(poly, x, primes)
        verdict, _ = jordan_verdict(types.keys(), n)
        passed = (expect == "S_18" and verdict.startswith("S_")) or \
                 (expect == "PROPER" and not verdict.startswith("S_") and "⊇A_" not in verdict)
        ok = ok and passed
        print(f"  {name}: {len(types)} cycle types; {verdict[:48]}  [{'PASS' if passed else 'FAIL'}]")
    print(f"  GATE: {'PASS' if ok else 'FAIL'}\n")


def emit_fd_literal(k: int, q0: int = 2):
    """Emit the H_B-mixed factor F_d(λ, q0) as exact integer (re, im) coefficient lists,
    lowest-first (the OcticCoefficientsAtQ2 convention) — the committed oracle the C# live
    Option-D witness is validated against (assert: live F_d == this literal). Run:
        python simulations/f89_pathk_galois.py emit 4
    """
    n_block = k + 1
    M = build_pathk_sym_over_qi(q0, n_block)
    cp = sp.expand(M.charpoly(LAM).as_expr())
    _, hb = rate_bucket_factors(cp)
    assert len(hb) == 1, f"path-{k}: expected one H_B-mixed factor, got {len(hb)}"
    fd = sp.Poly(hb[0], LAM)
    coeffs = list(reversed(fd.all_coeffs()))  # lowest-first
    re, im = [], []
    for c in coeffs:
        cr, ci = sp.nsimplify(sp.re(c)), sp.nsimplify(sp.im(c))
        assert cr.is_integer and ci.is_integer, f"non-integer F_d coefficient {c}"
        re.append(int(cr))
        im.append(int(ci))
    print(f"// path-{k} F_d(lambda, q={q0}); degree {fd.degree()}; lowest-first; monic; over Z[i]")
    print(f"re = {re}")
    print(f"im = {im}")
    return re, im


def emit_at_factors(k: int, q0: int = 2):
    """Print each AT-locked factor's degree and roots (Re, Im) — to study the F_a (rate 2γ,
    Re=−2) single-particle Bloch modes and the F_b (rate 6γ, Re=−6) 2-particle DE-Slater modes,
    the structure the live-D AT reconstruction must reproduce over Z[i]."""
    n_block = k + 1
    M = build_pathk_sym_over_qi(q0, n_block)
    cp = sp.expand(M.charpoly(LAM).as_expr())
    at, hb = rate_bucket_factors(cp)
    print(f"path-{k} (n_block={n_block}, q0={q0}): {len(at)} AT factor(s), {len(hb)} H_B-mixed factor(s)")
    for f in sorted(at, key=lambda x: sp.Poly(x, LAM).degree()):
        poly = sp.Poly(f, LAM)
        roots = np.roots([complex(c) for c in poly.all_coeffs()])
        coeffs = list(reversed(poly.all_coeffs()))
        intcoef = all(sp.nsimplify(sp.re(c)).is_integer and sp.nsimplify(sp.im(c)).is_integer for c in coeffs)
        rs = ", ".join(f"({r.real:+.2f},{r.imag:+.2f})" for r in sorted(roots, key=lambda z: (round(z.real, 2), z.imag)))
        print(f"  AT factor degree {poly.degree()} (integer coeffs: {intcoef}): {rs}")


def gen_csharp_literals():
    """Generate the committed C# F_d literal class F89PathKFdLiterals.cs (no hand-transcription;
    reproducible source-of-truth the semi-live witnesses verify against the live block). Run:
        python simulations/f89_pathk_galois.py gen-cs
    """
    out = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..",
        "compute", "RCPsiSquared.Core", "F89PathK", "F89PathKFdLiterals.cs"))
    lines = [
        "namespace RCPsiSquared.Core.F89PathK;",
        "",
        "// GENERATED by `python simulations/f89_pathk_galois.py gen-cs` — do not edit by hand.",
        "// F_d(lambda, 2) literals (lowest-first, monic, over Z[i]) for the path-k H_B-mixed factors;",
        "// the semi-live witnesses verify each divides the LIVE block characteristic polynomial.",
        "internal static class F89PathKFdLiterals",
        "{",
    ]
    for k in (4, 5, 6):
        re, im = emit_fd_literal(k)
        re_s = ", ".join('"%d"' % c for c in re)
        im_s = ", ".join('"%d"' % c for c in im)
        lines.append("    internal static readonly string[] Path%dRe = { %s };" % (k, re_s))
        lines.append("    internal static readonly string[] Path%dIm = { %s };" % (k, im_s))
    lines.append("}")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")


def reconstruct_at_invariant_subspace(k: int, q0: int = 2):
    """Full-D AT reconstruction (no F_d import): the AT-locked factor = product over rate sectors of
    charpoly(M | W), where W is the largest M-invariant subspace inside the sector (Re(λ) = −2 or −6).
    M = D + iK (D real diagonal rates, K real-symmetric hopping); W = nullspace of [P_Uc K^m]. This
    is the rate-confined-invariant-subspace construction the C# witness ports."""
    M = build_pathk_sym_over_qi(q0, k + 1)
    n = M.shape[0]
    rate = [sp.re(M[i, i]) for i in range(n)]
    K = M.applyfunc(sp.im)
    at = sp.Integer(1)
    for r0 in (-2, -6):
        mask = [rate[i] == r0 for i in range(n)]
        if not any(mask):
            continue
        p_u = sp.diag(*[0 if mask[i] else 1 for i in range(n)])
        blocks, mm = [], p_u
        for _ in range(n):
            blocks.append(mm)
            mm = mm * K
        ns = sp.Matrix.vstack(*blocks).nullspace()       # largest M-invariant subspace in the sector
        if not ns:
            continue
        W = sp.Matrix.hstack(*ns)
        B = (W.T * W).inv() * (W.T * (M * W))            # M | W
        at *= B.charpoly(LAM).as_expr()
    return sp.Poly(sp.expand(at), LAM).monic()


def gate_full_d():
    """Gate: the invariant-subspace AT reconstruction divides the full charpoly with quotient F_d,
    for path-4/5/6 (full Option D, F_d never imported). Run: f89_pathk_galois.py full-d."""
    print("=" * 78)
    print("GATE full-D: AT (rate-confined invariant subspace) divides C with quotient F_d, path-4/5/6")
    print("=" * 78)
    ok = True
    for k in (4, 5, 6):
        M = build_pathk_sym_over_qi(2, k + 1)
        cp = sp.expand(M.charpoly(LAM).as_expr())
        _, hb = rate_bucket_factors(cp)
        fd = sp.Poly(hb[0], LAM).monic()
        at = reconstruct_at_invariant_subspace(k)
        quotient, rem = sp.div(sp.Poly(cp, LAM).monic().as_expr(), at.as_expr(), LAM)
        match = sp.expand(rem) == 0 and \
            sp.expand(sp.Poly(quotient, LAM).monic().as_expr() - fd.as_expr()) == 0
        ok = ok and match
        print(f"  path-{k}: AT deg {at.degree()}, C / AT == F_{fd.degree()}: {match}")
    print(f"  GATE full-D: {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "full-d":
        gate_full_d()
    elif len(sys.argv) >= 2 and sys.argv[1] == "gen-cs":
        gen_csharp_literals()
    elif len(sys.argv) >= 3 and sys.argv[1] == "emit":
        emit_fd_literal(int(sys.argv[2]), q0=int(sys.argv[3]) if len(sys.argv) >= 4 else 2)
    elif len(sys.argv) >= 3 and sys.argv[1] == "emit_at":
        emit_at_factors(int(sys.argv[2]))
    else:
        run_gates()
        analyse_path(4, q0=2, prime_hi=1500)
        analyse_path(5, q0=2, prime_hi=2500)
        analyse_path(6, q0=2, prime_hi=4000)
