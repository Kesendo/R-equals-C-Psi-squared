"""One 2 x 2 in many settings: the rows of docs/THE_DOUBLE_ROOT.md's section of that name.

Every row is a 2 x 2, written as its own page writes it; where the page writes a characteristic polynomial or a
second-order equation, the companion matrix is this script's construction, and for it the "non-zero and nilpotent"
reading cannot fail (a companion with a double root is always one Jordan block): there the checks that can fail
are d and the roots, and the Jordan verdict is the physical block's. For each: M = (tr M / 2)·1 + M0 with M0
traceless. M0² = d·1 holds for every traceless 2 x 2 (Cayley-Hamilton), so the test that can fail is that d is the
row's stated value; then the page's eigenvalues are roots of the characteristic polynomial, and at the page's
point M0 is non-zero and nilpotent: one 2 x 2 Jordan block. Plus the facts the section states about single rows
(the edge block's part, the seed toy's self-orthogonal kernel, the F155 dimer's Liouvillian, the writing-dependence
of self-orthogonality). Everything exact sympy; parameters symbols.

Run: python simulations/two_by_two_ep_family.py
"""
import sys

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CLAIMS = []


def claim(label, ok):
    CLAIMS.append(bool(ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {label}")


I = sp.I
lam = sp.symbols("lambda")
J, g, gam, gam0, gam1, B, s, cN, h, x, eps, k, Jg, G, kap, cc = sp.symbols(
    "J g gamma gamma0 gamma1 B s kappa_N h x epsilon k Jg G kappa c", positive=True)
u = sp.symbols("u")
zs = sp.symbols("s_seed")


def zero(M):
    return M.applyfunc(sp.simplify) == sp.zeros(*M.shape)


def companion(a1, a0):
    """the companion matrix of lambda² + a1·lambda + a0"""
    return sp.Matrix([[0, 1], [-a0, -a1]])


def row(name, M, d_expected, eigen, point):
    M = sp.Matrix(M)
    M0 = (M - (M.trace() / 2) * sp.eye(2)).applyfunc(sp.simplify)
    d = sp.simplify((M0 * M0)[0, 0])
    chi = sp.expand((M - lam * sp.eye(2)).det())
    roots_ok = all(sp.simplify(chi.subs(lam, e)) == 0 for e in eigen)
    at = M0.subs(point).applyfunc(sp.simplify)
    nil = zero(at * at) and not zero(at)
    claim(f"{name}: d = {d_expected}, the page's eigenvalues are the roots, and at {point} the traceless part is "
          f"non-zero and nilpotent", sp.simplify(d - d_expected) == 0 and roots_ok and nil)
    return M


def main():
    print("closed blocks of a generator")
    # the two-site pair {|01>, |10>} in its four writings
    row("the horizon at N = 2, 3 (F2b corollary) lambda² + 4γλ + κ_N·J², companion (the nilpotency cannot fail)",
        companion(4 * gam, cN * J ** 2),
        4 * gam ** 2 - cN * J ** 2,
        [-2 * gam + sp.sqrt(4 * gam ** 2 - cN * J ** 2), -2 * gam - sp.sqrt(4 * gam ** 2 - cN * J ** 2)],
        {J: 2 * gam / sp.sqrt(cN)})
    row("the (1,1) pair at N = 2, Pauli book, lambda² + 4γλ + 16J² (PROOF_WEIGHT1_DEGENERACY), companion (the "
        "nilpotency cannot fail)",
        companion(4 * gam, 16 * J ** 2), 4 * gam ** 2 - 16 * J ** 2,
        [-2 * gam + 2 * I * sp.sqrt(4 * J ** 2 - gam ** 2), -2 * gam - 2 * I * sp.sqrt(4 * J ** 2 - gam ** 2)],
        {gam: 2 * J})
    row("the pair's own block in PROOF_MONOTONICITY_CPSI: with x = a - 1/2, ẋ = -4J·v, v̇ = 4J·x - 4γ·v, the matrix "
        "[[0, -4J], [4J, -4γ]] on (x, v), lambda = -2γ ± 2i√(4J² - γ²)",
        [[0, -4 * J], [4 * J, -4 * gam]], 4 * gam ** 2 - 16 * J ** 2,
        [-2 * gam + 2 * I * sp.sqrt(4 * J ** 2 - gam ** 2), -2 * gam - 2 * I * sp.sqrt(4 * J ** 2 - gam ** 2)],
        {gam: 2 * J})
    row("GAMMA0's z'' + 4γż + 4J²z = 0 (GAMMA0_IS_ALWAYS_THERE), companion (the nilpotency cannot fail)",
        companion(4 * gam, 4 * J ** 2),
        4 * gam ** 2 - 4 * J ** 2,
        [-2 * gam + 2 * sp.sqrt(gam ** 2 - J ** 2), -2 * gam - 2 * sp.sqrt(gam ** 2 - J ** 2)], {J: gam})
    row("the ring's {0,2} block [[0, B], [-B, -4γ]] (PROOF_RING_GAP_DOMINANCE)", [[0, B], [-B, -4 * gam]],
        4 * gam ** 2 - B ** 2,
        [-2 * gam + sp.sqrt(4 * gam ** 2 - B ** 2), -2 * gam - sp.sqrt(4 * gam ** 2 - B ** 2)], {B: 2 * gam})
    row("the proton doublet's (Z, Y) block [[0, 2J], [-2J, -2γ]] (PROTON_WIRE_CROSSING)",
        [[0, 2 * J], [-2 * J, -2 * gam]], gam ** 2 - 4 * J ** 2,
        [-gam + sp.sqrt(gam ** 2 - 4 * J ** 2), -gam - sp.sqrt(gam ** 2 - 4 * J ** 2)], {gam: 2 * J})
    Me = row("the N = 2 edge block [[-2γ₀, 2iJ], [2iJ, -2γ₁]] (PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE §(d))",
             [[-2 * gam0, 2 * I * J], [2 * I * J, -2 * gam1]], (gam1 - gam0) ** 2 - 4 * J ** 2,
             [-(gam0 + gam1) + sp.sqrt((gam1 - gam0) ** 2 - 4 * J ** 2),
              -(gam0 + gam1) - sp.sqrt((gam1 - gam0) ** 2 - 4 * J ** 2)], {gam1: gam0 + 2 * J})
    claim("the edge block's part at its point is [[2J, 2iJ], [2iJ, -2J]], as §(d) writes it",
          zero((Me - (Me.trace() / 2) * sp.eye(2)).subs(gam1, gam0 + 2 * J)
               - sp.Matrix([[2 * J, 2 * I * J], [2 * I * J, -2 * J]])))
    row("the sodium satellite pair -c·1 + [[-iG, κ], [κ, iG]] (THE_SODIUM_NUCLEUS_TWO_MIRRORS; G = 2q, κ = 3J₂; "
        "the frequency -iω dropped)", [[-cc - I * G, kap], [kap, -cc + I * G]], kap ** 2 - G ** 2,
        [-cc + sp.sqrt(kap ** 2 - G ** 2), -cc - sp.sqrt(kap ** 2 - G ** 2)], {G: kap})
    row("the neural pair -s·1 + [[x, -h], [h, -x]] (FLAVOR_RESOLVED_T2_INHERITANCE, x = Δ + c·a, h > 0)",
        [[-s + x, -h], [h, -s - x]], x ** 2 - h ** 2,
        [-s + sp.sqrt(x ** 2 - h ** 2), -s - sp.sqrt(x ** 2 - h ** 2)], {x: h})
    Mn = sp.Matrix([[x, -h], [h, -x]])
    kerp = sp.Matrix([1, 1])
    claim("the neural pair meets on both branches, x = h and x = -h (each a Jordan block, eigenvector weight one half on "
          "the two flavors), and at h = 0 = x its traceless part is zero: semisimple",
          zero(Mn.subs(x, h) * kerp) and zero(Mn.subs(x, -h) * sp.Matrix([1, -1]))
          and zero(Mn.subs(x, -h) ** 2) and not zero(Mn.subs(x, -h)) and zero(Mn.subs({x: 0, h: 0})))
    row("the sodium double-quantum pair -c·1 + [[-iG, -κ1], [-κ1, iG]] (G = 2q, κ1 = 3J₁; the frequency -2iω dropped)",
        [[-cc - I * G, -kap], [-kap, -cc + I * G]], kap ** 2 - G ** 2,
        [-cc + sp.sqrt(kap ** 2 - G ** 2), -cc - sp.sqrt(kap ** 2 - G ** 2)], {G: kap})

    print("\nreductions, toys and controls")
    row("the F86 toy [[-2γ₀(2k-1), iJg], [iJg, -2γ₀(2k+1)]] (PROOF_F86A_EP_MECHANISM; Jg = J·g_eff)",
        [[-2 * gam0 * (2 * k - 1), I * Jg], [I * Jg, -2 * gam0 * (2 * k + 1)]], 4 * gam0 ** 2 - Jg ** 2,
        [-4 * gam0 * k + sp.sqrt(4 * gam0 ** 2 - Jg ** 2), -4 * gam0 * k - sp.sqrt(4 * gam0 ** 2 - Jg ** 2)],
        {Jg: 2 * gam0})
    Ms = sp.Matrix([[zs, 1], [1, -zs]])
    for pt, ker in ((I, sp.Matrix([1, -I])), (-I, sp.Matrix([1, I]))):
        at = Ms.subs(zs, pt)
        claim(f"the seed toy [[s, 1], [1, -s]] (SEED_HOLONOMY_THIRD_CLOCK names s = i; -i by the same arithmetic): "
              f"eigenvalues ±√(s² + 1); at s = {pt} a "
              f"Jordan block at 0 with kernel vector {list(ker)}, self-orthogonal",
              sp.expand((Ms - lam * sp.eye(2)).det() - (lam ** 2 - zs ** 2 - 1)) == 0
              and zero(at * at) and not zero(at) and zero(at * ker) and sp.expand((ker.T * ker)[0]) == 0)
    Mr = -4 * sp.eye(2) + I * sp.Matrix([[u, eps], [eps, -u]])
    okr = True
    for pt in (I * eps, -I * eps):
        at = (Mr + 4 * sp.eye(2)).subs(u, pt)
        okr &= zero(at * at) and not zero(at) and sp.simplify((Mr.subs(u, pt) - lam * sp.eye(2)).det().subs(lam, -4)) == 0
    claim("the Route-B control -4·1 + i[[q-1, ε], [ε, -(q-1)]] (PROOF_ROUTE_B_N4_SELF_FOLD): discriminant "
          "-4((q-1)² + ε²); at q - 1 = ±iε a Jordan block at lambda = -4", okr
          and sp.expand(sp.discriminant((Mr - lam * sp.eye(2)).det(), lam) + 4 * (u ** 2 + eps ** 2)) == 0)

    print("\na Hamiltonian")
    X = sp.Matrix([[0, 1], [1, 0]]); Y = sp.Matrix([[0, -I], [I, 0]]); Zm = sp.diag(1, -1); Id = sp.eye(2)
    kron = sp.kronecker_product
    Hd = J * (kron(X, X) + kron(Y, Y)) / 2 + I * g * (kron(Zm, Id) - kron(Id, Zm)) / 2
    Hse = Hd.extract([1, 2], [1, 2])
    claim("F155's PT dimer H = J(XX + YY)/2 + ig(Z0 - Z1)/2: single-excitation block [[ig, J], [J, -ig]] on |01>, |10>, "
          "eigenvalues ±√(J² - g²), a Jordan block at g = J",
          zero(Hse - sp.Matrix([[I * g, J], [J, -I * g]])) and zero(Hse * Hse - (J ** 2 - g ** 2) * sp.eye(2))
          and zero(Hse.subs(g, J) ** 2) and not zero(Hse.subs(g, J)))
    Hn = Hd.subs({J: 1, g: 1})
    Ln = -I * (kron(Hn, sp.eye(4)) - kron(sp.eye(4), sp.conjugate(Hn)))
    null = [16 - (Ln ** p).rank() for p in (1, 2, 3, 4)]
    claim(f"F155's no-jump generator rho -> -i(H·rho - rho·H†), the superoperator -i(H⊗1 - 1⊗H̄), at g = J carries a "
          f"Jordan block of size 3 at 0 and none larger (nullities of G, G², G³, G⁴: {null})",
          null[2] > null[1] and null[3] == null[2])

    print("\nthe writing")
    a_, b_ = sp.symbols("a b")
    Dg = sp.diag(1, I)
    vr = sp.Matrix([1, -1])
    claim("the similarity diag(1, i) carries [[a, g], [-g, b]] to [[a, ig], [ig, b]]; the ring's merged vector (1, -1) "
          "at B = 2γ has vᵀv = 2, its image under diag(1, i)⁻¹, (1, i), has vᵀv = 0: self-orthogonality belongs to the "
          "complex-symmetric writing, not to the block",
          zero(Dg.inv() * sp.Matrix([[a_, g], [-g, b_]]) * Dg - sp.Matrix([[a_, I * g], [I * g, b_]]))
          and zero((sp.Matrix([[0, 2 * gam], [-2 * gam, -4 * gam]]) + 2 * gam * sp.eye(2)) * vr)
          and (vr.T * vr)[0] == 2 and sp.expand(((Dg.inv() * vr).T * (Dg.inv() * vr))[0]) == 0)

    held = sum(CLAIMS)
    print(f"\n{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS):
        sys.exit(1)


if __name__ == "__main__":
    main()
