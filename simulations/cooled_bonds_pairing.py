"""Which two-site bonds pair under cooling alone: exact checks with controls.

experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md. Every site carries cooling sigma^- = |0><1| at its own rate and nothing
else; H is a chain of one bond repeated, J_ab·sigma_a(x)sigma_b over a, b in {X, Y} plus Δ·ZZ, or a bond with a mixed
letter (XZ, YZ, ZX, ZY) among them; sections F and G also take the questions to rings. Pairing means the spectrum of L is symmetric about −c/2, c = Σ rates. Rows print PASS or FAIL; the
script exits nonzero on any FAIL; --deep adds the exact
polynomials named below as such. With --map it first prints the float survey the document quotes (a reading, not a
gate): random integer couplings over every subset of the nine two-letter bonds, N = 2, 3, 4.

  A  The turn lemma: for every real 2x2 J there is a rotation R with R·Jᵀ·Rᵀ = F·J·F, F = diag(1, −1) (Jᵀ and FJF share
     trace, antisymmetric part and the length of the traceless symmetric part). Checked on six J.
  B  The carrier at rates the chain's reversal P keeps: E = (Ad X^N ∘ Ad(U·P))⁻¹ ∘ S_IZ, S_IZ: rho -> rho·Z^N, U the
     uniform turn about Z realizing R, gives E·L·E⁻¹ = −L† − c exactly (F137's Ad_V∘S_XY with V = (U·P)†); controls:
     rates P does not keep.
  C  The Theodolite's elements (exact, read off the letters and verified): X/Y letters alone at any rates (a staggered
     A·rhoᵀ·B); with ZZ, a real or an imaginary X/Y part (A·rho·B), J symmetric (a = 0) or traceless (s = 0) at
     unequal rates, w = 0 (XXZ with an axial DM term), and any J on one bond
     (turned); silent where s, a and w are all nonzero at rates the reversal
     does not keep. J = s·1 + a·ε + T, ε = [[0, 1], [−1, 0]], T traceless symmetric, w = T₁₁ + i·T₁₂.
  D  Breaks, exact: an odd centred power sum Tr((L − tr L/dim)^k) that is not zero rules out every element; the first
     such k is asserted: 7 where s, a, w are all nonzero (XX + 2XY + ZZ at (2, 1, 1), (1, 2, 3) and, N = 4,
     (1, 2, 3, 4); 2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3)), 5 for a mixed letter (XZ; ZX + XX). Control: zero up to
     k = 7 at rates the reversal keeps.
  E  The seventh sum in closed form: Tr((L + c/2)^7) = −7·2¹⁴·4^(N−3)·Δ·s·a·|w|²·(γ₁γ₂ − γ_{N−1}γ_N), exact as a
     polynomial in the rates at N = 2 (zero, one J), 3 and 4 for three J (interpolated modulo primes, three surplus
     points asserted to fit, the CRT value unchanged by a fifth prime), at random rates at N = 5 and 6 (modulo a
     prime, every entry an integer); controls a = 0 and Δ = 0 give zero; where the seventh is silent off the mirror,
     rates (1, 6, 2, 3) and (1, 0, 2), the ninth is not.
  F  Bonds with a mixed letter: Tr((L + c/2)^3) = 0 for every bond of the nine letters (exact polynomial, N = 2, 3);
     Tr((L + c/2)^5) = 10·4^(N−2)·[3·S − 16·Im(μνw̄)·Σ γ_b γ_{b+1}], μ = c_XZ + i·c_YZ, ν = c_ZX + i·c_ZY,
     S = Σ_{m=2}^{N−1} γ_m·|μγ_{m+1} + νγ_{m−1}|² + |μ|²γ₁γ₂² + |ν|²γ_{N−1}²γ_N, exact as a polynomial in all nine
     coefficients and the rates at N = 2 (and 3, 4, 5 with --deep; the basis has the mixed letters in degree 2 or 4
     and 25 surplus points must fit), at random points modulo a prime at N = 3 to 6; where anisotropy cancels it
     (XZ + ZX − 3XY − 3YX at (8, 8)) the seventh is not zero; a site without cooling can empty S and the spectrum
     then can pair: XZ at (1, 0) and at (1, 0, 1), exact by the characteristic polynomial (control (1, 1, 0)).
     Every interpolation carries a guard prime.
  G  Rings: the seventh sum is zero at every rate on rings of 3 and 4 (and 5 with --deep) for two bonds, and the
     two bonds' (Im z³, Im z·|w|²) are independent, which with the page's parity count makes it zero on the triangle
     for every bond; at one random class bond on rings of 3 to 6 the third, fifth and seventh sums are zero and the
     ninth is not, modulo a prime; controls with a different bond on each edge: the triangle's fifth is not zero, on
     rings of 4 and 5 the fifth is and the seventh is not; the fifth
     sum of a bond with a mixed letter on the ring of four is the closed form without end terms (one random point
     modulo a prime), on the triangle for one bond it is not, and XZ − ZX at equal rates has the fifth and seventh exactly zero
     (zero modulo five primes, with a norm bound below their product) and the ninth not; on the ring of four
     Tr((L + c/2)^9) = −9·2¹²·Δ·s·a·|w|²·(γ₁ − γ₃)(γ₂ − γ₄)(Q(γ₁, γ₃) − Q(γ₂, γ₄)), Q(x, y) = 5x² + xy + 5y², for
     two bonds, and at (5, 3, 9, 10), where Q(5, 9) = Q(3, 10), the ninth is exactly zero and the eleventh not; on
     the triangle the ninth is a nonzero constant times (γ₁ − γ₂)(γ₁ − γ₃)(γ₂ − γ₃)(γ₁ + γ₂ + γ₃), the constants of
     the two bonds not in the ratio of Δ·s·a·|w|².
  H  The ninth cumulant on the chain: κ₉ = 144·Δ·s·a·|w|²·[Q + E(γ₁, γ₂, γ₃) − E(γ_N, γ_{N−1}, γ_{N−2})
     + 16(e(γ₁, γ₂, γ₃) − e(γ_N, γ_{N−1}, γ_{N−2}))] with the bulk Q, the quartic end E and the quadratic end e of
     the document, p₉ = 4^N·κ₉ + 36·p₇·p₂/4^N: exact as a polynomial in the rates (degrees 2 and 4, the only ones
     the phase count allows) at N = 3 and 4 for six bonds, with --deep also at N = 5 for four and at N = 6 for one,
     which with the cofactor count fixes every window; at one random bond with Δ, s, a, w nonzero and profile with
     unequal end products, modulo a prime, at N = 4, 5, 6 (and 7 with --deep); at
     N = 3 with no cooling in the middle it is a multiple of γ₁² − γ₃² whose other factor is positive, for every bond (symbolic); at
     two profiles with equal end products where it vanishes, the seventh and ninth are exactly zero and the eleventh
     is not.
"""
from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding='utf-8')

from framework import cube_tables as ct  # noqa: E402
from framework.diagnostics import cube_letters as cl  # noqa: E402

X = sp.Matrix([[0, 1], [1, 0]]); Y = sp.Matrix([[0, -sp.I], [sp.I, 0]]); Z = sp.diag(1, -1); I2 = sp.eye(2)
LOW = ct.LOWER
LET = {'X': X, 'Y': Y, 'Z': Z}
F = sp.diag(1, -1)
FAILURES = []


def check(name, ok, detail=''):
    print(f"  {'PASS' if ok else 'FAIL'}: {name}" + (f"  [{detail}]" if detail else ''))
    if not ok:
        FAILURES.append(name)


def kron(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = sp.kronecker_product(out, m)
    return out


def emb(ops, n):
    return kron(*[ops.get(l, I2) for l in range(n)])


def chain_terms(n, bond):
    """bond: {two-letter string: coefficient}, repeated on every bond of the open chain."""
    return [(s, (b, b + 1), c) for b in range(n - 1) for s, c in bond.items()]


def hamiltonian(n, bond):
    H = sp.zeros(2 ** n, 2 ** n)
    for s, (a, b), c in chain_terms(n, bond):
        H += c * emb({a: LET[s[0]], b: LET[s[1]]}, n)
    return H


def generator(H, rates):
    n = len(rates); d = 2 ** n; Id = sp.eye(d)
    L = -sp.I * (kron(H, Id) - kron(Id, H.T))
    for l, r in enumerate(rates):
        A = emb({l: LOW}, n); AdA = A.H * A
        L += r * (kron(A, A.conjugate()) - sp.Rational(1, 2) * (kron(AdA, Id) + kron(Id, AdA.T)))
    return L


def bond_of(J, zz):
    out = {a + b: J[i, j] for i, a in enumerate('XY') for j, b in enumerate('XY') if J[i, j] != 0}
    if zz:
        out['ZZ'] = zz
    return out


def turn_for(J):
    c, s = sp.symbols('c s', real=True)
    R = sp.Matrix([[c, -s], [s, c]])
    sols = sp.solve(list(R * J.T * R.T - F * J * F) + [c ** 2 + s ** 2 - 1], [c, s], dict=True)
    if not sols:
        return None
    sol = sols[0]
    if s not in sol:                       # a free angle (J commutes with every rotation's action): take s = 0
        sol = {c: sp.sympify(sol.get(c, 1)).subs(s, 0), s: 0}
        if sol[c] ** 2 != 1:
            sol[c] = 1
    return sol[c], sol[s]


def reversal(n):
    P = sp.zeros(2 ** n, 2 ** n)
    for k in range(2 ** n):
        bits = [(k >> (n - 1 - l)) & 1 for l in range(n)]
        P[sum(b << (n - 1 - l) for l, b in enumerate(reversed(bits))), k] = 1
    return P


def carrier_reflects(n, J, zz, rates):
    cth, sth = turn_for(J)
    W = emb({l: sp.diag(1, cth + sp.I * sth) for l in range(n)}, n) * reversal(n)
    Xn, Zn = emb({l: X for l in range(n)}, n), emb({l: Z for l in range(n)}, n)
    Ad = lambda A: kron(A, A.conjugate())                       # rho -> A·rho·A†, row-major
    E = (Ad(Xn) * Ad(W)).inv() * kron(sp.eye(2 ** n), Zn.T)    # S_IZ: rho -> rho·Z^N
    L = generator(hamiltonian(n, bond_of(J, zz)), rates)
    res = E * L + (L.H + sum(rates) * sp.eye(4 ** n)) * E
    return res.applyfunc(lambda e: sp.simplify(sp.expand(e))).is_zero_matrix


def setting(n, bond, rates):
    return ct.Setting(n, chain_terms(n, bond), [((l,), LOW, r) for l, r in enumerate(rates)])


def letters(s, rel):
    rd = cl.read(s, rel)
    if rd is not None and cl.verify(s, rel, rd.maps):
        return 'transposed' if 'ᵀ' in cl.describe(rd)[0] else 'plain'
    return None


def site_turns(s, rel):
    """The Theodolite's left turn x_l per site (site factor [[0, x], [1, 0]]), exact, or None."""
    rd = cl.read(s, rel)
    if rd is None or not cl.verify(s, rel, rd.maps):
        return None
    out = []
    for line in cl.describe(rd)[1:]:
        if line.startswith('site ') and 'A = [[0, x], [1, 0]] with x = ' in line:
            out.append(sp.sympify(line.split('with x = ')[1].split(';')[0], locals={'I': sp.I}))
    return out


def first_odd_trace(n, bond, rates, kmax=7):
    L = generator(hamiltonian(n, bond), rates)
    dim = L.shape[0]
    M = DomainMatrix.from_Matrix((dim * L - L.trace() * sp.eye(dim)).applyfunc(sp.expand)).convert_to(sp.QQ_I)
    P = M
    for k in range(1, kmax + 1):
        if k > 1:
            P = P * M
        if k % 2:
            tr = sum(P.to_Matrix()[i, i] for i in range(dim))
            if tr != 0:
                return k
    return None


def survey():
    """The float reading the document quotes: per (N, rates), how many random bonds pair, by letter class."""
    P = {k: np.array(v, dtype=complex) for k, v in ct.PAULI.items()}
    low = np.array(LOW, dtype=complex)
    letters9 = [a + b for a in 'XYZ' for b in 'XYZ']
    rng = np.random.default_rng(1)
    for n, rates in [(2, [1, 2]), (3, [1, 1, 1]), (3, [1, 2, 1]), (3, [1, 2, 3]), (4, [1, 2, 2, 1]), (4, [1, 2, 3, 5])]:
        d = 2 ** n; Id = np.eye(d); tally = {}; margin = {True: 0.0, False: np.inf}
        for k in range(1, 10):
            for sub in itertools.combinations(letters9, k):
                H = np.zeros((d, d), complex)
                coef = {s: rng.integers(1, 4) * rng.choice([-1, 1]) for s in sub}   # one bond, repeated
                for b in range(n - 1):
                    for s in sub:
                        ops = [np.eye(2)] * n
                        ops[b], ops[b + 1] = P[s[0]], P[s[1]]
                        term = ops[0]
                        for o in ops[1:]:
                            term = np.kron(term, o)
                        H += coef[s] * term
                L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
                for l, r in enumerate(rates):
                    ops = [np.eye(2)] * n; ops[l] = low
                    A = ops[0]
                    for o in ops[1:]:
                        A = np.kron(A, o)
                    AdA = A.conj().T @ A
                    L += r * (np.kron(A, A.conj()) - 0.5 * np.kron(AdA, Id) - 0.5 * np.kron(Id, AdA.T))
                ev = np.linalg.eigvals(L)
                dist = max(np.min(np.abs(ev + sum(rates) + e)) for e in ev)
                mixed = any(set(s) & {'Z'} and s != 'ZZ' for s in sub)
                Jc = {k: coef.get(k, 0) for k in ('XX', 'XY', 'YX', 'YY')}
                s_, a_ = (Jc['XX'] + Jc['YY']) / 2, (Jc['XY'] - Jc['YX']) / 2
                w_ = complex((Jc['XX'] - Jc['YY']) / 2, (Jc['XY'] + Jc['YX']) / 2)
                cls = ('mixed letter' if mixed else 'Δ = 0' if 'ZZ' not in sub
                       else 'Δ ≠ 0, s·a·w = 0' if s_ * a_ * abs(w_) == 0 else 'Δ ≠ 0, s, a, w ≠ 0')
                key = (cls, dist < 1e-6)
                tally[key] = tally.get(key, 0) + 1
                margin[key[1]] = (max if key[1] else min)(margin[key[1]], dist)
        print(f'  N = {n}, rates {rates}:', {f'{c} {"pairs" if p else "breaks"}': v for (c, p), v in sorted(tally.items())},
              f'(largest distance among the pairing {margin[True]:.1e}, smallest among the breaking {margin[False]:.2f})')


if '--map' in sys.argv:
    print("map (float reading, one random bond in ±{1, 2, 3} repeated along the chain, every subset of the nine letters)")
    survey()

print("A  the turn lemma")
for J in [sp.Matrix([[1, 2], [0, 0]]), sp.Matrix([[3, -1], [2, 5]]), sp.Matrix([[1, 0], [0, -4]]),
          sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[2, 3], [-3, 2]]), sp.Matrix([[0, 2], [3, 0]])]:
    t = turn_for(J)
    R = sp.Matrix([[t[0], -t[1]], [t[1], t[0]]]) if t else None
    check(f"J = {J.tolist()}: a rotation R with R·Jᵀ·Rᵀ = F·J·F exists",
          t is not None and sp.simplify(t[0] ** 2 + t[1] ** 2) == 1
          and (R * J.T * R.T - F * J * F).applyfunc(sp.simplify).is_zero_matrix, f'cos, sin = {t}')

print("B  the carrier: chain reversal, uniform turn, S_IZ")
J1, J2 = sp.Matrix([[1, 2], [0, 0]]), sp.Matrix([[3, -1], [2, 5]])
check("XX + 2XY + ZZ, N = 3, rates (1, 2, 1): E·L·E⁻¹ = −L† − c", carrier_reflects(3, J1, 1, [1, 2, 1]))
check("XX + 2XY + ZZ, N = 3, rates (1, 1, 1)", carrier_reflects(3, J1, 1, [1, 1, 1]))
check("J = [[3, −1], [2, 5]] − 2ZZ, N = 3, rates (2, 7, 2)", carrier_reflects(3, J2, -2, [2, 7, 2]))
check("XX + 2XY + ZZ, N = 4, rates (1, 2, 2, 1)", carrier_reflects(4, J1, 1, [1, 2, 2, 1]))
check("controls: rates the reversal does not keep, (2, 1, 1) and (2, 7, 3): it is not",
      not carrier_reflects(3, J1, 1, [2, 1, 1]) and not carrier_reflects(3, J2, -2, [2, 7, 3]))

print("C  the Theodolite's elements")
xy_only = {'XX': 1, 'XY': 2, 'YX': 3, 'YY': 5}
check("X/Y letters alone, XX + 2XY + 3YX + 5YY, rates (1, 2, 3): a transposed sandwich reflects ('refl')",
      letters(setting(3, xy_only, [1, 2, 3]), 'refl') == 'transposed')
check("real part with ZZ, XX + 3YY + ZZ, rates (1, 2, 3): a sandwich reflects",
      letters(setting(3, {'XX': 1, 'YY': 3, 'ZZ': 1}, [1, 2, 3]), 'refl') == 'plain')
check("imaginary part with ZZ, XY + 3YX + ZZ, rates (1, 2, 3): a sandwich reflects",
      letters(setting(3, {'XY': 1, 'YX': 3, 'ZZ': 1}, [1, 2, 3]), 'refl') == 'plain')
xs = site_turns(setting(3, {'XX': 1, 'XY': 2, 'YX': 2, 'ZZ': 1}, [1, 2, 3]), 'refl')
check("J symmetric (a = 0), XX + 2XY + 2YX + ZZ, rates (1, 2, 3): one turn on every site, x = ±(−1 + 4i)/√17",
      xs is not None and len(xs) == 3 and all(sp.simplify(x - xs[0]) == 0 for x in xs)
      and sp.simplify(xs[0] ** 2 - (-1 + 4 * sp.I) ** 2 / 17) == 0, f'x = {xs}')
xs = site_turns(setting(3, {'XY': 1, 'YX': 3, 'ZZ': 1}, [1, 2, 3]), 'refl')
check("imaginary part with ZZ, XY + 3YX + ZZ: the turn alternates by a half turn, x = (−1, 1, −1) up to sign",
      xs is not None and [sp.simplify(x / xs[0]) for x in xs] == [1, -1, 1], f'x = {xs}')
check("J traceless (s = 0), 2XX + 3XY − YX − 2YY + ZZ, rates (1, 2, 3): a turned sandwich reflects",
      letters(setting(3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': -2, 'ZZ': 1}, [1, 2, 3]), 'refl') is not None)
check("w = 0 (XXZ with an axial DM term), XX + YY + XY − YX + ZZ, rates (1, 2, 3): a turned sandwich reflects",
      letters(setting(3, {'XX': 1, 'YY': 1, 'XY': 1, 'YX': -1, 'ZZ': 1}, [1, 2, 3]), 'refl') is not None)
check("s, a, w all nonzero with ZZ on one bond, XX + 2XY + ZZ, N = 2, rates (1, 2): a turned sandwich reflects",
      letters(setting(2, {'XX': 1, 'XY': 2, 'ZZ': 1}, [1, 2]), 'refl') is not None)
check("s, a, w all nonzero with ZZ, N = 3, rates the reversal keeps (1, 2, 1): the letters are silent (no site map "
      "carries the reversal)", letters(setting(3, {'XX': 1, 'XY': 2, 'ZZ': 1}, [1, 2, 1]), 'refl') is None)
check("s, a, w all nonzero with ZZ, N = 3: the letters are silent (XX + 2XY + ZZ at (2, 1, 1), "
      "2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3))",
      letters(setting(3, {'XX': 1, 'XY': 2, 'ZZ': 1}, [2, 1, 1]), 'refl') is None
      and letters(setting(3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': 2, 'ZZ': 1}, [1, 2, 3]), 'refl') is None)

print("D  breaks, exact")
class_bond = {'XX': 1, 'XY': 2, 'ZZ': 1}
for n, bond, rates, name, first in [
        (3, class_bond, [2, 1, 1], 'XX + 2XY + ZZ, rates (2, 1, 1)', 7),
        (3, class_bond, [1, 2, 3], 'XX + 2XY + ZZ, rates (1, 2, 3)', 7),
        (4, class_bond, [1, 2, 3, 4], 'XX + 2XY + ZZ, N = 4, rates (1, 2, 3, 4)', 7),
        (3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': 2, 'ZZ': 1}, [1, 2, 3], '2XX + 3XY − YX + 2YY + ZZ, rates (1, 2, 3)', 7),
        (3, {'XZ': 1}, [1, 1, 1], 'XZ, rates (1, 1, 1)', 5),
        (3, {'ZX': 1, 'XX': 1}, [1, 2, 3], 'ZX + XX, rates (1, 2, 3)', 5)]:
    k = first_odd_trace(n, bond, rates)
    check(f"{name}: the first odd centred power sum that is not zero is the {first}th", k == first, f'k = {k}')
check("control: XX + 2XY + ZZ at rates (1, 2, 1) has every odd power sum zero up to k = 7",
      first_odd_trace(3, class_bond, [1, 2, 1]) is None)

print("E  the seventh sum in closed form")
PRIMES = [p_ for p_ in sp.primerange(900000, 1048576) if p_ % 4 == 1][-5:]   # float BLAS stays exact below 2^20


def modp_generator(n, bond, rates, p, ring=False):
    """2L + c over GF(p), i the square root of −1 (p = 1 mod 4); row-major vec, rho -> A rho B is A ⊗ Bᵀ; on a ring
    the bond also joins the last site to the first. A list of bonds gives each edge its own, in order."""
    i = next(r for r in (pow(g, (p - 1) // 4, p) for g in range(2, p)) if r * r % p == p - 1)
    let = {'X': np.array([[0, 1], [1, 0]]), 'Y': np.array([[0, p - i], [i, 0]]), 'Z': np.array([[1, 0], [0, p - 1]])}
    low = np.array([[0, 1], [0, 0]])

    def emb_(ops):
        out = np.array([[1]], dtype=np.int64)
        for l in range(n):
            out = np.kron(out, ops.get(l, np.eye(2, dtype=np.int64))) % p
        return out
    d = 2 ** n
    H = np.zeros((d, d), dtype=np.int64)
    for e, (b, c) in enumerate([(b, b + 1) for b in range(n - 1)] + ([(n - 1, 0)] if ring else [])):
        for lab, v in (bond[e] if isinstance(bond, list) else bond).items():
            H = (H + (v % p) * emb_({b: let[lab[0]], c: let[lab[1]]})) % p
    Id = np.eye(d, dtype=np.int64)
    M = (-2 * i) * ((np.kron(H, Id) - np.kron(Id, H.T)) % p) % p
    for l, g in enumerate(rates):
        A = emb_({l: low}); AdA = A.T @ A; g %= p
        M = (M + 2 * g * np.kron(A, A) - g * (np.kron(AdA, Id) + np.kron(Id, AdA.T)) + g * np.eye(d * d, dtype=np.int64)) % p
    return M


def modp_trace_power(M, k, p):
    R = M
    for _ in range(k - 1):
        R = (R.astype(np.float64) @ M.astype(np.float64) % p).astype(np.int64)
    return int(np.trace(R) % p)


def modp_solve(A, b, p):
    A = np.concatenate([np.array(A, dtype=np.int64) % p, np.array(b, dtype=np.int64)[:, None] % p], axis=1)
    n = A.shape[1] - 1
    for c in range(n):
        k = c + int(np.nonzero(A[c:, c])[0][0])
        A[[c, k]] = A[[k, c]]
        A[c] = A[c] * pow(int(A[c, c]), p - 2, p) % p
        f = A[:, c].copy(); f[c] = 0
        A = (A - (f[:, None] * A[c][None, :]) % p) % p
    assert not A[n:].any(), 'the surplus points do not fit the polynomial'
    return [int(x) for x in A[:n, n]]


def seventh_sum_polynomial(n, bond, k=7):
    """Tr((2L + c)^k) as a polynomial in the rates (degree at most k), interpolated modulo each prime; returns the
    CRT polynomial from four primes and whether a fifth prime changes it."""
    gs = sp.symbols(f'g1:{n + 1}')
    monos = sorted(e for e in itertools.product(range(k + 1), repeat=n) if sum(e) <= k)
    rng = np.random.default_rng(1)
    pts = [tuple(int(x) for x in rng.integers(1, 10 ** 5, n)) for _ in range(len(monos) + 3)]
    per_p = []
    for p in PRIMES:
        A = [[int(np.prod([pow(x, e_, p) for x, e_ in zip(pt, e)], dtype=object) % p) for e in monos] for pt in pts]
        per_p.append(modp_solve(A, [modp_trace_power(modp_generator(n, bond, pt, p), k, p) for pt in pts], p))

    def crt(m):
        out = 0
        for j, e in enumerate(monos):
            r, mod = sp.ntheory.modular.crt(PRIMES[:m], [c[j] for c in per_p[:m]])
            r = int(r) - (int(mod) if int(r) > int(mod) // 2 else 0)
            out += r * sp.Mul(*[g ** e_ for g, e_ in zip(gs, e)])
        return sp.expand(out)
    four = crt(4)
    return four, gs, sp.expand(four - crt(5)) == 0


def closed_form(n, bond, gs):
    J = {lab: sp.Integer(bond.get(lab, 0)) for lab in ('XX', 'XY', 'YX', 'YY', 'ZZ')}
    s_, a_ = (J['XX'] + J['YY']) / 2, (J['XY'] - J['YX']) / 2
    w2 = ((J['XX'] - J['YY']) / 2) ** 2 + ((J['XY'] + J['YX']) / 2) ** 2
    return -7 * 2 ** 14 * 4 ** (n - 3) * J['ZZ'] * s_ * a_ * w2 * (gs[0] * gs[1] - gs[-2] * gs[-1]) * 2 ** 7


bonds_E = [({'XX': 1, 'XY': 2, 'ZZ': 1}, 'XX + 2XY + ZZ'), ({'XX': 3, 'XY': -1, 'YX': 2, 'YY': 5, 'ZZ': -2},
           '3XX − XY + 2YX + 5YY − 2ZZ'), ({'XX': -2, 'XY': 1, 'YX': 3, 'YY': 1, 'ZZ': 3}, '−2XX + XY + 3YX + YY + 3ZZ')]
for n in (2, 3, 4):
    for bond, name in (bonds_E if n > 2 else bonds_E[:1]):
        P, gs, stable = seventh_sum_polynomial(n, bond)
        target = 0 if n == 2 else closed_form(n, bond, gs)
        check(f"N = {n}, {name}: the seventh sum is the closed form as a polynomial in the rates"
              + (" (zero)" if n == 2 else ""), stable and sp.expand(P - target) == 0, f'{sp.factor(P)}')
for bond, name in [({'XX': 1, 'XY': 2, 'YX': 2, 'ZZ': 1}, 'a = 0, XX + 2XY + 2YX + ZZ'),
                   ({'XX': 1, 'XY': 2}, 'Δ = 0, XX + 2XY')]:
    P, gs, stable = seventh_sum_polynomial(3, bond)
    check(f"control, N = 3, {name}: the seventh sum is zero at every rate", stable and P == 0)
rng = np.random.default_rng(3)
for n in (5, 6):
    for bond, name in bonds_E[:2]:
        g = [int(x) for x in rng.integers(1, 1000, n)]
        p = PRIMES[0]
        v = closed_form(n, bond, sp.symbols(f'g1:{n + 1}')).subs(dict(zip(sp.symbols(f'g1:{n + 1}'), g)))
        check(f"N = {n}, {name}, rates {tuple(g)}: the seventh sum equals the closed form modulo {p}",
              modp_trace_power(modp_generator(n, bond, g, p), 7, p) == int(v) % p)
M = modp_generator(4, bonds_E[0][0], [1, 6, 2, 3], PRIMES[0])
check("N = 4, XX + 2XY + ZZ, rates (1, 6, 2, 3) (γ₁γ₂ = γ₃γ₄, not mirror-symmetric): the seventh sum is zero and "
      "the ninth is not", modp_trace_power(M, 7, PRIMES[0]) == 0 and modp_trace_power(M, 9, PRIMES[0]) != 0)
M = modp_generator(3, bonds_E[0][0], [1, 0, 2], PRIMES[0])
check("N = 3, XX + 2XY + ZZ, rates (1, 0, 2) (no cooling on the middle site): the seventh sum is zero and the ninth "
      "is not", modp_trace_power(M, 7, PRIMES[0]) == 0 and modp_trace_power(M, 9, PRIMES[0]) != 0)


def interpolate(monos, syms, evaluate, nprimes=3, surplus=5, seed=1):
    """A polynomial in syms with exactly the monomials monos, from its values modulo primes (evaluate(point, p)) by
    elimination and CRT; the surplus points must fit, which they would not if a monomial outside monos were there.
    One guard prime more is computed, and the result is returned only if the CRT value does not change with it."""
    rng = np.random.default_rng(seed)
    pts = [tuple(int(x) for x in rng.integers(1, 10 ** 6, len(syms))) for _ in range(len(monos) + surplus)]
    per = []
    for p in PRIMES[:nprimes + 1]:
        A = [[int(np.prod([pow(x, e_, p) for x, e_ in zip(pt, e) if e_], dtype=object) % p) for e in monos]
             for pt in pts]
        per.append(modp_solve(A, [evaluate(pt, p) for pt in pts], p))

    def crt(m):
        out = 0
        for j, e in enumerate(monos):
            r, mod = sp.ntheory.modular.crt(PRIMES[:m], [c[j] for c in per[:m]])
            r = int(r) - (int(mod) if int(r) > int(mod) // 2 else 0)
            out += r * sp.Mul(*[v ** e_ for v, e_ in zip(syms, e)])
        return sp.expand(out)
    out = crt(nprimes)
    assert sp.expand(out - crt(nprimes + 1)) == 0, 'the guard prime changed the CRT value'
    return out


def rate_polynomial(n, bond, k, ring=False):
    """Tr((2L + c)^k) for one fixed bond as a polynomial in the rates (all degrees up to k)."""
    gs = sp.symbols(f'g1:{n + 1}', real=True)
    monos = [e for e in itertools.product(range(k + 1), repeat=n) if sum(e) <= k]
    return interpolate(monos, gs, lambda pt, p: modp_trace_power(modp_generator(n, bond, pt, p, ring), k, p)), gs


def float_generator(n, bond, rates, ring=False):
    """2L + c in double precision, for the norm bound only."""
    let = {'X': np.array([[0, 1], [1, 0]], complex), 'Y': np.array([[0, -1j], [1j, 0]]),
           'Z': np.array([[1, 0], [0, -1]], complex)}
    low = np.array([[0, 1], [0, 0]], complex)

    def emb_(ops):
        out = np.array([[1]], complex)
        for l in range(n):
            out = np.kron(out, ops.get(l, np.eye(2)))
        return out
    d = 2 ** n
    H = np.zeros((d, d), complex)
    for b, c in [(b, b + 1) for b in range(n - 1)] + ([(n - 1, 0)] if ring else []):
        for lab, v in bond.items():
            H += v * emb_({b: let[lab[0]], c: let[lab[1]]})
    Id = np.eye(d)
    M = -2j * (np.kron(H, Id) - np.kron(Id, H.T)) + sum(rates) * np.eye(d * d)
    for l, g in enumerate(rates):
        A = emb_({l: low}); AdA = A.T @ A
        M += g * (2 * np.kron(A, A) - np.kron(AdA, Id) - np.kron(Id, AdA.T))
    return M


def trace_is_exactly_zero(n, bond, rates, k, ring=False):
    """Tr((2L + c)^k) is an integer (a real Gaussian integer); it is zero if it vanishes modulo every prime and
    4^N·‖2L + c‖₂^k, which bounds its size, stays below half the primes' product."""
    if any(modp_trace_power(modp_generator(n, bond, rates, p, ring), k, p) for p in PRIMES):
        return False
    bound = 4 ** n * (1.001 * np.linalg.norm(float_generator(n, bond, rates, ring), 2) + 1) ** k
    return bound < np.prod([float(p) for p in PRIMES]) / 2


DEEP = '--deep' in sys.argv
MIXED, CLASS = ['XZ', 'YZ', 'ZX', 'ZY'], ['XX', 'XY', 'YX', 'YY', 'ZZ']


def fifth_closed_form(n, c, g, ring=False):
    """Tr((2L + c)^5) = 2⁵·10·4^(N−2)·[3·S − 16·Im(μνw̄)·Σ γ_b γ_{b+1}]; on the ring S has no end terms."""
    mu, nu = c['XZ'] + sp.I * c['YZ'], c['ZX'] + sp.I * c['ZY']
    w = sp.Rational(1, 2) * (c['XX'] - c['YY']) + sp.I * sp.Rational(1, 2) * (c['XY'] + c['YX'])
    at = (lambda m: g[m % n]) if ring else (lambda m: g[m])
    inner = range(n) if ring else range(1, n - 1)
    S = sum(at(m) * sp.expand((mu * at(m + 1) + nu * at(m - 1)) * sp.conjugate(mu * at(m + 1) + nu * at(m - 1)))
            for m in inner)
    if not ring:
        S += sp.expand(mu * sp.conjugate(mu)) * g[0] * g[1] ** 2 + sp.expand(nu * sp.conjugate(nu)) * g[-2] ** 2 * g[-1]
    cross = sum(at(m) * at(m + 1) for m in (range(n) if ring else range(n - 1)))
    return sp.expand(320 * 4 ** (n - 2) * (3 * S - 16 * sp.im(sp.expand(mu * nu * sp.conjugate(w))) * cross))


def ints(c):
    return {k_: sp.Integer(v) for k_, v in c.items()}


def exponents(nvars, total):
    """Every exponent tuple of nvars variables with the given total degree."""
    if nvars == 1:
        return [(total,)]
    return [(k,) + rest for k in range(total + 1) for rest in exponents(nvars - 1, total - k)]


def pairs_exactly(n, bond, rates):
    """The characteristic polynomial of L, shifted by c/2, has only terms of the degree's parity: the spectrum is
    symmetric about −c/2, exact over the Gaussian rationals."""
    L = generator(hamiltonian(n, bond), rates)
    x = sp.Symbol('x')
    P = DomainMatrix.from_Matrix(L.applyfunc(sp.expand)).convert_to(sp.QQ_I).charpoly()
    Q = sp.Poly(sp.expand(sum(sp.QQ_I.to_sympy(c_) * (x - sp.Rational(sum(rates), 2)) ** (len(P) - 1 - j)
                               for j, c_ in enumerate(P))), x)
    d = len(P) - 1
    return all(c_ == 0 for (k,), c_ in Q.terms() if (d - k) % 2)


print("F  bonds with a mixed letter: the fifth sum in closed form")
letters9 = sp.symbols(' '.join(MIXED + CLASS), real=True)
for n in (2, 3):
    gs = sp.symbols(f'g1:{n + 1}', real=True)
    monos = exponents(9 + n, 3)
    P = interpolate(monos, letters9 + gs, lambda pt, p: modp_trace_power(
        modp_generator(n, dict(zip(MIXED + CLASS, pt[:9])), pt[9:], p), 3, p), nprimes=2)
    check(f"N = {n}: the third sum is zero for every bond of the nine letters at every rate", P == 0)


def fifth_basis(n):
    """Mixed letters with degree 2 or 4 (a half turn on every site flips them alone), at least one rate."""
    return [m + c + r for dm in (2, 4) for dr in range(1, 6 - dm) for m in exponents(4, dm)
            for c in exponents(5, 5 - dm - dr) for r in exponents(n, dr)]


for n in ((2, 3, 4, 5) if DEEP else (2,)):
    gs = sp.symbols(f'g1:{n + 1}', real=True)
    P = interpolate(fifth_basis(n), letters9 + gs, lambda pt, p: modp_trace_power(
        modp_generator(n, dict(zip(MIXED + CLASS, pt[:9])), pt[9:], p), 5, p), nprimes=2, surplus=25)
    check(f"N = {n}: the fifth sum is the closed form as a polynomial in the nine coefficients and the rates",
          sp.expand(P - fifth_closed_form(n, dict(zip(MIXED + CLASS, letters9)), gs)) == 0)
rng = np.random.default_rng(11)
for n in (3, 4, 5, 6):
    c = {l: int(x) for l, x in zip(MIXED + CLASS, rng.integers(-5, 6, 9))}
    g = [int(x) for x in rng.integers(1, 50, n)]
    p = PRIMES[1]
    check(f"N = {n}, coefficients {tuple(c.values())}, rates {tuple(g)}: the fifth sum equals the closed form "
          f"modulo {p}",
          modp_trace_power(modp_generator(n, c, g, p), 5, p) == int(fifth_closed_form(n, ints(c), g)) % p)
c0 = {**{l: 0 for l in MIXED + CLASS}, 'XZ': 1, 'ZX': 1, 'XY': -3, 'YX': -3}
M = modp_generator(2, c0, [8, 8], PRIMES[0])
check("XZ + ZX − 3XY − 3YX, N = 2, rates (8, 8): the anisotropy cancels the fifth sum (the closed form gives zero), "
      "and the seventh is not zero", fifth_closed_form(2, ints(c0), [8, 8]) == 0
      and modp_trace_power(M, 5, PRIMES[0]) == 0 and modp_trace_power(M, 7, PRIMES[0]) != 0)

check("XZ alone, N = 2, rates (1, 0): a site without cooling empties S, the fifth sum is exactly zero, and the "
      "spectrum pairs (characteristic polynomial)",
      fifth_closed_form(2, ints({**{l: 0 for l in MIXED + CLASS}, 'XZ': 1}), [1, 0]) == 0
      and trace_is_exactly_zero(2, {'XZ': 1}, [1, 0], 5) and pairs_exactly(2, {'XZ': 1}, [1, 0]))
check("XZ alone, N = 3, rates (1, 0, 1): the spectrum pairs (characteristic polynomial)",
      pairs_exactly(3, {'XZ': 1}, [1, 0, 1]))
check("control: XZ alone, N = 3, rates (1, 1, 0): it does not", not pairs_exactly(3, {'XZ': 1}, [1, 1, 0]))

print("G  rings")
for n in ((3, 4, 5) if DEEP else (3, 4)):
    for bond, name in (bonds_E[1:] if n < 5 else bonds_E[1:2]):
        P, gs = rate_polynomial(n, bond, 7, ring=True)
        check(f"ring of {n}, {name}: the seventh sum is zero at every rate", P == 0)
c = {**{l: 0 for l in MIXED + CLASS}, 'XZ': 3, 'YZ': -1, 'ZX': 2, 'ZY': 4, 'XX': 1, 'XY': 2, 'YY': -1, 'ZZ': 1}
check("ring of 4, XZ, YZ, ZX, ZY with XX, XY, YY, ZZ, rates (3, 7, 2, 5): the fifth sum equals the ring's closed form "
      "(no end terms) modulo a prime",
      modp_trace_power(modp_generator(4, c, [3, 7, 2, 5], PRIMES[0], ring=True), 5, PRIMES[0])
      == int(fifth_closed_form(4, ints(c), [3, 7, 2, 5], ring=True)) % PRIMES[0])
check("ring of 4, XZ − ZX at equal rates (μ = −ν, every square zero): the fifth and seventh sums are exactly zero, "
      "the ninth is not", trace_is_exactly_zero(4, {'XZ': 1, 'ZX': -1}, [1, 1, 1, 1], 5, ring=True)
      and trace_is_exactly_zero(4, {'XZ': 1, 'ZX': -1}, [1, 1, 1, 1], 7, ring=True)
      and modp_trace_power(modp_generator(4, {'XZ': 1, 'ZX': -1}, [1, 1, 1, 1], PRIMES[0], ring=True), 9, PRIMES[0])
      != 0)
c = {**{l: 0 for l in MIXED + CLASS}, 'XZ': 3, 'YZ': -1, 'ZX': 2, 'ZY': 4, 'XX': 1}
check("triangle, 3XZ − YZ + 2ZX + 4ZY + XX, rates (3, 7, 2): its fifth sum is not the closed form without end terms "
      "(the windows wrap)", modp_trace_power(modp_generator(3, c, [3, 7, 2], PRIMES[0], ring=True), 5, PRIMES[0])
      != int(fifth_closed_form(3, ints(c), [3, 7, 2], ring=True)) % PRIMES[0])
check("ring of 4, XX + 2XY + ZZ, rates (5, 3, 9, 10) on the surface Q(γ₁, γ₃) = Q(γ₂, γ₄) = 575 that no reflection "
      "keeps: the ninth sum is exactly zero and the eleventh is not",
      trace_is_exactly_zero(4, bonds_E[0][0], [5, 3, 9, 10], 9, ring=True)
      and modp_trace_power(modp_generator(4, bonds_E[0][0], [5, 3, 9, 10], PRIMES[0], ring=True), 11, PRIMES[0]) != 0)

def wrap_invariants(bond):
    """(Im z³, Im z·|w|²), the two products that close round the triangle and change under H ↦ H̄."""
    J = {l: sp.Integer(bond.get(l, 0)) for l in CLASS}
    z = (J['XX'] + J['YY']) / 2 - sp.I * (J['XY'] - J['YX']) / 2
    w2 = ((J['XX'] - J['YY']) / 2) ** 2 + ((J['XY'] + J['YX']) / 2) ** 2
    return sp.im(sp.expand(z ** 3)), sp.im(z) * w2


(u1, v1), (u2, v2) = wrap_invariants(bonds_E[1][0]), wrap_invariants(bonds_E[2][0])
check("triangle: the two bonds whose seventh is zero above have independent (Im z³, Im z·|w|²), so both numbers of "
      "the triangle's seventh, α·Im z³ + β·Im z·|w|², are zero", u1 * v2 - u2 * v1 != 0,
      f'({u1}, {v1}) and ({u2}, {v2})')
rnd = random.Random(7)
for n in (3, 4, 5, 6):
    bond = {l: rnd.randint(-9, 9) for l in CLASS}
    rates = [rnd.randint(1, 60) for _ in range(n)]
    M, p = modp_generator(n, bond, rates, PRIMES[0], ring=True), PRIMES[0]
    tr = [modp_trace_power(M, k, p) for k in (3, 5, 7, 9)]
    check(f"ring of {n}, random class bond {bond}, rates {rates}: the third, fifth and seventh sums are zero and the "
          "ninth is not, modulo a prime", tr[:3] == [0, 0, 0] and tr[3] != 0, f'{tr}')
for n in (3, 4, 5):
    bonds = [{l: rnd.randint(-9, 9) for l in CLASS} for _ in range(n)]
    rates = [rnd.randint(1, 60) for _ in range(n)]
    M, p = modp_generator(n, bonds, rates, PRIMES[0], ring=True), PRIMES[0]
    t5, t7 = modp_trace_power(M, 5, p), modp_trace_power(M, 7, p)
    if n == 3:
        check("control, triangle with a different bond on each edge: no rotation, and the fifth sum is not zero "
              "(three factors close round it)", t5 != 0, f'{t5}')
    else:
        check(f"control, ring of {n} with a different bond on each edge: no rotation, the seventh sum is not zero, "
              "and the fifth still is (four factors and an odd power of the couplings leave no rate)",
              t5 == 0 and t7 != 0, f'{t5}, {t7}')

tri = []
for bond, name in bonds_E[:2]:
    J = {l: sp.Integer(bond.get(l, 0)) for l in CLASS}
    inv = J['ZZ'] * (J['XX'] + J['YY']) / 2 * (J['XY'] - J['YX']) / 2 * (
        ((J['XX'] - J['YY']) / 2) ** 2 + ((J['XY'] + J['YX']) / 2) ** 2)
    P, (g1, g2, g3, g4) = rate_polynomial(4, bond, 9, ring=True)
    target = -9 * 2 ** 12 * 2 ** 9 * inv * (g1 - g3) * (g2 - g4) * (
        5 * g1 ** 2 + g1 * g3 + 5 * g3 ** 2 - 5 * g2 ** 2 - g2 * g4 - 5 * g4 ** 2)
    check(f"ring of 4, {name}: Tr((L + c/2)^9) = −9·2¹²·Δ·s·a·|w|²·(γ₁ − γ₃)(γ₂ − γ₄)(Q(γ₁, γ₃) − Q(γ₂, γ₄)), "
          "Q(x, y) = 5x² + xy + 5y²", sp.expand(P - target) == 0)
    P, (g1, g2, g3) = rate_polynomial(3, bond, 9, ring=True)
    q = sp.cancel(P / ((g1 - g2) * (g1 - g3) * (g2 - g3) * (g1 + g2 + g3)))
    check(f"triangle, {name}: the ninth sum is a nonzero constant times (γ₁ − γ₂)(γ₁ − γ₃)(γ₂ − γ₃)(γ₁ + γ₂ + γ₃)",
          q.is_number and q != 0, f'constant {q / 2 ** 9}')
    tri.append((q, inv))
check("triangle: the two constants are not in the ratio of Δ·s·a·|w|² (the triangle is not bipartite)",
      tri[0][0] * tri[1][1] != tri[1][0] * tri[0][1], f'{tri[0][0] / tri[1][0]} against {tri[0][1] / tri[1][1]}')

print("H  the ninth cumulant on the chain in closed form")


def invariants(bond):
    """(Δ, s, a, |w|²) of a class bond, J = s·1 + a·ε + T, w = T₁₁ + i·T₁₂."""
    J = {l: sp.Integer(bond.get(l, 0)) for l in CLASS}
    return (J['ZZ'], (J['XX'] + J['YY']) / 2, (J['XY'] - J['YX']) / 2,
            ((J['XX'] - J['YY']) / 2) ** 2 + ((J['XY'] + J['YX']) / 2) ** 2)


def ninth_cumulant(g, D, s_, a, w2):
    """κ₉ of the spectrum of L + c/2 (p₉ = 4^N·κ₉ + 36·p₇·p₂/4^N): 144·Δ·s·a·|w|²·[Q + E_l − E_r + 16(e_l − e_r)]."""
    N = len(g)
    Q = sum(2 * g[l] * g[l + 1] * (g[l] ** 2 - g[l + 1] ** 2) for l in range(N - 1))
    Q += sum(g[l] * g[l + 2] * (11 * g[l + 1] * (g[l] - g[l + 2]) - 7 * (g[l] ** 2 - g[l + 2] ** 2)) for l in range(N - 2))

    def E(x, y, z):
        return x * y * (38 * x ** 2 - 3 * x * y + 40 * y ** 2) - 34 * x ** 2 * y * z

    def e(x, y, z):
        return (-6 * (s_ ** 2 + a ** 2 + w2) * x ** 2 + (16 * s_ ** 2 + 40 * a ** 2 - 40 * w2) * y ** 2
                - (111 * D ** 2 + 226 * s_ ** 2 + 374 * a ** 2 + 66 * w2) * x * y
                + (-23 * D ** 2 - 33 * s_ ** 2 - 181 * a ** 2 + 61 * w2) * y * z
                + (-7 * D ** 2 - 196 * s_ ** 2 + 28 * a ** 2 + 28 * w2) * x * z)
    l_, r_ = (g[0], g[1], g[2]), (g[N - 1], g[N - 2], g[N - 3])
    return 144 * D * s_ * a * w2 * (Q + E(*l_) - E(*r_) + 16 * (e(*l_) - e(*r_)))


def rate_polynomial_in(n, bond, k, degrees):
    """Tr((2L + c)^k) as a polynomial in the rates with only the given total degrees (the surplus points must fit)."""
    gs = sp.symbols(f'g1:{n + 1}', real=True)
    monos = [e for d_ in degrees for e in exponents(n, d_)]
    return interpolate(monos, gs, lambda pt, p: modp_trace_power(modp_generator(n, bond, pt, p), k, p)), gs


# four bonds whose (Δ², s², a², |w|²) are independent fix the quadratic part; two more check it
H_bonds = bonds_E + [({'XX': 2, 'XY': 1, 'YX': -1, 'YY': 1, 'ZZ': 1}, '2XX + XY − YX + YY + ZZ'),
                     ({'XX': 1, 'XY': 3, 'YX': 1, 'YY': -2, 'ZZ': 2}, 'XX + 3XY + YX − 2YY + 2ZZ'),
                     ({'XX': 4, 'XY': 1, 'YY': 1, 'ZZ': -1}, '4XX + XY + YY − ZZ')]
for n, bonds in ((3, H_bonds), (4, H_bonds), (5, H_bonds[:4] if DEEP else []), (6, H_bonds[:1] if DEEP else [])):
    for bond, name in bonds:
        P9, gs = rate_polynomial_in(n, bond, 9, (2, 4))
        P7, _ = rate_polynomial_in(n, bond, 7, (2,))
        P2, _ = rate_polynomial_in(n, bond, 2, (0, 1, 2))
        k9 = ninth_cumulant(list(gs), *invariants(bond))
        check(f"N = {n}, {name}: Tr((2L + c)^9) − 36·Tr((2L + c)^7)·Tr((2L + c)^2)/4^N = 2^9·4^N·κ₉, κ₉ the closed form, "
              "as a polynomial in the rates", sp.expand(P9 - sp.Rational(36, 4 ** n) * P7 * P2 - 512 * 4 ** n * k9) == 0)
rnd = random.Random(9)
for n in ((4, 5, 6, 7) if DEEP else (4, 5, 6)):
    while True:   # a draw where every part is there and the seventh is not silent, so the check can fail
        bond = {l: rnd.randint(-9, 9) for l in CLASS}
        rates = [rnd.randint(1, 99) for _ in range(n)]
        if 0 not in invariants(bond) and rates[0] * rates[1] != rates[-1] * rates[-2]:
            break
    p = PRIMES[1]
    M = modp_generator(n, bond, rates, p).astype(np.float64)
    R_, tr = M, {}
    for k in range(2, 10):
        R_ = (R_ @ M) % p
        if k in (2, 7, 9):
            tr[k] = int(np.trace(R_)) % p
    k9 = sp.Rational(ninth_cumulant(rates, *invariants(bond)))
    pred = (512 * pow(4, n, p) * int(k9.p) * pow(int(k9.q), p - 2, p)
            + 36 * tr[7] * tr[2] * pow(pow(4, n, p), p - 2, p)) % p
    check(f"N = {n}, random class bond {bond}, rates {rates}: the ninth sum is the closed form, modulo a prime",
          pred == tr[9] and tr[9] != 0, f'{tr[9]} against {pred}')
D_, s_, a_, w_, x_, z_ = sp.symbols('Delta s a w2 x z')
check("N = 3 with no cooling in the middle, every bond: κ₉ = −144·Δ·s·a·|w|²·(γ₁² − γ₃²)·(7γ₁γ₃ + 96(s² + a² + |w|²)), "
      "so the ninth breaks wherever γ₁ ≠ γ₃",
      sp.expand(ninth_cumulant([x_, 0, z_], D_, s_, a_, w_)
                + 144 * D_ * s_ * a_ * w_ * (x_ ** 2 - z_ ** 2) * (7 * x_ * z_ + 96 * (s_ ** 2 + a_ ** 2 + w_))) == 0)
for rates, bond, name in [([1, 6, 3, 2], {'XX': -6, 'XY': -2, 'YX': 1, 'YY': 1, 'ZZ': -2}, '−6XX − 2XY + YX + YY − 2ZZ'),
                          ([1, 8, 1, 8], {'XX': -6, 'XY': -6, 'YX': 2, 'YY': 2, 'ZZ': 3}, '−6XX − 6XY + 2YX + 2YY + 3ZZ')]:
    check(f"N = 4, {name}, rates {tuple(rates)} (equal end products, not kept by the reversal), where the closed form "
          "puts κ₉ = 0: the seventh and ninth sums are exactly zero and the eleventh is not",
          ninth_cumulant(rates, *invariants(bond)) == 0
          and trace_is_exactly_zero(4, bond, rates, 7) and trace_is_exactly_zero(4, bond, rates, 9)
          and modp_trace_power(modp_generator(4, bond, rates, PRIMES[0]), 11, PRIMES[0]) != 0)

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILED")
    sys.exit(1)
print("ALL PASS")
