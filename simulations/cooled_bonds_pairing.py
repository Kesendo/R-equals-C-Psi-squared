"""Which two-site bonds pair under cooling alone: exact checks with controls.

experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md. Every site carries cooling sigma^- = |0><1| at its own rate and nothing
else; H is a chain of one bond repeated, J_ab·sigma_a(x)sigma_b over a, b in {X, Y} plus Δ·ZZ, or a bond with a mixed
letter (XZ, ZX). Pairing means the spectrum of L is symmetric about −c/2, c = Σ rates. Rows print PASS or FAIL; the
script exits nonzero on any FAIL. With --map it first prints the float survey the document quotes (a reading, not a
gate): random integer couplings over every subset of the nine two-letter bonds, N = 2, 3, 4.

  A  The turn lemma: for every real 2x2 J there is a rotation R with R·Jᵀ·Rᵀ = F·J·F, F = diag(1, −1) (Jᵀ and FJF share
     trace, antisymmetric part and the length of the traceless symmetric part). Checked on six J.
  B  The carrier at rates the chain's reversal P keeps: E = (Ad X^N ∘ Ad(U·P))⁻¹ ∘ S_IZ, S_IZ: rho -> rho·Z^N, U the
     uniform turn about Z realizing R, gives E·L·E⁻¹ = −L† − c exactly (F137's Ad_V∘S_XY with V = (U·P)†); controls:
     rates P does not keep.
  C  The Theodolite's elements (exact, read off the letters and verified): X/Y letters alone at any rates (a staggered
     A·rhoᵀ·B); with ZZ, a real or an imaginary X/Y part (A·rho·B), J symmetric (a = 0) or traceless (s = 0) at
     unequal rates, and any J on one bond (turned); silent where s, a and w are all nonzero at rates the reversal
     does not keep. J = s·1 + a·ε + T, ε = [[0, 1], [−1, 0]], T traceless symmetric, w = T₁₁ + i·T₁₂.
  D  Breaks, exact: an odd centred power sum Tr((L − tr L/dim)^k) that is not zero rules out every element; the first
     such k is asserted: 7 where s, a, w are all nonzero (XX + 2XY + ZZ at (2, 1, 1), (1, 2, 3) and, N = 4,
     (1, 2, 3, 4); 2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3)), 5 for a mixed letter (XZ; ZX + XX). Control: zero up to
     k = 7 at rates the reversal keeps.
"""
from __future__ import annotations

import itertools
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
check("s, a, w all nonzero with ZZ on one bond, XX + 2XY + ZZ, N = 2, rates (1, 2): a turned sandwich reflects",
      letters(setting(2, {'XX': 1, 'XY': 2, 'ZZ': 1}, [1, 2]), 'refl') is not None)
check("s, a, w all nonzero with ZZ, N = 3, rates the reversal keeps (1, 2, 1): the letters are silent (no site map "
      "carries the reversal)", letters(setting(3, {'XX': 1, 'XY': 2, 'ZZ': 1}, [1, 2, 1]), 'refl') is None)
check("s, a, w all nonzero with ZZ, N = 3: the letters are silent (XX + 2XY + ZZ at (2, 1, 1), "
      "2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3))",
      letters(setting(3, {'XX': 1, 'XY': 2, 'ZZ': 1}, [2, 1, 1]), 'refl') is None
      and letters(setting(3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': 2, 'ZZ': 1}, [1, 2, 3]), 'refl') is None)

print("D  breaks, exact")
mixed = {'XX': 1, 'XY': 2, 'ZZ': 1}
for n, bond, rates, name, first in [
        (3, mixed, [2, 1, 1], 'XX + 2XY + ZZ, rates (2, 1, 1)', 7),
        (3, mixed, [1, 2, 3], 'XX + 2XY + ZZ, rates (1, 2, 3)', 7),
        (4, mixed, [1, 2, 3, 4], 'XX + 2XY + ZZ, N = 4, rates (1, 2, 3, 4)', 7),
        (3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': 2, 'ZZ': 1}, [1, 2, 3], '2XX + 3XY − YX + 2YY + ZZ, rates (1, 2, 3)', 7),
        (3, {'XZ': 1}, [1, 1, 1], 'XZ, rates (1, 1, 1)', 5),
        (3, {'ZX': 1, 'XX': 1}, [1, 2, 3], 'ZX + XX, rates (1, 2, 3)', 5)]:
    k = first_odd_trace(n, bond, rates)
    check(f"{name}: the first odd centred power sum that is not zero is the {first}th", k == first, f'k = {k}')
check("control: XX + 2XY + ZZ at rates (1, 2, 1) has every odd power sum zero up to k = 7",
      first_odd_trace(3, mixed, [1, 2, 1]) is None)

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILED")
    sys.exit(1)
print("ALL PASS")
