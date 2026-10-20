"""What the letter cube, applied site by site (Device B, a working name), found in the repo: exact checks with controls.

experiments/THE_CUBE_SITE_BY_SITE.md. Exact arithmetic (sympy over Q(i), or with sqrt where a point needs it).
Every row prints PASS or FAIL and the script exits nonzero on any FAIL. Conventions of framework.cube_tables: one
site's operators on the cells |i><j|, row-major; sigma^- = |0><1| cools, sigma^+ = |1><0| heats.

  A  The weight of the site-local Pi-shaped element is the one-site steady state. A site with cooling d, heating u and
     dephasing at rate p along any axis a·X + b·Y of the XY plane (a² + b² = 1) has the steady-state ratio
     q² = (u + p)/(d + p); the site map rho -> w·Y·rho^T·X·w, w = diag(1, q), solves the site equation exactly there.
     The identity row: the left eigenvector of the site generator at its middle rate is Z + κ·I with
     κ = (u − d)/(u + d + 2p) = (q² − 1)/(q² + 1), and with damping it spans that eigenspace.
  B  On whole chains: the matched element reflects L onto −L − c exactly (Heisenberg, cooling plus X dephasing; the
     flux-pi/2 triangle at Δ = 7/10 with three unaligned XY axes and one common q² = 1/4); a ratio that differs across
     a hopping bond leaves no product Pi-reflection (Device B's pre-check, exact, and the identity-row commutator
     ‖[H, ⊗(Z_l + κ_l)]‖_F = 8√803/85 on the Heisenberg chain with a temperature gradient) while the reflection
     −L† − c persists; raising one heating rate breaks the triangle's spectrum (Tr((L − μ)^5) = 192).
  C  The golden family. At c = 0 (windows XZY + YZX) one site map on every site, rho -> (X + Y)·rho·diag(1, −i),
     solves the site and window equations (every open chain, every ring); the bond XX + XY + YX has the uniform map
     rho -> X·(g·rho)^T·X with g = φX + Y, of order four, and its identity column is forced into the golden frame
     (α² − αβ − β² = 0 for g = αX + βY); it solves XIX + XIY + YIX bond by bond. Controls: c = 3/2, the weights
     (0, 1, 2), g = 2X + Y.
  D  Heisenberg + DM on a chain has no period-2 product of any invertible site maps (pure DM has one); Device B reads
     it as TURNED: the canonical Pi turned about Z by arctan(D/J) from site to site reflects it exactly, and twice that
     angle does not.
  E  The boundary of A, computed and not proven: with cooling alone (ratio 0) an isolated site has a site-local
     Pi-reflection, two sites joined by a Heisenberg bond have none.
  F  Where fields break M_beta: on F137's triangle with fields (3/10, −1/5, 1/2) at ratio 1/4, Δ = 0 and flux ±pi/2
     the spin flip with the whole thermal weight on one side, rho -> π·Z₀Z₂·X^N·rho·X^N·Z₁, is a site-local
     Pi-reflection; controls Δ = 7/10 and flux 0.
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding='utf-8')

from framework import cube_tables as ct  # noqa: E402
from framework.diagnostics import cube_reflector as cr  # noqa: E402

R = sp.Rational
X, Y, Z = ct.PAULI['X'], ct.PAULI['Y'], ct.PAULI['Z']
TRANSPOSE = sp.zeros(4, 4)
for _r, _c in ((0, 0), (1, 2), (2, 1), (3, 3)):
    TRANSPOSE[_r, _c] = 1
FAILURES = []


def check(name, ok, detail=''):
    print(f"  {'PASS' if ok else 'FAIL'}: {name}" + (f"  [{detail}]" if detail else ''))
    if not ok:
        FAILURES.append(name)


def matched_map(q):
    """rho -> w·Y·rho^T·X·w with w = diag(1, q), as a one-site cell map."""
    w = sp.diag(1, q)
    return ct.lr(w * Y, X * w) * TRANSPOSE


def traceless(M):
    return M - (M.trace() / M.shape[0]) * sp.eye(M.shape[0])


def full_generator(setting):
    """Independent route: the row-major Lindblad generator from the Lindblad formula, then per-site order."""
    n = setting.n
    d = 2 ** n

    def on(op_by_site):
        return sp.kronecker_product(*[op_by_site.get(l, sp.eye(2)) for l in range(n)])

    H = sp.zeros(d, d)
    for letters, sites, coef in setting.terms:
        H += coef * on({s: ct.PAULI[c] for s, c in zip(sites, letters)})
    Id = sp.eye(d)
    L = -sp.I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
    for sites, A, rate in setting.jumps:
        J = on({sites[0]: sp.Matrix(A)})
        JdJ = J.H * J
        L += rate * (sp.kronecker_product(J, J.conjugate())
                     - R(1, 2) * (sp.kronecker_product(JdJ, Id) + sp.kronecker_product(Id, JdJ.T)))
    P = ct._per_site_perm(n)
    return (P * L * P.T).applyfunc(sp.expand)


def reflects_pi(setting, E1s):
    E = sp.kronecker_product(*E1s)
    L = full_generator(setting)
    c = cr.constant_c(setting)
    return (E * L + (L + c * sp.eye(L.shape[0])) * E).applyfunc(sp.simplify).is_zero_matrix


def odd_sum(setting, k):
    L = full_generator(setting)
    n2 = L.shape[0]
    M = DomainMatrix.from_Matrix((L - (L.trace() / n2) * sp.eye(n2)).applyfunc(sp.expand)).convert_to(sp.QQ_I)
    P = M
    for _ in range(k - 1):
        P = P * M
    return sp.nsimplify(sum(P.to_Matrix()[i, i] for i in range(n2)))


def heisenberg(n):
    return [(p, (a, a + 1), 1) for a in range(n - 1) for p in ('XX', 'YY', 'ZZ')]


def bond(a, b, J, D, Dz):
    t = [('XX', (a, b), J), ('YY', (a, b), J)] if J else []
    t += [('XY', (a, b), D), ('YX', (a, b), -D)] if D else []
    return t + ([('ZZ', (a, b), Dz)] if Dz else [])


print("A  the weight is the one-site steady state")
d, u, p, q = sp.symbols('d u p q', positive=True)
tau = sp.symbols('tau', real=True)
a, b = (1 - tau ** 2) / (1 + tau ** 2), 2 * tau / (1 + tau ** 2)      # every XY axis but -X, rationally
axis = a * X + b * Y
site = ct.dissipator(ct.LOWER, d) + ct.dissipator(ct.RAISE, u) + ct.dissipator(axis, p)
L0 = traceless(site)
residual = (matched_map(q) * L0 + L0 * matched_map(q)).applyfunc(lambda e: sp.factor(sp.simplify(e)))
entries = {e for e in residual if e != 0}
check("site equation of w·Y·rho^T·X·w for any XY axis a·X + b·Y: every residual entry is ∝ q²(d + p) − (u + p)",
      bool(entries) and all(sp.simplify(e / (q ** 2 * (d + p) - (u + p))).free_symbols <= {sp.I} or
                            sp.simplify(e / (q ** 2 * (d + p) - (u + p))).is_number for e in entries),
      f"{len(entries)} nonzero entries")
check("so it vanishes exactly at q² = (u + p)/(d + p)",
      residual.subs(q, sp.sqrt((u + p) / (d + p))).applyfunc(sp.simplify).is_zero_matrix)
steady = sp.Matrix([1, 0, 0, sp.symbols('r')])
r_ss = sp.solve((site * steady)[0], sp.symbols('r'))
check("q² is the population ratio p₁/p₀ of the one-site steady state", r_ss == [(u + p) / (d + p)], f"{r_ss}")
kappa = sp.symbols('kappa')
left = (site.H * sp.Matrix([1 + kappa, 0, 0, kappa - 1]) - (site.trace() / 2) * sp.Matrix([1 + kappa, 0, 0, kappa - 1]))
k_sol = sp.solve([sp.simplify(e) for e in left if sp.simplify(e) != 0], kappa)
k_expect = (u - d) / (u + d + 2 * p)
check("identity row: Z + κ·I is a left eigenvector at tr L/2 with κ = (u − d)/(u + d + 2p) = (q² − 1)/(q² + 1)",
      bool(k_sol) and sp.simplify(k_sol[kappa] - k_expect) == 0
      and sp.simplify(k_expect - ((u + p) / (d + p) - 1) / ((u + p) / (d + p) + 1)) == 0, f"{k_sol}")
check("... and with damping (d + u > 0) it spans that eigenspace: its dimension is 1",
      len((site.H - (site.trace() / 2) * sp.eye(4)).nullspace()) == 1)

print("B  on whole chains")
s = ct.Setting(3, heisenberg(3), [((l,), ct.LOWER, 1) for l in range(3)] + [((l,), X, R(1, 2)) for l in range(3)])
check("Heisenberg N = 3, cooling 1 and X dephasing 1/2: the matched map (q² = 1/3) is a Pi-reflection",
      reflects_pi(s, [matched_map(sp.sqrt(R(1, 3)))] * 3))
check("control: q = 1 is not", not reflects_pi(s, [matched_map(1)] * 3))


def triangle(heat0):
    terms = bond(0, 1, 1, 0, R(7, 10)) + bond(1, 2, 2, 0, R(7, 10)) + bond(2, 0, 0, 3, R(7, 10))
    axes = [(3 * X + 4 * Y) / 5, X, (5 * X - 12 * Y) / 13]
    heat, deph = [heat0, R(1, 16), 0], [R(1, 5), R(1, 4), R(1, 3)]
    jumps = [((l,), ct.LOWER, 1) for l in range(3)] + [((l,), ct.RAISE, heat[l]) for l in range(3) if heat[l]]
    return ct.Setting(3, terms, jumps + [((l,), axes[l], deph[l]) for l in range(3)])


t = triangle(R(1, 10))
check("triangle, hoppings 1, 2, 3i (flux pi/2), Δ = 7/10, axes (3X+4Y)/5, X, (5X−12Y)/13, all q² = 1/4: "
      "the matched map is a Pi-reflection", reflects_pi(t, [matched_map(R(1, 2))] * 3))
check("control: q = 1 is not", not reflects_pi(t, [matched_map(1)] * 3))
check("raising the heating on site 0 by 1/10 breaks the spectrum: Tr((L − μ)^5) = 192",
      odd_sum(triangle(R(1, 5)), 5) == 192, f"{odd_sum(triangle(R(1, 5)), 5)}")
down, up = [R(1, 2), R(1), R(3, 2)], [R(1, 8), R(1, 3), R(1, 5)]
g = ct.Setting(3, heisenberg(3), [((l,), ct.LOWER, x) for l, x in enumerate(down)]
               + [((l,), ct.RAISE, x) for l, x in enumerate(up)])
check("temperature gradient: Device B's pre-check proves no product Pi-reflection", cr.precheck(g, 'pi')[0] == 'NONE')
check("... while the reflection −L† − c is found as cube moves", cr.cube_search(g, 'refl')[0] == 'FOUND')
kap = [(x - y) / (x + y) for y, x in zip(down, up)]
Hn = sp.zeros(8, 8)
for letters, sites, coef in g.terms:
    Hn += coef * sp.kronecker_product(*[ct.PAULI[letters[sites.index(l)]] if l in sites else sp.eye(2) for l in range(3)])
G = sp.kronecker_product(*[Z + k * sp.eye(2) for k in kap])
norm = sp.sqrt(sum(abs(v) ** 2 for v in (Hn * G - G * Hn)))
check("the identity-row commutator ‖[H, ⊗(Z_l + κ_l)]‖_F = 8√803/85", sp.simplify(norm - 8 * sp.sqrt(803) / 85) == 0,
      f"{sp.nsimplify(norm)}")

print("C  the golden family")
zsite = traceless(ct.site_table(ct.Setting(1, [], [((0,), Z, 1)]), 0))
router = ct.lr(X + Y, sp.diag(1, -sp.I))
check("c = 0: (X + Y)·rho·diag(1, −i) solves the Z-dephased site equation",
      (router * zsite + zsite * router).applyfunc(sp.expand).is_zero_matrix)
for w, expect in (((0, 1, 1), True), ((R(3, 2), 1, 1), False), ((0, 1, 2), False)):
    LW = ct.support_tables(ct.Setting(3, [('XZX', (0, 1, 2), w[0]), ('XZY', (0, 1, 2), w[1]),
                                          ('YZX', (0, 1, 2), w[2])]))[(0, 1, 2)]
    E3 = sp.kronecker_product(router, router, router)
    holds = (E3 * LW + LW * E3).applyfunc(sp.expand).is_zero_matrix
    check(f"window weights {w}: the uniform map {'solves' if expect else 'does not solve'} the window equation",
          holds == expect)
phi = (1 + sp.sqrt(5)) / 2
LB = ct.support_tables(ct.Setting(2, [('XX', (0, 1), 1), ('XY', (0, 1), 1), ('YX', (0, 1), 1)]))[(0, 1)]
for name, coef, expect in (('φ', phi, True), ('2', 2, False)):
    gg = coef * X + Y
    E1 = ct.lr(X, gg.T * X) * TRANSPOSE
    ok_site = (E1 * zsite + zsite * E1).applyfunc(sp.simplify).is_zero_matrix
    E2 = sp.kronecker_product(E1, E1)
    ok_bond = (E2 * LB + LB * E2).applyfunc(sp.simplify).is_zero_matrix
    check(f"bond XX + XY + YX: X·(g·rho)^T·X with g = {name}·X + Y {'is' if expect else 'is not'} a uniform Pi-map",
          ok_site and (ok_bond == expect))

print("D  Heisenberg + DM is turned, not alternating")


def period2_exists(setting, rel):
    fam = [cr.site_family(ct.site_table(setting, l), rel) for l in range(2)]
    E, gens, eqs = [], [], []
    for k, f in enumerate(fam):
        x = sp.symbols(f'x{k}_0:{len(f)}')
        tt = sp.Symbol(f't{k}')
        Ek = sum((xx * B for xx, B in zip(x, f)), sp.zeros(4, 4))
        E.append(Ek)
        gens += list(x) + [tt]
        eqs.append(tt * Ek.det() - 1)
    for S, LS in ct.support_tables(setting).items():
        ES = sp.kronecker_product(*[E[l % 2] for l in S])
        eqs += [e for e in (ES * LS - cr.target(LS, rel) * ES).applyfunc(sp.expand) if e != 0]
    return sp.groebner(eqs, *gens, order='grevlex', domain=sp.QQ_I).exprs != [1]


zj = [((l,), Z, x) for l, x in enumerate([R(1, 5), R(2, 7), R(1, 3)])]
heis_dm = ct.Setting(3, [t for a in range(2) for t in bond(a, a + 1, 1, R(1, 2), 1)], zj)
pure_dm = ct.Setting(3, [t for a in range(2) for t in bond(a, a + 1, 0, 1, 0)], zj)
check("Heisenberg + DM (D/J = 1/2), Z dephasing: no period-2 product of any invertible site maps (Pi-reflection)",
      not period2_exists(heis_dm, 'pi'))
check("control: pure DM has one", period2_exists(pure_dm, 'pi'))
cell = cr.Cell(jumps=[(Z, 1)], windows=[('XX', (0, 1), 1), ('YY', (0, 1), 1), ('ZZ', (0, 1), 1),
                                        ('XY', (0, 1), R(1, 2)), ('YX', (0, 1), -R(1, 2))])
check("Device B reads Heisenberg + DM as TURNED", cr.every_n(cell, 'pi') == 'TURNED')
PI_P = sp.zeros(4, 4)
for src, dst, ph in ((0, 1, 1), (1, 0, 1), (2, 3, sp.I), (3, 2, sp.I)):     # I -> X, X -> I, Y -> iZ, Z -> iY
    PI_P[dst, src] = ph
PI_CELLS = ct.to_cells(PI_P)
L_hd, c_hd = full_generator(heis_dm), cr.constant_c(heis_dm)


def turned_pi_reflects(z):
    maps = []
    for l in range(3):
        Ad = ct.lr(sp.diag(1, z ** l), sp.diag(1, sp.conjugate(z) ** l))
        Adinv = ct.lr(sp.diag(1, z ** -l), sp.diag(1, sp.conjugate(z) ** -l))
        maps.append((Ad * PI_CELLS * Adinv).applyfunc(sp.radsimp))
    E = sp.kronecker_product(*maps)
    res = (E * L_hd + (L_hd + c_hd * sp.eye(64)) * E).applyfunc(lambda e: sp.radsimp(sp.expand(e)))
    return all(sp.simplify(v) == 0 for v in res)


e1 = (2 + sp.I) / sp.sqrt(5)                                    # e^{i·arctan(1/2)}
check("the canonical Pi turned about Z by arctan(D/J) per site reflects the Heisenberg + DM chain exactly",
      turned_pi_reflects(e1))
check("controls: turned by −arctan(D/J) or by 2·arctan(D/J) it does not",
      not turned_pi_reflects(sp.conjugate(e1)) and not turned_pi_reflects(sp.expand(e1 ** 2)))

print("C' the golden bond inside a three-site window")
xix = ct.Setting(4, [(lab, (w, w + 1, w + 2), 1) for w in range(2) for lab in ('XIX', 'XIY', 'YIX')],
                 [((l,), Z, 1) for l in range(4)])
gmap = ct.lr(X, (phi * X + Y).T * X) * TRANSPOSE
tabs = ct.support_tables(xix)
ok = all((sp.kronecker_product(*[gmap] * len(S)) * LS + LS * sp.kronecker_product(*[gmap] * len(S)))
         .applyfunc(sp.simplify).is_zero_matrix for S, LS in tabs.items())
check("XIX + XIY + YIX: its support sets are the next-nearest bonds, and the golden map solves each of them",
      set(tabs) == {(0, 2), (1, 3)} and ok)
al, be = sp.symbols('alpha beta')
gab = al * X + be * Y
hb = sp.kronecker_product(X, X) + sp.kronecker_product(X, Y) + sp.kronecker_product(Y, X)
cm = [e for e in (hb * sp.kronecker_product(gab, gab) - sp.kronecker_product(gab, gab) * hb).applyfunc(sp.expand)
      if e != 0]
check("the identity column g = αX + βY: g⊗g commutes with XX + XY + YX only on α² − αβ − β² = 0, the golden frame",
      bool(cm) and all(sp.simplify(e / (al ** 2 - al * be - be ** 2)).is_number for e in cm))
G2, G4 = (gmap * gmap).applyfunc(sp.simplify), (gmap ** 4).applyfunc(sp.simplify)
check("the golden map has order four up to scale: its fourth power is a multiple of 1, its square is not",
      G4[0, 0] != 0 and (G4 - G4[0, 0] * sp.eye(4)).applyfunc(sp.simplify).is_zero_matrix
      and not (G2 - G2[0, 0] * sp.eye(4)).applyfunc(sp.simplify).is_zero_matrix)

print("E  the boundary of A (computed)")
lone = ct.Setting(1, [], [((0,), ct.LOWER, 1)])
pair = ct.Setting(2, heisenberg(2), [((l,), ct.LOWER, 1) for l in range(2)])
check("cooling alone, an isolated site: the general solve finds a site-local Pi-reflection",
      cr.machine(lone, 'pi', timeout=None) == 'FOUND')
check("cooling alone, two sites joined by a Heisenberg bond: the general solve finds none",
      cr.machine(pair, 'pi', timeout=None) == 'NONE')

print("F  where fields break M_beta")


def field_triangle(delta, J2, D2):
    """F137's triangle: unit hopping, real on (0, 1) and (1, 2), J + iD on (2, 0); fields (3/10, −1/5, 1/2), ratio 1/4."""
    terms = bond(0, 1, 1, 0, delta) + bond(1, 2, 1, 0, delta) + bond(2, 0, J2, D2, delta)
    terms += [('Z', (l,), h) for l, h in enumerate((R(3, 10), R(-1, 5), R(1, 2)))]
    return ct.Setting(3, terms, [((l,), ct.LOWER, x) for l, x in enumerate(down)]
                      + [((l,), ct.RAISE, x / 4) for l, x in enumerate(down)])


flip_end = ct.lr(sp.diag(1, R(1, 4)) * Z * X, X)         # rho -> diag(1, 1/4)·Z·X·rho·X
flip_mid = ct.lr(sp.diag(1, R(1, 4)) * X, X * Z)         # rho -> diag(1, 1/4)·X·rho·X·Z
flip = [flip_end, flip_mid, flip_end]
check("fields (3/10, −1/5, 1/2), ratio 1/4, Δ = 0, flux ±pi/2: rho -> π·Z₀Z₂·X^N·rho·X^N·Z₁ is a Pi-reflection",
      reflects_pi(field_triangle(0, 0, 1), flip) and reflects_pi(field_triangle(0, 0, -1), flip))
check("controls: at Δ = 7/10 or at flux 0 it is not",
      not reflects_pi(field_triangle(R(7, 10), 0, 1), flip) and not reflects_pi(field_triangle(0, 1, 0), flip))

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILED")
    sys.exit(1)
print("ALL PASS")
