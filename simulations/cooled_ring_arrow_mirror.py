"""The cooled XY triangle at flux pi/2: the mirror reverses the arrow, and where it cannot be undone.

Gates for experiments/THE_PALINDROME_WITHOUT_A_CARRIER.md. Under cooling alone (sigma^- on every site, no
heating, no dephasing) the Liouvillian L0 of an XY ring with longitudinal fields is block-triangular in the joint
popcount (p, q) of ket and bra: the feeding sigma^- rho sigma^+ steps (p, q) -> (p - 1, q - 1).

  Part A (exact, integers): a site-local sandwich P carries L0 onto -Lt - c, where Lt is L0 with the feeding
          reversed (sigma^+ rho sigma^-, same anticommutator, same diagonal blocks) and c = sum of the rates.
  Part B (exact, sympy over Q(i)): at a generic rational point no two blocks on one line p - q = const share an
          eigenvalue, so L0 and Lt are both similar to their block diagonal (Sylvester) and L0 ~ -L0 - c.
  Part C (error law, mpmath at 90 digits): at a resonance of the levels (2, 1) and (1, 0) where the uphill
          coupling vanishes and the downhill one does not (solved through the scale-free product c_down*c_up),
          L0 + eps*E splits its eigenvalue at lambda like eps^(1/2) (a Jordan block) and at the mirror point
          mu = -c - lambda like eps (two eigenvectors): the spectrum pairs, the Jordan structure does not, and
          no linear or antilinear map carries L0 onto -L0 - c there. Control: a resonance with both couplings
          nonzero splits like eps^(1/2) on both sides.
  Part D (free fermions, numerical): the levels are sums of single-particle energies of K_p = h(Phi + (p-1)pi)
          + diag(-2 h_l - i gamma_l / 2), and at Phi = pi/2 the two parity sectors have the same single-particle
          spectrum.
  Part E (exact modulo two primes): at two exceptional points, where L has Jordan blocks inside its levels,
          the Jordan structure is mirrored; at flux 0 the same rates give no palindrome.

Run: python simulations/cooled_ring_arrow_mirror.py            (the gates; numpy, scipy, sympy, mpmath)
     python simulations/cooled_ring_arrow_mirror.py --search   (a seeded search: which coupling vanishes)
"""
import itertools
import sys

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import least_squares

N = 3
D = 2 ** N
FAILS = []


def check(ok, label):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        FAILS.append(label)


I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
LO = np.array([[0, 1], [0, 0]], complex)   # sigma^-: |1> -> |0>, |0> the ground state (Z = +1)
RA = LO.T.copy()


def op(sites):
    m = np.array([[1.0]], complex)
    for l in range(N):
        m = np.kron(m, sites.get(l, I2))
    return m


def ad_h(H):          # rho -> -i[H, rho], row-major vec: vec(A rho B) = kron(A, B^T) vec(rho)
    return -1j * (np.kron(H, np.eye(D)) - np.kron(np.eye(D), H.T))


def dissipator_no_feed(J):
    JdJ = J.conj().T @ J
    return -0.5 * np.kron(JdJ, np.eye(D)) - 0.5 * np.kron(np.eye(D), JdJ.T)


def feed(J):
    return np.kron(J, J.conj())


def hamiltonian(fields, J, flux_pi_half=True):
    """XY triangle: J(XX + YY) on (0,1), (1,2); the seam (2,0) carries J(XY - YX), flux pi/2 (or J(XX + YY), flux 0)."""
    H = J * (op({0: X, 1: X}) + op({0: Y, 1: Y}) + op({1: X, 2: X}) + op({1: Y, 2: Y}))
    if flux_pi_half:
        H = H + J * (op({2: X, 0: Y}) - op({2: Y, 0: X}))
    else:
        H = H + J * (op({2: X, 0: X}) + op({2: Y, 0: Y}))
    for l, h in enumerate(fields):
        H = H + h * op({l: Z})
    return H


def generators(fields, rates, J, flux_pi_half=True):
    """(Lambda, Fdown): L0 = Lambda + Fdown, Lt = Lambda + Fdown^T (sigma^+ rho sigma^-; sigma^- is real)."""
    lam = ad_h(hamiltonian(fields, J, flux_pi_half))
    fd = np.zeros((D * D, D * D), complex)
    for l, g in enumerate(rates):
        sm = op({l: LO})
        lam = lam + g * dissipator_no_feed(sm)
        fd = fd + g * feed(sm)
    return lam, fd


def popcount_blocks():
    pc = lambda v: bin(v).count("1")
    idx = {}
    for a in range(D):
        for b in range(D):
            idx.setdefault((pc(a), pc(b)), []).append(a * D + b)
    return idx


IDX = popcount_blocks()


def search(seed=4, starts=120):
    """A seeded search for points where the resonance (2,1) <-> (1,0) meets a vanishing arrow coupling: solve
    gap = 0 and c_down*c_up = 0 (four real equations in seven parameters) from random starts around the page's
    point, then ask which of the two couplings vanished."""
    rng = np.random.default_rng(seed)
    x0 = np.array([-0.2569027905, -1.1561748484, -0.2671923618, 3.6080206858, 0.4111590516, 1.5774718924, 0.2086377941])
    ixA_, ixB_ = IDX[(2, 1)], IDX[(1, 0)]

    def pair(x, target):
        lam, fd = generators(x[:3], x[3:6], x[6])
        ea, ra = np.linalg.eig(lam[np.ix_(ixA_, ixA_)])
        eb, rb = np.linalg.eig(lam[np.ix_(ixB_, ixB_)])
        la, lb = np.linalg.inv(ra), np.linalg.inv(rb)
        i, j = np.argmin(np.abs(ea - target)), np.argmin(np.abs(eb - target))
        cd = lb[j] @ fd[np.ix_(ixB_, ixA_)] @ ra[:, i]
        cu = la[i] @ fd.T[np.ix_(ixA_, ixB_)] @ rb[:, j]
        nd = cd / (np.linalg.norm(lb[j]) * np.linalg.norm(ra[:, i]))
        nu = cu / (np.linalg.norm(la[i]) * np.linalg.norm(rb[:, j]))
        return ea[i], eb[j], cd * cu, nd, nu

    hits = {"c_up": [], "c_down": []}
    for _ in range(starts):
        xs = x0 + rng.normal(0, 0.3, 7) * np.array([1, 1, 1, 1, 1, 1, 0.3])
        xs[3:6] = np.abs(xs[3:6]) + 0.05
        xs[6] = abs(xs[6]) + 0.01
        target = pair(xs, -1.5378 - 0.4535j)[0]

        def res(x):
            a, b, prod, _, _ = pair(x, target)
            return [(a - b).real, (a - b).imag, prod.real, prod.imag]

        try:
            sol = least_squares(res, xs, bounds=([-8] * 3 + [0.01] * 3 + [0.005], [8] * 3 + [20] * 3 + [5]),
                                xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
        if np.max(np.abs(sol.fun)) < 1e-12:
            _, _, _, nd, nu = pair(sol.x, target)
            if abs(nu) < 1e-8 < abs(nd):
                hits["c_up"].append(sol.x[0] - sol.x[2])
            elif abs(nd) < 1e-8 < abs(nu):
                hits["c_down"].append(sol.x[0] - sol.x[2])
    print(f"search (seed {seed}, {starts} starts): converged with c_up = 0 at {len(hits['c_up'])} points, "
          f"with c_down = 0 at {len(hits['c_down'])}")
    for k, v in hits.items():
        print(f"  {k} = 0: h0 - h2 at the hits:", " ".join(f"{d:+.3f}" for d in sorted(v)))


if "--search" in sys.argv:
    search()
    sys.exit(0)

# ---------------------------------------------------------------- Part A
print("Part A: the mirror reverses the arrow (exact, integers)")
fields_q, rates_q, J_q = [1, -2, 3], [2, 5, 7], 1          # a rational point, unequal fields and rates
lam, fd = generators(fields_q, rates_q, J_q)
L0, Lt = lam + fd, lam + fd.T
c = sum(rates_q)


def as_int(M):
    R, Im = np.round(M.real), np.round(M.imag)
    assert np.abs(M - R - 1j * Im).max() == 0.0
    return R.astype(np.int64), Im.astype(np.int64)


def turns(ks):
    return np.diag([np.prod([1j ** (ks[l] * ((a >> (N - 1 - l)) & 1)) for l in range(N)]) for a in range(D)])


Xn, Yn = op({l: X for l in range(N)}), op({l: Y for l in range(N)})
found = None
Lr, Li = as_int(2 * L0)
Tr_, Ti = as_int(-2 * Lt - 2 * c * np.eye(D * D))
for kl in itertools.product(range(4), repeat=N):
    for kr in itertools.product(range(4), repeat=N):
        A, B = turns(kl) @ Xn, Yn @ turns(kr)
        E = np.kron(A, B.T)
        Er, Ei = as_int(E)
        # (Er + i Ei)(Lr + i Li) - (Tr + i Ti)(Er + i Ei) == 0, exactly in integers
        re = Er @ Lr - Ei @ Li - (Tr_ @ Er - Ti @ Ei)
        im = Er @ Li + Ei @ Lr - (Tr_ @ Ei + Ti @ Er)
        if not re.any() and not im.any():
            found = (kl, kr)
            break
    if found:
        break
check(found is not None, f"a sandwich (turns about Z)·X^N·rho·Y^N·(turns) carries L0 onto -Lt - c exactly: turns {found}")
lam0, fd0 = generators(fields_q, rates_q, J_q, flux_pi_half=False)
none_at_flux0 = True
for kl in itertools.product(range(4), repeat=N):
    for kr in itertools.product(range(4), repeat=N):
        E = np.kron(turns(kl) @ Xn, (Yn @ turns(kr)).T)
        if np.abs(E @ (lam0 + fd0) + (lam0 + fd0.T + c * np.eye(D * D)) @ E).max() == 0.0:
            none_at_flux0 = False
check(none_at_flux0, "at flux 0 on the odd triangle no such sandwich exists (F137's e^(2i Phi) = (-1)^N)")
check(not np.array_equal(L0, Lt), "L0 and Lt differ (the arrow is real): the identity is not a Pi-reflection of L0 itself")

# ---------------------------------------------------------------- Part B
print("Part B: the levels the arrow links share no eigenvalue at a generic rational point (exact, sympy over Q(i))")
x = sp.symbols("x")


def block_charpoly(M2, ix):
    sub = M2[np.ix_(ix, ix)]
    rows = [[sp.Integer(int(round(v.real))) + sp.I * sp.Integer(int(round(v.imag))) for v in r] for r in sub]
    return sp.Poly(sp.Matrix(rows).charpoly(x).as_expr(), x, domain="QQ_I")


lam2 = 2 * lam
coprime = True
pairs = 0
for (p, q), ixA in IDX.items():
    for (p2, q2), ixB in IDX.items():
        if p - q == p2 - q2 and p > p2:
            pairs += 1
            g = sp.gcd(block_charpoly(lam2, ixA), block_charpoly(lam2, ixB))
            if g.degree() > 0:
                coprime = False
check(coprime, f"all {pairs} pairs of blocks on a common line p - q have coprime characteristic polynomials")
print("    => L0 and Lt are each similar to their common block diagonal (Sylvester), so L0 ~ Lt ~ -L0 - c")

# ---------------------------------------------------------------- Part C
print("Part C: a resonance where only the downhill arrow chains (mpmath, 90 digits, error law)")
UNIT = {"bonds": ad_h(op({0: X, 1: X}) + op({0: Y, 1: Y}) + op({1: X, 2: X}) + op({1: Y, 2: Y})
                      + op({2: X, 0: Y}) - op({2: Y, 0: X})),
        "fields": [ad_h(op({l: Z})) for l in range(N)],
        "nofeed": [dissipator_no_feed(op({l: LO})) for l in range(N)],
        "feed": [feed(op({l: LO})) for l in range(N)]}


def tomp(M, ix, jx):
    S = M[np.ix_(ix, jx)]
    return mp.matrix([[mp.mpc(float(v.real), float(v.imag)) for v in r] for r in S])   # entries are halves: exact


def assemble(xx, ix, jx, with_feed=False, feed_up=False):
    h, g, J = xx[:3], xx[3:6], xx[6]
    M = tomp(UNIT["bonds"], ix, jx) * J
    for l in range(N):
        M += tomp(UNIT["fields"][l], ix, jx) * h[l] + tomp(UNIT["nofeed"][l], ix, jx) * g[l]
        if with_feed:
            F = UNIT["feed"][l].T if feed_up else UNIT["feed"][l]
            M += tomp(F, ix, jx) * g[l]
    return M


ixA, ixB = IDX[(2, 1)], IDX[(1, 0)]


def resonance(xx, lam_guess):
    """The resonant pair (lambda in A, its partner in B) and the two arrow couplings.

    c_down = <l_B|F_down|r_A> and c_up = <l_A|F_up|r_B> each depend on how the eigenvectors are scaled, the
    product c_down*c_up does not; the solver uses the product, the gates the normalised sizes."""
    LA, LB = assemble(xx, ixA, ixA), assemble(xx, ixB, ixB)
    ea, ra = mp.eig(LA)
    eb, rb = mp.eig(LB)
    i = min(range(len(ea)), key=lambda k: abs(ea[k] - lam_guess))
    j = min(range(len(eb)), key=lambda k: abs(eb[k] - lam_guess))
    la, lb = mp.inverse(ra), mp.inverse(rb)
    up = assemble(xx, ixA, ixB, with_feed=True, feed_up=True) - assemble(xx, ixA, ixB)
    down = assemble(xx, ixB, ixA, with_feed=True) - assemble(xx, ixB, ixA)
    cu = (la[i, :] * up * rb[:, j])[0]
    cd = (lb[j, :] * down * ra[:, i])[0]
    c_up = cu / (mp.norm(la[i, :]) * mp.norm(rb[:, j]))
    c_down = cd / (mp.norm(lb[j, :]) * mp.norm(ra[:, i]))
    return ea[i], eb[j], c_down, c_up, cd * cu, (ea, i), (eb, j)


def refine(start, lam_guess, free, kill_up):
    xx = [mp.mpf(v) for v in start]

    def F(*y):
        z = list(xx)
        for k, v in zip(free, y):
            z[k] = v
        a, b, cd, cu, prod, _, _ = resonance(z, lam_guess)
        out = [mp.re(a - b), mp.im(a - b)]
        return out + ([mp.re(prod), mp.im(prod)] if kill_up else [])

    sol = mp.findroot(F, [xx[k] for k in free], tol=mp.mpf(10) ** (-2 * mp.mp.dps + 10), maxsteps=100)
    sol = sol if isinstance(sol, mp.matrix) else mp.matrix([sol])
    for k, v in zip(free, sol):
        xx[k] = v
    return xx


def nearest(ev, z, k):
    return sorted(abs(e - z) for e in ev)[:k]


def splitting_law(xx, lam_guess, label):
    a, b, cd, cu, _, (ea, i), (eb, j) = resonance(xx, lam_guess)
    Lfull = assemble(xx, list(range(D * D)), list(range(D * D)), with_feed=True)
    ev0 = mp.eig(Lfull, left=False, right=False)
    centre = sum(ev0) / (D * D)
    mu = 2 * centre - a
    rng = np.random.default_rng(99)
    E = mp.matrix([[mp.mpc(rng.normal(), rng.normal()) for _ in range(D * D)] for _ in range(D * D)])
    rows = []
    for k in (20, 28, 36, 44):
        eps = mp.mpf(10) ** (-k)
        ev = mp.eig(Lfull + eps * E, left=False, right=False)
        sl = sorted(ev, key=lambda e: abs(e - a))[:2]
        sm = sorted(ev, key=lambda e: abs(e - mu))[:2]
        rows.append((eps, abs(sl[0] - sl[1]), abs(sm[0] - sm[1])))
    print(f"    {label}: |gap| {mp.nstr(abs(a - b), 3)}, |c_down| {mp.nstr(abs(cd), 4)}, |c_up| {mp.nstr(abs(cu), 4)},"
          f" centre {mp.nstr(centre, 10)}, lambda {mp.nstr(a, 10)}, mu {mp.nstr(mu, 10)}")
    for eps, s1, s2 in rows:
        print(f"      eps 1e-{int(round(-mp.log10(eps)))}: split at lambda {mp.nstr(s1, 4)}, at mu {mp.nstr(s2, 4)}")
    info = dict(a=a, b=b, cd=cd, cu=cu, ea=ea, i=i, eb=eb, j=j, ev0=ev0, centre=centre, mu=mu)
    return rows, info, xx


def constant(vals):
    vals = [float(v) for v in vals]
    return max(vals) / min(vals) - 1.0


mp.mp.dps = 90
TINY = mp.mpf(10) ** -60
# the carrierless point: start from a double-precision solve, keep three parameters, solve four conditions
# (the resonance (2,1) <-> (1,0), two real conditions, and c_down*c_up = 0, two more)
start = [-0.2569027905, -1.1561748484, -0.2671923618, 3.6080206858, 0.4111590516, 1.5774718924, 0.2086377941]
guess = mp.mpc(-1.53782344, -0.45353906)
xx = refine(start, guess, [0, 2, 4, 6], kill_up=True)
rows, info, xx = splitting_law(xx, guess, "c_up forced to 0")
check(abs(info["a"] - info["b"]) < TINY and abs(info["cu"]) < TINY and abs(info["cd"]) > 0.01,
      "Newton converges at 90 digits: resonance gap and c_up below 1e-60, c_down of order 0.1")
def conditions(z, lam_guess):
    a, b, _, _, prod, _, _ = resonance(z, lam_guess)
    return [mp.re(a - b), mp.im(a - b), mp.re(prod), mp.im(prod)]


step = mp.mpf(10) ** -40
jac = mp.matrix(4, 4)
for col, k in enumerate([0, 2, 4, 6]):
    zp, zm = list(xx), list(xx)
    zp[k] += step
    zm[k] -= step
    fp, fm = conditions(zp, guess), conditions(zm, guess)
    for row in range(4):
        jac[row, col] = (fp[row] - fm[row]) / (2 * step)
sv = sorted(mp.svd_r(jac, compute_uv=False))
print(f"    singular values of the 4x4 Jacobian of the conditions: {[mp.nstr(v, 4) for v in sv]}")
check(sv[0] > mp.mpf(10) ** -6 * sv[-1],
      "the Jacobian of the four conditions in the four solved parameters is nonsingular: near this point they "
      "leave a smooth set of dimension three in the seven parameters (implicit function theorem)")

lam = info["a"]
others = [e for blk in [(3, 2)] for e in mp.eig(assemble(xx, IDX[blk], IDX[blk]), left=False, right=False)]
simple = (nearest(info["ea"], lam, 2)[1] > 0.1 and nearest(info["eb"], lam, 2)[1] > 0.1
          and nearest(others, lam, 1)[0] > 0.1)
check(simple, "lambda is simple in each of the two blocks and absent from the third block on its line, (3, 2)")
pal_ok = True
for e in info["ev0"]:
    if nearest(info["ev0"], 2 * info["centre"] - e, 1)[0] > mp.mpf(10) ** -30:
        pal_ok = False
check(pal_ok, "the spectrum pairs about its centre: every mirror image lies within 1e-30 of an eigenvalue")
check(nearest(info["ev0"], info["mu"], 3)[1] < mp.mpf(10) ** -30 < nearest(info["ev0"], info["mu"], 3)[2],
      "mu is a double eigenvalue, like lambda")
check(constant([s1 / mp.sqrt(e) for e, s1, _ in rows]) < 1e-3,
      "split at lambda / eps^(1/2) constant over 24 decades of eps: a Jordan block of size 2")
check(constant([s2 / e for e, _, s2 in rows]) < 1e-3,
      "split at the mirror mu / eps constant: two eigenvectors, so L0 is not similar to -L0 - c at this point")
print("    parameters (fields, rates, J):", [mp.nstr(v, 20) for v in xx])

# control: a resonance with both couplings nonzero
ctrl = [-0.0122262, -0.97007525, -0.01769855, 4.22262647, 0.73046531, 2.24690822, 0.14247116]
xc = refine(ctrl, mp.mpc(-2, 0), [0, 2], kill_up=False)
rows_c, info_c, _ = splitting_law(xc, mp.mpc(-2, 0), "control, both arrows")
check(abs(info_c["a"] - info_c["b"]) < TINY and min(abs(info_c["cd"]), abs(info_c["cu"])) > 0.01,
      "control: resonance met to 1e-60 with both couplings of order 0.01 or more")
check(constant([s1 / mp.sqrt(e) for e, s1, _ in rows_c]) < 1e-3 and constant([s2 / mp.sqrt(e) for e, _, s2 in rows_c]) < 1e-3,
      "control: both sides split like eps^(1/2) (Jordan blocks mirrored), so the law tells the two cases apart")

# ---------------------------------------------------------------- Part D
print("Part D: the free-fermion reading (numerical)")


def single_particle(fields, rates, J, phi):
    M = np.diag([-2 * h - 0.5j * g for h, g in zip(fields, rates)]).astype(complex)
    M[0, 1] = M[1, 0] = M[1, 2] = M[2, 1] = 2 * J
    M[2, 0], M[0, 2] = 2 * J * np.exp(1j * phi), 2 * J * np.exp(-1j * phi)
    return M


fl = [float(v) for v in xx[:3]]
rt = [float(v) for v in xx[3:6]]
Jf = float(xx[6])
a_odd = np.linalg.eigvals(single_particle(fl, rt, Jf, np.pi / 2))
a_even = np.linalg.eigvals(single_particle(fl, rt, Jf, np.pi / 2 + np.pi))
best = min(max(abs(a_odd[i] - a_even[j]) for i, j in zip(range(N), perm)) for perm in itertools.permutations(range(N)))
check(best < 1e-12,
      "at flux pi/2 the odd and the even sector see the same single-particle spectrum (flux pi/2 and -pi/2 are transposes)")
lamf, fdf = generators(fl, rt, Jf)
ok = True
for (p, q), ix in IDX.items():
    ev = np.linalg.eigvals(lamf[np.ix_(ix, ix)])
    kp = a_odd if p % 2 == 1 else a_even
    kq = a_odd if q % 2 == 1 else a_even
    pred = [-1j * (sum(kp[list(A)]) - np.conj(sum(kq[list(B)])))
            for A in itertools.combinations(range(N), p) for B in itertools.combinations(range(N), q)]
    ok &= all(np.min(np.abs(np.array(pred) - e)) < 1e-9 for e in ev)
check(ok, "every block (p, q) has the eigenvalues -i(E_A - conj E_B), sums of single-particle energies of K_p and K_q")

# ---------------------------------------------------------------- Part E
print("Part E: Jordan blocks inside the levels are mirrored (exact modulo two primes, sympy over GF(p))")
from sympy.polys.matrices import DomainMatrix
from sympy import GF, Poly, legendre_symbol, nextprime, sqrt_mod


def modular_similarity(fields, rates_ab, rad, flux_pi_half, p):
    """2L over Z[i, sqrt(rad)] sent to GF(p) (i and sqrt(rad) exist there): is the characteristic polynomial
    palindromic, how many repeated factors are defective, and do the kernel dimensions of f(L)^k match those
    of the mirror factor f(2c - x) for every repeated factor f? A common reduction of an algebraic identity."""
    K = GF(p)
    ri, rr = int(sqrt_mod(p - 1, p)), int(sqrt_mod(rad, p))
    H = hamiltonian(fields, 1, flux_pi_half)
    parts = [(2 * ad_h(H), (1, 0))]
    for l, (a, b) in enumerate(rates_ab):
        sm = op({l: LO})
        parts.append((2 * (dissipator_no_feed(sm) + feed(sm)), (a, b)))
    acc = np.zeros((D * D, D * D), dtype=object)
    for M, (a, b) in parts:
        R, Im = as_int(M)
        acc = (acc + ((R + Im * ri) % p) * ((a + b * rr) % p)) % p
    A = DomainMatrix([[K(int(v)) for v in row] for row in acc], (D * D, D * D), K)
    n = D * D
    c2 = (sum(int(acc[i, i]) for i in range(n)) * pow(n, -1, p)) % p
    cp = Poly([int(v) for v in A.charpoly()], x, modulus=p)
    mirror = lambda f: f.compose(Poly([-1, 2 * c2], x, modulus=p))
    m = mirror(cp)
    palindromic = (m * pow(int(m.LC()), -1, p)) == (cp * pow(int(cp.LC()), -1, p))
    eye = DomainMatrix.eye(n, K)

    def at(f):
        R = DomainMatrix.zeros((n, n), K)
        for c in f.all_coeffs():
            R = R * A + eye * K(int(c))
        return R

    defective, mirrored = 0, True
    for f, mult in cp.factor_list()[1]:
        if mult < 2:
            continue
        F, G = at(f), at(mirror(f))
        Fk, Gk, kf, kg = F, G, [], []
        for k in range(mult):
            kf.append(n - Fk.rank())
            kg.append(n - Gk.rank())
            Fk, Gk = Fk * F, Gk * G
        defective += kf[0] < mult * f.degree()
        mirrored &= kf == kg
    return palindromic, defective, mirrored


def primes_with(rad, count=2, start=10 ** 6):
    out, q = [], start
    while len(out) < count:
        q = nextprime(q)
        if q % 4 == 1 and legendre_symbol(rad, q) == 1:
            out.append(q)
    return out


for name, fields, rates_ab, rad, flux in [
    ("triple exceptional point, uniform field 1, rates (8, 8 + 4 sqrt3, 8 - 4 sqrt3)", [1, 1, 1], [(8, 0), (8, 4), (8, -4)], 3, True),
    ("exceptional point, fields (0, 1/2, 1/2), rates (6, 6 + 2 sqrt7, 6 - 2 sqrt7)", [0, 0.5, 0.5], [(6, 0), (6, 2), (6, -2)], 7, True),
]:
    res = [modular_similarity(fields, rates_ab, rad, flux, q) for q in primes_with(rad)]
    check(all(r[0] and r[1] > 0 and r[2] for r in res),
          f"{name}: palindromic, {res[0][1]} defective repeated factors, every Jordan structure mirrored (two primes)")
res = [modular_similarity([0, 0.5, 0.5], [(6, 0), (6, 2), (6, -2)], 7, False, q) for q in primes_with(7)]
check(all(not r[0] for r in res), "control: the same rates at flux 0 give a characteristic polynomial that is not palindromic")

print()
print("FAILED:" if FAILS else "all gates pass", *FAILS, sep="\n  " if FAILS else "")
sys.exit(1 if FAILS else 0)
