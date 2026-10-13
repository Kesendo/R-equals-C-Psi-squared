"""The pairing condition read on the flow: gate for THE_PAIRING_CONDITION.md, subsection
"What U buys on the flow".

With U invertible, [U, H] = 0 and U A_l U^-1 = -A_l for every jump c_l = sqrt(g_l) A_l
(A_l Hermitian, A_l^2 = 1, sigma = sum g_l), the one-sided map Lu(rho) = U rho gives the
relation Lu Lc Lu^-1 = -Lc^dag for Lc = L + sigma. This gate checks:

  A  the relation, exactly (== 0.0), for U = X^N under Z dephasing and for the named row's
     U = SWAP02 Z0 Z1 Z2 at its stated values; controls that break it;
  B  the trace form Tr(B e^{t Lc} A) = Tr(e^{-t Lc}(B U^-1) . U A) against the rounding scale;
     for U = X^N the residual table of Pi L Pi^-1 + L + 2 sigma, L^T - L and
     Pi L^T Pi^-1 + L + 2 sigma (Pi(rho) = rho^T X^N), and the composite's trace form;
  C  the pairing K = Tr(rho_B(t)^dag U rho_A(t)), constant under Lc for two forward runs;
  D  N = 4, block (2, 2): neutrality off the centre line; the departures at small dephasing
     within Krein's bound; the first collision on the line tracked across it (two modes,
     opposite signs); the square-root closure of the lean;
  E  U = X^N: <rho(t), X^N rho(t)> = e^{-2 sigma t} <rho(0), X^N rho(0)> for any state, the cat
     pair's 2 cos(th) sin(th) e^{-2 sigma t}, and the anti-cat's lean -1.

Convention: D(rho) = sum_l g_l (A_l rho A_l - rho), so a Z-dephased coherence of XY weight k
decays at 2 sum over its sites of g_l, and the spectrum is centred at -sigma.

Exactness: part A compares exactly. Lu permutes matrix units up to sign (up to powers of two for
the non-unitary U), so the only place a
float can round is the diagonal rate sum; the rates here are dyadic, or (the named row) a single
jump, where -2g + g and 0 + g are exact. A non-dyadic profile of several rates leaves ~1e-16 on the
diagonal, a rounding of the rate sums; part A prints that reading on one profile and does not gate it.
Scales: parts B, C, E compare the holding rows against a heuristic rounding SCALE
(eps * dim * the norms of the propagators and operands), not a proven bound and not an error
law; a holding row passes within ten times it. For the unitary U the holding rows sit two to four
decades under it and the controls ten to thirteen above; with U = X^N 2^(Nhat^2) (part B2) the scale
carries ||U|| ||U^-1|| and the holding rows sit five to nine decades under it, its controls eight to
eleven above. The forward-only rewrite holds by Ad_U alone and tests nothing one-sided.

Run: python simulations/pairing_condition_flow.py   (writes simulations/results/pairing_condition_flow.txt)
"""
import itertools
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

EPS = np.finfo(float).eps
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)

OUT = []
FAILS = []


def say(line=""):
    print(line)
    OUT.append(line)


def check(name, ok, detail):
    say(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILS.append(name)


def kron(*ms):
    out = np.eye(1, dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


def op(o, l, n):
    return kron(*[o if k == l else I2 for k in range(n)])


def superop(f, d):
    """Matrix of a linear map on d x d matrices, row-major vec (orthonormal matrix units)."""
    m = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d * d):
        e = np.zeros(d * d, dtype=complex)
        e[a] = 1
        m[:, a] = f(e.reshape(d, d)).reshape(-1)
    return m


def lindblad(h, jumps):
    def f(r):
        out = -1j * (h @ r - r @ h)
        for a, g in jumps:
            out += g * (a @ r @ a - r)
        return out
    return f


def bonds_h(n, js, delta=1.0, flux=None):
    """sum_b J_b [2 (s+ s- + s- s+) + delta Z Z] = J (XX + YY + delta ZZ); flux puts a phase on a ring's wrap bond."""
    d = 2 ** n
    h = np.zeros((d, d), dtype=complex)
    bonds = [(l, l + 1) for l in range(n - 1)] + ([(n - 1, 0)] if flux is not None else [])
    for b, (i, j) in enumerate(bonds):
        jb = js[b % len(js)]
        ph = np.exp(1j * flux) if (flux is not None and b == len(bonds) - 1) else 1
        h += jb * 2 * (ph * op(SM.T, i, n) @ op(SM, j, n) + np.conj(ph) * op(SM, i, n) @ op(SM.T, j, n))
        h += jb * delta * op(Z, i, n) @ op(Z, j, n)
    return h


def tr(b, a):
    return np.trace(b @ a)


def apply(s, r):
    d = r.shape[0]
    return (s @ r.reshape(-1)).reshape(d, d)


def rand_state(d, rng):
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    r = a @ a.conj().T
    return r / np.trace(r)


def rand_herm(d, rng):
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    return (a + a.conj().T) / 2


# ---------------------------------------------------------------- the systems
def systems(n):
    """(name, H, jumps, U, expected-to-hold). Dyadic couplings and rates throughout."""
    g = [0.25, 0.5, 0.125, 0.375][:n]
    zj = [(op(Z, l, n), g[l]) for l in range(n)]
    xn = kron(*[X] * n)
    chain = bonds_h(n, [1.0, 0.5, 1.5])
    rows = [
        ("Heisenberg chain, uneven J and g", chain, zj, xn, True),
        ("XY chain + X field", bonds_h(n, [1.0], delta=0.0) + 0.75 * op(X, 1, n), zj, xn, True),
        ("chain + Y0 Z1 (Pi_Z no similarity, L^T != L)", chain + 0.5 * op(Y, 0, n) @ op(Z, 1, n), zj, xn, True),
        ("chain + X0 Y1 Z2 (Pi_Z no similarity, L^T != L)",
         chain + op(X, 0, n) @ op(Y, 1, n) @ op(Z, 2, n), zj, xn, True),
        ("control: chain + DM (XY - YX)",
         chain + sum(op(X, l, n) @ op(Y, l + 1, n) - op(Y, l, n) @ op(X, l + 1, n) for l in range(n - 1)), zj, xn, False),
        ("control: chain + Z field", chain + 0.75 * op(Z, 1, n), zj, xn, False),
        ("control: ring, quarter flux", bonds_h(n, [1.0], flux=np.pi / 2), zj, xn, False),
        ("control: chain + sigma- on site 0", chain, zj + [(None, 0.25)], xn, False),
    ]
    if n == 3:
        # the page's named row as stated there: P3, J(XX+YY+ZZ) at J = 1, X dephasing on site 1 at
        # gamma = 1/20, X field of magnitude 30/100 on sites 0 and 2 with opposite signs (U needs h2 = -h0)
        hn = bonds_h(3, [1.0]) + 0.3 * op(X, 0, 3) - 0.3 * op(X, 2, 3)
        swap02 = np.zeros((8, 8), dtype=complex)
        for i in range(8):
            b = [(i >> (2 - k)) & 1 for k in range(3)]
            b[0], b[2] = b[2], b[0]
            swap02[(b[0] << 2) | (b[1] << 1) | b[2], i] = 1
        un = swap02 @ op(Z, 0, 3) @ op(Z, 1, 3) @ op(Z, 2, 3)
        rows.insert(4, ("named row: X dephasing on site 1, U = SWAP02 Z0 Z1 Z2", hn, [(op(X, 1, 3), 1 / 20)], un, True))
    return rows


def build(h, jumps, n):
    d = 2 ** n
    std = [(a, g) for a, g in jumps if a is not None]
    am = [g for a, g in jumps if a is None]
    sig = sum(g for _, g in std)
    l = superop(lindblad(h, std), d)
    for g in am:  # amplitude damping on site 0, the non-Hermitian-jump control
        s = op(SM, 0, n)
        l = l + g * superop(lambda r: s @ r @ s.conj().T - 0.5 * (s.conj().T @ s @ r + r @ s.conj().T @ s), d)
    return l, sig


# ---------------------------------------------------------------- A: the metric relation, exactly
def part_a():
    say("A  metric relation Lu M Lu^-1 = -M^dag, compared exactly (== 0.0)")
    for n in (3, 4):
        d = 2 ** n
        for name, h, jumps, u, hold in systems(n):
            l, sig = build(h, jumps, n)
            m = l + sig * np.eye(d * d)
            lu = superop(lambda r: u @ r, d)
            lui = superop(lambda r: np.linalg.inv(u) @ r, d)
            res = np.abs(lu @ m @ lui + m.conj().T).max()
            ok = (res == 0.0) if hold else (res > 1e-3)
            check(f"N={n} {name}", ok, f"residual {res:.3g}" + ("" if hold else " (must break)"))
    # an admissible U that is neither unitary nor an involution: X^N 2^Nhat. f(Nhat) commutes with a
    # number-conserving H and with every Z_l, and 2^Nhat has a dyadic inverse, so the exact comparison
    # has a reason (a generic f(Nhat), e.g. 1 + 3 Nhat/8, has inverse entries like 8/11 and its residual
    # is 0.0 or ~1e-15 depending on the evaluation order)
    n, d = 4, 16
    nhat = sum(op((I2 - Z) / 2, l, n) for l in range(n))
    u2 = kron(*[X] * n) @ np.diag(2.0 ** np.diag(nhat).real)
    h = bonds_h(n, [1.0, 0.5, 1.5])
    jumps = [(op(Z, l, n), g) for l, g in enumerate([0.25, 0.5, 0.125, 0.375])]
    l, sig = build(h, jumps, n)
    m = l + sig * np.eye(d * d)
    lu = superop(lambda r: u2 @ r, d)
    lui = superop(lambda r: np.linalg.inv(u2) @ r, d)
    res = max(np.abs((lu @ m) @ lui + m.conj().T).max(), np.abs(lu @ (m @ lui) + m.conj().T).max())
    check("N=4 admissible U = X^N 2^Nhat, not unitary (both association orders)", res == 0.0, f"residual {res:.3g}")
    # Pi(rho) = U rho^T squares to rho -> U rho U^T, which commutes with L; it is Ad_U only for orthogonal U
    sq = superop(lambda r: u2 @ r @ u2.T, d)
    ad = superop(lambda r: u2 @ r @ np.linalg.inv(u2), d)
    c_sq = np.abs(sq @ l - l @ sq).max()
    check("N=4 Pi^2 = (rho -> U rho U^T) commutes with L, and differs from Ad_U",
          c_sq == 0.0 and np.abs(sq - ad).max() > 0.1,
          f"|[Pi^2, L]| = {c_sq:.2g}, |Pi^2 - Ad_U| = {np.abs(sq - ad).max():.3g}")
    # the non-dyadic reading (not gated)
    jumps = [(op(Z, l, n), g) for l, g in enumerate([0.3, 0.7, 0.1, 0.45])]
    l, sig = build(h, jumps, n)
    m = l + sig * np.eye(d * d)
    xn = kron(*[X] * n)
    lu = superop(lambda r: xn @ r, d)
    r_nd = lu @ m @ lu + m.conj().T
    off = np.abs(r_nd - np.diag(np.diag(r_nd))).max()
    say(f"  [read] N=4 non-dyadic rates (0.3, 0.7, 0.1, 0.45), U = X^N: residual {np.abs(r_nd).max():.2g}, "
        f"off the diagonal {off:.2g} (a rounding of the rate sums)")


# ---------------------------------------------------------------- B: the flow identity
def flow_scale(ef, eb, a, b):
    """Heuristic rounding scale of the two sides: eps * dim * ||e^{tM}||_2 ||e^{-tM}||_2 ||A||_F ||B||_F.
    A holding row must sit within ten times it; a breaking row far above it."""
    return EPS * ef.shape[0] * np.linalg.norm(ef, 2) * np.linalg.norm(eb, 2) * np.linalg.norm(a) * np.linalg.norm(b)


def part_b(rng):
    say("")
    say("B  flow identity Tr(B e^{tM} A) = Tr(e^{-tM}(B U^-1) . U A); Pi form for U = X^N")
    pi_map = {0: (1, 1), 1: (0, 1), 2: (3, 1j), 3: (2, 1j)}  # I->X, X->I, Y->iZ, Z->iY
    for n in (3, 4):
        d = 2 ** n
        strings = list(itertools.product(range(4), repeat=n))
        paul = [I2, X, Y, Z]
        mats = {s: kron(*[paul[k] for k in s]) for s in strings}

        def pi_f(r):
            out = np.zeros((d, d), dtype=complex)
            for s in strings:
                c = np.trace(mats[s] @ r) / d
                if c == 0:
                    continue
                ph = 1
                for k in s:
                    ph *= pi_map[k][1]
                out += c * ph * mats[tuple(pi_map[k][0] for k in s)]
            return out

        pi = superop(pi_f, d)
        pinv = np.linalg.inv(pi)
        tsw = superop(lambda r: r.T, d)
        pib = tsw @ pi.T @ tsw
        check(f"N={n} Pi is its own transpose under Tr(B X)", np.abs(pib - pi).max() == 0.0,
              f"max |Pi^b - Pi| = {np.abs(pib - pi).max():.3g}")
        for name, h, jumps, u, hold in systems(n):
            l, sig = build(h, jumps, n)
            m = l + sig * np.eye(d * d)
            ui = np.linalg.inv(u)
            if np.allclose(u, kron(*[X] * n)):
                r_pi0 = np.abs(pi @ l @ pinv + l + 2 * sig * np.eye(d * d)).max()
                r_t0 = np.abs(l.T - l).max()
                r_c0 = np.abs(pi @ l.T @ pinv + l + 2 * sig * np.eye(d * d)).max()
                say(f"  [table] N={n} {name}: |Pi L Pi^-1 + L + 2 sigma| = {r_pi0:.3g}, |L^T - L| = {r_t0:.3g}, "
                    f"|Pi L^T Pi^-1 + L + 2 sigma| = {r_c0:.3g}  (largest entry)")
                if not hold:
                    check(f"N={n} {name}: composite breaks", r_c0 > 0.1, f"{r_c0:.3g}")
                else:
                    check(f"N={n} {name}: composite exact", r_c0 == 0.0, f"{r_c0:.3g}")
                    if "no similarity" not in name:
                        check(f"N={n} {name}: both halves exact", r_pi0 == 0.0 and r_t0 == 0.0,
                              f"{r_pi0:.3g}, {r_t0:.3g}")
                if "Z field" in name:
                    check(f"N={n} {name}: keeps F119's transpose exactly", r_t0 == 0.0, f"|L^T - L| = {r_t0:.3g}")
            if "no similarity" in name:
                # the two factors each fail here, exactly; the composite (and the U form) holds
                r_pi = np.abs(pi @ l @ pinv + l + 2 * sig * np.eye(d * d)).max()
                r_t = np.abs(l.T - l).max()
                r_comp = np.abs(pi @ l.T @ pinv + l + 2 * sig * np.eye(d * d)).max()
                check(f"N={n} {name}: factors fail, composite holds", r_pi > 0.1 and r_t > 0.1 and r_comp == 0.0,
                      f"|Pi L Pi^-1 + L + 2 sigma| = {r_pi:.3g}, |L^T - L| = {r_t:.3g}, |Pi L^T Pi^-1 + L + 2 sigma| = {r_comp:.1e}")
            worst, worst_pi, worst_fwd = 0.0, 0.0, 0.0
            for t in (0.3, 1.1):
                a, b = rand_herm(d, rng), rand_herm(d, rng)
                ef, eb = expm(t * m), expm(-t * m)
                lhs = tr(b, apply(ef, a))
                rhs = tr(apply(eb, b @ ui), u @ a)
                scale = flow_scale(ef, eb, a, b)
                worst = max(worst, abs(lhs - rhs) / scale)
                # applied once more, the relation leaves forward runs only: the two-sided Ad_U
                fwd = tr(u @ b @ ui, apply(ef, u @ a @ ui))
                worst_fwd = max(worst_fwd, abs(lhs - fwd) / scale)
                if np.allclose(u, kron(*[X] * n)):
                    lhs_l = tr(b, apply(expm(t * l), a))
                    rhs_l = np.exp(-2 * sig * t) * tr(apply(pi, a).T, apply(expm(-t * l), apply(pinv, b).T))
                    worst_pi = max(worst_pi, abs(lhs_l - rhs_l) / scale)
            if hold:
                check(f"N={n} {name}", worst <= 10.0, f"deviation / rounding scale = {worst:.2e}")
                check(f"N={n} {name}: forward-only rewrite Tr(UBU^-1 . e^(tM)(UAU^-1))", worst_fwd <= 10.0,
                      f"deviation / rounding scale = {worst_fwd:.2e}")
            else:
                check(f"N={n} {name}", worst > 1e6, f"deviation / rounding scale = {worst:.2e} (must break)")
            if np.allclose(u, kron(*[X] * n)) and hold:
                # the composed Pi form; on the Y0 Z1 and X0 Y1 Z2 rows neither factor holds alone
                # (Pi_Z is no similarity, L^T != L), and the composition still does
                check(f"N={n} {name}, Pi form", worst_pi <= 10.0, f"deviation / rounding scale = {worst_pi:.2e}")


def part_b2(rng):
    say("")
    say("B2 the trace form and its forward-only rewrite with U = X^N 2^(Nhat^2), N = 4: U^2 is not a")
    say("   multiple of 1, so writing U where U^-1 belongs (or the reverse) would show")
    n, d = 4, 16
    nhat = np.diag(sum(op((I2 - Z) / 2, l, n) for l in range(n))).real
    u = kron(*[X] * n) @ np.diag(2.0 ** (nhat ** 2))
    ui = np.diag(2.0 ** (-(nhat ** 2))) @ kron(*[X] * n)
    check("U^2 is not a multiple of the identity", np.ptp(np.abs(np.diag(u @ u))) > 1.0,
          f"diagonal of U^2 spans {np.abs(np.diag(u @ u)).min():.3g} .. {np.abs(np.diag(u @ u)).max():.3g}")
    h = bonds_h(n, [1.0, 0.5, 1.5])
    l, sig = build(h, [(op(Z, k, n), g) for k, g in enumerate([0.25, 0.5, 0.125, 0.375])], n)
    m = l + sig * np.eye(d * d)
    w_right = w_fwd = w_mut1 = w_mut2 = 0.0
    for t in (0.3, 1.1):
        a, b = rand_herm(d, rng), rand_herm(d, rng)
        ef, eb = expm(t * m), expm(-t * m)
        lhs = tr(b, apply(ef, a))
        sc = flow_scale(ef, eb, a, b) * np.linalg.norm(u, 2) * np.linalg.norm(ui, 2)
        w_right = max(w_right, abs(lhs - tr(apply(eb, b @ ui), u @ a)) / sc)
        w_fwd = max(w_fwd, abs(lhs - tr(u @ b @ ui, apply(ef, u @ a @ ui))) / sc)
        w_mut1 = max(w_mut1, abs(lhs - tr(apply(eb, b @ u), u @ a)) / sc)
        w_mut2 = max(w_mut2, abs(lhs - tr(ui @ b @ u, apply(ef, u @ a @ ui))) / sc)
    check("trace form Tr(e^(-tM)(B U^-1) . U A)", w_right <= 10.0, f"deviation / rounding scale = {w_right:.2e}")
    check("forward-only rewrite Tr(U B U^-1 . e^(tM)(U A U^-1))", w_fwd <= 10.0,
          f"deviation / rounding scale = {w_fwd:.2e}")
    check("controls: B U in place of B U^-1, and U^-1 B U in place of U B U^-1, break",
          w_mut1 > 1e6 and w_mut2 > 1e6, f"{w_mut1:.2e}, {w_mut2:.2e}")
    ra, rb = rand_state(d, rng), rand_state(d, rng)
    k0 = tr(rb.conj().T @ u, ra)
    e = expm(2.0 * m)
    dev = abs(tr(apply(e, rb).conj().T @ u, apply(e, ra)) - k0)
    scale = EPS * d * d * np.linalg.norm(e, 2) ** 2 * np.linalg.norm(u, 2) * np.linalg.norm(ra) * np.linalg.norm(rb)
    check("the U-weighted product is constant under M at t = 2 with this U", dev / scale <= 10.0,
          f"drift / rounding scale = {dev / scale:.2e}")


# ---------------------------------------------------------------- C: the conserved pairing
def part_c(rng):
    say("")
    say("C  K(t) = Tr(rho_B(t)^dag U rho_A(t)) under the centred flow, both run FORWARD")
    n, d = 4, 16
    for name, h, jumps, u, hold in systems(n):
        if name.startswith("control: chain + sigma-") or name.startswith("control: ring"):
            continue
        l, sig = build(h, jumps, n)
        m = l + sig * np.eye(d * d)
        ra, rb = rand_state(d, rng), rand_state(d, rng)
        k0 = tr(rb.conj().T @ u, ra)
        devs = []
        for t in (0.5, 1.0, 2.0, 4.0):
            e = expm(t * m)
            a_t, b_t = apply(e, ra), apply(e, rb)
            scale = EPS * d * d * np.linalg.norm(e, 2) ** 2 * np.linalg.norm(ra) * np.linalg.norm(rb)
            devs.append(abs(tr(b_t.conj().T @ u, a_t) - k0) / scale)
        if hold:
            check(f"N=4 {name}", max(devs) <= 10.0,
                  f"K(0) = {k0:.6f}, drift to t = 4 / rounding scale = {max(devs):.2e}")
        else:
            check(f"N=4 {name}", max(devs) > 1e6, f"drift / rounding scale = {max(devs):.2e} (must drift)")


# ---------------------------------------------------------------- D: signatures and the first collision
def block_idx(n, p, q):
    d = 2 ** n
    pc = [bin(i).count("1") for i in range(d)]
    return [k * d + l for k in range(d) for l in range(d) if pc[k] == p and pc[l] == q]


def part_d():
    say("")
    say("D  N = 4, block (2, 2), which Lu maps to itself (Lu sends ket popcount p to N - p);")
    say("   K(r, r) on the modes of the centred block M = L + sigma")
    n, d = 4, 16
    u = kron(*[X] * n)
    idx = block_idx(n, 2, 2)
    lu = superop(lambda r: u @ r, d)[np.ix_(idx, idx)]
    cases = [("XY uniform", bonds_h(n, [1.0], delta=0.0)),
             ("Heisenberg uniform", bonds_h(n, [1.0])),
             ("XY uneven J = (1, 0.6, 1.4)", bonds_h(n, [1.0, 0.6, 1.4], delta=0.0))]

    def block_m(h, g):
        jumps = [(op(Z, l, n), g) for l in range(n)]
        return superop(lindblad(h, jumps), d)[np.ix_(idx, idx)] + n * g * np.eye(len(idx))

    def spectrum(h, g):
        w, v = np.linalg.eig(block_m(h, g))
        lean = np.array([np.vdot(v[:, k], lu @ v[:, k]).real / np.vdot(v[:, k], v[:, k]).real for k in range(len(w))])
        return w, lean, v

    def on_line(x):
        # eigenvalue rounding is ~ eps * dim * |M|; the slowest departure here is ~1e-2 * g
        return abs(x.real) < 1e3 * EPS * len(idx) * max(1.0, abs(x))

    for name, h in cases:
        say(f"  -- {name}")
        # (1) neutrality off the centre line, against the eigenvector conditioning
        w, lean, v = spectrum(h, 0.4)
        off = [abs(lean[k]) for k in range(len(w)) if not on_line(w[k])]
        scale = EPS * len(idx) * np.linalg.cond(v)
        check("off-centre modes neutral at g = 0.4", max(off) <= 10 * scale,
              f"{len(off)} modes, max |K(r,r)|/|r|^2 = {max(off):.1e}, scale eps*dim*cond(V) = {scale:.1e} (gate: 10x)")
        if name == "Heisenberg uniform":
            hz = h + 0.75 * op(Z, 1, n)
            w, lean, v = spectrum(hz, 0.4)
            offz = [abs(lean[k]) for k in range(len(w)) if not on_line(w[k])]
            check("control: + Z field (U no longer commutes with H), off-centre modes lean",
                  max(offz) > 0.1, f"max |K(r,r)|/|r|^2 = {max(offz):.2f}")
        # (2) at g = 0 the block is i times a Hermitian matrix; Krein's bound: of each eigenspace at most
        #     2 min(n+, n-) modes can leave the line at once
        m0 = block_m(h, 0.0)
        hw, hv = np.linalg.eigh(1j * m0)
        groups, cur = [], [0]
        for k in range(1, len(hw)):
            if abs(hw[k] - hw[cur[-1]]) < 1e-9:
                cur.append(k)
            else:
                groups.append(cur)
                cur = [k]
        groups.append(cur)
        bound = 0
        for gidx in groups:
            vs = hv[:, gidx]
            ev = np.linalg.eigvalsh(vs.conj().T @ lu @ vs)
            bound += 2 * min(int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
        left = sum(not on_line(x) for x in spectrum(h, 1e-6)[0])
        check("off the line at g = 1e-6 within Krein's bound from the g = 0 eigenspaces (an upper bound)", left <= bound,
              f"{left} modes off the line, sum of 2 min(n+, n-) = {bound}")
        # (3) the first collision among the modes still on the line: bisection on the count
        n0 = sum(on_line(x) for x in spectrum(h, 1e-6)[0])
        lo, hi = 1e-6, 0.01
        while sum(on_line(x) for x in spectrum(h, hi)[0]) == n0:
            lo, hi = hi, hi * 1.3
        for _ in range(45):
            mid = (lo + hi) / 2
            if sum(on_line(x) for x in spectrum(h, mid)[0]) == n0:
                lo = mid
            else:
                hi = mid
        gstar = lo
        # which centre frequencies are gone just after g*: track them across it
        wb, lb, _ = spectrum(h, gstar * (1 - 1e-7))
        wa, _, _ = spectrum(h, gstar * (1 + 1e-7))
        after = [x.imag for x in wa if on_line(x)]
        gone = sorted((wb[k].imag, lb[k]) for k in range(len(wb))
                      if on_line(wb[k]) and wb[k].imag > 0 and min(abs(wb[k].imag - y) for y in after) > 1e-5)
        ok = len(gone) == 2 and gone[0][1] * gone[1][1] < 0
        check(f"first detected collision (1.3x grid, then bisection), g* = {gstar:.8f}: the modes that leave", ok,
              ", ".join(f"omega {o:.4f} lean {s_:+.4f}" for o, s_ in gone) + "  (two, opposite signs, required)")
        if not ok:
            continue
        mid_om = (gone[0][0] + gone[1][0]) / 2

        # refine g* without the on-line threshold: just below g* the two colliding frequencies
        # approach as sqrt(g* - g), so their squared splitting is linear in g and vanishes at g*
        def split2(g):
            w, _, _ = spectrum(h, g)
            cand = sorted((x.imag for x in w if on_line(x) and abs(x.imag - mid_om) < 0.2 and x.imag > 0),
                          key=lambda y: abs(y - mid_om))[:2]
            return (cand[0] - cand[1]) ** 2

        gref = gstar
        for es in ((1e-3, 3e-4, 1e-4), (1e-4, 3e-5, 1e-5), (1e-5, 3e-6, 1e-6), (1e-6, 3e-7, 1e-7)):
            gs = np.array([gref * (1 - e_) for e_ in es])
            ys = np.array([split2(g) for g in gs])
            cq = np.polyfit(gs - gref, ys, 2)       # quadratic in g through three points
            roots = np.roots(cq)
            roots = roots[np.isreal(roots)].real
            gref = gref + roots[np.argmin(abs(roots))]
        say(f"  [read] g* refined from the squared splitting: {gref:.10f} (threshold bisection {gstar:.10f})")
        gstar = gref
        # (4) the lean closes as a square root: r(eps) = lean^2 / eps = c0 + c1 sqrt(eps) + c2 eps + ...;
        #     fit on the first three points, read the prediction of the fourth
        eps_list = [1e-2, 1e-3, 1e-4, 1e-5]
        rs = []
        for e_ in eps_list:
            w, lean, _ = spectrum(h, gstar * (1 - e_))
            pair = sorted((k for k in range(len(w)) if on_line(w[k]) and w[k].imag > 0),
                          key=lambda k: abs(w[k].imag - mid_om))[:2]
            rs.append(np.mean([lean[k] ** 2 for k in pair]) / e_)  # both members, averaged
        a = np.array([[1, np.sqrt(e_), e_] for e_ in eps_list[:3]])
        c = np.linalg.solve(a, rs[:3])
        pred = c[0] + c[1] * np.sqrt(eps_list[3]) + c[2] * eps_list[3]
        miss = abs(pred - rs[3])
        # the exponent: a lean ~ eps^(1/2) has a finite nonzero limit c0 that already dominates the
        # correction at the largest eps; a lean ~ eps^(3/4) would fit with c0 = 0 and fail here
        check("lean^2 / (1 - g/g*) -> c0, a limit of the series' own size: square-root exponent",
              c[0] > abs(c[1]) * np.sqrt(eps_list[0]) and c[0] > 0.1 * rs[0],
              f"c0 = {c[0]:.4f} against the correction c1 sqrt(1e-2) = {abs(c[1]) * 0.1:.2e} and 0.1 r(1e-2) = {0.1 * rs[0]:.3f}")
        say(f"  [read] lean^2/(1 - g/g*) at 1e-2..1e-5: " + ", ".join(f"{r:.5f}" for r in rs)
            + f"; limit c0 = {c[0]:.5f}; fourth decade predicted to {miss:.1e}")
        # convergence, not just a sign: each decade's step at most a third of the previous one
        # (a series that is already flat to rounding passes as converged)
        dif = [abs(rs[i + 1] - rs[i]) for i in range(3)]
        flat = 1e-10 * c[0]
        conv = all(dif[i + 1] <= dif[i] / 3 or dif[i + 1] <= flat for i in range(2))
        check("the series converges geometrically (each decade's step <= 1/3 of the last)", conv,
              "steps " + ", ".join(f"{x:.1e}" for x in dif))


# ---------------------------------------------------------------- E: the world and its one-sided reading
def part_e(rng):
    say("")
    say("E  U = X^N: the overlap <rho, X^N rho> between the world and its one-sided reading")
    n, d = 4, 16
    xn = kron(*[X] * n)
    rows = systems(n)
    for name, h, jumps, u, hold in [rows[0], rows[1], rows[5]]:
        l, sig = build(h, jumps, n)
        r0 = rand_state(d, rng)
        k0 = tr(r0.conj().T @ xn, r0)
        dev = 0.0
        for t in (0.5, 2.0, 8.0):
            e = expm(t * l)
            r = apply(e, r0)
            scale = EPS * d * d * np.linalg.norm(e, 2) ** 2 * np.linalg.norm(r0) ** 2
            dev = max(dev, abs(tr(r.conj().T @ xn, r) - np.exp(-2 * sig * t) * k0) / scale)
        if hold:
            check(f"N=4 {name}: <rho(t), X^N rho(t)> = e^(-2 sigma t) <rho(0), X^N rho(0)>", dev <= 10.0,
                  f"deviation / rounding scale = {dev:.1e}")
        else:
            check(f"N=4 {name}", dev > 1e6, f"deviation / rounding scale = {dev:.1e} (must break)")
    th = 0.3
    psi = np.zeros(d, dtype=complex)
    psi[0], psi[-1] = np.cos(th), np.sin(th)
    rc = np.outer(psi, psi.conj())
    worst = 0.0
    for js in ([1.0], [0.4, 1.7, 0.9]):
        l, sig = build(bonds_h(n, js), [(op(Z, k, n), [0.25, 0.5, 0.125, 0.375][k]) for k in range(n)], n)
        for t in (0.5, 2.0):
            e = expm(t * l)
            r = apply(e, rc)
            scale = EPS * d * d * np.linalg.norm(e, 2) ** 2
            worst = max(worst, abs(tr(r.conj().T @ xn, r) - 2 * np.cos(th) * np.sin(th) * np.exp(-2 * sig * t)) / scale)
    check("cat pair, theta = 0.3, two couplings: 2 cos sin e^(-2 sigma t)", worst <= 10.0,
          f"deviation / rounding scale = {worst:.1e}")
    ra = np.zeros((d, d), dtype=complex)
    ra[0, 0] = ra[-1, -1] = 0.5
    ra[0, -1] = ra[-1, 0] = -0.5
    la = tr(ra.conj().T @ xn, ra).real / np.linalg.norm(ra) ** 2
    # exact: the anti-cat's four entries are +-1/2, so every product and sum here is dyadic
    check("the anti-cat (|0..0> - |1..1>)/sqrt2 leans -1", la == -1.0, f"lean {float(la)}")
    # the two-copy reading of the product: Tr(rho_B U rho_A) = Tr[(U x 1) SWAP (rho_B x rho_A)]
    sw = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            sw[j * d + i, i * d + j] = 1
    q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    ra_, rb_ = rand_state(d, rng), rand_state(d, rng)
    worst2, miss2 = 0.0, np.inf
    for uu in (xn, q):
        lhs = np.trace(rb_ @ uu @ ra_)
        right = np.trace(np.kron(uu, np.eye(d)) @ sw @ np.kron(rb_, ra_))
        wrong = np.trace(np.kron(uu, np.eye(d)) @ sw @ np.kron(ra_, rb_))
        worst2 = max(worst2, abs(lhs - right))
        miss2 = min(miss2, abs(lhs - wrong))
    check("two copies: Tr(rho_B U rho_A) = Tr[(U x 1) SWAP (rho_B x rho_A)], X^N and a random unitary",
          worst2 < 1e3 * EPS and miss2 > 1e6 * max(worst2, EPS),
          f"max deviation {worst2:.1e}; the other order of the copies misses by at least {miss2:.2g}")


if __name__ == "__main__":
    rng = np.random.default_rng(20261013)
    part_a()
    part_b(rng)
    part_b2(rng)
    part_c(rng)
    part_d()
    part_e(rng)
    say("")
    say(f"{'ALL PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    out = Path(__file__).resolve().parent / "results" / "pairing_condition_flow.txt"
    out.write_text("\n".join(OUT) + "\n", encoding="utf-8", newline="\n")
    sys.exit(1 if FAILS else 0)
