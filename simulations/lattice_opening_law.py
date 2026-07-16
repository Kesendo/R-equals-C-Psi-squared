"""Gate for the lattice opening law (experiments/LATTICE_OPENING_LAW.md).

THE LAW (found playing the bridged lattice, 2026-07-16 evening): seed the cat pair
psi(theta) = cos(theta)|0..0> + sin(theta)|1..1> on the open XY chain with local Z
dephasing (rates gamma_l, any profile, H optional ZZ term), and run the lattice's e
world (normal rule) beside its L world (the one-sided reading X^N rho under the turned
rule). The lattice opening -- the entry-wise distance between the two worlds -- has the
closed form

    max_ij |e(t)[i,j] - L(t)[i,j]| = max(cos^2, sin^2) - cos*sin * exp(-2*Gamma*t),

with Gamma = sum_l gamma_l. In words: OPENING = the heavier sock's weight minus the
LIVING spook. Two separated contributions: the seed's chirality (timeless, lifts the
floor) and the spook (dies under the watching at the full k = N rate and closes exactly
the gap it owns). At theta = 45 deg the floor is 1/2 and the opening is purely the dead
fraction of the spook; at theta = 0 there is no spook and the lattice stands open at 1.

WHY IT IS EXACT AND J-FREE: the cat sector is H-dead. An excitation-conserving hop
needs an excitation next to a hole; |0..0> has no excitation and |1..1> has no hole, so
the XY handshake annihilates both ends (and the ZZ term gives both the SAME diagonal
energy, so the spook collects no phase). The e-world trajectory is pure dephasing:
populations frozen, spook = cos*sin*exp(-2*Gamma*t). The L world is the exact bridge
read-through L[i,j] = e[~i,j] (the site-resolved TURNED rule, gated here too), and the
entry-wise maximum of e - L is max(c^2 - cs*k, cs*k - s^2, ...) = max(c^2,s^2) - cs*k
because c^2 + s^2 = 1 and cs*k <= 1/2.

GATES (exit 0 iff all PASS; N = 2 and 3):
  T0  the site-resolved turned rule: L = X rho / R = rho X under mask
      -2*sum_l gamma_l*(1 - bit_l(i^j)) are exact bridges of the one-sided-watched
      e world (machine zero), and the WRONG (normal) rule breaks at O(1)
  T1  the opening law across theta in {0, 15, 30, 40, 45, 60, 75} deg, one-sided
      (gA = 0) AND two-sided watching, dev at the RK4/expm error floor
  T2  J-independence: two different J values give the SAME opening trajectory
      (the cat sector is H-dead); the ZZ term changes nothing
  T3  the Bell corner: theta = 45 -> all four vertices ONE matrix at t = 0; e = LR
      and L = R exactly for all t (the Bell lattice is TWO worlds, split only by
      the spook-death meter)

Runtime ~10 s. Scout history: a local scratch scout, superseded by (and folded into) this gate.
"""
import sys

import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
s0 = np.eye(2, dtype=complex)

FAILURES = []


def check(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


def run_level(N):
    D = 2 ** N
    print(f"\n===== N = {N} =====")

    def op(single, site):
        ops = [s0] * N
        ops[site] = single
        out = ops[0]
        for o in ops[1:]:
            out = np.kron(out, o)
        return out

    def bit(i, l):
        return (i >> (N - 1 - l)) & 1

    def hamiltonian(J, zz=0.0):
        H = np.zeros((D, D), dtype=complex)
        for b in range(N - 1):
            H += J * (op(sx, b) @ op(sx, b + 1) + op(sy, b) @ op(sy, b + 1))
            if zz:
                Z = np.diag([(1 - 2 * bit(i, b)) * (1 - 2 * bit(i, b + 1)) for i in range(D)])
                H += zz * Z.astype(complex)
        return H

    XN = op(sx, 0)
    for site in range(1, N):
        XN = XN @ op(sx, site)

    def masks(gammas):
        mn = np.zeros((D, D))
        mt = np.zeros((D, D))
        for i in range(D):
            for j in range(D):
                mn[i, j] = -2.0 * sum(g * (bit(i, l) ^ bit(j, l)) for l, g in enumerate(gammas))
                mt[i, j] = -2.0 * sum(g * (1 - (bit(i, l) ^ bit(j, l))) for l, g in enumerate(gammas))
        return mn, mt

    def rk4(rho, H, mask, dt, ticks):
        def rhs(r):
            return -1j * (H @ r - r @ H) + mask * r
        out = [rho.copy()]
        for _ in range(ticks):
            k1 = rhs(rho)
            k2 = rhs(rho + dt / 2 * k1)
            k3 = rhs(rho + dt / 2 * k2)
            k4 = rhs(rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            out.append(rho.copy())
        return out

    def maxdiff(a, b):
        return float(np.max(np.abs(a - b)))

    def cat(theta):
        p = np.zeros(D, dtype=complex)
        p[0], p[D - 1] = np.cos(theta), np.sin(theta)
        return np.outer(p, p.conj())

    dt, ticks = 0.05, 60
    J = 1.0

    print("T0  the site-resolved turned rule (one-sided watching)")
    gammas = [0.0] * N
    gammas[-1] = 0.5                                   # only the last site watched
    mn, mt = masks(gammas)
    H = hamiltonian(J)
    r0 = np.outer(*(lambda p: (p, p.conj()))(np.eye(D, dtype=complex)[1]))   # |0..01><0..01|
    e_tr = rk4(r0, H, mn, dt, ticks)
    l_tr = rk4(XN @ r0, H, mt, dt, ticks)
    r_tr = rk4(r0 @ XN, H, mt, dt, ticks)
    bL = max(maxdiff(l, XN @ e) for l, e in zip(l_tr, e_tr))
    bR = max(maxdiff(r, e @ XN) for r, e in zip(r_tr, e_tr))
    check("one-sided bridges exact", bL < 1e-12 and bR < 1e-12, f"L {bL:.1e}, R {bR:.1e}")
    r_wrong = rk4(r0 @ XN, H, mn, dt, ticks)
    bW = max(maxdiff(r, e @ XN) for r, e in zip(r_wrong, e_tr))
    check("wrong rule breaks O(1)", bW > 1e-2, f"break {bW:.3f}")

    print("T1  the opening law across theta (one- and two-sided)")
    for gam in ([0.0] * (N - 1) + [0.5], [0.3] + [0.0] * (N - 2) + [0.5]):
        G = sum(gam)
        mn, mt = masks(gam)
        worst = 0.0
        for deg in (0, 15, 30, 40, 45, 60, 75):
            th = np.deg2rad(deg)
            c, s = np.cos(th), np.sin(th)
            e_tr = rk4(cat(th), H, mn, dt, ticks)
            l_tr = rk4(XN @ cat(th), H, mt, dt, ticks)
            for k, (e, l) in enumerate(zip(e_tr, l_tr)):
                pred = max(c * c, s * s) - c * s * np.exp(-2 * G * dt * k)
                worst = max(worst, abs(maxdiff(e, l) - pred))
        check(f"law holds, gammas={gam}", worst < 1e-6, f"worst dev {worst:.1e}")

    print("T2  J-independence (the cat sector is H-dead) + ZZ blindness")
    gam = [0.2] * N
    mn, mt = masks(gam)
    th = np.deg2rad(30)
    open_a = [maxdiff(e, l) for e, l in zip(rk4(cat(th), hamiltonian(0.4), mn, dt, ticks),
                                            rk4(XN @ cat(th), hamiltonian(0.4), mt, dt, ticks))]
    open_b = [maxdiff(e, l) for e, l in zip(rk4(cat(th), hamiltonian(1.7), mn, dt, ticks),
                                            rk4(XN @ cat(th), hamiltonian(1.7), mt, dt, ticks))]
    open_c = [maxdiff(e, l) for e, l in zip(rk4(cat(th), hamiltonian(1.0, zz=0.8), mn, dt, ticks),
                                            rk4(XN @ cat(th), hamiltonian(1.0, zz=0.8), mt, dt, ticks))]
    dJ = max(abs(a - b) for a, b in zip(open_a, open_b))
    dZ = max(abs(a - b) for a, b in zip(open_a, open_c))
    check("opening blind to J", dJ < 1e-12, f"J=0.4 vs 1.7: {dJ:.1e}")
    check("opening blind to ZZ", dZ < 1e-12, f"zz=0 vs 0.8: {dZ:.1e}")

    print("T3  the Bell corner: one matrix at t=0, two worlds forever")
    gam = [0.0] * (N - 1) + [0.5]
    mn, mt = masks(gam)
    rb = cat(np.pi / 4)
    e_tr = rk4(rb, H, mn, dt, ticks)
    l_tr = rk4(XN @ rb, H, mt, dt, ticks)
    r_tr = rk4(rb @ XN, H, mt, dt, ticks)
    lr_tr = rk4(XN @ rb @ XN, H, mn, dt, ticks)
    sep0 = max(maxdiff(e_tr[0], w[0]) for w in (l_tr, r_tr, lr_tr))
    e_lr = max(maxdiff(e, lr) for e, lr in zip(e_tr, lr_tr))
    l_r = max(maxdiff(l, r) for l, r in zip(l_tr, r_tr))
    check("t=0: four vertices one matrix", sep0 < 1e-12, f"worst {sep0:.1e}")
    check("e = LR and L = R for all t", e_lr < 1e-12 and l_r < 1e-12,
          f"|e-LR| {e_lr:.1e}, |L-R| {l_r:.1e}")


def main():
    for N in (2, 3):
        run_level(N)
    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: the lattice opening law holds from below at N = 2, 3.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
