"""Gate for the h thread (experiments/LATTICE_H_THREAD.md): the lattice meets F131.

A longitudinal field h_l Z_l is X-ODD (X Z X = -Z), so X^N conjugation reflects the
field: X H(h) X = H(-h) exactly. That makes X^N a mirror in the F131 order-sorting
sense (sigma_op = -1 on the h axis, chi = +1 linear), the THIRD sighted instance after
the F71 site reversal R (Theorem A) and the antiunitary Floquet Theta (Theorem B).
Four findings, each gated below:

  T0  the mirror crosses the field: X H(h) X = H(-h) exactly; the XY handshake is
      X-even; the DOUBLE turn crosses worlds, X e_{+h}(t) X = the normal-rule world
      run with field -h from the conjugated seed, exactly.
  T1  the one-sided reading is a MIXED TWO-FIELD PENCIL: with a field, L = X^N rho(t)
      satisfies dL/dt = -i*(H(-h) L - L H(+h)) + turned mask -- ket leg in the mirror
      field, bra leg in the original -- and is NO single-field world at all (both
      single-field attempts break at O(1)). The naive "crossed bridge L = X e_{-h}"
      guess was falsified by this gate's scout run (0.819).
  T2  the opening law is h-ROBUST with the UNCHANGED closed form, for any h profile:
      both worlds ride the same cell detunes E_i - E_j, so every candidate magnitude
      of the opening is h-free (experiments/LATTICE_OPENING_LAW.md). NAMING CARE: the
      world measured here is the TURNED-WATCHING world run single-field (the opening
      law's L vertex), NOT T1's mixed-pencil reading X^N rho(t), whose distance to e
      IS phase-dependent under a field.
  T3  the four-cell trajectory face for M = X^N, swept over ALL 63 Pauli readouts on
      an operator-X-even population prep: every ALIVE readout sorts EXACTLY by its
      X-parity q (q = +1 reads even in h, q = -1 reads odd), the odd cell is
      non-vacuously populated, and the dead strings include every population readout:
      a SECOND mirror (entrywise conjugation composed with the chiral sublattice
      gauge, the same ingredients as Theorem B's Theta = T*K -- a rhyme, not an
      identification) forces all populations EVEN in h (gated exactly), so X-odd
      diagonal readouts like Z_0 are DOUBLY-MIRRORED ZEROS: identically zero at every
      h because the two mirrors disagree there. Which non-diagonal strings die is
      OBSERVED here (48/63 at N = 3), not derived; the derivation is open.

Runtime ~30 s. Scout history: a local scratch scout (falsified the crossed-bridge
guess, found the enumeration), folded into this gate.
"""
import sys
from itertools import product

import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
s0 = np.eye(2, dtype=complex)
PAULI = {"I": s0, "X": sx, "Y": sy, "Z": sz}

FAILURES = []


def check(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


N = 3
D = 2 ** N
J = 1.0


def op(single, site):
    ops = [s0] * N
    ops[site] = single
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def bit(i, l):
    return (i >> (N - 1 - l)) & 1


def hamiltonian(hs):
    H = np.zeros((D, D), dtype=complex)
    for b in range(N - 1):
        H += J * (op(sx, b) @ op(sx, b + 1) + op(sy, b) @ op(sy, b + 1))
    for l in range(N):
        H += hs[l] * np.diag([1 - 2 * bit(i, l) for i in range(D)]).astype(complex)
    return H


XN = op(sx, 0)
for site in range(1, N):
    XN = XN @ op(sx, site)


def masks(gammas):
    mn = np.zeros((D, D))
    mt = np.zeros((D, D))
    for i in range(D):
        for j in range(D):
            mn[i, j] = -2 * sum(g * (bit(i, l) ^ bit(j, l)) for l, g in enumerate(gammas))
            mt[i, j] = -2 * sum(g * (1 - (bit(i, l) ^ bit(j, l))) for l, g in enumerate(gammas))
    return mn, mt


def rk4(rho, Hket, mask, dt, ticks, Hbra=None):
    Hb = Hket if Hbra is None else Hbra

    def rhs(r):
        return -1j * (Hket @ r - r @ Hb) + mask * r

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


def main():
    gam = [0.0, 0.0, 0.5]                       # one-sided watching throughout
    mn, mt = masks(gam)
    dt, ticks = 0.05, 40
    hs = [0.4, 0.4, 0.4]
    hneg = [-v for v in hs]

    print("T0  the mirror crosses the field")
    check("X H(h) X = H(-h)", maxdiff(XN @ hamiltonian(hs) @ XN, hamiltonian(hneg)) < 1e-12,
          f"residual {maxdiff(XN @ hamiltonian(hs) @ XN, hamiltonian(hneg)):.1e}")
    check("the XY handshake is X-even",
          maxdiff(XN @ hamiltonian([0.0] * N) @ XN, hamiltonian([0.0] * N)) < 1e-12, "residual 0")
    psi = np.zeros(D, dtype=complex)
    psi[1] = 1.0
    r0 = np.outer(psi, psi.conj())
    e_plus = rk4(r0, hamiltonian(hs), mn, dt, ticks)
    lr_minus = rk4(XN @ r0 @ XN, hamiltonian(hneg), mn, dt, ticks)
    crossedLR = max(maxdiff(lr, XN @ e @ XN) for lr, e in zip(lr_minus, e_plus))
    check("the double turn crosses: world(-h) = X e_{+h} X", crossedLR < 1e-12, f"worst {crossedLR:.1e}")

    print("T1  the one-sided reading is a mixed two-field pencil")
    l_mixed = rk4(XN @ r0, hamiltonian(hneg), mt, dt, ticks, Hbra=hamiltonian(hs))
    mixed = max(maxdiff(l, XN @ e) for l, e in zip(l_mixed, e_plus))
    check("ket -h / bra +h pencil exact", mixed < 1e-12, f"worst {mixed:.1e}")
    l_same = rk4(XN @ r0, hamiltonian(hs), mt, dt, ticks)
    same = max(maxdiff(l, XN @ e) for l, e in zip(l_same, e_plus))
    e_minus = rk4(r0, hamiltonian(hneg), mn, dt, ticks)
    l_flip = rk4(XN @ r0, hamiltonian(hneg), mt, dt, ticks)
    flip = max(maxdiff(l, XN @ e) for l, e in zip(l_flip, e_minus))
    check("both single-field attempts break O(1)", same > 1e-2 and flip > 1e-2,
          f"same-field {same:.3f}, crossed-single {flip:.3f}")

    print("T2  the opening law is h-robust (unchanged closed form, any profile)")

    def cat(theta):
        p = np.zeros(D, dtype=complex)
        p[0], p[D - 1] = np.cos(theta), np.sin(theta)
        return np.outer(p, p.conj())

    G = sum(gam)
    worst = 0.0
    for hprof in ([0.0] * N, [0.4] * N, [0.7, -0.2, 0.3]):
        H = hamiltonian(hprof)
        for deg in (0, 30, 45, 75):
            th = np.deg2rad(deg)
            c, s = np.cos(th), np.sin(th)
            e_tr = rk4(cat(th), H, mn, dt, ticks)
            l_tr = rk4(XN @ cat(th), H, mt, dt, ticks)
            for k, (e, l) in enumerate(zip(e_tr, l_tr)):
                pred = max(c * c, s * s) - c * s * np.exp(-2 * G * dt * k)
                worst = max(worst, abs(maxdiff(e, l) - pred))
    check("opening = max(c^2,s^2) - cs*e^(-2Gt) under h", worst < 1e-6, f"worst dev {worst:.1e}")

    print("T3  the four cells for M = X^N, all 63 readouts")
    s_ = 1
    rb = np.zeros((D, D), dtype=complex)
    rb[s_, s_] = 0.5
    rb[D - 1 - s_, D - 1 - s_] = 0.5                    # operator-X-even population prep
    check("prep operator-X-even", maxdiff(XN @ rb @ XN, rb) < 1e-14, "X rho X = rho")

    t_scan, hdir = 0.4, [1.0, -0.6, 0.3]

    def traj(hprof, O):
        return [float(np.real(np.trace(r @ O)))
                for r in rk4(rb, hamiltonian(hprof), mn, dt, ticks)]

    p_plus = [np.real(np.diag(r)).copy() for r in rk4(rb, hamiltonian([t_scan * v for v in hdir]), mn, dt, ticks)]
    p_minus = [np.real(np.diag(r)).copy() for r in rk4(rb, hamiltonian([-t_scan * v for v in hdir]), mn, dt, ticks)]
    pop_even = max(float(np.max(np.abs(a - b))) for a, b in zip(p_plus, p_minus))
    check("the second mirror: ALL populations even in h", pop_even < 1e-12, f"worst dev {pop_even:.1e}")

    alive, sorted_ok, n_odd, n_dead = [], True, 0, 0
    for letters in product("IXYZ", repeat=N):
        name = "".join(letters)
        if name == "I" * N:
            continue
        O = op(PAULI[letters[0]], 0)
        for site in range(1, N):
            O = O @ op(PAULI[letters[site]], site)
        q = (-1) ** sum(1 for ch in name if ch in "YZ")
        fp = traj([t_scan * v for v in hdir], O)
        fm = traj([-t_scan * v for v in hdir], O)
        mag = max(abs(a) for a in fp)
        if mag < 1e-12:
            n_dead += 1
            continue
        odd = max(abs(a + b) for a, b in zip(fp, fm))
        even = max(abs(a - b) for a, b in zip(fp, fm))
        ok = (odd < 1e-10) if q < 0 else (even < 1e-10)
        sorted_ok = sorted_ok and ok
        alive.append((name, q, mag))
        if q < 0:
            n_odd += 1
    check("every alive readout sorts by its X-parity", sorted_ok,
          f"{len(alive)}/63 alive ({n_dead} dead; the diagonal deaths are the proven "
          "doubly-mirrored zeros), all sorted")
    top_odd = max((m for _, q, m in alive if q < 0), default=0.0)
    check("the odd cell is non-vacuous", n_odd > 0 and top_odd > 0.1,
          f"{n_odd} X-odd strings alive, top |<O>| = {top_odd:.3f}")
    zdead = all(m < 1e-12 for m in
                [max(abs(a) for a in traj([t_scan * v for v in hdir],
                                          np.diag([1 - 2 * bit(i, l) for i in range(D)]).astype(complex)))
                 for l in range(N)])
    check("every single-site Z is a doubly-mirrored zero", zdead, "Z_l identically 0 at every h")

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: X^N is a third F131 mirror; the double turn crosses fields; the one-sided "
          "reading is a mixed two-field pencil; the opening law is h-robust; the alive readouts "
          "sort by X-parity; the DIAGONAL dead ones are proven doubly-mirrored zeros (the "
          "off-diagonal deaths are observed, derivation open).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
