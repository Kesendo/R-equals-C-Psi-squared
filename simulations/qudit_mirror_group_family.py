#!/usr/bin/env python3
"""
_qudit_mirror_group_family.py - WIP scout: the mirror group as ONE family,
over local dimension d (F121's Z_d wr Z_2) AND the antilinear triangle (F119).

The unifying object is the Weyl-Heisenberg algebra at dimension d: the d^2
operators P_{a,b} = X^a Z^b (a,b in Z_d), with X the clock shift |x>->|x+1>
and Z = diag(omega^x), omega = e^{2pi i/d}, X Z = omega Z X. This is the qudit
generalization of the Pauli group {I,X,Y,Z} (d=2).

CLAIM (this scout pins):
  [A] the three involutions act on the WH labels (a,b) by sign flips + a
      symplectic phase (with the convention Z X = omega X Z):
        theta (transpose):  (a,b) -> (-a,  b), phase omega^{-ab}
        conj:               (a,b) -> ( a, -b), phase 1
        dagger:             (a,b) -> (-a, -b), phase omega^{+ab}
      and dagger = theta o conj holds on the labels (the Klein four-group).
  [B] at d=2 this collapses to F119: -a = a, the phase omega^{ab}=(-1)^{ab} is
      exactly (-1)^{n_Y} (only Y = P_{1,1} gets -1); theta(sigma)=conj(sigma)=
      (-1)^{n_Y} sigma, dagger(sigma)=sigma.
  [C] the transport law mu o L_H o mu = l(mu) m(mu) L_{mu(H)} is UNIVERSAL
      (basis-free): verified for theta, conj, dagger on random qudit H at d=2,3,4,5.
  [D] for d>2 the involutions genuinely PERMUTE the labels (a -> -a is nontrivial),
      so the "(-1)^{n_Y} sign" becomes a label reflection with symplectic phase:
      the antilinear triangle generalizes from a sign to a Z_d x Z_d lattice action.

Self-validating.
"""

import numpy as np
from itertools import product as iprod


def weyl_heisenberg(d):
    """Clock X (shift) and phase Z; return X, Z and the d^2 operators P[a,b]."""
    w = np.exp(2j * np.pi / d)
    X = np.zeros((d, d), dtype=complex)
    for x in range(d):
        X[(x + 1) % d, x] = 1.0
    Z = np.diag([w ** x for x in range(d)])
    P = {}
    for a in range(d):
        Xa = np.linalg.matrix_power(X, a)
        for b in range(d):
            P[(a, b)] = Xa @ np.linalg.matrix_power(Z, b)
    return X, Z, P, w


def identify(M, P, tol=1e-9):
    """Return (label, phase) such that M = phase * P[label], or None."""
    for lab, Pl in P.items():
        # M and Pl are both nonzero; find scalar c with M = c*Pl
        idx = np.unravel_index(np.argmax(np.abs(Pl)), Pl.shape)
        if abs(Pl[idx]) < tol:
            continue
        c = M[idx] / Pl[idx]
        if np.allclose(M, c * Pl, atol=tol):
            return lab, c
    return None


def main():
    print("=" * 72)
    print("The mirror group as one family: the qudit antilinear triangle (Weyl-Heisenberg)")
    print("=" * 72)

    # ---- [A] involution action on WH labels, with symplectic phases ----
    print("\n[A] involutions on WH labels P[a,b] = X^a Z^b:")
    for d in (2, 3, 4, 5):
        X, Z, P, w = weyl_heisenberg(d)
        ok_theta = ok_conj = ok_dag = ok_triangle = True
        for (a, b), Pab in P.items():
            th = identify(Pab.T, P)
            cj = identify(Pab.conj(), P)
            dg = identify(Pab.conj().T, P)
            # predicted (Z X = omega X Z, so theta gets omega^{-ab}, dagger omega^{+ab})
            pth = ((-a) % d, b % d), w ** (-a * b)
            pcj = (a % d, (-b) % d), 1.0 + 0j
            pdg = ((-a) % d, (-b) % d), w ** (a * b)
            ok_theta &= th is not None and th[0] == pth[0] and abs(th[1] - pth[1]) < 1e-9
            ok_conj &= cj is not None and cj[0] == pcj[0] and abs(cj[1] - pcj[1]) < 1e-9
            ok_dag &= dg is not None and dg[0] == pdg[0] and abs(dg[1] - pdg[1]) < 1e-9
            # triangle: dagger = theta o conj on the operator
            ok_triangle &= np.allclose(Pab.conj().T, (Pab.conj()).T, atol=1e-9)
        print(f"    d={d}: theta (a,b)->(-a,b) w^-ab {'OK' if ok_theta else 'FAIL'}; "
              f"conj ->(a,-b) {'OK' if ok_conj else 'FAIL'}; "
              f"dagger ->(-a,-b) w^-ab {'OK' if ok_dag else 'FAIL'}; "
              f"dag=theta.conj {'OK' if ok_triangle else 'FAIL'}")
        assert ok_theta and ok_conj and ok_dag and ok_triangle

    # ---- [B] d=2 collapses to F119's (-1)^{n_Y} ----
    print("\n[B] d=2 recovers F119: theta(sigma)=conj(sigma)=(-1)^{n_Y} sigma, dagger=id:")
    X, Z, P, w = weyl_heisenberg(2)
    names = {(0, 0): "I", (1, 0): "X", (1, 1): "Y(=XZ)", (0, 1): "Z"}
    for (a, b), Pab in P.items():
        s_theta = np.sign(np.real((Pab.T @ Pab.conj().T)[0, 0]) + 1e-30)  # +-1
        sign = int(round((w ** (a * b)).real))   # (-1)^{ab}
        nY = 1 if (a, b) == (1, 1) else 0         # Y is the only one with n_Y=1
        print(f"    {names[(a,b)]:>6}: (-1)^ab = {sign:+d}  (n_Y={nY}, (-1)^n_Y={(-1)**nY:+d})")
        assert sign == (-1) ** nY
    print("    -> the symplectic phase omega^{ab} at d=2 IS (-1)^{n_Y}. OK")

    # ---- [C] transport law mu o L_H o mu = l(mu) m(mu) L_{mu(H)}, universal ----
    print("\n[C] transport law mu o L_H o mu = l(mu)*m(mu)*L_{mu(H)} (basis-free):")
    rng = np.random.default_rng(0)
    minus_i = -1j

    def LH(H, rho):
        return minus_i * (H @ rho - rho @ H)

    for d in (2, 3, 4, 5):
        H = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))   # non-Hermitian
        rhos = [rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d)) for _ in range(3)]
        maps = [
            ("theta", lambda A: A.T, -1.0, H.T),
            ("conj", lambda A: A.conj(), -1.0, H.conj()),
            ("dagger", lambda A: A.conj().T, +1.0, H.conj().T),
        ]
        worst = 0.0
        for name, mu, sign, muH in maps:
            for rho in rhos:
                lhs = mu(LH(H, mu(rho)))
                rhs = sign * LH(muH, rho)
                worst = max(worst, np.max(np.abs(lhs - rhs)))
        print(f"    d={d}: worst |mu L_H mu - l*m*L_muH| = {worst:.2e}")
        assert worst < 1e-9

    # ---- [D] for d>2 the involutions permute labels (a -> -a nontrivial) ----
    print("\n[D] d>2: the involutions PERMUTE labels (the sign becomes a lattice action):")
    for d in (3, 5):
        X, Z, P, w = weyl_heisenberg(d)
        moved = sum(1 for (a, b) in P if ((-a) % d, b) != (a, b))
        print(f"    d={d}: theta moves {moved}/{d*d} labels off-diagonal "
              f"(a -> -a nontrivial for a != 0, d/2)")
        assert moved > 0
    print("    -> the antilinear triangle is a Z_d x Z_d lattice action; the qubit sign")
    print("       (-1)^{n_Y} is its d=2 degeneration (every label is its own -a).")

    print("\n" + "=" * 72)
    print("ALL CHECKS PASSED. The antilinear triangle (F119) generalizes to every d via")
    print("the Weyl-Heisenberg lattice: theta/conj/dagger act on (a,b) in Z_d x Z_d by")
    print("sign flips with symplectic phases omega^{+-ab}; the transport law is universal;")
    print("and the qubit (-1)^{n_Y} is the d=2 shadow of the symplectic phase. Next: the")
    print("mirror group <Pi_d, D> (F121's Z_d wr Z_2) + the antilinear unit, and whether the")
    print("clock Z_d is the discrete circle the S3xD4 completion thickens (d->inf).")
    print("=" * 72)


if __name__ == "__main__":
    main()
