"""The three-ladder bridge: Q is the hinge (gate-first).

Claim (read off PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §4): the three integer ladders the project
keeps separate -- the disagreement RUNG k = popcount(i^j), the GIRTH ell, the MOMENT j -- are not
three orthogonal axes. They are the two factors of one F87-hardness coefficient on M = A + gamma*Q:

    P_{m,1} = m * Tr(Q . A^{m-1})

where A = -i[H,.] carries the girth/moment side (A's closed walks on H's graph; t_j = Tr(Z_l H^j)),
and Q = sum_l Z_l (x) Z_l is the rung side (diagonal, Q_x = N - 2k(x)). Q's SPECTRUM is the rung
ladder; Q's ACTION projects A's walks onto the Z-weighted girth moments.

The cell-free m=3 (deg-1) face holds for EVERY Hermitian H (proof line 114):

    GATE A:  3 * Tr(Q . A^2)  ==  6 * 4^N * sum_l c_l^2        (rung-weighted walks == girth moments)
                              ==  3 * sum_x (N - 2k(x)) (A^2)_xx (the rung k literally weights the walks)

where c_l is the single-site-Z Pauli coefficient of H (the ell=1 girth moment t_1^(l) = 2^N c_l).

    GATE B (rung essential):  3 * Tr(A^2)  !=  the girth form   (remove the rung weighting -> the
                                                                  moment projection is gone)

A FIRING gate is the find: it would mean the bridge is mis-stated. Self-validating; run directly.
"""

import numpy as np
from math import comb

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def op(s):
    """Tensor product of a Pauli string, e.g. op('IZI')."""
    m = np.array([[1]], dtype=complex)
    for c in s:
        m = np.kron(m, PAULI[c])
    return m


def site_z(n, l):
    return op('I' * l + 'Z' + 'I' * (n - 1 - l))


def comm_super(h):
    """A = -i[H, .] on coherence space, row-stacking vec: kron(H, I) - kron(I, H^T)."""
    d = h.shape[0]
    idn = np.eye(d, dtype=complex)
    return -1j * (np.kron(h, idn) - np.kron(idn, h.T))


def q_super(n):
    """Q = sum_l Z_l (x) Z_l (diagonal; Q_x = N - 2*popcount(i^j))."""
    d = 2 ** n
    acc = np.zeros((d * d, d * d), dtype=complex)
    for l in range(n):
        zl = site_z(n, l)
        acc += np.kron(zl, zl)
    return acc


def rung_diag(n):
    """k(x) = popcount(i^j) for coherence x = i*d + j."""
    d = 2 ** n
    k = np.zeros(d * d, dtype=int)
    for i in range(d):
        for j in range(d):
            k[i * d + j] = bin(i ^ j).count('1')
    return k


def single_site_z_coeffs(h, n):
    """c_l = coefficient of the single-site string Z_l in H = Tr(H Z_l) / 2^N."""
    return np.array([np.trace(h @ site_z(n, l)).real / 2 ** n for l in range(n)])


def bridge_report(name, h, n):
    """Compute the m=3 deg-1 coefficient three ways and the rung-by-rung decomposition."""
    d = 2 ** n
    A = comm_super(h)
    Q = q_super(n)
    A2 = A @ A
    k = rung_diag(n)
    Qdiag = (n - 2 * k).astype(float)            # = diag(Q), the rung weight on each coherence
    assert np.allclose(np.diag(Q).real, Qdiag), "Q diagonal is not N - 2k"

    # three readings of P_{3,1}
    rung_weighted_walks = 3.0 * np.trace(Q @ A2).real                 # m * Tr(Q . A^2)
    explicit_rung_sum = 3.0 * np.sum(Qdiag * np.diag(A2).real)        # m * sum_x (N-2k) (A^2)_xx
    c = single_site_z_coeffs(h, n)
    girth_moment_form = 6.0 * (4 ** n) * np.sum(c ** 2)               # 6*4^N*sum c_l^2 (ell=1 face)

    # rung-essential control: drop Q (use identity weighting)
    rung_less = 3.0 * np.trace(A2).real                              # m * Tr(A^2), no rung

    # how each rung k contributes to the coefficient (the ladder, made concrete)
    walk_diag = np.diag(A2).real
    per_rung = {kk: 3.0 * np.sum(Qdiag[k == kk] * walk_diag[k == kk]) for kk in sorted(set(k))}

    print(f"\n[{name}]  N={n}")
    print(f"  rung-weighted walks  3*Tr(Q.A^2)         = {rung_weighted_walks:+.6f}")
    print(f"  explicit  3*sum_x (N-2k)(A^2)_xx         = {explicit_rung_sum:+.6f}")
    print(f"  girth-moment  6*4^N*sum c_l^2            = {girth_moment_form:+.6f}")
    print(f"  rung-less  3*Tr(A^2) (no rung weighting) = {rung_less:+.6f}")
    print(f"  per-rung contribution (k -> P)           = "
          + ", ".join(f"{kk}:{v:+.1f}" for kk, v in per_rung.items()))

    # GATE A: the bridge identity (all three readings agree)
    assert abs(rung_weighted_walks - explicit_rung_sum) < 1e-7, "Q-trace != explicit rung sum"
    assert abs(rung_weighted_walks - girth_moment_form) < 1e-7, \
        f"BRIDGE BROKEN: rung-weighted {rung_weighted_walks} != girth-moment {girth_moment_form}"
    assert abs(sum(per_rung.values()) - rung_weighted_walks) < 1e-7, "per-rung sum mismatch"
    return rung_weighted_walks, girth_moment_form, rung_less


def moment_factorization_check(h, n, m):
    """The bridge at a GENERAL rung: P_{m,1} = m*Tr(Q.A^{m-1}) equals the supertrace factorization
    m*(-1)^k * sum_l sum_j (-1)^j C(2k,j) t_j^l t_{2k-j}^l, k=(m-1)/2 (proof line 118). Holds for any
    Hermitian H -- so the rung-weights-the-girth-walks bridge is not an m=3 accident."""
    A = comm_super(h)
    Q = q_super(n)
    k2 = m - 1
    k = k2 // 2
    lhs = m * np.trace(Q @ np.linalg.matrix_power(A, k2)).real
    hpows = [np.linalg.matrix_power(h, j) for j in range(k2 + 1)]
    rhs = 0.0
    for l in range(n):
        zl = site_z(n, l)
        t = [np.trace(zl @ hpows[j]).real for j in range(k2 + 1)]
        rhs += sum((-1) ** j * comb(k2, j) * t[j] * t[k2 - j] for j in range(k2 + 1))
    rhs = m * ((-1) ** k) * rhs
    return lhs, rhs


def main():
    print("=" * 78)
    print("THE THREE-LADDER BRIDGE: Q is the hinge (rung x girth-walks -> girth moments)")
    print("=" * 78)

    # 1. The single-site-Z lift (ell=1, deg=1): the named case. c=(1,0) at N=2 -> P_{3,1}=6*16=96.
    rw, gm, rl = bridge_report("single-site-Z lift  H = Z_0", op('ZI'), 2)
    assert abs(gm - 96.0) < 1e-7, f"expected named P_{{3,1}}=96, got {gm}"
    assert abs(rl - rw) > 1.0, "rung-less should differ (the rung weighting is essential)"

    # 2. Z-component + hopping (ell=1, a real walk graph): c_0=1 -> 6*4^N still.
    bridge_report("Z_0 + 0.7 X_0X_1   H = Z_0 + 0.7 XX", op('ZI') + 0.7 * op('XX'), 2)

    # 3. Soft / bipartite control (pure hopping, no single-site-Z): c=0 -> coefficient 0.
    rw3, gm3, _ = bridge_report("bipartite control  H = X_0X_1 (soft)", op('XX'), 2)
    assert abs(gm3) < 1e-9 and abs(rw3) < 1e-9, "soft control should give 0 (no girth moment)"

    # 4. GENERAL: random Hermitian H (not cherry-picked) -- the m=3 face is cell-free, so the
    #    bridge identity must hold for ANY H. This is the strong gate.
    rng = np.random.default_rng(20260615)
    for n in (2, 3):
        d = 2 ** n
        for trial in range(3):
            re = rng.standard_normal((d, d))
            im = rng.standard_normal((d, d))
            raw = re + 1j * im
            h = raw + raw.conj().T                      # Hermitian
            h = h - np.trace(h) / d * np.eye(d)         # traceless (drop the identity component)
            bridge_report(f"random Hermitian H  trial {trial}", h, n)

    # 5. The bridge at a GENERAL rung (m=3 AND m=5): not an m=3 accident. Random Hermitian H.
    print("\n--- general-rung check: P_{m,1} = m*Tr(Q.A^{m-1}) == supertrace moment factorization ---")
    rng2 = np.random.default_rng(424242)
    for n in (2, 3):
        d = 2 ** n
        for m in (3, 5):
            raw = rng2.standard_normal((d, d)) + 1j * rng2.standard_normal((d, d))
            h = raw + raw.conj().T
            h = h - np.trace(h) / d * np.eye(d)
            lhs, rhs = moment_factorization_check(h, n, m)
            print(f"  N={n} m={m}:  m*Tr(Q.A^{m-1}) = {lhs:+.4f}   moment factorization = {rhs:+.4f}")
            assert abs(lhs - rhs) < 1e-6, f"factorization broke at N={n}, m={m}: {lhs} != {rhs}"

    print("\n" + "=" * 78)
    print("ALL GATES PASS: P_{3,1} = 3*Tr(Q.A^2) = 6*4^N*sum c_l^2 for every H tested.")
    print("Q's spectrum is the rung ladder k; Q weighting A's closed walks IS the F87 hardness,")
    print("and it factorizes into the girth moments t_1 = 2^N c_l. The rung is the hinge.")
    print("Removing Q (rung-less Tr(A^2)) destroys the moment projection -> the rung is essential.")
    print("=" * 78)


if __name__ == "__main__":
    main()
