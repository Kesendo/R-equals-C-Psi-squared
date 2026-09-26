"""The N=4 XXZ-Delta response of the F89 octic crossing: the producer's certificate and a recomputation.

simulations/f89_zz_break_gate.py locates the coalescences of the tracked pair by the discriminant
Newton and certifies each one (gap, isolation, alg/geo on the Riesz compression, the zero order of
the pair discriminant, and on the real axis a sign change of that discriminant). These tests read the
committed artifact it writes and hold it to the numbers: one double zero at Delta=0, two simple zeros
per sampled Delta > 0, each a Jordan EP2 at real q on Re lambda = -4. The references are the EP
locations as doubles; this file brackets each of them (the double it holds) to 1e-15 by a sign change of the real
discriminant in 40-digit arithmetic (mpmath), and recomputes one sign change and the crossing's double
zero in double precision, all from its own small builder; the producer is never imported (its result
file is written at the end of a run). No test requires or forbids verdict wording on a page.
"""
import re
from pathlib import Path

import mpmath as mp
import numpy as np


ROOT = Path(__file__).parents[2]
ARTIFACT = ROOT / "simulations" / "results" / "f89_zz_break_gate.txt"
HYPOTHESIS = ROOT / "hypotheses" / "DIABOLIC_BY_INTEGRABILITY.md"

Q_EP = np.sqrt((np.sqrt(13) - 1) / 6)
# Real-axis EP2 locations and their Im lambda; test_the_references_are_bracketed_in_40_digit_arithmetic
# certifies each location to 1e-15 by a sign change of the real discriminant.
REFERENCE = {
    0.02: ((0.656935638035895896, 1.303842206475765), (0.660248940906498777, 1.30790927659822)),
    0.10: ((0.650363027620131020, 1.243522164660881), (0.667060160649329616, 1.265942303756041)),
}
ROW = re.compile(r"^\s+(\d\.\d\d)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d)\s+(\d)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$")


def _rows():
    table = ARTIFACT.read_text(encoding="utf-8").split("COMPLEX-q LOCATOR", 1)[1].split("MECHANISM", 1)[0]
    rows = []
    for line in table.splitlines():
        m = ROW.match(line)
        if m:
            g = m.groups()
            rows.append(dict(delta=float(g[0]), q=float(g[1]), im_q=float(g[2]), re_lam_4=float(g[3]),
                             im_lam=float(g[4]), gap=float(g[5]), alg=int(g[6]), geo=int(g[7]), dep=float(g[8]),
                             order=float(g[9]), f_minus=g[10], f_plus=g[11], character=g[12].strip()))
    return rows


def test_the_crossing_is_a_double_zero_and_every_sampled_delta_splits_it_into_two_ep2s():
    rows = _rows()
    base = [r for r in rows if r["delta"] == 0.0]
    assert len(base) == 1
    assert base[0]["character"] == "DIABOLIC" and base[0]["geo"] == 2 and abs(base[0]["order"] - 2) < 1e-3
    assert abs(base[0]["q"] - Q_EP) < 1e-12
    deltas = sorted({r["delta"] for r in rows if r["delta"] > 0})
    assert deltas == [0.02, 0.05, 0.1, 0.2, 0.5]
    for d in deltas:
        pair = [r for r in rows if r["delta"] == d]
        assert len(pair) == 2
        for r in pair:
            assert r["character"] == "DEFECTIVE EP2", r
            assert (r["alg"], r["geo"]) == (2, 1)
            assert abs(r["order"] - 1) < 1e-3
            # located at the rounding floor: the square-root law keeps the gap near 2*sqrt(dep*u*|M|)
            assert r["gap"] <= 1e-6
            # on the real axis and on the line: rounding level, u*|M| ~ 1e-15
            assert r["im_q"] < 1e-12 and abs(r["re_lam_4"]) < 1e-12
            # the real discriminant changes sign across q* -+ 1e-9: an exact real-q zero lies between
            assert float(r["f_minus"]) * float(r["f_plus"]) < 0
        assert abs(pair[0]["q"] - pair[1]["q"]) > 1e-3


def test_the_located_ep2s_are_the_reference_zeros():
    rows = _rows()
    for d, refs in REFERENCE.items():
        located = sorted((r["q"], r["im_lam"]) for r in rows if r["delta"] == d)
        for (q, im), (q_ref, im_ref) in zip(located, sorted(refs)):
            # model: |dq| ~ (g/s)^2 + u|q| ~ 1e-15 at the reached gap g ~ 1e-8, slope s ~ 0.67 or 1.5
            assert abs(q - q_ref) < 1e-12, (d, q, q_ref)
            assert abs(im - im_ref) < 1e-8, (d, im, im_ref)   # the artifact prints nine decimals


def _sym_block(q, delta):
    """The R-even (SE,DE) block of the N=4 XXZ chain, gamma = 1, q = J (own copy, bit s = site s)."""
    n = 4
    kets = [m for m in range(1 << n) if bin(m).count("1") == 1]
    bras = [m for m in range(1 << n) if bin(m).count("1") == 2]
    basis = [(k, b) for k in kets for b in bras]
    index = {kb: i for i, kb in enumerate(basis)}

    def zz(c):
        return sum(1 if ((c >> b) & 1) == ((c >> (b + 1)) & 1) else -1 for b in range(n - 1))

    def refl(c):
        return sum(1 << (n - 1 - s) for s in range(n) if c & (1 << s))

    L = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, (k, b) in enumerate(basis):
        L[col, col] = -2.0 * bin(k ^ b).count("1") - 1j * q * delta * (zz(k) - zz(b))
        for s in range(n):
            for s2 in (s - 1, s + 1):
                if 0 <= s2 < n:
                    if k & (1 << s) and not k & (1 << s2):
                        L[index[((k & ~(1 << s)) | (1 << s2), b)], col] += -2j * q
                    if b & (1 << s) and not b & (1 << s2):
                        L[index[(k, (b & ~(1 << s)) | (1 << s2))], col] += 2j * q
    cols, done = [], set()
    for col, (k, b) in enumerate(basis):
        if col in done:
            continue
        mirror = index[(refl(k), refl(b))]
        v = np.zeros(len(basis), dtype=complex)
        v[col] = 1.0 if mirror == col else 1 / np.sqrt(2)
        if mirror != col:
            v[mirror] = 1 / np.sqrt(2)
        done.update({col, mirror})
        cols.append(v)
    P = np.column_stack(cols)
    return P.conj().T @ L @ P


def _real_discriminant(q, delta, lam):
    w = np.linalg.eigvals(_sym_block(q, delta))
    a, b = w[np.argsort(np.abs(w - lam))[:2]]
    return (a - b) ** 2


def test_the_real_discriminant_changes_sign_at_the_lower_ep2_and_not_at_the_crossing():
    # At Delta = 0.02 the real-axis spectrum is closed under lambda -> -8 - conj(lambda), so f is real;
    # a pair on the line gives f < 0, a mirror pair f > 0. |f| ~ 0.45 * 1e-9 here against a rounding
    # error of about 1e-15 * gap, so the sign is not a rounding effect.
    q_star, im_star = REFERENCE[0.02][0]
    lam = complex(-4.0, im_star)
    f_minus = _real_discriminant(q_star - 1e-9, 0.02, lam)
    f_plus = _real_discriminant(q_star + 1e-9, 0.02, lam)
    assert abs(f_minus.imag) < 1e-6 * abs(f_minus) and abs(f_plus.imag) < 1e-6 * abs(f_plus)
    assert f_minus.real * f_plus.real < 0
    # at the Delta=0 crossing the discriminant is a double zero: f = s^2 (q - q_EP)^2 >= 0 on both sides,
    # s = 11.5819 the linear gap slope
    lam0 = complex(-4.0, 2 * Q_EP)
    for dq in (-1e-6, 1e-6):
        f = _real_discriminant(Q_EP + dq, 0.0, lam0)
        assert f.real > 0 and abs(f.real / dq ** 2 - 11.5819 ** 2) < 0.1


def test_the_hypothesis_page_quotes_the_located_ep2s():
    text = HYPOTHESIS.read_text(encoding="utf-8")
    for d, refs in REFERENCE.items():
        for q_ref, _ in refs:
            assert f"{q_ref:.10f}" in text, (d, f"{q_ref:.10f}")


def _sym_block_mp(q, delta):
    """The same R-even block in mpmath arithmetic (q, delta as mpf)."""
    n = 4
    kets = [m for m in range(1 << n) if bin(m).count("1") == 1]
    bras = [m for m in range(1 << n) if bin(m).count("1") == 2]
    basis = [(k, b) for k in kets for b in bras]
    index = {kb: i for i, kb in enumerate(basis)}

    def zz(c):
        return sum(1 if ((c >> b) & 1) == ((c >> (b + 1)) & 1) else -1 for b in range(n - 1))

    def refl(c):
        return sum(1 << (n - 1 - s) for s in range(n) if c & (1 << s))

    size = len(basis)
    L = mp.zeros(size, size)
    for col, (k, b) in enumerate(basis):
        L[col, col] = -2 * bin(k ^ b).count("1") - mp.mpc(0, 1) * q * delta * (zz(k) - zz(b))
        for s in range(n):
            for s2 in (s - 1, s + 1):
                if 0 <= s2 < n:
                    if k & (1 << s) and not k & (1 << s2):
                        L[index[((k & ~(1 << s)) | (1 << s2), b)], col] += mp.mpc(0, -2) * q
                    if b & (1 << s) and not b & (1 << s2):
                        L[index[(k, (b & ~(1 << s)) | (1 << s2))], col] += mp.mpc(0, 2) * q
    orbits, done = [], set()
    for col, (k, b) in enumerate(basis):
        if col not in done:
            mirror = index[(refl(k), refl(b))]
            orbits.append((col, mirror))
            done.update({col, mirror})
    P = mp.zeros(size, len(orbits))
    for j, (col, mirror) in enumerate(orbits):
        if col == mirror:
            P[col, j] = 1
        else:
            P[col, j] = P[mirror, j] = 1 / mp.sqrt(2)
    return P.T * L * P


def test_the_references_are_bracketed_in_40_digit_arithmetic():
    # f = (lam_a - lam_b)^2 is real on the real axis (the spectrum is closed under lam -> -8 - conj(lam)).
    # At q_ref -+ 1e-15 it is about s^2 * 1e-15 (s^2 ~ 0.45 or 2.2), while 40-digit arithmetic leaves the
    # pair an error near sqrt(1e-40) = 1e-20 at an EP2: a sign change there puts an exact real-q zero within
    # 1e-15 of each reference.
    mp.mp.dps = 40
    for delta, refs in ((mp.mpf("0.02"), REFERENCE[0.02]), (mp.mpf("0.10"), REFERENCE[0.10])):
        for q_ref, im_ref in refs:
            q = mp.mpf(repr(q_ref))
            lam = mp.mpc(-4, im_ref)
            values = []
            for q_side in (q - mp.mpf("1e-15"), q + mp.mpf("1e-15")):
                w = sorted(mp.eig(_sym_block_mp(q_side, delta), left=False, right=False), key=lambda z: abs(z - lam))
                values.append((w[0] - w[1]) ** 2)
            for f in values:
                assert abs(mp.im(f)) < mp.mpf("1e-20") * abs(f)
            assert mp.re(values[0]) * mp.re(values[1]) < 0, (float(delta), q_ref)


def _block_generators(n):
    """The (SE,DE) block's q-independent parts: dephasing D (diagonal, -2 n_diff), the hopping pattern H
    (L's hopping is 2i*q*H: -1 on a ket hop, +1 on a bra hop) and the ZZ pattern Z (L's ZZ diagonal is
    -i*q*Delta*Z). All three are real integer matrices."""
    kets = [m for m in range(1 << n) if bin(m).count("1") == 1]
    bras = [m for m in range(1 << n) if bin(m).count("1") == 2]
    basis = [(k, b) for k in kets for b in bras]
    index = {kb: i for i, kb in enumerate(basis)}

    def zz(c):
        return sum(1 if ((c >> b) & 1) == ((c >> (b + 1)) & 1) else -1 for b in range(n - 1))

    D = np.diag([-2.0 * bin(k ^ b).count("1") for k, b in basis])
    Z = np.diag([float(zz(k) - zz(b)) for k, b in basis])
    H = np.zeros((len(basis), len(basis)))
    for col, (k, b) in enumerate(basis):
        for s in range(n):
            for s2 in (s - 1, s + 1):
                if 0 <= s2 < n:
                    if k & (1 << s) and not k & (1 << s2):
                        H[index[((k & ~(1 << s)) | (1 << s2), b)], col] -= 1
                    if b & (1 << s) and not b & (1 << s2):
                        H[index[(k, (b & ~(1 << s)) | (1 << s2))], col] += 1
    return D, H, Z


def _commutant_singular_values(D, others):
    """Singular values (ascending) of the linear map C -> ([C, A] for A in others), restricted to the C that
    commute with the diagonal integer D exactly (C[i, j] = 0 unless D_ii = D_jj)."""
    d = np.diag(D)
    size = len(d)
    eye = np.eye(size)
    keep = [i + size * j for j in range(size) for i in range(size) if d[i] == d[j]]
    M = np.vstack([(np.kron(eye, A) - np.kron(A.T, eye))[:, keep] for A in others])
    return np.linalg.svd(M, compute_uv=False)[::-1]


def test_the_commutant_of_the_block_is_10_and_19_at_delta_0_and_2_under_delta():
    # The q-independent commutant of the (SE,DE) block: C with [C, D] = [C, q-coefficient] = 0. At Delta=0 the
    # q-coefficient is 2i*H; at Delta = 1/7 it is i*(2H - Z/7); commuting with H and Z separately is the stronger
    # reading. The matrices are integers, so a null singular value is rounding, about u*|M| ~ 1e-15, while the
    # first nonzero one is 0.0197 or more (the smallest, at N=5 with Delta = 1/7): the rank is read across a gap
    # of more than ten decades.
    for n, dim0 in ((4, 10), (5, 19)):
        D, H, Z = _block_generators(n)
        for others, dim in (([H], dim0), ([2 * H - Z / 7], 2), ([H, Z], 2)):
            s = _commutant_singular_values(D, others)
            assert s[dim - 1] < 1e-12 and s[dim] > 1e-3, (n, dim, s[dim - 1], s[dim])
