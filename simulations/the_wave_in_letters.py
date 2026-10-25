"""The wave's look at itself, read in letters.

The letter-cube route behind the section "How far the look at itself reaches" of
reflections/ON_THE_WAVE_RECOGNIZING_THAT_IT_IS_NOT_ALONE.md (docs/THE_ONE_SQUARE.md §7 to §9).
Nothing here propagates: no e^{tL}, no eigensolver. A state is read as the weights it puts on
words (Pauli strings), each word at a height k_Z, and every sentence of the section is a
statement about heights, words and the two moves of the mirror, R (ρ ↦ ρ·X^N, the height flip
k ↦ N − k) and D (the transpose, the sign (−1)^{n_Y}, no word moves). The matrix script
simulations/the_wave_between.py is the second route of the handshake.

1. The letters: coordinates, heights, the two edges, R's flip of every word's height with a unit
   phase; D's sign and the rate 2γ·k_Z of every word checked against the word's own matrix
   (its transpose; Σ_l Z_l σ Z_l = (N − 2k_Z)·σ), N = 1..4.
2. One bit, watched: W_0 = (1 + z²)/2, W_1 = (x² + y²)/2 exactly (sympy); equal weights with the
   ball give z = 0 and x² + y² = 1, the equator.
3. Several bits, watched: the height profiles of |+⟩^N (binomial), the cat and the pair (half on
   each edge), a random product of equator bits (symmetric), the W state and a tilted product
   (not; the W state at N = 3, 4, since at N = 2 it is the pair |01⟩ + |10⟩ and pairs); W_k = W_{N−k}
   is the twin condition under watching alone.
4. The second move is invisible: a term with an even Y count joins two words of opposite Y
   parity; the step back is exactly minus the step there for EVERY Hermitian term (so
   for an even-Y term (−1)^{n_Y(σ)+n_Y(τ)} = −1 on every moved pair and D·L_H·D = −L_H = L_Hᵀ); a
   term with an odd Y count joins words of the same Y parity, so D does not reverse its steps
   (checked with commutator_coefficient, exact, N = 1..3, for XX, YY, ZZ, X, Z and the Y field). The letters
   say where the rule breaks; that the conjugate pair then parts is the matrix route's (and at
   N = 1 it does not part: the field steps X ↔ Z and leaves Y standing, cube_step).
5. The conjugate states: a random X^N eigenstate (either sign) is returned by R up to the sign
   (R·a = ± a, so |a_{Rσ}| = |a_σ| word for word), while D moves it (D·a ≠ ± a), so its mirror
   image Π·a = R·D·a is its conjugate D·a up to sign and it is not mirror-fixed; the cat, the pair
   and |+⟩^N (and at N = 2 the singlet and the triplet) are returned by both moves,
   mirror-fixed; the completely mixed state carries only I^N and nothing on X^N, where R would
   need it. A Y^N eigenstate: R·a = phase · (Z^N-shift of a) word for word, and the mirror image
   R(D·a) = phase · (Z^N-shift of D·a), the same heights as a.
6. The turn about Z keeps every height: the height profile of a product of equator bits at any
   angle is the binomial one; XX + YY and ZZ are carried to themselves by the turn, XX alone is
   not (symbolic, one bond).
7. A field along the watching: on one bit it steps X ↔ Y inside height 1 and leaves I and Z
   standing (cube_step); every normalized state has weight 1/2^N on I^N (so the heights cannot
   pair when height 1 holds only turning words); on two bits the field's moves on the words of
   |01⟩ ± |10⟩ cancel (commutator coefficients summed), and so do its moves on the bond's own
   words XX + YY and ZZ, so the cancellation holds along the run; on Bell+ they do not.
8. A field across the watching, along X: on one bit it steps Y ↔ Z between the heights and
   leaves I and X standing; R pairs I with X and Y with Z, standing with standing and turning
   with turning. Which weights must balance (x² = 1, y² = z²) is the matrix route's (section 7
   of the_wave_between.py); the letters give the pairs.
"""
from __future__ import annotations

import sys
from itertools import product
from math import comb
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding="utf-8")
import framework as fw  # noqa: E402
from framework.pauli import _pauli_label  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_all(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


def labels(N):
    return [''.join(t) for t in product('IXYZ', repeat=N)]


LETTER = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def mat(label):
    """The word's own matrix, for the checks of section 1 only."""
    return kron_all([LETTER[c] for c in label])


def coeffs(rho, N):
    return fw.pauli_basis_vector(rho, N)


def words(c, N, tol=1e-12):
    """The words a coefficient vector carries, label -> coefficient."""
    return {_pauli_label(k, N): c[k] for k in range(4 ** N) if abs(c[k]) > tol}


def state(psi):
    psi = np.asarray(psi, dtype=complex)
    return np.outer(psi, psi.conj()) / np.vdot(psi, psi).real


def letters():
    print("1. The letters: every word sits at (k_Z, k_X, k_Y); its height k_Z is the number of X or Y letters;")
    print("   R flips every height k -> N - k with a unit phase, D moves no word and signs the Y letters")
    for L, k in (("I", (0, 0, 0)), ("X", (1, 0, 1)), ("Y", (1, 1, 0)), ("Z", (0, 1, 1))):
        assert fw.cube_coords(L) == k
    for N in (1, 2, 3, 4):
        quiet = [s for s in labels(N) if fw.cube_coords(s)[0] == 0]
        far = [s for s in labels(N) if fw.cube_coords(s)[0] == N]
        assert len(quiet) == 2 ** N and len(far) == 2 ** N, "the two edges hold 2^N words each"
        assert all(set(s) <= set('IZ') for s in quiet) and all(set(s) <= set('XY') for s in far)
        for s in labels(N):
            t, ph = fw.cube_right_shift(s, 'X')
            assert fw.cube_coords(t)[0] == N - fw.cube_coords(s)[0] and abs(ph) == 1, "R: height k -> N - k"
            M = mat(s)
            assert np.array_equal(M.T, fw.cube_transpose_sign(s) * M), "D: the word's own transpose is its Y sign"
            conj_sum = sum(mat(s[:l] + 'Z' + s[l + 1:]) @ M @ mat(s[:l] + 'Z' + s[l + 1:]) for l in range(N))
            assert np.array_equal(conj_sum, (N - 2 * fw.cube_coords(s)[0]) * M), "the watching: sum_l Z_l sigma Z_l = (N - 2 k_Z) sigma, so the word fades at 2 gamma k_Z and moves nowhere"
            assert fw.cube_rate(s, 0.3) == 2 * 0.3 * fw.cube_coords(s)[0]
        assert {fw.cube_right_shift(s, 'X')[0] for s in quiet} == set(far), "R carries the quiet edge onto the far edge"
        print(f"   N = {N}: {4 ** N} words, quiet edge {len(quiet)}, far edge {len(far)}, R quiet -> far exactly, D signs {sum(1 for s in labels(N) if s.count('Y') % 2)} words")


def one_bit_watched():
    print("2. One bit, watched: W_0 = (1 + z^2)/2 (the word I alone is one half), W_1 = (x^2 + y^2)/2;")
    print("   equal with the ball only for z = 0, x^2 + y^2 = 1: the equator")
    x, y, z = sp.symbols('x y z', real=True)
    coeff = {'I': sp.Rational(1, 2), 'X': x / 2, 'Y': y / 2, 'Z': z / 2}       # rho = (I + xX + yY + zZ)/2
    W = {0: 0, 1: 0}
    for L, c in coeff.items():
        W[fw.cube_coords(L)[0]] += 2 * c ** 2                                  # 2^N |c|^2, N = 1
    assert sp.simplify(W[0] - (1 + z ** 2) / 2) == 0 and sp.simplify(W[1] - (x ** 2 + y ** 2) / 2) == 0
    assert sp.simplify(W[0] - W[1] - (1 + z ** 2 - x ** 2 - y ** 2) / 2) == 0
    # equal weights: x^2 + y^2 = 1 + z^2. With r^2 = x^2 + y^2 + z^2 the far edge's shortfall is
    # W_1 - W_0 = (r^2 - 1 - 2 z^2)/2 <= -z^2 <= 0 inside the ball, zero only for r = 1 and z = 0
    r2 = x ** 2 + y ** 2 + z ** 2
    shortfall = sp.simplify(W[1] - W[0])
    assert sp.simplify(shortfall - ((r2 - 1) / 2 - z ** 2)) == 0, "the far edge's shortfall is (r^2 - 1)/2 - z^2: <= 0 inside the ball, 0 only at r = 1, z = 0"
    on_surface = sp.expand(shortfall.subs(x ** 2, 1 - y ** 2 - z ** 2))                 # r = 1
    assert sp.simplify(on_surface + z ** 2) == 0 and sp.solve(on_surface, z) == [0], "on the surface the shortfall is -z^2: zero only at z = 0"
    print(f"   W_0 = {sp.simplify(W[0])}, W_1 = {sp.simplify(W[1])}; W_0 - W_1 = {sp.factor(W[0] - W[1])}")


def several_bits_watched():
    print("3. Several bits, watched: the twin condition is W_k = W_{N-k}; equator products are binomial, the pair and the")
    print("   cat sit half on each edge, the W state and a tilted product do not pair")
    rng = np.random.default_rng(3)
    for N in (2, 3, 4):
        d = 2 ** N
        plus = kron_all([np.array([[1], [1]])] * N).ravel()
        ghz = np.zeros(d); ghz[0] = ghz[-1] = 1
        w = np.zeros(d)
        for k in range(N):
            w[1 << k] = 1
        random_eq = kron_all([np.array([[1], [np.exp(1j * th)]]) for th in rng.uniform(0, 2 * np.pi, N)]).ravel()
        tilted = kron_all([np.array([[1], [1]])] * (N - 1) + [np.array([[np.cos(0.4)], [np.sin(0.4)]])]).ravel()
        cases = [("|+>^N", plus, True), ("cat", ghz, True), ("random equator product", random_eq, True), ("tilted product", tilted, False)]
        if N > 2:
            cases.append(("W", w, False))
        print(f"   N = {N}:")
        for name, psi, twin in cases:
            W = fw.cube_height_weights(coeffs(state(psi), N), N)
            sym = np.abs(W - W[::-1]).max()
            assert (sym < 1e-12) == twin, f"{name}: |W - reversed| = {sym:.1e}"
            if name == "|+>^N":
                assert np.abs(W - np.array([comb(N, k) for k in range(N + 1)]) / 2 ** N).max() < 1e-12, "binomial"
            if name == "cat":
                assert abs(W[0] - 0.5) < 1e-12 and abs(W[N] - 0.5) < 1e-12 and np.abs(W[1:N]).max() < 1e-12, "half on each edge"
            print(f"      {name:24s} W = {np.round(W, 4)}  |W - reversed| {sym:.1e}")


def second_move_invisible():
    print("4. The second move is invisible: a term of even Y count joins words of opposite Y parity and the step back is")
    print("   minus the step there, so D turns every step round; an odd-Y term joins words of equal parity and D does not")
    for N in (1, 2, 3):
        def placed(template):
            return ['I' * i + template + 'I' * (N - len(template) - i) for i in range(N - len(template) + 1)]
        terms = {"X field": placed("X"), "Z field": placed("Z"), "Y field": placed("Y")}
        if N > 1:
            terms.update({"XX": placed("XX"), "YY": placed("YY"), "ZZ": placed("ZZ")})
        for name, ts in terms.items():
            y_count = ts[0].count('Y')
            moved = 0
            for t in ts:
                for s in labels(N):
                    r = fw.commutator_coefficient(t, s)
                    if r is None:
                        continue
                    tau, c = r
                    back = fw.commutator_coefficient(t, tau)
                    assert back is not None and back[0] == s, "a term that moves sigma to tau moves tau back to sigma"
                    assert back[1] == -c, "the step back is minus the step there (Hermitian words: the phase is +-i)"
                    parity_flip = fw.cube_transpose_sign(s) * fw.cube_transpose_sign(tau) == -1
                    assert parity_flip == (y_count % 2 == 0), f"{name}: an even-Y term joins words of opposite Y parity, an odd-Y term words of equal parity"
                    moved += 1
            assert moved > 0
            verdict = "turned round by D on every pair" if y_count % 2 == 0 else "NOT turned round, on any pair"
            print(f"   N = {N} {name:8s} Y letters per term {y_count}: {moved} moved pairs, step back = -step there, {verdict}")
    assert fw.cube_step('Y', 'X')[0] == 'Z' and fw.cube_step('Y', 'Z')[0] == 'X'
    assert fw.cube_step('Y', 'Y') is None and fw.cube_step('Y', 'I') is None
    print("   one bit, Y field: X -> Z, Z -> X (the face keeping k_Y), I and Y stand: the sign of Y never enters a path")


def conjugate_states():
    print("5. The conjugate states: R returns an X^N eigenstate up to the sign, D moves it, so its mirror image R(D.a) is")
    print("   its conjugate D.a up to sign and it is not mirror-fixed;")
    print("   the cat, the pair, |+>^N are mirror-fixed; I/2^N is not paired; a Y^N eigenstate's image is its Z^N shift")
    rng = np.random.default_rng(5)

    def eigenstate(P, val):
        w, V = np.linalg.eigh(P)
        cols = V[:, np.abs(w - val) < 1e-9]
        psi = cols @ (rng.normal(size=cols.shape[1]) + 1j * rng.normal(size=cols.shape[1]))
        return state(psi)

    def up_to_sign(u, v):
        return min(np.abs(u - v).max(), np.abs(u + v).max())

    def up_to_phase(u, v):
        k = np.argmax(np.abs(v))
        return np.abs(u - v * (u[k] / v[k])).max() if abs(v[k]) > 1e-12 else np.abs(u).max()

    for N in (2, 3):
        d = 2 ** N
        XN, YN = kron_all([X] * N), kron_all([Y] * N)
        print(f"   N = {N}:")
        for val in (+1, -1):
            a = coeffs(eigenstate(XN, val), N)
            Ra, Da = fw.cube_right_shift_vector(a, N, 'X'), fw.cube_transpose_vector(a, N)
            word_pairing = max(abs(abs(Ra[k]) - abs(a[k])) for k in range(4 ** N))
            assert word_pairing < 1e-12, "paired word for word: |a_{R sigma}| = |a_sigma|"
            r_fixed = up_to_sign(Ra, a)
            assert r_fixed < 1e-12, "the first move returns the state up to the sign (rho X^N = +-rho)"
            mirror_image = fw.cube_right_shift_vector(Da, N, 'X')                   # Pi = R after D
            assert up_to_sign(mirror_image, Da) < 1e-12, "the mirror image R(D a) is the conjugate D a up to sign"
            d_moves = up_to_sign(Da, a)
            assert d_moves > 1e-2, f"the second move does not leave the state alone ({d_moves:.2f}), so it is not mirror-fixed"
            assert up_to_sign(mirror_image, a) > 1e-2, "not mirror-fixed"
            print(f"      X^N = {val:+d} eigenstate: word pairing {word_pairing:.1e}, R.a vs +-a {r_fixed:.1e}, D.a vs +-a {d_moves:.2f}, R(D.a) vs +-D.a {up_to_sign(mirror_image, Da):.1e}")
        ghz = np.zeros(d); ghz[0] = ghz[-1] = 1
        plus = kron_all([np.array([[1], [1]])] * N).ravel()
        fixed_cases = [("cat", ghz), ("|+>^N", plus)]
        if N == 2:
            fixed_cases += [("singlet", np.array([0, 1, -1, 0])), ("triplet", np.array([0, 1, 1, 0]))]
        for name, psi in fixed_cases:
            a = coeffs(state(psi), N)
            Da = fw.cube_transpose_vector(a, N)
            assert up_to_sign(Da, a) < 1e-12 and up_to_sign(fw.cube_right_shift_vector(Da, N, 'X'), a) < 1e-12, f"{name}: both moves return it, mirror-fixed"
        mixed = coeffs(np.eye(d) / d, N)
        iXN = next(k for k in range(4 ** N) if _pauli_label(k, N) == 'X' * N)
        assert mixed[0] == 1 / d and mixed[iXN] == 0.0, "I/2^N carries only I^N and nothing on X^N, where R would need it"
        print(f"      {', '.join(n for n, _ in fixed_cases)} mirror-fixed; I/2^N: amount on I^N {mixed[0].real:.3f}, on X^N {mixed[iXN].real:.1f}")
        for val in (+1, -1):
            a = coeffs(eigenstate(YN, val), N)
            Ra, Za = fw.cube_right_shift_vector(a, N, 'X'), fw.cube_right_shift_vector(a, N, 'Z')
            assert up_to_phase(Ra, Za) < 1e-12, "R.a = phase * (Z^N shift of a), word for word"
            Da = fw.cube_transpose_vector(a, N)
            image = fw.cube_right_shift_vector(Da, N, 'X')
            assert up_to_phase(image, fw.cube_right_shift_vector(Da, N, 'Z')) < 1e-12, "the mirror image R(D.a) is the Z^N shift of the conjugate D.a"
            W1, W2 = fw.cube_height_weights(a, N), fw.cube_height_weights(image, N)
            assert np.abs(W1 - W2).max() < 1e-12, "the Z^N shift keeps every height, so the image's profile is the state's: W_k = W_{N-k} follows"
            print(f"      Y^N = {val:+d} eigenstate: R.a vs phase * Z^N-shift {up_to_phase(Ra, Za):.1e}, R(D.a) vs Z^N-shift of D.a {up_to_phase(image, fw.cube_right_shift_vector(Da, N, 'Z')):.1e}, heights of the image equal {np.abs(W1 - W2).max():.1e}")


def turned_products():
    print("6. The turn about Z keeps every height: a product of equator bits at any angle has the binomial profile;")
    print("   XX + YY and ZZ are carried to themselves by the turn, XX alone is not")
    for N in (2, 3, 4):
        for phi in (0.7, 2.5):
            psi = kron_all([np.array([[1], [np.exp(1j * phi)]])] * N).ravel()
            W = fw.cube_height_weights(coeffs(state(psi), N), N)
            assert np.abs(W - np.array([comb(N, k) for k in range(N + 1)]) / 2 ** N).max() < 1e-12
        print(f"   N = {N}: products at angles 0.7 and 2.5 have W = binomial/2^N")
    t = sp.symbols('t', real=True)
    Xs, Ys, Zs = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])
    S = sp.diag(1, sp.exp(sp.I * t))
    turn = lambda M: sp.simplify((S * M * S.H))                                 # the turn about Z by t on one site
    kron = lambda A, B: sp.Matrix(sp.BlockMatrix([[A[i, j] * B for j in range(2)] for i in range(2)]).as_explicit())
    def turn2(M):
        S2 = kron(S, S)
        return sp.simplify(S2 * M * S2.H)
    xx, yy, zz = kron(Xs, Xs), kron(Ys, Ys), kron(Zs, Zs)
    assert turn2(xx + yy) == xx + yy and turn2(zz) == zz, "XX + YY and ZZ are fixed by the turn"
    assert sp.simplify(turn2(xx)[0, 3] - sp.exp(-2 * sp.I * t)) == 0 and xx[0, 3] == 1, "XX alone is not: its corner entry turns by e^(-2it)"
    assert turn(Zs) == Zs and sp.simplify(turn(Xs)[0, 1] - sp.exp(-sp.I * t)) == 0, "Z is fixed, X turns"
    print("   one bond: turn(XX + YY) = XX + YY, turn(ZZ) = ZZ, turn(XX) != XX; the turn keeps k_Z (Z and I fixed, X and Y mixed)")


def field_along_the_watching():
    print("7. A field along the watching: on one bit it steps X <-> Y inside height 1 and leaves I and Z standing; every")
    print("   state has weight 1/2^N on I^N; on two bits the field's moves on |01> +- |10> cancel, on Bell+ they do not")
    assert fw.cube_step('Z', 'X') == ('Y', (0, 1, -1)) and fw.cube_step('Z', 'Y') == ('X', (0, -1, 1)), "a face diagonal: k_Z kept"
    assert fw.cube_step('Z', 'I') is None and fw.cube_step('Z', 'Z') is None
    turning = {s for s in 'IXYZ' if fw.cube_step('Z', s) is not None}
    assert turning == {'X', 'Y'} and {fw.cube_coords(s)[0] for s in turning} == {1} and {fw.cube_coords(s)[0] for s in 'IXYZ' if s not in turning} == {0}, "height 1 holds only turning words, height 0 only standing ones"
    print("   one bit: X -> Y, Y -> X at height 1; I and Z stand at height 0; the word I carries (Tr rho)^2/2 = 1/2 in every state")
    N = 2
    index = {_pauli_label(k, N): k for k in range(16)}

    def field_moves(a):
        out = np.zeros(16, dtype=complex)
        for k in range(16):
            if abs(a[k]) < 1e-12:
                continue
            for term in ('ZI', 'IZ'):
                r = fw.commutator_coefficient(term, _pauli_label(k, N))
                if r is not None:
                    out[index[r[0]]] += r[1] * a[k]
        return np.abs(out).max()
    index_of = lambda lab: next(k for k in range(16) if _pauli_label(k, N) == lab)
    bond = np.zeros(16, dtype=complex)
    for lab in ('XX', 'YY'):
        bond[index_of(lab)] = 1.0
    bond[index_of('ZZ')] = 0.5
    assert field_moves(bond) == 0.0, "the field's moves on the bond's own words XX + YY and ZZ cancel too, so the cancellation holds along the run"
    print("   two bits: the field's moves on the bond's words XX + YY, ZZ cancel exactly")
    for name, psi, still in (("singlet |01>-|10>", [0, 1, -1, 0], True), ("triplet |01>+|10>", [0, 1, 1, 0], True), ("Bell+ |00>+|11>", [1, 0, 0, 1], False)):
        a = coeffs(state(psi), N)
        m = field_moves(a)
        assert (m == 0.0) == still, f"{name}: the field's moves sum to {m:.1e} (exact quarters, so exactly zero or not)"
        print(f"   two bits, {name}: words {sorted(words(a, N))}, field moves summed {m:.1e}")


def field_across_the_watching():
    print("8. A field across the watching, along X: on one bit it steps Y <-> Z between the heights and leaves I and X")
    print("   standing; R pairs I with X and Y with Z: standing with standing, turning with turning")
    assert fw.cube_step('X', 'Y') == ('Z', (-1, 0, 1)) and fw.cube_step('X', 'Z') == ('Y', (1, 0, -1))
    assert fw.cube_step('X', 'I') is None and fw.cube_step('X', 'X') is None
    for s, t in (('I', 'X'), ('Y', 'Z')):
        assert fw.cube_right_shift(s, 'X')[0] == t and fw.cube_right_shift(t, 'X')[0] == s
    standing = {s for s in 'IXYZ' if fw.cube_step('X', s) is None}
    assert standing == {'I', 'X'} and {fw.cube_right_shift(s, 'X')[0] for s in standing} == standing, "R keeps the standing pair and the turning pair"
    print("   Y -> Z (height 1 -> 0), Z -> Y; I and X stand; R: I <-> X, Y <-> Z. Which weights balance (x^2 = 1, y^2 = z^2)")
    print("   is the matrix route's, the_wave_between.py section 7; the letters say which words are paired")


if __name__ == "__main__":
    letters()
    one_bit_watched()
    several_bits_watched()
    second_move_invisible()
    conjugate_states()
    turned_products()
    field_along_the_watching()
    field_across_the_watching()
