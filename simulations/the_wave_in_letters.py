"""The wave's look at itself, read in letters.

The letter-cube route behind reflections/ON_THE_WAVE_RECOGNIZING_THAT_IT_IS_NOT_ALONE.md, section by
section in the page's order (docs/THE_ONE_SQUARE.md §7 to §9). Nothing here propagates: no e^{tL}, no
eigensolver. A state is read as the amounts it puts on words (Pauli strings), each word at a height
k_Z. The script works with the coefficients c of ρ = Σ c·σ (fw.pauli_basis_vector); the page's amount
on a word is √(2^N)·c, so that its square is the weight 2^N·c² (fw.cube_height_weights), the word of
all I letters weighing 1/2^N in every state. Every sentence of the page that is written in letters is
a statement about letters, words,
heights, the steps of a field or a coupling, and the two moves of the mirror, R (ρ ↦ ρ·X^N, the height
flip k ↦ N − k) and D (the transpose, the sign (−1)^{n_Y}, no word moves). The matrix script
simulations/the_wave_between.py is the second route of the handshake.

1. One bit in four letters ("Alone, nothing can be seen"): for a pure state at polar angle θ and phase
   φ the amount on I is 1/2 in every state, on Z the split cos θ / 2 between the two levels, on X and Y
   sin θ·(cos φ, sin φ)/2, so the angle φ between the two parts divides X from Y; a phase of the whole
   wave changes no amount; for "Known only through encounters", a probe along one letter reads that
   letter's amount and nothing else (the letters are orthogonal) and the four amounts rebuild the state
   (sympy, exact).
2. One bit, watched, with a field ("What you cannot see in yourself"): the generator written in the
   letters from the letters' heights and the field's steps alone, checked entry by entry against the
   matrices; the watching fades X and Y at one rate 2γ and leaves I and Z alone; a field along Z moves
   X onto Y and back and the turn about Z commutes with the generator; a field along X moves Y onto Z
   and back, leaves X at rest, and the turn does not commute (sympy, exact).
3. Words ("Everything can sit between"): the pair carries II, XX, YY and ZZ, a quarter of its weight on
   each, YY with a minus sign, and nothing on a word with one letter other than I; such words number
   3N and the others 4^N − 1 − 3N (6 and 9 at N = 2, 15 and 1008 at N = 5); sympy, exact.
4. One rule for the steps ("What happened between is one number", and the look blind to the second
   move in "How far the look at itself reaches"): a word's own transpose is its Y sign (checked against
   the word's matrix), so a matrix unchanged by the transpose carries no word with an odd Y count; the
   step back is exactly minus the step there for every Hermitian term, so −H walks the paths of H the
   other way (the page's "cannot tell a coupling run forwards from the same coupling run backwards";
   the look's side of it is checked in section 9 of the_wave_between.py); a term with an even Y count joins
   two words of opposite Y parity, so flipping every Y sign reverses each of its steps (D·L_H·D = −L_H;
   in the string basis this is the transpose of L_H's matrix, which the repository writes L_H†, while
   F119's transpose lives in the computational basis, where a real H gives L_H back); a term with an odd
   Y count, the Y field and the twisted coupling's XY and YX, joins words of equal parity and D leaves
   each of its steps as it is, so where every term is of that kind D commutes with the generator, a copy
   (commutator_coefficient, exact, N = 1..3); one bit under the Y field: X → Z,
   Z → X, I and Y at rest; the loss σ⁻: off the heights its only entry is the step I → Z, and nothing
   steps from Z back to I (framework.cube_tables, exact).
5. The mirror ("The mirror meets itself backwards"): the letters' coordinates, the heights, the two
   edges; R flips every height k → N − k with a unit phase, and the watching's charge 2γ·k_Z is checked
   against the word's own matrix (Σ_l Z_l σ Z_l = (N − 2k_Z)·σ); a word and its image are charged
   together 2γN, the charge of the full height (integer γ, exact); R carries X^N to I^N with phase 1
   and D leaves X^N alone, so the image of a state carries on I^N what the state carried on X^N; the
   image's weights stand at the opposite heights (sorted, exact); R keeps the steps of a term exactly
   when its Y count plus its Z count is even, so with item 4 the mirror reflects exactly the terms with
   an even Y count and an even Z count (XX, YY, ZZ, a field along X; not a field along Z or Y, not XY,
   YX, XZ, YZ); under watching alone the image's self-encounter is then e^{−2γNt}·S(−t), term by
   term, the algebra e^{−2γ(N−k)t} = e^{−2γNt}·e^{2γkt} on the weights the sorted check moved.
6. One bit, watched: W_0 = (1 + z²)/2, W_1 = (x² + y²)/2 exactly (sympy); equal weights with the ball
   give z = 0 and x² + y² = 1, the equator.
7. Several bits, watched: the height profiles of |+⟩^N (binomial), the cat and the pair (half on each
   edge), a random product of equator bits (symmetric), the W state and a tilted product (not; the W
   state at N = 3, 4, since at N = 2 it is the pair |01⟩ + |10⟩ and pairs); W_k = W_{N−k} is the twin
   condition under watching alone.
8. The conjugate states: a random X^N eigenstate (either sign) is returned by R up to the sign
   (R·a = ± a, so |a_{Rσ}| = |a_σ| word for word), while D moves it (D·a ≠ ± a), so its mirror image
   Π·a = R·D·a is its conjugate D·a up to sign and it is not mirror-fixed; the cat, the pair and |+⟩^N
   (and at N = 2 the singlet and the triplet) are returned by both moves, mirror-fixed; the completely
   mixed state carries only I^N and nothing on X^N, where R would need it. A Y^N eigenstate:
   R·a = phase · (Z^N-shift of a) word for word, and the mirror image R(D·a) = phase · (Z^N-shift of
   D·a), the same heights as a.
9. The turn about Z keeps every height: the height profile of a product of equator bits at any angle is
   the binomial one; XX + YY and ZZ are carried to themselves by the turn, XX alone is not (symbolic,
   one bond).
10. A field along the watching: on one bit it moves X ↔ Y inside height 1 and leaves I and Z at rest
    (cube_step); every normalized state has weight 1/2^N on I^N (so the heights cannot pair when height
    1 holds only moving words); on two bits the field's moves on the words of |01⟩ ± |10⟩ cancel
    (commutator coefficients summed), and so do its moves on the bond's own words XX and YY, while ZZ it
    does not move, so the cancellation holds along the run; on Bell+ they do not.
11. A field across the watching, along X: on one bit it moves Y ↔ Z between the heights and leaves I
    and X at rest; R pairs I with X and Y with Z, at rest with at rest and moving with moving. Which
    weights must balance (x² = 1, y² = z²) is the matrix route's (section 7 of the_wave_between.py); the
    letters give the pairs.
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
from framework import cube_tables  # noqa: E402
from framework.pauli import _pauli_label  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SYM = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]), 'Z': sp.diag(1, -1)}


def kron_all(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


def sym_word(label):
    """The word's own matrix in exact arithmetic."""
    out = sp.Matrix([[1]])
    for c in label:
        out = sp.kronecker_product(out, SYM[c])
    return out


def labels(N):
    return [''.join(t) for t in product('IXYZ', repeat=N)]


LETTER = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def mat(label):
    """The word's own matrix, for the checks against the matrices only."""
    return kron_all([LETTER[c] for c in label])


def coeffs(rho, N):
    return fw.pauli_basis_vector(rho, N)


def words(c, N, tol=1e-12):
    """The words a coefficient vector carries, label -> coefficient."""
    return {_pauli_label(k, N): c[k] for k in range(4 ** N) if abs(c[k]) > tol}


def state(psi):
    psi = np.asarray(psi, dtype=complex)
    return np.outer(psi, psi.conj()) / np.vdot(psi, psi).real


def placed(template, N):
    return ['I' * i + template + 'I' * (N - len(template) - i) for i in range(N - len(template) + 1)]


def one_bit_in_letters():
    print("1. One bit in four letters: I is one half in every state, Z the split between the two levels, X and Y the angle")
    print("   between the two parts; a phase of the whole wave is in no letter")
    theta, phi, alpha = sp.symbols('theta phi alpha', real=True)
    psi = sp.Matrix([sp.cos(theta / 2), sp.exp(sp.I * phi) * sp.sin(theta / 2)])
    rho = psi * psi.H
    c = {L: sp.simplify((rho * SYM[L]).trace() / 2) for L in 'IXYZ'}         # rho = sum_P c_P P
    assert c['I'] == sp.Rational(1, 2), "I: the same amount in every state, that it is there"
    assert sp.simplify(c['Z'] - sp.cos(theta) / 2) == 0, "Z: cos(theta)/2, the split p_0 - p_1 between the two levels"
    assert sp.simplify(c['X'] - sp.sin(theta) * sp.cos(phi) / 2) == 0 and sp.simplify(c['Y'] - sp.sin(theta) * sp.sin(phi) / 2) == 0, "X, Y: the angle phi between the two parts divides X from Y"
    whole = sp.exp(sp.I * alpha) * psi
    assert sp.simplify(whole * whole.H - rho) == sp.zeros(2, 2), "a phase of the whole wave changes no amount"
    assert [fw.cube_coords(L)[0] for L in 'IXYZ'] == [0, 1, 1, 0], "I and Z along the watching (height 0), X and Y across it (height 1)"
    for A in 'IXYZ':
        for B in 'IXYZ':
            assert (SYM[A] * SYM[B]).trace() / 2 == (1 if A == B else 0), "a probe along one letter reads that letter's amount and nothing else"
    rebuilt = sum((c[L] * SYM[L] for L in 'IXYZ'), sp.zeros(2, 2)) - rho
    assert rebuilt.applyfunc(lambda e: sp.simplify(sp.expand_complex(e))) == sp.zeros(2, 2), "the four amounts are the whole state"
    page = {L: sp.simplify(sp.sqrt(2) * c[L]) for L in 'IXYZ'}                # the page's amounts, sqrt(2^N) c
    assert page['I'] == 1 / sp.sqrt(2) and page['I'] ** 2 == sp.Rational(1, 2), "the word I weighs one half in every state"
    print(f"   coefficients of rho = sum c P: I {c['I']}, Z {c['Z']}, X {c['X']}, Y {c['Y']}; the whole wave's phase alpha appears in none;")
    print(f"   the page's amounts are sqrt(2) times these, I {page['I']}, so that the word I weighs {page['I'] ** 2};")
    print("   a probe along one letter reads that letter's amount alone, and the four amounts rebuild the state exactly")


def watching_and_two_fields():
    print("2. One bit, watched, with a field: the watching fades X and Y at one rate and leaves I and Z alone; a field along Z")
    print("   moves X onto Y and back, and the turn about Z is a move the dynamics do not notice; a field along X moves Y onto Z")
    print("   and back and leaves X at rest, and the turn is noticed")
    g, h, t = sp.symbols('gamma h t', positive=True)
    order = 'IXYZ'

    def in_letters(field):
        """The generator in the letters (I, X, Y, Z), from the heights and the field's steps alone."""
        M = sp.zeros(4, 4)
        for b, s in enumerate(order):
            M[b, b] -= 2 * g * fw.cube_coords(s)[0]
            r = fw.commutator_coefficient(field, s)
            if r is not None:
                tau, cc = r
                M[order.index(tau), b] += -sp.I * h * sp.nsimplify(cc)
        return M

    def from_matrices(field):
        """The same generator from the matrices: L(P) = -i h [F, P] + gamma (Z P Z - P), M[a, b] = Tr(P_a L(P_b))/2."""
        F = SYM[field]
        M = sp.zeros(4, 4)
        for b, s in enumerate(order):
            P = SYM[s]
            LP = -sp.I * h * (F * P - P * F) + g * (SYM['Z'] * P * SYM['Z'] - P)
            for a, r in enumerate(order):
                M[a, b] = sp.simplify((SYM[r] * LP).trace() / 2)
        return M

    turn = sp.Matrix([[1, 0, 0, 0], [0, sp.cos(t), -sp.sin(t), 0], [0, sp.sin(t), sp.cos(t), 0], [0, 0, 0, 1]])
    for field, noticed in (('Z', False), ('X', True)):
        M = in_letters(field)
        assert sp.simplify(M - from_matrices(field)) == sp.zeros(4, 4), f"field along {field}: the letters' generator is the matrices' generator"
        moved = {s: fw.cube_step(field, s)[0] for s in order if fw.cube_step(field, s) is not None}
        comm = sp.simplify(turn * M - M * turn)
        assert (comm != sp.zeros(4, 4)) == noticed, f"field along {field}: the turn about Z is {'noticed' if noticed else 'not noticed'}"
        print(f"   field along {field}: moves {moved}, at rest {[s for s in order if s not in moved]}; the turn about Z {'does NOT commute' if noticed else 'commutes'} with the generator")
    Mz = in_letters('Z')
    assert [Mz[k, k] for k in range(4)] == [0, -2 * g, -2 * g, 0], "the watching fades X and Y at 2 gamma, I and Z not at all"
    print(f"   the watching alone, on the diagonal: I {Mz[0, 0]}, X {Mz[1, 1]}, Y {Mz[2, 2]}, Z {Mz[3, 3]}")


def words_of_several_bits():
    print("3. Words: the pair carries II, XX, YY, ZZ, a quarter of its weight on each, and no word with a letter other than I")
    print("   at one bit only; words with one such letter number 3N, the others 4^N - 1 - 3N")
    N = 2
    psi = sp.Matrix([1, 0, 0, 1]) / sp.sqrt(2)
    rho = psi * psi.H
    amount = {s: sp.simplify((rho * sym_word(s)).trace() / 2 ** N) for s in labels(N)}
    carried = {s: a for s, a in amount.items() if a != 0}
    assert carried == {'II': sp.Rational(1, 4), 'XX': sp.Rational(1, 4), 'YY': -sp.Rational(1, 4), 'ZZ': sp.Rational(1, 4)}
    assert all(2 ** N * a ** 2 == sp.Rational(1, 4) for a in carried.values()), "a quarter of the weight on each word"
    assert not any(fw.cube_weight(s) == 1 for s in carried), "nothing on a word with a letter other than I at one bit only"
    print(f"   the pair: {carried}, weight 1/4 on each")
    for N in (1, 2, 3, 4, 5):
        one = sum(1 for s in labels(N) if fw.cube_weight(s) == 1)
        rest = sum(1 for s in labels(N) if fw.cube_weight(s) >= 2)
        assert one == 3 * N and rest == 4 ** N - 1 - 3 * N
        print(f"   N = {N}: words about one bit alone {one}, about two or more {rest}")


def one_rule_for_the_steps():
    print("4. One rule for the steps: the step back is minus the step there; a term with an even Y count joins words of")
    print("   opposite Y parity, so flipping every Y sign reverses its steps; an odd Y count (the Y field, the twisted coupling's")
    print("   XY and YX) joins words of equal parity; the loss steps from I to Z and nothing steps back")
    for N in (1, 2, 3):
        for s in labels(N):
            M = mat(s)
            assert np.array_equal(M.T, fw.cube_transpose_sign(s) * M), "a word's own transpose is its Y sign"
        terms = {"X field": placed("X", N), "Z field": placed("Z", N), "Y field": placed("Y", N)}
        if N > 1:
            terms.update({"XX": placed("XX", N), "YY": placed("YY", N), "ZZ": placed("ZZ", N), "XY": placed("XY", N), "YX": placed("YX", N)})
        for name, ts in terms.items():
            y_count = ts[0].count('Y')
            moved = 0
            for term in ts:
                for s in labels(N):
                    r = fw.commutator_coefficient(term, s)
                    if r is None:
                        continue
                    tau, c = r
                    back = fw.commutator_coefficient(term, tau)
                    assert back is not None and back[0] == s, "a term that moves sigma to tau moves tau back to sigma"
                    assert back[1] == -c, "the step back is minus the step there (Hermitian words: the phase is +-i)"
                    parity_flip = fw.cube_transpose_sign(s) * fw.cube_transpose_sign(tau) == -1
                    assert parity_flip == (y_count % 2 == 0), f"{name}: an even-Y term joins words of opposite Y parity, an odd-Y term words of equal parity"
                    moved += 1
            assert moved > 0
            verdict = "reversed by D on every pair" if y_count % 2 == 0 else "left as it is by D on every pair"
            print(f"   N = {N} {name:8s} Y letters per term {y_count}: {moved} moved pairs, step back = -step there, {verdict}")
    print("   so turning every term to its negative gives each step the value of its step back: -H walks the paths of H the")
    print("   other way, whatever the Y counts (that the look then cannot tell them apart under the watching is checked by")
    print("   the matrix route, the_wave_between.py section 9)")
    print("   so a matrix the transpose leaves alone carries no word with an odd Y count (every word's transpose is its Y sign)")
    assert fw.cube_step('Y', 'X')[0] == 'Z' and fw.cube_step('Y', 'Z')[0] == 'X'
    assert fw.cube_step('Y', 'Y') is None and fw.cube_step('Y', 'I') is None
    print("   one bit, Y field: X -> Z, Z -> X (the face keeping k_Y), I and Y at rest: the sign of Y never enters a path")
    g = sp.Symbol('gamma', positive=True)
    T = cube_tables.transfer_times_four(cube_tables.LOWER, g)            # 4 x the loss's transfer matrix, rows/columns I X Y Z
    off = {(a, b): T[a, b] for a in range(4) for b in range(4) if a != b and T[a, b] != 0}
    assert off == {(3, 0): 4 * g}, f"the loss: off the heights only the step I -> Z ({off})"
    assert [T[k, k] for k in range(4)] == [0, -2 * g, -2 * g, -4 * g], "its heights: X and Y at gamma/2, Z at gamma"
    print(f"   the loss sigma^-: 4 x its entries, heights {[T[k, k] for k in range(4)]}, step I -> Z {T[3, 0]}, step Z -> I {T[0, 3]}")


def the_mirror():
    print("5. The mirror: R flips every height k -> N - k; a word and its image are charged together the full height's charge;")
    print("   the image of a state carries on I^N what the state carried on X^N; R keeps the steps of a term exactly when its")
    print("   Y count plus its Z count is even; under watching alone the image meets itself backwards at the price")
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
            conj_sum = sum(mat(s[:l] + 'Z' + s[l + 1:]) @ M @ mat(s[:l] + 'Z' + s[l + 1:]) for l in range(N))
            assert np.array_equal(conj_sum, (N - 2 * fw.cube_coords(s)[0]) * M), "the watching: sum_l Z_l sigma Z_l = (N - 2 k_Z) sigma, so the word fades at 2 gamma k_Z and moves nowhere"
            assert fw.cube_rate(s, 0.3) == 2 * 0.3 * fw.cube_coords(s)[0]
            assert fw.cube_rate(s, 1) + fw.cube_rate(t, 1) == 2 * N, "a word and its image are charged together the full height's charge"
        assert {fw.cube_right_shift(s, 'X')[0] for s in quiet} == set(far), "R carries the quiet edge onto the far edge"
        assert fw.cube_right_shift('X' * N, 'X') == ('I' * N, 1) and fw.cube_transpose_sign('X' * N) == 1, "X^N -> I^N with phase 1, no Y sign"
        print(f"   N = {N}: {4 ** N} words, quiet edge {len(quiet)}, far edge {len(far)}, R quiet -> far, every word + its image charged 2N (gamma = 1)")
    rng = np.random.default_rng(7)
    for N in (2, 3):
        d = 2 ** N
        index = {_pauli_label(k, N): k for k in range(4 ** N)}
        for _ in range(3):
            psi = rng.normal(size=d) + 1j * rng.normal(size=d)
            a = coeffs(state(psi), N)
            image = fw.cube_right_shift_vector(fw.cube_transpose_vector(a, N), N, 'X')      # Pi = R after D
            assert image[index['I' * N]] == a[index['X' * N]], "the image carries on I^N exactly what the state carried on X^N"
            for k in range(N + 1):
                mine = sorted(abs(a[j]) ** 2 for j in range(4 ** N) if fw.cube_coords(_pauli_label(j, N))[0] == k)
                theirs = sorted(abs(image[j]) ** 2 for j in range(4 ** N) if fw.cube_coords(_pauli_label(j, N))[0] == N - k)
                assert mine == theirs, "the image's weights at height N - k are the state's at height k"
        print(f"   N = {N}: three random states; the image's amount on I^N is the state's on X^N, its weights stand at the opposite heights")

    def r_keeps(template, N):
        for term in placed(template, N):
            for s in labels(N):
                Rs, ps = fw.cube_right_shift(s, 'X')
                there, here = fw.commutator_coefficient(term, s), fw.commutator_coefficient(term, Rs)
                if there is None and here is None:
                    continue
                if (there is None) != (here is None):
                    return False
                tau, c = there
                Rtau, ptau = fw.cube_right_shift(tau, 'X')
                if here[0] != Rtau or ps * here[1] != c * ptau:
                    return False
        return True

    for template in ('XX', 'YY', 'ZZ', 'X', 'Z', 'Y', 'XY', 'YX', 'XZ', 'YZ'):
        ny, nz = template.count('Y'), template.count('Z')
        keeps = r_keeps(template, 3)
        assert keeps == ((ny + nz) % 2 == 0), template
        print(f"   N = 3 {template:2s}: R keeps its steps {keeps}, D reverses them {ny % 2 == 0} (section 4), the mirror reflects it {keeps and ny % 2 == 0}")
    # not a gate: with the weights moved as checked above, the backwards reading under watching alone is algebra
    print("   watching alone: the image's self-encounter sum_k W_k e^(-2 gamma (N - k) t) = e^(-2 gamma N t) S(-t), the weight")
    print("   of height k charged from the far edge (algebra on the weights the sorted check moved)")


def one_bit_watched():
    print("6. One bit, watched: W_0 = (1 + z^2)/2 (the word I alone is one half), W_1 = (x^2 + y^2)/2;")
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
    print("7. Several bits, watched: the twin condition is W_k = W_{N-k}; equator products are binomial, the pair and the")
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
            if name == "W" and N == 3:
                assert np.abs(W - np.array([1 / 3, 0, 2 / 3, 0])).max() < 1e-12, "the W state of three bits: heights 0 and 2, none at 1 or 3"
            print(f"      {name:24s} W = {np.round(W, 4)}  |W - reversed| {sym:.1e}")


def conjugate_states():
    print("8. The conjugate states: R returns an X^N eigenstate up to the sign, D moves it, so its mirror image R(D.a) is")
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
    print("9. The turn about Z keeps every height: a product of equator bits at any angle has the binomial profile;")
    print("   XX + YY and ZZ are carried to themselves by the turn, XX alone is not")
    for N in (2, 3, 4):
        for phi in (0.7, 2.5):
            psi = kron_all([np.array([[1], [np.exp(1j * phi)]])] * N).ravel()
            W = fw.cube_height_weights(coeffs(state(psi), N), N)
            assert np.abs(W - np.array([comb(N, k) for k in range(N + 1)]) / 2 ** N).max() < 1e-12
        print(f"   N = {N}: products at angles 0.7 and 2.5 have W = binomial/2^N")
    t = sp.symbols('t', real=True)
    Xs, Ys, Zs = SYM['X'], SYM['Y'], SYM['Z']
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
    print("10. A field along the watching: on one bit it moves X <-> Y inside height 1 and leaves I and Z at rest; every")
    print("    state has weight 1/2^N on I^N; on two bits the field's moves on |01> +- |10> cancel, on Bell+ they do not")
    assert fw.cube_step('Z', 'X') == ('Y', (0, 1, -1)) and fw.cube_step('Z', 'Y') == ('X', (0, -1, 1)), "a face diagonal: k_Z kept"
    assert fw.cube_step('Z', 'I') is None and fw.cube_step('Z', 'Z') is None
    moving = {s for s in 'IXYZ' if fw.cube_step('Z', s) is not None}
    assert moving == {'X', 'Y'} and {fw.cube_coords(s)[0] for s in moving} == {1} and {fw.cube_coords(s)[0] for s in 'IXYZ' if s not in moving} == {0}, "height 1 holds only moving words, height 0 only words at rest"
    rng = np.random.default_rng(10)
    for N_ in (1, 2, 3):
        for _ in range(3):
            psi = rng.normal(size=2 ** N_) + 1j * rng.normal(size=2 ** N_)
            c0 = coeffs(state(psi), N_)[0]                                   # the word of all I letters
            assert abs(2 ** N_ * abs(c0) ** 2 - 1 / 2 ** N_) < 1e-15, "the word of all I letters weighs 1/2^N in every state"
    print("    one bit: X -> Y, Y -> X at height 1; I and Z at rest at height 0; the word I weighs 1/2 in every state, I^N 1/2^N")
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
    bond = np.zeros(16, dtype=complex)
    for lab in ('XX', 'YY'):
        bond[index[lab]] = 1.0
    bond[index['ZZ']] = 0.5
    assert field_moves(bond) == 0.0, "the field's moves on the bond's own words XX and YY cancel too, and ZZ it does not move, so the cancellation holds along the run"
    assert all(fw.commutator_coefficient(term, 'ZZ') is None for term in ('ZI', 'IZ')), "the field does not move ZZ"
    print("    two bits: the field's moves on the bond's words XX and YY cancel exactly, and ZZ it does not move")
    for name, psi, still in (("singlet |01>-|10>", [0, 1, -1, 0], True), ("triplet |01>+|10>", [0, 1, 1, 0], True), ("Bell+ |00>+|11>", [1, 0, 0, 1], False)):
        a = coeffs(state(psi), N)
        m = field_moves(a)
        assert (m == 0.0) == still, f"{name}: the field's moves sum to {m:.1e} (exact quarters, so exactly zero or not)"
        print(f"    two bits, {name}: words {sorted(words(a, N))}, field moves summed {m:.1e}")


def field_across_the_watching():
    print("11. A field across the watching, along X: on one bit it moves Y <-> Z between the heights and leaves I and X")
    print("    at rest; R pairs I with X and Y with Z: at rest with at rest, moving with moving")
    assert fw.cube_step('X', 'Y') == ('Z', (-1, 0, 1)) and fw.cube_step('X', 'Z') == ('Y', (1, 0, -1))
    assert fw.cube_step('X', 'I') is None and fw.cube_step('X', 'X') is None
    for s, t in (('I', 'X'), ('Y', 'Z')):
        assert fw.cube_right_shift(s, 'X')[0] == t and fw.cube_right_shift(t, 'X')[0] == s
    at_rest = {s for s in 'IXYZ' if fw.cube_step('X', s) is None}
    assert at_rest == {'I', 'X'} and {fw.cube_right_shift(s, 'X')[0] for s in at_rest} == at_rest, "R keeps the pair at rest and the moving pair"
    print("    Y -> Z (height 1 -> 0), Z -> Y; I and X at rest; R: I <-> X, Y <-> Z. Which weights balance (x^2 = 1, y^2 = z^2)")
    print("    is the matrix route's, the_wave_between.py section 7; the letters say which words are paired")


if __name__ == "__main__":
    one_bit_in_letters()
    watching_and_two_fields()
    words_of_several_bits()
    one_rule_for_the_steps()
    the_mirror()
    one_bit_watched()
    several_bits_watched()
    conjugate_states()
    turned_products()
    field_along_the_watching()
    field_across_the_watching()
