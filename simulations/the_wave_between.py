"""What a state can see of itself, and what only a second can see.

The numbers behind reflections/ON_THE_WAVE_RECOGNIZING_THAT_IT_IS_NOT_ALONE.md.
Nothing here is new physics; the script lays known facts side by side so the
reflection can point at numbers instead of words.

1. The self-encounter S_A(t) = Tr(A† e^{tL} A) of a qubit cannot tell two
   states apart that differ only in their phase about Z, while an X probe can.
   The reason is exact: rotation about Z commutes with L (field along Z, Z
   dephasing), and the self-encounter is unchanged by any unitary symmetry of L.
   A field along X breaks that symmetry and the two self-encounters part.
2. The Bell pair: both one-site Bloch vectors are zero (each site alone is
   I/2), all content sits in the nine two-site coefficients.
3. The count of Pauli strings: 3N on one site, 4^N - 1 - 3N on two or more.
4. Reciprocity with content (F119 §7): for a real H and Z dephasing the
   generator equals its own transpose, L = L^T, exactly, so for symmetric A
   and B, Tr(B e^{tL} A) = Tr(A e^{tL} B) in the SAME system: preparation and
   readout trade places. A Dzyaloshinskii-Moriya term or amplitude damping
   breaks L = L^T, and then the swapped pairing differs. (The bare identity
   Tr(B e^{tL} A) = Tr(A e^{tL^T} B) for symmetric A, B holds for every L and
   is not tested here, since it cannot fail.)
5. The mirror: in F1's class (here an XXZ pair under Z dephasing) Π L Π⁻¹ =
   -L - 2σ with Π unitary, so for any operator A (ΠA is an operator, not a
   state) S_{ΠA}(t) = e^{-2σt} S_A(-t), backwards in the formal sense only.
6. The set S_A = S_{ΠA} for one qubit under dephasing alone (H = 0): exactly
   the pure states on the equator. Derived in sympy: S_A - S_{ΠA} factors as
   (1 - e^{-2γt})(1 + c² - a² - b²)/2, so the condition is a² + b² = 1 + c², and
   inside the ball a² + b² + c² <= 1 that forces c = 0 and a² + b² = 1;
   off-equator states as numeric controls.
7. One qubit with a field of its own, exact (sympy): the condition that
   e^{γt} S_A(t) is even in t, which is S_A = S_{MA} for the mirror M of the
   system wherever one exists. A field across the watching, along X, keeps Π
   as the mirror and leaves exactly |+⟩ and |−⟩ among all states, pure or
   mixed; along Y the mirror is Π turned by a quarter about Z and the answer
   is |+i⟩ and |−i⟩; at any other angle across, the turn about Z by that
   angle carries both, since it commutes with the watching (checked
   numerically at one angle). M must be unitary for S_{MA}(t) = e^{−2σt}S_A(−t);
   Π and its turns are. A field along the watching leaves no state: its
   spectrum {0, 0, −2γ ± 2ih} pairs about no centre (the real parts allow
   only −γ, and there the two non-turning modes at 0 lack non-turning
   partners at −2γ; L is real, so the same holds for −L† − c and for
   antilinear maps), so no map of any form reflects the generator at any
   price. Π still maps every pure equator state to an operator with the
   same self-encounter, but that is no mirror image.
8. Several bits watched at one rate. L = −2γ·k_Z on Pauli strings (k_Z = the
   number of X or Y letters), the mirror reflects k_Z to N − k_Z, so the
   condition is that the weight of the state sits symmetrically about the
   middle height N/2 (the only-if direction because distinct rates give
   linearly independent exponentials): a product of equatorial bits has it (binomial profile),
   the pair and the cat state have it with half on each edge, the W state
   does not. A state the mirror fixes up to sign, ΠA = ±A (the pair, the cat,
   every product of |+⟩ and |−⟩ bits), keeps it under every coupling that Π
   itself reflects, because L + σ maps each Π eigenspace to the opposite one
   and the odd moments vanish (the condition is that Π reflects the coupling,
   not that its spectrum pairs; couplings paired by another element are not
   tested here). A product of bits all at one other equator angle is the
   global Z-turn of |+⟩^N, and the turned Π reflects every U(1)-symmetric
   coupling Π reflects, so it keeps the property under XXZ at every angle; at
   the angles 0.7 and 2.5 it loses it under XX bonds alone, which break U(1)
   (at the quarter turns the turned Π reflects XX bonds too, so it would be
   kept there; not tested here).
   One seeded draw of random angles for each N loses it under both. N = 2, 3,
   4; asserted: the symmetry, the fixedness, the cat's halves, the gaps of the
   fixed states, of the one-angle products and of the random draw, and the
   watching-alone gaps of W and the tilted product; the coupled gaps of W, the
   tilted product and the random pure state are printed.
9. Coupled bits with a field, and the conjugate pair. A field along X, the
   same on every bit, is real and commutes with X^N, so Π reflects the full
   generator exactly when it reflected the coupling (XXZ, XX only). A random
   eigenstate of X^N (either sign) is NOT Π-fixed, but its mirror image is its
   complex conjugate up to sign (ρ X^N = ±ρ, then the transpose), and the
   conjugate meets itself alike because L = Lᵀ for real H; so it keeps the
   property with and without the field. Controls: a Y^N eigenstate under the
   X field fails, and so does the completely mixed state (unchanged by X^N
   under conjugation, not an eigenstate). Y^N eigenstates without a field keep
   it too (checked, both signs, both couplings). A field at another angle in
   the equator plane turns the string and the mirror with it: the turned Π
   reflects XXZ (U(1)) and not XX only; under XXZ the turned string's
   eigenstates keep it and the cat, an X^N eigenstate, does not. A random pure
   state and its conjugate meet themselves alike under the X field and part
   under a Y field once bits are coupled (L ≠ Lᵀ there, necessary and not by
   itself sufficient: one bit alone under a Y field keeps the pair, its Y
   coefficient only decays; asserted). Every two-site Pauli coupling that Π
   reflects is real (asserted), so the conjugate route covers every coupling
   that keeps the mirror. A field along the watching, equal on both bits
   of a pair: Π broken, spectrum paired about no centre, yet the singlet and
   the triplet |01⟩ ± |10⟩, which the field does not turn, keep it exactly;
   unequal fields on the two bits break it, and Bell+ fails. N = 2, 3 (the Z
   field at N = 2); all asserted.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding="utf-8")
from framework.lindblad import lindbladian_general, lindbladian_z_dephasing  # noqa: E402
from framework.symmetry import build_pi_full  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X, Y, Z]


def vec(M):
    return M.reshape(-1)


def kron_all(mats):
    out = np.array([[1.0 + 0j]])
    for m in mats:
        out = np.kron(out, m)
    return out


def unvec(v, d):
    return v.reshape(d, d)


def self_encounter_blind_to_phase():
    print("1. The self-encounter cannot see the phase; a second one can")
    omega, gamma = 1.0, 0.3
    L = lindbladian_z_dephasing(0.5 * omega * Z, [gamma])

    def bloch_state(theta, phi):
        a = (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta))
        return 0.5 * (I2 + a[0] * X + a[1] * Y + a[2] * Z)

    A1, A2 = bloch_state(1.0, 0.0), bloch_state(1.0, 2.0)
    print("   t    S(A1)      S(A2)      X probe on A1  X probe on A2")
    for t in (0.0, 0.5, 1.0, 2.0, 4.0):
        U = expm(L * t)
        r1, r2 = unvec(U @ vec(A1), 2), unvec(U @ vec(A2), 2)
        print(f"   {t:<4} {np.trace(A1 @ r1).real:.6f}   {np.trace(A2 @ r2).real:.6f}"
              f"   {np.trace(X @ r1).real:+.6f}      {np.trace(X @ r2).real:+.6f}")
    Lx = lindbladian_z_dephasing(0.5 * omega * X, [gamma])
    print("   control, field along X (no longer a symmetry of the phase):")
    for t in (0.5, 2.0):
        U = expm(Lx * t)
        r1, r2 = unvec(U @ vec(A1), 2), unvec(U @ vec(A2), 2)
        print(f"   {t:<4} {np.trace(A1 @ r1).real:.6f}   {np.trace(A2 @ r2).real:.6f}")


def bell_pair_is_all_between():
    print("2. The Bell pair: nothing on either site, everything between")
    b = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    R = np.outer(b, b.conj())

    def coeff(i, j):
        return np.trace(R @ np.kron(PAULI[i], PAULI[j])).real

    print("   A alone:", [float(round(coeff(i, 0), 6)) for i in (1, 2, 3)])
    print("   B alone:", [float(round(coeff(0, j), 6)) for j in (1, 2, 3)])
    print("   between:", [[float(round(coeff(i, j), 6)) for j in (1, 2, 3)] for i in (1, 2, 3)])


def count_selves_and_between():
    print("3. Pauli strings on one site against strings on two or more")
    for N in range(1, 6):
        print(f"   N = {N}: one site {3 * N}, between {4 ** N - 1 - 3 * N}")


def two_site_generator(gammas):
    H = 0.5 * (np.kron(X, X) + np.kron(Y, Y)) + 0.35 * np.kron(Z, Z)
    return lindbladian_z_dephasing(H, gammas)


def reciprocity():
    print("4. Reciprocity: preparation and readout trade places in the same system")
    rng = np.random.default_rng(1)
    Ma = rng.normal(size=(4, 4))
    Ma = Ma + Ma.T
    Mb = rng.normal(size=(4, 4))
    Mb = Mb + Mb.T
    H = 0.5 * (np.kron(X, X) + np.kron(Y, Y)) + 0.35 * np.kron(Z, Z)
    dm = 0.3 * (np.kron(X, Y) - np.kron(Y, X))
    sm = np.array([[0, 1], [0, 0]], dtype=complex)
    loss = 0.2 * (lindbladian_general(np.zeros((4, 4)), [np.kron(sm, I2), np.kron(I2, sm)]))
    cases = (("real H, Z dephasing", two_site_generator([0.2, 0.2])),
             ("with a DM term", lindbladian_z_dephasing(H + dm, [0.2, 0.2])),
             ("with amplitude damping", two_site_generator([0.2, 0.2]) + loss))
    for name, L in cases:
        U = expm(L * 0.7)
        forward = np.trace(Mb @ unvec(U @ vec(Ma), 4)).real
        swapped = np.trace(Ma @ unvec(U @ vec(Mb), 4)).real
        print(f"   {name}: max|L - L^T| = {np.abs(L - L.T).max()}, "
              f"Tr(B e^(tL) A) = {forward:.12f}, Tr(A e^(tL) B) = {swapped:.12f}")


def mirror_self_encounter():
    print("5. The mirror image meets itself as the original running backwards")
    gammas = [0.3, 0.5]
    sigma = sum(gammas)
    L = two_site_generator(gammas)
    strings = [np.kron(PAULI[a], PAULI[b]) for a, b in itertools.product(range(4), repeat=2)]
    T = np.array([vec(S) for S in strings]).T / 2
    LP = T.conj().T @ L @ T
    Pi = build_pi_full(2)
    residual = np.abs(Pi @ LP @ Pi.conj().T + LP + 2 * sigma * np.eye(16)).max()
    print(f"   palindrome residual in this basis: {residual}")
    rng = np.random.default_rng(3)
    a = rng.normal(size=16) + 1j * rng.normal(size=16)
    b = Pi @ a
    for t in (0.4, 1.3):
        s_mirror = b.conj() @ expm(LP * t) @ b
        s_back = np.exp(-2 * sigma * t) * (a.conj() @ expm(-LP * t) @ a)
        print(f"   t = {t}: S_PiA(t) = {s_mirror:.12f}, e^(-2σt) S_A(-t) = {s_back:.12f}")


def mirror_twin_exists():
    print("6. States whose self-encounter equals their mirror image's (H = 0, one qubit)")
    print("   closed form: S_A = S_PiA at every t iff a^2 + b^2 = 1 + c^2 for A = (I + aX + bY + cZ)/2,")
    print("   which among states means exactly the pure states on the equator")
    print("   (for gamma > 0; at gamma = 0 the generator vanishes and every state qualifies)")
    import sympy as sp
    g, t = sp.symbols("gamma t", positive=True)
    x, y, z = sp.symbols("x y z", real=True)
    letters = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)]
    rates = []
    for P in letters:   # the dephasing generator gamma (Z P Z - P) scales every letter by its own rate
        DP = g * (letters[3] * P * letters[3] - P)
        r = sp.simplify((DP * P).trace() / 2)
        assert sp.simplify(DP - r * P) == sp.zeros(2, 2)
        rates.append(r)
    Pi_num = build_pi_full(1)
    Pi_exact = sp.Matrix(4, 4, lambda i, j: sp.nsimplify(Pi_num[i, j].real) + sp.I * sp.nsimplify(Pi_num[i, j].imag))

    def exact_encounter(coeffs):
        A = sum((c * P for c, P in zip(coeffs, letters)), sp.zeros(2, 2))
        At = sum((c * sp.exp(r * t) * P for c, r, P in zip(coeffs, rates, letters)), sp.zeros(2, 2))
        return sp.expand((A.H * At).trace())

    a_exact = sp.Matrix([1, x, y, z]) / 2
    diff = sp.simplify(exact_encounter(list(a_exact)) - exact_encounter(list(Pi_exact * a_exact)))
    claimed = (1 - sp.exp(-2 * g * t)) * (1 + z**2 - x**2 - y**2) / 2
    assert sp.simplify(diff - claimed) == 0, diff
    print(f"   derived: S_A - S_PiA = {sp.factor(diff)}")
    print("   zero at every t > 0 iff x^2 + y^2 = 1 + z^2; with x^2 + y^2 + z^2 <= 1 this gives 1 + 2 z^2 <= 1,")
    print("   so z = 0 and x^2 + y^2 = 1: exactly the pure states on the equator")
    gamma = 0.3
    L = lindbladian_z_dephasing(np.zeros((2, 2), dtype=complex), [gamma])
    T = np.array([vec(S) for S in PAULI]).T / np.sqrt(2)
    LP = T.conj().T @ L @ T
    Pi = build_pi_full(1)
    cases = (("(I + Y)/2, equator", 0.5 * (I2 + Y)),
             ("(I + (X + Y)/sqrt2)/2, equator", 0.5 * (I2 + (X + Y) / np.sqrt(2))),
             ("(I + 0.5 X)/2, inside", 0.5 * (I2 + 0.5 * X)),
             ("(I + Z)/2, pole", 0.5 * (I2 + Z)))
    for name, A in cases:
        a = T.conj().T @ vec(A)
        b = Pi @ a
        diffs = [abs((a.conj() @ expm(LP * t) @ a) - (b.conj() @ expm(LP * t) @ b)) for t in (0.5, 1.7)]
        print(f"   {name}: |S_A - S_PiA| at t = 0.5, 1.7: "
              f"{diffs[0]:.3e}, {diffs[1]:.3e}")


def one_qubit_with_a_field():
    print("7. One qubit with a field of its own: watching along Z at rate gamma, field h along a Pauli letter (exact, sympy)")
    print("   the question of the page: S_A = S_{MA} for the mirror M of the system, which by section 5 is the same as")
    print("   e^(gamma t) S_A(t) even in t, i.e. all odd moments a^T (L + gamma)^(2m+1) a vanish, wherever such an M exists;")
    print("   A = (I + xX + yY + zZ)/2 with x^2 + y^2 + z^2 <= 1, pure on the sphere")
    import sympy as sp
    g, h = sp.symbols("gamma h", positive=True)
    x, y, z = sp.symbols("x y z", real=True)
    X2, Y2, Z2 = sp.symbols("X2 Y2 Z2", nonnegative=True)
    I2s = sp.eye(2); Xs = sp.Matrix([[0, 1], [1, 0]]); Ys = sp.Matrix([[0, -sp.I], [sp.I, 0]]); Zs = sp.diag(1, -1)
    paulis = (I2s, Xs, Ys, Zs)
    kron = sp.kronecker_product
    Ts = sp.Matrix.hstack(*[sp.Matrix(P.reshape(4, 1)) for P in paulis]) / sp.sqrt(2)   # row-major vec, as the framework's
    a = sp.Matrix([1, x, y, z]) / sp.sqrt(2)
    # the mirror Pi of F1 on the letters I, X, Y, Z: I <-> X, Y -> iZ, Z -> iY (column j is the image of letter j)
    Pi1 = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, sp.I], [0, 0, sp.I, 0]])
    assert np.abs(np.array(Pi1, dtype=complex) - build_pi_full(1)).max() < 1e-15, "the same matrix as the framework's Pi"
    # the quarter-turn about Z, S = diag(1, i), on the letters: X -> Y, Y -> -X, I and Z fixed. It commutes with the
    # watching, so it carries the mirror of one system to the mirror of the turned one: Pi turned by a quarter
    # reflects the Y-field generator as -L - 2 gamma. (The one-sided shift by the lit string Y, the colouring of
    # THE_ONE_SQUARE section 8, gives the other form, -L^dagger - 2 gamma; the two differ by the transpose.)
    Ad = sp.Matrix([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    PiY = Ad * Pi1 * Ad.inv()
    results = {}
    for name, Ps in (("X", Xs), ("Z", Zs), ("Y", Ys)):
        Hs = h * Ps
        L = -sp.I * (kron(Hs, I2s) - kron(I2s, Hs.T)) + g * (kron(Zs, Zs.T) - kron(I2s, I2s))
        LP = sp.simplify(Ts.H * L * Ts)
        assert LP == LP.conjugate(), "L is real in the Pauli basis"
        M = LP + g * sp.eye(4)
        odd = [sp.factor(sp.simplify((a.T * M ** (2 * m + 1) * a)[0])) for m in range(3)]
        sq = [sp.expand(o).subs({x ** 2: X2, y ** 2: Y2, z ** 2: Z2}) for o in odd]
        assert all(not (o.free_symbols & {x, y, z}) for o in sq), "the odd moments are polynomials in the squares"
        refl = {"Pi": sp.simplify(Pi1 * LP * Pi1.inv() + LP + 2 * g * sp.eye(4)) == sp.zeros(4),
                "quarter-turned Pi": sp.simplify(PiY * LP * PiY.inv() + LP + 2 * g * sp.eye(4)) == sp.zeros(4)}
        results[name] = (odd, LP, sq, refl)
        print(f"   field along {name}: mirror of the system: " + (", ".join(k for k, v in refl.items() if v) or "neither Pi nor its turned copy"))
        print(f"      first odd moment {odd[0]}")
        print(f"      third {odd[1]}")

    def two_states_exactly(sq, first, second, label):
        """The first two moments force first^2 = 1 and second^2 = z^2; the ball then forces z = second = 0."""
        sol = sp.solve([sq[0], sq[1]], [first, second], dict=True)
        assert sol == [{first: 1, second: Z2}], f"{label}: {sol}"
        assert sp.simplify((X2 + Y2 + Z2 - 1).subs(sol[0])) == 2 * Z2, f"{label}: on that set the ball reads 2 z^2 <= 0"

    # X field: Pi is the mirror; necessity from the first two moments and the ball, sufficiency from Pi-fixedness.
    oX, LPX, sqX, reflX = results["X"]
    assert reflX["Pi"], "Pi reflects the X-field generator"
    two_states_exactly(sqX, X2, Y2, "X field: x^2 = 1 and y^2 = z^2")
    for sgn in (1, -1):
        a_pm = sp.Matrix([1, sgn, 0, 0]) / sp.sqrt(2)
        assert Pi1 * a_pm == sgn * a_pm, "|+> and |-> are fixed by Pi up to sign"
        assert all(o.subs({x: sgn, y: 0, z: 0}) == 0 for o in oX)
    print("   field along X: exactly |+> and |->, among all states pure or mixed (x^2 = 1 and y^2 = z^2 from the first two moments,")
    print("      then the ball x^2 + y^2 + z^2 <= 1 forces y = z = 0; both are fixed by Pi, so every odd moment vanishes)")
    # Y field: the quarter-turned Pi is the mirror; the same two steps with X and Y exchanged; |+i> and |-i> are fixed by it
    oY, LPY, sqY, reflY = results["Y"]
    assert reflY["quarter-turned Pi"] and not reflY["Pi"]
    two_states_exactly(sqY, Y2, X2, "Y field: y^2 = 1 and x^2 = z^2")
    for sgn in (1, -1):
        a_pm = sp.Matrix([1, 0, sgn, 0]) / sp.sqrt(2)
        assert PiY * a_pm == sgn * a_pm, "|+i> and |-i> are fixed by the quarter-turned mirror"
        assert all(o.subs({x: 0, y: sgn, z: 0}) == 0 for o in oY)
    print("   field along Y: the mirror is Pi turned by a quarter about Z; exactly |+i> and |-i>, the same answer turned with the field")
    # Z field: the moments reduce to x^2 + y^2 = 0 and then 1 + z^2 = 0, so no state; and the spectrum pairs about no
    # centre at all, so no map of any form reflects this generator at any price.
    oZ, LPZ, sqZ, reflZ = results["Z"]
    assert not any(reflZ.values())
    x2_from_first = sp.solve(sqZ[0], X2)[0]                                   # x^2 = 1 + z^2 - y^2
    third_reduced = sp.factor(sqZ[1].subs(X2, x2_from_first))
    assert sp.simplify(third_reduced / (Z2 + 1)).free_symbols <= {g, h}, f"the third moment reduces to a multiple of (1 + z^2): {third_reduced}"
    assert sp.solve([sqZ[0], sqZ[1]], [X2, Y2], dict=True) == [], "no state meets the two moments"
    ev = LPZ.eigenvals()                                                       # {0: 2, -2 gamma + 2ih: 1, -2 gamma - 2ih: 1}
    assert ev == {0: 2, -2 * g + 2 * sp.I * h: 1, -2 * g - 2 * sp.I * h: 1}, ev
    for c in (0, g, 2 * g):                                                   # real parts allow a centre only at -gamma, i.e. c = 2 gamma
        paired = all(any(sp.simplify(-mu - c - nu) == 0 for nu in ev) for mu in ev)
        assert not paired, f"the Z-field spectrum pairs about -{c}/2"
    print("   field along Z: the moments force x^2 + y^2 = 0 and then 1 + z^2 = 0: NO state meets itself the same way forwards and")
    print(f"      backwards at the price; the spectrum {dict(ev)} (eigenvalue: multiplicity) pairs about no centre, so no mirror of any form is left.")
    print("      (Pi, no longer a mirror of this system, still maps every pure equator state to an operator with the same")
    print("      self-encounter: both give (1 + e^(-2 gamma t) cos 2ht)/2. That compares a state with an image, not with a mirror image.)")
    # the numeric handshake against the framework generator, with controls that must fail; and a field at a general angle
    # across the watching, where the turn about Z by that angle carries Pi to the mirror and |+> to the two states
    gamma_n, h_n, phi = 0.3, 0.7, 0.7
    T = np.array([vec(S) for S in PAULI]).T / np.sqrt(2)
    turn = np.diag([1, np.exp(1j * phi)])
    Ad_phi = T.conj().T @ np.kron(turn, turn.conj()) @ T                      # rho -> S rho S^dagger on the letters
    fields = (("X", X), ("Z", Z), ("Y", Y), (f"angle {phi} across", np.cos(phi) * X + np.sin(phi) * Y))
    for fname, Pn in fields:
        LPn = T.conj().T @ lindbladian_z_dephasing(h_n * Pn, [gamma_n]) @ T
        if fname in results:
            LPs = np.array(results[fname][1].subs({g: gamma_n, h: h_n}).evalf(), dtype=complex)
            assert np.abs(LPn - LPs).max() < 1e-12, "sympy and framework generators agree"
        else:
            Pi_phi = Ad_phi @ build_pi_full(1) @ np.linalg.inv(Ad_phi)
            assert np.abs(Pi_phi @ LPn @ np.linalg.inv(Pi_phi) + LPn + 2 * gamma_n * np.eye(4)).max() < 1e-12, "Pi turned by phi reflects the turned system"
        ts = (0.4, 1.3)

        def even_gap(A):
            a_n = T.conj().T @ vec(A)
            return max(abs(a_n.conj() @ expm(LPn * t) @ a_n - np.exp(-2 * gamma_n * t) * (a_n.conj() @ expm(-LPn * t) @ a_n)) for t in ts)
        along = np.cos(phi) * X + np.sin(phi) * Y
        gaps = {"|+>": even_gap(0.5 * (I2 + X)), "|->": even_gap(0.5 * (I2 - X)), "|+i>": even_gap(0.5 * (I2 + Y)), "|-i>": even_gap(0.5 * (I2 - Y)),
                "|+phi>": even_gap(0.5 * (I2 + along)), "|-phi>": even_gap(0.5 * (I2 - along)), "|0>": even_gap(0.5 * (I2 + Z)), "(I+0.5X)/2": even_gap(0.5 * (I2 + 0.5 * X))}
        expected = {"X": {"|+>", "|->"}, "Z": set(), "Y": {"|+i>", "|-i>"}, f"angle {phi} across": {"|+phi>", "|-phi>"}}[fname]
        for k, v in gaps.items():
            assert (v < 1e-12) == (k in expected), f"field {fname}, {k}: gap {v:.1e}"
        print(f"   field along {fname}, numeric (gamma {gamma_n}, h {h_n}): " + ", ".join(f"{k} {v:.1e}" for k, v in gaps.items()))


def several_bits_watched():
    print("8. Several bits watched along Z at one rate: forwards = backwards at the price iff the weights W_k on the heights")
    print("   k = number of X or Y letters are symmetric about N/2; mirror-fixed states keep it under every coupling Pi reflects")
    gamma, Delta, phis = 0.3, 0.5, (0.7, 2.5)
    rng = np.random.default_rng(3)
    for N in (2, 3, 4):
        d = 2 ** N
        idxs = list(itertools.product(range(4), repeat=N))
        T = np.array([vec(kron_all([PAULI[i] for i in idx])) for idx in idxs]).T / np.sqrt(d)
        kZ = np.array([sum(1 for i in idx if i in (1, 2)) for idx in idxs])
        L0 = T.conj().T @ lindbladian_z_dephasing(np.zeros((d, d), dtype=complex), [gamma] * N) @ T
        assert np.abs(L0 - np.diag(-2 * gamma * kZ)).max() < 1e-12, "watching alone is the height"
        Pi = build_pi_full(N)
        sigma = N * gamma
        assert np.abs(Pi @ L0 @ np.linalg.inv(Pi) + L0 + 2 * sigma * np.eye(d * d)).max() < 1e-12, "the mirror reflects the heights, k -> N - k"

        def bond(k, P):
            return kron_all([P if l in (k, k + 1) else I2 for l in range(N)])
        H_xxz = sum(bond(k, X) + bond(k, Y) + Delta * bond(k, Z) for k in range(N - 1))
        H_xx = sum(bond(k, X) for k in range(N - 1))
        gens = {}
        for cname, H in (("XXZ", H_xxz), ("XX only", H_xx)):
            L1 = T.conj().T @ lindbladian_z_dephasing(H, [gamma] * N) @ T
            assert np.abs(Pi @ L1 @ np.linalg.inv(Pi) + L1 + 2 * sigma * np.eye(d * d)).max() < 1e-12, f"{cname} is in the mirror's class"
            gens[cname] = L1
        ghz = np.zeros(d); ghz[0] = ghz[-1] = 1
        w = np.zeros(d)
        for k in range(N):
            w[1 << k] = 1
        plus = kron_all([np.array([1, 1])] * N)
        plus_minus = kron_all([np.array([1, (-1) ** k]) for k in range(N)])
        random_angles = kron_all([np.array([1, np.exp(1j * th)]) for th in rng.uniform(0, 2 * np.pi, N)])
        tilted = kron_all([np.array([1, 1])] * (N - 1) + [np.array([np.cos(0.4), np.sin(0.4)])])
        cases = [("cat (GHZ)", ghz), ("W", w), ("|+>^N", plus), ("|+>|->|+>...", plus_minus)]
        cases += [(f"one equator angle {phi} on every bit", kron_all([np.array([1, np.exp(1j * phi)])] * N)) for phi in phis]
        cases += [("random equator angles (one draw)", random_angles), ("one bit off the equator", tilted), ("random pure", rng.normal(size=d) + 1j * rng.normal(size=d))]
        mirror_fixed = {"cat (GHZ)", "|+>^N", "|+>|->|+>..."}
        print(f"   N = {N}:")
        for name, psi in cases:
            if name == "W" and N == 2:
                continue                                                   # W at N = 2 is the Bell state Psi+, mirror-fixed; not a W-type control
            rho = np.outer(psi, psi.conj()) / np.vdot(psi, psi)
            a = T.conj().T @ vec(rho)
            W = np.array([np.sum(np.abs(a[kZ == k]) ** 2) for k in range(N + 1)])
            sym = np.abs(W - W[::-1]).max()
            fixed = min(np.abs(Pi @ a - a).max(), np.abs(Pi @ a + a).max())

            def even_gap(L):
                return max(abs(a.conj() @ expm(L * t) @ a - np.exp(-2 * sigma * t) * (a.conj() @ expm(-L * t) @ a)) for t in (0.4, 1.3))
            g0, g1, g2 = even_gap(L0), even_gap(gens["XXZ"]), even_gap(gens["XX only"])
            assert (g0 < 1e-12) == (sym < 1e-12), f"{name}: watching alone, symmetric weights and the even self-encounter must agree"
            assert (fixed < 1e-12) == (name in mirror_fixed), f"{name}: mirror-fixed {fixed:.1e}"
            if name in mirror_fixed:
                assert g1 < 1e-12 and g2 < 1e-12, f"{name}: mirror-fixed, so it must hold under every coupling of the class"
            if name == "cat (GHZ)":
                assert abs(W[0] - 0.5) < 1e-12 and abs(W[N] - 0.5) < 1e-12 and np.abs(W[1:N]).max() < 1e-12, "the cat: half on each edge"
            if name.startswith("random equator angles"):
                assert sym < 1e-12 and g1 > 1e-3 and g2 > 1e-3, f"{name}: symmetric alone, this draw not under either coupling ({g1:.1e}, {g2:.1e})"
            if name.startswith("one equator angle"):
                assert sym < 1e-12 and g1 < 1e-12 and g2 > 1e-3, f"{name}: kept under XXZ (U(1), the turned Pi reflects it), lost under XX only ({g1:.1e}, {g2:.1e})"
            if name in ("W", "one bit off the equator"):
                assert sym > 1e-2 and g0 > 1e-2, f"{name} is the control that must fail"
            print(f"      {name:36s} W = {np.round(W, 3)}  |W - reversed| {sym:.1e}  mirror-fixed {fixed:.1e}  gap: alone {g0:.1e}, XXZ {g1:.1e}, XX only {g2:.1e}")


def coupled_bits_with_a_field():
    print("9. Coupled bits with a field, the conjugate pair, and the string eigenstates: the mirror image of an X^N eigenstate")
    print("   is its complex conjugate up to sign, which the self-encounter cannot see; under a field along X the mirror stays")
    gamma, Delta, h, phi = 0.3, 0.5, 0.7, 0.7
    rng = np.random.default_rng(5)
    ts = (0.4, 1.3)
    for N in (2, 3):
        d = 2 ** N
        sigma = N * gamma
        idxs = list(itertools.product(range(4), repeat=N))
        T = np.array([vec(kron_all([PAULI[i] for i in idx])) for idx in idxs]).T / np.sqrt(d)
        Pi = build_pi_full(N)
        Id = np.eye(d * d)

        def bond(k, P):
            return kron_all([P if l in (k, k + 1) else I2 for l in range(N)])

        def field(P):
            return sum(kron_all([P if l == k else I2 for l in range(N)]) for k in range(N))
        H_xxz = sum(bond(k, X) + bond(k, Y) + Delta * bond(k, Z) for k in range(N - 1))
        H_xx = sum(bond(k, X) for k in range(N - 1))
        XN, YN = kron_all([X] * N), kron_all([Y] * N)
        n_phi = np.cos(phi) * X + np.sin(phi) * Y
        turn = np.diag([1, np.exp(1j * phi)])
        Ad_phi = T.conj().T @ np.kron(kron_all([turn] * N), kron_all([turn] * N).conj()) @ T   # the global Z-turn by phi on the letters
        Pi_phi = Ad_phi @ Pi @ np.linalg.inv(Ad_phi)

        def gen(H, gammas=None):
            return T.conj().T @ lindbladian_z_dephasing(H, gammas or [gamma] * N) @ T

        def reflects(M, L):
            return np.abs(M @ L @ np.linalg.inv(M) + L + 2 * sigma * Id).max()

        def pairing(L):
            ev = np.linalg.eigvals(L)
            return max(min(abs(m - e) for e in ev) for m in -ev - 2 * sigma)

        def coords(rho):
            return T.conj().T @ vec(rho)

        def even_gap(L, a):
            return max(abs(a.conj() @ expm(L * t) @ a - np.exp(-2 * sigma * t) * (a.conj() @ expm(-L * t) @ a)) for t in ts)

        def encounter(L, a):
            return np.array([a.conj() @ expm(L * t) @ a for t in ts])

        def eigenstate(P, val):
            w, V = np.linalg.eigh(P)
            cols = V[:, np.abs(w - val) < 1e-9]
            psi = cols @ (rng.normal(size=cols.shape[1]) + 1j * rng.normal(size=cols.shape[1]))
            psi /= np.linalg.norm(psi)
            return np.outer(psi, psi.conj())

        print(f"   N = {N}:")
        # (a) a field along X, the same on every bit: Pi stays the mirror exactly when the coupling kept it
        for cname, H in (("XXZ", H_xxz), ("XX only", H_xx)):
            L_x = gen(H + h * field(X))
            assert reflects(Pi, L_x) < 1e-12, f"{cname} + X field: Pi still reflects (the field is real and commutes with X^N)"
            L_0 = gen(H)
            for val in (+1, -1):
                rho = eigenstate(XN, val)
                a, a_conj = coords(rho), coords(rho.conj())
                image = Pi @ a
                fixed = min(np.abs(image - a).max(), np.abs(image + a).max())
                conj_up_to_sign = min(np.abs(image - a_conj).max(), np.abs(image + a_conj).max())
                assert fixed > 1e-2, f"a random X^N = {val:+d} eigenstate is NOT mirror-fixed ({fixed:.2f})"
                assert conj_up_to_sign < 1e-12, "its mirror image is its complex conjugate up to sign (rho X^N = +-rho, then the transpose)"
                assert np.abs(encounter(L_x, a) - encounter(L_x, a_conj)).max() < 1e-12, "the conjugate meets itself alike (L = L^T, real H)"
                gx, g0 = even_gap(L_x, a), even_gap(L_0, a)
                assert gx < 1e-12 and g0 < 1e-12, f"X^N = {val:+d} eigenstate keeps the property with and without the X field ({gx:.1e}, {g0:.1e})"
                print(f"      {cname:8s} X field, X^N = {val:+d} eigenstate: mirror-fixed {fixed:.2f}, image = conjugate up to sign {conj_up_to_sign:.1e}, gap {gx:.1e} (no field {g0:.1e})")
            rho_y = eigenstate(YN, +1)
            gy = even_gap(L_x, coords(rho_y))
            assert gy > 1e-3, f"a Y^N eigenstate under the X field is the control that fails ({gy:.1e})"
            mixed = np.eye(d) / d
            gm = even_gap(L_x, coords(mixed))
            assert gm > 1e-2, f"the completely mixed state is unchanged by X^N under conjugation and fails ({gm:.1e}): eigenstate, not invariance"
            print(f"      {cname:8s} X field controls: Y^N eigenstate gap {gy:.1e}, completely mixed gap {gm:.1e}")
            # Y^N eigenstates without a field keep it too (checked, both signs; the reason is not the one above: L is still L^T,
            # but rho Y^N = +-rho makes the image rho^T X^N a shifted copy, not the conjugate)
            for val in (+1, -1):
                gyy = even_gap(L_0, coords(eigenstate(YN, val)))
                assert gyy < 1e-12, f"Y^N = {val:+d} eigenstate without field under {cname} ({gyy:.1e})"
        # (b) the field turned to the angle phi in the equator plane: the turned mirror reflects XXZ (U(1)) and not XX only;
        # under XXZ the eigenstates of the turned string keep it, the cat (an X^N eigenstate, not one of the turned string) does not
        L_phi_xxz, L_phi_xx = gen(H_xxz + h * field(n_phi)), gen(H_xx + h * field(n_phi))
        assert reflects(Pi_phi, L_phi_xxz) < 1e-12 and pairing(L_phi_xxz) < 1e-10, "the turned Pi reflects XXZ with the turned field"
        assert reflects(Pi_phi, L_phi_xx) > 1e-2 and pairing(L_phi_xx) > 1e-2, "with XX only, the turned field leaves no mirror and no pairing"
        rho_turned = eigenstate(kron_all([n_phi] * N), +1)
        ghz = np.zeros(d); ghz[0] = ghz[-1] = 1 / np.sqrt(2)
        g_turned, g_cat = even_gap(L_phi_xxz, coords(rho_turned)), even_gap(L_phi_xxz, coords(np.outer(ghz, ghz)))
        assert g_turned < 1e-12 and g_cat > 1e-3, f"turned-string eigenstate keeps it ({g_turned:.1e}), the cat does not ({g_cat:.1e})"
        print(f"      field at angle {phi} in the equator plane, XXZ: turned-string eigenstate gap {g_turned:.1e}, cat gap {g_cat:.1e}; XX only: no mirror ({pairing(L_phi_xx):.2f})")
        # (c) the conjugate pair parts under a Y field for a random pure state of coupled bits (H not real, L != L^T,
        #     necessary and not sufficient: one bit alone keeps the pair under a Y field)
        L_y = gen(H_xxz + h * field(Y))
        psi = rng.normal(size=d) + 1j * rng.normal(size=d); psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        part_y = np.abs(encounter(L_y, coords(rho)) - encounter(L_y, coords(rho.conj()))).max()
        part_x = np.abs(encounter(gen(H_xxz + h * field(X)), coords(rho)) - encounter(gen(H_xxz + h * field(X)), coords(rho.conj()))).max()
        assert part_x < 1e-12 and part_y > 1e-3, f"the conjugate pair: alike under the X field ({part_x:.1e}), parted under the Y field ({part_y:.1e})"
        print(f"      a random pure state and its conjugate: |S - S*| under X field {part_x:.1e}, under Y field {part_y:.1e}")
        # (c') one bit alone under the Y field keeps the pair (its Y coefficient only decays), and every two-site Pauli
        #      coupling that Pi reflects is real, so the conjugate sentence applies to every coupling that keeps the mirror
        if N == 2:
            T1 = np.array([vec(P) for P in PAULI]).T / np.sqrt(2)
            L1y = T1.conj().T @ lindbladian_z_dephasing(h * Y, [gamma]) @ T1
            b = rng.normal(size=3); b *= 0.8 / np.linalg.norm(b)
            rho1 = 0.5 * (I2 + b[0] * X + b[1] * Y + b[2] * Z)
            a1, a1c = T1.conj().T @ vec(rho1), T1.conj().T @ vec(rho1.conj())
            one_bit = max(abs(a1.conj() @ expm(L1y * t) @ a1 - a1c.conj() @ expm(L1y * t) @ a1c) for t in ts)
            assert np.abs(L1y - L1y.T).max() > 1e-2 and one_bit < 1e-12, f"one bit, Y field: L != L^T yet the conjugate pair holds ({one_bit:.1e})"
            reflected, non_real = [], []
            for idx in idxs:
                P = kron_all([PAULI[i] for i in idx])
                if np.abs(P - np.eye(d)).max() < 1e-12:
                    continue
                if reflects(Pi, gen(P)) < 1e-12:
                    reflected.append(idx)
                    if np.abs(P.imag).max() > 1e-12:
                        non_real.append(idx)
            assert reflected and not non_real, f"every two-site Pauli coupling Pi reflects is real ({len(reflected)} reflected, {len(non_real)} non-real)"
            print(f"      one bit under the Y field: conjugate pair gap {one_bit:.1e} (L != L^T: {np.abs(L1y - L1y.T).max():.2f}); two-site Pauli couplings Pi reflects: {len(reflected)}, all real")
        # (d) a field along the watching, the same on both bits of a pair: no mirror, yet the m = 0 pair keeps it; unequal fields break it
        if N == 2:
            L_z = gen(H_xxz + h * field(Z))
            assert reflects(Pi, L_z) > 1e-2 and pairing(L_z) > 0.1, f"Z field: Pi broken, spectrum pairs about no centre ({pairing(L_z):.2f})"
            singlet = np.array([0, 1, -1, 0]) / np.sqrt(2); triplet = np.array([0, 1, 1, 0]) / np.sqrt(2)
            bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
            L_z_unequal = gen(H_xxz + h * kron_all([Z, I2]) + 0.2 * kron_all([I2, Z]))
            for name, psi in (("singlet", singlet), ("triplet", triplet)):
                a = coords(np.outer(psi, psi.conj()))
                g, gu = even_gap(L_z, a), even_gap(L_z_unequal, a)
                assert g < 1e-12 and gu > 1e-2, f"{name}: kept under the equal Z field ({g:.1e}), lost under unequal fields ({gu:.1e})"
                print(f"      Z field on both bits, XXZ: {name} gap {g:.1e}; unequal fields {gu:.1e}; mirror pairing {pairing(L_z):.2f}")
            g_bell = even_gap(L_z, coords(np.outer(bell, bell)))
            assert g_bell > 1e-2, f"Bell+ under the Z field is turned by it and fails ({g_bell:.1e})"
            print(f"      Z field on both bits, XXZ: Bell+ (turned by the field) gap {g_bell:.1e}")


if __name__ == "__main__":
    self_encounter_blind_to_phase()
    bell_pair_is_all_between()
    count_selves_and_between()
    reciprocity()
    mirror_self_encounter()
    mirror_twin_exists()
    one_qubit_with_a_field()
    several_bits_watched()
    coupled_bits_with_a_field()
