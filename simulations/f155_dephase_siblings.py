"""F155 siblings gate: the X- and Y-dephase palindrome conventions.

THE QUESTION this file answers, and it is the OpenArc polarity_break_dephase_siblings.
F155 (docs/proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md) closed the
polarity break of the no-jump generator G_H : rho -> -i(H rho - rho H^dag), with
H = A + iB and both parts Hermitian, for the Z-dephasing palindrome convention:

    asymmetry_Z = 4^(N+1) * sum_{s : bit_b(s) odd} (-1)^(#Z(s)) * a_s * b_s .

Its scope note, now section (h), left the X- and Y-dephase versions unclaimed and
named the route: recompute the per-site swap rule for pi_X and pi_Y and reassemble.
That is done here, and the ANSWER IS NOT TWO NEW LAWS. Both siblings reduce to F155:

  * Y IS NOT A SIBLING AT ALL. Pi_Y = Pi_Z^-1 exactly (owned by the repo:
    docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md and, in the same directory,
    PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md), so Ad_{Pi_Y} = Ad_{Pi_Z}^3,
    and inverting an order-4 map EXCHANGES its +i and -i eigenprojections. Hence
    for EVERY superoperator M, not only for generators of the above form and with
    no bilinear argument at all, asymmetry_Y(M) = -asymmetry_Z(M). Gate S4 checks
    this bit-exactly on arbitrary complex M. The Y law is therefore F155 with a
    global minus, and it is the SAME operator that F155's section (a) already
    prices as the twisted-versus-untwisted pairing freedom: that freedom is
    conjugation by D = diag((-1)^#Y), and the repo owns D . Pi_Z . D = Pi_Y
    (PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md), so the Y convention in this
    pairing IS the Z convention in the untwisted one.

  * X IS THE ONE GENUINE SIBLING, and it is F155 in a rotated frame:

        asymmetry_X = 4^(N+1) * sum_{s : bit_a(s) odd} (-1)^(#X(s)) * a_s * b_s ,

    with bit_a(s) = #X(s) + #Y(s) mod 2. Two independent routes reach it and both
    are gated: the DIRECT one (recompute the swap table for pi_X, S2/S3, then
    reassemble, S6) and the TRANSPORT one (S10: asymmetry_X(A, B) equals
    asymmetry_Z(hAh, hBh) for h = H^(tensor N), the SAME Hadamard on every site).
    The transport is a route and not a coincidence because the repo owns the
    intertwiner Q_zx . Pi_Z . Q_zx^-1 = Pi_X
    (PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md section (d)), which is the same
    Hadamard route that proof uses for the F112 siblings. A Hadamard on a PROPER
    SUBSET of sites does NOT do this: it produces a mixed-letter law, not the X
    one, and S10 gates that too.

ONE LAW WITH THE LETTER AS A PARAMETER, in the pairing named below. Reading the
three swap tables together (S2) gives a single statement rather than three:

    asymmetry_l = 4^(N+1) * sum_{s : chi_l(s) odd} (-1)^(#l(s)) * a_s * b_s ,

where chi_l is the parity that Pi_l^2 grades by (F38/F88a: bit_a for l = X,
bit_b for l = Y and Z, which framework.symmetry.pi_squared_eigenvalue already
carries) and #l counts the occurrences of the dephasing letter ITSELF.

WHICH HALF OF THAT IS CONVENTION-FREE, because "the letter itself" is pretty and
only half true. The SUPPORT chi_l is convention-free. The WEIGHT is not: it is
the dephasing letter only in the twisted (pipeline) pairing this file computes
in. In the untwisted pairing the two per-letter signs exchange and the weights
become (-1)^#Y for BOTH Z and X and (-1)^#Z for Y, so there the weight letter is
the dephasing letter for none of them (S3 measures the exchange). This is the
same freedom F155 section (a) prices for its own law, inherited unchanged.

THE CONSEQUENCE WORTH QUOTING, and it is a SILENCE and not a balance. The
configuration that makes F155 subsume F113, a per-site Z-drive meeting amplitude
damping, reads exactly 0.0 under Pi_X (S8), and the two +-i norms are
INDIVIDUALLY and exactly zero: there is no cancellation, there is nothing there.
The repo has a name for this third answer beside balanced and broken (the SILENT
case; framework's is_structurally_degenerate). It is structural at every N and
every rate, and it is more general than the F113 family: ANY A, B whose Pauli
content is entirely bit_a-EVEN lies in the +1 eigenspace of Ad_{Pi_X}^2, so
M_anti is the exact zero array and both halves vanish. S8 gates the general
statement, not only the F113 instance. The asymmetry is therefore a property of
the physics PAIRED WITH a palindrome convention, which is the sharpest form of
the scope warning F155's gate G13 raised: G13 showed the value MOVES under a
per-site Hadamard, and S10 here shows WHERE it moves to.

WHY S2 CARRIES A HARD-CODED TABLE, and it is the lesson of this file's own audit.
A first version DERIVED the closed form's support and weight from the measured
swap table, and called that a safeguard against a hard-coded law masking a wrong
table. It is the reverse. If Pi itself is wrong, the measured table is wrong AND
the derived law is wrong in the MATCHING way, so the comparison succeeds on two
consistently wrong quantities. Demonstrated, not feared: flipping the phase sign
of the X convention inside framework.symmetry.pi_action left S2 and S6 both
PASSING, and only the independent routes S9 and S10 fired. So the three expected
tables are written out as literals below, S2 gates the MEASURED table against
them by exact set equality, and every closed form is built from the literals. The
measured table remains, as the cross-check it always should have been.

TOLERANCES, per the repo's no-rounding rule, with the error model stated here
rather than by reference. The asymmetry is a difference of two sums of squared
moduli, each O(||M||_F^2), so its rounding floor is C * eps * ||M||_F^2 with C an
O(1) constant depending on neither N nor the overall scale. That is a LAW and S6
gates it as one: it sweeps the input scale over nine decades at three values of
N and requires the ratio to STAY PUT (its spread bounded), not merely to come in
under a ceiling. S1, S2, S3, S4, S5, S7 and half of S8 and S9 compare EXACTLY
(==, np.array_equal, or a residual asserted to be 0.0), because an exact route
exists in each. S7's off-diagonal zeros are exact and are asserted as such rather
than gated, which is the repo's case (1) and not case (2). Every float gate that
asserts a MATCH also carries a signal guard, so that it can never be satisfied by
both sides being zero.

STAGE 0, what the repo already held (swept 2026-08-19 by agents over
docs/ANALYTICAL_FORMULAS.md, docs/proofs/, experiments/ including null results,
hypotheses/ and reflections/, docs/GLOSSARY.md, docs/CAUGHT_ERRORS.md,
compute/ including the OpenArcs registry and both Confirmations registries, and
simulations/framework/):
  * OWNED AND PROVEN, and it is what makes Y collapse: Pi_Y = Pi_Z^-1 at every N
    (PROOF_PI_FACTORS_AS_R_TIMES_D.md,
    PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md,
    and reflections/D_PI_Z_EQUALS_PI_Y.md), together with D . Pi_Z . D = Pi_Y and
    the Hadamard intertwiner Q_zx . Pi_Z . Q_zx^-1 = Pi_X. What is added here is
    that the first fact ALONE decides the Y asymmetry, universally in M and
    without touching the generator.
  * OWNED AND PROVEN: the Pi^2 grading per letter, F38/F88a, bit_a for X and
    bit_b for Y and Z (docs/ANALYTICAL_FORMULAS.md, which states the letter's
    Klein pair Z (0,1), X (1,0), Y (1,1) alongside it; typed at
    compute/RCPsiSquared.Core/Symmetry/PiOperator.cs; framework.symmetry.
    pi_squared_eigenvalue). S5 re-checks an owned fact, deliberately, because a
    typed claim in compute/ had contradicted it; see the note under S5.
  * OWNED, LETTER-PARAMETERISED AND ALREADY GATED, and an earlier draft of this
    file wrongly said otherwise: simulations/f112_klein_v4_cross_dephase_verify.py
    passes "X" and "Y" into polarity_coordinates_from_L (its
    asymmetry_for_dephase, called with those letters at lines 210, 240 and 284 and
    in the letter loops below them) and already implements the chi_l half of the
    law above, choosing bit_a for X and bit_b otherwise. That script is the
    verifier of PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md, which this file cites
    twice, so the earlier claim "no caller had ever passed a letter to the
    polarity primitive" was false with its counterexample in a file linked from
    the same page. What had never been asked in a non-Z convention is the VALUE of
    the asymmetry for THIS generator; the F112 work asks whether it VANISHES, and
    for the full Lindbladian.
  * A REAL GAP, the narrow half of that bullet: framework.diagnostics.
    polarity_coordinates hard-wires build_pi_full(N), i.e. Z, in its chain-facing
    entry points; only polarity_coordinates_from_L takes an injected Pi. Every
    gate below injects it explicitly.
  * NOTHING FOUND: fw.Confirmations and the C# ConfirmationsRegistry (no hardware
    confirmation touches this quantity), docs/GLOSSARY.md (no entry for the
    polarity asymmetry at all, in any convention), experiments/ (the X-dephase
    question appears as an OPEN question in PI_AS_TIME_REVERSAL.md and nowhere as
    an answer).
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("PYTHONIOENCODING", "utf-8")
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np

from framework.pauli import (
    _k_to_indices,
    _pauli_label,
    _vec_to_pauli_basis_transform,
    bit_a,
    bit_b,
    pauli_string,
)
from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L

EPS = float(np.finfo(float).eps)
EPS_RATIO_BOUND = 8.0
EPS_RATIO_SPREAD = 20.0    # how far C may WANDER across nine decades and three N
SIGNAL_GUARD = 1.0e6       # a matching gate must compare something non-trivial
LETTERS = ("Z", "X", "Y")
GATES = []
_CACHE = {}

# The law, written out. S2 gates the MEASURED swap table against these, and every
# closed form is built from these. See "WHY S2 CARRIES A HARD-CODED TABLE".
EXPECTED_SUPPORT = {"Z": frozenset("YZ"), "X": frozenset("XY"), "Y": frozenset("YZ")}
EXPECTED_WEIGHT = {"Z": frozenset("Z"), "X": frozenset("X"), "Y": frozenset("Y")}

I2 = np.eye(2, dtype=complex)
PAULI = {
    "I": I2,
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)
SIGMA_PLUS = np.array([[0, 0], [1, 0]], dtype=complex)
HADAMARD = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0)


def gate(name, passed, detail):
    GATES.append((name, bool(passed), detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}")
    return bool(passed)


def transform(N):
    """The Pauli transform in the order='F' convention, named at the call site.

    Not left to a default: flipping it to 'C' inverts the entire swap table (it
    conjugates the lifts by D = diag((-1)^#Y)), so the convention is written out.
    """
    if N not in _CACHE:
        _CACHE[N] = _vec_to_pauli_basis_transform(N, order="F")
    return _CACHE[N]


def to_pauli(L_vec, N):
    """Row-stack superoperator read against the order='F' Pauli transform.

    This is F155's twisted pairing, the one F113's registered numbers live in.
    The sibling laws inherit the same freedom, and S3 measures both spellings.
    """
    T = transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def ad_pi(M, N, letter):
    Pi = build_pi_full(N, letter)
    return Pi @ M @ Pi.conj().T


def lift(p_matrix, side, N):
    """The one-sided multiplication superoperator, in F155's twisted pairing."""
    Id = np.eye(2 ** N, dtype=complex)
    v = np.kron(p_matrix, Id) if side == "L" else np.kron(Id, p_matrix.T)
    return to_pauli(v, N)


def free_lift(p_matrix, side, N):
    """The same object in the untwisted (basis-free) pairing, for S3's mirror."""
    strings = [pauli_string(_k_to_indices(k, N)) for k in range(4 ** N)]
    d = 2 ** N
    if side == "L":
        def act(X):
            return p_matrix @ X
    else:
        def act(X):
            return X @ p_matrix
    return np.array([[np.trace(strings[a].conj().T @ act(strings[b])) / d
                      for b in range(4 ** N)] for a in range(4 ** N)])


def site_op(N, l, op2):
    ops = [I2] * N
    ops[l] = op2
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def generator(A, B, N):
    """-i(H rho - rho H^dag) in vec form, H = A + iB with A, B Hermitian."""
    Id = np.eye(2 ** N, dtype=complex)
    return (-1j * (np.kron(A, Id) - np.kron(Id, A.T))
            + (np.kron(B, Id) + np.kron(Id, B.T)))


IMAG_RESIDUE = {"worst_ratio": 0.0}


def pauli_coeffs(M, N):
    """Real Pauli coefficients of a Hermitian M, with the discard accounted for.

    Hermiticity is asserted EXACTLY, because an exact route exists: M + M^dag is
    bit-exactly self-adjoint. The imaginary part of each coefficient is then zero
    in exact arithmetic, but the trace sums 2^N products, so what is computed is a
    rounding RESIDUE, not a zero. It is therefore not asserted away and not
    silently dropped: it is measured against the floor eps*||M||_F it can be
    bounded by, and the worst ratio seen anywhere in the run is reported by S6.
    Silently taking the real part is what would turn a non-Hermitian input into a
    plausible number, and the Hermiticity assert above is what forbids that.
    """
    assert np.array_equal(M, M.conj().T), "pauli_coeffs: input is not Hermitian"
    d = 2 ** N
    floor = EPS * float(np.sqrt(np.sum(np.abs(M) ** 2)))
    out = []
    for k in range(4 ** N):
        c = np.trace(pauli_string(_k_to_indices(k, N)) @ M) / d
        if floor > 0.0:
            IMAG_RESIDUE["worst_ratio"] = max(IMAG_RESIDUE["worst_ratio"],
                                              abs(float(np.imag(c))) / floor)
        out.append(float(np.real(c)))
    return out


def polarity(M, N, letter):
    """The full polarity reading, so a gate can see the two halves separately."""
    return polarity_coordinates_from_L(M, N, 0.0, Pi=build_pi_full(N, letter))


def asym(M, N, letter):
    return float(polarity(M, N, letter)["asymmetry"])


def eps_ratio(residual, scale):
    return abs(residual) / (EPS * scale) if scale > 0.0 else 0.0


def hadamard_full(N):
    """h = H^(tensor N): the SAME Hadamard on EVERY site. A partial one is not this."""
    out = HADAMARD
    for _ in range(N - 1):
        out = np.kron(out, HADAMARD)
    return out


def hadamard_first_site_only(N):
    """The deliberate wrong version, so S10 can show the identity needs all sites."""
    out = HADAMARD
    for _ in range(N - 1):
        out = np.kron(out, I2)
    return out


def random_hermitian(rng, N):
    d = 2 ** N
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    return M + M.conj().T


def swap_table(letter, pairing="pipeline"):
    """(eps_L, eps_R) per single-qubit letter, or None where no exact sign exists.

    Ad_pi(L_p) = eps_L(p) * R_p and Ad_pi(R_p) = eps_R(p) * L_p, compared with
    np.array_equal: these are signed permutation matrices, so an exact route
    exists and a tolerance here would be an error code. This is the MEASUREMENT;
    EXPECTED_SUPPORT and EXPECTED_WEIGHT are the claim it is checked against.
    """
    make = lift if pairing == "pipeline" else free_lift
    table = {}
    for p in "IXYZ":
        L = make(PAULI[p], "L", 1)
        R = make(PAULI[p], "R", 1)
        e_l = next((s for s in (+1, -1)
                    if np.array_equal(ad_pi(L, 1, letter), s * R)), None)
        e_r = next((s for s in (+1, -1)
                    if np.array_equal(ad_pi(R, 1, letter), s * L)), None)
        table[p] = (e_l, e_r)
    return table


def support_letters(table):
    """The letters where the two swap signs DISAGREE: the law's support."""
    return frozenset(p for p in "IXYZ"
                     if table[p][0] is not None and table[p][1] is not None
                     and table[p][0] != table[p][1])


def weight_letters(table):
    """The letters carrying eps_R = -1: the exponent of the closed form's sign."""
    return frozenset(p for p in "IXYZ" if table[p][1] == -1)


def closed_form(A, B, N, letter):
    """4^(N+1) * sum over support-odd strings of (-1)^#(weight letters) a_s b_s.

    Built from the HARD-CODED tables, not from the measured one. A closed form
    derived from the measurement cannot disagree with it, and that blindness is
    exactly what this file's audit found and what S2 now closes.
    """
    supp, neg = EXPECTED_SUPPORT[letter], EXPECTED_WEIGHT[letter]
    ca, cb = pauli_coeffs(A, N), pauli_coeffs(B, N)
    total = 0.0
    for k in range(1, 4 ** N):
        label = _pauli_label(k, N)
        if sum(1 for c in label if c in supp) % 2 == 0:
            continue
        total += ((-1) ** sum(1 for c in label if c in neg)) * ca[k] * cb[k]
    return 4 ** (N + 1) * total


# ----------------------------------------------------------------------
print("=" * 78)
print("F155 dephase siblings: the X and Y palindrome conventions")
print("=" * 78)

# ---------------------------------------------------------------- S1
# Pi_Y = Pi_Z^-1, exactly. This single fact decides the whole Y case (S4).
ok, detail = True, []
for N in (1, 2, 3, 4):
    Pz, Py = build_pi_full(N, "Z"), build_pi_full(N, "Y")
    unitary = np.array_equal(Pz @ Pz.conj().T, np.eye(4 ** N, dtype=complex))
    order4 = np.array_equal(np.linalg.matrix_power(Pz, 4),
                            np.eye(4 ** N, dtype=complex))
    d_inv = float(np.max(np.abs(Py - Pz.conj().T)))
    d_sym = float(np.max(np.abs(Pz - Pz.T)))
    ok &= unitary and order4 and (d_inv == 0.0) and (d_sym == 0.0)
    detail.append(f"N={N}: unitary {unitary}, Pi^4 = I {order4}, "
                  f"|Pi_Y - Pi_Z^dag| = {d_inv!r}, |Pi_Z - Pi_Z^T| = {d_sym!r}")
gate("S1", ok,
     "Pi_Z is EXACTLY unitary (Pi Pi^dag is the identity entry for entry, so its "
     "adjoint IS its inverse) and Pi^4 is EXACTLY the identity (so Ad_Pi has order "
     "dividing 4 and Lemma 1's eigenprojection formula applies), and Pi_Y equals "
     "that adjoint entry for entry, and Pi_Z is exactly symmetric so the adjoint "
     "is also the entry-wise conjugate: Pi_Y, Pi_Z^-1 and conj(Pi_Z) are ONE "
     "matrix. " + "; ".join(detail)
     + ". The identity is owned by the repo (PROOF_PI_FACTORS_AS_R_TIMES_D.md, "
       "PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md); it is re-checked here "
       "because S4 rests entirely on it")

# ---------------------------------------------------------------- S2
# The swap tables per convention, MEASURED and checked against the written law.
measured = {letter: swap_table(letter) for letter in LETTERS}
ok, rows = True, []
for letter in LETTERS:
    t = measured[letter]
    well_formed = all(v[0] is not None and v[1] is not None for v in t.values())
    supp_ok = well_formed and support_letters(t) == EXPECTED_SUPPORT[letter]
    weight_ok = well_formed and weight_letters(t) == EXPECTED_WEIGHT[letter]
    ok &= well_formed and supp_ok and weight_ok
    rows.append(f"Pi_{letter}: " + " ".join(
        f"{p}({t[p][0]},{t[p][1]})" for p in "IXYZ")
        + f" -> support {''.join(sorted(support_letters(t)))} vs expected "
          f"{''.join(sorted(EXPECTED_SUPPORT[letter]))} [{supp_ok}], weight "
          f"(-1)^#{''.join(sorted(weight_letters(t)))} vs expected "
          f"#{''.join(sorted(EXPECTED_WEIGHT[letter]))} [{weight_ok}]")
gate("S2", ok,
     "the per-site swap rule Ad(L_p)=eps_L*R_p, Ad(R_p)=eps_R*L_p is bit-exact "
     "for all four letters in all three conventions, AND the support and weight "
     "it implies are checked by EXACT SET EQUALITY against the law as written, "
     "which is what makes this gate able to fail at all: " + " | ".join(rows)
     + ". The support is where the two signs DISAGREE (only their difference "
       "survives the assembly) and the weight is the eps_R=-1 letter. So the "
       "support parity is bit_b for Z and Y and bit_a for X, exactly the parity "
       "Pi^2 grades by (F38/F88a, gate S5), and the weight letter is the "
       "dephasing letter itself IN THIS PAIRING ONLY, see S3")

# ---------------------------------------------------------------- S3
# The string form, and both pairings, per convention.
ok, detail, free_rows = True, [], []
for letter in LETTERS:
    t_pipe, t_free = measured[letter], swap_table(letter, "free")
    if any(v[0] is None or v[1] is None for v in t_pipe.values()):
        # S2 has already failed; report that rather than crashing on a None sign
        ok = False
        detail.append(f"Pi_{letter}: no exact per-site sign exists, see S2")
        free_rows.append(f"Pi_{letter}: not evaluated")
        continue
    for N in (2, 3):
        worst = 0.0
        for k in range(4 ** N):
            label = _pauli_label(k, N)
            sig = pauli_string(_k_to_indices(k, N))
            s_l = int(np.prod([t_pipe[c][0] for c in label]))
            s_r = int(np.prod([t_pipe[c][1] for c in label]))
            L, R = lift(sig, "L", N), lift(sig, "R", N)
            worst = max(worst,
                        float(np.max(np.abs(ad_pi(L, N, letter) - s_l * R))),
                        float(np.max(np.abs(ad_pi(R, N, letter) - s_r * L))))
        ok &= (worst == 0.0)
        detail.append(f"Pi_{letter} N={N} dev {worst!r}")
    swapped = all(t_free[p] == (t_pipe[p][1], t_pipe[p][0]) for p in "IXYZ")
    ok &= swapped
    free_rows.append(f"Pi_{letter} untwisted weight (-1)^#"
                     f"{''.join(sorted(weight_letters(t_free)))}")
gate("S3", ok,
     "the string form Ad(L_s) = (-1)^#(eps_L letters) R_s and its partner hold "
     "entry-wise EXACTLY at N=2 AND N=3 in every convention, each deviation "
     "reported per N: " + ", ".join(detail)
     + ". And in the untwisted pairing the two per-letter signs EXCHANGE in every "
       "convention, which is where the pretty half of the law loses its frame: "
     + ", ".join(free_rows)
     + ". So the SUPPORT is convention-free but the WEIGHT is not, and 'the "
       "dephasing letter itself' is a statement about the twisted pairing this "
       "file computes in, inherited unchanged from F155 section (a)")

# ---------------------------------------------------------------- S4
# The Y case, decided by S1 alone, for ARBITRARY superoperators.
ok, detail = True, []
rng_s4 = np.random.default_rng(4155)
for N in (1, 2, 3):
    per_n, biggest, broken_at = True, 0.0, None
    for _ in range(5):
        d2 = 4 ** N
        M = rng_s4.normal(size=(d2, d2)) + 1j * rng_s4.normal(size=(d2, d2))
        a_z, a_y = asym(M, N, "Z"), asym(M, N, "Y")
        if a_y != -a_z:
            per_n, broken_at = False, (a_z, a_y)
        biggest = max(biggest, abs(a_z))
    ok &= per_n and biggest > 0.0
    detail.append(f"N={N} " + ("exact in all 5" if per_n
                               else f"BROKEN at {broken_at}")
                  + f", largest |asymmetry| seen {biggest:.4g} so the equality is "
                    f"not the trivial 0 = -0")
gate("S4", ok,
     "for FIFTEEN arbitrary complex superoperators, none of them of generator "
     "form, asymmetry_Y == -asymmetry_Z holds BIT-EXACTLY (==, not a tolerance): "
     + "; ".join(detail)
     + ". The mechanism: Ad_{Pi_Y} = Ad_{Pi_Z}^3 by S1, and inverting an order-4 "
       "map exchanges its +i and -i eigenprojections while leaving the anti part "
       "as the same array (it lies where Ad^2 = -1, so Ad^3 acts on it as -Ad "
       "does), whence the two squared norms are the SAME two floats reached in "
       "the opposite order. The Y convention adds no law: it is F155 with one "
       "global minus, and by D . Pi_Z . D = Pi_Y that minus is literally the "
       "twisted-versus-untwisted freedom of S3, not merely its analogue")

# ---------------------------------------------------------------- S5
# The Pi^2 grading per letter. A typed claim in compute/ had contradicted this.
ok, detail = True, []
for letter in LETTERS:
    for N in (1, 2, 3):
        P2 = build_pi_full(N, letter) @ build_pi_full(N, letter)
        off = float(np.max(np.abs(P2 - np.diag(np.diag(P2)))))
        diag = np.real(np.diag(P2))
        char_a = np.array([(-1.0) ** (sum(bit_a(i) for i in _k_to_indices(k, N)) % 2)
                           for k in range(4 ** N)])
        char_b = np.array([(-1.0) ** (sum(bit_b(i) for i in _k_to_indices(k, N)) % 2)
                           for k in range(4 ** N)])
        want = char_a if letter == "X" else char_b
        other = char_b if letter == "X" else char_a
        # the match is informative only if the two candidate parities differ
        ok &= (off == 0.0) and np.array_equal(diag, want) \
            and not np.array_equal(want, other)
    detail.append(f"Pi_{letter}^2 = (-1)^{'bit_a' if letter == 'X' else 'bit_b'}")
gate("S5", ok,
     "Pi^2 is exactly diagonal (off-diagonal exactly 0.0) and grades by ONE "
     "parity per letter: " + ", ".join(detail)
     + ", at N=1,2,3, bit-exact, with the two candidate parities confirmed "
       "DISTINCT at each N so a match is informative rather than automatic. This "
       "is F38/F88a and the repo owns it (PiOperator.cs, pi_squared_eigenvalue, "
       "and ANALYTICAL_FORMULAS, which states the letter's Klein pair beside it). "
       "The gate is here deliberately, because the typed claim "
       "DissipatorAxisSelectsPolarity.cs USED TO assert that Pi^2_Y activates "
       "bit_a and bit_b SIMULTANEOUSLY: that sentence confused the Klein index of "
       "the LETTER Y, which is indeed (1,1), with the grading of Pi_Y^2, which is "
       "bit_b alone. The claim did NOT invent the conflation; its anchor "
       "hypotheses/THE_POLARITY_LAYER.md carried it in the same words, in a list "
       "explicitly of polarity axes, and from there it reached a second typed "
       "claim and a live arc. All repaired together and recorded in "
       "docs/CAUGHT_ERRORS.md")

# ---------------------------------------------------------------- S6
# The closed form against the direct build, and the ERROR MODEL as a law.
rng_s6 = np.random.default_rng(6155)
ok, table_rows = True, []
for letter in LETTERS:
    cell, per_n, weakest_signal, exact_hits, total = [], {}, float("inf"), 0, 0
    for N in (2, 3, 4):
        here = []
        for decade in (-4, -2, 0, 2, 4):
            sc = 10.0 ** decade
            A = sc * random_hermitian(rng_s6, N)
            B = sc * random_hermitian(rng_s6, N)
            M = to_pauli(generator(A, B, N), N)
            scale = float(np.sum(np.abs(M) ** 2))
            predicted = closed_form(A, B, N, letter)
            residual = asym(M, N, letter) - predicted
            here.append(eps_ratio(residual, scale))
            exact_hits += 1 if residual == 0.0 else 0
            total += 1
            weakest_signal = min(weakest_signal, eps_ratio(predicted, scale))
        per_n[N] = max(here)
        cell.extend(here)
    hi = max(cell)
    nonzero = [v for v in cell if v > 0.0]
    typical = float(np.median(nonzero)) if nonzero else 0.0
    drift = (hi / typical) if typical > 0.0 else float("inf")
    letter_ok = (hi <= EPS_RATIO_BOUND and drift <= EPS_RATIO_SPREAD
                 and weakest_signal > SIGNAL_GUARD
                 and IMAG_RESIDUE["worst_ratio"] <= EPS_RATIO_BOUND)
    ok &= letter_ok
    table_rows.append(
        f"Pi_{letter}: worst per N "
        + ", ".join(f"N={n} {v:.2f}" for n, v in per_n.items())
        + f"; overall worst {hi:.2f}, typical (median of the nonzero) "
          f"{typical:.2f}, drift {drift:.1f}x, {exact_hits} of {total} instances "
          f"BIT-EXACT, weakest signal {weakest_signal:.2g} floors")
gate("S6", ok,
     "the closed form reproduces the measured asymmetry on dense random Hermitian "
     "A and B at N=2,3,4 with the input scaled over NINE DECADES (10^-4 to 10^4), "
     "fifteen instances per convention. THE ERROR MODEL, stated and not referred "
     "to: the asymmetry is a difference of two sums of squared moduli, each "
     "O(||M||_F^2), so the floor is C*eps*||M||_F^2 with C an O(1) constant "
     "independent of N and of the overall scale. What is gated is that C STAYS "
     "PUT: the worst ratio must not exceed the TYPICAL one (the median of the "
     "nonzero ratios) by more than a bounded factor, so a C that grew with N or "
     "with the scale would fail even while every instance sat under the ceiling. "
     "The per-N worsts are printed so the trend can be read directly, and the "
     "count of BIT-EXACT instances is reported rather than assumed, since some "
     "draws reproduce the closed form to the last bit: "
     + "; ".join(table_rows)
     + f" (bounds: worst <= {EPS_RATIO_BOUND:.0f}, drift <= "
       f"{EPS_RATIO_SPREAD:.0f}x, signal > {SIGNAL_GUARD:.0g} floors so the match "
       "can never be zero against zero). The laws confirmed: Z on bit_b-odd with "
       "(-1)^#Z (F155), X on bit_a-odd with (-1)^#X, Y on bit_b-odd with (-1)^#Y"
     + f". The one quantity DISCARDED anywhere in this file is the imaginary part "
       f"of each Pauli coefficient of a Hermitian operator, zero in exact "
       f"arithmetic and a rounding residue as computed; its worst value across "
       f"the whole run is {IMAG_RESIDUE['worst_ratio']:.2f} floors of eps*||M||_F, "
       f"so it is accounted for rather than assumed away")

# ---------------------------------------------------------------- S7
# Diagonality, and it is EXACT, so it is asserted and not gated.
ok, counts = True, []
for letter in LETTERS:
    N, worst_abs, pairs = 2, 0.0, 0
    for k1 in range(1, 4 ** N):
        for k2 in range(1, 4 ** N):
            if k1 == k2:
                continue
            A = pauli_string(_k_to_indices(k1, N))
            B = pauli_string(_k_to_indices(k2, N))
            M = to_pauli(generator(A, B, N), N)
            worst_abs = max(worst_abs, abs(asym(M, N, letter)))
            pairs += 1
    ok &= (worst_abs == 0.0)
    counts.append(f"Pi_{letter}: {pairs} pairs, worst |asymmetry| {worst_abs!r}")
gate("S7", ok,
     "DIAGONALITY holds in every convention and it is EXACT, so it is compared "
     "with == and not gated against a floor: A one Pauli string and B a DIFFERENT "
     "one gives asymmetry exactly 0.0, over all 210 distinct ordered pairs of "
     "NON-IDENTITY strings at N=2 (" + "; ".join(counts)
     + "). The identity string is excluded, as in the parent's gate G4 and for "
       "the same reason: L_I = R_I makes its commutator lift identically zero, so "
       "its vanishing is trivial for a different reason, and including it would "
       "report 240 pairs for a weaker statement. This is Lemma 3 of the parent "
       "proof, which is Pi-free and therefore letter-free, and it is what makes "
       "each sibling a DIAGONAL bilinear form rather than a general one")

# ---------------------------------------------------------------- S8
# The F113 configuration is SILENT under Pi_X, and so is every bit_a-even one.
rng_s8 = np.random.default_rng(8155)
ok, detail, silent_all, halves_all = True, [], True, True
for trial in range(8):
    N = 3
    omega = rng_s8.normal(size=N)
    g_t1 = np.abs(rng_s8.normal(size=N))
    g_pump = np.abs(rng_s8.normal(size=N))
    A = sum((omega[l] / 2.0) * site_op(N, l, PAULI["Z"]) for l in range(N))
    # B from the collapse operators themselves: iB = -(i/2) sum_k c_k^dag c_k
    c_ops = ([np.sqrt(g_t1[l]) * site_op(N, l, SIGMA_MINUS) for l in range(N)]
             + [np.sqrt(g_pump[l]) * site_op(N, l, SIGMA_PLUS) for l in range(N)])
    B = -0.5 * sum(c.conj().T @ c for c in c_ops)
    B = 0.5 * (B + B.conj().T)
    M = to_pauli(generator(A, B, N), N)
    scale = float(np.sum(np.abs(M) ** 2))
    f113 = (4 ** N / 2.0) * sum(omega[l] * (g_pump[l] - g_t1[l]) for l in range(N))
    a_z = asym(M, N, "Z")
    res_x = polarity(M, N, "X")
    r = eps_ratio(a_z - f113, scale)
    ok &= (r <= EPS_RATIO_BOUND) and eps_ratio(f113, scale) > SIGNAL_GUARD
    silent_all &= (float(res_x["asymmetry"]) == 0.0)
    halves_all &= (float(res_x["norm_sq"]["M_plus_half"]) == 0.0
                   and float(res_x["norm_sq"]["M_minus_half"]) == 0.0)
    if trial == 0:
        res_z = polarity(M, N, "Z")
        detail.append(
            f"first draw: Pi_Z reads {a_z:+.6f} against F113's {f113:+.6f} at "
            f"ratio {r:.2f}, its halves {float(res_z['norm_sq']['M_plus_half']):.2f} "
            f"and {float(res_z['norm_sq']['M_minus_half']):.2f}; Pi_X reads "
            f"{float(res_x['asymmetry'])!r} with halves "
            f"{float(res_x['norm_sq']['M_plus_half'])!r} and "
            f"{float(res_x['norm_sq']['M_minus_half'])!r}")
# the general statement: ANY bit_a-even content is silent, not just this family
general_silent = True
for _ in range(6):
    N = 3
    yy = site_op(N, 0, PAULI["Y"]) @ site_op(N, 1, PAULI["Y"])
    A = 0.7 * yy + rng_s8.normal() * site_op(N, 2, PAULI["Z"])
    B = 1.1 * yy - rng_s8.normal() * site_op(N, 2, PAULI["Z"])
    res = polarity(to_pauli(generator(A, B, N), N), N, "X")
    general_silent &= (float(res["asymmetry"]) == 0.0
                       and float(res["norm_sq"]["M_plus_half"]) == 0.0
                       and float(res["norm_sq"]["M_minus_half"]) == 0.0)
ok &= silent_all and halves_all and general_silent
gate("S8", ok,
     "the F113 configuration at N=3, with B built from the sigma-minus/sigma-plus "
     "COLLAPSE OPERATORS as B = -(1/2) sum_k c_k^dag c_k rather than asserted, "
     "over EIGHT random draws of the drive and of both rate sets: under Pi_Z it "
     "matches F113's (4^N/2)*sum omega*(gamma_pump - gamma_T1) every time. Under "
     "Pi_X it is not balanced, it is SILENT: the asymmetry is exactly 0.0 in all "
     "eight AND both +-i norms are individually exactly 0.0, so there is no "
     "cancellation, there is nothing there (the repo's third answer beside "
     "balanced and broken). " + "; ".join(detail)
     + ". And the silence is not about F113: six further draws of a bit_a-EVEN but "
       "NON-diagonal configuration (a YY bond against a single Z) are silent in "
       "exactly the same way, because any such generator lies in the +1 eigenspace "
       "of Ad_{Pi_X}^2, which makes the anti part the exact zero array at every N "
       "and every rate. The very configuration that makes F155 subsume F113 is "
       "invisible in the X convention")

# ---------------------------------------------------------------- S9
# Positive controls with their own negative half, one per convention.
ok, detail = True, []
for letter, carrier, off_carrier in (("Z", "Z", "X"), ("X", "X", "Z"), ("Y", "Z", "X")):
    N = 3
    M = to_pauli(generator(site_op(N, 0, PAULI[carrier]),
                           site_op(N, 0, PAULI[carrier]), N), N)
    value = asym(M, N, letter)
    expect = -(4 ** (N + 1)) if letter != "Y" else +(4 ** (N + 1))
    ok &= eps_ratio(value - expect, float(np.sum(np.abs(M) ** 2))) <= EPS_RATIO_BOUND
    wrong = asym(to_pauli(generator(site_op(N, 0, PAULI[off_carrier]),
                                    site_op(N, 0, PAULI[off_carrier]), N), N),
                 N, letter)
    ok &= (wrong == 0.0)
    detail.append(f"Pi_{letter}, A=B={carrier}_0: {value:+.1f} (expect {expect:+d}); "
                  f"off-support A=B={off_carrier}_0: {wrong!r}")
gate("S9", ok,
     "POSITIVE CONTROLS, each with its own negative half, and each convention "
     "needs its own carrier because the supports differ: " + "; ".join(detail)
     + ". A single site of the matched letter gives the full 4^(N+1) with the sign "
       "of its weight; the off-support letter gives exactly zero. The Y control "
       "sits on the SAME carrier as Z with the opposite sign, which is S4 again. "
       "The expected values here are literals and owe nothing to closed_form, so "
       "this gate stands even if the assembly is wrong, and it is one of the two "
       "that caught a corrupted Pi during this file's audit")

# ---------------------------------------------------------------- S10
# The transport route: X is Z in the Hadamard-rotated frame, all sites rotated.
rng_s10 = np.random.default_rng(10155)
ok, worst, weakest, partial_differs = True, 0.0, float("inf"), True
for N in (2, 3):
    h, h_partial = hadamard_full(N), hadamard_first_site_only(N)
    for _ in range(4):
        A, B = random_hermitian(rng_s10, N), random_hermitian(rng_s10, N)
        M = to_pauli(generator(A, B, N), N)
        M_rot = to_pauli(generator(h @ A @ h, h @ B @ h, N), N)
        M_part = to_pauli(generator(h_partial @ A @ h_partial,
                                    h_partial @ B @ h_partial, N), N)
        scale = float(np.sum(np.abs(M) ** 2))
        a_x = asym(M, N, "X")
        worst = max(worst, eps_ratio(a_x - asym(M_rot, N, "Z"), scale))
        weakest = min(weakest, eps_ratio(a_x, scale))
        partial_differs &= (eps_ratio(a_x - asym(M_part, N, "Z"), scale)
                            > SIGNAL_GUARD)
ok = (worst <= EPS_RATIO_BOUND and weakest > SIGNAL_GUARD and partial_differs)
gate("S10", ok,
     f"TRANSPORT, the second route to the X law and the only gate here that never "
     f"touches the swap table: asymmetry_X(A, B) equals asymmetry_Z(hAh, hBh) for "
     f"h = H^(tensor N), with no sign and no scale correction, worst ratio "
     f"{worst:.2f} over eight dense random Hermitian pairs at N=2,3, weakest "
     f"compared signal {weakest:.2g} floors so the match is never zero against "
     f"zero. It is a ROUTE and not a coincidence because the repo owns the "
     f"intertwiner Q_zx . Pi_Z . Q_zx^-1 = Pi_X "
     f"(PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md section (d), the same Hadamard "
     f"route it uses for the F112 siblings), so the two readings differ only in "
     f"which frame the Pauli content is expressed in. EVERY site must be rotated: "
     f"the same test with a Hadamard on site 0 alone disagrees macroscopically "
     f"(disagrees: {partial_differs}), because a partial rotation yields a "
     f"MIXED-letter law and not the X one. Together with S4 this is the file's "
     f"headline: neither sibling is a new law, one is D-conjugation and the other "
     f"Hadamard conjugation, and both operators were already in the repo")

# ----------------------------------------------------------------------
passed = sum(1 for _, p, _ in GATES if p)
print()
print("=" * 78)
print(f"{passed}/{len(GATES)} gates passed")
print("=" * 78)
if passed != len(GATES):
    raise AssertionError(
        "F155 sibling gates failed: "
        + ", ".join(n for n, p, _ in GATES if not p))
