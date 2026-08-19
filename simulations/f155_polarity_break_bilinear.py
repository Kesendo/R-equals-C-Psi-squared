"""F155 gate: the polarity break of the physical generator, in closed form.

THE CLAIM. Write the no-jump generator of a non-Hermitian H = A + iB, with A and
B Hermitian (the canonical split, always available and unique), as

    G_H : rho -> -i(H rho - rho H^dag) = -i[A, rho] + {B, rho}

and expand both parts over Pauli strings, A = sum_s a_s s, B = sum_s b_s s (real
coefficients). Then the F112 polarity asymmetry of G_H is

    asymmetry = 4^(N+1) * sum_{s : bit_b(s) odd} (-1)^(#Z(s)) * a_s * b_s

a BILINEAR form that is DIAGONAL in the Pauli basis: no pair of distinct strings
contributes, and only the bit_b-odd strings carry anything at all. bit_b(s) is
#Y(s) + #Z(s) mod 2, the F38 parity, which the repo also owns as the character of
conjugation by X^(tensor N), the global spin flip (docs/GLOSSARY.md, the polarity
cube). On the odd support (-1)^#Z is minus the transpose parity (-1)^#Y, the
cube's third axis, so the law reads: only the spin-flip-ODD content of A and B
pairs, each term weighted by its transpose parity; the sign of that weight is minus
in the pairing this file computes in and plus in the untwisted one, see WHICH
PAIRING below.

WHICH PAIRING, and this is not decoration. The framework reads a row-stack
generator against the order='F' Pauli transform. That pairing conjugates the bare
multiplication superoperators by D = diag((-1)^#Y), F118's transpose
superoperator, which is the same as conjugating Pi. F113's registered numbers live
in this pairing and so does everything below. In the untwisted pairing the swap
rule's two signs exchange and the law carries (-1)^#Y instead, i.e. the whole
value flips sign on the bit_b-odd support. Both spellings are gated (G9). What is
convention-free: the magnitude, the zero-versus-nonzero verdict, and the RATIO of
signs between any two odd strings; what is not: the attachment to #Z rather than
to #Y. F113's entry pins the same freedom for its own number.

WHAT IT SUBSUMES. F113 is the Z-diagonal special case of the no-jump face: for
H = sum_l (omega_l/2) Z_l and iB = -(i/2) sum_k c_k^dag c_k with sigma^- at
gamma_T1,l and sigma^+ at gamma_pump,l, one has a_{Z_l} = omega_l/2 and
b_{Z_l} = (gamma_T1,l - gamma_pump,l)/4, and the formula returns
(4^N/2) * sum_l omega_l * (gamma_pump,l - gamma_T1,l), F113 exactly (G7). F113's
own theorem is stated for the FULL Lindbladian; that the recycling term leaves the
value alone is a separate MEASURED fact (F113 entry, N = 2, 3, 4, sigma-minus/plus
family), not derived here. F112's balance for the commutator of an ARBITRARY
matrix also falls out, in one line (G12).

THE PROOF, and every step of it is gated here.

  Step 1 (G1). For any superoperator M and the order-4 palindrome operator Pi,
      asymmetry = ||P_{+i}M||^2 - ||P_{-i}M||^2 = Im<M, Pi M Pi^-1>,
  with P_{+-i} the eigenprojections of Ad_Pi. An identity, no hypothesis: Ad_Pi is
  Frobenius-unitary, so <M, Ad_Pi^3 M> is the conjugate of <M, Ad_Pi M> and
  P_{+i} - P_{-i} = (i/2)(Ad_Pi^3 - Ad_Pi).

  Step 2 (G2). Pi factorises over sites, Pi_N = pi^(tensor N).

  Step 3, THE SWAP RULE (G3a, G3b, G9). Ad_Pi EXCHANGES the two multiplications
  and pays a sign for it. In the pipeline pairing, per site and bit-exact:

      Ad(L_p) = eps_L(p) * R_p ,  Ad(R_p) = eps_R(p) * L_p
      (eps_L, eps_R) = (+1,+1) for I and X, (-1,+1) for Y, (+1,-1) for Z

  hence Ad(L_s) = (-1)^#Y(s) * R_s and Ad(R_s) = (-1)^#Z(s) * L_s. In the
  untwisted pairing the two exchange. Applying the rule twice gives
  eps_L*eps_R = (-1)^bit_b, the sign pattern of F38.

  PARENTAGE, because this is a refinement and not a discovery. F118 factors
  Pi_Z = R . D with D the transpose (docs/ANALYTICAL_FORMULAS.md F118); an
  antiautomorphism is exactly what exchanges left with right multiplication, and
  docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md (F119) records the sibling swap for
  F118's edge mirrors, K R K = F R, i.e. rho.F <-> F.rho. F114 (typed,
  compute/RCPsiSquared.Core/Symmetry/CommutatorDConjugationSign.cs) owns the
  commutator combination D L_s D = (-1)^(n_Y+1) L_s, and
  docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md section (f) owns the trichotomy
  "Pi_Z L_s Pi_Z^-1 = -L_s / +L_s / an anticommutator superoperator, neither".
  What is new here is the PER-SIDE form with both signs, which pins the sign that
  trichotomy leaves open, and what follows from it below.

  Step 4, ORTHOGONALITY AND ASSEMBLY (G3c, G4, G5). The lifts are orthogonal:
  <L_s, L_t> = <R_s, R_t> = 4^N delta_st, and <L_s, R_t> = Tr(s)Tr(t), which
  vanishes unless both strings are the identity. So only s = t survives. With the
  swap rule the same-string terms are

      commutator with commutator:      -a_s^2 * 4^N * (E_Z + E_Y)      REAL
      anticommutator with itself:      +b_s^2 * 4^N * (E_Z + E_Y)      REAL
      the two cross terms, summed:  2i * a_s b_s * 4^N * (E_Z - E_Y)

  with E_Z = (-1)^#Z and E_Y = (-1)^#Y. Taking imaginary parts drops both
  self-terms: the first IS F112's balance for the commutator part, the second is
  the anticommutator part's balance, measured on 2026-08-19 and derived here. The
  cross term is 0 for bit_b-even s (the signs agree) and 4^(N+1)*(-1)^#Z*a_s*b_s
  for bit_b-odd s (they are opposite). The identity string needs no special
  handling in the result, being bit_b-even, but note that L_I = R_I, so its
  commutator lift is identically zero and <L_I, R_I> = 4^N is the one nonvanishing
  L-R overlap; both facts move REAL terms only.

COROLLARIES, all gated:
  a) B = 0 (Hermitian H): balanced (G6).
  b) A carries only bit_b-EVEN strings: balanced against EVERY B (G6).
  c) B proportional to I: balanced against every A (G6).
  d) A and B carry their bit_b-odd content on disjoint strings: balanced (G6).
     So bit_b-odd content in A is NECESSARY and NOT SUFFICIENT.
  e) Adding a multiple of the identity to A or to B changes nothing (G10), so the
     law is BLIND to the gain/loss offset: balanced gain/loss and pure loss of the
     same imbalance read the same.
  f) The textbook PT dimer J(XX+YY)/2 + i g (Z_0 - Z_1)/2 is balanced at every g,
     because XX and YY are bit_b-even; g = J = 1 is where its single-excitation
     block sits at its exceptional point and the asymmetry does not notice. Add a
     0.4*Z_0 detuning and the value is exactly -12.8*g, a straight line with no
     feature at g = 1, and that detuned system has no EP on the real g axis at all
     (G11). The object is a detuning x gain-imbalance meter, not a PT meter.
  g) The transpose-paired lift K (x) I + I (x) K^T is balanced for ANY K, not only
     Hermitian: split K = K_1 + i K_2 into Hermitian parts; the self-terms are real
     as above and the two cross terms are symmetric in the pair, so they cancel
     (G12). This makes the drain's balance a theorem for every collapse operator,
     c^dag c being Hermitian anyway.
  h) The commutator of an ARBITRARY matrix is balanced (G12). That is the
     AGGREGATE face of F112's non-Hermitian extension, not that theorem: F112
     proves the per-pair identity for every pair of strings, and diagonality means
     the cross-pair information is exactly what drops out here.

SCOPE. Pi in the Z-dephasing convention. The X- and Y-dephase siblings are now
CLOSED, in section (g) of the proof and by simulations/f155_dephase_siblings.py,
and they are NOT two new laws: Y is this law with one global minus (Pi_Y = Pi_Z^-1,
so the +-i projections exchange, universally in M) and X is this law on the
bit_a-odd support with weight (-1)^#X, equivalently this law in the Hadamard
rotated frame. This block used to call them 'the Klein-V4 axis substitution',
which lumps Y in with X; only X moves the axis. The value is NOT basis-invariant: it
moves under a per-site Hadamard and under a local rotation (G13), so like F113's
number it must be quoted in the frame its Pauli content is defined in. The full
Lindbladian is a different operator, see WHAT IT SUBSUMES.

STAGE 0, what the repo already held (swept 2026-08-19 by agents over docs/proofs/,
docs/ANALYTICAL_FORMULAS.md, experiments/, hypotheses/ incl. archive, recovered/,
reflections/, compute/ incl. the OpenArcs and Confirmations registries,
docs/CAUGHT_ERRORS.md, simulations/framework/):
  * The parity half of Step 3 is OWNED and PROVEN: "L_s has zero Pi-conjugation
    +-i content whenever s is bit_b-even", section (b) of
    docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md, universal N, from F38.
  * The swap rule's PARENTS are owned and typed: F118, F119, F114 and
    PROOF_PI_FACTORS_AS_R_TIMES_D section (f). See PARENTAGE above. The first
    sweep of this session asked about the consequence and not about the mechanism,
    and so did not surface them; a second sweep did.
  * The hypothesis-free half of the bilinear expansion is written down in
    docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md.
  * The anticommutator lift's balance was MEASURED on 2026-08-19
    (hypotheses/THE_DRAIN_HAS_NO_CHIRALITY.md, gate G9b of
    simulations/f112_gain_loss_carrier_check.py) and never derived.
  * NOT the same object, and refuted: reflections/POLARITY_COORDINATES.md's
    conjecture that "A (x) B* + B (x) A* is automatically balanced" was killed at
    240 of 240 configurations. That is the CONJUGATE pairing; the lifts here are
    the TRANSPOSE pairing, and the conjugate arm is the positive control of gate
    G9b in the sibling checker.
  * fw.Confirmations and the C# ConfirmationsRegistry returned NOTHING: no
    hardware confirmation touches the polarity balance of the physical generator.

TOLERANCES. G2, G3a, G3b and G3c compare exactly, because an exact route exists
(signed permutations and integers). Everything that goes through the float
polarity projection is gated against the noise floor eps*||M||_F^2, the error
model established in f112_gain_loss_carrier_check.py, where the ratio and not the
residual is what stays put as the scale moves.
"""
from __future__ import annotations

import itertools
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
    _build_bilinear,
    _k_to_indices,
    _pauli_label,
    _vec_to_pauli_basis_transform,
    pauli_string,
)
from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L

EPS = float(np.finfo(float).eps)
EPS_RATIO_BOUND = 8.0
GATES = []
_CACHE = {}

I2 = np.eye(2, dtype=complex)
PAULI = {
    "I": I2,
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)
SIGMA_PLUS = np.array([[0, 0], [1, 0]], dtype=complex)


def gate(name, passed, detail):
    GATES.append((name, bool(passed), detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}")
    return bool(passed)


def basis(N):
    if N not in _CACHE:
        Pi = build_pi_full(N)
        _CACHE[N] = (_vec_to_pauli_basis_transform(N), Pi, np.linalg.inv(Pi))
    return _CACHE[N]


def to_pauli(L_vec, N):
    T, _, _ = basis(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def ad_pi(M, N):
    _, Pi, Pinv = basis(N)
    return Pi @ M @ Pinv


def im_form(M, N):
    """Im<M, Ad_Pi M>, the asymmetry by Step 1."""
    return float(np.imag(np.sum(np.conj(M) * ad_pi(M, N))))


def asym_and_scale(M, N):
    _, Pi, _ = basis(N)
    res = polarity_coordinates_from_L(M, N, 0.0, Pi=Pi)
    return float(res["asymmetry"]), float(np.sum(np.abs(M) ** 2))


def site_op(N, l, op2):
    ops = [I2] * N
    ops[l] = op2
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def generator(A, B, N):
    """-i(H rho - rho H^dag) in vec form, H = A + iB with A, B Hermitian."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    return (-1j * (np.kron(A, Id) - np.kron(Id, A.T))
            + (np.kron(B, Id) + np.kron(Id, B.T)))


def pauli_coeffs(M, N):
    d = 2 ** N
    return [float(np.real(np.trace(pauli_string(_k_to_indices(k, N)) @ M) / d))
            for k in range(4 ** N)]


def closed_form(A, B, N):
    ca, cb = pauli_coeffs(A, N), pauli_coeffs(B, N)
    total = 0.0
    for k in range(1, 4 ** N):
        label = _pauli_label(k, N)
        if sum(1 for c in label if c in "YZ") % 2 == 0:
            continue
        total += ((-1) ** sum(1 for c in label if c == "Z")) * ca[k] * cb[k]
    return 4 ** (N + 1) * total


def eps_ratio(residual, scale):
    return abs(residual) / (EPS * scale) if scale > 0.0 else 0.0


print("=" * 78)
print("F155 GATE: the polarity break of the physical generator is a DIAGONAL")
print("bilinear form on the bit_b-odd Pauli strings")
print("=" * 78)

# --- G1: the identity behind everything ------------------------------------
rng = np.random.default_rng(20260819)
worst = 0.0
for N in (1, 2, 3):
    D = 4 ** N
    for _ in range(3):
        M = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
        a, scale = asym_and_scale(M, N)
        worst = max(worst, eps_ratio(a - im_form(M, N), scale))
gate("G1 asymmetry == Im<M, Ad_Pi M> for arbitrary M (Step 1)",
     worst <= EPS_RATIO_BOUND,
     f"worst |difference|/(eps*||M||_F^2) = {worst:.2f} over 9 random superoperators, "
     f"N = 1, 2, 3 (an identity, not a property of the input)")

# --- G2: Pi factorises over sites ------------------------------------------
pi1 = build_pi_full(1)
ok = True
detail = []
for N in (2, 3, 4):
    prod = pi1
    for _ in range(N - 1):
        prod = np.kron(prod, pi1)
    diff = float(np.max(np.abs(build_pi_full(N) - prod)))
    ok = ok and diff == 0.0
    detail.append(f"N={N}: {diff!r}")
gate("G2 Pi_N == pi^(tensor N), bit for bit (Step 2)", ok,
     "max|Pi_N - pi^tensor| " + ", ".join(detail) + " (exact, both are signed permutations)")

# --- G3: the swap rule, per site and per string -----------------------------
T1, _, _ = basis(1)


def lift(p_matrix, side, N):
    """The PIPELINE's lift: a row-stack multiplication superoperator read against
    the order='F' Pauli transform. This is NOT the bare left/right multiplication
    superoperator: it is that object conjugated by D = diag((-1)^#Y), the
    transpose superoperator of F118 (gate G9 pins the relation). The framework's
    asymmetry, F113's registered numbers and this script all live in this pairing,
    which is why the law below carries (-1)^#Z; in the untwisted pairing the two
    spellings exchange and the sign flips on the odd support."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    v = np.kron(p_matrix, Id) if side == "L" else np.kron(Id, p_matrix.T)
    return to_pauli(v, N)


EPS_L = {"I": +1, "X": +1, "Y": -1, "Z": +1}
EPS_R = {"I": +1, "X": +1, "Y": +1, "Z": -1}

per_site_ok = True
detail = []
for p in "IXYZ":
    L = lift(PAULI[p], "L", 1)
    R = lift(PAULI[p], "R", 1)
    ok_l = np.array_equal(ad_pi(L, 1), EPS_L[p] * R)
    ok_r = np.array_equal(ad_pi(R, 1), EPS_R[p] * L)
    per_site_ok = per_site_ok and ok_l and ok_r
    detail.append(f"{p}:({EPS_L[p]:+d},{EPS_R[p]:+d})")
gate("G3a per-site swap rule in the PIPELINE pairing (Step 3)",
     per_site_ok,
     "bit-exact for all four letters, (eps_L, eps_R) = " + ", ".join(detail)
     + "; note eps_L*eps_R = (-1)^bit_b, which re-derives F38")

string_ok = True
worst_str = 0.0
for N in (2, 3):
    for k in range(4 ** N):
        label = _pauli_label(k, N)
        sig = pauli_string(_k_to_indices(k, N))
        L = lift(sig, "L", N)
        R = lift(sig, "R", N)
        sl = (-1) ** sum(1 for c in label if c == "Y")
        sr = (-1) ** sum(1 for c in label if c == "Z")
        worst_str = max(worst_str,
                        float(np.max(np.abs(ad_pi(L, N) - sl * R))),
                        float(np.max(np.abs(ad_pi(R, N) - sr * L))))
gate("G3b string swap rule in the PIPELINE pairing: Ad(L_s) = (-1)^#Y R_s, Ad(R_s) = (-1)^#Z L_s",
     worst_str == 0.0,
     f"worst entry-wise deviation over all 4^N strings at N = 2 and 3: {worst_str!r} "
     f"(exact: signed permutations tensored, so an exact route exists and is taken)")

# --- G3c: the lifts are orthogonal, which is what kills distinct strings -----
worst_ll = worst_rr = worst_lr = 0.0
N = 2
for ka in range(4 ** N):
    La = lift(pauli_string(_k_to_indices(ka, N)), "L", N)
    Ra = lift(pauli_string(_k_to_indices(ka, N)), "R", N)
    for kb in range(4 ** N):
        Lb = lift(pauli_string(_k_to_indices(kb, N)), "L", N)
        Rb = lift(pauli_string(_k_to_indices(kb, N)), "R", N)
        expect = (4 ** N) if ka == kb else 0.0
        worst_ll = max(worst_ll, abs(np.sum(np.conj(La) * Lb) - expect))
        worst_rr = max(worst_rr, abs(np.sum(np.conj(Ra) * Rb) - expect))
        expect_lr = (4 ** N) if (ka == 0 and kb == 0) else 0.0
        worst_lr = max(worst_lr, abs(np.sum(np.conj(La) * Rb) - expect_lr))
gate("G3c lift orthogonality <L_s,L_t> = <R_s,R_t> = 4^N delta, <L_s,R_t> = Tr(s)Tr(t)",
     max(worst_ll, worst_rr, worst_lr) == 0.0,
     f"worst deviations at N = 2 over all 256 string pairs: LL {worst_ll!r}, "
     f"RR {worst_rr!r}, LR {worst_lr!r} (exact; this is what makes the form diagonal)")

# --- G4: diagonality in the Pauli basis, string by string -------------------
for N in (2, 3):
    worst_off = 0.0
    worst_diag_dev = 0.0
    n_off = n_diag = 0
    for ka in range(1, 4 ** N):
        A = pauli_string(_k_to_indices(ka, N))
        for kb in range(1, 4 ** N):
            B = pauli_string(_k_to_indices(kb, N))
            M = to_pauli(generator(A, B, N), N)
            a, scale = asym_and_scale(M, N)
            if ka == kb:
                n_diag += 1
                worst_diag_dev = max(worst_diag_dev,
                                     eps_ratio(a - closed_form(A, B, N), scale))
            else:
                n_off += 1
                worst_off = max(worst_off, eps_ratio(a, scale))
    gate(f"G4 N={N} the form is DIAGONAL in the Pauli basis",
         worst_off <= EPS_RATIO_BOUND and worst_diag_dev <= EPS_RATIO_BOUND,
         f"{n_off} distinct-string pairs: worst ratio {worst_off:.2f}; "
         f"{n_diag} same-string pairs match the closed form to {worst_diag_dev:.2f}")

# --- G5: the closed form on dense random A, B -------------------------------
for N in (2, 3, 4):
    d = 2 ** N
    worst_dev = 0.0
    for _ in range(4):
        m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        A = (m + m.conj().T) / 2
        m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        B = (m + m.conj().T) / 2
        M = to_pauli(generator(A, B, N), N)
        a, scale = asym_and_scale(M, N)
        worst_dev = max(worst_dev, eps_ratio(a - closed_form(A, B, N), scale))
    gate(f"G5 N={N} closed form vs direct build, dense random Hermitian A and B",
         worst_dev <= EPS_RATIO_BOUND,
         f"worst |measured - formula|/(eps*||M||_F^2) = {worst_dev:.2f} over 4 pairs")

# --- G6: the corollaries ----------------------------------------------------
for N in (2, 3):
    d = 2 ** N
    m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    B_rand = (m + m.conj().T) / 2
    m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    A_rand = (m + m.conj().T) / 2

    a_herm, s_herm = asym_and_scale(to_pauli(generator(A_rand, np.zeros((d, d), complex), N), N), N)
    even = [k for k in range(1, 4 ** N)
            if sum(1 for c in _pauli_label(k, N) if c in "YZ") % 2 == 0]
    A_even = sum(rng.normal() * pauli_string(_k_to_indices(k, N)) for k in even)
    a_even, s_even = asym_and_scale(to_pauli(generator(A_even, B_rand, N), N), N)
    a_id, s_id = asym_and_scale(to_pauli(generator(A_rand, 0.7 * np.eye(d, dtype=complex), N), N), N)

    odd = [k for k in range(1, 4 ** N)
           if sum(1 for c in _pauli_label(k, N) if c in "YZ") % 2 == 1]
    half = len(odd) // 2
    A_dis = sum(rng.normal() * pauli_string(_k_to_indices(k, N)) for k in odd[:half])
    B_dis = sum(rng.normal() * pauli_string(_k_to_indices(k, N)) for k in odd[half:])
    a_dis, s_dis = asym_and_scale(to_pauli(generator(A_dis, B_dis, N), N), N)

    # the named instance of corollary 4, because a reflection cites it in prose:
    # Probe 1B's own pair, A on XY strings and B on XZ strings, disjoint odd supports
    bonds = [(i, i + 1) for i in range(N - 1)]
    A_xy = _build_bilinear(N, bonds, [("X", "Y", 1.0)])
    B_xz = _build_bilinear(N, bonds, [("X", "Z", 1.0)])
    a_probe1b, s_probe1b = asym_and_scale(to_pauli(generator(A_xy, B_xz, N), N), N)

    ratios = [eps_ratio(a_herm, s_herm), eps_ratio(a_even, s_even),
              eps_ratio(a_id, s_id), eps_ratio(a_dis, s_dis),
              eps_ratio(a_probe1b, s_probe1b)]
    gate(f"G6 N={N} the balanced corollaries (Hermitian H; bit_b-even A; B ~ I; "
         f"disjoint odd supports, twice)",
         max(ratios) <= EPS_RATIO_BOUND,
         "worst ratios " + ", ".join(f"{r:.2f}" for r in ratios)
         + " (each balanced; the third and fourth show that bit_b-odd content in A "
           "is NOT sufficient, and the fifth is the XY-chain against XZ-chain pair "
           "that reflections/POLARITY_COORDINATES.md names in prose)")

# --- G7: F113 is the Z-diagonal special case --------------------------------
for N in (2, 3):
    d = 2 ** N
    omega = rng.uniform(-1.0, 1.0, N)
    g_t1 = rng.uniform(0.0, 0.3, N)
    g_pump = rng.uniform(0.0, 0.3, N)
    A = sum((omega[l] / 2) * site_op(N, l, PAULI["Z"]) for l in range(N))
    c_ops = ([np.sqrt(g_t1[l]) * site_op(N, l, SIGMA_MINUS) for l in range(N)]
             + [np.sqrt(g_pump[l]) * site_op(N, l, SIGMA_PLUS) for l in range(N)])
    B = -0.5 * sum(c.conj().T @ c for c in c_ops)
    M = to_pauli(generator(A, B, N), N)
    a, scale = asym_and_scale(M, N)
    f113 = (4 ** N / 2) * sum(omega[l] * (g_pump[l] - g_t1[l]) for l in range(N))
    dev = max(eps_ratio(a - f113, scale), eps_ratio(closed_form(A, B, N) - f113, scale))
    gate(f"G7 N={N} F113 is the Z-diagonal special case",
         dev <= EPS_RATIO_BOUND,
         f"measured {a:+.10f}, F113 {f113:+.10f}, F155 {closed_form(A, B, N):+.10f}; "
         f"worst ratio {dev:.2f} (random per-site omega, gamma_T1 and gamma_pump)")

# --- G8: the positive control ----------------------------------------------
N = 3
d = 2 ** N
A = site_op(N, 0, PAULI["Z"])
B = site_op(N, 0, PAULI["Z"])
M = to_pauli(generator(A, B, N), N)
a, scale = asym_and_scale(M, N)
gate("G8 positive control: this gate can see a nonzero",
     abs(a) > 1e6 * EPS * scale and eps_ratio(a - closed_form(A, B, N), scale) <= EPS_RATIO_BOUND,
     f"A = B = Z_0 at N = 3 gives {a:+.1f}, the formula gives "
     f"{closed_form(A, B, N):+.1f} (= -4^(N+1) since #Z = 1)")

# --- G9: the two representations, and which one the law is written in --------
# The framework reads a row-stack generator against the order='F' Pauli
# transform. That pairing is the basis-free representation
# S[k,l] = Tr(sigma_k S(sigma_l))/2^N conjugated by C = diag((-1)^#Y), and
# conjugating Pi by C sends it to conj(Pi). So the swap rule and the law come in
# two mirrored spellings, #Y <-> #Z, differing by an overall sign on the
# bit_b-odd support. F113's entry pins the same freedom for its own number.
def free_rep(sig, side, N):
    S = [pauli_string(_k_to_indices(k, N)) for k in range(4 ** N)]
    d = 2 ** N
    act = (lambda X: sig @ X) if side == "L" else (lambda X: X @ sig)
    return np.array([[np.trace(S[a].conj().T @ act(S[b])) / d
                      for b in range(4 ** N)] for a in range(4 ** N)])


rep_ok = True
free_rule_ok = True
for N in (1, 2):
    Cmat = np.diag([(-1) ** _pauli_label(k, N).count("Y")
                    for k in range(4 ** N)]).astype(complex)
    Cinv = np.linalg.inv(Cmat)
    _, Pi, Pinv = basis(N)
    for k in range(4 ** N):
        sig = pauli_string(_k_to_indices(k, N))
        label = _pauli_label(k, N)
        for side in "LR":
            F = free_rep(sig, side, N)
            rep_ok = rep_ok and np.allclose(lift(sig, side, N), Cmat @ F @ Cinv)
        Lf, Rf = free_rep(sig, "L", N), free_rep(sig, "R", N)
        free_rule_ok = free_rule_ok and np.allclose(
            Pi @ Lf @ Pinv, ((-1) ** label.count("Z")) * Rf)
        free_rule_ok = free_rule_ok and np.allclose(
            Pi @ Rf @ Pinv, ((-1) ** label.count("Y")) * Lf)
gate("G9 the two pairings, and the #Y <-> #Z mirroring between them",
     rep_ok and free_rule_ok,
     "the pipeline representation is the basis-free one conjugated by "
     "C = diag((-1)^#Y), and in the basis-free one the swap rule reads "
     "Ad(L) = (-1)^#Z R, Ad(R) = (-1)^#Y L, the two spellings exchanged; both "
     "checked at N = 1, 2 (the law then carries (-1)^#Y there, i.e. the overall "
     "sign flips on the odd support, which is F113's registered convention freedom)")

# --- G10: invariance under identity shifts ----------------------------------
worst_shift = 0.0
for N in (2, 3):
    d = 2 ** N
    m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    A = (m + m.conj().T) / 2
    m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    B = (m + m.conj().T) / 2
    base, scale = asym_and_scale(to_pauli(generator(A, B, N), N), N)
    for shifted in (generator(A + 3.7 * np.eye(d), B, N),
                    generator(A, B + 3.7 * np.eye(d), N)):
        a, _ = asym_and_scale(to_pauli(shifted, N), N)
        worst_shift = max(worst_shift, eps_ratio(a - base, scale))
gate("G10 identity shifts of A or of B change nothing",
     worst_shift <= EPS_RATIO_BOUND,
     f"worst ratio {worst_shift:.2f} (the identity string is bit_b-even, so it "
     f"drops; physically the law is BLIND to the gain/loss offset, and balanced "
     f"gain/loss reads the same as pure loss of the same imbalance)")

# --- G11: the textbook PT dimer, across its exceptional point ----------------
N = 2
H_bond = _build_bilinear(N, [(0, 1)], [("X", "X", 0.5), ("Y", "Y", 0.5)])
Z0 = site_op(N, 0, PAULI["Z"])
Z1 = site_op(N, 1, PAULI["Z"])
pt_ratios = []
det_vals = []
gs = (0.0, 0.5, 0.99, 1.0, 1.5, 3.0)
for g in gs:
    B = (g / 2) * (Z0 - Z1)
    a_pt, s_pt = asym_and_scale(to_pauli(generator(H_bond, B, N), N), N)
    a_det, _ = asym_and_scale(to_pauli(generator(H_bond + 0.4 * Z0, B, N), N), N)
    pt_ratios.append(eps_ratio(a_pt, s_pt))
    det_vals.append(a_det)
linear = all(abs(det_vals[i] - g * (-12.8)) < 1e-9 for i, g in enumerate(gs))
gate("G11 the textbook PT dimer is balanced through its exceptional point",
     max(pt_ratios) <= EPS_RATIO_BOUND and linear,
     f"H = J(XX+YY)/2 + i*g*(Z_0 - Z_1)/2 gives 0 at g = 0 .. 3 (worst ratio "
     f"{max(pt_ratios):.2f}), because XX and YY are bit_b-EVEN; g = J = 1 is where "
     f"this dimer's single-excitation block sits at its exceptional point, and the "
     f"asymmetry does not notice. Add a 0.4*Z_0 detuning and the value is exactly "
     f"-12.8*g, a straight line with no feature at g = 1 (and that detuned system "
     f"has no EP on the real g axis at all). A detuning x gain-imbalance meter, "
     f"not a PT meter")

# --- G12: two statements beyond the theorem's own hypothesis ----------------
worst_k = 0.0
worst_comm = 0.0
for N in (2, 3):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    for _ in range(3):
        K = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))   # NOT Hermitian
        a, scale = asym_and_scale(to_pauli(np.kron(K, Id) + np.kron(Id, K.T), N), N)
        worst_k = max(worst_k, eps_ratio(a, scale))
        H = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))   # arbitrary matrix
        a, scale = asym_and_scale(to_pauli(-1j * (np.kron(H, Id) - np.kron(Id, H.T)), N), N)
        worst_comm = max(worst_comm, eps_ratio(a, scale))
gate("G12 the two corollaries beyond the Hermitian hypothesis",
     worst_k <= EPS_RATIO_BOUND and worst_comm <= EPS_RATIO_BOUND,
     f"transpose-paired lift with NON-Hermitian K: worst ratio {worst_k:.2f} "
     f"(the polarisation step); commutator of an ARBITRARY matrix: "
     f"{worst_comm:.2f} (the aggregate face of F112's non-Hermitian extension, "
     f"not its per-pair theorem, which diagonality cannot reach)")

# --- G13: the scope fact a reader most needs --------------------------------
N = 2
d = 4
m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
A = (m + m.conj().T) / 2
m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
B = (m + m.conj().T) / 2
base, _ = asym_and_scale(to_pauli(generator(A, B, N), N), N)
had = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
U = np.kron(had, had)
rot, _ = asym_and_scale(to_pauli(generator(U @ A @ U.conj().T, U @ B @ U.conj().T, N), N), N)
theta = 0.3
ry = np.array([[np.cos(theta / 2), -np.sin(theta / 2)],
               [np.sin(theta / 2), np.cos(theta / 2)]], dtype=complex)
V = np.kron(ry, np.eye(2, dtype=complex))
loc, _ = asym_and_scale(to_pauli(generator(V @ A @ V.conj().T, V @ B @ V.conj().T, N), N), N)
gate("G13 the value is NOT basis-invariant, and that is the scope fact",
     abs(rot - base) > 1e-6 and abs(loc - base) > 1e-6,
     f"base {base:+.4f}, after a per-site Hadamard (X <-> Z relabel) {rot:+.4f}, "
     f"after a local Y-rotation by 0.3 {loc:+.4f}: the law reads the Pauli content "
     f"of A and B in the frame the Z-palindrome convention fixes, so it must be "
     f"quoted in that frame, as F113's number must be quoted in the frame its "
     f"drive is defined in")

print("\n" + "=" * 78)
failed = [n for n, ok, _ in GATES if not ok]
print(f"{len(GATES) - len(failed)}/{len(GATES)} gates passed")
if failed:
    print("FIRED: " + ", ".join(failed))
print("=" * 78)
assert not failed, f"gates fired: {failed}"
