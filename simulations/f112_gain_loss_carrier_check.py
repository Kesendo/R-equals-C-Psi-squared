"""What decides the F112 polarity balance of the physical generator (checker)

WHY THIS EXISTS. `docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md` used to
name "the named PT/gain-loss spin systems" in its list of configurations whose
physical generator −i(Hρ − ρH†) is balanced (asym = 0). That was wrong, and the
first repair of it was wrong in the mirror-image way ("PT systems are not in the
balanced class"). Both errors name a CLASS OF SYSTEMS where the criterion is a
content condition. This script is the from-below check behind the repaired
sentences, and behind the "the recycling jump is asymmetry-inert" line added to
the F113 registry entry the same day.

THE LAW BEING CHECKED is F113 (Tier 1 derived, `docs/ANALYTICAL_FORMULAS.md`),
and it is SAME-SITE BILINEAR:

    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l),   H ⊃ Σ_l (ω_l/2)·Z_l

Each site contributes its own product and cross-site pairs contribute nothing
(registered 2026-05-26, `experiments/F113_BREAK_MAGNITUDE_FORMULA.md` item 5).
So the balance question is never about a class of systems: the sum is zero when
every term is zero (no Z anywhere, or detailed balance everywhere, or carrier and
imbalance simply sitting on different sites) and also when nonzero terms cancel
across sites. Gates G0-G2 and G6-G8 walk those cases.

SCOPE. The dissipator in G0-G8 is F113's own family: each c_k a per-site σ⁻ or
σ⁺. That family is OUTSIDE F112's bit_b-homogeneous hypothesis (σ⁻ = (X + iY)/2
has support on X, bit_b = 0, and Y, bit_b = 1, so it is bit_b-mixed; see
`docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md`), which is why F112 does not
already force these values to zero. G5 fails outside that family, measured: with
random per-site 2x2 jump operators the jump term stops being inert. G9 does NOT
fail there, and that is the point of its wider version below.

WHAT IT CHECKS (all gates must pass):

  G0   bond-only H (open XY chain), uniform σ⁻ damping        → asym == 0.0 exactly
  G0b  the hand-built FULL Lindbladian == the committed primitive
       `polarity_coordinates_from_hc`
  G1   bond-only H, gain/loss profile, no detuning            → asym == 0.0 exactly
  G2   bond H + 0.4·Z₀ detuning, no dissipation               → asym == 0.0 exactly
  G3   bond H + 0.4·Z₀ detuning + gain(0)/loss(N−1)           → F113's closed form
  G4   the same on −i(Hρ − ρH†) with H = A + i(γ/2)(n₀ − n_{N−1}); this is a
       CONVENTION check, not an independent configuration: that generator differs
       from G3's no-jump generator by γ·Id, and Id is Π-invariant, so it cannot
       move the asymmetry. What it verifies is that the two ways of writing a
       gain/loss system land on the same operator up to that shift.
  G5   the recycling jump c ⊗ c* is asymmetry-inert: asym(full) == asym(no-jump)
       over random Hermitian and non-Hermitian H with random per-site σ∓ rates,
       and asym(jump term alone) == 0.0 exactly
  G6   the same gain/loss profile with the field on EVERY site: the two ends carry
       +ω·γ and −ω·γ and cancel across sites
  G7   detailed balance (γ_pump,l = γ_T1,l per site) with the Z₀ detuning present
  G8   in G3's configuration the whole value sits on the gain at the DETUNED site;
       the loss at the far, undetuned site contributes nothing (the same-site law)
  G9   the drain term alone has asymmetry exactly 0.0, and not only for the σ∓
       family: also for random per-site 2x2 collapse operators and for full random
       and non-normal d×d ones. This is the title claim of
       hypotheses/THE_DRAIN_HAS_NO_CHIRALITY.md, and it is broader than that
       document scoped it. The older probe measured it on one uniform-loss
       profile only.
  G9b  what carries G9's zero is the SHAPE, not the Hermiticity of c†c: the
       transpose-paired lift −½(K ⊗ I + I ⊗ Kᵀ) sits at the noise floor for a
       NON-Hermitian K too, while the broken pairing K ⊗ I + I ⊗ conj(K) is
       macroscopic. That broken arm is this checker's positive control.
  G10  once the anti-Hermitian part B leaves the per-site n family, what sorts A
       is bit_b parity rather than F113's single-site Z, and one-directionally:
       bit_b-even A is balanced against every B tried, bit_b-odd A breaks against
       a GENERIC B but not against, say, a B proportional to the identity.
       DERIVED the same day as F155 (`simulations/f155_polarity_break_bilinear.py`,
       `docs/proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md`): this gate is
       the measured face of that law's corollary 5, and it is kept because a
       measured face of a derived law is a check on the derivation.
  GS   the sweep that keeps the rest honest: G0/G1/G2/G6/G7/G8 re-run across
       N = 2, 3, 4 × γ over six decades (1e-4 to 300) × four Z-coefficients (to 100,
       i.e. ω to 200) × two J, 144 configurations, reporting the residual against the error model below

TOLERANCES, and which zeros are bit-exact. Whether a zero comes out exactly 0.0
is an OBSERVATION here, not a taxonomy, and the sweep is what earns the words.

  * G0, G1, G2 are exactly 0.0 in all 144 swept configurations (GS-a), so they
    are compared EXACTLY and a nonzero there would be a finding about the
    construction, not a tolerance question.
  * G6, G7, G8 are exactly 0.0 in most swept configurations (122, 122 and 142 of
    144) and at the float noise level otherwise. Gating those on == 0.0 would be
    gating an input rather than the physics, so they are gated on the error model
    and the sweep prints the ratio.
  * The error model, measured rather than assumed: the noise floor is eps·‖L‖²_F,
    the squared Frobenius norm of the WHOLE Pauli-basis L. It has to be the whole
    L: under detailed balance the drain is the scalar −Nγ·Id, which sits entirely
    in M_zero, so a denominator built from ‖M₊ᵢ‖² + ‖M₋ᵢ‖² is blind to the largest
    thing the projection cancels, and the ratio then grows ∝ γ (it reached 400 at
    γ = 300 on G7, which is how this denominator was found). With ‖L‖²_F the worst
    observed ratio is 0.84, over γ from 1e-4 to 300, ω to 200, J to 2.7 and
    N = 2..4: it does not grow with any of them, and that constancy is the law.
    EPS_RATIO_BOUND = 8 is slack around it.
  * The same model gates G5, and there a fixed RELATIVE bound would be the wrong
    instrument: the relative deviation grows with N (2.8e-14, 5.4e-14, 1.1e-13 at N = 2, 3, 4,
    still inside F113's 1e-12 but on its way out of it) while
    the ratio to the noise floor stays below 1. The gate reports both.
  * The comparisons against F113's closed form (G3, G4) have a value of their own
    to divide by and sit at a few ulp (7e-16 relative). They keep the relative
    bound F113 itself is registered at, 1e-12
    (`experiments/F113_BREAK_MAGNITUDE_FORMULA.md`).

VEC CONVENTION: mirrors `framework.lindblad.lindbladian_general` exactly, as
`simulations/f112_nojump_cancellation_gate.py` does:

    L_H     = -i (H ⊗ I − I ⊗ H^T)
    jump    =    c ⊗ conj(c)
    drain   = -1/2 ( c†c ⊗ I + I ⊗ (c†c)^T )

The hand build exists only because the primitive cannot express the no-jump and
jump-only splits; G0b pins it against the primitive on the full Lindbladian.

The SIGN of the asymmetry is convention-attached (F113 pins this: the pipeline's
row-stack L read against the order='F' Pauli transform). The magnitude, the
zero-versus-nonzero verdict, and the relative sign between pump and loss are not.

NOTE ON GATE LABELS: `simulations/f112_nojump_cancellation_gate.py` also uses
G-labels, for different gates, and there a FIRED gate is the finding. Here every
gate must pass. The two label sets are unrelated.
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
from framework.diagnostics.polarity_coordinates import (
    polarity_coordinates_from_L,
    polarity_coordinates_from_hc,
)

SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)   # |0><1|, lowering
SIGMA_PLUS = np.array([[0, 0], [1, 0]], dtype=complex)    # |1><0|, raising
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)
N_UP = np.array([[0, 0], [0, 1]], dtype=complex)          # |1><1|

EPS = float(np.finfo(float).eps)
GAMMA = 0.1
DETUNING = 0.4          # coefficient of Z_0 in H, so omega_0 = 2 * DETUNING
NS = (2, 3)
REL_BOUND = 1e-12       # F113's own registered bound, for value-vs-formula only
EPS_RATIO_BOUND = 8.0   # slack over the observed worst ratio 0.84 (see TOLERANCES)

GATES = []
_PI_CACHE = {}


def gate(name, passed, detail):
    GATES.append((name, bool(passed), detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}")
    return bool(passed)


def site_op(N, l, op2):
    """N-qubit operator with a 2x2 op on site l, identity elsewhere (site 0 first)."""
    ops = [np.eye(2, dtype=complex)] * N
    ops[l] = op2
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def H_xy(N, J=1.0):
    """Uniform open XY chain, H = J * sum_l (1/2)(X_l X_{l+1} + Y_l Y_{l+1})."""
    bonds = [(i, i + 1) for i in range(N - 1)]
    return _build_bilinear(N, bonds, [("X", "X", 0.5 * J), ("Y", "Y", 0.5 * J)])


def pieces(H, c_ops):
    """(L_H, L_jump, L_drain) in vec form; c_ops already carry sqrt(rate)."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L_H = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    L_jump = np.zeros((d * d, d * d), dtype=complex)
    L_drain = np.zeros((d * d, d * d), dtype=complex)
    for c in c_ops:
        cdc = c.conj().T @ c
        L_jump += np.kron(c, c.conj())
        L_drain += -0.5 * np.kron(cdc, Id) - 0.5 * np.kron(Id, cdc.T)
    return L_H, L_jump, L_drain


def polarity(L_vec, N):
    """(asymmetry, ‖L‖²_F in the Pauli basis).

    The second value is the noise scale, and it must be the norm of the WHOLE
    Pauli-basis L, not of its ±i parts: under detailed balance the drain is the
    scalar −Nγ·Id, which lands entirely in M_zero, so a ±i-only scale is blind to
    the largest thing the projection had to cancel and the ratio then grows with γ.
    """
    if N not in _PI_CACHE:
        _PI_CACHE[N] = (_vec_to_pauli_basis_transform(N), build_pi_full(N))
    T, Pi = _PI_CACHE[N]
    L_pauli = (T.conj().T @ L_vec @ T) / (2 ** N)
    res = polarity_coordinates_from_L(L_pauli, N, 0.0, Pi=Pi)
    return float(res["asymmetry"]), float(np.sum(np.abs(L_pauli) ** 2))


def asymmetry(L_vec, N):
    return polarity(L_vec, N)[0]


def nonhermitian_generator(H):
    """-i (H rho - rho H^dagger) in vec form, for a genuinely non-Hermitian H."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.conj()))


def rel_dev(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale != 0.0 else abs(a - b)


def eps_ratio(residual, scale):
    """|residual| in units of the measured noise floor eps * ‖L‖²_F."""
    return abs(residual) / (EPS * scale) if scale > 0.0 else 0.0


def no_jump(H, c_ops, N):
    L_H, _, L_drain = pieces(H, c_ops)
    return polarity(L_H + L_drain, N)


def full(H, c_ops, N):
    L_H, L_jump, L_drain = pieces(H, c_ops)
    return polarity(L_H + L_jump + L_drain, N)


def profiles(N, gamma, det_coeff, J=1.0):
    """The seven configurations the gates and the sweep share."""
    H_bond = H_xy(N, J)
    H_det = H_bond + det_coeff * site_op(N, 0, PAULI_Z)
    H_field = H_bond + sum(det_coeff * site_op(N, l, PAULI_Z) for l in range(N))
    loss_all = [np.sqrt(gamma) * site_op(N, l, SIGMA_MINUS) for l in range(N)]
    gain_loss = [np.sqrt(gamma) * site_op(N, 0, SIGMA_PLUS),
                 np.sqrt(gamma) * site_op(N, N - 1, SIGMA_MINUS)]
    balanced = loss_all + [np.sqrt(gamma) * site_op(N, l, SIGMA_PLUS) for l in range(N)]
    return {
        "G0": (H_bond, loss_all),            # no carrier
        "G1": (H_bond, gain_loss),           # no carrier, gain/loss profile
        "G2": (H_det, []),                   # carrier, no dissipator
        "G6": (H_field, gain_loss),          # carrier everywhere, ends cancel
        "G7": (H_det, balanced),             # carrier, detailed balance
        "G8": (H_det, [gain_loss[1]]),       # carrier and damping on different sites
        "G3": (H_det, gain_loss),            # the one that fires
    }


# Which zeros come out bit-exact is an OBSERVATION over the sweep, not a
# taxonomy: G0/G1/G2 are exact in every configuration swept, G6/G7/G8 in most.
EXACT_ZEROS = ("G0", "G1", "G2")
NOISE_ZEROS = ("G6", "G7", "G8")

print("=" * 78)
print("F112 GAIN/LOSS CARRIER CHECK  (the label decides nothing, the sum does)")
print(f"gamma = {GAMMA}, detuning = {DETUNING}*Z_0 (omega_0 = {2 * DETUNING}), "
      f"J = 1, N = {NS}")
print("=" * 78)

for N in NS:
    print(f"\n{'-' * 78}\nN = {N}\n{'-' * 78}")
    cfg = profiles(N, GAMMA, DETUNING)
    H_det, gain_loss = cfg["G3"]

    a_hand, _ = full(H_det, gain_loss, N)
    res = polarity_coordinates_from_hc(
        H_det, [site_op(N, 0, SIGMA_PLUS), site_op(N, N - 1, SIGMA_MINUS)],
        [GAMMA, GAMMA], N, sigma=0.0, Pi=build_pi_full(N))
    a_prim = float(res["asymmetry"])
    gate(f"G0b N={N} hand-built FULL == polarity_coordinates_from_hc",
         rel_dev(a_hand, a_prim) <= REL_BOUND,
         f"hand {a_hand:+.12f}, primitive {a_prim:+.12f}, "
         f"rel dev {rel_dev(a_hand, a_prim):.2e}")

    for key, why in (("G0", "no Z-carrier, uniform loss"),
                     ("G1", "no Z-carrier, gain/loss profile"),
                     ("G2", "Z-carrier, no dissipator")):
        a, _ = no_jump(*cfg[key], N)
        gate(f"{key} N={N} {why}", a == 0.0,
             f"asym = {a!r} (exact zero required: exact in every swept config)")

    predicted = (4 ** N / 2) * (2 * DETUNING) * (GAMMA - 0.0)
    a_nojump, _ = no_jump(H_det, gain_loss, N)
    a_full, _ = full(H_det, gain_loss, N)
    dev = max(rel_dev(a_nojump, predicted), rel_dev(a_full, predicted))
    gate(f"G3 N={N} Z_0 detuning + gain/loss vs F113",
         dev <= REL_BOUND,
         f"no-jump {a_nojump:+.12f}, full {a_full:+.12f}, F113 {predicted:+.12f}, "
         f"rel dev {dev:.2e}")

    B = (GAMMA / 2.0) * (site_op(N, 0, N_UP) - site_op(N, N - 1, N_UP))
    L_nh = nonhermitian_generator(H_det + 1j * B)
    L_H, _, L_drain = pieces(H_det, gain_loss)
    shift = L_nh - (L_H + L_drain)
    scalar_resid = float(np.max(np.abs(shift - GAMMA * np.eye(4 ** N))))
    a_nh = asymmetry(L_nh, N)
    gate(f"G4 N={N} -i(H rho - rho H^dag) is that generator up to gamma*Id",
         rel_dev(a_nh, predicted) <= REL_BOUND and scalar_resid <= 1e-14,
         f"asym {a_nh:+.12f} (F113 {predicted:+.12f}); its difference to the "
         f"no-jump generator is gamma*Id to {scalar_resid:.1e} (SCALAR, not merely "
         f"diagonal, which is what makes it Pi-invariant), so this is a convention "
         f"check, not a second configuration")

    for key, why in (("G6", "field on EVERY site: the two ends cancel"),
                     ("G7", "detailed balance with the Z-carrier present"),
                     ("G8", "carrier and damping on DIFFERENT sites (same-site law)")):
        a, scale = no_jump(*cfg[key], N)
        r = eps_ratio(a, scale)
        gate(f"{key} N={N} {why}", r <= EPS_RATIO_BOUND,
             f"asym = {a!r}, |asym|/(eps*||L||_F^2) = {r:.2f} "
             f"(zero up to the noise floor, exact in most swept configs)")

# ---------------------------------------------------------------------------
# G5 / G9: random H, random per-site rates. The jump term and the drain term.
# ---------------------------------------------------------------------------
rng = np.random.default_rng(20260819)
worst_ratio = 0.0
worst_rel_by_N = {}
worst_jump = 0.0
worst_drain = 0.0
for N in (2, 3, 4):
    d = 2 ** N
    for _ in range(3):
        A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        for H in ((A + A.conj().T) / 2, A):
            rates_lo = rng.uniform(0.0, 0.3, size=N)
            rates_hi = rng.uniform(0.0, 0.3, size=N)
            c_ops = ([np.sqrt(r) * site_op(N, l, SIGMA_MINUS) for l, r in enumerate(rates_lo)]
                     + [np.sqrt(r) * site_op(N, l, SIGMA_PLUS) for l, r in enumerate(rates_hi)])
            L_H, L_jump, L_drain = pieces(H, c_ops)
            a_full, s_full = polarity(L_H + L_jump + L_drain, N)
            a_nojump, s_nojump = polarity(L_H + L_drain, N)
            worst_ratio = max(worst_ratio,
                              eps_ratio(a_full - a_nojump, max(s_full, s_nojump)))
            worst_rel_by_N[N] = max(worst_rel_by_N.get(N, 0.0),
                                    rel_dev(a_full, a_nojump))
            worst_jump = max(worst_jump, abs(asymmetry(L_jump, N)))
            worst_drain = max(worst_drain, abs(asymmetry(L_drain, N)))

print()
rel_by_N = ", ".join(f"N={n}: {v:.2e}" for n, v in sorted(worst_rel_by_N.items()))
gate("G5 recycling jump is asymmetry-inert (18 random H, N = 2, 3, 4, mixed rates)",
     worst_ratio <= EPS_RATIO_BOUND and worst_jump == 0.0,
     f"worst |full - no-jump|/(eps*||L||_F^2) = {worst_ratio:.2f}; "
     f"worst RELATIVE dev by N ({rel_by_N}) grows with N while the ratio does not, "
     f"which is why the relative bound is the wrong instrument here; "
     f"worst |asym(jump alone)| = {worst_jump!r}")

# G9: the drain's chirality-freedom is not about the sigma-minus/sigma-plus
# family. drain = -1/2 (K (x) I + I (x) K^T) with K = c^dag c HERMITIAN for ANY
# c, so the sweep runs three kinds of collapse operator, not one.
worst_drain_wide = worst_drain
kinds = []
for N in (2, 3):
    d = 2 ** N
    for kind in ("per-site 2x2", "full dxd", "non-normal full"):
        for _ in range(3):
            if kind == "per-site 2x2":
                c_ops = [site_op(N, l, rng.normal(size=(2, 2))
                                 + 1j * rng.normal(size=(2, 2))) for l in range(N)]
            else:
                m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
                c_ops = [m if kind == "full dxd" else np.triu(m)]
            _, _, L_drain = pieces(np.zeros((d, d), dtype=complex), c_ops)
            worst_drain_wide = max(worst_drain_wide, abs(asymmetry(L_drain, N)))
        kinds.append(f"{kind} N={N}")

gate("G9 the DRAIN alone has no chirality, for ANY collapse operator",
     worst_drain_wide == 0.0,
     f"worst |asym(drain alone)| = {worst_drain_wide!r} over the sigma-minus/plus "
     f"profiles above plus {len(kinds)} random families ({', '.join(kinds[:3])}, ...); "
     f"exact zero required; it is broader than the sigma-minus/plus family "
     f"because the drain's SHAPE, -1/2(K (x) I + I (x) K^T), has zero asymmetry "
     f"for every K, and K = c^dag c whatever c is")

# ---------------------------------------------------------------------------
# G9b: which SHAPE carries G9's zero. The drain is the transpose-paired lift of
# one operator to both sides; that shape, not the Hermiticity of c^dag c, is what
# kills the asymmetry. The broken pairing is the positive control: a checker that
# only ever reports zeros has not shown that it can see anything.
# ---------------------------------------------------------------------------
paired_worst = 0.0
broken_smallest = float("inf")
for N in (2, 3, 4):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    for _ in range(3):
        K = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))   # non-Hermitian
        a_paired, s_paired = polarity(-0.5 * (np.kron(K, Id) + np.kron(Id, K.T)), N)
        paired_worst = max(paired_worst, eps_ratio(a_paired, s_paired))
        a_broken, s_broken = polarity(np.kron(K, Id) + np.kron(Id, K.conj()), N)
        broken_smallest = min(broken_smallest, eps_ratio(a_broken, s_broken))

gate("G9b the transpose pairing is what carries it, not Hermiticity",
     paired_worst <= EPS_RATIO_BOUND and broken_smallest > 1e6,
     f"-1/2(K (x) I + I (x) K^T) with NON-Hermitian K: worst ratio {paired_worst:.2f} "
     f"(at the noise floor, so the shape alone suffices); break the pairing to "
     f"K (x) I + I (x) conj(K) and the smallest ratio is {broken_smallest:.1e}, "
     f"macroscopic (the positive control: this checker can see a nonzero)")

# ---------------------------------------------------------------------------
# G10: what sorts the general case, once the anti-Hermitian part B leaves the
# per-site n family. Not F113's single-site Z: bit_b parity of A. The sorting is
# ONE-directional, and the gate says so: bit_b-even A is balanced against every B
# tried, while bit_b-odd A breaks for a GENERIC B and can be balanced for special
# ones (B proportional to the identity, say).
# ---------------------------------------------------------------------------
even_worst = 0.0
odd_breaks = 0
odd_total = 0
for N in (2, 3):
    d = 2 ** N
    for _ in range(2):
        m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        B = (m + m.conj().T) / 2
        for k in range(1, 4 ** N):
            A = pauli_string(_k_to_indices(k, N))
            a, scale = polarity(nonhermitian_generator(A + 1j * B), N)
            r = eps_ratio(a, scale)
            if sum(1 for c in _pauli_label(k, N) if c in "YZ") % 2 == 0:
                even_worst = max(even_worst, r)
            else:
                odd_total += 1
                odd_breaks += (r > 1e6)

gate("G10 bit_b parity of A sorts the general-B case (one-directionally)",
     even_worst <= EPS_RATIO_BOUND and odd_breaks == odd_total,
     f"bit_b-EVEN A: worst ratio {even_worst:.2f} (balanced, exactly 0.0 in most "
     f"draws and at the noise floor otherwise); bit_b-ODD A: {odd_breaks}/{odd_total} "
     f"break against a random Hermitian B, and the converse does NOT hold, a special "
     f"B (proportional to the identity) leaves odd A balanced too")

# ---------------------------------------------------------------------------
# GS: the sweep. Does the exactness survive a change of gamma, omega, J, N?
# ---------------------------------------------------------------------------
sweep_worst = {k: 0.0 for k in EXACT_ZEROS + NOISE_ZEROS}
sweep_exact = {k: 0 for k in sweep_worst}
sweep_total = 0
for N in (2, 3, 4):
    for gamma, det_coeff, J in itertools.product(
            (1e-4, 0.01, 0.1, 3.0, 50.0, 300.0),
            (0.4, 0.1234567, float(np.sqrt(2)), 100.0), (1.0, 2.7)):
        cfg = profiles(N, gamma, det_coeff, J)
        sweep_total += 1
        for key in sweep_worst:
            a, scale = no_jump(*cfg[key], N)
            if a == 0.0:
                sweep_exact[key] += 1
            sweep_worst[key] = max(sweep_worst[key], eps_ratio(a, scale))

print()
gate(f"GS-a G0/G1/G2 stay EXACT across the sweep ({sweep_total} configs, N = 2..4, "
     f"gamma 1e-4..300, omega to 200, two J)",
     all(sweep_exact[k] == sweep_total for k in EXACT_ZEROS),
     ", ".join(f"{k}: {sweep_exact[k]}/{sweep_total} exact" for k in EXACT_ZEROS))

gate("GS-b G6/G7/G8 stay at the noise floor across the same sweep",
     all(sweep_worst[k] <= EPS_RATIO_BOUND for k in NOISE_ZEROS),
     ", ".join(f"{k}: worst ratio {sweep_worst[k]:.2f}, exact in "
               f"{sweep_exact[k]}/{sweep_total}" for k in NOISE_ZEROS)
     + f" (bound {EPS_RATIO_BOUND:.0f}; what stays put as the scale moves is the "
       f"ratio, not the residual)")

print("\n" + "=" * 78)
failed = [n for n, ok, _ in GATES if not ok]
print(f"{len(GATES) - len(failed)}/{len(GATES)} gates passed")
if failed:
    print("FIRED: " + ", ".join(failed))
print("=" * 78)
assert not failed, f"gates fired: {failed}"
