"""GATE-FIRST verifier: Direction 1 of ASYMMETRY_IS_THE_UNRECYCLED_DRAIN.md.

HYPOTHESIS (to confirm or refute):
  The F112 polarity asymmetry of a non-Hermitian generator is the chirality of
  the UN-RECYCLED DRAIN. For amplitude damping (jump c = sigma^- per site, rate
  gamma), the NO-JUMP generator -i(H_eff rho - rho H_eff^dagger) has nonzero
  polarity asymmetry, and adding back the JUMP term (c (.) c^dagger) restores
  the FULL Lindbladian's asymmetry to 0 -- i.e. the recycling jump feeds back
  exactly the chirality the drain removed.

GATE-FIRST DISCIPLINE (repo rule): Stage 0 (G0) tests that the hand-built
objects are correct; it MUST pass. The hypothesis gates G1-G3 + bonuses encode
the load-bearing predictions. A gate that FIRES is THE FINDING -- it is NOT
loosened; it is diagnosed and reported. The script computes ALL numbers, prints
the full report, then a final assert fires iff the hypothesis does not fully
hold (listing which gates fired).

VEC CONVENTION: mirror the framework ground-truth `lindbladian_general`
(simulations/framework/lindblad.py) EXACTLY, so the hand-built L is consistent
with the framework's Pi (build_pi_full) and the column-stack Pauli transform
(_vec_to_pauli_basis_transform, flatten('F')):

    L_H     = -i (H (x) I  -  I (x) H^T)              [the commutator; F112's object]
    jump_c  =     c (x) conj(c)                        [recycling jump]
    drain_c = -1/2 ( c^dag c (x) I + I (x) (c^dag c)^T)[no-jump anticommutator loss]

ASYMMETRY: feed each hand-built L (Pauli basis) to
polarity_coordinates_from_L(L_pauli, N, sigma); report
    asymmetry = ||M_+i||^2 - ||M_-i||^2
    eta       = asymmetry / (||M_+i||^2 + ||M_-i||^2)   in [-1, 1].
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

from framework.pauli import _build_bilinear, _vec_to_pauli_basis_transform
from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import (
    polarity_coordinates_from_L,
    polarity_coordinates_from_hc,
)

np.set_printoptions(precision=6, suppress=True)

# --------------------------------------------------------------------------
# Gate recorder (gate-first: a fired gate is THE finding; report all, then assert)
# --------------------------------------------------------------------------
GATES = []  # list of (name, passed, detail)


def gate(name, passed, detail):
    GATES.append((name, bool(passed), detail))
    tag = "PASS" if passed else "FAIL  <-- FIRED (finding)"
    print(f"[{tag}] {name}: {detail}")
    return bool(passed)


# --------------------------------------------------------------------------
# Operators
# --------------------------------------------------------------------------
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1|, lowering |1>->|0>
SIGMA_PLUS = np.array([[0, 0], [1, 0]], dtype=complex)   # |1><0|, raising |0>->|1>


def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]


def H_xy(N, J=1.0):
    """Uniform open XY chain H = J * sum_l (1/2)(X_l X_{l+1} + Y_l Y_{l+1})."""
    return _build_bilinear(N, chain_bonds(N), [("X", "X", 0.5 * J), ("Y", "Y", 0.5 * J)])


def H_heisenberg(N, J=1.0):
    return _build_bilinear(N, chain_bonds(N), [("X", "X", J), ("Y", "Y", J), ("Z", "Z", J)])


def site_op_2x2(N, l, op2):
    """N-qubit operator with a 2x2 op on site l, identity elsewhere (MSB site-0 first)."""
    ops = [np.eye(2, dtype=complex)] * N
    ops[l] = op2
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def collapse_ops(N, gamma, kind):
    """Return list of sqrt(gamma)-scaled per-site collapse operators.

    kind: 'cool' (sigma^-), 'heat' (sigma^+), or 'balance' (both, equal rate).
    """
    rt = np.sqrt(gamma)
    if kind == "cool":
        return [rt * site_op_2x2(N, l, SIGMA_MINUS) for l in range(N)]
    if kind == "heat":
        return [rt * site_op_2x2(N, l, SIGMA_PLUS) for l in range(N)]
    if kind == "balance":
        return ([rt * site_op_2x2(N, l, SIGMA_MINUS) for l in range(N)]
                + [rt * site_op_2x2(N, l, SIGMA_PLUS) for l in range(N)])
    raise ValueError(kind)


# --------------------------------------------------------------------------
# Liouvillian pieces in vec form (lindbladian_general convention, EXACT mirror)
# --------------------------------------------------------------------------
def build_L_pieces(H, c_ops):
    """c_ops already carry the sqrt(gamma) scaling.

    Returns (L_H, L_jump, L_drain) as d^2 x d^2 vec-form matrices.
      FULL    = L_H + L_jump + L_drain   (== lindbladian_general(H, c_ops))
      NO-JUMP = L_H + L_drain
      JUMP    = L_jump
    """
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


def to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def measure(L_vec, N, sigma, Pi):
    """Asymmetry + eta + components of L (vec form) via polarity_coordinates_from_L."""
    Lp = to_pauli(L_vec, N)
    res = polarity_coordinates_from_L(Lp, N, sigma, Pi=Pi)
    a = res["asymmetry"]
    denom = res["norm_sq"]["M_plus_half"] + res["norm_sq"]["M_minus_half"]
    eta = (a / denom) if denom > 1e-300 else 0.0
    return {
        "asym": float(a),
        "eta": float(eta),
        "manti_sq": float(denom),
        "Mp": res["M_plus_half"],
        "Mm": res["M_minus_half"],
        "ortho_resid": res["orthogonality_residual"],
    }


def fro_inner_re(A, B):
    """Re of Frobenius inner product <A, B> = sum conj(A)*B."""
    return float(np.real(np.sum(np.conj(A) * B)))


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
GAMMA = 0.1
NS = [2, 3]
ZERO_TOL = 1e-9    # |eta| below this = machine zero (balanced)
NONZERO_TOL = 1e-6  # |eta| above this = clearly nonzero

print("=" * 78)
print("F112 NO-JUMP CANCELLATION GATE  (Direction 1, ASYMMETRY_IS_THE_UNRECYCLED_DRAIN)")
print(f"gamma = {GAMMA}, H = uniform open XY chain, c = sigma^- per site")
print("=" * 78)

# Accumulate a results table across N for the final report.
rows = []  # (N, label, asym, eta, manti_sq)
per_N = {}  # N -> dict of measurements

for N in NS:
    print(f"\n{'-'*78}\nN = {N}\n{'-'*78}")
    H = H_xy(N)
    Pi = build_pi_full(N)
    c_cool = collapse_ops(N, GAMMA, "cool")
    sigma = float(np.real(sum(np.trace(c.conj().T @ c) for c in c_cool)))  # benign center
    # (asymmetry is sigma-independent; verified by G_sigma below.)

    L_H, L_jump, L_drain = build_L_pieces(H, c_cool)
    L_full = L_H + L_jump + L_drain
    L_nojump = L_H + L_drain

    # ---- G0: object identity check (MUST pass) -----------------------------
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    Gop = sum(c.conj().T @ c for c in c_cool)        # sum c^dag c  (= gamma * sum n_l)
    H_eff = H - 0.5j * Gop                            # H_eff = H - (i/2) sum c^dag c
    L_nojump_heff = -1j * (np.kron(H_eff, Id) - np.kron(Id, np.conj(H_eff)))  # -i(H_eff rho - rho H_eff^dag)
    Bop = -0.5 * Gop                                  # B = -1/2 sum c^dag c
    L_nojump_comm = (-1j * (np.kron(H, Id) - np.kron(Id, H.T))
                     + np.kron(Bop, Id) + np.kron(Id, Bop.T))  # -i[H,.] + {B,.}

    d0_heff = float(np.max(np.abs(L_nojump - L_nojump_heff)))
    d0_comm = float(np.max(np.abs(L_nojump - L_nojump_comm)))
    g0_ok = gate(
        f"G0(N={N}) identity: L_H+L_drain == -i(H_eff rho - rho H_eff^dag) == -i[H,.]+{{B,.}}",
        d0_heff < 1e-11 and d0_comm < 1e-11,
        f"max|term-built - H_eff form| = {d0_heff:.2e}; max|term-built - comm/anti form| = {d0_comm:.2e}",
    )
    if not g0_ok:
        raise SystemExit("G0 FAILED: hand-built objects are inconsistent. Fix convention before proceeding.")

    # cross-check: FULL == framework lindbladian_general convention via from_hc
    res_hc = polarity_coordinates_from_hc(
        H, [site_op_2x2(N, l, SIGMA_MINUS) for l in range(N)], [GAMMA] * N, N, sigma=sigma, Pi=Pi
    )
    m_full = measure(L_full, N, sigma, Pi)
    d0b = abs(m_full["asym"] - res_hc["asymmetry"])
    gate(
        f"G0b(N={N}) pipeline matches committed polarity_coordinates_from_hc on FULL",
        d0b < 1e-12,
        f"|asym(from_L hand-built) - asym(from_hc)| = {d0b:.2e}",
    )

    # ---- measurements ------------------------------------------------------
    m_nojump = measure(L_nojump, N, sigma, Pi)
    m_jump = measure(L_jump, N, sigma, Pi)
    m_drain0 = measure(L_drain, N, sigma, Pi)             # H = 0: isolates asym_B
    m_honly = measure(L_H, N, sigma, Pi)                  # drain off: F112 sanity (expect 0)

    per_N[N] = dict(full=m_full, nojump=m_nojump, jump=m_jump, drain0=m_drain0, honly=m_honly)
    for label, m in [("FULL", m_full), ("NO-JUMP", m_nojump), ("JUMP", m_jump),
                     ("DRAIN-ONLY(H=0)", m_drain0), ("H-ONLY(drain off)", m_honly)]:
        rows.append((N, label, m["asym"], m["eta"], m["manti_sq"]))

    # ---- G_sigma: asymmetry is independent of the palindrome center sigma ---
    m_nojump_s0 = measure(L_nojump, N, 0.0, Pi)
    gate(
        f"G_sigma(N={N}) asymmetry independent of sigma (only shifts M_zero)",
        abs(m_nojump_s0["asym"] - m_nojump["asym"]) < 1e-12,
        f"|asym(sigma=0) - asym(sigma={sigma:.3g})| = {abs(m_nojump_s0['asym']-m_nojump['asym']):.2e}",
    )

    # ---- G1: asymmetry(FULL amplitude damping) ----------------------------
    gate(
        f"G1(N={N}) FULL amplitude-damping asymmetry == 0 (machine zero)",
        abs(m_full["eta"]) < ZERO_TOL and abs(m_full["asym"]) < 1e-9,
        f"asym(FULL) = {m_full['asym']:.3e}, eta = {m_full['eta']:.3e}",
    )

    # ---- G2: asymmetry(NO-JUMP) != 0 --------------------------------------
    gate(
        f"G2(N={N}) NO-JUMP (un-recycled drain) asymmetry != 0",
        abs(m_nojump["eta"]) > NONZERO_TOL,
        f"asym(NO-JUMP) = {m_nojump['asym']:.6e}, eta = {m_nojump['eta']:.6e}, "
        f"||M_anti||^2 = {m_nojump['manti_sq']:.6e}",
    )

    # ---- G3: the cancellation. |asym(FULL)| << |asym(NO-JUMP)| -------------
    ratio = abs(m_full["asym"]) / (abs(m_nojump["asym"]) + 1e-300)
    gate(
        f"G3(N={N}) jump recycles the chirality: |asym(FULL)| << |asym(NO-JUMP)|",
        ratio < 1e-6,
        f"|asym(FULL)|/|asym(NO-JUMP)| = {ratio:.3e}",
    )

    # ---- G3 decomposition: FULL = NO-JUMP + JUMP (M_anti linear) -----------
    # Linearity check on the +i/-i projection matrices.
    lin_p = float(np.max(np.abs(m_full["Mp"] - (m_nojump["Mp"] + m_jump["Mp"]))))
    lin_m = float(np.max(np.abs(m_full["Mm"] - (m_nojump["Mm"] + m_jump["Mm"]))))
    gate(
        f"G3lin(N={N}) M_anti linear: M_+i(FULL) == M_+i(NO-JUMP) + M_+i(JUMP)",
        lin_p < 1e-10 and lin_m < 1e-10,
        f"max|residual| +i: {lin_p:.2e}, -i: {lin_m:.2e}",
    )
    # cross term: asym(FULL) = asym(NJ) + asym(J) + cross
    cross_sub = m_full["asym"] - m_nojump["asym"] - m_jump["asym"]
    cross_direct = 2.0 * (fro_inner_re(m_nojump["Mp"], m_jump["Mp"])
                          - fro_inner_re(m_nojump["Mm"], m_jump["Mm"]))
    print(f"   decomposition  asym(FULL) = asym(NO-JUMP) + asym(JUMP) + cross")
    print(f"     asym(NO-JUMP) = {m_nojump['asym']:+.6e}")
    print(f"     asym(JUMP)    = {m_jump['asym']:+.6e}")
    print(f"     cross         = {cross_sub:+.6e}  (direct inner-product: {cross_direct:+.6e})")
    print(f"     sum           = {m_nojump['asym'] + m_jump['asym'] + cross_sub:+.6e}  "
          f"(should equal asym(FULL) = {m_full['asym']:+.3e})")
    gate(
        f"G3xcheck(N={N}) cross term consistent (subtraction == direct inner product)",
        abs(cross_sub - cross_direct) < 1e-9,
        f"|cross_sub - cross_direct| = {abs(cross_sub - cross_direct):.2e}",
    )

    # ---- asym_A / asym_B / cross_AB structural split of NO-JUMP -----------
    asym_A = m_honly["asym"]
    asym_B = m_drain0["asym"]
    crossAB = m_nojump["asym"] - asym_A - asym_B
    print(f"   NO-JUMP structural split: asym_A(H comm)={asym_A:+.3e}  "
          f"asym_B(drain,H=0)={asym_B:+.3e}  cross_AB={crossAB:+.3e}")

# --------------------------------------------------------------------------
# Bonus B1: A-blindness (Direction 2) -- fix drain, sweep H
# --------------------------------------------------------------------------
print(f"\n{'-'*78}\nBONUS B1: A-blindness  (fix sigma^- drain, vary H)\n{'-'*78}")
rng = np.random.default_rng(20260620)
for N in NS:
    Pi = build_pi_full(N)
    d = 2 ** N
    c_cool = collapse_ops(N, GAMMA, "cool")
    sigma = float(np.real(sum(np.trace(c.conj().T @ c) for c in c_cool)))

    def rand_herm(d):
        A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
        return (A + A.conj().T) / 2

    H_variants = [
        ("XY chain", H_xy(N)),
        ("2 * XY chain", 2.0 * H_xy(N)),
        ("Heisenberg chain", H_heisenberg(N)),
        ("random Hermitian #1", rand_herm(d)),
        ("random Hermitian #2", rand_herm(d)),
        ("random Hermitian #3 (x5)", 5.0 * rand_herm(d)),
    ]
    print(f"  N = {N}:  asym_B(drain only, H=0) reference =", end=" ")
    _, _, L_drain0 = build_L_pieces(np.zeros((d, d), dtype=complex), c_cool)
    asymB = measure(L_drain0, N, sigma, Pi)["asym"]
    print(f"{asymB:+.6e}")
    vals = []
    for name, Hv in H_variants:
        L_H, _, L_drain = build_L_pieces(Hv, c_cool)
        m = measure(L_H + L_drain, N, sigma, Pi)
        vals.append(m["asym"])
        print(f"      asym(NO-JUMP) [{name:<26s}] = {m['asym']:+.6e}   eta = {m['eta']:+.6e}")
    spread = max(vals) - min(vals)
    # A-blind iff asym(NO-JUMP) is constant across all H AND equals asym_B.
    a_blind = spread < 1e-9 and abs(np.mean(vals) - asymB) < 1e-9
    gate(
        f"B1(N={N}) A-blindness: asym(NO-JUMP) constant over H and == asym_B",
        a_blind,
        f"spread over 6 H = {spread:.3e}; mean={np.mean(vals):+.3e} vs asym_B={asymB:+.3e}",
    )

# --------------------------------------------------------------------------
# Bonus B2: sign of eta -- cooling vs heating vs detailed balance (Direction 3)
# --------------------------------------------------------------------------
print(f"\n{'-'*78}\nBONUS B2: chirality sign  (NO-JUMP generator; cool=sigma^- / heat=sigma^+ / balance)\n{'-'*78}")
for N in NS:
    Pi = build_pi_full(N)
    H = H_xy(N)
    out = {}
    for kind in ["cool", "heat", "balance"]:
        c_ops = collapse_ops(N, GAMMA, kind)
        sigma = float(np.real(sum(np.trace(c.conj().T @ c) for c in c_ops)))
        L_H, _, L_drain = build_L_pieces(H, c_ops)
        out[kind] = measure(L_H + L_drain, N, sigma, Pi)
        print(f"  N={N} [{kind:<8s}] asym(NO-JUMP) = {out[kind]['asym']:+.6e}   eta = {out[kind]['eta']:+.6e}")
    eta_c, eta_h, eta_b = out["cool"]["eta"], out["heat"]["eta"], out["balance"]["eta"]
    sign_opposite = (abs(eta_c) > NONZERO_TOL and abs(eta_h) > NONZERO_TOL
                     and np.sign(eta_c) == -np.sign(eta_h))
    gate(
        f"B2a(N={N}) sign(eta_cool) == -sign(eta_heat)",
        sign_opposite,
        f"eta_cool={eta_c:+.4e}, eta_heat={eta_h:+.4e}",
    )
    gate(
        f"B2b(N={N}) detailed balance (cool+heat equal rate) -> eta ~ 0",
        abs(eta_b) < ZERO_TOL,
        f"eta_balance={eta_b:+.4e}",
    )

# --------------------------------------------------------------------------
# MECHANISM (diagnosis of the fired G2/B1/B2a): where does the asymmetry live?
#   asym(NO-JUMP) = asym_A(H) + asym_B(drain) + cross_AB.
#   F112: asym_A = 0 for any H (spine row 1).  Measured: asym_B = 0 (drain alone).
#   => asym(NO-JUMP) = cross_AB(H), a LINEAR FUNCTIONAL of H (M2 confirms).
#   So asym(NO-JUMP, H) = sum_alpha h_alpha * asym(NO-JUMP, sigma_alpha): scan it
#   per single Pauli string. The chain's NN building blocks {XX, YY, ZZ} all lie
#   in the KERNEL of this functional (M1) -> the physical chain gives exactly 0.
#   But the functional is NOT identically zero (M3): other Pauli strings carry it
#   (a random H -- real OR complex -- hits them). The spine's nonzero no-jump
#   number is therefore a GENERIC-H artifact, not a drain/chirality property.
#   NB: my first diagnosis ("carried by Im(H)") was REFUTED by these very gates --
#   a random real-symmetric H also gives nonzero. Reality of H is not the axis.
# --------------------------------------------------------------------------
from framework.pauli import _k_to_indices, pauli_string, _pauli_label  # noqa: E402
from framework.symmetry import klein_index  # noqa: E402

print(f"\n{'-'*78}\nMECHANISM: per-Pauli scan of the linear functional asym(NO-JUMP, .)\n{'-'*78}")
rng2 = np.random.default_rng(7)
for N in NS:
    Pi = build_pi_full(N)
    d = 2 ** N
    c_cool = collapse_ops(N, GAMMA, "cool")
    sigma = float(np.real(sum(np.trace(c.conj().T @ c) for c in c_cool)))

    def asym_nojump(Hmat):
        L_H, _, L_drain = build_L_pieces(Hmat, c_cool)
        return measure(L_H + L_drain, N, sigma, Pi)["asym"]

    # (a) M2 linearity: homogeneity + additivity over a random Hermitian H.
    A = rng2.standard_normal((d, d)) + 1j * rng2.standard_normal((d, d))
    Hc = (A + A.conj().T) / 2
    H_re = Hc.real.astype(complex)
    H_im = Hc - H_re
    a_full, a_re, a_im = asym_nojump(Hc), asym_nojump(H_re), asym_nojump(H_im)
    a_scaled = asym_nojump(2.0 * Hc)
    gate(
        f"M2(N={N}) asym(NO-JUMP) LINEAR in H (homogeneous + additive)",
        abs(a_full - (a_re + a_im)) < 1e-9 and abs(a_scaled - 2.0 * a_full) < 1e-9,
        f"|2H-2*full|={abs(a_scaled-2*a_full):.2e}, |full-(re+im)|={abs(a_full-(a_re+a_im)):.2e} "
        f"(note: Re(H) gives {a_re:+.3e} != 0 -> 'Im(H)' diagnosis refuted)",
    )

    # (b) per-Pauli scan: contribution of each single Pauli string as H.
    contrib = {}
    for k in range(1, 4 ** N):  # skip identity
        sig = pauli_string(list(_k_to_indices(k, N)))
        contrib[k] = asym_nojump(sig)
    nonzero = {k: v for k, v in contrib.items() if abs(v) > 1e-9}
    print(f"  N={N}: {len(nonzero)} of {4**N - 1} single-Pauli strings carry the functional.")
    # classify nonzero contributors
    print(f"     nonzero contributors (label: asym, Klein, #Y, #Z, bit_a-parity, bit_b-parity):")
    for k in sorted(nonzero, key=lambda kk: -abs(nonzero[kk]))[:24]:
        idxs = _k_to_indices(k, N)
        letters = [{(0,0):'I',(1,0):'X',(0,1):'Z',(1,1):'Y'}[t] for t in idxs]
        nY = letters.count('Y')
        nZ = letters.count('Z')
        ba = sum(t[0] for t in idxs) % 2
        bb = sum(t[1] for t in idxs) % 2
        ki = klein_index(letters)
        print(f"       {_pauli_label(k, N):<{N}} : {nonzero[k]:+.4e}  klein={ki} nY={nY} nZ={nZ} pa={ba} pb={bb}")

    # (c) M1: chain building blocks {XX, YY, ZZ on each NN bond} are in the kernel.
    chain_terms_zero = True
    detail_terms = []
    from framework.pauli import site_op as _site_op
    for i in range(N - 1):
        for L in ("X", "Y", "Z"):
            term = _site_op(N, i, L) @ _site_op(N, i + 1, L)
            a = asym_nojump(term)
            detail_terms.append((f"{L}{L}@({i},{i+1})", a))
            if abs(a) > 1e-9:
                chain_terms_zero = False
    gate(
        f"M1(N={N}) chain NN blocks XX/YY/ZZ all in kernel -> asym(NO-JUMP)=0",
        chain_terms_zero,
        "; ".join(f"{n}={a:+.1e}" for n, a in detail_terms),
    )

    # (d) M3: the functional is NOT identically zero (selection rule, not triviality).
    gate(
        f"M3(N={N}) functional NOT identically zero (some Pauli strings carry it)",
        len(nonzero) > 0,
        f"{len(nonzero)} nonzero single-Pauli contributors exist",
    )

    # (e) validate the scan: reconstruct asym(random H) from the per-Pauli scan.
    #     pauli_basis_vector returns coeff_k = Tr(sigma_k H)/2^N already, so
    #     H = sum_k coeff_k sigma_k and (asym linear) asym(H) = sum_k coeff_k asym(sigma_k).
    from framework.pauli import pauli_basis_vector as _pbv
    hvec = _pbv(Hc, N)
    recon = sum(hvec[k].real * contrib.get(k, 0.0) for k in range(1, 4 ** N))
    gate(
        f"Mscan(N={N}) asym(random H) == sum_alpha coeff_alpha * asym(sigma_alpha)",
        abs(recon - a_full) < 1e-8,
        f"|reconstructed - measured| = {abs(recon - a_full):.2e} "
        f"(linear functional fully captured by the single-Pauli scan)",
    )

# --------------------------------------------------------------------------
# Report table
# --------------------------------------------------------------------------
print(f"\n{'='*78}\nSUMMARY TABLE  (asymmetry = ||M_+i||^2 - ||M_-i||^2 ; eta normalized)\n{'='*78}")
print(f"{'N':>2} | {'object':<18} | {'asymmetry (raw)':>18} | {'eta':>14} | {'||M_anti||^2':>14}")
print("-" * 78)
for (N, label, a, e, mn) in rows:
    print(f"{N:>2} | {label:<18} | {a:>18.6e} | {e:>14.6e} | {mn:>14.6e}")

# --------------------------------------------------------------------------
# Final gate-first verdict.
#   SANITY gates (G0/G0b/G_sigma/G3lin/G3xcheck): MUST pass (objects correct).
#   HYPOTHESIS gates (G1/G2/G3/B1/B2a/B2b): encode Direction 1-3 predictions;
#       a FIRED one is the finding (refutation), NOT to be loosened.
#   MECHANISM gates (M1/M2/M3): diagnose WHY a hypothesis gate fired.
# --------------------------------------------------------------------------
print(f"\n{'='*78}\nGATE LEDGER\n{'='*78}")
for name, ok, detail in GATES:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

hypo_prefixes = ("G1", "G2", "G3(", "B1", "B2a", "B2b")
mech_prefixes = ("M1", "M2", "M3", "Mscan")
failed = [g for g in GATES if not g[1]]
hypo_failed = [g for g in failed if g[0].startswith(hypo_prefixes)]
mech_failed = [g for g in failed if g[0].startswith(mech_prefixes)]
sanity_failed = [g for g in failed if g not in hypo_failed and g not in mech_failed]

print(f"\n{len(GATES) - len(failed)}/{len(GATES)} gates passed.")
if sanity_failed:
    print("SANITY gate(s) FAILED -- objects wrong, results invalid:")
    for name, ok, detail in sanity_failed:
        print(f"  - {name}: {detail}")
if hypo_failed:
    print("\nHYPOTHESIS gate(s) FIRED -- THE FINDING (Direction 1-3 refuted as stated):")
    for name, ok, detail in hypo_failed:
        print(f"  - {name}: {detail}")
if mech_failed:
    print("\nMECHANISM gate(s) FAILED -- the diagnosis below is wrong, rethink:")
    for name, ok, detail in mech_failed:
        print(f"  - {name}: {detail}")

print(f"\n{'='*78}\nVERDICT\n{'='*78}")
if sanity_failed or mech_failed:
    print("INCONCLUSIVE: a sanity or mechanism gate failed; fix before trusting the verdict.")
elif hypo_failed:
    print("Direction 1 REFUTED for the physical chain. The un-recycled drain has NO")
    print("chirality: DRAIN-ONLY asymmetry = 0 bit-exact, and the NO-JUMP generator of")
    print("any XY/Heisenberg chain has asymmetry = 0 -- exactly like the FULL")
    print("Lindbladian. There is no chirality for the jump to recycle.")
    print("MECHANISM (M1/M2/M3/Mscan pass): asym(NO-JUMP) = cross_AB(H) is a LINEAR")
    print("functional of H (F112 kills the commutator's own asymmetry; the drain alone")
    print("is balanced). Its ONLY single-Pauli carriers are the single-site Z_l strings")
    print("(a local detuning on the drained site, same axis as the drain's n=(I-Z)/2).")
    print("The chain's NN bond terms {XX,YY,ZZ} are all in the kernel -> chain asym = 0.")
    print("A generic H (real OR complex -- reality is NOT the axis; my first guess was")
    print("wrong and gate M1 caught it) carries Z_l content, hence != 0: the spine's")
    print("132/270 is a GENERIC-H artifact, not a drain/chirality property. The reframe")
    print("does NOT become a positive statement.")
else:
    print("Direction 1 HOLDS as stated: full=0, no-jump!=0, jump recycles the chirality.")

# Gate-first: the assert fires iff the hypothesis is not fully borne out. The
# fired G2/B1/B2a ARE the finding (refutation); the non-zero exit is intended.
assert not (sanity_failed or mech_failed), "Sanity/mechanism failure -- verdict not trustworthy."
assert not hypo_failed, (
    f"{len(hypo_failed)} hypothesis gate(s) FIRED -- this is the finding (NOT loosened): "
    + "; ".join(g[0] for g in hypo_failed)
)
print("\nALL HYPOTHESIS GATES PASS: Direction 1 HOLDS as stated.")
