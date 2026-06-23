#!/usr/bin/env python3
r"""F112/F87 orthogonality derived: scope inclusion + mechanism separation + the F113 bridge (committed verifier).

The orthogonality between F112 (polarity balance) and F87 (dissipator-resonance trichotomy)
was an empirical line since 2026-05-26 (`polarity_probe_f87_connection.py`: all three classes
give asymmetry = 0 at N=3). This verifier upgrades it to a derived statement in three parts:

  (a) SCOPE INCLUSION. Every F87 input (Hermitian Pauli H-commutator + pure Z-dephasing) is an
      instance of PROOF_F112's hypotheses (Hermitian H, bit_b-homogeneous c: each Z_l is a
      single Pauli string). So F112's asymmetry is identically zero on F87's ENTIRE domain,
      all three trichotomy classes. The two functionals never co-vary because one is
      identically zero where the other lives.
  (b) MECHANISM SEPARATION. For bit_b-odd H (X^(x)N H X^(x)N = -H; covers the diagonal Klein
      cell that hosts every F87-hard pair, F110) the dagger involution IS the windowed
      converse's first reflection: M^dag = Fcal M Fcal with Fcal = X^(x)N (x) X^(x)N, diff
      0.00e+00. F112's functional exhausts that dagger structure at DEGREE 2 (Frobenius norms
      of Pi-eigenprojections, Step 5 of PROOF_F112); the F87 hardness decision lives in the
      second reflection R + the unsigned girth at ODD degree (p_3, p_5, ...), invisible to
      degree-2 norms. Same grading, different degree.
  (c) SCOPED ONE-WAY IMPLICATION (sigma-/sigma+ probe family). Off the shared domain, with
      amplitude damping added, both functionals read the SAME moment t_1^(l) = Tr(Z_l H)
      = 2^N c_l = 2^(N-1) omega_l: F113 reads it linearly paired with the net rate,
      asym = 2^N * Sum_l t_1^(l) * (gamma_pump,l - gamma_T1,l), while the girth ladder's
      ell=1 deg-1 face squares it, p_3 = 6 * gamma * Sum_l (t_1^(l))^2. Hence balance-broken
      (generic rates) => some t_1 != 0 => p_3 > 0 => F87-hard. One-way only: flux/K3 hard
      pairs have all t_1 = 0 and keep the balance under T1.

Convention reconciliation (2026-06-10): the hunt draft printed pred = -asym/2 because it used
the Z_l coefficient c_l where the documented F113 form (PROOF_F113, ANALYTICAL_FORMULAS F113,
LindbladBitBPiBreakMagnitude.cs) uses omega_l with H ⊇ Sum_l (omega_l/2) Z_l, i.e.
omega_l = 2 c_l, and the standard-physics sign sigma- = |0><1| (lowering; cooling with
omega > 0 gives NEGATIVE asymmetry): asym = (4^N/2) * Sum_l omega_l * (gamma_pump - gamma_T1).
With that chain the prediction matches the measured asymmetry to machine precision
(rel residual ~1e-15 at N=3; honest residual printed side by side below, not exactly equal
as floats).

Block ledger
------------
  Block A   shared involution       : M^dag = Fcal M Fcal diff 0.00e+00 for the two
                                      diagonal-cell pairs (XXZ+ZXX, IXY+XIY); off-cell H
                                      (Heisenberg, XY+Z-drive) fails the Fcal-form (diff 8, 4)
                                      but odd power sums stay real regardless.
  Block B+C scope inclusion         : pure Z-deph, F112 asym = 0.0 for truly/soft/hard while
                                      p_3 fires 6*4^N*Sum c_l^2*gamma for the Z-drive case.
  Block D   off-domain F113 channel : sigma- T1 added; documented-convention F113 prediction
                                      vs measured asymmetry side by side (machine precision).
  Block E   the shared moment       : t_1 scan at eps = 0.5/1.0/2.0; F87 squares it (p_3),
                                      F113 pairs it linearly with the rate (asym).
  Block F   N=4 scoping             : K3 XXZ+XZX and flux IXY+XIY are N=3-SOFT (their m* = 9
                                      first fires at N=4); at N=4 they are hard and F112
                                      asym is still 0. Scopes the statement honestly.

Run: python simulations/f112_f87_orthogonality.py   (~1 min; Block F's N=4 dominates).
Supersedes the hunt drafts _f112_f87_orthogonality_check.py + _f112_f87_orth_supplement.py.
2026-06-10.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw  # noqa: E402
from framework.lindblad import lindbladian_pauli_dephasing  # noqa: E402
from framework.pauli import _build_bilinear, _build_kbody_chain, site_op  # noqa: E402
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_hc  # noqa: E402

N = 3
GAMMA = 0.1
GT1 = 0.03
D = 2 ** N
chain = fw.ChainSystem(N=N, gamma_0=GAMMA)


def H_bilinear(terms, J=1.0):
    return _build_bilinear(N, chain.bonds, [(a, b, J) for a, b in terms])


def H_kbody(terms, J=1.0, n=N):
    return _build_kbody_chain(n, [tuple(t) + (J,) for t in terms])


def recentered_M_vec(H, gammas_z, n=N):
    """M = L + sigma*I in vec basis (the windowed-converse object)."""
    L = lindbladian_pauli_dephasing(H, gammas_z, dephase_letter='Z')
    sigma = sum(gammas_z)
    d = 2 ** n
    return L + sigma * np.eye(d * d, dtype=complex)


def odd_power_sums(M, ms=(3, 5)):
    out = {}
    P = M.copy()
    k = 1
    for m in ms:
        while k < m:
            P = P @ M
            k += 1
        out[m] = np.trace(P)
    return out


def sminus(l, n=N):
    sm = np.array([[0, 1], [0, 0]], dtype=complex)   # sigma- = |0><1| (lowering, standard physics)
    ops = [np.eye(2, dtype=complex)] * n
    ops[l] = sm
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def f113_prediction(H, gamma_t1, gamma_pump=0.0, n=N):
    """Documented-convention F113 closed form.

    H contains Sum_l (omega_l/2) Z_l, so omega_l = 2 * c_l with c_l = Tr(Z_l H) / 2^n.
    asym = (4^n / 2) * Sum_l omega_l * (gamma_pump,l - gamma_T1,l)  [uniform rates here].
    """
    d = 2 ** n
    omegas = [2.0 * np.trace(site_op(n, l, 'Z') @ H).real / d for l in range(n)]
    return (4 ** n / 2) * sum(o * (gamma_pump - gamma_t1) for o in omegas), omegas


X1 = np.array([[0, 1], [1, 0]], dtype=complex)
F = X1
for _ in range(N - 1):
    F = np.kron(F, X1)
Fcal = np.kron(F, F)           # the two-sided reflection (script-F, the windowed converse's first)

Z_OPS = [site_op(N, l, 'Z') for l in range(N)]

print("=" * 76)
print("BLOCK A: shared involution  M^dag = Fcal M Fcal  (dagger == first reflection)")
print("=" * 76)
cases = {
    'truly  Heisenberg':      (H_bilinear([('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]), False),
    'soft   XXZ+ZXX (k=3)':   (H_kbody([('X', 'X', 'Z'), ('Z', 'X', 'X')]), True),
    'hard   flux IXY+XIY':    (H_kbody([('I', 'X', 'Y'), ('X', 'I', 'Y')]), True),
    'hard   Zdrive XYchain':  (H_bilinear([('X', 'X'), ('Y', 'Y')], J=1.0)
                               + 1.0 * site_op(N, 0, 'Z') + 2.0 * site_op(N, 1, 'Z')
                               + 1.0 * site_op(N, 2, 'Z'), False),
}
for name, (H, bitb_odd) in cases.items():
    M = recentered_M_vec(H, [GAMMA] * N)
    diff = np.max(np.abs(M.conj().T - Fcal @ M @ Fcal))
    ps = odd_power_sums(M)
    tag = 'bit_b-odd H (diagonal cell)' if bitb_odd else 'mixed-cell H (Fcal-form n/a)'
    print(f"  {name:24s}  max|M^dag - Fcal M Fcal| = {diff:.3e}   "
          f"Im p3 = {ps[3].imag:+.3e}  Im p5 = {ps[5].imag:+.3e}   [{tag}]")
    if bitb_odd:
        assert diff == 0.0, f"Fcal dagger relation must be exact on bit_b-odd H, got {diff}"
    assert abs(ps[3].imag) < 1e-9 and abs(ps[5].imag) < 1e-9, "odd power sums must be real"

print()
print("=" * 76)
print("BLOCK B+C: scope inclusion (pure Z-deph): F112 asym vs F87 class / power sums")
print("=" * 76)
pair_cases = [
    ('truly', [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')], 'bilinear'),
    ('soft',  [('X', 'X', 'Z'), ('Z', 'X', 'X')],   'k3'),
    ('soft (N=3; hard from N=4)', [('I', 'X', 'Y'), ('X', 'I', 'Y')], 'k3'),
]
for expect, terms, kind in pair_cases:
    cls = fw.classify_pauli_pair(chain, terms)
    pol = fw.polarity_coordinates(chain, terms, strict=False)
    H = H_kbody(terms) if kind == 'k3' else H_bilinear(terms)
    M = recentered_M_vec(H, [GAMMA] * N)
    ps = odd_power_sums(M)
    print(f"  terms={str(terms):44s} class={cls:6s} (expect {expect})")
    print(f"      F112 asym = {pol['asymmetry']:+.3e}   "
          f"p3 = {ps[3].real:+.4e}  p5 = {ps[5].real:+.4e}")
    assert pol['asymmetry'] == 0.0, "F112 asymmetry must be exactly 0 on F87's domain"
# Z-drive case (single-site terms; classify_pauli_pair does not take them,
# hardness certified by the cell-free m=3 face p3 = 6*4^N*Sum c_l^2*gamma)
H_zd = cases['hard   Zdrive XYchain'][0]
M_zd = recentered_M_vec(H_zd, [GAMMA] * N)
ps = odd_power_sums(M_zd)
c_sq = 1.0 ** 2 + 2.0 ** 2 + 1.0 ** 2
pred_p3 = 6 * 4 ** N * c_sq * GAMMA
pol_zd = polarity_coordinates_from_hc(H_zd, Z_OPS, [GAMMA] * N, N)
print(f"  Zdrive XY chain (c=(1,2,1)):                 hard via p3")
print(f"      F112 asym = {pol_zd['asymmetry']:+.3e}   "
      f"p3 = {ps[3].real:+.4e}  (pred 6*4^N*Sum c^2*g = {pred_p3:.4e}, "
      f"rel diff {abs(ps[3].real - pred_p3) / pred_p3:.1e})")
assert pol_zd['asymmetry'] == 0.0, "F112 asymmetry must be exactly 0 on F87's domain"
assert abs(ps[3].real - pred_p3) / pred_p3 < 1e-12, "p3 must match the cell-free m=3 face"

print()
print("=" * 76)
print(f"BLOCK D: off-domain (sigma- T1 added, gamma_t1={GT1}): F113 break channel")
print("=" * 76)
t1_cases = {
    'hard flux IXY+XIY  + T1 (no Z-drive)': cases['hard   flux IXY+XIY'][0],
    'soft XXZ+ZXX       + T1 (no Z-drive)': cases['soft   XXZ+ZXX (k=3)'][0],
    'truly Heisenberg   + T1 (no Z-drive)': cases['truly  Heisenberg'][0],
    'XY chain + Z-drive + T1':              H_zd,
}
for name, H in t1_cases.items():
    c_ops = Z_OPS + [sminus(l) for l in range(N)]
    gs = [GAMMA] * N + [GT1] * N
    pol = polarity_coordinates_from_hc(H, c_ops, gs, N, sigma=N * GAMMA)
    meas = pol['asymmetry']
    pred, omegas = f113_prediction(H, GT1)
    resid = abs(meas - pred)
    scale = max(abs(pred), 1.0)
    print(f"  {name:38s} measured = {meas:+.15e}")
    print(f"  {'':38s} F113 pred= {pred:+.15e}   |diff| = {resid:.1e}   "
          f"omega = {[round(o, 3) for o in omegas]}")
    assert resid / scale < 1e-12, f"F113 prediction must match measurement, resid {resid}"

print()
print("=" * 76)
print("BLOCK E: the shared moment t1_l = Tr(Z_l H): F87 squares it, F113 pairs it")
print("=" * 76)
for eps in (0.5, 1.0, 2.0):
    H = H_bilinear([('X', 'X'), ('Y', 'Y')]) + eps * site_op(N, 1, 'Z')
    M = recentered_M_vec(H, [GAMMA] * N)
    p3 = odd_power_sums(M, ms=(3,))[3].real
    c_ops = Z_OPS + [sminus(l) for l in range(N)]
    pol = polarity_coordinates_from_hc(H, c_ops, [GAMMA] * N + [GT1] * N, N, sigma=N * GAMMA)
    meas = pol['asymmetry']
    t1m = np.trace(site_op(N, 1, 'Z') @ H).real           # t_1 = 2^N * eps
    pred_p3 = 6 * 4 ** N * eps ** 2 * GAMMA               # F87 deg-1 face (quadratic in t_1)
    pred_p3_moment = 6 * GAMMA * (t1m / 2 ** N) ** 2 * 4 ** N   # same face via t_1 = 2^N c
    pred_asym = 2 ** N * t1m * (0.0 - GT1)                # F113 via t_1 (linear, rate-paired)
    print(f"  eps={eps:4.1f}: t1=Tr(Z_1 H)={t1m:6.1f}")
    print(f"      p3   = {p3:+.6e}   pred 6*4^N*eps^2*g     = {pred_p3:+.6e}   "
          f"(via t1: {pred_p3_moment:+.6e})")
    print(f"      asym = {meas:+.15e}   pred 2^N*t1*(gp-gt1) = {pred_asym:+.15e}   "
          f"|diff| = {abs(meas - pred_asym):.1e}")
    assert abs(p3 - pred_p3) / pred_p3 < 1e-12
    assert abs(pred_p3 - pred_p3_moment) < 1e-9 * pred_p3
    assert abs(meas - pred_asym) / max(abs(pred_asym), 1.0) < 1e-12

print()
print("=" * 76)
print("BLOCK F: diagonal-cell hard pairs at their true m* vs F112 asym (N=3 vs N=4)")
print("=" * 76)
expected_class = {
    3: {'soft  XXZ+ZXX': 'soft', 'hard  K3 XXZ+XZX m*=9': 'soft',
        'hard  multiZ XXZ+ZZZ m*=5': 'hard', 'hard  flux IXY+XIY m*=9': 'soft'},
    4: {'soft  XXZ+ZXX': 'soft', 'hard  K3 XXZ+XZX m*=9': 'hard',
        'hard  multiZ XXZ+ZZZ m*=5': 'hard', 'hard  flux IXY+XIY m*=9': 'hard'},
}
for n in (3, 4):
    d = 2 ** n
    chain_n = fw.ChainSystem(N=n, gamma_0=GAMMA)
    print(f"  --- N={n} ---")
    for name, terms in [
        ('soft  XXZ+ZXX',             [('X', 'X', 'Z'), ('Z', 'X', 'X')]),
        ('hard  K3 XXZ+XZX m*=9',     [('X', 'X', 'Z'), ('X', 'Z', 'X')]),
        ('hard  multiZ XXZ+ZZZ m*=5', [('X', 'X', 'Z'), ('Z', 'Z', 'Z')]),
        ('hard  flux IXY+XIY m*=9',   [('I', 'X', 'Y'), ('X', 'I', 'Y')]),
    ]:
        H = _build_kbody_chain(n, [tuple(t) + (1.0,) for t in terms])
        M = recentered_M_vec(H, [GAMMA] * n, n=n)
        ps = {}
        P = M.copy()
        k = 1
        for m in range(2, 12):
            P = P @ M
            k += 1
            if k % 2 == 1:
                ps[k] = np.trace(P)
        cls = fw.classify_pauli_pair(chain_n, terms)
        pol = fw.polarity_coordinates(chain_n, terms, strict=False)
        pstr = "  ".join(f"p{m}={ps[m].real:+.3e}" for m in (3, 5, 7, 9, 11))
        print(f"    {name:30s} class={cls:6s} F112 asym={pol['asymmetry']:+.1e}")
        print(f"        {pstr}")
        assert cls == expected_class[n][name], \
            f"{name} at N={n}: expected {expected_class[n][name]}, got {cls}"
        assert pol['asymmetry'] == 0.0, "F112 asymmetry must be exactly 0 on F87's domain"

print()
print("All blocks passed: scope inclusion (asym = 0.0 exactly on F87's domain), mechanism")
print("separation (dagger = Fcal on bit_b-odd H, exact), and the F113 one-way bridge")
print("(documented convention, machine-precision match). 2026-06-10.")
