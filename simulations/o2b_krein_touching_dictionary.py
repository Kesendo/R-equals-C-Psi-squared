"""The Krein-touching dictionary: nu'' is the sign characteristic (defectiveness), s6 is
the orthogonal jet nu_q -- the two are DIFFERENT invariants (F89 beta-exotic arc).

Companion gate for the "second gauge-free equivalent" paragraph of
experiments/F89_BETA_EXOTIC_GENERICITY.md (the Hermitian Krein pencil
Q(lambda, q) = T(A - lambda) + iq*T*K = T(L(q) - lambda)). Asserts, from below:

  T1  L(q) is T-selfadjoint (T L T = L^dagger) and Q is Hermitian for real (lambda, q);
      exact 0.0 residuals.
  T2  at each forced N = 5 seed, nu = 0 is a SIMPLE eigenvalue of Q(lambda*, q*)
      (multiplicity 1 = the geometric multiplicity of L at the EP2).
  T3  the touching is T-neutral: nu'(lambda*) = -v^dag T v / ||v||^2 = 0
      (|v^dag T v| < 1e-13; finite-difference nu' < 1e-8).
  T4  nu''(lambda*) = -2 v^dag T w / ||v||^2 with (L - lambda*) w = v, w the min-norm
      generalized vector (Rellich reduced-resolvent second-order formula; Q' = -T,
      Q'' = 0). v^dag T w is REAL (the GLR sign characteristic of the size-2 Jordan
      block) and the formula matches the finite-difference curvature of the Krein
      eigencurve at every seed to < 1e-4. GLR canonical-form nondegeneracy of the
      root-subspace Gram forces v^dag T w != 0 for every genuine size-2 block: the
      lambda-touching is quadratic, never flatter -- nu'' pins DEFECTIVENESS.
  T5  nu_q(lambda*) = i v^dag T K v / ||v||^2 = -(4/q*) kappa_-2, three ways
      (formula, finite difference in q, the doc's banked relation) to < 1e-4.
  T6  nu'' and kappa_-2 (= |s6| in the antilinear gauge) are NOT proportional:
      the ratio nu''/kappa_-2 across the four N = 5 seeds spans sign and magnitude
      (asserted spread > 5, including one sign flip).
  T7  THE WITNESS (the moment-relation kill witness of the doc, 4x4): a genuine
      T-selfadjoint size-2 Jordan block with nu'' = -0.5 != 0 AND
      s6 = r^T C r = kappa_-2 = nu_q = 0 exactly. Hence nu'' != 0 (the sign
      characteristic) does NOT imply s6 != 0: the GLR machinery decides the
      touching order (H1-type defectiveness content) and is BLIND to O2b.

Conclusion the doc paragraph carries: the Krein pencil restates kappa_-2 as the
q-jet of the touching and supplies no new lever for the O2b nonvanishing; nu''
signs at seeds (this gate computes the four N = 5 ones; the 11-seed windows live in
o2b_krein_sign_law.py's arc) are a defectiveness/sign-characteristic diagnostic,
not evidence toward s6 != 0.

Axis: unit-hop q (twice the octic-book q), as in o2b_krein_sign_law.py.
Run:  python simulations/o2b_krein_touching_dictionary.py     (~1 min)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_existence_nullity_check import build                      # noqa: E402
from o2b_krein_sign_law import r_sectors, find_transitions          # noqa: E402


def refine_seed(A, C, q_lo, q_hi, lam0, re_w=0.6, im_w=0.6, iters=200):
    """Golden-section the closest-pair gap in the bracket -> (q*, lam*, gap)."""
    D = np.diag(A.astype(complex))

    def box_pair(qq):
        w = np.linalg.eigvals(D + qq * C)
        box = w[(np.abs(w.real - lam0) < re_w) & (np.abs(w.imag) < im_w)]
        if box.size < 2:
            return np.inf, None, None
        best = (np.inf, None, None)
        for a in range(box.size):
            for b in range(a + 1, box.size):
                dab = abs(box[a] - box[b])
                if dab < best[0]:
                    best = (dab, box[a], box[b])
        return best

    scan = np.linspace(q_lo, q_hi, 200)
    gv = np.array([box_pair(q)[0] for q in scan])
    j = int(np.argmin(gv))
    lo, hi = scan[max(0, j - 1)], scan[min(len(scan) - 1, j + 1)]
    gr = (np.sqrt(5) - 1) / 2
    c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    fc, fd = box_pair(c)[0], box_pair(d)[0]
    for _ in range(iters):
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - gr * (hi - lo)
            fc = box_pair(c)[0]
        else:
            lo, c, fc = c, d, fd
            d = lo + gr * (hi - lo)
            fd = box_pair(d)[0]
    q = 0.5 * (lo + hi)
    gap, la, lb = box_pair(q)
    lam = 0.5 * (la + lb).real if la is not None else lam0
    return q, lam, gap


def analyze(A, C, t, q, lam, tag):
    """All dictionary entries at one refined (q*, lam*); asserts T1-T5 locally."""
    n = len(A)
    D = np.diag(A.astype(complex))
    K = (C / 1j).real
    assert np.abs(C - 1j * K).max() < 1e-12
    T = np.diag(t.astype(float))
    L = D + q * C
    I = np.eye(n)

    assert np.max(np.abs(T @ L @ T - L.conj().T)) == 0.0, f"{tag}: not T-selfadjoint"  # T1
    Qmat = T @ (L - lam * I)
    assert np.max(np.abs(Qmat - Qmat.conj().T)) == 0.0, f"{tag}: Q not Hermitian"      # T1

    M = L - lam * I
    U, S, Vh = np.linalg.svd(M)
    v = Vh[-1].conj()
    v /= np.linalg.norm(v)
    assert S[-2] > 1e-2, f"{tag}: second singular value small (not geom mult 1)"

    evQ = np.linalg.eigvalsh(Qmat)
    assert int(np.sum(np.abs(evQ) < 1e-6)) == 1, f"{tag}: nu=0 not simple"             # T2

    vTv = np.vdot(v, T @ v).real
    assert abs(vTv) < 1e-13, f"{tag}: touching not T-neutral"                          # T3

    w = np.linalg.pinv(M, rcond=1e-8) @ v
    assert np.linalg.norm(M @ w - v) < 1e-12, f"{tag}: v not in range (not EP2)"
    vTw = np.vdot(v, T @ w)
    assert abs(vTw.imag) < 1e-12, f"{tag}: sign characteristic not real"               # T4
    nu2_formula = float(-2.0 * vTw.real)

    def branch_nu(lam_):
        ev = np.linalg.eigvalsh(T @ (L - lam_ * I))
        return ev[np.argmin(np.abs(ev))]

    h = 1e-4
    nu1_fd = (branch_nu(lam + h) - branch_nu(lam - h)) / (2 * h)
    nu2_fd = (branch_nu(lam + h) - 2 * branch_nu(lam) + branch_nu(lam - h)) / h ** 2
    assert abs(nu1_fd) < 1e-8, f"{tag}: nu' != 0"                                      # T3
    assert abs(nu2_formula - nu2_fd) < 1e-4, f"{tag}: nu'' formula vs FD"              # T4
    assert abs(nu2_formula) > 0.1, f"{tag}: nu'' unexpectedly small"                   # T4

    m2 = np.isclose(A, -2.0)
    kap2 = float(np.sum(t[m2] * np.abs(v[m2]) ** 2))
    nu_q_formula = float((1j * np.vdot(v, T @ (K @ v))).real)
    nu_q_doc = -(4.0 / q) * kap2

    def branch_nu_q(qq):
        ev = np.linalg.eigvalsh(T @ (D + qq * C - lam * I))
        return ev[np.argmin(np.abs(ev))]

    hq = 1e-4
    nu_q_fd = (branch_nu_q(q + hq) - branch_nu_q(q - hq)) / (2 * hq)
    assert abs(nu_q_formula - nu_q_doc) < 1e-10, f"{tag}: nu_q != -(4/q) kappa_-2"     # T5
    assert abs(nu_q_formula - nu_q_fd) < 1e-4, f"{tag}: nu_q formula vs FD"            # T5

    print(f"  {tag:>8} q*={q:8.5f} lam*={lam:8.4f}  kap2={kap2:+.4f}  "
          f"nu''={nu2_formula:+.4f} (FD {nu2_fd:+.4f})  nu_q={nu_q_formula:+.4f} "
          f"= -(4/q)kap2 {nu_q_doc:+.4f}  ratio nu''/kap2 = {nu2_formula / kap2:+.2f}")
    return nu2_formula, kap2


def run_seeds_n5():
    N = 5
    print(f"T1-T6 at the N = {N} forced seeds (unit-hop axis):")
    ratios = []
    n_seeds = 0
    for label, A_s, C_s, t_s in r_sectors(N):
        seen = []
        for q_lo, q_hi, lam0, born in sorted(find_transitions(A_s, C_s, 20.0, 8000)):
            q, lam, gap = refine_seed(A_s, C_s, q_lo, q_hi, lam0)
            if any(abs(q - s) < 1e-3 for s in seen) or gap > 1e-3:
                continue
            seen.append(q)
            nu2, kap2 = analyze(A_s, C_s, t_s, q, lam, label + ("+" if born else "-"))
            ratios.append(nu2 / kap2)
            n_seeds += 1
    assert n_seeds == 4, f"expected the 4 forced N=5 seeds, found {n_seeds}"
    spread = max(ratios) - min(ratios)
    assert spread > 5 and max(ratios) * min(ratios) < 0, \
        f"nu''/kap2 ratios unexpectedly uniform: {ratios}"                             # T6
    print(f"  T6: nu''/kappa_-2 spans {min(ratios):+.2f} .. {max(ratios):+.2f} "
          f"(sign flip included): nu'' and s6 are DIFFERENT invariants")


def run_witness():
    print("T7, the in-class witness (genuine size-2 Jordan, nu'' != 0, s6 = 0):")
    K = np.array([[0, 2, 0, 0], [2, 0, -4, 0], [0, -4, 0, -2], [0, 0, -2, 0]], float)
    A = np.array([-2.0, -2.0, -6.0, -6.0])
    t = np.array([1.0, -1.0, 1.0, -1.0])
    q, lam = 1.0, -4.0
    D = np.diag(A.astype(complex))
    C = 1j * K
    T = np.diag(t)
    L = D + q * C
    I = np.eye(4)
    assert np.max(np.abs(T @ L @ T - L.conj().T)) == 0.0
    assert np.max(np.abs(T @ K @ T + K)) == 0.0
    M = L - lam * I
    U, S, Vh = np.linalg.svd(M)
    assert int(np.sum(S < 1e-9)) == 1, "witness not geometric multiplicity 1"
    v = Vh[-1].conj()
    v /= np.linalg.norm(v)
    m6 = np.isclose(A, -6.0)
    s6 = complex(np.sum(v[m6] ** 2))
    rCr = complex(v @ (C @ v))
    assert abs(s6) < 1e-14 and abs(rCr) < 1e-14, "witness s6 not zero"
    w = np.linalg.pinv(M, rcond=1e-8) @ v
    assert np.linalg.norm(M @ w - v) < 1e-12, "witness v not in range (not EP2)"
    vTw = np.vdot(v, T @ w)
    nu2 = float(-2.0 * vTw.real)
    assert abs(nu2 + 0.5) < 1e-9, f"witness nu'' = {nu2}, expected -0.5"
    m2 = np.isclose(A, -2.0)
    kap2 = float(np.sum(t[m2] * np.abs(v[m2]) ** 2))
    nu_q = float((1j * np.vdot(v, T @ (K @ v))).real)
    assert abs(kap2) < 1e-14 and abs(nu_q) < 1e-14, "witness nu_q not zero"
    print(f"  size-2 Jordan (sing vals {np.round(S, 4)}), nu'' = {nu2:+.4f} != 0, "
          f"s6 = 0, nu_q = 0: the sign characteristic is blind to O2b")


if __name__ == "__main__":
    run_seeds_n5()
    run_witness()
    print("\nPASS: nu'' = -2 * (GLR sign characteristic) pins DEFECTIVENESS "
          "(forced != 0 at size-2 blocks); s6 lives in the orthogonal jet "
          "nu_q = -(4/q*) kappa_-2; the witness separates them. "
          "The Krein pencil restates kappa_-2 and adds no O2b lever.")
