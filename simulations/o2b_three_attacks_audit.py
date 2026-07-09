"""Self-asserting verifier for the F89 section "Three attacks on the sign law" (2026-07-09).

Companion to experiments/F89_SEED_EXISTENCE_REDUCTION.md. Recomputes from below, and ASSERTS,
every machine-checkable claim of that section (promoted from the empty-review auditor's script that
reproduced the section number for number):

  A. Kernel facts (baseline, N = 5, 7; NOT an all-N law -- the doc's own Piece-2 table has
     nullity 21 != 15 at N = 11): ker(K66) subset class O, baseline nullity 3(N-1)/2;
     ker(K22) subset E, nullity N-1 (this one IS the Piece-1 theorem);
     ker(K) baseline nullity 3(N-1)/2 splitting (N-1)/2 in E + (N-1) in O.
  B. Band sweep (simple real band eigenvectors, lam in (-6,-2), gauge conj(v) = T v):
     - slope identity  dlam/dq = -2 x^T K y / v^T v   (vs finite differences)
     - branch identities  kappa_-2 = sigma[(lam+6) - q dlam/dq]/4,
                          kappa_-6 = sigma[q dlam/dq - (lam+2)]/4     (machine precision)
     - separating law kappa_-2 < 0 => sigma < 0  (0 violations)
     - the STRONGER cell form b' >= b => a' >= a  (0 violations; strictly stronger than the
       separating law, see the doc)
     - cell identities  ||x2||^2 - ||y2||^2 = q^2 (a'-a)/(lam+2)^2  and
                        sigma ||v||^2 = q^2 [(a'-a)/(lam+2)^2 + (b'-b)/(lam+6)^2]
     - dead ends stay dead: monotonicity sigma*dlam/dq < 0 HAS violations; the J-cone law
       v^dag(iTK)v < 0 HAS violations, all on T-negative (sigma < 0) vectors.
  C. Seed-local claims at every defective seed (4 at N = 5, 7 at N = 7, via the committed
     o2b_krein_sign_law machinery):
     - kappa_-2 > 0 and kappa_-2 + kappa_-6 ~ 0
     - the Feshbach two-component system W x2 = q Ktil y2, W y2 = -q Ktil^T x2 with
       W = |lam+2| + Sigma_even, Sigma_even = -q^2 (lam+6) K26 R6 K62,
       R6 = ((lam+6)^2 + q^2 K66^2)^-1, Ktil = K22 - q^2 K26 K66 R6 K62;  W class-diagonal
     - the ratio identity ||x2||^2/||y2||^2 = (beta - |lam+2|)/(|lam+2| - alpha)
     - the binding-order FLIP alpha > |lam+2| > beta at exactly the three named seeds
       (N=5 q*=1.286074; N=7 q*=1.077392, 1.447107) -- the Semenoff-dominance falsifier
     - the surviving proximity form ||lam+2|-alpha| < ||lam+2|-beta| at every seed
     - W_E and W_O indefinite at every seed
     - adj(lam* I - L) has rank 1 and kappa_-2 = -S6/ST with S6 = tr(P_-6 adj), ST = tr(T adj),
       adj taken on the FULL block; N=5 full-block sign patterns, seeds sorted by q*:
       S6 = (+,-,-,+), ST = (-,+,+,-), c = (-,+,+,-)  (per-R-sector adjugates flip signs seed-wise,
       the product S6*ST < 0 is convention-independent)
     - Q-pencil Q = T(D - lam) + i q T K:  nu = nu' = 0 at the seed (tangential T-neutral
       touching), dnu/dq = -(4/q) kappa_-2 (the descending law), and the opening pattern:
       exits (complex pair born, count-drops) have nu'' < 0, re-entrant entries nu'' > 0.
  D. S6, ST are real and even in q at generic points, and S6 is NOT sign-definite on the strip.

Run:  python simulations/o2b_three_attacks_audit.py         # parts A (N=5,7), B (N=5), C (N=5), D
      python simulations/o2b_three_attacks_audit.py 7       # part C at N=7 (seed hunt, ~minutes)
      python simulations/o2b_three_attacks_audit.py b7      # part B at N=7
      python simulations/o2b_three_attacks_audit.py b9      # part B at N=9 (the separating-law leg)

Axis: unit-hop (this file's q = twice the octic census q*). Sample counts of band sweeps are
grid-dependent; the LAWS asserted here are grid-free, the counts printed are this grid's.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_existence_nullity_check import build                      # noqa: E402
from o2b_krein_sign_law import (bipartite_sign, r_sectors,          # noqa: E402
                                find_transitions, analyze_transition, classify)

FLIPS = {(5, 1.286074), (7, 1.077392), (7, 1.447107)}   # the three Semenoff-dominance flips


def nullity_split(M, Emask, tol=1e-9):
    """(total nullity, rank of kernel restricted to class E, to class O)."""
    if M.size == 0:
        return 0, 0, 0
    _, s, vt = np.linalg.svd(M)
    ker = vt[s < tol]
    tot = ker.shape[0]
    nE = np.linalg.matrix_rank(ker[:, Emask], tol=1e-8) if tot and Emask.any() else 0
    nO = np.linalg.matrix_rank(ker[:, ~Emask], tol=1e-8) if tot and (~Emask).any() else 0
    return tot, nE, nO


def part_A(N):
    """Kernel class facts. Baseline nullities asserted only at N = 5, 7 (N = 11 resonates: 21 != 15)."""
    A, C = build(N)
    K = (C / 1j).real
    t = bipartite_sign(N)
    m2, m6 = np.isclose(A, -2.0), np.isclose(A, -6.0)
    E = t > 0
    tot22, nE22, nO22 = nullity_split(K[np.ix_(m2, m2)], E[m2])
    tot66, nE66, nO66 = nullity_split(K[np.ix_(m6, m6)], E[m6])
    totK, nEK, nOK = nullity_split(K, E)
    print(f"[A] N={N}: ker(K22)={tot22} (E:{nE22},O:{nO22}); ker(K66)={tot66} (E:{nE66},O:{nO66}); "
          f"ker(K)={totK} (E:{nEK},O:{nOK})")
    assert tot22 == N - 1 and nO22 == 0, "Piece-1 theorem: ker(K22) = N-1, entirely class E"
    if N in (5, 7):
        assert tot66 == 3 * (N - 1) // 2 and nE66 == 0, "baseline: ker(K66) = 3(N-1)/2, entirely O"
        assert (totK, nEK, nOK) == (3 * (N - 1) // 2, (N - 1) // 2, N - 1), \
            "baseline: ker(K) = 3(N-1)/2 = (N-1)/2 E + (N-1) O"


def gauge_fix(v, t, tol=1e-8):
    """Scale v so conj(v) = T v (T = diag(t)); returns (v', residual)."""
    tv = t * v
    i = np.argmax(np.abs(v))
    theta = np.angle(np.conj(v[i]) / tv[i]) / 2
    vp = np.exp(1j * theta) * v
    res = np.linalg.norm(np.conj(vp) - t * vp) / np.linalg.norm(vp)
    if res > tol:
        vp2 = 1j * vp
        res2 = np.linalg.norm(np.conj(vp2) - t * vp2) / np.linalg.norm(vp2)
        if res2 < res:
            vp, res = vp2, res2
    return vp, res


def band_samples(N, qlist, gap_tol=1e-5):
    """Cell data for every simple real band eigenvector of D + qC over qlist."""
    A, C = build(N)
    D = np.diag(A)
    K = (C / 1j).real
    t = bipartite_sign(N)
    m2, m6 = np.isclose(A, -2.0), np.isclose(A, -6.0)
    out, h, n_gauge_drop = [], 1e-6, 0
    for q in qlist:
        w, V = np.linalg.eig(D + q * C)
        wp = np.linalg.eigvals(D + (q + h) * C)
        wm = np.linalg.eigvals(D + (q - h) * C)
        for k in range(len(w)):
            lam = w[k]
            if abs(lam.imag) > 1e-9 or not (-6 + 1e-6 < lam.real < -2 - 1e-6):
                continue
            lam = lam.real
            dists = np.abs(w - w[k]); dists[k] = np.inf
            if dists.min() < gap_tol:
                continue
            vp, gres = gauge_fix(V[:, k], t)
            if gres > 1e-6:
                n_gauge_drop += 1        # a simple real eigenvector must gauge-fix; report drops
                continue
            x, y = vp.real, vp.imag
            vTv = (vp @ vp).real
            n2 = np.vdot(vp, vp).real
            sigma = np.vdot(vp, t * vp).real / n2
            kap2 = float(np.sum(t[m2] * np.abs(vp[m2]) ** 2) / n2)
            kap6 = float(np.sum(t[m6] * np.abs(vp[m6]) ** 2) / n2)
            lamdot = -2 * (x @ K @ y) / vTv
            lp = wp[np.argmin(np.abs(wp - lam))]
            lm = wm[np.argmin(np.abs(wm - lam))]
            Kx, Ky = K @ x, K @ y
            out.append(dict(q=q, lam=lam, sigma=sigma, kap2=kap2, kap6=kap6, lamdot=lamdot,
                            lamdot_fd=((lp - lm) / (2 * h)).real,
                            a=np.sum(Kx[m2] ** 2), b=np.sum(Kx[m6] ** 2),
                            ap=np.sum(Ky[m2] ** 2), bp=np.sum(Ky[m6] ** 2),
                            x2n=np.sum(x[m2] ** 2), y2n=np.sum(y[m2] ** 2), nrm2=n2,
                            jval=np.vdot(vp, 1j * (t[:, None] * K) @ vp).real / n2))
    return out, n_gauge_drop


def part_B(N):
    S, dropped = band_samples(N, np.linspace(0.11, 20.0, 41))
    print(f"[B] N={N}: {len(S)} simple real band samples ({dropped} gauge-fix drops)")
    assert dropped <= max(2, len(S) // 50), "too many gauge-fix drops: the sweep is not clean"
    e_slope = max(abs(s['lamdot'] - s['lamdot_fd']) for s in S)
    e_br2 = max(abs(s['kap2'] - s['sigma'] * ((s['lam'] + 6) - s['q'] * s['lamdot']) / 4) for s in S)
    e_br6 = max(abs(s['kap6'] - s['sigma'] * (s['q'] * s['lamdot'] - (s['lam'] + 2)) / 4) for s in S)
    e_id1 = max(abs((s['x2n'] - s['y2n']) - s['q'] ** 2 * (s['ap'] - s['a']) / (s['lam'] + 2) ** 2)
                for s in S)
    e_id2 = max(abs(s['sigma'] * s['nrm2'] - s['q'] ** 2 * ((s['ap'] - s['a']) / (s['lam'] + 2) ** 2
                                                            + (s['bp'] - s['b']) / (s['lam'] + 6) ** 2))
                for s in S)
    sep = [s for s in S if s['kap2'] < -1e-12 and s['sigma'] >= 0]
    strong = [s for s in S if s['bp'] >= s['b'] - 1e-12 and s['ap'] < s['a'] - 1e-12]
    # NOTE: sigma*lamdot < 0 and v^dag(iTK)v < 0 are the SAME pointwise condition (both equal
    # -2 x^T K y / ||v||^2 up to a positive factor); both dead ends are kept for readability, but
    # they are one falsification, not two.
    mono = [s for s in S if s['sigma'] * s['lamdot'] >= 1e-10]
    jv = [s for s in S if s['jval'] >= 1e-10]
    print(f"    slope err {e_slope:.1e}; branch-identity errs {e_br2:.1e}/{e_br6:.1e}; "
          f"cell-identity errs {e_id1:.1e}/{e_id2:.1e}")
    print(f"    separating law: {len(sep)} violations; strong cell form: {len(strong)}; "
          f"monotonicity violations {len(mono)}; J-cone violations {len(jv)} (all sigma<0: "
          f"{all(s['sigma'] < 0 for s in jv)})")
    assert e_slope < 1e-4 and e_br2 < 1e-10 and e_br6 < 1e-10 and e_id1 < 1e-9 and e_id2 < 1e-9
    assert not sep, "separating law violated"
    assert not strong, "strong cell form b'>=b => a'>=a violated"
    assert mono, "monotonicity dead end unexpectedly clean (should have violations)"
    assert jv and all(s['sigma'] < 0 for s in jv), "J-cone violations must exist, all on sigma<0"


def seed_list(N, qmax=20.0, ngrid=8000):
    seeds = []
    for label, A_s, C_s, t_s in r_sectors(N):
        seen = []
        for q_lo, q_hi, lam0, born in sorted(find_transitions(A_s, C_s, qmax, ngrid)):
            r = analyze_transition(A_s, C_s, t_s, q_lo, q_hi, lam0, born)
            if any(abs(r['q'] - sq) < 1e-4 for sq in seen):
                continue
            seen.append(r['q'])
            if classify(r) == "defective":
                seeds.append(dict(label=label, A=A_s, C=C_s, t=t_s, **r))
    return seeds


def part_C(N, expected_seeds):
    print(f"[C] N={N}: seed hunt per R-parity sector ...")
    seeds = sorted(seed_list(N), key=lambda d: d['q'])
    print(f"    {len(seeds)} defective seeds (expected {expected_seeds})")
    assert len(seeds) == expected_seeds
    A_full, C_full = build(N)                     # S6, ST live on the FULL block (the doc's L)
    t_full = bipartite_sign(N)
    m6_full = np.isclose(A_full, -6.0)
    n_full = len(A_full)
    signs_S6, signs_ST, signs_c, flips = [], [], [], []
    for sd in seeds:
        A_s, C_s, t_s, q, lam = sd['A'], sd['C'], sd['t'], sd['q'], sd['lam0']
        n = len(A_s)
        assert sd['kappa2'] > 0.02 and abs(sd['kappa2'] + sd['kappa6']) < 1e-8
        D = np.diag(A_s)
        M = (D + q * C_s) - lam * np.eye(n)
        _, sv, vh = np.linalg.svd(M)
        vp, gres = gauge_fix(vh[-1, :].conj(), t_s)
        assert gres < 1e-6
        K = ((C_s / 1j) + (C_s / 1j).conj().T).real / 2
        m2, m6 = np.isclose(A_s, -2.0), np.isclose(A_s, -6.0)
        K22, K26, K66 = K[np.ix_(m2, m2)], K[np.ix_(m2, m6)], K[np.ix_(m6, m6)]
        R6 = np.linalg.inv((lam + 6) ** 2 * np.eye(K66.shape[0]) + q ** 2 * (K66 @ K66))
        Sigma = -q ** 2 * (lam + 6) * (K26 @ R6 @ K26.T)
        W = abs(lam + 2) * np.eye(Sigma.shape[0]) + Sigma
        Ktil = K22 - q ** 2 * (K26 @ K66 @ R6 @ K26.T)
        x2, y2 = vp.real[m2], vp.imag[m2]
        nx2, ny2 = x2 @ x2, y2 @ y2
        assert max(np.linalg.norm(W @ x2 - q * Ktil @ y2),
                   np.linalg.norm(W @ y2 + q * Ktil.T @ x2)) / np.linalg.norm(x2) < 1e-8
        Em2 = t_s[m2] > 0
        assert np.abs(W[np.ix_(Em2, ~Em2)]).max() < 1e-10, "W must be class-diagonal"
        alpha = x2 @ (-Sigma) @ x2 / nx2
        beta = y2 @ (-Sigma) @ y2 / ny2
        d2 = abs(lam + 2)
        assert abs(nx2 / ny2 - (beta - d2) / (d2 - alpha)) < 1e-6, "ratio identity"
        assert abs(d2 - alpha) < abs(d2 - beta), "proximity form"
        if alpha > d2 > beta:
            flips.append(round(q, 6))
        eE = np.linalg.eigvalsh(W[np.ix_(Em2, Em2)])
        eO = np.linalg.eigvalsh(W[np.ix_(~Em2, ~Em2)])
        assert eE.min() < 0 < eE.max() and eO.min() < 0 < eO.max(), "W_E, W_O indefinite"
        # adjugate on the FULL block: rank 1, kappa = -S6/ST, via SVD product
        Mf = lam * np.eye(n_full) - (np.diag(A_full) + q * C_full)
        U, s_all, Vh = np.linalg.svd(Mf)
        assert s_all[-1] < 1e-8 < s_all[-2], "adjugate rank 1 (geometric multiplicity 1)"
        phi = np.linalg.det(U @ Vh) * np.prod(s_all[:-1])
        vn, unc = Vh[-1, :].conj(), U[:, -1].conj()
        S6 = phi * np.sum(vn[m6_full] * unc[m6_full])
        ST = phi * np.sum(t_full * vn * unc)
        assert abs(-(S6 / ST).real - sd['kappa2']) < 1e-4, "kappa_-2 = -S6/ST"
        vf, _ = gauge_fix(vn, t_full)
        signs_S6.append(int(np.sign(S6.real)))
        signs_ST.append(int(np.sign(ST.real)))
        signs_c.append(int(np.sign((ST / np.vdot(vf, vf)).real)))
        # Q-pencil: tangential T-neutral touching, descending law, opening pattern
        Tm = np.diag(t_s)

        def nu(l, qq, u0):
            w, V = np.linalg.eigh(Tm @ (D - l * np.eye(n)) + 1j * qq * (Tm @ K))
            return w[np.argmax(np.abs(u0.conj() @ V))]

        w0, V0 = np.linalg.eigh(Tm @ (D - lam * np.eye(n)) + 1j * q * (Tm @ K))
        k0 = np.argmin(np.abs(w0))
        u0 = V0[:, k0]
        hl, hq = 2e-3, 1e-5
        nup, num = nu(lam + hl, q, u0), nu(lam - hl, q, u0)
        nu1 = (nup - num) / (2 * hl)
        nu2 = (nup - 2 * w0[k0] + num) / hl ** 2
        dnudq = (nu(lam, q + hq, u0) - nu(lam, q - hq, u0)) / (2 * hq)
        assert abs(w0[k0]) < 1e-4 and abs(nu1) < 5e-3, "tangential touching nu = nu' = 0"
        assert abs(dnudq - (-(4 / q) * sd['kappa2'])) < 1e-3, "descending law dnu/dq = -(4/q) kappa_-2"
        assert dnudq < 0, "T-neutral touchings descend"
        # sd['born']: True = complex pair born (exit of the real pair, a count-drop)
        assert (nu2 < 0) == bool(sd['born']), "exits open downward, re-entrant entries upward"
        print(f"  {sd['label']:>6} q*={q:.6f} lam*={lam:.6f} {'exit' if sd['born'] else 'entry'} "
              f"k2={sd['kappa2']:+.6f} alpha/|l+2|/beta {alpha:.4f}/{d2:.4f}/{beta:.4f} "
              f"flip={alpha > d2 > beta} nu''={nu2:+.3f} dnu/dq={dnudq:+.5f}")
    assert set((N, f) for f in flips) == {fl for fl in FLIPS if fl[0] == N}, \
        f"binding-order flips must be exactly the named ones; got {flips}"
    if N == 5:
        # FULL-block adjugate, seeds sorted by q*. (Per-R-sector adjugates give (+,+,-,+)/(-,-,+,-):
        # the co-sector characteristic factor flips S6 and ST together, so the product is
        # convention-independent while the individual signs are not.)
        assert signs_S6 == [1, -1, -1, 1] and signs_ST == [-1, 1, 1, -1] and signs_c == [-1, 1, 1, -1]
    assert all(a * b < 0 for a, b in zip(signs_S6, signs_ST)), "sign law <=> S6*ST < 0"
    print(f"    sign patterns (by q*): S6 {signs_S6}  ST {signs_ST}  c {signs_c}")


def part_D(N=5):
    A, C = build(N)
    t = bipartite_sign(N)
    m6 = np.isclose(A, -6.0)
    D = np.diag(A)
    n = len(A)
    for (lam, q) in [(-3.7, 0.9), (-4.9, 2.3), (-2.5, 5.0)]:
        vals = {}
        for qq in (q, -q):
            M = lam * np.eye(n) - (D + qq * C)
            adj = np.linalg.det(M) * np.linalg.inv(M)
            vals[qq] = (np.trace(adj[np.ix_(m6, m6)]), np.trace(np.diag(t) @ adj))
        (S6p, STp), (S6m, STm) = vals[q], vals[-q]
        assert abs(S6p.imag) < 1e-6 * abs(S6p) and abs(STp.imag) < 1e-6 * abs(STp), "S6, ST real"
        assert abs(S6p - S6m) < 1e-6 * abs(S6p) and abs(STp - STm) < 1e-6 * abs(STp), "even in q"
    signs = set()
    for lam in np.linspace(-5.9, -2.1, 30):
        for q in np.linspace(0.1, 6.0, 30):
            M = lam * np.eye(n) - (D + q * C)
            v = np.linalg.det(M) * np.trace(np.linalg.inv(M)[np.ix_(m6, m6)])
            if abs(v.real) > 1e-6 * abs(v):
                signs.add(int(np.sign(v.real)))
    assert signs == {-1, 1}, "S6 must be mixed-sign on the strip (dead end)"
    print(f"[D] N={N}: S6, ST real + even in q at generic points; S6 mixed-sign on strip. OK")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "5"
    if which == "5":
        part_A(5)
        part_A(7)
        part_B(5)
        part_D(5)
        part_C(5, 4)
        print("ALL ASSERTIONS PASSED (N=5 full + N=7 kernels)")
    elif which == "7":
        part_C(7, 7)
        print("ALL ASSERTIONS PASSED (N=7 seeds)")
    elif which == "b7":
        part_B(7)
        print("ALL ASSERTIONS PASSED (N=7 band)")
    elif which == "b9":
        part_B(9)
        print("ALL ASSERTIONS PASSED (N=9 band)")
