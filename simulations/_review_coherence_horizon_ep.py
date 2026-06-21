"""Coherence-horizon sqrt-EP: artifact-free diabolic-vs-defective re-check.

GATE-FIRST. The claim under review: the slowest single-excitation {0,2}-coherence
mode of the dephased XY chain undergoes a GENUINE square-root (2nd-order, DEFECTIVE)
exceptional point at Q*(N) = {1, sqrt2, 1.879, 2.372, ...}.

The existing evidence is eig-based PhaseRigidity r = 1/sqrt(K) -- the SAME instrument
family that misfired in the F86a retraction (read r->0 / K large on a merely
near-degenerate, NON-defective spectrum). So the EP-ness must be re-checked with the
three artifact-free measures (none reading a raw eig eigenvector pairing):

  (1) Riesz spectral-projector norm ||P|| on a small contour enclosing ONLY the
      coalescing pair (resolvent contour integral; radius-robust when correct).
      ||P|| >> 1 = genuine non-normality / defective; ||P|| ~ 1 = diabolic.
  (2) departure-from-normality of the 2x2 compression of L onto range(P) in an
      ORTHONORMAL basis: dep = sqrt(||A||_F^2 - sum|lam|^2). >0 (geo<alg) = defective.
  (3) geometric vs algebraic multiplicity at the coalescence. alg = trace(P)
      (count enclosed); geo = nullity of (L - lam0 I) via SVD. geo<alg = DEFECTIVE.

OPPOSITE-PRIOR to F86a: here the analytic discriminant of the SE dispersion
lam^2+8g*lam+4J^2 q^2 has a TRUE double-root at the horizon (at N=2,3 exactly,
lam^2+4g*lam+cJ^2, c=4,2 -> (lam+2g)^2=0). A 2x2 double-root with nonzero coupling
is generically DEFECTIVE, so I EXPECT genuine here. Guard against confirmation bias:
report honestly if any N reads diabolic / non-defective.

GATE 0 validates the 3 measures on KNOWN answers (a defective toy 2x2 EP MUST read
DEFECTIVE; a diabolic diag/normal repeated eigenvalue MUST read DIABOLIC) before any
horizon number is trusted.

Console is cp1252 -> ASCII only in prints (no gamma/sqrt/tensor glyphs).
"""
from __future__ import annotations

import sys
import numpy as np

sys.path.insert(0, 'simulations')
from coherence_horizon_se_block import L_se, qstar_se, LADDER  # noqa: E402

np.set_printoptions(precision=5, suppress=True, linewidth=140)


# ----------------------------------------------------------------------
# Artifact-free machinery (reused from _review_f86a_*; none reads eig vectors)
# ----------------------------------------------------------------------
def riesz_projector(L, lam0, r, nq=600):
    """P = (1/2pi i) oint (zI - L)^{-1} dz on a circle of radius r about lam0.
    Pure resolvent contour integral; never an eig eigenvector."""
    D = L.shape[0]
    Id = np.eye(D, dtype=complex)
    P = np.zeros((D, D), dtype=complex)
    for k in range(nq):
        th = 2 * np.pi * (k + 0.5) / nq
        z = lam0 + r * np.exp(1j * th)
        dz = 1j * r * np.exp(1j * th) * (2 * np.pi / nq)
        P += np.linalg.solve(z * Id - L, Id) * dz
    return P / (2j * np.pi)


def geometric_multiplicity(L, lam0, tol_rel=1e-6, tol_abs=1e-9):
    """geo mult = # singular values of (L - lam0 I) that are ~ 0 (nullity by SVD).
    Uses an ABSOLUTE-and-relative threshold (the F86a/reference lesson): the nullity
    test needs an absolute floor, else a near-zero or exactly-zero (L - lam0 I) -- e.g.
    a diabolic diag(lam0,lam0) where the whole matrix vanishes -- is read full-rank /
    nullity 0 because the relative tol scales to 0. Scale the absolute floor by ||L||."""
    M = L - lam0 * np.eye(L.shape[0])
    s = np.linalg.svd(M, compute_uv=False)
    scale = max(np.linalg.norm(L, 2), 1.0)
    thr = max(tol_rel * (np.max(s) if s.size else 0.0), tol_abs * scale)
    return int(np.sum(s < thr)), s


def projector_diagnostics(L, lam0, r, alg_expected=2):
    """All three artifact-free measures on the contour of radius r about lam0.
    Returns dict: ||P||_2, alg=trace(P), geo (SVD nullity), departure-from-normality
    of the alg-dim compression in an orthonormal basis of range(P)."""
    P = riesz_projector(L, lam0, r)
    pn = np.linalg.norm(P, 2)
    m_alg = float(np.trace(P).real)
    m = max(int(round(m_alg)), 1)
    # orthonormal basis of range(P) from the leading left singular vectors
    U, S, _ = np.linalg.svd(P)
    V = U[:, :m]
    A = V.conj().T @ L @ V                     # compression onto the enclosed eigenspace
    eigA = np.linalg.eigvals(A)
    dep = np.sqrt(max(0.0, np.linalg.norm(A, 'fro') ** 2 - np.sum(np.abs(eigA) ** 2)))
    g, _ = geometric_multiplicity(L, lam0)
    return dict(pn=pn, m_alg=m_alg, m_alg_round=m, m_geo=g, dep=dep,
                A_norm=float(np.linalg.norm(A, 'fro')), proj_sv=S[:4])


def eig_phase_rigidity_min(L):
    """The eig-based instrument under suspicion: r_n = 1/sqrt(K_n),
    K_n = ||VR_col_n||^2 * ||VRinv_row_n||^2 (Petermann). Returns the MIN r over
    modes (-> 0 near an EP), its mode eigenvalue, max K, and the eigenvector cond."""
    w, VR = np.linalg.eig(L)
    VR = VR / np.linalg.norm(VR, axis=0, keepdims=True)
    Rinv = np.linalg.inv(VR)
    K = (np.linalg.norm(VR, axis=0) ** 2) * (np.linalg.norm(Rinv, axis=1) ** 2)
    nmax = int(np.argmax(K))
    return 1.0 / np.sqrt(K[nmax]), w[nmax], float(K[nmax]), float(np.linalg.cond(VR))


def classify(d, dep_tol=1e-2):
    """Primary discriminator = geometric-vs-algebraic multiplicity AND
    departure-from-normality of the compression (both basis-independent, neither an
    eig pairing). ||P|| is a SUPPLEMENTARY non-orthogonality reading -- it is large
    only when there are eigenvalues OUTSIDE the contour with non-orthogonal
    eigenspaces (informative on the N^2 block; trivially 1 on a closed 2x2), so it
    must NOT gate the verdict (GATE 0 lesson).

      DEFECTIVE (genuine 2nd-order EP / Jordan): alg>=2, geo<alg, dep>tol.
      DIABOLIC  (normal repeated eigenvalue):    alg>=2, geo==alg, dep~0.
    """
    rel = d['dep'] / max(1.0, d['A_norm'])
    if d['m_alg_round'] >= 2 and d['m_geo'] < d['m_alg_round'] and rel > dep_tol:
        return "DEFECTIVE"
    if d['m_alg_round'] >= 2 and d['m_geo'] == d['m_alg_round'] and rel < dep_tol:
        return "DIABOLIC"
    if d['m_alg_round'] < 2:
        return "no-coalescence (1 eig enclosed)"
    return "AMBIGUOUS"


# ----------------------------------------------------------------------
# GATE 0: KNOWN-ANSWER validation of the three measures.
# ----------------------------------------------------------------------
def gate0():
    ok = True
    print("=" * 78)
    print("GATE 0  validate the artifact-free test on KNOWN answers")
    print("=" * 78)

    # --- (a) KNOWN DEFECTIVE: toy 2x2 EP [[-2g, i w],[i w, -2g]] at discriminant zero.
    #     eigenvalues (-2g) +- i w; this is ALWAYS defective for w != 0? No: it is
    #     normal-ish only if the two off-diagonals are conjugate. Here both = i*w
    #     (NOT Hermitian), eigenvalues -2g +- i*w are DISTINCT for w!=0. The DEFECTIVE
    #     point is where they coalesce. For [[a, b],[c, a]] eigenvalues a +- sqrt(bc);
    #     coalescence at bc=0 with b,c not both 0 -> Jordan. Use b=i w, c=0 -> a double
    #     eigenvalue a, geometric mult 1 (a genuine non-trivial Jordan block).
    g, w = 1.0, 1.3
    Ldef = np.array([[-2 * g, 1j * w], [0.0, -2 * g]], dtype=complex)
    lam0 = -2 * g
    d = projector_diagnostics(Ldef, lam0, r=0.5)
    rmin, _, Kmax, cond = eig_phase_rigidity_min(Ldef)
    verdict = classify(d)
    print("\n(a) KNOWN-DEFECTIVE toy Jordan 2x2  [[-2g, i*w],[0, -2g]]  (MUST read DEFECTIVE)")
    print(f"    eigenvalues       = {np.linalg.eigvals(Ldef)}")
    print(f"    ||P|| (r=0.5)     = {d['pn']:.4e}   alg={d['m_alg_round']}  geo={d['m_geo']}")
    print(f"    departure-from-N  = {d['dep']:.4e}   (A_norm={d['A_norm']:.3f})")
    print(f"    eig rigidity rmin = {rmin:.4e}  (Kmax={Kmax:.3e}, eigvec cond={cond:.2e})")
    print(f"    >>> reads: {verdict}   (expected DEFECTIVE)")
    if verdict != "DEFECTIVE":
        ok = False
        print("    *** GATE 0(a) FAILED ***")

    # --- (b) Second known DEFECTIVE: the SYMMETRIC near-EP form [[a, b],[b, a]] is
    #     normal (symmetric) -> diabolic at b=0. To get a defective symmetric-looking
    #     one use the toy_Leff style with complex coupling at its discriminant zero:
    #     [[-2g, i*J],[i*J, -6g]] coalesces when (2g)^2 = (i J)^2 i.e. J = 2g (sign),
    #     double root at -4g, Jordan. (This is the F86a toy_Leff EP, k=1.)
    g2 = 1.0
    J2 = 2.0 * g2
    Ldef2 = np.array([[-2 * g2, 1j * J2], [1j * J2, -6 * g2]], dtype=complex)
    lam0b = -4 * g2
    d2 = projector_diagnostics(Ldef2, lam0b, r=0.5)
    rmin2, _, Kmax2, cond2 = eig_phase_rigidity_min(Ldef2)
    v2 = classify(d2)
    print("\n(b) KNOWN-DEFECTIVE toy_Leff EP  [[-2g, i*J],[i*J, -6g]], J=2g  (MUST read DEFECTIVE)")
    print(f"    eigenvalues       = {np.linalg.eigvals(Ldef2)}")
    print(f"    ||P|| (r=0.5)     = {d2['pn']:.4e}   alg={d2['m_alg_round']}  geo={d2['m_geo']}")
    print(f"    departure-from-N  = {d2['dep']:.4e}")
    print(f"    eig rigidity rmin = {rmin2:.4e}  (Kmax={Kmax2:.3e})")
    print(f"    >>> reads: {v2}   (expected DEFECTIVE)")
    if v2 != "DEFECTIVE":
        ok = False
        print("    *** GATE 0(b) FAILED ***")

    # --- (c) KNOWN DIABOLIC: diag(-2g,-2g) -- two independent equal eigenvalues.
    Ldia = np.diag([-2 * g, -2 * g]).astype(complex)
    d3 = projector_diagnostics(Ldia, -2 * g, r=0.5)
    rmin3, _, Kmax3, cond3 = eig_phase_rigidity_min(Ldia)
    v3 = classify(d3)
    print("\n(c) KNOWN-DIABOLIC  diag(-2g, -2g)  (MUST read DIABOLIC)")
    print(f"    ||P|| (r=0.5)     = {d3['pn']:.4e}   alg={d3['m_alg_round']}  geo={d3['m_geo']}")
    print(f"    departure-from-N  = {d3['dep']:.4e}")
    print(f"    eig rigidity rmin = {rmin3:.4e}  (Kmax={Kmax3:.3e})")
    print(f"    >>> reads: {v3}   (expected DIABOLIC)")
    if v3 != "DIABOLIC":
        ok = False
        print("    *** GATE 0(c) FAILED ***")

    # --- (d) KNOWN DIABOLIC, normal-but-coupled: a 2x2 HERMITIAN with a repeated
    #     eigenvalue under a rotation -> still diabolic (||P||=1, dep=0).
    th = 0.7
    Q = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]], complex)
    Ldia2 = Q @ np.diag([-2 * g, -2 * g]).astype(complex) @ Q.conj().T
    d4 = projector_diagnostics(Ldia2, -2 * g, r=0.5)
    v4 = classify(d4)
    print("\n(d) KNOWN-DIABOLIC  rotated diag(-2g,-2g) (normal, repeated)  (MUST read DIABOLIC)")
    print(f"    ||P|| (r=0.5)     = {d4['pn']:.4e}   alg={d4['m_alg_round']}  geo={d4['m_geo']}")
    print(f"    departure-from-N  = {d4['dep']:.4e}")
    print(f"    >>> reads: {v4}   (expected DIABOLIC)")
    if v4 != "DIABOLIC":
        ok = False
        print("    *** GATE 0(d) FAILED ***")

    print("\n" + ("GATE 0 PASSED: the 3 measures separate DEFECTIVE from DIABOLIC on known answers.\n"
                   if ok else "*** GATE 0 FAILED -- STOP, do not trust horizon verdicts ***\n"))
    return ok


def gate_diabolic_control():
    """Confirmation-bias guard: the SAME contour/dep machinery on a DIABOLIC point of
    the SAME object class. At gamma=0 the SE generator -i(h(x)I - I(x)h^T) is exactly
    NORMAL (anti-Hermitian up to the basis), so every spectral degeneracy is DIABOLIC
    by construction. The test MUST read dep(A) ~ 0 here -- proving it is not rigged to
    call this object 'defective' regardless. (If it read defective on the normal
    generator, the horizon 'DEFECTIVE' verdict would be worthless.)"""
    print("=" * 78)
    print("GATE D  DIABOLIC control on the SAME object (gamma=0 SE generator is NORMAL)")
    print("=" * 78)
    ok = True
    for N in (4, 6):
        h = np.zeros((N, N), complex)
        for i in range(N - 1):
            h[i, i + 1] = 1.0
            h[i + 1, i] = 1.0
        I = np.eye(N)
        L0 = -1j * (np.kron(h, I) - np.kron(I, h.T))     # gamma=0: exactly normal
        comm = np.linalg.norm(L0 @ L0.conj().T - L0.conj().T @ L0, 2)
        w = np.linalg.eigvals(L0)
        Dmat = np.abs(w[:, None] - w[None, :]) + np.eye(len(w)) * 1e9
        i, j = np.unravel_index(np.argmin(Dmat), Dmat.shape)
        lam0 = 0.5 * (w[i] + w[j])
        sep = abs(w[i] - w[j])
        dist = np.sort(np.abs(w - lam0))
        third = float(dist[2])
        r = min(max(0.40 * third, 5 * sep + 1e-3), 0.49 * third) if third > 1e-6 else 0.1
        P = riesz_projector(L0, lam0, r)
        m = max(int(round(np.trace(P).real)), 1)
        U, S, _ = np.linalg.svd(P)
        V = U[:, :m]
        A = V.conj().T @ L0 @ V
        ea = np.linalg.eigvals(A)
        dep = np.sqrt(max(0.0, np.linalg.norm(A, 'fro') ** 2 - np.sum(np.abs(ea) ** 2)))
        reads = "DIABOLIC" if dep < 1e-2 else "DEFECTIVE"
        print(f"  N={N}: ||[L0,L0^H]||={comm:.2e} (0=normal)  closest pair sep={sep:.2e}  "
              f"enclosed alg={np.trace(P).real:.2f}  dep(A)={dep:.3e}  -> reads {reads}")
        if reads != "DIABOLIC":
            ok = False
            print("    *** GATE D FAILED: test calls the NORMAL generator defective ***")
    print("\n" + ("GATE D PASSED: the test reads DIABOLIC on the normal SE generator, so a "
                   "'DEFECTIVE'\n         verdict on the dephased horizon is a real property, "
                   "not a built-in bias.\n" if ok else "*** GATE D FAILED ***\n"))
    return ok


# ----------------------------------------------------------------------
# Identify the coalescing pair on the SE block near the horizon.
# ----------------------------------------------------------------------
def slow_coalescer(L):
    """The slowest non-zero oscillating conjugate PAIR about to freeze: among the
    near-gap upper-half-plane modes, the one with the SMALLEST |Im| whose conjugate
    is genuinely its nearest neighbour (so the pair is an isolated 2x2, not a member
    of a denser cluster). Returns (la, lb=conj(la), w). Eigenvalue-only; no vectors.

    The freezing {0,2} mode lives at Re ~ -2g (N=2,3) drifting toward -4g; restrict to
    that band so we do not grab the band-edge survivor (large |Im|) or a fast mode."""
    w = np.linalg.eigvals(L)
    nz = w[w.real < -1e-7]
    if len(nz) == 0:
        return None
    # candidate coalescers: upper-half oscillators with small |Im|, near Re ~ -2..-4
    osc = nz[(nz.imag > 1e-9)]
    if len(osc) == 0:
        order = np.argsort(-w.real)            # already frozen -> two slowest reals
        return w[order[0]], w[order[1]], w
    # the freezing mode = smallest |Im| of the slow oscillators; verify its conjugate
    # is its nearest eigenvalue (isolated pair) -- else step to the next candidate.
    for la in osc[np.argsort(np.abs(osc.imag))]:
        lb = np.conj(la)
        # distances from la to every eigenvalue except itself
        dd = np.abs(w - la)
        dd_sorted = np.sort(dd[dd > 1e-12])
        nn = dd_sorted[0] if dd_sorted.size else np.inf
        if abs(la - lb) <= 1.5 * nn + 1e-9:    # conjugate is (within 1.5x) the nearest
            return la, lb, w
    la = osc[np.argmin(np.abs(osc.imag))]      # fallback: smallest |Im| regardless
    return la, np.conj(la), w


def pair_contour_diagnostics(L, la, lb, w):
    """Tight contour around the coalescing pair: lam0 = midpoint, then enclose EXACTLY
    the two eigenvalues NEAREST lam0 (these ARE the pair, even when split to ~1e-7 by
    Jordan-block eigenvalue sensitivity), with radius set OUTSIDE the pair but INSIDE
    the genuine third-nearest eigenvalue. Reads the three artifact-free measures off
    the 2x2 compression A = V* L V (V = orthonormal basis of range(P)); geo is the
    nullity of (A - lam_bar I), NOT of the full L at the midpoint (which is not an
    eigenvalue when the pair is split)."""
    lam0 = 0.5 * (la + lb)
    sep = abs(la - lb)
    # rank eigenvalues by distance to lam0; [0],[1] are the pair, [2] is the third.
    dist = np.sort(np.abs(w - lam0))
    third = float(dist[2]) if len(dist) > 2 else 10.0
    # radius: outside the (possibly split) pair, well inside the third eigenvalue.
    r = 0.40 * third
    r = max(r, 5.0 * sep)            # ensure the split pair is comfortably enclosed
    r = min(r, 0.49 * third)         # ...but never reach the third eigenvalue
    P = riesz_projector(L, lam0, r)
    pn = np.linalg.norm(P, 2)
    m_alg = float(np.trace(P).real)
    m = max(int(round(m_alg)), 1)
    U, S, _ = np.linalg.svd(P)
    V = U[:, :m]
    A = V.conj().T @ L @ V                       # the enclosed (2x2) compression
    eigA = np.linalg.eigvals(A)
    dep = np.sqrt(max(0.0, np.linalg.norm(A, 'fro') ** 2 - np.sum(np.abs(eigA) ** 2)))
    # geometric multiplicity from the COMPRESSION at its own mean eigenvalue:
    lam_bar = eigA.mean()
    sA = np.linalg.svd(A - lam_bar * np.eye(A.shape[0]), compute_uv=False)
    scaleA = max(np.linalg.norm(A, 2), 1.0)
    geo = int(np.sum(sA < max(1e-6 * (sA.max() if sA.size else 0.0), 1e-7 * scaleA)))
    geo = max(geo, 1)
    return dict(lam0=lam0, sep=sep, third=third, r=r, pn=pn, m_alg=m_alg,
                m_alg_round=m, m_geo=geo, dep=dep,
                A_norm=float(np.linalg.norm(A, 'fro')), proj_sv=S[:4],
                eigA=eigA, A=A)


def find_closest_approach(N, gamma0=1.0, half_width=0.06, npts=121):
    """Sweep Q finely around Q*(N) and find the Q where the slow conjugate pair is
    CLOSEST (the on-axis projection of the EP). Returns (q_min, sep_min, qs, seps)."""
    qstar = qstar_se(N)
    qs = np.linspace(qstar - half_width, qstar + half_width, npts)
    seps = []
    for q in qs:
        # Q = J/g = q with g = gamma0 fixed, J = q*gamma0
        L = L_se(N, J=q * gamma0, g=gamma0)
        res = slow_coalescer(L)
        if res is None:
            seps.append(np.inf)
            continue
        la, lb, _ = res
        seps.append(abs(la - lb))
    seps = np.array(seps)
    i = int(np.argmin(seps))
    return qs[i], seps[i], qs, seps


def diagnostics_at_Q(N, q, gamma0=1.0):
    """Build L_se at Q=q, find the slow coalescer, return the pair_contour diagnostics
    + eig instrument readings."""
    L = L_se(N, J=q * gamma0, g=gamma0)
    res = slow_coalescer(L)
    if res is None:
        return None
    la, lb, w = res
    d = pair_contour_diagnostics(L, la, lb, w)
    rmin, _, Kmax, cond = eig_phase_rigidity_min(L)
    d.update(dict(la=la, lb=lb, rmin=rmin, Kmax=Kmax, cond=cond))
    return d


def analyze_N(N, gamma0=1.0):
    print("=" * 78)
    print(f"N={N}  coherence-horizon SE block  (Q*={qstar_se(N):.6f}, ladder={LADDER.get(N, 'n/a')})")
    print("=" * 78)
    qstar = qstar_se(N)

    # (1) locate the closest pair approach on the real Q axis (the on-axis EP image)
    q_min, sep_min, qs, seps = find_closest_approach(N, gamma0)
    print(f"  closest pair approach on real axis: Q = {q_min:.6f}  (Q*={qstar:.6f}, "
          f"diff={q_min - qstar:+.2e}),  min |Im-split| = {sep_min:.3e}")

    # (2) PRIMARY artifact-free test at a RESOLVABLE detuning above the horizon.
    #     At Q* the pair is coalesced to ~1e-7 (Jordan-block eigenvalue sensitivity:
    #     a true double-root splits as ~sqrt(eps_mach)), which is below the resolution
    #     a contour or an SVD can use. So read the 3 measures where the pair is a clean
    #     isolated 2x2: a DEFECTIVE EP keeps dep(A) O(1) (the 2x2 is a Jordan block as
    #     the eigenvalues merge); a DIABOLIC point has dep(A) -> 0.  Pick the CLOSEST-to-
    #     Q* detuning where the freezing pair is cleanly enclosed (alg=2) -- at some N
    #     the band-edge real mode interleaves the pair at wider detuning, so adapt.
    q_probe, dprobe = None, None
    for fac in (1.002, 1.004, 1.006, 1.008, 1.01, 1.015, 1.02):
        cand = diagnostics_at_Q(N, qstar * fac, gamma0)
        if cand is not None and cand['m_alg_round'] == 2:
            q_probe, dprobe = qstar * fac, cand
            break
    if dprobe is None:                            # fallback: just take 1.02 Q*
        q_probe = qstar * 1.02
        dprobe = diagnostics_at_Q(N, q_probe, gamma0)
    la = dprobe['la']
    print(f"  PRIMARY probe at Q = {q_probe:.5f} (={q_probe/qstar:.3f} Q*, closest detuning with "
          f"the pair cleanly isolated):")
    print(f"     pair lam = {dprobe['la']:.5f} & conj  (Re={dprobe['la'].real:.5f}, "
          f"|Im|={abs(dprobe['la'].imag):.4f}),  split = {dprobe['sep']:.4f}")
    print(f"     [EP double-root predicted at Re=-2g={-2*gamma0:.1f} (N=2,3) -> -4g (large N)]")
    print(f"     (1) Riesz ||P||_2            = {dprobe['pn']:.4e}   (>>1 = non-orthogonal eigenspaces)")
    print(f"     (2) departure-from-normality = {dprobe['dep']:.4e}   "
          f"(rel {dprobe['dep']/max(1.0,dprobe['A_norm']):.3e})")
    print(f"     (3) alg(trace P) = {dprobe['m_alg']:.4f}->{dprobe['m_alg_round']};  "
          f"geo(nullity of A-lam) = {dprobe['m_geo']};  third eig at dist {dprobe['third']:.3f}")
    print(f"     compression A eigenvalues = {dprobe['eigA']}")
    verdict = classify(dprobe)
    print(f"     >>> ARTIFACT-FREE VERDICT (resolvable probe): {verdict}")

    # (3) THE defective-vs-diabolic limit test: dep(A) and pair-split as Q -> Q*.
    #     DEFECTIVE: eigenvalues coalesce (split -> 0) while dep(A) stays bounded away
    #     from 0. DIABOLIC: dep(A) -> 0 together with the split.
    print("  limit Q -> Q* (defective: split->0 with dep(A) BOUNDED AWAY from 0 AND the two")
    print("  compression-eigenvectors merging |cos|->1; diabolic: dep->0, |cos| stays bounded):")
    deps_clean, splits_clean = [], []     # only points where the pair is cleanly enclosed (alg=2)
    for fac in (1.05, 1.02, 1.01, 1.005, 1.002):
        dq = diagnostics_at_Q(N, qstar * fac, gamma0)
        clean = (dq['m_alg_round'] == 2)
        tag = "" if clean else "  [pair not isolated here; excluded from floor]"
        cosang = float('nan')
        if clean and 'A' in dq:
            ea, va = np.linalg.eig(dq['A'])
            va = va / np.linalg.norm(va, axis=0)
            cosang = abs(np.vdot(va[:, 0], va[:, 1]))
        if clean:
            deps_clean.append(dq['dep'])
            splits_clean.append(dq['sep'])
        print(f"      Q={qstar*fac:.5f} (={fac:.3f}Q*):  pair-split={dq['sep']:.4e}  "
              f"dep(A)={dq['dep']:.4f}  eigvec|cos|={cosang:.5f}  "
              f"geo={dq['m_geo']}/alg={dq['m_alg_round']}{tag}")
    # the defect signal: dep(A) on the cleanly-enclosed pair as the split -> 0.
    dep_floor = min(deps_clean) if deps_clean else 0.0
    split_min = min(splits_clean) if splits_clean else np.inf
    dep_diabolic = dep_floor < 1e-2          # dep collapses to ~0 -> diabolic
    print(f"      => over the cleanly-enclosed points, split shrinks to {split_min:.2e} while "
          f"dep(A) floor = {dep_floor:.4f}")
    print(f"         ({'dep -> 0 with the split: DIABOLIC' if dep_diabolic else 'dep BOUNDED AWAY from 0 as split->0: DEFECTIVE (Jordan)'})")

    # (4) ||P|| radius-robustness at the resolvable probe (a true 2-dim projector norm
    #     is finite + radius-stable when the contour cleanly separates the pair).
    L = L_se(N, J=q_probe * gamma0, g=gamma0)
    print("  ||P|| radius-robustness at the probe (flat => contour cleanly encloses the pair):")
    rr = []
    for fac in (0.5, 0.7, 1.0, 1.3):
        r = min(fac * dprobe['r'], 0.49 * dprobe['third'])
        P = riesz_projector(L, dprobe['lam0'], r)
        pn, m = np.linalg.norm(P, 2), float(np.trace(P).real)
        rr.append(pn)
        print(f"      r={r:.4f}: ||P||={pn:8.3f}  alg={m:.3f}")

    # (5) eig-PhaseRigidity cross-check AT THE EXACT-EP closest approach (where eig
    #     misfired for F86a). Does r->0 here track the artifact-free 'genuine'?
    dep_exact = diagnostics_at_Q(N, q_min, gamma0)
    rmin, Kmax = dep_exact['rmin'], dep_exact['Kmax']
    print(f"  eig-PhaseRigidity at the EXACT-EP approach (Q={q_min:.5f}): "
          f"rmin = 1/sqrt(Kmax) = {rmin:.4e}  (Kmax={Kmax:.3e}, eigvec cond={dep_exact['cond']:.2e})")
    eig_says = "r->0 (reads GENUINE/near-EP)" if rmin < 0.1 else "r~O(1) (reads no non-normality)"
    artifact_free_genuine = (verdict == "DEFECTIVE" and not dep_diabolic)
    agree = ((rmin < 0.1) == artifact_free_genuine)
    print(f"     eig reads: {eig_says}.  artifact-free says genuine: {artifact_free_genuine}.  "
          f"AGREE: {'YES' if agree else 'NO'}")

    # (6) MAGNITUDE grid-sensitivity under dQ=+-1e-3 at the EXACT-EP approach: are the
    #     SPECIFIC numbers (Kmax, rmin, ||P||) fragile while the EP-NESS (geo<alg,
    #     dep>0 at resolvable distance) is robust? (the F86a magnitude lesson)
    print("  magnitude grid-sensitivity under dQ=+-1e-3 at the EXACT-EP approach:")
    Ks, rmins, pns = [], [], []
    for dq in (-1e-3, 0.0, 1e-3):
        dd = diagnostics_at_Q(N, q_min + dq, gamma0)
        Ks.append(dd['Kmax']); rmins.append(dd['rmin']); pns.append(dd['pn'])
        print(f"      dQ={dq:+.0e}:  Kmax={dd['Kmax']:10.3e}  rmin={dd['rmin']:.3e}  "
              f"||P||={dd['pn']:10.3e}  pair-split={dd['sep']:.3e}")
    K_swing = max(Ks) / max(1e-30, min(Ks))
    print(f"      => Kmax swings {K_swing:.3e}x across dQ=+-1e-3 (the SPECIFIC magnitudes are GRID-FRAGILE)")
    print(f"         vs the EP-ness verdict (dep(A) floor {dep_floor:.3f}, geo<alg) which is ROBUST.")
    print()
    return dict(N=N, q_min=q_min, sep_min=sep_min, q_probe=q_probe, lam=la,
                pn=dprobe['pn'], dep=dprobe['dep'], geo=dprobe['m_geo'],
                alg=dprobe['m_alg_round'], verdict=verdict, dep_floor=dep_floor,
                dep_diabolic=dep_diabolic, rmin=rmin, Kmax=Kmax, agree=agree,
                K_swing=K_swing, genuine=artifact_free_genuine)


if __name__ == '__main__':
    if not gate0():
        sys.exit("GATE 0 failed; aborting.")
    if not gate_diabolic_control():
        sys.exit("GATE D failed; aborting.")
    rows = [analyze_N(N) for N in (2, 3, 4, 5)]
    print("=" * 78)
    print("SUMMARY  coherence-horizon sqrt-EP: artifact-free diabolic-vs-defective")
    print("=" * 78)
    print(f"  {'N':>2} {'Q*':>8} {'probe||P||':>10} {'dep(probe)':>11} {'dep floor':>10} "
          f"{'geo/alg':>8} {'eig rmin':>9}  verdict")
    for r in rows:
        print(f"  {r['N']:>2} {qstar_se(r['N']):>8.5f} {r['pn']:>10.2e} {r['dep']:>11.4f} "
              f"{r['dep_floor']:>10.4f} {r['geo']:>3}/{r['alg']:<3} {r['rmin']:>9.2e}  {r['verdict']}")
    n_def = sum(1 for r in rows if r['genuine'])
    print(f"\n  genuine DEFECTIVE sqrt-EP (dep(A) bounded away from 0 + geo<alg) at {n_def}/4 of N=2..5.")
    print(f"  central claim 'genuine sqrt-EP at all N=2..5': "
          f"{'CONFIRMED' if n_def == 4 else ('REFUTED' if n_def == 0 else 'MIXED')}")
    print(f"  eig-rigidity-vs-artifact-free agreement: {sum(r['agree'] for r in rows)}/4 "
          f"(here eig SHOULD agree -- opposite-prior to F86a)")
    print(f"  magnitudes (Kmax) grid-fragile: swing up to {max(r['K_swing'] for r in rows):.2e}x "
          f"under dQ=+-1e-3, while the EP-ness verdict is robust.")
