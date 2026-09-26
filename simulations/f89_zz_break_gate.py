"""The N=4 XXZ-Delta response of the F89 path-3 octic crossing, located and certified.

At Delta=0 the (SE,DE) block has a semisimple crossing (diabolic) at
q_EP = sqrt((-1+sqrt(13))/6), lam_EP = -4g + 2iJ: g1 = g2 = 2, the 2x2 restriction
L|2D = lam*I. Turning on an XXZ anisotropy Delta (the ZZ term makes the two magnons
interact, so E_DE is no longer E_a + E_b) changes the pencil. This script asks what
happens to the crossing, and answers it with a certificate rather than a small gap.

THE LOCATOR (complex q, discriminant Newton)
  f(q) = (lam_a(q) - lam_b(q))^2 for the eigenvalue pair of the R-even block nearest the
  running midpoint. f is a symmetric function of the pair, hence holomorphic in q while the
  pair stays isolated, the coalescence included: a semisimple crossing is a DOUBLE zero of f
  (the gap closes linearly), and a SIMPLE zero forces a Jordan EP2 (the gap closes as sqrt|q - q*|);
  the kinds are told apart by f together with geo against alg, not by f alone.
  Newton runs on u = f/f' (step -u/u', u' = 1 - f f''/f'^2), quadratic at both orders; the
  second zero is found with the first divided out of f. A pattern search on the pair gap
  cannot do this job at an EP2: the square-root law leaves its gap near 1e-4 at |dq| ~ 1e-8.

THE CERTIFICATE (per located zero)
  - the pair gap after Newton, compared with the strict 1e-6 coincidence bound (its law: the
    rounding floor is 2*sqrt(dep*u*|M|) at an EP2, about u*|M| at a crossing);
  - alg = round(tr P), geo = nullity of A - mean(eig A) and dep on the Riesz compression A
    (gate-validated on known 2x2s in GATE 0); character by geo against alg, not by a
    relative departure threshold;
  - the zero order, the winding number of f on a circle of radius 1e-7 about the zero, closed on
    its opening value (so an integer to rounding) and read only if every step turns f by less than
    pi/4; a crossing needs order 2, an EP2 order 1 and a departure above 1e-6*max(|A|_F, 1);
  - on the real axis, a sign change of f itself across q* -+ 1e-9: f is real there because
    the spectrum at real q is closed under lam -> -8 - conj(lam), so a sign change proves an
    exact real-q coalescence between the two points.

CONSTRUCTION (robust, Pauli-built -- no hand-rolled magnon interaction)
  H(Delta) = J*sum_{i=0..2}(X_i X_{i+1} + Y_i Y_{i+1})
           + J*Delta*sum_{i=0..2} Z_i Z_{i+1}        (open 4-chain, 3 bonds)
  M_SE(Delta) = H restricted to the 1-excitation sector (4x4, incl. ZZ diagonal)
  M_DE(Delta) = H restricted to the 2-excitation sector (6x6, incl. ZZ diagonal)
  L = -i*M_SE(ket) + i*M_DE(bra) + diag(-2g if SE-site in DE-pair else -6g)
  Dephasing is Delta-INDEPENDENT.  g = 1, J = q.

  The (SE,DE) block is built two independent ways and CROSS-CHECKED at Stage 0:
   (i)  Pauli-sector assembly (above), reduced at Delta=0 to the committed reference build;
   (ii) the popcount-(1,2) sub-block of the FULL 256x256 XXZ Liouvillian.

THE MECHANISM TABLE restricts the dephasing D and the coherent part H to the Riesz 2-plane of
each located zero (HS-orthonormal R-even basis) and reports how far each is from a scalar.

Output: stdout, and the same text in simulations/results/f89_zz_break_gate.txt (written by
main(), not at import). ASCII console prints only. numpy + scipy.
"""
from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

np.set_printoptions(precision=6, suppress=True, linewidth=160)

# --------------------------------------------------------------------------
# Pauli operators and the 4-site XXZ Hamiltonian (site 0 = leftmost/MSB factor)
# --------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

N = 4
D = 2 ** N
STRICT_FULL_BLOCK_TOL = 1e-6
DE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
BASIS = [(i, jk) for i in range(4) for jk in DE_PAIRS]            # 24 (SE site, DE pair)
INDEX = {b: n for n, b in enumerate(BASIS)}
PERM = {0: 3, 1: 2, 2: 1, 3: 0}                                  # chain reflection s -> 3-s

# computational-basis index of a config: site s contributes bit 2^(N-1-s)
SE_IDX = [2 ** (N - 1 - i) for i in range(4)]                    # [8,4,2,1]
DE_IDX = [2 ** (N - 1 - j) + 2 ** (N - 1 - k) for (j, k) in DE_PAIRS]  # [12,10,9,6,5,3]


def kron_at(op, l):
    mats = [op if k == l else I2 for k in range(N)]
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def xxz_H(J, Delta):
    """H = J*sum(XX+YY) + J*Delta*sum(ZZ) on the open 4-chain."""
    H = np.zeros((D, D), dtype=complex)
    for l in range(N - 1):
        H += J * (kron_at(X, l) @ kron_at(X, l + 1) + kron_at(Y, l) @ kron_at(Y, l + 1))
        H += J * Delta * (kron_at(Z, l) @ kron_at(Z, l + 1))
    return H


def sectors(J, Delta):
    """Return (M_SE 4x4, M_DE 6x6), Pauli-built, including the ZZ diagonal.
    For real J the blocks are real-symmetric (asserted); complex J (used only by
    complex-q exploration) is the analytic continuation -- blocks stay complex."""
    H = xxz_H(J, Delta)
    M_SE = H[np.ix_(SE_IDX, SE_IDX)]
    M_DE = H[np.ix_(DE_IDX, DE_IDX)]
    if abs(np.imag(J)) < 1e-15 and abs(np.imag(Delta)) < 1e-15:
        assert np.max(np.abs(M_SE.imag)) < 1e-12 and np.max(np.abs(M_DE.imag)) < 1e-12
        return M_SE.real.astype(complex), M_DE.real.astype(complex)
    return M_SE.astype(complex), M_DE.astype(complex)


# --------------------------------------------------------------------------
# The 24x24 (SE,DE) coherence block, Pauli-derived sectors
# --------------------------------------------------------------------------
def build_L_block(J, gamma, Delta):
    M_SE, M_DE = sectors(J, Delta)
    L = np.zeros((24, 24), dtype=complex)
    for col, (i, jk) in enumerate(BASIS):
        for i2 in range(4):                                # SE side (ket): -i*M_SE
            if M_SE[i2, i] != 0:
                L[INDEX[(i2, jk)], col] += -1j * M_SE[i2, i]
        jk_idx = DE_PAIRS.index(jk)                        # DE side (bra): +i*M_DE
        for jk2 in range(6):
            if M_DE[jk_idx, jk2] != 0:
                L[INDEX[(i, DE_PAIRS[jk2])], col] += 1j * M_DE[jk_idx, jk2]
        L[col, col] += -2 * gamma if i in jk else -6 * gamma
    return L


def build_L_reference(J, gamma):
    """Committed Delta=0 reference build (pure hops), copied verbatim from
    f89_jordan_definitive.py -- the construction must reduce to this at Delta=0."""
    M_SE = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J
    M_DE = np.zeros((6, 6))
    for col, (j, k) in enumerate(DE_PAIRS):
        for nj in (j - 1, j + 1):
            if 0 <= nj <= 3 and nj != k:
                p = tuple(sorted((nj, k)))
                if p in DE_PAIRS:
                    M_DE[DE_PAIRS.index(p), col] += 2 * J
        for nk in (k - 1, k + 1):
            if 0 <= nk <= 3 and nk != j:
                p = tuple(sorted((j, nk)))
                if p in DE_PAIRS:
                    M_DE[DE_PAIRS.index(p), col] += 2 * J
    L = np.zeros((24, 24), dtype=complex)
    for col, (i, jk) in enumerate(BASIS):
        for i2 in range(4):
            if M_SE[i2, i] != 0:
                L[INDEX[(i2, jk)], col] += -1j * M_SE[i2, i]
        jk_idx = DE_PAIRS.index(jk)
        for jk2 in range(6):
            if M_DE[jk_idx, jk2] != 0:
                L[INDEX[(i, DE_PAIRS[jk2])], col] += 1j * M_DE[jk_idx, jk2]
        L[col, col] += -2 * gamma if i in jk else -6 * gamma
    return L


def build_L_full_block(J, gamma, Delta):
    """Extract the popcount-(1,2) coherence block from the FULL 256x256 XXZ
    Liouvillian (row-major vec: vec(A rho B) = (A kron B^T) vec(rho); coherence
    |i><jk| sits at vec-index i*D + jk).  Independent of build_L_block."""
    H = xxz_H(J, Delta)
    Id = np.eye(D, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = kron_at(Z, l)
        L += gamma * (np.kron(Zl, Zl) - np.kron(Id, Id))
    vidx = [SE_IDX[a] * D + DE_IDX[b] for a in range(4) for b in range(6)]
    return L[np.ix_(vidx, vidx)]


# --------------------------------------------------------------------------
# Reflection (R = +1) symmetric 12x12 sector -- the home of the octic
# --------------------------------------------------------------------------
def sym_projector():
    refl_pair = lambda jk: tuple(sorted((PERM[jk[0]], PERM[jk[1]])))
    cols, handled = [], set()
    for n, (i, jk) in enumerate(BASIS):
        if n in handled:
            continue
        m = INDEX[(PERM[i], refl_pair(jk))]
        v = np.zeros(24, dtype=complex)
        if n == m:
            v[n] = 1.0
        else:
            v[n] = v[m] = 1.0 / np.sqrt(2)
            handled.add(m)
        cols.append(v)
        handled.add(n)
    return np.column_stack(cols)                      # 24 x 12, orthonormal


P_SYM = sym_projector()


def build_L_sym(J, gamma, Delta):
    L = build_L_block(J, gamma, Delta)
    return P_SYM.conj().T @ L @ P_SYM


# --------------------------------------------------------------------------
# Artifact-free EP-character machinery (Riesz / departure-from-normality /
# geo-vs-alg / eigenvector-merge) -- reused from review_coherence_horizon_ep.py.
# NONE of these reads a raw eig eigenvector pairing (no eig-Petermann).
# --------------------------------------------------------------------------
def riesz_projector(L, lam0, r, nq=600):
    """P = (1/2pi i) oint (zI - L)^{-1} dz on a circle of radius r about lam0."""
    n = L.shape[0]
    Id = np.eye(n, dtype=complex)
    P = np.zeros((n, n), dtype=complex)
    for k in range(nq):
        th = 2 * np.pi * (k + 0.5) / nq
        z = lam0 + r * np.exp(1j * th)
        dz = 1j * r * np.exp(1j * th) * (2 * np.pi / nq)
        P += np.linalg.solve(z * Id - L, Id) * dz
    return P / (2j * np.pi)


def jordan_counts(L, lam, tol):
    """g1 = nullity(L-lamI), g2 = nullity((L-lamI)^2) by ABSOLUTE SVD tol."""
    n = L.shape[0]
    A = L - lam * np.eye(n)
    s1 = np.linalg.svd(A, compute_uv=False)
    s2 = np.linalg.svd(A @ A, compute_uv=False)
    return int(np.sum(s1 < tol)), int(np.sum(s2 < tol))


def closest_pair_near(L, anchor, radius=0.8):
    """Among eigenvalues of L within `radius` of anchor, return the closest pair
    (dist, la, lb, all_w).  Widens to the 4 nearest if fewer than 2 in range."""
    w = np.linalg.eigvals(L)
    near = w[np.abs(w - anchor) < radius]
    if len(near) < 2:
        near = w[np.argsort(np.abs(w - anchor))[:4]]
    best = (np.inf, None, None)
    for a in range(len(near)):
        for b in range(a + 1, len(near)):
            d = abs(near[a] - near[b])
            if d < best[0]:
                best = (d, near[a], near[b])
    return best[0], best[1], best[2], w


def characterize_pair(L, la, lb, w):
    """Artifact-free character of the (possibly split) pair (la,lb) of L.
    Contour encloses EXACTLY the 2 eigenvalues nearest the midpoint; the 2x2
    compression A = V^H L V (V = orthonormal basis of range P) carries dep and
    the enclosed alg/geo. Returns a dict."""
    lam0 = 0.5 * (la + lb)
    sep = abs(la - lb)
    dist = np.sort(np.abs(w - lam0))
    third = float(dist[2]) if len(dist) > 2 else 10.0
    isolated = dist[1] < third - 1e-12
    r = 0.5 * (dist[1] + third)                       # between 2nd- and 3rd-nearest
    P = riesz_projector(L, lam0, r)
    pn = float(np.linalg.norm(P, 2))
    m_alg = float(np.trace(P).real)
    m = max(int(round(m_alg)), 1)
    U, S, _ = np.linalg.svd(P)
    V = U[:, :m]
    A = V.conj().T @ L @ V
    eigA = np.linalg.eigvals(A)
    dep = float(np.sqrt(max(0.0, np.linalg.norm(A, 'fro') ** 2 - np.sum(np.abs(eigA) ** 2))))
    A_norm = float(np.linalg.norm(A, 'fro'))
    # geometric mult of the compression at its own mean eigenvalue
    lam_bar = complex(np.mean(eigA))
    sA = np.linalg.svd(A - lam_bar * np.eye(A.shape[0]), compute_uv=False)
    geo = int(np.sum(sA < max(1e-6 * (sA.max() if sA.size else 0.0),
                              1e-7 * max(A_norm, 1.0))))
    return dict(lam0=lam0, sep=sep, third=third, isolated=isolated, r=r, pn=pn,
                m_alg=m_alg, m=m, geo=max(geo, 1), dep=dep, A_norm=A_norm,
                A=A, eigA=eigA)


# --------------------------------------------------------------------------
# GATE 0 -- validate every measure on KNOWN-ANSWER toy 2x2 matrices
# --------------------------------------------------------------------------
def gate0():
    print("=" * 78)
    print("GATE 0  validate nullity + dep + ||P|| on KNOWN answers")
    print("=" * 78)
    ok = True
    lam = -4.0 + 1.3j

    Jb = np.array([[lam, 1.0], [0.0, lam]], dtype=complex)         # DEFECTIVE
    g1, g2 = jordan_counts(Jb, lam, 1e-8)
    d = characterize_pair(Jb, lam, lam, np.linalg.eigvals(Jb))
    okJ = (g1, g2) == (1, 2) and d['dep'] > 0.5
    ok = ok and okJ
    print(f"  defective [[l,1],[0,l]] : g1={g1} g2={g2} (want 1,2)  dep={d['dep']:.3e} (want >0)"
          f"  {'PASS' if okJ else 'FAIL'}")

    Db = np.diag([lam, lam]).astype(complex)                      # DIABOLIC
    g1, g2 = jordan_counts(Db, lam, 1e-8)
    d = characterize_pair(Db, lam, lam, np.linalg.eigvals(Db))
    okD = (g1, g2) == (2, 2) and d['dep'] < 1e-6
    ok = ok and okD
    print(f"  diabolic  diag(l,l)     : g1={g1} g2={g2} (want 2,2)  dep={d['dep']:.3e} (want ~0)"
          f"  {'PASS' if okD else 'FAIL'}")

    # DEFECTIVE embedded obliquely inside a non-normal 6x6 (the F86a failure mode):
    # the artifact-free dep/geo must still read DEFECTIVE even with ||P|| > 1.
    rng = np.random.default_rng(7)
    Jblk = np.array([[lam, 1.0], [0.0, lam]], dtype=complex)
    rest = np.diag([1 + 1j, -2 + 0.5j, 3 - 1j, 0.7j])
    big = np.zeros((6, 6), dtype=complex)
    big[:2, :2] = Jblk
    big[2:, 2:] = rest
    T = np.eye(6, dtype=complex) + 0.6 * rng.standard_normal((6, 6)) \
        + 0.6j * rng.standard_normal((6, 6))
    Lob = T @ big @ np.linalg.inv(T)
    w = np.linalg.eigvals(Lob)
    d = characterize_pair(Lob, lam, lam, w)
    okO = d['dep'] > 0.1 and d['geo'] < d['m'] and d['pn'] > 1.0
    ok = ok and okO
    print(f"  oblique-embedded Jordan : ||P||={d['pn']:.3f} (>1)  dep={d['dep']:.3e} (>0)"
          f"  geo={d['geo']}<alg={d['m']}  {'PASS' if okO else 'FAIL'}")

    print("\n" + ("GATE 0 PASSED.\n" if ok else "*** GATE 0 FAILED -- STOP ***\n"))
    return ok


# --------------------------------------------------------------------------
# STAGE 0 -- reproduce the settled Delta=0 DIABOLIC verdict + build cross-checks
# --------------------------------------------------------------------------
def stage0(q_ep, gamma):
    print("=" * 78)
    print("STAGE 0  Delta=0 must reproduce the settled DIABOLIC octic EP (or STOP)")
    print("=" * 78)
    J = q_ep * gamma
    ok = True

    # (i) Pauli build reduces to the committed reference build at Delta=0
    Lp = build_L_block(J, gamma, 0.0)
    Lref = build_L_reference(J, gamma)
    d_ref = float(np.linalg.norm(Lp - Lref))
    print(f"  (i)  ||L_pauli(Delta=0) - L_reference||           = {d_ref:.2e}  "
          f"{'PASS' if d_ref < 1e-12 else 'FAIL'}")
    ok = ok and d_ref < 1e-12

    # (ii) Pauli build == popcount-(1,2) block of the full 256x256 Liouvillian
    Lfb = build_L_full_block(J, gamma, 0.0)
    d_full = float(np.linalg.norm(Lp - Lfb))
    print(f"  (ii) ||L_pauli - full-256 (1,2)-coherence block||  = {d_full:.2e}  "
          f"{'PASS' if d_full < 1e-10 else 'FAIL'}")
    ok = ok and d_full < 1e-10

    # (iii) octic: 8 non-AT-locked modes in the symmetric sector; coalescing pair
    Lsym = build_L_sym(J, gamma, 0.0)
    ev = np.linalg.eigvals(Lsym)
    rates = -ev.real / gamma
    is_at = (np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6)
    octic = list(ev[~is_at])
    n_oct = len(octic)
    pr = sorted((abs(octic[a] - octic[b]), a, b)
                for a in range(n_oct) for b in range(a + 1, n_oct))
    dmin, a0, b0 = pr[0]
    lam_ep = 0.5 * (octic[a0] + octic[b0])
    lam_exact = complex(-4 * gamma, 2 * J)
    print(f"  (iii) octic modes (sym, non-AT)                    = {n_oct}  "
          f"{'PASS' if n_oct == 8 else 'FAIL'}")
    print(f"        min coalescing pair distance                = {dmin:.3e}")
    print(f"        lam_EP numeric  = {lam_ep:.8f}")
    print(f"        lam_EP analytic = {lam_exact:.8f}  "
          f"(|dRe|={abs(lam_ep.real - lam_exact.real):.1e} |dIm|={abs(lam_ep.imag - lam_exact.imag):.1e})")
    ok = ok and n_oct == 8 and dmin < 1e-5

    # (iv) character at the analytic lam_EP: g1=g2=2, dep~0, ||P|| finite => DIABOLIC
    g1, g2 = jordan_counts(Lsym, lam_exact, 1e-5)
    la, lb, w = octic[a0], octic[b0], np.linalg.eigvals(Lsym)
    ch = characterize_pair(Lsym, la, lb, w)
    verdict = ("DIABOLIC" if (g1 == 2 and g2 == 2 and ch['dep'] < 1e-2)
               else "DEFECTIVE" if (g1 == 1 and g2 == 2) else f"? (g1={g1},g2={g2})")
    print(f"  (iv)  g1={g1} g2={g2}  dep={ch['dep']:.3e}  ||P||={ch['pn']:.3f}  "
          f"geo={ch['geo']}/alg={ch['m']}")
    print(f"        => {verdict}  {'PASS (matches settled diabolic)' if verdict == 'DIABOLIC' else 'MISMATCH'}")
    ok = ok and verdict == "DIABOLIC"

    print("\n" + ("STAGE 0 PASSED: the right object is in hand.\n" if ok
                  else "*** STAGE 0 FAILED -- construction is wrong; STOP ***\n"))
    return ok, lam_ep


# --------------------------------------------------------------------------
# THE COMPLEX-q LOCATOR -- the discriminant Newton on the tracked pair
# --------------------------------------------------------------------------
NEWTON_H = 1e-6            # central-difference step for f', f'' (f holomorphic: the real direction gives f')
NEWTON_CAP = 0.05          # one step never moves q farther than this
NEWTON_STOP = 1e-12        # a step below this * max(1,|q|): the next is already at the rounding floor of q
WINDING_R = 1e-7           # zero-order circle
WINDING_PTS = 64
WINDING_SLACK = 1e-12      # closed circle: the sum is 2*pi*integer up to ~64*4u/2pi = 5e-15 turns of rounding
WINDING_MAX_STEP = np.pi / 4   # continuity: a simple zero turns f by 0.098 per step, a double zero by 0.196
DEP_FLOOR = 1e-6           # departure rounding on a semisimple 2x2 compression is <= 4.3e-8*|A|_F (k <= 17)
SIGN_H = 1e-9              # half-width of the real-axis sign-change bracket


def pair_discriminant(q, Delta, gamma, target, deflate=None):
    """f = (lam_a - lam_b)^2 for the R-even pair nearest target; f/(q - deflate) if a zero is divided out."""
    w = np.linalg.eigvals(build_L_sym(q * gamma, gamma, Delta))
    o = np.argsort(np.abs(w - target))
    a, b = w[o[0]], w[o[1]]
    f = (a - b) ** 2
    if deflate is not None:
        f = f / (q - deflate)
    return f, 0.5 * (a + b), abs(a - b)


def locate(Delta, gamma, q, target, deflate=None, maxit=60):
    """Newton on u = f/f' from (q, target); returns (q*, pair midpoint, pair gap, iterations)."""
    q = complex(q)
    f, mid, gap = pair_discriminant(q, Delta, gamma, target, deflate)
    for it in range(maxit):
        if f == 0:
            return q, mid, gap, it
        fp = pair_discriminant(q + NEWTON_H, Delta, gamma, mid, deflate)[0]
        fm = pair_discriminant(q - NEWTON_H, Delta, gamma, mid, deflate)[0]
        d1 = (fp - fm) / (2 * NEWTON_H)
        d2 = (fp - 2 * f + fm) / NEWTON_H ** 2
        if d1 == 0:
            return q, mid, gap, it
        u = f / d1
        du = 1 - f * d2 / (d1 * d1)
        step = -u / du if du != 0 else -u
        if abs(step) > NEWTON_CAP:
            step *= NEWTON_CAP / abs(step)
        q = q + step
        f, mid, gap = pair_discriminant(q, Delta, gamma, mid, deflate)
        if abs(step) <= NEWTON_STOP * max(1.0, abs(q)):
            return q, mid, gap, it + 1
    return q, mid, gap, maxit


def zero_order(Delta, gamma, q, lam):
    """Winding number of f about q on the WINDING_R circle (unrounded). The circle closes on its opening
    value, so the sum is an integer up to rounding; nan when a step fails the continuity guard."""
    turns = 0.0
    first = pair_discriminant(q + WINDING_R, Delta, gamma, lam)[0]
    prev = first
    for k in range(1, WINDING_PTS + 1):
        if k == WINDING_PTS:
            f = first
        else:
            f = pair_discriminant(q + WINDING_R * np.exp(2j * np.pi * k / WINDING_PTS), Delta, gamma, lam)[0]
        step = np.angle(f / prev)
        if not abs(step) < WINDING_MAX_STEP:
            return float("nan")
        turns += step
        prev = f
    return turns / (2 * np.pi)


def real_sign_change(Delta, gamma, q, lam):
    """f at q -+ SIGN_H on the real axis: (f_minus, f_plus, max |Im f|/|f|). A sign change of the
    real f proves an exact real-q zero between the two points."""
    fm = pair_discriminant(q.real - SIGN_H, Delta, gamma, lam)[0]
    fp = pair_discriminant(q.real + SIGN_H, Delta, gamma, lam)[0]
    return fm.real, fp.real, max(abs(fm.imag) / abs(fm), abs(fp.imag) / abs(fp))


def certify(Delta, gamma, q, lam):
    """The certificate at a located zero: gap, isolation, alg/geo/dep, zero order, real-axis sign change."""
    L = build_L_sym(q * gamma, gamma, Delta)
    dmin, la, lb, w = closest_pair_near(L, lam)
    ch = characterize_pair(L, la, lb, w)
    order = zero_order(Delta, gamma, q, ch['lam0'])
    sign = real_sign_change(Delta, gamma, q, ch['lam0']) if abs(q.imag) < 1e-12 else None
    coincident = dmin <= STRICT_FULL_BLOCK_TOL
    character = "Uncertified"
    integer = np.isfinite(order) and abs(order - round(order)) < WINDING_SLACK
    if coincident and ch['isolated'] and ch['m'] == 2 and integer:
        if ch['geo'] == 2 and round(order) == 2:
            character = "DIABOLIC"
        elif ch['geo'] == 1 and ch['dep'] > DEP_FLOOR * max(ch['A_norm'], 1.0) and round(order) == 1:
            character = "DEFECTIVE EP2"
    return dict(q=q, lam=ch['lam0'], gap=dmin, alg=ch['m'], geo=ch['geo'], dep=ch['dep'],
                order=order, sign=sign, character=character)


def scan(q_ep, gamma):
    print("=" * 78)
    print("COMPLEX-q LOCATOR  the discriminant Newton on the tracked pair, then the certificate")
    print("=" * 78)
    deltas = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5]
    lam_ep = complex(-4 * gamma, 2 * q_ep * gamma)
    rows = []
    print(f"  {'Delta':>5} {'q*':>20} {'|Im q*|':>8} {'Re lam+4':>9} {'Im lam':>12} {'gap':>8} "
          f"{'alg':>3} {'geo':>3} {'dep':>9} {'order':>6} {'f(q*-h)':>10} {'f(q*+h)':>10}  character")
    for Delta in deltas:
        z1 = locate(Delta, gamma, q_ep, lam_ep)
        zeros = [z1]
        if Delta > 0:
            zeros.append(locate(Delta, gamma, 2 * q_ep - z1[0], lam_ep, deflate=z1[0]))
        for q, lam, _, _ in sorted(zeros, key=lambda z: z[0].real):
            c = certify(Delta, gamma, q, lam)
            c['Delta'] = Delta
            rows.append(c)
            sm, sp = (f"{c['sign'][0]:>10.2e}", f"{c['sign'][1]:>10.2e}") if c['sign'] else (f"{'-':>10}", f"{'-':>10}")
            print(f"  {Delta:>5.2f} {c['q'].real:>20.15f} {abs(c['q'].imag):>8.1e} {c['lam'].real + 4:>9.1e} "
                  f"{c['lam'].imag:>12.9f} {c['gap']:>8.1e} {c['alg']:>3} {c['geo']:>3} {c['dep']:>9.7f} "
                  f"{c['order']:>6.3f} {sm} {sp}  {c['character']}")
    print()
    return rows


# --------------------------------------------------------------------------
# MECHANISM -- the two halves of L restricted to each located zero's 2-plane
# --------------------------------------------------------------------------
def mechanism(rows, gamma):
    print("=" * 78)
    print("MECHANISM  dephasing D|2D and coherent part H|2D on the Riesz 2-plane of each zero")
    print("=" * 78)
    print("  At Delta=0 both restrict to scalars (-4g*I and 2iJ*I): the twin-scalar restriction.")
    print("  scalar departure = ||M - (tr M / 2) I||_F / ||M||_F\n")
    bas_diag = np.array([-2.0 * gamma if i in jk else -6.0 * gamma for (i, jk) in BASIS], dtype=complex)
    D_sym = P_SYM.conj().T @ np.diag(bas_diag) @ P_SYM

    def departure(M):
        m = np.trace(M) / M.shape[0]
        return float(np.linalg.norm(M - m * np.eye(M.shape[0])) / max(np.linalg.norm(M), 1e-300)), m

    print(f"  {'Delta':>5} {'q*':>16} {'D|2D departure':>15} {'tr D|2D / 2':>13} {'H|2D departure':>15}")
    for c in rows:
        L = build_L_sym(c['q'] * gamma, gamma, c['Delta'])
        w = np.linalg.eigvals(L)
        dist = np.sort(np.abs(w - c['lam']))
        P = riesz_projector(L, c['lam'], min(0.4 * dist[2], 0.5))
        U, _, _ = np.linalg.svd(P)
        V = U[:, :2]
        D2 = V.conj().T @ D_sym @ V
        H2 = V.conj().T @ (L - D_sym) @ V
        dd, dm = departure(D2)
        hd, _ = departure(H2)
        print(f"  {c['Delta']:>5.2f} {c['q'].real:>16.12f} {dd:>15.2e} {dm.real:>13.9f} {hd:>15.2e}")
    print()


def report(rows):
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    base = [r for r in rows if r['Delta'] == 0.0]
    pert = [r for r in rows if r['Delta'] > 0.0]
    print(f"  Delta=0: {len(base)} zero, {base[0]['character']} (zero order {base[0]['order']:.3f}, "
          f"gap {base[0]['gap']:.1e})")
    by_delta = {}
    for r in pert:
        by_delta.setdefault(r['Delta'], []).append(r)
    for d, rs in by_delta.items():
        kinds = ", ".join(r['character'] for r in rs)
        on_line = max(abs(r['lam'].real + 4) for r in rs)
        real_q = max(abs(r['q'].imag) for r in rs)
        signs = all(r['sign'] is not None and r['sign'][0] * r['sign'][1] < 0 for r in rs)
        print(f"  Delta={d:.2f}: {len(rs)} zeros ({kinds}); max |Im q*| {real_q:.1e}, "
              f"max |Re lam + 4| {on_line:.1e}, real-axis sign change at both: {signs}")
    all_split = all(len(rs) == 2 and all(r['character'] == "DEFECTIVE EP2" for r in rs) for rs in by_delta.values())
    if base[0]['character'] == "DIABOLIC" and all_split:
        print("\n  The Delta=0 double zero splits, at every sampled Delta > 0, into two simple zeros of the pair")
        print("  discriminant, each a certified Jordan EP2 (alg 2, geo 1) at real q on Re lam = -4.")
        print("  A diabolic unfolds into EP2s under a generic perturbation, so this response alone does not")
        print("  isolate integrability as the cause.\n")
    else:
        print("\n  *** the certified pattern did not hold at every sampled Delta; read the table ***\n")


def main():
    lines = []

    class Tee:
        def write(self, text):
            lines.append(text)
            sys.__stdout__.write(text)

        def flush(self):
            sys.__stdout__.flush()

    sys.stdout = Tee()
    try:
        gamma = 1.0
        q_ep = np.sqrt((-1 + np.sqrt(13)) / 6)
        print(f"q_EP = sqrt((-1+sqrt13)/6) = {q_ep:.10f}")
        print(f"check 3q^4+q^2-1 = {3 * q_ep ** 4 + q_ep ** 2 - 1:.2e} (should be 0)")
        print(f"gamma={gamma}, J=q, lam_EP(Delta=0) = -4g+2iJ = {complex(-4 * gamma, 2 * q_ep * gamma):.6f}\n")

        if not gate0():
            sys.exit("GATE 0 failed; aborting.")
        ok, _ = stage0(q_ep, gamma)
        if not ok:
            sys.exit("STAGE 0 failed; construction broken, aborting (do not trust the scan).")

        rows = scan(q_ep, gamma)
        mechanism(rows, gamma)
        report(rows)
    finally:
        sys.stdout = sys.__stdout__
    out = Path(__file__).resolve().parent / "results" / "f89_zz_break_gate.txt"
    out.write_text("".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
