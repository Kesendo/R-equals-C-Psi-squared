"""Theorem D of "The count as a plane crossing" (docs/proofs/PROOF_WEIGHT1_DEGENERACY.md): the inertia identity.

In the diagonal block (p, p) of a real number-conserving H under uniform Z-dephasing (Pauli J, units gamma = 1,
x = J/gamma) write T(x) = B(x) + 2 with B = -i x A + D, A = H (x) I - I (x) H^T real symmetric, D = -2 Hamming, tau the
transpose. tau T is Hermitian (Lemma 1). Split the cells into populations P (Hamming 0) and coherences Q (Hamming >= 2),
and Q into the tau-even Q+ and tau-odd Q- (the symmetric and antisymmetric pairs); kappa_minus = dim Q-.

Theorem D.  For every x with T(x) nonsingular (every gamma off E_p),

    n(x) = #{ eigenvalues of B with Re > -2 }  =  #{ j : r_j(x) < 1 },

the r_j the eigenvalues of the Schur matrix R of Lemma 3, the c_p zero branches included. Proof: Haynsworth on the
Hermitian tau T, In(tau T) = In(tau_Q T_QQ) + In(S) with S = 2(I - R) the Schur complement onto P (tau_P = I);
tau_Q T_QQ is nonsingular at every x (Lemma 3) so its inertia is constant and equals its Zeno value, #pos = kappa_minus;
on the other side #pos(tau T) = n_{++} + n_{--} + n_c by the sign bookkeeping of the Krein form on T's root subspaces,
kappa_minus = n_{+-} + n_{--} + n_c, and n_{+-} = 0 by Theorem A. Hence n = n_{++} = #{r_j < 1}.
Corollaries: #E_p = C(N,p) - c_p - m_p EXACTLY (Theorem C is an equality), the converse of Theorem C, n non-decreasing
in gamma (no mode returns below the plane, on every graph), and the 2 gamma regime exists iff #E_p = C(N,p) - c_p in
every block.

Rows (an exact row compares to 0 or to an integer; a float row states its error model). A grid point is "at E" exactly
when det S = 0 (T_QQ is never singular, row I2), so a point is skipped when some Schur eigenvalue s_j lies within
1e3 eps ||S|| cond(T_QQ) of 0 (an upper bound on the solve's error, loose by four to five decades against a
high-precision Schur complement, and wide at large x where ||S|| grows as x^2 with the r_j(0+) = infinity branches;
nothing is counted inside it, and the crossings all sit at x < 2); the skip is exercised by the positive control of row I3.
Counts of eigenvalues of T use the window 1e3 eps ||T|| (rounding; the smallest genuine |Re| is printed beside it), and
counts of eigenvalues of the Hermitian tau T the same window (Weyl):
  I1  tau T - (tau T)^dagger = 0 EXACTLY at every grid point, every block (the entries are the same stored numbers,
      conjugated), and tau_P = I on the populations. Control: the DM term XY - YX breaks it (ring N = 5, a flux).
  I2  the inertia of tau_Q T_QQ is constant in x and equals kappa_minus = dim Q-; the law behind the count: every
      eigenvalue of tau_Q T_QQ has |eig| = a singular value of T_QQ >= 2 (tau_Q is a permutation; the Hermitian part of
      T_QQ is 2 - 2 Hamming <= -2), gated as min |eig| >= 2 - 8 eps ||tau_Q T_QQ|| at every grid point, tight at small x.
      Control: one coherence cell's diagonal replaced by 0 must break the bound.
  I3  Haynsworth: #pos(tau T) = kappa_minus + #pos(S) at every grid point; |S - S^dagger| / (eps ||T|| cond(T_QQ)) gated
      below 1e3 (S is real symmetric, the residual is the solve's). Positive control of the skip: on the N = 3 chain (1,1)
      the sign change of det S is bisected to x = 1/sqrt 3 (the point gamma/J = sqrt 3 of E(3, chain)) and the skip fires
      there and not at x (1 +- 1e-3).
  I4  the identity n(x) = #{r_j < 1} at every grid point of every block: chain N = 3, 4, 5, ring/star/K4 at N = 4,
      ring N = 5 (the blocks with m_p > 0 included: the ring's (1,1) and (2,2), the star's (1,1), K4's (1,1) and (2,2)).
      Control: the N = 5 ring with a DM flux (H complex Hermitian, Theorem A's hypothesis gone) must break the
      identity on a fine grid over x = [0.3, 0.8], where the flux pulls a complex mode below the plane.
  I5  Theorem C as an equality: #E_p + m_p = C(N,p) - c_p in every block. #E_p is counted as the CROSSINGS, the total
      rise of #{s_j < -noise} along a sweep x = 0.02..400 (Theorem B: each branch meets 1 at most once, transversally; a
      step holding two crossings raises the count by 2; a limit-1 branch, s_j -> 0+, never enters), independent of the
      identity; beside it #{r_j > 1 at x = 400} must
      agree (no crossing beyond x = 20, so min E > 1/20 is checked, not assumed), and along the same sweep n(x) must
      be non-increasing in x (Corollary D.2, no mode returns below the plane). m_p = n(x) - c_p is read at x = 100, 200,
      400 (constant, the law); the branches with r_j(0+) = 1 are counted in m_p: their modes sit inside at every finite x
      with Re(T-eigenvalue) ~ 1/x^2 (ratio 4 between x = 200 and 400, gated in [3.8, 4.2]; 7.8e-7 and 2.7e-6 at x = 400
      against the window 1e-9 and rounding eps ||T|| ~ 1e-12, three decades either side). The chain reads
      m_p = 0 and #E_p = C(N,p) - 1 in every block; the ring at N = 4 reads m_1 = 1 (the K = pi mode) and m_2 = 3
      (one mode at limiting height 0.845 and two at limiting height exactly 1, approached from below), K4 reads
      m_1 = 3 and m_2 = 2, the star m_1 = 1 and m_2 = 0.

Run: python simulations/f50_inertia_identity.py   (about one minute; prints "ALL GATES PASS" or a FAIL line)
"""
import sys
from math import comb
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
EPS = np.finfo(float).eps
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{(': ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


# ----------------------------------------------------------------------------------------------------------
# builders (the cell basis of f50_plane_crossing_count.py, copied: a producer is never imported)
# ----------------------------------------------------------------------------------------------------------
def bonds_chain(N):
    return [(i, i + 1) for i in range(N - 1)]


def bonds_ring(N):
    return [(i, (i + 1) % N) for i in range(N)]


def bonds_star(N):
    return [(0, i) for i in range(1, N)]


def bonds_complete(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


def ham_matrix(N, bonds, dm=0.0):
    """H = sum_bonds XX + YY + ZZ + dm (XY - YX); real symmetric when dm = 0."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex if dm else float)
    bit = lambda x, l: (x >> (N - 1 - l)) & 1
    for a, b in bonds:
        for x in range(d):
            if bit(x, a) != bit(x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                H[y, x] += 2
                if dm:
                    H[y, x] += 2j * dm * (1 if bit(x, a) else -1)
                H[x, x] -= 1
            else:
                H[x, x] += 1
    return H


def block(N, H, p):
    d = 2 ** N
    pc = lambda x: bin(x).count('1')
    idx = [(x, y) for x in range(d) for y in range(d) if pc(x) == p and pc(y) == p]
    pos = {e: i for i, e in enumerate(idx)}
    n = len(idx)
    A = np.zeros((n, n), dtype=H.dtype)
    for c, (x, y) in enumerate(idx):
        for xp in range(d):
            if H[xp, x]:
                A[pos[(xp, y)], c] += H[xp, x]
        for yp in range(d):
            if H[y, yp]:
                A[pos[(x, yp)], c] -= H[y, yp]
    ham = np.array([bin(x ^ y).count('1') for x, y in idx])
    return A, ham, idx


def tau_matrix(idx):
    pos = {e: i for i, e in enumerate(idx)}
    n = len(idx)
    Tau = np.zeros((n, n))
    for i, (x, y) in enumerate(idx):
        Tau[pos[(y, x)], i] = 1.0
    return Tau


def components(N, p, bonds):
    """c_p: components of the exclusion graph on the popcount-p configurations."""
    confs = [c for c in range(2 ** N) if bin(c).count('1') == p]
    parent = {c: c for c in confs}

    def find(c):
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    for c in confs:
        for a, b in bonds:
            if ((c >> a) & 1) != ((c >> b) & 1):
                parent[find(c)] = find(c ^ (1 << a) ^ (1 << b))
    return len({find(c) for c in confs})


class Block:
    def __init__(self, N, bonds, p, dm=0.0):
        self.N, self.p = N, p
        H = ham_matrix(N, bonds, dm=dm)
        self.A, self.ham, self.idx = block(N, H, p)
        self.n = len(self.idx)
        self.P = np.where(self.ham == 0)[0]
        self.Q = np.where(self.ham > 0)[0]
        self.Tau = tau_matrix(self.idx)
        self.kappa_minus = int(np.sum(np.linalg.eigvalsh(self.Tau) < 0))
        self.C = comb(N, p)
        self.c_p = components(N, p, bonds)

    def T(self, x):
        return -1j * x * self.A + np.diag(-2.0 * self.ham) + 2.0 * np.eye(self.n)

    def schur(self, T):
        P, Q = self.P, self.Q
        TQQ = T[np.ix_(Q, Q)]
        S = T[np.ix_(P, P)] - T[np.ix_(P, Q)] @ np.linalg.solve(TQQ, T[np.ix_(Q, P)])
        return S, TQQ


def herm_eigs(M):
    return np.linalg.eigvalsh((M + M.conj().T) / 2)


def schur_eigs(blk, x):
    T = blk.T(x)
    S, TQQ = blk.schur(T)
    return herm_eigs(S), T, S, TQQ


GRID = np.concatenate([np.geomspace(0.02, 30, 48), np.random.default_rng(11).uniform(0.05, 25, 16)])

BLOCKS = [(f"chain N={N}", bonds_chain(N), N, p) for N in (3, 4, 5) for p in range(1, N // 2 + 1)]
BLOCKS += [(f"{name} N=4", b, 4, p) for name, b in (("ring", bonds_ring(4)), ("star", bonds_star(4)), ("K4", bonds_complete(4))) for p in (1, 2)]
BLOCKS += [("ring N=5", bonds_ring(5), 5, p) for p in (1, 2)]

# ----------------------------------------------------------------------------------------------------------
print("I1  tau T Hermitian exactly; tau_P = I")
for label, bonds, N, p in BLOCKS:
    blk = Block(N, bonds, p)
    worst = 0.0
    for x in GRID:
        tT = blk.Tau @ blk.T(x)
        worst = max(worst, float(np.abs(tT - tT.conj().T).max()))
    tauP = blk.Tau[np.ix_(blk.P, blk.P)]
    check(f"{label} ({p},{p})", worst == 0.0 and np.array_equal(tauP, np.eye(len(blk.P))),
          f"max |tau T - (tau T)^dagger| = {worst} over {len(GRID)} points; tau_P = I {np.array_equal(tauP, np.eye(len(blk.P)))}")
blk = Block(5, bonds_ring(5), 1, dm=0.5)
tT = blk.Tau @ blk.T(1.0)
check("control: DM flux breaks the Hermiticity (ring N=5, p=1)", np.abs(tT - tT.conj().T).max() > 1.0,
      f"max |tau T - (tau T)^dagger| = {np.abs(tT - tT.conj().T).max():.3f}")

# ----------------------------------------------------------------------------------------------------------
print("I2  inertia of tau_Q T_QQ constant = kappa_minus; |eig| = a singular value of T_QQ >= 2 (Hermitian part 2 - 2 Hamming <= -2)")
for label, bonds, N, p in BLOCKS:
    blk = Block(N, bonds, p)
    inert = set()
    min_abs, worst_tol = np.inf, 0.0
    for x in GRID:
        T = blk.T(x)
        tQ = blk.Tau[np.ix_(blk.Q, blk.Q)] @ T[np.ix_(blk.Q, blk.Q)]
        e = herm_eigs(tQ)
        inert.add(int(np.sum(e > 0)))
        min_abs = min(min_abs, float(np.abs(e).min()))
        worst_tol = max(worst_tol, 8 * EPS * np.linalg.norm(tQ, 2))
    check(f"{label} ({p},{p})", inert == {blk.kappa_minus} and min_abs >= 2 - worst_tol,
          f"#pos = {sorted(inert)} vs kappa_minus = dim Q- = {blk.kappa_minus}; min |eig| = {min_abs:.15f} >= 2 - {worst_tol:.1e}")
blk = Block(4, bonds_chain(4), 1)
T = blk.T(1.0)
tQ = blk.Tau[np.ix_(blk.Q, blk.Q)] @ T[np.ix_(blk.Q, blk.Q)]
k0 = blk.Q[0]
Tm = T.copy(); Tm[k0, k0] = 0.0                   # one coherence cell's 2 - 2 Hamming replaced by 0
tQm = blk.Tau[np.ix_(blk.Q, blk.Q)] @ Tm[np.ix_(blk.Q, blk.Q)]
check("control: one coherence diagonal set to 0 breaks |eig| >= 2 (chain N=4, p=1, x=1)",
      np.abs(herm_eigs(tQm)).min() < 2 - 8 * EPS * np.linalg.norm(tQm, 2) <= np.abs(herm_eigs(tQ)).min(),
      f"min |eig| {np.abs(herm_eigs(tQ)).min():.6f} -> {np.abs(herm_eigs(tQm)).min():.6f}")

# ----------------------------------------------------------------------------------------------------------
print("I3/I4  Haynsworth #pos(tau T) = kappa_minus + #pos(S), and n(x) = #{r_j < 1}, every grid point")


def counts_at(blk, x):
    """Returns (n, r_below, pos_tT, im_ratio, near_E, smallest |Re| of an eigenvalue of T off the line).
    near_E: some Schur eigenvalue within the solve's noise of 0, where det S = 0 and both counts are undefined.
    A complex pair with Re exactly 0 (running along the line, Theorem A's boundary case) is counted on neither side."""
    s, T, S, TQQ = schur_eigs(blk, x)
    lam = np.linalg.eigvals(T)
    e = herm_eigs(blk.Tau @ T)
    normT = np.linalg.norm(T, 2)
    noise_S = 1e3 * EPS * np.linalg.norm(S, 2) * np.linalg.cond(TQQ)
    noise_n = 1e3 * EPS * normT
    im_ratio = float(np.abs(S - S.conj().T).max() / (EPS * normT * np.linalg.cond(TQQ)))
    near_E = bool(np.any(np.abs(s) < noise_S))
    re = np.abs(lam.real)
    off_line = re[re > noise_n]
    smallest = float(off_line.min()) if off_line.size else np.inf
    return int(np.sum(lam.real > noise_n)), int(np.sum(s > noise_S)), int(np.sum(e > noise_n)), im_ratio, near_E, smallest


for label, bonds, N, p in BLOCKS:
    blk = Block(N, bonds, p)
    bad_h = bad_id = skipped = 0
    worst_im, smallest_re = 0.0, np.inf
    for x in GRID:
        n, r_below, pos_tT, im_ratio, near_E, smallest = counts_at(blk, x)
        worst_im = max(worst_im, im_ratio)
        if near_E:
            skipped += 1
            continue
        smallest_re = min(smallest_re, smallest)
        bad_h += pos_tT != blk.kappa_minus + r_below
        bad_id += n != r_below
    check(f"{label} ({p},{p}) dim {blk.n}, kappa_minus {blk.kappa_minus}", bad_h == 0 and bad_id == 0 and worst_im < 1e3,
          f"Haynsworth fails at {bad_h}, identity fails at {bad_id} of {len(GRID) - skipped} points ({skipped} at E skipped); "
          f"max |S - S^dagger| / (eps ||T|| cond) = {worst_im:.1f}; smallest |Re| off the line {smallest_re:.1e} vs window ~{1e3 * EPS * 30:.0e}")
# positive control of the skip: the N = 3 chain's (1,1) point gamma/J = sqrt 3, x = 1/sqrt 3
blk = Block(3, bonds_chain(3), 1)
lo, hi = 0.55, 0.60
plo = int(np.sum(schur_eigs(blk, lo)[0] > 0))
for _ in range(60):
    mid = (lo + hi) / 2
    if int(np.sum(schur_eigs(blk, mid)[0] > 0)) == plo: lo = mid
    else: hi = mid
root = (lo + hi) / 2
fires = counts_at(blk, root)[4]
quiet = not counts_at(blk, root * (1 + 1e-3))[4] and not counts_at(blk, root * (1 - 1e-3))[4]
check("control: the skip fires at the bisected point of E (chain N=3, (1,1)) and not at x (1 +- 1e-3)",
      fires and quiet and abs(root - 1 / np.sqrt(3)) < 1e-9, f"root x = {root:.12f} vs 1/sqrt3 = {1/np.sqrt(3):.12f}; fires {fires}, quiet beside {quiet}")
for p in (1, 2):
    blk = Block(5, bonds_ring(5), p, dm=0.5)
    fails = 0
    fine = np.linspace(0.3, 0.8, 51)     # where the flux pulls a complex mode below the plane (gate R5 of the plane-crossing gate)
    for x in fine:
        T = blk.T(x)
        S, TQQ = blk.schur(T)
        lam = np.linalg.eigvals(T)
        fails += int(np.sum(lam.real > 1e-9)) != int(np.sum(herm_eigs(S) > 1e-9))
    check(f"control: the DM flux (ring N=5, p={p}) breaks n = #{{r_j < 1}}", fails > 0, f"the identity fails at {fails} of {len(fine)} points in x = [0.3, 0.8]")

# ----------------------------------------------------------------------------------------------------------
print("I5  Theorem C as an equality: #E_p (crossings along a sweep) + m_p = C(N,p) - c_p; n non-increasing in x; the 1/x^2 law")
EXPECT_M = {("ring N=4", 1): 1, ("ring N=4", 2): 3, ("K4 N=4", 1): 3, ("K4 N=4", 2): 2, ("star N=4", 1): 1, ("star N=4", 2): 0}
SWEEP = np.geomspace(0.02, 400.0, 1500)
for label, bonds, N, p in BLOCKS:
    blk = Block(N, bonds, p)
    neg_S, n_sweep = [], []
    for x in SWEEP:
        s, T, S, TQQ = schur_eigs(blk, x)
        noise_S = 1e3 * EPS * np.linalg.norm(S, 2) * np.linalg.cond(TQQ)
        neg_S.append(int(np.sum(s < -noise_S)))      # branches with r_j > 1; a limit-1 branch (s -> 0+) never enters
        n_sweep.append(int(np.sum(np.linalg.eigvals(T).real > 1e3 * EPS * np.linalg.norm(T, 2))))
    neg_S, n_sweep = np.array(neg_S), np.array(n_sweep)
    d = np.diff(neg_S)
    crossings = int(d[d > 0].sum())                  # #{r_j > 1} RISES as x rises (r_j rises through 1): the crossings
    rises = int(-d[d < 0].sum())                     # a branch going back below 1 would be a fall of #{r_j > 1}: none allowed
    last_cross_x = SWEEP[1:][d > 0].max() if np.any(d > 0) else 0.0
    n_mono = bool(np.all(np.diff(n_sweep) <= 0))
    ms, offs = [], {}
    for x in (100.0, 200.0, 400.0):
        lam = np.linalg.eigvals(blk.T(x))
        pos = np.sort(lam.real[lam.real > 1e-9])
        ms.append(len(pos) - blk.c_p)
        offs[x] = pos
    r_above_400 = int(np.sum(schur_eigs(blk, 400.0)[0] < -1e-9))
    m_p = ms[-1]
    small400 = offs[400.0][offs[400.0] < 1e-2]
    small200 = offs[200.0][offs[200.0] < 4e-2]
    ratio_ok = len(small400) == len(small200) and all(3.8 <= r <= 4.2 for r in small200 / small400)
    expect_m = 0 if label.startswith("chain") or label.startswith("ring N=5") else EXPECT_M[(label, p)]
    check(f"{label} ({p},{p})",
          crossings + m_p == blk.C - blk.c_p and rises == 0 and crossings == r_above_400 and last_cross_x < 20 and n_mono
          and len(set(ms)) == 1 and ratio_ok and m_p == expect_m,
          f"#E_p = {crossings} crossings (= #{{r_j > 1 at x=400}}: {r_above_400}; last at x = {last_cross_x:.3f}; rises {rises}), "
          f"n non-increasing in x: {n_mono}; m_p = {m_p} (constant over x: {len(set(ms)) == 1}; expected {expect_m}), "
          f"sum {crossings + m_p} = C - c_p = {blk.C - blk.c_p}; {len(small400)} limit-1 branches, Re at x=400 {np.array2string(small400, precision=2)}"
          + (f", ratio 200/400 {np.array2string(small200 / small400, precision=3)}" if len(small400) else ""))

print("\nALL GATES PASS" if not FAILS else "\nFAILED: " + ", ".join(FAILS))
