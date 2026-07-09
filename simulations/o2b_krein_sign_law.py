"""O2b Krein-index sign law (a CONJECTURE): sign(kappa_rung) = sign(|E_rung| - |O_rung|) at every seed.

Companion verifier for experiments/F89_SEED_EXISTENCE_REDUCTION.md, section "The beta-exotic,
sharpened". It tests, PROPERLY, the class-imbalance sign law that the ledger had recorded as
UNTESTED: the candidate for the one missing ingredient of the all-N beta-exotic exclusion (O2b).

The object (same as seed_existence_nullity_check.py): the (1,2) block of the Liouvillian of an
N-site XY chain under uniform Z-dephasing, the affine pencil L(q) = diag(A) + q*C with A = -2*n_diff
(rungs -2 and -6) and C = i*K, K real symmetric, so L is complex-symmetric (L^T = L). A seed is a
real q* where two real eigenvalues coalesce into a complex-conjugate pair (a defective EP2); N-1 of
them (net (N-1)/2 conjugate pairs) are forced at odd N.

The statement. Write the bipartite sign t_i = (-1)^(a0+a1+b0) (T = diag(t), the pencil's metric,
T L T = L^dag); it splits each rung into class E (t=+1) and class O (t=-1). The gauge-invariant
-rung Krein index is
    kappa_rung := v^dag T P_rung v / ||v||^2      (P_rung the rung projector).
Unlike the transpose sum s6 = sum_{rung -6} v_i^2 (whose SIGN flips with the antilinear branch and
is therefore meaningless), kappa is a normalized Hermitian form: its sign is content. At a defective
EP2 (lambda real, a 1-dimensional eigenspace) the coalescing vector is an eigenvector of A = T.conj
(A commutes with L since T L T = L^dag, and A^2 = 1), so conj(v) = c*T*v and v^dag T v = c*(v^T v)
= 0; hence kappa_-2 + kappa_-6 = 0 on this two-rung block. So the whole sign law is the ONE statement
kappa_-2 > 0 (equivalently kappa_-6 < 0), and AT A DEFECTIVE SEED it is the gauge-invariant signed
reading of s6 != 0: kappa_-6 < 0 => kappa_-6 != 0, and the same gauge gives |kappa_-6| = |s6|/||v||^2,
so kappa_-6 != 0 <=> s6 != 0. It is NOT a standalone beta-exotic exclusion: the gauge relation holds
only at geometric multiplicity 1, so kappa_-2 > 0 presupposes the defectiveness the exclusion must
establish. The direct beta-exotic test is that defectiveness itself (a 1-dim kernel, sigma[-2]=O(1)),
which the sweep checks separately (0 semisimple points).

The prediction (exact class counts, N=3,5,7 in the ledger): rung -2 has |E|-|O| = +(N-1), rung -6
has |E|-|O| = -3(N-1)/2, so the sign law predicts kappa_-2 > 0 and kappa_-6 < 0 at every seed.

RESULT (this script): the law holds at EVERY forced defective seed at N = 5, 7, 9, both R-parities,
exhaustively (net pairs = (N-1)/2 at each), and NO semisimple (beta-exotic) transition is found --
an independent re-confirmation of the beta-exotic exclusion at N=5,7 by a route orthogonal to the
disc-multiplicity certificate. Three chain lengths are not a proof (the arc's standing rule); this
identifies the candidate mechanism and the target sign, nothing more.

METHOD (the four traps this had to clear, all logged in the doc/ledger):
  1. Anchor at a real<->complex transition, never a global min gap: a min-gap finder drifts onto
     diabolic real-real crossings (gap->0 but NOT self-orthogonal, NOT a seed). Every transition is
     located as a change in the number of complex-conjugate pairs (npair), which diabolic crossings
     do not cause.
  2. |v^T v| ~ 0 certifies nothing (an isotropic vector exists in every 2-dim complex span). The
     certificate is the eigenvalue gap AND the 1-dim-kernel signature sigma[-2]=O(1) of L - lam* I.
  3. A genuine seed is an A-eigenvector, so it MUST have v^dag T v ~ 0; a stray near-real crossing
     the extractor latched onto has v^dag T v != 0 and is rejected (not read as a beta-exotic).
  4. At N>=9 the spectrum near a seed is a dense complex thicket. Split by R-parity (reflection
     i -> N-1-i, which commutes with A, C, and -- for odd N -- with T): each sector is invariant and
     T-graded, seeds live one per sector, and the halved density lets the extractor isolate the
     near-real born pair (further boxed to |Im|<0.6 to exclude the |Im|~1..5 thicket).

The coalescing vector is taken as the smallest right singular vector of L - lam*_I (SVD, robust;
eig() vectors are ill-conditioned at a near-defective point). If the sign law survives further
scrutiny it should become a C# witness (the repo rule); this is the scout that earned that.

Run:  python simulations/o2b_krein_sign_law.py         # asserts N=5,7 (exhaustive), ~3 min
      python simulations/o2b_krein_sign_law.py 9       # also N=9 (R-parity split), ~4 min more
"""
import sys
import os
import numpy as np
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_existence_nullity_check import build


# ---------------------------------------------------------------- basis, sign, imbalance
def basis(N):
    return [(a, b) for a in combinations(range(N), 2) for b in combinations(range(N), 1)]


def bipartite_sign(N):
    """t_i = (-1)^(a0+a1+b0); +1 -> class E, -1 -> class O."""
    return np.array([(-1.0) ** (a[0] + a[1] + b[0]) for (a, b) in basis(N)])


def class_imbalance(N, A, t):
    """(|E_rung| - |O_rung|) for rungs -2 and -6."""
    out = {}
    for rung in (-2.0, -6.0):
        m = np.isclose(A, rung)
        out[rung] = int(np.sum(t[m] > 0) - np.sum(t[m] < 0))
    return out


# ---------------------------------------------------------------- R-parity sectors
def r_sectors(N):
    """Split the (1,2) block into the two R-parity sectors (reflection i -> N-1-i).

    R commutes with A, C and (for odd N) with T -- under reflection the sign exponent maps
    a0+a1+b0 -> 3(N-1) - (a0+a1+b0), and 3(N-1) is even at odd N, so the T-class is preserved. Each
    sector is thus an invariant, T-graded subspace; seeds live one per sector and the spectral
    density halves. Returns [(label, A_sec, C_sec, t_sec), ...] with A_sec the (diagonal) rung
    vector, C_sec the projected coherent block, t_sec the sector bipartite sign."""
    A, C = build(N)
    t = bipartite_sign(N)
    states = basis(N)
    idx = {s: i for i, s in enumerate(states)}

    def refl(s):
        a, b = s
        return (tuple(sorted(N - 1 - x for x in a)), (N - 1 - b[0],))

    n = len(states)
    seen, even, odd = set(), [], []
    for i, s in enumerate(states):
        if i in seen:
            continue
        j = idx[refl(s)]
        seen.update((i, j))
        if i == j:
            v = np.zeros(n); v[i] = 1.0; even.append(v)
        else:
            ve = np.zeros(n); ve[i] = ve[j] = 1 / np.sqrt(2); even.append(ve)
            vo = np.zeros(n); vo[i] = 1 / np.sqrt(2); vo[j] = -1 / np.sqrt(2); odd.append(vo)
    out = []
    for label, cols in (("R-even", even), ("R-odd", odd)):
        Q = np.array(cols).T
        A_sec = np.diag(Q.T @ (A[:, None] * Q)).real.copy()    # diagonal: rung preserved on orbits
        C_sec = Q.T @ C @ Q
        t_sec = np.rint(np.diag(Q.T @ (t[:, None] * Q))).real.copy()
        out.append((label, A_sec, C_sec, t_sec))
    return out


# ---------------------------------------------------------------- transition finder + extractor
def find_transitions(A, C, qmax, ngrid, imtol=1e-6, progress=False):
    """Real q where the number of complex-conjugate pairs changes (a real<->complex transition).

    Returns (q_lo, q_hi, lam0, born) brackets; born=True for a birth (hi side complex). This filters
    OUT diabolic real crossings, which keep both eigenvalues real and so do not change npair."""
    D = np.diag(A)
    qs = np.linspace(qmax / ngrid, qmax, ngrid)
    npair = np.empty(ngrid, dtype=int)
    for i, q in enumerate(qs):
        npair[i] = int(np.sum(np.linalg.eigvals(D + q * C).imag > imtol))
        if progress and (i + 1) % 2000 == 0:
            print(f"  ... scan {i + 1}/{ngrid}", flush=True)
    brackets = []
    for i in range(1, ngrid):
        d = npair[i] - npair[i - 1]
        if d == 0:
            continue
        born = d > 0
        q_cplx = qs[i] if born else qs[i - 1]
        cplx = np.linalg.eigvals(D + q_cplx * C)
        cplx = cplx[cplx.imag > imtol]
        cplx = cplx[np.argsort(cplx.imag)]                     # smallest Im = freshest
        for k in range(min(abs(d), len(cplx))):
            brackets.append((qs[i - 1], qs[i], float(cplx[k].real), born))
    return brackets


def analyze_transition(A, C, t, q_lo, q_hi, lam0, born, iters=90, re_w=0.6, im_w=0.6):
    """Locate the coalescence in the tight cell, SVD-extract the vector, read the Krein indices."""
    D = np.diag(A)
    I = np.eye(len(A))

    def box_pair(q):
        """Min-gap pair among near-real eigenvalues in a box about (lam0, real axis). The coalescing
        pair is near-real close to q*, while the dense thicket has |Im| ~ 1..5; the box excludes it."""
        w = np.linalg.eigvals(D + q * C)
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

    scan = np.linspace(q_lo, q_hi, 120)
    gvals = np.array([box_pair(q)[0] for q in scan])
    j = int(np.argmin(gvals))
    lo, hi = scan[max(0, j - 1)], scan[min(len(scan) - 1, j + 1)]
    gr = (np.sqrt(5) - 1) / 2
    c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    fc, fd = box_pair(c)[0], box_pair(d)[0]
    for _ in range(iters):
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - gr * (hi - lo); fc = box_pair(c)[0]
        else:
            lo, c, fc = c, d, fd
            d = lo + gr * (hi - lo); fd = box_pair(d)[0]
    q = 0.5 * (lo + hi)
    gap, la, lb = box_pair(q)
    lam_mean = 0.5 * (la + lb) if la is not None else complex(lam0)

    _, S, Vh = np.linalg.svd((D + q * C) - lam_mean * I)
    v = Vh[-1, :].conj()                                       # coalescing vector (robust vs eig)
    nrm2 = np.vdot(v, v).real
    kappa = {}
    for rung in (-2.0, -6.0):
        m = np.isclose(A, rung)
        kappa[rung] = float(np.sum(t[m] * np.abs(v[m]) ** 2) / nrm2)
    return dict(q=q, gap=gap, lam0=lam_mean.real, sigma_min=S[-1], sigma_2=S[-2],
                vTtv=np.vdot(v, t * v).real, nrm2=nrm2,
                kappa2=kappa[-2.0], kappa6=kappa[-6.0], born=born)


def classify(r):
    """defective seed / SEMISIMPLE (beta-exotic) / not-A-sym (stray crossing) / unresolved.

    The A-symmetry gate is checked before the semisimple probe, so a genuine 2-dim kernel whose SVD
    vector came out a non-A-eigen mix would be tagged 'not-A-sym' and skipped, not flagged. The real
    exhaustiveness guard is the verify() assertion net_pairs == (N-1)/2 -- a missed or spurious birth
    breaks that count; classify alone is not the beta-exotic guard."""
    if r['gap'] >= 1e-4:
        return "unresolved"
    if abs(r['vTtv']) / r['nrm2'] >= 1e-5:
        return "not-A-sym"                                     # gap small but not the A-symmetric seed
    if r['sigma_2'] < 1e-4:
        return "SEMISIMPLE!"                                    # genuine 2-dim kernel = beta-exotic
    return "defective"


# ---------------------------------------------------------------- driver
def verify(N, qmax=20.0, ngrid=8000, verbose=True):
    """Run the finder inside each R-parity sector; return the aggregate reading."""
    A, C = build(N)
    imb = class_imbalance(N, A, bipartite_sign(N))
    assert imb[-2.0] == N - 1 and imb[-6.0] == -3 * (N - 1) // 2, \
        f"class imbalance {imb} != predicted (+(N-1)={N-1}, -3(N-1)/2={-3*(N-1)//2}) at N={N}"
    exp2, exp6 = np.sign(imb[-2.0]), np.sign(imb[-6.0])
    if verbose:
        print(f"N={N}: |E|-|O|  rung-2={imb[-2.0]:+d} (=+(N-1)), rung-6={imb[-6.0]:+d} (=-3(N-1)/2)"
              f"  => predict kappa_-2>0, kappa_-6<0")
        print(f"{'sector':>7} {'q*':>10} {'lam0':>9} {'dir':>5} {'gap':>9} {'sig[-2]':>9} "
              f"{'|v+Tv|':>9} {'kappa_-2':>10} {'kappa_-6':>10} {'class':>11} {'law?':>5}")
    defective, law_ok, n_semi, n_up, n_down = 0, True, 0, 0, 0
    for label, A_s, C_s, t_s in r_sectors(N):
        brackets = find_transitions(A_s, C_s, qmax, ngrid, progress=verbose)
        seen_q = []
        for q_lo, q_hi, lam0, born in sorted(brackets):
            r = analyze_transition(A_s, C_s, t_s, q_lo, q_hi, lam0, born)
            if any(abs(r['q'] - sq) < 1e-4 for sq in seen_q):
                continue
            seen_q.append(r['q'])
            cls = classify(r)
            law = "-"
            if cls == "defective":
                defective += 1
                n_up += born
                n_down += (not born)
                ok = (np.sign(r['kappa2']) == exp2 and np.sign(r['kappa6']) == exp6)
                law = str(ok)
                law_ok &= ok
            elif cls == "SEMISIMPLE!":
                n_semi += 1
            if verbose:
                print(f"{label:>7} {r['q']:>10.6f} {r['lam0']:>9.4f} {'up' if born else 'down':>5} "
                      f"{r['gap']:>9.1e} {r['sigma_2']:>9.1e} {abs(r['vTtv']):>9.1e} "
                      f"{r['kappa2']:>+10.4f} {r['kappa6']:>+10.4f} {cls:>11} {law:>5}")
    net = n_up - n_down
    res = dict(N=N, defective=defective, net_pairs=net, forced=(N - 1) // 2,
               semisimple=n_semi, law_ok=law_ok)
    if verbose:
        print(f"  => {defective} defective ({n_up} up, {n_down} down; net {net}), "
              f"forced (N-1)/2={res['forced']} [{net == res['forced']}], "
              f"semisimple/beta-exotic={n_semi}, sign law holds={law_ok}\n")
    return res


if __name__ == "__main__":
    todo = [int(a) for a in sys.argv[1:]] or [5, 7]
    for N in todo:
        r = verify(N)
        assert r['law_ok'], f"sign law FAILED at N={N}"
        assert r['net_pairs'] == r['forced'], \
            f"net pairs {r['net_pairs']} != forced (N-1)/2={r['forced']} at N={N} (not exhaustive)"
        assert r['semisimple'] == 0, f"a semisimple (beta-exotic) transition appeared at N={N}"
    print("All assertions passed: kappa_-2>0 and kappa_-6<0 at every forced defective seed, "
          f"exhaustively, no beta-exotic, at N in {todo}.")
