"""The rational gcd certificate that closes the O2b nonvanishing at N = 5 and N = 7 (F89 arc).

Companion verifier for experiments/F89_SEED_EXISTENCE_REDUCTION.md, sections "Three attacks on
the sign law" (the rational form kappa_-2 = -S6/ST) and "The beta-exotic is excluded at N = 5 and
N = 7". It recomputes, from scratch at runtime and per R-parity sector (spatial reflection
i -> N-1-i on the (1,2)-block states [(a,b) for a in C(N,2) for b in C(N,1)]):

  1. the real pencil L(q) = D + u K over Z (u = iq, K real symmetric integer, T K T = -K exact
     for the bipartite sign T), and the integer polynomials in (lam, Q2 = q^2)
         chi = det(lam I - L),  S6 = tr(P_-6 adj(lam I - L)),  ST = tr(T adj(lam I - L))
     at N = 5 by EXACT integer Faddeev-LeVerrier at dim//2 + 1 integer u-nodes + exact
     interpolation over Z (a proof: deg_u <= dim and evenness make those nodes determining;
     the CRT route runs alongside and is asserted equal), at N = 7 by mod-p FLV + Newton
     interpolation in u^2 + CRT (PROOF grade since 2026-07-16: primes consumed until
     prod(p) > 2*B for the rigorous permanent/row-sum l1 bound B of l1_bounds(), making the
     symmetric-range CRT provably exact; the stability stop with 3 slack primes and the
     2 extra u-nodes under 3 FRESH primes stay as cross-checks); both with sum_i adj_ii =
     d chi/d lam exact, evenness in u structural + checked at the u = -1 node, and a 50-digit
     mpmath full-block adjugate check incl. S6_full = S6_E chi_O + S6_O chi_E;
  2. the split chi = AT * F_res (AT = the q-independent invariant-subspace strand factor, exact
     division with zero remainder; F_res monic in lam, degrees 18/17 at N = 5, 53/52 at N = 7)
     and the S6 factor shape S6 = content * prefactors * G with every prefactor sign-fixed on
     the strip lam in (-6,-2), Q2 > 0 and the core G irreducible over Q;
  3. disc_Lam(F_res) = c * Q2^v * A1 * A2^2 with max multiplicity 2 off Q2 = 0 (the beta-exotic
     needs 3), layer degrees matching the certified table (N=5: A1 28/28, A2 16/13, v 77/69;
     N=7: A1 114/111, A2 210/195, v 772/753; E/O order, degrees in Q2), and A1 irreducible
     over Q; at N = 5 the discriminant is additionally reconstructed EXACTLY over Z
     (integer-node interpolation with exact resultants, cross-checked coefficient-wise
     against every certificate prime), its layers split exactly, A1 and A2 both factored
     exactly (each irreducible over Q), and the real-positive root INVENTORY asserted:
     A1's real positive roots are exactly the four seed loci w = q*^2 (two per sector,
     {1.541958, 4.645014} E / {1.653988, 31.469594} O) and A2 has none (R-even) resp.
     exactly one (R-odd, w ~ 5.100831, verified to carry a single REAL double lambda-root
     of F_res: the diabolic class, not a coincident EP2 pair);
  3b. COINCIDENT-EP2 EXCLUSION (both N, both sectors): with psc_1 the first principal
     subresultant coefficient of (F_res, dF_res/dlam) in lam, a w0 carrying TWO double
     lambda-roots (a coincident EP2 pair: an order-2 COUNT-DROPPING locus on A2 that the
     gcd certificate alone cannot see) or one triple root forces deg gcd(F_res, dF') >= 2,
     i.e. psc_0(w0) = psc_1(w0) = 0 (psc_0 = the discriminant). The asserted certificate
     gcd(disc/Q2^v, psc_1/Q2^v') = 1 therefore proves NO w0 != 0, real or complex, carries
     more than one double lambda-root of F_res; in particular no real positive A2 root is a
     coincident EP2 pair, and the one repeated root at every real disc point is REAL (a
     complex double root would bring its conjugate as a second one). This is stronger than
     the per-root A2 inventory the doc claim needs, and it needs no reconstruction of A2;
  4. THE CERTIFICATE:  gcd(Res_Lam(F_res, S6), A1) = 1  over Q, asserted at >= 3 good primes
     per sector, and the cross-sector gap gcd(Res_Lam(F_res_s, chi_other), A1_s) = 1 in both
     directions (so no OTHER-sector eigenvalue meets a seed locus: the kernel at a simple-layer
     collision is 1-dimensional on the FULL block, geometric multiplicity 1);
  5. seed spot-checks: each forced seed Newton-refined (exact Jacobian, dps 80) as a double
     root of F_res, then S6 != 0, S6*ST < 0 and kappa_-2 = -S6/ST > 0 asserted against the
     eigenvector-route reference values to 1e-5 (printed at 6 decimals).

Together: on the simple layer A1 of the discriminant (which carries every count-dropping locus,
by the disc-multiplicity table of the beta-exclusion section) S6 does not vanish, and no
cross-sector coincidence can raise the kernel dimension there. Hence at every count-drop the
geometric multiplicity is 1 full-block, adj(lam* I - L) = c*r r^T with c != 0, and
s6 = S6/c != 0 at every forced seed: the O2b nonvanishing, per N, at N = 5 and N = 7.

Proof status of each layer (read this before citing):
  * PROVED OVER Q (exact integer arithmetic), N = 5 end to end: the base polynomials
    chi/S6/ST themselves (exact integer Faddeev-LeVerrier + exact interpolation over Z, the
    determining-node argument of item 1), T K T = -K, tr adj = d chi/d lam, chi = AT * F_res,
    the S6 factor split (exact divisions), the strand-prefactor signs on the strip (each
    strand prefactor is written as a sum of nonnegative even-power monomials in (lam - r0)
    and Q2, the norm-form shape, so it is strictly positive on the OPEN strip; the linear
    factors (lam+2) < 0 < (lam+6) by inspection), and the whole N = 5 disc inventory of
    item 3 (exact D, A1, A2 over Z with A1 mod p asserted to BE the mult-1 Yun layer at
    every certificate prime).
  * BASE POLYNOMIALS AT N = 7: PROOF GRADE since 2026-07-16. chi/S6/S2/ST are reconstructed
    by CRT with the prime pool consumed until prod(p) > 2*B, B the rigorous a-priori
    permanent/row-sum l1 bound of l1_bounds() (validated from below at N = 5 against the
    exact-Z coefficients). The symmetric-range CRT is then provably exact -- no stability
    assumption left; the stability stop (+3 slack primes), fresh-prime extra-node checks and
    the 50-digit full-block adjugate check remain as independent cross-checks. Everything
    downstream at N = 7 (F_res, the disc layers, the psc_1 leg, both gcd certificates, the
    seed checks) now rests on proof-grade base polynomials. (Historical: before 2026-07-16
    this was stability-stopped only, verification grade -- the first of the two named N = 7
    premises in F89_BETA_EXOTIC_GENERICITY.md's bookkeeping; it is DISCHARGED. The layer
    identification premise below remains.)
  * PROVED OVER Q VIA A GOOD PRIME, given the base polynomials (a nontrivial Q-gcd of
    integer polynomials survives reduction mod any prime that preserves both degrees, by
    Gauss's lemma; F_res and chi are monic in lam, so Res_Lam commutes with reduction
    entrywise on the Sylvester matrix, and psc_1, a fixed determinant in F_res's integer
    coefficients, commutes unconditionally): both gcd certificates of item 4, and the
    coincident-EP2 exclusion of item 3b (which uses only D itself, whose true Q2-degree the
    DISCMULT engine certified by Hadamard-bounded prime sampling, never the layer
    identification). Asserted at >= 3 primes although ONE good prime already proves each.
    At N = 5 these are therefore unconditional; at N = 7 they now rest on the PROOF-grade
    base polynomials above (since 2026-07-16), so the psc_1 leg, being also free of the
    layer premise below, is unconditional at N = 7.
  * CERTIFIED MODULO A CONDITIONAL PREMISE (N = 7 only), ON TOP of the base grade: the
    identification "mod-p mult-1 Yun layer = A1 mod p (up to scalar)". A bad prime can merge
    two roots of A1 into a square that Yun then books into the mult-2 layer: such a merge
    preserves the total degree, the valuation and max-mult <= 2 ((A1, A2) = (116, 209) over
    Q with a d = 1 merge would imitate the table's (114, 210) exactly), and is invisible to
    the DISCMULT logic as well. Only finitely many primes are bad (those dividing
    disc(A1) * Res(A1, A2)), so the agreement of several independent primes is STRONG
    EVIDENCE, not a Q-proof. Consequently the N = 7 A1-irreducibility (subset-sum
    certificate: possible Q-factor degrees lie in every good prime's mod-p degree
    partition's subset sums, intersection {0, n}) and both N = 7 gcd conclusions of item 4
    rest on this layer-identification premise alone (since 2026-07-16 the base polynomials
    underneath are PROOF grade, so it is the single remaining N = 7 premise). At N = 5
    no premise remains: A1 and A2 are reconstructed exactly over Z and the mod-p layers are
    asserted to BE their reductions. The item-3b coincident-EP2 exclusion does not use the
    layer premise at either N.
  * MULTI-PRIME EVIDENCE ONLY (not asserted at N = 7): the A2^2-split of the resultant,
    Res_Lam(F_res, S6) = c * Q2^vR * A2^2 * B. At N = 5 the split was proved EXACTLY by the
    in-session sympy scouts, so deg gcd(Res, disc) = 2*deg A2 + v IS asserted there; at N = 7
    the same gcd degree is printed as evidence, never asserted (mod-p divisibility at finitely
    many primes does not lift).
  * Honest caveat (kept from the review): the final step s6 = S6/c != 0 uses the rank-1
    adjugate identity adj = c*r r^T in the antilinear gauge, standard linear algebra at
    geometric multiplicity 1; the certificate items above are what establish that geometric
    multiplicity 1 on the simple layer.

Axis: unit-hop (this file's q = twice the octic census q*; Q2 = q^2). All reference seed values
are unit-hop, from simulations/o2b_krein_sign_law.py runs.

Run:  python simulations/o2b_gcd_certificate.py         # N = 5 (~5 s measured, dev machine)
      python simulations/o2b_gcd_certificate.py 7       # N = 5 and N = 7 (~3 min measured)
"""
import os
import sys
import time
from fractions import Fraction

import numpy as np
import sympy as sp
from mpmath import mp, mpc, mpf, matrix as mpmat, det as mpdet

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_existence_nullity_check import build            # noqa: E402
from o2b_krein_sign_law import basis, bipartite_sign      # noqa: E402

lam, Q2 = sp.symbols("lam Q2")

# ------------------------------------------------------------------ certified reference tables
SECTOR_DIM = {5: {"E": 26, "O": 24}, 7: {"E": 75, "O": 72}}
FRES_DEG = {5: {"E": 18, "O": 17}, 7: {"E": 53, "O": 52}}
# disc_Lam(F_res) layers in Q2: simple layer A1, double layer A2, valuation v at Q2 = 0
DISC_TABLE = {5: {"E": dict(A1=28, A2=16, v=77), "O": dict(A1=28, A2=13, v=69)},
              7: {"E": dict(A1=114, A2=210, v=772), "O": dict(A1=111, A2=195, v=753)}}
# N = 5 only: the EXACT split gcd(Res(F,S6), disc) = Q2^v * A2^2 (proved by the sympy scouts)
GCD_R1_DISC_N5 = {"E": 109, "O": 95}          # = v + 2*deg A2
# forced seeds (sector, q*, lam*, kappa_-2 reference), unit-hop, from o2b_krein_sign_law.py
SEEDS = {
    5: [("E", 1.241756, -4.6189, 0.314463),
        ("O", 1.286074, -3.8196, 0.140184),
        ("E", 2.155229, -3.7917, 0.351821),
        ("O", 5.609777, -4.4882, 0.062939)],
    7: [("E", 1.077392, -3.9816, 0.090208),
        ("O", 1.200758, -4.9228, 0.238066),
        ("O", 1.447107, -4.8609, 0.184075),
        ("E", 1.848215, -4.6615, 0.060903),
        ("E", 3.029667, -4.8846, 0.025288),
        ("O", 3.181705, -3.8998, 0.377212),
        ("E", 11.355150, -4.8553, 0.050276)],
}
MIN_CERT_PRIMES = 3        # gcd certificates asserted at >= this many good primes


def make_primes(count, below=2 ** 25):
    """Primes just below 2^25: int64-safe in dim<=75 matmuls (75 * (2^25)^2 < 2^63)."""
    out, p = [], below
    while len(out) < count:
        p = int(sp.prevprime(p))
        out.append(p)
    return out


# ================================================================== 1. integer sector matrices
def integer_sectors(N):
    """R-parity sectors of the (1,2) block in an INTEGER (similarity) basis.

    Columns: fixed states e_i (R-even only); orbit pairs e_i +/- e_j. The inverse rows carry
    the 1/|orbit| weight, so K_s = Q^-1 K Q is asserted integer; D and t are constant on
    orbits, hence stay diagonal in the sector basis (rung projectors and T survive the split).
    """
    A, C = build(N)
    K = np.rint((C / 1j).real).astype(np.int64)
    assert np.abs(C - 1j * K).max() < 1e-12, "C = iK with integer K"
    t = bipartite_sign(N).astype(np.int64)
    # T-antisymmetry on the full block: T K T = -K exactly (hop flips the class parity)
    assert np.array_equal(t[:, None] * K * t[None, :], -K), "T K T != -K"
    states = basis(N)
    idx = {s: i for i, s in enumerate(states)}

    def refl(s):
        a, b = s
        return (tuple(sorted(N - 1 - x for x in a)), (N - 1 - b[0],))

    seen, orbits = set(), []
    for i, s in enumerate(states):
        if i in seen:
            continue
        j = idx[refl(s)]
        seen.update((i, j))
        orbits.append((i,) if i == j else (i, j))

    out = {}
    for label, sgn in (("E", +1), ("O", -1)):
        cols, reps = [], []
        for orb in orbits:
            if len(orb) == 1:
                if sgn == +1:
                    cols.append([(orb[0], 1)])
                    reps.append(orb[0])
            else:
                cols.append([(orb[0], 1), (orb[1], sgn)])
                reps.append(orb[0])
        dim = len(cols)
        Ks = np.zeros((dim, dim), dtype=np.int64)
        for r, orb_r in enumerate(cols):
            for c, orb_c in enumerate(cols):
                val = Fraction(0)
                for (ir, cr) in orb_r:
                    for (ic, cc) in orb_c:
                        val += Fraction(cr * cc * int(K[ir, ic]), len(orb_r))
                assert val.denominator == 1, f"non-integer K_s entry at {label} ({r},{c})"
                Ks[r, c] = int(val)
        Ds = np.array([int(round(A[reps[c]])) for c in range(dim)], dtype=np.int64)
        ts = np.array([t[reps[c]] for c in range(dim)], dtype=np.int64)
        for orb in orbits:
            assert len({int(round(A[i])) for i in orb}) == 1 and len({int(t[i]) for i in orb}) == 1
        out[label] = dict(K=Ks, D=Ds, t=ts, dim=dim)
    return out


# ================================================================== 2. exact chi / S6 / ST
def fl_modp(M, p):
    """Faddeev-LeVerrier mod p: chi coefficients (ascending, monic) and diag(B_j) with
    adj(lam I - M) = sum_j B_j lam^j. M int64 (n, n), entries reduced into [0, p)."""
    n = M.shape[0]
    c = np.zeros(n + 1, dtype=np.int64)
    c[n] = 1
    diagB = np.zeros((n, n), dtype=np.int64)
    B = np.eye(n, dtype=np.int64)
    diagB[n - 1] = 1
    for j in range(n - 1, 0, -1):
        MB = (M @ B) % p
        cj = (-int(np.trace(MB) % p) * pow(n - j, p - 2, p)) % p
        c[j] = cj
        B = MB
        B[np.diag_indices(n)] = (B[np.diag_indices(n)] + cj) % p
        diagB[j - 1] = B[np.diag_indices(n)]
    MB = (M @ B) % p
    c[0] = (-int(np.trace(MB) % p) * pow(n, p - 2, p)) % p
    return c, diagB


def _fl_selftest():
    """FL vs the sympy adjugate at dim 6 (gate-first: the engine is checked before use)."""
    rng = np.random.default_rng(7)
    p = 33554393
    n = 6
    A = rng.integers(-4, 5, size=(n, n))
    x = sp.Symbol("x")
    Ms = sp.Matrix(A.tolist())
    adj = (x * sp.eye(n) - Ms).adjugate()
    chi = (x * sp.eye(n) - Ms).det()
    c, diagB = fl_modp(A % p, p)
    cpoly = sp.Poly(chi, x).all_coeffs()[::-1]
    assert all(int(cpoly[j]) % p == c[j] % p for j in range(n + 1)), "FL chi self-test"
    for i in range(n):
        pi = sp.Poly(adj[i, i], x).all_coeffs()[::-1]
        pi += [0] * (n - len(pi))
        assert all(int(pi[j]) % p == diagB[j][i] % p for j in range(n)), "FL adj self-test"


def newton_interp_modp(wnodes, samples, p):
    """Newton interpolation mod p: samples (npts, m) at wnodes -> ascending monomial coeffs."""
    npts, m = samples.shape
    coef = samples.astype(np.int64).copy() % p
    wn = np.array(wnodes, dtype=np.int64)
    for j in range(1, npts):
        d = wn[j:] - wn[:-j]
        lut = {int(x): pow(int(x) % p, p - 2, p) for x in np.unique(d)}
        dinv = np.array([lut[int(x)] for x in d], dtype=np.int64)
        coef[j:] = ((coef[j:] - coef[j - 1:-1]) % p) * dinv[:, None] % p
    out = np.zeros((npts, m), dtype=np.int64)
    acc = np.zeros(npts, dtype=np.int64)
    acc[0] = 1
    deg_acc = 0
    for j in range(npts):
        out[:deg_acc + 1] = (out[:deg_acc + 1] + coef[j] * acc[:deg_acc + 1, None]) % p
        if j == npts - 1:
            break
        new = np.zeros(npts, dtype=np.int64)
        new[1:deg_acc + 2] = acc[:deg_acc + 1]
        new[:deg_acc + 1] = (new[:deg_acc + 1] - wn[j] * acc[:deg_acc + 1]) % p
        acc = new % p
        deg_acc += 1
    return out % p


def fl_exact_node(Kobj, Dvec, u):
    """Faddeev-LeVerrier over Z (exact Python ints): chi coefficients (ascending, monic) and
    diag(B_j) of adj(lam I - M) for M = u*K + diag(D). The divisions by n-j are exact over Z
    (a property of the FLV recurrence; asserted)."""
    n = len(Dvec)
    M = (u * Kobj).copy()
    M[np.diag_indices(n)] += Dvec
    c = [0] * (n + 1)
    c[n] = 1
    diagB = np.zeros((n, n), dtype=object)
    B = np.array([[1 if i == j else 0 for j in range(n)] for i in range(n)], dtype=object)
    diagB[n - 1] = 1
    for j in range(n - 1, 0, -1):
        MB = M @ B
        tr = int(np.trace(MB))
        assert (-tr) % (n - j) == 0, "FLV division not exact over Z"
        cj = -tr // (n - j)
        c[j] = cj
        B = MB
        B[np.diag_indices(n)] += cj
        diagB[j - 1] = B[np.diag_indices(n)]
    tr = int(np.trace(M @ B))
    assert (-tr) % n == 0, "FLV division not exact over Z"
    c[0] = -tr // n
    return c, diagB


def newton_interp_exact(wnodes, samples):
    """Exact Newton interpolation over Q at integer nodes; samples: object array (npts, m) of
    exact values. Returns an integer object array (npts, m) of ascending w-monomial
    coefficients (integrality asserted: the target polynomials live in Z[w])."""
    npts, m = samples.shape
    coef = [[Fraction(int(samples[i, c])) for c in range(m)] for i in range(npts)]
    for j in range(1, npts):
        for i in range(npts - 1, j - 1, -1):
            d = wnodes[i] - wnodes[i - j]
            for cix in range(m):
                coef[i][cix] = (coef[i][cix] - coef[i - 1][cix]) / d
    out = [[Fraction(0)] * m for _ in range(npts)]
    acc = [Fraction(0)] * npts
    acc[0] = Fraction(1)
    dega = 0
    for j in range(npts):
        for k in range(dega + 1):
            if acc[k]:
                for cix in range(m):
                    out[k][cix] += coef[j][cix] * acc[k]
        if j == npts - 1:
            break
        new = [Fraction(0)] * npts
        for k in range(dega + 1):
            new[k + 1] += acc[k]
            new[k] -= wnodes[j] * acc[k]
        acc = new
        dega += 1
    assert all(v.denominator == 1 for row in out for v in row), "interpolation not integral"
    return np.array([[int(v) for v in row] for row in out], dtype=object)


def crt_stack(residues, primes):
    """Symmetric-range CRT over a list of same-shape int64 arrays, one per prime."""
    M = 1
    x = np.zeros(residues[0].shape, dtype=object)
    for r, p in zip(residues, primes):
        Minv = pow(M % p, p - 2, p)
        fx, fr = x.reshape(-1), r.reshape(-1)
        for i in range(fx.size):
            fx[i] = int(fx[i]) + M * (((int(fr[i]) - int(fx[i])) % p) * Minv % p)
        M *= p
    half = M // 2
    f = x.reshape(-1)
    for i in range(f.size):
        if f[i] > half:
            f[i] -= M
    return x, M


def l1_bounds(sec):
    """Rigorous a-priori l1-coefficient bounds for chi, S6, S2, ST over Z[lam, u].

    chi = det(lam I - D - uK) is a signed sum over permutations of entry products. The
    l1 norm of coefficients (in the two variables lam, u) is submultiplicative on products
    and subadditive on sums, so ||chi||_1 <= perm(Nmat), where Nmat_ab = ||entry_ab||_1:
    the diagonal entry lam - D_a - u K_aa has l1 = 1 + |D_a| + |K_aa|, the off-diagonal
    -u K_ab has l1 = |K_ab|. The permanent of a nonnegative matrix is bounded by the
    product of its row sums, so ||chi||_1 <= prod_a rowsum_a with
    rowsum_a = 1 + |D_a| + sum_b |K_ab|. Each adj diagonal entry adj_aa is the (a,a)
    minor, so ||adj_aa||_1 <= prod_{a' != a} rowsum_{a'} (deleting a column only lowers
    row sums). Hence ||S6||_1 <= sum_{a on the -6 rung} prod_{a' != a} rowsum_{a'}, and
    likewise S2 (the -2 rung) and ST (all rows, |t_a| = 1). The u -> w = u^2 collection
    merges no terms (all polys are even in u, asserted structurally by T K T = -K and
    pinned by the u = -1 extra node), so the same numbers bound the (lam, w) coefficients,
    and every single coefficient obeys |c| <= ||.||_1. These are exact Python ints."""
    K, D = sec["K"], sec["D"]
    dim = sec["dim"]
    rowsums = [1 + abs(int(D[a])) + int(np.abs(K[a]).sum()) for a in range(dim)]
    prod_all = 1
    for s in rowsums:
        prod_all *= s
    minor = [prod_all // rowsums[a] for a in range(dim)]
    b = dict(chi=prod_all,
             S6=sum(minor[a] for a in range(dim) if int(D[a]) == -6),
             S2=sum(minor[a] for a in range(dim) if int(D[a]) == -2),
             ST=sum(minor))
    return b


def compute_exact_polys(label, sec, work_primes, check_primes, exact_over_Z=False):
    """chi, S6, S2, ST of the sector pencil D + uK as dicts {(j_desc, dw): int} in
    (lam, w = u^2).

    With exact_over_Z (N = 5): PROOF grade. deg_u(chi) <= dim and deg_u(adj entries) <= dim-1
    (every entry of lam I - D - uK is affine in u), and evenness in u makes deg_w <= dim//2,
    so the nw = dim//2 + 1 integer w-nodes u = 0..dim//2 DETERMINE the polynomials: exact
    integer FLV at those nodes + exact interpolation over Z is a proof, no CRT involved. The
    CRT route still runs and is asserted EQUAL, as a cross-check of the shared machinery.

    Without (N = 7): PROOF grade since 2026-07-16 -- the CRT prime pool is consumed until
    prod(used primes) > 2 * B for the rigorous a-priori l1 bound B of l1_bounds() (the
    permanent/row-sum bound), which makes the symmetric-range CRT reconstruction provably
    exact; the stability stop (+3 slack primes) and the fresh-prime extra u-nodes (+nw and
    -1, the latter pinning the odd-u part to zero) are kept as independent cross-checks.
    (Before 2026-07-16 this path was verification grade: stability-stopped only, no bound.)

    Both paths assert tr adj = d chi/d lam exactly."""
    dim, K, D, t = sec["dim"], sec["K"], sec["D"], sec["t"]
    nw = dim // 2 + 1
    unodes = list(range(nw))
    extra_u = [nw, -1]
    m6 = np.where(D == -6)[0]
    m2 = np.where(D == -2)[0]
    t0 = time.time()

    def samples_for(p):
        chi_s = np.zeros((nw, dim + 1), dtype=np.int64)
        s6_s = np.zeros((nw, dim), dtype=np.int64)
        s2_s = np.zeros((nw, dim), dtype=np.int64)
        st_s = np.zeros((nw, dim), dtype=np.int64)
        for k, u in enumerate(unodes):
            M = (u * K) % p
            M[np.diag_indices(dim)] = (M[np.diag_indices(dim)] + D) % p
            c, diagB = fl_modp(M, p)
            chi_s[k] = c
            s6_s[k] = diagB[:, m6].sum(axis=1) % p
            s2_s[k] = diagB[:, m2].sum(axis=1) % p
            st_s[k] = (diagB @ t) % p
        wnodes = [u * u for u in unodes]
        return tuple(newton_interp_modp(wnodes, arr, p)
                     for arr in (chi_s, s6_s, s2_s, st_s))

    bounds = l1_bounds(sec)
    B_max = max(bounds.values())

    res = {nm: [] for nm in ("chi", "S6", "S2", "ST")}
    used, exact = [], {}
    prodp = 1
    for p in work_primes:
        for nm, arr in zip(("chi", "S6", "S2", "ST"), samples_for(p)):
            res[nm].append(arr)
        used.append(p)
        prodp *= p
        # PROOF condition: symmetric-range CRT is exact once prod(primes) > 2*B (a-priori
        # l1 bound). The stability stop stays as an independent cross-check on top.
        if prodp > 2 * B_max and len(used) >= 8:
            ok = True
            for nm in ("chi", "S6", "S2", "ST"):
                x_all, _ = crt_stack(res[nm], used)
                x_less, _ = crt_stack(res[nm][:-3], used[:-3])
                if not np.array_equal(x_all, x_less):
                    ok = False
                    break
            if ok:
                for nm in ("chi", "S6", "S2", "ST"):
                    exact[nm], _ = crt_stack(res[nm], used)
                break
    assert exact, f"CRT did not close with {len(work_primes)} primes in sector {label}"
    assert prodp > 2 * B_max, "proof condition violated (unreachable: enforced in loop)"
    print(f"  [{label}] dim={dim}: CRT at {len(used)} primes; PROOF condition "
          f"prod(p) = 2^{prodp.bit_length() - 1} > 2*B = 2^{(2 * B_max).bit_length() - 1} "
          f"(l1 bounds: chi 2^{bounds['chi'].bit_length() - 1}, S6 2^{bounds['S6'].bit_length() - 1}, "
          f"S2 2^{bounds['S2'].bit_length() - 1}, ST 2^{bounds['ST'].bit_length() - 1}); "
          f"stability +3 slack also asserted  [{time.time() - t0:.1f}s]")

    if exact_over_Z:                              # N = 5: exact integer FLV, proof grade
        t1 = time.time()
        Kobj = K.astype(object)
        Dobj = D.astype(object)
        tobj = t.astype(object)
        raw = {nm: np.zeros((nw, cols), dtype=object)
               for nm, cols in (("chi", dim + 1), ("S6", dim), ("S2", dim), ("ST", dim))}
        for k, u in enumerate(unodes):
            c, diagB = fl_exact_node(Kobj, Dobj, u)
            raw["chi"][k] = c
            raw["S6"][k] = diagB[:, m6].sum(axis=1)
            raw["S2"][k] = diagB[:, m2].sum(axis=1)
            raw["ST"][k] = diagB @ tobj
        wnodes = [u * u for u in unodes]
        for nm in ("chi", "S6", "S2", "ST"):
            exZ = newton_interp_exact(wnodes, raw[nm])
            assert np.array_equal(exZ, exact[nm]), \
                f"exact-Z FLV != CRT result for {nm} in sector {label}"
            # from-below validation of the a-priori l1 bound against the exact coefficients
            mx = max(abs(int(v)) for v in exZ.ravel())
            assert mx <= bounds[nm], f"l1 bound violated by exact {nm} in sector {label}"
            exact[nm] = exZ
        print(f"  [{label}] exact integer FLV over Z at {nw} w-nodes: equals the CRT result "
              f"(N=5 base polynomials PROOF grade); exact max|coeff| <= l1 bound validated "
              f"for all four polys  [{time.time() - t1:.1f}s]")

    for p in check_primes:                        # fresh primes, extra u-nodes
        for u in extra_u:
            M = (u * K) % p
            M[np.diag_indices(dim)] = (M[np.diag_indices(dim)] + D) % p
            c, diagB = fl_modp(M, p)
            w = (u * u) % p
            wpow = np.array([pow(w, d, p) for d in range(nw)], dtype=object)
            for nm, direct in (("chi", c), ("S6", diagB[:, m6].sum(axis=1) % p),
                               ("S2", diagB[:, m2].sum(axis=1) % p),
                               ("ST", (diagB @ t) % p)):
                pred = (exact[nm].astype(object).T % p) @ wpow % p
                assert np.array_equal(pred % p, direct % p), \
                    f"extra-node check failed {label} {nm} u={u} p={p}"
    print(f"  [{label}] fresh-prime extra-node verification (u={extra_u}) OK")

    def to_dict(arr, deg_lam):
        d = {}
        for dw in range(arr.shape[0]):
            for jl in range(arr.shape[1]):
                v = int(arr[dw, jl])
                if v:
                    d[(deg_lam - jl, dw)] = v
        return d

    chi = to_dict(exact["chi"], dim)
    S6 = to_dict(exact["S6"], dim - 1)
    S2 = to_dict(exact["S2"], dim - 1)
    ST = to_dict(exact["ST"], dim - 1)
    # sum_i adj_ii = d chi / d lam exactly (rungs -2 and -6 exhaust the block)
    tradj = {}
    for d in (S6, S2):
        for k, v in d.items():
            tradj[k] = tradj.get(k, 0) + v
    tradj = {k: v for k, v in tradj.items() if v}
    dchi = {(j, dw): c * (dim - j) for (j, dw), c in chi.items() if dim - j > 0}
    dchi = {k: v for k, v in dchi.items() if v}
    assert tradj == dchi, f"tr adj != d chi/d lam in sector {label}"
    return dict(chi=chi, S6=S6, ST=ST, dim=dim)


def to_poly(d, deg_lam):
    """Stage dict (j descending in lam, dw in w = u^2) -> sympy Poly in (lam, Q2), u^2 -> -Q2."""
    expr = 0
    for (j, dw), c in d.items():
        expr += c * (-1) ** dw * lam ** (deg_lam - j) * Q2 ** dw
    return sp.Poly(expr, lam, Q2, domain="ZZ")


def full_block_50digit_check(N, polys, npts):
    """50-digit mpmath full-block adjugate check at rational points, incl. the composition
    S6_full = S6_E chi_O + S6_O chi_E (and the same for ST)."""
    mp.dps = 50
    A, C = build(N)
    t_full = bipartite_sign(N)
    m6_full = np.isclose(A, -6.0)
    n = len(A)
    K_full = np.rint((C / 1j).real).astype(np.int64)

    def ev(P, lam_f, w_f):
        tot = Fraction(0)
        for (dl, dw), c in zip(P.monoms(), P.coeffs()):
            tot += int(c) * lam_f ** dl * w_f ** dw
        return tot

    pts = [(Fraction(-43, 10), Fraction(31, 10)), (Fraction(-57, 10), Fraction(1, 2)),
           (Fraction(-25, 8), Fraction(9, 2))][:npts]
    for lam_f, q_f in pts:
        w = q_f * q_f                              # Q2 = q^2
        lam_mp = mpf(lam_f.numerator) / mpf(lam_f.denominator)
        q_mp = mpf(q_f.numerator) / mpf(q_f.denominator)
        Mm = mpmat(n, n)
        for i in range(n):
            for j in range(n):
                re = (lam_mp - mpf(float(A[i]))) if i == j else mpf(0)
                Mm[i, j] = mpc(re, -q_mp * int(K_full[i, j]))
        d = mpdet(Mm)
        inv = Mm ** -1
        S6_num = sum((d * inv[i, i]).real for i in range(n) if m6_full[i])
        ST_num = sum(t_full[i] * (d * inv[i, i]).real for i in range(n))
        for name, key in (("S6", "S6"), ("ST", "ST")):
            pv = (ev(polys["E"][key], lam_f, w) * ev(polys["O"]["chi"], lam_f, w)
                  + ev(polys["O"][key], lam_f, w) * ev(polys["E"]["chi"], lam_f, w))
            nv = S6_num if name == "S6" else ST_num
            pv_mp = mpf(pv.numerator) / mpf(pv.denominator)
            rel = abs(pv_mp - nv) / max(mpf(1), abs(pv_mp))
            assert rel < mpf("1e-30"), f"{name} full-block mismatch at {lam_f},{q_f}"
    print(f"  50-digit full-block adjugate check: {len(pts)}/{len(pts)} points OK")


# ================================================================== 3. AT split and F_res
def max_invariant_subspace(K, idx):
    """Maximal K-invariant subspace of span{e_i : i in idx}, exact over Q (sympy)."""
    dim = K.shape[0]
    U = sp.zeros(dim, len(idx))
    for c, i in enumerate(idx):
        U[i, c] = 1
    Km = sp.Matrix(K.tolist())
    while True:
        m = U.shape[1]
        if m == 0:
            return None
        S = U.row_join(-(Km * U))
        ns = S.nullspace()
        C = sp.Matrix([[v[m + j] for j in range(m)] for v in ns]).T if ns else sp.zeros(m, 0)
        if C.shape[1] == m:
            return U
        U = U * C
        for c in range(U.shape[1]):
            if any(U[r, c] != 0 for r in range(dim)):
                den = sp.ilcm(*[sp.fraction(sp.nsimplify(U[r, c]))[1] for r in range(dim)])
                if den != 1:
                    U[:, c] = U[:, c] * den


def at_factor(K, D, r0):
    """AT_{r0}(lam, Q2) = det((lam - r0) I - u K|_V) with u^2 -> -Q2 (even in u asserted)."""
    idx = [i for i in range(len(D)) if D[i] == r0]
    U = max_invariant_subspace(K, idx)
    if U is None or U.shape[1] == 0:
        return sp.Poly(1, lam, Q2), 0
    m = U.shape[1]
    Km = sp.Matrix(K.tolist())
    KV = (U.T * U).inv() * (U.T * (Km * U))
    x = sp.Symbol("x")
    coeffs = sp.Poly(KV.charpoly(x).as_expr(), x).all_coeffs()
    expr = 0
    for j, a in enumerate(coeffs):
        if a == 0:
            continue
        assert j % 2 == 0, f"AT_{r0}: odd-u coefficient a_{j} nonzero"
        expr += a * (lam - r0) ** (m - j) * (-Q2) ** (j // 2)
    P = sp.Poly(sp.expand(expr), lam, Q2)
    assert all(sp.Integer(c) == c for c in P.coeffs()), f"AT_{r0} non-integer"
    return P, m


# ================================================================== 4. S6 factor shape
def strip_sign_symbolic(P):
    """Sign of P on the OPEN strip lam in (-6,-2), Q2 > 0, PROVED symbolically. The linear
    factors by inspection; a strand factor is rewritten as sum c_{a,b} (lam - r0)^a Q2^b for
    r0 in {-2, -6} and every exponent a is checked even with every c >= 0 (the norm-form
    shape: each AT_{r0} irreducible factor is a Galois-orbit product of
    (lam - r0)^2 + Q2*mu_k^2 with mu_k^2 > 0). Such a form is strictly positive on the open
    strip: mu = lam - r0 != 0 there (r0 is an endpoint) and Q2 > 0, so every nonzero monomial
    is > 0. Returns +1 or -1; 0 means not provable this way (asserted against upstream)."""
    e = sp.expand(P.as_expr())
    if e == Q2:
        return +1
    if e == lam + 2:
        return -1
    if e == lam + 6:
        return +1
    mu = sp.Symbol("mu")
    for r0 in (-2, -6):
        Qp = sp.Poly(sp.expand(e.subs(lam, mu + r0)), mu, Q2)
        mc = list(zip(Qp.monoms(), Qp.coeffs()))
        if all(m[0] % 2 == 0 and c >= 0 for m, c in mc):
            return +1
        if all(m[0] % 2 == 0 and c <= 0 for m, c in mc):
            return -1
    return 0


def degree_partition_modp_sympy(coeffs_desc, p):
    """Degree partition of a univariate int poly mod p, or None if p is bad (lc drop or
    not squarefree)."""
    if coeffs_desc[0] % p == 0:
        return None
    x = sp.Symbol("_x")
    fp = sp.Poly(coeffs_desc, x, modulus=p, symmetric=False)
    if sp.gcd(fp, fp.diff()).degree() > 0:
        return None
    part = []
    for f, m in fp.factor_list()[1]:
        part += [f.degree()] * m
    return sorted(part)


def subset_sums(part, n):
    s = {0}
    for d in part:
        s |= {x + d for x in s if x + d <= n}
    return s


def certify_irreducible_univ(coeffs_desc, max_primes=25, p_start=10 ** 6):
    """Subset-sum irreducibility certificate over Q for a univariate integer polynomial:
    the possible Q-factor degrees lie in every good prime's partition subset sums; the
    intersection {0, n} proves irreducibility over Q. Returns (ok, primes_used)."""
    n = len(coeffs_desc) - 1
    possible, used, p = None, 0, p_start
    while used < max_primes:
        p = int(sp.nextprime(p))
        part = degree_partition_modp_sympy(coeffs_desc, p)
        if part is None:
            continue
        used += 1
        possible = subset_sums(part, n) if possible is None else possible & subset_sums(part, n)
        if possible == {0, n}:
            return True, used
    return False, used


def certify_bivariate_irreducible(G, tag):
    """G irreducible over Q: primitive in lam over Z[Q2] AND some full-degree specialization
    Q2 = w0 subset-sum certified irreducible (a nontrivial bivariate factorization would
    specialize to a nontrivial one, both lam-degrees preserved)."""
    n = G.degree(lam)
    coeff_polys = []
    for dl in range(n + 1):
        cexpr = sum(c * Q2 ** m[1] for m, c in zip(G.monoms(), G.coeffs()) if m[0] == dl)
        coeff_polys.append(sp.Poly(cexpr, Q2))
    cont = coeff_polys[-1]
    for cp in coeff_polys:
        cont = sp.Poly(sp.gcd(cont.as_expr(), cp.as_expr()), Q2)
        if cont.total_degree() == 0:
            break
    assert cont.total_degree() == 0, f"{tag}: not primitive in lam over Z[Q2]"
    for w0 in (2, 3, 5, 7, 11, 13):
        spec = sp.Poly(G.as_expr().subs(Q2, w0), lam)
        if spec.degree() != n:
            continue
        cs = [int(c) for c in spec.all_coeffs()]
        g = 0
        for c in cs:
            g = sp.igcd(g, c)
        cs = [c // g for c in cs]
        ok, used = certify_irreducible_univ(cs)
        if ok:
            print(f"    {tag}: IRREDUCIBLE over Q (specialization Q2={w0}, "
                  f"subset-sum certificate, {used} primes)")
            return
    raise AssertionError(f"{tag}: irreducibility certificate did not land")


def dissect_S6(S6, AT, tag):
    """S6 = content * Q2^a * (lam+2)^b * (lam+6)^c * prod(AT factors)^k * G, exact divisions;
    every prefactor sign-fixed on the strip (proved symbolically, strip_sign_symbolic);
    G irreducible over Q (certified)."""
    P = S6
    cs = [int(c) for c in P.coeffs()]
    content = 0
    for c in cs:
        content = sp.igcd(content, c)
    if content > 1:
        P = sp.Poly(P.as_expr() / content, lam, Q2)
    cands = [(sp.Poly(Q2, lam, Q2), "Q2"), (sp.Poly(lam + 2, lam, Q2), "(lam+2)"),
             (sp.Poly(lam + 6, lam, Q2), "(lam+6)")]
    for f, m in sp.factor_list(AT.as_expr(), lam, Q2)[1]:
        fp = sp.Poly(f, lam, Q2)
        cands.append((fp, f"ATfac(deg {fp.degree(lam)},{fp.degree(Q2)})"))
    found = []
    for cand, nm in cands:
        k = 0
        while True:
            q, r = sp.div(P, cand, lam, Q2)
            if not r.is_zero:
                break
            P = sp.Poly(q, lam, Q2)
            k += 1
        if k:
            sgn = strip_sign_symbolic(cand)
            assert sgn != 0, f"{tag}: prefactor {nm} sign on the strip not symbolically provable"
            found.append((nm, k, sgn))
    print(f"  {tag}: content {content}, prefactors {[(nm, k, '+' if s > 0 else '-') for nm, k, s in found]}, "
          f"core G deg ({P.degree(lam)},{P.degree(Q2)})")
    certify_bivariate_irreducible(P, f"{tag} core G")
    return P


# ================================================================== 5. mod-p univariate kit
def polytrim(f):
    i = 0
    while i < len(f) - 1 and f[i] == 0:
        i += 1
    return f[i:]


def is_zero(f):
    return len(f) == 1 and f[0] == 0


def polymod(f, g, p):
    f = f.copy() % p
    ginv = pow(int(g[0]), p - 2, p)
    dg = len(g) - 1
    while len(f) - 1 >= dg and len(f) > 0:
        c = f[0] % p
        if c:
            f[:dg + 1] = (f[:dg + 1] - ((c * ginv) % p) * g) % p
        f = f[1:]
        if len(f) == 0:
            break
    return polytrim(f) if len(f) else np.array([0], dtype=np.int64)


def polydivmod(f, g, p):
    f = f.copy() % p
    ginv = pow(int(g[0]), p - 2, p)
    dg = len(g) - 1
    qout = []
    while len(f) - 1 >= dg:
        q = (f[0] % p) * ginv % p
        qout.append(q)
        if q:
            f[:dg + 1] = (f[:dg + 1] - q * g) % p
        f = f[1:]
    r = polytrim(f) if len(f) else np.array([0], dtype=np.int64)
    return np.array(qout if qout else [0], dtype=np.int64), r


def polygcd(f, g, p):
    f, g = polytrim(f % p), polytrim(g % p)
    while not is_zero(g):
        f, g = g, polymod(f, g, p)
    return (f * pow(int(f[0]), p - 2, p)) % p


def polymul(f, g, p):
    ff = np.convolve(f.astype(object), g.astype(object))
    return np.array([int(x) % p for x in ff], dtype=np.int64)


def polyderiv(f, p):
    n = len(f) - 1
    if n == 0:
        return np.array([0], dtype=np.int64)
    return (f[:-1] * np.arange(n, 0, -1)) % p


def res_modp(f, g, p):
    """Resultant of two descending int64 coefficient arrays mod p (subresultant-free PRS)."""
    f, g = polytrim(f % p), polytrim(g % p)
    if is_zero(f) or is_zero(g):
        return 0
    res = 1
    while True:
        df, dg = len(f) - 1, len(g) - 1
        if dg == 0:
            return (res * pow(int(g[0]), df, p)) % p
        r = polymod(f, g, p)
        if is_zero(r):
            return 0
        res = (res * pow(int(g[0]), df - (len(r) - 1), p)) % p
        if (df * dg) % 2 == 1:
            res = (-res) % p
        f, g = g, r


def yun_sqf(f, p):
    """Yun's squarefree decomposition mod p (p >> deg): [(g_i, i)], f ~ prod g_i^i."""
    f = (f * pow(int(f[0]), p - 2, p)) % p
    df = polyderiv(f, p)
    g = polygcd(f, df, p)
    if len(g) == 1:
        return [(f, 1)]
    out = []
    c, _ = polydivmod(f, g, p)
    w = polydivmod(df, g, p)[0]
    i = 1
    while len(c) > 1:
        cp = polyderiv(c, p)
        m = max(len(w), len(cp))
        y = polytrim((np.pad(w, (m - len(w), 0)) - np.pad(cp, (m - len(cp), 0))) % p)
        if is_zero(y):
            gi = (c * pow(int(c[0]), p - 2, p)) % p
        else:
            gi = polygcd(c, y, p)
        if len(gi) > 1:
            out.append((gi, i))
        c, _ = polydivmod(c, gi, p)
        w = np.array([0], dtype=np.int64) if is_zero(y) else polydivmod(y, gi, p)[0]
        i += 1
    return out


def ddf_partition(f, p):
    """Distinct-degree factorization partition of a squarefree monic f mod p."""
    f = (f * pow(int(f[0]), p - 2, p)) % p
    part = []
    h = np.array([1, 0], dtype=np.int64)
    d = 0
    while len(f) - 1 > 0:
        d += 1
        if 2 * d > len(f) - 1:
            part.append(len(f) - 1)
            break
        hp = np.array([1], dtype=np.int64)
        base = h.copy()
        e = p
        while e:
            if e & 1:
                hp = polymod(polymul(hp, base, p), f, p)
            e >>= 1
            if e:
                base = polymod(polymul(base, base, p), f, p)
        h = hp
        m = max(len(h), 2)
        hpad = np.pad(h, (m - len(h), 0)).copy()
        hpad[-2] = (hpad[-2] - 1) % p
        diff = polytrim(hpad)
        g = f.copy() if is_zero(diff) else polygcd(diff, f, p)
        if len(g) > 1:
            deg = len(g) - 1
            assert deg % d == 0
            part += [d] * (deg // d)
            f, _ = polydivmod(f, g, p)
            h = polymod(h, f, p)
    return sorted(part)


def _modp_kit_selftest(p=33554393):
    """Resultant / Yun / DDF vs sympy on random polynomials (gate-first)."""
    rng = np.random.default_rng(3)
    x = sp.Symbol("x")
    for trial in range(20):
        f = sp.Poly([int(v) for v in rng.integers(-9, 10, 7)] + [1], x)
        g = sp.Poly([int(v) for v in rng.integers(-9, 10, 5)] + [3], x)
        if f.degree() < 1 or g.degree() < 1:
            continue
        want = int(sp.resultant(f, g)) % p
        got = res_modp(np.array([int(c) for c in f.all_coeffs()], dtype=np.int64),
                       np.array([int(c) for c in g.all_coeffs()], dtype=np.int64), p)
        assert got == want, f"res selftest {trial}"
    for trial in range(8):
        f1 = sp.Poly([1] + [int(v) for v in rng.integers(-5, 6, 3)], x)
        f2 = sp.Poly([1] + [int(v) for v in rng.integers(-5, 6, 4)], x)
        f = f1 * f2 ** 2
        arr = np.array([int(c) % p for c in f.all_coeffs()], dtype=np.int64)
        layers = yun_sqf(arr, p)
        degs = sorted((len(g) - 1, m) for g, m in layers)
        fp = sp.Poly(f.as_expr(), x, modulus=p, symmetric=False)
        want = sorted((sp.Poly(gg, x).degree(), mm) for gg, mm in fp.sqf_list()[1])
        assert degs == want, f"sqf selftest {trial}"
        for g, m in layers:
            pat = ddf_partition(g, p)
            fpg = sp.Poly([int(c) for c in g], x, modulus=p, symmetric=False)
            wantp = sorted(sum(([ff.degree()] * mm for ff, mm in fpg.factor_list()[1]), []))
            assert pat == wantp, f"ddf selftest {trial}"


# ================================================================== 6. bivariate mod-p engine
def poly_to_mat(P):
    """sympy Poly in (lam, Q2) -> int coefficient matrix C[dl][dw] (both ascending)."""
    dl_max = P.degree(lam)
    dw_max = max(P.degree(Q2), 0)
    C = [[0] * (dw_max + 1) for _ in range(dl_max + 1)]
    for m, c in zip(P.monoms(), P.coeffs()):
        C[m[0]][m[1]] = int(c)
    return C


def eval_nodes_modp(C, nodes, p):
    """Coefficient matrix -> (len(nodes), deg_lam+1) lam-coefficient arrays at Q2 = node."""
    Cm = np.array([[int(c) % p for c in row] for row in C], dtype=np.int64)
    vals = np.zeros((len(nodes), Cm.shape[0]), dtype=np.int64)
    wpow = np.ones(len(nodes), dtype=np.int64)
    nodes_arr = np.array(nodes, dtype=np.int64) % p
    for dw in range(Cm.shape[1]):
        vals = (vals + np.outer(wpow, Cm[:, dw])) % p
        wpow = (wpow * nodes_arr) % p
    return vals


def interp_modp(nodes, ys, p):
    """Newton interpolation of values ys at integer nodes -> descending trimmed array mod p."""
    npts = len(nodes)
    coef = np.array(ys, dtype=np.int64) % p
    wn = np.array(nodes, dtype=np.int64)
    for j in range(1, npts):
        d = (wn[j:] - wn[:-j]) % p
        lut = {int(x): pow(int(x), p - 2, p) for x in np.unique(d)}
        dinv = np.array([lut[int(x)] for x in d], dtype=np.int64)
        coef[j:] = ((coef[j:] - coef[j - 1:-1]) % p) * dinv % p
    out = np.zeros(npts, dtype=np.int64)
    acc = np.zeros(npts, dtype=np.int64)
    acc[0] = 1
    deg_acc = 0
    for j in range(npts):
        out[:deg_acc + 1] = (out[:deg_acc + 1] + coef[j] * acc[:deg_acc + 1]) % p
        if j == npts - 1:
            break
        new = np.zeros(npts, dtype=np.int64)
        new[1:deg_acc + 2] = acc[:deg_acc + 1]
        new[:deg_acc + 1] = (new[:deg_acc + 1] - wn[j] * acc[:deg_acc + 1]) % p
        acc = new % p
        deg_acc += 1
    return polytrim(out[::-1].copy())


def resultant_modp(CF, CG, nodes, p, derivative=False):
    """Res_lam(F, G)(Q2) mod p by evaluation at Q2-nodes + interpolation (descending array).
    With derivative=True, G is ignored and Res_lam(F, dF/dlam) is computed (the discriminant
    itself: F is monic in lam, so the lc factor is 1)."""
    Fv = eval_nodes_modp(CF, nodes, p)
    Gv = None if derivative else eval_nodes_modp(CG, nodes, p)
    vals = []
    for j in range(len(nodes)):
        fj = polytrim(Fv[j][::-1].copy())
        gj = polyderiv(fj, p) if derivative else polytrim(Gv[j][::-1].copy())
        vals.append(res_modp(fj, gj, p))
    return interp_modp(nodes, vals, p)


def det_modp(M, p):
    """Determinant of an int64 matrix mod p by Gaussian elimination (entries in [0, p))."""
    M = M % p
    n = M.shape[0]
    det = 1
    for k in range(n):
        piv = k + int(np.argmax(M[k:, k] != 0))
        if M[piv, k] == 0:
            return 0
        if piv != k:
            M[[k, piv]] = M[[piv, k]]
            det = (-det) % p
        det = (det * int(M[k, k])) % p
        inv = pow(int(M[k, k]), p - 2, p)
        if k + 1 < n:
            factors = (M[k + 1:, k] * inv) % p
            M[k + 1:, k:] = (M[k + 1:, k:] - factors[:, None] * M[k, k:]) % p
    return det


def psc_matrix(fdesc, gdesc, j):
    """The square matrix whose determinant is the j-th principal subresultant coefficient
    psc_j(f, g): rows x^{n-j-1}f..f and x^{m-j-1}g..g, columns of degrees m+n-j-1 down to j
    (m = deg f > n = deg g). psc_0 = Res(f, g); deg gcd(f, g) >= k <=> psc_0 = .. =
    psc_{k-1} = 0."""
    m, n = len(fdesc) - 1, len(gdesc) - 1
    size = m + n - 2 * j
    T = m + n - j - 1
    M = np.zeros((size, size), dtype=np.int64)
    row = 0
    for s in range(n - j - 1, -1, -1):
        lo = max(s, j)
        M[row, T - (s + m): T - lo + 1] = fdesc[: s + m - lo + 1]
        row += 1
    for s in range(m - j - 1, -1, -1):
        lo = max(s, j)
        M[row, T - (s + n): T - lo + 1] = gdesc[: s + n - lo + 1]
        row += 1
    return M


def psc1_poly_modp(CF, p):
    """psc_1(F_res, dF_res/dlam) as a polynomial in Q2 mod p (descending trimmed array),
    by nodewise determinants + interpolation. psc_1 is a FIXED determinant in F_res's
    integer coefficient polynomials, so its reduction mod p commutes unconditionally."""
    m = len(CF) - 1
    degw = max(len(row) for row in CF) - 1
    size = 2 * m - 3
    nodes = list(range(size * degw + 2))
    Fv = eval_nodes_modp(CF, nodes, p)
    deriv_mult = np.arange(m, 0, -1, dtype=np.int64)
    vals = []
    for jn in range(len(nodes)):
        fdesc = Fv[jn][::-1].copy() % p                # lc = 1 (F monic): degree stable
        gdesc = (fdesc[:-1] * deriv_mult) % p
        vals.append(det_modp(psc_matrix(fdesc, gdesc, 1), p))
    return interp_modp(nodes, vals, p)


def _psc_selftest(p=33554393, p2=33554383):
    """Gate-first: det vs sympy, and the psc vanishing criterion on pairs (f, f') with a
    known gcd degree k (f = f1 * h^(k... built as f1 * h^2 with h squarefree deg k)."""
    rng = np.random.default_rng(11)
    x = sp.Symbol("x")
    for trial in range(5):
        A = rng.integers(-9, 10, size=(6, 6))
        want = int(sp.Matrix(A.tolist()).det()) % p
        assert det_modp(A.astype(np.int64) % p, p) == want, "det selftest"
    for k in (1, 2, 3):
        while True:
            h = sp.Poly([1] + [int(v) for v in rng.integers(-4, 5, k)], x)
            f1 = sp.Poly([1] + [int(v) for v in rng.integers(-4, 5, 4)], x)
            f = f1 * h ** 2
            g = f.diff()
            if sp.gcd(f, g).degree() == k:              # generic: gcd = h exactly
                break
        fdesc = np.array([int(c) for c in f.all_coeffs()], dtype=np.int64)
        gdesc = np.array([int(c) for c in g.all_coeffs()], dtype=np.int64)
        for j in range(k + 1):
            v1 = det_modp(psc_matrix(fdesc % p, gdesc % p, j), p)
            v2 = det_modp(psc_matrix(fdesc % p2, gdesc % p2, j), p2)
            if j < k:
                assert v1 == 0 and v2 == 0, f"psc_{j} != 0 with gcd degree {k}"
            else:
                assert v1 != 0 or v2 != 0, f"psc_{k} = 0 at two primes with gcd degree {k}"


def valuation_at_zero(f):
    v = 0
    while v < len(f) and f[len(f) - 1 - v] == 0:
        v += 1
    return v


def disc_layers_modp(CF, nodes, p, table):
    """Disc mod p: assert deg, v, layer degrees vs the certified table, max multiplicity 2;
    return (Dp, A1p, A2p, vD)."""
    Dp = resultant_modp(CF, None, nodes, p, derivative=True)
    vD = valuation_at_zero(Dp)
    layers = yun_sqf(polytrim(Dp[:len(Dp) - vD].copy()), p)
    laydeg = {m: len(g) - 1 for g, m in layers}
    assert max(laydeg) <= 2, f"multiplicity >= 3 layer at p={p}: beta-exotic clash!"
    assert vD == table["v"], f"v_w(D) = {vD} != {table['v']} at p={p}"
    assert laydeg == {1: table["A1"], 2: table["A2"]}, \
        f"layer degrees {laydeg} != table {table} at p={p}"
    A1p = next(g for g, m in layers if m == 1)
    A2p = next(g for g, m in layers if m == 2)
    return Dp, A1p, A2p, vD


# ================================================================== 6b. exact N = 5 inventory
def exact_disc_inventory_n5(s, F):
    """N = 5 only: reconstruct disc_Lam(F_res) EXACTLY over Z (integer-node interpolation
    with exact sympy resultants), split it exactly, factor A1 and A2 (both irreducible over
    Q), and assert the real-positive root inventory: A1's real positive roots are exactly the
    sector's two seed loci w = q*^2; A2 has no real positive root (E) resp. exactly one (O,
    w ~ 5.100831). Returns dict(D, A1, A2) as exact sympy Polys in Q2."""
    t0 = time.time()
    table = DISC_TABLE[5][s]
    degD = table["v"] + table["A1"] + 2 * table["A2"]
    CF = poly_to_mat(F)
    nodes = list(range(degD + 3))
    vals = []
    for w0 in nodes:
        coeffs_desc = []
        for dl in range(len(CF) - 1, -1, -1):
            c = 0
            for dw in range(len(CF[dl]) - 1, -1, -1):        # Horner in w0, exact ints
                c = c * w0 + CF[dl][dw]
            coeffs_desc.append(c)
        fw = sp.Poly(coeffs_desc, lam)
        vals.append(int(sp.resultant(fw, fw.diff(lam))))
    # exact Newton interpolation at integer nodes 0..n-1
    n = len(nodes)
    coef = [Fraction(v) for v in vals]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / j
    poly = [Fraction(0)] * n                                 # ascending in w
    acc = [Fraction(0)] * n
    acc[0] = Fraction(1)
    dega = 0
    for j in range(n):
        for k in range(dega + 1):
            poly[k] += coef[j] * acc[k]
        if j == n - 1:
            break
        new = [Fraction(0)] * n
        for k in range(dega + 1):
            new[k + 1] += acc[k]
            new[k] -= j * acc[k]
        acc = new
        dega += 1
    assert all(c.denominator == 1 for c in poly), "exact disc interpolation not integral"
    ints = [int(c) for c in poly]
    while ints and ints[-1] == 0:
        ints.pop()
    assert len(ints) - 1 == degD, f"exact disc degree {len(ints) - 1} != {degD}"
    Dx = sp.Poly(ints[::-1], Q2)

    _, sq = Dx.sqf_list()
    v = 0
    A1x = A2x = None
    for f, mult in sq:
        fp = sp.Poly(f, Q2)
        if fp.as_expr() == Q2:
            v = mult
        elif mult == 1:
            A1x = fp
        elif mult == 2:
            A2x = fp
        else:
            raise AssertionError(f"unexpected exact layer multiplicity {mult}")
    assert v == table["v"] and A1x.degree() == table["A1"] and A2x.degree() == table["A2"], \
        f"exact layers (v={v}, A1={A1x.degree()}, A2={A2x.degree()}) != table {table}"
    assert len(sp.factor_list(A1x.as_expr())[1]) == 1, "exact A1 not irreducible over Q"
    assert len(sp.factor_list(A2x.as_expr())[1]) == 1, "exact A2 not irreducible over Q"

    # real-positive root inventory (exact isolation, then float comparison)
    r1 = sorted(float(r) for r in A1x.real_roots() if r > 0)
    expected = sorted(q0 * q0 for (sec, q0, _, _) in SEEDS[5] if sec == s)
    assert len(r1) == 2, f"A1 real positive roots: {r1} (expected 2 seed loci)"
    assert all(abs(a - b) < 1e-4 for a, b in zip(r1, expected)), \
        f"A1 real positive roots {r1} != seed loci {expected}"
    r2 = sorted(float(r) for r in A2x.real_roots() if r > 0)
    if s == "E":
        assert r2 == [], f"A2 (E) has unexpected real positive roots {r2}"
    else:
        assert len(r2) == 1 and abs(r2[0] - 5.100831) < 1e-4, \
            f"A2 (O) real positive roots {r2} != [5.100831]"
    print(f"  [{s}] EXACT disc over Z: deg {degD}, layers (v={v}, A1={A1x.degree()}, "
          f"A2={A2x.degree()}), A1 and A2 irreducible; A1 real>0 roots {[f'{x:.6f}' for x in r1]}"
          f" = seed loci; A2 real>0 roots {[f'{x:.6f}' for x in r2]}  [{time.time() - t0:.1f}s]")

    if s == "O":                                             # the single A2 point: diabolic
        mp.dps = 60
        w0 = mpf(str(sp.N([r for r in A2x.real_roots() if r > 0][0], 50)))
        coeffs = []
        for dl in range(len(CF) - 1, -1, -1):
            c = mpf(0)
            for dw in range(len(CF[dl]) - 1, -1, -1):
                c = c * w0 + CF[dl][dw]
            coeffs.append(c)
        rts = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
        close = [(i, j) for i in range(len(rts)) for j in range(i + 1, len(rts))
                 if abs(rts[i] - rts[j]) < mpf("1e-10")]
        assert len(close) == 1, f"A2 point carries {len(close)} double roots (expected 1)"
        i0, j0 = close[0]
        assert abs(rts[i0].imag) < mpf("1e-8") and abs(rts[j0].imag) < mpf("1e-8"), \
            "the A2-point double root is not real"
        print(f"  [O] A2 root w={mp.nstr(w0, 9)}: exactly ONE real double lambda-root of "
              f"F_res at lam={mp.nstr(rts[i0].real, 9)} (diabolic class, no coincident pair)")
    return dict(D=Dx, A1=A1x, A2=A2x)


# ================================================================== 7. seed spot-checks
def refine_and_check_seeds(N, polys, F_res):
    """Newton-refine each forced seed as a double root of its sector F_res (exact Jacobian,
    dps 80; mpmath findroot's absolute tolerance is useless here, the coefficients carry ~8
    digits of cancellation and the locus must be pinned far below that), then assert
    S6 != 0, S6*ST < 0, kappa_-2 = -S6/ST > 0 matching the eigenvector-route reference."""
    mp.dps = 80

    def mp_eval(P, l0, w0):
        tot = mpf(0)
        for (dl, dw), c in zip(P.monoms(), P.coeffs()):
            tot += int(c) * (l0 ** dl) * (w0 ** dw)
        return tot

    lamb = {}
    for s in ("E", "O"):
        Fe = F_res[s].as_expr()
        lamb[s] = dict(
            F=sp.lambdify((lam, Q2), Fe, "mpmath"),
            dF=sp.lambdify((lam, Q2), sp.diff(Fe, lam), "mpmath"),
            J=[sp.lambdify((lam, Q2), sp.diff(Fe, lam), "mpmath"),
               sp.lambdify((lam, Q2), sp.diff(Fe, Q2), "mpmath"),
               sp.lambdify((lam, Q2), sp.diff(Fe, lam, 2), "mpmath"),
               sp.lambdify((lam, Q2), sp.diff(Fe, lam, Q2), "mpmath")])

    print(f"  {'q*':>11} {'lam*':>11} sec {'S6_s':>12} {'ST_s':>12} "
          f"{'kappa_-2':>10} {'kappa_ref':>10}")
    for s, q0, l0, kref in SEEDS[N]:
        so = "O" if s == "E" else "E"
        l_s, w_s = mpf(l0), mpf(q0) ** 2
        dl_step = dw_step = mpf(1)
        for _ in range(80):
            f1 = lamb[s]["F"](l_s, w_s)
            f2 = lamb[s]["dF"](l_s, w_s)
            j11, j12, j21, j22 = (fn(l_s, w_s) for fn in lamb[s]["J"])
            det = j11 * j22 - j12 * j21
            dl_step = (f1 * j22 - f2 * j12) / det
            dw_step = (j11 * f2 - j21 * f1) / det
            l_s, w_s = l_s - dl_step, w_s - dw_step
            if abs(dl_step) < mpf("1e-45") and abs(dw_step) < mpf("1e-45"):
                break
        assert abs(dl_step) < mpf("1e-40") and abs(dw_step) < mpf("1e-40"), \
            f"Newton did not converge at q*={q0}"
        q_s = mp.sqrt(w_s)
        assert abs(l_s - l0) < 0.05 and abs(q_s - q0) < 0.01, f"Newton drifted at q*={q0}"
        S6_s = mp_eval(polys[s]["S6"], l_s, w_s)
        ST_s = mp_eval(polys[s]["ST"], l_s, w_s)
        chi_s = mp_eval(polys[s]["chi"], l_s, w_s)
        chi_o = mp_eval(polys[so]["chi"], l_s, w_s)
        S6_full = S6_s * chi_o + mp_eval(polys[so]["S6"], l_s, w_s) * chi_s
        ST_full = ST_s * chi_o + mp_eval(polys[so]["ST"], l_s, w_s) * chi_s
        kappa = float(-S6_full / ST_full)
        assert abs(S6_s) > mpf("1e-20"), f"S6 vanishes at seed q*={q0}!"
        assert S6_s * ST_s < 0 and S6_full * ST_full < 0, f"S6*ST not < 0 at q*={q0}"
        assert kappa > 0, f"kappa_-2 not positive at q*={q0}"
        assert abs(kappa - kref) < 1e-5, \
            f"kappa {kappa:.7f} != reference {kref} at q*={q0} (handshake break)"
        print(f"  {mp.nstr(q_s, 10):>11} {mp.nstr(l_s, 10):>11} R-{s} "
              f"{mp.nstr(S6_s, 4):>12} {mp.nstr(ST_s, 4):>12} {kappa:>+10.6f} {kref:>+10.6f}")


# ================================================================== driver
def run(N):
    t_start = time.time()
    print(f"\n================================ N = {N} ================================")
    all_primes = make_primes(60)
    check_primes, work_primes = all_primes[:3], all_primes[3:]

    # ---- 1+2: sectors and exact polynomials
    secs = integer_sectors(N)
    for s in ("E", "O"):
        assert secs[s]["dim"] == SECTOR_DIM[N][s], f"sector dim mismatch at {s}"
    exact = {s: compute_exact_polys(s, secs[s], work_primes, check_primes,
                                    exact_over_Z=(N == 5)) for s in ("E", "O")}
    polys = {}
    for s in ("E", "O"):
        dim = exact[s]["dim"]
        polys[s] = dict(chi=to_poly(exact[s]["chi"], dim),
                        S6=to_poly(exact[s]["S6"], dim - 1),
                        ST=to_poly(exact[s]["ST"], dim - 1))
    full_block_50digit_check(N, polys, npts=3 if N == 5 else 2)

    # ---- 2: AT split, F_res, S6 factor shape
    F_res, AT = {}, {}
    for s in ("E", "O"):
        t0 = time.time()
        A = sp.Poly(1, lam, Q2)
        for r0 in (-2, -6):
            P_r0, m = at_factor(secs[s]["K"], [int(x) for x in secs[s]["D"]], r0)
            A = sp.Poly(A.as_expr() * P_r0.as_expr(), lam, Q2)
        q_, r_ = sp.div(polys[s]["chi"], A, lam, Q2)
        assert r_.is_zero, f"chi not divisible by AT in sector {s}"
        F = sp.Poly(q_, lam, Q2)
        assert F.degree(lam) == FRES_DEG[N][s], \
            f"F_res degree {F.degree(lam)} != {FRES_DEG[N][s]} in sector {s}"
        assert sp.LC(F.as_expr(), lam) == 1, f"F_res not monic in lam in sector {s}"
        F_res[s], AT[s] = F, A
        print(f"  [{s}] chi = AT(deg {A.degree(lam)}) * F_res(deg {F.degree(lam)}, "
              f"w-deg {F.degree(Q2)}), F_res monic  [{time.time() - t0:.1f}s]")
        dissect_S6(polys[s]["S6"], A, f"[{s}] S6")

    # ---- 3+3b+4: disc layers (+ exact N=5 inventory), coincident-EP2 exclusion,
    #              A1 irreducibility, gcd certificates
    _modp_kit_selftest()
    _psc_selftest()
    print("  mod-p kit + psc self-tests OK")
    cert_primes = make_primes(8, below=2 ** 25 - 60)     # any good prime proves; 3 asserted
    for s in ("E", "O"):
        so = "O" if s == "E" else "E"
        table = DISC_TABLE[N][s]
        exact_inv = exact_disc_inventory_n5(s, F_res[s]) if N == 5 else None
        CF = poly_to_mat(F_res[s])
        CS = poly_to_mat(polys[s]["S6"])
        CX = poly_to_mat(polys[so]["chi"])
        degF_l, degF_w = F_res[s].degree(lam), F_res[s].degree(Q2)
        degS_l, degS_w = polys[s]["S6"].degree(lam), polys[s]["S6"].degree(Q2)
        degX_l, degX_w = polys[so]["chi"].degree(lam), polys[so]["chi"].degree(Q2)
        nodesD = list(range(degF_l * degF_w + (degF_l - 1) * degF_w + 2))
        nodesR = list(range(degF_l * degS_w + degS_l * degF_w + 2))
        nodesX = list(range(degF_l * degX_w + degX_l * degF_w + 2))
        # good-prime conditions, stated and checked:
        #   (i) F_res and chi_other are MONIC in lam (asserted above / by construction), so
        #       Res_lam commutes with reduction mod p once deg_lam of the second argument
        #       is preserved;
        #  (ii) p does not divide lc_lam(S6) (= the rung-6 state count, a tiny integer);
        # (iii) deg/valuation/layer degrees of D mod p match the certified table. At N = 5
        #       the identification A1p = A1 mod p is then ASSERTED against the exact A1; at
        #       N = 7 it is the conditional premise the module docstring spells out (a mod-p
        #       merge of two A1 roots is invisible here; finitely many primes can be bad).
        lcS6 = int(sp.LC(polys[s]["S6"].as_expr(), lam))
        r1_degs, gcd_r1_disc, psc_shapes = set(), set(), set()
        a1_images = []                               # (p, A1p) for the irreducibility cert
        n_ok = 0
        for p in cert_primes[:MIN_CERT_PRIMES]:
            t0 = time.time()
            assert lcS6 % p != 0, f"bad prime {p} for lc(S6)"
            Dp, A1p, A2p, vD = disc_layers_modp(CF, nodesD, p, table)
            a1_images.append((p, A1p))
            if exact_inv is not None:
                # N = 5: the mod-p objects are asserted to BE the exact ones' reductions
                Dx_p = polytrim(np.array([int(c) % p for c in exact_inv["D"].all_coeffs()],
                                         dtype=np.int64))
                assert np.array_equal(Dx_p, Dp), f"exact D mod {p} != interpolated Dp"
                for exact_poly, modp_layer, nm in ((exact_inv["A1"], A1p, "A1"),
                                                   (exact_inv["A2"], A2p, "A2")):
                    ax = polytrim(np.array([int(c) % p for c in exact_poly.all_coeffs()],
                                           dtype=np.int64))
                    ax = (ax * pow(int(ax[0]), p - 2, p)) % p        # monic-normalize
                    assert np.array_equal(ax, modp_layer), \
                        f"exact {nm} mod {p} != Yun mult layer (identification broken)"
            # coincident-EP2 exclusion: gcd(D/w^v, psc_1/w^v') = 1
            P1p = psc1_poly_modp(CF, p)
            vP = valuation_at_zero(P1p)
            Dt = polytrim(Dp[:len(Dp) - vD].copy())
            Pt = polytrim(P1p[:len(P1p) - vP].copy())
            g3 = polygcd(Dt, Pt, p)
            assert len(g3) - 1 == 0, \
                f"coincident-EP2 certificate failed: gcd(disc~, psc_1~) != 1 at p={p} sector {s}"
            psc_shapes.add((len(P1p) - 1, vP))
            R1p = resultant_modp(CF, CS, nodesR, p)
            Xp = resultant_modp(CF, CX, nodesX, p)
            g1 = polygcd(R1p, A1p, p)
            g2 = polygcd(Xp, A1p, p)
            assert len(g1) - 1 == 0, f"gcd(Res(F,S6), A1) != 1 at p={p} sector {s}"
            assert len(g2) - 1 == 0, f"gcd(Res(F,chi_{so}), A1) != 1 at p={p} sector {s}"
            r1_degs.add(len(R1p) - 1)
            gd = polygcd(R1p, Dp, p)
            gcd_r1_disc.add(len(gd) - 1)
            n_ok += 1
            print(f"  [{s}] p={p}: disc layers OK (max mult 2), gcd(Res(F,S6), A1) = 1, "
                  f"gcd(Res(F,chi_{so}), A1) = 1, gcd(disc~, psc_1~) = 1 (no coincident "
                  f"EP2 pair)  [{time.time() - t0:.1f}s]")
        assert n_ok >= MIN_CERT_PRIMES
        assert len(r1_degs) == 1, f"deg Res(F,S6) not prime-stable: {r1_degs}"
        assert len(psc_shapes) == 1, f"(deg, v) of psc_1 not prime-stable: {psc_shapes}"
        assert len(gcd_r1_disc) == 1, f"deg gcd(R1, disc) not prime-stable: {gcd_r1_disc}"
        pred = table["v"] + 2 * table["A2"]
        if N == 5:
            # the A2^2-split of the resultant is EXACT at N = 5 (sympy scouts): assert it
            assert gcd_r1_disc == {GCD_R1_DISC_N5[s]} and GCD_R1_DISC_N5[s] == pred, \
                f"gcd(R1, disc) degree {gcd_r1_disc} != exact split {GCD_R1_DISC_N5[s]}"
            print(f"  [{s}] deg gcd(Res(F,S6), disc) = {pred} = v + 2*deg A2 "
                  f"(the EXACT N=5 split, asserted)")
        else:
            print(f"  [{s}] deg gcd(Res(F,S6), disc) = {sorted(gcd_r1_disc)[0]} "
                  f"(prediction v + 2*deg A2 = {pred}; multi-prime EVIDENCE, not asserted)")
        # A1 irreducibility over Q: subset-sum across primes; start from the certificate
        # primes' images, extend with fresh good primes until {0, n} lands (budget 20 extra)
        possible, used, nA1 = None, 0, table["A1"]
        p_next = min(pp for pp, _ in a1_images)
        extra_budget = 20
        while True:
            if used < len(a1_images):
                p, A1p = a1_images[used]
            else:
                if extra_budget == 0:
                    break
                extra_budget -= 1
                p_next = int(sp.prevprime(p_next))
                p = p_next
                try:
                    _, A1p, _, _ = disc_layers_modp(CF, nodesD, p, table)
                except AssertionError:
                    continue                          # bad prime for the layer shape: skip
            part = ddf_partition(A1p, p)
            possible = subset_sums(part, nA1) if possible is None \
                else possible & subset_sums(part, nA1)
            used += 1
            if possible == {0, nA1}:
                break
        assert possible == {0, nA1}, \
            f"A1 irreducibility not certified in sector {s}: possible {sorted(possible)}"
        print(f"  [{s}] A1 (deg {nA1}): IRREDUCIBLE over Q (subset-sum certificate, "
              f"{used} primes)")

    # ---- 5: seed spot-checks
    refine_and_check_seeds(N, polys, F_res)

    print(f"  Coincident-EP2 pairs excluded at every w != 0 (both sectors, psc_1 certificate)"
          + (", N=5 inventory exact over Z" if N == 5 else "") + ".")
    qual = "" if N == 5 else (" (modulo the mod-p layer-identification premise; the base "
                              "polynomials are PROOF grade since 2026-07-16, see docstring)")
    print(f"\nO2B NONVANISHING CERTIFIED at N={N}: S6 != 0 on the entire simple disc layer "
          f"(all seeds), geometric multiplicity 1 full-block, s6 != 0 at every forced seed. "
          f"PASS{qual}  [{time.time() - t_start:.0f}s]")


if __name__ == "__main__":
    _fl_selftest()
    print("FL self-test OK")
    todo = [5] + ([7] if "7" in sys.argv[1:] else [])
    for n in todo:
        run(n)
    print(f"\nAll assertions passed at N in {todo}.")
