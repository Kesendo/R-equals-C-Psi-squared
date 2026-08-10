"""The (1,3) block at N=6: sector disc-reality and the real-q defective-locus census.

The sideways_spin_ladder arc's offered follow-up (experiments/F89_PATH_K_DIABOLIC.md, the
mechanism paragraph; OpenArcsRegistry entry sideways_spin_ladder): the self-fold is
block-scoped, fold legs fix (p, N/2) blocks at every even N, so the mechanism that gave
N=4 its four real defective (1,2) loci moves to (1,3) at N=6. This script runs that block
and gates what it finds. Everything is measured from below on the block builder of
eta_ladder_blocks.py (build_L, weight-general); gamma = 1 per site, octic q, J = 2q.

What is gated, in order:

  1. The antiunitary self-fold ON the block, entry-wise and exactly (no eigensolver):
     L + Pf conj(L) Pf^T + 2N I == 0.0 at real coupling, Pf the bra-complement permutation
     b -> ~b (f_P of the fold lattice, PROOF_CODIM1_BY_ADDITIVITY section 7; the antiunitary
     leg is F89CrossFoldSimilarityClaim, general in both weights). Exact 0.0: the relation
     is a permutation-signed rearrangement of identical floats.
  2. [Pf, Pr] == 0.0 and [L, Pr] == 0.0 (Pr the site reflection): the fold acts WITHIN each
     R-sector, so the sector discs inherit its reality.
  3. The bipartite gauge and the closure criterion, entry-wise: T = diag((-1)^(sigma(a)+
     sigma(b))), sigma the sum of occupied site indices, satisfies T L T == conj(L) exactly
     on EVERY block (at real J the composite of PROOF_CODIM1_BY_ADDITIVITY section-7
     ingredients (iii) conj L(q) = L(-conj q) and (iv) D L(q) D = L(-q); the block form is
     BetaExoticPerNExclusionClaim's metric, there at p+q = 3), and T commutes with the site reflection iff (p+q)(N-1) is
     even, anticommutes otherwise (both directions exact). So each R-sector of a block is
     conjugation-closed (a REAL matrix family) whenever p+q is even or N is odd (derived);
     on the other parity the gauge swaps the sectors, and non-closure is MEASURED, not
     derived (a second closing antiunitary is excluded by nothing on the page). Gated: the
     identity and the character at (1,3)/(1,2)/(2,3)/(0,3)@N=6, (1,2)@N=4, (1,3)@N=4, then
     the closure verdicts by nearest-neighbour residual: (1,3)@N=6 and (1,3)@N=4 closed at
     machine floor (the latter without being self-folded), (1,2)@N=4, (1,2)@N=6, (2,3)@N=6,
     (0,3)@N=6 NOT closed (residual O(1)) -- the last two being fold-fixed like (1,3), so
     the closure is the GAUGE's parity, not the self-fold's doing. Also gated here: the
     R-even kappa=0 AT triples at lambda = -4 and -8 are SEMISIMPLE (geometric = algebraic
     = 3, third singular value at floor, fourth O(1)), so the |delta| filter below hides no
     Jordan defect, and R-odd carries no exact degeneracy at all (60 distinct values). Consequence for (1,3)@N=6: with the holomorphic s = fold o conj
     (lambda -> -lambda - 2N) also acting within the sector, every defective locus is
     DOUBLED (conjugate-doubled on the fold-fixed line Re = -N, s-doubled off it), every
     observed disc zero has even order, and a sign-change (Sturm-style odd-order) detector
     is blind here, unlike at N=4 (gated in 7).
  4. Sector disc-reality, the arg-from-real device of even_n_literal_real_count.py with a
     degenerate-pair filter (|delta| > 1e-8; the (1,3) R-even sector carries two persistent
     kappa=0 triples at lambda = -4 and -8, whose intra-cluster noise angles would otherwise
     drown the read): (1,3)@N=6 both sectors at machine floor, (1,2)@N=6 both sectors O(1)
     (the failing control), (1,2)@N=4 at machine floor. The floor/fail split is 11 orders
     here; the gates sit at 1e-9 and 0.02 like the committed script's. NOT the committed
     gate's numbers re-read: even_n_literal_real_count.py probes the (2,1) block at the
     CARRIER knob (its q = J) with no pair filter, this script probes at octic q (J = 2q)
     filtered, so the two O(1) reads agree in kind, not in value.
  5. Detector validation on ground truth: the min-gap UNION log|disc|-dip channel on
     (1,2)@N=4 R-even must find EXACTLY the four known loci 0.460212 / 0.854438 / 0.857458
     / 1.738181 (ReferenceDefectiveLoci, the three-signal float scan; the exactness anchor
     is the committed factorization disc(F8) = c q^24 (3q^4+q^2-1)^2 P10(q^2), whose simple
     layer carries the four; no exact Sturm run on P10 exists in the repo), all with
     gap-scaling exponent 0.5.
     VALIDATION TRANSFERS ONLY TO THE SPECIES IT CONTAINS: the N=4 ground truth holds no
     real-lambda locus, and on (1,3)@N=6 this channel misses most of that species (a review
     round measured the miss); the census therefore runs TWO channels, and the second is
     the repo's own instrument for exactly this, the PT-break real-count scan of the real
     seed census ("immune to the closest-pair masking").
  6. THE CENSUS, (1,3)@N=6, octic q in [0.02, 8.0), two channels, AT LEAST THIRTEEN
     defective loci (a FLOOR, not a count: channel A rejects candidates whose refined gap
     stays above 1e-6, every rejection now printed, and a review round found a genuine
     channel-B locus in that bucket, so channel A alone under-counts; nothing bounds the
     census from above short of the exact layer read).
     Channel A (min-gap UNION dips, step 0.001) carries the fold-line species and one
     real-lambda locus:
       R-even (7): q* = 0.4630013 (lambda -6 +/- 0.3864821i), 0.6400121 (+/-1.2687899i),
                   0.7450837 (+/-1.4621529i), 0.9184132 (+/-7.7721940i),
                   1.5109958 (+/-0.9288088i), 3.0922209 (+/-1.5169893i),
                   3.3536850 (real s-pair -6.0787877 / -5.9212123)
       R-odd  (2): q* = 1.9567664 (+/-7.6439946i), 3.2608458 (+/-25.1047962i)
     Channel B (real-count transitions, step 0.001, bisected) carries the real-lambda
     species completely, four loci channel A misses plus 3.3536850 again:
       R-even: q* = 0.7026159 (12->8, s-pair -6.5624 / -5.4376),
               3.3536850 (8->12, s-pair -6.0788 / -5.9212)
       R-odd:  q* = 0.4635508 (12->8, -6.2681 / -5.7319), 0.8029541 (8->4, -6.6052 /
               -5.3948), 1.2273677 (4->0, -6.3761 / -5.6239)
     Union: 8 R-even + 5 R-odd = 13 found. Every locus gap-scaling exponent 0.5 (sqrt
     branch, defective); every fold-line locus conjugate-doubled (the partner pair
     coalesces at the SAME q*; given the exactly-gated closure this is derived, and the
     float check confirms it at the refinement floor, gated < 1e-6); every real-lambda
     locus s-paired, the two MEASURED merging pair means summing to -2N = -12 (gated
     1e-4). Diabolic crossings (exponent 1) are inventoried, printed, AND count-gated
     (10 R-even, 12 R-odd on channel A's grid), not gated as loci; a per-sector gate pins
     that no kept coalescence is unclassified. Edges [0.002, 0.02) and [8.0, 25.0)
     (section 6d): channel A empty at steps 0.0002/0.002, channel B empty at the same
     steps (its species cannot hide from a count except by exact cancellation inside one
     cell, the count-blind tangency the committed real-count gate already names).
  7. The sign-blindness, gated as the contrast it is: disc sign flips (parity of the
     filtered arg) -- ZERO for both (1,3)@N=6 sectors over [0.02, 8.0) at step 0.001
     despite the thirteen loci (every observed zero even-order), exactly FOUR for
     (1,2)@N=4 R-even at step 0.001, one at each known locus (single EPs there, odd-order
     zeros; the conjugate partner lives in the other sector). Both reads need the fine
     step: two odd-order zeros in one grid cell cancel their sign changes (the N=4 twins
     are 0.003 apart, the closest (1,3) loci 0.0028).

Exponent tolerances are a classification, not a fit: this verifier's defective exponents
read 0.500 and its diabolic band 1.0000..1.0006 (the scouting rounds' coarser probes read
0.492..0.508 vs 0.996..1.005), so any threshold in (0.51, 0.99) decides identically; the
gates use |exp - 0.5| < 0.05 and |exp - 1.0| < 0.05.

A linear gated script, not a module: importing it executes the whole suite; copy functions
rather than import. Known instrument seams, why channel B exists and why every count is a
floor: collapse(tol=1e-10) makes the collapsed min-gap discontinuous at a coalescence, the
trisection assumes local unimodality of exactly that function, and collapse also merges an
EXACTLY-converged coalescing pair, so features() then reads the next gap and the refiner
rejects the point (the known N=4 diabolic 0.6590 and the lambda=-6 real-real crossing
0.5439 land in the printed rejected bucket: the diabolic inventory is an inventory, not a
census). The arg-from-real device is one-sided (a floor proves reality via the exact
fold/gauge argument; an O(1) value is a mod-pi residue, only suggestive -- and on a sector
with exact degeneracies the UNFILTERED value is the angle of rounding noise, similarity-
dependent, the committed even_n_literal_real_count.py N=6 read included; its floor side
stands, its O(1) side needs this filter, a recorded follow-up).

Run:  python simulations/disc13_n6_census.py     (asserts everything; ~9 min quiet)
"""
import numpy as np
import sys, os

sys.path.insert(0, os.path.dirname(__file__))
from eta_ladder_blocks import build_L, block_basis, bit

FAIL = []


def gate(cond, msg):
    tag = "PASS" if cond else "FAIL"
    print(f"  [{tag}] {msg}")
    if not cond:
        FAIL.append(msg)


def refl_site(x, Nn):
    y = 0
    for l in range(Nn):
        if x & bit(l, Nn):
            y |= bit(Nn - 1 - l, Nn)
    return y


def reflection_perm(Nn, p, q):
    pairs, index = block_basis(Nn, p, q)
    P = np.zeros((len(pairs), len(pairs)))
    for col, (a, b) in enumerate(pairs):
        P[index[(refl_site(a, Nn), refl_site(b, Nn))], col] = 1.0
    return P


def bra_complement_perm(Nn, p, q):
    pairs, index = block_basis(Nn, p, q)
    mask = (1 << Nn) - 1
    P = np.zeros((len(pairs), len(pairs)))
    for col, (a, b) in enumerate(pairs):
        P[index[(a, (~b) & mask)], col] = 1.0
    return P


def sector_basis(P, even=True):
    wP, VP = np.linalg.eigh((P + P.T) / 2)
    return VP[:, wP > 0.5] if even else VP[:, wP < -0.5]


PCACHE = {}


def eigs(Nn, p, q, Jv, sector):
    Lb, _, _ = build_L(Nn, p, q, [Jv] * (Nn - 1), [1.0] * Nn)
    key = (Nn, p, q)
    if key not in PCACHE:
        PCACHE[key] = reflection_perm(Nn, p, q)
    B = sector_basis(PCACHE[key], even=(sector == "even"))
    return np.linalg.eigvals(B.T @ Lb @ B)


def collapse(w, tol=1e-10):
    w = list(w)
    out = []
    used = [False] * len(w)
    for i in range(len(w)):
        if used[i]:
            continue
        c = [w[i]]
        for j in range(i + 1, len(w)):
            if not used[j] and abs(w[j] - w[i]) < tol:
                c.append(w[j])
                used[j] = True
        out.append(np.mean(c))
    return np.array(out)


def features(w):
    w = collapse(w)
    g, pair, s = np.inf, None, 0.0
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            d = abs(w[i] - w[j])
            s += np.log(d)
            if d < g:
                g, pair = d, (w[i], w[j])
    return g, pair, s


def disc_arg_filtered(w, floor=1e-8):
    s = 0.0
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            d = w[i] - w[j]
            if abs(d) > floor:
                s += 2.0 * np.angle(d)
    return abs((s + np.pi / 2) % np.pi - np.pi / 2)


def conj_closure_residual(w):
    wc = np.conj(w)
    return max(min(abs(x - y) for y in wc) for x in w)


def census(Nn, p, q, sector, q_lo, q_hi, step):
    """Combined detector: min-gap local minima UNION log|disc| dips; refined; classified."""
    ef = lambda qq: eigs(Nn, p, q, 2.0 * qq, sector)
    qs = np.arange(q_lo, q_hi, step)
    feats = [features(ef(qq)) for qq in qs]
    gaps = np.array([f[0] for f in feats])
    F = np.array([f[2] for f in feats])
    cand = set()
    for k in range(1, len(qs) - 1):
        if gaps[k] < gaps[k - 1] and gaps[k] < gaps[k + 1] and gaps[k] < 0.2:
            cand.add(k)
        if F[k] < F[k - 1] and F[k] < F[k + 1] and min(F[k - 1], F[k + 1]) - F[k] > 0.4:
            cand.add(k)
    loci = []
    for k in sorted(cand):
        lo, hi = qs[k - 1], qs[k + 1]
        for _ in range(80):
            m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
            if features(ef(m1))[0] < features(ef(m2))[0]:
                hi = m2
            else:
                lo = m1
        qstar = (lo + hi) / 2
        g, pair, _ = features(ef(qstar))
        if g > 1e-6:
            # no silent bucket: every rejected candidate is printed (a review round found a
            # genuine channel-B locus in exactly this bucket, so channel A's count is a FLOOR)
            print(f"    q~{qstar:.7f}: refined gap {g:.3e} (avoided; not counted)")
            continue
        if any(abs(qstar - q0) < 1e-5 for q0, *_ in loci):
            continue
        lam = (pair[0] + pair[1]) / 2
        ds = (1e-3, 1e-4, 1e-5)
        slopes = []
        for sgn in (+1, -1):
            gs = []
            for d in ds:
                w = ef(qstar + sgn * d)
                near = w[np.argsort(np.abs(w - lam))[:2]]
                gs.append(abs(near[0] - near[1]))
            slopes.append(np.polyfit(np.log(ds), np.log(gs), 1)[0])
        exp = float(np.mean(slopes))
        w = ef(qstar)
        if abs(lam.imag) > 1e-6:
            partner = np.conj(lam)          # fold-line locus: conjugate partner
        else:
            partner = -2.0 * Nn - lam       # off-line real locus: s-partner
        near_p = w[np.argsort(np.abs(w - partner))[:2]]
        gap_p = abs(near_p[0] - near_p[1])
        loci.append((qstar, lam, exp, gap_p))
    return loci


# ---- 1. the self-fold, entry-wise ------------------------------------------------------
print("=== 1. antiunitary self-fold on (1,3)@N=6, entry-wise ===")
Pf = bra_complement_perm(6, 1, 3)
for Jv in (0.8, 2.0, 4.0):
    L, _, _ = build_L(6, 1, 3, [Jv] * 5, [1.0] * 6)
    res = np.abs(L + Pf @ np.conj(L) @ Pf.T + 12.0 * np.eye(120)).max()
    gate(res == 0.0, f"J={Jv}: |L + Pf conj(L) Pf^T + 12 I|_max == 0.0 (is {res:.1e})")

# ---- 2. fold within sectors ------------------------------------------------------------
print("\n=== 2. the fold acts within each R-sector ===")
Pr = reflection_perm(6, 1, 3)
gate(np.abs(Pf @ Pr - Pr @ Pf).max() == 0.0, "[Pf, Pr] == 0.0")
L, _, _ = build_L(6, 1, 3, [2.0] * 5, [1.0] * 6)
gate(np.abs(L @ Pr - Pr @ L).max() == 0.0, "[L, Pr] == 0.0 at J=2 (uniform coupling)")

# ---- 3. the bipartite gauge and the closure criterion ----------------------------------
print("\n=== 3. bipartite gauge T L T == conj(L), character, and sector closure ===")


def sigma_sum(x, Nn):
    return sum(l for l in range(Nn) if x & bit(l, Nn))


def stagger_metric(Nn, p, q):
    pairs, _ = block_basis(Nn, p, q)
    return np.diag([(-1.0) ** (sigma_sum(a, Nn) + sigma_sum(b, Nn)) for a, b in pairs])


for (Nn, p, q) in ((6, 1, 3), (6, 1, 2), (6, 2, 3), (6, 0, 3), (4, 1, 2), (4, 1, 3)):
    Lb, _, _ = build_L(Nn, p, q, [1.7] * (Nn - 1), [1.0] * Nn)
    T = stagger_metric(Nn, p, q)
    Prb = reflection_perm(Nn, p, q)
    rG = np.abs(T @ Lb @ T - np.conj(Lb)).max()
    sgn = (-1.0) ** ((p + q) * (Nn - 1))
    rC = np.abs(T @ Prb - sgn * Prb @ T).max()
    gate(rG == 0.0 and rC == 0.0,
         f"({p},{q})@N={Nn}: T L T == conj(L) and T Pr == {int(sgn):+d} Pr T, both exactly 0.0")

print("  --- closure verdicts (nearest-neighbour residual, J=2.6/7.4) ---")
for (Nn, p, q), expect_closed in (((6, 1, 3), True), ((4, 1, 3), True), ((4, 1, 2), False),
                                  ((6, 1, 2), False), ((6, 2, 3), False), ((6, 0, 3), False)):
    for sec in ("even", "odd"):
        r = max(conj_closure_residual(eigs(Nn, p, q, Jv, sec)) for Jv in (2.6, 7.4))
        if expect_closed:
            gate(r < 1e-9, f"({p},{q})@N={Nn} R-{sec} conj-closed (residual {r:.1e})")
        else:
            gate(r > 0.1, f"({p},{q})@N={Nn} R-{sec} NOT conj-closed (residual {r:.1e})")

# the filtered AT clusters are semisimple, and the filter is idle in R-odd: at the R-even
# kappa=0 triples (lambda = -4 and -8) geometric = algebraic = 3 (third singular value of
# M - lambda I at floor, fourth O(1)); R-odd holds no exact degeneracy at all
print("  --- AT-cluster semisimplicity (R-even) and R-odd non-degeneracy ---")
Be13 = sector_basis(reflection_perm(6, 1, 3), even=True)
for Jv in (0.6, 2.0, 5.2):
    Lb, _, _ = build_L(6, 1, 3, [Jv] * 5, [1.0] * 6)
    M = Be13.T @ Lb @ Be13
    ok = True
    wM = np.linalg.eigvals(M)
    for lam0 in (-4.0, -8.0):
        sv = np.linalg.svd(M - lam0 * np.eye(60), compute_uv=False)
        alg = int(np.sum(np.abs(wM - lam0) < 1e-8))
        ok = ok and sv[-3] < 1e-10 and sv[-4] > 1e-3 and alg == 3
    gate(ok, f"J={Jv}: geometric = 3 (nullity) AND algebraic = 3 (eigenvalue count) at -4 and -8")
for Jv in (0.6, 2.0, 5.2):
    w = eigs(6, 1, 3, Jv, "odd")
    gate(len(collapse(w)) == 60, f"J={Jv}: R-odd carries no exact degeneracy (60 distinct)")

# ---- 4. sector disc-reality ------------------------------------------------------------
print("\n=== 4. sector disc arg-from-real (filtered, |delta| > 1e-8), J = 0.8/2.0/4.0 ===")
for (Nn, p, q, sec), real in (
    ((4, 1, 2, "even"), True), ((4, 1, 2, "odd"), True),
    ((6, 1, 2, "even"), False), ((6, 1, 2, "odd"), False),
    ((6, 1, 3, "even"), True), ((6, 1, 3, "odd"), True),
):
    args = [disc_arg_filtered(eigs(Nn, p, q, Jv, sec)) for Jv in (0.8, 2.0, 4.0)]
    if real:
        gate(all(a < 1e-9 for a in args),
             f"({p},{q})@N={Nn} R-{sec} disc REAL (args {', '.join(f'{a:.1e}' for a in args)})")
    else:
        gate(all(a > 0.02 for a in args) and max(args) > 0.1,
             f"({p},{q})@N={Nn} R-{sec} disc COMPLEX (args {', '.join(f'{a:.1e}' for a in args)})")

# ---- 5. detector validation on N=4 ground truth ----------------------------------------
print("\n=== 5. detector validation: (1,2)@N=4 R-even, the known quartet ===")
KNOWN = [0.460212, 0.854438, 0.857458, 1.738181]
loci4 = census(4, 1, 2, "even", 0.02, 3.0, 0.001)
dfc = sorted(q for q, _, e, _ in loci4 if abs(e - 0.5) < 0.05)
gate(len(dfc) == 4 and all(abs(a - b) < 5e-5 for a, b in zip(dfc, KNOWN)),
     f"exactly the four known loci found as sqrt-branch: {[round(v, 6) for v in dfc]}")

# ---- 6. the census ----------------------------------------------------------------------
print("\n=== 6. the (1,3)@N=6 census, octic q ===")
EXPECT = {
    "even": [(0.4630013, 0.3864821), (0.6400121, 1.2687899), (0.7450837, 1.4621529),
             (0.9184132, 7.7721940), (1.5109958, 0.9288088), (3.0922209, 1.5169893),
             (3.3536850, 0.0)],
    "odd": [(1.9567664, 7.6439946), (3.2608458, 25.1047962)],
}
channel_a = {}
for sec in ("even", "odd"):
    loci = census(6, 1, 3, sec, 0.02, 8.0, 0.001)
    defective = [(q, lam, e, gp) for q, lam, e, gp in loci if abs(e - 0.5) < 0.05]
    diabolic = [(q, lam, e) for q, lam, e, _ in loci if abs(e - 1.0) < 0.05]
    channel_a[sec] = [q for q, *_ in defective]
    print(f"  R-{sec}: {len(defective)} defective, {len(diabolic)} diabolic crossings")
    for q, lam, e, gp in defective:
        dbl = "s-doubled" if abs(lam.imag) <= 1e-6 else "conj-doubled"
        print(f"    q* = {q:.7f}: lambda = {lam.real:+.7f}{lam.imag:+.7f}i, "
              f"exponent {e:.3f}, partner-pair gap {gp:.1e} ({dbl})")
    for q, lam, e in diabolic:
        print(f"    q* = {q:.7f}: lambda = {lam.real:+.7f}{lam.imag:+.7f}i, "
              f"exponent {e:.4f}, diabolic crossing")
    gate(len(loci) == len(defective) + len(diabolic),
         f"R-{sec}: no unclassified survivor (every kept coalescence is 0.5 or 1.0)")
    exp = EXPECT[sec]
    gate(len(defective) == len(exp), f"R-{sec}: {len(exp)} defective loci")
    n_dia = {"even": 10, "odd": 12}[sec]
    gate(len(diabolic) == n_dia, f"R-{sec}: {n_dia} diabolic crossings inventoried")
    for (q, lam, e, gp), (q0, y0) in zip(sorted(defective), exp):
        gate(abs(q - q0) < 1e-5 and abs(abs(lam.imag) - y0) < 1e-5,
             f"R-{sec} q* = {q0}: position and lambda pinned")
        gate(gp < 1e-6, f"R-{sec} q* = {q0}: partner pair coalesces at the same q* (doubling)")
        onl = abs(lam.real + 6.0) < 1e-6
        gate(onl if y0 > 0 else not onl,
             f"R-{sec} q* = {q0}: lambda {'on' if y0 > 0 else 'off'} the fold line Re = -6")

print("\n=== 6c. channel B: the real-lambda species via the PT-break real-count scan ===")


def real_count(Nn, p, q, sector, qoct, tol=1e-9):
    w = eigs(Nn, p, q, 2.0 * qoct, sector)
    return int(np.sum(np.abs(w.imag) < tol))


def count_transitions(Nn, p, q, sector, q_lo, q_hi, step):
    qs = np.arange(q_lo, q_hi, step)
    out = []
    prev = None
    prev_q = None
    for qq in qs:
        c = real_count(Nn, p, q, sector, qq)
        if prev is not None and c != prev:
            lo, hi = prev_q, qq
            for _ in range(45):
                m = (lo + hi) / 2
                if real_count(Nn, p, q, sector, m) == prev:
                    lo = m
                else:
                    hi = m
            out.append(((lo + hi) / 2, prev, c))
        prev, prev_q = c, qq
    return out


REAL_EXPECT = {
    "even": [(0.7026159, 12, 8, (-6.5624176, -5.4375824)),
             (3.3536850, 8, 12, (-6.0787877, -5.9212123))],
    "odd": [(0.4635508, 12, 8, (-6.2680671, -5.7319329)),
            (0.8029541, 8, 4, (-6.6051577, -5.3948423)),
            (1.2273677, 4, 0, (-6.3760926, -5.6239074))],
}
channel_b = {}
for sec in ("even", "odd"):
    trans = count_transitions(6, 1, 3, sec, 0.02, 8.0, 0.001)
    exp = REAL_EXPECT[sec]
    gate(len(trans) == len(exp), f"R-{sec}: {len(exp)} real-count transitions")
    channel_b[sec] = []
    for (qstar, c0, c1), (q0, e0, e1, lams) in zip(trans, exp):
        gate(abs(qstar - q0) < 1e-5 and (c0, c1) == (e0, e1),
             f"R-{sec} q* = {q0}: transition {e0}->{e1} pinned")
        # the merging s-pair: real eigenvalues just on the higher-count side, nearest to
        # each expected lambda, checked to coalesce with sqrt exponent and to sum to -2N
        side = -1e-6 if c0 > c1 else +1e-6
        w = eigs(6, 1, 3, 2.0 * (qstar + side), sec)
        ok_pair = True
        lam_meas = []
        for lam in lams:
            near = w[np.argsort(np.abs(w - lam))[:2]]
            lam_meas.append(np.mean(near).real)
            ds = (1e-3, 1e-4, 1e-5)
            gs = []
            for d in ds:
                w2 = eigs(6, 1, 3, 2.0 * (qstar + np.sign(side) * d), sec)
                n2 = w2[np.argsort(np.abs(w2 - lam))[:2]]
                gs.append(abs(n2[0] - n2[1]))
            e = np.polyfit(np.log(ds), np.log(gs), 1)[0]
            ok_pair = ok_pair and abs(e - 0.5) < 0.05
        gate(ok_pair, f"R-{sec} q* = {q0}: both s-pair members sqrt-branch (defective)")
        # the -2N sum is gated on the MEASURED pair means at qstar+side, not the literals
        gate(abs(lam_meas[0] + lam_meas[1] + 12.0) < 1e-4,
             f"R-{sec} q* = {q0}: measured s-pair sums to -2N = -12 "
             f"({lam_meas[0]:.7f} + {lam_meas[1]:.7f})")
        channel_b[sec].append(qstar)
    # edges empty for the count channel
    for (elo, ehi, estep) in ((0.002, 0.02, 0.0002), (8.0, 25.0, 0.002)):
        etr = count_transitions(6, 1, 3, sec, elo, ehi, estep)
        gate(len(etr) == 0, f"R-{sec}: no real-count transitions in [{elo}, {ehi})")

# the union count: 8 R-even + 5 R-odd = 13 distinct loci (tolerance dedup, not rounding:
# the two channels bisect 3.3536850 to endpoints that straddle a 5-decimal rounding edge)
def dedup(values, tol=1e-4):
    out = []
    for v in sorted(values):
        if not out or v - out[-1] > tol:
            out.append(v)
    return out


union_even = dedup(channel_a["even"] + list(channel_b["even"]))
union_odd = dedup(channel_a["odd"] + list(channel_b["odd"]))
gate(len(union_even) == 8 and len(union_odd) == 5,
     f"UNION: 8 R-even + 5 R-odd = AT LEAST 13 defective loci ({len(union_even)}/{len(union_odd)})")

print("\n=== 6d. the edges are empty (channel A) ===")
for sec in ("even", "odd"):
    lo_loci = census(6, 1, 3, sec, 0.002, 0.02, 0.0002)
    hi_loci = census(6, 1, 3, sec, 8.0, 25.0, 0.002)
    gate(len(lo_loci) == 0, f"R-{sec}: no channel-A loci in [0.002, 0.02)")
    gate(len(hi_loci) == 0, f"R-{sec}: no channel-A loci in [8.0, 25.0)")

print("\n=== 7. sign-blindness vs the N=4 sign changes ===")


def sign_flips(Nn, p, q, sector, q_lo, q_hi, step):
    ef = lambda qq: eigs(Nn, p, q, 2.0 * qq, sector)

    def dsign(w):
        s = 0.0
        for i in range(len(w)):
            for j in range(i + 1, len(w)):
                d = w[i] - w[j]
                if abs(d) > 1e-8:
                    s += 2.0 * np.angle(d)
        return np.cos(s)

    qs = np.arange(q_lo, q_hi, step)
    signs = np.array([dsign(ef(qq)) for qq in qs])
    return [(qs[k] + qs[k + 1]) / 2 for k in range(len(qs) - 1) if signs[k] * signs[k + 1] < 0]


# step 0.001 here too: R-even holds loci 0.0028 apart (0.6998 vs 0.7026), and two odd-order
# zeros inside one cell would cancel their flips, hiding a violation from a coarser grid
for sec in ("even", "odd"):
    flips = sign_flips(6, 1, 3, sec, 0.02, 8.0, 0.001)
    gate(len(flips) == 0, f"(1,3)@N=6 R-{sec}: ZERO disc sign flips despite the loci (even order)")
# step 0.001 at N=4: the twin loci 0.854438 / 0.857458 are 0.003 apart, and two odd-order
# zeros inside one grid cell cancel their sign changes (a 0.01 grid reads 2 flips, not 4)
flips4 = sign_flips(4, 1, 2, "even", 0.02, 3.0, 0.001)
gate(len(flips4) == 4 and all(abs(f - k) < 0.01 for f, k in zip(sorted(flips4), KNOWN)),
     f"(1,2)@N=4 R-even: FOUR sign flips, one at each known locus ({[round(f, 3) for f in sorted(flips4)]})")

print()
if FAIL:
    print(f"{len(FAIL)} GATE(S) FAILED")
    sys.exit(1)
print("ALL GATES PASS")
