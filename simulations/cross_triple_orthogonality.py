r"""The cross-triple orthogonality: B(tau, sigma) = 0 for distinct vanishing triples.

This is the item left open by `y_zero_and_level_law.py` (2026-07-10d).  The twinning of
the Heff levels at resonant N needs, beyond Y = 0 on a mirror pair, that the block
B(tau, sigma) := (K26 W_tau)^T (K26 W_sigma) vanish for every pair of DISTINCT vanishing
mode triples.  The open column was eps(tau,sigma) = -1, where the site reflection makes
F(tau,sigma) = 0 an empty identity.

SCOPE.  The chain is the OPEN XY chain (H1 is bare adjacency, no ZZ diagonal in any
excitation sector), the hop amplitude is fixed to 1, and gamma only sets the rung spacing
(-2 gamma, -6 gamma) and is divided out.  Heff = P_ker K62 K26 P_ker carries NO energy
denominator, so its levels (18/n, 6/n) and (24/n, 8/n) are shape numbers, not decay rates.
Twinning is an isospectrality and is invariant under both scalings.

VOCABULARY.  As in `y_zero_and_level_law.py` and in `experiments/F89_SEED_EXISTENCE_REDUCTION.md`,
section "The twinning is a selection rule, and the level law is one line" -- with ONE index
shift, stated here because it is a real stumbling block:

    that file uses sites z = -1..N with u_k(z) = sin(k (z+1) pi / n);
    this file uses sites x = z + 1, i.e. x = 0..n with u_k(x) = sin(k x pi / n).

So b, c run 1..N here and 0..N-1 there.  Everything else is shared:

    n = N + 1,  theta = pi/n,  lam_k = 2 cos(k theta),
    D_tau(x,y,z)  the 3-magnon Slater determinant of the mode triple tau = {k1<k2<k3},
    G_tau(b,c)    = D_tau(b-1,b,c) + D_tau(b,b+1,c),
    U^{+-}(tau,sigma) = sum over {c > b} / {c < b} of G_tau G_sigma / (||D_tau|| ||D_sigma||),
    F(tau,sigma)  = U+ + U-          (the ledger's F, the full-square sum),
    eps(tau,sigma) = (-1)^{k(tau)+k(sigma)},  k(tau) = the MODE SUM k1+k2+k3.

Lemma 3 of that section says U+ = eps * U-.  So eps = +1 makes F the informative half, and
eps = -1 makes the antisymmetric half U+ - U- the informative one; F = 0 is then forced by
the reflection and says nothing.  We never abbreviate U+ - U- to a single letter in prose: `A` in the
ledger is the diagonal rung matrix of the pencil L(q) = A + qC, a different object.  The one
local variable that holds it is spelled `Araw`.

For a complementary mode pair {p<q} = tau \ {k}, `M_pq` is the ledger's object, with the
closed form (its step (G))
    M_pq(x) = (sin q theta - sin p theta) sin((p+q) x theta)
            - (sin q theta + sin p theta) sin((q-p) x theta),
equivalently, in the half-angle form used below,
    M_pq(x) = 2 cos((p+q)theta/2) sin((q-p)theta/2) sin((p+q) x theta)
            - 2 sin((p+q)theta/2) cos((q-p)theta/2) sin((q-p) x theta).
We write `psum` = p+q and `pdif` = q-p in the code, never `S` and `D`: `D_tau` is the Slater
determinant and `S` will be the constraint polynomial's coefficient.

--------------------------------------------------------------------------------------
WHAT THIS FILE ESTABLISHES, and at what grade.  Every step below is asserted separately.

(A) AN OPERATOR FORM.  [PROVED, gated]
    On the (1,2) block with the FREE hop H = H2 (x) I + I (x) H1 -- the bra passes THROUGH
    the kets, H2 keeps the ket-ket hardcore constraint, H1 is bare adjacency -- the
    Liouvillian's coherent hop is
        K = -J H J,   J = diag((-1)^{bra site}),
    a one-line identity: J commutes with H2 (x) I and anticommutes with I (x) H1.  J is the
    bra-side Jordan-Wigner sign that turns "ket forward, bra backward" into a free model.
    The totally antisymmetric Slater lift Psi_tau(a0,a1,beta) := D_tau(a0,a1,beta) obeys
    H Psi_tau = Lam_tau Psi_tau with Lam_tau = sum lam_k (free fermions), and it is
    supported entirely on the -6 rung by Pauli exclusion (the determinant has a repeated
    column when the bra meets a ket).  With nu := #{kets strictly left of the bra},
        g_tau := -H(nu Psi_tau) restricted to the -2 rung  =  G_tau,
    and, since nu = 1{c<b} there,
        U+ = <(1-nu) g_tau, (1-nu) g_sigma>,   U- = <nu g_tau, nu g_sigma>.
    (The "-6 support" of g_tau needs Lam_tau = 0; the -2 identity g = G holds for every
    triple.  The gate below distinguishes the two.)

(B) THE DATUM.  [GATED at N = 11, 14]
    Over ALL triples with eps = -1, U+ = 0 exactly when BOTH triples vanish.  One vanishing
    triple is not enough (|U+| up to 0.70).  In particular U+ is NOT proportional to
    Lam_tau * Lam_sigma, which that bucket would force to zero.

(C) A CORRECTION TO THE 2026-07-10d HANDOVER.  [GATED at N = 8, 10, 12, 14]
    At EVEN N every mirror pair has eps = -1 and Y = 0 holds; the handover read that as
    "the eps=-1 obstruction belongs to the proof, not the phenomenon".  It is stronger than
    that.  At even N a NON-vanishing mode-disjoint mirror pair has |U+| up to 0.42.
    Mode-disjointness alone, which carries the odd-N theorem, does NOT carry the even-N one.
    Even N consumes the lambda-vanishing; it is not the same lemma minus a sign.

(D) THE CLOSED FORM: n CANCELS.  [each intermediate PROVED + gated; the endpoint gated]
    Four elementary identities, each asserted separately below:
      (D1) Laplace along the c column:  G_tau(b,c) = sum_i (-1)^{i+1} u_{k_i}(c) M_i(b),
           with M_i := M_pq for {p,q} = tau minus k_i.
      (D2) the half-angle closed form of M_pq quoted above.
      (D3) the geometric sum   Theta_P(b) := sum_c sgn(c-b) cos(P c theta)
                = [1-(-1)^P]/2 - cot(P theta/2) sin(P b theta)      (P not = 0 mod 2n),
           with the SEPARATE branch Theta_0(b) = n - 2b.
      (D4) the triple-sine sum  sum_b s_x s_y s_P
                = (1/4)[cot((P+x-y)th/2) + cot((P-x+y)th/2)
                        - cot((P+x+y)th/2) - cot((P-x-y)th/2)]   when x+y+P is ODD,
                = 0                                              when x+y+P is EVEN.
           When x+y+P is odd all four arguments are odd, hence never 0 mod 2n: no pole.
    The parity in (D4) is exactly eps: x+y+P = k(tau)+k(sigma) mod 2.  So eps=+1 kills every
    term (this re-derives Lemma 4), and eps=-1 leaves cotangents.  Assembling, EVERY cotangent
    argument is an integer combination of the MODE ANGLES alone; the n cancels, and
        (U+ - U-) * (n/2)^3  =  cross_form(a; b),      a_i := k_i theta,  b_j := l_j theta,
    a pure trigonometric function of six angles.  Its mirror specialisation sigma = tau'
    (where G_{tau'} = (-1)^c G_tau turns each cotangent into a tangent) is the single-triple
    `mirror_form(a)`, equal to -4 * sum_{c>b} (-1)^c G_tau(b,c)^2 * (n/2)^3.
    (D3)'s P=0 branch is hit exactly when tau contains an ANTIPODAL pair k_i + k_j = n; that
    never happens for a vanishing triple at odd n, since lam_{k}+lam_{n-k} = 0 would force the
    third mode to n/2.  So the branch is excluded, at no cost.  This is Lemma 5 recovered.
    THIS IS A CLOSED FORM, NOT A MECHANISM.  The full-domain Gram F dies by mode orthogonality;
    the half-domain U+ has no such symmetry, and none is claimed.  What kills it is (E)/(F),
    a divisibility identity.  The sibling file says the same: "the missing argument is not a
    symmetry we have overlooked."

(E) mirror_form(a) = 0 whenever cos a1 + cos a2 + cos a3 = 0.  [PROVED over Q]
    With z_j = exp(i a_j) the object is rational; its numerator is divisible, exactly, by
    z3^2 + S z3 + 1 with S = z1 + 1/z1 + z2 + 1/z2.  The numerator is not identically zero,
    the division is exact, and the denominator is coprime to the constraint, so the vanishing
    holds on the whole variety, not merely off a lower-dimensional set.  The symbolic object
    is pinned numerically against `mirror_form`, so what is proved is what (D) gates -- writing
    it twice by hand would otherwise risk proving a theorem about a different expression.
    And `mirror_form` is pinned to `cross_form`: exactly, cross_form(a; mirror a) =
    mirror_form(a)/2.  So (E) proves the vanishing of the SAME six-angle object of (F),
    restricted to the positive-dimensional mirror subvariety b_j = pi - a_{4-j}.
    CONSEQUENCE, riding on the gated (D) assembly: Y = 0 for every mirror pair at EVEN N.
    That consequence is therefore CERTIFICATE grade, not proof grade; only the identity above
    is proved over Q.  It is the handover's cheap lead, closed at that grade.

(F) cross_form(a; b) = 0 on {sum cos a = 0} x {sum cos b = 0}.  [CERTIFIED in GF(p), NOT
    PROVED over Q.  This is the one load-bearing hole.]
    The GF(p) run and the float run share ONE implementation, `_cross_form_generic`, with the
    coefficient field swapped; and that implementation is pinned against `cross_form` over C.
    Certified by exact arithmetic in GF(p) at random points OF the variety, three primes,
    with controls.  The constraint z3^2 + S z3 + 1 = 0 has the two roots z3 and 1/z3; BOTH
    are tested, which is what "the quadratic divides the numerator" requires.  This is a
    Schwartz-Zippel certificate on the variety.  No degree bound is stated, so the confidence
    is overwhelming but not rigorously quantified.
    THE CONCRETE NEXT ATTACK: do not expand the six-variable numerator (sympy does not leave
    cancel()).  Work in the quotient ring instead: z3^2 = -S z3 - 1 and 1/z3 = -(z3 + S) make
    every monomial a linear form in z3, inverses come from the norm, and the whole object is
    a rank-2 module over Q(i)(z1,z2,w1,w2,w3).  Then repeat in w3.

(G) The mode-sharing eps = -1 pairs (6 of the 10 at N = 11) are covered by (F).  [the
    removability PROVED + gated; the on-variety vanishing GATED, and inherited from (F)]
    Two separate facts, tested separately below, because conflating them is easy:
      (G1) REMOVABILITY.  cross_form has a pole at a_i = b_j, and it is removable: the true
           A_raw at a shared mode is its limit.  Reason: for eps = -1, x and y have opposite
           parity, so <M_i^tau, M_j^sigma> = 0 and the n * <M,M> term of Theta_0 dies; and in
           sum_b b cos((x +- y) b theta) = n/2 - csc^2(...)/2 the two n/2 cancel.  The n-free
           remainder is exactly lim_{mu->0} cot(mu/2) * Xsum(mu).  Gated against a NONZERO
           reference, which forces the test pairs OFF the variety -- that is what removability
           means, and it is all this gate shows.
      (G2) ON-VARIETY VANISHING.  Approaching the coincidence while STAYING on sum cos a = 0,
           cross_form -> 0.  Gated separately.  This is the step that inherits (F)'s grade.

(H) The eps = +1 mode-sharing column ("F = 0 only measured" in the 10d table) closes by hand,
    uniformly in N, and subsumes the old TRIV x TRIV "by hand" case.  [PROVED, uniform in N]
    Two distinct vanishing triples sharing one mode k have complementary pairs with the SAME
    two-magnon energy,
        lam_p + lam_q = 4 cos(psum*th/2) cos(pdif*th/2) = -lam_k = lam_p' + lam_q'.
    Write red(m) for the reduction of a sine index to 1..n-1 (s_m = +- s_red(m), and s_m = 0
    when m = 0 or n).  By Lemma 4, F = (n/2)(-1)^{i+j} <M_{tau\k}, M_{sigma\k}>, four sine
    overlaps.  Each dies unless red coincides, and no coincidence is possible:
      psum = n  <=>  lam_k = 0  <=>  k = n/2  <=>  psum' = n.  So both triples are TRIV or
        neither.  If both: red(psum) = 0, the psum term of each M vanishes identically, and
        pdif = pdif' would force the same pair.  (This is the old by-hand case.)
      Otherwise cos(psum*th/2) != 0, and:
        red(psum)=red(psum'):  psum=psum' gives cos(pdif/2)=cos(pdif'/2), same pair.
                               psum+psum'=2n gives cos(psum/2)[cos(pdif/2)+cos(pdif'/2)]=0,
                               impossible (both cosines positive, pdif < n).
        red(pdif)=red(pdif'):  forces psum=psum', same pair.
        red(psum)=red(pdif'):  psum=pdif' gives pdif=psum', then p+q = q'-p' and q-p = p'+q'
                               add to q=q' and subtract to p=-p'.  Impossible.
                               2n-psum=pdif' gives p=p' and q+q'=2n > 2n-2.  Impossible.
        red(pdif)=red(psum'):  the same, with tau and sigma exchanged.
    Half-angles stay in (0,pi) since psum <= 2n-3, so cos is injective on them.  The
    same-energy hypothesis is load-bearing: drop it and arbitrary shared-mode pairs at N=23
    collide 10130 times in 60000.

(I) THE HINGE.  [PROVED, gated]
    "B(tau,sigma)=0 for tau != sigma" gives block-diagonality of Heff ONLY because the W-basis
    is orthonormal, W^T W = I: otherwise the projector W (W^T W)^{-1} W^T re-couples the
    triples.  It is orthonormal because each ordering sector assigns the bra a fixed rank and
    so sees every sorted triple exactly once, making the sector Gram <D_tau, D_sigma> = 0.
    Gated below, together with the resulting block-diagonality of Heff and the twinning.

--------------------------------------------------------------------------------------
TOGETHER.  B(tau,sigma) = 0 for every pair of distinct vanishing triples splits into four
cells, and closes at these grades:

    eps=+1, mode-disjoint   PROVED     Lemma 4 (the Laplace sum is empty)
    eps=+1, shared mode     PROVED     (H), uniform in N, subsumes TRIV x TRIV
    eps=-1, mode-disjoint   CERTIFIED  (D) + (F): the six-angle identity on the variety
    eps=-1, shared mode     CERTIFIED  (G): (F) plus a proved removable limit

So the twinning at every resonant N -- and with it the cell law's protection against its
cheapest kill -- is CERTIFIED, not proved.  The single load-bearing hole is a Q-level proof
that cross_form vanishes on the variety, plus, one level down, a symbolic proof of the (D)
endpoint (its four intermediates are proved; only their assembly is checked numerically).
At any FIXED N the block-diagonality is directly and exactly checkable, and is checked here.

    python simulations/cross_triple_orthogonality.py          # ~2 min
    python simulations/cross_triple_orthogonality.py --slow   # adds the sympy proof (E)

Authors: Thomas Wicht and Claude, 2026-07-10.
"""
import itertools
import math
import random
import sys

import numpy as np

TOL = 1e-9
PRIMES = [998244353, 1004535809, 469762049]     # each = 1 mod 4, so sqrt(-1) exists


# --------------------------------------------------------------------------- primitives
def umat(N):
    """u[k-1, x] = sin(k x pi / n), x = 0..n; the two walls are exact zeros."""
    n = N + 1
    U = np.array([[math.sin(k * x * math.pi / n) for x in range(n + 1)] for k in range(1, N + 1)])
    U[:, 0] = 0.0
    U[:, n] = 0.0
    return U


def lam(k, n):
    return 2.0 * math.cos(k * math.pi / n)


def slater(U, tau, x, y, z):
    return np.linalg.det(np.array([[U[k - 1, x], U[k - 1, y], U[k - 1, z]] for k in tau]))


def slater_norm_sq(U, tau, N):
    return sum(slater(U, tau, *t) ** 2 for t in itertools.combinations(range(1, N + 1), 3))


def Ggrid(U, tau, N):
    """G_tau(b,c) / ||D_tau||, for b, c = 1..N."""
    nrm = math.sqrt(slater_norm_sq(U, tau, N))
    G = np.zeros((N + 2, N + 2))
    for b in range(1, N + 1):
        for c in range(1, N + 1):
            G[b, c] = (slater(U, tau, b - 1, b, c) + slater(U, tau, b, b + 1, c)) / nrm
    return G


def Upm(Gt, Gs, N):
    up = sum(Gt[b, c] * Gs[b, c] for b in range(1, N + 1) for c in range(b + 1, N + 1))
    um = sum(Gt[b, c] * Gs[b, c] for b in range(1, N + 1) for c in range(1, b))
    return up, um


def vanishing(N):
    n = N + 1
    return [t for t in itertools.combinations(range(1, N + 1), 3)
            if abs(sum(lam(k, n) for k in t)) < 1e-10]


def mirror(t, n):
    return tuple(sorted(n - k for k in t))


def eps(t, s):
    return (-1) ** ((sum(t) + sum(s)) % 2)


def red(m, n):
    """s_m = +- s_red(m); red(m) = 0 means s_m vanishes identically."""
    m %= 2 * n
    return 0 if m in (0, n) else (m if m < n else 2 * n - m)


# ------------------------------------------------------------------- (A) the operator form
def block_space(N):
    kets = list(itertools.combinations(range(N), 2))
    kidx = {a: i for i, a in enumerate(kets)}
    basis = [(a, b) for a in kets for b in range(N)]

    H2 = np.zeros((len(kets), len(kets)))          # keeps the ket-ket hardcore constraint
    for a in kets:
        for pos in (0, 1):
            for step in (-1, +1):
                new = list(a)
                new[pos] += step
                if new[0] < 0 or new[1] > N - 1 or new[0] == new[1]:
                    continue
                H2[kidx[a], kidx[tuple(sorted(new))]] += 1.0
    H1 = np.zeros((N, N))                          # bare adjacency: the bra passes through
    for i in range(N - 1):
        H1[i, i + 1] = H1[i + 1, i] = 1.0

    H = np.kron(H2, np.eye(N)) + np.kron(np.eye(len(kets)), H1)
    K = -(np.kron(H2, np.eye(N)) - np.kron(np.eye(len(kets)), H1))
    J = np.diag([(-1.0) ** b for (a, b) in basis])
    nu = np.array([sum(1 for x in a if x < b) for (a, b) in basis], float)
    six = np.array([b not in a for (a, b) in basis])
    return basis, H, K, J, nu, six


def check_operator_form(N, pairs):
    U = umat(N)
    basis, H, K, J, nu, six = block_space(N)
    n = N + 1
    assert np.abs(K + J @ H @ J).max() < 1e-12, "K = -J H J failed"

    def Psi(t):
        return np.array([slater(U, t, a[0] + 1, a[1] + 1, b + 1) for (a, b) in basis])

    worst = {"eig": 0.0, "supp": 0.0, "rung": 0.0, "U": 0.0}
    for t in {x for p in pairs for x in p}:
        P = Psi(t)
        L = sum(lam(k, n) for k in t)
        worst["eig"] = max(worst["eig"], np.abs(H @ P - L * P).max())
        worst["supp"] = max(worst["supp"], np.abs(P[~six]).max())
        if abs(L) < 1e-10:                 # -H(nu Psi) is -2-confined only for Lam = 0
            worst["rung"] = max(worst["rung"], np.abs((-(H @ (nu * P)))[six]).max())

    def g_of(t):
        v = -(H @ (nu * Psi(t)))
        v[six] = 0.0                       # on the -2 rung, g = G_tau for EVERY triple
        return v

    for t, s in pairs:
        gt, gs = g_of(t), g_of(s)
        up_op = float(((1 - nu) * gt) @ ((1 - nu) * gs))
        um_op = float((nu * gt) @ (nu * gs))
        d2 = math.sqrt(slater_norm_sq(U, t, N) * slater_norm_sq(U, s, N))
        up_d, um_d = Upm(Ggrid(U, t, N), Ggrid(U, s, N), N)
        worst["U"] = max(worst["U"], abs(up_op / d2 - up_d), abs(um_op / d2 - um_d))

    for key, v in worst.items():
        assert v < 1e-8, f"(A) {key} = {v:.2e}"
    return worst


# -------------------------------------------------------- (B)(C) the two data points
def check_datum(N):
    U = umat(N)
    triples = list(itertools.combinations(range(1, N + 1), 3))
    G = {t: Ggrid(U, t, N) for t in triples}
    n = N + 1
    L = {t: sum(lam(k, n) for k in t) for t in triples}
    buckets = {"both": 0.0, "one": 0.0, "none": 0.0}
    counts = {"both": 0, "one": 0, "none": 0}
    for t, s in itertools.combinations(triples, 2):
        if eps(t, s) != -1:
            continue
        up, um = Upm(G[t], G[s], N)
        assert abs(up + um) < 1e-8, "Lemma 3 (U+ = eps U-) broken"
        v = [abs(L[t]) < 1e-10, abs(L[s]) < 1e-10]
        key = "both" if all(v) else ("one" if any(v) else "none")
        buckets[key] = max(buckets[key], abs(up))
        counts[key] += 1
    if counts["both"]:
        assert buckets["both"] < TOL, f"(B) both-vanish U+ = {buckets['both']:.2e}"
    if counts["one"]:
        assert buckets["one"] > 0.1, "(B) one-vanishing pairs should NOT be orthogonal"
    return counts, buckets


def check_even_N_needs_lambda(N):
    assert N % 2 == 0
    n = N + 1
    U = umat(N)
    van_worst, non_best, nvan, nnon = 0.0, 0.0, 0, 0
    for t in itertools.combinations(range(1, N + 1), 3):
        tp = mirror(t, n)
        if tp <= t:
            continue
        assert eps(t, tp) == -1, "even N must give eps = -1 on every mirror pair"
        up, _ = Upm(Ggrid(U, t, N), Ggrid(U, tp, N), N)
        if abs(sum(lam(k, n) for k in t)) < 1e-10:
            van_worst = max(van_worst, abs(up)); nvan += 1
        elif not (set(t) & set(tp)):
            non_best = max(non_best, abs(up)); nnon += 1
    assert van_worst < TOL, f"(C) vanishing mirror pair at even N: {van_worst:.2e}"
    assert non_best > 0.1, "(C) mode-disjointness alone should FAIL at even N"
    return nvan, van_worst, nnon, non_best


# --------------------------------------------- (D1)-(D4) the four intermediate identities
def check_expansion_lemmas(ns):
    """Each of the four identities the closed form is assembled from, separately."""
    w1 = w2 = w3 = w4 = 0.0
    for n in ns:
        N, th = n - 1, math.pi / n
        U = umat(N)
        s = lambda m, x: math.sin(m * x * th)

        # (D2) the half-angle closed form of M_pq
        for p in range(1, N + 1):
            for q in range(p + 1, N + 1):
                psum, pdif = p + q, q - p
                for x in range(1, N + 1):
                    nab = lambda y: U[p - 1, y] * U[q - 1, y + 1] - U[q - 1, y] * U[p - 1, y + 1]
                    lhs = nab(x - 1) + nab(x)
                    rhs = (2 * math.cos(psum * th / 2) * math.sin(pdif * th / 2) * s(psum, x)
                           - 2 * math.sin(psum * th / 2) * math.cos(pdif * th / 2) * s(pdif, x))
                    w2 = max(w2, abs(lhs - rhs))

        # (D1) Laplace along the c column
        for tau in itertools.combinations(range(1, N + 1), 3):
            for b in range(1, N + 1):
                for c in range(1, N + 1):
                    lhs = slater(U, tau, b - 1, b, c) + slater(U, tau, b, b + 1, c)
                    rhs = 0.0
                    for i in range(3):
                        p, q = sorted(k for k in tau if k != tau[i])
                        nab = lambda y: U[p - 1, y] * U[q - 1, y + 1] - U[q - 1, y] * U[p - 1, y + 1]
                        rhs += ((-1) ** i) * U[tau[i] - 1, c] * (nab(b - 1) + nab(b))
                    w1 = max(w1, abs(lhs - rhs))
            if n > 8:
                break                       # (D1) is O(n^3) per triple; one triple suffices per n

        # (D3) the geometric sum, both branches
        for P in range(0, 2 * n):
            for b in range(1, N + 1):
                lhs = sum((1 if c > b else (-1 if c < b else 0)) * math.cos(P * c * th)
                          for c in range(1, N + 1))
                if P % (2 * n) == 0:
                    rhs = n - 2 * b
                else:
                    rhs = (1 - (-1) ** P) / 2 - (1 / math.tan(P * th / 2)) * math.sin(P * b * th)
                w3 = max(w3, abs(lhs - rhs))

        # (D4) the triple-sine sum, both parity branches, and the no-pole claim
        cot = lambda y: 1.0 / math.tan(y)
        for x in range(1, 2 * n, 3):
            for y in range(1, 2 * n, 3):
                for P in range(1, 2 * n, 3):
                    lhs = sum(s(x, b) * s(y, b) * s(P, b) for b in range(1, N + 1))
                    if (x + y + P) % 2 == 0:
                        w4 = max(w4, abs(lhs))
                    else:
                        for arg in (P + x - y, P - x + y, P + x + y, P - x - y):
                            assert arg % (2 * n) != 0, "(D4) a cotangent pole appeared"
                        rhs = 0.25 * (cot((P + x - y) * th / 2) + cot((P - x + y) * th / 2)
                                      - cot((P + x + y) * th / 2) - cot((P - x - y) * th / 2))
                        w4 = max(w4, abs(lhs - rhs))
    for name, v in (("D1", w1), ("D2", w2), ("D3", w3), ("D4", w4)):
        assert v < 1e-8, f"({name}) = {v:.2e}"
    return w1, w2, w3, w4


def check_lemma5(Ns):
    """A vanishing triple at ODD n never contains an antipodal pair -> (D3)'s P=0 branch is moot."""
    bad = 0
    for N in Ns:
        n = N + 1
        for t in vanishing(N):
            if any(t[i] + t[j] == n for i in range(3) for j in range(i + 1, 3)):
                bad += 1
    assert bad == 0, f"(D/Lemma 5) {bad} vanishing triples with an antipodal pair at odd n"
    return bad


# ---------------------------------------------- the two pure-angle forms, and the gate
def _pieces(ang):
    """For each i: the complementary pair's sum/difference and the two M coefficients."""
    P = []
    for i in range(3):
        j, l = [t for t in range(3) if t != i]
        P.append(dict(psum=ang[j] + ang[l], pdif=ang[l] - ang[j],
                      alpha=math.sin(ang[l]) - math.sin(ang[j]),
                      beta=-(math.sin(ang[l]) + math.sin(ang[j])),
                      lap=(-1) ** (i + 1)))
    return P


def _cot(x):
    return 1.0 / math.tan(x)


def mirror_form(ang):
    """The single-triple (mirror / even-N) object, in tangent form.  No n in it."""
    P = _pieces(ang)
    T = math.tan

    def Xi(mu, xi, up):
        return 0.25 * (T((mu - xi + up) / 2) + T((mu + xi - up) / 2)
                       - T((mu - xi - up) / 2) - T((mu + xi + up) / 2))

    tot = 0.0
    for i in range(3):
        pi = P[i]
        for j in range(3):
            pj = P[j]
            def Xs(mu, pi=pi, pj=pj):
                return (pi['alpha'] * pj['alpha'] * Xi(mu, pi['psum'], pj['psum'])
                        + pi['alpha'] * pj['beta'] * Xi(mu, pi['psum'], pj['pdif'])
                        + pi['beta'] * pj['alpha'] * Xi(mu, pi['pdif'], pj['psum'])
                        + pi['beta'] * pj['beta'] * Xi(mu, pi['pdif'], pj['pdif']))
            mp, mm = ang[i] + ang[j], ang[i] - ang[j]
            tot += pi['lap'] * pj['lap'] * (-T(mp / 2) * Xs(mp)
                                            + (T(mm / 2) * Xs(mm) if abs(mm) > 1e-14 else 0.0))
    return tot


def cross_form(ang_tau, ang_sigma):
    """The full cross-triple object, in cotangent form.  Pure function of six angles."""
    Pa, Pb = _pieces(ang_tau), _pieces(ang_sigma)

    def Xh(mu, xi, up):
        return 0.25 * (_cot((mu + xi - up) / 2) + _cot((mu - xi + up) / 2)
                       - _cot((mu + xi + up) / 2) - _cot((mu - xi - up) / 2))

    tot = 0.0
    for i in range(3):
        pi = Pa[i]
        for j in range(3):
            pj = Pb[j]
            def Xs(mu, pi=pi, pj=pj):
                return (pi['alpha'] * pj['alpha'] * Xh(mu, pi['psum'], pj['psum'])
                        + pi['alpha'] * pj['beta'] * Xh(mu, pi['psum'], pj['pdif'])
                        + pi['beta'] * pj['alpha'] * Xh(mu, pi['pdif'], pj['psum'])
                        + pi['beta'] * pj['beta'] * Xh(mu, pi['pdif'], pj['pdif']))
            mu_p, mu_m = ang_tau[i] + ang_sigma[j], ang_tau[i] - ang_sigma[j]
            tot += ((-1) ** (i + j)) * (_cot(mu_p / 2) * Xs(mu_p) - _cot(mu_m / 2) * Xs(mu_m))
    return 0.5 * tot


def check_mirror_specialisation(Ns):
    """cross_form(a; mirror a) = (1/2) * mirror_form(a), exactly.

    The mirror triple sends k -> n-k, i.e. the angle a -> pi - a.  The factor 1/2 is
    bookkeeping: cross_form is U+ - U- = 2 U+ (since eps = -1 there), while mirror_form is
    -4 * sum_{c>b} (-1)^c G^2 * (n/2)^3 = -4 * (-U+) * ... .  This gate is what lets (E),
    which proves a statement about mirror_form over Q, be read as a statement about
    cross_form on a positive-dimensional subvariety.  Without it, (E) and (F) would be
    theorems about two objects nobody had identified.
    """
    worst, cnt = 0.0, 0
    for N in Ns:
        n, th = N + 1, math.pi / (N + 1)
        for t in itertools.combinations(range(1, N + 1), 3):
            tp = mirror(t, n)
            if set(t) & set(tp):
                continue                       # cross_form's derivation assumes mode-disjoint
            a, b = tuple(k * th for k in t), tuple(k * th for k in tp)
            try:
                cf, mf = cross_form(a, b), mirror_form(a)
            except ZeroDivisionError:
                continue
            if abs(mf) < 1e-9:
                continue
            worst = max(worst, abs(cf - 0.5 * mf) / abs(mf))
            cnt += 1
    assert cnt and worst < 1e-9, f"mirror_form is NOT cross_form's mirror specialisation: {worst:.2e}"
    return cnt, worst


def gate_angle_forms(Ns_even, Ns_any):
    w3 = w6 = big = 0.0
    n3 = n6 = 0
    for N in Ns_even:
        n, th, U = N + 1, math.pi / (N + 1), umat(N)
        for t in itertools.combinations(range(1, N + 1), 3):
            if any(t[i] + t[j] == n for i in range(3) for j in range(i + 1, 3)):
                continue                                     # (D3)'s P = 0 branch
            G = Ggrid(U, t, N)
            Q = sum(((-1) ** c) * G[b, c] ** 2
                    for b in range(1, N + 1) for c in range(b + 1, N + 1)) * (n / 2) ** 3
            w3 = max(w3, abs(Q - (-0.25) * mirror_form(tuple(k * th for k in t))))
            n3 += 1
    for N in Ns_any:
        n, th, U = N + 1, math.pi / (N + 1), umat(N)
        G = {t: Ggrid(U, t, N) for t in itertools.combinations(range(1, N + 1), 3)}
        for t, s in itertools.combinations(G, 2):
            if eps(t, s) != -1 or (set(t) & set(s)):
                continue
            up, um = Upm(G[t], G[s], N)
            f = cross_form(tuple(k * th for k in t), tuple(k * th for k in s))
            w6 = max(w6, abs((up - um) * (n / 2) ** 3 - f)); big = max(big, abs(f)); n6 += 1
    assert w3 < 1e-8 and w6 < 1e-8, f"(D) gate: {w3:.2e} / {w6:.2e}"
    assert big > 1.0, "(D) the six-angle form must take large values somewhere"
    return n3, w3, n6, w6, big


# ------------------------------------------------------------- (E) the exact Q proof
def prove_mirror_form_over_Q():
    import sympy as sp
    z1, z2, z3 = sp.symbols('z1 z2 z3')
    Z, I = {0: z1, 1: z2, 2: z3}, sp.I
    sin_a = lambda i: (Z[i] - 1 / Z[i]) / (2 * I)
    mono = lambda d: sp.prod([Z[j] ** c for j, c in d.items()])

    def tan_half(d):
        if all(c == 0 for c in d.values()):
            return sp.Integer(0)
        W = mono(d)
        return -I * (W - 1) / (W + 1)

    def add(*ds):
        o = {0: 0, 1: 0, 2: 0}
        for d in ds:
            for j, c in d.items():
                o[j] += c
        return o
    neg = lambda d: {j: -c for j, c in d.items()}
    e = lambda j: {0: 0, 1: 0, 2: 0, **{j: 1}}

    P = []
    for i in range(3):
        j, l = [t for t in range(3) if t != i]
        P.append(dict(psum=add(e(j), e(l)), pdif=add(e(l), neg(e(j))),
                      alpha=sin_a(l) - sin_a(j), beta=-(sin_a(l) + sin_a(j)),
                      lap=(-1) ** (i + 1)))

    def Xi(mu, xi, up):
        return sp.Rational(1, 4) * (tan_half(add(mu, neg(xi), up)) + tan_half(add(mu, xi, neg(up)))
                                    - tan_half(add(mu, neg(xi), neg(up))) - tan_half(add(mu, xi, up)))

    expr = sp.Integer(0)
    for i in range(3):
        pi = P[i]
        for j in range(3):
            pj = P[j]
            def Xs(mu, pi=pi, pj=pj):
                return (pi['alpha'] * pj['alpha'] * Xi(mu, pi['psum'], pj['psum'])
                        + pi['alpha'] * pj['beta'] * Xi(mu, pi['psum'], pj['pdif'])
                        + pi['beta'] * pj['alpha'] * Xi(mu, pi['pdif'], pj['psum'])
                        + pi['beta'] * pj['beta'] * Xi(mu, pi['pdif'], pj['pdif']))
            mp, mm = add(e(i), e(j)), add(e(i), neg(e(j)))
            expr += pi['lap'] * pj['lap'] * (-tan_half(mp) * Xs(mp) + tan_half(mm) * Xs(mm))

    # The symbolic expr and the numeric mirror_form are two independent transcriptions of the
    # same object.  Pin them to each other, or (E) would prove a theorem about a DIFFERENT
    # expression than the one (D) gates against the discrete Gram.
    ang = (0.7, 1.9, 2.5)
    got = complex(sp.N(expr.subs({z1: sp.exp(I * sp.Float(ang[0])),
                                  z2: sp.exp(I * sp.Float(ang[1])),
                                  z3: sp.exp(I * sp.Float(ang[2]))}), 30))
    assert abs(got.imag) < 1e-12 and abs(got.real - mirror_form(ang)) < 1e-9, (
        f"(E) the symbolic expression is not mirror_form: {got} vs {mirror_form(ang)}")

    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    num = sp.expand(num)
    assert sp.simplify(num) != 0, "(E) numerator identically zero -- the claim would be vacuous"
    assert sp.simplify(den) != 0
    csum = z1 + 1 / z1 + z2 + 1 / z2
    constraint = sp.expand(sp.numer(sp.cancel(z3 ** 2 + csum * z3 + 1)))
    quo, rem = sp.div(sp.Poly(num, z3), sp.Poly(constraint, z3))
    assert sp.simplify(rem.as_expr()) == 0, "(E) numerator does not reduce"
    assert sp.simplify(sp.expand(num - quo.as_expr() * constraint)) == 0, "(E) division not exact"
    # the denominator must not vanish on the constraint, else "= 0 on the variety" is weaker
    g = sp.gcd(sp.Poly(sp.expand(den), z3), sp.Poly(constraint, z3))
    assert sp.degree(g, z3) == 0, "(E) denominator shares a factor with the constraint"
    return sp.degree(num, z3), (len(num.args) if num.is_Add else 1)


# ------------------------------------------------------- (F) the exact GF(p) certificate
def _tonelli(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r


class _ComplexField:
    """cross_form's coefficient field over C, for pinning the z-form to the angle-form."""
    i = 1j
    zero = 0j
    lift = staticmethod(complex)
    add = staticmethod(lambda a, b: a + b)
    sub = staticmethod(lambda a, b: a - b)
    mul = staticmethod(lambda a, b: a * b)

    @staticmethod
    def div(a, b):
        if abs(b) < 1e-300:
            raise ZeroDivisionError
        return a / b


class _ModPField:
    """The same field ops in GF(p), p = 1 mod 4 so that i = sqrt(-1) exists."""

    def __init__(self, p):
        self.p = p
        self.i = _tonelli(p - 1, p)
        assert self.i is not None and self.i * self.i % p == p - 1
        self.zero = 0

    def lift(self, k):
        return int(k) % self.p

    def add(self, a, b):
        return (a + b) % self.p

    def sub(self, a, b):
        return (a - b) % self.p

    def mul(self, a, b):
        return a * b % self.p

    def div(self, a, b):
        if b % self.p == 0:
            raise ZeroDivisionError
        return a * pow(b, self.p - 2, self.p) % self.p


def _cross_form_generic(F, zs):
    """ONE transcription of cross_form in the z-variables, z_j = exp(i a_j).

    Used twice: over C (to pin it against the angle-form `cross_form`) and over GF(p)
    (for the certificate (F)).  Writing it twice by hand would risk certifying a
    different function than the one gated in (D).
    """
    one = F.lift(1)

    def ipow(x, k):
        if k < 0:
            return ipow(F.div(one, x), -k)
        r = one
        for _ in range(k):
            r = F.mul(r, x)
        return r

    def sin_v(j):
        return F.div(F.sub(zs[j], F.div(one, zs[j])), F.mul(F.lift(2), F.i))

    def mono(d):
        r = one
        for j, c in d.items():
            if c:
                r = F.mul(r, ipow(zs[j], c))
        return r

    def cot_half(d):
        Z = mono(d)
        return F.mul(F.i, F.div(F.add(Z, one), F.sub(Z, one)))

    add_d = lambda *ds: {j: sum(d.get(j, 0) for d in ds) for j in set().union(*[d.keys() for d in ds])}
    neg = lambda d: {j: -c for j, c in d.items()}
    e = lambda j: {j: 1}

    def pieces(off):
        P = []
        for i in range(3):
            j, l = [t for t in range(3) if t != i]
            P.append(dict(psum=add_d(e(off + j), e(off + l)), pdif=add_d(e(off + l), neg(e(off + j))),
                          alpha=F.sub(sin_v(off + l), sin_v(off + j)),
                          beta=F.sub(F.zero, F.add(sin_v(off + l), sin_v(off + j))),
                          lap=one if i % 2 == 0 else F.sub(F.zero, one)))
        return P

    Pa, Pb, quarter = pieces(0), pieces(3), F.div(one, F.lift(4))

    def Xh(mu, xi, up):
        r = F.sub(F.add(cot_half(add_d(mu, xi, neg(up))), cot_half(add_d(mu, neg(xi), up))),
                  F.add(cot_half(add_d(mu, xi, up)), cot_half(add_d(mu, neg(xi), neg(up)))))
        return F.mul(quarter, r)

    tot = F.zero
    for i in range(3):
        pi = Pa[i]
        for j in range(3):
            pj = Pb[j]

            def Xs(mu, pi=pi, pj=pj):
                r = F.mul(F.mul(pi['alpha'], pj['alpha']), Xh(mu, pi['psum'], pj['psum']))
                r = F.add(r, F.mul(F.mul(pi['alpha'], pj['beta']), Xh(mu, pi['psum'], pj['pdif'])))
                r = F.add(r, F.mul(F.mul(pi['beta'], pj['alpha']), Xh(mu, pi['pdif'], pj['psum'])))
                r = F.add(r, F.mul(F.mul(pi['beta'], pj['beta']), Xh(mu, pi['pdif'], pj['pdif'])))
                return r

            mu_p, mu_m = add_d(e(i), e(3 + j)), add_d(e(i), neg(e(3 + j)))
            term = F.sub(F.mul(cot_half(mu_p), Xs(mu_p)), F.mul(cot_half(mu_m), Xs(mu_m)))
            tot = F.add(tot, F.mul(pi['lap'], F.mul(pj['lap'], term)))
    return F.div(tot, F.lift(2))


def check_transcription(ang_tau, ang_sigma):
    """Pin the z-form to the angle-form: they must be the SAME function."""
    zs = [complex(math.cos(a), math.sin(a)) for a in tuple(ang_tau) + tuple(ang_sigma)]
    got = _cross_form_generic(_ComplexField, zs)
    want = cross_form(ang_tau, ang_sigma)
    assert abs(got.imag) < 1e-8 and abs(got.real - want) < 1e-6 * max(1.0, abs(want)), (
        f"the z-form and the angle-form disagree: {got} vs {want}")
    return abs(got.real - want)


def certify_cross_form_gfp(trials=40):
    """Both roots of z3^2 + S z3 + 1 = 0 are tested; they are z3 and 1/z3 (product 1)."""
    out = []
    for p in PRIMES:
        F = _ModPField(p)
        inv = lambda x: pow(x % p, p - 2, p)
        rng = random.Random(p % 100003)

        def both_roots(z1, z2):
            S = (z1 + inv(z1) + z2 + inv(z2)) % p
            d = _tonelli((S * S - 4) % p, p)
            if d is None:
                return None
            r1, r2 = (p - S + d) % p * inv(2) % p, (p - S - d) % p * inv(2) % p
            if r1 == 0 or r2 == 0:
                return None
            assert r1 * r2 % p == 1 % p, "the two roots must multiply to 1"
            return (r1, r2)

        good = bad = ctrl = 0
        while good < trials:
            try:
                z1, z2 = rng.randrange(2, p), rng.randrange(2, p)
                w1, w2 = rng.randrange(2, p), rng.randrange(2, p)
                rz, rw = both_roots(z1, z2), both_roots(w1, w2)
                if rz is None or rw is None:
                    continue
                vals = [_cross_form_generic(F, [z1, z2, z3, w1, w2, w3]) for z3 in rz for w3 in rw]
                ctrl_val = _cross_form_generic(F, [z1, z2, rz[0], w1, w2, rng.randrange(2, p)])
            except ZeroDivisionError:
                continue
            good += 1
            bad += sum(v != 0 for v in vals)
            ctrl += (ctrl_val != 0)
        assert bad == 0, f"(F) cross_form != 0 at {bad}/{4 * good} variety points mod {p}"
        assert ctrl > 0.9 * good, f"(F) control too weak mod {p}"
        out.append((p, 4 * good, bad, ctrl, good))
    return out


# ---------------------------- (G1) removability  and  (G2) the on-variety approach
def check_removability(N, how_many=6):
    """cross_form's pole at a_i = b_j is removable: its limit is the true A_raw.

    Needs a NONZERO reference, so the pairs here are deliberately OFF the variety.
    This gate shows removability ONLY.  The on-variety vanishing is (G2).
    """
    n, th, U = N + 1, math.pi / (N + 1), umat(N)
    G = {t: Ggrid(U, t, N) for t in itertools.combinations(range(1, N + 1), 3)}
    worst, shown = 0.0, 0
    for t, s in itertools.combinations(G, 2):
        sh = sorted(set(t) & set(s))
        if eps(t, s) != -1 or len(sh) != 1:
            continue
        up, um = Upm(G[t], G[s], N)
        Araw = (up - um) * (n / 2) ** 3
        if abs(Araw) < 1.0:
            continue
        f = cross_form(tuple(k * th + (1e-6 if k == sh[0] else 0) for k in t),
                       tuple(k * th for k in s))
        worst = max(worst, abs(Araw - f) / max(1.0, abs(Araw)))
        shown += 1
        if shown == how_many:
            break
    assert shown and worst < 1e-4, f"(G1) removability: {worst:.2e} over {shown} pairs"
    return shown, worst


def check_on_variety_approach(N, deltas=(1e-2, 1e-4, 1e-6)):
    """Approach the mode coincidence WHILE STAYING on sum cos a = 0.  cross_form -> 0."""
    n, th = N + 1, math.pi / (N + 1)
    van = vanishing(N)
    worst, shown = 0.0, 0
    for t, s in itertools.combinations(van, 2):
        sh = sorted(set(t) & set(s))
        if eps(t, s) != -1 or len(sh) != 1:
            continue
        i = t.index(sh[0]); m = (i + 1) % 3; o = 3 - i - m
        for d in deltas:
            ang = [k * th for k in t]
            ang[i] += d
            tgt = -(math.cos(ang[i]) + math.cos(ang[o]))
            if abs(tgt) > 1:
                break
            ang[m] = math.acos(tgt)
            assert abs(sum(math.cos(v) for v in ang)) < 1e-12, "(G2) left the variety"
            worst = max(worst, abs(cross_form(tuple(ang), tuple(k * th for k in s))))
        shown += 1
    assert shown and worst < 1e-8, f"(G2) on-variety approach: {worst:.2e}"
    return shown, worst


# ------------------------------------ (H) the eps=+1 shared-mode column, uniform in N
def check_no_sine_coincidence(Ns):
    """TRIV x TRIV INCLUDED: red(psum) = 0 there, so that term of M vanishes identically."""
    total = hits = 0
    for N in Ns:
        n, th = N + 1, math.pi / (N + 1)
        for t, s in itertools.combinations(vanishing(N), 2):
            shq = set(t) & set(s)
            assert len(shq) != 2, "two distinct vanishing triples cannot share two modes"
            if eps(t, s) != 1 or not shq:
                continue
            k = shq.pop()
            total += 1
            p, q = sorted(x for x in t if x != k)
            p2, q2 = sorted(x for x in s if x != k)
            for (P, Q) in ((p, q), (p2, q2)):
                e = 4 * math.cos((P + Q) * th / 2) * math.cos((Q - P) * th / 2)
                assert abs(e + lam(k, n)) < 1e-9, "(H) the same-two-magnon-energy identity fails"
                assert 0 < P + Q <= 2 * n - 3 and 0 < Q - P < n, "(H) half-angle out of (0,pi)"
            for x in (p + q, q - p):
                for y in (p2 + q2, q2 - p2):
                    if red(x, n) != 0 and red(x, n) == red(y, n):
                        hits += 1
    assert hits == 0, f"(H) {hits} sine-mode coincidences -- the by-hand lemma fails"
    return total


def check_same_energy_is_load_bearing(N, cap=60000):
    """Drop the same-energy hypothesis and arbitrary shared-mode pairs DO collide."""
    n = N + 1
    coll = tot = 0
    for t in itertools.combinations(range(1, N + 1), 3):
        for s in itertools.combinations(range(1, N + 1), 3):
            if t >= s or len(set(t) & set(s)) != 1:
                continue
            tot += 1
            if tot > cap:
                return coll, tot - 1
            k = (set(t) & set(s)).pop()
            p, q = sorted(x for x in t if x != k)
            p2, q2 = sorted(x for x in s if x != k)
            if any(red(x, n) != 0 and red(x, n) == red(y, n)
                   for x in (p + q, q - p) for y in (p2 + q2, q2 - p2)):
                coll += 1
    return coll, tot


# ------------------------------------------ (I) the hinge: W^T W = I, and what it buys
def check_hinge(N):
    n, U = N + 1, umat(N)
    van = vanishing(N)
    six = [(a, b) for a in itertools.combinations(range(N), 2) for b in range(N) if b not in a]
    two = [(a, b) for a in itertools.combinations(range(N), 2) for b in range(N) if b in a]
    cols = []
    for t in van:
        nrm = math.sqrt(slater_norm_sq(U, t, N))
        for srank in range(3):
            v = np.zeros(len(six))
            for i, (a, b) in enumerate(six):
                z = sorted(a + (b,))
                if z.index(b) == srank:
                    v[i] = ((-1.0) ** b) * slater(U, t, z[0] + 1, z[1] + 1, z[2] + 1) / nrm
            cols.append(v)
    W = np.array(cols).T
    gram = W.T @ W
    orth = np.abs(gram - np.eye(W.shape[1])).max()
    assert orth < 1e-10, f"(I) W^T W != I: {orth:.2e}"

    _, H, K, _, _, sixmask = block_space(N)
    idx6 = np.where(sixmask)[0]
    idx2 = np.where(~sixmask)[0]
    K26 = K[np.ix_(idx2, idx6)]
    # our `six`/`two` orderings must match block_space's basis order
    KW = K26 @ W
    Heff = KW.T @ KW
    cross = 0.0
    for i in range(len(van)):
        for j in range(len(van)):
            if i == j:
                continue
            cross = max(cross, np.abs(Heff[3 * i:3 * i + 3, 3 * j:3 * j + 3]).max())
    assert cross < 1e-8, f"(I) Heff is NOT block-diagonal over the triples: {cross:.2e}"
    return W.shape[1], orth, cross


# -------------------------------------------------------------------------------- main
def main():
    slow = "--slow" in sys.argv
    print("The cross-triple orthogonality, step by separately asserted step.\n")

    w = check_operator_form(9, [((1, 5, 9), (2, 5, 8)), ((3, 5, 7), (4, 5, 6))])
    w2 = check_operator_form(11, [((1, 7, 9), (3, 5, 11)), ((1, 6, 11), (1, 7, 9))])
    print(f"(A) K = -J H J exact | H Psi = Lam Psi {max(w['eig'], w2['eig']):.1e}"
          f" | Psi on -6 only {max(w['supp'], w2['supp']):.1e}"
          f" | g = -H(nu Psi) on -2 only (Lam=0) {max(w['rung'], w2['rung']):.1e}"
          f" | U+- as half-Grams {max(w['U'], w2['U']):.1e}")

    for N in (11, 14):
        c, b = check_datum(N)
        print(f"(B) N={N:2d} eps=-1, ALL triples | both vanish {c['both']:3d} max|U+| {b['both']:.1e}"
              f" | one vanishes {c['one']:4d} max|U+| {b['one']:.3f}"
              f" | neither {c['none']:5d} max|U+| {b['none']:.3f}")
    print("    => U+ = 0 iff BOTH vanish; so U+ is not proportional to Lam_tau * Lam_sigma.")

    for N in (8, 10, 12, 14):
        nv, vw, nn, nb = check_even_N_needs_lambda(N)
        print(f"(C) N={N:2d} (even) mirror pairs | vanishing {nv} max|U+| {vw:.1e}"
              f" | NON-vanishing mode-disjoint {nn:3d} max|U+| {nb:.3f}")
    print("    => at even N, mode-disjointness alone does NOT give Y = 0.  Lambda is consumed.")

    d1, d2, d3, d4 = check_expansion_lemmas(range(5, 13))
    check_lemma5(range(8, 31, 2))
    print(f"(D) the four intermediates, n = 5..12 | (D1) Laplace {d1:.1e} | (D2) M_pq {d2:.1e}"
          f" | (D3) geometric sum + P=0 branch {d3:.1e} | (D4) triple sine + parity {d4:.1e}")
    n3, w3, n6, w6, big = gate_angle_forms((8, 10, 12, 14), (9, 11, 12, 13))
    print(f"    assembly | mirror_form on {n3} triples, err {w3:.1e}"
          f" | cross_form on {n6} pairs, err {w6:.1e} (values up to {big:.1f}).  n cancels.")
    nm, wm = check_mirror_specialisation((8, 10, 12, 14))
    rng = random.Random(5)
    wt = max(check_transcription(tuple(rng.uniform(0.3, 2.8) for _ in range(3)),
                                 tuple(rng.uniform(0.3, 2.8) for _ in range(3))) for _ in range(20))
    print(f"    pins | cross_form(a; mirror a) = mirror_form(a)/2 on {nm} pairs, rel err {wm:.1e}"
          f" | the z-form equals the angle-form, err {wt:.1e}")

    if slow:
        deg, terms = prove_mirror_form_over_Q()
        print(f"(E) mirror_form = 0 on cos a1+cos a2+cos a3 = 0: PROVED over Q"
              f" (numerator degree {deg} in z3, {terms} terms; exact division; den coprime).")
    else:
        print("(E) skipped (use --slow for the sympy proof over Q).")

    for p, pts, bad, ctrl, spec in certify_cross_form_gfp():
        print(f"(F) p = {p:>10} | cross_form vanishes at {pts - bad}/{pts} variety points"
              f" ({spec} specialisations x both roots of each quadratic)"
              f" | control nonzero {ctrl}/{spec}")
    print("    => CERTIFIED in GF(p), not proved over Q.  This is the one load-bearing hole.")

    for N in (11, 17):
        shown, worst = check_removability(N)
        shown2, worst2 = check_on_variety_approach(N)
        print(f"(G) N={N:2d} | (G1) removability, {shown} OFF-variety pairs with A != 0,"
              f" rel err {worst:.1e} | (G2) on-variety approach, {shown2} pairs,"
              f" max|cross_form| {worst2:.1e}")

    total = check_no_sine_coincidence(range(11, 42, 2))
    coll, tot = check_same_energy_is_load_bearing(23)
    print(f"(H) eps=+1 shared mode, odd N = 11..41: {total} pairs (TRIVxTRIV included),"
          f" 0 sine-mode coincidences")
    print(f"    control: WITHOUT the same-energy hypothesis, {coll}/{tot} arbitrary shared-mode"
          f" pairs at N=23 do collide => the hypothesis is load-bearing.")

    for N in (11, 17):
        ncol, orth, cross = check_hinge(N)
        print(f"(I) N={N:2d} | W is {ncol} columns, max|W^T W - I| = {orth:.1e}"
              f" | Heff cross-triple blocks {cross:.1e} => block-diagonal")

    print("\nAll steps verified.  B(tau,sigma) = 0 for every pair of DISTINCT vanishing triples:")
    print("  eps = +1, mode-disjoint   PROVED     Lemma 4 (the Laplace sum is empty)")
    print("  eps = +1, shared mode     PROVED     (H), uniform in N, subsumes TRIV x TRIV")
    print("  eps = -1, mode-disjoint   CERTIFIED  (D) + (F), the six-angle identity")
    print("  eps = -1, shared mode     CERTIFIED  (G), (F) plus a proved removable limit")
    print("With (I), Heff is block-diagonal over the triples, so the [[X,Y],[Y,X]] compression")
    print("is legitimate and the full-spectrum twinning follows -- at CERTIFICATE grade, not")
    print("proof grade.  THE ONE HOLE: a Q-level proof that cross_form vanishes on the variety")
    print("(and, one level down, a symbolic proof of the (D) assembly).  See (F) for the route.")


if __name__ == "__main__":
    main()
