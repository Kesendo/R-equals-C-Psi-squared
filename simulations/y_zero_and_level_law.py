"""The twinning orthogonality Y = 0 and the coupled-level law, step by step.

Companion to `experiments/F89_SEED_EXISTENCE_REDUCTION.md`, section
"The twinning is a selection rule, and the level law is one line (2026-07-10d)".
Sibling of `simulations/resonant_n_twinning.py`, which measures the two facts this
script's steps prove. Nothing here is a fit: every step is an identity, and the
script asserts each one separately rather than only the conclusion.

The objects (vocabulary in the F89 doc, section "The resonance criterion, closed"):
K is the real symmetric hop on the (1,2) coherence block, K_ab := P_{-a} K P_{-b};
n = N + 1; u_k(z) = sin(k(z+1)pi/n) the one-magnon modes; lam_k = 2 cos(k pi/n);
D_tau(x,y,z) = det[u_k(.)] over the mode triple tau the 3-magnon Slater determinant;
w_{tau,s} its lift to the -6 rung in ordering sector s, gauged by (-1)^{z_bra};
tau' = {n-k : k in tau} the mirror triple. On a mirror pair, Heff = P_ker K_62 K_26 P_ker
reads [[X, Y], [Y, X]], so the E-block X+Y and the O-block X-Y are twinned iff Y = 0.

The proof, and the step each check pins. B(tau,sigma) := (K26 W_tau)^T (K26 W_sigma) is the
block between two triples; X = B(tau,tau), Y = B(tau,tau'); k(tau) := k1+k2+k3 is the MODE SUM
(a number, not the hop matrix K).

  (A) K_26 w_{tau,s} (b,c) = -(-1)^b G_tau(b,c) 1{c in Omega_s} / ||D_tau||, on the -2 rung
      state (bra b, other ket site c), with G_tau(b,c) = D(b-1,b,c) + D(b,b+1,c) and
      Omega_0 = {c > b}, Omega_2 = {c < b}, Omega_1 = everything.  The sector indicator
      collapses because c always takes rank 0 (c < b) or rank 2 (c > b).
      Consequence, from the six region intersections alone: with
      U± := sum over {c >< b} of G_tau G_sigma / (||D_tau|| ||D_sigma||),
      B(tau,sigma) = [[U+, U+, 0], [U+, U+ + U-, U-], [0, U-, U-]] for ANY two triples.
  (B) G_{tau'}(b,c) = +(-1)^c G_tau(b,c).  Both terms of G carry the same site-sum parity,
      since (b-1)+b and b+(b+1) are both odd + 2b.  (The + is the doc's - in
      v_{tau'} = -(-1)^{sum z} v_tau composed with that odd parity.)  NOT needed for the
      theorem; it says Y is a scalar multiple of X.
  (C) The site reflection z -> N-1-z gives G_tau(rho b, rho c) = (-1)^{k(tau)} G_tau(b,c) and
      swaps {c > b} with {c < b}.  Hence for ANY two triples
      U+ = eps * U-,  eps := (-1)^{k(tau)+k(sigma)}.
      sigma = tau gives eps = +1 at every N, so X = a * [[1,1,0],[1,2,1],[0,1,1]].
      sigma = tau' gives eps = (-1)^{3n}, which is +1 exactly when n is even, i.e. N odd.
      THAT one sign is the whole content of the "N odd" assumption.
  (D) sum_{c=0}^{N-1} u_p(c) u_q(c) = (n/2) delta_{p,q}   (sine orthogonality, interior points).
  (E) A vanishing triple meets its own mirror (tau ∩ tau' != {}) IFF it is self-mirror, i.e.
      IFF it is TRIV.  So the extras, the triples that HAVE a mirror partner, are mode-disjoint
      from it.
  (F) Y = 0.  Laplace-expand G along its c column: G(b,c) = sum_i (-1)^{i+1} u_{k_i}(c) M_i(b).
      Because G(b,b) = 0, the FULL-SQUARE sum factorizes:
      F(tau,sigma) := sum_{b,c} G_tau G_sigma = (n/2) sum_{k_i = l_j} (-1)^{i+j} sum_b M_i M_j,
      an EMPTY sum when the triples share no mode (D).  For a mirror pair (E) gives
      mode-disjointness and (C) gives eps = +1, so U+ = U- and U+ + U- = F = 0, so both vanish
      and every entry of Y with them.  Uniform in N; no per-N check.
      What this does NOT use: the vanishing lam_k1 + lam_k2 + lam_k3 = 0 enters ONLY through (E).  Y = 0 holds
      for ANY triple mode-disjoint from its mirror, resonant or not (checked below, together
      with a witness showing the hypothesis is not removable).  The twinning is a
      mode-reflection selection rule, not an arithmetic accident of resonance.

  The level law falls out of the same computation:
  (G) M_pq(b) = m_pq(b-1,b) + m_pq(b,b+1) with m_pq the 2x2 Slater det of the mode pair.
      Writing zeta = b+1, the two-mode piece telescopes to
      M_pq(zeta) = (sin(q pi/n) - sin(p pi/n)) sin((p+q) zeta pi/n)
                 + (sin(q pi/n) + sin(p pi/n)) sin((p-q) zeta pi/n),
      and sine orthogonality gives  ||M_pq||^2 := sum_b M_pq(b)^2
                                              = n (sin^2(p pi/n) + sin^2(q pi/n))  exactly.
  (H) ||D_tau||^2 = (n/2)^3 for every triple of distinct modes (Cauchy-Binet).
  (I) X = a * [[1,1,0],[1,2,1],[0,1,1]] with a = U+(tau,tau), whose spectrum is {3a, a, 0}:
      the 3:1 ratio and the uncoupled level are structural, and summing (G) over the three
      mode pairs of tau, then dividing by (H),
      a = (4/n) sum_{k in tau} sin^2(k pi/n) = (12 - sum lam_k^2) / n,  for EVERY triple.
      For ROT3 sum lam^2 = 6 (the squares of a rotated cube-root triple are again one),
      for PENT sum lam^2 = 4 (the golden ratio cancels exactly), giving the measured
      coupled levels (18/n, 6/n) and (24/n, 8/n).

  (J) SCOPE, and one open item.  Y = 0 is the vanishing of the cross block of a mirror pair.
      The full-spectrum twinning that `resonant_n_twinning.py` measures (spec Heff|_E inside
      spec Heff|_O over the WHOLE kernel) needs more: the pair's 6-dimensional space must be
      Heff-invariant, i.e. B(tau, sigma) must vanish for DISTINCT triples tau != sigma as well,
      TRIV triples included.  That is measured to 1e-15 here and is NOT proved everywhere.
      Four populations, counted separately, never silently merged:
        lemma4  : no shared mode and reflection sign +1 -> (C) + (F) run verbatim.
        byhand  : two distinct TRIV triples             -> eps = +1, shared mode is exactly n/2,
                  and the two-mode overlap dies by sine orthogonality.  Proved.
        measured: reflection sign +1 with a shared mode, not two TRIV -> F = 0 is measured.
        open    : reflection sign -1                    -> (C) makes F = 0 vacuous.
      The doc's phrase "Heff restricted to the pair" tacitly assumed the invariance.
      Two facts that locate the gap, both checked below: the theorem's real hypothesis is
      mode-disjointness (it holds for NON-vanishing triples too, and fails when modes are
      shared), and at EVEN N every mirror pair has reflection sign -1 yet Y = 0 anyway.
      So the sign obstruction belongs to this proof, not to the phenomenon.

Run: python simulations/y_zero_and_level_law.py         (N = 11, 17; measured 3.4 s)
     python simulations/y_zero_and_level_law.py 23 29   (adds N = 23 and N = 29, the first
                                                         PENT pair; measured 30 s)

Authors: Thomas Wicht and Claude, 2026-07-10.
"""
import sys
from itertools import combinations

import numpy as np

sys.path.insert(0, "simulations")
from resonant_n_twinning import blocks, z3_exact  # noqa: E402

TOL = 1e-9


def sk(k, n):
    return np.sin(k * np.pi / n)


def lam(k, n):
    return 2.0 * np.cos(k * np.pi / n)


def modes(N):
    """u_k(z) for k = 1..N and z = -1..N; the walls z = -1 and z = N vanish identically."""
    n = N + 1
    return np.array([[np.sin(k * (j + 1) * np.pi / n) for j in range(-1, N + 1)]
                     for k in range(1, N + 1)])


def slater(U, tau, x, y, z):
    """D_tau(x,y,z), antisymmetric, zero on the walls and on coincidences."""
    return float(np.linalg.det(np.array([[U[k - 1, w + 1] for w in (x, y, z)] for k in tau])))


def G(U, tau, b, c):
    return slater(U, tau, b - 1, b, c) + slater(U, tau, b, b + 1, c)


def slater_norm(U, tau, N):
    return float(np.linalg.norm([slater(U, tau, *z) for z in combinations(range(N), 3)]))


def lift(U, tau, s, N, nrm):
    """w_{tau,s}: the gauged Slater determinant on the -6 rung, ordering sector s."""
    basis = [(a, b) for a in combinations(range(N), 2) for b in combinations(range(N), 1)]
    b6 = [(a, b[0]) for (a, b) in basis if b[0] not in a]
    w = np.zeros(len(b6))
    for i, (a, b) in enumerate(b6):
        z = tuple(sorted(a + (b,)))
        if list(z).index(b) == s:
            w[i] = ((-1.0) ** b) * slater(U, tau, *z) / nrm
    return w


# --------------------------------------------------------------- the analytic steps

def check_D_and_E(N):
    """(D) sine orthogonality: the plain form Lemma 4 uses, and the signed form that is
    Lemma 2's face of it; (E) a vanishing triple meets its own mirror iff it is self-mirror."""
    n, U = N + 1, modes(N)
    for p in range(1, N + 1):
        for q in range(1, N + 1):
            plain = sum(U[p - 1, c + 1] * U[q - 1, c + 1] for c in range(N))
            signed = sum((-1.0) ** c * U[p - 1, c + 1] * U[q - 1, c + 1] for c in range(N))
            assert abs(plain - ((n / 2) if p == q else 0.0)) < TOL, \
                f"(D) N={N} plain ({p},{q}) = {plain}"
            assert abs(signed - ((n / 2) if p + q == n else 0.0)) < TOL, \
                f"(D) N={N} signed ({p},{q}) = {signed}"

    (_, _, _), triples = z3_exact(N)
    for kind, tt in enumerate(triples):
        for tau in tt:
            meets = any(x + y == n for x, y in combinations(tau, 2)) or any(2 * k == n for k in tau)
            selfm = tuple(sorted(n - k for k in tau)) == tuple(tau)
            assert meets == selfm, f"(E) N={N} tau={tau}: meets-mirror {meets} != self-mirror {selfm}"
            assert selfm == (kind == 0), f"(E) N={N} tau={tau}: self-mirror is not exactly TRIV"


def check_G_and_H(N):
    """(G) ||M_pq||^2 = n(sin^2(p pi/n) + sin^2(q pi/n)) exactly; (H) ||D_tau||^2 = (n/2)^3."""
    n, U = N + 1, modes(N)
    for p, q in combinations(range(1, N + 1), 2):
        M = [(U[p - 1, b] * U[q - 1, b + 1] - U[p - 1, b + 1] * U[q - 1, b])
             + (U[p - 1, b + 1] * U[q - 1, b + 2] - U[p - 1, b + 2] * U[q - 1, b + 1])
             for b in range(N)]  # index shift: column b <-> site b-1
        R = float(np.dot(M, M))
        want = n * (sk(p, n) ** 2 + sk(q, n) ** 2)
        assert abs(R - want) < 1e-8 * max(1.0, want), f"(G) N={N} ({p},{q}): {R} != {want}"

    for tau in list(combinations(range(1, N + 1), 3))[:30]:
        nv2 = slater_norm(U, tau, N) ** 2
        assert abs(nv2 - (n / 2) ** 3) < 1e-6 * (n / 2) ** 3, f"(H) N={N} {tau}: {nv2} != {(n/2)**3}"


def check_level_law(N):
    """(I, arithmetic half) a*n = 12 - sum lam^2 = 6 on ROT3, 8 on PENT."""
    n = N + 1
    (_, _, _), triples = z3_exact(N)
    out = []
    for kind, want in ((1, 6.0), (2, 8.0)):
        for tau in triples[kind]:
            an = 12.0 - sum(lam(k, n) ** 2 for k in tau)
            an2 = 4.0 * sum(sk(k, n) ** 2 for k in tau)
            assert abs(an - want) < TOL and abs(an2 - want) < TOL, \
                f"(I) N={N} {tau}: a*n = {an} / {an2} != {want}"
        if triples[kind]:
            out.append(f"{['', 'ROT3', 'PENT'][kind]} a*n={want:.0f} ({len(triples[kind])})")
    return out


# --------------------------------------------------------------- the matrix steps

def check_matrix_steps(N):
    """(A) the elementwise formula, (B) the mirror law, (C) T+ = T- = 0, (F) Y = 0, (I) X's shape."""
    n, U = N + 1, modes(N)
    _, K26, _ = blocks(N)
    basis = [(a, b) for a in combinations(range(N), 2) for b in combinations(range(N), 1)]
    b2 = [(b[0], (set(a) - {b[0]}).pop()) for (a, b) in basis if b[0] in a]
    (_, _, _), triples = z3_exact(N)
    extras = triples[1] + triples[2]
    if not extras:
        print(f"  N={N:>3} | (D),(E),(G),(H) OK | not resonant, no mirror pair to test")
        return

    wa = wb = wc = wy = 0.0
    seen, rows = set(), []
    for tau in extras:
        taup = tuple(sorted(n - k for k in tau))
        if tau in seen:
            continue
        seen |= {tau, taup}
        nrm, nrmp = slater_norm(U, tau, N), slater_norm(U, taup, N)

        for s in range(3):
            got = K26 @ lift(U, tau, s, N, nrm)
            pred = np.array([-((-1.0) ** b) * G(U, tau, b, c) / nrm
                             * (1.0 if (s == 1 or (s == 0 and c > b) or (s == 2 and c < b)) else 0.0)
                             for (b, c) in b2])
            wa = max(wa, float(np.max(np.abs(got - pred))))

        wb = max(wb, max(abs(G(U, taup, b, c) / nrmp - ((-1.0) ** c) * G(U, tau, b, c) / nrm)
                         for b in range(N) for c in range(N)))

        g = np.array([[G(U, tau, b, c) / nrm for c in range(N)] for b in range(N)])
        sg = np.array([(-1.0) ** c for c in range(N)])
        Tp = sum(sg[c] * g[b, c] ** 2 for b in range(N) for c in range(N) if c > b)
        Tm = sum(sg[c] * g[b, c] ** 2 for b in range(N) for c in range(N) if c < b)
        wc = max(wc, abs(Tp - Tm), abs(Tp), abs(Tm))

        W = np.array([lift(U, tau, s, N, nrm) for s in range(3)]).T
        Wp = np.array([lift(U, taup, s, N, nrmp) for s in range(3)]).T
        wy = max(wy, float(np.linalg.norm((K26 @ W).T @ (K26 @ Wp))))

        X = (K26 @ W).T @ (K26 @ W)
        ap = sum(g[b, c] ** 2 for b in range(N) for c in range(N) if c > b)
        am = sum(g[b, c] ** 2 for b in range(N) for c in range(N) if c < b)
        assert abs(ap - am) < TOL, f"(C) N={N} {tau}: a+ != a-"
        assert np.max(np.abs(X - ap * np.array([[1.0, 1, 0], [1, 2, 1], [0, 1, 1]]))) < TOL, \
            f"(I) N={N} {tau}: X is not a*[[1,1,0],[1,2,1],[0,1,1]]"
        assert np.allclose(np.sort(np.linalg.eigvalsh(X))[::-1], [3 * ap, ap, 0.0], atol=TOL), \
            f"(I) N={N} {tau}: spec(X) != (3a, a, 0)"
        assert abs(ap * n - (12.0 - sum(lam(k, n) ** 2 for k in tau))) < 1e-7, \
            f"(I) N={N} {tau}: a*n != 12 - sum lam^2"
        rows.append((tau, "ROT3" if tau in triples[1] else "PENT", ap, ap * n))

    assert wa < TOL and wb < TOL and wc < TOL and wy < TOL, "a matrix step failed"
    print(f"  N={N:>3} | (A) {wa:.1e} | (B) {wb:.1e} | (C) T+=T-=0 {wc:.1e} | "
          f"(D),(E),(G),(H) OK | ||Y|| {wy:.1e} | pairs {len(seen)//2}")
    for tau, kind, a, an in rows:
        print(f"         tau={tau} {kind}: a*n = {an:.9f}, levels (3a, a, 0) = "
              f"({3*a:.9f}, {a:.9f}, 0) = ({3*an:.0f}/n, {an:.0f}/n, 0)")


def cross_block_class(tau, sigma, n, triv):
    """Which of the FOUR populations a distinct-triple block falls into.

    'lemma4' : Lemma 4 applies (no shared mode) AND the reflection sign is +1.
    'byhand' : two distinct TRIV triples; eps = +1 always, the single shared mode is n/2,
               and the overlap of the two two-mode factors vanishes by sine orthogonality.
    'measured': reflection sign +1 with a shared mode, but not two TRIV: F = 0 is measured.
    'open'   : reflection sign -1; there F = 0 is an identity and carries no information.
    """
    if (sum(tau) + sum(sigma)) % 2 == 1:
        return "open"
    if not set(tau) & set(sigma):
        return "lemma4"
    if tau in triv and sigma in triv:
        assert set(tau) & set(sigma) == {n // 2}, f"TRIVxTRIV {tau},{sigma} share more than n/2"
        return "byhand"
    return "measured"


def check_cross_triple_blocks(N):
    """(J) The invariance the full-spectrum twinning needs, and how much of it is proved."""
    n, U = N + 1, modes(N)
    _, K26, _ = blocks(N)
    (_, _, _), tri = z3_exact(N)
    allt = tri[0] + tri[1] + tri[2]
    nrm = {t: slater_norm(U, t, N) for t in allt}
    KW = {t: K26 @ np.array([lift(U, t, s, N, nrm[t]) for s in range(3)]).T for t in allt}
    g = {t: np.array([[G(U, t, b, c) / nrm[t] for c in range(N)] for b in range(N)]) for t in allt}

    triv = set(tri[0])
    kinds = ("lemma4", "byhand", "measured", "open")
    cnt = {k: 0 for k in kinds}
    worst = dict({k: 0.0 for k in kinds}, mirror=0.0)
    for t, s in combinations(allt, 2):
        blk = float(np.linalg.norm(KW[t].T @ KW[s]))
        if tuple(sorted(n - k for k in t)) == s:
            worst["mirror"] = max(worst["mirror"], blk)
            continue
        kind = cross_block_class(t, s, n, triv)
        cnt[kind] += 1
        worst[kind] = max(worst[kind], blk)
        if kind in ("byhand", "measured"):  # the premise those two rest on
            assert abs(float(np.sum(g[t] * g[s]))) < 1e-8, f"(J) N={N}: F({t},{s}) != 0"
    assert worst["mirror"] < TOL, f"(F) N={N}: the mirror block Y is nonzero"
    for k in ("lemma4", "byhand", "measured"):
        assert worst[k] < TOL, f"(J) N={N}: a '{k}' block is nonzero"
    # the by-hand column is a closed form, not a count: every unordered pair of distinct TRIV triples
    assert cnt["byhand"] == len(triv) * (len(triv) - 1) // 2 == (N - 1) // 2 * ((N - 1) // 2 - 1) // 2,         f"(J) N={N}: byhand {cnt['byhand']} != C(#TRIV, 2)"
    tot = sum(cnt.values())
    print(f"  N={N:>3} | Y (mirror, PROVED) {worst['mirror']:.1e} | {tot} distinct non-mirror pairs: "
          f"{cnt['lemma4']} by Lemma 4, {cnt['byhand']} by hand (TRIVxTRIV), "
          f"{cnt['measured']} F=0 only measured, {cnt['open']} OPEN (measured {worst['open']:.1e})")
    return cnt["open"]


def check_generalization(N=11):
    """The theorem's real hypothesis is mode-disjointness, not vanishing. And it is not removable."""
    n, U = N + 1, modes(N)
    _, K26, _ = blocks(N)

    def ynorm(tau):
        taup = tuple(sorted(n - k for k in tau))
        W = np.array([lift(U, tau, s, N, slater_norm(U, tau, N)) for s in range(3)]).T
        Wp = np.array([lift(U, taup, s, N, slater_norm(U, taup, N)) for s in range(3)]).T
        return taup, float(np.linalg.norm((K26 @ W).T @ (K26 @ Wp)))

    seen, worst, pairs = set(), 0.0, 0
    for tau in combinations(range(1, N + 1), 3):
        taup = tuple(sorted(n - k for k in tau))
        if taup == tau or set(tau) & set(taup) or tau in seen:
            continue
        if abs(sum(lam(k, n) for k in tau)) < 1e-9:      # keep only the NON-vanishing ones
            continue
        seen |= {tau, taup}
        pairs += 1
        worst = max(worst, ynorm(tau)[1])
    assert worst < TOL, f"N={N}: a non-vanishing mode-disjoint triple has Y != 0"

    # sharpness: drop mode-disjointness and Y is large
    tau = (5, 7, 11)
    taup, bad = ynorm(tau)
    assert set(tau) & set(taup) and bad > 0.1, "the sharpness witness stopped being sharp"
    print(f"  N={N:>3} | {pairs} NON-vanishing mode-disjoint mirror pairs: worst ||Y|| {worst:.1e} "
          f"| witness {tau} shares {sorted(set(tau) & set(taup))} with {taup}: ||Y|| = {bad:.3f}")


def check_even_N_lab(N=14):
    """N even => n odd => every mirror pair has reflection sign -1, so the proof says nothing.
    Y = 0 holds anyway, and the pair-block level law spec(X) = (3a, a, 0), a = (12 - sum lam^2)/n
    (6/n for ROT3, 8/n for PENT), extends to even N unchanged. At N = 8 (the smallest even lab,
    n = 9, found 2026-07-12 via
    MirrorWorld's seed table r(inf) = 6 = 3*Z3) the vanishing set is EXACTLY one mirror pair
    {2,4,8}/{1,5,7}, so the Q-proved mirror specialisation carries the entire vanishing
    population there; distinct (non-mirror) pairs first exist at N = 14."""
    n, U = N + 1, modes(N)
    _, K26, _ = blocks(N)
    seen, worst, worst_x, pairs = set(), 0.0, 0.0, 0
    for tau in combinations(range(1, N + 1), 3):
        if abs(sum(lam(k, n) for k in tau)) > 1e-9:
            continue
        taup = tuple(sorted(n - k for k in tau))
        assert taup != tau, f"N={N}: a self-mirror vanishing triple exists, but n/2 is not an integer"
        if tau in seen:
            continue
        seen |= {tau, taup}
        pairs += 1
        assert (sum(tau) + sum(taup)) % 2 == 1, f"N={N}: expected reflection sign -1"
        W = np.array([lift(U, tau, s, N, slater_norm(U, tau, N)) for s in range(3)]).T
        Wp = np.array([lift(U, taup, s, N, slater_norm(U, taup, N)) for s in range(3)]).T
        worst = max(worst, float(np.linalg.norm((K26 @ W).T @ (K26 @ Wp))))
        for t, Wt in ((tau, W), (taup, Wp)):
            a = (12.0 - sum(lam(k, n) ** 2 for k in t)) / n
            ev = sorted(np.linalg.eigvalsh((K26 @ Wt).T @ (K26 @ Wt)))
            worst_x = max(worst_x, max(abs(e - w) for e, w in zip(ev, sorted([0.0, a, 3 * a]))))
    if N == 8:
        assert sorted(seen) == [(1, 5, 7), (2, 4, 8)], f"N=8 inventory changed: {sorted(seen)}"
    assert pairs and worst < TOL, f"N={N}: the even-N laboratory does not show Y = 0"
    assert worst_x < TOL, f"N={N}: the level law spec(X) = (3a, a, 0) fails at even N"
    print(f"  N={N:>3} (n={n} odd) | {pairs} mirror pairs, 0 self-mirror, reflection sign -1 "
          f"everywhere | ||Y|| <= {worst:.1e}, level law spec(X) = (3a, a, 0) to {worst_x:.1e}: "
          f"the proof misses it, the phenomenon does not")


def main():
    extra = [int(v) for v in sys.argv[1:]]
    print("Y = 0 and the coupled-level law, proved step by step.\n"
          "(A) K26 formula  (B) mirror law  (C) reflection  (D) sine orthogonality\n"
          "(E) mode-disjoint from mirror    (F) Y = 0\n"
          "(G) ||M_pq||^2 = n(sin^2 p + sin^2 q)   (H) ||D_tau||^2 = (n/2)^3   (I) spec(X) = (3a, a, 0)\n")

    print("The analytic steps (D),(E),(G),(H) and the level arithmetic, all odd N = 3..29:")
    for N in range(3, 30, 2):
        check_D_and_E(N)
        check_G_and_H(N)
        fams = check_level_law(N)
        if fams:
            print(f"  N={N:>3} | (D),(E),(G),(H) OK | {', '.join(fams)}")
    print("  (all other odd N in range: OK, not resonant)\n")

    print("The matrix steps (A),(B),(C),(F),(I) at the resonant N:")
    resonant = sorted({11, 17} | {v for v in extra if (v + 1) % 3 == 0 and v >= 11})
    for N in sorted({5, 7} | set(resonant)):
        check_matrix_steps(N)

    print("\nThe theorem's real hypothesis (mode-disjointness), and that it is not removable:")
    check_generalization(11)

    print("\nThe even-N laboratory for the open cases (reflection sign -1, Y = 0 regardless):")
    for even_n in (8, 14):
        check_even_N_lab(even_n)

    print("\n(J) The scope: what the full-spectrum twinning needs beyond Y = 0.")
    still_open = sum(check_cross_triple_blocks(N) for N in resonant)

    print("\nAll steps verified. Y = 0 is a theorem, uniform in N, for every triple mode-disjoint from\n"
          "its mirror; the pair block is X = a*[[1,1,0],[1,2,1],[0,1,1]], spec (3a, a, 0), with\n"
          "a*n = 12 - sum lam^2 = 6 (ROT3) or 8 (PENT).\n"
          f"STILL OPEN: {still_open} cross-triple blocks (over the N run here) are measured zero,\n"
          "not proved; the full-spectrum twinning rests on them. See (J).")


if __name__ == "__main__":
    main()
