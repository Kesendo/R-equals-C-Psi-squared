"""Topology controls the radical-writability of open-chain relaxation.

For an N-site graph under uniform-J XY coupling + uniform Z-dephasing, the (SE,DE) Liouvillian
sector factors over Q(i) into an AT-locked half (rates 2γ/6γ, the free-fermion frequencies, radical-
closed) and an H_B-mixed residue. The chain's H_B-mixed factors are the full symmetric group S_n
(no radical closure; see f89_pathk_galois). This script asks whether NON-chain topologies differ.

Result (gate-first, N=4,5,6): topology controls writability. The COMPLETE graph K_N is the writable
extreme (every H_B-mixed factor is a quartic-or-less, definitively radical-solvable, at N=4,5,6),
while the chain (S_8/S_18/S_32), the star (S_9 from N=5), and the ring (S_15 by N=6) all scramble to
the full symmetric group. Conjectured mechanism: K_N's massively degenerate single-particle spectrum
(adjacency {N−1 once, −1 with multiplicity N−1}) shatters the (SE,DE) block into small factors.

Method: build the FULL (SE,DE) block per topology (no symmetry projection; the diagonal rate −2/−6 =
overlap/no-overlap is topology-independent, only the hopping changes with the adjacency), Berkowitz/
sympy-factor over Q(i), bucket by AT rate (Re ∈ {−2,−6}), and read each H_B-mixed factor's Galois
group: deg ≤4 trivially solvable; deg5 non-solvable iff a 3-cycle appears (A5/S5 vs the solvable
C5/D5/F20); deg6 iff a 5-cycle appears (order divisible by 5); deg ≥7 via the generalised Jordan
certificate (a prime cycle in (n/2,n−3] ⟹ ⊇A_n ⟹ non-solvable). The chain gate reproduces the landed
S_8/S_18/S_32 exactly. Reuses the committed f89_pathk_galois engine.

See experiments/TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md; open arc f89_galois_open_doors (door A).
"""
from __future__ import annotations

import os
import sys
import time
from itertools import combinations

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from f89_path3_octic_galois import frobenius_cycle_types, is_irreducible_over_qi  # noqa: E402
from f89_pathk_galois import rate_bucket_factors, jordan_verdict                   # noqa: E402

LAM = sp.symbols("lambda")


def build_se_de_block(q0: int, N: int, adj: dict):
    """Full (SE,DE) Liouvillian block (NO symmetry projection) at integer q0, gamma=1, over
    Q(i), for the topology given by adjacency `adj` (adj[i] = list of neighbours of site i).
    Same hop signs as f89_pathk_galois.build_pathk_sym_over_qi; only neighbours generalised."""
    q = sp.Integer(q0)
    de_pairs = list(combinations(range(N), 2))
    basis = [(i, jk) for i in range(N) for jk in de_pairs]
    idx = {b: t for t, b in enumerate(basis)}
    nb = len(basis)
    L = sp.zeros(nb, nb)
    for col, (i, jk) in enumerate(basis):
        for i2 in adj[i]:                                  # SE hop 2q (ket): -i
            L[idx[(i2, jk)], col] += -sp.I * 2 * q
        j, k = jk                                          # DE hop 2q (bra): +i
        for nj in adj[j]:
            if nj != k:
                L[idx[(i, tuple(sorted((nj, k))))], col] += sp.I * 2 * q
        for nk in adj[k]:
            if nk != j:
                L[idx[(i, tuple(sorted((j, nk))))], col] += sp.I * 2 * q
        L[col, col] += (-2 if i in jk else -6)             # overlap 2g / no-overlap 6g
    return L


def adj_chain(N):    return {i: [x for x in (i - 1, i + 1) if 0 <= x < N] for i in range(N)}
def adj_ring(N):     return {i: sorted({(i - 1) % N, (i + 1) % N}) for i in range(N)}
def adj_star(N):     return {0: list(range(1, N)), **{i: [0] for i in range(1, N)}}
def adj_complete(N): return {i: [x for x in range(N) if x != i] for i in range(N)}


def classify_factor(fd, primes):
    """Solvability of Gal(fd/Q(i)) from Frobenius cycle types over split primes.
    deg<=4 always solvable. deg5: non-solvable iff a (3,1,1) 3-cycle appears (only A5/S5
    have 3-cycles; C5/D5/F20 do not). deg6: non-solvable iff a 5-cycle (..,5,..) appears
    (only A5/S5/A6/S6 on 6 points have order div by 5; solvable transitive deg-6 groups
    do not). deg>=7: Jordan window (a prime p in (n/2,n-3] cycle => >=A_n => non-solvable)."""
    n = sp.Poly(fd, LAM).degree()
    if n <= 4:
        return n, "solvable (deg<=4)", "solvable"
    if not is_irreducible_over_qi(fd, LAM):
        return n, "REDUCIBLE (unexpected)", "solvable"
    types = list(frobenius_cycle_types(fd, LAM, primes).keys())
    if n == 5:
        ns = any(3 in ct for ct in types)
        return n, ("A5/S5 (3-cycle seen)" if ns else "solvable C5/D5/F20 (no 3-cycle)"), \
            ("nonsolvable" if ns else "solvable")
    if n == 6:
        ns = any(5 in ct for ct in types)
        return n, ("non-solvable (5-cycle seen)" if ns else "solvable (no order-5 element)"), \
            ("nonsolvable" if ns else "solvable")
    verdict, _ = jordan_verdict(types, n)
    if verdict.startswith("S_") or "⊇A_" in verdict:
        return n, verdict, "nonsolvable"
    # Jordan found no big-prime cycle: NOT proven solvable, just undetermined.
    return n, f"UNDETERMINED ({verdict[:40]})", "undetermined"


def analyse(name, N, adj, prime_hi=3000):
    print("=" * 74)
    print(f"{name}   N={N}")
    print("=" * 74)
    t0 = time.time()
    L = build_se_de_block(2, N, adj)
    cp = sp.expand(L.charpoly(LAM).as_expr())
    at, hb = rate_bucket_factors(cp)
    at_deg = sorted((sp.Poly(f, LAM).degree() for f in at), reverse=True)
    hb_deg = sorted((sp.Poly(f, LAM).degree() for f in hb), reverse=True)
    print(f"  full (SE,DE) block {L.shape[0]}x{L.shape[0]}, charpoly deg {sp.degree(cp, LAM)} "
          f"(built+factored {time.time()-t0:.1f}s)")
    print(f"  AT-locked deg {at_deg}")
    print(f"  H_B-mixed deg {hb_deg}")
    primes = list(sp.primerange(5, prime_hi))
    states = []
    for f in sorted(hb, key=lambda x: -sp.Poly(x, LAM).degree()):
        n, verdict, status = classify_factor(f, primes)
        states.append(status)
        if n >= 5:
            print(f"    H_B-mixed deg {n}: {verdict}   [{status.upper()}]")
    if "nonsolvable" in states:
        head, ret = "hits a NON-SOLVABLE factor (unwritable, like the chain)", "nonsolvable"
    elif "undetermined" in states:
        head, ret = "all-solvable EXCEPT undetermined large factor(s): needs a direct group check", "undetermined"
    else:
        head, ret = "ALL H_B-mixed factors SOLVABLE  ==>  fully writable", "solvable"
    print(f"  ==> {name} N={N}: {head}\n")
    return ret


def predicted_complete_hb_histogram(N):
    """The DERIVED H_B-mixed factor-degree histogram for K_N (N>=5): the S_N-irrep multiplicities in
    V = SE(x)DE = M^(N-2,1,1) (+) M^(N-3,2,1) cap at 4 (Kostka), so the H_B-mixed factors are exactly
    (N-1) quartics [standard rep [N-1,1], mult 4], N(N-3)/2 cubics [[N-2,2], mult 3], and
    (N-1)(N-2)/2 quadratics [[N-2,1,1], mult 2]. All degree <= 4 => radical-solvable, for all N."""
    return {4: N - 1, 3: N * (N - 3) // 2, 2: (N - 1) * (N - 2) // 2}


def verify_complete(N):
    """Gate the derivation: the actual K_N H_B-mixed factor-degree histogram must equal the
    rep-theory prediction (and so max degree = 4)."""
    L = build_se_de_block(2, N, adj_complete(N))
    cp = sp.expand(L.charpoly(LAM).as_expr())
    _, hb = rate_bucket_factors(cp)
    hist = {}
    for f in hb:
        d = sp.Poly(f, LAM).degree()
        hist[d] = hist.get(d, 0) + 1
    pred = predicted_complete_hb_histogram(N)
    ok = hist == pred
    print(f"  K_{N}: H_B-mixed degrees {dict(sorted(hist.items(), reverse=True))}  "
          f"vs rep-theory {dict(sorted(pred.items(), reverse=True))}   [{'MATCH' if ok else 'MISMATCH'}]")
    return ok


def ring_niven_split(N, ext, ext_label):
    """The ring's growth law is a GALOIS inflation by φ(N)/2 = [Q(cos(2π/N)):Q]. The dihedral D_N has
    an irrational character table (cos(2π/N)), so over Q(i) its Galois-conjugate 2-dim irreps merge and
    inflate the H_B-mixed factor degree by φ(N)/2 (> 1 except at the Niven-rational sizes N ∈ {1,2,3,4,6}
    where cos(2π/N) is rational). Gate: the max Q(i)-factor splits over Q(i, cos(2π/N)) into φ(N)/2 pieces."""
    from sympy import totient
    L = build_se_de_block(2, N, adj_ring(N))
    cp = sp.expand(L.charpoly(LAM).as_expr())
    _, hb = rate_bucket_factors(cp)
    f = max(hb, key=lambda p: sp.Poly(p, LAM).degree())
    dQi = sp.Poly(f, LAM).degree()
    fac = sp.factor(f, LAM, extension=ext)
    args = fac.args if fac.func == sp.Mul else [fac]
    degs = sorted((sp.Poly(a.base if a.func == sp.Pow else a, LAM).degree()
                   for a in args if (a.base if a.func == sp.Pow else a).has(LAM)), reverse=True)
    galois = totient(N) // 2 if N >= 3 else 1
    ok = degs == [dQi // galois] * galois
    print(f"  ring N={N}: max H_B factor deg over Q(i) = {dQi}; over {ext_label} -> {degs}; "
          f"φ(N)/2 = {galois}   [{'MATCH' if ok else 'see'}]")
    return ok


def verify_star(N):
    """Gate the star derivation: the star Liouvillian commutes with S_{N-1} (leaf permutations, hub
    fixed); the multiplicity of the standard irrep std_{N-1} in V is 9 (N-independent), so the max
    H_B-mixed factor degree is 9 and the number of degree-9 factors is dim(std_{N-1}) = N-2. The star
    scrambles to a FIXED S_9 for all N >= 5 (bounded, but > 4, so not writable)."""
    L = build_se_de_block(2, N, adj_star(N))
    cp = sp.expand(L.charpoly(LAM).as_expr())
    _, hb = rate_bucket_factors(cp)
    degs = sorted((sp.Poly(f, LAM).degree() for f in hb), reverse=True)
    maxd = degs[0] if degs else 0
    n9 = sum(1 for d in degs if d == 9)
    ok = maxd == 9 and n9 == N - 2
    print(f"  star N={N}: max H_B degree {maxd}, {n9} factor(s) of degree 9   "
          f"vs predicted max 9, {N-2} of degree 9   [{'MATCH' if ok else 'MISMATCH'}]")
    return ok


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "niven":
        print("Ring growth law = Galois inflation by φ(N)/2 (Niven: cos(2π/N) rational iff N in {1,2,3,4,6}).\n")
        ok6 = ring_niven_split(6, [sp.I], "Q(i) [cos(2π/6)=1/2 rational, no inflation]")
        ok5 = ring_niven_split(5, [sp.I, sp.sqrt(5)], "Q(i,√5) [cos(2π/5), φ/2=2]")
        print(f"\nNIVEN {'PASS:ring deg-16 (N=5) splits 8+8 over Q(i,√5); deg-15 (N=6) stays (rational)' if (ok5 and ok6) else 'FAIL'}")
        sys.exit(0 if (ok5 and ok6) else 1)
    if len(sys.argv) >= 2 and sys.argv[1] == "verify-star":
        print("Gate: star max H_B degree == 9 (mult of std_{N-1}), with N-2 degree-9 factors\n")
        targets = [int(a) for a in sys.argv[2:]] or [5, 6, 7, 8, 9]
        allok = all(verify_star(N) for N in targets)
        print(f"\nVERIFY-STAR {'PASS:star caps at S_9 (bounded, N-independent), all N' if allok else 'FAIL'}")
        sys.exit(0 if allok else 1)
    if len(sys.argv) >= 2 and sys.argv[1] == "verify":
        print("Gate: K_N H_B-mixed factor degrees == rep-theory prediction {4:N-1, 3:N(N-3)/2, 2:(N-1)(N-2)/2}\n")
        targets = [int(a) for a in sys.argv[2:]] or [5, 6, 7, 8]
        allok = all(verify_complete(N) for N in targets)
        print(f"\nVERIFY {'PASS:derivation confirmed (max degree 4, all N)' if allok else 'FAIL'}")
        sys.exit(0 if allok else 1)
    Ns = [int(a) for a in sys.argv[1:]] or [4, 5, 6]
    PRIME_HI = 4000    # enough for every clean S_n detection; raise (e.g. 20000) to probe a stubborn
    #                    'undetermined' factor (the ring N=5 deg-16 stays undetermined even at 20000)
    print("GATE: chain N=4 (= path-3) must show a degree-8 non-solvable (S_8) H_B-mixed factor.\n")
    gate = analyse("GATE chain", 4, adj_chain(4), prime_hi=PRIME_HI)
    print(f"GATE {'PASS' if gate == 'nonsolvable' else 'FAIL: builder wrong, stop'}\n")
    if gate != "nonsolvable":
        sys.exit(1)
    topo = {"complete": adj_complete, "star": adj_star, "ring": adj_ring, "chain": adj_chain}
    summary = {}
    for N in Ns:
        for name, fn in topo.items():
            summary[(name, N)] = analyse(f"{name.upper()}", N, fn(N), prime_hi=PRIME_HI)
    print("=" * 74)
    print("SUMMARY  (True = hits a non-solvable H_B-mixed factor, like the chain)")
    print("=" * 74)
    label = {"nonsolvable": "UNWRITABLE (non-solvable, like chain)",
             "undetermined": "undetermined large factor (needs direct group check)",
             "solvable": "WRITABLE (all factors solvable)"}
    for (name, N), v in summary.items():
        print(f"  {name:9s} N={N}: {label[v]}")
