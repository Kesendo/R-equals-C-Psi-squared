"""The break-input the no-rounding convention demands for "[L, S+] is exactly 0.0".

An exactness claim needs an input of the SAME KIND that makes it fail, or it is a claim about
the inputs tried rather than about the object. F142's Lemma 2.3 says the spin ladder commutes
iff Sigma h Sigma = -h, i.e. iff h connects only sites of OPPOSITE parity. So the break-inputs
are geometries that join sites of EQUAL parity, and each should break S+ while leaving Phi
untouched, because the eta ladder needs no condition on h at all.

  open chain, nearest neighbour   opposite parity only        S+ should hold
  + next-nearest-neighbour bond   joins l and l+2, same       S+ should break
  + on-site potential            h_ll, trivially same         S+ should break
  ring at ODD N                  bond (N-1, 0) joins same     S+ should break
  ring at EVEN N                 all bonds opposite           S+ should hold
  star                           centre-vs-leaves, not parity S+ should break (F142 Cor 2.4)

Phi is the control on every row: F125 Theorem B and F142 Lemma 2.2 both say it holds for ANY
quadratic particle-conserving H, so a Phi failure anywhere means the build is wrong.

Companion of eta_ladder_blocks.py. Recorded in the sideways_spin_ladder open arc.
"""
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from eta_ladder_blocks import block_basis, build_ladder, bit, popcount  # noqa: E402


def build_L_edges(N, edges, onsite, gamma, p, q):
    """L on block (p,q) for H = sum_(i,j,t) t (sigma+_i sigma-_j + h.c.) + sum_l v_l n_l."""
    pairs, index = block_basis(N, p, q)
    n = len(pairs)
    L = np.zeros((n, n), dtype=complex)

    def hops(a):
        """Fermionic c_i^dag c_j on a bit string, JW STRING INCLUDED.

        The first version of this used sigma+_i sigma-_j, which equals c_i^dag c_j only for
        ADJACENT i, j: for a longer bond the string between them is missing and the operator is
        not the quadratic fermionic H that F125 Theorem B and F142 Lemma 2.2 are about.  The Phi
        control caught it, failing on exactly the non-adjacent geometries; that is what the
        control is for.
        """
        out = []
        for (i, j, t) in edges:
            for (src, dst) in ((j, i), (i, j)):            # c^dag_dst c_src, both directions
                bs, bd = bit(src, N), bit(dst, N)
                if not (a & bs):
                    continue
                s1 = -1.0 if popcount(a >> (N - src)) % 2 else 1.0
                a1 = a & ~bs
                if a1 & bd:
                    continue
                s2 = -1.0 if popcount(a1 >> (N - dst)) % 2 else 1.0
                out.append((a1 | bd, t * s1 * s2))
        d = sum(v for l, v in enumerate(onsite) if a & bit(l, N))
        if d:
            out.append((a, d))
        return out

    for col, (a, b) in enumerate(pairs):
        diff = sum(gamma[l] for l in range(N) if bool(a & bit(l, N)) != bool(b & bit(l, N)))
        L[col, col] += -2.0 * diff
        for a2, amp in hops(a):
            L[index[(a2, b)], col] += -1j * amp
        for b2, amp in hops(b):
            L[index[(a, b2)], col] += +1j * amp
    return L


def worst_residual(N, edges, onsite, gamma, kind):
    dq = 1 if kind == "phi" else -1
    worst = 0.0
    for p in range(N + 1):
        for q in range(N + 1):
            if not (0 <= p + 1 <= N and 0 <= q + dq <= N):
                continue
            M = build_ladder(N, p, q, kind)
            if M.size == 0:
                continue
            Ls = build_L_edges(N, edges, onsite, gamma, p, q)
            Lt = build_L_edges(N, edges, onsite, gamma, p + 1, q + dq)
            worst = max(worst, np.linalg.norm(Lt @ M - M @ Ls))
    return worst


def chain(N, rng):
    return [(i, i + 1, float(0.4 + 1.2 * rng.random())) for i in range(N - 1)]


if __name__ == "__main__":
    rng = np.random.default_rng(4242)
    print(f"  {'geometry':<34} {'N':>2}  {'[L,Phi]':>11}  {'[L,S+]':>11}   verdict")
    print("  " + "-" * 78)
    cases = []
    for N in (5, 6):
        g = list(0.05 + 0.6 * rng.random(N))
        z = [0.0] * N
        cases.append((f"open chain, nearest neighbour", N, chain(N, rng), z, g, True))
        cases.append((f"chain + next-nearest bond (0,2)", N, chain(N, rng) + [(0, 2, 0.7)], z, g, False))
        cases.append((f"chain + on-site potential", N, chain(N, rng),
                      [0.0] * (N - 1) + [0.9], g, False))
        ring = chain(N, rng) + [(N - 1, 0, 0.83)]
        cases.append((f"ring (closes {N-1}-0)", N, ring, z, g, N % 2 == 0))
        star = [(0, j, float(0.4 + 1.2 * rng.random())) for j in range(1, N)]
        cases.append((f"star, centre 0", N, star, z, g, False))

    bad = 0
    phi_all_zero = True
    for (label, N, edges, onsite, gamma, expect_sp_holds) in cases:
        rp = worst_residual(N, edges, onsite, gamma, "phi")
        phi_all_zero = phi_all_zero and rp == 0.0
        rs = worst_residual(N, edges, onsite, gamma, "splus")
        sp_holds = rs == 0.0
        ok = (rp == 0.0) and (sp_holds == expect_sp_holds)
        bad += not ok
        verdict = ("S+ holds" if sp_holds else "S+ BREAKS") + ("" if ok else "   <-- UNEXPECTED")
        print(f"  {label:<34} {N:>2}  {rp:>11.3e}  {rs:>11.3e}   {verdict}")

    print()
    print(f"  Phi held on every row (it needs no condition on h): {phi_all_zero}")
    print(f"  rows behaving as F142 Lemma 2.3 predicts: {len(cases) - bad}/{len(cases)}")
