"""
The two ladders, block-locally, so N = 6 and 7 are cheap.

The intertwining is a statement PER RUNG, so it does not need the 4^N superoperator: build
L on a joint-popcount block (p,q), build the ladder as the rectangular map into its target
block, and check

    L_target . Ladder  -  Ladder . L_source  =  0

directly.  Blocks are at most C(7,3)^2 = 1225 at N = 7, so this reaches where the dense
construction cannot.

  Phi(|a><b|) = sum_l s_l(a) s_l(b) |a+e_l><b+e_l|          (p,q) -> (p+1, q+1)
  S+ (|a><b|) = sum_l (-1)^l s_l(a) t_l(b) |a+e_l><b-e_l|   (p,q) -> (p+1, q-1)

with the Jordan-Wigner string signs s_l(x) = (-1)^(occupied sites before l in x) for the
creation on either side, and t_l for the annihilation.

CONTROL: Phi's own commutator, on every rung.  Both F125 Theorem B and F142 Lemma 2.2 say it
is zero for any quadratic particle-conserving H, so a nonzero Phi residual means the build is
wrong and the S+ column may not be read.  (This file does not itself perform the dense-vs-block
cross-check; that was a one-off scout and is not committed.)

Companion of eta_ladder_chain.py (which gates numbers) and eta_ladder_breakinput.py (which
shows where the spin ladder stops holding).  Recorded in the sideways_spin_ladder open arc.
"""
import sys

import numpy as np


def popcount(x):
    return bin(x).count("1")


def bit(l, N):
    return 1 << (N - 1 - l)


def sign_before(x, l, N):
    """(-1)^(number of occupied sites m < l in x)."""
    return -1.0 if popcount(x >> (N - l)) % 2 else 1.0


def block_basis(N, p, q):
    A = [a for a in range(1 << N) if popcount(a) == p]
    B = [b for b in range(1 << N) if popcount(b) == q]
    pairs = [(a, b) for a in A for b in B]
    return pairs, {ab: i for i, ab in enumerate(pairs)}


def hamiltonian_hops(a, N, J):
    """H = sum_b J_b (sigma+_i sigma-_j + h.c.) on bonds (i, i+1): the XY chain's hopping.

    NEAREST NEIGHBOUR ONLY, and that restriction is load-bearing, not incidental: sigma+_i
    sigma-_j equals the fermionic c_i^dag c_j only when i and j are ADJACENT, because otherwise
    the Jordan-Wigner string between them is missing and the H is no longer the quadratic
    fermionic operator F125 Theorem B and F142 Lemma 2.2 are about.  Adding a longer bond here
    silently leaves that class.  eta_ladder_breakinput.py builds the general case correctly.
    Returns [(a', amplitude)] with H|a> = sum amplitude |a'>."""
    out = []
    for i in range(N - 1):
        j = i + 1
        bi, bj = bit(i, N), bit(j, N)
        if (a & bi) and not (a & bj):
            out.append((a ^ bi ^ bj, J[i]))
        elif (a & bj) and not (a & bi):
            out.append((a ^ bi ^ bj, J[i]))
    return out


def build_L(N, p, q, J, gamma):
    pairs, index = block_basis(N, p, q)
    n = len(pairs)
    L = np.zeros((n, n), dtype=complex)
    for col, (a, b) in enumerate(pairs):
        diff = sum(gamma[l] for l in range(N) if bool(a & bit(l, N)) != bool(b & bit(l, N)))
        L[col, col] += -2.0 * diff
        for a2, amp in hamiltonian_hops(a, N, J):          # -i H rho
            L[index[(a2, b)], col] += -1j * amp
        for b2, amp in hamiltonian_hops(b, N, J):          # +i rho H   (H real symmetric)
            L[index[(a, b2)], col] += +1j * amp
    return L, pairs, index


def build_ladder(N, p, q, kind):
    """kind 'phi': (p,q)->(p+1,q+1);  kind 'splus': (p,q)->(p+1,q-1)."""
    dq = 1 if kind == "phi" else -1
    src, _ = block_basis(N, p, q)
    tgt, tindex = block_basis(N, p + 1, q + dq)
    M = np.zeros((len(tgt), len(src)), dtype=complex)
    for col, (a, b) in enumerate(src):
        for l in range(N):
            bl = bit(l, N)
            if a & bl:
                continue                                    # creation on the left needs l empty
            if kind == "phi":
                if b & bl:
                    continue                                # creation on the right too
                nb = b | bl
            else:
                if not (b & bl):
                    continue                                # annihilation on the right needs l full
                nb = b & ~bl
            stagger = ((-1.0) ** l) if kind == "splus" else 1.0
            M[tindex[(a | bl, nb)], col] += stagger * sign_before(a, l, N) * sign_before(b, l, N)
    return M


def rung_report(N, J, gamma, kind, verbose=True):
    """Every rung of the given ladder: intertwining residual, rank, and spectral containment."""
    dq = 1 if kind == "phi" else -1
    rows = []
    for p in range(N + 1):
        for q in range(N + 1):
            tp, tq = p + 1, q + dq
            if not (0 <= tp <= N and 0 <= tq <= N):
                continue
            Ls, src, _ = build_L(N, p, q, J, gamma)
            Lt, tgt, _ = build_L(N, tp, tq, J, gamma)
            M = build_ladder(N, p, q, kind)
            if M.size == 0:
                continue
            resid = np.linalg.norm(Lt @ M - M @ Ls)
            sv = np.linalg.svd(M, compute_uv=False) if min(M.shape) else np.array([0.0])
            rank = int((sv > 1e-10 * max(sv[0], 1e-300)).sum())
            es, et = np.linalg.eigvals(Ls), np.linalg.eigvals(Lt)
            pool, shared = list(es), 0
            for lam in et:                      # is the TARGET spectrum inside the SOURCE one
                if not pool:
                    break
                k = int(np.argmin([abs(lam - m) for m in pool]))
                if abs(lam - pool[k]) < 1e-8 * max(1.0, abs(lam)):
                    shared += 1
                    pool.pop(k)
            rows.append((p, q, tp, tq, len(src), len(tgt), resid, rank, shared))
    if verbose:
        print(f"  {'rung':>17} {'dim':>11}  {'||L M - M L||':>13} {'rank':>5} "
              f"{'target spec in source':>22}")
        for (p, q, tp, tq, ns, nt, r, rk, sh) in rows:
            # The smaller block embeds in the larger, so the count to expect is min(ns, nt),
            # not nt.  Comparing against nt alone flags every INJECTIVE rung as a failure and
            # printed ten false alarms at N=5 while the arc quoting this tool said zero.
            flag = "" if sh == min(ns, nt) else "   <-- SMALLER SPECTRUM NOT CONTAINED"
            print(f"  ({p},{q})->({tp},{tq}) {ns:>5} ->{nt:>5}  {r:>13.3e} {rk:>5} "
                  f"{sh:>10}/{nt:<10}{flag}")
    return rows


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    rng = np.random.default_rng(20260807 + N)
    J = 0.4 + 1.2 * rng.random(N - 1)
    gamma = 0.05 + 0.6 * rng.random(N)
    print(f"N = {N}, non-uniform profile")
    print(f"  J     = {np.round(J, 6).tolist()}")
    print(f"  gamma = {np.round(gamma, 6).tolist()}")
    print()
    print("Phi rungs (p,q) -> (p+1,q+1):")
    rp = rung_report(N, J, gamma, "phi")
    print()
    print("S+ rungs (p,q) -> (p+1,q-1):")
    rs = rung_report(N, J, gamma, "splus")
    print()
    worst_phi = max((r[6] for r in rp), default=0.0)
    worst_sp = max((r[6] for r in rs), default=0.0)
    print(f"  worst Phi intertwining residual over all rungs: {worst_phi:.3e}")
    print(f"  worst S+  intertwining residual over all rungs: {worst_sp:.3e}")
    bad = [r for r in rs if r[8] != min(r[4], r[5])]
    print(f"  S+ rungs where the smaller spectrum does NOT embed in the larger: {len(bad)}")
    for r in bad:
        print(f"      ({r[0]},{r[1]})->({r[2]},{r[3]}) dims {r[4]}->{r[5]} rank {r[7]} shared {r[8]}")
