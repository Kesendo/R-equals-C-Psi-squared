"""The two AT eigenspaces of block (1,3) of the XY chain at N = 6 (lambda = -8 and -4, gamma = 1), read off the
heights on the letter cube: the second route beside F142's frozen band (experiments/F89_PATH_K_DIABOLIC.md, the
walk from a (1,3) seed, where -8 is l=2 under S+ and -4 is l=2 under Phi).

Both ladders keep every cell's height k = popcount(a XOR b) (PROOF_FROZEN_BAND_SO4 Lemma 2.1), and every cell of a
boundary block (0, q) has height q, so L there is -2q + iE with E an eigenvalue of H on q particles. The S+ chain
ends at (0,4) and (4,0), the Phi chain at (0,2) and (4,6). Each l=2 family in (1,3) is its ladder's image of the
boundary block (0,4) resp. (0,2). This gate checks, exactly, and fails unless every item holds:

  - the four chain ends sit at one height each, 4 for S+ and 2 for Phi;
  - the whole image of (0,4) under S+ lies on the height-4 cells of (1,3), the whole image of (0,2) under Phi on
    the height-2 cells;
  - the zero-energy states of H on 4 and on 2 particles (the rational nullspace of H's integer block, sympy) are
    three each;
  - their ladder images are eigenvectors of L on (1,3) at -8 resp. -4, checked in integers (L has Gaussian-integer
    entries: -2k on the diagonal, +-iJ for a hop), of rank three, on the height of their eigenvalue;
  - the nullity of L - lambda on (1,3) is three modulo two primes p = 1 mod 4 (i -> a square root of -1), an upper
    bound on the rational nullity, so the images span the whole eigenspace;
  - the bra flip P_f: rho -> rho X^N (F89_PATH_K_DIABOLIC's name) keeps (1,3) and carries one side onto the
    other: P_f Phi P_f = S+ on the rung
    (0,4) -> (1,3), and P_f carries the three images at -4 onto eigenvectors at -8, of rank three (P_f commutes with
    the Hamiltonian part and swaps the heights 2 and 4, so on the kernel of the Hamiltonian part it sends -4 to -8);
  - at the uniform chain (its kernels do not depend on the scale of J) and at the bond profile [1, 2, 3, 2, 5];
  - CONTROL: a potential n_0 + n_5 on the two end sites (the vacuum keeps energy 0, so a real eigenvalue on (0, q)
    still needs E = 0) leaves two zero-sum pairs (eps, -eps) among the single-particle energies, so two
    zero-energy states of two particles. On the -4 side,
    where Phi keeps intertwining for any particle-conserving quadratic H (Lemma 2.2), the eigenspace must follow the
    count down to two: two zero-energy states, two exact image eigenvectors, nullity two. The potential also breaks
    S+'s intertwining (Lemma 2.3 needs Sigma h Sigma = -h), so the -8 side is printed, not gated.

Run: python simulations/at_eigenspace_heights.py  (a few seconds; exit 0 iff every gate passes)."""
import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eta_ladder_blocks import block_basis, bit, sign_before, popcount

N = 6
PRIMES = (1000000009, 998244353)          # both 1 mod 4


def hops(a, J):
    out = []
    for i in range(N - 1):
        bi, bj = bit(i, N), bit(i + 1, N)
        if bool(a & bi) != bool(a & bj):
            out.append((a ^ bi ^ bj, J[i]))
    return out


def potential(x, mu):
    return sum(mu[l] for l in range(N) if x & bit(l, N))


def L_entries(p, q, J, mu):
    """L on block (p,q) as {(row, col): (re, im)} with integer parts, gamma = 1 on every site."""
    pairs, index = block_basis(N, p, q)
    E = {}

    def add(r, c, re, im):
        x = E.get((r, c), (0, 0))
        E[(r, c)] = (x[0] + re, x[1] + im)
    for col, (a, b) in enumerate(pairs):
        add(col, col, -2 * popcount(a ^ b), -(potential(a, mu) - potential(b, mu)))
        for a2, amp in hops(a, J):
            add(index[(a2, b)], col, 0, -amp)       # -i H rho
        for b2, amp in hops(b, J):
            add(index[(a, b2)], col, 0, +amp)       # +i rho H
    return pairs, E


def ladder(p, q, kind):
    """The map of eta_ladder_blocks.build_ladder, with integer entries."""
    dq = 1 if kind == "phi" else -1
    src, _ = block_basis(N, p, q)
    tgt, tindex = block_basis(N, p + 1, q + dq)
    M = {}
    for col, (a, b) in enumerate(src):
        for l in range(N):
            bl = bit(l, N)
            if a & bl:
                continue
            if kind == "phi":
                if b & bl:
                    continue
                nb = b | bl
            else:
                if not (b & bl):
                    continue
                nb = b & ~bl
            s = ((-1) ** l if kind == "splus" else 1) * int(sign_before(a, l, N)) * int(sign_before(b, l, N))
            key = (tindex[(a | bl, nb)], col)
            M[key] = M.get(key, 0) + s
    return M, tgt


def zero_energy_states(q, J, mu):
    """On (0,q) the ket is the vacuum (energy 0), L(|0><b|) = -2q |0><b| + i |0><b| H, so a real eigenvalue needs
    sum_b v_b <b|H = 0: the rational nullspace of H's (real symmetric, integer) q-particle block."""
    pairs, index = block_basis(N, 0, q)
    H = sp.zeros(len(pairs), len(pairs))
    for col, (a, b) in enumerate(pairs):
        H[col, col] += potential(b, mu)
        for b2, amp in hops(b, J):
            H[index[(a, b2)], col] += amp
    out = []
    for v in H.nullspace():
        den = sp.ilcm(*[sp.fraction(x)[1] for x in v])
        out.append([int(x * den) for x in v])
    return out


def apply(M, v, nt):
    w = [0] * nt
    for (r, c), s in M.items():
        w[r] += s * v[c]
    return w


def is_eigen(E, v, lam):
    re = [0] * len(v)
    im = [0] * len(v)
    for (r, c), (a, b) in E.items():
        re[r] += a * v[c]
        im[r] += b * v[c]
    return all(re[r] == lam * v[r] and im[r] == 0 for r in range(len(v)))


def nullity_mod_p(E, n, lam, p):
    r = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
    i_mod = pow(r, (p - 1) // 4, p)
    assert (i_mod * i_mod) % p == p - 1
    A = [[0] * n for _ in range(n)]
    for (rr, c), (a, b) in E.items():
        A[rr][c] = (a + b * i_mod) % p
    for d in range(n):
        A[d][d] = (A[d][d] - lam) % p
    rank, row = 0, 0
    for col in range(n):
        piv = next((x for x in range(row, n) if A[x][col]), None)
        if piv is None:
            continue
        A[row], A[piv] = A[piv], A[row]
        inv = pow(A[row][col], p - 2, p)
        A[row] = [(x * inv) % p for x in A[row]]
        for x in range(n):
            if x != row and A[x][col]:
                f = A[x][col]
                A[x] = [(A[x][y] - f * A[row][y]) % p for y in range(n)]
        row += 1
        rank += 1
    return n - rank


def bra_flip(v, cells):
    """P_f: |a><b| -> |a><b XOR 1...1|, on a vector over the cells of a block (p, q) with 2q = N, here (1,3)."""
    index = {c: i for i, c in enumerate(cells)}
    full = (1 << N) - 1
    w = [0] * len(v)
    for i, (a, b) in enumerate(cells):
        if v[i]:
            w[index[(a, b ^ full)]] += v[i]
    return w


def r_phi_r_is_splus():
    """P_f Phi P_f = S+ on the rung (0,4) -> (1,3), entry by entry in integers, the sites counted from 0."""
    src4, _ = block_basis(N, 0, 4)
    src2, idx2 = block_basis(N, 0, 2)
    Mphi, tgt = ladder(0, 2, "phi")
    Msp, tgt2 = ladder(0, 4, "splus")
    assert tgt == tgt2
    full = (1 << N) - 1
    for col, (a, b) in enumerate(src4):
        e = [0] * len(src2)
        e[idx2[(a, b ^ full)]] = 1
        lhs = bra_flip(apply(Mphi, e, len(tgt)), tgt)
        rhs = [Msp.get((r, col), 0) for r in range(len(tgt))]
        if lhs != rhs:
            return False
    return True


def heights_of(block):
    return sorted({popcount(a ^ b) for (a, b) in block_basis(N, *block)[0]})


SIDES = (  # (eigenvalue, boundary block, its ladder, the chain's two ends)
    (-8, (0, 4), "splus", ((0, 4), (4, 0))),
    (-4, (0, 2), "phi", ((0, 2), (4, 6))),
)


def read(J, mu):
    """One reading per side: (zero-energy count, exact eigenvectors, rank, image heights, family heights, ends'
    heights, nullities mod p)."""
    cells, E13 = L_entries(1, 3, J, mu)
    out = []
    flipped = None
    for lam, boundary, kind, ends in SIDES:
        zs = zero_energy_states(boundary[1], J, mu)
        M, tgt = ladder(*boundary, kind)
        assert tgt == cells
        images = [apply(M, v, len(cells)) for v in zs]
        out.append(dict(
            lam=lam, kind=kind, zero=len(zs),
            exact=all(is_eigen(E13, w, lam) for w in images),
            rank=sp.Matrix(images).rank() if images else 0,
            image_heights=sorted({popcount(cells[r][0] ^ cells[r][1]) for w in images for r in range(len(w)) if w[r]}),
            family_heights=sorted({popcount(cells[r][0] ^ cells[r][1]) for (r, c) in M}),
            ends={e: heights_of(e) for e in ends},
            nullity=[nullity_mod_p(E13, len(cells), lam, p) for p in PRIMES]))
        if kind == "phi":
            flipped = [bra_flip(w, cells) for w in images]
    # P_f carries the -4 images onto eigenvectors at -8
    out.append(dict(r_exact=all(is_eigen(E13, w, -8) for w in flipped) if flipped else False,
                    r_rank=sp.Matrix(flipped).rank() if flipped else 0))
    return out


def gate(J, label):
    print(f"--- {label}: J = {J}")
    ok = True
    *sides, rmap = read(J, [0] * N)
    for s in sides:
        h = -s["lam"] // 2
        good = (all(v == [h] for v in s["ends"].values()) and s["family_heights"] == [h]
                and s["zero"] == 3 and s["exact"] and s["rank"] == 3 and s["image_heights"] == [h]
                and s["nullity"] == [3, 3])
        print(f"  lambda = {s['lam']} ({s['kind']}): ends {s['ends']}, l=2 family on heights {s['family_heights']}; "
              f"zero-energy states {s['zero']}, images exact eigenvectors {s['exact']}, rank {s['rank']}, "
              f"on heights {s['image_heights']}; nullity of L - lambda on (1,3) mod p {s['nullity']}  "
              f"{'ok' if good else 'FAILED'}")
        ok &= good
    good = rmap["r_exact"] and rmap["r_rank"] == 3
    print(f"  the bra flip P_f: the -4 images carried to eigenvectors at -8 {rmap['r_exact']}, rank {rmap['r_rank']}  "
          f"{'ok' if good else 'FAILED'}")
    return ok and good


def control():
    print("--- CONTROL: the uniform chain with a potential n_0 + n_5 on the end sites, which leaves two zero-sum pairs (eps, -eps)")
    ok = True
    for s in read([1] * (N - 1), [1] + [0] * (N - 2) + [1])[:2]:
        line = (f"  lambda = {s['lam']} ({s['kind']}): zero-energy states {s['zero']}, images exact eigenvectors "
                f"{s['exact']}, rank {s['rank']}, nullity mod p {s['nullity']}")
        if s["kind"] == "phi":
            good = s["zero"] == 2 and s["exact"] and s["rank"] == 2 and s["nullity"] == [2, 2]
            print(line + f"  (Phi keeps intertwining: the eigenspace follows the count to two) {'ok' if good else 'FAILED'}")
            ok &= good
        else:
            print(line + "  (S+ no longer intertwines: printed only)")
    return ok


if __name__ == "__main__":
    good = gate([1] * (N - 1), "the uniform chain")
    good &= gate([1, 2, 3, 2, 5], "the bond profile [1, 2, 3, 2, 5]")
    good &= r_phi_r_is_splus()
    print(f"P_f Phi P_f = S+ on the rung (0,4) -> (1,3), exactly: {r_phi_r_is_splus()}")
    good &= control()
    print("at_eigenspace_heights: ALL GATES PASS" if good else "at_eigenspace_heights: FAILED")
    sys.exit(0 if good else 1)
