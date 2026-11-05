"""Gate for docs/THE_GREY_GLASS.md: rooms of operator space where the watching is one number.

Exact throughout: integer couplings and rates, generators built as full superoperators on the 4^N
cells from Kronecker products of integer matrices (complex entries with small integer or half-integer
parts, exact in floating point; every comparison is ==, no tolerance, no eigensolver). Self-contained;
imports nothing from the framework.

Run:
    python simulations/grey_glass_gate.py

What it checks:
  G1 the boundary blocks: under uniform Z dephasing a joint-popcount block (p, q) carries one Hamming
     distance exactly when min(p, q) = 0 or max(p, q) = N (N = 2..7, every block); on the full
     generator of an XXZ chain and ring with integer couplings and Z fields (N = 3, 4), every such
     block is kept by L and L - c Id equals L_H on it entry for entry, c = -2 gamma |p - q|, while on
     (1, 1) L - L_H is a diagonal of several values, so no c works (the control).
  G2 whole blocks need equal light: with the distinct rates 1..N the only blocks of one price are the
     four single cells (N = 2..7).
  G3 rooms under unequal light: on the uniform open XY chain at N = 5, phi = (1, 0, -1, 0, 1) (blind
     in F157's sense at the seats 1 and 3, counted from 0), the operator |phi><0...0| is annihilated by
     the full generator with those seats watched at the unequal rates 1 and 7, and moved when seat 0
     is watched (the control); at N = 3 with all rates positive and unequal, (1, 5, 1), |phi><000| with
     phi = (1, 0, -1) is an eigenvector of the full generator at -2, a grey room with mu != 0.
  G4 the shared field's price -2 u^T a u against the full correlated dissipator (N = 3, unequal rates);
     dark coherences inside one-price blocks: under one shared covariance kappa with the noise matrix
     positive semidefinite, every coherence that pays nothing lies in a block of one price (N = 2..7,
     every admissible integer kappa on a grid); control: a neighbour-only covariance on the ring at
     N = 4, a = M M^T exactly, leaves dark coherences in the block (2, 2), whose cells carry four prices;
     and one covariance at the unequal rates (1, 1, 3), N = 3, leaves a dark coherence in (1, 1), whose
     cells carry two prices: the containment needs equal rates.
  G5 amplitude damping: on the full generator with one rate gamma1 on every site and an XXZ chain,
     N = 3, the diagonal block (p, q) is L_H - (gamma1/2)(p + q) Id entry for entry, the block maps
     into (p, q) and (p - 1, q - 1) only, and it is kept by L exactly when min(p, q) = 0; control: the
     rates 1, 2, 3 make L - L_H on (1, 1) a diagonal of several values.
  G6 the parity under the height as a sign: every Hamming distance on (p, q) has the parity of p + q
     (N = 2..7), which is what makes the dissipator's own flow at gamma t = i pi / 2 one sign per block.
  G7 the grey half of the block (1, 1): on the full generator of an XY chain with integer couplings (N = 3, 4,
     5) and of the even ring (N = 4), uniform gamma, every operator odd under U v_(a,b) = (-1)^(a+b) v_(b,a)
     is mapped by L_H to an operator odd under U, and L = L_H - 4 gamma Id on them entry for entry;
     controls: the Heisenberg chain's L_H leaves the odd half, and unequal rates keep it but price it
     with several values. The size: the largest room at -4 gamma, the operators of (1, 1) whose whole
     motion under the single-excitation hopping h has no diagonal, is N^2 minus the rank of the
     ad_h-closure of the diagonal units; that rank is taken modulo two primes, which can only lower it,
     and the odd half already gives C(N, 2), so a count of C(N, 2) at one prime is exact. Uniform open XY
     chain, N = 3..12: exactly C(N, 2). Other chains, where the count differs, each against a lower
     bound: the open Heisenberg chain, N = 3..8, floor(N/2), the frozen modes PROOF_R90_FROZEN_DIVISOR
     section 7 counts (frozen modes have no diagonal and do not move, so they lie in the room); the odd
     ring, which has no U, N - 1 at N = 5, 7, 9, the frozen modes of PROOF_RING_FROZEN_BAND
     Proposition 6; the even ring, N = 4, 6, 8, C(N, 2) plus the U-even frozen modes without diagonal,
     counted exactly by a rational nullspace. The room at 0: the diagonal operators commuting with h
     are the multiples of the identity alone.
  G8 the ladders of PROOF_FROZEN_BAND_SO4 section 2 on the uniform open XY chain, N = 3..6: on every block
     (p, q) off the boundary and every k in its Hamming window, S-^s Phi^r intertwines the Hamiltonian
     part exactly at every step and carries the boundary block (k, 0) onto C(N, k) independent operators
     of Hamming distance k alone; the largest room at -2 gamma k, an upper bound modulo two primes, meets
     that lower bound everywhere except at k = N/2 in the blocks (2, 2) at N = 4 and (3, 2), (4, 3) at
     N = 6, where exact integer Krylov ranks give 8, 22 and 22, two more than the ladders; those two
     need the chain's mirror (couplings the reflection keeps hold them, others fall to the ladder count),
     and at N = 4 they are exhibited from Hamiltonian eigenstates on the mirror-pair configurations. The mirror
     pairs build them at N = 4, 6, 8: the filled mirror orbitals (one excitation on every mirror pair, energies
     +-J_c) against Lambda^dagger|0>, a staggered sum of full mirror pairs at energy 0, of which S+ of the mirror
     R is a multiple; their two products sit at Hamming distance N/2 alone and, with the ladder images, span the
     room of the block (N/2, 2), rank C(N, N/2) + 2, met from above on the uniform chain at N = 8. Only those
     orbitals live wholly on the one-per-pair kets, and only (Lambda^dagger)^m|0> on the empty-or-full ones; the
     products of all m, their transposes and their ladder images, beside the boundary ladder images, span the
     room at k = N/2 wherever it exceeds the ladders at N = 4, 6, and give at least 74 in (4, 4) at N = 8.
  G9 the right half is the bra: at even N = 2M, reading site l <= M as the ket and the mirror site N + 1 - l as
     the bra of site l of the half chain, the open XY chain's H acts on every basis state as
     rho -> H_L rho + rho H_L + J_c flip(rho), flip(rho) = s+_M rho s+_M + s-_M rho s-_M, and on rho Sigma
     (Sigma the product of Z on the odd half-chain sites) as ad_{H_L} + (-1)^M J_c flip, entry for entry at
     N = 4, 6, 8 with mirror-symmetric integer couplings; the even ring adds the same flip at the site 1 times
     the wrap bond (control: couplings the reflection does not keep break the identity). X^M and X^M Z^M, times Sigma, are eigenstates at +-J_c on the
     one-per-pair configurations, the popcount projectors P_j times Sigma are (Lambda^dagger)^j|0> / j!, and on
     the twisted cells Lambda^dagger acts as the eta ladder with the Jordan-Wigner strings counted from the
     centre end, exactly and with one sign, which is PROOF_FROZEN_BAND_SO4's Phi times (-1)^(p + q) on the half
     chain's block (p, q) (control: without the strings it does not). The largest subspace H keeps inside the
     one-per-pair configurations has dimension 2 and inside the empty-or-full ones M + 1, exact rational
     nullspaces, at non-uniform mirror-symmetric couplings; with zero bonds in the half chain (the central bond
     nonzero) the counts are 2 prod (|c| + 1) over the components c without the site M and prod (|c| + 1) over
     all, F4's count over components, met exactly.
  G10 the even ring, two bonds across the axis: the flip sends the letters X, Y, I, Z of its site to X, -Y, 0, 0;
     on rho Sigma the ring with the wrap bond J_w is ad_{H_L} + (-1)^M J_c flip_M - J_w flip_1 entry for entry
     (N = 4, 6, 8, mirror-symmetric, negative and equal-weight cases); X^M and X^M Z^M times Sigma are its
     eigenstates on the one-per-pair configurations, at +-((-1)^M J_c - J_w), and the only ones there except at
     N = 4 with J_w = -J_c, where X^2 (P_0 - P_2) Sigma is a third at 0; the popcount projectors are the only ones
     on the empty-or-full configurations; on the uniform ring the frequency is 0 at even M and 2J at odd M. The room at
     k = N/2 in (M, 2): at N = 6 the ladder images of (M, 0) are not kept and the room is the mirror products
     alone, 2 at the non-uniform mirror-symmetric coupling tried and 6 on the uniform ring (its three bond axes), upper bound modulo two
     primes meeting the exact eigen-products; at N = 4, 8 (non-uniform mirror-symmetric) the ladder images are kept,
     exactly, and with the products fill it, C(N, M) + 2. The uniform ring at N = 4 holds 11 against a span of 9
     from the ladder images and the products of every bond axis, both exact; at N = 8 at most 80 against
     at least 78, whether more not shown.
"""
import numpy as np

FAILS = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(label)


def pc(x):
    return bin(x).count("1")


def boundary(n, p, q):
    return min(p, q) == 0 or max(p, q) == n


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[1, 0], [0, 0]], dtype=complex) @ X        # sigma^- = |0><1|, |1> the excitation (Z|1> = -|1>)


def site(op, l, n):
    out = np.array([[1]], dtype=complex)
    for k in range(n):
        out = np.kron(out, op if k == l else I2)
    return out


def xxz(n, ring, Jxy, Jz, h):
    bonds = [(l, l + 1) for l in range(n - 1)] + ([(n - 1, 0)] if ring else [])
    H = sum(Jxy[b] * (site(X, l, n) @ site(X, m, n) + site(Y, l, n) @ site(Y, m, n))
            + Jz[b] * site(Zm, l, n) @ site(Zm, m, n) for b, (l, m) in enumerate(bonds))
    return H + sum(h[l] * site(Zm, l, n) for l in range(n))


def liouvillian(H, jumps):
    """row-major vec: vec(A rho B) = (A kron B^T) vec(rho)."""
    d = H.shape[0]; Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for r, A in jumps:
        AdA = A.conj().T @ A
        L = L + r * (np.kron(A, A.conj()) - 0.5 * np.kron(AdA, Id) - 0.5 * np.kron(Id, AdA.T))
    return L


def cells(n, p, q):
    return [(a, b) for a in range(2 ** n) for b in range(2 ** n) if pc(a) == p and pc(b) == q]


def block_index(n, p, q):
    return [a * 2 ** n + b for a, b in cells(n, p, q)]


# ---------------- G1 ----------------
print("G1 the boundary blocks")
ok = True
for n in range(2, 8):
    for p in range(n + 1):
        for q in range(n + 1):
            ks = {pc(a ^ b) for a, b in cells(n, p, q)}
            ok &= (len(ks) == 1) == boundary(n, p, q)
check("N = 2..7, every block (p, q): one Hamming distance exactly when min(p, q) = 0 or max(p, q) = N", ok)

rng = np.random.default_rng(7)
ok1, ctl1 = True, True
for n in (3, 4):
    for ring in (False, True):
        nb = n if ring else n - 1
        H = xxz(n, ring, rng.integers(1, 5, nb), rng.integers(-3, 4, nb), rng.integers(-2, 3, n))
        g = 3
        L = liouvillian(H, [(g, site(Zm, l, n)) for l in range(n)])
        LH = liouvillian(H, [])
        for p in range(n + 1):
            for q in range(n + 1):
                idx = block_index(n, p, q)
                rest = np.setdiff1d(np.arange(4 ** n), idx)
                Lb, LHb = L[np.ix_(idx, idx)], LH[np.ix_(idx, idx)]
                if boundary(n, p, q):
                    c = -2 * g * abs(p - q)
                    ok1 &= not np.any(L[np.ix_(rest, idx)])                      # kept by L
                    ok1 &= np.array_equal(Lb - c * np.eye(len(idx)), LHb)        # L - c Id = L_H there
                if (p, q) == (1, 1):
                    Dd = np.diag(Lb - LHb)
                    ctl1 &= np.array_equal(Lb - LHb, np.diag(Dd)) and len(set(Dd.tolist())) > 1
check("full generator, XXZ chain and ring with integer couplings and Z fields, gamma = 3, N = 3, 4: every boundary "
      "block is kept by L and L - c Id = L_H on it entry for entry, c = -2 gamma |p - q|; control: on (1, 1) "
      "L - L_H is diagonal with more than one value, so no c works", ok1 and ctl1)

# ---------------- G2 ----------------
print("G2 whole blocks need equal light")
ok2 = True
for n in range(2, 8):
    rates = list(range(1, n + 1))
    for p in range(n + 1):
        for q in range(n + 1):
            cl = cells(n, p, q)
            prices = {sum(rates[l] for l in range(n) if ((a ^ b) >> l) & 1) for a, b in cl}
            ok2 &= (len(prices) == 1) == (len(cl) == 1)
check("N = 2..7, rates 1..N: a block has one price exactly when it is a single cell, (0,0), (N,N), (0,N), (N,0)", ok2)

# ---------------- G3 ----------------
print("G3 rooms under unequal light")
n = 5
d = 2 ** n
H = xxz(n, False, [1] * (n - 1), [0] * (n - 1), [0] * n)
L = liouvillian(H, [(1, site(Zm, 1, n)), (7, site(Zm, 3, n))])


def room(vec):
    """sum_l vec_l |e_l><0...0|, e_l the one excitation on site l."""
    rho = np.zeros((d, d), dtype=complex)
    for l, v in enumerate(vec):
        rho[1 << (n - 1 - l), 0] = v
    return rho.reshape(-1)


Lc = liouvillian(H, [(1, site(Zm, 0, n))])              # control: a seat where phi does not vanish
check("uniform open XY chain, N = 5, seats counted from 0, dephasing on seats 1 and 3 at rates 1 and 7: "
      "L |phi><0...0| = 0 exactly for phi = (1, 0, -1, 0, 1); control: with seat 0 watched it is moved",
      not np.any(L @ room([1, 0, -1, 0, 1])) and np.any(Lc @ room([1, 0, -1, 0, 1])))

n = 3
d = 2 ** n
H = xxz(n, False, [1] * (n - 1), [0] * (n - 1), [0] * n)
L = liouvillian(H, [(r, site(Zm, l, n)) for l, r in enumerate([1, 5, 1])])
v = room([1, 0, -1])
check("rates positive and not all equal, (1, 5, 1), open XY chain N = 3: L |phi><000| = -2 |phi><000| exactly for "
      "phi = (1, 0, -1), a grey room with mu = -2, the light equal where the room lives",
      np.array_equal(L @ v, -2 * v) and np.any(v))

# ---------------- G4 ----------------
print("G4 dark coherences inside one-price blocks (one shared covariance)")
n = 3
a_k = np.array([[2, 1, 1], [1, 3, 1], [1, 1, 5]], dtype=np.int64)        # rates (2, 3, 5), one covariance 1
Zs = [site(Zm, l, n) for l in range(n)]
dd = 2 ** n
Dk = sum(a_k[k, l] * (np.kron(Zs[k], Zs[l].T) - 0.5 * np.kron(Zs[l] @ Zs[k], np.eye(dd))
                      - 0.5 * np.kron(np.eye(dd), (Zs[l] @ Zs[k]).T)) for k in range(n) for l in range(n))
want = np.zeros(dd * dd)
for x in range(dd):
    for y in range(dd):
        u = np.array([((x >> (n - 1 - l)) & 1) - ((y >> (n - 1 - l)) & 1) for l in range(n)])
        want[x * dd + y] = -2 * int(u @ a_k @ u)
check("the price -2 u^T a u against the full correlated dissipator sum_kl a_kl (Z_k rho Z_l - {Z_l Z_k, rho}/2), "
      "N = 3, rates (2, 3, 5), covariance 1: diagonal exactly, nothing off it", np.array_equal(Dk, np.diag(want)))

ok4, nd = True, 0
for n in range(2, 8):
    G = 6 * (n - 1)                                       # gamma; kappa = gamma eta on an integer grid
    for kap in range(-G // (n - 1), G + 1):               # a = gamma I + kappa (11^T - I) >= 0 iff -gamma/(n-1) <= kappa <= gamma
        for p in range(n + 1):
            for q in range(n + 1):
                cl = cells(n, p, q)
                pr = [(G - kap) * pc(a ^ b) + kap * (p - q) ** 2 for a, b in cl]     # u^T a u; the price is -2 u^T a u
                for (a, b), v in zip(cl, pr):
                    if a != b and v == 0:
                        nd += 1
                        ok4 &= len(set(pr)) == 1
                        ok4 &= kap in (G, -G // (n - 1))          # only at the two ends of the range
check("N = 2..7, gamma = 6(N-1), every integer kappa with a >= 0: every dark coherence lies in a block of one price, "
      "and they occur only at kappa = gamma and at the anti-correlated end",
      ok4 and nd > 0, f"{nd} dark coherences met")

n = 4
Msh = np.eye(n, dtype=np.int64) + np.roll(np.eye(n, dtype=np.int64), 1, axis=1)
a_ring = Msh @ Msh.T                                       # 2 I + neighbours: gamma = 2 kappa, a = M M^T >= 0 exactly


def uau(x, y):
    u = np.array([((x >> l) & 1) - ((y >> l) & 1) for l in range(n)])
    return int(u @ a_ring @ u)


cl = cells(n, 2, 2)
pr = [uau(a, b) for a, b in cl]
dark = [(a, b) for (a, b), v in zip(cl, pr) if a != b and v == 0]
ring_ok = all(a_ring[l, l] == 2 and a_ring[l, (l + 1) % n] == 1 and a_ring[l, (l + 2) % n] == 0 for l in range(n))
check("control: ring N = 4, neighbour-only covariance gamma = 2 kappa: dark coherences in (2, 2), whose cells carry "
      "four prices", ring_ok and len(dark) > 0 and len(set(pr)) == 4, f"dark {dark}, prices {sorted(set(pr))}")

n = 3
a_un = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 3]], dtype=np.int64)  # rates (1, 1, 3), one covariance 1: a >= 0
cl = cells(n, 1, 1)


def uau3(x, y):
    u = np.array([((x >> l) & 1) - ((y >> l) & 1) for l in range(n)])
    return int(u @ a_un @ u)


pr3 = [uau3(x, y) for x, y in cl]
dark3 = [(x, y) for (x, y), v in zip(cl, pr3) if x != y and v == 0]
psd3 = np.array_equal(a_un, np.ones((3, 3), dtype=np.int64) + np.diag([0, 0, 2]))   # 11^T + 2 e3 e3^T >= 0
check("control: the same covariance with unequal rates (1, 1, 3), N = 3: a dark coherence in (1, 1), whose cells "
      "carry two prices, so the containment needs equal rates", psd3 and len(dark3) > 0 and len(set(pr3)) > 1,
      f"dark {dark3}, prices {sorted(set(pr3))}")

# ---------------- G5 ----------------
print("G5 amplitude damping")
n = 3
H = xxz(n, False, [1, 2], [1, -1], [0] * n)
LH = liouvillian(H, [])
blk = {a * 2 ** n + b: (pc(a), pc(b)) for a in range(2 ** n) for b in range(2 ** n)}
ok5, ctl5 = True, True
for rates, uniform in (([2] * n, True), ([1, 2, 3], False)):
    L = liouvillian(H, [(r, site(SM, l, n)) for l, r in enumerate(rates)])
    for p in range(n + 1):
        for q in range(n + 1):
            idx = block_index(n, p, q)
            Lb, LHb = L[np.ix_(idx, idx)], LH[np.ix_(idx, idx)]
            if uniform:
                ok5 &= np.array_equal(Lb, LHb - rates[0] / 2 * (p + q) * np.eye(len(idx)))
                tgt = {blk[i] for i in np.nonzero(np.any(L[:, idx] != 0, axis=1))[0]}
                ok5 &= tgt <= {(p, q), (p - 1, q - 1)}
                ok5 &= ((p - 1, q - 1) not in tgt) == (min(p, q) == 0)      # kept exactly when one side is empty
            elif (p, q) == (1, 1):
                Dd = np.diag(Lb - LHb)
                ctl5 &= np.array_equal(Lb - LHb, np.diag(Dd)) and len(set(Dd.tolist())) > 1
check("full generator, XXZ chain N = 3, damping at gamma1 = 2 on every site: on every (p, q) the diagonal block is "
      "L_H - (gamma1/2)(p + q) Id entry for entry, the block maps into (p, q) and (p - 1, q - 1) only, and it is "
      "kept by L exactly when min(p, q) = 0; control: rates 1, 2, 3 make L - L_H on (1, 1) a diagonal of several "
      "values", ok5 and ctl5)

# ---------------- G6 ----------------
print("G6 the parity under the height as a sign")
ok6 = True
for n in range(2, 8):
    for a in range(2 ** n):
        for b in range(2 ** n):
            ok6 &= pc(a ^ b) % 2 == (pc(a) + pc(b)) % 2
check("N = 2..7: every Hamming distance on (p, q) has the parity of p + q (the parity THE_ONE_SQUARE's gate S8 "
      "meets as well)", ok6)

# ---------------- G7 ----------------
print("G7 the grey half of the block (1, 1)")


def one_exc(n):
    return [1 << (n - 1 - l) for l in range(n)]           # site l carries the excitation -> basis index


def u_odd_basis(n):
    e = one_exc(n); vs = []
    for a in range(n):
        for b in range(a + 1, n):
            v = np.zeros(4 ** n, dtype=complex)
            v[e[a] * 2 ** n + e[b]] = 1
            v[e[b] * 2 ** n + e[a]] = -(-1) ** (a + b)
            vs.append(v)
    return vs


def is_u_odd(n, w):
    e = one_exc(n); c = {(a, b): e[a] * 2 ** n + e[b] for a in range(n) for b in range(n)}
    inside = np.zeros(4 ** n, dtype=bool)
    for i in c.values():
        inside[i] = True
    ok = not np.any(w[~inside])                                     # nothing outside (1, 1)
    ok &= all(w[c[(a, a)]] == 0 for a in range(n))                 # no diagonal
    ok &= all(w[c[(b, a)]] == -(-1) ** (a + b) * w[c[(a, b)]] for a in range(n) for b in range(a + 1, n))
    return ok


rng = np.random.default_rng(11)
ok7, ctl7h, ctl7r = True, True, True
for n in (3, 4, 5):
    for ring in ((False, True) if n % 2 == 0 else (False,)):
        nb = n if ring else n - 1
        H = xxz(n, ring, rng.integers(1, 5, nb), np.zeros(nb, dtype=int), np.zeros(n, dtype=int))
        g = 3
        L = liouvillian(H, [(g, site(Zm, l, n)) for l in range(n)])
        LH = liouvillian(H, [])
        for v in u_odd_basis(n):
            ok7 &= is_u_odd(n, LH @ v) and np.array_equal(L @ v, LH @ v - 4 * g * v)
        if not ring:
            Hc, LHc = H, LH
    Hz = xxz(n, False, np.ones(n - 1, dtype=int), np.ones(n - 1, dtype=int), np.zeros(n, dtype=int))
    LHz = liouvillian(Hz, [])
    ctl7h &= not all(is_u_odd(n, LHz @ v) for v in u_odd_basis(n))
    Lr = liouvillian(Hc, [(r, site(Zm, l, n)) for l, r in enumerate(range(1, n + 1))])
    prices = set()
    for v in u_odd_basis(n):
        w = Lr @ v - LHc @ v
        ctl7r &= is_u_odd(n, Lr @ v)
        nz = np.nonzero(v)[0]
        prices |= {complex(w[i] / v[i]) for i in nz}
    ctl7r &= len(prices) > 1
check("full generator, XY chain (N = 3, 4, 5) and even ring (N = 4) with integer couplings, gamma = 3: L_H keeps "
      "the operators odd under U v_(a,b) = (-1)^(a+b) v_(b,a), and L = L_H - 4 gamma Id on them; controls: the "
      "Heisenberg chain's L_H leaves them at every N tried, and the rates 1..N keep them but with several prices",
      ok7 and ctl7h and ctl7r)


def hop(n, J, delta=0, ring=False):
    bonds = [(l, l + 1) for l in range(n - 1)] + ([(n - 1, 0)] if ring else [])
    h = np.zeros((n, n), dtype=np.int64)
    for (a, b), j in zip(bonds, J):
        h[a, b] += 2 * j; h[b, a] += 2 * j                          # J(XX + YY) on one excitation
        h[a, a] -= 2 * delta; h[b, b] -= 2 * delta                  # J_z ZZ, up to a constant
    return h


def closure_rank(h, p):
    n = len(h); Hm = h % p
    basis, pivots = [], []
    def add(rows):
        kept = []
        for v in rows:
            v = v % p
            for r, c in zip(basis, pivots):
                if v[c]:
                    v = (v - v[c] * r) % p
            nz = np.nonzero(v)[0]
            if len(nz) == 0:
                continue
            c = nz[0]; v = (v * pow(int(v[c]), p - 2, p)) % p
            for i, r in enumerate(basis):
                if r[c]:
                    basis[i] = (r - r[c] * v) % p
            basis.append(v); pivots.append(c); kept.append(v)
        return kept
    front = add([np.eye(n * n, dtype=np.int64)[l * n + l] for l in range(n)])
    while front:
        F = np.array(front).reshape(-1, n, n)
        A = (np.einsum('ik,mkj->mij', Hm, F) - np.einsum('mik,kj->mij', F, Hm)) % p
        front = add(list(A.reshape(len(front), n * n)))
    return len(pivots)


def grey_room(h):
    return min(len(h) ** 2 - closure_rank(h, p) for p in (32749, 32719))


ok7s = all(np.array_equal(hop(n, np.ones(n - 1, dtype=int)),
                          xxz(n, False, np.ones(n - 1), np.zeros(n - 1), np.zeros(n))[np.ix_(one_exc(n), one_exc(n))].real
                          ) for n in (3, 4, 5))                     # the hopping matrix is the full H's (1-excitation) block
sizes = {n: grey_room(hop(n, np.ones(n - 1, dtype=int))) for n in range(3, 13)}
ok7s &= all(sizes[n] == n * (n - 1) // 2 for n in sizes)
check("the largest room at -4 gamma in (1, 1), uniform open XY chain, N = 3..12: exactly C(N, 2), the odd half "
      "and nothing more", ok7s, f"{list(sizes.values())}")

import sympy as sp


def _dm_rank(A):
    from sympy.polys.matrices import DomainMatrix
    return DomainMatrix.from_Matrix(A).convert_to(sp.QQ).rank()      # exact rank over the rationals


def even_frozen(h):
    """U-even operators of (1, 1) with no diagonal that commute with h: an exact rational nullspace"""
    n = len(h); Hs = sp.Matrix(h); S = sp.diag(*[(-1) ** l for l in range(n)])
    Xs = sp.Matrix(n, n, sp.symbols(f"x0:{n * n}"))
    eqs = list(Hs * Xs - Xs * Hs) + [Xs[i, i] for i in range(n)] + list(S * Xs.T * S - Xs)
    A, _ = sp.linear_eq_to_matrix(eqs, list(Xs))
    return n * n - _dm_rank(A)


def frozen_nodiag(h):
    """operators of (1, 1) with no diagonal that commute with h: an exact rational nullspace, a lower bound on the
    room at -4 gamma (they do not move and the watching is -4 gamma on them)"""
    n = len(h); Hs = sp.Matrix(h)
    Xs = sp.Matrix(n, n, sp.symbols(f"y0:{n * n}"))
    A, _ = sp.linear_eq_to_matrix(list(Hs * Xs - Xs * Hs) + [Xs[i, i] for i in range(n)], list(Xs))
    return n * n - _dm_rank(A)


heis = {n: grey_room(hop(n, np.ones(n - 1, dtype=int), delta=1)) for n in range(3, 9)}
odd = {n: grey_room(hop(n, np.ones(n, dtype=int), ring=True)) for n in (5, 7, 9)}
ok7f = all(frozen_nodiag(hop(n, np.ones(n - 1, dtype=int), delta=1)) == heis[n] for n in heis)
ok7f &= all(frozen_nodiag(hop(n, np.ones(n, dtype=int), ring=True)) == odd[n] for n in odd)
even = {n: (grey_room(hop(n, np.ones(n, dtype=int), ring=True)), even_frozen(hop(n, np.ones(n, dtype=int), ring=True)))
        for n in (4, 6, 8)}
check("other chains, each upper bound meeting a lower bound: open Heisenberg chain, N = 3..8, floor(N/2), the frozen "
      "modes of PROOF_R90_FROZEN_DIVISOR section 7; odd XY ring, N = 5, 7, 9, N - 1, those of PROOF_RING_FROZEN_BAND "
      "Proposition 6; even XY ring, N = 4, 6, 8, C(N, 2) plus the U-even frozen modes without diagonal",
      ok7f and all(heis[n] == n // 2 for n in heis) and all(odd[n] == n - 1 for n in odd)
      and all(r == n * (n - 1) // 2 + e for n, (r, e) in even.items()),
      f"Heisenberg {list(heis.values())}, odd ring {odd}, even ring (room, U-even frozen) {even}")

ok70 = True
for n in range(3, 13):
    h = hop(n, np.ones(n - 1, dtype=int))
    # [h, D]_ab = (d_b - d_a) h_ab, so D commutes with h exactly when d is constant on every component of
    # the graph of h's nonzero off-diagonal entries: the dimension is the number of components
    comp = list(range(n))
    def root(x):
        while comp[x] != x:
            x = comp[x]
        return x
    for a in range(n):
        for b in range(n):
            if a != b and h[a, b]:
                comp[root(a)] = root(b)
    ok70 &= len({root(x) for x in range(n)}) == 1
check("the room at 0 in (1, 1), uniform open XY chain, N = 3..12: the diagonal operators commuting with the "
      "hopping are the multiples of the identity", ok70)

# ---------------- G8 ----------------
print("G8 the ladders carry the boundary blocks into the coloured blocks")
from math import comb


def subsets(n, p):
    return [A for A in range(2 ** n) if bin(A).count("1") == p]


def h_sector(n, p, J=1):
    """the uniform open XY chain on p excitations, cells as site sets (bit l = site l), integer"""
    bs = subsets(n, p); ix = {b: i for i, b in enumerate(bs)}
    Hs = np.zeros((len(bs), len(bs)), dtype=np.int64)
    for A in bs:
        for l in range(n - 1):
            if ((A >> l) & 1) != ((A >> (l + 1)) & 1):
                Hs[ix[A ^ (3 << l)], ix[A]] += J
    return Hs, bs


def turn(n, p, q):
    """K on the block (p, q): rho -> H_p rho - rho H_q, row-major, integer (L_H = -i K)"""
    Hp, _ = h_sector(n, p); Hq, _ = h_sector(n, q)
    return np.kron(Hp, np.eye(len(Hq), dtype=np.int64)) - np.kron(np.eye(len(Hp), dtype=np.int64), Hq.T)


def jw(A, l):
    return -1 if bin(A & ((1 << l) - 1)).count("1") % 2 else 1


def phi(n, p, q):
    """Phi(rho) = sum_l d+_l rho d_l : (p, q) -> (p + 1, q + 1)"""
    bp, bq, bP, bQ = subsets(n, p), subsets(n, q), subsets(n, p + 1), subsets(n, q + 1)
    iP = {b: i for i, b in enumerate(bP)}; iQ = {b: i for i, b in enumerate(bQ)}
    M = np.zeros((len(bP) * len(bQ), len(bp) * len(bq)), dtype=np.int64)
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(n):
                if not (A >> l) & 1 and not (B >> l) & 1:
                    M[iP[A | 1 << l] * len(bQ) + iQ[B | 1 << l], i * len(bq) + j] += jw(A, l) * jw(B, l)
    return M


def splus(n, p, q):
    """S+(rho) = sum_l (-1)^l d+_l rho d+_l : (p, q) -> (p + 1, q - 1); S- is its transpose"""
    bp, bq, bP, bQ = subsets(n, p), subsets(n, q), subsets(n, p + 1), subsets(n, q - 1)
    iP = {b: i for i, b in enumerate(bP)}; iQ = {b: i for i, b in enumerate(bQ)}
    M = np.zeros((len(bP) * len(bQ), len(bp) * len(bq)), dtype=np.int64)
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(n):
                if not (A >> l) & 1 and (B >> l) & 1:
                    M[iP[A | 1 << l] * len(bQ) + iQ[B ^ 1 << l], i * len(bq) + j] += jw(A, l) * jw(B ^ 1 << l, l) * (-1) ** l
    return M


def rank_mod(M, P):
    M = M.copy() % P; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), P - 2, P)) % P
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % P
        r += 1
    return r


def closure_room(n, p, q, k, P):
    Hp, bp = h_sector(n, p); Hq, bq = h_sector(n, q)
    T = turn(n, p, q).T % P                                 # dual action on functionals
    basis, piv = np.zeros((0, len(bp) * len(bq)), dtype=np.int64), []
    def add(rows):
        nonlocal basis
        kept = []
        for v in rows:
            v = v % P
            if piv:
                v = (v - v[piv] @ basis) % P
            nz = np.nonzero(v)[0]
            if not len(nz):
                continue
            c = nz[0]; v = (v * pow(int(v[c]), P - 2, P)) % P
            if piv:
                basis = (basis - np.outer(basis[:, c], v)) % P
            basis = np.vstack([basis, v]); piv.append(c); kept.append(v)
        return kept
    seeds = []
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            if bin(A ^ B).count("1") != k:
                e = np.zeros(len(bp) * len(bq), dtype=np.int64); e[i * len(bq) + j] = 1; seeds.append(e)
    front = add(seeds)
    while front:
        front = add([(T @ v) % P for v in front])
    return len(bp) * len(bq) - len(piv)


P1, P2 = 32749, 32719
EXCESS = {(4, 2, 2, 2): 8, (6, 3, 2, 3): 22, (6, 4, 3, 3): 22}             # the room above the ladders at k = N/2
ok8i, ok8b, ok8u, seen = True, True, True, {}
for n in (3, 4, 5, 6):
    for p in range(n + 1):
        for q in range(p + 1):
            if q == 0 or p == n:
                continue                                    # boundary blocks: G1
            for k in range(p - q, min(p + q, 2 * n - p - q) + 1, 2):
                s = (k - p + q) // 2; r = q - s             # S-^s then Phi^r from (k, 0)
                V = np.eye(comb(n, k), dtype=np.int64); a, b = k, 0
                Kk = turn(n, k, 0)
                for _ in range(s):
                    Sm = splus(n, a - 1, b + 1).T           # (a, b) -> (a - 1, b + 1)
                    ok8i &= np.array_equal(turn(n, a - 1, b + 1) @ Sm, Sm @ turn(n, a, b))
                    V = Sm @ V; a, b = a - 1, b + 1
                for _ in range(r):
                    F = phi(n, a, b)
                    ok8i &= np.array_equal(turn(n, a + 1, b + 1) @ F, F @ turn(n, a, b))
                    V = F @ V; a, b = a + 1, b + 1
                assert (a, b) == (p, q)
                bp, bq = subsets(n, p), subsets(n, q)
                ham = np.array([bin(A ^ B).count("1") for A in bp for B in bq])
                ok8b &= not np.any(V[ham != k]) and rank_mod(V.T, P1) == comb(n, k)   # one height, full rank
                room = min(closure_room(n, p, q, k, P) for P in (P1, P2))
                key = (n, max(p, q), min(p, q), k)
                ok8u &= room == EXCESS.get(key, comb(n, k))
                if key in EXCESS:
                    seen[key] = room
ok8u &= set(seen) == set(EXCESS)
check("uniform open XY chain, N = 3..6, every block (p, q) off the boundary and every k in its Hamming window: "
      "S-^s Phi^r intertwines K exactly on every step, and carries the boundary block (k, 0) onto C(N, k) "
      "independent operators of Hamming distance k alone, a grey room at -2 gamma k", ok8i and ok8b)
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ


def exact_room(n, p, q, k):
    """the room over the integers: the Krylov row space of the cell functionals off distance k under K, ranked
    exactly (DomainMatrix over ZZ) until it stops growing"""
    _, bp = h_sector(n, p); _, bq = h_sector(n, q)
    K = turn(n, p, q).astype(object)
    C = np.array([np.eye(len(bp) * len(bq), dtype=object)[i * len(bq) + j] for i, A in enumerate(bp)
                  for j, B in enumerate(bq) if bin(A ^ B).count("1") != k], dtype=object)
    blocks, r_old = [C], -1
    while True:
        M = np.vstack(blocks)
        r = DomainMatrix([[ZZ(int(x)) for x in row] for row in M], M.shape, ZZ).rank()
        if r == r_old:
            return C.shape[1] - r
        r_old = r
        blocks.append(blocks[-1] @ K)


exact = {key: exact_room(key[0], key[1], key[2], key[3]) for key in EXCESS}
check("the same blocks: the largest room at -2 gamma k is C(N, k) (an upper bound modulo two primes meeting the "
      "ladders' lower bound) everywhere except at k = N/2 strictly inside the window, where exact integer ranks give "
      "two more", ok8u and exact == EXCESS, f"exact {exact} against {[comb(k[0], k[3]) for k in exact]}")

def h_sector_bonds(n, p, J):
    bs = subsets(n, p); ix = {b: i for i, b in enumerate(bs)}
    Hs = np.zeros((len(bs), len(bs)), dtype=np.int64)
    for A in bs:
        for l in range(n - 1):
            if ((A >> l) & 1) != ((A >> (l + 1)) & 1):
                Hs[ix[A ^ (3 << l)], ix[A]] += J[l]
    return Hs, bs


def room_bonds(n, p, q, k, J, P):
    """closure_room with the couplings J (the ladders bound it from below by C(N, k) for any real J)"""
    Hp, bp = h_sector_bonds(n, p, J); Hq, bq = h_sector_bonds(n, q, J)
    T = (np.kron(Hp, np.eye(len(bq), dtype=np.int64)) - np.kron(np.eye(len(bp), dtype=np.int64), Hq.T)).T % P
    seeds = [np.eye(len(bp) * len(bq), dtype=np.int64)[i * len(bq) + j]
             for i, A in enumerate(bp) for j, B in enumerate(bq) if bin(A ^ B).count("1") != k]
    span = np.zeros((0, len(bp) * len(bq)), dtype=np.int64); piv = []
    front = seeds
    while front:
        new = []
        for v in front:
            v = v % P
            if piv:
                v = (v - v[piv] @ span) % P
            nz = np.nonzero(v)[0]
            if not len(nz):
                continue
            c = nz[0]; v = (v * pow(int(v[c]), P - 2, P)) % P
            if piv:
                span = (span - np.outer(span[:, c], v)) % P
            span = np.vstack([span, v]); piv.append(c); new.append(v)
        front = [(T @ v) % P for v in new]
    return len(bp) * len(bq) - len(piv)


mirror = {tuple(J): min(room_bonds(4, 2, 2, 2, J, P) for P in (P1, P2)) for J in ([1, 3, 1], [2, 1, 2], [1, 2, 3], [1, 1, 2])}
mirror6 = {tuple(J): min(room_bonds(6, 3, 2, 3, J, P) for P in (P1, P2)) for J in ([1, 3, 2, 3, 1], [1, 2, 1, 1, 1])}
ok8m = all(v == (8 if list(J) == list(J)[::-1] else 6) for J, v in mirror.items())
ok8m &= all(v == (22 if list(J) == list(J)[::-1] else 20) for J, v in mirror6.items())


def exact_room_bonds(n, p, q, k, J):
    _, bp = h_sector_bonds(n, p, J); _, bq = h_sector_bonds(n, q, J)
    Hp, _ = h_sector_bonds(n, p, J); Hq, _ = h_sector_bonds(n, q, J)
    K = (np.kron(Hp, np.eye(len(bq), dtype=np.int64)) - np.kron(np.eye(len(bp), dtype=np.int64), Hq.T)).astype(object)
    C = np.array([np.eye(len(bp) * len(bq), dtype=object)[i * len(bq) + j] for i, A in enumerate(bp)
                  for j, B in enumerate(bq) if bin(A ^ B).count("1") != k], dtype=object)
    blocks, r_old = [C], -1
    while True:
        M = np.vstack(blocks)
        r = DomainMatrix([[ZZ(int(x)) for x in row] for row in M], M.shape, ZZ).rank()
        if r == r_old:
            return C.shape[1] - r
        r_old = r
        blocks.append(blocks[-1] @ K)


# the counts above the ladders are upper bounds; the mirror-symmetric ones are ranked exactly as well
ok8m &= exact_room_bonds(4, 2, 2, 2, [1, 3, 1]) == 8 and exact_room_bonds(4, 2, 2, 2, [2, 1, 2]) == 8
ok8m &= exact_room_bonds(6, 3, 2, 3, [1, 3, 2, 3, 1]) == 22

# exact at N = 4: the six ladder images plus two operators |v><phi| built from Hamiltonian eigenvectors, v on the
# configurations with one excitation in each mirror pair (sites 0,3 and 1,2), phi = |0110> - |1001> on those with
# each pair empty or full; every such cell has Hamming distance 2. Rank 8 modulo a prime is a lower bound over Q.
ok8x = True
for J4 in ([1, 1, 1], [1, 3, 1]):
    H2, b2 = h_sector_bonds(4, 2, J4)
    pair_one = [i for i, A in enumerate(b2) if ((A >> 0) & 1) != ((A >> 3) & 1) and ((A >> 1) & 1) != ((A >> 2) & 1)]
    phi4 = np.zeros(len(b2), dtype=np.int64); phi4[b2.index(0b0110)] = 1; phi4[b2.index(0b1001)] = -1
    ok8x &= np.array_equal(H2 @ phi4, 0 * phi4)
    vs = []
    for e in (J4[1], -J4[1]):                             # H2 v = e v, e = +- the central bond, v on one-per-pair cells
        M = sp.Matrix(H2[:, pair_one]) - e * sp.Matrix(np.eye(len(b2), dtype=np.int64)[:, pair_one])
        ns = M.nullspace()
        ok8x &= len(ns) == 1
        w = ns[0] * sp.ilcm(*[x.q for x in ns[0]])
        v = np.zeros(len(b2), dtype=np.int64); v[pair_one] = [int(x) for x in w]
        ok8x &= np.array_equal(H2 @ v, e * v)
        vs.append(np.outer(v, phi4).ravel())
    ham4 = np.array([bin(A ^ B).count("1") for A in b2 for B in b2])
    ok8x &= all(not np.any(X[ham4 != 2]) for X in vs)
    Vlad = phi(4, 1, 1) @ (splus(4, 1, 1).T @ np.eye(comb(4, 2), dtype=np.int64))   # S- of (2, 0), then Phi
    def turn_b(p, q):
        Hp, _ = h_sector_bonds(4, p, J4); Hq, _ = h_sector_bonds(4, q, J4)
        return np.kron(Hp, np.eye(len(Hq), dtype=np.int64)) - np.kron(np.eye(len(Hp), dtype=np.int64), Hq.T)
    ok8x &= np.array_equal(turn_b(2, 2) @ Vlad, Vlad @ turn_b(2, 0))      # the ladder images are kept at these couplings
    ok8x &= rank_mod(np.vstack([Vlad.T] + [x[None, :] for x in vs]), P1) == 8
check("the excess at k = N/2 needs the chain's mirror: couplings the reflection keeps hold it ((1, 3, 1), (2, 1, 2) at "
      "N = 4: 8; (1, 3, 2, 3, 1) at N = 6: 22), couplings it does not keep fall to the ladder count ((1, 2, 3), (1, 1, 2): "
      "6; (1, 2, 1, 1, 1): 20), the mirror-symmetric counts also as exact integer ranks; at N = 4, couplings (1, 1, 1) "
      "and (1, 3, 1), the six ladder images and two operators |v><phi| with H v = +-J_c v (J_c the central bond) on "
      "one excitation per mirror pair and phi = |0110> - |1001>, every cell at Hamming distance 2, span 8",
      ok8m and ok8x, f"{mirror} {mirror6}")

# the excess built from the mirror pairs: the filled mirror orbitals Psi_e, Psi_o (one excitation on every mirror pair,
# energies +-J_c) against Lambda^dagger|0> = sum_l (-1)^l |{l, N-1-l}> (a whole mirror pair, staggered, energy 0)
def mirror_seeds(n, J):
    M = n // 2; Jc = J[M - 1]
    Hm, bm = h_sector_bonds(n, M, J); H2, b2 = h_sector_bonds(n, 2, J)
    lam = np.zeros(len(b2), dtype=np.int64)
    for l in range(M):
        lam[b2.index((1 << l) | (1 << (n - 1 - l)))] = (-1) ** l
    ok = np.array_equal(H2 @ lam, 0 * lam)
    O = [i for i, A in enumerate(bm) if all(((A >> l) & 1) != ((A >> (n - 1 - l)) & 1) for l in range(M))]
    seeds = []
    for e in (Jc, -Jc):
        Mx = sp.Matrix(Hm[:, O]) - e * sp.Matrix(np.eye(len(bm), dtype=np.int64)[:, O])
        nsp = Mx.nullspace()
        ok &= len(nsp) == 1                               # one eigenstate at +-J_c living on one-per-pair kets
        w = nsp[0] * sp.ilcm(*[x.q for x in nsp[0]])
        v = np.zeros(len(bm), dtype=np.int64); v[O] = [int(x) for x in w]
        ok &= np.array_equal(Hm @ v, e * v)
        seeds.append(np.outer(v, lam).ravel())
    ham = np.array([bin(A ^ B).count("1") for A in bm for B in b2])
    ok &= all(not np.any(X[ham != M]) for X in seeds)     # every mirror pair disagrees once: height N/2
    V = phi(n, M - 1, 1) @ (splus(n, M - 1, 1).T @ np.eye(comb(n, M), dtype=np.int64))   # ladder images in (M, 2)
    H1, _ = h_sector_bonds(n, 1, J); H0, _ = h_sector_bonds(n, 0, J)
    Kt = np.kron(Hm, np.eye(len(b2), dtype=np.int64)) - np.kron(np.eye(len(bm), dtype=np.int64), H2.T)
    Ks = np.kron(Hm, np.eye(1, dtype=np.int64)) - np.kron(np.eye(len(bm), dtype=np.int64), H0.T)
    ok &= np.array_equal(Kt @ V, V @ Ks)                  # the ladder images are kept at these couplings
    b1 = subsets(n, 1); Rm = np.zeros((n, n), dtype=np.int64)  # the mirror as an operator of (1, 1)
    for i, a in enumerate(b1):
        Rm[b1.index(1 << (n - 1 - (a.bit_length() - 1))), i] = 1
    ok &= np.array_equal(H1 @ Rm, Rm @ H1) and not np.any(np.diag(Rm))           # frozen in (1, 1), no diagonal
    ok &= np.array_equal(splus(n, 1, 1) @ Rm.ravel(), 2 * lam)                    # S+ carries it onto Lambda|0><0|
    rank = rank_mod(np.vstack([V.T] + [x[None, :] for x in seeds]), P1)
    return ok, rank


mirror_rk = {(n, tuple(J)): mirror_seeds(n, J) for n, J in [(4, [1, 1, 1]), (4, [1, 3, 1]), (6, [1] * 5), (6, [2, 1, 3, 1, 2]),
                                                     (8, [1] * 7), (8, [1, 2, 1, 3, 1, 2, 1])]}
ok8f = all(o and r == comb(n, n // 2) + 2 for (n, _), (o, r) in mirror_rk.items())
up8 = min(closure_room(8, 4, 2, 4, P) for P in (P1, P2))     # the room of (4, 2) at N = 8 from above
check("the excess built from the mirror pairs, N = 4, 6, 8 (uniform and mirror-symmetric couplings): Lambda^dagger|0>, a "
      "staggered sum of full mirror pairs, has energy 0, and the spin ladder carries the mirror R, frozen in (1, 1) with "
      "no diagonal, onto 2 Lambda^dagger|0><0|; one eigenstate at +J_c and one at -J_c lives wholly on the kets "
      "with one excitation per mirror pair; the two products sit at Hamming distance N/2 alone and are independent "
      "of the C(N, N/2) ladder images in the block (N/2, 2), rank C(N, N/2) + 2; on the uniform chain at N = 8 the upper bound meets it, 72",
      ok8f and up8 == comb(8, 4) + 2, f"ranks {[r for _, r in mirror_rk.values()]}, upper bound N = 8: {up8}")

def kept_inside(n, p, inset, P):
    """largest subspace the p-excitation hopping keeps inside span{configs in the set}: dimension minus the rank of
    the H-closure of the unit functionals off the set (mod P, an upper bound on the dimension)"""
    H, b = h_sector_bonds(n, p, [1] * (n - 1)); H = H % P
    piv, basis = [], np.zeros((0, len(b)), dtype=np.int64)
    front = [np.eye(len(b), dtype=np.int64)[i] for i, s in enumerate(b) if not inset(s)]
    while front:
        new = []
        for v in front:
            v = v % P
            if piv:
                v = (v - v[piv] @ basis) % P
            nz = np.nonzero(v)[0]
            if not len(nz):
                continue
            c = nz[0]; v = (v * pow(int(v[c]), P - 2, P)) % P
            if piv:
                basis = (basis - np.outer(basis[:, c], v)) % P
            basis = np.vstack([basis, v]); piv.append(c); new.append(v)
        front = [(H @ v) % P for v in new]
    return len(b) - len(piv)


def move(n, X, a, b, A, B):
    """carry the operators (columns of X) of the block (a, b) to (A, B) by powers of S+-, Phi and its adjoint"""
    ds, de = (A - B) - (a - b), (A + B) - (a + b)
    if ds % 2 or de % 2:
        return None
    for _ in range(abs(ds) // 2):
        if ds > 0:
            if b < 1:
                return None
            X = splus(n, a, b) @ X; a, b = a + 1, b - 1
        else:
            if a < 1:
                return None
            X = splus(n, a - 1, b + 1).T @ X; a, b = a - 1, b + 1
    for _ in range(abs(de) // 2):
        if de > 0:
            if a + 1 > n or b + 1 > n:
                return None
            X = phi(n, a, b) @ X; a, b = a + 1, b + 1
        else:
            if a < 1 or b < 1:
                return None
            X = phi(n, a - 1, b - 1).T @ X; a, b = a - 1, b - 1
    return X


def pair_products(n):
    """|Psi_eo><Lambda^m 0| and the transposes, m = 1..N/2 - 1, uniform chain, as (operator, block)"""
    M = n // 2
    one = lambda s: all(((s >> l) & 1) != ((s >> (n - 1 - l)) & 1) for l in range(M))
    sym = lambda s: all(((s >> l) & 1) == ((s >> (n - 1 - l)) & 1) for l in range(M))
    def eig_on(p, inset, e):
        H, b = h_sector_bonds(n, p, [1] * (n - 1)); S = [i for i, s in enumerate(b) if inset(s)]
        m = sp.Matrix(H[:, S]) - e * sp.Matrix(np.eye(len(b), dtype=np.int64)[:, S]); out = []
        for v in m.nullspace():
            w = v * sp.ilcm(*[x.q for x in v]); x = np.zeros(len(b), dtype=np.int64); x[S] = [int(t) for t in w]; out.append(x)
        return out
    psis = eig_on(M, one, 1) + eig_on(M, one, -1)
    out = []
    for m in range(1, M):
        lam = eig_on(2 * m, sym, 0)
        if len(lam) != 1:
            return None, psis
        for v in psis:
            out += [(np.outer(v, lam[0]).ravel(), (M, 2 * m)), (np.outer(lam[0], v).ravel(), (2 * m, M))]
    return out, psis


ok8u2, ok8b2 = True, True
for n in (4, 6, 8):
    _, ps = pair_products(n)
    for v in ps:                                          # Psi_e, Psi_o as integer amplitude tensors on 2^N
        _, bm = h_sector_bonds(n, n // 2, [1] * (n - 1))
        T = np.zeros(2 ** n, dtype=np.int64)
        for amp, s in zip(v, bm):
            T[s] = amp
        T = T.reshape([2] * n)                             # axis k is bit n - 1 - k
        for l in range(n // 2):
            ax = [n - 1 - l, l]
            A = np.moveaxis(T, ax, [0, 1]).reshape(4, -1)
            ok8b2 &= bool(np.any(A)) and all(not np.any(np.outer(A[i], A[j]) - np.outer(A[j], A[i]))
                                             for i in range(4) for j in range(i + 1, 4))   # rank 1 exactly: every 2x2 minor 0
            nz = np.abs(A[A != 0]); ok8b2 &= bool(np.all(nz == nz[0]))     # and equal moduli: Bell pairs
    M = n // 2
    ok8u2 &= kept_inside(n, M, lambda s: all(((s >> l) & 1) != ((s >> (n - 1 - l)) & 1) for l in range(M)), P1) == 2
    ok8u2 &= all(kept_inside(n, 2 * m, lambda s: all(((s >> l) & 1) == ((s >> (n - 1 - l)) & 1) for l in range(M)), P1) == 1
                 for m in range(M + 1))
check("only the two filled mirror orbitals live wholly on the kets with one excitation per mirror pair, and only "
      "(Lambda^dagger)^m|0> on those with every pair empty or full (uniform chain, N = 4, 6, 8: the largest subspace the "
      "hopping keeps there has dimension at most 2 and 1, modulo a prime; Psi_e, Psi_o and the staggered pairs, built "
      "in the checks around it, reach it), and Psi_e, Psi_o are products of one Bell state on every mirror pair "
      "(integer amplitudes, rank 1 exactly between each pair and the rest, equal moduli)", ok8u2 and ok8b2)

span_ok, spans = True, {}
for n, blocks in ((4, [(2, 2)]), (6, [(3, 2), (2, 3), (4, 3), (3, 4)]), (8, [(3, 3), (4, 4)])):
    prods, _ = pair_products(n)
    span_ok &= prods is not None
    for (p, q) in blocks:
        vecs = list(move(n, np.eye(comb(n, n // 2), dtype=np.int64), n // 2, 0, p, q).T)
        for X, (a, b) in prods:
            Y = move(n, X[:, None], a, b, p, q)
            if Y is not None and np.any(Y):
                vecs.append(Y[:, 0])
        V = np.array(vecs)
        _, bp = h_sector(n, p); _, bq = h_sector(n, q)
        ham = np.array([bin(A ^ B).count("1") for A in bp for B in bq])
        span_ok &= not np.any(V[:, ham != n // 2])                       # every vector at distance N/2 alone
        r = rank_mod(V, P1)
        span_ok &= rank_mod(np.vstack([V, (turn(n, p, q) @ V.T).T]), P1) == r   # and the span kept by the turning
        spans[(n, p, q)] = r
        target = EXCESS.get((n, max(p, q), min(p, q), n // 2)) if n < 8 else comb(8, 4) + (4 if p == q == 4 else 2)
        span_ok &= r == target
check("on the uniform chain the ladder images of the boundary block (N/2, 0), with the mirror-pair products and their "
      "ladder images, lie at distance N/2 alone, span a subspace consistent with being kept by the turning (modulo a prime; "
      "kept exactly by Lemmas 2.2, 2.3 and the exact eigenvectors), and fill the room at k = N/2 where it "
      "exceeds the ladders: every such block at N = 4, 6 (the room from above by the closure ranks of G8); at N = 8 at "
      "least C(8, 4) + 2 = 72 in (3, 3) and C(8, 4) + 4 = 74 in (4, 4)",
      span_ok, f"{spans}")


print("G9 the right half is the bra")
# spin language throughout, site 1 the most significant bit; the half chain has M = N/2 sites
_sp = np.array([[0, 0], [1, 0]], dtype=np.int64); _sm = _sp.T.copy()
_I2 = np.eye(2, dtype=np.int64); _Z = np.diag([1, -1]).astype(np.int64); _X = np.array([[0, 1], [1, 0]], dtype=np.int64)


def _kr(ms):
    out = np.ones((1, 1), dtype=np.int64)
    for m_ in ms:
        out = np.kron(out, m_)
    return out


def _site(n, l, op):
    return _kr([op if k == l else _I2 for k in range(1, n + 1)])


def _hop(n, J):
    H = np.zeros((2 ** n, 2 ** n), dtype=np.int64)
    for l in range(1, n):
        a = _site(n, l, _sp) @ _site(n, l + 1, _sm)
        H += J[l - 1] * (a + a.T)
    return H


def _cell(idx, n):
    """the basis state idx of N sites as the cell (a, b) of the half chain: ket = sites 1..M, bra = sites N..M+1"""
    m = n // 2
    bits = [(idx >> (n - k)) & 1 for k in range(1, n + 1)]
    return (int("".join(str(bits[l - 1]) for l in range(1, m + 1)), 2),
            int("".join(str(bits[n - l]) for l in range(1, m + 1)), 2))


def _to_op(v, n):
    m = n // 2; rho = np.zeros((2 ** m, 2 ** m), dtype=np.int64)
    for idx in range(2 ** n):
        a, b = _cell(idx, n); rho[a, b] = v[idx]
    return rho


def _to_vec(rho, n):
    return np.array([rho[_cell(idx, n)] for idx in range(2 ** n)], dtype=np.int64)


def _kept(H, S):
    """largest H-invariant subspace inside the coordinate span of S, exact rational arithmetic"""
    n = H.shape[0]; Hs = sp.Matrix(H)
    W = sp.zeros(n, len(S))
    for c, i in enumerate(S):
        W[i, c] = 1
    while W.shape[1]:
        A = (sp.eye(n) - W * (W.T * W).inv() * W.T) * Hs * W
        ns = A.nullspace()
        if len(ns) == W.shape[1]:
            break
        W = sp.Matrix.hstack(*[W * v for v in ns]) if ns else sp.zeros(n, 0)
    return W.shape[1]


def _edges(n):
    m = n // 2
    one = [i for i in range(2 ** n) if _cell(i, n)[0] ^ _cell(i, n)[1] == 2 ** m - 1]   # one excitation per mirror pair
    emf = [i for i in range(2 ** n) if _cell(i, n)[0] == _cell(i, n)[1]]               # every pair empty or full
    return one, emf


from math import factorial
ok9a = ok9b = ok9c = ok9d = True
ctl9a = ctl9d = True
for n, J in ((4, [1, 2, 1]), (6, [1, 3, 2, 3, 1]), (6, [1] * 5), (8, [2, 1, 3, 5, 3, 1, 2])):
    m = n // 2; HN = _hop(n, J); HL = _hop(m, J[:m - 1]); Jc = J[m - 1]
    Sig = _kr([_Z if l % 2 else _I2 for l in range(1, m + 1)])
    spM, smM = _site(m, m, _sp), _site(m, m, _sm)
    F = lambda r: spM @ r @ spM + smM @ r @ smM
    sgn = (-1) ** m
    for idx in range(2 ** n):
        e = np.zeros(2 ** n, dtype=np.int64); e[idx] = 1
        r = _to_op(e, n); img = _to_op(HN @ e, n)
        ok9a &= np.array_equal(img, HL @ r + r @ HL + Jc * F(r))
        rt = r @ Sig
        ok9a &= np.array_equal(img @ Sig, HL @ rt - rt @ HL + sgn * Jc * F(rt))
    XM = _kr([_X] * m); ZM = _kr([_Z] * m)
    for C, f in ((XM, 1), (XM @ ZM, -1)):
        ok9b &= not np.any(HL @ C - C @ HL)
        v = _to_vec(C @ Sig, n)
        ok9b &= np.array_equal(HN @ v, sgn * f * Jc * v)
        ok9b &= all((a ^ b) == 2 ** m - 1 for a, b in np.argwhere(C))
    Lam = np.zeros((2 ** n, 2 ** n), dtype=np.int64)
    for l in range(1, m + 1):
        r_ = n + 1 - l; s_ = np.eye(2 ** n, dtype=np.int64)
        for k in range(l + 1, r_):
            s_ = s_ @ _site(n, k, _Z)
        Lam += (-1) ** l * _site(n, l, _sp) @ s_ @ _site(n, r_, _sp)
    ok9c &= not np.any(HN @ Lam - Lam @ HN)
    T = [_kr([_Z if k > l else _I2 for k in range(1, m + 1)]) for l in range(1, m + 1)]   # strings from the centre end
    S = [_kr([_Z if k < l else _I2 for k in range(1, m + 1)]) for l in range(1, m + 1)]  # SO4's strings, k < l
    sgns, plain, so4 = set(), set(), True
    for idx in range(2 ** n):
        e = np.zeros(2 ** n, dtype=np.int64); e[idx] = 1
        rt = _to_op(e, n) @ Sig; img = _to_op(Lam @ e, n) @ Sig
        eta = sum(T[l - 1] @ _site(m, l, _sp) @ rt @ _site(m, l, _sm) @ T[l - 1] for l in range(1, m + 1))
        bare = sum(_site(m, l, _sp) @ rt @ _site(m, l, _sm) for l in range(1, m + 1))
        phi_so4 = sum(S[l - 1] @ _site(m, l, _sp) @ rt @ _site(m, l, _sm) @ S[l - 1] for l in range(1, m + 1))
        a_, b_ = _cell(idx, n)
        so4 &= np.array_equal(img, (-1) ** (bin(a_).count("1") + bin(b_).count("1")) * phi_so4)
        if np.any(img) or np.any(eta):
            sgns.add(1 if np.array_equal(img, eta) else -1 if np.array_equal(img, -eta) else 0)
        if np.any(img) or np.any(bare):
            plain.add(1 if np.array_equal(img, bare) else -1 if np.array_equal(img, -bare) else 0)
    ok9c &= len(sgns) == 1 and 0 not in sgns and so4         # and against SO4's Phi: (-1)^(p+q) on the block (p, q)
    ctl9a &= not (len(plain) == 1 and 0 not in plain)          # the bare ladder, without strings, is not it
    vac = np.zeros(2 ** n, dtype=np.int64); vac[0] = 1; w = vac.copy()
    for j in range(m + 1):
        Pj = np.diag([1 if bin(a).count("1") == j else 0 for a in range(2 ** m)]).astype(np.int64)
        pv = _to_vec(Pj @ Sig, n)
        ok9c &= np.array_equal(w, factorial(j) * pv)                 # (Lambda^dagger)^j|0> = j! P_j Sigma
        w = Lam @ w
    if J != [1] * 5:
        one, emf = _edges(n)
        ok9d &= _kept(HN, one) == 2 and _kept(HN, emf) == m + 1
# the even ring: the wrap bond joins the ket and the bra of the site 1, a second end flip
for n, J, Jw in ((4, [1, 2, 1], 3), (6, [1, 3, 2, 3, 1], 5), (8, [2, 1, 3, 5, 3, 1, 2], 4)):
    m = n // 2; wrap = _site(n, n, _sp) @ _site(n, 1, _sm); HR = _hop(n, J) + Jw * (wrap + wrap.T)
    HL = _hop(m, J[:m - 1]); Jc = J[m - 1]
    Fk = lambda r, k: _site(m, k, _sp) @ r @ _site(m, k, _sp) + _site(m, k, _sm) @ r @ _site(m, k, _sm)
    for e in np.eye(2 ** n, dtype=np.int64):
        r = _to_op(e, n)
        ok9a &= np.array_equal(_to_op(HR @ e, n), HL @ r + r @ HL + Jc * Fk(r, m) + Jw * Fk(r, 1))
# the control: couplings the reflection does not keep break the identity
HN = _hop(4, [1, 2, 3]); HL = _hop(2, [1]); spM, smM = _site(2, 2, _sp), _site(2, 2, _sm)
ctl9b = False
for i in range(16):
    r = _to_op(np.eye(16, dtype=np.int64)[i], 4)
    ctl9b |= not np.array_equal(_to_op(HN[:, i], 4), HL @ r + r @ HL + 2 * (spM @ r @ spM + smM @ r @ smM))
# zero bonds in the half chain: F4's count over the components
from math import prod
law9 = {}
for n, J in ((6, [1, 0, 2, 0, 1]), (8, [1, 0, 1, 2, 1, 0, 1]), (8, [0, 1, 1, 3, 1, 1, 0]), (8, [1, 1, 0, 2, 0, 1, 1])):
    m = n // 2; comps, cur = [], [1]
    for l in range(1, m):
        if J[l - 1] == 0:
            comps.append(cur); cur = [l + 1]
        else:
            cur.append(l + 1)
    comps.append(cur)
    one, emf = _edges(n)
    HN = _hop(n, J)
    got = (_kept(HN, one), _kept(HN, emf))
    pred = (2 * prod(len(c) + 1 for c in comps if m not in c), prod(len(c) + 1 for c in comps))
    law9[(n, tuple(J))] = got
    ctl9d &= got == pred
check("the right half is the bra, N = 4, 6, 8, mirror-symmetric integer couplings: on every basis state H_N acts as "
      "rho -> H_L rho + rho H_L + J_c flip(rho) on the half chain's cells, and on rho Sigma as ad_{H_L} + "
      "(-1)^M J_c flip, entry for entry; the even ring adds J_w times the same flip at the site 1; control: the couplings (1, 2, 3) at N = 4 break it", ok9a and ctl9b)
check("X^M and X^M Z^M commute with H_L, live on the far edge, and times Sigma are eigenstates of H_N at (-1)^M J_c and "
      "-(-1)^M J_c; the "
      "popcount projectors times Sigma are (Lambda^dagger)^j|0> / j!, j = 0..M; Lambda^dagger commutes with H_N "
      "and acts on the twisted cells as the eta ladder with Jordan-Wigner strings from the centre end, one sign "
      "throughout, which is PROOF_FROZEN_BAND_SO4's Phi times (-1)^(p + q) on the half chain's block (p, q); control: the ladder without strings does not", ok9b and ok9c and ctl9a)
check("at non-uniform mirror-symmetric couplings without zeros (N = 4, 6, 8) the largest subspace H_N keeps on the "
      "one-per-pair configurations is 2 and on the empty-or-full ones M + 1, exact rational nullspaces; with zero "
      "bonds in the half chain and the central bond nonzero, 2 prod(|c| + 1) over the components without the site M and prod(|c| + 1) over all, "
      "met exactly at four chains", ok9d and ctl9d, f"{law9}")


print("G10 the even ring: two bonds across the axis")
# the flip on one site's letters: X -> X, Y -> -Y, I and Z -> 0 (it transposes the site and forgets its quiet letters)
_Y = np.array([[0, -1j], [1j, 0]])
_fl = lambda P_: _sp @ P_ @ _sp + _sm @ P_ @ _sm
ok10a = (np.array_equal(_fl(_X), _X) and np.array_equal(_fl(_Y), -_Y) and not np.any(_fl(_I2)) and not np.any(_fl(_Z)))


def _ring(n, J, Jw):
    w = _site(n, n, _sp) @ _site(n, 1, _sm)
    return _hop(n, J) + Jw * (w + w.T)


en10 = {}
for n, J, Jw in ((4, [1, 2, 1], 3), (6, [1, 3, 2, 3, 1], 5), (8, [2, 1, 3, 5, 3, 1, 2], 4), (8, [2, 1, 3, 5, 3, 1, 2], 5),
                 (6, [1, 3, 2, 3, 1], -2), (4, [1, 1, 1], 1), (6, [1] * 5, 1), (8, [1] * 7, 1),
                 (6, [1, 3, 2, 3, 1], 2), (8, [2, 1, 3, 5, 3, 1, 2], -5)):           # the last two: equal weights, M >= 3
    m = n // 2; HN = _ring(n, J, Jw); HL = _hop(m, J[:m - 1]); Jc = J[m - 1]
    Sig = _kr([_Z if l % 2 else _I2 for l in range(1, m + 1)])
    Fk = lambda r, k: _site(m, k, _sp) @ r @ _site(m, k, _sp) + _site(m, k, _sm) @ r @ _site(m, k, _sm)
    for e in np.eye(2 ** n, dtype=np.int64):
        rt = _to_op(e, n) @ Sig
        ok10a &= np.array_equal(_to_op(HN @ e, n) @ Sig, HL @ rt - rt @ HL + (-1) ** m * Jc * Fk(rt, m) - Jw * Fk(rt, 1))
    tone = (-1) ** m * Jc - Jw
    XM = _kr([_X] * m); ZM = _kr([_Z] * m)
    for C, f in ((XM, 1), (XM @ ZM, -1)):
        v = _to_vec(C @ Sig, n)
        ok10a &= np.array_equal(HN @ v, f * tone * v)
    one, emf = _edges(n)
    ok10a &= _kept(HN, one) == 2 and _kept(HN, emf) == m + 1
    en10[(n, Jw, tuple(J))] = tone
# the exception: N = 4, J_w = -J_c, both half-chain sites carry a flip of equal weight
HN = _ring(4, [1, 2, 1], -2); one, emf = _edges(4)
v3 = _to_vec(_kr([_X, _X]) @ np.diag([1, 0, 0, -1]) @ _kr([_Z, _I2]), 4)
exc10 = _kept(HN, one) == 3 and _kept(HN, emf) == 3 and not np.any(HN @ v3) and np.any(v3)
check("the even ring, N = 4, 6, 8, mirror-symmetric integer couplings with the wrap bond J_w (negative and equal-"
      "weight cases included): the flip sends the letters X, Y, I, Z of its site to X, -Y, 0, 0; on rho Sigma the ring is "
      "ad_{H_L} + (-1)^M J_c flip_M - J_w flip_1, entry for entry; X^M and X^M Z^M times Sigma are its eigenstates at "
      "+-((-1)^M J_c - J_w), the only ones on the one-per-pair configurations but at N = 4 with J_w = -J_c, where "
      "X^2 (P_0 - P_2) Sigma is a third at 0, and the popcount projectors the only ones on the empty-or-full ones "
      "(exact); on the uniform ring the frequency is 0 at even M and 2J at odd M",
      ok10a and exc10 and en10[(4, 1, (1, 1, 1))] == 0 and en10[(8, 1, (1,) * 7)] == 0 and abs(en10[(6, 1, (1,) * 5)]) == 2,
      f"frequencies {en10}")


def _ring_sector(n, p, J, Jw):
    """gate convention (bit l = site l): the open bonds J and the wrap bond between the sites n - 1 and 0"""
    bs = subsets(n, p); ix = {b: i for i, b in enumerate(bs)}
    H = np.zeros((len(bs), len(bs)), dtype=np.int64)
    for A in bs:
        for l in range(n - 1):
            if ((A >> l) & 1) != ((A >> (l + 1)) & 1):
                H[ix[A ^ (3 << l)], ix[A]] += J[l]
        if (A & 1) != ((A >> (n - 1)) & 1):
            H[ix[A ^ (1 | 1 << (n - 1))], ix[A]] += Jw
    return H, bs


def _ring_room(n, p, q, k, J, Jw, P):
    """the largest room at -2 gamma k in (p, q) of the ring, from above (modulo P)"""
    Hp, bp = _ring_sector(n, p, J, Jw); Hq, bq = _ring_sector(n, q, J, Jw)
    K = (np.kron(Hp, np.eye(len(bq), dtype=np.int64)) - np.kron(np.eye(len(bp), dtype=np.int64), Hq)) % P
    ham = np.array([bin(A ^ B).count("1") for A in bp for B in bq]); dim = len(ham)
    piv, basis = [], np.zeros((0, dim), dtype=np.int64)
    front = [np.eye(dim, dtype=np.int64)[i] for i in np.nonzero(ham != k)[0]]
    while front:
        new = []
        for v in front:
            v = v % P
            if piv:
                v = (v - (v[piv] @ basis) % P) % P
            nz = np.nonzero(v)[0]
            if not len(nz):
                continue
            c = nz[0]; v = (v * pow(int(v[c]), P - 2, P)) % P
            if piv:
                basis = (basis - np.outer(basis[:, c], v)) % P
            basis = np.vstack([basis, v]); piv.append(c); new.append(v)
        front = [(K @ v) % P for v in new]
    return dim - len(piv)


def _ring_known(n, J, Jw, axes):
    """the ladder images of (M, 0) in (M, 2) and the mirror products of the bond axes r: their rank (modulo P1),
    whether they sit at height N/2 alone, whether the ladder images and the whole span are kept by the turning"""
    m = n // 2
    Hm, bm = _ring_sector(n, m, J, Jw); H2, b2 = _ring_sector(n, 2, J, Jw)
    K = np.kron(Hm, np.eye(len(b2), dtype=np.int64)) - np.kron(np.eye(len(bm), dtype=np.int64), H2)
    ham = np.array([bin(A ^ B).count("1") for A in bm for B in b2])
    lad = move(n, np.eye(comb(n, m), dtype=np.int64), m, 0, m, 2)
    _ex = lambda A_: DomainMatrix([[ZZ(int(x)) for x in row] for row in A_], A_.shape, ZZ).convert_to(sp.QQ).rank()
    lad_kept = _ex(np.vstack([lad.T, (K @ lad).T])) == _ex(lad.T)          # exact rational ranks
    tone = (-1) ** m * J[m - 1] - Jw
    prods = []
    for r in axes:
        mir = lambda l: (n - 1 - ((l - r) % n) + r) % n
        pairs = {tuple(sorted((l, mir(l)))) for l in range(n)}
        one = lambda s: all(((s >> a) & 1) != ((s >> b) & 1) for a, b in pairs)
        sym = lambda s: all(((s >> a) & 1) == ((s >> b) & 1) for a, b in pairs)
        def eig_on(H, bs, inset, e):
            S = [i for i, s in enumerate(bs) if inset(s)]
            Mx = sp.Matrix(H[:, S]) - e * sp.Matrix(np.eye(len(bs), dtype=np.int64)[:, S]); out = []
            for v in Mx.nullspace():
                w = v * sp.ilcm(*[x.q for x in v]); x = np.zeros(len(bs), dtype=np.int64); x[S] = [int(t) for t in w]
                out.append(x)
            return out
        psis = eig_on(Hm, bm, one, tone) + (eig_on(Hm, bm, one, -tone) if tone else [])
        for v in psis:
            for lv in eig_on(H2, b2, sym, 0):
                prods.append(np.outer(v, lv).ravel())
    Pr = np.array(prods)
    pr_ok = not np.any(Pr[:, ham != m])                     # at distance N/2 alone (kept: exact eigen-products)
    V = np.vstack([lad.T, Pr])
    return lad_kept, rank_mod(Pr, P1), (_ex(V) if n == 4 else rank_mod(V, P1)), pr_ok


r10 = {}
for n, J, Jw, axes in ((6, [1, 3, 2, 3, 1], 5, [0]), (6, [1] * 5, 1, range(6)), (4, [1, 2, 1], 3, [0]),
                       (8, [2, 1, 3, 5, 3, 1, 2], 4, [0]), (4, [1, 1, 1], 1, range(4)), (8, [1] * 7, 1, range(8))):
    m = n // 2
    room = min(_ring_room(n, m, 2, m, J, Jw, P) for P in (P1, P2))
    r10[(n, tuple(J), Jw)] = (room,) + _ring_known(n, J, Jw, axes)
odd_ok = all(r10[k][0] == r10[k][2] and not r10[k][1] and r10[k][4] for k in [(6, (1, 3, 2, 3, 1), 5), (6, (1,) * 5, 1)])
odd_ok &= r10[(6, (1, 3, 2, 3, 1), 5)][0] == 2 and r10[(6, (1,) * 5, 1)][0] == 6
even_ok = all(r10[k][1] and r10[k][4] and r10[k][0] == r10[k][3] == comb(k[0], k[0] // 2) + 2
              for k in [(4, (1, 2, 1), 3), (8, (2, 1, 3, 5, 3, 1, 2), 4)])
check("the room at k = N/2 in the block (M, 2) of the ring (from above modulo two primes, from below by exact "
      "eigen-products ranked modulo a prime; kept or not by exact ranks): at N = 6 the ladder images of (M, 0) are not kept and the room is "
      "the mirror products alone, 2 at mirror-symmetric couplings and 6 on the uniform ring, the products of its three "
      "bond axes; at N = 4 and 8 (non-uniform mirror-symmetric) the ladder images are kept and with the products "
      "fill the room, C(N, M) + 2", odd_ok and even_ok, f"{r10}")
H4, b4 = _ring_sector(4, 2, [1, 1, 1], 1)
K4 = np.kron(H4, np.eye(len(b4), dtype=np.int64)) - np.kron(np.eye(len(b4), dtype=np.int64), H4)
exact4 = _kept(K4, [i for i, (A, B) in enumerate((A, B) for A in b4 for B in b4) if bin(A ^ B).count("1") == 2])
gap_ok = (exact4, r10[(4, (1, 1, 1), 1)][0], r10[(4, (1, 1, 1), 1)][3]) == (11, 11, 9) and     (r10[(8, (1,) * 7, 1)][0], r10[(8, (1,) * 7, 1)][3]) == (80, 78) and r10[(4, (1, 1, 1), 1)][1] and r10[(8, (1,) * 7, 1)][1]
check("the uniform ring at even M: at N = 4 the room is 11 against a span of rank 9 from the ladder images and the "
      "products of all its bond axes, both exact over the rationals; at N = 8 at most 80 (modulo two primes) against "
      "a span of rank at least 78 (modulo a prime): whether more, not shown",
      gap_ok)


print("grey glass gate:", "ALL GREEN" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}")
