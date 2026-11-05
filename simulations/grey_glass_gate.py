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

print("grey glass gate:", "ALL GREEN" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}")
