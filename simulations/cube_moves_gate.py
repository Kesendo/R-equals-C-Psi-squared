"""The moves on the letter cube, and the two questions each move answers (docs/THE_ONE_SQUARE.md §9).

Coordinates of a Pauli string: (k_Z, k_X, k_Y) = the number of sites whose letter anticommutes with Z, X, Y. Under
Z-dephasing at rate γ the Hermitian part of L is −2γ·k_Z, so a rung −2kγ is the plane ⟨k_Z⟩ = k.

Every move is asked two questions. (1) Does it commute with L? Then it carries a COPY: an invariant subspace, a sector,
a reading that survives. (2) Is it the one-sided multiplication ρ ↦ ρ·F by an element F that commutes with H and
anticommutes with every jump? Then it carries a REFLECTION: an element of F158's far end, and the palindrome.
(Conjugation by such an F negates the jumps and still carries a copy, the dissipator of −Z being that of Z.) The gate below checks both answers for each row of the table, exactly where an exact route
exists (integer or dyadic entries, residuals compared to 0.0) and with a control that must fail.

Conventions: Pauli book H = J·Σ_bonds (XX + YY + ZZ) with J = 1 (integer entries); site l is bit N−1−l; row-major vec,
vec(AρB) = (A ⊗ Bᵀ) vec ρ; γ = 1/2 so every entry of L is a Gaussian integer or a half-integer and products are exact.

M1  the corner shift ρ ↦ ρ·F: F = Z^N (a dark corner, [H, F] = 0) commutes with L exactly; F = X^N (a lit corner,
    [H, F] = 0) carries S·L·S⁻¹ = −L† − 2σ exactly (σ = Σγ). Each field breaks exactly one: an X field (odd k_Z, it
    anticommutes with Z^N) breaks the copy and keeps the reflection, a Z field (odd k_X) the other way round.
M2  the three turns Ad_{P^N} = (−1)^{k_P}: commute with L exactly for the Heisenberg chain (every term even k_P for every
    P; a Pauli jump goes to ± itself and its dissipator is unchanged); an X field keeps only the turn by X^N.
M3  the eigenbasis on the centre line: the shift by Z^N is, in the computational basis, the parity of the bra's popcount,
    a diagonal sign, so L is block diagonal in that parity exactly when [H, Z^N] = 0 (Heisenberg, XY, XX alone), and
    EVERY vector supported on one bra-parity block has ⟨k_X⟩ = ⟨k_Y⟩ = N/2 exactly (the shift flips k_X and k_Y,
    k ↦ N − k, with |coefficient| kept): checked in Gaussian integers on random block vectors. So each eigenspace of L
    splits by parity and has a basis on the line (measured: eigenvectors computed per block sit on it to solver
    precision), while an eigenvector mixing the two parities need not (the identity in the kernel sits at the corner;
    at odd N every eigenvalue occurs in both sectors, Ad_{X^N} exchanging them, and at N = 4 several do). Control: an
    X field makes L's parity off-block nonzero exactly.
M4  the symmetry that is a reflection: P3 Heisenberg, Z jump on the middle site, fields X on site 0 and Y on site 2 of
    equal size. No single string is a colouring (exact over all 64), yet W = SWAP₀₂·(X+Y)^⊗3, the end swap composed
    with the half-turn about the bisector of X and Y, commutes with H exactly and anticommutes with the jump exactly:
    a symmetry of H that negates the light is a far element, and the full L pairs (measured, multiset matching); at
    h₂ = 2h₀ W no longer fixes H and the spectrum does not pair. Positive control for the colouring search: with the
    fields removed the same loop finds X^⊗3 and Y^⊗3.
M5  the quarter-turn about Z, S = diag(1, i) per site: commutes with L exactly for XXZ at Δ = 1 and Δ = 1/2; an X field
    breaks it.
M6  the letter swap that counts F85: among the Pauli strings of k letters, truly (#Y even and #Z even) = even non-truly
    (k_X even, not truly) + 1 at every k = 1..7, the extra one being X^k; the pairing is Y ↔ Z in the first non-X letter.
M7  the plane: at rung k the blocks (p, q) with p − q ≢ k (mod 2) whose Hamming distances straddle k hold cells on both
    sides of the plane ⟨k_Z⟩ = k and none on it, so det(B + 2kγ) has leading coefficient Π(2k − 2·Hamming) ≠ 0 and
    finitely many positive roots E_k; the blocks with p − q ≡ k may hold cells on the plane and are not fenced.
    Heisenberg chain N = 3, 4, k = 2, 3: leading coefficients exact; at γ/J = 1.3 no real eigenvalue at −2kγ in the
    wrong-parity blocks (0 is the fence) and some in the right-parity blocks (measured, not fenced); exactly one real
    eigenvalue per root; E_{N−k} = E_k block by block (the F1 mirror sends (p, q) to (p, N − q)) and
    E_2(4, chain) = {0.947023, 1.059767, √2} from the blocks (1, 2) and (2, 3), minimal polynomials asserted.
"""
import sys
from itertools import product
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ_I

sys.stdout.reconfigure(encoding="utf-8")
FAIL = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], dtype=complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
LET = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def kron_all(ms):
    r = np.eye(1, dtype=complex)
    for m in ms: r = np.kron(r, m)
    return r


def op(N, l, P):
    return kron_all([P if j == l else I2 for j in range(N)])


def string(s):
    return kron_all([LET[c] for c in s])


def heis(N, bonds, delta=1.0, hx=0.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for a, b in bonds:
        H += op(N, a, X) @ op(N, b, X) + op(N, a, Y) @ op(N, b, Y) + delta * op(N, a, Z) @ op(N, b, Z)
    for l in range(N): H += hx * op(N, l, X)
    return H


def liou(N, H, gam, jumps=None):
    """L = −i(H⊗I − I⊗Hᵀ) + Σ_l γ(Z_l ⊗ Z_lᵀ − I) over the jump sites (default: every site)."""
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in (range(N) if jumps is None else jumps):
        Zl = op(N, l, Z); L += gam * (np.kron(Zl, Zl.T) - np.eye(d * d))
    return L


def right_mult(N, F):
    return np.kron(np.eye(2 ** N), F.T)


def conj_super(N, U):
    return np.kron(U, U.conj())


def exact_zero(A):
    return float(np.max(np.abs(A))) == 0.0


GAM = 0.5
chain = lambda N: [(i, i + 1) for i in range(N - 1)]

# ---------------------------------------------------------------- M1 the corner shift
N = 3
H = heis(N, chain(N)); L = liou(N, H, GAM); sigma = N * GAM
SZ = right_mult(N, string('ZZZ')); SX = right_mult(N, string('XXX'))
check("M1 dark corner Z^N: ρ ↦ ρ·Z^N commutes with L exactly (P3 Heisenberg, Z jumps, γ = 1/2)", exact_zero(L @ SZ - SZ @ L))
res = SX @ L @ np.linalg.inv(SX) + L.conj().T + 2 * sigma * np.eye(L.shape[0])
check("M1 lit corner X^N: S·L·S⁻¹ = −L† − 2σ exactly (the one-sided factor of Π_Z, a far element)", exact_zero(res))
Hx = heis(N, chain(N), hx=1.0); Lx = liou(N, Hx, GAM)
Hz = heis(N, chain(N)) + sum(op(N, l, Z) for l in range(N)); Lz = liou(N, Hz, GAM)
resx = SX @ Lx @ np.linalg.inv(SX) + Lx.conj().T + 2 * sigma * np.eye(L.shape[0])
resz = SX @ Lz @ np.linalg.inv(SX) + Lz.conj().T + 2 * sigma * np.eye(L.shape[0])
Hy = heis(N, chain(N)) + sum(op(N, l, Y) for l in range(N)); Ly = liou(N, Hy, GAM)
resy = SX @ Ly @ np.linalg.inv(SX) + Ly.conj().T + 2 * sigma * np.eye(L.shape[0])
check("M1 control: an X field (odd k_Z) breaks the dark copy Z^N and keeps the lit reflection X^N; a Z field (odd k_X) keeps the copy and breaks the reflection; a Y field (odd in both) breaks both",
      not exact_zero(Lx @ SZ - SZ @ Lx) and exact_zero(resx) and exact_zero(Lz @ SZ - SZ @ Lz) and not exact_zero(resz)
      and not exact_zero(Ly @ SZ - SZ @ Ly) and not exact_zero(resy))

# ---------------------------------------------------------------- M2 the three turns
for P, name in ((X, 'X'), (Y, 'Y'), (Z, 'Z')):
    T = conj_super(N, string(name * N))
    check(f"M2 the turn by {name}^N commutes with L exactly (Heisenberg: every term has even k_{name})", exact_zero(L @ T - T @ L))
    ok_field = exact_zero(Lx @ T - T @ Lx)
    check(f"M2 control: under an X field the turn by {name}^N {'survives' if name == 'X' else 'breaks'}", ok_field == (name == 'X'))

# ---------------------------------------------------------------- M3 the centre line
def pauli_coeffs_times_2N(N, rho):
    """2^N · c_s = Tr(σ_s ρ) for every string s, exact for Gaussian-integer ρ (entries of σ_s are 0, ±1, ±i)."""
    out = {}
    for s in product('IXYZ', repeat=N):
        out[''.join(s)] = np.trace(string(''.join(s)) @ rho)
    return out


def kcoords(s):
    kZ = sum(c in 'XY' for c in s); kX = sum(c in 'YZ' for c in s); kY = sum(c in 'XZ' for c in s)
    return kZ, kX, kY


rng = np.random.default_rng(7)
for N in (3, 4):
    d = 2 ** N
    for label, Hc in (("Heisenberg", heis(N, chain(N))), ("XY", heis(N, chain(N), delta=0.0)),
                      ("XX alone", sum(op(N, a, X) @ op(N, b, X) for a, b in chain(N)))):
        Lc = liou(N, Hc, GAM); S = right_mult(N, string('Z' * N))
        check(f"M3 N={N} {label}: [L, ρ ↦ ρ·Z^N] = 0 exactly ([H, Z^N] = 0: the shift is the bra parity)", exact_zero(Lc @ S - S @ Lc))
    # every vector supported on one bra-parity block: <k_X> = <k_Y> = N/2 exactly, in Gaussian integers
    ok = True
    for parity in (0, 1):
        rho = np.zeros((d, d), dtype=complex)
        for x in range(d):
            for y in range(d):
                if bin(y).count('1') % 2 == parity:
                    rho[x, y] = complex(rng.integers(-5, 6), rng.integers(-5, 6))
        c = pauli_coeffs_times_2N(N, rho)
        m2 = lambda v: v.real ** 2 + v.imag ** 2          # exact for Gaussian integers (abs() goes through a root)
        sx = sum(m2(v) * (2 * kcoords(s)[1] - N) for s, v in c.items())
        sy = sum(m2(v) * (2 * kcoords(s)[2] - N) for s, v in c.items())
        tot = sum(m2(v) for v in c.values())
        ok &= sx == 0 and sy == 0 and tot > 0
    check(f"M3 N={N}: a random Gaussian-integer vector on one bra-parity block has Σ|c_s|²(2k_X − N) = 0 and the same for k_Y, exactly", ok)
    # the block structure, exactly: order the flat indices by bra parity; L's off-block must be 0.0
    even = [x * d + y for x in range(d) for y in range(d) if bin(y).count('1') % 2 == 0]
    odd = [x * d + y for x in range(d) for y in range(d) if bin(y).count('1') % 2 == 1]
    Lh = liou(N, heis(N, chain(N)), GAM)
    check(f"M3 N={N}: L is block diagonal in the bra parity exactly (Heisenberg; the off-block compared to 0.0)",
          exact_zero(Lh[np.ix_(even, odd)]) and exact_zero(Lh[np.ix_(odd, even)]))

    def line_deviation(vecs):
        out = 0.0
        for i in range(vecs.shape[1]):
            c = pauli_coeffs_times_2N(N, vecs[:, i].reshape(d, d)); tot = sum(abs(v) ** 2 for v in c.values())
            out = max(out, abs(sum(abs(v) ** 2 * kcoords(s)[1] for s, v in c.items()) / tot - N / 2))
        return out

    # measured beside it: eigenvectors computed per parity block sit on the line to solver precision, the full
    # solver's eigenvectors, free to mix degenerate sectors, do not have to
    per_block = []
    for idx in (even, odd):
        w, V = np.linalg.eig(Lh[np.ix_(idx, idx)])
        full = np.zeros((d * d, V.shape[1]), dtype=complex); full[idx, :] = V
        per_block.append(line_deviation(full))
    w, V = np.linalg.eig(Lh)
    mixed = line_deviation(V)
    print(f"      measured N={N}: eigenvectors per parity block |<k_X> − N/2| ≤ {max(per_block):.1e}; the full solver's, free to mix "
          f"degenerate sectors, up to {mixed:.2f}")
    check(f"M3 N={N}: the per-block eigenvectors sit on the line to solver precision (≤ 1e-10) while the mixed ones leave it by more than 0.1",
          max(per_block) < 1e-10 and mixed > 0.1)
    # control: an X field breaks the block structure exactly
    Lxf = liou(N, heis(N, chain(N), hx=1.0), GAM)
    check(f"M3 N={N} control: an X field makes the parity off-block of L nonzero (exactly)",
          not exact_zero(Lxf[np.ix_(even, odd)]))

# ---------------------------------------------------------------- M4 the symmetry that is a reflection
N = 3
H = heis(N, chain(N))
SWAP02 = np.zeros((8, 8))
for x in range(8):
    b = [(x >> 2) & 1, (x >> 1) & 1, x & 1]; y = (b[2] << 2) | (b[1] << 1) | b[0]
    SWAP02[y, x] = 1
# Heisenberg P3, Z jump on the middle site, fields h₀·X on site 0 and h₂·Y on site 2 (two different letters, so no
# single string is a colouring); the carrier is the half-turn about the bisector n = (X+Y)/√2 composed with the end swap,
# which sends the X field to the Y field exactly when h₀ = h₂. W = SWAP₀₂·(X+Y)^⊗3 is that carrier up to the scalar
# 2^(3/2), kept unnormalised so every entry is an integer.
Z1 = op(N, 1, Z)
W = SWAP02 @ string('XXX') + SWAP02 @ string('YYY') + SWAP02 @ (string('XXY') + string('XYX') + string('YXX') + string('XYY') + string('YXY') + string('YYX'))
Hs = heis(N, chain(N)) + op(N, 0, X) + op(N, 2, Y)          # h₀ = h₂ = 1
Hu = heis(N, chain(N)) + op(N, 0, X) + 2 * op(N, 2, Y)      # h₀ = 1, h₂ = 2
no_colouring = all(not (exact_zero(string(s) @ Hs - Hs @ string(s)) and exact_zero(string(s) @ Z1 + Z1 @ string(s))) for s in (''.join(t) for t in product('IXYZ', repeat=3)))
colourings_free = [s for s in (''.join(t) for t in product('IXYZ', repeat=3)) if exact_zero(string(s) @ heis(N, chain(N)) - heis(N, chain(N)) @ string(s)) and exact_zero(string(s) @ Z1 + Z1 @ string(s))]
check("M4 fields X on site 0 and Y on site 2 (equal), Z jump on site 1: no single Pauli string commutes with H and anticommutes with the jump (no colouring, exact over all 64); with the fields removed the same loop finds exactly XXX and YYY",
      no_colouring and sorted(colourings_free) == ['XXX', 'YYY'])
check("M4 the symmetry that is a reflection: W = SWAP₀₂·(X+Y)^⊗3 commutes with H exactly and anticommutes with the jump exactly (a far element)",
      exact_zero(W @ Hs - Hs @ W) and exact_zero(W @ Z1 + Z1 @ W))


def pairing_distance(L, sigma):
    """Assignment matching: every eigenvalue is paired with a distinct mirrored one by the minimum-sum assignment,
    and the largest matched distance is returned (an upper bound on the bottleneck matching distance)."""
    from scipy.optimize import linear_sum_assignment
    ev = np.linalg.eigvals(L)
    mirrored = -2 * sigma - ev
    cost = np.abs(ev[:, None] - mirrored[None, :])
    rows, cols = linear_sum_assignment(cost)
    return float(np.max(cost[rows, cols]))


ds, du = pairing_distance(liou(N, Hs, GAM, jumps=[1]), GAM), pairing_distance(liou(N, Hu, GAM, jumps=[1]), GAM)
check(f"M4 the spectrum pairs about −σ at equal fields (assignment pairing distance {ds:.1e}, measured) and not at h₂ = 2h₀ ({du:.3f}), where W no longer fixes H (nonzero exactly)",
      ds < 1e-10 and du > 1e-3 and not exact_zero(W @ Hu - Hu @ W))

# ---------------------------------------------------------------- M5 the quarter-turn
Sq = np.diag([1, 1j])
for N in (3, 4):
    Q = conj_super(N, kron_all([Sq] * N))
    for delta in (1.0, 0.5):
        Ld = liou(N, heis(N, chain(N), delta=delta), GAM)
        check(f"M5 N={N} XXZ Δ={delta}: the quarter-turn about Z commutes with L exactly", exact_zero(Ld @ Q - Q @ Ld))
    Lxf = liou(N, heis(N, chain(N), hx=1.0), GAM)
    check(f"M5 N={N} control: an X field breaks the quarter-turn", not exact_zero(Lxf @ Q - Q @ Lxf))

# ---------------------------------------------------------------- M6 the letter swap that counts F85
def partner(s):
    """Y ↔ Z in the first non-X letter."""
    for i, c in enumerate(s):
        if c != 'X':
            return s[:i] + ('Z' if c == 'Y' else 'Y') + s[i + 1:]
    return None


ok = True
rows = []
for k in range(1, 8):
    strings = [''.join(t) for t in product('XYZ', repeat=k)]
    truly = [s for s in strings if s.count('Y') % 2 == 0 and s.count('Z') % 2 == 0]
    even_nt = [s for s in strings if (s.count('Y') + s.count('Z')) % 2 == 0 and s not in truly]
    closed_t = (3 ** k + 2 + (-1) ** k) // 4; closed_e = (3 ** k - 2 + (-1) ** k) // 4
    # the pairing: partner is a bijection truly \ {X^k} -> even non-truly
    img = {partner(s) for s in truly if s != 'X' * k}
    ok &= len(truly) == closed_t and len(even_nt) == closed_e and len(truly) == len(even_nt) + 1 and img == set(even_nt) \
        and all(partner(partner(s)) == s for s in truly if s != 'X' * k)
    rows.append(f"k={k}: {len(truly)}/{len(even_nt)}")
check("M6 F85's strings: truly = (3^k+2+(−1)^k)/4, even non-truly = (3^k−2+(−1)^k)/4, truly = even non-truly + 1 with X^k the extra one and Y↔Z in the first non-X letter the pairing, k = 1..7",
      ok, "; ".join(rows))

# ---------------------------------------------------------------- M7 the plane
g = sp.symbols('g')
R = ZZ_I[g]


def ham_matrix_int(N, bonds):
    d = 2 ** N
    Hm = [[0] * d for _ in range(d)]
    bit = lambda x, l: (x >> (N - 1 - l)) & 1
    for a, b in bonds:
        for x in range(d):
            Hm[x][x] -= 1
            y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b)) if bit(x, a) != bit(x, b) else x
            Hm[y][x] += 2
    return Hm


def block_idx(N, p, q):
    d = 2 ** N
    pc = lambda x: bin(x).count('1')
    idx = [(x, y) for x in range(d) for y in range(d) if pc(x) == p and pc(y) == q]
    return idx, {e: i for i, e in enumerate(idx)}


def block_numeric(N, Hm, p, q, gam):
    d = 2 ** N
    idx, pos = block_idx(N, p, q)
    n = len(idx)
    B = np.zeros((n, n), dtype=complex)
    for c, (x, y) in enumerate(idx):
        for xp in range(d):
            if Hm[xp][x]: B[pos[(xp, y)], c] += -1j * Hm[xp][x]
        for yp in range(d):
            if Hm[y][yp]: B[pos[(x, yp)], c] += 1j * Hm[y][yp]
        B[c, c] += -2 * gam * bin(x ^ y).count('1')
    return B


def shifted_poly(N, Hm, p, q, k):
    d = 2 ** N
    idx, pos = block_idx(N, p, q)
    n = len(idx)
    rows_ = [[R.zero] * n for _ in range(n)]
    I = R.convert(sp.I); G = R.convert(g)
    hams = []
    for c, (x, y) in enumerate(idx):
        for xp in range(d):
            if Hm[xp][x]: rows_[pos[(xp, y)]][c] += -I * R.convert(Hm[xp][x])
        for yp in range(d):
            if Hm[y][yp]: rows_[pos[(x, yp)]][c] += I * R.convert(Hm[y][yp])
        h = bin(x ^ y).count('1'); hams.append(h)
        rows_[c][c] += G * R.convert(2 * k - 2 * h)
    return n, sp.expand(R.to_sympy(DomainMatrix(rows_, (n, n), R).det())), hams


def positive_roots(P):
    cs = sp.Poly(P, g).all_coeffs()
    A = sp.Poly([sp.re(c) for c in cs], g); B = sp.Poly([sp.im(c) for c in cs], g)
    Q = A if B.is_zero else (B if A.is_zero else sp.gcd(A, B))
    if Q.degree() <= 0: return []
    return sorted({r for r in Q.real_roots() if r > 0}, key=float)


x = sp.symbols('x')
E = {}
EBLOCK = {}
for N in (3, 4):
    Hm = ham_matrix_int(N, chain(N))
    for k in (2, 3):
        leads_ok = True; wrong_generic = 0; plane_generic = 0; roots_ok = True
        Ek = {}
        for p in range(N + 1):
            for q in range(N + 1):
                if abs(p - q) > k: continue
                idx, _ = block_idx(N, p, q)
                hams = {bin(xx ^ yy).count('1') for xx, yy in idx}
                B = block_numeric(N, Hm, p, q, 1.3)
                c = int(np.sum(np.abs(np.linalg.eigvals(B) + 2 * k * 1.3) < 1e-6))
                if (p - q - k) % 2 == 1:
                    wrong_generic += c
                    if q >= p and max(hams) > k and min(hams) < k:
                        n, P, hs = shifted_poly(N, Hm, p, q, k)
                        leads_ok &= sp.simplify(sp.Poly(P, g).LC() - sp.prod([2 * k - 2 * h for h in hs])) == 0
                        rts = positive_roots(P)
                        Ek[(p, q)] = rts
                        for r in rts:
                            Br = block_numeric(N, Hm, p, q, float(r))
                            roots_ok &= int(np.sum(np.abs(np.linalg.eigvals(Br) + 2 * k * float(r)) < 1e-6)) == 1
                else:
                    plane_generic += c
        E[(N, k)] = sorted({r for rts in Ek.values() for r in rts}, key=float)
        EBLOCK[(N, k)] = Ek
        check(f"M7 N={N} k={k}: straddling blocks {sorted(Ek)} have leading coefficient Π(2k − 2·Hamming) exactly; at γ/J = 1.3 the wrong-parity blocks hold {wrong_generic} real eigenvalues at −{2*k}γ (0 expected) and the right-parity blocks {plane_generic} (measured, not fenced); each root carries exactly one",
              leads_ok and wrong_generic == 0 and roots_ok,
              "E_k = {" + ", ".join(f"{float(r):.6f}" for r in E[(N, k)]) + "}")
# the F1 mirror: E_{N-k} = E_k; the self-mirror rung at N = 4 pinned
E1_3 = [sp.sqrt((sp.sqrt(17) - 1) / 2), sp.sqrt(3)]
check("M7 N=3: E_2 = E_1 = {√((√17−1)/2), √3} (rung 2 = rung N−1, the F1 mirror)",
      len(E[(3, 2)]) == 2 and all(sp.simplify(a - b) == 0 for a, b in zip(E[(3, 2)], E1_3)))
E1_4 = [0.745022, 0.745439, 0.910180, 1.545936, 1.572303, 1.867978, 2.088800, 2.197368]
check("M7 N=4: E_3 = E_1 (the eight points of f50_exceptional_couplings.py, 6 decimals)",
      len(E[(4, 3)]) == 8 and all(abs(float(a) - b) < 5e-7 for a, b in zip(E[(4, 3)], E1_4)))
# the F1 mirror block by block: right multiplication by X^N sends (p, q) to (p, N − q), so at N = 4 the rung-3 set of
# (1,3) is the rung-1 set of (1,1), and (2,2) is its own mirror; the rung-1 diagonal sets are recomputed here, not pinned
Hm4 = ham_matrix_int(4, chain(4))
def diag_roots(p):
    _, P, _ = shifted_poly(4, Hm4, p, p, 1)
    return positive_roots(P)
same = lambda a, b: len(a) == len(b) and all(sp.simplify(u - v) == 0 for u, v in zip(a, b))
check("M7 N=4: the F1 mirror block by block: E_3 of (1,3) = E_1 of (1,1), and E_3 of (2,2) = E_1 of (2,2), exactly",
      same(EBLOCK[(4, 3)][(1, 3)], diag_roots(1)) and same(EBLOCK[(4, 3)][(2, 2)], diag_roots(2)))
e2 = E[(4, 2)]
MINPOLY = [x**10 + 2*x**8 - 12*x**6 + 88*x**4 - 64, x**4 + 6*x**2 - 8, x**2 - 2]
check("M7 N=4: E_2 (the self-mirror rung) = {0.947023, 1.059767, √2} from the blocks (1,2) and (2,3), minimal polynomials x¹⁰+2x⁸−12x⁶+88x⁴−64, x⁴+6x²−8, x²−2 asserted",
      len(e2) == 3 and abs(float(e2[0]) - 0.947023) < 5e-7 and abs(float(e2[1]) - 1.059767) < 5e-7
      and all(sp.Poly(sp.minimal_polynomial(r, x), x) == sp.Poly(m, x) for r, m in zip(e2, MINPOLY)))

print()
print("FAILED: " + ", ".join(FAIL) if FAIL else "ALL GATES PASS")
sys.exit(1 if FAIL else 0)
