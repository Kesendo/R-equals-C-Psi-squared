"""The moves on the letter cube, and the two questions each move answers (docs/THE_ONE_SQUARE.md §9).

Coordinates of a Pauli string: (k_Z, k_X, k_Y) = the number of sites whose letter anticommutes with Z, X, Y. Under
Z-dephasing at rate γ the Hermitian part of L is −2γ·k_Z, so a rung −2kγ is the plane ⟨k_Z⟩ = k.

Every move is asked two questions. (1) Does it commute with L? Then it carries a COPY: an invariant subspace, a sector,
a reading that survives. (2) Is it the one-sided multiplication ρ ↦ ρ·F by an element F that commutes with H and
anticommutes with every jump (Hermitian jumps that square to 1)? Then it carries a REFLECTION: an element of F158's far end, and the palindrome.
(Conjugation by such an F negates the jumps and still carries a copy, the dissipator of −Z being that of Z.) The gate below checks both answers for each row of the table, exactly where an exact route
exists (integer or dyadic entries, residuals compared to 0.0) and with a control that must fail.

Conventions (M1 to M10; M11 states its own bonds, rates and fields and works over ℚ(i)): Pauli book H = J·Σ_bonds (XX + YY + ZZ) with J = 1 (integer entries); site l is bit N−1−l; row-major vec,
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
M8  amplitude damping, a height and an arrow: σ⁻ at one site is the two transverse lights (¼ along X and along Y) plus
    one move I → Z, the only nonzero Pauli-transfer entry of the rest. At N = 3 the moves that keep the arrow, I and Z
    fixed with their signs (the turn by Z^N, the quarter-turn about Z, the transpose), commute with the damping
    exactly; the turns by X^N and Y^N, the shift by Z^N, R and the letter swap X ↔ Z do not, and none of the eight reflects it.
    The turn by Z^N and the quarter-turn still commute with the full L (Heisenberg, Z dephasing and damping). At
    infinite temperature, loss and gain at equal rates, the arrows cancel and the lights remain exactly: the turns by
    X^N and Y^N carry copies again and the shift by Z^N reflects the full Heisenberg L exactly, R no longer.
M9  the arrow reversed: S(ρ) = X^N·ρ·Y^N (per site I ↦ iZ, Z ↦ −iI, X ↔ Y) reflects the thermal damping,
    S·D·S⁻¹ = −D† − (γ↓ + γ↑), symbolic in the two rates at one site, while ρ·Z alone does not; at N = 3 with unequal
    site rates it reflects the full L exactly for the Heisenberg chain, for XX bonds alone and at γ↑ = 0, centre
    −Σ(γ↓ + γ↑)/2; co-axial Z dephasing, a longitudinal Z field and DM bonds each break it.
M10 F137's composition decided: S_XY reflects one site's X or Y dephasing (symbolic, with and without heating; Z
    dephasing is the control that fails), so damping plus transverse
    dephasing (X, Y, or X, Y, X on the sites) is reflected exactly at N = 3 with c = Σγ_T1 + 2Σγ_deph; co-axial Z, a
    transverse field and a longitudinal field each give an exactly nonzero Tr((L + c/2)^3) (exact rationals), which
    rules out every reflecting element of these generators as written (the uniform longitudinal field is the lab
    frame of the field-free chain, M11); the trace is zero where S_XY reflects.
M11 outside S_XY's class: loss and gain exchanged. One site, symbolic: a turn about Z, g = diag(1, e^{iθ}), commutes
    with the thermal damping at any angle (σ± only take a phase) and with X dephasing exactly where sin θ = 0 (the
    commutator is sin θ times a fixed nonzero matrix); S_XY reflects dephasing along
    a·X + b·Y exactly where ab = 0. For H commuting with Z^N and jumps odd under Z^N with A†A + AA† ∝ 1 (σ±,
    dephasing along any axis in the XY plane), the shift S_IZ: ρ ↦ ρ·Z^N carries L onto −L′† − c, L′ with every A replaced
    by A†, cooling and heating exchanged (exact at N = 3 with a hopping phase, Z fields, a term XX − YY that does
    not conserve the excitation number, ZZ, and dephasing along three different axes; L itself is the control),
    so the palindrome is the spectrum unchanged by exchanging loss and gain; S_XY is the spin flip after it, and
    with dephasing along X or Y (the flip mirrors a tilted axis) S_XY·L(H)·S_XY⁻¹ = −L(X^N·H·X^N)† − c and
    Ad_{X^N}·L(X^N·H·X^N)·Ad_{X^N} = L(H) exchanged (exact, the fifth power sum changing sign), and at γ↑ = γ↓
    S_IZ reflects L itself (unequal rates the control); the decay rates of the triangle at flux π/2 miss their
    mirror by 2.53·10⁻² with numpy in double precision and with mpmath at 30 digits alike, at flux 0 only by each
    one's rounding (printed as a reading, no pass or fail);
    transverse DM YZ − ZY on the three-site chain, not commuting with Z^N, breaks under cooling alone (240) and at
    equal rates (1280). For bonds
    J(XX + YY) + D_z(XY − YX) + Δ·ZZ, X^N·H·X^N = H̄, every hopping phase reversed (exact on a triangle at flux
    arccos(3/5); L(H) itself in place of L(X^N·H·X^N) is the control). Carriers, each exact with a control: turns
    about Z on the chain (N = 3, Pythagorean J, D_z per bond, ZZ of both signs; at J = 0 the Sextant's DM chain
    under loss, the element X⊗Y⊗X·ρ·Y⊗X⊗Y); the site reflection of the triangle at rates symmetric under it
    (fluxes π/2 and arccos(3/5)); the spin flip with the shift by two sites on a 4-ring at ratios 3, 2, ⅓, ½;
    the general form V·(X^N·H·X^N)·V† = H with fields, the Heisenberg chain with fields (−½, 0, ½) under its
    reflection at rates the reflection keeps, with heating and under cooling alone;
    infinite temperature with ρ·Z^N;
    the bipartite ring at Δ = 0 with dephasing along two tilted axes and X; one ratio, γ↑ = γ↓/4 on every site, by the thermal transpose
    L(H)ᵀ = Γ⁻¹·L(Hᵀ)·Γ, Hᵀ = H̄ (F119 §7 at a complex H), palindromic characteristic polynomials at flux π/2 on
    two triangles and at ratio 4 (three ratios give a nonzero fifth power sum on the same H), and one ratio per
    component of the hopping bonds (a pendant joined by ZZ alone at its own ratio pairs, joined by hopping it does
    not); X dephasing beside the one-ratio damping breaks it at flux π/2 (a nonzero seventh power sum) and not
    with real hopping, and sees the gauge (flux 0 written with phases, a nonzero ninth power sum); a bipartite ring
    at Δ = 0 with fields and four ratios (N = 4), by complex conjugation and Z on one sublattice (Δ = 7⁄10 with the
    fields, and without them, the controls). On the triangle with three ratios Tr((L + c/2)^5) = −656·sin Φ at
    twelve points e^{iΦ} = (1 + it)²/(1 + t²), t rational, eleven of which fix a trigonometric polynomial of
    degree 5, so the power sum is −656·sin Φ at every flux and nothing reflects the spectrum off Φ = 0, π (+656
    with cooling and heating exchanged; at flux π S_XY itself reflects, at π/2 not); with cooling alone the odd
    power sums to 7 vanish at flux π/2; site-dependent fields at flux 0 under cooling alone give −2560 in that frame (the
    decay rates, which a common frame shift does not move, miss their mirror by 0.109 in both precisions, read);
    in the lab frame, M10's uniform field on the Heisenberg chain under cooling alone, L is block-diagonal in the
    coherence order and each order keeps its own mirror about −c/2 + iκ, κ growing with the order (exact,
    centred characteristic polynomials of the blocks; the triangle with those fields the control). The detuned triangle pairs in no common frame (the fifth power sum
    linear in the offset, zero once, where the seventh is not). Beside cooling on the Heisenberg chain one orthogonal pair of tilted
    axes pairs and a tilted axis beside X does not (exact). With longitudinal fields the
    three chained at one ratio (S_XY, the transpose, the conjugation) and turns with W·(−H̄)·W† = H carry the pairing,
    at Δ = 0 where e^{2iΦ} = (−1)^ℓ on every cycle: the triangle at flux ±π/2 pairs with fields (−½, 0, ½) under
    cooling alone and with fields (3⁄10, −⅕, ½) at ratio ¼ (exact characteristic polynomials; Δ = 7⁄10, flux
    arccos(3/5) and three ratios the controls). Contrast: under uniform Z dephasing the same triangle at flux π/2 is not
    palindromic and at flux π it is (exact characteristic polynomials), the one-sided ρ·X^N carrying L(H) onto a
    reflection only where X^N·H·X^N = H.
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

# ---------------------------------------------------------------- M8 amplitude damping: a height and an arrow
SM = np.array([[0, 1], [0, 0]], dtype=complex)        # sigma^- (|1> -> |0>), integer entries
def diss(N, Lk, rate):
    d = 2 ** N; Id = np.eye(d)
    LdL = Lk.conj().T @ Lk
    return rate * (np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T)))
# one site: the damping is the two transverse lights plus one move I -> Z (Pauli transfer: one entry)
D1 = diss(1, SM, 1.0)
LIGHTS1 = diss(1, X, 0.25) + diss(1, Y, 0.25)
A1 = D1 - LIGHTS1
PT = np.array([[np.trace(LET[a].conj().T @ (A1 @ LET[b].reshape(-1)).reshape(2, 2)) / 2 for b in 'IXYZ'] for a in 'IXYZ'])
nz = [(a, b) for i, a in enumerate('IXYZ') for j, b in enumerate('IXYZ') if PT[i, j] != 0]
check("M8 one site: σ⁻ at rate 1 = the lights ¼ along X and along Y (rates ½, ½, 1 on X, Y, Z) plus one arrow I → Z, the only nonzero Pauli-transfer entry of the rest",
      nz == [('Z', 'I')] and exact_zero(np.diag([0, 0.5, 0.5, 1.0]) + np.array([[np.trace(LET[a].conj().T @ (LIGHTS1 @ LET[b].reshape(-1)).reshape(2, 2)) / 2 for b in 'IXYZ'] for a in 'IXYZ'])),
      f"arrow entry {PT[3, 0].real:+.1f}")
N = 3
DT = sum(diss(N, op(N, l, SM), GAM) for l in range(N))
LT = liou(N, heis(N, chain(N)), GAM) + DT
Hd = X + Z                                             # the letter swap X <-> Z, unnormalised (integer entries)
TRANSPOSE = np.zeros((4 ** N, 4 ** N))
for a in range(2 ** N):
    for b in range(2 ** N):
        TRANSPOSE[b * 2 ** N + a, a * 2 ** N + b] = 1
MOVES = {
    "the turn by Z^N": conj_super(N, string('ZZZ')),
    "the quarter-turn about Z": conj_super(N, kron_all([np.diag([1, 1j])] * N)),
    "the transpose": TRANSPOSE,
    "the turn by X^N": conj_super(N, string('XXX')),
    "the turn by Y^N": conj_super(N, string('YYY')),
    "the shift by Z^N": right_mult(N, string('ZZZ')),
    "R, the shift by X^N": right_mult(N, string('XXX')),
    "the letter swap X ↔ Z": conj_super(N, kron_all([Hd] * N)),
}
KEEP = {"the turn by Z^N", "the quarter-turn about Z", "the transpose"}
got = {name: exact_zero(M @ DT - DT @ M) for name, M in MOVES.items()}
check("M8 N=3, damping σ⁻ at rate ½ on every site: the moves that keep the arrow (I and Z fixed, sign included) commute with the damping exactly, "
      "the turns by X^N and Y^N (Z ↦ −Z, loss ↦ gain), the shift by Z^N (the arrow's ends exchanged), R and the swap X ↔ Z (ends moved off the edge) do not",
      all(got[n] == (n in KEEP) for n in MOVES), "; ".join(f"{n}: {'keeps' if g else 'breaks'}" for n, g in got.items()))
CT = -2 * np.trace(DT).real / DT.shape[0]
refl = {name: exact_zero(M @ DT @ np.linalg.inv(M) + DT.conj().T + CT * np.eye(DT.shape[0])) for name, M in MOVES.items()}
check("M8 none of the eight gated moves carries a reflection of the damping (S·D·S⁻¹ + D† + c·1 nonzero for every move, c fixed by the trace)",
      not any(refl.values()))
check("M8 on the full L (Heisenberg chain, Z dephasing and damping at ½): the turn by Z^N and the quarter-turn about Z still carry copies exactly",
      exact_zero(LT @ MOVES["the turn by Z^N"] - MOVES["the turn by Z^N"] @ LT)
      and exact_zero(LT @ MOVES["the quarter-turn about Z"] - MOVES["the quarter-turn about Z"] @ LT))
# infinite temperature: loss and gain at equal rates, the two arrows cancel and the lights remain
DI = sum(diss(N, op(N, l, SM), GAM) + diss(N, op(N, l, SM.conj().T), GAM) for l in range(N))
LIGHTS = sum(diss(N, op(N, l, X), GAM / 2) + diss(N, op(N, l, Y), GAM / 2) for l in range(N))
check("M8 infinite temperature: loss and gain at ½ each are exactly the lights at ¼ along X and along Y (the arrows cancel)", exact_zero(DI - LIGHTS))
gotI = {name: exact_zero(M @ DI - DI @ M) for name, M in MOVES.items()}
check("M8 infinite temperature: the turns by X^N and Y^N carry copies again, with the turn by Z^N, the quarter-turn and the transpose",
      all(gotI[n] for n in ("the turn by Z^N", "the turn by X^N", "the turn by Y^N", "the quarter-turn about Z", "the transpose")))
LI = -1j * (np.kron(heis(N, chain(N)), np.eye(2 ** N)) - np.kron(np.eye(2 ** N), heis(N, chain(N)).T)) + DI
CI = -2 * np.trace(LI).real / LI.shape[0]
SZ3 = MOVES["the shift by Z^N"]; SX3 = MOVES["R, the shift by X^N"]
check("M8 infinite temperature on the Heisenberg chain: the shift by Z^N, the half-turn about the k_Z axis, reflects the full L exactly "
      "(S·L·S⁻¹ = −L† − c, c = 2σ of the lights), and R no longer does",
      exact_zero(SZ3 @ LI @ np.linalg.inv(SZ3) + LI.conj().T + CI * np.eye(LI.shape[0]))
      and not exact_zero(SX3 @ LI @ np.linalg.inv(SX3) + LI.conj().T + CI * np.eye(LI.shape[0])), f"c = {CI}")

# ---------------------------------------------------------------- M9 the arrow reversed: the thermal reflection
# S(ρ) = X^N·ρ·Y^N, the turn by X^N followed by the shift by Z^N up to a phase: per site I ↦ iZ, Z ↦ −iI, X ↔ Y.
gd, gu = sp.symbols('gamma_down gamma_up', positive=True)
sX = sp.Matrix([[0, 1], [1, 0]]); sY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
smm = sp.Matrix([[0, 1], [0, 0]]); spp = smm.T
def sdiss(Lk, rate):
    LdL = Lk.H * Lk
    I2s = sp.eye(2)
    return rate * (sp.kronecker_product(Lk, Lk.conjugate()) - sp.Rational(1, 2) * (sp.kronecker_product(LdL, I2s) + sp.kronecker_product(I2s, LdL.T)))
Ds = sdiss(smm, gd) + sdiss(spp, gu)
S1 = sp.kronecker_product(sX, sY.T)                    # vec(XρY) = (X ⊗ Yᵀ) vec ρ; S1 is its own inverse
check("M9 one site, symbolic rates: S·D·S⁻¹ = −D† − (γ↓ + γ↑) for the thermal damping γ↓·D[σ⁻] + γ↑·D[σ⁺], exactly (S = X·ρ·Y reverses the arrow and pairs the heights I with Z and X with Y)",
      sp.simplify(S1 * Ds * S1 + Ds.H + (gd + gu) * sp.eye(4)) == sp.zeros(4, 4) and S1 * S1 == sp.eye(4))
check("M9 one site control: the shift by Z alone (ρ·Z, the same letter pairs without the phases) does not reflect it",
      sp.simplify(sp.kronecker_product(sp.eye(2), sp.diag(1, -1)) * Ds * sp.kronecker_product(sp.eye(2), sp.diag(1, -1)) + Ds.H + (gd + gu) * sp.eye(4)) != sp.zeros(4, 4))
N = 3
SXY = np.kron(string('XXX'), string('YYY').T)
def thermal(N, H, down, up, zdeph=0.0):
    d = 2 ** N; Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        L += diss(N, op(N, l, SM), down[l]) + diss(N, op(N, l, SM.conj().T), up[l])
        if zdeph:
            Zl = op(N, l, Z); L += zdeph * (np.kron(Zl, Zl.T) - np.eye(d * d))
    return L
def xy_refl(L):
    c = -2 * np.trace(L).real / L.shape[0]
    return exact_zero(SXY @ L @ SXY + L.conj().T + c * np.eye(L.shape[0])), c
DOWN = [0.5, 1.0, 1.5]; UP = [0.125, 0.25, 0.375]
XXH = sum(op(N, a, X) @ op(N, b, X) for a, b in chain(N))
DMH = sum(op(N, a, X) @ op(N, b, Y) - op(N, a, Y) @ op(N, b, X) for a, b in chain(N))
okH, cH = xy_refl(thermal(N, heis(N, chain(N)), DOWN, UP))
okX, _ = xy_refl(thermal(N, XXH, DOWN, UP))
okC, _ = xy_refl(thermal(N, heis(N, chain(N)), DOWN, [0, 0, 0]))
check("M9 N=3, site rates γ↓ = (½, 1, 3⁄2), γ↑ = γ↓/4: X^N·ρ·Y^N reflects the full L exactly for the Heisenberg chain and for XX bonds alone "
      "(which do not conserve the excitation number), and at γ↑ = 0; c = Σ(γ↓ + γ↑), so the centre is −Σ(γ↓ + γ↑)/2",
      okH and okX and okC and cH == sum(DOWN) + sum(UP), f"c = {cH}")
badZ, _ = xy_refl(thermal(N, heis(N, chain(N)), DOWN, UP, zdeph=GAM))
badF, _ = xy_refl(thermal(N, heis(N, chain(N)) + sum(op(N, l, Z) for l in range(N)), DOWN, UP))
badD, _ = xy_refl(thermal(N, DMH, DOWN, UP))
check("M9 controls: co-axial Z dephasing, a longitudinal Z field (odd under X^N) and DM bonds XY − YX (odd under X^N) each break the reflection",
      not badZ and not badF and not badD)


# ---------------------------------------------------------------- M10 F137's composition, decided
# Transverse dephasing (X or Y on a site) is a height with no arrow, and S_XY pairs its heights (I with Z, X with Y,
# each pair summing to 2 gamma), so the one-site identity of M9 holds for it too and the composition with the
# damping is a reflection at every N; co-axial Z, a transverse field and a longitudinal field each make an odd
# centred trace Tr((L + c/2)^3) nonzero, which rules out every reflecting element (exact rationals, N = 3).
gx = sp.symbols('gamma_x', positive=True)
for name, P in (("X", sX), ("Y", sY)):
    DP = sdiss(P, gx)
    check(f"M10 one site, symbolic: S_XY·D[{name}]·S_XY⁻¹ = −D[{name}]† − 2γ for {name} dephasing at rate γ (I and Z, X and Y at rates summing to 2γ)",
          sp.simplify(S1 * DP * S1 + DP.H + 2 * gx * sp.eye(4)) == sp.zeros(4, 4))
DZs = sdiss(sp.diag(1, -1), gx)
check("M10 control: Z dephasing is not paired by S_XY (I and Z both at rate 0): S·D[Z]·S⁻¹ + D[Z]† + 2γ is not zero",
      sp.simplify(S1 * DZs * S1 + DZs.H + 2 * gx * sp.eye(4)) != sp.zeros(4, 4))
DT = Ds + sdiss(sX, gx)
check("M10 one site, symbolic: with a heating channel too, S·(D_thermal + D[X])·S⁻¹ = −(…)† − (γ↓ + γ↑ + 2γ)",
      sp.simplify(S1 * DT * S1 + DT.H + (gd + gu + 2 * gx) * sp.eye(4)) == sp.zeros(4, 4))
N = 3
HEIS = heis(N, chain(N))
def damped(H, extra):
    L = thermal(N, H, [1.0, 1.0, 1.0], [0, 0, 0])
    for l, P, r in extra:
        L += diss(N, op(N, l, P), r)
    return L
comp = {"X on every site": [(l, X, 0.5) for l in range(N)], "Y on every site": [(l, Y, 0.5) for l in range(N)],
        "X, Y, X": [(0, X, 0.5), (1, Y, 0.5), (2, X, 0.5)]}
for name, extra in comp.items():
    ok, c = xy_refl(damped(HEIS, extra))
    check(f"M10 N=3, Heisenberg, damping at 1 plus transverse dephasing ({name}) at ½: S_XY reflects the full L exactly, c = Σγ_T1 + 2Σγ_deph = {c} (F137's centre −(Σγ_deph + Σγ_T1/2))",
          ok and c == 3 * 1.0 + 2 * 3 * 0.5)
def odd_trace3(L):
    """Tr((L - Tr(L)/d)^3) in exact complex rationals: every entry of L here is a dyadic complex number, so
    Fraction(float) is exact; the trace is basis-free, so the computational basis serves."""
    from fractions import Fraction
    d = L.shape[0]
    C = lambda z: (Fraction(z.real), Fraction(z.imag))
    tr = (sum(Fraction(L[i, i].real) for i in range(d)), sum(Fraction(L[i, i].imag) for i in range(d)))
    shift = (tr[0] / d, tr[1] / d)
    M = {}
    for i in range(d):
        for j in range(d):
            re, im = C(L[i, j])
            if i == j: re, im = re - shift[0], im - shift[1]
            if re or im: M[(i, j)] = (re, im)
    rows = {}
    for (i, j), v in M.items(): rows.setdefault(i, []).append((j, v))
    def mul(A):
        out = {}
        for (i, k), (ar, ai) in A.items():
            for j, (br, bi) in rows.get(k, []):
                r0, i0 = out.get((i, j), (Fraction(0), Fraction(0)))
                out[(i, j)] = (r0 + ar * br - ai * bi, i0 + ar * bi + ai * br)
        return out
    M3 = mul(mul(M))
    return (sum(M3.get((i, i), (0, 0))[0] for i in range(d)), sum(M3.get((i, i), (0, 0))[1] for i in range(d)))
breaks = {"co-axial Z dephasing at ½": damped(HEIS, [(l, Z, 0.5) for l in range(N)]),
          "a transverse field X": damped(HEIS + sum(op(N, l, X) for l in range(N)), []),
          "a longitudinal field Z": damped(HEIS + sum(op(N, l, Z) for l in range(N)), [])}
pairs = {"damping alone": damped(HEIS, []), "transverse X dephasing": damped(HEIS, comp["X on every site"])}
t_break = {k: odd_trace3(L) for k, L in breaks.items()}
t_pair = {k: odd_trace3(L) for k, L in pairs.items()}
check("M10 N=3: Tr((L + c/2)^3) is exactly nonzero under co-axial Z, a transverse field and a longitudinal field (no element reflects these generators as written; the uniform longitudinal field is the lab frame of the field-free chain, M11), and exactly zero where S_XY reflects",
      all(v != (0, 0) for v in t_break.values()) and all(v == (0, 0) for v in t_pair.values()),
      "; ".join(f"{k}: {v[0]}{'' if v[1] == 0 else f' + {v[1]}i'}" for k, v in t_break.items()))


# ---------------------------------------------------------------- M11 outside S_XY's class: loss and gain exchanged
# Exact over Q(i) (sympy DomainMatrix): the J, D of the bonds and the rates are rational, the turns' phases Pythagorean.
from sympy import QQ_I
cs, sn = sp.symbols('c s', real=True)
def on_circle(e): return sp.simplify(sp.expand(e).subs(sn**2, 1 - cs**2))
g1 = sp.diag(1, cs + sp.I * sn)
def turn_cofactor_independent():
    """[Ad_g, D[X]] = sin θ·M with M a fixed nonzero matrix (free of the angle), so it vanishes exactly where sin θ = 0."""
    comm = (sp.kronecker_product(g1, g1.conjugate()) * sdiss(sX, gx) - sdiss(sX, gx) * sp.kronecker_product(g1, g1.conjugate()))
    M = comm.applyfunc(lambda e: sp.cancel(on_circle(e) / sn))
    return not (M.free_symbols - {gx}) and M != sp.zeros(4, 4) and (comm.applyfunc(on_circle) - sn * M).applyfunc(sp.expand) == sp.zeros(4, 4)
Adg1 = sp.kronecker_product(g1, g1.conjugate())
check("M11 one site, symbolic angle and rates: the turn about Z, g = diag(1, e^{iθ}), commutes with the thermal damping γ↓·D[σ⁻] + γ↑·D[σ⁺] at any angle, "
      "and with X dephasing exactly where sin θ = 0 (the commutator is sin θ times a fixed nonzero matrix)",
      (Adg1 * Ds - Ds * Adg1).applyfunc(on_circle) == sp.zeros(4, 4)
      and (Adg1 * sdiss(sX, gx) - sdiss(sX, gx) * Adg1).applyfunc(on_circle) != sp.zeros(4, 4)
      and turn_cofactor_independent())
ta, tb = sp.symbols('a b', real=True)
Dab = sdiss(ta * sX + tb * sY, gx)
res_ab = (S1 * Dab * S1 + Dab.H + 2 * gx * sp.eye(4)).applyfunc(lambda e: sp.factor(sp.expand(e).subs(tb**2, 1 - ta**2)))
check("M11 one site: S_XY reflects dephasing along a·X + b·Y (a² + b² = 1) only where ab = 0, so a turned element carries transverse dephasing only where the turn is a multiple of π/2",
      res_ab != sp.zeros(4, 4) and (res_ab / (ta * tb)).applyfunc(sp.simplify).free_symbols <= {gx}, f"residual entry {res_ab[1, 2]}")
sZm, sMm = sp.diag(1, -1), sp.Matrix([[0, 1], [0, 0]])
def s_op(n, l, P): return sp.kronecker_product(*[P if j == l else sp.eye(2) for j in range(n)]) if n > 1 else P
def s_string(n, P): return sp.kronecker_product(*[P] * n)
def u1_bond(n, a, b, J, D, Dl=0):
    """J(XX + YY) + D(XY − YX) + Δ·ZZ on (a, b): amplitude 2(J + iD) for the hop from a to b, ⟨01|H|10⟩."""
    return sp.expand(J * (s_op(n, a, sX) * s_op(n, b, sX) + s_op(n, a, sY) * s_op(n, b, sY))
                     + D * (s_op(n, a, sX) * s_op(n, b, sY) - s_op(n, a, sY) * s_op(n, b, sX)) + Dl * s_op(n, a, sZm) * s_op(n, b, sZm))
def s_thermal(n, H, down, up):
    d = 2 ** n; Id = sp.eye(d)
    L = -sp.I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
    for l in range(n):
        for Lk, r in ((s_op(n, l, sMm), down[l]), (s_op(n, l, sMm.T), up[l])):
            LdL = Lk.H * Lk
            L += r * (sp.kronecker_product(Lk, Lk.conjugate()) - sp.Rational(1, 2) * (sp.kronecker_product(LdL, Id) + sp.kronecker_product(Id, LdL.T)))
    return L
def qi(M): return DomainMatrix.from_Matrix(M.applyfunc(sp.expand)).convert_to(QQ_I)
def qi_zero(A): return all(v == 0 for v in A.to_Matrix())
def reflects(S, L, c):
    """S·L + (L† + c)·S = 0, i.e. S·L·S⁻¹ = −L† − c for invertible S."""
    n2 = L.shape[0]
    return qi_zero(qi(S) * qi(L) + (qi(L.H) + qi(c * sp.eye(n2))) * qi(S))
def centred_odd_trace(L, k):
    n2 = L.shape[0]
    M = qi(L) - qi((L.trace() / n2) * sp.eye(n2))
    P = M
    for _ in range(k - 1): P = P * M
    return sp.nsimplify(sum(P.to_Matrix()[i, i] for i in range(n2)))
def turns(n, bonds):
    """the turn about Z at each site that makes the bonds real: e^{iθ_b} = e^{iθ_a}·(J − iD)/|J + iD| along each (a, b)."""
    z = [sp.Integer(1)] * n
    for a, b, J, D, _ in bonds:
        z[b] = sp.expand(z[a] * (J - sp.I * D) / sp.sqrt(J**2 + D**2))
    return sp.diag(*[sp.expand(sp.prod([z[l] for l in range(n) if (i >> (n - 1 - l)) & 1])) for i in range(2 ** n)])
R_ = sp.Rational
N = 3
DOWN3 = [R_(1, 2), R_(1), R_(3, 2)]; UP3 = [R_(1, 8), R_(1, 3), R_(1, 5)]   # unequal, symmetric under no site permutation
C3 = sum(DOWN3) + sum(UP3)
XN3, YN3, ZN3 = s_string(N, sX), s_string(N, sY), s_string(N, sZm)
SXY3 = sp.kronecker_product(XN3, YN3.T)
def triangle(w):
    """unit hopping, real on (0, 1) and (1, 2), J + iD = w = e^{iΦ} on (2, 0), Δ = 0."""
    return u1_bond(N, 0, 1, 1, 0) + u1_bond(N, 1, 2, 1, 0) + u1_bond(N, 2, 0, sp.re(w), sp.im(w))
W35 = (3 + 4 * sp.I) / 5
HT = triangle(W35)
LT_ = s_thermal(N, HT, DOWN3, UP3)
LTc = s_thermal(N, (XN3 * HT * XN3).applyfunc(sp.expand), DOWN3, UP3)
check("M11 N=3 triangle at flux arccos(3/5), unequal rates: X^N·H·X^N = H̄ (every hopping phase reversed) and S_XY·L(H)·S_XY⁻¹ = −L(X^N·H·X^N)† − Σ(γ↓ + γ↑), exactly; "
      "with L(H) itself on the right it fails (control)",
      (XN3 * HT * XN3 - HT.conjugate()).applyfunc(sp.expand).is_zero_matrix
      and qi_zero(qi(SXY3) * qi(LT_) * qi(SXY3) + qi(LTc.H) + qi(C3 * sp.eye(64)))
      and not qi_zero(qi(SXY3) * qi(LT_) * qi(SXY3) + qi(LT_.H) + qi(C3 * sp.eye(64))))
# the chain: the turns make H real and carry S_XY to S' = Ad_G⁻¹∘S_XY∘Ad_G
BONDS = [(0, 1, 4, 3, 2), (1, 2, 12, 5, -1)]
HC = sum((u1_bond(N, a, b, J, D, Dl) for a, b, J, D, Dl in BONDS), sp.zeros(8))
G3 = turns(N, BONDS)
HG = (G3 * HC * G3.H).applyfunc(sp.expand)
LC_ = s_thermal(N, HC, DOWN3, UP3)
SP3 = sp.kronecker_product(G3.H, G3.T) * SXY3 * sp.kronecker_product(G3, G3.conjugate())
check("M11 N=3 chain, bonds (J, D, Δ) = (4, 3, 2) and (12, 5, −1), unequal rates: the turns e^{iθ_b} = e^{iθ_a}·(J − iD)/√(J² + D²) make G·H·G† real and commuting with X^N and Y^N, "
      "and S' = Ad_G⁻¹∘S_XY∘Ad_G, ρ ↦ ⊗(g†Xg)·ρ·⊗(g†Yg), reflects the full L exactly; the bare S_XY does not",
      all(sp.im(v) == 0 for v in HG) and (HG * XN3 - XN3 * HG).is_zero_matrix and (HG * YN3 - YN3 * HG).is_zero_matrix
      and reflects(SP3, LC_, C3) and not reflects(SXY3, LC_, C3))
# the Sextant's DM chain under loss (XY − YX on each bond, loss at rate 1): the turns are quarter turns, the element a pair of strings
BONDS_DM = [(0, 1, 0, 1, 0), (1, 2, 0, 1, 0)]
HD = sum((u1_bond(N, a, b, J, D, Dl) for a, b, J, D, Dl in BONDS_DM), sp.zeros(8))
GD = turns(N, BONDS_DM)
LD_ = s_thermal(N, HD, [1, 1, 1], [0, 0, 0])
AL, AR = (GD.H * XN3 * GD).applyfunc(sp.expand), (GD.H * YN3 * GD).applyfunc(sp.expand)
def pauli_string_up_to_sign(M):
    for word in product('IXYZ', repeat=N):
        P = sp.kronecker_product(*[{'I': sp.eye(2), 'X': sX, 'Y': sY, 'Z': sZm}[c] for c in word])
        if (M - P).is_zero_matrix or (M + P).is_zero_matrix: return ''.join(word)
    return None
SPD = sp.kronecker_product(AL, AR.T)
check("M11 the Sextant's DM chain under loss (N=3, XY − YX on each bond, σ⁻ at rate 1): the turns are quarter turns, the element is ρ ↦ A·ρ·B with A = XYX and B = YXY up to sign, "
      "and it reflects the full L exactly; the bare S_XY does not",
      pauli_string_up_to_sign(AL) == "XYX" and pauli_string_up_to_sign(AR) == "YXY"
      and reflects(SPD, LD_, 3) and not reflects(SXY3, LD_, 3), f"A = {pauli_string_up_to_sign(AL)}, B = {pauli_string_up_to_sign(AR)}")
# the triangle at flux pi/2: the reflection fixing site 1 reverses the flux; with rates symmetric under it, Ad_P∘S_XY reflects
PERM = sp.zeros(8)
for i in range(8):
    b0, b1, b2 = (i >> 2) & 1, (i >> 1) & 1, i & 1
    PERM[(b2 << 2) | (b1 << 1) | b0, i] = 1
HQ = triangle(sp.I)
DOWN_S = [R_(1, 2), R_(1), R_(1, 2)]; UP_S = [R_(1, 8), R_(1, 3), R_(1, 8)]
EP = sp.kronecker_product(PERM, PERM) * SXY3
check("M11 N=3 triangle at flux π/2 and arccos(3/5): the site reflection P fixing site 1 gives P·H·P† = X^N·H·X^N, and Ad_P∘S_XY reflects the full L exactly at rates symmetric under P; at flux π/2 not at the unequal rates",
      all((PERM * triangle(w) * PERM.T - XN3 * triangle(w) * XN3).applyfunc(sp.expand).is_zero_matrix
          and reflects(EP, s_thermal(N, triangle(w), DOWN_S, UP_S), sum(DOWN_S) + sum(UP_S)) for w in (sp.I, W35))
      and not reflects(EP, s_thermal(N, HQ, DOWN3, UP3), C3))
# the general second carrier, V·(X^N·H·X^N)·V† = H, with fields: the Heisenberg chain with a field profile the reflection negates
HHF = (u1_bond(N, 0, 1, 1, 0, 1) + u1_bond(N, 1, 2, 1, 0, 1) - R_(1, 2) * s_op(N, 0, sZm) + R_(1, 2) * s_op(N, 2, sZm)).applyfunc(sp.expand)
check("M11 N=3 Heisenberg chain with fields (−½, 0, ½): P·(X^N·H·X^N)·P† = H for the reflection P, and Ad_P∘S_XY reflects the full L exactly at rates P keeps "
      "(cooling (½, 1, ½) with heating (⅛, ⅓, ⅛) and under cooling alone), not at γ↓ = (½, 1, 3⁄2) under cooling alone",
      (PERM * XN3 * HHF * XN3 * PERM.T - HHF).applyfunc(sp.expand).is_zero_matrix
      and reflects(EP, s_thermal(N, HHF, DOWN_S, UP_S), sum(DOWN_S) + sum(UP_S))
      and reflects(EP, s_thermal(N, HHF, DOWN_S, [0, 0, 0]), sum(DOWN_S))
      and not reflects(EP, s_thermal(N, HHF, DOWN3, [0, 0, 0]), sum(DOWN3)))
# infinite temperature: gamma_up = gamma_down on every site, the spin flip is invisible, the element is rho·Z^N
SZ3 = sp.kronecker_product(sp.eye(8), ZN3.T)
inf_ok = all(reflects(SZ3, s_thermal(N, triangle(w), DOWN3, DOWN3), 2 * sum(DOWN3)) for w in (sp.I, W35))
off_tr = centred_odd_trace(s_thermal(N, HQ, DOWN3, [R_(1, 2), R_(1), R_(1, 2)]), 5)
check("M11 N=3 triangle, γ↑ = γ↓ = (½, 1, 3⁄2) on the sites: ρ ↦ ρ·Z^N (Ad_{X^N}∘S_XY up to a phase) reflects the full L exactly at flux π/2 and arccos(3/5); "
      "with site 2's heating moved to ½ the fifth trace at flux π/2 is nonzero",
      inf_ok and off_tr != 0, f"fifth trace off the end: {off_tr}")
# one ratio gamma_up/gamma_down on every site; the thermal transpose (F119 §7) holds for complex H too
xs = sp.symbols('x')
def palindromic_charpoly(L, c):
    pz = sp.Poly([QQ_I.to_sympy(cf) for cf in qi(L).charpoly()], xs).as_expr()
    q = sp.expand(pz.subs(xs, -c - xs))
    return sp.expand(q - pz) == 0   # even degree (4^N), so a palindrome is q = p
HW = u1_bond(N, 0, 1, 1, 0, R_(7, 10)) + u1_bond(N, 1, 2, 2, 0, R_(7, 10)) + u1_bond(N, 2, 0, 0, 3, R_(7, 10))
UP_Q = [x / 4 for x in DOWN3]
LW = s_thermal(N, HW, DOWN3, UP_Q)
SW = sp.diag(*[R_(1, 2) ** bin(i).count('1') for i in range(2 ** N)])   # pi_beta^(1/2), e^(-beta) = 1/4
XD3 = sum((R_(3, 10) * (sp.kronecker_product(s_op(N, l, sX), s_op(N, l, sX).T) - sp.eye(64)) for l in range(N)), sp.zeros(64))
GW = sp.kronecker_product(SW, SW)                                          # rho -> pi^(1/2) rho pi^(1/2)
check("M11 N=3 triangle, one ratio (γ↑ = γ↓/4 on every site, γ↓ = (½, 1, 3⁄2)), hoppings 1, 2 and 3i, Δ = 7⁄10, flux π/2: the thermal transpose holds for this complex H, "
      "L(H)ᵀ = Γ⁻¹·L(Hᵀ)·Γ with Γ: ρ ↦ π^½ρπ^½, Hᵀ = H̄ = X^N·H·X^N, and the characteristic polynomial is palindromic about −c/2 (exact); "
      "L(H) in place of L(Hᵀ) is the identity's control, and the same H with γ↑ = (⅛, ⅓, ⅕), three ratios, gives Tr((L + c/2)^5) = −3936 (control)",
      (LW.T - GW.inv() * s_thermal(N, HW.T, DOWN3, UP_Q) * GW).applyfunc(sp.expand).is_zero_matrix
      and not (LW.T - GW.inv() * s_thermal(N, HW, DOWN3, UP_Q) * GW).applyfunc(sp.expand).is_zero_matrix
      and (HW.T - XN3 * HW * XN3).applyfunc(sp.expand).is_zero_matrix
      and palindromic_charpoly(LW, sum(DOWN3) + sum(UP_Q))
      and centred_odd_trace(s_thermal(N, HW, DOWN3, UP3), 5) == -3936)
check("M11 N=3 unit triangle at flux π/2 (the H of the fifth-trace row): at one ratio, γ↑ = γ↓/4, the characteristic polynomial is palindromic, and at ratio 4, γ↑ = 4γ↓, too (exact); "
      "with X dephasing at 3⁄10 on every site beside the one-ratio damping it pairs at flux 0 and not at π/2, the seventh power sum about the mean −263277/400 (the weight does not carry dephasing)",
      palindromic_charpoly(s_thermal(N, HQ, DOWN3, UP_Q), sum(DOWN3) + sum(UP_Q))
      and palindromic_charpoly(s_thermal(N, HQ, DOWN3, [4 * x for x in DOWN3]), 5 * sum(DOWN3))
      and palindromic_charpoly(s_thermal(N, triangle(sp.Integer(1)), DOWN3, UP_Q) + XD3, sum(DOWN3) + sum(UP_Q) + R_(9, 5))
      and centred_odd_trace(s_thermal(N, HQ, DOWN3, UP_Q) + XD3, 7) == R_(-263277, 400))
# a bipartite ring at Delta = 0: complex conjugation sends L(H) to L(-Hbar); Z on one sublattice sends -hopping to hopping
N4 = 4
def u1_bond4(a, b, J, D, Dl=0):
    return sp.expand(J * (s_op(N4, a, sX) * s_op(N4, b, sX) + s_op(N4, a, sY) * s_op(N4, b, sY))
                     + D * (s_op(N4, a, sX) * s_op(N4, b, sY) - s_op(N4, a, sY) * s_op(N4, b, sX)) + Dl * s_op(N4, a, sZm) * s_op(N4, b, sZm))
H4 = (u1_bond4(0, 1, 1, 0) + u1_bond4(1, 2, 2, 0) + u1_bond4(2, 3, 3, 0) + u1_bond4(3, 0, 0, R_(17, 10))
      + R_(3, 10) * s_op(N4, 0, sZm) - R_(1, 2) * s_op(N4, 2, sZm))
D4 = [R_(1, 2), R_(1), R_(3, 2), R_(2)]; U4 = [R_(1, 8), R_(1, 3), R_(1, 5), R_(1, 7)]
XN4, W4 = s_string(N4, sX), s_op(N4, 1, sZm) * s_op(N4, 3, sZm)
H4z = H4 + R_(7, 10) * sum((s_op(N4, a, sZm) * s_op(N4, (a + 1) % N4, sZm) for a in range(N4)), sp.zeros(16))
check("M11 N=4 ring, Δ = 0, hoppings 1, 2, 3 and 1.7i (flux π/2), longitudinal fields on two sites, four ratios: the characteristic polynomial is palindromic (exact), "
      "carried by complex conjugation (L(H) entry by entry conjugated is L(−H̄), exact) and W = Z on sites 1 and 3 (W·(−conj(X^N·H·X^N))·W† = H, exact); with Δ = 7⁄10 it is not",
      palindromic_charpoly(s_thermal(N4, H4, D4, U4), sum(D4) + sum(U4))
      and not palindromic_charpoly(s_thermal(N4, H4z, D4, U4), sum(D4) + sum(U4))
      and s_thermal(N4, H4, D4, U4).conjugate() == s_thermal(N4, (-H4.conjugate()).applyfunc(sp.expand), D4, U4)
      and (W4 * (-(XN4 * H4 * XN4).conjugate()) * W4 - H4).applyfunc(sp.expand).is_zero_matrix)
# the arrow exchange: for any H commuting with Z^N, X or Y dephasing beside, X^N carries L(X^N·H·X^N) to L(H) with cooling and heating exchanged
HA = (triangle(W35) + u1_bond(N, 1, 2, 1, 0) + R_(1, 3) * s_op(N, 0, sZm) - R_(2, 5) * s_op(N, 2, sZm)
      + R_(1, 2) * (s_op(N, 0, sX) * s_op(N, 1, sX) - s_op(N, 0, sY) * s_op(N, 1, sY))).applyfunc(sp.expand)
XDA = sum((r * (sp.kronecker_product(s_op(N, l, sX), s_op(N, l, sX).T) - sp.eye(64)) for l, r in ((0, R_(3, 10)), (2, R_(1, 5)))), sp.zeros(64))
LA = s_thermal(N, HA, DOWN3, UP3) + XDA
LAx = s_thermal(N, (XN3 * HA * XN3).applyfunc(sp.expand), DOWN3, UP3) + XDA
LAswap = s_thermal(N, HA, UP3, DOWN3) + XDA
AdX3 = sp.kronecker_product(XN3, XN3)
CA = C3 + 2 * (R_(3, 10) + R_(1, 5))
tA, tAs = centred_odd_trace(LA, 5), centred_odd_trace(LAswap, 5)
check("M11 N=3, H with a hopping phase, Z fields on two sites and a term XX − YY that does not conserve the excitation number (commuting with Z^N), X dephasing on two sites: "
      "S_XY·L·S_XY⁻¹ = −L(X^N·H·X^N)† − c and Ad_{X^N}·L(X^N·H·X^N)·Ad_{X^N} = L(H) with cooling and heating exchanged, exactly; the fifth trace changes sign under the exchange",
      (HA * ZN3 - ZN3 * HA).is_zero_matrix
      and qi_zero(qi(SXY3) * qi(LA) * qi(SXY3) + qi(LAx.H) + qi(CA * sp.eye(64)))
      and (AdX3 * LAx * AdX3 - LAswap).applyfunc(sp.expand).is_zero_matrix
      and tA == -tAs and tA != 0, f"fifth trace {tA} against {tAs}")
check("M11 N=3, the same H and dephasing at γ↑ = γ↓ = (½, 1, 3⁄2), no arrow: ρ ↦ ρ·Z^N reflects the full L exactly, fields and the non-conserving term included; at the unequal rates it does not",
      reflects(SZ3, s_thermal(N, HA, DOWN3, DOWN3) + XDA, 2 * sum(DOWN3) + 2 * (R_(3, 10) + R_(1, 5)))
      and not reflects(SZ3, s_thermal(N, HA, DOWN3, UP3) + XDA, CA))
# the spin flip composed with a site permutation that keeps the flux and exchanges each site's cooling and heating
N4b = 4
def u1_ring4(ts, Dl):
    return sum((u1_bond4(a, (a + 1) % N4b, sp.re(t), sp.im(t), Dl) for a, t in enumerate(ts)), sp.zeros(16))
H4f = u1_ring4([W35, 2, W35, 2], R_(7, 10))
D4f = [R_(1, 2), R_(1), R_(3, 2), R_(2)]; U4f = [R_(3, 2), R_(2), R_(1, 2), R_(1)]
PERM4 = sp.zeros(16)
for i in range(16):
    bits = [(i >> (N4b - 1 - l)) & 1 for l in range(N4b)]
    nb = [0] * N4b
    for l in range(N4b): nb[(l + 2) % N4b] = bits[l]
    PERM4[sum(b << (N4b - 1 - l) for l, b in enumerate(nb)), i] = 1
XN4s, YN4s = s_string(N4b, sX), s_string(N4b, sY)
V4 = PERM4 * XN4s
E4 = sp.kronecker_product(V4, V4.conjugate()) * sp.kronecker_product(XN4s, YN4s.T)
U4bad = [R_(3, 2), R_(2), R_(1, 2), R_(3, 4)]
check("M11 N=4 ring, Δ = 7⁄10, hoppings (3 + 4i)/5, 2, (3 + 4i)/5, 2, ratios 3, 2, ⅓, ½: V = P·X^N with P the shift by two sites keeps H̄ onto H, and Ad_V∘S_XY reflects the full L exactly; "
      "with site 3's heating at ¾ it does not",
      (V4 * H4f.conjugate() * V4.H - H4f).applyfunc(sp.expand).is_zero_matrix
      and reflects(E4, s_thermal(N4b, H4f, D4f, U4f), sum(D4f) + sum(U4f))
      and not reflects(E4, s_thermal(N4b, H4f, D4f, U4bad), sum(D4f) + sum(U4bad)))
# the shift S_IZ: ρ ↦ ρ·Z^N: for H commuting with Z^N and jumps odd under Z^N with A†A + AA† ∝ 1, S_IZ·L·S_IZ⁻¹ = −L(J → J†)† − c
TILT = [(R_(3, 5), R_(4, 5)), (R_(1), R_(0)), (R_(5, 13), R_(-12, 13))]; GT = [R_(3, 10), R_(1, 4), R_(1, 5)]
HR = (HA + R_(3, 7) * s_op(N, 0, sZm) * s_op(N, 2, sZm)).applyfunc(sp.expand)
def tilted(n_down, n_up):
    L = s_thermal(N, HR, n_down, n_up)
    for l, ((a, b), g) in enumerate(zip(TILT, GT)):
        J = s_op(N, l, a * sX + b * sY); L += g * (sp.kronecker_product(J, J.conjugate()) - sp.eye(64))
    return L
LR, LRs = tilted(DOWN3, UP3), tilted(UP3, DOWN3)
CR = C3 + 2 * sum(GT)
check("M11 N=3, the same H with a ZZ bond added, damping at unequal rates and dephasing along (3X + 4Y)/5, X and (5X − 12Y)/13 on the three sites: "
      "S_IZ: ρ ↦ ρ·Z^N carries L onto −L′† − c with L′ the Liouvillian of cooling and heating exchanged, exactly; with L itself on the right it fails (control)",
      qi_zero(qi(SZ3) * qi(LR) * qi(SZ3) + qi(LRs.H) + qi(CR * sp.eye(64)))
      and not qi_zero(qi(SZ3) * qi(LR) * qi(SZ3) + qi(LR.H) + qi(CR * sp.eye(64))))
# one ratio per component of the hopping bonds: a pendant joined by ZZ alone may carry its own ratio
def pendant(join_by_hopping):
    H = (u1_bond4(0, 1, 1, 0) + u1_bond4(1, 2, 1, 0) + u1_bond4(2, 0, 0, 1)).applyfunc(sp.expand)
    H += u1_bond4(3, 0, 1, 0, 1) if join_by_hopping else u1_bond4(3, 0, 0, 0, 1)
    return H
DP = [R_(1, 2), R_(1), R_(3, 2), R_(1)]; UPp = [R_(1, 8), R_(1, 4), R_(3, 8), R_(1, 9)]
H4zz = (u1_bond4(0, 1, 1, 0, R_(7, 10)) + u1_bond4(1, 2, 2, 0, R_(7, 10)) + u1_bond4(2, 3, 3, 0, R_(7, 10)) + u1_bond4(3, 0, 0, R_(17, 10), R_(7, 10))).applyfunc(sp.expand)
check("M11 N=4: the unit triangle at flux π/2 at ratio ¼ with a fourth site at ratio 1⁄9 joined by a ZZ bond alone is palindromic, joined by a hopping bond it is not; "
      "the ring of the conjugation row at Δ = 7⁄10 without fields, four ratios, flux π/2, is not palindromic (exact characteristic polynomials)",
      palindromic_charpoly(s_thermal(N4, pendant(False), DP, UPp), sum(DP) + sum(UPp))
      and not palindromic_charpoly(s_thermal(N4, pendant(True), DP, UPp), sum(DP) + sum(UPp))
      and not palindromic_charpoly(s_thermal(N4, H4zz, D4, U4), sum(D4) + sum(U4)))
# transverse DM does not commute with Z^N: the three-site chain breaks under cooling alone, the axial chain does not
HTD = sum((s_op(N, a, sY) * s_op(N, b, sZm) - s_op(N, a, sZm) * s_op(N, b, sY) for a, b in ((0, 1), (1, 2))), sp.zeros(8)).applyfunc(sp.expand)
LTD = s_thermal(N, HTD, [1, 1, 1], [0, 0, 0])
check("M11 N=3 open chain, transverse DM YZ − ZY on both bonds (not commuting with Z^N), cooling alone at rate 1: Tr((L + c/2)^5) = 240, a break on a tree; "
      "the same chain at equal cooling and heating rates 1, where every H commuting with Z^N pairs, gives 1280; the axial DM chain under the same cooling has its odd power sums to 7 at zero",
      not (HTD * ZN3 - ZN3 * HTD).is_zero_matrix and centred_odd_trace(LTD, 5) == 240
      and centred_odd_trace(s_thermal(N, HTD, [1, 1, 1], [1, 1, 1]), 5) == 1280
      and all(centred_odd_trace(LD_, k) == 0 for k in (3, 5, 7)))
# the conjugation carrier at any dephasing axes in the XY plane: the bipartite ring at Δ = 0
def tilted4(H):
    L = s_thermal(N4, H, D4, U4)
    for l, (a, b), g in ((0, (R_(3, 5), R_(4, 5)), R_(3, 10)), (1, (R_(5, 13), R_(-12, 13)), R_(1, 4)), (2, (R_(1), R_(0)), R_(1, 5))):
        J = s_op(N4, l, a * sX + b * sY); L += g * (sp.kronecker_product(J, J.conjugate()) - sp.eye(256))
    return L
C4T = sum(D4) + sum(U4) + 2 * (R_(3, 10) + R_(1, 4) + R_(1, 5))
check("M11 N=4 ring of the conjugation row (Δ = 0, flux π/2, fields, four ratios) with dephasing along (3X + 4Y)/5, (5X − 12Y)/13 and X on three sites: "
      "palindromic (exact characteristic polynomial); at Δ = 7⁄10 not",
      palindromic_charpoly(tilted4(H4), C4T) and not palindromic_charpoly(tilted4(H4z), C4T))
# the triangle at unequal rates: the fifth odd trace as a function of the flux
TS = [R_(0), R_(1, 2), R_(1, 3), R_(2), R_(3), R_(1, 4), R_(4), R_(2, 3), R_(3, 2), R_(1, 5), R_(5), R_(1)]
tr5, tr3 = [], []
for t in TS:
    w = sp.expand((1 + sp.I * t) ** 2 / (1 + t ** 2))
    Lw = s_thermal(N, triangle(w), DOWN3, UP3)
    tr5.append((w, centred_odd_trace(Lw, 5))); tr3.append(centred_odd_trace(Lw, 3))
check("M11 N=3 triangle, unequal rates: Tr((L + c/2)^5) = −656·sin Φ at twelve points e^{iΦ} = (1 + it)²/(1 + t²), t rational, Φ in [0, π) (eleven fix a trigonometric polynomial of degree 5, "
      "L being linear in e^{±iΦ}), so the trace is −656·sin Φ at every flux and nothing reflects the spectrum off Φ = 0, π; Tr((L + c/2)^3) = 0 at the twelve",
      all(v == -656 * sp.im(w) for w, v in tr5) and all(v == 0 for v in tr3), f"Φ = π/2: {[v for w, v in tr5 if w == sp.I][0]}")
check("M11 N=3 triangle: with cooling and heating exchanged on every site the fifth trace at flux π/2 is +656; at flux π (the bond −1, real) S_XY itself reflects at the three ratios, at flux π/2 it does not",
      centred_odd_trace(s_thermal(N, HQ, UP3, DOWN3), 5) == 656 and reflects(SXY3, s_thermal(N, triangle(sp.Integer(-1)), DOWN3, UP3), C3)
      and not reflects(SXY3, s_thermal(N, HQ, DOWN3, UP3), C3))
WT = (3 + 4 * sp.I) / 5
HG0 = (u1_bond(N, 0, 1, sp.re(WT), sp.im(WT)) + u1_bond(N, 1, 2, sp.re(WT), -sp.im(WT)) + u1_bond(N, 2, 0, 1, 0)).applyfunc(sp.expand)
tg9 = centred_odd_trace(s_thermal(N, HG0, DOWN3, UP_Q) + XD3, 9)
check("M11 N=3 triangle at one ratio with X dephasing at 3⁄10: flux 0 written with hopping phases (θ, −θ, 0), cos θ = 3⁄5, does not pair (a nonzero ninth power sum), so X dephasing sees the gauge; with real hopping it pairs (control)",
      tg9 != 0 and palindromic_charpoly(s_thermal(N, triangle(sp.Integer(1)), DOWN3, UP_Q) + XD3, sum(DOWN3) + sum(UP_Q) + R_(9, 5)), f"ninth trace {tg9}")
# the decay rates themselves: a reading with two eigensolvers, numpy in double precision and mpmath at 30 digits
import mpmath as mp
def to_mp(z):
    a, b = sp.Rational(sp.re(z)), sp.Rational(sp.im(z))
    return mp.mpc(mp.mpf(a.p) / a.q, mp.mpf(b.p) / b.q)
def rate_miss(L, c):
    re = np.sort(np.linalg.eigvals(np.array(L.evalf(), dtype=complex)).real)
    miss = float(np.max(np.abs(re + re[::-1] + float(c))))
    with mp.workdps(30):
        ev = mp.eig(mp.matrix([[to_mp(L[i, j]) for j in range(L.shape[1])] for i in range(L.shape[0])]), left=False, right=False)
        r30 = sorted(mp.re(e) for e in ev)
        return miss, float(max(abs(x + y + to_mp(c).real) for x, y in zip(r30, reversed(r30))))
(m_flux, m_flux30), (m_zero, m_zero30) = rate_miss(s_thermal(N, triangle(sp.I), DOWN3, UP3), C3), rate_miss(s_thermal(N, triangle(sp.Integer(1)), DOWN3, UP3), C3)
print(f"READ M11 N=3 triangle, three ratios (two eigensolvers, no pass or fail): the sorted decay rates miss their mirror about −c/2 by "
      f"{m_flux:.6e} in double precision and {m_flux30:.6e} at 30 digits at flux π/2, by {m_zero:.1e} and {m_zero30:.1e} at flux 0", flush=True)
LDet = s_thermal(N, (triangle(sp.Integer(1)) - R_(1, 2) * s_op(N, 0, sZm) + R_(1, 2) * s_op(N, 2, sZm)).applyfunc(sp.expand), DOWN3, [0, 0, 0])
missD, missD30 = rate_miss(LDet, sum(DOWN3))
print(f"READ M11 N=3 triangle, flux 0, fields (−½, 0, ½), cooling alone (two eigensolvers, no pass or fail): the sorted decay rates miss their mirror by "
      f"{missD:.6f} in double precision and {missD30:.6f} at 30 digits; the rates are what no common frame shift moves", flush=True)
# the three chained, with fields: at one ratio S_XY, the transpose and the conjugation make the palindrome
# spec L(H + F) = spec L(−H̄ + F), and turns with W·(−H̄)·W† = H carry it (on a triangle at Δ = 0: flux ±π/2)
def centred_even(L):
    n2 = L.shape[0]; mu = sp.nsimplify(L.trace() / n2)
    cf = qi(L - mu * sp.eye(n2)).charpoly()
    return all(QQ_I.to_sympy(cf[j]) == 0 for j in range(1, len(cf), 2))   # a domain zero is not the integer 0
def fields3(h): return sum((h[l] * s_op(N, l, sZm) for l in range(N)), sp.zeros(2 ** N))
def triangle_zz(w, Dl): return u1_bond(N, 0, 1, 1, 0, Dl) + u1_bond(N, 1, 2, 1, 0, Dl) + u1_bond(N, 2, 0, sp.re(w), sp.im(w), Dl)
hF = [R_(3, 10), -R_(1, 5), R_(1, 2)]
pairs_c = [centred_even(s_thermal(N, (triangle(w) + fields3([-R_(1, 2), 0, R_(1, 2)])).applyfunc(sp.expand), DOWN3, [0, 0, 0])) for w in (sp.I, -sp.I)] + [
           centred_even(s_thermal(N, (triangle(w) + fields3(hF)).applyfunc(sp.expand), DOWN3, UP_Q)) for w in (sp.I, -sp.I)]
breaks_c = [centred_even(s_thermal(N, (triangle_zz(sp.I, R_(7, 10)) + fields3(hF)).applyfunc(sp.expand), DOWN3, UP_Q)),
            centred_even(s_thermal(N, (triangle(W35) + fields3(hF)).applyfunc(sp.expand), DOWN3, UP_Q)),
            centred_even(s_thermal(N, (triangle(sp.I) + fields3(hF)).applyfunc(sp.expand), DOWN3, UP3))]
check("M11 the three chained, with fields: on the triangle at flux ±π/2 the characteristic polynomial is palindromic with fields (−½, 0, ½) "
      "under cooling alone (which break at flux 0, −2560 below) and with fields (3⁄10, −⅕, ½) at one ratio ¼, at both signs of the flux; "
      "Δ = 7⁄10, flux arccos(3/5) and three ratios break it (controls)",
      all(pairs_c) and not any(breaks_c), f"pairs {pairs_c}, controls {breaks_c}")
# no common frame pairs the detuned triangle: Tr5 is a polynomial of degree at most 5 in the frame offset κ, fixed by six offsets
kk = sp.symbols('kk')
ZSUM3 = sum((s_op(N, l, sZm) for l in range(N)), sp.zeros(2 ** N))
def l_det(k): return s_thermal(N, (triangle(sp.Integer(1)) - R_(1, 2) * s_op(N, 0, sZm) + R_(1, 2) * s_op(N, 2, sZm) + k * ZSUM3).applyfunc(sp.expand), DOWN3, [0, 0, 0])
t5_line = sp.expand(sp.interpolate([(k, centred_odd_trace(l_det(k), 5)) for k in range(6)], kk))
t5_roots = sp.solve(t5_line, kk)
t7_at_roots = [centred_odd_trace(l_det(r), 7) for r in t5_roots]
check("M11 the detuned triangle (real hopping, fields (−½, 0, ½), cooling alone at γ↓ = (½, 1, 3⁄2)) pairs in no common frame: with κ·Σ Z_l added, "
      "Tr((L + c/2)^5) = −2560 − 15360κ (six offsets fix the degree-5 polynomial), zero only at κ = −1/6, where Tr((L + c/2)^7) = 1708448/9",
      t5_line == -2560 - 15360 * kk and t5_roots == [-R_(1, 6)] and t7_at_roots == [R_(1708448, 9)], f"Tr5 = {t5_line}, Tr7 at its zero {t7_at_roots}")
# transverse composition beside cooling: X or Y, or one orthogonal pair of axes (a common turn about Z), pairs; a tilted axis beside X can break
def deph3(axes):
    M = sp.zeros(4 ** N)
    for l, (n_ax, g) in enumerate(axes):
        if g:
            A = s_op(N, l, n_ax); M += g * (sp.kronecker_product(A, A.T) - sp.eye(4 ** N))
    return M
TI, TP = (3 * sX + 4 * sY) / 5, (-4 * sX + 3 * sY) / 5
LC1 = s_thermal(N, (u1_bond(N, 0, 1, 1, 0, 1) + u1_bond(N, 1, 2, 1, 0, 1)).applyfunc(sp.expand), [1, 1, 1], [0, 0, 0])
check("M11 transverse composition beside cooling, Heisenberg chain N=3: dephasing at ½ along the orthogonal pair (3X + 4Y)/5, (−4X + 3Y)/5, (3X + 4Y)/5 pairs "
      "(a common turn about Z brings it to X and Y), while (3X + 4Y)/5 on site 0 beside X on site 1 does not (an odd coefficient of the centred characteristic polynomial)",
      centred_even(LC1 + deph3([(TI, R_(1, 2)), (TP, R_(1, 2)), (TI, R_(1, 2))])) and not centred_even(LC1 + deph3([(TI, R_(1, 2)), (sX, R_(1, 2)), (sX, 0)])))
# the lab frame: a uniform field is a constant on each coherence order p − q
pcN = [bin(i).count('1') for i in range(2 ** N)]
order = [pcN[i // 2 ** N] - pcN[i % 2 ** N] for i in range(4 ** N)]
def order_mirrors(L):
    """L block-diagonal in the coherence order, and each block's characteristic polynomial, centred on the block's mean, even or odd about it; returns the means"""
    if any(order[i] != order[j] for (i, j) in L.todok()):
        return None
    cents = []
    for m in range(-N, N + 1):
        idx = [i for i in range(4 ** N) if order[i] == m]
        B = L.extract(idx, idx); mu = sp.nsimplify(B.trace() / len(idx)); cents.append(mu)
        cf = qi(B - mu * sp.eye(len(idx))).charpoly()
        if any(QQ_I.to_sympy(cf[j]) != 0 for j in range(1, len(cf), 2)):   # a domain zero is not the integer 0
            return None
    return cents
HLab = (u1_bond(N, 0, 1, 1, 0, 1) + u1_bond(N, 1, 2, 1, 0, 1) + sum((s_op(N, l, sZm) for l in range(N)), sp.zeros(2 ** N))).applyfunc(sp.expand)
cLab = order_mirrors(s_thermal(N, HLab, [1] * N, [0] * N))
cDet = order_mirrors(LDet)
check("M11 lab frame, N=3: the Heisenberg chain under cooling alone at rate 1 with the uniform field Σ Z_l (M10's uniform-field case, its lab frame): L is block-diagonal in the coherence order "
      "and each order's characteristic polynomial, centred on the order's mean −c/2 + iκ, is even or odd about it, κ growing with the order, so each order keeps its own mirror; "
      "the triangle with fields (−½, 0, ½) under cooling alone does not (control)",
      cLab is not None and all(sp.re(z) == -sp.Rational(N, 2) for z in cLab) and all(sp.im(cLab[k]) < sp.im(cLab[k + 1]) for k in range(2 * N))
      and cDet is None, f"means {cLab}")
check("M11 N=3 triangle with real hopping (flux 0) and fields Σ h_l·Z_l, h = (−½, 0, ½), cooling alone at γ↓ = (½, 1, 3⁄2): Tr((L + c/2)^5) = −2560 in this frame, so the complex spectrum does not pair here, and a uniform field ½ added (another common frame) moves that power sum; "
      "without the fields the odd power sums to 7 vanish",
      centred_odd_trace(LDet, 5) == -2560
      and centred_odd_trace(LDet - sp.I * (sp.kronecker_product(R_(1, 2) * sum((s_op(N, l, sZm) for l in range(N)), sp.zeros(8)), sp.eye(8))
                                            - sp.kronecker_product(sp.eye(8), (R_(1, 2) * sum((s_op(N, l, sZm) for l in range(N)), sp.zeros(8))).T)), 5) != -2560
      and all(centred_odd_trace(s_thermal(N, triangle(sp.Integer(1)), DOWN3, [0, 0, 0]), k) == 0 for k in (3, 5, 7)))
LQc = s_thermal(N, HQ, DOWN3, [0, 0, 0])
check("M11 N=3 triangle at flux π/2, γ↓ = (½, 1, 3⁄2), cooling alone: the odd traces to 7 vanish (a necessary reading; the proof at every flux is F137's route through H_eff), while with γ↑ = (⅛, ⅓, ⅕) the fifth is −656",
      all(centred_odd_trace(LQc, k) == 0 for k in (3, 5, 7)) and centred_odd_trace(s_thermal(N, HQ, DOWN3, UP3), 5) == -656)
def z_dephased(H):
    L = -sp.I * (sp.kronecker_product(H, sp.eye(8)) - sp.kronecker_product(sp.eye(8), H.T))
    for l in range(N):
        Zl = s_op(N, l, sZm); L += R_(1, 2) * (sp.kronecker_product(Zl, Zl.T) - sp.eye(64))
    return L
check("M11 contrast, N=3 triangle under uniform Z dephasing at ½: the characteristic polynomial is not palindromic about −σ at flux π/2 and is at flux π (exact), "
      "the one-sided ρ·X^N carrying L(H) onto a reflection only where X^N·H·X^N = H (derived; the rows compare characteristic polynomials)",
      not palindromic_charpoly(z_dephased(HQ), 3) and palindromic_charpoly(z_dephased(triangle(sp.Integer(-1))), 3))

print()
print("FAILED: " + ", ".join(FAIL) if FAIL else "ALL GATES PASS")
sys.exit(1 if FAIL else 0)
