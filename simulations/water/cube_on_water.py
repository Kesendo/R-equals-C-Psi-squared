"""The letter cube laid on the water wire: three rows of THE_ONE_SQUARE section 9 sighted on the wire, gated.

The site is a proton's position in its hydrogen bond, |L>, |R> the Z basis; the dephasing letter is Z, so a
coherence |i><j| decays at -2 gamma k_Z with k_Z = popcount(i xor j), the number of linkages in which the two
proton configurations it links disagree (the dipole mu = popcount is diagonal and sits on the face k_Z = 0; k_Z is
the disagreement of two configurations, not their dipole difference: |01><10| has mu = 1 on both sides and k_Z = 2).
Tunnelling -J X_l is a field (odd k_Z, even k_X), the Ising pair term K Z Z keeps k_Z, the bias Delta Z is a Z field
(odd k_X).

C1  Row 1, the lit corner shift S: rho -> rho X^N. On the two hydrogen-bond-qubit models (N = 2: -J(X1+X2) + K Z1Z2;
    N = 4: + J_inter X1X2 + K on (0,1), (2,3)) at Delta = 0 every term has even k_X, so [H, X^N] = 0 exactly and
    S L S^-1 = -L^dagger - 2 sigma entry for entry, at dyadic and at generic couplings alike (the global flip keeps
    every summation order). Control: Delta = 1/8 breaks both, exactly nonzero, and the spectrum no longer pairs.
C2  The turn row, Ad_{X^N} = (-1)^{k_X}: at Delta = 0 it commutes with L exactly (every term has even k_X); it is a
    copy and no reflection (conjugation negates no dissipator). This is the field section's first reading, the X^N
    reversal of a field along the wire, at zero field. Control: Delta = 1/8 breaks it.
C3  The half-turn row read on W8 (ii) of PROTON_WIRE_CROSSING: on a mirror-symmetric TFI (J = [j0, j1, j0],
    K = [k0, k0]) with a mirror-ODD bias [d, 0, -d] and no field, S = Rev X^N (Rev the site reversal) commutes with L,
    a copy, while the far kernel {U : [H, U] = 0, U anticommutes with every Z_l} is EMPTY and the spectrum does not
    pair, so no reflection of any kind exists. Both facts are exact in rationals (sympy) at generic rational couplings;
    in floats the copy's residual is a case-3 rounding residual (CLAUDE.md, no rounding): it sits on H's diagonal
    only, is the site-order summation of the bias and ZZ terms that Rev reverses, is a few eps under the site-by-site
    order and exactly 0.0 at every draw when each term is summed with its Rev image first; the order is the
    input the value does not contain, and the gate varies it. Controls: the mirror-EVEN bias [d, 0, d] breaks this copy
    (there Rev alone is one; under the odd bias Rev alone is not, the X^N is needed); at d = 0 X^N is lit, commutes, and the spectrum pairs.

THE SECOND AXIS (which letter the environment reads; the sections above use Z only):
C4  The letter swap. If the environment reads X (the delocalisation |L> +- |R>, the barrier) instead of Z (the
    position), the lit strings are {Y, Z}^N and Z^N is the one tested: the roles of tunnelling and bias swap
    exactly. Bias + ZZ with no tunnelling is palindromic under X reading (every term has even k_Z; the global Z flip
    keeps every summation order, so exact at generic floats) and the tunnelling breaks it; under Z reading it is
    the other way round. Both readings are printed for every wire, so the doc's table comes from one wire.
C5  The doublet discriminator, one proton, H = -J X + Delta Z. At Delta = 0 under X reading the doublet
    populations are dark (L(X) = 0 exactly) and the doublet coherence pays exactly 2 gamma; under Z reading the
    population difference X is an eigenvector at exactly -2 gamma (T1 = 1/(2 gamma)), while the coherence is the
    (Z, Y) block [[0, 2J], [-2J, -2 gamma]], pinned entry by entry, whose eigenvalues -gamma +- sqrt(gamma^2 - 4 J^2)
    give T2 = 1/gamma = 2 T1 only while gamma < 2J (overdamped above: slow rate gamma - sqrt(gamma^2 - 4 J^2)).
    Control: at Delta != 0 the upper level's population rate Tr(P L(P)) is nonzero under X reading. With both
    letters read the Pauli rates are X: 2 gamma_Z, Z: 2 gamma_X, Y: 2 gamma_X + 2 gamma_Z (the cube's rate
    2(gamma_Z k_Z + gamma_X k_X), THE_ONE_SQUARE section 7), so the coherence pair has real part -(2 gamma_X +
    gamma_Z) exactly while gamma_Z < 2J, at any gamma_X, and T2/T1 = 2 gamma_Z/(2 gamma_X + gamma_Z), a dial
    from 2 (pure Z) to 0 (pure X); gated against the eigensolver with gamma_Z > 2J as the control that must
    fail. The Y rate at generic rates is the sum of two dissipators and its float residual is a summation-order
    reading. The book is the unital Hermitian-jump one, D[rho] = gamma (P rho P - rho), which fixes no
    temperature; a thermal bath is F137's channel and outside these rows.
C6  The two-axis bath. With both Z and X jumps on every site the only string anticommuting with every jump is Y^N;
    on the wire the tunnelling and the bias have odd k_Y, so no lit string commutes with H and by F158 the
    spectrum does not pair, DEPOLARIZING_PALINDROME's "a field along either noise axis breaks the two-axis half"
    on the wire; ZZ alone (k_Y = 2) keeps Y^N. The pairing distance at a ten percent X admixture is read.

Exact rows compare to 0: the reflection and the turn in floats at dyadic and at generic couplings, the copy under the
odd bias in rationals at a generic point and in floats under the mirror-paired summation order. The eigensolver's pairing distance
is read beside them with its error model (eps times the spectral scale) and gated only on the controls, where it
must exceed that model by many decades.
"""
import sys
from itertools import product
import numpy as np
import sympy as sp
from scipy.optimize import linear_sum_assignment
sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], dtype=complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
LET = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
def kron_all(ms):
    r = np.eye(1, dtype=complex)
    for m in ms: r = np.kron(r, m)
    return r
def op(N, l, P): return kron_all([P if j == l else I2 for j in range(N)])
def string(s): return kron_all([LET[c] for c in s])
def liou(N, H, gam):
    """Row-stack Liouvillian of -i[H, rho] + gamma sum_l (Z_l rho Z_l - rho)."""
    d = 2 ** N; Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = op(N, l, Z); L += gam * (np.kron(Zl, Zl.T) - np.eye(d * d))
    return L
def right_mult(N, F): return np.kron(np.eye(2 ** N), F.T)
def conj_super(N, U): return np.kron(U, U.conj())
def maxabs(A): return float(np.max(np.abs(A)))
GAM = 0.5
EPS = np.finfo(float).eps
FAILS = []
def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}  {detail}")
    if not ok: FAILS.append(name)

def pairing_distance(L, sigma):
    """Maximum cost of the min-sum assignment between the spectrum and its mirror image about -sigma.
    If a pairing within eps existed the min-sum assignment would have every cost at most dim * eps."""
    ev = np.linalg.eigvals(L); mirrored = -2 * sigma - ev
    cost = np.abs(ev[:, None] - mirrored[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())
def scale(L): return float(np.max(np.abs(np.linalg.eigvals(L))))

# ---- exact (sympy) primitives, N = 3 -------------------------------------------------------------------------
sX = sp.Matrix([[0, 1], [1, 0]]); sZ = sp.diag(1, -1); sI = sp.eye(2)
def skron(ms):
    r = sp.Matrix([[1]])
    for m in ms: r = sp.kronecker_product(r, m)
    return r
def sop(N, l, P): return skron([P if j == l else sI for j in range(N)])
def sliou(N, H, gam):
    d = 2 ** N; Id = sp.eye(d)
    L = -sp.I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
    for l in range(N):
        Zl = sop(N, l, sZ); L += gam * (sp.kronecker_product(Zl, Zl.T) - sp.eye(d * d))
    return L
def far_kernel_dim(N, H):
    """dim {U : [H, U] = 0 and U Z_l + Z_l U = 0 for every l}, exact: U runs over the X/Y strings (the ones that
    anticommute with every Z_l, 2^N of them) and only [H, U] = 0 is imposed."""
    sY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    cols = [skron([{'X': sX, 'Y': sY}[c] for c in w]) for w in product('XY', repeat=N)]
    M = sp.Matrix.hstack(*[(H * c - c * H).reshape(4 ** N, 1) for c in cols])
    return len(cols) - M.rank()

print("=== C1: the hydrogen-bond-qubit models are exactly palindromic at Delta = 0 (row 1, the lit shift X^N) ===")
def hbq(N, delta, J=1.0, K=0.5, Jint=0.25):
    H = sum(-J * op(N, l, X) + delta * op(N, l, Z) for l in range(N))
    if N == 2: H = H + K * op(N, 0, Z) @ op(N, 1, Z)
    if N == 4: H = H + K * (op(N, 0, Z) @ op(N, 1, Z) + op(N, 2, Z) @ op(N, 3, Z)) + Jint * op(N, 1, X) @ op(N, 2, X)
    return H
for N in (2, 4):
    XN = string('X' * N); S = right_mult(N, XN); sigma = N * GAM
    for delta in (0.0, 0.125):
        H = hbq(N, delta); L = liou(N, H, GAM)
        comm = maxabs(H @ XN - XN @ H)
        res = maxabs(S @ L @ np.linalg.inv(S) + L.conj().T + 2 * sigma * np.eye(L.shape[0]))
        pd = pairing_distance(L, sigma); model = EPS * scale(L)
        if delta == 0.0:
            check(f"C1 N={N} Delta=0: [H, X^N] == 0.0", comm == 0.0, f"max {comm}")
            check(f"C1 N={N} Delta=0: S L S^-1 + L^dagger + 2 sigma == 0.0", res == 0.0, f"max {res}")
            print(f"       pairing distance {pd:.2e} = {pd / model:.1f} x (eps * spectral scale), read not gated")
        else:
            check(f"C1 control N={N} Delta=1/8: [H, X^N] != 0", comm != 0.0, f"max {comm}")
            check(f"C1 control N={N} Delta=1/8: reflection residual != 0", res != 0.0, f"max {res}")
            check(f"C1 control N={N} Delta=1/8: spectrum does not pair", pd > 1e6 * model, f"{pd:.2e} = {pd / model:.1e} x model")
    # generic (non-dyadic) couplings: the global flip keeps every summation order, so the identity stays exactly 0.0
    rng = np.random.default_rng(7); worst = 0.0
    for _ in range(50):
        J, K, Jint, g = rng.uniform(0.1, 2.0, 4)
        H = hbq(N, 0.0, J, K, Jint); L = liou(N, H, g)
        worst = max(worst, maxabs(right_mult(N, XN) @ L @ np.linalg.inv(right_mult(N, XN)) + L.conj().T + 2 * N * g * np.eye(L.shape[0])))
    check(f"C1 N={N} Delta=0, 50 generic (J, K, J_inter, gamma): reflection residual == 0.0 at every draw", worst == 0.0, f"worst {worst}")

print("\n=== C2: the turn Ad_{X^N} at Delta = 0 is a copy (the field section's first reading at zero field) ===")
for N in (2, 4):
    A = conj_super(N, string('X' * N))
    for delta, expect in ((0.0, True), (0.125, False)):
        L = liou(N, hbq(N, delta), GAM); comm = maxabs(L @ A - A @ L)
        check(f"C2 {'control ' if delta else ''}N={N} Delta={delta}: [L, Ad_X^N] {'==' if expect else '!='} 0.0", (comm == 0.0) == expect, f"max {comm}")

print("\n=== C3: W8 (ii), the mirror-odd bias: S = Rev X^N is a copy and no reflection exists ===")
N = 3
def tfi(J, K, dlt, Xm=X, Zm=Z, opf=op):
    H = sum(-J[l] * opf(N, l, Xm) + dlt[l] * opf(N, l, Zm) for l in range(N))
    return H + K[0] * opf(N, 0, Zm) @ opf(N, 1, Zm) + K[1] * opf(N, 1, Zm) @ opf(N, 2, Zm)
def stfi(J, K, dlt):
    H = sum((-J[l] * sop(N, l, sX) + dlt[l] * sop(N, l, sZ) for l in range(N)), sp.zeros(8, 8))
    return H + K[0] * sop(N, 0, sZ) * sop(N, 1, sZ) + K[1] * sop(N, 1, sZ) * sop(N, 2, sZ)
Rev = np.zeros((8, 8))
for x in range(8):
    b = [(x >> 2) & 1, (x >> 1) & 1, x & 1]; y = (b[2] << 2) | (b[1] << 1) | b[0]; Rev[y, x] = 1
Smat = Rev @ string('XXX'); sRev = sp.Matrix(Rev.astype(int)); sS = sRev * skron([sX] * 3)
jumps = [op(N, l, Z) for l in range(N)]
def lit_commuting_strings(H):
    return [s for s in (''.join(t) for t in product('IXYZ', repeat=N))
            if maxabs(string(s) @ H - H @ string(s)) == 0.0
            and all(maxabs(string(s) @ Zl + Zl @ string(s)) == 0.0 for Zl in jumps)]
for label, dlt, copy_expected, lit_expected, pairs in (
        ("mirror-odd bias [d,0,-d]", [0.25, 0, -0.25], True, [], False),
        ("control: mirror-even bias [d,0,d]", [0.25, 0, 0.25], False, [], False),
        ("control: no bias", [0, 0, 0], True, ['XXX'], True)):
    H = tfi([1.0, 0.5, 1.0], [0.25, 0.25], dlt); L = liou(N, H, GAM); A = conj_super(N, Smat)
    comm = maxabs(L @ A - A @ L); lit = lit_commuting_strings(H)
    pd = pairing_distance(L, N * GAM); model = EPS * scale(L)
    check(f"C3 {label} (dyadic): [L, Ad_S] {'==' if copy_expected else '!='} 0.0", (comm == 0.0) == copy_expected, f"max {comm}")
    check(f"C3 {label}: lit strings commuting with H = {lit_expected}", lit == lit_expected, f"found {lit or 'none'}")
    if pairs:
        print(f"       pairing distance {pd:.2e} = {pd / model:.1f} x (eps * spectral scale), read not gated")
    else:
        check(f"C3 {label}: spectrum does not pair (no reflection of any kind)", pd > 1e6 * model, f"{pd:.2e} = {pd / model:.1e} x model")

# exact in rationals at generic couplings: the copy and the empty far kernel
R = sp.Rational
J_r, K_r = [R(7, 10), R(13, 10), R(7, 10)], [R(3, 10), R(3, 10)]
for label, dlt, copy_expected, rev_expected, dim_expected in (
        ("mirror-odd bias", [R(1, 7), 0, R(-1, 7)], True, False, 0),
        ("control: mirror-even bias", [R(1, 7), 0, R(1, 7)], False, True, 0),
        ("control: no bias", [0, 0, 0], True, True, 1)):
    Hs = stfi(J_r, K_r, dlt); Ls = sliou(N, Hs, R(2, 9)); As = sp.kronecker_product(sS, sS)
    comm_exact = (Ls * As - As * Ls).is_zero_matrix
    check(f"C3 {label} (rational, generic): [L, Ad_S] = 0 exactly is {copy_expected}", comm_exact == copy_expected)
    rev_only = (sp.kronecker_product(sRev, sRev) * Ls - Ls * sp.kronecker_product(sRev, sRev)).is_zero_matrix
    check(f"C3 {label} (rational): the site reversal Rev alone is a copy is {rev_expected} (Rev sends the bias to its mirror image)", rev_only == rev_expected)
    dim = far_kernel_dim(N, Hs)
    check(f"C3 {label} (rational): far kernel dimension = {dim_expected}", dim == dim_expected, f"dim {dim}")

# the float residual of the copy at generic couplings is a function of the SUMMATION ORDER of H's terms and of
# nothing else in the physics (CLAUDE.md, no rounding, case 3): site by site it is a few eps on H's diagonal; with
# the mirror pairs summed together (Z_0 + Z_2 before scaling, Z_0 Z_1 + Z_1 Z_2 before scaling) it is 0.0 at every
# draw. The 0.0 setting is the ORDER, not a dyadic coupling.
def tfi_paired(J, K, dlt):
    """Same H, mirror-paired summation: every term and its image under Rev are added before anything else."""
    H = -J[1] * op(N, 1, X) - J[0] * (op(N, 0, X) + op(N, 2, X)) + dlt[1] * op(N, 1, Z)
    H = H + (dlt[0] * op(N, 0, Z) + dlt[2] * op(N, 2, Z))
    return H + K[0] * (op(N, 0, Z) @ op(N, 1, Z) + op(N, 1, Z) @ op(N, 2, Z))
rng = np.random.default_rng(1); nonzero = 0; worst = 0.0; offdiag_worst = 0.0; worst_paired = 0.0; same_H = 0.0
A = conj_super(N, Smat)
for _ in range(200):
    j0, j1, k0, d, g = rng.uniform(0.1, 2.0, 5)
    H = tfi([j0, j1, j0], [k0, k0], [d, 0, -d]); L = liou(N, H, g)
    r = maxabs(L @ A - A @ L); nonzero += r != 0.0; worst = max(worst, r)
    Hres = Smat @ H @ Smat.T - H; offdiag_worst = max(offdiag_worst, maxabs(Hres - np.diag(np.diag(Hres))))
    Hp = tfi_paired([j0, j1, j0], [k0, k0], [d, 0, -d]); Lp = liou(N, Hp, g)
    worst_paired = max(worst_paired, maxabs(Lp @ A - A @ Lp)); same_H = max(same_H, maxabs(Hp - H))
print(f"       200 generic float draws of the mirror-odd copy, site-by-site summation: residual nonzero in {nonzero}, "
      f"worst {worst:.2e} = {worst / EPS:.1f} eps, off-diagonal part of Rev H Rev - H exactly {offdiag_worst} (read)")
check("C3 generic floats: the copy residual is a function of the summation order: mirror-paired summation gives exactly 0.0 at all 200 draws",
      worst_paired == 0.0, f"worst {worst_paired}")
check("C3 generic floats: the site-by-site summation is NOT exactly 0.0 at every draw (the residual exists and the order is what moves it)",
      nonzero > 0, f"nonzero in {nonzero} of 200; the two H differ by at most {same_H:.1e}")

def liou_multi(N, H, rates):
    """Row-stack Liouvillian with one dephasing letter per (letter, rate) pair on every site."""
    d = 2 ** N; L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    for letter, g in rates:
        for l in range(N):
            P = op(N, l, letter); L += g * (np.kron(P, P.T) - np.eye(d * d))
    return L

print("\n=== C4: the letter swap: under X reading the bias keeps the palindrome and the tunnelling breaks it ===")
print("       wire J = 1 on every bond where present, K = 0.4, bias [0.3, 0.1, 0.2] where present, gamma = 0.5")
N = 3; ZN = string('Z' * N); SZ = right_mult(N, ZN)
def wire(J, K, dlt): return tfi([J, J, J], [K, K], dlt)
bias = [0.3, 0.1, 0.2]
for name, H, zn_commutes, xn_commutes in (("bias + ZZ, no tunnelling", wire(0.0, 0.4, bias), True, False),
                                          ("control: tunnelling + bias + ZZ", wire(1.0, 0.4, bias), False, False),
                                          ("control: tunnelling + ZZ, no bias", wire(1.0, 0.4, [0, 0, 0]), False, True)):
    comm = maxabs(H @ ZN - ZN @ H)
    LX = liou_multi(N, H, [(X, GAM)]); LZ = liou_multi(N, H, [(Z, GAM)]); sigma = N * GAM
    resX = maxabs(SZ @ LX @ np.linalg.inv(SZ) + LX.conj().T + 2 * sigma * np.eye(64))
    pdX, pdZ = pairing_distance(LX, sigma), pairing_distance(LZ, sigma); modelX, modelZ = EPS * scale(LX), EPS * scale(LZ)
    check(f"C4 {name}: [H, Z^N] {'==' if zn_commutes else '!='} 0.0", (comm == 0.0) == zn_commutes, f"max {comm}")
    check(f"C4 {name}: X reading, shift by Z^N is a reflection: residual {'==' if zn_commutes else '!='} 0.0", (resX == 0.0) == zn_commutes, f"max {resX}")
    print(f"       pairing distance, Z reading {pdZ:.2e} = {pdZ / modelZ:.1e} x model, X reading {pdX:.2e} = {pdX / modelX:.1e} x model")
    if zn_commutes: check(f"C4 {name}: Z reading does not pair (the bias breaks the Z palindrome)", pdZ > 1e6 * modelZ)
    else: check(f"C4 {name}: X reading does not pair", pdX > 1e6 * modelX)
    if xn_commutes: check(f"C4 {name}: Z reading pairs to rounding (read as a ratio to the model)", pdZ < 1e3 * modelZ)
rng = np.random.default_rng(3); worst = 0.0
for _ in range(50):
    J, K, g = 0.0, rng.uniform(0.1, 2.0), rng.uniform(0.1, 2.0); dl = list(rng.uniform(-1, 1, 3))
    H = wire(J, K, dl); LX = liou_multi(N, H, [(X, g)])
    worst = max(worst, maxabs(SZ @ LX @ np.linalg.inv(SZ) + LX.conj().T + 2 * N * g * np.eye(64)))
check("C4 bias + ZZ under X reading, 50 generic (K, bias profile, gamma): reflection residual == 0.0 at every draw", worst == 0.0, f"worst {worst}")

print("\n=== C5: the doublet discriminator, one proton, H = -J X + Delta Z ===")
J = 1.0; H1 = -J * X
plus = np.array([1, 1]) / np.sqrt(2); minus = np.array([1, -1]) / np.sqrt(2)
def apply(L, M): return (L @ M.reshape(-1)).reshape(2, 2)
def coef(M, P): return np.trace(M @ P) / 2          # Pauli coefficient
for g in (0.1, 0.7):
    LZ1 = liou_multi(1, H1, [(Z, g)]); LX1 = liou_multi(1, H1, [(X, g)])
    check(f"C5 gamma={g} X reading: L(X) == 0 exactly (the doublet populations are dark)", maxabs(apply(LX1, X)) == 0.0)
    check(f"C5 gamma={g} X reading: L(Z) == -2 gamma Z + 2J Y exactly and L(Y) == -2 gamma Y - 2J Z exactly (the coherence pays 2 gamma)",
          maxabs(apply(LX1, Z) + 2 * g * Z - 2 * J * Y) == 0.0 and maxabs(apply(LX1, Y) + 2 * g * Y + 2 * J * Z) == 0.0)
    check(f"C5 gamma={g} Z reading: L(X) == -2 gamma X exactly (the doublet relaxes at 2 gamma, T1 = 1/(2 gamma))", maxabs(apply(LZ1, X) + 2 * g * X) == 0.0)
    check(f"C5 gamma={g} Z reading: the (Z, Y) block is [[0, 2J], [-2J, -2 gamma]] exactly (L(Z) = 2J Y, L(Y) = -2J Z - 2 gamma Y)",
          maxabs(apply(LZ1, Z) - 2 * J * Y) == 0.0 and maxabs(apply(LZ1, Y) + 2 * J * Z + 2 * g * Y) == 0.0)
# the block's eigenvalues -gamma +- sqrt(gamma^2 - 4J^2) follow by hand; read them against the eigensolver in both regimes
for g in (0.1, 5.0):
    rates = np.sort(-np.linalg.eigvals(liou_multi(1, H1, [(Z, g)])).real)
    pred = np.sort(np.real(np.array([0, 2 * g, g - np.sqrt(complex(g * g - 4 * J * J)), g + np.sqrt(complex(g * g - 4 * J * J))])))
    print(f"       Z reading gamma={g} (J=1, {'under' if g < 2 * J else 'over'}damped): rates {np.round(rates, 4)}, by hand {np.round(pred, 4)}; "
          f"slowest coherence rate {rates[1]:.4f} vs gamma = {g} (read)")
# control: Delta != 0, the upper level's population rate is nonzero under X reading
Hd = -J * X + 0.5 * Z; w, v = np.linalg.eigh(Hd); up = v[:, 1]; Pup = np.outer(up, up.conj())
LXd = liou_multi(1, Hd, [(X, 0.1)]); rate = float(np.real(np.trace(Pup @ apply(LXd, Pup))))
check("C5 control Delta=0.5, X reading: population rate Tr(P L(P)) of the upper level != 0 (no longer dark)", rate != 0.0, f"rate {rate:.4f}")
check("C5 Delta=0, X reading: L(upper-level projector) == 0.0 exactly", maxabs(apply(liou_multi(1, H1, [(X, 0.1)]), np.outer(minus, minus))) == 0.0)
# both letters read: the dial, generic rates
rng = np.random.default_rng(11)
for gz, gx in ((0.1, 0.0), (0.1, 0.05), (0.1, 0.1), tuple(rng.uniform(0.05, 0.9, 2)), tuple(rng.uniform(0.05, 0.9, 2))):
    Lm = liou_multi(1, H1, [(Z, gz), (X, gx)])
    xres = maxabs(apply(Lm, X) + 2 * gz * X); zres = abs(coef(apply(Lm, Z), Z) + 2 * gx); yres = abs(coef(apply(Lm, Y), Y) + 2 * gx + 2 * gz)
    if gx == 0.0:
        check(f"C5 both letters gamma_Z={gz:.4f}, gamma_X=0: L(X) == -2 gamma_Z X exactly (one dissipator, no sum)", xres == 0.0, f"max {xres}")
    else:
        print(f"       gamma_Z={gz:.4f}, gamma_X={gx:.4f}: Pauli-rate residuals X {xres:.1e}, Z {zres:.1e}, Y {yres:.1e} against 2 gamma_Z, 2 gamma_X, "
              f"2 gamma_X + 2 gamma_Z: the two dissipators' diagonals are summed, a summation-order reading")
    ev = np.linalg.eigvals(Lm); coh = sorted(ev, key=lambda z: abs(z.imag))[-2:]      # the rotating pair
    re_res = max(abs(z.real + 2 * gx + gz) for z in coh); model = EPS * scale(Lm)
    check(f"C5 both letters gamma_Z={gz:.4f}, gamma_X={gx:.4f}: coherence pair real part == -(2 gamma_X + gamma_Z) to the eigensolver (gamma_Z < 2J)",
          re_res < 1e3 * model, f"residual {re_res:.1e} = {re_res / model:.1f} x model")
    t2_over_t1 = (2 * gz) / (-np.mean([z.real for z in coh]))      # T2 = 1/(-Re lambda), T1 = 1/(2 gamma_Z)
    print(f"       T2/T1 from L = {t2_over_t1:.3f}, formula 2 gamma_Z/(2 gamma_X + gamma_Z) = {2 * gz / (2 * gx + gz):.3f} (read)")
gz, gx = 3.0, 0.2
ev = np.linalg.eigvals(liou_multi(1, H1, [(Z, gz), (X, gx)])); reals = sorted(set(np.round(ev.real, 6)))
check("C5 control gamma_Z=3 > 2J: the coherence pair is overdamped, two distinct real rates, neither -(2 gamma_X + gamma_Z)",
      len(reals) == 4 and all(abs(r + 2 * gx + gz) > 1e-3 for r in reals), f"rates {[-r for r in reals]}")
from scipy.linalg import expm
g = 0.1; t = 5.0
for letter, name, pop_form, coh_form in ((Z, "Z (position)", f"(1 + e^(-2 gamma t))/2 = {(1 + np.exp(-2 * g * t)) / 2:.6f}", "a damped rotation, no single exponential"),
                                         (X, "X (barrier)", "1 exactly", f"e^(-2 gamma t)/2 = {np.exp(-2 * g * t) / 2:.6f}")):
    L1 = liou_multi(1, H1, [(letter, g)]); U = expm(L1 * t)
    rho = (U @ np.outer(plus, plus).reshape(-1)).reshape(2, 2); pop = float(np.real(plus @ rho @ plus))
    rho0 = np.outer(plus, plus) / 2 + np.outer(minus, minus) / 2 + np.outer(plus, minus) / 2 + np.outer(minus, plus) / 2
    coh = abs(plus @ (U @ rho0.reshape(-1)).reshape(2, 2) @ minus)
    print(f"       gamma t = {g * t}: bath reads {name:13s}: doublet population 1 -> {pop:.6f} ({pop_form}), doublet coherence 0.5 -> {coh:.6f} ({coh_form}), read")

print("\n=== C6: the mixed bath on the wire: both letters read, no lit element commutes with the tunnelling ===")
N = 3
jumps_zx = [op(N, l, P) for l in range(N) for P in (Z, X)]
anti_all = [w for w in (''.join(t) for t in product('IXYZ', repeat=N)) if all(maxabs(string(w) @ Q + Q @ string(w)) == 0.0 for Q in jumps_zx)]
check("C6 strings anticommuting with every Z_l and every X_l = ['YYY'] only (the lit span is one-dimensional, so [H, Y^N] != 0 means no colouring and, by F158, no pairing)", anti_all == ['YYY'], f"found {anti_all}")
YN = string('YYY')
for name, H, expect in (("tunnelling + ZZ", wire(1.0, 0.4, [0, 0, 0]), False), ("bias + ZZ", wire(0.0, 0.4, bias), False), ("ZZ only", wire(0.0, 0.4, [0, 0, 0]), True)):
    comm = maxabs(H @ YN - YN @ H)
    check(f"C6 {name}: [H, Y^N] {'==' if expect else '!='} 0.0", (comm == 0.0) == expect, f"max {comm}")
H = wire(1.0, 0.4, [0, 0, 0])
for gx, name in ((0.5, "Z and X equal"), (0.05, "Z with ten percent X")):
    L = liou_multi(N, H, [(Z, 0.5), (X, gx)]); sigma = N * (0.5 + gx); pd = pairing_distance(L, sigma); model = EPS * scale(L)
    if gx == 0.5:
        check(f"C6 tunnelling + ZZ, bath {name}: spectrum does not pair", pd > 1e6 * model, f"{pd:.2e} = {pd / model:.1e} x model")
    else:
        print(f"       tunnelling + ZZ, bath {name}: pairing distance {pd:.2e} = {pd / model:.1e} x model (read, not gated)")
L = liou_multi(N, wire(0.0, 0.4, [0, 0, 0]), [(Z, 0.5), (X, 0.5)]); pd = pairing_distance(L, 3.0)
print(f"       ZZ only, bath Z and X equal: pairing distance {pd:.2e} = {pd / (EPS * scale(L)):.1f} x model (Y^N lit and commuting, read)")

print(f"\n{'ALL OK' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
sys.exit(1 if FAILS else 0)
