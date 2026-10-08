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

print(f"\n{'ALL OK' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
sys.exit(1 if FAILS else 0)
