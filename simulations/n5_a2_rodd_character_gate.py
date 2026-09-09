"""Gate: the N=5 positive-real R-odd A2 locus is SEMISIMPLE, and the read that says so.

WHAT THIS IS BESIDE.  The verdict itself is not new. It is already gated in C#, by
compute/RCPsiSquared.Diagnostics.Tests/Foundation/RouteBA2PositiveAnchorTests.cs,
at these very numbers (WAnchor = 5.10083102, LambdaAnchor = -4.79196037, both q
lifts, R-odd), asserting Diabolic with alg = geo = 2. That gate is entirely
floating point: distances[1] < 1e-4, relative departure < 1e-6, dephase departure
< 1e-6. docs/CAUGHT_ERRORS.md (2026-09-06) caught an inference at this locus and
closed with the sentence that exact rank must determine its geometric
multiplicity. The companion proof supplies an exact argument in place of that
reading; this file is its numerical corroboration and nothing more. It does not
fill simulations/results/route_b_a2_n5.json's empty exactRankCertificates list.

THE ARGUMENT IT CORROBORATES (proved in docs/proofs/PROOF_N5_REAL_Q_DIABOLIC.md):
A2_O is irreducible over Q of degree 13 with 6 negative-real w roots, both
EXECUTED over Z in simulations/o2b_gcd_certificate.py (irreducibility at line 1336,
the A2 root inventory at lines 1342-1343).  In the R-odd sector, L_O = D_O + t*K_O
with D_O and K_O RATIONAL and K_O SYMMETRIC, and w = qUnitHop^2 = -t^2 with
t = i*qUnitHop.  At the positive-real root w0 the value t0 is nonzero purely
imaginary while Q(w0) is a subfield of R, so [Q(t0):Q] = 26 and all 26 t values
form a single Galois orbit.  That orbit contains the 12 REAL t coming from the 6
negative w, where L_O is real symmetric and geometric multiplicity equals
algebraic.  Both quantities are Galois invariant, so they agree at t0 as well.

WHAT EACH STEP IS SENSITIVE TO.

  STEP 0  ANCHOR       the inventory shape (13 = 6 + 1 + 6) and the q/J book.
                       The book is pinned twice, and the two catch different
                       mutations: a wrong BOOKS family dies here (no scaling
                       matches), while a wrong QPHYS halving survives here, because
                       the search simply relabels which scaling wins, and dies in
                       STEP 2 where all six source roots go red.
  STEP 1  STRUCTURE    Hermiticity at real t, entry-wise == 0.0, and the exact
                       R-invariance of the sector in the integer orbit basis,
                       also == 0.0.  Plus a CONSTRUCTION check on the commutator
                       sign, which is here because no spectral check can do it:
                       flipping the bra-side sign conjugates the whole R-odd
                       spectrum, and the A2 inventory is conjugation-closed, so
                       every anchor in it reads identically either way.
  STEP 2  SOURCE       nullity 2 at all six negative-real w roots (real t).  This
                       is the end of the Galois orbit the transfer carries from,
                       and it is the step that owns the factor-2 convention.
  STEP 3  CONCLUSION   nullity 2 at w0, against a defective control (must read 1)
                       and a generic-lambda control (must read 0), all three on
                       the same code path in the same sector.
  STEP 4  ERROR MODEL  the law behind the thresholds, measured rather than
                       asserted: at w0 the second singular value is proportional
                       to the coupling detune, sigma2 = c * (dJ/J) with c flat
                       across eight decades, while at the defective point sigma2 is
                       pinned at the coupling scale and does not move.  The
                       failure is one-sided AT THE DEFECTIVE POINT MEASURED, whose
                       Jordan coupling is 1.95e-2, eight orders above NULL_TOL.  That
                       is a measurement, not a bound: nothing here bounds the Jordan
                       coupling from below, so a defective point with a weak enough
                       coupling would read nullity 2 and be called semisimple.  The
                       field argument does not depend on this reading.

Run: python simulations/n5_a2_rodd_character_gate.py
"""
import json
import sys

import numpy as np
from scipy.optimize import minimize_scalar

# --- block builder, copied from simulations/eta_ladder_blocks.py (never import a producer) ---


def popcount(x):
    return bin(x).count("1")


def bit(l, N):
    return 1 << (N - 1 - l)


def block_basis(N, p, q):
    A = [a for a in range(1 << N) if popcount(a) == p]
    B = [b for b in range(1 << N) if popcount(b) == q]
    pairs = [(a, b) for a in A for b in B]
    return pairs, {ab: i for i, ab in enumerate(pairs)}


def hamiltonian_hops(a, N, J):
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
        for a2, amp in hamiltonian_hops(a, N, J):
            L[index[(a2, b)], col] += -1j * amp
        for b2, amp in hamiltonian_hops(b, N, J):
            L[index[(a, b2)], col] += +1j * amp
    return L, pairs, index


def reflect(x, N):
    y = 0
    for l in range(N):
        if x & bit(l, N):
            y |= bit(N - 1 - l, N)
    return y


def parity_columns(N, p, q, integer_basis=False):
    """Columns spanning the R-even and R-odd sectors of the (p,q) block.

    integer_basis=True returns the +-1 orbit basis (integer entries, so the
    invariance residual below is an exact 0.0 rather than a rounded one)."""
    pairs, index = block_basis(N, p, q)
    n = len(pairs)
    scale = 1.0 if integer_basis else 1 / np.sqrt(2)
    seen, cols_e, cols_o = set(), [], []
    for i, (a, b) in enumerate(pairs):
        if i in seen:
            continue
        j = index[(reflect(a, N), reflect(b, N))]
        seen.add(i)
        seen.add(j)
        v = np.zeros(n)
        if i == j:
            v[i] = 1.0
            cols_e.append(v)
        else:
            v[i] = v[j] = scale
            cols_e.append(v)
            u = np.zeros(n)
            u[i] = scale
            u[j] = -scale
            cols_o.append(u)
    return np.array(cols_e).T, np.array(cols_o).T


N, P, Q = 5, 1, 2
GAMMA = [1.0] * N
UE, UO = parity_columns(N, P, Q)
UE_INT, UO_INT = parity_columns(N, P, Q, integer_basis=True)

INVENTORY = "simulations/results/route_b_a2_n5.json"
LAM0 = -4.7919603651796410238703898211366803970127689805559
DEFECTIVE_Q_PHYS = 2.804888          # experiments/F89_PATH_K_DIABOLIC.md:146, six digits
DEFECTIVE_LAM_SEED = -4.4882
# The error model these two thresholds are read against is STEP 4's measured coupling
# law, sigma2 = c*(dJ/J) with c flat over eight decades, not a second model stated
# here. Inverting it turns any sigma2 into an implied coupling error in ULP of J,
# which is the unit both thresholds are reported in at the end of the run.
NULL_TOL = 1e-10                            # a null read must land far BELOW this
CLEAR_TOL = 1e-3                            # a clear read must land far ABOVE it


def sectors(j_scalar, bonds=None):
    b = bonds if bonds is not None else [j_scalar] * (N - 1)
    L, _, _ = build_L(N, P, Q, b, GAMMA)
    return L, UE.T @ L @ UE, UO.T @ L @ UO


def sigma3(M, lam):
    """The three smallest singular values of (M - lam*I), scaled by ||M||_2."""
    s = np.linalg.svd(M - lam * np.eye(M.shape[0]), compute_uv=False)
    return np.sort(s)[:3] / np.linalg.norm(M, 2)


def nullity(M, lam):
    s = sigma3(M, lam)
    return int(s[0] < NULL_TOL) + int(s[1] < NULL_TOL), s


FAIL = []


def check(name, ok, detail):
    print(("  PASS  " if ok else "  FAIL  ") + name + "   " + detail)
    if not ok:
        FAIL.append(name)


inv = json.load(open(INVENTORY))
odd = [s for s in inv["sectors"] if s["parity"] == "O"][0]
a2 = odd["a2Roots"]
pos = [r for r in a2 if r["rootKind"] == "positiveReal"]
negs = [r for r in a2 if r["rootKind"] == "negativeReal"]
nonreal = [r for r in a2 if r["rootKind"] == "nonreal"]

# ---------------------------------------------------------------- STEP 0
print("STEP 0  ANCHOR: the inventory shape the proof needs, and the q/J book")

check("block dims are (full 50, R-even 26, R-odd 24)",
      (UE.shape[0], UE.shape[1], UO.shape[1]) == (50, 26, 24),
      f"got ({UE.shape[0]}, {UE.shape[1]}, {UO.shape[1]})")
# the degree and the split both carry the proof: degree 13 with a negative root is
# what puts a real t in the orbit of the imaginary t0.
check("A2_O has degree 13, split 6 negative / 1 positive / 6 nonreal",
      (len(a2), len(negs), len(pos), len(nonreal)) == (13, 6, 1, 6),
      f"got total {len(a2)} = {len(negs)} + {len(pos)} + {len(nonreal)}")
if FAIL:
    print("\nINVENTORY SHAPE FAILED; nothing below is trustworthy.")
    sys.exit(1)

W0 = float(pos[0]["wSeed"]["real"])
QPHYS = np.sqrt(W0) / 2.0            # route_b_a2_n5.json: w = qUnitHop^2, qPhysical = qUnitHop/2
BOOKS = {"J = q_phys": 1.0, "J = 2*q_phys (unit hop)": 2.0, "J = q_phys/2": 0.5}
matched = [label for label, scale in BOOKS.items()
           if np.min(np.abs(np.linalg.eigvals(sectors(scale * QPHYS)[2]) - LAM0)) < 1e-12]
check("exactly one q/J book puts the stored lambda in the R-odd spectrum", len(matched) == 1,
      f"matched={matched}")
if len(matched) != 1:
    print("\nANCHOR FAILED: the convention is not pinned.")
    sys.exit(1)
SCALE = BOOKS[matched[0]]
print("        book LABEL pinned here: " + matched[0] + "; STEP 2 pins the factor-2 convention")

_, Le0, Lo0 = sectors(SCALE * QPHYS)
d_odd = np.min(np.abs(np.linalg.eigvals(Lo0) - LAM0))
d_even = np.min(np.abs(np.linalg.eigvals(Le0) - LAM0))
check("lambda0 is in the R-ODD sector", d_odd < 1e-12, f"|dlambda| = {d_odd:.2e}")
check("and the R-EVEN sector does NOT carry it (the parity label is measured)", d_even > 1e-6,
      f"|dlambda| = {d_even:.2e}")

# ---------------------------------------------------------------- STEP 1
print("\nSTEP 1  STRUCTURE: the three facts the proof reads off the matrix")

MAG = 0.7134912837465
L_realt, _, _ = sectors(1j * MAG)          # t real  <=>  the bond is purely imaginary
res_herm = np.max(np.abs(L_realt - L_realt.conj().T))
check("L is EXACTLY Hermitian at real t (== 0.0)", res_herm == 0.0,
      f"max|L - L^H| = {res_herm!r}")
# the pair-read: the same quantity on the physical axis is 2*MAG by construction, so it
# reports that the hopping is present, not that the Hermiticity check discriminates.
L_phys, _, _ = sectors(MAG)
res_phys = np.max(np.abs(L_phys - L_phys.conj().T))
# Not a threshold: 2*MAG is what this residual equals by construction, so it is compared
# to that value exactly. It reports that the hopping is present and discriminates nothing
# else; the Hermiticity verdict is the read above it.
check("and the physical axis is not Hermitian (a presence read, EXACTLY 2*MAG)",
      res_phys == 2 * MAG,
      f"max|L - L^H| = {res_phys:.4f}")

# The sector split is legitimate only if R commutes with L.  R is a permutation, so both
# sides of R L = L R are the same float values rearranged: the exact route, == 0.0.  A
# least-squares compression would only measure numpy's rounding here.
pairs_r, index_r = block_basis(N, P, Q)
Rperm = np.zeros((len(pairs_r), len(pairs_r)))
for i, (a, b) in enumerate(pairs_r):
    Rperm[index_r[(reflect(a, N), reflect(b, N))], i] = 1.0
L_at_w0, _, _ = sectors(SCALE * QPHYS)
res_inv = np.max(np.abs(Rperm @ L_at_w0 - L_at_w0 @ Rperm))
# Commuting alone is not enough: the identity permutation commutes with everything. Tie
# Rperm to the projectors actually used, so the check is about THIS reflection.
res_ro = np.max(np.abs(Rperm @ UO + UO))
res_re = np.max(np.abs(Rperm @ UE - UE))
check("R commutes with L AND is the reflection the projectors use (all == 0.0)",
      res_inv == 0.0 and res_ro == 0.0 and res_re == 0.0,
      f"max|RL - LR| = {res_inv!r}, max|R U_O + U_O| = {res_ro!r}, max|R U_E - U_E| = {res_re!r}")
# What steps 2 and 3 consume: in the +-1 orbit basis K_O is INTEGER and SYMMETRIC, so
# L_O = D_O + t*K_O is over Q and is real symmetric at real t.  Note what is NOT claimed:
# K_E's asymmetry in this basis is only the orbit-size normalisation (E mixes column norms
# 1 and sqrt 2, O is uniformly sqrt 2), and in an orthonormal basis BOTH restrictions are
# exactly symmetric.  The sector is chosen because A2_O is a per-sector object, not because
# the even sector fails to be symmetric.
K = (sectors(1.0)[0] - sectors(0.0)[0]) / (1j)   # L(J) = D + t*K with t = i*J
# U_O^T U_O = 2*I in the +-1 basis, so the restriction is an exact halving in binary;
# a least-squares solve here would only measure numpy's rounding.
Ko = 0.5 * (UO_INT.T @ K @ UO_INT)
res_ko = np.max(np.abs(Ko - Ko.T))
int_res = np.max(np.abs(Ko - np.round(Ko)))
check("K_O is EXACTLY symmetric and INTEGER in the orbit basis (== 0.0, on a nonzero K_O)",
      res_ko == 0.0 and int_res == 0.0 and np.max(np.abs(Ko)) > 0.5,
      f"max|K_O - K_O^T| = {res_ko!r}, max|K_O - round(K_O)| = {int_res!r}, "
      f"max|K_O| = {np.max(np.abs(Ko))!r}")

# The other half of "L_O = D_O + t*K_O over Q", which steps 2 and 3 consume just as much:
# D_O must restrict to a RATIONAL DIAGONAL on the odd sector, or the block would not be
# real symmetric at real t no matter what K_O does. Exact route, so exact comparison.
Do = 0.5 * (UO_INT.T @ sectors(0.0)[0] @ UO_INT)
off_do = np.max(np.abs(Do - np.diag(np.diag(Do))))
vals_do = sorted(set(np.round(np.diag(Do).real, 12)))
check("D_O is EXACTLY diagonal on the odd sector, with the AT dephasing values (== 0.0)",
      off_do == 0.0 and np.max(np.abs(np.diag(Do).imag)) == 0.0 and vals_do == [-6.0, -2.0],
      f"max|off-diagonal| = {off_do!r}, diagonal values = {vals_do}")

# The second q lift is not a second theorem: L_O(-t) = conj(L_O(t)) entry for entry and
# lambda0 is real, so the two shifted matrices are complex conjugates and share every rank.
# Exact route, read as 0.0; the nullity is then read at BOTH lifts rather than assumed.
Lp_lift, Lm_lift = sectors(SCALE * QPHYS)[2], sectors(-SCALE * QPHYS)[2]
res_lift = np.max(np.abs(Lm_lift - np.conj(Lp_lift)))
n_plus, s_plus = nullity(Lp_lift, LAM0)
n_minus, s_minus = nullity(Lm_lift, LAM0)
check("the SECOND q lift is the conjugate of the first (== 0.0) and reads the same nullity",
      res_lift == 0.0 and n_plus == n_minus == 2,
      f"max|L(-J) - conj(L(J))| = {res_lift!r}; nullity {n_plus} at +q and {n_minus} at -q, "
      f"sigma2 {s_plus[1]:.4e} vs {s_minus[1]:.4e}")

# CONSTRUCTION check, and it is here because no spectral check can do it: flipping the
# bra-side sign conjugates the whole R-odd spectrum, and the A2 inventory is closed under
# conjugation, so every anchor in it reads identically either way.
pairs, index = block_basis(N, P, Q)
Lc, _, _ = sectors(1.0)
ket_sign = bra_sign = None
for col, (a, b) in enumerate(pairs):
    for a2h, _ in hamiltonian_hops(a, N, [1.0] * (N - 1)):
        ket_sign = np.sign(Lc[index[(a2h, b)], col].imag)
        break
    if ket_sign is not None:
        for b2h, _ in hamiltonian_hops(b, N, [1.0] * (N - 1)):
            bra_sign = np.sign(Lc[index[(a, b2h)], col].imag)
            break
    if ket_sign is not None and bra_sign is not None:
        break
check("the commutator sign is written as -i on the ket and +i on the bra (CONSTRUCTION)",
      ket_sign == -1.0 and bra_sign == +1.0,
      f"ket {ket_sign:+.0f}i, bra {bra_sign:+.0f}i; no spectral anchor can see this")


def build_flipped(j_scalar):
    """The same block with the bra-side sign flipped: not a Lindbladian, the blindness probe."""
    pairs_f, index_f = block_basis(N, P, Q)
    Lf = np.zeros((len(pairs_f),) * 2, dtype=complex)
    bonds = [j_scalar] * (N - 1)
    for col, (a, b) in enumerate(pairs_f):
        Lf[col, col] += -2.0 * sum(GAMMA[l] for l in range(N)
                                   if bool(a & bit(l, N)) != bool(b & bit(l, N)))
        for a2h, amp in hamiltonian_hops(a, N, bonds):
            Lf[index_f[(a2h, b)], col] += -1j * amp
        for b2h, amp in hamiltonian_hops(b, N, bonds):
            Lf[index_f[(a, b2h)], col] += -1j * amp          # the flip
    return Lf


# why that check has to be a construction check: measured here, at the anchor coupling
J_ANCHOR = SCALE * QPHYS
Lc_a, _, Lo_a = sectors(J_ANCHOR)
Lf_a = build_flipped(J_ANCHOR)
Lof_a = UO.T @ Lf_a @ UO
block_diff = np.max(np.abs(Lc_a - Lf_a))
spec_c = np.sort_complex(np.linalg.eigvals(Lo_a))
spec_f = np.sort_complex(np.linalg.eigvals(Lof_a))
conj_res = np.max(np.abs(np.sort_complex(np.conj(spec_c)) - spec_f))
s_correct = sigma3(Lo_a, LAM0)
s_flipped = sigma3(Lof_a, LAM0)
read_diff = np.max(np.abs(s_correct - s_flipped))
check("the flip is invisible to the nullity read, which is why the check above is structural",
      block_diff > 1.0 and conj_res < 1e-12 and read_diff == 0.0,
      f"block differs by {block_diff:.4f}; R-odd spectrum is the conjugate to {conj_res:.2e}; "
      f"the sigma read differs by {read_diff!r}")

# ---------------------------------------------------------------- STEP 2
print("\nSTEP 2  SOURCE: nullity 2 at every negative-real w root (real t, Hermitian)")

for r in negs:
    w = float(r["wSeed"]["real"])
    lam = float(r["lambdaSeed"]["real"])
    j_unit = 1j * np.sqrt(-w)
    _, _, Lo = sectors(SCALE * (j_unit / 2.0))
    n, s = nullity(Lo, lam)
    check("  " + r["id"] + ": nullity 2", n == 2 and s[2] > CLEAR_TOL,
          f"nullity {n}, sigma = {s[0]:.2e} {s[1]:.2e} {s[2]:.2e}")

# ---------------------------------------------------------------- STEP 3
print("\nSTEP 3  CONCLUSION at w0, against a defective and a generic control")

J_DEF_SEED = SCALE * DEFECTIVE_Q_PHYS     # routed through the pinned book, like every other J


def pair_gap(j):
    ev = np.linalg.eigvals(sectors(j)[2])
    ev = ev[np.argsort(np.abs(ev - DEFECTIVE_LAM_SEED))][:2]
    return abs(ev[0] - ev[1])


try:
    res = minimize_scalar(pair_gap, bracket=(J_DEF_SEED - 8e-4, J_DEF_SEED, J_DEF_SEED + 8e-4),
                          method="brent", options={"xtol": 1e-14})
except ValueError as exc:
    # A wrong anchor (a mis-set QPHYS, a wrong book) moves the defective seed out of the
    # bracket and the optimiser raises. Report it as a failed check so the run ends on its
    # own verdict line rather than on a traceback from one screen earlier.
    check("the defective control's bracket still contains a minimum", False,
          f"minimize_scalar refused the bracket around J = {J_DEF_SEED:.6f}: {exc}")
    print("")
    print("RED: " + str(FAIL))
    sys.exit(1)
J_def = res.x
ev = np.linalg.eigvals(sectors(J_def)[2])
lam_def = float(np.mean(ev[np.argsort(np.abs(ev - DEFECTIVE_LAM_SEED))][:2]).real)
# The threshold on the residual gap is not a number but a consequence of the local law:
# at a defective coalescence the pair splits as C*sqrt|J - J*|, so a Brent xtol of 1e-14 in
# J floors the attainable gap near 1e-7. Measure the exponent instead of asserting a bound.
hs = [1e-10, 1e-8, 1e-6, 1e-4]
gaps = [pair_gap(J_def + h) for h in hs]
expo = np.polyfit(np.log(hs), np.log(gaps), 1)[0]
check("the defective control is a SQUARE-ROOT coalescence (exponent 1/2, so a Jordan block)",
      abs(expo - 0.5) < 0.02 and res.fun < 100 * np.sqrt(1e-14),
      f"J = {J_def:.12f}, gap exponent = {expo:.4f}, residual gap = {res.fun:.2e}")

# The exponent is the two-sided discriminator, so the SEMISIMPLE side is measured too and
# not merely named: a semisimple pair separates LINEARLY in the coupling, a defective one
# as a square root. A reading between the two would say the classification is wrong.
def w0_pair_gap(j):
    ev = np.linalg.eigvals(sectors(j)[2])
    two = ev[np.argsort(np.abs(ev - LAM0))][:2]
    return abs(two[0] - two[1])


J_W0 = SCALE * QPHYS
w0_hs = [1e-10, 1e-8, 1e-6, 1e-4]
w0_up = np.polyfit(np.log(w0_hs), np.log([w0_pair_gap(J_W0 + h) for h in w0_hs]), 1)[0]
w0_dn = np.polyfit(np.log(w0_hs), np.log([w0_pair_gap(J_W0 - h) for h in w0_hs]), 1)[0]
check("the w0 pair separates LINEARLY (exponent 1, against the defective control's 1/2)",
      abs(w0_up - 1.0) < 0.01 and abs(w0_dn - 1.0) < 0.01,
      f"exponent {w0_up:.7f} from above and {w0_dn:.7f} from below, against {expo:.4f} at "
      f"the defective point; the gap at J* itself is {w0_pair_gap(J_W0):.2e}")

for tag, j, lam, want in (("A2 locus w0        (claim: SEMISIMPLE)", SCALE * QPHYS, LAM0, 2),
                          ("defective EP control (must read 1)", J_def, lam_def, 1),
                          ("generic-lambda control (must read 0)", SCALE * QPHYS, LAM0 + 0.137, 0)):
    _, _, Lo = sectors(j)
    n, s = nullity(Lo, lam)
    check("  " + tag, n == want and s[2] > CLEAR_TOL,
          f"nullity {n}, sigma = {s[0]:.2e} {s[1]:.2e} {s[2]:.2e}")

# ---------------------------------------------------------------- STEP 4
print("\nSTEP 4  ERROR MODEL: what the thresholds rest on, measured over eight decades")

DELTAS = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4]
print(f"        {'dJ/J':>8}  {'sigma2 at w0':>14} {'ratio':>10}   {'sigma2 at the defective pt':>26}")
ratios = []
for d in DELTAS:
    bonds_w0 = [SCALE * QPHYS * (1 + d)] * (N - 1)
    bonds_df = [J_def * (1 + d)] * (N - 1)
    s_w0 = sigma3(sectors(0.0, bonds_w0)[2], LAM0)
    s_df = sigma3(sectors(0.0, bonds_df)[2], lam_def)
    ratios.append(s_w0[1] / d)
    print(f"        {d:8.0e}  {s_w0[1]:14.4e} {s_w0[1]/d:10.4e}   {s_df[1]:26.4e}")
spread = max(ratios) / min(ratios)
check("sigma2 at w0 is LINEAR in the UNIFORM coupling detune (ratio flat over eight decades)",
      spread < 1.02, f"max/min ratio = {spread:.5f}, c = {np.mean(ratios):.4e} (uniform detune)")

# Which directions the law above covers, and it is a theorem rather than a sample.
# L is affine in the bond couplings, so a bond detune gives V = sum_b w_b * J_b * dL/dJ_b.
# The chain reflection R is a real symmetric involution carrying dL/dJ_b to dL/dJ_{N-2-b},
# so an ANTI-palindromic delta pattern (w reversed = -w) obeys R V R = -V. For two R-odd
# vectors u, v one then has v^T V u = (Rv)^T (RVR) (Ru) = -v^T V u, hence U_O^T V U_O
# vanishes IDENTICALLY. That is a property of the DELTA alone: it holds at every N, on any
# base, symmetric or not, and the checks below use the exact integer orbit basis so the
# statement is read as == 0.0 rather than against a threshold.
#
# What the BASE's symmetry buys is different and is easy to mis-state: it makes the R-odd
# subspace L-INVARIANT, so that the restricted sigma2 is a spectral quantity of an actual
# invariant block rather than a compression of one. The control below separates the two.
L_base = sectors(0.0, [SCALE * QPHYS] * (N - 1))[0]
BASE_ODD_INT = UO_INT.T @ sectors(0.0, [SCALE * QPHYS] * (N - 1))[0] @ UO_INT
constants = {}
antis = []
for pattern, palindromic in (([1, 1, 1, 1], True), ([1, 0, 0, 1], True), ([1, -1, -1, 1], True),
                             ([1, -1, 1, -1], False), ([1, 0, 0, -1], False), ([0, 1, -1, 0], False)):
    V = sectors(0.0, [SCALE * QPHYS * (1 + w) for w in pattern])[0] - L_base
    reach_exact = np.max(np.abs(UO_INT.T @ V @ UO_INT))     # integer basis: an EXACT route
    if palindromic:
        dirs = [sigma3(sectors(0.0, [SCALE * QPHYS * (1 + d * w) for w in pattern])[2], LAM0)[1] / d
                for d in DELTAS]
        constants[tuple(pattern)] = float(np.mean(dirs))
        check(f"palindromic detune {pattern} reaches the R-odd sector and responds linearly",
              reach_exact > 0.0 and max(dirs) / min(dirs) < 1.05,
              f"exact-basis reach = {reach_exact:.4e}, c = {np.mean(dirs):.4e}, "
              f"ratio spread = {max(dirs) / min(dirs):.5f}")
    else:
        # The nonzero guard has an exact route too: max|V| is J*max|w| to the last bit, so it
        # is compared to that rather than to a threshold. Without it the pattern [0,0,0,0]
        # would pass on a perturbation that is not there.
        want_v = SCALE * QPHYS * max(abs(w) for w in pattern)
        check(f"anti-palindromic detune {pattern} misses the R-odd sector EXACTLY (theorem, == 0.0)",
              reach_exact == 0.0 and np.max(np.abs(V)) == want_v,
              f"max|U_O^T V U_O| in the integer orbit basis = {reach_exact!r}, on a V with "
              f"max|V| = {np.max(np.abs(V))!r} = J*max|w|")
        # The consequence, and it needs no threshold at all: the restricted block is not
        # merely hard to move, it is the SAME MATRIX, bit for bit, at every detune size.
        # Any function of it, sigma2 included, is therefore unchanged by construction.
        worst = max(np.max(np.abs(UO_INT.T @ sectors(0.0, [SCALE * QPHYS * (1 + d * w)
                                                           for w in pattern])[0] @ UO_INT
                                  - BASE_ODD_INT))
                    for d in DELTAS + [1e-3, 1e-2])
        antis += [sigma3(sectors(0.0, [SCALE * QPHYS * (1 + d * w) for w in pattern])[2], LAM0)[1]
                  for d in DELTAS + [1e-3, 1e-2]]
        check(f"and the R-odd block is then UNCHANGED bit for bit over ten decades",
              worst == 0.0,
              f"max|L_O(detuned) - L_O(base)| over 7 detunes = {worst!r}")

# The blindness is a PROJECTION, not the physics. On the full block the same anti-palindromic
# detune splits the pair quadratically, through the even sector. The law is the FLATNESS of
# split/delta^2 across the decades; its control is the palindromic direction read at the SAME
# power, where the split is linear so the quotient must run away as 1/delta, exactly 1e4 over
# these three decades.
def split_over(power, pattern):
    out = []
    for d in (1e-6, 1e-4, 1e-2):
        lf = np.linalg.eigvals(sectors(0.0, [SCALE * QPHYS * (1 + d * w) for w in pattern])[0])
        lo = np.argsort(np.abs(lf - LAM0))[:2]
        out.append(abs(lf[lo[0]] - lf[lo[1]]) / d ** power)
    return out


blind = split_over(2, [1, -1, 1, -1])
seeing = split_over(2, [1, 1, 1, 1])
check("the sector's blindness is a PROJECTION, not the physics: on the full block the same "
      "anti-palindromic detune splits the pair QUADRATICALLY",
      max(blind) / min(blind) < 1.02 and abs(max(seeing) / min(seeing) / 1e4 - 1) < 0.02,
      f"split/delta^2 = {[round(b, 3) for b in blind]}, flat to {max(blind) / min(blind):.5f}; "
      f"the palindromic control read at the same power runs away by "
      f"{max(seeing) / min(seeing):.4e} against the 1e4 its linear law predicts; at dJ/J = 1e-2 "
      f"the pair is already {np.mean(blind) * 1e-4:.2e} apart while the R-odd block has not "
      f"changed a bit")

# What the ORTHONORMAL basis reads instead, and why it is not a second finding: the same zero
# arrives there through a 1/sqrt(2) scaling, so sigma2 wanders over the SVD's last bits.
check("the orthonormal reading of that same zero is float noise, not a response",
      max(antis) / min(antis) < 3.0 and max(antis) < 1e-15,
      f"sigma2 over {len(antis)} anti-palindromic reads spans [{min(antis):.4e}, "
      f"{max(antis):.4e}], a factor {max(antis) / min(antis):.2f} with no trend in delta")

# The constants are NOT set by the parity: all three palindromic directions share a parity
# and give three different c. Parity decides only whether there is a response to be linear in.
spread_c = max(constants.values()) / min(constants.values())
check("the constant c is direction-dependent WITHIN one parity (so parity fixes existence, not size)",
      spread_c > 1.5, f"c ranges over {sorted(round(v, 6) for v in constants.values())}, "
      f"a factor {spread_c:.3f} across three directions of the SAME parity")

# The hypothesis, measured, and it is about INVARIANCE rather than about the theorem above:
# on a NON-palindromic base the reflection stops commuting with L, so the R-odd subspace is
# not L-invariant and a restricted sigma2 there is a compression rather than a sub-spectrum.
# The theorem itself survives: an anti-palindromic DELTA still projects to exactly zero.
bad_base = [SCALE * QPHYS * x for x in (1.0, 1.3, 0.7, 1.9)]
L_bad = sectors(0.0, bad_base)[0]
comm_bad = np.max(np.abs(Rperm @ L_bad - L_bad @ Rperm))
V_bad = sectors(0.0, [b + SCALE * QPHYS * w for b, w in zip(bad_base, [1, -1, 1, -1])])[0] - L_bad
reach_bad = np.max(np.abs(UO_INT.T @ V_bad @ UO_INT))
check("a NON-palindromic base breaks the sector's INVARIANCE while the theorem still holds",
      comm_bad > 1e-6 and reach_bad == 0.0,
      f"max|RL - LR| = {comm_bad:.4f} against 0.0 on the symmetric base, so the R-odd subspace "
      f"is no longer L-invariant; the anti-palindromic delta nevertheless projects to "
      f"{reach_bad!r}, because R V R = -V is a property of the DELTA and not of the base")

tight_s3 = min(sigma3(sectors(SCALE * (1j * np.sqrt(-float(r["wSeed"]["real"])) / 2.0))[2],
                      float(r["lambdaSeed"]["real"]))[2] for r in negs)
null_w0 = sigma3(sectors(SCALE * QPHYS)[2], LAM0)[1]
C_LAW = float(np.mean(ratios))       # STEP 4 measured this on the UNIFORM detune
# ULP of J is a count of absolute steps, so the relative dJ/J must be multiplied by J first.
ULP_PER_J = (SCALE * QPHYS) / np.spacing(SCALE * QPHYS)
print("")
print(f"        THRESHOLD WINDOW, in the coupling unit STEP 4 measured on the uniform detune")
print(f"        (sigma2 = {C_LAW:.4e}*dJ/J; the other two directions carry their own constant):")
print(f"        the null pair at w0 implies dJ/J = {null_w0 / C_LAW:.2e}, which is "
      f"{null_w0 / C_LAW * ULP_PER_J:.1f} units in the last place of J: no coupling error.")
print(f"        NULL_TOL = {NULL_TOL:.0e} implies {NULL_TOL / C_LAW * ULP_PER_J:.1e} ULP, so the "
      f"null read has {NULL_TOL / null_w0:.1e} of room against it.")
print(f"        CLEAR_TOL = {CLEAR_TOL:.0e} is not in that unit: it separates a nonzero singular "
      f"value from the pair, and the tightest one over the Hermitian anchors is {tight_s3:.2e}, a "
      f"factor {tight_s3 / CLEAR_TOL:.2f} above it. That factor is the narrow side, and it is a "
      f"separation of the operator rather than a measurement error.")

print("\n" + ("ALL GREEN" if not FAIL else "RED: " + str(FAIL)))
sys.exit(1 if FAIL else 0)
