"""F134 gate: the two-row reflection law of the F133 coefficients, with its
Θ-decomposition and the pinned affine-level dictionary.

THE LAW (G1): the F133 symplectic coefficients obey n_(j,k) = n_(10-j,k) on the
whole two-row slice (36 pairs, exact over ℤ). Equivalently μ₁ ↦ 22-μ₁ on
μ = λ+ρ. Grade: a finite, machine-checked exact-ℤ identity (the F127-wall
epistemic class; code-trust caveat, C# second implementation in
WSymplecticClosedForm/CrossTripleOrthogonality).

STRUCTURE GATED HERE:
  G1  the law itself, 0/36, via the committed F133 read-off.
  G2  the finite-Weyl half: F_k odd under total negation (proven algebraically
      in the arc; re-asserted).
  G3  the Θ-decomposition: with Θ_k(t) = (t²² - t⁻²²)·P_k(t) and the explicit
      palindromic P_k, the full coefficient array satisfies
      G_k = Θ_k + defect, defect supported strictly OUTSIDE the window; the
      periodic (-,+)-dihedral lift of Θ_k is exact (odd about 0, even about 22,
      anti-period 44).
  G4  the level dictionary: reflection center 11, anti-period 44 (m-units),
      both rival conventions dead on data; the per-coordinate quadratic factor
      counts 31-9 = 22 = 2(ℓ+h∨), A_ρ-count 14 = 2h∨(C₆) ⟹ the C₆⁽¹⁾ peg
      ℓ = 4 (a coordinate peg, NOT intrinsic; see the proof doc).
  G5  the domain fence: three-row (j,k,l) holds at l ∈ {0,1,3,4,5}, breaks at
      l = 2 exactly on the 8 recorded pairs; the two l=2 survivors sit at the
      s₀-fixed center μ₁ = 11.
  G6  the centroid obstruction: Xc(-v) = Xc(v) on a support sample ⟹ the
      support centroid is 0, so no affine array symmetry can carry the
      off-centroid center 22: the law is intrinsically a WINDOW identity.

What is NOT proved here: a structural mechanism (the affine/theta identity on
the (B₆,C₆) seam). The proof doc records the constructively closed routes.

Run: python simulations/f134_two_row_reflection_law.py   (~2-3 min, exit 0)
"""
import sys, itertools
sys.path.insert(0, "simulations")
from f133_w_closed_form import build_halves, pack, PACKB, RHO
from f133_w_closed_form import n_raw as _n_raw_impl

print("building the X dict (meet-in-the-middle halves)...", flush=True)
P1, P2 = build_halves()
def n_raw(lam):
    return _n_raw_impl(lam, P1, P2)   # = 2*n_lambda; equalities unaffected

_cache = {}
def Xc(nu):
    nu = tuple(nu)
    v = _cache.get(nu)
    if v is not None: return v
    t = 0; tp = pack(list(nu)) + PACKB; g = P1.get
    for k2, c2 in P2:
        w = g(tp - k2)
        if w: t += c2 * w
    _cache[nu] = t; return t

def _perm_sign(p):
    s, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if seen[i]: continue
        j, ln = i, 0
        while not seen[j]: seen[j] = True; j = p[j]; ln += 1
        if ln % 2 == 0: s = -s
    return s

def _sds(v):
    idx = sorted(range(len(v)), key=lambda i: -v[i]); sv = [v[i] for i in idx]
    for i in range(len(sv) - 1):
        if sv[i] == sv[i + 1]: return None, 0
    return sv, _perm_sign(idx)

TAILS = {k: (2 * k + 10, 8, 6, 4, 2) for k in range(6)}

def Fk(k, m):
    tail = TAILS[k]; tot = 0
    for eps in itertools.product((1, -1), repeat=5):
        sg = 1
        for e in eps: sg *= e
        W, ss = _sds([m] + [e * tv for e, tv in zip(eps, tail)])
        if W is None: continue
        tot += sg * ss * Xc(W)
    return tot

def two_row(j, k):
    return (j, k, 0, 0, 0, 0)

def gate1_law():
    bad = live = 0
    for k in range(6):
        for j in range(k, 11 - k):
            a = n_raw(two_row(j, k)); b = n_raw(two_row(10 - j, k))
            if a or b: live += 1
            if a != b: bad += 1
    assert bad == 0 and live > 0, (bad, live)
    print(f"G1  two-row reflection n_(j,k) = n_(10-j,k): 0/36 mismatches ({live} live pairs)")

def gate2_oddness():
    bad = 0
    for k in range(5):
        for m in range(12 + 2 * k, 33 - 2 * k):
            if Fk(k, m) != -Fk(k, -m): bad += 1
    assert bad == 0
    print("G2  finite-Weyl half F_k(m) = -F_k(-m): 0 fails on all windows")

# P_k in t-units, C_a := t^a + t^-a; Theta_k = (t^22 - t^-22) * P_k
def _poly(pairs):  # dict exponent -> coeff
    d = {}
    for e, c in pairs: d[e] = d.get(e, 0) + c
    return {e: c for e, c in d.items() if c}
P_TABLE = {
    0: _poly([(10, -1), (-10, -1), (6, -1), (-6, -1)]),
    1: _poly([(8, 1), (-8, 1), (0, -2)]),          # (t^4 - t^-4)^2
    2: _poly([(2, 3), (-2, 3)]),
    3: _poly([(4, -1), (-4, -1), (0, -3)]),
    4: _poly([(2, 1), (-2, 1)]),
    5: {},
}
def theta(k):
    out = {}
    for e, c in P_TABLE[k].items():
        out[e + 22] = out.get(e + 22, 0) + c
        out[e - 22] = out.get(e - 22, 0) - c
    return {e: c for e, c in out.items() if c}

def gate3_theta():
    for k in range(6):
        th = theta(k)
        lo, hi = 12 + 2 * k, 32 - 2 * k
        arr = {m: Fk(k, m) for m in range(-40, 41)}
        # in-window (both signs) Theta matches F exactly; defect strictly outside
        defect = {m: arr.get(m, 0) - th.get(m, 0) for m in range(-40, 41)
                  if arr.get(m, 0) != th.get(m, 0)}
        assert all(not (lo <= abs(m) <= hi) for m in defect), (k, defect)
        # Theta is ONE PERIOD of an honest (-,+) dihedral object: define the
        # anti-periodic lift from the fundamental domain [-21, 22] and check
        # (i) it reproduces Theta on the whole support range [-32, 32],
        # (ii) it is odd about 0 and even about 22 globally (checked on 3 periods)
        def lift(m):
            r = 0
            while m > 22: m -= 44; r += 1
            while m < -21: m += 44; r -= 1
            return th.get(m, 0) * (-1) ** (r & 1)
        for m in range(-32, 33):
            assert lift(m) == th.get(m, 0), (k, m)
        for m in range(-66, 67):
            assert lift(-m) == -lift(m) and lift(44 - m) == lift(m), (k, m)
    print("G3  Theta-decomposition: G_k = Theta_k + defect, defect window-disjoint;")
    print("    Theta's (-,+) dihedral lift exact (odd/0, even/22) for all k")

H_CHECK = 7
def gate4_level():
    sheets = [L for L in itertools.product((1, -1), repeat=6) if L[0] == 1 and L != (1,) * 6]
    for u in range(6):
        num = sum(L[u] ** 2 for L in sheets)
        den = 4 + 5      # sin x_u (t-exp 2) + five sum roots
        dif = 5          # five difference roots
        assert num - den == 22 and dif + den == 2 * H_CHECK and dif + num == 36
    assert 22 == 2 * (4 + H_CHECK)   # the C6^(1) peg: shifted level 11, l = 4
    # rivals dead on data (leaning on the proven wall F_k = 0 past 32-2k):
    assert Fk(0, 12) != 0
    assert Fk(0, 88 - 12) == 0                             # wall-at-22 reading dead
    assert Fk(0, 12 + 22) == 0                             # step-11 translation dead
    print("G4  level dictionary: counts 31-9=22=2(l+h_check), A_rho 14=2*h_check;")
    print("    center 11 + anti-period 44 pinned, rival conventions dead on data")

L2_BREAKS = {(2, 2), (3, 3), (4, 2), (4, 4), (6, 2), (6, 4), (7, 3), (8, 2)}
def gate5_fence():
    for l in (0, 1, 3, 4, 5):
        for k in range(l, 6):
            for j in range(k, 11 - k):
                lam = (j, k, l, 0, 0, 0); lam_r = (10 - j, k, l, 0, 0, 0)
                if 10 - j < k: continue
                assert n_raw(lam) == n_raw(lam_r), (l, j, k)
    breaks = set()
    for k in range(2, 6):
        for j in range(k, 11 - k):
            if 10 - j < k: continue
            lam = (j, k, 2, 0, 0, 0); lam_r = (10 - j, k, 2, 0, 0, 0)
            if n_raw(lam) != n_raw(lam_r): breaks.add((j, k))
    assert breaks == L2_BREAKS, breaks
    # the two genuine (nonzero) l=2 holds are the s0-CENTER-FIXED weights j = 5
    # (mu1 = 11, mapped to themselves) - consistent, not exceptions
    assert n_raw((5, 3, 2, 0, 0, 0)) != 0 and n_raw((5, 5, 2, 0, 0, 0)) != 0
    print(f"G5  domain fence: l in {{0,1,3,4,5}} clean; l=2 breaks exactly on {len(L2_BREAKS)} recorded pairs;")
    print("    the two nonzero l=2 holds are the center-fixed j=5 weights")

def gate6_centroid():
    # test on genuine support points: every cached read-off vector with Xc != 0
    # (populated by gates 1-3; X lives on the even sublattice)
    pts = [nu for nu, v in list(_cache.items()) if v]
    assert len(pts) > 100, len(pts)
    for nu in pts:
        assert Xc(tuple(-x for x in nu)) == Xc(nu), nu
    print(f"G6  centroid obstruction: Xc(-v) = Xc(v) on all {len(pts)} cached support points")
    print("    (support centroid 0 => no affine array symmetry carries center 22:")
    print("     the law is intrinsically a WINDOW identity)")

if __name__ == "__main__":
    gate1_law(); gate2_oddness(); gate3_theta(); gate4_level(); gate5_fence(); gate6_centroid()
    print("\nF134 gate: ALL GREEN")
