r"""core_grid.py  (grid+CRT engine for grid_proof_close / grid_proof_sweep).

Vectorized, RIGOROUS deterministic grid engine for the Q(i) residue/endpoint certificate.

Two levers over the original scalar `grid_proof_close.py`:

  LEVER 1 (monic quotient reduction).  Where every C-factor's leading w1-coefficient is a
  pure MONOMIAL (z3-free single term) -- the 6 linear singletons and the 3 inc-32 quads --
  the leading coefficient L of Pi_C is a UNIT at every nonzero grid point (its Qz-norm is a
  nonzero monomial), so we reduce P_C in S_p = R_p[w1]/(Pi_C) by inverting L per grid point
  (rp_inv) instead of the division-free pseudo-remainder.  The true remainder is then a genuine
  Laurent polynomial whose z1,z2,w2 span is NOT inflated by the L^k pseudo-remainder factor, so
  the grid shrinks (the pseudo-remainder inflates the certified degree ~2x).  For components with
  a MULTI-TERM leading coefficient (C06-C09 quads, the 12 pairs) L is not a unit, the monic
  remainder is a genuine rational function, and the polynomial grid lemma would not apply to the
  small grid; those keep the division-free pseudo-remainder at their (rigorous) pseudo span.

  SKIP-AND-ENLARGE (monic path).  A grid point where L (or a needed unit) vanishes is UNUSABLE.
  ADJUSTED TENSOR-GRID LEMMA:  a Laurent polynomial of z1,z2,w2-degree <= d_v per axis (after
  shifting min-exp to 0) that vanishes on a full product grid of d_v+1 DISTINCT GOOD residues per
  axis is identically 0.  Bad coordinates are simply not part of the grid, so replacing a bad axis
  value by a fresh good one is legitimate as long as every remaining coordinate is distinct and the
  count stays >= d_v+1.  For the monic components L's norm is a nonzero monomial, so NO grid point
  is ever bad; the enlarge path is a guarded safety net (asserted, essentially never taken).

  LEVER 2 (vectorization).  The whole (z1 x z2 x w2) tensor grid is evaluated at ONCE as int64
  numpy arrays (30-bit primes keep a*b < 2^63; every product is reduced mod p BEFORE the next
  multiply).  The original code vectorized only the z1 axis and looped z2 x w2 in Python; this
  removes that double loop.

RIGOROUS DEGREE BOUND.  The certified object is the reduction remainder as a Laurent polynomial in
z1,z2,w2.  We bound its per-axis exponent span by running the EXACT SAME eval + reduce over a
tropical "span ring" (each field element carries a (min,max) exponent interval per axis; multiply
adds intervals, add/subtract unions them, and the Sz = z1+1/z1+z2+1/z2 norm term contributes its
(-1,1) intervals exactly as in the arithmetic).  This is a guaranteed upper bound on the remainder
span with no analytic guesswork; the grid uses that bound + 1 per axis.

Two fixed bugs preserved: partial products are reduced mod p before combining (int64 overflow), and
Laurent divisors are shifted to min-exp 0 before w1-degree bookkeeping (clear_w1).

Authors: Thomas Wicht and Claude, 2026-07-10.
"""
import sys, time, math
from fractions import Fraction
import numpy as np
sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
from grid_proof_close import (build_structure, component_plan, gen_primes, tonelli_i,
                               factor_key, classify_factor)

W1 = 3  # w1 exponent index in the 6-tuple monomial key

# ===================================================================== coefficient -> GF(p)[i]
def gcoef(c, p, ii):
    cr = (int(c.re.numerator) * pow(int(c.re.denominator), p - 2, p)) % p
    ci = (int(c.im.numerator) * pow(int(c.im.denominator), p - 2, p)) % p
    return (cr + ii * ci) % p

# ===================================================================== ARRAY ring (mod p, 3D)
# A "field element" F is either None (zero) or an int64 numpy array (broadcastable over the grid).
# An R_p element is a pair (a0, a1) of F's:  value = a0 + a1*z3,   z3^2 = -Sz z3 - 1.
class ArrayRing:
    __slots__ = ("p", "ii", "Sz", "z1", "z2", "w2", "z1i", "z2i", "w2i", "w1val")
    def __init__(self, p, z1, z2, w2, Sz, ii=None):
        self.p = p; self.Sz = Sz; self.ii = ii; self.w1val = None
        self.z1 = z1; self.z2 = z2; self.w2 = w2
        self.z1i = pow_arr(z1, -1, p); self.z2i = pow_arr(z2, -1, p); self.w2i = pow_arr(w2, -1, p)
    # field ops
    def fmul(s, a, b):
        if a is None or b is None: return None
        return (a * b) % s.p
    def fadd(s, a, b):
        if a is None: return b
        if b is None: return a
        return (a + b) % s.p
    def fsub(s, a, b):
        if b is None: return a
        if a is None: return (-b) % s.p
        return (a - b) % s.p
    def fneg(s, a):
        return None if a is None else (-a) % s.p
    def from_mono(s, coef_val, e0, e1, e4):
        # coef_val is already a GF(p) scalar (int); build coef * z1^e0 z2^e1 w2^e4
        v = np.array(coef_val, dtype=np.int64)
        if e0: v = (v * axis_pow(s.z1, s.z1i, e0, s.p))
        v = v % s.p
        if e1: v = (v * axis_pow(s.z2, s.z2i, e1, s.p)) % s.p
        if e4: v = (v * axis_pow(s.w2, s.w2i, e4, s.p)) % s.p
        return v
    def fis_zero(s, a):
        return a is None or not np.any(a)

# ===================================================================== SPAN ring (tropical)
# A field element is None (zero) or a dict {0:(lo,hi),1:(lo,hi),4:(lo,hi)} of exponent intervals
# for axes z1(0), z2(1), w2(4).  mul adds intervals; add/sub unions; from_mono sets a point.
_AX = (0, 1, 4)
def _ivadd(x, y): return (x[0] + y[0], x[1] + y[1])
def _ivun(x, y): return (min(x[0], y[0]), max(x[1], y[1]))
class SpanRing:
    __slots__ = ("Sz",)
    def __init__(self):
        self.Sz = {0: (-1, 1), 1: (-1, 1), 4: (0, 0)}  # z1+1/z1+z2+1/z2
    def fmul(s, a, b):
        if a is None or b is None: return None
        return {ax: _ivadd(a[ax], b[ax]) for ax in _AX}
    def fadd(s, a, b):
        if a is None: return b
        if b is None: return a
        return {ax: _ivun(a[ax], b[ax]) for ax in _AX}
    def fsub(s, a, b): return s.fadd(a, b)
    def fneg(s, a): return a
    def from_mono(s, coef_val, e0, e1, e4):
        return {0: (e0, e0), 1: (e1, e1), 4: (e4, e4)}
    def fis_zero(s, a): return a is None

# ===================================================================== HEIGHT ring (coefficient L1)
# Rigorous coefficient-height bound, the third semiring beside ARRAY (values) and SPAN (exponents).
# A field element is None (zero) or a 4-tuple (I, e2, e3, e5):  the polynomial component it stands
# for, once its coefficients are cleared to Gaussian INTEGERS by the common denominator
# 2^e2 * 3^e3 * 5^e5, has integer L1-norm  sum_j (|Re num_j| + |Im num_j|)  <=  I.  (2,3,5 are the
# ONLY denominator primes: 2 from C = -1/8 and sin = -i/2, and 3,5 from the fixed w1 = 2,3,5 slices;
# every factor and every relation SZ,SW,Qz,Qw has Gaussian-integer coefficients.)  L1 is submulti-
# plicative for Gaussian integers ((|a|+|b|)(|c|+|d|) >= |ac-bd|+|ad+bc|) and subadditive, so:
#   mul  ->  (I1*I2, e+e)          (numerators multiply, denominators add exponents);
#   add  ->  common denom = max exponents; each side's numerator scaled up by the deficit, then L1
#            summed.  The monomial exponents (e0,e1,e4) never change a coefficient's magnitude, so
#            from_mono ignores them and reads |Re|+|Im| of the (w1-folded) Gaussian-rational coeff.
# The final bound H for an item = max over the remainder's (a0,a1) components of the numerator L1
# cleared to the item's GLOBAL common denominator (see height_bounds / endpoint_height below).
class HeightRing:
    __slots__ = ("Sz", "w1val")
    def __init__(self):
        self.Sz = (4, 0, 0, 0)   # z1 + 1/z1 + z2 + 1/z2 : four unit monomials, integer, L1 = 4
        self.w1val = None
    def fmul(s, a, b):
        if a is None or b is None: return None
        return (a[0] * b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])
    def fadd(s, a, b):
        if a is None: return b
        if b is None: return a
        e2, e3, e5 = max(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3])
        I = (a[0] * (2 ** (e2 - a[1])) * (3 ** (e3 - a[2])) * (5 ** (e5 - a[3]))
             + b[0] * (2 ** (e2 - b[1])) * (3 ** (e3 - b[2])) * (5 ** (e5 - b[3])))
        return (I, e2, e3, e5)
    def fsub(s, a, b): return s.fadd(a, b)   # |A - B| has the same L1 bound as |A| + |B|
    def fneg(s, a): return a
    def fis_zero(s, a): return a is None
    def from_mono(s, coef, e0, e1, e4):
        """coef is a (re, im) pair of Fractions (already w1-folded).  Returns the height element."""
        re, im = coef
        if re == 0 and im == 0:
            return None
        D = 1
        for fr in (re, im):
            if fr != 0:
                D = D * fr.denominator // math.gcd(D, fr.denominator)
        nr = abs(int(re * D)); ni = abs(int(im * D))
        e2 = e3 = e5 = 0; d = D
        while d % 2 == 0: d //= 2; e2 += 1
        while d % 3 == 0: d //= 3; e3 += 1
        while d % 5 == 0: d //= 5; e5 += 1
        assert d == 1, f"height denominator {D} has a prime factor outside {{2,3,5}}"
        return (nr + ni, e2, e3, e5)

def _height_clear(parts):
    """Max cleared-integer coefficient bound over a list of height elements (a common denominator)."""
    parts = [pt for pt in parts if pt is not None]
    if not parts:
        return 0
    E2 = max(pt[1] for pt in parts); E3 = max(pt[2] for pt in parts); E5 = max(pt[3] for pt in parts)
    H = 0
    for (I, e2, e3, e5) in parts:
        H = max(H, I * (2 ** (E2 - e2)) * (3 ** (E3 - e3)) * (5 ** (E5 - e5)))
    return H

# ===================================================================== axis power helpers (array)
def pow_arr(base, e, p):
    if e < 0:
        base = _invarr(base, p); e = -e
    r = np.ones_like(base); b = base % p
    while e:
        if e & 1: r = (r * b) % p
        b = (b * b) % p; e >>= 1
    return r
def _invarr(arr, p):
    r = np.ones_like(arr); b = arr % p; e = p - 2
    while e:
        if e & 1: r = (r * b) % p
        b = (b * b) % p; e >>= 1
    return r
def axis_pow(base, basei, e, p):
    if e == 0: return np.ones_like(base)
    if e > 0:  return pow_arr(base, e, p)
    return pow_arr(basei, -e, p)

# ===================================================================== R_p element ops (ring-generic)
def rp_mul(R, a, b):
    a0, a1 = a; b0, b1 = b
    a0b0 = R.fmul(a0, b0); a1b1 = R.fmul(a1, b1)
    a0b1 = R.fmul(a0, b1); a1b0 = R.fmul(a1, b0)
    r0 = R.fsub(a0b0, a1b1)
    r1 = R.fsub(R.fadd(a0b1, a1b0), R.fmul(a1b1, R.Sz))
    return (r0, r1)

def rp_inv_array(R, a):
    """Invert a=(a0,a1) in R_p (array ring).  Returns (inv, ok) where ok is a bool mask (all-True
    means every grid point is a unit).  norm = a0^2 - a0 a1 Sz + a1^2;  inv = ((a0 - a1 Sz) - a1 z3)/norm."""
    p = R.p; a0, a1 = a
    a0 = np.zeros_like(R.Sz) if a0 is None else a0
    a1 = np.zeros_like(R.Sz) if a1 is None else a1
    a0a0 = (a0 * a0) % p
    a1a1 = (a1 * a1) % p
    a0a1 = (a0 * a1) % p
    a0a1Sz = (a0a1 * R.Sz) % p
    norm = (a0a0 - a0a1Sz + a1a1) % p
    ok = norm != 0
    inv_norm = _invarr(np.where(ok, norm, 1), p)
    b0 = ((a0 - (a1 * R.Sz) % p) % p * inv_norm) % p
    b1 = ((-a1) % p * inv_norm) % p
    return (b0, b1), ok

# ===================================================================== w1-polynomial ops (dict)
def pmul_w1(R, A, B):
    out = {}
    for wa, ab in A.items():
        for wb, cd in B.items():
            pr = rp_mul(R, ab, cd); w = wa + wb
            cur = out.get(w)
            out[w] = pr if cur is None else (R.fadd(cur[0], pr[0]), R.fadd(cur[1], pr[1]))
    return out
def padd_w1(R, A, B):
    out = {w: (ab[0], ab[1]) for w, ab in A.items()}
    for w, cd in B.items():
        cur = out.get(w)
        out[w] = cd if cur is None else (R.fadd(cur[0], cd[0]), R.fadd(cur[1], cd[1]))
    return out
def clear_w1(R, poly):
    poly = {w: ab for w, ab in poly.items() if not (R.fis_zero(ab[0]) and R.fis_zero(ab[1]))}
    if not poly: return {}
    s = min(poly); return {w - s: ab for w, ab in poly.items()}

# ---- pseudo-remainder (division-free, rigorous for ANY leading coeff)
def pseudo_reduce(R, P, PHI):
    P = clear_w1(R, dict(P)); PHI = clear_w1(R, PHI)
    if not PHI: raise RuntimeError("empty divisor")
    m = max(PHI); L = PHI[m]; steps = 0
    while P:
        dP = max(P)
        if dP < m: break
        lead = P[dP]; newP = {}
        for w, ab in P.items(): newP[w] = rp_mul(R, ab, L)
        newP.pop(dP, None)  # leading term L*lead cancels PHI[m]*lead exactly (removed explicitly)
        for w, ab in PHI.items():
            if w == m: continue
            tgt = w + (dP - m); sub = rp_mul(R, ab, lead); cur = newP.get(tgt)
            newP[tgt] = (R.fneg(sub[0]), R.fneg(sub[1])) if cur is None else (R.fsub(cur[0], sub[0]), R.fsub(cur[1], sub[1]))
        P = {w: ab for w, ab in newP.items() if not (R.fis_zero(ab[0]) and R.fis_zero(ab[1]))}
        steps += 1
        if steps > 400: raise RuntimeError("pseudo-rem runaway")
    return P

# ---- monic remainder (inverts L; ARRAY ring only; requires L a unit -> ok mask all True)
def monic_reduce(R, P, PHI):
    P = clear_w1(R, dict(P)); PHI = clear_w1(R, PHI)
    if not PHI: raise RuntimeError("empty divisor")
    m = max(PHI); L = PHI[m]
    Linv, ok = rp_inv_array(R, L)
    if not np.all(ok):
        return None, ok  # bad grid point(s) present
    steps = 0
    while P:
        dP = max(P)
        if dP < m: break
        lead = P.pop(dP); f = rp_mul(R, lead, Linv)  # leading term cancels (f*L == lead)
        for w, ab in PHI.items():
            if w == m: continue
            tgt = w + (dP - m); sub = rp_mul(R, ab, f); cur = P.get(tgt)
            P[tgt] = (R.fneg(sub[0]), R.fneg(sub[1])) if cur is None else (R.fsub(cur[0], sub[0]), R.fsub(cur[1], sub[1]))
        P = {w: ab for w, ab in P.items() if not (R.fis_zero(ab[0]) and R.fis_zero(ab[1]))}
        steps += 1
        if steps > 400: raise RuntimeError("monic runaway")
    return P, ok

# ---- span-ring monic reduce (needs monomial L so Linv is a monomial)
def _mono_inv_span(L):
    """Inverse span of a MONOMIAL R_p element L=(a0,a1) with a1 zero and a0 a single monomial."""
    a0, a1 = L
    assert a1 is None, "monic span reduce requires z3-free (monomial) leading coeff"
    return ({ax: (-a0[ax][0], -a0[ax][1]) for ax in _AX}, None)
def monic_reduce_span(R, P, PHI):
    P = clear_w1(R, dict(P)); PHI = clear_w1(R, PHI)
    m = max(PHI); L = PHI[m]; Linv = _mono_inv_span(L)
    steps = 0
    while P:
        dP = max(P)
        if dP < m: break
        lead = P.pop(dP); f = rp_mul(R, lead, Linv)  # leading cancels by construction
        for w, ab in PHI.items():
            if w == m: continue
            tgt = w + (dP - m); sub = rp_mul(R, ab, f); cur = P.get(tgt)
            P[tgt] = (R.fneg(sub[0]), R.fneg(sub[1])) if cur is None else (R.fsub(cur[0], sub[0]), R.fsub(cur[1], sub[1]))
        P = {w: ab for w, ab in P.items() if not (R.fis_zero(ab[0]) and R.fis_zero(ab[1]))}
        steps += 1
        if steps > 400: raise RuntimeError("monic span runaway")
    return P

# ===================================================================== generic eval of P_C, Pi_C
def _eval_poly(R, poly, p=None, ii=None):
    """6-var Laurent dict -> {w1exp:(a0,a1)} over the ring R (z3 exp in {0,1})."""
    is_span = isinstance(R, SpanRing)
    is_height = isinstance(R, HeightRing)
    if not is_span and not is_height:
        p = R.p; ii = R.ii
    w1v = getattr(R, "w1val", None)  # if set (ArrayRing/HeightRing), fold w1 -> scalar (endpoints)
    out = {}
    for k, c in poly.items():
        if is_span:
            val = R.from_mono(None, k[0], k[1], k[4])  # w1 folds to a constant: no z-span effect
        elif is_height:
            re, im = c.re, c.im
            if w1v is not None and k[3]:
                f = Fraction(w1v) ** k[3]   # exact rational fold (negative w1-exp -> 1/w1v^|e|)
                re, im = re * f, im * f
            val = R.from_mono((re, im), k[0], k[1], k[4])
        else:
            cv = gcoef(c, p, ii)
            if w1v is not None and k[3]:
                e = k[3]
                cv = (cv * (pow(w1v, e, p) if e > 0 else pow(pow(w1v, p - 2, p), -e, p))) % p
            val = R.from_mono(cv, k[0], k[1], k[4])
        w = k[3]; z3e = k[2]
        cur = out.get(w)
        if cur is None:
            cur = (None, None); out[w] = cur
        if z3e == 0: out[w] = (R.fadd(cur[0], val), cur[1])
        else:        out[w] = (cur[0], R.fadd(cur[1], val))
    return out

def _ring_one(R):
    """The multiplicative identity as a field element for the given ring."""
    if isinstance(R, SpanRing): return {0: (0, 0), 1: (0, 0), 4: (0, 0)}
    if isinstance(R, HeightRing): return (1, 0, 0, 0)
    return np.ones_like(R.Sz)

def eval_PC_PHI(R, F, factors, C, which, inc2, distinct, p=None, ii=None, corrupt=False):
    Cset = set(C)
    Cfac = {key: _eval_poly(R, factors[key], p, ii) for key in C}
    dk = [factor_key(d) for d in distinct]
    drp = [_eval_poly(R, d, p, ii) for d in distinct]
    PC = {}; first = True
    for ti in inc2:
        A, B, dens = F[ti]; Num = A if which == 0 else B
        if not Num: continue
        k0, k1 = factor_key(dens[0]), factor_key(dens[1])
        inC0, inC1 = k0 in Cset, k1 in Cset
        used_C = set(x for x, inc in ((k0, inC0), (k1, inC1)) if inc)
        own_co = None if (inC0 and inC1) else (k1 if inC0 else k0)
        term = _eval_poly(R, Num, p, ii)
        if corrupt and first:
            one = _ring_one(R)
            term = padd_w1(R, term, {0: (one, None)}); first = False
        for key in C:
            if key in used_C: continue
            term = pmul_w1(R, term, Cfac[key])
        for d_k, d_rp in zip(dk, drp):
            if d_k == own_co: continue
            term = pmul_w1(R, term, d_rp)
        PC = padd_w1(R, PC, term)
    PHI = {0: (_ring_one(R), None)}
    for key in C: PHI = pmul_w1(R, PHI, Cfac[key])
    return PC, PHI

# ===================================================================== component analysis
def lead_is_monomial(factors, C):
    """True iff every C-factor's top-w1 coefficient is a single z3-free monomial."""
    for key in C:
        d = factors[key]
        wmax = max(k[3] for k in d)
        lead = [k for k in d if k[3] == wmax]
        if len(lead) != 1 or lead[0][2] != 0:   # >1 term, or carries z3
            return False
    return True

def span_bounds(F, factors, C, which, inc2, distinct, monic):
    """Rigorous per-axis (z1,z2,w2) span of the reduction remainder, via the tropical span ring."""
    R = SpanRing()
    PC, PHI = eval_PC_PHI(R, F, factors, C, which, inc2, distinct)
    rem = monic_reduce_span(R, PC, PHI) if monic else pseudo_reduce(R, PC, PHI)
    b = {0: (0, 0), 1: (0, 0), 4: (0, 0)}
    got = False
    for w, ab in rem.items():
        for part in ab:
            if part is None: continue
            got = True
            for ax in _AX: b[ax] = _ivun(b[ax], part[ax])
    return b, got

def component_widths(F, factors, incidence, C):
    """Return (monic:bool, widths dict axis->int) using the max span over which in {0,1}."""
    inc2, distinct = component_plan(F, factors, incidence, C)
    monic = lead_is_monomial(factors, C)
    W = {0: 0, 1: 0, 4: 0}
    for which in (0, 1):
        b, got = span_bounds(F, factors, C, which, inc2, distinct, monic)
        for ax in _AX: W[ax] = max(W[ax], b[ax][1] - b[ax][0])
    return monic, W, inc2, distinct

# ===================================================================== coefficient-height bounds
# The owed rigor item.  The same traversal as span_bounds/component_widths, but over the HEIGHT
# semiring: instead of exponent intervals we carry the cleared-integer L1 numerator bound.  The
# certified object is the DIVISION-FREE pseudo-remainder (never the monic remainder): for a
# monomial-leading (monic) component the monic remainder is the pseudo-remainder divided by the
# single monomial L^steps whose coefficient is a UNIT, so their coefficient magnitudes are equal
# and the pseudo bound is exact for it too; for a multi-leading component the pseudo-remainder IS
# the tested object, its leading-coefficient inflation counted step by step by pseudo_reduce.  So
# pseudo_reduce over the height ring is a rigorous upper bound on the tested object's integer
# coefficients for EVERY component class, with no reliance on the runtime monic/pseudo choice.
def height_bounds(F, factors, C, which, inc2, distinct):
    """Integer coefficient-height (cleared L1) bound of the pseudo-remainder for one which in {0,1}."""
    R = HeightRing()
    PC, PHI = eval_PC_PHI(R, F, factors, C, which, inc2, distinct)
    rem = pseudo_reduce(R, PC, PHI)
    parts = []
    for w, ab in rem.items():
        parts.extend(ab)
    return _height_clear(parts)

def component_height(F, factors, incidence, C, inc2=None, distinct=None):
    """Max cleared-integer coefficient-height bound H of the component's tested pseudo-remainder."""
    if inc2 is None or distinct is None:
        inc2, distinct = component_plan(F, factors, incidence, C)
    return max(height_bounds(F, factors, C, w, inc2, distinct) for w in (0, 1))

def _collapse_w1_height(R, poly):
    """Evaluate a w1-folded 6-var poly to a single height R_p element (a0,a1) (sums all w1 powers)."""
    d = _eval_poly(R, poly)
    a0 = a1 = None
    for w, ab in d.items():
        a0 = R.fadd(a0, ab[0]); a1 = R.fadd(a1, ab[1])
    return (a0, a1)

def endpoint_height(terms, distinct, w1val):
    """Integer coefficient-height bound H of an endpoint numerator N = sum_t Num_t * prod(distinct
    factors except the term's own two), all at the fixed integer w1 = w1val.  Same structured,
    division-free product the endpoint tester evaluates; distinct is {key: poly}, terms a list of
    (Num, key0, key1) as built by the runner's build_endpoint_terms."""
    R = HeightRing(); R.w1val = w1val
    dval = {k: _collapse_w1_height(R, d) for k, d in distinct.items()}
    N = (None, None)
    for (Num, k0, k1) in terms:
        term = _collapse_w1_height(R, Num)
        for k, dv in dval.items():
            if k == k0 or k == k1: continue
            term = rp_mul(R, term, dv)
        N = (R.fadd(N[0], term[0]), R.fadd(N[1], term[1]))
    return _height_clear([N[0], N[1]])

# ===================================================================== grid + run
def _distinct_good(n, p, seed):
    import random
    rng = random.Random(seed); s = []
    seen = set()
    while len(s) < n:
        v = rng.randrange(1, p)
        if v not in seen: seen.add(v); s.append(v)
    return s

def build_ring(p, gz1, gz2, gw2, seedbase, ii=None):
    z1 = np.array(_distinct_good(gz1, p, seedbase), dtype=np.int64).reshape(gz1, 1, 1)
    z2 = np.array(_distinct_good(gz2, p, seedbase + 1), dtype=np.int64).reshape(1, gz2, 1)
    w2 = np.array(_distinct_good(gw2, p, seedbase + 2), dtype=np.int64).reshape(1, 1, gw2)
    z1i = pow_arr(z1, -1, p); z2i = pow_arr(z2, -1, p)
    Sz = (z1 + z1i + z2 + z2i) % p            # shape (gz1,gz2,1) broadcast
    Sz = np.broadcast_to(Sz, (gz1, gz2, gw2)).copy()
    return ArrayRing(p, np.broadcast_to(z1, (gz1, gz2, gw2)).copy(),
                     np.broadcast_to(z2, (gz1, gz2, gw2)).copy(),
                     np.broadcast_to(w2, (gz1, gz2, gw2)).copy(), Sz, ii=ii)

def reduce_and_test(R, F, factors, C, which, inc2, distinct, monic, p, ii, corrupt=False):
    # ALWAYS reduce via monic (invert Pi_C's leading coeff per grid point): it touches only Pi_C's
    # few w1-terms per step (the pseudo-remainder multiplies the WHOLE P by L each step -> ~10x more
    # array multiplies) and never inflates coefficients.  At a GOOD grid point (leading-coeff norm
    # != 0) the monic remainder is 0  <=>  the pseudo-remainder is 0, so testing the monic remainder
    # on the pseudo-span grid (the rigorous certified degree) is an equivalent, cheaper certificate.
    # Bad points (norm == 0) are excluded and the axis enlarged by the caller (skip-and-enlarge).
    PC, PHI = eval_PC_PHI(R, F, factors, C, which, inc2, distinct, p, ii, corrupt=corrupt)
    rem, ok = monic_reduce(R, PC, PHI)
    if rem is None: return None  # bad point(s) present -> caller enlarges the grid
    return any((not R.fis_zero(ab[0])) or (not R.fis_zero(ab[1])) for ab in rem.values())

def run_component_prime(F, factors, incidence, C, p, monic, W, inc2, distinct, seedbase=100, do_control=False):
    """Full-tensor-grid reduction test at a single prime.  Returns (proved:bool, control_ok:bool|None, secs)."""
    t0 = time.time()
    ii = tonelli_i(p)
    gz1, gz2, gw2 = W[0] + 1, W[1] + 1, W[4] + 1
    # skip-and-enlarge: resample whole ring if any monic leading-coeff norm vanishes on the grid
    for attempt in range(6):
        R = build_ring(p, gz1, gz2, gw2, seedbase + 37 * attempt, ii)
        nz = None
        for which in (0, 1):
            nz = reduce_and_test(R, F, factors, C, which, inc2, distinct, monic, p, ii)
            if nz is None: break  # bad point -> enlarge
            if nz: break  # nonzero remainder -> FAIL (report immediately)
        if nz is None:
            gz1 += 1; gz2 += 1; gw2 += 1  # enlarge and resample
            continue
        proved = not nz
        control_ok = None
        if do_control:
            c = reduce_and_test(R, F, factors, C, 0, inc2, distinct, monic, p, ii, corrupt=True)
            control_ok = bool(c)  # True means corruption produced a nonzero remainder (test discriminates)
        return proved, control_ok, time.time() - t0
    raise RuntimeError("skip-and-enlarge exceeded attempts (unexpected for monomial-lead)")
