"""The two involution sources of F153, as one entry-wise identity, with controls.

WHAT THIS PINS. A joint-popcount block (p,q) of the Liouvillian carries an
antilinear involution about a REAL constant c when a permutation m of its cells
satisfies, entry for entry and with no gauge dressing,

    L[m i, m j] = -conj(L[i, j]) - c * delta_ij.                          (*)

(*) says P^T L P = -conj(L) - c*I with P a permutation matrix, a similarity, so
spec(L) = -conj(spec(L)) - c with algebraic multiplicities; no normality is
needed, and these blocks are decidedly not normal. Diagonalisability is not
needed either, and unlike normality it is not something these blocks lack: at
N = 4, Delta = 1, locus gamma, J = 1.3 the eigenvector-MATRIX condition numbers
kappa(V) of the blocks run 6 to 25 (per-eigenvalue condition numbers 1.0 to
10.4), and defective couplings are isolated points, catalogued for the corner
block by PROOF_R90_FROZEN_DIVISOR section 9, other blocks carrying their own
(e.g. the (1,2) real EP of PROOF_CODIM1_BY_ADDITIVITY section 7). With c real
the pairing is Re -> -Re - c and Im -> Im UNCHANGED. Source 1 below is F89d's
cross-fold (F89CrossFoldSimilarityClaim, live at `inspect --root crossfold`) in
the case where the fold's image block IS the block, i.e. the self-folded block of
PROOF_CODIM1_BY_ADDITIVITY section 7 consequence (b); source 2 is not F89d's map.

TWO MAPS SUPPLY (*), and they are F153's two sources.

  SOURCE 1. If p = N/2 or q = N/2, the bitwise complement of THAT leg maps the
  block to itself and (*) holds with c = 2*sigma, sigma = sum_l gamma_l, at any
  gamma profile and with no pinning. Neither half is new: F153 states the map and
  its profile generality ("sends rate -> Sum gamma - rate = 2*mean(rate), for ANY
  profile") and gates the profile blindness in PinnedBlockFloorClaimTests, while
  PROOF_CODIM1 section 7 (b) states the self-fold with no pinning hypothesis. What
  is run here is the entry-wise matrix form on the NON-pinned index-N/2 blocks.
  Source 1 is the global complement X^N on one leg, so it needs [H, X^N] = 0 and a
  LONGITUDINAL FIELD breaks it; that is its negative control.

  SOURCE 2, AND IT CARRIES A HYPOTHESIS ON H THAT SOURCE 1 DOES NOT. If gamma is
  anti-palindromic (the R90 locus) AND H is itself reflection-symmetric, the site
  reflection on both legs sends a cell's disagreement sum to 2*gbar*|S| - itself,
  S the disagreement set, so the CENTRE is fixed by the size class. The bond
  hypothesis is not decoration: with bonds [1, 1.25, 1.5, ...] the residual is
  1.0 / 2.0 / 3.0 at N = 4 / 5 / 6 on the pinned blocks and 2.0 / 4.0 on the
  gbar = 0 stratum, and it is proportional to the asymmetry with no threshold.
  That is this source's negative control and it is run below. Two regimes:

    gbar != 0 : one centre serves a block exactly when its disagreement window
                {|p-q|, ..., min(p+q, 2N-p-q)} holds ONE class, which is exactly
                min(p,q) = 0 or max(p,q) = N. Then c = 4*gbar*|p-q|.
    gbar == 0 : every class wants centre 0 at once, so EVERY block is served with
                c = 0. PROOF_R90_FROZEN_DIVISOR section 5 already states this.

The gbar != 0 regime is why F140's (1,1)-type corner block is not served: its
window is {0, 2}, the |S| = 2 cells force the constant 8*gbar (F140's frozen root
-4*gbar is that constant's axis) and the diagonal cells force 0. The mismatch is
not merely the same NUMBER as F140's defect, it is the same defect MATRIX reached
by a second and antilinear identity, and the gate below measures that: the whole
residual of the plain-c version sits on the D cells at exactly 8*gbar, is exactly
0.0 away from the block's own diagonal, and on the O cells is zero to within the
same 4*eps*max|L| the ZZ accumulation order buys everywhere else in this file
(0.0 at N = 4, 4.4e-16 to 8.9e-16 at N = 5 and 6). Writing that column as
"exactly 0.0" was itself a repair's overstatement and the run refutes it. And

    L[m i, m j] = -conj(L[i,j]) - 8*gbar*delta_ij + 8*gbar*(P_D)_ij

holds entry for entry, with P_D the projector onto the diagonal cells that
PROOF_R90_FROZEN_DIVISOR section 2 localises the LINEAR tauQ defect on. That is
an identification and not a rhyme. It reaches no COUNT: the tax of section 3.1 is
dim D_minus = floor(N/2), and the size classes say only that there IS a defect
on D, never how large the subtraction is.

WHAT (*) DOES NOT GIVE. A palindrome says nothing about flatness, and (*) is not
what makes F153's criterion true: at UNIFORM gamma a pinned block has one size
class, so its dissipator is the real scalar c = -2*gamma*|p-q| and L = c*Id + i*S
with S exactly Hermitian, which puts every Re lambda at c with no involution and
no locus in the argument. That is a statement about blocks whose dissipator is a
SCALAR and it does not generalise: on a block with two size classes a
block-preserving longitudinal field does move Re lambda (the N = 4 (2,2) block
goes from [-4.0000, 0.0000] to [-3.8680, 0.0000] at uniform gamma = 0.5,
Delta = 0, J = 1.3, h = (0.5, 0.25, 0.75, 0.125); the value is
parameter-dependent, J = 1.0 giving -3.805), so "the
Hamiltonian part is anti-Hermitian, therefore H cannot reach Re lambda" is false
as a general claim and true only with the scalar premise.

WHAT THE NON-PINNED INDEX-N/2 BLOCKS DO IS DELTA-DEPENDENT, and this file states
no mechanism for it. Measured at N = 6 on the locus, converged (identical at
J = 1e4, 1e5 and 1e6; that convergence check is a separate measurement, the
ladder this script runs stopping at 1e4), THREE behaviours and not two: at
Delta = 0 and 1 the
blocks (1,3) and (2,3) saturate on the interval spanned by their size-class
centres -2*gbar*|S|, endpoints included; at Delta = 2 they FLATTEN onto the trace
constant; at Delta = 0.5 they do NEITHER, settling at spreads 0.449 and 0.566
about that constant. The doubly-index-N/2 block (3,3) saturates at every Delta
measured, and every pinned block flattens at every Delta measured.

Three earlier readings died here and the table gates against all of them.
"Twice the window half-width" was a RATE-variable width, agreeing with the lambda
one only at gbar = 1/2, which is where the table used to be run. "A multi-class
block never flattens" was a Delta = 0 and 1 reading. And a J ladder stopping at
100 read Delta = 2 as broken when it had simply not converged: the spread there
is 4.29 at J = 100, 3.94 at 1e3 and flat at the eigensolver's backward error
(~1e-10) from 1e4 on. gbar, Delta and a
ladder that reaches the limit are all in the table now, and each row's reading is
COMPUTED and then compared with the committed one, so a drift either way fails.

THREE TRAPS THIS SCRIPT WALKED INTO AND NOW GATES AGAINST, recorded because each
made an earlier version pass for the wrong reason.

  (T1) An arithmetic ramp is automatically anti-palindromic, so linspace and every
       dyadic ramp lie ON the locus and cannot serve as an off-locus profile.
       Every profile's class is asserted, never commented.
  (T2) On the blocks (0, N/2) and (N/2, N) the two constants COINCIDE identically,
       4*gbar*|p-q| = 4*gbar*N/2 = 2*sigma, at any profile. Source 1 then supplies
       the spectral consequence off the locus as well, so the SPECTRAL column has
       no power on those rows however the reflection behaves. They are marked and
       excluded from the spectral assertion; the entry-wise verdict still stands.
  (T3) A "read, not gated" bucket with no bound cannot fail. An earlier version
       gated source 2 only at dyadic J and printed everything else; an injected
       O(1) break at the non-dyadic coupling passed. The bucket is now bounded by
       a stated error model and the sentence about it is COMPUTED, not printed.

EXACTNESS. The gammas are dyadic and each profile fixes sigma and gbar exactly, so
both constants are exactly representable and (*) has an exact route: compared with
== 0.0. Where a residual survives it is bounded, not tolerated: at a coupling whose
BOND PARTIAL SUMS are not representable, source 2 at Delta != 0 carries one
rounding of the ZZ diagonal. The mechanism is accumulation order, not the coupling
itself, which is why J = pi and J = 1.1 give exactly 0.0 while J = 1.3 does not,
and why N = 4 is exactly 0.0 at every J tried. That residual is a deterministic
function of an input the identity's VALUE does not contain and vanishes for
settings of it, which is this repo's test for reading a rounding rather than
gating it; it is nonetheless bounded here by 4*eps*max|L| so the bucket can fail.

The spectral half has no exact route (an eigensolver on a non-normal block), so it
is a ratio to the block's spectral radius, and the gate reads the gap between the
carried and the broken populations. Sorting for that comparison keys on rounded
(Re, Im) jointly: a naive complex sort scrambles near-tied real parts.

Pauli book: H = J * sum_bonds (XX + YY + Delta*ZZ), plus an optional longitudinal
field. |a><b| at a*d + b, site 0 the most significant bit.

Docs: experiments/WHAT_THE_R90_LOCUS_BUYS.md
Output: simulations/results/two_sources_gate.txt
"""

import os
import sys
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "results", "two_sources_gate.txt")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
_f = open(OUT, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(m=""):
    print(m, flush=True)
    _f.write(m + "\n")


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
EPS = np.finfo(float).eps
popcount = lambda x: bin(x).count("1")


def op_at(P, site, n):
    M = np.array([[1]], dtype=complex)
    for i in range(n):
        M = np.kron(M, P if i == site else I2)
    return M


def liouvillian(n, J, g, delta, h=None):
    """`J` is a scalar or a per-bond sequence of length n-1. The per-bond form
    exists for ONE purpose: source 2's reflection needs H itself reflection
    -symmetric, and a bond list that is not a palindrome is the only negative
    control for that hypothesis. Source 1 is blind to it."""
    d = 2 ** n
    Js = [J] * (n - 1) if np.isscalar(J) else list(J)
    assert len(Js) == n - 1, "a bond list must have n-1 entries"
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        H += Js[b] * (op_at(X, b, n) @ op_at(X, b + 1, n)
                      + op_at(Y, b, n) @ op_at(Y, b + 1, n))
        if delta:
            H += Js[b] * delta * (op_at(Z, b, n) @ op_at(Z, b + 1, n))
    if h is not None:
        for s in range(n):
            H += h[s] * op_at(Z, s, n)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for s in range(n):
        Zs = op_at(Z, s, n)
        L += g[s] * (np.kron(Zs, Zs.T) - np.kron(Id, Id))
    return L


def cells(n, p, q):
    d = 2 ** n
    return [(r, c) for r in range(d) for c in range(d)
            if popcount(r) == p and popcount(c) == q]


def block_of(L, n, p, q):
    d = 2 ** n
    cl = cells(n, p, q)
    ix = [r * d + c for (r, c) in cl]
    return L[np.ix_(ix, ix)], cl


def permutation_of(cl, image):
    """The cell permutation induced by `image`, as an index array."""
    pos = {x: i for i, x in enumerate(cl)}
    return np.array([pos[image(a, b)] for (a, b) in cl])


def antilinear_residual(sub, cl, image, c, m=None):
    """max |L[mi,mj] + conj(L[i,j]) + c*delta| over the block; 0 means (*) holds.

    Exact by construction: the permuted read and the conjugate are the SAME floats
    the builder produced, so a surviving nonzero is the identity failing and not a
    reassociation. Vectorised because the c scan below evaluates this hundreds of
    times per block."""
    if m is None:
        m = permutation_of(cl, image)
    d = sub[np.ix_(m, m)] + np.conj(sub)
    if c:
        d = d + c * np.eye(len(cl))
    return np.abs(d).max()


def best_constant(sub, cl, image, span):
    """Smallest residual over a scan of c wide enough to CONTAIN every candidate.

    An earlier version scanned a fixed [-8, 8] while the candidates 2*sigma and
    4*gbar*|p-q| sat outside it, so every printed minimum was the boundary value.
    The window is derived from the candidates now, and a boundary hit is flagged.
    """
    step = span / 256.0
    m = permutation_of(cl, image)
    base = sub[np.ix_(m, m)] + np.conj(sub)
    diag = np.diag(base)
    offdiag = np.abs(base - np.diag(diag)).max()      # constant in c, hoisted
    ks = np.arange(-256, 257)
    cs = ks * step
    diag_part = np.abs(diag[None, :] + cs[:, None]).max(axis=1)
    vals = np.maximum(offdiag, diag_part)
    i = int(np.argmin(vals))
    best, at = float(vals[i]), float(cs[i])
    return best, at, abs(abs(at) - span) < step


def linear_residual_with_gauge(sub, cl, image, c, n):
    """MirrorWorld Mirror's LINEAR fold leg with the bipartite gauge."""
    pos = {x: i for i, x in enumerate(cl)}
    m = [pos[image(a, b)] for (a, b) in cl]
    bit = lambda x, s: (x >> (n - 1 - s)) & 1
    e = [(-1) ** (sum(bit(a, s) for s in range(n) if s % 2)
                  + sum(bit(b, s) for s in range(n) if s % 2)) for (a, b) in cl]
    worst = 0.0
    for i in range(len(cl)):
        for j in range(len(cl)):
            v = sub[m[i], m[j]] + e[i] * e[j] * sub[i, j] + (c if i == j else 0.0)
            worst = max(worst, abs(v))
    return worst


def spectral_pairing(sub, c):
    """Ratio to the spectral radius of: Re -> -Re - c AND Im preserved."""
    ev = np.linalg.eigvals(sub)
    radius = max(np.abs(ev).max(), 1.0)
    partner = -np.conj(ev) - c
    key = lambda z: np.lexsort((np.round(z.imag, 9), np.round(z.real, 9)))
    return np.abs(ev[key(ev)] - partner[key(partner)]).max() / radius


def re_spread(sub):
    re = np.linalg.eigvals(sub).real
    return re.max() - re.min()


def window(n, p, q):
    return sorted({popcount(r ^ c) for (r, c) in cells(n, p, q)})


def pinned(n, p, q):
    return min(p, q) == 0 or max(p, q) == n


def constants_collide(n, p, q, g):
    """T2: on (0,N/2) and (N/2,N) the two constants are identically equal, so the
    spectral column is carried by source 1 whatever the reflection does."""
    return (2 * p == n or 2 * q == n) and \
        4 * g.mean() * abs(p - q) == 2 * g.sum()


def locus_profile(n, total=1.0):
    """Anti-palindromic and dyadic; each mirror pair sums to `total`, so
    gbar = total/2 exactly. `total` is swept so the gbar factor in source 2's
    constant is falsifiable: at a single gbar the formula cannot be told from
    2*|p-q|*total/2 with the gbar dropped."""
    g = np.empty(n)
    for l in range(n // 2):
        g[l] = total * (l + 1) / 16.0
        g[n - 1 - l] = total - total * (l + 1) / 16.0
    if n % 2 == 1:
        g[n // 2] = total / 2
    return g


def zero_mean_profile(n):
    """Anti-palindromic with gbar = 0. NOT a Lindblad profile: gbar = 0 with
    gamma >= 0 forces gamma = 0, so half the sites carry gain and the generator is
    not CP. Everything read on this stratum is algebraic, not physical."""
    return np.array([(n - 1 - 2 * l) / 8.0 for l in range(n)])


def free_profile(n):
    """Dyadic and OFF the locus. T1: an arithmetic ramp will not do, being
    automatically anti-palindromic. Squares break it."""
    return np.array([(l + 1) ** 2 for l in range(n)]) / 16.0


def on_locus(g):
    n, gb = len(g), g.mean()
    return max(abs(g[l] + g[n - 1 - l] - 2 * gb) for l in range(n)) == 0.0


comp_bra = lambda n: (lambda a, b: (a, ((1 << n) - 1) ^ b))
comp_ket = lambda n: (lambda a, b: (((1 << n) - 1) ^ a, b))


def reflection(n):
    rev = lambda x: int(format(x, f"0{n}b")[::-1], 2)
    return lambda a, b: (rev(a), rev(b))


def transpose_reflection(n):
    """tauQ, F140's cell mirror: the transpose dressed with the reflection."""
    rev = lambda x: int(format(x, f"0{n}b")[::-1], 2)
    return lambda a, b: (rev(b), rev(a))


fails = []
log("=" * 78)
log("THE TWO SOURCES OF F153, AS ONE ANTILINEAR ENTRY-WISE IDENTITY")
log("=" * 78)
log()
log("  (*)  L[m i, m j] = -conj(L[i,j]) - c*delta_ij     compared with == 0.0")
log()

# ------------------------------------------------------------------ source 1
log("SOURCE 1  complement of a leg whose index is N/2, c = 2*sigma. Both legs are")
log("  swept. It must hold OFF the locus and on NON-pinned blocks.")
log()
log(f"  {'N':>2} {'p':>2} {'q':>2} {'dim':>5} {'leg':>4} {'D':>4} {'J':>5} "
    f"{'prof':>6} {'pinned':>7} {'(*)':>6} {'spec':>9} {'lin+gauge':>10}")
for n in (4, 6):
    for tag, g in (("locus", locus_profile(n)), ("free", free_profile(n))):
        assert on_locus(g) == (tag == "locus"), "T1: a profile is misclassified"
        sigma = g.sum()
        for delta in (0.0, 1.0, 2.0):
            for Jc in (1.25, 1.3):
                L = liouvillian(n, Jc, g, delta)
                probe = [((p, n // 2), "bra") for p in range(n + 1)] + \
                        [((n // 2, q), "ket") for q in range(n + 1) if q != n // 2]
                for (p, q), leg in probe:
                    if len(cells(n, p, q)) < 2:
                        continue
                    sub, cl = block_of(L, n, p, q)
                    img = comp_bra(n) if leg == "bra" else comp_ket(n)
                    r = antilinear_residual(sub, cl, img, 2 * sigma)
                    sp = spectral_pairing(sub, 2 * sigma)
                    lin = linear_residual_with_gauge(sub, cl, img, 2 * sigma, n)
                    if r != 0.0:
                        fails.append(f"source 1 not exact: N={n} ({p},{q}) {leg} "
                                     f"D={delta} J={Jc} {tag} r={r:.3e}")
                    if sp > 1e-9:
                        fails.append(f"source 1 spectral pairing broken: N={n} "
                                     f"({p},{q}) {leg} D={delta} J={Jc} {tag} "
                                     f"ratio={sp:.3e}")
                    if delta == 0.0 and Jc == 1.25 and lin != 0.0:
                        fails.append(f"the linear+gauge leg must be exact at "
                                     f"Delta=0: N={n} ({p},{q}) {lin:.3e}")
                    if delta != 0.0 and lin == 0.0:
                        fails.append(f"the linear+gauge leg must FAIL at "
                                     f"Delta!=0: N={n} ({p},{q}) D={delta}")
                    log(f"  {n:>2} {p:>2} {q:>2} {len(cl):>5} {leg:>4} {delta:>4} "
                        f"{Jc:>5} {tag:>6} {str(pinned(n,p,q)):>7} "
                        f"{str(r == 0.0):>6} {sp:>9.1e} {lin:>10.3e}")

log()
log("  CONTROL, source 1 needs [H, X^N] = 0: a LONGITUDINAL FIELD must break it.")
log("  F153's CRITERION survives that field, so the control also shows what (*) is")
log("  NOT: it is not the mechanism of the criterion on the h axis.")
n = 4
log(f"    {'h':>22}  {'(*) on (1,2)':>13}  {'criterion on the pinned blocks':>32}")
for hname, h in (("none", [0.0] * n),
                 ("generic", [0.5, 0.25, 0.75, 0.125]),
                 ("reflection-symmetric", [0.5, 0.25, 0.25, 0.5])):
    gu = np.full(n, 0.5)                       # uniform, so the criterion applies
    Lh = liouvillian(n, 1.25, gu, 0.0, h=h)
    sub, cl = block_of(Lh, n, 1, 2)
    r = antilinear_residual(sub, cl, comp_bra(n), 2 * gu.sum())
    worst_floor = 0.0
    for (p, q) in [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4)]:
        s2, _ = block_of(Lh, n, p, q)
        re = np.linalg.eigvals(s2).real
        worst_floor = max(worst_floor, abs(re - (-2 * 0.5 * abs(p - q))).max())
    log(f"    {hname:>22}  {r:>13.3f}  {worst_floor:>32.2e}")
    if (hname == "none") != (r == 0.0):
        fails.append(f"source 1 longitudinal control wrong for h={hname}: {r:.3e}")
    if worst_floor > 1e-9:
        fails.append(f"F153's criterion should survive h={hname}: {worst_floor:.2e}")

# ------------------------------------------------------------------ source 2
log()
log("SOURCE 2  site reflection on both legs, on the R90 locus, c = 4*gbar*|p-q|.")
log("  TWO values of gbar are swept, so the gbar factor is falsifiable. Rows marked")
log("  T2 are the blocks where the two constants coincide identically and the")
log("  spectral column therefore has no power; their entry-wise verdict still counts.")
log()
log(f"  {'N':>2} {'p':>2} {'q':>2} {'dim':>5} {'D':>4} {'J':>5} {'gbar':>5} "
    f"{'prof':>6} {'window':>12} {'pin':>5} {'(*)':>6} {'spec':>9} {'best c scan':>16}")
read_rows = []
for n in (4, 5, 6):
    for total in (1.0, 1.5):
        profiles = [("locus", locus_profile(n, total))]
        # free_profile ignores `total`, so running it under both would print the
        # same rows twice, byte for byte, and a duplicated row reads as corroboration
        if total == 1.0:
            profiles.append(("free", free_profile(n)))
        for tag, g in profiles:
            assert on_locus(g) == (tag == "locus"), "T1: a profile is misclassified"
            gbar, sigma = g.mean(), g.sum()
            assert gbar != 0.0
            span = 4 * max(2 * sigma, 4 * gbar * n)
            for delta in (0.0, 1.0):
                for Jc in (1.25, 1.3):
                    L = liouvillian(n, Jc, g, delta)
                    probe = [(0, q) for q in range(1, n + 1)] + \
                            [(p, n) for p in range(0, n)] + [(1, 1), (1, 2), (2, 2)]
                    for (p, q) in probe:
                        if len(cells(n, p, q)) < 2:
                            continue
                        sub, cl = block_of(L, n, p, q)
                        c = 4 * gbar * abs(p - q)
                        r = antilinear_residual(sub, cl, reflection(n), c)
                        sp = spectral_pairing(sub, c)
                        expect = (tag == "locus") and pinned(n, p, q)
                        t2 = constants_collide(n, p, q, g)
                        scan = ""
                        if not expect and Jc == 1.25 and delta == 0.0:
                            b, bc, edge = best_constant(sub, cl, reflection(n), span)
                            scan = f"{b:.2f}@{bc:+.2f}{'!' if edge else ''}"
                            if b == 0.0:
                                fails.append(f"a constant DOES serve a control "
                                             f"block: N={n} ({p},{q}) {tag} c={bc}")
                            if edge:
                                fails.append(f"the c scan hit its boundary at "
                                             f"N={n} ({p},{q}); widen the span")
                        # the entry-wise verdict is gated at EVERY coupling; where a
                        # rounding survives it is BOUNDED, never merely printed
                        bound = 4 * EPS * np.abs(L).max()
                        if expect and r != 0.0:
                            read_rows.append((n, p, q, delta, Jc, r))
                            if r > bound:
                                fails.append(f"source 2 residual above the error "
                                             f"model: N={n} ({p},{q}) D={delta} "
                                             f"J={Jc} r={r:.3e} > {bound:.3e}")
                        if not expect and r == 0.0:
                            fails.append(f"source 2 should NOT hold: N={n} "
                                         f"({p},{q}) D={delta} J={Jc} {tag}")
                        if not t2 and sp > 1e-9 and expect:
                            fails.append(f"source 2 spectral pairing broken: N={n} "
                                         f"({p},{q}) D={delta} J={Jc} {tag}")
                        log(f"  {n:>2} {p:>2} {q:>2} {len(cl):>5} {delta:>4} "
                            f"{Jc:>5} {gbar:>5.2f} {tag:>6} "
                            f"{str(window(n,p,q)):>12} "
                            f"{('T2' if t2 else str(pinned(n,p,q))):>5} "
                            f"{str(r == 0.0):>6} {sp:>9.1e} {scan:>16}")

log()
log("  CONTROL, source 2 needs H REFLECTION-SYMMETRIC: a non-palindromic BOND LIST")
log("  must break it, while source 1 stays exactly 0.0 on the same chain. Without")
log("  this row the whole source-2 section runs at uniform J and cannot see the")
log("  hypothesis it depends on; source 2's every other row is uniform J.")
log(f"    {'N':>2} {'bonds':>22} {'palindrome':>11} {'source 2':>10} {'source 1':>10}")
for n in (4, 5, 6):
    g = locus_profile(n)
    gbar, sigma = g.mean(), g.sum()
    for Js in ([1.0] * (n - 1), [1.0 + 0.25 * b for b in range(n - 1)]):
        pal = Js == Js[::-1]
        L = liouvillian(n, Js, g, 1.0)
        # source 2 on the pinned blocks, where it is supposed to hold
        s2 = max(antilinear_residual(*block_of(L, n, 0, q),
                                     image=reflection(n), c=4 * gbar * q)
                 for q in range(1, n + 1))
        # source 1 on an index-N/2 block, blind to the bond profile either way
        s1 = (antilinear_residual(*block_of(L, n, 0, n // 2),
                                  image=comp_bra(n), c=2 * sigma)
              if n % 2 == 0 else float("nan"))
        log(f"    {n:>2} {str([round(j, 2) for j in Js]):>22} {str(pal):>11} "
            f"{s2:>10.3f} {s1:>10.3f}")
        if pal and s2 != 0.0:
            fails.append(f"source 2 broke on a SYMMETRIC bond list at N={n}")
        if not pal and s2 == 0.0:
            fails.append(f"source 2 survived an ASYMMETRIC bond list at N={n}; "
                         f"the reflection-symmetric-H hypothesis is not being tested")
        if n % 2 == 0 and s1 != 0.0:
            fails.append(f"source 1 is not bond-blind at N={n}: {s1:.3e}")

log()
log("  THE SURVIVING ROUNDING, bounded and then READ. The mechanism is the ORDER in")
log("  which the ZZ diagonal accumulates over bonds, not the coupling itself: J = pi")
log("  and J = 1.1 are just as unrepresentable as 1.3 and give exactly 0.0.")
for n in (5, 6):
    g = locus_profile(n)
    row = []
    for Jc in (1.25, 1.1, np.pi, 1.3, 0.9, 2.7):
        L = liouvillian(n, Jc, g, 1.0)
        sub, cl = block_of(L, n, 0, 1)
        row.append((Jc, antilinear_residual(sub, cl, reflection(n),
                                            4 * g.mean())))
    log(f"    N={n} (0,1) Delta=1: " +
        "  ".join(f"J={j:.4g}:{r:.1e}" for j, r in row))
# the N=4 claim is COMPUTED, not printed
n4_all_zero = True
for Jc in (1.25, 1.1, np.pi, 1.3, 0.9, 2.7, 1 / 3, 1.7):
    g = locus_profile(4)
    L = liouvillian(4, Jc, g, 1.0)
    for (p, q) in [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4)]:
        sub, cl = block_of(L, 4, p, q)
        if antilinear_residual(sub, cl, reflection(4), 4 * g.mean() * abs(p - q)) != 0.0:
            n4_all_zero = False
log(f"    N = 4 exactly 0.0 over eight couplings, computed: {n4_all_zero}")
if not n4_all_zero:
    fails.append("the N=4 exactness claim is false")

# ------------------------------------------------------------- the zero stratum
log()
log("  THE ZERO-MEAN STRATUM, gbar = 0: every class wants centre 0 at once, so EVERY")
log("  block is served with c = 0. NOT a Lindblad profile (gbar = 0 with gamma >= 0")
log("  forces gamma = 0), so this row is algebra about the generator, not physics.")
for n in (4, 5, 6):
    g = zero_mean_profile(n)
    assert on_locus(g) and g.mean() == 0.0 and g.min() < 0.0
    for delta in (0.0, 1.0):
        L = liouvillian(n, 1.25, g, delta)
        tot = served = nonpinned = 0
        for p in range(n + 1):
            for q in range(n + 1):
                if len(cells(n, p, q)) < 2:
                    continue
                sub, cl = block_of(L, n, p, q)
                tot += 1
                if antilinear_residual(sub, cl, reflection(n), 0.0) == 0.0:
                    served += 1
                    nonpinned += not pinned(n, p, q)
        log(f"    N={n} Delta={delta}: {served}/{tot} blocks exact at c = 0, "
            f"{nonpinned} of them NOT pinned")
        if served != tot:
            fails.append(f"zero-mean stratum: only {served}/{tot} at N={n}")
log("  The stratum inherits source 2's bond hypothesis, and it is OURS to state:")
log("  section 5 owns the observation, but its own entry-wise check runs tauQ on")
log("  q = p and q = N-p only. A non-palindromic bond list must break this too.")
for n in (4, 5):
    g = zero_mean_profile(n)
    L = liouvillian(n, [1.0 + 0.25 * b for b in range(n - 1)], g, 1.0)
    worst = max(antilinear_residual(*block_of(L, n, p, q), image=reflection(n), c=0.0)
                for p in range(n + 1) for q in range(n + 1)
                if len(cells(n, p, q)) > 1)
    log(f"    N={n} ASYMMETRIC bonds: worst over all dim>1 blocks {worst:.3f}")
    if worst == 0.0:
        fails.append(f"the gbar=0 stratum survived an asymmetric bond list at "
                     f"N={n}; its H hypothesis is not being tested")

# ---------------------------- the F140 defect, as a MATRIX and not as a number
log()
log("  THE F140 DEFECT AS A MATRIX. On the (1,1)-type corner the plain-c version of")
log("  (*) fails, and WHERE it fails is the content: the whole residual sits on the")
log("  DIAGONAL cells at exactly 8*gbar, is exactly 0.0 away from the diagonal of")
log("  the block, and on the off-diagonal cells is zero within the error model (the")
log("  'on O' column below is 0.0 at N = 4 and one ZZ rounding at N = 5, 6, not the")
log("  exact 0.0 an earlier version of this sentence claimed). So the size-class")
log("  mismatch is the same defect MATRIX 8*gbar*P_D that PROOF_R90_FROZEN_DIVISOR")
log("  section 2 localises for the LINEAR tauQ identity, reached antilinearly.")
log(f"    {'N':>2} {'J':>5} {'8*gbar':>7} {'on D':>8} {'on O':>9} {'off-diag':>9} "
    f"{'with P_D':>9}")
for n in (4, 5, 6):
    for Jc in (1.3, 0.7):
        g = locus_profile(n)
        gbar = g.mean()
        sub, cl = block_of(liouvillian(n, Jc, g, 1.0), n, 1, 1)
        m = permutation_of(cl, reflection(n))
        base = sub[np.ix_(m, m)] + np.conj(sub)
        plain = base + 8 * gbar * np.eye(len(cl))
        isD = np.array([1.0 if a == b else 0.0 for (a, b) in cl])
        onD = np.abs(np.diag(plain)[isD == 1.0]).max()
        onO = np.abs(np.diag(plain)[isD == 0.0]).max()
        off = np.abs(plain - np.diag(np.diag(plain))).max()
        withPD = np.abs(plain - 8 * gbar * np.diag(isD)).max()
        log(f"    {n:>2} {Jc:>5} {8*gbar:>7.2f} {onD:>8.3f} {onO:>9.1e} "
            f"{off:>9.1e} {withPD:>9.1e}")
        bound = 4 * EPS * np.abs(sub).max()
        if abs(onD - 8 * gbar) > 1e-12:
            fails.append(f"the D-cell residual is not 8*gbar at N={n} J={Jc}: {onD}")
        if off > bound or onO > bound:
            fails.append(f"the F140 defect is not confined to the D cells at "
                         f"N={n} J={Jc}: off={off:.3e} onO={onO:.3e}")
        if withPD > bound:
            fails.append(f"the defect-matrix identity failed at N={n} J={Jc}: "
                         f"{withPD:.3e} > {bound:.3e}")

# ------------------------------------------------- tauQ against the bare reflection
log()
log("  tauQ AGAINST THE BARE REFLECTION. F140's cell mirror is the transpose dressed")
log("  with the reflection; on the RATE DIAGONAL the transpose is invisible, since a")
log("  disagreement set is symmetric in its two legs: rate(a,b) = rate(b,a) by")
log("  arithmetic, so the two maps agree there as an IDENTITY, not a measurement.")
log("  An earlier version dressed that identity as a gated row; its failure branch")
log("  compared x^y with y^x and no input could reach it (the T3 shape). What lets")
log("  the two results share one reflection action is this identity; what the gate")
log("  MEASURES is where the maps part, below.")
log("  AND WHERE THEY PART, which is more than the hop part. At Delta != 0 the")
log("  transpose also moves the ZZ IMAGINARY diagonal; and the primary difference")
log("  is not a part of L at all but the KIND of identity each map satisfies:")
log("  tauQ's is LINEAR, the bare reflection's is ANTILINEAR, and on the two-class")
log("  corner each satisfies exactly one and fails the other.")
log("  Both identities are written WITH the 8*gbar*P_D defect, so each map is asked")
log("  the same question and only the linearity differs.")
log(f"    {'N':>2} {'map':>11} {'Im ZZ diag diff':>16} {'linear id':>10} "
    f"{'antilinear id':>14}")
for n in (4, 5):
    g = locus_profile(n)
    gbar = g.mean()
    sub, cl = block_of(liouvillian(n, 1.3, g, 1.0), n, 1, 1)
    isD = np.array([1.0 if a == b else 0.0 for (a, b) in cl])
    defect = 8 * gbar * (np.eye(len(cl)) - np.diag(isD))
    for name, img in (("tauQ", transpose_reflection(n)), ("reflection", reflection(n))):
        m = permutation_of(cl, img)
        imdiff = np.abs(np.diag(sub)[m].imag - np.diag(sub).imag).max()
        lin = np.abs(sub[np.ix_(m, m)] + sub + defect).max()
        anti = np.abs(sub[np.ix_(m, m)] + np.conj(sub) + defect).max()
        log(f"    {n:>2} {name:>11} {imdiff:>16.3f} {lin:>10.3f} {anti:>14.3f}")
        bound = 4 * EPS * np.abs(sub).max()
        want_lin = (name == "tauQ")
        if want_lin and lin > bound:
            fails.append(f"tauQ's LINEAR identity failed at N={n}: {lin:.3e}")
        if not want_lin and anti > bound:
            fails.append(f"the reflection's ANTILINEAR identity failed at "
                         f"N={n}: {anti:.3e}")
        if (anti if want_lin else lin) <= bound:
            fails.append(f"{name} satisfies BOTH identities at N={n}; the "
                         f"linear/antilinear split is not being measured")

# ------------------------------------ the palindrome extends, the flat set does not
log()
log("WHAT THE ENLARGED SET DOES NOT BUY. The pinned index-N/2 blocks flatten onto")
log("the trace constant at EVERY Delta; the non-pinned ones do something that")
log("DEPENDS ON Delta, and no mechanism for that dependence is claimed here.")
log("Two earlier readings of this table were wrong and both are gated against now.")
log("  (a) 'the saturation width is twice the window HALF-width' is the width in the")
log("      RATE variable and agrees with the lambda one only at gbar = 1/2, which is")
log("      where the table used to be run. gbar is swept.")
log("  (b) 'a multi-class block saturates on the size-class interval and NEVER")
log("      flattens' is a Delta = 0 and 1 reading. Delta is swept, and at Delta = 2")
log("      the blocks (1,3) and (2,3) flatten exactly onto the trace constant.")
log("What IS measured, and asserted per Delta rather than as a law: where a block")
log("saturates, it saturates on the interval spanned by its size-class centres")
log("-2*gbar*|S|, ENDPOINTS INCLUDED, which a width alone would not show.")
n = 6
# the ladder must REACH the strong-coupling limit at every Delta, and 100 does not:
# at Delta = 2 the spread is still 4.29 at J = 100, 3.94 at 1e3 and exactly 0 from
# 1e4 on. An earlier ladder stopped at 100 and read a converged phenomenon as broken.
couplings = (1.25, 100.0, 10000.0)
tracec = lambda p, q, gbar: -2 * gbar * (p + q - 2 * p * q / n)
# what the CONVERGED table says, block family by block family; the label is computed
# from the run and compared with this, so a drift in either direction fails
EXPECT = {0.0: "sat@centres", 0.5: "NEITHER", 1.0: "sat@centres", 2.0: "FLAT@trace"}
log(f"  {'gbar':>5} {'D':>4} {'block':>8} {'window':>12} {'pin':>5} " +
    "".join(f"{'J=' + str(j):>11}" for j in couplings) +
    f"{'reading':>12} {'as expected':>12}")
for total in (1.0, 1.5):
    g = locus_profile(n, total)
    gbar = g.mean()
    for delta in (0.0, 0.5, 1.0, 2.0):
        for p in range(n + 1):
            q = n // 2
            if len(cells(n, p, q)) < 2:
                continue
            w = window(n, p, q)
            row, radii, ends = [], [], None
            for J in couplings:
                sub, _ = block_of(liouvillian(n, J, g, delta), n, p, q)
                re = np.linalg.eigvals(sub).real
                row.append(re.max() - re.min())
                radii.append(max(np.abs(np.linalg.eigvals(sub)).max(), 1.0))
                if J == couplings[-1]:
                    ends = (re.max(), re.min())
            spread = row[-1] / radii[-1]
            hi, lo = -2 * gbar * w[0], -2 * gbar * w[-1]
            tc = tracec(p, q, gbar)
            if spread < 1e-9 and abs(ends[0] - tc) < 1e-6 * max(1.0, abs(tc)):
                tag = "FLAT@trace"
            elif (abs(row[-1] - 2 * gbar * (w[-1] - w[0])) < 1e-3 and
                  max(abs(ends[0] - hi), abs(ends[1] - lo)) < 1e-3 * max(1.0, abs(lo))):
                tag = "sat@centres"
            else:
                tag = "NEITHER"
            # the pinned blocks and (N/2,N/2) are Delta-robust; the rest follow EXPECT
            want = ("FLAT@trace" if pinned(n, p, q)
                    else "sat@centres" if p == q else EXPECT[delta])
            ok = (tag == want)
            log(f"  {gbar:>5.2f} {delta:>4} {'(' + str(p) + ',' + str(q) + ')':>8} "
                f"{str(w):>12} {str(pinned(n,p,q)):>5} " +
                "".join(f"{v:>11.5f}" for v in row) +
                f"{tag:>12} {str(ok):>12}")
            if not ok:
                fails.append(f"the converged reading changed at gbar={gbar} "
                             f"D={delta} ({p},{q}): got {tag}, the committed table "
                             f"says {want} (ends={ends}, centres={(hi, lo)}, "
                             f"trace={tc})")

# ------------------------------------- flatness at uniform gamma needs no (*)
log()
log("  AND THE CRITERION ITSELF NEEDS NEITHER (*) NOR THE LOCUS. At UNIFORM gamma a")
log("  pinned block has one size class, so its dissipator is the real scalar")
log("  c = -2*gamma*|p-q| and L = c*Id + i*S with S EXACTLY Hermitian. That is")
log("  F153's own derivation and it puts every Re lambda at c with no involution in")
log("  it. Worst non-Hermiticity of (L - c)/i over all pinned blocks, and worst")
log("  |Re lambda - c|, with and without a longitudinal field:")
for n_ in (4, 5):
    for hh in (None, [0.3 * (s + 1) for s in range(n_)]):
        gu = np.full(n_, 0.25)
        L = liouvillian(n_, 1.3, gu, 1.0, h=hh)
        wh = wr = 0.0
        for (p, q) in [(p, q) for p in range(n_ + 1) for q in range(n_ + 1)
                       if pinned(n_, p, q) and len(cells(n_, p, q)) > 1]:
            sub, _ = block_of(L, n_, p, q)
            c = -2 * 0.25 * abs(p - q)
            S = (sub - c * np.eye(sub.shape[0])) / 1j
            wh = max(wh, np.abs(S - S.conj().T).max())
            wr = max(wr, np.abs(np.linalg.eigvals(sub).real - c).max())
        log(f"    N={n_} field={'yes' if hh else 'no ':>3}: "
            f"non-Hermiticity {wh:.1e}, worst |Re - c| {wr:.1e}")
        if wh != 0.0:
            fails.append(f"S is not exactly Hermitian at N={n_}: {wh:.3e}")
        if wr > 1e-13:
            fails.append(f"a pinned block left its floor at uniform gamma, "
                         f"N={n_}: {wr:.3e}")

# ------------------------------------------------- the window is the criterion
log()
log("THE WINDOW COLLAPSE IS THE PINNING CRITERION (pure combinatorics; the repo owns")
log("this three times over, PROOF_CODIM1 section 6, DisagreementWindow,")
log("WindowShellLemmaTests, and it is restated here only to be used)")
bad = [(n, p, q) for n in range(2, 9) for p in range(n + 1) for q in range(n + 1)
       if (len(window(n, p, q)) == 1) != pinned(n, p, q)]
log(f"  single size class <=> pinned, N = 2..8, all (p,q): "
    f"{'no counterexample' if not bad else bad}")
if bad:
    fails.append(f"window criterion broken at {bad}")
for n in (4, 6, 8):
    log(f"  N={n}: F140's (1,1)-type corner block has window {window(n,1,1)}, two "
        f"classes, so at gbar != 0 no single centre serves it")

log()
log("=" * 78)
if read_rows:
    # every qualifier here is COMPUTED from the rows, not asserted in prose: an
    # earlier version printed a fixed sentence beside a bucket that could not fail
    all_aniso = all(d != 0.0 for (_, _, _, d, _, _) in read_rows)
    none_at_n4 = all(n != 4 for (n, _, _, _, _, _) in read_rows)
    worst = max(r for (_, _, _, _, _, r) in read_rows)
    log(f"READ, and bounded by 4*eps*max|L|: {len(read_rows)} rows carry one")
    log(f"rounding of the ZZ diagonal. All at Delta != 0: {all_aniso}. "
        f"None at N = 4: {none_at_n4}. Worst: {worst:.3e}.")
    if not (all_aniso and none_at_n4):
        fails.append("the read bucket's own description does not match its rows")
    for (n, p, q, d, j, r) in read_rows[:4]:
        log(f"  N={n} ({p},{q}) Delta={d} J={j}: {r:.3e}")
else:
    log("READ: no row carried a surviving rounding.")
log()
log(f"VERDICT: {'all checks pass' if not fails else str(len(fails)) + ' FAILURES'}")
for x in fails:
    log("  " + x)
log("=" * 78)
