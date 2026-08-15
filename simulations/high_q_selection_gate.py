"""The high-Q selection gate: flat is generic, spread is a resonance.

Companion gate of experiments/THE_SPREAD_IS_A_RESONANCE.md. Every number that
note asserts is computed here and compared against its committed value; a drift
either way fails. Writes simulations/results/high_q_selection_gate.txt; the
last line is the VERDICT line, and a file whose last line is not one was cut
off, not green.

WHAT THIS PINS, in five sections.

  (A) THE RULE, scope-extended. F122's registered Mechanism (high-Q degenerate
      PT, PROOF_STRUCTURAL_CEILING section 1) says: for each ad_H-eigenspace
      Omega, the strong-coupling rates are the eigenvalues of P_Omega D P_Omega.
      The repo gated it for XY (Delta = 0), uniform gamma, slowest mode only.
      Here it is gated on XXZ, a gamma PROFILE, the FULL Re-spectrum, the four
      distinct (p, 3) blocks at N = 6 (the mirrors (4,3)..(6,3) repeat value
      for value, the predecessor's record) and Delta in {0, 0.5, 1, 2}. The error model has
      TWO terms and both are laws, not pass marks:
        truncation:  c_t / J^2 with c_t <= 3e3 measured. The 1/J term is ABSENT
                     because A is anti-Hermitian: with omega_k the eigenvalues
                     OF A (purely imaginary; NOT the real w of eigh(iA) below),
                     the correction sum_m |D_mk|^2 / (J*(omega_k - omega_m))
                     has real numerators over imaginary denominators and is
                     purely imaginary; on a degenerate eigenspace the effective
                     operator sum_m Pi D|m><m|D Pi / (J*(omega - omega_m)) is
                     anti-Hermitian by the same split. Gated as a RATIO: every
                     row in the truncation window at J = 1e3 must shrink by
                     30..300x to J = 1e4 (the pure-1/J^2 value is 100).
        floor:       c_f * eps * J * ||A||_2 from the dense eigensolver's own
                     backward error (the reference, not the rule, degrades as
                     J grows), c_f <= 32 with measured ratios 1..11.
      Rows below the rule's validity scale J* ~ ||D|| / (min difference-gap)
      at J = 1e3 (the Delta = 2 rows) are printed and labeled, not gated: the
      predecessor's J-ladder trap IS this scale.
      Scope note: (A) certifies that degenerate PT reproduces the same L that
      block_L builds (prediction and reference share the builder); the physics
      is pinned by the hard-coded constants of (B)..(E).

  (B) FLAT IS GENERIC. On a 24-point Delta grid the predicted spread of (1,3)
      and (2,3) is zero except exactly at Delta in {0, 0.5, 1}; at those three
      it equals the committed values 2.0000 / 0.4494 / 2.0000 for (1,3) and
      4.0000 / 0.5663 / 4.0000 for (2,3). The committed four-point table of
      WHAT_THE_R90_LOCUS_BUYS.md sampled the three resonant Deltas plus one
      generic one. The mirror Delta -> -Delta (the sublattice rotation
      U = prod_{l odd} Z_l gives U H(Delta) U = -H(-Delta)) is spot-checked at
      Delta = -0.25 (flat), -0.5 (0.4494) and -1 (2.0000). The gamma-bar
      proportionality (A is H-only, D linear in gamma) is gated exactly:
      pred(1.5*gamma) == 1.5*pred(gamma) to 1e-12. The flatness threshold 1e-6
      is generous against the flat residual, which is eigenvector-conditioning
      noise ~ eps*||A||*||D||/mingap: it reads ~1e-8 at the grid edge
      Delta = 4 and would climb through 1e-6 near Delta ~ 6 as the difference
      spectrum crowds, so extending the grid needs the threshold rescaled by
      the measured min gap, not a wider net.

  (C) THE CANCELLATION, and it needs no locus. [H, X^N] = 0 exactly (compared
      == 0.0 entry-wise), and X^N maps the half-filling sector to itself, so
      every NON-degenerate half-filling eigenvector phi has
      |phi(a)|^2 = |phi(a-complement)|^2 and site density exactly 1/2. In the
      pair average sum_l gamma_l (m_l + n_l - 2 m_l n_l) the density n = 1/2
      cancels m COMPLETELY: the average is sigma/2, the mode sits at exactly
      -sigma, for ANY profile. The density check runs at Delta = 2 (all levels
      simple) AND Delta = 0.5 (where the 5/2 doublet must be exempted and the
      simple ones still read 1/2); flat-at--sigma is gated on (0,3), (1,3) and
      (2,3), on BOTH profiles, at three generic Delta.

  (D) THE RESONANCE CATALOGUE. Spread REQUIRES a collision of the difference
      spectrum {E_i - F_j} (necessary, not sufficient: on the locus at
      Delta = 0 the pinned (0,3) collapses although collisions are present).
      Gated:
        Delta = 1    SU(2): all 6 sector-1 and all 15 sector-2 levels reappear
                     exactly in sector 3 (0 of each at Delta = 2);
        Delta = 0.5  the half-filling multiplicity HISTOGRAM is 6x1 at N = 4,
                     18x1 + 1x2 at N = 6, 54x1 + 8x2 at N = 8 (a histogram,
                     not a pair count, so a triple cannot masquerade); the
                     N = 6 doublet at E = 5J/2 is EXACT over Z (Faddeev-
                     LeVerrier charpoly of the integer matrix 2H, built from
                     bit logic in exact ints and compared == to the float
                     route; p(5) = 0, p'(5) = 0, p''(5) != 0); the shared-level
                     census at Delta = 0.5 is 0 / 0 / 8 at N = 4 / 6 / 8, so
                     at N = 8 the resonance is CROSS-sector (the eight doublets
                     are the eight sector-1 levels), unlike N = 6's
                     within-sector doublet;
        (3,3)        ker(ad_H) on a (p,p) block is the commutant, dimension
                     exactly sum(m_i^2) over the sector multiplicities (>= 20
                     for ANY H, so the >= is not gateable; the runtime equality
                     alone is not either, both sides moving with H); gated as
                     the runtime equality AND the hard-coded values 22 / 20 /
                     20 at Delta = 0.5 / 1 / 2, the resonance count's pin.

  (E) THE FAMILY TIES AND THE FOURTH BEHAVIOUR. The rule reproduces F153's
      whole off-locus (0, N/2) family: at Delta = 0.5 the spreads 0.0000 /
      0.6456 / 1.2903 at N = 4 / 6 / 8 and at Delta = 0 the spreads 0.2236 /
      0.5325 / 1.0756 (F153 quotes the first five; 1.0756 is this rule's
      prediction, measured dense 2026-08-15). On a block with BOTH indices off
      N/2 the cancellation is absent: (1,2) at N = 6, Delta = 2 is a stable
      density-overlap band, spread 0.9982, range [-2.8510, -1.8528], the same
      at J = 1e4 and 1e6.

EXACTNESS POLICY. [H, X^N], the integer charpoly, its float-int comparison and
the gamma-bar proportionality are exact routes, compared with == 0. Everything
through an eigensolver states its tolerance next to the check; the rule check
(A) gates its two-term error model as scaling laws.

Pauli book throughout: H = J * sum_bonds (XX + YY + Delta*ZZ), site l = bit l,
cell |a><b| rate 2 * sum_{l: a_l != b_l} gamma_l, sigma = sum gamma_l.
"""
import numpy as np
from collections import Counter
from fractions import Fraction

OUT = "simulations/results/high_q_selection_gate.txt"
lines = []
fails = []

def log(s=""):
    lines.append(s)
    print(s)

def popcount(x):
    return bin(x).count("1")

# ----------------------------------------------------------------- primitives
def build_H(n, J, D):
    d = 1 << n
    H = np.zeros((d, d))
    for a in range(d):
        for l in range(n - 1):
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a, a] += J * D * za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2 * J
    return H

def block_L(n, p, q, J, D, gamma):
    """The (p,q) joint-popcount block of L, cells (a,b) ket-major."""
    H = build_H(n, J, D)
    kets = [a for a in range(1 << n) if popcount(a) == p]
    bras = [b for b in range(1 << n) if popcount(b) == q]
    cells = [(a, b) for a in kets for b in bras]
    idx = {c: i for i, c in enumerate(cells)}
    m = len(cells)
    L = np.zeros((m, m), dtype=complex)
    for (a, b), i in idx.items():
        for c in kets:
            if H[a, c] != 0:
                L[i, idx[(c, b)]] += -1j * H[a, c]
        for d_ in bras:
            if H[b, d_] != 0:
                L[i, idx[(a, d_)]] += 1j * H[b, d_]
        rate = sum(g for l, g in enumerate(gamma) if ((a >> l) & 1) != ((b >> l) & 1))
        L[i, i] += -2 * rate
    return L, cells

def split_A_D(n, p, q, D, gamma):
    """L(J) = J*A + Dm on the block: A the J = 1 Hamiltonian part, Dm diagonal."""
    LJ1, cells = block_L(n, p, q, 1.0, D, gamma)
    Dm, _ = block_L(n, p, q, 0.0, D, gamma)
    return LJ1 - Dm, Dm, cells

def predicted_res(n, p, q, D, gamma, gtol=1e-9, want_norm=False):
    """The high-Q degenerate-PT Re-spectrum: union over ad_H eigenspaces of
    spec(P_Omega Dm P_Omega) (Hermitian, so eigvalsh), sorted. gtol sits in the
    measured gap void [3e-14, 4.3e-6] of the difference spectrum on the grid."""
    A, Dm, _ = split_A_D(n, p, q, D, gamma)
    w, V = np.linalg.eigh(1j * A)
    order = np.argsort(w)
    w, V = w[order], V[:, order]
    res = []
    start = 0
    for i in range(1, len(w) + 1):
        if i == len(w) or abs(w[i] - w[start]) > gtol:
            P = V[:, start:i]
            res.extend(np.linalg.eigvalsh(P.conj().T @ Dm @ P))
            start = i
    out = np.sort(np.array(res))
    return (out, np.abs(w).max()) if want_norm else out

def measured_res(n, p, q, D, gamma, J):
    L, _ = block_L(n, p, q, J, D, gamma)
    return np.sort(np.linalg.eigvals(L).real)

def sector_spectrum(n, p, D, J=1.0, vectors=False):
    H = build_H(n, J, D)
    idx = [a for a in range(1 << n) if popcount(a) == p]
    Hs = H[np.ix_(idx, idx)]
    if vectors:
        E, V = np.linalg.eigh(Hs)
        return E, V, idx
    return np.linalg.eigvalsh(Hs)

def mult_histogram(E, tol=1e-9):
    """Multiplicity histogram of a sorted spectrum: {multiplicity: count}."""
    groups, start = [], 0
    for i in range(1, len(E) + 1):
        if i == len(E) or E[i] - E[start] > tol:
            groups.append(i - start)
            start = i
    return dict(sorted(Counter(groups).items()))

def locus_profile(n):
    if n % 2:
        raise ValueError("even N only")
    g = np.array([(l + 1) / 16 for l in range(n // 2)])
    return np.concatenate([g, (1 - g)[::-1]])

def free_profile(n):
    return np.array([(l + 1) ** 2 / 16 for l in range(n)])

def on_locus(g):
    return bool(np.all(g + g[::-1] == g[0] + g[-1]))

log("=" * 78)
log("THE HIGH-Q SELECTION GATE: flat is generic, spread is a resonance")
log("=" * 78)

n = 6
gl = locus_profile(6)
gf = free_profile(6)
log(f"locus profile   {np.array2string(gl, precision=4)}  on locus: {on_locus(gl)}  sigma = {gl.sum()}")
log(f"off-locus prof. {np.array2string(gf, precision=4)}  on locus: {on_locus(gf)}  sigma = {gf.sum()}")
if not on_locus(gl) or on_locus(gf):
    fails.append("profile locus classes wrong")
sigma_l, sigma_f = gl.sum(), gf.sum()
eps = np.finfo(float).eps

# ------------------------------------------------------------- (A) the rule
log()
log("(A) THE RULE. dev(J) = max|predicted - dense| per row. Bound: each dev")
log("    <= 3e3/J^2 + 32*eps*J*||A||; truncation-window rows (1e-7 < dev(1e3)")
log("    < 1e-2) must shrink 30..300x to J = 1e4 (pure 1/J^2 = 100x); rows")
log("    above the window at J = 1e3 are below the validity scale J*, printed")
log("    and not gated there.")
for D in (0.0, 0.5, 1.0, 2.0):
    for p in range(4):
        pred, normA = predicted_res(n, p, 3, D, gl, want_norm=True)
        devs = {J: np.max(np.abs(pred - measured_res(n, p, 3, D, gl, J)))
                for J in (1e3, 1e4, 1e6)}
        tag = ""
        for J in (1e4, 1e6):
            b = 3e3 / J ** 2 + 32 * eps * J * normA
            if devs[J] > b:
                fails.append(f"(A) dev over model at Delta={D} ({p},3) J={J:.0e}: "
                             f"{devs[J]:.2e} > {b:.2e}")
        if 1e-7 < devs[1e3] < 1e-2:
            r = devs[1e3] / devs[1e4]
            tag = f"  ratio(1e3/1e4)={r:.0f}"
            if not (30 <= r <= 300):
                fails.append(f"(A) 1/J^2 ratio off at Delta={D} ({p},3): {r:.0f}")
        elif devs[1e3] >= 1e-2:
            tag = "  [below J* at 1e3, not gated there]"
        log(f"    Delta={D:<4} ({p},3): dev(1e3)={devs[1e3]:.2e} dev(1e4)={devs[1e4]:.2e} "
            f"dev(1e6)={devs[1e6]:.2e} ||A||={normA:.1f}{tag}")

# ------------------------------------------------------ (B) flat is generic
log()
log("(B) FLAT IS GENERIC: predicted spread on the Delta grid (locus profile).")
log("    Zero (< 1e-6) off {0, 0.5, 1}; the committed values at the three.")
GRID = [0.0, 0.1, 0.2, 0.25, 1 / 3, 0.4, 0.5, 0.6, 2 / 3, 0.75, 0.8, 0.9, 1.0,
        1.1, 1.25, 4 / 3, 1.5, 5 / 3, 1.75, 2.0, 2.25, 2.5, 3.0, 4.0]
RESONANT = {0.0: {(1, 3): 2.0000, (2, 3): 4.0000},
            0.5: {(1, 3): 0.4494, (2, 3): 0.5663},
            1.0: {(1, 3): 2.0000, (2, 3): 4.0000}}
for D in GRID:
    sp = {}
    for (p, q) in ((1, 3), (2, 3)):
        r = predicted_res(n, p, q, D, gl)
        sp[(p, q)] = r.max() - r.min()
    tag = "resonant" if D in RESONANT else "generic"
    log(f"    Delta={D:<7.4f} spread (1,3)={sp[(1, 3)]:.4f}  (2,3)={sp[(2, 3)]:.4f}  [{tag}]")
    if D in RESONANT:
        for k, want in RESONANT[D].items():
            if abs(sp[k] - want) > 1e-4:
                fails.append(f"(B) resonant spread drifted at Delta={D} {k}: {sp[k]:.4f} != {want}")
    else:
        for k, v in sp.items():
            if v > 1e-6:
                fails.append(f"(B) generic Delta={D} not flat on {k}: spread {v:.2e}")
log("    the Delta -> -Delta mirror (sublattice rotation), spot-checked:")
for (D, want) in ((-0.25, 0.0), (-0.5, 0.4494), (-1.0, 2.0000)):
    r = predicted_res(n, 1, 3, D, gl)
    v = r.max() - r.min()
    log(f"      Delta={D}: (1,3) spread {v:.4f} (want {want})")
    if abs(v - want) > 1e-4:
        fails.append(f"(B) mirror value at Delta={D}: {v:.4f} != {want}")
log("    gamma-bar proportionality, exact: max|pred(1.5g) - 1.5*pred(g)|:")
for D in (0.5, 2.0):
    r1 = predicted_res(n, 1, 3, D, gl)
    r2 = predicted_res(n, 1, 3, D, 1.5 * gl)
    d = np.max(np.abs(r2 - 1.5 * r1))
    log(f"      Delta={D}: {d:.2e}")
    if d > 1e-12:
        fails.append(f"(B) gamma-bar proportionality broken at Delta={D}: {d:.2e}")

# ---------------------------------------------------- (C) the cancellation
log()
log("(C) THE CANCELLATION. [H, X^N] compared == 0.0 entry-wise; nondegenerate")
log("    half-filling densities at 1/2 (eigh backward error, tol 1e-10), at")
log("    Delta = 2 (all simple) and 0.5 (the doublet exempted); flat value")
log("    exactly -sigma on THREE blocks and BOTH profiles at generic Delta.")
for D in (0.5, 1.0, 2.0):
    H = build_H(6, 1.0, D)
    X = np.zeros_like(H)
    for a in range(64):
        X[a ^ 63, a] = 1.0
    comm = np.abs(H @ X - X @ H).max()
    log(f"    Delta={D}: max|[H, X^N]| = {comm}")
    if comm != 0.0:
        fails.append(f"(C) [H, X^N] != 0 at Delta={D}: {comm}")
for D in (0.5, 2.0):
    E, V, idx = sector_spectrum(6, 3, D, vectors=True)
    occ = np.array([[(a >> l) & 1 for l in range(6)] for a in idx], float)
    dens = (np.abs(V) ** 2).T @ occ
    simple = [i for i in range(len(E))
              if (i == 0 or E[i] - E[i - 1] > 1e-9) and (i == len(E) - 1 or E[i + 1] - E[i] > 1e-9)]
    dev = np.abs(dens[simple] - 0.5).max()
    devall = np.abs(dens - 0.5).max()
    log(f"    Delta={D}: {len(simple)}/{len(E)} half-filling levels simple; "
        f"max|n_l - 1/2| = {dev:.2e} on them ({devall:.2e} over all)")
    if dev > 1e-10:
        fails.append(f"(C) half-filling density off 1/2 at Delta={D}: {dev:.2e}")
for (gname, g, sig) in (("locus", gl, sigma_l), ("off-locus", gf, sigma_f)):
    for D in (0.75, 1.5, 2.0):
        worst_sp = worst_off = 0.0
        for (p, q) in ((0, 3), (1, 3), (2, 3)):
            r = predicted_res(6, p, q, D, g)
            worst_sp = max(worst_sp, r.max() - r.min())
            worst_off = max(worst_off, np.abs(r + sig).max())
        log(f"    {gname:<9} Delta={D:<4}: worst spread {worst_sp:.2e}, worst |Re + sigma| = {worst_off:.2e}"
            f"  (sigma = {sig})")
        if worst_sp > 1e-6 or worst_off > 1e-6:
            fails.append(f"(C) flat-at--sigma fails on {gname} at Delta={D}")

# ---------------------------------------------------------- (D) the catalogue
log()
log("(D) THE RESONANCE CATALOGUE.")
log("    SU(2) at Delta = 1: sector-1 AND sector-2 levels reappearing exactly")
log("    in sector 3 (both legs of the reported resonances).")
for D in (1.0, 2.0):
    E3 = sector_spectrum(6, 3, D)
    for src in (1, 2):
        Es = sector_spectrum(6, src, D)
        shared = sum(1 for e in Es if np.min(np.abs(E3 - e)) < 1e-9)
        want = len(Es) if D == 1.0 else 0
        log(f"      Delta={D}: sector {src}: {shared} of {len(Es)} shared (want {want})")
        if shared != want:
            fails.append(f"(D) SU(2) shared-level count at Delta={D} sector {src}: {shared} != {want}")

log("    Delta = 0.5 half-filling multiplicity HISTOGRAMS (a triple cannot")
log("    masquerade as a double).")
HIST = {4: {1: 6}, 6: {1: 18, 2: 1}, 8: {1: 54, 2: 8}}
for (nn, want) in HIST.items():
    h = mult_histogram(sector_spectrum(nn, nn // 2, 0.5))
    log(f"      N={nn}: {h} (want {want})")
    if h != want:
        fails.append(f"(D) Delta=0.5 multiplicity histogram at N={nn}: {h} != {want}")

log("    The Delta = 0.5 shared-level census (sector 1 in sector N/2):")
log("    0 / 0 / 8 at N = 4 / 6 / 8 - at N = 8 the resonance is CROSS-sector.")
for (nn, want) in ((4, 0), (6, 0), (8, 8)):
    E1 = sector_spectrum(nn, 1, 0.5)
    Eh = sector_spectrum(nn, nn // 2, 0.5)
    shared = sum(1 for e in E1 if np.min(np.abs(Eh - e)) < 1e-9)
    log(f"      N={nn}: {shared} of {len(E1)} shared (want {want})")
    if shared != want:
        fails.append(f"(D) Delta=0.5 shared census at N={nn}: {shared} != {want}")

log("    The N = 6 doublet is EXACT: 2H built in exact ints from bit logic,")
log("    compared == to the float route; Faddeev-LeVerrier charpoly over Z.")
idx6 = [a for a in range(64) if popcount(a) == 3]
m = len(idx6)
pos = {a: i for i, a in enumerate(idx6)}
M = [[0] * m for _ in range(m)]
for a in idx6:
    diag = 0
    for l in range(5):
        za = 1 - 2 * ((a >> l) & 1)
        zb = 1 - 2 * ((a >> (l + 1)) & 1)
        diag += za * zb                      # 2H diagonal: 2*J*Delta*sum = sum at Delta=1/2, J=1
        if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
            M[pos[a]][pos[a ^ (1 << l) ^ (1 << (l + 1))]] += 4   # 2 * hop 2J
    M[pos[a]][pos[a]] = diag
H6f = build_H(6, 1.0, 0.5)
mismatch = max(abs(2 * H6f[a, b] - M[pos[a]][pos[b]]) for a in idx6 for b in idx6)
log(f"      max|2H_float - 2H_int| = {mismatch}")
if mismatch != 0.0:
    fails.append(f"(D) integer 2H does not match the float route: {mismatch}")
MF = [[Fraction(x) for x in row] for row in M]
c = [Fraction(1)]
Mk = [[Fraction(int(i == j)) for j in range(m)] for i in range(m)]
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(m)] for i in range(m)]
for k in range(1, m + 1):
    Ak = matmul(MF, Mk)
    tr = sum(Ak[i][i] for i in range(m))
    ck = -tr / k
    c.append(ck)
    Mk = [row[:] for row in Ak]
    for i in range(m):
        Mk[i][i] += ck
def poly_eval(coeffs, x):
    v = Fraction(0)
    for a in coeffs:
        v = v * x + a
    return v
def poly_deriv(coeffs):
    deg = len(coeffs) - 1
    return [a * (deg - i) for i, a in enumerate(coeffs[:-1])]
p0 = poly_eval(c, Fraction(5))
d1 = poly_deriv(c)
p1 = poly_eval(d1, Fraction(5))
p2 = poly_eval(poly_deriv(d1), Fraction(5))
log(f"      p(5) = {p0}, p'(5) = {p1}, p''(5) = {p2}")
if not (p0 == 0 and p1 == 0 and p2 != 0):
    fails.append(f"(D) exact multiplicity of 5 in charpoly(2H) is not 2")

log("    (3,3): ker(ad_H) is the commutant, dim = sum(m_i^2) exactly (>= 20")
log("    holds for ANY H; the equality plus the hard-coded 22/20/20 is the pin).")
for (D, pin) in ((0.5, 22), (1.0, 20), (2.0, 20)):
    A, _, _ = split_A_D(6, 3, 3, D, gl)
    w = np.linalg.eigvalsh(1j * A)
    dim0 = int(np.sum(np.abs(w) < 1e-9))
    h = mult_histogram(sector_spectrum(6, 3, D))
    want = sum(mult * mult * cnt for mult, cnt in h.items())
    log(f"      Delta={D}: dim ker = {dim0}, sum(m^2) from sector = {want}, pinned = {pin}")
    if dim0 != want or dim0 != pin:
        fails.append(f"(D) (3,3) commutant dim at Delta={D}: {dim0} vs sum(m^2)={want}, pin={pin}")

# ------------------------------------------ (E) family ties, fourth behaviour
log()
log("(E) F153's off-locus (0, N/2) family, reproduced by the rule:")
FAM = {0.5: {4: 0.0000, 6: 0.6456, 8: 1.2903},
       0.0: {4: 0.2236, 6: 0.5325, 8: 1.0756}}
for D, row in FAM.items():
    for nn, want in row.items():
        r = predicted_res(nn, 0, nn // 2, D, free_profile(nn))
        v = r.max() - r.min()
        log(f"    Delta={D} N={nn} (0,{nn//2}): predicted spread {v:.4f} (want {want})")
        if abs(v - want) > 2e-4:
            fails.append(f"(E) family spread at Delta={D} N={nn}: {v:.4f} != {want}")

log("    THE FOURTH BEHAVIOUR: (1,2) at Delta = 2, locus profile; both indices")
log("    off N/2, no cancellation: a stable density-overlap band.")
pred = predicted_res(6, 1, 2, 2.0, gl)
m4 = measured_res(6, 1, 2, 2.0, gl, 1e4)
m6 = measured_res(6, 1, 2, 2.0, gl, 1e6)
sp_p, sp4, sp6 = pred.max() - pred.min(), m4.max() - m4.min(), m6.max() - m6.min()
log(f"    predicted spread {sp_p:.4f} range [{pred.min():.4f},{pred.max():.4f}]")
log(f"    measured  spread {sp4:.4f} (J=1e4), {sp6:.4f} (J=1e6), max|pred-meas(1e6)| = {np.max(np.abs(pred - m6)):.2e}")
for (name, v) in (("predicted", sp_p), ("J=1e4", sp4), ("J=1e6", sp6)):
    if abs(v - 0.9982) > 1e-4:
        fails.append(f"(E) (1,2) Delta=2 spread {name} drifted: {v:.4f} != 0.9982")
if abs(pred.min() - (-2.8510)) > 1e-4 or abs(pred.max() - (-1.8528)) > 1e-4:
    fails.append(f"(E) (1,2) range drifted: [{pred.min():.4f},{pred.max():.4f}]")

# --------------------------------------------------------------- the verdict
log()
if fails:
    log("FAILURES:")
    for f in fails:
        log("  " + f)
    log("=" * 78)
    log("VERDICT: FAIL")
else:
    log("=" * 78)
    log("VERDICT: all checks pass")

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
