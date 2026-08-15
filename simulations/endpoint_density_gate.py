"""THE ENDPOINT DENSITY GATE: why the resonances saturate exactly on the
size-class-centre interval, endpoints attained.

Answers the first open item of experiments/THE_SPREAD_IS_A_RESONANCE.md. The
objects: a degenerate eigenspace Omega of ad_H on a coherence block (p, q); N_l
the diagonal indicator that a cell's ket and bra disagree at site l (the
site-resolved face of the dissipator, D = -2 Sum_l gamma_l N_l, and of F122's
size operator N_XY = Sum_l N_l); the site reflection R acting on cells as
(a, b) -> (rev a, rev b), a symmetry of the block exactly when H is
reflection-symmetric, preserving every Omega.

Layer 1, the THEOREM (needs only [H, R] = 0): C_l = Pi_Omega (N_l - N_{N-1-l})
Pi_Omega is R-odd, so on any Omega where R acts as a SCALAR +-1 it satisfies
C_l = -C_l = 0, and on a parity-MIXED Omega its two diagonal parity blocks
vanish, the cross block being unconstrained. Section (A) measures the parity
census: at Delta = 0 every colliding Omega is scalar up to N = 7, and at
Delta = 0.5 up to N = 8, on every block in the census, so there the vanishing
is FORCED.

Layer 2, the LAW (free-fermion content beyond the symmetry): at Delta = 0 the
mixed Omegas that first appear at N = 8 STILL satisfy C_l = 0, at the
eigensolver floor, while a generic R-odd diagonal compressed on the same
Omega has cross norm ~0.1 (section B). The bond control there separates what
can be separated: the mixed collisions are themselves a uniform-chain
resonance (any generic bond profile kills them, palindromic included), and a
non-palindromic J breaks even the forced layer. At Delta = 1 the
mixed Omegas are exactly where the law CAN break (section C): every breaking
eigenspace at N = 6 is mixed, every scalar one is forced, and the omega = 0
carriers of the extremes are scalar.

The WEAK law (Delta = 1): comp(N_l) = (s/N) Id on Omega_s = Omega intersect
(pure size-class-s span), uniform hence profile-blind; FALSE at Delta = 0
(gated from below).

Consequence (section D): on the R90 locus (gamma_l + gamma_{N-1-l} = 2*gbar),
wherever C_l = 0, pairing the sites gives
    Pi_Omega D Pi_Omega = -2*gbar * Pi_Omega N_XY Pi_Omega,
so the compressed spectrum is CONTAINED in [-2*gbar*smax, -2*gbar*smin]
(Rayleigh) and every pure-size-class vector pins an endpoint ON its centre.
Off the locus the consequence is gone (section E); on the self-folded block
(1,3) the two overshoot ends are one fact, the undressed one-sided X^N fold
F: (a,b) -> (a, ~b) with F D F = -D - 2*Sigma (sibling full-block statement
at uniform gamma: PROOF_CODIM1_BY_ADDITIVITY section 7 clause (b)), and the
genuinely two-sided case is (2,4), whose lower centre is an arithmetic
endpoint.

Eigenvalue grouping uses gtol = 1e-9, which sits in the measured gap void of
the difference spectrum on this grid ([3e-14, 4.3e-6], documented at
high_q_selection_gate.predicted_res); the class-purity cut 1 - 1e-9 and the
support cut 1e-12 sit five decades above every floor printed here.

Error model, not a chosen number: every zero here is exact in exact
arithmetic; the float route is eigh's invariant-subspace error,
res <= C * eps * ||A|| / gap(Omega), gap the distance to the nearest ad_H
eigenvalue outside Omega. Gates compare spectral-norm residuals (basis-
invariant; a max-entry residual would depend on LAPACK's arbitrary basis
inside a degenerate space) against 64 * that floor, per eigenspace. Breaks
are gated from below at 0.05, four decades above every floor printed here.

Writes simulations/results/endpoint_density_gate.txt; last line is the VERDICT.
"""

import numpy as np

OUT = "simulations/results/endpoint_density_gate.txt"
lines = []
fails = []
eps = np.finfo(float).eps

def log(s=""):
    lines.append(s)
    print(s)

def popcount(x):
    return bin(x).count("1")

def revbits(a, n):
    r = 0
    for l in range(n):
        if (a >> l) & 1:
            r |= 1 << (n - 1 - l)
    return r

# ----------------------------------------------------------------- primitives
def build_H(n, J, D):
    """Pauli book H = Sum_b J_b (XX + YY + D*ZZ), J scalar or per-bond array."""
    Jb = np.full(n - 1, float(J)) if np.isscalar(J) else np.asarray(J, float)
    d = 1 << n
    H = np.zeros((d, d))
    for a in range(d):
        for l in range(n - 1):
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a, a] += Jb[l] * D * za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2 * Jb[l]
    return H

def block_A(n, p, q, D, J=1.0):
    """-i ad_H on the (p,q) block, cells ket-major (high_q_selection_gate's
    split_A_D at gamma = 0, generalized to a bond profile)."""
    H = build_H(n, J, D)
    kets = [a for a in range(1 << n) if popcount(a) == p]
    bras = [b for b in range(1 << n) if popcount(b) == q]
    cells = [(a, b) for a in kets for b in bras]
    idx = {c: i for i, c in enumerate(cells)}
    m = len(cells)
    A = np.zeros((m, m), dtype=complex)
    for (a, b), i in idx.items():
        for c in kets:
            if H[a, c] != 0:
                A[i, idx[(c, b)]] += -1j * H[a, c]
        for d_ in bras:
            if H[b, d_] != 0:
                A[i, idx[(a, d_)]] += 1j * H[b, d_]
    return A, cells, idx

def site_indicators(n, cells):
    return [np.array([1.0 if (a ^ b) & (1 << l) else 0.0 for (a, b) in cells])
            for l in range(n)]

def reflection_matrix(n, cells, idx):
    m = len(cells)
    R = np.zeros((m, m))
    for (a, b), i in idx.items():
        R[idx[(revbits(a, n), revbits(b, n))], i] = 1.0
    return R

def eigenspaces(A, gtol=1e-9):
    """Orthonormal bases of the ad_H eigenspaces, with the outside gap."""
    w, V = np.linalg.eigh(1j * A)
    nrm = np.abs(w).max()
    out = []
    start = 0
    for i in range(1, len(w) + 1):
        if i == len(w) or abs(w[i] - w[start]) > gtol:
            gap = np.inf
            if start > 0:
                gap = min(gap, w[start] - w[start - 1])
            if i < len(w):
                gap = min(gap, w[i] - w[i - 1])
            out.append((w[start], V[:, start:i], gap))
            start = i
    return out, nrm

def locus_profile(n):
    g = np.array([(l + 1) / 16 for l in range(n // 2)])
    return np.concatenate([g, (1 - g)[::-1]])

def free_profile(n):
    return np.array([(l + 1) ** 2 / 16 for l in range(n)])

def parity_of(P, R):
    """(n_plus, n_minus) of R restricted to the eigenspace (an involution when
    [H, R] = 0); mixed <=> both nonzero."""
    dim = P.shape[1]
    RW = P.conj().T @ R @ P
    inv_res = float(np.max(np.abs(RW @ RW - np.eye(dim))))
    ev = np.linalg.eigvalsh((RW + RW.conj().T) / 2)
    nplus = int(np.sum(ev > 0))
    return nplus, dim - nplus, inv_res

def strong_rows(n, p, q, D, J=1.0):
    """Per-eigenspace rows (omega, dim, nplus, nminus, res2, floor); res2 the
    SPECTRAL norm of the worst C_l."""
    A, cells, idx = block_A(n, p, q, D, J)
    Nl = site_indicators(n, cells)
    R = reflection_matrix(n, cells, idx)
    spaces, nrm = eigenspaces(A)
    rows = []
    for wv, P, gap in spaces:
        if P.shape[1] < 2:
            continue
        npl, nmi, inv_res = parity_of(P, R)
        r = max(float(np.linalg.norm(
            P.conj().T @ ((Nl[l] - Nl[n - 1 - l])[:, None] * P), 2))
            for l in range(n // 2))
        rows.append((wv, P.shape[1], npl, nmi, r, eps * nrm / gap, inv_res))
    return rows

def weak_rows(n, p, q, D):
    """Per-(eigenspace, size class) rows (omega, s, dimWs, res2, floor)."""
    A, cells, idx = block_A(n, p, q, D)
    S = np.array([popcount(a ^ b) for (a, b) in cells])
    Nl = site_indicators(n, cells)
    spaces, nrm = eigenspaces(A)
    rows = []
    for wv, P, gap in spaces:
        if P.shape[1] < 2:
            continue
        for s in sorted(set(S)):
            mask = (S == s).astype(float)
            Ps = P * mask[:, None]
            evg, Ug = np.linalg.eigh(Ps.conj().T @ Ps)
            keep = evg > 1 - 1e-9
            dWs = int(keep.sum())
            if dWs == 0:
                continue
            B = P @ Ug[:, keep]
            r = max(float(np.linalg.norm(
                B.conj().T @ (Nl[l][:, None] * B) - (s / n) * np.eye(dWs), 2))
                for l in range(n))
            rows.append((wv, s, dWs, r, eps * nrm / gap))
    return rows

log("=" * 78)
log("THE ENDPOINT DENSITY GATE: the theorem, the law, and the consequence")
log("=" * 78)

# ------------------------------------ (A) parity census + the forced theorem
log()
log("(A) PARITY CENSUS + THEOREM LAYER. Per case: multi-dim eigenspaces, how")
log("    many are parity-MIXED, and the strong-law spectral-norm residual on")
log("    the SCALAR (single-parity) ones, where C_l = 0 is forced by [H,R]=0")
log("    (gated at 64 * eps*||A||/gap as a consistency check). At Delta = 0")
log("    every colliding eigenspace is scalar up to N = 7; at Delta = 0.5 up")
log("    to N = 8: there the whole census is the theorem, not the law.")
census_cases = (
    [(4, p, q, D) for D in (0.0, 0.5) for (p, q) in [(1, 2), (1, 3)]]
    + [(5, p, q, 0.0) for (p, q) in [(1, 3), (2, 3)]]
    + [(6, p, q, D) for D in (0.0, 0.5)
       for (p, q) in [(1, 2), (1, 3), (2, 3), (2, 4), (1, 4), (3, 3), (1, 5), (2, 5)]]
    + [(7, p, q, 0.0) for (p, q) in [(1, 3), (1, 4), (2, 4)]]
    + [(8, p, q, 0.0) for (p, q) in [(1, 3), (1, 4), (2, 4)]]
    + [(8, p, q, 0.5) for (p, q) in [(1, 3), (1, 4)]]
)
mixed_by_case = {}
maxdim_census = 0
for (n, p, q, D) in census_cases:
    rows = strong_rows(n, p, q, D)
    if not rows:
        log(f"    N={n} ({p},{q}) Delta={D}: no multi-dim eigenspaces (vacuous)")
        mixed_by_case[(n, p, q, D)] = (0, 0)
        continue
    maxdim_census = max(maxdim_census, max(d for (_, d, *_ ) in rows))
    mixed = [row for row in rows if row[2] > 0 and row[3] > 0]
    scalar = [row for row in rows if not (row[2] > 0 and row[3] > 0)]
    mixed_by_case[(n, p, q, D)] = (len(rows), len(mixed))
    ok = all(r <= 64 * f for (_, _, _, _, r, f, _) in scalar)
    ok_inv = all(iv <= 64 * f for (_, _, _, _, _, f, iv) in rows)
    res_s = max((r for (_, _, _, _, r, _, _) in scalar), default=0.0)
    log(f"    N={n} ({p},{q}) Delta={D}: {len(rows):3d} eigenspaces "
        f"({len(mixed):2d} mixed), scalar-side max res {res_s:.2e}  "
        f"{'ok' if ok and ok_inv else 'FAIL'}")
    if not (ok and ok_inv):
        fails.append(f"theorem layer N={n} ({p},{q}) Delta={D}")
log(f"    max eigenspace dim in the census: {maxdim_census}")
n_mixed_small = sum(m for ((n, p, q, D), (t, m)) in mixed_by_case.items()
                    if n <= 7 or D == 0.5)
n_mixed_n8d0 = sum(m for ((n, p, q, D), (t, m)) in mixed_by_case.items()
                   if n == 8 and D == 0.0)
log(f"    mixed eigenspaces at N <= 7 or Delta = 0.5: {n_mixed_small} "
    f"(the theorem covers all of those census rows)")
log(f"    mixed eigenspaces at N = 8, Delta = 0: {n_mixed_n8d0}")
if n_mixed_small != 0:
    fails.append("unexpected mixed eigenspace below the N=8 Delta=0 frontier")
if n_mixed_n8d0 == 0:
    fails.append("no mixed eigenspaces at N=8 Delta=0 (the law has no content)")

# --------------------------------------- (B) the LAW at Delta = 0 on N = 8
log()
log("(B) THE LAW: at Delta = 0 the parity-MIXED eigenspaces of N = 8 still")
log("    satisfy C_l = 0 at the floor, where a generic R-odd diagonal,")
log("    compressed on the same eigenspace, has cross norm ~0.1. Bond")
log("    control on (1,3): the mixed collisions are a UNIFORM-chain")
log("    resonance (generic bonds kill them, palindromic or not), and a")
log("    non-palindromic J breaks even the forced layer, R being no")
log("    symmetry there; the two generic profiles share collision dims.")
rng = np.random.default_rng(20260815)
for (p, q) in [(1, 3), (1, 4), (2, 4)]:
    rows = strong_rows(8, p, q, 0.0)
    mixed = [row for row in rows if row[2] > 0 and row[3] > 0]
    ok = all(r <= 64 * f for (_, _, _, _, r, f, _) in mixed)
    log(f"    N=8 ({p},{q}): {len(mixed)} mixed eigenspaces, law max res "
        f"{max(r for (_, _, _, _, r, _, _) in mixed):.2e}  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"the law on mixed spaces N=8 ({p},{q})")
# generic-R-odd contrast on the largest mixed space of (1,3)
A, cells, idx = block_A(8, 1, 3, 0.0)
R = reflection_matrix(8, cells, idx)
spaces, nrm = eigenspaces(A)
contrast = 0.0
for wv, P, gap in spaces:
    if P.shape[1] < 2:
        continue
    npl, nmi, _ = parity_of(P, R)
    if npl == 0 or nmi == 0:
        continue
    d = rng.uniform(-1, 1, len(cells))
    d = (d - d[[idx[(revbits(a, 8), revbits(b, 8))] for (a, b) in cells]]) / 2
    contrast = max(contrast, float(np.linalg.norm(P.conj().T @ (d[:, None] * P), 2)))
log(f"    generic R-odd diagonal on the same mixed spaces: max cross norm "
    f"{contrast:.3f} (gate > 0.01; the law's zeros are not support artifacts)")
if not contrast > 0.01:
    fails.append("generic R-odd contrast missing")
J_pal = np.array([1.3, 0.7, 2.1, 1.4, 2.1, 0.7, 1.3])
J_non = np.array([1.3, 0.7, 2.1, 1.4, 1.9, 0.9, 1.6])
rows_p = strong_rows(8, 1, 3, 0.0, J_pal)
rows_n = strong_rows(8, 1, 3, 0.0, J_non)
mixed_p = [row for row in rows_p if row[2] > 0 and row[3] > 0]
res_n = max(r for (_, _, _, _, r, _, _) in rows_n)
dims_p = sorted(d for (_, d, *_ ) in rows_p)
dims_n = sorted(d for (_, d, *_ ) in rows_n)
ok_p = all(r <= 64 * f for (_, _, _, _, r, f, _) in rows_p)
log(f"    palindromic non-uniform J: {len(mixed_p)} mixed of {len(rows_p)} "
    f"eigenspaces (the mixed collisions are a UNIFORM-chain resonance and die "
    f"under any generic bond profile, palindromic included, so only the "
    f"theorem layer remains and it holds: max res "
    f"{max(r for (_, _, _, _, r, _, _) in rows_p):.2e})  "
    f"{'ok' if ok_p else 'FAIL'}")
log(f"    non-palindromic J: max res {res_n:.3f} (gate > 0.05; R is no longer "
    f"a symmetry, so even the forced layer breaks); "
    f"same collision dims as palindromic: {dims_p == dims_n}")
if mixed_p:
    fails.append("mixed spaces unexpectedly survive a generic palindromic J")
if not ok_p:
    fails.append("palindromic-J theorem layer broken")
if not res_n > 0.05:
    fails.append("non-palindromic-J break missing")
if dims_p != dims_n:
    fails.append("collision structure changed between the two generic bond profiles")
# one profile is not "generic": random sweep over palindromic bond profiles
n_sweep, sweep_mixed, sweep_res, sweep_min_spaces = 10, 0, 0.0, 10 ** 9
for k in range(n_sweep):
    half = rng.uniform(0.5, 2.0, 4)
    Jr = np.concatenate([half[:3], [half[3]], half[:3][::-1]])
    rows_r = strong_rows(8, 1, 3, 0.0, Jr)
    sweep_min_spaces = min(sweep_min_spaces, len(rows_r))
    sweep_mixed += sum(1 for row in rows_r if row[2] > 0 and row[3] > 0)
    sweep_res = max(sweep_res, max((r / f for (_, _, _, _, r, f, _) in rows_r),
                                   default=0.0))
log(f"    random palindromic sweep ({n_sweep} profiles on (1,3)): "
    f"{sweep_mixed} mixed spaces total, >= {sweep_min_spaces} eigenspaces per "
    f"profile, worst res/floor {sweep_res:.1f} "
    f"(gate: 0 mixed, res <= 64*floor everywhere, eigenspaces present)")
if sweep_mixed != 0 or sweep_res > 64 or sweep_min_spaces < 10:
    fails.append("random palindromic sweep")

# ------------------------------- (C) Delta = 1: parity characterization + weak
log()
log("(C) Delta = 1. Every eigenspace that breaks C_l = 0 at N = 6 is parity-")
log("    MIXED (the scalar ones are forced); mixedness is NOT sufficient on")
log("    (1,3), where two mixed spaces hold at the floor. The omega = 0")
log("    carriers of the extremes are scalar. At N = 8 (1,3) and (2,4) every")
log("    eigenspace is scalar, so the Delta = 1 strong law there is entirely")
log("    the theorem. WEAK law comp(N_l) = (s/N) Id on Omega_s: holds at")
log("    Delta = 1 everywhere tested, FAILS at Delta = 0 (gated from below).")
for (n, p, q) in [(6, 1, 3), (6, 2, 4)]:
    rows = strong_rows(n, p, q, 1.0)
    bad_scalar = [row for row in rows
                  if not (row[2] > 0 and row[3] > 0) and row[4] > 64 * row[5]]
    breaks = [row for row in rows if row[4] > 0.05]
    all_breaks_mixed = all(row[2] > 0 and row[3] > 0 for row in breaks)
    zero_rows = [row for row in rows if abs(row[0]) < 1e-9]
    zero_scalar = all(not (row[2] > 0 and row[3] > 0) for row in zero_rows)
    holds_mixed = [row for row in rows
                   if row[2] > 0 and row[3] > 0 and row[4] <= 64 * row[5]]
    log(f"    N={n} ({p},{q}): breaks: "
        + ", ".join(f"omega={row[0]:+.0f} dim={row[1]} ({row[2]}+,{row[3]}-) "
                    f"res={row[4]:.3f}" for row in breaks)
        + "; mixed-that-hold: "
        + (", ".join(f"omega={row[0]:+.0f} dim={row[1]} ({row[2]}+,{row[3]}-)"
                     for row in holds_mixed) or "none")
        + f"; omega=0 scalar: {zero_scalar}")
    if bad_scalar or not all_breaks_mixed or not zero_scalar or not breaks:
        fails.append(f"Delta=1 parity characterization N={n} ({p},{q})")
for (n, p, q) in [(8, 1, 3), (8, 2, 4)]:
    rows = strong_rows(n, p, q, 1.0)
    n_mixed = sum(1 for row in rows if row[2] > 0 and row[3] > 0)
    ok = n_mixed == 0 and all(r <= 64 * f for (_, _, _, _, r, f, _) in rows)
    log(f"    N={n} ({p},{q}): {len(rows)} eigenspaces, {n_mixed} mixed, "
        f"max res {max(r for (_, _, _, _, r, _, _) in rows):.2e}  "
        f"{'ok (all scalar, all forced)' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"Delta=1 N=8 scalar census ({p},{q})")
for (n, p, q) in [(6, 1, 2), (6, 1, 3), (6, 2, 3), (6, 2, 4), (6, 1, 4),
                  (6, 3, 3), (6, 1, 5), (6, 2, 5), (8, 1, 3), (8, 1, 4)]:
    rows = weak_rows(n, p, q, 1.0)
    ok = all(r <= 64 * f for (_, _, _, r, f) in rows)
    res = max((r for (_, _, _, r, _) in rows), default=0.0)
    log(f"    weak law N={n} ({p},{q}) Delta=1: {len(rows):2d} Omega_s, "
        f"max res {res:.2e}  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"weak law N={n} ({p},{q})")
weak_d0 = max(max((r for (_, _, _, r, _) in weak_rows(6, p, q, 0.0)), default=0.0)
              for (p, q) in [(1, 3), (2, 4)])
log(f"    weak law at Delta = 0 (N=6, (1,3) and (2,4)): max res {weak_d0:.3f} "
    f"(gate > 0.1: the two resonances carry different density structure)")
if not weak_d0 > 0.1:
    fails.append("Delta=0 weak-law failure missing")

# ------------------------------- (D) the consequence on the locus, and inside
log()
log("(D) On the R90 locus, wherever C_l = 0 (forced or lawful): the")
log("    compression IS -2*gbar * comp(N_XY); containment holds for EVERY")
log("    multi-dim eigenspace (breaking ones included, with slack); both")
log("    endpoints are attained by pure-class vectors where they exist; and")
log("    the dim-12 interior of (1,3) at Delta = 0 is a rational ladder.")
def consequence(n, p, q, D, gamma, skip_broken=False):
    gbar = gamma.mean()
    A, cells, idx = block_A(n, p, q, D)
    S = np.array([popcount(a ^ b) for (a, b) in cells])
    Nl = site_indicators(n, cells)
    rate = np.array([sum(gamma[l] for l in range(n) if ((a ^ b) >> l) & 1)
                     for (a, b) in cells])
    spaces, nrm = eigenspaces(A)
    worst_id = worst_cont = worst_att = 0.0
    attained = 0
    for wv, P, gap in spaces:
        dim = P.shape[1]
        if dim < 2:
            continue
        floor = eps * nrm / gap
        tol_c = 64 * floor * max(1.0, 2 * gbar * S.max())
        M = P.conj().T @ ((-2 * rate)[:, None] * P)
        ev = np.linalg.eigvalsh(M)
        # containment for EVERY eigenspace, against its own class support
        supp = (np.abs(P) ** 2).sum(axis=1) > 1e-12
        allS = sorted(set(S[supp]))
        worst_cont = max(worst_cont,
                         max(0.0, (-2 * gbar * max(allS)) - ev.min(),
                             ev.max() - (-2 * gbar * min(allS))) / tol_c)
        broken = max(float(np.linalg.norm(
            P.conj().T @ ((Nl[l] - Nl[n - 1 - l])[:, None] * P), 2))
            for l in range(n // 2)) > 0.05
        if skip_broken and broken:
            continue
        CS = P.conj().T @ (S.astype(float)[:, None] * P)
        worst_id = max(worst_id,
                       float(np.linalg.norm(M + 2 * gbar * CS, 2)) / tol_c)
        # attainment at BOTH ends where the extreme class has a pure vector
        sW = []
        for s in sorted(set(S)):
            mask = (S == s).astype(float)
            evg = np.linalg.eigvalsh((P * mask[:, None]).conj().T @ (P * mask[:, None]))
            if evg.max() > 1 - 1e-9:
                sW.append(s)
        if sW:
            if max(sW) == max(allS):
                worst_att = max(worst_att, abs(ev.min() + 2 * gbar * max(sW)) / tol_c)
            if min(sW) == min(allS):
                worst_att = max(worst_att, abs(ev.max() + 2 * gbar * min(sW)) / tol_c)
            attained += 1
    return worst_id, worst_cont, worst_att, attained

for (n, p, q, D, skip, expect_pure) in [
        (6, 1, 3, 0.0, False, None), (6, 2, 4, 0.0, False, None),
        (6, 1, 3, 0.5, False, 0), (6, 1, 3, 1.0, True, None),
        (6, 2, 4, 1.0, True, None), (8, 1, 4, 0.0, False, None)]:
    gamma = locus_profile(n)
    wid, wcont, watt, att = consequence(n, p, q, D, gamma, skip_broken=skip)
    ok = wid <= 1.0 and wcont <= 1.0 and watt <= 1.0
    if expect_pure is not None and att != expect_pure:
        ok = False
    log(f"    N={n} ({p},{q}) Delta={D}{' [identity on law-abiding spaces]' if skip else ''}: "
        f"identity ratio {wid:.3f}, containment ratio {wcont:.3f} (ALL spaces), "
        f"attainment ratio {watt:.3f}, {att} eigenspaces with pure vectors"
        f"{' (gated == 0: no attainment possible, the small resonance)' if expect_pure == 0 else ''}  "
        f"{'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"consequence N={n} ({p},{q}) Delta={D}")

def block_extremes(n, p, q, D, gamma):
    gbar = gamma.mean()
    A, cells, idx = block_A(n, p, q, D)
    rate = np.array([sum(gamma[l] for l in range(n) if ((a ^ b) >> l) & 1)
                     for (a, b) in cells])
    spaces, nrm = eigenspaces(A)
    lo, hi = np.inf, -np.inf
    gap_min = np.inf
    for wv, P, gap in spaces:
        if P.shape[1] < 2:
            continue
        ev = np.linalg.eigvalsh(P.conj().T @ ((-2 * rate)[:, None] * P))
        lo, hi = min(lo, ev.min()), max(hi, ev.max())
        gap_min = min(gap_min, gap)
    return lo, hi, 64 * eps * nrm / gap_min

for (n, p, q, D) in [(6, 1, 3, 0.0), (6, 1, 3, 1.0), (6, 2, 4, 0.0), (8, 1, 4, 0.0)]:
    gamma = locus_profile(n)
    gbar = gamma.mean()
    smin, smax = abs(p - q), min(p + q, 2 * n - p - q)
    lo, hi, floor_blk = block_extremes(n, p, q, D, gamma)
    tol = floor_blk * max(1.0, 2 * gbar * smax)
    dev = max(abs(lo + 2 * gbar * smax), abs(hi + 2 * gbar * smin))
    log(f"    block saturation N={n} ({p},{q}) Delta={D}: range "
        f"[{lo:.6f},{hi:.6f}] = centres to {dev:.2e} vs floor {tol:.2e} "
        f"(spread {hi - lo:.4f})  {'ok' if dev <= tol else 'FAIL'}")
    if not dev <= tol:
        fails.append(f"block saturation N={n} ({p},{q}) Delta={D}")

# the rational interior: 7 * spec(comp N_XY) on the dim-12 space = {14,20,22,28} x3
A, cells, idx = block_A(6, 1, 3, 0.0)
S = np.array([popcount(a ^ b) for (a, b) in cells])
spaces, nrm = eigenspaces(A)
big = [(wv, P, gap) for (wv, P, gap) in spaces if P.shape[1] == 12]
if len(big) != 1:
    fails.append("dim-12 eigenspace not unique on (1,3) Delta=0")
else:
    wv, P, gap = big[0]
    ev7 = 7 * np.linalg.eigvalsh(P.conj().T @ (S.astype(float)[:, None] * P))
    ints = np.round(ev7).astype(int)
    frac = float(np.max(np.abs(ev7 - ints)))
    want = sorted([14] * 3 + [20] * 3 + [22] * 3 + [28] * 3)
    got = sorted(ints.tolist())
    ok = got == want and frac <= 64 * 7 * 4 * eps * nrm / gap and abs(wv) < 1e-9
    log(f"    dim-12 interior (at omega = {wv:+.2e}): 7*spec(comp N_XY) = "
        f"{sorted(set(got))} x3, integrality defect {frac:.2e}  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append("rational interior of the dim-12 space")

# ------------------------------------------------------ (E) off-locus fences
log()
log("(E) Off the locus the consequence is gone. On the SELF-FOLDED block (1,3)")
log("    (bra index N/2) the two overshoot ends are ONE fact: the ket/bra fold")
log("    F: (a,b) -> (a, ~b) is an exact involution, commutes with ad_H (from")
log("    [H, X^N] = 0) and satisfies F D F = -D - 2*Sigma entry-wise (compared")
log("    == 0.0; the dyadic profile makes the float sums exact), so every")
log("    compressed spectrum is symmetric about -Sigma for ANY profile (the")
log("    sibling full-block statement at uniform gamma is the s-symmetric")
log("    self-folded block, PROOF_CODIM1_BY_ADDITIVITY section 7 clause (b);")
log("    F here is the undressed one-sided X^N, any Delta, any profile).")
log("    The genuinely two-sided case is (2,4): its lower")
log("    centre -2*gbar*N is the global minimum of the rate diagonal for any")
log("    profile (an arithmetic endpoint), and only its upper end escapes.")
gf = free_profile(6)
Sig = gf.sum()
# fold identities on (1,3), exact
for D in (0.0, 1.0):
    A, cells, idx = block_A(6, 1, 3, D)
    m = len(cells)
    mask_full = (1 << 6) - 1
    F = np.zeros((m, m))
    for (a, b), i in idx.items():
        F[idx[(a, b ^ mask_full)], i] = 1.0
    rate = np.array([sum(gf[l] for l in range(6) if ((a ^ b) >> l) & 1)
                     for (a, b) in cells])
    Dm = np.diag(-2 * rate)
    r_inv = float(np.max(np.abs(F @ F - np.eye(m))))
    r_com = float(np.max(np.abs(A @ F - F @ A)))
    r_fold = float(np.max(np.abs(F @ Dm @ F + Dm + 2 * Sig * np.eye(m))))
    ok = r_inv == 0.0 and r_com == 0.0 and r_fold == 0.0
    log(f"    fold identities on (1,3) Delta={D}: F^2-I = {r_inv}, [ad_H,F] = "
        f"{r_com}, F D F + D + 2*Sigma = {r_fold}  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"fold identity Delta={D}")

def overshoots(n, p, q, D, gamma):
    gbar_ = gamma.mean()
    A, cells, idx = block_A(n, p, q, D)
    S = np.array([popcount(a ^ b) for (a, b) in cells])
    rate = np.array([sum(gamma[l] for l in range(n) if ((a ^ b) >> l) & 1)
                     for (a, b) in cells])
    spaces, nrm = eigenspaces(A)
    lo, hi = np.inf, -np.inf
    gap_min = np.inf
    for wv, P, gap in spaces:
        if P.shape[1] < 2:
            continue
        ev = np.linalg.eigvalsh(P.conj().T @ ((-2 * rate)[:, None] * P))
        lo, hi = min(lo, ev.min()), max(hi, ev.max())
        gap_min = min(gap_min, gap)
    floor = 64 * eps * nrm / gap_min * max(1.0, 2 * gbar_ * S.max())
    return (-2 * gbar_ * S.max()) - lo, hi - (-2 * gbar_ * S.min()), floor

for D, gate_lo in [(0.0, 0.1), (1.0, 1e-3)]:
    over_lo, over_hi, floor = overshoots(6, 1, 3, D, gf)
    ok = over_lo > gate_lo and over_hi > gate_lo and abs(over_lo - over_hi) <= floor
    log(f"    off-locus Delta={D} (1,3): overshoot lower {over_lo:.10f}, upper "
        f"{over_hi:.10f}, |diff| {abs(over_lo - over_hi):.1e} <= fold floor "
        f"{floor:.1e} (both gated > {gate_lo})  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"off-locus Delta={D} overshoot/fold")
over_lo, over_hi, floor = overshoots(6, 2, 4, 0.0, gf)
ok = over_lo <= floor and over_hi > 0.1
log(f"    off-locus Delta=0.0 (2,4): lower {over_lo:.2e} (arithmetic endpoint, "
    f"gate <= floor {floor:.1e}), upper {over_hi:.4f} (gate > 0.1)  "
    f"{'ok' if ok else 'FAIL'}")
if not ok:
    fails.append("off-locus (2,4) asymmetric overshoot")

log()
log("=" * 78)
if fails:
    log("VERDICT: FAIL: " + "; ".join(fails))
else:
    log("VERDICT: all checks pass")
with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
