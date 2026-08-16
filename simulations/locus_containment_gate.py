"""THE LOCUS CONTAINMENT GATE: on the mirror-balanced locus the breaking
mixed eigenspaces' compressed dissipator spectra stay inside the
size-class-centre interval for EVERY physical profile, by a mirror-transversal
certificate in closed form.

Answers the second half of arc `compressed_density_laws` NextStep (3). The
question (endpoint note, F154): the containment
spec(comp D) in [-2*gbar*smax, -2*gbar*smin] was DERIVED wherever C_l = 0
(theorem + locus pairing, any locus profile) but only MEASURED, with slack,
on the Delta = 1, N = 6 eigenspaces that BREAK C_l = 0, and only at one
interior profile. Three steps close it.

STEP 1, the REDUCTION (an argument, not a computation): comp(D) =
-2 Sum_l gamma_l comp(N_l) is affine in the profile (the eigenspace Omega is
built from ad_H alone, so comp(N_l) is profile-independent), and the
physical locus set {gamma_l >= 0, gamma_l + gamma_{N-1-l} = 2*gbar} is a
product of floor(N/2) segments whose corners are the MIRROR TRANSVERSALS A
(one site of each pair lit at 2*gbar, its partner dark). lambda_max of a
Hermitian affine family is CONVEX in the parameters and lambda_min is
CONCAVE (the variational characterization), so over the box both extremes
are attained at corners: containment for every physical locus profile is
EQUIVALENT to containment at the 2^floor(N/2) transversals. At EVEN N
(this gate's case) comp(D) = -4*gbar*comp(N_A) there; at odd N the middle
site is its own mirror, stays pinned at gbar, and adds a fixed term
-2*gbar*comp(N_mid) to every corner (and N_A + N_A^c = N_XY - N_mid), so
the reduction survives verbatim while the corner FORM carries the middle
term. Section (B) cross-checks the convexity consequence on random
profiles.

STEP 2, the SPLIT: comp(N_A) = (comp(N_XY) + Sum_l sigma_l C_l)/2, sigma
the transversal's signs, comp(N_XY) R-parity-even (gated), each C_l
R-parity-odd, i.e. supported on cross-parity dyad pairs only (gated;
THE_TWO_SPIN_ZEROS' object). comp(N_XY) is DIAGONAL in the multiplet-dyad
eigenbasis with /24-rational entries (gated, values printed); it is SCALAR
on all 24 dim-3 spaces and on the 4 dim-2 spaces of the m = 0 blocks, and
genuinely non-scalar on the remaining 20 dim-2 spaces (gaps 1/6, 1/12,
1/24). The interval itself is exact combinatorics, no eigenvector needed:
a (p, q) cell has popcount(a XOR b) == p + q (mod 2), so the block's class
set is {|p-q|, |p-q|+2, ..., min(p+q, 2N-p-q)}; section (E) gates that
every breaking eigenspace carries mass >= 0.1 on every one of those
classes and EXACTLY 0.0 on every other (the support is the full block
ladder, endpoints live). For the m = 0 blocks (q = N/2) the closed form
derives the centre: (|p-3| + p+3)/2 = 3.

STEP 3, the CERTIFICATE, closed-form on all 48: with delta =
ptp(diag comp N_XY)/2 (zero on scalar spaces), w = max_sigma
||Sum_l sigma_l C_l||, and bound = min(smax - mbar, mbar - smin),

    sqrt(delta^2 + w^2) <= bound,

equivalent to spec(comp N_A) in [smin/2, smax/2] at every transversal
(dim-2: exact 2x2 eigenvalues; dim-3: all scalar, delta = 0, spectrum
symmetric about mbar/2). The measured worst transversal slack EQUALS
(bound - sqrt(delta^2 + w^2))/2 on every space (gated at the floor; that
equality also implies the split identity, which section (C) additionally
gates entry-wise at one transversal), so
tight/non-tight is an exact split, no band needed: EQUALITY holds on the
six dim-3 corner spaces ((1,1), (1,5), (5,5) at omega = +-2, w = 1/3 =
bound), where a transversal eigenvalue reaches an interval edge exactly,
the UPPER edge smax/2 on (1,1) and (5,5), the LOWER edge smin/2 on (1,5)
(edge sides gated). Everywhere else the slack is the closed form's
margin:
0.020141 on the six corner dim-2 spaces (where the cruder Weyl bound
delta + w = 0.275783 would EXCEED the 0.25 bound: the exact 2x2 form is
load-bearing there, gated), >= 0.122515 on the remaining 36. Both parts
are rational: 2160 * w^2 is an integer on every space (the same measured
denominator as THE_TWO_SPIN_ZEROS' sizes, a fit, not derived) and
2304 * delta^2 is an integer (delta in {0, 1/12, 1/24, 1/48}), both
gated. Together: containment on the locus box is derived (step 1 + the
split) plus finitely certified in closed form (step 3) on exactly the
spaces where it was only measured; the endpoint note's slack is the
margin the certificate leaves, and the corner class shows it cannot be
improved.

Error model: the certificate quantities are exact in exact arithmetic; the
float route is eigh backward error, floors eps*||H_sector||/gap as in
mixed_break_rule_gate.py; gates compare at 64 * floor (the rationality
gates at den * 4 * tol, den the claimed denominator). The census cut
carries the predecessor's VOID branch (resC between 64*floor and 0.04
FAILS rather than reclassifies), the dyad grouping void is gated >= 1e-6
per block, sector labels are gated ((E, S, r) nondegenerate, R = +-1,
S(S+1) integer), and section (F) bridges every space against an
independent cells-route construction.

Writes simulations/results/locus_containment_gate.txt; VERDICT last line.
"""

import numpy as np

OUT = "simulations/results/locus_containment_gate.txt"
lines = []
fails = []
eps = np.finfo(float).eps

n = 6
GBAR = 0.5


def log(s=""):
    lines.append(s)
    print(s)


def popcount(x):
    return bin(x).count("1")


def revbits(a, nn):
    r = 0
    for l in range(nn):
        if (a >> l) & 1:
            r |= 1 << (nn - 1 - l)
    return r


def build_H(nn, D=1.0, J=1.0):
    d = 1 << nn
    H = np.zeros((d, d))
    for a in range(d):
        for l in range(nn - 1):
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a, a] += J * D * za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2 * J
    return H


def build_S2(nn):
    d = 1 << nn
    S2 = np.zeros((d, d))
    for a in range(d):
        for i in range(nn):
            for j in range(i + 1, nn):
                za = 1 - 2 * ((a >> i) & 1)
                zb = 1 - 2 * ((a >> j) & 1)
                S2[a, a] += 0.5 * za * zb
                if ((a >> i) & 1) != ((a >> j) & 1):
                    S2[a, a ^ (1 << i) ^ (1 << j)] += 1.0
    S2 += (3 * nn / 4) * np.eye(d)
    return S2


def sector_vectors(nn, p, H, S2):
    cfg = [a for a in range(1 << nn) if popcount(a) == p]
    idx = {a: k for k, a in enumerate(cfg)}
    Hp = H[np.ix_(cfg, cfg)]
    S2p = S2[np.ix_(cfg, cfg)]
    Rp = np.zeros((len(cfg), len(cfg)))
    for a in cfg:
        Rp[idx[revbits(a, nn)], idx[a]] = 1.0
    out = []
    gap_min = np.inf
    s2e, s2v = np.linalg.eigh(S2p)
    k = 0
    while k < len(s2e):
        m2 = k
        while m2 + 1 < len(s2e) and s2e[m2 + 1] - s2e[k] < 0.5:
            m2 += 1
        W = s2v[:, k:m2 + 1]
        Sval = (-1 + np.sqrt(1 + 4 * float(np.mean(s2e[k:m2 + 1])))) / 2
        if abs(Sval - round(2 * Sval) / 2) > 1e-9:
            fails.append(f"sector {p}: S^2 eigenvalue not of the form S(S+1)")
        Sval = round(2 * Sval) / 2
        He, Hv = np.linalg.eigh(W.T @ Hp @ W)
        j = 0
        while j < len(He):
            m3 = j
            while m3 + 1 < len(He) and He[m3 + 1] - He[j] < 1e-8:
                m3 += 1
            if j > 0:
                gap_min = min(gap_min, He[j] - He[j - 1])
            U = W @ Hv[:, j:m3 + 1]
            rb = U.T @ Rp @ U
            re, rv = np.linalg.eigh(0.5 * (rb + rb.T))
            Ur = U @ rv
            for c in range(Ur.shape[1]):
                if abs(abs(re[c]) - 1) > 1e-9:
                    fails.append(f"sector {p}: R eigenvalue not +-1")
                out.append(dict(E=float(np.mean(He[j:m3 + 1])), S=Sval,
                                r=int(round(re[c])), v=Ur[:, c]))
            j = m3 + 1
        k = m2 + 1
    labels = [(round(d["E"], 8), d["S"], d["r"]) for d in out]
    if len(labels) != len(set(labels)):
        fails.append(f"sector {p}: degenerate (E, S, r) label (the dyad "
                     f"basis is not well-defined)")
    floor = eps * float(np.linalg.norm(Hp, 2)) / min(gap_min, 2.0)
    return out, cfg, floor


def occ_all(cfg, vecs, nn):
    B = np.column_stack([d["v"] for d in vecs])
    occ = []
    for l in range(nn):
        w = np.array([(a >> l) & 1 for a in cfg], dtype=float)
        occ.append(B.T @ (w[:, None] * B))
    return occ


def mixed_spaces(nn, vp, vq, gtol=1e-9):
    dyads = [(vp[i]["E"] - vq[j]["E"], i, j)
             for i in range(len(vp)) for j in range(len(vq))]
    dyads.sort(key=lambda t: t[0])
    groups = []
    start = 0
    void = np.inf
    for k in range(1, len(dyads) + 1):
        if k == len(dyads) or dyads[k][0] - dyads[k - 1][0] > gtol:
            om = float(np.mean([t[0] for t in dyads[start:k]]))
            groups.append((om, [(i, j) for (_, i, j) in dyads[start:k]]))
            if k < len(dyads):
                void = min(void, dyads[k][0] - dyads[k - 1][0])
            start = k
    mixed = []
    for om, dy in groups:
        if abs(om) < gtol or len(dy) < 2:
            continue
        pars = [vp[i]["r"] * vq[j]["r"] for (i, j) in dy]
        if pars.count(1) > 0 and pars.count(-1) > 0:
            mixed.append((om, dy, pars))
    return mixed, void


def comp_matrices(nn, dy, op, oq):
    dim = len(dy)
    comp = []
    for l in range(nn):
        M = np.zeros((dim, dim))
        for a, (i, j) in enumerate(dy):
            for b, (i2, j2) in enumerate(dy):
                npl = op[l][i, i2]
                nql = oq[l][j, j2]
                dp = 1.0 if i == i2 else 0.0
                dq = 1.0 if j == j2 else 0.0
                M[a, b] = npl * (dq - nql) + (dp - npl) * nql
        comp.append(M)
    return comp


log("=" * 78)
log("THE LOCUS CONTAINMENT GATE: corner reduction + parity split + the")
log("closed-form certificate on the breaking mixed eigenspaces (Delta=1, N=6)")
log("=" * 78)

# ------------------------------------------------- (A) the breaking census
log()
log("(A) SETUP. Sector eigenbases with (E, r) labels; the 48 breaking mixed")
log("    eigenspaces over the fifteen blocks p <= q (THE_TWO_SPIN_ZEROS'")
log("    census; gate: exactly 48, and nothing in the VOID band between")
log("    64*floor and 0.04).")
H = build_H(n)
S2 = build_S2(n)
sector = {}
for p in range(1, 6):
    vecs, cfg, floor = sector_vectors(n, p, H, S2)
    sector[p] = (vecs, cfg, occ_all(cfg, vecs, n), floor)
spaces = []
for p in range(1, 6):
    for q in range(p, 6):
        vp, cfg_p, op, fp = sector[p]
        vq, cfg_q, oq, fq = sector[q]
        mixed, void = mixed_spaces(n, vp, vq)
        if void < 1e-6:
            fails.append(f"({p},{q}): dyad grouping void {void:.1e} < 1e-6")
        for om, dy, pars in mixed:
            comp = comp_matrices(n, dy, op, oq)
            resC = max(np.linalg.norm(comp[l] - comp[n - 1 - l], 2)
                       for l in range(n // 2))
            if resC <= 64 * (fp + fq):
                continue
            if resC < 0.04:
                fails.append(f"resC in the VOID band ({p},{q}) "
                             f"om={om:+.0f}: {resC:.2e}")
                continue
            spaces.append(dict(p=p, q=q, om=om, dy=dy, pars=pars, comp=comp,
                               floor=fp + fq))
log(f"    breaking spaces: {len(spaces)}  "
    f"{'ok' if len(spaces) == 48 else 'FAIL'}")
if len(spaces) != 48:
    fails.append(f"breaking census {len(spaces)} != 48")

# ------------------------------------- (E) exact class support, no cutoff
log()
log("(E) THE INTERVAL IS EXACT COMBINATORICS. popcount(a^b) == p+q (mod 2),")
log("    so the block's class set is {|p-q|, |p-q|+2, ..., min(p+q, 2N-p-q)}.")
log("    Gate per breaking space: mass EXACTLY 0.0 on every class outside")
log("    that set, and >= 0.1 (an enumerated finite minimum) on every class")
log("    inside it: the support is the full block ladder, endpoints live,")
log("    smin = |p-q|, smax = min(p+q, 2N-p-q).")
live_min = np.inf
for sp in spaces:
    p, q, dy = sp["p"], sp["q"], sp["dy"]
    vp, cfg_p, _, _ = sector[p]
    vq, cfg_q, _, _ = sector[q]
    s_cell = np.array([[popcount(a ^ b) for b in cfg_q] for a in cfg_p])
    masses = {s: 0.0 for s in range(n + 1)}
    for (i, j) in dy:
        amp2 = np.outer(vp[i]["v"], vq[j]["v"]) ** 2
        for s in range(n + 1):
            masses[s] += float(amp2[s_cell == s].sum())
    smin, smax = abs(p - q), min(p + q, 2 * n - p - q)
    ladder = set(range(smin, smax + 1, 2))
    for s in range(n + 1):
        if s in ladder:
            live_min = min(live_min, masses[s])
            if masses[s] < 0.1:
                fails.append(f"ladder class {s} underweight ({p},{q}) "
                             f"om={sp['om']:+.0f}: {masses[s]:.3f}")
        elif masses[s] != 0.0:
            fails.append(f"off-ladder class {s} carries mass ({p},{q}) "
                         f"om={sp['om']:+.0f}: {masses[s]:.1e}")
    sp["smin"], sp["smax"] = smin, smax
log(f"    smallest ladder-class mass across all 48: {live_min:.3f}; every")
log("    off-ladder mass exactly 0.0  "
    + ("ok" if not any("ladder" in f for f in fails) else "FAIL"))

# ------------------------------------ (B) corner reduction, numerical check
log()
log("(B) CORNER REDUCTION. comp(D) is affine in gamma (Omega is built from")
log("    ad_H alone); lambda_max is convex and lambda_min concave, so the")
log("    box extremes sit at the 8 mirror transversals. Consequence check:")
log("    for 20 random locus profiles per space, the spectrum stays inside")
log("    the hull of the 8 transversal spectra (worst signed margin printed;")
log("    gate: no overshoot beyond 64 * floor).")
rng = np.random.default_rng(20260816)
worst_margin = np.inf   # positive = strictly inside
worst_over = 0.0
for sp in spaces:
    comp = sp["comp"]
    verts = []
    for bits in range(8):
        A = [l if (bits >> l) & 1 else n - 1 - l for l in range(3)]
        NA = sum(comp[l] for l in A)
        ev = np.linalg.eigvalsh(-4 * GBAR * NA)
        verts.append((ev.min(), ev.max()))
    vlo = min(v[0] for v in verts)
    vhi = max(v[1] for v in verts)
    for _ in range(20):
        half = rng.uniform(0.0, 2 * GBAR, n // 2)
        gamma = np.concatenate([half, (2 * GBAR - half)[::-1]])
        ev = np.linalg.eigvalsh(sum(-2 * gamma[l] * comp[l]
                                    for l in range(n)))
        worst_margin = min(worst_margin, ev.min() - vlo, vhi - ev.max())
        worst_over = max(worst_over, vlo - ev.min(), ev.max() - vhi)
floor_all = max(sp["floor"] for sp in spaces)
ok = worst_over <= 64 * floor_all
log(f"    worst overshoot {worst_over:.2e} (gate <= {64 * floor_all:.1e}); "
    f"tightest interior margin {worst_margin:.2e}  {'ok' if ok else 'FAIL'}")
if not ok:
    fails.append("random profile escapes the transversal hull")

# ---------------- (C) parity split: M0 even+diagonal+rational, C_l odd
log()
log("(C) PARITY SPLIT. Gates: comp(N_XY) commutes with the compressed")
log("    reflection R~ = diag(dyad parities) (even); every C_l anticommutes")
log("    with it (odd, cross-parity support only); comp(N_XY) is diagonal")
log("    in the dyad eigenbasis, 24 * entries integer (values printed);")
log("    scalar exactly on the 24 dim-3 spaces and the 4 dim-2 spaces of")
log("    the m = 0 blocks (ptp of the diagonal <= 64*floor there, and")
log("    > 1/25 on the other 20).")
m0_values = set()
gap_ints = set()
n_scalar = 0
n_dim3 = 0
for sp in spaces:
    comp = sp["comp"]
    dim = len(sp["dy"])
    if dim == 3:
        n_dim3 += 1
    M0 = sum(comp[l] for l in range(n))
    Rt = np.diag([float(x) for x in sp["pars"]])
    r_even = float(np.max(np.abs(Rt @ M0 @ Rt - M0)))
    r_odd = max(float(np.max(np.abs(Rt @ (comp[l] - comp[n - 1 - l]) @ Rt
                                    + (comp[l] - comp[n - 1 - l]))))
                for l in range(n // 2))
    if r_even > 64 * sp["floor"] or r_odd > 64 * sp["floor"]:
        fails.append(f"parity split broken ({sp['p']},{sp['q']}) "
                     f"om={sp['om']:+.0f}")
    off = max(abs(M0[a, b]) for a in range(dim) for b in range(dim)
              if a != b)
    if off > 64 * sp["floor"]:
        fails.append(f"M0 not diagonal ({sp['p']},{sp['q']}) "
                     f"om={sp['om']:+.0f}: {off:.1e}")
    ints = []
    for a in range(dim):
        v24 = 24 * M0[a, a]
        if abs(v24 - round(v24)) > 64 * sp["floor"] * 24:
            fails.append(f"M0 diagonal not /24-rational "
                         f"({sp['p']},{sp['q']}) om={sp['om']:+.0f}")
        ints.append(round(v24))
        m0_values.add(round(v24) / 24)
    # scalar vs non-scalar decided on the INTEGERS, no float fence
    int_gap = max(ints) - min(ints)
    want_scalar = (dim == 3) or (sp["p"] == 3 or sp["q"] == 3)
    if want_scalar != (int_gap == 0):
        fails.append(f"scalarity mismatch ({sp['p']},{sp['q']}) "
                     f"om={sp['om']:+.0f}: 24-integer gap {int_gap}")
    if int_gap == 0:
        n_scalar += 1
    else:
        gap_ints.add(int_gap)
    # the m = 0 blocks measure scalar AT THE CENTRE, mbar = 3 (gated)
    if (sp["p"] == 3 or sp["q"] == 3) and             abs(float(np.mean(np.diag(M0))) - 3.0) > 64 * sp["floor"]:
        fails.append(f"m=0 block not scalar at 3 ({sp['p']},{sp['q']}) "
                     f"om={sp['om']:+.0f}")
    # the split identity itself, entry-wise, at one transversal
    NA0 = sum(comp[l] for l in [0, 1, 2])
    W0 = sum(comp[l] - comp[n - 1 - l] for l in range(n // 2))
    r_split = float(np.max(np.abs(NA0 - (M0 + W0) / 2)))
    if r_split > 64 * sp["floor"]:
        fails.append(f"split identity broken ({sp['p']},{sp['q']}) "
                     f"om={sp['om']:+.0f}: {r_split:.1e}")
    sp["M0"] = M0
want_vals = sorted([5 / 3, 11 / 6, 7 / 3, 29 / 12, 3.0, 8 / 3, 65 / 24,
                    10 / 3, 79 / 24, 11 / 3, 43 / 12, 13 / 3, 25 / 6])
log(f"    scalar spaces: {n_scalar}/48 (gate == 28; 24 dim-3 + 4 m=0 "
    f"dim-2); dims: {n_dim3} dim-3 (gate == 24); /24-rational diagonal "
    f"values: {sorted(m0_values)} (gate: this 13-value set); non-scalar "
    f"24-integer gaps: {sorted(gap_ints)} (gate: {{1, 2, 4}})")
if n_scalar != 28:
    fails.append(f"scalar count {n_scalar} != 28")
if n_dim3 != 24:
    fails.append(f"dim-3 count {n_dim3} != 24")
if sorted(m0_values) != want_vals:
    fails.append("M0 diagonal value set changed")
if gap_ints != {1, 2, 4}:
    fails.append(f"non-scalar gap set {sorted(gap_ints)} != {{1, 2, 4}}")

# ------------------------------------------------- (D) the certificate
log()
log("(D) THE CERTIFICATE, closed-form on all 48. Per space: delta =")
log("    ptp(diag M0)/2, w = max_sigma ||Sum sigma_l C_l||, bound =")
log("    min(smax - mbar, mbar - smin). Gates: sqrt(delta^2 + w^2) <=")
log("    bound (containment); the measured worst transversal slack equals")
log("    (bound - sqrt(delta^2 + w^2))/2 to the floor (the closed form IS")
log("    the spectrum); 2160 * w^2 and 2304 * delta^2 integers;")
log("    equality exactly on the six corner dim-3 spaces, upper edge on")
log("    (1,1)/(5,5), lower on (1,5); on the corner dim-2 spaces the")
log("    cruder Weyl form delta + w exceeds the bound (the 2x2 form is")
log("    load-bearing, gated).")
tight = []
all_slacks = []
for sp in spaces:
    comp = sp["comp"]
    dim = len(sp["dy"])
    smin, smax = sp["smin"], sp["smax"]
    worst_slack = np.inf
    edge = None
    for bits in range(8):
        A = [l if (bits >> l) & 1 else n - 1 - l for l in range(3)]
        NA = sum(comp[l] for l in A)
        ev = np.linalg.eigvalsh(NA)
        lo_sl = ev.min() - smin / 2
        hi_sl = smax / 2 - ev.max()
        if min(lo_sl, hi_sl) < worst_slack:
            worst_slack = min(lo_sl, hi_sl)
            edge = "lower" if lo_sl < hi_sl else "upper"
    Cs = [comp[l] - comp[n - 1 - l] for l in range(n // 2)]
    w = max(np.linalg.norm(sum(s * C for s, C in zip(sg, Cs)), 2)
            for sg in [[1 if (b >> k) & 1 else -1 for k in range(3)]
                       for b in range(8)])
    mbar = float(np.mean(np.diag(sp["M0"])))
    delta = float(np.ptp(np.diag(sp["M0"]))) / 2
    bound = min(smax - mbar, mbar - smin)
    d2 = np.sqrt(delta ** 2 + w ** 2)
    tol = 64 * sp["floor"]
    if d2 > bound + tol:
        fails.append(f"certificate violated ({sp['p']},{sp['q']}) "
                     f"om={sp['om']:+.0f}")
    pred_slack = (bound - d2) / 2
    if abs(pred_slack - worst_slack) > tol:
        fails.append(f"closed form misses the slack ({sp['p']},{sp['q']}) "
                     f"om={sp['om']:+.0f}: {pred_slack:.6f} vs "
                     f"{worst_slack:.6f}")
    for val, name, den in [(2160 * w ** 2, "w^2", 2160),
                           (2304 * delta ** 2, "delta^2", 2304)]:
        if abs(val - round(val)) > den * 4 * tol:
            fails.append(f"{den}*{name} not integer ({sp['p']},{sp['q']}) "
                         f"om={sp['om']:+.0f}: {val:.6f}")
    is_tight = abs(bound - d2) <= tol
    if is_tight:
        tight.append((sp["p"], sp["q"], int(round(sp["om"])), edge))
    else:
        all_slacks.append(worst_slack)
    corner2 = dim == 2 and sp["p"] in (1, 5) and sp["q"] in (1, 5)
    if corner2 and not delta + w > bound:
        fails.append(f"Weyl form unexpectedly suffices ({sp['p']},{sp['q']})"
                     f" om={sp['om']:+.0f}")
    log(f"    ({sp['p']},{sp['q']}) om={sp['om']:+.0f} dim={dim} "
        f"[{smin},{smax}] slack {worst_slack:+.6f} = closed form "
        f"{pred_slack:+.6f}  w={w:.6f} delta={delta:.6f} "
        f"2160w2={round(2160 * w ** 2):3d} 2304del2={round(2304 * delta ** 2):3d} "
        f"bound={bound:.4f}"
        f"{'  TIGHT(' + edge + ')' if is_tight else ''}")
want_tight = sorted([(1, 1, -2, "upper"), (1, 1, 2, "upper"),
                     (1, 5, -2, "lower"), (1, 5, 2, "lower"),
                     (5, 5, -2, "upper"), (5, 5, 2, "upper")])
if sorted(tight) != want_tight:
    fails.append(f"tight set {sorted(tight)} != corner six with edges")
log(f"    tight set: {sorted(tight)}")
log("    (gate: exactly the corner six, upper on (1,1)/(5,5), lower on (1,5))")
# the slack landscape, gated as counts and closed-form values
corner_slack = (0.25 - np.sqrt(1 / 144 + 1 / 27)) / 2
n_corner2 = sum(1 for sl in all_slacks if abs(sl - corner_slack) < 1e-9)
n_rest = sum(1 for sl in all_slacks if sl > 0.1)
rest_min = min((sl for sl in all_slacks if sl > 0.1), default=0.0)
log(f"    slack landscape: 6 tight, {n_corner2} at the corner dim-2 value "
    f"(0.25 - sqrt(1/144 + 1/27))/2 = {corner_slack:.6f}, {n_rest} at >= "
    f"0.12 (min {rest_min:.6f}); gate: 6 + {n_corner2} + {n_rest} == 48")
if n_corner2 != 6 or n_rest != 36 or rest_min < 0.12:
    fails.append(f"slack landscape changed: {n_corner2} corner, {n_rest} "
                 f"rest, min {rest_min:.4f}")

# ------------------------- (F) independent cells-route bridge per space
log()
log("(F) BRIDGE. Independent route (ad_H on the raw (p,q) cell block, no")
log("    dyads): the eigenspace at each omega gives eig(comp N_A) for the")
log("    transversal A = {0,1,2}; gate: matches the dyad route to 64*floor")
log("    on every breaking space.")
worst_bridge = 0.0
for pq in sorted(set((sp["p"], sp["q"]) for sp in spaces)):
    p, q = pq
    kets = [a for a in range(1 << n) if popcount(a) == p]
    bras = [b for b in range(1 << n) if popcount(b) == q]
    cells = [(a, b) for a in kets for b in bras]
    idx = {c: i for i, c in enumerate(cells)}
    m = len(cells)
    Acell = np.zeros((m, m))
    for (a, b), i in idx.items():
        for c in kets:
            if H[a, c] != 0:
                Acell[i, idx[(c, b)]] += H[a, c]
        for d_ in bras:
            if H[b, d_] != 0:
                Acell[i, idx[(a, d_)]] -= H[b, d_]
    NAv = np.array([sum(1.0 for l in [0, 1, 2] if (a ^ b) & (1 << l))
                    for (a, b) in cells])
    wv, V = np.linalg.eigh(Acell)
    for sp in [s2 for s2 in spaces if (s2["p"], s2["q"]) == pq]:
        sel = np.abs(wv - sp["om"]) < 1e-9
        P = V[:, sel]
        ev_cells = np.linalg.eigvalsh(P.T @ (NAv[:, None] * P))
        NA_dyad = sum(sp["comp"][l] for l in [0, 1, 2])
        ev_dyad = np.linalg.eigvalsh(NA_dyad)
        if P.shape[1] != len(sp["dy"]):
            fails.append(f"bridge dimension mismatch ({p},{q}) "
                         f"om={sp['om']:+.0f}")
            continue
        worst_bridge = max(worst_bridge,
                           float(np.max(np.abs(ev_cells - ev_dyad))))
        if np.max(np.abs(ev_cells - ev_dyad)) > 64 * sp["floor"]:
            fails.append(f"cells-route bridge broken ({p},{q}) "
                         f"om={sp['om']:+.0f}")
log(f"    worst spectrum deviation across all 48: {worst_bridge:.1e}  "
    + ("ok" if not any("bridge" in f for f in fails) else "FAIL"))

log()
log("=" * 78)
if fails:
    log("VERDICT: FAIL: " + "; ".join(fails))
else:
    log("VERDICT: all checks pass")
with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
