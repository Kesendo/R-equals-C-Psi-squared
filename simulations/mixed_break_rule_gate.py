"""THE TWO-ZERO RULE GATE: which parity-mixed eigenspaces at Delta = 1 break
C_l = 0, decided by two Clebsch-Gordan zeros, with closed-form residual sizes.

Answers the first half of arc `compressed_density_laws` NextStep (3). The
objects: on a coherence block (p, q) (convention p <= q throughout) the
eigenspace of ad_H at omega != 0 is spanned by multiplet dyads |i><j| with
E_i - F_j = omega, i and j labeled by (energy, total spin S, reflection
parity r), all (S, E, r) labels nondegenerate at N = 6 (gated in section A,
so every leg below is basis-unique up to sign). Every compressed matrix
element of N_l between dyads factors into one-sector density legs: for
(i, j), (i', j') with i != i' and j != j' it is -2 <i|n_l|i'> <j|n_l|j'>,
and for i = i' it is <j|n_l|j'> (1 - 2 <i|n_l|i>). Site occupation
n_l = (1 - Z_l)/2 in this file's convention (Z = +1 on an empty site).

THE RULE (two zeros): a leg <i|n_l|i'> is DEAD iff
  (a) |S_i - S_i'| >= 2: Z_l is the q = 0 component of a rank-1 spherical
      tensor (the site spin), so the Wigner-Eckart triangle kills the
      element (the identity half of n_l needs i = i'); or
  (b) the sector is half filling (magnetization m = 0) and S_i = S_i': the
      Clebsch-Gordan factor <S m; 1 0 | S m> = m/sqrt(S(S+1)) vanishes at
      m = 0. An equivalent route for this zero only: Z_l anticommutes with
      X^N, X^N maps the m = 0 sector to itself with eigenvalue
      i^N * (-1)^S on a spin-S state ((-1)^(S+1) at N = 6), so same-S
      states share the parity and an X^N-odd operator cannot connect them
      (that route kills every even Delta S at m = 0, more than the rule
      needs). The same zero makes the diagonal exact, <i|n_l|i> = 1/2,
      killing the i = i' factor (1 - 2 <i|n_l|i>).

TWO DIRECTIONS, only one a theorem. HOLD: if every cross-parity dyad pair
has a dead leg, C_l = 0 exactly (R-conjugation makes the C_l element of a
cross-parity pair 2x the comp(N_l) element and cancels same-parity pairs,
so the dead legs kill everything that survives the cancellation). BREAK:
"both legs alive" means "not forced zero by the two symmetries", and the
converse (some alive cross-parity pair => the space breaks) is GENERICITY,
measured 60/60 here, not derived: a leg alive by the rule could still
vanish at every common site, and in this census none does.

Verified: 60/60 mixed spaces over ALL fifteen N = 6 blocks p <= q,
element-exact (every rule-dead element measures at the floor). Four of the
60 break BELOW the endpoint gate's inherited 0.05 threshold
(resC = sqrt(1/432) ~ 0.0481: (1,2) and (1,4) at omega = +4, (2,5) and
(4,5) at omega = -4). No committed claim classified anything in that band
(the earlier census's blocks carry no break in (0.04, 0.05)); the honest
fragility is (2,4)'s smallest break 0.0609, 22% above 0.05. Classification
here is by the FLOOR: holds <= 64 * floor (measured <= 5.6e-16), breaks
>= 0.04, and anything in the nearly-fourteen-decade band between FAILS (the
VOID branch is the load-bearing fence; the 64 is slack the census never
uses).

SIZES (measured, not derived): every residual satisfies 2160 * resC^2
integer (gated), values {5, 8, 9, 20, 32, 36, 80, 128, 144}; the legs
entering these particular mixed spaces measure as square roots of
rationals, and the breaking elements are their products. Aliveness does
NOT imply rationality (generic alive legs have irrational squares), and
the denominator 2160 is a fit, not a derivation.

THE PROTECTION MAP (gate section C): each half-filled SIDE protects one
omega sign, bra side (q = N/2) the sign omega > 0, ket side (p = N/2) the
sign omega < 0; (3,3) has both sides and holds at all four omega; blocks
without an m = 0 sector hold nowhere. Mechanism, the derived half: at
m = 0 zero (b) leaves only Delta S = 1 legs alive, and every alive m = 0
leg in the census touches the S = N/2 ferromagnetic level (gated), whose
energy 5 is the global maximum, so its dyads reach only one omega sign per
side. That no other Delta S = 1 pair appears is census, not derivation.

Frontier: NO mixed spaces at Delta = 1 at any N in 3..8 except N = 6, on
ALL blocks p <= q per size: among N <= 8 the phenomenon lives at N = 6
alone (N = 2 has no omega != 0 dyad pairs to mix), and it
is a UNIFORM-chain resonance there (a generic palindromic bond profile
kills all 60 mixed spaces with R still a symmetry, asserted; section E).

Error model, not a chosen number: every rule-dead element and every hold
is exact in exact arithmetic; the float route is eigh's backward error,
element error <= C * eps * ||A|| / gap. The break fence 0.04 sits
thirteen decades above every floor printed here and a factor 1.2
below the smallest break. Eigenvalue grouping uses gtol = 1e-9; the dyad
difference spectrum's gap void is gated >= 1e-6 in every section, the
null-result sections D and E included (the thinnest margin is N = 8,
void 6.3e-6).

Writes simulations/results/mixed_break_rule_gate.txt; last line is the VERDICT.
"""

import numpy as np

OUT = "simulations/results/mixed_break_rule_gate.txt"
lines = []
fails = []
eps = np.finfo(float).eps

BLOCKS6 = [(p, q) for p in range(1, 6) for q in range(p, 6)]


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


def build_H(n, D=1.0, J=1.0):
    """Pauli book H = Sum_bonds J (XX + YY + D*ZZ) on the open chain."""
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


def build_S2(n):
    """S^2 = (3N/4) Id + (1/2) Sum_{i<j} sigma_i . sigma_j."""
    d = 1 << n
    S2 = np.zeros((d, d))
    for a in range(d):
        for i in range(n):
            for j in range(i + 1, n):
                za = 1 - 2 * ((a >> i) & 1)
                zb = 1 - 2 * ((a >> j) & 1)
                S2[a, a] += 0.5 * za * zb
                if ((a >> i) & 1) != ((a >> j) & 1):
                    S2[a, a ^ (1 << i) ^ (1 << j)] += 1.0
    S2 += (3 * n / 4) * np.eye(d)
    return S2


def sector_vectors(n, p, H, S2):
    """Orthonormal sector basis with sharp labels: diagonalize S^2 first
    (integer spectrum S(S+1), gaps >= 2), then H within each S-eigenspace,
    then R within each (S, E) level. Returns (list of dicts (E, S, r, v),
    cfg, floor) with floor = eps * ||Hp|| / (min level gap)."""
    cfg = [a for a in range(1 << n) if popcount(a) == p]
    idx = {a: k for k, a in enumerate(cfg)}
    Hp = H[np.ix_(cfg, cfg)]
    S2p = S2[np.ix_(cfg, cfg)]
    Rp = np.zeros((len(cfg), len(cfg)))
    for a in cfg:
        Rp[idx[revbits(a, n)], idx[a]] = 1.0
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
    floor = eps * float(np.linalg.norm(Hp, 2)) / min(gap_min, 2.0)
    return out, cfg, floor


def occ_all(cfg, vecs, n):
    """occ[l][i, i'] = <v_i| n_l |v_i'> over the labeled sector basis."""
    B = np.column_stack([d["v"] for d in vecs])
    occ = []
    for l in range(n):
        w = np.array([(a >> l) & 1 for a in cfg], dtype=float)
        occ.append(B.T @ (w[:, None] * B))
    return occ


def mixed_spaces(n, vp, vq, gtol=1e-9):
    """Group dyads by omega = E_i - F_j (sorted, split at gaps > gtol) and
    return the parity-mixed groups with omega != 0, plus the smallest
    inter-group gap (the void the grouping lives in)."""
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


def leg_dead(vecs, m, i, i2):
    """The two zeros. m is the sector magnetization (N - 2p)/2."""
    if i == i2:
        return m == 0
    dS = abs(vecs[i]["S"] - vecs[i2]["S"])
    return dS >= 2 or (m == 0 and dS == 0)


def comp_matrices(n, dy, op, oq):
    """Compressed N_l on the dyad space, all l, in the dyad basis."""
    dim = len(dy)
    comp = []
    for l in range(n):
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


def cells_route_res(n, p, q, H, gtol=1e-9):
    """Independent route: ad_H on the (p, q) cell block (endpoint-gate
    style), per-eigenspace worst spectral-norm residual of C_l, keyed by
    omega. Guards the dyad construction against being a different object."""
    kets = [a for a in range(1 << n) if popcount(a) == p]
    bras = [b for b in range(1 << n) if popcount(b) == q]
    cells = [(a, b) for a in kets for b in bras]
    idx = {c: i for i, c in enumerate(cells)}
    m = len(cells)
    A = np.zeros((m, m))
    for (a, b), i in idx.items():
        for c in kets:
            if H[a, c] != 0:
                A[i, idx[(c, b)]] += H[a, c]
        for d_ in bras:
            if H[b, d_] != 0:
                A[i, idx[(a, d_)]] -= H[b, d_]
    Nl = [np.array([1.0 if (a ^ b) & (1 << l) else 0.0 for (a, b) in cells])
          for l in range(n)]
    w, V = np.linalg.eigh(A)
    res = {}
    start = 0
    for k in range(1, len(w) + 1):
        if k == len(w) or w[k] - w[k - 1] > gtol:
            P = V[:, start:k]
            om = float(np.mean(w[start:k]))
            if P.shape[1] >= 2 and abs(om) > gtol:
                res[round(om, 6)] = max(float(np.linalg.norm(
                    P.T @ ((Nl[l] - Nl[n - 1 - l])[:, None] * P), 2))
                    for l in range(n // 2))
            start = k
    return res


log("=" * 78)
log("THE TWO-ZERO RULE GATE: mixed-space break/hold at Delta = 1, from two")
log("Clebsch-Gordan zeros, with closed-form residual sizes")
log("=" * 78)

# ---------------------------------------------- (A) sector spectroscopy, N=6
log()
log("(A) SECTOR SPECTROSCOPY, N = 6. Sharp (E, S, r) labels per sector; gate:")
log("    every label nondegenerate (legs basis-unique up to sign), S half-")
log("    integer-exact, R parities +-1 (asserted inside the builder).")
n = 6
H6 = build_H(n)
S26 = build_S2(n)
sector = {}
for p in range(1, 6):
    vecs, cfg, floor = sector_vectors(n, p, H6, S26)
    occ = occ_all(cfg, vecs, n)
    sector[p] = (vecs, cfg, occ, floor)
    labels = [(round(d["E"], 8), d["S"], d["r"]) for d in vecs]
    ndeg = len(labels) - len(set(labels))
    log(f"    sector {p}: {len(vecs)} states, degenerate (E,S,r) labels: "
        f"{ndeg}, floor {floor:.1e}  {'ok' if ndeg == 0 else 'FAIL'}")
    if ndeg != 0:
        fails.append(f"degenerate (E,S,r) label in sector {p}")

# -------------------------------------- (B) census + rule + element exactness
log()
log("(B) THE CENSUS AND THE RULE. Per mixed space: rule prediction vs the")
log("    floor-classified measurement, worst rule-dead element, and the")
log("    rational size 2160 * resC^2. Gate: 60 mixed spaces, 60/60 predicted,")
log("    dead elements <= 64 * floor, holds <= 64 * floor, breaks >= 0.04,")
log("    2160 * resC^2 integer for every break.")
n_mixed_total = 0
n_correct = 0
hold_map = {}
m0_alive_pairs = []
subfence = []
for (p, q) in BLOCKS6:
    vp, cfg_p, op, fp = sector[p]
    vq, cfg_q, oq, fq = sector[q]
    mp = (n - 2 * p) / 2
    mq = (n - 2 * q) / 2
    floor = fp + fq
    mixed, void = mixed_spaces(n, vp, vq)
    cells_res = cells_route_res(n, p, q, H6)
    if void < 1e-6:
        fails.append(f"({p},{q}): dyad difference void {void:.1e} < 1e-6")
    holds_at = []
    for om, dy, pars in mixed:
        n_mixed_total += 1
        comp = comp_matrices(n, dy, op, oq)
        resC = max(np.linalg.norm(comp[l] - comp[n - 1 - l], 2)
                   for l in range(n // 2))
        pred_break = False
        dead_max = 0.0
        for a, (i, j) in enumerate(dy):
            for b, (i2, j2) in enumerate(dy):
                if b <= a:
                    continue
                pdead = leg_dead(vp, mp, i, i2)
                qdead = leg_dead(vq, mq, j, j2)
                dead = pdead or qdead
                if mp == 0 and not pdead and i != i2:
                    m0_alive_pairs.append(tuple(sorted((vp[i]["S"],
                                                        vp[i2]["S"]))))
                if mq == 0 and not qdead and j != j2:
                    m0_alive_pairs.append(tuple(sorted((vq[j]["S"],
                                                        vq[j2]["S"]))))
                elem = max(abs(comp[l][a, b]) for l in range(n))
                if dead:
                    dead_max = max(dead_max, elem)
                elif pars[a] * pars[b] == -1:
                    pred_break = True
        if resC <= 64 * floor:
            meas = "hold"
            holds_at.append(om)
        elif resC >= 0.04:
            meas = "BREAK"
        else:
            meas = "VOID"
            fails.append(f"({p},{q}) omega={om:+.0f}: resC {resC:.2e} in the "
                         f"forbidden band (floor {floor:.1e})")
        ok = (pred_break and meas == "BREAK") or \
             (not pred_break and meas == "hold")
        if ok:
            n_correct += 1
        else:
            fails.append(f"rule mismatch ({p},{q}) omega={om:+.0f}")
        if dead_max > 64 * floor:
            fails.append(f"dead element alive ({p},{q}) omega={om:+.0f}: "
                         f"{dead_max:.2e}")
        # closed-form size: 2160 * resC^2 integer (breaks; holds give 0)
        k2160 = 2160 * resC ** 2
        tol_rat = 2160 * 2 * max(resC, 1.0) * 64 * floor
        rat_ok = abs(k2160 - round(k2160)) <= tol_rat
        if not rat_ok:
            fails.append(f"size not rational ({p},{q}) omega={om:+.0f}: "
                         f"2160*resC^2 = {k2160:.6f}")
        # bridge to the independent cells route
        cres = cells_res.get(round(om, 6))
        bridge_ok = cres is not None and abs(cres - resC) <= 64 * floor
        if not bridge_ok:
            fails.append(f"cells-route bridge ({p},{q}) omega={om:+.0f}")
        if 0.04 <= resC < 0.05:
            subfence.append((p, q, om))
        pairpar = f"({pars.count(1)}+,{pars.count(-1)}-)"
        log(f"    ({p},{q}) omega={om:+.0f} dim={len(dy)} {pairpar} "
            f"pred={'BREAK' if pred_break else 'hold '} meas={meas:5s} "
            f"resC={resC:.3e} 2160*resC^2={round(k2160):3d} "
            f"dead_max={dead_max:.1e} bridge={'ok' if bridge_ok else 'FAIL'}")
    hold_map[(p, q)] = holds_at
log(f"    census: {n_mixed_total} mixed spaces, rule correct on "
    f"{n_correct}/{n_mixed_total}")
log(f"    breaks below the inherited 0.05 fence: "
    f"{[(f'({p},{q})', f'{om:+.0f}') for (p, q, om) in subfence]} "
    f"(gate: exactly four)")
if len(subfence) != 4:
    fails.append(f"sub-fence break count {len(subfence)} != 4")
if n_mixed_total != 60:
    fails.append(f"census count {n_mixed_total} != 60")
if n_correct != n_mixed_total:
    fails.append("rule not 60/60")

# --------------------------------------------- (C) the m = 0 sign protection
log()
log("(C) THE HALF-FILLING PROTECTION. Each half-filled side protects one")
log("    omega sign; (3,3) has both sides and holds everywhere; blocks")
log("    without an m = 0 sector hold nowhere. Gate: the map.")
for (p, q) in BLOCKS6:
    holds = hold_map[(p, q)]
    half = n // 2
    if p == half and q == half:
        ok = len(holds) == 4
        desc = "both sides m=0, both signs protected"
    elif q == half:
        ok = len(holds) == 2 and all(om > 0 for om in holds)
        desc = "bra side m=0, protects omega > 0"
    elif p == half:
        ok = len(holds) == 2 and all(om < 0 for om in holds)
        desc = "ket side m=0, protects omega < 0"
    else:
        ok = holds == []
        desc = "no m=0 sector, no protection"
    log(f"    ({p},{q}): holds at {[f'{om:+.0f}' for om in holds]} "
        f"({desc})  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"protection map ({p},{q})")
# every alive leg recorded in an m = 0 sector touches the ferromagnetic
# S = N/2 level (the omega-sign mechanism's measured half)
bad_pairs = [sp for sp in m0_alive_pairs if max(sp) != n / 2]
log(f"    alive S-pairs recorded at m = 0 across the census: "
    f"{sorted(set(m0_alive_pairs))} (gate: every pair touches S = N/2)  "
    f"{'ok' if not bad_pairs else 'FAIL'}")
if bad_pairs:
    fails.append(f"m=0 alive pair without the S=N/2 level: {bad_pairs}")

# -------------------------------------------------- (D) frontier: N = 7, 8
log()
log("(D) FRONTIER. At Delta = 1 there are NO parity-mixed omega != 0")
log("    eigenspaces at any N in 3..8 except N = 6, on ANY block p <= q")
log("    (gate == 0): among N <= 8 the phenomenon lives at N = 6 alone.")
for nn in (3, 4, 5, 7, 8):
    Hn = build_H(nn)
    S2n = build_S2(nn)
    secn = {s: sector_vectors(nn, s, Hn, S2n) for s in range(1, nn)}
    blocks_n = [(p, q) for p in range(1, nn) for q in range(p, nn)]
    cnt = 0
    void_min = np.inf
    for (p, q) in blocks_n:
        mixed, vd = mixed_spaces(nn, secn[p][0], secn[q][0])
        cnt += len(mixed)
        void_min = min(void_min, vd)
    ok_v = void_min >= 1e-6
    log(f"    N={nn}: {cnt} mixed spaces across all {len(blocks_n)} blocks "
        f"p <= q, grouping void {void_min:.1e} (gated >= 1e-6)  "
        f"{'ok' if cnt == 0 and ok_v else 'FAIL'}")
    if cnt != 0:
        fails.append(f"unexpected mixed space at N={nn}")
    if not ok_v:
        fails.append(f"grouping void too small at N={nn}")

# ------------------------------- (E) uniform-chain resonance (bond control)
log()
log("(E) BOND CONTROL. The 60 mixed spaces are a UNIFORM-chain resonance,")
log("    as at Delta = 0 (endpoint gate, section B): a generic palindromic")
log("    bond profile kills every one (gate == 0), R staying a symmetry.")
Jb = [1.25, 0.75, 2.25, 0.75, 1.25]  # dyadic: the float sums are exact
d6 = 1 << n
Hpal = np.zeros((d6, d6))
for a in range(d6):
    for l in range(n - 1):
        za = 1 - 2 * ((a >> l) & 1)
        zb = 1 - 2 * ((a >> (l + 1)) & 1)
        Hpal[a, a] += Jb[l] * za * zb
        if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
            Hpal[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2 * Jb[l]
Rfull = np.zeros((d6, d6))
for a in range(d6):
    Rfull[revbits(a, n), a] = 1.0
r_com_pal = float(np.max(np.abs(Rfull @ Hpal - Hpal @ Rfull)))
if r_com_pal != 0.0:
    fails.append(f"[R, H_pal] != 0: {r_com_pal:.1e}")
cnt_pal = 0
void_pal = np.inf
sec_pal = {}
for (p, q) in BLOCKS6:
    for s in (p, q):
        if s not in sec_pal:
            sec_pal[s] = sector_vectors(n, s, Hpal, S26)
    mixed_pal, vd = mixed_spaces(n, sec_pal[p][0], sec_pal[q][0])
    cnt_pal += len(mixed_pal)
    void_pal = min(void_pal, vd)
if void_pal < 1e-6:
    fails.append(f"grouping void too small on the palindromic profile")
log(f"    palindromic J = {Jb}: {cnt_pal} mixed spaces over the fifteen "
    f"blocks, [R, H_pal] max entry {r_com_pal} (gated == 0.0), grouping "
    f"void {void_pal:.1e} (gated >= 1e-6)  "
    f"{'ok' if cnt_pal == 0 else 'FAIL'}")
if cnt_pal != 0:
    fails.append("mixed spaces survive a generic palindromic J at Delta = 1")
# one profile is not "generic": a seeded sweep of random DYADIC palindromic
# profiles (multiples of 1/16 in [0.5, 2], so the float sums stay exact and
# [R, H] == 0.0 holds exactly per profile)
rng = np.random.default_rng(20260816)
sweep_mixed = 0
for _ in range(6):
    half = rng.integers(8, 33, 3) / 16.0
    Jr = [half[0], half[1], half[2], half[1], half[0]]
    Hr = np.zeros((d6, d6))
    for a in range(d6):
        for l in range(n - 1):
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            Hr[a, a] += Jr[l] * za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                Hr[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2 * Jr[l]
    if float(np.max(np.abs(Rfull @ Hr - Hr @ Rfull))) != 0.0:
        fails.append("[R, H] != 0 in the random palindromic sweep")
    sec_r = {ss: sector_vectors(n, ss, Hr, S26) for ss in range(1, 6)}
    for (pp, qq) in BLOCKS6:
        mixed_r, vd = mixed_spaces(n, sec_r[pp][0], sec_r[qq][0])
        sweep_mixed += len(mixed_r)
        if vd < 1e-6:
            fails.append("grouping void too small in the random sweep")
log(f"    random dyadic palindromic sweep (6 profiles): {sweep_mixed} mixed "
    f"spaces total (gate == 0)  {'ok' if sweep_mixed == 0 else 'FAIL'}")
if sweep_mixed != 0:
    fails.append("mixed spaces survive the random palindromic sweep")

log()
log("=" * 78)
if fails:
    log("VERDICT: FAIL: " + "; ".join(fails))
else:
    log("VERDICT: all checks pass")
with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
