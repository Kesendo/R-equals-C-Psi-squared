"""Gate for docs/proofs/PROOF_UNIFORM_LAW.md: the Delta = 1 uniform law
comp(N_l) = (s/N) Id on every class-pure slice Omega_s, derived in two layers.

Sections, matching the proof:
  A. Lemma A2 as an INTEGER identity (no eigensolver, no tolerance): on every
     block (p,q) and class s, N * #(class-s cells with bit l in a XOR b)
     == s * #(class-s cells), for every l.  N = 2..10.
  B. Layer B0 exact facts, N = 5..8: the sector-1 restriction of the chain H
     equals (N-1)J*Id - 2J*L_path ENTRY-EXACTLY (both sides are small exact
     floats, compared with ==); the cosine-mode densities obey
     d_m + d_{N-m} = 2/N for m = 1..N-1 and d_0 = 1/N (plus d_{N/2} = 1/N at
     even N), gated against the error-model floor eps*||M||/gap. Plus the B2
     identity itself at N = 6 and 8: the proof's expansion of <y, Zcal_l x>
     (diagonal bracket + factored cross term) against the direct trace, on
     GENERIC dyad-space vectors, not only on the kernel.
  C. The census, N = 5..8, all blocks of the uniform chain: every Omega_s
     row contains T_s (projection norm 1 within floor), every row is
     span{T_s} (dim 1) except the three corner rows (1,1) s=2, (1,N-1) s=N-2,
     (N-1,N-1) s=2 with dim floor(N/2); law residual within 64 * floor on
     every row; no class-pure vector inside any degenerate eigenspace at
     omega != 0; the enhanced (1,1) space equals the predicted pair-equality
     span. The purity cutoff 1e-9 is reported as a measured separation (the
     nearest kept/dropped Gram eigenvalues per N), not trusted blindly.
  D. Controls: a random SWAP graph keeps the law and loses the enhancement
     (ALL blocks at N = 6, corner blocks at N = 8, every row dim 1);
     Delta = 0.5 violates on (1,5) s=4 (break gated from below at 0.05).

Error model (case 2 of the no-rounding rule): eigenvector-dependent residuals
are gated at 64 * eps * ||A|| / gap per eigenspace, gap the distance to the
nearest distinct ad_H eigenvalue. Exact routes (sections A and the matrix
identity in B) are compared exactly.

Writes simulations/results/uniform_law_gate.txt; last line is the VERDICT.
"""

import numpy as np
from math import comb

eps = np.finfo(float).eps
TOL_GROUP = 1e-9
OUT = "simulations/results/uniform_law_gate.txt"
_lines = []
fails = []


def log(msg=""):
    print(msg)
    _lines.append(msg)


def popcount(x):
    return bin(x).count("1")


def build_H_swap_graph(n, bonds):
    """H = Sum w * (XX + YY + ZZ) on the given bonds (Pauli book, J inside w)."""
    d = 1 << n
    H = np.zeros((d, d))
    for (i, j, w) in bonds:
        for a in range(d):
            za = 1 - 2 * ((a >> i) & 1)
            zb = 1 - 2 * ((a >> j) & 1)
            H[a, a] += w * za * zb
            if ((a >> i) & 1) != ((a >> j) & 1):
                H[a, a ^ (1 << i) ^ (1 << j)] += 2 * w
    return H


def build_H_xxz(n, J, Delta):
    d = 1 << n
    H = np.zeros((d, d))
    for a in range(d):
        for l in range(n - 1):
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a, a] += J * Delta * za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2 * J
    return H


def chain_bonds(n, J=1.0):
    return [(l, l + 1, J) for l in range(n - 1)]


# ---------------------------------------------------------------- section A
log("A. Lemma A2 as an integer identity, N = 2..10, all blocks and classes")
worst_n = 0
for n in range(2, 11):
    ok_all = True
    for p in range(n + 1):
        for q in range(p, n + 1):
            counts = {}
            for a in range(1 << n):
                if popcount(a) != p:
                    continue
                for b in range(1 << n):
                    if popcount(b) != q:
                        continue
                    s = popcount(a ^ b)
                    tot, per_l = counts.setdefault(s, [0, [0] * n])
                    counts[s][0] += 1
                    for l in range(n):
                        if (a ^ b) & (1 << l):
                            counts[s][1][l] += 1
            for s, (tot, per_l) in counts.items():
                for l in range(n):
                    if n * per_l[l] != s * tot:
                        ok_all = False
    log(f"   N={n}: N*count_l == s*total on every (p,q,s,l): "
        f"{'exact' if ok_all else 'FAIL'}")
    if not ok_all:
        fails.append(f"integer double count N={n}")

# ---------------------------------------------------------------- section B
log()
log("B. Layer B0: one-magnon exact facts of the uniform chain")
for n in (5, 6, 7, 8):
    J = 1.0
    H = build_H_swap_graph(n, chain_bonds(n, J))
    cfg1 = [a for a in range(1 << n) if popcount(a) == 1]
    site = [int(np.log2(a)) for a in cfg1]
    ordk = np.argsort(site)
    M1 = H[np.ix_(cfg1, cfg1)][np.ix_(ordk, ordk)]
    L = np.zeros((n, n))
    for l in range(n - 1):
        L[l, l] += 1
        L[l + 1, l + 1] += 1
        L[l, l + 1] -= 1
        L[l + 1, l] -= 1
    target = (n - 1) * J * np.eye(n) - 2 * J * L
    exact = bool(np.all(M1 == target))
    log(f"   N={n}: sector-1 H == (N-1)J*Id - 2J*L_path entry-exactly: "
        f"{exact}")
    if not exact:
        fails.append(f"sector-1 identity N={n}")
    E1, V1 = np.linalg.eigh(M1)
    gap = np.diff(np.sort(E1)).min()
    floor = eps * np.abs(E1).max() / gap
    # label columns by cosine mode m
    labels = np.empty(n, dtype=int)
    for col in range(n):
        best, bm = 0.0, -1
        for m in range(n):
            c = np.cos(np.pi * m * (2 * np.arange(n) + 1) / (2 * n))
            ov = abs(np.dot(c / np.linalg.norm(c), V1[:, col]))
            if ov > best:
                best, bm = ov, m
        labels[col] = bm
    d = V1 ** 2
    col_of = {m: int(np.where(labels == m)[0][0]) for m in range(n)}
    pair_dev = max(np.abs(d[:, col_of[m]] + d[:, col_of[n - m]] - 2 / n).max()
                   for m in range(1, (n + 1) // 2))
    flat_devs = [np.abs(d[:, col_of[0]] - 1 / n).max()]
    if n % 2 == 0:
        flat_devs.append(np.abs(d[:, col_of[n // 2]] - 1 / n).max())
    flat_dev = max(flat_devs)
    ok = pair_dev <= 64 * floor and flat_dev <= 64 * floor
    log(f"        d_m + d_(N-m) - 2/N: {pair_dev:.2e}, flat modes - 1/N: "
        f"{flat_dev:.2e} (floor {floor:.2e})  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"cosine pairing N={n}")

# B2 identity on GENERIC dyad-space vectors (not only the kernel): the
# expansion  <y, Zcal_l x> = sum_m ybar_m x_m (2 psi_m(l)^2 + 2 chi_m(l)^2 - 1)
#            - 4 (sum_m ybar_m psi_m(l) chi_m(l)) (sum_m x_m psi_m(l) chi_m(l))
# against the direct trace tr(Y^dag Z_l X Z_l) on block (1, N-1).
log()
log("B2. identity check on generic dyad vectors, block (1,N-1)")
for n in (6, 8):
    H = build_H_swap_graph(n, chain_bonds(n))
    cfg1 = [a for a in range(1 << n) if popcount(a) == 1]
    cfgq = [b for b in range(1 << n) if popcount(b) == n - 1]
    site1 = np.argsort([int(np.log2(a)) for a in cfg1])
    siteq = np.argsort([int(np.log2(((1 << n) - 1) ^ b)) for b in cfgq])
    E1, V1 = np.linalg.eigh(H[np.ix_(cfg1, cfg1)])
    Eq, Vq = np.linalg.eigh(H[np.ix_(cfgq, cfgq)])
    psi = V1[site1, :]
    chi = np.zeros_like(psi)
    for m in range(n):
        j = int(np.argmin(np.abs(Eq - E1[m])))
        chi[:, m] = Vq[siteq, j]
    rng = np.random.default_rng(1)
    x = rng.normal(size=n)
    y = rng.normal(size=n)
    worst = 0.0
    for l in range(n):
        lhs = sum(y[m] * x[m] * (2 * psi[l, m] ** 2 + 2 * chi[l, m] ** 2 - 1)
                  for m in range(n))
        lhs -= 4 * (sum(y[m] * psi[l, m] * chi[l, m] for m in range(n))
                    * sum(x[m] * psi[l, m] * chi[l, m] for m in range(n)))
        # direct: X = sum_m x_m |psi_m><chi_m| as a cells matrix, Z_l diag
        X = sum(x[m] * np.outer(psi[:, m], chi[:, m]) for m in range(n))
        Y = sum(y[m] * np.outer(psi[:, m], chi[:, m]) for m in range(n))
        # repo convention: Z = -1 on an occupied site. ket |e_i>: occupied at
        # i, so Z_l eigenvalue = -1 if l == i else +1
        zket = np.array([-1.0 if i == l else 1.0 for i in range(n)])
        # bra |1bar - e_j>: all occupied except j: Z_l = +1 if l == j else -1
        zbra = np.array([1.0 if j == l else -1.0 for j in range(n)])
        rhs = float(np.sum(Y * (zket[:, None] * X * zbra[None, :])))
        worst = max(worst, abs(lhs - rhs))
    ok = worst <= 1e-12
    log(f"   N={n}: max |expansion - direct| over l = {worst:.2e}  "
        f"{'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"B2 identity N={n}")


# ---------------------------------------------------------------- section C
def sector(H, n, p):
    cfg = [a for a in range(1 << n) if popcount(a) == p]
    E, V = np.linalg.eigh(H[np.ix_(cfg, cfg)])
    return cfg, E, V


def block_rows(n, H, p, q):
    """Per degenerate ad_H group and class: (omega, s, dim, law_res, T_in,
    floor, basis B). Also returns count of class-pure rows at omega != 0."""
    kets, Ek, Vk = sector(H, n, p)
    bras, Eb, Vb = sector(H, n, q)
    nb = len(bras)
    cls = np.array([popcount(a ^ b) for a in kets for b in bras])
    Nl = [np.array([1.0 if (a ^ b) & (1 << l) else 0.0
                    for a in kets for b in bras]) for l in range(n)]
    om = (Ek[:, None] - Eb[None, :]).ravel()
    order = np.argsort(om)
    oms = om[order]
    diffs = np.diff(oms)
    if diffs.size:
        within = diffs[diffs <= TOL_GROUP]
        between = diffs[diffs > TOL_GROUP]
        if within.size:
            GSEP[0] = max(GSEP[0], float(within.max()))
        if between.size:
            GSEP[1] = min(GSEP[1], float(between.min()))
    rows, off0 = [], 0
    start = 0
    for t in range(1, len(order) + 1):
        if t == len(order) or oms[t] - oms[start] > TOL_GROUP:
            g = order[start:t]
            if len(g) >= 2:
                gap = np.inf
                if start > 0:
                    gap = min(gap, oms[start] - oms[start - 1])
                if t < len(order):
                    gap = min(gap, oms[t] - oms[t - 1])
                floor = eps * np.abs(om).max() / gap
                Vg = np.empty((len(cls), len(g)))
                for c, idx in enumerate(g):
                    i, j = divmod(idx, nb)
                    Vg[:, c] = np.kron(Vk[:, i], Vb[:, j])
                for s in sorted(set(cls)):
                    mask = (cls == s).astype(float)
                    W = Vg * mask[:, None]
                    evg, Ug = np.linalg.eigh(W.T @ W)
                    keep = evg > 1 - 1e-9
                    if keep.any():
                        SEP[0] = min(SEP[0], float(evg[keep].min()))
                    if (~keep).any():
                        SEP[1] = max(SEP[1], float(evg[~keep].max()))
                    dim = int(keep.sum())
                    if dim == 0:
                        continue
                    if abs(oms[start]) > TOL_GROUP:
                        off0 += 1
                        continue
                    B = Vg @ Ug[:, keep]
                    res = max(float(np.linalg.norm(
                        B.T @ (Nl[l][:, None] * B) - (s / n) * np.eye(dim), 2))
                        for l in range(n))
                    t_s = mask / np.sqrt(mask.sum())
                    t_in = float(np.linalg.norm(B.T @ t_s))
                    rows.append((oms[start], int(s), dim, res, t_in, floor, B))
            start = t
    return rows, off0


def predicted_enhanced_11(n, V1, labels):
    """Predicted (1,1) s=2 span: dyad coords with x_m = x_{N-m}, sum x = 0."""
    col_of = {m: int(np.where(labels == m)[0][0]) for m in range(n)}
    P_cols = []
    for m in range(1, (n + 1) // 2):
        v = np.zeros(n)
        v[col_of[m]] = 1
        v[col_of[n - m]] = 1
        P_cols.append(v)
    singles = ([n // 2, 0] if n % 2 == 0 else [0])
    for m in singles:
        v = np.zeros(n)
        v[col_of[m]] = 1
        P_cols.append(v)
    P = np.stack(P_cols, 1)
    ns_ = np.linalg.svd((np.ones(n) @ P)[None, :])[2][1:, :].T
    K, _ = np.linalg.qr(P @ ns_)
    D = np.empty((n * n, n))
    for A in range(n):
        D[:, A] = np.kron(V1[:, A], V1[:, A])
    return K, D


log()
log("C. Census, uniform chain, all blocks: T_s in every row, dims, law")
SEP = [np.inf, 0.0]     # [min kept, max dropped] class-purity Gram eigenvalue
GSEP = [0.0, np.inf]    # [max within-group, min between-group] omega spacing
for n in (5, 6, 7, 8):
    H = build_H_swap_graph(n, chain_bonds(n))
    corner = {(1, 1, 2), (1, n - 1, n - 2), (n - 1, n - 1, 2)}
    n_rows = 0
    worst_res_f = worst_tin = 0.0
    dims_ok = True
    off0_total = 0
    for p in range(n + 1):
        for q in range(p, n + 1):
            if comb(n, p) * comb(n, q) < 2:
                continue
            rows, off0 = block_rows(n, H, p, q)
            off0_total += off0
            for (wv, s, dim, res, t_in, floor, B) in rows:
                n_rows += 1
                worst_res_f = max(worst_res_f, res / (64 * floor))
                worst_tin = max(worst_tin, abs(t_in - 1.0))
                want = n // 2 if (p, q, s) in corner else 1
                if dim != want:
                    dims_ok = False
                    log(f"   N={n} ({p},{q}) s={s}: dim {dim} != expected "
                        f"{want}  FAIL")
                if (p, q, s) == (1, 1, 2):
                    cfg1 = [a for a in range(1 << n) if popcount(a) == 1]
                    site = [int(np.log2(a)) for a in cfg1]
                    M1 = H[np.ix_(cfg1, cfg1)][np.ix_(np.argsort(site),
                                                      np.argsort(site))]
                    E1, V1 = np.linalg.eigh(M1)
                    labels = np.empty(n, dtype=int)
                    for col in range(n):
                        best, bm = 0.0, -1
                        for m in range(n):
                            c = np.cos(np.pi * m * (2 * np.arange(n) + 1)
                                       / (2 * n))
                            ov = abs(np.dot(c / np.linalg.norm(c), V1[:, col]))
                            if ov > best:
                                best, bm = ov, m
                        labels[col] = bm
                    K, D = predicted_enhanced_11(n, V1, labels)
                    coords = D.T @ B
                    Q1, _ = np.linalg.qr(coords)
                    dist = float(np.linalg.norm(Q1 @ (Q1.T @ K) - K, 2))
                    ok = dist <= 64 * floor * 10
                    log(f"   N={n} (1,1) s=2: predicted pair-equality span == "
                        f"measured, distance {dist:.2e}  "
                        f"{'ok' if ok else 'FAIL'}")
                    if not ok:
                        fails.append(f"enhanced span match N={n}")
    ok = dims_ok and worst_res_f <= 1.0 and worst_tin <= 1e-12 \
        and off0_total == 0
    log(f"   N={n}: {n_rows} Omega_s rows, worst law res/floor(64) "
        f"{worst_res_f:.3f}, worst |T_s proj - 1| {worst_tin:.1e}, "
        f"class-pure rows off omega=0: {off0_total}, dims as predicted: "
        f"{dims_ok}  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"census N={n}")

log(f"   class-purity cutoff separation across the census: kept Gram "
    f"eigenvalues >= 1 - {1 - SEP[0]:.1e}, dropped <= {SEP[1]:.6f}")
log(f"   omega-grouping cutoff {TOL_GROUP:.0e} separation: within-group "
    f"spacings <= {GSEP[0]:.1e}, between-group >= {GSEP[1]:.1e}")

# ---------------------------------------------------------------- section D
log()
log("D. Controls")
for n in (6, 8):
    rng = np.random.default_rng(20260815 + n)
    bonds = [(i, j, rng.uniform(0.2, 2.0))
             for i in range(n) for j in range(i + 1, n) if rng.random() < 0.7]
    H = build_H_swap_graph(n, bonds)
    ok = True
    n_rows_d = 0
    if n == 6:
        block_list = [(p, q) for p in range(n + 1) for q in range(p, n + 1)
                      if comb(n, p) * comb(n, q) >= 2]
    else:
        block_list = [(1, 1), (1, n - 1), (n - 1, n - 1)]
    for (p, q) in block_list:
        rows, off0 = block_rows(n, H, p, q)
        for (wv, s, dim, res, t_in, floor, B) in rows:
            n_rows_d += 1
            if dim != 1 or res > 64 * floor or abs(t_in - 1) > 1e-12:
                ok = False
    scope = "all blocks" if n == 6 else "corner blocks"
    log(f"   N={n} random SWAP graph, {scope} ({n_rows_d} rows): all rows "
        f"dim 1 = span(T_s), law at floor: {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"random graph control N={n}")

H_ctl = build_H_xxz(6, 1.0, 0.5)
rows, _ = block_rows(6, H_ctl, 1, 5)
viol = max((res for (_, s, dim, res, _, _, _) in rows if s == 4), default=0.0)
ok = viol > 0.05
log(f"   N=6 Delta=0.5 control, (1,5) s=4: law residual {viol:.3f} "
    f"(gate > 0.05, break present): {'ok' if ok else 'FAIL'}")
if not ok:
    fails.append("Delta=0.5 control")

log()
if fails:
    log("VERDICT: FAIL: " + "; ".join(fails))
else:
    log("VERDICT: PASS: the uniform law's two layers hold as derived: the "
        "integer double count on the T_s line (N=2..10), the exact one-magnon "
        "facts and the B2 identity, the N=5..8 census (T_s in every row, "
        "corner dims floor(N/2), law at floor, nothing off omega=0), the "
        "enhanced-span match, and both controls")

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(_lines) + "\n")
