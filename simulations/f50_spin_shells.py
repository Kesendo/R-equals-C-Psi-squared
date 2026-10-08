"""F50, the sector's spin structure: W_p is one reduced matrix per N.

The object (docs/proofs/PROOF_WEIGHT1_DEGENERACY.md, the paragraph "The
sector's spin structure"): on the uniform Heisenberg chain, sector of p
flipped spins, eigenstates |E_a> of H_p, the matrix W_p with
W_ab = sum_l |<E_a|n_l|E_b>|^2. Its eigenvalues mu give the limiting heights
2(p - mu); the open statement asks mu_2 < p - 1/2 at every N and p.

Every eigenstate of the SU(2)-invariant chain carries a total spin S, and the
populations sort into spin shells. T_s is PROOF_UNIFORM_LAW's orbit sum, the
sum of the cells of Hamming s in block (p,p), s = 0, 2, ..., 2 min(p, N-p).

  G0  (law)   the zeros: W_ab = 0 at |S_a - S_b| >= 2 in every sector, and
      W_ab = 0 for a != b in one shell at half filling (THE_TWO_SPIN_ZEROS).
  G1  (exact, integers)   4[H_p, T_s] == 0 entry for entry (Lemma A1 on the
      chain). Controls: T_s with one symmetric pair of cells removed does not
      commute (s >= 2); at Delta = 2 every T_s with s >= 2 stops commuting
      except T_N at half filling, which is X^N and commutes with the XXZ chain.
  G2  (law)   Wigner-Eckart across sectors: two multiplets present in sectors
      p and p' satisfy W_ab(p) CG^2(p') = W_ab(p') CG^2(p), CG^2 the squared
      Clebsch-Gordan factor of the q = 0 vector component at weight
      m = p - N/2. Every multiplet of the smaller sector is matched exactly once
      (asserted). Error model: an eigenvector of a simple level is off by
      eps |H| / (its local gap) at first order, so the residual is held against
      eps |H| sum(1/gap) over the four levels of the pair and the ratio must
      stay below 8; it is 0.01 .. 0.1 in the large-gap cells and the bound is
      loose by about 100x in the near-degenerate ones (the per-decade maxima
      are printed), so a defect confined to those cells below ~1e-10 would
      pass. The mutations "swap up and down", "drop m^2" and "m -> m + 1" all
      fail by 1e12.
  G3  (law)   the shell sums 1_S (the spin-S projector in the sector, Lemma
      A1's M_S) span a W_p-invariant subspace; on it W_p has the eigenvalues
      p - s/2 (heights 0, 2, ..., 2 min(p, N-p)); mu_2 > p - 1, so the mu_2
      eigenspace lies in the shell-traceless complement (zero sum on every
      shell), which is then a corollary and is read as such.
  G4  (law + reading)   at half filling the diagonal is N/4, the spectrum is
      mirrored mu <-> N/2 - mu, (-1)^S is in the kernel, and the second-
      smallest eigenvalue (the open statement's face there) is read, above 1/2.
  G5  (reading)   at fixed N the smallest non-stationary height falls with p
      toward half filling.

Rounding model for G3 and G4: eps sqrt(dim) p, the backward error of a dense
symmetric eigensolver on a matrix of norm p; the first-order eigenvector error
eps |H| / gap does not enter these rows, because a rotation inside a
near-degenerate pair of ONE spin leaves the shell sums, the X^N parity and so
the diagonal 1/2 unchanged, and a pair of different spins is caught by the
sharpness check. The ratios are printed and must stay below 8 (read 0.04 .. 7.2
to N = 12, a scale read off the runs, not derived). G0's zeros are squared
matrix elements and carry the squared eigenvector error, (eps |H| sum 1/gap)^2 N.
The simple-spectrum check asks every gap > 1000 eps sqrt(dim) |H|, so that the
first-order eigenvector error is below 1e-3 everywhere.

Run:  python simulations/f50_spin_shells.py [--nmax 12]
Must print "ALL GATES PASS". About 50 s at --nmax 12. Prints ASCII only.
"""
import argparse
import sys
from itertools import combinations

import numpy as np

EPS = np.finfo(float).eps
RATIO_CAP = 8.0
FAIL = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


# ----------------------------------------------------------------- sector
def sector_basis(N, p):
    states = sorted(sum(1 << l for l in occ) for occ in combinations(range(N), p))
    return states, {x: i for i, x in enumerate(states)}


def four_h_int(N, p, states, idx, delta=1):
    """4 H_p over the integers, H = sum_bonds (XX+YY+ZZ)/4 with ZZ -> delta ZZ."""
    d = len(states)
    H4 = np.zeros((d, d), dtype=np.int64)
    for i, x in enumerate(states):
        for l in range(N - 1):
            bl, bm = (x >> l) & 1, (x >> (l + 1)) & 1
            H4[i, i] += delta if bl == bm else -delta
            if bl != bm:
                H4[idx[x ^ (1 << l) ^ (1 << (l + 1))], i] += 2
    return H4


def s2_sector(N, states, idx):
    d = len(states)
    S2 = np.zeros((d, d))
    for i, x in enumerate(states):
        S2[i, i] += 0.75 * N
        for a in range(N):
            for b in range(a + 1, N):
                if ((x >> a) & 1) == ((x >> b) & 1):
                    S2[i, i] += 0.5
                else:
                    S2[i, i] -= 0.5
                    S2[idx[x ^ (1 << a) ^ (1 << b)], i] += 1.0
    return S2


def hamming_matrix(states):
    x = np.array(states)
    return np.array([[bin(a ^ b).count("1") for b in x] for a in x])


def orbit_sum(states, s):
    """T_s: 1 on every cell of Hamming s, else 0 (integers)."""
    return (hamming_matrix(states) == s).astype(np.int64)


def commutes(H4, T):
    return bool(np.all(H4 @ T - T @ H4 == 0))


def spin_resolved(N, p):
    states, idx = sector_basis(N, p)
    H = four_h_int(N, p, states, idx) / 4.0
    S2 = s2_sector(N, states, idx)
    E, V = np.linalg.eigh(H)
    Sval = np.einsum("ia,ij,ja->a", V, S2, V)
    S = np.round((np.sqrt(4 * Sval + 1) - 1) / 2, 6)
    srt = np.sort(E)
    lgap = np.array([np.min(np.abs(srt[srt != e] - e)) for e in E])
    resid = np.linalg.norm(S2 @ V - V * Sval, axis=0)
    hnorm = np.linalg.norm(H, 2)
    sharp = np.max(resid * lgap) / (EPS * np.sqrt(len(E)) * np.linalg.norm(S2, 2) * hnorm)
    return states, idx, H, E, V, S, sharp, lgap, hnorm


def w_matrix(N, V, states):
    W = np.zeros((V.shape[1],) * 2)
    for l in range(N):
        nl = np.array([(x >> l) & 1 for x in states], dtype=float)
        M = V.T @ (nl[:, None] * V)
        W += M ** 2
    return W


def cg_sq(Sa, Sb, m):
    """|<S_b m; 1 0 | S_a m>|^2, the q = 0 vector-operator factor (S_b the ket)."""
    if abs(Sa - Sb) > 1 + 1e-9:
        return 0.0
    if abs(Sa - Sb) < 1e-9:
        return 0.0 if Sb < 1e-9 else m * m / (Sb * (Sb + 1))
    if Sa > Sb:
        S = Sb
        return ((S + 1) ** 2 - m * m) / ((S + 1) * (2 * S + 1))
    S = Sb
    return (S * S - m * m) / (S * (2 * S + 1))


def model(d, p):
    return EPS * np.sqrt(d) * p


# ------------------------------------------------------------------ gates
def gate_g0(N, p, S, W, lgap, hnorm):
    # a zero of W_ab is a squared matrix element: its rounding is the square of the
    # eigenvector error, (eps |H| (1/gap_a + 1/gap_b))^2, summed over the N sites
    d = W.shape[0]
    q = (EPS * hnorm * (1 / lgap[:, None] + 1 / lgap[None, :])) ** 2 * N
    far = np.abs(S[:, None] - S[None, :]) >= 2 - 1e-9
    worst = np.max(np.abs(W[far]) / q[far]) if far.any() else 0.0
    check(f"G0 N={N} p={p}: W_ab = 0 at |S_a - S_b| >= 2 ({far.sum()} cells), max/model = {worst:.3f}",
          worst < RATIO_CAP)
    if 2 * p == N:
        same = (np.abs(S[:, None] - S[None, :]) < 1e-9) & ~np.eye(d, dtype=bool)
        worst = np.max(np.abs(W[same]) / q[same])
        check(f"G0 N={N} half filling: W_ab = 0 within a shell ({same.sum()} cells), max/model = {worst:.3f}",
              worst < RATIO_CAP)


def gate_g1(N, p):
    states, idx = sector_basis(N, p)
    H4 = four_h_int(N, p, states, idx)
    smax = 2 * min(p, N - p)
    classes = list(range(0, smax + 1, 2))
    check(f"G1 N={N} p={p}: 4[H_p,T_s] == 0 for s in {classes} (integers)",
          all(commutes(H4, orbit_sum(states, s)) for s in classes))
    # control 1: remove one symmetric pair of cells from T_s (s >= 2)
    broken = []
    for s in classes[1:]:
        T = orbit_sum(states, s)
        i, j = np.argwhere(T == 1)[0]
        T[i, j] = T[j, i] = 0
        broken.append(not commutes(H4, T))
    check(f"G1 control N={N} p={p}: T_s minus one symmetric cell pair does not commute, s >= 2", all(broken))
    # control 2: Delta = 2 breaks every s >= 2 except T_N at half filling
    H4d = four_h_int(N, p, states, idx, delta=2)
    actual = [s for s in classes if not commutes(H4d, orbit_sum(states, s))]
    expect = [s for s in classes if s >= 2 and not (s == N)]
    check(f"G1 control N={N} p={p}: at Delta = 2 the non-commuting classes are {expect}", actual == expect,
          f"found {actual}")


def gate_g2(N, res):
    worst = 0.0
    count = skipped = 0
    decades = {}
    ps = sorted(res)
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            p1, p2 = ps[i], ps[j]
            E1, S1, W1, m1, g1v, h1 = res[p1]
            E2, S2_, W2, m2, g2v, h2 = res[p2]
            lookup = {}
            for k in range(len(E2)):
                lookup.setdefault((round(E2[k], 6), S2_[k]), []).append(k)
            match = []
            for a in range(len(E1)):
                ka = lookup.get((round(E1[a], 6), S1[a]), [])
                if len(ka) != 1:
                    skipped += 1
                    match.append(None)
                else:
                    match.append(ka[0])
            for a in range(len(E1)):
                if match[a] is None:
                    continue
                for b in range(a + 1, len(E1)):
                    if match[b] is None:
                        continue
                    g1, g2 = cg_sq(S1[a], S1[b], m1), cg_sq(S1[a], S1[b], m2)
                    if g1 == 0.0 or g2 == 0.0:
                        continue
                    inv = 1 / g1v[a] + 1 / g1v[b] + 1 / g2v[match[a]] + 1 / g2v[match[b]]
                    r = abs(W1[a, b] * g2 - W2[match[a], match[b]] * g1) / (EPS * max(h1, h2) * inv)
                    worst = max(worst, r)
                    dec = int(np.floor(np.log10(inv)))
                    decades[dec] = max(decades.get(dec, 0.0), r)
                    count += 1
    check(f"G2 N={N}: every multiplet of the smaller sector matched exactly once (skipped {skipped})", skipped == 0)
    per = ", ".join(f"1e{k}: {v:.3f}" for k, v in sorted(decades.items()))
    check(f"G2 N={N}: Wigner-Eckart across sectors, {count} cell pairs, max residual/model = {worst:.3f}",
          count > 0 and worst < RATIO_CAP, f"per decade of sum(1/gap): {per}")


def gate_g3(N, p, E, S, W):
    d = len(E)
    shells = sorted(set(S))
    B = np.stack([(S == s).astype(float) for s in shells], axis=1)
    Q, _ = np.linalg.qr(B)
    resid = np.linalg.norm(W @ Q - Q @ (Q.T @ (W @ Q)), 2)
    check(f"G3 N={N} p={p}: shell sums span a W-invariant subspace, residual/model = {resid/model(d,p):.2f}",
          resid / model(d, p) < RATIO_CAP)
    ev = np.sort(np.linalg.eigvalsh(Q.T @ W @ Q))
    expect = np.sort([p - s / 2 for s in range(0, 2 * min(p, N - p) + 1, 2)])
    dev = np.max(np.abs(ev - expect))
    check(f"G3 N={N} p={p}: block eigenvalues are p - s/2, s = 0, 2, ..., {2*min(p,N-p)}, deviation/model = {dev/model(d,p):.2f}",
          dev / model(d, p) < RATIO_CAP)
    mu, C = np.linalg.eigh(W)
    o = np.argsort(-mu)
    mu2 = mu[o[1]]
    check(f"G3 N={N} p={p}: mu_2 = {mu2:.6f} > p - 1 (margin {mu2-(p-1):.4f}), so the mu_2 eigenspace is shell-traceless",
          mu2 - (p - 1) > RATIO_CAP * model(d, p))
    space = C[:, np.abs(mu - mu2) < RATIO_CAP * model(d, p)]
    sums = np.max(np.abs(B.T @ space))
    check(f"G3 N={N} p={p}: the mu_2 eigenspace (dim {space.shape[1]}) is shell-traceless, max shell sum/model = {sums/model(d,p):.2f}",
          sums / model(d, p) < RATIO_CAP)
    return 2 * (p - mu2)


def gate_g4(N, S, W):
    d = W.shape[0]
    p = N // 2
    mu = np.sort(np.linalg.eigvalsh(W))
    diag_dev = np.max(np.abs(W.diagonal() - N / 4))
    check(f"G4 N={N} half filling: diagonal N/4, deviation/model = {diag_dev/model(d,p):.2f}",
          diag_dev / model(d, p) < RATIO_CAP)
    mirror = np.max(np.abs(mu - (N / 2 - mu[::-1])))
    check(f"G4 N={N} half filling: spectrum mirrored mu <-> N/2 - mu, residual/model = {mirror/model(d,p):.2f}",
          mirror / model(d, p) < RATIO_CAP)
    par = (-1.0) ** np.round(S)
    par /= np.linalg.norm(par)
    kern = np.linalg.norm(W @ par)
    check(f"G4 N={N} half filling: (-1)^S is in the kernel, |W v|/model = {kern/model(d,p):.2f}",
          kern / model(d, p) < RATIO_CAP)
    check(f"G4 N={N} half filling: second-smallest eigenvalue {mu[1]:.6f} > 1/2 (reading)", mu[1] > 0.5)


def main(nmax):
    for N in range(4, nmax + 1):
        res = {}
        heights = {}
        for p in range(1, N // 2 + 1):
            states, idx, H, E, V, S, sharp, lgap, hnorm = spin_resolved(N, p)
            d = len(states)
            gap = np.min(lgap)
            floor = 1000 * EPS * np.sqrt(d) * hnorm
            check(f"sector N={N} p={p}: simple spectrum (dim {d}, smallest gap {gap:.2e} > {floor:.1e}), "
                  f"S sharp on every eigenvector (residual gap_a / (eps sqrt(dim) |S^2| |H|) = {sharp:.2f})",
                  gap > floor and sharp < RATIO_CAP)
            W = w_matrix(N, V, states)
            res[p] = (E, S, W, p - N / 2, lgap, hnorm)
            gate_g0(N, p, S, W, lgap, hnorm)
            gate_g1(N, p)
            heights[p] = gate_g3(N, p, E, S, W)
            if 2 * p == N:
                gate_g4(N, S, W)
        gate_g2(N, res)
        hs = [heights[p] for p in sorted(heights)]
        check(f"G5 N={N}: smallest height falls with p toward half filling (reading) " +
              ", ".join(f"{h:.4f}" for h in hs), all(a > b for a, b in zip(hs, hs[1:])))
    print()
    if FAIL:
        print(f"{len(FAIL)} GATE(S) FAILED:")
        for f in FAIL:
            print("   " + f)
        sys.exit(1)
    print("ALL GATES PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=12)
    main(ap.parse_args().nmax)
