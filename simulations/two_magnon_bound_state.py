#!/usr/bin/env python3
"""Isolate the two-magnon bound state: the new object born when two excitations interact.

Free particles (XY, Delta=0) only scatter: the two excitations spread apart, no new object.
Turn on the ZZ interaction (Heisenberg/XXZ, Delta != 0) and a genuinely new eigenstate can
appear, one where the two excitations STAY TOGETHER, a bound pair that travels as a single
composite, split off in energy from the two-magnon scattering continuum.

We build H = sum_b (X_b X_{b+1} + Y_b Y_{b+1}) + Delta * sum_b Z_b Z_{b+1}, diagonalize, and in
the popcount-2 (two-excitation) sector measure for every eigenstate:
  - mean separation  <|i - j|>  of the two excited sites (small = bound, large = spread)
  - P(adjacent)      probability the two sit on neighbouring sites
  - energy, and its gap from the scattering continuum

The most-bound eigenstate per Delta is the new element. Its separation distribution (printed
as a sparkline over distance d=1..N-1) shows the pair clinging together.
"""
import sys

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
BLOCKS = "▁▂▃▄▅▆▇█"


def bond(N, b, P, Q):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == b else (Q if i == b + 1 else I2))
    return o


def H_xxz(N, delta):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        H += bond(N, b, X, X) + bond(N, b, Y, Y) + delta * bond(N, b, Z, Z)
    return H


def occupied_sites(idx, N):
    """Site positions (0 = leftmost / MSB) occupied in basis index idx."""
    return [s for s in range(N) if (idx >> (N - 1 - s)) & 1]


def sparkline(ys):
    ys = np.asarray(ys, float)
    hi = ys.max()
    if hi < 1e-12:
        return BLOCKS[0] * len(ys)
    idx = np.clip((ys / hi * 7).round().astype(int), 0, 7)
    return "".join(BLOCKS[i] for i in idx)


def analyse(N, delta):
    d = 2 ** N
    # Project onto the popcount-2 sector directly (H conserves number), so every eigenstate is
    # purely two-excitation. This avoids cross-sector degeneracy mixing at small Delta.
    idxs = [b for b in range(d) if bin(b).count("1") == 2]
    seps = [abs(occupied_sites(b, N)[0] - occupied_sites(b, N)[1]) for b in idxs]
    H = H_xxz(N, delta).real  # XX+YY+ZZ is real-symmetric in the computational basis
    Hsub = H[np.ix_(idxs, idxs)]
    w, V = np.linalg.eigh(Hsub)

    rows = []  # (energy, mean_sep, p_adjacent, sep_distribution)
    for k in range(len(w)):
        p = np.abs(V[:, k]) ** 2
        sep_dist = np.zeros(N, float)  # index by distance 1..N-1
        for m, s in enumerate(seps):
            sep_dist[s] += p[m]
        mean_sep = sum(dd * sep_dist[dd] for dd in range(1, N))
        rows.append((w[k], mean_sep, sep_dist[1], sep_dist))
    rows.sort(key=lambda r: r[1])  # most bound (smallest mean separation) first
    return rows


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def bound_pair_vec(N, delta):
    """The most-bound popcount-2 eigenstate of H_xxz(N, delta), embedded in full 2^N space."""
    d = 2 ** N
    idxs = [b for b in range(d) if bin(b).count("1") == 2]
    seps = [abs(occupied_sites(b, N)[0] - occupied_sites(b, N)[1]) for b in idxs]
    Hsub = H_xxz(N, delta).real[np.ix_(idxs, idxs)]
    w, V = np.linalg.eigh(Hsub)
    best_k, best = 0, 1e9
    for k in range(len(w)):
        p = np.abs(V[:, k]) ** 2
        msep = sum(seps[m] * p[m] for m in range(len(idxs)))
        if msep < best:
            best, best_k = msep, k
    full = np.zeros(d, complex)
    for m, idx in enumerate(idxs):
        full[idx] = V[m, best_k]
    return full


def P_adjacent(N):
    """Diagonal projector onto two-excitation states whose excitations sit on adjacent sites."""
    d = 2 ** N
    diag = np.zeros(d, complex)
    for b in range(d):
        occ = occupied_sites(b, N)
        if len(occ) == 2 and abs(occ[0] - occ[1]) == 1:
            diag[b] = 1.0
    return np.diag(diag)


def liouvillian(N, H, gamma):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += gamma * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def alive(N, delta, gammas, ts):
    """The convergence: prepare the bound pair (the new element born from interaction), let it
    live in the loop under the carrier gamma_0, and watch how long the carrier keeps it bound.
    Q = J/gamma with J = 1 (unit XY coupling), so Q = 1/gamma."""
    d = 2 ** N
    Padj = P_adjacent(N)
    floor = 2.0 / N  # equipartition adjacency: (N-1) adjacent pairs / C(N,2) = 2/N
    v0 = bound_pair_vec(N, delta)
    rho0 = np.outer(v0, v0.conj())
    vec0 = rho0.flatten(order="F")
    H = H_xxz(N, delta)
    print(f"Bound pair ALIVE in the loop, N={N}, Delta={delta}, J=1 (so Q=1/gamma)")
    print(f"  prepare the bound pair, evolve open-system; P(adjacent) starts ~1, relaxes to floor 2/N={floor:.3f}")
    print(f"  the carrier gamma_0 gives it time AND erodes it: gentler carrier (smaller gamma, higher Q) = longer life\n")
    print(f"  {'gamma':>6}  {'Q=1/g':>6}  {'P_adj(0)':>8}  {'P_adj(end)':>10}  {'half-life tau':>13}  P(adjacent)(t)")
    for g in gammas:
        if g == 0.0:
            # closed: the bound pair is an eigenstate, P_adj is frozen (eternally bound, no decay)
            w, V = np.linalg.eigh(H.real)
            c = V.T @ v0
            series = []
            for t in ts:
                psit = V @ (np.exp(-1j * w * t) * c)
                series.append(float((psit.conj() @ (Padj @ psit)).real))
            p0, pend = series[0], series[-1]
            hl = "frozen"
        else:
            L = liouvillian(N, H, g)
            wL, VL = np.linalg.eig(L)
            cL = np.linalg.solve(VL, vec0)
            series = []
            for t in ts:
                rho = (VL @ (np.exp(wL * t) * cL)).reshape(d, d, order="F")
                series.append(float(np.trace(Padj @ rho).real))
            p0, pend = series[0], series[-1]
            half = p0 - 0.5 * (p0 - floor)
            hl = next((f"{ts[i]:.2f}" for i in range(len(ts)) if series[i] <= half), ">tmax")
        spark = sparkline(np.array(series) - floor)  # relative to the floor
        print(f"  {g:>6.2f}  {1.0/g if g>0 else float('inf'):>6.1f}  {p0:>8.3f}  {pend:>10.3f}  {hl:>13}  {spark}")
    print(f"\n  reading: at gamma=0 the pair is frozen-bound (an eigenstate, eternal but timeless);")
    print(f"  the carrier brings time and slowly unbinds it; the half-life grows with Q=1/gamma.")
    print(f"  birth (interaction) + carrier (gamma_0) + the loop, in one frame.")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "alive":
        N = int(sys.argv[2]) if len(sys.argv) > 2 else 5
        ts = np.linspace(0.0, 60.0, 240)
        alive(N, 2.0, [0.0, 0.05, 0.2, 0.5, 1.0], ts)
        return
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    deltas = [0.0, 0.5, 1.0, 2.0, 4.0, -1.0, -4.0]
    print(f"Two-magnon bound state, N={N} chain, popcount-2 sector ({N*(N-1)//2} states)")
    print(f"  mean separation <|i-j|>: small = bound pair, large = spread/scattering")
    print(f"  a free chain (Delta=0) has no bound state; watch one split off as |Delta| grows\n")
    print(f"  {'Delta':>6}  {'min<sep>':>8}  {'P(adj)':>7}  {'E_bound':>9}  {'E_cont_lo':>9}  {'gap':>7}  sep-distribution d=1..{N-1}")
    for delta in deltas:
        rows = analyse(N, delta)
        if not rows:
            print(f"  {delta:>6.2f}   (no clean 2-excitation eigenstates)")
            continue
        e_bound, msep, padj, sep_dist = rows[0]          # most bound
        energies = sorted(r[0] for r in rows)
        # continuum = the bulk of 2-excitation levels; bound state sits at an extreme.
        cont_lo = energies[1] if len(energies) > 1 else e_bound
        cont_hi = energies[-2] if len(energies) > 1 else e_bound
        gap = min(abs(e_bound - cont_lo), abs(e_bound - cont_hi))
        spark = sparkline(sep_dist[1:N])
        tag = "BOUND" if msep < 1.6 and padj > 0.5 else ("partial" if msep < (N / 3) else "spread/free")
        print(f"  {delta:>6.2f}  {msep:>8.3f}  {padj:>7.3f}  {e_bound:>9.3f}  {cont_lo:>9.3f}  {gap:>7.3f}  {spark}  [{tag}]")
    print(f"\n  reading: at Delta=0 the tightest state is still spread (free, no binding); as |Delta|")
    print(f"  grows the most-bound state collapses to <sep>~1, P(adj)~1, and splits off in energy")
    print(f"  (the gap) -- a new composite object. Note the Delta -> -Delta symmetry of <sep>/P(adj).")


if __name__ == "__main__":
    main()
