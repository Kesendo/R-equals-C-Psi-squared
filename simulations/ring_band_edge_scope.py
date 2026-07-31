"""Gate for the SCOPE of the XY floor band-edge law, the part PROOF_RING_GAP_DOMINANCE adds.

`ring_gap_dominance.py` and `ring_gap_completeness.py` gate the ring numbers. This gates the
sentences around them: who else obeys the law, what the name does not name, and how thin the
N=4 exception is.

Convention, the proof's own: H = (J/2) sum_bonds (X_iX_j + Y_iY_j), uniform Z-dephasing at
rate gamma, L = -i[H,.] + sum_l gamma (Z_l rho Z_l - rho). "The floor" is the set of
Liouvillian eigenvalues with Re = -2*gamma exactly (tolerance 1e-8), and the quantity is
max|Im| over that set. It is NOT the spectral max|Im| (F148); at N=6 the two read 2.0 and
7.996 against DeltaE_max = 8.0.

What this gate pins:

  STAGE 1  the law is not the ring's: over ALL 38 connected labelled graphs at N=4, exactly
           three violate max|Im|_floor = J*rho, and they are the three labellings of the
           4-cycle. Named graphs at N=5 (star, complete, path, an asymmetric one) all obey it.
  STAGE 2  it is not even about equal couplings: random bond weights obey it against their own
           WEIGHTED adjacency radius, and the N=4 exception disappears when the weights differ.
  STAGE 3  the name names nothing: at even N the maximiser is NOT unique (m=0 and m=N/2 both
           attain |E| = 2), so no cyclic-symmetry argument singles out the uniform mode there.
  STAGE 4  the four scope fences: uniform gamma, the full single-site {Z_l} set, the XY form,
           and the irrelevance of sign(J).
  STAGE 5  the exception is fine-tuned in COUPLING space, not just in gamma: a wrap-bond
           detuning of 1e-4 removes it, while 1e-6 does not.
  STAGE 6  the N=4 floor is the numerically thin one (nearest off-floor eigenvalue 4.7e-5,
           against ~3e-2 at N=3, 5, 6).

  --slow adds the N=6 checks: the isolated 28-mode crossing at gamma = 0.5 exactly, and the
  floor-vs-spectrum comparison against DeltaE_max.

Usage:  python simulations/ring_band_edge_scope.py [--slow]
Importable: all checks live inside main().
"""
from __future__ import annotations

import itertools
import sys

import numpy as np
from scipy.linalg import eig

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

CHECKS: list[tuple[str, bool]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, bool(ok)))
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")


def op(P, k, N):
    M = np.array([[1]], dtype=complex)
    for s in range(N):
        M = np.kron(M, P if s == k else I2)
    return M


def xy(N, bonds, weights=None, delta=0.0):
    """H = (1/2) sum_bonds w*(XX + YY) + (delta/2) sum_bonds w*ZZ, w = 1 unless given."""
    w = [1.0] * len(bonds) if weights is None else list(weights)
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for (i, j), wij in zip(bonds, w):
        H += (wij / 2.0) * (op(X, i, N) @ op(X, j, N) + op(Y, i, N) @ op(Y, j, N))
        if delta:
            H += (wij / 2.0) * delta * (op(Z, i, N) @ op(Z, j, N))
    return H


def liou(N, H, gammas, sites=None):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    gs = [gammas] * N if np.isscalar(gammas) else list(gammas)
    for l in (range(N) if sites is None else sites):
        Zl = op(Z, l, N)
        L += gs[l] * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def floor_of(N, H, g, sites=None, gammas=None):
    """(count, max|Im|) over the eigenvalues at Re = -2*g exactly."""
    w = eig(liou(N, H, gammas if gammas is not None else g, sites), right=False)
    fl = w[np.abs(w.real + 2 * g) < 1e-8]
    return len(fl), (float(np.abs(fl.imag).max()) if len(fl) else float("nan")), w


def rho(N, bonds, weights=None):
    w = [1.0] * len(bonds) if weights is None else list(weights)
    A = np.zeros((N, N))
    for (i, j), wij in zip(bonds, w):
        A[i, j] = A[j, i] = wij
    return float(np.abs(np.linalg.eigvalsh(A)).max())


def connected(N, bonds):
    seen, st = {0}, [0]
    adj = {v: set() for v in range(N)}
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    while st:
        v = st.pop()
        for u in adj[v]:
            if u not in seen:
                seen.add(u)
                st.append(u)
    return len(seen) == N


def ring(N):
    return [(i, (i + 1) % N) for i in range(N)]


def degseq(N, bonds):
    d = [0] * N
    for i, j in bonds:
        d[i] += 1
        d[j] += 1
    return tuple(sorted(d))


def main() -> int:
    slow = "--slow" in sys.argv
    g = 0.05

    # -------------------------------------------------- 1. the law is not the ring's
    print("\n1. Over ALL connected graphs at N=4, who violates max|Im|_floor = J*rho?")
    edges = list(itertools.combinations(range(4), 2))
    tot, viol = 0, []
    for m in range(1, 1 << len(edges)):
        b = [edges[k] for k in range(len(edges)) if m >> k & 1]
        if not connected(4, b):
            continue
        tot += 1
        _, fm, _ = floor_of(4, xy(4, b), g)
        if abs(fm - rho(4, b)) > 1e-6:
            viol.append(tuple(b))
    check("38 connected labelled graphs on 4 vertices", tot == 38, f"{tot}")
    check("exactly three violate the law", len(viol) == 3, f"{len(viol)}")
    check("and all three are 4-cycles (degree sequence 2,2,2,2)",
          all(degseq(4, b) == (2, 2, 2, 2) for b in viol),
          f"{[degseq(4, b) for b in viol]}")

    print("\n1b. named graphs at N=5, including one with no symmetry")
    for nm, b in (("ring", ring(5)), ("path", [(i, i + 1) for i in range(4)]),
                  ("star", [(0, k) for k in range(1, 5)]),
                  ("complete", list(itertools.combinations(range(5), 2))),
                  ("bull", [(0, 1), (1, 2), (2, 0), (1, 3), (2, 4)]),
                  ("paw+tail", [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)])):
        cnt, fm, _ = floor_of(5, xy(5, b), g)
        r = rho(5, b)
        check(f"{nm}: floor max = J*rho", abs(fm - r) < 1e-7,
              f"max|Im| {fm:.7f}  rho {r:.7f}  (floor modes {cnt})")

    # ------------------------------------------- 2. not about equal couplings either
    print("\n2. Random bond weights: still J*rho, against the WEIGHTED radius")
    rng = np.random.default_rng(7)
    for nm, N, b in (("ring4", 4, ring(4)), ("ring5", 5, ring(5)),
                     ("path4", 4, [(0, 1), (1, 2), (2, 3)])):
        for trial in range(2):
            w = rng.uniform(0.4, 1.6, len(b))
            cnt, fm, _ = floor_of(N, xy(N, b, weights=w), g)
            r = rho(N, b, weights=w)
            check(f"{nm} trial{trial}: floor max = J*rho(weighted)", abs(fm - r) < 1e-7,
                  f"{fm:.9f} vs {r:.9f}  (floor modes {cnt})")
    cnt, _, _ = floor_of(4, xy(4, ring(4), weights=[1.0, 1.1, 0.9, 1.05]), g)
    check("and the N=4 exception is gone once the weights differ", cnt == 16,
          f"floor modes {cnt} (16 = the pure V_1 count)")

    # ------------------------------------------------ 3. the name names nothing
    print("\n3. At even N the band-top maximiser is NOT unique, so C_N singles out nothing")
    for N in (4, 5, 6, 7, 8):
        a = np.array([2 * np.cos(2 * np.pi * m / N) for m in range(N)])
        hits = int(np.sum(np.abs(np.abs(a) - np.abs(a).max()) < 1e-12))
        check(f"N={N}: |E| max attained {hits}x", hits == (2 if N % 2 == 0 else 1),
              f"eigenvalues {np.round(a, 6)}")

    # ----------------------------------------------------- 4. the four fences
    print("\n4. The scope fences")
    H4, H5 = xy(4, ring(4)), xy(5, ring(5))
    for N, H, prof in ((4, H4, [0.1, 0.05, 0.05, 0.05]), (5, H5, [0.1, 0.05, 0.05, 0.05, 0.05])):
        cnt, _, _ = floor_of(N, H, float(np.mean(prof)), gammas=prof)
        check(f"non-uniform gamma at N={N}: the exact floor set is EMPTY", cnt == 0,
              f"modes at Re=-2*mean(gamma): {cnt}")
    # N=4 is confounded here: its own exception already breaks the law, so the clean
    # demonstration is N=5 (and N=6 under --slow), where the law holds with the full set.
    cnt, fm, _ = floor_of(5, H5, g, sites=range(1, 5))
    check("dropping Z_0 at N=5 breaks the law: floor max is no longer J*rho",
          abs(fm - 1.6180340) < 1e-6 and abs(fm - rho(5, ring(5))) > 1e-6,
          f"floor {cnt}  max|Im| {fm:.7f}  vs rho {rho(5, ring(5)):.7f}")
    cnt, fm, w = floor_of(4, H4, g, sites=[0, 1])
    check("keeping only {Z_0,Z_1} makes the statement vacuous",
          cnt == 126 and abs(fm - np.abs(w.imag).max()) < 1e-9,
          f"floor {cnt} modes, floor max {fm:.7f} = whole-spectrum max {np.abs(w.imag).max():.7f}")
    for N, want in ((4, 2.0200000), (5, 1.9800000)):
        cnt, fm, _ = floor_of(N, xy(N, ring(N), delta=0.01), g)
        check(f"XXZ delta=0.01 at N={N} moves the value off J*rho", abs(fm - want) < 1e-6,
              f"floor {cnt}  max|Im| {fm:.7f}")
    for N, want in ((3, 2.0), (4, 2.8266588), (5, 2.0)):
        _, fm, _ = floor_of(N, xy(N, ring(N), weights=[-1.0] * N), g)
        check(f"sign(J) is irrelevant at N={N} (rho is a radius)", abs(fm - want) < 1e-6,
              f"max|Im| at J=-1: {fm:.7f}")

    # --------------------------------- 5. the exception is fine-tuned in coupling space
    print("\n5. The N=4 exception under a wrap-bond detuning J -> J(1-delta)")
    for d, want_ex in ((0.0, True), (1e-6, True), (1e-4, False), (1e-3, False)):
        cnt, fm, _ = floor_of(4, xy(4, ring(4), weights=[1, 1, 1, 1 - d]), g)
        check(f"delta={d:g}: exception {'survives' if want_ex else 'is gone'}",
              (fm > 2.0 + 1e-9) == want_ex,
              f"floor {cnt}  max|Im| {fm:.9f}  rho {rho(4, ring(4), [1,1,1,1-d]):.9f}")

    # -------------------------------------------- 6. how thin the N=4 floor really is
    print("\n6. Distance from the floor to the nearest NON-floor eigenvalue")
    seps = {}
    for N in (3, 4, 5):
        _, _, w = floor_of(N, xy(N, ring(N)), g)
        d = np.abs(w.real + 2 * g)
        seps[N] = float(np.sort(d[d > 1e-8])[0])
        check(f"N={N}: nearest off-floor Re-distance", True, f"{seps[N]:.3e}")
    check("N=4 is ~1000x thinner than its neighbours",
          seps[4] < seps[3] / 100 and seps[4] < seps[5] / 100,
          f"N=4 {seps[4]:.2e} vs N=3 {seps[3]:.2e}, N=5 {seps[5]:.2e}")

    # ------------------------------------------------------------ slow: N=6
    if slow:
        print("\n7. (--slow) N=6: the isolated crossing, and floor vs whole spectrum")
        H6 = xy(6, ring(6))
        dE = float(np.linalg.eigvalsh(H6).max() - np.linalg.eigvalsh(H6).min())
        counts = {}
        for gg in (0.49, 0.5, 0.51):
            cnt, fm, _ = floor_of(6, H6, gg)
            counts[gg] = cnt
            check(f"gamma={gg}: floor modes", True, f"{cnt}  max|Im| {fm:.7f}")
        check("the 28 at gamma=0.5 is an ISOLATED crossing, not a trend",
              counts[0.5] == 28 and counts[0.49] == 24 and counts[0.51] == 24, f"{counts}")
        cnt, fm, w = floor_of(6, H6, 0.05)
        check("floor max is far below the spectral max (F148 is a different quantity)",
              abs(fm - 2.0) < 1e-9 and np.abs(w.imag).max() > 7.9,
              f"floor {fm:.7f}  spectrum {np.abs(w.imag).max():.7f}  DeltaE_max {dE:.7f}")
    else:
        print("\n7. (N=6 checks need --slow: three 4096x4096 eigendecompositions.)")

    bad = [n for n, ok in CHECKS if not ok]
    print(f"\n{'=' * 70}")
    print(f"{len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed")
    for n in bad:
        print(f"  FAILED: {n}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
