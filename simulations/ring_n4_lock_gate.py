"""Gate for the ring N=4 lock and its scope.

Runs the checks behind `docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md` and the typed
claim `compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs`. Every
number those two state about the 4-cycle is produced here.

The law that survives:

    Delta_max(H_C4) = 3*J = (3/4)*J*N            (Casimir on K_{2,2}, exact)
    max|Im(lambda_L)| = Delta_max(H)             (uniform gamma; EVERY graph)

What is the 4-cycle's own is the CLOSED FORM. The value 3*J is NOT the ring's
alone: it is the MAXIMUM of Delta_max over connected graphs at N=4, and ten of
the 38 connected labelled graphs attain it (three labellings of C_4, six of the
diamond K_4 minus an edge, and K_4 itself). K_4 is not bipartite, so
bipartite-completeness cannot be what produces the value.

The realiser is |0000><singlet|, not a generic product of the two extremal
eigenspaces: the E = +J quintuplet has five members and only the two Z-product
ones (|0000> and |1111>) give a Liouvillian eigenoperator. For the other three
the residual is order 1.

Sign convention: H = J * sum_bonds S_i.S_j with J > 0, so |0000> is the MAXIMUM
energy eigenstate (E = +J = J*bonds/4) and the (S_A=1, S_B=1, S_tot=0) singlet
is the minimum (E = -2J). Hence lambda = -4*gamma - 3i*J for |0000><singlet| and
its conjugate for the reversed ordering; both are in the spectrum.

Usage:  python simulations/ring_n4_lock_gate.py
Importable: all checks live inside main().
"""
from __future__ import annotations

import collections
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

from star_saturation_gate import (
    Z,
    chain,
    complete,
    connected,
    dephasing_liouvillian,
    heisenberg,
    max_imag,
    op,
    ring,
    spread,
    star,
)

REPO = Path(__file__).resolve().parent.parent

CHECKS: list[tuple[str, bool]] = []

# The observed Im/sigma values hard-coded in
# compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs (_observedImOverSigma),
# read off the stored sweep as max(MaxImag, |MinImag|) / sigma. Section 8 rebuilds them
# from the JSON and compares bit for bit, so the typed surface cannot drift from the data.
CSHARP_OBSERVED = {
    "0.5000": 0.37500000000000044,
    "1.0000": 0.7499999999999999,
    "1.5000": 1.1250000000000009,
    "1.7321": 1.2990381056766593,
    "2.0000": 1.5000000000000049,
    "2.5000": 1.8750000000000095,
}


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, bool(ok)))
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")


def levels(H, nd=9):
    vals, cnts = np.unique(np.round(np.linalg.eigvalsh(H), nd), return_counts=True)
    return list(zip(vals.tolist(), cnts.tolist()))


def cluster(values, tol=1e-6):
    """Distinct values under a tolerance, so that a spread-out numerical zero
    counts once. Plain rounding splits it into several near-zero entries."""
    out = []
    for v in sorted(np.asarray(values).ravel().tolist()):
        if not out or abs(v - out[-1]) > tol:
            out.append(v)
    return np.array(out)


def eigop_residual(L, M):
    """Relative residual of the vectorised operator M as an eigenvector of L."""
    v = M.reshape(-1)
    Lv = L @ v
    lam = (v.conj() @ Lv) / (v.conj() @ v)
    return float(np.linalg.norm(Lv - lam * v) / np.linalg.norm(v)), complex(lam)


def main() -> int:
    N = 4
    C4 = ring(N)

    # --------------------------------------------- 1. the Casimir spectrum
    print("\n1. C_4 = K_{2,2} Casimir spectrum {-2J, -J^3, 0^7, +J^5}")
    H = heisenberg(N, C4, 1.0)
    lv = levels(H)
    check("four distinct levels", len(lv) == 4, f"{lv}")
    check("levels are -2J, -J, 0, +J",
          [v for v, _ in lv] == [-2.0, -1.0, 0.0, 1.0], f"{[v for v, _ in lv]}")
    check("multiplicities are 1, 3, 7, 5",
          [c for _, c in lv] == [1, 3, 7, 5], f"{[c for _, c in lv]}")
    check("multiplicities sum to 16 = 2^4", sum(c for _, c in lv) == 16)
    s = spread(H)
    check("Delta_max = 3J = (3/4)*J*N", abs(s - 3.0) < 1e-12,
          f"{s:.12f}; (3/4)*J*N = {0.75 * 1.0 * N}")

    print("\n1b. the |0000> end is E_max, the singlet end is E_min")
    check("E_max = +J = J*bonds/4", abs(np.linalg.eigvalsh(H).max() - len(C4) / 4) < 1e-12)
    check("|0000> IS that maximum", abs(H[0, 0].real - 1.0) < 1e-12, f"E={H[0,0].real:+.6f}")
    check("E_min = -2J is non-degenerate", levels(H)[0][1] == 1)

    # ------------------------------------------------ 2. the lock, uniform gamma
    print("\n2. max|Im(lambda_L)| = 3J, over J and over gamma")
    for J, g in ((1.0, 0.05), (1.0, 0.5), (0.1, 0.05), (2.0, 0.5)):
        im = max_imag(heisenberg(N, C4, J), N, [g] * N)
        check(f"J={J} gamma={g}: max|Im| = (3/4)*J*N", abs(im - 0.75 * J * N) < 1e-11,
              f"max|Im|={im:.12f}  Im/sigma={im/(N*g):.6f}  (3/4)Q={0.75*J/g:.6f}")

    print("\n2b. gamma-independence at fixed J, four decades")
    for g in (0.005, 0.05, 0.5, 5.0, 50.0):
        im = max_imag(H, N, [g] * N)
        check(f"gamma={g}: still 3J", abs(im - 3.0) < 1e-9, f"max|Im|={im:.12f}")

    print("\n2c. the six Q-sweep anchors of the proof table (gamma_0 = 0.05)")
    g0 = 0.05
    for label, Q in (("0.5", 0.5), ("1.0", 1.0), ("1.5", 1.5), ("sqrt3", np.sqrt(3.0)),
                     ("2.0", 2.0), ("2.5", 2.5)):
        J = Q * g0
        im = max_imag(heisenberg(N, C4, J), N, [g0] * N)
        obs = im / (N * g0)
        check(f"Q={label}: Im/sigma = (3/4)*Q", abs(obs - 0.75 * Q) < 1e-11,
              f"observed {obs:.12f} vs predicted {0.75*Q:.12f}  rel.err {abs(obs-0.75*Q)/(0.75*Q):.1e}")

    # ---------------------------------------------------- 3. the realising mode
    print("\n3. The realiser is |0000><singlet|, lambda = -4*gamma - 3i*J")
    g = 0.5
    L = dephasing_liouvillian(H, N, [g] * N)
    ev, evec = np.linalg.eigh(H)
    ferro = np.zeros(2 ** N, dtype=complex)
    ferro[0] = 1.0
    sing = evec[:, int(np.argmin(ev))]
    wts = sorted({bin(k).count("1") for k in range(2 ** N) if abs(sing[k]) > 1e-10})
    check("the E=-2J singlet has definite Hamming weight 2", wts == [2], f"weights {wts}")
    r, lam = eigop_residual(L, np.outer(ferro, sing.conj()))
    check("|0000><singlet| IS an L-eigenoperator", r < 1e-12, f"residual {r:.2e}")
    check("its eigenvalue is -4*gamma - 3i*J",
          abs(lam - (-4 * g - 3j)) < 1e-10, f"lambda = {lam:.12f}")
    check("the dephasing rate is the scalar -2*gamma*popcount = -4*gamma",
          abs(lam.real - (-2 * g * 2)) < 1e-10, f"Re = {lam.real:.12f}")
    for gg in (0.005, 0.05, 0.5, 5.0, 50.0):
        rg, lg = eigop_residual(dephasing_liouvillian(H, N, [gg] * N), np.outer(ferro, sing.conj()))
        check(f"still an eigenoperator at gamma={gg}, lambda = -4*gamma - 3i*J",
              rg < 1e-12 and abs(lg - (-4 * gg - 3j)) < 1e-10,
              f"residual {rg:.2e}  lambda = {lg:.10f}")
    r2, lam2 = eigop_residual(L, np.outer(sing, ferro.conj()))
    check("the reversed ordering carries the CONJUGATE eigenvalue",
          r2 < 1e-12 and abs(lam2 - (-4 * g + 3j)) < 1e-10, f"lambda = {lam2:.12f}")

    print("\n3b. 'a member of the S_tot=2 multiplet' is NOT enough")
    top = [k for k in range(2 ** N) if abs(ev[k] - 1.0) < 1e-9]
    check("the E=+J quintuplet has five members", len(top) == 5, f"{len(top)}")
    good, bad, rayleigh = [], [], []
    for k in top:
        rk, lk = eigop_residual(L, np.outer(evec[:, k], sing.conj()))
        (good if rk < 1e-12 else bad).append((k, rk))
        rayleigh.append(lk)
    check("only two of the five give an eigenoperator", len(good) == 2 and len(bad) == 3,
          f"eigen: {len(good)}, not: {[f'{r:.2f}' for _, r in bad]}")
    check("the two are the Z-product members |0000> and |1111>",
          len(good) == 2 and all(sorted({bin(b).count("1") for b in range(2 ** N)
                                         if abs(evec[b, k]) > 1e-10}) in ([0], [4])
                                 for k, _ in good),
          f"{len(good)} eigenoperator members")
    check("all five share the Rayleigh quotient -4*gamma - 3i*J, so <v,Lv> cannot see the difference",
          all(abs(lk - (-4 * g - 3j)) < 1e-10 for lk in rayleigh),
          f"max deviation {max(abs(lk - (-4*g - 3j)) for lk in rayleigh):.2e}")
    check("the failing residuals are the exact values 1, 2/sqrt(3), 1, not numerical noise",
          sorted(round(r, 12) for _, r in bad) == sorted(
              [1.0, 1.0, round(2 / np.sqrt(3.0), 12)]),
          f"{[repr(r) for _, r in bad]}  vs 2/sqrt(3) = {2/np.sqrt(3.0)!r}")

    # -------------------------------------------- 4. what is NOT the ring's
    print("\n4. 3J is the MAXIMUM over connected N=4 graphs, and it is not unique")
    edges = list(combinations(range(N), 2))
    graphs = []
    for mask in range(1, 1 << len(edges)):
        bonds = [edges[k] for k in range(len(edges)) if mask >> k & 1]
        if connected(N, bonds):
            graphs.append((spread(heisenberg(N, bonds)), tuple(bonds)))
    check("38 connected labelled graphs on 4 vertices", len(graphs) == 38, f"{len(graphs)}")
    mx = max(d for d, _ in graphs)
    at_max = [b for d, b in graphs if abs(d - mx) < 1e-9]
    mn = min(d for d, _ in graphs)
    at_min = [b for d, b in graphs if abs(d - mn) < 1e-9]
    check("the maximum spread is 3J", abs(mx - 3.0) < 1e-9, f"{mx:.12f}")
    check("the 4-cycle attains it", abs(spread(heisenberg(N, C4)) - mx) < 1e-9)
    check("but it is NOT the unique maximiser: ten labelled graphs tie",
          len(at_max) == 10, f"{len(at_max)} graphs at 3J")
    check("K_4 ties, and K_4 is not bipartite",
          abs(spread(heisenberg(N, complete(N))) - 3.0) < 1e-9,
          f"spread(K_4) = {spread(heisenberg(N, complete(N))):.12f}")
    k4lv = levels(heisenberg(N, complete(N)))
    check("K_4 reaches 3J with a DIFFERENT multiplet structure {-1.5^2, -0.5^9, +1.5^5}",
          [v for v, _ in k4lv] == [-1.5, -0.5, 1.5] and [c for _, c in k4lv] == [2, 9, 5],
          f"{k4lv}")
    check("the diamond K_4 minus an edge ties too",
          abs(spread(heisenberg(N, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)])) - 3.0) < 1e-9)
    degseq = collections.Counter(
        tuple(sorted(collections.Counter(v for e in b for v in e).values())) for b in at_max)
    check("the ten split 3 x C_4 + 6 x diamond + 1 x K_4, by degree sequence",
          degseq == {(2, 2, 2, 2): 3, (2, 2, 3, 3): 6, (3, 3, 3, 3): 1}, f"{dict(degseq)}")
    check("the minimum is the star's 2J = J*N/2", abs(mn - 2.0) < 1e-9, f"{mn:.12f}")
    check("four labelled stars attain the minimum", len(at_min) == 4, f"{len(at_min)}")
    check("the chain sits strictly between", 2.0 < spread(heisenberg(N, chain(N))) < 3.0,
          f"spread(P_4) = {spread(heisenberg(N, chain(N))):.12f}")

    print("\n4a. and the maximality does NOT extend past N=4")
    for Nb in (5, 6):
        rs, ks = spread(heisenberg(Nb, ring(Nb))), spread(heisenberg(Nb, complete(Nb)))
        check(f"N={Nb}: the ring is strictly below the complete graph", rs < ks - 1e-9,
              f"ring={rs:.10f}  K_{Nb}={ks:.10f}")

    print("\n4c. the complete graph is a THIRD elementary Casimir family")
    for Nb in range(3, 9):
        sk = spread(heisenberg(Nb, complete(Nb)))
        closed = Nb * (Nb + 2) / 8 if Nb % 2 == 0 else (Nb - 1) * (Nb + 3) / 8
        check(f"Delta_max(K_{Nb}) = " + ("N(N+2)/8" if Nb % 2 == 0 else "(N-1)(N+3)/8"),
              abs(sk - closed) < 1e-9, f"{sk:.10f} vs {closed:.10f}")
    print("\n4b. every N=4 topology saturates its OWN spread (the value differs, the law does not)")
    for name, bonds in (("ring C_4", C4), ("chain P_4", chain(N)),
                        ("star S_4", star(N)), ("complete K_4", complete(N))):
        Hb = heisenberg(N, bonds)
        d, im = spread(Hb), max_imag(Hb, N, [0.5] * N)
        check(f"{name}: max|Im| = spread(H)", abs(im - d) < 1e-11,
              f"spread={d:.10f}  max|Im|={im:.10f}  |diff|={abs(im-d):.1e}")

    # -------------------------------------------------------- 5. the fences
    print("\n5. Uniform gamma is load-bearing for the EQUALITY (the bound survives)")
    for prof in ([1.0, 0.5, 0.5, 0.5], [1.0, 0.1, 1.0, 0.1], [2.0, 0.0, 0.0, 0.0]):
        im = max_imag(H, N, prof)
        check(f"gamma={prof}: strictly below 3J, and still under the bound",
              im < 3.0 - 1e-6, f"max|Im|={im:.10f}  deficit {3.0-im:.3e}")

    print("\n5a. The FULL single-site jump set is load-bearing, not just Hermiticity")
    d2 = 2 ** N
    IdN = np.eye(d2, dtype=complex)

    def max_im_jumps(Hb, jumps):
        Lj = -1j * (np.kron(Hb, IdN) - np.kron(IdN, Hb.T))
        for c_op, rate in jumps:
            Lj += rate * (np.kron(c_op, c_op.conj())
                          - 0.5 * np.kron(c_op.conj().T @ c_op, IdN)
                          - 0.5 * np.kron(IdN, (c_op.conj().T @ c_op).T))
        return float(np.max(np.abs(np.linalg.eigvals(Lj).imag)))

    z0 = op(Z, 0, N)
    check("Z_0 alone at rate 0.5 (Hermitian) does NOT saturate",
          abs(max_im_jumps(H, [(z0, 0.5)]) - 2.8961214544) < 1e-8,
          f"max|Im| = {max_im_jumps(H, [(z0, 0.5)]):.10f} vs 3J")
    z0z1 = z0 @ op(Z, 1, N)
    check("the two-body Z_0Z_1 at rate 0.5 (Hermitian) does NOT saturate either",
          abs(max_im_jumps(H, [(z0z1, 0.5)]) - 2.9433995088) < 1e-8,
          f"max|Im| = {max_im_jumps(H, [(z0z1, 0.5)]):.10f} vs 3J")
    check("and those values are RATE-dependent, unlike the saturating case",
          abs(max_im_jumps(H, [(z0, 0.05)]) - 2.9990042309) < 1e-8,
          f"Z_0 alone at rate 0.05: {max_im_jumps(H, [(z0, 0.05)]):.10f}")
    check("gamma=(2,0,0,0) IS Z_0 alone at rate 2, the same system read the other way",
          abs(max_imag(H, N, [2.0, 0.0, 0.0, 0.0]) - max_im_jumps(H, [(z0, 2.0)])) < 1e-10,
          f"{max_imag(H, N, [2.0,0,0,0]):.10f} = {max_im_jumps(H, [(z0, 2.0)]):.10f}")
    check("the full set at one rate DOES", abs(max_im_jumps(H, [(op(Z, l, N), 0.5) for l in range(N)]) - 3.0) < 1e-10)

    print("\n5b. A fully polarised H extreme is load-bearing")
    xy = np.zeros((d2, d2), dtype=complex)
    for (i, j) in C4:
        xy += 0.25 * (op(np.array([[0, 1], [1, 0]], dtype=complex), i, N)
                      @ op(np.array([[0, 1], [1, 0]], dtype=complex), j, N)
                      + op(np.array([[0, -1j], [1j, 0]], dtype=complex), i, N)
                      @ op(np.array([[0, -1j], [1j, 0]], dtype=complex), j, N))
    im_xy = max_imag(xy, N, [0.5] * N)
    check("XY on the same ring does NOT saturate its spread",
          im_xy < spread(xy) - 1e-6, f"max|Im|={im_xy:.10f} vs spread={spread(xy):.10f}")

    print("\n5b'. Isotropy is NOT load-bearing; the polarised extreme is (F148 past Heisenberg)")

    def xxz(Nb, bonds, J=1.0, delta=1.0):
        dd = 2 ** Nb
        M = np.zeros((dd, dd), dtype=complex)
        XM = np.array([[0, 1], [1, 0]], dtype=complex)
        YM = np.array([[0, -1j], [1j, 0]], dtype=complex)
        for (i, j) in bonds:
            M += (J / 4) * (op(XM, i, Nb) @ op(XM, j, Nb) + op(YM, i, Nb) @ op(YM, j, Nb))
            M += (J / 4) * delta * (op(Z, i, Nb) @ op(Z, j, Nb))
        return M

    saturating = {}
    for label, delta, expect in (("XXZ delta=-2", -2.0, True), ("XXZ delta=-1.5", -1.5, True),
                                 ("XXZ delta=-1", -1.0, True), ("XXZ delta=-0.5", -0.5, False),
                                 ("XY (delta=0)", 0.0, False), ("XXZ delta=0.5", 0.5, False),
                                 ("XXZ delta=1 (isotropic)", 1.0, True), ("XXZ delta=1.5", 1.5, True),
                                 ("XXZ delta=2", 2.0, True)):
        Hd = xxz(N, C4, delta=delta)
        sd = spread(Hd)
        polarised = (abs(Hd[0, 0].real - np.linalg.eigvalsh(Hd).max()) < 1e-12
                     or abs(Hd[0, 0].real - np.linalg.eigvalsh(Hd).min()) < 1e-12)
        mis = [max_imag(Hd, N, [gg] * N) for gg in (0.05, 0.5, 2.0)]
        sats = [abs(mi - sd) < 1e-10 for mi in mis]
        check(f"{label}: |0...0> extremal <=> saturates at every gamma",
              polarised == expect and all(s == expect for s in sats),
              f"extremal={polarised}  spread={sd:.10f}  "
              f"max|Im| at gamma=(0.05,0.5,2) = " + " ".join(f"{m:.10f}" for m in mis))
        saturating[delta] = all(sats)
    check("so for XXZ the dividing line is exactly |delta| >= 1",
          {dd for dd, ok in saturating.items() if ok} == {dd for dd in saturating if abs(dd) >= 1.0},
          f"saturating deltas: {sorted(dd for dd, ok in saturating.items() if ok)}")

    tfi = np.zeros((d2, d2), dtype=complex)
    for (i, j) in C4:
        tfi += 0.25 * (op(Z, i, N) @ op(Z, j, N))
    for l in range(N):
        tfi += op(np.array([[0, 1], [1, 0]], dtype=complex), l, N)
    hxxz = xxz(N, C4, delta=0.5)
    check("the two failing XXZ delta=0.5 values the docs quote",
          abs(max_imag(hxxz, N, [0.05] * N) - 2.8672194968) < 5e-10
          and abs(max_imag(hxxz, N, [2.0] * N) - 2.1861406616) < 5e-10,
          f"gamma=0.05: {max_imag(hxxz, N, [0.05]*N):.10f}   gamma=2: {max_imag(hxxz, N, [2.0]*N):.10f}")
    hising = np.zeros((d2, d2), dtype=complex)
    for (i, j) in C4:
        hising += 0.25 * (op(Z, i, N) @ op(Z, j, N))
    check("pure Ising: spread 2.0000000000, saturating at every gamma",
          abs(spread(hising) - 2.0) < 1e-12
          and all(abs(max_imag(hising, N, [gg] * N) - 2.0) < 1e-10 for gg in (0.05, 0.5, 2.0)),
          f"spread={spread(hising):.10f}")
    check("XXZ delta=2 spread is exactly 3 + sqrt(3)",
          abs(spread(xxz(N, C4, delta=2.0)) - (3 + np.sqrt(3.0))) < 1e-12,
          f"{spread(xxz(N, C4, delta=2.0)):.12f} vs {3+np.sqrt(3.0):.12f}")
    check("transverse-field Ising falls SHORT of its spread (not over it)",
          max_imag(tfi, N, [0.5] * N) < spread(tfi) - 1e-3,
          f"max|Im|={max_imag(tfi, N, [0.5]*N):.7f} vs spread={spread(tfi):.7f}")

    print("\n5c. The 'adds only real decay' justification is false; Hermitian jumps carry it")
    c = np.eye(2, dtype=complex) + 1j * np.array([[0, -1j], [1j, 0]], dtype=complex)
    D = (np.kron(c, c.conj())
         - 0.5 * np.kron(c.conj().T @ c, np.eye(2))
         - 0.5 * np.kron(np.eye(2), (c.conj().T @ c).T))
    w = np.linalg.eigvals(D)
    want = sorted([0j, 0j, -2 + 2j, -2 - 2j], key=lambda z: (z.real, z.imag))
    got = sorted(w.tolist(), key=lambda z: (round(z.real, 9), round(z.imag, 9)))
    check("H=0 with c=I+iY has spectrum {0, 0, -2+2i, -2-2i}: max|Im| = 2 against a bound of 0",
          all(abs(a - b) < 1e-9 for a, b in zip(want, got)) and abs(max(abs(w.imag)) - 2.0) < 1e-9,
          f"spectrum {np.round(np.sort_complex(w), 6)}")
    check("what it lacks is self-adjointness: ||D - D*||_F = 4*sqrt(2)",
          abs(np.linalg.norm(D - D.conj().T) - 4 * np.sqrt(2.0)) < 1e-12,
          f"||D - D*||_F = {np.linalg.norm(D - D.conj().T):.4f} (spectral norm "
          f"{np.linalg.norm(D - D.conj().T, 2):.4f})")

    # -------------------------------------------- 6. the numerical-range reading
    print("\n6. The bound is an interval, not a fixed point")
    w = np.linalg.eigvals(L)
    ims = cluster(w.imag)
    diffs = cluster(np.subtract.outer(np.linalg.eigvalsh(H), np.linalg.eigvalsh(H)).ravel())
    outside = [x for x in ims if min(abs(diffs - x)) > 1e-6]
    check("35 distinct imaginary parts at gamma=0.5", len(ims) == 35, f"{len(ims)}")
    check("the signed difference set has 7 members", len(diffs) == 7, f"{sorted(diffs)}")
    check("28 of the 35 are NOT in the difference set", len(outside) == 28, f"{len(outside)}")
    check("none of them leaves the interval", max(abs(ims)) <= 3.0 + 1e-9,
          f"max|Im| = {max(abs(ims)):.12f}")

    # ------------------------------------------- 7. the larger rings, c_N
    print("\n7. c_N = 1/4 - E_0/(J*N): the closed forms and the ln 2 limit")
    c4 = 0.25 - min(np.linalg.eigvalsh(H)) / N
    check("c_4 = 3/4", abs(c4 - 0.75) < 1e-12, f"{c4:.12f}")
    H6 = heisenberg(6, ring(6))
    c6 = 0.25 - min(np.linalg.eigvalsh(H6)) / 6
    check("c_6 = (5 + sqrt(13))/12", abs(c6 - (5 + np.sqrt(13)) / 12) < 1e-11,
          f"{c6:.12f} vs {(5+np.sqrt(13))/12:.12f}")
    check("c_6 reproduces the empirical 0.717129", abs(c6 - 0.717129) < 1e-6, f"{c6:.6f}")
    e0_6 = min(np.linalg.eigvalsh(H6))
    check("4*E_0(N=6) is a root of x^2 + 8x - 36",
          abs((4 * e0_6) ** 2 + 8 * (4 * e0_6) - 36) < 1e-9, f"4*E_0 = {4*e0_6:.12f}")
    c8 = 0.25 - min(np.linalg.eigvalsh(heisenberg(8, ring(8)))) / 8
    check("c_8 is a root of 512c^3 - 640c^2 + 232c - 25",
          abs(512 * c8 ** 3 - 640 * c8 ** 2 + 232 * c8 - 25) < 1e-9, f"c_8 = {c8:.12f}")
    check("the EVEN sequence decreases and crosses 1/sqrt(2) between N=6 and N=8",
          c4 > c6 > 1 / np.sqrt(2) > c8 > np.log(2),
          f"{c4:.6f} > {c6:.6f} > {1/np.sqrt(2):.6f} > {c8:.6f} > {np.log(2):.6f}")

    print("\n7a. and the ODD rings approach ln 2 from BELOW, rising: two branches")
    odd = {Nb: 0.25 - min(np.linalg.eigvalsh(heisenberg(Nb, ring(Nb)))) / Nb for Nb in (5, 7, 9, 11)}
    for Nb, expect in ((5, 0.623607), (7, 0.657883), (9, 0.671922), (11, 0.678994)):
        check(f"c_{Nb} = {expect} and sits BELOW ln 2",
              abs(odd[Nb] - expect) < 1e-6 and odd[Nb] < np.log(2), f"c_{Nb} = {odd[Nb]:.6f}")
    check("the odd branch rises", odd[5] < odd[7] < odd[9] < odd[11],
          " < ".join(f"{odd[k]:.6f}" for k in (5, 7, 9, 11)))
    check("so c_N over ALL N is not monotone", not (c4 > odd[5] > c6),
          f"c_4={c4:.6f}, c_5={odd[5]:.6f}, c_6={c6:.6f}")

    print("\n7b. saturation at N=6: the |ferro><ground| mode, without a 4096^2 eigensolve")
    for Nb in (6, 8):
        Hb = heisenberg(Nb, ring(Nb))
        evb, evcb = np.linalg.eigh(Hb)
        gs = evcb[:, int(np.argmin(evb))]
        w_gs = sorted({bin(k).count("1") for k in range(2 ** Nb) if abs(gs[k]) > 1e-9})
        check(f"N={Nb}: the ring ground state has definite Hamming weight {Nb//2}",
              w_gs == [Nb // 2], f"weights {w_gs}")
        check(f"N={Nb}: |ferro> is a Z-product state at E_max",
              abs(Hb[0, 0].real - evb.max()) < 1e-12)
    for Nb in (6, 8, 10):
        Hb = heisenberg(Nb, ring(Nb))
        evb, evcb = np.linalg.eigh(Hb)
        gs = evcb[:, int(np.argmin(evb))]
        M = np.zeros((2 ** Nb, 2 ** Nb), dtype=complex)
        M[0, :] = gs.conj()
        comm = Hb @ M - M @ Hb
        dE = evb.max() - evb.min()
        xors = {bin(j).count("1") for j in range(2 ** Nb) if abs(M[0, j]) > 1e-9}
        check(f"N={Nb}: [H, |ferro><gs|] = +Delta_max * M with a constant XOR weight",
              np.linalg.norm(comm - dE * M) < 1e-9 * dE and xors == {Nb // 2},
              f"||[H,M] - Delta*M|| = {np.linalg.norm(comm - dE * M):.2e}, XOR weights {sorted(xors)}"
              f" => D acts as the scalar -2*gamma*{Nb // 2}")
    Hr6 = heisenberg(6, ring(6))
    L6 = dephasing_liouvillian(Hr6, 6, [0.05] * 6)
    ev6, evc6 = np.linalg.eigh(Hr6)
    f6 = np.zeros(2 ** 6, dtype=complex)
    f6[0] = 1.0
    r6, lam6 = eigop_residual(L6, np.outer(f6, evc6[:, int(np.argmin(ev6))].conj()))
    check("N=6: residual of |ferro><gs| is machine zero", r6 < 1e-11, f"{r6:.2e}")
    check("N=6: |Im| = Delta_max = 4.302776", abs(abs(lam6.imag) - spread(Hr6)) < 1e-10,
          f"|Im|={abs(lam6.imag):.6f} spread={spread(Hr6):.6f}")
    check("N=6: Re = -2*gamma*3 = -0.3", abs(lam6.real + 0.3) < 1e-10, f"{lam6.real:.12f}")

    # -------------------------------------------------- 8. the stored anchors
    print("\n8. The stored Q-sweep JSON anchors")
    adir = REPO / "simulations" / "results" / "q_sweep_anchor"
    for label, Q in (("0.5000", 0.5), ("1.0000", 1.0), ("1.5000", 1.5),
                     ("1.7321", np.sqrt(3.0)), ("2.0000", 2.0), ("2.5000", 2.5)):
        p = adir / f"ring_N4_Q{label}.json"
        if not p.exists():
            check(f"anchor file {p.name} exists", False, "MISSING")
            continue
        d = json.loads(p.read_text())
        check(f"anchor file {p.name} exists", True)
        mi = d.get("MaxImag", d.get("max_imag"))
        g = d.get("GammaValue", 0.05)
        obs = mi / (4 * g)
        check(f"{p.name}: stored Im/sigma = (3/4)*Q", abs(obs - 0.75 * Q) / (0.75 * Q) < 1e-13,
              f"stored MaxImag={mi} -> Im/sigma={obs:.12f} vs {0.75*Q:.12f}")
        two = max(mi, abs(d["MinImag"])) / (4 * g)
        stored = CSHARP_OBSERVED[label]
        check(f"{p.name}: two-sided max|Im|/sigma reproduces RingN4DihedralLockClaim's stored value",
              two == stored and abs(two - 0.75 * Q) / (0.75 * Q) < 1e-14,
              f"two-sided {two!r} vs C# {stored!r}  rel.err {abs(two - 0.75*Q)/(0.75*Q):.2e}"
              f"  {'|MinImag| is the larger' if abs(d['MinImag']) > mi else 'MaxImag is the larger'}")

    # ------------------------------------------------------------ summary
    bad_checks = [n for n, ok in CHECKS if not ok]
    print(f"\n{'=' * 70}")
    print(f"{len(CHECKS) - len(bad_checks)}/{len(CHECKS)} checks passed")
    if bad_checks:
        print("FAILED:")
        for n in bad_checks:
            print(f"  - {n}")
    return 1 if bad_checks else 0


if __name__ == "__main__":
    sys.exit(main())
