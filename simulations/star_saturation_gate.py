"""Gate for the star Im-max law and its scope.

Runs the checks behind `experiments/STAR_CONFOCAL_LIMIT.md` and
`docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md`. Every number those two
documents state about the Liouvillian spectrum is produced here.

The law that survives:

    Delta_max(H_star) = J*N/2                       (Casimir, exact, every N)
    max|Im(lambda_L)| = Delta_max(H)                (EVERY graph, uniform gamma)

so the star's content is the CLOSED FORM and the fact that J*N/2 is the MINIMUM
over connected graphs (searched at N <= 6), not the saturation, which is universal.

The realising modes are |beta_k><ferro| coherences, with
lambda = -2*gamma*k + i*Delta_max: the fully polarised extreme is a Z-product
state, so uniform dephasing acts on that block as the SCALAR -2*gamma*k.
Three hypotheses are load-bearing and each has a fence below:
  - drop the fully polarised extreme (XY, transverse-field Ising): saturation fails
  - drop uniform gamma (site-dependent gamma_l): saturation fails, bound survives
  - drop Hermitian jump operators: the BOUND itself fails (H = 0, c = I + iY)

Sign convention: H = J * sum_bonds S_i.S_j with J > 0, so |0...0> is the
MAXIMUM-energy eigenstate (each bond simultaneously at its triplet maximum +J/4).

Usage:  python simulations/star_saturation_gate.py
Importable: all checks live inside main(), so `from star_saturation_gate import
heisenberg` does not run the suite.
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent

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


def heisenberg(N, bonds, J=1.0):
    """H = J * sum_bonds S_i.S_j = (J/4) * sum_bonds (XX + YY + ZZ)."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for P in (X, Y, Z):
            H += (J / 4.0) * (op(P, i, N) @ op(P, j, N))
    return H


def dephasing_liouvillian(H, N, gammas):
    """L = -i[H, .] + sum_l gamma_l (Z_l . Z_l - .)."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = op(Z, k, N)
        L += gammas[k] * (np.kron(Zk, Zk.T) - np.kron(Id, Id))
    return L


def star(N):
    return [(0, k) for k in range(1, N)]


def chain(N):
    return [(i, i + 1) for i in range(N - 1)]


def ring(N):
    return [(i, (i + 1) % N) for i in range(N)]


def complete(N):
    return list(combinations(range(N), 2))


def max_imag(H, N, gammas):
    return float(np.abs(np.linalg.eigvals(dephasing_liouvillian(H, N, gammas)).imag).max())


def spread(H):
    w = np.linalg.eigvalsh(H)
    return float(w.max() - w.min())


def connected(N, bonds):
    adj = {i: set() for i in range(N)}
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    seen, stack = {0}, [0]
    while stack:
        for nb in adj[stack.pop()]:
            if nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return len(seen) == N


def automorphism_count(N, bonds):
    E = {frozenset(b) for b in bonds}
    return sum(1 for p in permutations(range(N))
               if {frozenset((p[i], p[j])) for (i, j) in bonds} == E)


def emin_rungs(N, bonds):
    """Which Hamming-weight sectors contain the global H minimum."""
    H = heisenberg(N, bonds)
    d = 2 ** N
    Emin = np.linalg.eigvalsh(H).min()
    out = []
    for k in range(N + 1):
        idx = [i for i in range(d) if bin(i).count("1") == k]
        if abs(np.linalg.eigvalsh(H[np.ix_(idx, idx)]).min() - Emin) < 1e-10:
            out.append(k)
    return out


def main() -> int:
    # ------------------------------------------------------------ 1. Casimir
    print("\n1. Hub-leaf Casimir: Delta_max(H_star) = J*N/2")
    for N in range(3, 8):
        s = spread(heisenberg(N, star(N)))
        check(f"N={N} spread(H_star) = J*N/2", abs(s - N / 2) < 1e-12, f"{s:.12f} vs {N/2}")

    print("\n1b. The S_L = 0 sector carries ONE level, not two (odd N)")
    for N in (3, 5, 7):
        w = sorted({round(v, 9) for v in np.linalg.eigvalsh(heisenberg(N, star(N)))})
        check(f"N={N}: E = -0.5 absent from the star spectrum", -0.5 not in w, f"spectrum {w}")

    print("\n1c. |0...0> is the H MAXIMUM (the J > 0 sign convention)")
    for N in (3, 4, 5, 6):
        for nm, bonds in (("star", star(N)), ("chain", chain(N)), ("ring", ring(N))):
            H = heisenberg(N, bonds)
            check(f"N={N} {nm}: |0...0> at E_max = J*B/4",
                  abs(H[0, 0].real - np.linalg.eigvalsh(H).max()) < 1e-12
                  and abs(H[0, 0].real - len(bonds) / 4) < 1e-12,
                  f"E={H[0,0].real:+.6f}, B={len(bonds)}")

    # ----------------------------------------- 2. the star law, uniform gamma
    print("\n2. max|Im(lambda_L)| = J*N/2 for the star, over J and over gamma")
    for N in (3, 4, 5, 6):
        for J, g in ((1.0, 0.5), (1.0, 0.05), (0.1, 0.05), (0.125, 0.05)):
            im = max_imag(heisenberg(N, star(N), J), N, [g] * N)
            check(f"N={N} J={J} gamma={g}", abs(im - J * N / 2) < 1e-11,
                  f"max|Im|={im:.12f} J*N/2={J*N/2:.12f} Im/sigma={im/(N*g):.6f} Q/2={J/(2*g):.6f}")

    print("\n2b. gamma-independence at FIXED J, over four decades of gamma")
    for N in (3, 4, 5):
        J = 1.0
        for g in (0.005, 0.05, 0.5, 5.0, 50.0):
            im = max_imag(heisenberg(N, star(N), J), N, [g] * N)
            check(f"N={N} J=1 gamma={g}", abs(im - J * N / 2) < 1e-9,
                  f"max|Im|={im:.12f} (sigma = {N*g})")

    # ---------------------------- 3. the saturation is UNIVERSAL, not the star's
    print("\n3. Every Heisenberg topology saturates its OWN spread")
    asym6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (3, 5)]
    for N in (3, 4, 5, 6):
        tops = [("star", star(N)), ("chain", chain(N)), ("ring", ring(N)), ("complete", complete(N))]
        if N == 5:
            tops.append(("T-shape", [(0, 1), (1, 2), (2, 3), (1, 4)]))
        if N == 6:
            tops.append(("K3+K3", [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)]))
            tops.append(("asymmetric", asym6))
        for name, bonds in tops:
            H = heisenberg(N, bonds)
            s, im = spread(H), max_imag(H, N, [0.5] * N)
            check(f"N={N} {name}: max|Im| = spread(H)", abs(im - s) < 1e-11,
                  f"{im:.10f} vs {s:.10f}  |diff|={abs(im-s):.1e}  Im/sigma={im/(N*0.5):.6f}")
    check("the 'asymmetric' graph really has no symmetry",
          automorphism_count(6, asym6) == 1, f"|Aut| = {automorphism_count(6, asym6)}")
    check("(for contrast) the N=5 T-shape DOES have a symmetry",
          automorphism_count(5, [(0, 1), (1, 2), (2, 3), (1, 4)]) == 2,
          f"|Aut| = {automorphism_count(5, [(0,1),(1,2),(2,3),(1,4)])}")

    # -------------------------------------------------- 4. the realising modes
    print("\n4. The realisers are |beta_k><ferro|: lambda = -2*gamma*k + i*Delta")
    for N in (3, 4, 5, 6):
        J, g = 1.0, 0.5
        d = 2 ** N
        H = heisenberg(N, star(N), J)
        L = dephasing_liouvillian(H, N, [g] * N)
        w = np.linalg.eigvalsh(H)
        Emax, Emin = w.max(), w.min()
        ferro = np.zeros(d, dtype=complex)
        ferro[0] = 1.0
        exact = 0
        for k in range(1, N):
            idx = [i for i in range(d) if bin(i).count("1") == k]
            ws, Vs = np.linalg.eigh(H[np.ix_(idx, idx)])
            if abs(ws.min() - Emin) > 1e-10:
                continue
            beta = np.zeros(d, dtype=complex)
            beta[idx] = Vs[:, 0]
            v = np.kron(beta, np.conj(ferro))
            v /= np.linalg.norm(v)
            lam = np.vdot(v, L @ v)
            resid = float(np.linalg.norm(L @ v - lam * v))
            # sign matters: [H, |beta><ferro|] = -Delta*rho, so -i[H,.] gives +i*Delta
            if (resid < 1e-10 and abs(lam.real + 2 * g * k) < 1e-10
                    and abs(lam.imag - J * N / 2) < 1e-10):
                exact += 1
        check(f"N={N}: one exact realiser per Hamming rung k = 1..N-1", exact == N - 1,
              f"{exact}/{N-1}, each with lambda = -2*gamma*k + i*J*N/2")
        ev = np.linalg.eigvals(L)
        at_max = ev[np.abs(np.abs(ev.imag) - J * N / 2) < 1e-9]
        check(f"N={N}: realising-mode count = 4(N-1) with multiplicity", len(at_max) == 4 * (N - 1),
              f"{len(at_max)} = 2 polarised extremes x (N-1) rungs x 2 signs of Im")
        distinct = {(round(z.real, 8), round(z.imag, 8)) for z in at_max}
        check(f"N={N}: distinct values = 2(N-1), each of multiplicity 2",
              len(distinct) == 2 * (N - 1) and len(at_max) == 2 * len(distinct),
              f"{len(distinct)} distinct, {len(at_max)} with multiplicity")

    print("\n4b. Why 4(N-1): the star's E_min sits in EVERY rung; other graphs' does not")
    for N in (4, 5, 6):
        r = emin_rungs(N, star(N))
        check(f"N={N} star: E_min occupies rungs 1..N-1", r == list(range(1, N)), f"rungs {r}")
        for nm, bonds in (("chain", chain(N)), ("ring", ring(N)), ("complete", complete(N))):
            r = emin_rungs(N, bonds)
            check(f"N={N} {nm}: E_min occupies only the middle rung(s)",
                  r != list(range(1, N)) and len(r) < N - 1, f"rungs {r}")

    print("\n4c. NEGATIVE CONTROL: extremal-multiplet products are mostly NOT eigenoperators")
    for N in (3, 4, 5):
        J, g = 1.0, 0.5
        H = heisenberg(N, star(N), J)
        L = dephasing_liouvillian(H, N, [g] * N)
        w, V = np.linalg.eigh(H)
        Vmax = V[:, np.abs(w - w.max()) < 1e-9]
        Vmin = V[:, np.abs(w - w.min()) < 1e-9]
        worst, total, bad = 0.0, 0, 0
        for a in range(Vmax.shape[1]):
            for b in range(Vmin.shape[1]):
                v = np.kron(Vmax[:, a], np.conj(Vmin[:, b]))
                v /= np.linalg.norm(v)
                r = float(np.linalg.norm(L @ v - np.vdot(v, L @ v) * v))
                total += 1
                worst = max(worst, r)
                if r > 1e-9:
                    bad += 1
        check(f"N={N}: most extremal products fail to be eigenoperators", bad >= total // 2,
              f"{bad}/{total} fail, worst residual {worst:.4f}")

    # --------------------------------------- 5. scope: uniform gamma is required
    print("\n5. SCOPE: site-dependent gamma breaks the saturation (the bound survives)")
    for N, prof in ((3, [2.0, 0.5, 0.5]), (4, [2.0, 0.5, 0.5, 0.5]), (4, [0.1, 0.2, 0.3, 0.4]),
                    (5, [0.11, 0.29, 0.47, 0.65, 0.83])):
        im = max_imag(heisenberg(N, star(N)), N, prof)
        check(f"N={N} gamma={prof}: strictly below J*N/2, bound intact",
              im < N / 2 - 1e-6, f"max|Im|={im:.9f} vs J*N/2={N/2}")

    print("\n5a. CONTROLLED: at equal perturbation the HUB costs more than a leaf")
    for N in (4, 5):
        for delta in (0.05, 0.2):
            gh = [0.5 + delta] + [0.5] * (N - 1)
            gl = [0.5] * (N - 1) + [0.5 + delta]
            ih = max_imag(heisenberg(N, star(N)), N, gh)
            il = max_imag(heisenberg(N, star(N)), N, gl)
            check(f"N={N} delta={delta}: hub-perturbed falls further than leaf-perturbed",
                  ih < il, f"hub {ih:.9f} < leaf {il:.9f}  (J*N/2 = {N/2})")

    print("\n5b. SCOPE: the polarised extreme is required (XY and transverse-field Ising fail)")
    for N in (4, 5):
        d = 2 ** N
        Hxy = np.zeros((d, d), dtype=complex)
        for (i, j) in chain(N):
            Hxy += 0.25 * (op(X, i, N) @ op(X, j, N) + op(Y, i, N) @ op(Y, j, N))
        im, s = max_imag(Hxy, N, [0.5] * N), spread(Hxy)
        check(f"N={N} XY chain does NOT saturate", im < s - 1e-6, f"{im:.7f} vs spread {s:.7f}")

        Htfim = np.zeros((d, d), dtype=complex)
        for (i, j) in chain(N):
            Htfim += 0.25 * (op(Z, i, N) @ op(Z, j, N))
        for i in range(N):
            Htfim += 0.5 * op(X, i, N)
        im, s = max_imag(Htfim, N, [0.5] * N), spread(Htfim)
        check(f"N={N} transverse-field Ising does NOT saturate", im < s - 1e-6,
              f"{im:.7f} vs spread {s:.7f}")

    print("\n5c. SCOPE: Hermitian jumps are required -- otherwise even the BOUND fails")
    c = I2 + 1j * Y
    cd = c.conj().T
    Lc = np.kron(c, c.conj()) - 0.5 * (np.kron(cd @ c, I2) + np.kron(I2, (cd @ c).T))
    evc = np.linalg.eigvals(Lc)
    check("H = 0, c = I + iY: max|Im| > 0 = max H gap", np.abs(evc.imag).max() > 1.9,
          f"spectrum {np.round(np.sort_complex(evc), 9)}")
    check("  ... and it is a legitimate trace-preserving generator",
          np.linalg.norm(Lc.conj().T @ I2.reshape(-1)) < 1e-12 and (evc.real <= 1e-12).all())
    check("  ... its dissipator is NOT Hilbert-Schmidt self-adjoint",
          np.linalg.norm(Lc - Lc.conj().T) > 1.0,
          f"||D-D*||={np.linalg.norm(Lc - Lc.conj().T):.4f}")

    # ---------------------- 6. J*N/2 is the MINIMUM over connected graphs, N <= 6
    print("\n6. The star is the unique minimiser of Delta_max over connected graphs")
    for N in (4, 5, 6):
        edges = list(combinations(range(N), 2))
        best, winners, total = None, [], 0
        for r in range(N - 1, len(edges) + 1):
            for bonds in combinations(edges, r):
                if not connected(N, bonds):
                    continue
                total += 1
                g = spread(heisenberg(N, list(bonds)))
                if best is None or g < best - 1e-9:
                    best, winners = g, [bonds]
                elif abs(g - best) < 1e-9:
                    winners.append(bonds)
        stars = {frozenset(frozenset(e) for e in [(h, k) for k in range(N) if k != h])
                 for h in range(N)}
        won = {frozenset(frozenset(e) for e in wgraph) for wgraph in winners}
        check(f"N={N}: min over {total} connected labelled graphs = J*N/2, minimisers = the {N} stars",
              abs(best - N / 2) < 1e-9 and won == stars, f"min={best:.10f}, {len(winners)} minimisers")

    # ------------------------------- 7. closed forms in the cross-topology table
    print("\n7. The cross-topology closed forms")
    for N in (4, 5, 6):
        kn = spread(heisenberg(N, complete(N)))
        pred = N * (N + 2) / 8 if N % 2 == 0 else (N - 1) * (N + 3) / 8
        check(f"K_{N} spread = J*N(N+2)/8 (even) / J*(N-1)(N+3)/8 (odd)",
              abs(kn - pred) < 1e-10, f"{kn:.10f} vs {pred:.10f}")
    for N, cf, nm in ((4, 0.75, "3/4"), (6, (5 + np.sqrt(13)) / 12, "(5+sqrt13)/12")):
        c_N = spread(heisenberg(N, ring(N))) / N
        check(f"ring c_{N} = {nm}", abs(c_N - cf) < 1e-10, f"{c_N:.12f} vs {cf:.12f}")

    print("\n7b. The ring coefficient approaches ln 2 in two branches")
    cs = {N: spread(heisenberg(N, ring(N))) / N for N in range(3, 13)}
    ln2 = float(np.log(2))
    evens = [cs[N] for N in (4, 6, 8, 10, 12)]
    odds = [cs[N] for N in (3, 5, 7, 9, 11)]
    check("even-N branch lies ABOVE ln 2 and decreases",
          all(c > ln2 for c in evens) and all(a > b for a, b in zip(evens, evens[1:])),
          f"{[f'{c:.6f}' for c in evens]} vs ln2={ln2:.6f}")
    check("odd-N branch lies BELOW ln 2 and increases",
          all(c < ln2 for c in odds) and all(a < b for a, b in zip(odds, odds[1:])),
          f"{[f'{c:.6f}' for c in odds]}")

    print("\n7c. Disconnected components: the spread ADDS")
    N = 8
    k4_p4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (4, 5), (5, 6), (6, 7)]
    tot = spread(heisenberg(N, k4_p4))
    g_k4 = spread(heisenberg(4, complete(4)))
    g_p4 = spread(heisenberg(4, chain(4)))
    check("K_4 (+) P_4 spread = sum of component spreads", abs(tot - (g_k4 + g_p4)) < 1e-9,
          f"{tot:.10f} = {g_k4:.10f} + {g_p4:.10f}  (product would be {g_k4*g_p4:.10f})")
    check("  ... and dividing by sigma reproduces the repo's Im/sigma = 1.3415063509",
          abs(tot / (N * 0.5) - 1.3415063509) < 1e-9, f"{tot/(N*0.5):.10f}")

    # ----------------------------------------- 8. precision of the stored anchors
    print("\n8. Stored-anchor precision (machine precision, NOT bit-exact)")
    worst = 0.0
    files = sorted((REPO / "simulations/results/q_sweep_anchor").glob("star_N*_Q*.json"))
    check("24 star Q-sweep files present", len(files) == 24, f"{len(files)} files")
    for f in files:
        d = json.loads(f.read_text(encoding="utf-8-sig"))
        worst = max(worst, abs(d["MaxImag"] - d["JValue"] * d["N"] / 2) / (d["JValue"] * d["N"] / 2))
    for name in ("star_N3_python", "star_N4_python", "star_N5_python", "star_N6_python", "star_N8"):
        d = json.loads((REPO / f"simulations/results/f1_n8_n9_metrics/{name}.json")
                       .read_text(encoding="utf-8-sig"))
        worst = max(worst, abs(d["MaxImag"] - d["JValue"] * d["N"] / 2) / (d["JValue"] * d["N"] / 2))
    check("no stored anchor is exactly J*N/2 (so 'bit-exact' is the wrong word)", worst > 0.0)
    check("all 29 stored anchors agree to < 5e-14", worst < 5e-14,
          f"worst relative deviation {worst:.3e}")
    print(f"\n  -> worst relative deviation across all 29 stored runs: {worst:.3e}")

    passed = sum(1 for _, ok in CHECKS if ok)
    print(f"\n{'=' * 70}\n{passed}/{len(CHECKS)} checks passed")
    if passed != len(CHECKS):
        print("FAILURES:")
        for name, ok in CHECKS:
            if not ok:
                print(f"  - {name}")
    return 0 if passed == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
