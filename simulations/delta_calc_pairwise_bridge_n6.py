"""
The N = 6 ring of Dynamic Entanglement §9.1, in the retired tool's reading and exactly
=====================================================================================
experiments/DYNAMIC_ENTANGLEMENT.md §9.1 reports the retired delta_calc tool's
simulate_subsystem_crossing on the alternating state |0+0+0+> on the six-site Heisenberg ring,
gamma 0.05, t_max 5 (the Claude Desktop log of the chat's MCP holds the call, 2026-02-18 11:29 UTC,
and keeps its head and its tail, crossing_count 0, but not the maxima). The tool returns each
maximum rounded to four digits (round(max_C_psi, 4) in its source); the transcription reproduces the
page's maxima of this run, 0.248, 0.131, 0.113 and 0.111, from that four-digit output read at three
digits, which ties it to that run (the tool's source read here is dated after the call). The tool's source is kept outside the repo and is
transcribed here, not imported: an Euler step of dt 0.01, then the Hermitian part, the negative
eigenvalues set to zero and the trace renormalized; every pair (i, j) read in the pairwise
correlation bridge C_corr = clip((P_AB - P_A P_B)/(1 - P_A P_B), 0, 1) times Psi = l1/3 of the
pair's reduced state, at every step. simulations/delta_calc_pairwise_bridge.py does the same at
N = 4; this script adds N = 6, and the exact N = 4 maxima and crossing times of pair (0,2) that
§9.1 and §11.2 set beside the tool's.

Conventions: H = sum over ring bonds of (XX + YY + ZZ) (J = 1, the Pauli book), local jumps Z_k at
rate gamma, D[rho] = gamma sum_k (Z_k rho Z_k - rho); qubit 0 the most significant bit; the
alternating state puts |0> on the even sites and |+> on the odd ones. Exact means DOP853 at rtol
1e-10 on the density matrix, read every 0.001 (the claims over all 15 pairs every 0.01).

Import-inert; prints only, and exits with status 1 if any check fails.
"""

from itertools import combinations
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
ZERO = np.array([1, 0], dtype=complex)
PLUS = np.array([1, 1], dtype=complex) / np.sqrt(2)
CLAIMS, CHECKS = [], []


def site(P, k, n):
    m = np.array([[1.0]], dtype=complex)
    for i in range(n):
        m = np.kron(m, P if i == k else I2)
    return m


def ring(n):
    H = np.zeros((2**n, 2**n), dtype=complex)
    for i in range(n):
        for P in (X, Y, Z):
            H += site(P, i, n) @ site(P, (i + 1) % n, n)
    return H


def alternating(n):
    psi = np.array([1.0], dtype=complex)
    for k in range(n):
        psi = np.kron(psi, ZERO if k % 2 == 0 else PLUS)
    return np.outer(psi, psi.conj())


def reduced(rho, keep, n):
    t = rho.reshape([2] * (2 * n))
    for i in sorted((q for q in range(n) if q not in keep), reverse=True):
        t = np.trace(t, axis1=i, axis2=i + t.ndim // 2)
    return t.reshape(2 ** len(keep), 2 ** len(keep))


def purity(r):
    return float(np.real(np.trace(r @ r)))


def cpsi(rho, i, j, n):
    rab = reduced(rho, [i, j], n)
    pa, pb, pab = purity(reduced(rho, [i], n)), purity(reduced(rho, [j], n)), purity(rab)
    den = 1.0 - pa * pb
    c = min(1.0, max(0.0, (pab - pa * pb) / den)) if den > 1e-12 else 0.0
    return c * float(np.abs(rab).sum() - np.abs(np.diag(rab)).sum()) / 3.0


def rhs(H, n, gamma):
    Ls = [site(Z, k, n) for k in range(n)]
    d = 2**n

    def f(_, v):
        r = v.reshape(d, d)
        dr = -1j * (H @ r - r @ H)
        for L in Ls:
            dr += gamma * (L @ r @ L - r)
        return dr.reshape(-1)

    return f, Ls


def tool_maxima(n, gamma, t_max, pairs, dt=0.01):
    H = ring(n)
    _, Ls = rhs(H, n, gamma)
    rho = alternating(n)
    best = {p: 0.0 for p in pairs}
    count = 0
    prev = {p: cpsi(rho, *p, n) for p in pairs}
    for step in range(int(t_max / dt) + 1):
        for p in pairs:
            v = cpsi(rho, *p, n)
            best[p] = max(best[p], v)
            if step and (prev[p] - 0.25) * (v - 0.25) < 0:
                count += 1
            prev[p] = v
        if step == int(t_max / dt):
            break
        d = -1j * (H @ rho - rho @ H)
        for L in Ls:
            d += gamma * (L @ rho @ L - rho)
        rho = rho + dt * d
        rho = 0.5 * (rho + rho.conj().T)
        w, v = np.linalg.eigh(rho)
        rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
        rho /= np.real(np.trace(rho))
    return best, count


def exact_solution(n, gamma, t_max):
    f, _ = rhs(ring(n), n, gamma)
    return solve_ivp(f, (0, t_max), alternating(n).reshape(-1), method="DOP853", rtol=1e-10, atol=1e-12,
                     dense_output=True)


def claim(label, ok, detail=""):
    CLAIMS.append((label, ok))
    print(f"  {label} {detail}  {'holds' if ok else 'FAILS'}")


def check(label, page, value, digits=3):
    ok = round(value, digits) == round(page, digits)
    CHECKS.append((label, ok))
    print(f"  {label:<58s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")


def main():
    n, gamma, t_max = 6, 0.05, 5.0
    shown = [(0, 1), (0, 2), (1, 3), (0, 3)]
    print("N = 6 ring, |0+0+0+>, gamma 0.05, t <= 5: largest C Psi of each pair")
    best, count = tool_maxima(n, gamma, t_max, list(combinations(range(n), 2)))
    for p, former in zip(shown, (0.248, 0.131, 0.113, 0.111)):
        check(f"the tool, pair {p}: its four-digit output, at three", former, round(best[p], 4))
    claim("the tool: no pair crosses 1/4 (no sign change of C Psi - 1/4 at any step, largest below 1/4)",
          count == 0 and max(best.values()) < 0.25, f"(largest {max(best.values()):.4f})")
    claim("the tool ranks the non-neighbour classes (0,2) > (1,3) > (0,3), the antipodal last",
          best[(0, 2)] > best[(1, 3)] > best[(0, 3)])
    classes = {(0, 1): [(1, 2), (2, 3), (3, 4), (4, 5), (0, 5)], (0, 2): [(2, 4), (0, 4)], (1, 3): [(3, 5), (1, 5)],
               (0, 3): [(1, 4), (2, 5)]}
    covered = sorted([p for p in classes] + [q for qs in classes.values() for q in qs])
    claim("the four classes cover the 15 pairs once each", covered == list(combinations(range(n), 2)))
    tsym = max(abs(best[q] - best[p]) for p, qs in classes.items() for q in qs)
    claim("the tool: the pairs of each class carry the same maximum (<= 1e-9)", tsym <= 1e-9, f"({tsym:.1e})")
    sol = exact_solution(n, gamma, t_max)
    ts = np.arange(0, t_max + 1e-9, 0.001)
    ex = {p: max(cpsi(sol.sol(t).reshape(64, 64), *p, n) for t in ts) for p in shown}
    for p, page in zip(shown, (0.248, 0.122, 0.102, 0.129)):
        check(f"exactly, pair {p}", page, ex[p])
    allex = max(max(cpsi(sol.sol(t).reshape(64, 64), *p, n) for t in ts[::10]) for p in combinations(range(n), 2))
    claim("exactly, no pair reaches 1/4 (sampled every 0.01 for all 15 pairs)", allex < 0.25, f"(largest {allex:.4f})")
    sym = max(abs(max(cpsi(sol.sol(t).reshape(64, 64), *q, n) for t in ts[::10])
                  - max(cpsi(sol.sol(t).reshape(64, 64), *p, n) for t in ts[::10]))
              for p, qs in classes.items() for q in qs)
    claim("exactly: the pairs of each class carry the same maximum (sampled every 0.01, <= 1e-9)", sym <= 1e-9, f"({sym:.1e})")
    claim("exactly, the antipodal class (0,3) has the largest C Psi of the three non-neighbour classes",
          ex[(0, 3)] > ex[(0, 2)] > ex[(1, 3)])

    print("\nN = 4 ring, |0+0+>, gamma 0.05, exactly (Sec. 9.1's N = 4 column and Sec. 11.2)")
    sol4 = exact_solution(4, 0.05, t_max)
    for p, page in (((0, 1), 0.247), ((0, 2), 0.320), ((1, 3), 0.224)):
        check(f"pair {p}, largest C Psi", page, max(cpsi(sol4.sol(t).reshape(16, 16), *p, 4) for t in ts))
    g = lambda t: cpsi(sol4.sol(t).reshape(16, 16), 0, 2, 4) - 0.25
    check("pair (0, 2), t_cross_up", 0.285, brentq(g, 0.20, 0.35))
    check("pair (0, 2), t_cross_down", 0.453, brentq(g, 0.40, 0.60))

    held = sum(ok for _, ok in CLAIMS)
    matched = sum(ok for _, ok in CHECKS)
    print(f"\n{matched} of {len(CHECKS)} printed values match at the page's precision; "
          f"{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS) or matched < len(CHECKS):
        sys.exit(1)


if __name__ == "__main__":
    main()
