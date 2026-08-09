#!/usr/bin/env python3
"""
Gamma_0 as Eigenvalue: Exact-Hit Positions
==========================================
Correction scan for the gamma_0-as-eigenvalue feature of the dissipation
interval [0, 2*gamma_0] (hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md, the
"dissipation interval" section).

The original scan (structure_points_large_n.py) ran at gamma_0 = 1e-4 with
bin tolerance gamma_0 * 1e-3, which cannot separate an exact eigenvalue hit
from an O((gamma_0/J)^2) perturbative near-miss: both bin to alpha/gamma_0 = 1.
This script asks the sharp question instead: is gamma_0 an EXACT decay rate
(double-precision hit, |alpha - gamma_0| < 1e-12 * max(1, gamma_0)) of the
full Liouvillian, and how close is the nearest miss otherwise? Run at three
well-separated gamma_0 values so a hit must be gamma-robust to count.

Found 2026-08-09 (review of PRIMORDIAL_GAMMA_CONSTANT): the original doc
sentence had the chain map INVERTED. Measured here: at the odd sizes N=3,5
the exact hits sit at the odd interior positions (N=3: B=1; N=5: B=1,3),
never at the endpoints; at the even sizes N=4,6 every position hits, with
gamma-DEPENDENT multiplicity (the presence is robust, the counts are not).
Star hub at N=5 collapses the whole spectrum to alpha/gamma_0 in {0, 1, 2}
at gamma_0/J = 0.1 and 0.3; at 0.7 the collapse opens palindromically while
{0, 1, 2} stay exact. The Y-junction, reported absent in the April doc,
hits at the central node, the long-arm joint, and the long-arm end,
gamma-robustly; the short-arm ends, ring, and star leaf carry no hit.

Conventions: H = sum_bonds (J/2)(XX+YY) (single-excitation hopping = J),
dephasing on the single site B via jump operator sqrt(gamma_0) * Z_B,
matching structure_points_large_n.py. The Y-junction is that script's
shape at N=5: bonds (0,1),(0,2),(0,3),(3,4), center 0, long-arm end 4.
"""

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)
OUT_PATH = RESULTS_DIR / "gamma0_eigenvalue_positions.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, site, n):
    out = np.array([[1]], dtype=complex)
    for s in range(n):
        out = np.kron(out, op if s == site else I2)
    return out


def build_hamiltonian(n, bonds, J=1.0):
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        H += J * 0.5 * (site_op(X, i, n) @ site_op(X, j, n)
                        + site_op(Y, i, n) @ site_op(Y, j, n))
    return H


def liouvillian_single_site_dephasing(H, n, b, gamma):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    Zb = site_op(Z, b, n)
    # row-stacking vec: rho -> A rho B^T maps to kron(A, B) on vec(rho)
    L = -1j * (np.kron(H, Idd) - np.kron(Idd, H.T))
    L += gamma * (np.kron(Zb, Zb.T) - np.kron(Idd, Idd))
    return L


def chain_bonds(n):
    return [(i, i + 1) for i in range(n - 1)]


def star_bonds(n):
    return [(0, i) for i in range(1, n)]


def ring_bonds(n):
    return [(i, (i + 1) % n) for i in range(n)]


def y_junction_bonds_n5():
    return [(0, 1), (0, 2), (0, 3), (3, 4)]


def exact_hit_report(n, bonds, b, gamma):
    H = build_hamiltonian(n, bonds)
    ev = np.linalg.eigvals(liouvillian_single_site_dephasing(H, n, b, gamma))
    alpha = -ev.real
    tol = 1e-12 * max(1.0, gamma)
    hits = int(np.sum(np.abs(alpha - gamma) < tol))
    misses = np.abs(alpha - gamma)
    misses = misses[misses >= tol]
    nearest = float(np.min(misses)) / gamma if misses.size else float("nan")
    return hits, nearest


def main():
    gammas = (0.1, 0.3, 0.7)
    log("gamma_0 exact-eigenvalue positions (hit = |alpha - gamma_0| < 1e-12*max(1,gamma_0))")
    log(f"gamma_0 values tested: {gammas}; J = 1 on every bond")
    log()

    log("== Chains N = 3..6, every dephasing position B ==")
    for n in range(3, 7):
        for b in range(n):
            per_gamma = [exact_hit_report(n, chain_bonds(n), b, g) for g in gammas]
            hits = {h for (h, _) in per_gamma}
            robust = hits if len(hits) == 1 else None
            nearest = max(m for (_, m) in per_gamma)
            log(f"N={n} chain B={b}: exact hits per gamma = "
                f"{[h for (h, _) in per_gamma]}"
                f"{' (gamma-robust)' if robust else ' (NOT gamma-robust!)'}"
                f", nearest relative miss <= {nearest:.3e}")
        log()

    log("== N=5 non-chain topologies ==")
    cases = [
        ("star B=hub", star_bonds(5), 0),
        ("star B=leaf", star_bonds(5), 1),
        ("ring B=0 (all sites equivalent)", ring_bonds(5), 0),
        ("Y-junction B=4 (long-arm end)", y_junction_bonds_n5(), 4),
        ("Y-junction B=3 (long-arm joint)", y_junction_bonds_n5(), 3),
        ("Y-junction B=1 (short-arm end; B=2 equivalent)", y_junction_bonds_n5(), 1),
        ("Y-junction B=0 (central node)", y_junction_bonds_n5(), 0),
    ]
    for name, bonds, b in cases:
        per_gamma = [exact_hit_report(5, bonds, b, g) for g in gammas]
        nearest = max(m for (_, m) in per_gamma)
        log(f"N=5 {name}: exact hits per gamma = {[h for (h, _) in per_gamma]}"
            f", nearest relative miss <= {nearest:.3e}")

    log()
    log("== Chain endpoint near-misses per gamma (B=0; the nearest rate to gamma_0) ==")
    for n in (3, 5):
        parts = []
        for g in gammas:
            _, near = exact_hit_report(n, chain_bonds(n), 0, g)
            parts.append(f"gamma_0={g}: {near:.4e}")
        log(f"N={n} chain B=0: " + " | ".join(parts))

    log()
    log("== N=5 chain, gamma_0=0.1: distinct eigenvalue counts and minimum gaps ==")
    for b in (0, 1, 2):
        H = build_hamiltonian(5, chain_bonds(5))
        ev = np.linalg.eigvals(liouvillian_single_site_dephasing(H, 5, b, 0.1))
        alpha = np.sort(-ev.real)
        distinct = [alpha[0]]
        for a in alpha[1:]:
            if a - distinct[-1] > 1e-9:
                distinct.append(a)
        gaps = np.diff(distinct)
        log(f"B={b}: {len(distinct)} distinct alpha (clustering tol 1e-9), "
            f"min gap {np.min(gaps):.3e}")

    log()
    log("== N=5 star hub: full distinct alpha/gamma_0 sets ==")
    for g in gammas:
        H = build_hamiltonian(5, star_bonds(5))
        ev = np.linalg.eigvals(liouvillian_single_site_dephasing(H, 5, 0, g))
        alpha = np.sort(-ev.real) / g
        distinct = []
        for a in alpha:
            if not distinct or a - distinct[-1] > 1e-10:
                distinct.append(a)
        log(f"gamma_0={g}: distinct alpha/gamma_0 = "
            + ", ".join(f"{a:.12f}" for a in distinct))

    log()
    log("== N=3 chain, B=2 (endpoint): first-order vs measured S-coherence rates ==")
    log("(the (1exc, vac) coherence block is M = -i*H1 - 2*gamma_0*P_B;")
    log(" first order: alpha = 2*gamma_0*|a_B|^2 with a_B the H1 eigenvector;")
    log(" exact: alpha = -Re eig(M) = 2*gamma_0*|v(B)|^2 with v the M eigenvector, EQ-015)")
    g = 0.1
    H1 = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=complex)  # hopping J=1
    PB = np.diag([0, 0, 1]).astype(complex)
    evals, evecs = np.linalg.eigh(H1)
    M = -1j * H1 - 2 * g * PB
    mevals, mvecs = np.linalg.eig(M)
    order = np.argsort(evals)
    for idx in order:
        aB2 = abs(evecs[2, idx]) ** 2
        overlaps = np.abs(mvecs.conj().T @ evecs[:, idx]) ** 2
        k = int(np.argmax(overlaps))
        exact = -mevals[k].real
        vB2 = abs(mvecs[2, k]) ** 2 / np.sum(np.abs(mvecs[:, k]) ** 2)
        log(f"H1 mode E={evals[idx]:+.6f}: |a_B|^2={aB2:.6f}, "
            f"first-order alpha={2 * g * aB2:.6f}, measured alpha={exact:.12f}, "
            f"alpha/(2*gamma_0)={exact / (2 * g):.12f} (=|v(B)|^2 up to norm: {vB2:.12f})")

    log()
    log("== N=3 chain, B=2: first-order error law (max rel err of 2*gamma_0*|a_B|^2")
    log("   against -Re eig(M), relative to the first-order prediction) ==")
    for g in (1e-4, 1e-3, 1e-2, 0.1, 0.5):
        M = -1j * H1 - 2 * g * PB
        mevals2, mvecs2 = np.linalg.eig(M)
        worst = 0.0
        for idx in range(3):
            aB2 = abs(evecs[2, idx]) ** 2
            overlaps = np.abs(mvecs2.conj().T @ evecs[:, idx]) ** 2
            exact = -mevals2[int(np.argmax(overlaps))].real
            pred = 2 * g * aB2
            worst = max(worst, abs(exact - pred) / pred)
        log(f"gamma_0/J={g}: max rel err = {worst:.4e}")


if __name__ == "__main__":
    main()
