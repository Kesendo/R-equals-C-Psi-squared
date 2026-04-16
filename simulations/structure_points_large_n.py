#!/usr/bin/env python3
"""
Structure Points at Large N
============================
Extend the structure-points scan from N=3,4,5 to N=6,7 (chains)
and test non-chain topologies (ring, star, Y-junction).

Question: are {0, gamma_0, 2*gamma_0} the only universal anchors?

Date: 2026-04-16
"""

import numpy as np
import sys
from scipy.linalg import eig
from fractions import Fraction
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)
OUT_PATH = RESULTS_DIR / "structure_points_large_n.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

np.set_printoptions(precision=8, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron(*factors)


def liouvillian(H, jump_ops):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Idd, LdL)
              - 0.5 * np.kron(LdL.T, Idd))
    return L


def build_hamiltonian(N, bonds, J=1.0):
    """Build XX+YY Hamiltonian from a list of (i, j) bonds."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        H += J * 0.5 * (site_op(X, i, N) @ site_op(X, j, N)
                        + site_op(Y, i, N) @ site_op(Y, j, N))
    return H


def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]


def ring_bonds(N):
    return [(i, (i + 1) % N) for i in range(N)]


def star_bonds(N):
    """Site 0 is the hub, 1..N-1 are leaves."""
    return [(0, i) for i in range(1, N)]


def y_junction_bonds(N):
    """Three arms of roughly equal length meeting at site 0."""
    bonds = []
    arm_len = (N - 1) // 3
    site = 1
    for arm in range(3):
        prev = 0
        for _ in range(arm_len):
            bonds.append((prev, site))
            prev = site
            site += 1
    # remaining sites on last arm
    while site < N:
        bonds.append((prev, site))
        prev = site
        site += 1
    return bonds


def distinct_alphas(L_super, tol=1e-8):
    eigs = np.linalg.eigvals(L_super)
    alphas = -eigs.real
    sorted_alphas = np.sort(alphas)
    distinct = []
    counts = []
    for a in sorted_alphas:
        if distinct and abs(a - distinct[-1]) < tol:
            counts[-1] += 1
        else:
            distinct.append(a)
            counts.append(1)
    return distinct, counts


def best_fraction(x, max_denom=24):
    if abs(x) < 1e-9:
        return Fraction(0)
    return Fraction(x).limit_denominator(max_denom)


def scan_structure(N, bonds, topology_name, gamma_0, B_site, J=1.0):
    """Run structure-point scan for given topology."""
    log(f"\n{'=' * 70}")
    log(f"N={N}, topology={topology_name}, B={B_site}, "
        f"dim={4**N}, gamma_0={gamma_0}")
    log(f"{'=' * 70}")

    H = build_hamiltonian(N, bonds, J)
    L = liouvillian(H, [np.sqrt(gamma_0) * site_op(Z, B_site, N)])
    distinct, counts = distinct_alphas(L, tol=gamma_0 * 1e-3)
    normalized = [a / gamma_0 for a in distinct]

    log(f"Distinct rates: {len(distinct)}")
    log(f"{'alpha/g0':>12} {'mult':>5} {'rational':>10} {'err':>10}")
    log("-" * 42)
    fracs = set()
    for a, c, x in zip(distinct, counts, normalized):
        f = best_fraction(x, max_denom=24)
        f_val = float(f)
        err = abs(x - f_val)
        if err < 1e-3:
            fracs.add(f)
        log(f"{x:12.6f} {c:5d} {str(f):>10} {err:10.2e}")

    return fracs, normalized


# =====================================================================
if __name__ == "__main__":
    log("STRUCTURE POINTS AT LARGE N")
    log("Are {0, gamma_0, 2*gamma_0} the only universal anchors?")

    J = 1.0
    gamma_0 = 1e-4  # smaller than original scan for cleaner rationals

    # === Part 1: Chains N=3..7 ===
    log("\n\n" + "#" * 70)
    log("PART 1: CHAINS (B = last site)")
    log("#" * 70)

    chain_fracs = {}
    for N in range(3, 8):
        if 4**N > 20000:
            log(f"\nN={N}: dim={4**N}, using scipy.linalg.eig (may be slow)...")
        fracs, _ = scan_structure(N, chain_bonds(N), "chain", gamma_0, N - 1)
        chain_fracs[N] = fracs

    # cross-N universal
    log("\n\n" + "=" * 70)
    log("CHAIN: Cross-N universal fractions (N=3..7)")
    log("=" * 70)

    all_chain_fracs = sorted(set().union(*chain_fracs.values()))
    universal_chain = set.intersection(*chain_fracs.values())

    header = f"{'fraction':>12} {'value':>8} "
    for N in sorted(chain_fracs.keys()):
        header += f"{'N='+str(N):>5}"
    header += f"  {'UNIV':>5}"
    log(header)

    for f in all_chain_fracs:
        val = float(f)
        if val < -0.01 or val > 2.01:
            continue
        in_all = f in universal_chain
        row = f"{str(f):>12} {val:8.4f} "
        for N in sorted(chain_fracs.keys()):
            present = "y" if f in chain_fracs[N] else "-"
            row += f"{present:>5}"
        row += f"  {'YES' if in_all else '':>5}"
        log(row)

    log(f"\nUniversal: {sorted(universal_chain)}")

    # === Part 2: Non-chain topologies at N=5 ===
    log("\n\n" + "#" * 70)
    log("PART 2: NON-CHAIN TOPOLOGIES (N=5, B=last site)")
    log("#" * 70)

    topo_fracs = {}
    topologies = [
        ("chain", chain_bonds(5)),
        ("ring", ring_bonds(5)),
        ("star", star_bonds(5)),
        ("Y-junction", y_junction_bonds(5)),
    ]

    for name, bonds in topologies:
        fracs, _ = scan_structure(5, bonds, name, gamma_0, 4)
        topo_fracs[name] = fracs

    log("\n\n" + "=" * 70)
    log("TOPOLOGY: Cross-topology universal fractions (N=5)")
    log("=" * 70)

    all_topo_fracs = sorted(set().union(*topo_fracs.values()))
    universal_topo = set.intersection(*topo_fracs.values())

    header = f"{'fraction':>12} {'value':>8} "
    for name in topo_fracs:
        header += f"{name[:6]:>7}"
    header += f"  {'UNIV':>5}"
    log(header)

    for f in all_topo_fracs:
        val = float(f)
        if val < -0.01 or val > 2.01:
            continue
        in_all = f in universal_topo
        row = f"{str(f):>12} {val:8.4f} "
        for name in topo_fracs:
            present = "y" if f in topo_fracs[name] else "-"
            row += f"{present:>7}"
        row += f"  {'YES' if in_all else '':>5}"
        log(row)

    log(f"\nUniversal across topologies: {sorted(universal_topo)}")

    # === Part 3: Different B-site positions ===
    log("\n\n" + "#" * 70)
    log("PART 3: B-SITE POSITION (N=5 chain)")
    log("#" * 70)

    bsite_fracs = {}
    for b in range(5):
        fracs, _ = scan_structure(5, chain_bonds(5), f"chain,B={b}", gamma_0, b)
        bsite_fracs[b] = fracs

    universal_bsite = set.intersection(*bsite_fracs.values())
    log(f"\nUniversal across B-sites: {sorted(universal_bsite)}")

    # === Verdict ===
    log("\n\n" + "=" * 70)
    log("VERDICT")
    log("=" * 70)
    log(f"Universal across chains N=3..7: {sorted(universal_chain)}")
    log(f"Universal across topologies:     {sorted(universal_topo)}")
    log(f"Universal across B-sites:        {sorted(universal_bsite)}")

    _outf.close()
    print(f"\nResults written to {OUT_PATH}")
