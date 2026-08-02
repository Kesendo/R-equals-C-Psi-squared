#!/usr/bin/env python3
"""The cavity finesse law of experiments/VEFFECT_CAVITY_MODES.md, from below and against its dump.

On the single-excitation sector the Heisenberg Hamiltonian is exactly H = J(N-1)*I - 2J*L_graph:
the ZZ term supplies -2J*deg(j) and turns the adjacency matrix into the graph LAPLACIAN. That
block is invariant (H conserves magnetization) and the dephasing is flat at -2gamma across it, so

    lambda_m = -2 gamma - 2i J mu_m        and        Q_max = J * mu_max / gamma

with mu the Laplacian spectrum of the coupling graph. For the path P_N that is
mu_max = 2(1 + cos(pi/N)); for the star it is N.

This script builds the Laplacians itself and compares against simulations/results/
veffect_cavity_modes.txt, which is the measurement. Two things it deliberately does NOT do:
it does not take the document's closed form on trust (it derives Q from mu_max, so a wrong
closed form in the prose cannot make the gate pass), and it does not assume the highest-Q mode
sits at weight 1 (the dump's own claim). Where an expected N or topology is missing from the
dump it FAILS rather than dropping the row, which an earlier version of this file did silently.

The dump prints Q rounded to one decimal but carries the eigenvalue at full precision on the
same line, so the comparison uses Im(lambda) where it can: at N=5 the printed Q alone spends
more than half of a 1e-3 tolerance on display rounding.
"""
from __future__ import annotations
import os
import re
import sys
import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = os.path.join(HERE, "results", "veffect_cavity_modes.txt")
J, GAMMA = 1.0, 0.05


def path_laplacian(n):
    L = np.zeros((n, n))
    for i in range(n - 1):
        L[i, i + 1] = L[i + 1, i] = -1.0
        L[i, i] += 1.0
        L[i + 1, i + 1] += 1.0
    return L


def star_laplacian(n):
    L = np.zeros((n, n))
    for i in range(1, n):
        L[0, i] = L[i, 0] = -1.0
        L[0, 0] += 1.0
        L[i, i] += 1.0
    return L


def ring_laplacian(n):
    L = path_laplacian(n)
    L[0, n - 1] = L[n - 1, 0] = -1.0
    L[0, 0] += 1.0
    L[n - 1, n - 1] += 1.0
    return L


LAPLACIAN = {"chain": path_laplacian, "star": star_laplacian, "ring": ring_laplacian}


def q_from_graph(kind, n, J=J, gamma=GAMMA):
    return J * np.linalg.eigvalsh(LAPLACIAN[kind](n))[-1] / gamma


def read_chain_highest_q(path=DUMP):
    """Per N: (Q as printed, |Im lambda| at full precision, the weight shell it sat in)."""
    out, current = {}, None
    for line in open(path, encoding="utf-8"):
        m = re.match(r"\s*N=(\d+):", line)
        if m:
            current = int(m.group(1))
        m = re.search(r"Highest-Q mode: .*?([-+][\d.]+)([-+][\d.]+)j, Q=([\d.eE+-]+), weight shell k=(\d+)",
                      line)
        if m and current is not None:
            out[current] = (float(m.group(3)), abs(float(m.group(2))), int(m.group(4)))
    return out


def read_topology_table(path=DUMP):
    """The Result 5 rows: (N, topology) -> Q_max as printed."""
    out = {}
    for line in open(path, encoding="utf-8"):
        m = re.match(r"\s*(\d+)\s+(chain|star|ring|complete)\s+\d+\s+\d+\s+([\d.]+)\s+[\d.]+\s*$", line)
        if m:
            out[(int(m.group(1)), m.group(2))] = float(m.group(3))
    return out


if __name__ == "__main__":
    if not os.path.exists(DUMP):
        sys.exit(f"missing {DUMP}")
    ok = True
    highest = read_chain_highest_q()
    expected_ns = [2, 3, 4, 5, 6]
    missing = [n for n in expected_ns if n not in highest]
    if missing:
        print(f"  FAIL  the dump is missing N={missing} for the chain")
        ok = False

    print("chain: Q_max from the path Laplacian against the dump")
    print(f"{'N':>3} {'law':>12} {'|Im| based':>12} {'printed Q':>10} {'shell':>6} {'rel dev':>10}")
    for n in sorted(set(expected_ns) & set(highest)):
        printed_q, im_abs, shell = highest[n]
        q_law = q_from_graph("chain", n)
        q_meas = im_abs / (2 * GAMMA)          # Q = |Im| / |Re|, and |Re| = 2 gamma at weight 1
        # The dump prints Im(lambda) to four decimals, so its own resolution is about 7e-6
        # relative at these magnitudes. 2e-5 sits just above that and is still fifty times
        # tighter than comparing against the one-decimal Q column would allow.
        dev = abs(q_meas - q_law) / q_law
        if dev > 2e-5:
            ok = False
        if shell != 1:
            print(f"  NOTE  N={n}: the highest-Q mode sits at weight shell {shell}, not 1;"
                  f" the |Re| = 2 gamma step below does not apply")
            ok = False
        print(f"{n:>3} {q_law:>12.6f} {q_meas:>12.6f} {printed_q:>10.1f} {shell:>6} {dev:>10.1e}")

    print("\ntopology: the same law on other graphs, against the dump's own Result 5 table")
    topo = read_topology_table()
    checked = 0
    for (n, kind), q_dump in sorted(topo.items()):
        if kind not in LAPLACIAN:
            continue
        q_law = q_from_graph(kind, n)
        good = abs(q_dump - q_law) <= 0.05 + 1e-9      # the table prints one decimal
        ok &= good
        checked += 1
        print(f"  {'PASS' if good else 'FAIL'}  N={n} {kind:<5} law {q_law:8.4f}  dump {q_dump:8.1f}")
    if checked == 0:
        print("  FAIL  no topology rows parsed")
        ok = False

    print(f"\nthe chain's finesse limit is 4J/gamma = {4*J/GAMMA:.0f}, since mu_max < 4 on a path.")
    star_over = [(n, q_from_graph('star', n)) for n in (4, 5) if q_from_graph('star', n) >= 4 * J / GAMMA]
    print("  it does NOT bound other graphs: " +
          ", ".join(f"star N={n} reaches {q:.0f}" for n, q in star_over))
    print("\nall checks pass" if ok else "\nFAILURES above")
    sys.exit(0 if ok else 1)
