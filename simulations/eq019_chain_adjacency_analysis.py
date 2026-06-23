#!/usr/bin/env python3
"""EQ-019 chain-adjacency analysis.

The Dicke c_1 bond-position scan (c1_bond_scan_multi_N) shows c_1(|S_n>, b)
varies with bond b. At N=5 the factor between endpoint and interior is ~26
for |S_2>, |S_3>.

Hypothesis tested today: chain-adjacency layer (V-Effect / project_v_effect_
combinatorial) — endpoint bonds have 1-sided chain-adjacency, interior bonds
have 2-sided. If the layer story extends from V-Effect categories to the c_1
closure-breaking coefficient, we should see a clean monotone dependence of
|c_1| on chain-adjacency degree.

Reads:
  results/c1_bond_scan/bond_scan.json                    (N=5)
  results/c1_bond_scan_multi_N/bond_scan_multi_N_3_4.json (N=3, 4)
  results/c1_bond_scan_multi_N/bond_scan_multi_N_6.json   (N=6)

For each (N, n), reports:
  - c_1 at each bond
  - bond's chain-adjacency degree (0, 1, 2)
  - mirror-pair averaging where applicable
  - endpoint/interior ratio
  - whether the data supports chain-adjacency monotonicity

Output: print-only (no new compute, no new JSON).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
RESULTS = SCRIPT_DIR / "results"


def load_data():
    """Returns dict: N -> {bond_label: {S_n: c_1}}."""
    data = {}

    p5 = RESULTS / "c1_bond_scan" / "bond_scan.json"
    with open(p5) as f:
        d5 = json.load(f)
    data[5] = d5["bond_results"]

    p34 = RESULTS / "c1_bond_scan_multi_N" / "bond_scan_multi_N_3_4.json"
    with open(p34) as f:
        d34 = json.load(f)
    data[3] = d34["3"]["bond_results"]
    data[4] = d34["4"]["bond_results"]

    p6 = RESULTS / "c1_bond_scan_multi_N" / "bond_scan_multi_N_6.json"
    with open(p6) as f:
        d6 = json.load(f)
    data[6] = d6["6"]["bond_results"]

    return data


def chain_adjacency(b, N):
    """Number of chain-adjacent neighbors for bond b in N-chain.

    Bond b connects sites b and b+1. Its adjacent bonds are (b-1,b) and (b+1,b+2).
    Endpoint bonds (b=0 or b=N-2) have 1 chain-adjacent neighbor.
    Interior bonds have 2.
    For N=2 (only 1 bond), it has 0 chain-adjacent neighbors.
    """
    n_bonds = N - 1
    left = (b - 1 >= 0)
    right = (b + 1 < n_bonds)
    return int(left) + int(right)


def parse_bond_label(label):
    """'bond_0_1' -> (0, 1)."""
    parts = label.split("_")
    return int(parts[1]), int(parts[2])


def main():
    data = load_data()

    print("=" * 80)
    print("EQ-019: Chain-adjacency hypothesis for c_1 bond-position dependence")
    print("=" * 80)
    print()
    print("Hypothesis: c_1(|S_n>, bond b) decreases monotonically with bond's")
    print("chain-adjacency degree (V-Effect compatibility-layer logic).")
    print("  endpoint bonds: 1-sided adjacency")
    print("  interior bonds: 2-sided adjacency")
    print()
    print("If hypothesis holds, expect |c_1| to be larger at endpoints than")
    print("at interior, monotonically across all (N, n).")
    print()

    # For each N, build (b, adj_degree, n) -> c_1 table
    for N in sorted(data.keys()):
        print("-" * 80)
        print(f"N = {N}")
        print("-" * 80)
        print()

        bonds = data[N]
        sorted_bonds = sorted(bonds.keys(), key=parse_bond_label)
        bond_indices = [parse_bond_label(lbl)[0] for lbl in sorted_bonds]
        adj = [chain_adjacency(b, N) for b in bond_indices]

        # Print header
        header = f"  {'bond':<10} {'adj':>4}   "
        states = sorted(bonds[sorted_bonds[0]].keys(),
                        key=lambda s: int(s.split("_")[1]))
        for s in states:
            header += f"{s:>10}"
        print(header)
        for lbl, b, a in zip(sorted_bonds, bond_indices, adj):
            row = f"  {lbl:<10} {a:>4}   "
            for s in states:
                row += f"{bonds[lbl][s]:>+10.4f}"
            print(row)
        print()

        # For each n, group bonds by adjacency degree, report mean |c_1|
        print(f"  Mean |c_1| by chain-adjacency degree (N={N}):")
        adj_to_bonds = {}
        for lbl, a in zip(sorted_bonds, adj):
            adj_to_bonds.setdefault(a, []).append(lbl)
        adj_keys = sorted(adj_to_bonds.keys())
        header2 = f"    {'state':>6}   "
        for a in adj_keys:
            header2 += f" adj={a} (n={len(adj_to_bonds[a])})  "
        print(header2)
        for s in states:
            row = f"    {s:>6}   "
            for a in adj_keys:
                vals = [abs(bonds[lbl][s]) for lbl in adj_to_bonds[a]]
                mean = sum(vals) / len(vals)
                row += f"  {mean:>10.4f}    "
            print(row)
        print()

        # Endpoint / interior ratio per state, only if both adjacency
        # degrees are present
        if 1 in adj_to_bonds and 2 in adj_to_bonds:
            print(f"  Endpoint/Interior ratio (|c_1|; 1-sided / 2-sided):")
            print(f"    {'state':>6}   {'mean endpoint':>16}   "
                  f"{'mean interior':>16}   {'ratio':>10}")
            for s in states:
                ep = sum(abs(bonds[lbl][s])
                          for lbl in adj_to_bonds[1]) / len(adj_to_bonds[1])
                it = sum(abs(bonds[lbl][s])
                          for lbl in adj_to_bonds[2]) / len(adj_to_bonds[2])
                if it < 1e-10:
                    ratio_str = "   inf"
                else:
                    ratio_str = f"{ep / it:>10.2f}"
                # Filter unreliable states (where c_1 is exactly 0 or huge)
                if abs(ep) < 1e-9 and abs(it) < 1e-9:
                    continue
                print(f"    {s:>6}   {ep:>16.4f}   {it:>16.4f}   {ratio_str}")
        print()

    # Cross-N pattern: pull |S_1> (single-excitation symmetric) across all N,
    # endpoint vs interior
    print("=" * 80)
    print("Cross-N: |S_1> Dicke state, endpoint-vs-interior scaling")
    print("=" * 80)
    print()
    print(f"  {'N':>3}   {'|c_1|@endpoint':>16}   {'|c_1|@interior':>16}   "
          f"{'EP/IT ratio':>14}")
    for N in sorted(data.keys()):
        bonds = data[N]
        sorted_bonds = sorted(bonds.keys(), key=parse_bond_label)
        bond_indices = [parse_bond_label(lbl)[0] for lbl in sorted_bonds]
        adj = [chain_adjacency(b, N) for b in bond_indices]
        adj_to_bonds = {}
        for lbl, a in zip(sorted_bonds, adj):
            adj_to_bonds.setdefault(a, []).append(lbl)
        if 1 not in adj_to_bonds or 2 not in adj_to_bonds:
            continue
        s = "S_1"
        if s not in bonds[sorted_bonds[0]]:
            continue
        ep = sum(abs(bonds[lbl][s]) for lbl in adj_to_bonds[1]) / len(adj_to_bonds[1])
        it = sum(abs(bonds[lbl][s]) for lbl in adj_to_bonds[2]) / len(adj_to_bonds[2])
        ratio = ep / it if it > 1e-10 else float('inf')
        print(f"  {N:>3d}   {ep:>16.4f}   {it:>16.4f}   {ratio:>14.2f}")
    print()

    # Same for |S_2>
    print("Cross-N: |S_2> Dicke state, endpoint-vs-interior scaling")
    print()
    print(f"  {'N':>3}   {'|c_1|@endpoint':>16}   {'|c_1|@interior':>16}   "
          f"{'EP/IT ratio':>14}")
    for N in sorted(data.keys()):
        bonds = data[N]
        sorted_bonds = sorted(bonds.keys(), key=parse_bond_label)
        bond_indices = [parse_bond_label(lbl)[0] for lbl in sorted_bonds]
        adj = [chain_adjacency(b, N) for b in bond_indices]
        adj_to_bonds = {}
        for lbl, a in zip(sorted_bonds, adj):
            adj_to_bonds.setdefault(a, []).append(lbl)
        if 1 not in adj_to_bonds or 2 not in adj_to_bonds:
            continue
        s = "S_2"
        if s not in bonds[sorted_bonds[0]]:
            continue
        ep = sum(abs(bonds[lbl][s]) for lbl in adj_to_bonds[1]) / len(adj_to_bonds[1])
        it = sum(abs(bonds[lbl][s]) for lbl in adj_to_bonds[2]) / len(adj_to_bonds[2])
        ratio = ep / it if it > 1e-10 else float('inf')
        print(f"  {N:>3d}   {ep:>16.4f}   {it:>16.4f}   {ratio:>14.2f}")
    print()

    # Verdict
    print("=" * 80)
    print("Verdict")
    print("=" * 80)
    print()
    print("Read the EP/IT ratio columns above:")
    print("  - If ratio > 1 monotonically across N for all states:")
    print("    chain-adjacency hypothesis HOLDS (1-sided > 2-sided in |c_1|).")
    print("  - If ratio swings around 1 or inverts state-by-state:")
    print("    bond-position dependence is NOT explained by chain-adjacency")
    print("    layer alone. Additional structure (state-dependent sine-mode")
    print("    overlap) is at play.")


if __name__ == "__main__":
    main()
