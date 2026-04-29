"""Systematic data sweep for Π²-odd universality at higher N.

Hypothesis: chain (X,Y) ≡ chain (X,Z) ≡ chain (Y,X) ≡ chain (Z,X) at M-spectrum
level for all N. We've verified at N=3, 4, 5. Sweep N=3..7 across topologies
to gather more data, look for emergent patterns, and check if universality
holds or breaks at higher N.

Output: CSV with cluster patterns. Pattern analysis: cluster values as
functions of N, whether they follow a fittable formula.

Run from repo root: python simulations/_pi2_odd_universality_data_sweep.py
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np
import time
from pathlib import Path

import framework as fw
from framework.pauli import _build_bilinear, bit_b
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual


def chain_bonds(N): return [(i, i + 1) for i in range(N - 1)]
def star_bonds(N): return [(0, i) for i in range(1, N)]
def ring_bonds(N): return chain_bonds(N) + [(N - 1, 0)]


def build_M(N, terms, bonds, gamma=1.0):
    H = _build_bilinear(N, bonds, [(t[0], t[1], 1.0) for t in terms])
    L = lindbladian_z_dephasing(H, [gamma] * N)
    return palindrome_residual(L, N * gamma, N)


def cluster_svs(M, tol=1e-6):
    svs = np.linalg.svd(M, compute_uv=False)
    out = []
    for s in svs:
        placed = False
        for i, (v, c) in enumerate(out):
            if abs(s - v) < tol:
                out[i] = (v, c + 1); placed = True; break
        if not placed:
            out.append((s, 1))
    return sorted(out, key=lambda x: -x[0])


def main():
    pi2_odd_pairs = [('X', 'Y'), ('X', 'Z'), ('Y', 'X'), ('Z', 'X')]
    pi2_even_non_truly = [('Y', 'Z'), ('Z', 'Y')]  # for contrast

    topologies = {
        'chain': chain_bonds,
        'star': star_bonds,
        'ring': ring_bonds,
    }

    out_dir = Path('simulations/results/pi2_odd_universality_sweep')
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("Π²-odd 2-body universality sweep")
    print("=" * 78)
    print()

    rows = []  # (N, topology, P, Q, parity, cluster_pattern_string, num_distinct, max_sv, frob_sq)

    # Maximum N — depends on memory. 4^N × 4^N matrix for SVD.
    # N=6: 4^6 = 4096, matrix 16.7M entries × 16 bytes = 268 MB. ~5s SVD.
    # N=7: 4^7 = 16384, matrix 268M entries × 16 bytes = 4.3 GB. ~3 min SVD.
    # For first sweep: stick to N=3..6.
    N_max = 6

    for N in range(3, N_max + 1):
        for topo_name, topo_func in topologies.items():
            bonds = topo_func(N)
            print(f"--- N={N}, topology={topo_name}, {len(bonds)} bonds ---")
            for label, pair_list, parity in [('odd', pi2_odd_pairs, 1),
                                              ('even', pi2_even_non_truly, 0)]:
                results = []
                for (P, Q) in pair_list:
                    t0 = time.time()
                    M = build_M(N, [(P, Q)], bonds)
                    cl = cluster_svs(M)
                    elapsed = time.time() - t0
                    cluster_str = ', '.join(f'({v:.4f},{m})' for v, m in cl)
                    frob_sq = float(np.linalg.norm(M) ** 2)
                    rows.append({
                        'N': N, 'topology': topo_name, 'P': P, 'Q': Q,
                        'parity': label, 'num_distinct': len(cl),
                        'max_sv': cl[0][0] if cl else 0.0,
                        'frob_sq': frob_sq, 'cluster': cluster_str,
                    })
                    results.append((P, Q, cluster_str))
                    print(f"  ({P},{Q}) {label}: {cluster_str} ({elapsed:.1f}s)")

                # Check universality within this parity class
                if len(set(r[2] for r in results)) == 1:
                    print(f"  → {label} universality: YES")
                else:
                    print(f"  → {label} universality: NO ({len(set(r[2] for r in results))} distinct patterns)")
            print()

    # Write CSV
    csv_path = out_dir / 'sweep_results.csv'
    with open(csv_path, 'w', encoding='utf-8') as f:
        f.write('N,topology,P,Q,parity,num_distinct,max_sv,frob_sq,cluster\n')
        for r in rows:
            f.write(f"{r['N']},{r['topology']},{r['P']},{r['Q']},{r['parity']},"
                    f"{r['num_distinct']},{r['max_sv']:.6f},{r['frob_sq']:.4f},"
                    f"\"{r['cluster']}\"\n")
    print(f"Written: {csv_path}")

    # Print summary
    print()
    print("=" * 78)
    print("Summary: where does universality hold/fail?")
    print("=" * 78)
    from collections import defaultdict
    by_class = defaultdict(list)
    for r in rows:
        key = (r['N'], r['topology'], r['parity'])
        by_class[key].append((r['P'], r['Q'], r['cluster']))

    for (N, topo, par), results in sorted(by_class.items()):
        n_unique_clusters = len(set(r[2] for r in results))
        verdict = "UNIVERSAL" if n_unique_clusters == 1 else f"{n_unique_clusters} CLASSES"
        print(f"  N={N} {topo} {par}: {verdict} ({len(results)} pairs tested)")


if __name__ == '__main__':
    main()
