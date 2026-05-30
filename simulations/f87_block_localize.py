#!/usr/bin/env python3
"""Localize the F87 hard break: which degenerate-frequency block carries the asymmetry?

The windowed converse (non-bipartite => hard) is reduced (PROOF_F103 §7.3,
BIPARTITE_CHIRALITY_DIAGONAL_CELL) to a first-order-in-γ statement: the degenerate dephasing
block D̂ (diagonalized within each degenerate-ω subspace of the H-eigenbasis coherences
|E_a><E_b|) reproduces the break bit-exact, and c = 0 ⟺ bipartite.

This probe sharpens the reduction analytically and checks it. Because H is REAL (real hopping in
the dephasing basis), the (−ω)-block matrix is the transpose of the ω-block matrix, so they share
a spectrum. The palindrome partner of μ = −iω + γs is −μ − 2σ = −i(−ω) + γ(−s − 2N) (σ = Nγ).
With s_{−ω} = s_ω as sets, the pairing collapses to a PER-BLOCK condition:

    {s} = {−s − 2N}   for each degenerate-ω block  ⟺  the block spectrum is symmetric about −N
                                                   ⟺  Zovl_ω := M_ω + N·I is symmetric about 0.

So "soft" ⟺ EVERY block's shift-spectrum is symmetric about −N. A hard pair must have at least
one block asymmetric about −N. This probe prints, per case, the per-block asymmetry, confirming
(a) it sums to the measured c, and (b) WHERE (which ω, which block) the odd cycle breaks it.
"""
from __future__ import annotations
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain, site_op


def blocks(pair, N=4, tol=1e-6):
    """Per-ω block matrices M_ω of the degenerate first-order dephasing shift operator."""
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    groups = defaultdict(list)
    for a in range(d):
        for b in range(d):
            groups[round(E[a] - E[b], 6)].append((a, b))
    out = {}
    for omega, modes in groups.items():
        n = len(modes)
        M = np.zeros((n, n), dtype=complex)
        for i, (a, b) in enumerate(modes):
            for j, (ap, bp) in enumerate(modes):
                val = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
                if (a, b) == (ap, bp):
                    val -= N
                M[i, j] = val
        out[omega] = (modes, M)
    return out, H


def asym_about(vals, center):
    """Optimal-transport asymmetry of the multiset {vals} about `center`:
    mean |v_i − reflect(matched partner)| pairing v ↔ 2·center − v."""
    vals = np.asarray(vals)
    tgt = 2.0 * center - vals
    cost = np.abs(vals[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].mean())


def hopping_graph_bipartite(H, tol=1e-9):
    """2-colour H's off-diagonal hopping graph in the (computational) dephasing basis.
    Returns (is_bipartite, n_nodes_in_largest_component, has_diag)."""
    d = H.shape[0]
    has_diag = any(abs(H[i, i]) > tol for i in range(d))
    adj = defaultdict(list)
    for a in range(d):
        for b in range(d):
            if a != b and abs(H[a, b]) > tol:
                adj[a].append(b)
    color = {}
    bip = True
    for s in range(d):
        if s in color or not adj[s]:
            continue
        color[s] = 0
        stack = [s]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in color:
                    color[w] = color[u] ^ 1
                    stack.append(w)
                elif color[w] == color[u]:
                    bip = False
    return bip, len(color), has_diag


def main():
    N = 4
    cases = [
        ("SOFT  XXZ+ZXX", [('X', 'X', 'Z'), ('Z', 'X', 'X')]),
        ("FLUX  IXY+XIY", [('I', 'X', 'Y'), ('X', 'I', 'Y')]),
        ("REAL  XXZ+XZX", [('X', 'X', 'Z'), ('X', 'Z', 'X')]),
    ]
    for label, pair in cases:
        blk, H = blocks(pair, N)
        bip, nc, has_diag = hopping_graph_bipartite(H)
        # measured c for cross-check
        g = 1e-4
        L = lindbladian_pauli_dephasing(H, [g] * N, dephase_letter='Z')
        ev = np.linalg.eigvals(L)
        c_meas = asym_about(ev, -N * g) / g

        print(f"\n{label}   bipartite={bip}  (graph nodes={nc}, diagonal-lift={has_diag})")
        print(f"  measured c (about −σ) = {c_meas:.4f}")
        print(f"  {'ω':>8}  {'block':>5}  {'asym about −N':>13}")
        total = 0.0
        worst = (None, 0.0, 0)
        for omega in sorted(blk, key=lambda w: (abs(w), w)):
            modes, M = blk[omega]
            s = np.linalg.eigvals(M).real
            a = asym_about(s, -N)
            total += a * len(s)              # weight by block size (each eigenvalue is a mode)
            if a > worst[1]:
                worst = (omega, a, len(s))
            if a > 1e-9:
                print(f"  {omega:8.3f}  {len(s):5d}  {a:13.4f}   <- breaks")
        # normalise the per-block sum the same way opt_break averages (per eigenvalue)
        n_modes = sum(len(m) for m, _ in blk.values())
        print(f"  Σ(per-block asym, size-weighted)/d² = {total / n_modes:.4f}   "
              f"(worst block: ω={worst[0]}, size {worst[2]}, asym {worst[1]:.4f})")

        # The analytic backbone: M(−ω) = M(ω)^T (H real, Z symmetric), so the +ω and −ω blocks
        # share a spectrum. That is why the palindrome partner (ω,s) ↔ (−ω, −s−2N) collapses to a
        # PER-BLOCK condition {s} = {−s−2N} (symmetry about −N), rather than the +ω block pairing
        # with a genuinely different −ω partner. The equal-asymmetry ± pairs above are this fact.
        tmis = 0.0
        for omega in blk:
            if omega <= 1e-9:
                continue
            key_m = round(-omega, 6)
            if key_m in blk:
                sp = np.sort(np.linalg.eigvals(blk[omega][1]).real)
                sm = np.sort(np.linalg.eigvals(blk[key_m][1]).real)
                if len(sp) == len(sm):
                    tmis = max(tmis, float(np.max(np.abs(sp - sm))))
        print(f"  transpose M(−ω)=M(ω)^T: max ±ω spectrum mismatch = {tmis:.1e}")


if __name__ == "__main__":
    main()
