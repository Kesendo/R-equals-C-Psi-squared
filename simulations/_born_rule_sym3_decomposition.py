#!/usr/bin/env python3
"""Decompose the Dyson sym3 = 8 integer for F94 by (bond1, bond2, site, ordering, component).

F94 (Tier 1, 2026-05-16) states Δ_|00⟩ = (4/3)·Q²·K³ for the dominant outcome of |0+0+⟩
N=4 Heisenberg ring + Z-dephasing on pair (0,2). The "8" in the coefficient 4/3 = 8/6
comes from

    ⟨00|_pair Tr_{1,3}[ sym3·ρ_0 ] |00⟩_pair = 8

where sym3 = L_H² L'_dis + L_H L'_dis L_H + L'_dis L_H². The proof PROOF_F94 closes the
bit-exact derivation but leaves open the *structural* decomposition of 8 by bonds,
sites, orderings, and Heisenberg components (XX/YY/ZZ).

This script enumerates all (b1, b2, s, ordering, c1, c2) contributions and reports
which non-vanishing combinations sum to 8. We want to see the structure: do specific
bonds dominate? do XX/YY/ZZ contribute equally? what is the site-ordering pattern?
"""
from __future__ import annotations

import sys
from pathlib import Path
from collections import defaultdict
from fractions import Fraction

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
EYE = np.eye(2, dtype=complex)
PAULI = {"X": SX, "Y": SY, "Z": SZ, "I": EYE}


def single_site_op(N, i, op):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        out = np.kron(out, op if k == i else EYE)
    return out


def two_site_op(N, i, j, op_a, op_b):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        if k == i:
            out = np.kron(out, op_a)
        elif k == j:
            out = np.kron(out, op_b)
        else:
            out = np.kron(out, EYE)
    return out


def commutator(A, rho):
    return -1j * (A @ rho - rho @ A)


def Z_dis_site(rho, N, s):
    """Single-site Z-dephasing dissipator: Z_s ρ Z_s − ρ."""
    Zs = single_site_op(N, s, SZ)
    return Zs @ rho @ Zs - rho


def reduced_density(rho, N, keep):
    n_keep = len(keep)
    trace = [i for i in range(N) if i not in keep]
    rho_tensor = rho.reshape([2] * N + [2] * N)
    for q in sorted(trace, reverse=True):
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + (N - sum(1 for t in trace if t > q)))
    return rho_tensor.reshape((2 ** n_keep, 2 ** n_keep))


def initial_state(N):
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi_0 = np.kron(zero, np.kron(plus, np.kron(zero, plus)))
    return np.outer(psi_0, psi_0.conj())


def heisenberg_bond_component(N, a, b, component):
    """Single component of bond Hamiltonian: (J/4)·σ^c_a σ^c_b with J=1."""
    op = PAULI[component]
    return 0.25 * two_site_op(N, a, b, op, op)


def evaluate_triple(N, rho_0, b1, b2, s, ordering, c1, c2):
    """Evaluate ⟨00|_pair Tr_{1,3}[ ordering(L_H(H_{b1,c1}), L_H(H_{b2,c2}), L_s) ρ_0 ] |00⟩_pair.

    ordering ∈ {1, 2, 3}:
        1 = L_H² L'_dis  →  apply L_s, then L_H_{b2,c2}, then L_H_{b1,c1}
        2 = L_H L'_dis L_H  →  apply L_H_{b2,c2}, then L_s, then L_H_{b1,c1}
        3 = L'_dis L_H²  →  apply L_H_{b2,c2}, then L_H_{b1,c1}, then L_s
    """
    bonds = [(i, (i + 1) % N) for i in range(N)]
    a1, ap1 = bonds[b1]
    a2, ap2 = bonds[b2]
    H_b1c1 = heisenberg_bond_component(N, a1, ap1, c1)
    H_b2c2 = heisenberg_bond_component(N, a2, ap2, c2)

    if ordering == 1:
        # L_H L_H L_s ρ
        x = Z_dis_site(rho_0, N, s)
        x = commutator(H_b2c2, x)
        x = commutator(H_b1c1, x)
    elif ordering == 2:
        # L_H L_s L_H ρ
        x = commutator(H_b2c2, rho_0)
        x = Z_dis_site(x, N, s)
        x = commutator(H_b1c1, x)
    elif ordering == 3:
        # L_s L_H L_H ρ
        x = commutator(H_b2c2, rho_0)
        x = commutator(H_b1c1, x)
        x = Z_dis_site(x, N, s)
    else:
        raise ValueError(ordering)

    x_pair = reduced_density(x, N, keep=[0, 2])
    val = x_pair[0, 0]
    if abs(val.imag) > 1e-12:
        raise RuntimeError(f"Non-real |00⟩ element: {val}")
    return float(val.real)


def main():
    N = 4
    rho_0 = initial_state(N)
    bonds = [(i, (i + 1) % N) for i in range(N)]
    sites = list(range(N))
    components = ["X", "Y", "Z"]
    orderings = [1, 2, 3]

    print("Decomposing F94 sym3 = 8")
    print(f"  Setup: |0+0+⟩ N=4 Heisenberg ring, pair (0,2), |00⟩ outcome")
    print(f"  Bonds: {bonds}")
    print(f"  Sites: {sites}  (|0⟩ on 0, 2; |+⟩ on 1, 3)")
    print()

    contributions = defaultdict(float)  # key = (b1, b2, s, ordering, c1, c2)
    total = 0.0
    nonzero_terms = 0
    sampled_terms = 0

    for ordering in orderings:
        for b1 in range(len(bonds)):
            for b2 in range(len(bonds)):
                for s in sites:
                    for c1 in components:
                        for c2 in components:
                            v = evaluate_triple(N, rho_0, b1, b2, s, ordering, c1, c2)
                            sampled_terms += 1
                            if abs(v) > 1e-12:
                                key = (b1, b2, s, ordering, c1, c2)
                                contributions[key] = v
                                total += v
                                nonzero_terms += 1

    print(f"  Sampled {sampled_terms} (b1, b2, s, ord, c1, c2) terms.")
    print(f"  Non-vanishing terms: {nonzero_terms}")
    print(f"  Total = {total:.6f}  (expect 8.0)")
    print()

    # Aggregate by various axes
    print("=" * 70)
    print("Sum by ordering")
    print("=" * 70)
    by_ord = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_ord[ordering] += v
    for ordering in orderings:
        print(f"  ordering {ordering}: {by_ord[ordering]:+.6f}")
    print()

    print("=" * 70)
    print("Sum by site")
    print("=" * 70)
    by_site = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_site[s] += v
    for s in sites:
        marker = "  |0⟩" if s in (0, 2) else "  |+⟩"
        print(f"  site {s}: {by_site[s]:+.6f}  {marker}")
    print()

    print("=" * 70)
    print("Sum by component pair (c1, c2)")
    print("=" * 70)
    by_cc = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_cc[(c1, c2)] += v
    for c1 in components:
        for c2 in components:
            print(f"  ({c1}, {c2}): {by_cc[(c1, c2)]:+.6f}")
    print()

    print("=" * 70)
    print("Sum by bond pair (b1, b2)")
    print("=" * 70)
    by_bb = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_bb[(b1, b2)] += v
    for b1 in range(len(bonds)):
        for b2 in range(len(bonds)):
            v = by_bb[(b1, b2)]
            if abs(v) > 1e-9:
                print(f"  bond_pair (b1={b1}={bonds[b1]}, b2={b2}={bonds[b2]}): {v:+.6f}")
    print()

    print("=" * 70)
    print("Sum by (component_pair, site) — the structural cell")
    print("=" * 70)
    by_cc_s = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_cc_s[(c1, c2, s)] += v
    for c1 in components:
        for c2 in components:
            row = []
            for s in sites:
                row.append(f"s={s}: {by_cc_s[(c1, c2, s)]:+.4f}")
            print(f"  ({c1}, {c2})  " + "  ".join(row))
    print()

    print("=" * 70)
    print("Sum by (ordering, component_pair) — the dynamic cell")
    print("=" * 70)
    by_ord_cc = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_ord_cc[(ordering, c1, c2)] += v
    for ordering in orderings:
        for c1 in components:
            for c2 in components:
                v = by_ord_cc[(ordering, c1, c2)]
                if abs(v) > 1e-9:
                    print(f"  ord={ordering}  ({c1}, {c2}): {v:+.6f}")
    print()

    # Now: scaled by (1/4)^2 since each H_b carries J/4 = 1/4 at J=1
    # The Heisenberg coupling we ought to use is (J/4)·(XX+YY+ZZ), so each H_b has 3
    # components each with coefficient 1/4. The (1/4)² is already baked in via
    # heisenberg_bond_component returning 0.25·op.  Let's verify the total matches.

    # Equivalent count after dividing out the (1/4)² absorbs:
    raw_count = total * 16  # undo two (1/4) factors
    print(f"  After undoing (1/4)² factor from heisenberg_bond_component:")
    print(f"  raw integer count = {raw_count:.3f}  (expect a clean integer)")
    print()

    # Now do the same per (ordering, c1, c2)
    print("=" * 70)
    print("Raw integer count (×16 to undo the 1/16 = (1/4)² coupling normalization)")
    print("=" * 70)
    for ordering in orderings:
        for c1 in components:
            for c2 in components:
                v = by_ord_cc[(ordering, c1, c2)] * 16
                if abs(v) > 1e-9:
                    print(f"  ord={ordering}  ({c1}, {c2}): {v:+.3f}")
    print()

    print("=" * 70)
    print("Component-pair (c1, c2) totals × 16")
    print("=" * 70)
    for c1 in components:
        for c2 in components:
            v = by_cc[(c1, c2)] * 16
            print(f"  ({c1}, {c2}): {v:+.3f}")
    print()

    print("=" * 70)
    print("(ordering, c1, c2) × bond_pair breakdown")
    print("=" * 70)
    by_ord_cc_bb = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_ord_cc_bb[(ordering, c1, c2, b1, b2)] += v
    for ordering in orderings:
        for c1 in components:
            for c2 in components:
                cell_total = by_ord_cc[(ordering, c1, c2)]
                if abs(cell_total) < 1e-9:
                    continue
                print(f"  ord={ordering}  ({c1}, {c2}) total = {cell_total:+.4f}")
                for b1 in range(len(bonds)):
                    for b2 in range(len(bonds)):
                        v = by_ord_cc_bb[(ordering, c1, c2, b1, b2)]
                        if abs(v) > 1e-9:
                            shared = "self" if b1 == b2 else "adj@" + str(set(bonds[b1]) & set(bonds[b2]))
                            print(f"    b1={b1}={bonds[b1]}, b2={b2}={bonds[b2]}  ({shared}): {v:+.4f}")
                print()

    print("=" * 70)
    print("(ordering, c1, c2) × site breakdown")
    print("=" * 70)
    by_ord_cc_s = defaultdict(float)
    for (b1, b2, s, ordering, c1, c2), v in contributions.items():
        by_ord_cc_s[(ordering, c1, c2, s)] += v
    for ordering in orderings:
        for c1 in components:
            for c2 in components:
                cell_total = by_ord_cc[(ordering, c1, c2)]
                if abs(cell_total) < 1e-9:
                    continue
                print(f"  ord={ordering}  ({c1}, {c2}) total = {cell_total:+.4f}")
                for s in sites:
                    v = by_ord_cc_s[(ordering, c1, c2, s)]
                    marker = "(kept |0⟩)" if s in (0, 2) else "(traced |+⟩)"
                    print(f"    s={s}  {marker}: {v:+.4f}")
                print()

    print("=" * 70)
    print("Per-term enumeration of all 32 surviving (b1, b2, s, ord, c1, c2) terms")
    print("=" * 70)
    sorted_terms = sorted(contributions.items())
    for (b1, b2, s, ordering, c1, c2), v in sorted_terms:
        v_raw = v * 16  # × 16 to undo the (1/4)² Heisenberg normalization
        print(f"  ord={ordering}  ({c1},{c2})  b1={b1}={bonds[b1]} b2={b2}={bonds[b2]}  s={s}  "
              f": v={v:+.4f}  (raw ×16: {v_raw:+.1f})")
    print()
    print(f"  Sum of per-term raw values: {sum(v * 16 for v in contributions.values()):+.4f}")
    print(f"  Sum of per-term values:     {sum(contributions.values()):+.4f}")
    print()

    print("=" * 70)
    print("Bucket by absolute value (raw ×16)")
    print("=" * 70)
    bucket = defaultdict(int)
    for v in contributions.values():
        v_raw = round(v * 16, 1)
        bucket[v_raw] += 1
    for val, count in sorted(bucket.items()):
        print(f"  raw value {val:+.1f}: {count} terms, total = {val * count:+.1f}")
    print()


if __name__ == "__main__":
    main()
