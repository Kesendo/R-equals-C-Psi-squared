#!/usr/bin/env python3
"""F94 topology visibility probe: when does F94 see K_4 vs ring vs chain?

Investigates the question: F50's central-weight excess sees K_4's extra
bonds (over ring) as a +23 weight-2 ker excess at N=4. F94's sym3 matrix
element on |0+0+⟩ pair (0,2) is a different observable on the same algebra.
Does F94 see the same K_4 vs ring topology dependence?

Findings (2026-05-17):

  (1) For F94's canonical lens (|0+0+⟩, pair (0,2)), K_4 and Ring give
      sym3 matrix elements agreeing in this finite 1e-9 comparison. That
      numerical agreement does not establish exact equality and does not
      identify a cancellation mechanism for the two extra K_4 bonds.

  (2) Ring |10> slope -16/9, while chain |10> slope -4/3. This refutes
      topology universality. Ring and K_4 agree only for the named lens;
      that agreement does not extend to the chain or to arbitrary states.

  (3) For asymmetric initial states (|++00⟩, |10+0⟩): F94's sym3 IS
      sensitive to K_4 vs ring. The K_4 - Ring difference is an
      antisymmetric integer shift between outcomes that differ by one bit
      flip. Magnitudes are small (1-3) — much smaller than F50's +23
      central-weight excess at K_4 N=4 weight-2.

Conclusion: F94/F96 (Dyson series at fixed initial state) and F50 (operator
algebra independent of initial state) are complementary lenses on
Heisenberg + Z-deph topology dependence. F94 sees observable-specific
integer shifts when the state breaks the relevant symmetry; F50 sees the
full centralizer dimension which integrates over all initial-state choices.
The K_4-vs-ring topology contains structure visible at both levels with
different magnitude relationships.
"""
from __future__ import annotations

import sys
import numpy as np

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def two_site(N, i, j, a, b):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        if k == i:
            out = np.kron(out, a)
        elif k == j:
            out = np.kron(out, b)
        else:
            out = np.kron(out, PAULI["I"])
    return out


def single_site(N, i, op):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        out = np.kron(out, op if k == i else PAULI["I"])
    return out


def heisenberg(N, bonds, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for (a, b) in bonds:
        H += (J / 4) * (
            two_site(N, a, b, PAULI["X"], PAULI["X"])
            + two_site(N, a, b, PAULI["Y"], PAULI["Y"])
            + two_site(N, a, b, PAULI["Z"], PAULI["Z"]))
    return H


def L_H_apply(rho, H):
    return -1j * (H @ rho - rho @ H)


def L_dis_apply(rho, N):
    out = np.zeros_like(rho)
    for l in range(N):
        Zl = single_site(N, l, PAULI["Z"])
        out += Zl @ rho @ Zl - rho
    return out


def reduced(rho, N, keep):
    n_keep = len(keep)
    trace = [i for i in range(N) if i not in keep]
    t = rho.reshape([2] * N + [2] * N)
    for q in sorted(trace, reverse=True):
        t = np.trace(t, axis1=q, axis2=q + (N - sum(1 for tt in trace if tt > q)))
    return t.reshape((2 ** n_keep, 2 ** n_keep))


def sym3_matrix_elements(N, bonds, state, pair):
    """Compute the sym3·ρ_0 pair-element for all 4 outcomes."""
    H = heisenberg(N, bonds)
    rho_0 = np.outer(state, state.conj())
    a = L_dis_apply(rho_0, N)
    a = L_H_apply(a, H)
    a = L_H_apply(a, H)
    b = L_H_apply(rho_0, H)
    b = L_dis_apply(b, N)
    b = L_H_apply(b, H)
    c = L_H_apply(rho_0, H)
    c = L_H_apply(c, H)
    c = L_dis_apply(c, N)
    sym3 = a + b + c
    r = reduced(sym3, N, pair)
    return [r[i, i].real for i in range(4)]


def unitary_second_order_elements(N, bonds, state, pair):
    """Return the raw pair diagonal of L_H^2 rho_0."""
    H = heisenberg(N, bonds)
    rho_0 = np.outer(state, state.conj())
    lh2 = L_H_apply(L_H_apply(rho_0, H), H)
    r = reduced(lh2, N, pair)
    return [r[i, i].real for i in range(4)]


def normalized_dyson_slope(dyson_element, order, unitary_element):
    return dyson_element / (order * unitary_element)


def main():
    if sys.platform == "win32":
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    one = np.array([0, 1], dtype=complex)

    test_states = {
        "|0+0+⟩ pair=(0,2) [F94 canonical]": (
            np.kron(zero, np.kron(plus, np.kron(zero, plus))), [0, 2]),
        "|0+0+⟩ pair=(0,1) [pair asym]": (
            np.kron(zero, np.kron(plus, np.kron(zero, plus))), [0, 1]),
        "|++00⟩ pair=(0,2) [state asym]": (
            np.kron(plus, np.kron(plus, np.kron(zero, zero))), [0, 2]),
        "|+0+0⟩ pair=(0,2) [|0+0+⟩ rotated]": (
            np.kron(plus, np.kron(zero, np.kron(plus, zero))), [0, 2]),
        "|10+0⟩ pair=(0,1) [fully asym]": (
            np.kron(one, np.kron(zero, np.kron(plus, zero))), [0, 1]),
    }

    topologies = {
        "Chain": [(0, 1), (1, 2), (2, 3)],
        "Ring":  [(0, 1), (1, 2), (2, 3), (3, 0)],
        "K_4":   [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
    }
    comparison_tolerance = 1e-9

    print("=" * 78)
    print("F94 sym3 visibility of K_4 vs Ring topology, across initial-state choices")
    print("=" * 78)
    print()

    for s_name, (state, pair) in test_states.items():
        print(f"--- {s_name} ---")
        print(f"  {'topology':>10} {'sym3@|00⟩':>10} {'sym3@|01⟩':>10} {'sym3@|10⟩':>10} {'sym3@|11⟩':>10}")
        rows = {}
        for t_name, bonds in topologies.items():
            vals = sym3_matrix_elements(4, bonds, state, pair)
            rows[t_name] = vals
            print(f"  {t_name:>10} " + " ".join(f"{v:>10.4f}" for v in vals))
        diff_RvK = [rows["K_4"][i] - rows["Ring"][i] for i in range(4)]
        visible = any(abs(d) > comparison_tolerance for d in diff_RvK)
        if visible:
            print(f"  → K_4 vs Ring DIFFERENCE: {diff_RvK}  ★ topology visible")
        else:
            print("  → no difference resolved by this finite 1e-9 comparison")
        print()

    print("=" * 78)
    print("Key finding:")
    print("=" * 78)
    print()
    canonical_state, canonical_pair = test_states[
        "|0+0+⟩ pair=(0,2) [F94 canonical]"
    ]
    canonical_ring = sym3_matrix_elements(
        4, topologies["Ring"], canonical_state, canonical_pair
    )
    canonical_k4 = sym3_matrix_elements(
        4, topologies["K_4"], canonical_state, canonical_pair
    )
    canonical_difference = [
        canonical_k4[index] - canonical_ring[index] for index in range(4)
    ]
    max_abs_difference = max(abs(value) for value in canonical_difference)
    print(
        "Named-lens finite 1e-9 comparison: "
        f"max |K_4 - Ring| = {max_abs_difference:.3e}"
    )
    print("This does not establish exact equality and does not identify a cancellation mechanism.")
    print()
    print("F96 topology control for the named lens (raw M3/(3 U2)):")
    for topology_name in ("Chain", "Ring", "K_4"):
        bonds = topologies[topology_name]
        m3 = sym3_matrix_elements(4, bonds, canonical_state, canonical_pair)
        u2 = unitary_second_order_elements(4, bonds, canonical_state, canonical_pair)
        slope_10 = normalized_dyson_slope(m3[2], 3, u2[2])
        print(f"  {topology_name:>5}: |10> slope = {slope_10:+.9f}")
    print("  Chain -4/3 versus ring -16/9 is the explicit topology counterexample.")
    print()
    print("Asymmetric initial states (|++00⟩, |10+0⟩) break the symmetry, and K_4's")
    print("extra bonds become visible as ANTISYMMETRIC INTEGER SHIFTS between outcomes")
    print("that differ by one bit flip. Magnitudes are small (1-3), much smaller than")
    print("F50's +23 central-weight excess at K_4 N=4 weight-2 — different lenses on")
    print("the same underlying topology dependence.")


if __name__ == "__main__":
    main()
