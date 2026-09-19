"""Finite checks that keep several appearances of the numeral four separate.

The operator-space identity (d^2)^N, a Hamming-distance dephasing ladder, and a
constructed flavour-hopping model are three different objects. Their values at
d=2 do not give one typed genealogy. In particular, the sourced F94 qudit
continuation gives 40/27 at d=3 rather than 3: the qutrit continuation refutes
typed a_-1 ancestry. Its qubit coefficient 4/3 comes from the 8/6
setup-specific diagram count. This is a sourced counterexample, not a
derivation in this script, and no qudit fold is inferred.

The executable rows below are finite sanity checks for the definitions used in
this file. They do not identify F94 with F56, F97, F98, or a quarter boundary.
Importing the module is silent.
"""

import itertools
from fractions import Fraction

import numpy as np


GAMMA = 0.05
J = 1.0


def all_states(d, N):
    return list(itertools.product(range(d), repeat=N))


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def dephasing_rate_ladder(d, N):
    """Distinct rates in the stipulated equidistant full-Cartan model."""
    states = all_states(d, N)
    rates = {
        round(2.0 * GAMMA * hamming(a, b), 12)
        for a in states
        for b in states
    }
    return sorted(rates)


def se_band(d, N):
    """Distinct energies of the stipulated one-flavour-at-a-time hopping model."""
    states = [
        state
        for state in all_states(d, N)
        if sum(value != 0 for value in state) == 1
    ]
    idx = {state: index for index, state in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for a, state in enumerate(states):
        for left in range(N - 1):
            right = left + 1
            if (state[left] == 0) != (state[right] == 0):
                moved = list(state)
                moved[left], moved[right] = moved[right], moved[left]
                s2 = tuple(moved)
                H[idx[s2], a] += J
    return sorted({round(value, 9) for value in np.linalg.eigvalsh(H)})


def se_band_multiplicity(d, N):
    """Flavour multiplicity in the constructed hopping model."""
    if N < 1 or d < 2:
        raise ValueError("requires N>=1 and d>=2")
    return d - 1


def f94_qudit_coefficient(d):
    """Sourced F94 continuation used here solely as an ancestry counterexample."""
    if d < 2:
        raise ValueError("requires d>=2")
    return Fraction(4 * (d + 2) * (d - 1), 3 * d * d)


def main():
    print("=== FINITE FOUR-COEFFICIENT CLASSIFICATION CHECKS ===")
    print("Distinct objects can share the value four at d=2.")

    print("\nStage A: exact operator-space dimensions")
    for N in (1, 2, 3, 4):
        qubit_dimension = 2 ** (2 * N)
        qutrit_dimension = 3 ** (2 * N)
        assert qubit_dimension == 4 ** N
        assert qutrit_dimension == 9 ** N
        print(f"  N={N}: d=2 -> {qubit_dimension}; d=3 -> {qutrit_dimension}")

    print("\nStage C: finite rate-ladder sanity check")
    stage_c_rows = []
    for N in (2, 3, 4):
        lad2 = dephasing_rate_ladder(2, N)
        lad3 = dephasing_rate_ladder(3, N)
        expected = [round(2.0 * GAMMA * k, 12) for k in range(N + 1)]
        assert lad2 == expected and lad3 == expected
        stage_c_rows.append((N, lad2, lad3, expected))
    for N, lad2, lad3, expected in stage_c_rows:
        print(f"  N={N}: d2={lad2}; d3={lad3}; stipulated={expected}")
    print("This confirms only the leading model's defined Hamming-rate ladder.")

    print("\nStage D: constructed flavour-hopping model")
    stage_d_rows = []
    for N in (3, 4, 5):
        b2 = se_band(2, N)
        b3 = se_band(3, N)
        multiplicity2 = se_band_multiplicity(2, N)
        multiplicity3 = se_band_multiplicity(3, N)
        rho = 2.0 * np.cos(np.pi / (N + 1))
        assert b2 == b3
        assert abs(max(abs(b2[0]), abs(b2[-1])) - rho) < 1e-7
        stage_d_rows.append((N, b2, b3, multiplicity2, multiplicity3, rho))
    for N, b2, b3, multiplicity2, multiplicity3, rho in stage_d_rows:
        print(
            f"  N={N}: distinct bands equal={b2 == b3}; "
            f"flavour multiplicities={multiplicity2},{multiplicity3}; radius={rho:.9f}"
        )

    f94_qutrit = f94_qudit_coefficient(3)
    assert f94_qutrit == Fraction(40, 27)
    print("\nF94 qutrit continuation: "
          f"{f94_qutrit} = 40/27; refutes typed a_-1 ancestry.")
    print("This is a sourced counterexample, not a derivation in this script.")
    print("The qubit 4/3 is the named 8/6 setup-specific diagram count.")
if __name__ == "__main__":
    main()
