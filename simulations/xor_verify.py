"""Retired verifier for the historical XOR state-coordinate diagnostic.

The former checks treated squared coordinates in a non-normal right-eigenvector
basis as state probabilities. They were also mutation-insensitive when one of
the compared coordinate totals vanished. They are intentionally not retained
as callable gates.

Run ``python simulations/f22_operator_charge.py`` for the independent
operator-level F22 check.
"""


if __name__ == "__main__":
    print(
        "RETIRED: right-eigenvector coordinate squares are not invariant state "
        "weights. Run f22_operator_charge.py for the F22 operator check."
    )
