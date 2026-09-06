"""Operator-level scope check for the former error-correction palindrome study.

This producer intentionally does not decompose states into right Liouvillian
eigenvectors. For a non-normal generator, squared right-eigenvector coordinates
are not invariant probabilities and cannot rank state protection.

The surviving result is F22: a computational-basis coherence |a><b| under
local Z-dephasing is charged at twice the sum of the rates on sites where a
and b differ. The GHZ off-diagonal operator has disagreement on every site.

Output: simulations/results/error_correction_palindrome_f22.txt
"""

from pathlib import Path

import numpy as np


OUT = Path(__file__).resolve().parent / "results" / "error_correction_palindrome_f22.txt"
I2 = np.eye(2, dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)


def basis_operator(n: int, bra: int, ket: int) -> np.ndarray:
    """Return |bra><ket| in the N-qubit computational basis."""
    op = np.zeros((2**n, 2**n), dtype=complex)
    op[bra, ket] = 1.0
    return op


def site_operator(n: int, site: int, local: np.ndarray) -> np.ndarray:
    factors = [I2] * n
    factors[site] = local
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor)
    return result


def dephasing_action(op: np.ndarray, gammas: tuple[float, ...]) -> np.ndarray:
    """Apply sum_l gamma_l (Z_l op Z_l - op)."""
    n = len(gammas)
    result = np.zeros_like(op)
    for site, gamma in enumerate(gammas):
        z_site = site_operator(n, site, Z)
        result += gamma * (z_site @ op @ z_site - op)
    return result


def hamming_charge(bra: int, ket: int, gammas: tuple[float, ...]) -> float:
    """Return 2 sum gamma_l over the differing computational-basis bits."""
    n = len(gammas)
    return 2.0 * sum(
        gamma
        for site, gamma in enumerate(gammas)
        if ((bra >> (n - 1 - site)) & 1) != ((ket >> (n - 1 - site)) & 1)
    )


def verify_operator_charge(
    label: str,
    bra: int,
    ket: int,
    gammas: tuple[float, ...],
) -> tuple[str, int, float]:
    op = basis_operator(len(gammas), bra, ket)
    charge = hamming_charge(bra, ket, gammas)
    actual = dephasing_action(op, gammas)
    expected = -charge * op
    if not np.array_equal(actual, expected):
        raise AssertionError(f"{label}: direct dissipator action != -charge * operator")
    distance = (bra ^ ket).bit_count()
    return label, distance, charge


def main() -> None:
    # Dyadic rates make the two independent routes exactly comparable in IEEE
    # arithmetic: explicit Z conjugation versus the Hamming-disagreement sum.
    gammas = (0.125, 0.25, 0.375)
    n = len(gammas)
    maximum_charge = 2.0 * sum(gammas)

    rows = [
        verify_operator_charge("population |000><000|", 0b000, 0b000, gammas),
        verify_operator_charge("one-site coherence |000><001|", 0b000, 0b001, gammas),
        verify_operator_charge("W/Bell-type coherence |001><100|", 0b001, 0b100, gammas),
        verify_operator_charge("GHZ coherence |000><111|", 0b000, 0b111, gammas),
        verify_operator_charge("GHZ adjoint |111><000|", 0b111, 0b000, gammas),
    ]

    ghz_charges = [rows[-2][2], rows[-1][2]]
    if ghz_charges != [maximum_charge, maximum_charge]:
        raise AssertionError("GHZ off-diagonal operators must carry the maximal charge")
    if rows[1][2] == maximum_charge:
        raise AssertionError("one-site control must not be mistaken for maximal disagreement")
    if rows[0][2] != 0.0:
        raise AssertionError("population control must remain uncharged")

    verified_sizes = []
    for checked_n in range(2, 6):
        checked_gammas = tuple(0.125 * (site + 1) for site in range(checked_n))
        checked_maximum = 2.0 * sum(checked_gammas)
        top = (1 << checked_n) - 1
        ghz = verify_operator_charge(
            f"GHZ N={checked_n}", 0, top, checked_gammas
        )[2]
        adjoint = verify_operator_charge(
            f"GHZ adjoint N={checked_n}", top, 0, checked_gammas
        )[2]
        population = verify_operator_charge(
            f"population N={checked_n}", 0, 0, checked_gammas
        )[2]
        one_site = verify_operator_charge(
            f"one-site N={checked_n}", 0, 1, checked_gammas
        )[2]
        if [ghz, adjoint] != [checked_maximum, checked_maximum]:
            raise AssertionError(f"N={checked_n}: GHZ charge is not maximal")
        if population != 0.0 or one_site == checked_maximum:
            raise AssertionError(f"N={checked_n}: controls do not separate the maximum")
        verified_sizes.append(checked_n)

    lines = [
        "OPERATOR-LEVEL DEPHASING SUPPORT CHECK",
        "=======================================",
        f"N={n}; local rates gamma={gammas}",
        "",
        "Current verdict:",
        "  F22 concerns computational-basis coherence operators, not state weights.",
        "  |0...0><1...1| and its adjoint differ on every site and therefore",
        "  receive the maximal local Z-dephasing charge 2*sum(gamma_l).",
        "  No GHZ/W protection ranking, optimal state, standing-wave syndrome,",
        "  or information lifetime is inferred here.",
        "",
        f"{'operator':<42} {'distance':>8} {'charge':>10}",
        "-" * 62,
    ]
    for label, distance, charge in rows:
        lines.append(f"{label:<42} {distance:>8d} {charge:>10.6f}")
    lines.extend(
        [
            "",
            "Gates:",
            "  PASS direct Z-conjugation action equals the Hamming-charge formula.",
            "  PASS population control has zero charge.",
            "  PASS one-site coherence control is nonzero but not maximal.",
            "  PASS both GHZ off-diagonal operators have maximal charge.",
            f"  PASS GHZ, adjoint, population, and one-site controls at N={verified_sizes}.",
        ]
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    report = "\n".join(lines) + "\n"
    OUT.write_text(report, encoding="utf-8")
    print(report, end="")
    print(f"Results written to {OUT}")


if __name__ == "__main__":
    main()
