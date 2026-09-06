#!/usr/bin/env python3
"""What a longitudinal field costs the F1 rate pairing, and what decides the cost.

The Heisenberg chain under Z-dephasing satisfies F1 exactly. Adding a
longitudinal field sum_k h_k Z_k breaks the full complex multiset at any
nonzero field, and the question here is the weaker RATE pairing: whether the
decay-rate multiset stays symmetric about Sigma_gamma.

Two things decide it, and neither is the field's strength.

  WHETHER it breaks: commutation. A constant profile is a multiple of the
  conserved total spin sum_k Z_k, commutes with the bond Hamiltonian exactly,
  and cannot move a rate at any amplitude. Every non-constant profile can,
  ordered or not: a disordered all-positive profile breaks the pairing as
  readily as a staggered one, and a global offset is free.

  AT WHAT ORDER: the reflection parity of the deviation from constant. Along an
  R-odd direction (R the chain reversal k -> N-1-k) the response is quadratic;
  a direction with no definite parity responds linearly. The staggered profile
  is R-odd at even N, which is why it alone gives a clean square there.

This is the repository's existing distinction (MIRROR_SYMMETRY_PROOF.md, the
uniform-versus-non-uniform rows, and the OpenArcs arc two_coast_classifier_repair,
which states the commutation mechanism) with the response order added.

Script: simulations/afm_field_palindrome.py
Output: simulations/results/afm_field_palindrome.txt
"""

from pathlib import Path
import math
import numpy as np

RNG = np.random.default_rng(7)

RESULTS_DIR = Path(__file__).parent / "results"
J, GAMMA = 1.0, 0.05
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def tensor(ops):
    result = ops[0]
    for operator in ops[1:]:
        result = np.kron(result, operator)
    return result


def site(pauli, k, n):
    ops = [I2] * n
    ops[k] = pauli
    return tensor(ops)


def liouvillian(n, profile):
    dimension = 2**n
    identity = np.eye(dimension, dtype=complex)
    hamiltonian = bond_hamiltonian(n)
    for k in range(n):
        hamiltonian += profile[k] * site(Z, k, n)
    generator = -1j * (np.kron(hamiltonian, identity)
                       - np.kron(identity, hamiltonian.T))
    for k in range(n):
        z_k = site(Z, k, n)
        generator += GAMMA * (np.kron(z_k, z_k.conj()) - np.eye(dimension**2))
    return generator


def rate_defect(n, profile):
    """Worst mismatch of the DECAY-RATE multiset under r -> 2*sigma - r.

    The rates are sorted and compared with their sorted reflections about
    sigma, which is a genuine multiset distance: every rate must find a
    partner, so no mode can satisfy it by sitting still.

    Only the rate reading is computed here. A defect for the full complex
    multiset would need a bottleneck matching, and a greedy nearest-unused
    pass is not one: it saturates at the diameter 2*sigma and is not even
    monotone in the field, so it reports the matcher rather than the spectrum.
    The repository's canonical object for that is
    F1SpectrumStatistics.MaxF1RatePairingDistance.
    """
    values = np.linalg.eigvals(liouvillian(n, profile))
    sigma = n * GAMMA
    rates = np.sort(-values.real)
    return float(np.max(np.abs(rates - np.sort(2 * sigma - rates))))


def field_operator(n, profile):
    return sum(profile[k] * site(Z, k, n) for k in range(n))


def bond_hamiltonian(n):
    dimension = 2**n
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for bond in range(n - 1):
        for pauli in (X, Y, Z):
            ops = [I2] * n
            ops[bond], ops[bond + 1] = pauli, pauli
            hamiltonian += J * tensor(ops)
    return hamiltonian


def commutator_norm(n, profile):
    """||[H_bond, sum h_k Z_k]||, the quantity that decides the whole question."""
    bonds = bond_hamiltonian(n)
    field = field_operator(n, profile)
    return float(np.linalg.norm(bonds @ field - field @ bonds))


def render_report():
    lines = [
        "A LONGITUDINAL FIELD AGAINST THE F1 RATE PAIRING",
        f"Heisenberg chain, J={J}, Z-dephasing gamma={GAMMA}; field sum_k h_k Z_k.",
        "Defect = worst mismatch of the decay-rate multiset under r -> 2*Sigma_gamma - r.",
        "",
        "The discriminator is commutation, not shape. A profile that commutes with the",
        "bond Hamiltonian cannot move a decay rate; every other profile can, whether or",
        "not it is ordered.",
        "",
        f"{'N':>3} {'profile':>26} {'||[H_bond,F]||':>15} {'rate defect':>13}",
    ]

    def row(n, label, profile):
        lines.append(f"{n:3d} {label:>26} {commutator_norm(n, profile):15.3e} "
                     f"{rate_defect(n, profile):13.3e}")

    families = {}
    for n in (2, 3, 4):
        uniform = [0.3] * n
        staggered = [0.3 * (-1.0) ** k for k in range(n)]
        disordered = list(np.round(RNG.uniform(0.2, 0.4, n), 4))
        row(n, "uniform (0.3, ...)", uniform)
        row(n, "staggered (+/-0.3)", staggered)
        row(n, "disordered, all positive", disordered)
        lines.append("")
        families[n] = (uniform, staggered, disordered)

    lines += [
        "The commutator column is the mechanism, read directly. At amplitude 0.3 it is a",
        "few ulp rather than 0.0 because 0.3 is not dyadic; the exact assertion below runs",
        "at a dyadic amplitude, where it is 0.0 to the bit.",
        "",
        "The disordered profile carries no order at all and breaks the pairing anyway,",
        "as strongly as the staggered one. What the uniform profile has is not order but",
        "a vanishing commutator: sum_k Z_k is the conserved total spin of the Heisenberg",
        "bond, so its superoperator is diagonal with purely imaginary spectrum and commutes",
        "with the rest of the generator. It can only move frequencies. A global offset is",
        "therefore free, and the defect depends on the profile's deviation from constant.",
        "",
        "Global offset, added to a fixed deviation:",
        f"{'N':>3} {'offset':>10} {'rate defect':>15}",
    ]
    offset_values = {}
    for n in (2, 3, 4):
        base = [0.001 * (-1.0) ** k for k in range(n)]
        defects = []
        for offset in (0.0, 0.5, 3.0):
            defect = rate_defect(n, [b + offset for b in base])
            defects.append(defect)
            lines.append(f"{n:3d} {offset:10.1f} {defect:15.6e}")
        offset_values[n] = defects
        lines.append("")

    lines += [
        "The response ORDER is set by the reflection parity of the deviation, not by its",
        "amplitude. R is the chain reversal k -> N-1-k. Exponent is read from the defect",
        "at two field scales a decade apart, so a wrong exponent cannot hide in a constant.",
        "",
        f"{'N':>3} {'deviation':>22} {'R-parity':>9} {'defect(1e-3)':>14} {'defect(1e-4)':>14} {'exponent':>9}",
    ]
    directions = [
        (4, "staggered (+,-,+,-)", [1.0, -1.0, 1.0, -1.0]),
        (4, "(+2,-1,+1,-2)", [2.0, -1.0, 1.0, -2.0]),
        (4, "(+1,+1,-1,-1)", [1.0, 1.0, -1.0, -1.0]),
        (4, "single site (0,0,0,1)", [0.0, 0.0, 0.0, 1.0]),
    ]
    odd_exponents, unparity_exponents = [], []
    for n, label, raw in directions:
        vector = np.array(raw, dtype=float)
        vector = vector - vector.mean()
        reversed_vector = vector[::-1]
        if np.allclose(reversed_vector, -vector):
            parity = "odd"
        elif np.allclose(reversed_vector, vector):
            parity = "even"
        else:
            parity = "none"
        big = rate_defect(n, list(1e-3 * vector))
        small = rate_defect(n, list(1e-4 * vector))
        exponent = math.log10(big / small)
        lines.append(f"{n:3d} {label:>22} {parity:>9} {big:14.3e} {small:14.3e} {exponent:9.2f}")
        if parity == "odd":
            odd_exponents.append((label, exponent))
        elif parity == "none":
            unparity_exponents.append((label, exponent))
    lines.append("")

    # Gates, each able to fail.
    # (1) The theorem's hypothesis, exact and eigensolver-free: a constant profile
    #     commutes with the bond Hamiltonian. This is the reason the rates hold,
    #     so it is checked exactly rather than inferred from the rate column.
    for n in (2, 3, 4):
        uniform, staggered, disordered = families[n]
        exact_commutator = commutator_norm(n, uniform)
        if exact_commutator > 1e-12:
            raise RuntimeError(
                f"a constant profile failed to commute at N={n}: {exact_commutator:.3e}")
        if rate_defect(n, uniform) > 1e-12:
            raise RuntimeError(
                f"a commuting profile moved the rates at N={n}: {rate_defect(n, uniform):.3e}")
        # (2) Order is not the discriminator: a disordered profile must break too.
        if rate_defect(n, disordered) < 1e-4:
            raise RuntimeError(
                f"a disordered non-constant profile did not break the pairing at N={n}")
    # (3) A global offset is inert. The exact route is the reason, not the
    #     consequence: the offset part of any profile is a multiple of sum_k Z_k,
    #     which commutes with the bond Hamiltonian exactly, so it cannot reach a
    #     decay rate. That is checked at 0.0; the eigensolver's own spread across
    #     offsets is then reported rather than gated, and it is relative to a
    #     defect four decades larger.
    for n in (2, 3, 4):
        offset_commutator = commutator_norm(n, [3.0] * n)
        if offset_commutator != 0.0:
            raise RuntimeError(
                f"the offset part does not commute exactly at N={n}: {offset_commutator:.3e}")
    for n, defects in offset_values.items():
        spread = max(defects) - min(defects)
        relative = spread / max(defects)
        lines.append(f"  N={n}: offsets agree to {relative:.1e} relative "
                     f"(absolute spread {spread:.1e}, eigensolver noise on an exact identity)")
        if relative > 1e-6:
            raise RuntimeError(
                f"a global offset moved the defect at N={n} beyond the solver's reach: "
                f"{relative:.3e} relative")
    lines.append("")
    # (4) Every R-odd direction responds at order 2; a direction with no parity does not.
    for label, exponent in odd_exponents:
        if abs(exponent - 2.0) > 0.02:
            raise RuntimeError(f"R-odd direction {label} responded at order {exponent:.2f}, not 2")
    for label, exponent in unparity_exponents:
        if abs(exponent - 2.0) <= 0.02:
            raise RuntimeError(
                f"the parity-free direction {label} responded at order 2; the parity law "
                "no longer separates")

    lines += [
        "Every R-odd deviation responds at order exactly 2, and the single-site bump, which",
        "is neither R-odd nor R-even, responds at order 1. That is F131's order-sorting law",
        "read as a response order: along an R-odd direction the spectrum is even in the",
        "parameter, so the linear term is absent and the leading response is quadratic.",
        "The staggered profile is R-odd at even N, which is the whole reason it gave a clean",
        "square; at odd N it is not, and no clean power should be expected there.",
        "",
        "So the field costs the F1 rate pairing exactly its non-commuting part, at an order",
        "the deviation's reflection parity decides. Order in the antiferromagnetic sense is",
        "neither necessary nor sufficient: a disordered profile breaks the pairing, and a",
        "perfectly ordered uniform one does not.",
        "",
        "Scope: finite chains N=2..4, one coupling, one dephasing rate, and the RATE",
        "multiset only. The full complex F1 multiset breaks under every non-zero field,",
        "uniform included, and is a separate statement not measured here.",
    ]
    return "\n".join(lines) + "\n"


def main():
    output = render_report()
    print(output, end="")
    RESULTS_DIR.mkdir(exist_ok=True)
    (RESULTS_DIR / "afm_field_palindrome.txt").write_text(output, encoding="utf-8")
    print("Results: simulations/results/afm_field_palindrome.txt")


if __name__ == "__main__":
    main()
