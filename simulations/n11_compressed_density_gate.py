"""Exact N=11 witness for the (1,1) compressed-density closure.

For the uniform open XY chain H = J sum(XX+YY), J > 0, the one-excitation
energies are 4J cos(k*pi/12). This gate enumerates the complete ad_H frequency
space of |psi_1><psi_6| and uses the *physical cell indicator*

    N_l |a><b| = [a == l XOR b == l] |a><b|

to compute its compressed site densities and dissipator. All comparisons are
exact SymPy equalities; no eigensolver or frequency grouping tolerance enters.
The interval claim belongs to the compression, not to finite-J Liouvillian
spectra. The zero-frequency room is a scalar-parity room, where F154's identity
is a theorem; the gate reads its two ends, 0 on the one-excitation identity and
-4*gbar on a chiral difference of mode projectors, off the computed action of
the physical cell dissipator. Both halves run through the same injectable cell
rule, so a wrong cell action fails the endpoint half on its own.

Run: python simulations/n11_compressed_density_gate.py
"""

import sympy as sp


def sine_modes(n):
    """Normalized open-chain sine modes, one-indexed in k and zero-indexed in site."""
    m = n + 1
    return {
        k: tuple(sp.sqrt(sp.Rational(2, m)) * sp.sin(sp.pi * k * (z + 1) / m)
                 for z in range(n))
        for k in range(1, n + 1)
    }


def frequency_dyads(n, reference):
    """Every (k,b) with cos(k*pi/M)-cos(b*pi/M) equal to the reference."""
    m = n + 1
    energies = {k: sp.cos(sp.pi * k / m) for k in range(1, n + 1)}
    target = energies[reference[0]] - energies[reference[1]]
    return tuple(
        (k, b) for k in range(1, n + 1) for b in range(1, n + 1)
        if sp.simplify(energies[k] - energies[b] - target) == 0
    )


def reflection_parity(dyad):
    """Parity of |psi_k><psi_b| under site reflection."""
    k, b = dyad
    return (-1) ** (k + b)


def physical_disagreement(ket_site, bra_site, site):
    """Local Z-dephasing disagreement on the physical |ket><bra| cell."""
    return int((ket_site == site) != (bra_site == site))


def _cell_amplitude(modes, dyad, ket_site, bra_site):
    k, b = dyad
    return modes[k][ket_site] * modes[b][bra_site]


def _compressed_cell_diagonal(modes, dyads, cell_value):
    """Direct V.T diag(cell_value) V on all n*n physical (1,1) cells."""
    n = len(next(iter(modes.values())))
    rows = []
    for outgoing in dyads:
        row = []
        for incoming in dyads:
            value = sum(
                cell_value(a, b)
                * _cell_amplitude(modes, outgoing, a, b)
                * _cell_amplitude(modes, incoming, a, b)
                for a in range(n) for b in range(n)
            )
            row.append(sp.simplify(value))
        rows.append(row)
    return sp.Matrix(rows)


def compressed_indicator(modes, dyads, site, cell_rule=None):
    """Compression of N_l; an injected cell rule follows exactly this path."""
    n = len(next(iter(modes.values())))
    assert 0 <= site < n
    rule = physical_disagreement if cell_rule is None else cell_rule
    return _compressed_cell_diagonal(modes, dyads, lambda a, b: rule(a, b, site))


def compressed_dissipator(modes, dyads, rates, cell_rule=None):
    """Direct physical D=-2 sum_l gamma_l N_l, before compression."""
    n = len(next(iter(modes.values())))
    assert len(rates) == n
    rule = physical_disagreement if cell_rule is None else cell_rule
    return _compressed_cell_diagonal(
        modes, dyads,
        lambda a, b: -2 * sum(rates[site] * rule(a, b, site) for site in range(n))
    )


def compressed_diagonal_overlap(modes, dyads, rates):
    """T=sum_l gamma_l |ll><ll|, compressed independently of N_l."""
    n = len(next(iter(modes.values())))
    assert len(rates) == n
    return _compressed_cell_diagonal(
        modes, dyads, lambda a, b: rates[a] if a == b else 0
    )


def compressed_gamma_legs(modes, dyads, rates):
    """Compression of Gamma on the ket and bra legs, from physical cells."""
    n = len(next(iter(modes.values())))
    assert len(rates) == n
    return _compressed_cell_diagonal(modes, dyads, lambda a, b: rates[a] + rates[b])


def compressed_dissipator_action(modes, dyads, rates, coefficients, cell_rule=None):
    """Project physical D acting on a single superposition of mode dyads."""
    n = len(next(iter(modes.values())))
    assert len(rates) == n and len(coefficients) == len(dyads)
    rule = physical_disagreement if cell_rule is None else cell_rule
    active = [(coefficient, dyad) for coefficient, dyad in zip(coefficients, dyads)
              if coefficient != 0]
    state = {(a, b): sum(c * _cell_amplitude(modes, d, a, b) for c, d in active)
             for a in range(n) for b in range(n)}
    return sp.Matrix([
        sp.simplify(sum(
            -2 * sum(rates[site] * rule(a, b, site)
                     for site in range(n))
            * _cell_amplitude(modes, outgoing, a, b) * state[a, b]
            for a in range(n) for b in range(n)
        ))
        for outgoing in dyads
    ])


def matrix_equal(left, right):
    if left.shape != right.shape:
        return False
    return all(sp.simplify(a - b) == 0 for a, b in zip(left, right))


def rayleigh_on_eigenvector(action, vector):
    """The eigenvalue read off a computed action, after checking it is one."""
    pivot = next(i for i, entry in enumerate(vector) if entry != 0)
    value = sp.simplify(action[pivot] / vector[pivot])
    assert matrix_equal(action, value * vector), "not an eigenvector of the compression"
    return value


def verify_zero_frequency_room(cell_rule=None):
    """The block-wide endpoints, computed from physical cells through cell_rule."""
    n = 11
    modes = sine_modes(n)
    balanced = (2,) + (1,) * 9 + (0,)
    zero_dyads = frequency_dyads(n, (1, 1))
    assert zero_dyads == tuple((k, k) for k in range(1, n + 1))
    assert all(
        sp.simplify(sum(modes[k][a] * modes[k][b] for k in range(1, n + 1))
                    - int(a == b)) == 0
        for a in range(n) for b in range(n)
    ), "zero-frequency identity is not physical-cell diagonal"
    identity = sp.ones(n, 1)
    upper = rayleigh_on_eigenvector(
        compressed_dissipator_action(modes, zero_dyads, balanced, identity, cell_rule),
        identity)
    assert upper == 0, "zero-frequency identity misses the upper endpoint"
    chiral = sp.Matrix([1] + [0] * (n - 2) + [-1])
    assert all(sp.simplify(modes[1][a] ** 2 - modes[n][a] ** 2) == 0
               for a in range(n))
    lower = rayleigh_on_eigenvector(
        compressed_dissipator_action(modes, zero_dyads, balanced, chiral, cell_rule),
        chiral)
    assert lower == -4, "zero-frequency chiral difference misses the lower endpoint"
    return {"zero_frequency_upper": upper, "zero_frequency_lower": lower}


def verify_n11(cell_rule=None):
    """Run the literal-value checks; wrong physical cell actions raise here."""
    n = 11
    modes = sine_modes(n)
    dyads = frequency_dyads(n, (1, 6))
    expected_dyads = ((1, 6), (3, 7), (5, 9), (6, 11))
    assert dyads == expected_dyads, f"incomplete frequency space: {dyads}"
    assert tuple(reflection_parity(d) for d in dyads) == (-1, 1, 1, -1)

    left = compressed_indicator(modes, dyads, 0, cell_rule)
    right = compressed_indicator(modes, dyads, 10, cell_rule)
    contrast = left - right
    expected_contrast = -sp.sqrt(2) / 72 * sp.Matrix(
        [[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 1], [0, 1, 1, 0]]
    )
    assert matrix_equal(contrast, expected_contrast), "contrast matrix differs from literal"
    assert left[2, 0] == -sp.sqrt(2) / 144  # opposite-parity positive control
    assert right[2, 0] == sp.sqrt(2) / 144
    assert left[2, 1] == -sp.Rational(1, 72) - sp.sqrt(3) / 144
    assert contrast[2, 1] == 0  # same parity, but N_0 itself is nonzero
    assert contrast.rank() == 2
    assert sp.trace(contrast * contrast) == sp.Rational(1, 324)
    expected_eigenvalues = {-sp.sqrt(2) / 36: 1, sp.Integer(0): 2, sp.sqrt(2) / 36: 1}
    assert contrast.eigenvals() == expected_eigenvalues
    # The one-legged mirror contrast vanishes on a complete fixed-frequency
    # room. The physical C_0 is the signed agreement/double-occupancy term.
    signed_ends = (1,) + (0,) * 9 + (-1,)
    assert matrix_equal(
        contrast, -2 * compressed_diagonal_overlap(modes, dyads, signed_ends)
    )

    balanced = (2,) + (1,) * 9 + (0,)
    diss = compressed_dissipator(modes, dyads, balanced, cell_rule)
    gamma_legs = compressed_gamma_legs(modes, dyads, balanced)
    diagonal_overlap = compressed_diagonal_overlap(modes, dyads, balanced)
    assert matrix_equal(gamma_legs, 2 * sp.eye(4)), "balanced Gamma compression"
    assert matrix_equal(diss, -4 * sp.eye(4) + 4 * diagonal_overlap), "diagonal-overlap identity"
    assert diss[2, 0] == sp.sqrt(2) / 36, "balanced F154 identity failure"
    # The old conditional identity would give -2*N_XY and a zero cross entry.
    size_cross = sum(compressed_indicator(modes, dyads, site, cell_rule)[2, 0]
                     for site in range(n))
    assert size_cross == 0
    expected_room_spectrum = {
        sp.Integer(-4): 2,
        -sp.Rational(10, 3) - sp.sqrt(2) / 18: 1,
        -sp.Rational(10, 3) + sp.sqrt(2) / 18: 1,
    }
    assert diss.eigenvals() == expected_room_spectrum
    contact = sp.Matrix([[modes[k][z] * modes[b][z] for k, b in dyads]
                         for z in range(n)])
    assert contact[:, 0] == contact[:, 3]
    assert contact[:, 1] == contact[:, 2]
    assert contact.rank() == 2

    # At omega=0, I is physical-cell diagonal, so D I=0. The chiral
    # difference Q=P_1-P_11 has no physical diagonal and its direct
    # compressed physical-D action is -4Q for this balanced profile.
    endpoints = verify_zero_frequency_room(cell_rule)

    off_locus = (1,) + (0,) * 10
    off_diss = compressed_dissipator(modes, dyads, off_locus, cell_rule)
    off_gamma_legs = compressed_gamma_legs(modes, dyads, off_locus)
    assert off_gamma_legs[2, 2] == (4 + sp.sqrt(3)) / 24
    assert off_gamma_legs[2, 2] != sp.Rational(2, 11), "off-locus Gamma contrast"
    assert off_diss[2, 2] == -(22 + 5 * sp.sqrt(3)) / 72
    assert sp.simplify(off_diss[2, 2] + sp.Rational(4, 11)) < 0
    return {
        "dyads": dyads,
        "contrast": contrast,
        "balanced_cross": diss[2, 0],
        "nonzero_room_spectrum": diss.eigenvals(),
        "zero_frequency_upper": endpoints["zero_frequency_upper"],
        "zero_frequency_lower": endpoints["zero_frequency_lower"],
        "off_locus_gamma_legs": off_gamma_legs[2, 2],
        "off_locus_rayleigh": off_diss[2, 2],
    }


if __name__ == "__main__":
    result = verify_n11()
    print("N=11 (1,1) complete frequency space:", result["dyads"])
    print("C0 =", result["contrast"])
    print("spec(C0) =", result["contrast"].eigenvals())
    print("balanced physical D cross =", result["balanced_cross"])
    print("balanced nonzero-room spec(D) =", result["nonzero_room_spectrum"])
    print("zero-frequency block endpoints =",
          (result["zero_frequency_lower"], result["zero_frequency_upper"]))
    print("off-locus Gamma legs =", result["off_locus_gamma_legs"])
    print("off-locus Rayleigh =", result["off_locus_rayleigh"])
    print("VERDICT: PASS (exact N=11 compression witness; finite-J spectrum not tested)")
