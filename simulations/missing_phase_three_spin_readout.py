"""Exact check for PROOF_MISSING_PHASE_SLOW_READOUT.md, section 7.1.

Uniform XY chain, h hopping 2, centre-only Z-dephasing. This follows the
owned slow-readout quadrature algebra. It does not simulate a hardware run
or extend the finite-defect witness. Run from any directory; writes
results/missing_phase_three_spin_readout.json. Bit j is physical site j.
"""

import json
from itertools import combinations
from pathlib import Path

import sympy as sp


N = 7
DIM = 1 << N
FLIP = DIM - 1
HALF = sp.Rational(1, 2)


def zero(value):
    entries = value.values() if isinstance(value, sp.SparseMatrix) else (
        list(value) if isinstance(value, sp.MatrixBase) else [value])
    return all(sp.simplify(entry) == 0 for entry in entries)


def require(label, condition):
    if not condition:
        raise RuntimeError(label)


def physical(a, b=None):
    """Direct computational-basis assembly of both spin-flipped copies."""
    if b is None:
        b = sp.zeros(N)
    rho = sp.MutableSparseMatrix(DIM, DIM, {})
    for i in range(N):
        for j in range(N):
            left_i, left_j = 1 << i, 1 << j
            right_i, right_j = FLIP ^ left_i, FLIP ^ left_j
            rho[left_i, left_j] = a[i, j] / 2
            rho[right_i, right_j] = a[i, j] / 2
            rho[left_i, right_j] = b[i, j] / 2
            rho[right_i, left_j] = b[i, j] / 2
    return rho


def pauli_action(word, ops):
    target, factor = word, sp.S.One
    for site, letter in ops.items():
        bit = (word >> site) & 1
        if letter == "X":
            target ^= 1 << site
        elif letter == "Y":
            factor *= sp.I * (-1) ** bit
            target ^= 1 << site
        elif letter == "Z":
            factor *= (-1) ** bit
        else:
            raise ValueError(letter)
    return target, factor


def pauli_matrix(ops):
    return sp.MutableSparseMatrix(DIM, DIM, {
        (pauli_action(word, ops)[0], word): pauli_action(word, ops)[1]
        for word in range(DIM)})


def expectation(rho, ops):
    value = 0
    for (row, col), coefficient in rho.todok().items():
        target, factor = pauli_action(row, ops)
        if target == col:
            value += factor * coefficient
    return sp.simplify(value)


def partial_trace(rho, kept):
    mask = sum(1 << site for site in kept)
    out = sp.zeros(1 << len(kept))
    for (row, col), coefficient in rho.todok().items():
        if (row & (FLIP ^ mask)) != (col & (FLIP ^ mask)):
            continue
        r = sum(((row >> site) & 1) << k for k, site in enumerate(kept))
        c = sum(((col >> site) & 1) << k for k, site in enumerate(kept))
        out[r, c] += coefficient
    return out.applyfunc(sp.simplify)


def main():
    eye = sp.eye(N)
    b0, b1, b2 = [(eye[:, j] - eye[:, 6 - j]) / sp.sqrt(2)
                   for j in range(3)]
    v = (b0 - b2) / sp.sqrt(2)
    u = (b0 + b2) / sp.sqrt(2)
    e_plus, e_minus = (u + b1) / sp.sqrt(2), (u - b1) / sp.sqrt(2)
    m = (e_plus * v.T + v * e_minus.T) / (2 * sp.sqrt(2))
    cosine, sine = m + m.H, sp.I * (m.H - m)
    c, s = sp.symbols("c s", real=True)
    eta, eta_plus, eta_minus = sp.symbols("eta eta_plus eta_minus", real=True)
    gamma = sp.symbols("gamma", positive=True)
    omega = 2 * sp.sqrt(2)
    psi = (v + c * u - sp.I * s * b1) / sp.sqrt(2)
    psi_theta = (-s * u - sp.I * c * b1) / sp.sqrt(2)
    h = sp.zeros(N)
    for j in range(N - 1):
        h[j, j + 1] = h[j + 1, j] = 2
    z = sp.eye(N)
    z[3, 3] = -1
    require("exact Schrodinger trajectory", zero(sp.I * omega * psi_theta - h * psi))
    require("centre is empty within the moving copy", zero(z * psi - psi))
    require("initial state", zero(psi.subs({c: 1, s: 0}) - b0))
    require("phase-circle normalization", zero((psi.H * psi)[0] - (1 + c*c + s*s)/2))
    a = psi * psi.H
    rho = physical(a, eta * a)
    da = -sp.I * (h * a - a * h)
    require("A has no dephasing cost", zero(gamma * (z * a * z - a)))
    require("B follows eta=e^(-2 gamma t)",
            zero(-sp.I * (h * (eta*a) - (eta*a) * h)
                 - gamma * (z * (eta*a) * z + eta*a) - eta * (da - 2*gamma*a)))

    q = {0: "Y", 1: "X", 3: "Z"}
    q_other = {0: "X", 1: "Y", 3: "Z"}
    untagged = {0: "Y", 1: "X"}
    residue = expectation(physical(m), q)
    sine_read = expectation(physical(sine), q)
    cosine_read = expectation(physical(cosine), q)
    full_read = expectation(rho, q)
    expected_residue = sp.I / (4 * sp.sqrt(2))
    require("three-spin slow residue", zero(residue - expected_residue))
    require("sine readout", zero(sine_read - 1/(2*sp.sqrt(2))))
    require("cosine null", cosine_read == 0)
    require("full signal has slow and doubled frequency", zero(full_read - s*(1+c)/(2*sp.sqrt(2))))
    require("untagged signal cancels", expectation(rho, untagged) == 0)
    require("opposite current term", zero(expectation(rho, q_other) + full_read))
    require("tag current alone does not read intercopy B", not full_read.has(eta))
    require("removing the tag kills sine", expectation(physical(sine), untagged) == 0)
    require("changing Z tag to X is not a copy reading",
            expectation(physical(sine), {0: "Y", 1: "X", 3: "X"}) == 0)
    single_copy = sp.MutableSparseMatrix(DIM, DIM, {
        (1 << i, 1 << j): sine[i, j] for i in range(N) for j in range(N) if sine[i, j] != 0})
    require("one-copy control exposes current", zero(expectation(single_copy, untagged) - sine_read))

    # Physical states at two positive times of the SAME forward trajectory:
    # theta=pi/2 and theta=3pi/2. B has different positive decay factors there.
    ap = a.subs({c: 0, s: 1})
    am = a.subs({c: 0, s: -1})
    rp, rm = physical(ap, eta_plus * ap), physical(am, eta_minus * am)
    pair_count = 0
    for pair in combinations(range(N), 2):
        require(f"same pair snapshot {pair}", zero(partial_trace(rp-rm, pair)))
        require(f"isolated sine invisible at {pair}", zero(partial_trace(physical(sine), pair)))
        pair_count += 1
    require("physical opposite phases differ at three sites",
            not zero(partial_trace(rp-rm, (0, 1, 3))))
    require("positive-phase signal", zero(expectation(rp, q) - 1/(2*sp.sqrt(2))))
    require("negative-phase signal", zero(expectation(rm, q) + 1/(2*sp.sqrt(2))))
    require("endpoint ZZ equal", expectation(rp, {0: "Z", 6: "Z"}) == HALF
            and expectation(rm, {0: "Z", 6: "Z"}) == HALF)
    for sign, e in ((1, eta_plus), (-1, eta_minus)):
        d = (v - sign*sp.I*b1) / sp.sqrt(2)
        left = sp.zeros(DIM, 1)
        right = sp.zeros(DIM, 1)
        for j in range(N):
            left[1 << j], right[FLIP ^ (1 << j)] = d[j], d[j]
        even, odd = (left+right)/sp.sqrt(2), (left-right)/sp.sqrt(2)
        # Convex density decomposition for 0 <= eta <= 1, no positivity clipping.
        convex = (1+e)/2 * (even*even.H) + (1-e)/2 * (odd*odd.H)
        require("positive physical preparation", zero(physical(d*d.H, e*d*d.H) - convex))
        require("normalization", (d.H*d)[0] == 1)

    # Independent full-spin current convention from the continuity equation.
    h_full = sum((pauli_matrix({j: "X", j+1: "X"})
                  + pauli_matrix({j: "Y", j+1: "Y"}) for j in range(N-1)),
                 sp.MutableSparseMatrix(DIM, DIM, {}))
    n0 = (sp.SparseMatrix.eye(DIM) - pauli_matrix({0: "Z"}))/2
    j01 = -(pauli_matrix({0: "X", 1: "Y"}) - pauli_matrix({0: "Y", 1: "X"}))
    require("outward current is -d n0/dt", zero(j01 + sp.I*(h_full*n0 - n0*h_full)))
    z3 = pauli_matrix({3: "Z"})
    current_tag = z3 * j01
    require("current and single Pauli read have stated normalization",
            zero(sp.trace(current_tag*rp) - 2*expectation(rp, q)))
    require("full Pauli matrix agrees with direct bit action",
            zero(sp.trace(pauli_matrix(q)*rp) - expectation(rp, q)))
    zz = pauli_matrix({0: "Z", 6: "Z"})
    derivative_zz = sp.simplify(sp.trace(sp.I * (h_full*zz - zz*h_full) * rho))
    require("tagged motion equals the endpoint slope on this trajectory",
            zero(derivative_zz - 8*full_read))
    zz_read = expectation(rho, {0: "Z", 6: "Z"})
    chain_rule_zz = omega * (-s*sp.diff(zz_read, c) + c*sp.diff(zz_read, s))
    require("direct commutator agrees with trajectory derivative",
            zero(derivative_zz - chain_rule_zz))

    result = {
        "scope": "N=7 uniform XY, J=1 (h hopping 2), centre-only gamma; internal A motion. Exact symbolic arithmetic.",
        "observable": "Y_0 X_1 Z_3",
        "slow_complex_residue": str(residue),
        "isolated_sine_expectation": str(sine_read),
        "isolated_cosine_expectation": str(cosine_read),
        "full_trajectory_expectation": str(full_read),
        "phase_definition": "c=cos(theta), s=sin(theta), theta=2*sqrt(2)*t",
        "harmonic_form": "sin(theta)/(2*sqrt(2)) + sin(2*theta)/(4*sqrt(2))",
        "untagged_Y0X1": "0",
        "paired_snapshot_times": ["pi/(4*sqrt(2))", "3*pi/(4*sqrt(2))"],
        "all_pair_marginals_equal": pair_count,
        "snapshot_Z0Z6": ["1/2", "1/2"],
        "snapshot_tagged_readouts": [str(expectation(rp, q)), str(expectation(rm, q))],
        "intercopy_coherence_independent": True,
        "minimal_support_for_uniform_sine": 3,
        "controls": ["remove centre tag: zero", "replace tag Z by X: zero", "remove flipped copy: nonzero", "cosine component: zero", "opposite phase: sign reversal"],
        "physical_current_convention": "j_0_to_1=-(X0*Y1-Y0*X1)=-d n0/dt at J=1; <Z3*j>=2*<Y0*X1*Z3> on this trajectory",
        "endpoint_slope_identity": "d<Z0 Z6>/dt = 8 <Y0 X1 Z3> on this uniform trajectory at J=1",
        "fences": ["The full signal includes a second harmonic.", "Minimum support three is for a single-time readout before decoding; a two-spin time series can recover direction from its slope.", "This is internal A motion, not the B intercopy phase with minimum support 5.", "No finite-defect minimal-support theorem, all-N result, or hardware qualification."]
    }
    out = Path(__file__).resolve().parent / "results" / "missing_phase_three_spin_readout.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
