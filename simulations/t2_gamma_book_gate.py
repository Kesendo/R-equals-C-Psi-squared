"""Gate for the T2 -> gamma conversion book (the GLOSSARY anchor).

The repository carries several forms of the conversion from a hardware
coherence time T2 into a Lindblad Z-dephasing rate gamma. This gate decides
which of them are right, and under which model, by measurement rather than by
citation. Every number the GLOSSARY section "The T2 -> gamma conversion"
states is produced here.

The three forms under test, on a chain whose only dissipator is Z-dephasing
(the framework's D[Z]) unless stated otherwise:

    A   gamma = 1/T2                     the factor-2 trap
    B   gamma = 1/(2*T2)                 exact for a dephasing-only model
    C   gamma = (1/T2 - 1/(2*T1))/2      the pure-dephasing part, exact for a
                                         model that ALSO carries sigma- at 1/T1

What this gate is sensitive to: it compares each form's predicted coherence
against the measured decay envelope exp(-t/T2) directly, so a wrong factor of
two shows up as a factor of two in the residual, not as a threshold miss. Every
machine-zero assertion is paired with a read on the same code path whose right
answer is nonzero, so a broken computation cannot pass by returning zeros.

What it cannot see, stated because a review found it the hard way: the
coherence envelope CANNOT distinguish the two models. Any jump operator whose
c-dagger-c is a population projector contributes the same -Gamma/2 to rho_01,
so replacing sigma- by sigma+ (amplitude heating, the opposite population flow)
leaves check 2 identical to the last bit. Section 5 measures that blindness
rather than leaving it implicit, and section 6 reads the populations, where the
two models are not alike at all. A green section 2 is evidence about ONE matrix
element and nothing more.

It also says nothing about which form any given producer script uses. That is a
separate sweep.
"""

import csv
import os
import sys

import numpy as np
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)

# Real IBM Torino calibration pairs (T1, T2) in microseconds, from the two
# chains of experiments/CHAIN_SELECTION_TEST.md, plus two synthetic ends:
# a T1-unlimited qubit and a strongly T1-limited one.
CALIBRATION_PAIRS = [
    ("Q80", 103.3, 27.7),
    ("Q8", 204.1, 218.7),
    ("Q79", 77.5, 91.9),
    ("Q53", 22.4, 62.4),   # T2 > 2*T1: form C is undefined here, form B is not
    ("Q85", 2.9, 5.0),
    ("Q18", 267.6, 211.3),
    ("Q90", 267.0, 275.4),
    ("T1-unlimited", 1.0e6, 50.0),
    ("T1-limited", 30.0, 55.0),
]

TORINO_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "data", "ibm_history", "ibm_torino_history.csv",
)

_failures = []


def check(label, passed, detail):
    status = "PASS" if passed else "FAIL"
    print(f"  [{status}] {label}: {detail}")
    if not passed:
        _failures.append(label)


def dissipator(c, rate):
    """Lindblad dissipator rate * D[c] in the row-stack (order='C') vec
    convention the framework uses: vec(A rho B) = (A (x) B^T) vec(rho)."""
    cdc = c.conj().T @ c
    return rate * (
        np.kron(c, c.conj())
        - 0.5 * np.kron(cdc, I2)
        - 0.5 * np.kron(I2, cdc.T)
    )


def coherence(L, t):
    """|rho_01|(t) starting from the equal superposition."""
    v = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex).flatten()
    return abs((expm(L * t) @ v).reshape(2, 2)[0, 1])


def gamma_A(t1, t2):
    return 1.0 / t2


def gamma_B(t1, t2):
    return 1.0 / (2.0 * t2)


def gamma_C(t1, t2):
    return (1.0 / t2 - 1.0 / (2.0 * t1)) / 2.0


def check_dephasing_decay_rate():
    """D[Z] at rate gamma decays coherences at 2*gamma, not gamma.

    This is the fact the whole book rests on; measure it, do not cite it.
    """
    print("\n1. The D[Z] decay rate (the origin of the factor 2)")
    worst = 0.0
    for gamma in (0.05, 0.3, 1.7):
        L = dissipator(Z, gamma)
        for t in (0.5, 2.0, 7.0):
            worst = max(worst, abs(coherence(L, t) - 0.5 * np.exp(-2.0 * gamma * t)))
    # Exact route, so compare exactly: dissipator(Z, gamma) is diagonal
    # (0, -2g, -2g, 0), expm is element-wise exp, and the reference is the same
    # exp of the same product. A nonzero here would be a finding about the
    # construction, not a tolerance to widen.
    check(
        "coherence decays as exp(-2*gamma*t)",
        worst == 0.0,
        f"worst |deviation| = {worst:.3e} over gamma in (0.05, 0.3, 1.7), "
        "compared exactly (the dissipator is diagonal here)",
    )
    # Paired nonzero read: the ONE-gamma envelope must NOT fit, and must miss
    # by a definite amount. A broken expm returning zeros would fail this.
    L = dissipator(Z, 0.3)
    miss = abs(coherence(L, 2.0) - 0.5 * np.exp(-0.3 * 2.0))
    check(
        "the exp(-gamma*t) envelope is wrong by a measurable amount",
        miss > 0.1,
        f"|deviation| = {miss:.6f} at gamma = 0.3, t = 2 (a zero here would mean "
        "the check cannot fail)",
    )


def check_both_books_reproduce_T2():
    """Forms B and C both reproduce the measured envelope exp(-t/T2) exactly.

    They are not canonical-versus-shorthand. They are the correct D[Z] rate
    for two different models: B when the model carries dephasing alone, C when
    it also carries amplitude damping at 1/T1. Form A reproduces neither.
    """
    print("\n2. Forms B and C against the measured envelope exp(-t/T2)")
    worst_b = worst_c = 0.0
    worst_a = 0.0
    tested_c = 0
    for name, t1, t2 in CALIBRATION_PAIRS:
        grid = np.linspace(1.0, 3.0 * t2, 25)
        L_b = dissipator(Z, gamma_B(t1, t2))
        L_a = dissipator(Z, gamma_A(t1, t2))
        run_c = t2 < 2.0 * t1
        if run_c:
            L_c = dissipator(Z, gamma_C(t1, t2)) + dissipator(SIGMA_MINUS, 1.0 / t1)
            tested_c += 1
        for t in grid:
            ref = 0.5 * np.exp(-t / t2)
            worst_b = max(worst_b, abs(coherence(L_b, t) - ref))
            worst_a = max(worst_a, abs(coherence(L_a, t) - ref))
            if run_c:
                worst_c = max(worst_c, abs(coherence(L_c, t) - ref))
    check(
        "form B (dephasing-only model) is exact",
        worst_b < 1e-14,
        f"worst |deviation| = {worst_b:.3e} over {len(CALIBRATION_PAIRS)} (T1, T2) pairs",
    )
    check(
        "form C (D[Z] + sigma- at 1/T1) is exact where it is defined",
        worst_c < 1e-14,
        f"worst |deviation| = {worst_c:.3e} over the {tested_c} pairs with T2 < 2*T1",
    )
    check(
        "form C's domain exclusion actually fires on this pair list",
        tested_c < len(CALIBRATION_PAIRS),
        f"{len(CALIBRATION_PAIRS) - tested_c} of {len(CALIBRATION_PAIRS)} pairs are "
        "outside form C's domain (a zero here would mean the exclusion is untested)",
    )
    check(
        "form A misses the envelope (the paired nonzero read)",
        worst_a > 0.1,
        f"worst |deviation| = {worst_a:.6f}",
    )


def check_form_A_is_exactly_double():
    """The A-form error is a clean factor of two in the rate, at every T2.

    Stated, not gated: gamma_A / gamma_B is 2.0 for every positive T2 because
    division by two is exact in binary floating point. Asserting it would be a
    check that cannot fail, so it is printed as an identity and the gating is
    left to check 2, where form A is measured against the envelope.
    """
    print("\n3. The size of the A-form error")
    ratios = {gamma_A(t1, t2) / gamma_B(t1, t2) for _, t1, t2 in CALIBRATION_PAIRS}
    print(
        f"      identity (not a check): gamma_A / gamma_B = {ratios.pop():.1f} "
        f"for all {len(CALIBRATION_PAIRS)} pairs, exactly, by construction"
    )


def check_form_C_domain():
    """Form C is undefined where T2 > 2*T1, and that is not a corner case.

    The bound T2 <= 2*T1 is physical, but reported calibrations cross it, and
    one such qubit sits in a live experiment table (Q53 of
    experiments/CHAIN_SELECTION_TEST.md, T1 = 22.4, T2 = 62.4).
    """
    print("\n4. The domain of form C on real calibration data")
    t1_q53, t2_q53 = 22.4, 62.4
    g = gamma_C(t1_q53, t2_q53)
    check(
        "form C returns a negative rate for Q53 of CHAIN_SELECTION_TEST",
        g < 0.0,
        f"gamma_C(T1 = {t1_q53}, T2 = {t2_q53}) = {g:+.6f} per us",
    )
    if not os.path.exists(TORINO_CSV):
        check("Torino history CSV present", False, f"not found at {TORINO_CSV}")
        return
    total = crossing = unparsed = seen = 0
    per_qubit_t2 = {}
    per_qubit_t1 = {}
    q53_crossings = 0
    with open(TORINO_CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            seen += 1
            try:
                t1 = float(row["T1_us"])
                t2 = float(row["T2_us"])
                qubit = int(row["qubit"])
            except (TypeError, ValueError, KeyError):
                unparsed += 1
                continue
            if t1 <= 0.0 or t2 <= 0.0:
                unparsed += 1
                continue
            total += 1
            per_qubit_t2.setdefault(qubit, set()).add(round(t2, 4))
            per_qubit_t1.setdefault(qubit, []).append(t1)
            if t2 > 2.0 * t1:
                crossing += 1
                if qubit == 53:
                    q53_crossings += 1
    pct = 100.0 * crossing / total
    check(
        "T2 > 2*T1 records exist in the Torino history",
        crossing > 0,
        f"{crossing} of {total} records ({pct:.1f}%) give a negative gamma under "
        f"form C; {unparsed} rows carried no usable (T1, T2) pair",
    )
    # The crossings are not spread evenly, and reading the 1.4% as a physical
    # regime would be wrong: a frozen calibration field supplies a large share.
    frozen = sorted(
        q for q, vals in per_qubit_t2.items()
        if len(vals) == 1 and len(per_qubit_t1[q]) > 10
    )
    check(
        "the crossings concentrate on one frozen calibration record",
        frozen == [53] and q53_crossings > 0.3 * crossing,
        f"qubits with a single distinct T2 over >10 dates: {frozen}; Q53 holds "
        f"T2 = {sorted(per_qubit_t2[53])[0]} across {len(per_qubit_t1[53])} dates "
        f"while its T1 ranges over {min(per_qubit_t1[53]):.1f} to "
        f"{max(per_qubit_t1[53]):.1f} us, and supplies {q53_crossings} of the "
        f"{crossing} crossings",
    )
    print(
        f"      residue: {seen} rows read, {total} carried a usable (T1, T2), "
        f"{unparsed} skipped as unparsable or non-positive"
    )


def check_the_coherence_read_is_blind():
    """Section 2's observable cannot tell the two models apart. Measure that.

    Swapping sigma- for sigma+ reverses the population flow and the steady
    state, and leaves the coherence identical, because both jump operators have
    c-dagger-c equal to a population projector and neither transfers 0<->1
    coherence. So a green section 2 is not evidence that the model is right.
    """
    print("\n5. What the coherence read is blind to (the gate's own scope)")
    worst = 0.0
    for name, t1, t2 in CALIBRATION_PAIRS:
        if t2 >= 2.0 * t1:
            continue
        gz = gamma_C(t1, t2)
        L_down = dissipator(Z, gz) + dissipator(SIGMA_MINUS, 1.0 / t1)
        L_up = dissipator(Z, gz) + dissipator(SIGMA_MINUS.conj().T, 1.0 / t1)
        for t in np.linspace(1.0, 3.0 * t2, 15):
            worst = max(worst, abs(coherence(L_down, t) - coherence(L_up, t)))
    check(
        "relaxation and heating are indistinguishable in the coherence",
        worst == 0.0,
        f"worst |coherence(sigma-) - coherence(sigma+)| = {worst:.3e}, exactly "
        "(so section 2 constrains one matrix element, not the model)",
    )


def check_the_models_differ_in_populations():
    """Where the two rows part company, and by how much.

    This is the scope statement the GLOSSARY table needs: the rows agree on the
    coherence envelope and on nothing else. Row 1 has no T1 channel at all.
    """
    print("\n6. Where the two rows differ (populations, and the rate itself)")
    name, t1, t2 = "Q80", 103.3, 27.7
    L_b = dissipator(Z, gamma_B(t1, t2))
    L_c = dissipator(Z, gamma_C(t1, t2)) + dissipator(SIGMA_MINUS, 1.0 / t1)
    v = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex).flatten()
    rho_b = (expm(L_b * 200.0) @ v).reshape(2, 2)
    rho_c = (expm(L_c * 200.0) @ v).reshape(2, 2)
    check(
        f"the two rows separate in the populations on {name}",
        abs(rho_b[0, 0].real - rho_c[0, 0].real) > 0.1,
        f"rho_00 at t = 200 us: row 1 = {rho_b[0, 0].real:.6f} (no T1 channel, "
        f"stays mixed), row 2 = {rho_c[0, 0].real:.6f}",
    )
    # Report the real calibration pairs separately from the two synthetic ends,
    # because the GLOSSARY quotes a real qubit and a reader comparing the two
    # must not meet a number produced by an invented (T1, T2).
    real = []
    for q, a, b in CALIBRATION_PAIRS:
        if b >= 2.0 * a or q.startswith("T1-"):
            continue
        real.append((gamma_B(a, b) / gamma_C(a, b), q))
    worst_ratio, worst_q = max(real)
    check(
        "the rows disagree on gamma itself by more than the factor-2 trap",
        worst_ratio > 2.0,
        f"largest gamma_B / gamma_C on a real Torino pair = {worst_ratio:.2f} "
        f"on {worst_q}; spread over the {len(real)} real pairs is "
        f"{min(s for s, _ in real):.2f} to {worst_ratio:.2f} (the synthetic "
        "T1-limited end reaches 12.00 and is excluded from this read)",
    )


def main():
    print("T2 -> gamma conversion book gate")
    print("=" * 70)
    check_dephasing_decay_rate()
    check_both_books_reproduce_T2()
    check_form_A_is_exactly_double()
    check_form_C_domain()
    check_the_coherence_read_is_blind()
    check_the_models_differ_in_populations()
    print("\n" + "=" * 70)
    if _failures:
        print(f"FAILED ({len(_failures)}): " + "; ".join(_failures))
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
