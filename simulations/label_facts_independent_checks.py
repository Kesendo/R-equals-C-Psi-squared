"""From-below checks that own rules in simulations/tests/test_label_facts.py.

Definitions follow the repository's own: C = purity Tr(rho^2), Psi = l1 coherence
(sum of |off-diagonal| in the computational basis) / (d - 1), CPsi = C * Psi.

Each check_* function returns the values its CHECK line reports. Where an exact
route exists (rationals, surds, integers, an exact matrix identity) the value is
computed in exact arithmetic with sympy; the numerical checks name their error
model where they print a float. A CHECK id is a stable name that the rules cite, so
the ids are not renumbered when a check is removed.

    python simulations/label_facts_independent_checks.py          one line per CHECK
    python simulations/label_facts_independent_checks.py --json   the values that
        simulations/tests/test_label_facts.py asserts

Read-only: the script writes nothing. The repository files it reads are named in
the check that reads them (CHECK-10, CHECK-12).
"""
from __future__ import annotations

import json
import sys
from math import comb
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]

# ---------------------------------------------------------------- exact helpers
S_I2 = sp.eye(2)
S_X = sp.Matrix([[0, 1], [1, 0]])
S_Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
S_Z = sp.diag(1, -1)


def s_kron(*factors):
    out = factors[0]
    for factor in factors[1:]:
        out = sp.kronecker_product(out, factor)
    return out


def s_ket(bits):
    vector = sp.zeros(2 ** len(bits), 1)
    vector[int("".join(map(str, bits)), 2)] = 1
    return vector


def s_outer(ket, bra=None):
    bra = ket if bra is None else bra
    return ket * bra.conjugate().T


def s_l1(rho):
    return sp.simplify(sum(sp.Abs(rho[i, j]) for i in range(rho.rows)
                           for j in range(rho.cols) if i != j))


def s_cpsi(rho):
    return sp.simplify(sp.trace(rho * rho) * s_l1(rho) / (rho.rows - 1))


def s_dissipator(jump, rho):
    adjoint = jump.conjugate().T
    return sp.simplify(jump * rho * adjoint - (adjoint * jump * rho + rho * adjoint * jump) / 2)


def s_rhs(rho, hamiltonian, jumps):
    out = -sp.I * (hamiltonian * rho - rho * hamiltonian)
    for jump in jumps:
        out += s_dissipator(jump, rho)
    return sp.simplify(out)


def s_bell_plus_concurrence():
    """Wootters concurrence of the dephased Bell+ state, exact in f (0 < f < 1).

    rho(f) has 1/2 on the |00>,|11> diagonal and f/2 on the |00><11| coherence.
    The spin-flipped product rho * (YY rho* YY) has eigenvalues ((1+f)/2)^2,
    ((1-f)/2)^2, 0, 0, so the concurrence is (1+f)/2 - (1-f)/2 = f."""
    f = sp.symbols("f", positive=True)
    rho = sp.Matrix([[sp.Rational(1, 2), 0, 0, f / 2], [0, 0, 0, 0],
                     [0, 0, 0, 0], [f / 2, 0, 0, sp.Rational(1, 2)]])
    yy = s_kron(S_Y, S_Y)
    product = rho * (yy * rho.conjugate() * yy)
    eigenvalues = product.eigenvals()
    expected = {((1 + f) / 2) ** 2, ((1 - f) / 2) ** 2, 0}
    assert {sp.factor(value) for value in eigenvalues} == {sp.factor(value) for value in expected}
    assert eigenvalues[0] == 2
    return f, sp.simplify((1 + f) / 2 - (1 - f) / 2)


def s_directional_cpsi(rho, rho_dot):
    """Exact one-sided d/dt CPsi at t=0 along rho_dot (d = 4).

    Differentiable because no coherence that is zero at t=0 is created at first
    order; the assertion keeps that premise checked rather than assumed."""
    purity = sp.simplify(sp.trace(rho * rho))
    l1 = s_l1(rho)
    purity_dot = sp.simplify(2 * sp.re(sp.trace(rho * rho_dot)))
    l1_dot = 0
    for i in range(rho.rows):
        for j in range(rho.cols):
            if i == j:
                continue
            if rho[i, j] == 0:
                assert sp.simplify(rho_dot[i, j]) == 0, "a new coherence would make l1 non-differentiable"
            else:
                l1_dot += sp.re(sp.conjugate(rho[i, j]) * rho_dot[i, j]) / sp.Abs(rho[i, j])
    return sp.simplify((purity_dot * l1 + purity * sp.simplify(l1_dot)) / (rho.rows - 1))


# ------------------------------------------------------------ numeric helpers
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)


def cpsi(rho):
    d = rho.shape[0]
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    return float(np.real(np.trace(rho @ rho)) * l1 / (d - 1))


def liouvillian(generator, d):
    """Column k is generator(E_k) for the row-major basis matrix E_k."""
    columns = []
    for k in range(d * d):
        basis = np.zeros(d * d, dtype=complex)
        basis[k] = 1.0
        columns.append(generator(basis.reshape(d, d)).reshape(-1))
    return np.column_stack(columns)


# ----------------------------------------------------------------------- checks
def check_1():
    """CPsi = 1/4 is not a separability or quantum/classical boundary."""
    plus = (s_ket([0]) + s_ket([1])) / sp.sqrt(2)
    ppp = s_kron(plus, plus, plus)
    ghz = (s_ket([0, 0, 0]) + s_ket([1, 1, 1])) / sp.sqrt(2)
    return {"cpsi_plus3_separable": str(s_cpsi(s_outer(ppp))),
            "cpsi_ghz3_entangled": str(s_cpsi(s_outer(ghz)))}


def check_2():
    """A local non-Pauli unitary moves CPsi upward."""
    rho00 = s_outer(s_ket([0, 0]))
    hadamard = (S_X + S_Z) / sp.sqrt(2)
    local = s_kron(hadamard, S_I2)
    return {"cpsi_00": str(s_cpsi(rho00)),
            "cpsi_after_hadamard_x_identity": str(s_cpsi(sp.simplify(local * rho00 * local.T)))}


def check_3():
    """Upward crossing of 1/4 under a fixed Markovian semigroup.

    |01>, Heisenberg XX+YY+ZZ (Pauli, J = 1), local Z-dephasing gamma = 0.05 on
    both sites. The trajectory is the matrix exponential of the 16x16 generator
    (no integrator error); the crossing times are brentq roots of that trajectory."""
    hamiltonian = np.kron(X, X) + np.kron(Y, Y) + np.kron(Z, Z)
    gamma = 0.05
    z1, z2 = np.kron(Z, I2), np.kron(I2, Z)

    def generator(rho):
        return (-1j * (hamiltonian @ rho - rho @ hamiltonian)
                + gamma * (z1 @ rho @ z1 - rho) + gamma * (z2 @ rho @ z2 - rho))

    lindblad = liouvillian(generator, 4)
    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[1, 1] = 1.0

    def trajectory(t):
        return cpsi((expm(t * lindblad) @ rho0.reshape(-1)).reshape(4, 4))

    grid = np.arange(0.0, 3.0 + 1e-12, 0.005)
    values = np.array([trajectory(t) for t in grid])
    peak = int(np.argmax(values))
    refined = minimize_scalar(lambda t: -trajectory(t), bounds=(grid[peak - 1], grid[peak + 1]),
                              method="bounded", options={"xatol": 1e-12})
    first_above = int(np.argmax(values > 0.25))
    t_up = brentq(lambda t: trajectory(t) - 0.25, grid[first_above - 1], grid[first_above], xtol=1e-14)
    first_below = peak + int(np.argmax(values[peak:] < 0.25))
    t_down = brentq(lambda t: trajectory(t) - 0.25, grid[first_below - 1], grid[first_below], xtol=1e-14)
    return {"cpsi_at_0": values[0], "cpsi_max": -refined.fun, "t_max": float(refined.x),
            "t_up": t_up, "t_down": t_down}


def check_4b():
    """The proof's exact Z-pulse counterexample (PROOF_MONOTONICITY_CPSI.md Part 5).

    rho0 = ((I + X/2 + Z/2)/2) (x) |0><0|, H = Y (x) I, D(rho) = (Z1 rho Z1 - rho)/4."""
    rho0 = s_kron((S_I2 + S_X / 2 + S_Z / 2) / 2, sp.diag(1, 0))
    hamiltonian = s_kron(S_Y, S_I2)
    z_first = s_kron(S_Z, S_I2)
    jump = z_first / 2
    pulsed = sp.simplify(z_first * rho0 * z_first)
    passive = sp.simplify(z_first * hamiltonian * z_first)
    zero = sp.zeros(4)
    return {"cpsi_rho0": str(s_cpsi(rho0)),
            "cpsi_after_pulse": str(s_cpsi(pulsed)),
            "slope_rho0": str(s_directional_cpsi(rho0, s_rhs(rho0, hamiltonian, (jump,)))),
            "slope_rho0_h_zero": str(s_directional_cpsi(rho0, s_rhs(rho0, zero, (jump,)))),
            "slope_after_pulse_same_h": str(s_directional_cpsi(pulsed, s_rhs(pulsed, hamiltonian, (jump,)))),
            "slope_after_pulse_passive_frame": str(s_directional_cpsi(pulsed, s_rhs(pulsed, passive, (jump,))))}


def check_5():
    """The collective Z1+Z2 dissipator is not the sum of two local Z dissipators.

    The superoperators have integer entries, so the Frobenius norm is exact."""
    z1, z2 = s_kron(S_Z, S_I2), s_kron(S_I2, S_Z)
    columns_local, columns_collective = [], []
    for k in range(16):
        basis = sp.zeros(4)
        basis[k // 4, k % 4] = 1
        columns_local.append(list(s_dissipator(z1, basis) + s_dissipator(z2, basis)))
        columns_collective.append(list(s_dissipator(z1 + z2, basis)))
    difference = sp.Matrix(columns_local).T - sp.Matrix(columns_collective).T
    norm_squared = sum(entry ** 2 for entry in difference)
    return {"norm_squared": str(norm_squared), "norm": str(sp.sqrt(norm_squared))}


def check_6():
    """Bell+ under local Z-dephasing, purity book CPsi = f(1+f^2)/6, f = e^(-4 gamma t).

    The crossing root is the real root of f^3 + f - 3/2 = 0; the concurrence of the
    dephased Bell+ state is f (s_bell_plus_concurrence), so it is f* > 0 there."""
    f, concurrence = s_bell_plus_concurrence()
    root = sp.real_roots(f ** 3 + f - sp.Rational(3, 2))[0]
    return {"f_star": float(sp.N(root, 30)), "k_fold": float(sp.N(-sp.log(root) / 4, 30)),
            "concurrence_formula": str(concurrence),
            "concurrence_at_crossing": float(sp.N(concurrence.subs(f, root), 30))}


def check_7():
    """No chain ratio of the V-Effect census tends to phi.

    F6 gain V(N) = 1 + cos(pi/N); cold single-excitation frequencies
    omega_m = 4J(1 - cos(pi m/N)) (simulations/n5_optimal_cavity_size.py)."""
    n = sp.symbols("N", positive=True)
    phi = (1 + sp.sqrt(5)) / 2
    gain = 1 + sp.cos(sp.pi / n)
    ratio_21 = (1 - sp.cos(2 * sp.pi / n)) / (1 - sp.cos(sp.pi / n))
    return {"V_limit": str(sp.limit(gain, n, sp.oo)),
            "V5_minus_1_plus_phi_half": str(sp.simplify(gain.subs(n, 5) - (1 + phi / 2))),
            "omega21_N5_minus_2_plus_phi": str(sp.simplify(sp.nsimplify(ratio_21.subs(n, 5)) - (2 + phi))),
            "omega21_limit": str(sp.limit(ratio_21, n, sp.oo)),
            "phi": float(sp.N(phi, 30)),
            "V5": float(sp.N(gain.subs(n, 5), 30))}


def check_8():
    """The free-|+> crossing and the Bell+ fold are two different numbers.

    Free |+> in the pure-dephasing limit: x^3 + x = 1/2 with x = e^(-t*/T2)."""
    x = sp.symbols("x", positive=True)
    root = sp.real_roots(x ** 3 + x - sp.Rational(1, 2))[0]
    return {"x_root": float(sp.N(root, 30)), "t_star_over_T2": float(sp.N(-sp.log(root), 30)),
            "bell_k_fold": check_6()["k_fold"]}


def check_9():
    """The historical 120-pair catalogue is C(16, 2), not C(6, 2)."""
    return {"C_6_2": comb(6, 2), "C_16_2": comb(16, 2)}


def check_10():
    """The April-26 cusp-precision payloads hold no CPsi row exactly at 1/4.

    Reads data/ibm_cusp_precision_april2026/*.json and counts the stored cpsi
    rows, not the delay factors (one delay factor is 0.25)."""
    rows = {}
    for path in sorted((ROOT / "data/ibm_cusp_precision_april2026").glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        values = [row["cpsi"] for row in payload["cpsi_data"]]
        rows[path.name] = {"rows": len(values),
                           "exactly_quarter": sum(value == 0.25 for value in values),
                           "within_0_01_of_quarter": sum(abs(value - 0.25) < 0.01 for value in values),
                           "phase_keys": path.read_text(encoding="utf-8").count('"phase')}
    return {"files": len(rows), "per_file": rows,
            "exactly_quarter_total": sum(entry["exactly_quarter"] for entry in rows.values())}


def check_12():
    """The Q52 '115.0 vs 114.7' pair is one hardware record read twice.

    Reads data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json.
    crossing_us is the run's stored crossing of the MEASURED CPsi; the second
    reading is the linear interpolation of the stored density matrices that
    simulations/cockpit_ibm_hardware.py:309-317 performs and prints at :319-324
    as 'Measured' against the stored value printed as 'Predicted'."""
    path = ROOT / "data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    times, values = [], []
    for point in record["raw_tomography"]:
        rho = np.array(point["density_matrix_real"]) + 1j * np.array(point["density_matrix_imag"])
        times.append(point["delay_us"])
        values.append(cpsi(rho))
    interpolated = None
    for i in range(1, len(values)):
        if values[i - 1] >= 0.25 > values[i]:
            fraction = (0.25 - values[i]) / (values[i - 1] - values[i])
            interpolated = times[i] * (1 - fraction) + times[i - 1] * fraction
            break
    stored = record["crossing_us"]
    return {"stored_crossing_us": stored, "interpolated_crossing_us": interpolated,
            "relative_gap": abs(interpolated - stored) / stored,
            "stored_pure_dephasing_prediction_t_over_T2": record["analytical_prediction_pure_dephasing"]}


def check_13():
    """One Bell+ trajectory, three CPsi books, three crossing doses.

    f = e^(-4 gamma t); concurrence book f^2/3, purity book f(1+f^2)/6, constant
    book f/3 (simulations/subsystem_crossing_pairs.py:18-23). The dephased Bell+
    concurrence is f, positive at every finite t, so the concurrence-book crossing
    at f = sqrt(3)/2 is not where the pair stops being entangled."""
    f, concurrence = s_bell_plus_concurrence()
    f_concurrence = sp.sqrt(3) / 2
    assert sp.simplify(f_concurrence ** 2 / 3 - sp.Rational(1, 4)) == 0
    f_purity = sp.real_roots(f ** 3 + f - sp.Rational(3, 2))[0]
    f_constant = sp.Rational(3, 4)
    doses = {name: sp.simplify(-sp.log(value) / 4)
             for name, value in (("concurrence", f_concurrence), ("purity", f_purity),
                                 ("constant", f_constant))}
    return {"f_concurrence_book": str(f_concurrence),
            "k_concurrence_book": str(doses["concurrence"]),
            "k_concurrence_book_value": float(sp.N(doses["concurrence"], 30)),
            "k_purity_book_value": float(sp.N(doses["purity"], 30)),
            "k_constant_book_value": float(sp.N(doses["constant"], 30)),
            "t_at_gamma_0_05": [round(float(sp.N(doses[name], 30)) / 0.05, 4)
                                for name in ("concurrence", "purity", "constant")],
            "concurrence_formula": str(concurrence),
            "concurrence_at_concurrence_crossing": str(sp.simplify(concurrence.subs(f, f_concurrence)))}


def check_14():
    """Whether a trajectory crosses 1/4 depends on the book, not on the physics alone.

    |++> under Heisenberg XX+YY+ZZ with local Z-dephasing gamma on both sites. The
    product family rho(a) = rho1(a) (x) rho1(a), rho1 = (I + aX)/2, a = e^(-2 gamma t),
    solves the Lindblad equation exactly ([H, rho(a)] = 0 and D(rho(a)) =
    -2 gamma a d rho/da), so the state stays a product state. Its spin-flipped
    product rho * (YY rho* YY) is ((1-a^2)/4)^2 times the identity, so the four
    Wootters roots are equal and the concurrence max(0, s - 3s) is 0 for 0 < a <= 1:
    the concurrence book stays at 0, while the purity book starts at 1 and crosses 1/4."""
    a, gamma = sp.symbols("a gamma", positive=True)
    rho1 = (S_I2 + a * S_X) / 2
    rho = s_kron(rho1, rho1)
    hamiltonian = s_kron(S_X, S_X) + s_kron(S_Y, S_Y) + s_kron(S_Z, S_Z)
    z1, z2 = s_kron(S_Z, S_I2), s_kron(S_I2, S_Z)
    commutator = sp.simplify(hamiltonian * rho - rho * hamiltonian)
    dephasing = sp.simplify(gamma * (z1 * rho * z1 - rho) + gamma * (z2 * rho * z2 - rho))
    flow = sp.simplify(dephasing + 2 * gamma * a * rho.diff(a))
    yy = s_kron(S_Y, S_Y)
    root = (1 - a ** 2) / 4
    spin_flip_is_scalar = sp.simplify(rho * (yy * rho.conjugate() * yy) - root ** 2 * sp.eye(4)) == sp.zeros(4)
    concurrence_samples = [str(sp.Max(0, root - 3 * root).subs(a, sp.Rational(k, 10))) for k in range(1, 11)]
    purity_book = sp.factor(sp.simplify(sp.trace(rho * rho) * s_l1(rho) / 3))
    crossing = [value for value in sp.real_roots(sp.numer(sp.together(purity_book - sp.Rational(1, 4))))
                if 0 < value < 1]
    return {"commutator_is_zero": commutator == sp.zeros(4),
            "product_family_solves_lindblad": flow == sp.zeros(4),
            "spin_flip_product_is_scalar": spin_flip_is_scalar,
            "concurrence_at_a_0_1_to_1": concurrence_samples,
            "purity_book": str(purity_book),
            "purity_book_at_a_1": str(purity_book.subs(a, 1)),
            "purity_book_crossing_a": float(sp.N(crossing[0], 30)),
            "purity_book_crossing_count_in_unit_interval": len(crossing)}


CHECKS = {
    "CHECK-1": check_1, "CHECK-2": check_2, "CHECK-3": check_3, "CHECK-4b": check_4b,
    "CHECK-5": check_5, "CHECK-6": check_6, "CHECK-7": check_7, "CHECK-8": check_8,
    "CHECK-9": check_9, "CHECK-10": check_10, "CHECK-12": check_12, "CHECK-13": check_13,
    "CHECK-14": check_14,
}


def report(values):
    lines = {
        "CHECK-1": "CPsi(|+++>, separable) = {cpsi_plus3_separable} ; CPsi(GHZ3, entangled) = {cpsi_ghz3_entangled}",
        "CHECK-2": "CPsi(|00>) = {cpsi_00} -> after H x I: {cpsi_after_hadamard_x_identity}",
        "CHECK-3": "|01>: CPsi(0) = {cpsi_at_0:.3f}, max CPsi = {cpsi_max:.3f} at t = {t_max:.3f}; "
                   "crosses 1/4 upward at t = {t_up:.3f}, back down at t = {t_down:.3f}",
        "CHECK-4b": "CPsi(rho0) = {cpsi_rho0}; CPsi'(0) = {slope_rho0}; H->0: {slope_rho0_h_zero}; "
                    "after Z1 pulse CPsi = {cpsi_after_pulse}, same H: {slope_after_pulse_same_h}; "
                    "passive frame Z1HZ1: {slope_after_pulse_passive_frame}",
        "CHECK-5": "||D[Z1]+D[Z2] - D[Z1+Z2]||^2 = {norm_squared} -> norm {norm} (0 would mean equivalent)",
        "CHECK-6": "Bell+ purity book f* = {f_star:.6f} at CPsi = 1/4; concurrence there = {concurrence_at_crossing:.6f} > 0; "
                   "K_fold = -ln(f*)/4 = {k_fold:.5f}",
        "CHECK-7": "V(N) -> {V_limit}; V(5) - (1 + phi/2) = {V5_minus_1_plus_phi_half}; omega2/omega1 at N=5 minus (2 + phi) = "
                   "{omega21_N5_minus_2_plus_phi}; omega2/omega1 -> {omega21_limit}; phi = {phi:.4f}",
        "CHECK-8": "free-|+> root x = {x_root:.6f}, t*/T2 = {t_star_over_T2:.6f}; Bell+ K_fold = {bell_k_fold:.5f} -> distinct crossings",
        "CHECK-9": "C(6,2) = {C_6_2} ; C(16,2) = {C_16_2} = pairs with repetition from 15 words",
        "CHECK-10": "cusp payload files = {files} ; cpsi rows exactly 1/4 = {exactly_quarter_total} ; per file {per_file}",
        "CHECK-12": "Q52 stored crossing_us = {stored_crossing_us:.3f} us ; interpolated from the same stored density "
                    "matrices = {interpolated_crossing_us:.3f} us ; relative gap {relative_gap:.4f}",
        "CHECK-13": "Bell+ books: K = {k_concurrence_book_value:.6f} (concurrence), {k_purity_book_value:.6f} (purity), "
                    "{k_constant_book_value:.6f} (constant); t at gamma=0.05 = {t_at_gamma_0_05}; concurrence at the "
                    "concurrence-book crossing = {concurrence_at_concurrence_crossing} > 0",
        "CHECK-14": "|++> Heisenberg + local Z: [H, rho] = 0: {commutator_is_zero}; product family solves Lindblad: "
                    "{product_family_solves_lindblad}; spin-flip product scalar: {spin_flip_product_is_scalar}; "
                    "purity book starts at {purity_book_at_a_1} and crosses at a = {purity_book_crossing_a:.6f}; "
                    "concurrence at a = 0.1..1: {concurrence_at_a_0_1_to_1}",
    }
    return [f"{name} " + lines[name].format(**values[name]) for name in CHECKS]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    values = {name: check() for name, check in CHECKS.items()}
    if "--json" in argv:
        print(json.dumps(values, sort_keys=True))
    else:
        for line in report(values):
            print(line)


if __name__ == "__main__":
    main()
