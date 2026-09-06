#!/usr/bin/env python3
"""Basis-independent N=3 observable time traces under Z dephasing.

This replaces eigenvector-coordinate 'state weights', which are not invariant
for a non-normal Liouvillian and become basis-dependent in degenerate spaces.
The output reports only direct expectation-value ranges.
"""

from itertools import product
from pathlib import Path
import numpy as np
from scipy.linalg import expm

RESULTS_DIR = Path(__file__).parent / "results"
N, J, GAMMA = 3, 1.0, 0.05
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


def tensor(ops):
    result = ops[0]
    for operator in ops[1:]:
        result = np.kron(result, operator)
    return result


def operator(label):
    return tensor([PAULI[letter] for letter in label])


def liouvillian():
    dimension = 2**N
    identity = np.eye(dimension, dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for site in range(N - 1):
        for pauli in (X, Y, Z):
            ops = [I2] * N
            ops[site], ops[site + 1] = pauli, pauli
            hamiltonian += J * tensor(ops)
    result = -1j * (np.kron(hamiltonian, identity) - np.kron(identity, hamiltonian.T))
    for site in range(N):
        ops = [I2] * N
        ops[site] = Z
        jump = np.sqrt(GAMMA) * tensor(ops)
        norm = jump.conj().T @ jump
        result += np.kron(jump, jump.conj()) - 0.5 * (
            np.kron(norm, identity) + np.kron(identity, norm.T)
        )
    return result


def ket(bits):
    vector = np.zeros(2**N, dtype=complex)
    vector[int(bits, 2)] = 1
    return vector


def states():
    zero, one = np.array([1, 0]), np.array([0, 1])
    plus = (zero + one) / np.sqrt(2)
    return {
        "GHZ": (ket("000") + ket("111")) / np.sqrt(2),
        "W": (ket("100") + ket("010") + ket("001")) / np.sqrt(3),
        "Bell01": (ket("000") + ket("110")) / np.sqrt(2),
        "+++": tensor([plus, plus, plus]),
    }


def direct_trace_gate(generator=None):
    """Return fixed N=3 anchors, raising if the direct dynamics are corrupted.

    The signed W/IYY derivative makes this gate sensitive to a zero generator
    and to a reversed dissipative direction.  The trace derivative separately
    pins trace preservation on the same moving state and generator path.
    """
    if generator is None:
        generator = liouvillian()
    psi = states()["W"]
    rho0_matrix = np.outer(psi, psi.conj())
    rho0 = rho0_matrix.reshape(-1)
    derivative = (generator @ rho0).reshape(2**N, 2**N)
    observable = operator("IYY")
    anchors = {
        "w_iyy_t0": float(np.trace(observable @ rho0_matrix).real),
        "w_iyy_dt0": float(np.trace(observable @ derivative).real),
        "w_state_dt0_norm": float(np.linalg.norm(derivative)),
        "trace_dt0": float(abs(np.trace(derivative))),
    }
    scale = max(1.0, float(np.linalg.norm(generator)))
    tolerance = 256 * np.finfo(float).eps * scale
    failures = []
    if abs(anchors["w_iyy_t0"] - 2.0 / 3.0) > tolerance:
        failures.append("W/IYY t=0 anchor")
    if abs(anchors["w_iyy_dt0"] - (-2.0 / 15.0)) > tolerance:
        failures.append("W/IYY signed derivative")
    if anchors["w_state_dt0_norm"] <= 0.1:
        failures.append("moving-state derivative norm")
    if anchors["trace_dt0"] > tolerance:
        failures.append("trace derivative")
    if failures:
        raise RuntimeError("direct Pauli trace gate failed: " + ", ".join(failures))
    return anchors


def trace_rows(times=np.linspace(0, 10, 101), generator=None):
    if generator is None:
        generator = liouvillian()
    observables = ("ZZZ", "IYY", "XXZ", "ZXX", "YYI", "XZX", "YIY")
    propagators = [expm(generator * time) for time in times]
    rows = []
    for state_name, psi in states().items():
        rho0 = np.outer(psi, psi.conj()).reshape(-1)
        for label in observables:
            observable = operator(label)
            values = np.array([
                np.trace(observable @ (propagator @ rho0).reshape(2**N, 2**N)).real
                for propagator in propagators
            ])
            rows.append((state_name, label, values[0], values.min(), values.max(), (values.max()-values.min())/2))
    return rows


def main():
    generator = liouvillian()
    anchors = direct_trace_gate(generator)
    lines = [
        "N=3 DIRECT PAULI-OBSERVABLE TIME TRACES",
        f"Heisenberg chain, J={J}, gamma={GAMMA}, t=0..10 (101 samples)",
        "Direct expectations are basis-independent; no eigenmode state-weight percentage is reported.",
        "Amplitude is half the sampled max-min range, not a standing-wave certificate.",
        (
            "Dynamic gate: W/IYY t0="
            f"{anchors['w_iyy_t0']:.6f}, dt0={anchors['w_iyy_dt0']:.6f}; "
            f"||drho/dt||={anchors['w_state_dt0_norm']:.6f}, "
            f"|Tr drho/dt|={anchors['trace_dt0']:.2e}."
        ),
        "",
        f"{'state':<9} {'Pauli':<5} {'t0':>11} {'min':>11} {'max':>11} {'half-range':>11}",
    ]
    for row in trace_rows(generator=generator):
        lines.append(f"{row[0]:<9} {row[1]:<5} {row[2]:11.6f} {row[3]:11.6f} {row[4]:11.6f} {row[5]:11.6f}")
    lines += ["", "Scope: observable traces only; spatial counter-propagation and mode-pair phase relations were not tested."]
    output = "\n".join(lines) + "\n"
    print(output, end="")
    path = RESULTS_DIR / "standing_wave_analysis.txt"
    path.write_text(output, encoding="utf-8")
    print(f"\nResults saved to: {path}")


if __name__ == "__main__":
    main()
