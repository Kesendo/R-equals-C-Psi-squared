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

    The signed W/IYY derivative pins the dissipative direction.  The signed
    Bell01/IXY derivative is zero when the Hamiltonian term is deleted and
    therefore pins the coherent path independently.  The trace derivative
    separately pins trace preservation on the same generator path.
    """
    if generator is None:
        generator = liouvillian()
    psi = states()["W"]
    rho0_matrix = np.outer(psi, psi.conj())
    rho0 = rho0_matrix.reshape(-1)
    derivative = (generator @ rho0).reshape(2**N, 2**N)
    observable = operator("IYY")
    bell = states()["Bell01"]
    bell_rho0 = np.outer(bell, bell.conj()).reshape(-1)
    bell_derivative = (generator @ bell_rho0).reshape(2**N, 2**N)
    anchors = {
        "w_iyy_t0": float(np.trace(observable @ rho0_matrix).real),
        "w_iyy_dt0": float(np.trace(observable @ derivative).real),
        "w_state_dt0_norm": float(np.linalg.norm(derivative)),
        "trace_dt0": float(abs(np.trace(derivative))),
        "bell01_ixy_dt0": float(np.trace(operator("IXY") @ bell_derivative).real),
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
    if abs(anchors["bell01_ixy_dt0"] - (-2.0)) > tolerance:
        failures.append("Hamiltonian-sensitive Bell01/IXY derivative")
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


def spectral_census(generator):
    """The eigenvalue half of the reading, which the state-weight defect never touched.

    Eigenvalues are invariant under any change of eigenvector basis and under
    degeneracy, so the pairing, the centred classification and the count of
    standing-wave candidates survive the withdrawal of eigenvector-coordinate
    weights. What made those weights unusable was the choice of coordinates;
    none of that reaches the spectrum.

    The classification never rests on a bare threshold. It reports the gap that
    separates the two sides, so a reader can see whether the verdict is a
    measurement or an artefact of where a line was drawn.
    """
    sigma = N * GAMMA
    values = np.linalg.eigvals(generator)
    mu = values + sigma                       # centred: F1 sends mu -> -mu
    scale = float(np.linalg.norm(generator))
    solver_noise = 256 * np.finfo(float).eps * scale

    # F1 pairing on the centred spectrum, matched greedily with multiplicity.
    # The residual is only meaningful alongside the count: an unmatched mode
    # contributes no distance, so a run that pairs nothing would otherwise
    # report a perfect residual. Both are gated together below.
    remaining = list(range(len(mu)))
    pairs, worst_residual = 0, 0.0
    unmatched = 0
    while remaining:
        i = remaining.pop(0)
        best, best_distance = None, float("inf")
        for j in remaining:
            distance = abs(mu[i] + mu[j])
            if distance < best_distance:
                best, best_distance = j, distance
        # The match window is the noise floor with room to spare, not a free
        # number: F1 is exact here, so a genuine partner sits at the floor and
        # anything a thousand floors away is not one.
        if best is not None and best_distance <= 1000 * solver_noise:
            remaining.remove(best)
            pairs += 1
            worst_residual = max(worst_residual, best_distance)
        else:
            unmatched += 1

    real_parts = np.abs(mu.real)
    imag_parts = np.abs(mu.imag)
    oscillatory = int(np.count_nonzero((real_parts <= solver_noise) & (imag_parts > solver_noise)))
    purely_real = int(np.count_nonzero(imag_parts <= solver_noise))
    mixed = len(mu) - oscillatory - purely_real

    off_axis = real_parts[real_parts > solver_noise]
    closest = float(np.min(off_axis)) if off_axis.size else float("nan")

    lines = [
        "CENTRED SPECTRAL CENSUS",
        f"mu = lambda + Sigma_gamma, Sigma_gamma = {sigma:.6f}; F1 acts as mu -> -mu.",
        f"{len(mu)} eigenvalues; eigensolver noise floor 256*eps*||L|| = {solver_noise:.2e}.",
        "",
        f"  F1 pairing: {2 * pairs}/{len(mu)} matched = {pairs} pairs, "
        f"{unmatched} unmatched; worst |mu_k + mu_partner| = {worst_residual:.2e}",
        f"  purely imaginary mu (standing-wave candidates): {oscillatory}",
        f"  purely real mu (no oscillation):                {purely_real}",
        f"  mixed decay and oscillation:                    {mixed}",
        "",
    ]
    # Gates. Each states something the run could contradict.
    # (1) F1 is exact here, so every mode must find a partner. A count that
    #     drops is the finding; the residual alone cannot report it, because an
    #     unmatched mode contributes no distance.
    if unmatched != 0 or 2 * pairs != len(mu):
        raise RuntimeError(
            f"F1 pairing incomplete: {2 * pairs}/{len(mu)} matched, {unmatched} unmatched")
    if worst_residual > 1000 * solver_noise:
        raise RuntimeError(
            f"F1 pairing residual {worst_residual:.2e} exceeds the match window "
            f"{1000 * solver_noise:.2e}")
    # (2) The three classes must exhaust the spectrum.
    if oscillatory + purely_real + mixed != len(mu):
        raise RuntimeError("the centred classification does not partition the spectrum")
    # (3) The standing-wave count is only readable if nothing sits near the axis.
    #     A mode inside a decade of the floor would make the count a threshold effect.
    if oscillatory == 0 and closest <= 10 * solver_noise:
        raise RuntimeError(
            f"the nearest off-axis mode is {closest:.3e}, within a decade of the "
            f"{solver_noise:.3e} floor; the zero count cannot be read")

    if oscillatory == 0:
        lines += [
            f"  No mode sits on the imaginary axis. The nearest one is {closest:.6f} away,",
            f"  which is {closest / solver_noise:.1e} times the noise floor, so the count is a",
            "  measurement and not a rounding: at N=3 no centred pair has a FLAT envelope.",
            "",
            "  The 20 mixed pairs do share an envelope, exactly: each is a conjugate pair",
            "  with a common real part, and the distinct values are the whole story",
            f"  ({', '.join(f'{v:+.7f}' for v in sorted(set(np.round(mu[(np.abs(mu.real) > solver_noise) & (np.abs(mu.imag) > solver_noise)].real, 7))))}).",
            "  What condition 1 asks for is not a shared envelope but a flat one, Re mu = 0,",
            "  so the two members neither grow nor decay relative to each other. That is what",
            "  is absent, and conditions 2 and 3 are never reached.",
        ]
    else:
        lines.append(f"  {oscillatory} modes sit within the noise floor of the imaginary axis.")
    lines.append("")
    return "\n".join(lines)


def render_report():
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
            f"|Tr drho/dt|={anchors['trace_dt0']:.2e}; "
            f"Bell01/IXY dt0={anchors['bell01_ixy_dt0']:.6f}."
        ),
        "",
        spectral_census(generator),
        f"{'state':<9} {'Pauli':<5} {'t0':>11} {'min':>11} {'max':>11} {'half-range':>11}",
    ]
    for row in trace_rows(generator=generator):
        lines.append(f"{row[0]:<9} {row[1]:<5} {row[2]:11.6f} {row[3]:11.6f} {row[4]:11.6f} {row[5]:11.6f}")
    lines += ["", "Scope: observable traces only; spatial counter-propagation and mode-pair phase relations were not tested."]
    return "\n".join(lines) + "\n"


def main():
    output = render_report()
    print(output, end="")
    path = RESULTS_DIR / "standing_wave_analysis.txt"
    path.write_text(output, encoding="utf-8")
    print(f"\nResults saved to: {path}")


if __name__ == "__main__":
    main()
