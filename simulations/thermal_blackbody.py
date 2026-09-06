#!/usr/bin/env python3
"""Finite thermal-channel census; no blackbody or wave mechanism is inferred."""

from pathlib import Path
import os
import sys
import numpy as np

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

RESULTS_DIR = Path(__file__).parent / "results"
J, GAMMA, TOL_FREQ = 1.0, 0.05, 1e-6
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Sm = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma- = |0><1|: emission
Sp = np.array([[0, 0], [1, 0]], dtype=complex)  # sigma+ = |1><0|: absorption


def kron_chain(ops):
    result = ops[0]
    for operator in ops[1:]:
        result = np.kron(result, operator)
    return result


def dissipator_term(jump, dimension):
    identity = np.eye(dimension, dtype=complex)
    norm = jump.conj().T @ jump
    return np.kron(jump, jump.conj()) - 0.5 * (
        np.kron(norm, identity) + np.kron(identity, norm.T)
    )


def build_thermal_liouvillian(n, gammas, n_bar, gamma_thermal=None):
    """Build Z dephasing plus emission/absorption for n_bar >= 0."""
    if n_bar < 0:
        raise ValueError("n_bar must be non-negative")
    if gamma_thermal is None:
        gamma_thermal = gammas
    if len(gammas) != n or len(gamma_thermal) != n:
        raise ValueError("one dephasing and thermal rate is required per site")
    dimension = 2**n
    identity = np.eye(dimension, dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for site in range(n - 1):
        for pauli in (Xm, Ym, Zm):
            ops = [I2] * n
            ops[site], ops[site + 1] = pauli, pauli
            hamiltonian += J * kron_chain(ops)
    generator = -1j * (np.kron(hamiltonian, identity) - np.kron(identity, hamiltonian.T))
    for site in range(n):
        ops = [I2] * n
        ops[site] = Zm
        generator += dissipator_term(np.sqrt(gammas[site]) * kron_chain(ops), dimension)
        ops = [I2] * n
        ops[site] = Sm
        generator += dissipator_term(
            np.sqrt(gamma_thermal[site] * (n_bar + 1.0)) * kron_chain(ops), dimension
        )
        if n_bar > 0:
            ops = [I2] * n
            ops[site] = Sp
            generator += dissipator_term(
                np.sqrt(gamma_thermal[site] * n_bar) * kron_chain(ops), dimension
            )
    return generator


def spectral_row(values):
    oscillating = np.abs(values.imag) > TOL_FREQ
    rates = -values.real
    q = np.abs(values.imag[oscillating]) / np.maximum(rates[oscillating], 1e-15)
    return int(np.count_nonzero(oscillating)), float(np.mean(rates)), float(np.max(q))


def main():
    lines = []
    def log(message=""):
        print(message)
        lines.append(message)

    n = 4
    gammas = [GAMMA] * n
    log("THERMAL-CHANNEL NUMERICAL CENSUS")
    log("sigma- = |0><1| emission; sigma+ = |1><0| absorption")
    log("L(nbar)=L_H+L_Z+gamma_T[(nbar+1)D[sigma-]+nbar D[sigma+]]")
    log(f"N={n}, J={J}, gamma_Z=gamma_T={GAMMA}, imag tolerance={TOL_FREQ:.0e}")
    log("Counts are numerical at the stated tolerance; they are not EP certificates.")
    log()
    log(f"{'nbar':>10} {'osc':>6} {'osc%':>8} {'mean decay':>12} {'Qmax':>10}")
    for n_bar in (0.0, 1e-9, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0):
        values = np.linalg.eigvals(build_thermal_liouvillian(n, gammas, n_bar, gammas))
        osc, mean_decay, qmax = spectral_row(values)
        log(f"{n_bar:10.3g} {osc:6d} {100*osc/len(values):7.3f}% {mean_decay:12.6f} {qmax:10.3f}")
    cold = build_thermal_liouvillian(1, [0.0], 0.0, [GAMMA])
    near = build_thermal_liouvillian(1, [0.0], 1e-9, [GAMMA])
    log()
    log(f"negative control: ||L(0)||_F={np.linalg.norm(cold):.6e} > 0 (spontaneous emission present)")
    log(f"continuity: ||L(1e-9)-L(0)||_F={np.linalg.norm(near-cold):.6e}")
    log("Scope: finite spectral counts only; no thermal-photon, standing-wave, or phase-transition mechanism is inferred.")
    path = RESULTS_DIR / "thermal_blackbody.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nResults saved to: {path}")


if __name__ == "__main__":
    main()
