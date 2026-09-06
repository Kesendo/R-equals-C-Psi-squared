"""Audit what the local finite-occupation channel actually fixes.

The earlier fixed-point construction compared the steady state of local
sigma-/sigma+ jumps with a Gibbs state of the interacting Heisenberg
Hamiltonian. Those are different bath models, so their energy difference
cannot be interpreted as heat produced by the simulated system. This
producer keeps only the channel statement that can be checked from below.
"""

from pathlib import Path

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1|
SP = np.array([[0, 0], [1, 0]], dtype=complex)  # |1><0|


def kron_at(op, target, n):
    result = np.eye(1, dtype=complex)
    for site in range(n):
        result = np.kron(result, op if site == target else I2)
    return result


def build_h(n, coupling=1.0):
    dim = 2**n
    hamiltonian = np.zeros((dim, dim), dtype=complex)
    for site in range(n - 1):
        for pauli in (X, Y, Z):
            hamiltonian += (coupling * kron_at(pauli, site, n)
                            @ kron_at(pauli, site + 1, n))
    return hamiltonian


def add_dissipator(generator, jump):
    dim = jump.shape[0]
    identity = np.eye(dim, dtype=complex)
    jdj = jump.conj().T @ jump
    generator += np.kron(jump.conj(), jump)
    generator -= 0.5 * np.kron(identity, jdj)
    generator -= 0.5 * np.kron(jdj.T, identity)


def build_generator(hamiltonian, gamma_z, gamma_amp, n_bar):
    dim = hamiltonian.shape[0]
    n = int(np.log2(dim))
    identity = np.eye(dim, dtype=complex)
    generator = -1j * (np.kron(identity, hamiltonian)
                       - np.kron(hamiltonian.T, identity))
    for site in range(n):
        gz = gamma_z[site] if hasattr(gamma_z, "__len__") else gamma_z
        if gz > 0:
            add_dissipator(generator, np.sqrt(gz) * kron_at(Z, site, n))
        add_dissipator(
            generator,
            np.sqrt(gamma_amp * (1 + n_bar)) * kron_at(SM, site, n),
        )
        if n_bar > 0:
            add_dissipator(
                generator,
                np.sqrt(gamma_amp * n_bar) * kron_at(SP, site, n),
            )
    return generator


def local_target(n_bar):
    """Stationary one-site state fixed by the chosen jump-rate ratio."""
    p_excited = n_bar / (2 * n_bar + 1)
    return np.diag([1 - p_excited, p_excited]).astype(complex)


def tensor_power(matrix, n):
    result = np.eye(1, dtype=complex)
    for _ in range(n):
        result = np.kron(result, matrix)
    return result


def stationarity_residual(n, gamma_z, gamma_amp, n_bar, target=None):
    hamiltonian = build_h(n)
    generator = build_generator(hamiltonian, gamma_z, gamma_amp, n_bar)
    rho = tensor_power(local_target(n_bar) if target is None else target, n)
    return np.linalg.norm(generator @ rho.reshape(-1, order="F"))


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "self_heating_fixpoint.txt"
    lines = []

    def out(message=""):
        print(message)
        lines.append(message)

    out("LOCAL FINITE-OCCUPATION CHANNEL AUDIT")
    out("sigma- = |0><1|, sigma+ = |1><0|")
    out("n_bar is supplied externally; this generator contains no feedback law.")
    out()
    out("One-qubit directional endpoint gate")
    out("n_bar=0: |0><0| fixed; |1><1| flows toward |0><0|.")

    gamma_amp = 0.2
    cold = build_generator(np.zeros((2, 2), dtype=complex), 0.0,
                           gamma_amp, 0.0)
    rho0 = np.diag([1.0, 0.0]).astype(complex)
    rho1 = np.diag([0.0, 1.0]).astype(complex)
    d0 = (cold @ rho0.reshape(-1, order="F")).reshape(2, 2, order="F")
    d1 = (cold @ rho1.reshape(-1, order="F")).reshape(2, 2, order="F")
    endpoint_residual = max(abs(d0[1, 1]), abs(d1[1, 1] + gamma_amp))
    out(f"endpoint residual = {endpoint_residual:.2e}")

    wrong_target = np.diag([0.25, 0.75]).astype(complex)
    correct_residual = stationarity_residual(1, 0.0, 0.1, 0.5)
    wrong_residual = stationarity_residual(
        1, 0.0, 0.1, 0.5, target=wrong_target
    )
    out(f"one-site target residual at n_bar=0.5 = {correct_residual:.2e}")
    out(f"wrong-population control residual = {wrong_residual:.2e}")
    if endpoint_residual > 1e-14 or correct_residual > 1e-14 \
            or wrong_residual < 1e-2:
        raise AssertionError("direction/stationarity gate failed")

    out()
    out("Interacting-chain stationarity of the local-channel target")
    out(f"{'N':>2} {'n_bar':>7} {'gamma_z':>9} {'||L rho_target||':>18}")
    for n, n_bar, gamma_z in (
        (3, 0.0, 0.0),
        (3, 0.5, 0.1),
        (5, 0.5, 0.1),
        (5, 2.0, [0.5, 0.01, 0.01, 0.01, 0.01]),
    ):
        residual = stationarity_residual(n, gamma_z, 0.05, n_bar)
        out(f"{n:2d} {n_bar:7.2f} {str(gamma_z):>9} {residual:18.2e}")
        if residual > 1e-12:
            raise AssertionError("local-channel product target is not stationary")

    out()
    out("This checks a local bath target, not a Gibbs state of the interacting H.")
    out("No self-heating, heat-production, cooling-rate, or biological claim follows.")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
