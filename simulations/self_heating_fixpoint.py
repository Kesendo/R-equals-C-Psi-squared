"""
Self-Heating Fixed Point: Does the resonator find its own temperature?

The feedback loop:
  decay -> heat -> new modes -> more decay -> more heat

This script finds the self-consistent n_bar where the steady-state
energy of the Liouvillian matches the thermal energy at that n_bar.

Method: fixed-point iteration on n_bar.
  1. Start with n_bar = 0
  2. Build Liouvillian L(gamma_z, gamma_amp, n_bar)
  3. Find steady state rho_ss (null eigenvector of L)
  4. Compute E_ss = Tr(H * rho_ss)
  5. Compute E_thermal(n_bar) = Tr(H * rho_thermal(n_bar))
  6. Adjust n_bar until E_ss = E_thermal
"""

import numpy as np
from scipy.linalg import null_space
from pathlib import Path


# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)
SM = np.array([[0, 0], [1, 0]], dtype=complex)


def kron_at(op, t, n):
    r = np.eye(1, dtype=complex)
    for k in range(n):
        r = np.kron(r, op if k == t else I2)
    return r


def build_H(n, J=1.0):
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i in range(n - 1):
        for p in [X, Y, Z]:
            H += J * kron_at(p, i, n) @ kron_at(p, i + 1, n)
    return H


def build_L_thermal(H, gamma_z, gamma_amp, n_bar):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n_qubits = int(np.log2(d))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(n_qubits):
        gz = gamma_z[k] if hasattr(gamma_z, '__len__') else gamma_z
        ga = gamma_amp[k] if hasattr(gamma_amp, '__len__') else gamma_amp
        for Lop in ([np.sqrt(gz) * kron_at(Z, k, n_qubits)] if gz > 0 else []) + \
                   ([np.sqrt(ga * (1 + n_bar)) * kron_at(SM, k, n_qubits)] if ga > 0 else []) + \
                   ([np.sqrt(ga * n_bar) * kron_at(SP, k, n_qubits)] if ga > 0 and n_bar > 0 else []):
            LdL = Lop.conj().T @ Lop
            L += np.kron(Lop.conj(), Lop)
            L -= 0.5 * np.kron(Id, LdL)
            L -= 0.5 * np.kron(LdL.T, Id)
    return L


def find_steady_state(L, d):
    """Find rho_ss such that L @ rho_ss = 0."""
    # The steady state is the null eigenvector of L
    evals, evecs = np.linalg.eig(L)
    # Find eigenvalue closest to 0
    idx = np.argmin(np.abs(evals))
    rho_vec = evecs[:, idx]
    rho = rho_vec.reshape(d, d)
    rho = (rho + rho.conj().T) / 2  # enforce hermiticity
    rho /= np.trace(rho)  # normalize
    return rho


def thermal_state(H, n_bar):
    """Compute thermal state rho = exp(-beta*H) / Z."""
    d = H.shape[0]
    n_qubits = int(np.log2(d))
    if n_bar < 1e-12:
        # Zero temperature: ground state
        evals, evecs = np.linalg.eigh(H)
        rho = np.outer(evecs[:, 0], evecs[:, 0].conj())
        return rho
    # beta from n_bar: for a qubit with splitting Delta,
    # n_bar = 1/(exp(beta*Delta) - 1), so beta = ln(1+1/n_bar) / Delta
    # Use the mean level spacing as Delta
    evals = np.linalg.eigvalsh(H)
    Delta = (max(evals) - min(evals)) / (d - 1)
    if Delta < 1e-15:
        return np.eye(d, dtype=complex) / d
    beta = np.log(1 + 1 / n_bar) / Delta
    rho = np.zeros((d, d), dtype=complex)
    evals_H, evecs_H = np.linalg.eigh(H)
    boltz = np.exp(-beta * (evals_H - min(evals_H)))
    Z = np.sum(boltz)
    for i in range(d):
        rho += (boltz[i] / Z) * np.outer(evecs_H[:, i], evecs_H[:, i].conj())
    return rho


def spectral_metrics(L, d):
    """Q-factor and frequency count from Liouvillian eigenvalues."""
    evals = np.linalg.eigvals(L)
    osc = []
    for ev in evals:
        rate = -ev.real
        freq = abs(ev.imag)
        if freq > 1e-10 and rate > 1e-15:
            osc.append((rate, freq, freq / rate))
    if not osc:
        return 0, 0, 0
    Qs = [q for _, _, q in osc]
    freqs = sorted(set(round(f, 4) for _, f, _ in osc if f > 1e-8))
    return max(Qs), len(freqs), len([q for q in Qs if q > 5])


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "self_heating_fixpoint.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 70)
    out("SELF-HEATING FIXED POINT")
    out("Does the resonator find its own temperature?")
    out("=" * 70)

    J = 1.0

    configs = [
        ("N=3, pure amp", 3, 0.0, 0.1),
        ("N=3, Z+amp", 3, 0.1, 0.05),
        ("N=5, pure amp", 5, 0.0, 0.1),
        ("N=5, Z+amp", 5, 0.1, 0.05),
        ("N=5, Z+amp (weak)", 5, 0.1, 0.01),
        ("N=5, sacrifice+amp", 5, [0.5, 0.01, 0.01, 0.01, 0.01], 0.05),
    ]

    for label, N, gz, ga in configs:
        out(f"\n{'='*70}")
        out(f"CONFIG: {label}")
        out(f"{'='*70}")

        d = 2**N
        H = build_H(N, J)
        E_ground = min(np.linalg.eigvalsh(H))
        E_max = max(np.linalg.eigvalsh(H))

        out(f"N={N}, J={J}, gamma_z={gz}, gamma_amp={ga}")
        out(f"H spectrum: [{E_ground:.2f}, {E_max:.2f}]")

        # Fixed-point iteration
        out(f"\n{'iter':>4} {'n_bar':>8} {'E_ss':>10} {'E_th':>10}"
            f" {'gap':>10} {'Q_max':>8} {'freqs':>6}")

        n_bar = 0.001  # start slightly above 0 to avoid division issues
        converged = False

        for iteration in range(30):
            # Build Liouvillian at current n_bar
            L = build_L_thermal(H, gz, ga, n_bar)

            # Steady state
            rho_ss = find_steady_state(L, d)
            E_ss = np.real(np.trace(H @ rho_ss))

            # Thermal state at same n_bar
            rho_th = thermal_state(H, n_bar)
            E_th = np.real(np.trace(H @ rho_th))

            # Spectral metrics
            Q_max, n_freq, n_hiQ = spectral_metrics(L, d)

            gap = E_ss - E_th

            out(f"{iteration:4d} {n_bar:8.4f} {E_ss:10.4f} {E_th:10.4f}"
                f" {gap:+10.4f} {Q_max:8.1f} {n_freq:6d}")

            # Check convergence
            if abs(gap) < 0.001 * abs(E_max - E_ground):
                out(f"\n  CONVERGED at n_bar = {n_bar:.4f}")
                out(f"  E_ss = E_th = {E_ss:.4f}")
                out(f"  Q_max = {Q_max:.1f}, {n_freq} frequencies")
                converged = True
                break

            # Update n_bar: if E_ss > E_th, system is hotter, raise n_bar
            # Use bisection-like update
            if gap > 0:
                n_bar *= 1.5  # system hotter than bath: raise bath temp
            else:
                n_bar *= 0.7  # bath hotter than system: lower bath temp

            # Clamp
            n_bar = max(1e-6, min(n_bar, 50.0))

        if not converged:
            out(f"\n  NOT CONVERGED after 30 iterations")
            out(f"  Final n_bar = {n_bar:.4f}, gap = {gap:+.4f}")

    # Summary
    out(f"\n{'='*70}")
    out("SUMMARY")
    out(f"{'='*70}")
    out("\nThe fixed-point iteration finds n_bar* where the steady-state")
    out("energy of the Liouvillian equals the thermal energy at n_bar*.")
    out("At this point, the system is in thermal equilibrium with its")
    out("own decay products.")

    out("\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
