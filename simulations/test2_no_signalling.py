#!/usr/bin/env python3
"""
No-signalling with a separate scalar quarter-band readout.

For a Bell+ pair, B is dephased in the Z basis and outcomes are averaged.
A's reduced state rho_A is unchanged. The separately defined scalar
Tr(rho_AB^2) times the largest eigenvalue of rho_A moves from 0.500 to 0.250.

Script:  simulations/test2_no_signalling.py
Output:  simulations/results/test2_no_signalling.txt
Docs:    experiments/NO_SIGNALLING_BOUNDARY.md
"""

import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


OUT_PATH = Path(__file__).parent / "results" / "test2_no_signalling.txt"
_lines = None


def log(msg=""):
    if _lines is None:
        raise RuntimeError("log() is available only during main().")
    print(msg, flush=True)
    _lines.append(str(msg))


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)

N = 2
d = 4
d2 = 16


def site_op(op, k, nq=N):
    ops = [I2] * nq
    ops[k] = op
    result = ops[0]
    for other in ops[1:]:
        result = np.kron(result, other)
    return result


def build_H(J=1.0):
    hamiltonian = np.zeros((d, d), dtype=complex)
    for pauli in (sx, sy, sz):
        hamiltonian += J * site_op(pauli, 0) @ site_op(pauli, 1)
    return hamiltonian


def build_L(hamiltonian, gamma):
    identity = np.eye(d)
    generator = -1j * (
        np.kron(hamiltonian, identity) - np.kron(identity, hamiltonian.T)
    )
    for k in range(N):
        z_k = site_op(sz, k)
        generator += gamma * (np.kron(z_k, z_k.conj()) - np.eye(d2))
    return generator


def evolve(generator, rho, time):
    vector = expm(generator * time) @ rho.flatten()
    evolved = vector.reshape(d, d)
    return (evolved + evolved.conj().T) / 2


def ptrace_A(rho):
    """Partial trace over B, keeping A."""
    reshaped = rho.reshape(2, 2, 2, 2)
    return np.trace(reshaped, axis1=1, axis2=3)


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


def von_neumann_entropy(rho):
    eigenvalues = np.real(np.linalg.eigvalsh(rho))
    eigenvalues = eigenvalues[eigenvalues > 1e-15]
    return float(-np.sum(eigenvalues * np.log2(eigenvalues)))


def apply_B_measurement_Z(rho):
    """Apply the outcome-averaged local Z measurement channel on B."""
    p0 = site_op(np.outer(up, up.conj()), 1)
    p1 = site_op(np.outer(dn, dn.conj()), 1)
    return p0 @ rho @ p0.conj().T + p1 @ rho @ p1.conj().T


def cpsi(rho_AB):
    """Separately defined scalar for the fixed A|B subsystem split."""
    rho_A = ptrace_A(rho_AB)
    psi_readout = np.max(np.real(np.linalg.eigvalsh(rho_A)))
    return purity(rho_AB) * psi_readout


def quarter_band(value, tolerance=1e-4):
    if value > 0.25 + tolerance:
        return "above 1/4"
    if value < 0.25 - tolerance:
        return "below 1/4"
    return "equal to 1/4"


def no_signalling_reading(rho_A_before, rho_A_after, tolerance=1e-12):
    """Return the reduced-state distance and its tolerance-scoped equality."""
    delta = float(np.linalg.norm(rho_A_before - rho_A_after))
    unchanged = bool(np.allclose(
        rho_A_before, rho_A_after, atol=tolerance, rtol=0.0
    ))
    return delta, unchanged


def main(output_path=OUT_PATH):
    global _lines
    if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    _lines = []

    log("QUARTER-CURRENT")
    log("Current reading: reduced-state equality carries no-signalling; the scalar bin is separate.")
    log()
    log("=" * 70)
    log("No-signalling and a separate scalar quarter-band readout")
    log("=" * 70)
    log()

    bell_plus = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho_before = ket2dm(bell_plus)
    rho_after = apply_B_measurement_Z(rho_before)

    rho_A_before = ptrace_A(rho_before)
    rho_A_after = ptrace_A(rho_after)
    delta, no_signalling = no_signalling_reading(rho_A_before, rho_A_after)

    log("Test 1: Bell+ pair, averaged local Z measurement channel on B")
    log("-" * 70)
    log(f"  ||delta rho_A|| = {delta:.10f}")
    log()
    log(f"  {'Quantity':>20}  {'Before':>10}  {'After':>10}  {'Changed?':>9}")
    log("  " + "-" * 55)

    for name, pauli in (("<sigma_x>_A", sx), ("<sigma_y>_A", sy), ("<sigma_z>_A", sz)):
        before = float(np.real(np.trace(rho_A_before @ pauli)))
        after = float(np.real(np.trace(rho_A_after @ pauli)))
        changed = abs(before - after) > 1e-12
        log(f"  {name:>20}  {before:10.3f}  {after:10.3f}  {'YES' if changed else 'NO':>9}")

    purity_A_before = purity(rho_A_before)
    purity_A_after = purity(rho_A_after)
    log(f"  {'Purity(rho_A)':>20}  {purity_A_before:10.3f}  {purity_A_after:10.3f}  "
        f"{'YES' if abs(purity_A_before - purity_A_after) > 1e-12 else 'NO':>9}")

    entropy_before = von_neumann_entropy(rho_A_before)
    entropy_after = von_neumann_entropy(rho_A_after)
    log(f"  {'S(rho_A)':>20}  {entropy_before:10.3f}  {entropy_after:10.3f}  "
        f"{'YES' if abs(entropy_before - entropy_after) > 1e-12 else 'NO':>9}")

    psi_before = float(np.max(np.real(np.linalg.eigvalsh(rho_A_before))))
    psi_after = float(np.max(np.real(np.linalg.eigvalsh(rho_A_after))))
    global_purity_before = purity(rho_before)
    global_purity_after = purity(rho_after)
    scalar_before = global_purity_before * psi_before
    scalar_after = global_purity_after * psi_after

    log(f"  {'max eig(rho_A)':>20}  {psi_before:10.3f}  {psi_after:10.3f}  "
        f"{'YES' if abs(psi_before - psi_after) > 1e-12 else 'NO':>9}")
    log(f"  {'Tr(rho_AB^2)':>20}  {global_purity_before:10.3f}  {global_purity_after:10.3f}  "
        f"{'YES' if abs(global_purity_before - global_purity_after) > 1e-12 else 'NO':>9}")
    log(f"  {'scalar readout':>20}  {scalar_before:10.3f}  {scalar_after:10.3f}  "
        f"{'YES' if abs(scalar_before - scalar_after) > 1e-12 else 'NO':>9}")
    log(f"  {'quarter band':>20}  {quarter_band(scalar_before):>10}  "
        f"{quarter_band(scalar_after):>10}  {'YES':>9}")

    if not no_signalling:
        raise AssertionError(f"rho_A changed: ||delta rho_A||={delta:.3e}")
    log()
    log("No-signalling verdict: rho_A is unchanged under the averaged local channel on B.")
    log(f"Global purity: {global_purity_before:.3f} -> {global_purity_after:.3f}")
    log(f"Scalar readout: {scalar_before:.3f} -> {scalar_after:.3f} "
        f"({quarter_band(scalar_before)} -> {quarter_band(scalar_after)}).")

    log()
    log("Test 2: Bell+ under local Z dephasing; apply the same channel at each sampled time")
    log("-" * 70)
    gamma = 0.05
    hamiltonian = build_H(1.0)
    generator = build_L(hamiltonian, gamma)
    rho0 = ket2dm(bell_plus)
    log(f"  {'t':>6}  {'scalar(no)':>12}  {'scalar(B->Z)':>12}  "
        f"{'||delta rho_A||':>16}  {'band(no)':>14}  {'band(Z)':>14}")
    for time in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0):
        rho_t = evolve(generator, rho0, time)
        rho_bz = apply_B_measurement_Z(rho_t)
        local_delta = np.linalg.norm(ptrace_A(rho_t) - ptrace_A(rho_bz))
        scalar_no = cpsi(rho_t)
        scalar_bz = cpsi(rho_bz)
        log(f"  {time:6.1f}  {scalar_no:12.4f}  {scalar_bz:12.4f}  "
            f"{local_delta:16.2e}  {quarter_band(scalar_no):>14}  "
            f"{quarter_band(scalar_bz):>14}")

    log()
    log("Test 3: finite scalar trajectories after applying the channel at t=2")
    log("-" * 70)
    rho_t2 = evolve(generator, rho0, 2.0)
    rho_t2_bz = apply_B_measurement_Z(rho_t2)
    log(f"  {'t':>6}  {'scalar(no)':>12}  {'scalar(channel@2)':>18}  {'difference':>12}")
    for after in (0.0, 0.5, 1.0, 2.0, 5.0):
        scalar_no = cpsi(evolve(generator, rho_t2, after))
        scalar_bz = cpsi(evolve(generator, rho_t2_bz, after))
        log(f"  {2.0 + after:6.1f}  {scalar_no:12.4f}  {scalar_bz:18.4f}  "
            f"{scalar_bz - scalar_no:12.4f}")

    log()
    log("=" * 70)
    log("The exact reduced-state equality carries the no-signalling statement.")
    log("The scalar quarter-band change is not a physical boundary.")
    log("Completed: deterministic rerun")
    log("Results: simulations/results/test2_no_signalling.txt")

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(("\n".join(_lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=OUT_PATH)
    main(output_path=parser.parse_args().output)
