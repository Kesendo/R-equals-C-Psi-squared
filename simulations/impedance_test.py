#!/usr/bin/env python3
"""
Phase 2: Impedance Test - Does ||Z rho Z - rho|| peak at CΨ = 1/4?

The Lindblad dissipator acts as: gamma * (Z rho Z - rho).
The magnitude of this term depends on rho. If it peaks at CΨ = 1/4,
then the fold catastrophe is literally the point of maximum coupling
between the system and the external signal.

Usage: python impedance_test.py
"""
import numpy as np
from scipy.linalg import expm
import warnings
warnings.filterwarnings('ignore')

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def tensor(*ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result


def partial_trace(rho, n_qubits, keep):
    """Trace out all qubits except those in keep list."""
    d = 2 ** n_qubits
    d_keep = 2 ** len(keep)
    d_trace = d // d_keep
    rho_reduced = np.zeros((d_keep, d_keep), dtype=complex)

    trace_qubits = [q for q in range(n_qubits) if q not in keep]

    for i in range(d_keep):
        for j in range(d_keep):
            for k in range(d_trace):
                # Build full indices
                row = 0
                col = 0
                ki, kj = i, j
                kk = k
                for q in range(n_qubits - 1, -1, -1):
                    if q in keep:
                        idx = keep.index(q)
                        bit_i = (ki >> (len(keep) - 1 - idx)) & 1
                        bit_j = (kj >> (len(keep) - 1 - idx)) & 1
                    else:
                        tidx = trace_qubits.index(q)
                        bit_i = (kk >> (len(trace_qubits) - 1 - tidx)) & 1
                        bit_j = bit_i  # trace: same index
                    row |= (bit_i << q)
                    col |= (bit_j << q)
                rho_reduced[i, j] += rho[row, col]

    return rho_reduced


def cpsi(rho):
    """Compute CΨ = Tr(rho^2) * L1/(d-1)."""
    d = rho.shape[0]
    purity = np.real(np.trace(rho @ rho))
    l1 = np.sum(np.abs(rho)) - np.real(np.trace(np.abs(rho)))
    if d <= 1:
        return 0
    return purity * l1 / (d - 1)


def impedance(rho, n_qubits):
    """
    Compute ||Z_i rho Z_i - rho||_F for each qubit, summed.
    This is the Frobenius norm of the dissipator action.
    """
    d = 2 ** n_qubits
    total = 0
    for q in range(n_qubits):
        # Build Z_q (Z on qubit q, I on others)
        ops = [I2] * n_qubits
        ops[q] = Z
        Zq = tensor(*ops)
        diff = Zq @ rho @ Zq - rho
        total += np.linalg.norm(diff, 'fro')
    return total / np.linalg.norm(rho, 'fro')


def build_hamiltonian(n_qubits, J=1.0):
    """Heisenberg chain Hamiltonian."""
    d = 2 ** n_qubits
    H = np.zeros((d, d), dtype=complex)
    for i in range(n_qubits - 1):
        for P in [X, Y, Z]:
            ops = [I2] * n_qubits
            ops[i] = P
            ops[i + 1] = P
            H += J * tensor(*ops)
    return H


def lindblad_step(rho, H, gammas, dt, n_qubits):
    """One Euler step of Lindblad evolution."""
    d = 2 ** n_qubits
    # Hamiltonian part
    drho = -1j * (H @ rho - rho @ H)
    # Dissipator
    for q in range(n_qubits):
        ops = [I2] * n_qubits
        ops[q] = Z
        Zq = tensor(*ops)
        drho += gammas[q] * (Zq @ rho @ Zq - rho)
    rho_new = rho + dt * drho
    # Enforce Hermiticity
    rho_new = (rho_new + rho_new.conj().T) / 2
    return rho_new


print("=" * 60)
print("IMPEDANCE TEST: Does ||Z rho Z - rho|| peak at CΨ = 1/4?")
print("=" * 60)

# Test 1: N=2 Bell pair under uniform dephasing
print("\n### N=2, Bell pair, uniform gamma=0.05")

n = 2
d = 4
H = build_hamiltonian(n, J=1.0)
gammas = [0.05, 0.05]

# Bell state
psi = np.zeros(d, dtype=complex)
psi[0] = 1 / np.sqrt(2)
psi[3] = 1 / np.sqrt(2)
rho = np.outer(psi, psi.conj())

dt = 0.01
print(f"\n{'T':>6s} {'CΨ':>8s} {'Impedance':>10s} {'D=1-4CΨ':>8s}")

for step in range(2000):
    t = step * dt
    if step % 20 == 0:
        c = cpsi(rho)
        imp = impedance(rho, n)
        disc = 1 - 4 * c
        marker = " <<<" if abs(c - 0.25) < 0.02 else ""
        print(f"{t:6.2f} {c:8.4f} {imp:10.4f} {disc:8.4f}{marker}")
    rho = lindblad_step(rho, H, gammas, dt, n)

# Test 2: N=3 with structured bath (the heartbeat config)
print("\n\n### N=3, Bell(0,1)+|+>_bath, J=5.0, gamma=[0.0001, 0.0001, 0.005]")

n = 3
d = 8
H = build_hamiltonian(n, J=5.0)
gammas = [0.0001, 0.0001, 0.005]

# Bell(0,1) x |+>_2
psi = np.zeros(d, dtype=complex)
norm = 1 / (np.sqrt(2) * np.sqrt(2))
for b in range(2):  # bath qubit in |+>
    psi[0 * 4 + 0 * 2 + b] += norm  # |00b>
    psi[1 * 4 + 1 * 2 + b] += norm  # |11b>
rho = np.outer(psi, psi.conj())

dt = 0.01
print(f"\n{'T':>6s} {'CΨ(01)':>8s} {'Imp(01)':>10s} {'Imp(full)':>10s} {'D':>8s}")

prev_cpsi = 1.0
crossings = 0

for step in range(3000):
    t = step * dt
    if step % 50 == 0:
        rho01 = partial_trace(rho, n, [0, 1])
        c01 = cpsi(rho01)

        # Impedance of subsystem (qubits 0,1 only)
        imp_sub = 0
        for q in [0, 1]:
            ops = [I2] * n
            ops[q] = Z
            Zq = tensor(*ops)
            diff = Zq @ rho @ Zq - rho
            imp_sub += np.linalg.norm(diff, 'fro')
        imp_sub /= np.linalg.norm(rho, 'fro')

        imp_full = impedance(rho, n)
        disc = 1 - 4 * c01

        marker = ""
        if prev_cpsi >= 0.25 and c01 < 0.25:
            crossings += 1
            marker = " ↓"
        elif prev_cpsi < 0.25 and c01 >= 0.25:
            crossings += 1
            marker = " ↑"

        if marker or step % 200 == 0 or step < 200:
            print(f"{t:6.2f} {c01:8.4f} {imp_sub:10.4f} {imp_full:10.4f} {disc:8.4f}{marker}")

        prev_cpsi = c01

    rho = lindblad_step(rho, H, gammas, dt, n)

print(f"\nTotal crossings: {crossings}")

# Test 3: Impedance vs CΨ scatter (N=2, many time points)
print("\n\n### Impedance vs CΨ scatter (N=2)")
print("Collecting 200 data points across CΨ trajectory...")

n = 2
d = 4
H = build_hamiltonian(n, J=1.0)
gammas = [0.05, 0.05]

psi = np.zeros(d, dtype=complex)
psi[0] = 1 / np.sqrt(2)
psi[3] = 1 / np.sqrt(2)
rho = np.outer(psi, psi.conj())

cpsi_vals = []
imp_vals = []

for step in range(2000):
    if step % 10 == 0:
        c = cpsi(rho)
        imp = impedance(rho, n)
        cpsi_vals.append(c)
        imp_vals.append(imp)
    rho = lindblad_step(rho, H, gammas, dt, n)

# Find where impedance peaks
max_imp_idx = np.argmax(imp_vals)
cpsi_at_max = cpsi_vals[max_imp_idx]
print(f"Max impedance at CΨ = {cpsi_at_max:.4f} (distance from 1/4: {abs(cpsi_at_max - 0.25):.4f})")
print(f"Impedance at max: {imp_vals[max_imp_idx]:.4f}")

# Bin by CΨ range
bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.6), (0.6, 1.0)]
print(f"\n{'CΨ range':>12s} {'Mean Imp':>10s} {'Count':>6s}")
for lo, hi in bins:
    mask = [(lo <= c < hi) for c in cpsi_vals]
    if any(mask):
        vals = [imp_vals[i] for i, m in enumerate(mask) if m]
        print(f"  {lo:.1f}-{hi:.1f}     {np.mean(vals):10.4f} {len(vals):6d}")

print("\nDone.")
