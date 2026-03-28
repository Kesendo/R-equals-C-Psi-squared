#!/usr/bin/env python3
"""
The hydrogen bond as a qubit.

The proton in O-H...O tunnels between two wells: |L> (donor side)
and |R> (acceptor side). This IS a qubit. d=2. The palindrome is
PROVEN. CΨ MUST cross 1/4. We compute the timescales.

Phase 1: Single proton qubit (N=1)
Phase 2: One water molecule (N=2, two proton qubits)
Phase 3: Two water molecules + hydrogen bond (N=4)
"""
import numpy as np
from scipy.linalg import expm
import os

# Physical constants
hbar = 6.582e-16  # eV*s
kB = 8.617e-5     # eV/K
meV = 1e-3        # eV


def build_liouvillian(H, c_ops):
    """Build Liouvillian superoperator for density matrix evolution."""
    d = H.shape[0]
    eye = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, eye) - np.kron(eye, H.T))
    for c in c_ops:
        cd = c.conj().T
        cdc = cd @ c
        L += (np.kron(c, c.conj())
              - 0.5 * np.kron(cdc, eye)
              - 0.5 * np.kron(eye, cdc.T))
    return L


def compute_cpsi(rho):
    """CΨ = Tr(rho^2) * L1(rho) / (d-1)."""
    d = rho.shape[0]
    purity = np.real(np.trace(rho @ rho))
    # L1 coherence: sum of absolute off-diagonal elements
    L1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    return purity * L1 / (d - 1)


def evolve_and_track(L, rho0, times):
    """Evolve density matrix and track CΨ."""
    d2 = L.shape[0]
    d = int(np.sqrt(d2))
    rho_vec = rho0.flatten()
    cpsi_values = []
    for t in times:
        rho_t_vec = expm(L * t) @ rho_vec
        rho_t = rho_t_vec.reshape(d, d)
        cpsi_values.append(compute_cpsi(rho_t))
    return np.array(cpsi_values)


def check_palindrome(L, tol=1e-8):
    """Check palindromic eigenvalue pairing."""
    ev = np.linalg.eigvals(L)
    sum_gamma = -np.min(ev.real)  # steady state eigenvalue ~ 0
    # Find the pairing constant
    rates = sorted(ev.real)
    # Pair: for each eigenvalue, find partner at -ev - 2*sum_gamma...
    # Actually just check all pair sums
    n = len(ev)
    pair_sums = []
    used = set()
    for i in range(n):
        if i in used:
            continue
        target = -ev[i].real  # approximate partner
        best_j, best_d = -1, np.inf
        for j in range(n):
            if j != i and j not in used:
                # Check if ev[i] + ev[j] is constant
                d = abs((ev[i] + ev[j]).real - (ev[0] + ev[1]).real) if len(pair_sums) == 0 else abs((ev[i] + ev[j]).real - pair_sums[0])
                if len(pair_sums) == 0:
                    d = 0  # first pair defines the constant
                if d < best_d:
                    best_j = j
                    best_d = d
        if best_j >= 0:
            pair_sums.append((ev[i] + ev[best_j]).real)
            used.add(i)
            used.add(best_j)

    if pair_sums:
        mean_sum = np.mean(pair_sums)
        std_sum = np.std(pair_sums)
        return mean_sum, std_sum, len(pair_sums)
    return 0, 0, 0


# Pauli matrices
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result


# ================================================================
# PHASE 1: Single Proton Qubit (N=1)
# ================================================================
print("=" * 65)
print("PHASE 1: Single Proton Qubit")
print("The proton in O-H...O tunnels between |L> and |R>.")
print("=" * 65)

# Sweep J/gamma ratio
print(f"\n  {'J (meV)':>8s}  {'gamma':>10s}  {'J/gamma':>8s}  "
      f"{'CΨ(0)':>6s}  {'t_cross':>10s}  {'regime':>10s}")
print(f"  {'-'*60}")

for J_meV in [0.5, 1.0, 5.0, 10.0]:
    for gamma_ratio in [0.01, 0.1, 1.0, 10.0]:
        J = J_meV * meV  # eV
        gamma = J / gamma_ratio if gamma_ratio > 0 else 0.001 * J

        # Hamiltonian: tunneling
        H = -J * sx

        # Dephasing
        c_ops = [np.sqrt(gamma / hbar) * sz]

        # Liouvillian
        L = build_liouvillian(H / hbar, c_ops)  # work in angular frequency

        # Initial state: |L> = |0>
        rho0 = np.array([[1, 0], [0, 0]], dtype=complex)

        # Time evolution
        tau = hbar / J  # natural timescale
        times = np.linspace(0, 50 * tau, 500)
        cpsi = evolve_and_track(L, rho0, times)

        # Find crossing of 1/4
        t_cross = None
        for k in range(len(cpsi) - 1):
            if cpsi[k] > 0.25 and cpsi[k+1] <= 0.25:
                frac = (0.25 - cpsi[k]) / (cpsi[k+1] - cpsi[k])
                t_cross = times[k] + frac * (times[k+1] - times[k])
                break

        jg = gamma_ratio
        regime = "quantum" if jg > 1 else ("fold" if 0.3 < jg < 3 else "classical")
        cross_str = f"{t_cross/1e-12:.2f} ps" if t_cross else "no cross"

        if jg in [0.01, 1.0, 10.0] or (jg == 0.1 and J_meV in [0.5, 5.0]):
            print(f"  {J_meV:8.1f}  {gamma/meV:10.4f}  {jg:8.2f}  "
                  f"{cpsi[0]:6.4f}  {cross_str:>10s}  {regime:>10s}")


# ================================================================
# PHASE 2: One Water Molecule (N=2)
# ================================================================
print("\n" + "=" * 65)
print("PHASE 2: One Water Molecule (2 proton qubits)")
print("H-O-H: two O-H bonds, coupled through shared oxygen")
print("=" * 65)

J_intra = 1.0 * meV  # intramolecular tunneling
K_OO = 0.1 * meV     # proton-proton coupling through O
gamma_val = 1.0 * meV  # dephasing (J/gamma = 1, near fold)

# Hamiltonian: two qubits
H2 = (-J_intra * kron_list([sx, I2])        # qubit 1 tunneling
      - J_intra * kron_list([I2, sx])        # qubit 2 tunneling
      + K_OO * kron_list([sz, sz]))          # coupling through O

# Dephasing on both qubits
c_ops2 = [np.sqrt(gamma_val / hbar) * kron_list([sz, I2]),
          np.sqrt(gamma_val / hbar) * kron_list([I2, sz])]

L2 = build_liouvillian(H2 / hbar, c_ops2)

# Check palindrome
ev2 = np.linalg.eigvals(L2)
n_osc = np.sum(np.abs(ev2.imag) > 1e-3)
pair_sum, pair_std, n_pairs = check_palindrome(L2)

# Count distinct frequencies
freqs = sorted(set(np.round(np.abs(ev2.imag), 4)))
freqs = [f for f in freqs if f > 0.01]

print(f"\n  Parameters: J_intra={J_intra/meV:.1f} meV, K={K_OO/meV:.1f} meV, "
      f"gamma={gamma_val/meV:.1f} meV")
print(f"  Liouvillian: {L2.shape[0]}x{L2.shape[0]}")
print(f"  Oscillatory eigenvalues: {n_osc}")
print(f"  Distinct frequencies: {len(freqs)}")
print(f"  Palindromic pair sums: mean={pair_sum:.6f}, std={pair_std:.2e}")
print(f"  Palindrome: {'EXACT' if pair_std < 1e-6 else 'approximate'}")

# CΨ evolution
rho0_2 = np.zeros((4, 4), dtype=complex)
rho0_2[0, 0] = 1.0  # |LL>
tau2 = hbar / J_intra
times2 = np.linspace(0, 100 * tau2, 1000)
cpsi2 = evolve_and_track(L2, rho0_2, times2)

t_cross2 = None
for k in range(len(cpsi2) - 1):
    if cpsi2[k] > 0.25 and cpsi2[k+1] <= 0.25:
        frac = (0.25 - cpsi2[k]) / (cpsi2[k+1] - cpsi2[k])
        t_cross2 = times2[k] + frac * (times2[k+1] - times2[k])
        break

print(f"  CΨ(0) = {cpsi2[0]:.4f}")
print(f"  CΨ crosses 1/4: {'t = ' + f'{t_cross2/1e-12:.2f} ps' if t_cross2 else 'no crossing'}")


# ================================================================
# PHASE 3: Two Water Molecules + Hydrogen Bond (N=4)
# ================================================================
print("\n" + "=" * 65)
print("PHASE 3: Two Water Molecules + H-bond (4 proton qubits)")
print("H-O-H ... O-H-H: mediator topology")
print("=" * 65)

J_intra = 1.0 * meV    # intramolecular tunneling (strong)
J_inter = 0.1 * meV     # intermolecular tunneling through H-bond (weak)
K_intra = 0.1 * meV     # intramolecular proton-proton coupling
gamma_val = 1.0 * meV   # dephasing

# 4 qubits: [q1, q2] - molecule 1, [q3, q4] - molecule 2
# H-bond: q2 (donor proton) couples to q3 (acceptor side)
d4 = 16
I4 = [I2, I2, I2, I2]

def op4(pauli, site):
    """4-qubit operator: pauli on site, identity elsewhere."""
    ops = [I2, I2, I2, I2]
    ops[site] = pauli
    return kron_list(ops)

def op4_2(p1, s1, p2, s2):
    """4-qubit two-body operator."""
    ops = [I2, I2, I2, I2]
    ops[s1] = p1
    ops[s2] = p2
    return kron_list(ops)

H4 = (# Intramolecular tunneling
      - J_intra * op4(sx, 0)
      - J_intra * op4(sx, 1)
      - J_intra * op4(sx, 2)
      - J_intra * op4(sx, 3)
      # Intramolecular coupling (within each molecule)
      + K_intra * op4_2(sz, 0, sz, 1)
      + K_intra * op4_2(sz, 2, sz, 3)
      # INTERMOLECULAR: H-bond coupling q2-q3
      + J_inter * op4_2(sx, 1, sx, 2))

# Dephasing on all 4 qubits
c_ops4 = [np.sqrt(gamma_val / hbar) * op4(sz, i) for i in range(4)]

L4 = build_liouvillian(H4 / hbar, c_ops4)

ev4 = np.linalg.eigvals(L4)
n_osc4 = np.sum(np.abs(ev4.imag) > 1e-3)
freqs4 = sorted(set(np.round(np.abs(ev4.imag), 2)))
freqs4 = [f for f in freqs4 if f > 0.1]
pair_sum4, pair_std4, n_pairs4 = check_palindrome(L4)

print(f"\n  Parameters: J_intra={J_intra/meV:.1f}, J_inter={J_inter/meV:.1f}, "
      f"K={K_intra/meV:.1f}, gamma={gamma_val/meV:.1f} meV")
print(f"  Liouvillian: {L4.shape[0]}x{L4.shape[0]}")
print(f"  Oscillatory eigenvalues: {n_osc4}")
print(f"  Distinct frequencies: {len(freqs4)}")
print(f"  Palindromic pair sums: mean={pair_sum4:.4f}, std={pair_std4:.2e}")
print(f"  Palindrome: {'EXACT' if pair_std4 < 1e-4 else 'approximate'}")

# V-Effect: compare with single molecule
print(f"\n  V-Effect check:")
print(f"    Single molecule (N=2): {len(freqs)} distinct frequencies")
print(f"    Two molecules (N=4):   {len(freqs4)} distinct frequencies")
print(f"    Ratio: {len(freqs4)} / (2 x {len(freqs)}) = "
      f"{len(freqs4) / (2 * len(freqs)) if len(freqs) > 0 else 'inf':.2f}")
if len(freqs4) > 2 * len(freqs):
    print(f"    >>> V-EFFECT: {len(freqs4) - 2*len(freqs)} new frequencies!")

# CΨ across the hydrogen bond (qubits 2 and 3)
rho0_4 = np.zeros((d4, d4), dtype=complex)
rho0_4[0, 0] = 1.0  # |LLLL>
tau4 = hbar / J_intra
times4 = np.linspace(0, 100 * tau4, 1000)
cpsi4 = evolve_and_track(L4, rho0_4, times4)

# Count CΨ crossings of 1/4
crossings = 0
for k in range(len(cpsi4) - 1):
    if (cpsi4[k] - 0.25) * (cpsi4[k+1] - 0.25) < 0:
        crossings += 1

print(f"\n  CΨ(0) = {cpsi4[0]:.4f}")
print(f"  CΨ crossings of 1/4: {crossings}")
if crossings > 0:
    # Find first crossing time
    for k in range(len(cpsi4) - 1):
        if cpsi4[k] > 0.25 and cpsi4[k+1] <= 0.25:
            frac = (0.25 - cpsi4[k]) / (cpsi4[k+1] - cpsi4[k])
            t_first = times4[k] + frac * (times4[k+1] - times4[k])
            print(f"  First crossing: t = {t_first/1e-12:.4f} ps")
            break
    print(f"  Q-factor (crossings/2): {crossings // 2}")


# ================================================================
# Summary
# ================================================================
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"""
The proton in a hydrogen bond IS a qubit (d=2).
The palindrome MUST hold (proven for d=2 with Z-dephasing).
CΨ MUST cross 1/4 (proven for all Markovian channels).

The question was not WHETHER but WHEN and HOW FAST.
""")
