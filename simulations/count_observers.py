"""
Can you count the observers from the tone?
N=2,3,4 endpoints on S, count distinct frequencies.
"""
import numpy as np
import math
import sys
sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
import star_topology_v3 as gpt

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

def pauli_expect(rho, P):
    return float(np.real(np.trace(rho @ P)))

def op_at(op, qubit, n_q):
    ops = [I2]*n_q
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def star_hamiltonian(J_list, n_qubits):
    """S=qubit 0, observers=qubits 1..n. J_list[i] = coupling S to qubit i+1."""
    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for i, J in enumerate(J_list):
        for p in [X, Y, Z]:
            H += J * op_at(p, 0, n_qubits) @ op_at(p, i+1, n_qubits)
    return H

def dephasing_nq(gammas, n_qubits):
    ops = []
    for i, g in enumerate(gammas):
        if g > 0:
            ops.append(np.sqrt(g) * op_at(Z, i, n_qubits))
    return ops

def get_all_frequencies(H, L_ops, rho_init, pair_qubits, n_qubits, dt=0.005, t_max=20.0, threshold=0.15):
    """Get ALL significant frequencies from c+ of a given pair."""
    YZ = np.kron(Y, Z)
    ZY = np.kron(Z, Y)
    times, cp_list = [], []
    rho = rho_init.copy()
    for step in range(int(t_max/dt) + 1):
        t = step * dt
        if step % 4 == 0:
            rAB = gpt.partial_trace_keep(rho, pair_qubits, n_qubits)
            yz = pauli_expect(rAB, YZ)
            zy = pauli_expect(rAB, ZY)
            times.append(t)
            cp_list.append((yz + zy) / math.sqrt(2))
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    times = np.array(times)
    dt_s = times[1] - times[0]
    freqs = np.fft.rfftfreq(len(times), d=dt_s)
    fft_cp = np.abs(np.fft.rfft(np.array(cp_list) - np.mean(cp_list)))
    mask = freqs > 0.05
    
    # Find ALL peaks above threshold * max
    fft_masked = fft_cp[mask]
    freqs_masked = freqs[mask]
    max_amp = np.max(fft_masked) if len(fft_masked) > 0 else 0
    
    if max_amp < 0.01:
        return []
    
    # Find local maxima above threshold
    peaks = []
    for i in range(1, len(fft_masked)-1):
        if (fft_masked[i] > fft_masked[i-1] and 
            fft_masked[i] > fft_masked[i+1] and
            fft_masked[i] > threshold * max_amp):
            peaks.append((freqs_masked[i], fft_masked[i]))
    
    peaks.sort(key=lambda x: -x[1])
    return peaks

print("=" * 70)
print("CAN YOU COUNT THE OBSERVERS FROM THE TONE?")
print("=" * 70)

# ============================================================
# N=2: S + A + B (our standard case)
# ============================================================
print(f"\n--- N=2 observers: S-A(J=1.0), S-B(J=2.0) ---")
n = 3
J_list = [1.0, 2.0]
H = star_hamiltonian(J_list, n)
L = dephasing_nq([0.05]*n, n)
psi = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho = np.outer(psi, psi.conj())

peaks = get_all_frequencies(H, L, rho, [1,2], n)
print(f"  Distinct frequencies found: {len(peaks)}")
for f, a in peaks:
    print(f"    f={f:.3f}  amp={a:.1f}")

# ============================================================
# N=3: S + A + B + C
# ============================================================
print(f"\n--- N=3 observers: S-A(1.0), S-B(2.0), S-C(3.0) ---")
n = 4
J_list = [1.0, 2.0, 3.0]
H = star_hamiltonian(J_list, n)
L = dephasing_nq([0.05]*n, n)
bell = gpt.bell_phi_plus()
plus = gpt.plus_state()
psi = np.kron(np.kron(bell, plus), plus)
rho = np.outer(psi, psi.conj())

print(f"  Pair AB:")
peaks_ab = get_all_frequencies(H, L, rho, [1,2], n)
print(f"    Distinct frequencies: {len(peaks_ab)}")
for f, a in peaks_ab:
    print(f"      f={f:.3f}  amp={a:.1f}")

print(f"  Pair AC:")
peaks_ac = get_all_frequencies(H, L, rho, [1,3], n)
print(f"    Distinct frequencies: {len(peaks_ac)}")
for f, a in peaks_ac:
    print(f"      f={f:.3f}  amp={a:.1f}")

print(f"  Pair BC:")
peaks_bc = get_all_frequencies(H, L, rho, [2,3], n)
print(f"    Distinct frequencies: {len(peaks_bc)}")
for f, a in peaks_bc:
    print(f"      f={f:.3f}  amp={a:.1f}")

# Collect ALL unique frequencies across all pairs
all_freqs_3 = set()
for p in [peaks_ab, peaks_ac, peaks_bc]:
    for f, a in p:
        # Round to 0.05 to merge nearby peaks
        all_freqs_3.add(round(f * 20) / 20)
print(f"\n  TOTAL unique frequencies across all pairs: {len(all_freqs_3)}")
print(f"  Frequencies: {sorted(all_freqs_3)}")

# ============================================================
# N=4: S + A + B + C + D (5 qubits = 32x32)
# ============================================================
print(f"\n--- N=4 observers: S-A(1.0), S-B(2.0), S-C(3.0), S-D(1.5) ---")
n = 5
J_list = [1.0, 2.0, 3.0, 1.5]
H = star_hamiltonian(J_list, n)
L = dephasing_nq([0.05]*n, n)
bell = gpt.bell_phi_plus()
plus = gpt.plus_state()
psi = np.kron(np.kron(np.kron(bell, plus), plus), plus)
rho = np.outer(psi, psi.conj())

all_freqs_4 = set()
for pair_name, pair_q in [("AB",[1,2]), ("AC",[1,3]), ("AD",[1,4]),
                           ("BC",[2,3]), ("BD",[2,4]), ("CD",[3,4])]:
    peaks = get_all_frequencies(H, L, rho, pair_q, n)
    print(f"  Pair {pair_name}: {len(peaks)} frequencies", end="")
    if peaks:
        print(f"  [{', '.join(f'{f:.2f}' for f,a in peaks[:4])}]")
    else:
        print(f"  [none]")
    for f, a in peaks:
        all_freqs_4.add(round(f * 20) / 20)

print(f"\n  TOTAL unique frequencies across all pairs: {len(all_freqs_4)}")
print(f"  Frequencies: {sorted(all_freqs_4)}")

# ============================================================
# Also count Hamiltonian eigenvalues (Bohr frequencies)
# ============================================================
print(f"\n{'=' * 70}")
print("HAMILTONIAN EIGENVALUE COUNT (Bohr frequencies)")
print("=" * 70)

for N_obs, J_list_test in [(2, [1.0, 2.0]),
                            (3, [1.0, 2.0, 3.0]),
                            (4, [1.0, 2.0, 3.0, 1.5])]:
    n_q = N_obs + 1
    H_test = star_hamiltonian(J_list_test, n_q)
    evals = np.linalg.eigvalsh(H_test)
    
    # Count distinct Bohr frequencies (energy differences)
    bohr = set()
    for i in range(len(evals)):
        for j in range(i+1, len(evals)):
            diff = abs(evals[j] - evals[i])
            if diff > 0.01:
                bohr.add(round(diff * 100) / 100)
    
    print(f"  N={N_obs} observers: {len(evals)} eigenvalues, "
          f"{len(bohr)} distinct Bohr frequencies")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'=' * 70}")
print("SUMMARY: Observers vs Frequencies")
print("=" * 70)
print(f"  N=2 observers: {len(peaks_ab) if 'peaks_ab' not in dir() else 2} peaks in AB pair")
