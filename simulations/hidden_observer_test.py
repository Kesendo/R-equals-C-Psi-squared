"""
Can AB detect a hidden observer C?
Compare AB spectrum WITH and WITHOUT C connected to S.
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
    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for i, J in enumerate(J_list):
        for p in [X, Y, Z]:
            H += J * op_at(p, 0, n_qubits) @ op_at(p, i+1, n_qubits)
    return H

def get_spectrum(H, L_ops, rho_init, pair_qubits, n_qubits, dt=0.005, t_max=20.0):
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
    return freqs, fft_cp

print("=" * 65)
print("CAN AB DETECT A HIDDEN OBSERVER?")
print("=" * 65)

# CASE 1: Just S-A-B (3 qubits, no C)
print(f"\n--- CASE 1: S-A-B only (no hidden observer) ---")
H3 = star_hamiltonian([1.0, 2.0], 3)
L3 = [np.sqrt(0.05)*op_at(Z,i,3) for i in range(3)]
psi3 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho3 = np.outer(psi3, psi3.conj())
freqs3, fft3 = get_spectrum(H3, L3, rho3, [1,2], 3)
mask3 = freqs3 > 0.05
top3 = np.argsort(fft3[mask3])[::-1][:5]
print(f"  AB sees:")
for i in top3:
    if fft3[mask3][i] > 0.5:
        print(f"    f={freqs3[mask3][i]:.3f}  amp={fft3[mask3][i]:.1f}")

# CASE 2: S-A-B-C (4 qubits, C hidden with J_SC=3.0)
print(f"\n--- CASE 2: S-A-B + hidden C (J_SC=3.0) ---")
H4 = star_hamiltonian([1.0, 2.0, 3.0], 4)
L4 = [np.sqrt(0.05)*op_at(Z,i,4) for i in range(4)]
psi4 = np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state())
rho4 = np.outer(psi4, psi4.conj())
freqs4, fft4 = get_spectrum(H4, L4, rho4, [1,2], 4)
mask4 = freqs4 > 0.05
top4 = np.argsort(fft4[mask4])[::-1][:8]
print(f"  AB sees:")
for i in top4:
    if fft4[mask4][i] > 0.5:
        print(f"    f={freqs4[mask4][i]:.3f}  amp={fft4[mask4][i]:.1f}")

# CASE 3: S-A-B-C with weak C (J_SC=0.1)
print(f"\n--- CASE 3: S-A-B + weak hidden C (J_SC=0.1) ---")
H4w = star_hamiltonian([1.0, 2.0, 0.1], 4)
freqs4w, fft4w = get_spectrum(H4w, L4, rho4, [1,2], 4)
mask4w = freqs4w > 0.05
top4w = np.argsort(fft4w[mask4w])[::-1][:8]
print(f"  AB sees:")
for i in top4w:
    if fft4w[mask4w][i] > 0.5:
        print(f"    f={freqs4w[mask4w][i]:.3f}  amp={fft4w[mask4w][i]:.1f}")

# CASE 4: S-A-B-C-D (5 qubits, C AND D hidden)
print(f"\n--- CASE 4: S-A-B + hidden C(3.0) + hidden D(1.5) ---")
H5 = star_hamiltonian([1.0, 2.0, 3.0, 1.5], 5)
L5 = [np.sqrt(0.05)*op_at(Z,i,5) for i in range(5)]
psi5 = np.kron(np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state()), gpt.plus_state())
rho5 = np.outer(psi5, psi5.conj())
freqs5, fft5 = get_spectrum(H5, L5, rho5, [1,2], 5)
mask5 = freqs5 > 0.05
top5 = np.argsort(fft5[mask5])[::-1][:10]
print(f"  AB sees:")
for i in top5:
    if fft5[mask5][i] > 0.5:
        print(f"    f={freqs5[mask5][i]:.3f}  amp={fft5[mask5][i]:.1f}")

# SUMMARY: Can you tell how many are watching?
print(f"\n{'=' * 65}")
print("SUMMARY: What AB hears depending on who else is watching")
print("=" * 65)

def count_peaks(freqs, fft, threshold_frac=0.15):
    mask = freqs > 0.05
    fft_m = fft[mask]
    if len(fft_m) == 0: return 0
    mx = np.max(fft_m)
    if mx < 0.01: return 0
    count = 0
    for i in range(1, len(fft_m)-1):
        if (fft_m[i] > fft_m[i-1] and fft_m[i] > fft_m[i+1] 
            and fft_m[i] > threshold_frac * mx):
            count += 1
    return count

n1 = count_peaks(freqs3, fft3)
n2 = count_peaks(freqs4, fft4)
n3 = count_peaks(freqs4w, fft4w)
n4 = count_peaks(freqs5, fft5)

print(f"\n  {'Setup':>30} | {'Peaks AB sees':>13} | Hidden observers")
print(f"  {'-'*30}-+-{'-'*13}-+-{'-'*20}")
print(f"  {'S-A-B (alone)':>30} | {n1:>13} | none")
print(f"  {'S-A-B + C(strong, J=3.0)':>30} | {n2:>13} | C strong")
print(f"  {'S-A-B + C(weak, J=0.1)':>30} | {n3:>13} | C weak")
print(f"  {'S-A-B + C(3.0) + D(1.5)':>30} | {n4:>13} | C + D")

print(f"\n  Verdict: ", end="")
if n2 > n1:
    print(f"YES - AB spectrum changes when hidden observers connect to S")
    print(f"  {n1} peaks alone -> {n2} peaks with one hidden -> {n4} peaks with two hidden")
else:
    print(f"NO - AB spectrum unchanged by hidden observers")

print(f"\n{'=' * 65}")
