"""
QUANTUM SONAR: Detection threshold for hidden observers.
At what coupling strength does a hidden C become visible in AB's spectrum?
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

def get_spectrum_vector(H, L_ops, rho_init, pair_qubits, n_qubits, dt=0.005, t_max=20.0):
    """Return full FFT spectrum of c+ for a given pair."""
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

# Baseline: AB alone (3 qubits)
print("=" * 65)
print("QUANTUM SONAR: Detection threshold for hidden observers")
print("=" * 65)

print("\nComputing baseline (AB alone)...")
H_base = star_hamiltonian([1.0, 2.0], 3)
L_base = [np.sqrt(0.05)*op_at(Z,i,3) for i in range(3)]
psi_base = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho_base = np.outer(psi_base, psi_base.conj())
freqs_base, fft_base = get_spectrum_vector(H_base, L_base, rho_base, [1,2], 3)

# Sweep J_SC from 0 to 5 and measure spectral distance from baseline
print(f"\nSweeping J_SC (hidden observer coupling strength)...")
print(f"\n{'J_SC':>5} | {'Spectral dist':>13} | {'New peaks':>9} | {'Dominant f':>10} | Detectable?")
print("-" * 65)

# Need to interpolate baseline to 4-qubit frequency grid
# Actually easier: compute spectral distance as correlation between FFT vectors

results = []
for j_sc in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]:
    H4 = star_hamiltonian([1.0, 2.0, j_sc], 4)
    L4 = [np.sqrt(0.05)*op_at(Z,i,4) for i in range(4)]
    psi4 = np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state())
    rho4 = np.outer(psi4, psi4.conj())
    freqs4, fft4 = get_spectrum_vector(H4, L4, rho4, [1,2], 4)
    
    # Spectral distance: use only overlapping frequency range
    # Both have same dt_sample so same frequency grid
    min_len = min(len(fft_base), len(fft4))
    fb = fft_base[:min_len]
    f4 = fft4[:min_len]
    
    # Normalize
    fb_n = fb / (np.linalg.norm(fb) + 1e-15)
    f4_n = f4 / (np.linalg.norm(f4) + 1e-15)
    
    # Cosine similarity (1 = identical, 0 = orthogonal)
    cos_sim = float(np.dot(fb_n, f4_n))
    spectral_dist = 1 - cos_sim
    
    # Count peaks in 4-qubit spectrum
    mask = freqs4 > 0.05
    fft_m = fft4[mask]
    freqs_m = freqs4[mask]
    mx = np.max(fft_m) if len(fft_m) > 0 else 0
    n_peaks = 0
    if mx > 0.01:
        for i in range(1, len(fft_m)-1):
            if (fft_m[i] > fft_m[i-1] and fft_m[i] > fft_m[i+1] 
                and fft_m[i] > 0.15 * mx):
                n_peaks += 1
    
    dom_f = freqs_m[np.argmax(fft_m)] if mx > 0.01 else 0
    
    detectable = "YES" if spectral_dist > 0.05 else "no"
    if spectral_dist > 0.3: detectable = "STRONG"
    
    results.append((j_sc, spectral_dist, n_peaks, dom_f, detectable))
    print(f"{j_sc:>5.2f} | {spectral_dist:>13.4f} | {n_peaks:>9} | {dom_f:>10.3f} | {detectable}")

# Find threshold
print(f"\n--- Detection threshold ---")
for i in range(len(results)-1):
    if results[i][4] == "no" and results[i+1][4] != "no":
        print(f"  Hidden observer becomes detectable between J_SC={results[i][0]} and J_SC={results[i+1][0]}")
        print(f"  Spectral distance jumps from {results[i][1]:.4f} to {results[i+1][1]:.4f}")
        break

# Can AB distinguish ONE hidden observer from TWO?
print(f"\n--- Can AB count hidden observers? ---")
print(f"{'Setup':>30} | {'Spec dist':>9} | {'Peaks':>5} | {'Dominant':>8}")
print("-" * 60)

configs = [
    ("AB alone", [1.0, 2.0], 3),
    ("AB + C(2.0)", [1.0, 2.0, 2.0], 4),
    ("AB + C(3.0)", [1.0, 2.0, 3.0], 4),
    ("AB + C(2.0) + D(2.0)", [1.0, 2.0, 2.0, 2.0], 5),
    ("AB + C(2.0) + D(3.0)", [1.0, 2.0, 2.0, 3.0], 5),
    ("AB + C(2.0)+D(3.0)+E(1.5)", [1.0, 2.0, 2.0, 3.0, 1.5], 6),
]

for name, J_list, nq in configs:
    H = star_hamiltonian(J_list, nq)
    L = [np.sqrt(0.05)*op_at(Z,i,nq) for i in range(nq)]
    # Build initial state: Bell_SA x |+> x |+> x ...
    psi = gpt.bell_phi_plus()
    for _ in range(nq - 2):
        psi = np.kron(psi, gpt.plus_state())
    rho = np.outer(psi, psi.conj())
    
    freqs_t, fft_t = get_spectrum_vector(H, L, rho, [1,2], nq)
    
    min_len = min(len(fft_base), len(fft_t))
    fb = fft_base[:min_len] / (np.linalg.norm(fft_base[:min_len]) + 1e-15)
    ft = fft_t[:min_len] / (np.linalg.norm(fft_t[:min_len]) + 1e-15)
    sd = 1 - float(np.dot(fb, ft))
    
    mask = freqs_t > 0.05
    fft_m = fft_t[mask]
    freqs_m = freqs_t[mask]
    mx = np.max(fft_m) if len(fft_m) > 0 else 0
    n_peaks = 0
    if mx > 0.01:
        for i in range(1, len(fft_m)-1):
            if (fft_m[i] > fft_m[i-1] and fft_m[i] > fft_m[i+1] 
                and fft_m[i] > 0.15 * mx):
                n_peaks += 1
    dom_f = freqs_m[np.argmax(fft_m)] if mx > 0.01 else 0
    
    print(f"{name:>30} | {sd:>9.4f} | {n_peaks:>5} | {dom_f:>8.3f}")

print(f"\n{'=' * 65}")
print("SONAR TEST COMPLETE")
print("=" * 65)
