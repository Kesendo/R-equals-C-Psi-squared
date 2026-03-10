"""
PARAMETER SWEEP: Does the two-sector structure survive?
Tests: varying J_SB, varying gamma, varying initial states.
For each: find dominant frequencies via FFT on c+ and c- coordinates.
"""
import numpy as np
import math
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

I2 = np.eye(2, dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)

def pauli_expect(rho, P):
    return float(np.real(np.trace(rho @ P)))

def run_and_extract_frequencies(J_SA, J_SB, gamma, initial_state_fn, dt=0.005, t_max=20.0):
    """Run simulation, extract c+ and c- time series, find dominant frequencies."""
    H = gpt.star_hamiltonian_n(n_observers=2, J_SA=J_SA, J_SB=J_SB)
    L_ops = gpt.dephasing_ops_n([gamma]*3)
    rho = initial_state_fn()

    YZ = np.kron(Y, Z)
    ZY = np.kron(Z, Y)
    XX_op = np.kron(X, X)
    
    times, cplus, cminus, cpsi_trace = [], [], [], []
    
    steps = int(t_max / dt)
    for step in range(steps + 1):
        t = step * dt
        if step % 4 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], 3)
            yz = pauli_expect(rAB, YZ)
            zy = pauli_expect(rAB, ZY)
            cp = (yz + zy) / math.sqrt(2)
            cm = (yz - zy) / math.sqrt(2)
            c = gpt.concurrence_two_qubit(rAB)
            p = gpt.psi_norm(rAB)
            
            # Check XX symmetry
            xx_comm = np.linalg.norm(rAB @ XX_op - XX_op @ rAB)
            
            times.append(t)
            cplus.append(cp)
            cminus.append(cm)
            cpsi_trace.append(c * p)
        if step < steps:
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    times = np.array(times)
    cplus = np.array(cplus)
    cminus = np.array(cminus)
    dt_s = times[1] - times[0]

    # FFT on c+ and c-
    freqs = np.fft.rfftfreq(len(times), d=dt_s)
    mask = freqs > 0.05
    
    fft_cp = np.abs(np.fft.rfft(cplus - np.mean(cplus)))
    fft_cm = np.abs(np.fft.rfft(cminus - np.mean(cminus)))
    
    f_cp = freqs[mask][np.argmax(fft_cp[mask])] if np.any(fft_cp[mask] > 0) else 0
    f_cm = freqs[mask][np.argmax(fft_cm[mask])] if np.any(fft_cm[mask] > 0) else 0
    
    amp_cp = np.max(fft_cp[mask]) if np.any(mask) else 0
    amp_cm = np.max(fft_cm[mask]) if np.any(mask) else 0
    
    # Check XX symmetry at last sampled point
    rAB_final = gpt.partial_trace_keep(rho, [1,2], 3)
    xx_comm_final = np.linalg.norm(rAB_final @ XX_op - XX_op @ rAB_final)
    
    # Ratio of frequencies
    ratio = f_cp / f_cm if f_cm > 0.01 else float('inf')
    
    return {
        'f_cplus': f_cp, 'f_cminus': f_cm,
        'amp_cplus': amp_cp, 'amp_cminus': amp_cm,
        'ratio': ratio,
        'xx_sym': xx_comm_final,
        'f_predicted': (J_SA + J_SB) / 2,
    }

# Initial state functions
def bell_plus_state():
    return gpt.density_from_statevector(
        np.kron(gpt.bell_phi_plus(), gpt.plus_state()))

def w_state():
    w = (np.array([0,1,1,0,1,0,0,0], dtype=complex) / math.sqrt(3))
    return np.outer(w, w.conj())

def product_state():
    return gpt.density_from_statevector(
        np.kron(np.kron(gpt.plus_state(), gpt.plus_state()), gpt.plus_state()))

# ============================================================
print("=" * 75)
print("PARAMETER SWEEP: Does the two-sector structure survive?")
print("=" * 75)

# SWEEP 1: Varying J_SB
print(f"\n--- SWEEP 1: Varying J_SB (J_SA=1.0, gamma=0.05, Bell+ initial) ---")
print(f"{'J_SB':>5} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'f_pred':>7} {'err%':>5} | {'XX_sym':>8} | note")
print("-" * 75)

for jsb in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
    r = run_and_extract_frequencies(1.0, jsb, 0.05, bell_plus_state)
    err = abs(r['f_cplus'] - r['f_predicted']) / r['f_predicted'] * 100 if r['f_predicted'] > 0 else 0
    sym = "EXACT" if r['xx_sym'] < 1e-10 else f"{r['xx_sym']:.1e}"
    note = ""
    if 3.0 < r['ratio'] < 5.0:
        note = "two sectors clear"
    elif r['ratio'] > 5.0:
        note = "sectors well separated"
    elif r['ratio'] < 2.0:
        note = "sectors merging?"
    print(f"{jsb:>5.1f} | {r['f_cplus']:>7.3f} {r['f_cminus']:>7.3f} {r['ratio']:>6.2f} | "
          f"{r['f_predicted']:>7.3f} {err:>5.1f} | {sym:>8} | {note}")

# SWEEP 2: Varying gamma
print(f"\n--- SWEEP 2: Varying gamma (J_SA=1.0, J_SB=2.0, Bell+ initial) ---")
print(f"{'gamma':>6} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'amp_c+':>7} {'amp_c-':>7} | note")
print("-" * 70)

for g in [0.001, 0.01, 0.05, 0.10, 0.20, 0.50, 1.00]:
    r = run_and_extract_frequencies(1.0, 2.0, g, bell_plus_state)
    note = ""
    if r['amp_cplus'] < 0.5:
        note = "c+ damped out"
    if r['amp_cminus'] < 0.1:
        note += " c- very weak"
    print(f"{g:>6.3f} | {r['f_cplus']:>7.3f} {r['f_cminus']:>7.3f} {r['ratio']:>6.2f} | "
          f"{r['amp_cplus']:>7.2f} {r['amp_cminus']:>7.2f} | {note}")

# SWEEP 3: Different initial states
print(f"\n--- SWEEP 3: Different initial states (J_SA=1.0, J_SB=2.0, gamma=0.05) ---")
print(f"{'State':>12} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'amp_c+':>7} {'amp_c-':>7} | {'XX_sym':>8}")
print("-" * 70)

for name, fn in [("Bell+", bell_plus_state), ("W-state", w_state), ("|+++>", product_state)]:
    r = run_and_extract_frequencies(1.0, 2.0, 0.05, fn)
    sym = "EXACT" if r['xx_sym'] < 1e-10 else f"{r['xx_sym']:.1e}"
    print(f"{name:>12} | {r['f_cplus']:>7.3f} {r['f_cminus']:>7.3f} {r['ratio']:>6.2f} | "
          f"{r['amp_cplus']:>7.2f} {r['amp_cminus']:>7.2f} | {sym:>8}")

# SWEEP 4: Break the symmetry - asymmetric gamma
print(f"\n--- SWEEP 4: Asymmetric gamma (J_SA=1.0, J_SB=2.0, Bell+ initial) ---")
print(f"{'gA':>5} {'gB':>5} {'gS':>5} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'XX_sym':>8} | note")
print("-" * 70)

for ga, gb, gs in [(0.05, 0.05, 0.05), (0.01, 0.10, 0.05), (0.10, 0.01, 0.05),
                    (0.05, 0.05, 0.20), (0.05, 0.05, 0.00)]:
    H = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
    L_ops = gpt.dephasing_ops_n([gs, ga, gb])
    psi = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
    rho_init = gpt.density_from_statevector(psi)
    
    YZ = np.kron(Y, Z)
    ZY = np.kron(Z, Y)
    XX_op = np.kron(X, X)
    
    times, cp_list, cm_list = [], [], []
    rho = rho_init.copy()
    dt = 0.005
    for step in range(4001):
        t = step * dt
        if step % 4 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], 3)
            yz = pauli_expect(rAB, YZ)
            zy = pauli_expect(rAB, ZY)
            times.append(t)
            cp_list.append((yz + zy) / math.sqrt(2))
            cm_list.append((yz - zy) / math.sqrt(2))
        if step < 4000:
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    times = np.array(times)
    dt_s = times[1] - times[0]
    freqs = np.fft.rfftfreq(len(times), d=dt_s)
    mask = freqs > 0.05
    
    fft_cp = np.abs(np.fft.rfft(np.array(cp_list) - np.mean(cp_list)))
    fft_cm = np.abs(np.fft.rfft(np.array(cm_list) - np.mean(cm_list)))
    f_cp = freqs[mask][np.argmax(fft_cp[mask])]
    f_cm = freqs[mask][np.argmax(fft_cm[mask])]
    
    rAB_f = gpt.partial_trace_keep(rho, [1,2], 3)
    xx_c = np.linalg.norm(rAB_f @ XX_op - XX_op @ rAB_f)
    sym = "EXACT" if xx_c < 1e-10 else f"{xx_c:.1e}"
    ratio = f_cp / f_cm if f_cm > 0.01 else float('inf')
    
    note = ""
    if xx_c > 1e-6:
        note = "XX BROKEN"
    
    print(f"{ga:>5.2f} {gb:>5.2f} {gs:>5.2f} | {f_cp:>7.3f} {f_cm:>7.3f} {ratio:>6.2f} | {sym:>8} | {note}")

print(f"\n{'=' * 75}")
print("PARAMETER SWEEP COMPLETE")
print("=" * 75)
