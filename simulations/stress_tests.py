"""
STRESS TESTS: Where does the two-sector structure break?
Tests that deliberately break symmetry to find the boundaries.

Test 1: Local fields (h_A * Z_A, h_B * Z_B) - breaks X*X symmetry
Test 2: Anisotropic coupling (XXZ instead of Heisenberg)
Test 3: Direct A-B coupling (J_AB != 0)
Test 4: Different jump operators (sigma_x instead of sigma_z)
"""
import numpy as np
import math
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

def pauli_expect(rho, P):
    return float(np.real(np.trace(rho @ P)))

def op_on_qubit(op, qubit, n_qubits):
    ops = [I2]*n_qubits
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def heisenberg_term(i, j, n_qubits, Jx=1.0, Jy=1.0, Jz=1.0):
    H = Jx * op_on_qubit(X, i, n_qubits) @ op_on_qubit(X, j, n_qubits)
    H += Jy * op_on_qubit(Y, i, n_qubits) @ op_on_qubit(Y, j, n_qubits)
    H += Jz * op_on_qubit(Z, i, n_qubits) @ op_on_qubit(Z, j, n_qubits)
    return H

def run_freq_analysis(H, L_ops, rho_init, dt=0.005, t_max=20.0):
    """Run simulation and extract c+/c- frequencies + XX symmetry."""
    YZ = np.kron(Y, Z)
    ZY = np.kron(Z, Y)
    XX_op = np.kron(X, X)
    
    times, cp_list, cm_list = [], [], []
    rho = rho_init.copy()
    
    for step in range(int(t_max/dt) + 1):
        t = step * dt
        if step % 4 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], 3)
            yz = pauli_expect(rAB, YZ)
            zy = pauli_expect(rAB, ZY)
            times.append(t)
            cp_list.append((yz + zy) / math.sqrt(2))
            cm_list.append((yz - zy) / math.sqrt(2))
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    times = np.array(times)
    dt_s = times[1] - times[0]
    freqs = np.fft.rfftfreq(len(times), d=dt_s)
    mask = freqs > 0.05
    
    fft_cp = np.abs(np.fft.rfft(np.array(cp_list) - np.mean(cp_list)))
    fft_cm = np.abs(np.fft.rfft(np.array(cm_list) - np.mean(cm_list)))
    
    f_cp = freqs[mask][np.argmax(fft_cp[mask])] if np.any(fft_cp[mask] > 0.01) else 0
    f_cm = freqs[mask][np.argmax(fft_cm[mask])] if np.any(fft_cm[mask] > 0.01) else 0
    amp_cp = np.max(fft_cp[mask]) if np.any(mask) else 0
    amp_cm = np.max(fft_cm[mask]) if np.any(mask) else 0
    
    # XX symmetry at midpoint
    rAB_mid = gpt.partial_trace_keep(rho, [1,2], 3)
    xx_comm = np.linalg.norm(rAB_mid @ XX_op - XX_op @ rAB_mid)
    
    ratio = f_cp / f_cm if f_cm > 0.01 else float('inf')
    return f_cp, f_cm, ratio, amp_cp, amp_cm, xx_comm

# Standard initial state
psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho0 = gpt.density_from_statevector(psi0)

# ============================================================
# BASELINE
# ============================================================
print("=" * 75)
print("STRESS TESTS: Where does the two-sector structure break?")
print("=" * 75)

H_base = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
L_base = gpt.dephasing_ops_n([0.05]*3)
f_cp, f_cm, ratio, amp_cp, amp_cm, xx = run_freq_analysis(H_base, L_base, rho0)
print(f"\nBASELINE (Heisenberg, sigma_z, no fields):")
print(f"  f(c+)={f_cp:.3f} f(c-)={f_cm:.3f} ratio={ratio:.2f} XX={xx:.1e}")

# ============================================================
# TEST 1: Local fields (break X*X symmetry)
# ============================================================
print(f"\n{'=' * 75}")
print("TEST 1: Local Z-fields on A and/or B")
print("Prediction: X*X symmetry should break, sectors may merge or shift")
print("=" * 75)
print(f"\n{'hA':>5} {'hB':>5} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'XX':>8} | note")
print("-" * 65)

for hA, hB in [(0.0, 0.0), (0.01, 0.0), (0.05, 0.0), (0.1, 0.0),
               (0.5, 0.0), (1.0, 0.0), (0.1, 0.1), (0.5, 0.5),
               (0.1, -0.1), (1.0, -1.0)]:
    H = H_base.copy()
    H += hA * op_on_qubit(Z, 1, 3)  # field on A
    H += hB * op_on_qubit(Z, 2, 3)  # field on B
    f_cp, f_cm, ratio, amp_cp, amp_cm, xx = run_freq_analysis(H, L_base, rho0)
    sym = "EXACT" if xx < 1e-10 else f"{xx:.1e}"
    note = ""
    if xx > 1e-6 and xx < 0.01: note = "XX slightly broken"
    elif xx > 0.01: note = "XX BROKEN"
    if abs(ratio - 3.75) < 0.3 and xx < 1e-10: note = "structure intact"
    print(f"{hA:>5.2f} {hB:>5.2f} | {f_cp:>7.3f} {f_cm:>7.3f} {ratio:>6.2f} | {sym:>8} | {note}")

# ============================================================
# TEST 2: Anisotropic coupling (XXZ)
# ============================================================
print(f"\n{'=' * 75}")
print("TEST 2: Anisotropic coupling (Jz != Jx=Jy)")
print("Prediction: XX symmetry may survive (depends on anisotropy axis)")
print("=" * 75)
print(f"\n{'Jz/Jxy':>7} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'XX':>8} | note")
print("-" * 60)

for jz_ratio in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]:
    # SA coupling: Jx=Jy=1.0, Jz=jz_ratio
    H = heisenberg_term(0, 1, 3, Jx=1.0, Jy=1.0, Jz=jz_ratio)
    # SB coupling: Jx=Jy=2.0, Jz=2*jz_ratio (keep ratio same)
    H += heisenberg_term(0, 2, 3, Jx=2.0, Jy=2.0, Jz=2.0*jz_ratio)
    
    f_cp, f_cm, ratio, amp_cp, amp_cm, xx = run_freq_analysis(H, L_base, rho0)
    sym = "EXACT" if xx < 1e-10 else f"{xx:.1e}"
    note = ""
    if xx < 1e-10 and ratio > 2: note = "structure survives"
    elif xx > 1e-6: note = "XX broken"
    print(f"{jz_ratio:>7.2f} | {f_cp:>7.3f} {f_cm:>7.3f} {ratio:>6.2f} | {sym:>8} | {note}")

# ============================================================
# TEST 3: Direct A-B coupling
# ============================================================
print(f"\n{'=' * 75}")
print("TEST 3: Direct A-B coupling (J_AB != 0)")
print("Prediction: changes spectrum but may preserve symmetry")
print("=" * 75)
print(f"\n{'J_AB':>5} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'XX':>8} | note")
print("-" * 55)

for jab in [0.0, 0.1, 0.3, 0.5, 1.0, 2.0]:
    H = H_base.copy()
    H += heisenberg_term(1, 2, 3, Jx=jab, Jy=jab, Jz=jab)
    
    f_cp, f_cm, ratio, amp_cp, amp_cm, xx = run_freq_analysis(H, L_base, rho0)
    sym = "EXACT" if xx < 1e-10 else f"{xx:.1e}"
    note = ""
    if xx < 1e-10 and ratio > 2: note = "structure survives"
    print(f"{jab:>5.2f} | {f_cp:>7.3f} {f_cm:>7.3f} {ratio:>6.2f} | {sym:>8} | {note}")

# ============================================================
# TEST 4: Different jump operators
# ============================================================
print(f"\n{'=' * 75}")
print("TEST 4: Different jump operators (not just sigma_z)")
print("Prediction: sigma_x/sigma_y may break X*X symmetry")
print("=" * 75)
print(f"\n{'Jump':>10} | {'f(c+)':>7} {'f(c-)':>7} {'ratio':>6} | {'XX':>8} | note")
print("-" * 60)

for jump_name, jump_op in [("sigma_z", Z), ("sigma_x", X), ("sigma_y", Y),
                             ("mixed_xz", None)]:
    if jump_name == "mixed_xz":
        # sigma_x on A, sigma_z on B, sigma_z on S
        L_ops = [
            np.sqrt(0.05) * op_on_qubit(Z, 0, 3),
            np.sqrt(0.05) * op_on_qubit(X, 1, 3),
            np.sqrt(0.05) * op_on_qubit(Z, 2, 3),
        ]
    else:
        L_ops = [np.sqrt(0.05) * op_on_qubit(jump_op, i, 3) for i in range(3)]
    
    f_cp, f_cm, ratio, amp_cp, amp_cm, xx = run_freq_analysis(H_base, L_ops, rho0)
    sym = "EXACT" if xx < 1e-10 else f"{xx:.1e}"
    note = ""
    if xx < 1e-10 and ratio > 2: note = "structure survives"
    elif xx > 1e-6: note = "XX BROKEN"
    print(f"{jump_name:>10} | {f_cp:>7.3f} {f_cm:>7.3f} {ratio:>6.2f} | {sym:>8} | {note}")

print(f"\n{'=' * 75}")
print("STRESS TESTS COMPLETE")
print("=" * 75)
