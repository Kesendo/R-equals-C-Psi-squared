"""
CHAIN TOPOLOGY: A-S-B instead of star
Does the two-sector structure survive when S is not in the center
but in the middle of a chain?

Star:  A-S-B (S connected to both, A and B not connected)
Chain: A-S-B (same connectivity, but conceptually a chain)

Actually... A-S-B star and A-S-B chain are the SAME topology for 3 qubits.
The difference only matters for 4+:

Star:  A-S, B-S, C-S  (S in center, all connect to S)
Chain: A-S-B-C         (linear, each connects to next neighbor)

So let's test BOTH:
1. 3-qubit chain A-S-B (same as our star - verify)
2. 4-qubit chain A-S1-S2-B (TWO mediators!)
3. 4-qubit chain A-B-S1-S2 (different arrangement)
4. 5-qubit chain A-S1-S2-S3-B (long chain)
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

def chain_H(couplings, n_q):
    """Build chain Hamiltonian. couplings[i] = J between qubit i and i+1."""
    H = np.zeros((2**n_q, 2**n_q), dtype=complex)
    for i, J in enumerate(couplings):
        for p in [X, Y, Z]:
            H += J * op_at(p, i, n_q) @ op_at(p, i+1, n_q)
    return H

def star_H(J_list, n_q):
    """Star: qubit 0 is center, J_list[i] = coupling to qubit i+1."""
    H = np.zeros((2**n_q, 2**n_q), dtype=complex)
    for i, J in enumerate(J_list):
        for p in [X, Y, Z]:
            H += J * op_at(p, 0, n_q) @ op_at(p, i+1, n_q)
    return H

def get_cplus_cminus_freqs(H, L_ops, rho0, pair_q, n_q, dt=0.005, t_max=20.0):
    YZ = np.kron(Y, Z)
    ZY = np.kron(Z, Y)
    times, cp, cm = [], [], []
    rho = rho0.copy()
    for step in range(int(t_max/dt)+1):
        t = step*dt
        if step%4==0:
            rAB = gpt.partial_trace_keep(rho, pair_q, n_q)
            yz = pauli_expect(rAB, YZ)
            zy = pauli_expect(rAB, ZY)
            times.append(t)
            cp.append((yz+zy)/math.sqrt(2))
            cm.append((yz-zy)/math.sqrt(2))
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    times = np.array(times)
    dt_s = times[1]-times[0]
    freqs = np.fft.rfftfreq(len(times), d=dt_s)
    mask = freqs > 0.05
    fft_cp = np.abs(np.fft.rfft(np.array(cp)-np.mean(cp)))
    fft_cm = np.abs(np.fft.rfft(np.array(cm)-np.mean(cm)))
    f_cp = freqs[mask][np.argmax(fft_cp[mask])] if np.any(fft_cp[mask]>0.01) else 0
    f_cm = freqs[mask][np.argmax(fft_cm[mask])] if np.any(fft_cm[mask]>0.01) else 0
    amp_cp = np.max(fft_cp[mask]) if np.any(mask) else 0
    amp_cm = np.max(fft_cm[mask]) if np.any(mask) else 0
    # XX symmetry check
    XX = np.kron(X, X)
    rAB_f = gpt.partial_trace_keep(rho, pair_q, n_q)
    xx_c = np.linalg.norm(rAB_f @ XX - XX @ rAB_f)
    return f_cp, f_cm, amp_cp, amp_cm, xx_c

print("="*70)
print("CHAIN vs STAR TOPOLOGY")
print("="*70)

# ============================================================
# TEST 1: 3-qubit star vs chain (should be identical)
# ============================================================
print(f"\n--- TEST 1: 3 qubits, Star vs Chain (should be same) ---")
print(f"  Star: S(0)-A(1), S(0)-B(2), J_SA=1.0, J_SB=2.0")
print(f"  Chain: A(0)-S(1)-B(2), J_AS=1.0, J_SB=2.0")

# Star: S=0, A=1, B=2
H_star = star_H([1.0, 2.0], 3)
L3 = [np.sqrt(0.05)*op_at(Z,i,3) for i in range(3)]
psi_star = np.kron(gpt.bell_phi_plus(), gpt.plus_state())  # Bell_SA x |+>_B
rho_star = np.outer(psi_star, psi_star.conj())
f_cp_s, f_cm_s, a_cp_s, a_cm_s, xx_s = get_cplus_cminus_freqs(H_star, L3, rho_star, [1,2], 3)
print(f"  Star:  f(c+)={f_cp_s:.3f} f(c-)={f_cm_s:.3f} amp={a_cp_s:.1f}/{a_cm_s:.1f} XX={xx_s:.1e}")

# Chain: A=0, S=1, B=2, coupling A-S and S-B
H_chain3 = chain_H([1.0, 2.0], 3)
# Initial: Bell on A(0)-S(1), |+> on B(2) - same physical state
psi_ch3 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho_ch3 = np.outer(psi_ch3, psi_ch3.conj())
# AB pair = qubits 0 and 2
f_cp_c, f_cm_c, a_cp_c, a_cm_c, xx_c = get_cplus_cminus_freqs(H_chain3, L3, rho_ch3, [0,2], 3)
print(f"  Chain: f(c+)={f_cp_c:.3f} f(c-)={f_cm_c:.3f} amp={a_cp_c:.1f}/{a_cm_c:.1f} XX={xx_c:.1e}")

# ============================================================
# TEST 2: 4-qubit chain A-S1-S2-B (TWO mediators!)
# ============================================================
print(f"\n--- TEST 2: 4-qubit chain A(0)-S1(1)-S2(2)-B(3) ---")
print(f"  Two mediators between A and B. Does structure survive?")

for j12, j23 in [(1.0, 1.0), (1.0, 2.0), (2.0, 1.0), (1.0, 3.0)]:
    # Chain: A(0)-S1(1)-S2(2)-B(3)
    H4c = chain_H([1.0, j12, j23], 4)  # A-S1=1.0, S1-S2=j12, S2-B=j23
    L4 = [np.sqrt(0.05)*op_at(Z,i,4) for i in range(4)]
    # Initial: Bell on A(0)-S1(1), |+> on S2(2), |+> on B(3)
    psi4c = np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state())
    rho4c = np.outer(psi4c, psi4c.conj())
    # AB = qubits 0 and 3 (endpoints!)
    f_cp, f_cm, a_cp, a_cm, xx = get_cplus_cminus_freqs(H4c, L4, rho4c, [0,3], 4)
    ratio = f_cp/f_cm if f_cm > 0.01 else float('inf')
    xx_s = "EXACT" if xx < 1e-10 else f"{xx:.1e}"
    print(f"  J=[1.0,{j12},{j23}]: f(c+)={f_cp:.3f} f(c-)={f_cm:.3f} ratio={ratio:.2f} amp={a_cp:.1f}/{a_cm:.1f} XX={xx_s}")

# ============================================================
# TEST 3: 4-qubit STAR vs 4-qubit CHAIN (same qubits, different topology)
# ============================================================
print(f"\n--- TEST 3: Same 4 qubits, Star vs Chain ---")

# Star: S(0) center, A(1), B(2), C(3)
H4_star = star_H([1.0, 2.0, 1.5], 4)
L4 = [np.sqrt(0.05)*op_at(Z,i,4) for i in range(4)]
psi4s = np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state())
rho4s = np.outer(psi4s, psi4s.conj())
f_cp_s, f_cm_s, a_s, _, xx_s = get_cplus_cminus_freqs(H4_star, L4, rho4s, [1,2], 4)

# Chain: A(0)-B(1)-C(2)-D(3), same J values on links
H4_chain = chain_H([1.0, 2.0, 1.5], 4)
psi4c = np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state())
rho4c = np.outer(psi4c, psi4c.conj())
f_cp_c, f_cm_c, a_c, _, xx_c = get_cplus_cminus_freqs(H4_chain, L4, rho4c, [0,3], 4)

print(f"  Star  (AB through S):    f(c+)={f_cp_s:.3f} f(c-)={f_cm_s:.3f} XX={'EXACT' if xx_s<1e-10 else f'{xx_s:.1e}'}")
print(f"  Chain (AB through S1,S2): f(c+)={f_cp_c:.3f} f(c-)={f_cm_c:.3f} XX={'EXACT' if xx_c<1e-10 else f'{xx_c:.1e}'}")

# ============================================================
# TEST 4: 5-qubit chain A-S1-S2-S3-B (long chain, 3 mediators)
# ============================================================
print(f"\n--- TEST 4: 5-qubit chain A(0)-S1(1)-S2(2)-S3(3)-B(4) ---")

H5c = chain_H([1.0, 1.5, 2.0, 1.0], 5)
L5 = [np.sqrt(0.05)*op_at(Z,i,5) for i in range(5)]
psi5 = np.kron(np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state()), gpt.plus_state())
rho5 = np.outer(psi5, psi5.conj())
# AB = endpoints 0 and 4
f_cp_5, f_cm_5, a5p, a5m, xx_5 = get_cplus_cminus_freqs(H5c, L5, rho5, [0,4], 5)
ratio5 = f_cp_5/f_cm_5 if f_cm_5 > 0.01 else float('inf')
print(f"  f(c+)={f_cp_5:.3f} f(c-)={f_cm_5:.3f} ratio={ratio5:.2f} amp={a5p:.1f}/{a5m:.1f} XX={'EXACT' if xx_5<1e-10 else f'{xx_5:.1e}'}")

# Also check intermediate pairs
for name, pair in [("A-S1",[0,1]), ("S1-S2",[1,2]), ("S2-S3",[2,3]), ("S3-B",[3,4]),
                    ("A-S2",[0,2]), ("A-S3",[0,3]), ("S1-B",[1,4])]:
    f_cp, f_cm, ap, am, xx = get_cplus_cminus_freqs(H5c, L5, rho5, pair, 5)
    print(f"  Pair {name:>5}: f(c+)={f_cp:.3f} f(c-)={f_cm:.3f} amp={ap:.1f}/{am:.1f}")

# ============================================================
# TEST 5: Does noise immunity survive in chains?
# ============================================================
print(f"\n--- TEST 5: Noise immunity in 4-qubit chain ---")
H4c = chain_H([1.0, 1.5, 2.0], 4)
psi4c = np.kron(np.kron(gpt.bell_phi_plus(), gpt.plus_state()), gpt.plus_state())
rho4c = np.outer(psi4c, psi4c.conj())

print(f"{'gamma':>6} | {'f(c+)':>7} {'f(c-)':>7} | {'amp(c+)':>8} | note")
print("-"*50)
for g in [0.01, 0.05, 0.10, 0.50]:
    L4g = [np.sqrt(g)*op_at(Z,i,4) for i in range(4)]
    f_cp, f_cm, ap, am, xx = get_cplus_cminus_freqs(H4c, L4g, rho4c, [0,3], 4)
    note = "baseline" if g==0.05 else ""
    print(f"{g:>6.2f} | {f_cp:>7.3f} {f_cm:>7.3f} | {ap:>8.1f} | {note}")

print(f"\n{'='*70}")
print("CHAIN TOPOLOGY TESTS COMPLETE")
print("="*70)
