"""
N_eff: Effective neighborhood size
GPT's formula: N_eff = (sum w_j)^2 / (sum w_j^2)
Test: does N_eff predict phase complexity better than raw degree?

Simulate star topologies with varying coupling patterns,
compute N_eff and phase std, check correlation.
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

def op_at(op, qubit, n_q):
    ops = [I2]*n_q
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def star_H(J_list, n_q):
    H = np.zeros((2**n_q, 2**n_q), dtype=complex)
    for i, J in enumerate(J_list):
        for p in [X, Y, Z]:
            H += J * op_at(p, 0, n_q) @ op_at(p, i+1, n_q)
    return H

def n_eff(weights):
    """Participation ratio: N_eff = (sum w)^2 / sum(w^2)"""
    w = np.array(weights)
    w = w[w > 0]  # only nonzero
    if len(w) == 0: return 0
    return (np.sum(w))**2 / np.sum(w**2)

def phase_complexity(H, L_ops, rho0, target_qubit, n_q, dt=0.005, t_max=10.0):
    """Measure phase std of target qubit's rho_01 over time."""
    phases = []
    rho = rho0.copy()
    for step in range(int(t_max/dt)+1):
        if step % 10 == 0 and step > 0:
            # Get single-qubit reduced state
            rho_q = gpt.partial_trace_keep(rho, [target_qubit], n_q)
            re = float(np.real(rho_q[0,1]))
            im = float(np.imag(rho_q[0,1]))
            if abs(re) + abs(im) > 0.01:
                phases.append(math.degrees(math.atan2(im, re)))
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    return np.std(phases) if len(phases) > 3 else 0

print("="*70)
print("N_eff: EFFECTIVE NEIGHBORHOOD SIZE")
print("Does N_eff predict phase complexity better than raw degree?")
print("="*70)

# Test configurations: vary degree AND coupling strengths
# Each config: (name, J_list for star around qubit 0, description)
configs = [
    # Degree 1 (one neighbor)
    ("1x strong",    [2.0],              "1 neighbor at J=2.0"),
    ("1x weak",      [0.3],              "1 neighbor at J=0.3"),
    
    # Degree 2 (two neighbors) - like Q102
    ("2x equal",     [1.0, 1.0],         "2 equal neighbors"),
    ("2x close",     [1.8, 2.0],         "2 CLOSE neighbors (Q102-like)"),
    ("2x spread",    [0.3, 2.0],         "2 spread neighbors"),
    ("2x very close",[1.9, 2.0],         "2 VERY close neighbors"),
    
    # Degree 3 (three neighbors) - like Q80
    ("3x equal",     [1.0, 1.0, 1.0],    "3 equal neighbors"),
    ("3x one dom",   [2.0, 0.1, 0.1],    "3 but 1 dominates (Q80-like)"),
    ("3x spread",    [0.5, 1.0, 2.0],    "3 spread neighbors"),
    ("3x all close", [1.8, 1.9, 2.0],    "3 ALL close neighbors"),
    
    # Degree 4
    ("4x equal",     [1.0, 1.0, 1.0, 1.0], "4 equal neighbors"),
    ("4x one dom",   [3.0, 0.1, 0.1, 0.1], "4 but 1 dominates"),
]

print(f"\n{'Name':>14} | {'Degree':>6} {'N_eff':>6} | {'Phase std':>9} | Description")
print("-"*70)

results = []
for name, J_list, desc in configs:
    n_q = 1 + len(J_list)
    degree = len(J_list)
    ne = n_eff(J_list)
    
    H = star_H(J_list, n_q)
    L = [np.sqrt(0.05)*op_at(Z,i,n_q) for i in range(n_q)]
    
    # Initial: |+> on all qubits
    psi = gpt.plus_state()
    for _ in range(n_q - 1):
        psi = np.kron(psi, gpt.plus_state())
    rho0 = np.outer(psi, psi.conj())
    
    pc = phase_complexity(H, L, rho0, 0, n_q)
    results.append((name, degree, ne, pc, desc))
    
    print(f"{name:>14} | {degree:>6} {ne:>6.2f} | {pc:>9.1f} | {desc}")

# Correlation analysis
degrees = [r[1] for r in results]
n_effs = [r[2] for r in results]
phase_stds = [r[3] for r in results]

# Pearson correlation
from numpy import corrcoef
r_degree = corrcoef(degrees, phase_stds)[0,1]
r_neff = corrcoef(n_effs, phase_stds)[0,1]

print(f"\n--- CORRELATION ANALYSIS ---")
print(f"  Degree vs Phase complexity:  r = {r_degree:.3f}")
print(f"  N_eff  vs Phase complexity:  r = {r_neff:.3f}")
print(f"  Winner: {'N_eff' if abs(r_neff) > abs(r_degree) else 'Degree'}")

# IBM comparison
print(f"\n--- IBM COMPARISON ---")
print(f"  Q80:  Degree=3, Phase std=12.4 deg (smooth)")
print(f"  Q102: Degree=2, Phase std=108.8 deg (chaotic)")
print(f"")
print(f"  If Q80's neighbors are [2.0, 0.1, 0.1] (one dominates):")
print(f"    N_eff = {n_eff([2.0, 0.1, 0.1]):.2f}")
print(f"  If Q102's neighbors are [1.8, 2.0] (both close):")
print(f"    N_eff = {n_eff([1.8, 2.0]):.2f}")
print(f"")
print(f"  Degree says: Q80 should be MORE complex (3 > 2). WRONG.")
print(f"  N_eff says:  Q80 N_eff={n_eff([2.0, 0.1, 0.1]):.2f} vs Q102 N_eff={n_eff([1.8, 2.0]):.2f}")
ne_q80 = n_eff([2.0, 0.1, 0.1])
ne_q102 = n_eff([1.8, 2.0])
if ne_q102 > ne_q80:
    print(f"  N_eff says:  Q102 should be MORE complex ({ne_q102:.2f} > {ne_q80:.2f}). CORRECT.")
else:
    print(f"  N_eff says:  Q80 more complex. Same as degree - still wrong.")

print(f"\n{'='*70}")
print("N_eff TEST COMPLETE")
print("="*70)
