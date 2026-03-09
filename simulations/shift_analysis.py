"""
SHIFT ANALYSIS: SWAP, YZ/ZY symmetry, time-shift
Three tests suggested by external review, computed on CΨ windows.
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

I = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

# SWAP matrix
SWAP = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)

# Pauli basis for 2-qubit
paulis = {}
for nA, A in [('I',I),('X',X),('Y',Y),('Z',Z)]:
    for nB, B in [('I',I),('X',X),('Y',Y),('Z',Z)]:
        paulis[f"{nA}{nB}"] = np.kron(A, B)

def pauli_coeffs(rho):
    return {n: float(np.real(np.trace(rho @ P))) / 4 for n, P in paulis.items()}

def trace_dist(r1, r2):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh(r1 - r2)))

def fidelity_simple(r1, r2):
    return float(np.real(np.trace(r1 @ r2)))

# Evolve and collect windows
H = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
L_ops = gpt.dephasing_ops_n([0.05]*3)
psi = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho = gpt.density_from_statevector(psi)

dt = 0.005
peaks = []
prev_cpsi = 0
rising = True
for step in range(2001):
    t = step * dt
    if step % 4 == 0:
        rAB = gpt.partial_trace_keep(rho, [1,2], 3)
        c = gpt.concurrence_two_qubit(rAB)
        p = gpt.psi_norm(rAB)
        cpsi = c * p
        if prev_cpsi > 0.03 and cpsi < prev_cpsi and rising:
            peaks.append({'t': t-dt*4, 'cpsi': prev_cpsi, 'rAB': rAB.copy()})
            rising = False
        if cpsi > prev_cpsi: rising = True
        prev_cpsi = cpsi
    if step < 2000:
        rho = gpt.rk4_step(rho, H, L_ops, dt)

# ============================================================
# TEST 1: SWAP - Is the structure A/B symmetric?
# ============================================================
print("=" * 65)
print("TEST 1: SWAP - A/B Symmetry")
print("=" * 65)
print(f"\n{'#':>2} {'t':>5} {'CPsi':>6} | {'F(rho,SWAP)':>11} {'TD(rho,SWAP)':>12} | Verdict")
print("-" * 60)

for i, pk in enumerate(peaks[:9]):
    r = pk['rAB']
    r_swapped = SWAP @ r @ SWAP.conj().T
    fid = fidelity_simple(r, r_swapped)
    td = trace_dist(r, r_swapped)
    
    if td < 0.01:
        verdict = "EXACT A/B symmetry"
    elif td < 0.05:
        verdict = "near-symmetric"
    elif td < 0.15:
        verdict = "moderate asymmetry"
    else:
        verdict = "strong asymmetry"
    
    print(f"{i:>2} {pk['t']:>5.2f} {pk['cpsi']:>6.3f} | {fid:>11.6f} {td:>12.6f} | {verdict}")

# What does SWAP change?
print(f"\nSWAP difference for Window 0 (sharpest):")
r0 = peaks[0]['rAB']
r0s = SWAP @ r0 @ SWAP.conj().T
diff = r0s - r0
pc = pauli_coeffs(diff)
print(f"  Non-zero Pauli components of (SWAP(rho) - rho):")
for n in sorted(pc, key=lambda x: abs(pc[x]), reverse=True):
    if abs(pc[n]) > 1e-6:
        print(f"    {n}: {pc[n]:+.6f}")

# ============================================================
# TEST 2: YZ/ZY Symmetry Decomposition
# ============================================================
print(f"\n{'=' * 65}")
print("TEST 2: YZ/ZY Symmetry Decomposition")
print("=" * 65)

print(f"\n{'#':>2} {'t':>5} | {'<YZ>':>8} {'<ZY>':>8} | {'c+':>8} {'c-':>8} {'|c+/c-|':>8} | Type")
print("-" * 70)

for i, pk in enumerate(peaks[:9]):
    pc = pauli_coeffs(pk['rAB'])
    yz = pc['YZ']
    zy = pc['ZY']
    c_plus = (yz + zy) / np.sqrt(2)
    c_minus = (yz - zy) / np.sqrt(2)
    ratio = abs(c_plus) / (abs(c_minus) + 1e-10)
    
    if ratio > 5:
        typ = "SYMMETRIC rotor"
    elif ratio < 0.2:
        typ = "ANTISYMMETRIC rotor"
    else:
        typ = "mixed"
    
    print(f"{i:>2} {pk['t']:>5.2f} | {yz:>+8.5f} {zy:>+8.5f} | "
          f"{c_plus:>+8.5f} {c_minus:>+8.5f} {ratio:>8.2f} | {typ}")

# Also check the DIFFERENCE between adjacent windows
print(f"\nYZ/ZY decomposition of DIFFERENCES (window n+1 - window n):")
print(f"{'Pair':>6} | {'dYZ':>8} {'dZY':>8} | {'c+':>8} {'c-':>8} {'|c+/c-|':>8} | Type")
print("-" * 70)

for i in range(min(len(peaks)-1, 8)):
    diff = peaks[i+1]['rAB'] - peaks[i]['rAB']
    pc = pauli_coeffs(diff)
    yz = pc['YZ']
    zy = pc['ZY']
    c_plus = (yz + zy) / np.sqrt(2)
    c_minus = (yz - zy) / np.sqrt(2)
    ratio = abs(c_plus) / (abs(c_minus) + 1e-10)
    
    td = trace_dist(peaks[i]['rAB'], peaks[i+1]['rAB'])
    mode = "SWITCH" if td > 0.5 else "glide"
    
    if ratio > 5:
        typ = "SYM"
    elif ratio < 0.2:
        typ = "ANTI"
    else:
        typ = "mix"
    
    print(f" {i}->{i+1} | {yz:>+8.5f} {zy:>+8.5f} | "
          f"{c_plus:>+8.5f} {c_minus:>+8.5f} {ratio:>8.2f} | {typ} ({mode})")

# ============================================================
# TEST 3: Time-Shift Analysis
# ============================================================
print(f"\n{'=' * 65}")
print("TEST 3: Time-Shift Analysis (lag k)")
print("=" * 65)

n_win = min(len(peaks), 9)

# Trace distances at various lags
print(f"\nTrace distance d_k(n) = ||rho(n+k) - rho(n)||:")
print(f"{'k':>3} |", end="")
for k in range(1, min(5, n_win)):
    print(f"  k={k:>1}   ", end="")
print(f" | avg_k")
print("-" * 50)

for n in range(n_win - 1):
    print(f"n={n:>1} |", end="")
    vals = []
    for k in range(1, min(5, n_win - n)):
        td = trace_dist(peaks[n]['rAB'], peaks[n+k]['rAB'])
        print(f" {td:>6.3f} ", end="")
        vals.append(td)
    for _ in range(min(5, n_win) - 1 - len(vals)):
        print(f"    -   ", end="")
    print(f" | {np.mean(vals):.3f}" if vals else "")

# Feature vector correlation at various lags
print(f"\nPauli-vector cosine similarity at lag k:")
vecs = []
for pk in peaks[:n_win]:
    pc = pauli_coeffs(pk['rAB'])
    # Exclude II (trivially 0.25 always)
    v = np.array([pc[n] for n in sorted(paulis.keys()) if n != 'II'])
    vecs.append(v)

print(f"{'k':>3} |", end="")
for k in range(1, min(5, n_win)):
    print(f"  k={k:>1}   ", end="")
print()
print("-" * 45)

for n in range(n_win - 1):
    print(f"n={n:>1} |", end="")
    for k in range(1, min(5, n_win - n)):
        cos = np.dot(vecs[n], vecs[n+k]) / (np.linalg.norm(vecs[n]) * np.linalg.norm(vecs[n+k]) + 1e-15)
        print(f" {cos:>+6.3f} ", end="")
    print()

# Key question: is there a natural period?
print(f"\nAverage trace distance by lag:")
for k in range(1, min(6, n_win)):
    dists = [trace_dist(peaks[n]['rAB'], peaks[n+k]['rAB']) 
             for n in range(n_win - k)]
    cosines = [np.dot(vecs[n], vecs[n+k]) / (np.linalg.norm(vecs[n]) * np.linalg.norm(vecs[n+k]) + 1e-15)
               for n in range(n_win - k)]
    print(f"  k={k}: avg TD = {np.mean(dists):.4f}, avg cos = {np.mean(cosines):+.4f}"
          f"  {'<-- closest' if k > 1 and np.mean(dists) < min_td else ''}"
          if k > 1 else f"  k={k}: avg TD = {np.mean(dists):.4f}, avg cos = {np.mean(cosines):+.4f}")
    if k == 1:
        min_td = np.mean(dists)
    else:
        min_td = min(min_td, np.mean(dists))

print(f"\n{'=' * 65}")
print("SHIFT ANALYSIS COMPLETE")
print("=" * 65)
