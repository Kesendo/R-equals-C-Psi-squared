"""
Which operation transforms Window 0 into Window 1?
Test all single-qubit and two-qubit Pauli operations.
Find which one maximizes the overlap after transformation.
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

# Pauli matrices
I = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

# All single and two-qubit Pauli operations on AB
ops = {}
for nA, A in [('I',I),('X',X),('Y',Y),('Z',Z)]:
    for nB, B in [('I',I),('X',X),('Y',Y),('Z',Z)]:
        name = f"{nA}x{nB}"
        ops[name] = np.kron(A, B)

def apply_unitary(rho, U):
    return U @ rho @ U.conj().T

def overlap_fraction(r1, r2):
    m1 = np.abs(r1)
    m2 = np.abs(r2)
    shared = np.sum(np.minimum(m1, m2))
    total = max(np.sum(m1), np.sum(m2))
    return shared / total if total > 0 else 0

def trace_dist(r1, r2):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh(r1 - r2)))

def fidelity_dm(r1, r2):
    return float(np.real(np.trace(r1 @ r2))) + 2*np.real(
        np.sqrt(max(0, np.linalg.det(r1) * np.linalg.det(r2))))

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

print("=" * 70)
print("OPERATION SEARCH: Which transform maps Window_i to Window_j?")
print(f"Found {len(peaks)} windows")
print("=" * 70)

# Test 1: For the two sharpest windows, which operation maximizes overlap?
print(f"\n--- Test 1: Transform Window 0 -> Window 1 ---")
print(f"Window 0: t={peaks[0]['t']:.2f}, CPsi={peaks[0]['cpsi']:.3f}")
print(f"Window 1: t={peaks[1]['t']:.2f}, CPsi={peaks[1]['cpsi']:.3f}")
print(f"\nNo transform (raw):  overlap={overlap_fraction(peaks[0]['rAB'], peaks[1]['rAB'])*100:.1f}%"
      f"  TD={trace_dist(peaks[0]['rAB'], peaks[1]['rAB']):.4f}")

print(f"\n{'Operation':>6} | {'overlap%':>8} {'TD':>8} {'dTD':>8} | note")
print("-" * 55)

results = []
r0 = peaks[0]['rAB']
r1 = peaks[1]['rAB']
raw_td = trace_dist(r0, r1)

for name, U in sorted(ops.items()):
    r0_transformed = apply_unitary(r0, U)
    ov = overlap_fraction(r0_transformed, r1)
    td = trace_dist(r0_transformed, r1)
    dtd = td - raw_td
    results.append((name, ov, td, dtd))

results.sort(key=lambda x: x[2])  # sort by trace distance (lower = better match)

for name, ov, td, dtd in results:
    note = ""
    if td < raw_td * 0.5:
        note = "<-- IMPROVES"
    elif td < raw_td * 0.8:
        note = "<-- helps"
    elif td > raw_td * 1.2:
        note = "worsens"
    print(f"{name:>6} | {ov*100:>7.1f}% {td:>8.4f} {dtd:>+8.4f} | {note}")

# Test 2: Does the SAME operation work for ALL adjacent pairs?
print(f"\n{'=' * 70}")
print("Test 2: Best operation for EACH adjacent window pair")
print("=" * 70)
print(f"\n{'Pair':>6} | {'Raw TD':>7} | {'Best Op':>6} {'Best TD':>8} {'Improv%':>8} | {'2nd Op':>6} {'2nd TD':>8}")
print("-" * 75)

best_ops_per_pair = []
for k in range(min(len(peaks)-1, 8)):
    ra = peaks[k]['rAB']
    rb = peaks[k+1]['rAB']
    raw = trace_dist(ra, rb)
    
    pair_results = []
    for name, U in ops.items():
        ra_t = apply_unitary(ra, U)
        td = trace_dist(ra_t, rb)
        pair_results.append((name, td))
    
    pair_results.sort(key=lambda x: x[1])
    best = pair_results[0]
    second = pair_results[1]
    improv = (raw - best[1]) / raw * 100 if raw > 0 else 0
    
    best_ops_per_pair.append(best[0])
    
    print(f" {k}->{k+1} | {raw:>7.4f} | {best[0]:>6} {best[1]:>8.4f} {improv:>+7.1f}% | "
          f"{second[0]:>6} {second[1]:>8.4f}")

# Test 3: Is there one UNIVERSAL operation?
print(f"\n{'=' * 70}")
print("Test 3: Is there ONE operation that works for all pairs?")
print("=" * 70)

# For each operation, sum the trace distances across all pairs
total_scores = {}
for name, U in ops.items():
    total = 0
    for k in range(min(len(peaks)-1, 8)):
        ra_t = apply_unitary(peaks[k]['rAB'], U)
        total += trace_dist(ra_t, peaks[k+1]['rAB'])
    total_scores[name] = total

raw_total = sum(trace_dist(peaks[k]['rAB'], peaks[k+1]['rAB']) 
                for k in range(min(len(peaks)-1, 8)))

print(f"\nRaw total TD (no operation): {raw_total:.4f}")
print(f"\n{'Op':>6} | {'Total TD':>9} | {'vs raw':>8} | note")
print("-" * 45)
for name in sorted(total_scores, key=total_scores.get):
    t = total_scores[name]
    vs = (t - raw_total) / raw_total * 100
    note = "BEST" if t == min(total_scores.values()) else ""
    if t < raw_total * 0.9:
        note = "universal help"
    print(f"{name:>6} | {t:>9.4f} | {vs:>+7.1f}% | {note}")

# Test 4: What about continuous rotations? Find the optimal angle
print(f"\n{'=' * 70}")
print("Test 4: Continuous rotations - optimal angle for W0->W1")
print("=" * 70)

r0 = peaks[0]['rAB']
r1 = peaks[1]['rAB']

# Test Rz(theta) on A, Rz(theta) on B, and both
for axis_name, axis_gen in [("Rz_A", lambda t: np.kron(np.array([[np.exp(-1j*t/2),0],[0,np.exp(1j*t/2)]]), I)),
                             ("Rz_B", lambda t: np.kron(I, np.array([[np.exp(-1j*t/2),0],[0,np.exp(1j*t/2)]]))),
                             ("Ry_A", lambda t: np.kron(np.array([[np.cos(t/2),-np.sin(t/2)],[np.sin(t/2),np.cos(t/2)]]), I)),
                             ("Ry_B", lambda t: np.kron(I, np.array([[np.cos(t/2),-np.sin(t/2)],[np.sin(t/2),np.cos(t/2)]]))),
                             ("Rz_AB", lambda t: np.kron(np.array([[np.exp(-1j*t/2),0],[0,np.exp(1j*t/2)]]),
                                                          np.array([[np.exp(-1j*t/2),0],[0,np.exp(1j*t/2)]])))]:
    
    best_td = 999
    best_angle = 0
    for angle_deg in range(-180, 181, 5):
        theta = angle_deg * np.pi / 180
        U = axis_gen(theta)
        r0_t = apply_unitary(r0, U)
        td = trace_dist(r0_t, r1)
        if td < best_td:
            best_td = td
            best_angle = angle_deg
    
    # Fine-tune around best
    for angle_deg_fine in np.linspace(best_angle - 5, best_angle + 5, 101):
        theta = angle_deg_fine * np.pi / 180
        U = axis_gen(theta)
        r0_t = apply_unitary(r0, U)
        td = trace_dist(r0_t, r1)
        if td < best_td:
            best_td = td
            best_angle = angle_deg_fine
    
    raw = trace_dist(r0, r1)
    improv = (raw - best_td) / raw * 100
    print(f"  {axis_name:>6}: best angle = {best_angle:>+7.1f} deg, TD = {best_td:.4f} "
          f"(was {raw:.4f}, {improv:>+.1f}% better)")

# Test 5: Repeat for ALL adjacent pairs - same optimal axis?
print(f"\n{'=' * 70}")
print("Test 5: Optimal Rz_A angle for each adjacent pair")
print("=" * 70)

for k in range(min(len(peaks)-1, 8)):
    ra = peaks[k]['rAB']
    rb = peaks[k+1]['rAB']
    
    best_td = 999
    best_angle = 0
    for angle_deg in range(-180, 181, 2):
        theta = angle_deg * np.pi / 180
        U = np.kron(np.array([[np.exp(-1j*theta/2),0],[0,np.exp(1j*theta/2)]]), I)
        td = trace_dist(apply_unitary(ra, U), rb)
        if td < best_td:
            best_td = td
            best_angle = angle_deg
    
    raw = trace_dist(ra, rb)
    improv = (raw - best_td) / raw * 100
    mode = "SWITCH" if raw > 0.5 else "glide"
    print(f"  {k}->{k+1} ({mode:>6}): Rz_A = {best_angle:>+4d} deg, "
          f"TD {raw:.3f} -> {best_td:.3f} ({improv:>+.0f}%)")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
