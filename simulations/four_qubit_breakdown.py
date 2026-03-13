"""
4-QUBIT BREAKDOWN: How exactly does frequency-decay orthogonality break?
========================================================================
In the 3-qubit star, decay rates {2g, 8g/3, 10g/3} are EXACT and
topology-independent. At 4 qubits, some rates start moving with J.

Questions:
1. Which rates stay fixed? Which move?
2. How do they move as a function of J?
3. Is there a pattern to the breakdown?
4. Does the topology (star vs chain) matter?
"""
import numpy as np
import math
import sys

sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")

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

def build_liouvillian(n_q, bonds, gamma):
    """Build Liouvillian for arbitrary qubit network.
    bonds: list of (i, j, J) tuples."""
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for i, j, J in bonds:
        for p in [X, Y, Z]:
            H += J * op_at(p,i,n_q) @ op_at(p,j,n_q)
    L_ops = [np.sqrt(gamma)*op_at(Z,k,n_q) for k in range(n_q)]
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in L_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    return L_mat

def extract_decay_rates(L_mat, gamma, threshold=0.1):
    """Extract unique decay rates (in units of gamma) from oscillatory modes."""
    evals = np.linalg.eigvals(L_mat)
    rates = []
    for ev in evals:
        if abs(np.imag(ev)) > threshold:
            rate = -np.real(ev)
            if rate > 0.001:
                rates.append(round(rate / gamma, 4))
    return sorted(set(rates))

# ============================================================
# PART 1: 3-qubit reference (should be topology-independent)
# ============================================================
print("=" * 75)
print("4-QUBIT BREAKDOWN ANALYSIS")
print("How does frequency-decay orthogonality break?")
print("=" * 75)

gamma = 0.05

print(f"\n--- PART 1: 3-qubit reference (gamma={gamma}) ---")
print(f"  Expected: all rates are {2}, {8/3:.4f}, {10/3:.4f} regardless of J")
print(f"\n  {'J_SA':>5} {'J_SB':>5} | decay rates (in gamma units)")
print("-" * 60)
for J_SA, J_SB in [(1,1), (1,2), (0.5,3), (2,2), (1,5), (0.3,0.7)]:
    bonds = [(0,1,J_SA), (0,2,J_SB)]  # S=0, A=1, B=2
    L = build_liouvillian(3, bonds, gamma)
    rates = extract_decay_rates(L, gamma)
    print(f"  {J_SA:>5.1f} {J_SB:>5.1f} | {rates}")

# ============================================================
# PART 2: 4-qubit star - sweep one coupling
# ============================================================
print(f"\n--- PART 2: 4-qubit STAR S(0)-A(1), S(0)-B(2), S(0)-C(3) ---")
print(f"  Fix J_SA=1, J_SB=2, sweep J_SC from 0.1 to 5.0")
print(f"\n  {'J_SC':>5} | decay rates (gamma units)")
print("-" * 70)

j_sc_values = [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
star_4q_data = []
for J_SC in j_sc_values:
    bonds = [(0,1,1.0), (0,2,2.0), (0,3,J_SC)]
    L = build_liouvillian(4, bonds, gamma)
    rates = extract_decay_rates(L, gamma)
    star_4q_data.append((J_SC, rates))
    rates_str = ", ".join(f"{r:.4f}" for r in rates)
    print(f"  {J_SC:>5.2f} | {rates_str}")

# Identify which rates are FIXED and which MOVE
print(f"\n--- FIXED vs MOVING rates (4-qubit star) ---")
if star_4q_data:
    # Find rates that appear at every J_SC value
    all_rate_sets = [set(r for r in rates) for _, rates in star_4q_data]
    # Check each rate from first config
    first_rates = star_4q_data[0][1]
    for ref_rate in first_rates:
        appearances = 0
        max_dev = 0
        for _, rates in star_4q_data:
            closest = min(rates, key=lambda r: abs(r - ref_rate))
            if abs(closest - ref_rate) < 0.05:
                appearances += 1
                max_dev = max(max_dev, abs(closest - ref_rate))
        status = "FIXED" if max_dev < 0.001 else f"MOVES (max dev={max_dev:.4f})"
        print(f"  rate ~ {ref_rate:.4f}g: {status} (appeared in {appearances}/{len(star_4q_data)})")

# ============================================================
# PART 3: 4-qubit chain - same sweep
# ============================================================
print(f"\n--- PART 3: 4-qubit CHAIN 0-1-2-3 ---")
print(f"  Fix J_01=1, J_12=1.5, sweep J_23 from 0.1 to 5.0")
print(f"\n  {'J_23':>5} | decay rates (gamma units)")
print("-" * 70)

chain_4q_data = []
for J_23 in j_sc_values:
    bonds = [(0,1,1.0), (1,2,1.5), (2,3,J_23)]
    L = build_liouvillian(4, bonds, gamma)
    rates = extract_decay_rates(L, gamma)
    chain_4q_data.append((J_23, rates))
    rates_str = ", ".join(f"{r:.4f}" for r in rates)
    print(f"  {J_23:>5.2f} | {rates_str}")

# ============================================================
# PART 4: Do the MOVING rates at least scale linearly with gamma?
# ============================================================
print(f"\n--- PART 4: Gamma scaling test (4-qubit star, J=[1,2,1.5]) ---")
print(f"  If rates still scale with gamma, they are 'weakly broken'")
print(f"\n  {'gamma':>7} | decay rates (gamma units) | {'stable?':>7}")
print("-" * 70)

ref_rates = None
for g in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
    bonds = [(0,1,1.0), (0,2,2.0), (0,3,1.5)]
    L = build_liouvillian(4, bonds, g)
    rates = extract_decay_rates(L, g)
    rates_str = ", ".join(f"{r:.4f}" for r in rates)
    
    if ref_rates is None:
        ref_rates = rates
        print(f"  {g:>7.3f} | {rates_str} | REF")
    else:
        # Compare with reference
        max_dev = 0
        for r_ref, r_now in zip(ref_rates[:min(len(ref_rates),len(rates))], rates):
            max_dev = max(max_dev, abs(r_ref - r_now))
        stable = "YES" if max_dev < 0.01 else f"NO ({max_dev:.4f})"
        print(f"  {g:>7.3f} | {rates_str} | {stable}")

# ============================================================
# PART 5: Boundary values - what stays at every N?
# ============================================================
print(f"\n--- PART 5: N-qubit comparison (star topology, uniform J=1) ---")
print(f"  Looking for rates that survive at every system size")
print(f"\n  {'N':>3} | decay rates (gamma units)")
print("-" * 70)

for n_q in [2, 3, 4, 5]:
    if n_q == 2:
        bonds = [(0,1,1.0)]
    else:
        bonds = [(0,i,1.0) for i in range(1, n_q)]
    L = build_liouvillian(n_q, bonds, gamma)
    rates = extract_decay_rates(L, gamma)
    rates_str = ", ".join(f"{r:.4f}" for r in rates)
    print(f"  {n_q:>3} | {rates_str}")

# ============================================================
# PART 6: High-resolution J sweep to see rate trajectories
# ============================================================
print(f"\n--- PART 6: Rate trajectories (4-qubit star, fine J_SC sweep) ---")
print(f"  J_SA=1, J_SB=2, J_SC sweeps 0.01 to 3.0 in 30 steps")
print(f"  Tracking individual rate bands")
print(f"\n  {'J_SC':>6} | rates (gamma units, sorted)")
print("-" * 75)

fine_data = []
for J_SC in np.linspace(0.01, 3.0, 30):
    bonds = [(0,1,1.0), (0,2,2.0), (0,3,J_SC)]
    L = build_liouvillian(4, bonds, gamma)
    rates = extract_decay_rates(L, gamma)
    fine_data.append((J_SC, rates))
    rates_str = ", ".join(f"{r:.3f}" for r in rates[:8])
    print(f"  {J_SC:>6.3f} | {rates_str}")

# Identify the boundary rates (min and max)
print(f"\n--- BOUNDARY ANALYSIS ---")
all_min = min(min(r) for _, r in fine_data)
all_max = max(max(r) for _, r in fine_data)
print(f"  Smallest rate seen:  {all_min:.4f}g")
print(f"  Largest rate seen:   {all_max:.4f}g")
print(f"  3-qubit boundaries:  2.0000g and {10/3:.4f}g")
print(f"  Rate 2g present at 4q? {any(any(abs(r-2.0)<0.01 for r in rates) for _,rates in fine_data)}")

print(f"\n{'='*75}")
print("4-QUBIT BREAKDOWN ANALYSIS COMPLETE")
print("="*75)
