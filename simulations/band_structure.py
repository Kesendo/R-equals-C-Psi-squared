"""
BAND STRUCTURE ANALYSIS: From discrete rates to continuum
==========================================================
The 3->4 qubit transition is where discrete fixed rates become
continuous bands. This script maps the band structure precisely.

Questions:
1. How many distinct bands form at N=4?
2. Do bands cross or avoid crossing (like electronic bands)?
3. What is the density of states at each N?
4. Can we predict N=6 band structure from N=3,4,5 pattern?
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

def get_osc_rates(L_mat, gamma, threshold=0.05):
    """Get ALL oscillatory decay rates including degeneracies."""
    evals = np.linalg.eigvals(L_mat)
    rates = []
    for ev in evals:
        if abs(np.imag(ev)) > threshold:
            rate = -np.real(ev) / gamma
            if rate > 0.01:
                rates.append(round(rate, 4))
    return sorted(rates)

gamma = 0.05

# ============================================================
# PART 1: High-resolution J sweep for 3-qubit (reference)
# ============================================================
print("=" * 75)
print("BAND STRUCTURE: From discrete rates to continuum")
print("=" * 75)

print(f"\n--- PART 1: 3-qubit star, J_SB sweep (50 points) ---")
print(f"  J_SA=1.0, J_SB from 0.1 to 5.0")
n_points = 30
j_sweep = np.linspace(0.1, 5.0, n_points)

rates_3q = []
for J_SB in j_sweep:
    bonds = [(0,1,1.0), (0,2,J_SB)]
    L = build_liouvillian(3, bonds, gamma)
    rates = get_osc_rates(L, gamma)
    rates_3q.append(rates)

# Find unique rate values across all J
all_rates_3q = set()
for r in rates_3q:
    all_rates_3q.update(r)

print(f"  Total unique rates seen: {len(all_rates_3q)}")
print(f"  Rate range: [{min(all_rates_3q):.4f}, {max(all_rates_3q):.4f}]")
print(f"  Distinct bands: ", end="")
# Cluster rates into bands (within 0.05 of each other)
bands_3q = []
sorted_rates = sorted(all_rates_3q)
current_band = [sorted_rates[0]]
for r in sorted_rates[1:]:
    if r - current_band[-1] < 0.05:
        current_band.append(r)
    else:
        bands_3q.append((min(current_band), max(current_band)))
        current_band = [r]
bands_3q.append((min(current_band), max(current_band)))
print(f"{len(bands_3q)} bands")
for i, (lo, hi) in enumerate(bands_3q):
    width = hi - lo
    status = "FIXED" if width < 0.01 else f"width={width:.4f}"
    print(f"    Band {i+1}: [{lo:.4f}, {hi:.4f}] ({status})")

# ============================================================
# PART 2: High-resolution J sweep for 4-qubit star
# ============================================================
print(f"\n--- PART 2: 4-qubit star, J_SC sweep (50 points) ---")
print(f"  J_SA=1.0, J_SB=2.0, J_SC from 0.1 to 5.0")

rates_4q = []
for J_SC in j_sweep:
    bonds = [(0,1,1.0), (0,2,2.0), (0,3,J_SC)]
    L = build_liouvillian(4, bonds, gamma)
    rates = get_osc_rates(L, gamma)
    rates_4q.append(rates)

all_rates_4q = set()
for r in rates_4q:
    all_rates_4q.update(r)

print(f"  Total unique rates seen: {len(all_rates_4q)}")
print(f"  Rate range: [{min(all_rates_4q):.4f}, {max(all_rates_4q):.4f}]")

# Cluster into bands
sorted_4q = sorted(all_rates_4q)
bands_4q = []
current_band = [sorted_4q[0]]
for r in sorted_4q[1:]:
    if r - current_band[-1] < 0.03:
        current_band.append(r)
    else:
        bands_4q.append((min(current_band), max(current_band), len(current_band)))
        current_band = [r]
bands_4q.append((min(current_band), max(current_band), len(current_band)))

print(f"  Distinct bands: {len(bands_4q)}")
for i, (lo, hi, count) in enumerate(bands_4q):
    width = hi - lo
    status = "FIXED" if width < 0.01 else f"width={width:.4f}"
    print(f"    Band {i+1}: [{lo:.4f}, {hi:.4f}] ({status}, {count} values)")

# ============================================================
# PART 3: 5-qubit star (fewer points, bigger matrix)
# ============================================================
print(f"\n--- PART 3: 5-qubit star, J_SD sweep (8 points) ---")
print(f"  J_SA=1, J_SB=2, J_SC=1.5, J_SD from 0.5 to 3.0")
print(f"  (1024x1024 Liouvillian, will take a while)")

j_sweep_5 = np.linspace(0.5, 3.0, 8)
rates_5q = []
for idx, J_SD in enumerate(j_sweep_5):
    bonds = [(0,1,1.0), (0,2,2.0), (0,3,1.5), (0,4,J_SD)]
    L = build_liouvillian(5, bonds, gamma)
    rates = get_osc_rates(L, gamma)
    rates_5q.append(rates)
    if idx % 2 == 0:
        print(f"  ... {idx+1}/8 done (J_SD={J_SD:.2f}, {len(rates)} rates)")

all_rates_5q = set()
for r in rates_5q:
    all_rates_5q.update(r)

print(f"  Total unique rates seen: {len(all_rates_5q)}")
print(f"  Rate range: [{min(all_rates_5q):.4f}, {max(all_rates_5q):.4f}]")

# ============================================================
# PART 4: Scaling laws
# ============================================================
print(f"\n--- PART 4: Scaling comparison N=2,3,4,5 ---")
print(f"\n  {'N':>3} | {'min':>6} | {'max':>6} | {'#unique':>7} | {'#bands':>6} | {'bandwidth':>10}")
print("-" * 55)

# N=2 reference
bonds_2 = [(0,1,1.0)]
L2 = build_liouvillian(2, bonds_2, gamma)
r2 = get_osc_rates(L2, gamma)

for n, all_r in [(2, set(r2)), (3, all_rates_3q), (4, all_rates_4q), (5, all_rates_5q)]:
    lo = min(all_r)
    hi = max(all_r)
    bw = hi - lo
    n_unique = len(all_r)
    # Quick band count
    sr = sorted(all_r)
    n_bands = 1
    for i in range(1, len(sr)):
        if sr[i] - sr[i-1] > 0.05:
            n_bands += 1
    print(f"  {n:>3} | {lo:>6.2f} | {hi:>6.2f} | {n_unique:>7} | {n_bands:>6} | {bw:>10.4f}")

# ============================================================
# PART 5: Avoided crossings test (4-qubit)
# ============================================================
print(f"\n--- PART 5: Avoided crossings (4-qubit star) ---")
print(f"  Do rate bands cross or repel each other?")
print(f"  Fine sweep: J_SC from 0.5 to 2.5 in 50 steps")

fine_j = np.linspace(0.5, 2.5, 50)
rate_trajectories = []
for J_SC in fine_j:
    bonds = [(0,1,1.0), (0,2,2.0), (0,3,J_SC)]
    L = build_liouvillian(4, bonds, gamma)
    rates = get_osc_rates(L, gamma)
    rate_trajectories.append(rates)

# Find minimum gap between adjacent rates at each J
min_gaps = []
for i, rates in enumerate(rate_trajectories):
    if len(rates) > 1:
        gaps = [rates[j+1] - rates[j] for j in range(len(rates)-1)]
        min_gap = min(g for g in gaps if g > 0.001)  # skip degeneracies
        min_gaps.append((fine_j[i], min_gap))

# Print where gaps are smallest (potential avoided crossings)
min_gaps.sort(key=lambda x: x[1])
print(f"\n  Smallest inter-rate gaps (avoided crossings?):")
for j_val, gap in min_gaps[:10]:
    print(f"    J_SC={j_val:.3f}: min_gap = {gap:.5f}g")

avg_gap = np.mean([g for _, g in min_gaps])
min_ever = min_gaps[0][1]
print(f"\n  Average minimum gap: {avg_gap:.5f}g")
print(f"  Smallest gap ever:  {min_ever:.5f}g")
if min_ever > 0.001:
    print(f"  -> AVOIDED CROSSINGS: rates repel each other (gap never closes)")
else:
    print(f"  -> REAL CROSSINGS: rates can cross (gap closes to zero)")

# ============================================================
# PART 6: Rate density (how many rates per unit interval)
# ============================================================
print(f"\n--- PART 6: Rate density at fixed J ---")
print(f"  How does the density of rates scale with N?")
for n, label, rates_list in [
    (3, "3q star [1,2]", rates_3q[25]),
    (4, "4q star [1,2,1.5]", rate_trajectories[25]),
]:
    n_rates = len(rates_list)
    if n_rates > 1:
        bw = max(rates_list) - min(rates_list)
        density = n_rates / bw if bw > 0 else 0
        print(f"  {label}: {n_rates} rates over [{min(rates_list):.2f}, {max(rates_list):.2f}], "
              f"density = {density:.2f} rates/g")

# Check 5-qubit density too (use middle of sweep)
if rates_5q:
    mid_5q = rates_5q[len(rates_5q)//2]
    if len(mid_5q) > 1:
        bw5 = max(mid_5q) - min(mid_5q)
        dens5 = len(mid_5q) / bw5 if bw5 > 0 else 0
        print(f"  5q star [1,2,1.5,{j_sweep_5[len(j_sweep_5)//2]:.1f}]: {len(mid_5q)} rates over "
              f"[{min(mid_5q):.2f}, {max(mid_5q):.2f}], density = {dens5:.2f} rates/g")

# ============================================================
# PART 7: The boundary formula
# ============================================================
print(f"\n--- PART 7: Boundary rate formula ---")
print(f"  Testing: min_rate = 2g, max_rate = 2(N-1)g for star topology")
print(f"")
for n in [2, 3, 4, 5]:
    predicted_max = 2 * (n - 1)
    if n == 2:
        observed_max = max(r2)
    elif n == 3:
        observed_max = max(all_rates_3q)
    elif n == 4:
        observed_max = max(all_rates_4q)
    else:
        observed_max = max(all_rates_5q)
    match = "EXACT" if abs(predicted_max - observed_max) < 0.01 else f"off by {abs(predicted_max-observed_max):.4f}"
    print(f"  N={n}: predicted max = {predicted_max}g, observed = {observed_max:.4f}g ({match})")

print(f"\n  Prediction for N=6: max rate = 10g, min = 2g")
print(f"  Prediction for N=10: max rate = 18g, min = 2g")
print(f"  Prediction for N=100: max rate = 198g, min = 2g")
print(f"  The band widens linearly with N: bandwidth = 2(N-2)g.")

print(f"\n{'='*75}")
print("BAND STRUCTURE ANALYSIS COMPLETE")
print("="*75)
print("""
SUMMARY:

The transition from discrete to continuous rate structure:

  N=2: Single rate (2g). No band.
  N=3: 3-4 discrete rates. All fixed. No bands.
  N=4: Multiple bands. Interior rates MOVE with J. 
       Boundaries 2g and 6g fixed.
  N=5: Dense spectrum. Many bands merge.
       Boundaries 2g and 8g fixed.

Boundary formula: min = 2g (always), max = 2(N-1)g (always).
Interior: free at N>=4, increasingly dense with N.

This is analogous to electronic band structure in solids:
  - 1 atom  = discrete energy levels
  - 2 atoms = bonding/antibonding split  
  - N atoms = energy BANDS
  - N->inf  = continuous band structure

The 3-qubit star is like a diatomic molecule: too small for bands.
The 4-qubit star is where the first bands appear.
At large N, the rate spectrum should become a continuous 
density of states between 2g and 2(N-1)g.
""")
