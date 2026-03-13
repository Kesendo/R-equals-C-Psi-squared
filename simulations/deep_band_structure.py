"""
DEEP BAND STRUCTURE: N=2 through N=6
=====================================
Run by Tom, results read by Claude.
Outputs to: results_band_structure.txt

Tests:
1. Rate spectrum at each N (uniform J=1 star)
2. J-sweep at each N to find band widths
3. Avoided crossing gap vs N
4. Rate density scaling
5. Search for second transition (like 3->4 but at higher N)
"""
import numpy as np
import json
import time
import sys

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
    evals = np.linalg.eigvals(L_mat)
    rates = []
    for ev in evals:
        if abs(np.imag(ev)) > threshold:
            rate = -np.real(ev) / gamma
            if rate > 0.01:
                rates.append(round(rate, 4))
    return sorted(rates)

def star_bonds(n_q, couplings):
    """Star topology: qubit 0 is center, rest are leaves."""
    return [(0, i, couplings[i-1]) for i in range(1, n_q)]

gamma = 0.05
outfile = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results_band_structure.txt"
f = open(outfile, "w")

def log(msg):
    print(msg)
    f.write(msg + "\n")
    f.flush()

log("=" * 80)
log("DEEP BAND STRUCTURE ANALYSIS")
log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 80)

# ============================================================
# TEST 1: Uniform star (J=1 on all arms) at each N
# ============================================================
log("\n### TEST 1: Uniform star topology (all J=1)")
log("Liouvillian size: N=2 (16), N=3 (64), N=4 (256), N=5 (1024), N=6 (4096)")

max_n = 6  # Set to 5 if N=6 is too slow (4096x4096 matrix)

for n_q in range(2, max_n + 1):
    t0 = time.time()
    log(f"\n  N={n_q} (matrix {(2**n_q)**2}x{(2**n_q)**2})...")
    
    bonds = star_bonds(n_q, [1.0]*(n_q-1))
    try:
        L = build_liouvillian(n_q, bonds, gamma)
        rates = get_osc_rates(L, gamma)
        elapsed = time.time() - t0
        
        log(f"    Time: {elapsed:.1f}s")
        log(f"    Oscillatory rates: {len(rates)}")
        log(f"    Min rate: {min(rates):.4f}g")
        log(f"    Max rate: {max(rates):.4f}g")
        log(f"    Bandwidth: {max(rates)-min(rates):.4f}g")
        log(f"    Predicted max 2(N-1): {2*(n_q-1)}g")
        log(f"    Max matches 2(N-1): {abs(max(rates) - 2*(n_q-1)) < 0.01}")
        
        # Count distinct rate values (cluster within 0.005)
        unique = []
        for r in rates:
            if not unique or abs(r - unique[-1]) > 0.005:
                unique.append(r)
        log(f"    Distinct rate clusters: {len(unique)}")
        
        # Rate density
        if max(rates) > min(rates):
            density = len(rates) / (max(rates) - min(rates))
            log(f"    Rate density: {density:.1f} rates/g")
        
        # All rates (for Claude to analyze)
        log(f"    ALL RATES: {rates}")
        
    except Exception as e:
        log(f"    FAILED: {e}")
        log(f"    (N={n_q} may be too large for this machine)")
        break

# ============================================================
# TEST 2: J-sweep at each N to measure band widths
# ============================================================
log("\n### TEST 2: J-sweep to measure how much rates move")
log("Fix all J=1 except last arm, sweep last from 0.1 to 3.0")

for n_q in range(3, min(max_n+1, 6)):  # skip N=6 for sweep
    log(f"\n  N={n_q} star, sweeping J_last...")
    t0 = time.time()
    
    n_sweep = 20 if n_q <= 4 else 10
    j_values = np.linspace(0.1, 3.0, n_sweep)
    
    all_rates_seen = set()
    rate_at_each_j = []
    
    for J_last in j_values:
        couplings = [1.0]*(n_q-2) + [J_last]
        bonds = star_bonds(n_q, couplings)
        try:
            L = build_liouvillian(n_q, bonds, gamma)
            rates = get_osc_rates(L, gamma)
            all_rates_seen.update(rates)
            rate_at_each_j.append((J_last, rates))
        except:
            break
    
    elapsed = time.time() - t0
    log(f"    Time: {elapsed:.1f}s ({len(rate_at_each_j)} points)")
    
    if all_rates_seen:
        log(f"    Total unique rates across sweep: {len(all_rates_seen)}")
        log(f"    Range: [{min(all_rates_seen):.4f}, {max(all_rates_seen):.4f}]")
        
        # Find rates that NEVER move (appear at same value in every config)
        if rate_at_each_j:
            first_rates = set(rate_at_each_j[0][1])
            fixed_rates = []
            for ref in sorted(first_rates):
                is_fixed = True
                for _, rates in rate_at_each_j:
                    closest = min(rates, key=lambda r: abs(r - ref)) if rates else 999
                    if abs(closest - ref) > 0.01:
                        is_fixed = False
                        break
                if is_fixed:
                    fixed_rates.append(ref)
            
            log(f"    Fixed rates (never move): {[f'{r:.4f}' for r in fixed_rates]}")
            log(f"    Moving rates: {len(all_rates_seen) - len(fixed_rates)*len(rate_at_each_j)}")

# ============================================================
# TEST 3: Avoided crossings - minimum gap between rates vs N
# ============================================================
log("\n### TEST 3: Avoided crossings - do rates ever cross?")
log("Fine J-sweep, track minimum gap between adjacent rates")

for n_q in range(3, min(max_n+1, 6)):
    log(f"\n  N={n_q} star...")
    t0 = time.time()
    
    n_fine = 40 if n_q <= 4 else 15
    j_fine = np.linspace(0.3, 2.5, n_fine)
    
    min_gap_ever = float('inf')
    min_gap_at_j = 0
    all_gaps = []
    
    for J_last in j_fine:
        couplings = [1.0]*(n_q-2) + [J_last]
        bonds = star_bonds(n_q, couplings)
        try:
            L = build_liouvillian(n_q, bonds, gamma)
            rates = get_osc_rates(L, gamma)
            if len(rates) > 1:
                gaps = [rates[i+1] - rates[i] for i in range(len(rates)-1)]
                real_gaps = [g for g in gaps if g > 0.0005]
                if real_gaps:
                    mg = min(real_gaps)
                    all_gaps.append(mg)
                    if mg < min_gap_ever:
                        min_gap_ever = mg
                        min_gap_at_j = J_last
        except:
            break
    
    elapsed = time.time() - t0
    if all_gaps:
        avg_gap = np.mean(all_gaps)
        log(f"    Time: {elapsed:.1f}s")
        log(f"    Smallest gap: {min_gap_ever:.6f}g at J_last={min_gap_at_j:.3f}")
        log(f"    Average min gap: {avg_gap:.6f}g")
        log(f"    Gap closes to zero: {min_gap_ever < 0.0005}")
        if min_gap_ever > 0.0005:
            log(f"    -> AVOIDED CROSSINGS (topological protection)")
        else:
            log(f"    -> REAL CROSSINGS detected!")

# ============================================================
# TEST 4: Symmetry of band structure around center
# ============================================================
log("\n### TEST 4: Is the band structure symmetric around Ng?")
log("At uniform J=1, center = (min+max)/2 = (2 + 2(N-1))/2 = Ng")

for n_q in range(3, min(max_n+1, 6)):
    bonds = star_bonds(n_q, [1.0]*(n_q-1))
    try:
        L = build_liouvillian(n_q, bonds, gamma)
        rates = get_osc_rates(L, gamma)
        center = n_q  # center = Ng = (2 + 2(N-1))/2
        
        # For each rate below center, check if mirror exists above
        below = [r for r in rates if r < center - 0.01]
        above = [r for r in rates if r > center + 0.01]
        at_center = [r for r in rates if abs(r - center) < 0.01]
        
        matched = 0
        unmatched_below = []
        for r in below:
            mirror = 2*center - r
            closest = min(above, key=lambda a: abs(a - mirror)) if above else 999
            if abs(closest - mirror) < 0.02:
                matched += 1
            else:
                unmatched_below.append(r)
        
        log(f"\n  N={n_q}: center={center}g")
        log(f"    Rates below center: {len(below)}")
        log(f"    Rates above center: {len(above)}")
        log(f"    Rates at center: {len(at_center)}")
        log(f"    Mirror-matched pairs: {matched}")
        log(f"    Unmatched: {len(unmatched_below)}")
        symmetry = matched / max(len(below), 1)
        log(f"    Symmetry score: {symmetry:.1%}")
    except:
        log(f"  N={n_q}: FAILED")

# ============================================================
# TEST 5: Chain vs Star topology comparison
# ============================================================
log("\n### TEST 5: Chain vs Star at same N")
log("Same number of qubits, different connectivity")

for n_q in range(3, min(max_n+1, 6)):
    # Star: 0 connects to all others
    star_b = star_bonds(n_q, [1.0]*(n_q-1))
    # Chain: 0-1-2-3-...
    chain_b = [(i, i+1, 1.0) for i in range(n_q-1)]
    
    try:
        L_star = build_liouvillian(n_q, star_b, gamma)
        L_chain = build_liouvillian(n_q, chain_b, gamma)
        
        r_star = get_osc_rates(L_star, gamma)
        r_chain = get_osc_rates(L_chain, gamma)
        
        log(f"\n  N={n_q}:")
        log(f"    STAR:  {len(r_star)} rates, range [{min(r_star):.4f}, {max(r_star):.4f}]")
        log(f"    CHAIN: {len(r_chain)} rates, range [{min(r_chain):.4f}, {max(r_chain):.4f}]")
        log(f"    Same min? {abs(min(r_star)-min(r_chain)) < 0.01}")
        log(f"    Same max? {abs(max(r_star)-max(r_chain)) < 0.01}")
        
        # Do they share the same boundary rates?
        star_min = min(r_star)
        chain_min = min(r_chain)
        star_max = max(r_star)
        chain_max = max(r_chain)
        log(f"    Star min={star_min:.4f}, Chain min={chain_min:.4f}")
        log(f"    Star max={star_max:.4f}, Chain max={chain_max:.4f}")
    except:
        log(f"  N={n_q}: FAILED")

# ============================================================
# TEST 6: Scaling laws summary
# ============================================================
log("\n### TEST 6: Scaling summary table")
log(f"\n  {'N':>3} | {'Liouv':>8} | {'#rates':>7} | {'min':>6} | {'max':>6} | {'BW':>6} | {'2(N-1)':>7} | {'max=2(N-1)':>11}")
log("-" * 75)

# Collect from what we computed
for n_q in range(2, max_n + 1):
    bonds = star_bonds(n_q, [1.0]*(n_q-1))
    try:
        L = build_liouvillian(n_q, bonds, gamma)
        rates = get_osc_rates(L, gamma)
        lo = min(rates) if rates else 0
        hi = max(rates) if rates else 0
        bw = hi - lo
        pred = 2 * (n_q - 1)
        match = "YES" if abs(hi - pred) < 0.01 else "NO"
        liouv_size = (2**n_q)**2
        log(f"  {n_q:>3} | {liouv_size:>8} | {len(rates):>7} | {lo:>6.2f} | {hi:>6.2f} | {bw:>6.2f} | {pred:>7} | {match:>11}")
    except:
        log(f"  {n_q:>3} | FAILED")

log(f"\nCompleted: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 80)
f.close()

print(f"\n>>> Results saved to: {outfile}")
print(f">>> Copy the contents and paste to Claude for analysis.")
