"""
MIRROR BREAKING TRANSITION: Where does 1/4 appear?
====================================================
Tom's prediction: the transition from perfect mirrors
to broken mirrors has a critical point at 1/4.

We blend: L = (1-alpha)*L_dephasing + alpha*L_amplitude_damping
alpha=0: pure dephasing (mirrors perfect)
alpha=1: pure amplitude damping (mirrors broken)
Somewhere in between: the mirrors break. WHERE?

Also test: does the symmetry score pass through 0.25 or 0.75?
Does the critical alpha equal 0.25?

Results -> results_mirror_transition.txt
"""
import numpy as np
import math
import time

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
sigma_minus = np.array([[0,0],[1,0]], dtype=complex)

def op_at(op, qubit, n_q):
    ops = [I2]*n_q
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def build_liouvillian_general(n_q, bonds, jump_ops):
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for i, j, J, paulis in bonds:
        for p in paulis:
            H += J * op_at(p,i,n_q) @ op_at(p,j,n_q)
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in jump_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    return L_mat

def heisenberg_bonds(n_q, couplings):
    return [(0, i, couplings[i-1], [X,Y,Z]) for i in range(1, n_q)]

def get_osc_rates_raw(L_mat, threshold=0.05):
    evals = np.linalg.eigvals(L_mat)
    rates = []
    for ev in evals:
        if abs(np.imag(ev)) > threshold:
            rate = -np.real(ev)
            if rate > 0.0001:
                rates.append(round(rate, 6))
    return sorted(rates)

def check_mirror_symmetry(rates, center, tol=0.005):
    below = [r for r in rates if r < center - tol]
    above = [r for r in rates if r > center + tol]
    matched = 0
    for r in below:
        mirror = 2*center - r
        closest = min(above, key=lambda a: abs(a - mirror)) if above else 999
        if abs(closest - mirror) < tol:
            matched += 1
    score = matched / max(len(below), 1) if below else 1.0
    return score

def best_symmetry_center(rates):
    """Find the center that maximizes mirror symmetry."""
    if len(rates) < 2:
        return rates[0] if rates else 0, 1.0
    mid = (min(rates) + max(rates)) / 2
    best_c = mid
    best_s = check_mirror_symmetry(rates, mid)
    # Search around midpoint
    for offset in np.linspace(-0.5*(max(rates)-min(rates)), 0.5*(max(rates)-min(rates)), 200):
        c = mid + offset
        s = check_mirror_symmetry(rates, c)
        if s > best_s:
            best_s = s
            best_c = c
    return best_c, best_s

outfile = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results_mirror_transition.txt"
f = open(outfile, "w")

def log(msg):
    print(msg)
    f.write(msg + "\n")
    f.flush()

log("=" * 80)
log("MIRROR BREAKING TRANSITION")
log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("Does 1/4 appear at the transition?")
log("=" * 80)

# ============================================================
# TEST 1: Blend dephasing -> amplitude damping
# ============================================================
log("\n### TEST 1: Blend dephasing -> amplitude damping")
log("L_jump = sqrt((1-a)*g)*Z + sqrt(a*g)*sigma_minus per qubit")
log("alpha=0: pure dephasing. alpha=1: pure amplitude damping.")
log("High resolution: 100 steps")

gamma = 0.05
alpha_values = np.linspace(0.0, 1.0, 101)

for n_q in [3, 4, 5]:
    log(f"\n  N={n_q}:")
    log(f"  {'alpha':>7} | {'sym_Ng':>7} | {'sym_mid':>7} | {'best_sym':>8} | "
        f"{'best_ctr':>8} | {'min_r':>7} | {'max_r':>7} | {'mid':>7}")
    log("-" * 85)
    
    transition_data = []
    for alpha in alpha_values:
        bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
        jumps = []
        for k in range(n_q):
            if alpha < 1.0:
                jumps.append(np.sqrt((1-alpha)*gamma) * op_at(Z, k, n_q))
            if alpha > 0.0:
                jumps.append(np.sqrt(alpha*gamma) * op_at(sigma_minus, k, n_q))
        
        L = build_liouvillian_general(n_q, bonds, jumps)
        rates = get_osc_rates_raw(L)
        
        if not rates:
            continue
        
        center_Ng = n_q * gamma
        sym_Ng = check_mirror_symmetry(rates, center_Ng)
        mid = (min(rates) + max(rates)) / 2
        sym_mid = check_mirror_symmetry(rates, mid)
        best_c, best_s = best_symmetry_center(rates)
        
        transition_data.append({
            'alpha': alpha, 'sym_Ng': sym_Ng, 'sym_mid': sym_mid,
            'best_sym': best_s, 'best_center': best_c,
            'min_r': min(rates), 'max_r': max(rates), 'mid': mid
        })
        
        if alpha * 100 % 5 < 0.5:  # print every 5%
            log(f"  {alpha:>7.3f} | {sym_Ng:>7.1%} | {sym_mid:>7.1%} | {best_s:>8.1%} | "
                f"{best_c:>8.5f} | {min(rates):>7.5f} | {max(rates):>7.5f} | {mid:>7.5f}")
    
    # Find the critical alpha where Ng-symmetry drops below 100%
    if transition_data:
        critical_alpha = None
        for d in transition_data:
            if d['sym_Ng'] < 0.99:
                critical_alpha = d['alpha']
                break
        
        if critical_alpha is not None:
            log(f"\n  CRITICAL ALPHA (Ng-symmetry first breaks): {critical_alpha:.4f}")
            log(f"  Is it 0.25? Difference: {abs(critical_alpha - 0.25):.4f}")
        else:
            log(f"\n  Ng-symmetry never broke!")
        
        # Find where best_symmetry drops below various thresholds
        for threshold in [0.99, 0.95, 0.90, 0.75, 0.50, 0.25]:
            for d in transition_data:
                if d['best_sym'] < threshold:
                    log(f"  Best symmetry drops below {threshold:.0%} at alpha={d['alpha']:.4f}")
                    break
        
        # Check: at alpha=0.25, what is the symmetry?
        for d in transition_data:
            if abs(d['alpha'] - 0.25) < 0.006:
                log(f"\n  AT ALPHA=0.25: Ng_sym={d['sym_Ng']:.1%}, mid_sym={d['sym_mid']:.1%}, "
                    f"best_sym={d['best_sym']:.1%}")
                break
        
        # Check: does any rate or center value equal N*gamma/4 or gamma/4?
        for d in transition_data:
            if abs(d['alpha'] - 0.25) < 0.006:
                log(f"  Rates at alpha=0.25: [{d['min_r']:.5f}, {d['max_r']:.5f}]")
                log(f"  Center at alpha=0.25: {d['best_center']:.5f}")
                log(f"  N*gamma/4 = {n_q*gamma/4:.5f}")
                log(f"  gamma*N*(1-0.25) = {gamma*n_q*0.75:.5f}")

# ============================================================
# TEST 2: Different blend: dephasing -> depolarizing
# ============================================================
log("\n### TEST 2: Blend dephasing -> depolarizing")
log("alpha=0: Z dephasing only. alpha=1: full depolarizing (X+Y+Z)")

n_q = 4
gamma = 0.05

log(f"\n  N={n_q}:")
log(f"  {'alpha':>7} | {'sym_Ng':>7} | {'sym_mid':>7} | {'best_sym':>8}")
log("-" * 45)

for alpha in np.linspace(0.0, 1.0, 51):
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = []
    for k in range(n_q):
        # Z component always present
        jumps.append(np.sqrt(gamma) * op_at(Z, k, n_q))
        # X and Y components grow with alpha
        if alpha > 0:
            jumps.append(np.sqrt(alpha*gamma) * op_at(X, k, n_q))
            jumps.append(np.sqrt(alpha*gamma) * op_at(Y, k, n_q))
    
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    
    if rates:
        sym_Ng = check_mirror_symmetry(rates, n_q*gamma)
        mid = (min(rates)+max(rates))/2
        sym_mid = check_mirror_symmetry(rates, mid)
        best_c, best_s = best_symmetry_center(rates)
        
        if alpha * 50 % 5 < 0.5:
            log(f"  {alpha:>7.3f} | {sym_Ng:>7.1%} | {sym_mid:>7.1%} | {best_s:>8.1%}")

# ============================================================
# TEST 3: Look for 1/4 in the rates themselves
# ============================================================
log("\n### TEST 3: Does 1/4 appear in the rate values?")
log("At the transition, do any rates equal exactly 1/4 of something?")

n_q = 4
gamma = 0.05

for alpha in [0.0, 0.10, 0.20, 0.25, 0.30, 0.50, 0.75, 1.0]:
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = []
    for k in range(n_q):
        if alpha < 1.0:
            jumps.append(np.sqrt((1-alpha)*gamma) * op_at(Z, k, n_q))
        if alpha > 0.0:
            jumps.append(np.sqrt(alpha*gamma) * op_at(sigma_minus, k, n_q))
    
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    
    if rates:
        unique_r = sorted(set(round(r, 4) for r in rates))
        max_r = max(rates)
        log(f"\n  alpha={alpha:.2f}: {len(unique_r)} unique rates")
        log(f"    Range: [{min(rates):.5f}, {max_r:.5f}]")
        log(f"    max/4 = {max_r/4:.5f}")
        
        # Check ratios between rates
        ratios_quarter = []
        for r in unique_r:
            if abs(r/max_r - 0.25) < 0.01:
                ratios_quarter.append(r)
        if ratios_quarter:
            log(f"    Rates at 1/4 of max: {ratios_quarter}")
        
        # Check: does min_rate / max_rate = 1/4?
        ratio = min(rates) / max_r
        log(f"    min/max = {ratio:.5f} (1/4 = 0.25, diff = {abs(ratio-0.25):.5f})")

# ============================================================
# TEST 4: Bernoulli connection - is sym_score = 4*alpha*(1-alpha)?
# ============================================================
log("\n### TEST 4: Is the symmetry score related to Bernoulli variance?")
log("If sym = 1 - 4*alpha*(1-alpha), the 1/4 boundary appears naturally")
log("Or: sym = (1-2*alpha)^2, which equals 0 at alpha=1/2 and 1 at alpha=0,1")

n_q = 4
gamma = 0.05

log(f"\n  {'alpha':>7} | {'sym_best':>8} | {'1-4a(1-a)':>10} | {'(1-2a)^2':>8} | {'4a(1-a)':>8}")
log("-" * 60)

for alpha in np.linspace(0.0, 1.0, 41):
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = []
    for k in range(n_q):
        if alpha < 1.0:
            jumps.append(np.sqrt((1-alpha)*gamma) * op_at(Z, k, n_q))
        if alpha > 0.0:
            jumps.append(np.sqrt(alpha*gamma) * op_at(sigma_minus, k, n_q))
    
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    
    if rates:
        best_c, best_s = best_symmetry_center(rates)
        bernoulli = 4*alpha*(1-alpha)
        comp_bernoulli = 1 - bernoulli
        sq = (1-2*alpha)**2
        log(f"  {alpha:>7.3f} | {best_s:>8.4f} | {comp_bernoulli:>10.4f} | {sq:>8.4f} | {bernoulli:>8.4f}")

# ============================================================
# TEST 5: Center drift - how does the symmetry center move?
# ============================================================
log("\n### TEST 5: How does the symmetry center move during transition?")
log("Does it cross N*gamma/4 or N*gamma*(1-alpha)?")

n_q = 4
gamma = 0.05

log(f"\n  {'alpha':>7} | {'best_ctr':>8} | {'Ng':>8} | {'ctr/Ng':>8} | {'1-a':>6}")
log("-" * 50)

for alpha in np.linspace(0.0, 1.0, 41):
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = []
    for k in range(n_q):
        if alpha < 1.0:
            jumps.append(np.sqrt((1-alpha)*gamma) * op_at(Z, k, n_q))
        if alpha > 0.0:
            jumps.append(np.sqrt(alpha*gamma) * op_at(sigma_minus, k, n_q))
    
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    
    if rates:
        best_c, best_s = best_symmetry_center(rates)
        Ng = n_q * gamma
        ratio = best_c / Ng if Ng > 0 else 0
        log(f"  {alpha:>7.3f} | {best_c:>8.5f} | {Ng:>8.4f} | {ratio:>8.4f} | {1-alpha:>6.3f}")

# ============================================================
# TEST 6: N-dependence of critical alpha
# ============================================================
log("\n### TEST 6: Does the critical alpha depend on N?")
log("If 1/4 is universal, it should be the same at every N")

gamma = 0.05

for n_q in [3, 4, 5]:
    log(f"\n  N={n_q}:")
    for alpha in np.linspace(0.0, 0.5, 51):
        bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
        jumps = []
        for k in range(n_q):
            if alpha < 1.0:
                jumps.append(np.sqrt((1-alpha)*gamma) * op_at(Z, k, n_q))
            if alpha > 0.0:
                jumps.append(np.sqrt(alpha*gamma) * op_at(sigma_minus, k, n_q))
        
        L = build_liouvillian_general(n_q, bonds, jumps)
        rates = get_osc_rates_raw(L)
        
        if rates:
            sym_Ng = check_mirror_symmetry(rates, n_q*gamma)
            if sym_Ng < 0.99:
                log(f"    Ng-symmetry first breaks at alpha = {alpha:.4f}")
                log(f"    Is it 1/4? diff = {abs(alpha - 0.25):.4f}")
                break
    else:
        log(f"    Ng-symmetry never broke in range [0, 0.5]")

log("\n" + "=" * 80)
log("COMPLETE")
log("=" * 80)
f.close()
print(f"\n>>> Results saved to: {outfile}")
