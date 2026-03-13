"""
WHAT IS z* IN THE DENSITY MATRIX?
=================================
z* = (1 - sqrt(1-4*CPsi))/2 satisfies z*(1-z*) = CPsi
Correlates with purity (r=0.917) but not exactly.

Goal: find the EXACT relationship between z* and rho_AB properties.
Compute z* alongside every known quantum quantity at each timestep.
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

def concurrence(rho_2q):
    """Wootters concurrence for 2-qubit state."""
    YY = np.kron(Y, Y)
    rho_t = YY @ rho_2q.conj() @ YY
    R = rho_2q @ rho_t
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, float(ev[0] - ev[1] - ev[2] - ev[3]))

def l1_coherence(rho, normalized=True):
    """l1-norm of coherence."""
    d = rho.shape[0]
    l1 = sum(abs(rho[i,j]) for i in range(d) for j in range(d) if i!=j)
    return float(l1) / (d-1) if normalized else float(l1)

def quantum_discord_simple(rho_2q):
    """Simplified quantum discord estimate via measurement on B."""
    # Total mutual information
    rho_A = np.array([[rho_2q[0,0]+rho_2q[1,1], rho_2q[0,2]+rho_2q[1,3]],
                      [rho_2q[2,0]+rho_2q[3,1], rho_2q[2,2]+rho_2q[3,3]]])
    rho_B = np.array([[rho_2q[0,0]+rho_2q[2,2], rho_2q[0,1]+rho_2q[2,3]],
                      [rho_2q[1,0]+rho_2q[3,2], rho_2q[1,1]+rho_2q[3,3]]])
    S_A = von_neumann_entropy(rho_A)
    S_B = von_neumann_entropy(rho_B)
    S_AB = von_neumann_entropy(rho_2q)
    mutual_info = S_A + S_B - S_AB
    # Classical correlation (Z-measurement on B)
    p0 = float(np.real(rho_2q[0,0] + rho_2q[2,2]))
    p1 = float(np.real(rho_2q[1,1] + rho_2q[3,3]))
    if p0 > 1e-10:
        rho_A_0 = np.array([[rho_2q[0,0], rho_2q[0,2]],
                            [rho_2q[2,0], rho_2q[2,2]]]) / p0
    else:
        rho_A_0 = np.eye(2)/2
    if p1 > 1e-10:
        rho_A_1 = np.array([[rho_2q[1,1], rho_2q[1,3]],
                            [rho_2q[3,1], rho_2q[3,3]]]) / p1
    else:
        rho_A_1 = np.eye(2)/2
    S_cond = p0*von_neumann_entropy(rho_A_0) + p1*von_neumann_entropy(rho_A_1)
    classical_corr = S_A - S_cond
    return max(0, mutual_info - classical_corr)

def von_neumann_entropy(rho):
    ev = np.real(np.linalg.eigvalsh(rho))
    ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def run_trajectory(J_SA=1.0, J_SB=2.0, gamma=0.05, dt=0.01, t_max=5.0):
    """Run dynamics, compute everything at each step."""
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
    rho = np.outer(psi0, psi0.conj())
    
    data = []
    for step in range(int(t_max/dt)+1):
        t = step * dt
        if step % 10 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], n_q)
            C = concurrence(rAB)
            psi_n = l1_coherence(rAB, normalized=True)
            cpsi = C * psi_n
            
            # z* from CPsi
            if cpsi > 0.001 and cpsi <= 0.25:
                zstar = (1 - math.sqrt(1-4*cpsi)) / 2
            elif cpsi > 0.25:
                zstar = 0.5
            else:
                zstar = 0.0
            
            # Eigenvalues of rho_AB (sorted descending)
            eigvals = np.sort(np.real(np.linalg.eigvalsh(rAB)))[::-1]
            lam1, lam2, lam3, lam4 = eigvals
            
            # Many candidate quantities
            purity = float(np.real(np.trace(rAB @ rAB)))
            lin_ent = (4/3)*(1 - purity)  # normalized linear entropy
            S_vN = von_neumann_entropy(rAB)
            
            # Eigenvalue-based candidates
            lam_max = lam1
            lam_gap = lam1 - lam2  # gap between two largest
            lam_spread = lam1 - lam4
            lam_prod12 = lam1 * lam2
            
            # Concurrence-related
            sqrt_cpsi = math.sqrt(cpsi) if cpsi > 0 else 0
            
            # Purity-based candidates for z*
            # If z* ~ f(purity), what f?
            # purity ranges from 0.25 (maximally mixed) to 1.0 (pure)
            # z* ranges from 0 to 0.5
            pur_shifted = purity - 0.25  # 0 to 0.75
            pur_norm = (purity - 0.25) / 0.75  # 0 to 1
            
            # Participation ratio
            part_ratio = 1.0 / sum(e**2 for e in eigvals if e > 1e-15)
            
            # Negativity (entanglement monotone)
            # Partial transpose of rho_AB w.r.t. B
            rAB_pt = rAB.copy().reshape(2,2,2,2)
            rAB_pt = rAB_pt.transpose(0,3,2,1).reshape(4,4)
            pt_eigvals = np.real(np.linalg.eigvalsh(rAB_pt))
            negativity = abs(sum(e for e in pt_eigvals if e < 0))
            log_neg = math.log2(2*negativity + 1) if negativity > 0 else 0
            
            data.append({
                't': t, 'C': C, 'psi': psi_n, 'cpsi': cpsi, 'zstar': zstar,
                'purity': purity, 'lin_ent': lin_ent, 'S_vN': S_vN,
                'lam1': lam1, 'lam2': lam2, 'lam3': lam3, 'lam4': lam4,
                'lam_max': lam_max, 'lam_gap': lam_gap, 'lam_spread': lam_spread,
                'lam_prod12': lam_prod12, 'sqrt_cpsi': sqrt_cpsi,
                'pur_shifted': pur_shifted, 'pur_norm': pur_norm,
                'part_ratio': part_ratio, 'negativity': negativity,
                'log_neg': log_neg,
            })
        
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    return data

# ============================================================
# MAIN ANALYSIS
# ============================================================
if __name__ == "__main__":
    print("=" * 75)
    print("WHAT IS z* IN THE DENSITY MATRIX?")
    print("=" * 75)
    
    data = run_trajectory()
    
    # Filter to points where z* > 0 (CPsi > 0)
    valid = [d for d in data if d['zstar'] > 0.001]
    print(f"\n  Total points: {len(data)}, with z* > 0: {len(valid)}")
    
    if len(valid) < 3:
        print("  Not enough valid points!")
        sys.exit(1)
    
    z_arr = np.array([d['zstar'] for d in valid])
    
    # Test every candidate
    candidates = {
        'purity': [d['purity'] for d in valid],
        'lin_ent': [d['lin_ent'] for d in valid],
        'S_vN': [d['S_vN'] for d in valid],
        'lam1': [d['lam1'] for d in valid],
        'lam2': [d['lam2'] for d in valid],
        'lam_gap': [d['lam_gap'] for d in valid],
        'lam_spread': [d['lam_spread'] for d in valid],
        'lam_prod12': [d['lam_prod12'] for d in valid],
        'sqrt(CPsi)': [d['sqrt_cpsi'] for d in valid],
        'pur_shifted': [d['pur_shifted'] for d in valid],
        'pur_norm': [d['pur_norm'] for d in valid],
        'part_ratio': [d['part_ratio'] for d in valid],
        'negativity': [d['negativity'] for d in valid],
        'log_neg': [d['log_neg'] for d in valid],
        'C (concurrence)': [d['C'] for d in valid],
        'Psi (coherence)': [d['psi'] for d in valid],
    }
    
    # Linear correlations
    print(f"\n--- Linear correlation with z* ---")
    for name, vals in sorted(candidates.items(), key=lambda x: -abs(np.corrcoef(z_arr, x[1])[0,1])):
        r = np.corrcoef(z_arr, vals)[0,1]
        print(f"  {name:>20}: r = {r:+.6f}")
    
    # Now try NONLINEAR candidates: functions of eigenvalues
    print(f"\n--- Nonlinear candidates ---")
    nonlin = {}
    for d in valid:
        l1, l2, l3, l4 = d['lam1'], d['lam2'], d['lam3'], d['lam4']
        p = d['purity']
        c = d['C']
        psi = d['psi']
        z = d['zstar']
        cp = d['cpsi']
        
        # Try many functional forms
        nonlin.setdefault('sqrt(purity)', []).append(math.sqrt(p))
        nonlin.setdefault('(lam1-0.25)/0.75', []).append((l1-0.25)/0.75 if l1 > 0.25 else 0)
        nonlin.setdefault('sqrt(lam1-0.25)', []).append(math.sqrt(max(0,l1-0.25)))
        nonlin.setdefault('lam1-lam4', []).append(l1-l4)
        nonlin.setdefault('(lam1-lam4)/2', []).append((l1-l4)/2)
        nonlin.setdefault('sqrt(lam1*lam2)', []).append(math.sqrt(max(0,l1*l2)))
        nonlin.setdefault('C/2', []).append(c/2)
        nonlin.setdefault('C*psi/(C+psi)', []).append(c*psi/(c+psi) if (c+psi)>0.001 else 0)
        nonlin.setdefault('negativity', []).append(d['negativity'])
        nonlin.setdefault('negativity/2', []).append(d['negativity']/2)
        nonlin.setdefault('sqrt(negativity)', []).append(math.sqrt(d['negativity']) if d['negativity']>0 else 0)
        
        # Key test: is z* = negativity? 
        # Negativity N = (||rho^TB||_1 - 1)/2
        # For pure 2-qubit states, negativity = C/2
        
        # Another key test: z* related to conditional entropy?
        # z* related to Holevo quantity?
        
        # Try: z* = f(C, Psi) directly
        nonlin.setdefault('C*Psi', []).append(c*psi)
        nonlin.setdefault('(C+Psi)/4', []).append((c+psi)/4)
        nonlin.setdefault('sqrt(C*Psi)', []).append(math.sqrt(max(0,c*psi)))
    
    for name, vals in sorted(nonlin.items(), key=lambda x: -abs(np.corrcoef(z_arr, x[1])[0,1])):
        r = np.corrcoef(z_arr, vals)[0,1]
        # Also compute max absolute error for near-exact matches
        err = max(abs(z_arr[i] - vals[i]) for i in range(len(z_arr)))
        print(f"  {name:>25}: r = {r:+.6f}, max_err = {err:.6f}")
    
    # DIRECT COMPARISON TABLE for best candidates
    print(f"\n--- Point-by-point comparison (best candidates vs z*) ---")
    print(f"  {'t':>5} | {'z*':>7} | {'neg':>7} | {'C/2':>7} | {'lam_gap':>7} | {'sqrt_cp':>7} | {'(l1-l4)/2':>7}")
    print("-" * 70)
    for d in valid:
        z = d['zstar']
        print(f"  {d['t']:>5.2f} | {z:>7.4f} | {d['negativity']:>7.4f} | "
              f"{d['C']/2:>7.4f} | {d['lam_gap']:>7.4f} | "
              f"{d['sqrt_cpsi']:>7.4f} | {(d['lam1']-d['lam4'])/2:>7.4f}")
    
    # KEY TEST: Is z*(1-z*) EXACTLY equal to C*Psi_norm?
    print(f"\n--- Algebraic identity check: z*(1-z*) = CPsi ---")
    max_err = 0
    for d in valid:
        z = d['zstar']
        lhs = z * (1-z)
        rhs = d['cpsi']
        err = abs(lhs - rhs)
        if err > max_err:
            max_err = err
    print(f"  max |z*(1-z*) - CPsi| = {max_err:.2e}")
    print(f"  Identity holds: {'YES (exact)' if max_err < 1e-10 else 'APPROXIMATE'}")
    
    # MULTI-PARAMETER TEST: does the relationship change with J, gamma?
    print(f"\n{'='*75}")
    print("MULTI-PARAMETER VERIFICATION")
    print("Does z* = f(density matrix) hold universally?")
    print("="*75)
    
    configs = [
        (1.0, 1.0, 0.05, "symmetric J"),
        (1.0, 2.0, 0.05, "asymmetric J"),
        (1.0, 2.0, 0.20, "strong noise"),
        (1.0, 5.0, 0.02, "strong asymmetry"),
        (0.5, 0.5, 0.01, "weak coupling"),
    ]
    
    print(f"\n  {'config':>20} | {'r(z*,neg)':>10} | {'r(z*,C/2)':>10} | {'r(z*,pur)':>10} | {'r(z*,lam1)':>10}")
    print("-" * 75)
    
    for J_SA, J_SB, gamma, label in configs:
        d2 = run_trajectory(J_SA, J_SB, gamma, dt=0.01, t_max=5.0)
        v2 = [d for d in d2 if d['zstar'] > 0.001]
        if len(v2) < 3:
            print(f"  {label:>20} | too few points with z*>0")
            continue
        z2 = np.array([d['zstar'] for d in v2])
        r_neg = np.corrcoef(z2, [d['negativity'] for d in v2])[0,1]
        r_c2 = np.corrcoef(z2, [d['C']/2 for d in v2])[0,1]
        r_pur = np.corrcoef(z2, [d['purity'] for d in v2])[0,1]
        r_l1 = np.corrcoef(z2, [d['lam1'] for d in v2])[0,1]
        print(f"  {label:>20} | {r_neg:>10.4f} | {r_c2:>10.4f} | {r_pur:>10.4f} | {r_l1:>10.4f}")
    
    # FINAL: What is the BEST identity?
    print(f"\n{'='*75}")
    print("CONCLUSION")
    print("="*75)
    print("""
  z* is defined by z*(1-z*) = CPsi (algebraically exact).
  
  The question is whether z* equals some KNOWN quantum quantity.
  
  If z* = negativity (or C/2 for pure states), then CPsi = N(1-N)
  and the 1/4 boundary means: negativity cannot exceed 1/2.
  
  If z* = f(eigenvalues), then CPsi has a spectral interpretation.
  
  If z* does NOT match any standard quantity, then it is genuinely
  new: a composite diagnostic that combines entanglement (C) and
  coherence (Psi) into a single number with Bernoulli structure.
""")
