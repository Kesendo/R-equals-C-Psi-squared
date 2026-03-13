"""
What IS z* in the density matrix?
==================================
z* = (1 - sqrt(1 - 4*CPsi)) / 2
z*(1-z*) = CPsi
Correlation with purity: r = 0.917 (close but not exact)

Strategy: compute z* alongside ALL density matrix properties
and look for exact matches.
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

def run_and_analyze():
    J_SA, J_SB, gamma = 1.0, 2.0, 0.05
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
    rho = np.outer(psi0, psi0.conj())
    YY = np.kron(Y, Y)
    dt = 0.02
    t_max = 5.0  # short, we only need early dynamics where CPsi > 0
    
    data = []
    for step in range(int(t_max/dt)+1):
        t = step * dt
        if step % 2 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], n_q)
            
            # Eigenvalues of rho_AB (sorted descending)
            eigvals = np.sort(np.real(np.linalg.eigvalsh(rAB)))[::-1]
            lam1, lam2, lam3, lam4 = eigvals
            
            # Concurrence
            rho_t = YY @ rAB.conj() @ YY
            R_mat = rAB @ rho_t
            ev_r = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R_mat), 0))))[::-1]
            C = max(0, float(ev_r[0] - ev_r[1] - ev_r[2] - ev_r[3]))
            
            # l1-coherence normalized
            l1 = sum(abs(rAB[i,j]) for i in range(4) for j in range(4) if i!=j)
            psi_n = float(l1) / 3.0
            
            cpsi = C * psi_n
            
            # z*
            if cpsi > 0.001:
                z_star = (1 - math.sqrt(max(0, 1-4*cpsi))) / 2
            else:
                z_star = 0
            
            # Purity and linear entropy
            purity = float(np.real(np.trace(rAB @ rAB)))
            lin_ent = (4/3)*(1 - purity)
            
            # Von Neumann entropy
            S_vN = -sum(e*math.log2(e) for e in eigvals if e > 1e-12)
            
            # Candidate expressions for z*
            candidates = {}
            candidates['lam1'] = lam1
            candidates['lam1-lam2'] = lam1 - lam2
            candidates['lam1-lam4'] = lam1 - lam4
            candidates['lam1*lam2'] = lam1 * lam2
            candidates['sqrt(lam1)'] = math.sqrt(max(0, lam1))
            candidates['sqrt(lam1*lam2)'] = math.sqrt(max(0, lam1*lam2))
            candidates['(lam1-0.25)*4'] = (lam1 - 0.25) * 4
            candidates['1-2*lam4'] = 1 - 2*lam4
            candidates['2*lam1-0.5'] = 2*lam1 - 0.5
            candidates['lam1-lam2-lam3+lam4'] = lam1-lam2-lam3+lam4
            candidates['(lam1-lam3)'] = lam1 - lam3
            candidates['(lam2-lam4)'] = lam2 - lam4
            candidates['purity'] = purity
            candidates['2*purity-0.5'] = 2*purity - 0.5
            candidates['sqrt(purity-0.25)'] = math.sqrt(max(0, purity-0.25))
            candidates['(purity-0.25)*4'] = (purity - 0.25) * 4
            
            # Concurrence-related
            candidates['C/2'] = C / 2
            candidates['C*C'] = C * C
            candidates['sqrt(C)'] = math.sqrt(max(0, C))
            
            # R-matrix eigenvalues (sqrt of eigenvalues of rho*rho_tilde)
            r1, r2, r3, r4 = ev_r
            candidates['r1'] = float(r1)
            candidates['r1/2'] = float(r1) / 2
            candidates['r1*r1'] = float(r1*r1)
            candidates['(r1-r2)/2'] = float(r1-r2) / 2
            
            # Partial transpose eigenvalues (for negativity)
            rAB_pt = rAB.copy().reshape(2,2,2,2).transpose(0,3,2,1).reshape(4,4)
            pt_eigvals = np.sort(np.real(np.linalg.eigvalsh(rAB_pt)))
            negativity = -sum(e for e in pt_eigvals if e < 0)
            candidates['negativity'] = negativity
            candidates['negativity/2'] = negativity / 2
            candidates['2*negativity'] = 2 * negativity
            
            if cpsi > 0.001:
                data.append({
                    't': t, 'z_star': z_star, 'cpsi': cpsi,
                    'C': C, 'psi': psi_n,
                    'lam': (lam1, lam2, lam3, lam4),
                    'candidates': candidates
                })
        
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    return data

# ============================================================
# MAIN
# ============================================================
print("=" * 75)
print("WHAT IS z* IN THE DENSITY MATRIX?")
print("=" * 75)

data = run_and_analyze()
print(f"\nCollected {len(data)} data points where CPsi > 0.001")

# For each candidate, compute correlation with z*
print(f"\n--- Correlation of each candidate with z* ---")
print(f"  r = 1.0000 means PERFECT match")
print(f"{'candidate':>30} | {'r':>8} | {'max_err':>8} | {'mean_err':>8}")
print("-" * 65)

z_vals = [d['z_star'] for d in data]
cand_names = list(data[0]['candidates'].keys())

results = []
for name in cand_names:
    c_vals = [d['candidates'][name] for d in data]
    if len(set(c_vals)) < 2:
        continue
    r = np.corrcoef(z_vals, c_vals)[0,1]
    errors = [abs(z - c) for z, c in zip(z_vals, c_vals)]
    max_err = max(errors)
    mean_err = sum(errors) / len(errors)
    results.append((name, r, max_err, mean_err))

# Sort by absolute correlation
results.sort(key=lambda x: -abs(x[1]))
for name, r, mx, mn in results:
    marker = " <== EXACT!" if mx < 0.001 else (" <-- CLOSE" if mx < 0.01 else "")
    print(f"{name:>30} | {r:>8.4f} | {mx:>8.5f} | {mn:>8.5f}{marker}")

# Show the best candidates in detail
print(f"\n--- TOP 3 CANDIDATES (detail) ---")
for name, r, mx, mn in results[:3]:
    print(f"\n  {name}: r={r:.6f}, max_error={mx:.6f}")
    c_vals = [d['candidates'][name] for d in data]
    for d, c in zip(data[:6], c_vals[:6]):
        print(f"    t={d['t']:.2f}: z*={d['z_star']:.6f}, {name}={c:.6f}, "
              f"diff={abs(d['z_star']-c):.6f}")

# BRUTE FORCE: try linear combinations a*X + b for top candidates
print(f"\n--- LINEAR FIT: z* = a*X + b ---")
for name, r, mx, mn in results[:5]:
    c_vals = np.array([d['candidates'][name] for d in data])
    z_arr = np.array(z_vals)
    # Least squares: z* = a * candidate + b
    A = np.vstack([c_vals, np.ones(len(c_vals))]).T
    coef, residuals, _, _ = np.linalg.lstsq(A, z_arr, rcond=None)
    a, b = coef
    fitted = a * c_vals + b
    fit_err = np.max(np.abs(z_arr - fitted))
    print(f"  z* = {a:.4f} * {name} + {b:.4f}  (max_err={fit_err:.6f})")

# KEY TEST: Is z* exactly some function of the eigenvalues?
print(f"\n--- EXACT EIGENVALUE RELATIONSHIPS ---")
print(f"  Checking if z*(1-z*) = f(lam1, lam2, lam3, lam4)")
print(f"")
for d in data[:8]:
    z = d['z_star']
    l1, l2, l3, l4 = d['lam']
    cpsi = d['cpsi']
    
    # z*(1-z*) = CPsi. What is CPsi in terms of eigenvalues?
    # Try various eigenvalue combinations
    tests = {
        'l1*l2 - l3*l4': l1*l2 - l3*l4,
        'l1*l4 - l2*l3': l1*l4 - l2*l3,
        '(l1-l4)*(l2-l3)': (l1-l4)*(l2-l3),
        '(l1-l3)*(l2-l4)': (l1-l3)*(l2-l4),
        '(l1+l4)*(l2+l3)-0.25': (l1+l4)*(l2+l3)-0.25,
        'l1*l2+l3*l4': l1*l2+l3*l4,
        'purity-sum(li^2)': sum(li**2 for li in [l1,l2,l3,l4]),
    }
    
    print(f"  t={d['t']:.2f}: CPsi={cpsi:.6f}, eigenvals=[{l1:.4f},{l2:.4f},{l3:.4f},{l4:.4f}]")
    for name, val in tests.items():
        match = "MATCH!" if abs(val - cpsi) < 0.0001 else ""
        if d == data[0]:
            print(f"    {name:>30} = {val:.6f} vs CPsi={cpsi:.6f} {match}")

# Check specifically: is CPsi a function of purity alone?
print(f"\n--- CPsi vs simple functions of purity ---")
print(f"  {'t':>5} | {'CPsi':>8} | {'purity':>8} | {'pur-0.25':>8} | {'4*(p-.25)':>8} | ratio")
print("-" * 65)
for d in data[:10]:
    p = sum(li**2 for li in d['lam'])
    ratio = d['cpsi'] / (p - 0.25) if p > 0.251 else 0
    print(f"  {d['t']:>5.2f} | {d['cpsi']:>8.5f} | {p:>8.5f} | {p-0.25:>8.5f} | "
          f"{4*(p-0.25):>8.5f} | {ratio:>8.4f}")

print(f"\n{'='*75}")
print("DONE")
print("="*75)
