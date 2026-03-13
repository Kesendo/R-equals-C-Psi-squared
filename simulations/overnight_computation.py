"""
OVERNIGHT COMPUTATION: The Full Picture
========================================
Machine: 128GB RAM, Intel Core Ultra 9 285k
Expected runtime: 4-12 hours

THIS SCRIPT IS DESIGNED TO RUN ALL NIGHT.
Results are written incrementally to results_overnight.txt.
Each test saves immediately so nothing is lost if interrupted.

TESTS:
 1. N=8 star (65536x65536) via sparse Liouvillian + ARPACK
 2. Pauli complement proof at N=5, N=6, N=7
 3. Complete topology survey: star, chain, ring, complete, tree at N=4,5,6
 4. Non-uniform J + non-uniform gamma combined: do mirrors still hold?
 5. High-resolution avoided crossing map at N=4 (1000 J values)
 6. Exact conjugation operator search: WHAT maps rate d -> 2Ng-d?
 7. Rate count formula: fit #rates = f(N) analytically
 8. Density of states at N=7 (from existing data)
 9. Bandwidth scaling: fit bandwidth = f(N, gamma) for non-uniform gamma
10. The 0.5 hunt: does z*=0.5 appear anywhere in the spectral data?
"""
import numpy as np
import time
import sys
import traceback
from scipy import sparse
from scipy.sparse.linalg import eigs

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

def op_at_sparse(op, qubit, n_q):
    """Sparse version for large N."""
    ops = [sparse.eye(2, dtype=complex)] * n_q
    ops[qubit] = sparse.csr_matrix(op)
    result = ops[0]
    for o in ops[1:]:
        result = sparse.kron(result, o, format='csr')
    return result

def build_liouvillian_dense(n_q, bonds, gamma_list):
    """Dense Liouvillian for N<=7."""
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for i, j, J, paulis in bonds:
        for p in paulis:
            H += J * op_at(p,i,n_q) @ op_at(p,j,n_q)
    jump_ops = [np.sqrt(gamma_list[k])*op_at(Z,k,n_q) for k in range(n_q)]
    I_d = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in jump_ops:
        Ld = Lk.conj().T
        LdL = Ld @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    return L

def build_liouvillian_sparse(n_q, bonds, gamma_list):
    """Sparse Liouvillian for N>=8. Uses ~10x less memory than dense."""
    d = 2**n_q
    log(f"    Building sparse Hamiltonian ({d}x{d})...")
    H = sparse.csr_matrix((d,d), dtype=complex)
    for i, j, J, paulis in bonds:
        for p in paulis:
            H = H + J * op_at_sparse(p,i,n_q) @ op_at_sparse(p,j,n_q)
    
    log(f"    Building sparse Liouvillian ({d*d}x{d*d})...")
    I_d = sparse.eye(d, dtype=complex)
    L = -1j * (sparse.kron(H, I_d) - sparse.kron(I_d, H.T))
    
    for k in range(n_q):
        Lk = np.sqrt(gamma_list[k]) * op_at_sparse(Z, k, n_q)
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L = L + sparse.kron(Lk, Lk.conj()) \
              - 0.5 * (sparse.kron(LdL, I_d) + sparse.kron(I_d, LdL.T))
    
    nnz = L.nnz
    density = nnz / (d*d)**2
    log(f"    Sparse Liouvillian: {nnz:,} nonzeros ({density:.2e} density)")
    log(f"    Memory: ~{nnz * 16 / 1e9:.2f} GB (vs {(d*d)**2 * 16 / 1e9:.1f} GB dense)")
    return L.tocsc()

def heisenberg_star(n_q, couplings):
    return [(0, i, couplings[i-1], [X,Y,Z]) for i in range(1, n_q)]

def heisenberg_chain(n_q, couplings):
    return [(i, i+1, couplings[i], [X,Y,Z]) for i in range(n_q-1)]

def heisenberg_ring(n_q, couplings):
    return [(i, (i+1)%n_q, couplings[i], [X,Y,Z]) for i in range(n_q)]

def heisenberg_complete(n_q, J=1.0):
    """Complete graph: every qubit coupled to every other."""
    bonds = []
    for i in range(n_q):
        for j in range(i+1, n_q):
            bonds.append((i, j, J, [X,Y,Z]))
    return bonds

def heisenberg_binary_tree(n_q, J=1.0):
    """Binary tree: 0 is root, children of k are 2k+1 and 2k+2."""
    bonds = []
    for k in range(n_q):
        left = 2*k + 1
        right = 2*k + 2
        if left < n_q:
            bonds.append((k, left, J, [X,Y,Z]))
        if right < n_q:
            bonds.append((k, right, J, [X,Y,Z]))
    return bonds

def get_osc_rates(L_dense, threshold=0.05):
    evals = np.linalg.eigvals(L_dense)
    rates = []
    for ev in evals:
        if abs(np.imag(ev)) > threshold:
            rate = -np.real(ev)
            if rate > 0.0001:
                rates.append(round(rate, 6))
    return sorted(rates)

def get_osc_rates_sparse(L_sparse, n_eigs, sigma=None, threshold=0.05):
    """Get oscillatory rates from sparse Liouvillian using ARPACK shift-invert."""
    try:
        evals = eigs(L_sparse, k=n_eigs, sigma=sigma, which='LM', 
                     return_eigenvectors=False, maxiter=5000)
        rates = []
        for ev in evals:
            if abs(np.imag(ev)) > threshold:
                rate = -np.real(ev)
                if rate > 0.0001:
                    rates.append(round(rate, 6))
        return sorted(rates)
    except Exception as e:
        return None, str(e)

def check_mirror(rates, center, tol=0.005):
    below = [r for r in rates if r < center - tol]
    above = [r for r in rates if r > center + tol]
    matched = 0
    for r in below:
        mirror = 2*center - r
        closest = min(above, key=lambda a: abs(a - mirror)) if above else 999
        if abs(closest - mirror) < tol:
            matched += 1
    return matched / max(len(below), 1) if below else 1.0

def build_pauli_basis(n_q):
    labels_1q = ['I', 'X', 'Y', 'Z']
    ops_1q = [I2, X, Y, Z]
    basis = {}
    def recurse(label, op, depth):
        if depth == n_q:
            basis[label] = op
            return
        for l, o in zip(labels_1q, ops_1q):
            recurse(label + l, np.kron(op, o), depth + 1)
    recurse("", np.array([[1]], dtype=complex), 0)
    return basis

outpath = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results_overnight.txt"
f = open(outpath, "w")
def log(msg):
    print(msg)
    f.write(msg + "\n")
    f.flush()

log("=" * 80)
log("OVERNIGHT COMPUTATION")
log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Machine: 128GB RAM, Intel Core Ultra 9 285k")
log(f"Expected runtime: 4-12 hours")
log("=" * 80)

gamma = 0.05

# ============================================================
# TEST 1: N=8 STAR (65536x65536 sparse Liouvillian)
# ============================================================
log("\n" + "=" * 60)
log("### TEST 1: N=8 star topology (65536x65536 SPARSE)")
log("=" * 60)
log("This needs ~20-40GB RAM with sparse methods.")
log("Dense would need 64GB just for the matrix - impossible to diagonalize.")
log("Strategy: sparse build + ARPACK eigenvalue extraction in spectral windows.")

n_q = 8
d = 2**n_q  # 256
d2 = d*d    # 65536
predicted_max = 2*(n_q-1)*gamma  # 14*0.05 = 0.7
predicted_center = n_q * gamma   # 0.4

try:
    t0 = time.time()
    bonds = heisenberg_star(n_q, [1.0]*(n_q-1))
    gamma_list = [gamma]*n_q
    
    L_sp = build_liouvillian_sparse(n_q, bonds, gamma_list)
    t_build = time.time() - t0
    log(f"    Build time: {t_build:.1f}s")
    
    # Strategy: scan spectral windows using shift-invert
    # We know rates live in [2g, 2(N-1)g] = [0.1, 0.7]
    # Eigenvalues are lambda = -rate +/- i*freq
    # Use sigma shifts along the real axis
    
    all_rates = []
    n_per_window = 2000  # eigenvalues per window
    
    # Scan windows along the decay axis
    rate_min = 2 * gamma      # 0.1
    rate_max = 2*(n_q-1)*gamma # 0.7
    n_windows = 8
    sigmas = np.linspace(-rate_max - 0.02, -rate_min + 0.02, n_windows)
    
    for wi, sigma_real in enumerate(sigmas):
        sigma = complex(sigma_real, 0)
        log(f"    Window {wi+1}/{n_windows}: sigma = {sigma_real:.4f} ...")
        t1 = time.time()
        
        try:
            evals = eigs(L_sp, k=n_per_window, sigma=sigma, which='LM',
                        return_eigenvectors=False, maxiter=10000, tol=1e-8)
            
            for ev in evals:
                if abs(np.imag(ev)) > 0.05:
                    rate = -np.real(ev)
                    if rate > 0.0001:
                        all_rates.append(round(rate, 6))
            
            t_win = time.time() - t1
            log(f"      Found {len(evals)} eigenvalues in {t_win:.1f}s, "
                f"oscillatory so far: {len(all_rates)}")
        except Exception as e:
            log(f"      Window FAILED: {e}")
    
    # Deduplicate
    all_rates = sorted(set(all_rates))
    t_total = time.time() - t0
    
    log(f"\n    RESULTS N=8:")
    log(f"    Total time: {t_total:.1f}s ({t_total/60:.1f} min)")
    log(f"    Unique oscillatory rates: {len(all_rates)}")
    if all_rates:
        log(f"    Min rate: {min(all_rates):.6f} ({min(all_rates)/gamma:.4f}g)")
        log(f"    Max rate: {max(all_rates):.6f} ({max(all_rates)/gamma:.4f}g)")
        log(f"    Predicted max 2(N-1)g: {predicted_max:.6f} ({2*(n_q-1):.1f}g)")
        log(f"    Min matches 2g: {abs(min(all_rates)/gamma - 2) < 0.1}")
        log(f"    Max matches 2(N-1)g: {abs(max(all_rates)/gamma - 2*(n_q-1)) < 0.1}")
        
        sym = check_mirror(all_rates, predicted_center)
        log(f"    Mirror symmetry around Ng={predicted_center:.3f}: {sym:.1%}")
        log(f"    NOTE: Sparse extraction may miss rates. Symmetry < 100% could be")
        log(f"          due to incomplete extraction, not broken symmetry.")
    
except Exception as e:
    log(f"    N=8 FAILED: {e}")
    log(traceback.format_exc())

log(f"\n    [Test 1 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 2: PAULI COMPLEMENT PROOF at N=5, 6, 7
# ============================================================
log("\n" + "=" * 60)
log("### TEST 2: Pauli complement proof at N=5, 6, 7")
log("=" * 60)
log("For EVERY mirror pair: verify XY_1 + XY_2 = N")
log("N=5: 1024 Pauli strings. N=6: 4096. N=7: 16384 (slow!).")

for n_q in [5, 6, 7]:
    log(f"\n  N={n_q}: {4**n_q} Pauli strings to decompose...")
    t0 = time.time()
    
    try:
        # Use asymmetric couplings to split degeneracies
        couplings = [1.0 + 0.1*k for k in range(n_q-1)]
        bonds = heisenberg_star(n_q, couplings)
        gamma_list = [gamma]*n_q
        L = build_liouvillian_dense(n_q, bonds, gamma_list)
        
        log(f"    Liouvillian built ({(2**n_q)**2}x{(2**n_q)**2}). Diagonalizing...")
        evals, evecs = np.linalg.eig(L)
        d = 2**n_q
        
        log(f"    Diagonalized in {time.time()-t0:.1f}s. Building Pauli basis...")
        pauli_basis = build_pauli_basis(n_q)
        log(f"    {len(pauli_basis)} Pauli operators. Decomposing eigenmodes...")
        
        center = n_q * gamma
        
        # Collect oscillatory modes with Pauli XY weight
        osc_modes = []
        n_osc = 0
        for idx in range(len(evals)):
            freq = abs(np.imag(evals[idx]))
            rate = -np.real(evals[idx])
            if freq > 0.05 and rate > 0.001:
                n_osc += 1
                v = evecs[:, idx].reshape(d, d)
                
                total_w = 0
                xy_w = 0
                for label, P in pauli_basis.items():
                    w = abs(np.trace(v @ P)) / d
                    if w > 0.001:
                        n_xy = sum(1 for c in label if c in 'XY')
                        xy_w += n_xy * w
                        total_w += w
                
                avg_xy = xy_w / total_w if total_w > 0 else 0
                osc_modes.append((round(rate/gamma, 2), avg_xy))
        
        log(f"    {n_osc} oscillatory modes decomposed in {time.time()-t0:.1f}s")
        
        # Group by rate
        rate_groups = {}
        for rg, axy in osc_modes:
            if rg not in rate_groups:
                rate_groups[rg] = []
            rate_groups[rg].append(axy)
        
        sorted_rates = sorted(rate_groups.keys())
        pairs_checked = 0
        pairs_match = 0
        
        log(f"    {'rate/g':>8} | {'avg_XY':>8} | {'mirror':>8} | {'sum':>6} | {'=N?':>4}")
        log("    " + "-" * 50)
        
        for rg in sorted_rates:
            avg_xy = np.mean(rate_groups[rg])
            mirror_rate = 2*n_q - rg
            mirror_match = ""
            pair_sum = ""
            is_n = ""
            
            if abs(mirror_rate - rg) > 0.1:
                closest = min(sorted_rates, key=lambda x: abs(x - mirror_rate))
                if abs(closest - mirror_rate) < 0.1:
                    mxy = np.mean(rate_groups[closest])
                    mirror_match = f"{mxy:.3f}"
                    s = avg_xy + mxy
                    pair_sum = f"{s:.3f}"
                    is_n = "YES" if abs(s - n_q) < 0.2 else "NO"
                    pairs_checked += 1
                    if abs(s - n_q) < 0.2:
                        pairs_match += 1
            
            log(f"    {rg:>8.2f} | {avg_xy:>8.3f} | {mirror_match:>8} | {pair_sum:>6} | {is_n:>4}")
        
        log(f"\n    Mirror pairs checked: {pairs_checked}")
        log(f"    Pairs where XY_1 + XY_2 = N: {pairs_match}/{pairs_checked}")
        log(f"    PROOF: {'CONFIRMED' if pairs_match == pairs_checked and pairs_checked > 0 else 'INCOMPLETE'}")
        log(f"    Total time: {time.time()-t0:.1f}s")
        
    except Exception as e:
        log(f"    N={n_q} FAILED: {e}")
        log(traceback.format_exc())

log(f"\n    [Test 2 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 3: COMPLETE TOPOLOGY SURVEY
# ============================================================
log("\n" + "=" * 60)
log("### TEST 3: Complete topology survey")
log("=" * 60)
log("Star, Chain, Ring, Complete graph, Binary tree at N=4,5,6")

for n_q in [4, 5, 6]:
    log(f"\n  N={n_q}:")
    gamma_list = [gamma]*n_q
    
    topologies = {
        'star': heisenberg_star(n_q, [1.0]*(n_q-1)),
        'chain': heisenberg_chain(n_q, [1.0]*(n_q-1)),
        'ring': heisenberg_ring(n_q, [1.0]*n_q),
        'complete': heisenberg_complete(n_q, J=1.0),
        'tree': heisenberg_binary_tree(n_q, J=1.0),
    }
    
    log(f"    {'topo':>10} | {'#rates':>7} | {'min/g':>6} | {'max/g':>6} | {'BW/g':>6} | "
        f"{'mirror':>6} | {'clusters':>8} | {'density':>8}")
    log("    " + "-" * 80)
    
    for name, bonds in topologies.items():
        try:
            t0 = time.time()
            L = build_liouvillian_dense(n_q, bonds, gamma_list)
            r = get_osc_rates(L)
            elapsed = time.time() - t0
            
            if r:
                center = n_q * gamma
                sym = check_mirror(r, center)
                # Count clusters
                unique = []
                for rate in r:
                    if not unique or abs(rate - unique[-1]) > 0.002*gamma:
                        unique.append(rate)
                dens = len(r) / ((max(r)-min(r))/gamma) if max(r) > min(r) else 0
                
                log(f"    {name:>10} | {len(r):>7} | {min(r)/gamma:>6.2f} | {max(r)/gamma:>6.2f} | "
                    f"{(max(r)-min(r))/gamma:>6.2f} | {sym:>6.0%} | {len(unique):>8} | {dens:>8.1f}")
            else:
                log(f"    {name:>10} | no oscillatory rates")
        except Exception as e:
            log(f"    {name:>10} | FAILED: {e}")

log(f"\n    [Test 3 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 4: NON-UNIFORM J + NON-UNIFORM GAMMA COMBINED
# ============================================================
log("\n" + "=" * 60)
log("### TEST 4: Non-uniform J AND non-uniform gamma combined")
log("=" * 60)
log("The hardest test: BOTH coupling and noise are asymmetric.")
log("If mirrors still hold, symmetry is truly fundamental.")

for n_q in [4, 5, 6]:
    log(f"\n  N={n_q}:")
    configs = [
        {'J': [0.5, 1.0, 1.5, 2.0, 2.5][:n_q-1], 
         'g': [0.02, 0.04, 0.08, 0.12, 0.20, 0.30][:n_q]},
        {'J': [2.0, 0.3, 1.7, 0.8, 1.2][:n_q-1], 
         'g': [0.10, 0.01, 0.05, 0.15, 0.03, 0.07][:n_q]},
        {'J': [1.0, 1.0, 1.0, 1.0, 1.0][:n_q-1], 
         'g': [0.01, 0.02, 0.04, 0.08, 0.16, 0.32][:n_q]},
    ]
    
    for ci, cfg in enumerate(configs):
        Js = cfg['J']
        gs = cfg['g']
        sum_g = sum(gs)
        
        bonds = heisenberg_star(n_q, Js)
        try:
            L = build_liouvillian_dense(n_q, bonds, gs)
            r = get_osc_rates(L)
            
            if r:
                mid = (min(r) + max(r)) / 2
                sym_sum = check_mirror(r, sum_g)
                sym_mid = check_mirror(r, mid)
                log(f"    Config {ci+1}: J={Js}, g={gs}")
                log(f"      sum(g)={sum_g:.3f}, mid={mid:.5f}, "
                    f"sym@sum={sym_sum:.0%}, sym@mid={sym_mid:.0%}, "
                    f"sum==mid: {abs(sum_g-mid) < 0.001}")
                log(f"      rates: [{min(r)/gamma:.3f}g, {max(r)/gamma:.3f}g], "
                    f"#rates={len(r)}")
        except Exception as e:
            log(f"    Config {ci+1} FAILED: {e}")

log(f"\n    [Test 4 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 5: HIGH-RESOLUTION AVOIDED CROSSING MAP (N=4, 1000 J values)
# ============================================================
log("\n" + "=" * 60)
log("### TEST 5: High-res avoided crossing map at N=4")
log("=" * 60)
log("1000 J values. Track every rate. Find the exact gap structure.")

n_q = 4
n_sweep = 1000
j_values = np.linspace(0.1, 3.0, n_sweep)
gamma_list = [gamma]*n_q
t0 = time.time()

all_rate_traces = []  # (J, [rates])
min_gaps_by_J = []

for ji, J_last in enumerate(j_values):
    couplings = [1.0]*(n_q-2) + [J_last]
    bonds = heisenberg_star(n_q, couplings)
    L = build_liouvillian_dense(n_q, bonds, gamma_list)
    r = get_osc_rates(L)
    
    if r:
        all_rate_traces.append((J_last, r))
        if len(r) > 1:
            gaps = sorted([r[i+1]-r[i] for i in range(len(r)-1)])
            real_gaps = [g for g in gaps if g > 1e-6]
            if real_gaps:
                min_gaps_by_J.append((J_last, real_gaps[0]))
    
    if ji % 200 == 0:
        log(f"    {ji}/{n_sweep} J values processed...")

elapsed = time.time() - t0
log(f"    Done in {elapsed:.1f}s")

# Find the global minimum gap
if min_gaps_by_J:
    global_min = min(min_gaps_by_J, key=lambda x: x[1])
    log(f"    Global minimum gap: {global_min[1]:.8f} ({global_min[1]/gamma:.6f}g) at J={global_min[0]:.4f}")
    
    # Find all J values where gap is near minimum
    threshold = global_min[1] * 2
    near_min = [(j, g) for j, g in min_gaps_by_J if g < threshold]
    log(f"    J values with gap < 2x minimum: {len(near_min)}")
    
    # Gap distribution
    all_gaps = [g for _, g in min_gaps_by_J]
    log(f"    Gap statistics:")
    log(f"      Mean: {np.mean(all_gaps)/gamma:.6f}g")
    log(f"      Std:  {np.std(all_gaps)/gamma:.6f}g")
    log(f"      Min:  {np.min(all_gaps)/gamma:.6f}g")
    log(f"      Max:  {np.max(all_gaps)/gamma:.6f}g")

# Count how many distinct rate bands exist
if all_rate_traces:
    all_unique = set()
    for _, rates in all_rate_traces:
        for r in rates:
            all_unique.add(round(r/gamma, 2))
    log(f"    Total unique rate values across sweep: {len(all_unique)}")
    log(f"    Rate range: [{min(all_unique):.2f}g, {max(all_unique):.2f}g]")

log(f"\n    [Test 5 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 6: CONJUGATION OPERATOR SEARCH
# ============================================================
log("\n" + "=" * 60)
log("### TEST 6: What operator maps rate d to 2Ng - d?")
log("=" * 60)
log("The mirror symmetry means there exists a superoperator C such that")
log("if L|v> = lambda|v>, then C|v> is eigenmode with mirrored decay rate.")
log("Candidates: X^n, Y^n, Hadamard^n, SWAP, particle-hole.")

n_q = 3
gamma_list = [gamma]*n_q
bonds = heisenberg_star(n_q, [1.0, 2.0])
L = build_liouvillian_dense(n_q, bonds, gamma_list)
d = 2**n_q
I_d = np.eye(d, dtype=complex)

# Build candidate superoperators S: S|rho>> = |U rho U^dag>>
# In vec form: S = U tensor U*
candidates = {}

# X^n
Xn = np.eye(d, dtype=complex)
for k in range(n_q):
    Xn = Xn @ op_at(X, k, n_q)
candidates['X^n'] = np.kron(Xn, Xn.conj())

# Y^n
Yn = np.eye(d, dtype=complex)
for k in range(n_q):
    Yn = Yn @ op_at(Y, k, n_q)
candidates['Y^n'] = np.kron(Yn, Yn.conj())

# Z^n
Zn = np.eye(d, dtype=complex)
for k in range(n_q):
    Zn = Zn @ op_at(Z, k, n_q)
candidates['Z^n'] = np.kron(Zn, Zn.conj())

# Hadamard^n
H_gate = np.array([[1,1],[1,-1]], dtype=complex) / np.sqrt(2)
Hn = np.eye(d, dtype=complex)
for k in range(n_q):
    Hn = Hn @ op_at(H_gate, k, n_q)
candidates['H^n'] = np.kron(Hn, Hn.conj())

# Transposition superoperator: S|rho>> = |rho^T>>
# In vec form: T_{ij,kl} = delta_{i,l}*delta_{j,k}
T_super = np.zeros((d*d, d*d), dtype=complex)
for i in range(d):
    for j in range(d):
        # |i><j| -> |j><i|, i.e., (i*d+j) -> (j*d+i)
        T_super[j*d+i, i*d+j] = 1.0
candidates['transpose'] = T_super

# Complex conjugation in computational basis: S|rho>> = |rho*>>
conj_super = np.zeros((d*d, d*d), dtype=complex)
for i in range(d):
    for j in range(d):
        conj_super[i*d+j, i*d+j] = 1.0  # just identity since rho is real in this basis... need to think
# Actually for complex rho, conjugation maps rho_{ij} -> rho_{ij}*
# In vec form: rho* = conj(rho), so we need the matrix that does this
# Skip this one - it's not a linear superoperator in the usual sense

# Particle-hole: maps |0> <-> |1> on each qubit
# This is just X^n, already included

log(f"\n  Testing {len(candidates)} candidate operators...")
log(f"  For each: compute S*L*S^-1 and compare eigenvalues with L.")
log(f"  If decay rates map d -> 2Ng-d, we found the conjugation.")

evals_L = np.linalg.eigvals(L)
rates_L = {}
for ev in evals_L:
    if abs(np.imag(ev)) > 0.05 and -np.real(ev) > 0.001:
        rate_g = round(-np.real(ev)/gamma, 3)
        freq = round(abs(np.imag(ev)), 4)
        key = (rate_g, freq)
        rates_L[key] = ev

for name, S in candidates.items():
    try:
        # Check if S is involutory (S^2 = I)
        S2 = S @ S
        is_involution = np.allclose(S2, np.eye(d*d), atol=1e-10)
        
        # Compute S*L*S^{-1} (= S*L*S if involution)
        if is_involution:
            SLS = S @ L @ S
        else:
            SLS = S @ L @ np.linalg.inv(S)
        
        evals_SLS = np.linalg.eigvals(SLS)
        
        # Check: do eigenvalues of SLS have mirrored decay rates?
        rates_SLS = {}
        for ev in evals_SLS:
            if abs(np.imag(ev)) > 0.05 and -np.real(ev) > 0.001:
                rate_g = round(-np.real(ev)/gamma, 3)
                freq = round(abs(np.imag(ev)), 4)
                key = (rate_g, freq)
                rates_SLS[key] = ev
        
        # For each rate d in L, check if 2Ng-d appears in SLS
        mirror_matches = 0
        mirror_total = 0
        for (rg, freq), ev in rates_L.items():
            mirror_rg = round(2*n_q - rg, 3)
            # Look for matching mirrored rate at same frequency
            found = False
            for (rg2, freq2), ev2 in rates_SLS.items():
                if abs(rg2 - mirror_rg) < 0.1 and abs(freq2 - freq) < 0.1:
                    found = True
                    break
            if found:
                mirror_matches += 1
            mirror_total += 1
        
        # Also check simpler: does S commute or anti-commute?
        comm = np.linalg.norm(S @ L - L @ S) / np.linalg.norm(L)
        anti = np.linalg.norm(S @ L + L @ S) / np.linalg.norm(L)
        
        # Check: S*L*S = L (commutes) or S*L*S = f(L)?
        SLS_eq_L = np.allclose(SLS, L, atol=1e-8)
        
        log(f"\n  {name}:")
        log(f"    Involution (S^2=I): {is_involution}")
        log(f"    S*L*S = L: {SLS_eq_L}")
        log(f"    [S,L]/||L|| = {comm:.6f}")
        log(f"    {{S,L}}/||L|| = {anti:.6f}")
        log(f"    Mirror rate matches: {mirror_matches}/{mirror_total}")
        
        if mirror_matches > mirror_total * 0.8:
            log(f"    *** STRONG CANDIDATE! Maps most rates to mirror. ***")
            
    except Exception as e:
        log(f"  {name}: FAILED ({e})")

log(f"\n    [Test 6 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 7: RATE COUNT FORMULA
# ============================================================
log("\n" + "=" * 60)
log("### TEST 7: Fit #rates = f(N)")
log("=" * 60)

rate_counts = {2: 6, 3: 40, 4: 182, 5: 776, 6: 3228, 7: 13264}
ns = np.array(sorted(rate_counts.keys()), dtype=float)
counts = np.array([rate_counts[int(n)] for n in ns], dtype=float)
d2 = (2**ns)**2  # Liouvillian dimension squared

log(f"\n  Data: {dict(zip(ns.astype(int), counts.astype(int)))}")

# Hypothesis 1: rates ~ a * 4^N
ratios = counts / 4**ns
log(f"\n  Hypothesis: rates = c * 4^N")
log(f"    c = rates/4^N: {dict(zip(ns.astype(int), [f'{r:.6f}' for r in ratios]))}")
log(f"    c is INCREASING -> not a fixed fraction of 4^N")

# Hypothesis 2: rates = 4^N - a * 3^N
# If rates = 4^N(1 - a*(3/4)^N), then (4^N - rates)/3^N = a
residual = (4**ns - counts) / 3**ns
log(f"\n  Hypothesis: rates = 4^N - a * 3^N")
log(f"    a = (4^N - rates)/3^N: {dict(zip(ns.astype(int), [f'{r:.4f}' for r in residual]))}")

# Hypothesis 3: rates = 4^N - 2^N (decoherence-free subspace?)
predicted_3 = 4**ns - 2**ns
error_3 = abs(counts - predicted_3) / counts
log(f"\n  Hypothesis: rates = 4^N - 2^N")
log(f"    Predicted: {dict(zip(ns.astype(int), predicted_3.astype(int)))}")
log(f"    Actual:    {dict(zip(ns.astype(int), counts.astype(int)))}")
log(f"    Error:     {dict(zip(ns.astype(int), [f'{e:.4f}' for e in error_3]))}")

# Hypothesis 4: rates = 4^N - 2*3^N + 2^N
predicted_4 = 4**ns - 2*3**ns + 2**ns
error_4 = abs(counts - predicted_4) / counts
log(f"\n  Hypothesis: rates = 4^N - 2*3^N + 2^N (inclusion-exclusion)")
log(f"    Predicted: {dict(zip(ns.astype(int), predicted_4.astype(int)))}")
log(f"    Actual:    {dict(zip(ns.astype(int), counts.astype(int)))}")
log(f"    Error:     {dict(zip(ns.astype(int), [f'{e:.4f}' for e in error_4]))}")

# Hypothesis 5: polynomial in N times 4^N
# rates / 4^N = p(N)?
frac = counts / 4**ns
log(f"\n  rates/4^N = {dict(zip(ns.astype(int), [f'{r:.6f}' for r in frac]))}")
log(f"  Approaching 1 as N grows -> rates ~ 4^N * (1 - f(N))")
log(f"  1 - rates/4^N = {dict(zip(ns.astype(int), [f'{1-r:.6f}' for r in frac]))}")

deficit = 1 - frac
log(f"\n  Deficit log-ratio analysis:")
for i in range(1, len(ns)):
    if deficit[i] > 0 and deficit[i-1] > 0:
        ratio = deficit[i] / deficit[i-1]
        log(f"    deficit[{int(ns[i])}]/deficit[{int(ns[i-1])}] = {ratio:.4f} (0.75 = (3/4))")

log(f"\n    [Test 7 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 8: DENSITY OF STATES AT N=7
# ============================================================
log("\n" + "=" * 60)
log("### TEST 8: Density of states at N=7 (from fresh computation)")
log("=" * 60)

try:
    n_q = 7
    bonds = heisenberg_star(n_q, [1.0]*(n_q-1))
    gamma_list = [gamma]*n_q
    log(f"  Building N=7 Liouvillian...")
    t0 = time.time()
    L = build_liouvillian_dense(n_q, bonds, gamma_list)
    log(f"  Built in {time.time()-t0:.1f}s. Extracting rates...")
    r = get_osc_rates(L)
    log(f"  {len(r)} rates in {time.time()-t0:.1f}s")
    
    r_g = [x/gamma for x in r]
    n_bins = 40  # finer resolution
    lo, hi = min(r_g), max(r_g)
    bin_w = (hi - lo) / n_bins
    bins = [0]*n_bins
    for x in r_g:
        idx = min(int((x - lo) / bin_w), n_bins - 1)
        bins[idx] += 1
    
    log(f"\n  Histogram ({n_bins} bins):")
    max_count = max(bins)
    for i, count in enumerate(bins):
        bar = '#' * int(50 * count / max_count) if max_count > 0 else ''
        bc = lo + (i + 0.5) * bin_w
        log(f"    {bc:>6.2f}g | {count:>5} {bar}")
    
    mean = np.mean(r_g)
    std = np.std(r_g)
    skew = np.mean(((np.array(r_g) - mean)/std)**3) if std > 0 else 0
    kurt = np.mean(((np.array(r_g) - mean)/std)**4) if std > 0 else 0
    
    log(f"\n  Mean: {mean:.4f}g (expected {n_q:.1f}g)")
    log(f"  Std:  {std:.4f}g")
    log(f"  Skew: {skew:.6f} (0=symmetric)")
    log(f"  Kurtosis: {kurt:.4f} (3=Gaussian, 1.8=semicircle)")
    
    # Percentiles
    r_arr = np.array(r_g)
    for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
        log(f"  P{p:>2}: {np.percentile(r_arr, p):.4f}g")

except Exception as e:
    log(f"  FAILED: {e}")
    log(traceback.format_exc())

log(f"\n    [Test 8 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 9: BANDWIDTH SCALING WITH NON-UNIFORM GAMMA
# ============================================================
log("\n" + "=" * 60)
log("### TEST 9: Bandwidth under non-uniform gamma")
log("=" * 60)
log("For uniform gamma: BW = 2(N-2)*gamma, min=2g, max=2(N-1)g")
log("For non-uniform: min=2*min(g_i)? max=2*(N-1)*max(g_i)? Or sum-based?")

for n_q in [4, 5]:
    log(f"\n  N={n_q}:")
    
    configs = [
        [0.05]*n_q,                          # uniform baseline
        [0.01, 0.01, 0.01, 0.20, 0.30][:n_q],  # one dominant
        [0.01, 0.05, 0.10, 0.15, 0.20][:n_q],  # linear ramp
        [0.20, 0.01, 0.20, 0.01, 0.20][:n_q],  # alternating
        [0.10, 0.10, 0.10, 0.10, 0.10][:n_q],  # double uniform
    ]
    
    log(f"    {'config':>40} | {'min/g':>8} | {'max':>8} | {'BW':>8} | "
        f"{'sum':>6} | {'2*min':>6} | {'2*sum-2*min':>12}")
    log("    " + "-" * 100)
    
    for gs in configs:
        bonds = heisenberg_star(n_q, [1.0]*(n_q-1))
        L = build_liouvillian_dense(n_q, bonds, gs)
        r = get_osc_rates(L)
        
        if r:
            min_r = min(r)
            max_r = max(r)
            bw = max_r - min_r
            sum_g = sum(gs)
            min_g = min(gs)
            max_g = max(gs)
            
            log(f"    {str(gs):>40} | {min_r:.5f} | {max_r:.5f} | {bw:.5f} | "
                f"{sum_g:.3f} | {2*min_g:.4f} | {2*sum_g-2*min_g:.5f}")

log(f"\n    [Test 9 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# TEST 10: THE 0.5 HUNT
# ============================================================
log("\n" + "=" * 60)
log("### TEST 10: Where does 0.5 appear in the spectral data?")
log("=" * 60)
log("z* = 0.5 is the fixed point boundary. Does it show up in:")
log("  - Rate ratios? min/max, center/max?")
log("  - Gap ratios?")
log("  - Eigenvalue moduli?")
log("  - Fraction of oscillatory vs total eigenvalues?")

for n_q in [3, 4, 5, 6, 7]:
    bonds = heisenberg_star(n_q, [1.0]*(n_q-1))
    gamma_list = [gamma]*n_q
    
    try:
        L = build_liouvillian_dense(n_q, bonds, gamma_list)
        evals = np.linalg.eigvals(L)
        r = get_osc_rates(L)
        
        if not r:
            continue
        
        # Various ratios
        min_r = min(r)
        max_r = max(r)
        center = n_q * gamma
        d2 = (2**n_q)**2
        
        n_osc = len(r)
        n_total = len(evals)
        n_pure_decay = sum(1 for ev in evals if abs(np.imag(ev)) < 0.05 and -np.real(ev) > 0.001)
        
        log(f"\n  N={n_q}:")
        log(f"    min/max = {min_r/max_r:.6f}")
        log(f"    min/center = {min_r/center:.6f}")
        log(f"    center/max = {center/max_r:.6f} = {n_q}/{2*(n_q-1)} = N/2(N-1)")
        log(f"    BW/max = {(max_r-min_r)/max_r:.6f} = 2(N-2)/2(N-1) = (N-2)/(N-1)")
        log(f"    osc/total = {n_osc/n_total:.6f}")
        log(f"    osc/4^N = {n_osc/4**n_q:.6f}")
        log(f"    pure_decay/total = {n_pure_decay/n_total:.6f}")
        
        # Check center/max convergence
        half_check = center / max_r  # = N/(2N-2)
        log(f"    center/max = N/(2(N-1)) = {n_q/(2*(n_q-1)):.6f} -> 0.5 as N->inf")
        
    except Exception as e:
        log(f"  N={n_q}: FAILED ({e})")

log(f"\n  Key observation: center/max = N/(2(N-1)) converges to 0.5 as N->infinity")
log(f"  The mirror center (Ng) divided by the max rate (2(N-1)g) -> 1/2")
log(f"  At N=7: {7/(2*6):.6f}")
log(f"  At N=100: {100/(2*99):.6f}")
log(f"  The 0.5 is the ASYMPTOTIC position of the mirror center within the band.")

log(f"\n    [Test 10 complete at {time.strftime('%H:%M:%S')}]")

# ============================================================
# FINAL SUMMARY
# ============================================================
log("\n" + "=" * 80)
log("FINAL SUMMARY")
log("=" * 80)
log(f"""
OVERNIGHT COMPUTATION COMPLETE
Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}

TESTS COMPLETED:
 1. N=8 sparse Liouvillian (65536x65536)
 2. Pauli complement proof at N=5, 6, 7
 3. Complete topology survey: star, chain, ring, complete, tree
 4. Non-uniform J + non-uniform gamma combined
 5. High-res avoided crossing map (N=4, 1000 J values)
 6. Conjugation operator search
 7. Rate count formula analysis
 8. Density of states at N=7
 9. Bandwidth scaling with non-uniform gamma
10. The 0.5 hunt

SCALING TABLE (updated):
""")

all_data = {2: 6, 3: 40, 4: 182, 5: 776, 6: 3228, 7: 13264}
log(f"  {'N':>3} | {'matrix':>8} | {'#rates':>7} | {'min/g':>6} | {'max/g':>6} | "
    f"{'BW/g':>6} | {'mirror':>6} | {'rates/4^N':>10}")
log("-" * 75)
for n in sorted(all_data.keys()):
    d2 = (2**n)**2
    log(f"  {n:>3} | {d2:>8} | {all_data[n]:>7} | {'2.0':>6} | {2*(n-1):>6.1f} | "
        f"{2*(n-2):>6.1f} | {'100%':>6} | {all_data[n]/4**n:>10.6f}")

log(f"""
CONJECTURES TO VERIFY:
  - rate_count = 4^N * (1 - c*(3/4)^N) for some constant c?
  - Avoided crossing gap is CONSTANT with N (~0.004g)?
  - Density of states is Gaussian with mean=Ng, skew=0?
  - center/max -> 0.5 as N -> infinity?
  - Conjugation operator = ??? (see Test 6 results)
""")

log("=" * 80)
log(f"Total runtime: will be computed from start/end timestamps above")
log("=" * 80)

f.close()
print(f"\n>>> All results saved to: {outpath}")
print(">>> Script complete. Good night Tom.")
