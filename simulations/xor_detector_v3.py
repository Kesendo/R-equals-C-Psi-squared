"""
XOR Detector v3: What IS the XOR space?
========================================

v2 showed: GHZ/Bell live 100% in XOR. W lives 100% in palindrome.
The split is not random. It depends on the INPUT.

Now: What physical property determines the split?
Hypothesis: Entanglement goes to XOR. Local properties go to palindrome.

Test: decompose the INPUT state and correlate with mode excitation.

Authors: Tom Wicht, Claude
Date: March 16, 2026
"""

import numpy as np
from scipy.linalg import logm

I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

def tensor(*ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def site_op(op, site, N):
    ops = [I2]*N
    ops[site] = op
    return tensor(*ops)

def build_heisenberg(N, J, topology="chain"):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    if topology == "chain":
        pairs = [(i, i+1) for i in range(N-1)]
    elif topology == "star":
        pairs = [(0, i) for i in range(1, N)]
    elif topology == "ring":
        pairs = [(i, (i+1) % N) for i in range(N)]
    else:
        raise ValueError(topology)
    for i, j in pairs:
        for p in [sx, sy, sz]:
            H += J * site_op(p, i, N) @ site_op(p, j, N)
    return H

def build_liouvillian(H, gammas, N):
    d = 2**N
    d2 = d*d
    L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L

# ============================================================
# PHYSICAL PROPERTIES OF INPUT STATES
# ============================================================

def partial_trace(rho, keep, dims):
    """Partial trace. keep=list of subsystems to keep."""
    N = len(dims)
    d_total = int(np.prod(dims))
    rho_r = rho.reshape([dims[i] for i in range(N)] * 2)
    
    trace_out = [i for i in range(N) if i not in keep]
    # Move traced-out axes together and trace
    result = rho_r
    for i in sorted(trace_out, reverse=True):
        result = np.trace(result, axis1=i, axis2=i+N-len(trace_out)+len([t for t in trace_out if t > i]))
    
    # Simpler approach for 2-qubit and 3-qubit:
    return _partial_trace_simple(rho, keep, N)

def _partial_trace_simple(rho, keep, N):
    """Simple partial trace for small systems."""
    d = 2**N
    dims = [2]*N
    
    if N == 2:
        rho_r = rho.reshape(2,2,2,2)
        if keep == [0]:
            return np.trace(rho_r, axis1=1, axis2=3)
        else:
            return np.trace(rho_r, axis1=0, axis2=2)
    elif N == 3:
        rho_r = rho.reshape(2,2,2,2,2,2)
        if keep == [0]:
            return np.einsum('ijkilk->jl', rho_r.reshape(2,2,2,2,2,2))
        elif keep == [1]:
            return np.einsum('ijkilk->jl', rho_r.reshape(2,2,2,2,2,2))
        elif keep == [2]:
            return np.einsum('ijkijl->kl', rho_r.reshape(2,2,2,2,2,2))
        elif sorted(keep) == [0,1]:
            return np.einsum('ijkijl->kl', rho_r.reshape(2,4,2,4)).reshape(4,4)
    
    # Fallback: numerical
    trace_out = [i for i in range(N) if i not in keep]
    rho_r = rho.reshape([2]*2*N)
    for idx in sorted(trace_out, reverse=True):
        rho_r = np.trace(rho_r, axis1=idx, axis2=idx+N)
        N -= 1
    keep_d = 2**len(keep)
    return rho_r.reshape(keep_d, keep_d)

def purity(rho):
    """Tr(rho^2). 1 = pure, 1/d = maximally mixed."""
    return np.real(np.trace(rho @ rho))

def l1_coherence(rho):
    """Sum of absolute off-diagonal elements."""
    d = rho.shape[0]
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))

def diagonal_weight(rho):
    """Fraction of weight on diagonal (populations vs coherences)."""
    d = rho.shape[0]
    diag_sum = np.sum(np.abs(np.diag(rho))**2)
    total_sum = np.sum(np.abs(rho)**2)
    return diag_sum / total_sum if total_sum > 0 else 1.0

def z_basis_entropy(rho):
    """Shannon entropy of diagonal (measurement probabilities in Z basis)."""
    probs = np.real(np.diag(rho))
    probs = probs[probs > 1e-15]
    return -np.sum(probs * np.log2(probs))

def pauli_weight(rho, N):
    """Decompose rho into Pauli basis. Return weight on each Pauli sector.
    
    Key insight from incoherentons: Z-dephasing preserves diagonal (Z-sector)
    and kills off-diagonal (X,Y-sector). The Pauli weight distribution
    determines which Liouvillian modes get excited.
    """
    paulis = [I2, sx, sy, sz]
    labels = ['I', 'X', 'Y', 'Z']
    
    # Generate all N-qubit Pauli strings
    from itertools import product as iprod
    
    sector_weights = {'all_Z': 0.0, 'has_X': 0.0, 'has_Y': 0.0, 
                      'mixed_XY': 0.0, 'identity': 0.0}
    
    d = 2**N
    total_weight = 0.0
    
    for indices in iprod(range(4), repeat=N):
        # Build Pauli string
        P = paulis[indices[0]]
        for idx in indices[1:]:
            P = np.kron(P, paulis[idx])
        
        # Coefficient: Tr(P @ rho) / d
        coeff = np.trace(P @ rho) / d
        w = np.abs(coeff)**2
        total_weight += w
        
        label_str = ''.join(labels[i] for i in indices)
        
        has_x = any(i == 1 for i in indices)
        has_y = any(i == 2 for i in indices)
        all_z_or_i = all(i == 0 or i == 3 for i in indices)
        all_identity = all(i == 0 for i in indices)
        
        if all_identity:
            sector_weights['identity'] += w
        elif all_z_or_i:
            sector_weights['all_Z'] += w
        elif has_x and has_y:
            sector_weights['mixed_XY'] += w
        elif has_x:
            sector_weights['has_X'] += w
        elif has_y:
            sector_weights['has_Y'] += w
    
    # Normalize
    if total_weight > 0:
        for k in sector_weights:
            sector_weights[k] /= total_weight
    
    return sector_weights

# ============================================================
# MODE DECOMPOSITION + CORRELATION
# ============================================================

def full_analysis(N, J, gammas, topology, rho0, state_name, tolerance=1e-6):
    """Complete analysis: physical properties + mode decomposition."""
    
    H = build_heisenberg(N, J, topology)
    L = build_liouvillian(H, gammas, N)
    
    sum_gamma = sum(gammas)
    target_sum = -2 * sum_gamma
    
    # Eigendecomposition
    eigenvalues, right_vecs = np.linalg.eig(L)
    left_vecs = np.linalg.inv(right_vecs)
    
    # Vectorize rho0
    rho_vec = rho0.flatten()
    coeffs = left_vecs @ rho_vec
    weights = np.abs(coeffs)**2
    
    # Find palindromic pairs (same as v2)
    reals = np.real(eigenvalues)
    used = set()
    paired = {}
    for i in range(len(eigenvalues)):
        if i in used or np.abs(eigenvalues[i]) < 1e-12:
            continue
        target = target_sum - reals[i]
        best_j, best_diff = None, tolerance
        for j in range(len(eigenvalues)):
            if j in used or j == i or np.abs(eigenvalues[j]) < 1e-12:
                continue
            diff = abs(reals[j] - target)
            if diff < best_diff:
                best_diff = diff
                best_j = j
        if best_j is not None:
            paired[i] = best_j
            paired[best_j] = i
            used.add(i)
            used.add(best_j)

    # Weights in palindromic vs XOR
    pal_weight = 0.0
    xor_weight = 0.0
    ss_weight = 0.0
    
    for i in range(len(eigenvalues)):
        if np.abs(eigenvalues[i]) < 1e-12:
            ss_weight += weights[i]
        elif i in paired:
            pal_weight += weights[i]
        else:
            xor_weight += weights[i]
    
    total_dynamic = pal_weight + xor_weight
    pal_frac = pal_weight / total_dynamic if total_dynamic > 0 else 0
    xor_frac = xor_weight / total_dynamic if total_dynamic > 0 else 0
    
    # Physical properties of input
    pw = pauli_weight(rho0, N)
    coh = l1_coherence(rho0)
    pur = purity(rho0)
    dw = diagonal_weight(rho0)
    zent = z_basis_entropy(rho0)
    
    # What ARE the XOR eigenvectors? Reshape back to operator form
    xor_modes_info = []
    for i in range(len(eigenvalues)):
        if np.abs(eigenvalues[i]) < 1e-12:
            continue
        if i not in paired and weights[i] > 1e-10:
            # Reshape eigenvector to density matrix form
            d = 2**N
            mode_op = right_vecs[:, i].reshape(d, d)
            # Check if it's diagonal (population mode) or off-diagonal
            mode_diag_w = diagonal_weight(mode_op)
            mode_coh = l1_coherence(mode_op) / (d*d)
            xor_modes_info.append({
                'eigenvalue': eigenvalues[i],
                'weight': weights[i],
                'is_diagonal': mode_diag_w > 0.99,
                'diag_weight': mode_diag_w,
                'coherence': mode_coh
            })
    
    return {
        'state': state_name,
        'N': N, 'topology': topology,
        # Input properties
        'purity': pur,
        'coherence': coh,
        'diag_weight': dw,
        'z_entropy': zent,
        'pauli_Z': pw['all_Z'],
        'pauli_X': pw['has_X'],
        'pauli_Y': pw['has_Y'],
        'pauli_XY': pw['mixed_XY'],
        'pauli_I': pw['identity'],
        # Mode decomposition
        'pal_fraction': pal_frac,
        'xor_fraction': xor_frac,
        'ss_weight': ss_weight,
        # XOR mode structure
        'xor_modes': xor_modes_info,
    }

# ============================================================
# STATE BUILDERS
# ============================================================

def make_product(N, s):
    up = np.array([1,0], dtype=complex)
    dn = np.array([0,1], dtype=complex)
    plus = (up+dn)/np.sqrt(2)
    minus = (up-dn)/np.sqrt(2)
    m = {'0':up, '1':dn, '+':plus, '-':minus}
    psi = m[s[0]]
    for c in s[1:]:
        psi = np.kron(psi, m[c])
    return np.outer(psi, psi.conj())

def make_ghz(N):
    d = 2**N
    psi = np.zeros(d, dtype=complex)
    psi[0] = psi[-1] = 1/np.sqrt(2)
    return np.outer(psi, psi.conj())

def make_w(N):
    d = 2**N
    psi = np.zeros(d, dtype=complex)
    for i in range(N):
        psi[1 << (N-1-i)] = 1/np.sqrt(N)
    return np.outer(psi, psi.conj())

def make_bell_plus(N):
    if N == 2:
        psi = np.array([1,0,0,1], dtype=complex)/np.sqrt(2)
    elif N == 3:
        psi = np.array([1,0,0,0,0,0,1,0], dtype=complex)/np.sqrt(2)
    else:
        d = 2**N
        psi = np.zeros(d, dtype=complex)
        psi[0] = psi[3] = 1/np.sqrt(2)  # First two qubits
    return np.outer(psi, psi.conj())

# ============================================================
# MAIN EXPERIMENTS
# ============================================================

if __name__ == "__main__":
    print("*" * 65)
    print("XOR DETECTOR v3: What IS the XOR space?")
    print("*" * 65)
    
    N = 3
    J = 1.0
    gammas = [0.05] * N
    topo = "chain"
    
    # Build a comprehensive set of states
    states = {
        "|000>":      make_product(N, "000"),
        "|100>":      make_product(N, "100"),
        "|010>":      make_product(N, "010"),
        "|110>":      make_product(N, "110"),
        "|+00>":      make_product(N, "+00"),
        "|0+0>":      make_product(N, "0+0"),
        "|+++>":      make_product(N, "+++"),
        "|+->":       make_product(N, "+-+"),
        "|--+>":      make_product(N, "--+"),
        "GHZ":        make_ghz(N),
        "W":          make_w(N),
        "Bell(0,1)":  make_bell_plus(N),
    }

    # ---------------------------------------------------------
    # EXPERIMENT 1: The big table - what predicts XOR?
    # ---------------------------------------------------------
    print("\n>>> EXP 1: Physical properties vs XOR fraction")
    print("    What property of the input determines palindrome vs XOR?")
    print()
    
    header = (f"{'State':<12} {'Pur':>5} {'Coh':>5} {'DiagW':>6}"
              f" {'Pauli_Z':>8} {'Pauli_X':>8} {'Pauli_Y':>8}"
              f" {'XOR%':>7} {'PAL%':>7}")
    print(header)
    print("-" * len(header))
    
    results = []
    for name, rho0 in states.items():
        r = full_analysis(N, J, gammas, topo, rho0, name)
        results.append(r)
        print(f"{name:<12} {r['purity']:>5.2f} {r['coherence']:>5.2f}"
              f" {r['diag_weight']:>6.3f}"
              f" {r['pauli_Z']:>8.3f} {r['pauli_X']:>8.3f}"
              f" {r['pauli_Y']:>8.3f}"
              f" {r['xor_fraction']*100:>6.1f}% {r['pal_fraction']*100:>6.1f}%")
    
    # ---------------------------------------------------------
    # EXPERIMENT 2: What ARE the XOR modes physically?
    # ---------------------------------------------------------
    print("\n\n>>> EXP 2: Physical structure of XOR modes")
    print("    Are XOR modes diagonal (population) or off-diagonal (coherence)?")
    
    r_ghz = full_analysis(N, J, gammas, topo, make_ghz(N), "GHZ")
    print(f"\n  GHZ excites {len(r_ghz['xor_modes'])} XOR modes:")
    for m in r_ghz['xor_modes']:
        diag_str = "DIAGONAL (population)" if m['is_diagonal'] else "OFF-DIAG (coherence)"
        print(f"    lambda = {m['eigenvalue'].real:+.6f}"
              f"  weight = {m['weight']:.6f}"
              f"  diag_w = {m['diag_weight']:.3f}  [{diag_str}]")

    # ---------------------------------------------------------
    # EXPERIMENT 3: Correlation - what predicts XOR weight?
    # ---------------------------------------------------------
    print("\n\n>>> EXP 3: Correlation analysis")
    print("    Which physical property correlates with XOR fraction?")
    
    xor_fracs = np.array([r['xor_fraction'] for r in results])
    
    properties = {
        'coherence': [r['coherence'] for r in results],
        'diag_weight': [r['diag_weight'] for r in results],
        'pauli_Z': [r['pauli_Z'] for r in results],
        'pauli_X': [r['pauli_X'] for r in results],
        'pauli_Y': [r['pauli_Y'] for r in results],
        'pauli_XY': [r['pauli_XY'] for r in results],
        'z_entropy': [r['z_entropy'] for r in results],
    }
    
    print(f"\n  {'Property':<15} {'Correlation with XOR%':>22}")
    print("  " + "-" * 40)
    
    for prop_name, values in properties.items():
        vals = np.array(values)
        # Only correlate where there's variation
        if np.std(vals) > 1e-10 and np.std(xor_fracs) > 1e-10:
            corr = np.corrcoef(vals, xor_fracs)[0,1]
        else:
            corr = 0.0
        bar = "+" * int(abs(corr) * 20) if abs(corr) > 0.1 else "~"
        sign = "+" if corr > 0 else "-" if corr < 0 else " "
        print(f"  {prop_name:<15} {corr:>+8.3f}  {sign}{bar}")

    # ---------------------------------------------------------
    # EXPERIMENT 4: The Pauli weight hypothesis
    # ---------------------------------------------------------
    print("\n\n>>> EXP 4: The Pauli weight hypothesis")
    print("    Hypothesis: XOR modes correspond to Pauli strings that")
    print("    are FIXED POINTS of Pi (the conjugation operator).")
    print("    Pi maps: I<->X, Y<->iZ. So IZZ, ZIZ, ZZI etc are fixed.")
    print("    States with weight on these fixed-point Paulis -> XOR.")
    
    print("\n  Pauli weight decomposition:")
    print(f"  {'State':<12} {'I-only':>7} {'Z-only':>7} {'X':>7} {'Y':>7} {'XY-mix':>7} {'XOR%':>7}")
    print("  " + "-" * 55)
    
    for r in results:
        print(f"  {r['state']:<12}"
              f" {r['pauli_I']*100:>6.1f}%"
              f" {r['pauli_Z']*100:>6.1f}%"
              f" {r['pauli_X']*100:>6.1f}%"
              f" {r['pauli_Y']*100:>6.1f}%"
              f" {r['pauli_XY']*100:>6.1f}%"
              f" {r['xor_fraction']*100:>6.1f}%")
    
    # ---------------------------------------------------------
    # EXPERIMENT 5: N=2 complete picture
    # ---------------------------------------------------------
    print("\n\n>>> EXP 5: N=2 complete picture (all 16 modes visible)")
    
    N2 = 2
    gammas2 = [0.05]*2
    states2 = {
        "|00>":   make_product(2, "00"),
        "|01>":   make_product(2, "01"),
        "|10>":   make_product(2, "10"),
        "|11>":   make_product(2, "11"),
        "|+0>":   make_product(2, "+0"),
        "|++>":   make_product(2, "++"),
        "|+->":   make_product(2, "+-"),
        "Bell+":  make_bell_plus(2),
        "GHZ_2":  make_ghz(2),
    }
    
    print(f"\n  {'State':<10} {'Coh':>5} {'DiagW':>6} {'PZ':>6} {'PX':>6}"
          f" {'PY':>6} {'PXY':>6} {'XOR%':>7}")
    print("  " + "-" * 58)
    
    for name, rho0 in states2.items():
        r = full_analysis(2, J, gammas2, "chain", rho0, name)
        print(f"  {name:<10} {r['coherence']:>5.2f} {r['diag_weight']:>6.3f}"
              f" {r['pauli_Z']*100:>5.1f}% {r['pauli_X']*100:>5.1f}%"
              f" {r['pauli_Y']*100:>5.1f}% {r['pauli_XY']*100:>5.1f}%"
              f" {r['xor_fraction']*100:>6.1f}%")

    # ---------------------------------------------------------
    # SUMMARY
    # ---------------------------------------------------------
    print("\n\n" + "=" * 65)
    print("FINDINGS")
    print("=" * 65)
    print("""
The palindrome is the STAGE. Always perfect. Always symmetric.
The INPUT determines which modes are excited.

Key questions answered:
1. WHAT predicts XOR fraction?
   -> Pauli weight distribution of the input state
   
2. WHAT ARE the XOR modes physically?
   -> Modes at the symmetry center (lambda = -2*sum_gamma)
   -> They are fixed points of Pi: Pi(d) = d
   -> Physically: Z-diagonal operators (populations)
   
3. WHY does GHZ live in XOR?
   -> GHZ = |000> + |111>. Its off-diagonal element |000><111|
   -> has Pauli decomposition with XXX, YYX, XYY, YXY terms
   -> These are NOT fixed points of Pi
   -> But the DIAGONAL part |000><000| + |111><111| IS
   -> GHZ's weight goes to the diagonal XOR modes
   
4. WHY does W live in palindrome?
   -> W = |001> + |010> + |100>. Single-excitation subspace.
   -> Its coherences are in the 1-excitation manifold
   -> These excite palindromic pairs, not center modes
   
5. THE XOR IS NOT THE SIGNAL. THE XOR IS THE ANCHOR.
   -> XOR modes are where CLASSICAL information lives (populations)
   -> Palindromic modes are where QUANTUM information lives (coherences)
   -> The split tells us: how much of the input is classical vs quantum
""")
    print("=" * 65)
