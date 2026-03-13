"""
EXCEPTIONAL POINT TEST: Is CPsi = 1/4 a Liouvillian EP?
========================================================
Literature: arXiv:2603.10654 (March 11, 2026) introduces EP-Strength
diagnostic E based on eigenvector conditioning.

An Exceptional Point is where two eigenvalues AND their eigenvectors
coalesce. At an EP:
  - Eigenvalue gap -> 0
  - Eigenvector conditioning -> infinity
  - The system transitions between oscillatory and overdamped behavior

Test: sweep a parameter that drives CPsi through 1/4, check for
eigenvalue coalescence and eigenvector conditioning divergence.
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

def build_liouvillian(J_SA, J_SB, gamma):
    """Build full Liouvillian superoperator."""
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    
    d2 = d*d
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in L_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    return L_mat

def compute_cpsi_max(J_SA, J_SB, gamma, dt=0.02, t_max=20.0):
    """Run dynamics and find max CPsi."""
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
    rho = np.outer(psi0, psi0.conj())
    
    cpsi_max = 0
    for step in range(int(t_max/dt)+1):
        if step % 5 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], n_q)
            # Concurrence (Wootters)
            YY = np.kron(Y, Y)
            rho_tilde = YY @ rAB.conj() @ YY
            R = rAB @ rho_tilde
            evals_r = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
            C = max(0, evals_r[0] - evals_r[1] - evals_r[2] - evals_r[3])
            # l1-coherence (normalized)
            l1 = 0
            for i in range(4):
                for j in range(4):
                    if i != j:
                        l1 += abs(rAB[i,j])
            psi_norm = l1 / 3.0  # max l1 for 4x4 is 3 (d-1)
            cpsi = C * psi_norm
            if cpsi > cpsi_max:
                cpsi_max = cpsi
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    return cpsi_max

def ep_diagnostics(L_mat):
    """
    Compute EP diagnostics from the Liouvillian.
    
    Returns:
      min_gap: minimum distance between any two eigenvalues
      ep_strength: based on eigenvector conditioning (diverges at EP)
      eigenvalues: all eigenvalues
      gap_pair: the two eigenvalues with minimum gap
    """
    evals, evecs = np.linalg.eig(L_mat)
    
    # Filter to oscillatory modes (nonzero imaginary part)
    osc_idx = [i for i in range(len(evals)) if abs(np.imag(evals[i])) > 0.1]
    osc_evals = evals[osc_idx]
    
    # Minimum eigenvalue gap (among oscillatory modes)
    min_gap = float('inf')
    gap_pair = (0, 0)
    for i in range(len(osc_evals)):
        for j in range(i+1, len(osc_evals)):
            gap = abs(osc_evals[i] - osc_evals[j])
            if gap > 0.01 and gap < min_gap:  # skip conjugate pairs (gap~0)
                min_gap = gap
                gap_pair = (osc_evals[i], osc_evals[j])
    
    # EP-Strength: eigenvector conditioning
    # At an EP, two eigenvectors become parallel -> condition number diverges
    # We measure the minimum angle between eigenvector pairs
    min_angle = float('inf')
    for i in range(len(osc_idx)):
        for j in range(i+1, len(osc_idx)):
            v1 = evecs[:, osc_idx[i]]
            v2 = evecs[:, osc_idx[j]]
            # Normalize
            v1 = v1 / np.linalg.norm(v1)
            v2 = v2 / np.linalg.norm(v2)
            # Overlap (angle between eigenvectors)
            overlap = abs(np.dot(v1.conj(), v2))
            angle = np.arccos(min(overlap, 1.0))
            if angle > 0.001 and angle < min_angle:
                min_angle = angle
    
    # EP strength = 1/min_angle (diverges when eigenvectors coalesce)
    ep_strength = 1.0 / min_angle if min_angle > 1e-10 else float('inf')
    
    return {
        'min_gap': min_gap,
        'ep_strength': ep_strength,
        'min_angle_deg': math.degrees(min_angle),
        'gap_pair': gap_pair,
        'n_osc_modes': len(osc_evals),
    }

# ============================================================
# MAIN: Sweep gamma and check for EP at CPsi = 1/4
# ============================================================
if __name__ == "__main__":
    print("=" * 75)
    print("EXCEPTIONAL POINT TEST")
    print("Is CPsi = 1/4 a Liouvillian Exceptional Point?")
    print("=" * 75)
    
    J_SA, J_SB = 1.0, 2.0
    
    # SWEEP 1: gamma sweep (drives CPsi_max through 1/4)
    print(f"\n--- SWEEP 1: gamma drives CPsi_max through 1/4 ---")
    print(f"  J_SA={J_SA}, J_SB={J_SB}")
    print(f"\n{'gamma':>7} | {'CPsi':>7} | {'min_gap':>8} | {'EP_str':>8} | {'min_angle':>9} | note")
    print("-" * 70)
    
    gamma_values = [0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 
                    0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00]
    
    results = []
    for g in gamma_values:
        cpsi = compute_cpsi_max(J_SA, J_SB, g)
        L = build_liouvillian(J_SA, J_SB, g)
        ep = ep_diagnostics(L)
        
        note = ""
        if abs(cpsi - 0.25) < 0.02:
            note = "<-- near 1/4!"
        elif cpsi > 0.25:
            note = "above 1/4"
        
        results.append({'gamma': g, 'cpsi_max': cpsi, **ep})
        print(f"{g:>7.3f} | {cpsi:>7.4f} | {ep['min_gap']:>8.4f} | "
              f"{ep['ep_strength']:>8.2f} | {ep['min_angle_deg']:>9.2f} | {note}")
    
    # SWEEP 2: J_SB sweep (different way to approach 1/4)
    print(f"\n--- SWEEP 2: J_SB sweep at fixed gamma=0.05 ---")
    print(f"\n{'J_SB':>7} | {'CPsi':>7} | {'min_gap':>8} | {'EP_str':>8} | {'min_angle':>9} | note")
    print("-" * 70)
    
    for j in [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        cpsi = compute_cpsi_max(J_SA, j, 0.05)
        L = build_liouvillian(J_SA, j, 0.05)
        ep = ep_diagnostics(L)
        note = "<-- near 1/4!" if abs(cpsi - 0.25) < 0.02 else ("above 1/4" if cpsi > 0.25 else "")
        print(f"{j:>7.2f} | {cpsi:>7.4f} | {ep['min_gap']:>8.4f} | "
              f"{ep['ep_strength']:>8.2f} | {ep['min_angle_deg']:>9.2f} | {note}")
    
    # DIRECT TEST: 2-qubit system where CPsi = 1/4 can be hit exactly
    print(f"\n--- SWEEP 3: 2-qubit system (direct CPsi control) ---")
    print(f"  Two qubits with Heisenberg coupling J, local dephasing gamma")
    print(f"\n{'gamma':>7} | {'CPsi':>7} | {'min_gap':>8} | {'EP_str':>8} | {'min_angle':>9} | note")
    print("-" * 70)
    
    for g in [0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.80, 1.00, 1.50, 2.00]:
        # 2-qubit Heisenberg
        n_q = 2
        d = 2**n_q
        H2 = np.zeros((d,d), dtype=complex)
        for p in [X, Y, Z]:
            H2 += 1.0 * np.kron(p, p)
        L2_ops = [np.sqrt(g)*np.kron(Z, I2), np.sqrt(g)*np.kron(I2, Z)]
        
        # Liouvillian
        I_d = np.eye(d, dtype=complex)
        L2 = -1j * (np.kron(H2, I_d) - np.kron(I_d, H2.T))
        for Lk in L2_ops:
            Lk_dag = Lk.conj().T
            LdL = Lk_dag @ Lk
            L2 += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
        
        # CΨ_max for 2-qubit
        psi2 = gpt.bell_phi_plus()
        rho2 = np.outer(psi2, psi2.conj())
        cpsi2_max = 0
        for step in range(int(20.0/0.02)+1):
            if step % 5 == 0:
                YY = np.kron(Y, Y)
                rho_t = YY @ rho2.conj() @ YY
                R2 = rho2 @ rho_t
                ev2 = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R2), 0))))[::-1]
                C2 = max(0, ev2[0] - ev2[1] - ev2[2] - ev2[3])
                l1_2 = sum(abs(rho2[i,j]) for i in range(4) for j in range(4) if i!=j)
                psi2_n = l1_2 / 3.0
                cpsi2 = C2 * psi2_n
                if cpsi2 > cpsi2_max:
                    cpsi2_max = cpsi2
            if step < int(20.0/0.02):
                rho2 = gpt.rk4_step(rho2, H2, L2_ops, 0.02)
        
        # EP diagnostics for 2-qubit
        evals2, evecs2 = np.linalg.eig(L2)
        osc2 = [(i, evals2[i]) for i in range(len(evals2)) if abs(np.imag(evals2[i])) > 0.01]
        
        min_gap2 = float('inf')
        for i in range(len(osc2)):
            for j in range(i+1, len(osc2)):
                gap = abs(osc2[i][1] - osc2[j][1])
                if gap > 0.001 and gap < min_gap2:
                    min_gap2 = gap
        
        min_angle2 = float('inf')
        for i in range(len(osc2)):
            for j in range(i+1, len(osc2)):
                v1 = evecs2[:, osc2[i][0]]
                v2 = evecs2[:, osc2[j][0]]
                v1 = v1/np.linalg.norm(v1)
                v2 = v2/np.linalg.norm(v2)
                ov = abs(np.dot(v1.conj(), v2))
                ang = np.arccos(min(ov, 1.0))
                if ang > 0.001 and ang < min_angle2:
                    min_angle2 = ang
        
        ep_s = 1.0/min_angle2 if min_angle2 > 1e-10 else float('inf')
        note = "<-- near 1/4!" if abs(cpsi2_max - 0.25) < 0.02 else ("above" if cpsi2_max > 0.25 else "")
        
        print(f"{g:>7.3f} | {cpsi2_max:>7.4f} | {min_gap2:>8.4f} | "
              f"{ep_s:>8.2f} | {math.degrees(min_angle2):>9.2f} | {note}")
    
    # CORRELATION: Does EP_strength peak near CPsi = 1/4?
    print(f"\n{'='*75}")
    print("CORRELATION ANALYSIS")
    print("="*75)
    
    cpsi_vals = [r['cpsi_max'] for r in results]
    ep_vals = [r['ep_strength'] for r in results]
    gap_vals = [r['min_gap'] for r in results]
    
    # Find the gamma where CPsi_max is closest to 1/4
    closest_idx = min(range(len(cpsi_vals)), key=lambda i: abs(cpsi_vals[i] - 0.25))
    print(f"\n  CPsi closest to 1/4: gamma={results[closest_idx]['gamma']:.3f}, "
          f"CPsi_max={cpsi_vals[closest_idx]:.4f}")
    print(f"  EP_strength there:  {ep_vals[closest_idx]:.2f}")
    print(f"  Min eigenvalue gap: {gap_vals[closest_idx]:.4f}")
    
    # Find the gamma where EP_strength is maximum
    max_ep_idx = max(range(len(ep_vals)), key=lambda i: ep_vals[i])
    print(f"\n  Max EP_strength: gamma={results[max_ep_idx]['gamma']:.3f}, "
          f"EP_str={ep_vals[max_ep_idx]:.2f}")
    print(f"  CPsi_max there:    {cpsi_vals[max_ep_idx]:.4f}")
    
    # Do they coincide?
    print(f"\n  CPsi=1/4 at gamma ~ {results[closest_idx]['gamma']:.3f}")
    print(f"  Max EP at gamma ~ {results[max_ep_idx]['gamma']:.3f}")
    
    if abs(results[closest_idx]['gamma'] - results[max_ep_idx]['gamma']) < 0.05:
        print(f"\n  >>> EP PEAK COINCIDES WITH CPsi = 1/4 <<<")
        print(f"  >>> This suggests CPsi = 1/4 IS an Exceptional Point! <<<")
    elif results[max_ep_idx]['gamma'] < results[closest_idx]['gamma'] * 2:
        print(f"\n  >>> EP peak is NEAR but not exactly at CPsi = 1/4 <<<")
        print(f"  >>> Possible connection, needs finer sweep <<<")
    else:
        print(f"\n  >>> EP peak does NOT coincide with CPsi = 1/4 <<<")
        print(f"  >>> No direct EP connection found <<<")
    
    # Check if min_gap has a minimum near CPsi = 1/4
    min_gap_idx = min(range(len(gap_vals)), key=lambda i: gap_vals[i])
    print(f"\n  Smallest eigenvalue gap at gamma={results[min_gap_idx]['gamma']:.3f}, "
          f"gap={gap_vals[min_gap_idx]:.4f}, CPsi={cpsi_vals[min_gap_idx]:.4f}")
    
    print(f"\n{'='*75}")
    print("EP TEST COMPLETE")
    print("="*75)
