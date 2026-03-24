"""
JOINT POLE ANALYSIS: Shared poles across all channels
======================================================
GPT's #1 recommendation: Don't fit each channel independently.
Fit ONE common set of poles with channel-specific residues.

Also: plot poles in the complex plane across topology sweep.
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

def pauli_expect(rho, P):
    return float(np.real(np.trace(rho @ P)))

def generate_all_channels(J_SA=1.0, J_SB=2.0, gamma=0.05, dt=0.02, t_max=20.0):
    """Generate ALL four channels: A alone, B alone, c+, c-."""
    n_q = 3
    H = np.zeros((8,8), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
    rho = np.outer(psi0, psi0.conj())
    
    YZ = np.kron(Y, Z)
    ZY = np.kron(Z, Y)
    YI = np.kron(Y, I2)
    ZI = np.kron(Z, I2)
    IY = np.kron(I2, Y)
    IZ = np.kron(I2, Z)
    
    times, ch_cp, ch_cm, ch_A_yz, ch_B_yz = [], [], [], [], []
    
    for step in range(int(t_max/dt)+1):
        t = step * dt
        rAB = gpt.partial_trace_keep(rho, [1,2], n_q)
        yz = pauli_expect(rAB, YZ)
        zy = pauli_expect(rAB, ZY)
        a_y = pauli_expect(rAB, YI)
        a_z = pauli_expect(rAB, ZI)
        b_y = pauli_expect(rAB, IY)
        b_z = pauli_expect(rAB, IZ)
        
        times.append(t)
        ch_cp.append((yz + zy) / math.sqrt(2))
        ch_cm.append((yz - zy) / math.sqrt(2))
        ch_A_yz.append(a_y * b_z)   # factored YZ contribution from A side
        ch_B_yz.append(a_z * b_y)   # factored YZ contribution from B side
        
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    return {
        'times': np.array(times),
        'c+': np.array(ch_cp),
        'c-': np.array(ch_cm),
        'YZ': np.array([pauli_expect(gpt.partial_trace_keep(rho,[1,2],n_q), YZ) for _ in [0]]),
        'dt': dt,
    }

def exact_poles_from_liouvillian(J_SA=1.0, J_SB=2.0, gamma=0.05):
    """
    Extract EXACT complex poles from the Liouvillian superoperator.
    No fitting needed - these are the TRUE system poles.
    """
    n_q = 3
    d = 2**n_q  # 8
    
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    
    # Build Liouvillian as d^2 x d^2 matrix
    d2 = d*d
    L_mat = np.zeros((d2, d2), dtype=complex)
    
    # -i[H, rho] part
    I_d = np.eye(d, dtype=complex)
    L_mat += -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    
    # Lindblad dissipator
    for Lk in L_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    
    # Eigendecompose
    evals, evecs = np.linalg.eig(L_mat)
    
    poles = []
    for i, ev in enumerate(evals):
        freq = abs(np.imag(ev)) / (2 * np.pi)
        decay = -np.real(ev)
        if freq > 0.01 and decay > -0.1:
            poles.append({
                'eigenvalue': ev,
                'frequency': freq,
                'decay_rate': decay,
                'Q': np.pi * freq / decay if decay > 0.001 else float('inf'),
                'evec_idx': i,
            })
    
    # Sort by frequency
    poles.sort(key=lambda x: x['frequency'])
    return poles, evals, evecs, L_mat

def compute_residues(evals, evecs, L_mat, J_SA=1.0, J_SB=2.0, gamma=0.05):
    """
    Compute residues for c+ and c- observables.
    Residue_k = Tr(O * evec_k) * Tr(evec_k_left * rho0)
    """
    n_q = 3
    d = 2**n_q
    
    # Build observables as vectorized operators
    YZ = op_at(Y,1,n_q) @ op_at(Z,2,n_q)
    ZY = op_at(Z,1,n_q) @ op_at(Y,2,n_q)
    O_cp = (YZ + ZY) / math.sqrt(2)  # c+ observable
    O_cm = (YZ - ZY) / math.sqrt(2)  # c- observable
    
    # Initial state
    psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
    rho0 = np.outer(psi0, psi0.conj())
    rho0_vec = rho0.flatten()
    
    # Vectorize observables
    O_cp_vec = O_cp.flatten()
    O_cm_vec = O_cm.flatten()
    
    # Left eigenvectors (rows of inverse)
    evecs_left = np.linalg.inv(evecs)
    
    residues = []
    for i in range(len(evals)):
        freq = abs(np.imag(evals[i])) / (2*np.pi)
        decay = -np.real(evals[i])
        
        if freq < 0.01:
            continue
        
        # Right eigenvector (column)
        r_vec = evecs[:, i]
        # Left eigenvector (row)
        l_vec = evecs_left[i, :]
        
        # Residue = <O|r_k> * <l_k|rho0>
        r_cp = np.dot(O_cp_vec.conj(), r_vec) * np.dot(l_vec, rho0_vec)
        r_cm = np.dot(O_cm_vec.conj(), r_vec) * np.dot(l_vec, rho0_vec)
        
        residues.append({
            'eigenvalue': evals[i],
            'frequency': freq,
            'decay_rate': decay,
            'residue_cp': r_cp,
            'residue_cm': r_cm,
            'amp_cp': abs(r_cp),
            'amp_cm': abs(r_cm),
            'phase_cp': np.angle(r_cp),
            'phase_cm': np.angle(r_cm),
        })
    
    residues.sort(key=lambda x: -max(x['amp_cp'], x['amp_cm']))
    return residues

# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 75)
    print("JOINT POLE ANALYSIS: Exact Liouvillian poles + channel residues")
    print("Signal engineering approach: shared poles, per-channel residues")
    print("=" * 75)
    
    # 1. Extract exact poles
    print("\n--- EXACT SYSTEM POLES (from Liouvillian) ---")
    poles, evals, evecs, L_mat = exact_poles_from_liouvillian()
    
    # Group by unique frequency (merge conjugate pairs)
    seen = set()
    unique_poles = []
    for p in poles:
        f_round = round(p['frequency'], 3)
        d_round = round(p['decay_rate'], 3)
        key = (f_round, d_round)
        if key not in seen:
            seen.add(key)
            unique_poles.append(p)
    
    print(f"\n{'#':>3} | {'Frequency':>10} | {'Decay':>10} | {'Q-factor':>10} | {'Pole (s)':>30}")
    print("-" * 75)
    for i, p in enumerate(unique_poles[:15]):
        s = p['eigenvalue']
        print(f"{i+1:>3} | {p['frequency']:>10.4f} | {p['decay_rate']:>10.4f} | "
              f"{p['Q']:>10.1f} | {np.real(s):+.4f} {np.imag(s):+.4f}j")
    print(f"\nTotal unique poles with f>0.01: {len(unique_poles)}")

    # 2. Compute residues per channel
    print(f"\n--- RESIDUES: Which poles are bright in which channel ---")
    residues = compute_residues(evals, evecs, L_mat)
    
    print(f"\n{'#':>3} | {'Freq':>7} | {'Decay':>7} | {'|r_c+|':>8} {'ph_c+':>7} | {'|r_c-|':>8} {'ph_c-':>7} | Bright in")
    print("-" * 80)
    for i, r in enumerate(residues[:12]):
        bright = []
        if r['amp_cp'] > 0.001: bright.append('c+')
        if r['amp_cm'] > 0.001: bright.append('c-')
        bright_str = '+'.join(bright) if bright else 'dark'
        print(f"{i+1:>3} | {r['frequency']:>7.4f} | {r['decay_rate']:>7.4f} | "
              f"{r['amp_cp']:>8.5f} {math.degrees(r['phase_cp']):>+7.1f} | "
              f"{r['amp_cm']:>8.5f} {math.degrees(r['phase_cm']):>+7.1f} | {bright_str}")
    
    # Key question: do c+ and c- share poles or have different poles?
    print(f"\n--- KEY QUESTION: Shared poles or separate poles? ---")
    cp_poles = [(r['frequency'], r['decay_rate']) for r in residues if r['amp_cp'] > 0.005]
    cm_poles = [(r['frequency'], r['decay_rate']) for r in residues if r['amp_cm'] > 0.005]
    
    shared = []
    cp_only = []
    cm_only = []
    for f, d in cp_poles:
        if any(abs(f-f2)<0.02 and abs(d-d2)<0.02 for f2,d2 in cm_poles):
            shared.append((f,d))
        else:
            cp_only.append((f,d))
    for f, d in cm_poles:
        if not any(abs(f-f2)<0.02 and abs(d-d2)<0.02 for f2,d2 in cp_poles):
            cm_only.append((f,d))
    
    print(f"  Shared poles (bright in both): {len(shared)}")
    for f,d in shared: print(f"    f={f:.4f} decay={d:.4f}")
    print(f"  c+ only poles: {len(cp_only)}")
    for f,d in cp_only: print(f"    f={f:.4f} decay={d:.4f}")
    print(f"  c- only poles: {len(cm_only)}")
    for f,d in cm_only: print(f"    f={f:.4f} decay={d:.4f}")

    # 3. Pole trajectory across topology
    print(f"\n{'='*75}")
    print(f"POLE MAP: Tracking poles across J_SB sweep")
    print(f"If 'decay fixed, freq moves' is true, poles move horizontally in complex plane")
    print(f"{'='*75}")
    
    print(f"\n{'J_SB':>5} | Top poles (freq, decay) sorted by brightness")
    print("-" * 75)
    
    pole_tracks = {}
    for J_SB in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        poles_j, ev_j, evec_j, L_j = exact_poles_from_liouvillian(J_SB=J_SB)
        res_j = compute_residues(ev_j, evec_j, L_j, J_SB=J_SB)
        
        # Top 4 brightest poles
        top = res_j[:4]
        desc = "  ".join(f"f={r['frequency']:.3f}(d={r['decay_rate']:.3f})" for r in top)
        print(f"{J_SB:>5.1f} | {desc}")
        
        # Collect for tracking
        for r in res_j:
            f_bin = round(r['frequency'], 1)
            if f_bin not in pole_tracks:
                pole_tracks[f_bin] = []
            pole_tracks[f_bin].append({
                'J_SB': J_SB,
                'freq': r['frequency'],
                'decay': r['decay_rate'],
                'amp_cp': r['amp_cp'],
                'amp_cm': r['amp_cm'],
            })

    # 4. Decay rate stability check
    print(f"\n--- DECAY RATE STABILITY across J_SB ---")
    print(f"If decay is truly topology-independent, all decay values should cluster")
    
    all_decays = []
    for J_SB in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
        _, ev_j, evec_j, L_j = exact_poles_from_liouvillian(J_SB=J_SB)
        res_j = compute_residues(ev_j, evec_j, L_j, J_SB=J_SB)
        bright = [r for r in res_j if max(r['amp_cp'], r['amp_cm']) > 0.005]
        decays = sorted(set(round(r['decay_rate'], 4) for r in bright))
        all_decays.append((J_SB, decays))
        print(f"  J_SB={J_SB:.1f}: decay values = {decays}")
    
    # Collect all unique decay rates
    flat_decays = []
    for _, ds in all_decays:
        flat_decays.extend(ds)
    if flat_decays:
        print(f"\n  All decay rates: min={min(flat_decays):.4f} max={max(flat_decays):.4f}")
        print(f"  Spread: {max(flat_decays)-min(flat_decays):.4f}")
        if max(flat_decays) - min(flat_decays) < 0.05:
            print(f"  >>> DECAY RATES ARE CLUSTERED -- topology-independence CONFIRMED <<<")
        else:
            print(f"  >>> DECAY RATES SPREAD -- multiple distinct decay channels <<<")
    
    # 5. The GPT test: is the Q difference real?
    print(f"\n--- GPT's KEY QUESTION: Is the Q-factor difference REAL? ---")
    res_base = compute_residues(evals, evecs, L_mat)
    cp_dominant = max(res_base, key=lambda r: r['amp_cp'])
    cm_dominant = max(res_base, key=lambda r: r['amp_cm'])
    
    print(f"  Pole dominating c+: f={cp_dominant['frequency']:.4f} decay={cp_dominant['decay_rate']:.4f}")
    print(f"  Pole dominating c-: f={cm_dominant['frequency']:.4f} decay={cm_dominant['decay_rate']:.4f}")
    
    if abs(cp_dominant['decay_rate'] - cm_dominant['decay_rate']) > 0.01:
        print(f"  >>> DIFFERENT POLES with different decay rates <<<")
        print(f"  >>> GPT correct: sector-specific damping, not channel-specific <<<")
    else:
        same_pole = abs(cp_dominant['frequency'] - cm_dominant['frequency']) < 0.02
        if same_pole:
            print(f"  >>> SAME POLE dominates both channels <<<")
        else:
            print(f"  >>> Different frequencies but similar decay <<<")
    
    print(f"\n{'='*75}")
    print("JOINT POLE ANALYSIS COMPLETE")
    print("="*75)
