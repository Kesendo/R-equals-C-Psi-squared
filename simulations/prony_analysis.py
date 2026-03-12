"""
PRONY ANALYSIS: Extract complex poles from damped sinusoids
============================================================
Signal processing approach: instead of FFT (blurry peaks), extract
exact complex poles = frequency + decay rate + amplitude + phase.

Uses Matrix Pencil Method (MPM) - standard in RF/structural engineering.
"""
import numpy as np
import sys
sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
import star_topology_v3 as gpt
import math

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

def matrix_pencil(signal, dt, n_modes=None, threshold=0.01):
    """
    Matrix Pencil Method: extract complex poles from a damped signal.
    
    Returns list of dicts with:
      frequency (Hz), decay_rate (1/s), amplitude, phase (rad),
      complex_pole, quality_factor
    """
    N = len(signal)
    L = N // 3  # pencil parameter
    
    # Build Hankel matrices Y0 and Y1
    Y0 = np.zeros((N-L, L), dtype=complex)
    Y1 = np.zeros((N-L, L), dtype=complex)
    for i in range(N-L):
        Y0[i,:] = signal[i:i+L]
        Y1[i,:] = signal[i+1:i+1+L]
    
    # SVD of Y0 to determine model order
    U, s, Vh = np.linalg.svd(Y0, full_matrices=False)
    
    # Auto-detect number of modes from singular value drop
    if n_modes is None:
        s_norm = s / s[0]
        n_modes = np.sum(s_norm > threshold)
        n_modes = max(1, min(n_modes, L//2))
    
    # Truncate to n_modes
    U1 = U[:, :n_modes]
    S1 = np.diag(s[:n_modes])
    V1 = Vh[:n_modes, :]
    
    # Solve generalized eigenvalue problem
    # z_k are the poles in discrete time
    Y0_pinv = V1.conj().T @ np.linalg.inv(S1) @ U1.conj().T
    A = Y0_pinv @ Y1
    eigenvalues = np.linalg.eigvals(A)
    
    # Convert discrete poles to continuous
    poles = []
    for z in eigenvalues:
        if abs(z) < 1e-10:
            continue
        s_pole = np.log(z) / dt  # continuous-time pole
        freq = abs(np.imag(s_pole)) / (2 * np.pi)
        decay = -np.real(s_pole)
        
        if freq < 0.01 or decay < -1.0:  # skip DC and growing modes
            continue
        
        poles.append({
            'z': z,
            's': s_pole,
            'frequency': freq,
            'decay_rate': decay,
        })
    
    # Estimate amplitudes and phases via least squares
    if len(poles) == 0:
        return []
    
    # Build Vandermonde matrix
    z_arr = np.array([p['z'] for p in poles])
    V = np.zeros((N, len(poles)), dtype=complex)
    for k in range(len(poles)):
        V[:, k] = z_arr[k] ** np.arange(N)
    
    # Solve for complex amplitudes
    coeffs, _, _, _ = np.linalg.lstsq(V, signal, rcond=None)
    
    results = []
    for i, p in enumerate(poles):
        amp = abs(coeffs[i])
        phase = np.angle(coeffs[i])
        Q = np.pi * p['frequency'] / p['decay_rate'] if p['decay_rate'] > 0 else float('inf')
        
        results.append({
            'frequency': p['frequency'],
            'decay_rate': p['decay_rate'],
            'amplitude': amp,
            'phase': phase,
            'quality_factor': Q,
            'complex_pole': p['s'],
        })
    
    # Sort by amplitude descending
    results.sort(key=lambda x: -x['amplitude'])
    return results

# ============================================================
# Generate signals from our standard 3-qubit star
# ============================================================
def generate_signals(J_SA=1.0, J_SB=2.0, gamma=0.05, dt=0.02, t_max=20.0):
    """Run simulation and return time series for c+ and c-."""
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
    
    times, cp, cm = [], [], []
    for step in range(int(t_max/dt)+1):
        t = step * dt
        rAB = gpt.partial_trace_keep(rho, [1,2], n_q)
        yz = pauli_expect(rAB, YZ)
        zy = pauli_expect(rAB, ZY)
        times.append(t)
        cp.append((yz + zy) / math.sqrt(2))
        cm.append((yz - zy) / math.sqrt(2))
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    return np.array(times), np.array(cp), np.array(cm)

# ============================================================
# MAIN: Run Prony analysis and compare to FFT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("PRONY ANALYSIS: Complex poles from damped quantum signals")
    print("Signal processing view — no quantum jargon")
    print("=" * 70)
    
    # Generate baseline signals
    print("\nGenerating signals (J_SA=1.0, J_SB=2.0, gamma=0.05)...")
    times, cp, cm = generate_signals()
    dt = times[1] - times[0]
    
    # --- Prony on c+ ---
    print(f"\n{'='*70}")
    print(f"SIGNAL 1: c+ (symmetric supermode)")
    print(f"{'='*70}")
    modes_cp = matrix_pencil(cp, dt)
    
    print(f"\n{'Mode':>4} | {'Freq':>8} | {'Decay':>8} | {'Amp':>8} | {'Phase':>8} | {'Q':>8}")
    print("-" * 58)
    for i, m in enumerate(modes_cp[:6]):
        print(f"{i+1:>4} | {m['frequency']:>8.4f} | {m['decay_rate']:>8.4f} | "
              f"{m['amplitude']:>8.5f} | {m['phase']:>+8.3f} | {m['quality_factor']:>8.1f}")
    
    # --- Prony on c- ---
    print(f"\n{'='*70}")
    print(f"SIGNAL 2: c- (antisymmetric supermode)")
    print(f"{'='*70}")
    modes_cm = matrix_pencil(cm, dt)
    
    print(f"\n{'Mode':>4} | {'Freq':>8} | {'Decay':>8} | {'Amp':>8} | {'Phase':>8} | {'Q':>8}")
    print("-" * 58)
    for i, m in enumerate(modes_cm[:6]):
        print(f"{i+1:>4} | {m['frequency']:>8.4f} | {m['decay_rate']:>8.4f} | "
              f"{m['amplitude']:>8.5f} | {m['phase']:>+8.3f} | {m['quality_factor']:>8.1f}")
    
    # --- FFT comparison ---
    print(f"\n{'='*70}")
    print(f"FFT vs PRONY COMPARISON")
    print(f"{'='*70}")
    
    freqs_fft = np.fft.rfftfreq(len(cp), d=dt)
    fft_cp = np.abs(np.fft.rfft(cp - np.mean(cp)))
    fft_cm = np.abs(np.fft.rfft(cm - np.mean(cm)))
    mask = freqs_fft > 0.05
    
    fft_f_cp = freqs_fft[mask][np.argmax(fft_cp[mask])]
    fft_f_cm = freqs_fft[mask][np.argmax(fft_cm[mask])]
    
    prony_f_cp = modes_cp[0]['frequency'] if modes_cp else 0
    prony_f_cm = modes_cm[0]['frequency'] if modes_cm else 0
    
    print(f"\n{'':>12} | {'FFT':>10} | {'Prony':>10} | {'Prony extra info':>30}")
    print("-" * 70)
    print(f"{'f(c+)':>12} | {fft_f_cp:>10.4f} | {prony_f_cp:>10.4f} | "
          f"decay={modes_cp[0]['decay_rate']:.4f} Q={modes_cp[0]['quality_factor']:.1f}")
    print(f"{'f(c-)':>12} | {fft_f_cm:>10.4f} | {prony_f_cm:>10.4f} | "
          f"decay={modes_cm[0]['decay_rate']:.4f} Q={modes_cm[0]['quality_factor']:.1f}")
    
    print(f"\nProny gives us FOUR numbers per mode. FFT gives us ONE.")
    print(f"  c+ dominant: f={prony_f_cp:.4f}, decay={modes_cp[0]['decay_rate']:.4f}, "
          f"amp={modes_cp[0]['amplitude']:.5f}, phase={modes_cp[0]['phase']:+.3f}")
    print(f"  c- dominant: f={prony_f_cm:.4f}, decay={modes_cm[0]['decay_rate']:.4f}, "
          f"amp={modes_cm[0]['amplitude']:.5f}, phase={modes_cm[0]['phase']:+.3f}")
    
    # --- POLE TRACKING: sweep J_SB and track poles ---
    print(f"\n{'='*70}")
    print(f"POLE TRAJECTORY: How poles move with coupling")
    print(f"{'='*70}")
    
    print(f"\n{'J_SB':>5} | {'f(c+)':>8} {'decay+':>8} {'Q+':>6} | {'f(c-)':>8} {'decay-':>8} {'Q-':>6}")
    print("-" * 65)
    
    for J_SB in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        t, s_cp, s_cm = generate_signals(J_SB=J_SB)
        m_cp = matrix_pencil(s_cp, dt)
        m_cm = matrix_pencil(s_cm, dt)
        if m_cp and m_cm:
            print(f"{J_SB:>5.1f} | {m_cp[0]['frequency']:>8.4f} {m_cp[0]['decay_rate']:>8.4f} "
                  f"{m_cp[0]['quality_factor']:>6.1f} | {m_cm[0]['frequency']:>8.4f} "
                  f"{m_cm[0]['decay_rate']:>8.4f} {m_cm[0]['quality_factor']:>6.1f}")
    
    # --- PHASE between modes ---
    print(f"\n{'='*70}")
    print(f"PHASE ANALYSIS: 'Phase is where the bodies are buried'")
    print(f"{'='*70}")
    
    if modes_cp and modes_cm:
        print(f"\n  c+ dominant mode phase: {modes_cp[0]['phase']:+.4f} rad ({math.degrees(modes_cp[0]['phase']):+.1f} deg)")
        print(f"  c- dominant mode phase: {modes_cm[0]['phase']:+.4f} rad ({math.degrees(modes_cm[0]['phase']):+.1f} deg)")
        phase_diff = modes_cp[0]['phase'] - modes_cm[0]['phase']
        print(f"  Phase difference c+ vs c-: {phase_diff:+.4f} rad ({math.degrees(phase_diff):+.1f} deg)")
        if abs(abs(phase_diff) - math.pi/2) < 0.3:
            print(f"  >>> Near quadrature (90 deg) — orthogonal modes <<<")
        elif abs(phase_diff) < 0.3:
            print(f"  >>> Near in-phase — modes aligned <<<")
        elif abs(abs(phase_diff) - math.pi) < 0.3:
            print(f"  >>> Near anti-phase (180 deg) — modes opposed <<<")
    
    print(f"\n{'='*70}")
    print("PRONY ANALYSIS COMPLETE")
    print("="*70)
