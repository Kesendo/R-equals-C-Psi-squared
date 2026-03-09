"""
Apply the window-XOR / skeleton analysis to IBM hardware data.
Same question: what's shared, what rotates?
"""
import json
import numpy as np

# Load hardware data
with open(r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\results\shadow_march\shadow_hardware_combined_20260309_181852.json") as f:
    data = json.load(f)

for qr in data['qubit_results']:
    qubit = qr['qubit']
    points = qr['points']
    
    print("=" * 65)
    print(f"QUBIT {qubit} - Skeleton Analysis")
    print("=" * 65)
    
    # Reconstruct density matrices
    rhos = []
    times = []
    for pt in points:
        rho = np.array([[pt['rho01_re'] * 2, pt['rho01_re'] + 1j * pt['rho01_im']],
                        [pt['rho01_re'] - 1j * pt['rho01_im'], 0]])
        # Actually, let's use the raw data properly
        # From the metrics: C = purity, rho01 = off-diagonal
        # For a 2x2: rho = [[rho00, rho01], [rho01*, rho11]]
        # We have rho01_re, rho01_im, and can get populations from C (purity)
        # C = rho00^2 + rho11^2 + 2|rho01|^2
        # But we also know rho00 + rho11 = 1
        # From CΨ data: cpsi, C (purity), psi
        # psi = 2*|rho01| for single qubit
        # So |rho01| = psi/2 = rho01_abs
        
        rho01 = pt['rho01_re'] + 1j * pt['rho01_im']
        rho01_abs = abs(rho01)
        C_purity = pt['C']
        
        # From purity: C = rho00^2 + rho11^2 + 2|rho01|^2
        # And rho00 + rho11 = 1, so rho00 = 1 - rho11
        # C = (1-rho11)^2 + rho11^2 + 2|rho01|^2
        # C = 1 - 2*rho11 + 2*rho11^2 + 2|rho01|^2
        # 2*rho11^2 - 2*rho11 + (1 + 2|rho01|^2 - C) = 0
        
        a_coeff = 2
        b_coeff = -2
        c_coeff = 1 + 2*rho01_abs**2 - C_purity
        disc = b_coeff**2 - 4*a_coeff*c_coeff
        if disc < 0: disc = 0
        rho11 = (-b_coeff - np.sqrt(disc)) / (2*a_coeff)  # smaller root = |1> population
        rho00 = 1 - rho11
        
        rho_mat = np.array([[rho00, rho01], [np.conj(rho01), rho11]], dtype=complex)
        rhos.append(rho_mat)
        times.append(pt['delay_us'])
    
    # Show all density matrices
    print(f"\n{'#':>2} {'t(us)':>7} | {'rho00':>7} {'rho11':>7} | {'Re01':>8} {'Im01':>8} {'|01|':>7} {'ph01':>7}")
    print("-" * 70)
    for i, (t, rho) in enumerate(zip(times, rhos)):
        r00 = np.real(rho[0,0])
        r11 = np.real(rho[1,1])
        re01 = np.real(rho[0,1])
        im01 = np.imag(rho[0,1])
        abs01 = abs(rho[0,1])
        ph01 = np.angle(rho[0,1]) / np.pi
        print(f"{i:>2} {t:>7.1f} | {r00:>7.4f} {r11:>7.4f} | {re01:>+8.4f} {im01:>+8.4f} {abs01:>7.4f} {ph01:>+6.3f}pi")

    # XOR analysis: overlay adjacent time points
    print(f"\nXOR ANALYSIS: What changes between adjacent measurements?")
    print(f"{'Pair':>6} | {'shared%':>8} | {'Re chg':>8} {'Im chg':>8} {'pop chg':>8} | {'ph_advance':>10}")
    print("-" * 70)
    
    for i in range(len(rhos) - 1):
        r_a = rhos[i]
        r_b = rhos[i+1]
        
        # Magnitudes
        mag_a = np.abs(r_a)
        mag_b = np.abs(r_b)
        shared = np.sum(np.minimum(mag_a, mag_b))
        total = max(np.sum(mag_a), np.sum(mag_b))
        shared_pct = shared / total * 100 if total > 0 else 0
        
        # What changed?
        d_re = np.real(r_b[0,1]) - np.real(r_a[0,1])
        d_im = np.imag(r_b[0,1]) - np.imag(r_a[0,1])
        d_pop = np.real(r_b[0,0]) - np.real(r_a[0,0])
        
        # Phase advance
        ph_a = np.angle(r_a[0,1])
        ph_b = np.angle(r_b[0,1])
        d_ph = (ph_b - ph_a) / np.pi
        while d_ph > 1: d_ph -= 2
        while d_ph < -1: d_ph += 2
        
        print(f" {i}->{i+1} | {shared_pct:>7.1f}% | {d_re:>+8.4f} {d_im:>+8.4f} {d_pop:>+8.4f} | {d_ph:>+9.3f}pi")

    # The key question: separate skeleton from rotation
    print(f"\nSKELETON vs ROTATION decomposition:")
    print(f"  Skeleton = populations (diagonal)")
    print(f"  Rotation = phase of off-diagonal")
    
    # Skeleton: how much do populations change?
    pop_changes = []
    phase_changes = []
    for i in range(len(rhos) - 1):
        dp = abs(np.real(rhos[i+1][0,0]) - np.real(rhos[i][0,0]))
        pop_changes.append(dp)
        
        ph_a = np.angle(rhos[i][0,1])
        ph_b = np.angle(rhos[i+1][0,1])
        dph = abs(ph_b - ph_a)
        if dph > np.pi: dph = 2*np.pi - dph
        phase_changes.append(dph)
    
    avg_pop = np.mean(pop_changes)
    avg_phase = np.mean(phase_changes)
    
    print(f"\n  Avg population change per step:  {avg_pop:.4f}")
    print(f"  Avg phase change per step:       {avg_phase/np.pi:.3f} pi")
    print(f"  Ratio (phase/population):        {avg_phase/np.pi / (avg_pop + 1e-10):.1f}x")
    
    if avg_phase/np.pi > 3 * avg_pop:
        print(f"  -> PHASE DOMINATES: the skeleton is stable, what rotates is the phase")
    elif avg_pop > 3 * avg_phase/np.pi:
        print(f"  -> POPULATIONS DOMINATE: the skeleton is changing, phase is secondary")
    else:
        print(f"  -> MIXED: both populations and phase change comparably")

    # Compare to Lindblad prediction
    print(f"\n  Comparison to Lindblad (T1/T2 model):")
    T1 = qr['verdict']['T1_us']
    T2 = qr['verdict']['T2_us']
    print(f"  T1 = {T1:.1f} us, T2 = {T2:.1f} us")
    
    for i, (t, rho) in enumerate(zip(times, rhos)):
        if t == 0: continue
        # Lindblad prediction
        rho11_theory = 0.5 * np.exp(-t / T1)
        rho00_theory = 1 - rho11_theory
        rho01_theory = 0.5 * np.exp(-t / T2)  # real, no phase rotation
        
        # Actual
        rho01_actual = rho[0,1]
        phase_actual = np.angle(rho01_actual) / np.pi
        
        # The Lindblad model predicts rho01 stays REAL (phase = 0)
        # Any phase != 0 is structure the model doesn't predict
        print(f"  t={t:>6.1f}: Lindblad phase=0.000pi, actual phase={phase_actual:>+.3f}pi"
              f"  excess_phase={phase_actual:>+.3f}pi")

print(f"\n{'=' * 65}")
print("DONE")
print("=" * 65)
