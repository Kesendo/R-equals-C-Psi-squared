"""
CORRELATED BATH SWEEP: Can bath geometry select which sector is visible?
=========================================================================
GPT's experiment suggestion (March 12, 2026):

All our noise tests used LOCAL dephasing (each qubit independent).
What happens with CORRELATED noise where A and B share a common bath?

Two knobs:
  eta in [0,1]: local -> fully correlated bath
  phi in [0,pi/2]: common ZZ bath -> common XX bath

L_corr(phi) = sqrt(eta*gamma) * (cos(phi)*Z + sin(phi)*X) on A+B combined
L_local = sqrt((1-eta)*gamma) * Z on each qubit individually

Key question: Can the bath geometry itself choose which sector is visible?
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

def build_jump_ops(eta, phi, gamma, n_q=3):
    """
    Build jump operators for correlated + local noise.
    n_q=3: S(0), A(1), B(2)
    
    Local part: sqrt((1-eta)*gamma) * Z on each qubit
    Correlated part: sqrt(eta*gamma) * (cos(phi)*Z + sin(phi)*X) on A+B combined
    S always gets local dephasing at full gamma.
    """
    L_ops = []
    
    # S always gets local dephasing (mediator noise unchanged)
    L_ops.append(np.sqrt(gamma) * op_at(Z, 0, n_q))
    
    # Local noise on A and B (scaled by 1-eta)
    if eta < 1.0:
        L_ops.append(np.sqrt((1-eta)*gamma) * op_at(Z, 1, n_q))
        L_ops.append(np.sqrt((1-eta)*gamma) * op_at(Z, 2, n_q))
    
    # Correlated noise on A+B (scaled by eta)
    if eta > 0.0:
        # Correlated operator: (cos(phi)*Z + sin(phi)*X) applied to BOTH A and B
        op_A = math.cos(phi) * op_at(Z, 1, n_q) + math.sin(phi) * op_at(X, 1, n_q)
        op_B = math.cos(phi) * op_at(Z, 2, n_q) + math.sin(phi) * op_at(X, 2, n_q)
        L_corr = np.sqrt(eta * gamma / 2) * (op_A + op_B)
        L_ops.append(L_corr)
    
    return L_ops

def measure_structure(H, L_ops, rho0, dt=0.005, t_max=20.0):
    """Measure all structural indicators."""
    YZ = np.kron(Y, Z)
    ZY = np.kron(Z, Y)
    XX = np.kron(X, X)
    
    times, cp_vals, cm_vals = [], [], []
    rho = rho0.copy()
    
    for step in range(int(t_max/dt)+1):
        t = step * dt
        if step % 4 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], 3)
            yz = pauli_expect(rAB, YZ)
            zy = pauli_expect(rAB, ZY)
            times.append(t)
            cp_vals.append((yz + zy) / math.sqrt(2))
            cm_vals.append((yz - zy) / math.sqrt(2))
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    # Final state measurements
    rAB_final = gpt.partial_trace_keep(rho, [1,2], 3)
    xx_comm = np.linalg.norm(rAB_final @ XX - XX @ rAB_final)
    
    # Skeleton/rotation: compare to time-average
    times = np.array(times)
    dt_s = times[1] - times[0]
    
    # FFT for frequencies
    cp = np.array(cp_vals) - np.mean(cp_vals)
    cm = np.array(cm_vals) - np.mean(cm_vals)
    freqs = np.fft.rfftfreq(len(cp), d=dt_s)
    fft_cp = np.abs(np.fft.rfft(cp))
    fft_cm = np.abs(np.fft.rfft(cm))
    
    mask = freqs > 0.05
    f_cp = freqs[mask][np.argmax(fft_cp[mask])] if np.any(fft_cp[mask] > 0.01) else 0
    f_cm = freqs[mask][np.argmax(fft_cm[mask])] if np.any(fft_cm[mask] > 0.01) else 0
    a_cp = np.max(fft_cp[mask]) if np.any(mask) else 0
    a_cm = np.max(fft_cm[mask]) if np.any(mask) else 0
    
    return {
        'f_cp': f_cp, 'f_cm': f_cm,
        'a_cp': a_cp, 'a_cm': a_cm,
        'xx_comm': xx_comm,
        'ratio': a_cp/a_cm if a_cm > 0.01 else float('inf'),
    }

# Setup: standard 3-qubit star
J_SA, J_SB, gamma = 1.0, 2.0, 0.05

def star_H_3q():
    H = np.zeros((8,8), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,3) @ op_at(p,1,3)
        H += J_SB * op_at(p,0,3) @ op_at(p,2,3)
    return H

H = star_H_3q()
psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho0 = np.outer(psi0, psi0.conj())

# Baseline: pure local noise
L_local = [np.sqrt(gamma)*op_at(Z,i,3) for i in range(3)]
baseline = measure_structure(H, L_local, rho0)

print("="*70)
print("CORRELATED BATH SWEEP")
print("Can bath geometry select which sector is visible?")
print("="*70)
print(f"\nBaseline (local only): f(c+)={baseline['f_cp']:.3f} f(c-)={baseline['f_cm']:.3f}")
print(f"  amp ratio={baseline['ratio']:.2f} XX={baseline['xx_comm']:.1e}")

# ============================================================
# SWEEP 1: eta sweep at fixed phi=0 (correlated ZZ bath)
# ============================================================
print(f"\n--- SWEEP 1: eta (local -> correlated), phi=0 (ZZ bath) ---")
print(f"{'eta':>5} | {'f(c+)':>7} {'f(c-)':>7} | {'A+/A-':>6} | {'XX':>8} | Note")
print("-"*60)

for eta in [0.0, 0.1, 0.2, 0.5, 0.8, 1.0]:
    L = build_jump_ops(eta, 0.0, gamma)
    r = measure_structure(H, L, rho0)
    note = ""
    if abs(r['f_cp'] - baseline['f_cp']) > 0.02: note += "FREQ SHIFT! "
    if r['xx_comm'] > 1e-6: note += "XX BROKEN! "
    if abs(r['ratio'] - baseline['ratio']) > 0.5: note += "RATIO SHIFT! "
    print(f"{eta:>5.1f} | {r['f_cp']:>7.3f} {r['f_cm']:>7.3f} | {r['ratio']:>6.2f} | {r['xx_comm']:>8.1e} | {note}")

# ============================================================
# SWEEP 2: phi sweep at fixed eta=1.0 (fully correlated, vary bath type)
# ============================================================
print(f"\n--- SWEEP 2: phi (ZZ bath -> XX bath), eta=1.0 (fully correlated) ---")
print(f"{'phi/pi':>6} | {'f(c+)':>7} {'f(c-)':>7} | {'A+/A-':>6} | {'XX':>8} | Bath type")
print("-"*65)

for phi_frac in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    phi = phi_frac * math.pi / 2
    L = build_jump_ops(1.0, phi, gamma)
    r = measure_structure(H, L, rho0)
    bath = "ZZ" if phi_frac < 0.1 else "XX" if phi_frac > 0.9 else "mixed"
    note = ""
    if abs(r['f_cp'] - baseline['f_cp']) > 0.02: note += "FREQ! "
    if r['xx_comm'] > 1e-6: note += "XX BROKEN! "
    print(f"{phi_frac:>6.1f} | {r['f_cp']:>7.3f} {r['f_cm']:>7.3f} | {r['ratio']:>6.2f} | {r['xx_comm']:>8.1e} | {bath} {note}")

# ============================================================
# SWEEP 3: Full 2D grid (coarse)
# ============================================================
print(f"\n--- SWEEP 3: 2D grid (eta x phi) ---")
print(f"Values = amplitude ratio A+/A- (baseline={baseline['ratio']:.2f})")
print(f"{'':>6}", end="")
phi_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
for pf in phi_vals:
    print(f" | phi={pf:.2f}", end="")
print()
print("-" * (8 + 10*len(phi_vals)))

eta_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
xx_broken = []

for eta in eta_vals:
    print(f"e={eta:.2f}", end="")
    for pf in phi_vals:
        phi = pf * math.pi / 2
        L = build_jump_ops(eta, phi, gamma)
        r = measure_structure(H, L, rho0)
        marker = "*" if r['xx_comm'] > 1e-6 else " "
        print(f" | {r['ratio']:>5.2f}{marker}  ", end="")
        if r['xx_comm'] > 1e-6:
            xx_broken.append((eta, pf, r['xx_comm']))
    print()

print(f"\n* = XX symmetry broken")
if xx_broken:
    print(f"XX broken at: {[(f'eta={e:.2f},phi={p:.2f}',f'{x:.1e}') for e,p,x in xx_broken]}")
else:
    print(f"XX survived EVERYWHERE")

# ============================================================
# VERDICT
# ============================================================
print(f"\n{'='*70}")
print("VERDICT")
print("="*70)
print("""
Three possible outcomes (GPT's prediction):

1. Only amplitudes move, frequencies stay fixed
   -> "Noise damps" survives, upgraded: bath geometry selects visibility

2. Frequencies or mode identities move
   -> "Noise is irrelevant to structure" is too strong, locality did the work

3. Nothing changes at all
   -> Robustness claim becomes much harder to dismiss
""")
