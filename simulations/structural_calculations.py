"""
STRUCTURAL CALCULATIONS - Reconstructed from GPT's March 2026 analysis
Verifies: X⊗X symmetry, Bell-sector structure, S-coherence gating
"""
import numpy as np
import sys
sys.path.insert(0, '.')
import gpt_code as gpt

def bell_fidelities(rho_AB):
    bells = {
        "Phi+": np.array([1,0,0,1], dtype=complex)/np.sqrt(2),
        "Phi-": np.array([1,0,0,-1], dtype=complex)/np.sqrt(2),
        "Psi+": np.array([0,1,1,0], dtype=complex)/np.sqrt(2),
        "Psi-": np.array([0,1,-1,0], dtype=complex)/np.sqrt(2),
    }
    return {n: float(np.real(np.trace(rho_AB @ np.outer(b, b.conj())))) 
            for n, b in bells.items()}

def check_xx_symmetry(rho_AB):
    """Check if rho_AB commutes with X⊗X."""
    sx = np.array([[0,1],[1,0]], dtype=complex)
    XX = np.kron(sx, sx)
    comm = rho_AB @ XX - XX @ rho_AB
    return np.linalg.norm(comm)


# Setup
H = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
L_ops = gpt.dephasing_ops_n([0.05]*3)
psi = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho = gpt.density_from_statevector(psi)

dt = 0.005
peaks = []
prev_cpsi = 0
rising = True

# Also run unitary reference
rho_unitary = rho.copy()
L_ops_zero = gpt.dephasing_ops_n([0.0]*3)

all_cpsi = []
all_s_coh = []
max_comm = 0.0

for step in range(2001):
    t = step * dt
    if step % 4 == 0:
        rho_AB = gpt.partial_trace_keep(rho, keep=[1,2], n_qubits=3)
        rho_S = gpt.partial_trace_keep(rho, keep=[0], n_qubits=3)
        c = gpt.concurrence_two_qubit(rho_AB)
        p = gpt.psi_norm(rho_AB)
        cpsi = c * p
        s_coh = gpt.l1_coherence(rho_S)
        
        # Check X⊗X symmetry
        comm_norm = check_xx_symmetry(rho_AB)
        max_comm = max(max_comm, comm_norm)
        
        all_cpsi.append(cpsi)
        all_s_coh.append(s_coh)
        
        # Detect peaks
        if prev_cpsi > 0.05 and cpsi < prev_cpsi and rising:
            # Get unitary reference at same time
            rho_AB_u = gpt.partial_trace_keep(rho_unitary, keep=[1,2], n_qubits=3)
            peaks.append({
                't': t - dt*4,
                'cpsi': prev_cpsi,
                'rho_AB': gpt.partial_trace_keep(rho, keep=[1,2], n_qubits=3),
                'rho_AB_unitary': rho_AB_u,
            })
            rising = False
        if cpsi > prev_cpsi:
            rising = True
        prev_cpsi = cpsi
    
    if step < 2000:
        rho = gpt.rk4_step(rho, H, L_ops, dt)
        rho_unitary = gpt.rk4_step(rho_unitary, H, L_ops_zero, dt)


# Results
print("=" * 70)
print("STRUCTURAL CALCULATIONS")
print("=" * 70)

print(f"\n1. X⊗X SYMMETRY: max commutator norm = {max_comm:.2e}")
if max_comm < 1e-10:
    print("   -> EXACT symmetry confirmed. rho_AB block-diagonal in Bell basis.")
    print("   -> Two parity sectors: {Phi+, Psi+} and {Phi-, Psi-}")
    print("   -> 7-parameter family, not generic 15-parameter")

print(f"\n2. WINDOW METRICS (first {min(8, len(peaks))} peaks):")
print(f"{'#':>2} {'t':>5} {'CPsi':>6} | {'Phi+':>5} {'Phi-':>5} {'Psi+':>5} {'Psi-':>5} | {'H_Bell':>6} | {'dom':>4}")
print("-" * 65)

for i, pk in enumerate(peaks[:8]):
    fids = bell_fidelities(pk['rho_AB'])
    F = [fids['Phi+'], fids['Phi-'], fids['Psi+'], fids['Psi-']]
    H_bell = -sum(f * np.log2(f) if f > 0.001 else 0 for f in F)
    dom = max(fids, key=fids.get)
    print(f"{i:>2} {pk['t']:>5.2f} {pk['cpsi']:>6.3f} | "
          f"{F[0]:>5.3f} {F[1]:>5.3f} {F[2]:>5.3f} {F[3]:>5.3f} | "
          f"{H_bell:>6.3f} | {dom:>4}")

# Population symmetry
print(f"\n3. POPULATION SYMMETRY:")
for i, pk in enumerate(peaks[:3]):
    r = pk['rho_AB']
    d = np.real(np.diag(r))
    print(f"   Window {i}: |rho00-rho11| = {abs(d[0]-d[3]):.2e}, "
          f"|rho01-rho10| = {abs(d[1]-d[2]):.2e}")

# S-coherence correlation
from scipy.stats import pearsonr
cpsi_arr = np.array(all_cpsi)
s_coh_arr = np.array(all_s_coh)
mask = (cpsi_arr > 0) | (s_coh_arr > 0)
r_sc, p_sc = pearsonr(cpsi_arr[mask], s_coh_arr[mask])
print(f"\n4. S-COHERENCE GATING:")
print(f"   Pearson(CΨ_AB, S_l1_coh) = {r_sc:.3f} (p = {p_sc:.2e})")

# Noisy vs unitary comparison
print(f"\n5. NOISY vs UNITARY at peaks:")
for i, pk in enumerate(peaks[:4]):
    r_noisy = pk['rho_AB']
    r_unit = pk['rho_AB_unitary']
    td = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(r_noisy - r_unit)))
    pur_loss = float(np.real(np.trace(r_unit @ r_unit) - np.trace(r_noisy @ r_noisy)))
    phase_n = np.angle(r_noisy[0,3]) / np.pi
    phase_u = np.angle(r_unit[0,3]) / np.pi
    print(f"   Window {i}: TD={td:.3f}, pur_loss={pur_loss:.3f}, "
          f"phase(noisy)={phase_n:+.3f}pi, phase(unit)={phase_u:+.3f}pi")

print(f"\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)
