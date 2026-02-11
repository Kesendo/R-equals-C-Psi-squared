"""
Bridge Fingerprint Experiment — Reconstructed Simulation

Two coupled 2-qubit systems (A = receiver, B = sender) connected through a
Heisenberg bridge. System A starts classical (|00⟩, C·Ψ = 0). System B starts
in various quantum states. A's C·Ψ trajectory constitutes a unique "fingerprint"
of B's initial state.

Original simulation: 2026-02-09 (lost)
Reconstruction: 2026-02-11 (Guardian review session)

Physical setup:
  4 qubits: A = {0,1}, B = {2,3}
  H = H_A + H_B + H_bridge
  H_A: Heisenberg 0↔1, J = 1.0
  H_B: Heisenberg 2↔3, J = 1.0
  H_bridge: Heisenberg 1↔2, J = variable
  Decoherence: Local dephasing σ_z on all 4 qubits, γ = 0.1
  Evolution: First-order Lindblad, dt = 0.005, t_max = 5.0

Usage:
  python bridge_fingerprints.py
"""

import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import json

# ============================================================
# Pauli matrices and tensor products
# ============================================================

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_n(*matrices):
    """Kronecker product of N matrices."""
    result = matrices[0]
    for m in matrices[1:]:
        result = np.kron(result, m)
    return result


def pauli_on_qubit(pauli, qubit, n_qubits=4):
    """Place a Pauli operator on a specific qubit in an n-qubit system."""
    ops = [I2] * n_qubits
    ops[qubit] = pauli
    return kron_n(*ops)


# ============================================================
# Hamiltonian construction
# ============================================================

def heisenberg_coupling(q1, q2, J, n_qubits=4):
    """Heisenberg coupling H = J(σx⊗σx + σy⊗σy + σz⊗σz) between qubits q1, q2."""
    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for pauli in [sx, sy, sz]:
        ops = [I2] * n_qubits
        ops[q1] = pauli
        ops[q2] = pauli
        H += J * kron_n(*ops)
    return H


def build_hamiltonian(J_internal=1.0, J_bridge=0.5, n_qubits=4):
    """Build full Hamiltonian: H_A(0↔1) + H_B(2↔3) + H_bridge(1↔2)."""
    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    H += heisenberg_coupling(0, 1, J_internal, n_qubits)  # H_A
    H += heisenberg_coupling(2, 3, J_internal, n_qubits)  # H_B
    H += heisenberg_coupling(1, 2, J_bridge, n_qubits)    # H_bridge
    return H


# ============================================================
# Lindblad operators (local dephasing)
# ============================================================

def build_dephasing_operators(gamma, n_qubits=4):
    """Local dephasing: L_k = √γ σ_z on each qubit."""
    ops = []
    for q in range(n_qubits):
        L = np.sqrt(gamma) * pauli_on_qubit(sz, q, n_qubits)
        ops.append(L)
    return ops


# ============================================================
# Lindblad evolution
# ============================================================

def lindblad_rhs(rho, H, L_ops):
    """Compute dρ/dt = -i[H,ρ] + Σ_k (L_k ρ L_k† - ½{L_k†L_k, ρ})."""
    drho = -1j * (H @ rho - rho @ H)
    for L in L_ops:
        Ld = L.conj().T
        LdL = Ld @ L
        drho += L @ rho @ Ld - 0.5 * (LdL @ rho + rho @ LdL)
    return drho


def evolve_lindblad(rho0, H, L_ops, dt, n_steps):
    """First-order Euler Lindblad evolution. Returns density matrix at each step."""
    dim = rho0.shape[0]
    rho = rho0.copy()
    trajectory = [rho.copy()]

    for _ in range(n_steps):
        drho = lindblad_rhs(rho, H, L_ops)
        rho = rho + dt * drho
        # Enforce Hermiticity (numerical drift)
        rho = 0.5 * (rho + rho.conj().T)
        trajectory.append(rho.copy())

    return trajectory


# ============================================================
# Partial trace
# ============================================================

def partial_trace_B(rho, d_A=4, d_B=4):
    """Trace out system B, return ρ_A. Assumes ρ is (d_A*d_B × d_A*d_B)."""
    rho_A = np.zeros((d_A, d_A), dtype=complex)
    for j in range(d_B):
        proj = np.zeros((d_B, 1), dtype=complex)
        proj[j, 0] = 1.0
        bra = np.kron(np.eye(d_A, dtype=complex), proj.conj().T)
        ket = np.kron(np.eye(d_A, dtype=complex), proj)
        rho_A += bra @ rho @ ket
    return rho_A


def partial_trace_A(rho, d_A=4, d_B=4):
    """Trace out system A, return ρ_B."""
    rho_B = np.zeros((d_B, d_B), dtype=complex)
    for i in range(d_A):
        proj = np.zeros((d_A, 1), dtype=complex)
        proj[i, 0] = 1.0
        bra = np.kron(proj.conj().T, np.eye(d_B, dtype=complex))
        ket = np.kron(proj, np.eye(d_B, dtype=complex))
        rho_B += bra @ rho @ ket
    return rho_B


# ============================================================
# C·Ψ metrics
# ============================================================

def purity(rho):
    """C = Tr(ρ²)."""
    return np.real(np.trace(rho @ rho))


def l1_coherence(rho):
    """L1 coherence = Σ_{i≠j} |ρ_{ij}|."""
    d = rho.shape[0]
    total = 0.0
    for i in range(d):
        for j in range(d):
            if i != j:
                total += abs(rho[i, j])
    return total


def cpsi(rho):
    """C·Ψ where C = purity, Ψ = L1/(d-1)."""
    d = rho.shape[0]
    C = purity(rho)
    Psi = l1_coherence(rho) / (d - 1)
    return C * Psi


# ============================================================
# Initial states
# ============================================================

def ket(bits, n_qubits=4):
    """Computational basis state from bit string, e.g. '0110'."""
    dim = 2**n_qubits
    index = int(bits, 2)
    v = np.zeros((dim, 1), dtype=complex)
    v[index, 0] = 1.0
    return v


def plus_state():
    """Single-qubit |+⟩ = (|0⟩+|1⟩)/√2."""
    return np.array([[1], [1]], dtype=complex) / np.sqrt(2)


def minus_state():
    """Single-qubit |−⟩ = (|0⟩−|1⟩)/√2."""
    return np.array([[1], [-1]], dtype=complex) / np.sqrt(2)


def zero_state():
    """Single-qubit |0⟩."""
    return np.array([[1], [0]], dtype=complex)


def one_state():
    """Single-qubit |1⟩."""
    return np.array([[0], [1]], dtype=complex)


def bell_plus_2q():
    """2-qubit Bell+ = (|00⟩+|11⟩)/√2."""
    v = np.zeros((4, 1), dtype=complex)
    v[0, 0] = 1 / np.sqrt(2)  # |00⟩
    v[3, 0] = 1 / np.sqrt(2)  # |11⟩
    return v


def bell_minus_2q():
    """2-qubit Bell- = (|00⟩-|11⟩)/√2."""
    v = np.zeros((4, 1), dtype=complex)
    v[0, 0] = 1 / np.sqrt(2)
    v[3, 0] = -1 / np.sqrt(2)
    return v


def psi_plus_2q():
    """2-qubit |Ψ+⟩ = (|01⟩+|10⟩)/√2."""
    v = np.zeros((4, 1), dtype=complex)
    v[1, 0] = 1 / np.sqrt(2)  # |01⟩
    v[2, 0] = 1 / np.sqrt(2)  # |10⟩
    return v


def make_initial_state(label):
    """
    Build 4-qubit initial state: A = |00⟩ always, B = specified state.
    Returns 16×16 density matrix.
    """
    # System A always starts as |00⟩ (classical)
    psi_A = np.kron(zero_state(), zero_state())  # 4×1

    if label == '|++>':
        psi_B = np.kron(plus_state(), plus_state())
    elif label == '|+0>':
        psi_B = np.kron(plus_state(), zero_state())
    elif label == 'Bell+':
        psi_B = bell_plus_2q()
    elif label == 'Bell-':
        psi_B = bell_minus_2q()
    elif label == '|Psi+>':
        psi_B = psi_plus_2q()
    elif label == '|+->':
        psi_B = np.kron(plus_state(), minus_state())
    elif label == '|01>':
        psi_B = np.kron(zero_state(), one_state())
    elif label == '|11>':
        psi_B = np.kron(one_state(), one_state())
    else:
        raise ValueError(f"Unknown state label: {label}")

    psi_full = np.kron(psi_A, psi_B)  # 16×1
    rho0 = psi_full @ psi_full.conj().T  # 16×16
    return rho0


# ============================================================
# Fingerprint extraction
# ============================================================

def extract_fingerprint(times, cpsi_values):
    """Extract fingerprint metrics from a C·Ψ trajectory."""
    threshold = 0.25
    cpsi_arr = np.array(cpsi_values)
    t_arr = np.array(times)

    max_cpsi = np.max(cpsi_arr)
    t_peak = t_arr[np.argmax(cpsi_arr)]

    # Crossing detection
    above = cpsi_arr > threshold
    crossings_up = []
    crossings_down = []
    for i in range(1, len(above)):
        if above[i] and not above[i-1]:
            crossings_up.append(t_arr[i])
        if not above[i] and above[i-1]:
            crossings_down.append(t_arr[i])

    t_up = crossings_up[0] if crossings_up else None
    t_down = crossings_down[0] if crossings_down else None

    # Time above threshold
    dt = t_arr[1] - t_arr[0] if len(t_arr) > 1 else 0
    above_time = np.sum(above) * dt

    # Integrated quantum content
    integral = np.trapezoid(cpsi_arr, t_arr)

    # Rise/fall rates (simple finite difference at peak)
    peak_idx = np.argmax(cpsi_arr)
    if peak_idx > 0:
        rise_rate = (cpsi_arr[peak_idx] - cpsi_arr[peak_idx-1]) / (t_arr[peak_idx] - t_arr[peak_idx-1])
    else:
        rise_rate = 0.0
    if peak_idx < len(cpsi_arr) - 1:
        fall_rate = (cpsi_arr[peak_idx+1] - cpsi_arr[peak_idx]) / (t_arr[peak_idx+1] - t_arr[peak_idx])
    else:
        fall_rate = 0.0

    return {
        'max_cpsi': float(max_cpsi),
        't_peak': float(t_peak),
        't_up': float(t_up) if t_up is not None else None,
        't_down': float(t_down) if t_down is not None else None,
        'above_time': float(above_time),
        'integral': float(integral),
        'rise_rate': float(rise_rate),
        'fall_rate': float(fall_rate),
        'crosses_quarter': bool(max_cpsi > threshold),
    }


# ============================================================
# Main simulation
# ============================================================

def run_single(state_label, J_bridge=0.5, gamma=0.1, dt=0.005, t_max=5.0):
    """Run simulation for one initial B state. Returns times, cpsi_A, cpsi_B, fingerprint."""
    n_steps = int(t_max / dt)
    H = build_hamiltonian(J_internal=1.0, J_bridge=J_bridge)
    L_ops = build_dephasing_operators(gamma)
    rho0 = make_initial_state(state_label)

    # Also compute initial C·Ψ of B
    rho_B0 = partial_trace_A(rho0)
    cpsi_B0 = cpsi(rho_B0)

    print(f"  Running {state_label:8s} (B: C*Psi_0 = {cpsi_B0:.3f}) ...", end=" ", flush=True)

    trajectory = evolve_lindblad(rho0, H, L_ops, dt, n_steps)

    times = []
    cpsi_A_list = []
    cpsi_B_list = []

    # Sample every Nth step to keep output manageable
    sample_every = max(1, n_steps // 1000)

    for i, rho in enumerate(trajectory):
        if i % sample_every == 0:
            t = i * dt
            rho_A = partial_trace_B(rho)
            rho_B = partial_trace_A(rho)
            times.append(t)
            cpsi_A_list.append(cpsi(rho_A))
            cpsi_B_list.append(cpsi(rho_B))

    fp = extract_fingerprint(times, cpsi_A_list)
    fp['state'] = state_label
    fp['cpsi_B0'] = float(cpsi_B0)
    fp['J_bridge'] = J_bridge
    fp['J_over_gamma'] = J_bridge / gamma

    print(f"max(C*Psi_A) = {fp['max_cpsi']:.3f}, crosses 1/4: {fp['crosses_quarter']}")

    return times, cpsi_A_list, cpsi_B_list, fp


def run_all_states(J_bridge=0.5, gamma=0.1, dt=0.005, t_max=5.0):
    """Run all 8 B-states at a given coupling strength."""
    states = ['|++>', '|+0>', 'Bell+', 'Bell-', '|Psi+>', '|+->', '|01>', '|11>']

    results = {}
    for label in states:
        times, cpsi_A, cpsi_B, fp = run_single(label, J_bridge, gamma, dt, t_max)
        results[label] = {
            'times': times,
            'cpsi_A': cpsi_A,
            'cpsi_B': cpsi_B,
            'fingerprint': fp,
        }
    return results


# ============================================================
# Visualization
# ============================================================

def plot_fingerprints_grid(all_coupling_results, output_path):
    """6-panel grid: C*Psi_A trajectories across coupling strengths."""
    j_values = sorted(all_coupling_results.keys())
    n_panels = len(j_values)
    cols = 3
    rows = (n_panels + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(18, 5*rows), sharex=True, sharey=True)
    if rows == 1:
        axes = [axes]

    colors = {
        '|++>': '#e41a1c', '|+0>': '#377eb8', 'Bell+': '#4daf4a',
        'Bell-': '#984ea3', '|Psi+>': '#ff7f00', '|+->': '#a65628',
        '|01>': '#999999', '|11>': '#666666',
    }

    for idx, j_val in enumerate(j_values):
        ax = axes[idx // cols][idx % cols] if rows > 1 else axes[0][idx % cols]
        res = all_coupling_results[j_val]
        gamma = 0.1
        j_over_gamma = j_val / gamma

        for label, data in res.items():
            ax.plot(data['times'], data['cpsi_A'], label=label,
                    color=colors.get(label, 'black'), linewidth=1.2)

        ax.axhline(y=0.25, color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax.set_title(f'J/gamma = {j_over_gamma:.0f} (J_bridge = {j_val})', fontsize=11)
        ax.set_ylim(-0.02, 0.40)
        ax.set_xlabel('Time')
        ax.set_ylabel('C*Psi (System A)')
        ax.legend(fontsize=7, ncol=2, loc='upper right')

    # Hide unused panels
    for idx in range(n_panels, rows * cols):
        axes[idx // cols][idx % cols].set_visible(False)

    plt.suptitle('Bridge Fingerprints: A responds to B across coupling strengths',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_dual_sender_receiver(results, J_bridge, output_path):
    """Sender (B) and Receiver (A) C*Psi comparison at one coupling strength."""
    fig, axes = plt.subplots(2, 4, figsize=(20, 8), sharex=True, sharey=True)

    states_order = ['|++>', '|+0>', 'Bell+', 'Bell-', '|Psi+>', '|+->', '|01>', '|11>']

    for idx, label in enumerate(states_order):
        row = idx // 4
        col = idx % 4
        ax = axes[row][col]
        data = results[label]

        ax.plot(data['times'], data['cpsi_A'], label='A (receiver)', color='blue', linewidth=1.5)
        ax.plot(data['times'], data['cpsi_B'], label='B (sender)', color='red', linewidth=1.5, alpha=0.7)
        ax.axhline(y=0.25, color='gray', linestyle='--', alpha=0.5)
        ax.set_title(f'B0 = {label}', fontsize=11)
        ax.set_ylim(-0.02, 0.50)
        ax.legend(fontsize=8)
        if row == 1:
            ax.set_xlabel('Time')
        if col == 0:
            ax.set_ylabel('C*Psi')

    gamma = 0.1
    plt.suptitle(f'Sender (B) vs Receiver (A) -- J/gamma = {J_bridge/gamma:.0f}', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_phase_portrait(all_coupling_results, output_path):
    """Phase portrait: detection speed (1/t_peak) vs signal strength (max C*Psi_A)."""
    fig, ax = plt.subplots(figsize=(10, 8))

    markers = {
        '|++>': 'o', '|+0>': 's', 'Bell+': '^', 'Bell-': 'v',
        '|Psi+>': 'D', '|+->': 'p', '|01>': '*', '|11>': 'h',
    }
    colors_j = plt.cm.viridis(np.linspace(0.2, 0.9, len(all_coupling_results)))

    for j_idx, (j_val, results) in enumerate(sorted(all_coupling_results.items())):
        for label, data in results.items():
            fp = data['fingerprint']
            speed = 1.0 / fp['t_peak'] if fp['t_peak'] > 0 else 0
            ax.scatter(fp['max_cpsi'], speed, marker=markers.get(label, 'o'),
                       color=colors_j[j_idx], s=80, alpha=0.8, edgecolors='black', linewidth=0.5)

    ax.axvline(x=0.25, color='red', linestyle='--', alpha=0.5, label='1/4 boundary')
    ax.set_xlabel('Signal Strength: max(C*Psi_A)', fontsize=12)
    ax.set_ylabel('Detection Speed: 1/t_peak', fontsize=12)
    ax.set_title('Phase Portrait -- Fingerprint Space', fontsize=14)

    # Custom legend for states
    from matplotlib.lines import Line2D
    state_handles = [Line2D([0], [0], marker=m, color='gray', linestyle='None',
                            markersize=8, label=s) for s, m in markers.items()]
    ax.legend(handles=state_handles, fontsize=9, loc='upper left')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_entanglement_barrier(all_coupling_results, output_path):
    """Entanglement barrier: product vs entangled max signal across coupling strengths."""
    product_states = ['|++>', '|+0>', '|+->', '|01>', '|11>']
    entangled_states = ['Bell+', 'Bell-', '|Psi+>']

    j_values = sorted(all_coupling_results.keys())
    product_max = []
    entangled_max = []

    for j_val in j_values:
        res = all_coupling_results[j_val]
        p_vals = [res[s]['fingerprint']['max_cpsi'] for s in product_states if s in res]
        e_vals = [res[s]['fingerprint']['max_cpsi'] for s in entangled_states if s in res]
        product_max.append(np.mean(p_vals))
        entangled_max.append(np.mean(e_vals))

    j_over_gamma = [j / 0.1 for j in j_values]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(j_over_gamma, product_max, 'o-', color='blue', linewidth=2, markersize=8,
            label='Product states (mean)')
    ax.plot(j_over_gamma, entangled_max, 's-', color='red', linewidth=2, markersize=8,
            label='Entangled states (mean)')
    ax.axhline(y=0.25, color='gray', linestyle='--', alpha=0.5, label='1/4 boundary')
    ax.set_xlabel('J/gamma (coupling / decoherence)', fontsize=12)
    ax.set_ylabel('Mean max(C*Psi_A)', fontsize=12)
    ax.set_title('Entanglement Barrier: Product vs Entangled Signal', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0, max(j_over_gamma) + 1)

    # Add ratio annotation
    for i, j in enumerate(j_over_gamma):
        if entangled_max[i] > 0.001:
            ratio = product_max[i] / entangled_max[i]
            ax.annotate(f'{ratio:.1f}x', xy=(j, product_max[i]),
                        xytext=(0, 10), textcoords='offset points',
                        fontsize=8, ha='center', color='blue')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    gamma = 0.1
    dt = 0.005
    t_max = 5.0

    # Coupling strengths to sweep (J/gamma = 1, 3, 5, 7, 10, 15)
    j_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5]

    output_dir = os.path.dirname(os.path.abspath(__file__))
    viz_dir = os.path.join(output_dir, '..', 'visualizations', 'fingerprints')
    os.makedirs(viz_dir, exist_ok=True)

    all_results = {}

    for j_bridge in j_values:
        j_over_gamma = j_bridge / gamma
        print(f"\n{'='*60}")
        print(f"J_bridge = {j_bridge}, J/gamma = {j_over_gamma:.0f}")
        print(f"{'='*60}")
        all_results[j_bridge] = run_all_states(j_bridge, gamma, dt, t_max)

    # Print summary table (matching doc format)
    print(f"\n{'='*80}")
    print(f"FINGERPRINT SUMMARY TABLE (J/gamma = 5, J_bridge = 0.5)")
    print(f"{'='*80}")
    ref_results = all_results[0.5]
    print(f"{'B State':10s} {'B: C*Psi_0':>10s} {'A: max':>10s} {'Crosses 1/4':>12s} {'Above time':>12s}")
    print(f"{'-'*60}")
    for label in ['|++>', '|+0>', 'Bell+', 'Bell-', '|Psi+>', '|+->', '|01>', '|11>']:
        fp = ref_results[label]['fingerprint']
        crosses = "YES" if fp['crosses_quarter'] else "NEVER"
        above = f"{fp['above_time']:.2f}s" if fp['crosses_quarter'] else "---"
        print(f"{label:10s} {fp['cpsi_B0']:10.3f} {fp['max_cpsi']:10.3f} {crosses:>12s} {above:>12s}")

    # Save fingerprint data as JSON
    fp_data = {}
    for j_val, results in all_results.items():
        fp_data[str(j_val)] = {
            label: data['fingerprint'] for label, data in results.items()
        }
    json_path = os.path.join(viz_dir, 'fingerprints_data.json')
    with open(json_path, 'w') as f:
        json.dump(fp_data, f, indent=2)
    print(f"\nSaved fingerprint data: {json_path}")

    # Generate plots
    print("\nGenerating visualizations...")
    plot_fingerprints_grid(all_results, os.path.join(viz_dir, 'fingerprints_grid.png'))
    plot_dual_sender_receiver(all_results[0.5], 0.5, os.path.join(viz_dir, 'fingerprints_dual.png'))
    plot_phase_portrait(all_results, os.path.join(viz_dir, 'fingerprints_phase.png'))
    plot_entanglement_barrier(all_results, os.path.join(viz_dir, 'fingerprints_barrier.png'))

    print("\nBridge fingerprint simulation complete.")
