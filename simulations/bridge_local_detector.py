"""
Bridge Local Detector Test (Simulation Test #1 from BRIDGE_PROTOCOL.md)

Question: What local observable on system A can detect the 1/4-crossing,
and can it distinguish Bell+ from product states?

Strategy: Evolve Bell+, |++>, and |00> through the Heisenberg bridge.
At each timestep, compute a battery of local observables on rho_A.
Look for ANY observable that shows different behavior at the crossing point.

Observables on rho_A (4x4 reduced density matrix of 2-qubit subsystem A):
  1. Purity:           Tr(rho_A^2)
  2. Von Neumann S:    -Tr(rho_A log rho_A)
  3. C*Psi:            purity * normalized_coherence
  4. L1 coherence:     sum |rho_ij| for i!=j
  5. Concurrence:      entanglement within A
  6. Single-qubit:     <sigma_x>, <sigma_y>, <sigma_z> on qubit 0 and 1
  7. MAP theta:        arctan(sqrt(4*C*Psi - 1)) -- real above 1/4, imaginary below
  8. Qubit 0 purity:   Tr(rho_q0^2) -- single qubit purity

2026-02-24  Test derived from v033 agent experiment insight
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

# Import shared infrastructure from bridge_fingerprints
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bridge_fingerprints import (
    build_hamiltonian, build_dephasing_operators, evolve_lindblad,
    partial_trace_B, partial_trace_A, purity, l1_coherence, cpsi,
    make_initial_state, sx, sy, sz, I2
)


# ============================================================
# Additional observables not in bridge_fingerprints
# ============================================================

def von_neumann_entropy(rho):
    """S = -Tr(rho log rho). Uses eigenvalues to avoid log(0)."""
    eigvals = np.linalg.eigvalsh(rho)
    eigvals = eigvals[eigvals > 1e-15]  # filter numerical zeros
    return -np.sum(eigvals * np.log2(eigvals))


def concurrence_2qubit(rho):
    """Concurrence for a 2-qubit density matrix (Wootters formula)."""
    # sigma_y tensor sigma_y
    sy_sy = np.kron(sy, sy)
    rho_tilde = sy_sy @ rho.conj() @ sy_sy
    R = rho @ rho_tilde
    eigvals = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0.0, eigvals[0] - eigvals[1] - eigvals[2] - eigvals[3])


def partial_trace_qubit1(rho_2q):
    """Trace out qubit 1 from a 2-qubit system, return rho of qubit 0."""
    rho_q0 = np.zeros((2, 2), dtype=complex)
    rho_q0[0, 0] = rho_2q[0, 0] + rho_2q[1, 1]
    rho_q0[0, 1] = rho_2q[0, 2] + rho_2q[1, 3]
    rho_q0[1, 0] = rho_2q[2, 0] + rho_2q[3, 1]
    rho_q0[1, 1] = rho_2q[2, 2] + rho_2q[3, 3]
    return rho_q0


def partial_trace_qubit0(rho_2q):
    """Trace out qubit 0 from a 2-qubit system, return rho of qubit 1."""
    rho_q1 = np.zeros((2, 2), dtype=complex)
    rho_q1[0, 0] = rho_2q[0, 0] + rho_2q[2, 2]
    rho_q1[0, 1] = rho_2q[0, 1] + rho_2q[2, 3]
    rho_q1[1, 0] = rho_2q[1, 0] + rho_2q[3, 2]
    rho_q1[1, 1] = rho_2q[1, 1] + rho_2q[3, 3]
    return rho_q1


def expectation_pauli(rho_1q, pauli):
    """<sigma> for a single-qubit density matrix."""
    return np.real(np.trace(rho_1q @ pauli))


def map_theta(cpsi_val):
    """MAP theta = arctan(sqrt(4*C*Psi - 1)).
    Returns (theta_real, theta_imag):
      Above 1/4: theta_real > 0, theta_imag = 0
      Below 1/4: theta_real = 0, theta_imag = atanh(sqrt(1 - 4*C*Psi))
    """
    arg = 4.0 * cpsi_val - 1.0
    if arg >= 0:
        return np.arctan(np.sqrt(arg)), 0.0
    else:
        return 0.0, np.arctanh(np.sqrt(-arg)) if abs(arg) < 1 else np.arctanh(np.sqrt(min(-arg, 0.9999)))


# ============================================================
# Compute full observable battery at each timestep
# ============================================================

def compute_observables(rho_A):
    """Compute all local observables on the 2-qubit subsystem A."""
    # Subsystem-level
    pur = purity(rho_A)
    ent = von_neumann_entropy(rho_A)
    cpsi_val = cpsi(rho_A)
    l1 = l1_coherence(rho_A)
    conc = concurrence_2qubit(rho_A)
    theta_r, theta_i = map_theta(cpsi_val)

    # Single-qubit level
    rho_q0 = partial_trace_qubit1(rho_A)
    rho_q1 = partial_trace_qubit0(rho_A)

    q0_x = expectation_pauli(rho_q0, sx)
    q0_y = expectation_pauli(rho_q0, sy)
    q0_z = expectation_pauli(rho_q0, sz)
    q0_pur = np.real(np.trace(rho_q0 @ rho_q0))

    q1_x = expectation_pauli(rho_q1, sx)
    q1_y = expectation_pauli(rho_q1, sy)
    q1_z = expectation_pauli(rho_q1, sz)
    q1_pur = np.real(np.trace(rho_q1 @ rho_q1))

    return {
        'purity': pur,
        'entropy': ent,
        'cpsi': cpsi_val,
        'l1_coherence': l1,
        'concurrence': conc,
        'theta_real': theta_r,
        'theta_imag': theta_i,
        'q0_sx': q0_x, 'q0_sy': q0_y, 'q0_sz': q0_z, 'q0_purity': q0_pur,
        'q1_sx': q1_x, 'q1_sy': q1_y, 'q1_sz': q1_z, 'q1_purity': q1_pur,
    }


# ============================================================
# Run simulation for one state
# ============================================================

def run_detector_test(state_label, J_bridge=0.5, gamma=0.1, dt=0.002, t_max=5.0):
    """Evolve and collect all observables on rho_A at each sampled step."""
    n_steps = int(t_max / dt)
    H = build_hamiltonian(J_internal=1.0, J_bridge=J_bridge)
    L_ops = build_dephasing_operators(gamma)
    rho0 = make_initial_state(state_label)

    print(f"  Evolving {state_label:8s} ({n_steps} steps, dt={dt}) ...", flush=True)
    trajectory = evolve_lindblad(rho0, H, L_ops, dt, n_steps)

    sample_every = max(1, n_steps // 1000)
    times = []
    obs_series = {key: [] for key in [
        'purity', 'entropy', 'cpsi', 'l1_coherence', 'concurrence',
        'theta_real', 'theta_imag',
        'q0_sx', 'q0_sy', 'q0_sz', 'q0_purity',
        'q1_sx', 'q1_sy', 'q1_sz', 'q1_purity',
    ]}

    for i, rho in enumerate(trajectory):
        if i % sample_every == 0:
            rho_A = partial_trace_B(rho)
            obs = compute_observables(rho_A)
            times.append(i * dt)
            for key in obs_series:
                obs_series[key].append(obs[key])

    # Convert to numpy
    for key in obs_series:
        obs_series[key] = np.array(obs_series[key])

    return np.array(times), obs_series


# ============================================================
# Visualization: Observable comparison
# ============================================================

def find_crossing_time(times, cpsi_vals, threshold=0.25):
    """Find first upward crossing of threshold."""
    for i in range(1, len(cpsi_vals)):
        if cpsi_vals[i] >= threshold and cpsi_vals[i-1] < threshold:
            # Linear interpolation
            frac = (threshold - cpsi_vals[i-1]) / (cpsi_vals[i] - cpsi_vals[i-1])
            return times[i-1] + frac * (times[i] - times[i-1])
    return None


def plot_observable_comparison(all_data, output_dir):
    """Multi-panel plot comparing observables across states."""
    states = list(all_data.keys())
    colors = {'Bell+': '#e41a1c', '|++>': '#377eb8', '|00>': '#4daf4a',
              '|+0>': '#984ea3', '|+->': '#ff7f00'}

    # Find crossing times
    crossing_times = {}
    for label, (times, obs) in all_data.items():
        t_cross = find_crossing_time(times, obs['cpsi'])
        crossing_times[label] = t_cross
        if t_cross:
            print(f"  {label:8s} crosses 1/4 at t = {t_cross:.3f}")
        else:
            print(f"  {label:8s} NEVER crosses 1/4")

    # ---- Panel 1: Key observables overview (2x3 grid) ----
    obs_to_plot = [
        ('cpsi', 'C*Psi (MAP metric)', True),
        ('purity', 'Purity Tr(rho_A^2)', False),
        ('entropy', 'Von Neumann Entropy S', False),
        ('concurrence', 'Concurrence (A internal)', False),
        ('l1_coherence', 'L1 Coherence', False),
        ('q0_purity', 'Qubit 0 Purity', False),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    for idx, (obs_key, title, show_quarter) in enumerate(obs_to_plot):
        ax = axes[idx // 3][idx % 3]
        for label in states:
            times, obs = all_data[label]
            ax.plot(times, obs[obs_key], label=label, color=colors.get(label, 'gray'), linewidth=1.5)
            # Mark crossing time
            t_cross = crossing_times[label]
            if t_cross and obs_key == 'cpsi':
                ax.axvline(x=t_cross, color=colors.get(label, 'gray'),
                          linestyle=':', alpha=0.5, linewidth=1)
        if show_quarter:
            ax.axhline(y=0.25, color='red', linestyle='--', alpha=0.5, linewidth=1, label='1/4 threshold')
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(title.split('(')[0].strip())
        if idx >= 3:
            ax.set_xlabel('Time')
        ax.legend(fontsize=8)

    plt.suptitle('Local Detector Test: Observable Battery on System A', fontsize=14, y=1.01)
    plt.tight_layout()
    path1 = os.path.join(output_dir, 'local_detector_overview.png')
    plt.savefig(path1, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path1}")

    # ---- Panel 2: Single-qubit Bloch components ----
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    bloch_obs = [
        ('q0_sx', '<sigma_x> qubit 0'), ('q0_sy', '<sigma_y> qubit 0'), ('q0_sz', '<sigma_z> qubit 0'),
        ('q1_sx', '<sigma_x> qubit 1'), ('q1_sy', '<sigma_y> qubit 1'), ('q1_sz', '<sigma_z> qubit 1'),
    ]
    for idx, (obs_key, title) in enumerate(bloch_obs):
        ax = axes[idx // 3][idx % 3]
        for label in states:
            times, obs = all_data[label]
            ax.plot(times, obs[obs_key], label=label, color=colors.get(label, 'gray'), linewidth=1.5)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(obs_key)
        if idx >= 3:
            ax.set_xlabel('Time')
        ax.legend(fontsize=8)

    plt.suptitle('Local Detector Test: Single-Qubit Bloch Components', fontsize=14, y=1.01)
    plt.tight_layout()
    path2 = os.path.join(output_dir, 'local_detector_bloch.png')
    plt.savefig(path2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path2}")

    # ---- Panel 3: MAP theta (real vs imaginary parts) ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for label in states:
        times, obs = all_data[label]
        axes[0].plot(times, obs['theta_real'], label=label, color=colors.get(label, 'gray'), linewidth=1.5)
        axes[1].plot(times, obs['theta_imag'], label=label, color=colors.get(label, 'gray'), linewidth=1.5)
    axes[0].set_title('MAP theta (real part) -- above 1/4', fontsize=11)
    axes[0].set_xlabel('Time'); axes[0].set_ylabel('theta_real')
    axes[0].legend(fontsize=9)
    axes[1].set_title('MAP theta (imaginary part) -- below 1/4', fontsize=11)
    axes[1].set_xlabel('Time'); axes[1].set_ylabel('theta_imag')
    axes[1].legend(fontsize=9)
    plt.suptitle('MAP Phase Transition: theta = arctan(sqrt(4*C*Psi - 1))', fontsize=13, y=1.02)
    plt.tight_layout()
    path3 = os.path.join(output_dir, 'local_detector_theta.png')
    plt.savefig(path3, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path3}")

    # ---- Panel 4: DIFFERENCE plots (Bell+ minus |++>) ----
    if 'Bell+' in all_data and '|++>' in all_data:
        times_b, obs_b = all_data['Bell+']
        times_p, obs_p = all_data['|++>']

        # Ensure same length (should be if same params)
        n = min(len(times_b), len(times_p))
        diff_obs = ['cpsi', 'purity', 'entropy', 'concurrence', 'l1_coherence', 'q0_purity']
        diff_labels = ['Delta C*Psi', 'Delta Purity', 'Delta Entropy',
                       'Delta Concurrence', 'Delta L1', 'Delta Q0 Purity']

        fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
        for idx, (obs_key, dlabel) in enumerate(zip(diff_obs, diff_labels)):
            ax = axes[idx // 3][idx % 3]
            diff = obs_b[obs_key][:n] - obs_p[obs_key][:n]
            ax.plot(times_b[:n], diff, color='purple', linewidth=1.5)
            ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
            ax.set_title(dlabel + ' (Bell+ minus |++>)', fontsize=11)
            ax.set_ylabel(dlabel)
            if idx >= 3:
                ax.set_xlabel('Time')

            # Mark crossing times
            for label, col in [('Bell+', '#e41a1c'), ('|++>', '#377eb8')]:
                t_cross = crossing_times.get(label)
                if t_cross:
                    ax.axvline(x=t_cross, color=col, linestyle=':', alpha=0.5,
                              linewidth=1, label=f'{label} crosses')
            ax.legend(fontsize=8)

        plt.suptitle('CRITICAL: Bell+ vs |++> Observable Differences', fontsize=14, y=1.01)
        plt.tight_layout()
        path4 = os.path.join(output_dir, 'local_detector_differences.png')
        plt.savefig(path4, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {path4}")

    # ---- Summary: Values AT crossing point ----
    print(f"\n{'='*80}")
    print(f"VALUES AT CROSSING POINT (or nearest if no crossing)")
    print(f"{'='*80}")

    for label in states:
        times, obs = all_data[label]
        t_cross = crossing_times[label]
        if t_cross:
            idx = np.argmin(np.abs(times - t_cross))
        else:
            # Use time of max cpsi instead
            idx = np.argmax(obs['cpsi'])
            t_cross = times[idx]

        print(f"\n  {label} at t = {t_cross:.3f}:")
        print(f"    C*Psi       = {obs['cpsi'][idx]:.6f}")
        print(f"    Purity      = {obs['purity'][idx]:.6f}")
        print(f"    Entropy     = {obs['entropy'][idx]:.6f}")
        print(f"    Concurrence = {obs['concurrence'][idx]:.6f}")
        print(f"    L1 Coh      = {obs['l1_coherence'][idx]:.6f}")
        print(f"    Q0 Purity   = {obs['q0_purity'][idx]:.6f}")
        print(f"    Q0 <sx>     = {obs['q0_sx'][idx]:.6f}")
        print(f"    Q0 <sy>     = {obs['q0_sy'][idx]:.6f}")
        print(f"    Q0 <sz>     = {obs['q0_sz'][idx]:.6f}")
        print(f"    theta_real  = {obs['theta_real'][idx]:.6f}")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    # Parameters matching bridge_fingerprints reference case
    J_bridge = 0.5
    gamma = 0.1
    dt = 0.002    # Slightly coarser than fingerprints for speed, still accurate
    t_max = 5.0

    # States to test: Bell+ (entangled), |++> (product, similar energy),
    # |00> (classical control), |+0> and |+-> (asymmetric products)
    test_states = ['Bell+', '|++>', '|+0>', '|+->', '|01>']

    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'visualizations', 'local_detector')
    os.makedirs(output_dir, exist_ok=True)

    print(f"Bridge Local Detector Test")
    print(f"J_bridge={J_bridge}, gamma={gamma}, J/gamma={J_bridge/gamma:.0f}")
    print(f"dt={dt}, t_max={t_max}")
    print(f"States: {test_states}")
    print(f"{'='*60}")

    all_data = {}
    for label in test_states:
        times, obs = run_detector_test(label, J_bridge, gamma, dt, t_max)
        all_data[label] = (times, obs)

    print(f"\nGenerating plots...")
    plot_observable_comparison(all_data, output_dir)

    print(f"\n{'='*60}")
    print(f"Local detector test complete.")
    print(f"Output: {output_dir}")
