#!/usr/bin/env python3
"""
Find CΨ in Wilson-Cowan: The Exploratory Phase

Five CΨ candidates tested across four sweeps.
We do not know which is right. We test all of them.

Candidates:
  1. E_freq / (E_freq + E_decay)   - oscillation vs decay energy ratio
  2. Palindromic fraction           - paired oscillation / total oscillation
  3. Distance from equilibrium      - purity analog (dynamic)
  4. Cross-type coupling fraction   - coherence analog (static)
  5. R = CΨ^2 direct               - amplitude at E/I boundary

Sweeps:
  1. Coupling strength alpha        - does CΨ cross 1/4?
  2. E/I ratio                      - does balance matter?
  3. Time evolution                  - monotonic decrease?
  4. Two coupled networks           - V-Effect analog?

Criteria for a valid CΨ_neural:
  1. Decreases over time
  2. Crosses a parameter-independent value (like 1/4)
  3. Observable change at the crossing
  4. Connected to palindromic structure
"""
import numpy as np
from scipy.integrate import solve_ivp
import os
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# === Wilson-Cowan infrastructure ===

def sigmoid(x, a=1.3, theta=4.0):
    arg = -a * (x - theta)
    arg = np.clip(arg, -500, 500)
    return 1.0 / (1.0 + np.exp(arg))


def dsigmoid(x, a=1.3, theta=4.0):
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)


def build_wc_jacobian(W, signs, tau_E, tau_I, alpha,
                       a_E=1.3, theta_E=4.0, a_I=2.0, theta_I=3.7,
                       x_star=None):
    """Build Wilson-Cowan Jacobian at steady state.

    Each neuron i has activity x_i. The dynamics are:
      dx_i/dt = (-x_i + S_i(alpha * sum_j W[i,j]*x_j + P)) / tau_i

    Linearized around x_star:
      J[i,j] = alpha * W[i,j] * S'_i(input_i) / tau_i    (i != j)
      J[i,i] = (-1 + alpha * W[i,i] * S'_i(input_i)) / tau_i
    """
    n = len(signs)

    if x_star is None:
        # Approximate steady state: iterate
        x = np.ones(n) * 0.3
        P = 1.5
        for _ in range(200):
            inputs = alpha * W @ x + P
            for i in range(n):
                a_i = a_E if signs[i] > 0 else a_I
                th_i = theta_E if signs[i] > 0 else theta_I
                tau_i = tau_E if signs[i] > 0 else tau_I
                x[i] = sigmoid(inputs[i], a_i, th_i)
        x_star = x

    # Compute sigmoid derivatives at steady state
    inputs = alpha * W @ x_star + 1.5
    J = np.zeros((n, n))
    for i in range(n):
        a_i = a_E if signs[i] > 0 else a_I
        th_i = theta_E if signs[i] > 0 else theta_I
        tau_i = tau_E if signs[i] > 0 else tau_I
        dS_i = dsigmoid(inputs[i], a_i, th_i)

        J[i, i] = (-1.0 + alpha * W[i, i] * dS_i) / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] * dS_i / tau_i

    return J, x_star


def build_linear_jacobian(W, signs, tau_E, tau_I, alpha):
    """Simple linear Jacobian (no sigmoid, direct coupling)."""
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_E if signs[i] > 0 else tau_I
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def make_balanced_network(N, n_exc, density=0.3, seed=42):
    """Create a random balanced network with Dale's law."""
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    inh_idx = rng.choice(N, N - n_exc, replace=False)
    signs[inh_idx] = -1

    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (N, N))

    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if mask[i, j]:
                W[i, j] = signs[j] * weights[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


# === CΨ candidate functions ===

def cpsi_1_energy_ratio(eigenvalues):
    """Candidate 1: E_freq / (E_freq + E_decay).
    Measures oscillation vs decay balance."""
    E_freq = np.sum(np.abs(np.imag(eigenvalues)))
    E_decay = np.sum(np.abs(np.real(eigenvalues)))
    if E_freq + E_decay == 0:
        return 0
    return E_freq / (E_freq + E_decay)


def cpsi_2_palindromic_fraction(eigenvalues, tol_frac=0.03):
    """Candidate 2: fraction of oscillation energy in palindromic pairs."""
    rates = sorted(-ev.real for ev in eigenvalues if abs(ev.real) > 1e-10)
    if len(rates) < 2:
        return 0
    center = (min(rates) + max(rates)) / 2.0
    tol = tol_frac * (max(rates) - min(rates))

    paired_osc = 0
    total_osc = np.sum(np.abs(np.imag(eigenvalues)))
    if total_osc == 0:
        return 0

    used = [False] * len(rates)
    paired_indices = []
    for i in range(len(rates)):
        if used[i]:
            continue
        target = 2 * center - rates[i]
        for j in range(len(rates)):
            if i != j and not used[j] and abs(rates[j] - target) < tol:
                used[i] = True
                used[j] = True
                paired_indices.extend([i, j])
                break

    # Sum oscillation energy of paired eigenvalues
    evs_sorted = sorted(eigenvalues, key=lambda x: -x.real)
    for idx in paired_indices:
        if idx < len(evs_sorted):
            paired_osc += abs(evs_sorted[idx].imag)

    return paired_osc / total_osc if total_osc > 0 else 0


def cpsi_3_distance(x, x_star, x_max_norm=1.0):
    """Candidate 3: ||x - x*|| / ||x_max - x*||. Purity analog."""
    dist = np.linalg.norm(x - x_star)
    return min(dist / x_max_norm, 1.0) if x_max_norm > 0 else 0


def cpsi_4_crosstype_fraction(W, signs):
    """Candidate 4: cross-type coupling / total coupling. Static."""
    n = len(signs)
    cross = 0
    same = 0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            t_i = 'E' if signs[i] > 0 else 'I'
            t_j = 'E' if signs[j] > 0 else 'I'
            w = abs(W[i, j])
            if t_i != t_j:
                cross += w
            else:
                same += w
    total = cross + same
    return cross / total if total > 0 else 0


def cpsi_5_boundary_amplitude(eigenvalues, eigenvectors, signs):
    """Candidate 5: amplitude of oscillation at E/I boundary.

    For each eigenmode, compute how much amplitude sits at the
    boundary between E and I neurons. The "boundary" modes are
    those with significant amplitude at both E and I sites.
    """
    n = len(signs)
    e_idx = np.where(signs > 0)[0]
    i_idx = np.where(signs < 0)[0]

    boundary_energy = 0
    total_energy = 0

    for k in range(len(eigenvalues)):
        ev = eigenvalues[k]
        vec = eigenvectors[:, k]
        osc = abs(ev.imag)

        # Amplitude at E and I sites
        amp_E = np.sum(np.abs(vec[e_idx])**2)
        amp_I = np.sum(np.abs(vec[i_idx])**2)
        total_amp = amp_E + amp_I

        if total_amp > 0:
            # Boundary character: high when both E and I contribute
            # Maximum at equal contribution (amp_E = amp_I = 0.5)
            boundary = 4 * amp_E * amp_I / (total_amp**2)  # range [0, 1]
            boundary_energy += osc * boundary
            total_energy += osc

    return boundary_energy / total_energy if total_energy > 0 else 0


# === Sweep functions ===

def sweep_coupling(N=50, n_exc=25, density=0.3, tau_E=5.0, tau_I=10.0):
    """Sweep 1: coupling strength alpha."""
    W, signs = make_balanced_network(N, n_exc, density=density, seed=42)
    alphas = [0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]

    print(f"{'alpha':>6s}  {'CΨ_1':>6s}  {'CΨ_2':>6s}  {'CΨ_4':>6s}  {'CΨ_5':>6s}  "
          f"{'cross?':>6s}")
    print("-" * 50)

    for alpha in alphas:
        J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)
        ev = np.linalg.eigvals(J)
        _, evec = np.linalg.eig(J)

        c1 = cpsi_1_energy_ratio(ev)
        c2 = cpsi_2_palindromic_fraction(ev)
        c4 = cpsi_4_crosstype_fraction(W, signs)
        c5 = cpsi_5_boundary_amplitude(ev, evec, signs)

        # Does any candidate cross 1/4?
        crosses = []
        if abs(c1 - 0.25) < 0.02:
            crosses.append("C1")
        if abs(c5 - 0.25) < 0.02:
            crosses.append("C5")

        cross_str = ",".join(crosses) if crosses else ""
        print(f"  {alpha:4.2f}  {c1:6.4f}  {c2:6.4f}  {c4:6.4f}  {c5:6.4f}  {cross_str}")

    return W, signs


def sweep_ei_ratio(N=50, density=0.3, alpha=0.3, tau_E=5.0, tau_I=10.0):
    """Sweep 2: E/I ratio."""
    ratios = [(50, 0), (40, 10), (30, 20), (25, 25), (20, 30), (10, 40)]

    print(f"{'E:I':>6s}  {'CΨ_1':>6s}  {'CΨ_2':>6s}  {'CΨ_4':>6s}  {'CΨ_5':>6s}")
    print("-" * 40)

    for n_e, n_i in ratios:
        W, signs = make_balanced_network(N, n_e, density=density, seed=42)
        J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)
        ev = np.linalg.eigvals(J)
        _, evec = np.linalg.eig(J)

        c1 = cpsi_1_energy_ratio(ev)
        c2 = cpsi_2_palindromic_fraction(ev)
        c4 = cpsi_4_crosstype_fraction(W, signs)
        c5 = cpsi_5_boundary_amplitude(ev, evec, signs)

        print(f" {n_e:2d}:{n_i:<2d}  {c1:6.4f}  {c2:6.4f}  {c4:6.4f}  {c5:6.4f}")


def sweep_time_evolution(N=50, n_exc=25, density=0.3, alpha=0.3,
                         tau_E=5.0, tau_I=10.0):
    """Sweep 3: time evolution of CΨ candidates."""
    W, signs = make_balanced_network(N, n_exc, density=density, seed=42)
    J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)

    # Steady state
    x_star = np.zeros(N)

    # Initial perturbation
    rng = np.random.RandomState(123)
    delta = rng.randn(N) * 0.1
    x0 = x_star + delta
    x_max_norm = np.linalg.norm(delta)

    # Linear time evolution: x(t) = x* + expm(J*t) @ delta
    from scipy.linalg import expm

    times = np.linspace(0, 50, 100)

    print(f"{'t':>6s}  {'CΨ_1':>6s}  {'CΨ_3':>6s}  {'CΨ_5':>6s}  {'monotonic?':>10s}")
    print("-" * 45)

    prev_c3 = None
    monotonic_violations = 0

    for t in times[::5]:  # every 5th point for display
        x_t = x_star + expm(J * t) @ delta
        ev_t = np.linalg.eigvals(J)  # J doesn't change (linear)
        _, evec_t = np.linalg.eig(J)

        c1 = cpsi_1_energy_ratio(ev_t)
        c3 = cpsi_3_distance(x_t, x_star, x_max_norm)
        c5 = cpsi_5_boundary_amplitude(ev_t, evec_t, signs)

        mono = ""
        if prev_c3 is not None:
            if c3 > prev_c3 + 1e-10:
                mono = "VIOLATION"
                monotonic_violations += 1
        prev_c3 = c3

        if t < 5 or t > 45 or abs(t - 25) < 3 or mono:
            print(f"  {t:4.1f}  {c1:6.4f}  {c3:6.4f}  {c5:6.4f}  {mono}")

    print(f"\n  CΨ_3 monotonicity violations: {monotonic_violations}")
    print(f"  (CΨ_1 and CΨ_5 are Jacobian properties, constant in linear regime)")

    # Check if CΨ_3 crosses 1/4
    c3_values = []
    for t in times:
        x_t = x_star + expm(J * t) @ delta
        c3_values.append(cpsi_3_distance(x_t, x_star, x_max_norm))

    c3_arr = np.array(c3_values)
    crosses_quarter = np.any((c3_arr[:-1] > 0.25) & (c3_arr[1:] <= 0.25))
    if crosses_quarter:
        idx = np.where((c3_arr[:-1] > 0.25) & (c3_arr[1:] <= 0.25))[0][0]
        print(f"  CΨ_3 crosses 1/4 at t ~ {times[idx]:.1f}")
    else:
        print(f"  CΨ_3 does NOT cross 1/4 (range: [{c3_arr.min():.4f}, {c3_arr.max():.4f}])")


def sweep_v_effect(N_single=25, n_exc_single=12, density=0.3, alpha=0.3,
                   tau_E=5.0, tau_I=10.0):
    """Sweep 4: Two coupled networks (V-Effect analog)."""
    # Single network
    W1, signs1 = make_balanced_network(N_single, n_exc_single, density, seed=42)
    J1 = build_linear_jacobian(W1, signs1, tau_E, tau_I, alpha)
    ev1 = np.linalg.eigvals(J1)

    # Count distinct frequencies
    freqs_single = sorted(set(np.round(np.abs(np.imag(ev1)), 4)))
    freqs_single = [f for f in freqs_single if f > 0.001]

    # Two coupled networks through a mediator
    N_total = 2 * N_single + 1  # +1 mediator
    W_coupled = np.zeros((N_total, N_total))
    signs_coupled = np.zeros(N_total)

    # Place network 1 at indices 0..N_single-1
    W_coupled[:N_single, :N_single] = W1
    signs_coupled[:N_single] = signs1

    # Place network 2 at indices N_single..2*N_single-1
    W2, signs2 = make_balanced_network(N_single, n_exc_single, density, seed=99)
    W_coupled[N_single:2*N_single, N_single:2*N_single] = W2
    signs_coupled[N_single:2*N_single] = signs2

    # Mediator at index 2*N_single (excitatory)
    med = 2 * N_single
    signs_coupled[med] = 1.0

    # Couple mediator to both networks (edge nodes)
    coupling_strength = 0.3
    for net_offset in [0, N_single]:
        # Mediator connects to first and last node of each network
        W_coupled[med, net_offset] = coupling_strength
        W_coupled[net_offset, med] = coupling_strength
        W_coupled[med, net_offset + N_single - 1] = coupling_strength
        W_coupled[net_offset + N_single - 1, med] = coupling_strength

    J_coupled = build_linear_jacobian(W_coupled, signs_coupled, tau_E, tau_I, alpha)
    ev_coupled = np.linalg.eigvals(J_coupled)

    freqs_coupled = sorted(set(np.round(np.abs(np.imag(ev_coupled)), 4)))
    freqs_coupled = [f for f in freqs_coupled if f > 0.001]

    # CΨ candidates for both
    _, evec1 = np.linalg.eig(J1)
    _, evec_c = np.linalg.eig(J_coupled)

    c1_single = cpsi_1_energy_ratio(ev1)
    c1_coupled = cpsi_1_energy_ratio(ev_coupled)
    c5_single = cpsi_5_boundary_amplitude(ev1, evec1, signs1)
    c5_coupled = cpsi_5_boundary_amplitude(ev_coupled, evec_c, signs_coupled)

    print(f"  Single network (N={N_single}):  {len(freqs_single)} distinct frequencies")
    print(f"  Coupled (N={N_total}):           {len(freqs_coupled)} distinct frequencies")
    print(f"  Ratio: {len(freqs_coupled)} / (2 x {len(freqs_single)}) = "
          f"{len(freqs_coupled) / (2 * len(freqs_single)):.2f}")
    print(f"\n  CΨ_1 single: {c1_single:.4f}  coupled: {c1_coupled:.4f}")
    print(f"  CΨ_5 single: {c5_single:.4f}  coupled: {c5_coupled:.4f}")

    if len(freqs_coupled) > 2 * len(freqs_single):
        v_effect = len(freqs_coupled) - 2 * len(freqs_single)
        print(f"\n  >>> V-EFFECT DETECTED: {v_effect} new frequencies from coupling")
    else:
        print(f"\n  No V-Effect (coupled frequencies <= 2 x single)")


# === Main ===

if __name__ == '__main__':
    print("=" * 60)
    print("CΨ NEURAL EXPLORATION")
    print("Five candidates, four sweeps, four criteria")
    print("=" * 60)

    N = 50
    n_exc = 25
    tau_E = 5.0
    tau_I = 10.0

    print(f"\nNetwork: N={N}, E={n_exc}, I={N-n_exc}")
    print(f"tau_E={tau_E}ms, tau_I={tau_I}ms, ratio={tau_I/tau_E:.1f}")

    # Sweep 1
    print("\n" + "=" * 60)
    print("SWEEP 1: Coupling Strength")
    print("=" * 60)
    W, signs = sweep_coupling(N, n_exc, tau_E=tau_E, tau_I=tau_I)

    # Sweep 2
    print("\n" + "=" * 60)
    print("SWEEP 2: E/I Ratio")
    print("=" * 60)
    sweep_ei_ratio(N, alpha=0.3, tau_E=tau_E, tau_I=tau_I)

    # Sweep 3
    print("\n" + "=" * 60)
    print("SWEEP 3: Time Evolution")
    print("=" * 60)
    sweep_time_evolution(N, n_exc, alpha=0.3, tau_E=tau_E, tau_I=tau_I)

    # Sweep 4
    print("\n" + "=" * 60)
    print("SWEEP 4: V-Effect (Two Coupled Networks)")
    print("=" * 60)
    sweep_v_effect(tau_E=tau_E, tau_I=tau_I)

    # === Verdict ===
    print("\n" + "=" * 60)
    print("CANDIDATE ASSESSMENT")
    print("=" * 60)
    print("""
Criteria:
  1. Decreases over time (like quantum CΨ)
  2. Crosses a parameter-independent value
  3. Observable change at the crossing
  4. Connected to palindromic structure
""")
