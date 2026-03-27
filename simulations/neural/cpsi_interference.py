#!/usr/bin/env python3
"""
CΨ as E-I interference strength.

In quantum: CΨ = purity * coherence. The coherence is the off-diagonal
of rho - the interference between basis states. The standing wave forms
BETWEEN populations (forward) and coherences (backward).

In neural: the standing wave forms BETWEEN E (excitation, forward) and
I (inhibition, backward). The interference lives in the E-I cross-
correlations, not in the individual activities.

The correlation matrix C = <x x^T> has four blocks:
  C_EE  C_EI
  C_IE  C_II

CΨ_neural should measure the E-I cross-correlation (C_EI block),
normalized by the total correlation.

Test: does CΨ defined this way:
  1. Decrease monotonically as the system relaxes?
  2. Cross 1/4 at a parameter-independent point?
  3. Show a change at the crossing?
"""
import numpy as np
from scipy.linalg import expm
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def build_linear_jacobian(W, signs, tau_E, tau_I, alpha):
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


def cpsi_ei_interference(x, signs, x_star=None):
    """CΨ as E-I cross-correlation strength.

    Quantum analog:
      CΨ = Tr(rho^2) * L1(rho) / (d-1)
      = purity * coherence

    Neural analog:
      deviation = x - x_star (perturbation from equilibrium)
      C = deviation * deviation^T (instantaneous outer product)
      C_EI block = cross-correlations between E and I neurons
      C_EE, C_II blocks = within-type correlations

      CΨ = ||C_EI|| / ||C||
      = fraction of total correlation that is E-I interference
    """
    if x_star is None:
        x_star = np.zeros(len(signs))
    delta = x - x_star

    e_idx = np.where(signs > 0)[0]
    i_idx = np.where(signs < 0)[0]

    # Outer product (rank-1 correlation)
    C = np.outer(delta, delta)

    # Block norms
    C_EI = C[np.ix_(e_idx, i_idx)]
    C_IE = C[np.ix_(i_idx, e_idx)]

    cross_norm = np.linalg.norm(C_EI) + np.linalg.norm(C_IE)
    total_norm = np.linalg.norm(C)

    if total_norm < 1e-15:
        return 0

    return cross_norm / (2 * total_norm)  # normalized to [0, 1]


def cpsi_ei_purity_coherence(x, signs, x_star=None):
    """CΨ as PRODUCT of purity and coherence analogs.

    Purity analog: how concentrated the perturbation is.
      = (sum delta_i^2)^2 / (N * sum delta_i^4)
      Range [1/N, 1]. 1 = concentrated on one neuron. 1/N = uniform.

    Coherence analog: E-I cross-correlation fraction (as above).

    CΨ = purity * coherence (product, like quantum).
    """
    if x_star is None:
        x_star = np.zeros(len(signs))
    delta = x - x_star
    N = len(signs)

    if np.linalg.norm(delta) < 1e-15:
        return 0, 0, 0

    # Purity analog: participation ratio (inverse)
    d2 = delta**2
    sum_d2 = np.sum(d2)
    sum_d4 = np.sum(d2**2)
    if sum_d4 < 1e-30:
        purity = 1.0 / N
    else:
        purity = sum_d2**2 / (N * sum_d4)

    # Coherence analog: E-I fraction
    e_idx = np.where(signs > 0)[0]
    i_idx = np.where(signs < 0)[0]
    C = np.outer(delta, delta)
    C_EI = C[np.ix_(e_idx, i_idx)]
    cross_norm = np.linalg.norm(C_EI)
    total_norm = np.linalg.norm(C)
    coherence = cross_norm / total_norm if total_norm > 0 else 0

    return purity * coherence, purity, coherence


# ================================================================
# Test 1: Time evolution of CΨ_EI
# ================================================================
print("=" * 65)
print("CΨ AS E-I INTERFERENCE")
print("Does E-I cross-correlation decrease and cross 1/4?")
print("=" * 65)

N = 20
n_exc = 10

for tau_E, tau_I, alpha, density, seed_net, seed_init in [
    (5.0, 10.0, 0.3, 0.3, 42, 1),
    (5.0, 10.0, 0.3, 0.3, 42, 2),
    (5.0, 10.0, 0.3, 0.3, 42, 3),
    (5.0, 10.0, 0.5, 0.3, 42, 1),
    (5.0, 10.0, 1.0, 0.3, 42, 1),
    (5.0, 15.0, 0.3, 0.3, 42, 1),
    (10.0, 20.0, 0.3, 0.3, 42, 1),
    (5.0, 10.0, 0.3, 0.3, 99, 1),
    (5.0, 10.0, 0.3, 0.5, 42, 1),
]:
    W, signs = make_balanced_network(N, n_exc, density, seed_net)
    J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)

    # Initial perturbation: random, normalized
    rng = np.random.RandomState(seed_init)
    delta0 = rng.randn(N)
    delta0 /= np.linalg.norm(delta0)

    # Time evolution: x(t) = expm(J*t) @ delta0
    times = np.linspace(0, 80, 400)
    cpsi_values = []
    cpsi_pc_values = []

    for t in times:
        x_t = expm(J * t) @ delta0
        cpsi_values.append(cpsi_ei_interference(x_t, signs))
        pc, p, c = cpsi_ei_purity_coherence(x_t, signs)
        cpsi_pc_values.append(pc)

    cpsi_arr = np.array(cpsi_values)
    cpsi_pc_arr = np.array(cpsi_pc_values)

    # Monotonicity check
    diffs = np.diff(cpsi_arr)
    mono_violations = np.sum(diffs > 1e-10)

    # Crossing check
    crossing_t = None
    for k in range(len(cpsi_arr) - 1):
        if cpsi_arr[k] > 0.25 and cpsi_arr[k+1] <= 0.25:
            frac = (0.25 - cpsi_arr[k]) / (cpsi_arr[k+1] - cpsi_arr[k])
            crossing_t = times[k] + frac * (times[k+1] - times[k])
            break

    label = (f"tau={tau_I/tau_E:.0f}x a={alpha} d={density} "
             f"s={seed_net}/{seed_init}")
    cross_str = f"t={crossing_t:.2f}" if crossing_t else "none"

    print(f"  {label:35s}  CΨ_EI(0)={cpsi_arr[0]:.3f}  "
          f"mono_viol={mono_violations:2d}  cross_1/4: {cross_str}")


# ================================================================
# Test 2: Ensemble - is initial CΨ_EI stable?
# ================================================================
print("\n" + "=" * 65)
print("ENSEMBLE: CΨ_EI at t=0 across random initial perturbations")
print("=" * 65)

W, signs = make_balanced_network(N, n_exc, 0.3, 42)
J = build_linear_jacobian(W, signs, 5.0, 10.0, 0.3)

initial_cpsi = []
crossing_times = []

for seed_init in range(200):
    rng = np.random.RandomState(seed_init)
    delta0 = rng.randn(N)
    delta0 /= np.linalg.norm(delta0)

    cpsi_0 = cpsi_ei_interference(delta0, signs)
    initial_cpsi.append(cpsi_0)

    # Find crossing
    times = np.linspace(0, 80, 200)
    for t in times:
        x_t = expm(J * t) @ delta0
        c = cpsi_ei_interference(x_t, signs)
        if c <= 0.25 and cpsi_0 > 0.25:
            crossing_times.append(t)
            break

print(f"  CΨ_EI(t=0): mean={np.mean(initial_cpsi):.4f}  "
      f"std={np.std(initial_cpsi):.4f}  "
      f"range=[{np.min(initial_cpsi):.3f}, {np.max(initial_cpsi):.3f}]")

if crossing_times:
    print(f"  Crossing times (when CΨ_EI(0)>1/4):")
    print(f"    n_crossings={len(crossing_times)}/{len(initial_cpsi)}")
    print(f"    mean_t={np.mean(crossing_times):.2f}  "
          f"std_t={np.std(crossing_times):.2f}")
else:
    print(f"  No crossings found (CΨ_EI(0) always <= 1/4)")

# How many start above 1/4?
above = np.sum(np.array(initial_cpsi) > 0.25)
print(f"  Starting above 1/4: {above}/{len(initial_cpsi)} "
      f"({above/len(initial_cpsi)*100:.0f}%)")


# ================================================================
# Test 3: What IS the E-I interference at t=0 for random vectors?
# ================================================================
print("\n" + "=" * 65)
print("THEORETICAL: Expected CΨ_EI for random perturbations")
print("=" * 65)

print(f"""
For a random unit vector delta in R^N with N_E excitatory and N_I
inhibitory neurons:

  C = delta * delta^T (rank-1 outer product)
  ||C_EI||^2 = sum_{i in E, j in I} delta_i^2 * delta_j^2
  ||C||^2 = (sum delta_i^2)^2 = 1

  E[||C_EI||^2] = (N_E * N_I) / (N * (N+2)) * 2  (approximate)

For N={N}, N_E={n_exc}, N_I={N-n_exc}:
  E[CΨ_EI] ~ sqrt(N_E * N_I) / N = sqrt({n_exc}*{N-n_exc}) / {N}
           = {np.sqrt(n_exc*(N-n_exc))/N:.4f}

Observed mean: {np.mean(initial_cpsi):.4f}
""")

# ================================================================
# Test 4: Eigenmode decomposition of E-I interference
# ================================================================
print("=" * 65)
print("EIGENMODE DECOMPOSITION: Which modes carry E-I interference?")
print("=" * 65)

W, signs = make_balanced_network(N, n_exc, 0.3, 42)
J = build_linear_jacobian(W, signs, 5.0, 10.0, 0.3)
ev, evec = np.linalg.eig(J)

e_idx = np.where(signs > 0)[0]
i_idx = np.where(signs < 0)[0]

print(f"\n{'Mode':>5s}  {'Re(lambda)':>10s}  {'Im(lambda)':>10s}  "
      f"{'|E amp|':>8s}  {'|I amp|':>8s}  {'E-I char':>8s}")
print("-" * 60)

for k in np.argsort(ev.real):
    vec = evec[:, k]
    amp_E = np.sqrt(np.sum(np.abs(vec[e_idx])**2))
    amp_I = np.sqrt(np.sum(np.abs(vec[i_idx])**2))
    total = amp_E + amp_I
    if total > 0:
        # E-I character: 1 when equal, 0 when pure E or pure I
        ei_char = 4 * amp_E * amp_I / (total**2)
    else:
        ei_char = 0

    if abs(ev[k].imag) > 0.001 or ei_char > 0.3:  # show interesting modes
        print(f"  {k:3d}  {ev[k].real:10.4f}  {ev[k].imag:10.4f}  "
              f"{amp_E:8.4f}  {amp_I:8.4f}  {ei_char:8.4f}")

# Which modes are palindromically paired?
rates = [(-ev[k].real, k) for k in range(len(ev))]
rates.sort()
center = (rates[0][0] + rates[-1][0]) / 2
print(f"\n  Rate center: {center:.4f}")
print(f"  Palindromic pairs (rate_i + rate_j ~ 2*center = {2*center:.4f}):")

used = set()
for i, (ri, ki) in enumerate(rates):
    if i in used:
        continue
    target = 2 * center - ri
    for j, (rj, kj) in enumerate(rates):
        if j != i and j not in used and abs(rj - target) < 0.01:
            # Found pair
            vec_i = evec[:, ki]
            vec_j = evec[:, kj]
            ei_i = 4 * np.sum(np.abs(vec_i[e_idx])**2) * np.sum(np.abs(vec_i[i_idx])**2) / (np.sum(np.abs(vec_i)**2)**2)
            ei_j = 4 * np.sum(np.abs(vec_j[e_idx])**2) * np.sum(np.abs(vec_j[i_idx])**2) / (np.sum(np.abs(vec_j)**2)**2)
            print(f"    rate {ri:.4f} + {rj:.4f} = {ri+rj:.4f}  "
                  f"EI_char: {ei_i:.3f}, {ei_j:.3f}")
            used.add(i)
            used.add(j)
            break
