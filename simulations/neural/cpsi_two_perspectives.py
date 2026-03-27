#!/usr/bin/env python3
"""
Two-Perspective CΨ: the E-side and the I-side.

Quantum: CΨ and CΨ_Pi are two perspectives. Pi maps one to the other.
  The palindrome couples them: lambda + lambda' = -2Sg.
  Each mode has XY-weight (how much "coherence" vs "population").
  Pi swaps the character: low XY-weight <-> high XY-weight.

Neural: CΨ_E and CΨ_I are two perspectives. Q maps one to the other.
  The palindrome couples them: Q*J*Q + J + 2S = 0.
  Each mode has E-character and I-character.
  Q should swap the character: high E-char <-> high I-char.

Test:
  1. For each palindromic pair (k, k'):
     Does E-character of k ≈ I-character of k'? (character swap)
  2. Is there a mode at the CENTER where E-char = I-char = 0.5?
  3. Does the center mode correspond to the "fold"?
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


# ================================================================
# Analysis: E/I character of eigenmodes and palindromic pairing
# ================================================================

N = 20
n_exc = 10
tau_E, tau_I = 5.0, 10.0

print("=" * 65)
print("TWO-PERSPECTIVE CΨ: E-side and I-side")
print("=" * 65)

for alpha in [0.3, 0.5, 1.0]:
    W, signs = make_balanced_network(N, n_exc, 0.3, 42)
    J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)
    ev, evec = np.linalg.eig(J)

    e_idx = np.where(signs > 0)[0]
    i_idx = np.where(signs < 0)[0]

    # E/I character for each mode
    e_char = np.zeros(N)
    i_char = np.zeros(N)
    for k in range(N):
        vec = evec[:, k]
        norm2 = np.sum(np.abs(vec)**2)
        if norm2 > 0:
            e_char[k] = np.sum(np.abs(vec[e_idx])**2) / norm2
            i_char[k] = np.sum(np.abs(vec[i_idx])**2) / norm2

    # Find palindromic pairs by rate matching
    rates = [-ev[k].real for k in range(N)]
    center = (min(rates) + max(rates)) / 2
    tol = 0.03 * (max(rates) - min(rates))

    pairs = []
    used = set()
    for k in range(N):
        if k in used:
            continue
        target = 2 * center - rates[k]
        best_j, best_d = -1, tol
        for j in range(N):
            if j != k and j not in used:
                d = abs(rates[j] - target)
                if d < best_d:
                    best_j = j
                    best_d = d
        if best_j >= 0:
            pairs.append((k, best_j))
            used.add(k)
            used.add(best_j)

    print(f"\n--- alpha = {alpha} ({len(pairs)} palindromic pairs) ---")
    print(f"{'Pair':>5s}  {'rate_k':>8s}  {'rate_k_':>8s}  {'sum':>8s}  "
          f"{'E(k)':>6s}  {'I(k)':>6s}  {'E(k_)':>6s}  {'I(k_)':>6s}  "
          f"{'E(k)~I(k_)?':>11s}")
    print("-" * 85)

    char_swap_errors = []
    for k, kp in pairs:
        r_k = rates[k]
        r_kp = rates[kp]
        # Check character swap: E-char(k) should ≈ I-char(k')
        swap_err = abs(e_char[k] - i_char[kp]) + abs(i_char[k] - e_char[kp])
        char_swap_errors.append(swap_err)

        swap_ok = "YES" if swap_err < 0.3 else "no"
        print(f"  {k:2d},{kp:<2d}  {r_k:8.4f}  {r_kp:8.4f}  {r_k+r_kp:8.4f}  "
              f"{e_char[k]:6.3f}  {i_char[k]:6.3f}  {e_char[kp]:6.3f}  "
              f"{i_char[kp]:6.3f}  {swap_ok:>11s} ({swap_err:.3f})")

    mean_swap = np.mean(char_swap_errors)
    print(f"\n  Mean character-swap error: {mean_swap:.4f}")
    print(f"  (0 = perfect swap, like quantum Pi; >1 = no swap)")

    # Find the center mode (closest to E-char = I-char = 0.5)
    center_dists = np.abs(e_char - 0.5)
    center_mode = np.argmin(center_dists)
    print(f"\n  Center mode (E-char closest to 0.5):")
    print(f"    Mode {center_mode}: E-char={e_char[center_mode]:.3f}  "
          f"I-char={i_char[center_mode]:.3f}  "
          f"rate={rates[center_mode]:.4f}  "
          f"osc={abs(ev[center_mode].imag):.4f}")
    print(f"    Rate center: {center:.4f}")


# ================================================================
# Time evolution: track CΨ_E(t) and CΨ_I(t) separately
# ================================================================
print("\n" + "=" * 65)
print("TIME EVOLUTION: CΨ_E(t) and CΨ_I(t)")
print("Does CΨ_E + CΨ_I stay constant? Do they cross?")
print("=" * 65)

W, signs = make_balanced_network(N, n_exc, 0.3, 42)
J = build_linear_jacobian(W, signs, tau_E, tau_I, 0.3)
e_idx = np.where(signs > 0)[0]
i_idx = np.where(signs < 0)[0]

rng = np.random.RandomState(1)
delta0 = rng.randn(N)
delta0 /= np.linalg.norm(delta0)

times = np.linspace(0, 60, 300)
cpsi_E_vals = []
cpsi_I_vals = []

for t in times:
    x_t = expm(J * t) @ delta0
    # CΨ_E = squared energy of E-component
    norm_E = np.sum(x_t[e_idx]**2)
    norm_I = np.sum(x_t[i_idx]**2)
    total = norm_E + norm_I
    if total > 1e-30:
        cpsi_E_vals.append(norm_E / total)
        cpsi_I_vals.append(norm_I / total)
    else:
        cpsi_E_vals.append(0.5)
        cpsi_I_vals.append(0.5)

cpsi_E = np.array(cpsi_E_vals)
cpsi_I = np.array(cpsi_I_vals)
cpsi_sum = cpsi_E + cpsi_I

# Find crossing point
crossing_t = None
for k in range(len(times) - 1):
    if (cpsi_E[k] > cpsi_I[k] and cpsi_E[k+1] <= cpsi_I[k+1]) or \
       (cpsi_E[k] < cpsi_I[k] and cpsi_E[k+1] >= cpsi_I[k+1]):
        frac = (0.5 - cpsi_E[k]) / (cpsi_E[k+1] - cpsi_E[k]) if abs(cpsi_E[k+1] - cpsi_E[k]) > 1e-15 else 0
        crossing_t = times[k] + frac * (times[k+1] - times[k])
        crossing_val = (cpsi_E[k] + cpsi_E[k+1]) / 2
        break

print(f"\n{'t':>6s}  {'CΨ_E':>8s}  {'CΨ_I':>8s}  {'sum':>8s}  {'E-I':>8s}")
print("-" * 45)
for k in range(0, len(times), 15):
    print(f"  {times[k]:4.1f}  {cpsi_E[k]:8.4f}  {cpsi_I[k]:8.4f}  "
          f"{cpsi_sum[k]:8.4f}  {cpsi_E[k]-cpsi_I[k]:+8.4f}")

print(f"\n  Sum CΨ_E + CΨ_I: mean={np.mean(cpsi_sum):.4f}  "
      f"std={np.std(cpsi_sum):.6f}  "
      f"{'CONSTANT' if np.std(cpsi_sum) < 0.01 else 'NOT constant'}")

if crossing_t:
    print(f"\n  >>> CROSSING at t = {crossing_t:.2f}")
    print(f"  >>> CΨ_E = CΨ_I = {crossing_val:.4f} at crossing")
    print(f"  >>> (quantum analog would be CΨ = 1/4 = 0.250)")
else:
    print(f"\n  No crossing found (CΨ_E {'>' if cpsi_E[0] > cpsi_I[0] else '<'} CΨ_I throughout)")


# ================================================================
# Multiple initial conditions: is the crossing value stable?
# ================================================================
print("\n" + "=" * 65)
print("STABILITY: Is the crossing value CΨ_E = CΨ_I always the same?")
print("=" * 65)

crossing_values = []
crossing_times_list = []

for seed_init in range(100):
    rng = np.random.RandomState(seed_init)
    delta0 = rng.randn(N)
    delta0 /= np.linalg.norm(delta0)

    times_fine = np.linspace(0, 60, 600)
    for k in range(len(times_fine) - 1):
        x_k = expm(J * times_fine[k]) @ delta0
        x_k1 = expm(J * times_fine[k+1]) @ delta0

        eE_k = np.sum(x_k[e_idx]**2)
        eI_k = np.sum(x_k[i_idx]**2)
        eE_k1 = np.sum(x_k1[e_idx]**2)
        eI_k1 = np.sum(x_k1[i_idx]**2)

        tot_k = eE_k + eI_k
        tot_k1 = eE_k1 + eI_k1
        if tot_k < 1e-30 or tot_k1 < 1e-30:
            continue

        cE_k = eE_k / tot_k
        cI_k = eI_k / tot_k
        cE_k1 = eE_k1 / tot_k1
        cI_k1 = eI_k1 / tot_k1

        if (cE_k - cI_k) * (cE_k1 - cI_k1) < 0:  # sign change
            # Crossing found
            frac = abs(cE_k - cI_k) / (abs(cE_k - cI_k) + abs(cE_k1 - cI_k1))
            cross_t = times_fine[k] + frac * (times_fine[k+1] - times_fine[k])
            cross_val = (cE_k + cE_k1) / 2  # approximately 0.5 at crossing
            crossing_values.append(cross_val)
            crossing_times_list.append(cross_t)
            break

if crossing_values:
    print(f"\n  Crossings found: {len(crossing_values)} / 100")
    print(f"  Crossing value CΨ_E = CΨ_I:")
    print(f"    mean = {np.mean(crossing_values):.6f}")
    print(f"    std  = {np.std(crossing_values):.6f}")
    print(f"    range = [{np.min(crossing_values):.4f}, {np.max(crossing_values):.4f}]")
    print(f"  Crossing time:")
    print(f"    mean = {np.mean(crossing_times_list):.2f}")
    print(f"    std  = {np.std(crossing_times_list):.2f}")

    print(f"\n  >>> Crossing value = {np.mean(crossing_values):.4f}")
    if abs(np.mean(crossing_values) - 0.5) < 0.01:
        print(f"  >>> This is 0.5 (trivial: E and I contribute equally)")
        print(f"  >>> The FRACTIONS CΨ_E/(CΨ_E+CΨ_I) always sum to 1")
        print(f"  >>> Crossing at 0.5 is geometric, not dynamical")
    if abs(np.mean(crossing_values) - 0.25) < 0.05:
        print(f"  >>> Close to 1/4! The neural fold?")
else:
    print(f"\n  No crossings found in any of 100 trials")


# ================================================================
# Alternative: CΨ_E and CΨ_I as ABSOLUTE energies (not fractions)
# ================================================================
print("\n" + "=" * 65)
print("ALTERNATIVE: Absolute CΨ_E(t) and CΨ_I(t)")
print("Not normalized by total. Do they individually cross 1/4?")
print("=" * 65)

delta0 = np.random.RandomState(1).randn(N)
norm0 = np.linalg.norm(delta0)

print(f"\n{'t':>6s}  {'CΨ_E':>8s}  {'CΨ_I':>8s}  {'sum':>8s}")
print("-" * 38)

for t in np.linspace(0, 60, 20):
    x_t = expm(J * t) @ delta0
    cE = np.sum(x_t[e_idx]**2) / norm0**2
    cI = np.sum(x_t[i_idx]**2) / norm0**2
    marker = ""
    if abs(cE - 0.25) < 0.02 or abs(cI - 0.25) < 0.02:
        marker = " <-- near 1/4"
    print(f"  {t:4.1f}  {cE:8.4f}  {cI:8.4f}  {cE+cI:8.4f}{marker}")
