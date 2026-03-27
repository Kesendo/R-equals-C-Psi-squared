#!/usr/bin/env python3
"""
Test palindromic pairing in LARGER, DENSER balanced subnetworks.

N=10 at density 0.02 is too sparse - the palindrome is trivially true.
Here we test N=20, 30, 50 with 1:1 balance and measure:
1. The actual internal density of the subnetwork
2. Whether pair sums are closer to constant for C. elegans vs random
3. The antisymmetry condition quality for real vs random coupling
"""
import json
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(SCRIPT_DIR, "celegans_connectome.json")) as f:
    data = json.load(f)
W_chem = np.array(data["chemical"])
signs_full = np.array(data["chemical_sign"])
N_full = len(signs_full)

W_signed_full = np.zeros((N_full, N_full))
for i in range(N_full):
    for j in range(N_full):
        W_signed_full[i, j] = signs_full[j] * W_chem[j, i]
max_w = np.max(np.abs(W_signed_full))
W_norm_full = W_signed_full / max_w

exc_idx = np.where(signs_full > 0)[0]
inh_idx = np.where(signs_full < 0)[0]


def build_jacobian(W, tau_exc, tau_inh, signs, alpha=0.3):
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_exc if signs[i] > 0 else tau_inh
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def pair_sum_stats(eigenvalues):
    """Compute pair sum statistics (pairing smallest with largest real part)."""
    evs = sorted(eigenvalues, key=lambda x: x.real)
    n = len(evs)
    sums = []
    for k in range(n // 2):
        sums.append(evs[k].real + evs[n - 1 - k].real)
    return np.array(sums)


def antisymmetry_quality(W, signs, perm, tau_E, tau_I):
    """Measure how well W satisfies the antisymmetry condition.

    Condition: (1/tau_{Q(i)}) * W[Q(i),Q(j)] + (1/tau_i) * W[i,j] = 0
    Returns relative residual.
    """
    n = len(signs)
    num, denom = 0.0, 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            tau_i = tau_E if signs[i] > 0 else tau_I
            qi, qj = perm[i], perm[j]
            tau_qi = tau_E if signs[qi] > 0 else tau_I
            lhs = W[qi, qj] / tau_qi + W[i, j] / tau_i
            rhs = W[i, j] / tau_i
            num += lhs**2
            denom += rhs**2
    return np.sqrt(num / denom) if denom > 0 else 0.0


tau_E, tau_I = 10.0, 20.0
predicted_sum = -(1.0 / tau_E + 1.0 / tau_I)

# === Test balanced subnetworks at increasing size ===
print("=" * 70)
print("PALINDROMIC PAIRING IN DENSE BALANCED SUBNETWORKS")
print("=" * 70)
print(f"tau_E={tau_E}, tau_I={tau_I}, alpha=0.3")
print(f"Predicted pair sum: {predicted_sum:.4f}")
print(f"C. elegans: {len(exc_idx)} E neurons, {len(inh_idx)} I neurons")
print()

for n_half in [5, 10, 13]:
    n_total = 2 * n_half
    # C. elegans: limited by I-neurons (only 26)
    if n_half > len(inh_idx):
        print(f"N={n_total}: skipped (only {len(inh_idx)} I-neurons available)")
        continue

    print(f"\n--- N={n_total} ({n_half}E + {n_half}I) ---")

    ce_pair_stds = []
    rand_pair_stds = []
    ce_densities = []
    ce_antisym = []
    rand_antisym = []

    n_trials = min(200, len(inh_idx))  # can't sample more I than exist

    for trial in range(n_trials):
        rng = np.random.RandomState(trial + 100)

        # C. elegans subnetwork
        e_pick = rng.choice(exc_idx, n_half, replace=False)
        i_pick = rng.choice(inh_idx, n_half, replace=False)
        idx = np.concatenate([e_pick, i_pick])
        W_sub = W_norm_full[np.ix_(idx, idx)]
        signs_sub = signs_full[idx]

        # Density
        n_nonzero = np.count_nonzero(W_sub) - np.count_nonzero(np.diag(W_sub))
        density = n_nonzero / (n_total * (n_total - 1))
        ce_densities.append(density)

        # Build E-I swap permutation for this subnetwork
        e_local = np.where(signs_sub > 0)[0]
        i_local = np.where(signs_sub < 0)[0]
        perm_local = np.arange(n_total)
        for k in range(min(len(e_local), len(i_local))):
            perm_local[e_local[k]] = i_local[k]
            perm_local[i_local[k]] = e_local[k]

        # Eigenvalues and pairing
        J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, alpha=0.3)
        ev_ce = np.linalg.eigvals(J_ce)
        sums_ce = pair_sum_stats(ev_ce)
        ce_pair_stds.append(np.std(sums_ce))

        # Antisymmetry quality
        ce_antisym.append(antisymmetry_quality(
            W_sub, signs_sub, perm_local, tau_E, tau_I))

        # Random control: same density, same signs, random weights
        W_rand = np.zeros((n_total, n_total))
        for i in range(n_total):
            for j in range(n_total):
                if i != j and rng.random() < density:
                    W_rand[i, j] = signs_sub[j] * rng.exponential(0.3)
        mx = np.max(np.abs(W_rand))
        if mx > 0:
            W_rand /= mx

        J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, alpha=0.3)
        ev_rand = np.linalg.eigvals(J_rand)
        sums_rand = pair_sum_stats(ev_rand)
        rand_pair_stds.append(np.std(sums_rand))
        rand_antisym.append(antisymmetry_quality(
            W_rand, signs_sub, perm_local, tau_E, tau_I))

    print(f"  Internal density: {np.mean(ce_densities):.3f} "
          f"(range {np.min(ce_densities):.3f}-{np.max(ce_densities):.3f})")
    print()
    print(f"  Pair-sum std (lower = more palindromic):")
    print(f"    C. elegans:  mean={np.mean(ce_pair_stds):.6f}  "
          f"std={np.std(ce_pair_stds):.6f}")
    print(f"    Random:      mean={np.mean(rand_pair_stds):.6f}  "
          f"std={np.std(rand_pair_stds):.6f}")
    ratio = np.mean(ce_pair_stds) / np.mean(rand_pair_stds) if np.mean(rand_pair_stds) > 0 else 0
    print(f"    Ratio (C.e./random): {ratio:.3f}")

    print()
    print(f"  Antisymmetry quality (lower = more palindromic):")
    print(f"    C. elegans:  mean={np.mean(ce_antisym):.4f}  "
          f"std={np.std(ce_antisym):.4f}")
    print(f"    Random:      mean={np.mean(rand_antisym):.4f}  "
          f"std={np.std(rand_antisym):.4f}")
    a_ratio = np.mean(ce_antisym) / np.mean(rand_antisym) if np.mean(rand_antisym) > 0 else 0
    print(f"    Ratio (C.e./random): {a_ratio:.3f}")

    if ratio < 0.8:
        print(f"  >>> C. elegans is MORE palindromic than random (ratio {ratio:.2f})")
    elif ratio > 1.2:
        print(f"  >>> C. elegans is LESS palindromic than random (ratio {ratio:.2f})")
    else:
        print(f"  >>> C. elegans and random are COMPARABLE (ratio {ratio:.2f})")


# === Full connectome with forced balance ===
print("\n" + "=" * 70)
print("FULL CONNECTOME WITH FORCED E/I BALANCE")
print("=" * 70)

# Take all 300 neurons but relabel half as inhibitory
# (keeping original connectivity, just changing signs)
for balance in [(274, 26), (200, 100), (150, 150)]:
    n_e, n_i = balance
    signs_forced = np.ones(N_full)
    # Randomly assign n_i neurons as inhibitory
    np.random.seed(42)
    i_assign = np.random.choice(N_full, n_i, replace=False)
    signs_forced[i_assign] = -1

    # Rebuild signed weights with new signs
    W_forced = np.zeros((N_full, N_full))
    for i in range(N_full):
        for j in range(N_full):
            W_forced[i, j] = signs_forced[j] * W_chem[j, i]
    mx = np.max(np.abs(W_forced))
    W_forced /= mx

    J_full = build_jacobian(W_forced, tau_E, tau_I, signs_forced, alpha=0.05)
    ev_full = np.linalg.eigvals(J_full)
    sums_full = pair_sum_stats(ev_full)

    n_nonzero = np.count_nonzero(W_forced) - np.count_nonzero(np.diag(W_forced))
    density = n_nonzero / (N_full * (N_full - 1))

    print(f"\n  N=300, E={n_e}, I={n_i}, density={density:.3f}, alpha=0.05")
    print(f"  Pair sums: mean={np.mean(sums_full):.6f}, std={np.std(sums_full):.6f}")
    print(f"  Predicted: {predicted_sum:.6f}")
    print(f"  Pair-sum range: [{np.min(sums_full):.4f}, {np.max(sums_full):.4f}]")

    # Traditional palindrome check
    rates = sorted(-ev.real for ev in ev_full if abs(ev.real) > 1e-10)
    if len(rates) >= 2:
        center = (min(rates) + max(rates)) / 2.0
        tol = 0.03 * (max(rates) - min(rates))
        paired = 0
        used = [False] * len(rates)
        for i_r in range(len(rates)):
            if used[i_r]:
                continue
            target = 2 * center - rates[i_r]
            for j_r in range(len(rates)):
                if i_r != j_r and not used[j_r] and abs(rates[j_r] - target) < tol:
                    used[i_r] = True
                    used[j_r] = True
                    paired += 2
                    break
        print(f"  Traditional pairing: {paired}/{len(rates)} = "
              f"{paired/len(rates)*100:.1f}%")
