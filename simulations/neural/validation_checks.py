#!/usr/bin/env python3
"""
Validation checks for the neural algebraic palindrome claim.

Before publishing biological results, validate robustness against:
1. Parameter sensitivity (tau ratio, alpha)
2. Pairing sensitivity (sequential vs random vs optimized)
3. Degree-preserving null model (sharper than Erdos-Renyi)
4. Effect size (is 5-8x within normal network variation?)
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


def build_jacobian(W, tau_E, tau_I, signs, alpha=0.3):
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_E if signs[i] > 0 else tau_I
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def build_swap(signs, pairing='sequential'):
    """Build E-I swap permutation with different pairing strategies."""
    n = len(signs)
    e_local = list(np.where(signs > 0)[0])
    i_local = list(np.where(signs < 0)[0])
    n_pairs = min(len(e_local), len(i_local))

    if pairing == 'random':
        np.random.shuffle(e_local)
        np.random.shuffle(i_local)

    perm = np.arange(n)
    for k in range(n_pairs):
        perm[e_local[k]] = i_local[k]
        perm[i_local[k]] = e_local[k]

    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, perm[i]] = 1.0
    return perm, Q


def palindrome_residual(J, Q):
    """Compute coupling residual ||R_off|| / ||J||."""
    n = J.shape[0]
    QJQ = Q @ J @ Q.T
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2 * np.diag(S_diag)
    R_off = R - np.diag(np.diag(R))
    norm_J = np.linalg.norm(J)
    return np.linalg.norm(R_off) / norm_J if norm_J > 0 else 0


def sample_subnetwork(rng, n_half=5):
    """Sample balanced subnetwork from C. elegans."""
    e_pick = rng.choice(exc_idx, n_half, replace=False)
    i_pick = rng.choice(inh_idx, n_half, replace=False)
    idx = np.concatenate([e_pick, i_pick])
    return W_norm_full[np.ix_(idx, idx)], signs_full[idx]


def random_dale_network(rng, n_total, signs, density):
    """Random network with Dale's law signs, matched density."""
    W = np.zeros((n_total, n_total))
    for i in range(n_total):
        for j in range(n_total):
            if i != j and rng.random() < max(density, 0.01):
                W[i, j] = signs[j] * rng.exponential(0.3)
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W


def degree_preserving_rewire(W, signs, n_swaps=None, rng=None):
    """Degree-preserving randomization: rewire edges keeping degree sequence.

    For each swap attempt: pick two edges (i->j) and (k->l), swap to
    (i->l) and (k->j), but only if the new edges don't already exist
    and signs are preserved (Dale's law).
    """
    if rng is None:
        rng = np.random.RandomState()
    W_new = W.copy()
    n = W.shape[0]

    # Find all existing edges
    edges = [(i, j) for i in range(n) for j in range(n)
             if i != j and W[i, j] != 0]
    if len(edges) < 2:
        return W_new

    if n_swaps is None:
        n_swaps = len(edges) * 10

    for _ in range(n_swaps):
        # Pick two random edges
        idx1, idx2 = rng.choice(len(edges), 2, replace=False)
        i, j = edges[idx1]
        k, l = edges[idx2]

        # Skip if same source or target
        if i == k or j == l or i == l or k == j:
            continue

        # Check new edges don't exist
        if W_new[i, l] != 0 or W_new[k, j] != 0:
            continue

        # Check Dale's law: source sign determines edge sign
        # Edge (i->j) has sign of neuron j. After swap: (i->l) has sign of l.
        # This changes the sign structure. Only swap if signs match.
        if signs[j] != signs[l]:
            continue

        # Swap
        W_new[i, l] = W_new[i, j]
        W_new[k, j] = W_new[k, l]
        W_new[i, j] = 0
        W_new[k, l] = 0

        # Update edge list
        edges[idx1] = (i, l)
        edges[idx2] = (k, j)

    return W_new


# ================================================================
# CHECK 1: Parameter sensitivity
# ================================================================
print("=" * 70)
print("CHECK 1: Parameter Sensitivity")
print("Is the 5-8x ratio robust across tau ratios and alpha values?")
print("=" * 70)

n_trials = 200
n_half = 5

print(f"\n{'tau_I/tau_E':>12s} {'alpha':>6s}  {'C.eleg':>8s}  {'Random':>8s}  {'Ratio':>7s}")
print("-" * 50)

for tau_ratio in [1.5, 2.0, 2.5, 3.0]:
    for alpha in [0.1, 0.3, 0.5]:
        tau_E = 10.0
        tau_I = tau_E * tau_ratio

        ce_res = []
        rand_res = []

        for trial in range(n_trials):
            rng = np.random.RandomState(trial + 100)
            W_sub, signs_sub = sample_subnetwork(rng, n_half)
            density = np.count_nonzero(W_sub) / (10 * 9)

            _, Q = build_swap(signs_sub)
            J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, alpha)
            ce_res.append(palindrome_residual(J_ce, Q))

            W_rand = random_dale_network(rng, 10, signs_sub, density)
            J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, alpha)
            rand_res.append(palindrome_residual(J_rand, Q))

        ce_m = np.mean(ce_res)
        rand_m = np.mean(rand_res)
        ratio = ce_m / rand_m if rand_m > 0 else 0
        print(f"  {tau_ratio:10.1f}  {alpha:5.1f}  {ce_m:8.4f}  {rand_m:8.4f}  {ratio:7.3f}")


# ================================================================
# CHECK 2: Pairing sensitivity
# ================================================================
print("\n" + "=" * 70)
print("CHECK 2: Pairing Sensitivity")
print("Does the E-I pairing choice affect the ratio?")
print("=" * 70)

tau_E, tau_I = 10.0, 20.0

for pairing_strategy in ['sequential', 'random_best']:
    ce_res = []
    rand_res = []

    for trial in range(n_trials):
        rng = np.random.RandomState(trial + 100)
        W_sub, signs_sub = sample_subnetwork(rng, n_half)
        density = np.count_nonzero(W_sub) / (10 * 9)

        if pairing_strategy == 'sequential':
            _, Q = build_swap(signs_sub, 'sequential')
            J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
            ce_res.append(palindrome_residual(J_ce, Q))

            W_rand = random_dale_network(rng, 10, signs_sub, density)
            J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, 0.3)
            rand_res.append(palindrome_residual(J_rand, Q))

        elif pairing_strategy == 'random_best':
            # Try 20 random pairings, take the best (lowest residual)
            J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
            best_ce = np.inf
            for _ in range(20):
                _, Q_try = build_swap(signs_sub, 'random')
                r = palindrome_residual(J_ce, Q_try)
                if r < best_ce:
                    best_ce = r
            ce_res.append(best_ce)

            W_rand = random_dale_network(rng, 10, signs_sub, density)
            J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, 0.3)
            best_rand = np.inf
            for _ in range(20):
                _, Q_try = build_swap(signs_sub, 'random')
                r = palindrome_residual(J_rand, Q_try)
                if r < best_rand:
                    best_rand = r
            rand_res.append(best_rand)

    ce_m = np.mean(ce_res)
    rand_m = np.mean(rand_res)
    ratio = ce_m / rand_m if rand_m > 0 else 0
    print(f"\n  {pairing_strategy:20s}:  C.eleg={ce_m:.4f}  Random={rand_m:.4f}  Ratio={ratio:.3f}")


# ================================================================
# CHECK 3: Degree-preserving null model
# ================================================================
print("\n" + "=" * 70)
print("CHECK 3: Degree-Preserving Null Model")
print("Sharper control: rewire edges but keep degree sequence")
print("=" * 70)

ce_res = []
dp_res = []

for trial in range(n_trials):
    rng = np.random.RandomState(trial + 100)
    W_sub, signs_sub = sample_subnetwork(rng, n_half)

    _, Q = build_swap(signs_sub)
    J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
    ce_res.append(palindrome_residual(J_ce, Q))

    # Degree-preserving rewiring
    W_dp = degree_preserving_rewire(W_sub, signs_sub, rng=rng)
    J_dp = build_jacobian(W_dp, tau_E, tau_I, signs_sub, 0.3)
    dp_res.append(palindrome_residual(J_dp, Q))

ce_m = np.mean(ce_res)
dp_m = np.mean(dp_res)
ratio = ce_m / dp_m if dp_m > 0 else 0

print(f"\n  C. elegans:              mean={ce_m:.4f}  std={np.std(ce_res):.4f}")
print(f"  Degree-preserving rand:  mean={dp_m:.4f}  std={np.std(dp_res):.4f}")
print(f"  Ratio (C.e./deg-pres):   {ratio:.3f}")

if ratio < 0.5:
    print(f"  >>> C. elegans STILL more palindromic (beyond degree sequence)")
elif ratio > 0.8 and ratio < 1.2:
    print(f"  >>> Comparable: degree sequence EXPLAINS the palindrome advantage")
else:
    print(f"  >>> Ratio {ratio:.2f}")


# ================================================================
# CHECK 4: Effect size and overlap
# ================================================================
print("\n" + "=" * 70)
print("CHECK 4: Effect Size and Distribution Overlap")
print("Is 5-8x within normal variation or clearly separated?")
print("=" * 70)

# Reuse check 1 data at default params
ce_res_full = []
rand_res_full = []
for trial in range(n_trials):
    rng = np.random.RandomState(trial + 100)
    W_sub, signs_sub = sample_subnetwork(rng, n_half)
    density = np.count_nonzero(W_sub) / (10 * 9)

    _, Q = build_swap(signs_sub)
    J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
    ce_res_full.append(palindrome_residual(J_ce, Q))

    W_rand = random_dale_network(rng, 10, signs_sub, density)
    J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, 0.3)
    rand_res_full.append(palindrome_residual(J_rand, Q))

ce_arr = np.array(ce_res_full)
rand_arr = np.array(rand_res_full)

# Cohen's d
pooled_std = np.sqrt((np.var(ce_arr) + np.var(rand_arr)) / 2)
cohens_d = (np.mean(rand_arr) - np.mean(ce_arr)) / pooled_std if pooled_std > 0 else 0

# Overlap: fraction of C. elegans samples that fall within random distribution range
ce_in_rand_range = np.sum((ce_arr >= np.percentile(rand_arr, 5)) &
                          (ce_arr <= np.percentile(rand_arr, 95)))
overlap_pct = ce_in_rand_range / len(ce_arr) * 100

# Rank-based: what fraction of random samples is C. elegans better than?
better_than = np.mean([np.mean(ce_arr[i] < rand_arr) for i in range(len(ce_arr))])

print(f"\n  C. elegans:  mean={np.mean(ce_arr):.4f}  "
      f"median={np.median(ce_arr):.4f}  "
      f"[{np.percentile(ce_arr, 5):.4f}, {np.percentile(ce_arr, 95):.4f}]")
print(f"  Random:      mean={np.mean(rand_arr):.4f}  "
      f"median={np.median(rand_arr):.4f}  "
      f"[{np.percentile(rand_arr, 5):.4f}, {np.percentile(rand_arr, 95):.4f}]")
print(f"\n  Cohen's d:   {cohens_d:.2f}  "
      f"({'large' if cohens_d > 0.8 else 'medium' if cohens_d > 0.5 else 'small'})")
print(f"  C.e. samples within random 5-95% range: {overlap_pct:.0f}%")
print(f"  C.e. better than random (rank-based):   {better_than*100:.0f}%")


# ================================================================
# SUMMARY
# ================================================================
print("\n" + "=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)
