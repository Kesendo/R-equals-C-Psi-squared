#!/usr/bin/env python3
"""
TASK: Inhibitory Position Effect on Palindromic Pairing (C. elegans)

Does WHERE the inhibitory neurons sit in the network affect palindromic
pairing? The qubit analogy predicts:
- I at edges (peripheral): highest pairing
- I at center: different character
- I random: lowest pairing

Usage: python celegans_inhibitory_position.py
"""
import json
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(SCRIPT_DIR, "celegans_connectome.json")) as f:
    data = json.load(f)

W_chem = np.array(data["chemical"])
signs = np.array(data["chemical_sign"])

with open(os.path.join(SCRIPT_DIR, "celegans_neuron_ids.txt")) as f:
    neuron_ids = [line.strip() for line in f.readlines()]

N = len(neuron_ids)

# Build signed weight matrix
W_signed = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        W_signed[i, j] = signs[j] * W_chem[j, i]

max_w = np.max(np.abs(W_signed))
W_norm = W_signed / max_w if max_w > 0 else W_signed


def build_jacobian(W, tau_exc, tau_inh, neuron_signs, alpha=0.3):
    n = len(neuron_signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_exc if neuron_signs[i] > 0 else tau_inh
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j and W[i, j] != 0:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def check_palindrome(eigenvalues, tol_frac=0.03):
    rates = sorted(-ev.real for ev in eigenvalues if abs(ev.real) > 1e-10)
    if len(rates) < 2:
        return 0, len(rates), 0
    center = (min(rates) + max(rates)) / 2.0
    tol = tol_frac * (max(rates) - min(rates))
    paired = 0
    used = [False] * len(rates)
    for i in range(len(rates)):
        if used[i]: continue
        target = 2 * center - rates[i]
        best_j, best_dist = -1, tol
        for j in range(len(rates)):
            if i != j and not used[j]:
                dist = abs(rates[j] - target)
                if dist < best_dist:
                    best_dist = dist
                    best_j = j
        if best_j >= 0:
            used[i] = True
            used[best_j] = True
            paired += 2
    return paired, len(rates), center


def degree_centrality(W_sub):
    """Total connection strength per neuron."""
    return np.sum(np.abs(W_sub), axis=0) + np.sum(np.abs(W_sub), axis=1)


# ============================================================
# Step 1: Reproduce the 80% baseline
# ============================================================
print("=" * 60)
print("STEP 1: Reproduce Phase 2 baseline (N=10, E=5, I=5)")
print("=" * 60)

connectivity = np.sum(np.abs(W_norm), axis=0) + np.sum(np.abs(W_norm), axis=1)
top10 = np.argsort(connectivity)[-10:].tolist()
W_sub = W_norm[np.ix_(top10, top10)]
signs_sub = signs[top10]
n_e = np.sum(signs_sub > 0)
n_i = np.sum(signs_sub < 0)

J = build_jacobian(W_sub, 10.0, 20.0, signs_sub)
ev = np.linalg.eigvals(J)
p, t, c = check_palindrome(ev)
print(f"Top-10 connected: E={n_e}, I={n_i}, pairing={p}/{t} = {p/t*100:.1f}%")

# ============================================================
# Step 2: Sample 200 balanced subnetworks, classify I-position
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Sample balanced subnetworks (N=10, E=5, I=5)")
print("=" * 60)

exc_idx = np.where(signs > 0)[0]
inh_idx = np.where(signs < 0)[0]

np.random.seed(42)
results = []

for trial in range(200):
    # Pick 5 random E and 5 random I neurons
    e_pick = np.random.choice(exc_idx, 5, replace=False)
    i_pick = np.random.choice(inh_idx, 5, replace=False)
    idx = np.concatenate([e_pick, i_pick])

    W_sub = W_norm[np.ix_(idx, idx)]
    signs_sub = signs[idx]

    # Compute centrality within subnetwork
    cent = degree_centrality(W_sub)

    # Mean centrality of I-neurons (indices 5-9 in subnetwork)
    mean_cent_i = np.mean(cent[5:])
    mean_cent_e = np.mean(cent[:5])
    i_centrality_ratio = mean_cent_i / (mean_cent_e + 1e-10)

    # Palindromic pairing
    J = build_jacobian(W_sub, 10.0, 20.0, signs_sub)
    ev = np.linalg.eigvals(J)
    p, t, c = check_palindrome(ev)
    score = p / t * 100 if t > 0 else 0

    results.append((score, i_centrality_ratio, mean_cent_i, mean_cent_e))

scores = [r[0] for r in results]
ratios = [r[1] for r in results]

print(f"Samples: {len(results)}")
print(f"Pairing: mean={np.mean(scores):.1f}%, std={np.std(scores):.1f}%, "
      f"min={np.min(scores):.0f}%, max={np.max(scores):.0f}%")

# Correlation
corr = np.corrcoef(scores, ratios)[0, 1]
print(f"\nCorrelation (pairing vs I-centrality ratio): r = {corr:.3f}")

# Bin by I-centrality
print("\nBinned by I-neuron centrality ratio:")
print(f"{'Bin':15s} {'Count':>5s} {'Mean Pairing':>12s} {'Std':>8s}")
sorted_results = sorted(results, key=lambda x: x[1])
bin_size = len(sorted_results) // 4
for i, label in enumerate(["I-peripheral", "I-moderate-low", "I-moderate-high", "I-central"]):
    chunk = sorted_results[i*bin_size:(i+1)*bin_size]
    s = [c[0] for c in chunk]
    r = [c[1] for c in chunk]
    print(f"{label:15s} {len(chunk):5d} {np.mean(s):12.1f}% {np.std(s):7.1f}%  "
          f"(ratio {np.mean(r):.2f})")


# ============================================================
# Step 3: Direct sacrifice-zone test on one subnetwork
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Fix neurons, vary I-assignment by position")
print("=" * 60)

# Pick 10 well-connected neurons (mixed E/I from original)
top20 = np.argsort(connectivity)[-20:]
# Use all of them, assign I based on position
W_20 = W_norm[np.ix_(top20, top20)]
cent_20 = degree_centrality(W_20)
sorted_by_cent = np.argsort(cent_20)  # low centrality = peripheral

print(f"\nUsing top-20 connected neurons, assigning 10 as E and 10 as I:")

for label, i_indices in [
    ("I-PERIPHERAL (5 least central as I)",
     sorted_by_cent[:5]),
    ("I-CENTRAL (5 most central as I)",
     sorted_by_cent[-5:]),
    ("I-MIXED (5 medium central as I)",
     sorted_by_cent[7:12]),
]:
    # Assign signs: +1 for E, -1 for I
    test_signs = np.ones(20)
    test_signs[i_indices] = -1

    n_e = np.sum(test_signs > 0)
    n_i = np.sum(test_signs < 0)

    J = build_jacobian(W_20, 10.0, 20.0, test_signs)
    ev = np.linalg.eigvals(J)
    p, t, c = check_palindrome(ev)
    score = p / t * 100 if t > 0 else 0
    print(f"  {label}: E={n_e:.0f} I={n_i:.0f} -> {p}/{t} = {score:.1f}%")

# Also try N=10 version
top10_neurons = np.argsort(connectivity)[-10:]
W_10 = W_norm[np.ix_(top10_neurons, top10_neurons)]
cent_10 = degree_centrality(W_10)
sorted_10 = np.argsort(cent_10)

print(f"\nUsing top-10 connected neurons, assigning 5 as E and 5 as I:")

for label, i_indices in [
    ("I-PERIPHERAL", sorted_10[:5]),
    ("I-CENTRAL", sorted_10[-5:]),
]:
    test_signs = np.ones(10)
    test_signs[i_indices] = -1

    J = build_jacobian(W_10, 10.0, 20.0, test_signs)
    ev = np.linalg.eigvals(J)
    p, t, c = check_palindrome(ev)
    score = p / t * 100 if t > 0 else 0
    print(f"  {label}: {p}/{t} = {score:.1f}%")

# Random assignment baseline (20 trials)
random_scores = []
for _ in range(20):
    perm = np.random.permutation(10)
    test_signs = np.ones(10)
    test_signs[perm[:5]] = -1
    J = build_jacobian(W_10, 10.0, 20.0, test_signs)
    ev = np.linalg.eigvals(J)
    p, t, c = check_palindrome(ev)
    random_scores.append(p / t * 100 if t > 0 else 0)

print(f"  I-RANDOM (20 trials): mean={np.mean(random_scores):.1f}%, "
      f"std={np.std(random_scores):.1f}%, range=[{np.min(random_scores):.0f}%, {np.max(random_scores):.0f}%]")

print("\nDone.")
