#!/usr/bin/env python3
"""
Phase 2: C. elegans connectome palindrome test.

300 neurons, real synaptic weights, real E/I classification.
Build Wilson-Cowan-style Jacobian from the connectome and test
for palindromic eigenvalue pairing.

Source: WormNeuroAtlas (francescorandi/wormneuroatlas)
Data: Cook et al. 2019 updated connectome

Usage: python celegans_palindrome.py
"""
import json
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Load connectome
with open(os.path.join(SCRIPT_DIR, "celegans_connectome.json")) as f:
    data = json.load(f)

W_chem = np.array(data["chemical"])       # 300x300 directed synapse counts
W_gap = np.array(data["electrical"])       # 300x300 symmetric gap junctions
signs = np.array(data["chemical_sign"])    # +1 excitatory, -1 inhibitory

# Load neuron IDs
with open(os.path.join(SCRIPT_DIR, "celegans_neuron_ids.txt")) as f:
    neuron_ids = [line.strip() for line in f.readlines()]

N = len(neuron_ids)
n_exc = np.sum(signs > 0)
n_inh = np.sum(signs < 0)

print(f"=== C. ELEGANS CONNECTOME PALINDROME TEST ===")
print(f"Neurons: {N}")
print(f"Excitatory: {n_exc}, Inhibitory: {n_inh}")
print(f"Chemical synapses (nonzero): {np.count_nonzero(W_chem)}")
print(f"Gap junctions (nonzero): {np.count_nonzero(W_gap)}")
print()

# Build signed weight matrix: W_ij = sign_j * synapses_ji
# (presynaptic j with sign, postsynaptic i)
W_signed = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        W_signed[i, j] = signs[j] * W_chem[j, i]  # j -> i, signed by j's type

# Normalize weights
max_w = np.max(np.abs(W_signed))
if max_w > 0:
    W_norm = W_signed / max_w


def build_jacobian(W, tau_exc, tau_inh, signs, alpha=0.1):
    """
    Build linearized dynamics matrix.

    dX_i/dt = (-X_i + alpha * sum_j W_ij X_j) / tau_i

    where tau_i = tau_exc if neuron i is excitatory, tau_inh if inhibitory.
    alpha scales the effective coupling strength.
    """
    N = len(signs)
    J = np.zeros((N, N))

    for i in range(N):
        tau_i = tau_exc if signs[i] > 0 else tau_inh

        # Self-decay
        J[i, i] = -1.0 / tau_i

        # Coupling from other neurons
        for j in range(N):
            if i != j and W[i, j] != 0:
                J[i, j] = alpha * W[i, j] / tau_i

    return J


def check_palindrome(eigenvalues, tol_frac=0.03):
    """Check palindromic pairing of decay rates."""
    rates = sorted(-ev.real for ev in eigenvalues if abs(ev.real) > 1e-10)

    if len(rates) < 2:
        return 0, len(rates), 0

    center = (min(rates) + max(rates)) / 2.0
    tol = tol_frac * (max(rates) - min(rates))

    paired = 0
    used = [False] * len(rates)

    for i in range(len(rates)):
        if used[i]:
            continue
        target = 2 * center - rates[i]
        best_j = -1
        best_dist = tol
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


# Test 1: Full connectome with selective damping
print("=" * 60)
print("TEST 1: Full C. elegans (300 neurons), selective damping")
print("=" * 60)

for tau_ratio in [1.0, 1.5, 2.0, 3.0, 5.0]:
    tau_exc = 10.0
    tau_inh = tau_exc * tau_ratio

    J = build_jacobian(W_norm, tau_exc, tau_inh, signs, alpha=0.3)
    eigenvalues = np.linalg.eigvals(J)

    paired, total, center = check_palindrome(eigenvalues)
    score = paired / total * 100 if total > 0 else 0

    # Rate distribution
    rates = sorted(-ev.real for ev in eigenvalues if abs(ev.real) > 1e-10)
    r_min, r_max = min(rates), max(rates)

    print(f"\n  tau_E={tau_exc}, tau_I={tau_inh} (ratio={tau_ratio})")
    print(f"  Rates: {total} distinct, range [{r_min:.4f}, {r_max:.4f}]")
    print(f"  Center: {center:.4f}")
    print(f"  Paired: {paired}/{total} = {score:.1f}%")


# Test 2: Subnetwork - motor neurons only
print("\n\n" + "=" * 60)
print("TEST 2: Motor neuron subnetwork")
print("=" * 60)

# Simple heuristic: neurons with names starting with typical motor neuron prefixes
motor_prefixes = ['DA', 'DB', 'DD', 'VA', 'VB', 'VD', 'AS']
motor_idx = [i for i, name in enumerate(neuron_ids)
             if any(name.startswith(p) for p in motor_prefixes)]
print(f"Motor neurons found: {len(motor_idx)}")

if len(motor_idx) > 5:
    W_motor = W_norm[np.ix_(motor_idx, motor_idx)]
    signs_motor = signs[motor_idx]

    for tau_ratio in [1.0, 2.0, 3.0]:
        tau_exc = 10.0
        tau_inh = tau_exc * tau_ratio
        J = build_jacobian(W_motor, tau_exc, tau_inh, signs_motor, alpha=0.3)
        eigenvalues = np.linalg.eigvals(J)
        paired, total, center = check_palindrome(eigenvalues)
        score = paired / total * 100 if total > 0 else 0
        print(f"  ratio={tau_ratio}: {paired}/{total} = {score:.1f}%")


# Test 3: Small subnetworks (N=10, 20, 50) for scaling
print("\n\n" + "=" * 60)
print("TEST 3: Scaling with subnetwork size")
print("=" * 60)

np.random.seed(42)
tau_ratio = 2.0
tau_exc = 10.0
tau_inh = tau_exc * tau_ratio

for sub_n in [10, 20, 50, 100, 200, 300]:
    if sub_n >= N:
        idx = list(range(N))
    else:
        # Take most connected neurons
        connectivity = np.sum(np.abs(W_norm), axis=0) + np.sum(np.abs(W_norm), axis=1)
        idx = np.argsort(connectivity)[-sub_n:].tolist()

    W_sub = W_norm[np.ix_(idx, idx)]
    signs_sub = signs[idx]

    J = build_jacobian(W_sub, tau_exc, tau_inh, signs_sub, alpha=0.3)
    eigenvalues = np.linalg.eigvals(J)
    paired, total, center = check_palindrome(eigenvalues)
    score = paired / total * 100 if total > 0 else 0

    rates = sorted(-ev.real for ev in eigenvalues if abs(ev.real) > 1e-10)
    bandwidth = max(rates) - min(rates) if rates else 0

    print(f"  N={sub_n:3d}: {paired:3d}/{total:3d} = {score:5.1f}% paired, "
          f"bandwidth={bandwidth:.4f}, center={center:.4f}")

print("\n\nDone.")
