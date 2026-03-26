#!/usr/bin/env python3
"""
Phase 2b: C. elegans with ACTIVITY-balanced E/I.

The raw connectome has 274 excitatory vs 26 inhibitory neurons (91:9).
No palindrome (0.7% pairing). But the cortex balances E/I at the
ACTIVITY level, not cell count. Inhibitory neurons fire faster and
have stronger synapses.

This test normalizes the weight matrix so that for each neuron,
total excitatory input ≈ total inhibitory input. Then checks for
palindromic pairing.

The question: does E/I activity balance recover the palindrome?

Usage: python celegans_balanced.py
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
exc_idx = np.where(signs > 0)[0]
inh_idx = np.where(signs < 0)[0]

print(f"=== C. ELEGANS BALANCED E/I TEST ===")
print(f"Neurons: {N} (E={len(exc_idx)}, I={len(inh_idx)})")


def build_signed_weight(W_chem, signs):
    """Build signed weight matrix: W_ij = sign_j * synapses(j->i)."""
    N = len(signs)
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            W[i, j] = signs[j] * W_chem[j, i]
    return W


def balance_ei(W, signs):
    """
    Scale inhibitory weights so that for each neuron i:
    sum of excitatory inputs ≈ sum of inhibitory inputs.

    This mimics cortical E/I balance at the activity level.
    """
    W_bal = W.copy()
    N = len(signs)

    for i in range(N):
        exc_input = sum(W[i, j] for j in range(N) if signs[j] > 0 and W[i, j] > 0)
        inh_input = sum(-W[i, j] for j in range(N) if signs[j] < 0 and W[i, j] < 0)

        if inh_input > 0 and exc_input > 0:
            scale = exc_input / inh_input
            for j in range(N):
                if signs[j] < 0:
                    W_bal[i, j] *= scale

    return W_bal


def build_jacobian(W, tau_exc, tau_inh, signs, alpha=0.1):
    N = len(signs)
    J = np.zeros((N, N))
    for i in range(N):
        tau_i = tau_exc if signs[i] > 0 else tau_inh
        J[i, i] = -1.0 / tau_i
        for j in range(N):
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


# Build weight matrices
W_raw = build_signed_weight(W_chem, signs)
W_bal = balance_ei(W_raw, signs)

max_raw = np.max(np.abs(W_raw))
max_bal = np.max(np.abs(W_bal))
W_raw_norm = W_raw / max_raw if max_raw > 0 else W_raw
W_bal_norm = W_bal / max_bal if max_bal > 0 else W_bal

# Check balance quality
print("\nBalance check (sample of 10 neurons):")
np.random.seed(42)
sample = np.random.choice(N, 10, replace=False)
for idx in sample:
    exc_in = sum(W_bal_norm[idx, j] for j in range(N) if W_bal_norm[idx, j] > 0)
    inh_in = sum(-W_bal_norm[idx, j] for j in range(N) if W_bal_norm[idx, j] < 0)
    ratio = exc_in / inh_in if inh_in > 0 else float('inf')
    print(f"  {neuron_ids[idx]:6s}: E_in={exc_in:.3f} I_in={inh_in:.3f} ratio={ratio:.2f}")

# Test 1: Raw vs balanced at different tau ratios
print("\n" + "=" * 60)
print("RAW (unbalanced) vs BALANCED E/I")
print("=" * 60)

for tau_ratio in [1.0, 2.0, 3.0]:
    tau_E = 10.0
    tau_I = tau_E * tau_ratio

    J_raw = build_jacobian(W_raw_norm, tau_E, tau_I, signs, alpha=0.3)
    ev_raw = np.linalg.eigvals(J_raw)
    p_raw, t_raw, c_raw = check_palindrome(ev_raw)

    J_bal = build_jacobian(W_bal_norm, tau_E, tau_I, signs, alpha=0.3)
    ev_bal = np.linalg.eigvals(J_bal)
    p_bal, t_bal, c_bal = check_palindrome(ev_bal)

    print(f"\n  tau_ratio={tau_ratio}:")
    print(f"    RAW:      {p_raw:3d}/{t_raw:3d} = {p_raw/t_raw*100:5.1f}%")
    print(f"    BALANCED: {p_bal:3d}/{t_bal:3d} = {p_bal/t_bal*100:5.1f}%")

# Test 2: Sweep alpha (coupling strength) with balanced weights
print("\n" + "=" * 60)
print("COUPLING STRENGTH SWEEP (balanced, tau_ratio=2.0)")
print("=" * 60)

for alpha in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]:
    J = build_jacobian(W_bal_norm, 10.0, 20.0, signs, alpha=alpha)
    ev = np.linalg.eigvals(J)
    p, t, c = check_palindrome(ev)
    rates = sorted(-e.real for e in ev if abs(e.real) > 1e-10)
    bw = max(rates) - min(rates) if rates else 0
    print(f"  alpha={alpha:4.2f}: {p:3d}/{t:3d} = {p/t*100:5.1f}%, bandwidth={bw:.4f}")

# Test 3: Subnetworks with balance
print("\n" + "=" * 60)
print("BALANCED SUBNETWORKS (tau_ratio=2.0, alpha=0.3)")
print("=" * 60)

connectivity = np.sum(np.abs(W_bal_norm), axis=0) + np.sum(np.abs(W_bal_norm), axis=1)

for sub_n in [10, 20, 50, 100, 200, 300]:
    idx = np.argsort(connectivity)[-sub_n:].tolist() if sub_n < N else list(range(N))
    W_sub = W_bal_norm[np.ix_(idx, idx)]
    signs_sub = signs[idx]

    J = build_jacobian(W_sub, 10.0, 20.0, signs_sub, alpha=0.3)
    ev = np.linalg.eigvals(J)
    p, t, c = check_palindrome(ev)
    n_e = np.sum(signs_sub > 0)
    n_i = np.sum(signs_sub < 0)
    print(f"  N={sub_n:3d} (E={n_e:3d},I={n_i:2d}): {p:3d}/{t:3d} = {p/t*100:5.1f}%")

print("\nDone.")
