#!/usr/bin/env python3
"""
Q6: Does the 98.2% palindromic pairing require C. elegans topology?

Control experiment: generate random networks with the same parameters
(N=10, E=5, I=5, tau_E=10, tau_I=20, alpha=0.3) but random connectivity.
If random networks also show ~98%, the palindrome is GENERIC to any
balanced oscillatory network. If they show significantly less, the
C. elegans topology is essential.

Controls:
  A. Erdos-Renyi random graphs (same density as C. elegans subnetworks)
  B. Fully connected with random Gaussian weights
  C. Uniform weights (all connections equal strength)
  D. Unbalanced E/I (8E:2I) with random weights
  E. Uniform damping (tau_E = tau_I) with C. elegans topology
"""
import json
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# === Shared infrastructure ===

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
        if used[i]:
            continue
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


def run_ensemble(label, W_func, n_trials=200, **kw):
    """Run n_trials networks, report pairing statistics."""
    scores = []
    for trial in range(n_trials):
        W, s = W_func(seed=trial, **kw)
        tau_e = kw.get('tau_exc', 10.0)
        tau_i = kw.get('tau_inh', 20.0)
        alpha = kw.get('alpha', 0.3)
        J = build_jacobian(W, tau_e, tau_i, s, alpha)
        ev = np.linalg.eigvals(J)
        p, t, c = check_palindrome(ev)
        scores.append(p / t * 100 if t > 0 else 0)
    print(f"  {label:45s}  mean={np.mean(scores):5.1f}%  "
          f"std={np.std(scores):5.1f}%  "
          f"range=[{np.min(scores):.0f}%, {np.max(scores):.0f}%]")
    return scores


# === Load C. elegans baseline ===

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
W_norm_full = W_signed_full / max_w if max_w > 0 else W_signed_full

exc_idx = np.where(signs_full > 0)[0]
inh_idx = np.where(signs_full < 0)[0]


# === Network generators ===

def celegans_balanced_subnetwork(seed=0, **kw):
    """Random balanced subnetwork from C. elegans (the original test)."""
    rng = np.random.RandomState(seed + 42)
    e_pick = rng.choice(exc_idx, 5, replace=False)
    i_pick = rng.choice(inh_idx, 5, replace=False)
    idx = np.concatenate([e_pick, i_pick])
    W_sub = W_norm_full[np.ix_(idx, idx)]
    signs_sub = signs_full[idx]
    return W_sub, signs_sub


def erdos_renyi_balanced(seed=0, N=10, n_exc=5, density=0.5, **kw):
    """Erdos-Renyi random graph with balanced E/I."""
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    signs[n_exc:] = -1
    rng.shuffle(signs)

    # Random connectivity (directed)
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (N, N))
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if mask[i, j]:
                W[i, j] = signs[j] * weights[i, j]
    # Normalize
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def fully_connected_gaussian(seed=0, N=10, n_exc=5, **kw):
    """Fully connected network with Gaussian weights."""
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    signs[n_exc:] = -1
    rng.shuffle(signs)

    W_raw = rng.randn(N, N) * 0.3
    np.fill_diagonal(W_raw, 0)
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            W[i, j] = signs[j] * abs(W_raw[i, j])
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def uniform_weights_balanced(seed=0, N=10, n_exc=5, density=0.5, **kw):
    """Random graph with uniform (equal) nonzero weights."""
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    signs[n_exc:] = -1
    rng.shuffle(signs)

    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if mask[i, j]:
                W[i, j] = signs[j] * 1.0
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def small_world_balanced(seed=0, N=10, n_exc=5, k=4, p_rewire=0.3, **kw):
    """Watts-Strogatz small-world network with balanced E/I."""
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    signs[n_exc:] = -1
    rng.shuffle(signs)

    # Build ring lattice
    mask = np.zeros((N, N), dtype=bool)
    for i in range(N):
        for j in range(1, k // 2 + 1):
            mask[i, (i + j) % N] = True
            mask[i, (i - j) % N] = True

    # Rewire
    for i in range(N):
        for j in range(N):
            if mask[i, j] and rng.random() < p_rewire:
                mask[i, j] = False
                new_j = rng.randint(N)
                while new_j == i or mask[i, new_j]:
                    new_j = rng.randint(N)
                mask[i, new_j] = True

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


def erdos_renyi_unbalanced(seed=0, N=10, n_exc=8, density=0.5, **kw):
    """Erdos-Renyi with UNbalanced E/I (8E:2I like C. elegans ratio)."""
    return erdos_renyi_balanced(seed=seed, N=N, n_exc=n_exc,
                                density=density, **kw)


def celegans_uniform_damping(seed=0, **kw):
    """C. elegans topology but with UNIFORM damping (tau_E = tau_I)."""
    W, s = celegans_balanced_subnetwork(seed=seed)
    return W, s


# === Compute C. elegans subnetwork density for matched controls ===

densities = []
for trial in range(200):
    W, s = celegans_balanced_subnetwork(seed=trial)
    n_nonzero = np.count_nonzero(W) - np.count_nonzero(np.diag(W))
    densities.append(n_nonzero / (10 * 9))
mean_density = np.mean(densities)


# === Main ===

print("=" * 70)
print("Q6: Random Network Controls for Palindromic Pairing")
print("=" * 70)
print(f"\nC. elegans subnetwork density: {mean_density:.2f} (mean over 200 samples)")
print(f"Parameters: N=10, tau_E=10, tau_I=20, alpha=0.3")
print(f"Trials per condition: 200\n")

# Baseline: C. elegans balanced subnetworks
print("--- BASELINE ---")
scores_ce = run_ensemble(
    "C. elegans balanced subnetworks (5E:5I)",
    celegans_balanced_subnetwork)

# Control A: Erdos-Renyi, matched density
print("\n--- RANDOM TOPOLOGY CONTROLS (balanced E/I) ---")
scores_er = run_ensemble(
    "Erdos-Renyi (density matched)",
    erdos_renyi_balanced, density=mean_density)

# Control B: Fully connected Gaussian
scores_fc = run_ensemble(
    "Fully connected, Gaussian weights",
    fully_connected_gaussian)

# Control C: Uniform weights
scores_uw = run_ensemble(
    "Random graph, uniform weights",
    uniform_weights_balanced, density=mean_density)

# Control D: Small-world (Watts-Strogatz)
scores_sw = run_ensemble(
    "Watts-Strogatz small-world (k=4, p=0.3)",
    small_world_balanced)

# Control E: Unbalanced E/I
print("\n--- E/I BALANCE CONTROLS ---")
scores_ub = run_ensemble(
    "Erdos-Renyi UNBALANCED (8E:2I)",
    erdos_renyi_unbalanced, density=mean_density)

# Control F: Uniform damping (tau_E = tau_I)
scores_ud = run_ensemble(
    "C. elegans topology, UNIFORM damping (tau_E=tau_I)",
    celegans_uniform_damping, tau_exc=15.0, tau_inh=15.0)

# Control G: Stronger selectivity
print("\n--- DAMPING RATIO CONTROLS ---")
scores_s3 = run_ensemble(
    "Erdos-Renyi, tau_I/tau_E = 3.0",
    erdos_renyi_balanced, density=mean_density,
    tau_exc=10.0, tau_inh=30.0)

scores_s4 = run_ensemble(
    "Erdos-Renyi, tau_I/tau_E = 4.0",
    erdos_renyi_balanced, density=mean_density,
    tau_exc=10.0, tau_inh=40.0)

# === Verdict ===
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

ce_mean = np.mean(scores_ce)
er_mean = np.mean(scores_er)
fc_mean = np.mean(scores_fc)
uw_mean = np.mean(scores_uw)
sw_mean = np.mean(scores_sw)

print(f"""
  C. elegans balanced subnetworks:  {ce_mean:.1f}%
  Erdos-Renyi (density matched):    {er_mean:.1f}%
  Fully connected Gaussian:         {fc_mean:.1f}%
  Uniform weights random graph:     {uw_mean:.1f}%
  Small-world (Watts-Strogatz):     {sw_mean:.1f}%
""")

# === Density sweep ===
print("\n--- DENSITY SWEEP (Erdos-Renyi, balanced 5E:5I) ---")
print(f"{'density':>8s}  {'mean':>6s}  {'std':>6s}  {'range':>15s}")
for d in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
    sc = []
    for trial in range(200):
        W, s = erdos_renyi_balanced(seed=trial, density=d)
        J = build_jacobian(W, 10.0, 20.0, s, alpha=0.3)
        ev = np.linalg.eigvals(J)
        p, t, c = check_palindrome(ev)
        sc.append(p / t * 100 if t > 0 else 0)
    print(f"  {d:6.2f}    {np.mean(sc):5.1f}%  {np.std(sc):5.1f}%  "
          f"[{np.min(sc):.0f}%, {np.max(sc):.0f}%]")

# === Alpha (coupling strength) sweep at fixed density ===
print("\n--- COUPLING STRENGTH SWEEP (Erdos-Renyi, density=0.3, balanced) ---")
print(f"{'alpha':>8s}  {'mean':>6s}  {'std':>6s}  {'range':>15s}")
for a in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0]:
    sc = []
    for trial in range(200):
        W, s = erdos_renyi_balanced(seed=trial, density=0.3)
        J = build_jacobian(W, 10.0, 20.0, s, alpha=a)
        ev = np.linalg.eigvals(J)
        p, t, c = check_palindrome(ev)
        sc.append(p / t * 100 if t > 0 else 0)
    print(f"  {a:6.2f}    {np.mean(sc):5.1f}%  {np.std(sc):5.1f}%  "
          f"[{np.min(sc):.0f}%, {np.max(sc):.0f}%]")

if abs(er_mean - ce_mean) < 10 and abs(fc_mean - ce_mean) < 10:
    print("  >>> The palindromic pairing is GENERIC.")
    print("  >>> Any balanced oscillatory network shows it.")
    print("  >>> C. elegans topology is NOT required.")
    print("  >>> The mechanism is purely: selective damping (tau_E != tau_I).")
elif er_mean < ce_mean - 20:
    print("  >>> The palindromic pairing DEPENDS on topology.")
    print("  >>> C. elegans structure contributes beyond E/I balance.")
    print("  >>> Random networks show significantly less pairing.")
else:
    print("  >>> Intermediate result. Topology contributes but is not essential.")
    print("  >>> Detailed analysis needed.")
