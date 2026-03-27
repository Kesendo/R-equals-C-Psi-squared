#!/usr/bin/env python3
"""
Algebraic Palindrome Condition for Wilson-Cowan

NOT tolerance-based eigenvalue matching. ALGEBRAIC residual:
  R = Q*J*Q + J + 2*S
where Q is the E-I swap permutation.

If R = 0: exact palindrome. ||R|| measures deviation.

Decomposed:
  R_self = 0 (exact by construction, self-decay swaps tau_E <-> tau_I)
  R_coupling = T'*W' + T*W (antisymmetry residual of coupling)

For exact palindrome, coupling must satisfy:
  W[Q(i),Q(j)] = -(tau_{Q(i)}/tau_i) * W[i,j]

Dale's law gives the signs. The magnitude condition for tau_I/tau_E = 2:
  |W[I,I]| = 2 * |W[E,E]|
  |W[I,E]| = 2 * |W[E,I]|
"""
import json
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# === Load C. elegans ===
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


def build_swap(signs):
    """Build E-I swap permutation. Pair E neurons with I neurons.

    Returns permutation array and permutation matrix.
    """
    n = len(signs)
    e_local = np.where(signs > 0)[0]
    i_local = np.where(signs < 0)[0]
    n_pairs = min(len(e_local), len(i_local))

    perm = np.arange(n)
    for k in range(n_pairs):
        perm[e_local[k]] = i_local[k]
        perm[i_local[k]] = e_local[k]

    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, perm[i]] = 1.0

    return perm, Q


def palindrome_residual(J, Q, tau_E, tau_I, signs):
    """Compute palindrome residual R = Q*J*Q + J + 2*S.

    Returns total residual, self-decay residual, coupling residual.
    """
    n = len(signs)
    # S_k = (1/tau_E - 1/tau_I) / 2 for all k
    # Actually S depends on which neuron: for E-neuron, we need
    # Q*D*Q diagonal to match. Let me compute directly.

    QJQ = Q @ J @ Q.T
    # 2S: chosen so that diagonal part is zero
    # QJQ[i,i] + J[i,i] + 2*S[i] = 0
    # S[i] = -(QJQ[i,i] + J[i,i]) / 2
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    S = np.diag(S_diag)

    R = QJQ + J + 2 * S
    R_diag = np.diag(np.diag(R))
    R_offdiag = R - R_diag

    norm_J = np.linalg.norm(J)
    return (np.linalg.norm(R) / norm_J,
            np.linalg.norm(R_diag) / norm_J,
            np.linalg.norm(R_offdiag) / norm_J,
            R)


def coupling_block_analysis(W, signs, perm, tau_E, tau_I):
    """Analyze the antisymmetry condition by connection type.

    Condition: W[Q(i),Q(j)] = -(tau_{Q(i)}/tau_i) * W[i,j]
    Decompose by type: EE, EI, IE, II.
    """
    n = len(signs)
    blocks = {'EE': [], 'EI': [], 'IE': [], 'II': []}
    residuals = {'EE': [], 'EI': [], 'IE': [], 'II': []}

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            qi, qj = perm[i], perm[j]
            tau_i = tau_E if signs[i] > 0 else tau_I
            tau_qi = tau_E if signs[qi] > 0 else tau_I

            # Connection type
            t_i = 'E' if signs[i] > 0 else 'I'
            t_j = 'E' if signs[j] > 0 else 'I'
            block = t_i + t_j

            w_ij = W[i, j]
            w_qi_qj = W[qi, qj]
            predicted = -(tau_qi / tau_i) * w_ij
            residual = w_qi_qj - predicted

            if abs(w_ij) > 1e-10:
                blocks[block].append((w_ij, w_qi_qj, predicted))
                residuals[block].append(residual)

    return blocks, residuals


def magnitude_ratio_test(W, signs, perm, tau_E, tau_I):
    """Test the magnitude condition: |W[Q(i),Q(j)]| / |W[i,j]| = tau_{Q(i)}/tau_i.

    For EE connections: ratio should be tau_I/tau_E = 2.0
    For II connections: ratio should be tau_E/tau_I = 0.5
    For EI connections: ratio should be tau_I/tau_I = 1.0 (I stays I)...

    Wait, need to be more careful about what Q does to each pair.
    """
    n = len(signs)
    ratios_by_type = {'EE': [], 'EI': [], 'IE': [], 'II': []}

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            qi, qj = perm[i], perm[j]
            tau_i = tau_E if signs[i] > 0 else tau_I
            tau_qi = tau_E if signs[qi] > 0 else tau_I
            expected_ratio = tau_qi / tau_i

            w_ij = W[i, j]
            w_qi_qj = W[qi, qj]

            if abs(w_ij) > 1e-10:
                actual_ratio = abs(w_qi_qj / w_ij) if abs(w_ij) > 1e-10 else 0
                t_i = 'E' if signs[i] > 0 else 'I'
                t_j = 'E' if signs[j] > 0 else 'I'
                ratios_by_type[t_i + t_j].append(
                    (actual_ratio, expected_ratio))

    return ratios_by_type


# === Phase 1: Synthetic networks ===
print("=" * 70)
print("PHASE 1: Synthetic Networks - Does Dale's Law Give Exact Palindrome?")
print("=" * 70)

tau_E, tau_I = 10.0, 20.0
N = 10
n_half = 5

np.random.seed(42)
signs_syn = np.ones(N)
signs_syn[n_half:] = -1

# (a) Perfect Dale's law with symmetric magnitudes (scaled by tau ratio)
print("\n(a) Perfect Dale's law + magnitude condition:")
W_perfect = np.zeros((N, N))
perm_syn, Q_syn = build_swap(signs_syn)
for i in range(N):
    for j in range(N):
        if i == j:
            continue
        qi, qj = perm_syn[i], perm_syn[j]
        tau_i = tau_E if signs_syn[i] > 0 else tau_I
        tau_qi = tau_E if signs_syn[qi] > 0 else tau_I
        if i < j:  # set upper triangle, derive lower from condition
            base = np.random.randn() * 0.3
            W_perfect[i, j] = signs_syn[j] * abs(base)
            W_perfect[qi, qj] = -(tau_qi / tau_i) * W_perfect[i, j]
mx = np.max(np.abs(W_perfect))
if mx > 0:
    W_perfect /= mx

J_perfect = build_jacobian(W_perfect, tau_E, tau_I, signs_syn, alpha=1.0)
r_tot, r_diag, r_off, R = palindrome_residual(J_perfect, Q_syn, tau_E, tau_I, signs_syn)
print(f"  ||R||/||J||     = {r_tot:.2e}")
print(f"  ||R_diag||/||J|| = {r_diag:.2e}")
print(f"  ||R_off||/||J||  = {r_off:.2e}")
if r_tot < 1e-10:
    print(f"  >>> EXACT PALINDROME (R = 0 to machine precision)")

# (b) Dale's law signs, random magnitudes
print("\n(b) Dale's law signs, random magnitudes:")
W_dale = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if i == j:
            continue
        W_dale[i, j] = signs_syn[j] * abs(np.random.randn() * 0.3)
mx = np.max(np.abs(W_dale))
if mx > 0:
    W_dale /= mx

J_dale = build_jacobian(W_dale, tau_E, tau_I, signs_syn, alpha=1.0)
r_tot_d, r_diag_d, r_off_d, _ = palindrome_residual(J_dale, Q_syn, tau_E, tau_I, signs_syn)
print(f"  ||R||/||J||     = {r_tot_d:.4f}")
print(f"  ||R_off||/||J||  = {r_off_d:.4f}")

# (c) Random signs and magnitudes
print("\n(c) Random signs, random magnitudes:")
W_random = np.random.randn(N, N) * 0.3
np.fill_diagonal(W_random, 0)
mx = np.max(np.abs(W_random))
W_random /= mx

J_random = build_jacobian(W_random, tau_E, tau_I, signs_syn, alpha=1.0)
r_tot_r, r_diag_r, r_off_r, _ = palindrome_residual(J_random, Q_syn, tau_E, tau_I, signs_syn)
print(f"  ||R||/||J||     = {r_tot_r:.4f}")
print(f"  ||R_off||/||J||  = {r_off_r:.4f}")


# === Phase 2: C. elegans subnetworks ===
print("\n" + "=" * 70)
print("PHASE 2: C. elegans Balanced Subnetworks vs Random")
print("=" * 70)

for n_half in [5, 10, 13]:
    n_total = 2 * n_half
    if n_half > len(inh_idx):
        continue

    print(f"\n--- N={n_total} ({n_half}E + {n_half}I) ---")

    ce_residuals = []
    dale_residuals = []

    for trial in range(200):
        rng = np.random.RandomState(trial + 100)

        # C. elegans subnetwork
        e_pick = rng.choice(exc_idx, n_half, replace=False)
        i_pick = rng.choice(inh_idx, n_half, replace=False)
        idx = np.concatenate([e_pick, i_pick])
        W_sub = W_norm_full[np.ix_(idx, idx)]
        signs_sub = signs_full[idx]

        perm_sub, Q_sub = build_swap(signs_sub)
        J_sub = build_jacobian(W_sub, tau_E, tau_I, signs_sub, alpha=0.3)
        r_tot, _, r_off, _ = palindrome_residual(
            J_sub, Q_sub, tau_E, tau_I, signs_sub)
        ce_residuals.append(r_off)

        # Random Dale's law control (same density)
        density = np.count_nonzero(W_sub) / (n_total * (n_total - 1))
        W_rand_dale = np.zeros((n_total, n_total))
        for i in range(n_total):
            for j in range(n_total):
                if i != j and rng.random() < max(density, 0.01):
                    W_rand_dale[i, j] = signs_sub[j] * rng.exponential(0.3)
        mx = np.max(np.abs(W_rand_dale))
        if mx > 0:
            W_rand_dale /= mx

        J_rand = build_jacobian(W_rand_dale, tau_E, tau_I, signs_sub, alpha=0.3)
        r_tot_r, _, r_off_r, _ = palindrome_residual(
            J_rand, Q_sub, tau_E, tau_I, signs_sub)
        dale_residuals.append(r_off_r)

    ce_mean = np.mean(ce_residuals)
    dale_mean = np.mean(dale_residuals)
    ratio = ce_mean / dale_mean if dale_mean > 0 else 0

    print(f"  Coupling residual ||R_off||/||J|| (lower = more palindromic):")
    print(f"    C. elegans:       mean={ce_mean:.6f}  std={np.std(ce_residuals):.6f}")
    print(f"    Random (Dale):    mean={dale_mean:.6f}  std={np.std(dale_residuals):.6f}")
    print(f"    Ratio (C.e./random): {ratio:.3f}")

    if ratio < 0.8:
        print(f"    >>> C. elegans is SIGNIFICANTLY more palindromic")
    elif ratio > 1.2:
        print(f"    >>> C. elegans is LESS palindromic than random Dale's law")
    else:
        print(f"    >>> Comparable (ratio {ratio:.2f})")


# === Phase 3: Magnitude ratio test on C. elegans ===
print("\n" + "=" * 70)
print("PHASE 3: Magnitude Ratio Test")
print("=" * 70)
print(f"Predicted ratio for tau_I/tau_E = {tau_I/tau_E:.1f}:")
print(f"  EE connections: |W[I-partner,I-partner]| / |W[E,E]| = {tau_I/tau_E:.1f}")
print(f"  II connections: |W[E-partner,E-partner]| / |W[I,I]| = {tau_E/tau_I:.1f}")
print(f"  EI connections: |W[I-partner,E-partner]| / |W[E,I]| = varies")

# Collect magnitude ratios across many subnetworks
all_ratios = {'EE': [], 'EI': [], 'IE': [], 'II': []}
for trial in range(200):
    rng = np.random.RandomState(trial + 200)
    e_pick = rng.choice(exc_idx, 5, replace=False)
    i_pick = rng.choice(inh_idx, 5, replace=False)
    idx = np.concatenate([e_pick, i_pick])
    W_sub = W_norm_full[np.ix_(idx, idx)]
    signs_sub = signs_full[idx]
    perm_sub, _ = build_swap(signs_sub)

    ratios = magnitude_ratio_test(W_sub, signs_sub, perm_sub, tau_E, tau_I)
    for block in ratios:
        for actual, expected in ratios[block]:
            all_ratios[block].append((actual, expected))

print(f"\n  Observed magnitude ratios (from 200 subnetworks):")
for block in ['EE', 'EI', 'IE', 'II']:
    if all_ratios[block]:
        actuals = [a for a, e in all_ratios[block]]
        expecteds = [e for a, e in all_ratios[block]]
        expected_val = np.mean(expecteds)
        print(f"    {block}: observed={np.mean(actuals):.3f} +/- {np.std(actuals):.3f}"
              f"  predicted={expected_val:.1f}"
              f"  n={len(actuals)}")
    else:
        print(f"    {block}: no data")

# === Summary ===
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
The algebraic palindrome condition Q*J*Q + J + 2S = 0 requires:
  1. Selective damping (tau_E != tau_I)     -- UNIVERSAL in biology
  2. Dale's law sign structure               -- UNIVERSAL in neurons
  3. Magnitude condition: |W[Q(i),Q(j)]| = (tau_Q(i)/tau_i) * |W[i,j]|

Condition 3 predicts specific coupling strength ratios that depend
on tau_I/tau_E. For tau_I/tau_E = {tau_I/tau_E:.1f}:
  |inhibition-to-inhibition| = {tau_I/tau_E:.1f}x |excitation-to-excitation|

This is testable on any connectome with known E/I labels and
synaptic weights.
""")
