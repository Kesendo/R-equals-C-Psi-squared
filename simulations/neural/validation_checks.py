#!/usr/bin/env python3
"""Instrument controls for the fitted off-diagonal residual.

This fitted diagonal condition does not test scalar-s F36. The random arm
uses its own maximum normalization, unlike the connectome's global scale.
On sparse blocks without Q-partner edges the residual sees the weight norm;
degree-preserving rewiring can therefore leave it unchanged. Neither those
ratios nor the distribution separation establish a biological advantage.
"""
import numpy as np

from celegans_trichotomy import load_worm
from neural_palindrome import fitted_offdiagonal_residual


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


def sample_subnetwork(rng, n_half, W_norm_full, signs_full):
    """Sample balanced subnetwork from the supplied C. elegans arrays."""
    exc_idx = np.flatnonzero(signs_full > 0)
    inh_idx = np.flatnonzero(signs_full < 0)
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


def main():
    W_norm_full, signs_full = load_worm()
    print("Instrument: fitted off-diagonal residual; this does not test scalar-s F36.")
    # ================================================================
    # CHECK 1: Parameter sensitivity
    # ================================================================
    print("=" * 70)
    print("CHECK 1: Parameter Sensitivity")
    print("Is the 5-8x ratio robust across tau ratios and alpha values?")
    print("WITHDRAWN 2026-08-26: the 5-8x is a difference of normalisation")
    print("constants, so its robustness is the robustness of an artifact.")
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
                W_sub, signs_sub = sample_subnetwork(rng, n_half, W_norm_full, signs_full)
                density = np.count_nonzero(W_sub) / (10 * 9)

                _, Q = build_swap(signs_sub)
                J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, alpha)
                ce_res.append(fitted_offdiagonal_residual(J_ce, np.argmax(Q, axis=1)))

                W_rand = random_dale_network(rng, 10, signs_sub, density)
                J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, alpha)
                rand_res.append(fitted_offdiagonal_residual(J_rand, np.argmax(Q, axis=1)))

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
    print("The ratio itself is WITHDRAWN 2026-08-26; what this sweep shows is")
    print("that the pairing choice is not the reason it was wrong.")
    print("=" * 70)

    tau_E, tau_I = 10.0, 20.0

    for pairing_strategy in ['sequential', 'random_best']:
        ce_res = []
        rand_res = []

        for trial in range(n_trials):
            rng = np.random.RandomState(trial + 100)
            W_sub, signs_sub = sample_subnetwork(rng, n_half, W_norm_full, signs_full)
            density = np.count_nonzero(W_sub) / (10 * 9)

            if pairing_strategy == 'sequential':
                _, Q = build_swap(signs_sub, 'sequential')
                J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
                ce_res.append(fitted_offdiagonal_residual(J_ce, np.argmax(Q, axis=1)))

                W_rand = random_dale_network(rng, 10, signs_sub, density)
                J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, 0.3)
                rand_res.append(fitted_offdiagonal_residual(J_rand, np.argmax(Q, axis=1)))

            elif pairing_strategy == 'random_best':
                # Try 20 random pairings, take the best (lowest residual)
                J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
                best_ce = np.inf
                for _ in range(20):
                    _, Q_try = build_swap(signs_sub, 'random')
                    r = fitted_offdiagonal_residual(J_ce, np.argmax(Q_try, axis=1))
                    if r < best_ce:
                        best_ce = r
                ce_res.append(best_ce)

                W_rand = random_dale_network(rng, 10, signs_sub, density)
                J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, 0.3)
                best_rand = np.inf
                for _ in range(20):
                    _, Q_try = build_swap(signs_sub, 'random')
                    r = fitted_offdiagonal_residual(J_rand, np.argmax(Q_try, axis=1))
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
    print("Rewire edges while keeping degree sequence; equality does not identify a cause.")
    print("=" * 70)

    ce_res = []
    dp_res = []

    for trial in range(n_trials):
        rng = np.random.RandomState(trial + 100)
        W_sub, signs_sub = sample_subnetwork(rng, n_half, W_norm_full, signs_full)

        _, Q = build_swap(signs_sub)
        J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
        ce_res.append(fitted_offdiagonal_residual(J_ce, np.argmax(Q, axis=1)))

        # Degree-preserving rewiring
        W_dp = degree_preserving_rewire(W_sub, signs_sub, rng=rng)
        J_dp = build_jacobian(W_dp, tau_E, tau_I, signs_sub, 0.3)
        dp_res.append(fitted_offdiagonal_residual(J_dp, np.argmax(Q, axis=1)))

    ce_m = np.mean(ce_res)
    dp_m = np.mean(dp_res)
    ratio = ce_m / dp_m if dp_m > 0 else 0

    print(f"\n  C. elegans:              mean={ce_m:.4f}  std={np.std(ce_res):.4f}")
    print(f"  Degree-preserving rand:  mean={dp_m:.4f}  std={np.std(dp_res):.4f}")
    print(f"  Ratio (C.e./deg-pres):   {ratio:.3f}")

    print(f"  >>> ratio {ratio:.2f}, and NO verdict is drawn from it: WITHDRAWN")
    print("      Without Q-partner edges the fitted residual only reads weight norms.")
    print("      Rewiring then leaves it unchanged unless it creates a Q-partner pair;")
    print("      on these sparse 5E+5I blocks valid swaps are also often absent.")
    print("      See docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md.")


    # ================================================================
    # CHECK 4: Effect size and overlap
    # ================================================================
    print("\n" + "=" * 70)
    print("CHECK 4: Effect Size and Distribution Overlap")
    print("Is 5-8x within normal variation or clearly separated?")
    print("WITHDRAWN 2026-08-26: the separation below is between two arms")
    print("normalised by different constants; it is not an effect size.")
    print("=" * 70)

    # Reuse check 1 data at default params
    ce_res_full = []
    rand_res_full = []
    for trial in range(n_trials):
        rng = np.random.RandomState(trial + 100)
        W_sub, signs_sub = sample_subnetwork(rng, n_half, W_norm_full, signs_full)
        density = np.count_nonzero(W_sub) / (10 * 9)

        _, Q = build_swap(signs_sub)
        J_ce = build_jacobian(W_sub, tau_E, tau_I, signs_sub, 0.3)
        ce_res_full.append(fitted_offdiagonal_residual(J_ce, np.argmax(Q, axis=1)))

        W_rand = random_dale_network(rng, 10, signs_sub, density)
        J_rand = build_jacobian(W_rand, tau_E, tau_I, signs_sub, 0.3)
        rand_res_full.append(fitted_offdiagonal_residual(J_rand, np.argmax(Q, axis=1)))

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
    print(f"  C.e. smaller fitted residual (rank-based):   {better_than*100:.0f}%")


    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)



if __name__ == "__main__":
    main()
