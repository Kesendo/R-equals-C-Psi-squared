#!/usr/bin/env python3
"""Thresholded scalar-center and complex spectral readings of subcircuits.

The labels truly/soft/hard are numerical bins, not exact F36 verdicts.
The scalar residual uses the declared s=(1/tau_E+1/tau_I)/2 and a chosen
E/I swap; complex spectral assignment preserves multiplicities. Controls
use unequal weight normalizations, so biological comparisons are not valid.
Equal worm/degree-preserved scores do not establish degree-distribution
causality: the withdrawn fitted metric is insensitive to these rewires.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from neural_palindrome import scalar_center_residual, spectral_pairing_error

SCRIPT_DIR = Path(__file__).parent


def load_worm():
    """Load C. elegans connectome and signed weight matrix."""
    with open(SCRIPT_DIR / 'celegans_connectome.json') as f:
        data = json.load(f)
    W_chem = np.array(data['chemical'])
    signs_full = np.array(data['chemical_sign'])
    N_full = len(signs_full)

    W_signed_full = np.zeros((N_full, N_full))
    for i in range(N_full):
        for j in range(N_full):
            W_signed_full[i, j] = signs_full[j] * W_chem[j, i]
    max_w = np.max(np.abs(W_signed_full))
    W_norm_full = W_signed_full / max_w

    return W_norm_full, signs_full


def build_jacobian(W, tau_E, tau_I, signs, alpha=0.3):
    """Wilson-Cowan-style linearized dynamics matrix."""
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
    """Build E-I swap permutation matrix Q with Q^2 = I."""
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


def classify_subcircuit(W, signs, tau_E, tau_I, alpha=0.3,
                        truly_threshold=0.01, pairing_threshold=0.05):
    """Return ('truly'|'soft'|'hard', residual_norm, max_pair_deviation)."""
    perm, _ = build_swap(signs)
    J = build_jacobian(W, tau_E, tau_I, signs, alpha=alpha)
    s = 0.5 * (1 / tau_E + 1 / tau_I)
    r_tot = scalar_center_residual(J, perm, s)
    values = np.linalg.eigvals(J)
    scale = float(np.max(np.abs(values)))
    pair_dev = spectral_pairing_error(values, s) / (scale if scale else 1.0)

    if r_tot < truly_threshold:
        return 'truly', r_tot, pair_dev
    if pair_dev < pairing_threshold:
        return 'soft', r_tot, pair_dev
    return 'hard', r_tot, pair_dev


def degree_preserving_rewire(W, signs, n_swaps=None, rng=None):
    """Rewire same-sign edges, keeping row weights and in/out degree counts.

    On the sparse blocks used here, the fitted residual is often unchanged
    or no valid swap occurs. Equality cannot establish degree causality.
    """
    if rng is None:
        rng = np.random.RandomState()
    W_new = W.copy()
    n = W.shape[0]
    edges = [(i, j) for i in range(n) for j in range(n)
             if i != j and W[i, j] != 0]
    if len(edges) < 2:
        return W_new
    if n_swaps is None:
        n_swaps = len(edges) * 10

    for _ in range(n_swaps):
        idx1, idx2 = rng.choice(len(edges), 2, replace=False)
        i, j = edges[idx1]
        k, l = edges[idx2]
        if i == k or j == l or i == l or k == j:
            continue
        if W_new[i, l] != 0 or W_new[k, j] != 0:
            continue
        if signs[j] != signs[l]:
            continue
        W_new[i, l] = W_new[i, j]
        W_new[k, j] = W_new[k, l]
        W_new[i, j] = 0
        W_new[k, l] = 0
        edges[idx1] = (i, l)
        edges[idx2] = (k, j)
    return W_new


def main():
    W_norm_full, signs_full = load_worm()
    exc_idx = np.where(signs_full > 0)[0]
    inh_idx = np.where(signs_full < 0)[0]

    tau_E, tau_I = 10.0, 20.0
    n_half = 5
    n_total = 2 * n_half
    n_trials = 200

    print(f'Thresholded scalar-center / spectral readings on C. elegans subcircuits ({n_half}E + {n_half}I)')
    print(f'tau_E={tau_E}, tau_I={tau_I}, alpha=0.3, n_trials={n_trials}')
    print(f'Threshold labels only: truly scalar_r<0.01, soft assigned_complex_error/max|mu|<0.05')
    print('=' * 78)

    sources = ('worm', 'random_dale', 'degree_preserved')
    counts = {s: {'truly': 0, 'soft': 0, 'hard': 0} for s in sources}
    residuals = {s: [] for s in sources}

    for trial in range(n_trials):
        rng = np.random.RandomState(trial + 1000)

        e_pick = rng.choice(exc_idx, n_half, replace=False)
        i_pick = rng.choice(inh_idx, n_half, replace=False)
        idx = np.concatenate([e_pick, i_pick])
        W_sub = W_norm_full[np.ix_(idx, idx)]
        signs_sub = signs_full[idx]

        # Worm subcircuit
        klass, r_tot, _ = classify_subcircuit(
            W_sub, signs_sub, tau_E, tau_I, alpha=0.3)
        counts['worm'][klass] += 1
        residuals['worm'].append(r_tot)

        # Random Dale's-law control of matched density (Erdős-Rényi-Dale)
        density = np.count_nonzero(W_sub) / (n_total * (n_total - 1))
        W_rand = np.zeros((n_total, n_total))
        for i in range(n_total):
            for j in range(n_total):
                if i != j and rng.random() < max(density, 0.01):
                    W_rand[i, j] = signs_sub[j] * rng.exponential(0.3)
        mx = np.max(np.abs(W_rand))
        if mx > 0:
            W_rand /= mx
        klass_r, r_tot_r, _ = classify_subcircuit(
            W_rand, signs_sub, tau_E, tau_I, alpha=0.3)
        counts['random_dale'][klass_r] += 1
        residuals['random_dale'].append(r_tot_r)

        # Degree-preserving null: inspect its score without inferring causality.
        W_dp = degree_preserving_rewire(W_sub, signs_sub, rng=rng)
        klass_dp, r_tot_dp, _ = classify_subcircuit(
            W_dp, signs_sub, tau_E, tau_I, alpha=0.3)
        counts['degree_preserved'][klass_dp] += 1
        residuals['degree_preserved'].append(r_tot_dp)

    print(f'{"":<18s} | {"truly":>8s} | {"soft":>8s} | {"hard":>8s} | scalar residual: median')
    print('-' * 78)
    for src in sources:
        c = counts[src]
        n = sum(c.values())
        med_r = np.median(residuals[src])
        print(f'{src:<18s} | {c["truly"]:>5d}/{n} | {c["soft"]:>5d}/{n} | {c["hard"]:>5d}/{n} | {med_r:.4f}')

    print()
    p_truly = {s: counts[s]['truly'] / n_trials for s in sources}
    print(f'truly fractions: worm={p_truly["worm"]:.1%}, '
          f'erdos_dale={p_truly["random_dale"]:.1%}, '
          f'degree_preserved={p_truly["degree_preserved"]:.1%}')
    print()
    print("Equal worm/degree-preserved scores do not establish degree-distribution causality.")
    print("The withdrawn fitted metric is insensitive to these rewires; sparse nulls often do not move.")
    print("Random controls retain a different weight normalization; these fractions are not biological evidence.")


if __name__ == '__main__':
    main()
