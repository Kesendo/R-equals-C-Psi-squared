#!/usr/bin/env python3
"""F77-style trichotomy classification of C. elegans subcircuits.

The existing palindrome test (algebraic_palindrome.py, celegans_palindrome.py)
gives a binary reading: how palindromic is the subcircuit, with a worm/
random ratio of ~8x. This script refines that by classifying each
subcircuit into one of three classes — the neural analogue of the F77
trichotomy on the quantum side:

  truly:  Q*J*Q + J + 2S*I = 0 holds exactly (algebraic equation closes,
          residual ~ 0). Strict palindrome, no residual operator.

  soft:   eigenvalues of J still come in palindromic pairs (lambda <-> -2S - lambda),
          but the algebraic equation has nonzero residual. Spectral
          symmetry without operator symmetry.

  hard:   eigenvalues do NOT come in palindromic pairs. Both spectral
          and operator symmetry fail.

Thresholds:
  truly_residual < 0.01 * ||J||  -> algebraic equation holds to 1%
  paired_max_dev < 0.05 * |lambda|  -> every eigenvalue has a partner within 5%

Tested over 200 random C. elegans 5E+5I subcircuits and 200 random
Dale's-law controls of matched density.

Self-contained (no imports from sibling neural scripts to avoid
top-level-script side effects).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

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


def palindrome_residual_norm(J, Q):
    """Return ||Q*J*Q + J + 2S*I||_F / ||J||_F where S is chosen to zero
    the diagonal (matches the algebraic_palindrome.py convention)."""
    QJQ = Q @ J @ Q.T
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2.0 * np.diag(S_diag)
    norm_J = np.linalg.norm(J)
    return float(np.linalg.norm(R) / norm_J), float(S_diag.mean())


def palindromic_pairing_max_deviation(J, S):
    """For each eigenvalue lambda of J, find the closest partner candidate
    among {-2S - lambda_j} and return the worst (max) deviation as a fraction
    of |lambda|.

    Small => paired. Large => some eigenvalues are unpartnered.
    """
    evals = np.linalg.eigvals(J)
    n = len(evals)
    max_rel_dev = 0.0
    for i in range(n):
        target = -2.0 * S - evals[i]
        deviations = np.abs(evals - target)
        nearest = float(np.min(deviations))
        scale = max(abs(evals[i]), 1e-10)
        rel_dev = nearest / scale
        if rel_dev > max_rel_dev:
            max_rel_dev = rel_dev
    return max_rel_dev


def classify_subcircuit(W, signs, tau_E, tau_I, alpha=0.3,
                        truly_threshold=0.01, pairing_threshold=0.05):
    """Return ('truly'|'soft'|'hard', residual_norm, max_pair_deviation)."""
    _, Q = build_swap(signs)
    J = build_jacobian(W, tau_E, tau_I, signs, alpha=alpha)
    r_tot, S = palindrome_residual_norm(J, Q)
    pair_dev = palindromic_pairing_max_deviation(J, S)

    if r_tot < truly_threshold:
        return 'truly', r_tot, pair_dev
    if pair_dev < pairing_threshold:
        return 'soft', r_tot, pair_dev
    return 'hard', r_tot, pair_dev


def main():
    W_norm_full, signs_full = load_worm()
    exc_idx = np.where(signs_full > 0)[0]
    inh_idx = np.where(signs_full < 0)[0]

    tau_E, tau_I = 10.0, 20.0
    n_half = 5
    n_total = 2 * n_half
    n_trials = 200

    print(f'F77 trichotomy on C. elegans subcircuits ({n_half}E + {n_half}I)')
    print(f'tau_E={tau_E}, tau_I={tau_I}, alpha=0.3, n_trials={n_trials}')
    print(f'Thresholds: truly r_tot<0.01, soft pair_dev<0.05')
    print('=' * 78)

    counts = {'worm': {'truly': 0, 'soft': 0, 'hard': 0},
              'random_dale': {'truly': 0, 'soft': 0, 'hard': 0}}
    residuals = {'worm': [], 'random_dale': []}
    pair_devs = {'worm': [], 'random_dale': []}

    for trial in range(n_trials):
        rng = np.random.RandomState(trial + 1000)

        e_pick = rng.choice(exc_idx, n_half, replace=False)
        i_pick = rng.choice(inh_idx, n_half, replace=False)
        idx = np.concatenate([e_pick, i_pick])
        W_sub = W_norm_full[np.ix_(idx, idx)]
        signs_sub = signs_full[idx]

        # Worm subcircuit
        klass, r_tot, pair_dev = classify_subcircuit(
            W_sub, signs_sub, tau_E, tau_I, alpha=0.3)
        counts['worm'][klass] += 1
        residuals['worm'].append(r_tot)
        pair_devs['worm'].append(pair_dev)

        # Random Dale's-law control of matched density
        density = np.count_nonzero(W_sub) / (n_total * (n_total - 1))
        W_rand = np.zeros((n_total, n_total))
        for i in range(n_total):
            for j in range(n_total):
                if i != j and rng.random() < max(density, 0.01):
                    W_rand[i, j] = signs_sub[j] * rng.exponential(0.3)
        mx = np.max(np.abs(W_rand))
        if mx > 0:
            W_rand /= mx
        klass_r, r_tot_r, pair_dev_r = classify_subcircuit(
            W_rand, signs_sub, tau_E, tau_I, alpha=0.3)
        counts['random_dale'][klass_r] += 1
        residuals['random_dale'].append(r_tot_r)
        pair_devs['random_dale'].append(pair_dev_r)

    print(f'{"":<14s} | {"truly":>8s} | {"soft":>8s} | {"hard":>8s} | residual: median  | pair_dev: median')
    print('-' * 78)
    for src in ('worm', 'random_dale'):
        c = counts[src]
        n = sum(c.values())
        med_r = np.median(residuals[src])
        med_p = np.median(pair_devs[src])
        print(f'{src:<14s} | {c["truly"]:>5d}/{n} | {c["soft"]:>5d}/{n} | {c["hard"]:>5d}/{n} | {med_r:.4f}           | {med_p:.4f}')

    print()
    p_soft_worm = counts['worm']['soft'] / n_trials
    p_soft_rand = counts['random_dale']['soft'] / n_trials
    p_hard_worm = counts['worm']['hard'] / n_trials
    p_hard_rand = counts['random_dale']['hard'] / n_trials
    print(f'soft fraction: worm={p_soft_worm:.2%}, random={p_soft_rand:.2%}, '
          f'enrichment={p_soft_worm/max(p_soft_rand, 1e-3):.2f}x')
    print(f'hard fraction: worm={p_hard_worm:.2%}, random={p_hard_rand:.2%}')


if __name__ == '__main__':
    main()
