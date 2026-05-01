#!/usr/bin/env python3
"""F77-style trichotomy classification of C. elegans subcircuits.

The existing palindrome test (algebraic_palindrome.py, celegans_palindrome.py)
gives a binary reading: how palindromic is the subcircuit, with a worm/
random ratio of ~8x. This script refines that by classifying each
subcircuit into one of three classes, the neural analogue of the F77
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


def degree_preserving_rewire(W, signs, n_swaps=None, rng=None):
    """Degree-preserving randomization: keep each neuron's in/out-degree
    sequence, randomly reroute edges. Inlined from validation_checks.py
    to avoid top-level-script side effects on import.

    This is the proper null model: it controls for degree distribution
    so any palindrome advantage seen against this null is purely a
    wiring effect, not a degree-distribution artifact.
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

    print(f'F77 trichotomy on C. elegans subcircuits ({n_half}E + {n_half}I)')
    print(f'tau_E={tau_E}, tau_I={tau_I}, alpha=0.3, n_trials={n_trials}')
    print(f'Thresholds: truly r_tot<0.01, soft pair_dev<0.05')
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

        # Degree-preserving null: keeps the worm's degree sequence,
        # randomly reroutes edges. Tests whether the trichotomy
        # enrichment is due to specific wiring or just degree distribution
        # (per the existing caveat in docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md
        # that the binary palindrome advantage is fully explained by
        # degree distribution).
        W_dp = degree_preserving_rewire(W_sub, signs_sub, rng=rng)
        klass_dp, r_tot_dp, _ = classify_subcircuit(
            W_dp, signs_sub, tau_E, tau_I, alpha=0.3)
        counts['degree_preserved'][klass_dp] += 1
        residuals['degree_preserved'].append(r_tot_dp)

    print(f'{"":<18s} | {"truly":>8s} | {"soft":>8s} | {"hard":>8s} | residual: median')
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
    print('Reading: if "worm" and "degree_preserved" are similar but both '
          'much higher than "random_dale", the trichotomy enrichment is due '
          'to degree distribution alone (per existing caveat). If "worm" is '
          'distinctly higher than "degree_preserved", specific wiring contributes.')


if __name__ == '__main__':
    main()
