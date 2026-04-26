#!/usr/bin/env python3
"""Neural networks viewed through the framework's V-Effect lens.

Earlier draft of this script tested "bit_a-even (E-E + I-I) only vs full"
and was naive: removing the E-I cross-couplings disconnects the two
populations, and two disconnected damped networks are trivially
palindromic. That comparison says nothing about the framework.

The framework's actual claim is the V-Effect's: TWO palindromic
populations coupled through a BRIDGE produce a palindromic joint
network if and only if the bridge satisfies a specific algebraic
condition. The condition (March 2026) is:

    |W[Q(i), Q(j)]|  =  (τ_{Q(i)}/τ_i)  ·  |W[i,j]|

When this magnitude relation holds across the E-I bridge, the joint
Wilson-Cowan Jacobian satisfies Q · J · Q⁻¹ + J + 2·S = 0 exactly.
When it fails, the bridge is the source of palindrome-breaking.

This script tests the bridge condition on C. elegans cross-couplings.

Run: `python simulations/_neural_framework_lens.py`
"""
import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))


def palindrome_residual_generic(M, A, S):
    if np.allclose(A @ A.conj().T, np.eye(A.shape[0]), atol=1e-9):
        A_inv = A.conj().T
    else:
        A_inv = np.linalg.inv(A)
    n = M.shape[0]
    if np.isscalar(S):
        S_mat = S * np.eye(n)
    else:
        S_mat = np.diag(S) if S.ndim == 1 else S
    return A @ M @ A_inv + M + 2 * S_mat


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


def relative_residual(J, Q):
    QJQ = Q @ J @ Q.T
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = palindrome_residual_generic(J, Q, S_diag)
    norm_J = np.linalg.norm(J)
    return float(np.linalg.norm(R) / norm_J) if norm_J > 0 else 0.0


def bridge_magnitude_deviation(W, signs, perm, tau_E, tau_I):
    """For each E-I cross-coupling, compute deviation from V-Effect's
    magnitude condition |W[Q(i),Q(j)]| = (τ_{Q(i)}/τ_i) · |W[i,j]|.

    Returns:
        mean_deviation: average |observed_ratio − predicted_ratio| / predicted
        n_bridge: number of E-I + I-E couplings tested
    """
    n = len(signs)
    deviations = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # Cross-coupling only (the bridge)
            if signs[i] * signs[j] >= 0:
                continue
            qi, qj = perm[i], perm[j]
            tau_i = tau_E if signs[i] > 0 else tau_I
            tau_qi = tau_E if signs[qi] > 0 else tau_I
            predicted_ratio = tau_qi / tau_i

            w_ij = W[i, j]
            w_qi_qj = W[qi, qj]
            if abs(w_ij) > 1e-10:
                actual_ratio = abs(w_qi_qj / w_ij)
                deviations.append(abs(actual_ratio - predicted_ratio) / predicted_ratio)
    if not deviations:
        return 0.0, 0
    return float(np.mean(deviations)), len(deviations)


def make_corrected_W(W, signs, perm, tau_E, tau_I):
    """Return a copy of W where each E-I + I-E coupling is forced to satisfy
    the V-Effect magnitude condition. The signs and the (i,j) magnitudes
    are kept; the (Q(i), Q(j)) magnitudes are overwritten so that the
    ratio matches the prediction. This is the "what if the bridge were
    perfect" counterfactual.
    """
    n = len(signs)
    W_corrected = W.copy()
    seen = set()
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if signs[i] * signs[j] >= 0:
                continue
            qi, qj = perm[i], perm[j]
            if (i, j) in seen or (qi, qj) in seen:
                continue
            tau_i = tau_E if signs[i] > 0 else tau_I
            tau_qi = tau_E if signs[qi] > 0 else tau_I
            ratio = tau_qi / tau_i
            # Force |W[Q(i), Q(j)]| = ratio · |W[i, j]|, keep its sign
            w_ij = W[i, j]
            sign_q = np.sign(W[qi, qj]) if W[qi, qj] != 0 else (1 if signs[qj] > 0 else -1)
            W_corrected[qi, qj] = sign_q * abs(w_ij) * ratio
            seen.add((i, j))
            seen.add((qi, qj))
    return W_corrected


# ════════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════════

def main():
    print()
    print("█" * 78)
    print("█" + "  Neural networks through the framework's V-Effect lens  ".center(76) + "█")
    print("█" * 78)
    print()
    print("Framework's V-Effect prediction: two palindromic populations")
    print("(E and I) coupled through a bridge satisfy Q·J·Q⁻¹ + J + 2·S = 0")
    print("if and only if the bridge satisfies the magnitude condition")
    print("    |W[Q(i),Q(j)]| = (τ_{Q(i)}/τ_i) · |W[i,j]|")
    print("for each cross-coupling.")
    print()
    print("Test: on C. elegans subnetworks, measure (a) the actual residual")
    print("and (b) what the residual WOULD BE if the bridge were forced to")
    print("satisfy the condition. The gap between (a) and (b) tells us how")
    print("much palindrome the worm's bridge already produces vs how much")
    print("would be possible.")
    print()

    connectome_path = SCRIPT_DIR / "neural" / "celegans_connectome.json"
    with open(connectome_path) as f:
        data = json.load(f)
    W_chem = np.array(data["chemical"])
    signs_full = np.array(data["chemical_sign"])
    N_full = len(signs_full)

    W_signed = np.zeros((N_full, N_full))
    for i in range(N_full):
        for j in range(N_full):
            W_signed[i, j] = signs_full[j] * W_chem[j, i]
    W_norm = W_signed / np.max(np.abs(W_signed))

    exc_idx = np.where(signs_full > 0)[0]
    inh_idx = np.where(signs_full < 0)[0]
    tau_E, tau_I = 10.0, 20.0

    print(f"Worm: N={N_full} neurons, {len(exc_idx)} E + {len(inh_idx)} I")
    print(f"τ_E={tau_E}, τ_I={tau_I}, predicted bridge ratio = {tau_I/tau_E:.1f}")
    print(f"α=0.3, 200 trials per size")
    print()
    print(f"  {'size':>6s}    {'‖R‖/‖J‖ actual':>16s}    {'‖R‖/‖J‖ corrected':>18s}    "
          f"{'bridge dev':>12s}    {'gap closed':>12s}")
    print(f"  {'─' * 6}    {'─' * 16}    {'─' * 18}    {'─' * 12}    {'─' * 12}")

    n_trials = 200
    for n_half in [5, 10, 13]:
        if n_half > len(inh_idx):
            continue
        n_total = 2 * n_half
        actual_residuals = []
        corrected_residuals = []
        bridge_deviations = []
        for trial in range(n_trials):
            rng = np.random.RandomState(trial + 7000)
            e_pick = rng.choice(exc_idx, n_half, replace=False)
            i_pick = rng.choice(inh_idx, n_half, replace=False)
            idx = np.concatenate([e_pick, i_pick])
            W_sub = W_norm[np.ix_(idx, idx)]
            signs_sub = signs_full[idx]
            perm, Q = build_swap(signs_sub)

            J_actual = build_jacobian(W_sub, tau_E, tau_I, signs_sub)
            actual_residuals.append(relative_residual(J_actual, Q))

            W_corrected = make_corrected_W(W_sub, signs_sub, perm, tau_E, tau_I)
            J_corrected = build_jacobian(W_corrected, tau_E, tau_I, signs_sub)
            corrected_residuals.append(relative_residual(J_corrected, Q))

            dev, _ = bridge_magnitude_deviation(W_sub, signs_sub, perm, tau_E, tau_I)
            bridge_deviations.append(dev)

        actual_mean = np.mean(actual_residuals)
        corrected_mean = np.mean(corrected_residuals)
        bridge_dev_mean = np.mean(bridge_deviations)
        # "Gap closed" = how much of the actual residual is removed by correcting the bridge
        gap_closed = (actual_mean - corrected_mean) / actual_mean if actual_mean > 0 else 0
        print(f"  {n_total:>6d}    {actual_mean:>16.4e}    {corrected_mean:>18.4e}    "
              f"{bridge_dev_mean:>12.4f}    {gap_closed:>12.1%}")

    print()
    print("Reading the table:")
    print("  ‖R‖ actual       = residual on the worm's bridge as-is.")
    print("  ‖R‖ corrected    = residual after forcing the V-Effect magnitude")
    print("                     condition on each cross-coupling.")
    print("  bridge dev       = mean relative deviation of |W[Q(i),Q(j)]|/|W[i,j]|")
    print("                     from the predicted τ_I/τ_E = 2.0.")
    print("  gap closed       = fraction of residual removed by correcting bridge.")
    print()
    print("If gap closed ≈ 100%: the bridge is essentially the only source")
    print("of palindrome-breaking, and the V-Effect framework explains it all.")
    print()
    print("If gap closed << 100%: the worm's palindrome breaking has sources")
    print("beyond the bridge — within-population magnitude or sign mismatches.")


if __name__ == "__main__":
    main()
