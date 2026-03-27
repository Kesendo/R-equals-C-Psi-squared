#!/usr/bin/env python3
"""
Is the Wilson-Cowan palindrome EXACT or just approximate?

In the quantum Lindblad case:
  lambda_k + lambda_{paired} = -2*sum_gamma  (EXACT, for all k)

This holds for ANY coupling strength. It's an algebraic identity, not
a numerical coincidence.

Question: does the Wilson-Cowan Jacobian have the same property?
  mu_k + mu_{paired} = constant?

If yes: there must be a hidden algebraic structure (a classical Pi).
If no: the palindrome is approximate and the 98.2% is a sparsity artifact.
"""
import numpy as np
import json
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def build_jacobian(W, tau_exc, tau_inh, signs, alpha=0.3):
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_exc if signs[i] > 0 else tau_inh
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def best_pairing(eigenvalues):
    """Find the best pairing of eigenvalues and report pair sums.

    For N eigenvalues (N even), find the pairing that minimizes
    the variance of pair sums. Returns pair sums and the pairing.
    """
    evs = np.array(sorted(eigenvalues, key=lambda x: x.real))
    n = len(evs)

    # Strategy: pair smallest-real with largest-real (greedy by real part)
    pairs = []
    used = set()
    order = np.argsort([e.real for e in evs])

    left = 0
    right = n - 1
    while left < right:
        while left in used:
            left += 1
        while right in used:
            right -= 1
        if left < right:
            pairs.append((order[left], order[right]))
            used.add(left)
            used.add(right)
            left += 1
            right -= 1

    pair_sums = [evs[i].real + evs[j].real for i, j in pairs]
    return pair_sums, pairs, evs


# === Test 1: Pure structure (no connectome, synthetic network) ===
print("=" * 65)
print("TEST 1: Synthetic balanced network, varying coupling strength")
print("=" * 65)

N = 10
n_exc = 5
np.random.seed(42)
signs = np.ones(N)
signs[n_exc:] = -1
np.random.shuffle(signs)

# Random weight matrix (moderate density)
W = np.random.randn(N, N) * 0.3
np.fill_diagonal(W, 0)
for i in range(N):
    for j in range(N):
        W[i, j] = signs[j] * abs(W[i, j])
W /= np.max(np.abs(W))

tau_E, tau_I = 10.0, 20.0
predicted_sum = -(1.0 / tau_E + 1.0 / tau_I)

print(f"\nN={N}, E={n_exc}, I={N-n_exc}")
print(f"tau_E={tau_E}, tau_I={tau_I}")
print(f"Predicted pairing sum (from tau): {predicted_sum:.4f}")
print(f"  (if palindrome is exact: all pair sums = {predicted_sum:.4f})")
print()
print(f"{'alpha':>8s}  {'mean_sum':>10s}  {'std_sum':>10s}  "
      f"{'min_sum':>10s}  {'max_sum':>10s}  {'exact?':>8s}")
print("-" * 65)

for alpha in [0.0, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0]:
    J = build_jacobian(W, tau_E, tau_I, signs, alpha=alpha)
    ev = np.linalg.eigvals(J)
    sums, pairs, evs = best_pairing(ev)

    mean_s = np.mean(sums)
    std_s = np.std(sums)
    exact = "YES" if std_s < 1e-10 else ("~" if std_s < 0.001 else "no")

    print(f"  {alpha:6.3f}  {mean_s:10.6f}  {std_s:10.6f}  "
          f"{np.min(sums):10.6f}  {np.max(sums):10.6f}  {exact:>8s}")


# === Test 2: Algebraic check - does Q exist? ===
print("\n" + "=" * 65)
print("TEST 2: Does a conjugation Q exist for Wilson-Cowan?")
print("=" * 65)

print("""
For Lindblad: Pi L Pi^{-1} = -L - 2Sg*I (EXACT)
For Wilson-Cowan: Q J Q^{-1} = -J - c*I ?

If J = T*(-I + alpha*W) where T = diag(1/tau_i):
  Q must satisfy: Q*T*Q^{-1} maps tau_E <-> tau_I
  AND: Q*(T*W)*Q^{-1} = -(T*W)

The first condition is a permutation swapping E and I sites.
The second requires T'*W' = -(T*W) where T', W' are permuted.
""")

# Try swapping E and I neurons
e_idx = np.where(signs > 0)[0]
i_idx = np.where(signs < 0)[0]

# Build swap permutation (pair each E with an I)
n_pairs = min(len(e_idx), len(i_idx))
perm = np.arange(N)
for k in range(n_pairs):
    perm[e_idx[k]] = i_idx[k]
    perm[i_idx[k]] = e_idx[k]

# Build permutation matrix
Q_perm = np.zeros((N, N))
for i in range(N):
    Q_perm[i, perm[i]] = 1.0

alpha_test = 0.3
J_test = build_jacobian(W, tau_E, tau_I, signs, alpha=alpha_test)

# Check Q J Q^{-1} + J + c*I = 0
QJQ = Q_perm @ J_test @ Q_perm.T  # Q^{-1} = Q^T for permutation

# What is c?
# At alpha=0: J = diag(-1/tau_i). QJQ = diag(-1/tau_{perm(i)}).
# QJQ + J = diag(-1/tau_{perm(i)} - 1/tau_i)
# For E->I swap: -1/tau_I - 1/tau_E = -(1/tau_E + 1/tau_I) = c
c = -(1.0 / tau_E + 1.0 / tau_I)

residual = QJQ + J_test - c * np.eye(N)
err = np.linalg.norm(residual)
err_diag = np.linalg.norm(np.diag(residual))
err_offdiag = np.linalg.norm(residual - np.diag(np.diag(residual)))

print(f"Permutation swap (E <-> I):")
print(f"  ||Q*J*Q^T + J - c*I|| = {err:.6f}")
print(f"  ||diagonal residual||  = {err_diag:.6f}")
print(f"  ||off-diag residual||  = {err_offdiag:.6f}")
print(f"  c = {c:.4f}")

if err_diag < 1e-10:
    print(f"\n  Diagonal is EXACT (self-decay swaps perfectly)")
    print(f"  Off-diagonal error = {err_offdiag:.6f}")
    print(f"  This means: the palindrome is exact for self-decay")
    print(f"  but broken by coupling. The coupling error is O(alpha).")

# Check how off-diagonal error scales with alpha
print(f"\n  Off-diagonal error vs alpha:")
for alpha in [0.0, 0.01, 0.1, 0.3, 1.0, 3.0]:
    J_a = build_jacobian(W, tau_E, tau_I, signs, alpha=alpha)
    QJQ_a = Q_perm @ J_a @ Q_perm.T
    res_a = QJQ_a + J_a - c * np.eye(N)
    err_od = np.linalg.norm(res_a - np.diag(np.diag(res_a)))
    print(f"    alpha={alpha:5.2f}:  ||off-diag|| = {err_od:.6f}"
          f"  ratio to alpha: {err_od/alpha:.4f}" if alpha > 0 else
          f"    alpha={alpha:5.2f}:  ||off-diag|| = {err_od:.6f}")


# === Test 3: What WOULD make it exact? ===
print("\n" + "=" * 65)
print("TEST 3: Condition for exact palindrome in Wilson-Cowan")
print("=" * 65)

print("""
Q*J*Q^{-1} + J + c*I = 0 requires (off-diagonal):
  (1/tau_{Q(i)}) * W[Q(i),Q(j)] + (1/tau_i) * W[i,j] = 0

For E-I swap Q: tau_{Q(i)} swaps E<->I. So:
  For (i=E, j=E): (1/tau_I)*W[I_partner, I_partner] + (1/tau_E)*W[i,j] = 0
  This requires: W[I_partner, I_partner] = -(tau_I/tau_E) * W[i,j]

The coupling matrix must be ANTISYMMETRIC under the E-I swap,
scaled by the tau ratio. Random weights violate this.

But in the QUANTUM case, the analogous condition is satisfied
automatically by the Pauli algebra. The commutator [H, rho] has
a built-in antisymmetry: H*rho - rho*H.

Wilson-Cowan dynamics dX/dt = JX has no such built-in antisymmetry.
The palindrome is exact only at zero coupling and approximate otherwise.
""")

# Verify: construct a W that satisfies the exact condition
print("Constructing a W that satisfies the exact palindrome condition...")
W_exact = np.zeros((N, N))
ratio = tau_I / tau_E

for i in range(N):
    for j in range(N):
        if i == j:
            continue
        qi = perm[i]
        qj = perm[j]
        # Set W[qi,qj] = -(tau_{Q(i)}/tau_i) * W[i,j]
        # Only set upper triangle of original, compute lower from constraint
        if i < j:
            W_exact[i, j] = np.random.randn() * 0.3
            tau_i = tau_E if signs[i] > 0 else tau_I
            tau_qi = tau_E if signs[qi] > 0 else tau_I
            W_exact[qi, qj] = -(tau_qi / tau_i) * W_exact[i, j]

W_exact /= np.max(np.abs(W_exact)) if np.max(np.abs(W_exact)) > 0 else 1

J_exact = build_jacobian(W_exact, tau_E, tau_I, signs, alpha=1.0)
ev_exact = np.linalg.eigvals(J_exact)
sums_ex, _, _ = best_pairing(ev_exact)

QJQ_ex = Q_perm @ J_exact @ Q_perm.T
res_ex = QJQ_ex + J_exact - c * np.eye(N)
err_ex = np.linalg.norm(res_ex)

print(f"  Constructed W satisfying antisymmetry condition")
print(f"  ||Q*J*Q^T + J - c*I|| = {err_ex:.2e}")
print(f"  Pair sums: mean={np.mean(sums_ex):.6f}, std={np.std(sums_ex):.6f}")
print(f"  Predicted: {c:.6f}")
if np.std(sums_ex) < 1e-8:
    print(f"  >>> EXACT palindrome achieved with constrained W!")
