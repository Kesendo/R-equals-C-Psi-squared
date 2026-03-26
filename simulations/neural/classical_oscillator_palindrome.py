#!/usr/bin/env python3
"""
Phase 1, Test 1: Do classical damped oscillator chains have palindromic spectra?

Three coupled masses with springs (coupling k) and friction (damping γ).
State vector: [x1, v1, x2, v2, x3, v3] (positions and velocities).

The dynamics matrix A has the same form as a Lindblad system:
    coupling (springs) = Hamiltonian
    dissipation (friction) = dephasing

Question: does the eigenvalue spectrum of A exhibit palindromic pairing?
If yes: the palindrome is not quantum-specific but a property of any
coupled-dissipative system. The pattern propagates upward.

Usage: python classical_oscillator_palindrome.py
"""
import numpy as np
from itertools import combinations

def build_oscillator_chain(N, k=1.0, gamma=0.05, m=1.0):
    """
    Build dynamics matrix for N coupled damped oscillators.

    dx_i/dt = v_i
    dv_i/dt = -(gamma/m)*v_i + (k/m)*(x_{i-1} - 2*x_i + x_{i+1})

    State vector: [x1, v1, x2, v2, ..., xN, vN]
    Returns: 2N x 2N dynamics matrix A where d(state)/dt = A @ state
    """
    dim = 2 * N
    A = np.zeros((dim, dim))

    for i in range(N):
        xi = 2 * i      # position index
        vi = 2 * i + 1  # velocity index

        # dx_i/dt = v_i
        A[xi, vi] = 1.0

        # dv_i/dt = -(gamma/m)*v_i + (k/m)*(spring forces)
        A[vi, vi] = -gamma / m  # damping

        # Spring coupling to self (from both neighbors)
        n_neighbors = 0
        if i > 0:
            n_neighbors += 1
            A[vi, 2*(i-1)] = k / m  # pull from left neighbor
        if i < N - 1:
            n_neighbors += 1
            A[vi, 2*(i+1)] = k / m  # pull from right neighbor
        A[vi, xi] = -n_neighbors * k / m  # restoring force

    return A


def check_palindrome(eigenvalues, tol=0.01):
    """
    Check if decay rates (negative real parts) are palindromically paired.
    For each rate d, look for a partner at 2*center - d.
    """
    rates = sorted(-ev.real for ev in eigenvalues if abs(ev.real) > 1e-10)

    if len(rates) == 0:
        return 0, 0, 0, []

    center = (min(rates) + max(rates)) / 2.0

    paired = 0
    pairs = []
    used = [False] * len(rates)

    for i in range(len(rates)):
        if used[i]:
            continue
        target = 2 * center - rates[i]
        for j in range(len(rates)):
            if i != j and not used[j] and abs(rates[j] - target) < tol * center:
                pairs.append((rates[i], rates[j]))
                used[i] = True
                used[j] = True
                paired += 2
                break

    return paired, len(rates), center, pairs


def run_test(N, k, gamma, label=""):
    """Run palindrome test for one configuration."""
    A = build_oscillator_chain(N, k=k, gamma=gamma)
    eigenvalues = np.linalg.eigvals(A)

    paired, total, center, pairs = check_palindrome(eigenvalues)
    score = paired / total * 100 if total > 0 else 0

    print(f"\n{'='*60}")
    print(f"{label}N={N}, k={k}, gamma={gamma}")
    print(f"{'='*60}")
    print(f"Matrix size: {2*N}x{2*N}")
    print(f"Eigenvalues ({len(eigenvalues)}):")

    for ev in sorted(eigenvalues, key=lambda x: -x.real):
        rate = -ev.real
        freq = ev.imag
        if abs(ev) > 1e-10:
            print(f"  rate={rate:8.4f}  freq={freq:8.4f}")

    print(f"\nDecay rates: {total}")
    print(f"Center: {center:.4f}")
    print(f"Paired: {paired}/{total} = {score:.1f}%")

    if pairs:
        print(f"\nPalindromic pairs (rate_a + rate_b = 2*center = {2*center:.4f}):")
        for a, b in pairs:
            print(f"  {a:.4f} + {b:.4f} = {a+b:.4f}")

    return score, center, pairs


print("=" * 60)
print("PHASE 1: CLASSICAL OSCILLATOR PALINDROME TEST")
print("=" * 60)

# Test 1: N=3, uniform coupling and damping
run_test(3, k=1.0, gamma=0.1, label="Uniform: ")

# Test 2: N=5
run_test(5, k=1.0, gamma=0.1, label="Uniform: ")

# Test 3: N=7
run_test(7, k=1.0, gamma=0.1, label="Uniform: ")

# Test 4: Vary damping
print("\n\n" + "=" * 60)
print("DAMPING SWEEP (N=3, k=1.0)")
print("=" * 60)
for gamma in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]:
    score, center, _ = run_test(3, k=1.0, gamma=gamma, label=f"gamma={gamma}: ")

# Test 5: Non-uniform damping (sacrifice zone analogy)
print("\n\n" + "=" * 60)
print("NON-UNIFORM DAMPING (N=3, sacrifice zone)")
print("=" * 60)

A = np.zeros((6, 6))
gammas = [0.5, 0.01, 0.01]  # edge sacrifice
k = 1.0
for i in range(3):
    xi, vi = 2*i, 2*i+1
    A[xi, vi] = 1.0
    A[vi, vi] = -gammas[i]
    n_nb = 0
    if i > 0: A[vi, 2*(i-1)] = k; n_nb += 1
    if i < 2: A[vi, 2*(i+1)] = k; n_nb += 1
    A[vi, xi] = -n_nb * k

eigenvalues = np.linalg.eigvals(A)
print(f"\nNon-uniform gammas: {gammas}")
print(f"Eigenvalues:")
for ev in sorted(eigenvalues, key=lambda x: -x.real):
    print(f"  rate={-ev.real:8.4f}  freq={ev.imag:8.4f}")

paired, total, center, pairs = check_palindrome(eigenvalues)
print(f"\nPaired: {paired}/{total} = {paired/total*100:.1f}%")
print(f"Center: {center:.4f}, Sum gammas/2: {sum(gammas)/2:.4f}")

# Test 6: Random networks for comparison
print("\n\n" + "=" * 60)
print("RANDOM NETWORKS (N=5, 10 trials)")
print("=" * 60)
np.random.seed(42)
scores = []
for trial in range(10):
    W = np.random.randn(5, 5) * 0.5
    W = (W + W.T) / 2  # symmetric coupling
    gamma_rand = np.abs(np.random.randn(5)) * 0.1

    dim = 10
    A = np.zeros((dim, dim))
    for i in range(5):
        xi, vi = 2*i, 2*i+1
        A[xi, vi] = 1.0
        A[vi, vi] = -gamma_rand[i]
        for j in range(5):
            if i != j:
                A[vi, 2*j] = W[i, j]
        A[vi, xi] -= sum(abs(W[i, j]) for j in range(5) if i != j)

    eigenvalues = np.linalg.eigvals(A)
    paired, total, center, _ = check_palindrome(eigenvalues, tol=0.05)
    score = paired / total * 100 if total > 0 else 0
    scores.append(score)
    print(f"  Trial {trial+1}: {score:.0f}% paired")

print(f"\nRandom average: {np.mean(scores):.1f}% (vs chain: should be higher if palindrome is structural)")
