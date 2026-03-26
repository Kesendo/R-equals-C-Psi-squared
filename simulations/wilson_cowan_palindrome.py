#!/usr/bin/env python3
"""
Phase 1, Test 2: Wilson-Cowan neural population model - palindrome test.

Wilson-Cowan has SELECTIVE damping: excitatory (tau_E) and inhibitory (tau_I)
populations decay at different rates. This is the biological analogue of
the 2:2 Pauli split that makes the quantum palindrome possible.

Test: does the linearized Jacobian of a Wilson-Cowan network have
palindromic eigenvalue pairing?

Usage: python wilson_cowan_palindrome.py
"""
import numpy as np

def sigmoid(x, a=1.3, theta=4.0):
    arg = -a * (x - theta)
    arg = np.clip(arg, -500, 500)
    return 1.0 / (1.0 + np.exp(arg))

def dsigmoid(x, a=1.3, theta=4.0):
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)


def build_jacobian_chain(N, tau_E=8.0, tau_I=18.0,
                          c_EE=16.0, c_EI=12.0, c_IE=15.0, c_II=3.0,
                          c_lat=2.0, P=1.5,
                          a_E=1.3, a_I=2.0, theta_E=4.0, theta_I=3.7):
    """
    Build Jacobian for chain of N Wilson-Cowan E-I nodes.
    Uses approximate fixed point from literature (E*~0.2, I*~0.4).
    """
    # Find single-node fixed point by iteration
    E, I = 0.1, 0.1
    for _ in range(5000):
        inp_E = c_EE * E - c_EI * I + P
        inp_I = c_IE * E - c_II * I
        E_new = sigmoid(inp_E, a_E, theta_E)
        I_new = sigmoid(inp_I, a_I, theta_I)
        E = 0.99 * E + 0.01 * E_new
        I = 0.99 * I + 0.01 * I_new
    E_star, I_star = E, I

    # Sigmoid derivatives at fixed point
    inp_E = c_EE * E_star - c_EI * I_star + P
    inp_I = c_IE * E_star - c_II * I_star
    dSE = dsigmoid(inp_E, a_E, theta_E)
    dSI = dsigmoid(inp_I, a_I, theta_I)

    dim = 2 * N
    J = np.zeros((dim, dim))

    for i in range(N):
        ei = 2 * i
        ii = 2 * i + 1

        # Local: dE/dt = (-E + S_E(input)) / tau_E
        J[ei, ei] = (-1.0 + dSE * c_EE) / tau_E
        J[ei, ii] = (-dSE * c_EI) / tau_E
        J[ii, ei] = (dSI * c_IE) / tau_I
        J[ii, ii] = (-1.0 - dSI * c_II) / tau_I

        # Lateral E-E coupling
        if i > 0:
            J[ei, 2*(i-1)] += dSE * c_lat / tau_E
        if i < N - 1:
            J[ei, 2*(i+1)] += dSE * c_lat / tau_E

    return (E_star, I_star), J


def check_palindrome(eigenvalues, tol_frac=0.05):
    rates = sorted(set(round(-ev.real, 6) for ev in eigenvalues if abs(ev.real) > 1e-10))
    if len(rates) < 2:
        return 0, len(rates), 0, []

    center = (min(rates) + max(rates)) / 2.0
    paired = 0
    pairs = []
    used = [False] * len(rates)
    tol = max(tol_frac * abs(center), 0.001)

    for i in range(len(rates)):
        if used[i]: continue
        target = 2 * center - rates[i]
        for j in range(len(rates)):
            if i != j and not used[j] and abs(rates[j] - target) < tol:
                pairs.append((rates[i], rates[j]))
                used[i] = True
                used[j] = True
                paired += 2
                break

    return paired, len(rates), center, pairs


def run_test(N, tau_E, tau_I, label="", **kwargs):
    (E_star, I_star), J = build_jacobian_chain(N, tau_E=tau_E, tau_I=tau_I, **kwargs)
    eigenvalues = np.linalg.eigvals(J)
    paired, total, center, pairs = check_palindrome(eigenvalues)
    score = paired / total * 100 if total > 0 else 0

    print(f"\n{'='*60}")
    print(f"{label}N={N}, tau_E={tau_E}, tau_I={tau_I}, ratio={tau_I/tau_E:.1f}")
    print(f"Fixed point: E*={E_star:.4f}, I*={I_star:.4f}")
    print(f"Eigenvalues:")

    for ev in sorted(eigenvalues, key=lambda x: x.real):
        rate = -ev.real
        freq = ev.imag
        print(f"  rate={rate:8.4f}  freq={freq:8.4f}")

    print(f"\nDistinct rates: {total}, Center: {center:.4f}")
    print(f"Paired: {paired}/{total} = {score:.1f}%")
    if pairs:
        for a, b in pairs:
            print(f"  {a:.4f} + {b:.4f} = {a+b:.4f} (2*center={2*center:.4f})")

    return score, center, eigenvalues


print("=" * 60)
print("PHASE 1, TEST 2: WILSON-COWAN PALINDROME TEST")
print("=" * 60)

# Single node
run_test(1, 8.0, 18.0, label="Single node: ")

# Chain N=3
run_test(3, 8.0, 18.0, label="Chain: ")

# Chain N=5
run_test(5, 8.0, 18.0, label="Chain: ")

# Tau sweep
print("\n\n" + "=" * 60)
print("TAU RATIO SWEEP (N=3)")
print("=" * 60)
for tau_I in [8.0, 12.0, 18.0, 30.0, 50.0]:
    run_test(3, 8.0, tau_I, label=f"")

# Uniform vs selective
print("\n\n" + "=" * 60)
print("UNIFORM (tau_E = tau_I) VS SELECTIVE")
print("=" * 60)
run_test(3, 8.0, 8.0, label="UNIFORM: ")
run_test(3, 8.0, 18.0, label="SELECTIVE: ")
