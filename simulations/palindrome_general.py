#!/usr/bin/env python3
"""The general palindrome condition, formulated and verified across substrates.

Claim. A linear open-system generator L has a spectrum palindromic about -c iff there is an
involution Pi (Pi^2 = I) with

    Pi L Pi = -L - 2c I        (equivalently: the Pi-even part of L is the scalar -c I).

Split L = A + B (A = coupling/conservative, B = dissipative diagonal). Two sufficient sub-conditions:
    (i)  Pi A Pi = -A          (the coupling is mirror-ODD: contributes 0 to the centre -> topology-blind)
    (ii) B + Pi B Pi = -2c I   (the bath rates pair to a constant sum under Pi)
Then Pi L Pi = (A + Pi A Pi) + (B + Pi B Pi) = 0 + (-2c I) = -2c I, and for L v = lam v one has
L(Pi v) = (-2c - lam)(Pi v): the partner mode, palindromic about -c. The centre -c is set by B (the
bath) alone, never by A (the coupling): the topology-blindness, substrate-independent.

We verify on three substrates: (1) a fresh abstract system built from the two sub-conditions (the
generalization itself, neither quantum nor neural); (2) the Wilson-Cowan E-I Jacobian (exact
palindromic network, Q the E-I swap); (3) the quantum Z-dephasing Liouvillian (the original F1).
"""
import numpy as np


def spectrum_palindromic(L, c, tol=1e-7):
    """True if every eigenvalue lam has a partner -2c - lam in the spectrum."""
    ev = np.linalg.eigvals(L)
    partner = -2 * c - ev
    used = np.zeros(len(ev), bool)
    worst = 0.0
    for p in partner:
        d = np.abs(ev - p) + np.where(used, np.inf, 0.0)
        k = int(np.argmin(d))
        worst = max(worst, float(d[k]))
        used[k] = True
    return worst < tol, worst


# ---- 1. the generalization: a fresh abstract system from the two sub-conditions ----

def abstract_demo(n=12, seed=0):
    rng = np.random.RandomState(seed)
    # Pi = a random involution permutation (pairs of indices, allow fixed points)
    order = rng.permutation(n)
    perm = np.arange(n)
    k = 0
    while k + 1 < n:
        if rng.rand() < 0.85:                      # pair them
            perm[order[k]] = order[k + 1]
            perm[order[k + 1]] = order[k]
            k += 2
        else:                                      # leave order[k] a fixed point
            k += 1
    Pi = np.zeros((n, n))
    for i in range(n):
        Pi[i, perm[i]] = 1.0

    c = 1.0
    # A: any matrix, projected to be Pi-ODD:  A = (M - Pi M Pi)/2  =>  Pi A Pi = -A
    M = rng.randn(n, n) + 1j * rng.randn(n, n)
    A = (M - Pi @ M @ Pi) / 2.0
    # B: diagonal, paired to sum -2c under Pi; fixed points get -c
    b = np.zeros(n, complex)
    done = np.zeros(n, bool)
    for i in range(n):
        if done[i]:
            continue
        j = perm[i]
        if i == j:
            b[i] = -c
            done[i] = True
        else:
            bi = -(0.2 + rng.rand())               # a generic negative rate
            b[i], b[j] = bi, -2 * c - bi
            done[i] = done[j] = True
    L = A + np.diag(b)

    cond = np.linalg.norm(L + Pi @ L @ Pi + 2 * c * np.eye(n))
    pal, worst = spectrum_palindromic(L, c)
    return cond, pal, worst, c


# ---- 2. the neural Wilson-Cowan Jacobian (builders from neural_clock_two_hands.py) ----

def build_exact_palindromic_network(N, tau_E, tau_I, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    signs[rng.choice(N, N // 2, replace=False)] = -1
    e_idx = list(np.where(signs > 0)[0]); i_idx = list(np.where(signs < 0)[0])
    perm = np.arange(N)
    for k in range(min(len(e_idx), len(i_idx))):
        perm[e_idx[k]] = i_idx[k]; perm[i_idx[k]] = e_idx[k]
    W = np.zeros((N, N))
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    for i in range(N):
        for j in range(N):
            if i == j or i >= j or not mask[i, j]:
                continue
            qi, qj = perm[i], perm[j]
            base = rng.exponential(0.3)
            W[i, j] = signs[j] * base
            tau_i = tau_E if signs[i] > 0 else tau_I
            tau_qi = tau_E if signs[qi] > 0 else tau_I
            W[qi, qj] = -(tau_qi / tau_i) * W[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def build_linear_jacobian(W, signs, tau_E, tau_I, alpha):
    n = len(signs); J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_E if signs[i] > 0 else tau_I
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def swap_Q(signs):
    n = len(signs)
    e_idx = list(np.where(signs > 0)[0]); i_idx = list(np.where(signs < 0)[0])
    perm = np.arange(n)
    for k in range(min(len(e_idx), len(i_idx))):
        perm[e_idx[k]] = i_idx[k]; perm[i_idx[k]] = e_idx[k]
    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, perm[i]] = 1.0
    return Q


def neural_demo(N=20, tau_E=5.0, tau_I=10.0, alpha=0.5):
    W, signs = build_exact_palindromic_network(N, tau_E, tau_I)
    J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)
    Q = swap_Q(signs)
    S = (1.0 / tau_E + 1.0 / tau_I) / 2.0
    cond = np.linalg.norm(J + Q @ J @ Q + 2 * S * np.eye(N))
    pal, worst = spectrum_palindromic(J, S)
    return cond, pal, worst, S, np.linalg.norm(J)


# ---- 3. the quantum Z-dephasing Liouvillian (the original F1) ----

def quantum_demo(N=2, gamma=0.7, J=1.0):
    d = 2 ** N
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]], complex)
    Z = np.array([[1, 0], [0, -1]], complex)

    def op(P, k):
        o = np.array([[1]], complex)
        for i in range(N):
            o = np.kron(o, P if i == k else I2)
        return o

    H = np.zeros((d, d), complex)
    for b in range(N - 1):
        H += J * (op(X, b) @ op(X, b + 1) + op(Y, b) @ op(Y, b + 1))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = op(Z, k)
        L += gamma * (np.kron(Zk, Zk) - np.kron(Id, Id))
    c = N * gamma                                  # palindrome centre = Sum gamma = N*gamma (uniform)
    pal, worst = spectrum_palindromic(L, c)
    return pal, worst, c


def main():
    print("THE GENERAL PALINDROME CONDITION:  exists involution Pi with Pi L Pi = -L - 2c I\n")

    print("1. FRESH abstract system (Pi-odd coupling A + Pi-paired bath B; not quantum, not neural):")
    for seed in range(4):
        cond, pal, worst, c = abstract_demo(n=12, seed=seed)
        print(f"   seed {seed}:  ||L + Pi L Pi + 2cI|| = {cond:.2e}   "
              f"spectrum palindromic about -c={-c}:  {pal} (worst {worst:.2e})")
    print("   -> the two sub-conditions ALONE force the palindrome, on a fresh abstract generator.\n")

    print("2. Wilson-Cowan E-I Jacobian (exact palindromic network, Q = the E-I swap):")
    for alpha in [0.2, 0.5, 1.0]:
        cond, pal, worst, S, nJ = neural_demo(alpha=alpha)
        print(f"   alpha={alpha}:  ||J + Q J Q + 2S I||/||J|| = {cond/nJ:.2e}   "
              f"palindromic about -S={-S:.3f}: {pal} (worst {worst:.2e})")
    print("   -> the classical neural Jacobian is an instance: Dale's Law makes W mirror-odd.\n")

    print("3. Quantum Z-dephasing Liouvillian (the original F1, Pi = the proven conjugation):")
    for N in [2, 3]:
        pal, worst, c = quantum_demo(N=N)
        print(f"   N={N}:  spectrum palindromic about -Sum_gamma={-c:.2f}:  {pal} (worst {worst:.2e})")
    print("   -> the original quantum palindrome is the d=2 instance of the same condition.")


if __name__ == "__main__":
    main()
