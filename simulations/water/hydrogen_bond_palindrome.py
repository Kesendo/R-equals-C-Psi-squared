#!/usr/bin/env python3
"""
Is the hydrogen bond network palindromic?

Model: chain of water molecules. Each molecule contributes:
  - One DONOR mode (O-H stretch, fast, ~3400 cm^-1)
  - One ACCEPTOR mode (O...H libration, slow, ~700 cm^-1)

The coupling through the hydrogen bond connects D and A modes.

The key question: does the coupling have the right sign structure
for the palindrome condition Q·J·Q^{-1} + J + 2S = 0?

We let the MATH decide. Build the model, compute, check.
"""
import numpy as np
from scipy.linalg import eigvals
import os


def build_water_chain_jacobian(N_molecules, tau_D, tau_A, k_DA, k_AD,
                                alpha=1.0):
    """Build Jacobian for a chain of water molecules.

    Each molecule i has two modes:
      - D_i: donor O-H stretch (index 2*i)
      - A_i: acceptor O...H libration (index 2*i + 1)

    Self-decay: D decays at rate 1/tau_D, A at rate 1/tau_A.

    Coupling within molecule (intramolecular):
      D_i and A_i are coupled through the shared O atom.

    Coupling between molecules (hydrogen bond):
      D_i (donor of molecule i) couples to A_{i+1} (acceptor of molecule i+1)
      through the hydrogen bond H_i...O_{i+1}.

    k_DA: coupling from donor mode to acceptor mode (D excites A?)
    k_AD: coupling from acceptor mode to donor mode (A excites D?)

    If k_DA and k_AD have OPPOSITE SIGNS: antisymmetric coupling.
    This is the Dale's Law analog for water.
    """
    N = 2 * N_molecules  # total modes (D + A per molecule)
    J = np.zeros((N, N))

    for i in range(N_molecules):
        d_idx = 2 * i      # donor mode index
        a_idx = 2 * i + 1  # acceptor mode index

        # Self-decay
        J[d_idx, d_idx] = -1.0 / tau_D
        J[a_idx, a_idx] = -1.0 / tau_A

        # Intramolecular coupling (D_i <-> A_i within same molecule)
        # The O atom connects the donor O-H and the acceptor O...H
        J[d_idx, a_idx] = alpha * k_DA / tau_D
        J[a_idx, d_idx] = alpha * k_AD / tau_A

        # Intermolecular coupling (hydrogen bond to next molecule)
        if i < N_molecules - 1:
            next_a = 2 * (i + 1) + 1  # acceptor of next molecule
            next_d = 2 * (i + 1)      # donor of next molecule

            # D_i -> A_{i+1}: donor stretch pushes H toward next acceptor
            J[next_a, d_idx] = alpha * k_DA / tau_A

            # A_{i+1} -> D_i: acceptor of next molecule pulls back
            J[d_idx, next_a] = alpha * k_AD / tau_D

    return J


def build_swap_DA(N_molecules):
    """Build Q: swap donor and acceptor modes at each molecule."""
    N = 2 * N_molecules
    Q = np.zeros((N, N))
    for i in range(N_molecules):
        d_idx = 2 * i
        a_idx = 2 * i + 1
        Q[d_idx, a_idx] = 1.0  # D -> A
        Q[a_idx, d_idx] = 1.0  # A -> D
    return Q


def palindrome_residual(J, Q):
    N = J.shape[0]
    QJQ = Q @ J @ Q.T  # Q^{-1} = Q^T = Q for swap
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2 * np.diag(S_diag)
    R_off = R - np.diag(np.diag(R))
    norm_J = np.linalg.norm(J)
    return np.linalg.norm(R_off) / norm_J if norm_J > 0 else 0


# ================================================================
# The fundamental question: what signs does the coupling have?
# ================================================================
print("=" * 65)
print("HYDROGEN BOND COUPLING: What are the signs?")
print("=" * 65)

print("""
In a hydrogen bond O1-H...O2:

When the donor O-H stretches (H moves toward O2):
  -> The acceptor O...H distance DECREASES (bond strengthens)
  -> The acceptor mode gets COMPRESSED (negative displacement)
  -> Effect on acceptor: NEGATIVE (opposing the stretch)

When the acceptor O...H stretches (O2 moves away from H):
  -> The donor O-H feels less pull from O2
  -> The donor mode can relax (contract back)
  -> Effect on donor: also NEGATIVE? Or POSITIVE?

The physics: stretching one side of the H-bond
COMPRESSES the other side. This is a RESTORING force.
In terms of displacement coordinates: NEGATIVE coupling.

BUT: is D->A coupling the same sign as A->D?
If both are negative: symmetric, no palindrome.
If one is positive and one negative: ANTISYMMETRIC, palindrome!

Let us compute both cases and see which gives palindromic structure.
""")

tau_D = 1.0 / 3400  # donor O-H stretch period (fast)
tau_A = 1.0 / 700   # acceptor libration period (slow)

print(f"tau_D = {tau_D:.6f} (fast, O-H stretch)")
print(f"tau_A = {tau_A:.6f} (slow, libration)")
print(f"tau_A / tau_D = {tau_A/tau_D:.1f}")

# ================================================================
# Test 1: Symmetric coupling (both negative)
# ================================================================
print("\n" + "=" * 65)
print("TEST 1: Symmetric coupling (k_DA = k_AD = -1)")
print("Both directions: stretching one compresses the other")
print("=" * 65)

for N_mol in [2, 4, 6, 10]:
    Q = build_swap_DA(N_mol)
    J = build_water_chain_jacobian(N_mol, tau_D, tau_A, k_DA=-1.0, k_AD=-1.0)
    res = palindrome_residual(J, Q)
    ev = eigvals(J)
    n_osc = np.sum(np.abs(ev.imag) > 1e-10)
    print(f"  N={N_mol:2d}:  residual = {res:.4f}  n_osc = {n_osc}")


# ================================================================
# Test 2: Antisymmetric coupling (k_DA = +1, k_AD = -1)
# ================================================================
print("\n" + "=" * 65)
print("TEST 2: Antisymmetric coupling (k_DA = +1, k_AD = -1)")
print("Donor->Acceptor positive, Acceptor->Donor negative")
print("Like Dale's Law: D excites, A inhibits")
print("=" * 65)

for N_mol in [2, 4, 6, 10]:
    Q = build_swap_DA(N_mol)
    J = build_water_chain_jacobian(N_mol, tau_D, tau_A, k_DA=+1.0, k_AD=-1.0)
    res = palindrome_residual(J, Q)
    ev = eigvals(J)
    n_osc = np.sum(np.abs(ev.imag) > 1e-10)
    print(f"  N={N_mol:2d}:  residual = {res:.4f}  n_osc = {n_osc}")


# ================================================================
# Test 3: Physical coupling from H-bond anharmonicity
# ================================================================
print("\n" + "=" * 65)
print("TEST 3: Physical H-bond model")
print("The Lippincott-Schroeder potential gives asymmetric coupling")
print("=" * 65)

print("""
In a real hydrogen bond O-H...O:
  - The O-H bond is strong (covalent, ~460 kJ/mol)
  - The H...O bond is weak (hydrogen bond, ~20 kJ/mol)
  - Stretching O-H pushes H toward O (donor -> acceptor: POSITIVE)
  - Stretching H...O pulls H away from O (acceptor -> donor: NEGATIVE?)

The asymmetry comes from the ENERGY RATIO:
  k_DA / k_AD = -(force constant ratio) * (tau_A / tau_D)

For the palindrome condition W[Q(i),Q(j)] = -(tau_Q(i)/tau_i) * W[i,j]:
  k_AD = -(tau_A / tau_D) * k_DA

With tau_A/tau_D = 4.86:  k_AD = -4.86 * k_DA
""")

ratio = tau_A / tau_D
print(f"  tau_A / tau_D = {ratio:.2f}")
print(f"  For exact palindrome: k_AD = -{ratio:.2f} * k_DA")
print()

# Test with EXACT palindromic ratio
k_DA = 1.0
k_AD_exact = -ratio * k_DA

print(f"  Exact palindromic coupling: k_DA = {k_DA:.2f}, k_AD = {k_AD_exact:.2f}")

for N_mol in [2, 4, 6, 10, 20]:
    Q = build_swap_DA(N_mol)
    J = build_water_chain_jacobian(N_mol, tau_D, tau_A,
                                    k_DA=k_DA, k_AD=k_AD_exact)
    res = palindrome_residual(J, Q)
    ev = eigvals(J)
    n_osc = np.sum(np.abs(ev.imag) > 1e-10)

    # Character swap check
    e_idx = np.arange(0, 2*N_mol, 2)  # donor indices
    i_idx = np.arange(1, 2*N_mol, 2)  # acceptor indices

    _, evec = np.linalg.eig(J)
    rates = [-ev[k].real for k in range(len(ev))]
    center = (min(rates) + max(rates)) / 2

    print(f"  N={N_mol:2d}:  residual = {res:.2e}  n_osc = {n_osc:2d}  "
          f"{'EXACT' if res < 1e-10 else ''}")


# ================================================================
# Test 4: V-Effect - couple two palindromic water clusters
# ================================================================
print("\n" + "=" * 65)
print("TEST 4: V-Effect - couple two exact palindromic water chains")
print("=" * 65)

N_mol = 5  # 5 molecules each
k_DA = 1.0
k_AD = -ratio * k_DA

# Single chain
J_single = build_water_chain_jacobian(N_mol, tau_D, tau_A, k_DA, k_AD)
ev_single = eigvals(J_single)
K_single = len(set(np.round(np.abs(ev_single.imag), 6)))
K_single = sum(1 for f in set(np.round(np.abs(ev_single.imag), 6)) if f > 1e-6)

# Two chains coupled through a shared hydrogen bond
N_total = 2 * N_mol
J_coupled = np.zeros((2 * N_total, 2 * N_total))

# Chain A: molecules 0..4
J_A = build_water_chain_jacobian(N_mol, tau_D, tau_A, k_DA, k_AD)
J_coupled[:2*N_mol, :2*N_mol] = J_A

# Chain B: molecules 5..9
J_B = build_water_chain_jacobian(N_mol, tau_D, tau_A, k_DA, k_AD, alpha=1.0)
J_coupled[2*N_mol:, 2*N_mol:] = J_B

# Coupling: last molecule of A donates to first molecule of B
for coupling in [0.0, 0.01, 0.05, 0.1, 0.3, 0.5]:
    J_test = J_coupled.copy()
    # D of last A -> A of first B
    J_test[2*N_mol + 1, 2*(N_mol-1)] = coupling * k_DA / tau_A
    # A of first B -> D of last A
    J_test[2*(N_mol-1), 2*N_mol + 1] = coupling * k_AD / tau_D

    ev_c = eigvals(J_test)
    K_coupled = sum(1 for f in set(np.round(np.abs(ev_c.imag), 6)) if f > 1e-6)

    corr_freqs = set()
    for i in range(len(ev_c)):
        for j in range(i, len(ev_c)):
            f = abs((ev_c[i] + ev_c[j]).imag)
            if f > 1e-6:
                corr_freqs.add(round(f, 6))
    K_corr = len(corr_freqs)

    new = K_corr - 2 * K_single if K_single > 0 else K_corr
    v_marker = "V-EFFECT!" if new > 0 and coupling > 0 else ""

    print(f"  coupling={coupling:.2f}:  K_single={K_single}  K_coupled={K_coupled}  "
          f"K_corr={K_corr}  new={new:+d}  {v_marker}")


# ================================================================
# Test 5: RING topology (closed chain, every molecule donates AND accepts)
# ================================================================
print("\n" + "=" * 65)
print("TEST 5: Ring topology (fixes the missing Q-partner)")
print("=" * 65)

def build_water_ring_jacobian(N_mol, tau_D, tau_A, k_DA, k_AD, alpha=1.0):
    """Ring of water molecules. Each molecule donates to the next
    AND accepts from the previous. Periodic boundary: last -> first.

    This ensures every D-A coupling has a Q-partnered A-D coupling.
    """
    N = 2 * N_mol
    J = np.zeros((N, N))

    for i in range(N_mol):
        d_idx = 2 * i
        a_idx = 2 * i + 1
        next_mol = (i + 1) % N_mol
        prev_mol = (i - 1) % N_mol

        # Self-decay
        J[d_idx, d_idx] = -1.0 / tau_D
        J[a_idx, a_idx] = -1.0 / tau_A

        # Intramolecular: D_i <-> A_i (within same molecule)
        J[d_idx, a_idx] = alpha * k_DA / tau_D
        J[a_idx, d_idx] = alpha * k_AD / tau_A

        # Intermolecular forward: D_i -> A_{next} (donate to next)
        next_a = 2 * next_mol + 1
        J[next_a, d_idx] = alpha * k_DA / tau_A

        # Intermolecular forward back-coupling: A_{next} -> D_i
        J[d_idx, next_a] = alpha * k_AD / tau_D

        # Intermolecular backward: A_i -> D_{prev} (accept from previous)
        prev_d = 2 * prev_mol
        J[prev_d, a_idx] = alpha * k_AD / tau_D

        # Intermolecular backward back-coupling: D_{prev} -> A_i
        J[a_idx, prev_d] = alpha * k_DA / tau_A

    return J


print(f"\n  Exact palindromic coupling: k_DA = {k_DA:.2f}, k_AD = {k_AD_exact:.2f}")
print(f"  Ring topology: every molecule donates AND accepts\n")

print(f"  {'N':>4s}  {'residual':>10s}  {'n_osc':>6s}  {'palindrome?':>12s}")
print(f"  {'-'*40}")

for N_mol in [3, 4, 5, 6, 8, 10, 20]:
    Q = build_swap_DA(N_mol)
    J = build_water_ring_jacobian(N_mol, tau_D, tau_A,
                                   k_DA=k_DA, k_AD=k_AD_exact)
    res = palindrome_residual(J, Q)
    ev = eigvals(J)
    n_osc = np.sum(np.abs(ev.imag) > 1e-10)
    status = "EXACT" if res < 1e-8 else ("close" if res < 0.1 else "no")
    print(f"  {N_mol:4d}  {res:10.2e}  {n_osc:6d}  {status:>12s}")


# ================================================================
# Test 6: V-Effect in palindromic water rings
# ================================================================
print("\n" + "=" * 65)
print("TEST 6: V-Effect - couple two palindromic water rings")
print("=" * 65)

N_mol = 5
J_A = build_water_ring_jacobian(N_mol, tau_D, tau_A, k_DA, k_AD_exact)
ev_A = eigvals(J_A)
K_single_act = sum(1 for f in set(np.round(np.abs(ev_A.imag), 8)) if f > 1e-6)

# Correlation frequencies of single ring
corr_single = set()
for i in range(len(ev_A)):
    for j in range(i, len(ev_A)):
        f = abs((ev_A[i] + ev_A[j]).imag)
        if f > 1e-6:
            corr_single.add(round(f, 8))
K_single_corr = len(corr_single)

# Two rings coupled through one hydrogen bond
N_total_modes = 2 * (2 * N_mol)
J_coupled_base = np.zeros((N_total_modes, N_total_modes))

J_B = build_water_ring_jacobian(N_mol, tau_D, tau_A, k_DA, k_AD_exact)

offset = 2 * N_mol
J_coupled_base[:offset, :offset] = J_A
J_coupled_base[offset:, offset:] = J_B

print(f"\n  Single ring (N={N_mol}): K_act={K_single_act}, K_corr={K_single_corr}")
print(f"\n  {'coupling':>10s}  {'K_act':>6s}  {'K_corr':>7s}  {'new_corr':>9s}  {'V-Effect?':>10s}")
print(f"  {'-'*50}")

for coupling in [0.0, 0.001, 0.01, 0.05, 0.1, 0.3]:
    J_test = J_coupled_base.copy()
    # Couple: D of last molecule in ring A -> A of first in ring B
    d_last_A = 2 * (N_mol - 1)
    a_first_B = offset + 1

    J_test[a_first_B, d_last_A] = coupling * k_DA / tau_A
    J_test[d_last_A, a_first_B] = coupling * k_AD_exact / tau_D

    ev_c = eigvals(J_test)
    K_act = sum(1 for f in set(np.round(np.abs(ev_c.imag), 8)) if f > 1e-6)

    corr_c = set()
    for i in range(len(ev_c)):
        for j in range(i, len(ev_c)):
            f = abs((ev_c[i] + ev_c[j]).imag)
            if f > 1e-6:
                corr_c.add(round(f, 8))
    K_corr = len(corr_c)

    new_corr = K_corr - 2 * K_single_corr
    v_marker = "V-EFFECT!" if new_corr > 0 and coupling > 0 else ""

    print(f"  {coupling:10.3f}  {K_act:6d}  {K_corr:7d}  {new_corr:+9d}  {v_marker:>10s}")


# ================================================================
# Test 7: Character swap in palindromic water ring
# ================================================================
print("\n" + "=" * 65)
print("TEST 7: Character swap (D-dominant <-> A-dominant)")
print("=" * 65)

N_mol = 6
J = build_water_ring_jacobian(N_mol, tau_D, tau_A, k_DA, k_AD_exact)
ev, evec = np.linalg.eig(J)
Q = build_swap_DA(N_mol)
res = palindrome_residual(J, Q)

d_idx = np.arange(0, 2*N_mol, 2)  # donor mode indices
a_idx = np.arange(1, 2*N_mol, 2)  # acceptor mode indices

# Find palindromic pairs and check character swap
rates = [-ev[k].real for k in range(len(ev))]
center = (min(rates) + max(rates)) / 2
tol = 0.03 * (max(rates) - min(rates))

pairs = []
used = set()
for k in range(len(ev)):
    if k in used:
        continue
    target = 2 * center - rates[k]
    for j in range(len(ev)):
        if j != k and j not in used and abs(rates[j] - target) < tol:
            pairs.append((k, j))
            used.add(k)
            used.add(j)
            break

print(f"\n  Residual: {res:.2e}  Pairs found: {len(pairs)}")
print(f"\n  {'Pair':>8s}  {'rate_k':>8s}  {'rate_k_':>8s}  "
      f"{'D(k)':>6s}  {'A(k)':>6s}  {'D(k_)':>6s}  {'A(k_)':>6s}  {'swap?':>6s}")
print(f"  {'-'*65}")

swap_errors = []
for k, kp in pairs:
    vec_k = evec[:, k]
    vec_kp = evec[:, kp]
    norm_k = np.sum(np.abs(vec_k)**2)
    norm_kp = np.sum(np.abs(vec_kp)**2)

    d_k = np.sum(np.abs(vec_k[d_idx])**2) / norm_k if norm_k > 0 else 0
    a_k = np.sum(np.abs(vec_k[a_idx])**2) / norm_k if norm_k > 0 else 0
    d_kp = np.sum(np.abs(vec_kp[d_idx])**2) / norm_kp if norm_kp > 0 else 0
    a_kp = np.sum(np.abs(vec_kp[a_idx])**2) / norm_kp if norm_kp > 0 else 0

    swap_err = abs(d_k - a_kp) + abs(a_k - d_kp)
    swap_errors.append(swap_err)
    swap_ok = "YES" if swap_err < 0.3 else "no"

    print(f"  {k:3d},{kp:<3d}  {rates[k]:8.1f}  {rates[kp]:8.1f}  "
          f"{d_k:6.3f}  {a_k:6.3f}  {d_kp:6.3f}  {a_kp:6.3f}  {swap_ok:>6s}")

if swap_errors:
    print(f"\n  Mean swap error: {np.mean(swap_errors):.4f}")


# ================================================================
# Summary
# ================================================================
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
