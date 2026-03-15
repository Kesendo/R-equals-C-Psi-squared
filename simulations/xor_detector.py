"""
XOR Detector: Finding what breaks palindromic symmetry
======================================================

The palindrome (Pi) is always there. Universal. Boring.
What's NOT palindromic is the signal. The information.
The dynamic space.

Method:
1. Compute full Liouvillian spectrum
2. Identify palindromic pairs (d + d' = 2*sum_gamma)
3. Extract the RESIDUAL: modes that don't pair
4. These unpaired modes ARE the information

Authors: Tom Wicht, Claude
Date: March 16, 2026
"""

import numpy as np
from itertools import product
from dataclasses import dataclass
from typing import List, Tuple

# Pauli matrices
I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

def tensor(*ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def site_op(op, site, N):
    ops = [I2]*N
    ops[site] = op
    return tensor(*ops)

def build_heisenberg_hamiltonian(N, J, topology="chain"):
    """Build Heisenberg XXX Hamiltonian."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    
    if topology == "chain":
        pairs = [(i, i+1) for i in range(N-1)]
    elif topology == "ring":
        pairs = [(i, (i+1) % N) for i in range(N)]
    elif topology == "star":
        pairs = [(0, i) for i in range(1, N)]
    else:
        raise ValueError(f"Unknown topology: {topology}")
    
    for i, j in pairs:
        for pauli in [sx, sy, sz]:
            H += J * site_op(pauli, i, N) @ site_op(pauli, j, N)
    return H


def build_liouvillian(H, gammas, N):
    """Build full Liouvillian superoperator with Z-dephasing."""
    d = 2**N
    d2 = d*d
    
    # Hamiltonian part: -i[H, rho] = -i(H x I - I x H^T)
    L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    
    # Dephasing: gamma_k * (Zk rho Zk - rho)
    for k in range(N):
        Zk = site_op(sz, k, N)
        ZkZ = np.kron(Zk, Zk.conj())
        L += gammas[k] * (ZkZ - np.eye(d2))
    
    return L

@dataclass
class PalindromePair:
    rate_d: float       # Decay rate d
    rate_partner: float # Decay rate 2*sum_gamma - d
    imag_d: float       # Imaginary part of d
    imag_partner: float # Imaginary part of partner
    quality: float      # How well they match (0=bad, 1=perfect)
    idx_d: int
    idx_partner: int

@dataclass 
class XORResult:
    """The residual: what's NOT palindromic."""
    N: int
    topology: str
    gammas: list
    sum_gamma: float
    total_modes: int
    paired_modes: int
    unpaired_modes: int
    pairs: List[PalindromePair]
    unpaired_eigenvalues: List[complex]  # THE XOR
    palindrome_fraction: float
    xor_fraction: float


def find_palindromic_pairs(eigenvalues, sum_gamma, tolerance=1e-6):
    """Find pairs where Re(d) + Re(d') = 2*sum_gamma."""
    n = len(eigenvalues)
    reals = np.real(eigenvalues)
    imags = np.imag(eigenvalues)
    
    used = set()
    pairs = []
    
    # Sort by real part for efficient matching
    sorted_idx = np.argsort(reals)
    
    for i in range(n):
        if i in used:
            continue
        
        d = reals[sorted_idx[i]]
        target = 2 * sum_gamma - d  # Expected partner
        
        # Find closest unused match
        best_j = None
        best_diff = tolerance
        
        for j in range(n):
            if j in used or j == i:
                continue
            diff = abs(reals[sorted_idx[j]] - target)
            if diff < best_diff:
                best_diff = diff
                best_j = j
        
        if best_j is not None:
            idx_i = sorted_idx[i]
            idx_j = sorted_idx[best_j]
            quality = 1.0 - best_diff / (abs(sum_gamma) + 1e-10)
            
            pairs.append(PalindromePair(
                rate_d=reals[idx_i],
                rate_partner=reals[idx_j],
                imag_d=imags[idx_i],
                imag_partner=imags[idx_j],
                quality=max(0, quality),
                idx_d=idx_i,
                idx_partner=idx_j
            ))
            used.add(i)
            used.add(best_j)
    
    # Unpaired modes = THE XOR
    unpaired_idx = [sorted_idx[i] for i in range(n) if i not in used]
    unpaired = eigenvalues[unpaired_idx] if unpaired_idx else np.array([])
    
    return pairs, unpaired

def analyze_xor(N, J, gammas, topology="chain", tolerance=1e-6):
    """Full XOR analysis: build system, find palindrome, extract residual."""
    
    H = build_heisenberg_hamiltonian(N, J, topology)
    L = build_liouvillian(H, gammas, N)
    eigenvalues = np.linalg.eigvals(L)
    
    sum_gamma = sum(gammas)
    
    # Only analyze non-zero eigenvalues (skip steady state)
    nonzero_mask = np.abs(eigenvalues) > 1e-10
    nonzero_eigs = eigenvalues[nonzero_mask]
    
    pairs, unpaired = find_palindromic_pairs(nonzero_eigs, -sum_gamma, tolerance)
    
    total = len(nonzero_eigs)
    paired = len(pairs) * 2
    n_unpaired = len(unpaired)
    
    return XORResult(
        N=N,
        topology=topology,
        gammas=list(gammas),
        sum_gamma=sum_gamma,
        total_modes=total,
        paired_modes=paired,
        unpaired_modes=n_unpaired,
        pairs=pairs,
        unpaired_eigenvalues=list(unpaired),
        palindrome_fraction=paired / total if total > 0 else 0,
        xor_fraction=n_unpaired / total if total > 0 else 0
    )


def print_xor_result(result):
    """Pretty print the XOR analysis."""
    print(f"\n{'='*60}")
    print(f"XOR ANALYSIS: N={result.N}, {result.topology}")
    print(f"gammas = {result.gammas}, sum_gamma = {result.sum_gamma}")
    print(f"{'='*60}")
    print(f"Total non-zero modes: {result.total_modes}")
    print(f"Palindromic pairs:    {result.paired_modes} ({result.palindrome_fraction:.1%})")
    print(f"UNPAIRED (XOR):       {result.unpaired_modes} ({result.xor_fraction:.1%})")
    
    if result.unpaired_modes > 0:
        print(f"\n--- THE XOR: Unpaired eigenvalues ---")
        for ev in result.unpaired_eigenvalues:
            print(f"  lambda = {ev.real:.6f} + {ev.imag:.6f}i")
    
    if len(result.pairs) <= 20:
        print(f"\n--- Palindromic pairs (d + d' = {2*(-result.sum_gamma):.4f}) ---")
        for p in sorted(result.pairs, key=lambda x: x.rate_d):
            summ = p.rate_d + p.rate_partner
            print(f"  {p.rate_d:+.6f} + {p.rate_partner:+.6f} = {summ:.6f}"
                  f"  (target: {2*(-result.sum_gamma):.6f})"
                  f"  quality: {p.quality:.4f}")

# ============================================================
# EXPERIMENTS
# ============================================================

if __name__ == "__main__":
    
    print("*" * 60)
    print("XOR DETECTOR: What breaks palindromic symmetry?")
    print("*" * 60)
    
    # ---------------------------------------------------------
    # EXP 1: Uniform dephasing (should be fully palindromic)
    # ---------------------------------------------------------
    print("\n>>> EXP 1: Uniform dephasing (symmetric system)")
    print("    Expectation: 100% palindromic, 0% XOR")
    
    result = analyze_xor(N=3, J=1.0, 
                         gammas=[0.05, 0.05, 0.05],
                         topology="chain")
    print_xor_result(result)
    
    # ---------------------------------------------------------
    # EXP 2: Non-uniform dephasing (broken symmetry)
    # ---------------------------------------------------------
    print("\n>>> EXP 2: Non-uniform dephasing (asymmetric noise)")
    print("    One qubit dephases faster. Does this break palindrome?")
    
    result = analyze_xor(N=3, J=1.0,
                         gammas=[0.05, 0.05, 0.20],
                         topology="chain")
    print_xor_result(result)

    # ---------------------------------------------------------
    # EXP 3: Different topologies, same parameters
    # ---------------------------------------------------------
    print("\n>>> EXP 3: Topology comparison (same gammas)")
    print("    Does topology change the XOR?")
    
    for topo in ["chain", "ring", "star"]:
        result = analyze_xor(N=3, J=1.0,
                             gammas=[0.05, 0.05, 0.05],
                             topology=topo)
        print(f"\n  {topo:6s}: {result.palindrome_fraction:.1%} palindromic, "
              f"{result.unpaired_modes} unpaired")
    
    # ---------------------------------------------------------
    # EXP 4: Star with asymmetric coupling (2:1 optimal)
    # ---------------------------------------------------------
    print("\n>>> EXP 4: Star topology - symmetric vs asymmetric")
    print("    The 2:1 ratio that optimizes QST")
    
    # Symmetric: all same gamma
    r_sym = analyze_xor(N=3, J=1.0,
                        gammas=[0.05, 0.05, 0.05],
                        topology="star")
    
    # Asymmetric: mediator different
    r_asym = analyze_xor(N=3, J=1.0,
                         gammas=[0.10, 0.05, 0.05],
                         topology="star")
    
    print(f"  Symmetric gammas:  XOR = {r_sym.unpaired_modes} modes")
    print(f"  Asymmetric gammas: XOR = {r_asym.unpaired_modes} modes")
    print_xor_result(r_asym)

    # ---------------------------------------------------------
    # EXP 5: Gamma sweep - how does XOR change with noise?
    # ---------------------------------------------------------
    print("\n>>> EXP 5: Gamma sweep (N=3 chain)")
    print("    How does XOR fraction change with dephasing strength?")
    print(f"  {'gamma':>8} {'paired':>8} {'unpaired':>10} {'XOR%':>8}")
    print("  " + "-" * 40)
    
    for g in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0]:
        result = analyze_xor(N=3, J=1.0,
                             gammas=[g, g, g],
                             topology="chain")
        print(f"  {g:>8.3f} {result.paired_modes:>8} {result.unpaired_modes:>10}"
              f" {result.xor_fraction:>8.1%}")
    
    # ---------------------------------------------------------
    # EXP 6: System size scaling
    # ---------------------------------------------------------
    print("\n>>> EXP 6: System size scaling")
    print("    Does XOR grow with N?")
    print(f"  {'N':>4} {'total':>8} {'paired':>8} {'unpaired':>10} {'XOR%':>8}")
    print("  " + "-" * 44)
    
    for N in [2, 3, 4, 5]:
        gammas = [0.05] * N
        result = analyze_xor(N=N, J=1.0, gammas=gammas, topology="chain")
        print(f"  {N:>4} {result.total_modes:>8} {result.paired_modes:>8}"
              f" {result.unpaired_modes:>10} {result.xor_fraction:>8.1%}")
    
    # ---------------------------------------------------------
    # SUMMARY
    # ---------------------------------------------------------
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("If XOR = 0 everywhere: palindrome is perfect, no dynamic space")
    print("If XOR > 0 somewhere: THAT is where information lives")
    print("If XOR depends on topology: structure creates information")
    print("If XOR depends on gamma: noise creates information")
    print("=" * 60)
