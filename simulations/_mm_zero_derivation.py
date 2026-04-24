#!/usr/bin/env python3
"""F75 Mirror-Pair MI closed form verified against C# brecher simulation peaks.

Derivation: for single-excitation pure state |psi> = sum_j c_j |1_j> with
mirror-symmetric amplitudes c_{N-1-j} = eta c_j (eta = +/- 1), the mutual
information between mirror-pair sites (ell, N-1-ell) at t=0 is

    MI = 2 h(p_ell) - h(2 p_ell),   p_ell = |c_ell|^2

where h is binary Shannon entropy (base 2).

For the bonding modes |psi_k> = sqrt(2/(N+1)) sum_j sin(pi k (j+1)/(N+1)) |1_j>
(F65), p_ell(k, N) = (2/(N+1)) sin^2(pi k (ell+1)/(N+1)), and

    MM(0) = sum_{ell=0..floor(N/2)-1} [2 h(p_ell) - h(2 p_ell)]

Comparison to PeakMM from the C# brecher scan at uniform J, gamma_0 = 0.05:
MM(0) matches simulation peak within 7% (decay envelope from dephasing +
Heisenberg same-parity mixing over t = 0 to 0.1).
"""
import math


def h(p: float) -> float:
    """Binary Shannon entropy in bits."""
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def bonding_mm_zero(N: int, k: int) -> float:
    """MM(0) closed form for bonding:k on an N-site uniform chain."""
    def p_ell(ell: int) -> float:
        return (2.0 / (N + 1)) * math.sin(math.pi * k * (ell + 1) / (N + 1)) ** 2
    total = 0.0
    for ell in range(N // 2):
        p = p_ell(ell)
        total += 2 * h(p) - h(2 * p)
    return total


def main() -> None:
    print("F75 Mirror-Pair MI(0) closed form for bonding:k modes")
    print("    MM(0) = sum_ell [2 h(p_ell) - h(2 p_ell)]")
    print("    p_ell = (2/(N+1)) sin^2(pi k (ell+1)/(N+1))")
    print()

    for N in (5, 7, 9):
        print(f"=== N={N} ===")
        for k in (1, 2, 3):
            mm0 = bonding_mm_zero(N, k)
            print(f"  bonding:{k}  MM(0) = {mm0:.4f}")
        print()

    print("=== Comparison to C# brecher PeakMM (uniform J, gamma_0 = 0.05) ===")
    sim_peaks = {
        (5, 1): 0.789, (5, 2): 1.241, (5, 3): 0.865,
        (7, 1): 0.801, (7, 2): 1.090, (7, 3): 0.819,
        (9, 1): 0.830, (9, 2): 1.049, (9, 3): 0.829,
    }
    print(f"{'N':>3} {'k':>3} {'MM(0) analytic':>15} {'PeakMM sim':>12} {'sim/analytic':>14}")
    for (N, k), sim_mm in sorted(sim_peaks.items()):
        mm0 = bonding_mm_zero(N, k)
        ratio = sim_mm / mm0 if mm0 > 0 else float("nan")
        print(f"{N:>3} {k:>3} {mm0:>15.4f} {sim_mm:>12.4f} {ratio:>14.4f}")


if __name__ == "__main__":
    main()
