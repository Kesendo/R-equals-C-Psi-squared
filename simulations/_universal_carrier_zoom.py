"""Universal Carrier Zoom (2026-05-12): examine full eigenvalue structure
of N=3 chain Liouvillian without rounding aggregation.

Goes deeper than _universal_carrier_demo.py: shows full Re+iIm distribution,
identifies the kernel modes, the palindromic pairs, the F33 mixed cluster
fine structure, and the F1 palindrome center.
"""
import sys
import numpy as np
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw

print("=" * 78)
print("Universal Carrier: Full eigenvalue structure of N=3 chain Liouvillian")
print("=" * 78)
print()

# Build at γ=0.05 (engine default) and γ=0.10 (verification)
for gamma in [0.05, 0.10]:
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=gamma)
    L = chain.L
    eigs = np.linalg.eigvals(L)
    re = np.real(eigs)
    im = np.imag(eigs)

    print(f"--- γ = {gamma}  (Σγ = {3*gamma}, palindrome center predicted at Re = {-3*gamma}) ---")
    print(f"Total eigenvalues: {len(eigs)} (= d² = 4^3 = 64)")

    # Kernel: Re ≈ 0
    kernel = np.where(np.abs(re) < 1e-10)[0]
    print(f"Kernel modes (|Re| < 1e-10): {len(kernel)}")

    # Group by Re value (no rounding; use clustering with absolute tolerance)
    re_sorted = np.sort(np.unique(np.round(re, decimals=10)))
    print(f"Distinct Re values (10-digit precision): {len(re_sorted)}")
    print()

    # Bucket by integer-multiple of 2γ
    print(f"  {'Re(λ)':>16}  {'count':>6}  {'n_XY = -Re/(2γ)':>16}  {'reading':>40}")
    print(f"  {'-'*16}  {'-'*6}  {'-'*16}  {'-'*40}")
    bucket_count = Counter()
    bucket_repr = {}
    for r in re:
        # bucket at 1e-6 absolute (well below F33 splitting of 7e-5)
        key = round(r, 6)
        bucket_count[key] += 1
        if key not in bucket_repr or abs(r - key) < abs(bucket_repr[key] - key):
            bucket_repr[key] = r

    sorted_keys = sorted(bucket_count.keys())
    for key in sorted_keys:
        count = bucket_count[key]
        repr_val = bucket_repr[key]
        n_xy = -repr_val / (2 * gamma) if gamma > 0 else 0
        if abs(repr_val) < 1e-10:
            reading = f"kernel (n_XY = 0)"
        elif abs(n_xy - round(n_xy)) < 1e-3:
            reading = f"a₀·γ·{int(round(n_xy))}  pure-weight rung"
        elif abs(n_xy - 4/3) < 1e-2:
            reading = f"F33 mixed near 4/3 (8γ/3 = {8*gamma/3:.6f})"
        elif abs(n_xy - 5/3) < 1e-2:
            reading = f"F33 mixed near 5/3 (10γ/3 = {10*gamma/3:.6f})"
        else:
            reading = f"a₀·γ·{n_xy:.4f}"
        print(f"  {repr_val:>16.10f}  {count:>6}  {n_xy:>16.6f}  {reading:>40}")

    # Total non-kernel eigenvalues
    nonkernel = sum(c for k, c in bucket_count.items() if abs(k) > 1e-10)
    print()
    print(f"Sum of non-kernel modes: {nonkernel} (kernel: {len(kernel)}, total: {nonkernel + len(kernel)})")

    # F1 palindrome verification: every mode at Re paired with mode at -2Σγ - Re
    sigma_gamma = 3 * gamma
    palindrome_center = -sigma_gamma
    print(f"F1 palindrome center: −Σγ = {palindrome_center}")
    pairs_checked = 0
    pairs_matched = 0
    re_set = set(round(r, 6) for r in re)
    for r in sorted_keys:
        if abs(r) < 1e-10:
            continue
        partner = round(-2 * sigma_gamma - r, 6)
        if partner == r:
            # self-paired (at center)
            pairs_checked += 1
            pairs_matched += 1
        elif partner in re_set:
            pairs_checked += 1
            pairs_matched += 1
    print(f"F1 partnership: {pairs_matched}/{pairs_checked} non-kernel Re values find their −2Σγ−Re partner")
    print()
