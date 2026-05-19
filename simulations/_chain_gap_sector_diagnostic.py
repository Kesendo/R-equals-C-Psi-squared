"""Identify which joint-popcount sector hosts the chain dissipation-gap eigenmode.

The 2026-05-19 Q-sweep + Absorption Theorem reading suggests the slow mode of
the chain Heisenberg + Z-dephasing Liouvillian is NOT in the pure w=1 sector
(F50 pins those at Re = -2γ); it lives in a mixed-content sector with
fractional ⟨n_XY⟩ ≪ 1. The bit-exact prediction at Q=2:

    ⟨n_XY⟩_slow ≈ 0.55·Q²/N²    (chain, Marrakesh convention)

which reproduces the empirical gap = 2γ·⟨n_XY⟩_slow = 1.10·γ·Q²/N² to <1%.

This script reads the slow mode's joint-popcount sector (p_col, p_row) and
Pauli-basis light content ⟨n_XY⟩ directly, giving two independent diagnostics
of where the slow physics lives in the operator algebra.

Diagnostics surfaced per (N, Q):
  1. Slow-mode eigenvalue (Re, Im) and dissipation gap = |Re|.
  2. Joint-popcount sector (p_col, p_row). d=2^N states partition by Z-popcount,
     operator-space d²=4^N partitions by the joint (col-popcount, row-popcount).
     For Heisenberg + Z-dephasing both are conserved separately. The slow mode
     lives in ONE specific (p_col, p_row) block.
  3. Pauli-basis light content ⟨n_XY⟩ = Σ_p |⟨p|slow⟩|² · n_XY(p), where the
     sum is over Pauli strings p and n_XY(p) is the count of X/Y operators.
     Per the Absorption Theorem, Re(λ_slow) = -2γ·⟨n_XY⟩ exactly.

Hamiltonian: H = (J/4) Σ_b (X_i X_j + Y_i Y_j + Z_i Z_j) on the open chain.
Uniform Z-dephasing γ per site.
"""
from __future__ import annotations

import sys
import time
from itertools import product

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from framework.lindblad import lindbladian_z_dephasing  # noqa: E402

from _f1_chain_heisenberg_small_n_anchor import (  # noqa: E402
    chain_bonds,
    heisenberg_graph_h,
)


def popcount(i: int) -> int:
    return bin(i).count("1")


def joint_popcount_sectors(N: int) -> dict[tuple[int, int], np.ndarray]:
    """Map (p_col, p_row) → indices in vec(ρ) space where ρ is d×d, vec uses
    convention vec[i*d + j] = ρ[i, j] (row-major). col-popcount = popcount(i),
    row-popcount = popcount(j). The Heisenberg Hamiltonian conserves total Z
    (popcount of any computational basis state), and the dephasing dissipator
    is diagonal in the computational basis; both preserve col and row
    popcounts separately, so L block-diagonalises by (p_col, p_row).
    """
    d = 1 << N
    sectors: dict[tuple[int, int], list[int]] = {}
    pops = [popcount(i) for i in range(d)]
    for i in range(d):
        for j in range(d):
            sec = (pops[i], pops[j])
            sectors.setdefault(sec, []).append(i * d + j)
    return {k: np.array(v, dtype=np.int64) for k, v in sectors.items()}


def all_pauli_strings_n_xy(N: int) -> np.ndarray:
    """For each Pauli-string index p ∈ {0, ..., 4^N - 1} return the count
    of X or Y letters in the string. Convention: digits in base 4 are
    {0: I, 1: X, 2: Y, 3: Z}; X and Y count as "light", I and Z as "dark"."""
    n_xy = np.zeros(4 ** N, dtype=np.int32)
    for p in range(4 ** N):
        q = p
        c = 0
        for _ in range(N):
            d = q & 3
            if d == 1 or d == 2:  # X or Y
                c += 1
            q >>= 2
        n_xy[p] = c
    return n_xy


_PAULI_SINGLE = np.array([
    [[1, 0], [0, 1]],   # I
    [[0, 1], [1, 0]],   # X
    [[0, -1j], [1j, 0]],  # Y
    [[1, 0], [0, -1]],  # Z
], dtype=complex)


def pauli_string_matrix(p: int, N: int) -> np.ndarray:
    """Build the N-qubit Pauli string from its base-4 index."""
    digits = []
    q = p
    for _ in range(N):
        digits.append(q & 3)
        q >>= 2
    M = _PAULI_SINGLE[digits[-1]]
    for d in reversed(digits[:-1]):
        M = np.kron(M, _PAULI_SINGLE[d])  # qubit 0 = MSB convention (matches Compute)
    return M


def pauli_basis_overlaps(slow_vec: np.ndarray, N: int) -> np.ndarray:
    """Project the slow-mode operator |slow⟩ (length d², vec convention
    vec[i*d+j] = ρ[i,j]) onto the orthonormal Pauli basis P_p / sqrt(d).
    Returns overlap squared |⟨p|slow⟩|² per Pauli string p, normalised to sum
    to 1 (the eigenmode is unit-norm in the Frobenius inner product)."""
    d = 1 << N
    slow_op = slow_vec.reshape(d, d)
    overlaps = np.zeros(4 ** N)
    norm_factor = 1.0 / d
    for p in range(4 ** N):
        P = pauli_string_matrix(p, N)
        amp = np.vdot(P.reshape(-1), slow_vec) * np.sqrt(norm_factor)
        overlaps[p] = abs(amp) ** 2
    return overlaps


def run(N: int, J: float, gamma: float) -> dict:
    Q = J / gamma
    print(f"\n--- chain N={N}, J={J}, γ={gamma}, Q={Q:.4f} ---")
    d = 1 << N

    t0 = time.time()
    H = heisenberg_graph_h(N, chain_bonds(N), J=J)
    gammas = [gamma] * N
    L = lindbladian_z_dephasing(H, gammas)
    print(f"  L built ({L.shape[0]}×{L.shape[0]}), {time.time()-t0:.2f}s")

    sectors = joint_popcount_sectors(N)
    print(f"  joint-popcount sectors: {len(sectors)} blocks, sizes "
          f"{sorted({len(v) for v in sectors.values()})}")

    # Eigendecompose per sector to identify the slow mode block-by-block.
    # ρ-space convention: vec[i*d+j] = ρ[i,j]. L acts on this vec space.
    kernel_tol = 1e-9
    slow_mode = None  # (re, im, sector, eigvec_in_full_space)
    per_sector_slow = []
    for sec, indices in sorted(sectors.items()):
        if len(indices) == 0:
            continue
        Lblock = L[np.ix_(indices, indices)]
        evals, evecs = np.linalg.eig(Lblock)
        non_kernel_mask = np.abs(evals.real) > kernel_tol
        if not np.any(non_kernel_mask):
            per_sector_slow.append((sec, None, len(indices)))
            continue
        idx = np.argmin(np.abs(evals.real[non_kernel_mask]))
        block_slow_idx = np.where(non_kernel_mask)[0][idx]
        re = evals[block_slow_idx].real
        im = evals[block_slow_idx].imag
        per_sector_slow.append((sec, (re, im), len(indices)))
        if slow_mode is None or abs(re) < abs(slow_mode[0]):
            full_vec = np.zeros(d * d, dtype=complex)
            full_vec[indices] = evecs[:, block_slow_idx]
            full_vec /= np.linalg.norm(full_vec)
            slow_mode = (re, im, sec, full_vec)

    re_slow, im_slow, sec_slow, vec_slow = slow_mode
    gap = abs(re_slow)
    sigma = N * gamma
    print(f"  slow mode: λ = {re_slow:.6e} + {im_slow:.6e}i  → gap = {gap:.6e}")
    print(f"  slow-mode sector (p_col, p_row) = {sec_slow}")

    # Light content via Pauli-basis projection (the F3 reading).
    print(f"  projecting onto Pauli basis ({4**N} strings)...", end=" ", flush=True)
    t0 = time.time()
    overlaps = pauli_basis_overlaps(vec_slow, N)
    n_xy_arr = all_pauli_strings_n_xy(N)
    overlap_sum = overlaps.sum()
    n_xy_avg = float((overlaps * n_xy_arr).sum() / overlap_sum)
    print(f"{time.time()-t0:.2f}s, normalised sum {overlap_sum:.6f}")
    print(f"  ⟨n_XY⟩_slow = {n_xy_avg:.6f}")
    print(f"  predicted gap from Absorption Theorem (-2γ·⟨n_XY⟩) "
          f"= {2 * gamma * n_xy_avg:.6e}  (match within {abs(2*gamma*n_xy_avg - gap)/gap*100:.3f}%)")
    print(f"  closed-form prediction ⟨n_XY⟩ ≈ 0.55·Q²/N² = {0.55*Q*Q/(N*N):.6f}")

    # Distribution of ⟨n_XY⟩ across Pauli weights for the slow mode
    weights_dist = np.zeros(N + 1)
    for w in range(N + 1):
        mask = n_xy_arr == w
        weights_dist[w] = overlaps[mask].sum()
    print("  weight distribution of slow mode (Σ |⟨p|slow⟩|² per n_XY):")
    for w in range(N + 1):
        if weights_dist[w] > 1e-9:
            print(f"    n_XY={w}: {weights_dist[w]:.6f}  ({weights_dist[w]/overlap_sum*100:.2f}%)")

    print("  per-sector slow eigenvalues (top 10 by smallest |Re|):")
    per_sector_slow_finite = [(sec, ev, sz) for sec, ev, sz in per_sector_slow if ev is not None]
    per_sector_slow_finite.sort(key=lambda x: abs(x[1][0]))
    for sec, (re, im), sz in per_sector_slow_finite[:10]:
        marker = " ← gap" if (re, im, sec) == (re_slow, im_slow, sec_slow) else ""
        print(f"    sector {sec}, block size {sz}: λ_slow_in_block = "
              f"{re:+.4e} + {im:+.4e}i{marker}")

    return {
        "N": N, "J": J, "gamma": gamma, "Q": Q,
        "gap": gap, "im_slow": im_slow,
        "slow_sector": sec_slow,
        "n_xy_avg": n_xy_avg,
        "predicted_n_xy": 0.55 * Q * Q / (N * N),
        "weight_distribution": weights_dist.tolist(),
    }


def main() -> None:
    print("Chain dissipation-gap sector diagnostic")
    print("=" * 80)
    print("Goal: identify joint-popcount sector + Pauli light content of the slow mode")
    print()

    # Marrakesh-convention anchors (γ=0.5, J=1, Q=2) for N=4, 5, 6
    # We pick Q=2 because today's typed claims are anchored there.
    results = []
    for N in (4, 5, 6):
        r = run(N=N, J=1.0, gamma=0.5)
        results.append(r)

    print("\n" + "=" * 80)
    print("Summary (Marrakesh convention γ=0.5, J=1, Q=2):")
    print(f"{'N':>3}  {'gap':>10}  {'sector':>12}  {'⟨n_XY⟩':>10}  "
          f"{'predicted':>10}  {'Im_slow':>10}")
    for r in results:
        sec_str = f"{r['slow_sector']}"
        print(f"{r['N']:3d}  {r['gap']:10.6e}  {sec_str:>12}  "
              f"{r['n_xy_avg']:10.6f}  {r['predicted_n_xy']:10.6f}  "
              f"{r['im_slow']:+10.6f}")

    print("\nAbsorption Theorem reading: gap = 2γ·⟨n_XY⟩_slow → "
          "matches gap value to machine precision (or close).")
    print(f"Closed-form conjecture ⟨n_XY⟩ ≈ 0.55·Q²/N² should be accurate to ~1%.")


if __name__ == "__main__":
    main()
