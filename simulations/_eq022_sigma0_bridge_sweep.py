"""sigma_0 bridge sweep across (c, N): Item 3 asymptote + g_eff_E·sqrt(3/8) bridge.

Numerical verification recommended by the J-synthesis (2026-05-07-j-as-coupling-strength.md):
sigma_0(c, N) is the inter-channel SVD-top singular value of V_inter, the structural bridge
between the F2b OBC dispersion E_k = 2J cos(pi k/(N+1)) and the per-bond g_eff(N, b).

Two empirical hypotheses to test:

1. **Asymptotic (Item 3 closed form):** sigma_0(c, N -> infty) -> 2*sqrt(2*(c-1)).
   The dimensionless ratio sigma_0/sqrt(2(c-1)) should converge monotonically from below to
   2.0 with growing N, for each c. Q_EP -> 1/sqrt(2(c-1)) follows.

2. **Bridge to g_eff_E (c=2):** g_eff_Endpoint(N) ~ sigma_0(c=2, N) * sqrt(3/8).
   Today's Direction (alpha) attempt found Delta <= 0.01 for N >= 6, but Delta = 0.063 at
   N=5. Verify across more N, check higher c. Empirical g_eff_E from pinned PolarityInheritanceLink
   witnesses: g_eff_E = 4.39382 / Q_peak_Endpoint (composition reading).

Usage: python simulations/_eq022_sigma0_bridge_sweep.py

Output: machine-parseable table on stdout (N, c, sigma_0, sigma_0/sqrt(2(c-1)),
sigma_0*sqrt(3/8), g_eff_E_empirical, Delta_sqrt38_bridge, Delta_asymptote).

This is a research-attempt script per simulations/_*.py convention; produces verification
data, the synthesis writeup is the place for interpretation. The script computes every number
directly via framework primitives (no hardcoded sigma_0 values).
"""
from __future__ import annotations

import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import math
from pathlib import Path

import numpy as np
from numpy.linalg import svd

# Make framework importable when run from repo root.
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))

import framework as fw  # noqa: E402  (validates the framework cockpit is reachable)
from framework.coherence_block import block_L_split_xy, hd_channel_basis  # noqa: E402


# Empirical anchors from PolarityInheritanceLink.cs witnesses (c=2, N=5..8, gamma_0=0.05):
# Q_peak_Endpoint values (4-decimal precision, pinned from the live C2HwhmRatio pipeline).
# g_eff_E = 4.39382 / Q_peak_Endpoint follows from BareDoubledPtfXPeak (= 2.196910) universal
# composition: Q_peak = 2.196910 * Q_EP = 2.196910 * 2 / g_eff = 4.39382 / g_eff.
EMPIRICAL_Q_PEAK_ENDPOINT = {
    5: 2.5008,
    6: 2.5470,
    7: 2.5299,
    8: 2.5145,
}
BARE_DOUBLED_PTF_X_PEAK = 2.196910  # Tier 1 derived universal (C2HwhmRatio)


def empirical_g_eff_endpoint(N):
    """g_eff_Endpoint(N) extracted by inverting Q_peak via BareDoubledPtfXPeak composition."""
    Q_peak = EMPIRICAL_Q_PEAK_ENDPOINT.get(N)
    if Q_peak is None:
        return None
    return 2.0 * BARE_DOUBLED_PTF_X_PEAK / Q_peak


def build_v_inter(N, n, gamma_0, hd1, hd2):
    """Compute V_inter = P_HD1^dagger * M_H_total * P_HD2 for the (n, n+1) coherence block.

    M_H_total = sum_b M_H_per_bond[b] is the J-coefficient block of L. The HD-channel
    projectors P_HD1, P_HD2 isolate the popcount(p ^ q) = HD1 and = HD2 subspaces.

    Returns the rectangular V_inter matrix (n_HD1, n_HD2) and a list of singular values
    sorted descending.
    """
    D, M_H_per_bond, P_n, P_np1 = block_L_split_xy(N, n, gamma_0)
    M_H_total = np.zeros_like(D)
    for M_b in M_H_per_bond:
        M_H_total += M_b

    P, HDs = hd_channel_basis(N, n)
    # We need full HD-subspace projectors (not the channel-uniform 1D projectors from
    # hd_channel_basis). Reconstruct: project onto every (p,q) basis pair with HD = hd_target.
    Mtot = D.shape[0]
    Mnp1 = len(P_np1)

    def hd_subspace_proj(hd_target):
        """Build an (Mtot, n_hd) orthonormal projector onto basis pairs (p,q) with
        Hamming-distance hd_target. Each column is a single computational-basis vector
        of the block."""
        cols = []
        for ip, p in enumerate(P_n):
            for iq, q in enumerate(P_np1):
                if bin(p ^ q).count("1") == hd_target:
                    e = np.zeros(Mtot, dtype=complex)
                    e[ip * Mnp1 + iq] = 1.0
                    cols.append(e)
        if not cols:
            return None
        return np.column_stack(cols)

    P_hd1 = hd_subspace_proj(hd1)
    P_hd2 = hd_subspace_proj(hd2)
    if P_hd1 is None or P_hd2 is None:
        return None, None

    V_inter = P_hd1.conj().T @ M_H_total @ P_hd2
    s = svd(V_inter, compute_uv=False)
    return V_inter, np.array(s, dtype=float)


def sweep_row(c, N, gamma_0=0.05, hd1=1, hd2=3):
    """Compute one row of the sweep: sigma_0 + derived ratios for (c, N)."""
    n = c - 1  # smallest n giving chromaticity c (matches WitnessCache convention)
    if N < 2 * c - 1:
        return None
    V_inter, s = build_v_inter(N, n, gamma_0, hd1=hd1, hd2=hd2)
    if s is None or len(s) == 0:
        return None
    sigma_0 = s[0]
    asymptote = 2.0 * math.sqrt(2.0 * (c - 1))
    normalised_ratio = sigma_0 / math.sqrt(2.0 * (c - 1))
    sqrt_3_8 = math.sqrt(3.0 / 8.0)
    bridge_predicted = sigma_0 * sqrt_3_8
    g_eff_E_empirical = empirical_g_eff_endpoint(N) if c == 2 else None
    delta_bridge = (bridge_predicted - g_eff_E_empirical) if g_eff_E_empirical is not None else None
    delta_asymptote = sigma_0 - asymptote
    block_dim = V_inter.shape[0] * V_inter.shape[1] // (V_inter.shape[1])  # purely informational
    return dict(
        c=c,
        N=N,
        block_dim_hd1=V_inter.shape[0],
        block_dim_hd2=V_inter.shape[1],
        sigma_0=sigma_0,
        asymptote=asymptote,
        normalised_ratio=normalised_ratio,
        bridge_predicted=bridge_predicted,
        g_eff_E_empirical=g_eff_E_empirical,
        delta_bridge=delta_bridge,
        delta_asymptote=delta_asymptote,
    )


def main():
    print(f"# sigma_0 bridge sweep, gamma_0 = 0.05, sqrt(3/8) = {math.sqrt(3.0/8.0):.6f}")
    print()
    print(f"# (HD pair convention: HD1=1, HD2=3, slowest-pair coupling, matches "
          f"InterChannelSvd.Build(block, hd1: 1, hd2: 3) and SigmaZeroChromaticityScaling.cs)")
    print(f"# framework loaded: ChainSystem? {hasattr(fw, 'ChainSystem')}")
    print()

    # (c, N) grid: c=2 N=5..9, c=3 N=5..8, c=4 N=7..9.
    grid = []
    for N in (5, 6, 7, 8, 9):
        grid.append((2, N))
    for N in (5, 6, 7, 8):
        grid.append((3, N))
    for N in (7, 8, 9):
        grid.append((4, N))

    rows = []
    for c, N in grid:
        try:
            row = sweep_row(c, N)
            if row is None:
                continue
            rows.append(row)
        except Exception as exc:  # noqa: BLE001  (research script: log + continue)
            print(f"# SKIP c={c} N={N}: {type(exc).__name__}: {exc}")
            continue

    # Print machine-parseable table.
    header = ["c", "N", "n_HD1", "n_HD2", "sigma_0", "sigma_0/sqrt(2(c-1))",
              "sigma_0*sqrt(3/8)", "g_eff_E_emp", "Delta_bridge", "Delta_asymptote"]
    print("\t".join(header))
    for r in rows:
        line = [
            str(r["c"]),
            str(r["N"]),
            str(r["block_dim_hd1"]),
            str(r["block_dim_hd2"]),
            f"{r['sigma_0']:.6f}",
            f"{r['normalised_ratio']:.6f}",
            f"{r['bridge_predicted']:.6f}",
            f"{r['g_eff_E_empirical']:.6f}" if r["g_eff_E_empirical"] is not None else "N/A",
            f"{r['delta_bridge']:+.6f}" if r["delta_bridge"] is not None else "N/A",
            f"{r['delta_asymptote']:+.6f}",
        ]
        print("\t".join(line))

    print()
    print("# === Hypothesis check (A): sigma_0/sqrt(2(c-1)) -> 2.0 monotone in N ===")
    for c in (2, 3, 4):
        sub = [r for r in rows if r["c"] == c]
        if not sub:
            continue
        print(f"#   c={c}:")
        prev = None
        for r in sub:
            tag = "OK" if r["normalised_ratio"] >= (prev or 0) else "NON-MONOTONE"
            within_001 = abs(r["normalised_ratio"] - 2.0) <= 0.01
            print(f"#     N={r['N']}: ratio={r['normalised_ratio']:.4f} "
                  f"(Delta_to_2.0={r['normalised_ratio']-2.0:+.4f}, |Delta|<=0.01? "
                  f"{'YES' if within_001 else 'no'}, monotone? {tag})")
            prev = r["normalised_ratio"]

    print()
    print("# === Hypothesis check (B): sigma_0*sqrt(3/8) =? g_eff_E (c=2 only) ===")
    sub = [r for r in rows if r["c"] == 2 and r["g_eff_E_empirical"] is not None]
    for r in sub:
        within_001 = abs(r["delta_bridge"]) <= 0.01
        within_005 = abs(r["delta_bridge"]) <= 0.005
        print(f"#   N={r['N']}: predicted={r['bridge_predicted']:.4f}, "
              f"empirical={r['g_eff_E_empirical']:.4f}, "
              f"Delta={r['delta_bridge']:+.4f} "
              f"(|Delta|<=0.01? {'YES' if within_001 else 'no'}, "
              f"|Delta|<=0.005? {'YES' if within_005 else 'no'})")

    print()
    print("# === Hypothesis check (C): does sqrt(3/8) generalise to c=3, c=4? ===")
    print("#   No empirical g_eff_E table for c >= 3 in PolarityInheritanceLink.")
    print("#   Sanity: sigma_0*sqrt(3/8) at higher c lands in the asymptotic 2*sqrt(2(c-1))*sqrt(3/8) =")
    for c in (2, 3, 4):
        ref = 2.0 * math.sqrt(2.0 * (c - 1)) * math.sqrt(3.0 / 8.0)
        print(f"#     c={c}: 2*sqrt(2*{c-1})*sqrt(3/8) = sqrt({3*(c-1)}) = {ref:.4f}")
        for r in rows:
            if r["c"] == c:
                gap_to_asym_bridge = r["bridge_predicted"] - ref
                print(f"#       N={r['N']}: sigma_0*sqrt(3/8) = {r['bridge_predicted']:.4f} "
                      f"(gap to c-asymptote-bridge {ref:.4f}: {gap_to_asym_bridge:+.4f})")


if __name__ == "__main__":
    main()
