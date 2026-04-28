#!/usr/bin/env python3
"""EQ-024 follow-up: why is F71-symmetry capacity-optimal?

Sub-question 2 closure showed F71-symmetric receivers achieve higher max
capacity than F71-breaking ones (12.41 vs 11.99 product; 10.26 vs 8.80
non-product), even though F71-symmetric is a strict subset of pure states.

Hypothesis: F71-symmetry structurally constrains the J-Jacobian into a
2+2 block decomposition (R-symmetric J-modes J_0+J_3, J_1+J_2 form one
block; R-antisymmetric J_0-J_3, J_1-J_2 form the other). This decomposition
allows one block's singular values to be large while the other's are
small — bimodal SV spectrum, favorable for waterfilling capacity.

F71-breaking receivers have no block structure; SVs are generic, more
uniform.

Empirical test: compute full SV spectra for ~30 samples per mode and
compare distribution of (sv_max, sv_min, sv_ratio = sv_max/sv_min,
"peakedness" measures). If F71-symmetric is systematically more bimodal
than F71-breaking, hypothesis confirmed.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

from _eq024_three_class_completeness import (
    N, d, compute_jacobian_J, waterfilling_capacity, bloch_state,
)
from _eq024_f71_breaking_capacity import (
    random_product_asymmetric, f71_asymmetry,
)


def random_f71_sym_product(rng):
    """F71-symmetric product state ψ = |a⟩|b⟩|c⟩|b⟩|a⟩."""
    th_a = rng.uniform(0, np.pi); ph_a = rng.uniform(0, 2 * np.pi)
    th_b = rng.uniform(0, np.pi); ph_b = rng.uniform(0, 2 * np.pi)
    th_c = rng.uniform(0, np.pi); ph_c = rng.uniform(0, 2 * np.pi)
    a = bloch_state(th_a, ph_a)
    b = bloch_state(th_b, ph_b)
    c = bloch_state(th_c, ph_c)
    psi = a
    for v in [b, c, b, a]:
        psi = np.kron(psi, v)
    return psi


def measure(psi):
    rho0 = np.outer(psi, psi.conj())
    rho0 = (rho0 + rho0.conj().T) / 2.0
    rho0 /= np.trace(rho0).real
    A = compute_jacobian_J(rho0)
    sv = np.linalg.svd(A, compute_uv=False)
    C = waterfilling_capacity(sv)
    return sv, C


def main():
    n_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 42

    rng_sym = np.random.default_rng(seed)
    rng_brk = np.random.default_rng(seed + 1)

    print(f"EQ-024 F71-optimality SVD analysis ({n_samples} samples per mode)")
    print(f"  Hypothesis: F71-sym SVDs bimodal (one big sv, one small),")
    print(f"  F71-breaking SVDs flat (all 4 sv comparable).")
    print()

    sym_data = []
    brk_data = []

    t0 = time.time()
    for i in range(n_samples):
        psi_s = random_f71_sym_product(rng_sym)
        sv_s, C_s = measure(psi_s)
        sym_data.append({'sv': sv_s.tolist(), 'C': C_s,
                         'asym': f71_asymmetry(psi_s)})

        psi_b = random_product_asymmetric(rng_brk)
        sv_b, C_b = measure(psi_b)
        brk_data.append({'sv': sv_b.tolist(), 'C': C_b,
                         'asym': f71_asymmetry(psi_b)})

        if (i + 1) % 5 == 0:
            print(f"  [{i + 1}/{n_samples}] elapsed {time.time() - t0:.1f}s")

    # Analysis
    sym_sv = np.array([d['sv'] for d in sym_data])  # (n, 4)
    brk_sv = np.array([d['sv'] for d in brk_data])

    sym_C = np.array([d['C'] for d in sym_data])
    brk_C = np.array([d['C'] for d in brk_data])

    sym_asym = np.array([d['asym'] for d in sym_data])
    brk_asym = np.array([d['asym'] for d in brk_data])

    print()
    print("=== F71-asymmetry sanity ===")
    print(f"  F71-sym samples:      asym range [{sym_asym.min():.4f}, {sym_asym.max():.4f}]  (should be ~0)")
    print(f"  F71-breaking samples: asym range [{brk_asym.min():.4f}, {brk_asym.max():.4f}]  (should be > 0)")

    print()
    print("=== Capacity ===")
    print(f"  F71-sym:      C range [{sym_C.min():.3f}, {sym_C.max():.3f}], mean {sym_C.mean():.3f}")
    print(f"  F71-breaking: C range [{brk_C.min():.3f}, {brk_C.max():.3f}], mean {brk_C.mean():.3f}")

    # SV statistics: per-sample sv_max, sv_min, ratio
    sym_sv_max = sym_sv[:, 0]
    sym_sv_min = sym_sv[:, -1]
    sym_ratio = sym_sv_max / np.maximum(sym_sv_min, 1e-12)

    brk_sv_max = brk_sv[:, 0]
    brk_sv_min = brk_sv[:, -1]
    brk_ratio = brk_sv_max / np.maximum(brk_sv_min, 1e-12)

    print()
    print("=== SV spectrum (sorted high to low: sv₁ ≥ sv₂ ≥ sv₃ ≥ sv₄) ===")
    print(f"  {'mode':>15s}  {'sv₁ mean':>9s}  {'sv₂ mean':>9s}  {'sv₃ mean':>9s}  {'sv₄ mean':>9s}  {'sv₁/sv₄ mean':>13s}")
    print(f"  {'F71-sym':>15s}  "
          f"{sym_sv[:, 0].mean():>9.3f}  {sym_sv[:, 1].mean():>9.3f}  "
          f"{sym_sv[:, 2].mean():>9.3f}  {sym_sv[:, 3].mean():>9.3f}  "
          f"{sym_ratio.mean():>13.2f}")
    print(f"  {'F71-breaking':>15s}  "
          f"{brk_sv[:, 0].mean():>9.3f}  {brk_sv[:, 1].mean():>9.3f}  "
          f"{brk_sv[:, 2].mean():>9.3f}  {brk_sv[:, 3].mean():>9.3f}  "
          f"{brk_ratio.mean():>13.2f}")

    print()
    print("=== SV ratio distribution (peakedness measure) ===")
    print(f"  F71-sym sv₁/sv₄:      min {sym_ratio.min():.2f}, max {sym_ratio.max():.2f}, "
          f"mean {sym_ratio.mean():.2f}, median {np.median(sym_ratio):.2f}")
    print(f"  F71-breaking sv₁/sv₄: min {brk_ratio.min():.2f}, max {brk_ratio.max():.2f}, "
          f"mean {brk_ratio.mean():.2f}, median {np.median(brk_ratio):.2f}")

    # Test the hypothesis quantitatively: is sv₁/sv₄ systematically larger for F71-sym?
    # Mann-Whitney style: count how many F71-sym ratios exceed F71-breaking ratios pairwise.
    pairs = 0
    sym_wins = 0
    for r_s in sym_ratio:
        for r_b in brk_ratio:
            pairs += 1
            if r_s > r_b:
                sym_wins += 1
    win_frac = sym_wins / pairs
    print()
    print(f"=== Hypothesis test ===")
    print(f"  P(F71-sym sv₁/sv₄ > F71-breaking sv₁/sv₄) ≈ {win_frac:.3f} "
          f"(0.5 = no effect, 1.0 = always sym more bimodal)")

    # Also: SV-skewness measure (sv₁² / Σ sv²)
    sym_skew = sym_sv[:, 0] ** 2 / np.sum(sym_sv ** 2, axis=1)
    brk_skew = brk_sv[:, 0] ** 2 / np.sum(brk_sv ** 2, axis=1)
    print()
    print(f"=== Concentration of leading SV (sv₁² / Σ sv²) ===")
    print(f"  F71-sym:      min {sym_skew.min():.3f}, max {sym_skew.max():.3f}, "
          f"mean {sym_skew.mean():.3f}")
    print(f"  F71-breaking: min {brk_skew.min():.3f}, max {brk_skew.max():.3f}, "
          f"mean {brk_skew.mean():.3f}")

    pairs = 0
    sym_wins = 0
    for s in sym_skew:
        for b in brk_skew:
            pairs += 1
            if s > b:
                sym_wins += 1
    print(f"  P(F71-sym sv₁²-fraction > F71-breaking) ≈ {sym_wins / pairs:.3f}")

    # Save raw
    out = SCRIPT_DIR / "results" / "eq024_f71_optimality_svd.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {'n_samples': n_samples, 'seed': seed, 'N': N},
            'f71_sym': sym_data,
            'f71_breaking': brk_data,
        }, f, indent=1)
    print()
    print(f"Saved: {out}")
    print(f"Total runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
