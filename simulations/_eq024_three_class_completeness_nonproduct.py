#!/usr/bin/env python3
"""EQ-024 sub-question 1b: three-class completeness in NON-PRODUCT F71-sym states.

Companion to _eq024_three_class_completeness.py. The product-state version
(2026-04-28, 100 samples, all C ≥ 7.54 bits, no fourth-class candidate)
established empirical support for completeness within F71-symmetric *product*
states (6-parameter Bloch family).

Non-product F71-symmetric pure states span a much larger subspace at N=5: the
+1 eigenspace of the chain-mirror operator R (|a₀a₁a₂a₃a₄⟩ → |a₄a₃a₂a₁a₀⟩),
which has dimension 20 (8 self-mirror computational basis states + 12 mirror-
pair symmetric combinations). This gives 38 real parameters modulo norm and
phase — much more room for a fourth blind mechanism to hide than the product
slice.

Sampling strategy: Haar-random pure state in C^32, project onto F71-symmetric
+1-eigenspace via P_sym = (I + R)/2, renormalise. This gives a Haar-uniform
random sample in the 20-dim symmetric subspace.

Same class-membership tests and J-Jacobian pipeline as the product version.
Reuses helper functions from _eq024_three_class_completeness.py.
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

# Import shared infrastructure from the product-version script
from _eq024_three_class_completeness import (
    N, d,
    classify, is_class_1_or_2, is_class_3,
    compute_jacobian_J, waterfilling_capacity,
    CLASS_TOL, CAPACITY_BLIND,
    GAMMA_0, J_REF, DJ, T_POINTS, SPREAD, SIGMA_NOISE,
)


def build_R_mirror(N):
    """Chain-mirror operator R: |b₀…b_{N-1}⟩ → |b_{N-1}…b₀⟩ on C^{2^N}."""
    dim = 2 ** N
    R = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        rev = 0
        x = i
        for _ in range(N):
            rev = (rev << 1) | (x & 1)
            x >>= 1
        R[rev, i] = 1.0
    return R


R = build_R_mirror(N)
P_SYM = (np.eye(d, dtype=complex) + R) / 2.0


def random_f71_symmetric(rng):
    """Haar-random pure state in the F71-symmetric (+1 R-eigenspace) subspace.

    Method: draw Gaussian random vector in C^d, project onto +1 R-eigenspace,
    normalise. Equivalent to Haar-uniform on the 20-dim symmetric subspace.
    """
    g = rng.standard_normal(d) + 1j * rng.standard_normal(d)
    psi = P_SYM @ g
    psi /= np.linalg.norm(psi)
    # Verify R-symmetry numerically
    err = np.linalg.norm(R @ psi - psi)
    assert err < 1e-10, f"R-symmetry violated: ‖Rψ − ψ‖ = {err}"
    return psi


def main():
    n_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = np.random.default_rng(seed)

    print(f"EQ-024 Three-Class Completeness — NON-PRODUCT F71-symmetric states")
    print(f"  N={N}, Heisenberg, γ₀={GAMMA_0}, J={J_REF}, dJ={DJ}")
    print(f"  Sample space: 20-dim +1 eigenspace of chain-mirror R")
    print(f"  Random samples: {n_samples}, seed={seed}")
    print(f"  Tolerance: class membership {CLASS_TOL:.0e}, "
          f"J-blind threshold {CAPACITY_BLIND} bits")
    print()

    in_class = {'1or2': 0, '3': 0}
    outside_capacities = []
    blind_outside = []  # candidate fourth-class states

    t0 = time.time()
    for i in range(n_samples):
        psi = random_f71_symmetric(rng)

        verdict, label = classify(psi)
        if verdict == 'in-class':
            in_class[label] += 1
            continue

        rho0 = np.outer(psi, psi.conj())
        rho0 = (rho0 + rho0.conj().T) / 2.0
        rho0 /= np.trace(rho0).real
        A = compute_jacobian_J(rho0)
        sv = np.linalg.svd(A, compute_uv=False)
        C = waterfilling_capacity(sv)
        outside_capacities.append(C)

        if C < CAPACITY_BLIND:
            blind_outside.append({
                'sample_idx': i,
                'capacity_bits': C,
                'sv_max': float(sv[0]),
                'sv': [float(s) for s in sv],
                'psi_real': psi.real.tolist(),
                'psi_imag': psi.imag.tolist(),
            })

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  [{i + 1}/{n_samples}] elapsed {elapsed:.1f}s, "
                  f"in-class {sum(in_class.values())}, "
                  f"outside {len(outside_capacities)}, "
                  f"blind-outside {len(blind_outside)}")

    print()
    print("=== Summary ===")
    print(f"Total samples: {n_samples}")
    print(f"In-class breakdown: Class 1or2: {in_class['1or2']}, Class 3: {in_class['3']}")
    print(f"Outside-class samples: {len(outside_capacities)}")
    if outside_capacities:
        arr = np.array(outside_capacities)
        print(f"  capacity range: {arr.min():.4f} – {arr.max():.4f} bits")
        print(f"  capacity mean: {arr.mean():.4f}, std: {arr.std():.4f}")
        print(f"  capacity median: {np.median(arr):.4f}")
    print()
    print(f"Blind outside (C < {CAPACITY_BLIND} bits, fourth-class candidates): "
          f"{len(blind_outside)}")
    if blind_outside:
        print("  *** WARNING: fourth-class candidate(s) found ***")
        for entry in blind_outside[:5]:
            print(f"    sample {entry['sample_idx']}  C={entry['capacity_bits']:.6f}  "
                  f"sv_max={entry['sv_max']:.3e}")
    else:
        print("  None found — strong empirical support for three-class completeness")
        print("  within F71-symmetric pure states (full 20-dim symmetric subspace).")

    print()
    print(f"Total runtime: {time.time() - t0:.1f}s")

    out = SCRIPT_DIR / "results" / "eq024_three_class_completeness_nonproduct.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {
                'N': N, 'gamma_0': GAMMA_0, 'J_ref': J_REF, 'dJ': DJ,
                't_points': T_POINTS, 'spread': SPREAD, 'sigma_noise': SIGMA_NOISE,
                'n_samples': n_samples, 'seed': seed,
                'sample_space': '20-dim F71-symmetric subspace (Haar-random)',
                'class_tol': CLASS_TOL, 'capacity_blind_threshold': CAPACITY_BLIND,
            },
            'in_class_counts': in_class,
            'outside_capacities': outside_capacities,
            'blind_outside_candidates': blind_outside,
        }, f, indent=1)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
