#!/usr/bin/env python3
"""EQ-024 sub-question 2: F71-breaking receiver capacity distribution.

Direction 3 (April 23) swept only F71-symmetric receivers and saturated at
12.07 bits. The today's sub-question-1 closure showed F71-symmetric pure
states have C ∈ [7.54, 12.41] bits across 100 product + 100 non-product
samples. None blind.

This script tests F71-BREAKING receivers (no chain-mirror symmetry):

  Mode 'product': random product state ⊗_i |a_i⟩ with 10 independent
                  Bloch angles (5 sites × 2 angles), no symmetry constraint.
  Mode 'nonproduct': Haar-random pure state in full C^32 (62 real params).

For each sample, classify in/out of Classes 1-3 (Classes 1-3 are all
F71-symmetric, so F71-breaking samples are AUTOMATICALLY outside; we
nonetheless run the test to flag any edge case where Haar randomness
accidentally lands very close to a symmetric class). Compute J-Jacobian
and Shannon capacity.

Comparison to F71-symmetric distributions answers:
  - Does F71-symmetry constrain max capacity? (i.e., is there a
    F71-breaking receiver with C > 12.41 bits?)
  - Is the F71-breaking distribution shifted up, down, or comparable?
  - Are there F71-breaking states with C ≈ 0 (would be a structural
    surprise — F71-breaking should preserve all kinematic responsiveness)?
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
    N, d,
    classify, compute_jacobian_J, waterfilling_capacity,
    CLASS_TOL, CAPACITY_BLIND,
    GAMMA_0, J_REF, DJ, T_POINTS, SPREAD, SIGMA_NOISE,
    bloch_state,
)
from _eq024_three_class_completeness_nonproduct import (
    R as R_mirror,
)


def random_product_asymmetric(rng):
    """Random product state with 10 independent Bloch angles (no symmetry).

    Generically F71-breaking by sampling each site independently.
    """
    psi = np.array([1.0], dtype=complex)
    for _ in range(N):
        theta = rng.uniform(0, np.pi)
        phi = rng.uniform(0, 2 * np.pi)
        psi = np.kron(psi, bloch_state(theta, phi))
    return psi


def random_nonproduct_haar(rng):
    """Haar-random pure state in full C^32 (62 real params)."""
    g = rng.standard_normal(d) + 1j * rng.standard_normal(d)
    return g / np.linalg.norm(g)


def f71_asymmetry(psi):
    """Quantify how non-F71-symmetric a state is.

    Returns ‖ψ − R ψ‖ (zero if F71-symmetric, sqrt(2) if maximally asym).
    """
    return float(np.linalg.norm(psi - R_mirror @ psi))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'product'
    n_samples = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2

    if mode == 'product':
        sampler = random_product_asymmetric
        label = 'F71-breaking product (10 indep Bloch angles)'
    elif mode == 'nonproduct':
        sampler = random_nonproduct_haar
        label = 'F71-breaking non-product (Haar on C^32, 62 real params)'
    else:
        raise ValueError(f"unknown mode: {mode!r} (use 'product' or 'nonproduct')")

    rng = np.random.default_rng(seed)

    print(f"EQ-024 sub-question 2: F71-breaking receiver capacity")
    print(f"  Mode: {label}")
    print(f"  N={N}, Heisenberg, γ₀={GAMMA_0}, J={J_REF}, dJ={DJ}")
    print(f"  Random samples: {n_samples}, seed={seed}")
    print()

    in_class = {'1or2': 0, '3': 0}
    capacities = []
    asymmetries = []
    blind_outside = []

    t0 = time.time()
    for i in range(n_samples):
        psi = sampler(rng)
        asym = f71_asymmetry(psi)
        asymmetries.append(asym)

        verdict, lbl = classify(psi)
        if verdict == 'in-class':
            in_class[lbl] += 1
            capacities.append(0.0)  # in-class is by construction J-blind
            continue

        rho0 = np.outer(psi, psi.conj())
        rho0 = (rho0 + rho0.conj().T) / 2.0
        rho0 /= np.trace(rho0).real
        A = compute_jacobian_J(rho0)
        sv = np.linalg.svd(A, compute_uv=False)
        C = waterfilling_capacity(sv)
        capacities.append(C)

        if C < CAPACITY_BLIND:
            blind_outside.append({
                'sample_idx': i,
                'capacity_bits': C,
                'asymmetry': asym,
                'sv_max': float(sv[0]),
                'sv': [float(s) for s in sv],
            })

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  [{i + 1}/{n_samples}] elapsed {elapsed:.1f}s, "
                  f"in-class {sum(in_class.values())}, "
                  f"min asym {min(asymmetries):.3f}, "
                  f"blind-outside {len(blind_outside)}")

    print()
    print("=== Summary ===")
    arr_C = np.array(capacities)
    arr_A = np.array(asymmetries)
    print(f"Total samples: {n_samples}")
    print(f"In-class breakdown: Class 1or2: {in_class['1or2']}, Class 3: {in_class['3']}")
    print(f"F71-asymmetry ‖ψ − Rψ‖: min={arr_A.min():.4f}, "
          f"max={arr_A.max():.4f}, mean={arr_A.mean():.4f}")
    print(f"Capacity range: {arr_C.min():.4f} – {arr_C.max():.4f} bits")
    print(f"Capacity mean / median / std: {arr_C.mean():.4f} / "
          f"{np.median(arr_C):.4f} / {arr_C.std():.4f}")
    print()
    print(f"Blind candidates (C < {CAPACITY_BLIND}): {len(blind_outside)}")
    if blind_outside:
        print("  *** WARNING: F71-breaking blind state found ***")
        for entry in blind_outside[:5]:
            print(f"    sample {entry['sample_idx']}: C={entry['capacity_bits']:.6f}, "
                  f"asym={entry['asymmetry']:.3f}")
    print()
    print(f"Total runtime: {time.time() - t0:.1f}s")

    out = SCRIPT_DIR / "results" / f"eq024_f71_breaking_{mode}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {
                'mode': mode, 'sampler_label': label,
                'N': N, 'gamma_0': GAMMA_0, 'J_ref': J_REF, 'dJ': DJ,
                't_points': T_POINTS, 'spread': SPREAD, 'sigma_noise': SIGMA_NOISE,
                'n_samples': n_samples, 'seed': seed,
            },
            'in_class_counts': in_class,
            'capacities': capacities,
            'asymmetries': asymmetries,
            'blind_outside_candidates': blind_outside,
        }, f, indent=1)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
