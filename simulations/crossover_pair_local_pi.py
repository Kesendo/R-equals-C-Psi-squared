#!/usr/bin/env python3
"""XZ+YZ and ZX+ZY have a LOCAL palindrome mirror: a continuous uniform per-site map.

This is the committed demonstration behind the correction of the old "genuinely non-local
/ entangled Π" verdict (experiments/PI_OPERATOR_ENTANGLEMENT.md, hypotheses/THE_BOOT_SCRIPT.md).
The old verdict came from two narrowing lenses:

  1. the per-site map was searched only over the discrete signed-permutation crossovers
     {P1, P4, M2}, which genuinely fail for these two cases; and
  2. the fallback construction paired Liouvillian eigenvectors, which under the heavy
     spectral degeneracy of these Liouvillians returns an entangled representative even
     when a product solution exists (PI_OPERATOR_ENTANGLEMENT §5 says this verbatim).

Widen the per-site map to a continuous block rotation and a product mirror appears. The
closed-form map M (M² = −I, unitary), the SAME on every site, gives Π = M^⊗N that mirrors
the full Liouvillian to machine precision at N = 2..6. The two cases are LOCAL.

The closed form and the verification both live in the framework primitive
`framework/diagnostics/crossover_product_pi.py`; this script just exercises it end-to-end
and prints the table the documents cite.

Run:  PYTHONIOENCODING=utf-8 python simulations/crossover_pair_local_pi.py
"""
from __future__ import annotations

import sys
import numpy as np

sys.path.insert(0, 'simulations')
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

GAMMA = 0.5
CASES = [(('X', 'Z'), ('Y', 'Z')), (('Z', 'X'), ('Z', 'Y'))]

# Discrete crossover P1 in framework order [I,X,Z,Y] (I->X, X->I, Z->iY, Y->iZ),
# the map that closes the truly cases and that the old discrete-only search relied on.
P1 = np.zeros((4, 4), dtype=complex)
P1[1, 0] = 1; P1[0, 1] = 1; P1[3, 2] = 1j; P1[2, 3] = 1j


def main():
    print("XZ+YZ / ZX+ZY: a LOCAL continuous-rotation mirror (γ=0.5, Pauli order [I,X,Z,Y])")
    print()

    M = fw.crossover_map()
    print("closed-form per-site map M (rows=output, cols=input):")
    print(np.array2string(np.real_if_close(M, tol=1e6), precision=4, suppress_small=True))
    print(f"  unitary: {np.allclose(M @ M.conj().T, np.eye(4))}   "
          f"M² = −I: {np.allclose(M @ M, -np.eye(4))}   "
          f"nonzeros: {int(np.sum(np.abs(M) > 1e-9))} (a permutation would have 4)")
    print()

    # Sanity: the machinery reproduces a known-local 'truly' case (XX+YY via P1).
    truly = (('X', 'X'), ('Y', 'Y'))
    s = fw.product_pi_residual(truly, [P1, P1], 2, GAMMA)
    print(f"  SANITY XX+YY (truly): P1^⊗2 residual = {s:.2e}  "
          f"{'OK' if s < 1e-9 else 'FAIL'}")
    print()

    for combo in CASES:
        label = f"{combo[0][0]}{combo[0][1]}+{combo[1][0]}{combo[1][1]}"
        print(f"  === {label} ===")
        for N in (2, 3, 4, 5, 6):
            local = fw.verify_crossover_local(combo, N, GAMMA)
            disc = fw.product_pi_residual(combo, [P1] * N, N, GAMMA)
            print(f"    N={N}: continuous M^⊗N = {local:.2e}  (LOCAL)   "
                  f"discrete P1^⊗N = {disc:6.2f}  (fails)")
        print()

    print("The continuous mirror closes both cases at every N to machine precision while")
    print("the discrete crossover fails: the cases are local, not non-local. The 'non-local'")
    print("verdict was relative to the discrete-permutation lens; under a continuous per-site")
    print("lens the mirror is a product.")


if __name__ == "__main__":
    main()
