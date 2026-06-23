"""F112 non-Hermitian extension: N=5 enumeration runner.

Calls enumerate_F_on_basis(5) from f112_open_identity_basis_enum.py.

Estimated cost: ~16 GB working memory (1024 pre-computed L_alpha,-i matrices
of size 1024x1024 complex = 16 MB each); ~30 min on multi-core OpenBLAS
(524,800 distinct upper-triangular pairs, each requiring one Frobenius
inner product).

If F = 0 bit-exact across all 524,800 pairs: F112 non-Hermitian extension
extends from Tier1Derived N <= 4 to Tier1Derived N <= 5 via basis-spanning
argument (per LindbladBitBPiBalance.NonHermitianExtension docstring).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from f112_open_identity_basis_enum import enumerate_F_on_basis

print("F112 non-Hermitian extension: N=5 enumeration")
print("=" * 80)
print("Welle 10a (2026-05-26): extending N <= 4 Tier1Derived to N = 5.")
print("Expected: ~16 GB working memory, several minutes runtime.")
print()

t0 = time.time()
max_im, mean_im, count_nonzero, total = enumerate_F_on_basis(5, tol=1e-10)
elapsed = time.time() - t0

print()
print("=" * 80)
print(f"WELLE 10a RESULT (elapsed {elapsed:.1f} s = {elapsed/60:.1f} min):")
print(f"  max |Im|: {max_im:.4e}")
print(f"  mean |Im|: {mean_im:.4e}")
print(f"  pairs with |Im| > 1e-10: {count_nonzero} / {total}")
if count_nonzero == 0:
    print(f"  *** PASS: F = 0 bit-exact across all {total} N=5 pairs ***")
    print("  *** F112 non-Hermitian Tier1Derived extends to N <= 5 ***")
else:
    print(f"  *** FAIL: {count_nonzero} non-zero entries at N=5 ***")
