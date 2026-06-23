"""Polarity probe #4: k >= 3 body Hamiltonians.

Tests whether the +/-1/2 polarity balance survives k-body Hamiltonians beyond
2-body bilinears. The polarity-coordinates wave found bit-exact balance for all
6 tested bilinear (k=2) families under various dissipators. This probe extends
to k=3 chain Hamiltonians via existing pi_decompose_M API which already handles
k-body via _build_kbody_chain.

If balance still holds at k=3: Hermitian-conjugate symmetry explanation is
unaffected by k (which makes structural sense).
If balance breaks at k=3: the explanation is finer than pure
Hermitian-conjugate symmetry; possibly tied to 2-body specific structure.
"""

import sys
sys.path.insert(0, 'simulations')
import framework as fw

N = 4
gamma_z = 0.05

# k=3 chain Hamiltonian families. Each tuple is (a, b, c) = a_l b_(l+1) c_(l+2)
# at every chain triplet (sliding window).
kbody_cases = [
    ("k=3: XYZ pure",                 [('X', 'Y', 'Z')]),
    ("k=3: ZYX (Hermitian conjugate)", [('Z', 'Y', 'X')]),
    ("k=3: XYZ + ZYX (Hermitian pair)", [('X', 'Y', 'Z'), ('Z', 'Y', 'X')]),
    ("k=3: YYY pure",                 [('Y', 'Y', 'Y')]),
    ("k=3: XXX + YYY + ZZZ",          [('X', 'X', 'X'), ('Y', 'Y', 'Y'), ('Z', 'Z', 'Z')]),
    ("k=3: XYI + IXY",                [('X', 'Y', 'I'), ('I', 'X', 'Y')]),
]

chain = fw.ChainSystem(N=N, gamma_0=gamma_z)

CASE_W = 50
TABLE_W = 110

print("=" * TABLE_W)
print(f"Polarity probe #4: k=3 chain Hamiltonians at N={N}, gz={gamma_z} (pure Z-deph)")
print("=" * TABLE_W)
print(f"{'Case':<{CASE_W}}  {'||M||^2':>12}  {'0%':>7}  {'+1/2%':>7}  {'-1/2%':>7}  {'asym':>12}")
print("-" * TABLE_W)

for label, terms in kbody_cases:
    try:
        # strict=False because F81 identity M_anti = L_{H_odd} was proven for
        # 2-body H only; at k>=3 the identity has a non-zero residual which is
        # itself a discovery channel.
        result = fw.polarity_coordinates(chain, terms, gamma_z=gamma_z, strict=False)
        ns_M = result['norm_sq']['M']
        if ns_M > 1e-12:
            p0 = 100.0 * result['norm_sq']['M_zero'] / ns_M
            pp = 100.0 * result['norm_sq']['M_plus_half'] / ns_M
            pm = 100.0 * result['norm_sq']['M_minus_half'] / ns_M
            print(f"{label:<{CASE_W}}  {ns_M:>12.4f}  {p0:>6.2f}%  {pp:>6.2f}%  {pm:>6.2f}%  {result['asymmetry']:+.3e}")
        else:
            print(f"{label:<{CASE_W}}  {ns_M:>12.2e}    ---      ---      ---  {result['asymmetry']:+.3e}")
    except Exception as e:
        # Encode error message to ASCII so Windows cp1252 console can print it.
        msg = repr(str(e)).encode('ascii', 'backslashreplace').decode('ascii')
        print(f"{label:<{CASE_W}}  ERROR: {type(e).__name__}: {msg}")

print()
print("=" * TABLE_W)
print("Reading:")
print("=" * TABLE_W)
print("  asym = 0 (bit-exact) -> +/-1/2 balance preserved at k=3 (consistent with Hermitian-conjugate explanation)")
print("  asym != 0            -> balance broken at k=3 -> new structural axis beyond Hermitian symmetry")
