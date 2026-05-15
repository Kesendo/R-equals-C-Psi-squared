"""Extend F89 path-k closed-form tabulation from k=24 to k=25+ via the Chebyshev pipeline.

Reuses extract_path_polynomial from f89_pathk_symbolic_derivation.py. Runs for
k = 25, 26, ... until either sympy timing or a manual cutoff. Emits coefficient
arrays in the C# tabulation format ready to paste into F89UnifiedFaClosedFormClaim.

Run: python simulations/_f89_pathk_extend_k25_plus.py
"""
import sys
import time

sys.path.insert(0, "simulations")
from f89_pathk_symbolic_derivation import extract_path_polynomial, predicted_D

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def fmt_csharp_array(coefs):
    """Format coefficient array low-to-high as C# initializer with .0 doubles."""
    parts = [f"{int(c)}.0" for c in coefs]
    return "{ " + ", ".join(parts) + " }"


print("# F89 path-k tabulation extension k=25+")
print(f"# {'k':>3} {'time(s)':>8} {'D_ext':>12} {'D_pred':>12} {'match':>6}  P_k(y)")

MIN_K = 33   # first k not yet in F89UnifiedFaClosedFormClaim tabulation
MAX_K = 50   # adjust based on tolerance; sympy time grows ~O(k^2)
TIME_BUDGET_SECONDS = 600  # per-k upper bound; abort if exceeded
total_t0 = time.time()

new_polys = {}
for k in range(MIN_K, MAX_K + 1):
    t0 = time.time()
    try:
        coefs, D_ext, orbit = extract_path_polynomial(k)
        t1 = time.time()
        elapsed = t1 - t0
        D_pred = predicted_D(k)
        ok = (D_ext == D_pred)
        new_polys[k] = (coefs, D_ext)
        # Format polynomial display
        terms = []
        for i, c in enumerate(coefs):
            if c == 0:
                continue
            if i == 0:
                terms.append(f"{c}")
            elif i == 1:
                terms.append(f"{c}*y")
            else:
                terms.append(f"{c}*y^{i}")
        poly_str = " + ".join(reversed(terms))
        print(f"  {k:>3} {elapsed:>8.2f} {D_ext:>12} {D_pred:>12} {('OK' if ok else 'FAIL'):>6}  P_{k}(y) = {poly_str}")
        if elapsed > TIME_BUDGET_SECONDS:
            print(f"# WARNING: k={k} took {elapsed:.1f}s, near budget; stopping")
            break
    except Exception as e:
        print(f"  k={k}: error {type(e).__name__}: {e}")
        break

total_elapsed = time.time() - total_t0
print(f"\n# Total time: {total_elapsed:.1f}s for k=25..{max(new_polys.keys()) if new_polys else '24'}")

# Emit C# tabulation strings for F89UnifiedFaClosedFormClaim
print("\n# === C# coefficient arrays (paste into F89UnifiedFaClosedFormClaim.cs) ===\n")
for k in sorted(new_polys.keys()):
    coefs, D = new_polys[k]
    print(f"    private static readonly double[] _path{k}Coefs = {fmt_csharp_array(coefs)};")

print("\n# === Switch entries for PathPolynomial ===\n")
for k in sorted(new_polys.keys()):
    coefs, D = new_polys[k]
    print(f"            {k} => (_path{k}Coefs, {D}),")
