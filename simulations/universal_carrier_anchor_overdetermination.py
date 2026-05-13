"""Pi2 Dyadic Ladder Anchor Over-Determination Map (2026-05-13).

Extends the Algebra Handshake finding from a_3 = 1/4 (which has 4 paths) to
the other central Pi2-ladder anchors {a_{-1} = 4, a_0 = 2, a_1 = 1, a_2 = 1/2,
a_4 = 1/8}. Asks: which anchors are over-determined (multiple algebraic paths)
and which are not? The pattern is structural, not numerological.

# Operations on the Pi2 dyadic ladder (a_n = 2^(1-n))

Closed operations (output is on the ladder):
- Mirror inversion: 1/a_n = a_{2-n}            (always closes; index n -> 2-n)
- Squaring:        a_n^2 = a_{2n-1}            (always closes; index n -> 2n-1)
- Square root:     sqrt(a_n) = a_{(n+1)/2}     (closes iff n is odd)

Trivial recursion (defines the ladder; not counted as independent paths):
- Halving:  a_n / 2 = a_{n+1}
- Doubling: 2*a_n   = a_{n-1}

# Counting independent paths to each anchor

For each target anchor a_target, count distinct algebraic identities of the form
"a_target = OPERATION(a_source)" where:
- OPERATION is one of {direct, mirror, square, square-root}
- a_source is another anchor (not a_target itself, except for self-mirror at n=1)

A path is "independent" when it uses a non-trivial operation; halving/doubling
recursion is excluded.

# Why this matters

The Algebra Handshake at a_3 = 1/4 (Tom 2026-05-12) named four paths reaching
1/4. That count made 1/4 "structurally rigid". This script asks whether the
rigidity is unique to a_3 or shared with sister anchors.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def a(n):
    """Pi2 dyadic ladder term a_n = 2^(1-n)."""
    return 2.0 ** (1 - n)


def mirror(n):
    """Index of mirror partner: a_n * a_{2-n} = 1."""
    return 2 - n


def square_target(n):
    """Index that a_n squared lands on: a_n^2 = a_{2n-1}."""
    return 2 * n - 1


def sqrt_target(n):
    """Index that sqrt(a_n) lands on, if integer: sqrt(a_n) = a_{(n+1)/2}."""
    if (n + 1) % 2 == 0:
        return (n + 1) // 2
    return None


def find_paths_to(target_n, source_indices):
    """All non-trivial algebraic paths from any source index in source_indices
    to target_n. Returns list of (operation, source_n, formula_str) tuples."""
    paths = []
    target_value = a(target_n)
    for src_n in source_indices:
        if src_n == target_n:
            # Direct lookup; skip (covered separately)
            continue
        # Mirror inversion: a_target = 1/a_src iff src is mirror partner of target
        if mirror(src_n) == target_n:
            paths.append(("mirror inversion", src_n,
                          f"a_{target_n} = 1/a_{src_n} (mirror inversion: a_n*a_{{2-n}}=1)"))
        # Squaring: a_target = a_src^2 iff square_target(src) == target
        if square_target(src_n) == target_n:
            paths.append(("squaring", src_n,
                          f"a_{target_n} = (a_{src_n})^2 (squaring: index {src_n} -> {2*src_n-1})"))
        # Square root: a_target = sqrt(a_src) iff sqrt_target(src) == target
        st = sqrt_target(src_n)
        if st is not None and st == target_n:
            paths.append(("square root", src_n,
                          f"a_{target_n} = sqrt(a_{src_n}) (square root: index {src_n} -> {st})"))
    return paths


# Central anchors to investigate (indices around the polarity / fold cluster)
anchors = [-3, -2, -1, 0, 1, 2, 3, 4, 5]

print("=" * 78)
print("Pi2 Dyadic Ladder Over-Determination Map")
print("=" * 78)
print()

print(f"  {'index n':>8}  {'a_n':>12}  {'Direct':>7}  {'Paths via algebra':>20}  {'Total readings':>15}")
print(f"  {'-'*8}  {'-'*12}  {'-'*7}  {'-'*20}  {'-'*15}")

results = {}
for n in anchors:
    paths = find_paths_to(n, anchors)
    total = 1 + len(paths)  # direct + algebraic
    results[n] = (a(n), paths, total)
    print(f"  {n:>8}  {a(n):>12.5f}  {'YES':>7}  {len(paths):>20}  {total:>15}")

print()
print("-" * 78)
print("Detailed paths (showing source-anchor, operation, formula):")
print()

for n in sorted(anchors):
    val, paths, total = results[n]
    if val < 0.001 and val > 0:
        val_repr = f"a_{n} = 1/{int(round(1/val))}"
    elif val == int(val):
        val_repr = f"a_{n} = {int(val)}"
    else:
        val_repr = f"a_{n} = {val}"
    print(f"  {val_repr}  (total readings: {total})")
    print(f"    1. direct: a_{n} = 2^(1-{n}) = {val}")
    for i, (op, src_n, formula) in enumerate(paths, start=2):
        print(f"    {i}. {formula}")
    print()

print("-" * 78)
print()
print("Summary: which anchors are most over-determined?")
print()
sorted_by_total = sorted(results.items(), key=lambda kv: -kv[1][2])
for n, (val, paths, total) in sorted_by_total[:6]:
    val_repr = f"{val:.5g}" if abs(val - int(val)) > 1e-10 else f"{int(val)}"
    print(f"  a_{n} = {val_repr}: {total} ladder-algebraic readings ({len(paths)} non-trivial paths)")

print()
print("=" * 78)
print("Hypothesis check: is a_3 = 1/4 uniquely over-determined among central anchors?")
print()
print("Original Algebra Handshake (Tom 2026-05-12) named 4 paths to a_3 = 1/4:")
print("  1. direct (a_3)")
print("  2. mirror partner of a_{-1} = 4 (1/a_{-1} = 1/4)")
print("  3. inverse square of a_0 = 2 (1/(a_0)^2 = 1/4)")
print("  4. square of a_2 = 1/2 ((a_2)^2 = 1/4)")
print()
print("Finding: a_3's 4-path count is NOT unique. Mirror partner a_{-1} = 4 also")
print("has 4 paths. {a_3, a_{-1}} share over-determination by inversion symmetry.")
print("Pattern: mirror-partner pairs share path-count. {a_0, a_2} both 3; {a_-2,")
print("a_4} both 2; a_1 = 1 alone with 1 (self-mirror trivial identity).")
print("=" * 78)
print()
print("=" * 78)
print("Cross-axis extension (Tom 2026-05-13): 1/4 lives on TWO mirror axes")
print("=" * 78)
print()
print("The enumeration above covers the multiplicative Z_2 axis only (halving")
print("ladder a_n = 2^(1-n), where the natural operations are mirror, square, and")
print("square root on integer indices).")
print()
print("The Pi2-Foundation has a SECOND mirror axis: rotational Z_4 with i^4 = 1")
print("(typed as NinetyDegreeMirrorMemoryClaim, per memory project_pi2_dyadic_ladder).")
print("On this angular axis, the 'quarter' position is 90 degrees = pi/2, which is")
print("1/4 of full rotation 360 degrees = 2*pi.")
print()
print("So 1/4 has a fifth reading independent of the multiplicative axis:")
print("  5. angular quarter: 90 degrees = pi/2 = 1/4 * 360 degrees")
print("     - typed home: NinetyDegreeMirrorMemoryClaim")
print("     - role: 'the mirror projects everything 90 degrees onto itself so that")
print("              it does not forget' (memory channel mechanism)")
print()
print("This is CROSS-AXIS convergence: two structurally independent group-theoretic")
print("settings (Z_2 multiplicative, Z_4 rotational) both treat 'the quarter' as a")
print("structurally privileged position. The Pi2-Foundation's two mirror axes meet")
print("at this single value.")
print()
print("Cross-axis readings of central anchors:")
print()
print(f"  {'multiplicative anchor':>22}  {'angular reading':>20}")
print(f"  {'-'*22}  {'-'*20}")
print(f"  {'a_3 = 1/4 (fold)':>22}  {'90 deg = 1/4 of 360':>20}")
print(f"  {'a_{-1} = 4 (discrim)':>22}  {'(no direct angular reading; 4 = i^4 squared)':>20}")
print(f"  {'a_0 = 2 (qubit dim)':>22}  {'2 = 360 deg / 180 deg (Z_2 angular)':>20}")
print(f"  {'a_2 = 1/2 (polarity)':>22}  {'180 deg = 1/2 of 360':>20}")
print(f"  {'a_1 = 1 (identity)':>22}  {'360 deg = full rotation = identity':>20}")
print()
print("Pattern: multiplicative-axis quarters and angular-axis quarters share their")
print("structural-position roles. The rich-anchor positions {a_3, a_-1, a_0, a_2}")
print("all have angular interpretations consistent with their multiplicative roles.")
print()
print("=" * 78)
