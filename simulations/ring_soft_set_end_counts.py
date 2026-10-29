"""The soft set of t1*XZX + t2*XZY + t3*YZX on the periodic ring, read through F158's exact end counts.

The question: on a ring of N sites (windows w = 0..N-1, sites mod N), Z dephasing on every site, which real
triples pair the Liouvillian spectrum about -sigma? docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md §8 ("the soft set
on rings") states: every triple at N = 4; t1 = 0 or t2 = t3 at 4 | N >= 8; t1 = 0 alone at 4 ∤ N.

This script is the second route beside the hand argument: for every integer triple of a grid it computes
dim ker L (near end) and dim ker(L + 2 sigma) (far end) exactly over Q(i) through
framework.diagnostics.f158_end_count, and checks pairing (near == far) against the statement. The live C# twin
is GoldenRouterClaim's ring battery (compute/RCPsiSquared.Diagnostics/F87/GoldenRouterClaim.cs).

It also reads both ends at N = 4, where the proof states six at generic points, nine on the line t2 = t3 and
ten on the plane t1 = 0 (each at a generic point of it), and six again on t2 + t3 = 0.

Run:  python simulations/ring_soft_set_end_counts.py            # N = 3..7 on {-2..2}^3, N = 8 on {-1..2}^3
Exits nonzero on any mismatch.
"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from framework.diagnostics.f158_end_count import end_count  # noqa: E402


def ring_terms(n, t1, t2, t3):
    """The weighted template on the periodic ring as {Pauli string: coefficient}."""
    h = {}
    for w in range(n):
        for (a, b), t in ((('X', 'X'), t1), (('X', 'Y'), t2), (('Y', 'X'), t3)):
            if t == 0:
                continue
            s = ['I'] * n
            s[w], s[(w + 1) % n], s[(w + 2) % n] = a, 'Z', b
            key = ''.join(s)
            h[key] = h.get(key, 0) + t
    return {k: v for k, v in h.items() if v != 0}


def pairs_by_statement(n, t1, t2, t3):
    if n == 4 or t1 == 0:
        return True
    return t2 == t3 and n % 4 == 0


def sweep(n, values):
    jumps = ['I' * l + 'Z' + 'I' * (n - l - 1) for l in range(n)]
    mismatches = []
    for t in itertools.product(values, repeat=3):
        if t == (0, 0, 0):
            continue
        ec = end_count(n, ring_terms(n, *t), jumps)
        if ec.pairs != pairs_by_statement(n, *t):
            mismatches.append((t, ec.near, ec.far))
    return mismatches


# (t1, t2, t3) -> the dimension both ends hold on the ring N = 4
N4_ENDS = {(1, 2, 3): 6, (2, -1, 5): 6, (1, 1, 1): 9, (5, 2, 2): 9, (0, 2, 7): 10, (1, 2, -2): 6}


def main():
    plan = [(n, range(-2, 3)) for n in range(3, 8)] + [(8, range(-1, 3))]
    failed = False
    jumps4 = ['ZIII', 'IZII', 'IIZI', 'IIIZ']
    for t, dim in N4_ENDS.items():
        ec = end_count(4, ring_terms(4, *t), jumps4)
        ok = ec.near == ec.far == dim
        print(f"ring N=4, {t}: near {ec.near}, far {ec.far}, expected {dim}: {'ok' if ok else 'MISMATCH'}", flush=True)
        failed |= not ok
    for n, values in plan:
        bad = sweep(n, values)
        print(f"ring N={n}, grid {{{values.start}..{values.stop - 1}}}^3: mismatches {len(bad)} {bad[:5]}", flush=True)
        failed |= bool(bad)
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
