"""Exact modular machinery for deciding a spectral palindrome, plus its guard.

The palindrome is spec(L) = {-lambda - s}. With p(x) = det(xI - L) and the
dimension 4^N even,

    prod_i (x + s + lambda_i) = (-1)^n p(-x - s) = p(-x - s),

so the spectrum is palindromic IFF p(x) == p(-x - s) identically. Every entry
of L is a Gaussian rational whenever gamma, J, delta and the field magnitudes
are rationals, so L maps into GF(p) once i is sent to a square root of -1, and
the identity is tested by comparing determinants. No eigensolver, no rounding.

WHAT THIS IS AND IS NOT. A 'not palindromic' verdict is a PROOF: a single
witness x with p(x) != p(-x-s) mod p rules the identity out over Q(i) too. A
'palindromic' verdict is NOT a decision, and saying so was the first version's
error. It has two escape channels:

  (a) Schwartz-Zippel. Two distinct degree-n polynomials agree at at most n of
      the p field elements, so one random x lets a difference hide with
      probability <= n/p.
  (b) THE MODULUS ITSELF. q = p(x) - p(-x-s) can be nonzero over Q(i) and
      vanish identically mod p. No number of evaluation points sees this.
      Minimal witness, checked in main() below: L = diag(1, p-1), s = 0, whose
      true spectrum {1, p-1} is not palindromic while mod p it reads {1, -1},
      which is.

Channel (b) is why this module works over SEVERAL independent primes. A
spurious 'palindromic' then needs every coefficient of q divisible by all of
them at once. The honest label for the result is a one-sided certificate:
breaks are proved, holds are certified to a stated bound.

Run this file directly to execute its guards.
"""

import random

import numpy as np

# Primes 1 mod 4, so that -1 is a square and i has an image in each.
PRIMES = (998244353, 1004535809, 469762049)


def sqrt_minus_one(p):
    """An element i with i*i == -1 mod p, found rather than assumed."""
    for g in range(2, 200):
        cand = pow(g, (p - 1) // 4, p)
        if (cand * cand) % p == p - 1:
            return cand
    raise ValueError(f'no square root of -1 mod {p}')


def inv_p(a, p):
    return pow(int(a) % p, p - 2, p)


def det_mod_np_p(mat, p):
    """Determinant mod p by vectorized Gaussian elimination.

    No accumulation anywhere: the inner operation is one product subtracted
    from one entry, so with p < 2^30 every intermediate stays inside int64.
    Guarded against the silent-overflow failure mode by main() below, which
    compares this against an arbitrary-precision object-dtype replica.
    """
    a = np.asarray(mat, dtype=np.int64) % p
    n = a.shape[0]
    det = 1
    for c in range(n):
        nz = np.nonzero(a[c:, c])[0]
        if nz.size == 0:
            return 0
        piv = c + int(nz[0])
        if piv != c:
            a[[c, piv]] = a[[piv, c]]
            det = (-det) % p
        det = (det * int(a[c, c])) % p
        ic = inv_p(int(a[c, c]), p)
        rows = np.nonzero(a[c + 1:, c])[0]
        if rows.size:
            rows = rows + c + 1
            f = (a[rows, c] * ic) % p
            a[rows, c:] = (a[rows, c:] - f[:, None] * a[c, c:]) % p
    return det % p


def det_mod_object(mat, p):
    """The same determinant in Python ints. Slow, cannot overflow, the guard."""
    a = [[int(v) % p for v in row] for row in mat]
    n = len(a)
    det = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if a[r][c]), None)
        if piv is None:
            return 0
        if piv != c:
            a[c], a[piv] = a[piv], a[c]
            det = (-det) % p
        det = (det * a[c][c]) % p
        ic = inv_p(a[c][c], p)
        for r in range(c + 1, n):
            if a[r][c]:
                f = (a[r][c] * ic) % p
                for k in range(c, n):
                    a[r][k] = (a[r][k] - f * a[c][k]) % p
    return det % p


def _main():
    rng = random.Random(20260803)
    print('guard 1: sqrt(-1) exists in every prime used')
    for p in PRIMES:
        i = sqrt_minus_one(p)
        assert (i * i) % p == p - 1
        assert p % 4 == 1
        print(f'  p={p:11d}  i={i:11d}  i^2 = -1 mod p: OK')

    print('guard 2: the vectorized determinant equals the overflow-proof one')
    for trial in range(4):
        p = PRIMES[trial % len(PRIMES)]
        n = 24
        mat = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
        fast = det_mod_np_p(np.array(mat, dtype=np.int64), p)
        slow = det_mod_object(mat, p)
        assert fast == slow, (trial, fast, slow)
        print(f'  trial {trial}: p={p:11d}  det={fast:11d}  match: OK')

    print('guard 3: a single modulus really can certify a false palindrome')
    p = PRIMES[0]
    L = np.array([[1, 0], [0, p - 1]], dtype=np.int64)
    eye = np.eye(2, dtype=np.int64)
    agree = all(det_mod_np_p((x * eye - L) % p, p)
                == det_mod_np_p(((-x) % p * eye - L) % p, p)
                for x in (rng.randrange(p) for _ in range(8)))
    print(f'  L = diag(1, p-1), s = 0: true spectrum {{1, {p - 1}}} over Z is '
          f'NOT palindromic,')
    print(f'  yet the single-prime test reports palindromic at every point: '
          f'{agree}')
    assert agree, 'the documented escape channel did not reproduce'
    print('  -> this is why the sweep runs several primes, and why a HOLDS')
    print('     verdict is a certificate rather than a decision.')
    print()
    print('all guards passed')


if __name__ == '__main__':
    _main()
