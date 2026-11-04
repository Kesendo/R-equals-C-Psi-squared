"""Even-parity sector of the dephased XY chain: is it point-symmetric about its mean?

A Liouvillian quadratic in the Jordan-Wigner fermions (free-fermion, Nambu) carries on its even-parity operator sector a
constant plus the even-cardinality subset sums of 2N rapidities. An even subset of 2N
elements has an even complement, so that sector is point-symmetric about its mean. This
script builds the even sector (popcount(ket) + popcount(bra) even) of the open XY chain,
H = (J/2) sum (X_l X_{l+1} + Y_l Y_{l+1}), under uniform Z dephasing D = gamma (Z rho Z - rho),
and prints the worst mismatch between the sector and its reflection through its mean.
A nonzero mismatch excludes a generator quadratic in the Jordan-Wigner fermions at that N.

Recorded in docs/CAUGHT_ERRORS.md, entry of 2026-11-04 (the Nambu reduction).
"""
import numpy as np


def even_sector(N, J=1.0, gamma=0.3):
    d = 2 ** N

    def bits(x):
        return [(x >> (N - 1 - i)) & 1 for i in range(N)]

    pc = [sum(bits(x)) for x in range(d)]
    H = np.zeros((d, d))
    for x in range(d):
        b = bits(x)
        for i in range(N - 1):
            if b[i] != b[i + 1]:
                H[x ^ (1 << (N - 1 - i)) ^ (1 << (N - 2 - i)), x] += J
    cells = [(a, b) for a in range(d) for b in range(d) if (pc[a] + pc[b]) % 2 == 0]
    pos = {t: i for i, t in enumerate(cells)}
    M = np.zeros((len(cells), len(cells)), complex)
    for (a, b), i in pos.items():
        M[i, i] -= 2 * gamma * bin(a ^ b).count("1")
        for a2 in range(d):
            if H[a2, a]:
                M[pos[(a2, b)], i] += -1j * H[a2, a]
        for b2 in range(d):
            if H[b, b2]:
                M[pos[(a, b2)], i] += 1j * H[b, b2]
    return np.linalg.eigvals(M)


if __name__ == "__main__":
    for N in (2, 3, 4, 5):
        ev = even_sector(N)
        mean = ev.mean()
        worst = max(np.min(np.abs(ev - (2 * mean - z))) for z in ev)
        print(f"N={N}: even-sector dim {len(ev)}, mean {mean.real:+.6f}{mean.imag:+.6f}i, "
              f"worst point-reflection mismatch {worst:.3e}")
