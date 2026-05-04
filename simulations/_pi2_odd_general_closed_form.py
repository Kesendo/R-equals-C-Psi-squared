"""F88 Π²-odd/memory generalised closed form: numerical verifier.

Sweeps all popcount-coherence pair states (|p⟩+|q⟩)/√2 at N = 2..7 (any
popcount pair n_p, n_q, any HD compatible with popcount difference) and
verifies the closed-form prediction bit-exactly against the numerical
state-level reading (kernel projection + Pauli-basis Π²-projection).

Closed form, anchor structure (mirror, K-intermediate, generic, HD = N
Π²-classical), and connection to F60 documented in
docs/proofs/PROOF_F86_QPEAK.md §Structural inheritance from F88. C# port:
compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs.
"""
from __future__ import annotations

import sys
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def krawtchouk(n: int, s: int, N: int) -> int:
    k_lo = max(0, n - (N - s))
    k_hi = min(s, n)
    return sum((-1) ** k * comb(s, k) * comb(N - s, n - k) for k in range(k_lo, k_hi + 1))


def static_fraction_general(N: int, n_p: int, n_q: int) -> float:
    if n_p == n_q:
        return 1.0 / comb(N, n_p)
    return 0.25 / comb(N, n_p) + 0.25 / comb(N, n_q)


def alpha_krawtchouk_general(N: int, n_p: int, n_q: int) -> float:
    """Π²-odd-fraction-of-static via Krawtchouk closed form.

    Works for both inter-sector (n_p ≠ n_q) and intra-sector (n_p = n_q).
    For intra-sector A_s = B_s, so (A_s + B_s)² = 4·A_s² and the formula
    simplifies to Σ_{s odd} C(N,s) K_n(s;N)² / Σ_s C(N,s) K_n(s;N)².
    """
    Cnp = comb(N, n_p)
    Cnq = comb(N, n_q)
    odd_sum = 0.0
    total = 0.0
    for s in range(N + 1):
        A_s = krawtchouk(n_p, s, N) / Cnp
        B_s = krawtchouk(n_q, s, N) / Cnq
        term = comb(N, s) * (A_s + B_s) ** 2
        total += term
        if s & 1:
            odd_sum += term
    return odd_sum / total if total > 0 else 0.0


def pi2_odd_in_memory_general(N: int, n_p: int, n_q: int, hd: int) -> float:
    """Closed-form prediction Π²-odd/memory for popcount-coherence pair state."""
    if hd == N:
        return 0.0
    s = static_fraction_general(N, n_p, n_q)
    a = alpha_krawtchouk_general(N, n_p, n_q)
    return (0.5 - a * s) / (1.0 - s)


# ---------- numerical verification ----------

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, SX, SY, SZ]
BIT_B = [0, 0, 1, 1]


def kron_n(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def coherence_state(N, p_bits, q_bits, sign=+1):
    psi = np.zeros(2**N, dtype=complex)
    if p_bits == q_bits:
        # degenerate (would not be a coherence pair); skip
        return None
    psi[p_bits] = 1.0 / np.sqrt(2.0)
    psi[q_bits] = sign / np.sqrt(2.0)
    return np.outer(psi, psi.conj())


def kernel_projection_popcount(rho, N):
    """ρ_d0 = Σ_n (Tr(P_n · ρ) / Tr(P_n²)) · P_n with P_n the popcount-n sector
    projector. For Heisenberg + Z-dephasing, kernel of L = span{P_n}, so this is
    the canonical state-level kernel projector. Shared helper for the F88-Lens
    sweep + multi-state probe."""
    d = 2**N
    rho_d0 = np.zeros_like(rho)
    for n in range(N + 1):
        P_n = np.zeros((d, d), dtype=complex)
        for b in range(d):
            if bin(b).count('1') == n:
                P_n[b, b] = 1.0
        rank_n = int(np.real(np.trace(P_n @ P_n)))
        if rank_n == 0:
            continue
        coeff = np.real(np.trace(P_n @ rho)) / rank_n
        rho_d0 = rho_d0 + coeff * P_n
    return rho_d0


def numerical_pi2_odd_in_memory(rho, N):
    rho_d0 = kernel_projection_popcount(rho, N)
    rho_d2 = rho - rho_d0
    d = 2**N

    inv = 1.0 / d
    odd_memory_sq = 0.0
    for k in range(4**N):
        kk = k
        idxs = []
        for _ in range(N):
            idxs.append(kk & 3)
            kk >>= 2
        if (sum(BIT_B[i] for i in idxs) & 1) == 0:
            continue
        sigma = kron_n([PAULIS[i] for i in idxs])
        c_memory = np.trace(sigma @ rho_d2) * inv
        odd_memory_sq += abs(c_memory) ** 2 * d
    norm_memory_sq = np.linalg.norm(rho_d2, 'fro') ** 2
    return odd_memory_sq / norm_memory_sq if norm_memory_sq > 1e-14 else 0.0


def first_pair(N, n_p, n_q, hd):
    for p in range(2**N):
        if bin(p).count('1') != n_p:
            continue
        for q in range(2**N):
            if bin(q).count('1') != n_q:
                continue
            if p == q:
                continue
            if bin(p ^ q).count('1') == hd:
                return p, q
    return None, None


def main():
    print("F88 generalised closed form vs numerical verification")
    print("=" * 110)
    print(f"{'N':>3} {'n_p':>4} {'n_q':>4} {'HD':>3} {'category':<22} {'static':>9} {'α (closed)':>11} {'predict':>10} {'numeric':>10} {'dev':>10}")
    print("-" * 110)

    cases = []
    for N in range(2, 8):
        for n_p in range(N + 1):
            for n_q in range(n_p, N + 1):
                # HD compatible: |n_q − n_p| ≤ HD ≤ n_p + n_q, same parity as n_q − n_p
                hd_min = abs(n_q - n_p) if n_q != n_p else 2  # intra needs HD ≥ 2 (must differ)
                hd_max = min(n_p + n_q, N)
                # parity: HD ≡ |n_q − n_p| (mod 2)
                parity = (n_q - n_p) % 2
                for hd in range(hd_min, hd_max + 1):
                    if hd % 2 != parity:
                        continue
                    if hd == 0:
                        continue
                    cases.append((N, n_p, n_q, hd))

    max_dev = 0.0
    fail_count = 0
    for N, n_p, n_q, hd in cases:
        p, q = first_pair(N, n_p, n_q, hd)
        if p is None:
            continue
        rho = coherence_state(N, p, q, sign=+1)
        if rho is None:
            continue

        s = static_fraction_general(N, n_p, n_q)
        a = alpha_krawtchouk_general(N, n_p, n_q)
        pred = pi2_odd_in_memory_general(N, n_p, n_q, hd)
        num = numerical_pi2_odd_in_memory(rho, N)
        dev = abs(pred - num)
        max_dev = max(max_dev, dev)

        if hd == N:
            cat = "HD=N (Π²-classical)"
        elif n_p == n_q and 2 * n_p == N:
            cat = "intra-mirror"
        elif n_p == n_q:
            cat = "intra generic"
        elif n_p + n_q == N:
            cat = "inter popcount-mirror"
        elif (n_p == N // 2 or n_q == N // 2) and N % 2 == 0:
            cat = "inter K-intermediate"
        else:
            cat = "inter generic"

        marker = " " if dev < 1e-9 else " ✗"
        if dev >= 1e-9:
            fail_count += 1
        print(f"{N:>3} {n_p:>4} {n_q:>4} {hd:>3} {cat:<22} {s:>9.6f} {a:>11.6f} {pred:>10.6f} {num:>10.6f} {dev:>10.2e}{marker}")

    print()
    print(f"Total cases: {len(cases)}")
    print(f"Max deviation (predict − numeric): {max_dev:.2e}")
    print(f"Failures (dev ≥ 1e−9): {fail_count}")


if __name__ == "__main__":
    main()
