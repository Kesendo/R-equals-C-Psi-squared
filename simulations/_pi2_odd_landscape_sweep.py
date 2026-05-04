"""F88 Π²-odd/memory landscape across popcount-coherence configurations.

For each popcount-coherence state |ψ⟩ = (|p⟩+|q⟩)/√2 with popcount(p) = n_p,
n_q = n_p+1, sweep N = 3..10 and verify the closed form

    Π²-odd / memory  =  (1/2 − α · s) / (1 − s)

against numerical MemoryAxisRho-equivalent at N = 3..7 (max deviation
machine precision). α is the Π²-odd-Frobenius²-fraction of the kernel
projection; computed both via the explicit Krawtchouk sum and via a
three-anchor closed form, cross-checked against each other.

Canonical statement, derivation, and verified table: see
docs/proofs/PROOF_F86_QPEAK.md §Structural inheritance from F88. C#
counterpart: compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs.
"""
from __future__ import annotations

import sys
from math import comb
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def krawtchouk(n: int, s: int, N: int) -> int:
    """K_n(s; N) = Σ_k (−1)^k · C(s, k) · C(N−s, n−k). Integer-valued.

    Valid k range: max(0, n − (N − s)) ≤ k ≤ min(s, n) so that both binomials
    have non-negative arguments.
    """
    k_lo = max(0, n - (N - s))
    k_hi = min(s, n)
    return sum((-1) ** k * comb(s, k) * comb(N - s, n - k) for k in range(k_lo, k_hi + 1))


def static_frac_closed(N: int, np_: int, nq: int) -> float:
    """static_frac = 1/(4·C(N, n_p)) + 1/(4·C(N, n_q))."""
    return 1.0 / (4 * comb(N, np_)) + 1.0 / (4 * comb(N, nq))


def alpha_pi2odd_static_closed(N: int, np_: int, nq: int) -> float:
    """Π²-odd-fraction-of-static via Krawtchouk closed form.

    α = (Σ_{s odd} C(N,s) (A_s+B_s)²) / (Σ_s C(N,s) (A_s+B_s)²)
    """
    Cn_p = comb(N, np_)
    Cn_q = comb(N, nq)
    coeffs_sq = []
    for s in range(N + 1):
        A_s = krawtchouk(np_, s, N) / Cn_p
        B_s = krawtchouk(nq, s, N) / Cn_q
        coeffs_sq.append(comb(N, s) * (A_s + B_s) ** 2)
    total = sum(coeffs_sq)
    if total == 0:
        return 0.0
    odd_sum = sum(coeffs_sq[s] for s in range(1, N + 1, 2))
    return odd_sum / total


def predict_pi2odd_in_memory(N: int, np_: int, nq: int) -> float:
    s = static_frac_closed(N, np_, nq)
    a = alpha_pi2odd_static_closed(N, np_, nq)
    return (0.5 - a * s) / (1.0 - s)


def alpha_three_anchor_closed(N: int, np_: int, nq: int) -> float:
    """Three-anchor closed form for α(N, n_p, n_q) at popcount-coherence (n_p, n_p+1):

    - popcount-mirror (n_p + n_q = N, odd N central pair): α = 0
    - near-mirror near-half (even N, n_p ∈ {N/2−1, N/2}): α = (N+2)/(4(N+1))
    - elsewhere: α = 1/2

    Empirical fit across N = 3..10. The two structured cases are derived; the
    α = 1/2 generic case is verified bit-exact on the sweep but not fully
    proven analytically here (one boundary case (0, 1) is proven via Krawtchouk
    moment identities; general proof is open).
    """
    if np_ + nq == N:
        return 0.0
    near_mirror_half_pairs = {(N // 2 - 1, N // 2), (N // 2, N // 2 + 1)}
    if N % 2 == 0 and (np_, nq) in near_mirror_half_pairs:
        return (N + 2) / (4 * (N + 1))
    return 0.5


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


def coherence_state(N, p_bits, q_bits):
    psi = np.zeros(2**N, dtype=complex)
    psi[p_bits] = 1.0 / np.sqrt(2.0)
    psi[q_bits] = 1.0 / np.sqrt(2.0)
    return np.outer(psi, psi.conj())


def numerical_diagnostic(rho, N):
    """Returns (static_frac, memory_frac, alpha, pi2odd_in_memory)."""
    d = 2**N
    rho_d0 = np.zeros_like(rho)
    sectors_used = 0
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
        sectors_used += 1
    rho_d2 = rho - rho_d0

    # Π²-odd Frobenius² of static and memory parts
    inv = 1.0 / d
    odd_static_norm_sq = 0.0
    odd_memory_norm_sq = 0.0
    for k in range(4**N):
        kk = k
        idxs = []
        for _ in range(N):
            idxs.append(kk & 3)
            kk >>= 2
        if (sum(BIT_B[i] for i in idxs) & 1) == 0:
            continue
        sigma = kron_n([PAULIS[i] for i in idxs])
        c_static = np.trace(sigma @ rho_d0) * inv
        c_memory = np.trace(sigma @ rho_d2) * inv
        odd_static_norm_sq += abs(c_static) ** 2 * d
        odd_memory_norm_sq += abs(c_memory) ** 2 * d

    norm_static_sq = np.linalg.norm(rho_d0, 'fro') ** 2
    norm_memory_sq = np.linalg.norm(rho_d2, 'fro') ** 2
    norm_total_sq = np.linalg.norm(rho, 'fro') ** 2

    return (
        norm_static_sq / norm_total_sq,
        norm_memory_sq / norm_total_sq,
        odd_static_norm_sq / norm_static_sq if norm_static_sq > 1e-14 else 0.0,
        odd_memory_norm_sq / norm_memory_sq if norm_memory_sq > 1e-14 else 0.0,
    )


def first_hd_pair(N, np_, nq, hd_target):
    """Find first (p_bits, q_bits) with popcount(p)=n_p, popcount(q)=n_q, HD(p,q)=hd_target."""
    for p in range(2**N):
        if bin(p).count('1') != np_:
            continue
        for q in range(2**N):
            if bin(q).count('1') != nq:
                continue
            if bin(p ^ q).count('1') == hd_target:
                return p, q
    return None, None


def main():
    print("F88 Π²-odd/memory landscape across popcount-coherence (n_p, n_p+1) configurations")
    print("=" * 96)
    print()
    print(f"{'N':>2} {'(n_p, n_q)':>10} {'mirror':>7} {'static':>10} {'α (closed)':>11} {'α (numeric)':>12} {'odd/mem (predict)':>18} {'odd/mem (numeric)':>18}")
    print("-" * 96)

    rows = []
    N_NUMERIC_MAX = 7  # numerical verification only up to N=7 (4^N Pauli loop)
    for N in range(3, 11):
        for np_ in range(N):
            nq = np_ + 1
            mirror = (np_ + nq == N)
            s_closed = static_frac_closed(N, np_, nq)
            a_closed = alpha_pi2odd_static_closed(N, np_, nq)
            pred = predict_pi2odd_in_memory(N, np_, nq)

            if N <= N_NUMERIC_MAX:
                p, q = first_hd_pair(N, np_, nq, hd_target=1)
                rho = coherence_state(N, p, q)
                s_num, m_num, a_num, om_num = numerical_diagnostic(rho, N)
            else:
                a_num, om_num = float('nan'), float('nan')

            rows.append((N, np_, nq, mirror, s_closed, a_closed, a_num, pred, om_num))
            mark = "*" if mirror else " "
            num_a = f"{a_num:>12.6f}" if a_num == a_num else f"{'(skip)':>12}"
            num_om = f"{om_num:>18.6f}" if om_num == om_num else f"{'(skip)':>18}"
            print(f"{N:>2} ({np_:>2}, {nq:>2})  {mark:>4}    {s_closed:>10.6f} {a_closed:>11.6f} {num_a} {pred:>18.6f} {num_om}")

    print()
    print("Legend: * = popcount-mirror (n_p + n_q = N → α = 0 by Krawtchouk reflection)")
    print()

    print(f"Formula verification (predict − numeric, only N ≤ {N_NUMERIC_MAX}):")
    max_dev = 0.0
    for (N, np_, nq, mirror, s, a_c, a_n, pred, om) in rows:
        if om != om:
            continue  # skipped row
        dev = abs(pred - om)
        if dev > max_dev:
            max_dev = dev
        if dev > 1e-9:
            print(f"  N={N} ({np_},{nq}): predict={pred:.10f} numeric={om:.10f} dev={dev:.2e}")
    print(f"Max deviation across verified rows: {max_dev:.2e}")

    print()
    print("α-closed vs α-numerical (Krawtchouk validation):")
    max_a_dev = 0.0
    for (N, np_, nq, mirror, s, a_c, a_n, pred, om) in rows:
        if a_n != a_n:
            continue
        dev = abs(a_c - a_n)
        if dev > max_a_dev:
            max_a_dev = dev
    print(f"Max α deviation: {max_a_dev:.2e}")

    # Identify special cases and intermediate-α rationals
    print()
    print("Anchor-α cases (popcount-mirror, half-half) + intermediate rationals:")
    from fractions import Fraction
    for (N, np_, nq, mirror, s, a_c, a_n, pred, om) in rows:
        if abs(a_c) < 1e-12:
            tag = "α = 0       (mirror, full Π²-odd cancellation in static)"
        elif abs(a_c - 0.5) < 1e-12:
            tag = "α = 1/2     (half-half static)"
        else:
            f = Fraction(a_c).limit_denominator(10000)
            tag = f"α = {f}     (intermediate; near-mirror at even N)"
        print(f"  N={N} ({np_:>2}, {nq:>2})  s = {Fraction(s).limit_denominator(10000)}   {tag}")

    # Three-anchor closed-form vs Krawtchouk closed-form
    print()
    print("Three-anchor formula vs Krawtchouk closed-form:")
    max_3anchor_dev = 0.0
    for (N, np_, nq, mirror, s, a_c, a_n, pred, om) in rows:
        a_3 = alpha_three_anchor_closed(N, np_, nq)
        dev = abs(a_c - a_3)
        if dev > max_3anchor_dev:
            max_3anchor_dev = dev
        if dev > 1e-12:
            print(f"  N={N} ({np_}, {nq}): Krawtchouk={a_c:.10f}, 3-anchor={a_3:.10f}, dev={dev:.2e}")
    print(f"Max 3-anchor deviation across all (N, n_p, n_p+1) pairs N=3..10: {max_3anchor_dev:.2e}")

    print()
    print("2n+1 vs N pattern (anchor categories per row):")
    for (N, np_, nq, mirror, s, a_c, a_n, pred, om) in rows:
        gap = N - (np_ + nq)
        if abs(a_c) < 1e-12:
            cat = "MIRROR  (n+nq = N)"
        elif abs(a_c - 0.5) < 1e-12:
            cat = "HALF    "
        else:
            cat = "INTER   "
        print(f"  N={N:>2} (n,n+1)=({np_:>2},{nq:>2})  n+nq = {np_+nq:>2}  N−(n+nq) = {gap:>+2}  α = {a_c:.6f}  [{cat}]")


if __name__ == "__main__":
    main()
