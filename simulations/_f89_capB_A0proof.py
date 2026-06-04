#!/usr/bin/env python3
"""F89 capB: prove A(0)=m (k=2 mod4), A(0)=0 (4|k); B(0) closed form (2026-06-04).

A(0) = sum_{i=0}^{k/2} (-1)^i (k-4i),  k even, let K=k/2.
  = sum_{i=0}^{K} (-1)^i (k-4i).
This is an alternating arithmetic series.  Closed form:
  Let f(i)=k-4i.  sum_{i=0}^{K}(-1)^i f(i).  Pairing i=2t,2t+1: f(2t)-f(2t+1)=(k-8t)-(k-8t-4)=4.
  - K even (4|k): K+1 odd terms; (K/2) pairs each +4, plus last term f(K)=k-4K=k-2k=-k... wait K=k/2 so 4K=2k, f(K)=k-2k=-k.
    Actually with K even, terms i=0..K, last i=K even (+). pairs (0,1),(2,3),...,(K-2,K-1) = K/2 pairs *4, + f(K).
    f(K)=k-4*(k/2)=k-2k=-k.  So A0 = 4*(K/2) + (-k) = 2K - k = k - k = 0.  => A(0)=0 for 4|k. PROVEN.
  - K odd (k=2 mod4): terms i=0..K, K odd so last i=K is '-'. pairs (0,1),...,(K-1,K) = (K+1)/2 pairs each +4.
    A0 = 4*(K+1)/2 = 2(K+1) = 2(k/2+1) = k+2 = m.  => A(0)=m for k=2 mod4. PROVEN.

B(0) = sum_{i=0}^{K}(k-4i)^2.  Closed form via sum of squares of arithmetic progression.

This script verifies both A(0) closed forms and gives B(0) closed form + v2, confirming the
k=2 mod4 cap ingredient A(0)=m (=> 2-content 2 v2(m) from A(0)^2).
"""
from __future__ import annotations

import sys

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def v2(n) -> int:
    n = abs(int(n))
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def A0(k):
    return sum((-1) ** i * (k - 4 * i) for i in range(k // 2 + 1))


def B0(k):
    return sum((k - 4 * i) ** 2 for i in range(k // 2 + 1))


def main():
    print("=" * 96)
    print("F89 capB: A(0) closed forms PROVEN; B(0) closed form")
    print("=" * 96)
    print()
    print("Claim:  A(0) = m  (k=2 mod4),   A(0) = 0  (4|k).   [alternating arithmetic series]")
    print(f"  {'k':>3} {'kmod4':>5} {'A0_direct':>9} {'A0_claim':>8} {'ok':>3}")
    okall = True
    for k in range(6, 70, 2):
        a = A0(k)
        claim = (k + 2) if k % 4 == 2 else 0
        ok = (a == claim)
        okall = okall and ok
        if k <= 40 or not ok:
            print(f"  {k:>3} {k%4:>5} {a:>9} {claim:>8} {('OK' if ok else 'X'):>3}")
    print(f"  ... A(0) closed form holds for all even k in [6,68]: {okall}")
    print()

    # B(0) closed form: sum_{i=0}^{K}(k-4i)^2, K=k/2.  Get sympy closed form.
    K = sp.Symbol("K", integer=True, nonnegative=True)
    i = sp.Symbol("i", integer=True)
    kk = 2 * K
    expr = sp.summation((kk - 4 * i) ** 2, (i, 0, K))
    expr = sp.simplify(expr)
    print(f"B(0) = sum_{{i=0}}^{{K}}(2K-4i)^2 = {sp.factor(expr)}   (K=k/2)")
    print()
    # verify
    print(f"  {'k':>3} {'B0_direct':>9} {'B0_closed':>9} {'v2(B0)':>6} {'ok':>3}")
    okall2 = True
    for k in [6, 8, 10, 12, 14, 16, 20, 24, 32]:
        Kv = k // 2
        b = B0(k)
        bc = int(expr.subs(K, Kv))
        ok = (b == bc)
        okall2 = okall2 and ok
        print(f"  {k:>3} {b:>9} {bc:>9} {v2(b):>6} {('OK' if ok else 'X'):>3}")
    print(f"  B(0) closed form ok: {okall2}")
    print()
    print("  => For k=2 mod4:  A(0)=m, so A(0)^2 contributes 2 v2(m) to v2(Nint(0)).")
    print("     This is the chain-length 2-content carrier for the k=2 mod4 cap.")
    print("     (The full v2(Rm_top)=2v2(m)+2 still needs the trace's leading behavior,")
    print("      but A(0)=m is the elementary algebraic source of the v2(m) dependence.)")
    print()


if __name__ == "__main__":
    main()
