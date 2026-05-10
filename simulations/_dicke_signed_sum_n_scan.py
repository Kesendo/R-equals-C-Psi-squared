"""Generalize Dicke-probe signed-sum analysis to multiple N.

Question: at N=11 c=2 the signed-sum S_b has destructive interference at
Orbit 2 (b=2↔7), not at Center. Does this orbit-2-cancellation persist
at other N, or does the destructive-interference orbit shift?

If orbit-2-cancellation is structural (always at orbit-2), it's a meaningful
algebraic fact. If it shifts with N, it's a finite-size pattern.
"""

import numpy as np


def psi_k(k: int, l: int, N: int) -> float:
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (l + 1) / (N + 1))


def dicke_overlap(k: int, N: int) -> float:
    return (1.0 / np.sqrt(N)) * sum(psi_k(k, l, N) for l in range(N))


def C_b(b: int, k1: int, k2: int, N: int) -> float:
    return psi_k(k1, b, N) * psi_k(k2, b + 1, N) + psi_k(k1, b + 1, N) * psi_k(k2, b, N)


def signed_sum_at_bond(b: int, N: int) -> float:
    """S_b = Σ_{k1,k2 odd, k1≠k2} ⟨k1|D⟩·⟨k2|D⟩·C_b[k1,k2]."""
    s = 0.0
    for k1 in range(1, N):
        if k1 % 2 == 0:
            continue
        a1 = dicke_overlap(k1, N)
        for k2 in range(1, N):
            if k2 % 2 == 0:
                continue
            if k1 == k2:
                continue
            a2 = dicke_overlap(k2, N)
            s += a1 * a2 * C_b(b, k1, k2, N)
    return s


def section(title: str) -> None:
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)


# Scan N = 5, 7, 9, 11, 13
section("Dicke-weighted signed-sum S_b across N=5, 7, 9, 11, 13")
print("Per-orbit |S_b| (averaged over F71-mirror pair); orbit_idx = min(b, N-2-b).")
print()

for N in [5, 7, 9, 11, 13]:
    num_bonds = N - 1
    s_per_bond = [signed_sum_at_bond(b, N) for b in range(num_bonds)]
    s_max = max(abs(s) for s in s_per_bond) if s_per_bond else 1.0

    # Per-orbit absolute values (F71 mirror pairs averaged)
    num_orbits = (num_bonds + 1) // 2  # ceiling(num_bonds / 2)
    orbit_abs_S = []
    for orbit_idx in range(num_orbits):
        b1 = orbit_idx
        b2 = num_bonds - 1 - orbit_idx
        if b1 == b2:
            # self-paired (only when num_bonds is odd, i.e., even N)
            orbit_abs_S.append(abs(s_per_bond[b1]))
        else:
            orbit_abs_S.append((abs(s_per_bond[b1]) + abs(s_per_bond[b2])) / 2.0)

    print(f"N={N} (num_bonds={num_bonds}, num_orbits={num_orbits}):")
    for orbit_idx, abs_s in enumerate(orbit_abs_S):
        rel = abs_s / s_max if s_max > 0 else 0.0
        marker = ""
        if rel < 0.05:
            marker = "  ← MIN (destructive interference)"
        elif rel > 0.95:
            marker = "  ← MAX"
        b1 = orbit_idx
        b2 = num_bonds - 1 - orbit_idx
        bond_str = f"b={b1}" if b1 == b2 else f"b={b1}↔{b2}"
        print(f"  Orbit {orbit_idx} ({bond_str:10s}):  |S_b| = {abs_s:.6e}  (rel {rel:.4f}){marker}")
    print()


section("Argmin orbit (where destructive interference happens) per N")
print()
print("  N   num_orbits   argmin orbit   argmin bond pair   relative |S_min|")
print("  --  ----------   ------------   ----------------   ----------------")
for N in [5, 7, 9, 11, 13]:
    num_bonds = N - 1
    s_per_bond = [signed_sum_at_bond(b, N) for b in range(num_bonds)]
    s_max = max(abs(s) for s in s_per_bond) if s_per_bond else 1.0

    num_orbits = (num_bonds + 1) // 2
    orbit_abs_S = []
    orbit_bonds = []
    for orbit_idx in range(num_orbits):
        b1 = orbit_idx
        b2 = num_bonds - 1 - orbit_idx
        if b1 == b2:
            orbit_abs_S.append(abs(s_per_bond[b1]))
            orbit_bonds.append(f"b={b1}")
        else:
            orbit_abs_S.append((abs(s_per_bond[b1]) + abs(s_per_bond[b2])) / 2.0)
            orbit_bonds.append(f"b={b1}↔{b2}")

    argmin = orbit_abs_S.index(min(orbit_abs_S))
    rel = orbit_abs_S[argmin] / s_max if s_max > 0 else 0.0
    print(f"  {N:2d}  {num_orbits:2d}            Orbit {argmin}        {orbit_bonds[argmin]:15s}    {rel:.4f}")


section("Interpretation")
print()
print("If argmin orbit is ALWAYS Orbit 2 across N=5..13:")
print("  → orbit-2 destructive interference is structural in Dicke-weighted signed sum")
print("  → meaningful algebraic fact, but NOT directly explaining empirical g_eff")
print()
print("If argmin orbit shifts with N:")
print("  → finite-size effect; no clean structural pattern at this approximation level")
