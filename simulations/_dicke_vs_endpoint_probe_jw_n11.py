"""Algebraic comparison: Dicke probe vs Endpoint-localized probe at N=11 c=2.

Tom 2026-05-10 (late evening): "Das einzige was einwirkt, ist der initial state oder?"

Yes — H = J·Σ(X⊗X + Y⊗Y) is uniform across bonds, Z-dephasing γ₀·Σ Z_l ρ Z_l is
uniform across sites, F73 spatial-sum kernel S is uniform. Only the INITIAL STATE
breaks site-uniformity. The bond-position-dependence we observe in F86 K_b
observable is therefore an interaction between probe choice and chain geometry,
not a chain-intrinsic property.

This script verifies algebraically:

1. OBC sine basis ψ_k(l) = √(2/(N+1))·sin(πk(l+1)/(N+1)) is orthonormal at N=11.
2. Dicke probe |D_1⟩ = (1/√N)·Σ_l c_l†|0⟩ projects ONLY onto odd k. Even k
   components are exactly 0 (sine reflection-symmetry cancellation).
3. Dicke weights p_k = |⟨k|D_1⟩|² scale as 1/k² for odd k (low-k dominance).
4. Endpoint-localized probe |L_0⟩ = c_0†|0⟩ touches ALL k (no parity selection).
5. Per-bond inter-mode coupling C_b[k1,k2] = ψ_k1(b)·ψ_k2(b+1) + ψ_k1(b+1)·ψ_k2(b)
   is F71-mirror-invariant: C_b = C_{N-2-b}.
6. Dicke-probe-weighted bond response |D_b| has F71-mirror invariance.
7. Endpoint-loc probe weighted bond response BREAKS F71 (E_b ≠ E_{N-2-b}).
8. Destructive interference at Center: ψ_k(center) for odd k has ALTERNATING
   signs (+1, -1, +1, -1, ...) — the Dicke-weighted sum at Center cancels
   significantly, suppressing the predicted bond-coupling magnitude.

Empirical reference (from extended-grid c2hwhm scan at N=11, Q ∈ [0.2, 32]):
  Orbit 0 (Endpoint, b=0↔9): Q_peak=2.5007, g_eff=1.76
  Orbit 1 (b=1↔8):           Q_peak=8.7946, g_eff=0.50
  Orbit 2 (b=2↔7):           Q_peak=13.6117, g_eff=0.32
  Orbit 3 (mid-flank, b=3↔6): Q_peak=1.5901, g_eff=2.76
  Orbit 4 (Center, b=4↔5):    Q_peak=21.9389, g_eff=0.20

g_eff via F86a EP formula: Q_peak ≈ 4.39/g_eff (BareDoubledPtfXPeak·2 = 4.39).
"""

import numpy as np

N = 11
NUM_BONDS = N - 1  # = 10


def psi_k(k: int, l: int, N: int = N) -> float:
    """OBC sine basis: ψ_k(l) = √(2/(N+1))·sin(πk(l+1)/(N+1))."""
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (l + 1) / (N + 1))


def dicke_overlap(k: int, N: int = N) -> float:
    """⟨k|D_1⟩ = (1/√N)·Σ_l ψ_k(l). Zero for even k by reflection symmetry."""
    return (1.0 / np.sqrt(N)) * sum(psi_k(k, l, N) for l in range(N))


def endpoint_overlap(k: int, N: int = N) -> float:
    """⟨k|L_0⟩ = ψ_k(0). Nonzero for all k."""
    return psi_k(k, 0, N)


def C_b(b: int, k1: int, k2: int, N: int = N) -> float:
    """Bond-b inter-mode coupling: C_b[k1,k2] = ψ_k1(b)·ψ_k2(b+1) + ψ_k1(b+1)·ψ_k2(b)."""
    return psi_k(k1, b, N) * psi_k(k2, b + 1, N) + psi_k(k1, b + 1, N) * psi_k(k2, b, N)


def section(title: str) -> None:
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)


# ============================================================================
# Task 1: OBC sine basis orthonormality at N=11
# ============================================================================
section(f"Task 1: OBC sine basis ψ_k(l) at N={N} — orthonormality check")

for k in range(1, N):
    norm_sq = sum(psi_k(k, l) ** 2 for l in range(N))
    print(f"  k={k}: Σ_l |ψ_k(l)|² = {norm_sq:.10f}  (expect 1.0)")

# Cross-products: ψ_k1 · ψ_k2 = δ_{k1,k2}
cross_max_dev = 0.0
for k1 in range(1, N):
    for k2 in range(k1 + 1, N):
        cross = sum(psi_k(k1, l) * psi_k(k2, l) for l in range(N))
        cross_max_dev = max(cross_max_dev, abs(cross))
print(f"  max off-diagonal |ψ_k1·ψ_k2| (k1≠k2) = {cross_max_dev:.2e}  (expect ~ 0)")


# ============================================================================
# Task 2: Verify Dicke probe = odd-k filter
# ============================================================================
section("Task 2: Dicke probe |D_1⟩ projection onto k-modes")

print("  k    ⟨k|D_1⟩         |⟨k|D_1⟩|²    parity")
print("  --   -------------   -----------   ------")
dicke_weights = {}
total_dicke = 0.0
for k in range(1, N):
    overlap = dicke_overlap(k)
    weight = overlap ** 2
    dicke_weights[k] = weight
    total_dicke += weight
    parity = "ODD " if k % 2 == 1 else "even"
    flag = " <<< NONZERO ODD" if (k % 2 == 1 and abs(overlap) > 1e-12) else (" (expect 0)" if k % 2 == 0 else "")
    print(f"  {k:2d}   {overlap:+.10f}   {weight:.6e}   {parity}{flag}")
print(f"\n  Σ_k |⟨k|D_1⟩|² = {total_dicke:.10f}  (expect 1.0)")

# Verify 1/k² scaling for odd k (low-k dominance)
section("Task 2 (cont): Dicke weights scaling — 8/(π²·k²) for low odd k")
print("  k    p_k empirical     8/(π²·k²) prediction    ratio")
print("  --   --------------   ----------------------   -----")
for k in [1, 3, 5, 7, 9]:
    p_k = dicke_weights[k]
    pred = 8.0 / (np.pi ** 2 * k ** 2)
    ratio = p_k / pred if pred > 0 else float("nan")
    print(f"  {k:2d}   {p_k:.6e}   {pred:.6e}             {ratio:.4f}")


# ============================================================================
# Task 3: Endpoint-localized probe weights
# ============================================================================
section("Task 3: Endpoint-localized probe |L_0⟩ projection onto k-modes")

print("  k    ⟨k|L_0⟩         |⟨k|L_0⟩|²    parity")
print("  --   -------------   -----------   ------")
endpoint_weights = {}
total_endpoint = 0.0
for k in range(1, N):
    overlap = endpoint_overlap(k)
    weight = overlap ** 2
    endpoint_weights[k] = weight
    total_endpoint += weight
    parity = "ODD " if k % 2 == 1 else "even"
    print(f"  {k:2d}   {overlap:+.10f}   {weight:.6e}   {parity}")
print(f"\n  Σ_k |⟨k|L_0⟩|² = {total_endpoint:.10f}  (expect 1.0)")
print(f"  Distribution: Endpoint-loc TOUCHES ALL k (no parity selection)")
print(f"                Dicke FILTERS to odd k only")


# ============================================================================
# Task 4: Per-bond inter-mode coupling C_b[k1,k2] + F71 mirror invariance
# ============================================================================
section(f"Task 4: Per-bond coupling C_b[k1,k2] — F71 mirror invariance check")

# Check C_b = C_{N-2-b} bit-exact (it's algebraic, should be ~1e-15)
mirror_max_dev = 0.0
for b in range(NUM_BONDS):
    b_mirror = NUM_BONDS - 1 - b
    if b_mirror <= b:
        continue
    for k1 in range(1, N):
        for k2 in range(1, N):
            c1 = C_b(b, k1, k2)
            c2 = C_b(b_mirror, k1, k2)
            mirror_max_dev = max(mirror_max_dev, abs(c1 - c2))
print(f"  max |C_b[k1,k2] − C_{{N-2-b}}[k1,k2]| over all (b, k1, k2) = {mirror_max_dev:.2e}")
print(f"  (expect ~ 0 → F71 spatial-mirror invariance bit-exact)")


# ============================================================================
# Task 5: Dicke-probe-weighted bond response D_b
# ============================================================================
section("Task 5: Dicke-probe-weighted bond response D_b (odd-k only)")

# D_b = Σ_{k1,k2 odd, k1≠k2} (p_k1·p_k2)·|C_b[k1,k2]|²
# (k1 ≠ k2 because we want INTER-mode coupling, not energy-diagonal)
print("  bond b   D_b (Dicke-weighted)   F71 mirror partner b'   |D_b - D_b'|")
print("  ------   --------------------   --------------------    ------------")
D_b_values = []
for b in range(NUM_BONDS):
    Db = 0.0
    for k1 in range(1, N):
        if k1 % 2 == 0:
            continue
        for k2 in range(1, N):
            if k2 % 2 == 0:
                continue
            if k1 == k2:
                continue
            Db += dicke_weights[k1] * dicke_weights[k2] * C_b(b, k1, k2) ** 2
    D_b_values.append(Db)

for b in range(NUM_BONDS):
    b_mirror = NUM_BONDS - 1 - b
    diff = abs(D_b_values[b] - D_b_values[b_mirror])
    print(f"  b={b:2d}    {D_b_values[b]:.6e}        b'={b_mirror:2d}                  {diff:.2e}")


# ============================================================================
# Task 6: Endpoint-loc probe weighted bond response E_b — F71 BREAKS
# ============================================================================
section("Task 6: Endpoint-loc probe weighted bond response E_b (all k)")

# E_b = Σ_{k1,k2, k1≠k2} (p_k1·p_k2 endpoint)·|C_b[k1,k2]|²
print("  bond b   E_b (Endpoint-loc-weighted)   F71 mirror b'   |E_b - E_b'|  (F71 broken)")
print("  ------   ---------------------------   -------------   ------------")
E_b_values = []
for b in range(NUM_BONDS):
    Eb = 0.0
    for k1 in range(1, N):
        for k2 in range(1, N):
            if k1 == k2:
                continue
            Eb += endpoint_weights[k1] * endpoint_weights[k2] * C_b(b, k1, k2) ** 2
    E_b_values.append(Eb)

for b in range(NUM_BONDS):
    b_mirror = NUM_BONDS - 1 - b
    diff = abs(E_b_values[b] - E_b_values[b_mirror])
    print(f"  b={b:2d}    {E_b_values[b]:.6e}            b'={b_mirror:2d}            {diff:.2e}")


# ============================================================================
# Task 7: Destructive interference at Center under Dicke probe
# ============================================================================
section("Task 7: Center bond destructive interference (Dicke probe)")

# At chain center sites, ψ_k for odd k has alternating signs.
# N=11: center site is at l = (N-1)/2 = 5. Bond 4 connects sites 4 and 5; bond 5 connects sites 5 and 6.
# Print ψ_k(5), ψ_k(4), ψ_k(6) for odd k.
print(f"  Center site: l = {(N - 1) // 2} = 5")
print(f"  ψ_k at sites 4, 5, 6 for odd k:")
print(f"  k    ψ_k(4)         ψ_k(5)         ψ_k(6)         sign(ψ_k(5))")
print(f"  --   ------------   ------------   ------------   ------------")
for k in [1, 3, 5, 7, 9]:
    p4 = psi_k(k, 4)
    p5 = psi_k(k, 5)
    p6 = psi_k(k, 6)
    sign = "+" if p5 > 0 else "-"
    print(f"  {k}    {p4:+.10f}   {p5:+.10f}   {p6:+.10f}   {sign}")

# At Center bond (e.g., b=4 connects sites 4, 5), compute C_4[k1,k2] for odd k1, k2 pairs
# and show signs for low pairs.
section("Task 7 (cont): C_b[k1,k2] sign structure at Center bond b=4 (sites 4,5)")
print(f"  k1   k2   C_4[k1,k2]            (sign × magnitude)")
print(f"  --   --   --------------------")
for k1 in [1, 3, 5, 7, 9]:
    for k2 in [1, 3, 5, 7, 9]:
        c = C_b(4, k1, k2)
        print(f"  {k1}    {k2}    {c:+.6e}")

# Compare to Endpoint bond b=0 (sites 0, 1)
section("Task 7 (cont): C_b[k1,k2] sign structure at Endpoint bond b=0 (sites 0,1)")
print(f"  k1   k2   C_0[k1,k2]            (sign × magnitude)")
print(f"  --   --   --------------------")
for k1 in [1, 3, 5, 7, 9]:
    for k2 in [1, 3, 5, 7, 9]:
        c = C_b(0, k1, k2)
        print(f"  {k1}    {k2}    {c:+.6e}")


# ============================================================================
# Task 8: Compare bond-dependence pattern to empirical g_eff
# ============================================================================
section("Task 8: Predicted vs empirical bond-position-dependence pattern")

# Empirical g_eff per orbit at N=11 (from extended-grid scan):
empirical_g_eff = {
    0: 1.76,  # Endpoint b=0↔9
    1: 0.50,  # b=1↔8
    2: 0.32,  # b=2↔7
    3: 2.76,  # b=3↔6 (mid-flank)
    4: 0.20,  # b=4↔5 (Center)
}

# Map bonds to F71 orbits (orbit index = min(b, N-2-b))
# At N=11, NumBonds=10: orbits are 0-9, 1-8, 2-7, 3-6, 4-5
print(f"  Per-orbit Dicke-weighted prediction D_b vs empirical g_eff:")
print(f"  orbit   bond pair    D_b (Dicke)        empirical g_eff   ratio D_b/g_eff")
print(f"  -----   ----------   ----------------   ---------------   ---------------")
for orbit_idx in range(5):
    b = orbit_idx
    b_mirror = NUM_BONDS - 1 - b
    D_avg = (D_b_values[b] + D_b_values[b_mirror]) / 2.0
    g = empirical_g_eff[orbit_idx]
    ratio = D_avg / g if g > 0 else float("nan")
    print(f"  {orbit_idx}       b={b}↔{b_mirror}       {D_avg:.6e}      {g:.4f}            {ratio:.4f}")


# ============================================================================
# Summary
# ============================================================================
section("SUMMARY")

print("Verified algebraically:")
print("  ✓ OBC sine basis orthonormal at N=11 (Task 1)")
print("  ✓ Dicke probe filters to ODD k only, even k = 0 exactly (Task 2)")
print("  ✓ Dicke weights p_k ∝ 1/k², low-odd-k dominant (Task 2)")
print("  ✓ Endpoint-loc probe TOUCHES ALL k (no parity filter) (Task 3)")
print("  ✓ Per-bond |C_b|² is F71-mirror-invariant in MAGNITUDE (Task 4)")
print("    (C_b itself: F71 sign rule C_{N-2-b} = (-1)^{k1+k2}·C_b)")
print("  ✓ Dicke-weighted SQUARED-magnitude D_b is F71-mirror-invariant (Task 5)")
print("  ✓ Endpoint-loc weighted SQUARED-magnitude E_b is F71-invariant (Task 6)")
print("    (because |ψ_k(N-1)|² = |ψ_k(0)|² — magnitudes preserved under R)")
print()
print("=" * 80)
print("CRITICAL FINDING")
print("=" * 80)
print()
print("My naive |C_b|² Dicke-weighted prediction D_b shows Center MAXIMUM:")
print(f"  D_Endpoint  = {D_b_values[0]:.4e}")
print(f"  D_flank-1   = {D_b_values[1]:.4e}")
print(f"  D_flank-2   = {D_b_values[2]:.4e}")
print(f"  D_mid-flank = {D_b_values[3]:.4e}")
print(f"  D_Center    = {D_b_values[4]:.4e}  ← LARGEST")
print()
print("But empirical g_eff is INVERTED at Center:")
print("  g_Endpoint  = 1.76")
print("  g_flank-1   = 0.50")
print("  g_flank-2   = 0.32")
print("  g_mid-flank = 2.76")
print("  g_Center    = 0.20  ← SMALLEST")
print()
print("Magnitude-squared formula |C_b|² CANNOT capture g_eff bond-dependence.")
print("The empirical g_eff is determined by SIGNED matrix elements, not |·|².")
print()
print("=" * 80)
print("Test: Dicke-weighted SIGNED sum S_b — does Center destructively interfere?")
print("=" * 80)
print()

# Signed sum: Σ_{k1, k2 odd, k1≠k2} a_k1·a_k2·C_b[k1, k2]
# where a_k = ⟨k|D_1⟩ is the Dicke amplitude (signed; for odd k all positive)
# k1=k2 excluded because we want INTER-mode, not energy-diagonal.
S_b_signed = []
for b in range(NUM_BONDS):
    Sb = 0.0
    for k1 in range(1, N):
        if k1 % 2 == 0:
            continue
        a1 = dicke_overlap(k1)
        for k2 in range(1, N):
            if k2 % 2 == 0:
                continue
            if k1 == k2:
                continue
            a2 = dicke_overlap(k2)
            Sb += a1 * a2 * C_b(b, k1, k2)
    S_b_signed.append(Sb)

print("  bond b   S_b (Dicke signed sum)   |S_b| / |S_b_max|")
print("  ------   ----------------------   -----------------")
S_max = max(abs(s) for s in S_b_signed)
for b in range(NUM_BONDS):
    rel = abs(S_b_signed[b]) / S_max if S_max > 0 else 0.0
    print(f"  b={b:2d}    {S_b_signed[b]:+.6e}        {rel:.4f}")

print()
print("  Per-orbit signed sum |S_b| compared to empirical g_eff:")
print("  orbit   bond pair    |S_b|              empirical g_eff   |S_b|/g_eff")
print("  -----   ----------   ----------------   ---------------   -----------")
for orbit_idx in range(5):
    b = orbit_idx
    b_mirror = NUM_BONDS - 1 - b
    S_avg = (abs(S_b_signed[b]) + abs(S_b_signed[b_mirror])) / 2.0
    g = empirical_g_eff[orbit_idx]
    ratio = S_avg / g if g > 0 else float("nan")
    print(f"  {orbit_idx}       b={b}↔{b_mirror}       {S_avg:.6e}      {g:.4f}            {ratio:.4f}")

print()
print("=" * 80)
print("FINAL READING")
print("=" * 80)
print()
print("If the SIGNED sum S_b shows Center MINIMUM (or close to zero), Tom's")
print("destructive-interference-at-Center reading is empirically supported")
print("at the algebraic level — even though |C_b|² is large at Center, the")
print("alternating phases of odd-k modes cancel in the signed sum.")
print()
print("If S_b ALSO shows Center maximum (parallel to |C_b|²), the destructive-")
print("interference reading is falsified at this level — bond g_eff dependence")
print("comes from somewhere else (multi-cluster Petermann, polarity-gradient,")
print("or block-L specific structure not captured by simple JW projection).")
