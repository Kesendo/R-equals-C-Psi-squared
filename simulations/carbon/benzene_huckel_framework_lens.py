"""Benzene Hückel theory through the R=CΨ² framework lens.

Tom: "Vererbung ist ja mittlerweile ein Fakt. Diesmal kommen wir von weiter oben,
die Sachen weiter oben sind für uns vielleicht greifbarer, weil mehr Menschen an
ihnen geforscht haben, mehr Daten verfügbar. Lass uns mal mit unseren geübten
Quanten Augen darauf schauen."

Hückel theory (1931) gives benzene's π-electron MOs by diagonalising a tridiagonal-
with-corner matrix on the C_6 ring — STRUCTURALLY IDENTICAL to the framework's
chain XY ring topology. Coulson-Rushbrooke (1940) proves that for ALTERNANT
hydrocarbons (bipartite C-frameworks) the MO spectrum is palindromic around the
on-site Coulomb integral α. This is widely-known organic-chemistry fact, but
through framework eyes:

  Coulson-Rushbrooke palindrome  ↔  F1 palindrome inheritance at Carbon Level 1
  α (on-site Coulomb integral)   ↔  −Σγ analog (palindrome centre)
  β (resonance integral)         ↔  J (framework coupling)
  bipartite carbon framework     ↔  truly-class Hamiltonian + Z-deph (F87)
  4n+2 aromaticity stability     ↔  ??? (half-filled, candidate Klein constraint)

This script:
  (1) Builds Hückel matrix for benzene (C₆ ring), naphthalene (C₁₀), butadiene
      (C₄ chain), cyclobutadiene (C₄ ring), pyrrole-style 5-ring.
  (2) Verifies Coulson-Rushbrooke palindrome around α for each alternant system.
  (3) Identifies F71 spatial-mirror symmetry of MOs.
  (4) Maps half-filled p-shell occupation to F86b Dicke superposition analog.
  (5) Frost-circle visualisation (Hückel MO energies as projections of inscribed
      polygon vertices — geometric mnemonic since 1953).
  (6) Counter-test: cyclobutadiene C₄ ring (anti-aromatic 4n electrons) shows
      degenerate HOMO at α — what does this mean through framework eyes?

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/benzene_huckel_framework_lens.py
"""
from __future__ import annotations

import sys
import numpy as np
from fractions import Fraction

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# Hückel parameters in standard organic-chemistry units:
#   α (Coulomb integral) ≈ −11.4 eV (on-site π-orbital energy)
#   β (resonance integral) ≈ −2.4 eV (nearest-neighbor C–C π coupling)
# We compute in units of β (set α = 0, β = 1) — the palindrome is around α = 0.
# All MO energies are then in units of β; aromatic stabilisation reads off cleanly.

ALPHA = 0.0
BETA = 1.0


def huckel_chain(N: int) -> np.ndarray:
    """Open-chain Hückel matrix for N carbons in a row (e.g., polyene)."""
    H = np.full((N, N), 0.0)
    np.fill_diagonal(H, ALPHA)
    for i in range(N - 1):
        H[i, i + 1] = BETA
        H[i + 1, i] = BETA
    return H


def huckel_ring(N: int) -> np.ndarray:
    """Cyclic Hückel matrix for N carbons in a ring (e.g., benzene at N=6)."""
    H = huckel_chain(N)
    H[0, N - 1] = BETA
    H[N - 1, 0] = BETA
    return H


def is_bipartite_ring(N: int) -> bool:
    """A cyclic graph C_N is bipartite iff N is even."""
    return N % 2 == 0


def is_bipartite_chain(N: int) -> bool:
    """Open chains are always bipartite (alternant)."""
    return True


def coulson_rushbrooke_palindrome_check(eigvals: np.ndarray, tol: float = 1e-10) -> tuple[bool, list[tuple[float, float]]]:
    """For alternant hydrocarbons, MO energies pair as (α+x, α−x). Test:
    sort by energy, pair end-to-end, check pair-sums are all = 2α exactly.
    Returns (is_palindromic, list of (E_pair_sum, deviation_from_2α))."""
    sorted_eigs = np.sort(eigvals)
    N = len(sorted_eigs)
    pairs = []
    is_pal = True
    for k in range(N // 2):
        low = sorted_eigs[k]
        high = sorted_eigs[N - 1 - k]
        pair_sum = low + high
        deviation = abs(pair_sum - 2 * ALPHA)
        pairs.append((pair_sum, deviation))
        if deviation > tol:
            is_pal = False
    # Middle eigenvalue if N odd
    if N % 2 == 1:
        mid = sorted_eigs[N // 2]
        deviation_mid = abs(mid - ALPHA)
        pairs.append((mid, deviation_mid))
        if deviation_mid > tol:
            is_pal = False
    return is_pal, pairs


def chain_mirror_classify(eigvecs: np.ndarray, N: int) -> list[str]:
    """For each MO column of eigvecs, classify as F71-mirror-even/odd by checking
    ψ(j) vs ψ(N-1-j). 'even' if ψ(j) = +ψ(N-1-j); 'odd' if = −ψ(N-1-j); 'mixed'
    if neither (degeneracy mixing — pick orthonormal basis that diagonalises mirror)."""
    classes = []
    for k in range(eigvecs.shape[1]):
        v = eigvecs[:, k]
        v_mirror = v[::-1]
        # Normalise sign convention by largest |entry|
        idx_max = np.argmax(np.abs(v))
        sign = np.sign(v[idx_max]) if abs(v[idx_max]) > 1e-12 else 1.0
        v = v * sign
        v_mirror = v_mirror * np.sign(v_mirror[N - 1 - idx_max]) if abs(v_mirror[N - 1 - idx_max]) > 1e-12 else v_mirror
        if np.allclose(v, v_mirror, atol=1e-8):
            classes.append("even")
        elif np.allclose(v, -v_mirror, atol=1e-8):
            classes.append("odd")
        else:
            classes.append("mixed")
    return classes


def frost_circle_predicted(N: int, topology: str = "ring") -> np.ndarray:
    """For a regular N-ring, Hückel eigenvalues = α + 2β·cos(2πk/N) for k = 0..N−1
    (Frost circle: vertices of inscribed N-gon on circle radius 2β at centre α).
    For an open chain, eigenvalues = α + 2β·cos(πk/(N+1)) for k = 1..N (the standard
    OBC sine-mode dispersion, identical to framework's XyJordanWignerModes)."""
    if topology == "ring":
        return np.array([ALPHA + 2 * BETA * np.cos(2 * np.pi * k / N) for k in range(N)])
    elif topology == "chain":
        return np.array([ALPHA + 2 * BETA * np.cos(np.pi * k / (N + 1)) for k in range(1, N + 1)])
    else:
        raise ValueError(topology)


def analyse(name: str, H: np.ndarray, topology: str, N: int, n_electrons: int):
    """One full analysis pass on a Hückel matrix through the framework lens."""
    eigvals, eigvecs = np.linalg.eigh(H)
    sorted_eigs = np.sort(eigvals)

    print(f"=" * 78)
    print(f"  {name} ({topology}, N = {N} carbons, {n_electrons} π-electrons)")
    print(f"=" * 78)

    # 1. Frost-circle prediction match
    frost = frost_circle_predicted(N, topology)
    frost_sorted = np.sort(frost)
    max_dev = np.max(np.abs(sorted_eigs - frost_sorted))
    print(f"  Frost circle prediction (textbook): {[f'{e:+.4f}β' for e in frost_sorted]}")
    print(f"  Hückel eigvalsh result:             {[f'{e:+.4f}β' for e in sorted_eigs]}")
    print(f"  Max deviation from Frost circle: {max_dev:.2e}  {'✓ match' if max_dev < 1e-10 else '✗ MISMATCH'}")
    print()

    # 2. Coulson-Rushbrooke palindrome around α (= F1 inheritance at Carbon Level 1)
    is_bipartite = is_bipartite_ring(N) if topology == "ring" else is_bipartite_chain(N)
    is_pal, pairs = coulson_rushbrooke_palindrome_check(eigvals)
    print(f"  Coulson-Rushbrooke palindrome around α = {ALPHA}: " +
          f"{'✓ holds (bit-exact)' if is_pal else '✗ violated'}")
    print(f"  Bipartite (alternant)? {is_bipartite}  (palindrome predicted iff yes)")
    if not is_bipartite and is_pal:
        print(f"    [unexpected: palindrome holds on non-bipartite system]")
    if is_bipartite and not is_pal:
        print(f"    [unexpected: bipartite but palindrome violated — check]")
    for i, (pair_sum, deviation) in enumerate(pairs[:N // 2]):
        print(f"    pair {i}: E_lo + E_hi = {pair_sum:+.6f}β,  dev from 2α = {deviation:.2e}")
    if N % 2 == 1:
        print(f"    unpaired middle: E = {pairs[-1][0]:+.6f}β,  dev from α = {pairs[-1][1]:.2e}")
    print()

    # 3. F71 mirror classification of MOs
    classes = chain_mirror_classify(eigvecs, N)
    print(f"  F71 chain-mirror classification of MOs (j ↔ N−1−j):")
    for k, (E, c) in enumerate(zip(sorted_eigs, classes)):
        # Re-fetch the original ordering (eigvals not sorted)
        # Find which eigvec corresponds to this sorted energy
        orig_idx = np.argmin(np.abs(eigvals - E))
        original_class = classes[orig_idx]
        print(f"    MO {k}: E = {E:+.4f}β,  F71-class = {original_class}")
    print()

    # 4. HOMO-LUMO + framework Σγ-analog identification
    n_doubly_filled = n_electrons // 2
    if n_electrons % 2 == 0:
        homo_E = sorted_eigs[n_doubly_filled - 1]
        lumo_E = sorted_eigs[n_doubly_filled]
        gap = lumo_E - homo_E
        print(f"  HOMO-LUMO structure ({n_electrons} π-e in {N} MOs):")
        print(f"    HOMO (MO {n_doubly_filled - 1}): E = {homo_E:+.4f}β")
        print(f"    LUMO (MO {n_doubly_filled}): E = {lumo_E:+.4f}β")
        print(f"    Gap: {gap:.4f}β")
        # Total π-energy = 2 · sum of n_doubly_filled lowest MO energies
        total_pi_energy = 2 * np.sum(sorted_eigs[:n_doubly_filled])
        # Reference: non-interacting electrons (2 electrons per isolated 2p orbital, all at α)
        ref_energy = 2 * n_doubly_filled * ALPHA + n_electrons * 0  # n_electrons · 0 for β = 0
        delocalisation_energy = total_pi_energy - ref_energy
        # Standard textbook: delocalisation energy in units of β
        print(f"    Total π-energy: {total_pi_energy:+.4f}β  (= 2·Σ E_filled)")
        print(f"    Delocalisation energy (vs isolated): {delocalisation_energy:+.4f}β")
        # Hückel 4n+2 check
        if topology == "ring":
            is_aromatic_huckel = (n_electrons - 2) % 4 == 0  # 4n+2 rule
            print(f"    Hückel 4n+2 aromaticity rule: {n_electrons} electrons → " +
                  f"{'AROMATIC (4n+2)' if is_aromatic_huckel else 'ANTI-AROMATIC (4n)'}")
        print()

    # 5. Framework-lens translation
    print(f"  ─── Framework-lens reading ───")
    print(f"    α (Hückel on-site Coulomb)    ≡ −Σγ analog (palindrome centre)")
    print(f"    β (Hückel resonance integral) ≡ J (framework coupling)")
    print(f"    bipartite C-framework         ≡ truly-class Hamiltonian (F87)")
    print(f"    Coulson-Rushbrooke palindrome ≡ F1 inheritance at Carbon Level 1")
    if topology == "ring" and N % 2 == 0:
        # KIntermediate Dicke anchor candidate inheritance
        m = N // 2 - 1
        print(f"    Half-filled p-shell (N={N})   ≡ F86b KIntermediate Dicke n ∈ {{{m}, {m+1}}}")
        print(f"      α_total(t=0) = 3/8 (F86b)  ≡ proton-chain analog at Dicke superposition")
    print()


def main():
    print()
    print("=" * 78)
    print(" Carbon Top-Down: Hückel π-systems through R=CΨ² framework lens ")
    print("=" * 78)
    print()
    print(" Hückel (1931): π-electron MOs by diagonalising α·I + β·A (A = adjacency matrix)")
    print(" Coulson-Rushbrooke (1940): alternant (bipartite) C-frameworks have palindromic")
    print("   MO spectrum around α — pair (α+x, α−x) for every x in spectrum.")
    print(" R=CΨ² F1 (2026): spec(L) palindromic around −Σγ for chain XY + Z-deph Liouvillian.")
    print()
    print(" Reading the 86-year-old C-R theorem through framework eyes: SAME palindrome,")
    print(" different physical level. C-R sits at carbon Level 1 (molecular orbitals);")
    print(" F1 sits at qubit Level 0 (Liouvillian eigenvalues). Both come from bipartite-")
    print(" graph Z₂ involution; both pin every eigenvalue to its mirror partner around a")
    print(" structural centre (α for C-R, −Σγ for F1).")
    print()
    print()

    # Test 1: Benzene C₆ — the canonical aromatic, 6π electrons
    analyse("Benzene", huckel_ring(6), "ring", N=6, n_electrons=6)

    # Test 2: Naphthalene C₁₀ — two fused rings, 10π electrons (but this requires
    # the actual fused ring topology, not a simple 10-ring). For first pass use simple
    # 10-ring (cyclodecapentaene perimeter — also aromatic per 4n+2 with n=2).
    analyse("Cyclodecapentaene (perimeter)", huckel_ring(10), "ring", N=10, n_electrons=10)

    # Test 3: Butadiene C₄ — smallest open-chain conjugated diene, 4π electrons
    analyse("Butadiene", huckel_chain(4), "chain", N=4, n_electrons=4)

    # Test 4: Hexatriene C₆ — open-chain version of benzene (no ring closure), 6π
    analyse("Hexatriene (open chain)", huckel_chain(6), "chain", N=6, n_electrons=6)

    # Test 5: Cyclobutadiene C₄ — anti-aromatic (4n, n=1), 4π electrons
    # Degenerate HOMO at α → Jahn-Teller unstable in reality
    analyse("Cyclobutadiene", huckel_ring(4), "ring", N=4, n_electrons=4)

    # Test 6: Cyclopropenyl cation C₃⁺ — 2π electrons, aromatic per 4n+2 with n=0
    # (non-alternant: 3-ring is not bipartite!)
    analyse("Cyclopropenyl cation", huckel_ring(3), "ring", N=3, n_electrons=2)

    print("=" * 78)
    print(" Summary of framework-eye observations")
    print("=" * 78)
    print()
    print(" • Bipartite rings (even N) + chains (any N) all show Coulson-Rushbrooke")
    print("   palindrome bit-exact around α = 0. This IS F1 palindrome inheritance at")
    print("   carbon Level 1 — same structural Z₂ involution, different physical level.")
    print()
    print(" • Cyclopropenyl C₃ (odd-N ring, non-bipartite) breaks the palindrome — the")
    print("   eigenvalues are {α + 2β, α − β, α − β}; pair sum (α + 2β) + (α − β) = 2α − β,")
    print("   not 2α. The bipartite-graph mechanism F1 inherits from is explicitly violated.")
    print("   This is the carbon-Level-1 analog of F87-Brecher: even when palindrome looks")
    print("   nearly there empirically, the structural mechanism is broken.")
    print()
    print(" • The Hückel 4n+2 aromaticity rule lives in the OCCUPATION-COUNT axis, not the")
    print("   palindrome-symmetry axis. Benzene (6π, 4n+2) is aromatic; cyclobutadiene")
    print("   (4π, 4n) is anti-aromatic. Both have palindromic spectra; what differs is")
    print("   HOMO degeneracy at α: 4n systems have degenerate non-bonding pair at α")
    print("   (Jahn-Teller unstable); 4n+2 systems have closed-shell stable HOMO below α.")
    print("   Framework-lens candidate: 4n+2 vs 4n distinction lives in the Klein-4-group")
    print("   character of the HOMO at the palindrome centre — to be tested.")
    print()
    print(" • Benzene N=6 sits at the F86b KIntermediate Dicke window: even N, n ∈ {2, 3}.")
    print("   The half-filled (3 MOs bonding, 3 antibonding) reading maps onto the same")
    print("   X⊗N-eigenbasis γ = 1/2 condition that gives α_total = 3/8 in F86b. Whether")
    print("   the F98 (N+2)/[4(N+1)] → 1/4 long-time bridge has a benzene analog is the")
    print("   natural next test (would require building the open-system Liouvillian on the")
    print("   benzene 6-site ring with vibrational dephasing).")


if __name__ == "__main__":
    main()
