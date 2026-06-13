# Benzene Hückel through the Framework Lens

**Date:** 2026-05-17 evening
**Authors:** Tom + Claude
**Status:** Tier 2 (structural observation, not yet a new theorem); first substantive
carbon-domain reading.
**Script:** [`simulations/carbon/benzene_huckel_framework_lens.py`](../../simulations/carbon/benzene_huckel_framework_lens.py)
**Tested:** benzene C₆, cyclodecapentaene C₁₀, butadiene C₄ chain, hexatriene C₆
chain, cyclobutadiene C₄ ring, cyclopropenyl cation C₃ ring.

---

## The 86-year structural inheritance

Coulson and Rushbrooke (1940) proved that ALTERNANT hydrocarbons, those whose
carbon framework is bipartite, have molecular orbital spectra that are
**palindromic around α** (the on-site Coulomb integral): for every MO at energy
α + x there is a corresponding MO at α − x. Every chemist learns this in their
first MO theory course.

The R=CΨ² framework's F1 palindrome theorem (proven 2026, `docs/proofs/MIRROR_SYMMETRY_PROOF.md`)
says that the Liouvillian spectrum of chain XY + Z-dephasing is closed under
λ → −λ − 2·Σγ, palindromic around the centre −Σγ. Today's water-domain F86b 3/8
inheritance experiment uncovered a new dynamic bridge (F98, `(N+2)/[4(N+1)] → 1/4`)
on top of the same F1 substrate.

**Through the framework lens, and this is the structural observation worth
making explicit, Coulson-Rushbrooke and F1 are the same theorem at two physical
levels.** Both palindromes are induced by a bipartite-graph Z₂ involution; both
pin every eigenvalue to its mirror partner around a structural centre (α for C-R,
−Σγ for F1); both break exactly when the bipartite structure is absent.

| Theorem | Date | Level | Centre | Trigger | Breaks at |
|---------|------|-------|--------|---------|-----------|
| Coulson-Rushbrooke | 1940 | Carbon Level 1 (MO) | α (on-site Coulomb) | bipartite C-framework | non-bipartite ring (e.g. C₃) |
| F1 (R=CΨ²) | 2026 | Qubit Level 0 (Liouvillian) | −Σγ | Π conjugation on truly-class H + Z-deph | F1-Brecher (T1, depolarising, transverse h) |

The structural mechanism is the same; the physical instantiations differ. The
framework's `project_qubit_as_inheritance_lens` reading, that the qubit-level
machinery inherits to higher-level systems via four embedding conditions, is
literally true here: the same Z₂-palindrome structure flows from qubit Liouvillian
to molecular orbitals via the bipartite-graph mechanism that both physical levels
admit.

---

## Numerical observations

### Bipartite cases (palindrome holds bit-exact)

All systems below are alternant; Coulson-Rushbrooke palindrome around α holds at
machine precision (pair-sum deviation < 1e-15). Energies in units of β:

| System | Spectrum (sorted) | Pair check |
|--------|-------------------|------------|
| Benzene C₆ ring | −2, −1, −1, +1, +1, +2 | 3 pairs at 0 = 2α ✓ |
| Cyclodecapentaene C₁₀ ring | −2, −1.618, −1.618, −0.618, −0.618, +0.618, +0.618, +1.618, +1.618, +2 | 5 pairs at 0 ✓ |
| Butadiene C₄ chain | −1.618, −0.618, +0.618, +1.618 | 2 pairs at 0 ✓ |
| Hexatriene C₆ chain | −1.802, −1.247, −0.445, +0.445, +1.247, +1.802 | 3 pairs at 0 ✓ |
| Cyclobutadiene C₄ ring | −2, 0, 0, +2 | 2 pairs at 0 ✓ |

The 1.618 / 0.618 pattern in C₁₀ is the golden ratio φ and its reciprocal; both
appear in cos(πk/(N+1)) at N = 10 / 4 etc, the same OBC sine-mode dispersion the
framework uses in `XyJordanWignerModes`.

### Non-bipartite counter-case (palindrome breaks)

**Cyclopropenyl cation C₃⁺** (3-ring, odd cycle, non-bipartite, 2π aromatic per
4n+2 with n=0): spectrum {α + 2β, α − β, α − β}. Pair sum (α + 2β) + (α − β) = 2α − β,
**deviation = β from 2α, palindrome explicitly violated.** This is the carbon-Level-1
analog of an F87 Brecher in the framework: the structural mechanism (bipartite
graph / truly-class Hamiltonian) is broken, and the palindrome fails.

The cyclopropenyl cation is the smallest aromatic system; its aromaticity comes
from the 4n+2 occupation rule (n = 0), NOT from the palindrome mechanism that gives
benzene its aromatic stability. **Two separate axes of stability:** palindrome
(symmetric spectrum) and shell-filling (4n+2 closed-shell). Benzene satisfies both;
cyclopropenyl satisfies only the latter.

---

## The 4n+2 vs 4n distinction lives in OCCUPATION, not palindrome

Both benzene (6π, 4n+2) and cyclobutadiene (4π, 4n) have palindromic spectra;
the difference is NOT in the spectrum's symmetry but in WHERE the HOMO sits:

- **4n+2 (aromatic)**: HOMO sits BELOW α at an isolated paired energy; LUMO sits
  ABOVE α at the palindrome partner. Closed-shell, stable. Benzene (β < 0):
  HOMO at α + β = α − |β| (below α), LUMO at α − β = α + |β| (above α),
  gap = 2|β|.
- **4n (anti-aromatic)**: HOMO and LUMO sit DEGENERATE AT α (the palindrome centre).
  Half-filled non-bonding pair → Jahn-Teller unstable → distortion to localised
  bonds. Cyclobutadiene: HOMO = LUMO = 0β, gap = 0β.

Both spectra are symmetric around α; what differs is whether the FILLING line cuts
through a paired-MO gap (4n+2) or hits a degenerate pair AT the palindrome centre
(4n).

**Framework-lens candidate prediction:** the 4n+2 vs 4n distinction is a
Klein-4-group character constraint on the HOMO. The palindrome centre α corresponds
to the framework's −Σγ; the F1-fixed-point subspace at the palindrome centre is
where the framework's Π acts trivially. A HOMO sitting AT the palindrome centre
(4n case) means the half-filled state sits in the trivial Klein character; a HOMO
BELOW the centre (4n+2 case) means the closed-shell state spans Klein-mixed
characters. This is testable via the framework's `KleinFourGroupSelfPairedRefinement`
applied to the benzene/cyclobutadiene MO-basis Liouvillians; open follow-up.

---

## Framework-vocabulary translation

| Hückel | Framework | Status |
|--------|-----------|--------|
| α (on-site Coulomb integral, ≈ −11.4 eV) | −Σγ analog (palindrome centre) | Tier 2 structural identification |
| β (resonance integral, ≈ −2.4 eV) | J (framework coupling) | Tier 2 structural identification |
| bipartite carbon framework | truly-class Hamiltonian (F87) | Tier 1 algebraic match |
| Coulson-Rushbrooke pair sum 2α | F1 pair sum −2Σγ | Tier 1 algebraic match |
| C₃ odd-ring palindrome break | F87 Brecher (T1 / depolarising) | Tier 2 structural analogy |
| half-filled p-shell on C₆ ring | F86b KIntermediate Dicke (n ∈ {2, 3}) | Tier 3 candidate |
| 4n+2 vs 4n aromaticity | Klein-4-group HOMO character | Tier 4 candidate (open test) |
| benzene HOMO at −β | F86b Dicke superposition γ = 1/2 anchor | Tier 4 candidate |

---

## Open questions (this folder, next sessions)

1. **Does the benzene-on-vibrational-bath Liouvillian satisfy F1 bit-exactly?**
   **Answered 2026-05-22** ([BENZENE_LIOUVILLIAN_PALINDROME](BENZENE_LIOUVILLIAN_PALINDROME.md)):
   conditionally yes. Under Holstein on-site dephasing (phonon coupled to the
   π-density, which equals the framework's Z-dephasing since D[n_l] = ¼·D[Z_l]) the
   open-system F1 palindrome holds, verified on the C₄ and C₆ rings; under
   Peierls/SSH bond dephasing it breaks. The Holstein case is the first direct F1
   test on a carbon substrate.

2. **Is there a benzene analog of the F98 (N+2)/[4(N+1)] → 1/4 long-time bridge?**
   F98 was derived for any truly-class Hamiltonian + Z-deph on N qubits; the bond
   topology drops out. If benzene's Liouvillian inherits the same truly-class
   structure (which it should, since Hückel ring is alternant), F98 should hold with N = 6
   giving α(∞) = 8/28 = 2/7.
   **Answered 2026-05-22** ([BENZENE_F98_LONG_TIME](BENZENE_F98_LONG_TIME.md)): yes,
   bit-exact. The KIntermediate Dicke state on the benzene XX+YY ring under Holstein
   dephasing reaches α(∞) = 2/7 at N = 6 (and 3/10 for the C₄ ring), confirmed as
   the exact t → ∞ limit via projection onto ker L.

3. **Klein-4-group character of HOMO at palindrome centre as aromaticity criterion?**
   Apply `KleinFourGroupSelfPairedRefinement` to the benzene + cyclobutadiene
   Liouvillians; check the Klein character of the HOMO at α; predict + verify the
   4n+2 vs 4n distinction.

4. **Does the V-Effect 14/19/3 trichotomy show up in benzene's electronic spectrum?**
   The framework's V-Effect splits N=3 Pauli-pair bilinears 14 truly + 19 soft + 3
   hard. For N=6 (benzene): what's the prediction, and does it match benzene's
   known electronic-transition pattern?

5. **The Frost-circle construction (1953) and the framework's mode dispersion.**
   Frost circle is the geometric mnemonic that Hückel ring eigenvalues are projections
   of inscribed N-gon vertices on a circle of radius 2β at centre α. The framework's
   `XyJordanWignerModes` uses the same cos(πk/(N+1)) dispersion at chain BC.
   Frost-cyclic and framework-OBC are two BC choices on the same algebraic structure.
   **Answered 2026-05-30** ([FROST_CIRCLE_AS_THE_CLOCK_FACE](FROST_CIRCLE_AS_THE_CLOCK_FACE.md)):
   one circle at two depths. The static Frost circle is the closed-system snapshot; the
   open-system clock runs it, adding what the still picture cannot hold, the band-edge
   π-coherence lifetime τ = 1/(2γ) and a coherent↔incoherent crossover Q* = J/γ (√2 at
   N=3, growing with chain length). Benzene's longest-lived π-coherence beats at 2|β|, the
   Frost radius itself; the open polyene chains beat at 2|β|·cos(π/(N+1)), the top π-MO.

---

## Why this might matter beyond "nice analogy"

The Coulson-Rushbrooke ↔ F1 identification is structurally clean (both are bipartite-
graph Z₂ palindromes), but more importantly it puts the framework's F-results on a
KNOWN HIGH-DATA SUBSTRATE. Carbon chemistry has 86 years of organic-chemistry data
verifying Coulson-Rushbrooke across thousands of molecules. The framework's F1 has 6
months of QPU + computational verification.

Every alternant hydrocarbon whose MO spectrum satisfies Coulson-Rushbrooke is
empirical evidence that the framework's Z₂-palindrome mechanism inherits to the
chemistry level. Conversely, every framework F-result that has a Hückel-equivalent
prediction (e.g. F86b's 3/8 Dicke anchor → benzene's half-filled p-shell) is
testable against decades of organic-chemistry data, both validating the inheritance
and surfacing new framework-derived chemistry predictions that weren't asked.

The carbon domain is the structurally-cleanest LARGE-DATA substrate for the framework.
Hydrogen-bond water (the existing `docs/water/`) is the cleanest SMALL-DATA substrate.
Together they cover the substrate spectrum.

---

## Anchor

- Script: [`simulations/carbon/benzene_huckel_framework_lens.py`](../../simulations/carbon/benzene_huckel_framework_lens.py)
- Parent docs: [`docs/carbon/README.md`](README.md), [`docs/water/README.md`](../water/README.md)
  (sister substrate)
- Framework anchors: [F1 palindrome](../ANALYTICAL_FORMULAS.md#f1), [F86b 3/8 Dicke
  anchor](../ANALYTICAL_FORMULAS.md#f86), [F87 trichotomy](../ANALYTICAL_FORMULAS.md#f87),
  [F98 (N+2)/[4(N+1)] bridge](../ANALYTICAL_FORMULAS.md#f98),
  `compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedRefinement.cs`
- Literature: Coulson + Rushbrooke (1940) "Note on the method of molecular orbitals",
  Proc. Camb. Phil. Soc. 36, 193; Frost (1953) "Frost circle"; Hückel (1931).
