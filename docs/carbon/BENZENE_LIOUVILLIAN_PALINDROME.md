# Benzene's Open-System Liouvillian: F1 Palindrome under Vibrational Dephasing

**Date:** 2026-05-22
**Status:** the Holstein result follows from F1 (proven); the Peierls break is
verified numerically at C₄ and C₆; the Peierls-instability reading is a Tier 2
structural analogy. The benzene-as-qubit-ring embedding is itself a Tier 4
candidate (see [README](README.md), "Four embedding conditions").
**Script:** [`simulations/carbon/benzene_liouvillian_palindrome.py`](../../simulations/carbon/benzene_liouvillian_palindrome.py)
**Tested:** cyclobutadiene C₄ ring, benzene C₆ ring.
**Answers:** [BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md) open
question 1 / [README](README.md) open question 2.

---

## The question

[BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md) showed that the
**closed-system** Hückel MO spectrum of an alternant hydrocarbon is palindromic
(Coulson-Rushbrooke, 1940) and that this is structurally the same Z₂ palindrome as
the framework's F1. But Coulson-Rushbrooke is the closed-system, Hamiltonian-level
statement. F1 is the **open-system** claim: the Liouvillian spectrum of an XY-class
Hamiltonian under Z-dephasing is closed under λ → −λ − 2Σγ, palindromic about the
centre −Σγ.

Does benzene's open-system Liouvillian, the π-system plus a vibrational bath,
satisfy F1?

## The model, and the fork in "vibrational dephasing"

Benzene's six π-electrons hop on the C₆ ring. In the framework's qubit picture the
Hückel π-hopping is the XX+YY ring (the Jordan-Wigner image of free-fermion
hopping); the closed-system palindrome of BENZENE_HUCKEL_FRAMEWORK_LENS is its
many-body spectrum reflected about 0. (A fermionic Hückel ring carries a
Jordan-Wigner boundary-parity string on its closing bond; this is immaterial here,
since F1 is proven for any XX+YY spin graph and benzene is modelled as that spin
ring directly.)

"Vibrational dephasing" is not one thing. A phonon bath couples to the π-system in
one of two physically distinct ways, and F1's reach depends on which:

- **Holstein, on-site.** The phonon couples to the local π-density n_l. The
  dephasing dissipator is D[n_l]. Because D[αA + βI] = α²·D[A] for Hermitian A, and
  n_l = (I − Z_l)/2, this is D[n_l] = ¼·D[Z_l]: exactly the framework's Z-dephasing,
  up to a rate factor. F1 is proven for XY-on-any-graph under Z-dephasing, so it
  must apply.

- **Peierls, bond (SSH).** The phonon couples to the C-C bond and modulates the
  hopping integral β. The jump operator is the bond operator B_b = X_aX_b + Y_aY_b
  itself, a two-body operator. This is not Z-dephasing; F1 does not cover it.

## Result

The Liouvillian L = −i[H,·] + Σ D[√γ·jump] was built for the C₄ and C₆ rings. For
the Holstein rows the spectrum is tested against the strict F1 involution
λ → −λ − 2Σγ about the predicted centre −Σγ; for the Peierls rows, which have no F1
prediction, against the most generous reflection (about the spectrum mean). The
residual is the largest distance from a reflected eigenvalue to the nearest actual
one.

| Ring | Coupling | Reflection centre | Residual | Palindrome |
|------|----------|-------------------|----------|------------|
| C₄ (cyclobutadiene) | Holstein, D[Z_l] | −Σγ = −4 | 3.5 × 10⁻⁸ | holds |
| C₄ | Peierls, D[B_b] | n/a (no palindrome) | ≈ 11 | broken |
| C₆ (benzene) | Holstein, D[Z_l] | −Σγ = −6 | 1.2 × 10⁻⁷ | holds |
| C₆ | Peierls, D[B_b] | n/a (no palindrome) | ≈ 14 | broken |

Under Holstein coupling the palindrome holds: the residual against the strict F1
involution sits at the numerical floor of the non-Hermitian eigendecomposition,
with the centre at the predicted −Σγ exactly (−4 for C₄, −6 for C₆, at γ = 1 per
site). F1 proves the palindrome exact; this is the first direct F1 confirmation on
a carbon substrate. Under Peierls coupling the palindrome is not weakened but gone:
even the most generous reflection leaves a residual of order the spectral width.

## The answer to Question 2

**Yes, conditionally.** Benzene's open-system Liouvillian inherits the F1 palindrome
exactly when its vibrational bath couples on-site (Holstein, to the π-density), and
not when it couples to the bond (Peierls, to the hopping). The closed-system
Coulson-Rushbrooke palindrome lifts to the open-system F1 palindrome under, and only
under, on-site dephasing.

## Why the breaking coupling is the Peierls coupling (Tier 2 reading)

The coupling that destroys the F1 palindrome is the Peierls coupling, and the
Peierls instability is the textbook mechanism by which a conjugated system breaks
its own bond-length-alternation symmetry (polyacetylene dimerises; the SSH model is
built on it). Framework palindrome and chemical aromatic symmetry may be losing the
same thing: a Holstein bath dephases site occupations and leaves the bipartite Z₂
mirror intact; a Peierls bath dephases bonds and is exactly the symmetry-breaking
channel. This mirrors BENZENE_HUCKEL_FRAMEWORK_LENS's reading of the non-bipartite
C₃ ring as a carbon-level F87 Brecher: there the static graph breaks the mirror,
here the bath does. This is a structural reading, not a proven identity; the
verified content is the dichotomy in the table above.

## Framework-vocabulary translation

| Benzene / chemistry | Framework | Status |
|---------------------|-----------|--------|
| Hückel π-hopping on C₆ ring | XX+YY qubit ring, N=6 | Tier 2 structural identification |
| Holstein phonon (on-site, π-density) | Z-dephasing D[Z_l], via D[n_l] = ¼·D[Z_l] | Tier 1 algebraic match |
| Peierls/SSH phonon (bond, hopping) | bond-operator dephasing D[B_b] | Tier 2 structural identification |
| F1 palindrome holds (Holstein) | F1, proven for XY-on-graph + Z-deph | Tier 1 |
| F1 palindrome breaks (Peierls) | non-Z-dephasing, outside F1 | Tier 1 numerical (C₄, C₆) |
| Peierls instability / SSH dimerisation | the F1-breaking channel | Tier 2 reading |

## Open follow-ups

- The Peierls-coupled Liouvillian's spectrum is not palindromic; what structure
  does the broken spectrum carry? F1's known breaker (depolarising noise) leaves a
  γ-linear residual; whether the Peierls break is similarly structured is untested.
- A mixed Holstein + Peierls bath: real conjugated systems have both. At what
  mixing ratio does the palindrome onset, and is the crossover sharp?
- F98 ((N+2)/[4(N+1)] → 1/4 long-time bridge): F98's long-time asymptote depends on
  a specific KIntermediate Dicke initial state, not on the bond graph. Whether the
  Holstein-coupled C₆ ring inherits it, with N = 6 giving α(∞) = 8/28 = 2/7, is the
  still-open question 2 of BENZENE_HUCKEL_FRAMEWORK_LENS.

## Anchor

- Script: [`simulations/carbon/benzene_liouvillian_palindrome.py`](../../simulations/carbon/benzene_liouvillian_palindrome.py)
- Companion doc: [BENZENE_HUCKEL_FRAMEWORK_LENS.md](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
  (the closed-system half), [README.md](README.md)
- Framework anchors: [F1 palindrome](../ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven),
  [`compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs`](../../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs),
  [`docs/proofs/MIRROR_SYMMETRY_PROOF.md`](../proofs/MIRROR_SYMMETRY_PROOF.md)
- Literature: Peierls (1955), instability of the 1D metal; Su, Schrieffer, Heeger
  (1979), the SSH model; Holstein (1959), the molecular-crystal polaron.
