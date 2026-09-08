# Random Matrix Theory Analysis of the Palindromic Liouvillian

<!-- Keywords: random matrix theory Liouvillian, spacing ratio palindromic
spectrum, Poisson level statistics open quantum, Liouvillian
spectral statistics, shifted-generator sectorwise P symmetry, XY-weight
light-content band universality class, Heisenberg dephasing RMT analysis,
R=CPsi2 random matrix theory -->

**Status:** Finite-size computational study (N=2-7, 21,840 eigenvalues, Heisenberg chain)
**Date:** April 1, 2026; sector-resolved comparison June 30, 2026; band-hint analysis June 27, 2026.
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** compute/RCPsiSquared.Compute (C# eigenvalue export),
[`simulations/rmt_analysis.py`](../simulations/rmt_analysis.py) (Python spacing analysis),
[`simulations/rmt_goe_hint_verdict.py`](../simulations/rmt_goe_hint_verdict.py) (raw band-frequency multiset diagnostics)

---

## What this document is about

This document reports direct consecutive-gap ratios of the raw full decay-rate
multiset, including zero gaps. The mean uses defined ratios only and has no
standard Poisson/GOE/GUE calibration for this degenerate unresolved population.
It neither establishes level clustering inside irreducible blocks nor proves integrability.
The companion frequency-band diagnostic preserves the selected multiset too;
neither analysis supplies a within-band GOE verdict (Result 3).

---

## Abstract

The C# exports contain all 4^N eigenvalues at N=2 through N=7. Direct
consecutive-gap ratios retain each original adjacent pair, including zero gaps;
0/0 is undefined and counted separately. At N=7 the defined-ratio mean is
0.2021456120489688 over 16,366 ratios, with 16 undefined pairs. These parsed-double
statistics are sensitive to near-degenerate eigensolver splittings and do not
certify physical level resolution or a universality class. The independent
**joint-popcount grading** is exact: the Hamiltonian conserves
excitation number and Z-dephasing does not move it either, so L is
block-diagonal across the (N+1)² blocks indexed by the popcounts of
ρ's row and column, and eigenvalues from different blocks never see
each other. The palindrome is a separate fact and not this one: it
pairs λ with −λ−2Σγ across the whole spectrum rather than splitting it
into non-interacting parts. The F1 palindrome (the centered spectrum has
exact ± pairing) is confirmed to machine precision. The separate raw
absolute-frequency multisets in average-light bands do not establish a
within-band onset of GOE-like repulsion or an integrability theorem.

---

## Background

### Why RMT?

Random-matrix universality tests require a specified spectral population and
appropriate symmetry resolution. The pooled real parts here include exact
numeric ties, near-degeneracies, and multiple invariant blocks. Standard
nondegenerate ensemble means are not a calibration for this population.

The adjacent spacing ratio r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1}) is a
standard diagnostic. It is exactly invariant under a common affine rescaling
of the spectrum and is comparatively insensitive to a slowly varying local
mean density because neighbouring gaps share approximately the same local
scale. It is not invariant under an arbitrary smooth nonlinear transformation,
and one mean value does not by itself classify a mixed or unresolved spectrum.

### What we knew before this analysis

The palindromic Liouvillian has:
- Exact eigenvalue pairing: λ + λ' = −2Σγ
  ([Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md))
- Light-content structure: the decay rate is 2γ·⟨n_XY⟩ exactly (the
  Absorption Theorem), so a mode of PURE XY-weight w sits at 2wγ
  ([F3](../docs/ANALYTICAL_FORMULAS.md)) and a mode of mixed content
  sits at its average, which need not be an integer
- Weight-parity mixing: the Hamiltonian mixes sectors w with w±2
  (discovered during [derivation verification](../docs/proofs/derivations/D05_DYNAMIC_MODE_COUNT.md))
- A symmetry-classification question distinct from random-matrix universality
  ([KMS analysis](../docs/KMS_DETAILED_BALANCE.md))

What we did NOT know: whether the spectrum shows level repulsion
(a chaotic-statistics signature) or level clustering.

---

## Method

### Eigenvalue computation (C# engine)

Heisenberg chain, J=1.0, uniform Z-dephasing γ=0.05, open
boundaries. All complex Liouvillian eigenvalues computed and
exported as CSV.

| N | Matrix | Eigenvalues | Engine | Time |
|---|--------|-------------|--------|------|
| 2 | 16x16 | 16 | MKL z_eigen | <1s |
| 3 | 64x64 | 64 | MKL z_eigen | <1s |
| 4 | 256x256 | 256 | MKL z_eigen | <1s |
| 5 | 1024x1024 | 1,024 | MKL z_eigen | 2s |
| 6 | 4096x4096 | 4,096 | MKL z_eigen | 58s |
| 7 | 16384x16384 | 16,384 | MKL z_eigen | 95min |

Command: `dotnet run -c Release -- rmt`
(in compute/RCPsiSquared.Compute/)

### Spacing ratio analysis (Python)

For decay rates d = −Re(λ):

1. Sort the complete parsed-double multiset: all 4^N rates, including zero and
   small negative numerical rates; no cutoff or deduplication.
2. For every original adjacent gap pair compute min/max. A zero next to a positive
   gap gives zero; 0/0 is undefined and remains in its original ratio position.
3. Average only defined ratios. Report level, gap, exact-zero-gap, adjacent-pair,
   defined-ratio and undefined-0/0 counts. Nearzero gaps are retained, not certified
   as physically resolved. There is no standard Poisson/GOE/GUE calibration.

No explicit spectral unfolding was applied in this analysis. Adjacent-gap
ratios cancel a common local scale exactly, hence are affine-invariant, and
are only approximately robust when the mean density varies slowly across two
neighbouring gaps. They are not invariant under general smooth reparameterizations,
so this property does not remove the unresolved-population limitation.

### Band analysis

For N=2-7 (the same C# CSVs), select rate windows |d−2wγ| < 0.3γ,
then analyze the complete absolute-frequency multiset in each window.
Rates within a finite window need not be equal.

The bin is a **band** around average light content ⟨n_XY⟩ = w, not a
weight sector. Modes of pure weight w land in it, and so does any mixed
mode within the tolerance, including mixed modes exactly on the rung.
The distinction affects which reference population could be justified;
these bands are not a symmetry-resolved universality test.

---

## Result 1: Direct Raw-Multiset Consecutive-Gap Ratios

Executed with `python simulations/rmt_analysis.py` from the committed CSVs.
Here L is the level count: there are L−1 gaps and L−2 original adjacent pairs.
The last two count columns sum to L−2; zero gaps count exact parsed-double ties.

| N | L | Zero gaps | Defined ratios | Undefined 0/0 | ⟨r⟩ defined |
|---|---:|---:|---:|---:|---:|
| 2 | 16 | 4 | 14 | 0 | 0.09635955831608055 |
| 3 | 64 | 7 | 62 | 0 | 0.22260363264213059 |
| 4 | 256 | 35 | 244 | 10 | 0.27789182063152246 |
| 5 | 1024 | 14 | 1021 | 1 | 0.21627032685910053 |
| 6 | 4096 | 1340 | 3107 | 987 | 0.18988807786718778 |
| 7 | 16384 | 136 | 16366 | 16 | 0.2021456120489688 |

The optional subset d < Nγ has no lower cutoff. It is not an irreducible sector
and does not certify removal of correlations or numerical degeneracies.

| N | Subset L | Zero gaps | Defined ratios | Undefined 0/0 | ⟨r⟩ defined |
|---|---:|---:|---:|---:|---:|
| 2 | 9 | 3 | 7 | 0 | 0.03571428571428635 |
| 3 | 32 | 3 | 30 | 0 | 0.19280201821404228 |
| 4 | 114 | 9 | 111 | 1 | 0.27863805704900285 |
| 5 | 512 | 6 | 509 | 1 | 0.21681497117520818 |
| 6 | 1973 | 608 | 1540 | 431 | 0.18227295486195011 |
| 7 | 8192 | 65 | 8184 | 6 | 0.20163224832233276 |

These are descriptions of the numeric multiset, not estimates of a resolved
nondegenerate spacing process. They do not establish an irreducible-sector
repulsion law, integrability, or a large-N limit.

### Why the pooled population does not classify a sector

ONE grading cuts the Liouvillian into pieces that never interact: the
**joint popcount**. H conserves excitation number and Z-dephasing acts
diagonally on ρ's indices, so L is block-diagonal across the (N+1)²
blocks indexed by (popcount(row), popcount(col)). Levels in different
blocks need not repel each other. Mixing these sectors can suppress
repulsion in the pooled statistic; block diagonality alone neither forces
Poisson statistics nor proves integrability inside a block.

XY-weight parity does not add a second cut, though it is easy to think
it does. (−1)^n_XY does commute with L, but it is not independent of
the grading above: inside the block (p, q) every Pauli component has
n_XY ≡ p − q (mod 2), so the parity is a function of the block and
splits nothing. Checked at N=3, 4, 5: no joint-popcount block carries
both parities.

The palindrome is not one of these. Π conjugation sends λ to −λ−2Σγ
and so relates the two halves of the spectrum; it pairs levels rather
than separating them, and a pairing does not by itself suppress
repulsion. The d < Nγ restriction is only a numerical subset diagnostic.

The spacing statistics are consistent with what the analytical
formulas already allowed. They do not establish that the spectrum is
fully determined: only the (0,1) coherence block has a closed form
(F2), and the interior blocks do not. What is measured is a finite-N pooled
mean ratio, which is much weaker than an unfolded irreducible-sector point
process or knowledge of the levels.

---

## Result 2: Centered Decay-Rate Reflection Check

The exact F1 theorem says that the centered Liouvillian L_c = L + Σγ·I
has ± pairing of the full complex eigenvalue multiset. The diagnostic here
checks only the induced reflection of the centered decay rates
d_c = -Re(λ_c). It performs a multiplicity-preserving perfect matching of the
full real-part projection, including the modes at the centre; on the real line
this is the bottleneck distance between sorted d_c and sorted -d_c. It does not
match the imaginary parts or independently re-establish the full complex
multiset bijection.

| N | Positive | Negative | Central | Multiplicity bottleneck error |
|---|----------|----------|---------|-------------------------------|
| 2 | 3 | 3 | 10 | 2.78e-16 |
| 3 | 32 | 32 | 0 | 4.11e-15 |
| 4 | 52 | 52 | 152 | 2.38e-14 |
| 5 | 512 | 512 | 0 | 1.64e-14 |
| 6 | 1,130 | 1,130 | 1,836 | 3.10e-14 |
| 7 | 8,192 | 8,192 | 0 | 6.20e-14 |

The projected decay-rate reflection is satisfied to machine precision at every
tested N and is consistent with the independently proved full palindrome
([Mirror Symmetry](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)). This projected
check alone neither confirms the imaginary-part pairing nor proves the
operator identity. It does not assign
a random-matrix class by itself. In the Sá-Ribeiro-Prosen construction,
Hermiticity preservation also supplies T₊, and the global order-4 Π must be
restricted and phase-normalized inside the Π² parity sectors before the
remaining symmetry algebra is classified. That irreducible-sector calculation
has not been completed here. The raw-multiset finite-size statistic and the
exact spectral pairing are therefore separate observations; neither proves
integrability or licenses a class label.

---

## Result 3: Raw Absolute-Frequency Multisets in Average-Light Bands

The companion [producer](../simulations/rmt_goe_hint_verdict.py) uses the same
C# CSVs and direct adjacent-gap rule. Its legacy filename does not name a
current ensemble verdict. It selects |d−2wγ| < 0.3γ, an average-light window
|⟨n_XY⟩−w| < 0.15, and sorts every selected |Im λ|, including zero frequencies
and multiplicities. There is no rounding, deduplication, or frequency cutoff.
These are not invariant fixed-XY-weight sectors, nor independent nondegenerate
frequency samples. No standard-ensemble calibration applies.

Executed with `python simulations/rmt_goe_hint_verdict.py`;
[all band populations](../simulations/results/rmt_band_multiset.txt) are printed.
Representative central bands:

| N | w | L | Zero gaps | Defined ratios | Undefined 0/0 | ⟨r⟩ defined |
|---|---|---:|---:|---:|---:|---:|
| 4 | 2 | 152 | 1 | 150 | 0 | 0.18169077812422024 |
| 5 | 2 | 80 | 1 | 78 | 0 | 0.24815201187882771 |
| 5 | 3 | 80 | 0 | 78 | 0 | 0.27828836818199409 |
| 6 | 3 | 1948 | 14 | 1946 | 0 | 0.1659385190868666 |
| 7 | 3 | 1710 | 17 | 1708 | 0 | 0.19686572210722628 |
| 7 | 4 | 1710 | 7 | 1708 | 0 | 0.21256889000702267 |

F43 pairs exact reflected bands with a multiplicity-preserving frequency-sign
bijection. It does not make raw floating-point tiny gaps reliable: reflected
bands can have different numeric mean ratios despite the exact identity.
The endpoint bands likewise contain tiny numerical frequencies, not a
certified nonzero-frequency population. This sensitivity is a limitation of
these descriptive summaries, not physical breaking of the palindrome.
No within-band GOE repulsion or integrability verdict is established.

---

## Result 4: All Eigenvalues in the Left Half-Plane

Every nonzero eigenvalue has Re(λ) < 0; the fraction with Re < 0 is
1.0000 at every N tested. This is a finite-size spectral-stability
check only. The location of the eigenvalues does not establish complete
positivity or trace preservation. The Lindblad construction supplies
those properties independently; this spectrum merely shows no growing
mode in the sampled generators.

---

## Result 5: Finite-Size Filling-Associated Crossover Evidence (June 2026)

The sector-resolved comparison gives finite executed CSR evidence for a
**filling dependence** at canonical Delta=1 plus disorder. Working with
the complex spacing ratio (CSR, Sá-Ribeiro-Prosen) on coherence blocks
(wKet, wBra) of the Z-dephased XXZ Liouvillian:

This is filling-associated crossover evidence at N=6..8: movement toward GinUE,
not a causal or thermodynamic threshold theorem.
Each spectrum is first reduced to one 1e-9 finite-precision cluster representative per rounded
coordinate pair. The resulting counts and CSR values are tolerance-dependent and not an exact degeneracy census;
sufficiently close nondegenerate levels can merge at this resolution.

- The **dilute** (SE,DE) = (1,2) block, the Door-C block, where the non-solvable
  Galois group S_d lives, stays Poisson-like / non-GinUE over the sampled
  anisotropy and random-field sweeps (`inspect --root galoischaos`, the Δ=0
  control; the two Door-C sweep stages).
- A **dense** block (p, p+1) near half-filling at the same model parameters and disorder
  distribution, but from a separately sampled realization ensemble,
  the **same** disorder + interactions, **moves toward GinUE**: its radial CSR ⟨|z|⟩ is
  near the GinUE reference and its angular repulsion ⟨cos θ⟩ goes negative and climbs
  toward GinUE with the block size (≈ −0.09 → −0.13 → −0.16 at N = 6/7/8 = 43% →
  56% → 67% of the size-matched GinUE angle), while the dilute block stays flat at
  ⟨cos θ⟩ ≈ 0 (~14–23%).

The knobs have three distinct Hamiltonian meanings: nonzero Delta breaks free-fermion
additivity, but uniform XXZ remains Bethe-integrable;
at Delta=0 the random-field XY Hamiltonian remains quadratic (Anderson/free fermions),
which does not classify the Z-dephasing Liouvillian as a quadratic generator;
generic random field plus Delta!=0 is the interacting disordered nonintegrable test.
Reflection/conjugation/cross-fold breaking need not break Hamiltonian integrability.
The canonical Delta=1 plus disorder comparison supports a finite-size filling dependence,
not a universal threshold or proof of thermalization. Galois structure over the coupling
and spectral statistics at fixed coupling are distinct. GinUE is used here only
as the class-A comparison ensemble. Unequal weight (p,p+1) sends Π to the conjugate
(p+1,p) block rather than furnishing an internal symmetry; the near-zero disordered
conjugation-match fraction excludes that sampled pairing but does not prove that no
other antiunitary symmetry survives. The full irreducible class remains **OPEN**
until the shifted generator's sectorwise P symmetry and all unitary/strong sectors
have been resolved. Live:
`inspect --root fillcsr` (`FillingThresholdWitness`); full writeup in
[FILLING_THRESHOLD_CHAOS.md](FILLING_THRESHOLD_CHAOS.md).

## What This Does Not Answer

Result 3 reports the tested within-band statistics, not their large-N limit.
Result 5 settles only the executed N=6..8 fixed-parameter comparison:
dense coherence blocks move toward GinUE more strongly than dilute blocks
at canonical Delta=1 plus disorder. It does not settle convergence,
thermalization, or a thermodynamic threshold. Two open items remain.

1. **Comparison with Denisov lemon shape.** The complex-plane density
   of random Lindbladians (Denisov et al., PRL 2019) has a specific
   "lemon" shape (the characteristic boundary curve of eigenvalue density for structureless random Lindbladians). Our palindromic constraint modifies this. A
   quantitative comparison needs the 2D density, not just 1D rates.

2. **Topological dependence.** All results above use chain topology. Star,
   ring, and complete topologies were since surveyed with the complex
   spacing ratio in `simulations/rmt_topology_csr.py` (chain reads clean
   2D-Poisson; at the producer's declared 1e-9 rounding, the symmetric
   topologies fragment into too few finite-precision cluster representatives
   for global non-Hermitian RMT). These cluster counts are tolerance-dependent,
   not exact degeneracy counts. A full sector-resolved comparison across
   topologies is still open.

---

## Connection to the Framework

The pooled mean-ratio result is consistent with the picture the
[analytical formulas](../docs/ANALYTICAL_FORMULAS.md) draw, and it is
worth stating the implication in the direction it actually runs.
The measurement does not run backward to integrability or a universality
class. In fact the repository has a
closed form for ONE block, the (0,1) coherence block (F2); the interior
blocks have none.

Joint popcount is one concrete reason the pooled population is not an
irreducible-level statistic: it cuts the Liouvillian into blocks that never
interact, and levels from different blocks need not repel. The present pooled
analysis does not compare against a symmetry-resolved control, so it does not
isolate block mixing as the cause of the measured ratios. The palindromic
constraint (F1) pairs the two halves of the spectrum instead of defining this
block decomposition, but this run likewise does not measure whether that
pairing changes the pooled statistic.

The exact P-type anticommutation and raw-multiset spacing ratios are different
statements. The spacing statistic does not assign the full irreducible-sector
symmetry class and does not identify Hamiltonian or Liouvillian integrability.

The raw non-unfolded frequency spectral form factor (SFF) shows sampled modulation,
compatible with structured spectra but not proof of integrability. It reports
sampled modulation near ω_min = 4J(1-cos(π/N)) and
palindrome-paired decay-rate bands in the time domain. Under uniform dephasing
their centres can be labeled by average light w ↔ N-w, but they are not
invariant fixed-XY-weight eigenvalue sectors. The spacing ratio describes local adjacent-gap correlations in the chosen ordering;
the SFF describes global spectral structure. These diagnostics cannot classify universality
from the raw mixed spectrum. The SFF's multiplicity-dependent raw multiset density
scale defines descriptive bins, not physical time regimes. See [Spectral Form Factor](SPECTRAL_FORM_FACTOR.md).

---

## Reproducibility

| Component | Location |
|-----------|----------|
| C# eigenvalue export | compute/RCPsiSquared.Compute/ (`dotnet run -c Release -- rmt`) |
| Python RMT analysis | [`simulations/rmt_analysis.py`](../simulations/rmt_analysis.py) |
| Eigenvalue CSVs | `simulations/results/rmt_eigenvalues_N{2..7}.csv` |
| Analysis output | [`simulations/results/rmt_analysis.txt`](../simulations/results/rmt_analysis.txt) |

---

## References

- Sa, L., Ribeiro, P., Prosen, T. (2023). "Symmetry Classification
  of Many-Body Lindbladians." PRX 13, 031019.
  (38 symmetry classes of Lindbladians)
- Denisov, S. et al. (2019). "Universal Spectra of Random Lindblad
  Operators." PRL 123, 140403. (Lemon-shaped spectral density)
- Oganesyan, V., Huse, D.A. (2007). "Localization of interacting
  fermions at high temperature." PRB 75, 155111.
  (Spacing ratio diagnostic, original paper)
- Atas, Y.Y. et al. (2013). "Distribution of the Ratio of Consecutive
  Level Spacings." PRL 110, 084101. (Spacing ratio reference values)
- Verbaarschot, J. (1994). "Spectrum of the QCD Dirac operator and
  chiral random matrix theory." PRL 72, 2531.
  (chiral RMT reference; not a class assignment for these pooled Liouvillian data)

---

*The 21,840 eigenvalues supply finite-N pooled spacing-ratio data, not an integrability proof.
Joint-popcount grading separates invariant blocks; the palindrome pairs levels
rather than separating them. One block has a closed form, while the measured
mean ratios do not determine the interior spectra or a universality class.*
