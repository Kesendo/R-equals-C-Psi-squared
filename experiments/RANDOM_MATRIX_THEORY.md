# Random Matrix Theory Analysis of the Palindromic Liouvillian

<!-- Keywords: random matrix theory Liouvillian, spacing ratio palindromic
spectrum, Poisson level statistics open quantum, integrable Lindbladian
spectral statistics, chiral symmetry Liouvillian eigenvalues, XY-weight
light-content band universality class, Heisenberg dephasing RMT analysis,
R=CPsi2 random matrix theory -->

**Status:** Computationally verified (N=2-7, 21,840 eigenvalues, Heisenberg chain)
**Date:** April 1, 2026 (updated June 30, 2026: dissipative chaos located as a FILLING threshold, Result 5; the within-band GOE hint resolved as a small-sample artifact June 27, Result 3)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** compute/RCPsiSquared.Compute (C# eigenvalue export),
[`simulations/rmt_analysis.py`](../simulations/rmt_analysis.py) (Python spacing analysis),
[`simulations/rmt_goe_hint_verdict.py`](../simulations/rmt_goe_hint_verdict.py) (band GOE-hint verdict: bootstrap + larger-N)

---

## What this document is about

This document asks whether the palindromic Liouvillian spectrum looks
chaotic or integrable, using the standard diagnostic from random matrix
theory: the spacing ratio between consecutive eigenvalues. The answer:
the spectrum is Poisson (integrable, eigenvalues cluster rather than
repel), which is what a Liouvillian cut into non-interacting blocks by
its conserved quantities should look like. An early small-N read hinted
at GOE-like repulsion inside individual decay-rate bands, but driving it
to a verdict (bootstrap + larger samples) shows it was small-sample
noise: the bands are integrable too, with no chaotic transition
(Result 3).

---

## Abstract

We perform the first random matrix theory (RMT) analysis of the
palindromic Liouvillian spectrum. Using spacing ratios (robust,
unfolding-free) on eigenvalues computed by the C# engine with
MKL/OpenBLAS (N=2 through N=7, up to 16,384 eigenvalues), we find
that the decay rate spectrum is **Poisson** (integrable, no level
repulsion) at every system size tested. The mean spacing ratio
converges to ⟨r⟩ = 0.36-0.39, consistent with the Poisson value
0.386 and far from GOE (0.536) or GUE (0.603). The decomposition that
does this is the **joint-popcount grading**: the Hamiltonian conserves
excitation number and Z-dephasing does not move it either, so L is
block-diagonal across the (N+1)² blocks indexed by the popcounts of
ρ's row and column, and eigenvalues from different blocks never see
each other. The palindrome is a separate fact and not this one: it
pairs λ with −λ−2Σγ across the whole spectrum rather than splitting it
into non-interacting parts. The chiral symmetry (centered spectrum has
exact ± pairing) is confirmed to machine precision. An early
band-resolved read at N=5 showed an apparent GOE-like ⟨r⟩ = 0.513
inside individual decay-rate bands; driven to a verdict it is a
small-sample artifact (the n=15 value is a 1.5σ Poisson fluctuation, and
the same bands read Poisson/sub-Poisson at N=6-7 with hundreds of
frequencies). There is no integrable-to-chaotic transition: the spectrum
is integrable in every band and at every tested N.

---

## Background

### Why RMT?

Random Matrix Theory classifies quantum spectra by their level
statistics. The three standard universality classes are:

| Ensemble | Spacing ratio ⟨r⟩ | Level repulsion | Physics |
|----------|-------------------|-----------------|---------|
| Poisson  | 0.386 | None (clustering) | Integrable, conserved quantities |
| GOE      | 0.536 | Linear (s^1) | Time-reversal invariant, real |
| GUE      | 0.603 | Quadratic (s^2) | Time-reversal broken, complex |

The spacing ratio r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1}) is the
modern standard diagnostic: it requires no spectral unfolding and gives
a single number that distinguishes the three classes.

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
- Does not fit any of the 38 Sa-Ribeiro-Prosen classes
  ([KMS analysis](../docs/KMS_DETAILED_BALANCE.md))

What we did NOT know: whether the spectrum shows level repulsion
(chaotic) or level clustering (integrable).

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

For decay rates (Re parts of eigenvalues):
1. Sort all nonzero rates
2. Compute consecutive spacing ratios r_n
3. Average ⟨r⟩ and compare with reference values

No spectral unfolding needed. The spacing ratio is invariant under
smooth transformations of the spectrum.

### Band analysis

For N=2-5 (Python eigenvalues): bin modes by nearest integer w
(rate ≈ 2·w·γ), then analyze FREQUENCIES (imaginary parts) within each
bin. Rates are degenerate inside a bin by construction; frequencies
carry the spectral information.

The bin is a **band** around average light content ⟨n_XY⟩ = w, not a
weight sector. Modes of pure weight w land in it, and so does any mixed
mode within the tolerance, including mixed modes exactly on the rung.
Everything below reads bands; the distinction does not affect the
statistics but it does decide what the sets are called.

---

## Result 1: The Spectrum is Poisson (Integrable)

| N | Eigenvalues | ⟨r⟩ (all rates) | ⟨r⟩ (lower half) | Classification |
|---|-------------|-----------------|-------------------|----------------|
| 3 | 64 | 0.220 | 0.301 | Poisson |
| 4 | 256 | 0.408 | 0.385 | Poisson |
| 5 | 1,024 | 0.369 | 0.368 | Poisson |
| 6 | 4,096 | 0.364 | 0.363 | Poisson |
| 7 | 16,384 | 0.383 | 0.383 | Poisson |

(N=2 has too few eigenvalues for statistics.)

**The palindromic Liouvillian is integrable.** The mean spacing ratio
converges to ~0.37, consistent with the Poisson value 0.386. There is
no level repulsion. Eigenvalues cluster rather than repel.

The "lower half" analysis (only rates below Σγ, removing the
palindromic pairing correlation) gives the same result. The Poisson
statistics are intrinsic to each half of the spectrum, not an artifact
of the palindromic doubling.

### Why Poisson?

ONE grading cuts the Liouvillian into pieces that never interact: the
**joint popcount**. H conserves excitation number and Z-dephasing acts
diagonally on ρ's indices, so L is block-diagonal across the (N+1)²
blocks indexed by (popcount(row), popcount(col)). Levels in different
blocks are unrelated by construction and have no reason to repel; that
is what produces Poisson, and it is the same reason integrable systems
are Poisson.

XY-weight parity does not add a second cut, though it is easy to think
it does. (−1)^n_XY does commute with L, but it is not independent of
the grading above: inside the block (p, q) every Pauli component has
n_XY ≡ p − q (mod 2), so the parity is a function of the block and
splits nothing. Checked at N=3, 4, 5: no joint-popcount block carries
both parities.

The palindrome is not one of these. Π conjugation sends λ to −λ−2Σγ
and so relates the two halves of the spectrum; it pairs levels rather
than separating them, and a pairing does not by itself suppress
repulsion. That the two halves are individually Poisson (measured
above) is evidence for exactly this reading.

The spacing statistics are consistent with what the analytical
formulas already implied. They do not establish that the spectrum is
fully determined: only the (0,1) coherence block has a closed form
(F2), and the interior blocks do not. What is Poisson is the level
STATISTICS, which is a much weaker statement than knowing the levels.

---

## Result 2: Perfect Chiral Symmetry

The centered Liouvillian L_c = L + Σγ·I has exact ±
eigenvalue pairing (the palindromic constraint becomes
λ_c + λ_c' = 0).

| N | ± pairs | Mean pairing error |
|---|-----------|-------------------|
| 3 | 32 | 8.3e-16 |
| 4 | 52 | 1.5e-15 |
| 5 | 512 | 2.8e-15 |
| 6 | 1,130 | 4.7e-15 |
| 7 | 8,192 | 6.6e-15 |

Pairing is exact to machine precision at every N. This confirms the
algebraic proof ([Mirror Symmetry](../docs/proofs/MIRROR_SYMMETRY_PROOF.md))
and places the centered Liouvillian in the chiral symmetry class.

In the Altland-Zirnbauer classification (the tenfold taxonomy of symmetry classes for random matrices, extending Wigner-Dyson's three classes to include particle-hole and chiral symmetries), this is class AIII (chiral
unitary). However, the Poisson level statistics show that the system
does not exhibit the level repulsion expected for a random chiral
GUE ensemble. The palindromic Liouvillian is chiral but integrable:
it has the symmetry of class AIII but the statistics of an integrable
system.

---

## Result 3: The Within-Band GOE Hint, Resolved (Artifact)

Within individual decay-rate bands, frequencies (not rates) were measured
with the spacing ratio ⟨r⟩. An early small-N read looked GOE-like:

| N | Band | Unique freq | ⟨r⟩ | Class |
|---|--------|-------------|-----|-------|
| 4 | w=2 | 41 | 0.130 | sub-Poisson |
| 5 | w=2 | 15 | 0.513 | GOE |
| 5 | w=3 | 15 | 0.513 | GOE |

**This hint is an artifact of small sample size.** It was driven to a verdict
in [`simulations/rmt_goe_hint_verdict.py`](../simulations/rmt_goe_hint_verdict.py)
(reproduce + bootstrap + extend to N=6-7):

1. **The two "GOE" rows are one sample, not two.** At N=5, w=2 and w=N−2=3 are
   palindromic partners with identical frequency content (the F43 sector pairing
   K_freq(w,t) = K_freq(N−w,t)), so they read identically (0.513 on the same 15
   frequencies). There was never an independent second band.

2. **0.513 on 15 frequencies is a Poisson fluctuation.** The Poisson sampling
   band for ⟨r⟩ at n=15 (Monte Carlo over a homogeneous Poisson process, same
   spacing-ratio statistic) is 0.386 ± 0.087, with [5%, 95%] = [0.245, 0.533].
   The observed 0.513 sits inside the band; the one-sided p(Poisson ≥ 0.513) =
   0.076, a ~1.5σ upward fluctuation, not significant.

3. **Larger samples read Poisson, not GOE.** The decay rate can be read off the
   full Liouvillian spectrum at N=6-7 with no eigenvectors and no extra memory
   (this is the obstacle the original write-up wrongly thought blocked the
   check), so the binning is cheap. What it selects needs saying precisely,
   because the earlier wording here got it wrong twice over. The Absorption
   Theorem gives rate = 2γ·⟨n_XY⟩, the AVERAGE light content, not the weight,
   so a rate bin does not select a weight sector. And the bin is a BAND, not a
   rung: `rmt_analysis.py` admits |rate − 2wγ| < 0.3γ, so it collects every
   mode whose ⟨n_XY⟩ falls within 0.3 of w, and most of those are not at w.
   So these are **bands** around ⟨n_XY⟩ = w. Sitting exactly ON the rung while
   being mixed is possible and does happen: at N=4, γ=0.05 the 4γ rung carries
   modes with histogram {1: ½, 3: ½} and ⟨n_XY⟩ = 2.000000 exactly, beside the
   pure weight-2 ones. The measurement is unaffected: a band is still a
   well-defined set of genuine eigenvalues of the full Liouvillian, and it is
   still large. N=6 band 3 (546 freq) → ⟨r⟩ = 0.272, N=7 bands 3 and 4
   (414 freq) → 0.283. The reading does not approach GOE (0.536); it converges
   to Poisson and below (sub-Poisson = level clustering, the opposite of
   repulsion, the signature of a strongly degenerate integrable additive
   spectrum, consistent with the N=4 band-2 row). The tiny bands 1 and N−1
   (5-6 frequencies) throw a spurious "GUE" ⟨r⟩ > 0.79, plainly small-sample
   noise.

**Verdict: no within-band chaos.** The system is integrable in every band and
at every tested N; the earlier GOE hint was small-sample noise. This matches
the global Poisson result above and the sector-resolved non-Hermitian test (the
`galoischaos` witness, `inspect --root galoischaos`), which independently reads
the Galois-S_n half of the (SE,DE) block Poisson-like / sub-Poisson, not
Ginibre.

---

## Result 4: All Eigenvalues in the Left Half-Plane

Every nonzero eigenvalue has Re(λ) < 0, confirming that the
Liouvillian is a proper generator of a completely positive trace-
preserving (CPTP) semigroup. The fraction of eigenvalues with
Re < 0 is 1.0000 at every N tested.

This is a consistency check, not a new result. But it confirms that
the C# eigenvalue export is producing physically valid spectra.

---

## Result 5: Dissipative Chaos is a Filling Threshold (June 2026)

The deeper "does any sector reach dissipative quantum chaos?" question is now
answered, and the answer is about **filling**, not integrability. Working with
the complex spacing ratio (CSR, Sá-Ribeiro-Prosen) on coherence blocks
(wKet, wBra) of the Z-dephased XXZ Liouvillian:

- The **dilute** (SE,DE) = (1,2) block, the Door-C block, where the non-solvable
  Galois group S_d lives, stays Poisson-like / non-GinUE under **every**
  integrability-breaking knob (XXZ Δ, a random Z-field, with or without
  interactions). That null is robust (`inspect --root galoischaos`, the Δ=0
  control; the two Door-C sweep stages).
- A **dense** block (p, p+1) near half-filling of the **same** Liouvillian, under
  the **same** disorder + interactions, **is chaotic**: its radial CSR ⟨|z|⟩ sits
  at the GinUE value and its angular repulsion ⟨cos θ⟩ goes negative and climbs
  toward GinUE with the block size (≈ −0.09 → −0.13 → −0.16 at N = 6/7/8 = 43% →
  56% → 67% of the size-matched GinUE angle), while the dilute block stays flat at
  ⟨cos θ⟩ ≈ 0 (~14–23%).

So fixed-q dissipative chaos switches on with **extensive excitation content**,
not with breaking the Galois/Hamiltonian integrability. Galois chaos (over the
coupling q) and spectral chaos (GinUE at fixed q) merge only at extensive filling;
the dilute (SE,DE) sector that carries S_d is too dilute to thermalize, and its
persistent Poisson statistics are the kinematic shadow of that. Class A is licensed
by the unequal weight (p,p+1) (Π maps it to the conjugate (p+1,p) block, not itself;
the disordered conjugation-match fraction is ≈ 0). Live:
`inspect --root fillcsr` (`FillingThresholdWitness`); full writeup in
[FILLING_THRESHOLD_CHAOS.md](FILLING_THRESHOLD_CHAOS.md).

## What This Does Not Answer

The GOE question that earlier topped this list (does the within-band
⟨r⟩ approach GOE/GUE as N grows?) is resolved in
Result 3: it does not, the dilute bands stay integrable. The dissipative-chaos
question is resolved in Result 5: a dense (extensive-filling) coherence sector of
the same Liouvillian does reach toward GinUE; chaos is a filling threshold. Two
genuinely open items remain.

1. **Comparison with Denisov lemon shape.** The complex-plane density
   of random Lindbladians (Denisov et al., PRL 2019) has a specific
   "lemon" shape (the characteristic boundary curve of eigenvalue density for structureless random Lindbladians). Our palindromic constraint modifies this. A
   quantitative comparison needs the 2D density, not just 1D rates.

2. **Topological dependence.** All results above use chain topology. Star,
   ring, and complete topologies were since surveyed with the complex
   spacing ratio in `simulations/rmt_topology_csr.py` (chain reads clean
   2D-Poisson; the symmetric topologies fragment into too few distinct
   levels for non-Hermitian RMT). A full sector-resolved comparison across
   topologies is still open.

---

## Connection to the Framework

The Poisson result is consistent with the picture the
[analytical formulas](../docs/ANALYTICAL_FORMULAS.md) draw, and it is
worth stating the implication in the direction it actually runs.
Closed forms would imply integrability, and integrability implies
Poisson; the measurement gives Poisson, which is the weakest of the
three and does not run back up the chain. In fact the repository has a
closed form for ONE block, the (0,1) coherence block (F2); the interior
blocks have none.

What does the work is the conserved structure, and the palindrome is
not it. The joint popcount cuts the Liouvillian into blocks that never
interact, and levels in different blocks have no reason to repel. The
palindromic constraint (F1) pairs the two halves of the spectrum
instead; it is a strong symmetry but it is not the one that suppresses
repulsion here.

This explains why the system does not fit any of the 38
Sa-Ribeiro-Prosen classes: those classes assume random matrix
statistics within each symmetry sector. Our system has the symmetry
of class AIII (chiral) but Poisson statistics instead of GUE. It is
an **integrable chiral Lindbladian**, which is not one of the 38
standard cases.

**Update (April 2026):** The spectral form factor (SFF) confirms
integrability (no dip-ramp-plateau) but reveals richer structure than
Poisson alone: palindromic modulation at ω_min = 4J(1-cos(π/N)),
w ↔ N-w band pairing in the time domain, and visibility scaling
as ~1/4^N. The spacing ratio (Poisson) describes local correlations;
the SFF describes global spectral structure. Both are consistent:
integrable + chiral = unique palindromic signature, neither standard
Poisson nor GUE. See [Spectral Form Factor](SPECTRAL_FORM_FACTOR.md).

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
  (Chiral RMT, class AIII)

---

*The 21,840 eigenvalues answer the question nobody asked: the
palindromic Liouvillian is integrable. What prevents chaos is the
joint-popcount grading, which cuts the spectrum into blocks that never
interact; the palindrome pairs those levels rather than separating
them. The spectrum is not random, and one of its blocks is in closed
form.*
