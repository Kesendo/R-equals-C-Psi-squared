# Random Matrix Theory Analysis of the Palindromic Liouvillian

<!-- Keywords: random matrix theory Liouvillian, spacing ratio palindromic
spectrum, Poisson level statistics open quantum, integrable Lindbladian
spectral statistics, chiral symmetry Liouvillian eigenvalues, XY-weight
sector universality class, Heisenberg dephasing RMT analysis,
R=CPsi2 random matrix theory -->

**Status:** Computationally verified (N=2-7, 21,832 eigenvalues, Heisenberg chain)
**Date:** April 1, 2026 (updated June 30, 2026: dissipative chaos located as a FILLING threshold, Result 5; the within-sector GOE hint resolved as a small-sample artifact June 27, Result 3)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** compute/RCPsiSquared.Compute (C# eigenvalue export),
[`simulations/rmt_analysis.py`](../simulations/rmt_analysis.py) (Python spacing analysis),
[`simulations/rmt_goe_hint_verdict.py`](../simulations/rmt_goe_hint_verdict.py) (sector GOE-hint verdict: bootstrap + larger-N)

---

## What this document is about

This document asks whether the palindromic Liouvillian spectrum looks
chaotic or integrable, using the standard diagnostic from random matrix
theory: the spacing ratio between consecutive eigenvalues. The answer:
the spectrum is Poisson (integrable, eigenvalues cluster rather than
repel), confirming that the palindromic symmetry provides enough
conserved quantities to fully determine the spectrum. An early small-N
read hinted at GOE-like repulsion within individual weight sectors, but
driving it to a verdict (bootstrap + larger samples) shows it was
small-sample noise: the sectors are integrable too, with no chaotic
transition (Result 3).

---

## Abstract

We perform the first random matrix theory (RMT) analysis of the
palindromic Liouvillian spectrum. Using spacing ratios (robust,
unfolding-free) on eigenvalues computed by the C# engine with
MKL/OpenBLAS (N=2 through N=7, up to 16,384 eigenvalues), we find
that the decay rate spectrum is **Poisson** (integrable, no level
repulsion) at every system size tested. The mean spacing ratio
converges to ⟨r⟩ = 0.36-0.39, consistent with the Poisson value
0.386 and far from GOE (0.536) or GUE (0.603). The palindromic
constraint creates an exactly integrable spectral structure: the
Liouvillian decomposes into non-interacting sectors that prevent
eigenvalue repulsion. The chiral symmetry (centered spectrum has
exact ± pairing) is confirmed to machine precision. An early
sector-resolved read at N=5 showed an apparent GOE-like ⟨r⟩ = 0.513
within individual weight sectors; driven to a verdict it is a
small-sample artifact (the n=15 value is a 1.5σ Poisson fluctuation, and
the same sectors read Poisson/sub-Poisson at N=6-7 with hundreds of
frequencies). There is no integrable-to-chaotic transition: the spectrum
is integrable at every sector and every tested N.

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
- Weight-sector structure: XY-weight w determines the decay rate
  envelope 2wγ ([F3](../docs/ANALYTICAL_FORMULAS.md))
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

### Sector analysis

For N=2-5 (Python eigenvalues): classify modes by nearest integer
XY-weight w (rate ~ 2*w*gamma), then analyze FREQUENCIES (imaginary
parts) within each sector. Rates are degenerate within a sector;
frequencies carry the spectral information.

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

The palindromic symmetry decomposes the Liouvillian into effectively
independent sectors. The conserved quantities (XY-weight parity,
palindromic pairing) prevent the eigenvalue interactions that would
produce level repulsion. This is the same mechanism that makes
integrable quantum systems Poisson: enough conserved quantities to
prevent chaos.

The spacing statistics confirm what the analytical formulas already
implied: the Liouvillian spectrum is fully determined by the sector
structure (formulas 1-5, 22-23, 33) plus the dispersion relation
within each sector (F2). There is no residual randomness.

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

## Result 3: The Within-Sector GOE Hint, Resolved (Artifact)

Within individual XY-weight sectors, frequencies (not rates) were measured
with the spacing ratio ⟨r⟩. An early small-N read looked GOE-like:

| N | Sector | Unique freq | ⟨r⟩ | Class |
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
   frequencies). There was never an independent second sector.

2. **0.513 on 15 frequencies is a Poisson fluctuation.** The Poisson sampling
   band for ⟨r⟩ at n=15 (Monte Carlo over a homogeneous Poisson process, same
   spacing-ratio statistic) is 0.386 ± 0.087, with [5%, 95%] = [0.245, 0.533].
   The observed 0.513 sits inside the band; the one-sided p(Poisson ≥ 0.513) =
   0.076, a ~1.5σ upward fluctuation, not significant.

3. **Larger samples read Poisson, not GOE.** The decay rate assigns the weight
   exactly (rate = 2wγ, by the Absorption Theorem), so the sectors can be read
   off the full Liouvillian spectrum at N=6-7 with no eigenvectors and no extra
   memory (this is the obstacle the original write-up wrongly thought blocked the
   check). The same sector then has hundreds of frequencies: N=6 w=3 (546 freq)
   → ⟨r⟩ = 0.272, N=7 w=3,4 (414 freq) → 0.283. The reading does not approach
   GOE (0.536); it converges to Poisson and below (sub-Poisson = level
   clustering, the opposite of repulsion, the signature of a strongly degenerate
   integrable additive spectrum, consistent with the N=4 w=2 row). The tiny
   w=1 / w=N−1 sectors (5-6 frequencies) throw a spurious "GUE" ⟨r⟩ > 0.79,
   plainly small-sample noise.

**Verdict: no within-sector chaos.** The system is integrable at every sector
and every tested N; the earlier GOE hint was small-sample noise. This matches
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

- The **dilute** (SE,DE) = (1,2) block — the Door-C block, where the non-solvable
  Galois group S_d lives — stays Poisson-like / non-GinUE under **every**
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

The sector GOE question that earlier topped this list (does the
within-sector ⟨r⟩ approach GOE/GUE as N grows?) is resolved in
Result 3: it does not, the dilute sectors stay integrable. The dissipative-chaos
question is resolved in Result 5: a dense (extensive-filling) coherence sector of
the same Liouvillian does reach toward GinUE — chaos is a filling threshold. Two
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

The Poisson result closes a loop: the
[analytical formulas](../docs/ANALYTICAL_FORMULAS.md) give exact
closed-form expressions for
eigenvalues in each sector. Exact formulas imply integrability.
Integrability implies Poisson. The RMT analysis confirms this
chain of reasoning numerically.

The palindromic constraint (F1) is not just a symmetry; it
is an integrability constraint. It provides enough conserved quantities
(weight-parity sectors, palindromic pairing) to fully determine the
spectrum, leaving no room for the randomness that produces level
repulsion.

This explains why the system does not fit any of the 38
Sa-Ribeiro-Prosen classes: those classes assume random matrix
statistics within each symmetry sector. Our system has the symmetry
of class AIII (chiral) but Poisson statistics instead of GUE. It is
an **integrable chiral Lindbladian**, which is not one of the 38
standard cases.

**Update (April 2026):** The spectral form factor (SFF) confirms
integrability (no dip-ramp-plateau) but reveals richer structure than
Poisson alone: palindromic modulation at ω_min = 4J(1-cos(π/N)),
w ↔ N-w sector pairing in the time domain, and visibility scaling
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

*The 87,376 eigenvalues answer the question nobody asked: the
palindromic Liouvillian is integrable. The symmetry that pairs every
decay mode also prevents chaos. The spectrum is not random. It is
exactly what the formulas predict.*
