# Spectral Form Factor of the Palindromic Liouvillian

<!-- Keywords: spectral form factor Liouvillian, palindromic modulation SFF,
raw frequency spectral statistics, dip ramp plateau integrable, Heisenberg
chain Lindblad SFF, average-light decay-rate bands, raw density scale,
palindromic-pair period, R=CPsi2 spectral form factor -->

**Status:** Sampled raw-frequency diagnostic (Tier 2). N=2-7 (21,840 eigenvalues); not a universality-class proof.
**Date:** April 1, 2026
**Script:** [spectral_form_factor.py](../simulations/spectral_form_factor.py)
**Data:** [spectral_form_factor.txt](../simulations/results/spectral_form_factor.txt)
**Input:** Eigenvalue CSVs from [RMT analysis](RANDOM_MATRIX_THEORY.md)
(`simulations/results/rmt_eigenvalues_N{2..7}.csv`)
**Depends on:**
- [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) (raw spacing statistics)
- [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) (palindromic pairing)
- [Analytical Spectrum](ANALYTICAL_SPECTRUM.md) ((0,1) coherence block dispersion: ω_k = 4J(1-cos(πk/N)))

---

## What this document is about

We compute a raw non-unfolded frequency SFF of the palindromic Liouvillian
across N=2 to N=7. The reached window shows modulation;
this is compatible with integrability or symmetry/block fragmentation, not a
unique diagnosis of either. FFT components near the predicted F2/D10
(0,1)-block k=1 reference are identified at N=2-4 and N=6; at N=5 and N=7 that component is not cleanly
identified among the reported candidates. The exact pairing is an algebraic result;
the sampled modulation and its visibility are separate numerical readings.

---

## Abstract

The spectral form factor (SFF) is the standard time-domain diagnostic
for spectral statistics. Using the 21,840 complex Liouvillian eigenvalues
from the RMT analysis (N=2 through N=7), we compute the frequency SFF
K_freq(t) = (1/M²)|Σ exp(i·Im(λ_k)·t)|², where M=4^N is the eigenvalue count, and find:

1. **Sampled modulation candidates.** An FFT candidate lies within 1% of
   ω_F2 = 4J(1-cos(π/N)), the (0,1)-block k=1 reference, at N=2-4 and N=6; the producer does not identify
   it at N=5 or N=7. Each palindromic pair (μ, -μ)
   contributes 2cos(Im(μ)·t) to the trace amplitude; its square has doubled
   and cross frequencies. Such modulation alone is not a unique
   universality-class fingerprint.

2. **Modulation fades exponentially with N.** Visibility drops from
   0.72 (N=2) to 0.002 (N=7), roughly as ~1/4^N over this finite range.

3. **Period versus raw multiset density scale.** The modulation period
   t_Π = 2π/ω_F2 and the multiplicity-dependent raw scale are distinct
   readings. Their ratio is descriptive, not a separation of physical regimes.
   The beyond-scale bins at N=5-7 are not sampled.

4. **Reflected-band SFF matches palindromic pairing.** The sampled
   decay-rate bins centred at 2wγ and 2(N−w)γ have identical reported
   SFF statistics. By the Absorption Theorem they select average-light
   bands, not invariant fixed-XY-weight eigenspaces. At strictly positive
   uniform dephasing, the exact endpoint
   contains N+1 zero-frequency palindrome partners, so its normalized
   frequency SFF is the constant 1 at every time. At γ=0 the endpoints
   collapse into the larger commutator kernel, so the N+1 count is not used.

5. **Bounded diagnostic.** Modulation in the reached window is
   not a universality-class proof. The raw frequencies are
   not unfolded or reduced to irreducible symmetry blocks.

---

## Method

### Frequency SFF

For Liouvillian eigenvalues {λ_k}, the frequency SFF uses only
imaginary parts (oscillation frequencies):

    K_freq(t) = (1/M²) |Σ_k exp(i · Im(λ_k) · t)|², M = 4^N

This measures raw oscillation-frequency correlations, ignoring decay rates.
With this normalization K_freq(0)=1, and the independent-phase reference
is 1/M, not 1. This reference is not an asserted plateau: exact frequency
degeneracies retain additional cross terms in the long-time average.

Define the **raw multiset density scale** S = 2π / mean(diff(sort(f))),
where f contains every |Im λ| > 1e-10, retaining duplicates. It is
multiplicity-dependent: duplicating frequencies [1,2] to [1,1,2,2] changes
S from 2π to 6π, while the M²-normalized SFF is unchanged at every t.
Thus S is not a physical spectral-resolution invariant; no unfolding is performed.

The grid ends at min(3S, 50t_Π, 200). Its descriptive bins are below
(t < 0.1S), intermediate (0.1S < t < S), and beyond (t > S).
They are not physical time regimes or ramp/plateau evidence. Empty bins are
reported as `not sampled`, not zero. An intermediate-bin slope requires more
than 5 samples; an absent slope makes the heuristic not classifiable.
The helper tests distinguish a sampled zero from an empty bin and execute
the producer against a temporary output without changing the tracked report.

The producer computes only this raw oscillation-frequency SFF. It does not
define a decay-weighted SFF: physical Liouvillian propagation would involve
`exp(lambda*t)`, whereas `exp(i*lambda*t)` would weight `Im(lambda)`
exponentially and is not a decay observable. It also does not report a
connected SFF, which requires a specified ensemble or averaging prescription
before subtracting the disconnected contribution.

### Eigenvalue source

All eigenvalues from the C# engine (MKL/OpenBLAS), exported as CSV
by the RMT analysis: J=1.0, γ=0.05, Heisenberg chain, open boundaries.

---

## Result 1: Sampled Modulation Candidates

Each palindromic pair (μ, -μ) in the centered spectrum contributes:

    exp(i·Im(μ)·t) + exp(-i·Im(μ)·t) = 2·cos(Im(μ)·t)

The paired trace amplitude therefore contains cosine components. Squaring
that sum also introduces frequency sums and differences. The candidate
ω_F2 = 4J(1-cos(π/N)) is the k=1 mode of the (0,1) coherence block (F2),
not the smallest nonzero frequency of the full spectrum;
it need not be the largest FFT component of the SFF.

### FFT verification

| N | ω_F2 (block reference) | FFT peak | Match | Visibility |
|---|-------------------|----------|-------|------------|
| 2 | 4.000 | 3.998 | 0.1% | 0.718 |
| 3 | 2.000 | 1.999 | 0.1% | 0.232 |
| 4 | 1.172 | 1.162 | 0.8% | 0.042 |
| 5 | 0.764 | 2.010 | 163% (no match) | 0.013 |
| 6 | 0.536 | 0.534 | 0.4% | 0.007 |
| 7 | 0.396 | 0.848 | 114% (no match among candidates) | 0.002 |

The matched FFT amplitude rankings are N=2 first, N=3 first, N=4 third,
and N=6 second; each is within 1% of ω_F2. The producer compares the five
largest non-DC FFT amplitudes and prints the top three separately.
N=5 has no match among those five candidates (163% nearest error).
N=7 has no match either (114% nearest error); the output does not establish
that its ω_F2 component is present. The visibility column measures the
largest FFT amplitude relative to the sum, not the visibility of ω_F2.

### Visibility scaling

The modulation visibility decreases roughly as the inverse of the
number of eigenvalues (~1/4^N). At N=2 (16 eigenvalues), 6 oscillating
modes create a strong modulation. At N=7 (16,384 eigenvalues), 16,040
oscillating modes dephase and average out the modulation.

This is not a loss of structure: the palindromic pairing is exact at
every N. It is a loss of VISIBILITY: as more modes contribute, the
coherent modulation from any single pair is drowned in the sum.

---

## Result 2: Period and Descriptive Raw Density Scale

| N | t_Π | raw scale S | t_Π/S | ω_F2 | raw mean gap |
|---|-----|------|---------|-------|-------------|
| 2 | 1.57 | 25129 | 0.0001 | 4.000 | 0.00025 |
| 3 | 3.14 | 61.3 | 0.051 | 2.000 | 0.103 |
| 4 | 5.36 | 138.8 | 0.039 | 1.172 | 0.045 |
| 5 | 8.22 | 497.3 | 0.017 | 0.764 | 0.013 |
| 6 | 11.72 | 1649 | 0.007 | 0.536 | 0.0038 |
| 7 | 15.86 | 5810 | 0.003 | 0.396 | 0.0011 |

**F1-paired F2 reference period** t_Π = 2π/ω_F2 grows polynomially (~N² for large N,
since ω_F2 ~ 2Jπ²/N², or 2π²/N² at J=1, from the block dispersion relation).

The raw scale increases with the multiset density; it can increase merely
by duplicating modes without changing the SFF. Its ratio to t_Π is therefore
not a physical short/long-time separation.

The grid reaches at most t=200, below S=497.3, 1649, and 5810 at
N=5, 6, and 7. Their beyond bins are not sampled; N=7 does not reach
the intermediate bin. These facts describe coverage of arbitrary raw-scale
bins, not physical asymptotic behavior. Exact pairing has an independent
algebraic owner and is not restricted by this grid.

---

## Result 3: Decay-Rate-Band SFF

The producer groups full-spectrum eigenvalues into finite-width rate bins
`|d-2wγ|<γ`. These are decay-rate bands; under uniform γ their labels are
centres in average-light coordinates. They are not invariant integer-
XY-weight eigenspaces because the Hamiltonian mixes Pauli weights by ±2.
The sampled reflected bands have matching reported SFF statistics:

### N=5 band analysis

| Band centre w | Eigenvalues | <K_freq> | std(K) |
|----------|-------------|----------|--------|
| w=1 | 28 | 0.227 | 0.186 |
| w=2 | 478 | 0.023 | 0.064 |
| w=3 | 478 | 0.023 | 0.064 |
| w=4 | 28 | 0.227 | 0.186 |
| w=5 (endpoint) | 6 | 1.000 | 0.000 |

**Palindromic pairing:** the w=1 and w=4 bins have identical reported SFF
statistics, as do w=2 and w=3. D09 supplies the exact reason for consistently
reflected intervals: Π maps `d → 2Nγ-d` and negates every frequency while
preserving algebraic multiplicity. The labels w and N-w name average-light
band centres here, not eigenvalue sectors.

**Exact endpoint:** at the sampled γ=0.05>0, a connected N=5 chain has N+1=6 stationary modes and six
palindrome partners at `λ=-2Nγ`. Every endpoint frequency is zero, hence the
trace amplitude is constantly 6, the unnormalized frequency SFF constantly
36, and the normalized frequency SFF constantly 1. This is not an impulse at
`t=0`, nor a claim about the entire 2^N-dimensional XY-weight-N Pauli space.
The finite-width w=5 bin contains exactly those six modes in this run.

**Interior bands (w=2,3):** the sampled mean K is 0.023, reflecting a broader
frequency multiset in these bins. No universality classification follows from
that raw, non-unfolded number.

---

## Result 4: Modulation in the Reached Window

The raw non-unfolded frequency SFF shows modulation
in the reached window. This is compatible with integrability or symmetry/block
fragmentation, not a universality-class proof. An irreducible-block, unfolded
and suitably averaged statistic would be needed for a controlled RMT comparison.

The producer's below-bin mean and intermediate-bin slope thresholds describe only their
implemented numerical conditions. They do not certify Poisson, GUE, or an
physical asymptotic behavior. In particular the N=5-7 beyond-bin means
are `not sampled`, and the N=7 slope heuristic is not classifiable.

---

## Connection to Previous Results

**RMT (raw global spacing ratio):** The direct raw-multiset consecutive-gap
ratio retains zero gaps and counts undefined 0/0 separately (N=7: mean
0.2021456120489688 over 16,366 defined ratios, 16 undefined). It has no
standard-ensemble calibration for this degenerate unresolved population. It and
modulation in the reached window are separate finite-size observations. Their
combination does not prove integrability or establish a symmetry class. A
non-Hermitian class assignment requires reduction to irreducible strong-
symmetry sectors and the full sectorwise symmetry algebra.

**Palindrome identity:** each exact pair contributes
`2cos(ωt)` to the frequency trace amplitude. The SFF squares the complete
amplitude, producing doubled and cross frequencies. The raw-multiset spacing
statistic does not calibrate a universality class; spectral reflection alone neither assigns a
global AIII class nor implies Poisson statistics.

**Analytical spectrum (F2):** a reported FFT candidate matches
4J(1-cos(π/N)) to within 1% for N=2-4 and N=6. The producer does not identify
that candidate at N=5 or N=7; D10, not this sampled FFT association, proves
the (0,1) coherence-block dispersion.

**Topological analysis (geometric):** The SFF is independent of spatial
localization (it measures spectral correlations, not mode profiles).
The geometric localization and the spectral form factor are orthogonal
diagnostics. A symmetry-class label would require the irreducible-sector
reduction and its full symmetry algebra; the standing-wave structure and the
measured modulation require their own evidence.

---

## What This Does Not Answer

- A decay-sensitive or ensemble-connected non-Hermitian SFF; no such kernel or estimator is defined in this run
- SFF for non-chain topologies (ring, star, complete)
- SFF under non-uniform dephasing (sacrifice zone)
- Finite-size scaling of modulation visibility (limited to N=2-7)
- Whether the modulation peak at 2ω_F2 (second harmonic) carries
  independent information

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [`simulations/spectral_form_factor.py`](../simulations/spectral_form_factor.py) |
| Output | [`simulations/results/spectral_form_factor.txt`](../simulations/results/spectral_form_factor.txt) |
| Input CSVs | `simulations/results/rmt_eigenvalues_N{2..7}.csv` |

---

## References

- Sa, L., Ribeiro, P., Prosen, T. (2020). "Spectral and Steady-State
  Properties of Random Liouvillians." PRX 10, 021019.
- Cotler, J. et al. (2017). "Black Holes and Random Matrices."
  JHEP 2017, 118.
- Chan, A., De Luca, A., Chalker, J. (2018). "Spectral Statistics in
  Spatially Extended Chaotic Quantum Many-Body Systems." PRX 8, 041019.
- Gharibyan, H. et al. (2018). "Onset of Random Matrix Behavior in
  Scrambling Systems." JHEP 2018, 124.

---

*Exact pairing supplies cosine pairs; the raw SFF reads their combined
modulation over a finite window. The measured decay-rate bands and
visibility do not determine an unsampled large-N or late-time limit.*
