# Spectral Form Factor of the Palindromic Liouvillian

<!-- Keywords: spectral form factor Liouvillian, palindromic modulation SFF,
dissipative spectral statistics, dip ramp plateau integrable, Heisenberg
chain Lindblad SFF, XY-weight sector spectral, Heisenberg time palindromic
time, R=CPsi2 spectral form factor -->

**Status:** Computationally verified (Tier 2). N=2-7 (21,832 eigenvalues).
**Date:** April 1, 2026
**Script:** [spectral_form_factor.py](../simulations/spectral_form_factor.py)
**Data:** [spectral_form_factor.txt](../simulations/results/spectral_form_factor.txt)
**Input:** Eigenvalue CSVs from [RMT analysis](RANDOM_MATRIX_THEORY.md)
(simulations/results/rmt_eigenvalues_N{2..7}.csv)
**Depends on:**
- [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) (Poisson, class AIII)
- [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) (Pi chiral, palindromic pairing)
- [Analytical Spectrum](ANALYTICAL_SPECTRUM.md) (w=1 dispersion: ω_k = 4J(1-cos(πk/N)))

---

## Abstract

The spectral form factor (SFF) is the standard time-domain diagnostic
for spectral statistics. Using the 21,832 complex Liouvillian eigenvalues
from the RMT analysis (N=2 through N=7), we compute the frequency SFF
K_freq(t) = (1/N²)|Σ exp(i·Im(λ_k)·t)|² and find:

1. **Palindromic modulation confirmed.** The SFF contains periodic
   modulation at the predicted frequency ω_min = 4J(1-cos(π/N)),
   matching to <1% for N=2-4 and N=6. Each palindromic pair (μ, -μ)
   contributes cos(Im(μ)·t), creating a fingerprint unique to the
   palindromic symmetry.

2. **Modulation fades exponentially with N.** Visibility drops from
   0.72 (N=2) to 0.002 (N=7), roughly as ~1/4^N. At large N, the
   many modes dephase and the modulation averages out.

3. **Two timescales separate.** The palindromic time t_Π = 2π/ω_min
   (modulation period) grows as ~N², while the Heisenberg time
   t_H = 2π/Δ (spectral resolution) grows as ~4^N. Their ratio
   t_Π/t_H → 0: the palindromic structure lives in the short-time
   regime, invisible to long-time diagnostics.

4. **Sector SFF confirms palindromic pairing.** Sectors w and N-w have
   identical SFF statistics (same mean, same variance). The XOR sector
   (w=N) has K=1.0 exactly (all eigenvalues degenerate at rate 2Nγ).

5. **Neither Poisson nor GUE.** The SFF does not match the Poisson
   prediction (K=1/N, flat) nor the GUE prediction (dip-ramp-plateau).
   It is a distinct palindromic signature: modulated, decaying, with
   sector structure.

---

## Method

### Frequency SFF

For Liouvillian eigenvalues {λ_k}, the frequency SFF uses only
imaginary parts (oscillation frequencies):

    K_freq(t) = (1/N²) |Σ_k exp(i · Im(λ_k) · t)|²

This measures correlations between oscillation frequencies, ignoring
decay rates. It is well-defined for all t (no overflow from decay).

The dissipative SFF K_diss(t) = (1/N²) Σ_{j,k} exp((λ_j+conj(λ_k))t)
(Sa & Prosen 2020) was not computed because eigenvalues with large
negative imaginary parts cause numerical overflow at moderate t. A
proper implementation would require filtering to the slow-decay subspace.

### Eigenvalue source

All eigenvalues from the C# engine (MKL/OpenBLAS), exported as CSV
by the RMT analysis: J=1.0, γ=0.05, Heisenberg chain, open boundaries.

---

## Result 1: Palindromic Modulation

Each palindromic pair (μ, -μ) in the centered spectrum contributes:

    exp(i·Im(μ)·t) + exp(-i·Im(μ)·t) = 2·cos(Im(μ)·t)

The SFF therefore contains cosine components at every palindromic
frequency, with the dominant modulation at the lowest frequency
ω_min = 4J(1-cos(π/N)) (the k=1 mode of the w=1 sector, formula 2).

### FFT verification

| N | ω_min (predicted) | FFT peak | Match | Visibility |
|---|-------------------|----------|-------|------------|
| 2 | 4.000 | 3.998 | 0.1% | 0.718 |
| 3 | 2.000 | 1.999 | 0.1% | 0.232 |
| 4 | 1.172 | 1.162 | 0.8% | 0.042 |
| 5 | 0.764 | 2.010 | 163% (no match) | 0.013 |
| 6 | 0.536 | 0.534 | 0.4% | 0.007 |
| 7 | 0.396 | 3.391 | (no match) | 0.002 |

For N=2-4 and N=6: the dominant FFT peak matches ω_min to within 1%.
The palindromic modulation is confirmed in the time domain.

For N=5 and N=7: the ω_min peak is present but subdominant. Higher
harmonics and inter-sector beating dominate the FFT. The visibility
is too low (<1.5%) for clean identification.

### Visibility scaling

The modulation visibility decreases roughly as the inverse of the
number of eigenvalues (~1/4^N). At N=2 (16 eigenvalues), 6 oscillating
modes create a strong modulation. At N=7 (16,384 eigenvalues), 16,040
oscillating modes dephase and average out the modulation.

This is not a loss of structure: the palindromic pairing is exact at
every N. It is a loss of VISIBILITY: as more modes contribute, the
coherent modulation from any single pair is drowned in the sum.

---

## Result 2: Two Timescales Separate

| N | t_Π | t_H | t_Π/t_H | ω_min | Δ (spacing) |
|---|-----|------|---------|-------|-------------|
| 2 | 1.57 | 25129 | 0.0001 | 4.000 | 0.00025 |
| 3 | 3.14 | 61.3 | 0.051 | 2.000 | 0.103 |
| 4 | 5.36 | 138.8 | 0.039 | 1.172 | 0.045 |
| 5 | 8.22 | 497.3 | 0.017 | 0.764 | 0.013 |
| 6 | 11.72 | 1649 | 0.007 | 0.536 | 0.0038 |
| 7 | 15.86 | 5810 | 0.003 | 0.396 | 0.0011 |

**Palindromic time** t_Π = 2π/ω_min grows polynomially (~N² for large N,
since ω_min ~ π²/N² from the dispersion relation).

**Heisenberg time** t_H = 2π/Δ grows exponentially (~4^N, set by the
total number of eigenvalues).

Their ratio t_Π/t_H → 0: the palindromic modulation lives in an
exponentially shrinking fraction of the SFF time window. At N=7,
the palindromic period is 0.3% of the Heisenberg time.

Physical meaning: the palindromic structure (paired modes, mirrored
decay rates) is a SHORT-TIME phenomenon. The long-time spectral
correlations (level repulsion, spectral rigidity) are set by the
Poisson statistics of the integrable system. The palindrome organizes
the fast dynamics; Poisson governs the slow dynamics.

---

## Result 3: Sector-Resolved SFF

The SFF computed separately for each XY-weight sector confirms the
palindromic pairing at the sector level:

### N=5 sector analysis

| Sector w | Eigenvalues | <K_freq> | std(K) |
|----------|-------------|----------|--------|
| w=1 | 28 | 0.227 | 0.186 |
| w=2 | 478 | 0.023 | 0.064 |
| w=3 | 478 | 0.023 | 0.064 |
| w=4 | 28 | 0.227 | 0.186 |
| w=5 (XOR) | 6 | 1.000 | 0.000 |

**Palindromic pairing:** w=1 and w=4 have identical SFF statistics.
w=2 and w=3 have identical SFF statistics. This is the palindromic
symmetry in the time domain: Π maps w → N-w, so sectors w and N-w
must have the same spectral structure.

**XOR sector (w=N):** K=1.000 exactly, with zero variance. All XOR
eigenvalues are degenerate at rate 2Nγ (the maximum decay rate). They
oscillate at exactly the same frequency → perfect correlation → K=1.

**Interior sectors (w=2,3):** Much lower K (0.023), indicating more
spectral diversity (many distinct frequencies, less correlation).

---

## Result 4: Neither Poisson Nor GUE

The SFF does not match either standard universality class:

- **Poisson predicts** K(t>0) = 1/N (flat, no temporal structure).
  Observed: K has strong time dependence with palindromic modulation.

- **GUE predicts** dip at t=0+, linear ramp to K=1/N at t_H, plateau.
  Observed: no ramp (K decreases or stays flat, no linear rise).

The palindromic Liouvillian has a UNIQUE SFF signature: modulated
oscillation at short times (t < t_Π), decay to a low plateau at
intermediate times, no ramp. This is consistent with integrability
(Poisson spacing) combined with the palindromic modulation (chiral
symmetry).

---

## Connection to Previous Results

**RMT (Poisson, class AIII):** The SFF confirms integrability: no
dip-ramp-plateau, no level repulsion signature. The Poisson spacing
ratio (⟨r⟩ = 0.383) manifests as absence of ramp in the SFF.

**PT analysis (Π chiral, order 4):** The palindromic modulation in the
SFF is the TIME-DOMAIN fingerprint of the same λ ↔ -(λ+2Σγ) pairing
that the PT analysis proved algebraically. The SFF sees it as a cosine;
the level statistics sees it as Poisson; the Π operator explains both.

**Analytical spectrum (formula 2):** The modulation frequency ω_min
matches the predicted 4J(1-cos(π/N)) to within 1% for N=2-4, 6. This
confirms the w=1 dispersion relation in the time domain.

**Topological analysis (geometric):** The SFF is independent of spatial
localization (it measures spectral correlations, not mode profiles).
The geometric localization and the spectral form factor are orthogonal
diagnostics that agree on the same underlying structure: the palindromic
Liouvillian is integrable, chiral (AIII), and organized by standing
wave modes.

---

## What This Does Not Answer

- Dissipative SFF K_diss (overflow issue; requires spectral filtering)
- SFF for non-chain topologies (ring, star, complete)
- SFF under non-uniform dephasing (sacrifice zone)
- Finite-size scaling of modulation visibility (limited to N=2-7)
- Whether the modulation peak at 2ω_min (second harmonic) carries
  independent information

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | simulations/spectral_form_factor.py |
| Output | simulations/results/spectral_form_factor.txt |
| Input CSVs | simulations/results/rmt_eigenvalues_N{2..7}.csv |

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

*The SFF is the Fourier transform of "how are the eigenvalues
distributed?" The palindromic answer: in cosine pairs, organized by
XY-weight sectors, integrable at every N, with a modulation that fades
as 1/4^N but never disappears. The structure is exact; the visibility
is finite.*
