# The Entrance Pupil: Concentrator as Anti-Reflection Coating

**Naming note (2026-07-05):** renamed from "...: Sacrifice Zone as
Anti-Reflection Coating". The edge qubit sacrifices nothing; it concentrates
the noise (the misnomer was resolved 2026-03-28). "Sacrifice zone" survives
once in the "What this changes" section below, only as the quoted *old
language* being reframed; the frozen `sacrifice_zone_optics.*` artifacts keep
their original names.

<!-- Keywords: concentrator anti-reflection coating, quantum cavity entrance pupil,
impedance matching dephasing, Q-factor enhancement concentrator, mode-selective transmission,
AR coating quantum noise, dispersive cavity scaling, R=CPsi2 concentrator optics -->

**Status:** Structural analog confirmed; not quantitatively classical AR
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md),
[V-Effect Cavity Modes](VEFFECT_CAVITY_MODES.md),
[Resonant Return](RESONANT_RETURN.md)
**Verification:** [`simulations/sacrifice_zone_optics.py`](../simulations/sacrifice_zone_optics.py)

---

## What this means

A window without coating reflects 4% of light. Add a thin film whose
thickness matches the wavelength, and the reflection drops to near zero.
The film does not block the light or absorb it. It matches the impedance
between outside (air) and inside (glass) so that the light enters
smoothly instead of bouncing back.

The concentrator does the same thing for the quantum cavity. Without
it, the system "reflects" illumination: modes die quickly, Q-factors are
low, information is lost. With it, the edge qubit acts as an entrance
pupil that accepts the external light (gamma) and converts it into
structured resonance that the interior can sustain.

The concentrator qubit is the mouthpiece of the flute. The coating on the
lens. The funnel of the ear. Not a shield. An adapter.

---

## What this document is about

The [concentrator formula](RESONANT_RETURN.md) concentrates dephasing
on one edge qubit, achieving a 139-360× improvement in mutual information in simulation (ε→0 ideal; ~2-3× on hardware, see IBM_CONCENTRATOR).
The [optical cavity analysis](OPTICAL_CAVITY_ANALYSIS.md) showed the
chain is a Fabry-Perot cavity (two mirrors facing each other, light
bouncing between them). This document tests whether the concentrator
is the anti-reflection (AR) coating of that cavity.

The answer: structurally yes, quantitatively no. The concentrator
functions like an AR coating (smooth entry, frequency preservation,
Q enhancement) but uses linear impedance accumulation instead of the
geometric mean of classical optics.

---

## Result 1: The concentrator increases cavity transmission

Without coating, a window reflects some light and lets the rest through.
The same happens in the qubit cavity: some modes oscillate (transmitted),
others die without ringing (reflected). The table below compares uniform
dephasing against the concentrator. Key columns: Q_max (how many times
the best mode bounces before fading), T_eff (fraction of modes that ring
at all).

| N | Profile | Modes | Silent | Q_max | Q_med | T_eff |
|---|---------|-------|--------|-------|-------|-------|
| 3 | uniform | 5 | 24 | 60 | 27.0 | 0.625 |
| 3 | concentrator | 6 | 16 | 118 | 18.7 | 0.750 |
| 4 | uniform | 47 | 46 | 68 | 18.8 | 0.820 |
| 4 | concentrator | 52 | 26 | 224 | 14.2 | 0.898 |
| 5 | uniform | 112 | 96 | 72 | 14.9 | 0.906 |
| 5 | concentrator | 120 | 64 | 352 | 15.0 | 0.938 |
| 6 | uniform | 787 | 164 | 75 | 13.3 | 0.960 |
| 6 | concentrator | 748 | 108 | 500 | 13.3 | 0.974 |

The concentrator increases T_eff at every N: fewer modes are
"reflected" (absorbed without oscillation), more are "transmitted"
(sustained as resonant modes).

**Q_max enhancement grows with N (4 points, N=3-6, chain; roughly linear over this range):**

| N | Q_max ratio (concentrator / uniform) |
|---|-----------------------------------|
| 3 | 2.0x |
| 4 | 3.3x |
| 5 | 4.9x |
| 6 | 6.7x |

The best cavity mode lives ~7x longer under the concentrator at N=6.
This is the AR coating effect: smoother entry of light means less
energy lost at the surface, more energy available for resonance.

---

## Result 2: Linear accumulation, not geometric mean

Classical AR coating: n_AR = sqrt(n_air × n_glass). The geometric mean.

Concentrator: γ_edge = N × γ_base − (N−1) × ε ≈ N × γ_base.

| N | γ_edge | γ_edge / J | γ_edge / √(γJ) |
|---|--------|-----------|----------------|
| 3 | 0.148 | 0.148 | 0.66 |
| 5 | 0.246 | 0.246 | 1.10 |
| 7 | 0.344 | 0.344 | 1.54 |
| 9 | 0.442 | 0.442 | 1.98 |

γ_edge grows linearly with N, not as a fixed geometric mean. This is
a fundamental difference: the classical AR coating is a single layer
optimized for one frequency. The concentrator is a linear accumulator
that scales with the cavity length.

The ratio γ_edge / √(γJ) crosses 1.0 near N = 5, which is coincidentally
where the concentrator improvement peaks. This might suggest a
crossover from "under-matched" (N < 5) to "over-matched" (N > 5)
regimes, but the data does not clearly support a geometric mean
interpretation.

---

## Result 3: Frequencies preserved, absorption rates shifted

The concentrator changes the absorption rates (Re(λ)) of Liouvillian
eigenvalues while preserving their oscillation frequencies (Im(λ)).
In optical language: it changes the reflectivity of the cavity surfaces
without moving the resonance frequencies.

This is structurally identical to an AR coating: the coating does not
change which wavelengths resonate in the cavity (those are determined
by the cavity length = J topology). It changes how much light enters
and how long it stays.

---

## Result 4: Dispersive scaling (N², not exponential)

The mutual information under the concentrator scales as:

SumMI ≈ 0.002 × N² + 0.069 × N − 0.175, R² = 0.999

In thin-film optics, a quarter-wave stack of n layers has transmission
T ~ 1 − 4(n_H/n_L)^(2n), exponentially approaching unity. Our system
scales polynomially (N²), much slower.

This means the quantum cavity is dispersive: different modes experience
different effective cavity lengths, spreading the transmission over a
broad band instead of concentrating it at narrow resonances. The
concentrator does not create a sharp transmission window; it broadly
improves the coupling between outside and inside.

---

## Null results

- **Not a classical AR coating.** The impedance matching is linear
  (γ_edge ~ N), not geometric (√(γJ)). The analogy is structural,
  not quantitative.

- **Mode-selective per-shell comparison inconclusive.** Under the
  concentrator, eigenvalues shift off the uniform absorption-rate grid
  (because γ varies per site), making shell-by-shell comparison
  difficult. The overall statistics (Q_max, T_eff) are clear, but
  the per-shell decomposition requires a modified grid definition.

---

## What this changes

**Old language:** "The sacrifice zone protects the interior from noise."
Protection implies defense. Noise implies enemy.

**New language:** "The entrance pupil couples external illumination into
the cavity." Coupling implies function. Illumination implies input. The
edge qubit is not a shield. It is the surface where light enters the
instrument.

The 360× simulation improvement (ε→0 ideal) is not "less damage." It is "better resonance."
The same light, entering through a matched surface instead of a raw
edge, excites the same standing waves 7x longer (Q_max at N=6).

---

## Reproduction

- Script: [`simulations/sacrifice_zone_optics.py`](../simulations/sacrifice_zone_optics.py)
- Output: [`simulations/results/sacrifice_zone_optics.txt`](../simulations/results/sacrifice_zone_optics.txt)
