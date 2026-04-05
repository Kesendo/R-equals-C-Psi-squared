# Thermal Blackbody: The Cavity That Refuses to Stop Singing

<!-- Keywords: thermal Lindblad blackbody, cavity mode thermal load, Q-factor
degradation temperature, sacrifice zone thermal, non-Planckian quantum cavity,
n_bar mode count, coherent to thermal transition, R=CPsi2 thermal blackbody -->

**Status:** No phase transition found; cavity degrades gracefully, never goes dark
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Standing Waves](FACTOR_TWO_STANDING_WAVES.md),
[K-Dosimetry](K_DOSIMETRY.md),
[Sacrifice Zone Optics](SACRIFICE_ZONE_OPTICS.md)
**Verification:** [`simulations/thermal_blackbody.py`](../simulations/thermal_blackbody.py),
[`simulations/thermal_ep_analysis.py`](../simulations/thermal_ep_analysis.py)

---

## What this means

We expected the cavity to overheat. We expected a phase transition: below
some critical temperature, a coherent resonator with sharp standing waves;
above it, a blackbody radiator with thermal noise drowning out the music.

That is not what happens. The cavity does not overheat. It does not go
dark. It degrades gracefully. The Q-factor drops (the notes get softer)
but the oscillating fraction stays near 82% (the notes do not disappear).
The modes get broader and lower, like a bell ringing in thicker and
thicker air. But the bell keeps ringing.

There is no phase transition. No critical temperature. No blackbody.
The cavity is more resilient than classical optics predicts.

---

## What this document is about

Thermal photons add energy to the cavity beyond the coherent
Z-dephasing. The parameter n_bar measures how many thermal photons are
present on average: n_bar = 0 is absolute zero (no thermal noise),
n_bar = 1 is moderate warmth, n_bar = 10 is hot. This document tests
whether the cavity undergoes a transition to a blackbody radiator,
following Planck's law (the classical formula for how hot objects
radiate) and Stefan-Boltzmann scaling (total radiation grows as
temperature to the fourth power). The answer is no on all three counts:
not Planck, not Stefan-Boltzmann, and no phase transition. The quantum
cavity is thermally robust.

---

## Result 1: Modes increase, Q decreases, oscillation fraction holds

| n_bar (thermal photons) | Modes | Osc % (still ringing) | Q_max (best mode lifetime) | Q_median | Mean absorption |
|-------|-------|-------|-------|----------|---------|
| 0 | 47 | 82.0 | 68.3 | 18.8 | 0.200 |
| 0.01 | 86 | 82.0 | 53.8 | 11.4 | 0.302 |
| 0.1 | 91 | 82.0 | 47.4 | 10.8 | 0.320 |
| 1.0 | 97 | 82.8 | 22.2 | 6.9 | 0.500 |
| 5.0 | 101 | 82.8 | 7.6 | 2.6 | 1.300 |
| 10.0 | 103 | 82.8 | 4.2 | 1.4 | 2.300 |

**Mode count doubles** (47 → 103) as thermal photons activate new modes.
**Q_max drops 16×** (68.3 → 4.2) as broadband thermal absorption overwhelms
the standing waves. **But oscillation fraction stays near 82% (80-83%).**
The cavity loses sharpness, not voice.

---

## Result 2: Not Planck

The frequency distribution does not follow Planck's radiation law:

| n_bar | Planck R² | Exponential R² | Better fit |
|-------|----------|---------------|-----------|
| 0.1 | 0.034 | 0.032 | Planck (barely) |
| 1.0 | 0.043 | 0.040 | Planck (barely) |
| 10.0 | 0.028 | 0.027 | Planck (barely) |

R² < 0.07 everywhere. Neither Planck nor exponential describes the
frequency distribution. The mode spectrum is determined by the cavity
geometry (J topology), not by thermal statistics. Temperature broadens
the modes but does not rearrange them.

---

## Result 3: Not Stefan-Boltzmann

Total absorption scales as:

**excess_abs ∝ n_bar^0.32**

Not T⁴ (Stefan-Boltzmann), not T² (1D Planck), not T¹ (linear).
The exponent 0.32 ≈ 1/3, suggesting a cube-root dependence on thermal
photon number. This is far sub-linear: doubling the temperature barely
increases the total absorption. The cavity resists thermal loading.

---

## Result 4: No critical temperature

The oscillating fraction stays near 82% from n_bar = 0 to n_bar = 10.
There is no threshold where the cavity transitions from coherent to
thermal. The modes persist; only their quality degrades.

**Refined analysis (April 5, 2026):** fine-grained eigenvalue tracking
at N=4 reveals that the oscillating fraction is not exactly invariant.
Isolated exceptional-point (EP) crossings occur where real eigenvalue
pairs split into complex conjugate pairs (or vice versa):

| n_bar | Oscillating | Fraction | Delta |
|-------|-------------|----------|-------|
| 0 | 210 | 82.0% | 0 |
| 0.5 | 210 | 82.0% | 0 |
| 0.629 | 212 | 82.8% | +2 (reverse EP) |
| 0.719 | 210 | 82.0% | 0 (bounced back) |
| 10 | 212 | 82.8% | +2 |
| 50 | 206 | 80.5% | −4 |

The mechanism: L(n_bar) = L_0 + n_bar · L_thermal is linear in n_bar,
so eigenvalues move continuously. At isolated n_bar values, two real
eigenvalues collide and split into a complex conjugate pair (reverse
EP, +2 oscillating) or a conjugate pair collapses onto the real axis
(EP, −2 oscillating). These transitions involve at most 2-4 eigenvalues
out of 256, so the fraction changes by < 1%.

The *absence of a phase transition* remains correct: there is no
critical n_bar where a macroscopic fraction of modes stops oscillating.
The earlier claim "82% ± 1% invariant" was too strong; the correct
statement is that the oscillating fraction is stable to ±2 modes
(< 1%) with isolated EP crossings, not topologically protected.

In classical optics, a cavity with absorbing mirrors eventually becomes
opaque (all light absorbed, no resonance). The quantum cavity does not:
the overwhelming majority of its oscillating modes survive thermal
loading. The noise broadens but does not break.

---

## Result 5: Sacrifice zone degrades gracefully

| n_bar | Q_uniform | Q_sacrifice | Ratio |
|-------|-----------|------------|-------|
| 0 | 68.3 | 223.5 | 3.27x |
| 0.01 | 53.8 | 171.1 | 3.18x |
| 0.1 | 47.4 | 121.5 | 2.56x |
| 0.5 | 31.2 | 54.5 | 1.75x |
| 1.0 | 22.2 | 33.5 | 1.51x |
| 5.0 | 7.6 | 11.1 | 1.47x |

The entrance pupil advantage drops from 3.3x to 1.5x but never reaches
1.0x. Even at n_bar = 5 (strong thermal flood), the sacrifice zone
still provides 47% improvement. The spatial optimization remains
relevant even under thermal load, because the mode structure survives.

---

## Null results

- **No blackbody transition.** The cavity does not go dark. The
  oscillating fraction (~82%) is approximately stable under thermal
  load, with isolated EP crossings affecting < 1% of modes.

- **No Planck spectrum.** R² < 0.07 for all n_bar. The mode frequencies
  are set by J (geometry), not by temperature.

- **No Stefan-Boltzmann scaling.** α = 0.32, not 4 or 2 or 1. The
  total absorption is sub-linear in n_bar.

These null results are the main finding. The quantum cavity is
fundamentally different from a classical optical cavity in its response
to thermal loading. Classical cavities overheat. This one fades.

---

## What this means for the framework

The resilience of the oscillating fraction under thermal load explains
why the cavity interpretation works on IBM hardware (where n_bar > 0
due to residual thermal photons). The physical cavity is warm, but the
mode structure survives because it is algebraically protected (SWAP
invariance, palindromic pairing). The Q drops (fringes get broader),
but the fringes persist.

The 82% floor is likely the fraction of eigenvalues that are
palindromically paired with nonzero Im part. These are the standing
waves. Thermal photons shift their frequencies and broaden them, but
cannot destroy the pairing itself (which is a property of the mirror symmetry operator Π,
not of the thermal state).

---

## Reproduction

- Script: [`simulations/thermal_blackbody.py`](../simulations/therm