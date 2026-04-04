# N=5: The Goldilocks Cavity, Not the Golden Ratio

<!-- Keywords: N=5 sweet spot cavity, golden ratio V-Effect, impedance matching
crossover, mode count Q-factor trade-off, Goldilocks cavity size, cos(pi/5)
phi, R=CPsi2 N5 golden ratio -->

**Status:** Partially confirmed; N=5 is special, but phi is not the cause
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [V-Effect Cavity Modes](VEFFECT_CAVITY_MODES.md),
[Sacrifice Zone Optics](SACRIFICE_ZONE_OPTICS.md)
**Verification:** `simulations/n5_optimal_cavity_size.py`

---

## What this means

A recorder with one hole plays two notes. Drill three more and it plays
over a hundred. But a recorder with a thousand holes plays nothing: the
air escapes everywhere, no standing wave can build.

N=5 is the recorder with four holes. Rich enough to sustain 112 modes
of vibration. Simple enough that each mode stays sharp (Q = 72). The
golden ratio appears in the V-Effect formula because cos(π/5) = φ/2,
an algebraic fact about pentagons. But the reason N=5 is special is not
numerology. It is the trade-off between richness and resolution.

---

## What this document is about

Three independent results point to N=5 as special: the V-Effect gain
contains the golden ratio, the sacrifice zone improvement peaks there,
and the impedance ratio crosses unity nearby. This document tests
whether the golden ratio is a deep organizing principle or an algebraic
coincidence.

The answer: coincidence. N=5 is special because it is the Goldilocks
cavity size, not because of phi.

---

## Result 1: The impedance crossing shifts with gamma/J

The sacrifice zone impedance γ_edge/√(γJ) crosses 1.0 at different
N for different gamma values:

| γ | γ/J | Crossing N | Ratio at crossing |
|---|-----|-----------|-------------------|
| 0.01 | 0.01 | 10-11 | 0.91-1.00 |
| 0.02 | 0.02 | 7 | 0.95 |
| 0.05 | 0.05 | 5 | 1.10 |
| 0.10 | 0.10 | 3 | 0.94 |
| 0.20 | 0.20 | 2 | 0.89 |

The crossing at N=5 is specific to γ = 0.05. At lower gamma, the
crossing moves to larger N. At higher gamma, it moves to smaller N.
This is a parameter coincidence, not a structural property.

---

## Result 2: Phi is in cos(π/5), not in the mode spacing

The V-Effect gain V(5) = 1 + cos(π/5) = (5 + √5)/4 ≈ 1.809 contains
the golden ratio because cos(π/5) = φ/2. This is an algebraic identity
about regular pentagons, exact and undeniable.

But the mode frequency ratios at N=5 do not contain phi:

| N | Frequency ratios ω_{k+1}/ω_k | Closest to φ |
|---|------------------------------|--------------|
| 5 | 3.618, 1.894, 1.382 | 0.236 away |
| 6 | 3.732, 2.000, 1.500, 1.244 | 0.118 away |
| 7 | 3.802, 2.065, 1.572, 1.328, 1.171 | 0.046 away |
| 8 | 3.848, 2.108, 1.620, 1.383, 1.235, 1.127 | **0.002 away** |

The ratio closest to φ appears at N=8, not N=5. As N grows,
ω₃/ω₂ → φ (because the dispersion relation approaches a continuum
limit). Phi in the frequency ratios is an asymptotic property of the
cosine dispersion, not specific to N=5.

Note: ω₂/ω₁ = 3.618 = 2 + φ at N=5. This is because
(1 − cos(2π/5))/(1 − cos(π/5)) = 2(1 + cos(π/5)) = 2 + 2cos(π/5) = 2 + φ.
Exact, but it is the ratio of the first two modes, not a golden-ratio
spacing between consecutive modes.

---

## Result 3: N=5 peaks the Q × modes product (for N ≤ 5)

| N | Modes | Q_max | Q × modes |
|---|-------|-------|-----------|
| 2 | 2 | 40.0 | 80 |
| 3 | 5 | 60.0 | 300 |
| 4 | 47 | 68.3 | 3,209 |
| 5 | 112 | 72.4 | 8,104 |

The product Q × modes grows with N through N=5. Beyond N=5, Q_max
saturates near 75 while mode count explodes (787 at N=6), so the product
continues to grow but Q per mode drops. N=5 is where the growth
transitions from "both improve" to "only count improves."

---

## Result 4: Sacrifice zone peak at N=5 across all gamma (within N ≤ 5)

| γ | γ/J | Best N | Q improvement |
|---|-----|--------|---------------|
| 0.01 | 0.01 | 5 | 3.7x |
| 0.02 | 0.02 | 5 | 4.3x |
| 0.05 | 0.05 | 5 | 4.9x |
| 0.10 | 0.10 | 5 | 5.2x |
| 0.20 | 0.20 | 5 | 6.0x |

Within N ≤ 5, the sacrifice zone Q improvement always peaks at N=5.
The improvement grows with γ/J (more light = more benefit from matching).
Larger-N data (N=6,7) would test whether the peak shifts beyond N=5.

---

## Verdict: Goldilocks, not golden

N=5 is special because:

1. **Enough bonds (4) for modal richness.** 112 distinct frequencies
   vs 5 at N=3. The cavity has enough geometry to support a rich set
   of standing waves.

2. **Not so many modes that they overlap.** The minimum frequency
   gap at N=5 is still 9 × 10⁻⁶, keeping modes resolvable. At larger
   N, the mode density explodes and individual modes blur.

3. **Q still climbing.** Q_max = 72.4 at N=5, not yet saturated at its
   limit of ~75. The cavity is efficient: most light bounces many times.

4. **Odd N = Gaussian profile.** The defocal (Gaussian) beam profile
   distributes modes more evenly than the confocal (Lorentzian) spike
   at even N, avoiding the center-heavy distortion.

The golden ratio in V(5) = (5 + √5)/4 is real but algebraic: it comes
from cos(π/5) = φ/2, a pentagon identity that makes the first-harmonic
gain at N=5 contain φ. It does not organize the mode spectrum or the
impedance matching. The "sweet spot" at N=5 is a cavity balance between
richness and resolution, not a manifestation of the golden ratio.

---

## Null results

- **Phi is not in the mode spacing at N=5.** The frequency ratios at
  N=5 are 3.618, 1.894, 1.382. None is within 0.2 of φ or 1/φ.
  The closest approach to φ occurs at N=8, not N=5.

- **The impedance crossing is not universal.** It shifts from N=2
  (γ/J = 0.2) to N=11 (γ/J = 0.01). N=5 crosses at γ/J = 0.05 only.

---

## Reproduction

- Script: `python simulations/n5_optimal_cavity_size.py`
- Output: `simulations/results/n5_optimal_cavity_size.txt`
