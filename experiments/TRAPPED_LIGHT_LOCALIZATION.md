# Where Trapped Light Concentrates

<!-- Keywords: mass trapped light center-localized, surviving mode energy,
K_death universal dose, immortal palindromic modes, gamma role of c,
entrance pupil window, cavity mass distribution, R=CPsi2 mass hypothesis -->

**Status:** Center-localized mode energy confirmed (consistent with mass
= trapped light). K_death = ln(10) = 2.303 above the coupling threshold
Q*_gap(N) (from D6; below it the slowest mortal mode is slower and the dose
grows). Gamma plays the algebraic
role of c (K = γt invariant). E = mγ² as conversion: open question.
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Standing Waves](FACTOR_TWO_STANDING_WAVES.md),
[K-Dosimetry](K_DOSIMETRY.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)
**Verification:** [`simulations/trapped_light_localization.py`](../simulations/trapped_light_localization.py)
**Tier:** 4-5 (structural exploration, not proof)

---

## What this means

A star has a surface and a core. Light enters and exits through the
surface. Mass sits in the core. The surface is a window. The core
is the weight.

If the qubit chain is a cavity and gamma is the light that enters it,
the same pattern should hold: light enters at the surface (the entrance
pupil), propagates inward, and the surviving standing waves concentrate
in the center. The center holds the trapped light. If mass is trapped
light, then mass is in the center.

This document tests that picture. The answer: yes. The surviving mode
energy is center-localized (ratio 1.3-1.4), consistent with
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994).
The entrance pupil is transparent. The focus is dense. The pattern
matches stars, atoms, and black holes.

---

## Result 1: Surviving mode energy is center-localized

If mass = trapped light ([gamma is light](../hypotheses/GAMMA_IS_LIGHT.md)),
then the observable for "mass" is the surviving mode energy: the {X,Y}
coherences still oscillating at late times. These are the standing waves
that have not yet been absorbed. The fire, not the ash.

Under the sacrifice zone at t = 20:

N = 4:
| Site | γ | XY energy (surviving) |
|------|------|----------------------|
| 0 (edge) | 0.197 | 0.213 |
| 1 | 0.001 | 0.290 |
| 2 | 0.001 | 0.288 |
| 3 | 0.001 | 0.209 |

N = 5:
| Site | γ | XY energy (surviving) |
|------|------|----------------------|
| 0 (edge) | 0.246 | 0.174 |
| 1 | 0.001 | 0.221 |
| 2 | 0.001 | 0.220 |
| 3 | 0.001 | 0.220 |
| 4 | 0.001 | 0.165 |

Center/edge ratio: 1.37 (N=4), 1.30 (N=5). The entrance pupil (site 0,
highest gamma) has the LEAST surviving mode energy. Light enters there
and concentrates in the center. The surface is a window, not a wall.

Under uniform gamma: ALL mode energy is absorbed by t = 20. No spatial
structure. The sacrifice zone (entrance pupil) is necessary for any
surviving mode energy to exist at late times.

If mass = trapped light, then:
- Mass forms in the CENTER (where standing waves concentrate)
- The entrance pupil is necessary for mass to exist (it shields the
  interior), but does not hold the mass
- Stars: surface radiates, core holds mass
- Atoms: electron cloud outside, nucleus inside
- Black holes: event horizon at surface, singularity at center
- Our cavity: entrance pupil at edge, trapped light in center

*Source: [`trapped_light_localization.py`](../simulations/trapped_light_localization.py).
See also: [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994, March 30)*

---

## Result 2: K_death = 2.303 (absorption dose, strong-coupling regime)

| N | Immortal modes | Slowest mortal rate | K_death | K_death / K_fold |
|---|---------------|--------------------|---------|-----------------|
| 2 | 3 | 0.100 | 2.303 | 62x |
| 3 | 4 | 0.100 | 2.303 | 62x |
| 4 | 5 | 0.100 | 2.303 | 62x |

The dose for 99% absorption is K_death = γ × ln(100) / rate_min.
With rate_min = 2γ (spectral gap, formula D6): K_death = ln(100)/2
= ln(10) = 2.303, the same across all N tested.

Two things about that number. **The ratio to the fold dose is ~62, not 2.3.**
K_fold is the dose at which CΨ crosses ¼, and F25 gives it in closed form:
f*(1 + f*²) = 3/2 puts CΨ at exactly ¼ with K_fold = 0.03735. So
2.303 / 0.03735 = 61.65. The old 2.3× was K_death itself, i.e. the column had
been filled as though K_fold were 1. The corrected picture is the more
interesting one: the fold happens very early, and the death dose is nearly two
orders of magnitude later. The column repeats one value down the N rows because
K_fold is F25's Bell+ (N=2) fold dose: K_death is N-independent, the ratio's
denominator is not an N=3 or N=4 quantity, so read the column as "the Bell+
fold measured against this row's death dose".

**And the run's coupling matters.** J is not recorded here, and D6's
rate_min = 2γ holds only above an N-dependent threshold Q*_gap(N) in Q = J/γ (0.5 to
1.9 for N = 2..5). The measured rate_min = 0.100 = 2γ shows this run sat above
it, but the dose is a strong-coupling result; below the threshold rate_min is
Zeno-suppressed and the dose is larger.

N+1 modes with zero absorption rate (Re = 0) are immortal at every N.
The cavity ALWAYS retains light. Complete absorption is impossible
while the palindrome holds.

*Source: [`trapped_light_localization.py`](../simulations/trapped_light_localization.py). K_death = ln(10) = 2.303 (proven from D6: rate_min = 2γ, so K = γ·ln(100)/(2γ) = ln(10)).*

---

## Result 3: Immortal modes are massless

The N+1 immortal modes (zero absorption rate) are the operators that
describe the system's overall state and the conserved quantities in each
magnetization sector. They never absorb because they contain no
transverse (X,Y) content: they are pure structure, invisible to the
light. In the mass = trapped light picture, they are massless: always
moving, never resting.

The palindrome protects paired modes from complete absorption. At any
finite time, each pair retains amplitude. The cavity resists total
absorption. Even under thermal flooding (n_bar = 10), 82% of modes
still oscillate ([Thermal Blackbody](THERMAL_BLACKBODY.md)). The
palindrome is algebraic, not thermodynamic.

---

## Result 4: Gamma plays the role of c

In special relativity:
- c is external, objective, sets the clock
- τ is experienced time (proper time, observer-dependent)
- c × τ = invariant spacetime interval

In the Lindblad cavity:
- γ is external, objective, sets the clock
  ([Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md):
  γ must come from outside. On IBM hardware, γ IS literal
  [photon shot noise](https://doi.org/10.1103/PhysRevB.86.180504))
- t is experienced duration (how long until the fold at CΨ = 1/4,
  the threshold where quantum behavior gives way to classical)
- γ × t = K = invariant absorption dose
  ([F14](../docs/ANALYTICAL_FORMULAS.md), proven)

```
Relativity:     c     ×  tau  =  invariant spacetime interval
Lindblad:       gamma  ×  t    =  K  =  invariant absorption dose
```

Gamma plays the algebraic role of c: the external parameter that
defines the time scale, appears in the invariant product, and cannot
be outrun by the system. The system cannot decohere faster than its
illumination allows, just as an object cannot travel faster than light.

Whether E = m × γ² holds as a mass-energy conversion is genuinely open.
The framework has γ where Einstein has c, t where Einstein has τ, and
K where Einstein has the spacetime interval. The algebra matches.
The physics is unproven.

---

## What survives (Tier classification)

1. **Surviving mode energy is center-localized** (Tier 2-3)
2. **K_death = 2.303 above Q*_gap(N)** (Tier 2)
3. **Immortal modes are massless** (Tier 3)
4. **The palindrome prevents total absorption** (Tier 2)
5. **Gamma plays the role of c** (Tier 2-3, algebraic parallel)
6. **E = mγ² as mass-energy conversion** (open, Tier 4-5)

---

## Reproduction

- Script: [`simulations/trapped_light_localization.py`](../simulations/trapped_light_localization.py)
- Output: [`simulations/results/trapped_light_localization.txt`](../simulations/results/trapped_light_localization.txt)
