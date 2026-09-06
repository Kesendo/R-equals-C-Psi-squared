# Where Trapped Light Concentrates

<!-- Keywords: mass trapped light center-localized, surviving mode energy,
K_death universal dose, immortal palindromic modes, gamma role of c,
entrance pupil window, cavity mass distribution, R=CPsi2 mass hypothesis -->

**Status:** Center-weighted late-time XY diagnostic observed for N=4,5 in the
declared sacrifice-zone runs. K_death = ln(10) = 2.303 above the coupling threshold
Q*_gap(N) (from D6; below it the slowest mortal mode is slower and the dose
grows). Gamma plays the algebraic
role of an inverse-time scale in K = γt. No mass, light, or energy-conversion
mechanism is established.
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
pupil), and the late-time transverse diagnostic might be centre-weighted. The
simulation tests that spatial diagnostic; it does not measure light or mass.

The declared runs give a centre/edge ratio of 1.3-1.4 for the late-time XY
diagnostic, consistent with
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994).
This finite numerical resemblance does not establish a match to stars, atoms,
or black holes.

---

## Result 1: Surviving mode energy is center-localized

The computed observable is the site-resolved late-time {X,Y} weight in the
declared simulation. Calling it mass, trapped light, or standing-wave energy
would require independent energetic and interference observables that are not
present here.

Under the sacrifice zone at t = 20:

N = 4:
| Site | γ | XY energy (surviving) |
|------|------|----------------------|
| 0 (edge) | 0.197 | 0.209 |
| 1 | 0.001 | 0.290 |
| 2 | 0.001 | 0.288 |
| 3 | 0.001 | 0.213 |

N = 5:
| Site | γ | XY energy (surviving) |
|------|------|----------------------|
| 0 (edge) | 0.246 | 0.165 |
| 1 | 0.001 | 0.221 |
| 2 | 0.001 | 0.220 |
| 3 | 0.001 | 0.220 |
| 4 | 0.001 | 0.174 |

Center/edge ratio: 1.37 (N=4), 1.30 (N=5), taken against the mean of the
two edges. Both edges hold less surviving mode energy than the centre, and
the entrance pupil is the lower of the two: 0.209 against site 3's 0.213 at
N=4, 0.165 against site 4's 0.174 at N=5. The sampled XY weight is lower
at both edges than in the centre. The entrance-pupil language is an analogy.

Under uniform gamma in these runs, the sampled late-time XY diagnostic is
numerically absent at t=20. The comparison shows that the chosen sacrifice
profile preserves more of that diagnostic; it is not a general necessity
theorem.

The calculation establishes only the two finite centre/edge ratios above. It
does not infer where mass forms or a causal relation to stellar, atomic, or
black-hole structure.

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
Complete absorption is impossible: the Liouvillian always has a kernel.
What the cavity always retains is lens content, never light.

*Source: [`trapped_light_localization.py`](../simulations/trapped_light_localization.py). K_death = ln(10) = 2.303 (proven from D6: rate_min = 2γ, so K = γ·ln(100)/(2γ) = ln(10)).*

---

## Result 3: Zero-rate modes are structural, not an energy census

The N+1 immortal modes (zero absorption rate) are the operators that
describe the system's overall state and the conserved quantities in each
magnetization sector. They never absorb because they contain no
transverse (X,Y) content: they are pure structure, invisible to the
declared dephasing channel. This says nothing about their mass or energy.

The nontrivial kernel blocks complete convergence to zero. That fact follows
from stationary modes, not protection of every paired mode or a physical
absorption-energy balance. The separate [Thermal Blackbody](THERMAL_BLACKBODY.md)
census is algebraic, not thermodynamic.

---

## Result 4: Gamma supplies an inverse-time scale

As a dimensional analogy only, special relativity supplies a velocity scale
and proper time whose product has dimensions of length. No identification of
those quantities with the open-system variables below is established.

In the Lindblad cavity:
- gamma is a chosen decay rate that sets the model's dissipative timescale.
  The [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) certifies
  openness, not an external origin; IBM dephasing can have several microscopic
  contributions and is not identified here with photon shot noise.
- t is the evolution parameter
- γ × t = K = invariant absorption dose
  ([F14](../docs/ANALYTICAL_FORMULAS.md), proven)

```
Relativity:     c     ×  tau  =  invariant spacetime interval
Lindblad:       gamma  ×  t    =  K  =  invariant absorption dose
```

Gamma defines a decay timescale and appears in the dimensionless product K.
There is no speed-limit theorem here: Liouvillian decay rates also depend on
the channel, system size, Hamiltonian and mode content.

No computed quantity in this experiment supports `E=m*gamma^2` or an
energy-conversion causality. The dimensional resemblance to `c*tau` is not a
physical equivalence.

---

## What survives (Tier classification)

1. **Surviving mode energy is center-localized** (Tier 2-3)
2. **K_death = 2.303 above Q*_gap(N)** (Tier 2)
3. **The zero-rate sector is insensitive to the declared dephasing channel** (Tier 1 within scope)
4. **The kernel prevents convergence of the full operator space to zero** (Tier 1 within scope)
5. **Gamma has units of inverse time and forms K=gamma*t** (dimensional identity)
6. **No mass-energy conversion is established**

---

## Reproduction

- Script: [`simulations/trapped_light_localization.py`](../simulations/trapped_light_localization.py)
- Output: [`simulations/results/trapped_light_localization.txt`](../simulations/results/trapped_light_localization.txt)
