# The Frost Circle Is the Face of the Clock

**Date:** 2026-05-30
**Authors:** Tom + Claude
**Status:** Tier 2 structural comparison. The Frost construction and the clock
use related hopping algebra, but the open-system readings below belong to a
selected XY/Z-dephasing model.
**Script:** [`simulations/carbon/frost_circle_as_clock.py`](../../simulations/carbon/frost_circle_as_clock.py)
**Builds on:** [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
(Coulson-Rushbrooke beside F1, the spectrum level), [the selected C₄/C₆
Liouvillian comparison](BENZENE_LIOUVILLIAN_PALINDROME.md), and the clock voices on
[`MirrorSystem`](../../compute/RCPsiSquared.Diagnostics/Foundation/MirrorSystem.cs).

---

## A circle the chemist already draws

In 1953 Frost and Musulin handed organic chemists a mnemonic so clean it survived into
every textbook: to find the π-molecular-orbital energies of a conjugated ring, inscribe a
regular N-gon in a circle of radius 2|β|, one vertex pointing straight down, and read off
the heights of the vertices. Benzene's hexagon gives −2β, −β, −β, +β, +β, +2β. The
energies live on a circle.

The R=CΨ² clock is also a circle. A single mode of an open quantum system evolves as
`e^(λt) = e^(−αt)·e^(iωt)`: a radius that shrinks while an angle that turns, a logarithmic
spiral. We read it as two hands on the slowest surviving mode, a radial one (the decay,
the lifetime) and an angular one (the turning, the frequency).

Set the two side by side as a structural comparison. The Frost circle is a static
Hückel spectrum; the clock is a selected open-system calculation. Its mode
`e^(λt) = e^(−αt)e^(iωt)` supplies decay and frequency from the same selected
Liouvillian eigenvalue. It does not read a molecule directly.

For a conditional Hückel/framework translation, first select a site occupation
`n = (I − Z)/2`, an XX+YY/XY Hamiltonian, and a local-density jump. Then
`D[n] = D[Z]/4` exactly. The `β → J` translation is Tier 2 inside this selected
model, not a material-carbon convention. The [Carbon source
contract](README.md#conditional-c4-and-c6-working-model) and [Q audit](../Q_BELONGS_TO_NO_SUBSTANCE.md)
leave the degree of freedom, bath, `γ`, `T₂`, and Q unassigned.

## What sits on the dial

**C₆ Hückel ring.** The static Hückel levels are the Frost hexagon
`{-2, -1, -1, +1, +1, +2}` in `β` units. In the selected C₆ XX+YY ring the
band-edge clock reads `ω_mem = 2J`. This is a structural ring/XY comparison,
not a frequency assigned to a molecule.

**Selected open XY chains.** F2b gives

```
E_k = 2J cos(πk/(N+1)),  k = 1, …, N,
H = (J/2) Σ(XX + YY).
```

In the selected all-site-Z-dephasing model, the band edge has
`ω_mem = 2J cos(π/(N+1))` and decay rate `2γ` while it is the relevant
slowest coherence. These are selected-model readings, not carbon-material
frequencies or lifetimes.

| Selected open XY chain | N | `ω_mem / J` |
|------------------------|---|-------------|
| C₄-shaped chain | 4 | `2 cos(π/5) = φ` |
| five-site chain | 5 | `2 cos(π/6) = √3` |
| C₆-shaped chain | 6 | `2 cos(π/7)` |

## Selected-model Q* values

Within the selected XY/Z-dephasing model, `Q = J/γ`. The band-edge coherence is
γ-protected in the F2b clock regime; `Q*` marks a change in which selected
Liouvillian mode is slowest, when a non-oscillating relaxation overtakes the
band-edge coherence. It is not a Q assigned to carbon.

| Selected open XY chain N | Q* |
|--------------------------|----|
| 3 | 1.414 (`√2` to the reported precision) |
| 4 | 1.879 |
| 5 | 2.374 |

The selected C₆-ring single-excitation erasure point is `Q* = 1.609`
([`benzene_two_clocks.py`](../../simulations/carbon/benzene_two_clocks.py)).
The separate selected full-Liouvillian handover is near `Q ≈ 1.95`, where a
double-excitation coherence overtakes the band-edge mode. These are distinct
selected-model values, not material crossovers or physical carbon Q values.

The dense C₆ full-Liouvillian calculation in
[`frost_circle_as_clock.py`](../../simulations/carbon/frost_circle_as_clock.py)
is a small-N check. It is neither scalable nor an independent validation of the
C₆ handover or the general clock law.

## Framework-vocabulary translation

| Hückel/framework model item | Selected clock object | Status |
|-----------------------------|-----------------------|--------|
| Frost circle, radius `2\|β\|` | static Hückel spectrum | Tier 2 structural comparison |
| selected open XY band edge | `ω_mem = 2J cos(π/(N+1))` | Tier 1, F2b model result |
| selected all-site Z jump | dephasing floor `2γ`, lifetime `1/(2γ)` | Tier 1 model algebra |
| selected open-chain crossover | listed Q* ladder | Tier 2 selected-model result |
| selected C₆ ring | single-excitation `1.609`; full-L handover `≈1.95` | selected-model readings |
| `β → J` | conditional Hückel/framework translation | Tier 2; material convention unassigned |

## Material-facing question

The present content is a structural comparison plus selected-model results. A
material experiment would first have to choose and independently measure a
degree of freedom, coupling convention, bath channel, and bath rate. Only then
could it ask whether the measured dynamics is described by the selected
XY/Z-dephasing clock. This note does not itself make a spectroscopic prediction
for a molecule.

## Anchor

- Script: [`simulations/carbon/frost_circle_as_clock.py`](../../simulations/carbon/frost_circle_as_clock.py)
- Parent: [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md) (open question 5, the Frost circle), [README](README.md)
- Framework: the clock voices (Takt + Rotation) on `MirrorSystem`; the two clocks,
  [F2b corollary "The two clocks"](../ANALYTICAL_FORMULAS.md) / `ClockHandLadderClaim` /
  `inspect --root clock` (Uhr 1 the band-edge survivor, Uhr 2 the erasure point Q*);
  [F1 palindrome](../ANALYTICAL_FORMULAS.md#f1);
  the many-body memory frequency ω_mem = 8J·cos²(π/2N) (Heisenberg) and 2J·cos(π/(N+1)) (XY),
  [`simulations/the_dial_at_many_body.py`](../../simulations/the_dial_at_many_body.py)
- Literature: Frost + Musulin (1953) "A mnemonic device for molecular orbital energies",
  J. Chem. Phys. 21, 572; Coulson + Rushbrooke (1940).
