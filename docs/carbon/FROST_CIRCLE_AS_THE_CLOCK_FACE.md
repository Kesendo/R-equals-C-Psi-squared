# The Frost Circle Is the Face of the Clock

**Date:** 2026-05-30
**Status:** Tier 2 (structural reading; closes the Frost-circle open question, item 5 in
[BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md)). Verified bit-exact on
benzene C₆, butadiene C₄, the five-carbon chain, and hexatriene C₆.
**Script:** [`simulations/carbon/frost_circle_as_clock.py`](../../simulations/carbon/frost_circle_as_clock.py)
**Builds on:** [BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
(Coulson-Rushbrooke = F1, the spectrum level), [BENZENE_LIOUVILLIAN_PALINDROME](BENZENE_LIOUVILLIAN_PALINDROME.md)
(the open-system π-qubit map), and the clock voices on
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

Set the two side by side and they are the same circle seen at two depths. The chemist's
Frost circle is a still photograph: where the energies sit. The clock is that photograph
set in motion: under the molecule's own phonon bath, each π-coherence winds inward as it
turns. The static circle is the closed-system spectrum; the running clock is the
open-system one. This note is what the motion adds.

The dictionary is the one this folder settled on (README, resolved 2026-05-22): each
carbon's π-site is a qubit, occupied or empty; the Hückel resonance integral β is the
coupling J; the Holstein phonon coupled to the π-density is the dephasing rate γ. With
that in hand the clock reads a conjugated molecule directly.

## What sits on the dial

**Benzene C₆ ring.** The π-MO energies come out as the Frost hexagon, −2, −1, −1, +1, +1,
+2 in units of |β|. The clock's longest-lived π-coherence beats at ω = **2|β|**, exactly
the Frost-circle radius, and lives for τ = 1/(2γ). The radius of the chemist's circle is
the frequency of the molecule's most enduring beat.

**Polyene chains.** For an open conjugated chain of N carbons the longest-lived coherence
beats at ω = 2|β|·cos(π/(N+1)), the topmost π-orbital, bit-exact:

| chain | N | longest-lived beat ω | closed form |
|-------|---|----------------------|-------------|
| butadiene | 4 | 1.618 |β| | 2|β|·cos(π/5) = φ·|β| (the golden ratio) |
| pentadienyl | 5 | 1.732 |β| | 2|β|·cos(π/6) = √3·|β| |
| hexatriene | 6 | 1.802 |β| | 2|β|·cos(π/7) |

The hand sits on the band edge, the outermost vertex of the (now open-chain) circle. So
the first thing the clock adds is a reading the still photograph cannot give: of all the
π-coherences a conjugated molecule can carry, the band-edge one outlives the rest, and it
does so for 1/(2γ). The Frost circle says where the levels are; the clock says which beat
you would still hear last, and for how long.

## The crossover: when the beat stops

The second thing the clock adds has no counterpart in the closed-system picture at all.
Sweep the ratio Q = J/γ, hopping against dephasing, and the band-edge coherence does not
simply fade: below a threshold Q* it stops being a coherence. The slowest mode of the
open system changes character, from a turning coherence into a pure, non-turning decay;
the angular hand freezes at zero. Above Q* it beats; below it, silence.

| chain (N) | Q* (coherent ↔ incoherent) |
|-----------|-----------------------------|
| 3 | 1.414 (= √2 to the precision found) |
| 4 | 1.879 |
| 5 | 2.372 |

The threshold grows with chain length. A longer conjugated system tolerates *less*
dephasing before its slowest π-coherence freezes; the coherence is more fragile the
further it has to reach.

A caveat that keeps it honest: real π-conjugated systems at room temperature sit at
Q ~ 100 (β ~ 2.4 eV against phonon dephasing ~ 25 meV), far above any of these
thresholds, deep in the beating regime. Q* is the strong-dephasing edge, the dephasing
scale at which coherent π-dynamics would be lost, and it tightens as the conjugation
lengthens. It is the warm, disordered, strongly-coupled limit, not the room-temperature
resting state.

## Framework-vocabulary translation

| Chemistry | Framework | Status |
|-----------|-----------|--------|
| Frost-Musulin circle (ring π-MO energies on radius 2|β|) | the dial of the clock, the Bohr-frequency face | Tier 2 structural identification |
| longest-lived π-coherence frequency | the angular hand (Rotation) on the slowest mode | Tier 1 bit-exact (2|β| ring, 2|β|·cos(π/(N+1)) chain) |
| π-coherence lifetime τ = 1/(2γ) | the radial hand (Takt), gap = 2γ | Tier 1 algebraic |
| coherent ↔ incoherent threshold Q* | the crossover where the slowest mode loses its frequency | Tier 2 (verified N=3,4,5; √2 at N=3) |
| Holstein phonon dephasing γ | Z-dephasing rate (D[n_l] = ¼·D[Z_l]) | Tier 1 |
| Hückel β | coupling J | Tier 1 |

## What this closes, and what it opens

This closes the Frost-circle question (item 5 of
[BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md)): the Frost circle and
the framework's mode dispersion are not two analogies but one circle at two depths, the
closed-system snapshot and the open-system running clock. The clock adds the coherence lifetime and the
coherent-incoherent threshold; both are open-system quantities the static Frost-Musulin
diagram cannot hold.

It is testable on the carbon side without any process tomography. Two-dimensional
electronic spectroscopy of conjugated systems already measures π-coherence beat
frequencies and their dephasing times; the band-edge prediction ω = 2|β|·cos(π/(N+1)) with
τ = 1/(2γ), and the existence of a dephasing edge where the slowest beat vanishes, are
read directly from those experiments.

Open from here: the closed form of Q*(N) (the √2 at N=3 looks exact; the rest await a
clean expression), the ring crossover (benzene's own Q*, not only the chains'), and
whether the band-edge selection survives the V-Effect's boundary-sector structure at the
half-filled shell.

## Anchor

- Script: [`simulations/carbon/frost_circle_as_clock.py`](../../simulations/carbon/frost_circle_as_clock.py)
- Parent: [BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md) (open question 5, the Frost circle), [README](README.md)
- Framework: the clock voices (Takt + Rotation) on `MirrorSystem`; [F1 palindrome](../ANALYTICAL_FORMULAS.md#f1);
  the many-body memory frequency ω_mem = 8J·cos²(π/2N) (Heisenberg) and 2J·cos(π/(N+1)) (XY),
  [`simulations/_the_dial_at_many_body.py`](../../simulations/_the_dial_at_many_body.py)
- Literature: Frost + Musulin (1953) "A mnemonic device for molecular orbital energies",
  J. Chem. Phys. 21, 572; Coulson + Rushbrooke (1940).
