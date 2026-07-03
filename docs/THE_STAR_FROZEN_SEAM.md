# The star's frozen survivor is the structural ceiling, read dynamically

**Status:** Tier 1 candidate (gate-verified N = 4..8, numerical). The star's longest-lived survivor
never un-freezes for N ≥ 5; the frozen seam is the structural ceiling g2 = 4/(N−1) ≤ 1, read dynamically.
**Date:** 2026-06-18
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Typed claim:** [`StarFrozenSeamClaim`](../compute/RCPsiSquared.Core/Symmetry/StarFrozenSeamClaim.cs)
(Tier 1 candidate) + the live witness `StarFrozenSeamWitness` (`inspect --root starseam`).
**Verifier:** [`star_frozen_seam.py`](../simulations/star_frozen_seam.py) (self-validating, N = 4..8).
**Builds on:** [`PROOF_STRUCTURAL_CEILING.md`](proofs/PROOF_STRUCTURAL_CEILING.md) §7 (the star has no
coherence horizon, the structural ceiling g2 = 4/(N−1)); `SecondClockRegimeClaim`, `CoherenceHorizonClaim`.

## What this is about

Watch the longest-lived coherence of a dephased spin network as you weaken the dephasing (raise
`Q = J/γ`, the coupling-to-observation ratio). Does that survivor stay *frozen* (a pure decay, no
oscillation, `|Im λ| = 0`) or does it *un-freeze* (start ringing)? The answer separates the topologies,
and it gives the star its own place in a trichotomy the typed layer had named only for the chain and the
ring. (`|Im λ|` is the imaginary part of the survivor's Liouvillian eigenvalue, its oscillation
frequency. `g2 = 4/(N−1)` is the star's structural ceiling, the darkest reachable decay-gap.)

## Abstract

A dephased spin network's longest-lived coherence is either frozen (overdamped, `|Im λ| = 0`) or
oscillating, and which one it is as `Q = J/γ` grows separates the topologies. The chain and ring both
*un-freeze*: the chain's `(p,p)` interior survivor yields to an oscillating `(0,1)` band edge above its
coherence horizon `Q*(N)` (a square-root exceptional point of the dispersive band), and the ring's
`(2,2)` frozen level crossing yields to its band edge above the handover. The **star does not.** Its
survivor is the darkest `[H,A] = 0` commutant `(1,1)` coherence, which lies in the `ad_H` kernel
(`−i[H,ρ] = 0`) and so cannot oscillate; it is frozen at *every* `Q`. But it is the survivor (the slowest
mode) only when it undercuts the `−2γ` Absorption floor, i.e. exactly when the structural ceiling
`g2 = 4/(N−1) ≤ 1`, which holds for `N ≥ 5`. `N = 4` (`g2 = 4/3 > 1`) is the outlier: the commutant mode
is not dark enough, an oscillating band-edge mode wins, and the star un-freezes (max `|Im| = √3·J ≈ 1.73`,
the `(0,1)` band edge oscillating at the star's single-excitation top `√(N−1)·J`). So the star's frozen seam *is* the structural ceiling read
dynamically, and it is the third member of the trichotomy chain (un-freezing square-root exceptional point, SE-EP) / ring (level
crossing) / star (commutant). Gate-verified `N = 4..8`; Tier 1 candidate.

## Model scope: the XY ceiling vs the Heisenberg survivor

The darkness value `g2 = 4/(N−1)` here is the **XY network** (hopping-only) ceiling, by design:
`StructuralCeilingClaim` is "the high-Q gap rate of an XY network," and both `StructuralCeilingWitness`
and `StarFrozenSeamWitness` build the star's XY hopping Hamiltonian with no ZZ diagonal. For the XY star
this is exactly the survivor darkness, `⟨n_XY⟩(Q→∞) = 4/(N−1)`.

The project's canonical survivor model is **Heisenberg** (`XX+YY+ZZ`: the chain survivor,
[`SurvivalIncompletenessMirrorClaim`](../compute/RCPsiSquared.Diagnostics/Foundation/SurvivalIncompletenessMirrorClaim.cs),
[`SURVIVOR_FLIP_AND_REFLECTION_ODD`](../experiments/SURVIVOR_FLIP_AND_REFLECTION_ODD.md)). The Heisenberg
star survivor darkens *further*, to **`⟨n_XY⟩(Q→∞) = 4/N`** (verified N = 5..8, full-`4^N` cross-checked
at N = 6): the single-excitation ZZ potential (hub `−(N−1)`, leaves `N−3`) shifts which `ad_H`-kernel
commutant is darkest, from `4/(N−1)` (XY) to `4/N` (Heisenberg). Companion verifier:
[`star_survivor_heisenberg.py`](../simulations/star_survivor_heisenberg.py).

What is **model-robust** (holds in both XY and Heisenberg): the survivor is frozen (`|Im λ| = 0`) at
every Q for N ≥ 5; it is the `[H,ρ] = 0` commutant only in the high-Q *limit* (`‖[H,ρ]‖ ∝ 1/Q`, not zero
at finite Q); it sits at the `(1,1)/(N−1,N−1)` popcount boundary; the star has no coherence horizon. Only
the darkness *value* is model-specific. (The XY threshold logic, frozen iff `g2 = 4/(N−1) ≤ 1` with the N=4
outlier, is likewise XY-specific; the Heisenberg `4/N` analog is not characterized below N = 5 here.)

## The finding

A single `Q` does not tell the topologies apart: *below* its coherence horizon every topology's slowest
mode is overdamped (real, `|Im| = 0`). At `Q = 1.5` the chain `(2,2)`, ring `(2,2)`, and star `(1,1)`
survivors are all frozen. The signature is the `|Im|(Q)` *curve* (gate-verified, `simulations/star_frozen_seam.py`):

| topology | survivor as `Q` grows | un-freezes? |
|---|---|---|
| **chain** | `(p,p)` interior, frozen, then the `(0,1)` band edge takes over above `Q*(N)` | **yes** (the coherence horizon) |
| **ring** | `(2,2)` frozen seam, then the oscillating band edge above the handover | **yes** (the handover) |
| **star** | `(1,1)` boundary survivor, frozen at **every** `Q` (for `N ≥ 5`) | **no** |

So the star has its own frozen survivor, and it **never un-freezes**: the third case of the trichotomy.

## Why, and the threshold

The star's survivor is the darkest `[H,A] = 0` **commutant** coherence (the structural-ceiling mode).
A coherence that commutes with `H` has no coherent evolution, `−i[H,ρ] = 0`, so it cannot oscillate: it
is frozen by construction. But it is only the *survivor* (the slowest mode) when it is darker than the
`−2γ` Absorption floor, i.e. exactly when the ceiling `g2 = 4/(N−1) ≤ 1`:

- `N ≥ 5` (`g2 ≤ 1`): the commutant coherence undercuts the floor, is the survivor, and is frozen at all
  `Q`. Gate-verified frozen at `N = 5, 6, 7, 8` (max `|Im|` over the `Q` sweep ≈ 1e-15, machine zero).
- `N = 4` (`g2 = 4/3 > 1`): the commutant mode is not dark enough; an oscillating band-edge mode wins and
  the star **un-freezes** (max `|Im| = √3·J ≈ 1.73`). This is the known star outlier: the `(0,1)` band edge
  oscillating at the star's single-excitation top `√(N−1)·J`.

**So the star's frozen seam IS the structural ceiling `g2 = 4/(N−1) ≤ 1`, read dynamically.** The high-Q
ceiling and the all-Q frozenness of the survivor are the same fact: the commutant coherence sits below the
floor. This is the survivor-level reading of `PROOF_STRUCTURAL_CEILING.md` §7's "the star has no coherence
horizon": the global slowest mode never acquires a frequency.

## What this closes

The typed layer characterized *frozen vs oscillating* for the chain (oscillates, the SE-EP horizon) and
the ring (frozen `(2,2)` level crossing), but the star was absent from those statements (only "flat band /
no horizon / boundary" was typed). This names the star's case and ties it to an already-typed quantity
(the ceiling), so it is a consolidation, not a new mechanism. The star's freeze is a *different route* to
`|Im| = 0` than the ring's: the ring freezes by a **level crossing** (two real eigenvalues coincide), the
star by a **commutant** (the mode commutes with `H`). Both are frozen; the chain is neither.

## See also

- Verifier (self-validating): `simulations/star_frozen_seam.py` (the pin at `Q=1.5`; the `|Im|(Q)` sweep;
  the `g2 = 4/(N−1) ≤ 1` threshold across `N = 4..8`). Reuses `simulations/carbon/incompleteness_survivor.py`
  (the sector-projected survivor, validated bit-for-bit vs the full `4^N` `L` at `N = 4`).
- The ceiling: `docs/proofs/PROOF_STRUCTURAL_CEILING.md` (`g2 = 4/(N−1)`, §7 no-horizon); `StructuralCeilingClaim`,
  `SecondClockRegimeClaim` (the regimes), `CoherenceHorizonClaim` (the chain's un-freezing).
- Typed: `StarFrozenSeamClaim` (Core/Symmetry, Tier1Candidate, parents `StructuralCeilingClaim` +
  `SecondClockRegimeClaim`) + the live witness `StarFrozenSeamWitness` (`inspect --root starseam`, recomputes the
  `g2 = 4/(N−1)` commutant ceiling and the `|Im|(Q)` frozen/un-freeze sweep, N≤5).
- The seam picture this extends: `docs/STERILE_BIRTHCANAL_AND_THE_JUNCTION.md` (chain/ring frozen vs
  oscillating); the sibling negative result `experiments/THE_HUB_KILLS_THE_HORIZON.md` (a hub removes the
  horizon however dispersive the leaves).
- The unified view: [`docs/THE_TRICHOTOMY_SEEN.md`](THE_TRICHOTOMY_SEEN.md), the chain/ring/star trichotomy
  rendered as one figure + the live `inspect --root trichotomy` tree; this star frozen seam sits there
  alongside the chain's un-freeze and the ring's level crossing.
