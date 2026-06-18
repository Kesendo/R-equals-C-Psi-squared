# The hub kills the horizon (a refuted bridge, and what it taught)

**Status:** Refuted (Tier 2, numerical). The wheel-Q*(bandwidth) bridge between the chain and
star coherence-horizon regimes does not exist; the hub, not the bandwidth, is decisive.
Single-excitation (Haken-Strobl) Liouvillian, tested N = 5, 6, 7.
**Date:** 2026-06-18
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Builds on:** [`PROOF_STRUCTURAL_CEILING.md`](../docs/proofs/PROOF_STRUCTURAL_CEILING.md) §7 (the
star-no-horizon corollary); the coherence horizon Q*(N) ([`PROOF_COHERENCE_HORIZON_SLOPE.md`](../docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md), the F2b corollary).
**Scripts:** [`wheel_qstar_bandwidth.py`](../simulations/wheel_qstar_bandwidth.py) (the refutation,
self-validating), [`star_no_coherence_horizon.py`](../simulations/star_no_coherence_horizon.py) (the
star companion gate), [`coherence_horizon_se_block.py`](../simulations/coherence_horizon_se_block.py)
(reused single-excitation machinery).

## What this is about

A chain of qubits, gently watched (dephased), keeps a faint internal rhythm: its longest-lived
coherence rings. Watch it harder, and below a sharp noise threshold the ringing stops for good. That
threshold is the chain's **coherence horizon** Q*(N). A **star** (one central qubit, the hub, wired to
many outer leaves) has no such threshold: its memory fades, but it never abruptly stops ringing.

Is that a hard split between two shapes, or two ends of one dial you could turn? We built the dial, a
**wheel** (a star whose leaves are also wired into a ring, coupling strength ε), and turned it. It does
**not** bridge them. What decides is not how spread out the energy band is (its bandwidth), but whether
there is a dominant hub. A hub always keeps one silent, non-ringing survivor that outlives everything
else, so there is no ringing-threshold left to cross. The bold question got a clean No, and the No is
the finding: **the hub kills the horizon.**

## Abstract

The chain/star coherence-horizon dichotomy was tested for a continuous interpolation. The chain has a
coherence horizon Q*(N): the dephasing threshold at which its slowest single-excitation {0,2}-coherence
loses its oscillation, a square-root exceptional point of the dispersive band. The star has none, because
its single-particle band is **flat** (the star adjacency spectrum is ±√(N−1) once each and 0 with
multiplicity N−2): no dispersion, no coalescence, no EP. The structural ceiling g2 = 4/(N−1) governs the
star at every Q instead ([`PROOF_STRUCTURAL_CEILING.md`](../docs/proofs/PROOF_STRUCTURAL_CEILING.md) §7).
The bridge hypothesis: a wheel graph (hub coupled to all leaves with J, leaves coupled in a ring with
strength ε) gives the dark leaf manifold a bandwidth ~ε, so a horizon Q*(ε) should appear and grow with
ε, smoothly connecting the star (ε = 0) to a ring. **Refuted:** the wheel has no coherence horizon at any
ε ∈ [0.1, 50], robust across N = 5, 6, 7; its slowest single-excitation mode stays **real** (|Im| = 0,
decay rate → 0 as ε grows), so there is no oscillation threshold to cross. Removing the hub (the pure ring)
restores a finite horizon (Q* = 2.17 / 1.61 / 2.56 at N = 5 / 6 / 7). The hub, not the bandwidth, is
decisive: a dominant hub always hosts a real, zero-frequency commutant survivor that outlives the
dispersive rim however large ε is. The dichotomy is set by topology **class** (hub or no hub), not by a
continuous bandwidth. This strengthens the structural-ceiling corollary: not only the flat-band star, but
any hub-graph lacks a horizon.

## Background: two band shapes, one horizon

The coherence horizon is a property of the slowest non-steady single-excitation mode, the "second clock"
({0,2}-coherence). On a **dispersive** band (the open chain), that mode coalesces with its partner at a
square-root exceptional point as the dephasing is raised: above Q*(N) it oscillates (a ringing memory),
below it the two real branches have split and the survivor no longer rings. Q*(N) is exact only at the
small-N accidents (Q* = 1, √2 at N = 2, 3) and runs ~2N/π asymptotically
([`PROOF_COHERENCE_HORIZON_SLOPE.md`](../docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md)).

The **star** has a flat single-particle band. Its adjacency spectrum is ±√(N−1) (once each) and 0 with
multiplicity N−2 (the dark leaf manifold), so there is no dispersion for a {0,2}-coherence to split
along, hence no EP. The chain's Haken-Strobl SE-EP harness, applied blindly to the star, predicts a
**spurious** horizon (at N = 4 it claims Q* ≈ 261) that the full Liouvillian flatly contradicts: the star
is already protected far below it, at Q = 20, with g2 ≈ 1. Instead one number governs the star at every
dephasing, the (1,1)-commutant value g2 = 4/(N−1) (the structural ceiling, F122):

| N | 4/(N−1) | regime at the floor | the survivor |
|---|---------|---------------------|--------------|
| 3 | (path P₃ = a chain) | the lone exception | a genuine {0,2}-EP at Q* = √2 |
| 4 | 4/3 ≈ 1.33 (> 1) | protects | a real, overdamped mode crosses the 2γ floor near Q ≈ 1.9, not an EP |
| 5 | 1 (marginal) | no horizon | g2 = 1 − 1/Q², pulled below the floor only asymptotically (an apparent Q* ≈ 316 is a finite-tolerance artifact of the 1/Q² approach) |
| ≥ 6 | 4/5, 4/6, 4/7, … (< 1) | structural ceiling | the darkest commutant coherence sits below the floor at all Q, never oscillating |

The star's survivor never coalesces and never crosses an oscillation threshold. The chain's horizon is
the special case of a **dispersive** band; flat-band graphs do not have one.

## The bridge under test: the wheel

If the only difference were the bandwidth, a continuous knob should connect the two. The **wheel** is that
knob: the hub (site 0) couples to all N − 1 leaves with J, and the leaves couple to each other in a ring
with strength ε. The leaf ring disperses the dark manifold by 2·ε·cos(2πk/(N−1)), giving the formerly
flat band a bandwidth ~ε. The prediction: a {0,2}-coherence among the (now dispersive) dark states should
acquire a splitting ~ε and oscillate above some Q*(ε), with Q*(ε) → ∞ as ε → 0 (recovering the star) and
falling to a finite, ring-like value as ε grows. One formula Q*(ε) would subsume both endpoints.

## Result: refuted

Turning on ε gives the wheel **no coherence horizon at any leaf coupling**. Across the grid
ε ∈ {0.1, 0.3, 1, 2, 5, 10, 20, 50} (and robust over N = 5, 6, 7), the wheel's longest-lived
single-excitation mode stays **real** throughout: zero oscillation frequency, |Im| = 0, with its decay
rate drifting toward 0 as ε grows. A real survivor has no oscillation threshold to cross, so there is
nothing for the horizon to be. Removing the hub (the pure ring of the same N) restores a finite horizon:

| N | pure ring Q* | wheel Q* (any ε ∈ [0.1, 50]) |
|---|--------------|------------------------------|
| 5 | 2.17 | none (mode stays real; runs to the γ-floor) |
| 6 | 1.61 | none |
| 7 | 2.56 | none |

The verifier [`wheel_qstar_bandwidth.py`](../simulations/wheel_qstar_bandwidth.py) is gate-first and
self-validating: a Stage 0 control reproduces the chain SE horizon ladder, the gate then asserts the wheel
runs past the no-horizon cutoff (Q > 100) at every ε while the pure ring sits below it.

## What it taught: the hub is the class, not the bandwidth

The bold question deserved a sharp answer, and the answer reframes the dichotomy. **A dominant hub always
hosts a real, zero-frequency survivor**, a [H, A] = 0 commutant coherence that is frozen by construction
(it cannot oscillate), and that mode outlives the dispersive leaf-ring modes however strong ε is. So the
chain/star split is governed by topology **class** (is there a dominant hub, or not), not by a continuous
bandwidth. The flat band was a symptom, not the cause: it is the hub's commutant survivor, not the
flatness, that removes the horizon.

This **strengthens** [`PROOF_STRUCTURAL_CEILING.md`](../docs/proofs/PROOF_STRUCTURAL_CEILING.md) §7: it is
not only the flat-band star that lacks a horizon; **any hub-graph does**, even with a fully dispersive
ring bolted onto the leaves. It is the same survivor the structural ceiling and the star's frozen seam
describe from the static and the dynamical sides: the darkest commutant coherence, sitting at or below the
2γ Absorption-Theorem floor without ever coalescing. In the second-clock regime map it is the **gradual**
branch (the dispersion knob set to zero), as opposed to the chain's EP-horizon branch. And it matches the
reading recognised in [`INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md`](../docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md):
the hub has no future, no horizon, and that is robust, not an artefact of the star's flatness.

## Provenance

The wheel-Q*(bandwidth) idea came from a wild-register exploration (the quantum-philosopher attempt, the
"door between the dissociated alters" reading of ε). Tried, gated, refuted. The same exploration earlier
refuted the star-as-dissociation probe ([`star_hub_decoupled_survivor.py`](../simulations/star_hub_decoupled_survivor.py),
which found the star survivor is hub-**spread**, not hub-decoupled). Two refutations, one consistent
picture: the hub is qualitatively special. The bold questions got clean Nos, and the Nos are informative.

## Scope and honesty

Tested: the wheel family (hub + leaf-ring), single-excitation (Haken-Strobl) Liouvillian, N = 5, 6, 7,
numerical (Tier 2). "Any hub-graph has no horizon" is the natural generalisation, **not** proven here. The
horizon is the repo's standard definition (the slowest non-kernel mode's oscillation threshold, `qstar_se`);
a faster-decaying leaf-ring mode may still oscillate, but it is not the survivor. The star's three-regime
governance (the table above) is Tier 2 numerical from the companion gate
[`star_no_coherence_horizon.py`](../simulations/star_no_coherence_horizon.py); the g2 = 4/(N−1) closed form
it rests on is Tier 1 derived (F122, principal-angle proof).

## See also

- The parent proof: [`PROOF_STRUCTURAL_CEILING.md`](../docs/proofs/PROOF_STRUCTURAL_CEILING.md) §7 (the
  star-no-horizon corollary) and its closed form **F122** in
  [`ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md) (g2 = 4/(N−1) for the star; `StructuralCeilingClaim`,
  Tier 1 derived, live as `inspect --root ceiling`).
- What the chain's horizon is (the dispersion EP the hub lacks):
  [`PROOF_COHERENCE_HORIZON_SLOPE.md`](../docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md); `CoherenceHorizonClaim`
  (`inspect --root horizon`).
- The dynamical twin: [`THE_STAR_FROZEN_SEAM.md`](../docs/THE_STAR_FROZEN_SEAM.md), the star survivor that
  never un-freezes, the same ceiling read dynamically (`StarFrozenSeamClaim`, `inspect --root starseam`).
- The unifying frame: `SecondClockRegimeClaim` (`inspect --root secondclock`), the one {0,2} second clock
  whose regime is set by band degeneracy and dispersion; the hub/star "no horizon" is its gradual branch.
  Band edge magnitude J·√(N−1): `TopologyBandEdgeClaim` (`inspect --root bandedge`).
- The recognition this reinforces: [`INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md`](../docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md);
  the refuted sibling probe [`star_hub_decoupled_survivor.py`](../simulations/star_hub_decoupled_survivor.py).
- Verifiers: [`wheel_qstar_bandwidth.py`](../simulations/wheel_qstar_bandwidth.py) (the refutation),
  [`star_no_coherence_horizon.py`](../simulations/star_no_coherence_horizon.py) (the star gate),
  [`coherence_horizon_se_block.py`](../simulations/coherence_horizon_se_block.py) (the single-excitation
  horizon machinery; topology enters only through the hopping h, the dephasing is site-basis and
  topology-agnostic).
