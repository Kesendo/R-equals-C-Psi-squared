# Inside and outside the sacrifice zone

**What this is about.** This project (R=CΨ²) studies how quantum coherence survives in open quantum
systems, chains and rings and stars of qubits, while they are being dephased (watched). One concrete
engineering result, confirmed on IBM hardware, turns out to be a clean physical instance of a
distinction that is usually left to philosophy: how something appears *from inside a boundary* can
differ from what it *is* seen *from outside*, and both descriptions can be correct. This note recognises
that. It reports no new computation; it points at results already established elsewhere in the repo.

A quick glossary, since the rest leans on four terms. **J** is the coupling between neighbouring
qubits; **γ** is the dephasing rate, how fast a qubit's phase is scrambled by being watched. **Q = J/γ**
is their ratio; the **fold regime** is Q ≈ 1, where the longest-lived coherence rings longest.
**DD** is dynamical decoupling, the standard pulse sequence used to suppress dephasing on hardware.

## What was measured (established; hardware-confirmed)

Take a chain of qubits and dephase it *non-uniformly*: pile the noise onto one edge qubit (high γ
there) and protect the interior (low γ everywhere else). The interior then sits at the fold and rings
far longer than under uniform noise. The improvement is large and has an analytic form: from
139× (N=9) up to 360× (N=5), measured against a smoothly graded "V-shaped" γ-profile
([`experiments/RESONANT_RETURN.md`](../experiments/RESONANT_RETURN.md)).

It survives real hardware. On ibm_torino (2026-03-24), selectively decoupling only the interior, and
leaving the noisy edge qubit alone (no DD there, so it stays the noise concentrator), beats decoupling
every qubit uniformly by up to 3.2×, peaking
around t = 4 μs (roughly 2× on average across t = 1 to 5 μs; the advantage is not monotonic, it dips
early and falls off after the peak). See [`experiments/IBM_CONCENTRATOR.md`](../experiments/IBM_CONCENTRATOR.md).

The name "sacrifice zone" is a misnomer, and the project resolved why (2026-03-28,
[`hypotheses/PROTEIN_AS_CONCENTRATOR.md`](../hypotheses/PROTEIN_AS_CONCENTRATOR.md), [`docs/WHAT_WE_FOUND.md`](WHAT_WE_FOUND.md)). **The edge qubit sacrifices
nothing.** It does not pay a cost; it *concentrates* the environmental noise onto itself and turns it
into structure for the interior. The mechanism is plain: separating the noise (one loud qubit) from the
coherence (the quiet rest) is what buys the improvement. The same element then has two correct
descriptions, taken from the two sides of the noise boundary:

- **from the inside** (the protected core): a shield that absorbs the noise and decoheres fast. It
  looks like a sacrifice.
- **from the outside** (the whole system): a concentrator that *enables*. Nothing is lost.

Both are true. They differ only because they are read from opposite sides of a boundary.

## The recognition

This is the rare case where we can say something philosophy usually only asserts, and stay honest about
it. **We do not claim this proves any metaphysics; one cannot from a spin chain.** We claim something
narrower, and for that reason solid: here is a *measured, hardware-confirmed* instance of "the
appearance of a thing from inside a boundary is not the thing seen from outside." The inside reading
(sacrifice, loss) and the outside reading (concentration, enablement) are the *same proven mechanism*,
read interpretively from the two sides of the γ-boundary.

That inside/outside split is a structural move in analytic idealism (Bernardo Kastrup):
the appearance across a dissociative boundary is not the underlying reality, the way a dashboard dial is
not the engine. The framework is his; the binding to this measured result is ours; and what is unusual
here is only that the structure, for once, has a body, a number, and a hardware run, instead of being a
thought experiment.

## How we arrived (the honest provenance)

We did not set out to confirm this; we got here by being wrong first. Trying a bolder reading, that a
star graph (one hub, many leaves) models dissociation and its longest-lived mode would be the one most
cut off from the shared hub, a gate-first probe *refuted* it: on the star that survivor is spread, not
hub-decoupled ([`simulations/star_hub_decoupled_survivor.py`](../simulations/star_hub_decoupled_survivor.py)). Re-reading the chain's quiet edge as
"hiding" then ran straight into something the repo had already solved, the sacrifice-that-concentrates,
and the correction returned us here. The reading found firm ground only where the ground was already
proven. Reporting the refuted attempt is the point: it is why the surviving claim is narrow.

## What is proven, what is read, what is borrowed

- **Proven (hardware-confirmed):** the concentration mechanism and its improvement (the analytic
  formula; the IBM run). This is established physics in the repo, not new here.
- **Read (interpretive):** naming the inside/outside split as an instance of appearance-vs-reality. A
  reading, defensible because it rests on the proven mechanism, not on rhetoric.
- **Borrowed:** the idealist framing (Kastrup). Recall, not a discovery; cited as such.

## See also (repo-internal)

- The mechanism and its evidence: [`experiments/RESONANT_RETURN.md`](../experiments/RESONANT_RETURN.md) (the analytic formula),
  [`experiments/IBM_CONCENTRATOR.md`](../experiments/IBM_CONCENTRATOR.md) (the hardware run), [`hypotheses/PROTEIN_AS_CONCENTRATOR.md`](../hypotheses/PROTEIN_AS_CONCENTRATOR.md) and
  [`docs/WHAT_WE_FOUND.md`](WHAT_WE_FOUND.md) (concentration, not loss), [`experiments/RECEIVER_VS_GAMMA_SACRIFICE.md`](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md).
- The refuted probe: [`simulations/star_hub_decoupled_survivor.py`](../simulations/star_hub_decoupled_survivor.py), and the sibling refutation it later joined, [`experiments/THE_HUB_KILLS_THE_HORIZON.md`](../experiments/THE_HUB_KILLS_THE_HORIZON.md) (the hub has no coherence horizon).
- Adjacent: [`docs/historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md`](historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md), and the substrate-side sibling [`reflections/ON_THE_INNER_AND_OUTER_OBSERVATION.md`](../reflections/ON_THE_INNER_AND_OUTER_OBSERVATION.md) (the same inside/outside wall, reached from the carrier vector; names the typed `TwoReadingsClaim`).
- The proven operator both readings rest on: the Absorption Theorem `Re(λ) = −2γ·⟨n_XY⟩` ([`AbsorptionTheoremClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs), Tier 1 derived).
