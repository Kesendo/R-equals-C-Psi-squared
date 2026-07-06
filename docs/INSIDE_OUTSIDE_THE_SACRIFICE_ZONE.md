# Inside and outside the sacrifice zone

**What this is about.** This project (R=CΨ²) studies how quantum coherence survives in open quantum
systems, chains and rings and stars of qubits, while they are being dephased (watched). One concrete
engineering result, measured on IBM hardware, turns out to be a clean physical instance of a
distinction that is usually left to philosophy: how something appears *from inside a boundary* can
differ from what it *is* seen *from outside*, and both descriptions can be correct. This note recognises
that. It reports no new computation; it points at results already established elsewhere in the repo.

A quick glossary, since the rest leans on four terms. **J** is the coupling between neighbouring
qubits; **γ** is the dephasing rate, how fast a qubit's phase is scrambled by being watched. **Q = J/γ**
is their ratio; the **fold regime** is Q ≈ 1, where the longest-lived coherence rings longest.
**DD** is dynamical decoupling, the standard pulse sequence used to suppress dephasing on hardware.

## What was measured (simulation + a hardware improvement; the mechanism attribution open)

Take a chain of qubits and dephase it *non-uniformly*: pile the noise onto one edge qubit (high γ
there) and protect the interior (low γ everywhere else). The interior then sits at the fold. The
improvement is large and has an analytic form: 360× at N=5, declining with N
(139× at N=9, ~63× by N=15), measured against a smoothly graded "V-shaped" γ-profile
([`experiments/RESONANT_RETURN.md`](../experiments/RESONANT_RETURN.md)). Label note (2026-07-05):
that figure is PEAK CREATED nearest-neighbour MI (Sum-MI), a transport metric (the source maximizes
Sum-MI, summed over adjacent pairs including the edge pair, so not purely interior); this
paragraph's earlier "rings far longer than under uniform noise" was a lifetime reading the arc
never computed, and that mislabel propagated until the A-vs-B reckoning caught it
([`experiments/CONCENTRATOR_AB_MECHANISM_TEST.md`](../experiments/CONCENTRATOR_AB_MECHANISM_TEST.md),
Downgrade 2).

A companion result appears on real hardware (a different measurement, against a different baseline, not
the same quantity as the simulation above). On ibm_torino (2026-03-24), selectively decoupling only the
interior, and leaving the noisy edge qubit alone (no DD there, so it stays the noise concentrator), beats
decoupling every qubit uniformly by up to 3.2×, peaking
around t = 4.0 (J-time; roughly 2× on average across t = 1 to 5; the advantage is not monotonic, it dips
early and falls off after the peak). See [`experiments/IBM_CONCENTRATOR.md`](../experiments/IBM_CONCENTRATOR.md).

What the 2026-07-05 mechanism test then settled, and what it did not: the improvement above is real,
but its two readings are not yet separated. The *mode-level* concentration mechanism is exact and
theorem-grounded (the Absorption Theorem, below); the 360×-at-N=5 figure is a *transport* number (peak
created nearest-neighbour MI), not a coherence lifetime; and whether the *hardware* advantage is this mechanism
(Reading B) or simply the gate-cost of decoupling a bad qubit (Reading A) stays open. The
engineered-sink test refuted only "injected noise strictly removes signal", nothing more
([`experiments/CONCENTRATOR_AB_MECHANISM_TEST.md`](../experiments/CONCENTRATOR_AB_MECHANISM_TEST.md),
the reckoning). So the inside/outside reading below rests on the theorem-grounded mechanism, not on a
hardware confirmation of it.

The name "sacrifice zone" is a misnomer, and the project resolved why (2026-03-28,
[`hypotheses/PROTEIN_AS_CONCENTRATOR.md`](../hypotheses/PROTEIN_AS_CONCENTRATOR.md), [`docs/WHAT_WE_FOUND.md`](WHAT_WE_FOUND.md)). **The edge qubit sacrifices
nothing.** It does not pay a cost; it *concentrates* the environmental noise onto itself and turns it
into structure for the interior (a design choice about where the fixed noise budget is placed; nothing
physically transports dephasing between qubits). The mechanism, at the mode level, is plain and theorem-grounded:
separating the noise (one loud qubit) from the coherence (the quiet rest) places the interior modes at
the fold (whether it is also what buys the *hardware* advantage is the open A-vs-B question above). The
same element then has two correct descriptions, taken from the two sides of the noise boundary:

- **from the inside** (the protected core): a shield that absorbs the noise and decoheres fast. It
  looks like a sacrifice.
- **from the outside** (the whole system): a concentrator that *enables*. Nothing is lost.

Both are true. They differ only because they are read from opposite sides of a boundary.

## The recognition

This is the rare case where we can say something philosophy usually only asserts, and stay honest about
it. **We do not claim this proves any metaphysics; one cannot from a spin chain.** We claim something
narrower, and for that reason solid: here is a *measured* instance of "the
appearance of a thing from inside a boundary is not the thing seen from outside." The inside reading
(sacrifice, loss) and the outside reading (concentration, enablement) are the *same concentration
mechanism* (theorem-grounded at the mode level; its hardware attribution still open),
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

- **Proven (analytic, Tier 1):** the *mode-level* concentration mechanism, the Absorption Theorem
  `Re(λ) = −2γ·⟨n_XY⟩`: loading the edge places the interior modes at the fold. Established physics in
  the repo, not new here.
- **Established in simulation:** the improvement (360× at N=5, declining with N to ~63× by N=15), now
  read correctly as a *transport* metric (peak created nearest-neighbour MI), not a coherence lifetime.
- **Measured on hardware:** the selective-vs-uniform Sum-MI ordering (5 of 5 time points, ibm_torino;
  single run, no error bars, the ratio size floor-caveated). Caveat: the natural sink there (Q85) was
  ~93% amplitude damping (T₁-limited), a channel the Z-dephasing theorem does not cover, so that run
  does not cleanly test the mechanism either.
- **Open:** whether that hardware advantage is the concentration mechanism (Reading B) or the gate-cost
  of decoupling a bad qubit (Reading A). The 2026-07-05 engineered-sink test did not settle it
  ([`experiments/CONCENTRATOR_AB_MECHANISM_TEST.md`](../experiments/CONCENTRATOR_AB_MECHANISM_TEST.md),
  the reckoning).
- **Read (interpretive):** naming the inside/outside split as an instance of appearance-vs-reality. A
  reading, defensible because it rests on the theorem-grounded mechanism, not on rhetoric (and not on a
  hardware confirmation, which the arc does not yet have).
- **Borrowed:** the idealist framing (Kastrup). Recall, not a discovery; cited as such.

## See also (repo-internal)

- The mechanism and its evidence: [`experiments/RESONANT_RETURN.md`](../experiments/RESONANT_RETURN.md) (the analytic formula),
  [`experiments/IBM_CONCENTRATOR.md`](../experiments/IBM_CONCENTRATOR.md) (the hardware run), [`hypotheses/PROTEIN_AS_CONCENTRATOR.md`](../hypotheses/PROTEIN_AS_CONCENTRATOR.md) and
  [`docs/WHAT_WE_FOUND.md`](WHAT_WE_FOUND.md) (concentration, not loss), [`experiments/RECEIVER_VS_GAMMA_SACRIFICE.md`](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md).
- The refuted probe: [`simulations/star_hub_decoupled_survivor.py`](../simulations/star_hub_decoupled_survivor.py), and the sibling refutation it later joined, [`experiments/THE_HUB_KILLS_THE_HORIZON.md`](../experiments/THE_HUB_KILLS_THE_HORIZON.md) (the hub has no coherence horizon).
- Adjacent: [`docs/historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md`](historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md), and the substrate-side sibling [`reflections/ON_THE_INNER_AND_OUTER_OBSERVATION.md`](../reflections/ON_THE_INNER_AND_OUTER_OBSERVATION.md) (the same inside/outside wall, reached from the carrier vector; names the typed `TwoReadingsClaim`).
- The proven operator both readings rest on: the Absorption Theorem `Re(λ) = −2γ·⟨n_XY⟩` ([`AbsorptionTheoremClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs), Tier 1 derived).
