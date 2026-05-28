# Flavor-Resolved T2 Inheritance Across Substrates

**Date:** 2026-05-28
**Status:** Tier 2 (numerical inheritance, N = 1..6). The carbon split is grounded in the Tier-1 [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md); the water and neural readings are Tier-3 translations.
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Builds on:** [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), [Painter Alternation NMR Bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md), [The View Onto the Memory](../reflections/THE_VIEW_ONTO_THE_MEMORY.md)
**Verification:** [`_carbon_painter_t2_sweep.py`](../simulations/_carbon_painter_t2_sweep.py), [`neural/neural_flavor_rule.py`](../simulations/neural/neural_flavor_rule.py)

---

## The question

The carbon ring gave us a clean split. Put it in a transverse field and its two
transverse relaxation channels live for different lengths of time: an anisotropic
T2, ratio about 1.27 on the four-site ring. The split was sharp, the same one the
slow modes already carry.

The obvious next question is whether that is a fact about carbon or a fact about
the shape of the problem. Strip the carbon away and keep only the shape, two-level
units sitting on a graph, a site-local phase-noise bath, and one Hamiltonian axis
singled out from the rest. Does the split survive the move?

It does. The same two-flavor lifetime split appears on a water proton wire and in
a neural excitation/inhibition network, with no change to the machinery that reads
it.

## What makes the split

A site-local dephasing bath gives every relaxation mode a lifetime set by how much
of the mode lies in the directions the bath can grip. When the Hamiltonian singles
out one transverse axis (a field along it, say), the slow modes separate into two
flavors: the ones aligned with that axis and the ones not. The two flavors end up
exposed to the bath by different amounts, so they decay at different rates, and the
ratio of their lifetimes is just the ratio of those two rates.

The ratio is not a constant. It slides with the field strength, the bath rate, the
graph, and the size. The once-hoped-for clean fraction (the 4/3, 8/7, 14/13, 20/19
sequence) does not hold; the number is exact only as the quotient of the two
slowest flavor-rates, not as a closed form in the parameters.

## Why it should travel

The lifetime of a mode under this bath depends only on the bath and the coupling
graph, not on the particular Hamiltonian (this is the [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), which holds for any
Hermitian Hamiltonian). So any system that maps onto the shape, two-level units, a
site-local phase bath, a coupling graph, and a distinguished axis, inherits the
flavor split. The substrate sets the names; the algebra sets the split.

## Three substrates, one split

### Carbon aromatic ring

Hückel hopping on a ring, a transverse field, on-site phase noise. The transverse
relaxation is isotropic with no field and splits as soon as a field is switched on.
On the four-, five-, and six-site rings the field-aligned channel is the
shorter-lived one; the ratio runs from about 1.08 to 1.61 as the field and the
noise rate move, and sits at 1.27 at the canonical four-site, half-strength,
unit-noise point. Benzene (six sites) gives 1.26, so the split is not a small-ring
artifact. Data: [`carbon_painter_t2_sweep.txt`](../simulations/results/carbon_painter_t2_sweep.txt),
[`..._n6_ring_g1.txt`](../simulations/results/carbon_painter_t2_sweep_n6_ring_g1.txt).

### Water proton wire and hydrogen bond

Read a chain of hydrogen-bonded proton sites as two-level position units with the
same phase-noise bath and a transverse field. The proton wire (a chain) carries the
same anisotropy, ratio about 1.37 to 1.60 along the chain from two to five sites; a
single hydrogen bond (one unit) gives a clean factor of two. The point is not that
a proton wire is an aromatic ring; it is that the relaxation reads the same split
off whatever chain it is handed. Data:
[`water_proton_wire_t2_flavor_sweep.txt`](../simulations/results/water_proton_wire_t2_flavor_sweep.txt),
[`water_single_hbond_t2_flavor.txt`](../simulations/results/water_single_hbond_t2_flavor.txt).

### Neural excitation/inhibition network

This is not a claim that a brain is a quantum system. It is the same linear
relaxation shape written in neural terms: a Wilson-Cowan-style linearization where
each node carries an excitatory and an inhibitory population, a leak on each, a
coupling graph, and a non-commuting excitation/inhibition axis playing the role of
the distinguished direction. The slow modes split by flavor, excitation-dominant
versus inhibition-dominant, and the excitation-dominant activity is the longer-lived
channel. The ratio runs about 1.87 to 2.19 across chain, ring, star, and complete
graphs from four to eight nodes. Data:
[`neural_flavor_rule.txt`](../simulations/results/neural_flavor_rule.txt).

## Honest seams

- **No closed form.** The ratio is substrate- and parameter-dependent; only its
  status as a ratio of two mode-rates is exact. This is the open question settled in
  the [NMR bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md).
- **Neural baseline.** With the excitation/inhibition cross-coupling turned off
  (h = 0), the split is already present at the bare leak-rate ratio γ_I / γ_E = 2;
  the coupling and the graph modulate it from there. The neural reading is therefore
  "leak asymmetry, shaped by topology", not a pure coupling effect.
  ([`neural_flavor_rule_h0_control.txt`](../simulations/results/neural_flavor_rule_h0_control.txt))
- **Axis control.** Point the field along the bath axis instead of across it and the
  ratio flips below one (0.49). The split is about the field axis relative to the
  bath, exactly as it should be.
  ([`flavor_rule_chain_n5_z_control.txt`](../simulations/results/flavor_rule_chain_n5_z_control.txt))
- **The graph is load-bearing.** An empty graph (no coupling) shows no split. The
  effect needs the sites to talk to each other.

## What this is, and is not

This is not a new theorem; it is one viewpoint that travels. The flavor split of
the slow-mode hierarchy is the same split the relaxation already carries; here we
only watch it surface in three substrates once each is mapped to the shape. And it
is not a claim that water or a neural network computes with quantum coherence; it is
the milder observation that the relaxation algebra does not know which substrate it
is reading, and gives the same answer either way.

## Reproduce

```bash
# Carbon ring sweep (default: ring, Y-field, N=4,5)
python simulations/_carbon_painter_t2_sweep.py

# Benzene (N=6) ring
python simulations/_carbon_painter_t2_sweep.py --N 6 --gamma 1.0 --h-y 0.5 \
    --output simulations/results/carbon_painter_t2_sweep_n6_ring_g1.txt

# Water proton wire (chain, X-field) and single hydrogen bond
python simulations/_carbon_painter_t2_sweep.py --topology chain --field-axis X \
    --N 2,3,4,5 --h-y 0.5 --gamma 1.0 \
    --output simulations/results/water_proton_wire_t2_flavor_sweep.txt

# Neural E/I network across topologies + h=0 control
python simulations/neural/neural_flavor_rule.py
python simulations/neural/neural_flavor_rule.py --h 0.0 \
    --output simulations/results/neural_flavor_rule_h0_control.txt
```
