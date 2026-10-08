# Flavor-Resolved T2 Inheritance Across Substrates

**Date:** 2026-05-28
**Status:** Tier 2 (numerical inheritance, N = 1..6) for carbon and water: the carbon split is grounded in the Tier-1 [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), the water reading is a Tier-3 translation. The neural network is not an inheritance: its split is [F36](../docs/ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra)'s pairing with a closed form (Tier 1 algebra, see its section).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Builds on:** [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), [Painter Alternation NMR Bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md), [The View Onto the Memory](../reflections/THE_VIEW_ONTO_THE_MEMORY.md)
**Verification:** [`carbon_painter_t2_sweep.py`](../simulations/carbon_painter_t2_sweep.py), [`neural/neural_flavor_rule.py`](../simulations/neural/neural_flavor_rule.py), [`neural/flavor_split_is_f36.py`](../simulations/neural/flavor_split_is_f36.py)

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

On the water wire it does. In a neural excitation/inhibition network a two-flavor
lifetime split appears too, read by the same machinery, but it is not this split
carried over: that network has no bath, and its split is the neural palindrome's own
pairing, in closed form (its section below).

## What makes the split

A site-local dephasing bath gives every relaxation mode a lifetime set by how much
of the mode lies in the directions the bath can grip. When the Hamiltonian singles
out one transverse axis (a field along it, say), the slow modes separate into two
flavors: the ones aligned with that axis and the ones not. The two flavors end up
exposed to the bath by different amounts, so they decay at different rates, and the
ratio of their lifetimes is just the ratio of those two rates.

The ratio is not a constant. It slides with the field strength, the bath rate, the
graph, and the size. The once-hoped-for clean fraction (the 4/3, 8/7, 14/13, 20/19
sequence) does not hold; for carbon and water the number is exact only as the
quotient of the two slowest flavor-rates, not as a closed form in the parameters (the
neural look-alike has one, below).

## Why it should travel

The lifetime of a mode under this bath depends only on the bath and the coupling
graph, not on the particular Hamiltonian (this is the [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), which holds for any
Hermitian Hamiltonian). So any system that maps onto the shape, two-level units, a
site-local phase bath, a coupling graph, and a distinguished axis, inherits the
flavor split. The substrate sets the names; the algebra sets the split. The neural
network below does not map onto the shape, and gets its split from elsewhere.

## Three substrates, one shape of split

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

This is not a claim that a brain is a quantum system. It is a linear relaxation of
the same outward shape written in neural terms, with leaks in place of a bath: a Wilson-Cowan-style linearization where
each node carries an excitatory and an inhibitory population, a leak on each, a
coupling graph, and a non-commuting excitation/inhibition axis playing the role of
the distinguished direction. The slow modes split by flavor, excitation-dominant
versus inhibition-dominant, and the excitation-dominant activity is the longer-lived
channel. The ratio runs about 1.83 to 2.19 across chain, ring, star, and complete
graphs from four to eight nodes. Data:
[`neural_flavor_rule.txt`](../simulations/results/neural_flavor_rule.txt).

This split is not an inheritance of the Absorption split. The generator has no density
matrix and no dephasing dissipator, only leaks, and what it carries instead is the neural
palindrome, [F36](../docs/ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra),
exactly: with Q the per-node excitation/inhibition swap and `s = (γ_E + γ_I)/2`, the generator
`[[−γ_E I + cA, −hI], [hI, −γ_I I − cA]]` (c the graph coupling, h the excitation/inhibition
cross-coupling, A the adjacency scaled to spectral radius one) satisfies `Q J Q = −J − 2sI` at
every h, every coupling, every graph and every pair of leaks: entry for entry a sign flip plus
one sum, `γ_E + γ_I` against 2s, the same rounded sum on both sides, so the float residual is
exactly 0.0 at generic couplings; +cA on both blocks, or a wrong s, is the control that fails.
And it has a closed form. Per eigenvalue a of A, with eigenvector w, the pair (w, 0), (0, w)
spans an invariant block `−sI + [[Δ + ca, −h], [h, −Δ − ca]]`, `Δ = (γ_I − γ_E)/2`, with
eigenvalues `λ = −s ± √((Δ + ca)² − h²)`. A block is overdamped, two real rates mirrored
about s, exactly when `|Δ + ca| > h`, and its slow mode is excitation-dominant exactly when
`Δ + ca > 0`; an underdamped block sits on `Re λ = −s` with excitation weight exactly one
half, the classifier's "mixed"; a block with `Δ + ca = ±h` is defective, a Jordan block at
−s whose eigenvector has weight one half, and at this page's parameters that is every a = 0
block (the empty graph, the star, odd chains, rings with 4 | N). So the flavor split is the
mirror pairing read on the overdamped blocks: every excitation-dominant mode at rate r has
an inhibition-dominant partner at 2s − r. The ratio this page reports pairs the slowest
excitation-dominant rate, `s − √((Δ + c)² − h²)` from the block a = 1, the same 0.940983 on
every connected graph, with the slowest inhibition-dominant rate, which is the faster rate
`s + √(·)` of the weakest overdamped block, so the graph enters through which blocks are
overdamped: chain N = 4 reads `(1.5 + 0.3235)/(1.5 − 0.5590) = 1.9378`; the star and the
complete graph at every N and the ring at N = 4 have a single overdamped block and read
2.1882; every row of the data file and of the h = 0 control reproduces from the closed form,
with the classifier's own condition, the weakest block's fast mode carrying inhibition
weight above 0.65, read beside each. At h = 0 every block with `Δ + ca ≠ 0` is overdamped
and the ratio is `(γ_I + c·a_min)/(γ_E − c·a_max)`: 2.333 on chain, ring and star, 2.6 on the
complete graph, 2 on the empty graph. The equal-leak control
([`neural_flavor_rule_equal_gamma_control.txt`](../simulations/results/neural_flavor_rule_equal_gamma_control.txt))
reads "classifier-failed" on every topology because at Δ = 0 with h > c every block is
underdamped and every weight is one half; with c > h overdamped blocks appear at equal leaks
too, but there the blocks a and −a share their eigenvalues, so on a bipartite graph each
overdamped rate is carried by one mode of each flavor at once and only a non-bipartite graph
splits genuinely. Which flavor is slow is the sign of Δ + ca and not F36: at γ_E = 1,
γ_I = 2, h = 0.1, c = 1 on the chain, inhibition-dominant modes sit below s. Gate
[`neural/flavor_split_is_f36.py`](../simulations/neural/flavor_split_is_f36.py), 65 checks,
on the producer's own generator and its data files. [The view onto the memory](../reflections/THE_VIEW_ONTO_THE_MEMORY.md)
says the neural translation carries the palindrome rather than the split; this is that
sentence with the data files' own rows under it.

## Honest seams

- **No closed form for carbon and water.** There the ratio is substrate- and
  parameter-dependent; only its status as a ratio of two mode-rates is exact. This is
  the open question settled in the [NMR bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md).
  The neural ratio has one, above.
- **Neural baseline.** With the excitation/inhibition cross-coupling turned off
  (h = 0), the split is already present: `(γ_I + c·a_min)/(γ_E − c·a_max)`, which is the
  bare leak ratio γ_I / γ_E = 2 on the empty graph and 2.333 on chain, ring and star, 2.6
  on the complete graph at c = 0.25. The neural reading is therefore "leak asymmetry,
  shaped by topology", not a pure coupling effect; by the closed form above it is the
  mirror pairing about the mean leak, the graph choosing which blocks are overdamped.
  ([`neural_flavor_rule_h0_control.txt`](../simulations/results/neural_flavor_rule_h0_control.txt))
- **Axis control.** Point the field along the bath axis instead of across it and the
  ratio flips below one (0.49). The split is about the field axis relative to the
  bath, exactly as it should be.
  ([`flavor_rule_chain_n5_z_control.txt`](../simulations/results/flavor_rule_chain_n5_z_control.txt))
- **The graph is load-bearing for carbon and water.** There an empty graph (no
  coupling) shows no split; the effect needs the sites to talk to each other. The neural
  empty graph at the default parameters sits exactly at an exceptional point (Δ = h), a
  defective block at −s whose eigenvector has weight one half, and at h = 0 it splits at
  γ_I/γ_E; what the neural graph does is decide which blocks are overdamped.

## What this is, and is not

This is not a new theorem; it is one viewpoint that travels. The flavor split of
the slow-mode hierarchy is the same split the relaxation already carries; here we
watch it surface in two substrates once each is mapped to the shape, and meet its
look-alike in a third, the neural network, where the same reading machinery finds a split
that the neural palindrome makes rather than the bath. For the two substrates the
relaxation algebra does not know which one it is reading; the third is not read by it. And it
is not a claim that water or a neural network computes with quantum coherence; it is
the milder observation that the relaxation algebra does not know which substrate it
is reading, and gives the same answer either way.

## Reproduce

```bash
# Carbon ring sweep (default: ring, Y-field, N=4,5)
python simulations/carbon_painter_t2_sweep.py

# Neural: the split is F36's pairing, closed form, on the producer's own generator
python simulations/neural/flavor_split_is_f36.py
python simulations/neural/neural_flavor_rule.py --gamma-e 1.0 --gamma-i 1.0 --N 6 \
    --output simulations/results/neural_flavor_rule_equal_gamma_control.txt   # the equal-leak control

# Benzene (N=6) ring
python simulations/carbon_painter_t2_sweep.py --N 6 --gamma 1.0 --h-y 0.5 \
    --output simulations/results/carbon_painter_t2_sweep_n6_ring_g1.txt

# Water proton wire (chain, X-field) and single hydrogen bond
python simulations/carbon_painter_t2_sweep.py --topology chain --field-axis X \
    --N 2,3,4,5 --h-y 0.5 --gamma 1.0 \
    --output simulations/results/water_proton_wire_t2_flavor_sweep.txt

# Neural E/I network across topologies + h=0 control
python simulations/neural/neural_flavor_rule.py
python simulations/neural/neural_flavor_rule.py --h 0.0 --N 6 \
    --output simulations/results/neural_flavor_rule_h0_control.txt
```
