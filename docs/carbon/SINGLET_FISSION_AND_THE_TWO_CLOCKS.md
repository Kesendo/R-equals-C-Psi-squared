# Singlet Fission and the Two Clocks: A Conditional Bridge

**Date:** 2026-05-30
**Authors:** Tom + Claude
**Status:** Tier 3 conditional analogy. The V-Effect exchange `J_eff = (3/8)·α²/J` is a Tier 1-2 framework/model result; the proposed connection to singlet fission is neither a physical mapping nor a confirmed mechanism.
**Builds on:** [exchange from the V-Effect](../../experiments/EXCHANGE_FROM_V_EFFECT.md), [the Frost circle as the clock face](FROST_CIRCLE_AS_THE_CLOCK_FACE.md), the F87 trichotomy, and [Heisenberg reloaded](../../hypotheses/HEISENBERG_RELOADED.md).

---

## Two distinct framework readings

The repository has two separately defined model readings. The Frost-circle work reads a selected XX+YY plus dephasing model through radial/decay and angular/frequency coordinates. The V-Effect work studies a selected direct-Heisenberg system of two singlet pairs joined by a bridge `α`; F87 is a separate framework classification.

They are not thereby two physical descriptions of one carbon degree of freedom. Treating one as a charge/orbital model and the other as a spin model is a conditional translation that needs a specified molecular Hamiltonian, bath, state preparation, and producer before it can be applied to a molecule.

## V-Effect result retained

In the named direct-Heisenberg calculation, each pair has intra-pair coupling `J` and the pairs are joined by a bridge `α`. The bridge couples the singlet-singlet ground state to the subspace in which both pairs are triplet. The latter is separated by

    8J = 2 × 4J,

and the small-`α` second-order ground-state shift is

    δE_GS^(2) = −(3/8) · α² / J,
    J_eff = (3/8) · α² / J.

The factor `3/8` is the named Pauli-algebra/direct-Heisenberg-bridge result, with its N=4 numerical check, in [exchange from the V-Effect](../../experiments/EXCHANGE_FROM_V_EFFECT.md). It is neither a molecular Hamiltonian nor a parameter assignment for a material. The producing calculation explicitly records that mapping `α` and `J` to a specific physical system requires independently identifying those parameters.

## The singlet-fission analogy

The structural analogy is limited: a selected model contains a singlet-singlet state, a both-pairs-triplet subspace, a bridge, and a second-order energy scale. Those named ingredients resemble labels used in singlet-fission discussions. The repository has not shown that the V-Effect both-pairs-triplet state is a carotenoid `2Ag` state, that either has the same microscopic content, or that the V-Effect bridge is a physical mixing operator for singlet fission.

Consequently, the relation between a molecular dark state, a charge-transfer sector, and a triplet-pair sector remains unassigned. The V-Effect calculation cannot select a molecular preparation, specify a molecular bath, establish an optical transition, or supply a physical `J`, `α`, `γ`, `T₂`, or `Q`.

## What the two clocks do not establish

The XX+YY/dephasing clock and the direct-Heisenberg V-Effect calculation are separate selected models. Their co-presence does not establish a common molecule, a common state space, a direct physical channel between their states, or a mechanism that mixes charge and spin in a material.

The framework depth/absorption and F87 classification remain framework statements on their respective operators. They do not by themselves identify a molecular bright state, dark state, optical selection rule, triplet-pair state, or singlet-fission pathway.

## The depth rail: a separate framework result

The [Absorption Theorem](../../experiments/ABSORPTION_THEOREM_DISCOVERY.md) sorts framework coherences by drain depth, the number of sites at which a bra and ket differ. Its parity relation is `n_diff ≡ Δpopcount (mod 2)`: odd `Δpopcount` has odd `n_diff`, while number-conserving (`Δpopcount = 0`) coherences have even `n_diff`.

That classification is not an optical assignment for a molecule. It does not identify a molecular bright or dark state, make a triplet pair a particular framework coherence, or show that the V-Effect bridge connects the two clock models.

## Open model questions

1. Specify a molecular Hamiltonian and state/preparation map, then determine whether its relevant subspaces admit either selected framework reading.
2. Specify a bath channel and rate model before comparing a molecular time scale with the selected dephasing clock or defining a material `Q`.
3. Construct a producer for the proposed bridge in that named model and test whether it connects the proposed charge-like and spin-like subspaces. A direct-Heisenberg V-Effect result alone does not do this.
4. Only after those choices, compare the model's calculated states and observables with a defined molecular target.

## Anchor

- **V-Effect producer/model:** [exchange from the V-Effect](../../experiments/EXCHANGE_FROM_V_EFFECT.md)
  and its [`level1_emergent_exchange.py` calculation](../../simulations/level1_emergent_exchange.py):
  the named four-site direct-Heisenberg system with intra-pair `J` and bridge `α`.
- **Clock producer/model:** [the Frost circle as the clock face](FROST_CIRCLE_AS_THE_CLOCK_FACE.md)
  and [`frost_circle_as_clock.py`](../../simulations/carbon/frost_circle_as_clock.py): the
  selected XX+YY/local-`D[n]` clock model. F87 remains a framework classification in the
  [formula registry](../ANALYTICAL_FORMULAS.md).
- **Scope:** the [Carbon source contract](README.md#carbon-source-contract) and the
  [Q audit](../Q_BELONGS_TO_NO_SUBSTANCE.md) keep these model results distinct from a
  material mapping.
- **Molecular mapping:** remains open for a named target. It requires a selected molecular
  degree of freedom, Hamiltonian/coupling convention, state preparation, bridge/channel,
  bath channel, and rate before the Tier 3 analogy can be tested against that target.
