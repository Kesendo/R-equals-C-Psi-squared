# Water-domain map: neutral proton-position wire

This directory maps selected framework statements onto a deliberately narrow
theoretical model: a **chosen two-level proton-position coordinate** (`|L⟩` /
`|R⟩`) in a neutral, oriented, single-file, one-dimensional proton wire. It is
not a model of a whole water molecule, bulk liquid water, ice, a branched
network, an excess proton, a charged Grotthuss carrier, or a transport
calculation.

The local Z-dephasing channel and the selected Hamiltonian (including its
coupling `J`) are modelling assumptions. They are not measured facts about
water. The canonical scope, dipole map, and exclusions are in [the proton-wire
crossing pass](PROTON_WIRE_CROSSING.md).

The crossing pass controls the scope of this directory. Linked sibling pages
record finite calculations for the selected coordinate; they do not enlarge
this neutral-wire model into a molecular, charged-carrier, or transport claim.

## Map of this directory

- [HYDROGEN_BOND_QUBIT.md](HYDROGEN_BOND_QUBIT.md) introduces the chosen
  two-position coordinate for one O-H···O linkage.
- [PROTON_WATER_CHAIN.md](PROTON_WATER_CHAIN.md) records the finite-chain
  model calculations.
- [PROTON_WIRE_CROSSING.md](PROTON_WIRE_CROSSING.md) supplies the domain map:
  the neutral-wire dipole identity, the fixed-dipole boundary, the bias
  exception, and the excluded physical situations.

The associated model scripts live in
[`simulations/water/`](../../simulations/water/), including
`proton_wire_crossing.py`, `proton_chain_dicke_anchor.py`, and
`proton_water_chain.py`.

## Embedding conditions, kept separate

The following are different kinds of statement; none follows merely from
calling the coordinate a proton qubit.

1. **Exact model algebra.** [F1](../ANALYTICAL_FORMULAS.md) is exact only
   under its stated Hamiltonian and local-Z-dephasing hypotheses. In this
   water mapping, the un-biased model must exclude the displayed longitudinal
   double-well bias `Δ Σ_l Z_l`; see [the bias exception](PROTON_WIRE_CROSSING.md#the-bias-field-is-the-same-operator).
   The un-biased transverse-field-Ising calculation has a machine-zero F1
   residual in this selected model. That is a model result, not a property
   measured for water.

2. **Fixed-dipole structural results.** These require their own premises,
   rather than one blanket “fixed-dipole” qualifier.

   - [F4](../ANALYTICAL_FORMULAS.md)'s principal stationary-mode formula is
     for the Heisenberg Hamiltonian at `Σγ = 0`. Its separate
     **dephased-kernel extension** concerns a Heisenberg connected component
     with local Z-dephasing on every site; it is this extension that yields
     the `N + 1` kernel count under its stated conditions.
   - F86b's `3/8` K-intermediate anchor, stated in [F88b's
     multi-state Dicke extension](../ANALYTICAL_FORMULAS.md), is a static
     formula for its specified even-`N` Dicke superposition.
   - [F88b](../ANALYTICAL_FORMULAS.md) is the stated Krawtchouk formula
     for its popcount-coherence inputs; its dynamical memory reading retains
     F88b's own conditions.
   - [F98](../ANALYTICAL_FORMULAS.md) additionally requires the
     K-intermediate Dicke state, a magnetization-conserving Hamiltonian
     (`[H, Ŵ] = 0`), all-site Z-dephasing, and a sector-uniform initial
     diagonal. It is not implied by the F1 palindrome class.

3. **Physical parameter assumptions.** The position coordinate, the choice
   of `J`, and the selected Hamiltonian and decoherence channel are inputs to
   the model. The repository does not identify a measured water value of
   `Q = J/γ` from these calculations.

4. **Excluded situations.** Bulk water, ice, branches, an excess proton, a
   charged Grotthuss carrier, and a transport rate are outside this state
   space. The [crossing pass's scope](PROTON_WIRE_CROSSING.md#scope-stated-plainly)
   explains why.

## Comparison to hardware

The pure-Z model is unital and fixes every Z-diagonal density matrix; it does
not alone select `I/d`, a Gibbs β, or a unique steady state. In the
crossing-pass calculations at `N = 3, 4, 5`, the selected un-biased TFI model
has a one-dimensional kernel and reaches `I/d`; that finite-model calculation
is not a general TFI theorem. The all-site-dephased Heisenberg branch instead
retains the F4 sector structure. Hardware T1 amplitude damping introduces a
preferred direction and therefore changes the comparison. See
[`simulations/memory_reading_ibm_torino.py`](../../simulations/memory_reading_ibm_torino.py).

Repository-recorded hardware flights test framework qubit models; they do not
confirm water chemistry. This is not a statement about what hardware evidence
may exist outside the repository.

## F86b and F98 in the selected model

For even `N`, the K-intermediate Dicke state

```
ψ = (|D_{N/2−1}⟩ + |D_{N/2}⟩) / √2
```

has the F86b static output

```
α(t = 0) = 3/8.
```

With the selected Heisenberg Hamiltonian and all-site Z-dephasing,
[`proton_chain_dicke_anchor.py`](../../simulations/water/proton_chain_dicke_anchor.py)
reproduces that structural identity and the F98 model output

```
‖P_{N/2−1, odd}‖² = C(N, N/2−1) / 2
α(∞)_KIntermediate = (N + 2) / [4(N + 1)] → 1/4.
```

These are formula and model outputs, not water measurements. More generally,
F98 is the result with the K-intermediate Dicke state, a
magnetization-conserving `H`, all-site Z-dephasing, and a sector-uniform
initial diagonal; its conclusion does not cover arbitrary initial states or a
Hamiltonian that changes `Ŵ`. The [F98 scope gate](../../simulations/f98_scope.py)
and [crossing pass](PROTON_WIRE_CROSSING.md) give the contrasting cases.

## Dipole crossing and the bias boundary

Only for a neutral point-charge wire with equal O···O spacing and a chosen
orientation in one dimension does

```
μ = Ŵ = Σ_l (I − Z_l)/2
```

hold, in units of elementary charge times the spacing. Under those conditions
`[Ŵ, H] = 0` fixes the **total** dipole; it does not prohibit movement.
`XX+YY` can move protons in correlated opposite directions while their net
dipole displacement remains zero. The full derivation and its limits are in
[the crossing pass](PROTON_WIRE_CROSSING.md#the-answer-popcount-is-the-wires-dipole-moment).

The un-biased TFI F1 result is therefore only a result of this selected model.
The longitudinal `Δ Σ_l Z_l` bias is the dipole coupled to a uniform field and
is the relevant F1 exception here; see [the bias-field discussion](PROTON_WIRE_CROSSING.md#the-bias-field-is-the-same-operator).

## Open follow-ups

- **Temperature and parameters.** A chemistry-grounded relation between a
  specified environment and the model inputs `J` and `γ` remains to be
  established. No liquid-water `Q` is asserted here.
- **Spectroscopy and pump-probe.** Whether a specified observable of a
  confined neutral wire could read this model's dipole coordinate is an open
  experimental and literature question. No observed water signature is
  claimed.
- **DNA tautomer coordinates.** A DNA application would need its own
  coordinate, Hamiltonian, channel, and scope check; it does not inherit this
  water-wire mapping automatically. The existing context is
  [`experiments/DNA_BASE_PAIRING.md`](../../experiments/DNA_BASE_PAIRING.md).
- **Charged and branched extensions.** Extending the coordinate to an excess
  proton, a branch, or a transport observable is separate work, not a
  continuation inside this neutral one-dimensional state space.
