# Proton Water Chain: Two Selected Models

<!-- Keywords: neutral 1-D proton-coordinate chain, proton qubit water chain, hydrogen-bond
network, palindromic Lindblad model, water-chain model, R=CPsi2 proton chain -->

**Status:** Tier 2 selected-model calculations
**Date:** April 1, 2026
**Script:** [proton_water_chain.py](../../simulations/water/proton_water_chain.py)
**Data:** [proton_water_chain.txt](../../simulations/results/proton_water_chain.txt)
**Related model scope:** [Hydrogen Bond Qubit](HYDROGEN_BOND_QUBIT.md),
[Analytical Formulas](../ANALYTICAL_FORMULAS.md), and
[Mirror Symmetry Proof](../proofs/MIRROR_SYMMETRY_PROOF.md)

---

## Scope

This page records two different N = 1..5 spin-model calculations. They can
motivate questions about a chosen proton-transfer coordinate, but neither
calculation validates a physical proton wire, Grotthuss transport, or liquid
water.

- **Heisenberg + all-site local Z dephasing** is the formula/model-validation
  branch. Its Hamiltonian is the XX+YY+ZZ chain and its dissipator has one
  local Z-dephasing operator at every site. The formula comparisons below are
  statements about this generator.
- **Unbiased transverse-field Ising (TFI) + all-site local Z dephasing** is a
  dipole-moving toy mapping: `-J Σ X_i + K Σ Z_i Z_{i+1}` represents chosen
  tunnelling and ZZ-coupling terms. It is not the Heisenberg calculation and
  does not inherit its popcount, fixed-dipole, or stationary-kernel results
  merely because it is called a proton-wire model.

For the second, un-biased TFI generator, the F1 palindrome can hold under the
named local-Z model conditions. That is a property of the specified
Liouvillian, not a consequence of the words “water” or “proton wire.” A
longitudinal bias, a different channel, or a different generator must be
checked against F1's own premises before making a pairing claim.

## Q provenance and model labels

The displayed TFI sweeps are selected illustrative model runs, not a
measurement or regime verdict for ordinary liquid water. [Q Belongs to No
Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md) gives only the illustrative
conditional ceiling `Q ≲ 4.6`: it requires a selected proton coordinate, the
ice-derived `J = 0.5 meV` convention, and a 1–3 ps H-bond-lifetime proxy for
the missing proton-coordinate `T₂`. It has no lower endpoint, is unvalidated
for liquid water, and does not calibrate the model's `J`, `γ`, or `K` for
liquid water.

Two different quantities are called `Q` here, and they must not be
interchanged.

- **`Q_max` in the spectrum tables** is the spectral ratio
  `max |Im λ| / (−Re λ)` over oscillatory Liouvillian eigenvalues. It is an
  output of the displayed Heisenberg or TFI generator.
- **`Q = J / γ₀` in the provenance range and F86 scan** is an input control
  ratio: `J` is the XY hopping coefficient in
  `H_XY = Σ_b (J/2)(X_b X_(b+1) + Y_b Y_(b+1))`, and `γ₀` is the uniform
  per-site local-Z-dephasing rate. `Q_peak` and its HWHM are response features
  along this control axis, not `Q_max` values.

The TFI Hamiltonian uses its separately defined tunnelling coefficient in
`−J Σ_i X_i + K Σ_b Z_b Z_(b+1)`; its `Q_max` cells remain spectral outputs.
Likewise, the all-Pauli Heisenberg normalization used in the full-state F88b
diagnostic below is stated there explicitly. A numerical equality of any of
these ratios would not identify their Hamiltonians or their observables.

The script's labels such as “Water,” “Enhanced,” and “Zundel” name parameter
rows. They do not classify water, ice, enzymes, or individual
proton-transfer events.

---

## Heisenberg branch: formula/model-validation calculation

The following output is from the open, all-site-Z-dephased Heisenberg chain at
`J = γ = 1`. The `V(N)` and `Q_max` columns compare the selected generator to
the corresponding formula outputs; they are not measurements of a molecular
system.

| N | Eigenvalues | Frequencies | Q_max | V(N) (run) | V(N) (formula) |
|---|-------------|-------------|-------|------------|----------------|
| 1 | 4 | 0 | 0.000 | 0.000 | 0.000 |
| 2 | 16 | 2 | 2.000 | 1.000 | 1.000 |
| 3 | 64 | 5 | 3.000 | 1.500 | 1.500 |
| 4 | 256 | 34 | 3.414 | 1.707 | 1.707 |
| 5 | 1024 | 109 | 3.618 | 1.809 | 1.809 |

The same output records the following whole-spectrum minimum-rate diagnostic
at this selected `J = γ = 1` point. It must not be read as a violation of a
pure-sector bound: mixed Hamiltonian sectors can have rates below the
weight-one value. The relevant mode-count discussion is
[D05 Dynamic Mode Count](../proofs/derivations/D05_DYNAMIC_MODE_COUNT.md).

| N | Pure weight-one reference | Whole-spectrum minimum in the run |
|---|---------------------------|-----------------------------------|
| 2 | 2.000 | 2.000 |
| 3 | 2.000 | 2.000 |
| 4 | 2.000 | 0.978 |
| 5 | 2.000 | 0.617 |

### Formula premises that stay on this branch

F4 has two distinct statements and premises.

- Its principal Clebsch–Gordan stationary-mode formula is for the Heisenberg
  Hamiltonian with `Σγ = 0`; it is not the dephased run above.
- Its dephased-kernel extension requires Z dephasing on **every** site of the
  Heisenberg graph. For graph components `c`, it gives
  `dim ker L = Π_c (|c| + 1)`; the familiar `N + 1` count is the connected
  case. It is not a fixed-dipole assertion for the TFI model. See
  [the F4 component proof](../proofs/PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md).

F86b's selected Dicke/popcount response and F88b's popcount-coherence result
are mathematical statements, not water properties, but they belong to
different generators. The F86 Q-peak/HWHM response is a **Tier-1 candidate**
from the uniform XY `(n, n+1)` coherence block under all-site local Z
dephasing, scanned in the control `Q = J / γ₀`. F88b separately classifies the
Π²-odd memory fraction of specified popcount-coherence states for the
Heisenberg-plus-all-site-local-Z kernel. Their derivations are in [the F86
Q-peak proof](../proofs/PROOF_F86_QPEAK.md) and the [formula
registry](../ANALYTICAL_FORMULAS.md).

F98 has a still narrower premise: even `N`, the K-intermediate Dicke initial
state `( |D_(N/2-1)> + |D_(N/2)> ) / √2`, whose initial Z-basis diagonal is
uniform within its two popcount sectors, a magnetization-conserving Hamiltonian
`[H, W] = 0`, and Z dephasing on every site. Under those conditions its
long-time Π²-odd asymptote is `(N + 2) / [4(N + 1)]`; the long-time
sector-uniform diagonal is a consequence. It does not apply to arbitrary
initial states or to the dipole-moving TFI branch.

---

## TFI branch: selected dipole-moving toy-model outputs

For `J = 50`, `K = 20`, and `γ = 50` in the script's units, the two models
produce the following output. The comparison holds the displayed model inputs
fixed; it does not identify either parameter set with a chemical system.

| N | Model | Frequencies | Q_max | Minimum rate |
|---|-------|-------------|-------|--------------|
| 2 | Heisenberg | 2 | 2.000 | 100.0 |
| 2 | TFI | 3 | 1.792 | 59.0 |
| 3 | Heisenberg | 5 | 3.000 | 100.0 |
| 3 | TFI | 15 | 1.811 | 58.9 |
| 4 | Heisenberg | 36 | 3.414 | 48.9 |
| 4 | TFI | 47 | 1.820 | 58.9 |
| 5 | Heisenberg | 109 | 3.618 | 30.8 |
| 5 | TFI | 222 | 1.826 | 58.9 |

The high-Q parameter row in the same TFI output is retained below as a model
row, without treating its label as a property of a Zundel cation or of water.

| N | Frequencies | Q_max | Minimum rate |
|---|-------------|-------|--------------|
| 1 | 1 | 9.95 | 50.0 |
| 2 | 3 | 9.96 | 54.3 |
| 3 | 15 | 9.96 | 54.2 |
| 4 | 46 | 9.96 | 54.2 |
| 5 | 228 | 9.96 | 54.2 |

### Selected thermal and dephasing-profile runs

The `N = 3` thermal table is a selected TFI calculation at `T = 300 K` in the
script. Its “warm” column adds amplitude damping to local Z dephasing, so it
is outside the local-Z-only F1 model condition and is not a thermal prediction
for liquid water.

| Property | Local Z only | Local Z plus amplitude damping |
|----------|--------------|-------------------------------|
| Frequencies | 15 | 23 |
| Q_max | 1.81 | 0.43 |
| Rate range | 59–300 | 226–697 |

The `N = 5` non-uniform-dephasing sweep is also a TFI model calculation. Its
values select long-lived modes under the prescribed `γ` profiles; they do not
establish a water dephasing profile or an experimentally protected proton
mode.

| Prescribed γ profile | Q_max | Improvement over displayed uniform row |
|----------------------|-------|----------------------------------------|
| Uniform `[50]*5` | 1.83 | baseline |
| Edge `[100, 10, 10, 10, 10]` | 9.24 | 5.1× |
| Both edges `[100, 10, 10, 10, 100]` | 7.96 | 4.4× |
| Center `[10, 10, 100, 10, 10]` | 8.43 | 4.6× |

For the N = 5 Heisenberg chain with all-site Z-dephasing profiles, the
[cavity-mode calculation](../../experiments/CAVITY_MODE_LOCALIZATION.md)
reports a `0.994` correlation between prescribed sacrifice weight and decay
rate. Separately, at N = 5 with uniform `γ₀ = 0.05` on every site and uniform
`J = 1`, the selected alt-z-bits receiver gives an `11.5×` Peak Sum-MI advantage
over that source's γ-sacrifice reference. See
[the fixed-receiver comparison](../../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md).
Neither is a water observable.

### Parameterized-model comparison

The script also compares named parameterizations at identical selected inputs;
equal rows mean that the input Hamiltonian and dissipator are the same. This
is a parameterized-model comparison, not an equivalence between water and DNA.

| Parameterization label in output | N | Frequencies | Q_max |
|----------------------------------|---|-------------|-------|
| Symmetric two-site TFI row | 2 | 3 | 1.79 |
| A–T parameter row | 2 | 3 | 1.79 |
| Symmetric three-site TFI row | 3 | 15 | 1.81 |
| G–C parameter row with selected central asymmetry | 3 | 15 | 1.95 |

---

## F86 XY-block Q-peak scan and separate F88b Heisenberg diagnostic

The separate [EP-resonance script](../../simulations/water/proton_chain_ep_resonance.py)
first runs the **XY-only** F86 block scan through `fw.block_L_split_xy`, not
the Heisenberg branch: `N = 5`, a popcount-`(2, 3)` coherence block, all-site
local Z dephasing at `γ₀ = 0.05`, and a selected control scan
`Q = J / γ₀` with spacing `0.025`. Its bond Hamiltonian is
`(J/2)(XX + YY)`; it has no `ZZ` term. This is a model witness for the
block-resolved response, not spectroscopy and not a water measurement.

| Bond class | Q_peak (control, selected run) | HWHM−/Q* (selected run) | F86b Tier-1-candidate reference |
|------------|--------------------------------|--------------------------|---------------------------------|
| Interior (b = 1, 2) | 1.566 | 0.7458 | 0.756 ± 0.005 |
| Endpoint (b = 0, 3) | 2.400 | 0.7663 | 0.770 |

The same script later performs a separate **full state-level Heisenberg +
local-Z-dephasing** diagnostic. There it uses
`H_Heis = (J/4) Σ_b (X_b X_(b+1) + Y_b Y_(b+1) + Z_b Z_(b+1))` with `J = 1`
and the selected local-Z rate `γ = 1 / Q_peak(interior)`. This Pauli
normalization and generator are not those of the XY F86 scan; the copied
numerical control ratio selects the later diagnostic but does not turn its
table into an F86 XY result.

The selected state-level readout has
`(|00011> + |00111>) / √2`: two popcount labels `(2, 3)`, Hamming distance
one, and `n_p + n_q = N`. In the Heisenberg plus all-site-Z-dephasing kernel,
F88b gives its time-zero Π²-odd/memory value `10/19 = 0.5263`; the later rows
are the selected Heisenberg state-evolution reading. These are mathematical
consequences of the named state and generator, not proton signatures.

| t | Static / total | Memory / total | Π²-odd / memory | Per-site \|r\| (0..4) |
|---|----------------|----------------|------------------|--------------------|
| 0 | 0.05 | 0.95 | 0.5263 | 1.00, 1.00, 1.00, 1.00, 1.00 |
| 1 | 0.11 | 0.89 | 0.561 | 1.00, 0.90, 0.22, 0.90, 1.00 |
| 5 | 0.27 | 0.73 | 0.677 | 0.83, 0.56, 0.00, 0.56, 0.83 |
| 20 | 0.82 | 0.18 | 0.957 | 0.29, 0.19, 0.00, 0.19, 0.29 |
| 50 | 1.00 | 0.00 | 1.00 | ~0 |

---

## What these calculations do not establish

- No room-temperature coherent oscillation in ordinary liquid water.
- No Grotthuss transport rate or charged-carrier prediction.
- No calibrated inter-coordinate coupling `K` for water.
- No observed water `CΨ` crossing.
- No direct proton-coordinate `T₂` and no validation that a physical liquid
  water coordinate realizes the selected local-Z channel.

## Reproducibility

| Component | Location |
|-----------|----------|
| Two-branch sweep | [proton_water_chain.py](../../simulations/water/proton_water_chain.py) |
| Stored output | [proton_water_chain.txt](../../simulations/results/proton_water_chain.txt) |
| F86 XY-block scan and F88b Heisenberg diagnostic | [proton_chain_ep_resonance.py](../../simulations/water/proton_chain_ep_resonance.py) |

*From one to five chosen two-level coordinates, these calculations compare an
all-site-dephased Heisenberg formula model with an all-site-dephased,
unbiased-TFI dipole-moving toy model. They report properties of those named
generators and selected parameter rows, not a calibrated account of liquid
water or proton transport.*
