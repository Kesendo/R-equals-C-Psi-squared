# The Hydrogen Bond as a Qubit

<!-- Keywords: hydrogen bond proton qubit tunneling, double-well potential
palindromic spectral symmetry, CΨ crossing water enzyme, V-Effect hydrogen
bond network, proton transfer fold catastrophe, R=CPsi2 hydrogen bond -->

**Status:** Tier 2 (computed from proven framework)
**Date:** March 28, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../proofs/MIRROR_SYMMETRY_PROOF.md),
[CΨ Monotonicity](../proofs/PROOF_MONOTONICITY_CPSI.md)

---

## What this document is about

This document makes a modelling map, not an identity claim about matter.
For a selected proton-position coordinate, choose two localized basis states
`|L⟩` and `|R⟩` and use

`H = -J σ_X + Δ σ_Z`.

Here `J` is the selected tunnelling matrix element (coupling) and `Δ` is the
longitudinal bias of that coordinate. At `Δ = 0`, this Hamiltonian has energy
splitting `2|J|`; a source-reported tunnelling splitting must therefore not be
inserted automatically as the model parameter `J`. A selected local
Z-dephasing channel is then an effective environmental model. Physical protons
and abstract qubits are not thereby identical, and whether this two-level,
Markovian channel is useful is an empirical modelling question.

The F1 palindrome calculations discussed below use the un-biased mapping
`Δ = 0` together with their stated Hamiltonian and local-dephasing premises.
A longitudinal `Δ σ_Z` term is outside that un-biased mapping. That does not
license the false shortcut that every field breaks F1: field and channel
directions must be assessed against the named F1 premises.

For liquid water, the repository gives no Q. It records only an illustrative
conditional ceiling, `Q ≲ 4.6`, when a proton coordinate, the ice-derived
`J = 0.5 meV` convention, and a 1–3 ps H-bond-lifetime proxy for its missing
`T₂` are all stipulated. It has no lower endpoint and is not a diagnosis of
ordinary water or a statement about every proton-transfer event.

---

## Abstract

For a chosen O-H...O coordinate, the model keeps `|L⟩` (donor-side) and
`|R⟩` (acceptor-side) as a two-level basis and starts from
`H = -J σ_X + Δ σ_Z`. It selects local Z-dephasing as an effective channel;
this is a modelling choice, not a claim that a physical proton is identical
to an abstract qubit.

The F1 palindrome is exact for the theorem's named generator and channel
premises. The 87,376-eigenvalue calculation is a finite verification of that
mathematical model, not a transfer of the theorem to arbitrary hydrogen bonds.
The displayed water-inspired rows are calculations at selected parameters.

**Results:**
- Selected single-coordinate runs: CΨ crosses 1/4 at 0.07-1.32 ps when
  `J/γ = 1` in the displayed model.
- Selected two-coordinate run: a finite numerical pairing check is small and
  CΨ crosses at 0.46 ps.
- Selected four-coordinate run: coupling changes the calculated frequency
  count from 11 per isolated model molecule to 126 for the coupled model.
- Illustrative conditional ceiling: `Q ≲ 4.6` under the stated coordinate,
  ice-derived-`J`, and lifetime-proxy assumptions; it assigns no Q to liquid
  water and does not classify it as classical or quantum.

---

## Limits of an Earlier Classical Model

A previous attempt (V17 negative result, local analysis not published)
modeled water classically: hydrogen bonds as springs with friction,
donor/acceptor modes as coupled oscillators. Palindrome residual: 1.33.
No palindromic structure found.

This only diagnoses a limitation of that selected classical model: its
edge-based donor/acceptor variables did not supply the coordinate used by its
palindrome test. It neither proves that water is classical nor establishes
that the two-level map is universally valid.

The alternative model selects a proton-position coordinate with `|L⟩` and
`|R⟩`, uses `H = -J σ_X + Δ σ_Z`, and chooses a local Z-dephasing channel.
With `Δ = 0` and the other F1 premises, the theorem supplies the palindrome;
outside those premises the model must be checked on its own terms.

---

## The Model

### Single selected coordinate (N=1)

Hamiltonian: H = -J · σ_X + Δ · σ_Z

- J: selected tunnelling matrix element; at `Δ = 0`, the energy splitting is
  `2|J|`
- Δ: asymmetry of the double well (0 for symmetric bonds)

Dephasing: L_k = √γ_eff · σ_Z (where γ_eff = γ/ℏ in angular frequency)

- γ: selected dephasing energy scale in the model (in eV); its relation to
  surrounding molecules requires a specified bath model or measurement

Initial state: |L⟩ = |0⟩ (proton on donor side, no coherence)

The repository observable is `CΨ = Tr(ρ²) · L1(ρ) / (d - 1)`: purity times
normalized L1 coherence. In the selected initial-state runs below it starts at
zero, can rise as the model Hamiltonian creates coherence, and may or may not
reach the reference value 1/4 before the selected dephasing suppresses it.

The F1 calculations listed below take `Δ = 0`. A nonzero longitudinal bias is
not part of that un-biased palindrome mapping; it requires a separate
generator-level analysis rather than a blanket conclusion about fields.

### Three regimes

| Selected model J/γ | Model-run label | CΨ crosses 1/4 in this model? | Interpretation |
|--------------------|-----------------|-------------------------------|----------------|
| << 1 | classical | No (overdamped) | selected model run, not a water diagnosis |
| ~ 1 | **fold** | **Yes (sub-ps)** | selected model run |
| >> 1 | quantum | Yes (slower) | selected model run |

### Two selected coordinates (N=2): a one-molecule model

This model represents two chosen O-H coordinates coupled through a shared
oxygen. It is not a literal assertion that a water molecule contains two
independent qubits. Its selected Hamiltonian is:

```
H = -J · σ_X(1) - J · σ_X(2) + K · σ_Z(1) · σ_Z(2)
```

Dephasing on both qubits: γ · σ_Z(1) and γ · σ_Z(2).

### Four selected coordinates (N=4): a two-molecule model

```
Molecule 1       H-bond       Molecule 2
H(1)-O ... H(2)---O ... H(3)-O-H(4)
            donor  M  acceptor
```

Coordinate 2 is assigned to the donating side of the selected intermolecular
bond. Coordinates 1,2 and 3,4 are coupled intramolecularly by `J_intra`; the
selected intermolecular coupling is `J_inter`. These are model assignments.

---

## Results

The following are selected model runs. Their regime labels describe the
displayed `J/γ` rows, not the state of water, ice, enzymes, or chemistry.

### Phase 1: Single-coordinate model

| J (meV) | J/γ | CΨ crossing time | Regime |
|---------|-----|------------------|--------|
| 0.5 | 0.01 | no crossing | classical |
| 0.5 | 1.0 | 1.32 ps | fold |
| 1.0 | 1.0 | 0.66 ps | fold |
| 5.0 | 1.0 | 0.13 ps | fold |
| 10.0 | 1.0 | 0.07 ps | fold |

### Phase 2: One-molecule model (N=2)

Parameters: J_intra = 1.0 meV, K = 0.1 meV, γ = 1.0 meV (J/γ = 1).
This γ is chosen to place the run at the fold, not derived from water; it is
neither the 25 meV row below nor a measured rate.

- Finite numerical pairing check: pair-sum std = 5.4e-3 relative to a mean
  of ~6e12. This small residual is a model check, not a proof or a replacement
  for the exact F1 theorem under its named premises.
- Distinct frequencies: 11
- CΨ crosses 1/4 at **0.46 ps**

### Phase 3: Two-molecule hydrogen-bond model (N=4)

Parameters: J_intra = 1.0 meV, J_inter = 0.1 meV, K = 0.1 meV,
γ = 1.0 meV.

- Finite numerical pairing check: pair-sum std = 3.5e-2; it is not an
  assertion of theorem-level exactness.
- Distinct frequencies: **126**
- V-Effect: 11 per molecule → 126 coupled = **104 new frequencies**

Full-system CΨ does not cross 1/4 (N-scaling suppression, d-1 = 15).
Subsystem CΨ for the selected pair (q2, q3) across the H-bond would need a
separate computation; no local subsystem-crossing document is currently
available.

---

## The V-Effect in Hydrogen Bonds

**Continued in:** [Proton Water Chain](PROTON_WATER_CHAIN.md) (selected neutral
1-D proton-coordinate chain N=1-5, formula validation, frequency explosion
0→222) and
[DNA Base Pairing](../../experiments/DNA_BASE_PAIRING.md) (A-T N=2, G-C N=3, sacrifice
zone in base pairs). At biological temperature (310 K) those selected DNA
model runs put the displayed H-bonds at J/γ ~ 0.01 with all modes overdamped.
That 0.01 is a floor from a γ whose stated source and conversion are
insufficient; it is not directly comparable with a water Q (see
[Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md)).

| System | Distinct frequencies |
|--------|---------------------|
| Single molecule (N=2) | 11 |
| Two molecules coupled (N=4) | 126 |
| New from coupling | **104** |
| Ratio | 5.73 |

For this N=4 model and its displayed parameters, coupling changes the
calculated frequency count by 104 relative to the two selected isolated-model
counts. This is the model's V-Effect reading, not a count established for
physical hydrogen bonds generally.

Note: the abstract qubit V-Effect (N=5 MediatorBridge) produces 109
total frequencies. The numbers differ: 104 NEW at N=4 H-bond vs 109
at N=5 abstract.

---

## Physical Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Tunnelling splitting (ice) | 0.2-1 meV | Bove et al 2009 |
| Tunneling splitting (strong H-bond) | 1-10 meV | Cleland & Kreevoy 1994 |
| O...O distance (normal) | 2.7-3.0 A | Steiner 2002 |
| O...O distance (strong) | 2.4-2.6 A | Cleland & Kreevoy 1994 |
| H-bond lifetime (liquid water) | 1-3 ps | Luzar & Chandler 1996 |
| Thermal energy (300 K) | k_B T ≈ 25 meV | Standard |

The 0.2-1 meV row is an ice result. The `J = 0.5 meV` row used in the model
runs is an ice-derived convention, not a direct identification of a reported
splitting with the matrix element `J`; its mapping and validity for liquid
water are unverified. `k_B T ≈ 25 meV` is a thermal energy, not a measured
proton-coordinate dephasing rate or a universal bath-correlation time. No bath
spectral density, cutoff, or system-bath coupling is specified here, so this
energy cannot be converted into the model rate `γ`.

**Conditional estimate only.** The provenance derivation in
[Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md) gives the
illustrative conditional ceiling `Q ≲ 4.6` if a chosen proton coordinate, the
ice-derived `J = 0.5 meV` convention, and a local-decoherence channel support
`Q = 2JT₂/ℏ`, with the 1–3 ps H-bond lifetime used as a proxy for the missing
proton-coordinate `T₂`. A direct proton-coordinate `T₂` measurement and a
microscopic bath model are absent, so this has no lower endpoint and assigns no
Q to ordinary liquid water.

Whether an enzyme environment realizes a specified coordinate, coupling, and
channel near a selected model-run ratio remains an empirical question.

---

## Connection to the Framework

The table records the modelling correspondence; it does not identify a
hydrogen-bond proton with an abstract qubit or establish these parameters in
liquid water:

| Framework concept | H-bond realization |
|-------------------|-------------------|
| d = 2 | \|L⟩ (donor) and \|R⟩ (acceptor) |
| σ_X coupling (J) | selected tunnelling coupling in the coordinate model |
| σ_Z dephasing (γ) | selected effective environmental channel |
| CΨ = 1/4 crossing | the repository observable reaches its reference value in a model trajectory |
| P(L) = P(R) | an additional symmetric-population condition of a selected state, not the definition of CΨ |
| V-Effect | selected-model frequency-count change under coupling |
| Concentrator | Tier 4 hypothesis outside this document |

---

## Zundel Cation: an Open Parameter Question

The repository has not established the older 124-meV value as a Zundel
tunnelling splitting. It may instead be a shared-proton vibrational
fundamental. Primary literature is needed to identify the relevant coordinate,
the appropriate two-level reduction if any, its coupling, and its decoherence
channel.

Until that evidence is assessed, this repository cannot derive a Zundel `Q`, a
fold-crossing count, a per-molecule or per-drop rate, or a chemistry conclusion
from that value. The honest open question is whether a literature-supported
Zundel coordinate admits a useful two-level local-dephasing model and, if so,
which F1 premises it satisfies.

---

## Open Questions

1. Does a measured proton-position coordinate support a two-level reduction
   and a local Z-dephasing channel in liquid water?
2. What is the proton-coordinate `T₂`, rather than a bond-lifetime proxy?
3. Which literature-supported Zundel energy is a tunnelling splitting versus a
   shared-proton vibrational fundamental?
4. For a selected molecular generator, which F1 premises hold, including the
   treatment of longitudinal bias?
5. How do the selected-model V-Effect counts scale with the number of
   coordinates?

---

## Scripts

| Script | What it computes |
|--------|-----------------|
| [hydrogen_bond_qubit.py](../../simulations/water/hydrogen_bond_qubit.py) | Phases 1-3 selected-model calculations; its greedy pairing output is not a proof |
| [hydrogen_bond_palindrome.py](../../simulations/water/hydrogen_bond_palindrome.py) | V17 selected classical-model calculation (negative result) |

---

*See also:* [Mirror Symmetry Proof](../proofs/MIRROR_SYMMETRY_PROOF.md)
(the F1 theorem and its stated premises). The selected-model V-Effect and the
open subsystem-crossing question are described above; no corresponding local
legacy documents are currently present.
