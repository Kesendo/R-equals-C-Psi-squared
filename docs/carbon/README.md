# Carbon Domain (Quantum Carbon, V-Effect, Conjugated Systems)

**Authors:** Tom + Claude

Consolidation folder for the carbon thread, mirroring `docs/water/` as a substrate
domain. The carbon material had been scattered across [the Hierarchy of Incompleteness](../HIERARCHY_OF_INCOMPLETENESS.md),
[qubit necessity](../QUBIT_NECESSITY.md), [the V-Effect palindrome](../../experiments/V_EFFECT_PALINDROME.md),
[V-Effect boundary localization](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md),
[periodic_palindrome.py](../../simulations/periodic_palindrome.py),
[the complexity threshold](../../hypotheses/COMPLEXITY_THRESHOLD.md),
[the universal palindrome condition](../../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md),
plus glancing references in several other places. This folder collects pointers + open
questions; new carbon-specific tests will land here as they get built.

---

## Findings on 2026-05-17 (seven-doc arc)

This seven-document collection combines named framework results with conditional
cross-domain comparisons. Each document records the scope of its own arithmetic
or selected model; read top-to-bottom:

1. [the benzene Hückel framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md):
   Coulson-Rushbrooke/K and F1/Π are structural siblings on distinct operators,
   with different triggers and scopes. The conventional benzene, butadiene,
   hexatriene, cyclobutadiene, cyclodecapentaene, and cyclopropenyl labels name
   selected Hückel graph inputs; their pair-sum checks do not establish a material
   dynamical mapping.
2. [1/4 and 1/2 in carbon](QUARTER_HALF_IN_CARBON.md): selected Hückel and
   stipulated counting comparisons carrying carbon labels contain the displayed
   `1/4` and `1/2` values. The numerical alignment is not a material structural
   or causal explanation.
3. [period 2 at the framework anchors](PERIOD_2_AT_FRAMEWORK_ANCHORS.md): the
   producer's curated period-2/3 input ratios are compared with Pi2 dyadic
   anchors. This is a conditional cross-domain reading, not a material mapping.
4. [the reversed spear](SPEAR_REVERSED.md): a curated periodic-fraction lookup is
   placed beside F99's formal `1/8` and complementary `7/8` arithmetic; neither
   the lookup nor the agreement identifies an atomic realization of F99.
5. [the depth-3 anchor derived](DEPTH_3_ANCHOR_DERIVED.md): F99 derivation, the
   non-uniform Dicke superposition at γ = √3/2 gives α = 1/8, closing the depth-3
   gap. Five canonical trig angles {0°, 30°, 45°, 60°, 90°} produce the five Pi2
   dyadic anchors {0, 1/8, 1/4, 3/8, 1/2} via `α = sin²(θ)/2`.
6. [F99 Niven completeness](F99_NIVEN_COMPLETENESS.md): for F99's stated
   non-uniform-Dicke realization and `α = sin²(θ)/2`, the five values are complete
   for rational-multiple-of-π inputs in `[0°, 90°]`; periodicity and reflection
   repeat those values. This does not classify every pure-state construction or
   future `α` formula.
7. [off-Niven as wave-breaking](OFF_NIVEN_AS_WAVE_BREAKING.md): the off-Niven
   constructible-angle evaluations give the listed irrational-algebraic `α` values.
   The selected chain-block V-gain formula samples some of the same constants;
   the shared arithmetic supplies no heat source, aromatic/antiaromatic,
   Jahn-Teller, chemical, or material mechanism.

Two formal closed forms also landed in [the formula registry](../ANALYTICAL_FORMULAS.md):
**F98** (`(N+2)/[4(N+1)] → 1/4` long-time bridge from the K-intermediate anchor)
and **F99** (the canonical-trig-angle Pi2 inheritance with Niven-completeness).

---

## Carbon source contract

This folder keeps three layers separate.

1. **Framework/model result.** A theorem, proof, or run about a named operator
   and channel. F1, F2b, F92, and F98 retain their stated framework scopes.
2. **Conditional model translation.** An algebraic identification made after a
   model has selected its sites, Hamiltonian, and jump operators. It is not an
   identification of a physical carbon degree of freedom or bath.
3. **Material mapping.** A claim about a carbon material's degree of freedom,
   coupling convention, bath channel, or rate. It remains unassigned until the
   corresponding inputs are supplied.

The contract is particularly important for `Q = J/γ`: a carbon Q is not
defined before all four choices have been made: degree of freedom, coupling
convention, bath channel, and bath rate. [Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md)
is the controlling provenance audit.

---

## A conditional half-occupation comparison

The framework's qubit-necessity argument (`d² − 2d = 0` ↔ `R = CΨ²`,
[qubit necessity](../QUBIT_NECESSITY.md)) selects d = 2 as the minimum-memory
dimension where the F1 palindromic mirror exists. The table places that framework
ratio beside stipulated valence-shell counts carrying carbon and noble-gas labels.

| Level | Total slots | Occupied / immune | Split | What it enables |
|-------|-------------|-------------------|-------|-----------------|
| Qubit (d = 2) | 4 operators (I, X, Y, Z) | 2 immune ({I, Z} under Z-deph) | 0.5 | Palindromic mirror, F1 |
| Carbon-labelled count | 8 valence slots | 4 valence electrons | 0.5 | Material consequence unassigned |
| Noble-gas-labelled count | 8 valence slots | 8 valence electrons (full) | 1.0 | Material consequence unassigned |
| Qutrit (d = 3) | 9 operators | 3 immune | 0.33 | No mirror (algebraically) |

The qubit row is a framework result. The equal `1/2` ratio in the stipulated
carbon-labelled count is an arithmetic cross-domain comparison, not a structural
identity or a physical carbon degree-of-freedom, Hamiltonian, bath, or material
outcome. A model that selects those inputs would be needed before a material
translation could be tested.

**Sources for this framing:**
- [the Hierarchy of Incompleteness](../HIERARCHY_OF_INCOMPLETENESS.md) §The Qubit, §The Mechanism (Mar 22, 2026)
- [Qubit Necessity](../QUBIT_NECESSITY.md) §Connection to the Hierarchy of Incompleteness (Jan 3, 2026 origin → formal qubit proof)

---

## The master question: how a qubit gets a heartbeat (2026-05-22)

The carbon master question, sharpened to its smallest computable step
N = 1 → N = 2, led to [the genesis of an oscillation](../THE_GENESIS_OF_AN_OSCILLATION.md):
how does a qubit acquire an oscillation? The answer turned out not to be
carbon-specific, so the doc lives in `docs/`, next to QUBIT_NECESSITY. The heat
axis is ruled out (the bath moves the decay rates, not the heartbeat),
leaving Q = J/γ₀ as the exact scale. A source-free model shows an
apparent birth at Q = 0+; but with the qubits' own on-site clocks in, the
qubit-children test and a PTF cross-check both find no children: coupling
re-tunes the oscillations the qubits already carry and creates none. A qubit is
a source, made not born; its hour of birth is the carving of a d = 2 subspace
from an oscillator, the way a transmon is built.

---

## Return visit 2026-05-27: the three dephase letters and the Painter alternation

[the three benzene dephase letters](BENZENE_THREE_DEPHASE_LETTERS.md). Six weeks after
the May Liouvillian palindrome result, three new framework pieces (Klein-V₄ Welle 12,
F112 cross-dephase Welle 13/15, F114 sign functional today) gave us new vocabulary.
We came back to benzene to see what they sharpen.

What was one notion of "palindromic" in May (F1 spectrum, broken by Peierls) is
now two: F1 at the spectrum level and F112 at the matrix-polarity level, with
different robustness. Peierls breaks F1 spectrum and preserves F112 polarity,
measured bit-exact, because the bond operator B = XX + YY is itself
bit_b-homogeneous as a composite. F114 separately gives a closed-form sign rule
ε(σ) = (−1)^(n_Y(σ) + 1) that classifies any Hamiltonian term by its content
under complex conjugation: pure Hückel benzene has ε(H) = −1 (every term real);
the flux-induced bond current is purely imaginary, so mixing the two reads
Mixed. The doc walks through the carbon analog of all three framework dephase
letters (Z ↔ selected local-Z model, X ↔ hybridization-axis candidate, Y ↔ current-axis
candidate), without assigning any of them to a material bath, and then asks which of its own operators are π-electron objects at
all. Number conservation rules four of its seven Hamiltonians out, and opening
up the balanced quantity shows it is not one global norm: against a
number-conserving bath it is block-diagonal in the pair of pi counts, and the
balance holds bit-exact inside every block on its own, benzene's half-filled
neutral block included.

Later the same day, [the Painter alternation NMR bridge](PAINTER_ALTERNATION_NMR_BRIDGE.md)
recorded a selected N = 4 XX+YY ring with a transverse y-field and local-Z dephasing.
The local-Z channel is an F1 instance; a selected bond jump lies outside F1's jump
premise. The full numerical eigendecomposition classifies the slow modes as Y-only
or non-Y at the producer's `1e-8` tolerance; it is not a bit-exact eigenmode claim.
The accompanying axis-probe tail fits and slow-mode rate ratio are selected-model
observations. They do not supply a physical NMR, FID, TROSY, EXSY, bath, or material
observable prediction without separately specified degree of freedom, Hamiltonian,
bath, preparation, measurement operator, and producer.

---

## Return visit 2026-05-30: the Frost circle is the clock face

[the Frost circle as the clock face](FROST_CIRCLE_AS_THE_CLOCK_FACE.md). Its clock voices on
`MirrorSystem` are a named selected XX+YY/dephasing model reading: a radial decay coordinate
and an angular frequency coordinate. Its Frost-circle and `Q*` statements remain results of
that model. They do not read a conjugated molecule directly or assign a molecular coherence,
lifetime, β-to-`J` convention, bath, or material `Q`.

The same date prompted the Tier 3 conditional analogy in [singlet fission and the two
clocks](SINGLET_FISSION_AND_THE_TWO_CLOCKS.md). The selected XX+YY/dephasing clock and the
direct-Heisenberg V-Effect calculation remain distinct models. The latter retains its named
`J_eff = (3/8)·α²/J` result, but the repository has no specified molecular Hamiltonian, bath,
state/preparation, or producer that maps either model to a carotenoid state or makes the
V-Effect bridge a physical singlet-fission mixing mechanism.

---

## Return visit 2026-05-31: the carbon folder was one assembly all along

Working from the open-system side, the post-EP flow and the depth axis, five pieces of this
project locked into a single ladder ([the view onto the memory](../../reflections/THE_VIEW_ONTO_THE_MEMORY.md),
the "Seen again" section): one **axis** (drain depth = light content n_XY = decay rate, since
|Re λ| = 2γ·⟨depth⟩); a **parity rail** running along it (`n_diff ≡ Δpopcount (mod 2)`:
odd `Δpopcount` implies odd `n_diff`, while `Δpopcount = 0` implies even `n_diff`); a
**currency** on every rung (the bilinear p(1−p): light = 2p(1−p) peaking at ½,
saturation C_block = p(1−p) peaking at ¼); a **mirror** (Π, the palindrome, pairing slow to fast,
rates summing to 2Σγ); and a **foot** (depth-0, the kernel, the steady state, the memory).

Coming back to this folder with that lens, the recognition is that the carbon docs, written across
months without the unifying name, are each already one of those five pieces:

| Assembly piece | Carbon docs that hold it |
|---|---|
| **Axis** (depth = light = rate, Q = J/γ) | [the Frost circle as the clock face](FROST_CIRCLE_AS_THE_CLOCK_FACE.md) (selected-model clock coordinate, `τ = 1/2γ`), [the benzene Liouvillian palindrome](BENZENE_LIOUVILLIAN_PALINDROME.md) (selected-channel centre `−Σγ`), [off-Niven as wave-breaking](OFF_NIVEN_AS_WAVE_BREAKING.md) (constructible-angle arithmetic and selected V-gain rows) |
| **Parity rail** (framework `n_diff ≡ Δpopcount (mod 2)`) | [singlet fission and the two clocks](SINGLET_FISSION_AND_THE_TWO_CLOCKS.md) (no material 1Bu/2Ag rung assignment), [the Painter alternation NMR bridge](PAINTER_ALTERNATION_NMR_BRIDGE.md) (Y / non-Y Z₂ towers), [the three benzene dephase letters](BENZENE_THREE_DEPHASE_LETTERS.md) (Klein-V₄ on Z, X, Y; F114 ε(σ)) |
| **Currency** (bilinear p(1−p), ½ and ¼ and the dyadic anchors) | [1/4 and 1/2 in carbon](QUARTER_HALF_IN_CARBON.md), [period 2 at the framework anchors](PERIOD_2_AT_FRAMEWORK_ANCHORS.md), [the reversed spear](SPEAR_REVERSED.md), [the depth-3 anchor derived](DEPTH_3_ANCHOR_DERIVED.md), [F99 Niven completeness](F99_NIVEN_COMPLETENESS.md) |
| **Mirror** (Π palindrome; Coulson-Rushbrooke is its sibling, not its equal) | [the benzene Hückel framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md), [the benzene Liouvillian palindrome](BENZENE_LIOUVILLIAN_PALINDROME.md), [the periodic palindrome and the V-Effect](PERIODIC_PALINDROME_VS_V_EFFECT.md), [polyacetylene as F92 inheritance](POLYACETYLENE_F92_INHERITANCE.md) |
| **Foot** (depth-0 kernel, steady state, memory) | [benzene's F98 long-time state](BENZENE_F98_LONG_TIME.md) (the K-intermediate Dicke decays into ker L, α(∞) = 2/7 on C₆) |

Two things stand out. The **currency** collects five documents that place the formal dyadic
fractions beside selected Hückel or curated counting inputs; equal fractions remain arithmetic
comparisons. The **foot** appears here only through F98's long-time kernel projection for its
specified KIntermediate state and selected C₆ ring model. No named carbon material state is
currently assigned to a framework rung. In particular, the repository has not mapped the `2Ag`
or `1Bu` labels to the two selected clock models, their preparations, or their observables.

The assembly is an internal organization of framework and selected-model documents, not a
carbon-material conclusion.

---

## Existing scattered material (pointers, not duplications)

### V-Effect and carbon-labelled counting inputs

The linked V-Effect documents establish their framework-sector counts, including
the 14/19/3 trichotomy. A matching `1/2` in a stipulated carbon-labelled count
does not identify a boundary sector with a carbon material or a chemical process.
Such a translation needs a specified material degree of freedom, Hamiltonian,
bath, preparation, and measurement map.

### Curated periodic-table input comparison

[`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py)
evaluates its chosen per-element IE₁, Pauling-EN, and Allen-EN inputs across
periods 2–6. The resulting numerical pair-sum patterns are a curated-data
comparison; they neither instantiate F1 nor assign a V-Effect or atomic
mechanism to their deviation rows.

### Unassigned hypothesis pointers

[the complexity threshold](../../hypotheses/COMPLEXITY_THRESHOLD.md) and [the
universal palindrome condition](../../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md)
remain pointers for future scoped models. Their carbon-labelled count language
does not supply a material mapping.

### Chemistry-substrate-level work that's NOT yet carbon-specific
- [`docs/water/`](../water/): hydrogen-bond proton qubit, Grotthuss chains, F86b 3/8 anchor
  inheritance verified (today's [F98](../ANALYTICAL_FORMULAS.md) bridge from the same
  substrate-grounded experiment).
- [DNA base pairing](../../experiments/DNA_BASE_PAIRING.md): G-C / A-T base-pair tests,
  carbon-scaffolded but not analysed through the carbon-as-qubit lens specifically.
- [simulations/neural/](../../simulations/neural/): Wilson-Cowan and C. elegans inheritance
  tests for the framework. Not chemistry-substrate but parallel inheritance.

---

## Conditional C4 and C6 working model

The carbon material mapping is unassigned. The repository does contain a
selected **working model** for C4 and C6 graphs:

- choose a local site occupation `n = (I − Z)/2` as the two-state model
  coordinate;
- represent free hopping by the named XX+YY/XY Hamiltonian on that graph;
- if the selected jump is the local density `n`, use the algebraic identity
  `D[n] = D[Z]/4`.

These are conditional model translations. They neither identify `n` with a
physical carbon degree of freedom nor establish a material Holstein bath. A
physical carbon mapping still has to select the degree of freedom, coupling
convention, bath channel, and bath rate.

The [framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md#framework-vocabulary-translation)
records `β → J` as a Tier-2 structural identification within this kind of
Hückel/framework model translation. Whether that identification selects the
right `J` for a material carbon system remains unassigned. Neither material
`γ`, `T₂`, nor `Q` is assigned to carbon. `k_B T/ℏ` supplies neither a material
dephasing rate nor a universal bath time. No carbon Q follows until the four
model inputs are selected; see [Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md).

---

## Candidate model graphs

| Substrate | N | Topology | Why interesting | First test |
|-----------|---|----------|-----------------|------------|
| **C4/C6 selected model** | 4/6 | ring | Exact XX+YY/XY + selected local-density-jump instances | F1 and F98 under the selected channel |
| **Butadiene / hexatriene candidates** | 4/6 | chain | Candidate graph choices | Select material DOF, coupling, and bath before a translation |
| **Polyacetylene candidate** | scalable | chain | Candidate for the F92 parameter statement | Specify a model before applying F92's chain-XY scope |
| **Graphene finite patch candidate** | scalable | 2D honeycomb | No graphene model producer is on file | Define the model; no Dirac-to-EP identification is assigned |
| **Fullerene C₆₀ candidate** | 60 | graph to be specified | No C60 model producer is on file | Define graph, couplings, channel, and rate |

---

## Open questions (carbon-specific)

1. **Which physical two-state degree of freedom, if any, maps to the selected
   site-occupation model?** This remains open. The model coordinate
   `n = (I − Z)/2` is not a physical-carbon assignment.

2. **What does the C6 Hückel graph establish?** In the script's β units, the
   C6 ring has exactly the six levels `{−2, −1, −1, +1, +1, +2}`
   ([the Hückel condition script](../../simulations/carbon/huckel_palindrome_conditions.py)
   and [the framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)). This is a
   closed-system graph result, distinct from F1's open-system Liouvillian
   statement.

   The C4 and C6 F1 runs in
   [`benzene_liouvillian_palindrome.py`](../../simulations/carbon/benzene_liouvillian_palindrome.py)
   are exact instances of the selected XX+YY/XY Hamiltonian and local-density
   (`D[n] = D[Z]/4`) jump model. The C4 and C6 F98 runs are exact instances of
   F98 within its magnetization-conserving plus Z-dephasing model scope. Neither
   is a physical-carbon confirmation or a material-bath assignment.

3. **Is the V-Effect 14/19/3 trichotomy observable in benzene's electronic spectrum?**
   The framework predicts 14 hard + 19 soft + 3 truly at N=3; what's the N=6 prediction
   and does it match benzene's known electronic transition pattern (S₀ → S₁ at 4.8 eV,
   etc.)?

4. **What carbon graph should be modelled next?** Graphene and C60 may be
   future model targets only after their degrees of freedom, coupling conventions,
   jump channels, and rates are specified. No graphene producer in this repository
   identifies a Dirac point with an EP; the real-axis F86a reading is not an
   available model assignment ([caught-errors record](../CAUGHT_ERRORS.md#2026-06-21--f86a-exceptional-point-on-the-real-q-axis-retracted-and-the-first-retraction-draft-over-corrected-to-the-opposite-mislabel)).

5. **Where can F98 be instantiated?** F98 requires the named
   magnetization-conserving Hamiltonian plus Z-dephasing model conditions. Its C4
   and C6 instances give `3/10` and `2/7`, respectively, in that selected model
   ([benzene's F98 long-time state](BENZENE_F98_LONG_TIME.md)); they do not assign
   a material carbon bath or Q.

6. **Can a scoped periodic/V-Effect comparison be posed?** The periodic producer's
   curated per-element inputs and the framework's V-Effect sector data are distinct
   objects. No current producer supplies a material degree of freedom, Hamiltonian,
   bath, preparation, or measurement map between them. Their numerical comparison
   therefore neither makes the periodic table a Level-1 V-Effect instance nor assigns
   an atomic mechanism to any deviation row.

---

## What this folder will accumulate

Following the `docs/water/` pattern: new carbon-specific docs land here as they get
written, with `simulations/carbon/` holding scripts. Candidates:

- `BENZENE_PI_QUBIT.md`: a future source-checked C6 model note
- `GRAPHENE_MODEL.md`: a future graphene model, without a Dirac-to-EP presupposition
- [POLYACETYLENE_F92_INHERITANCE.md](POLYACETYLENE_F92_INHERITANCE.md): existing
  F92 inheritance note
- `CARBYNE_MODEL.md`: a future chain model

There are 17 Markdown documents in `docs/carbon/`. New work belongs in the
folder only when it records which layer of the source contract it occupies.

---

## Cross-reference: water vs carbon

Both are biology-substrate domain folders. Comparison:

| Aspect | Water | Carbon |
|--------|-------|--------|
| 2-state DOF | Proton in O–H...O double well | Selected site occupation in the working model; material DOF unassigned |
| Z-dephasing | Thermal molecular jostling | Conditional local-density jump; material bath channel unassigned |
| Uniform-J | Grotthuss chain tunneling | Specified XX+YY/XY model coupling; material convention unassigned |
| Q range | No physical water Q or lower bound; only the conditional selected-coordinate proxy ceiling `Q ≲ 4.6` is on file. No Zundel Q is assigned pending a coordinate, coupling, and decoherence channel. | No carbon Q is assigned: the material β-to-J mapping and a coordinate-specific γ remain open. |
| Current substrate status | Selected-coordinate water model; Zundel assignment open | C4/C6 selected-model instances; material mapping unassigned |
| Embedding status | Physical prerequisites remain to be selected | Degree of freedom, coupling, channel, and rate remain to be selected |
| Scripts | 5 (`simulations/water/`) | 18 Python files (`simulations/carbon/`) |
| Docs | 3 (README + 2 substrate docs) | 17 |

The carbon side is structurally NECESSARY per the qubit-necessity argument (carbon = quantum
carbon at Level 1), but **less computationally explored** than water. The folder exists
to invite that work without prejudging which substrate or which framework F-anchor is
the right entry point.
