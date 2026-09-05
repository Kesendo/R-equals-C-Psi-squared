# Q Belongs to No Substance

*last refreshed 2026-08-22 (the change history lives in git)*

A substance does not hand over a Q. Q = J/γ₀ is defined only once a two-level
degree of freedom, a coupling, and a dephasing channel have been chosen, and a
real material offers several of each. "The Q of water" is therefore
under-determined by water: it is a property of a system in a condition and a
channel, not of a substance.

Within a fixed homogeneous model, γ₀ supplies the unit in which the
Liouvillian is read. At fixed Q, changing γ₀ rescales physical time and energy
while leaving the dimensionless/rescaled dynamics Q-controlled. The scale
result is therefore not an assertion that γ₀ has no physical role. The
repository owned, proved, and gated that distinction in three lineages that
were never in one room. This document puts them in one room and adds what
follows for the substrate folders, where the choice of J and γ was never made
explicitly and two denominators were imported without a source and reused.

## What the sweep returned

`docs/ANALYTICAL_FORMULAS.md` returned the bit-exact γ₀-invariance of Q_peak
and HWHM_left/Q_peak across γ₀ ∈ {0.025, 0.05, 0.10} (`:2870`), the measurement
half of the scale claim. `experiments/` returned the theorem half twice:
[`Q_SCALE_THREE_BANDS.md`](../experiments/Q_SCALE_THREE_BANDS.md) (Tier 1,
2026-04-22) owns the rescaling identity `L(λ·J, λ·γ₀) = λ·L(J, γ₀)` at `:83`,
with a λ-sweep spanning a factor 100 in γ₀ at N=5 (`:66-81`) in which two of the
four columns agree to three decimals; the source explains the γ₀ = 0.01 scatter
as a tolerance artefact and leaves the γ₀ = 0.05 column, which deviates by up to
about 8 %, unexplained. The identity itself is at
[`GAMMA0_IS_ALWAYS_THERE.md`](../experiments/GAMMA0_IS_ALWAYS_THERE.md) states
it as "it is the unit, and in any ratio the unit drops out" (`:26`), scoped
there to a fixed chip rather than stated as the general theorem.
`docs/proofs/` returned [`INCOMPLETENESS_PROOF.md`](proofs/INCOMPLETENESS_PROOF.md)
`:272`: γ defines the scale against which t is counted, which makes the unit
circular by construction. `docs/` itself returned
[`THE_GENESIS_OF_AN_OSCILLATION.md`](THE_GENESIS_OF_AN_OSCILLATION.md) `:49`,
a section headed "Q is the scale, exactly", and two more an earlier draft of
this sweep missed although the same session was repairing their index entries:
[`Q_REGIME_ANCHORS.md`](Q_REGIME_ANCHORS.md) `:98-107`, a section headed
"γ₀ as time-tick", and
[`THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md`](THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md)
`:75`, "γ₀-independent by the very definition of γ₀ being the carrier: the
uniform part cancels in any ratio". `hypotheses/` returned
[`GAMMA_IS_LIGHT.md`](../hypotheses/GAMMA_IS_LIGHT.md) `:298` (the absolute
scale always cancels; only Q shows from within) and
[`Q_AS_THE_EXCHANGE_RATE.md`](../hypotheses/Q_AS_THE_EXCHANGE_RATE.md) `:16`,
whose primary clause is object 1 exactly ("how many H-clock ticks correspond to
one γ₀-clock tick"); only its "equivalently" restatement in **periods** is loose,
and its own table at `:29-37` lists the count as Q, so the looseness is in the
word and not in the number. The same wording sits at `Q_REGIME_ANCHORS.md:107`. `reflections/` returned
[`ON_HOW_GAMMA_BECAME_THE_TICK.md`](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md),
where one tick = 1/γ₀.

`docs/GLOSSARY.md` returned the T₂ → γ conversion in full (`:208`), including
the structural discriminator between the dephasing-only and the T₁-aware form,
and the whole q-versus-Q section at `:317-334`, which already called J/γ₀ "the
central control parameter" (`:319`) and already carried the q = Q/2 trap
(`:326`). What it did
**not** hold was the factorisation L = γ₀·L₁(Q) or any consequence of it for a
substance; the paragraph now at `:336` was written in the same pass and points
here.

`docs/CAUGHT_ERRORS.md` returned three normalization traps (`:566`, `:761`,
`:807`) and, at `:785`, an entry that is entirely about a dephasing rate's
provenance: "γ = 0.5 is not an operating point, it is a test constant that
inherited its way into three projects". That is this document's shape one book
down, and an earlier draft of this sweep reported the store as returning
nothing.

The OpenArcs registry returned `gamma_book_enforced_nowhere` (`:7415`, open),
whose finding is the same shape one book down: the knowledge is complete, the
enforcement is absent. That arc **is** the units arc, so the claim an earlier
draft of this sweep made, that no arc on units existed, was contradicted by the
sentence before it. It also returned `whirlpool_carbon_layers` (`:2596`, open
since 2026-06-03), parked at "water adaptation done; carbon layers and a water
prose note never written", which is the substrate-mapping arc for exactly the
two folders audited here. What was genuinely absent is an arc on Q as a scale
and on substrate provenance; it exists now as `substrate_q_provenance`, opened
with this document's third pass, and the open items below are carried there in
full rather than only here. `fw.Confirmations` returned `gamma_0_marrakesh_calibration`,
a fitted γ_Z that is model-degenerate (0.05 continuous, 0.1 with Trotter n=3
modeled, the same data through two physics models), and
`gamma0_off_the_lever_kingston_may2026` (`confirmations.py:242`), the hardware
read-off of γ₀ from its only lever J, which is the measured face of the unit
claim below. It returned **no** entry in which a substrate Q was measured, and
that narrower statement is the one this document rests on.

`recovered/` returned nothing. `hypotheses/archive/` returned nothing beyond a
single line in the retired `GAMMA_TO_GRAVITY.md`.

## γ₀ is the unit

For the uniform Z-dephasing model with an H homogeneous of degree 1 in J, the
formalism carries two rates and one dimensionless ratio. H sets J, the
dissipator sets γ, and the Liouvillian factors into a unit times a function of
their ratio:

> L(J, γ₀) = γ₀ · L₁(Q),  Q = J/γ₀

(`THE_GENESIS_OF_AN_OSCILLATION.md:53`; the rescaling identity it rests on is
`Q_SCALE_THREE_BANDS.md:83`, and the numerical gate is
`ANALYTICAL_FORMULAS.md:2883`.) At fixed Q every eigenvalue is proportional to
γ₀: γ₀ sets the physical time and energy scale, while the rescaled spectral
shape and dimensionless dynamics are Q-controlled. At fixed J, changing γ₀
also moves Q and reshapes L₁(Q). So `docs/Q_REGIME_ANCHORS.md:136` is stating a
property of a unit when it says γ₀ = 0.05 "was chosen as a convenient round
number, not a physical constant".

In the repository's convention a `D[Z]` channel at rate γ decays coherences at
2γ, so γ = 1/(2T₂) (`GLOSSARY.md:208`). Written with times and an energy J:

> Q = 2 · J · T₂ / ℏ

**Watch the constants.** Q is a ratio of rates, not a count of turns of any
particular mode: a mode of the (0,1) block turns at 2Jμ_m, so its radians per
tick are 2μ_m·Q and no single conversion factor exists. Where the repository
says "Q rotation periods per tick" the word is loose and the number is Q. And
the factor 2
above is the Lindblad book; in the coherence book (γ as the coherence-decay
rate itself) it is absent. The open arc `gamma_book_enforced_nowhere` exists
because that seam is unenforced repo-wide.

**Language note.** γ sits at the emitting end. The chain is not looked at, it
stands in light (`docs/quantum/THE_LABEL_MAP.md:120`; working record: the arc
`gamma_is_the_sender_not_the_watching`). Q counts turning per tick, not turning
per glance.

## A thermal energy is not a dephasing rate

Pure Z-dephasing is unital. It leaves every density matrix diagonal in the Z
basis fixed and removes Z-basis off-diagonals; it therefore does not by itself
select `I/d`, a Gibbs inverse temperature β, or a unique steady state. Once a
Hamiltonian is added, the stationary set is a property of the full generator.
For example, the selected finite unbiased-TFI runs in `docs/water/` reach
`I/d`, whereas the all-site-dephased Heisenberg branch retains the F4 sector
structure. Neither finite model assigns a bath temperature to the pure-Z
channel.

No microscopic bath model, spectral density, cutoff, or system-bath coupling is
specified for the selected proton coordinate. Thus `k_B T/ℏ` is neither its
dephasing rate nor a universal bath-correlation time. A substrate γ may depend
on temperature in a derived model, but it must be measured for that coordinate
or calculated from a specified spectral density and coupling.

The substrate documents historically divided by a temperature-attributed
quantity from a single table row in the pre-repair baseline
`ec7cc619fa075d82137698cedee27b742b7dd6fc:docs/water/HYDROGEN_BOND_QUBIT.md:209`:

> `| Thermal decoherence (300K, upper bound) | γ ~ kT/ℏ ~ 25 meV | Standard |`

Its source column reads "Standard". A real pure-dephasing rate depends on the
system-bath coupling strength and can vanish as that coupling vanishes; a
coupling-independent thermal energy cannot supply that rate. Nor can the
Born-Markov time-scale separation be inferred from temperature alone.

The second temperature-attributed rate is not in a document at all. It is a
printed line in a runner, `simulations/dna_base_pairing.py:195`, giving γ_deph
as 10 to 100 cm⁻¹ for the "molecular environment at 310 K", with no functional
form behind the attribution and no source. It is audited two sections down,
because what it turns out to be is not a second estimate.

## Where each substrate Q came from

| System | Q as written | Numerator | Denominator | Site |
|---|---|---|---|---|
| liquid water, 300 K | 0.02 (historical direct-liquid-water value) | J ~ 0.5 meV: a round value inside the ice-scoped cited row (Bove 2009, pre-repair baseline `ec7cc619fa075d82137698cedee27b742b7dd6fc:docs/water/HYDROGEN_BOND_QUBIT.md:204`), applied to liquid water without comment | the 25 meV estimate | pre-repair baseline `ec7cc619fa075d82137698cedee27b742b7dd6fc:docs/water/HYDROGEN_BOND_QUBIT.md:211` |
| Zundel cation | 4.8 (historical; no current Zundel-Q assignment) | J = 124 meV, no source. The pre-repair baseline `ec7cc619fa075d82137698cedee27b742b7dd6fc:docs/water/PROTON_WIRE_CROSSING.md:244-251` identifies 124 meV = 1000.1 cm⁻¹ as the H₅O₂⁺ shared-proton stretch fundamental, i.e. a vibrational quantum rather than an inter-site coupling | the 25 meV estimate | pre-repair baseline `ec7cc619fa075d82137698cedee27b742b7dd6fc:docs/water/HYDROGEN_BOND_QUBIT.md:255` |
| DNA base pair | 0.01 | J = 0.5 cm⁻¹ | γ = 50 cm⁻¹ ≈ 6.2 meV, and **not** the 25 meV estimate. Both come from the cm⁻¹ set below | `experiments/DNA_BASE_PAIRING.md:96` (J), `:97` (γ) |
| liquid water again, 300 K | 0.01 | the same J = 0.5 cm⁻¹ | the same γ = 50 cm⁻¹ | `simulations/water/proton_water_chain.py:280-285` |
| enzyme active site | ~1 | J ~ 0.5 meV, water's numerator | γ ≈ 0.5 meV, a stipulated 50× reduction of the 25 meV estimate | `hypotheses/PROTEIN_AS_CONCENTRATOR.md:71` |
| π-conjugated carbon, 300 K | ~100 (historical energy quotient; no current carbon-Q assignment) | Hückel β ≈ 2.4 eV, written without a source citation in [the framework lens](carbon/BENZENE_HUCKEL_FRAMEWORK_LENS.md#framework-vocabulary-translation) | the 25 meV thermal energy, not an established γ | [carbon's conditional C4 and C6 working model](carbon/README.md#conditional-c4-and-c6-working-model) |

Dividing two energies is not itself a defect: ℏ cancels out of a ratio of
rates, so E_J/E_γ is the same number as (E_J/ℏ)/(E_γ/ℏ). The defect is
provenance, and there are **two** denominators, not one. Four rows divide by the
25 meV estimate. The other two divide by 50 cm⁻¹, which is a second borrowed
denominator, borrowed the same way.

## The second denominator, and the source it does not have

`simulations/dna_base_pairing.py:193` heads its parameter list "Single H-bond
parameters (from `HYDROGEN_BOND_QUBIT.md`)" and gives J_tunnel as 0.01 to
100 cm⁻¹ and γ_deph as 10 to 100 cm⁻¹, "molecular environment at 310 K".
`HYDROGEN_BOND_QUBIT.md` contains no cm⁻¹ at all. What it contains is
0.2 to 1 meV and 25 meV, which are **1.6 to 8.1 cm⁻¹** and **201.6 cm⁻¹**.
`simulations/water/proton_water_chain.py:280-285` then carries the identical
four-parameter set, value for value.

| | the meV table | the cm⁻¹ set | the meV value in cm⁻¹ | apart by |
|---|---|---|---|---|
| J, ordinary H-bond | 0.5 meV | 0.5 cm⁻¹ | 4.03 cm⁻¹ | 8.1 |
| γ | 25 meV | 50 cm⁻¹ | 201.6 cm⁻¹ | 4.0 |
| J, Zundel | 124 meV | 250 cm⁻¹ | 1000 cm⁻¹ | 4.0 |
| K, inter-bond | 0.1 meV (pre-repair baseline `ec7cc619fa075d82137698cedee27b742b7dd6fc:docs/water/HYDROGEN_BOND_QUBIT.md:149`) | 20 cm⁻¹ | 0.81 cm⁻¹ | 24.8 |

**What the table shows, and what it does not.** The matching displacement of the
γ and Zundel rows can make historical quotients numerically similar after a unit
conversion. It does not identify a two-level Zundel coupling: 124 meV is the
shared-proton stretch fundamental, and the 250 cm⁻¹ candidate is unsourced.
Neither quotient therefore creates a current physical Zundel-Q assignment. For
the ordinary-water rows, converting the entries gives the historical energy
quotient 4.03/201.6 ≈ 0.02; that arithmetic likewise leaves its coordinate,
coupling, and decoherence channel unresolved.

What the table does **not** show is a mechanism, and an earlier version of this
section claimed one: values carried across a unit change without converting. The
fourth row is why that reading fails. K is displaced by 24.8 and fits neither
factor, and it was omitted from the first version of this table, which is how a
story survives three rows. Two further pieces of evidence point away from it.
`simulations/dna_base_pairing.py:39-42` computes and prints its own kT in cm⁻¹
(215 at 310 K) six lines above the γ line, so transcribing 25 meV as 50 cm⁻¹
would contradict the author's own output; and the lines in question give
**ranges** (0.01 to 100, 10 to 100 cm⁻¹) where the meV table gives values, with
50 sitting mid-range. That file also separates sourced from unsourced explicitly
at `:198-199` ("estimated, NOT from literature"), which is not how blind copying
behaves. 50 and 250 cm⁻¹ are round numbers with obvious independent homes in
molecular spectroscopy.

So the finding is the attribution and the arithmetic, not the story: the set
names a source that does not carry it, and it is not the meV table in other
units. Where 10 to 100 cm⁻¹ actually comes from is open.

In the pre-repair baseline `ec7cc619fa075d82137698cedee27b742b7dd6fc`, the
0.01 appeared in ten places: four water readings
(`docs/water/HYDROGEN_BOND_QUBIT.md:141`,
`docs/water/PROTON_WATER_CHAIN.md:55,270`,
`ec7cc619fa075d82137698cedee27b742b7dd6fc:docs/carbon/README.md:386`) and six
DNA readings (`docs/water/HYDROGEN_BOND_QUBIT.md:178`,
`experiments/DNA_BASE_PAIRING.md:47,99,141,219`, and
`experiments/README.md:138`). The baseline's two bare sweep rows were
`docs/water/HYDROGEN_BOND_QUBIT.md:141` and
`experiments/DNA_BASE_PAIRING.md:99`.

**The rule this leaves.** γ is a rate, obtained either from a measured time as
1/(2T₂) or computed from a bath spectral density together with a coupling
strength. A bath energy scale on its own is neither.

That rule is already met on hardware, and has been for months: γ = 1/(2·T₂) in
µs⁻¹, structurally discriminated from the T₁-aware form (`GLOSSARY.md:208`),
implemented in `compute/RCPsiSquared.Core/Calibration/CalibrationChain.cs`,
gated by `simulations/t2_gamma_book_gate.py`, its failure modes logged twice in
`CAUGHT_ERRORS.md`. Nothing in the substrate folders meets it.

A second, unaudited question sits on the numerator side and is the same shape
as the Zundel finding: Hückel β is a single-particle π-orbital hopping, while
the framework's J is an XX+YY exchange between two-level sites. [The framework
lens](carbon/BENZENE_HUCKEL_FRAMEWORK_LENS.md#framework-vocabulary-translation)
records β → J as a Tier-2 structural identification for that model translation.
Whether it selects the right J for a material carbon system is not settled;
[carbon's conditional C4 and C6 working model](carbon/README.md#conditional-c4-and-c6-working-model)
records that the underlying material two-state degree of freedom is unassigned.

## A conditional ceiling in the selected-coordinate model

With ℏ = 0.6582 meV·ps and the repository's `D[Z]` convention,
`Q = 2·J·T₂/ℏ`. No proton-coordinate `T₂` or microscopic bath model is on
file, so this relation supplies no water-model lower endpoint.

One illustrative ceiling can be stated only after three additional
stipulations: a chosen proton coordinate, the ice-derived `J = 0.5 meV`
convention for that coordinate, and the use of the liquid-water H-bond lifetime
as a proxy upper bound for its unavailable `T₂`. [Physical
parameters](water/HYDROGEN_BOND_QUBIT.md#physical-parameters) records that
lifetime as 1 to 3 ps (Luzar & Chandler 1996). Taking `T₂ ≲ 3 ps` gives

> Q ≲ 2 · 0.5 meV · 3 ps / 0.6582 meV·ps ≈ **4.6**

This is an illustrative conditional ceiling, `Q ≲ 4.6`, not a measurement or a
bound for ordinary liquid water. The 3 ps value is the top of a mean bond
lifetime rather than a maximum or a proton-coordinate coherence time. The
`J = 0.5 meV` input is an ice-derived convention rather than a direct
identification of a measured splitting with the matrix element `J`; its mapping
and validity for liquid water are unverified. Without those stipulations, and
without a coordinate `T₂` or bath model, the repository assigns no water Q.

## A pre-repair water-folder mismatch

The pre-repair baseline `ec7cc619fa075d82137698cedee27b742b7dd6fc` placed the
following two readings side by side. At `docs/water/README.md:55`, embedding
condition 4 read:

> "**Decoherence ~ J** (proton tunneling rates and bath fluctuations on the
> same picosecond scale) → Q is in the framework's testable range."

At `docs/water/HYDROGEN_BOND_QUBIT.md:211`, it read:

> "For liquid water at 300K: J ~ 0.5 meV, γ ~ 25 meV. J/γ ~ 0.02.
> Classical regime. The palindrome exists but is overdamped."

Same substance, same question. The README does not estimate a Q, it asserts one
as an **embedding condition** ("under these conditions the F-chain inherits"),
so the two are an assumption and an estimate rather than two estimates. Read its
picosecond scale as T₂ ~ 1 ps and the same Lindblad book used throughout gives
Q = 2 · 0.5 meV · 1 ps / 0.6582 meV·ps ≈ 1.5; at the top of its own 1 to 3 ps
range it gives 4.6, the band's upper end. The parameter table said 0.02. That is
roughly two orders of magnitude, in one folder, from May 2026 until now, with
neither page citing the other. The README reasons in times, the parameter table
reasons from a bath energy, and nothing in the folder made them meet.

## The letter Q, and how many things wear it

`GLOSSARY.md:317-334` already lays this out and names more senses than the three
the substrate and neural material actually uses. The three below are the ones
that collide in that material; the glossary's own list adds the coherence
horizon Q*(N), the gap threshold Q*_gap(N), and the right-popcount index q of a
(p,q)-block.

**The dangerous one is not in the list at all**, because it is not a different
quantity but a second normalization of the same one: **q = Q/2**, the F89 octic
book (`GLOSSARY.md:326`, worked at `:328`, where q_EP ≈ 0.659 octic is Q = 1.318
carrier). A visibly different object gets noticed; a silently halved one does
not.

1. **Q = J/γ₀**, the scale. This document, `Q_REGIME_ANCHORS.md`,
   `Q_SCALE_THREE_BANDS.md`.
2. **Q_max = |Im λ| / |Re λ|**, a resonator quality factor, defined at
   `experiments/VEFFECT_CAVITY_MODES.md:189` as `J·μ_max/γ` and computed at
   `simulations/neural/neural_gamma_cavity.py:247`. Both numbers in the
   comparison it draws (C. elegans 0.1 against the qubit chain's 68 to 75) are
   this object, not object 1. That comparison was withdrawn on 2026-08-25: see
   `experiments/NEURAL_GAMMA_CAVITY.md`, Result 2b, which supplies the
   numerator and denominator provenance item (5)(b) of the arc
   `substrate_q_provenance` asked for. The arc stays open; (5)(b)'s third leg,
   the thermal-regime explanation, was deleted rather than answered. The 68 to
   75 is not a property of any qubit: on the N=5 chain
   `VEFFECT_CAVITY_MODES.md:170` gives 72.4 = μ_max·J/γ with
   μ_max = 2(1+cos 36°) = 3.618, evaluated at the sweep values J = 1 and
   γ = 0.05 that the source names twice (`:234`, `:239`). Since J/γ is object 1,
   **object 2 is object 1 times a graph invariant**: 3.618 × 20. It is not
   independent of object 1, and 72.4 is fixed by Q, not by γ. Attributing it to
   γ alone is the error this document exists to name, and an earlier draft of
   this line did exactly that.
3. **Q = the count of CΨ = ¼ crossings**, in
   `hypotheses/COMPLEXITY_THRESHOLD.md:33`.

## What this does not claim

- Not that ordinary liquid water has an assigned Q, or that it is quantum at
  room temperature, or that it is classical. The sole water-adjacent numerical
  statement is the illustrative conditional ceiling `Q ≲ 4.6`, which requires a
  chosen coordinate, the ice-derived `J = 0.5 meV` convention, and the H-bond
  lifetime proxy for the absent coordinate `T₂`; it has no lower endpoint.
- Not a value for T₂ of the proton coordinate in a hydrogen bond. That number
  is not in this repository and this document does not supply one.
- Not that the substrate Q values are arithmetically wrong. Given their inputs
  they are correct divisions. The inputs are the finding.
- Not a correction to `Q_SCALE_THREE_BANDS.md`,
  `THE_GENESIS_OF_AN_OSCILLATION.md`, or the tick reflection. They are right,
  they came first, and the scale result is theirs. What is new here is the
  provenance audit, the band, and the join.
- Not that one ratio always suffices. With J, γ_A and γ_B the unit still
  divides out, since L is homogeneous of degree 1 in all its rates jointly, but
  the shape then depends on J/γ_A and γ_B/γ_A separately rather than on a single
  Q (`experiments/OBSERVER_DEPENDENT_VISIBILITY.md`, its "Important Correction"
  section, which resolves the same document's earlier "depends on the absolute
  noise rate" wording into exactly these ratios; that wording has since been
  withdrawn from the prose there and from the linked script, though the stored
  output of the earlier run still carries it, and an earlier draft of this
  section had adopted it as written).

## Open

- **Missing inputs for a water-model Q.** A proton coordinate, its associated
  coupling, and a decoherence channel must all be selected. A
  proton-coordinate T₂ would constrain the last input, but cannot alone
  collapse a band. [The README's open follow-ups](water/README.md#open-follow-ups)
  scopes the needed γ_Z(T) estimate and pump-probe IR data.
- **The Zundel coupling.** 124 meV is unsourced and suspected to be a stretch
  fundamental; `experiments/DNA_BASE_PAIRING.md:96` carries the same cation at
  250 cm⁻¹ ≈ 31 meV. The historical 4.8 and 5.0 figures assign no current
  Zundel Q; neither candidate coupling is sourced. [Open item 1](water/PROTON_WIRE_CROSSING.md#open)
  retains the primary-source question.
- **Where 10 to 100 cm⁻¹ came from.** The cm⁻¹ set names a source that does not
  carry it. Whether the range has a real origin elsewhere, or is the meV numbers
  transcribed, is not settled here. Two of the repository's substrate Q values
  rest on it.
- **Hückel β.** 2.4 eV appears in five carbon documents with no citation
  anywhere in the repository, and whether a hopping integral belongs in the J
  slot is a separate and unasked question. The on-site α ≈ 11.4 eV appears in
   one document (`BENZENE_HUCKEL_FRAMEWORK_LENS.md:119-120`, which writes both with
   their signs, α ≈ −11.4 eV and β ≈ −2.4 eV; the audit above quotes magnitudes,
   which is harmless inside |J|/γ but not when α is quoted as a value). No separately
   sourced material bridge-`α` value is assigned in this repository.
- **Which J.** A substrate does not hand over one coupling. Superconducting
  hardware has a gate coupling in the MHz range and a static residual coupling
  orders of magnitude below it, and both are that substrate's J. Choosing
  between them is choosing the reading, not resolving an ambiguity. This is the
  substance of the opening claim.
- **The shipped defaults read at Q = 20.** `Propagate/Program.cs:1254` and
  `framework/chain_system.py:75` both default to J = 1.0 at γ₀ = 0.05.
  `Q_REGIME_ANCHORS.md:113` names this the pre-Q-band baseline and not a
  framework anchor. The canonical point is J = 0.075, Q = 1.5.

## Related

- [`experiments/Q_SCALE_THREE_BANDS.md`](../experiments/Q_SCALE_THREE_BANDS.md): the scale result and the three bands, with the rescaling identity.
- [`docs/Q_REGIME_ANCHORS.md`](Q_REGIME_ANCHORS.md): the ten named anchors on the Q axis.
- [`docs/THE_GENESIS_OF_AN_OSCILLATION.md`](THE_GENESIS_OF_AN_OSCILLATION.md): the exact factorisation L = γ₀·L₁(Q).
- [`docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md`](THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md): what cancelling the unit costs.
- [`reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md`](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md): one tick = 1/γ₀.
- [`docs/GLOSSARY.md`](GLOSSARY.md) §T₂ → γ: the measured route, and its two books.
- [`docs/GAMMA_TIME_DISTINCTION.md`](GAMMA_TIME_DISTINCTION.md): γ as the source of experienced time; its Part 3 verdict on τ = γt is sharpened by the scale reading here.
- [`docs/water/README.md`](water/README.md), [`docs/carbon/README.md`](carbon/README.md): the substrate folders whose Q values this document audits.
- The model dissipator this all refers to is built in [`simulations/water/hydrogen_bond_qubit.py`](../simulations/water/hydrogen_bond_qubit.py) `:132` as `√(γ/ℏ)·σ_z`, which does honour the rate conversion its own prose then bypasses.
