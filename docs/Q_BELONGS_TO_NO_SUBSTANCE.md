# Q Belongs to No Substance

*last refreshed 2026-08-22 (the change history lives in git)*

A substance does not hand over a Q. Q = J/γ₀ is defined only once a two-level
degree of freedom, a coupling, and a dephasing channel have been chosen, and a
real material offers several of each. "The Q of water" is therefore
under-determined by water: it is a property of a system in a condition and a
channel, not of a substance.

What γ₀ is, separately, is the unit. The Liouvillian factors exactly through
the ratio, so the absolute rate carries no physics and only Q does. That part
the repository already owned, proved and gated, in three lineages that were
never in one room. This document puts them in one room and adds what follows
for the substrate folders, where the choice of J and γ was never made
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

The formalism carries two rates and one number. H sets J, the dissipator sets
γ, and the Liouvillian factors into a unit times a function of their ratio:

> L(J, γ₀) = γ₀ · L₁(Q),  Q = J/γ₀

(`THE_GENESIS_OF_AN_OSCILLATION.md:53`; the rescaling identity it rests on is
`Q_SCALE_THREE_BANDS.md:83`, and the numerical gate is
`ANALYTICAL_FORMULAS.md:2870`.) At fixed Q every eigenvalue is proportional to
γ₀, and the shape as a function of Q does not move. Read at fixed J instead the
statement is false, because changing γ₀ then moves Q and reshapes L₁(Q). The
scope both sources carry, uniform Z-dephasing and an H homogeneous of degree 1
in J, travels with it. So `docs/Q_REGIME_ANCHORS.md:136` is stating
a property of a unit, not conceding a weakness, when it says γ₀ = 0.05 "was
chosen as a convenient round number, not a physical constant".

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

## Temperature enters only through the rate, and one rate was guessed

Pure Z-dephasing exchanges no energy: it leaves the populations in the Z basis
standing and removes only the off-diagonals. Its steady state is I/d, which is
the β = 0 Gibbs state, so the repository is right to call it an
infinite-temperature bath ([`KMS_DETAILED_BALANCE.md`](KMS_DETAILED_BALANCE.md)
`:250`) rather than a bath with no temperature. What follows is sharper than
"the formalism has no temperature": the channel **form** is T-independent, and
the **rate** γ is the only place temperature can enter.

So a substrate γ that depends on T is not the error. The error is the
functional form that was used for it. Two temperature-attributed rates exist in
the repository and neither carries a derivation. The one the substrate documents
divide by is a single table row,
[`docs/water/HYDROGEN_BOND_QUBIT.md`](water/HYDROGEN_BOND_QUBIT.md) `:209`:

> `| Thermal decoherence (300K, upper bound) | γ ~ kT/ℏ ~ 25 meV | Standard |`

Its source column reads "Standard". A real pure-dephasing rate depends on the
system-bath coupling strength and vanishes as that coupling vanishes; kT/ℏ is
coupling-independent, so it cannot be a dephasing rate. It is the inverse
thermal correlation time of the bath. The one defensible reading of it as a
ceiling is a validity limit rather than a physical bound: at ħγ ≳ kT the
Born-Markov assumption behind the Lindblad description fails, so a faster γ
could not be described by this model in the first place. Everything below
inherits that caveat, including the conclusion that survives.

The second temperature-attributed rate is not in a document at all. It is a
printed line in a runner, `simulations/dna_base_pairing.py:195`, giving γ_deph
as 10 to 100 cm⁻¹ for the "molecular environment at 310 K", with no functional
form behind the attribution and no source. It is audited two sections down,
because what it turns out to be is not a second estimate.

## Where each substrate Q came from

| System | Q as written | Numerator | Denominator | Site |
|---|---|---|---|---|
| liquid water, 300 K | 0.02 | J ~ 0.5 meV: a round value inside a cited row (Bove 2009, `HYDROGEN_BOND_QUBIT.md:204`) that is scoped to **ice**, applied here to liquid water without comment | the 25 meV estimate | `water/HYDROGEN_BOND_QUBIT.md:211` |
| Zundel cation | 4.8 | J = 124 meV, no source. `water/PROTON_WIRE_CROSSING.md:246` (its own Open item 1) identifies 124 meV = 1000.1 cm⁻¹ as the H₅O₂⁺ shared-proton stretch fundamental, i.e. a vibrational quantum rather than an inter-site coupling | the 25 meV estimate | `water/HYDROGEN_BOND_QUBIT.md:255` |
| DNA base pair | 0.01 | J = 0.5 cm⁻¹ | γ = 50 cm⁻¹ ≈ 6.2 meV, and **not** the 25 meV estimate. Both come from the cm⁻¹ set below | `experiments/DNA_BASE_PAIRING.md:96` (J), `:97` (γ) |
| liquid water again, 300 K | 0.01 | the same J = 0.5 cm⁻¹ | the same γ = 50 cm⁻¹ | `water/PROTON_WATER_CHAIN.md:270`, `simulations/water/proton_water_chain.py:280-285` |
| enzyme active site | ~1 | J ~ 0.5 meV, water's numerator | γ ≈ 0.5 meV, a stipulated 50× reduction of the 25 meV estimate | `hypotheses/PROTEIN_AS_CONCENTRATOR.md:71` |
| π-conjugated carbon, 300 K | ~100 | Hückel β ≈ 2.4 eV, carried by five documents and cited in none | the same 25 meV, now carrying its caveats at `carbon/README.md:263` and `FROST_CIRCLE_AS_THE_CLOCK_FACE.md:109` | `carbon/README.md:263` |

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
| K, inter-bond | 0.1 meV (`HYDROGEN_BOND_QUBIT.md:149`) | 20 cm⁻¹ | 0.81 cm⁻¹ | 24.8 |

**What the table shows, and what it does not.** Two consequences are solid. The
γ row and the Zundel row are displaced by the same factor, so the Zundel Q
survives the crossing: 250/50 = 5.0 against 124/25.85 = 4.8, using kT at 300 K
rather than the table's rounded 25. The ordinary-water numerator is displaced by
twice as much, and that extra factor of two is the entire difference between the
two Q values the repository carries for one substance; converted rather than
reused, the cm⁻¹ set gives 4.03/201.6 = **0.02**, the meV table's own number.

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

The 0.01 it produces stands in ten places, four as water's verdict
(`HYDROGEN_BOND_QUBIT.md:141`, `PROTON_WATER_CHAIN.md:55`, `:270`,
`carbon/README.md:386`) and six as DNA's (`HYDROGEN_BOND_QUBIT.md:178`,
`DNA_BASE_PAIRING.md:47`, `:99`, `:141`, `:219`, and the experiments index at
`experiments/README.md:138`). Eight of the ten now carry their scope; the two
bare sweep rows, `HYDROGEN_BOND_QUBIT.md:141` and `DNA_BASE_PAIRING.md:99`, are
covered by the note in their own section rather than in the row.

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
the framework's J is an XX+YY exchange between two-level sites. Whether β is
the right object for that slot at all is not settled anywhere, and
`carbon/README.md:281` concedes the underlying two-state degree of freedom is
not identified.

## Read as times, and the estimate points the other way

With ℏ = 0.6582 meV·ps, and Q = 2·J·T₂/ℏ:

| Statement | As a time |
|---|---|
| the imported estimate, read as a validity limit on the coherence rate: 1/T₂ ≤ kT/ℏ, at kT = 25.85 meV | T₂ ≥ **25.5 fs** |
| Q = 1 at J = 0.5 meV, the ice-row splitting | needs T₂ = **658 fs** |
| Q = 1 at J = 10 meV, the top of the strong-H-bond range (Cleland & Kreevoy 1994, 1 to 10 meV; at the bottom of that range, 329 fs) | needs T₂ = **33 fs** |
| Q = 1 at J = 2.4 eV, Hückel β | needs T₂ = **0.14 fs** |

**Where the floor comes from, and why it is 25.5 fs and not 13.** The estimate's
one defensible reading is a Markov-validity limit: the bath correlation time
ℏ/kT must be short against the time the model has the coherence surviving. In
the Lindblad book a `D[Z]` channel at rate γ decays coherences at 2γ, so the
condition binds 1/T₂ and not γ, and the floor is T₂ ≳ ℏ/kT = 25.5 fs at
kT = 25.85 meV. Setting γ ≤ kT/ℏ instead permits T₂ = 13 fs, which is **half**
the bath correlation time, i.e. exactly the regime where the description being
used is guaranteed not to apply. The looser floor is self-undermining and the
correct one is twice as high.

The doubling carries all the way to Q, and every number below is the corrected
one. Q ≥ 2·J·T₂/ℏ with T₂ ≥ ℏ/kT is simply **Q ≥ 2J/kT**, twice the bare quotient
J/γ the substrate documents wrote down. For water that is 2·0.5/25.85 = **0.04**,
not 0.02, and an earlier version of this document applied the doubling to carbon
while leaving water's endpoint at the discredited 13 fs floor. The bare quotients
in the provenance table are what those documents computed; the floors on Q are
twice them.

An upper estimate on γ is a **lower** bound on Q. It can place a system above a
scale; it can never place one below. The two substrate conclusions therefore
come apart:

- **Carbon.** The requirement (0.14 fs) sits **186× below** the floor
  (25.5 fs), i.e. Q ≳ 186. "Far above the framework window" is a lower-bound
  claim and has a lower bound, so it holds. But a validity limit cannot bound a
  physical rate, so what is supported is a disjunction and not a milder verdict:
  **either** Q ≳ 186 for π-conjugated carbon at 300 K, **or** the Lindblad
  description does not apply there, in which case the framework Q is not defined
  for it. Both branches also inherit an uncited β.
- **Water.** The requirement (658 fs for an ordinary hydrogen bond) sits
  **26× above** the floor, i.e. Q ≳ 0.04. A floor twenty-six times under the
  requirement decides nothing. "Classical at room temperature, overdamped" is an
  upper-bound claim on Q, and the estimate supplies no upper bound on Q at all.

**What does bound water's Q from above** is on the same page and was not used.
A coherence carried by the proton in a hydrogen bond cannot outlive the bond,
and `HYDROGEN_BOND_QUBIT.md:208` records the H-bond lifetime in liquid water as
1 to 3 ps (Luzar & Chandler 1996). With T₂ ≲ 3 ps and J = 0.5 meV,

> Q ≲ 2 · 0.5 meV · 3 ps / 0.6582 meV·ps ≈ **4.6**

Two caveats travel with that number and neither is small. The 3 ps is the top of
a **mean** bond lifetime, not a maximum, so the second digit is not decidable and
the bound is an order of magnitude. And the J is the ice splitting the table
above flags as applied to liquid water without comment: both endpoints of this
band inherit it, so the band is conditional on the very input this document
criticises. For a strong hydrogen bond at J = 10 meV the same lifetime gives
Q ≲ 91 and bounds nothing useful.

So ordinary liquid water sits somewhere in **0.04 ≲ Q ≲ 4.6**. That excludes
anything carbon-like and it does not confirm the classical verdict. What it does
is worse for deciding: the band **contains the whole framework-anchor range**
[0.2, 2.0] (`Q_REGIME_ANCHORS.md:113`, the anchor table at `:7-18`), with room left over
at both ends. Every anchor the framework names, onset, Balance, the F86 peak, is
inside the band, so nothing in it can be ruled in or out. The band, not either
endpoint, is what the repository currently knows.

## The water folder had answered the same question twice

[`docs/water/README.md`](water/README.md) `:55`, embedding condition 4:

> "**Decoherence ~ J** (proton tunneling rates and bath fluctuations on the
> same picosecond scale) → Q is in the framework's testable range."

[`docs/water/HYDROGEN_BOND_QUBIT.md`](water/HYDROGEN_BOND_QUBIT.md) `:211`:

> "For liquid water at 300K: J ~ 0.5 meV, γ ~ 25 meV. J/γ ~ 0.02.
> Classical regime. The palindrome exists but is overdamped."

That is the text as it stood until this audit; the second sentence has since
been replaced there, and that page now carries the band and links here, so the
two no longer disagree. What is worth keeping is how long they did.

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
   `simulations/neural/neural_gamma_cavity.py:241`. Both numbers in that
   comparison (C. elegans 0.1 against the qubit chain's 68 to 75) are this
   object, not object 1. The comparison itself was withdrawn on 2026-08-25:
   see `experiments/NEURAL_GAMMA_CAVITY.md`, Result 2b, which also closes the
   neural half of item (5)(b) of the arc `substrate_q_provenance`. The 68 to
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

- Not that water is quantum at room temperature, and not that it is classical.
  The band 0.02 ≲ Q ≲ 4.6 is what the sources support; both verdicts in the
  repository were stated outside it.
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
  Q (`experiments/OBSERVER_DEPENDENT_VISIBILITY.md:154-160`, which resolves its
  own earlier "depends on the absolute noise rate" wording at `:38` and `:114`
  into exactly these ratios; that wording is never withdrawn there, and an
  earlier draft of this section had adopted it as written).

## Open

- **The one missing number.** T₂ of the proton coordinate in a confined
  hydrogen bond. With J already in the literature it collapses the band above
  to a value. `docs/water/README.md:179` already scoped the search (a γ_Z(T)
  estimate plus pump-probe IR data) and deferred it.
- **The Zundel J.** 124 meV is unsourced and suspected to be a stretch
  fundamental; `experiments/DNA_BASE_PAIRING.md:96` carries the same cation at
  250 cm⁻¹ ≈ 31 meV. The factor of four between them is the same displacement the
  section above finds on the denominator, which is why the two Zundel Q values
  agree; it does not tell us which numerator is right, and neither is sourced.
  Open item 1 of `water/PROTON_WIRE_CROSSING.md:242`.
- **Where 10 to 100 cm⁻¹ came from.** The cm⁻¹ set names a source that does not
  carry it. Whether the range has a real origin elsewhere, or is the meV numbers
  transcribed, is not settled here. Two of the repository's substrate Q values
  rest on it.
- **Hückel β.** 2.4 eV appears in five carbon documents with no citation
  anywhere in the repository, and whether a hopping integral belongs in the J
  slot is a separate and unasked question. The on-site α ≈ 11.4 eV appears in
  one document (`BENZENE_HUCKEL_FRAMEWORK_LENS.md:123`, which writes both with
  their signs, α ≈ −11.4 eV and β ≈ −2.4 eV; the audit above quotes magnitudes,
  which is harmless inside |J|/γ but not when α is quoted as a value); a second,
  unrelated
  α ≈ 0.2 to 0.4 eV is the bridge coupling at
  `SINGLET_FISSION_AND_THE_TWO_CLOCKS.md:58`.
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
