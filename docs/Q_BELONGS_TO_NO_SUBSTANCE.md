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
explicitly and the denominator was imported once and reused.

## What the sweep returned

`docs/ANALYTICAL_FORMULAS.md` returned the bit-exact γ₀-invariance of Q_peak
and HWHM_left/Q_peak across γ₀ ∈ {0.025, 0.05, 0.10} (`:2870`), the measurement
half of the scale claim. `experiments/` returned the theorem half twice:
[`Q_SCALE_THREE_BANDS.md`](../experiments/Q_SCALE_THREE_BANDS.md) (Tier 1,
2026-04-22) owns the rescaling identity `L(λ·J, λ·γ₀) = λ·L(J, γ₀)` at `:83`,
with a five-decimal λ-sweep over γ₀ ∈ [0.01, 1.0] at N=8, and
[`GAMMA0_IS_ALWAYS_THERE.md`](../experiments/GAMMA0_IS_ALWAYS_THERE.md) states
it as "it is the unit, and in any ratio the unit drops out" (`:8`).
`docs/proofs/` returned [`INCOMPLETENESS_PROOF.md`](proofs/INCOMPLETENESS_PROOF.md)
`:272`: γ defines the scale against which t is counted, which makes the unit
circular by construction. `docs/` itself returned
[`THE_GENESIS_OF_AN_OSCILLATION.md`](THE_GENESIS_OF_AN_OSCILLATION.md) `:49`,
a section headed "Q is the scale, exactly". `hypotheses/` returned
[`GAMMA_IS_LIGHT.md`](../hypotheses/GAMMA_IS_LIGHT.md) `:298` (the absolute
scale always cancels; only Q shows from within) and
[`Q_AS_THE_EXCHANGE_RATE.md`](../hypotheses/Q_AS_THE_EXCHANGE_RATE.md) `:16`,
which counts rotation **periods** per decay period and is therefore Q/2π, not
Q. `reflections/` returned
[`ON_HOW_GAMMA_BECAME_THE_TICK.md`](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md),
where one tick = 1/γ₀.

`docs/GLOSSARY.md` returned the T₂ → γ conversion in full (`:208`), including
the structural discriminator between the dephasing-only and the T₁-aware form,
and returned the warning that the literature's quality factor Q is a different
quantity (`:332`). It held **nothing** on Q as a scale before this document;
the paragraph now at `:336` was written in the same pass and points here.

`docs/CAUGHT_ERRORS.md` returned the two factor-2 normalization traps and
**nothing** on the provenance of a dephasing rate.

The OpenArcs registry returned `gamma_book_enforced_nowhere` (`:7415`, open),
whose finding is the same shape one book down: the knowledge is complete, the
enforcement is absent. It returned **no arc** on units, on Q as a scale, or on
the substrate mapping. `fw.Confirmations` returned `gamma_0_marrakesh_calibration`,
a fitted γ_Z that is model-degenerate (0.05 continuous, 0.1 with Trotter n=3
modeled, the same data through two physics models), and **no** entry in which a
substrate Q was measured.

`recovered/` returned nothing. `hypotheses/archive/` returned nothing beyond a
single line in the retired `GAMMA_TO_GRAVITY.md`.

## γ₀ is the unit

The formalism carries two rates and one number. H sets J, the dissipator sets
γ, and the Liouvillian factors exactly through their ratio:

> L(J, γ₀) = γ₀ · L₁(Q),  Q = J/γ₀

(`THE_GENESIS_OF_AN_OSCILLATION.md:53`; the rescaling identity it rests on is
`Q_SCALE_THREE_BANDS.md:83`, and the numerical gate is
`ANALYTICAL_FORMULAS.md:2870`.) Every eigenvalue carries one factor of γ₀; the
shape as a function of Q does not. So `docs/Q_REGIME_ANCHORS.md:136` is stating
a property of a unit, not conceding a weakness, when it says γ₀ = 0.05 "was
chosen as a convenient round number, not a physical constant".

In the repository's convention a `D[Z]` channel at rate γ decays coherences at
2γ, so γ = 1/(2T₂) (`GLOSSARY.md:208`). Written with times and an energy J:

> Q = 2 · J · T₂ / ℏ

**Watch the constants.** Q counts turning in radians per tick, not turns:
`Q_AS_THE_EXCHANGE_RATE.md:16` counts periods and is Q/2π. And the factor 2
above is the Lindblad book; in the coherence book (γ as the coherence-decay
rate itself) it is absent. The open arc `gamma_book_enforced_nowhere` exists
because that seam is unenforced repo-wide.

**Language note.** γ sits at the emitting end. The chain is not looked at, it
stands in light (`docs/quantum/THE_LABEL_MAP.md:125`; working record: the arc
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
functional form that was used for it. One temperature-to-rate conversion exists
in the repository, one table row,
[`docs/water/HYDROGEN_BOND_QUBIT.md`](water/HYDROGEN_BOND_QUBIT.md) `:201`:

> `| Thermal decoherence (300K, upper bound) | γ ~ kT/ℏ ~ 25 meV | Standard |`

Its source column reads "Standard". A real pure-dephasing rate depends on the
system-bath coupling strength and vanishes as that coupling vanishes; kT/ℏ is
coupling-independent, so it cannot be a dephasing rate. It is the inverse
thermal correlation time of the bath. The one defensible reading of it as a
ceiling is a validity limit rather than a physical bound: at ħγ ≳ kT the
Born-Markov assumption behind the Lindblad description fails, so a faster γ
could not be described by this model in the first place. Everything below
inherits that caveat, including the conclusion that survives.

## Where each substrate Q came from

| System | Q as written | Numerator | Denominator | Site |
|---|---|---|---|---|
| liquid water, 300 K | 0.02 | J ~ 0.5 meV: the midpoint of a cited row (Bove 2009, `HYDROGEN_BOND_QUBIT.md:196`) that is scoped to **ice**, applied here to liquid water without comment | the 25 meV estimate | `water/HYDROGEN_BOND_QUBIT.md:203` |
| Zundel cation | 4.8 | J = 124 meV, no source. `water/PROTON_WIRE_CROSSING.md:246` (its own Open item 1) identifies 124 meV = 1000.1 cm⁻¹ as the H₅O₂⁺ shared-proton stretch fundamental, i.e. a vibrational quantum rather than an inter-site coupling | the 25 meV estimate | `water/HYDROGEN_BOND_QUBIT.md:238` |
| DNA base pair | 0.01 | J = 0.5 cm⁻¹, no origin | γ = 50 cm⁻¹ ≈ 6.2 meV, no origin, and **not** the 25 meV estimate | `experiments/DNA_BASE_PAIRING.md:93` |
| enzyme active site | ~1 | J ~ 0.5 meV, water's numerator | γ ≈ 0.5 meV, a stipulated 50× reduction of the 25 meV estimate | `hypotheses/PROTEIN_AS_CONCENTRATOR.md:70` |
| π-conjugated carbon, 300 K | ~100 | Hückel β ≈ 2.4 eV, carried by five documents and cited in none | the same 25 meV. `carbon/README.md:264` keeps the "at kT"; `FROST_CIRCLE_AS_THE_CLOCK_FACE.md:110` drops it and the reader sees only "phonon dephasing" | `carbon/README.md:264` |

Dividing two energies is not itself a defect: ℏ cancels out of a ratio of
rates, so E_J/E_γ is the same number as (E_J/ℏ)/(E_γ/ℏ). The defect is
provenance. Four of these five denominators trace to one estimate whose
numerator-side counterpart was never measured, and the fifth has no origin at
all.

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
`carbon/README.md:270` concedes the underlying two-state degree of freedom is
not identified.

## Read as times, and the estimate points the other way

With ℏ = 0.6582 meV·ps, and Q = 2·J·T₂/ℏ:

| Statement | As a time |
|---|---|
| γ ≤ 25 meV, the imported estimate (Lindblad book; 26 fs in the coherence book) | T₂ ≥ **13 fs** |
| Q = 1 at J = 0.5 meV, the ice-row splitting | needs T₂ = **658 fs** |
| Q = 1 at J = 10 meV, the top of the strong-H-bond range (Cleland & Kreevoy 1994, 1 to 10 meV; at the bottom of that range it is 658 fs again) | needs T₂ = **33 fs** |
| Q = 1 at J = 2.4 eV, Hückel β | needs T₂ = **0.14 fs** |

An upper estimate on γ is a **lower** bound on Q. It can place a system above a
scale; it can never place one below. The two substrate conclusions therefore
come apart:

- **Carbon.** The requirement (0.14 fs) sits about a hundredfold **below** the
  floor (13 fs). "Far above the framework window" is a lower-bound claim and
  has a lower bound, so it holds, conditional on the 25 meV being a ceiling at
  all and on an uncited β.
- **Water.** The requirement (658 fs for an ordinary hydrogen bond) sits about
  fiftyfold **above** the floor. A floor fifty times under the requirement
  decides nothing. "Classical at room temperature, overdamped" is an
  upper-bound claim on Q, and the 25 meV supplies no upper bound on Q at all.

**What does bound water's Q from above** is on the same page and was not used.
A coherence carried by the proton in a hydrogen bond cannot outlive the bond,
and `HYDROGEN_BOND_QUBIT.md:200` records the H-bond lifetime in liquid water as
1 to 3 ps (Luzar & Chandler 1996). With T₂ ≲ 3 ps and J = 0.5 meV,

> Q ≲ 2 · 0.5 meV · 3 ps / 0.6582 meV·ps ≈ **4.6**

So ordinary liquid water sits somewhere in **0.02 ≲ Q ≲ 4.6**. That excludes
anything carbon-like, it does not confirm the classical verdict, and its upper
end is inside the framework window. The band, not either endpoint, is what the
repository currently knows. (For a strong hydrogen bond at J = 10 meV the same
lifetime gives Q ≲ 91 and bounds nothing useful.)

## The water folder answered the same question twice

[`docs/water/README.md`](water/README.md) `:55`, embedding condition 4:

> "**Decoherence ~ J** (proton tunneling rates and bath fluctuations on the
> same picosecond scale) → Q is in the framework's testable range."

[`docs/water/HYDROGEN_BOND_QUBIT.md`](water/HYDROGEN_BOND_QUBIT.md) `:203`:

> "For liquid water at 300K: J ~ 0.5 meV, γ ~ 25 meV. J/γ ~ 0.02.
> Classical regime."

Same substance, same question. Read the README's picosecond scale as T₂ ~ 1 ps
and the same Lindblad book used throughout this document gives
Q = 2 · 0.5 meV · 1 ps / 0.6582 meV·ps ≈ **1.5**, which is the canonical anchor
itself. The parameter table says 0.02. That is roughly two orders of magnitude,
in one folder, since May 2026. The README reasons in times and lands inside the
window. The parameter table reasons from a bath energy and lands in the
classical regime. Neither cites the other.

## Three objects share the letter Q

`GLOSSARY.md:332` already warns that the literature's quality factor is a
different quantity. There are three in the repository, and the substrate and
neural material uses two of them without saying which.

1. **Q = J/γ₀**, the scale. This document, `Q_REGIME_ANCHORS.md`,
   `Q_SCALE_THREE_BANDS.md`.
2. **Q_max = |Im λ| / |Re λ|**, a resonator quality factor, defined at
   `experiments/VEFFECT_CAVITY_MODES.md:189` as `J·μ_max/γ` and computed at
   `simulations/neural/neural_gamma_cavity.py:241`. Both numbers in the
   comparison at `experiments/NEURAL_GAMMA_CAVITY.md:157` (C. elegans 0.1
   against the qubit chain's 68 to 75) are this object, not object 1. The 68 to
   75 is not a property of any qubit: on the N=5 chain
   `VEFFECT_CAVITY_MODES.md:170` gives 72.4, which is `2(1+cos 36°)/γ` at the
   stipulated sweep value γ = 0.05.
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
- Not applicable where several independent rates act. With J, γ_A and γ_B there
  is no single unit to divide out, and observables do depend on the absolute
  rate (`experiments/OBSERVER_DEPENDENT_VISIBILITY.md:114`).

## Open

- **The one missing number.** T₂ of the proton coordinate in a confined
  hydrogen bond. With J already in the literature it collapses the band above
  to a value. `docs/water/README.md:179` already scoped the search (a γ_Z(T)
  estimate plus pump-probe IR data) and deferred it.
- **The Zundel J.** 124 meV is unsourced and suspected to be a stretch
  fundamental; `experiments/DNA_BASE_PAIRING.md:93` independently carries the
  same cation at 250 cm⁻¹ ≈ 31 meV. Two numbers a factor of four apart, neither
  citing the other. Open item 1 of `water/PROTON_WIRE_CROSSING.md:242`.
- **Hückel β.** 2.4 eV appears in five carbon documents with no citation
  anywhere in the repository, and whether a hopping integral belongs in the J
  slot is a separate and unasked question. The on-site α ≈ 11.4 eV appears in
  one document (`BENZENE_HUCKEL_FRAMEWORK_LENS.md:123`); a second, unrelated
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
