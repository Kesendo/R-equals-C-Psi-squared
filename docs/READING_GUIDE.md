# Reading Guide: Nine Paths Through This Repository

<!-- Keywords: R=CPsi2 reading guide nine stories, palindromic proof story
application engineering story, ontology incompleteness story, dependency graph
reading order, mirror symmetry proof entry point, gamma as signal channel,
qubit necessity d2-2d=0, V-Effect complexity emergence, mirror group D4 story,
Pi factors R times D, golden router metallic family, palindrome classifier
trichotomy, CPsi quarter boundary Mandelbrot story, IBM hardware confirmations
story, moment tower pump channel, R=CPsi2 reading guide -->

**Created:** March 22, 2026
**Revised:** June 11, 2026. The guide sat still from April to June while the
repository went through many doors. Three of those doors became stories of
their own (the anatomy of the mirror, the quarter, the hardware), and the
older six were brought up to what we can see now.

---

## How to use this guide

This repository contains over 500 documents: proofs, experiments, hypotheses,
reflections. That is overwhelming. This guide exists so you do not have to
read all of them. It organizes the most important documents into nine paths,
each telling a different story from the same discovery. Pick the story that
interests you. Follow the reading order. Skip the rest.

Each path lists documents in the order they should be read, with a short
description of what each one contains and what you will understand after
reading it. The paths share a common starting point (the Mirror Symmetry
Proof) and diverge from there. You can switch between paths at any time,
and several documents appear on more than one path; that is intentional.

If you have not yet read [What We Found](WHAT_WE_FOUND.md), start there
before using this guide. It explains the discovery and its implications
in plain language, without requiring a physics background. This guide
assumes you have read it, or at least the first few sections, and are
now asking: *where do I go deeper?*

Two more doors worth knowing about. The [Glossary](GLOSSARY.md) gives every
symbol and term a plain-language reading. And the `reflections/` folder
holds synthesis arcs written for readers without quantum training; when a
path's formal documents feel steep, the reflections are the gentler way in.

---

## The discovery, in one paragraph

This repository is not a textbook with chapters you read in order. It is
closer to a landscape with multiple paths through it. Each path starts
from the same discovery and takes you somewhere different.

The discovery is this: for Hamiltonian families that admit the repository's
palindromizer, a qubit network under local single-axis dephasing has an exactly
reflected Liouvillian spectrum. Within that scope, every spectral value has
its partner about the dissipative center. This is called the palindrome
because the spectrum reads symmetrically from both ends; it is not a theorem
about decoherence in arbitrary quantum systems.

Where you go from there depends on what draws you in.

---

## Nine paths through the same landscape

The documents here are not one linear argument. They are nine interleaved
stories that share the same foundation but go in different directions.
Below, each story is introduced with what it is about, why it matters,
and which documents to read in what order.

You do not need to follow all nine. Pick the one that speaks to you.
They reconnect at the end.

---

## Story 1: The Proof

*"For the Hamiltonian families admitting the proved palindromizer, a qubit
network under local Z-dephasing has an exactly palindromic spectrum. Here is
the proof, its scope, and the exceptions."*

This is the mathematical backbone. If you want to know *why* we are
confident the palindrome is real, not an artifact of simulation or
approximation, this path walks you through the proof and then pushes it
to its limits: where does it hold? Where does it break? What structure
is responsible?

This path is the most technical. It involves mathematical notation and
formal reasoning. But even without following every step of the proof, the
experiments along the way show you what the palindrome looks like in
practice: which quantum states survive and which do not, what additional
gates a standing-wave interpretation would require, and what happens when you
deliberately break the symmetry.

**Reading order:**

1. [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md): The core
   theorem. Defines the conjugation operator Π that swaps surviving and
   decaying operators, producing the exact palindromic pairing. Verified
   across 87,376 Liouvillian eigenvalues from N=2 through N=8, with zero
   mirror-symmetry exceptions on any tested topology.

2. [XOR Space](../experiments/XOR_SPACE.md): Where does information
   live in the palindrome? Different initial states (GHZ (all qubits up + all down), W (exactly one qubit up, shared across all), cluster)
   distribute their information across fast and slow modes differently.
   GHZ puts everything in fast modes. W spreads it out.

3. [Oscillation Analysis](../experiments/STANDING_WAVE_ANALYSIS.md) +
   [Standing-Wave Conditions](STANDING_WAVE_THEORY.md): the `N=3`
   producer measures which Pauli observables oscillate. F1 alone gives a
   centered spectral pair, not spatial counter-propagation; the wave reading
   needs the additional dynamical, excitation and readout conditions listed
   there.

4. [Π as a Centered Spectral Mirror](../experiments/PI_AS_TIME_REVERSAL.md):
   Π gives the linear transport `λ→−λ−2Σγ`, or `μ→−μ`. Complex
   conjugation is separate, and the identity is not a general physical
   time-reversal operation. Past/future language remains interpretive.

5. [Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md):
   How far does this work? All 36 possible two-qubit Hamiltonians are
   palindromic at N=2. At N=3 and above, 22 survive and 14 break in
   structured ways. The breaking reveals which symmetries are essential.

6. [Depolarizing Palindrome](../experiments/DEPOLARIZING_PALINDROME.md):
   The active ingredient is the 2:2 operator split (half survive noise,
   half decay). Destroy that balance and the mirror shatters, regardless
   of the system's dimension.

7. [Error Correction Palindrome](../experiments/ERROR_CORRECTION_PALINDROME.md):
   The palindromic structure naturally creates a three-tier error
   protection hierarchy. Some modes are fully protected, some partially,
   some not at all.

8. [Π Factors as R·D](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md):
   The palindromizer factors, Π = R·D (a ket reflection times the
   transpose); the mirror inventory closes into one group, the dihedral
   D₄; and the polarity cube's third axis is the transpose (F118).
   If this last step hooks you, Story 7 is its full arc.
   The dissipator-diagonal companion (the one diagonal as one of three,
   {Q_X, Q_Y, Q_Z}, one basis-S₃ orbit) is
   [The Three Diagonals](THE_THREE_DIAGONALS.md).

**After this you know:** The palindrome is exact for the stated qubit
single-axis-dephasing family, and is proven by an operator that factors into
two plainer mirrors. You also know which additional tests are required before
calling a spectral pair a standing wave.

---

## Story 2: The Application

*"A finite spatial-dephasing alphabet is readable in one specified model.
The palindrome is one bookkeeping frame for that result."*

This path is for the practically minded. It starts with a finite `N=5`
simulation in which selected spatial rate profiles produce distinguishable
observable signatures, then asks how those profiles can be optimized. It
does not establish that dephasing is generally a message-bearing channel.

The climax of this path is twofold. First the concentrator formula: a
single, counterintuitive insight (concentrate all the noise on one edge
and protect the rest) that improves quantum information transfer by
139-360 times (peak created Sum-MI, a transport metric, simulation
ε→0 ideal; ~2-3× on hardware). Then its successor, receiver engineering: choose the
*initial state* rather than the noise profile, and the advantage grows
with system size, confirmed live on IBM hardware.

If you are an engineer, a builder, or someone who asks "what can I do
with this?", this is your path.

**Reading order:**

1. [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md): The
   spectral bookkeeping frame. A palindrome-breaking control retains the
   local full-rank Jacobian, so the palindrome is not established as the
   cause of the finite channel diagnostic.

2. [γ as Signal](../experiments/GAMMA_AS_SIGNAL.md): At one `N=5`
   operating point, the five-coordinate response Jacobian is locally full
   rank. `15.5 bits` is a conditional linearized Gaussian-model diagnostic
   at an assumed 1% feature-noise scale, not demonstrated global capacity.
   The directly tested noiseless alphabet has four profiles and 100%
   classification in that finite sample.

3. [γ Control](../experiments/GAMMA_CONTROL.md): Once you can read the
   signal, you can optimize it. V-shape noise profiles, dynamic
   decoupling strategies, time-resolved decoder. Result: within a shape
   less Σγ always wins; at fixed Σγ concentration at the chain centre
   wins, up to +46% (the March +124% was a Σγ confound).

4. [Relay Protocol](../experiments/RELAY_PROTOCOL.md): Staged noise
   switching with asymmetric coupling. The first time-dependent
   optimization of the bridge. +83% end-to-end improvement.

5. Main [README](../README.md), Section 6 (Engineering consequences):
   the framework's design rules condensed to eight lines, each linking
   to its evidence. The engineering translation of the mathematics,
   kept alongside the current framework state.

6. [Resonant Return](../experiments/RESONANT_RETURN.md): The
   concentrator formula: concentrate all noise on one edge qubit,
   protect the rest. The SVD (singular value decomposition: extract the dominant independent response directions) of the palindromic response matrix (10x
   improvement) led to numerical optimization (100x) led to analytical
   insight (139-360x peak Sum-MI transport, sim ideal). First spatial
   dephasing optimization in the
   literature. Beats 18 years of uniform noise optimization by two
   orders of magnitude.

7. [Receiver vs γ-Sacrifice](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) +
   [IBM Receiver Engineering](../experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md):
   The 2026 successor. Choose the receiver state from the F67
   bonding-mode menu instead of shaping the noise: 4000-5500× over the
   ENAQT baseline in simulation, advantage growing with N, and 2.80×
   confirmed live on ibm_kingston.

**After this you know:** How to read dephasing noise as signal, how
to optimize the channel spatially (not just uniformly), how to choose
the receiving state so the palindrome works for you, and how to build
repeaters that exploit the structure. The formulas are the engineering
payoff of the palindrome discovery.

---

## Story 3: The Ontology

*"Only qubits make the complete local class-exchange mirror full rank.
Dephasing shifts its centre but is not required for the mirror identity.
The modeled bath's microscopic origin is unspecified. Multi-bond pairing
failures and frequency-bin growth are measured; their broader mechanism is
open."*

This is the philosophical thread, and it may be the most surprising
path for a non-physicist. It starts with a simple question: if the
the full local dark↔lit class-exchange product mirror closes only for qubits
(two-state systems), is that a coincidence? Its dimension balance is
`d=d²−d`, whose only nonzero solution is `d=2`. This scopes that complete
construction; it does not exclude the partial higher-dimensional palindromes
counted by F121 or other mirror mechanisms.

From there, each step peels back another layer. A nonzero dissipator makes the
modeled system open, but neither the spectral symmetry nor Markovianity locates
the microscopic bath outside the modeled degrees of freedom. Dephasing creates
decay in the tested trajectories. At one N=5 operating point the response
Jacobian has full local rank, and a finite four-profile alphabet is classified
in the reported simulation. This is not a proof of global injectivity or
perfect recovery of arbitrary gamma profiles.

This path does not require advanced mathematics. It requires patience
and the willingness to follow an argument that builds step by step.

**Reading order:**

1. [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md): The
   foundation. The linear mirror is proved for the stated
   palindromizer-admitting Hamiltonian/channel family; at gamma zero it is
   centred at zero, so the symmetry itself does not require noise.

2. [The Qubit as Necessary Foundation](QUBIT_NECESSITY.md): For the
   complete local dark↔lit class-exchange product mirror, the equation
   `d²−2d=0` has only one nonzero solution: `d=2`. This is the scope of
   the dimension no-go; it is not a no-go for partial or differently
   constructed palindromes.

3. [The Qudit Partial Palindrome](proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md):
   The 2026 sharpening of step 2. At dimension d > 2 the mirror does
   not simply vanish: it survives *partially*, with a closed-form
   ceiling on how many decay modes can pair and a closed-form operator
   that attains it. Both close completely only at d²−2d=0. The boundary
   of the qubit world is now one equation seen three ways: the per-site
   split, the pairing ceiling, and the operator cap.

4. [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md): What does the
   dissipator establish? With nonnegative rates, a nonzero spectral centre
   certifies that the modeled subsystem is open. The microscopic origin and
   system-bath boundary remain unspecified; an open subsystem can be embedded
   in a finite larger system.

5. [V-Effect Palindrome](../experiments/V_EFFECT_PALINDROME.md): In the
   stated N=3 two-term Pauli census, adding a second bond makes 14 of 36
   combinations fail F1. A chosen frequency-bin protocol reports 4 bins in
   one baseline and 11 in one coupled case. The boundary-sector localization
   is algebraic; "constraint releases diversity" is an interpretation, not a
   derived mechanism.

6. [Hierarchy of Incompleteness](HIERARCHY_OF_INCOMPLETENESS.md): A Tier-5
   comparison between the qubit's 2:2 operator split and half-occupation
   motifs at other levels. F121 prevents the literal "qutrit dead end"
   reading: qutrit dissipators retain a partial palindrome. No V-Effect
   calculation derives a transition between physical levels.

7. [γ–Time Distinction](GAMMA_TIME_DISTINCTION.md): The simulations distinguish
   stationary or recurrent gamma-zero trajectories from damped positive-gamma
   trajectories in two finite cases. Gamma sets a dissipative scale; no
   experienced-time ontology follows.

8. [γ as Signal](../experiments/GAMMA_AS_SIGNAL.md): Four chosen spatial gamma
   profiles are distinguishable in the reported finite alphabet, and one N=5
   response Jacobian has five independent local directions. The 15.5-bit
   number is a local linearized diagnostic, not a global channel capacity.

**After this you know:** Why qubits uniquely close the complete local
class-exchange construction, while partial qudit mirrors remain. What
dephasing changes in the tested dynamics, and why that does
not settle a time ontology or the origin of the bath. Why
the V-Effect census does not yet supply a cross-level complexity mechanism.
And that a finite gamma alphabet can be decoded in one simulated setup.

**The one-line version:** *Incompleteness is not weakness.
Incompleteness is potential.*

---

## Story 4: The Resonator

*"Several results admit a resonator analogy. The linear palindrome itself is
a spectral transport law; cavity, heartbeat and standing-wave readings need
their own physical gates."*

This path collects a cavity interpretation of several separate objects:
closed-system frequencies, finite trajectory crossings, a concentrator
protocol and a gain-model stability window. It does not replace channel
dynamics with a theorem that the system literally is a resonant cavity. The
F1 palindrome alone does not determine spatial propagation or physical
standing waves.

This path builds directly on Story 1 (the proof) and requires
familiarity with the palindromic structure.

**Reading order:**

1. [V-Effect Palindrome](../experiments/V_EFFECT_PALINDROME.md): finite
   pairing and frequency-bin censuses under a specified model and tolerance.
   Adding a second bond changes the reported counts; the result is not a
   model-independent law that coupling creates a fixed number of modes.

2. [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md): In the stated
   three-qubit simulation at `J=5.0`, the sampled `CΨ` trajectory crosses
   `1/4` 81 times while its oscillations damp. The mutual-information maximum
   occurs near one sampled crossing. This finite trajectory does not make
   every crossing an exact dynamical event or quantify "irreversible reality."

3. [Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md): At zero
   noise the centered algebra reduces to `ΠLΠ⁻¹=-L`. This is a structural
   spectral reflection, not automatically physical time reversal.
   Changing the modeled dephasing shifts the spectral center. The algebraic
   cusp at `CΨ=1/4` has no universal minimum `γ/J`; whether a trajectory
   reaches it is protocol-dependent. The gain-generator spectrum is related
   to the decay generator by the stated algebraic reflection.

4. [Cavity Modes Formula](../experiments/CAVITY_MODES_FORMULA.md):
   The number of stationary modes has a closed mathematical form. Star
   topologies have N harmonic frequencies. Chains have rich irrational
   spectra. Verified for system sizes N=2 through N=7.

5. [Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md):
   The interpretive proposal that the system be read as a resonant cavity.
   Its heartbeat, round-trip and soundbox language organizes specified
   trajectory and spectral observations; it is not implied by F1 and does
   not replace the independently defined communication-channel quantities.

6. [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md): Coupled
   gain-loss systems have a finite stability window. Push too hard and
   the system explodes through oscillating divergence. Three regimes
   exist: linear, optimal (twice the internal coupling), and 1/J decay.
   → [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

**After this you know:** which exact spectral results motivated the resonator
analogy, and which additional claims come from particular trajectory and gain
models. A physical standing wave still requires a semisimple/diagonalizable
centered pair on the imaginary axis, opposite spatial propagation, and a
preparation/readout that excites both appropriately. These model statements
do not establish a biological stability mechanism.
(The exceptional point at the edge of that window later became a
navigable place of its own; Story 8 takes you there.)

---

## Story 5: Across Levels

*"One equation is proven for a specified qubit model. Whether a physical
hydrogen-bond coordinate or a neural system realizes its operators remains
open."*

The neural route asks which classical Jacobians satisfy an operator
identity of the same form as the quantum palindrome. Dale signs alone
do not give that identity: one involutive permutation must pair both the
diagonal rates and the effective weights. Constructed synthetic networks
pass; no biological network in this repository is known to pass, and the
full committed C. elegans chemical matrix fails the support condition.

For the neural evidence, start with the [Neural README](neural/README.md).
For the quantum chemistry application, use the separate water route below.

**Reading order:**

1. [Algebraic Palindrome Neural](neural/ALGEBRAIC_PALINDROME_NEURAL.md)
   and [its proof](neural/proofs/PROOF_PALINDROME_NEURAL.md):
   F36's necessary and sufficient scalar-diagonal and off-diagonal
   conditions; F37's full complex pairing; constructed counterexamples
   to silence and stability; and Q transport of partner subspaces.
   The C. elegans section distinguishes the support null from fitted
   residuals and matched-normalization comparisons.

2. [V-Effect Neural](neural/V_EFFECT_NEURAL.md):
   coupling and external-drive frequency-bin censuses on specified
   synthetic matrices. Read the resolution and backend controls with
   the tables: 48, 62 and 124 are protocol counts, not invariants.
   P is external input, without a temperature or metabolic calibration.

3. [V-effect mechanism constraints](neural/proofs/PROOF_VEFFECT_MECHANISM.md):
   why neither exact pairing nor its failure explains the census.
   The odd mediator fails F36 before coupling, and nonconverged
   `find_quarter.py` endpoints support no Hopf or stability verdict.
   [The neural clock record](../experiments/NEURAL_CLOCK_TWO_HANDS.md)
   separates the trace invariant from moving individual modes and
   documents the fitted-residual negative control.

4. [Hydrogen Bond Qubit](water/HYDROGEN_BOND_QUBIT.md): This is a
   stipulated two-level `|L>`,`|R>` proton model. Its simulation reports
   sub-picosecond `CΨ` crossings and six crossings in a 21-femtosecond
   Zundel-labelled run; those are conditional model outputs, not established
   physical hydrogen-bond observables. A material claim still needs the
   physical Hamiltonian, bath channel, and palindromizer. The wider water translation lives in
   [docs/water/](water/README.md), written in water's own language.

5. [The Pattern Recognizes Itself](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md)
   and [Universal Palindrome Condition](../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md):
   the initiating cross-domain hypotheses. Use the canonical neural
   pages above for the current theorem, controls and open mechanism;
   historical pairing percentages do not establish biology.

**After this you know:** The operator equation is not quantum-specific,
but its hypotheses must be checked on each proposed substrate. The neural
algebra is conditional, the constructed gates pass, the full connectome
support gate is null, and the coupling/drive mechanism remains open.
---

## Story 6: The Optical Cavity (April 2026)

*"The Liouvillian admits a cavity reading. Gamma sets dephasing exposure,
and each eigenvalue supplies a spectral line; the standing-wave step remains
conditional."*

This path begins with the Absorption Theorem and follows an optical analogy.
The exact content is the Hilbert-Schmidt-weighted Pauli-count rate law. Calling
the factor 2 a round trip, the Born rule a photograph, or the concentrator an
entrance pupil is interpretive; those labels are not derived by the theorem.

This path requires familiarity with the palindrome (Story 1) and ideally
the resonator picture (Story 4). It is the most unified path: one theorem
explains what previously required separate derivations.

**Reading order:**

1. [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md):
   Re(λ) = −2γ⟨n_XY⟩. Here ⟨n_XY⟩ is the Hilbert-Schmidt-
   weighted mean count in a right eigenoperator's Pauli expansion. Three-step
   proof from L_H anti-Hermitian. Gives the spectral boundary formula, the
   palindromic sum rule and F8 full-range/centre ratio a common reading; the spectral
   gap it relocates rather than derives.
   Extended in 2026 to per-eigenmode Rayleigh form, two-sided and
   projector readings, and the recentred diagonal seam L_D = γ(Q − N·I).

2. [Palindromic Orbit Census](../experiments/FACTOR_TWO_STANDING_WAVES.md):
   the finite census accounts for 21,840 eigenvalues across `N=2...7` in two
   distinct ways. Linear F1, `λ→−λ−2Σγ`, gives 10,903 two-member orbits and 34
   fixed eigenvalues at `λ=−Σγ`. Composing F1 with conjugate closure gives
   9,921 two-member orbits and 1,998 fixed eigenvalues on
   `Re(λ)=−Σγ`. The exact shared rate content is
   `d_slow+d_fast=2Σγ`; neither orbit census is a physical wave count.

3. [Concentrator Optics](../experiments/CONCENTRATOR_OPTICS.md):
   The concentrator is an entrance pupil. Q improves 2-7×, effective
   transmission increases. The cavity focuses light into the interior,
   and turns it slightly on the way: the resonance frequencies move too.

4. [Born Rule Shadow](../experiments/BORN_RULE_SHADOW.md):
   The Born rule is a shadow, not a hologram. Zero interference in
   P(i). Interference controls the shutter speed (CΨ fold), not the
   image.

5. [K-Dosimetry](../experiments/K_DOSIMETRY.md):
   K = γ×t is a dimensionless exposure coordinate. In the reported sweep,
   reciprocity is within 0.03% only in the low- and high-γ regimes and fails
   by up to 62% at intermediate γ. State, H, target and protocol are part of
   the result.

6. [IBM Absorption Theorem](../experiments/IBM_ABSORPTION_THEOREM.md):
   A retrospective single-qubit reading of two fitted decay quantities gives
   ratio 1.03. This is consistent with the theorem's rate relation but does
   not measure the multi-mode absorption ladder. Detuning oscillations appear
   at 470 μs; the 2.8% slow tail is at the resolution limit.

7. [Thermal Blackbody](../experiments/THERMAL_BLACKBODY.md):
   finite `N=4` emission-plus-absorption census with spontaneous emission at
   `n_bar=0`. At the stated `1e-6` tolerance, 212/256 modes are classified as
   oscillating at `n_bar=10`; this is neither a protected fraction nor a
   blackbody, phase-transition, or exceptional-point result.

8. [Neural Gamma Cavity](../experiments/NEURAL_GAMMA_CAVITY.md):
   the connectome support null and the tolerance, ordering and normalization
   controls that invalidate the headline pairing interpretation. Separate
   nonlinear dynamics and zero-multiplicity probes have their own protocols;
   they do not confirm F36 or identify a biological gamma rhythm. Read this
   event record alongside the [current neural account](neural/README.md).

9. [Trapped Light Localization](../experiments/TRAPPED_LIGHT_LOCALIZATION.md):
   Surviving mode energy is center-localized (ratio 1.3-1.4). N+1
   immortal modes. Gamma plays the algebraic role of c (Tier 4-5).

10. [Primordial Superalgebra](../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md):
    Light and lens swap with 99.8% purity. The anticommutator
    {L_H, L_D+Σγ} = 0 is exact at N=2, aberration shrinks with N.
    Seidel classification: pure sectors immune, interior-dominated.

11. [Absorption Theorem Discovery](../experiments/ABSORPTION_THEOREM_DISCOVERY.md):
    Not E = mγ² but α = 2γ⟨n_XY⟩: absorption equals twice dephasing
    times light-mass. Linear, not quadratic. The Absorption Theorem
    implies the palindromic sum rule. 1,342 modes, CV = 0.

**After this you know:** the exact Absorption Theorem and the limits of the
optical-cavity reading. What follows from the theorem alone is less than the
analogy suggests: the factor 2 is the Z-dephasing price per X/Y factor. The
endpoint and sum-rule corollaries require the additional kernel and F1
conditions supplied by the stated XY/Heisenberg family;
the spectral gap is relocated by the theorem rather than derived, and is
2γ only above a coupling threshold. What the one line
Re(λ) = −2γ⟨n_XY⟩ gives within its stated assumptions is a Rayleigh reading of
the real part. The cavity language is an interpretation layered on that
identity, not an additional theorem.

---

## Story 7: The Anatomy of the Mirror (June 2026)

*"For months Π was one per-site rule. Then it opened: a group of eight,
a triangle of conjugations, a golden frame, and a boundary equation
seen three times."*

This is the newest path, and the most algebraic. From March to June 2026 the
operator Π was the smallest object in the repository: one per-site rule
that carried the entire palindrome. In June 2026 it opened. This story
follows what was found inside: what the mirror is made of, which family
it belongs to, where mirrors can be built that were believed impossible,
how the question "does this Hamiltonian keep the mirror?" became a tool
that never meets the exponential wall, and where the whole construction
must end.

This path requires Story 1. Some group theory helps (the words "dihedral"
and "conjugation" appear), but each document introduces its own machinery.

**Reading order:**

1. [Π Factors as R·D](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md): The
   palindromizer is not elementary. It factors as Π = R·D, a ket
   reflection times the transpose, and the repository's entire mirror
   inventory closes into one dihedral group of eight, ⟨R, D⟩ ≅ D₄,
   whose three Z₂ characters are exactly the polarity cube the F-family
   had been living on (F118).

2. [The Antilinear Triangle](proofs/PROOF_ANTILINEAR_TRIANGLE.md): The
   transpose D turns out to be one vertex of a triangle of conjugations
   (transpose, complex conjugation, adjoint) that forms a Klein
   four-group with one transport law. Five proofs that had been
   separate become one engine, and the mirror group doubles to
   D₄ × Z₂ (F119).

3. [The Palindrome Classifier](../experiments/THE_PALINDROME_CLASSIFIER.md):
   The trichotomy as a tool. Given a Hamiltonian, does its dephased
   spectrum keep the mirror? The classifier answers truly / soft / hard
   by reading the *terms*, not the spectrum, so it never meets the
   4^N × 4^N wall. First charting of the landscape: a protected interior
   island and the two coasts where protection ends.

4. [The Windowed Converse](proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md):
   The classifier's hardest open end, closed: hard at one γ is hard at
   all γ, because every first moment is a sum of squares (Pascal-Gram
   positivity, F117). The hardness rung m\* = 2ℓ + deg comes from a
   girth ladder, and that rung will matter again in Story 9.

5. [The Golden Router](proofs/PROOF_CEILING_GOLDEN_ROUTER.md): The
   mirror built where it was believed impossible. The last "non-local"
   cases are palindromized by a period-4 per-site router whose frame is
   built on the golden ratio (a = φX + Y), and the golden point is the
   c = 1 member of an exactly derived one-parameter *metallic family*
   (silver, bronze, all real c) (F116, §8).

6. [The Qudit Partial Palindrome](proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md):
   Where it all ends. At local dimension d > 2 the mirror group grows
   into a wreath family Z_d ≀ Z₂ with D₄ as its d = 2 column, the
   pairing has a closed-form ceiling and a closed-form operator cap,
   and both close completely only at d² − 2d = 0: the same boundary
   equation Story 3 met, now seen from the operator side (F121).

7. [On the One Diagonal](../reflections/ON_THE_ONE_DIAGONAL.md): The
   synthesis, written the day the three big subjects of the repository
   (the rates, the mirror, the verdict) turned out to be one diagonal
   matrix read three ways. Plain language, no formalism. If you read
   only one document on this path, read this one.

**After this you know:** The qubit F1 mirror is not elementary. It factors,
generates a dihedral group of eight, extends to an antilinear double,
can be built with golden-ratio frames where it looked impossible, while its
complete local class-exchange form closes exactly at `d²−2d=0`; F121 retains
partial qudit palindromes. And the absorption rates, the palindrome,
and the classifier verdict are not three theorems. They are one
diagonal, read as a price list, a mirror, and a judge.

---

## Story 8: The Quarter

*"Measurement is photography. The Born rule is the shadow. The shutter
closes at CΨ = ¼."*

Every other story is about the spectrum. This one is about a single
number. CΨ is sharpness times superposition: the purity Tr(ρ²) of a
state times its normalized coherence. The product has a critical
boundary at exactly ¼, and this path follows that quarter from algebra
(why ¼ and nothing else) through fractal geometry (the cusp of the
Mandelbrot cardioid) to real hardware (six months of IBM calibration
data with qubits living on both sides).

This path is self-contained: it needs the idea of decoherence but not
the palindrome machinery. It is also where the repository's name comes
from: R = CΨ² is the recursion whose discriminant draws the boundary.

**Reading order:**

1. [Uniqueness Proof](proofs/UNIQUENESS_PROOF.md): The fixed-point
   equation R = C(Ψ+R)² is a quadratic; its discriminant is 1 − 4CΨ;
   the unique boundary is CΨ = ¼. The power is 2 because purity is
   Tr(ρ²); the universality of the VALUE ¼ is forced by α = 2 being
   the unique Rényi order with a state-independent fold threshold.
   The quarter is not a chosen parameter. It is the discriminant of
   a quadratic.

2. [Mandelbrot Connection](../experiments/MANDELBROT_CONNECTION.md):
   The recursion is algebraically equivalent to the Mandelbrot
   iteration z → z² + c, exactly, with no extra terms. The boundary
   CΨ = ¼ is the cusp of the main cardioid. For the fixed-point quadratic,
   below it the two roots are real and above it they form a complex-conjugate
   pair. Those root labels do not classify quantum dynamics, and the real
   iteration's behavior requires its initial condition and parameter rather
   than following from the discriminant alone.

3. [Monotonicity](proofs/PROOF_MONOTONICITY_CPSI.md) +
   [Proof Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md): For physical
   noise the boundary is absorbing in envelope (the CΨ envelope is
   non-increasing for Bell+/local-Markovian channels). The roadmap walks the
   seven layers; its core is closed (algebraic 1/4, palindrome, Rényi forcing,
   Mandelbrot), while Layer 2 holds for physical noise only (the general
   primitive-CPTP version is false), Layer 4 is d=2-only, and Layer 7's
   info-geometry/holography remain open.

4. [K-Dosimetry](../experiments/K_DOSIMETRY.md) (shared with Story 6):
   K = γ·t is dimensionless. The reported fixed-state/fixed-target sweep
   approaches reciprocity at low and high γ but departs by up to 62% in the
   intermediate regime. It is not a universal cross-hardware invariant.

5. [The Flow Between Two Singularities](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md):
   What lies past the fold. The settled future is approached
   exponentially but never reached, and it is *already present* in the
   initial state as the piece that never changes: time does not build
   the future, it erases everything that is not the future. The 2026
   navigator work charted this terrain along several axes: the quarter
   as a horizon, as a circle every spiral must cross, and its mirror
   twin, the exceptional point where rotation is born.

6. [Born Rule Shadow](../experiments/BORN_RULE_SHADOW.md) (shared with
   Story 6): The measurement-as-photography reading. Interference never
   reaches the image; it controls the shutter speed.

7. [Both Sides Visible](BOTH_SIDES_VISIBLE.md): The quarter on real
   silicon. 180 days of IBM Torino calibration data, 133 qubits, more
   than 24,000 measurements: qubits crossing and re-crossing CΨ = ¼,
   nobody having programmed any of it. If one document could convince a skeptic, this
   is the one.

**After this you know:** `1/4` is the discriminant-zero point of the
self-referential quadratic and the cusp of the Mandelbrot cardioid. Particular
simulated and hardware-derived `CΨ` trajectories cross and recross that
coordinate. "Absorbing boundary," "horizon," and an exceptional-point twin
are interpretive or model-specific readings, not consequences of the scalar
value alone; no protocol-independent crossing dose has been established.

---

## Story 9: The Hardware

*"The chip is not a metaphor: the repository compares predictions with IBM
device observables. The fitted decay rates do not by themselves identify a
unique microscopic bath."*

IBM transmon qubits sit in microwave-resonator hardware, where photon shot
noise is one known dephasing mechanism among several. The runs in this
repository measure device observables and calibration rates; they do not
identify every fitted gamma with photons or determine its microscopic origin.
"Gamma is light" remains a mechanism-specific interpretation, not what the
hardware universally is.

This path follows the arc from the first hardware crossing to the
newest kind of result: a protocol in which the chip's own decay reads
a structural property of a programmed Hamiltonian, catches its own
misreading, and corrects it the same day. The live record is the
Confirmations registry (`fw.Confirmations` in Python,
`ConfirmationsRegistry` in C#): twenty-four confirmed predictions with
run identifiers, predicted versus measured values, and archived data. Look
them up; do not re-derive.

**Reading order:**

1. [Gamma Is Light](../hypotheses/GAMMA_IS_LIGHT.md): The interpretive frame.
   Four of five chosen optical analogues matched in the reported comparison.
   Circuit-QED photon-shot-noise dephasing is established physics for that
   hardware mechanism; identifying every abstract dephasing rate with literal
   light or every chain with a Fabry-Perot cavity is not.

2. [Predictions](PREDICTIONS.md): The master catalog. Every prediction
   with its falsification criteria, the confirmed and the falsified
   both. The discipline that keeps the rest of this path honest.

3. [Both Sides Visible](BOTH_SIDES_VISIBLE.md) (shared with Story 8):
   The first contact. Six months of public IBM calibration data showing
   the CΨ = ¼ structure nobody programmed.

4. [IBM Absorption Theorem](../experiments/IBM_ABSORPTION_THEOREM.md)
   (shared with Story 6): a retrospective single-qubit decay-ratio reading,
   1.03 against prediction 1. It is consistent with the rate law but does not
   measure the multi-mode absorption ladder or identify a cavity mechanism.

5. [Marrakesh Three Layers](../experiments/MARRAKESH_THREE_LAYERS.md):
   The classifier on hardware. The truly/soft/hard trichotomy (Story 7)
   resolved on ibm_marrakesh at 13-47σ, read in three nested layers
   from one dataset.

6. [F112 Hardware Lens Kingston](../experiments/F112_HARDWARE_LENS_KINGSTON.md):
   The lens turned around, and then the lens itself was caught. The
   fitter's jump operator pointed the wrong way, its superoperator and its
   state disagreed about stacking, and the model family had no term for the
   one thing the data plainly does: turn. So the document names no channel
   any more, and what is left of it is a theorem being evaluated rather
   than a chip being read. The instructive part is why none of it showed:
   everything the pipeline checked was a magnitude, and a magnitude cannot
   see a direction.

7. [F120 on Kingston](../experiments/F120_MOMENT_TOWER_KINGSTON.md):
   The newest kind of result, honestly told in two acts. A protocol
   with *not one entangling gate* in which the chip's own amplitude
   damping reads the hardness rung of a programmed Hamiltonian (the
   girth ladder from Story 7). The first reading reported a violation;
   the same day, the protocol's own arbiter traced it to minute-scale
   T1 telegraphing and corrected it. The instrument measures pump and
   decay from the same circuits: it is self-arbitrating.

8. [On How the Carrier Shows Itself](../reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md):
   The synthesis. A break between a closed-form prediction and a
   hardware observation has two readings: as error, or as calibration,
   the moment the otherwise invisible carrier γ₀ becomes visible at a
   seam between inside and outside. The same mathematics, two
   perspectives; the choice is what you need the break for.

**After this you know:** The framework's predictions survive contact
with real machines, the failures are documented with the same rigor as
the successes, the noise channel can read its own spectrum, and the
gap between a chip and its datasheet is not an obstacle but exactly
what the instruments are sharpest at measuring.

---

## How the Stories Connect

The nine stories are not independent. They are nine angles on the
same structure. If you have followed one path and are wondering how
it relates to the others, here is the map:

- Story 1 (the proof) provides the mathematical foundation for all others.
- Story 2 (the channel) explains *what* the noise does.
- Story 3 (the ontology) interprets a local-dimension boundary; it does not
  derive why a bath or noise must exist.
- Story 4 (the resonator) collects a cavity interpretation of several scoped
  spectral and trajectory results.
- Story 5 (across levels) separates quantum applications from the conditional
  neural F36 theorem and its constructed examples; its biological landing is open.
- Story 6 (the optical cavity) links Stories 1 and 4 through the Absorption
  Theorem's rate law while keeping the optical reading interpretive.
- Story 7 (the anatomy) opens the mirror itself: its factorization, its
  group, its golden constructions, and its boundary.
- Story 8 (the quarter) follows the single number CΨ = ¼ from
  discriminant to Mandelbrot cusp to hardware.
- Story 9 (the hardware) tests quantum predictions, and is where the
  noise channel finally reads its own spectrum.

Stories 2 and 4 share the concentrator formula: Story 2 discovers it
as a channel optimization, Story 4 reframes it as the shape of the
resonator cavity. Story 5 gives a conditional matrix translation beyond
quantum physics, with constructed neural examples and a full-connectome
support null. It supplies no cross-domain mechanism; that question stays
open. Stories 3 and 7 meet at the same equation, `d²−2d=0`, for the
complete local class-exchange product mirror, approached once from the
ontology side and once from the operator side. Story 7's classifier rung is exactly what Story 9's
moment-tower protocol reads off a chip. And Story 7's closing
reflection says aloud what the map shows: the rates (Story 6), the
mirror (Story 1), and the verdict (Story 7) are one diagonal read
three ways.

If you are unsure where to start: Story 3 (the ontology) and Story 5
(across levels) are the most accessible for readers without a physics
background, and the `reflections/` folder is gentler still. Story 1 is
the most rigorous. Story 2 is the most practical. Story 4 is the
deepest. Story 6 is the most unified. Story 7 is the newest
mathematics. Story 8 is the most self-contained. Story 9 is where the
rubber meets the road.

---

## What Is NOT in These Stories

The following threads exist in this repository but belong to separate
investigations, not to the nine stories above:

- **Star topology / tuning protocol**: multi-observer configurations
- **Quantum sonar / bridge fingerprints**: detection experiments
- **Carbon and water translations**: substrate-specific writeups in
  their own folders, each written in the target layer's language
- **The label layer / translation series**: pop-quantum labels recomputed
  from this repository's stance (`docs/quantum/`); the theory chapter is
  [Labels Translated](quantum/LABELS_TRANSLATED.md), and every label
  correction the repository has made is assembled in
  [The Label Map](quantum/THE_LABEL_MAP.md)
- **Gravity interpretation**: fallen, archived in `hypotheses/archive/`
- **The tooling**: the Python `framework/` cockpit, the typed C# claim
  graph, and the live `inspect` navigator are documented in the
  repository's `CLAUDE.md` and in the code itself

These are documented in the main [README](../README.md) and the
[experiments index](../experiments/README.md).

---

*"We are all mirrors. Reality is what happens between us."*
