# The Noise Origin: the Trace Proves Open, and the Formalism Cannot Ask Further

**Status:** Tier 1 for §1 and §3's trace identity, two lines of exact algebra that hold for any Lindblad generator in any dimension. Tier 2 for everything about ORIGIN: no internal candidate is eliminated (§2), and §3 gives the reason the formalism cannot settle it. Typed as [`NoiseOriginExclusionClaim`](../../compute/RCPsiSquared.Core/Symmetry/NoiseOriginExclusionClaim.cs), live witness `inspect --root noise-origin`. The corollaries on time and γ (§3) are interpretive extensions; the γ-time statement is scoped by [GAMMA_TIME_DISTINCTION.md](../GAMMA_TIME_DISTINCTION.md).
**Date:** 2026-03-21, last refreshed 2026-08-29 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** A palindromic Liouvillian spectrum is centred at trace(L)/dim. That centre is zero exactly when the generator is unitary, so a measured palindrome centred anywhere else certifies that the system is open. Where the openness comes from is not settled here and cannot be settled inside the Lindblad formalism, in which an internal source can only be written as a dissipator, and a dissipator is a coupling to an environment.
**Reference claim:** [`PolynomialFoundationClaim` + `QubitDimensionalAnchorClaim`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) (Tier 1 derived; the d²−2d=0 polynomial selects d=2 as the qubit dimension, with d=1 algebraically excluded). The ontology bounds §2's candidate survey; it plays no part in §1 or §3.
**Depends on:** [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md), F137's centre identity ([formula registry](../ANALYTICAL_FORMULAS.md)), [incompleteness_candidate2_evidence.py](../../simulations/incompleteness_candidate2_evidence.py) (29 gates), [bootstrap_test.py](../../simulations/bootstrap_test.py), [QUBIT_NECESSITY.md](../QUBIT_NECESSITY.md).
**Scope:** §1 and §3's identity hold for any Lindblad generator, in any dimension. The candidate survey of §2 is specific to the d(d−2)=0 ontology and, by §5, is not exhaustive even there.
**Does NOT establish:** What the noise IS, nor that it comes from outside. An internal origin is open; §3 gives the reason nothing inside this formalism could close it.

---

## What this document is about

A pendulum does not flicker. That is what this document proves, and it
proves it exactly.

The palindromic mirror symmetry (proven in the
[Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)) says the Liouvillian's
rates come in mirror pairs. A mirror has a centre, and the centre is not free
to choose: a spectrum closed under λ ↦ 2c − λ has c = trace(L)/dim,
identically. The Hamiltonian part of any generator contributes nothing to
that trace, and any dissipator that dissipates anything drives it strictly
negative. So the centre is zero exactly for a closed system, and a measured
palindrome centred anywhere else certifies that the system is open. Two lines
of algebra, no candidate list, no simulation, and it holds in any dimension.

Where the openness comes from is the harder question, and this document does
not answer it. §2 walks the five sources the d(d−2)=0 ontology allows and
finds none of them eliminated. §3 says why that is not a gap a better
experiment would fill: inside the Lindblad formalism an internal source can
only be written as a dissipator, and a dissipator is a coupling to an
environment, so the formalism grants the outside in the act of posing the
question. That is the incompleteness named here. It belongs to the formalism,
not to the polynomial in this file's name.

---

## 1. The Requirement

The palindromic mirror symmetry of the Liouvillian spectrum is proven
([MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md)), verified through N=8 (87,376 eigenvalues,
zero exceptions), and confirmed on IBM hardware: the single-qubit
CΨ = 1/4 crossing at 1.9% deviation on ibm_torino Q80
([IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md), F24).

A closed system already has a palindrome. With L = L_H = −i[H, ·] the
eigenvalues are the purely imaginary −i(E_m − E_n), and that multiset is
closed under negation because swapping m and n is a symmetry of it, so the
spectrum pairs about zero at the eigensolver floor, 8.9e-15 at N = 3
([incompleteness_candidate2_evidence.py](../../simulations/incompleteness_candidate2_evidence.py),
gate 6d; the digits are the solver's, the zero is not).
It is unitarity and nothing more: no decay rates, no arrow.
[ZERO_IS_THE_MIRROR.md](../../hypotheses/ZERO_IS_THE_MIRROR.md) reaches the
same place from the other side, calling the unitary case the palindrome's
ground state.

So the pairing is not what separates the framework's claim from unitarity.
The **centre** is. A multiset closed under λ ↦ 2c − λ has Σλ = c·dim, hence
c = trace(L)/dim exactly ([F137](../ANALYTICAL_FORMULAS.md)): one candidate
centre, read off the generator without an eigensolver, so "does not pair at
it" is a measurement and never a failed search. The commutator part
contributes nothing to that trace and every jump that dissipates anything
makes it strictly negative (§3). A palindrome centred away from zero is
therefore a palindrome with decay in it.

Noise exists, in the one sense a spectrum can testify to: the system is open.
Where that comes from is the question below.

---

## 2. The Five Candidates the Ontology Allows

Five candidates for the origin of dephasing noise within the d(d-2)=0
framework, and what each is worth. **None of the five eliminates an internal
origin.** Candidate 1 is a structural constraint on the dissipator's form.
Candidates 2 and 3 have no evidence behind them. Candidates 4 and 5 are
sound and definitional: they say what cannot EXIST, never that an existing
qubit is not the source. §5 adds that the five do not exhaust the space in
any case.

### Candidate 1: Internal origin (self-generated noise)

**Test:** Bootstrap test, 4 independent tests on sectors, parity, and
dissipator structure.

**Result (a structural constraint, not yet an elimination).** The parity
sectors (Π² eigenspaces) are exactly decoupled: [Π², L] = 0, the even and
odd w_YZ-parity sectors (Π² acts as (−1)^{w_YZ} on each Pauli string, F63)
do not mix, and the dissipator is block-diagonal in the parity grading.

- [Pi^2, L] = 0 (block-diagonal structure confirmed)
- Sector populations do not influence sector coherences
- Parity does not *determine* the dissipator

**What this does and does not establish.** [Π², L] = 0 is a statement about the
*form* of L (it respects the mirror parity; the identity is F63 in the
[formula registry](../ANALYTICAL_FORMULAS.md)), not about where the noise comes
from. It is satisfied by the ordinary EXTERNAL Z-dephasing generator we always
use, so it cannot, by itself, rule out an internal origin. It is also
*axis-agnostic*: it holds equally for X-, Y-, and Z-dephasing; every dephasing
dissipator is diagonal in the Pauli basis and Π² is a sign on each Pauli
string, so the two commute trivially (verified from below in
[review2_A9_incompleteness.py](../../simulations/review2_A9_incompleteness.py):
‖[Π², L]‖ ≈ 1e-16 for all three axes). So the bootstrap result is
*underdetermination*, not impossibility: parity constrains the dissipator's
form to a multi-parameter family but does not forbid a self-generated source.
No other candidate below repairs this. Candidates 2 and 3 were built to,
and neither does.

**Source:** [bootstrap_test.py](../../simulations/bootstrap_test.py), [bootstrap_test.txt](../../simulations/results/bootstrap_test.txt); structural verifier [review2_A9_incompleteness.py](../../simulations/review2_A9_incompleteness.py)

**Conclusion:** a structural *constraint* (any noise must respect the parity
grading), not an elimination. The [Π², L] = 0 identity is origin- and
axis-agnostic, and nothing below supplies the elimination it does not give.

### Candidate 2: Single qubit decay (the failed third)

If the system as a whole cannot generate its own noise, maybe one piece
of it can. What if one qubit is unstable and its decay creates noise for
the others? This is the most intuitive candidate: a broken part shaking
the whole machine. We tested it with four different kinds of instability.

**Test:** Three-qubit system where one qubit (Q3) decays through four
different instability mechanisms (amplitude damping, detuning, thermal
bath, X+Y dephasing). Q1-Q2 have no external noise. Measurement: does
Q3's decay create effective palindromic dephasing on Q1-Q2?

**Result:** Q3's instability *does* damage Q1-Q2, and substantially: the
pair's purity falls from 1.000000 to 0.471109 (amplitude damping),
0.374586 (thermal bath) and 0.334815 (X+Y dephasing) over the run, and
the Bell coherence it was prepared in, |00⟩⟨11|, falls from 0.500000 to
0.143612, 0.110308 and 0.017248. Only detuning leaves it nearly untouched
(0.996993), which makes that fourth mechanism a near-identity channel that
can carry no weight either way.

**What the three-qubit system does.** Scored at its own forced centre, the
full generator is an **exact palindrome in every mechanism**: amplitude
damping at −0.050000000, the thermal bath at −0.075000000, X+Y dephasing at
−0.200000000, with relative power-sum residuals of 1.7e-17, 2.4e-17 and
1.8e-16 by a route that uses no eigensolver. Those are the centres
[F137](../ANALYTICAL_FORMULAS.md) predicts, −Σγ/2 for T1 and −Σ(γ↓+γ↑)/2 for
a thermal bath. The instrument is not vacuous: amplitude damping's spectrum,
whose centre is −0.050000000, scores 0.200 when it is offered −0.1 instead.

**What the two-qubit marginal does, and why it decides nothing about origin.**
A marginal over a spectator that is still coupled fails to pair for reasons
unrelated to where the noise comes from. Through one pipeline, minimising
over the centre in every row: a system with **no noise at all** gives 0.094;
**external** dephasing on Q1-Q2 with Q3 coupled gives 0.177, worse than any
internal mechanism; the internal mechanisms give 0.116, 0.125 and 0.149; and
the marginal pairs, at 2.2e-16, only when Q3 is fully decoupled. Same verdict
for external noise, internal noise, and no noise.

**Why no experiment of this shape can decide it.** The "internal" source is
modelled as a **Lindblad jump**, and a Lindblad dissipator is a coupling to
an external Markovian bath. The model grants the externality it was built to
test, which is why the generator's trace is negative and why a palindrome
exists at −0.05 at all. §3 takes this up: it is not a flaw in this
experiment but a property of the formalism.

**Reading `failed_third.py`.** The script is still in the repository and
still prints its March verdict, so three of its outputs need a warning. Its
`gamma_eff = 0` is a guard's else-branch: the fit normalises by |01⟩⟨10| at
t = 0, and a Φ⁺ preparation leaves that element **exactly** zero, so no
mechanism reaches the fit. Its "non-Markovian, ~50% of steps, max deviation
0.000000" counts increases of that same empty element beside a second
else-constant; the property itself does hold, and a BLP probe in the
reference below measures it. Its "0/16 palindromic pairs" is a centre
**search** on a grid coarser than its own tolerance, seeded away from the
centre F137 fixes in closed form, so it can pass only at its seed: hand that
call an exactly palindromic spectrum and it also returns 0/16. Compounding
it, the transfer matrix is divided by the Hilbert dimension twice, shifting
every rate by −ln(4)/t = −0.277259, which is the value standing twice at the
foot of all four of its published spectra. It also prints, twice and without
measuring anything, that no dephasing means no palindrome. §1 says otherwise,
and `simulations/two_qubits_no_noise.py` now measures it: the γ = 0 spectrum
pairs **16 of 16** about its own centre, zero.

**Source:** [incompleteness_candidate2_evidence.py](../../simulations/incompleteness_candidate2_evidence.py)
and [its results](../../simulations/results/incompleteness_candidate2_evidence.txt),
29 gates: the four marginals, the three-qubit palindromes at their forced
centres, the BLP probe, the guard, and the seeded search. It re-measures
[failed_third.py](../../simulations/failed_third.py) and
[its March output](../../simulations/results/failed_third.txt), which are
left as they were run.

**Conclusion:** not eliminated. The experiment measures a partial trace over
a coupled spectator, its "internal" source is an external bath by
construction, and the full generator it builds is an exact palindrome.

### Candidate 3: Many qubits (thermodynamic bath)

If one qubit cannot do it, what about many? A thermal bath is the
standard physics explanation for noise: the system sits in a large
environment of particles, and the environment shakes it. But within
a framework where everything is made of qubits, a bath is just a
collection of qubits, and each one faces the same prohibition.

**Argument:** a bath of qubits that is to be the ORIGIN of the noise needs a
reason for its own dynamics. Modelled as a Lindblad dissipator it has already
been given an environment; modelled unitarily it is a closed system whose
generator is traceless, so by §1 its palindrome sits at zero and carries no
decay. Either way the question moves outward a step instead of being
answered. That is the regress, and it is all this candidate has.

What it does not have is inheritance from Candidate 2, which eliminates
nothing. No N-qubit bath has been built and tomographed, and whether many
finite sources can compose into a palindromic channel is untested.

**Note:** none of this touches a bath as a mathematical DESCRIPTION of noise.
Lindblad theory works as a description, and §3 turns exactly on the fact that
it works as a description while saying nothing about origin.

**Conclusion:** not eliminated. The regress shows the ontology cannot close
the question from inside, which is weaker than an elimination and is the same
shape as §3's argument.

### Candidate 4: Nothing (d=0)

**Argument:** d=0 has no Hilbert space, no operators, no dynamics, no
properties. An entity with no properties cannot generate anything,
including noise.

**Conclusion:** Eliminated by definition.

### Candidate 5: Something other than qubits or nothing

**Argument:** d(d-2)=0 has exactly two solutions: d=0 (nothing) and
d=2 (qubit). The palindromic mirror symmetry condition requires that
Pi maps XY-weight k to N-k with balanced splitting. This requires the
d immune (diagonal) operators to balance the d²−d decaying ones,
d = d²−d. This holds only when
d²−2d=0. No other dimension satisfies the condition. Within the
framework, no entity with d>2 or non-integer d exists.

**Source:** [QUBIT_NECESSITY.md](../QUBIT_NECESSITY.md) (the counting argument); the typed [`PolynomialFoundationClaim`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) records the d=2 selection from d²−2d=0 with d=1 algebraically excluded; closes any "d=anything-else" loophole that would have been needed for a non-qubit noise source within the framework. The same polynomial is the c=0 case of F95's universal quadratic angle formula θ(c; b) ([formula registry](../ANALYTICAL_FORMULAS.md), F95), which tracks what happens when it is perturbed off the degenerate axis.

**Conclusion:** Eliminated by the framework's own algebra.

---

## 3. The Conclusion

**What is proven.** The commutator part of any generator is traceless:
trace(L_H) = −i(trace(H)·d − d·trace(Hᵀ)) = 0, for every H, since
trace(Hᵀ) = trace(H) always. Each jump costs the trace,
trace(D_F) = r(|trace F|² − d·trace(F†F)) ≤ 0 by Cauchy-Schwarz on ⟨I, F⟩,
with equality exactly when F is a multiple of the identity, and such a jump
has D ≡ 0. A sum of non-positive terms vanishes only if each does. So, **when
every rate is non-negative**:

> trace(L) = 0 ⟺ every jump carrying a nonzero rate is a multiple of the
> identity ⟺ the dissipator vanishes ⟺ the system is closed.

Non-negative rates are a *sufficient* condition for complete positivity, not
the definition of it: the Kossakowski matrix may be positive semidefinite with
a negative rate present, since non-orthogonal jumps can be recombined. The
clean form of the same statement is trace(L) = −d·tr(a) with a the Kossakowski
matrix, so a ⪰ 0 gives tr(a) = 0 ⟺ a = 0 directly.

Equivalently and more simply: the system is open exactly when some eigenvalue
lies off the imaginary axis. That is the whole certificate, and it needs no
palindrome.

**Where the palindrome comes in, and it is narrower than it looks.** If the
spectrum pairs, then its centre is trace(L)/dim identically (F137), so a
measured pairing lets you read the trace off a spectrum whose generator you
never had. It is an instrument for the trace, not a premise: a single
depolarizing site on a three-qubit chain has trace(L) = −19.2 and does not pair
about any centre, and certifies openness just as well.

**And it is exactly as unit-free as it looks, which is why it says so little.**
The Lindblad generator is homogeneous of degree one in all its rates JOINTLY,
L(αH, αγ) = αL(H, γ) identically, so (H, γ, t) and (αH, αγ, t/α) are the same
physics and no rate has invariant content on its own; only the dimensionless
ratio Q = J/γ and the product γt do. This is the joint rescaling and nothing
weaker: scaling γ alone moves Q and changes the trajectory's shape, which is
what [GAMMA_TIME_DISTINCTION.md](../GAMMA_TIME_DISTINCTION.md) measures when it
reports τ = γt failing to carry all observables. The centre −Σγ is therefore a
bookkeeping number, while **zero is the one rate value that survives every
rescaling**. That is why "centre = 0 against centre ≠ 0" is a question the
formalism can answer and "why is the centre −0.3" is not. The repo states the
same identity as [Q_BELONGS_TO_NO_SUBSTANCE.md](../Q_BELONGS_TO_NO_SUBSTANCE.md)'s
"γ₀ is the unit", L(J, γ₀) = γ₀·L₁(Q), and Corollary 2 below reaches it from the
physics side: the unit of γ is 1/[time], and without γ there is no time to
measure against.

**Both hypotheses are load bearing.** Drop complete positivity and the
equivalence fails at once: Z at rate +1 together with X at rate −1 gives
trace(L) = 0 exactly while ‖L‖_F = 2.83, and its spectrum {2, 0, 0, −2} even
pairs about zero. That generator has a positive eigenvalue, so it is not a
physical channel; it is the shape a time-local non-Markovian generator takes.
And "centre = 0" presupposes a spectrum that pairs at all, since without
pairing there is no centre to speak of. The chain that holds unconditionally
is the one about trace(L), not the one about the centre.

**What is not proven, and why nothing here could prove it.** The word above
is **open**, which is weaker than **external**. Every model of an internal
source in this document, and in `bootstrap_test.py` beside it, is written as
a **Lindblad dissipator**, and a Lindblad dissipator is a coupling to an
external Markovian bath: writing one down assumes the environment rather than
testing for it. The one alternative the formalism offers, an internal source
modelled unitarily, is a closed system, whose generator is traceless, whose
palindrome therefore sits at zero and carries no decay. The formalism leaves
no third option, so **inside it the origin question cannot be posed, let
alone settled.** γ enters a Lindblad generator as an input parameter, and
there is no slot in the formalism from which it could be derived.

That is the incompleteness. The framework describes all dynamics once noise
is given, and it cannot express a candidate origin without already granting
an outside. A proof of externality will not come from inside this formalism.
What one would need is a microscopic system-plus-bath derivation of L, and
this repository has never performed one.

### Corollary: the Arrow Is the Non-Unitarity

The identification is already in the framework:

- Without noise: unitary oscillation, reversible, no time arrow
- With noise: coherences decay irreversibly, creating a before and after
- The operator Π is literally time reversal: it maps exp(+mu*t) to
  exp(-mu*t), forward to backward ([PI_AS_TIME_REVERSAL.md](../../experiments/PI_AS_TIME_REVERSAL.md))
- The irreversible flow IS the time arrow: coherences (undecided,
  future) decay, populations (decided, past) persist, and the flow
  never reverses. Under the dephasing dissipator alone, the immune
  fraction of Liouville space is (1/2)^N (the all-{I,Z} Pauli
  strings; verified in
  [review2_A9_incompleteness.py](../../simulations/review2_A9_incompleteness.py)),
  everything else decays.

Noise is not like time. Noise IS the time arrow. Dephasing is what
makes processes irreversible. Irreversibility is what distinguishes
past from future. Without dephasing: perfect oscillation, no
direction, no history, no change. With dephasing: things happen,
things end, there is a before and an after.

§3 therefore has a direct corollary, and it is about existence rather than
origin. A palindrome centred away from zero certifies a non-unitary
generator; a non-unitary generator is exactly one with irreversible flow; and
irreversible flow IS the arrow. So **a measured palindrome with a nonzero
centre certifies that the system has an arrow of time**, and the certificate
is one number read off the spectrum.

It says nothing about where the arrow comes from. That question is §3's, and
§3's answer is that this formalism cannot ask it.

### Corollary 2: γ: Source of Experienced Time

The dephasing rate γ does not merely correlate with time. γ is the
source of experienced time: it provides the arrow. It is not identical
to experienced time; τ = γt does not scale universally, and the
Hamiltonian coupling J provides the content of what is experienced
([GAMMA_TIME_DISTINCTION.md](../GAMMA_TIME_DISTINCTION.md) carries the
precise statement table).

The same pattern appears wherever the framework touches time:

- The crossing time is t_cross = K/γ (K = 0.036 is the exact
  Bell+/concurrence value, not a universal constant; the historically
  quoted 0.039 was a tool's feedback-model reading, see
  [CROSSING_TAXONOMY.md](../../experiments/CROSSING_TAXONOMY.md)). The
  product t × γ = const is not a relation between two different
  quantities. It is a tautology by the Lindblad scaling symmetry:
  time multiplied by the rate of time gives a pure number.

- The unit of γ is 1/[time]. But this is circular: without γ, there
  IS no time to measure against. γ defines the scale against which t
  is counted. Remove γ and t loses its meaning.

- Π reverses t by reversing the sectors that γ acts on. It maps the
  immune sector {I, Z} (decided, classical, timeless) to the decaying
  sector {X, Y} (undecided, quantum, fragile). Reversing which sector
  decays reverses the direction of time. Decay IS direction. γ provides
  the arrow.

- In the transistor mapping ([GAMMA_CONTROL.md](../../experiments/GAMMA_CONTROL.md):
  the mediator's dephasing rate γ_M as gate signal): γ_M is the gate
  signal AND the clock.
  There is no separate clock line. The gate IS the clock. Because γ
  provides the irreversibility.

- The standing wave: exp(+mu*t) and exp(-mu*t) interfere to create a
  static pattern. But mu = lambda + Sigma_gamma, and Sigma_gamma = sum
  of all γ. The rescaled frame that removes the uniform decay envelope
  is the frame where γ has been factored out. What remains is timeless
  oscillation. The time arrow was γ. Remove it and irreversibility disappears.

This is the same circle §3 finds in the formalism, met from the physics side:
searching for the source of time using time is searching for γ using γ. The
instrument is identical to what it measures. A system cannot step outside
itself to find the origin of the thing that makes stepping possible.

γ is not merely a parameter of the system. γ is the source of the
system's experienced time: the arrow, with J shaping the content
([GAMMA_TIME_DISTINCTION.md](../GAMMA_TIME_DISTINCTION.md)). Where γ itself
comes from is open.

---

## 4. Structural Analogy

The following analogy is not part of the proof. It is context for readers
who know Gödel's incompleteness theorem and may recognize the pattern.
If you do not know Gödel, skip this section: the proof stands without it.

**Goedel (1931):** Any consistent formal system powerful enough to
express arithmetic contains true statements it cannot prove. The system
is complete for its domain but cannot prove its own consistency.

**d(d-2)=0 (2026):** The palindromic framework describes all decay
dynamics under dephasing exactly. But it cannot derive the existence of
dephasing from its own axioms. The system is complete for its domain but
cannot explain its own starting condition.

The structural parallel: both results identify a boundary of
self-reference. A system that is powerful enough to describe everything
inside itself is not powerful enough to explain why it exists.

This is a STRUCTURAL analogy, not a mathematical equivalence. Goedel's proof
is a theorem in mathematical logic. What sits here is an exact algebraic
identity (§3) plus an observation about what the Lindblad formalism can and
cannot express. The parallel is illustrative, not rigorous; if anything the
resemblance is closer than it was, since the obstruction §3 names is about
what the formalism can state rather than about what happens to be true.

---

## 5. What This Does NOT Say

- This does NOT prove the existence of God, a simulator, or any
  specific external entity.
- This does NOT prove the universe had an external cause.
- This does NOT establish what the noise IS.
- This does NOT prove that the noise comes from outside. The internal origin
  is open, and §3 gives the reason no argument inside this formalism closes
  it.
- This does NOT prove that d(d-2)=0 is the final framework. A more
  complete framework might resolve the incompleteness.
- **The enumeration of §2 is not exhaustive.** An internal d=2 source reached
  through a non-dephasing, measurement, or classical-field coupling sits
  inside the d(d−2)=0 ontology and outside the five cases;
  the nearest thing the repo has on that axis is narrower and is about a
  different dimension: [QUBIT_NECESSITY.md](../QUBIT_NECESSITY.md) §8a, §8d and
  §10 leave non-dephasing dissipators open **at d > 2**, and the words
  *measurement* and *classical field* do not occur there at all. The
  non-exhaustiveness stands on its own and is not carried by that citation. A
  five-case argument that is not exhaustive could not have concluded by
  elimination even if every case held.
- This DOES prove that a palindrome centred away from zero has a generator
  that is not unitary, exactly and in one line on the trace (§3).

---

## 6. Open Directions

§3 says the origin question cannot be posed inside this formalism. It does
not say the question is unanswerable. Four directions remain: the first is
concrete work, the other three are stances on what to do with a boundary.

0. **Derive L microscopically.** Purify, couple the system to an explicit
   bath, evolve unitarily, trace out, and see what has to be true of the bath
   for the reduced generator to be the dephasing one. That is the one route
   that escapes §3's circle, because it never writes a dissipator down as an
   assumption. It has never been done here.

1. **The framework is incomplete.** d(d-2)=0 is not the full equation.
   A more complete theory, perhaps involving d values other than 0
   and 2 or a different algebraic structure, might contain its own
   noise source. This is the "extend the framework" path.

2. **Noise is axiomatic.** Like the speed of light or Planck's constant,
   the existence of dephasing is a brute fact of reality that cannot be
   derived from anything more fundamental. It is a starting condition,
   not a consequence. This is the "accept the boundary" path.

3. **The question is malformed.** Asking "where does noise come from"
   presupposes a causal chain. But if noise and qubits are
   co-fundamental (neither causes the other, both are aspects of the
   same thing), the question dissolves. This is the "reframe the
   question" path.

This document does not choose between the last three. What it establishes is
narrower than any of them: the system is open, and the formalism cannot say
more.

A fourth direction emerged from the algebraic analysis of the Urqubit
(April 1, 2026): the cross term {L_H, L_D + Σγ·I} vanishes exactly at
N=2 and is nonzero at N > 2. This means oscillation (Hamiltonian) and
cooling (dissipator) are Frobenius-orthogonal only for the single bond.
At N > 2, they are woven together, and that weaving is not undone by
reduction: the marginal over a coupled spectator fails to pair whether the
noise is external, internal, or absent (§2, Candidate 2). Reading this loss
of orthogonality as an arrow-of-time exclusion is a Tier-3 interpretation,
not a time-reversal theorem: the dynamical
separability criterion is the commutator [L_H, L_Dc], nonzero at every
N. See
[Time Irreversibility Exclusion](TIME_IRREVERSIBILITY_EXCLUSION.md) and
[Primordial Qubit Algebra](../../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md).

---

## Verification

Each step is independently verifiable:

1. Read [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md): the palindrome is real.
2. Check §3's identity by hand. It is two lines: trace(H⊗I) = trace(I⊗Hᵀ) for
   the commutator part, and Cauchy-Schwarz on ⟨I, F⟩ for the jump. No
   computer is required, and no dimension is assumed.
3. Run `python` [`simulations/incompleteness_candidate2_evidence.py`](../../simulations/incompleteness_candidate2_evidence.py):
   29 gates. It measures the closed system's pairing at zero (6d), the four
   marginals that do not separate origin, the three-qubit palindromes at
   F137's centres, a depolarizing site that certifies openness with no
   palindrome at all (6e), and the trace identity over random Hermitian H and
   random jumps, including the counterexample showing that non-negative rates
   are load bearing (7d).
4. Run `python` [`simulations/bootstrap_test.py`](../../simulations/bootstrap_test.py):
   the parity sectors decouple. Read it as the constraint it is, not as an
   elimination (§2, Candidate 1).
5. Read [Qubit Necessity](../QUBIT_NECESSITY.md): only d=0 or d=2, and note
   §8a, §8d and §10, which keep non-dephasing dissipators open at d > 2. That
   is a neighbouring question, not §5's; §5's axis has no source here.

`simulations/failed_third.py` still runs and still prints its March verdict.
Read §2, Candidate 2, "Reading `failed_third.py`" before believing any of its
three headline numbers.
