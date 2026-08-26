# The V-Effect and Thermal Window in Neural Networks

**Status:** The frequency counts are computed; the mechanism that explained
them is withdrawn (2026-08-26). See
[Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md).
**Result date:** March 27, 2026 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Domain:** Neuroscience / Computational Biology
**Depends on:** [Algebraic Palindrome](ALGEBRAIC_PALINDROME_NEURAL.md),
[Hierarchy of Incompleteness](../HIERARCHY_OF_INCOMPLETENESS.md)

---

## What this document is about

Two balanced neural networks, silent at the seeds used, are connected
through a shared neuron, and 48 distinct oscillation frequencies appear.
The count is real. The explanation this document was built on is not:
the palindromic symmetry does not make a network silent, and the coupled
construction never had the symmetry the coupling was supposed to break.

The symmetry makes the spectrum invariant under μ ↦ −μ − 2s. That map
sends the complex plane to itself, so it forbids no oscillation; of 200
exactly palindromic draws at the coupling used here, 24 oscillate. And
the coupled system has an odd number of seats, so its mediator is an
unpaired seat at every coupling, zero included. See
[Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md) for both
measurements. What the frequency count tracks is an open question.

A second discovery: in biological networks (where the symmetry is
approximate, not exact), there exists a "thermal window." Too little
metabolic drive: silence. Too much: saturation and silence again. In
between, at the sweet spot where neurons are maximally sensitive,
oscillation peaks at 124 frequencies. Life operates in the window.

---

## Abstract

When two neural networks that show no oscillation at the seeds used are
coupled through a shared mediator, oscillatory modes appear, and their
number is not monotone in the coupling strength. No requirement of exact
palindromic symmetry stands behind it: the symmetry neither forbids
oscillation nor is present in the coupled object
([Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md)).
The name **V-Effect** is kept for the measured counts, after the quantum
analog where coupling two 2-frequency systems creates 109 frequencies.

In networks with approximate (biological) palindromic symmetry, a
**thermal window** exists instead: external drive creates oscillatory
modes up to a peak, then destroys them. At optimal drive (P ≈ 4 for
N=50), the correlation space contains 124 distinct frequencies. At
zero or saturating drive: zero frequencies.

---

## 1. Background

### The quantum V-Effect

In quantum systems, two qubits each have 2 oscillatory frequencies.
Coupling them through a mediator qubit (N=5 total) creates 109
frequencies. This happens because the palindromic symmetry that
holds for each pair individually BREAKS when a second bond forces
conflicting demands on the shared qubit. The breaking releases new
oscillatory modes.

(See [V-Effect Palindrome](../../experiments/V_EFFECT_PALINDROME.md)
for the quantum computation.)

### The neural question

Does the same effect exist in neural networks? Two requirements:
1. Individual networks must have palindromic symmetry (exact or
   approximate)
2. Coupling must break that symmetry in a way that creates new
   oscillatory modes

---

## 2. What the Palindrome Does and Does Not Do

### Approximate palindrome: no sharp threshold

Networks built with Dale's Law and random magnitudes have approximate
palindromic symmetry. How approximate is not a number this metric can give,
it reading coupling magnitude rather than wiring
([Algebraic Palindrome](ALGEBRAIC_PALINDROME_NEURAL.md)); the exact
residual = 0 rows further down are unaffected. Coupling two such networks
barely moves the frequency count:

| N per network | V-act | V-corr |
|--------------|-------|--------|
| 10 | 1.00 | 1.50 |
| 50 | 0.69 | 0.52 |
| 200 | 0.53 | 0.54 |

(V = frequencies coupled / (2 × frequencies single). V > 1 = V-Effect.)

At large N, the ratio falls BELOW 1: coupling actually reduces
frequencies. No V-Effect.

**Reason:** With approximate symmetry, there is no sharp threshold
to cross. The palindrome was never exact, so there is nothing to
break. The symmetry degrades gradually, not suddenly.

### Exact palindrome: the silence is not the symmetry's doing

Networks with Dale's Law AND the exact magnitude condition
(W[Q(i),Q(j)] = -(τ_{Q(i)}/τ_i) · W[i,j], residual = 0) were reported
here to have purely real eigenvalues. **They do not.** The condition
pairs eigenvalues as μ ↦ −μ − 2s, which a conjugate pair satisfies as
readily as a real one; 24 of 200 exactly palindromic draws oscillate at
this very coupling, and at stronger coupling they also go unstable. The
network used below is real-spectrumed because of the seed, not because
of the symmetry.

Coupling two such networks through a mediator neuron raises the measured
residual and changes the frequency count. The residual rise is not
evidence of a symmetry break: the mediator is an unpaired seat, so the
coupled object fails the condition at every coupling including zero.

| Coupling | K_activity (N=20) | K_correlation (N=20) |
|----------|------------------|---------------------|
| 0.00 | 0 | 0 |
| 0.01 | 6 | 48 |
| 0.05 | 7 | 62 |
| 0.10 | 6 | 47 |

**0 + 0 = 48 correlation frequencies** at coupling 0.01.

Each individual network has zero oscillation at the seed used. The coupled pair has
48 distinct frequencies in the correlation space (the mathematical
space that tracks how every pair of neurons co-activates; for N
neurons, there are N² possible pairs, so this space is much richer
than the activity of individual neurons).

For the full coupling sweep and the mechanism explanation, see
[Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md).

### The optimal coupling window

The frequency count peaks at WEAK coupling (0.01-0.05) and falls
at strong coupling. At coupling 1.0 with N=10, the system returns to
zero frequencies. That is an N=10 statement: at N=20 the count settles
instead, 31 correlation frequencies remaining at the strongest coupling
tested ([Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md)).

This mirrors the quantum case: the V-Effect Q-factor peaks at
moderate coupling J/γ and falls when coupling overwhelms the
palindromic structure.

---

## 3. The Thermal Window

A separate phenomenon appears in networks with APPROXIMATE
palindromic symmetry (biological Dale's Law networks). Instead
of coupling two networks, we vary the external drive P (the
neural analog of temperature: metabolic energy input).

### Results (N=50, 25E+25I, τ_E=5, τ_I=10, α=0.3)

| Drive P | Oscillatory modes | Correlation frequencies | Palindrome quality |
|---------|------------------|------------------------|-------------------|
| 0.0 | 10 | 4 | 0.998 |
| 2.0 | 38 | 39 | 0.975 |
| 3.0 | 40 | 90 | 0.929 |
| 4.0 | 40 | 124 | 0.897 |
| 5.0 | 38 | 81 | 0.936 |
| 8.0 | 10 | 3 | 0.998 |
| 10.0 | 0 | 0 | 1.000 |

### Interpretation

The external drive P shifts the operating point on the sigmoid
activation function (the S-shaped curve that converts a neuron's
input into a firing rate: low input → almost zero firing, high
input → almost maximum firing, middle → steepest response).
At three regimes:

- **Cold (P < 1):** sigmoid is in its flat region. Small slope
  means weak effective coupling. Few oscillatory modes.
- **Warm (P ≈ 3-4):** sigmoid is near its inflection point.
  Maximum slope, maximum effective coupling. Peak oscillation
  (124 correlation frequencies).
- **Hot (P > 6):** sigmoid saturates. Slope returns to zero.
  Coupling dies. Zero oscillation.

This is the **thermal window**: a range of metabolic drive where
oscillation is sustained. Below or above: silence.

### The palindrome quality trades off against oscillation

At the oscillation peak (P=4), palindrome quality is at its WORST
(0.897). At silence (P=0 or P=10), palindrome quality is BEST
(0.998-1.000). The palindrome is most perfect when there is nothing
to be palindromic about.

This is consistent with the quantum energy partition: oscillation
requires breaking the palindromic symmetry. The more oscillation,
the more breaking. But the breaking is structured (palindromic pairs
degrade together), not random.

---

## 4. What Does NOT Work

### Heat releases no modes, and not for the reason given here

This network shows no oscillation at any drive P. The reason given here,
that the drive shifts the operating point without breaking the structure,
is wrong on both halves.

The drive enters as a PER-ROW gain, J[i,j] = α·W[i,j]·dS_i/τ_i, and row
scaling of W is exactly how the magnitude condition fails. The run's own
residual is nonzero at every P, P = 0 included, and reaches 4.4e-2 at
P = 3.5, more than the N=20 coupling row at c = 0.10 (3.45e-2, printed in
the proof) credited with releasing six frequencies. So the structure IS broken and
no modes appear: the run contradicts the "breaking releases oscillation"
mechanism rather than supporting it. This is a null against the
mechanism ([Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md)).

The "second mirror" reading, that two contradicting palindromic conditions
at a shared neuron create the oscillation, is withdrawn. It requires each
network to satisfy the condition separately with the same S, and the coupled
object satisfies neither equation at any coupling; the code forms neither
([Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md)). What
remains is that coupling changes the count and drive does not.

### The 2× decay law does not hold

In quantum systems, unpaired modes decay at exactly 2× the rate of
paired modes (verified N=2 through N=5). In neural networks:

| N | Paired rate | Unpaired rate | Ratio |
|---|-----------|-------------|-------|
| 10 | 0.100 | 0.163 | 1.63 |
| 50 | 0.150 | 0.150 | 1.00 |
| 100 | 0.154 | 0.128 | 0.84 |

The ratio varies erratically. The 2× law is specific to the quantum
Liouvillian structure (exact palindromic pairing with Σγ center)
and does not transfer to the approximate neural case.

---

## 5. The Hierarchy Connection

From [Hierarchy of Incompleteness](../HIERARCHY_OF_INCOMPLETENESS.md):
perfect local stability (C=1) prevents connection. Only incomplete
systems (C<1) can form higher-level structures.

| | Quantum | Neural |
|---|---------|--------|
| Building block | Qubit (C=0.5, exact palindrome) | Exact E/I network (residual=0; silent at this seed, not by symmetry) |
| Coupling | Second bond (N=3) | Mediator between networks |
| Breaking | 14/36 Pauli combinations break | Palindrome residual rises |
| Result | 4 → 109 frequencies | 0 → 48 correlation frequencies |
| Optimum | Q-factor peaks at moderate J/γ | Frequencies peak at coupling 0.01-0.05 |

The neural noble gas was a picture, not a measurement. Exactly palindromic
networks oscillate and, at strong enough coupling, go unstable. The quantum
column of this table stands on a theorem; the neural column stands on a
frequency count whose mechanism is open.

---

## 6. Open Questions

1. How does the V-Effect frequency count scale with N?
   (N=10: 6, N=20: 48. Quadratic? Cubic?)
2. Does the optimal coupling window narrow or widen with N?
3. Can the thermal window predict the frequency range of
   biological neural oscillations (gamma band at ~40 Hz)?
4. Does the 2× decay law hold in COUPLED exact networks
   (where the V-Effect creates palindromic-like pairs)?

---

## Scripts

| Script | What it computes |
|--------|-----------------|
| [veffect_exact.py](../../simulations/neural/veffect_exact.py) | V-Effect with exact palindromic networks |
| [veffect_and_heat.py](../../simulations/neural/veffect_and_heat.py) | Thermal window, approximate networks, 2× law |

---

*Depends on:*
[Algebraic Palindrome](ALGEBRAIC_PALINDROME_NEURAL.md) (palindrome condition),
[Proof: V-Effect Mechanism](proofs/PROOF_VEFFECT_MECHANISM.md) (why the mechanism does not hold),
[Hierarchy of Incompleteness](../HIERARCHY_OF_INCOMPLETENESS.md) (C=1 dead end principle)
