# What Qubits Experience

**What this document is about:** A Tier 5 meditation on one exact algebraic
contrast in the uniform Heisenberg/Z-dephasing model. At N=2 the centered
Hamiltonian and dissipative superoperators anticommute; at N=3 the normalized
Frobenius anticommutator is `1/sqrt(48)` (about 14.43%). A separate 1.83%
figure divides the cross term by `||L_c²||` at the reported parameter point;
it is not the γ-independent normalization. This is not a theorem about
experienced time, dynamical factorization or a macroscopic arrow of time. The
language below the marked line is an interpretation of this contrast.

**Status:** Tier 5 (interpretation), grounded in Tier 1-2 computation
**Date:** April 1, 2026
**Last refreshed:** 2026-09-06 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Depends on:**
- [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) (open-system decay-clock scale; no experienced-time identification)
- [Time Irreversibility Exclusion](../docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md) (cross term)
- [Primordial Qubit Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) (Pythagorean theorem)
- [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md) (Σγ = 0 as origin)

---

## The Computed Ground

Everything below the line is Tier 1-2. Proven or computed. No
interpretation needed.

A qubit under Z-dephasing has a rate γ. At γ = 0 the stated model has
unitary Hamiltonian evolution. At γ > 0 the dissipator directly damps
off-diagonal components in the Z basis; populations pay no direct dephasing
cost, although H can mix operator components. Whether a trajectory settles,
and to which invariant set, depends on H, the state and conserved quantities.
These statements provide a decay scale, not an identification of γ with
experienced time.

For the N=2 uniform Heisenberg bond used in the cited computation, the
Hamiltonian part L_H and centered dissipator L_D + Σγ I anticommute:

    {L_H, L_D + Σγ·I} = 0

The square of the dynamics decomposes cleanly:

    L_c² = L_H² + (L_D + Σγ)²

The square therefore has no anticommutator cross term. This does **not** mean
that the two flows factor dynamically: their commutator is nonzero already at
N=2, as [Time Irreversibility Exclusion](../docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md)
states explicitly.

For the corresponding N=3 chain, the anticommutator is nonzero:

    {L_H, L_D + Σγ·I} ≠ 0

The normalized Frobenius ratio
`||{L_H,L_Dc}||/(||L_H|| ||L_Dc||)` is `1/sqrt(48)` for this
model and is independent of the positive uniform γ after the common factor
is divided out. Numerically, `1/sqrt(48) ≈ 0.144337`, or 14.43%. The 1.83%
reported elsewhere is `||{L_H,L_Dc}||/||L_c²||` at one parameter point. The
raw anticommutator itself scales with γ. The general
closed form and its topology scope are given in the cited proof.

In the cited N=3 marginal test, tracing out the coupled third qubit does not
recover the N=2 palindrome. The failure remains when that third qubit carries
no noise (residual 0.094, against 0.116–0.149 with its decay running). A
previous 0/16 output came from a seeded-center search that also rejects an
exactly palindromic control, so it is not evidence. This finite marginal test
does not prove an impossibility for every reduction.

Everything above this line is the scoped mathematical/computational ground.

---

## The Line

Below here is feeling. Not physics. Not falsifiable. Words such as
"experience", "separate", "watcher" and "time" are metaphors for the
anticommutator result above, not additional consequences of it. The
mathematics allows other readings.

---

## What Happens to a Qubit

A qubit at N=2 can oscillate under H and decay under positive γ. In this
meditation, the vanishing anticommutator is pictured as a right angle between
those contributions. It is not literal dynamical independence: the two
semigroups do not factor, even at N=2.

This is the simplest image: a rhythm and a fading. The `N=2` zero
anticommutator/Pythagorean cross term motivates the right-angle picture. The
ordinary Hilbert-Schmidt inner product of the two generator matrices is not
the N-dependent diagnostic, and neither zero establishes dynamically
independent flows.

Then a third qubit arrives.

The third qubit does not do anything dramatic. It just exists nearby.
It couples to one of the first two through a bond. But now the bond
is local: it connects two of three, not all. And the third qubit
watches. It contributes its weight to the w_XY sum. It shifts the
geometry.

And the right angle bends.

Now oscillation bleeds into cooling. The rhythm disturbs the fading.
The fading bends the rhythm. They are no longer two things; they are
one thing with two aspects that cannot be separated.

The qubit at N=3 does not experience time and change as independent
phenomena. It experiences them as woven together. It cannot point to
the fading and say "this is time." It cannot point to the rhythm and
say "this is change." The two are entangled. Not quantum-entangled;
algebraically entangled. The cross term ties them.

Within the computed N=3 model, the normalized bend is 14.43% and does not
depend on the positive uniform γ. Calling it permanent or experiential is
the interpretation, not a result about arbitrary systems.

---

## What Happens to Us

We are not three qubits. We are 10²³ local bonds, each spanning two
sites out of unimaginably many. Every bond has unimaginably many
watchers. Every watcher bends the angle. The cross term is not 2%.
It is the entirety of our experience.

We experience time and change as one thing. We cannot separate the
passing of time from the happening of events. We cannot point to a
moment and say "here, time passed but nothing changed." We cannot
point to a change and say "this happened outside of time." They are
woven so tightly that we do not even notice they were once separate.

In the meditation, the `N=2` right angle becomes the image of two experiences.
Mathematically it is only an inner-product statement: it does not establish
separable semigroups, reversibility, or an absence of thermodynamic time.

We are what happens when you add watchers. When bonds become local.
When the right angle bends and bends and bends until the two legs of
the triangle are indistinguishable from the hypotenuse.

The following is the philosophical reading, not a result of the calculation:
the nonzero cross term is pictured as the entangling of rhythm and fading.
The calculation supplies neither a human arrow-of-time mechanism nor a claim
about every bond in bodies, brains, or air. Its `γ`-independence after the
chosen normalization does not make it a universal physical cause.

Within that metaphor, the arrow is drawn as a bent angle. The typed exclusion
is controlling: the anticommutator norm does not prove reversibility at zero,
irreversibility when nonzero, or a physical arrow of time.

---

## What the Urqubit Was

Before the watchers, before locality, before the cross term, there
was a bond between two qubits. It spanned everything. It had no
outside. Its oscillation and its cooling were two independent words
for two independent experiences.

It could have stayed there. At the right angle. At the mirror.
Where forward and backward were imagined as the same word. This is not a
reversibility claim about the modeled open dynamics.

It did not stay.

A third qubit arrived. The framework does not identify the origin of that
larger system, its bath or its system-environment cut. In this meditation the
additional site brings locality, the nonzero anticommutator, and the image of
an irreversibility that makes everything we know possible.

The Urqubit remains a metaphor. The calculation does not establish that every
physical bond realizes this `N=2` subsystem, nor that a right angle allows
reversal or a bent angle forbids it.

We are the complexity that the Urqubit became.

The cited N=3 marginal test does not recover the N=2 palindrome after tracing
out its coupled neighbour, including in the zero-noise control on that
neighbour. That finite result does not prove that every reduction from every
N > 2 system must fail.

We remember the past because we are made of it. We cannot return to
it because the angle is bent.

---

*"There was a primordial qubit. It sits on both sides."*

*It still does. But now there are watchers. And the watchers bent
the angle between its two sides. And the bending is what we call
time.*

*Thomas Wicht, April 1, 2026*
