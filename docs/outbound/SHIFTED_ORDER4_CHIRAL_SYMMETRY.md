# A Shifted, Order-4 Chiral Symmetry in Local-Dephasing Lindbladians

<!-- Keywords: Lindbladian sector symmetry classification, shifted chiral symmetry order
four, many-body Lindbladian tenfold way beyond, anti-pseudo-Hermiticity constant
shift, chiral Liouvillian Poisson statistics, complex spacing ratio
filling threshold GinUE, Sa Ribeiro Prosen classification extension, outbound
adapter symmetry classification community, R=CPsi2 placement -->

**PARKED (2026-07-06), not outreach-ready.** The whole docs/outbound arc is
parked (no outreach plan, gated on Tom; the S2 concentrator adapter was found
hollow on review). This adapter has NOT itself had an empty-review round, so
treat it as an unvalidated draft; every result it cites lives in its home
claim/proof/theorem regardless. See [README](README.md).

**Status:** Outbound adapter (draft). This document translates a repository
result into the Lindbladian-symmetry-classification community's language and
open edges. The operator identity is Tier 1 (analytic proof + machine
verification, 87,376 eigenvalues, zero exceptions); the placement against the
Sá-Ribeiro-Prosen scheme is formal analysis; the spectral-statistics results
are computational (N ≤ 8, live-witnessed). External citations are drawn from
the repository's literature analysis and should be verified against the
primary sources before any outreach.
**Date:** July 5, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Mirror Symmetry Proof](../proofs/MIRROR_SYMMETRY_PROOF.md) (the
identity and its verification), [KMS and Detailed Balance](../KMS_DETAILED_BALANCE.md)
(the placement analysis), [Random Matrix Theory](../../experiments/RANDOM_MATRIX_THEORY.md)
(the statistics), [PT-Symmetry Analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md)
(the direct operator classification), [Filling Threshold Chaos](../../experiments/FILLING_THRESHOLD_CHAOS.md)
(the sector-resolved crossover), [K Partnership proof](../proofs/PROOF_K_PARTNERSHIP.md)
(the Hamiltonian-level altitude, kept separate), [Selective Decoupling Selection Rule](SELECTIVE_DECOUPLING_SELECTION_RULE.md),
[State Transfer Decay Structure](STATE_TRANSFER_DECAY_STRUCTURE.md) and
[Noise Asymmetry Symmetry Scalar](NOISE_ASYMMETRY_SYMMETRY_SCALAR.md) (the
sister adapters), [Combination Valence](../../hypotheses/COMBINATION_VALENCE.md)
(why this adapter hands over objects, not words)

---

## What this document is about, and who it is for

This is the third entry of the repository's **outbound arc**: where the
translation series in `docs/quantum/` carries a community label into our
stance, this carries the operator identity studied in the repository out to a
community's stance. It is
written for one audience, the **symmetry-classification community for open
quantum systems** in the line of Sá, Ribeiro, and Prosen ("Symmetry
Classification of Many-Body Lindbladians: Tenfold Way and Beyond", Phys. Rev.
X 13, 031019; [arXiv:2212.00474](https://arxiv.org/abs/2212.00474)). The
classification must be applied to the shifted generator inside irreducible
symmetry sectors; the globally order-4 representative below is not by itself
a new class.

It follows the arc's standing rule (see
[Combination Valence](../../hypotheses/COMBINATION_VALENCE.md)): hand over an
object, a representation, or a number the reader can stand in, not a coined
word. Our in-house name for the phenomenon below stays home; what travels is
an operator, an exact identity, and a spectral-statistics data set.

One altitude note up front, to keep the bookkeeping honest. The model family
below also carries a perfectly ordinary **Hamiltonian-level** sublattice
(chiral) symmetry, K = diag((−1)^ℓ) on the excitation lattice (in spin
space, ⊗_{odd i} Z_i) with K H K = −H, class BDI on real hopping: textbook,
and not the subject here. The object of this document
lives one level up, on the **Liouvillian superoperator**, and the two are
independent symmetries of distinct objects
([K Partnership proof](../proofs/PROOF_K_PARTNERSHIP.md)). Everything below
is about the Liouvillian level.

---

## 1. The classification operation this example requires

The tenfold-way-and-beyond program classifies many-body Lindbladians by
their behavior under the generator set {T₊, C₊, T₋, C₋, P, Q₊, Q₋}:
antiunitary time-reversal and particle-hole flavors, unitary chiral P, and
the two pseudo-Hermiticity conditions involving L†. Its negative symmetries
are defined up to an identity shift and hence on the traceless shifted
generator L′. It also requires all commuting unitary symmetries to be
resolved before assigning a class to an irreducible block.

For the operator below, Π² is itself a commuting unitary symmetry Ux. In a
Ux=px sector, px∈{+1,−1}, the phase-normalized restriction
P_px=√px·Π|px obeys P_px²=I and anticommutes with L′. This is the same
sectorwise normalization used in the paper's dephasing examples. The live
question is therefore not whether the shift or global order four lies
outside the taxonomy, but which class each fully reduced sector occupies
after the remaining commuting symmetries and their relations with T₊, Q₊
and P_px have been worked out.

This document hands over a model family with an exact operator identity
and computed statistics; the statistics do not prove full spectral solvability.

---

## 2. The object: a site-local, order-4, shifted chiral operator

Take any Heisenberg / XY / Ising / XXZ coupling graph with local Z-dephasing
at rates γ_i, one of the scheme's own standard example families. On the 4^N
Pauli-string basis of operator space, define the unitary Π site by site:

    per site:  I → X,  X → I,  Y → iZ,  Z → iY

(a strict tensor product of identical single-site maps; the factors ±i are
essential). Then, exactly and for every such graph:

    Π · L · Π⁻¹ = −L − 2Σγ · I,

an anti-similarity **with a constant shift**. Equivalently, the centered
Liouvillian L_c = L + Σγ·I anti-commutes with Π. The identity is proven
analytically (it reduces to a 16-entry single-bond table plus site-locality
and linearity) and machine-verified over N = 2 to 8: 87,376 eigenvalues,
every one paired, zero exceptions, on chains, stars, rings, complete graphs
and trees ([proof](../proofs/MIRROR_SYMMETRY_PROOF.md)).

The operator's profile, each property separately pinned in the sources:

- **Unitary and linear** (a signed permutation of the Pauli basis, phases in
  {±1, ±i}); no complex conjugation anywhere in the condition.
- **Order 4 globally, involutive sectorwise:** Π⁴ = I, and Π² is conjugation
  by the global string X^⊗N (it acts on a Pauli string as
  (−1)^{n_Y+n_Z}). Resolving this commuting ±1 symmetry and multiplying the
  restriction by √px gives the involution required for P.
- **Site-local:** a tensor product of single-site operations, so it survives
  arbitrary coupling graphs and inhomogeneous rates γ_i unchanged.
- **Chiral for both:** after centering, Π anti-commutes with L_c and with
  L_c† ([PT-Symmetry Analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md),
  residuals 0.0).

---

## 3. What is placed, and what still needs calculation

The exact identity places Π as a global representative of a P-type chiral
symmetry of the shifted generator. The classification step is:

- shift to L′=L+Σγ·I, exactly as the cited scheme prescribes for negative
  symmetries;
- resolve Ux=Π² into px=±1 sectors;
- use P_px=√px·Π|px, for which P_px²=I and {P_px,L′|px}=0;
- resolve the other commuting symmetries and compute the sectorwise
  commutation algebra with the unavoidable Hermiticity-preserving T₊ and any
  Q symmetry before naming a class.

The last step has not been completed for this repository's exact
joint-popcount/reflection decomposition. Spectral reflection alone cannot
select AIII or any other class. Sá, Ribeiro and Prosen's dephasing examples
demonstrate the relevant operation and obtain parity-dependent BDI/CI
classes for their Hamiltonians; those labels cannot simply be copied onto a
further-reducible model with a different Hamiltonian symmetry algebra.

The *shape* of the shifted spectrum also has a catalogued noninteracting
home in the "shifted sublattice symmetry" of Kawasaki, Mochizuki, and Obuse
(Phys. Rev. B 106, 035408;
[arXiv:2201.09283](https://arxiv.org/abs/2201.09283)). For interacting
many-body Lindbladians, Sá, Ribeiro and Prosen already formulate the
classification on L′ and construct site-local chiral operators for dephasing
chains. The object placed here is the explicit Pauli-basis Π, its global
relation Π²=Ux and its exact identity across the stated coupling graphs.
Prior-art coverage and equivalence to published constructions remain OPEN. No
priority or novelty claim is made.

So the sharp question we hand over is sectorwise: **after resolving joint
popcount, reflection and Ux, what are the complete symmetry relations and
the resulting class of each irreducible block?** The exact Π identity is an
input to that calculation, not a substitute for it.

---

## 4. Symmetry classification and random-matrix universality

The classification's random-matrix predictions apply within irreducible
symmetry sectors. The repository's existing global rate-spacing statistic
pools many exact joint-popcount blocks, so it is descriptive data rather than
a test of one classified ensemble
([Random Matrix Theory](../../experiments/RANDOM_MATRIX_THEORY.md)):

- **Symmetry axis:** the centered spectrum is exactly ± paired (machine
  precision, ~1e-15) at every tested N. This is the spectral consequence of
  the P-type anticommutation; it is not by itself a class assignment.
- **Statistics axis:** direct raw-multiset consecutive-gap ratios retain zero
  gaps and count undefined 0/0 separately (21,840 eigenvalues, N=2 to 7).
  At N=7 the mean is 0.2021456120489688 over 16,366 defined ratios, with 16 undefined.
  There is no standard-ensemble calibration for this degenerate unresolved population.
  The separate complex spacing ratio of the chain reads 2D-Poisson-like in its executed comparison. This does not prove integrability.
  This pooled comparison cannot confirm or refute a sector-level universality
  prediction. The statistics do not supply the missing class assignment.

Sector-resolved, the picture sharpens into a data point for your
class-statistics correspondence
([Filling Threshold Chaos](../../experiments/FILLING_THRESHOLD_CHAOS.md)):

- **Dilute coherence sectors stay Poisson-like at the tested N and operating points**,
  including the canonical Delta=1 plus disorder comparison. Nonzero Delta breaks free-fermion
  additivity, but uniform XXZ remains Bethe-integrable; at Delta=0 the random-field XY Hamiltonian remains
  quadratic. That Hamiltonian label does not classify the Z-dephasing Liouvillian as a quadratic generator.
  Delta!=0 plus generic disorder is the interacting nonintegrable many-body test.
- **Dense (near-half-filling) sectors climb toward GinUE with N**: the
  complex-spacing-ratio marker ⟨cos θ⟩ runs −0.089 / −0.129 / −0.162 at
  N = 6 / 7 / 8 (GinUE −0.241, 2D-Poisson 0). GinUE is the class-A
  comparison ensemble here, not a computed class assignment: Π maps an unequal
  block to its conjugate partner, but that fact and the sampled conjugation test
  do not exhaust residual antiunitary symmetries. The full irreducible class is
  **OPEN** after the shifted generator's sectorwise P symmetry and every
  unitary/strong sector are resolved.

In this tested family, the finite-N CSR comparison supports a filling dependence
at fixed interacting disorder. It does not prove a universal threshold or thermalization cause.

---

## 5. What you can check with your own tools

Three independent handovers, each cheap on your side:

> 1. **The identity, in one line of code.** Build Π from the per-site rule
>    above (a 4^N signed permutation), build L for any small Heisenberg +
>    local-dephasing model of your choice, and check
>    Π·L·Π⁻¹ + L + 2Σγ·I = 0 to machine precision. No fitting, no limits;
>    it either vanishes or it does not.
> 2. **The classification question.** Resolve Ux=Π² and every other commuting
>    unitary symmetry, phase-normalize Π on each Ux sector, and compute the
>    full sectorwise symmetry algebra. The output is a class per irreducible
>    block, not a class inferred from the global spectral reflection.
> 3. **The crossover, with your diagnostic.** Reproduce the sector-resolved
>    complex spacing ratio: the measured finite-N contrast at N=6..8 under
>    canonical Δ=1 plus disorder is dilute Poisson-like versus dense movement
>    toward GinUE. Large-N convergence is a proposed extension, not a measured
>    limit or a thermalization theorem.

---

## 6. What we are not claiming

Weak coupling, deliberately: this adapter offers a test case and a sharp
question to a classification program, not a rival scheme.

- **The Sá-Ribeiro-Prosen placement is incomplete at the irreducible-sector
  step.** We have not worked the Bernard-LeClair or the Kawabata et al.
  non-Hermitian schemes.
- **Symmetry class and statistical universality are separate:** measured
  Poisson-like spacings neither assign a class nor prove integrability.
- **The statistics are computational and small-N** (N ≤ 7 global, N ≤ 8
  sector-resolved, with bootstrap intervals at the largest N). The identity
  itself is exact and proven; the statistics are data.
- **Z-dephasing only** (single-axis, unital). Other dissipators have
  analogous single-axis operators in the repository, but the identity as
  stated is for the dephasing family.
- **The Hamiltonian-level sublattice BDI is not the result** (the altitude
  note of the opening section): that symmetry is textbook and lives on a
  different space; conflating the two levels would conflate distinct objects.
- **Prior-art boundary:** the dated repository record is provenance for this
  implementation, not evidence of priority. Prior-art coverage and equivalence
  to published constructions remain OPEN. No priority or novelty claim is made.

---

## 7. The stance-objects (what to take with you)

If one thing survives this document, let it be the objects, not the phrasing:

- **The operator:** per site I → X, X → I, Y → iZ, Z → iY; unitary; order 4;
  Π² = conjugation by X^⊗N; a strict site-local product.
- **The identity:** Π·L·Π⁻¹ = −L − 2Σγ·I, exact on every tested graph;
  equivalently {Π, L_c} = 0 for the centered Liouvillian.
- **The placement:** P-type anticommutation of the shifted generator; globally
  Π²=Ux, sectorwise P_px=√px·Π is involutive. The final irreducible-sector
  classes have not been computed.
- **The data:** exact centered spectral pairing plus raw-multiset gap-ratio summaries,
  with no standard-ensemble calibration; separately, sector-resolved ⟨cos θ⟩ → −0.089 / −0.129 / −0.162 (N = 6/7/8)
  in dense blocks; dilute blocks are Poisson-like in the measured finite-N comparison.
- **The question:** which class does each fully reduced sector occupy?

The repository's in-house name for the spectral reflection is left behind on purpose; it
was painted true at our stance (2026) and is not needed to use any object
above. They are yours to place, rename, and test in your own language; that
is what an adapter is for.
