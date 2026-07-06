# A Shifted, Order-4 Chiral Symmetry in Local-Dephasing Lindbladians

<!-- Keywords: Lindbladian symmetry classification gap, shifted chiral symmetry order
four, many-body Lindbladian tenfold way beyond, anti-pseudo-Hermiticity constant
shift, integrable chiral Lindbladian Poisson statistics, complex spacing ratio
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
stance, this carries a result of ours out to a community's stance. It is
written for one audience, the **symmetry-classification community for open
quantum systems** in the line of Sá, Ribeiro, and Prosen ("Symmetry
Classification of Many-Body Lindbladians: Tenfold Way and Beyond", Phys. Rev.
X 13, 031019; [arXiv:2212.00474](https://arxiv.org/abs/2212.00474)), and it
concerns an object living in one of that scheme's own margins.

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

## 1. Your open edges, in your words

The tenfold-way-and-beyond program classifies many-body Lindbladians by
their behavior under the generator set {T₊, C₊, T₋, C₋, P, Q₊, Q₋}:
antiunitary time-reversal and particle-hole flavors, unitary chiral P, and
the two pseudo-Hermiticity conditions involving L†. Two edges of that
program are live:

1. **Symmetry conditions the generator set does not span.** The scheme
   classifies conditions of the form S·L·S⁻¹ = ±L or ±L†, without constant
   shifts and with involutive chiral operators (P² = I). Physical
   Lindbladians that satisfy a symmetry condition *with* a shift, or with a
   chiral operator of higher order, sit outside the net by construction.

2. **The class-statistics correspondence when integrability intervenes.**
   The classification's diagnostic arm predicts random-matrix statistics
   (per class) within each symmetry sector. An integrable Lindbladian can
   carry the full symmetry structure of a class while showing none of its
   level repulsion, and where the crossover to the predicted statistics
   happens in physical models is a data question.

This document hands over one concrete, exactly-solvable model family that
sits on both edges at once, with the operator in closed form and the
statistics computed.

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
- **Order 4, not 2:** Π⁴ = I, and Π² is conjugation by the global string
  X^⊗N (it acts on a Pauli string as (−1)^{n_Y + n_Z}). A chiral-type
  operator with P⁴ = I where the scheme's P requires P² = I.
- **Site-local:** a tensor product of single-site operations, so it survives
  arbitrary coupling graphs and inhomogeneous rates γ_i unchanged.
- **Chiral for both:** after centering, Π anti-commutes with L_c and with
  L_c† ([PT-Symmetry Analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md),
  residuals 0.0).

---

## 3. Where it sits in your scheme

We placed Π against the generator set ourselves
([KMS and Detailed Balance](../KMS_DETAILED_BALANCE.md), the repository's
formal analysis):

- The **closest slot is Q₋** (anti-pseudo-Hermiticity, Q·L·Q⁻¹ = −L†), but
  the condition above has **L, not L†**, on the right side, and L ≠ L† as
  soon as H ≠ 0.
- The **shift 2Σγ has no slot**: the scheme's symmetry conditions carry no
  constant term. Centering absorbs the shift, but then the classification
  applies to L_c, not to L.
- After centering, the honest reading is a **generalized P**: a unitary
  chiral operator of order 4 instead of 2, anti-commuting with both L_c and
  L_c†, with all of P's spectral consequences (the exact ± pairing of the
  centered spectrum) intact.

Two placement notes we carry along honestly. First, the *shape* of the
shifted spectrum (pairs at ±λ plus a constant) has a catalogued
noninteracting home: the "shifted sublattice symmetry" of Kawasaki,
Mochizuki, and Obuse (Phys. Rev. B 106, 035408;
[arXiv:2201.09283](https://arxiv.org/abs/2201.09283)), for quadratic
systems. What we did not find catalogued is the **interacting** case with a
strictly **site-local product** operator on arbitrary graphs, which is what
Π is. Second, as we read the tenfold-way paper, its own local-dephasing
examples note the resulting chiral symmetry of the spectrum; the closed-form
site-local operator that generates it, its order-4 structure, and the shift
bookkeeping are what this document adds (this reading of the paper should be
re-verified against it before anything is claimed in public).

So the sharp question we hand over, in your terms: **is "unitary chiral of
order 4, with constant shift" a natural extension slot of the
classification, or is centering plus generalized-P the canonical resolution?**
Either answer organizes a physical family your scheme already uses as an
example class.

---

## 4. The statistics twist: the symmetry of a class, the statistics of none

The classification's diagnostic arm assumes random-matrix statistics within
symmetry sectors. This family declines, in an instructive way
([Random Matrix Theory](../../experiments/RANDOM_MATRIX_THEORY.md)):

- **Symmetry axis:** the centered spectrum is exactly ± paired (machine
  precision, ~1e-15) at every tested N: the structure of class AIII
  (chiral unitary).
- **Statistics axis:** the level statistics are **Poisson**, not chiral-RMT:
  the spacing ratio converges to ⟨r⟩ = 0.36 to 0.39 with N (Poisson 0.386,
  GOE 0.536, GUE 0.603; 21,840 eigenvalues total, N = 2 to 7 chain), and the
  complex spacing ratio of the chain reads clean 2D-Poisson. An **integrable
  chiral Lindbladian**: the symmetry of a class, the statistics of none of
  the 38.

Sector-resolved, the picture sharpens into a data point for your
class-statistics correspondence
([Filling Threshold Chaos](../../experiments/FILLING_THRESHOLD_CHAOS.md)):

- **Dilute coherence sectors stay Poisson at every tested N**, and stay
  Poisson under integrability breaking (XXZ anisotropy Δ, added disorder):
  the null is structural (too few interacting excitations to thermalize),
  not fine-tuned.
- **Dense (near-half-filling) sectors climb toward GinUE with N**: the
  complex-spacing-ratio marker ⟨cos θ⟩ runs −0.089 / −0.129 / −0.162 at
  N = 6 / 7 / 8 (GinUE −0.241, 2D-Poisson 0), with class A the licensed
  comparison inside such a block (Π maps the block to its conjugate partner,
  so no residual antiunitary constrains it).

In one sentence: in this family, **dissipative chaos is a filling
threshold**, not an integrability-breaking threshold, and the crossover is
visible in exactly the diagnostic your community uses.

---

## 5. What you can check with your own tools

Three independent handovers, each cheap on your side:

> 1. **The identity, in one line of code.** Build Π from the per-site rule
>    above (a 4^N signed permutation), build L for any small Heisenberg +
>    local-dephasing model of your choice, and check
>    Π·L·Π⁻¹ + L + 2Σγ·I = 0 to machine precision. No fitting, no limits;
>    it either vanishes or it does not.
> 2. **The classification question.** Decide where "unitary chiral, order 4,
>    constant shift" belongs: a new slot, or generalized-P after centering.
>    The operator, its algebra (Π² = X^⊗N conjugation, Π⁴ = I), and our own
>    placement walk are laid out in the sources; the family is analytically
>    tame enough to be a clean test case.
> 3. **The crossover, with your diagnostic.** Reproduce the sector-resolved
>    complex spacing ratio: the prediction is GinUE convergence in the dense
>    blocks as N grows and a structurally pinned Poisson in the dilute
>    blocks, robust against Δ and disorder.

---

## 6. What we are not claiming

Weak coupling, deliberately: this adapter offers a test case and a sharp
question to a classification program, not a rival scheme.

- **We placed Π against the Sá-Ribeiro-Prosen scheme only.** We have not
  worked the Bernard-LeClair or the Kawabata et al. non-Hermitian schemes;
  any placement there is yours to make, and may well close the gap.
- **"Outside all 38" is a statistics statement, not an algebra complaint:**
  the classes assume random-matrix statistics within sectors; this family is
  integrable and declines that assumption while carrying the AIII symmetry
  structure. Nothing in the scheme is wrong about it.
- **The statistics are computational and small-N** (N ≤ 7 global, N ≤ 8
  sector-resolved, with bootstrap intervals at the largest N). The identity
  itself is exact and proven; the statistics are data.
- **Z-dephasing only** (single-axis, unital). Other dissipators have
  analogous single-axis operators in the repository, but the identity as
  stated is for the dephasing family.
- **The Hamiltonian-level sublattice BDI is not the result** (the altitude
  note of the opening section): that symmetry is textbook and lives on a
  different space; conflating the two levels would overstate what is new.
- **Priority bookkeeping:** the identity is dated 2026-03-14 in the
  repository's history. The nearest published neighbour we track
  ([arXiv:2605.20930](https://arxiv.org/abs/2605.20930), May 2026) works in
  the same model family and does not state the operator or the spectral
  reflection.

---

## 7. The stance-objects (what to take with you)

If one thing survives this document, let it be the objects, not the phrasing:

- **The operator:** per site I → X, X → I, Y → iZ, Z → iY; unitary; order 4;
  Π² = conjugation by X^⊗N; a strict site-local product.
- **The identity:** Π·L·Π⁻¹ = −L − 2Σγ·I, exact on every tested graph;
  equivalently {Π, L_c} = 0 for the centered Liouvillian.
- **The placement:** nearest slot Q₋ (fails: L, not L†); shift unslotted;
  generalized-P (order 4) after centering; noninteracting shape catalogued
  (shifted sublattice symmetry), interacting site-local case not.
- **The data:** AIII structure + Poisson statistics (⟨r⟩ = 0.36 to 0.39,
  N ≤ 7); sector-resolved ⟨cos θ⟩ → −0.089 / −0.129 / −0.162 (N = 6/7/8)
  in dense blocks; dilute blocks structurally Poisson.
- **The question:** new slot, or generalized-P after centering?

Our in-house name for the spectral reflection is left behind on purpose; it
was painted true at our stance (2026) and is not needed to use any object
above. They are yours to place, rename, and test in your own language; that
is what an adapter is for.
