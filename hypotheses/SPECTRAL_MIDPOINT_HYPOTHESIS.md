<!-- QUARTER-CURRENT -->
# Spectral-midpoint question after the basis-census check

Current reading: the retained N=3 and N=5 tables form a basis-dependent coefficient census
in a fixed, non-orthogonal eigenvector normalization. They are not spectral
projector weights.  A centered eigenvalue does not imply that an individual
eigenvector is Pi-fixed; an invariant two-perspective test remains open.

<!-- QUARTER-HISTORICAL -->
**Historical reading:** the original midpoint hypothesis and constructed
SLOW/FAST perspective are retained below as research context, not as a measured
Pi-side observation.

# Spectral Midpoint Hypothesis: Both Sides See the Center

<!-- CROSSING-CURRENT -->

<!-- Keywords: spectral midpoint palindromic Liouvillian, CΨ quarter boundary spectral decomposition,
fold catastrophe eigenvalue center, geometric mean palindromic perspectives,
dual perspective spectral analysis, R=CPsi2 spectral midpoint confirmation,
OOP polymorphism quantum palindrome analogy, AM-GM inequality spectral bands,
palindromic mirror both sides see center, quantum classical boundary spectral resonance -->

**Status:** Historical hypothesis; finite basis-dependent census, invariant test still open
**Tier:** 5 (analogy, not derivable; see the postscript)
**Date:** March 25, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Temporal Sacrifice (fold observation)](../experiments/TEMPORAL_SACRIFICE.md), [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md), [Crossing Taxonomy](../experiments/CROSSING_TAXONOMY.md)

---

## What this is about

Is there a connection between two independently proven structures in
this project?

- The scalar value **CΨ = ¼**, which is the discriminant zero of the
  separate recurrence R = C(Ψ+R)² after its variables are identified
- The **palindromic midpoint** Σγ in spectral space: the center of
  the Liouvillian's mirror-symmetric decay spectrum (every rate d has
  a partner at 2Σγ - d)

This hypothesis arose from the [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md)
experiment, where one finite grid put an endpoint-mutual-information maximum
near an endpoint CΨ crossing. The question was whether an invariant relation
exists; this document's coefficient census does not establish one.

---

## The Problem in Plain Language

Imagine two people standing on opposite sides of a glass wall. Each
can see through the glass, but what looks like "left" to one looks
like "right" to the other. If you ask just one person what they see,
you get a biased picture. Only by combining both perspectives do you
see the wall itself.

The palindromic spectrum works the same way. It pairs every decay
mode with a mirror partner:

```
Fast mode (rate d)  <-->  Slow mode (rate 2Σγ - d)
```

The midpoint rate Σγ is the glass in the analogy: its scalar rate coordinate
is fixed by d → 2Σγ-d. A centered eigenvalue does not imply that its
eigenvector is fixed by Π, especially in a degenerate eigenspace.

We asked: at the sampled time where this CΨ readout equals ¼, does the
midpoint bin dominate a chosen coefficient census?

From one side alone: no. The spectral weight is spread across bands.
Under a constructed SLOW↔FAST relabelling: the displayed geometric mean makes
the MID bin largest for these two tables. That is not a second observation.

---

## The historical conjecture

The conjecture asked whether midpoint modes dominate at t\* where CΨ(t\*)=¼
under an invariant two-perspective measure. The calculation below instead uses
basis- and normalization-dependent coefficient magnitudes, so that conjecture
remains open.

The historical proposal used the **geometric mean** of a coefficient census and
its constructed rate-bin reversal:

```
w_combined(band) = √( w_our_side(band) × w_Π_side(band) )
```

Why the geometric mean? Because it is the simplest combination that
treats both perspectives equally and preserves multiplicative structure.
An arithmetic mean would mask the asymmetry. A geometric mean exposes it.

---

## Why One Side Is Not Enough

### Step 1: Eigendecomposition (what we computed)

We built the full Liouvillian matrix for N=3 (64×64) and N=5 (1024×1024),
decomposed the initial state |+⟩ᴺ in a numerically normalized eigenbasis, and
tracked coefficient-magnitude sums over time. These are not invariant spectral
projector weights. Three bins were classified by distance
from the palindromic midpoint Σγ:

- **SLOW:** modes with decay rate d < Σγ - γ (slow from our perspective)
- **MID:** modes with |d - Σγ| < γ (near the midpoint, the "glass")
- **FAST:** modes with d > Σγ + γ (fast from our perspective)

(The rate≈0 bin was excluded from the percentages; calling it a classical
floor was an interpretation, not an eigendecomposition result.)

At the CΨ = ¼ crossing, from our side alone:

| | SLOW | MID | FAST |
|-----|------|-----|------|
| N=3 | 13% | 53% | 34% |
| N=5 | 45% | 47% | 8% |

MID is the largest band, but not dominant. At N=5, SLOW nearly ties
MID (45% vs 47%). From one perspective alone, there is no clear
midpoint concentration. **The simple hypothesis appears falsified.**

### Step 2: The product problem (why it seemed broken)

CΨ = Tr(ρ²) × L₁/(d-1) is a product of purity (diagonal in the
Pauli basis) and coherence (off-diagonal in the computational basis).
These live in different bases. The cross-terms between them mix all
spectral bands, preventing clean midpoint concentration in any
single-basis decomposition.

This is real. The product structure is why one perspective cannot
see the midpoint clearly. But the solution is not to fix the math.
It is to **add the missing perspective**.

---

## Why Both Sides Together See the Midpoint

### Step 3: a constructed rate-bin reversal swaps the edges

The Π conjugation (the proven operator from the
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md))
maps each decay rate to its palindromic partner:

```
d  -->  2Σγ - d
```

A spectral partner has the complementary rate 2Σγ-d. The table below was
constructed by swapping the aggregate SLOW and FAST numbers; it did not apply
Π to the state and remeasure coefficients. At the scalar-bin level:
- What is SLOW from our side (d small) is FAST from the Π side (2Σγ - d large)
- What is FAST from our side is SLOW from the Π side
- the MID rate bin maps to the MID rate bin because d ≈ Σγ implies 2Σγ-d ≈ Σγ

Thus the historical construction assigned a second column with **SLOW and FAST swapped**:

| | SLOW (Π) | MID (Π) | FAST (Π) |
|-----|----------|---------|----------|
| N=3 | 34% | 53% | 13% |
| N=5 | 8% | 47% | 45% |

Each side sees a lopsided spectrum. Each side thinks one edge band is
large. But they disagree about WHICH edge band is large.

### Step 4: The geometric mean reveals the midpoint

How do we combine two perspectives that each see a biased picture?
The geometric mean: √(perspective_A × perspective_B). It punishes
disagreement and rewards consensus.

**N=3 at CΨ = ¼ crossing:**

| Band | Our side | Π side | Geometric mean √(ours × Π) |
|------|----------|--------|-----------------------------|
| SLOW | 13% | 34% | **21%** |
| **MID** | **53%** | **53%** | **53%** |
| FAST | 34% | 13% | **21%** |

**N=5 at CΨ = ¼ crossing:**

| Band | Our side | Π side | Geometric mean |
|------|----------|--------|----------------|
| SLOW | 45% | 8% | **19%** |
| **MID** | **47%** | **47%** | **47%** |
| FAST | 8% | 45% | **19%** |

**In this constructed coefficient summary, the MID number is largest.** At N=3:
MID=53% vs SLOW=FAST=21%. At N=5: MID=47% vs SLOW=FAST=19%.

### Why this works (the math)

The exact palindrome guarantees the rate pairing. The following arithmetic then
holds for the deliberately swapped aggregate table:

1. **SLOW and FAST labels swap under the constructed rate reversal.** This does
   not establish equality of state coefficients under a physical second view.

2. **The MID scalar rate is centered.** The arithmetic
   2Σγ-Σγ=Σγ fixes the eigenvalue coordinate. It neither fixes an
   individual eigenvector nor supplies an invariant weight.

3. **Asymmetry shrinks the edges.** At the crossing, SLOW and FAST
   are unequal (one large, one small). Their geometric mean
   √(large × small) is always less than either alone. (For any two
   positive numbers, the geometric mean is less than or equal to the
   arithmetic mean. When they are unequal, strictly less.)

In formulas:
```
MID_combined  = √(MID × MID) = MID              (unchanged)
SLOW_combined = √(SLOW_us × FAST_us)             (mixed)
FAST_combined = √(FAST_us × SLOW_us)             (same as SLOW_combined)
```

**Important caveat:** This arithmetic does not guarantee
`MID > sqrt(SLOW*FAST)` in general. It only makes the two constructed edge
numbers equal and no larger than their original maximum. In the N=3 and N=5
tables the midpoint coefficient sum is numerically larger than that constructed
geometric mean. Because the coefficients are basis- and normalization-dependent,
the comparison establishes neither spectral dominance nor a physical mechanism
at the selected time.

---

## Current disposition of the hypothesis

### The two centered scalars are not one condition

CΨ=¼ can be inserted into a recurrence whose discriminant is 1-4CΨ.
The spectral midpoint Σγ is the center of a separate rate involution. The
shared visual language of "meeting" does not make these the same condition.

- **Recurrence:** two fixed-point roots merge at the specified parameter.
- **Spectrum:** partner eigenvalue coordinates are centered at Σγ.

No map between these objects has been proved here.

### What would be needed

An invariant test would use projectors or another normalization-independent
object, specify how Π acts on the state and degenerate subspaces, and compare it
against controls away from the selected CΨ time. The present tables do none of
those things.

The facing-mirrors image remains an invitation for designing that test, not a
result extracted from the coefficient census.

### Historical prediction, not yet supported

The asymmetry between SLOW and FAST grows with N:
- N=3: SLOW/FAST ratio = 2.6×, √(SLOW × FAST) = 21% vs MID = 53%
- N=5: SLOW/FAST ratio = 5.6×, √(SLOW × FAST) = 19% vs MID = 47%

The constructed geometric mean shrinks unequal edge numbers. Extrapolating that
arithmetic to near-complete midpoint dominance at larger N is not licensed
without an invariant measure. N=7 and N=9 remain possible controls.

---

## The Failed Path (documented for honesty)

### Simple form (one-sided, falsified)

The original hypothesis predicted >80% spectral weight in the MID
band at the crossing, viewed from our side alone. This fails because
CΨ = C × Ψ is a product that mixes spectral bands via cross-terms
between two different bases (Pauli for purity, computational for
coherence).

### The attempted two-sided construction

The cross-term problem remains. The second column was produced by swapping
the aggregate SLOW and FAST numbers; it was not obtained by applying Π to a
state, eigenspace, projector, or normalized coefficient vector. The exact F1
eigenvalue pairing therefore does not certify this geometric-mean construction.

Looking from both sides remains a useful invitation for the missing invariant
test. It is not a resolved result of the present tables.

---

## Computational Details

**Tool:** `dotnet run -c Release -- spectral <N>` in
compute/RCPsiSquared.Propagate.Test

**Method:** Build Liouvillian from LindbladPropagator.EvalRHS applied
to 64 (N=3) or 1024 (N=5) basis density matrices. Eigendecompose
via MathNet.Numerics Evd(). Decompose |+⟩ᴺ in eigenbasis. Track
band weights |c_k · exp(-d_k · t)| at 0.1 time intervals.

**Band definition:** |d - Σγ| < γ for MID, d < Σγ - γ for SLOW,
d > Σγ + γ for FAST, d < 0.001 for IMMUNE. IMMUNE excluded from
dynamic analysis.

---

## Postscript: The Palindrome as Polymorphism

*For those who think in code. Not physics. Not proof.
A translation, for those who build systems for a living.*

*Tier 5 (analogy, not derivable).*

---

**Interpretive software analogy.** A palindromic eigenvalue has a partner at the
complementary rate. The two "interfaces" below are a programming metaphor, not
two physical observations of one mode.

This is polymorphism.

The Liouvillian is the base class. It defines the eigenvalue spectrum,
the pairing rule, the midpoint. Every observer inherits from it. But
each observer implements `observe()` differently: our side measures
populations and calls them "real." The Π side measures coherences
and calls them "real." Both are correct. Both are incomplete.

The older analogy cast CΨ=¼ as an abstract method. Current mathematics does not
connect that scalar readout to collapse, classicality, or the spectral rate
Σγ; those lines are story, not an interface contract of the model.

And the midpoint, the spectral center where d = Σγ, is the point
where the two implementations are identical. The mode that returns
the same value regardless of which interface you call it through.
The base class made visible. The only place where inheritance
collapses to identity.

You cannot see this from one side. A caller holding a reference to
`OurSide` sees SLOW=45%, MID=47%, FAST=8%. A caller holding
`PiSide` sees SLOW=8%, MID=47%, FAST=45%. Each thinks the spectrum
is lopsided. Only by holding both references and computing the
geometric mean does the midpoint emerge: MID=47%, edges=19%.

The single-perspective view is like reading one class in isolation.
You see the methods, the fields, the behavior. But you miss the
contract. You miss what the system *is*, because what it *is* lives
in the interface between implementations, not in any one of them.

The programmers who first struggled with OOP remember the moment
it clicked: the object is not the class. The object is the behavior
that emerges when multiple classes interact through a shared
interface. The class is just a perspective. The interface is the
reality.

The palindrome is the interface. The two perspectives are the
classes. The midpoint is where `@Override` returns `super`.

In the metaphor, CΨ=¼ is when the program runs. In the calculation, it is only
the time selected for the finite coefficient table.

---

## References

- [Temporal Sacrifice (fold observation)](../experiments/TEMPORAL_SACRIFICE.md)
- [Mirror Symmetry Proof (Π conjugation)](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
- [CΨ monotonicity: named trajectories, not a universal Hamiltonian-live envelope](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)
- [Crossing Taxonomy: fixed-readout, Hamiltonian-dead Bell+ gamma sweep in two books](../experiments/CROSSING_TAXONOMY.md)
- [Mathematical Connections (fold, Mandelbrot)](../docs/MATHEMATICAL_CONNECTIONS.md)
