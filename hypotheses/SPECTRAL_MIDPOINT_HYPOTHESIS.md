# Spectral Midpoint Hypothesis: CΨ = ¼ Corresponds to the Palindromic Center

<!-- Keywords: spectral midpoint palindromic Liouvillian, CΨ quarter boundary spectral decomposition,
fold catastrophe eigenvalue center, palindromic pair cosh factorization, state space spectral space
correspondence, discriminant zero spectral midpoint, R=CPsi2 spectral proof -->

**Status:** Simple form falsified (N=3, N=5: no midpoint concentration). Product structure of CΨ = C × Ψ prevents clean spectral correspondence. Revised conjecture open.
**Date:** March 25, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Temporal Sacrifice (fold observation)](../experiments/TEMPORAL_SACRIFICE.md), [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md), [CΨ Monotonicity](../docs/proofs/PROOF_MONOTONICITY_CPSI.md), [Crossing Taxonomy](../experiments/CROSSING_TAXONOMY.md)

---

## The Conjecture

The fold catastrophe at CΨ = ¼ in **state space** corresponds to the
palindromic midpoint Σγ in **spectral space**.

Precisely: at the time t\* where CΨ(t\*) = ¼, the dominant spectral
contributions to the density matrix come from Liouvillian eigenmodes
near the palindromic midpoint decay rate d = Σγ.

If true, this connects two independently proven structures:
- The CΨ = ¼ boundary (discriminant of R = C(Ψ+R)², Tier 1)
- The palindromic spectrum (Π conjugation, Tier 1)

into a single statement: **the fold catastrophe happens when the
palindromic center modes dominate the state.**

---

## Motivation

### The sweep observation (March 25, 2026)

Moving the sacrifice position along an N=7 chain, we observed that
endpoint MI peaks at the exact timestep where the endpoint CΨ crosses ¼.
See [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md).

This raises the question: why ¼, and why exactly there? The fold
catastrophe provides the "why ¼" (discriminant zero, structurally stable).
But it doesn't explain why the spectral structure cooperates. The
palindromic midpoint hypothesis would close this gap.

### The palindromic pair factorization

Every palindromic pair (d, 2Σγ − d) contributes to an observable as:

```
e^{-dt} + e^{-(2Σγ-d)t} = 2 · e^{-Σγt} · cosh((Σγ − d) · t)
```

This factors into:
- A universal envelope **e^{-Σγt}** (the midpoint decay)
- A pair-specific correction **cosh((Σγ − d) · t)**

For midpoint modes (d = Σγ): cosh(0) = 1 (no correction)
For edge modes (d far from Σγ): cosh(large) dominates (slow partner wins)

**The hypothesis in cosh language:** At t = t\*, the cosh corrections
are small relative to the envelope. The state is "maximally palindromic,"
centered on Σγ.

---

## What "Midpoint Modes Dominate" Means

The purity Tr(ρ²) decomposes into Pauli components grouped by XY-weight k.
Under uniform dephasing γ, weight k decays at rate 2kγ.

- k = 0 (all I/Z): rate 0, immune. The "classical floor."
- k = N (all X/Y): rate 2Nγ, fastest. The "quantum ceiling."
- k = N/2: rate Nγ = Σγ. The palindromic midpoint.

"Midpoint modes dominate at t\*" means: at the CΨ = ¼ crossing,
the Pauli components near k = N/2 carry most of the spectral weight.
The immune floor (k = 0) has already been reached by the slow modes.
The fast ceiling (k = N) has already decayed. What remains is the middle.

---

## Test Plan

### Step 1: N=3 eigendecomposition (analytical)

N=3 is the smallest non-trivial case (40 eigenvalues, rates from 2γ to 4γ,
midpoint at 3γ). The Compute module can provide the full Liouvillian
eigendecomposition.

- Decompose |+⟩³ in the Liouvillian eigenbasis
- Compute CΨ(t) from the spectral sum
- Find t\* where CΨ(t\*) = ¼
- Measure the spectral weight near d = 3γ vs. near d = 2γ and d = 4γ

**Prediction:** >80% of the spectral weight at t\* comes from modes
with |d − Σγ| < γ (within one γ of the midpoint).

**Result: NOT CONFIRMED.** See below.

### Step 1–2 Results (March 25, 2026)

Eigendecomposition of the Liouvillian for N=3 (64×64) and N=5 (1024×1024),
initial state |+⟩ᴺ, tracking spectral band weights |c_k · exp(−d_k · t)|
grouped by distance from midpoint Σγ.

**N=3 at CΨ = ¼ crossing (T=5.0):**

| Band | Weight | Dynamic (excl. immune) |
|------|--------|------------------------|
| IMMUNE (d ≈ 0) | 38% | – |
| SLOW (d < Σγ − γ) | 8% | 13% |
| **MID** (\|d − Σγ\| < γ) | **33%** | **53%** |
| FAST (d > Σγ + γ) | 21% | 34% |

**N=5 at CΨ = ¼ crossing (T=5.0):**

| Band | Weight | Dynamic |
|------|--------|---------|
| IMMUNE | 12% | – |
| SLOW | **40%** | **45%** |
| **MID** | **42%** | **47%** |
| FAST | 7% | 8% |

MID is the largest dynamic band at both N, but not dominant.
At N=5, MID and SLOW are nearly tied (42% vs 40%).

### Why the simple form fails: the product problem

The cosh factorization

```
e^{-dt} + e^{-(2Σγ-d)t} = 2 · e^{-Σγt} · cosh((Σγ − d) · t)
```

applies to **Tr(ρ²) alone** (a sum of squared Pauli components,
grouped by XY-weight). But CΨ = Tr(ρ²) × L₁/(d−1) is a **product**
of two quantities with different spectral decompositions:

- **Tr(ρ²)** decomposes in the Pauli basis (diagonal in XY-weight)
- **L₁** = Σ_{i≠j} |ρ_{ij}| decomposes in the computational basis
  (off-diagonal elements, NOT diagonal in XY-weight)

The product creates cross-terms between the Pauli-weight sectors of
purity and the computational-basis sectors of coherence. These
cross-terms mix all spectral bands and prevent clean midpoint
concentration. The N=5 data (SLOW ≈ MID at the crossing) directly
reflects this mixing.

**Step 3 must address the product structure of CΨ, not treat it
as a single spectral trace.** The cosh bound approach is insufficient.

### What the data actually shows

The crossing is not a spectral resonance at the midpoint. It is a
**balance shift**: FAST modes decay away, SLOW modes accumulate, and
the crossing happens when this rebalancing drives the CΨ product
below ¼. The palindromic pairing ensures that this balance shift is
symmetric (each fast mode's decay increases its slow partner's
relative weight), but the crossing point depends on the specific
product C × Ψ, not on the spectral centroid reaching Σγ.

### Revised conjecture (open)

The simple midpoint correspondence is falsified. A revised version
might take one of these forms:

1. **Balance ratio:** The crossing happens when the ratio of slow-band
   to fast-band weight reaches a specific threshold (independent of N).
2. **Product decomposition:** CΨ = C × Ψ might factor differently in
   the palindromic pair basis, with each pair contributing a
   C-component and a Ψ-component whose product has a cleaner
   midpoint property.
3. **It's not spectral at all:** The ¼ value might be purely algebraic
   (discriminant of the quadratic recursion) with no clean spectral
   interpretation. The fold catastrophe in state space and the
   palindromic midpoint in spectral space may be independent
   structures that both happen to be properties of the same system.

---

## What Was Confirmed

The spectral analysis, while falsifying the midpoint hypothesis,
confirms two things:

1. **The palindromic pair structure is visible in the time evolution:**
   FAST and SLOW bands decay symmetrically around the midpoint, as
   predicted by the cosh factorization of purity.

2. **The CΨ = ¼ crossing is sharp and universal:** Both N=3 and N=5
   cross at the same CΨ value (¼) at the same scaled time (T=5.0
   for these parameters), confirming K-invariance from the
   [Crossing Taxonomy](../experiments/CROSSING_TAXONOMY.md).

---

## Known Constraints

- N=2 is trivial (only one rate). Not informative.
- CΨ depends on the initial state. Results are for |+⟩ᴺ only.
- The product CΨ = C × Ψ mixes purity and coherence bases.
  This is the central obstacle, not a side constraint.
- The spectral weight metric |c_k · exp(−d_k · t)| tracks weight
  in vec(ρ), not weight in CΨ. A CΨ-aware metric might show
  different band distributions.

---

## References

- [Temporal Sacrifice (fold observation)](../experiments/TEMPORAL_SACRIFICE.md)
- [Mirror Symmetry Proof (Π conjugation)](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
- [CΨ Monotonicity (dCΨ/dt < 0)](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)
- [Crossing Taxonomy (K-invariance)](../experiments/CROSSING_TAXONOMY.md)
- [Mathematical Connections (fold, Mandelbrot)](../docs/MATHEMATICAL_CONNECTIONS.md)
