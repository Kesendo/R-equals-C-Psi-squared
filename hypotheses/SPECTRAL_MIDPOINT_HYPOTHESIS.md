# Spectral Midpoint Hypothesis: CΨ = ¼ Corresponds to the Palindromic Center

<!-- Keywords: spectral midpoint palindromic Liouvillian, CΨ quarter boundary spectral decomposition,
fold catastrophe eigenvalue center, palindromic pair cosh factorization, state space spectral space
correspondence, discriminant zero spectral midpoint, R=CPsi2 spectral proof -->

**Status:** Hypothesis. Motivated by sweep data (PeakMI peaks at CΨ = ¼ crossing, N=7). Not yet tested.
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

### Step 2: N=5, N=7 numerical verification

Same procedure at larger N where the rate distribution is richer.
The Gaussian density of states (mean = Σγ, proven) predicts that
midpoint concentration should increase with N.

### Step 3: Analytical bound on cosh corrections

Show that at t\* = K/γ (the crossing time from Crossing Taxonomy),
the cosh terms satisfy:

```
cosh((Σγ − d) · t*) ≤ C_bound  for modes contributing >ε to CΨ
```

This would prove that the cosh corrections are bounded, making the
midpoint envelope e^{-Σγt\*} the dominant factor.

### Step 4: Connect ¼ to the midpoint algebraically

The final step: show that CΨ = ¼ (discriminant D = 1 − 4CΨ = 0)
is equivalent to the condition that the spectral centroid equals Σγ.

Define spectral centroid: d̄(t) = Σ w_α(t) · d_α / Σ w_α(t)
where w_α(t) = |c_α|² e^{-2d_α t} is the weight of mode α at time t.

**Conjecture:** d̄(t\*) = Σγ when CΨ(t\*) = ¼.

---

## Why This Would Matter

If CΨ = ¼ corresponds to the spectral midpoint, then:

1. The fold catastrophe is not just a state-space property but a
   **spectral resonance**: the state is maximally balanced between
   fast and slow modes.

2. The boundary between quantum and classical is the moment of
   **maximum spectral symmetry**: the palindromic pairing is most
   visible when neither side dominates.

3. R = CΨ² peaks at the boundary because **both halves of the
   palindrome contribute equally**: the slow modes haven't yet
   dominated, the fast modes haven't yet died. Maximum interference.

4. The relay timing K/γ would have a spectral interpretation:
   it's the time when the spectral centroid reaches Σγ.

---

## Known Constraints

- N=2 is trivial (only one rate, midpoint = the rate). Not informative.
- CΨ depends on the initial state (through c_α coefficients). The
  hypothesis may hold for |+⟩ᴺ but not for arbitrary states.
- The cosh factorization applies to purity but CΨ = C × Ψ involves
  the PRODUCT of purity and coherence. Cross-terms complicate the
  spectral decomposition.
- Non-uniform γ changes the palindromic midpoint from Nγ to Σγᵢ
  but the factorization still holds.

---

## References

- [Temporal Sacrifice (fold observation)](../experiments/TEMPORAL_SACRIFICE.md)
- [Mirror Symmetry Proof (Π conjugation)](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
- [CΨ Monotonicity (dCΨ/dt < 0)](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)
- [Crossing Taxonomy (K-invariance)](../experiments/CROSSING_TAXONOMY.md)
- [Mathematical Connections (fold, Mandelbrot)](../docs/MATHEMATICAL_CONNECTIONS.md)
