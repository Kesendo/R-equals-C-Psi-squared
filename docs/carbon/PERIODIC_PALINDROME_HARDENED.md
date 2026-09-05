# The Periodic Palindrome, Hardened: Mostly Ramp, Mixed Beyond It

**Date:** 2026-06-28
**Authors:** Tom + Claude
**Status:** Tier 1 computational result for the stated decomposition and seeded
sign-flip null on the repository's curated input arrays. The calculation qualifies
the shuffle-null result in
[`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py).
It does not evaluate the formal F1 Liouvillian identity or establish an atomic,
material, or causal correspondence; no degree of freedom, Hamiltonian, bath,
preparation, measurement, or material producer is selected here. It complements
[The Periodic Palindrome and the V-Effect](PERIODIC_PALINDROME_VS_V_EFFECT.md).
**Script:** [`simulations/periodic_palindrome_gate.py`](../../simulations/periodic_palindrome_gate.py)

---

## The question

`periodic_palindrome.py` stores curated numeric rows under period/property labels and
reports reflection-pair sums `v_k + v_{N+1-k}`. Its shuffle null returns small values
for its pair-sum score. **A linear ramp satisfies pair-sum constancy exactly**, while
that shuffle null rejects only an unordered row. The question here is therefore
narrow: after subtracting the fitted linear antisymmetric component, how does the
remaining finite input vector compare with a sign-flip null?

That vector calculation is not the formal F1 statement
`Π L Π⁻¹ = −L − 2Σγ I`: the latter concerns a specified Liouvillian with its
Hamiltonian and local dephasing channel. The present rows define neither such an
operator nor a physical measurement.

## The hardened test (gate-first)

For each supplied row, center the vector (`w = v − mean`) and split it about the
array centre into an antisymmetric part `a_k = (w_k − w_{N−1−k})/2` and a symmetric
part `s_k = (w_k + w_{N−1−k})/2` (pair-sum constancy ⟺ `s = 0`). Remove the fitted
linear antisymmetric component `a_lin`; the residual is `r = a_non + s`.

**The null.** `R = E_non / (E_non + E_sym)` is the antisymmetric fraction of the
post-ramp residual. A magnitude-preserving **sign-flip null** flips each array
partner independently: it preserves every `|r_k|` and swaps that pair's
antisymmetric and symmetric energies. Its mean is exactly 0.5; the script's
seeded 40,000-draw distribution supplies the reported one-sided p-values. Its
synthetic controls give no residual for a pure ramp, `R > 0.5` for an
antisymmetric bump, and `R < 0.5` for a symmetric bump. These are properties of
the calculation, not a physical observer or mechanism.

## Result

| Curated row | reflection pairs | R (antisymmetric residual fraction) | seeded sign-flip readout |
|-------------|------------------|-------------------------------------|--------------------------|
| 2 (Li–Ne) | 4 | 0.380 | symmetric-residual direction, p = 0.496 |
| 3 (Na–Ar) | 4 | 0.247 | symmetric-residual direction, p = 0.249 |
| 4 (K–Kr) | 9 | 0.588 | antisymmetric-residual direction, p = 0.219 |
| 5 (Rb–Xe) | 9 | 0.729 | antisymmetric-residual direction, p = 0.045 |
| 6 (Cs–Rn) | 16 | 0.583 | antisymmetric-residual direction, p = 0.120 |
| **pool rows 2,3** | 8 | **0.327** | symmetric-residual direction, p = 0.199 |
| **pool rows 4,5,6** | 34 | **0.622** | antisymmetric-residual direction, p = 0.010 |

For the five supplied IE rows, the fitted-ramp fraction
`1 − E_resid / ||w||²` is 93 %, 90 %, 68 %, 67 %, and 70 %. The two supplied
Allen-labelled rows give post-ramp residual energies 0.017 and 0.007 in squared
input units. These are reproducible properties of the embedded arrays.

The direct readout is limited:

1. The original shuffle score cannot separate the pair-sum regularity of these
   inputs from their fitted-ramp content.
2. Pooling the supplied rows 2 and 3 gives `R = 0.327`, `p = 0.199`; that is a
   symmetric-residual direction under this null, with eight pairs.
3. Pooling rows 4 through 6 gives `R = 0.622`, `p = 0.010`; that is an
   antisymmetric-residual direction under this null, with 34 pairs.

The producer does not assign either direction to shell filling, a band model, an
atomic mechanism, or a material process. Those translations remain open until a
separately specified model and producer supply them.

## Verdict

For these curated arrays, a small shuffle-null score is not evidence for the formal
F1 identity beyond fitted-ramp content. Conversely, the post-ramp `R` values are
not counterexamples to F1: they are not Liouvillian spectra or operator residuals.
They neither validate nor refute the formal result.

Any physical bridge would have to name a degree of freedom, Hamiltonian and coupling
convention, bath/jump channel, preparation, measurement operator, and a producer
that calculates the stated observable. This page supplies none of those inputs and
makes no periodic-table, atomic, material-observer, or causal claim.

## Anchor

- Script: [`simulations/periodic_palindrome_gate.py`](../../simulations/periodic_palindrome_gate.py),
  built on [`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py).
- Complements [The Periodic Palindrome and the V-Effect](PERIODIC_PALINDROME_VS_V_EFFECT.md),
  which runs a distinct per-pair test on the same curated-input genre.
- Parent: [README.md](README.md).
