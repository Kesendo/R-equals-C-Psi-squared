# τ_max vs the Spectral Gap: γ Is the Timekeeper

<!-- Keywords: tau_max spectral gap formula rejected, gamma is the timekeeper, coherence time
1/(2gamma), inverse spectral gap clock, palindromic Liouvillian 2gamma floor, decay rate
J-independent, gamma=0 clock stops, R=CPsi2 gamma time, auto-extracted formula consistency test -->

**Status:** Rejected: the formula is wrong; the clock is τ = 1/(2γ)
**Date:** May 28, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** [tau_max_spectral_gap.py](../simulations/tau_max_spectral_gap.py),
[tau_max_spectral_gap_verify.py](../simulations/tau_max_spectral_gap_verify.py)
**Raw data:** [tau_max_spectral_gap.txt](../simulations/results/tau_max_spectral_gap.txt)

---

## Abstract

An auto-extracted, fragmentary formula,

    τ_max = ħ / √(λ₂(L) · J²)

surfaced from `D:\new_formulas.txt`, framed as a "universal clock ticked by topology" with
claimed NMR / ion-trap validation that **appears nowhere in this repo**. Tested against the
repo's proven facts about the palindromic Liouvillian (Heisenberg chain + local Z-dephasing),
it is **rejected**. The relaxation timescale is the **inverse spectral gap**,

    τ = 1 / λ₂ = 1 / (2γ),

set by **γ alone; γ is the timekeeper**. At γ = 0 the gap is exactly 0 and the clock stops
(τ → ∞: no dephasing, only Hamiltonian oscillation). The coupling J sets oscillation
*frequencies*, never decay *rates*. The formula is wrong twice over: it puts the gap under a
square root (1/√λ₂ where the timescale is 1/λ₂ → γ-power −½ instead of −1), and it injects a
spurious 1/J (J-power −1 instead of 0). Both errors are exhibited as clean fitted exponents.

---

## Claim under test

τ_max = ħ/√(λ₂(L)·J²), with λ₂(L) the Liouvillian spectral gap (smallest nonzero decay rate)
and J the Heisenberg coupling. Since √(λ₂·J²) = √λ₂·J, the formula asserts **τ_max ∝ 1/J** and
**τ_max ∝ 1/√λ₂**. With the palindrome floor λ₂ = 2γ (proven N-independent) it also implies
**τ_max = const(N)**.

---

## Dimensional pre-check: the square root is the tell

In natural units (ħ = 1), J and λ₂ are both rates [1/t]:

    √(λ₂·J²) = √( [1/t]·[1/t²] ) = t^(−3/2)   ⟹   τ_max ~ t^(+3/2), not a time.

A relaxation time is **1 / rate**, i.e. **1/λ₂**, not 1/√λ₂. The square root in the formula is
already dimensionally wrong; the experiment below confirms the consequence quantitatively.

---

## What the repo already proves (re-confirmed here, all checks pass)

- **Palindrome floor:** above an N-dependent coupling threshold Q*_gap(N) in the coupling ratio Q = J/γ the smallest
  nonzero Liouvillian decay rate is exactly **2γ**, so λ₂(L) = 2γ. Below it the gap is
  Zeno-suppressed and smaller; see
  [the Absorption Theorem proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §4.3.
  Every number in this document sits above the threshold.
- **Gap is N-independent above Q*_gap(N):** λ₂ = 2γ for every N (verified N = 2…5). Q*_gap(N) itself
  grows with N: 0.50, 0.80, 1.34, 1.82 at N = 2…5 on the Heisenberg chain.
- **Gap is J-independent above the coupling threshold:** verified gap = 2γ for J ∈ 0.5…4 at
  γ = 0.05, i.e. Q = J/γ ∈ 10…80, one to two orders above the threshold Q*_gap(N) ≲ 1.9. The sweep
  never crossed the boundary, so it could not see it. Below Q*_gap(N) the gap is Zeno-suppressed
  and J-dependent; see [the Absorption Theorem proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)
  §4.3. The τ_max rejection below is unaffected: it argues from the rate scale, which is γ.
- **Gap is γ-linear:** gap = 2γ for γ ∈ 0.02…0.20.

The only rate scale in the problem is γ. The `ħ/√(λ₂·J²)` form appears nowhere in the repo.

**This positive law is the Absorption Theorem / D6.** τ = 1/λ₂ = 1/(2γ) is the slowest-mode
(⟨popcount(i⊕j)⟩ = 1) reading of the **Absorption Theorem** `Re(λ) = −2γ⟨n_XY⟩`
([`PROOF_ABSORPTION_THEOREM.md`](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)); the gap = 2γ and
the mixing time are **D6** in [`ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md). What this
experiment adds is the *rejection* of the extracted formula, plus an H-independence check that
closes the Absorption Theorem's former complex-Hamiltonian caveat.

---

## γ is the timekeeper

The relaxation clock is the inverse gap, τ = 1/λ₂ = 1/(2γ). Three facts pin it down:

1. **At γ = 0 the clock stops.** With no dephasing the Liouvillian is anti-Hermitian; every
   eigenvalue is purely imaginary, the gap is exactly 0, and τ = ∞. There is no decay clock:
   only Hamiltonian oscillation at frequencies set by J. γ is what makes time *tick*.
   (This is the repo's own thesis: [γ = source of experienced time](../docs/GAMMA_TIME_DISTINCTION.md).)
2. **J is frequency, not rate, in this regime.** Above Q*_gap(N) the gap does not move with J, so
   no relaxation timescale set by the gap can depend on J there, and the formula's 1/J is
   spurious. (Below Q*_gap(N) the spectral gap does depend on J, approaching
   2(1−cos(π/N))·2J²/γ as J/γ → 0, the prefactor being N-dependent and equal to 1 only at
   N=3. That regime does not rescue τ_max, since the dependence has the wrong sign and
   power. It also does not touch the clock: the modes that fall below 2γ are population
   modes, and every coherence keeps its 2γ floor at any coupling, so τ = 1/(2γ) as a
   coherence lifetime stands regardless.)
3. **The scale is purely γ.** τ = 1/(2γ): halve γ, double the clock.

---

## Results

Sweeps over N, J, γ comparing the true clock (1/λ₂) to the formula
([raw data](../simulations/results/tau_max_spectral_gap.txt)):

```
 sweep  N     J  gamma  lambda2  clock=1/gap   formula
     N  2  1.00  0.050   0.1000       10.000     3.162
     N  3  1.00  0.050   0.1000       10.000     3.162
     N  4  1.00  0.050   0.1000       10.000     3.162
     N  5  1.00  0.050   0.1000       10.000     3.162
     J  3  0.50  0.050   0.1000       10.000     6.325
     J  3  1.00  0.050   0.1000       10.000     3.162
     J  3  2.00  0.050   0.1000       10.000     1.581
     J  3  4.00  0.050   0.1000       10.000     0.791
 gamma  3  1.00  0.020   0.0400       25.000     5.000
 gamma  3  1.00  0.050   0.1000       10.000     3.162
 gamma  3  1.00  0.100   0.2000        5.000     2.236
 gamma  3  1.00  0.200   0.4000        2.500     1.581
```

Fitted power-law exponents (log-log slope):

| quantity | vs J | vs γ | reading |
|----------|------|------|---------|
| **clock** = 1/λ₂ | **0.000** | **−1.000** | depends only on γ; τ = 1/(2γ) |
| formula = ħ/√(λ₂J²) | −1.000 | −0.500 | spurious 1/J; wrong γ-power (−½ not −1) |

γ = 0: spectral gap = 0.0000 → clock stops (τ = ∞).

**VERDICT: REJECTED.** The formula fails on **both** required exponents: it must have J-power 0
and γ-power −1 to match the clock; it has −1 and −½.

---

## Verdict

τ_max = ħ/√(λ₂·J²) is **rejected as a law**. The relaxation timescale is the inverse spectral
gap, **τ = 1/λ₂ = 1/(2γ)** above Q*_gap(N), set by γ alone there, N-independent, J-independent,
and divergent at γ = 0 (the clock stops without dephasing). The formula's two errors are quantitative and clean:

- **Wrong functional form:** 1/√λ₂ gives a γ-power of −½; the true timescale 1/λ₂ has γ-power −1.
- **Spurious coupling dependence:** a J-power of −1, where the true clock is J-independent
  (J sets oscillation frequency, never the decay rate).

The only correct kernel one might salvage, "the gap matters", is true but trivial here, since
λ₂ = 2γ; stripped of the square root and the J, the statement is simply **the clock is γ**.

---

## Note on method (why no coherence-envelope here)

An earlier version of this experiment operationalized τ_max as the 1/e decay time of a
state's coherence (an RK4 envelope), with a per-axis classifier. That approach was **dropped**:
the L₁-coherence observable is non-monotonic (the Heisenberg coupling delocalizes the
excitation, so total coherence *grows* before decaying), which manufactured a spurious
"J-independence", and a purity-based 1/e time has a system-dependent floor that the threshold
cannot always reach. The gap-based result needs **no** such observable: τ = 1/λ₂ follows
directly from the proven spectrum, with no metric to bias it. Avoiding that fragility is the
point.

---

## What this does and does not claim

- **Does** claim: the formula is rejected (wrong γ-power *and* spurious J); the relaxation
  clock is τ = 1/(2γ), set by γ; at γ = 0 the clock stops; J is frequency, not rate.
- **Does not** claim: any reproduction of the (unrecoverable) historical "τ_max ~ N" figures;
  a quantum-speed-limit interpretation (not tested here).

---

## Reproduce

```bash
python simulations/tau_max_spectral_gap.py          # sweep + power-law fits + verdict + results file
python simulations/tau_max_spectral_gap_verify.py   # asserts the invariants and the gamma-clock
```

The verify script confirms, at γ = 0.05 and J ∈ 0.5…4, i.e. Q = 10…80, one to two orders
above every Q*_gap(N): gap = 2γ exactly; gap J-independent over that range and γ-linear
(γ 0.02…0.20); gap = 0 at γ = 0 (clock stops); clock = 1/(2γ) is N- and J-independent;
clock γ-power = −1 while the formula's is −½ and its J-power is −1; and the palindrome pairing
rᵢ + r₍ₙ₋₁₋ᵢ₎ = 2Nγ holds across the full rate multiset.

---

## References

- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md): `Re(λ) = −2γ⟨n_XY⟩`, the master rate identity; τ = 1/(2γ) is its slowest-mode case. (Generality to complex Hermitian H closed here.)
- [Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md): **D6** (spectral gap = 2γ, mixing time); F3 / F8 / F74 are decay-rate corollaries of the Absorption Theorem.
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the 2γ floor / palindromic spectrum.
- [γ–Time Distinction](../docs/GAMMA_TIME_DISTINCTION.md): γ as the source of experienced time.
- `simulations/decay_derivation.py`: decay spectrum is J-independent rational multiples of γ.
- Lindblad, G. (1976). "On the generators of quantum dynamical semigroups." Commun. Math. Phys. 48, 119–130.
