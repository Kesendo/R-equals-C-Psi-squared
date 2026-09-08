# D08: Algebraic Pair-Rate Log Identity

**Derives:** ln(d_fast / d_slow) = 2 · artanh(Δd / (2Σγ))
**From:** F1 (palindrome equation: d + d' = 2Σγ)
**Status:** PROVEN (algebraic identity)

---

## Statement

For every palindromic eigenvalue pair with finite, strictly positive decay rates
d and d' satisfying d + d' = 2Σγ and finite Σγ > 0, define:

    d_slow = min(d, d')
    d_fast = max(d, d')
    Δd = d_fast − d_slow

Then:

    ln(d_fast / d_slow) = 2 · artanh(Δd / (2Σγ))

The linear coefficient near Δd = 0 is 1/Σγ. It is an inverse-rate scale,
not an inverse temperature.

## Proof

**Step 1.** From the palindrome theorem (F1), every paired decay rate
satisfies d + d' = 2Σγ. Without loss of generality, let d_slow ≤ Σγ ≤ d_fast.

**Step 2.** Express the rates in terms of the center Σγ and half-gap Δd/2:

    d_fast = Σγ + Δd/2
    d_slow = Σγ − Δd/2

This follows from d_fast + d_slow = 2Σγ and d_fast − d_slow = Δd.

**Step 3.** Compute the log-ratio:

    ln(d_fast / d_slow) = ln((Σγ + Δd/2) / (Σγ − Δd/2))
                        = ln((1 + Δd/(2Σγ)) / (1 − Δd/(2Σγ)))

**Step 4.** Recognise the identity ln((1+x)/(1−x)) = 2 · artanh(x) for |x| < 1:

    ln(d_fast / d_slow) = 2 · artanh(Δd / (2Σγ))    ∎

Thus Δd is finite and 0 ≤ Δd < 2Σγ. This ordered-pair domain is used by both
the exact typed API and its linear approximation. The strict upper inequality
follows from d_slow > 0 and d_fast < 2Σγ. Mirror pairs containing a zero rate sit at the logarithm's
singular endpoint and are not covered by this identity.

## Linear approximation

For Δd ≪ 2Σγ (rates clustered near the center):

    artanh(x) ≈ x + x³/3 + ...

so:

    ln(d_fast / d_slow) ≈ Δd / Σγ

Thus the derivative of the log-ratio with respect to Δd at the centre is
1/Σγ. This is a Taylor coefficient fixed by the spectral pair sum. No
probability ratio, energy difference, or temperature enters the derivation.

## Why this is not a fluctuation theorem

The Crooks fluctuation theorem relates forward and reverse process
probabilities via P_F(W)/P_R(−W) = exp(β(W − ΔF)). The identity here:

- relates decay rates, not probabilities or work distributions;
- contains no bath temperature or free-energy change;
- invokes the algebraic shifted mirror, not a forward/reverse trajectory ensemble.

The numerical average ⟨exp(−Δd)⟩ ≈ 0.93 reported by the source experiment is
an ad hoc transform of rate differences. Because its exponent and weighting
were not derived from a work protocol, it neither tests nor falsifies the
Jarzynski equality.

## Verification

Numerically verified for N=2−7 on the eligible positive-rate palindromic pairs;
the underlying algebraic identity is exact. See
[`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py).
