# D08: Crooks-like Rate Identity

**Derives:** ln(d_fast / d_slow) = 2 · artanh(Δd / (2Σγ))
**From:** F1 (palindrome equation: d + d' = 2Σγ)
**Status:** PROVEN (algebraic identity)

---

## Statement

For every palindromic eigenvalue pair with decay rates d and d' satisfying
d + d' = 2Σγ, define:

    d_slow = min(d, d')
    d_fast = max(d, d')
    Δd = d_fast − d_slow

Then:

    ln(d_fast / d_slow) = 2 · artanh(Δd / (2Σγ))

with effective inverse temperature β_eff ~ 1/Σγ in the linear regime.

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

The condition |Δd/(2Σγ)| < 1 is guaranteed because d_slow > 0 (from the spectral
gap D06: d_min = 2γ > 0) and d_fast < 2Σγ (since d_slow > 0).

## Linear approximation

For Δd ≪ 2Σγ (rates clustered near the center):

    artanh(x) ≈ x + x³/3 + ...

so:

    ln(d_fast / d_slow) ≈ Δd / Σγ

This has the form of a Boltzmann factor ln(p₁/p₂) = −β · ΔE with
β_eff = 1/Σγ. The palindrome mimics detailed balance structurally,
but this is NOT a thermodynamic identity: ⟨exp(−Δd)⟩ ≈ 0.93 ≠ 1,
so no Jarzynski equality holds.

## Why this is NOT Crooks

The Crooks fluctuation theorem relates forward and reverse process
probabilities via P_F(W)/P_R(−W) = exp(β(W − ΔF)). Our identity:

- Has the same functional form (log-ratio = function of difference)
- But involves decay rates, not work distributions
- The "temperature" β_eff = 1/Σγ is set by the dissipation strength, not a thermal bath
- No time-reversal symmetry is invoked; the palindrome is a spatial symmetry (Π)

The analogy is structural, not physical.

## Verification

Numerically verified for N=2−7 (all palindromic pairs satisfy the identity
to machine precision). See [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py).
