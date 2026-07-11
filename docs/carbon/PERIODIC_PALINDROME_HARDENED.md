# The Periodic Palindrome, Hardened: Mostly Ramp, Mixed Beyond It

**Date:** 2026-06-28
**Authors:** Tom + Claude
**Status:** Tier 1 computational result. A hardened, gate-first re-examination of the
periodic palindrome reported by
[`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py). It
qualifies that script's shuffle-null "p < 10⁻⁴" significance and complements
[The Periodic Palindrome and the V-Effect](PERIODIC_PALINDROME_VS_V_EFFECT.md).
**Script:** [`simulations/periodic_palindrome_gate.py`](../../simulations/periodic_palindrome_gate.py)

---

## The question

`periodic_palindrome.py` finds that across a period, element first ionization
energies form pairs `v_k + v_{N+1-k} ≈ const` (the F1 palindrome shape,
λ → −λ − 2σ), and calls it "significant" against a shuffle null at p < 10⁻⁴. But
**pair-sum-constant is satisfied exactly by any linear ramp**, and a shuffle null
only rejects "unordered". Any property that merely rises smoothly across a period
passes by smoothness alone. So: is the periodic palindrome a real F1 signature, or
a monotonic-smoothness artifact?

This is the [validation-vs-retrodiction](../../reflections/README.md) knife-edge in
concrete form: the elements are the framework's one free external anchor that can
say *no*, but only if the test is sharp enough to let them.

## The hardened test (gate-first)

Center each period (`w = v − mean`). Split about the period centre into the
F1-respecting antisymmetric part `a_k = (w_k − w_{N−1−k})/2` and the F1-breaking
symmetric part `s_k = (w_k + w_{N−1−k})/2` (pair-sum-constant ⟺ s = 0). Remove the
**linear ramp** `a_lin` (the trivial part any monotone property carries). The
discriminating residual is `r = a_non + s`.

**The null.** `R = E_non / (E_non + E_sym)` is the fraction of the non-ramp residual
that respects the mirror. A magnitude-preserving **sign-flip null** flips each
mirror partner's sign independently: it preserves every `|r_k|` exactly and
randomizes only the mirror relationship. Flipping a pair swaps its antisym↔symm
energy, so the null mean of R is exactly 0.5. R significantly above 0.5 means
mirror-respecting beyond chance; below 0.5 means mirror-breaking; ≈ 0.5 means F1
adds nothing beyond the ramp. Controls behave: a pure ramp leaves no residual; an
antisymmetric bump gives R > 0.5; a symmetric bump gives R < 0.5.

## Result

| Period | mirror pairs | R (mirror-respecting fraction) | verdict |
|--------|--------------|--------------------------------|---------|
| 2 (Li–Ne) | 4 | 0.380 | anti-F1, p = 0.50 (n.s.) |
| 3 (Na–Ar) | 4 | 0.247 | anti-F1, p = 0.25 (n.s.) |
| 4 (K–Kr) | 9 | 0.588 | F1-leaning, p = 0.22 |
| 5 (Rb–Xe) | 9 | 0.729 | **F1-leaning, p = 0.045** ✱ |
| 6 (Cs–Rn) | 16 | 0.583 | F1-leaning, p = 0.12 |
| **pool light (2,3)** | 8 | **0.327** | anti-F1, p = 0.20 (n.s.) |
| **pool heavy (4,5,6)** | 34 | **0.622** | **F1-leaning, p = 0.010** ✱ |

Ramp fraction of total period variance: 93 %, 90 %, 68 %, 67 %, 70 % (periods 2–6).
Allen 1989 electronegativity, the original's "strictest palindrome", has a post-ramp
residual energy of ~0.01 eV², pure smoothness.

Three things, honestly:

1. **The palindrome is mostly the trivial ramp.** Two-thirds to nine-tenths of each
   IE period is the linear ramp every rising property carries; Allen EN is ~100 %
   ramp. The shuffle-null "p < 10⁻⁴" was largely a monotonicity artifact.

2. **At the light elements (periods 2–3, carbon's home) the residual leans anti-F1.**
   Pooled R = 0.327: the non-ramp structure *breaks* the mirror more than it respects
   it, and the dominant breakers are exactly the half-filled-region shell anomalies
   (B/O, Al/S = p¹/p⁴). But with only 8 mirror pairs the lean is **not significant**
   (p = 0.20). The periodic table is too small here to settle it.

3. **At the heavy elements (periods 4–6) the residual is significantly
   mirror-respecting** (pooled R = 0.622, p = 0.010; period 5 alone p = 0.045). This
   is a real signal, but it plausibly reflects the point-symmetric (sigmoid)
   band-filling profile of a d/f period, which is antisymmetric-beyond-linear and
   which a generic band model produces too. F1 is consistent with it but does not
   *uniquely* predict it; F1-specific vs generic-band attribution is open.

## Verdict

F1-on-ionization-energy is **too loose a test**: pair-sum-constant ≈
antisymmetric-about-centre ≈ "rises across the period", which almost any property
does trivially. The elements are an honest but mostly-weak anchor for F1 here. The
cheap p < 10⁻⁴ is a smoothness artifact; the one sharp thing the light elements say
is that the shell anomalies *break* the mirror; and the significant heavy-period
signal cannot be pinned to F1 specifically. This neither validates nor refutes F1;
it shows the periodic table cannot, at this resolution, do either.

A real element-anchor for the framework needs a prediction a smooth ramp **cannot
fake**: the particle-hole-*odd* parity of the survivor mode against measured spectra
(see [The Survivor is Spin-Flip-Odd and Reflection-Odd](../../experiments/SURVIVOR_FLIP_AND_REFLECTION_ODD.md)),
the Π-vs-fermionic-particle-hole operator identity (the open analytical seam in
[the Majorana axis-modes lens](../../experiments/MAJORANA_AXIS_MODES.md)), or the T2 anisotropy
of [`carbon_painter_t2_anisotropy.py`](../../simulations/carbon_painter_t2_anisotropy.py).
Those are where F1 claims something non-trivial that the elements could genuinely
confirm or deny.

## Anchor

- Script: [`simulations/periodic_palindrome_gate.py`](../../simulations/periodic_palindrome_gate.py),
  built on [`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py).
- Complements [The Periodic Palindrome and the V-Effect](PERIODIC_PALINDROME_VS_V_EFFECT.md): that
  doc showed the *deviations* are textbook anomalies (the mechanism does not transfer);
  this one adds that even the *presence* of the palindrome is mostly smoothness, and
  what survives the ramp is mixed.
- Parent: [README.md](README.md).
