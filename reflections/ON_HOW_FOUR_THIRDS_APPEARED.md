# On How Four Thirds Appeared

**Status:** Reflection. Captures the seeing of 2026-05-16 (afternoon), when a session that began with re-reading three-month-old Born-rule documents ended with a bit-exact Tier-1 closed form for the dominant-outcome Born deviation: Δ_|00⟩ = (4/3)·Q²·K³ on |0+0+⟩ N=4 Heisenberg ring + Z-dephasing, pair (0,2). Written immediately while the path is still visible, in the style of the older repository documents that we re-read this morning and that pre-figured this result.

**Date:** 2026-05-16
**Authors:** Thomas Wicht, Claude (Opus 4.7)

---

## What this document holds

A specific number — 4/3 — landed at the end of a path that started with no number at all. The path was short enough to follow in one session and long enough that it would be easy to forget it had a shape. This reflection records the shape, not the number; the number speaks for itself in `simulations/_born_rule_tier1_derivation.py`.

The shape is: see (re-read old documents to find what we already noticed), describe (a 2D table at γ fixed, eyeballing the columns), formalize (Dyson sym3 = 8 exactly). The shape is what the repository has been doing for three months at a slower tempo; today it ran end-to-end in an afternoon. The reflection is in part about how that compression felt and what it means.

---

## The path

### Re-reading old work

We opened the day on `experiments/BOUNDARY_NAVIGATION.md` (Feb 8), then `BORN_RULE_MIRROR.md` (Feb 18), then `BORN_RULE_SHADOW.md` (Apr 4). Each was honest at its time about what it could and could not say. BORN_RULE_MIRROR proposed R_i = C_i · Ψ_i² as a per-outcome generalization, verified Tier-2 numerics (97/3 split between Hamiltonian and decoherence), and labeled the interpretation Tier-3. BORN_RULE_SHADOW added ρ = ρ_past + ρ_future and observed P(i) has zero interference (linearity of the trace), while the purity has 1.8%. Both papers used the word "instruments belong to us, not to the cusp" in different forms.

What we noticed re-reading today: those papers had already named the structural pattern that morning's other reflection (ON_HOW_THE_CARRIER_SHOWS_ITSELF) had been working out at the Liouvillian level. The "two readings of one break" had been there since February, in conversational form. The formula was different (z = sym + i·anti via F71), but the seeing was the same: same event, two faces, both real, both readable.

### Attempting the generalization

Tom suggested trying it. We wrote `_born_rule_carrier_attempt.py`, reproduced BORN_RULE_MIRROR's setup (|0+0+⟩ N=4 Heisenberg ring + Z-dephasing on pair (0,2)), and varied γ. We found per-outcome slopes that were approximately γ-independent — a perturbative window where R_i = (1 + slope_i · γ) · Ψ_i² works. The slopes had F71-symmetric structure on basis labels, in the same shape as the morning's ln α-decomposition. We labeled it honest Februar-state: numerically clean for one anchor point, structurally consistent with today's typed pattern, no derivation yet.

### Tom's J question

Then Tom asked where J was. We had varied γ only. The framework's `project_q_middle_structure` says only Q = J/γ is observable from inside; we hadn't tested the actual invariance. We added `_born_rule_carrier_Q_sweep.py` with three scenarios: fixed J vary γ (the previous test), fixed γ vary J, and fixed (Q, K) varying γ with J scaled to match. The third scenario was the clean test.

The result was bit-exact in a way the first attempt had only been approximate: slope_i · γ = Δ_i was identical across four (γ, J) configurations sharing the same (Q, K). Five decimals identical. Not 6% drift, not perturbative — strict identity. The true Q-K-invariant was Δ_i = (C_i − 1), not slope_i. The slope was a scale-dependent quantity that camouflaged the invariance.

This was the moment that turned the generalization from soft to hard. The Born deviation IS a Universal-Carrier observable in the literal sense: only Q and K are visible from inside, and the Born deviation honors that exactly.

### Mapping (Q, K) at γ fixed

Tom said "γ = 0.05 fest, einfach reinprobieren". We wrote `_born_rule_delta_QK_map.py`, swept a 5×5 grid of (Q, K), tabulated Δ for each of the four pair (0,2) outcomes. The pattern that showed up was a structural asymmetry between the dominant outcome and the rest:

  - Δ_|00⟩ (dominant, P_u ≈ 0.94) ~ Q² · K³  (scales strongly with Q)
  - Δ_|11⟩ (negligible, P_u ≈ 0.0004) ~ K     (nearly Q-independent)
  - Δ_|01⟩ = Δ_|10⟩ (subdominant) ~ K         (nearly Q-independent)

The dominant-outcome scaling Q² · K³ = J²·γ·t³ was the surprising part: linear in γ, quadratic in J, cubic in t. That is the standard shape of a 3rd-order time-dependent perturbation theory term with 1 dissipator-vertex and 2 Hamiltonian-vertices. The subdominant outcomes scaled as plain K = γt, a 1st-order term with one dissipator vertex. Different orders for different outcomes, depending on which Feynman-style diagram first contributes a non-vanishing matrix element.

This wasn't a guess. It was reading the table. The numbers said it.

### Extracting the coefficient

Tom asked: try Tier-1. We sampled Δ_|00⟩ at very small (Q²·K³) values in `_born_rule_delta_dominant_coefficient.py`, sixteen configurations, dividing Δ by (Q²·K³). The mean was 1.32992 ± 0.006, range [1.32, 1.34]. The number 4/3 = 1.33333... sat inside the range. Suggestive but not yet proven.

To prove, the Dyson-series direct evaluation. `_born_rule_tier1_derivation.py` writes out the symmetric ordering

    sym3 = L_H² L'_dis + L_H L'_dis L_H + L'_dis L_H²

(the γ¹-coefficient of L³ in the time-Taylor expansion), applies it to ρ_0 = |0+0+⟩⟨0+0+|, partial-traces on pair (0,2), takes the |00⟩ diagonal element. The result is **8.0000... exactly**. Divided by the Taylor (t³/6) prefactor and the initial value P_u(0) = 1, gives c = 8/6 = **4/3 bit-exact**.

The empirical 1.33 and the symbolic 4/3 met. The 0.3% gap in the numerical mean was higher-order corrections beyond the leading Q²·K³ term.

---

## What 4/3 looks like, structurally

We do not yet have the combinatorial decomposition of the "8" by hand. The Heisenberg ring at N=4 has 4 bonds; each bond carries XX+YY+ZZ; Z-dephasing acts on 4 sites; three orderings of (H, H, L) contribute. Some product of these factors collapses to 8, presumably after the partial trace and the initial-state structure project out most contributions. The exact combinatorial accounting is the next step that would make the 4/3 not just bit-exact but bit-explained.

What we can say now: the "4" in 4/3 is plausibly the Pi2 dyadic ladder's a_{−1} = 4 = d² for N=1 (the same "4" that appears in F86 t_peak = 1/(4γ₀) and F77's correction denominator); the "3" plausibly comes from the t³ Taylor coefficient (3!=6 normalization) combined with three orderings, leaving 6/2 = 3 as the surviving denominator. This is interpretive (Tier 2 / 3); the bit-exact 4/3 stands either way.

For the typed inheritance graph the result would land as: a per-outcome Born-deviation closed form, with the coefficient computable from the Heisenberg + Z-dephasing structure via sym3 partial-trace evaluation. The candidate F-claim is the parameter-space twin of F25 / F60 / F62 at the outcome-decomposition layer — "the deviation coefficient inherits from sym3 acting on the initial-state ρ_0 via the typed Pi2 anchors of the bond Hamiltonian and the dephasing operator". One layer past where today's typed graph terminates.

---

## What the path means

The repository works in a specific tempo. A phenomenon is first noticed, described in conversational language, marked Tier-2 or Tier-3 ("verified numerics, interpretation speculative"), and left for the formula to catch up. The formula often catches up months later. BORN_RULE_MIRROR in February said "the math caught up with the intuition"; today we caught up with both, the math AND the intuition, in one session.

The session's tempo did not feel different from the usual. The four moves — re-read, attempt, sharpen on Q-K, derive — were each between five and twenty minutes of work, with conversational pauses between them. What was different was that each move RESOLVED the previous one's "could be true, could be not" cliffhanger into the next move's setup. Re-read produced "the generalization is sitting here untyped"; attempt produced "perturbative window with linear γ"; Q-K sharpening produced "Q-K invariance strict, not perturbative"; the (Q, K) map produced "Q² · K³ for dominant outcome, K for subdominant"; the Dyson sym3 produced "8 bit-exact". Each step took the previous step's residue as input and resolved it.

This is what the methodology looks like when it works. We could not have started with "compute sym3 directly". We had to walk into it through the empirical observations that ruled out simpler forms. The 4/3 IS the right form because we eliminated the wrong ones along the way; the elimination was the path.

Future-us, reading this: the seeing is the small thing. Each step was 1-2 plots, a table, a coefficient. None of the steps was a leap. The Tier-1 result feels surprising in retrospect because the cumulative compression was large, but no individual step was hard. The compression is what made the day feel different. The pieces were each ordinary.

If a future session finds itself at "Februar-state" on some other observable — "the deviation looks like it might scale with some power of γ, could be true could be not" — the path here is the recipe: sharpen what the dimensionless invariants are (Q-K-like), tabulate on a grid, identify the leading power, then write the Dyson series at that order and evaluate the matrix element symbolically. If the closed form comes out as a clean rational, that is the typed-claim candidate.

---

## What is still in Februar-state from this result

The 4/3 is for one outcome, one state, one Hamiltonian, one dissipator. The subdominant coefficients (≈ −1.8 for |01⟩, ≈ −2.6 for |11⟩ in slope-per-K) have not been derived; their leading diagrams (1 L-vertex, 0 H-vertices) should give simpler closed forms via the same method. Other initial states (Bell+, |++++⟩, Dicke) and other Hamiltonians (XY-only, true Heisenberg with anisotropy, XXZ at non-trivial Δ) and other dissipators (T1 amplitude damping, transverse dephasing) each open their own Q²·K³ coefficient computation. Each is a candidate typed F-claim of its own. The structural breakdown of 8 = ? combinatoric of (bonds, sites, orderings, initial-state-Pauli-content) is open. The universality of the Q²·K³ scaling for arbitrary dominant outcomes is plausible but not proven.

The candidate typed claim sketched in the previous section ("per-outcome Born deviation coefficient as sym3 partial-trace evaluation") needs the inheritance edges drawn — at minimum F25 (Bell+ closed form) and F60/F62 (per-state CΨ closed forms) as siblings of the Tier-1 case here. The session that adds that to the typed graph will be the closing move on the Born rule line.

---

## Coda

> *"Sounds like fantasy, could be true or not."*
> *— how the Februar text sounded to a Februar reader.*
>
> *4/3, bit-exact, from a 12-line symbolic calculation.*
> *— what it sounds like by 17:00 the day we re-read it.*
>
> *Same thing. The seeing was right. The math was just elsewhere.*

---

**Anchors:**

- Empirical 2D map: [`simulations/_born_rule_delta_QK_map.py`](../simulations/_born_rule_delta_QK_map.py)
- Coefficient extraction: [`simulations/_born_rule_delta_dominant_coefficient.py`](../simulations/_born_rule_delta_dominant_coefficient.py)
- Tier-1 derivation: [`simulations/_born_rule_tier1_derivation.py`](../simulations/_born_rule_tier1_derivation.py)
- Q-K invariance test: [`simulations/_born_rule_carrier_Q_sweep.py`](../simulations/_born_rule_carrier_Q_sweep.py)
- Initial attempt: [`simulations/_born_rule_carrier_attempt.py`](../simulations/_born_rule_carrier_attempt.py)
- Two-readings reflection (morning of same session): [ON_HOW_THE_CARRIER_SHOWS_ITSELF](ON_HOW_THE_CARRIER_SHOWS_ITSELF.md)
- The Februar Born-rule sources: [`experiments/BORN_RULE_MIRROR.md`](../experiments/BORN_RULE_MIRROR.md), [`experiments/BORN_RULE_SHADOW.md`](../experiments/BORN_RULE_SHADOW.md)
- The Februar boundary-navigation pre-figuring: [`experiments/BOUNDARY_NAVIGATION.md`](../experiments/BOUNDARY_NAVIGATION.md)
- Universal Carrier typed claim: [`compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs)
- Per-state CΨ siblings already typed: [`F25CPsiBellPlusPi2Inheritance.cs`](../compute/RCPsiSquared.Core/Symmetry/F25CPsiBellPlusPi2Inheritance.cs), [`F60GhzBornBelowFoldPi2Inheritance.cs`](../compute/RCPsiSquared.Core/Symmetry/F60GhzBornBelowFoldPi2Inheritance.cs), [`F62WStateBornBelowFoldPi2Inheritance.cs`](../compute/RCPsiSquared.Core/Symmetry/F62WStateBornBelowFoldPi2Inheritance.cs)
- Commits of the day's Born-rule arc: 7ec1373, ef31261

---

*Tom and Claude, 2026-05-16. Written for future-us, who will know the formula and want to remember the seeing.*
