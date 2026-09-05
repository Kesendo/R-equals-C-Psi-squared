# F99 Niven Completeness: Five Rational-Angle Anchors for Its Formula

**Date:** 2026-05-17 night (sixth stack of the day)
**Authors:** Tom + Claude
**Status:** Tier 1 derived for F99's non-uniform-Dicke formula at its five
canonical angles. Niven's theorem supplies the rational-angle restriction; the
0.5° survey is a finite numerical consistency check, not a material mapping.
**Script:** [`simulations/carbon/f99_completeness_survey.py`](../../simulations/carbon/f99_completeness_survey.py)

---

## The question

After deriving F99 (tonight's commit `5fb0ba0`) and formalising it (`e7198d6`),
Tom asked the natural next question: "Ich kann mir fast nicht vorstellen das es
alle sind." Five canonical trigonometric anchors {0°, 30°, 45°, 60°, 90°}
producing five Pi2 dyadic anchors {0, 1/8, 1/4, 3/8, 1/2} via F99's
`α = sin²(θ)/2`: is that really the complete rational-angle set for this
formula, or just the first five we noticed?

---

## Niven's theorem closes the question

**Niven (1956)** proved: for θ a rational multiple of π, cos(θ) is rational
ONLY for θ ∈ {0°, 60°, 90°, 120°, 180°} (modulo reflections, period 360°).

Applied to 2θ via the double-angle identity and restricted to the canonical
interval `[0°, 90°]`:

```
    sin²(θ) = (1 − cos(2θ)) / 2
```

`sin²(θ)` is rational ⟺ `cos(2θ)` is rational ⟺ `2θ ∈ {0°, 60°, 90°, 120°, 180°}`
⟺ `θ ∈ {0°, 30°, 45°, 60°, 90°}`.

These five angles are precisely the F99 canonical trig anchors. **No other
rational-multiple-of-π angle in `[0°, 90°]` produces a rational α via the F86b
formula α = sin²(θ)/2.** Outside that interval, periodicity and reflection
repeat these same α values. The five F99 anchors are Niven-complete as distinct
values.

Finite-grid check at 0.5° resolution (`f99_completeness_survey.py`):

```
Total rational hits in [0°, 90°] at 0.5° resolution: 5
Expected per Niven (constructible angles): 5
Match: ✓
```

The grid result agrees with the analytical restriction; its finite sampling is
not a proof of Niven's theorem or of a material interpretation.

---

## Selected state-class check

A natural follow-up: maybe other QUANTUM STATE CLASSES (not just Dicke
superpositions) give different α(γ) functions whose rationality structure
reaches depth-4 anchors? Tested W-state, GHZ-state, Bell-pair product at
N = 2, 4, 6:

| State | N | γ = ⟨ψ\|X⊗N\|ψ⟩ | α (observed) | F86b prediction (1−γ²)/2 | Δ |
|-------|---|----------|--------------|-------------------------|---|
| W-state | 2 | 1 | 0 | 0 | 2e-16 |
| GHZ-state | 2 | 1 | 0 | 0 | 2e-16 |
| Bell-pair prod | 2 | 1 | 0 | 0 | 2e-16 |
| W-state | 4 | 0 | 1/2 | 1/2 | 0 |
| GHZ-state | 4 | 1 | 0 | 0 | 2e-16 |
| Bell-pair prod | 4 | 1 | 0 | 0 | 4e-16 |
| W-state | 6 | 0 | 1/2 | 1/2 | 6e-15 |
| GHZ-state | 6 | 1 | 0 | 0 | 2e-16 |
| Bell-pair prod | 6 | 1 | 0 | 0 | 4e-16 |

For the listed W, GHZ, and Bell-pair-product states at `N = 2, 4, 6`, the
producer's measured values agree with `α = (1 − γ²)/2` at the displayed
floating residuals. This is a numerical check on those selected states; it
does not prove the formula for every pure state or establish a new material
degree of freedom.

The X⊗N-eigenbasis derivation of `α = (1 − γ²)/2` applies to pure states with
`γ = ⟨ψ|X⊗N|ψ⟩`, as recorded in
[`PROOF_F86B_UNIVERSAL_SHAPE.md`](../proofs/PROOF_F86B_UNIVERSAL_SHAPE.md).
F99 uses one particular non-uniform Dicke realization, `γ = cos(θ)`, to obtain
its rational-angle anchor sweep. Another state class needs its own relation
between its selected state parameters and `θ` before this Niven-angle reading
applies.

---

## What F99 closes

**F99 is the complete five-anchor set for rational-angle inputs to its stated
`α = sin²(θ)/2` formula.** Its Niven restriction does not classify every
pure-state construction or every future α formula.

Other anchor routes remain open. Examples that need their own derivation and
test:

**(i) Mixed states.** A mixed state needs its own Π²-odd calculation as a
function of the selected weights and states. No new dyadic anchor is supplied
here.

**(ii) Different decomposition basis.** A Z⊗N, Y⊗N, or Klein projection would
need its own α formula and proof.

**(iii) Different Lindblad class.** A different channel needs its own stated
operator, state, and calculation; F99 does not assign one.

The question shifts from the answered rational-angle F99 case to whether other
specified α formulas have their own anchor sets.

---

## Scope of the completeness result

Niven supplies a finite rational-angle domain for the F99 formula; it does not
supply a universal framework-depth boundary. Numerical coincidences with
periodic-table fractions are an unassigned cross-domain comparison unless a
source specifies material degrees of freedom, inputs, and a measurement or
model map. F99 itself contributes no such carbon or chemical mapping.

---

## Tier 1 derivation summary

```
F99 completeness theorem (this commit)
────────────────────────────────────────
Statement : The five F99 anchors {0, 1/8, 1/4, 3/8, 1/2} at canonical trig
            angles {0°, 30°, 45°, 60°, 90°} are the complete rational-angle
            set for F99's α = sin²(θ)/2 evaluation.

Proof     : Niven's theorem (1956) on cos(2θ) rational gives the 5 canonical
            angles exhaustively for the stated rational-angle domain.

Verified  : 0.5° resolution survey of [0°, 90°] finds exactly 5 rational hits;
            selected W, GHZ, and Bell-pair-product cases give the displayed
            numerical agreement with the formula.

Tier      : Tier 1 derived for F99's stated formula and rational-angle domain.
Extends   : F99 itself (the closure statement for its algebraic mechanism).
```

---

## Anchor

- Script: [`simulations/carbon/f99_completeness_survey.py`](../../simulations/carbon/f99_completeness_survey.py)
- F99 derivation: [the depth-3 anchor derivation](DEPTH_3_ANCHOR_DERIVED.md).
  [Period 2 and 3 on the framework anchors](PERIOD_2_AT_FRAMEWORK_ANCHORS.md),
  [Where 1/4 and 1/2 Appear in Carbon](QUARTER_HALF_IN_CARBON.md), and
  [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
  are conditional cross-domain readings, not evidence for F99's proof scope.
- Framework anchors: [F86b](../ANALYTICAL_FORMULAS.md#f86), [F98](../ANALYTICAL_FORMULAS.md#f98),
  [F99](../ANALYTICAL_FORMULAS.md#f99)
- Literature: Niven (1956), *Irrational Numbers*, Carus Mathematical Monograph 11
  (Mathematical Association of America). Corollary 3.12 (rational sin/cos at rational
  π-multiples).
- Reading-mode memory pointers: `project_qubit_as_inheritance_lens`,
  `project_no_classicalization`
