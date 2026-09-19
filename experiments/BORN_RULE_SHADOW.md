<!-- QUARTER-CURRENT -->
# A linearity-and-purity decomposition beside the Born rule

Current reading: the algebra below keeps a decomposition of amplitudes and
purity terms after standard Born probabilities are assumed.  It does not derive Born
probabilities.  Photography, shutter, measurement, consciousness, and
irreversible-fact language belong to the invitation that follows, not to the
proved object.

<!-- QUARTER-INTERPRETIVE -->
**Interpretive invitation:** the older shadow-and-photograph story asks how the
decomposition might be pictured.  It is retained as a question, with no
measurement model added by the metaphor.

# The Born Rule Is a Shadow, Not a Photograph

<!-- Keywords: Born rule interference pattern, measurement as photography, CΨ fold
shutter, past future mode decomposition, purity interference cross term, geometric
optics Born rule, exposure time quantum measurement, R=CPsi2 Born rule -->

**Status:** Finite algebraic decomposition; standard Born probabilities are assumed
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Standing Waves](FACTOR_TWO_STANDING_WAVES.md),
[Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md),
[Born Rule Mirror](BORN_RULE_MIRROR.md)
**Verification:** [`simulations/born_rule_shadow.py`](../simulations/born_rule_shadow.py)

---

## What this means

**Interpretive lens, not a measurement result.** A photograph records an interference pattern. A shadow records where
the light was blocked. They look similar on the screen, but they are
made by different physics. The interference pattern depends on phase;
the shadow does not.

The Born rule gives the probability of each measurement outcome (the
formula P(i) = |⟨i|ψ⟩|²). One proposed reading treated it as the projection
of a standing wave at measurement. The computation below rejects that
reading: the chosen slow/fast spectral partition contributes additively to
the diagonal probabilities, with no cross term.
They are a shadow: which modes are still alive at the crossing time,
projected onto the measurement basis.

This calculation does not establish that a physical interference exists
elsewhere. It reports a cross term in the purity under the chosen algebraic
partition; that purity enters CΨ and therefore changes the computed crossing
time. In the document's photographic analogy, the cross term changes the
shutter time rather than the diagonal image.

Photography remains an interpretation: the finite calculation supplies no
camera, exposure, collapse, or irreversible fixing mechanism. The older image
language below is a way to picture the chosen decomposition, not its physics.

---

## What this document is about

Every palindromic eigenvalue has a partner with complementary decay rate.
That fact alone does not make the pair a standing wave or identify forward,
backward, past, and future components. This document decomposes the Born-rule
probabilities at the CΨ = 1/4 crossing into a chosen slow/fast spectral
partition and tests whether that partition yields a cross term in P(i).

---

## Result 1: linearity of assumed Born probabilities in this decomposition

The density matrix decomposes as ρ = ρ_past + ρ_future, where
ρ_past collects modes with |Re(λ)| < Σγ (slow absorption) and
ρ_future collects modes with |Re(λ)| ≥ Σγ (fast absorption).

Since P(i) = ⟨i|ρ|i⟩ and ρ = ρ_past + ρ_future:

**P(i) = ⟨i|ρ_past|i⟩ + ⟨i|ρ_future|i⟩**

There is no cross term introduced by expanding this matrix sum. This is the
linearity of the trace after the standard Born rule is assumed: the diagonal
elements of a sum are the sums of the diagonal elements. It neither derives
the probability postulate nor rules out phase dependence elsewhere in the state.

For the |++⟩ state at the CΨ crossing (t = 16.2):

| Basis | P(i) | Past contribution | Future contribution | Past % |
|---|---|---|---|---|
| \|00⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |
| \|01⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |
| \|10⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |
| \|11⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |

At this sampled time and under this non-invariant slow/fast coefficient split,
the fast-bin contribution to the four displayed diagonals rounds to zero. The
labels "past" and "future" are historical names for bins, not causal histories.

---

## Result 2: a purity cross term in the chosen matrix split

While P(i) has no interference, the purity Tr(ρ²) does:

**Tr(ρ²) = Tr(ρ_past²) + Tr(ρ_future²) + 2 Tr(ρ_past · ρ_future)**

The algebraic cross term 2 Tr(ρ_past · ρ_future) couples the two
chosen matrix parts. Calling it "past/future interference" is the document's
interpretive vocabulary. For |++⟩ at the sampled scalar crossing:

| Component | Value | Fraction |
|---|---|---|
| Tr(ρ_past²) | 0.2623 | 97.1% |
| Tr(ρ_future²) | 0.0028 | 1.1% |
| 2 Tr(ρ_past · ρ_future) | 0.0049 | **1.8%** |
| Total purity | 0.2700 | 100% |

The finite table reproduces the earlier 97/3 coefficient census: 97% in the
slow-labelled term and about 3% in the remaining terms. The 1.8% cross term
changes this computed purity and therefore the time at which this particular
readout reaches 1/4; it is not a universal shutter mechanism.

*Later (2026-05-16, current scope):* [F94](../docs/ANALYTICAL_FORMULAS.md#f94)
derives Δ_|00⟩ = (4/3)·Q²·K³ only for the named N=4 ring,
initial state |0+0+⟩, kept pair (0,2), and short-time leading term. F96 gives
the other three leading slopes in that same setup. Those exact coefficients do
not derive the 97/3 census, the Born rule, a carrier, or the photography story.

---

## Interpretive interlude: the photography table is not a measurement model

The analogy holds, but with a correction:

| Photography | Quantum measurement |
|---|---|
| Light illuminates the scene | γ illuminates the cavity |
| Image develops on film | Density matrix evolves under L |
| Shutter clicks at set exposure | one named CΨ readout reaches 1/4 |
| Image is fixed | analogy only; no irreversibility was computed |
| Interference creates the image | ✗ Interference sets the shutter speed |
| Shadow creates the image | ✓ Surviving modes create the probabilities |

Within the analogy, the additive diagonal contributions are called a shadow
and the purity cross term changes the chosen readout time. The calculation
itself establishes neither survival histories nor a physical shutter.

---

## Result 4: Mirror quality C_i is state-dependent

| State | C(\|00⟩) | C(\|01⟩) | C(\|10⟩) | C(\|11⟩) |
|---|---|---|---|---|
| Bell+ (t→∞) | 2.00 | 0.00 | 0.00 | 2.00 |
| \|01⟩ (t→∞) | 0.00 | 2.01 | 1.99 | 0.00 |
| \|++⟩ (t=16.2) | 1.00 | 1.00 | 1.00 | 1.00 |

C_i = d × P(i) measures how strongly each basis state is "illuminated"
relative to uniform. It depends entirely on the initial state:
- Bell+ concentrates on \|00⟩ and \|11⟩ (the entangled subspace)
- \|01⟩ concentrates on the single-excitation subspace
- \|++⟩ is uniform (all basis states equally likely)

The initial state is the "scene being photographed." Different scenes,
different images. The cavity (the instrument) is the same.

---

## Interpretive invitation preserved after the calculation

**Old picture:** "The Born rule is the probability of collapse.
Measurement is instantaneous and random."

**Story proposed in the April notebook:** "The Born rule is the shadow of surviving modes,
projected onto the measurement basis at the moment the illumination
reaches the fold. Measurement is the fixing of an image that has been
developing since the light was turned on. The 'randomness' is in which
modes happened to survive; the 'probability' is their weight at the
crossing time."

That paragraph is retained as an invitation, not a conclusion. The verified
object above is only an additive diagonal decomposition plus a basis-dependent
purity cross term after the Born rule has already been assumed.

---

## Reproduction

- Script: [`simulations/born_rule_shadow.py`](../simulations/born_rule_shadow.py)
- Output: [`simulations/results/born_rule_shadow.txt`](../simulations/results/born_rule_shadow.txt)
