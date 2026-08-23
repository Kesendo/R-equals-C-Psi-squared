# Shaping the Dephasing Profile: Practical γ Control for Quantum State Transfer

<!-- Keywords: dephasing profile optimization, quantum state transfer noise shaping,
total dephasing budget sum gamma confound, dynamical decoupling mediator chain,
V-shape dephasing gradient, quantum transistor gate noise, Lindblad noise engineering,
open quantum system control, gamma profile quantum channel,
spin chain mutual information optimization, assignment vs multiset -->

**Status:** Computationally verified (all simulations reproducible)
**Date:** March 22, 2026, last refreshed 2026-08-23 (the change history lives in git)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** [using_gamma.py](../simulations/using_gamma.py) (March); [gamma_control_two_lever_check.py](../simulations/gamma_control_two_lever_check.py) (the 2026-08-23 gate over [mediator_bridge.py](../simulations/mediator_bridge.py))
**Data:** [using_gamma.txt](../simulations/results/using_gamma.txt)

---

## Abstract

In a 5-qubit mediator chain under Heisenberg coupling and local dephasing,
we test strategies for shaping the spatial dephasing profile to maximize
mutual information between the end pairs. Two levers govern the outcome.
Within a fixed shape, less total dephasing Σγ is always better,
monotonically in every family tested. At fixed Σγ, the ASSIGNMENT is a
large effect with a simple law: **concentrating the budget beats spreading
it, and the chain centre is the cheapest seat**, worth up to +46% over
uniform for the whole budget on the mediator site; the March "V-shape" is a
partial concentration worth +6%. The March headline of "+124% from the
V-shape" compared arms with different Σγ (0.13 against 0.25) and so
measured mostly the budget, and its conclusion that gate noise HELPS is
refuted by the one-factor controls, including two rows the March run
computed and never reported: the mediator's own noise always harms; it
merely harms least, because sensitivity to noise falls from the chain ends
toward the centre. Dynamical decoupling follows the same positional law,
the two ends equally first, the centre last. AC modulation of the gate and
state-dependent feedback both fail. The March "time-resolved decoder"
loses its resolution figure to its missing control, run here: a doubled
γ₂ imprints at most 0.013 in absolute MI at any lag, systematic but
small, and the readable form of that fact is the full channel analysis
([γ as Signal](GAMMA_AS_SIGNAL.md): 15.5 bits capacity at 1% readout
noise, 100% classification noiseless).

**Scope on "always harms", added 2026-08-23.** That clause is measured on the
preparations swept here, and it turns out to be preparation-dependent rather than
general. The counterexample is on a **different system**, so read it as a limit
on the word "always" and not as a correction of any number on this page:
[The Blind Site](THE_BLIND_SITE.md) §7 runs an **eleven**-site chain (not the
five-site one used here), observes I(A:B) between five-site halves (not the
end-pair MI used here), and prepares a reflection-odd single excitation, for
which the centre is an exact node of every mode the state occupies. There the
mediator's own γ eventually **helps**: the response to γ_M is not monotone, and
past a crossover the correlation goes back **up**, reaching a 10.8 % gain at
γ_M = 100 while the window-mean occupation of the centre falls from 0.0663 to
0.0310. Strong dephasing at the seat every path crosses blocks the transport, and
the state keeps more of the correlation it began with.

**Where that crossover sits, so the limit is not overstated.** At a per-site
background of 0.05, this page's own baseline, the crossover needs γ_M between 5
and 20, which is one to two decades above anything swept here; at a background of
0.2, which is above this page's range, γ_M = 0.5 already helps by 0.024 %. So the
counterexample does not reach into the rates this page measured and no number
above changes. What it limits is the word **always**: the sign of the mediator's
own noise is a fact about the prepared state and the rate, not about the seat.

---

## Background: total noise first, then the seat

### The standard approach

In quantum computing, dephasing noise (γ) is treated as a budget to be
minimized. This experiment asks the follow-up question: once the total
budget Σγ is fixed, does its spatial ASSIGNMENT across the chain matter,
and which seats should carry it?

### What this experiment shows

Two levers, and neither dominates the other. The BUDGET: within any fixed
shape, less Σγ means more end-to-end mutual information, monotonically
(uniform: 0.709 at Σγ = 0.13 falls to 0.338 at 0.25; all-on-mediator:
1.039 at 0.13 falls through 0.747 at 0.25 to 0.577 at 0.40). A comparison
whose arms carry different Σγ inside one shape therefore measures the
budget. The ASSIGNMENT: at fixed Σγ the spread across arrangements is
large, +51% between the best and worst at Σγ = 0.13, and the law is
concentration: putting the whole budget on one site beats spreading it, at
every seat tried, and the chain centre is by far the cheapest seat. The
two levers cross, so neither is "the dominant" one: all of 0.25 on the
centre beats uniform 0.13 at nearly twice the budget (0.747 vs 0.709 at
t = 5.0, +16.7% on the trajectory integrated over t = 0-10). The mediator's
own noise still harms, always and monotonically; concentration at the
centre wins because sensitivity to noise falls from the ends toward the
centre, so the centre is where each unavoidable unit does the least
damage. The March version of this page instead concluded that gate noise
helps ("the gate needs to be the noisiest node"); that conclusion was the
Σγ confound, and the one-factor rows in
[the March data file](../simulations/results/using_gamma.txt) already
contradicted it (see the results table below).

**What the repo already held on this, swept 2026-08-23:** the Σγ confound
was first written down on 2026-05-31 in a local design spec that never
became a tracked page, so no committed document carried it until now.
`hypotheses/THE_OTHER_SIDE.md` §22 measured the one-factor γ_M sweep on the
N=5 bridge in March (0.44 → 0.12) and concluded "the mediator must be quiet
for the bridge to be wide", consistent with the controls here and in direct
contradiction to this page's March conclusion.
[Mediator as Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md)
prescribes the same polarity (reduce γ_M to open the channel, DD on the
pairs and NOT on M). The formula registry's F91 states, for its own
object (the F71-refined block spectra over the pair-sum profile), that the
dependence is on the ASSIGNMENT and never on the bare multiset; the
comparison repaired here did not even hold the multiset's sum fixed. `docs/CAUGHT_ERRORS.md` (2026-08-23 entry) records
the error shape.

### Connection to the dephasing channel

If the spatial γ profile affects internal observables strongly, then the
profile itself is **readable from inside** the system. This led to the
[γ as Signal](GAMMA_AS_SIGNAL.md) experiment, which proved that γ profiles
encode information with 15.5 bits of theoretical channel capacity.

---

## System Setup

| Parameter | Value |
|-----------|-------|
| System | 5-qubit linear chain 0-1-2-3-4; pair A = {0,1}, mediator = 2, pair B = {3,4} |
| Coupling | Heisenberg (J = 1.0 between nearest neighbors) |
| Baseline dephasing | γ = 0.05 per qubit (uniform), Σγ = 0.25 |
| Initial state | Bell(0,1) ⊗ \|0⟩₂ ⊗ \|+⟩₃ ⊗ \|+⟩₄ (`bell_A_0M_pp_B`) |
| Observable | MI({0,1} : {3,4}) at t = 5.0 |
| Noise model | Local Z-dephasing (σ_z per qubit) |
| Master equation | Lindblad: dρ/dt = −i[H,ρ] + Σᵢ γᵢ(σ_z⁽ⁱ⁾ρσ_z⁽ⁱ⁾ − ρ) |

---

## Results

### The March profiles, with the column that decides

All eight static profiles of the March run, five of which the March page
never reported, the two `gate-*` rows among them. Sorted by Σγ; within this particular set of
eight, MI rises as Σγ falls (that is the set, not a law: the matched-Σγ
table below contains profiles that break the ordering across shapes):

| Profile | γ per site | Σγ | MI(A:B) |
|---------|-----------|-----|---------|
| gate-closed | [.05,.05,.20,.05,.05] | 0.40 | 0.192 |
| uniform (baseline) | [.05,.05,.05,.05,.05] | 0.25 | 0.338 |
| gate-open | [.05,.05,.01,.05,.05] | 0.21 | 0.416 |
| inverse V | [.05,.03,.01,.03,.05] | 0.17 | 0.530 |
| quiet receiver | [.05,.05,.05,.01,.01] | 0.17 | 0.569 |
| ramp-down | [.05,.04,.03,.02,.01] | 0.15 | 0.628 |
| ramp-up | [.01,.02,.03,.04,.05] | 0.15 | 0.631 |
| V-shape | [.01,.03,.05,.03,.01] | 0.13 | 0.755 |

Three readings, all from this one table. The `gate-open`/`gate-closed` pair
is the direct one-factor test of the mediator's own noise: only γ₂ moves,
and MI falls monotonically as it rises (0.416 → 0.338 → 0.192). A
one-factor move is of course also a Σγ move; what the pair shows is that
ADDING a unit at the centre harms, which is all the refutation of "gate
noise helps" needs. The two
ramps share Σγ = 0.15 and differ by 0.4%, though little should be read
into that: the two are spatial reversals of each other and the observable
is nearly reversal-symmetric, so their near-equality is close to forced
and the 0.4% is the initial state's asymmetry. And the second matched
pair, across different shapes, inverse V against quiet receiver at
Σγ = 0.17, is won by the quiet receiver; its arms differ in more than one
seat (relative to uniform the inverse V removes noise from sites 1, 2 and
3, the quiet-receiver profile from sites 3 and 4), so the clean per-seat
ranking is the matched-removal DD table below, whose positional law this
pair agrees with.

### Σγ matched: the assignment axis in full

All arms at Σγ = 0.13 (the inverse V rescaled by 13/17 exactly; the
single-site and sharp-V rows computed 2026-08-23, twice independently):

| Profile (Σγ = 0.13) | γ per site | MI(A:B) |
|---------------------|-----------|---------|
| all on the mediator | [0, 0, .13, 0, 0] | 1.039 |
| sharp V | [0, .015, .10, .015, 0] | 0.888 |
| all on site 3 | [0, 0, 0, .13, 0] | 0.809 |
| all on site 1 | [0, .13, 0, 0, 0] | 0.789 |
| V-shape | [.010,.030,.050,.030,.010] | 0.755 |
| all on site 4 | [0, 0, 0, 0, .13] | 0.742 |
| all on site 0 | [.13, 0, 0, 0, 0] | 0.740 |
| uniform | [.026,.026,.026,.026,.026] | 0.709 |
| inverse V rescaled | (13/17) × [.05,.03,.01,.03,.05] | 0.690 |

At matched total noise the V-shape's advantage over uniform shrinks from
the March "+124%" to **+6%**, but the assignment axis itself is wide: the
full spread at this budget is +51%, and its shape is simple. Every
single-site concentration beats uniform, whichever seat carries it; the
centre seat wins by far (+46% over uniform); and the V-shape sits
mid-field because it is only a partial concentration, 38% of its budget on
the centre. None of this supports "gate noise helps": the one-factor rows
above show added noise always hurts, and concentration wins by placing the
unavoidable budget where each unit does the least damage.

### One-factor controls (edges held fixed, only the mediator moves)

| Edges | mediator γ = 0.01 | mediator γ = 0.05 |
|-------|-------------------|-------------------|
| quiet [.01,.03,·,.03,.01] | 0.936 | 0.755 |
| loud [.05,.03,·,.03,.05] | 0.530 | 0.428 |

Monotone in both backgrounds: the mediator's own noise harms.

---

## Dynamical decoupling: it removes Σγ, and where it removes it matters

DD is simulated as a ×0.1 multiplier on the local γ of the treated sites,
so every DD row removes Σγ; the rows differ in how much and where. The
March "receiver is 3× more important than the mediator" compared a two-site
removal against a one-site removal. At MATCHED removal (one site, −0.045
of Σγ) the fair ranking:

| DD applied to | Σγ | MI | vs baseline |
|---------------|-----|------|------|
| none (baseline) | 0.250 | 0.338 | 0% |
| site 0 alone (sender end) | 0.205 | 0.459 | +36% |
| site 1 alone | 0.205 | 0.448 | +33% |
| mediator (site 2) | 0.205 | 0.428 | +27% |
| site 3 alone | 0.205 | 0.441 | +31% |
| site 4 alone (receiver end) | 0.205 | 0.459 | +36% |
| receiver pair (3+4) | 0.160 | 0.612 | +81% |
| mediator + receiver | 0.115 | 0.784 | +132% |
| everywhere | 0.025 | 1.535 | +355% |

At matched one-site removal the ranking is positional, not functional: it
falls from the ends toward the centre, the two ends tie to within 0.1%,
and the small left-right asymmetry at sites 1 and 3 is the initial state's
Bell-pair-versus-|+⟩|+⟩ imbalance. Sender and receiver are equally
profitable places to remove noise; the March "receiver first" reading came
from never measuring the sender side. The consistent picture with the
static profiles: removing a unit helps most at the ends and least at the
centre, which is the same fact as the centre being the cheapest seat to
LEAVE a unit on.

### Hardware translation (IBM Torino)

On real hardware the assignment can be implemented passively:
- Select **high-T2\* qubits** for the chain ENDS; sender and receiver
  sides profit equally
- The centre tolerates a **low-T2\* qubit** best, but its noise still
  costs; nothing is gained by making it noisy
- Apply **CPMG pulse sequences** (a standard dynamical decoupling
  technique using evenly spaced spin-echo pulses) on the end qubits first
- Toggle DD on/off for the **relay protocol** (staged γ switching,
  see [Relay Protocol](RELAY_PROTOCOL.md))

---

## What Does Not Work

### AC modulation on the gate (no effect within 0.3%)

Sinusoidal modulation of γ_M(t) = γ₀(1 + A·sin(2πft)) was tested at
frequencies f = 0.1 to 8.0 (spanning the palindromic mode frequencies
and the system's natural oscillation frequencies, called Bohr
frequencies). All frequencies produce MI within 0.3%
of the unmodulated baseline.

**Why it fails:** The palindromic mode structure does not couple to γ
modulation. The decay rates are set by the *time-averaged* γ, not by its
instantaneous value. AC modulation averages to the same mean γ and produces
no net effect. Only the **DC component** (the spatial profile) matters.

### State-dependent feedback (−3%, harmful)

Making γ_M depend on the system's own coherence:
γ_eff = γ_base × (1 + κ|⟨Z_i Z_j⟩|), tested at κ = 0.1 to 1.0.

All κ values slightly degrade performance. The feedback increases γ_M when
the system is coherent, the opposite of what helps. Negative feedback
(reducing γ_M when coherent) was not tested and might work. This remains
open.

---

## The Key Discovery: Time-Resolved γ Detection

### The first sign that γ is readable

The March run doubled γ on qubit 2 at t=5 and reported a detectable change
in the MI trajectory with "~0.5 time-unit resolution", from a single
trajectory with no control. The control, run 2026-08-23 (identical
evolution, γ₂ left at 0.05), with a lag scan past the switch: the two MI
trajectories differ by 0.0017 at half a unit, 0.0085 at one, and at most
0.0127 at any lag out to ten, against an intrinsic sample-to-sample swing
of ~0.1 on the same half-unit grid in that window. No lag gives clean separation from a single
trajectory; the March resolution figure was the run's own sampling grid
read back as a result, and no resolution figure survives the control.

What stands is weaker and still true: doubling γ₂ leaves a small,
systematic imprint, at most 0.013 in absolute MI at any lag (up to ~14%
of the MI value at the same lag), readable with a precise measurement
against a known baseline. That is
exactly the properly instrumented form the question took in
[γ as Signal](GAMMA_AS_SIGNAL.md): template matching over the whole
trajectory, 15.5 bits, 100% classification.

This observation directly motivated the
[γ as Signal](GAMMA_AS_SIGNAL.md) experiment, which formalized the question:
"If γ changes are detectable, how much information does the γ profile
carry?" The answer: **15.5 bits** of theoretical channel capacity at 1%
measurement noise, with 5 independent spatial modes. The dephasing rate is
not noise but a high-bandwidth information channel from outside to
inside.

---

## Design Rules

1. **Within a shape, minimize Σγ.** Monotone in every family tested.
   Across shapes it is NOT the ranking variable: a concentrated 0.25
   beats a uniform 0.13.

2. **At fixed Σγ, concentrate at the centre.** Concentration beats
   spreading at every seat tried (+4% to +46% over uniform at Σγ = 0.13),
   and the centre seat is the cheapest because sensitivity to noise falls
   from the ends toward the centre. The V-shape is a partial form of this
   rule, worth +6%. Noise on the centre still harms; it harms least
   there.

3. **DD at the ends first, centre last.** At matched Σγ removal both
   chain ends return +36% against the centre's +27%; sender and receiver
   ends tie.

4. **No AC modulation.** The palindromic structure does not couple to
   time-varying γ. Only DC spatial shaping matters.

5. **No positive feedback.** Self-regulating γ based on correlations is
   harmful. The system cannot improve its own channel.

6. **Single-trajectory γ monitoring is weak.** A doubled γ₂ imprints at
   most 0.013 in absolute MI at any lag; for channel estimation use the
   template-matching channel analysis of
   [γ as Signal](GAMMA_AS_SIGNAL.md), not a single-trajectory decay read.

### A naming caution

An inversion, not an ambiguity: plotted over
position, γ = [.01,.03,.05,.03,.01] rises to the center and falls again, a
peak (Λ), yet this page's March label for it is "V-shape", and its label
for the center-dipping [.05,.03,.01,.03,.05], which IS a geometric V, is
"inverse V". `docs/GLOSSARY.md`'s "V-shape profile" (edges higher than
center) is the geometrically correct usage, so the glossary and this page
use the same word for opposite objects and this page's labels are the
backwards ones. They are kept here unflipped because every citing document
refers to them; a profile is stated as its five numbers wherever it
matters, and the repo-wide unification is the follow-up arc's job (the one
pre-existing accurate cross-reading is the correction note in
[the Resonant Return experiment](RESONANT_RETURN.md), before its tests).

---

## Connection to the Broader Framework

**Upstream:** F91 of [`docs/ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md)
states that spectra depend on the ASSIGNMENT of a parameter profile, never
on its bare multiset; the assignment residual measured here is an effect of
the same kind on the γ axis, a dynamical observable rather than F91's block
spectra, so kinship and not derivation, and the repaired comparison is the
cautionary case where not even the multiset's sum was held fixed.
`hypotheses/THE_OTHER_SIDE.md` §22 and the arc `the_gate_that_does_not_gate`
(`compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs`) carry the
mediator-noise polarity question at N=5 and N=11; the controls here
confirm the N=5 polarity datum that arc already carries (quiet mediator
opens, loud closes), an input to its NextStep (1) rather than a
resolution of it.

**Downstream:** The γ-as-signal result
([GAMMA_AS_SIGNAL.md](GAMMA_AS_SIGNAL.md)) shows the profile is readable
from inside; the V-shape's role there is as one of the distinguishable
antenna shapes, a statement about reading the channel, not about transfer.

**The word (translation added 2026-08-23):** "noise" is the engineering
stance's word, painted true where γ degrades an intended computation, and
this page keeps it because every document citing it lives on that canvas.
The repo's own label for γ is different: γ is light, the sending; the
chain is not looked at, it stands in illumination
([the label map](../docs/quantum/THE_LABEL_MAP.md), arc
`gamma_is_the_sender_not_the_watching`), and this page's own sibling
experiment, published the same March day, had already repainted the label
from below: "This is not noise. This is a channel."
([γ as Signal](GAMMA_AS_SIGNAL.md)). Read on that canvas,
the two levers are light and dark as ingredients: keep the total
illumination low, and where light must fall, let it fall on the centre and
keep the ends dark; a γ profile is a pattern of light, and this page's
assignment axis is the shape of that pattern. The translation matters for
search as much as for meaning: a sweep for "noise" and a sweep for "light"
return different halves of this repository, the same mechanism as the
confound this page repairs, two vocabularies that cannot find each other.

---

## Reproducibility

| Script | What it computes | Runtime |
|--------|-----------------|---------|
| [using_gamma.py](../simulations/using_gamma.py) | All strategies incl. the gate-open/closed one-factor rows + time-resolved decoder | ~10 min |
| [using_gamma.txt](../simulations/results/using_gamma.txt) | Full numerical results | - |

The Σγ-matched and matched-removal control numbers in this page were
computed 2026-08-23 through the committed primitives in
[mediator_bridge.py](../simulations/mediator_bridge.py)
(`make_initial_state('bell_A_0M_pp_B')`, `build_mediator_system`,
`evolve_rho` to t = 5.0, `mutual_info(ρ, 5, [0,1], [3,4])`), independently
in two sessions with matching results (one rounding-route deviation,
resolved; the catch entry records it), and are gated by
[gamma_control_two_lever_check.py](../simulations/gamma_control_two_lever_check.py)
(34 gates: every table row, the orderings, the lever crossing, and the
switch-control bound; ~2 min).

All scripts use NumPy and SciPy. No proprietary dependencies.

---

## References

- [γ as Signal](GAMMA_AS_SIGNAL.md): the γ profile is not just shapeable --
  it is **readable**. 100% classification, 15.5 bits channel capacity.
  This experiment was the precursor to that discovery.
- [Relay Protocol](RELAY_PROTOCOL.md): staged γ switching
- [Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md):
  CΨ = 1/4 as threshold voltage, γ_M as gate signal; its polarity (quiet
  gate = open channel) is the one the controls here confirm
- [Scaling Curve](SCALING_CURVE.md): MI vs chain length baseline
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the
  palindromic spectral structure of the chain
- `docs/CAUGHT_ERRORS.md`, 2026-08-23 entry: the Σγ confound, its shape,
  and how it was caught
