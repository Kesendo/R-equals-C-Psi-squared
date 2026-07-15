# F130 on Hardware: The Beat Protocol Is Not Flyable

**Status:** Verified infeasibility (design-stage null result; no QPU spent)
**Date:** 2026-07-15
**Authors:** Thomas Wicht & Claude
**Gate:** [`simulations/f130_hw_infeasibility.py`](../simulations/f130_hw_infeasibility.py)
**Relation:** the law this protocol tried to fly is [F130 in the time domain](F130_TIME_DOMAIN_DECOUPLING.md); what survives is the F129 Ramsey-fringe design in [IBM_F129_RAMSEY_FRINGE.md](IBM_F129_RAMSEY_FRINGE.md).

## What this document is about

F130 says that two mode triples with equal levels never talk to each other under
dephasing: the second-order coupling between their coherence multiplets is exactly
zero, and the time-domain gate measured what replaces it, an effective coupling of
5.9·10⁻⁵ q⁴ (protected) against 1.9·10⁻² q² (generic), q = J/γ. The obvious hardware
experiment is "the beat": prepare coherence weight in one multiplet and watch the q²
beat appear on a generic pair and stay absent on the collision pair. This document
records, before any runner was built and before any QPU was spent, that the beat
protocol cannot fly on Heron-class hardware. The verdict is a design-stage null
result: negative results matter, and this one redirects the hardware question to a
first-order observable that does work.

## What this document settles

1. **Trotterization does not kill the degeneracy.** The exact level collision
   survives a first-order Trotter step (odd/even Givens layering) to third order:
   the Floquet splitting of the n = 12 pair (1,2,10) ~ (3,5,6) is 0.13·θ³ per step
   (measured 0.1326-0.1336 over θ ∈ [0.02, 0.32], gate part A). At the safe angle
   θ = 0.32 that is 4.4·10⁻³ rad per step, under 9 % of the engineered dephasing
   linewidth γ_step = 0.05 the concentrator pipeline certifies. The protocol does
   not die on the degeneracy.

2. **It dies on amplitude, everywhere the beat exists.** An upper bound on the beat
   signal that is optimistic in every disputable factor (gate part B) is

   A(q, M) ≤ ½ · M·g_step · e^(−6·M·γ_step) · (1 − p₂)^CX(M),

   with g_step = 1.9·10⁻² q² · γ_step the generic coupling per step (the committed
   coefficient from [`f130_time_domain_decoupling.py`](../simulations/f130_time_domain_decoupling.py) G2),
   ½ the maximal initial (1,2)-coherence weight, e^(−6γ_step) the rung −6 decay per
   step, and CX(M) the two-qubit count of the 3-magnon preparation network, M Trotter
   steps, and a cross-sector readout network at its theoretical minimum. (γ_step in
   this bound is the scan's free per-step dephasing knob, a different number from the
   concentrator pipeline's engineered 0.05 quoted in item 1.)

   The domain logic matters, because the bound has NO interior maximum: at the
   optimal M it reduces to amp ∝ q·e^(−1.92/q), monotonically growing in q, so its
   ceiling is set by the coupling law's validity, not by an optimum the scan finds.
   The q² law is second-order perturbation theory in q = J/γ (fit window q ≤ 0.32),
   meaningful up to q ≈ 1; beyond that the chain is hopping-dominated, the
   dephasing-induced beat B is no longer the object being measured, and extrapolating
   the law there does not describe any experiment. Inside the domain the best point
   (q = 1, γ_step = 0.32, M = 1) gives amplitude 1.5·10⁻⁴ and SNR 0.014 at 8192
   shots: **360× short of the 5σ discrimination bar**, roughly 2.6 orders in
   amplitude and over 5 orders in shots. Even the generous over-extrapolation to
   q = 2 (printed separately by the gate; 6× past the fit window) reaches only
   SNR 0.073, still 69× short. The protected arm is smaller by another factor
   q²·(5.9·10⁻⁵/1.9·10⁻²) ≈ 1.2·10⁻² and is invisible a fortiori. And the large-q
   region where the extrapolated bound would eventually cross the bar is exactly
   where the beat picture has dissolved into hopping-dominated first-order physics,
   i.e. the regime of the Ramsey pivot below.

The optimism list, spelled out because each item strengthens the verdict (the bound
keeps only the two unavoidable physical losses, the rung decay and the gate
attenuation; every disputable choice is resolved in the protocol's favor):
linear-in-t coherent growth with no detuning suppression (the gate's G5 says the
self-blocks split at q², which only reduces transfer); no T1/T2* idle decay; no
readout error; no Trotter-splitting dephasing; unit fringe contrast; and the q²
coupling law taken at face value up to the perturbative edge q = 1 (fit window
q ≤ 0.32), with the q = 2 over-extrapolation shown besides. One caution the list
does NOT contain: extrapolating the growing law to large q is the FEASIBLE
direction, not the conservative one; the domain cap is what carries the verdict,
and it is physics (second-order perturbation theory), not a knob. Constants:
p₂ = 0.3 % median two-qubit error and the coherent NN ZZ ≈ 4 kHz noise model of
the `price_pair_locality_marrakesh_july2026` confirmation.

## The reading: the smallness is the law

This is not an instrument falling short of a signal; the absence of a signal IS the
claim. F130's content is that the q² beat, the generic second-order channel, is
switched off exactly on the collision locus. A protocol that needs to resolve the q⁴
remainder is asking the law to violate itself loudly enough to measure. The
experiment's honest form is therefore not "see the q⁴", it is "see a first-order
quantity that distinguishes equal levels from detuned ones", and that observable
exists: the collision as a standing Ramsey fringe, where the phase slope of a
two-triple superposition is first-order in the level difference and the F129
collision pins it to zero.

## What survives

- **The q⁴-vs-q² law itself** is untouched; it stays verified from below
  ([F130_TIME_DOMAIN_DECOUPLING.md](F130_TIME_DOMAIN_DECOUPLING.md), gate G1-G5).
- **The Trotter-survival fact** (0.13·θ³ per step) is load-bearing for the successor:
  the collision remains a collision on the stepped device, so a fringe experiment
  inherits an exactly standing fringe up to a computable third-order drift.
- **The pivot:** the F129 Ramsey-fringe design,
  [IBM_F129_RAMSEY_FRINGE.md](IBM_F129_RAMSEY_FRINGE.md), a first-order protocol at
  the n = 9 clean mirror collision (impostor rejection 21.7σ, verdict control 50.2σ,
  6.6-7.7 QPU minutes projected). F130's own n = 12 pair does NOT fly there: its nearest impostor
  sits at dS = 0.018 and is inseparable at this budget; the pivot carries F129's law,
  not F130's pair.

## Reproduction

```bash
python simulations/f130_hw_infeasibility.py   # part A: Trotter splitting; part B: the optimistic bound
```

Part A output (per step): splitting 1.06·10⁻⁶ at θ = 0.02 rising as θ³ to 4.38·10⁻³
at θ = 0.32, ratio split/θ³ stable at 0.133; the detuned control (3,4,9) winds at
first order (1.9·10⁻² at θ = 0.02). Part B output: best in-domain point q = 1,
γ_step = 0.167, M = 1, amplitude 2.005·10⁻⁴, SNR 0.0181 (276× short of 5σ); the
q = 2 over-extrapolation at amplitude 8.012·10⁻⁴, SNR 0.073 (69× short); verdict
INFEASIBLE.
