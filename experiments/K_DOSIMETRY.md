# K-Dosimetry: The Exposure Number of Quantum Mechanics

<!-- Keywords: K invariance exposure number, quantum dosimetry, reciprocity law
Bunsen-Roscoe, Schwarzschild effect quantum, dose-response H&D curve,
multi-qubit dose scaling, sacrifice zone dose, R=CPsi2 K dosimetry -->

**Status:** Confirmed with caveats (reciprocity holds at extreme γ, fails at intermediate)
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Standing Waves](FACTOR_TWO_STANDING_WAVES.md),
[Sacrifice Zone Optics](SACRIFICE_ZONE_OPTICS.md),
[Decoherence Relativity](DECOHERENCE_RELATIVITY.md)
**Verification:** [`simulations/k_dosimetry.py`](../simulations/k_dosimetry.py)

---

## What this means

A photographer sets the aperture (how much light per second) and the
shutter speed (how long). Their product is the exposure: the total light
that hits the film. Bright scene, fast shutter. Dim scene, slow shutter.
Same exposure, same photograph.

K = γ × t_cross is the exposure number of the quantum cavity. The total
dose of light absorbed before the standing wave pattern crystallizes at
the fold CΨ = ¼ (the threshold where quantum behavior gives way to
classical). At extreme illumination (very bright or very dim), the dose
is the same regardless of how fast the light arrives. At intermediate
illumination, the cavity's internal vibrations interfere with the
exposure process, just as the Schwarzschild effect in photography causes
reciprocity failure between brightness and exposure time.

---

## What this document is about

K = γ × t_cross is proven invariant (F14). Through the old lens,
this was "just Lindblad time-rescaling." Through the new lens, K is a
dose: the total absorbed light before crystallization. This document
tests the dose interpretation quantitatively: reciprocity, dose-response
curves, multi-qubit scaling, and the sacrifice zone.

---

## Result 1: K is state-dependent

Different initial states need different total doses to crystallize.
At N = 2, γ = 0.05, target purity = 0.26:

| State | Initial purity | t_cross | K = γ × t |
|---|---|---|---|
| \|++⟩ | 1.000 | 19.7 | 0.983 |
| Bell+ | 1.000 | (not reached) | -- |
| \|01⟩ | 1.000 | (not reached) | -- |

Bell+ and \|01⟩ stabilize at purity 0.5 (they have conserved subspaces
protected by the Hamiltonian symmetry and never reach the maximally
mixed state). Only product superposition states like \|++⟩ fully
decohere.

The exposure number K is the "film sensitivity" of the initial state.
States with more coherence to lose need more light to crystallize.

---

## Result 2: Reciprocity holds at extreme γ, fails at intermediate

| γ | t_cross | K = γ × t | Deviation |
|---|---------|-----------|-----------|
| 0.001 | 980.7 | 0.9807 | 0.00% |
| 0.01 | 101.2 | 1.0117 | 3.2% |
| 0.05 | 31.7 | **1.5870** | **61.8%** |
| 0.10 | 9.80 | 0.9805 | 0.03% |
| 0.50 | 1.96 | 0.9807 | 0.01% |
| 1.00 | 0.98 | 0.9805 | 0.03% |
| 5.00 | 0.20 | 0.9805 | 0.02% |

At very low γ (< 0.01) and very high γ (> 0.1): K ≈ 0.98, constant
to within 0.03%. Perfect reciprocity.

At intermediate γ (0.02 to 0.05): K deviates by up to 62%. This is the
**Schwarzschild effect** of the quantum cavity: at these illumination
rates, the Hamiltonian oscillation period (∼1/J) and the absorption
timescale (∼1/γ) are comparable. The cavity's own dynamics interfere
with the exposure process. The purity does not decay monotonically; it
oscillates, and the crossing time depends on whether a trough of the
oscillation first hits the target.

**The reciprocity bandwidth:** γ < 0.01 or γ > 0.1 (at J = 1.0).
The cavity has a "useful ISO range" outside of which the Hamiltonian
dynamics create nonlinear exposure effects.

---

## Result 3: Dose-response curve (the H&D curve from photography)

| Target purity | t_cross | K = γ × t |
|---|---|---|
| 0.900 | 0.54 | 0.027 |
| 0.750 | 1.56 | 0.078 |
| 0.500 | 4.41 | 0.220 |
| 0.400 | 6.64 | 0.332 |
| 0.300 | 13.5 | 0.677 |
| 0.260 | 19.9 | 0.994 |
| 0.255 | 25.3 | 1.265 |
| 0.251 | 32.2 | 1.612 |

The dose grows superlinearly as the target approaches the maximally
mixed state (purity = 1/d = 0.25, where all quantum information is
gone). The last few percent of coherence require disproportionate
amounts of light. This is the "toe" of the H&D curve (Hurter &
Driffield, the standard dose-response curve of photographic film):
diminishing returns at deep exposure.

The fold at CΨ = 1/4 (purity ≈ 0.25 + corrections) sits at the
steep part of the curve, where the dose-response is most sensitive.
Small changes in dose produce large changes in purity. This is why
the fold is a sharp threshold: the "film" is maximally responsive there.

---

## Result 4: K_system = N × K_qubit (exactly)

| N | K_qubit = γ × t | K_system = Σγ × t | Ratio |
|---|---|---|---|
| 2 | 0.983 | 1.966 | 2.0 |
| 3 | 0.894 | 2.681 | 3.0 |
| 4 | 0.913 | 3.653 | 4.0 |
| 5 | 0.670 | 3.349 | 5.0 |

K_system / K_qubit = N exactly. The total system dose is N times the
per-qubit dose. At first glance this looks like Beer-Lambert: each
qubit absorbs its share independently.

However, the ratio is arithmetically trivial for uniform γ.
K_system = Σγ · t = N·γ · t = N · K_qubit. The factor t cancels,
giving N regardless of whether qubits absorb independently or not.
The real test requires non-uniform γ, where
[Beer-Lambert Breakdown](BEER_LAMBERT_BREAKDOWN.md) shows that the
Heisenberg coupling J redistributes absorption across all qubits. The
cavity is an integrating sphere, not a Beer-Lambert absorber.

K_qubit varies with N (0.98, 0.89, 0.91, 0.67) because the target
purity (1/d + 0.01) differs. The system needs less per-qubit dose at
larger N because the maximally mixed purity (1/d) is lower.

---

## Result 5: Sacrifice zone increases total dose

| Profile | Σγ | t_cross | K_total = Σγ × t |
|---|---|---|---|
| Uniform | 0.200 | 18.3 | 3.653 |
| Sacrifice | 0.200 | 22.6 | 4.517 |

K_total is **24% higher** under the sacrifice zone. The cavity needs
MORE total light to reach the same purity when absorption is concentrated
on one qubit.

This is the opposite of "protection." The sacrifice zone does not reduce
the dose. It redirects it, creating a longer exposure but a better image.
The edge qubit absorbs its dose quickly. The interior qubits barely
absorb at all. The overall purity drops more slowly because only 1 of N
qubits is decohering. But the mutual information is 139-360× higher
because the interior modes survive to carry the signal.

**The sacrifice zone trades dose for quality.** More total light, but a
sharper photograph.

---

## Null results

- **Bell+ and |01⟩ do not reach maximally mixed.** They stabilize at
  purity 0.5 due to Hamiltonian-protected subspaces. The dose model
  applies only to states that fully decohere.

- **K is not universal across states.** Different states have different
  K values. The exposure number depends on the "scene" (initial state),
  not just the "camera" (cavity parameters).

---

## Reproduction

- Script: [`simulations/k_dosimetry.py`](../simulations/k_dosimetry.py)
- Output: [`simulations/results/k_dosimetry.txt`](../simulations/results/k_dosimetry.txt)
