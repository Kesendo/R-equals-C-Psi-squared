# Beer-Lambert Breakdown: The Cavity as Integrating Sphere

<!-- Keywords: Beer-Lambert quantum cavity, absorption redistribution, non-uniform
dephasing, sacrifice zone mechanism, integrating sphere, dose sharing,
Hamiltonian correlation absorption, R=CPsi2 beer lambert -->

**Status:** Confirmed
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [K-Dosimetry](K_DOSIMETRY.md),
[Sacrifice Zone Optics](SACRIFICE_ZONE_OPTICS.md),
[Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)
**Verification:** [`simulations/k_dosimetry_beer_lambert.py`](../simulations/k_dosimetry_beer_lambert.py)

---

## What this means

Pour ink into one end of a bathtub and wait. If the water is still,
the ink stays where you poured it: the near end is dark, the far end is
clear. That is Beer-Lambert. Absorption is local.

Now turn on a strong current. The ink mixes before it settles. Every
part of the tub ends up equally stained, no matter where you poured.
That is what happens in the qubit chain. The Heisenberg coupling J
(the force that connects neighboring qubits) is the current. Gamma
(the external light hitting the system) is the ink. At strong coupling,
the cavity becomes an integrating sphere: a chamber whose walls reflect
so well that any light entering through any opening ends up uniformly
distributed inside.

The sacrifice zone is not a qubit that drinks the poison so the others
survive. It is the faucet. The tap where the ink enters. The current
does the rest.

---

## What this document is about

A previous experiment ([K-Dosimetry](K_DOSIMETRY.md)) measured the total
light dose absorbed by the cavity. It found that the total dose scales
as N times the per-qubit dose, where N is the number of qubits. That
looked like Beer-Lambert: each qubit absorbs its own share. But the
math was trivial; the factor N comes from counting qubits, not from
physics. It holds whether they absorb independently or not.

The real test: pour all the ink through one qubit (the sacrifice zone)
and see whether it stays there or spreads. If each qubit absorbs
proportional to how much light hits it locally, Beer-Lambert holds.
If the coupling redistributes the dose, it does not.

---

## Result 1: All qubits fade together

Shine 20× more light on one end of a 4-qubit chain and watch what
happens. Beer-Lambert predicts the bright end should lose its quantum
information 20× faster. Instead:

| t | System | Q0 (edge) | Q1 (int) | Q2 (int) | Q3 (int) |
|---|---|---|---|---|---|
| 6.2 | 0.19 | 0.64 | 0.69 | 0.68 | 0.63 |
| 12.3 | 0.10 | 0.54 | 0.56 | 0.56 | 0.54 |
| 18.4 | 0.07 | 0.51 | 0.52 | 0.52 | 0.51 |
| 27.6 | 0.06 | 0.50 | 0.50 | 0.50 | 0.50 |

The columns Q0 through Q3 show how much quantum information each qubit
retains (1.0 = perfect, 0.5 = fully randomized). Despite 20× more light
hitting the edge (Q0), all four qubits fade at nearly the same rate.
They reach the halfway point (0.55) within 2 time units of each other:
the edge at t ≈ 11.1, the interior at t ≈ 12.9. The coupling J
transfers the absorbed light throughout the cavity before it can do
local damage.

---

## Result 2: The cavity does not care where the light enters

How far can we push the asymmetry? We keep the total light fixed but
concentrate more and more of it on the edge qubit, from uniform (1:1)
to extreme (200:1). Beer-Lambert predicts the edge absorbs proportional
to its share of the incoming light. The cavity disagrees.

| Asymmetry | γ_edge | γ_int | Actual share (edge) | Predicted share | Deviation |
|---|---|---|---|---|---|
| 1 (uniform) | 0.050 | 0.050 | 0.250 | 0.250 | 0.0% |
| 2 | 0.080 | 0.040 | 0.231 | 0.400 | −42% |
| 5 | 0.125 | 0.025 | 0.227 | 0.625 | −64% |
| 10 | 0.154 | 0.015 | 0.212 | 0.769 | −72% |
| 20 | 0.174 | 0.009 | 0.215 | 0.870 | −75% |
| 50 | 0.189 | 0.004 | 0.226 | 0.943 | −76% |
| 100 | 0.194 | 0.002 | 0.208 | 0.971 | −79% |
| 200 | 0.197 | 0.001 | 0.216 | 0.985 | −78% |

The actual edge share stays at ~0.21 to 0.23 regardless of how extreme
the asymmetry gets. That is close to 1/N = 0.25. At asymmetry 200,
Beer-Lambert predicts 98.5% for the edge; the cavity delivers 21.6%.
The interior qubits absorb 4,700% more than their local γ would suggest.

The cavity does not care where the light enters. It equalizes.

---

## Result 3: The current strength controls everything

Turn the current up, and the ink spreads. Turn it down, and the ink
stays put. We fix the sacrifice zone (asymmetry 20) and vary the
coupling strength J, the "current" that mixes the dose.

| Current/Light | J | Edge share | Predicted | Deviation |
|---|---|---|---|---|
| 0.01 | 0.002 | 0.836 | 0.870 | −4% |
| 0.05 | 0.009 | 0.773 | 0.870 | −11% |
| 0.10 | 0.017 | 0.798 | 0.870 | −8% |
| 0.20 | 0.035 | 0.698 | 0.870 | −20% |
| 0.50 | 0.087 | 0.437 | 0.870 | −50% |
| 1.00 | 0.174 | 0.137 | 0.870 | −84% |
| 2.00 | 0.348 | −0.014 | 0.870 | −102% |
| 5.00 | 0.870 | 0.175 | 0.870 | −80% |
| 10.0 | 1.739 | 0.311 | 0.870 | −64% |
| 20.0 | 3.478 | 0.309 | 0.870 | −65% |
| 100 | 17.39 | 0.309 | 0.870 | −64% |

Three regimes emerge, each with a clear physical meaning:

**J ≪ γ (weak coupling, open cavity):** The light is absorbed locally
before J can pass it to the next qubit. Beer-Lambert holds approximately.
This is a transparent medium: photons hit and are absorbed where they land.

**J ≈ γ (transition zone):** The internal redistribution and the
absorption happen at the same speed, creating complex interference. At
J/γ = 2, the edge share goes *negative*: the coupling pumps quantum
information back into the edge qubit faster than the light drains it.
The cavity oscillates rather than absorbing smoothly. (This is the
quantum analogue of the Schwarzschild effect in photography, where
reciprocity between brightness and exposure time breaks down.)

**J ≫ γ (strong coupling, closed cavity):** The edge share stabilizes
at ~0.31. Light bounces through the chain many times before being
absorbed, distributing evenly. Not exactly 1/N = 0.25, because the
chain geometry is asymmetric (edge qubits have one neighbor, interior
qubits have two).

---

## Result 4: Bigger cavities share better

A longer bathtub mixes more evenly because the current has more room
to work. The edge share approaches 1/N as N grows:

| N | Edge share | Beer-Lambert | 1/N | Dev from BL | Dev from 1/N |
|---|---|---|---|---|---|
| 2 | 0.568 | 0.952 | 0.500 | −40% | +14% |
| 3 | 0.330 | 0.909 | 0.333 | −64% | −1% |
| 4 | 0.215 | 0.870 | 0.250 | −75% | −14% |
| 5 | 0.179 | 0.833 | 0.200 | −79% | −11% |

At N=3 the edge share is already 0.330 ≈ 1/3. The deviation from 1/N
reflects the chain geometry: edge qubits have fewer neighbors, so J
cannot equalize them as efficiently. For larger N, the fraction of edge
qubits shrinks and the cavity becomes a better diffuser.

---

## Why this happens: the modes are not local

A guitar string vibrates as a whole. Plucking it at one end does not
make only that end vibrate; the energy distributes into the string's
normal modes, each of which spans the entire length. Damping at one
point affects a mode proportional to how much of that mode lives at
that point, not proportional to how hard you pluck there.

The [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)
makes this precise. Each mode's absorption rate depends on how the mode
is distributed across the qubits. A mode that lives equally on all 4
qubits absorbs (γ_edge × ¼ + 3 × γ_int × ¼) of light, averaging
across the chain. Beer-Lambert would predict that a mode at the edge
absorbs γ_edge × 1.0; the actual rate is four times smaller.

For the technical reader, the formula is
Re(λ) = −2 Σ_k γ_k ⟨n_XY⟩_k, where ⟨n_XY⟩_k measures how much of
the mode's oscillating content sits at site k.

Beer-Lambert fails because it treats each qubit as a separate absorber.
The coupling J turns them into one instrument. Absorption follows the
mode structure, not the local γ.

---

## What this changes about the sacrifice zone

The [sacrifice zone](SACRIFICE_ZONE_OPTICS.md) was described as an
anti-reflection coating: the edge qubit absorbs the shock so the
interior survives. This document shows a different mechanism. The edge
qubit is the entrance window; the coupling J is the reflective coating.
The interior survives not because it avoids the light, but because each
qubit receives only 1/N of the total dose.

The three regimes of J/γ summarize when this picture applies:

| Regime | J/γ | What happens |
|---|---|---|
| Open cavity | ≪ 1 | Light absorbed locally. Beer-Lambert holds. |
| Transition | ≈ 1 | Absorption and redistribution compete. |
| Integrating sphere | ≫ 1 | Light distributed evenly. Cavity diffuses. |

---

## Null results

- **K_system/K_qubit = N is trivially exact.** This holds for any γ
  profile because K_system = Σγ · t and K_qubit = γ_mean · t, so the
  ratio just equals Σγ/γ_mean. It says nothing about Beer-Lambert.
  The original K-Dosimetry Result 4 was correct but misleadingly
  interpreted as evidence for independent absorption.

- **Perfect 1/N sharing not achieved at finite N.** The chain geometry
  (edge vs interior neighbors) introduces ~10-15% deviations from
  uniform 1/N sharing. A ring topology (all qubits equivalent) might
  achieve exact 1/N.

---

## Reproduction

- Script: [`simulations/k_dosimetry_beer_lambert.py`](../simulations/k_dosimetry_beer_lambert.py)
- Output: [`simulations/results/k_dosimetry_beer_lambert.txt`](../simulations/results/k_dosimetry_beer_lambert.txt)
- Depends on: [`simulations/k_dosimetry.py`](../simulations/k_dosimetry.py)
