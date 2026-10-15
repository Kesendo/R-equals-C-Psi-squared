# Beer-Lambert Breakdown: The Cavity as Integrating Sphere

<!-- Keywords: Beer-Lambert quantum cavity, absorption redistribution, non-uniform
dephasing, sacrifice zone mechanism, integrating sphere, dose sharing,
Hamiltonian correlation absorption, R=CPsi2 beer lambert -->

**Status:** Confirmed (Beer-Lambert fails; at strong coupling the lit end shares its excess with its mirror image, and the share is not even)
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [K-Dosimetry](K_DOSIMETRY.md),
[Concentrator Optics](CONCENTRATOR_OPTICS.md),
[Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)
**Verification:** [`simulations/k_dosimetry_beer_lambert.py`](../simulations/k_dosimetry_beer_lambert.py)

---

## What this means

Pour ink into one end of a bathtub and wait. If the water is still,
the ink stays where you poured it: the near end is dark, the far end is
clear. That is Beer-Lambert. Absorption is local.

Now turn on a current. The ink mixes as it settles, and the stronger the
current, the less the stain remembers where it was poured. That is what
happens in the qubit chain. The Heisenberg coupling J (the force that
connects neighboring qubits) is the current. Gamma (the external light
hitting the system) is the ink. At strong coupling the cavity acts like an
integrating sphere only in part: a chamber whose walls reflect so well that
the light no longer stays at its opening. It does not spread the light
evenly, though. The lit end keeps more than its even share, and its mirror
image, the dark far end, keeps exactly as much; the middle keeps less. At
moderate coupling the lit end alone still stands out.

The sacrifice zone is not a qubit that drinks the poison so the others
survive. It is the faucet. The tap where the ink enters. The current
carries part of it on.

---

## What this document is about

A previous experiment ([K-Dosimetry](K_DOSIMETRY.md)) measured the total
light dose absorbed by the cavity. It found that the total dose scales
as N times the per-qubit dose, where N is the number of qubits. That
looked like Beer-Lambert: each qubit absorbs its own share. But the
math was trivial; the factor N comes from counting qubits, not from
physics. It holds whether they absorb independently or not.

The real test: pour most of the ink through one qubit (the sacrifice
zone) and see whether it stays there or spreads. If each qubit absorbs
proportional to how much light hits it locally, Beer-Lambert holds.
If the coupling redistributes the dose, it does not.

Every run starts from |+⟩^⊗N on the Heisenberg chain and keeps the total
light fixed, Σγ = 0.2. A qubit's **share** is its purity loss over the
first time unit (t from 0.01 to about 0.97), divided by the sum of all
qubits' losses; Beer-Lambert predicts γ_k/Σγ.

---

## Result 1: All qubits fade over a similar time

Shine 20× more light on one end of a 4-qubit chain (J = 1) and watch
what happens. Beer-Lambert predicts the bright end should lose its
quantum information 20× faster. Instead:

| t | System | Q0 (edge) | Q1 (int) | Q2 (int) | Q3 (int) |
|---|---|---|---|---|---|
| 6.2 | 0.19 | 0.63 | 0.67 | 0.68 | 0.64 |
| 12.3 | 0.10 | 0.54 | 0.56 | 0.56 | 0.54 |
| 18.4 | 0.07 | 0.51 | 0.52 | 0.52 | 0.51 |
| 27.6 | 0.06 | 0.50 | 0.50 | 0.50 | 0.50 |

The columns Q0 through Q3 show how much quantum information each qubit
retains (1.0 = perfect, 0.5 = fully randomized). Despite 20× more light
hitting the edge (Q0), the four qubits reach purity 0.55 within 2.3 time
units of each other: the edge at t = 11.03, the far end Q3 at 11.36, the
two middle qubits at 13.13 and 13.30. The two ends fade first, the
lit one and the dark one alike; the coupling J carries the absorbed light
along the chain before it can do purely local damage.

---

## Result 2: The edge keeps a large share at moderate coupling

How far can we push the asymmetry? We keep the total light fixed but
concentrate more and more of it on the edge qubit, from uniform (1:1)
to extreme (200:1), at J = 1. Beer-Lambert predicts the edge absorbs
proportional to its share of the incoming light. The cavity disagrees.

| Asymmetry | γ_edge | γ_int | Actual share (edge) | Predicted share | Deviation |
|---|---|---|---|---|---|
| 1 (uniform) | 0.050 | 0.050 | 0.250 | 0.250 | 0.0% |
| 2 | 0.080 | 0.040 | 0.292 | 0.400 | −27% |
| 5 | 0.125 | 0.025 | 0.355 | 0.625 | −43% |
| 10 | 0.154 | 0.015 | 0.397 | 0.769 | −48% |
| 20 | 0.174 | 0.009 | 0.427 | 0.870 | −51% |
| 50 | 0.189 | 0.004 | 0.451 | 0.943 | −52% |
| 100 | 0.194 | 0.002 | 0.460 | 0.971 | −53% |
| 200 | 0.197 | 0.001 | 0.464 | 0.985 | −53% |

The edge share rises with the asymmetry and levels off near 0.46: about
half of what Beer-Lambert predicts, and still far above the even share
1/N = 0.25. At asymmetry 200, Beer-Lambert predicts 98.5% for the edge;
the cavity delivers 46%. Interior qubit Q1's share of the purity loss is
then 0.29, about 5,800% above what its local γ would suggest; the light
itself is still taken up only where γ is (the Absorption Theorem below),
and the coupling J carries the loss inward.

At this coupling the cavity redistributes, but it does not equalize.

---

## Result 3: The current strength controls the share

Turn the current up, and the ink spreads. Turn it down, and the ink
stays put. We fix the sacrifice zone (asymmetry 20) and vary the
coupling strength J, the "current" that mixes the dose.

| Current/Light | J | Edge share | Predicted | Deviation |
|---|---|---|---|---|
| 0.01 | 0.002 | 0.831 | 0.870 | −4% |
| 0.05 | 0.009 | 0.831 | 0.870 | −4% |
| 0.10 | 0.017 | 0.831 | 0.870 | −4% |
| 0.20 | 0.035 | 0.830 | 0.870 | −5% |
| 0.50 | 0.087 | 0.823 | 0.870 | −5% |
| 1.00 | 0.174 | 0.799 | 0.870 | −8% |
| 2.00 | 0.348 | 0.718 | 0.870 | −17% |
| 5.00 | 0.870 | 0.465 | 0.870 | −47% |
| 10.0 | 1.739 | 0.331 | 0.870 | −62% |
| 20.0 | 3.478 | 0.262 | 0.870 | −70% |
| 50.0 | 8.696 | 0.275 | 0.870 | −68% |
| 100 | 17.39 | 0.290 | 0.870 | −67% |

The share falls smoothly from the weak-coupling value to its lowest scanned value, 0.262
at J/γ = 20, and rises after it to 0.29, where it stays up to J/γ = 3000, in
three regimes:

**J ≪ γ (weak coupling, open cavity):** The light is absorbed locally
before J can pass it to the next qubit. The edge share is 0.831, close to
Beer-Lambert's 0.870; the remaining 4% is the share's own nonlinearity,
since purity loss over a finite time is not proportional to γ.

**J ≈ γ (transition zone):** Absorption and redistribution compete, and
the edge share falls steadily, from 0.80 at J/γ = 1 to 0.33 at J/γ = 10.

**J ≫ γ (strong coupling, closed cavity):** From J/γ ≈ 100 on the four
shares are near their limit, and by J/γ = 1000 they settle at [0.293, 0.207, 0.207, 0.293] (J/γ = 1000 and 3000): the
dark far end, the lit end's mirror image, takes the same share as the lit
end to three digits, and the profile is not the even 1/N = 0.25; those two
take 17% more than that, the middle 17% less. The pair is the lit site and
its image, not the chain's ends as such: on a ring of four, where every site
has two neighbours, the lit site and the site opposite it take 0.34 each and
the two between them 0.16. With equal light on every site each share is
exactly 1/N, since the Hamiltonian then never touches |+⟩^⊗N
([K-Dosimetry](K_DOSIMETRY.md)).

---

## Result 4: Cavity size does not equalize at moderate coupling

At J = 1 and asymmetry 20, a longer chain does not bring the edge share
toward 1/N:

| N | Edge share | Beer-Lambert | 1/N | Dev from BL | Dev from 1/N |
|---|---|---|---|---|---|
| 2 | 0.428 | 0.952 | 0.500 | −55% | −14% |
| 3 | 0.466 | 0.909 | 0.333 | −49% | +40% |
| 4 | 0.427 | 0.870 | 0.250 | −51% | +71% |
| 5 | 0.408 | 0.833 | 0.200 | −51% | +104% |

The edge keeps about 0.41 to 0.47 at every N from 2 to 5, so its distance
from 1/N grows with N while its distance from Beer-Lambert stays near
50%. At this coupling the length of the tub does not do the mixing; the
current's strength does (Result 3).

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
The coupling J turns them into one instrument, and the stronger the
coupling, the more the absorption follows the mode structure rather than
the local γ.

---

## What this changes about the sacrifice zone

The [sacrifice zone](CONCENTRATOR_OPTICS.md) was described as an
anti-reflection coating: the edge qubit absorbs the shock so the
interior survives. This document shows a different mechanism. The edge
qubit is the entrance window; the coupling J is the reflective coating.
The interior survives not because it avoids the light, but because the
coupling spreads the dose: at strong coupling into a pattern where the lit
edge and its mirror image keep 0.29 each and the middle 0.21 (four sites),
at moderate coupling with the lit edge still keeping nearly half.

The three regimes of J/γ summarize when this picture applies:

| Regime | J/γ | What happens |
|---|---|---|
| Open cavity | ≪ 1 | Light absorbed locally. Beer-Lambert holds up to the share's 4% nonlinearity. |
| Transition | ≈ 1 to 10 | Absorption and redistribution compete; the edge share falls from 0.80 to 0.33. |
| Integrating sphere | ≳ 100 | The shares approach, and by J/γ = 1000 reach, 0.29 for the lit edge and its mirror image, 0.21 for the middle (N = 4); J/γ ≈ 20 is the dip on the way. |

---

## Null results

- **K_system/K_qubit = N is trivially exact.** This holds for any γ
  profile because K_system = Σγ · t and K_qubit = γ_mean · t, so the
  ratio just equals Σγ/γ_mean. It says nothing about Beer-Lambert.
  K-Dosimetry's K_system = N·K_qubit row is correct and is no
  evidence for independent absorption.

- **No even 1/N sharing.** At strong coupling the shares settle at 0.29
  and 0.21 on the chain (N = 4) and at 0.34 and 0.16 on the ring, where
  all sites are equivalent; at J = 1 the edge keeps 0.41 to 0.47 for
  every N from 2 to 5.

---

## Reproduction

- Script: [`simulations/k_dosimetry_beer_lambert.py`](../simulations/k_dosimetry_beer_lambert.py)
- Output: [`simulations/results/k_dosimetry_beer_lambert.txt`](../simulations/results/k_dosimetry_beer_lambert.txt)
