# K-Dosimetry: the exposure number of a state the Hamiltonian cannot touch

**Status:** Closed form for every |+⟩^⊗N row (exact, checked against the full Lindblad propagation); the edge-concentrator row is numerical, on the exact propagator.
**Date:** April 4, 2026
**Source:** [k_dosimetry.py](../simulations/k_dosimetry.py),
[stored output](../simulations/results/k_dosimetry.txt).

A photographer trades aperture against shutter time: bright and brief, dim
and long, the same picture. Film obeys this reciprocity (Bunsen-Roscoe) as
long as only the product of brightness and time matters. The question here
is whether K = γt, illumination times time, plays that part for a fading
quantum state: does the same K bring the state to the same place, whatever
the rate?

For the preparation this experiment uses, |+⟩^⊗N on the isotropic Heisenberg
chain with the same Z-dephasing rate γ on every site, the answer is yes,
exactly, at every γ and at every coupling J. The reason is that the
Hamiltonian never touches the state, and the letter cube says what the
dephasing then charges.

## Why the answer is exact

Each one-site state (I + X)/2 is a sum of two Pauli letters. Under
Z-dephasing (jump √γ·Z, the convention of the script) the letter I sits at
height 0 and pays nothing, the letter X sits at height 1 and fades as
g = e^(−2γt). So without the chain the state at time t is

  ρ(t) = ⊗ (I + g·X)/2.

The chain does not change this. Every Heisenberg bond is XX + YY + ZZ =
2·SWAP − I, and a product of N identical one-site states commutes with
every SWAP, so [H, ρ(t)] = 0 along the whole trajectory: the state is
Hamiltonian-dead, as [Observer-Gravity Bridge](OBSERVER_GRAVITY_BRIDGE.md)
derives for |++⟩ at N = 2;
every further bond is one more SWAP. With equal rates the one-site states
stay identical, so it stays dead. The purity is a product of one-site
purities,

  Tr ρ² = ((1 + g²)/2)^N,

and the exposure needed to reach a purity target P is

  **K = γt = −¼ · ln(2·P^(1/N) − 1).**

J does not appear, and γ enters only through γt. That is reciprocity, exact,
for this state.

## The γ sweep at N = 2, target purity 0.26

K = −¼·ln(2√0.26 − 1) = 0.980469. The full Lindblad propagation (dense
expm, crossing found by root search) gives the same K at every rate tried:

| γ | t_cross | K = γt |
|---|---|---|
| 0.001 | 980.47 | 0.980469 |
| 0.01 | 98.047 | 0.980469 |
| 0.05 | 19.609 | 0.980469 |
| 0.1 | 9.8047 | 0.980469 |
| 0.5 | 1.9609 | 0.980469 |
| 1.0 | 0.98047 | 0.980469 |
| 5.0 | 0.19609 | 0.980469 |

The largest relative distance from the closed form over the eleven rates the
script runs (0.001 to 5) is 4·10⁻¹³, at the longest time.

## Bell+ and |01⟩ never get there

Both stop at purity ½, and the cube says why. Bell+ is (II + XX − YY + ZZ)/4:
XX and YY sit at height 2 and fade, II and ZZ sit at height 0 and stay, and
what stays has purity (1 + 1)/4 = ½. |01⟩ is (II + ZI − IZ − ZZ)/4: the
chain turns ZI − IZ into the spin current XY − YX, which sits at height 2
and fades; II − ZZ stays, again purity ½. A target of 0.26 is below what
the populations keep.

## The purity-target sweep (N = 2, γ = 0.05)

| Target purity | K = γt | t_cross |
|---|---|---|
| 0.900 | 0.02707 | 0.541 |
| 0.750 | 0.07798 | 1.560 |
| 0.500 | 0.22034 | 4.407 |
| 0.400 | 0.33209 | 6.642 |
| 0.300 | 0.58730 | 11.746 |
| 0.270 | 0.80958 | 16.192 |
| 0.260 | 0.98047 | 19.609 |
| 0.255 | 1.15253 | 23.051 |
| 0.251 | 1.55390 | 31.078 |

This is the film's response curve for this state. Near the floor ¼ (the
maximally mixed value of two qubits) the last step costs a logarithm: with
P = ¼ + δ the closed form gives K ≈ −¼·ln(2δ), so each tenfold step
toward ¼ adds the same ¼·ln 10 ≈ 0.576 to the exposure. Purity ¼ is a
purity, not the CΨ quarter boundary; the two are different readings.

## Across N (γ = 0.05, target 1/2^N + 0.01)

| N | K_qubit = γt | K_system = Σγ·t |
|---|---|---|
| 2 | 0.98047 | 1.96094 |
| 3 | 0.91255 | 2.73766 |
| 4 | 0.81885 | 3.27539 |
| 5 | 0.71575 | 3.57876 |

K_system = Nγt = N·K_qubit is arithmetic for a uniform rate profile, so the
ratio N says nothing about independent absorption
([Beer-Lambert](BEER_LAMBERT_BREAKDOWN.md) takes that question up). K_qubit
falls with N because the target moves with N, not because the sites
interact: every row is the closed form with that target.

## The edge-concentrator profile (N = 4)

Put almost all of the light on one end, γ = (0.197, 0.001, 0.001, 0.001),
with the same total Σγ = 0.2 as the uniform γ = 0.05, at J = 1. Now the one-site
states differ, the product no longer commutes with the chain, and the
closed form no longer applies. On the exact propagator:

| Profile | Σγ | t_cross | K_total = Σγ·t |
|---|---|---|---|
| Uniform | 0.200 | 16.377 | 3.2754 |
| Edge concentrator | 0.200 | 19.583 | 3.9166 |

The edge concentrator needs 1.196 times the total exposure to reach the same
purity. The chain is what spreads the edge's dephasing to the rest: with the
coupling switched off the three dim sites would take t ≈ 745 to get there.
These rows compare total exposure at one purity threshold; they say nothing
about mode selection or the mutual-information benefit the concentrator was
built for ([Concentrator Geometry](CONCENTRATOR_GEOMETRY.md)).

## What the picture is, and is not

The photographer's image stays a picture: there is no photon absorption
model and no crystallization event behind these thresholds. What it gets
right is narrow and exact: for a state the Hamiltonian cannot touch, the
exposure γt alone decides how far the state has faded. For a state the
Hamiltonian moves, a γ sweep at fixed J also changes J/γ, and the dose can
move with it; [Gamma-Time Distinction](../docs/GAMMA_TIME_DISTINCTION.md)
holds that negative control.

## Reproduction

`python simulations/k_dosimetry.py` checks the closed form against the full
propagation at every row, reads the Bell+ and |01⟩ floors, and runs the
edge-concentrator profile, in about half a minute.
