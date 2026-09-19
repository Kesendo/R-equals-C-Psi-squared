# Minimum quarter-crossing energy in a finite two-qubit scan

<!-- F14-CURRENT -->

CΨ_max is the sampled maximum for t>=0.01; excludes t=0.
The separately printed initial value is not included. Thus Bell+ has
CΨ(0)=0.3333 but the first post-initial maximum is 0.3320.

<!-- F14-CURRENT -->

**Status:** Finite numerical scan with an exact equal-energy counterexample to
an energy-only crossing criterion. The readout is Wootters tangle divided by
three, CΨ = C_conc²/3. It is not the purity-times-l1 book.

**Source:** [minimum_energy.py](../simulations/minimum_energy.py).
**Reproduced output:** [minimum_energy.txt](../simulations/results/minimum_energy.txt).

## What this asks

Can energy alone tell us whether a named scalar readout reaches one quarter?
The family cos(α)|00⟩ + sin(α)|11⟩ gives a useful answer: every member has
⟨H⟩ = J for H = J(XX+YY+ZZ), yet its initial tangle depends on α.
Energy alone does not decide this crossing. That is a statement about this
family and readout, not an energy-independent criterion for every quantum system.

The dynamics is the clean linear Lindblad equation with equal local
Z-dephasing. On this Hamiltonian-dead family,

    CΨ(t) = sin²(2α) exp(-8γt)/3
    K_conc(α) = γ t_cross = ln(4 sin²(2α)/3)/8.

For 30° < α ≤ 45° the downward crossing is at positive time. At α = 30°
equality is initial; below 30° this trajectory has no downward crossing.
The retained floating-point table prints NO at 30°; it does not override
the exact initial equality.

A readout crossing is not a physical measurement event or a quantum/classical
classifier. The scalar equation comes before any story about a clock.

## What the current producer says

The six alpha rows, eight product preparations and nine J/γ rows in the linked
output are the finite catalogue. For |01⟩, the retained sweep is below one
quarter at J/γ = 5 and above it at 10. That brackets the sampled change;
it is not an optimized critical ratio.

For |01⟩ and |10⟩ the current producer prints 0.308/YES at J/γ=20,
while the historical document table prints 0.309. For |+,0⟩, |0,+⟩ and |+,1⟩ the producer reports CΨ_max = 0.077 and NO.
The old document table below instead says 0.295 and YES. Those rows are kept
as a historical record, not current evidence. No numerical solver or trajectory
was changed in repairing the labels.

## Hamiltonian eigenstates are not stationary states of every open-system law

[H,ρ₀] = 0 removes the Hamiltonian contribution initially; it does not imply
D[ρ₀] = 0, nor that subsequent states commute with H. Here |00⟩ and |11⟩
are stationary for the full generator. In contrast, |++⟩ follows the exact
separable trajectory ρ_A(t)⊗ρ_A(t), with

    ρ_A(t) = (I + exp(-2γt) X)/2
    Tr(ρ(t)²) = ((1 + exp(-4γt))/2)².

Its purity changes, its tangle remains zero, and its Hamiltonian contribution
vanishes throughout because H = J(2 SWAP-I) commutes with ρ_A⊗ρ_A.
No crossing in this readout does not mean no dynamics or no time.

## Three sampled regimes, with the clock specified

1. A trajectory starting above one quarter can cross downward. The law
   t_cross = K/γ belongs to a named Hamiltonian-dead, fixed-bridge trajectory
   such as the alpha family above, or to joint (J,γ) scaling at fixed Q.
   Initial position above a threshold alone does not imply gamma-only scaling.
2. A state starting below can rise above and return, as the sampled |01⟩
   trajectories do. Their J/γ sweep is a local negative control against a
   state-independent gamma-only clock.
3. A trajectory whose maximum stays below one quarter does not cross this
   threshold. It can still evolve and carry other observables.

<!-- F14-HISTORICAL -->

**Historical record:** The following four tables preserve the original
document's numerical catalogue. They are not a second independently verified
producer. In particular the product-state YES entries just discussed disagree
with the committed source; the old 45° maximum-at-zero convention also differs
from the producer's positive-time grid.

### Equal-energy family

| α     | ⟨H⟩  | CΨ(0)  | Crosses? |
|-------|------|--------|----------|
| 45°   | J    | 0.3333 | YES      |
| 35°   | J    | 0.2943 | YES      |
| 31°   | J    | 0.2599 | YES (barely) |
| 30°   | J    | 0.2500 | NO       |
| 25°   | J    | 0.1956 | NO       |
| 15°   | J    | 0.0833 | NO       |

### Single-excitation family

| α     | CΨ(0)  | CΨ_max | t(max) | Crosses? |
|-------|--------|--------|--------|----------|
| 45°   | 0.3333 | 0.3333 | 0.00   | YES      |
| 25°   | 0.1956 | 0.2956 | 0.37   | YES      |
| 15°   | 0.0833 | 0.3028 | 0.38   | YES      |
| 5°    | 0.0101 | 0.3078 | 0.39   | YES      |

### Product preparations

| State   | CΨ(0) | CΨ_max | Crosses? |
|---------|--------|--------|----------|
| \|0,1⟩  | 0.000  | 0.309  | YES      |
| \|1,0⟩  | 0.000  | 0.309  | YES      |
| \|+,0⟩  | 0.000  | 0.295  | YES      |
| \|0,+⟩  | 0.000  | 0.295  | YES      |
| \|+,1⟩  | 0.000  | 0.295  | YES      |
| \|+,+⟩  | 0.000  | 0.000  | NO       |
| \|0,0⟩  | 0.000  | 0.000  | NO       |
| \|1,1⟩  | 0.000  | 0.000  | NO       |

### J/γ sweep for |01⟩

| J/γ   | CΨ_max | Crosses? |
|-------|--------|----------|
| 0.1   | 0.003  | NO       |
| 0.5   | 0.045  | NO       |
| 1.0   | 0.100  | NO       |
| 2.0   | 0.169  | NO       |
| 5.0   | 0.248  | NO       |
| 10.0  | 0.286  | YES      |
| 20.0  | 0.309  | YES      |
| 50.0  | 0.323  | YES      |
| 100.0 | 0.328  | YES      |

<!-- F14-INTERPRETIVE -->

**Interpretive invitation — not a result:** “Pump and drain” is a way to picture
the competition between Hamiltonian-generated entanglement and dephasing in
these examples. One can also ask what a clock would mean if its ticks were
defined by a chosen threshold. The Wheeler-DeWitt comparison belongs to that
question, not to a deduction that an H-eigenstate has no Lindblad dynamics,
or that a missing crossing removes physical or experienced time.

<!-- F14-CURRENT -->

## Reproduction and next question

Run `python simulations/minimum_energy.py` for stdout only; use an explicit
`--output` path to retain a new result. The numerical grids and state updates
are unchanged. A general crossing criterion would have to specify the initial
state, Hamiltonian, channel, readout, threshold and crossing convention.

Related: [Observer-Gravity Bridge](OBSERVER_GRAVITY_BRIDGE.md),
[Coherence Density](COHERENCE_DENSITY.md),
[Time as Crossing Rate](../hypotheses/TIME_AS_CROSSING_RATE.md).
