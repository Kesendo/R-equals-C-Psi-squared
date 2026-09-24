<!-- QUARTER-CURRENT -->

<!-- CROSSING-CURRENT -->

# Quarter-crossing taxonomy: finite scalar equations in two named books

Current reading: choose the scalar response `C(f)` and the evolution book before
solving `C(f)f/3=1/4`.  The resulting finite/never categories are properties of
those equations, not observer, measurement, or phase categories.

<!-- QUARTER-HISTORICAL -->
**Historical record:** the detailed repaired taxonomy follows under its own
crossing-scope marker.

<!-- CROSSING-CURRENT -->

# Quarter-crossing taxonomy: three scalar-response classes in two dynamics books

**Status:** Finite Bell+ readout taxonomy; clean equations and reconstructed
feedback equations, with the retired February tables preserved below.
**Date:** February 18, 2026.
**Depends on:** [Quarter crossings in two books](OBSERVER_DEPENDENT_CROSSING.md),
[Metric Discrimination](METRIC_DISCRIMINATION.md).

Five definitions of C read one Bell+ local-Z-dephasing trajectory in the clean
book. With Ψ = f/3, the readouts are mutual information
[2−h₂((1+f)/2)]/2, Wootters concurrence f, the correlation bridge min(1, ½ + f²), which is 1
through the crossing, the mutual purity √(P_A·P_B) = ½, and the overlap
|Tr(ρ_A·ρ_B)|² = ¼.
They place the adopted scalar boundary differently because they are different
functions of the same coherence factor f.

The clean Lindblad book uses df/dt = −4γf. The retired feedback book uses
df/dt = −4γC(f)f, so the bridge also changes the dynamics. Mutual information
and concurrence make that law nonlinear; correlation keeps it linear while it
sits on its cap, which lasts past the crossing, and the two constant bridges
give linear constant-rate scalar decay. The feedback family is not one linear
Lindblad generator.

In each book the question is **C(f)f/3 = ¼**. Three bridges cross and two
never cross, giving six finite crossings over the two books.
A scalar readout is not a physical measurement operation; the crossing is
not a quantum/classical transition. Historical Type A/B/C names below describe
scalar-response shapes, not physical observer classes.

## The finite crossing values

| Bridge | K (standard Lindblad) | t at γ = 0.05 | K (tool, feedback) |
|--------|----------------------|---------------|--------------------|
| mutual_info | 0.02966 | 0.593 | 0.033 |
| concurrence | 0.03596 = ln(4/3)/8 | 0.719 | 0.039 |
| correlation | 0.07192 | 1.438 | 0.072 (identical: C = 1 through the crossing makes the feedback inert there) |

The table displays rounded values. Independent regression references in
[test_crossing_taxonomy_books.py](../simulations/tests/test_crossing_taxonomy_books.py)
pin the six K values without rounding the reference constants. For clean
concurrence K = ln(4/3)/8; feedback concurrence gives
K = (2/√3−1)/4. Correlation gives K = ln(4/3)/4 in both books.
Mutual information uses its stated entropy function, not Shannon information
from a chosen classical measurement. The two below-threshold bridges give
CΨ(0) = 1/6 and 1/12 respectively and decrease thereafter.

## Why the gamma sweep has a fixed K here

K = γt_cross is fixed during a γ sweep **within a fixed bridge and the named
Hamiltonian-dead Bell+ book**. Isotropic Heisenberg coupling does not move any
density matrix on this trajectory, so J drops out. In dimensionless time
the clean law is df/dK = −4f and the feedback law is df/dK = −4C(f)f.

For a general state at fixed J, varying γ changes Q = J/γ and need not leave
the dimensionless generator fixed. Joint scaling obeys
L(λJ,λγ) = λL(J,γ); all dimensionful couplings must scale with it.
The [gamma-time negative control](../docs/GAMMA_TIME_DISTINCTION.md) uses
|01⟩ at fixed J = 1 over a 20× γ range and reports deviations up to
0.861 in concurrence at matched γt. The
[scaling gate](../simulations/gamma_unit_scaling_gate.py) independently
compares Bell+ and |01⟩. An initial H eigenstate alone is insufficient:
the dissipator must also keep the whole trajectory Hamiltonian-dead.

Purity is a different bridge from concurrence. The framework F25 reading
CΨ = f(1+f²)/6 has f* ≈ 0.8612 and K ≈ 0.03735; its t* ≈ 0.747
at γ = 0.05 is the [Boundary Navigation](BOUNDARY_NAVIGATION.md) book.
The values 0.03596, 0.03735 and 0.0387 therefore name two readouts and a
changed evolution law, not competing estimates of one universal constant.

## February sweep and scalar classes

<!-- CROSSING-HISTORICAL -->

**Historical nomenclature:** Type A, B and C label the three response shapes
recorded by the retired delta_calc tool, whose source is kept outside the repo
and was read there.
The tables preserve the numerical record, including finite-step discrepancies;
the current reconstruction does not certify every printed late-time entry.

### Setup

| Parameter | Value |
|-----------|-------|
| State | Bell+ (maximally entangled, (\|00⟩+\|11⟩)/√2) |
| Hamiltonian | Heisenberg (J = 1, h = 0) |
| Noise | Local Z-dephasing (σ_z per qubit) |
| Time step | dt = 0.01 |
| γ values | 0.01, 0.05, 0.10, 0.20 |
| Bridge metrics | Concurrence, mutual information, correlation, mutual purity, overlap |

### Crossing coefficients in the feedback book

| Bridge | γ = 0.01 | γ = 0.05 | γ = 0.10 | γ = 0.20 | K |
|--------|----------|----------|----------|----------|----------|
| mutual_info | t=3.264, K=0.03265 | t=0.653, K=0.03265 | t=0.326, K=0.03265 | t=0.163, K=0.03265 | **0.03265** |
| concurrence | t=3.868, K=0.03868 | t=0.773, K=0.03868 | t=0.387, K=0.03868 | t=0.193, K=0.03868 | **0.03868** |
| correlation | t=7.192, K=0.07192 | t=1.438, K=0.07192 | t=0.719, K=0.07192 | t=0.360, K=0.07192 | **0.07192** |

The twelve finite rows are three bridges over four γ values. The committed
default producer uses γ = 0.05; the other γ columns follow from the fixed K
by t = K/γ. These are not independent evidence for a universal scaling law.

### Matched dimensionless time

| τ = γt | Bridge | γ = 0.01 | γ = 0.05 |
|--------|--------|----------|----------|
| 0.005 | concurrence C | 0.980385 | 0.980354 |
| 0.005 | concurrence Ψ | 0.326795 | 0.326785 |

The finite-dt readings sit near the concurrence feedback value
f = 1/(1+4τ) = 0.980392 at τ = 0.005; the clean value is
exp(−4τ) = 0.980199. The small difference between columns does not
change which book was used.

### Type A: constant C, decaying Ψ

| t | C(t) | Ψ(t) | CΨ |
|-----|------|-------|---------|
| 0.0 | 1.000 | 0.333 | 0.333 |
| 0.5 | 1.000 | 0.302 | 0.302 |
| 1.0 | 1.000 | 0.273 | 0.273 |
| 1.437 | 1.000 | 0.250 | **0.250** |
| 1.8 | 0.986 | 0.232 | 0.229 |

The bridge function the taxonomy ran (bridge_type = correlation) is an excess
purity, doubled and capped:
C = min(1, 2(P_AB − P_A·P_B)), which is min(1, ½ + f²) on this trajectory.
It holds at 1 until f = 1/√2, t = ln 2/(8γ) ≈ 1.73, and then slides as
½ + f², so the stored 0.986 at t = 1.8 is the ceiling releasing, not a
sudden physical loss of a protected correlation.

### Type B: both factors decrease

| t | C(t) | Ψ(t) | CΨ |
|-----|-------|-------|---------|
| 0.0 | 1.000 | 0.333 | 0.333 |
| 0.5 | 0.909 | 0.303 | 0.275 |
| 0.773 | 0.866 | 0.289 | **0.250** |

Concurrence and mutual information change with f. Their scalar crossings
precede the constant-C crossing in these books. The feedback spread is about
2.2× (0.0327 to 0.0719); the clean spread is about 2.4×
(0.02966 to 0.07192).

### Type C: initially below the adopted threshold

Mutual purity is the geometric mean √(P_A·P_B) of the subsystem purities, 0.5
here, not their product 0.25. Overlap is |Tr(ρ_A·ρ_B)|², ¼ on this trajectory,
not fidelity to the initial state. They give
the stored rounded initial products 0.167 and 0.083. Neither has a downward
quarter crossing; no quantum/classical conclusion follows from that absence.

| Class | Mechanism | C at crossing | Bridges | K (tool) | K (exact) |
|-------|-----------|---------------|---------|------|------|
| **Type A** | C stable, only Ψ decays | 1.000 | correlation | 0.072 | 0.07192 |
| **Type B** | C and Ψ both decay | 0.84-0.87 | concurrence, mutual_info | 0.039, 0.033 | 0.03596, 0.02966 |
| **Type C** | CΨ(0) < 1/4 already | n/a | mutual_purity, overlap | never | never |

<!-- CROSSING-CURRENT -->

## What the classification tells us

The chosen C(f) determines whether and where the clean scalar crosses.
In the feedback book it also enters the rate. The Liouvillian palindrome
does not derive the taxonomy: eigenvalue pairing and a nonlinear scalar
readout are separate objects.

The [noise-channel record](NOISE_ROBUSTNESS.md) holds one σ_z run: for local noise the retired
tool applied σ_z on every site whatever jump operator it was given, so its σ_x and σ_y
columns repeat it. Run for real, σ_x and σ_y keep every bridge's C curve but
hold Ψ at 1/3 and move every crossing, and depolarizing noise turns the
correlation bridge Type B.

The GHZ/W results in [N-Scaling Barrier](N_SCALING_BARRIER.md) concern
normalization and named scalar readouts. A full-state or reduced-pair
quarter equality does not identify where entanglement lives.

For the listed Bell+ Z-dephasing functions the decay and crossing directions
follow directly from the equations above. Other channels and states require
their own dynamics; the [envelope work](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)
addresses that separate boundary.

<!-- CROSSING-INTERPRETIVE -->

**Interpretive invitation, not a result:** The old observer picture asks how a
relation determines the landmarks one sees. A robust grip, a fading grip and
an unseen landmark are useful images for the three shapes. Turning those
images into experienced time or a detector event needs a physical model.
The [time-as-crossing-rate hypothesis](../hypotheses/TIME_AS_CROSSING_RATE.md)
keeps that question open.

<!-- CROSSING-CURRENT -->

## Reproduction

Run `python simulations/crossing_taxonomy_books.py`, followed by
`python -m pytest simulations/tests/test_crossing_taxonomy_books.py -q`.
The tests compare six fixed independent references, both never-crossing
bridges, gamma rescalings within each book and same-door bridge mutations.
The tighter feedback solve is an empirical convergence check, not a
certified error bound.

[Quarter crossings](OBSERVER_DEPENDENT_CROSSING.md) ·
[Metric Discrimination](METRIC_DISCRIMINATION.md) ·
[Bridge Fingerprints](BRIDGE_FINGERPRINTS.md) ·
[γ as Signal](GAMMA_AS_SIGNAL.md)
