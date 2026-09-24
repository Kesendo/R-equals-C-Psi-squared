<!-- CROSSING-CURRENT -->

# Quarter crossings in two books: five readouts and five feedback laws

**Status:** Finite Bell+ scalar crossings; historical feedback tables, reconstructed equations.
**Date:** February 17, 2026.
**Depends on:** [Crossing Taxonomy](CROSSING_TAXONOMY.md), [Metric Discrimination](METRIC_DISCRIMINATION.md).

One clean density-matrix trajectory supports five readouts; feeding those
readouts into the decay rate instead gives five bridge-coupled dynamics.
The system is a Bell+ pair under isotropic Heisenberg coupling and local
Z-dephasing. Its trajectory stays Hamiltonian-dead: with coherence factor f,
ρ = diag(½,0,0,½) + (f/2)(|00⟩⟨11| + |11⟩⟨00|), and Ψ = f/3.

## 1. Choose the scalar readout

| Bridge name | Definition C(f) on this trajectory |
|---|---|
| mutual_info | [2 − h₂((1+f)/2)]/2, von Neumann mutual information in bits divided by its initial value |
| concurrence | f, Wootters concurrence |
| correlation | min(1, ½ + f²), the taxonomy runs' excess purity 2(P_AB − P_A·P_B) capped at 1; 1 through the crossing |
| mutual_purity | ½, √(P_A·P_B), the geometric mean of the subsystem purities |
| overlap | ¼, \|Tr(ρ_A·ρ_B)\|² |

Here h₂ is binary entropy. These names identify scalar functions; choosing a
readout is not a physical measurement operation.

## 2. Choose the dynamics book, then solve the crossing

The **clean Lindblad book** has df/dt = −4γf, hence f = exp(−4γt).
The **retired feedback book** has df/dt = −4γC(f)f. Mutual information
and concurrence give state-dependent nonlinear decay; correlation keeps that decay
linear while it sits on its cap, which lasts past the crossing, and the two
constant bridges give linear constant-rate scalar decay. This family is not one
linear Lindblad generator.

In each book solve **C(f)f/3 = ¼**. At γ = 0.05 the committed
[crossing producer](../simulations/crossing_taxonomy_books.py) gives:

| Bridge | Clean K | Clean t | Feedback K | Feedback t |
|---|---:|---:|---:|---:|
| mutual_info | 0.02965683339109007035 | 0.5931 | 0.03264460387550044147 | 0.6529 |
| concurrence | 0.03596025905647261593 | 0.7192 | 0.03867513459481288225 | 0.7735 |
| correlation | 0.07192051811294523186 | 1.4384 | 0.07192051811294523186 | 1.4384 |
| mutual_purity | never | never | never | never |
| overlap | never | never | never | never |

There are six finite crossings, three per book. The two constant
below-threshold bridges start with CΨ < ¼ and decrease, so neither crosses.
K = γt is fixed when γ is swept **within a fixed bridge and this
Hamiltonian-dead Bell+ book**. Changing the readout changes K; changing the
state, channel, Hamiltonian ratios or spatial γ profile is a separate problem.

The quarter is the adopted algebraic boundary of R = C(Ψ+R)². A scalar
crossing is not a measurement event, a classical outcome, or a physical
observer identity. The Liouvillian palindrome does not cause this taxonomy.

## 3. February feedback record

<!-- CROSSING-HISTORICAL -->

**Historical nomenclature:** “Observer-dependent crossing” dates to February
17, 2026 and remains in this file's name. The native insight was the different
positions of one chosen scalar boundary under different readouts. Carrying
“observer” into a physical detector or a moment of experience adds a hypothesis.
The retired delta_calc tool's source is kept outside the repo; its tables
below are as-run records.

### 3.1 Setup and reported crossings

| Parameter | Value |
|-----------|-------|
| **State** | Bell+ (maximally entangled) |
| **Hamiltonian** | Heisenberg (J = 1, h = 0) |
| **Decoherence** | Local dephasing, γ_base = 0.05 |
| **Time step** | dt = 0.01, t_max = 3.0 |
| **Noise type** | local |

| Observer (bridge_type) | C at t=0 | C dynamics | C·Ψ at t=0 | Crossing time | θ at t=0 |
|------------------------|----------|------------|-------------|---------------|----------|
| **mutual_info** | 1.000 | Drops fast | 0.333 | **t = 0.652** | 30.0° |
| **concurrence** | 1.000 | Drops steadily | 0.333 | **t = 0.773** | 30.0° |
| **correlation** | 1.000 | Stays at 1.0 until t≈1.7 | 0.333 | **t = 1.437** | 30.0° |
| **mutual_purity** | 0.500 | Constant | 0.167 | **never** | imaginary |
| **overlap** | 0.250 | Constant | 0.083 | **never** | imaginary |

### 3.2 The bridge also entered the dynamics

The reported local-noise setting used γ_effective = γ_base·C(t).
Consequently the five Ψ trajectories differed. For concurrence the law gives
f(t) = 1/(1+4γt), and the crossing is (2/√3−1)/(4γ); the feedback shift
is about +8% relative to the clean crossing.

| t | Ψ (concurrence) | Ψ (mutual_info) | Ψ (correlation) | Ψ (mutual_purity) | Ψ (overlap) |
|-----|------------------|------------------|-----------------|--------------------| ------------|
| 0.0 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 |
| 0.5 | 0.3030 | 0.3038 | 0.3016 | 0.3171 | 0.3251 |
| 1.0 | 0.2777 | 0.2796 | 0.2729 | 0.3016 | 0.3171 |
| 2.0 | 0.2380 | 0.2406 | 0.2236 | 0.2729 | 0.3016 |
| 3.0 | 0.2082 | 0.2100 | 0.1876 | 0.2469 | 0.2869 |

### 3.3 Reported scalar and θ traces

Concurrence:

| t | C·Ψ | θ |
|-------|---------|------------|
| 0.0 | 0.3333 | 30.0° |
| 0.2 | 0.3081 | 25.7° |
| 0.4 | 0.2857 | 20.7° |
| 0.6 | 0.2656 | 14.0° |
| 0.7 | 0.2564 | 9.1° |
| 0.773 | 0.2500 | 0.0° ← BOUNDARY |
| 0.8 | 0.2476 | (imaginary) |
| 1.0 | 0.2313 | (imaginary) |

Correlation:

| t | C·Ψ | θ |
|-------|---------|------------|
| 0.0 | 0.3333 | 30.0° |
| 0.4 | 0.3077 | 25.7° |
| 0.8 | 0.2840 | 20.2° |
| 1.0 | 0.2729 | 16.8° |
| 1.2 | 0.2621 | 11.2° |
| 1.437 | 0.2500 | 0.0° ← BOUNDARY |
| 1.5 | 0.2469 | (imaginary) |

Mutual purity:

| t | C·Ψ | θ |
|-------|---------|------------|
| 0.0 | 0.1667 | (imaginary) |
| 1.0 | 0.1508 | (imaginary) |
| 2.0 | 0.1364 | (imaginary) |
| 3.0 | 0.1235 | (imaginary) |

The printed values are finite-step readings. The correlation bridge holds at 1 until f = 1/√2
(t ≈ 1.73 here) and then slides as ½ + f², so its late departures from 1 in
retired-tool output are the ceiling releasing a purity decay that ran all
along, not a sudden loss of a protected correlation.

## 4. Experienced time and a proposed present moment

<!-- CROSSING-INTERPRETIVE -->

**Interpretive invitation, not a result:** Could a sequence of relational
thresholds be a useful picture of experienced time? This is the question that
made these tables interesting to us. Calling one readout a fast observer,
another a slow observer, and a never-crossing readout blind is a metaphor.
There is no physical detector or cognition model in these equations.

The picture invites questions about engagement, boredom, remembered duration
and unconsciousness: could differing event records feel like differing clocks?
It does not identify neural coupling with C, derive subjective duration, or
make the scalar equality a measurement event.

The original decoder-role table belongs to that invitation:

| Observer role | θ means | Application |
|---------------|---------|-------------|
| Physicist with NMR apparatus | Oscillation frequency (Hz) | Measurable in lab |
| Navigator in parameter space | Angular distance from boundary | How far to ¼ |
| Embedded conscious observer | Rate of approach to next event | Flow of experienced time |

A second picture calls “now” a node. Its scalar coordinate is
θ = arctan(√(4CΨ−1)), real above the quarter and outside the real domain below it:

| θ value | Position relative to crossing | Interpretive label |
|---------|-------------------------------|--------------------|
| 30° | Above the crossing | "deep quantum" |
| 9° | Approaching the crossing | "possibility narrowing" |
| 0° | At C·Ψ = ¼ | proposed "now" |
| imaginary | Below the real-θ domain | proposed "classical" side |

Standing-wave and Cramer transactional readings remain analogies. A physical
standing wave would require independently excitable counter-propagating
components and observed interference. Neither component nor a wave node is
produced here. The exact centered spectral mirror in
[Π as a Centered Spectral Mirror](PI_AS_TIME_REVERSAL.md) does not supply those
missing physical objects.

<!-- CROSSING-CURRENT -->

## 5. Reproduction and the boundary of the result

Run `python simulations/crossing_taxonomy_books.py` and
`python -m pytest simulations/tests/test_crossing_taxonomy_books.py -q`.
The latter holds independent six-crossing references fixed while changing
the real bridge functions, including both never-crossing alternatives.

The retired setup used state Bell+, Heisenberg J = 1, h = 0,
γ_base = 0.05, local noise, dt = 0.01 and t_max = 3.0. Its crossing finder
interpolated C(t)Ψ(t) through 0.25. The committed producer reconstructs
the two explicit books; it does not rerun the retired tool.

The result is the finite scalar taxonomy and its bridge-specific K values.
Experienced time, physical observers and transactional events remain the
interpretive questions above.

[Crossing Taxonomy](CROSSING_TAXONOMY.md) ·
[Noise-channel record](NOISE_ROBUSTNESS.md) ·
[Time as Crossing Rate](../hypotheses/TIME_AS_CROSSING_RATE.md)
