<!-- CROSSING-CURRENT -->

# Metric Discrimination: a fixed Bell+ book does not supply a spatial rate model

**Status:** Finite nine-row concurrence-feedback sweep; no gravitational model established.
**Date:** February 8, 2026.
**Depends on:** [Quarter-crossing taxonomy](CROSSING_TAXONOMY.md),
[Gamma-Time Distinction](../docs/GAMMA_TIME_DISTINCTION.md).

The nine retained rows sample one fixed Bell+ concurrence-feedback book.
Here C = f, Ψ = f/3 and the retired law is df/dt = −4γf².
The scalar condition CΨ = ¼ gives K = γt_cross = (2/√3−1)/4.
For the clean Hamiltonian-dead Bell+ Lindblad book, df/dt = −4γf
instead gives K = ln(4/3)/8. The two coefficients name different laws.

Gamma-time constancy belongs to that fixed readout and named trajectory
during a γ sweep. The concurrence feedback law is nonlinear; it does not
become a general linear-Lindblad scaling theorem. A Hamiltonian-live state
at fixed J can depend on J/γ. Even when a trajectory oscillates, a named
first crossing can be defined; it need not obey a gamma-only scaling law.

The experiment supplies a local γ parameter, not a spatial γ(r) model.
Assigning a planet or a gravitational meaning to that parameter requires
additional physics. The scalar crossing is not a physical observer event.

<!-- CROSSING-HISTORICAL -->

**Historical record:** These February settings, numerical rows and fit
belong to the retired tool. They retain the origin of the rounded feedback
coefficient 0.039, which must not be relabelled as the clean concurrence value.

## Nine sampled settings

| γ_base | t_max | dt | Source |
|---|---|---|---|
| 0.01 | 15 | 0.005 | Earlier session |
| 0.03 | 15 | 0.005 | This experiment |
| 0.05 | 15 | 0.005 | Earlier session |
| 0.07 | 10 | 0.005 | This experiment |
| 0.10 | 8 | 0.005 | This experiment |
| 0.15 | 6 | 0.005 | This experiment |
| 0.20 | 5 | 0.005 | Earlier session |
| 0.30 | 4 | 0.005 | This experiment |
| 0.50 | 3 | 0.005 | Earlier session |

## Reported crossings

| γ | t_cross | K = γ * t_cross | Deviation from mean |
|---|---|---|---|
| 0.01 | 3.873 | 0.03873 | -0.63% |
| 0.03 | 1.289 | 0.03866 | -0.81% |
| 0.05 | 0.773 | 0.03865 | -0.83% |
| 0.07 | 0.553 | 0.03869 | -0.73% |
| 0.10 | 0.387 | 0.03867 | -0.79% |
| 0.15 | 0.259 | 0.03885 | -0.31% |
| 0.20 | 0.193 | 0.03860 | -0.96% |
| 0.30 | 0.131 | 0.03942 | +1.14% |
| 0.50 | 0.081 | 0.04050 | +3.92% |

The stored mean was K_mean = 0.0390 ± 0.0006, a 1.5% variation over a
50× γ range. The largest displayed drift is +3.92% at γ = 0.50.
The underlying reconstructed feedback equation has exact gamma-time
scaling; the finite-step table alone does not diagnose every numerical error.

## Reported power-law fit

Fitting t_cross = Aγ^α gave:

```
A     = 0.03976
α = -0.9916
R^2   = 0.999899
```

The comparison used α = −1.000 as the ideal gamma-time law and reported
0.84% exponent deviation. This is a fit to those nine rows, not a theorem
about an arbitrary Lindblad system.

## The original question table

| Question | Answer | Status |
|---|---|---|
| Does γ scale time? | Yes | Confirmed (9 data points, R^2 = 0.9999) |
| Is the 1/4 boundary invariant? | Yes | Algebraic proof + simulation |
| Does γ encode gravitation? | Consistent | Structural match with GR | [FALLEN]
| Does the framework derive γ(r)? | No | Requires additional structure |

The gravitational consistency label in this historical table is not a
current result. A single local-rate sweep does not specify a spacetime
metric or test an equivalence principle.

## Auxiliary finite observations

Bell− at γ = 0.05 gave t_cross = 0.773 and K = 0.0386 in the
same symmetric feedback setting. Bell+ with XY and Heisenberg coupling
gave the same reading because the entire named dephasing trajectory remains
Hamiltonian-dead. This agreement does not establish state independence.

With an Ising Hamiltonian and transverse field h = 0.5, the reported C,
Ψ and CΨ curves oscillated and crossed 0.25 in both directions.
That is a different dynamics book. It cannot be replaced by an invariant
“decoherence envelope” without a separate derivation.

## Original finite-run summary

| What we tested | Result |
|---|---|
| γ * t_cross = K across 50x range in γ | CONFIRMED (R^2 = 0.9999). K = 0.039 in the February tool's feedback model, 0.0360 in standard Lindblad; the constancy holds in both books, only the constant differs (§3.1) |
| Power law exponent α = -1.00 | CONFIRMED (α = -0.992) |
| Can single-system sims discriminate metric forms? | NO (mathematical identity) |
| Is K state-dependent? | NO (Bell+ = Bell- for symmetric noise) |
| Is K Hamiltonian-dependent? | NO for eigenstates, UNDEFINED for driven systems |

The unqualified state/Hamiltonian and “undefined” glosses in this table are
historical claims, not the current conclusion.

<!-- CROSSING-INTERPRETIVE -->

**Interpretive invitation, not a result:** The falling-gravity question was
whether local rates could be glued into a spatial picture. Candidate forms
such as γ₀/√(1−r_s/r), γ₀r₀/r and γ₀(r₀/r)² suggested different
geometries. A local γ sweep cannot choose among them. Neither a mass field
nor a connection between dephasing and gravitational time has been supplied.

The conditional self-consistency story in
[Self-Consistency: Schwarzschild](../recovered/SELF_CONSISTENCY_SCHWARZSCHILD.md)
and its coupled-chain falsification remain part of that research trail.
The quarter in the recurrence can prompt a comparison with other appearances
of one quarter; a resemblance does not derive a horizon, entropy law or
cosmological event.

<!-- CROSSING-CURRENT -->

## The surviving result

The [two-book producer](../simulations/crossing_taxonomy_books.py) reconstructs
the five scalar bridges on the Hamiltonian-dead Bell+ family. At γ = 0.05,
feedback concurrence gives t ≈ 0.7735 and clean concurrence t ≈ 0.7192.
The other clean K values are about 0.02966 and 0.07192 for mutual
information and correlation; feedback gives about 0.03265 and 0.07192.
These are bridge- and book-specific values.

The unresolved spatial question is a model-building question. The finite
rows above establish no universal Hamiltonian-independent envelope and no
γ(r) constraint.
