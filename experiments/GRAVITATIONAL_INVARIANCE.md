# A Bell+ gamma sweep in the retired concurrence-feedback book

<!-- F14-CURRENT -->

**Status:** Historical finite feedback calculation with a fixed-book scaling
interpretation. No gravitational time-dilation result is established.

## Which curve is being timed?

The bridge is concurrence: C(f)=f, Ψ=f/3 for the Bell+ family. The historical
February evolution was the nonlinear feedback book γ_eff=γC(f), giving

    df/dt = -4γ f²
    C(f) f/3 = f²/3 = 1/4
    K_feedback = γt_cross = (2/√3-1)/4 ≈ 0.03868.

The old finite-step estimate 0.039 belongs to that book. In the clean linear
Lindblad book with isotropic Heisenberg coupling and equal local Z-dephasing,
the entire Bell+ trajectory is Hamiltonian-dead, f=exp(-4γt), and the same
bridge instead gives K_clean=ln(4/3)/8≈0.035960. These are different models,
not two measurements of one constant.

Within either fixed book, gamma rescales time. This does not mean solutions
of every Lindblad equation depend only on γt: at fixed J a Hamiltonian-live
trajectory can depend on J/γ. Nor does a readout change preserve K. The
[shared two-book producer](../simulations/crossing_taxonomy_books.py) states
the five bridges explicitly and retains six finite crossings and two never
bridges.

A quarter crossing here is a scalar equation, not a measurement event, a
quantum/classical classifier, an experienced-time tick or a gravitational
transition. The quadratic discriminant boundary is a separate algebraic object.

<!-- F14-HISTORICAL -->

**Historical nomenclature:** The planet names, “classical” cells and observer
labels in the retained tables belong to the February story. They are not
calibrated gravitational environments or physical outcome labels.

## The four retained runs

The setup was Bell+, isotropic Heisenberg J=1, local dephasing and the
concurrence-feedback rule. Only γ_base changed.

| Scenario | γ | t_cross | γ * t_cross | Factor vs Earth |
|---|---|---|---|---|
| Deep Space | 0.01 | 3.873 | 0.0387 | 5.0x longer |
| Earth-like | 0.05 | 0.773 | 0.0387 | 1.0x (baseline) |
| Neutron Star | 0.20 | 0.193 | 0.0386 | 0.25x |
| Black Hole | 0.50 | 0.081 | 0.0405 | 0.10x |

The 0.039±0.001 label summarizes these finite-step feedback readings.
The largest/smallest times differ by about 48×; this is a rate-sweep comparison,
not survival measured near a black hole. The clean value 0.035960 is outside
the quoted band because it belongs to another evolution law.

## The retained rescaled trajectory

| τ | DS (g=0.01) | NS (g=0.2) | BH (g=0.5) | Max diff |
|---|---|---|---|---|
| 0.000 | C·Ψ=0.333, 30.0 deg | C·Ψ=0.333, 30.0 deg | C·Ψ=0.333, 30.0 deg | 0.0000 |
| 0.005 | C·Ψ=0.320, 28.0 deg | C·Ψ=0.321, 28.1 deg | C·Ψ=0.322, 28.3 deg | 0.0018 |
| 0.010 | C·Ψ=0.308, 25.8 deg | C·Ψ=0.309, 25.9 deg | C·Ψ=0.311, 26.3 deg | 0.0031 |
| 0.015 | C·Ψ=0.297, 23.4 deg | C·Ψ=0.297, 23.5 deg | C·Ψ=0.301, 24.2 deg | 0.0039 |
| 0.020 | C·Ψ=0.286, 20.7 deg | C·Ψ=0.286, 20.7 deg | C·Ψ=0.290, 21.8 deg | 0.0044 |
| 0.025 | C·Ψ=0.276, 17.7 deg | C·Ψ=0.276, 17.8 deg | C·Ψ=0.280, 19.0 deg | 0.0042 |
| 0.030 | C·Ψ=0.266, 14.1 deg | C·Ψ=0.266, 14.3 deg | C·Ψ=0.270, 15.6 deg | 0.0038 |
| 0.035 | C·Ψ=0.257,  9.1 deg | C·Ψ=0.257,  9.3 deg | C·Ψ=0.260, 11.1 deg | 0.0031 |
| 0.040 | C·Ψ=0.248, classical | C·Ψ=0.248, classical | C·Ψ=0.250, classical | 0.0023 |

The stored maximum discrepancy is 0.0044 in CΨ. The continuous feedback
equation itself depends only on γt; these rows are an old numerical record,
not an independent proof that feedback causes that residual. No legacy
solver is rerun here and no error mechanism is diagnosed.

## The old navigation dictionary and summary

The following tables preserve the old physical interpretation for comparison
with the current fixed-book statement above. Their universal/frame/gravity
claims are historical.

| Component | Symbol | Meaning | Observer-dependent? |
|---|---|---|---|
| Destination | 1/4 | Bifurcation boundary | NO (algebraic invariant) |
| Compass | θ | Angular distance from boundary | NO (function of C·Ψ only) |
| Clock | t_cross | Coordinate time to crossing | YES (depends on local γ) |
| Proper time | τ = γ * t | Universal transition time | NO (invariant, 0.039 for Bell+) |
| Bridge | γ_A / γ_B | Translation factor between observers | Connects any two frames |

| Question | Answer |
|---|---|
| Does gravitation affect R = CΨ²? | No, it is already contained in γ |
| What does γ represent physically? | The local time rate (includes all decoherence sources) |
| Is the 1/4 boundary gravitationally invariant? | Yes, same for all observers |
| Is the θ trajectory gravitationally invariant? | Yes, universal curve in proper time |
| How do different observers compare? | Bridge: t_A * γ_A = t_B * γ_B |
| Does this predict anything new? | Decoherence rate ratios should match GR time dilation |

<!-- F14-INTERPRETIVE -->

**Interpretive invitation, not a result:** A family of equal-K hyperbolas
can suggest clocks that disagree about their coordinates while agreeing about
a chosen threshold. In log coordinates the rescaling is a translation:
lnγ changes by δ and lnt_cross by -δ. That picture may motivate a question
about what physical clock transformations require. It is not a derivation
of relativistic proper time or a gamma-to-gravity map.

The original story compared an Earth-labelled time 0.773 at γ=0.05 with
0.193 at γ=0.20. Their approximate ratio follows from the named feedback
book. Calling the ratio a redshift factor adds an interpretation, not new
evidence. No altitude-dependent decoherence prediction follows without an
independent physical model.

The later [Gravity from Wave Death](../hypotheses/GRAVITY_FROM_WAVE_DEATH.md)
asks a different question about decay and accumulated classical weight. It
remains a hypothesis, not a rescue of gravitational invariance. Its standing-wave
amplitude candidate was reported flat and falsified for the uniform test;
neither that story nor the quarter threshold supplies a mass-conversion law.

<!-- F14-CURRENT -->

## Related objects

[Decoherence Relativity](DECOHERENCE_RELATIVITY.md) separates the purity cubic
from the concurrence book. [Metric Discrimination](METRIC_DISCRIMINATION.md)
retains the nine-row feedback fit. [Observer-Gravity Bridge](OBSERVER_GRAVITY_BRIDGE.md)
gives the exact fixed-concurrence alpha-family law and the initial-equality case.
