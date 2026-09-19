<!-- F14-CURRENT -->

# Decoherence clocks: one Bell+ cubic and two evolution books

**Status:** Exact fixed-family scalar crossing, with retained finite numerical
checks. The relativity and experienced-time readings are interpretations.

## The object before the clock

Choose the purity bridge C(f)=(1+f²)/2 and normalized l1 coherence Ψ=f/3
for the Bell+ density-matrix family. Under clean linear Lindblad evolution
with isotropic Heisenberg H and equal local Z-dephasing, the entire family
commutes with H and f(t)=exp(-4γt). Thus

    C(f) f/3 = 1/4
    f³ + f = 3/2
    K_purity = γt_cross = -ln(f_cross)/4 ≈ 0.03735.

This is a scalar threshold for one readout. It is not a physical measurement
event, a quantum/classical boundary, or a phase-transition theorem. The
discriminant of the separate quadratic fixed-point equation does not change
that distinction.

The same state family with Wootters concurrence C(f)=f gives a different
clean-book value, K_conc=ln(4/3)/8≈0.035960. The retired February feedback law
γ_eff=γC(f) instead gives K_conc=(2/√3-1)/4≈0.03868, historically rounded
to 0.039. For purity with that feedback, K≈0.04013. Choose the bridge first,
then the evolution book, then solve the finite equation.

## What is fixed in a gamma sweep

When the dimensionless generator, initial state, readout, target and crossing
convention remain fixed, writing K=γt removes the overall rate unit. The
Hamiltonian-dead Bell+ trajectory meets that condition during a gamma-only
sweep. A general fixed-J trajectory instead depends on J/γ. Unitary evolution
preserves purity but can redistribute basis-fixed l1 coherence; it need not
leave CΨ or a threshold time unchanged.

The standard joint scaling L(aH,aγ)=aL(H,γ) preserves the dimensionless
generator. It does not make gamma-only scaling universal across states,
Hamiltonians, channels or spatial profiles. The negative control in
[GAMMA_TIME_DISTINCTION](../docs/GAMMA_TIME_DISTINCTION.md) and the independently
rerunnable [gamma-unit gate](../simulations/gamma_unit_scaling_gate.py) keep
that distinction executable.

## Retained numerical checks, with their book attached

The original clean Bell+/purity comparison reports K=0.037350 analytically
against 0.037345 numerically, a 0.014% difference. The rounded crossing
coordinates remain

    f_cross = 0.8612
    Ψ_cross = 0.2871
    C_cross = 0.8709
    C Ψ = 0.2500.

The run reported no printed K change between J=0 and J=10 and spread below
0.1% across γ=0.01 to 1.0. That is consistent with this H-dead family,
not evidence that every Hamiltonian leaves a crossing unchanged.

For two independent unit-prefactor Z dephasers the Bell coherence decays
at 4γ. For the additive collective jump Z⊗I+I⊗Z it decays at 8γ; the same
purity-family cubic then gives K≈0.01868. A product jump Z⊗Z is a different
operator and leaves this Bell coherence untouched. Sharing a scalar equation
does not make noise models equivalent.

The positive root above is from the chosen state-family relation and threshold;
it is not derived from R=CΨ² alone. Ψ=l1/(d-1) is the declared normalization.
A GHZ state at N≥3 has initial Ψ=1/(2^N-1) below the quarter threshold in
this purity book. That does not classify its entanglement or reality.

## Log coordinates within this book

For a positive finite K, lnγ+lnt_cross=lnK is a straight line. Relabelling
the rate changes the coordinates along that line. This does not itself
construct the Lorentz group or a gravitational metric.

The retained Bell+ logarithmic coordinates are

    ξ₀ = ln(1/3) = -1.099
    ξ_cross = ln(0.2871) = -1.248
    Δξ = -0.149.

They describe the same fixed-family threshold as above. Other states and
readouts need not traverse this interval; Hamiltonian-live l1 coherence need
not even be log-linear. A Markovian generator can curve ln(l1), so curvature
alone does not certify memory.

<!-- F14-HISTORICAL -->

**Historical reading:** This comparison table preserves the original
relativity dictionary, including its gravitational identifications. The
gamma-to-gravity, “proper time” and observer entries are not current results.

| Property | Einstein (GR) | Decoherence Relativity |
|----------|--------------|----------------------|
| Invariant | ds² | K = γ · t_cross |
| Coordinates | (ct, x) | (ln γ, ln t_cross) |
| Transformation | Lorentz boost | Translation δ |
| Parameter | φ = arctanh(v/c) | δ = ln(g_B/g_A) |
| Source | velocity / gravity | gravity + initial state |
| γ formula | Φ/c² = g·R/c² | 2·m·g·Δx/ℏ |
| K for Bell+ | N/A | 0.03735 (local dephasing, b³+b=3/2) |
| Crossing | N/A | C = 0.87, Ψ = 0.29 |
| Geometry | hyperbolic rotation | translation on line |

The original gravity construction used γ=2m g Δx/ℏ, the ratio
γ_Earth/γ_Mars=g_Earth/g_Mars=2.6371, and cancellation of g in
γ/(Φ/c²)=2mΔx c²/(ℏR). These are consequences of the proposed substitution;
the fixed-book quantum calculation does not establish that substitution.

<!-- F14-INTERPRETIVE -->

**Interpretive invitation — not a result:** Two readers can draw the same
hyperbola and tell different stories about its coordinates. “Proper
decoherence time” is one such story; a comparison to relativity can help us ask
what an invariant clock would require. Neither γ nor K here is a gravitational
field, physical proper time or experienced duration.

The standing-wave picture asks another question: what becomes visible when
two contributions are added before squaring? Its cross-term is ordinary
algebra. Calling that term a relation between observers does not derive
information conservation, physical interference or a measurement event.

The notebook can keep these questions without asking the Bell+ cubic to
answer them. To connect to gravity would require an independently specified
physical coupling and an observation that distinguishes it from other noise.

<!-- F14-CURRENT -->

## Where to continue

[Crossing Taxonomy](CROSSING_TAXONOMY.md) gives the five bridges in two books.
[Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md) retains the old
feedback scan under its own label. A new crossing comparison should specify
the full state, channel, Hamiltonian, readout and target before comparing K.
