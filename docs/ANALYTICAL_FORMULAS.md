# Analytical Formulas Reference

**Purpose:** Before building a Liouvillian, check here. Every formula
below replaces a matrix computation.

**Convention:** J = coupling strength. γ = dephasing rate per qubit.
N = number of qubits. w = XY-weight (count of X/Y in a Pauli string).
Formulas in ASCII; prose uses Unicode (Ψ, Π, Σ, γ).

---

## Spectral Structure (replace eigenvalue computation)

### 1. Palindrome equation (Tier 1, proven)

    Π · L · Π⁻¹ = -L - 2Σγ · I

Every Liouvillian eigenvalue λ has a partner at -λ - 2Σγ.
Every decay rate d pairs with 2Σγ - d.

**Valid for:** Heisenberg, XY, Ising, XXZ, DM; Z-dephasing; any graph;
any N; non-uniform γ per qubit. Two Π families (P1, P4).
**Breaks for:** depolarizing noise (error = (2/3)Σγ, linear in γ and N).
**Replaces:** palindrome verification (54,118 eigenvalues at N=8).
**Source:** [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md)

### 2. w=1 dispersion relation (Tier 1-2, verified N=2-6)

    omega_k = 4J * (1 - cos(pi*k/N)),    k = 1, ..., N-1

N-1 distinct frequencies for the Heisenberg chain. Machine-precision
match for 15 frequencies (N=2-6). Tight-binding model with hopping 2J.

**Valid for:** Heisenberg chain, open boundaries, all N (verified N=2-6).
**Replaces:** full Liouvillian diagonalization for w=1 frequencies.
O(N) instead of O(4^{3N}).
**Source:** [Analytical Spectrum](../experiments/ANALYTICAL_SPECTRUM.md)

### 3. Decay rate bounds (Tier 1, from palindrome proof)

    min rate = 2*gamma       (w=1 modes)
    max rate = 2*(N-1)*gamma (w=N-1 modes)
    bandwidth = 2*(N-2)*gamma

**Valid for:** Heisenberg chain, uniform Z-dephasing, all N.
**Replaces:** eigenvalue range computation.
**Source:** [README](../README.md)

### 4. Stationary mode count (Tier 1, Clebsch-Gordan)

    Stat(N) = Sum_J m(J,N) * (2J+1)^2

m(J,N) = multiplicity of total spin J in N spin-1/2 particles.
Exact for chain topology, lower bound for higher-symmetry topologies.

**Valid for:** Heisenberg Hamiltonian, Σγ = 0, all N.
**Replaces:** null-space computation of Liouvillian.
**Source:** [Cavity Modes Formula](../experiments/CAVITY_MODES_FORMULA.md)

### 5. Depolarizing error (Tier 1, proven)

    error = gamma * 2*(N-2)/3

Linear in γ and N. Hamiltonian-independent.

**Valid for:** any Hamiltonian under depolarizing noise.
**Replaces:** numerical palindrome check for depolarizing channels.
At γ ~ 0.001 (typical IBM): error < 0.1%.
**Source:** [Depolarizing Palindrome](../experiments/DEPOLARIZING_PALINDROME.md)

### 22. GHZ XOR-drain (Tier 2, verified N=2-5)

    GHZ  -> 100% weight in XOR modes (N+1 modes at rate 2*Sigma_gamma)
    W    -> 0% XOR, 100% palindromic  (N >= 3)
    Bell -> 0% XOR, 100% palindromic  (N >= 3, Hamming distance 2)

Predictor: mixed XY Pauli weight, r = 0.976 correlation with XOR
fraction (N >= 3). GHZ fragility is not vague "delicateness" but
exact projection onto the fastest-decaying spectral sector.

**Valid for:** Heisenberg, XY, Ising, XXZ, DM; Z-dephasing; N=2-5.
**Replaces:** GHZ fragility analysis; explains why W-encoding
outperforms GHZ for state transfer (0% vs 100% drain weight).
**Source:** [XOR Space](../experiments/XOR_SPACE.md)

### 23. XOR drain vanishing fraction (Tier 1, combinatorial proof)

    fraction(XOR) = (N+1) / 4^N

N=3: 6.25%. N=5: 0.59%. N=8: 0.014%. N=20: ~10^-11.
GHZ fragility is a small-N phenomenon; at macroscopic N the XOR
sector has measure zero.

**Valid for:** any N, Z-dephasing.
**Replaces:** large-N XOR mode counting; confirms the drain is
irrelevant at macroscopic scale.
**Source:** [N->infinity Palindrome](../experiments/N_INFINITY_PALINDROME.md)

### 33. N=3 exact intermediate decay rates (Tier 1, exact rational)

    rate_1 = 2*gamma       (w=1, pure single-site coherence)
    rate_2 = 8*gamma/3     (w=2 mixed, Hamiltonian superposition)
    rate_3 = 10*gamma/3    (w=2 mixed, Hamiltonian superposition)

Three distinct rates for the N=3 Heisenberg chain. The Hamiltonian
mixes w=1 and w=2 Pauli strings into supermodes with exact rational
decay rates. At N=3 these three rates plus the extremes (0 and
6*gamma) fully determine the spectrum. At N >= 4, internal rates
become topology-dependent (only the boundary rates 2*gamma and
2*(N-1)*gamma remain universal).

**Valid for:** N=3 Heisenberg chain, Z-dephasing.
**Replaces:** Liouvillian diagonalization for N=3 decay rates.
Two independent information channels (frequency vs decay) are
perfectly orthogonal at N=3.
**Source:** [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md)

---

## Q-Factor and V-Effect (replace resonator analysis)

### 6. V-Effect gain (Tier 1-2, verified N=2-6)

    V(N) = 1 + cos(pi/N) = 2*cos^2(pi/(2N))

Q-factor amplification from coupling. γ-independent (cancels in ratio).
For N=5: (5+sqrt(5))/4 = 1.80902. For N→∞: V = 2 (saturation).
Under non-uniform γ: applies only to the extremal (best-Q) mode.

**Valid for:** Heisenberg chain, Z-dephasing, all N.
**Replaces:** paired Liouvillian diagonalization for V-Effect measurement.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

### 7. Q-factor spectrum (Tier 1-2, follows from formula 2)

    Q_k = 2*J/gamma * (1 - cos(pi*k/N))
    Q_max = 2*J/gamma * (1 + cos(pi/N))
    Q_min = 2*J/gamma * (1 - cos(pi/N))
    Q_mean = 2*J/gamma  (exactly, from sum of cos = 0)
    Q_spread = Q_max/Q_min = cot^2(pi/(2N))  (~N^2/pi^2 for large N)

**Valid for:** Heisenberg chain, uniform Z-dephasing, w=1 sector.
**Replaces:** Q-factor computation from eigenvalues.
**Source:** [Analytical Spectrum](../experiments/ANALYTICAL_SPECTRUM.md)

### 8. 2x universal decay law (Tier 2, verified N=2-5)

    rate(unpaired) = 2*N*gamma
    rate(paired mean) = N*gamma
    ratio = 2.00 exactly

**Valid for:** Heisenberg chain, Z-dephasing, N=2-5.
**Replaces:** unpaired mode rate computation.
**Source:** [Energy Partition](../hypotheses/ENERGY_PARTITION.md)

---

## Sacrifice Zone (replace numerical optimization)

### 9. Sacrifice-zone formula (Tier 2, verified N=2-15)

    gamma_edge = N * gamma_base - (N-1) * epsilon
    gamma_other = epsilon

One-line formula. 360x improvement at N=5, 63.5x at N=15 vs V-shape.

**Valid for:** Heisenberg chain, Z-dephasing, SumMI objective.
**Replaces:** numerical optimization (Nelder-Mead, DE).
**Source:** [Resonant Return](../experiments/RESONANT_RETURN.md)

### 10. SumMI quadratic scaling (Tier 2, verified N=2-15)

    SumMI ~ 0.0053 * N^2 + 0.028 * N - 0.062

Under the sacrifice-zone formula. Quadratic, not exponential.

**Valid for:** Heisenberg chain, sacrifice-zone profile, |+>^N initial state.
**Replaces:** time evolution simulation for SumMI estimation.
**Source:** [Signal Analysis Scaling](../experiments/SIGNAL_ANALYSIS_SCALING.md)

### 11. Mode localization profile, N=5 (Tier 2, topological)

    slowest modes: [0.52, 0.63, 0.70, 0.63, 0.52]
    fastest modes: [0.98, 0.87, 0.80, 0.87, 0.98]

Profile is identical under all noise profiles (uniform, sacrifice, IBM).
Correlation edge-weight vs rate: r = 0.994.

**Valid for:** N=5 Heisenberg chain, ALL Z-dephasing profiles.
**Replaces:** eigenvector decomposition of Liouvillian.
**Source:** [Cavity Mode Localization](../experiments/CAVITY_MODE_LOCALIZATION.md)

---

## CΨ Crossing (replace trajectory computation)

### 12. Single-qubit universal crossing fraction (Tier 2)

    t*/T2 = 0.858367
    from: x^3 + x = 1/2,  x = e^{-t/T2}

Platform-independent. Bell states: ~10x entanglement penalty.

**Valid for:** single qubit, maximal superposition, pure dephasing.
**Replaces:** CΨ(t) trajectory simulation for crossing time.
**Source:** [Universal Quantum Lifetime](../experiments/UNIVERSAL_QUANTUM_LIFETIME.md)

### 13. r* threshold (Tier 2-3, 24,073 records)

    r* = T2 / (2*T1) = 0.2128

Separates crossers from non-crossers. Precision 0.000014. Zero false
positives across 133 qubits, 181 days.

**Valid for:** single qubit, amplitude damping + dephasing, T2echo basis.
**Replaces:** CΨ(t) simulation for crossing prediction per qubit.
**Source:** [IBM Hardware Synthesis](../experiments/IBM_HARDWARE_SYNTHESIS.md)

### 14. K-invariance (Tier 2, Lindblad scaling)

    K = gamma * t_cross = constant per bridge type

K_concurrence = 0.03596. K_MI = 0.033. K_correlation = 0.072.
Standard Lindblad time-rescaling (τ = γt), not deep physics.

**Valid for:** any Lindblad system, any bridge metric.
**Replaces:** multi-γ parameter sweeps for crossing time.
**Source:** [Crossing Taxonomy](../experiments/CROSSING_TAXONOMY.md)

### 15. θ compass (Tier 2)

    theta = arctan(sqrt(4*C*Psi - 1))

Angular distance from CΨ = 1/4 boundary. θ = 0 at crossing.

**Valid for:** any system where CΨ is defined, CΨ > 1/4.
**Replaces:** nothing directly, but provides geometric intuition.
**Source:** [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md)

### 24. Generalized crossing equation (Tier 1, algebraic + hardware-validated)

    C   = 1 - b^r + b^{2r}/2 + b^2/2
    Psi = b
    Crossing:  [1 - b^r + b^{2r}/2 + b^2/2] * b = 1/4

    b = e^{-t/T2*},   r = T2*/T1

Extends formula 12 (pure dephasing, r -> 0) to finite T1.
Polynomial approximation (max error < 0.001):

    t*(r)/T2* ~ 0.858 + 0.012*r + 0.375*r^2 - 0.019*r^3 - 0.084*r^4

| r = T2*/T1 | t*/T2* | Regime           |
|------------|--------|------------------|
| r -> 0     | 0.858  | Pure dephasing   |
| r = 0.5    | 0.950  | T1 = 2*T2        |
| r = 1      | 1.141  | T1 = T2          |

Hardware validated: IBM Torino qubit 52 (MAE = 0.053 with fitted T2*).
The 1/4 crossing emerges from a global fit (1/4 was not a fit target).

**Valid for:** single qubit |+> under combined T1 + T2* decay.
**Replaces:** formula 12 when T1 is finite; numerical CPsi(t)
simulation for superconducting qubits.
**Source:** [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md)

### 25. CPsi closed form, Bell+ Z-dephasing (Tier 1, proven)

    CPsi(t) = f * (1 + f^2) / 6,       f = e^{-4*gamma*t}

    dCPsi/dt = -2*gamma*f*(1 + 3*f^2) / 3

Crossing at f* = 0.8612 (from f*(1 + f*^2) = 3/2).
K = gamma * t_cross = 0.0374.

**Valid for:** Bell+ initial state, Z-dephasing, 2 qubits.
**Replaces:** numerical integration for CPsi(t) trajectory.
O(1) evaluation instead of ODE solver.
**Source:** [CPsi Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### 26. CPsi closed form, general Pauli channels (Tier 1, proven)

    CPsi = u * (1 + u^2 + v^2 + w^2) / 12

    u = e^{-alpha*t},   v = e^{-beta*t},   w = e^{-delta*t}
    alpha = 4*(gamma_y + gamma_z)
    beta  = 4*(gamma_x + gamma_z)
    delta = 4*(gamma_x + gamma_y)

Derivative proven < 0 for any nonzero noise (all coefficients in
the bracket are non-negative; product with -u is strictly negative).

**Valid for:** Bell+ initial state, any combination of
(gamma_x, gamma_y, gamma_z), 2 qubits.
**Replaces:** Lindblad master equation solver for multi-axis
noise on Bell states. O(1) instead of matrix exponentiation.
**Source:** [CPsi Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### 27. K values per noise channel (Tier 1, from formula 26)

    K_Z     = 0.0374    (pure Z-dephasing)
    K_X     = 0.0867    (pure X-noise)
    K_Y     = 0.0867    (pure Y-noise)
    K_depol = 0.0440    (depolarizing, gamma/3 each axis)

Complements formula 14 (K per bridge metric). These are K per
noise TYPE, all measured with CPsi on Bell+ state.

**Valid for:** Bell+ state, single-axis or depolarizing noise.
**Replaces:** per-channel crossing time derivation.
**Source:** [CPsi Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### 28. Fixed-point absorber theorem (Tier 1-2)

    CPsi(rho*) < 1/4    for all primitive CPTP maps

Proven analytically:
- Case A: unital maps (rho* = I/d, CPsi = 0)
- Case B: local channels (rho* = product state, CPsi < 1/4)
Verified numerically:
- Case C: 100 random primitive maps, max CPsi(rho*) = 0.138

Consequence: CPsi = 1/4 is an eventual absorber. Every initial
state with CPsi > 1/4 must eventually cross below 1/4.

**Valid for:** any primitive (unique fixed-point) quantum channel,
2 qubits.
**Replaces:** per-channel verification that the fixed point sits
below 1/4.
**Source:** [Subsystem Crossing Proof](proofs/PROOF_SUBSYSTEM_CROSSING.md)

---

## Fold and Recursion (foundational)

### 16. Fold normal form (Tier 1, proven)

    R = C * (Psi + R)^2

Equivalent to Mandelbrot: u -> u^2 + c with c = C*Psi.
Boundary at CΨ = 1/4 (discriminant of fixed-point equation).

**Source:** [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md)

### 17. CΨ monotonicity (Tier 1, proven)

    dCPsi/dt < 0  for all local Markovian channels

Envelope theorem for arbitrary states. 300 random CPTP maps, 0
exceptions. CΨ is Pauli-invariant (DD cannot change it).

**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### 18. Fold threshold (Tier 2, N-independent)

    Sigma_gamma_crit / J = 0.00249 (Bell state)
    Sigma_gamma_crit / J = 0.00497 (product state)

Below: no fold, CΨ oscillates forever. Above: CΨ crosses 1/4
irreversibly. Max/min ratio across N=2-5: 1.015 (1.5% variation).

**Valid for:** Heisenberg chain, Z-dephasing, N=2-5.
**Replaces:** γ sweep to find fold onset.
**Source:** [Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)

### 19. Fragile bridge asymptotic constant (Tier 2)

    gamma_crit * J_bridge -> 0.50  (strong bridge limit)

Instability is Hopf bifurcation, not PT breaking.
Linear regime: gamma_crit = 0.19 * J_bridge.
Optimal: J_bridge ~ 2J, gamma_crit = 0.41.

**Valid for:** coupled gain-loss Heisenberg chains.
**Replaces:** stability analysis for large J_bridge.
**Source:** [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md)

---

## Thermal (replace n_bar sweep)

### 20. Thermal V-Effect: gain decreases, diversity increases (Tier 2)

    V(N, n_bar=0) = 1 + cos(pi/N)  (exact)
    V(N, n_bar=0.5) ~ 1.44  (N=5)
    V(N, n_bar->inf) ~ 1.29  (saturates)

    Frequencies(N=5, n_bar=0) = 111
    Frequencies(N=5, n_bar=5) = 445  (4x diversity gain)

**Valid for:** Heisenberg chain, Z-dephasing + amplitude damping.
**Replaces:** thermal Liouvillian sweep.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

### 21. Self-heating divergence (Tier 2)

    Fixed point: n_bar -> infinity (all 6 configs tested)

Without external cooling, the system thermalizes to maximum entropy.
No self-consistent operating point exists.

**Valid for:** closed Lindblad systems with amplitude damping.
**Replaces:** convergence search for thermal equilibrium.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

---

## Topology and Protocols (replace parameter sweeps)

### 29. Star-topology coupling threshold (Tier 2, N=3)

    J_SB / J_SA >= 1.466    (at gamma = 0.05)

Threshold for observer-observer (AB) crossing through shared
mediator S in star topology. Below: no AB crossing. Above: AB
crosses 1/4. Receiver noise is fatal (gamma_A > 0.2 kills the
connection); sender noise is tolerable (gamma_B <= 0.5).

**Valid for:** 3-qubit star, Heisenberg, Z-dephasing,
Bell_SA x |+>_B initial state.
**Replaces:** coupling sweep for star-topology crossing threshold.
**Source:** [Star Topology Observers](../experiments/STAR_TOPOLOGY_OBSERVERS.md)

### 30. Gamma-as-signal channel capacity (Tier 2, SVD + Shannon)

    Capacity = 15.5 bits    (at 1% measurement noise, sigma = 0.01)
    Independent channels: 5  (full rank, condition number 14.8)

Spatial dephasing profile is a readable information channel.
100% classification accuracy with 4-symbol alphabet (2 bits
empirical). GHZ is completely blind (d_min = 0); |+>^N is the
optimal receiver (phased array, not omnidirectional).

**Valid for:** N=5 Heisenberg chain, Z-dephasing, |+>^5 initial.
**Replaces:** assumption that dephasing is unstructured noise;
channel capacity analysis from scratch.
**Source:** [Gamma as Signal](../experiments/GAMMA_AS_SIGNAL.md)

### 31. Relay protocol MI bound (Tier 2, N=11)

    MI improvement = +83%    (relay + 2:1 coupling vs passive)

Six relay stages, each t_stage = K/gamma. Receiving qubits get
10x noise reduction during their reception phase. Combines three
results: K/gamma timing (formula 14), quiet receiver (formula 29),
and 2:1 impedance matching.

**Valid for:** N=11 Heisenberg chain, Z-dephasing, Bell pair initial.
**Replaces:** passive propagation baseline for long chains.
**Source:** [Relay Protocol](../experiments/RELAY_PROTOCOL.md)

### 32. Optimal protection state (Tier 2, N=3)

    Slow-mode weight: 90%
    Concurrence:      0.364
    XOR weight:       0.02%

Outperforms GHZ (0% slow-mode, 100% XOR drain), W (0% slow-mode),
Bell (7% slow-mode). Composed mainly of |010>, |000>, |100>, |001>.
Loads boundary-tier palindromic pairs (rates 0.10/0.20) that decay
slowest among the dynamic modes.

**Valid for:** N=3 Heisenberg chain, Z-dephasing.
**Replaces:** trial-and-error state selection for dephasing survival.
**Source:** [Error Correction Palindrome](../experiments/ERROR_CORRECTION_PALINDROME.md)

---

## Structural (replace dimensionality arguments)

### 34. Qubit necessity equation (Tier 1, proven)

    d^2 - 2*d = 0    -->    d = 0 (nothing)  or  d = 2 (qubit)

Palindromic dephasing requires exactly 2 immune Pauli choices (I, Z)
and 2 decaying choices (X, Y) per site. This fixes d = 2 algebraically.
236 qutrit dissipators tested: 0/236 palindromic.

**Valid for:** any local dephasing model.
**Replaces:** dimensional search for palindrome-compatible systems.
**Source:** [Qubit Necessity](QUBIT_NECESSITY.md)

### 35. Dual-perspective lifetime ratio (Tier 2, hardware-validated)

    t_cross(Pi side) / t_cross(direct) ~ T1 / T2

IBM Torino qubit 52: CPsi_A crosses 1/4 at ~140 us, CPsi_B (Pi
perspective) at ~895 us. Factor ~6x. The palindromic partner
decays at the T1 rate, not T2.

**Valid for:** single qubit under T1 + T2 decay, both CPsi
perspectives computed from the same density matrix.
**Replaces:** dual-perspective CPsi simulation.
**Source:** [Both Sides Visible](BOTH_SIDES_VISIBLE.md)

---

## Neural Analog (replace neural symmetry analysis)

### 36. Neural palindrome condition (Tier 1-2, proven + verified)

    Q * J * Q + J + 2*S = 0

    Q = E-I neuron swap operator
    J = Jacobian of Wilson-Cowan dynamics
    S = (1/tau_E + 1/tau_I) / 2 * I

Exact structural analog of quantum palindrome (Pi * L * Pi^-1 =
-L - 2*Sigma_gamma * I). Derived algebraically from quantum proof
via E-I swap mapping. C. elegans connectome: residual 0.013 vs
random 0.108 (8x more palindromic than chance).

**Valid for:** Wilson-Cowan neural networks with Dale's Law.
**Replaces:** ad-hoc neural symmetry analysis; connectome
palindromic quality assessment.
**Source:** [Algebraic Palindrome Neural](neural/ALGEBRAIC_PALINDROME_NEURAL.md)

### 37. Neural eigenvalue pairing (Tier 1, from formula 36)

    mu_k + mu_k' = -(1/tau_E + 1/tau_I)

Analog of lambda + lambda' = -2*Sigma_gamma. Every neural mode
pairs with a partner; decay rates sum to the E-I time constant sum.
Verified: mean sum = -0.3012, predicted = -0.300 (1.6% max deviation).

**Valid for:** linearized Wilson-Cowan networks satisfying formula 36.
**Replaces:** neural eigenvalue computation for pairing verification.
**Source:** [Proof Palindrome Neural](neural/proofs/PROOF_PALINDROME_NEURAL.md)

---

## Derived Relations (follow from combining formulas above)

All derivations verified numerically (N=2-5) against Liouvillian
eigenvalues. Full proofs: [docs/proofs/derivations/](proofs/derivations/).

### D1. Bandwidth and mode density (from formula 2) [VERIFIED]

    BW = omega_{N-1} - omega_1 = 8*J * cos(pi/N) --> 8*J

    rho(omega) = N / (pi * sqrt(omega * (8*J - omega)))

Van Hove singularities at band edges. Exact 1D tight-binding
density of states. Max frequency error < 5e-9.

**Valid for:** Heisenberg chain, w=1 sector, large N.
**Replaces:** numerical mode density estimation.

### D2. V-Effect = Q_max / Q_mean (from formulas 6 + 7) [VERIFIED]

    V(N) = Q_max / Q_mean = (1 + cos(pi/N)) / 1

Q_mean = 2*J/gamma exactly. Proof: Sum_{k=1}^{N-1} cos(pi*k/N) = 0
(geometric series identity). Deviation < 3e-15.

### D3. Crossing time ratios (from formula 27) [VERIFIED]

    t_X / t_Z = K_X / K_Z = ln(2) / (2*ln(1/f*)) = 2.320
    t_depol / t_Z = K_depol / K_Z = 1.177

Verified by Lindblad propagation (Bell+ N=2). Ratio deviation < 3e-6.

### D4. Dimensional factor in crossing (from formulas 12 + 25) [VERIFIED]

    Single qubit (d=2):    f*(1 + f*^2) = 1/2
    Bell+ 2-qubit (d=4):   f*(1 + f*^2) = 3/2 = (d-1) * 1/2

The crossing condition scales with Hilbert space dimension as (d-1)/2.
Exact to machine precision.

### D5. Dynamic palindromic mode count (from formulas 4 + 22 + 23) [VERIFIED]

    Oscillating modes = 4^N - (N+1) - Stat(N)

N=3: 36.  N=5: 898.  Fraction --> 1 exponentially.
**Caveat:** exact at gamma -> 0 only. At finite gamma, the
Hamiltonian mixes weight-parity sectors (w with w +/- 2).

### D6. Spectral gap and mixing time (from formulas 1 + 3) [VERIFIED]

    Spectral gap = 2*gamma    (minimum non-zero decay rate)
    Mixing time  <= N*ln(4) / (2*gamma)

Deviation < 1e-14 for all N tested.

### D7. Q-factor distribution (from formulas 2 + 7) [VERIFIED]

    rho(Q) = 1 / (pi * sqrt((Q - Q_min) * (Q_max - Q)))

    Q_min = 2*J/gamma * (1 - cos(pi/N))
    Q_max = 2*J/gamma * (1 + cos(pi/N))
    Variance = (Q_max - Q_min)^2 / 8

Arcsine distribution (U-shaped). Modes cluster at band edges,
not at the mean. Converges with N (variance rel. error < 2% at N=50).

**Valid for:** Heisenberg chain, w=1 sector, large N.
**Replaces:** numerical Q-factor histogram.
**Source:** [D07 Derivation](proofs/derivations/D07_Q_DISTRIBUTION.md)

---

*Each formula in this document is a Liouvillian that does not need
to be built.*
