# Analytical Formulas Reference

**Purpose:** Before building a Liouvillian, check here. Every formula
below replaces a matrix computation.

**Convention:** J = coupling strength. γ = dephasing rate per qubit.
N = number of qubits. w = XY-weight (count of X/Y in a Pauli string).
Formulas in ASCII; prose uses Unicode (Ψ, Π, Σ, γ).

---

## Spectral Structure (replace eigenvalue computation)

### F1. Palindrome equation (Tier 1, proven)

    Π · L · Π⁻¹ = -L - 2Σγ · I

Every Liouvillian eigenvalue λ has a partner at -λ - 2Σγ.
Every decay rate d pairs with 2Σγ - d.

**Valid for:** Heisenberg, XY, Ising, XXZ, DM; Z-dephasing; any graph;
any N; non-uniform γ per qubit. Two Π families (P1, P4).
**Breaks for:** depolarizing noise (error = (2/3)Σγ, linear in γ and N).
**Replaces:** palindrome verification (54,118 eigenvalues at N=8).
**Source:** [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md)

### F2. w=1 dispersion relation (Tier 1, proven D10)

    omega_k = 4J * (1 - cos(pi*k/N)),    k = 1, ..., N-1

N-1 distinct frequencies for the Heisenberg chain. Machine-precision
match for 15 frequencies (N=2-6). Tight-binding model with hopping 2J.
Three independent validations: (1) eigenvalue match < 1e-12, (2) Poisson
spacing in w=1 sector (RMT), (3) SFF modulation peak at omega_1 matches
to <1% for N=2-4, 6 ([Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)).

**Valid for:** Heisenberg chain, open boundaries, all N (verified N=2-6).
**Replaces:** full Liouvillian diagonalization for w=1 frequencies.
O(N) instead of O(4^{3N}).
**Source:** [Analytical Spectrum](../experiments/ANALYTICAL_SPECTRUM.md)

### AT. Absorption Theorem (Tier 1, proven)

    Re(lambda_k) = -2*gamma * <n_XY>_k       (exact)

The absorption rate of any Liouvillian eigenmode equals twice the
dephasing rate times the mode's mean light content. Here ⟨n_XY⟩ is the
expectation of the X/Y Pauli factor count in the eigenvector's Pauli
decomposition: ⟨n_XY⟩ = Σ_P |c_P|² n_XY(P) / ||v||².

**Proof sketch:** L_H is anti-Hermitian (v†L_Hv purely imaginary). L_D
is diagonal in the Pauli basis with eigenvalues -2γ n_XY(P). For any
right eigenvector v: Re(λ) = v†L_D v/||v||² = -2γ⟨n_XY⟩. Full proof
in three steps.

**The absorption quantum is 2γ.** Each X/Y Pauli factor costs exactly 2γ
in absorption rate. The spectrum is a ladder with rung spacing 2γ. The
Hamiltonian smooths the ladder (⟨n_XY⟩ can be non-integer) but cannot
change the endpoints or the fundamental quantum.

**F3, F8, F33, and D6 are corollaries** of this theorem.
The palindromic sum rule (α_fast + α_slow = 2Σγ) follows from combining
this theorem with the palindromic weight swap (⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N).

**Valid for:** any real Hermitian Hamiltonian, Z-dephasing, any graph,
any N, non-uniform γ_k per site (replace 2γ with 2Σ_k γ_k × [σ_k ∈ {X,Y}]).
**Breaks for:** complex Hermitian Hamiltonians (DM interactions), where
L_H is not anti-Hermitian.
**Replaces:** eigenvalue range computation; palindromic sum rule verification;
spectral gap derivation; unpaired mode rate identification.
**Verified:** 1,343 modes, N=2-5, γ=0.01-1.0, J=0.1-5.0, CV = 0.0000.
**Source:** [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md),
[Absorption Theorem Discovery](../experiments/ABSORPTION_THEOREM_DISCOVERY.md)

### F3. Decay rate bounds (Tier 1, corollary of Absorption Theorem)

    min rate = 2*gamma       (w=1 modes, pure sector)
    max rate = 2*(N-1)*gamma (w=N-1 modes)
    bandwidth = 2*(N-2)*gamma

**Now a corollary of the Absorption Theorem (AT):** min = 2γ because
the smallest nonzero ⟨n_XY⟩ ≈ 1 (pure weight-1 modes). Max = 2(N-1)γ
for the fastest paired modes (⟨n_XY⟩ ≈ N-1). The XOR drain at 2Nγ
(⟨n_XY⟩ = N) sits above this range.

**Caveat resolved:** At N ≥ 4, Hamiltonian mixing creates hybrid modes
with rates below 2γ (N=4: 0.98γ, N=5: 0.62γ). These are NOT exceptions:
they are mixed-sector modes with fractional ⟨n_XY⟩ < 1. The Absorption
Theorem holds exactly for these modes; the rate 2γ⟨n_XY⟩ is correct
for non-integer ⟨n_XY⟩ ([Proton Water Chain](../experiments/PROTON_WATER_CHAIN.md)).

**Valid for:** Heisenberg chain, uniform Z-dephasing, all N.
**Replaces:** eigenvalue range computation.
**Source:** [README](../README.md),
[Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md)

### F4. Stationary mode count (Tier 1, Clebsch-Gordan decomposition)

    Stat(N) = Sum_J m(J,N) * (2J+1)^2

m(J,N) = multiplicity of total spin J in N spin-1/2 particles.
Exact for chain topology, lower bound for higher-symmetry topologies.

**Valid for:** Heisenberg Hamiltonian, Σγ = 0, all N.
**Replaces:** null-space computation of Liouvillian.
**Source:** [Cavity Modes Formula](../experiments/CAVITY_MODES_FORMULA.md)

### F5. Depolarizing error (Tier 1, proven)

    error = gamma * 2*(N-2)/3

Linear in γ and N. Hamiltonian-independent.

**Valid for:** any Hamiltonian under depolarizing noise.
**Replaces:** numerical palindrome check for depolarizing channels.
At γ ~ 0.001 (typical IBM): error < 0.1%.
**Source:** [Depolarizing Palindrome](../experiments/DEPOLARIZING_PALINDROME.md)

### F22. GHZ XOR-drain (Tier 2, verified N=2-5)

    GHZ  -> 100% weight in XOR modes (N+1 modes at rate 2*Σγ)
    W    -> 0% XOR, 100% palindromic  (N >= 3)
    Bell -> 0% XOR, 100% palindromic  (N >= 3, Hamming distance 2)

Predictor: mixed XY Pauli weight, r = 0.976 correlation with XOR
fraction (N >= 3). GHZ fragility is not vague "delicateness" but
exact projection onto the fastest-decaying spectral sector.

**Valid for:** Heisenberg, XY, Ising, XXZ, DM; Z-dephasing; N=2-5.
**Replaces:** GHZ fragility analysis; explains why W-encoding
outperforms GHZ for state transfer (0% vs 100% drain weight).
**Source:** [XOR Space](../experiments/XOR_SPACE.md)

### F23. XOR drain vanishing fraction (Tier 1, combinatorial proof)

    fraction(XOR) = (N+1) / 4^N

N=3: 6.25%. N=5: 0.59%. N=8: 0.014%. N=20: ~10^-11.
GHZ fragility is a small-N phenomenon; at macroscopic N the XOR
sector has measure zero.

**Valid for:** any N, Z-dephasing.
**Replaces:** large-N XOR mode counting; confirms the drain is
irrelevant at macroscopic scale.
**Source:** [N->infinity Palindrome](../experiments/N_INFINITY_PALINDROME.md)

### F33. N=3 exact intermediate decay rates (Tier 1, exact rational)

    rate_1 = 2*gamma       (w=1, pure single-site coherence)
    rate_2 = 8*gamma/3     (w=2 mixed, Hamiltonian superposition)
    rate_3 = 10*gamma/3    (w=2 mixed, Hamiltonian superposition)

Three distinct rates for the N=3 Heisenberg chain. The Hamiltonian
mixes w=1 and w=2 Pauli strings into supermodes with exact rational
decay rates. At N=3 these three rates plus the extremes (0 and
6*gamma) fully determine the spectrum. At N >= 4, internal rates
become topology-dependent (only the boundary rates 2*gamma and
2*(N-1)*gamma remain universal).

**Absorption Theorem interpretation:** The fractional rates correspond
to fractional ⟨n_XY⟩ from Hamiltonian mixing:

    rate_1 = 2*gamma     -> <n_XY> = 1    (pure weight-1)
    rate_2 = 8*gamma/3   -> <n_XY> = 4/3  (mix of w=1 and w=2)
    rate_3 = 10*gamma/3  -> <n_XY> = 5/3  (mix of w=1 and w=2)

The theorem α = 2γ⟨n_XY⟩ holds exactly, including for non-integer ⟨n_XY⟩.

**Valid for:** N=3 Heisenberg chain, Z-dephasing.
**Replaces:** Liouvillian diagonalization for N=3 decay rates.
Two independent information channels (frequency vs decay) are
perfectly orthogonal at N=3.
**Source:** [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md)

### F50. Weight-1 degeneracy / conserved operator count (Tier 1, proven + verified N=2-7)

    d_real(Re = -2*gamma) = 2N

Exactly 2N purely-real Liouvillian eigenvalues at the first non-zero
grid position. These are the SWAP-invariant conserved operators:

    T_c^{(a)} = Sum_j Sum_{S subset complement(j), |S|=c} sigma_a^{(j)} x Z_S x I_rest

for a in {X, Y} and c = 0, 1, ..., N-1. Each T_c^{(a)} commutes with H
because the Heisenberg Hamiltonian is a sum of SWAPs, and SWAP preserves
both the active Pauli type (X or Y) and the Z-count c. The 2N operators
are linearly independent (disjoint Pauli string support).

Special cases:
- T_0^{(X)} = 2*S_x, T_0^{(Y)} = 2*S_y (global SU(2) generators)
- T_{N-1}^{(a)} = Sum_j sigma_a^{(j)} x Z_{all others} (Jordan-Wigner-type)

Lower bound dim(ker) >= 2N: proven (SWAP invariance constructs 2N kernel vectors).
Upper bound dim(ker) <= 2N: proven (triangle inequality forces each SWAP to
fix v individually; adjacent transpositions generate S_N; one invariant vector
per transitive orbit). Numerically verified for N = 2, ..., 7.

**Valid for:** ANY connected graph with isotropic Heisenberg coupling
+ uniform Z-dephasing. Not only chains -- also star, ring, complete, tree.
**Breaks for:** anisotropic XXZ (Delta != 1), where ZZ term mixes X/Y types.
**Caveat:** This universality is UNIQUE to k=0 and k=1. For k >= 2,
d_real(k) is topology-dependent (Chain < Star < Ring < Complete).
See [Weight-2 Kernel](../experiments/WEIGHT2_KERNEL.md).
**Replaces:** eigenvector analysis at the first grid position;
numerical counting of purely-real eigenvalues.
**Source:** [Weight-1 Degeneracy Proof](proofs/PROOF_WEIGHT1_DEGENERACY.md)

---

## Q-Factor and V-Effect (replace resonator analysis)

### F6. V-Effect gain (Tier 1-2, verified N=2-6)

    V(N) = 1 + cos(pi/N) = 2*cos^2(pi/(2N))

Q-factor amplification from coupling. γ-independent (cancels in ratio).
For N=5: (5+sqrt(5))/4 = 1.80902. For N→∞: V = 2 (saturation).
Under non-uniform γ: applies only to the extremal (best-Q) mode.

**Valid for:** Heisenberg chain, Z-dephasing, all N.
**Physically validated:** proton water chain N=1-5, machine-precision
match ([Proton Water Chain](../experiments/PROTON_WATER_CHAIN.md)).
**Replaces:** paired Liouvillian diagonalization for V-Effect measurement.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

### F7. Q-factor spectrum (Tier 1, corollary of D10)

    Q_k = 2*J/gamma * (1 - cos(pi*k/N))
    Q_max = 2*J/gamma * (1 + cos(pi/N))
    Q_min = 2*J/gamma * (1 - cos(pi/N))
    Q_mean = 2*J/gamma  (exactly, from sum of cos = 0)
    Q_spread = Q_max/Q_min = cot^2(pi/(2N))  (~N^2/pi^2 for large N)

**Valid for:** Heisenberg chain, uniform Z-dephasing, w=1 sector.
**Replaces:** Q-factor computation from eigenvalues.
**Source:** [Analytical Spectrum](../experiments/ANALYTICAL_SPECTRUM.md)

### F8. 2× universal decay law (Tier 1, corollary of Absorption Theorem)

    rate(unpaired) = 2*N*gamma       (<n_XY> = N, pure {X,Y}^N)
    rate(paired mean) = N*gamma      (<n_XY>_fast + <n_XY>_slow = N)
    ratio = 2.00 exactly

**Now a corollary:** The palindromic sum rule α_fast + α_slow = 2Σγ
follows from the Absorption Theorem (α = 2γ⟨n_XY⟩) combined with the
palindromic weight swap (⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N, proven in
[Primordial Superalgebra](../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)).
The "ratio 2" is the ratio of the full range (0 to 2Σγ) to the center
(Σγ) of a symmetric interval; it is a definition, not a separate law.

**Valid for:** any real Hermitian Hamiltonian, Z-dephasing, any graph, any N.
**Replaces:** unpaired mode rate computation; palindromic sum verification.
**Source:** [Energy Partition](../hypotheses/ENERGY_PARTITION.md),
[Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md)

---

## Sacrifice Zone (replace numerical optimization)

### F9. Sacrifice-zone formula (Tier 2, verified N=2-15)

    gamma_edge = N * gamma_base - (N-1) * epsilon
    gamma_other = epsilon

One-line formula. 360x improvement at N=5, 63.5x at N=15 vs V-shape.

**Valid for:** Heisenberg chain, Z-dephasing, SumMI objective.
**Replaces:** numerical optimization (Nelder-Mead, DE).
**Source:** [Resonant Return](../experiments/RESONANT_RETURN.md)

### F10. SumMI quadratic scaling (Tier 2, verified N=2-15)

    SumMI ~ 0.0053 * N^2 + 0.028 * N - 0.062

Under the sacrifice-zone formula. Quadratic, not exponential.

**Valid for:** Heisenberg chain, sacrifice-zone profile, |+>^N initial state.
**Replaces:** time evolution simulation for SumMI estimation.
**Source:** [Signal Analysis Scaling](../experiments/SIGNAL_ANALYSIS_SCALING.md)

### F11. Mode localization profile, N=5 (Tier 2, geometric)

    slowest modes (edge sacrifice): [0.52, 0.63, 0.70, 0.63, 0.52]
    fastest modes (edge sacrifice): [0.98, 0.87, 0.80, 0.87, 0.98]

Profile of each individual mode is identical under all noise profiles
(geometric, from standing wave patterns sin(πkj/N)). But which modes
are slowest depends on the noise profile: edge sacrifice selects
center-heavy modes, center sacrifice selects edge-heavy modes. Not
topologically protected (winding number = 0, Berry phase not quantized).
Correlation edge-weight vs rate: r = 0.994.

**Valid for:** N=5 Heisenberg chain, ALL Z-dephasing profiles.
**Replaces:** eigenvector decomposition of Liouvillian.
**Source:** [Cavity Mode Localization](../experiments/CAVITY_MODE_LOCALIZATION.md),
[Topological Edge Modes](../experiments/TOPOLOGICAL_EDGE_MODES.md)

---

## CΨ Crossing (replace trajectory computation)

### F12. Single-qubit universal crossing fraction (Tier 2)

    t*/T2 = 0.858367
    from: x^3 + x = 1/2,  x = e^{-t/T2}

Platform-independent. Bell states: ~10x entanglement penalty.

**Valid for:** single qubit, maximal superposition, pure dephasing.
**Replaces:** CΨ(t) trajectory simulation for crossing time.
**Source:** [Universal Quantum Lifetime](../experiments/UNIVERSAL_QUANTUM_LIFETIME.md)

### F13. r* threshold (Tier 2-3, 24,073 records)

    r* = T2 / (2*T1) = 0.2128

Separates crossers from non-crossers. Precision 0.000014. Zero false
positives across 133 qubits, 181 days.

**Valid for:** single qubit, amplitude damping + dephasing, T2echo basis.
**Replaces:** CΨ(t) simulation for crossing prediction per qubit.
**Source:** [IBM Hardware Synthesis](../experiments/IBM_HARDWARE_SYNTHESIS.md)

### F14. K-invariance (Tier 2, Lindblad scaling)

    K = gamma * t_cross = constant per bridge type

K_concurrence = 0.03596. K_MI = 0.033. K_correlation = 0.072.
Standard Lindblad time-rescaling (τ = γt), not deep physics.

**Interpretation:** In the [optical cavity analogy](../experiments/OPTICAL_CAVITY_ANALYSIS.md),
γ plays the role of external illumination and t is the system's experienced
duration. K = γ × t is the invariant decoherence dose: more light means
shorter experience, less light means longer, the product unchanged.
Structural parallel to c × τ in relativity (observation, not identification).

**Connection to Absorption Theorem:** γ is the absorption quantum rate
(2γ per X/Y Pauli factor). K = γ × t is the total absorbed dose. The
invariance of K means the total dose is state-dependent but
parameter-independent. See [K-Dosimetry](../experiments/K_DOSIMETRY.md).

**Valid for:** any Lindblad system, any bridge metric.
**Replaces:** multi-γ parameter sweeps for crossing time.
**Source:** [Crossing Taxonomy](../experiments/CROSSING_TAXONOMY.md)

### F15. θ compass (Tier 2)

    theta = arctan(sqrt(4*C*Psi - 1))

Angular distance from CΨ = 1/4 boundary. θ = 0 at crossing.

**Valid for:** any system where CΨ is defined, CΨ > 1/4.
**Replaces:** nothing directly, but provides geometric intuition.
**Source:** [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md)

### F24. Generalized crossing equation (Tier 1, algebraic + hardware-validated)

    C   = 1 - b^r + b^{2r}/2 + b^2/2
    Psi = b
    Crossing:  [1 - b^r + b^{2r}/2 + b^2/2] * b = 1/4

    b = e^{-t/T2*},   r = T2*/T1

Extends F12 (pure dephasing, r -> 0) to finite T1.
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
**Replaces:** F12 when T1 is finite; numerical CΨ(t)
simulation for superconducting qubits.
**Source:** [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md)

### F25. CΨ closed form, Bell+ Z-dephasing (Tier 1, proven)

    CPsi(t) = f * (1 + f^2) / 6,       f = e^{-4*gamma*t}

    dCPsi/dt = -2*gamma*f*(1 + 3*f^2) / 3

Crossing at f* = 0.8612 (from f*(1 + f*^2) = 3/2).
K = gamma * t_cross = 0.0374.

**Valid for:** Bell+ initial state, Z-dephasing, 2 qubits.
**Replaces:** numerical integration for CΨ(t) trajectory.
O(1) evaluation instead of ODE solver.
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F26. CΨ closed form, general Pauli channels (Tier 1, proven)

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
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F27. K values per noise channel (Tier 1, from F26)

    K_Z     = 0.0374    (pure Z-dephasing)
    K_X     = 0.0867    (pure X-noise)
    K_Y     = 0.0867    (pure Y-noise)
    K_depol = 0.0440    (depolarizing, gamma/3 each axis)

Complements F14 (K per bridge metric). These are K per
noise TYPE, all measured with CΨ on Bell+ state.

**Valid for:** Bell+ state, single-axis or depolarizing noise.
**Replaces:** per-channel crossing time derivation.
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F28. Fixed-point absorber theorem (Tier 1-2)

    CPsi(rho*) < 1/4    for all primitive CPTP maps (completely positive trace-preserving: the most general physically allowed quantum operations)

Proven analytically:
- Case A: unital maps (rho* = I/d, CΨ = 0)
- Case B: local channels (rho* = product state, CΨ < 1/4)
Verified numerically:
- Case C: 100 random primitive maps, max CΨ(rho*) = 0.138

Consequence: CΨ = 1/4 is an eventual absorber. Every initial
state with CΨ > 1/4 must eventually cross below 1/4.

**Valid for:** any primitive (unique fixed-point) quantum channel,
2 qubits.
**Replaces:** per-channel verification that the fixed point sits
below 1/4.
**Source:** [Subsystem Crossing Proof](proofs/PROOF_SUBSYSTEM_CROSSING.md)

---

## Fold and Recursion (foundational)

### F16. Fold normal form (Tier 1, proven)

    R = C * (Psi + R)^2

Equivalent to Mandelbrot: u -> u^2 + c with c = C*Psi.
Boundary at CΨ = 1/4 (discriminant of fixed-point equation).

**Source:** [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md)

### F17. CΨ monotonicity (Tier 1, proven)

    dCPsi/dt < 0  for all local Markovian channels

Envelope theorem for arbitrary states. 300 random CPTP maps, 0
exceptions. CΨ is Pauli-invariant (DD cannot change it).

**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F18. Fold threshold (Tier 2, N-independent)

    Σγ_crit / J = 0.00249 (Bell state)
    Σγ_crit / J = 0.00497 (product state)

Below: no fold, CΨ oscillates forever. Above: CΨ crosses 1/4
irreversibly. Max/min ratio across N=2-5: 1.015 (1.5% variation).

**Valid for:** Heisenberg chain, Z-dephasing, N=2-5.
**Replaces:** γ sweep to find fold onset.
**Source:** [Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)

### F19. Fragile bridge asymptotic constant (Tier 2)

    gamma_crit * J_bridge -> 0.50  (strong bridge limit)

Instability is a Hopf bifurcation (eigenvalues leave the imaginary
axis as a conjugate pair, creating oscillatory instability), now
identified as Liouvillian chiral
symmetry breaking (Π forces λ ↔ −λ at Σγ = 0; eigenvalues leave
the imaginary axis at γ_crit).
Linear regime: gamma_crit = 0.19 * J_bridge.
Optimal: J_bridge ~ 2J, gamma_crit = 0.41.

**Valid for:** coupled gain-loss Heisenberg chains.
**Replaces:** stability analysis for large J_bridge.
**Source:** [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md),
[PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

---

## Thermal (replace n_bar sweep)

### F20. Thermal V-Effect: gain decreases, diversity increases (Tier 2)

    V(N, n_bar=0) = 1 + cos(pi/N)  (exact)
    V(N, n_bar=0.5) ~ 1.44  (N=5)
    V(N, n_bar->inf) ~ 1.29  (saturates)

    Frequencies(N=5, n_bar=0) = 111
    Frequencies(N=5, n_bar=5) = 445  (4x diversity gain)

**Valid for:** Heisenberg chain, Z-dephasing + amplitude damping.
**Replaces:** thermal Liouvillian sweep.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

### F21. Self-heating divergence (Tier 2)

    Fixed point: n_bar -> infinity (all 6 configs tested)

Without external cooling, the system thermalizes to maximum entropy.
No self-consistent operating point exists.

**Valid for:** closed Lindblad systems with amplitude damping.
**Replaces:** convergence search for thermal equilibrium.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

---

## Topology and Protocols (replace parameter sweeps)

### F29. Star-topology coupling threshold (Tier 2, N=3)

    J_SB / J_SA >= 1.466    (at gamma = 0.05)

Threshold for observer-observer (AB) crossing through shared
mediator S in star topology. Below: no AB crossing. Above: AB
crosses 1/4. Receiver noise is fatal (gamma_A > 0.2 kills the
connection); sender noise is tolerable (gamma_B <= 0.5).

**Valid for:** 3-qubit star, Heisenberg, Z-dephasing,
Bell_SA x |+>_B initial state.
**Replaces:** coupling sweep for star-topology crossing threshold.
**Source:** [Star Topology Observers](../experiments/STAR_TOPOLOGY_OBSERVERS.md)

### F30. Gamma-as-signal channel capacity (Tier 2, SVD + Shannon)

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

### F31. Relay protocol MI bound (Tier 2, N=11)

    MI improvement = +83%    (relay + 2:1 coupling vs passive)

Six relay stages, each t_stage = K/gamma. Receiving qubits get
10x noise reduction during their reception phase. Combines three
results: K/gamma timing (F14), quiet receiver (F29),
and 2:1 impedance matching.

**Valid for:** N=11 Heisenberg chain, Z-dephasing, Bell pair initial.
**Replaces:** passive propagation baseline for long chains.
**Source:** [Relay Protocol](../experiments/RELAY_PROTOCOL.md)

### F32. Optimal protection state (Tier 2, N=3)

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

### F34. Qubit necessity equation (Tier 1, proven)

    d^2 - 2*d = 0    -->    d = 0 (nothing)  or  d = 2 (qubit)

Palindromic dephasing requires exactly 2 immune Pauli choices (I, Z)
and 2 decaying choices (X, Y) per site. This fixes d = 2 algebraically.
236 qutrit dissipators tested: 0/236 palindromic.

**Valid for:** any local dephasing model.
**Replaces:** dimensional search for palindrome-compatible systems.
**Source:** [Qubit Necessity](QUBIT_NECESSITY.md)

### F35. Dual-perspective lifetime ratio (Tier 2, hardware-validated)

    t_cross(Pi side) / t_cross(direct) ~ T1 / T2

IBM Torino qubit 52: CΨ_A crosses 1/4 at ~140 us, CΨ_B (Π
perspective) at ~895 us. Factor ~6x. The palindromic partner
decays at the T1 rate, not T2.

**Valid for:** single qubit under T1 + T2 decay, both CΨ
perspectives computed from the same density matrix.
**Replaces:** dual-perspective CΨ simulation.
**Source:** [Both Sides Visible](BOTH_SIDES_VISIBLE.md)

---

## Neural Analog (replace neural symmetry analysis)

### F36. Neural palindrome condition (Tier 1-2, proven + verified)

    Q * J * Q + J + 2*S = 0

    Q = E-I neuron swap operator
    J = Jacobian of Wilson-Cowan dynamics
    S = (1/tau_E + 1/tau_I) / 2 * I

Exact structural analog of quantum palindrome (Pi * L * Pi^-1 =
-L - 2*Σγ * I). Derived algebraically from quantum proof
via E-I swap mapping. C. elegans connectome: residual 0.013 vs
random 0.108 (8x more palindromic than chance).

**Valid for:** Wilson-Cowan neural networks with Dale's Law.
**Replaces:** ad-hoc neural symmetry analysis; connectome
palindromic quality assessment.
**Source:** [Algebraic Palindrome Neural](neural/ALGEBRAIC_PALINDROME_NEURAL.md)

### F37. Neural eigenvalue pairing (Tier 1, from F36)

    mu_k + mu_k' = -(1/tau_E + 1/tau_I)

Analog of lambda + lambda' = -2*Σγ. Every neural mode
pairs with a partner; decay rates sum to the E-I time constant sum.
Verified: mean sum = -0.3012, predicted = -0.300 (1.6% max deviation).

**Valid for:** linearized Wilson-Cowan networks satisfying F36.
**Replaces:** neural eigenvalue computation for pairing verification.
**Source:** [Proof Palindrome Neural](neural/proofs/PROOF_PALINDROME_NEURAL.md)

---

## Derived Relations (follow from combining formulas above)

All derivations verified numerically (N=2-5) against Liouvillian
eigenvalues. Full proofs: [docs/proofs/derivations/](proofs/derivations/).

### D1. Bandwidth and mode density (from F2) [VERIFIED]

    BW = omega_{N-1} - omega_1 = 8*J * cos(pi/N) --> 8*J

    rho(omega) = N / (pi * sqrt(omega * (8*J - omega)))

Van Hove singularities (peaks in the density of states where the
dispersion curve is flat) at band edges. Exact 1D tight-binding
density of states. Max frequency error < 5e-9.

**Valid for:** Heisenberg chain, w=1 sector, large N.
**Replaces:** numerical mode density estimation.

### D2. V-Effect = Q_max / Q_mean (from F6 + F7) [VERIFIED]

    V(N) = Q_max / Q_mean = (1 + cos(pi/N)) / 1

Q_mean = 2*J/gamma exactly. Proof: Sum_{k=1}^{N-1} cos(pi*k/N) = 0
(geometric series identity). Deviation < 3e-15.

### D3. Crossing time ratios (from F27) [VERIFIED]

    t_X / t_Z = K_X / K_Z = ln(2) / (2*ln(1/f*)) = 2.320
    t_depol / t_Z = K_depol / K_Z = 1.177

Verified by Lindblad propagation (Bell+ N=2). Ratio deviation < 3e-6.

### D4. Dimensional factor in crossing (from F12 + F25) [VERIFIED]

    Single qubit (d=2):    f*(1 + f*^2) = 1/2
    Bell+ 2-qubit (d=4):   f*(1 + f*^2) = 3/2 = (d-1) * 1/2

The crossing condition scales with Hilbert space dimension as (d-1)/2.
Exact to machine precision.

### D5. Dynamic palindromic mode count (from F4 + F22 + F23) [VERIFIED]

    Oscillating modes = 4^N - (N+1) - Stat(N)

N=3: 36.  N=5: 898.  Fraction --> 1 exponentially.
**Caveat:** exact at gamma -> 0 only. At finite gamma, the
Hamiltonian mixes weight-parity sectors (w with w +/- 2).

### D6. Spectral gap and mixing time (from Absorption Theorem) [VERIFIED]

    Spectral gap = 2*gamma    (minimum non-zero decay rate)
    Mixing time  <= N*ln(4) / (2*gamma)

**Immediate from Absorption Theorem:** the smallest nonzero ⟨n_XY⟩ for
any eigenmode is bounded below by the weight-1 contribution. For modes
dominated by single-site coherences: ⟨n_XY⟩ → 1, giving gap = 2γ × 1.
The spectral gap is the cost of one X/Y Pauli factor: one absorption
quantum.

Deviation < 1e-14 for all N tested.

### D7. Q-factor distribution (from F2 + F7) [VERIFIED]

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

## Π Operator Properties (from PT-symmetry classification)

### F38. Π squared (Tier 1, proven + verified N=2,3)

    Pi^2 = (-1)^{w_YZ}

Diagonal parity operator in Pauli basis. w_YZ = count of Y,Z entries
in the Pauli string. Pi has order 4 (Pi^4 = I), NOT order 2.
Eigenvalues of Pi^2: +1 (half) and -1 (half), equally split.

In Hilbert space, Pi^2 is realized as conjugation by U = X^{tensor N}
(the global bit-flip): U sigma U = (-1)^{w_YZ} sigma for any Pauli string.
The two definitions agree on the per-site map (I -> I, X -> X, Y -> -Y, Z -> -Z).

**Companion result (F63):** [L, Pi^2] = 0 exactly for all N (proven analytically).
Pi^2 is therefore a conserved quantum number of every Liouvillian eigenmode.

**Valid for:** any N, Z-dephasing Π (P1 family).
**Replaces:** assumption that Π is involutory.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md), [PROOF_BIT_B_PARITY_SYMMETRY](proofs/PROOF_BIT_B_PARITY_SYMMETRY.md)

### F39. det(Π) (Tier 1, proven + verified N=1-4)

    det(Pi) = (-1)^{N * 4^{N-1}}

    N=1: -1.  N >= 2: +1  (since 4^{N-1} is even).

**Valid for:** any N.
**Replaces:** manual determinant computation.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

### F40. Fragile bridge gamma_crit at J_bridge = J (Tier 2, verified)

    gamma_crit = 0.1873  (N=2 per chain, J = J_bridge = 1.0)

Below gamma_crit: all eigenvalues on the imaginary axis (chiral phase).
Above: eigenvalue pairs leave the axis (chiral symmetry breaking, Hopf).
Petermann factor (a measure of eigenstate non-orthogonality that
diverges at exceptional points) peaks at K = 403 at gamma/gamma_crit ~ 1.46 (near-EP).

**Valid for:** N=2 per chain Heisenberg, J_bridge = J = 1.0.
**Replaces:** bisection search at this specific parameter set.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

### F41. Palindromic time (Tier 1, corollary of D10)

    t_Pi = 2*pi / omega_min = pi / (2*J * sin^2(pi/(2*N)))

Period of the slowest palindromic modulation in the SFF. Grows as
~N^2/pi^2 for large N. Confirmed by FFT peak matching (<1% for N=2-4, 6).

**Valid for:** Heisenberg chain, w=1 sector.
**Replaces:** numerical FFT of SFF for modulation period.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### F42. Timescale separation (Tier 2, verified N=2-7)

    t_Pi / t_H ~ (Delta * N^2) / (2 * pi^2 * J)  -->  0  for N -> inf

t_Pi ~ N^2 (polynomial), t_H = 2*pi/Delta ~ 4^N (exponential).
Palindromic modulation is a short-time effect; long-time behavior
is Poisson (integrable). Visibility of modulation ~1/4^N.

| N | t_Pi | t_H | t_Pi/t_H |
|---|------|------|----------|
| 3 | 3.14 | 61.3 | 0.051 |
| 5 | 8.22 | 497 | 0.017 |
| 7 | 15.9 | 5810 | 0.003 |

**Valid for:** Heisenberg chain, Z-dephasing.
**Replaces:** numerical SFF timescale extraction.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### F43. Sector SFF pairing (Tier 1, proven D09)

    K_freq(w, t) = K_freq(N-w, t)    (identical SFF for paired sectors)

Palindromic symmetry Pi maps w -> N-w, so sectors w and N-w have
identical spectral statistics. XOR sector (w=N): K=1.000 (all eigenvalues
degenerate at rate 2*N*gamma).

**Valid for:** Heisenberg chain, Z-dephasing, all N.
**Replaces:** sector-by-sector SFF comparison.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### F44. Crooks-like rate identity (Tier 1, proven D08)

    ln(d_fast / d_slow) = 2 * artanh(Delta_d / (2*Σγ))

For each palindromic pair (d_fast, d_slow) with d_fast + d_slow = 2*Σγ.
Linear approximation: ln(d_fast/d_slow) ~ Delta_d / Σγ,
giving β_eff ~ 1/Σγ (effective inverse temperature).

This is ALGEBRAIC (follows from d_fast + d_slow = 2*Σγ),
NOT a Crooks fluctuation theorem (the thermodynamic identity
relating forward and reverse process probabilities). No Jarzynski
equality holds
(<exp(-Delta_d)> ~ 0.93, not 1). The palindrome has the FORM of
detailed balance without BEING detailed balance.

**Valid for:** any palindromic Liouvillian, all N.
**Replaces:** ad-hoc thermodynamic analogies for the palindrome.
**Source:** [Entropy Production](../experiments/ENTROPY_PRODUCTION.md)

---

## Information Geometry (from Bures metric analysis)

### F45. Bures metric at the fold (Tier 2, N=2 Bell state)

    g(CΨ = 1/4) = 3.36    (Bures metric [the natural Riemannian distance between quantum states based on fidelity], finite, no singularity)

The fold at CΨ = 1/4 has no Riemannian singularity. CΨ is a
smooth coordinate everywhere along the Lindblad trajectory.

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### F46. Geodesic decoherence (Tier 2, N=2 Bell state)

    Geodesic deviation = 9.1e-4    (Lindblad ~ shortest Bures path)

The Lindblad trajectory is approximately geodesic in the Bures metric.
Decoherence follows the geometrically shortest path to equilibrium.
Geometric interpretation of dCPsi/dt < 0 (proven in
[Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)).

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### F47. Gaussian curvature at the fold (Tier 2, N=2)

    K(CΨ = 1/4) = -25    (negative, hyperbolic, finite)

Strong negative curvature at the fold (states diverge quickly).
Finite: no geometric singularity. Decays toward the maximally mixed
state (K → -15 at CΨ ~ 0.2).

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### F48. Pythagorean decomposition (Tier 2, exact at N=2)

    L_c² = L_H² + (L_D + Σγ·I)²

where L_c = L + Σγ·I (centered Liouvillian). The cross term
{L_H, L_D + Σγ·I} vanishes exactly at N=2 because all nonzero L_H
entries connect Pauli strings with w_XY(a) + w_XY(b) = N.

The decomposition: (time evolution)² = (oscillation)² + (cooling)².

**Valid for:** N=2 Heisenberg chain, Z-dephasing, any γ. Exact.
At N≥3: cross term ~2% of ||L_c²||, γ-independent.
**Replaces:** Nothing (new structural insight, not a shortcut).
**Source:** [Primordial Qubit Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md),
[Time Irreversibility Exclusion](proofs/TIME_IRREVERSIBILITY_EXCLUSION.md)

### F49. Cross-term formula (Tier 1, proven)

    R(N) = ||{L_H, L_Dc}|| / (||L_H|| · ||L_Dc||) = √((N-2) / (N · 4^(N-1)))

Equivalently: R(N)² = 4(N-2) / (N · 4^N). At N=2: R = 0 (exact
Pythagorean decomposition). At N=3: R = 1/√48. At N=4: R = 1/√128.

Follows from the key identity ||{L_H, L_Dc}||² = 4γ²(N-2)||L_H||²
(bond-sum rule + spectator variance + disjoint bond supports) and
||L_Dc||² = γ² · 4^N · N.

γ-independent, J-independent, topology-independent. Depends only on N.

**Valid for:** Any shadow-balanced bond coupling (both Paulis in {X,Y}
or both in {I,Z}): Heisenberg XXX, XXZ, XY model, Ising, DM interaction.
Uniform Z-dephasing, any graph, all N >= 2.
**Breaks for:** Shadow-crossing couplings (X_iZ_j, Y_iZ_j).
**Replaces:** per-N numerical computation of the cross-term magnitude.
**Verified:** N=2-6, 4 topologies, 5 gamma values, 10 coupling types.
**Source:** [Proof](proofs/PROOF_CROSS_TERM_FORMULA.md),
[Experiment](../experiments/CROSS_TERM_FORMULA.md)

### F49b. Centered dissipator norm (Tier 1, proven)

    ||L_Dc||² = γ² · 4^N · N

where L_Dc = L_D + Nγ·I is the centered dissipator for uniform
Z-dephasing. Auxiliary lemma used in the proof of F49.

**Valid for:** Uniform Z-dephasing at rate γ per site, any N.
**Source:** [Proof](proofs/PROOF_CROSS_TERM_FORMULA.md), Lemma 1

### F49c. Cross-term formula for shadow-crossing couplings (Tier 1, proven)

    R(N) = √((N-1) / (N · 4^(N-1)))

For couplings where one bond Pauli is in {X,Y} and the other in {I,Z}
(e.g., X_iZ_j, Y_iZ_j). The bond-site variance is 1 (not 0), so
N-2 becomes N-1. Companion to F49.

**Valid for:** Shadow-crossing bond couplings, uniform Z-dephasing,
any graph, all N >= 2.
**Verified:** N=3-6, 5 coupling types, 2 topologies.
**Source:** [Proof](proofs/PROOF_CROSS_TERM_CROSSING.md),
[Experiment](../experiments/CROSS_TERM_CROSSING.md)

---

## Cockpit and Diagnostics (replace full tomography)

### F51. Decoherence cockpit: 3-observable reduction (Tier 2, verified N=2-5, IBM-validated)

    n_eff = 3    (Purity, Concurrence, Ψ-norm)
    coverage = 88-96%    (of trajectory variance, via PCA)
    cost = 3 measurements per pair    (vs 4^N for full tomography)

PCA selects automatically which observable is PC1: Concurrence
dominates in sparse topologies (chains, stars at small N), Purity
dominates in dense topologies (rings, complete, chains at large N).
The effective dimensionality n95 grows as ~N, but the first 3 PCs
always capture 88-96%. Two practical regimes: monitoring (3 PCs,
88-96%) and full diagnostics (~N PCs, 95%).

θ is the most sensitive instrument near the ¼ boundary: 1.68×
amplification over CΨ under sacrifice-zone optimization, because
the arctan mapping amplifies small CΨ changes near ¼.

**Hardware validation:** IBM Torino Q52, CΨ = ¼ crossing predicted
at 114.7 μs, measured at 115.0 μs (0.3% error). Selective DD beats
uniform DD by 3.2× in mutual information on 5-qubit chain.
**Valid for:** Heisenberg chain, Z-dephasing and depolarizing, N=2-5,
9 topologies tested.
**Replaces:** full quantum state tomography for decoherence monitoring.
**Caveat:** Concurrence (PC1 proxy, 57% variance) never validated on
a qubit pair. Single-qubit instruments consistent, 2-qubit untested.
**Source:** [Cockpit Universality](../experiments/COCKPIT_UNIVERSALITY.md)

### F52. Thermal oscillation resilience (Tier 2, verified N=4)

    f_osc(n_bar) ~ 82%    (stable to ±2 modes out of 256, N=4)
    |Delta_osc| <= 4       (for n_bar in [0, 50])

The fraction of oscillating modes is approximately stable under
thermal load. Q_max degrades 16× (68.3 → 4.2), mode count doubles
(47 → 103 distinct frequencies), but the oscillating fraction stays
near 82%. The cavity loses sharpness, not voice.

**Mechanism (April 5, 2026):** L(n_bar) = L_0 + n_bar · L_thermal
is linear in n_bar, so eigenvalues move continuously. The oscillating
count can only change at isolated exceptional-point (EP) crossings
where a real pair splits into a complex conjugate pair (or vice
versa). At N=4: four EP crossings in [0, 2], each affecting exactly
2 eigenvalues. No macroscopic fraction ever changes character.

**Not invariant:** the earlier claim "82% ± 1% invariant" was too
strong. The correct statement: oscillating fraction is stable to
< 1% with isolated EP crossings, not topologically protected. At
n_bar = 50 the fraction drops to 80.5% (−4 modes).

**No phase transition.** Frequency distribution follows neither
Planck nor Stefan-Boltzmann. The mode spectrum is set by the cavity
geometry (J topology), not by thermal statistics.

**Valid for:** N=4 Heisenberg chain, Z-dephasing + amplitude damping.
**Replaces:** thermal Liouvillian sweep to check cavity survival.
The answer is: the cavity survives; no sweep needed.
**Source:** [Thermal Blackbody](../experiments/THERMAL_BLACKBODY.md),
[`thermal_ep_analysis.py`](../simulations/thermal_ep_analysis.py)

---

## Absorption Doses (replace time-evolution for dose estimates)

### F55. Universal absorption dose K_death (Tier 1, proven from D6)

    K_death = ln(10) = 2.303    (dose for 99% absorption)
    K_death / K_fold ~ 2.3      (ratio to CΨ = ¼ crossing dose)
    Immortal modes = N + 1      (zero absorption rate, all N)

**Derivation:** 99% absorption of the slowest mortal mode means
e^{−rate_min · t} = 0.01, so rate_min · t = ln(100). By formula D6:
rate_min = 2γ (spectral gap). Therefore t = ln(100)/(2γ), and
K = γ · t = ln(100)/2 = ln(10) = 2.303. Independent of N, γ,
topology. QED.

N+1 modes have exactly zero absorption rate (pure {I,Z} content,
invisible to the light). Complete absorption is impossible while the
palindrome holds. The cavity always retains light.

**Valid for:** any Heisenberg chain, Z-dephasing, all N.
**Replaces:** time evolution to find "when does the system die."
**Source:** [Trapped Light Localization](../experiments/TRAPPED_LIGHT_LOCALIZATION.md)

---

## Cusp Dynamics (replace iteration counting and trajectory integration)

### F56. Critical slowing iteration count (Tier 1, closed-form, zero fit parameters)

    K(eps, tol) = (1/2)*ln(4*eps/tol) + alpha(tol)*sqrt(eps)

    alpha(tol) = -4 + (1/2)*ln(16*tol)

K = n*sqrt(eps) is the rescaled iteration count of u_{n+1} = u^2 + c
near the cardioid cusp (c = 1/4 - eps). The leading logarithm comes
from saddle-node passage (ODE integral). The -4 comes from the
starting-transient integral (eta_0 = -1/4). The ln(16*tol) term comes
from the Modified Equation correction (Euler discretization error).

Verified: 0.5-2% accuracy over 5 tol decades (10^-8 to 10^-16) and
10 eps decades (10^-1 to 10^-10). Modified Equation slope 0.504 vs
predicted 0.500 (0.8% deviation).

**Valid for:** Mandelbrot iteration near cardioid cusp, any eps > 0,
any tol > 0. Equivalent to CΨ recursion near the 1/4 boundary.
**Replaces:** numerical iteration counting near the saddle-node.
O(1) evaluation instead of O(1/sqrt(eps)) iterations.
**Source:** [Critical Slowing at the Cusp](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md)

### F57. Trajectory dwell time at CΨ = 1/4 (Tier 1, analytical)

    t_dwell(delta) = 2*delta / |dCPsi/dt|_{t_cross}

For Bell+ Z-dephasing (using F25 derivative):

    t_dwell = 1.080088 * delta / gamma    (Bell+ specific)
    K_dwell = gamma * t_dwell = 1.080088 * delta    (gamma-independent)

The prefactor 1.080088 = 2/1.851701 is state-specific (depends on
f_cross and |dCΨ/dt| at the crossing). K_dwell is independent of γ
to machine precision (std < 2 × 10^-17 across γ = 0.1 to 10.0).

**Valid for:** any state with a CΨ = 1/4 crossing under Z-dephasing.
Prefactor is state-specific; γ-invariance of K_dwell is universal.
**Replaces:** trajectory integration for dwell-time estimation.
**Hardware verified:** ibm_kingston (Heron r2), 2026-04-16. Two Bell+ pairs with 2.55x gamma ratio (qubits 124-125, T2=[150,310] us; qubits 14-15, T2=[537,381] us). K_dwell/delta = 0.649 (pair A) and 0.694 (pair B), spread 6.3% despite 2.55x gamma difference. Gamma-invariance of K_dwell confirmed on open quantum hardware. Absolute prefactor 0.67 vs theoretical 1.08 (difference from T1 amplitude damping; the F57 formula assumes pure Z-dephasing, Kingston has T1 comparable to T2). Both CΨ(t) trajectories cross 1/4 monotonically. First two-qubit observation of the CΨ = 1/4 boundary crossing on a quantum computer; the single-qubit case was validated separately on ibm_torino Q80 at 1.9% deviation (F24, IBM Run 3).
**Data:** [data/ibm_cusp_slowing_april2026/](../data/ibm_cusp_slowing_april2026/README.md) (full JSON, PNG, and reanalysis scripts).
**Related:** [CPSI_COMPLEX_PLANE](../experiments/CPSI_COMPLEX_PLANE.md) (the saved density matrices additionally reveal a 2D-spiral structure in the complex c-plane, extending the 1D real-axis picture of BOUNDARY_NAVIGATION).
**Source:** [Critical Slowing at the Cusp](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) (Section 6)

### F58. Weight-based dwell prefactor (Tier 2, even-weight states only)

    prefactor = (2 + 4*W2) / (1 + 6*W2)

    dCPsi/dt = -2*gamma*Psi*(1 + 6*W2)    (weight-based derivative)

W2 is the light-face ({X,Y}) sector weight at the crossing moment.
For Bell+: W2 = 0.3709 (from f_cross = 0.8612), giving prefactor =
1.080088, matching F57 exactly. The static Pauli decomposition
and the dynamic cusp passage are algebraically identical for even-weight
states. For states with odd-weight Pauli content (k = 1 terms), the
prefactor additionally requires individual coefficient magnitudes
(via sqrt(W1)), not just sector weights.

**Valid for:** states with only even-weight Pauli content (Bell+, GHZ,
basis-flip superpositions). Partial for odd-weight states (2.6% error
for |+⟩^{⊗2}).
**Generalized:** F59 extends this to any (W_0, k) pair under the same
two-sector assumption. Bell+ is the W_0 = 1/2, k = 2 special case.
**Replaces:** dynamics solution when sector weights are known at crossing.
**Source:** [Dwell Prefactor from Weights](../experiments/DWELL_PREFACTOR_FROM_WEIGHTS.md)

---

### F59. Generalized dwell prefactor (Tier 1, analytical, verified Bell+ and W_3)

    prefactor = (4/k) * (W_0 + W_k) / (W_0 + 3*W_k)

For any state whose Pauli content lives in exactly two sectors (a stationary
sector at weight W_0, immune to dephasing, and a single coherent sector at
XY-weight k > 0 with weight W_k), the dwell-time prefactor at the CPsi = 1/4
crossing is a pure algebraic function of the sector weights at the crossing
moment. Independent of N and d.

Reduces to F58 (Bell+) when W_0 = 1/2, k = 2:

    (4/2) * (1/2 + W_2) / (1/2 + 3*W_2) = (2 + 4*W_2) / (1 + 6*W_2)

The factor 3 in the denominator is structural (product rule: 2*W_k from the
dC/dt contribution plus W_k from the C*dPsi/dt contribution), NOT from d - 1
as initially conjectured in the April 5 dwell-prefactor work. The d - 1
factor cancels because both Psi and dPsi/dt carry it.

Verified on two states with different stationary weights:

| State | N | W_0 | k | Prefactor (formula) | Prefactor (direct) |
|-------|---|-----|---|--------------------:|-------------------:|
| Bell+ | 2 | 1/2 | 2 | 1.080088            | 1.080097           |
| W_3   | 3 | 1/3 | 2 | 0.876832            | 0.876839           |

Both methods agree to within 10^-5 relative error (limited by the
finite-difference resolution of the direct measurement, not by the formula).

**Valid for:** states with Pauli content in exactly two sectors, with linear
Psi-to-f relation. Includes W states, basis-flip superpositions between
computational basis states at even Hamming distance, and any state
constructed from a single off-diagonal coherence pattern.
**Breaks for:** states populating multiple coherent sectors simultaneously
(|+⟩^{⊗2} with k=1 and k=2 content); states with odd-weight Pauli content
where the linear Psi-to-f relation fails.
**Replaces:** Lindblad propagation for dwell-time estimation in the two-sector
class. O(1) algebraic evaluation given W_0 and W_k at the crossing.
**Source:** [Generalized Dwell Prefactor](../experiments/DWELL_PREFACTOR_GENERALIZED.md)

### F60. GHZ_N born below the fold (Tier 1, geometric corollary)

    CPsi(0) for GHZ_N = 1 / (2^N - 1)

For GHZ_N = (|0...0> + |1...1>) / sqrt(2), the initial value of CPsi drops
below 1/4 as soon as N >= 3, regardless of any dynamics:

| N | CPsi(0)        | Above 1/4?              |
|---|----------------|-------------------------|
| 2 | 1/3 = 0.3333   | Yes (Bell+ crosses)     |
| 3 | 1/7 = 0.1429   | No (already classical)  |
| 4 | 1/15 = 0.0667  | No                      |
| 5 | 1/31 = 0.0323  | No                      |

Derivation: C(0) = 1 (pure state), and the only nonzero off-diagonal matrix
elements are rho[0...0, 1...1] = 1/2 and its conjugate, giving L1 coherence
= 1 and Psi(0) = L1/(d - 1) = 1/(2^N - 1). Therefore CPsi(0) = 1/(2^N - 1).

For N >= 3, GHZ_N starts in the classical regime of the R = CPsi^2 framework
before any dephasing acts on it. This is gamma-independent: reducing the
dephasing rate cannot fix the geometric deficit. The only escape is to
change the state.

This geometric statement complements F22 (GHZ projects 100% of its coherent
weight onto the fastest-decaying XOR modes at rate 2*Sigma*gamma). The two
arguments converge: GHZ encoding is structurally unsuitable for state
transfer, regardless of how slowly it would decay.

**Valid for:** GHZ_N for any N >= 2, gamma-independent.
**Replaces:** trajectory simulation to confirm GHZ_N (N >= 3) never crosses 1/4.
**Source:** [Generalized Dwell Prefactor](../experiments/DWELL_PREFACTOR_GENERALIZED.md) Section 4,
[Engineering Blueprint Rule 1](../publications/ENGINEERING_BLUEPRINT.md) Note (April 6, 2026)

### F61. n_XY Parity Selection Rule (Tier 1, proven, verified 64 configs N=2-6)

The Liouvillian of the isotropic Heisenberg model under Z-dephasing
preserves n_XY parity. Consequences:

1. Every eigenmode has definite n_XY parity (even or odd).
2. Every single-excitation density matrix has purely even n_XY content.
3. No SE state can excite an odd-n_XY eigenmode. The overlap is exactly zero.

**Corollary (Accessibility Boundary).** If the slowest SE-accessible mode
has rate alpha_1, and a slower mode exists at rate alpha_2 < alpha_1 with
odd n_XY parity, then alpha_1 is an exact ceiling for SE protection. No
optimizer within the SE family can reach alpha_2.

**Companion symmetry (F63):** L also commutes with Pi^2 (the w_YZ-parity
operator). Together these are the two independent Z2 symmetries of L
admitted by the Pauli algebra of d=2: bit_a (n_XY parity, this rule) and
bit_b (w_YZ parity, F63). Combined, they decompose the operator space into
4 sectors, each of dimension 4^(N-1).

**Valid for:** any isotropic Heisenberg (XX+YY+ZZ) or XY (XX+YY) model,
any graph topology, any site-dependent Z-dephasing gamma_k, any N.
**Breaks for:** amplitude damping (T1), transverse fields (odd-n_XY terms).
**Verified:** 64 configurations (N=2-6, Chain/Star/Ring/Complete, 4 gamma
profiles). Second slow mode SE Frobenius ratio < 1e-3 in all 64. Data:
`simulations/results/lens_survey/lens_survey_scaling.txt`.
**Source:** [Proof](proofs/PROOF_PARITY_SELECTION_RULE.md)

### F62. CΨ(0) for W_N (Tier 1, analytical, verified N=2-10)

    CΨ(0) = 2(N^2 - 4N + 8) / (3N^3)

For the W_N state (equal superposition of single excitations) on any
N-qubit system, the initial CΨ on any pair is given by the formula above.
The reduced density matrix for any pair (a,b) is:

    rho_ab = diag((N-2)/N, 1/N, 1/N, 0) + (1/N)|01><10| + (1/N)|10><01|

from which Tr(rho_ab^2) = (N^2 - 4N + 8)/N^2, L1 = 2/N, Psi = 2/(3N).

| N | CΨ(0) | Fraction | Above 1/4? |
|---|-------|----------|------------|
| 2 | 1/3 = 0.3333 | 8/24 | Yes (W_2 = Bell+) |
| 3 | 10/81 = 0.1235 | | No |
| 5 | 26/375 = 0.0693 | | No |
| 10 | 68/1500 = 0.0453 | | No |

**Corollary (W_N born below the fold).** For N >= 3, CΨ(0) < 1/4.
Proof: 2(N^2 - 4N + 8)/(3N^3) < 1/4 iff 3N^3 - 8N^2 + 32N - 64 > 0.
At N=3 this evaluates to 41 > 0, and the cubic is monotonically increasing
for N >= 3. Combined with the Parity Selection Rule (F61), this proves
that single-excitation states on Heisenberg chains under Z-dephasing
never cross CΨ = 1/4. The cusp exit is structurally inaccessible to them.

**Valid for:** W_N on any N-qubit system, pair-independent (permutation symmetry).
**Verified:** numerically at N=2-10, all within machine precision.
**Source:** `simulations/cpsi_wn_analytical.py`,
[Cusp-Lens Connection](../experiments/CUSP_LENS_CONNECTION.md)

### F63. w_YZ Parity Symmetry, [L, Pi^2] = 0 (Tier 1, proven analytically, verified N=2-5)

    [L, Pi^2_super] = 0     (exactly, for all N)

The Liouvillian commutes exactly with the Pi^2-parity superoperator
Pi^2_super(rho) = U rho U where U = X^{tensor N} (global bit-flip).
Equivalently in Pauli basis: Pi^2 acts as multiplication by (-1)^{w_YZ},
and L preserves w_YZ parity.

**Six-line proof.** U conjugation acts per site as I -> I, X -> X,
Y -> -Y, Z -> -Z. Heisenberg bonds XX/YY/ZZ each contain two Y's or two
Z's (or none), so signs cancel and U H U = H. Z-dephasing dissipator is
quadratic in Z (term Z rho Z), where the two minus signs from U Z U = -Z
also cancel, and the anti-commutator term {Z dagger Z, rho} = {I, rho}
is trivially U-invariant.

**Two independent Z2 symmetries.** Together with F61 (n_XY parity, bit_a),
this gives L two independent Z2 symmetries proven for all N. Per-site Pauli
{I, X, Y, Z} factorizes as C2 x C2 indexed by (bit_a, bit_b) =
(n_XY, w_YZ). The 4-block decomposition has dimension 4^(N-1) per block.
This is the maximal symmetry decomposition admitted by the Pauli algebra
of d=2 (no third independent Z2 classification exists per F34/QUBIT_NECESSITY).

**Per-sector mode count (closed form).** For Heisenberg coupling with Z-dephasing on the boundary qubit B, modes split within each Pi^2-sector by the Absorption Theorem applied at site B:

    conserved per sector:  even = floor(N/2) + 1,  odd = ceil(N/2)
    correlation per sector: same as conserved (palindrome symmetry)
    mirror per sector:     2^(2N-1) - 2 * (conserved per sector)

Mechanism: conserved modes are the (N+1) elementary symmetric polynomials e_d(Z_1, ..., Z_N) for d=0..N (functions of S_z, commuting with both H and Z_B). Each e_d has w_YZ-parity = d mod 2. The asymmetry for even N (one extra even-parity conserved) comes from e_N having even parity when N is even.

| N | sector | cons (e, o) | mirror (e, o) |
|---|--------|-------------|---------------|
| 2 | 8      | (2, 1)      | (4, 6)        |
| 3 | 32     | (2, 2)      | (28, 28)      |
| 4 | 128    | (3, 2)      | (122, 124)    |
| 5 | 512    | (3, 3)      | (506, 506)    |

**Valid for:** Heisenberg (XX+YY+ZZ), XY (XX+YY); Z-dephasing on any subset
of sites; any N; any graph; uniform or non-uniform gamma_k.
**Breaks for:** single-site Y or Z terms in H (transverse field, magnetic
field along Z); Y or X jump operators (no two-factor cancellation).
**Verified:** ||[L, Pi^2]|| = 0.000000e+00 (identically zero, not numerically
small) at N=2, 3, 4, 5. Also for Heisenberg XXX with uniform gamma at N=3.
Per-sector mode count formula verified at N=2-5; conserved-modes-as-e_d(Z) verified at N=2-4.
Data: `simulations/primordial_bit_a_bit_b_N_scaling.py`,
`simulations/mirror_mode_split_formula.py`.
**Source:** [PROOF_BIT_B_PARITY_SYMMETRY](proofs/PROOF_BIT_B_PARITY_SYMMETRY.md),
[PRIMORDIAL_QUBIT](../hypotheses/PRIMORDIAL_QUBIT.md) Section 9

### F64. Effective gamma from cavity mode exposure (Tier 1-2, analytical + verified N=3,4)

For an N-qubit chain with XX+YY coupling and Z-dephasing only on the outermost site B at rate γ_B, the effective dephasing rate of the slowest eigenmode contributing to inner-site S coherence is:

    γ_eff = γ_B · |a_B|²        (decoherence rate, Lorentzian half-width)

equivalently in Liouvillian-eigenvalue units:

    α = 2γ_B · |a_B|²           (Liouvillian decay constant, α = -Re(λ))

where a_B is the B-site amplitude of the single-excitation Hamiltonian eigenvector. The factor of 2 between the two forms is the standard QM convention: ρ_{ij}(t) ∝ exp(-γ_eff·t) corresponds to a Liouvillian eigenvalue λ = -2γ_eff. Both express the same content; choose the convention that fits the surrounding context. This is the Absorption Theorem (F1/AT) applied to the single-excitation sector: α = 2γ_B · ⟨n_XY⟩_B, with ⟨n_XY⟩_B = |a_B|² verified to machine precision for these modes ([`factor_two_clarification.py`](../simulations/factor_two_clarification.py)).

γ_B appears as a constant prefactor. It is not diminished by intervening sites.

**Closed form at N=3** (chain S-M-B, couplings J_SM, J_MB, r = J_SM/J_MB):

                 ⎧ r² / (r² + 1)       for r < 1/√2    [zero mode]
    g(r) =       ⎨
                 ⎩ 1 / (2(r² + 1))     for r ≥ 1/√2    [bonding mode]

Crossover at r = 1/√2, g = ⅓. Special value: g(1) = ¼.

Derived from the 3×3 single-excitation Hamiltonian eigenvalues {0, ±√(J_SM² + J_MB²)} and eigenvectors.

**General N.** Diagonalize the N×N tridiagonal single-excitation Hamiltonian. Find the eigenvector with the smallest |a_B|² among those with nonzero |a_S|². No layered composition: the formula is a global eigenvector property, not a product of per-layer factors. Multiplicative stacking fails at N=4 (ratio 0.04 to 62); eigenvector formula exact (ratio 1.0000 ± 0.0003).

**Replaces:** time-domain exponential fit for γ_eff extraction.
**Valid for:** XX+YY chains with Z-dephasing on one end site; good-cavity regime (γ_B ≪ J_MB). Breaks when γ_B ≥ J_MB (bad cavity: B decoheres before transmitting).
**Verified:** N=3 (max relative error 1.8% vs 64×64 Liouvillian), N=4 (9 configs, ratio 1.0000 ± 0.0003 vs 256×256 Liouvillian).
**Scripts:** [`primordial_gamma_analytical.py`](../simulations/primordial_gamma_analytical.py), [`primordial_gamma_stacking_4qubit.py`](../simulations/primordial_gamma_stacking_4qubit.py), [`factor_two_clarification.py`](../simulations/factor_two_clarification.py)
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md)

### F65. Single-excitation spectrum of uniform open XX chain (Tier 1, proven, verified N=3..30)

For the uniform open XX chain (all couplings J, N sites) with Z-dephasing at rate γ₀ on endpoint B = site N-1, the single-excitation dissipation rates are:

    α_k / γ₀ = (4 / (N+1)) · sin²(kπ / (N+1)),    k = 1, ..., N

This is F64 evaluated on the analytically known eigenvectors ψ_k(i) = √(2/(N+1)) · sin(πk(i+1)/(N+1)) of the N×N tridiagonal single-excitation Hamiltonian. The endpoint amplitude |ψ_k(N-1)|² = (2/(N+1)) · sin²(kπ/(N+1)), and the Absorption Theorem gives α_k = 2γ₀ · |a_B|².

**Properties:**
- All α_k lie in [0, 2γ₀].
- Internal symmetry: α_k = α_{N+1-k}, from sin²(kπ/(N+1)) = sin²((N+1-k)π/(N+1)). This mirror is within the single-excitation spectrum; the palindromic pairing α_a + α_b = 2γ₀ of F1 in general maps single-excitation modes to multi-excitation sectors.
- Maximum rate: α_max / γ₀ = 4/(N+1) when N is odd (then k = (N+1)/2 is integer and sin² = 1 is attained exactly); when N is even the maximum lies strictly below 4/(N+1). The single-excitation sector never reaches 2γ₀ for N ≥ 2; its maximum decays as 4/(N+1) → 0 for growing N.

**Niven rationality.** All α_k/γ₀ are rational if and only if N+1 ∈ {1, 2, 3, 4, 6}, i.e., N ∈ {0, 1, 2, 3, 5}. For all other N the values are algebraic irrationals (golden-ratio family at N=4,9; √2 family at N=7; general cyclotomic otherwise). This follows from Niven's theorem: sin²(rπ) is rational only for sin(rπ) ∈ {0, ±1/2, ±1}.

**Verified values:**
- N=3: α/γ₀ ∈ {1/2, 1, 1/2}
- N=4: α/γ₀ ∈ {0.276393, 0.723607, 0.723607, 0.276393} (algebraic irrationals from sin²(π/5), sin²(2π/5); golden-ratio family)
- N=5: α/γ₀ ∈ {1/6, 1/2, 2/3, 1/2, 1/6}

**Verified:** Formula matches the tridiagonal N×N single-excitation eigendecomposition to machine precision (max error 1.2 · 10⁻¹⁵) for N=3..30. All single-excitation rates confirmed present (to within O((γ₀/J)²) perturbative corrections, see below) in the full 4^N Liouvillian spectrum for N=3..7. Dynamical check at γ₀ = 0.01, where second-order shifts are ~10⁻⁶: formula predicts the decay rate of coherence operators ρ_k = |ψ_k⟩⟨0| under full Liouvillian propagation to within 10⁻⁴ relative error for all k at N=5. Asymptotic 1/(N+1)³ scaling of α_min verified; ratio to 4π²/(N+1)³ rises monotonically from 0.81 at N=3 to 0.99 at N=15.

**Perturbative nature.** The formula is derived by applying the Absorption Theorem (AT) to single-excitation coherence operators |ψ_k⟩⟨vac|, treating them as decoupled Liouvillian right eigenvectors. This is exact to first order in γ₀/J. At finite γ₀ the Lindblad dissipator mixes |ψ_k⟩⟨vac| with other sectors, and the full-Liouvillian eigenvalue shifts by O((γ₀/J)²) relative to the formula. For γ₀/J = 0.05 and N=5, the relative shift is ≈ 4·10⁻³ (verified via full eigendecomposition in `palindromic_partner_f67.py`). The palindromic pairing F1 survives this shift exactly: α_b + α_p = 2γ₀ to machine precision, even as each individual rate deviates from its first-order value (see F68).
**Scripts:** [`single_excitation_spectrum.py`](../simulations/single_excitation_spectrum.py), [`f65_dynamic_verification.py`](../simulations/f65_dynamic_verification.py)
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md)

### F66. Pole modes at α = 0 and α = 2γ₀ (Tier 1, verified XY chain with B at endpoint, N=3..7)

For the uniform open XY chain with single-site Z-dephasing at the endpoint site B = N-1, the dissipation interval [0, 2γ₀] has exact eigenvalues at both endpoints:

- **α = 0 modes:** ⟨n_XY⟩_B = 0 (forced by the Absorption Theorem α = 2γ₀·⟨n_XY⟩_B). Dominant Pauli strings have I or Z at every site (total XY-weight = 0). Z-basis population content, completely shielded from Z-dephasing at B.

- **α = 2γ₀ modes:** ⟨n_XY⟩_B = 1 (same mechanism). Dominant Pauli strings have X or Y at every site (total XY-weight = N). Maximally off-diagonal at B, fully exposed to Z-dephasing.

The two poles are palindromic partners under the conjugation Π, which maps total XY-weight w ↔ N-w (see F43). The single-excitation sector (F65) never reaches either pole for N ≥ 3; both poles live in the extreme XY-weight sectors (w = 0 and w = N).

**Multiplicity:** exactly N+1 at each pole, verified for N=3..7. Each α = 0 mode corresponds to one of the N+1 elementary symmetric polynomials e_d(Z_1, ..., Z_N) in F63 (commuting with both H and Z_B). The α = 2γ₀ sector has matching multiplicity by Π-symmetry.

**Scope.** Verified only for the uniform XY chain with B at the endpoint. Whether the same structure (existence of both poles, multiplicity N+1) persists for other topologies (ring, star, Y-junction) or for interior B-positions is open. Indirect evidence from the structure-points scan: at B = center of N=5 chain, α = 0 has multiplicity 64 (not 6), so the N+1 count is endpoint-specific.

**Verified:** ⟨n_XY⟩_B = 1.000000 exact for all α = 2γ₀ modes (N=3..5, from Pauli basis projection). Dominant Pauli strings have total XY-weight N for α = 2γ₀ modes and total XY-weight 0 for α = 0 modes (N=3, N=4 explicit). Multiplicity N+1 at each pole verified for N=3..7. Dynamical check of F63 conservation: all N+1 elementary symmetric polynomials e_d(Z_1,...,Z_N) drift by < 10⁻¹⁴ under Lindblad evolution for N=4 over 80 time units, while the non-symmetric control Z_0 Z_2 drifts by 3 × 10⁻². Confirms the conserved observables at the α = 0 pole are precisely the e_d, not arbitrary Z-products.
**Scripts:** [`two_gamma_pole.py`](../simulations/two_gamma_pole.py), [`f65_dynamic_verification.py`](../simulations/f65_dynamic_verification.py)
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) section "The dissipation interval [0, 2γ₀]", [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md)

### F67. Bonding-mode encoding is the optimal dephasing-protected Bell pair (Tier 1, verified N=3, N=5)

For an isolated external reference qubit R entangled with an N-site uniform XY chain C under endpoint Z-dephasing γ₀ at Q_{N-1}, the delocalized encoding

    |Ψ⟩ = (|0⟩_R |vac⟩_C + |1⟩_R |ψ_1⟩_C) / √2

with |ψ_1⟩ = √(2/(N+1)) · Σᵢ sin(π(i+1)/(N+1)) |1ᵢ⟩ (the k=1 bonding mode of F65) decays as a pure exponential from t=0:

    N(R:C)(t) = N(0) · exp(-α_1 · t),   α_1 = (4γ₀/(N+1)) sin²(π/(N+1))

because |vac⟩⟨ψ_1| is a Liouvillian right eigenvector with eigenvalue -α_1 (Absorption Theorem applied to F65, ⟨n_XY⟩_B = |a_B(ψ_1)|² = (2/(N+1)) sin²(π/(N+1))). The R half of the Bell pair is isolated, so the only dissipation channel is the chain's slowest mode.

**Comparison with localized encodings.** Bell pairs |Φ+⟩_{R,Q_j} ⊗ |0…0⟩_rest at any chain site j are superpositions of all k-modes, |1_j⟩ = Σ_k U_{jk} |ψ_k⟩, and decay multi-exponentially: an initial fast transient from high-α_k components, then an asymptotic tail at α_1 (the slowest nonzero Liouvillian eigenvalue). The tail survives with amplitude |U_{j,1}|² = (2/(N+1)) sin²(π(j+1)/(N+1)).

**Counterintuitive equivalence A ≡ C.** By palindromic symmetry of the uniform chain, |U_{0,k}|² = |U_{N-1,k}|², so inner-localized (j=0) and outer-localized (j=N-1) encodings have **identical** decay dynamics despite their very different spatial relationships to the dephased site B=N-1. Spatial distance from noise is not the protecting mechanism; spectral encoding is. The bonding mode is optimal because it has least ⟨n_XY⟩_B of all single-excitation modes.

**Scaling.** T_2 ≡ 1/α_1 = (N+1)/(4γ₀ · sin²(π/(N+1))) → (N+1)³/(4π²γ₀) for large N. Cubic improvement with chain length, no saturation regime identified.

**Verified:** Variant B (bonding-mode) α_fit/α_1 = 0.9989 (N=3) and 0.9963 (N=5), both within 10⁻³. Variant A (inner-localized) long-time tail α_fit/α_1 = 1.046 (N=3) and 1.015 (N=5), within 5% as expected for multi-exponential. Variant C (outer-localized) yields 1.023 (N=3) and 1.008 (N=5), differing from A only by fit noise (confirms palindromic equivalence). At fixed decay rate, bonding-mode preserves ~2.2× more entanglement than either localized variant at t = 0.4 · T_2 in the N=5 run.
**Scripts:** [`bell_pair_chain_protection.py`](../simulations/bell_pair_chain_protection.py)
**Source:** F65, F66, [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md), [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)

### F68. Palindromic partner of the bonding mode (Tier 1, verified N=3, 4, 5)

By F1 (palindromic spectrum), the bonding-mode eigenvalue -α_b of the full Liouvillian has a partner at -α_p with α_b + α_p = 2γ₀ exactly. Direct eigendecomposition of the 4^N-dim Liouvillian (N=3, 4, 5; γ₀ = 0.05, J = 1.0) confirms this pairing to machine precision (|α_b + α_p - 2γ₀| < 2·10⁻¹⁵), and identifies the partner's structure:

**Sector.** For N ≥ 4, the partner mode lives entirely in the XY-weight-(N-1) sector, while the bonding mode lives entirely in w=1. These are Π-mirror sectors under the total XY-weight conjugation w ↔ N-w (F43). ⟨n_XY⟩_B = α_p / (2γ₀) ≈ 0.86 (N=4), 0.92 (N=5), close to but distinct from the F66 pole at ⟨n_XY⟩_B = 1.

**Rank.** For N ≥ 4, the partner's right eigenvector V_p is a rank-1 operator (SVD: σ₁/σ₀ = 0 to numerical precision). The partner therefore admits a pure-state factorization V_p = |u⟩⟨v|, and can be operationally prepared as the off-diagonal coherence of an R-C Bell-pair-like state: ρ₀ has off-diagonal |0⟩⟨1|_R ⊗ |u⟩⟨v|, decaying at rate α_p = 2γ₀ - α_b.

**N=3 special case.** At N=3, bonding and partner are both rank-2 and mix adjacent w-sectors (bonding: w=0 + w=2; partner: w=1 + w=3). This is a degeneracy artifact: multiplicities are 4 each, and any orthonormal basis of a 4-dim degenerate eigenspace is a valid eigenvector set. The palindromic pairing α_b + α_p = 2γ₀ nonetheless holds exactly.

**Palindrom robustness.** F65 predicts α_b from first-order perturbation in γ₀/J; at finite γ₀ the full-Liouvillian eigenvalue shifts by O((γ₀/J)²). Both α_b and α_p track this shift, but in opposite directions: their sum is pinned to 2γ₀ by F1 regardless of perturbative order. This is a concrete demonstration that F1 is algebraically exact while F65 is perturbative; the palindrome is the more robust structural feature.

| N | α_b (F65 formula) | α_b (full L) | α_p (full L) | α_b + α_p | error from 2γ₀ |
|---|-------------------|--------------|--------------|-----------|----------------|
| 3 | 0.025000 | 0.025003 | 0.074997 | 0.100000 | 4.2·10⁻¹⁶ |
| 4 | 0.013820 | 0.013784 | 0.086216 | 0.100000 | 2.4·10⁻¹⁶ |
| 5 | 0.008333 | 0.008303 | 0.091697 | 0.100000 | 1.1·10⁻¹⁵ |

**Operational consequence.** The partner mode is a Bell-pair encoding that decays at 2γ₀ - α_b, the fastest non-pole rate. It is the "dark twin" of the F67 bonding-mode encoding: identical spatial structure under Π-mirror, opposite exposure. Where F67 gives cubic T₂ scaling (N+1)³/(4π²γ₀), the partner encoding gives T₂ → 1/(2γ₀) for large N (maximum exposure, saturating at the F66 pole in the limit).

**Verified:** α_b + α_p + 2γ₀ = 0 to <2·10⁻¹⁵ for N=3, 4, 5. Rank-1 verified by SVD (σ₁/σ₀ < 10⁻¹⁵) for N=4, 5. Pauli-basis decomposition shows bonding w=1 dominance (100% for N≥4) and partner w=N-1 dominance (100% for N≥4); mixed-sector at N=3 (degeneracy artifact).
**Scripts:** [`palindromic_partner_f67.py`](../simulations/palindromic_partner_f67.py)
**Source:** F1, F43, F65, F66, F67, [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md)

---

*Each formula in this document is a Liouvillian that does not need
to be built.*
