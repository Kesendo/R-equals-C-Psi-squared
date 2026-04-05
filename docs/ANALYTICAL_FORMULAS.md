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

**Formulas 3, 8, 33, and D6 are corollaries** of this theorem.
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

### 3. Decay rate bounds (Tier 1, corollary of Absorption Theorem)

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

### 4. Stationary mode count (Tier 1, Clebsch-Gordan decomposition)

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

### 50. Weight-1 degeneracy / conserved operator count (Tier 1, proven + verified N=2-7)

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

### 6. V-Effect gain (Tier 1-2, verified N=2-6)

    V(N) = 1 + cos(pi/N) = 2*cos^2(pi/(2N))

Q-factor amplification from coupling. γ-independent (cancels in ratio).
For N=5: (5+sqrt(5))/4 = 1.80902. For N→∞: V = 2 (saturation).
Under non-uniform γ: applies only to the extremal (best-Q) mode.

**Valid for:** Heisenberg chain, Z-dephasing, all N.
**Physically validated:** proton water chain N=1-5, machine-precision
match ([Proton Water Chain](../experiments/PROTON_WATER_CHAIN.md)).
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

### 8. 2× universal decay law (Tier 1, corollary of Absorption Theorem)

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

### 11. Mode localization profile, N=5 (Tier 2, geometric)

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
**Replaces:** formula 12 when T1 is finite; numerical CΨ(t)
simulation for superconducting qubits.
**Source:** [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md)

### 25. CΨ closed form, Bell+ Z-dephasing (Tier 1, proven)

    CPsi(t) = f * (1 + f^2) / 6,       f = e^{-4*gamma*t}

    dCPsi/dt = -2*gamma*f*(1 + 3*f^2) / 3

Crossing at f* = 0.8612 (from f*(1 + f*^2) = 3/2).
K = gamma * t_cross = 0.0374.

**Valid for:** Bell+ initial state, Z-dephasing, 2 qubits.
**Replaces:** numerical integration for CΨ(t) trajectory.
O(1) evaluation instead of ODE solver.
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### 26. CΨ closed form, general Pauli channels (Tier 1, proven)

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

### 27. K values per noise channel (Tier 1, from formula 26)

    K_Z     = 0.0374    (pure Z-dephasing)
    K_X     = 0.0867    (pure X-noise)
    K_Y     = 0.0867    (pure Y-noise)
    K_depol = 0.0440    (depolarizing, gamma/3 each axis)

Complements formula 14 (K per bridge metric). These are K per
noise TYPE, all measured with CΨ on Bell+ state.

**Valid for:** Bell+ state, single-axis or depolarizing noise.
**Replaces:** per-channel crossing time derivation.
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### 28. Fixed-point absorber theorem (Tier 1-2)

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

    Σγ_crit / J = 0.00249 (Bell state)
    Σγ_crit / J = 0.00497 (product state)

Below: no fold, CΨ oscillates forever. Above: CΨ crosses 1/4
irreversibly. Max/min ratio across N=2-5: 1.015 (1.5% variation).

**Valid for:** Heisenberg chain, Z-dephasing, N=2-5.
**Replaces:** γ sweep to find fold onset.
**Source:** [Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)

### 19. Fragile bridge asymptotic constant (Tier 2)

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

IBM Torino qubit 52: CΨ_A crosses 1/4 at ~140 us, CΨ_B (Π
perspective) at ~895 us. Factor ~6x. The palindromic partner
decays at the T1 rate, not T2.

**Valid for:** single qubit under T1 + T2 decay, both CΨ
perspectives computed from the same density matrix.
**Replaces:** dual-perspective CΨ simulation.
**Source:** [Both Sides Visible](BOTH_SIDES_VISIBLE.md)

---

## Neural Analog (replace neural symmetry analysis)

### 36. Neural palindrome condition (Tier 1-2, proven + verified)

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

### 37. Neural eigenvalue pairing (Tier 1, from formula 36)

    mu_k + mu_k' = -(1/tau_E + 1/tau_I)

Analog of lambda + lambda' = -2*Σγ. Every neural mode
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

Van Hove singularities (peaks in the density of states where the
dispersion curve is flat) at band edges. Exact 1D tight-binding
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

### D6. Spectral gap and mixing time (from Absorption Theorem) [VERIFIED]

    Spectral gap = 2*gamma    (minimum non-zero decay rate)
    Mixing time  <= N*ln(4) / (2*gamma)

**Immediate from Absorption Theorem:** the smallest nonzero ⟨n_XY⟩ for
any eigenmode is bounded below by the weight-1 contribution. For modes
dominated by single-site coherences: ⟨n_XY⟩ → 1, giving gap = 2γ × 1.
The spectral gap is the cost of one X/Y Pauli factor: one absorption
quantum.

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

## Π Operator Properties (from PT-symmetry classification)

### 38. Π squared (Tier 1, proven + verified N=2,3)

    Pi^2 = (-1)^{w_YZ}

Diagonal parity operator in Pauli basis. w_YZ = count of Y,Z entries
in the Pauli string. Pi has order 4 (Pi^4 = I), NOT order 2.
Eigenvalues of Pi^2: +1 (half) and -1 (half), equally split.

**Valid for:** any N, Z-dephasing Π (P1 family).
**Replaces:** assumption that Π is involutory.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

### 39. det(Π) (Tier 1, proven + verified N=1-4)

    det(Pi) = (-1)^{N * 4^{N-1}}

    N=1: -1.  N >= 2: +1  (since 4^{N-1} is even).

**Valid for:** any N.
**Replaces:** manual determinant computation.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

### 40. Fragile bridge gamma_crit at J_bridge = J (Tier 2, verified)

    gamma_crit = 0.1873  (N=2 per chain, J = J_bridge = 1.0)

Below gamma_crit: all eigenvalues on the imaginary axis (chiral phase).
Above: eigenvalue pairs leave the axis (chiral symmetry breaking, Hopf).
Petermann factor (a measure of eigenstate non-orthogonality that
diverges at exceptional points) peaks at K = 403 at gamma/gamma_crit ~ 1.46 (near-EP).

**Valid for:** N=2 per chain Heisenberg, J_bridge = J = 1.0.
**Replaces:** bisection search at this specific parameter set.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

### 41. Palindromic time (Tier 1-2, from formula 2)

    t_Pi = 2*pi / omega_min = pi / (2*J * sin^2(pi/(2*N)))

Period of the slowest palindromic modulation in the SFF. Grows as
~N^2/pi^2 for large N. Confirmed by FFT peak matching (<1% for N=2-4, 6).

**Valid for:** Heisenberg chain, w=1 sector.
**Replaces:** numerical FFT of SFF for modulation period.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### 42. Timescale separation (Tier 2, verified N=2-7)

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

### 43. Sector SFF pairing (Tier 2, verified N=3-5)

    K_freq(w, t) = K_freq(N-w, t)    (identical SFF for paired sectors)

Palindromic symmetry Pi maps w -> N-w, so sectors w and N-w have
identical spectral statistics. XOR sector (w=N): K=1.000 (all eigenvalues
degenerate at rate 2*N*gamma).

**Valid for:** Heisenberg chain, Z-dephasing, all N.
**Replaces:** sector-by-sector SFF comparison.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### 44. Crooks-like rate identity (Tier 2, algebraic)

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

### 45. Bures metric at the fold (Tier 2, N=2 Bell state)

    g(CΨ = 1/4) = 3.36    (Bures metric [the natural Riemannian distance between quantum states based on fidelity], finite, no singularity)

The fold at CΨ = 1/4 has no Riemannian singularity. CΨ is a
smooth coordinate everywhere along the Lindblad trajectory.

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### 46. Geodesic decoherence (Tier 2, N=2 Bell state)

    Geodesic deviation = 9.1e-4    (Lindblad ~ shortest Bures path)

The Lindblad trajectory is approximately geodesic in the Bures metric.
Decoherence follows the geometrically shortest path to equilibrium.
Geometric interpretation of dCPsi/dt < 0 (proven in
[Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)).

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### 47. Gaussian curvature at the fold (Tier 2, N=2)

    K(CΨ = 1/4) = -25    (negative, hyperbolic, finite)

Strong negative curvature at the fold (states diverge quickly).
Finite: no geometric singularity. Decays toward the maximally mixed
state (K → -15 at CΨ ~ 0.2).

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### 48. Pythagorean decomposition (Tier 2, exact at N=2)

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

### 49. Orthogonality constant at N=3 (Tier 2, verified)

    ||{L_H, L_D + Σγ·I}|| / (||L_H|| · ||L_D + Σγ·I||) = 1/√48

γ-independent. Geometric constant of the Heisenberg chain at N=3.

**Valid for:** N=3 Heisenberg chain, Z-dephasing, all γ.
**Source:** [Primordial Qubit Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md)

---

## Cockpit and Diagnostics (replace full tomography)

### 51. Decoherence cockpit: 3-observable reduction (Tier 2, verified N=2-5, IBM-validated)

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

### 52. Thermal oscillation resilience (Tier 2, verified N=4)

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

### 55. Universal absorption dose K_death (Tier 1, proven from D6)

    K_death = ln(10) = 2.303    (dose for 99% absorption)
    K_death / K_fold ~ 2.3      (ratio to CΨ = ¼ crossing dose)
    Immortal modes = N + 1      (zero absorption rate, all N)

**Derivation:** 99% absorption of the slowest mortal mode means
e^{−rate_min · t} = 0.01, so rate_min · t = ln(100). By formula D6:
rate_min = 2γ (spectral gap). Therefore t = ln(100)/(2γ), and
K = γ · t = ln(100)/2 = ln(10) = 2.303. Independent of N, γ,
topology. QED.

**Note:** The source document (Trapped Light Localization) writes
"K_death = ln(100) = 2.303," which is a typo. The value 2.303 is
ln(10), not ln(100) = 4.605. The derivation above is clean: the
factor 2 in rate_min = 2γ halves the exponent.

N+1 modes have exactly zero absorption rate (pure {I,Z} content,
invisible to the light). Complete absorption is impossible while the
palindrome holds. The cavity always retains light.

**Valid for:** any Heisenberg chain, Z-dephasing, all N.
**Replaces:** time evolution to find "when does the system die."
**Source:** [Trapped Light Localization](../experiments/TRAPPED_LIGHT_LOCALIZATION.md)

---

*Each formula in this document is a Liouvillian that does not need
to be built.*
