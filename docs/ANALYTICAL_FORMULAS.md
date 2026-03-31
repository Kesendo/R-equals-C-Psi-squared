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

*Each formula in this document is a Liouvillian that does not need
to be built.*
