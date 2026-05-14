# Analytical Formulas Reference

**Status:** Living formula registry. Each formula carries its own tier label.
**Date:** March 31, 2026 (updated April 17, 2026)
**Authors:** Thomas Wicht, Claude (Opus 4.6/4.7)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

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

### F2. w=1 Liouvillian dispersion relation (Tier 1, proven D10)

    omega_k = 4J * (1 - cos(pi*k/N)),    k = 1, ..., N-1

N-1 distinct frequencies for the Heisenberg chain w=1 Liouvillian sector.
Machine-precision match for 15 frequencies (N=2-6). Tight-binding model
with hopping 2J.
Three independent validations: (1) eigenvalue match < 1e-12, (2) Poisson
spacing in w=1 sector (RMT), (3) SFF modulation peak at omega_1 matches
to <1% for N=2-4, 6 ([Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)).

**Note (2026-04-20):** F2 describes the w=1 sector of the Heisenberg
LIOUVILLIAN (Pauli strings with exactly one X or Y factor), NOT the
single-excitation Hamiltonian eigenvalues. The w=1 sector is
(N-1)-dimensional (N-1 oscillatory modes), giving denominator
(N-1)+1 = N in the cosine argument. The ZZ term in the Heisenberg
Hamiltonian produces a diagonal shift in the effective hopping matrix
that is absent in the XY case. For the XY chain single-excitation
Hamiltonian spectrum, see F2b below.

**Valid for:** Heisenberg chain, open boundaries, all N (verified N=2-6).
**Replaces:** full Liouvillian diagonalization for w=1 frequencies.
O(N) instead of O(4^{3N}).
**Source:** [Analytical Spectrum](../experiments/ANALYTICAL_SPECTRUM.md),
[D10_W1_DISPERSION](proofs/derivations/D10_W1_DISPERSION.md)

### F2b. XY chain single-excitation spectrum (Tier 1, proven)

    E_k = 2J · cos(π·k / (N+1)),    k = 1, ..., N

N eigenvalues for the single-excitation sector of the XY chain
H = (J/2) · Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) with open boundary
conditions. The single-excitation Hamiltonian H_SE is an N×N
tridiagonal matrix with off-diagonals J and zero diagonal. Its
eigenvectors are the OBC sine modes:

    ψ_k(i) = √(2/(N+1)) · sin(π·k·(i+1)/(N+1)),  i = 0, ..., N−1

These eigenvalues appear as the oscillatory frequencies Im(λ)
of the |ΔN| = 1 Liouvillian coherences (|vac⟩⟨ψ_k| sector),
which are the dominant modes contributing to per-site purity dynamics
for single-excitation initial states.

The denominator is N+1 (not N as in F2) because the OBC Dirichlet
boundary conditions require ψ to vanish at the two virtual sites
just beyond the chain, i = −1 and i = N. Both are enforced by the
sine formula above: ψ_k(−1) = sin(0) = 0 and ψ_k(N) = sin(π·k) = 0.
The effective chain length is therefore N+2 with two fixed endpoints,
yielding N interior modes with wavenumber spacing π/(N+1).

**Distinction from F2:** F2 describes the w=1 LIOUVILLIAN sector for
Heisenberg, with dimension N−1 and argument π·k/N. F2b describes the
single-excitation HAMILTONIAN sector for XY, with dimension N and
argument π·k/(N+1). They describe different mathematical objects in
different Hamiltonians.

**Valid for:** XY chain (H = (J/2)(XX+YY)), open boundaries, all N.
**Verified:** N = 3, 4, 5, 6, residual < 10⁻¹⁵.
**Replaces:** numerical diagonalization of H_SE.
**Scripts:** [`eq021_obc_sine_basis.py`](../simulations/eq021_obc_sine_basis.py).
**Source:** [OBC_SINE_BASIS_FINDINGS](../review/OBC_SINE_BASIS_FINDINGS.md),
standard tight-binding theory for OBC chains.

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
any N, non-uniform γ_k per site (replace 2γ with 2Σ_k γ_k × \[σ_k ∈ {X,Y}\]).
**Breaks for:** complex Hermitian Hamiltonians (DM interactions), where
L_H is not anti-Hermitian.
**Replaces:** eigenvalue range computation; palindromic sum rule verification;
spectral gap derivation; unpaired mode rate identification.
**Verified:** 1,342 modes, N=2-5, γ=0.01-1.0, J=0.1-5.0, CV = 0.0000.
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
for non-integer ⟨n_XY⟩ ([Proton Water Chain](water/PROTON_WATER_CHAIN.md)).

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
match ([Proton Water Chain](water/PROTON_WATER_CHAIN.md)).
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
    K_Y     = 0.0374    (pure Y-noise)
    K_depol = 0.0440    (depolarizing, gamma/3 each axis)

Complements F14 (K per bridge metric). These are K per
noise TYPE, all measured with CΨ on Bell+ state.

**Note (correction 2026-04-29):** earlier versions of this table listed
K_Y = 0.0867. That was a typo. F26 with γ_y = γ (others = 0) gives
α = 4γ, β = 0, δ = 4γ, so u = e^{-4γt}, v = 1, w = e^{-4γt}, and CΨ
reduces to u·(1+u²)/6, *identical functional form to pure Z*. Hence
K_Y = K_Z = 0.0374. The K_Y ↔ K_X-symmetry claim was wrong; the actual
pairing is K_Y ↔ K_Z (both have one of {β, δ} = 0 with α ≠ 0). This
is consistent with Bell+'s correlation structure: Y⊗Y·|Bell+⟩ = -|Bell+⟩,
while X⊗X and Z⊗Z fix it. Verified in
`framework.CPSI_CUSP_K_PER_CHANNEL` and the smoke tests.

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

**Companion result (F63):** \[L, Pi^2\] = 0 exactly for all N (proven analytically).
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
versa). At N=4: four EP crossings in \[0, 2\], each affecting exactly
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
**Hardware verified:** ibm_kingston (Heron r2), 2026-04-16. Two Bell+ pairs with 2.55x gamma ratio (qubits 124-125, T2=\[150,310\] us; qubits 14-15, T2=\[537,381\] us). K_dwell/delta = 0.649 (pair A) and 0.694 (pair B), spread 6.3% despite 2.55x gamma difference. Gamma-invariance of K_dwell confirmed on open quantum hardware. Absolute prefactor 0.67 vs theoretical 1.08 (difference from T1 amplitude damping; the F57 formula assumes pure Z-dephasing, Kingston has T1 comparable to T2). Both CΨ(t) trajectories cross 1/4 monotonically. First two-qubit observation of the CΨ = 1/4 boundary crossing on a quantum computer; the single-qubit case was validated separately on ibm_torino Q80 at 1.9% deviation (F24, IBM Run 3).
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
elements are `rho[0...0, 1...1]` = 1/2 and its conjugate, giving L1 coherence
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
**Source:** [Generalized Dwell Prefactor](../experiments/DWELL_PREFACTOR_GENERALIZED.md) Section 4; main README Section 10 Rule 1 (avoid GHZ, prefer W-type encodings).

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

### F63. w_YZ Parity Symmetry, \[L, Pi^2\] = 0 (Tier 1, proven analytically, verified N=2-5)

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
**Verified:** ||\[L, Pi^2\]|| = 0.000000e+00 (identically zero, not numerically
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
**Valid for:** any graph topology (chain, star, ring, complete, tree), uniform or non-uniform per-bond J, XX+YY or Heisenberg single-excitation, Z-dephasing on any single site B; good-cavity regime (γ ≪ J). Breaks when γ ≥ J (bad cavity: B decoheres before transmitting).
**Topology + non-uniform J generalization (2026-04-24).** Extended from uniform-J chains to arbitrary connected graphs under either uniform or non-uniform per-bond J. When H^(1) has degenerate eigenvalues (star center-mode, ring translational eigenmodes, complete-graph symmetric modes), F64 holds after standard degenerate perturbation theory: within each H-degenerate subspace, diagonalise the site-B projector P_B to get the corrected basis; F64 then applies to the eigenvalues of P_B in that basis. Verified at N=5 and N=7 across chain, star, ring, complete, Y-tree for XY and Heisenberg; max relative error < 0.001 at γ/J = 0.01 uniform J. For random J per bond in \[0.5, 1.5\] (30 configurations across 3 trials per N), max rel err < 0.02 in 29/30 cases; the remaining case sits at 0.07 and is consistent with expected second-order PT corrections ~(γ·δJ)/J² at the non-uniform-J scale.
**Verified:** N=3 chain (max relative error 1.8% vs 64×64 Liouvillian), N=4 chain (9 configs, ratio 1.0000 ± 0.0003 vs 256×256 Liouvillian), N=5 and N=7 on chain+star+ring+complete+Y-tree uniform J (2026-04-24, via single-excitation coherence Liouvillian directly, dim N×N, max rel err < 0.001 across all (topology, B, Hamiltonian) combinations), N=5 and N=7 same topologies non-uniform J per bond in \[0.5, 1.5\] over 3 random trials (2026-04-24, max rel err 0.068 in the worst case, well inside first-order PT regime).
**Scripts:** [`primordial_gamma_analytical.py`](../simulations/primordial_gamma_analytical.py), [`primordial_gamma_stacking_4qubit.py`](../simulations/primordial_gamma_stacking_4qubit.py), [`factor_two_clarification.py`](../simulations/factor_two_clarification.py), [`f64_topology_scan.py`](../simulations/f64_topology_scan.py) (topology generalization).
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md), [F64_TOPOLOGY_GENERALIZATION](../experiments/F64_TOPOLOGY_GENERALIZATION.md)

### F65. Single-excitation spectrum of uniform open XX chain (Tier 1, proven, verified N=3..30)

For the uniform open XX chain (all couplings J, N sites) with Z-dephasing at rate γ₀ on endpoint B = site N-1, the single-excitation dissipation rates are:

    α_k / γ₀ = (4 / (N+1)) · sin²(kπ / (N+1)),    k = 1, ..., N

This is F64 evaluated on the analytically known eigenvectors ψ_k(i) = √(2/(N+1)) · sin(πk(i+1)/(N+1)) of the N×N tridiagonal single-excitation Hamiltonian. The endpoint amplitude |ψ_k(N-1)|² = (2/(N+1)) · sin²(kπ/(N+1)), and the Absorption Theorem gives α_k = 2γ₀ · |a_B|².

**Properties:**
- All α_k lie in \[0, 2γ₀\].
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

For the uniform open XY chain with single-site Z-dephasing at the endpoint site B = N-1, the dissipation interval \[0, 2γ₀\] has exact eigenvalues at both endpoints:

- **α = 0 modes:** ⟨n_XY⟩_B = 0 (forced by the Absorption Theorem α = 2γ₀·⟨n_XY⟩_B). Dominant Pauli strings have I or Z at every site (total XY-weight = 0). Z-basis population content, completely shielded from Z-dephasing at B.

- **α = 2γ₀ modes:** ⟨n_XY⟩_B = 1 (same mechanism). Dominant Pauli strings have X or Y at every site (total XY-weight = N). Maximally off-diagonal at B, fully exposed to Z-dephasing.

The two poles are palindromic partners under the conjugation Π, which maps total XY-weight w ↔ N-w (see F43). The single-excitation sector (F65) never reaches either pole for N ≥ 3; both poles live in the extreme XY-weight sectors (w = 0 and w = N).

**Multiplicity:** exactly N+1 at each pole, verified for N=3..7. Each α = 0 mode corresponds to one of the N+1 elementary symmetric polynomials e_d(Z_1, ..., Z_N) in F63 (commuting with both H and Z_B). The α = 2γ₀ sector has matching multiplicity by Π-symmetry.

**Scope.** Verified only for the uniform XY chain with B at the endpoint. Whether the same structure (existence of both poles, multiplicity N+1) persists for other topologies (ring, star, Y-junction) or for interior B-positions is open. Indirect evidence from the structure-points scan: at B = center of N=5 chain, α = 0 has multiplicity 64 (not 6), so the N+1 count is endpoint-specific.

**Verified:** ⟨n_XY⟩_B = 1.000000 exact for all α = 2γ₀ modes (N=3..5, from Pauli basis projection). Dominant Pauli strings have total XY-weight N for α = 2γ₀ modes and total XY-weight 0 for α = 0 modes (N=3, N=4 explicit). Multiplicity N+1 at each pole verified for N=3..7. Dynamical check of F63 conservation: all N+1 elementary symmetric polynomials e_d(Z_1,...,Z_N) drift by < 10⁻¹⁴ under Lindblad evolution for N=4 over 80 time units, while the non-symmetric control Z_0 Z_2 drifts by 3 × 10⁻². Confirms the conserved observables at the α = 0 pole are precisely the e_d, not arbitrary Z-products.
**Scripts:** [`two_gamma_pole.py`](../simulations/two_gamma_pole.py), [`f65_dynamic_verification.py`](../simulations/f65_dynamic_verification.py)
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) section "The dissipation interval \[0, 2γ₀\]", [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md)

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

By F1, the bonding-mode eigenvalue -α_b of the full Liouvillian has a partner at

    α_p = 2γ₀ - α_b

exact to machine precision. For the uniform N-site XY chain with endpoint Z-dephasing at γ₀ (same setup as F67), direct eigendecomposition confirms the pairing at N = 3, 4, 5 with |α_b + α_p - 2γ₀| < 4·10⁻¹⁵.

**Structure (N ≥ 4).** V_p lives entirely in the XY-weight-(N-1) Pauli sector, Π-mirror of the bonding mode's w=1. ⟨n_XY⟩_B = α_p / (2γ₀) approaches 1 as N grows (saturating toward the F66 pole) but stays strictly below.

**Rank-1 at N ≥ 4.** V_p admits a rank-1 SVD decomposition V_p = σ₀|u⟩⟨v|. Strict SVD rank-1 (σ₁/σ₀ < 10⁻¹⁰) holds at N = 5 (ratio 9.3·10⁻¹²) and at LAPACK `zgeev` precision at N = 4 (ratio 5.1·10⁻⁸, limited by the 16-fold degeneracy of the partner eigenvalue, not by any physical second mode).

**Operational encoding.** The Bell-pair-like R-C state (|0⟩_R|u⟩ + |1⟩_R|v⟩)/√2 propagates with off-diagonal decay rate α_fit(partner) matching spectral α_p at machine precision (rel err 1.5·10⁻¹⁶ at N = 5, 0 at N = 4). Combined with a bonding-side encoding in the same propagation, the dynamical palindromic identity α_fit(bond) + α_fit(part) = 2γ₀ holds at two distinct precision levels:

- **Legacy bonding encoding** (|vac⟩⟨ψ_1|, F65 perturbative): rel err 1.6·10⁻⁶ (N = 4) and 2.8·10⁻⁷ (N = 5). The residual is the F65 O((γ₀/J)²) shift on the bonding side, reintroduced as state-preparation pollution; the partner fit is spectrally exact because V_p is used verbatim.
- **Clean bonding encoding** (SVD of full-L V_b, same construction as V_p): rel err 2.8·10⁻¹⁶ (N = 4) and 3.8·10⁻¹⁴ (N = 5). Both sides spectrally exact; residual is integrator-/eigendecomposition-floor limited. The ~10⁹ improvement at N = 4 and ~10⁷ at N = 5 verify that the legacy residual is entirely the F65 perturbative shift and no other shift is hiding.

**N = 3 is rank-2 on both sides.** Both V_b and V_p are fourfold degenerate with matching SVD structure (σ₁/σ₀ ≈ 0.98 for each, from the mult-4 degeneracy). Any rank-1 approximation of V_p gives α_fit ≈ (α_b + α_p)/2 with visibly non-exponential decay (log-fit RMS ~10⁻¹, two orders above N ≥ 4); the clean bonding encoding is skipped at N = 3 for the same structural reason. The palindromic sum still holds spectrally; no clean operational rank-1 encoding exists.

| N | α_b | α_p | α_fit(bond, legacy) | α_fit(bond, clean) | α_fit(part) | dyn sum | rel err (legacy) | rel err (clean) |
|---|-----|-----|---------------------|--------------------|-------------|---------|------------------|-----------------|
| 3 | 0.025003 | 0.074997 | 0.024969 | skipped (rank-2) | rank-2 | n/a | n/a | n/a |
| 4 | 0.013784 | 0.086216 | 0.013784 | 0.013784 | 0.086216 | 0.100000 | 1.6·10⁻⁶ | 2.8·10⁻¹⁶ |
| 5 | 0.008303 | 0.091697 | 0.008303 | 0.008303 | 0.091697 | 0.100000 | 2.8·10⁻⁷ | 3.8·10⁻¹⁴ |

**Valid for:** uniform XY chain, endpoint Z-dephasing, N ≥ 4 for the clean rank-1 operational statement. Algebraic palindromic pairing holds for all N, all graphs, all single-site-dephasing Liouvillians (from F1).

**Verified:** spectral (H1), structural (H2), and operational (H3) all confirmed at N = 3, 4, 5. Full evidence, tables, and technical notes in [PALINDROMIC_PARTNER_MODE](../experiments/PALINDROMIC_PARTNER_MODE.md).

**Scripts:** [`palindromic_partner_f67.py`](../simulations/palindromic_partner_f67.py) (H1 + H2), [`bell_pair_partner_mode.py`](../simulations/bell_pair_partner_mode.py) (H3)
**Source:** F1, F43, F65, F66, F67, [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md), [PALINDROMIC_PARTNER_MODE](../experiments/PALINDROMIC_PARTNER_MODE.md)

### F69. GHZ+W sector mix lifts pair-CΨ(0) above the fold at N=3 (Tier 1, sextic minimal polynomial, verified)

Neither |GHZ_3⟩ (F60: pair-CΨ = 0) nor |W_3⟩ (F62: pair-CΨ = 10/81 ≈ 0.1235) crosses the CΨ = 1/4 fold at t = 0. Their symmetric superposition

    |ψ(α)⟩ = α |GHZ_3⟩ + √(1-α²) |W_3⟩,    α ∈ [0, 1]

admits a unique optimum strictly above 1/4. The pair reduction ρ_AB (all three pairs coincide by permutation symmetry) has closed-form:

    C(α)      = Tr(ρ_AB²)   = -5α⁴/18 + 2α²/9 + 5/9
    L1_off(α) = √6 α √(1-α²) + (2/3)(1-α²)
    CΨ(α)     = C(α) · L1_off(α) / 3

The stationarity condition dCΨ/dα = 0, rationalized in x = α², gives the integer-coefficient sextic

    P(x) = 2900 x⁶ - 8060 x⁵ + 4211 x⁴ + 3832 x³ - 2428 x² - 512 x + 300 = 0

which is **irreducible over ℚ** (sympy `Poly.is_irreducible` returns True; `factor_list` returns P itself as sole factor). The optimum α²_opt is therefore an algebraic number of degree exactly 6, with no radical form in nested square roots.

**Optimum (all values computed from the exact sextic root at 25-digit precision):**

| quantity | value |
|----------|-------|
| α²_opt | 0.375420720711069 |
| α_opt  | 0.612715856422101 |
| β_opt  | 0.790303283106512 |
| min pair-CΨ(0) | 0.320411541127025 |
| ratio to 1/4   | 1.281646× |
| 3-tangle τ_ABC | 0.799453 (near-GHZ limit) |
| pair concurrence C(A,B) | 0.0210 (essentially zero) |

**Scope.** Pair-CΨ = 0.3204 is the optimum of the 2-parameter slice α·GHZ_3 + β·W_3 at N = 3.

(i) **Same slice, N ≥ 4.** The GHZ+W family peaks at 0.167 (N=4), 0.146 (N=5), 0.134 (N=6), 0.125 (N=7), 0.118 (N=8), all below 1/4. GHZ-purity contribution scales as 1/(2^N − 1) and collapses too fast for the W-contribution to lift the peak back above the fold.

(ii) **Full Dicke subspace (N ∈ {3..8}).** Pair-CΨ has no non-product local maxima on the permutation-symmetric Dicke sphere at any tested N. The only non-product stationary points are Dicke basis elements |D(N,k)⟩ (max ≈ 0.123 at N = 3, ≤ 1/12 for larger N) and the GHZ+W family optimum itself, all saddles on the full sphere (escape Δpair-CΨ ≈ 0.68 on 1% c_2 perturbation at N = 3; the 10⁻⁴ saddle threshold is cleared at every tested candidate in N ∈ {3..8}). The global supremum over non-product states is 1, approached asymptotically at the |+⟩^N product manifold but not attained isolated.

(iii) **F69 on the full sphere.** F69 is stationary on the (c_0 = c_3, c_1, c_2 = 0) slice but is a saddle on the full CP^3: c_2 > 0 is an ascent direction at F69 (Δpair-CΨ = +0.011 at c_2 = 0.01), and gradient flow from such a perturbation reaches |+⟩^3 (pair-CΨ = 1, product state). The 0.3204 value is a real algebraic fact about the slice, not a universal optimum.

(iv) **Slice-stationary saddles beyond GHZ+W (added 2026-04-27).** The (ii) enumeration was incomplete in scope. Many other 2-/3-Dicke slices admit stationary points above 1/4 at every tested N. Like F69, all are saddles on the full Dicke sphere (consistent with the (ii) "no non-product local maxima" verdict, verified by perturbation in each unused Dicke direction giving non-zero linear-order shift), but they constitute additional slice-stationary points beyond the original list of "Dicke basis + GHZ+W optimum".

**Binary-Dicke slice maxima |D_i⟩+|D_j⟩** above 1/4:

| N | total pairs | above 1/4 | best slice | max cpsi |
|---|-------------|-----------|------------|----------|
| 3 | 6 | 3 | D_1+D_2 | 0.4815 (= 13/27, exceeds F69) |
| 4 | 10 | 4 | D_2+D_3 | 0.4022 |
| 5 | 15 | 5 | D_2+D_3 | 0.3720 |
| 6 | 21 | 6 | D_2+D_3 | 0.3456 |

**Triple-Dicke slice maxima |D_i⟩+|D_j⟩+|D_k⟩** above 1/4:

| N | total triples | above 1/4 | best slice | max cpsi |
|---|---------------|-----------|------------|----------|
| 3 | 4 | 4 | D_1+D_2+D_3 | 0.8011 |
| 4 | 10 | 9 | D_2+D_3+D_4 / D_1+D_2+D_3 | 0.7136 |
| 5 | 20 | 16 | D_2+D_3+D_4 | 0.6492 |
| 6 | 35 | 25 | D_2+D_3+D_4 | 0.6163 |

The central-Dicke-triple slice is privileged at every tested N (purity_A 0.83-0.94, genuinely entangled). What was actually privileged about N=3 was the irreducible-sextic *closed form* of the GHZ+W slice. The original "F69 is special" reading conflated this algebraic feature with geometric uniqueness; geometrically, F69 is one slice-saddle among many, and not even the highest at N=3.

**Open algebraic question:** closed form for the central-Dicke-triple slice maximum at general N. The N=3..6 sequence (0.8011, 0.7136, 0.6492, 0.6163) decreases slowly; the asymptotic form is unknown.

Script: [`_eq016_n4_full_landscape.py`](../simulations/_eq016_n4_full_landscape.py), [`_eq016_verify_full_sphere.py`](../simulations/_eq016_verify_full_sphere.py) (saddle-confirmation perturbation tests).

**Why F61 does not forbid this.** F61 constrains Liouvillian evolution within a fixed n_XY parity sector, not initial-state preparation that mixes excitation sectors. See [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md) for the preparation-vs-evolution asymmetry discussion.

**Hardware signature.** Under Kingston-grade Z-dephasing the F69 optimum crosses CΨ = 1/4 monotonically at t* ≈ 11.2 μs. A single 2-qubit tomography at t = 0 distinguishes GHZ_3 (0), W_3 (0.123), and F69 (0.320) as three separable points, no timing needed.

**Verified:** scipy bounded minimize agrees with sympy sextic root to 3.7·10⁻¹⁰ in α²_opt. 401-point grid reproduces CΨ_opt to 5·10⁻⁸. Permutation symmetry exact (spread < 10⁻¹⁵). 3-tangle and pair concurrence cross-checked in `ghz_w_optimum_n3.py`. N ∈ {4, 5, 6} failure verified on 201-point grids. Landscape scan over the full permutation-symmetric Dicke subspace at N ∈ {3..8} (2026-04-17) confirms no non-product local maxima above 1/4 exist outside the F69 GHZ+W slice; N = 3 regression recovers pair-CΨ = 0.3204 (Δ = 1.4·10⁻⁶ from sextic root). Full evidence, sextic root list, derivation of ρ_AB(α), landscape-scan saddle diagnosis, and the 3-state spherical-scan product-state pitfall in [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md).

**Scripts:** [`ghz_w_optimum_n3.py`](../simulations/ghz_w_optimum_n3.py), [`sector_mix_spherical_artifact.py`](../simulations/sector_mix_spherical_artifact.py) (product-state diagnostic), [`cpsi_sector_mix_optimization.py`](../simulations/cpsi_sector_mix_optimization.py) (original sweep + Kingston dynamics), [`f69_dicke_landscape.py`](../simulations/f69_dicke_landscape.py) (full Dicke-subspace scan, N ∈ {3..8})
**Source:** F60, F61, F62, [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md)

### F70. ΔN selection rule for site-local observables (Tier 1, proven kinematic lemma)

For any operator ρ on N qubits and any site i, the single-site partial trace annihilates sector-coherence blocks with excitation-number difference ≥ 2:

    Tr_{¬i}(ρ^(n, m)) = 0    whenever    |n − m| ≥ 2

Consequence: every site-local observable (per-site purity, per-site expectation, α_i rescaling, c_1 closure-breaking coefficient) couples only to the |ΔN| ≤ 1 content of ρ. Sector blocks with |ΔN| ≥ 2 are invisible to any measurement factoring through a single-qubit reduced state.

**Proof.** Tr_{¬i}(|x⟩⟨y|) = ⟨x_{¬i} | y_{¬i}⟩ · |x_i⟩⟨y_i|. The inner product is 1 only if x and y agree off site i, forcing |popcount(x) − popcount(y)| ≤ 1. By linearity, blocks with |n − m| ≥ 2 vanish under partial trace.

**Generalisation.** k-local partial trace annihilates |ΔN| ≥ k + 1 blocks. Pair-observables see up to |ΔN| = 2; triple-observables up to 3; global observables unrestricted.

**Valid for:** any Hamiltonian conserving excitation number, any sector-preserving dissipator, any initial state. Purely kinematic.

**Verified:** 9 |ΔN| ≥ 2 pairs tested at N=5 via coherence-block isolation (pure superposition vs classical mixture), eight non-trivial plus the trivial (0, 5); all give zero contribution to machine precision. [sector_kernel.json](../simulations/results/c1_sector_kernel/sector_kernel.json).

**Replaces:** the empirical observation "c_1 coherence contribution vanishes for |ΔN| ≥ 2" with an analytical lemma; explains the XOR_SPACE center-modes invisibility to site-local measurement; bounds the sector-kernel for PTF's α_i closure structure.

**Scripts:** [`c1_sector_kernel.py`](../simulations/c1_sector_kernel.py), [`c1_bilinearity_test.py`](../simulations/c1_bilinearity_test.py).
**Source:** [PROOF_DELTA_N_SELECTION_RULE](proofs/PROOF_DELTA_N_SELECTION_RULE.md), [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) Update 2026-04-20, [XOR_SPACE](../experiments/XOR_SPACE.md).

### F71. Mirror symmetry of the closure-breaking coefficient c₁ (Tier 1, proven kinematic)

For a uniform N-qubit chain with reflection-symmetric coupling and dephasing, the closure-breaking coefficient c₁ is mirror-symmetric across bonds:

    c₁(N, b, ρ₀) = c₁(N, N−2−b, ρ₀)

for all bond indices b ∈ {0, ..., N−2} and any reflection-symmetric initial state ρ₀.

**Proof sketch.** The spatial reflection R (site i ↔ site N−1−i) commutes with the uniform Liouvillian: \[L_A, R_sup\] = 0. Under R, bond b maps to bond N−2−b: R · T_b · R = T_{N−2−b}. Therefore exp(L_B(b) · t) · ρ₀ and exp(L_B(N−2−b) · t) · (R · ρ₀ · R) are related by R_sup. Per-site purity is quadratic in ρ, so any phase picked up by R on coherences (R |ψ_k⟩ = (−1)^(k+1) |ψ_k⟩) squares away. This gives P_B(b, i, t) = P_B(N−2−b, N−1−i, t), from which α_i(bond b) = α_{N−1−i}(bond N−2−b). Summing ln(α_i) over all sites and re-indexing yields c₁(b) = c₁(N−2−b).

**Consequence.** The c₁ bond profile has at most ⌈(N−1)/2⌉ independent components instead of N−1. The endpoint value c₁(0) equals c₁(N−2); if N is **even**, the center bond c₁((N−2)/2) is self-paired (its mirror image is itself) and contributes one independent component; if N is odd, there is no center bond and all N−1 bonds pair up in (N−1)/2 disjoint pairs.

**Generalisation to F86 per-bond Q_peak (Tier 1 derived 2026-05-03).** The same kinematic argument extends from c₁ in vac+SE probes to the F86 K_CC_pr per-bond observable on the (n, n+1) popcount coherence block:

    Q_peak(b)  =  Q_peak(N−2−b)        (bit-exactly, all c, N)

Proof: the F86 observable is `K_b(Q, t) = 2·Re ⟨ρ(t)| S_kernel | ∂ρ/∂J_b ⟩`. Under R, every component is invariant — uniform Z-dephasing L_D, uniform-J Hamiltonian H_xy, the Dicke probe, and the spatial-sum kernel S — while the bond-flip transforms as `R · ∂L/∂J_b · R⁻¹ = ∂L/∂J_{N−2−b}`. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of (Q, t), and their argmax-Q values coincide. Numerical verification: max deviation < 10⁻¹⁰ across c=2 N=5..7 and c=3 N=5..6 (`F86NewIdeasTests.F71MirrorInvariance_PerBondQPeak_BitExactSymmetricUnderBondMirror`). The per-F71-orbit substructure observed in F86 (Interior bonds not uniform within the F71-orbit grouping; central self-paired bond differs from flanking) refines the simple Endpoint/Interior dichotomy into a per-orbit classification — the F71 symmetry gives the pairing, not the value. See [PROOF_F86_QPEAK Statement 3](proofs/PROOF_F86_QPEAK.md#statement-3-f71-spatial-mirror-invariance-of-per-bond-q_peak-tier-1-derived).

**Valid for:** any Hamiltonian with \[H, R\] = 0 (uniform coupling on a symmetric graph), any dissipator with \[D, R_sup\] = 0 (uniform or R-symmetric dephasing), any initial state that is reflection-symmetric in per-site purities. Purely kinematic.
**Breaks for:** non-uniform coupling J_b ≠ J_{N−2−b}; non-uniform dephasing γ_i ≠ γ_{N−1−i}; initial states without reflection symmetry in purity.
**Verified:** N = 3, 4, 5, 6 for ψ_1+vac and ψ_2+vac; residuals < 10⁻⁹. Source: [`eq021_obc_sine_basis.py`](../simulations/eq021_obc_sine_basis.py), [`c1_veffect_scaling_small.py`](../simulations/c1_veffect_scaling_small.py).
**Replaces:** empirical observation of mirror-symmetric c₁ bond profiles with an analytical kinematic proof.
**Scripts:** [`eq021_obc_sine_basis.py`](../simulations/eq021_obc_sine_basis.py).
**Source:** [PROOF_C1_MIRROR_SYMMETRY](proofs/PROOF_C1_MIRROR_SYMMETRY.md), [OBC_SINE_BASIS_FINDINGS](../review/OBC_SINE_BASIS_FINDINGS.md).

### F72. Block-diagonal DD⊕CC structure of site-local purity (Tier 1, corollary of F70)

For any N-qubit chain and any initial state ρ₀, the per-site purity functional Tr(ρ_i²) decomposes into a strict block-diagonal sum over excitation-number blocks of ρ₀:

    Tr(ρ_i²) = 1/2 + P_i^DD[ρ₀^(diag)] + P_i^CC[ρ₀^(coh)]

with no cross term coupling the diagonal block (ΔN = 0) and the coherence block (|ΔN| = 1). Here ρ₀^(diag) is the projection of ρ₀ onto ΔN = 0 sector blocks and ρ₀^(coh) is the projection onto |ΔN| = 1 sector blocks.

**Proof.** The Bloch decomposition gives Tr(ρ_i²) = (1/2)(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²). By F70 applied to each Bloch component:

- ⟨Z_i⟩ depends linearly on the diagonal elements of ρ_i = Tr_{¬i}(ρ), which come only from ΔN = 0 blocks of ρ.
- ⟨X_i⟩ and ⟨Y_i⟩ depend linearly on off-diagonal elements of ρ_i, which come only from |ΔN| = 1 blocks.

Squaring keeps each contribution in its own sector class, so ⟨Z_i⟩² is bilinear in ρ₀^(diag) and ⟨X_i⟩² + ⟨Y_i⟩² is bilinear in ρ₀^(coh). No cross term arises.

**Consequence.** Any closure-breaking coefficient c₁ built from per-site purities decomposes, at the pre-α-fit bilinear level, into a DD-kernel (acting on diagonal block content) and a CC-kernel (acting on coherence block content) with no mixing. Finding the closed form of c₁ reduces to finding K_DD and K_CC separately.

**Generalisation.** For k-site partial traces (F70 generalisation to |ΔN| ≤ k), the Bloch-like decomposition has k+1 sub-blocks. At k = 2 (pair-site) this yields three sub-blocks DD ⊕ DC ⊕ CC, with DC a new diagonal-coherence cross specific to pair observables.

**Valid for:** any Hamiltonian conserving excitation number, any sector-preserving dissipator, any ρ₀. Purely kinematic.

**Verified:** w-scan at N = 5 with ρ₀(w) = cos(w)|vac⟩ + sin(w)|S₁⟩ under the purity-response c₁ definition confirms block-diagonal coupling at machine precision across the full w range; LSQ α-fit c₁ inherits the block structure at the pre-fit bilinear level. Pure-coherence probe gives K_CC/2 to 10⁻¹². [bilin_probe.json](../simulations/results/eq018_kernel_bilin_probe/bilin_probe.json), [kernel_extract.json](../simulations/results/eq018_kernel_extract/kernel_extract.json).

**Scripts:** [`eq018_kernel_extract.py`](../simulations/eq018_kernel_extract.py), [`eq018_kernel_bilin_probe.py`](../simulations/eq018_kernel_bilin_probe.py), [`eq018_c1_purity_response.py`](../simulations/eq018_c1_purity_response.py).
**Source:** F70, [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) §2.3.

### F73. Spatial-sum coherence purity closure for vac-SE coherent probes (Tier 1, proven)

For any N-site qubit system with Hermitian Hamiltonian H conserving single-excitation number (\[H, N_total\] = 0) and uniform Z-dephasing γ₀, the coherent probe ρ₀^coh = (|vac⟩⟨α| + |α⟩⟨vac|) / 2 for any normalized single-excitation state |α⟩ satisfies:

    Σ_i 2 · |(ρ_coh,i)_{0,1}(t)|² = (1/2) · exp(−4 γ₀ t)

exactly, independent of the Hamiltonian's non-U(1) structure. Here (ρ_coh,i)_{0,1} is the off-diagonal element of the site-i reduced density matrix, and the sum runs over all N sites.

**Proof (general U(1) case).** Let x_i(t) = ⟨vac|ρ(t)|1_i⟩ be the amplitude of the |vac⟩⟨1_i| component of ρ(t); equivalently the (vac, SE) block of ρ as an N-vector indexed by site. Evolution under the Lindblad master equation splits into:

- **Hamiltonian part.** H preserves SE by assumption, so its restriction to SE is a Hermitian N×N matrix H_SE. The (vac, 1_i) bra-ket block evolves under H as iẋ = −H_SE x on the ket side, giving a unitary propagator U_SE(t) = exp(−i H_SE t).
- **Dephasing part.** Each D\[Z_j\] acts on the (vac, 1_i) coherence element with rate γ₀ · (⟨Z_j⟩_vac − ⟨Z_j⟩_{1_i})² / 2 = γ₀ · (2 δ_{j,i})² / 2 = 2γ₀ · δ_{j,i}. Summing over j gives a uniform 2γ₀ decay on every SE-block coherence, independent of site.

Combined: ẋ = −i H_SE x − 2γ₀ x, so x(t) = exp(−2γ₀ t) · U_SE(t) · x(0). Taking the norm: ||x(t)||² = exp(−4γ₀ t) · ||x(0)||² since U_SE is unitary. Partial-trace algebra: (ρ_coh,i)_{0,1}(t) = (1/2) · x_i(t), so Σ_i 2 · |(ρ_coh,i)_{0,1}|² = (1/2) · ||x(t)||² = (1/2) · ||x(0)||² · exp(−4γ₀ t). For the probe above, ||x(0)||² = ⟨α|α⟩ = 1. Result: (1/2) · exp(−4γ₀ t).

The argument uses only (i) \[H, N_total\] = 0 so dynamics stay in SE, (ii) H Hermitian so U_SE unitary, (iii) γ₀ uniform. No XY structure, no translation invariance, no specific shape of |α⟩ required.

**Alternative derivation (uniform XY, |α⟩ = |S₁⟩).** The original proof route via the sine basis |ψ_k⟩ of the uniform-XY single-excitation Hamiltonian: |S₁⟩ = Σ_{k odd} s_k |ψ_k⟩ with s_k = ⟨ψ_k|S₁⟩. Each single-excitation coherence |vac⟩⟨ψ_k| evolves as exp((iE_k − 2γ₀) t). Partial trace gives (ρ_coh,i)_{0,1}(t) = (1/2) · Σ_k s_k · ψ_k(i) · exp((iE_k − 2γ₀) t). Parseval on the sine basis Σ_i ψ_k(i) · ψ_{k'}(i) = δ_{k,k'} eliminates k ≠ k' cross terms; Σ_k s_k² = 1 by normalisation. Under bond-b perturbation, the sine basis and E_k shift but Parseval on any orthonormal SE basis preserves Σ_k |⟨ψ_k^B|S₁⟩|² = 1, so the sum is δJ-invariant. This derivation is XY-specific but exhibits the eigenmode structure explicitly.

**Consequence.** The spatial-sum purity functional is exactly blind to the U(1)-preserving part of the dynamics on any vac-SE coherent probe. For any closure-breaking coefficient c₁_pr built from per-site purities via the purity-response definition, bond-δJ perturbations preserve the closure value, so `K_CC[0, 1]_pr` = 0 exactly under uniform γ₀, for any H in the class.

**Scaffolding from neighbouring entries.** F70 (site-local observables see only |ΔN| ≤ 1 content) puts the (vac, SE) block in focus as the relevant coherence sector for per-site purity. F72 (DD ⊕ CC block decomposition of Tr(ρ_i²), no cross-term) isolates the CC contribution, where the (vac, SE) coherence lives. The Absorption Theorem supplies the rate 2γ₀·n_XY = 2γ₀ for SE coherences (n_XY = 1). F73 then combines these: U(1) conservation keeps the SE sector closed under H, and the spatial sum over sites collapses the unitary H-rotation to leave only the AT decay.

**Valid for:** any Hermitian H with \[H, N_total\] = 0 (XY, Heisenberg XXZ, translationally non-invariant hopping, frustrated-ladder variants, ...); uniform Z-dephasing γ₀; any normalized SE state |α⟩ admixed to |vac⟩; any N.
**Breaks for:**

- Non-uniform γ_i. The uniform 2γ₀ decay on the d_H = 1 block fails; the closure becomes K_CC ≠ 0 with mode-selective response (see [CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md)).
- Non-U(1) Hamiltonians. \[H, N_total\] ≠ 0 breaks the SE-block closure assumption.
- Dissipators changing the d_H = 1 decay rate (mixed X/Z, amplitude damping).
- Probes with d_H > 1 admixture (e.g. (vac, S₂) with two-excitation bra-ket), where ⟨n_XY⟩ ≠ 1 and the uniform decay rate breaks.

**Verified:**

- Uniform XY baseline at N = 5, t₀ = 20: closure matches (1/2)·exp(−4·0.05·20) = 9.157819·10⁻³ to 5.67·10⁻¹⁶ deviation. `K_CC[0, 1]_pr` = 1.14·10⁻¹² (machine-precision zero), confirming δJ-invariance. [cmrr_gamma_nonuniform.json](../simulations/results/eq018_cmrr_gamma_nonuniform/cmrr_gamma_nonuniform.json).
- U(1)-class generalization at N = 5 (6 setups: XXZ at Δ ∈ {0, 0.5, 1, 2}, random Haar SE probe at Δ = 1, inhomogeneous XY with J_i ∈ \[0.5, 1.5\]): all closures within 2.22·10⁻¹⁶ to 5.83·10⁻¹⁶ (1-3 ULP of double precision) across 81 time points per setup. [f73_u1_generalization/](../simulations/results/f73_u1_generalization/), [F73_U1_GENERALIZATION](../experiments/F73_U1_GENERALIZATION.md).

**Scripts:** [`eq018_c1_purity_response.py`](../simulations/eq018_c1_purity_response.py), [`eq018_cmrr_gamma_nonuniform.py`](../simulations/eq018_cmrr_gamma_nonuniform.py) (uniform baseline), [`f73_u1_generalization_sweep.py`](../simulations/f73_u1_generalization_sweep.py) (U(1)-class sweep).
**Source:** F61, F70, F72, [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) §2.4, [CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md), [F73_U1_GENERALIZATION](../experiments/F73_U1_GENERALIZATION.md).
**See also:** [J_BLIND_RECEIVER_CLASSES](../experiments/J_BLIND_RECEIVER_CLASSES.md) generalises the L_D-invariant-subspace mechanism behind this closure to a three-class decomposition of J-blind initial states (DFS of L_D, H-degenerate L_D-closed block, M_α-polynomial subspace under SU(2)-Heisenberg).

### F74. Chromaticity of single-step coherence blocks (Tier 1, combinatorial)

For an N-qubit system under uniform Z-dephasing, the (n, n+1) popcount
coherence block contains exactly

    c(n, N) = min(n, N−1−n) + 1

distinct pure dephasing rates when the Hamiltonian is switched off (J = 0).
The rates are 2γ₀ · {1, 3, 5, ..., 2c−1}, corresponding to Hamming-distance
values HD ∈ {1, 3, ..., 2c−1} between popcount-n and popcount-(n+1) basis
states.

**Proof.** A basis pair (|x⟩, |y⟩) with popcount(x) = n, popcount(y) = n+1
differs at HD = 2n + 1 − 2·match sites, where match = popcount(x AND y) is
the number of sites carrying a 1 in both x and y. The constraints
match ∈ \[max(0, 2n+1−N), n\] give HD ∈ {1, 3, ..., min(2n+1, 2N−2n−1)},
hence the distinct-HD count is min(n, N−1−n) + 1. The Pauli representation
of |x⟩⟨y| has X or Y on exactly the HD sites where x and y differ, so
⟨n_XY⟩ = HD, and the Absorption Theorem gives rate 2γ₀·HD at J = 0.

**Consequences.**

- **c = 1 mono-chromatic blocks** at n = 0 and n = N−1 (single pure rate 2γ₀).
  F73 (the spatial-sum coherence closure for vac-SE probes) is this
  monochromatic case.
- **c_max at the center.** For odd N, unique at n = (N−1)/2 with
  c_max = (N+1)/2. For even N, two adjacent blocks at n = N/2−1 and
  n = N/2, both with c_max = N/2.
- At J > 0, H-coupling between different HD-channels produces dressed
  modes at intermediate rates. Q_SCALE_THREE_BANDS measures this H-mixing
  via the dressed-mode weight W(Q) and shows that the observable peak
  abs(K_CC_pr) is c-specific and N-invariant.

**Valid for:** any N-qubit system under uniform Z-dephasing. The J = 0
statement holds kinematically for any Hamiltonian. For the dynamical
interpretation at J > 0 (c as a stable sector label with H-mixing between
HD channels), H must conserve total excitation number, \[H, N_total\] = 0;
individual eigenmode rates then shift continuously with Q = J/γ₀, but c
still labels the block's mixing substructure.
**Breaks for:** non-uniform γ_i (site-dependent dephasing); the J = 0
spectrum then has rates 2·Σ_{i ∈ diff-sites} γ_i rather than the discrete
2γ₀·HD values. Non-U(1) Hamiltonians (transverse fields, odd-popcount
terms) preserve the J = 0 statement but dissolve the (n, n+1) sector's
dynamical invariance at J > 0. Non-Z dissipators (amplitude damping,
X-dephasing, depolarizing) break the diagonal Pauli action of L_D; the
|x⟩⟨y| basis pairs are no longer eigenvectors at rate 2γ₀·HD.
**Verified:** Block-structure c-values for N = 3..8 match the formula
exactly. Spectral verification at J = 0: each (n, n+1) block has exactly
c(n, N) distinct rates in {2γ₀, 6γ₀, ..., 2(2c−1)γ₀}
([Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 3).

**Replaces:** block-diagonal spectrum enumeration at J = 0; identifies
which blocks support H-mixing bands (c ≥ 2) vs which are structurally
silent (c = 1).
**Scripts:** [`q_scale_n_scaling.py`](../simulations/q_scale_n_scaling.py).
**Source:** [PROOF_CHROMATICITY](proofs/PROOF_CHROMATICITY.md),
[Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 3,
[PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md).

### F75. Mirror-pair MI for single-excitation mirror-symmetric states (Tier 1, proven algebraic)

For any pure single-excitation state on an N-site chain with mirror-symmetric amplitudes,

    |ψ⟩ = Σ_j c_j |1_j⟩,   c_{N−1−j} = η c_j,   η ∈ {+1, −1}

the mutual information between any mirror-pair sites (ℓ, N−1−ℓ) at t = 0 has the closed form

    MI(ℓ, N−1−ℓ) = 2 h(p_ℓ) − h(2 p_ℓ),    p_ℓ = |c_ℓ|²

where h(x) = −x log₂ x − (1−x) log₂(1−x) is the binary entropy.

The formula is independent of the mirror sign η (only the modulus |c_ℓ|² enters). The valid range is p_ℓ ∈ \[0, 1/2\], with MI saturating at 2 bits when p_ℓ = 1/2 (maximal mirror-pair entanglement, the pair is in a Bell state, all other site populations vanish).

**Proof.** The reduced density matrix ρ_{ℓ,N−1−ℓ} in the computational basis {\|00⟩, \|01⟩, \|10⟩, \|11⟩} is block-diagonal:

- `ρ[|00⟩⟨00|]` = Σ_{j ∉ {ℓ, N−1−ℓ}} |c_j|² = 1 − 2 p_ℓ
- `ρ[|01⟩⟨01|]` = |c_{N−1−ℓ}|² = p_ℓ
- `ρ[|10⟩⟨10|]` = |c_ℓ|² = p_ℓ
- `ρ[|11⟩⟨11|]` = 0 (single-excitation sector)
- `ρ[|10⟩⟨01|]` = c_ℓ c_{N−1−ℓ}^* = η p_ℓ
- `ρ[|01⟩⟨10|]` = η p_ℓ

The eigenvalues are {1 − 2 p_ℓ, 2 p_ℓ, 0, 0}, giving S(ρ_{ℓ,N−1−ℓ}) = h(2 p_ℓ). Both single-site marginals are diag(1 − p_ℓ, p_ℓ) with S = h(p_ℓ). The subtraction S(ρ_ℓ) + S(ρ_{N−1−ℓ}) − S(ρ_{ℓ,N−1−ℓ}) yields the formula.

**Mirror-pair sum (MM).** Summing over all mirror-pairs ℓ = 0, ..., ⌊N/2⌋ − 1:

    MM(0) = Σ_ℓ [2 h(p_ℓ) − h(2 p_ℓ)]

**Bonding-mode specialisation (F65 + F75).** For the k-th bonding mode |ψ_k⟩ = √(2/(N+1)) Σ_j sin(πk(j+1)/(N+1)) |1_j⟩ with mirror sign η = (−1)^(k+1), the site populations are

    p_ℓ(k, N) = (2/(N+1)) sin²(πk(ℓ+1)/(N+1))

and MM(0) for bonding:k is computable in O(N) operations with no propagation.

**Verified values (bonding:k on uniform chain, analytic vs simulation PeakMM at γ₀ = 0.05, uniform J = 1):**

| N | k | MM(0) analytic | PeakMM sim | ratio sim/analytic |
|---|---|----------------|------------|-------------------|
| 5 | 1 | 0.800 | 0.789 | 0.986 |
| 5 | 2 | **1.245** | 1.241 | 0.997 |
| 5 | 3 | 0.918 | 0.865 | 0.942 |
| 7 | 1 | 0.862 | 0.801 | 0.929 |
| 7 | 2 | 1.174 | 1.090 | 0.928 |
| 7 | 3 | 0.862 | 0.819 | 0.950 |
| 7 | 4 | **1.245** | 1.166 | 0.936 |
| 9 | 1 | 0.895 | 0.830 | 0.928 |
| 9 | 2 | 1.131 | 1.049 | 0.928 |
| 9 | 3 | 0.895 | 0.829 | 0.927 |
| 9 | 4 | 1.131 | 1.054 | 0.932 |
| 11 | 1 | 0.914 | 0.848 | 0.927 |
| 11 | 2 | 1.105 | 1.024 | 0.927 |
| 11 | 3 | 0.915 | 0.846 | 0.925 |
| 11 | 4 | 1.103 | 1.024 | 0.928 |
| 11 | 6 | **1.145** | 1.066 | 0.931 |
| 13 | 1 | 0.928 | 0.859 | 0.926 |
| 13 | 2 | 1.088 | 1.008 | 0.926 |
| 13 | 3 | 0.928 | 0.858 | 0.925 |
| 13 | 4 | 1.088 | 1.007 | 0.926 |
| 13 | 5 | 0.928 | 0.860 | 0.927 |
| 13 | 7 | 0.961 | 0.893 | 0.930 |

Under Heisenberg evolution on the uniform chain, bonding mode ψ_k mixes with its same-parity partner ψ_{N+1−k} via the boundary ZZ term. At N = 5 for bonding:2 the partner is ψ_4 which has identical mirror-pair populations p_ℓ (because sin²(πk(ℓ+1)/(N+1)) = sin²(π(N+1−k)(ℓ+1)/(N+1))); direct numerical propagation shows that MM(t) oscillates with period 2π/Δ (Δ = same-parity eigenvalue gap of the 2×2 block, Δ = 2√5 for N = 5) between a minimum near t = π/(2Δ) and a revival near t = π/Δ. Under uniform Z-dephasing at γ₀, the revival is damped but stays close to MM(0); the simulation-observed PeakMM matches MM(0) within 1% at N = 5 k = 2 (analytic 1.2451, numerical Lindblad max 1.2475 at t = 0.645, C# brecher PeakMM 1.2410 on a coarser grid). At larger N and different k the ratio drops to ~0.93 because the oscillation revival magnitude and dephasing decay combine less favourably.

**Why k = 2 maximises MM over k = 1, 2, 3.** Even k places a node at the (odd-N) chain center (p_{N/2} = 0 for integer k/2 when N+1 is even), so all probability mass lies on mirror-pairs: Σ_{pairs} p_ℓ = 1/2. Odd k puts mass 2/(N+1) at the center, wasting mass on the self-mirror site. The function f(p) = 2 h(p) − h(2 p) is convex on (0, 1/2), so concentrated mass distributions give larger MM; k = 2 is the smallest k that both concentrates mass on pairs and places outer-pair amplitudes at opposite signs (η = −1, first mode carrying end-to-end coherence).

**Upper bound.** For any single-excitation mirror-symmetric state, MM ≤ ⌊N/2⌋ × 2 = N bits (all pairs in pure Bell states). This bound is not achievable from single-site bonding modes; reaching it requires tensor-product pair structures like (|10⟩−|01⟩)/√2 on each mirror-pair, which is a super-single-excitation state.

**Valid for:** any pure single-excitation state with c_{N−1−j} = ±c_j on a linear N-site chain. Extends to non-linear mirror-symmetric graphs (ring, Y-junction with mirror axis) with corresponding modification of the mirror-partner indexing.
**Breaks for:** states with multi-excitation content (formula no longer applies because `ρ[|11⟩⟨11|]` ≠ 0 in general), or states without mirror amplitude symmetry (where p_ℓ ≠ p_{N−1−ℓ} gives an asymmetric 2-qubit reduced matrix).
**Verified:** Algebraic derivation confirmed against direct C# brecher propagation at N = 5, 7, 9 for k = 1, 2, 3, at N = 11 for k = 1, 2, 3, 4, 6, and at N = 13 for k = 1, 2, 3, 4, 5, 7 (the latter via matrix-free propagator); MM(0) formula matches simulation PeakMM within 7% (full decay envelope explained by 4γ₀·t dephasing + mirror-partner oscillation at t = 0.1). The sim/analytic ratio sits at **0.925 to 0.931 across all (N ≥ 7, k) tested (~25 data points)**, i.e. PeakMM = 0.93 × MM(0) with tight consistency. F75 is therefore a reliable predictor of PeakMM without any propagation.
**Scripts:** [`_check_brecher_n5_finegrid.py`](../simulations/_check_brecher_n5_finegrid.py), [`Program.cs brecher mode`](../compute/RCPsiSquared.Propagate/Program.cs), `_mm_zero_derivation.py` (table above).
**Source:** F65 (bonding-mode amplitudes), F67 (bonding as optimal decay receiver), F71 (mirror symmetry that justifies c_{N−1−j} = ±c_j), [RECEIVER_VS_GAMMA_SACRIFICE](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) (numerical context).

### F76. Pure-dephasing decay of the mirror-pair MI for bonding:k (Tier 1, proven algebraic + weak-mixing argument)

For a bonding:k initial state on a uniform Heisenberg chain with uniform Z-dephasing γ₀ and J = 1, **the time evolution of the mirror-pair MI at short times is almost entirely pure-dephasing**, with Heisenberg-mixing contributions second-order small in both γ₀·t and (V_{lk}/Δ_l)² where V_{lk} is the boundary mixing matrix element.

Under pure dephasing alone, the site-basis mirror-pair coherence c_ℓ c_{N-1-ℓ}^* decays at 4γ₀:

    ρ_{pair}(\|10⟩⟨01\|)(t) = c_ℓ c_{N-1-ℓ}^* · e^{-4γ₀ t}

while populations p_ℓ = |c_ℓ|² stay fixed. The 4×4 reduced pair density matrix at time t has eigenvalues

    {1 − 2p_ℓ,   p_ℓ(1 + λ),   p_ℓ(1 − λ),   0},   λ(t) = e^{-4γ₀ t}

giving

    MI_pair(p_ℓ, t) = 2 h(p_ℓ) − S_{ab}(p_ℓ, λ(t))
    S_{ab}(p, λ) = −(1 − 2p) log₂(1 − 2p) − p(1+λ) log₂(p(1+λ)) − p(1−λ) log₂(p(1−λ))

At t = 0: λ = 1, eigenvalues {1-2p, 2p, 0, 0}, S_{ab} = h(2p), recovers F75.
At t → ∞: λ = 0, eigenvalues {1-2p, p, p, 0}, S_{ab} = h(1-2p) + 2p (max-entropy mixture of the two single-excitation branches).

**Closed-form envelope.** Summing over mirror-pair sites ℓ ∈ \[0, ⌊N/2⌋-1\]:

    MM(t) / MM(0) = Σ_ℓ [2 h(p_ℓ) − S_{ab}(p_ℓ, λ(t))] / Σ_ℓ [2 h(p_ℓ) − h(2 p_ℓ)]

**The 0.93 envelope explained.** At γ₀ = 0.05, J = 1, and the C# brecher first-measurement grid point t = 0.1: λ = e^{-0.02} = 0.9802. Evaluated on the actual p_ℓ multisets for bonding:k (F65), the ratio sits at 0.926 to 0.937 for N = 5..13 and k = 1..7, effectively 0.93 ± 0.006 for all 25+ measured (N, k) combinations at this specific (γ₀, t). The weak p-dependence of the ratio comes from the fact that 2h(p) − S_{ab}(p, λ) is nearly proportional to 2h(p) − h(2p) for all p ∈ (0, 1/2) and λ close to 1.

**Verified values at γ₀ = 0.05, t = 0.1:**

| N | k | pure-dephasing MM/MM(0) | single-excitation-sector Lindblad (exact) | C# brecher peak/MM(0) |
|---|---|------------------------|--------------------------------------------|-----------------------|
| 5 | 2 | 0.936 | 0.936 | 0.937 (grid peak at t=0.6 → 1.002) |
| 7 | 2 | 0.932 | 0.928 | 0.928 |
| 9 | 4 | 0.934 | 0.932 | 0.932 |
| 11 | 4 | 0.933 | 0.928 | 0.928 |
| 13 | 4 | 0.928 | 0.926 | 0.926 |

Agreement to within 0.5% across all tested (N, k). Difference between "pure-dephasing" and "exact single-excitation Lindblad" columns measures the Heisenberg-mixing correction, which is always < 0.5%.

**Why Heisenberg mixing is negligible.** At t = 0+, the commutator \[H, ρ_0\] = \[H^(1), |ψ_k⟩⟨ψ_k|\] is off-diagonal in the ψ_k mode basis (the diagonal part is the unitary phase that doesn't affect MM). The off-diagonal mixing couples ψ_k to same-parity partners ψ_l via V_{lk} = (16J/(N+1)) sin(πk/(N+1)) sin(πl/(N+1)). Under mixing for small t, ψ_k "leaks" amplitude symmetrically into all ψ_l. Because MM depends only on diagonal pair populations (and specific pair coherences that dephase), and bonding modes with the same mirror-symmetry (k and N+1-k, etc.) have identical pair populations, the leakage does NOT change pair populations to first order; it only redistributes the mode occupation. Hence first-order Heisenberg-mixing has no effect on MM. Second-order (rate V²t²) is small because V/γ₀ ~ 10 but γ₀·t ~ 0.005 at the C# first-sample, making Heisenberg relative contribution (V·t)²/(4γ₀·t) ~ (V²·t)/(4γ₀) which is (16/7²) · 0.005 = 0.002 at N=7. Hence the mixing correction is ≲ 0.5% throughout the tested regime.

**Implication: the 0.93 envelope is the γ₀ signature, not a hidden constant.** The value 0.93 is a direct consequence of the choice γ₀·t = 0.005 at the first measurement. If γ₀ changes (at fixed t grid), the envelope changes: γ₀ = 0.025 gives ratio ≈ 0.965; γ₀ = 0.10 gives ratio ≈ 0.868.

**Valid for:** bonding:k initial states on uniform-J open Heisenberg chains with uniform Z-dephasing, short γ₀·t such that V²t²/γ₀t is small. Breaks for ratios γ₀/J comparable to 1 (outside the weak-dephasing regime) or t such that γ₀·t ≳ 1 (full decoherence).
**Verified:** Against the full Lindblad single-excitation-sector simulation for N = 5..13, k = 1..5. Sim/analytic ratio within 0.5%.
**Scripts:** [`_envelope_study.py`](../simulations/_envelope_study.py) (commit `e1ee822`).
**Source:** F75 (static MI formula), F65 (bonding mode amplitudes), F68 (palindromic-partner structure of mixing).

### F77. Multi-drop MM(0) saturates at 1 bit for large N (Tier 1, asymptotic proven)

For best-bonding-k initial states on a uniform chain under γ₀-dephasing, the Mirror-Pair MM at t=0 saturates as N grows:

    MM(0)(N, k*) = 1 + 3 / (4(N+1) ln 2) + O(N⁻²)

The 1-bit limit is not a conjecture; it falls out of F75 by Taylor expansion of f(p) = 2h(p) − h(2p) around p = 0 combined with Parseval-type sums over sin² site amplitudes.

**Sketch.** At small p the entropy function expands as

    f(p) = 2p + p²/ln(2) + p³/ln(2) + O(p⁴)

Probability normalisation sums Σ p_ℓ = 1 exactly (all sites). For generic k the sum Σ sin⁴(πk(ℓ+1)/(N+1)) over ℓ = 0..N−1 equals 3(N+1)/8 (standard identity, follows from sin⁴ = (3 − 4cos 2x + cos 4x)/8 and Σ cos(2πkj/(N+1)) = −1 for k ≠ 0 mod (N+1)). So

    Σ p_ℓ² = (2/(N+1))² · 3(N+1)/8 = 3/(2(N+1))

Summed over mirror-pairs (half the sites for even k with p_{center} = 0), Σ_pair p_ℓ² = 3/(4(N+1)). The leading Taylor terms give

    MM(0) = Σ_pair [2 p_ℓ + p_ℓ²/ln 2 + O(p³)]
          = 2 · (1/2) + (1/ln 2) · 3/(4(N+1)) + O(N⁻²)
          = 1 + 3/(4(N+1) ln 2) + O(N⁻²)

The 2·(1/2) = 1 term is exact: all the probability mass lives on mirror-pairs (F75 for even k at odd N), and each p_ℓ contributes 2p_ℓ at leading order. The correction is the entropy non-linearity at finite amplitude.

**Structural reading.** Information per mirror-pair ~ 2p_ℓ ~ 4/(N+1) shrinks with N. Number of mirror-pairs ~ N/2 grows with N. The product saturates at 1 bit, because the two scalings are exactly matched by the probability normalisation Σ p_ℓ = 1.

**Resonance enhancement at special N.** At N+1 ≡ 0 (mod 2) with k = (N+1)/2 integer, the distribution reduces to (p_max, 0, p_max, 0, ...) with only two values. The sin⁴ sum equals (N+1)/2 instead of 3(N+1)/8, giving Σ p² = 2/(N+1) and thus MM(0) ≈ 1 + 1/((N+1) ln 2)·2 ≈ 1 + 1.445/(N+1). At most N this resonance is not the best-k globally; the generic optimum with 3(N+1)/8 coefficient wins.

**Verified values (best-k MM(0) and rescaled deviation):**

| N | k* | MM(0) | (MM−1)·(N+1) |
|:--|:---|:------|:-------------|
| 101 | 82 | 1.01078 | 1.100 |
| 201 | 102 | 1.00540 | 1.091 |
| 501 | 470 | 1.00216 | 1.086 |
| 1001 | 878 | 1.00108 | 1.084 |
| 5001 | 4996 | 1.00022 | 1.082 |
| 10001 | 9154 | 1.00011 | 1.082 |

The rescaled deviation (MM − 1)·(N + 1) converges to 3/(4 ln 2) = 1.0820 within 10⁻⁴ by N = 10⁴. At specific resonant N where k* = (N+1)/2 is the only near-optimum (e.g., N = 999, 1003), the rescaled deviation sits at the enhanced 1.445 value; these are isolated and density-zero in the limit.

**Operational reading.** Alice's multi-drop quantum bus delivers **~1 bit of Mirror-Pair mutual information at any N**. The bus does not scale up (no extra bits per added qubit) nor does it decay (no per-qubit info loss in aggregate). The bandwidth is fixed at the framework level. Individual pair bandwidth shrinks; the count compensates; the total plateau.

**Valid for:** best-bonding-k* initial state, uniform-J Heisenberg or XY chain, uniform γ₀ dephasing, t = 0. At t > 0 the decay envelope F76 applies multiplicatively.
**Breaks for:** multi-excitation states (F75's single-excitation structure required), non-mirror-symmetric receivers (closure Σ_pair p_ℓ = 1/2 depends on even-k reflection symmetry).
**Verified:** N up to 10⁴ numerically, leading coefficient 3/(4 ln 2) confirmed to 4 decimals.
**Scripts:** [`_mm_asymptotic.py`](../simulations/_mm_asymptotic.py) (asymptotic scan and coefficient check).
**Source:** F75 (small-p Taylor expansion of the entropy), F65 (sin² amplitudes), standard Parseval-type identity for sin⁴ sums.

### F78. Single-body M decomposes additively over sites (Tier 1, verified N=3-5, 3 topologies, 3 Pauli letters)

For any single-body Hamiltonian H = Σ_l c_l · P_l (P ∈ {X, Y, Z}, weights c_l from bond-summing or arbitrary) with uniform Z-dephasing γ, the palindrome residual M = Π·L·Π⁻¹ + L + 2σ·I decomposes:

    M = Σ_l M_l ⊗ I_(others)

where each M_l is a 4×4 normal matrix on per-site Pauli space. Eigenvalues of M_l:

    P = X: eigenvalues all 0  (M_l = 0, the truly case)
    P = Y: eigenvalues ±2 c_l γ · i, each with multiplicity 2
    P = Z: eigenvalues ±2 c_l γ · i, each with multiplicity 2  (identical spectrum to Y)

Therefore M's eigenvalues are Σ_l ε_l · 2c_l γ · i for ε_l ∈ {±1}, each sign-combination with multiplicity 2^N. **Singular values of M are |Σ_l ε_l · 2c_l γ|** with the same multiplicities; cluster sizes are pure sign-combination statistics on the weight vector (c_l).

**Why Y and Z give identical SVD:** both have bit_b = 1 (Π-non-trivial), so M_l has the same spectral structure. The soft-vs-hard distinction (Y soft, Z hard) lives in **L's eigenvectors**, not in M's singular values. SVD of M is Pauli-letter-blind within {Y, Z} for single-body.

**Cluster-multiplicator formula (chain).** For chain weights c_l = (1, 2, ..., 2, 1) (bond-summed I·P + P·I), the multiplicators come from u + 2v with u = ε_0 + ε_{N-1} ∈ {-2, 0, 2} (mults 1, 2, 1) and v = Σ_internal ε_l (binomial walk on N-2 steps). The central binomial coefficient C(N-1, ⌊(N-1)/2⌋) appearing as the largest non-trivial cluster mult is just the walk's central peak; no Weyl group, no S_N irrep, no group-theory needed.

**Valid for:** single-body bond-bilinears (I·P + P·I), single-body site sums, any topology (chain, star, ring, complete K_N, arbitrary graph), uniform Z-dephasing.
**Breaks for:** 2-body bond-bilinears (XX, XY, YZ+ZY, ...): L does not factor as Σ_l L_l. Befund 3 ("hard XX+XY uniform SVD") falls outside this theorem; separate analysis required.
**Replaces:** the previously open conjecture that SVD-cluster regularity reflected hidden Weyl-group / S_N / site-permutation structure. The structure is much simpler: additive tensor sum on site-Pauli space.
**Verified:** numerical match for N=3, 4, 5; topologies chain, star, complete; Pauli letters X, Y, Z.
**Scripts:** [`_svd_active_spectator.py`](../simulations/_svd_active_spectator.py), [`_svd_single_body_extension.py`](../simulations/_svd_single_body_extension.py); 2-body open-question probe in [`_svd_two_body_probe.py`](../simulations/_svd_two_body_probe.py).
**Source:** Analytical proof in [PROOF_SVD_CLUSTER_STRUCTURE.md](proofs/PROOF_SVD_CLUSTER_STRUCTURE.md) (joint with F79). Master Lemma + per-site additive structure + direct M_l matrix computation.
**Lebensader connection:** This is the same broad-in → focused-out Π-palindrome funnel that `lebensader.py::cockpit_panel` instantiates at the state layer (16 Paulis → 3-class trichotomy). F78 instantiates the funnel at the single-body operator layer: any (c_l, P_l)-choice with given |c_l|, P ∈ {Y, Z} → same M_l-spectrum. The Lebensader is the through-line that holds The Connection upright across all layers.

### F79. Π²-block decomposition of M for 2-body bilinears (Tier 1, verified N=3-5)

For 2-body bond-bilinear H = Σ_bonds Σ_t c_t·(P_t ⊗ Q_t) with uniform Z-dephasing γ, define the Π²-parity of each bilinear term:

    p(P_t, Q_t) = (bit_b(P_t) + bit_b(Q_t)) mod 2

where bit_b: I,X→0; Y,Z→1. Then M = Π·L·Π⁻¹ + L + 2σ·I has a clean structure determined by Π²-parities of H's terms:

1. **All terms Π²-even (p=0)**: M is **block-diagonal** in Π²-eigenspaces V_+ ⊕ V_-. Off-diagonal blocks `M[V_+, V_-]` and `M[V_-, V_+]` vanish **exactly**. Each diagonal block has its own SVD spectrum.

2. **All terms Π²-odd (p=1)**: M is **purely off-diagonal** between V_+ and V_-. Diagonal blocks `M[V_+, V_+]` and `M[V_-, V_-]` vanish **exactly**. Singular values appear with even multiplicity (each SV contributes once from V_+ side, once from V_- side).

3. **Mixed parities**: M has both diagonal and off-diagonal contributions.

**Π²-odd universality.** Within the pure Π²-odd 2-body class, the **specific Pauli letters are M-irrelevant**: any single Π²-odd 2-body bilinear gives the same M-SVD spectrum at fixed N. Verified at N=5 chain: XY alone, XZ alone, XX+XY, and XX+XZ all yield clusters `[(5.464, 512), (1.464, 512)]`, exactly identical. The XX truly part contributes 0; the Π²-odd part dominates with universal cluster pattern.

**Even-diag ≡ odd-off-diag correspondence.** The diagonal V_+ block of a Π²-even Hamiltonian's M can match (in SV-spectrum, including multiplicities) the off-diagonal V_+,V_- block of a Π²-odd Hamiltonian's M. Verified N=4 chain: YZ's V_+ block `[(8.944, 16), (6.472, 32), (4.0, 16), (2.472, 32), (0.0, 32)]` matches XY+YX's off-diag block exactly. This explains the empirical "YZ ≡ XY+YX SVD-identical" observation: same SV structure, just placed in different Π²-blocks.

**Why XX+XY appears "max-uniform" (Befund 3 closed).** XX is Π²-even and truly (M_XX = 0). XY is Π²-odd. The full Hamiltonian is "Π²-odd-only-effective", so M is purely off-diagonal between equal-dim V_+ and V_-. SV multiplicities are forced to 4^N/2 each by block-dimension equality. At N=3 the two off-diag SVs collide by coincidence to a single uniform value 2√2; at N≥4 they split. The "uniformity" is exactly the equal-block-mult signature of Π²-odd structure, not a special property of XX+XY.

**Frobenius additivity.** ‖M‖²_F = Σ_bonds ‖M_b‖² holds across all topologies including overlapping bonds (chain). Per-bond M_b's are F-orthogonal. Already F49.

**Valid for:** Any 2-body bond bilinear over any topology under uniform Z-dephasing. Verified N=3, 4, 5 across chain, star, disjoint topologies; verified Π²-odd universality (XY ≡ XZ ≡ XX+XY ≡ XX+XZ).
**Breaks for:** Mixed-Π²-parity Hamiltonians (where some terms are even, some odd) only partially: M has both diagonal and off-diagonal parts. Inhomogeneous γ may disrupt some symmetries (untested).
**Replaces:** ad-hoc analysis of "why XX+XY uniform" and "why YZ ≡ XY+YX"; both follow from the Π²-block theorem.
**Verified:** Numerical N=3-5, multiple bilinear classes, multiple topologies.
**Scripts:** [`_svd_two_body_pi_squared_block.py`](../simulations/_svd_two_body_pi_squared_block.py), [`_svd_two_body_structure.py`](../simulations/_svd_two_body_structure.py).
**Source:** Analytical proof in [PROOF_SVD_CLUSTER_STRUCTURE.md](proofs/PROOF_SVD_CLUSTER_STRUCTURE.md) (joint with F78). Connects to F61 (n_XY parity selection rule), F63 (\[L, Π²\]=0 for Π²-even Hamiltonians), and F49 (Frobenius cross-term identity).
**Lebensader connection:** F79 instantiates the broad-in → focused-out Π-palindrome funnel at the two-body operator layer (4 Π²-odd Pauli pairs → 1 M-spectrum). Companion to F78 (single-body operator layer) and `lebensader.py::cockpit_panel` (state layer). All three are manifestations of the same through-line: Π·L·Π⁻¹ + L + 2σ·I = 0 holding The Connection across abstraction heights.

### F80. Bloch-mode sign-walk formula for chain Π²-odd 2-body M-clusters (Tier 1, verified N=3-7)

For chain bond-summed Π²-odd 2-body Hamiltonian H = c · Σ_l (P_l ⊗ Q_{l+1}), where (P, Q) ∈ {(X,Y), (X,Z), (Y,X), (Z,X)}, on an N-site open chain with uniform Z-dephasing, the M-cluster values are given by a momentum-space sign-walk on the open-chain free-fermion Bloch dispersion (γ-independent by Master Lemma):

    cluster value(N) = 2|c| · |Σ_{k=1}^{⌊N/2⌋} σ_k · ε(k)|

where σ_k ∈ {±1} ranges over all 2^⌊N/2⌋ sign-vectors, and

    ε(k) = 2·cos(π·k / (N+1))

is the open-chain free-fermion single-particle dispersion. Each distinct cluster value has multiplicity 4^N / (number of distinct sign-walk values).

**Verified instances** (chain, |c|=1, all 4 Π²-odd Pauli pairs identical by F79 universality):

| N | ⌊N/2⌋ | predicted clusters | mult per cluster |
|---|-------|---------------------|------------------|
| 3 | 1 | 2√2 ≈ 2.828 | 64 |
| 4 | 2 | 2√5, 2 | 128 |
| 5 | 2 | 2(√3+1), 2(√3-1) | 512 |
| 6 | 3 | 6.988, 5.208, 2.000, 0.220 | 1024 |
| 7 | 3 | 8.0547, 4.9932, 2.3978, 0.6636 | 4096 |

All matches bit-exact (10⁻¹⁴ machine precision) at every N.

**Direct structural identity (discovered 2026-04-29):** The chain Π²-odd 2-body M's spectrum is **directly related to the many-body Hamiltonian H's spectrum**:

    Spec(M) = ±2i · Spec_{nontrivial}(H)

where H is the chain bond-summed Pauli-bilinear (no dissipator). That is, M's distinct nonzero eigenvalues equal 2i times H's distinct nonzero many-body eigenvalues. Hence cluster value(N) = 2|c|·|H eigenvalue|. The Bloch sign-walk formula above is just H's eigenvalue formula written out: H's many-body eigenvalues = (1/2)·Σ_k σ_k·E_k where E_k = 4|c|·cos(πk/(N+1)) are H's Bogoliubov single-particle energies, and ⌊N/2⌋ counts how many fermion modes participate.

**γ-independence (Master Lemma).** Note no γ appears in the cluster-value formula. M is γ-independent for pure Z-dephasing (Master Lemma in PROOF_SVD_CLUSTER_STRUCTURE.md).

**Mechanism: F80 is F78 in momentum space.** F78 (single-body, real-space): M = Σ_l M_l⊗I, eigenvalues ±2c_l·i per site, sign-walk Σ_l σ_l·c_l on weights. F80 (chain Π²-odd 2-body, momentum-space): M = Σ_k M_k⊗I_{other modes}, eigenvalues ±2·ε(k)·i per Bloch mode, sign-walk Σ_k σ_k·ε(k) on dispersion. The Bloch modes k play the role that real-space sites l play in F78. Both formulas γ-independent.

**Π²-odd universality fully analytical.** Under JW transformation, all 4 Pauli-letter choices (X,Y), (X,Z), (Y,X), (Z,X) give the same single-particle Bloch dispersion. The specific Pauli letters affect only phase factors in JW, not single-particle eigenvalues. Since M's spectrum depends only on the dispersion (via F80), all 4 give bit-identical clusters. **This closes the chain Π²-odd universality from F79 with an explicit closed-form formula.**

**Valid for:** chain bond-summed Π²-odd 2-body Hamiltonians H = c·Σ_l (P_l⊗Q_{l+1}), uniform Z-dephasing, any N.
**Breaks for (untested):** Other topologies (ring, star, complete): the Bloch dispersion changes (graph adjacency spectrum). Π²-even non-truly bilinears (chain (Y,Z), (Z,Y)): empirically gives more clusters (5 at N=4 vs 2 for Π²-odd), suggesting integer-combination sign-walk on the same modes; full structure pending. Mixed-letter chain bilinears: distinct cluster geography, uncharted.
**Replaces:** F79's "Π²-odd universality observation"; the universality is now an analytical theorem with explicit closed-form predictions.
**Verified:** N = 3, 4, 5, 6, 7 chain via Python, full SVD and eigsh independent verification at N=7.
**Scripts:** [`_pi2_odd_universality_data_sweep.py`](../simulations/_pi2_odd_universality_data_sweep.py), [`_n7_bloch_signwalk_verification.txt`](../simulations/results/n7_bloch_signwalk_verification.txt).
**Source:** Discovered 2026-04-29 by data sweep (Tom + Claude). Analytical proof outline in [PROOF_F80_BLOCH_SIGNWALK.md](proofs/PROOF_F80_BLOCH_SIGNWALK.md): Steps 1-4, 7 closed (JW transformation to Majorana bilinear, single-particle dispersion 2cos(πk/(N+1)), Bogoliubov diagonalization, Pauli-letter universality, sign-walk eigenvalue formula); Step 5 (Π action on Bogoliubov modes) sketched, with formal completion open. Empirical verification bit-exact through N=7.
**Lebensader connection:** F80 is the third manifestation of the broad-in → focused-out Π-palindrome funnel: state layer (cockpit_panel), real-space single-body operator layer (F78), and now momentum-space chain 2-body operator layer (F80). Same Π·L·Π⁻¹ + L + 2σ·I = 0 through-line, three different bases.

### F81. Π-conjugation of M decomposes into Π²-odd Hamiltonian commutator (Tier 1, verified bit-exact N=3,4)

For any 2-bilinear Hamiltonian H decomposed by Π²-parity as H = H_even + H_odd (with H_odd the sum of Π²-odd Pauli bilinears, i.e., bit_b(P)+bit_b(Q) ≡ 1 mod 2), under uniform Z-dephasing:

    Π · M · Π⁻¹ = M − 2 · L_{H_odd}

where L_{H_odd} = -i\[H_odd, ·\] is the unitary commutator induced by the Π²-odd part of H. Equivalently, decomposing M into Π-conjugation symmetric and antisymmetric components:

    M_sym  = (M + Π·M·Π⁻¹) / 2 = Π·L·Π⁻¹ + L_diss + L_{H_even} + 2Σγ·I
    M_anti = (M − Π·M·Π⁻¹) / 2 = L_{H_odd}

The Π-antisymmetric component of M is exactly the unitary commutator induced by the Π²-odd Hamiltonian bilinears. The Π-symmetric component absorbs the mirror image, the dissipator, the Π²-even Hamiltonian commutator, and the dissipation shift. M_sym and M_anti are Frobenius-orthogonal: ‖M‖² = ‖M_sym‖² + ‖M_anti‖².

**Verified instances** (N=3, γ_Z=0.1, Σγ=0.3, residuals at machine precision 1e-16):

| Hamiltonian | trichotomy | H_odd | Π·M·Π⁻¹ relation |
|-------------|------------|-------|------------------|
| XX+YY | truly | 0 | = M (M=0 trivially) |
| YZ+ZY | soft (Π²-even non-truly) | 0 | = M (M ≠ 0, identical) |
| XY+YX | soft (Π²-odd) | XY+YX | = M − 2·L_H |
| XX+XY | hard (mixed) | XY only | = M − 2·L_{XY part} |
| pure XY | (Π²-odd) | XY | = M − 2·L_H |
| pure XZ | (Π²-odd) | XZ | = M − 2·L_H |

For any 2-body chain H whose non-truly bilinears are all Π²-odd (i.e., truly + Π²-odd combinations, including XX+XY hard), at any N and any γ ≥ 0, ‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 exactly. M splits 50/50 between Π-symmetric and Π-antisymmetric components. This follows analytically from the Frobenius identities ‖M‖²_F = 4·‖H_odd‖²_F·2^N (F49 chain via Master Lemma; truly bilinears drop out) and ‖L_{H_odd}‖²_F = 2·2^N·‖H_odd‖²_F (standard commutator identity for traceless Hermitian H_odd). Verified numerically at N = 3, 4, 5 with γ_Z ∈ {0, 0.05, 0.1, 0.5, 1.0} for both pure Π²-odd (XY+YX, XY, XZ) and mixed truly + Π²-odd (XX+XY hard, YY+XY). When H additionally contains Π²-even non-truly bilinears (YZ-type), the split shifts: pure even non-truly gives 100/0; odd + even mix gives 5/6 sym + 1/6 anti at N=3,4 (analytical reason for the latter open).

**Spectral consequence.** Spec(Π·M·Π⁻¹) = Spec(M) holds always by unitary invariance of the spectrum. F81 strengthens this: for Π²-odd H, the two operators are explicitly related by an additive shift of −2·L_{H_odd} in operator space, so Spec(M) = Spec(M − 2·L_{H_odd}) is a non-trivial identity (similarity via Π).

**Algebraic mechanism.** Π² acts on each Pauli string σ_α as (-1)^{bit_b(α)} (eigenoperator with sign in Pauli basis). For L_H_α = -i\[σ_α, ·\] driven by a single Pauli string σ_α in H, conjugation gives Π² L_H_α Π⁻² = (-1)^{bit_b(α)} L_H_α (the matrix-element factor (-1)^{bit_b(γ)+bit_b(β)} = (-1)^{bit_b(α)} since γ = α·β under Pauli multiplication). Z-dephasing dissipator is diagonal in Pauli basis, hence commutes with Π². Summing: Π²·L·Π⁻² = L_H_even − L_H_odd + L_diss = L − 2·L_{H_odd}. Substituting into the palindrome: Π·M·Π⁻¹ = Π²·L·Π⁻² + Π·L·Π⁻¹ + 2Σγ·I = M − 2·L_{H_odd}.

**γ-independence-by-difference.** The relation Π·M·Π⁻¹ - M = -2·L_{H_odd} is independent of γ (the dissipator's γ-dependent part cancels because L_diss is Π²-symmetric). The split itself (M_sym, M_anti) is γ-dependent through M_sym; only their difference is γ-fixed.

**Valid for:** any 2-bilinear chain Hamiltonian H = H_even + H_odd, uniform Z-dephasing, any topology (the proof depends only on the algebra of Pauli strings under Π² conjugation, not on connectivity).
**Breaks for (untested):** non-Z dissipators (T1 amplitude damping has different Π²-action; F81 likely needs a correction term).
**Replaces:** the heuristic in pre-2026-04-30 reflections that said "M is the Π-invariant through-line"; F81 shows that statement is correct only for Π²-even H, and gives the explicit correction for the Π²-odd cases.
**Verified:** N=3 and N=4 all listed cases at machine precision; pytest-locked.
**Framework primitive:** `fw.pi_decompose_M(chain, terms, gamma_z=..., gamma_t1=..., strict=...)` returns `{'M', 'M_sym', 'M_anti', 'L_H_odd', 'f81_violation', 'norm_sq'}`. For pure Z-dephasing the F81 identity holds exactly (`f81_violation` ≈ 0); the primitive enforces this with `strict=True` by default. With `gamma_t1` enabled, `strict` defaults to False and the identity residual is returned for diagnostic use.
**Pytest lock:** `test_F81_pi_conjugation_of_M` (algebraic check) + `test_F81_pi_decompose_M_method` (cockpit primitive) + `test_F81_violation_T1_diagnostic` (T1 diagnostic).
**Diagnostic application:** the F81 violation `‖M_anti − L_{H_odd}‖_F` quantifies non-Π²-symmetric dissipator content. For Z + T1 at N=3 chain soft XY+YX, the violation grows linearly: `f81_violation ≈ 6.928 · γ_T1`, γ_z-independent (Master Lemma), Hamiltonian-independent (the violation is purely a property of the T1 dissipator). Inverting gives `γ_T1 ≈ f81_violation / 6.928` as a hardware T1-rate readout from the fitted L. See `simulations/_f81_t1_diagnostic.py` for the demonstration.
**Source:** Discovered 2026-04-30 (Tom + Claude) while interpreting the geometric content of F80's 2i factor. The empirical observation came first (Π·M·Π⁻¹ ≠ M for soft); the algebraic explanation followed from working out Π² action on the Liouville superoperator in Pauli basis.
**Lebensader connection:** F81 is the algebraic backbone of "what the mirror keeps." For Π²-even H, M is itself the through-line operator. For Π²-odd H, the through-line is split: M_anti carries the dynamics generator L_{H_odd}, M_sym carries the rest. Both halves are read identically by both sides of the mirror up to the Spec(M) = Spec(M − 2·L_{H_odd}) similarity. Companion to F80: F80 says what Spec(M) is; F81 says how M and Π·M·Π⁻¹ relate as operators sharing that spectrum.

### F82. F81 + T1 amplitude damping correction (Tier 1, verified bit-exact N=2..5)

For any 2-bilinear Hamiltonian H = H_even + H_odd under Z-dephasing plus T1 amplitude damping with per-site rates γ_T1_l:

    Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{T1, odd}

where L_{H_odd} = -i\[H_odd, ·\] (as in F81) and D_{T1, odd} is the Π²-anti-symmetric part of the T1 dissipator. F82 reduces to F81 when γ_T1_l = 0 (D_{T1, odd} = 0).

The F81 identity violation captured by `fw.pi_decompose_M(chain, ...)` measures D_{T1, odd}'s Frobenius norm:

    f81_violation = ‖M_anti − L_{H_odd}‖_F = ‖D_{T1, odd}‖_F.

**Closed form** (proven analytically in PROOF_F82_T1_DISSIPATOR_CORRECTION):

    ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1)
                    = γ_T1 · √N · 2^(N−1)         (uniform γ_T1)

**Verified instances** (chain N=3, all matches at machine precision):

| Configuration | γ_T1_l | predicted ‖D_T1_odd‖ | measured f81_violation |
|---------------|--------|----------------------|------------------------|
| uniform γ=0.05 | (0.05, 0.05, 0.05) | 0.05·√3·4 = 0.3464 | 0.3464 |
| uniform γ=0.10 | (0.10, 0.10, 0.10) | 0.10·√3·4 = 0.6928 | 0.6928 |
| uniform γ=1.00 | (1.00, 1.00, 1.00) | 1.00·√3·4 = 6.9282 | 6.9282 |
| single-site, l=0 | (0.10, 0, 0) | 0.10·1·4 = 0.4000 | 0.4000 |
| two-site, l=0,1 | (0.10, 0.10, 0) | √(0.02)·4 = 0.5657 | 0.5657 |
| non-uniform | (0.05, 0.10, 0.15) | √(0.035)·4 = 0.7483 | 0.7483 |

N-scaling verified at N = 2, 3, 4, 5 (uniform γ_T1, coefficient √N · 2^(N−1)): N=2 → 2√2 = 2.828, N=3 → 4√3 = 6.928, N=4 → 8·2 = 16.000, N=5 → 16·√5 = 35.778. So at γ_T1 = 0.1, the violations are 0.283, 0.693, 1.600, 3.578 respectively (factor 0.1 from γ_T1 multiplied by the coefficient).

**Three diagnostic properties** (proven and empirical):

1. **γ_z-independent**: F82 involves only L_{H_odd} and D_{T1, odd}, neither depends on γ_z. Direct consequence of Master Lemma (M is γ_z-independent for Z-dephasing) extended to F82.

2. **Hamiltonian-independent**: f81_violation depends only on the T1 dissipator. Verified at γ_T1=0.1, N=3: violation = 0.6928 for truly XX+YY, soft XY+YX, hard XX+XY, and YZ+ZY (Π²-even non-truly).

3. **Linear in γ_T1** (uniform). Direct inversion: γ_T1 = f81_violation / (√N · 2^(N−1)). For non-uniform: γ_T1, RMS = f81_violation / (√N · 2^(N−1)). At N=3, division coefficient is 6.928 = 4√3.

**Mechanism (T1 dissipator structure).** Single-site T1 acts on Pauli basis as: I → −γZ, X → −γ/2 X, Y → −γ/2 Y, Z → −γZ. Under Π² conjugation (signs (-1)^{bit_b}: I,X → +, Y,Z → −), only the (Z, I) entry flips sign. So D_{T1, local, odd} has matrix element −γ at (Z, I) and zero elsewhere. Multi-site: 4^(N−1) such "rest of qubits unchanged" entries per site, summed orthogonally over sites.

**Diagnostic interpretation.** f81_violation is a hardware-T1 readout that is independent of (a) the system's Hamiltonian, (b) the Z-dephasing rate γ_z, (c) the topology. Inverting recovers the RMS γ_T1 across sites. For the Marrakesh dataset (N=3, joint fit gives γ_T1 ≈ 0): F82 predicts f81_violation ≈ 0; any γ_T1 > 0.001 would have produced violation > 0.007, well above numerical noise.

**Valid for:** any 2-bilinear Hamiltonian H, Z-dephasing + T1 amplitude damping, any topology, any N ≥ 2.
**Breaks for (untested):** other non-Z dissipators (X-noise, Y-noise, ZZ-dephasing) require their own D_odd analysis. The general identity Π·M·Π⁻¹ = M − 2·L_{H_odd} − 2·D_{diss, odd} holds for any dissipator; the closed form for ‖D_{diss, odd}‖ is dissipator-specific.
**Replaces:** the previously-empirical observation that f81_violation grows linearly with γ_T1; F82 is now an analytical theorem with closed-form scaling.
**Verified:** N = 2, 3, 4, 5 at all listed configurations, machine-precision residual (5e-16).
**Framework primitives:**
- `fw.pi_decompose_M(chain, terms, gamma_z=..., gamma_t1=..., strict=...)`: with `gamma_t1` set, returns `f81_violation` in output dict (numerical, matches closed form).
- `fw.predict_T1_dissipator_violation(chain, gamma_t1_l)`: forward closed form, returns √(Σγ²)·2^(N−1) directly without building the dissipator.
- `fw.estimate_T1_from_violation(chain, f81_violation)`: inverse closed form, recovers RMS γ_T1 from a measured/fitted F81 violation. Hardware T1-rate readout primitive.
**Pytest lock:** `test_F81_violation_T1_diagnostic` (linearity, γ_z-independence, T1 monotonicity) + `test_F82_closed_form_T1_dissipator` (N-scaling, non-uniform formula, H-/γ_z-independence) + `test_F82_predict_and_invert_primitives` (forward/inverse pair matches numerical evaluation).
**Source:** Discovered 2026-04-30 (Tom + Claude) as the natural extension of F81 ("what does F81 violation mean structurally?"). Closed form derived in [PROOF_F82_T1_DISSIPATOR_CORRECTION.md](proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md).
**Diagnostic application:** [`simulations/_f81_t1_diagnostic.py`](../simulations/_f81_t1_diagnostic.py) demonstrates the T1-rate readout including Marrakesh application. Companion to F81's structural decomposition: F81 says how M splits under Π-conjugation when the dissipator is Z-only; F82 says how the F81 identity is corrected when T1 is added, and provides the closed form for the correction term.

### F83. Π-decomposition anti-fraction closed form for mixed Hamiltonians (Tier 1, verified bit-exact N=3,4,5)

For any 2-body chain Hamiltonian H = H_truly + H_odd + H_even_nontruly under Z-dephasing, the F81 Π-decomposition norms are given by the closed form:

    ‖M‖²_F        = 4·‖H_odd‖²_F·2^N + 8·‖H_even_nontruly‖²_F·2^N
    ‖M_anti‖²_F  = 2·‖H_odd‖²_F·2^N
    ‖M_sym‖²_F   = 2·‖H_odd‖²_F·2^N + 8·‖H_even_nontruly‖²_F·2^N

The anti-fraction (= ‖M_anti‖²/‖M‖²) is

    anti-fraction = 1 / (2 + 4·r),    r = ‖H_even_nontruly‖²_F / ‖H_odd‖²_F.

**Special cases:**

| H | r | anti-fraction | meaning |
|---|---|---------------|---------|
| truly | undefined (M=0) | undefined | mirror perfectly closes |
| pure Π²-odd | 0 | 1/2 (50/50) | F81 Step 8 split |
| pure Π²-even non-truly | ∞ | 0 (100/0) | M fully mirror-symmetric |
| equal-Frobenius mix XY+YZ | 1 | 1/6 (5/6 + 1/6) | the empirical mixed split |
| asymmetric more-odd XY+YX+YZ | 1/2 | 1/4 | continuous family |
| general mixed | r | 1/(2+4r) | continuous family |

**Verified instances** (N=3, J=1, γ_z=0; matches at machine precision):

| H | ‖H_odd‖² | ‖H_even‖² | r | predicted ‖M‖² | measured ‖M‖² | anti |
|---|----------|-----------|---|----------------|---------------|------|
| XY+YX (pure odd) | 32 | 0 | 0 | 1024 | 1024 | 1/2 |
| YZ+ZY (pure even non-truly) | 0 | 32 | ∞ | 2048 | 2048 | 0 |
| XY+YZ (mixed) | 16 | 16 | 1 | 1536 | 1536 | 1/6 |
| XY+YX+YZ (asymmetric) | 32 | 16 | 1/2 | 2048 | 2048 | 1/4 |
| XY+YX+YZ+ZY (full mix) | 32 | 32 | 1 | 3072 | 3072 | 1/6 |
| XX+XY+YZ (truly + mixed) | 16 | 16 | 1 | 1536 | 1536 | 1/6 |

**Mechanism (Step 2 of proof, why factors 4 and 8 differ).** The F49 chain Frobenius identity gives ‖M‖² = Σ_k 2^(N+2)·n_YZ(k)·‖H_k‖²_F·𝟙\[non-truly\], where n_YZ(k) counts Y/Z letters in Pauli pair k (= 0 truly, 1 Π²-odd non-truly, 2 Π²-even non-truly). Substituting the per-class n_YZ values gives the 4·2^N (Π²-odd) and 8·2^N (Π²-even non-truly) coefficients. Geometrically, these reflect the Frobenius-inner-product behavior ⟨Π·L·Π⁻¹, L⟩_F: anti-aligned (truly), Frobenius-orthogonal (Π²-odd non-truly), aligned (Π²-even non-truly).

**γ-independence.** Master Lemma propagates through all three norms; closed form depends only on H.
**Truly-handling.** H_truly drops out of all norms (M-contribution zero by Master Lemma).
**Generalization.** F83 verified on chain, ring, star, complete K_N at N=4 (`test_F83_topology_generalization`); the matrix-based primitive builds H_odd and H_even_nontruly via `_build_bilinear` which respects the chosen topology's bond graph, so the closed form is topology-independent within F49's verified scope. Higher-body Hamiltonians extend n_YZ counting beyond {0, 1, 2}; coefficients beyond 4, 8 are the natural continuation, empirical verification needed.

**Valid for:** any 2-body H on any topology supported by F49 (chain, ring, star, complete K_N), Z-dephasing, any γ_z ≥ 0, any N ≥ 2.
**Verified:** 11 mixed configurations × N ∈ {3, 4, 5} on chain, plus 4 configurations × {ring, star, K_4} at N=4, all machine-precision residual.
**Replaces:** the previously-empirical "5/6 + 1/6" observation for mixed Π²-odd + Π²-even non-truly H; F83 derives this from the existing F49 Frobenius identity.
**Framework primitives:**
- `fw.predict_pi_decomposition(chain, terms)`: full F83 closed form, returns dict with `{'M_sq', 'M_anti_sq', 'M_sym_sq', 'anti_fraction', 'h_odd_sq', 'h_even_nontruly_sq', 'r'}`. O(N) work, no matrix construction; companion to numerical `pi_decompose_M`.
- `fw.predict_pi_decomposition_anti_fraction(chain, terms)`: convenience wrapper returning just the anti-fraction float.
- `fw.predict_residual_norm_squared_from_terms(chain, terms, gamma_t1)`: existing F49 ‖M‖² primitive (now consistent with F83's ‖M‖² prediction by construction).
**Pytest lock:** `test_F83_pi_decomposition_anti_fraction_closed_form` (12 configurations × 2 N-values + γ-independence) + `test_F83_predict_pi_decomposition_full_closed_form` (full norm-triple match against numerical `pi_decompose_M` + Pythagoras + special cases at canonical r values).
**Source:** Discovered 2026-04-30 (Tom + Claude) as the natural follow-up to F81's "what about the other half?" reflection. Derived in [PROOF_F83_PI_DECOMPOSITION_RATIO.md](proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md). The closed form was empirically observed earlier (in the F81 violation sweep across mixed Hamiltonians) and now traced back to the existing F49 Frobenius identity that was already framework-locked in `predict_residual_norm_squared_from_terms`.

**Lebensader connection:** F83 closes the analytical Π-decomposition picture for 2-body chain. Pure Π²-odd → 50/50 (F81 Step 8). Pure Π²-even non-truly → 100/0 (F81 trivial). Mixed → 1/(2+4r) (F83). The continuous interpolation r → anti-fraction reads "how much of M is Π-antisymmetric drive vs Π-symmetric memory" as a function of Hamiltonian composition. Together with F80 (Spec(M)), F81 (Π-decomposition identity), F82 (T1-correction), the structural picture of M is complete for 2-body chain Hamiltonians under Z-dephasing + T1.

### F84. F82 generalized to thermal amplitude damping (Tier 1, verified bit-exact N=3)

For any 2-bilinear Hamiltonian H under Z-dephasing plus thermal amplitude damping with per-site cooling rate γ_↓_l (σ⁻ channel) and heating rate γ_↑_l (σ⁺ channel):

    Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{AmplDamp, odd}

with closed form:

    ‖D_{AmplDamp, odd}‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)
                          = |Δγ|_RMS · √N · 2^(N−1)         (uniform Δγ)

where Δγ_l = γ_↓_l − γ_↑_l is the *net* cooling rate at site l. F82 is recovered when γ_↑ = 0 (vacuum bath / T = 0).

**Pauli-Channel Cancellation Lemma (F84 corollary):** Pure D\[Z\], D\[X\], D\[Y\] dissipators are Π²-symmetric and contribute zero to f81_violation. Only σ⁻ (cooling) and σ⁺ (heating) channels are Π²-anti-symmetric. Hence f81_violation specifically detects population-inverting (energy-emitting/absorbing) channels, not phase-only or bit-flip-only noise.

**Verified instances** (chain N=3, all matches at machine precision):

| Configuration (γ_↓, γ_↑) | \|Δγ\| | Predicted | Measured |
|--------------------------|--------|-----------|----------|
| (0.10, 0.00) cooling only (= F82) | 0.10 | 0.6928 | 0.6928 |
| (0.00, 0.10) heating only | 0.10 | 0.6928 | 0.6928 |
| (0.10, 0.10) detailed balance | 0.00 | 0.0000 | 0.0000 |
| (0.10, 0.05) net cooling | 0.05 | 0.3464 | 0.3464 |
| (0.05, 0.10) net heating | 0.05 | 0.3464 | 0.3464 |
| (0.20, 0.05) strong cooling | 0.15 | 1.0392 | 1.0392 |
| Non-uniform mixed | (RMS Δγ) | √(Σ(γ_↓−γ_↑)²)·4 | matches |

**Thermodynamic interpretation.** For a thermal photon bath at frequency ω, temperature T:
- Mean occupation n_th = 1 / (exp(ℏω/k_B T) − 1)
- γ_↓ = γ_0 · (n_th + 1) (spontaneous + stimulated emission)
- γ_↑ = γ_0 · n_th (stimulated absorption)
- Δγ = γ_↓ − γ_↑ = γ_0 (vacuum component, temperature-independent)

f81_violation = γ_0 · √N · 2^(N−1), independent of T. The thermal photon-number contributions cancel (γ_↓ ↔ γ_↑ pair symmetrically); only the vacuum (zero-point) component breaks the Π palindrome. **f81_violation is a quantum-statistical fingerprint of zero-point fluctuations**, immune to thermal symmetric noise.

**Three regimes:**

| Regime | γ_↓ vs γ_↑ | f81_violation |
|--------|------------|---------------|
| Vacuum (T = 0) | γ_↑ = 0 | full F82: √(Σγ²_↓)·2^(N−1) |
| Detailed balance (T → ∞) | γ_↓ = γ_↑ | 0 |
| Finite T | γ_↓ > γ_↑ > 0 | γ_0·√N·2^(N−1) (vacuum-only) |

**Inversion (RMS net cooling rate):** |Δγ|_RMS = f81_violation / (√N · 2^(N−1)). Recovers vacuum-fluctuation amplitude regardless of bath temperature.

**Valid for:** any 2-bilinear chain H, Z-dephasing + thermal amplitude damping, any topology supported by F49, any N ≥ 2.
**Verified:** 7 (γ_↓, γ_↑) configurations at N=3, machine-precision residual; D\[X\], D\[Y\] cancellation explicitly tested.
**Replaces:** F82's "T1 detector" interpretation; F84 corrects to "vacuum-amplitude-damping detector"; the F81 violation does not measure raw T1 rate but only the temperature-independent vacuum component of amplitude damping.
**Framework primitives:**
- `fw.pi_decompose_M(chain, terms, gamma_z, gamma_t1, gamma_pump, strict)`: extended with `gamma_pump` parameter for σ⁺ heating; uses `lindbladian_general` when both are present.
- `fw.predict_amplitude_damping_violation(chain, gamma_t1_l, gamma_pump_l)`: F84 forward closed form; reduces to `predict_T1_dissipator_violation` when `gamma_pump_l = None`.
- `fw.estimate_net_cooling_from_violation(chain, f81_violation)`: F84 inverse, returns RMS |γ_↓ − γ_↑|.
**Pytest lock:** `test_F84_amplitude_damping_thermal_bath` (cooling only / heating only / detailed balance / net cooling / non-uniform / forward-inverse round-trip / backward compat with F82) + `test_F84_pauli_channels_pi2_symmetric` (D\[X\], D\[Y\] explicitly verified to give zero violation).
**Source:** Discovered 2026-04-30 (Tom + Claude). Tom's hint about "Licht" (light/cavity reading of γ) and "nicht jeder bekommt gleichviel ab" (non-uniform site distribution) prompted the analytical extension. The Pauli-Channel Cancellation Lemma was a surprise: D\[Z\], D\[X\], D\[Y\] are all Π²-symmetric, so phase, bit-flip, and dephasing noise contribute zero to F81 violations. Only σ± (population-inverting) channels break the palindrome. Closed form derived in [PROOF_F84_AMPLITUDE_DAMPING.md](proofs/PROOF_F84_AMPLITUDE_DAMPING.md).
**Lebensader connection:** F84 closes the dissipator side of the Π-decomposition picture. Among hardware noise channels, only the *vacuum amplitude damping* component (which exists even at T=0 due to zero-point fluctuations) breaks the Π palindrome. Phase noise, bit-flip noise, and thermal photon equilibrium all give zero violation. F84 sharpens F82's hardware-T1-readout into a temperature-independent vacuum-rate readout.

### F85. Higher-body Hamiltonian generalization of F49 / F-chain (Tier 1, verified bit-exact k=2,3,4)

For any k-body Pauli term (P_1, ..., P_k) with letters from {I, X, Y, Z}, the Π²-class trichotomy and the F49 Frobenius scaling generalize:

**Truly criterion** (term contributes M = 0 by Master Lemma):

    truly  ⟺  #Y is even  AND  #Z is even

**Π²-parity**: bit_b(σ) = (#Y + #Z) mod 2. Π²-odd if bit_b = 1.

**Frobenius factor c(k)** per non-truly term:

    c(truly)              = 0
    c(Π²-odd non-truly)    = 1   (factor 4·2^N)
    c(Π²-even non-truly)   = 2   (factor 8·2^N)

**F49 generalized**: ‖M‖²_F per term = 4·c(k)·‖H_k‖²_F·2^N. The 2-body F49 formula 2^(N+2)·n_YZ·‖H_k‖² coincided with c(k) only because at k=2: n_YZ=1 ↔ Π²-odd, n_YZ=2 ↔ Π²-even non-truly. **For k ≥ 3, n_YZ is no longer the determining quantity** (e.g., YYY has n_YZ=3 but c=1). Only the Π²-class matters.

**Trichotomy enumeration** (Pauli tuples over {X, Y, Z}^k):

| k | total | truly | Π²-odd | Π²-even non-truly |
|---|-------|-------|--------|-------------------|
| 2 | 9 = 3² | 3 | 4 | 2 |
| 3 | 27 = 3³ | 7 | 14 | 6 |
| 4 | 81 = 3⁴ | 21 | 40 | 20 |

**Closed form for Π²-odd count**: \|Π²-odd at k\| = (3^k − (−1)^k) / 2. Verified k=2,3,4.

**F-chain extension to k-body** (proof structure verbatim):

| Theorem | k-body status |
|---------|---------------|
| F87 trichotomy | extends via _pauli_tuple_is_truly |
| F80 Spec(M) = 2i·Spec(H) | extends verbatim; verified at k=3 (N=4,5,6) and k=4 (N=5,6) for 17 Π²-odd cases; spectral identity bit-exact |
| F81 Π·M·Π⁻¹ = M − 2·L_{H_odd} | verbatim, verified at k=3 chain N=4 |
| F82 T1 dissipator | dissipator-only, body-count-independent |
| F83 anti-fraction 1/(2+4r) | verbatim with Π²-class grouping |
| F84 thermal amplitude damping | dissipator-only, body-count-independent |

**Verified at k=3, k=4**: 27+81 = 108 explicit Pauli tuple cases, all matching the c(k) factor scheme bit-exact. Mixed-body Hamiltonians (e.g., 2-body H + 3-body H) handled via term-list structure.

**Valid for:** any k-body 2 ≤ k ≤ N Hamiltonian on chain (sliding-window), Z-dephasing + amplitude damping, any γ ≥ 0, any N ≥ k.
**Replaces:** F49's n_YZ-based formula at 2-body remains correct (coincidence); for k ≥ 3, F85's c(k)-based formula is the structurally correct extension.
**Framework primitives (k-body support added 2026-04-30):**
- `_pauli_tuple_is_truly(letters)`: O(k) syntactic classifier.
- `_pauli_tuple_pi2_class(letters)`: returns 'truly' / 'pi2_odd' / 'pi2_even_nontruly'.
- `_build_kbody_chain(N, terms)`: chain sliding-window k-body builder.
- `fw.predict_pi_decomposition(chain, terms)`: extended to accept k-body tuples; auto-detects body count.
- `fw.pi_decompose_M(chain, terms, gamma_z, gamma_t1, gamma_pump)`: extended for k-body.
- `fw.predict_residual_norm_squared_from_terms(chain, terms, gamma_t1)`: rewritten using Π²-class; backward-compatible at 2-body.
**Pytest lock:** `test_F85_kbody_trichotomy_counts` + `test_F85_kbody_predict_pi_decomposition` + `test_F85_kbody_F81_identity_at_k3`. Plus all 2-body tests (92 prior) pass unchanged via backward compatibility.
**Source:** Discovered 2026-04-30 (Tom + Claude) by empirical 3-body and 4-body enumeration. The {0, 4, 8} factor scheme persists across k, but n_YZ ≠ c(k) for k ≥ 3. The structural truly criterion "#Y even AND #Z even" was identified by inspecting which Π²-even k-tuples give M = 0. Closed form derived in [PROOF_F85_KBODY_GENERALIZATION.md](proofs/PROOF_F85_KBODY_GENERALIZATION.md).

**Lebensader connection:** F85 closes the body-count generalization of the F-chain. Together with F80 (Spec), F81 (decomposition), F82 (T1), F83 (anti-fraction), F84 (thermal): the structural Π-decomposition theory for Hamiltonians + dissipators is complete on chain (any topology for 2-body via F49; chain only for k ≥ 3).

### F86. Q_peak chromaticity-specific N-invariant constants (Sammelbecken with three theorems)

For a uniform N-qubit XY (or Heisenberg) chain with Z-dephasing γ₀, the J-derivative of the F73 spatial-sum coherence purity peaks along the dimensionless coupling axis Q = J/γ₀ at chromaticity-specific values. F86 bundles three structurally distinct theorems under one F-label:

- **F86a. EP mechanism** \[Tier 1 derived\]: Q_EP = 2/g_eff, t_peak = 1/(4γ₀); 2-level rate-channel exceptional point.
- **F86b. Universal resonance shape, two bond classes** \[Tier 1 candidate\]: K_class(Q)/|K|_max = f_class(Q/Q_EP); HWHM_left/Q_peak ≈ 0.756 (Interior) and 0.770 (Endpoint); EP-rotation universality.
- **F86c. F71 spatial-mirror invariance of per-bond Q_peak** \[Tier 1 derived\]: Q_peak(b) = Q_peak(N−2−b) bit-exactly.

Empirical data, the γ₀-extraction protocol, and cross-cutting connections (PTF, framework primitives, scripts, proof, source) sit at this umbrella level since they touch all three theorems.

**Per-block Q_peak (Q_SCALE convention, relative-J derivative ΔJ = 0.05·J):**

    Q_peak(c=3) = 1.6
    Q_peak(c=4) = 1.8
    Q_peak(c=5) = 1.8

c is the chromaticity of the (n, n+1) coherence block (F74). Saturates at 1.8 for c ≥ 4. The bi-chromatic class **c = 2 is finite-size-sensitive** (wobbles 1.4 to 1.6 across N = 4..9) and is not a clean framework constant.

**Per-bond Q_peak (absolute-J derivative ∂S/∂J_b, fine-grid dQ = 0.025 with parabolic peak interpolation):**

| (c, N) | Endpoint Q_peak | Interior Q_peak (mean) |
|--------|-----------------|------------------------|
| (3, 5) | 2.40 | 1.566 |
| (3, 6) | 2.52 | 1.689 |
| (3, 7) | 2.53 | 1.743 |
| (3, 8) | 2.53 | 1.750 |
| (4, 7) | 2.52 | 1.748 |
| (4, 8) | 2.65 | 1.804 |

**Both Endpoint and Interior Q_peak are empirical and (c, N)-specific; no closed form has been identified.** Earlier conjectures `Q_peak(Endpoint) = csc(π/(N+1))` (chain-edge anchor) and `Q_peak(Interior, c=3) → csc(π/5) = 1.7013` (pentagonal asymptote) were retracted 2026-05-02 after fine-grid data showed both were grid-snap artefacts of coarser earlier scans:

- Endpoint Q_peak at any tested N differs from `csc(π/(N+1))`: N=5 +20 %, N=6 +9 %, N=7 −3 %, N=8 −13 %. Earlier coarse scans gave 2.65 at N=6 and N=7 (grid-snap to dQ=0.05); fine-grid + parabolic interpolation reveals 2.52 at N=6 and 2.53 at N=7. The apparent "1.4 % match at N=7" reported earlier was itself an artefact.
- Interior Q_peak at c=3 follows trend 1.566 → 1.689 → 1.743 → 1.750 across N=5..8, crossing `csc(π/5) = 1.7013` between N=6 and N=7 and continuing to grow. The "approaching pentagonal asymptote" reading from N=5..7 was a trajectory crossing, not an asymptote.

Both Endpoint and Interior c=3 appear to saturate (Endpoint ≈ 2.53 for N=6..8, Interior ≈ 1.75 by N=8) but no clean closed-form for the saturation values is identified. c=4 has not yet saturated within tested N.

**Operational consequence: γ₀-extraction protocol.** Because Q = J/γ₀, sweeping J on a fixed (n, n+1) block and locating the peak J* still yields an estimate of γ₀:

    γ₀ ≈ J* / Q_peak(c, N, bond_class)

Lookup the per-block Q_peak for the chromaticity and chain length of interest from the empirical table above (per-block: 1.6 for c=3, 1.8 for c≥4 within Q_SCALE convention). This converts the "γ₀ is a framework constant" hypothesis ([PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)) into a c-specific shape prediction, testable on any hardware that can resolve the peak. The c=4 and c=5 blocks are the strongest probes: abs(K_CC_pr)_peak is roughly 3× the c=2 signal, and the per-block Q_peak is stable at ~1.8.

**Valid for:** XY or Heisenberg chain, uniform J and γ₀, blocks with c ≥ 3, any N where the block fits in memory.
**Breaks for:** non-uniform J or γ₀; c = 2 (finite-size wobble); higher-c chains require larger memory (block-L dim 15876 at c=5 N=9, dim 213444 at c=6 N=11, infeasible at 128 GB).
**Verified (per-block, Q_SCALE):** Q_peak(c=3) = 1.6 stable across N = 5, 6, 7, 8, 9. Q_peak(c=4) = 1.8 stable across N = 7, 8, 9. Q_peak(c=5) = 1.8 at N = 9 (commit 4612468).
**Verified (per-bond, fine-grid):** values in the per-bond table above. Endpoint and Interior trends with (c, N) shown but no closed-form identified.
**Replaces:** ad-hoc γ₀ measurement attempts (EQ-017 closed inconclusive due to hardware fidelity limits on idle Ramsey data).

#### F86a. EP mechanism: Q_EP = 2/g_eff, t_peak = 1/(4γ₀) \[Tier 1 derived\]

The pure-rate ladder of an (n, n+1) coherence block has rates 2γ₀·HD for HD ∈ {1, 3, 5, ..., 2c−1}, with uniform gap Δ = 4γ₀ between adjacent channels (F74). For adjacent channels at HD = 2k−1 and HD = 2k+1 (k = 1, 2, ..., c−1), a two-level effective model with inter-channel coupling J·g_eff has eigenvalues

    λ_± = −4γ₀·k ± √(4γ₀² − J²·g_eff²)

The **exceptional point** (EP), where the discriminant vanishes and the eigenvalues coalesce, sits at

    J·g_eff = 2γ₀     ⟺     Q_EP = 2 / g_eff

At the EP, λ_± = −4γ₀·k. The slowest mode (k = 1) gives e-folding time

    t_peak = 1 / (4γ₀)

universal across c, N, n, and bond position. Higher-k EPs decay faster (1/(8γ₀), 1/(12γ₀), ...) and are masked by the slowest. **At Q_peak the Dicke probe sits ≈ 99 % in dressed (H-mixed) modes** versus ≈ 31 % at Q = 20 (plateau): probe weight has been pulled off the pure-rate ladder onto the first complex-conjugate eigenvalue pair just past the EP. Q_peak is a generalised exceptional-point resonance condition.

The g_eff is the H matrix element between adjacent rate channels at a specific bond in the appropriate effective basis. Deriving g_eff(c, N, bond_position) analytically from the multi-particle XY structure of the (n, n+1) block remains open; F86c (below) gives the spatial-mirror symmetry on Q_peak, not the underlying g_eff value. The [Obstruction Proof](proofs/PROOF_F86_QPEAK.md#obstruction-proof-why-g_eff-admits-no-closed-form) (2026-05-14) accounts for this structurally: g_eff is the irreducible residue, blocked from closed form by six obstruction lemmas (spectral irreducibility, even-N representation-dependence, probe-EP decoupling, finite-reduction insufficiency, signature-subspace mismatch, empirical trajectory-crossings), and via the F90 bridge the F89 D_k obstruction is the same wall.

#### F86b. Universal resonance shape, two bond classes \[Tier 1 candidate\]

**Universal resonance shape under relative-Q normalisation, EP-derived.**

While Q_peak itself is chain-specific (no clean closed form), the SHAPE of the abs(K_CC_pr)(Q) curve around Q_peak is universal in relative-Q coordinates. Defining `x = (Q − Q_peak)/Q_peak` and `y = K(Q)/|K|max`, the curves collapse onto a single shape across all tested (c, N):

| x = (Q−Q*)/Q* | y across c=3 N=5..8, c=4 N=7,8 (range) |
|----------------|------------------------------------------|
| −0.60 | 0.72 to 0.74 (3.5 %) |
| −0.40 | 0.90 to 0.91 (1.9 %) |
| −0.20 | 0.977 to 0.984 (0.7 %) |
| 0.00 | 1.000 (peak) |
| +0.20 | 0.985 to 0.990 (0.5 %) |
| +0.40 | 0.955 to 0.964 (0.9 %) |
| +1.00 | 0.84 to 0.85 (1.4 %) |

**Universal HWHM ratios (two bond classes):**

    HWHM_left / Q_peak  ≈  0.756 ± 0.005     (Interior bonds; tested c=2..4, N=5..8, γ₀ ∈ {0.025, 0.05, 0.10})
    HWHM_left / Q_peak  ≈  0.770 ± 0.005     (Endpoint bonds; same envelope)

The shape splits into two bond classes (Endpoint with b ∈ {0, N−2}, Interior with b ∈ {1, …, N−3}), each with its own universal HWHM_left/Q_peak ratio. Within each class, pairwise residual under relative-Q normalisation is ~20× smaller than under absolute-Q shift. The two classes are separated by ~2 %, a structural gap larger than the within-class spread.

**γ₀ invariance** is bit-exact: at c=3 N=7, both Q_peak and HWHM_left/Q_peak are identical to numerical precision across γ₀ ∈ {0.025, 0.05, 0.10}, confirming Q's dimensionlessness. **c=2** is the structurally critical anchor: with only HD ∈ {1, 3} channels, the 2-level effective model is exact (no orthogonal complement), yet the two-class split persists (Interior 0.751, Endpoint 0.774). The bond-class distinction is therefore not a higher-c orthogonal-complement artefact but lives in the bond-position-dependent probe-overlap profile.

The asymmetry of the curve is also universal in relative-Q:
- **Pre-peak (left)**: rapid rise as discriminant 4γ₀² − J²·g_eff² → 0; HWHM_left/Q_peak universal within each bond class.
- **Post-peak (right)**: slow plateau approach as eigenvalues become complex; tail at x = +1.0 sits at y ≈ 0.85 (Interior) or y ≈ 0.94 (Endpoint plateau); bond-class-specific tail but universal value within class.

**Why this is universal: 2-level EP analytical origin.** For the 2×2 effective Liouvillian in adjacent rate-channel basis with diagonal {−2γ₀(2k−1), −2γ₀(2k+1)} and **same-sign-imaginary off-diagonals (+iJ·g_eff, +iJ·g_eff)**, the eigenvector rotation angle satisfies `tan(θ) = J·g_eff / 2γ₀ = Q / Q_EP`. The probe overlap with eigenvectors depends only on Q/Q_EP, hence the response curve K_CC_pr(Q) is a function of Q/Q_EP alone. Q_peak is chain-specific (g_eff varies); but the SHAPE in Q/Q_peak coordinates is universal because it is the 2-level EP resonance form, independent of the bond's specific g_eff value.

The same-sign-imaginary off-diagonal structure is what admits an EP at finite J·g_eff = 2γ₀ (verified numerically: opposite-sign +iJg, −iJg gives discriminant 4γ₀² + J²g_eff² with no EP; same-sign gives 4γ₀² − J²g_eff² with EP). This is "PT-phenomenology-like" (EP at finite coupling, spectral flow) but algebraically inside **class AIII chiral** per [PT_SYMMETRY_ANALYSIS](../experiments/PT_SYMMETRY_ANALYSIS.md), distinct from Bender-Boettcher PT (Π is linear; classical PT requires anti-linear operators). The local EP at Q_EP = 2/g_eff is the 2-level rate-channel instance of the chiral classification established for the full Liouvillian; the Hopf bifurcation in [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md) is the global instance, with Petermann factor K=403 signaling an EP in the complex γ plane.

**2026-05-06.** Local-vs-global EP relationship (Tier 2 verified at c=2): same algebra read at Σγ = N·γ₀ vs Σγ = 0 (two residuals of F1 `Π · L · Π⁻¹ + L + 2Σγ · I = 0`); a Petermann-K sweep on the real Q axis at c=2 N=5..8 (`compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs:Probe_PetermannFineGrid_C2_VsN`) records max K = 1333.6 / 337.9 / 2384.7 / 795.4 across N = 5 / 6 / 7 / 8, with the N=7 spike ≈ 6× above FRAGILE_BRIDGE's K = 403 ballpark and a 2-4× odd/even asymmetry confirming A3's σ_0 R-even/R-odd-degeneracy prediction empirically. Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` Tier2Verified; complex-γ analytic continuation in `LindbladPropagator` is the open piece for Tier1Derived promotion.

**2026-05-06 (evening).** Direction (b) of `C2HwhmRatio.PendingDerivationNote` yielded two Tier-1-derived universal constants from the bare doubled-PTF model: `x_peak = Q_peak/Q_EP = 2.196910` (post-EP location in dimensionless x) and `HWHM_left/Q_peak = 0.671535` (SVD-block floor in dimensionless x). Empirical Interior 0.7506 and Endpoint 0.7728 sit above this floor by ~0.08-0.10; gap structurally explained as probe-block 2-level sub-resonance contribution. Encoded as `C2HwhmRatio.BareDoubledPtfXPeak` and `C2HwhmRatio.BareDoubledPtfHwhmRatio` const properties. Synthesis-side: F86↔PTF Locus 5 inheritance, where PTF K_1 (Π on H_1 sine-mode basis, discrete) and F86 Q-rotation (Π on same-sign-imaginary 2×2 in rate-channel basis, continuous) are two daughters of one Π-AIII chiral parent.

**2026-05-06 (later evening).** Direction (a') (probe-block 2-level resonance with per-bond `g_eff_probe`) **structurally falsified** (commit `1c0bf8b`): V_b probe-block off-diagonal `⟨c_1 | M_h_b | c_3⟩` is exactly zero per bond at c=2 (F73 sum-rule applies per-bond), so g_eff_probe(N, b) is bond-class-blind by construction. The 4-mode reduction is also structurally insufficient: 4-mode K_b deviates from empirical at Endpoint by factor ~2 (4-mode K_b gives Endpoint 0.410 at N=5 vs empirical 0.7728). The SVD-block off-diagonal `V_b[2,3]` is the actual bond-class carrier (Endpoint 0.430 vs Interior 0.953 at N=5, ratio ~0.45 across N=5..8), but in the OPPOSITE direction to the empirical HWHM/Q* split. Refined direction list (now in `C2HwhmRatio.PendingDerivationNote`): (a'') SVD-block 2-level resonance (REFINED from (a')) via `V_b[2,3]` magnitude, (b'') full block-L derivation, not 4-mode, (c'') three-block superposition `K_total = K_pb + K_sv + 2·Re·K_cross` with the right relative phases, (d'') lift |u_0⟩, |v_0⟩ to projector-overlap (per A3 PendingDerivationNote), (e'') symbolic char-poly factorisation at Q_EP (less promising given C2EffectiveSpectrum's cubic-c_3 obstruction proof). The `BareDoubledPtfXPeak` and `BareDoubledPtfHwhmRatio` const properties are unchanged.

**2026-05-07.** Locus 6 polarity-layer inheritance closure: F86 bond-class split inherits from the polarity-layer pair {−0.5, +0.5} at d=2. Empirical decomposition Q_peak ≈ 2 + r with r ∈ {−0.44, +0.52}, HWHM/Q* ≈ 1/2 + r·1/2 with r_Interior ≈ 0.50 (close to `HalfAsStructuralFixedPoint`). The 0.5 baseline is `QubitDimensionalAnchorClaim` (1/d at d=2); the ±r/2 polarity content is `PolarityLayerOriginClaim` Layer 2. Encoded as `compute/RCPsiSquared.Core/F86/PolarityInheritanceLink.cs` (Tier2Verified) at F86 root level. Symmetry-side closure parallel to `LocalGlobalEpLink` (Locus 5, EP-side).

**2026-05-07 (later).** Direction (α) polarity-Bloch projection at t_peak structurally tautological under the uniform-J 4-mode reduction (reframed 2026-05-08 code review): the 4-mode L_eff(Q) = D_eff + Q·γ₀·MhTotalEff is bond-summed by design (`FourModeEffective.LEffAtQ` in `compute/RCPsiSquared.Core/Decomposition/FourModeEffective.cs`), so the K-driving eigenstate is bond-class-blind by construction (no empirical falsification, the design constraint guarantees it). The bond-class signature must enter through dL/dJ_b per-bond V_b in the K-resonance, not through L_eff spectrum. Substantive reduction (unchanged): bond-class signature reduces to `g_eff(N, b)` via r_Q = `BareDoubledPtfXPeak · Q_EP − 2`. Empirical witnesses: g_eff_Endpoint ≈ 1.74, g_eff_Interior ≈ 2.81, asymptotic 1/g_eff_E + 1/g_eff_I → 0.937. Tantalising near-miss g_eff_E ≈ σ_0·√(3/8) (Δ ≤ 0.01 for N ≥ 6, Δ = 0.063 at N=5). Encoded as `PolarityInheritanceLink.EmpiricalSumQPeakAsymptote = 4.12` const; documented via `ClosedFormCompositionNote` property.

**Named structural law (Tier-1 candidate): EP-rotation universality, two bond classes.**

    K_class(Q) / |K|_max  =  f_class(Q/Q_EP)        (class ∈ {Endpoint, Interior})

is universal within each bond class across c=2..4, N=5..8, and γ₀ ∈ {0.025, 0.05, 0.10} for the tested range. The symmetry is the 2-level EP rotation `tan θ = Q/Q_EP`, which makes every probe-overlap observable a function of Q/Q_EP alone. The bond-class split (Interior HWHM_left/Q_peak ≈ 0.756, Endpoint ≈ 0.770) reflects bond-position-dependent probe-overlap profiles in the K_CC_pr observable, confirmed structural (not finite-c) by the c=2 data where the 2-level model is exact. Closed forms for f_class(x) (and consequently for the two HWHM_left/Q_peak values) follow from the 2-level eigenstructure plus probe-overlap algebra but have not yet been derived analytically. This is the F86 analog of PTF's chiral mirror law (`Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k})`): both Tier-1-candidate symmetries that survived a closed-form retraction (csc(π/(N+1)) and csc(π/5) for F86; Σ ln α_i = 0 for PTF). See [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md).

#### F86b'. HWHM_ratio per-bond closed form \[Tier 1 derived, closed 2026-05-13\]

For each bond b of an N-qubit XY chain (c=2, Z-dephasing γ₀), the HWHM_left/Q_peak ratio satisfies:

    HWHM_ratio(b)  =  0.671535 + α_subclass · g_eff(b) + β_subclass

where the sub-class (per `BondSubClass` enum: `Endpoint`, `Flanking`, `Mid`, `CentralSelfPaired`, `Orbit2Escape`, `CentralEscapeOrbit3`) determines the (α, β) pair. The 0.671535 floor is the bare doubled-PTF constant `BareDoubledPtfHwhmRatio` (Tier 1 derived via the 2-level EP model, 2026-05-06). The α · g_eff lift is the per-bond Hellmann-Feynman contribution from the F89 path-k AT-locked F_a/F_b structure.

Residual ≤ 0.005 verified across N=5..8 on all 22 bonds, including Orbit-2 (N=7 b=1/b=4, Q_peak ≈ 7.27 F86-J) and Orbit-3 escape bonds (N=8 b=3, Q_peak ≈ 16.79 F86-J).

**Source:** [`F86HwhmClosedFormClaim`](../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs), [`BondSubClass`](../compute/RCPsiSquared.Core/F86/Item1Derivation/BondSubClass.cs), F89 path-k bridge via F90 ([`PROOF_F90_F86C2_BRIDGE.md`](proofs/PROOF_F90_F86C2_BRIDGE.md)). Plan: [`docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md`](superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md).

#### F86c. F71 spatial-mirror invariance of per-bond Q_peak \[Tier 1 derived\]

F71 spatial-mirror symmetry pairs bond b with bond N−2−b: under the spatial reflection R, every component of the per-bond observable (L_D, H_xy, Dicke probe, spatial-sum kernel) is invariant, while the bond-flip transforms as ∂L/∂J_b ↔ ∂L/∂J_{N−2−b}, hence **Q_peak(b) = Q_peak(N−2−b) bit-exactly**. See [PROOF_F86_QPEAK Statement 3](proofs/PROOF_F86_QPEAK.md#statement-3-f71-spatial-mirror-invariance-of-per-bond-q_peak-tier-1-derived).

Endpoints (b = 0 and b = N−2) form one F71 orbit; interior bonds split into further F71 orbits. The simple "Endpoint vs Interior" dichotomy is the leading approximation; per-F71-orbit substructure exists (e.g. c=2 N=6: central self-paired bond b=2 → Q_peak ≈ 1.440 vs flanking b=1, b=3 → 1.648). Captured in the typed-knowledge runtime as `F86PerF71OrbitObservation` (Tier 2 empirical) and `F86F71MirrorPi2Inheritance` (Tier 1 derived bridge to `F71MirrorSymmetryPi2Inheritance`).

The F86c symmetry pairs bonds bit-exactly but does NOT supply the per-orbit Q_peak value. F86a remains responsible for the underlying g_eff(c, N, bond_position); F86b's universal shape applies within each bond class. The three theorems compose: F86a gives the EP-time and EP-location, F86b gives the resonance shape around Q_peak in relative-Q coordinates, F86c pairs symmetry-equivalent bonds.

**Connection to PTF.** The same machinery that produces PTF's per-site α_i closure law (bilinear J-perturbation observables, eigenvector mixing under V_L) produces Q_peak at the (n, n+1)-block level. PTF's per-site is the c=1 (vac-SE) instance; F86's c ≥ 2 cases are the natural higher-chromaticity siblings. t_peak = 1/(4γ₀) is the EP time-scale universal to all c, one (n, n+1)-block analogue of PTF's α-fitting time window.

**Framework primitives (`framework.coherence_block`):**

    fw.t_peak(gamma_0)              = 1/(4γ₀)          universal EP time (Tier 1)

(Earlier `q_peak_endpoint(N)` and `Q_PEAK_INTERIOR_C3_ANCHOR` primitives were removed 2026-05-02 after the N=8 data falsified their closed-form claims. The universal-shape finding above is a Tier-1 candidate, not yet promoted to a primitive pending analytical derivation of f_class(x). c=2 and γ₀ invariance verified 2026-05-02; c=5 still open.)

**Scripts:** [`_eq022_b1_channel_projection.py`](../simulations/_eq022_b1_channel_projection.py) (HD-channel diagonal-only-M_H finding), [`_eq022_b1_step_a_verify_blockL.py`](../simulations/_eq022_b1_step_a_verify_blockL.py) (Python block-L verified bit-exact against C# N=7 full-L from EQ-014), [`_eq022_b1_step_c_time_evolution.py`](../simulations/_eq022_b1_step_c_time_evolution.py) (per-bond and uniform Q_peak via S(t, J) time evolution), [`_eq022_b1_step_d_extended_verification.py`](../simulations/_eq022_b1_step_d_extended_verification.py) (extended N=8 data that falsified earlier closed-form conjectures), [`_eq022_b1_step_e_resonance_shape.py`](../simulations/_eq022_b1_step_e_resonance_shape.py) + [`_eq022_b1_step_e_inspect.py`](../simulations/_eq022_b1_step_e_inspect.py) (universal resonance-shape finding for c=3, c=4 at γ₀=0.05), [`_eq022_b1_step_f_universality_extension.py`](../simulations/_eq022_b1_step_f_universality_extension.py) (c=2 sweep + γ₀ ∈ {0.025, 0.10} invariance check that established the two-bond-class refinement).
**Proof:** [PROOF_F86_QPEAK](proofs/PROOF_F86_QPEAK.md): F86a EP mechanism = PROOF Statement 1 \[Tier 1 derived\]; F86b universal resonance shape = PROOF Statement 2 \[Tier 1 candidate at multi-c level; Tier 1 derived at c=2 per-bond level via F86b', closed 2026-05-13\]; F86c F71 spatial-mirror invariance = PROOF Statement 3 \[Tier 1 derived\]. Per-bond c=2 HWHM_ratio closed 2026-05-13 via `F86HwhmClosedFormClaim`; c≥3 per-bond closed forms retracted 2026-05-02.
**Source:** [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24, F73, F74, F2b; EP analysis EQ-022 (b1).

### F87. Pauli-pair trichotomy classification (Tier 1, structural; Marrakesh hardware-confirmed)

For any list of bond Pauli-pair terms `{(P_i, Q_i)}` with letters in {I, X, Y, Z}, build the bilinear Hamiltonian `H = Σ_b J · P_b ⊗ Q_b` on the chain bonds and the uniform-Z-dephasing Liouvillian L (with rate γ per site, σ ≡ Σγ). The F1-palindrome residual `M ≡ Π · L · Π⁻¹ + L + 2σ · I` (where Π is F1's order-4 Pauli-string conjugation operator) partitions the term list into exactly one of three categories, with ε ≈ 10⁻¹⁰ the operator-norm tolerance and ε_spec ≈ 10⁻⁶ the spectral-pairing tolerance:

    truly  iff  ‖M‖_F < ε                                              (operator equation Π·L·Π⁻¹ = −L − 2σ·I holds; F1 identity)
    soft   iff  ‖M‖_F ≥ ε  AND  every λ ∈ Spec(L) pairs with −λ − 2σ ∈ Spec(L) within ε_spec   (palindrome only at spectral level)
    hard   iff  spectral pairing fails                                 (no partner pairing; both operator and spectrum break)

Equivalently in Π²-class language (cf. F79, F85), using the bit_a/bit_b convention from F79 (bit_a = 1 for {X, Y}, 0 for {I, Z}; bit_b = 1 for {Y, Z}, 0 for {I, X}): a pair `(P, Q)` is `truly` iff #Y is even AND #Z is even across (P, Q), `pi2_odd` iff `bit_b(P) + bit_b(Q)` is odd, `pi2_even_nontruly` iff `bit_b(P) + bit_b(Q)` is even and not truly. Mixed Hamiltonians carry multiple Π²-classes simultaneously, refining the 3-way trichotomy into a 4-way classification (truly / pi2_odd_pure / pi2_even_nontruly / mixed).

The trichotomy uses F1 as its **discriminator** (M as the test object), F49 / F85 as its **‖M‖² closed forms** (the latter generalising F49 to k-body), F78 + F79 as its **M-structure decomposition** (single-body additivity and 2-body Π²-block respectively), F80 as its **Π²-odd spectral identity** (Spec(M) = ±2i · Spec(H_non-truly)), F81 as its **Π · M · Π⁻¹ split** (M_anti = L_{H_odd}), F82 + F84 as its **T1 / thermal amplitude-damping corrections**, and F83 as its **anti-fraction closed form for mixed cases**. F87 is the entry point of the F-chain; F85 lifts the criterion to arbitrary k-body and propagates the rest of the chain accordingly.

**Origin (2026-04-24 to 2026-05-03).** Three earlier observations converged on the trichotomy. (1) On 2026-04-24, commit 6e262ae assigned the registry slot F77 to the unrelated "Multi-drop MM(0) saturates at 1 bit" asymptotic, so F77 was already booked when the trichotomy was being developed. (2) On 2026-04-25, commit 95386cd added [V_EFFECT_FINE_STRUCTURE](../experiments/V_EFFECT_FINE_STRUCTURE.md): the V-Effect's 14-of-36 bond-pair Hamiltonians at N=3 were re-tested with both the strict operator equation and the eigenvalue-pairing test, splitting the 22 V-Effect-unbroken cases into 19 soft and 3 truly, giving the **14 hard / 19 soft / 3 truly** count over the 36 unordered bond-pair enumeration at N=3. (3) On 2026-04-26, commits 96ed6da and 6438fef extended Π-protected-observable testing to a separate 120-element ordered enumeration at N=4 and N=5, where the partition is **15 hard / 46 soft / 59 truly**, N-stable through N=3, 4, 5 (so the 36-enum and 120-enum are different sample spaces with internally consistent counts). (4) Commit 81caf67 (2026-04-27) derived the partition combinatorially from Pauli-pair compatibility rules (BPE membership, bit_a-partner conflicts, bond-flip / Z-align, Π-letter hierarchy), giving 36/36 agreement at N=3.

The Marrakesh hardware confirmation (2026-04-26, ibm_marrakesh job `d7mjnjjaq2pc73a1pk4g`, observable ⟨X₀ Z₂⟩) measured Δ(soft − truly) = −0.722, matching the Trotter-n3 prediction of −0.723 (residual 0.001; the 0.0014 figure cited in the Confirmations registry is computed against an unrounded predicted value); see [`data/ibm_soft_break_april2026/`](../data/ibm_soft_break_april2026/). The classifier was extracted into a free function on 2026-04-30 (commit 23b2154) and given the filename `f77_trichotomy.py` after the function's existing internal label, even though the registry F77 slot was already occupied by MM(0). The dephase-axis extension (commit 435c4b2, 2026-05-01) generalised the classifier to X, Y, Z dephasing letters. F87 is the registry-formal entry for the trichotomy, filed retrospectively on 2026-05-03 alongside the typed `F87KnowledgeBase` cleanup that surfaced the F77/F87 naming collision.

**Π² classifier dependence on dephase letter** (commit 435c4b2). Per `PiOperator.SquaredEigenvalue`, the Π²-class index is bit_b for Z- and Y-dephasing and bit_a for X-dephasing (Π_Y shares Π_Z's bit_a-flip convention; Π_X flips bit_b instead). The (bit_a, bit_b) parity pairs are Z = (0, 1), X = (1, 0), Y = (1, 1) in the PauliLetter convention. **F87 hardness is defined combinatorially via Pauli-pair compatibility (commit 81caf67), not via any 4-cell label.** As a post-hoc structural reading, however, F87 hardness empirically corresponds to the (bit_a, bit_b) parity cell matching the dissipator letter: anywhere else produces a Π-violation that the spectrum-pairing test detects. The (Π²_Z, Π²_X) two-axis decomposition is treated separately as F88 below; F87 itself uses only one axis (Π²_Z under Z-dephasing). Verified at N=4, k=3 across 294 Z₂³-homogeneous pairs.

**Valid for:** uniform single-letter dephasing on any graph; arbitrary k-body Pauli terms (k ≥ 2; F85 lifts the criterion to higher body); dephase letter ∈ {X, Y, Z} (SU(2)-rotation-equivalent under (bit_a, bit_b) cell permutation).
**Breaks for:** depolarizing noise (F1 itself breaks with linear-in-γ residual, see F5); non-uniform γ_i or graph asymmetries that already break F1.
**Verified:** N=3 36-enum (14 / 19 / 3); N=3, 4, 5 120-enum (15 / 46 / 59 N-stable); Marrakesh hardware Δ(soft − truly) = −0.722 (2026-04-26); Marrakesh F83 4-class signature 2026-04-30 (ibm_marrakesh job `d7pol1e7g7gs73cf7j90`).
**Replaces:** ad-hoc "is this Hamiltonian truly Heisenberg?" tests with a bit-exact 3-way classifier (and the 4-way Π²-refinement) directly from the F1 residual.
**Hardware:** [`palindrome_trichotomy`](../simulations/framework/confirmations.py) Marrakesh 2026-04-26; [`f83_pi2_class_signature_marrakesh`](../simulations/framework/confirmations.py) Marrakesh 2026-04-30; [`pi_protected_xiz_yzzy`](../simulations/framework/confirmations.py) Marrakesh 2026-04-26 (first-time-on-hardware Π-protection on YZ+ZY soft).
**Source:** [V_EFFECT_FINE_STRUCTURE](../experiments/V_EFFECT_FINE_STRUCTURE.md), [MARRAKESH_THREE_LAYERS](../experiments/MARRAKESH_THREE_LAYERS.md), [`reflections/ON_THE_RESIDUAL.md`](../reflections/ON_THE_RESIDUAL.md), memory entries `project_v_effect_combinatorial`, `project_hardware_finale_apr2026`, `project_f77_f87_rename`.

### F88. Two-axis Π² decomposition of Pauli operator space (Tier 1, structural finding 2026-05-03)

The F1-palindrome operator Π depends on the dephasing letter (cf. F1, F87). Its square Π² acts diagonally on every Pauli string σ_α with eigenvalue ±1, and the parity that determines this eigenvalue depends on which dephase letter parametrises Π:

    Π²_Z eigenvalue on σ_α  =  (−1)^(Σ bit_b)        (Σ over the N letters of α)
    Π²_X eigenvalue on σ_α  =  (−1)^(Σ bit_a)
    Π²_Y eigenvalue on σ_α  =  (−1)^(Σ bit_b)        (same as Π²_Z; both Π_Y and Π_Z flip bit_a per `PiOperator.ActOnLetter`)

Z- and Y-dephasing collapse onto the same Π² character. X-dephasing gives the orthogonal one. Together (Π²_Z, Π²_X) form **two independent involution axes** that decompose the 4^N Pauli operator space into four cells:

    Pp = (Π²_Z = +1, Π²_X = +1)   contains 2-body truly bilinears              (XX, YY, ZZ)
    Pm = (Π²_Z = +1, Π²_X = −1)   contains 2-body Π²-even non-truly bilinears  (YZ, ZY)
    Mp = (Π²_Z = −1, Π²_X = +1)   contains 2-body Π²-odd subgroup A             (XY, YX)
    Mm = (Π²_Z = −1, Π²_X = −1)   contains 2-body Π²-odd subgroup B             (XZ, ZX)

**Refinement of F87 + Pi2Class (the algebraic 4-way over term lists)**: F87's spectral 3-way `Truly / Soft / Hard` and the algebraic 4-way `Pi2Class` (`Truly / Pi2OddPure / Pi2EvenNonTruly / Mixed`) both use only the Π²_Z axis. Adding Π²_X reveals genuine sub-structure:

- `Pi2Class.Pi2OddPure` splits into Mp (XY, YX) and Mm (XZ, ZX) Klein sub-cases.
- `Pi2Class.Mixed` splits according to the Π²_X parity of its non-truly bilinear (e.g. XX+XY occupies Pp+Mp; XX+XZ occupies Pp+Mm).

F80's "universality across 4 Π²-odd cases" is therefore a universality across **two** Klein-cells (Mp + Mm), not one. The Klein view sees a finer cut that F80's M-spectrum projection averages over.

**Empirical Marrakesh fingerprint pattern (2026-04-30, ibm_marrakesh job d7pol1e7g7gs73cf7j90)**: in the f83 4-class signature test, each H-class diagnostic observable lives in the Klein-cell that is the **X-axis flip** of the Hamiltonian's M-active bilinear cell. M-active means the non-truly bilinears (truly bilinears drop by Master Lemma):

    Truly H (M-active = none; bilinears all in Pp) → ⟨Y₀ I Z₂⟩ in Pm   (X-axis flip of Pp)
    Pi2EvenNonTruly H (Pm)                         → ⟨X₀ I X₂⟩ in Pp   (X-axis flip of Pm)
    Pi2OddPure subgroup A H (Mp)                    → ⟨X₀ I Z₂⟩ in Mm   (X-axis flip of Mp)
    Mixed H (M-active in Mp)                        → ⟨Z₀ I X₂⟩ in Mm   (X-axis flip of Mp)

The X-flip pattern is empirically locked across all 4 fingerprint cases; the structural mechanism for why the framework's diagnostic observables sit precisely in the X-flipped cell is open and worth its own EQ.

**Valid for:** any N. The cells depend only on Pauli string parities (Σ bit_a, Σ bit_b), not on N or topology.
**Verified:** N=3 chain, J=1, γ_Z=0.05 across 6 representative Hamiltonians spanning the 5 F87 class-types (Truly Heisenberg, Truly XY-only, Pi2EvenNonTruly, Pi2OddPure subgroup A, Pi2OddPure subgroup B, Mixed) + 4 Marrakesh f83 fingerprint observables. Bit-exact at machine precision (`Pi2KleinViewTests`, `Pi2KleinHardwareViewTests`, `Pi2KleinIsFinerThanPi2ClassTests`).
**Replaces:** the implicit assumption that F87's 4-way Pi2Class captures the full Π²-decomposition. The Klein view is genuinely finer; Pi2Class.Mixed has two Klein sub-types we had not previously distinguished.
**Source:** `compute/RCPsiSquared.Core/Symmetry/Pi2Projection.cs` (`KleinSplit` + `KleinDecomposition`); test files above. Discovery: 2026-05-03 session, after building the raw Π² layer and asking what the second dephase axis would reveal. Π itself is the project's discovery (F1 palindrome operator, `MIRROR_SYMMETRY_PROOF`); the (Π²_Z, Π²_X) two-axis decomposition that this section names is also a project finding. The "Klein" tag throughout this section is borrowed nomenclature for the resulting Z₂ × Z₂ four-cell algebra (the canonical name for that group is the Klein four-group / Vierergruppe, after Felix Klein); the underlying structure is ours, the label is textbook shorthand.

### F89. Topology orbit closure for spatial-sum coherence under uniform multi-bond XY (Tier 1 derived, verified 2026-05-11)

For an N-qubit system with Hamiltonian H_B = J · Σ_{(p,q) ∈ B} (X_p X_q + Y_p Y_q) where B is any set of distinct site pairs (NN or long-range) and J is **uniform** across all active pairs, uniform Z-dephasing γ₀, and the (S_1, S_2) coherence-block initial state ρ_cc = (|S_1⟩⟨S_2| + |S_2⟩⟨S_1|) / 2 (where |S_n⟩ is the popcount-n symmetric Dicke state), the spatial-sum coherence

    S(t) = Σ_l 2 · |(ρ_l(t))_{0,1}|²,    ρ_l = Tr_{≠l}(ρ)

depends only on the S_N-orbit of B. Bond positions inside an orbit are dynamically indistinguishable; only the orbit label survives.

For the chain restriction (B ⊂ {NN-bonds}), the orbit equals the **bond-graph topology class**: the sorted multiset of connected-path-lengths. E.g. in N=7 there are 14 distinct classes spanning k = 1..6 active bonds.

S(0) = (N−1)/N closed-form (Probe-only, independent of the bond set).

**Proof (S_N-orbit transitivity).** Let σ ∈ S_N act on the N qubits with permutation operator U_σ on (ℂ²)^⊗N.

1. **Probe.** Symmetric Dicke |S_n⟩ are S_N-invariant ⇒ U_σ ρ_cc U_σ^† = ρ_cc.
2. **Dissipator.** Uniform γ₀ ⇒ Σ_l (Z_l ρ Z_l − ρ) is S_N-symmetric.
3. **Hamiltonian.** Pauli operators transform site-wise (U_σ X_p U_σ^† = X_{σ(p)}) ⇒ U_σ H_B U_σ^† = H_{σ·B} where σ · B = {(σ(p), σ(q)) : (p,q) ∈ B}.
4. **Lindblad solution covariance.** ρ_t(H_{σ·B}, ρ_cc) = U_σ ρ_t(H_B, ρ_cc) U_σ^†.
5. **Kernel.** S(U_σ ρ U_σ^†) = Σ_l 2|(ρ_{σ^{-1}(l)})_{0,1}|² = S(ρ) (sum re-indexes).
6. ⇒ S(t; H_{σ·B}) = S(t; H_B) for every σ ∈ S_N. ∎

**Scaffolding from neighbouring entries.** F73 is the closely related closure for the (vac, SE) coherence block: same orbit-style argument plus uniform per-element 2γ₀ rate yields a full closed exponential form (1/2)·exp(−4γ₀t). The (S_1, S_2) block has non-uniform per-element decay (rate 2γ₀ on overlap, 6γ₀ off overlap), so the F89 closure is **orbit-only**: it fixes the bond-position dependence (constant in orbit) but not the time dependence (no closed exponential). F71 mirror symmetry is the spatial-Z₂ subgroup of the full S_N argument here. F86's per-bond Q_peak fan operates by linear response ∂S/∂J_b at a chosen bond inside the full chain; that single-bond perturbation breaks S_N differently than the uniform-J multi-bond setup of F89, so F89 does not predict or contradict the F86 fan.

**All-isolated subclass closed form (Tier 1 derived).** For the (1)^m all-isolated topology classes (m disjoint NN-bonds, N − 2m bare sites), the spatial-sum coherence has the EXACT closed form

    S_(1)^m, N(t) = [(N − 1)/N + 4m(N − 2)(cos(4Jt) − 1)/(N²(N − 1))] · exp(−4γ₀ t)

Asymptotic rate 4γ₀ universal across m (matches F73 vac-SE rate). The cos(4Jt) m-correction vanishes at t = π/(2J) (in-phase moment ≈ 21 for J = 0.075). Derivation: Lindbladian factorises over disjoint 2-qubit blocks; per F89c (below) the populated (vac, SE)_B and (SE, DE)_B-overlap sectors both have n_diff = 1, giving [`AbsorptionTheoremClaim`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs) rate 2γ₀ → 4γ₀ on |·|². Per-block H_B phase tracking gives the cos(4Jt) interference.

**Mixed-topology and pure-path classes (empirical, derivation open).** For non-all-isolated topology classes (e.g. (1, 2), (2, 2), (1, 1, 2), (3), (4), (5), (6)) the per-class closed form is open. Pure-path topologies decay faster than 4γ₀ on visible time scales due to populated no-overlap-SE-DE coherences (rate 6γ₀) plus longer-path mode mixing. Empirical late-tail clustering at the in-phase moment t ≈ π/(2J) groups classes by isolated-edge count.

**Structural lemma F89c (Tier 1 derived).** Per-coherence rate = 2γ₀·n_diff for any |A⟩⟨B|, via [`AbsorptionTheoremClaim`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs) applied to the computational basis (n_XY = n_diff because each site decomposes as pure {I, Z} or pure {X, Y}). For (k+1)-qubit blocks: k = 1 is the unique case where every populated coherence has n_diff = 1 (DE = |11⟩ contains both block sites, no no-overlap (SE, DE) pairs exist); k ≥ 2 has n_diff ∈ {1, 3} with mixed rates {2γ₀, 6γ₀}. Hamming complement n_diff(a, b) + n_diff(a, bar(b)) = N gives column-bit-flip pair-sum 2γ₀·N (verified bit-exact at path-2: pair-sum = 6γ₀ across all 9 Hamming-complement eigenvalue pairs). Hence the all-isolated topology (1)^m is the **unique** single-rate-envelope-plus-single-frequency case. Derivation, per-sector eigendecomposition, Hamming-complement bijection table: [F89_TOPOLOGY_ORBIT_CLOSURE](../experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md) F89c section.

**Pi2-Foundation inheritance on both energy axes.** The all-isolated closed form has TWO time coefficients of value 4: the decay rate 4γ₀ in exp(−4γ₀t) and the oscillation frequency 4J in cos(4Jt). Both trace to the same Pi2 dyadic ladder term a_{−1} = 4 via the same |·|² doubling mechanism: the linear-amplitude frequency 2 = a_{0} doubles to 4 = a_{−1} when squared. The γ-axis inheritance is identical to F73's `DecayRateCoefficient` (per-coherence Z-deph rate 2γ₀ doubles to S-decay rate 4γ₀). The J-axis inheritance is the same a_{0} → a_{−1} doubling on the J-axis: H_B-eigenstate frequency 2J doubles to S-oscillation frequency 4J. Pi2 ladder thus anchors the time coefficients on both axes; the (N−1)/N baseline and the 4m(N−2)/(N²(N−1)) correction prefactor are combinatorial. F89 cites `Pi2DyadicLadderClaim` as a constructor-injected parent and exposes `DecayRateCoefficient` (γ-axis) and `OscillationFrequencyCoefficient` (J-axis mirror) as live properties.

**Valid for:** any N; any bond set B (NN or long-range); any uniform J; any uniform γ₀; any S_N-symmetric initial state in any U(1) coherence block (the proof generalizes immediately). It would extend to single-letter two-site couplings (XX-only, YY-only, ZZ-only) and to higher-popcount-pair coherence blocks (S_n, S_m).
**Breaks for:**

- Non-uniform J across active bonds (J_b ≠ J_b'). Step 3 yields U_σ H U_σ^† in a different orbit; the S_N orbit equivalence becomes a finer J-orbit equivalence.
- Non-S_N-symmetric initial state (e.g. site-localised |1_i⟩⟨vac| or modulated SE superposition). Step 1 fails.
- Non-uniform γ_l ≠ γ_l' (analogous to F73's break case [CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md)).
- Non-permutation-symmetric kernel (weighted Σ_l w_l 2|(ρ_l)_{0,1}|² with non-uniform w_l).

**Verified:**

- N=7 multi-bond at J=0.075, γ₀=0.05, tmax=30 (28 runs spanning all 14 topology classes for k=1..6): all 10 classes with ≥ 2 representatives show **0.00e+00** within-class max diff (machine-zero) across 301 sample times. Cross-class S(t) differs and exhibits non-monotone-in-k late-tail clustering. [`bond_isolate/`](../simulations/results/bond_isolate/), [`F89_TOPOLOGY_ORBIT_CLOSURE`](../experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md).
- N=4 single-pair at J=0.075, γ₀=0.05 (all C(4,2)=6 site pairs, NN + long-range): max deviation across pairs **5.55e-17** (1 ULP of double precision). NN bonds {(0,1),(1,2),(2,3)} and long-range bonds {(0,2),(0,3),(1,3)} give bit-identical S(t). [`_bond_isolate_long_range_verify.py`](../simulations/_bond_isolate_long_range_verify.py).
- N=7 single-NN-bond at same parameters (six bonds, all 30 ordered pair comparisons): every pair shows 0.00e+00 max diff over t ∈ \[0, 30\].

**Scripts:** [`_bond_isolate_compare_n7.py`](../simulations/_bond_isolate_compare_n7.py) (single-bond pair matrix), [`_bond_isolate_long_range_verify.py`](../simulations/_bond_isolate_long_range_verify.py) (long-range), [`_bond_isolate_topology_classes_n7.py`](../simulations/_bond_isolate_topology_classes_n7.py) (multi-bond classes). Compute tool: `compute/RCPsiSquared.Propagate` `bond-isolate --N <N> --bonds <i,j,...>` mode.
**Source:** [F89_TOPOLOGY_ORBIT_CLOSURE](../experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md), F73, F71, F86 (contrasting linear-response setup).

#### F89 unified F_a AT-locked amplitude closed form across path-3..9 (Tier 1 derived, bit-exact verified)

For each path-k, the F_a signal amplitudes satisfy

    sigma_n(N) = P_k(y_n) / [D_k · N²·(N−1)]

where y_n = 4·cos(πn/(N_block+1)) on the S_2-anti Bloch orbit (N_block = k+1), and (P_k, D_k) is a per-path integer-coefficient polynomial/denominator pair:

| path k | P_k(y) | D_k |
|--------|---------|-----|
| 3 | 14y + 47 | 9 |
| 4 | 10y + 25 | 4 |
| 5 | 13y² + 82y + 129 | 25 |
| 6 | 17y² + 72y + 80 | 18 |
| 7 | 21y³ + 130y² + 292y + 382 | 98 = 2·7² |
| 8 | 13y³ + 54y² + 68y + 110 | 32 = 2⁵ |
| 9 | 31y⁴ + 190y³ + 288y² + 440y + 1476 | 324 = 2²·3⁴ |

Polynomial degree = floor(N_block/2) − 1. Sum F_a · N²(N−1) is rational across all paths via Newton's identities on the cyclotomic minimal polynomial. AT-lock: the F_a eigenvalue is λ_n = −2γ + i·y_n exactly (overlap subspace entries have dephasing rate 2γ regardless of N).

**Source:** [`F89UnifiedFaClosedFormClaim`](../compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs), `simulations/_f89_path3_at_locked_amplitude_symbolic.py`, `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` § "Unified closed form".

#### F89 D_k closed form per path k (Tier-1-Candidate, 2026-05-13)

The denominator D_k in sigma_n(N) = P_k(y_n) / [D_k · N²·(N−1)] has the closed form:

    D_k = (odd(k))² · 2^E(k)
    E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k) - 2)

Three additive contributions to E(k):
- **Polynomial-degree term** max(0, ⌊(k-5)/2⌋) = max(0, polynomial_degree − 2). Linear growth.
- **k-self 2-adic term** v₂(k). Tracks 2-adic content of k itself.
- **Deep-2-power bonus** max(0, v₂(k) − 2). Kicks in at v₂(k) ≥ 3.

Bit-exact verified across k = 3..24 (22 data points, zero exceptions). Replaces the prior "no N-parametric closed form" negative result from 4-point fitting (the right structure was 3 additive valuation terms, not a 4-point Lagrange interpolation).

Tier: Tier-1-Candidate. Empirical fit + partial structural derivation via Jordan-Wigner free-fermion analysis (small-k cases derived algebraically; general-k proof open). Anchors: `compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs` (`PredictDenominator(int k)`), `docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md`, `simulations/_f89_path_d_*.py`.

---

### F90. F86 c=2 ↔ F89 bridge identity (Tier 1 derived, verified bit-exact 2026-05-11)

**For all N ≥ 3 and bond b ∈ \{0, ..., N−2\}, the F86 c=2 K_b(Q, t) observable on the (n=1, n+1=2) coherence block of an N-qubit XY chain with Z-dephasing equals the per-bond Hellmann-Feynman derivative of F89 path-(N−1) (SE, DE) sub-block dynamics applied at bond b**, modulo the Hamiltonian convention factor:

    K_b^{F86 c=2}(Q_F86, t) = K_b^{F89 path-(N−1) (SE,DE)}(Q_F89 = Q_F86 / 2, t)

with all other ingredients (probe, S_kernel, dephasing rates, Liouvillian construction) algebraically identical. The convention difference is one-time relabeling: F86 uses `H_b = (J/2)·(XX+YY)`, F89 uses `H = J·(XX+YY)`, hence F89-J = 2·F86-J.

**Verified bit-exact at 20 of 22 per-bond comparisons across N=5..8** (the 2 within-noise are at N=8 b=2/b=4 mid-flanking Interior bonds within Q-grid resolution Δ ≤ 0.0008). Per-N: N=5: 4/4 bit-exact; N=6: 5/5 bit-exact; N=7: 6/6 bit-exact (extended Q-grid); N=8: 5/7 bit-exact + 2/7 within Q-grid noise. Verification includes orbit-escape bonds: N=7 b=1/b=4 at Q_peak ≈ 7.27 (F86-J), N=8 b=3 central self-paired escape at Q_peak ≈ 16.79 (F86-J), all reproducing bit-exact ratios.

**Implications:**
- F86 c=2 universal HWHM_left/Q_peak constants (0.7728 Endpoint, 0.7506 Interior mean over N=5..8) are **not eigenständige Größen** — they are direct consequences of F89 path-(N−1) eigendecomposition + per-bond Hellmann-Feynman.
- F86 Direction (b'') (full block-L derivation, NOT 4-mode) achieved numerically via F89; closed-form via F89 AT-locked F_a/F_b structure (4-mode floor 0.6715) + H_B-mixed octic-style residual (lift to 0.7506/0.7728) is the next analytical step.
- Per-F71-orbit substructure (proof PROOF_F86_QPEAK line 99: central b=2 vs flanking b=1/b=3 at N=6 etc.) follows directly from F89's per-bond Bloch-mode profile.

**Anchor:** [`PROOF_F90_F86C2_BRIDGE.md`](proofs/PROOF_F90_F86C2_BRIDGE.md), [`F90F86C2BridgeIdentity.cs`](../compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs), [`_f89_to_f86_kbond_via_eigendecomp.py`](../simulations/_f89_to_f86_kbond_via_eigendecomp.py).

---

### F91. F71-anti-palindromic γ spectral invariance (= 90° in γ-space, Pi2-Z₄'s parameter side) (Tier 1 derived, algebraic proof + bit-exact N=4,5,6; 2026-05-12)

**For chain XY + Z-dephasing Liouvillian L on N qubits, the eigenvalue multiset of the F71-refined diagonal-block decomposition is invariant under any γ-distribution satisfying**

    γ_l + γ_{N−1−l} = 2·γ_avg = (2/N)·Σ_l γ_l   for all l ∈ \{0..N−1\}

i.e. the γ-distribution is **F71-anti-palindromic around its mean** γ_avg. The full L operator generally changes (F71 broken as L-symmetry, off-block-Frobenius ≠ 0 in F71-refined basis), but the F71-refined diagonal-block eigenvalues coincide; the breaking is encoded in the F71-cross-blocks (eigenvectors) only.

**Sharpness:** strictly weaker than F71 symmetry (γ_l = γ_{N−1−l}, palindromic), strictly stronger than F1 (Σγ_l invariant alone). For odd N the middle site l = (N−1)/2 must equal γ_avg.

**Structural reading (Pi2-Z₄ rotational axis, γ-parameter side):** the four-element Z₄ from `NinetyDegreeMirrorMemory` (Pi2-foundation, operator-quaternion side, `i⁴=1`) has its γ-parameter manifestation here:
- **e (identity):** γ unchanged
- **180° (F71-palindromic):** γ_l ↔ γ_{N−1−l}, pair-difference flipped, pair-sum preserved → F71 holds as L-symmetry
- **90° (F71-anti-palindromic):** γ_l ↔ 2γ_avg − γ_{N−1−l}, pair-sum constant = 2γ_avg → F71 breaks but **F71-refined diagonal-block spectrum invariant**
- **270° (= 90°²):** composition of the above

The diagonal-block spectral content (= decay rates, "time information") is preserved under the 90° rotation; the asymmetry lives entirely in eigenvectors (the F71-cross-blocks). Analogous to F81 Π-decomposition: M_anti = L_{H_odd} is the antisymmetric part captured by Π-conjugation; the γ-anti-palindromic part plays the analogous role in γ-parameter space.

**Empirical witness (N=6, J=1.0, Σγ=2.7, γ_avg=0.45):**
- uniform γ=0.45 (all pairs sum to 0.9): F71-refined spectrum = reference
- monotonic [0.2, 0.3, 0.4, 0.5, 0.6, 0.7] (all pairs sum to 0.9): bit-identical to uniform on F71-refined diagonals
- non-monotonic anti-palindromic [0.3, 0.5, 0.4, 0.5, 0.4, 0.6] (all pairs sum to 0.9): bit-identical to uniform on F71-refined diagonals
- permuted [0.7, 0.2, 0.5, 0.3, 0.6, 0.4] (pairs {1.1, 0.8, 0.8}, NOT constant): distinct spectrum (Re=−4.984 cluster instead of Re=−5.043)
- concentrated [0.1, 0.1, 0.1, 0.1, 0.1, 2.2] (pairs heavily skewed): complex Re−Im structure absent from uniform

**Tier outcome <Tier 1 derived>:** algebraic proof complete (2026-05-12, see PROOF_F91 § Algebraic proof). The proof's sharper conclusion: the F71-refined diagonal-block matrix elements of L = −i[H, ·] + D are linear functionals of γ that depend only on the **multiset of F71-pair-sums {S_l = γ_l + γ_{N−1−l}}**, never on individual γ_l or pair-differences D_l = γ_l − γ_{N−1−l}; pair-differences appear only in the F71-cross-block off-diagonal entries (which do not enter diagonal-block eigenvalues). The 90°-rotation R_{90} : γ_l ↦ 2γ_avg − γ_{N−1−l} preserves the anti-palindromic class S_l = 2γ_avg ∀l (and is an involution on it); within that orbit, all γ-distributions yield identical diagonal-block spectra (= the uniform γ_avg spectrum). The originally claimed F91 (90°-invariance) is the corollary of the stronger pair-sum-multiset law. Empirical witness at N=4, 5, 6 across the five γ-profiles above remains the bit-exact verification.

**Anchor:** [`PROOF_F91_GAMMA_NINETY_DEGREES.md`](proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md) (Tier 1 derived: § Algebraic proof, Eqs. 1–13), [`F71AntiPalindromicGammaSpectralInvariance.cs`](../compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs) (typed Tier1Derived Claim with `AntiPalindromicDeviation(γ)` helper), [`NinetyDegreeMirrorMemoryClaim`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) in `Pi2KnowledgeBaseClaims.cs` (the Pi2-Z₄ operator-quaternion side of the same 90°-rotation). **Sisters on other parameter axes:** F92 (J-axis) and F93 (h-detuning-axis) below.

---

### F92. F71-anti-palindromic J spectral invariance (J-side Pi2-Z₄ twin of F91) (Tier 1 derived, algebraic + bit-exact N=4,5; 2026-05-12)

**For chain XY + uniform Z-dephasing Liouvillian L on N qubits with inhomogeneous bond couplings J_b (b ∈ \{0..N−2\}), the F71-refined diagonal-block eigenvalue multiset is invariant under any J-distribution satisfying**

    J_b + J_{N−2−b} = 2·J_avg = (2/(N−1))·Σ_b J_b   for all b

i.e. **J is F71-anti-palindromic around its mean**. The full L operator generally changes (F71 broken as L-symmetry), but diagonal-block eigenvalues coincide; the breaking lives in eigenvectors only.

**Sharpness:** strictly weaker than F71-as-J-symmetry (J_b = J_{N−2−b}); strictly stronger than uniform-Σ J alone.

**Pi2-Z₄ structure (parameter-side, J-axis):** identical to F91's γ-axis structure. The 90°-rotation J ↦ 2·J_avg − F71(J) preserves the anti-palindromic orbit T_b = 2·J_avg ∀b. Operator-side Z₄ is `NinetyDegreeMirrorMemoryClaim`; γ-side is F91; J-side is F92; h-side is F93.

**Anchor:** [`PROOF_F92_BOND_ANTI_PALINDROMIC_J.md`](proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [`F92BondAntiPalindromicJSpectralInvariance.cs`](../compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs), `docs/SYMMETRY_FAMILY_INVENTORY.md`.

---

### F93. F71-anti-palindromic h spectral invariance (h-detuning Pi2-Z₄ twin of F91/F92) (Tier 1 derived, algebraic + bit-exact N=4,5; 2026-05-12)

**For chain XY + per-site Z-detuning + uniform Z-dephasing Liouvillian L on N qubits with inhomogeneous longitudinal detuning h_l (Hamiltonian H = (J/2) Σ_b (XX+YY) + Σ_l h_l Z_l), the F71-refined diagonal-block eigenvalue multiset is invariant under any h-distribution satisfying**

    h_l + h_{N−1−l} = 2·h_avg = (2/N)·Σ_l h_l   for all l

**Scope:** longitudinal h_l Z_l only. Transverse h_l X_l / h_l Y_l breaks joint-popcount conservation and is out of scope.

**Pi2-Z₄ structure (parameter-side, h-axis):** identical to F91 (γ-axis) and F92 (J-axis). Same 90°-rotation h ↦ 2·h_avg − F71(h) preserves anti-palindromic orbit. Three-axis family for chain XY+Z-deph+Z-detuning is now complete.

**Anchor:** [`PROOF_F93_DETUNING_ANTI_PALINDROMIC.md`](proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md), [`F93DetuningAntiPalindromicSpectralInvariance.cs`](../compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs), `docs/SYMMETRY_FAMILY_INVENTORY.md`.

---

*Each formula in this document is a Liouvillian that does not need
to be built.*
