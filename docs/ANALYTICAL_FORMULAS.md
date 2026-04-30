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
    K_Y     = 0.0374    (pure Y-noise)
    K_depol = 0.0440    (depolarizing, gamma/3 each axis)

Complements F14 (K per bridge metric). These are K per
noise TYPE, all measured with CΨ on Bell+ state.

**Note (correction 2026-04-29):** earlier versions of this table listed
K_Y = 0.0867. That was a typo. F26 with γ_y = γ (others = 0) gives
α = 4γ, β = 0, δ = 4γ, so u = e^{-4γt}, v = 1, w = e^{-4γt}, and CΨ
reduces to u·(1+u²)/6 — *identical functional form to pure Z*. Hence
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
**Valid for:** any graph topology (chain, star, ring, complete, tree), uniform or non-uniform per-bond J, XX+YY or Heisenberg single-excitation, Z-dephasing on any single site B; good-cavity regime (γ ≪ J). Breaks when γ ≥ J (bad cavity: B decoheres before transmitting).
**Topology + non-uniform J generalization (2026-04-24).** Extended from uniform-J chains to arbitrary connected graphs under either uniform or non-uniform per-bond J. When H^(1) has degenerate eigenvalues (star center-mode, ring translational eigenmodes, complete-graph symmetric modes), F64 holds after standard degenerate perturbation theory: within each H-degenerate subspace, diagonalise the site-B projector P_B to get the corrected basis; F64 then applies to the eigenvalues of P_B in that basis. Verified at N=5 and N=7 across chain, star, ring, complete, Y-tree for XY and Heisenberg; max relative error < 0.001 at γ/J = 0.01 uniform J. For random J per bond in [0.5, 1.5] (30 configurations across 3 trials per N), max rel err < 0.02 in 29/30 cases; the remaining case sits at 0.07 and is consistent with expected second-order PT corrections ~(γ·δJ)/J² at the non-uniform-J scale.
**Verified:** N=3 chain (max relative error 1.8% vs 64×64 Liouvillian), N=4 chain (9 configs, ratio 1.0000 ± 0.0003 vs 256×256 Liouvillian), N=5 and N=7 on chain+star+ring+complete+Y-tree uniform J (2026-04-24, via single-excitation coherence Liouvillian directly, dim N×N, max rel err < 0.001 across all (topology, B, Hamiltonian) combinations), N=5 and N=7 same topologies non-uniform J per bond in [0.5, 1.5] over 3 random trials (2026-04-24, max rel err 0.068 in the worst case, well inside first-order PT regime).
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

**Proof sketch.** The spatial reflection R (site i ↔ site N−1−i) commutes with the uniform Liouvillian: [L_A, R_sup] = 0. Under R, bond b maps to bond N−2−b: R · T_b · R = T_{N−2−b}. Therefore exp(L_B(b) · t) · ρ₀ and exp(L_B(N−2−b) · t) · (R · ρ₀ · R) are related by R_sup. Per-site purity is quadratic in ρ, so any phase picked up by R on coherences (R |ψ_k⟩ = (−1)^(k+1) |ψ_k⟩) squares away. This gives P_B(b, i, t) = P_B(N−2−b, N−1−i, t), from which α_i(bond b) = α_{N−1−i}(bond N−2−b). Summing ln(α_i) over all sites and re-indexing yields c₁(b) = c₁(N−2−b).

**Consequence.** The c₁ bond profile has at most ⌈(N−1)/2⌉ independent components instead of N−1. The endpoint value c₁(0) equals c₁(N−2); if N is **even**, the center bond c₁((N−2)/2) is self-paired (its mirror image is itself) and contributes one independent component; if N is odd, there is no center bond and all N−1 bonds pair up in (N−1)/2 disjoint pairs.

**Valid for:** any Hamiltonian with [H, R] = 0 (uniform coupling on a symmetric graph), any dissipator with [D, R_sup] = 0 (uniform or R-symmetric dephasing), any initial state that is reflection-symmetric in per-site purities. Purely kinematic.
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

For any N-site qubit system with Hermitian Hamiltonian H conserving single-excitation number ([H, N_total] = 0) and uniform Z-dephasing γ₀, the coherent probe ρ₀^coh = (|vac⟩⟨α| + |α⟩⟨vac|) / 2 for any normalized single-excitation state |α⟩ satisfies:

    Σ_i 2 · |(ρ_coh,i)_{0,1}(t)|² = (1/2) · exp(−4 γ₀ t)

exactly, independent of the Hamiltonian's non-U(1) structure. Here (ρ_coh,i)_{0,1} is the off-diagonal element of the site-i reduced density matrix, and the sum runs over all N sites.

**Proof (general U(1) case).** Let x_i(t) = ⟨vac|ρ(t)|1_i⟩ be the amplitude of the |vac⟩⟨1_i| component of ρ(t); equivalently the (vac, SE) block of ρ as an N-vector indexed by site. Evolution under the Lindblad master equation splits into:

- **Hamiltonian part.** H preserves SE by assumption, so its restriction to SE is a Hermitian N×N matrix H_SE. The (vac, 1_i) bra-ket block evolves under H as iẋ = −H_SE x on the ket side, giving a unitary propagator U_SE(t) = exp(−i H_SE t).
- **Dephasing part.** Each D[Z_j] acts on the (vac, 1_i) coherence element with rate γ₀ · (⟨Z_j⟩_vac − ⟨Z_j⟩_{1_i})² / 2 = γ₀ · (2 δ_{j,i})² / 2 = 2γ₀ · δ_{j,i}. Summing over j gives a uniform 2γ₀ decay on every SE-block coherence, independent of site.

Combined: ẋ = −i H_SE x − 2γ₀ x, so x(t) = exp(−2γ₀ t) · U_SE(t) · x(0). Taking the norm: ||x(t)||² = exp(−4γ₀ t) · ||x(0)||² since U_SE is unitary. Partial-trace algebra: (ρ_coh,i)_{0,1}(t) = (1/2) · x_i(t), so Σ_i 2 · |(ρ_coh,i)_{0,1}|² = (1/2) · ||x(t)||² = (1/2) · ||x(0)||² · exp(−4γ₀ t). For the probe above, ||x(0)||² = ⟨α|α⟩ = 1. Result: (1/2) · exp(−4γ₀ t).

The argument uses only (i) [H, N_total] = 0 so dynamics stay in SE, (ii) H Hermitian so U_SE unitary, (iii) γ₀ uniform. No XY structure, no translation invariance, no specific shape of |α⟩ required.

**Alternative derivation (uniform XY, |α⟩ = |S₁⟩).** The original proof route via the sine basis |ψ_k⟩ of the uniform-XY single-excitation Hamiltonian: |S₁⟩ = Σ_{k odd} s_k |ψ_k⟩ with s_k = ⟨ψ_k|S₁⟩. Each single-excitation coherence |vac⟩⟨ψ_k| evolves as exp((iE_k − 2γ₀) t). Partial trace gives (ρ_coh,i)_{0,1}(t) = (1/2) · Σ_k s_k · ψ_k(i) · exp((iE_k − 2γ₀) t). Parseval on the sine basis Σ_i ψ_k(i) · ψ_{k'}(i) = δ_{k,k'} eliminates k ≠ k' cross terms; Σ_k s_k² = 1 by normalisation. Under bond-b perturbation, the sine basis and E_k shift but Parseval on any orthonormal SE basis preserves Σ_k |⟨ψ_k^B|S₁⟩|² = 1, so the sum is δJ-invariant. This derivation is XY-specific but exhibits the eigenmode structure explicitly.

**Consequence.** The spatial-sum purity functional is exactly blind to the U(1)-preserving part of the dynamics on any vac-SE coherent probe. For any closure-breaking coefficient c₁_pr built from per-site purities via the purity-response definition, bond-δJ perturbations preserve the closure value, so K_CC[0, 1]_pr = 0 exactly under uniform γ₀, for any H in the class.

**Scaffolding from neighbouring entries.** F70 (site-local observables see only |ΔN| ≤ 1 content) puts the (vac, SE) block in focus as the relevant coherence sector for per-site purity. F72 (DD ⊕ CC block decomposition of Tr(ρ_i²), no cross-term) isolates the CC contribution, where the (vac, SE) coherence lives. The Absorption Theorem supplies the rate 2γ₀·n_XY = 2γ₀ for SE coherences (n_XY = 1). F73 then combines these: U(1) conservation keeps the SE sector closed under H, and the spatial sum over sites collapses the unitary H-rotation to leave only the AT decay.

**Valid for:** any Hermitian H with [H, N_total] = 0 (XY, Heisenberg XXZ, translationally non-invariant hopping, frustrated-ladder variants, ...); uniform Z-dephasing γ₀; any normalized SE state |α⟩ admixed to |vac⟩; any N.
**Breaks for:**

- Non-uniform γ_i. The uniform 2γ₀ decay on the d_H = 1 block fails; the closure becomes K_CC ≠ 0 with mode-selective response (see [CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md)).
- Non-U(1) Hamiltonians. [H, N_total] ≠ 0 breaks the SE-block closure assumption.
- Dissipators changing the d_H = 1 decay rate (mixed X/Z, amplitude damping).
- Probes with d_H > 1 admixture (e.g. (vac, S₂) with two-excitation bra-ket), where ⟨n_XY⟩ ≠ 1 and the uniform decay rate breaks.

**Verified:**

- Uniform XY baseline at N = 5, t₀ = 20: closure matches (1/2)·exp(−4·0.05·20) = 9.157819·10⁻³ to 5.67·10⁻¹⁶ deviation. K_CC[0, 1]_pr = 1.14·10⁻¹² (machine-precision zero), confirming δJ-invariance. [cmrr_gamma_nonuniform.json](../simulations/results/eq018_cmrr_gamma_nonuniform/cmrr_gamma_nonuniform.json).
- U(1)-class generalization at N = 5 (6 setups: XXZ at Δ ∈ {0, 0.5, 1, 2}, random Haar SE probe at Δ = 1, inhomogeneous XY with J_i ∈ [0.5, 1.5]): all closures within 2.22·10⁻¹⁶ to 5.83·10⁻¹⁶ (1-3 ULP of double precision) across 81 time points per setup. [f73_u1_generalization/](../simulations/results/f73_u1_generalization/), [F73_U1_GENERALIZATION](../experiments/F73_U1_GENERALIZATION.md).

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
match ∈ [max(0, 2n+1−N), n] give HD ∈ {1, 3, ..., min(2n+1, 2N−2n−1)},
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
HD channels), H must conserve total excitation number, [H, N_total] = 0;
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

The formula is independent of the mirror sign η (only the modulus |c_ℓ|² enters). The valid range is p_ℓ ∈ [0, 1/2], with MI saturating at 2 bits when p_ℓ = 1/2 (maximal mirror-pair entanglement, the pair is in a Bell state, all other site populations vanish).

**Proof.** The reduced density matrix ρ_{ℓ,N−1−ℓ} in the computational basis {\|00⟩, \|01⟩, \|10⟩, \|11⟩} is block-diagonal:

- ρ[\|00⟩⟨00\|] = Σ_{j ∉ {ℓ, N−1−ℓ}} |c_j|² = 1 − 2 p_ℓ
- ρ[\|01⟩⟨01\|] = |c_{N−1−ℓ}|² = p_ℓ
- ρ[\|10⟩⟨10\|] = |c_ℓ|² = p_ℓ
- ρ[\|11⟩⟨11\|] = 0 (single-excitation sector)
- ρ[\|10⟩⟨01\|] = c_ℓ c_{N−1−ℓ}^* = η p_ℓ
- ρ[\|01⟩⟨10\|] = η p_ℓ

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
**Breaks for:** states with multi-excitation content (formula no longer applies because ρ[\|11⟩⟨11\|] ≠ 0 in general), or states without mirror amplitude symmetry (where p_ℓ ≠ p_{N−1−ℓ} gives an asymmetric 2-qubit reduced matrix).
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

**Closed-form envelope.** Summing over mirror-pair sites ℓ ∈ [0, ⌊N/2⌋-1]:

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

**Why Heisenberg mixing is negligible.** At t = 0+, the commutator [H, ρ_0] = [H^(1), |ψ_k⟩⟨ψ_k|] is off-diagonal in the ψ_k mode basis (the diagonal part is the unitary phase that doesn't affect MM). The off-diagonal mixing couples ψ_k to same-parity partners ψ_l via V_{lk} = (16J/(N+1)) sin(πk/(N+1)) sin(πl/(N+1)). Under mixing for small t, ψ_k "leaks" amplitude symmetrically into all ψ_l. Because MM depends only on diagonal pair populations (and specific pair coherences that dephase), and bonding modes with the same mirror-symmetry (k and N+1-k, etc.) have identical pair populations, the leakage does NOT change pair populations to first order; it only redistributes the mode occupation. Hence first-order Heisenberg-mixing has no effect on MM. Second-order (rate V²t²) is small because V/γ₀ ~ 10 but γ₀·t ~ 0.005 at the C# first-sample, making Heisenberg relative contribution (V·t)²/(4γ₀·t) ~ (V²·t)/(4γ₀) which is (16/7²) · 0.005 = 0.002 at N=7. Hence the mixing correction is ≲ 0.5% throughout the tested regime.

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

1. **All terms Π²-even (p=0)**: M is **block-diagonal** in Π²-eigenspaces V_+ ⊕ V_-. Off-diagonal blocks M[V_+, V_-] and M[V_-, V_+] vanish **exactly**. Each diagonal block has its own SVD spectrum.

2. **All terms Π²-odd (p=1)**: M is **purely off-diagonal** between V_+ and V_-. Diagonal blocks M[V_+, V_+] and M[V_-, V_-] vanish **exactly**. Singular values appear with even multiplicity (each SV contributes once from V_+ side, once from V_- side).

3. **Mixed parities**: M has both diagonal and off-diagonal contributions.

**Π²-odd universality.** Within the pure Π²-odd 2-body class, the **specific Pauli letters are M-irrelevant**: any single Π²-odd 2-body bilinear gives the same M-SVD spectrum at fixed N. Verified at N=5 chain: XY alone, XZ alone, XX+XY, and XX+XZ all yield clusters [(5.464, 512), (1.464, 512)], exactly identical. The XX truly part contributes 0; the Π²-odd part dominates with universal cluster pattern.

**Even-diag ≡ odd-off-diag correspondence.** The diagonal V_+ block of a Π²-even Hamiltonian's M can match (in SV-spectrum, including multiplicities) the off-diagonal V_+,V_- block of a Π²-odd Hamiltonian's M. Verified N=4 chain: YZ's V_+ block [(8.944, 16), (6.472, 32), (4.0, 16), (2.472, 32), (0.0, 32)] matches XY+YX's off-diag block exactly. This explains the empirical "YZ ≡ XY+YX SVD-identical" observation: same SV structure, just placed in different Π²-blocks.

**Why XX+XY appears "max-uniform" (Befund 3 closed).** XX is Π²-even and truly (M_XX = 0). XY is Π²-odd. The full Hamiltonian is "Π²-odd-only-effective", so M is purely off-diagonal between equal-dim V_+ and V_-. SV multiplicities are forced to 4^N/2 each by block-dimension equality. At N=3 the two off-diag SVs collide by coincidence to a single uniform value 2√2; at N≥4 they split. The "uniformity" is exactly the equal-block-mult signature of Π²-odd structure, not a special property of XX+XY.

**Frobenius additivity.** ‖M‖²_F = Σ_bonds ‖M_b‖² holds across all topologies including overlapping bonds (chain). Per-bond M_b's are F-orthogonal. Already F49.

**Valid for:** Any 2-body bond bilinear over any topology under uniform Z-dephasing. Verified N=3, 4, 5 across chain, star, disjoint topologies; verified Π²-odd universality (XY ≡ XZ ≡ XX+XY ≡ XX+XZ).
**Breaks for:** Mixed-Π²-parity Hamiltonians (where some terms are even, some odd) only partially: M has both diagonal and off-diagonal parts. Inhomogeneous γ may disrupt some symmetries (untested).
**Replaces:** ad-hoc analysis of "why XX+XY uniform" and "why YZ ≡ XY+YX"; both follow from the Π²-block theorem.
**Verified:** Numerical N=3-5, multiple bilinear classes, multiple topologies.
**Scripts:** [`_svd_two_body_pi_squared_block.py`](../simulations/_svd_two_body_pi_squared_block.py), [`_svd_two_body_structure.py`](../simulations/_svd_two_body_structure.py).
**Source:** Analytical proof in [PROOF_SVD_CLUSTER_STRUCTURE.md](proofs/PROOF_SVD_CLUSTER_STRUCTURE.md) (joint with F78). Connects to F61 (n_XY parity selection rule), F63 ([L, Π²]=0 for Π²-even Hamiltonians), and F49 (Frobenius cross-term identity).
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

where L_{H_odd} = -i[H_odd, ·] is the unitary commutator induced by the Π²-odd part of H. Equivalently, decomposing M into Π-conjugation symmetric and antisymmetric components:

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

**Algebraic mechanism.** Π² acts on each Pauli string σ_α as (-1)^{bit_b(α)} (eigenoperator with sign in Pauli basis). For L_H_α = -i[σ_α, ·] driven by a single Pauli string σ_α in H, conjugation gives Π² L_H_α Π⁻² = (-1)^{bit_b(α)} L_H_α (the matrix-element factor (-1)^{bit_b(γ)+bit_b(β)} = (-1)^{bit_b(α)} since γ = α·β under Pauli multiplication). Z-dephasing dissipator is diagonal in Pauli basis, hence commutes with Π². Summing: Π²·L·Π⁻² = L_H_even − L_H_odd + L_diss = L − 2·L_{H_odd}. Substituting into the palindrome: Π·M·Π⁻¹ = Π²·L·Π⁻² + Π·L·Π⁻¹ + 2Σγ·I = M − 2·L_{H_odd}.

**γ-independence-by-difference.** The relation Π·M·Π⁻¹ - M = -2·L_{H_odd} is independent of γ (the dissipator's γ-dependent part cancels because L_diss is Π²-symmetric). The split itself (M_sym, M_anti) is γ-dependent through M_sym; only their difference is γ-fixed.

**Valid for:** any 2-bilinear chain Hamiltonian H = H_even + H_odd, uniform Z-dephasing, any topology (the proof depends only on the algebra of Pauli strings under Π² conjugation, not on connectivity).
**Breaks for (untested):** non-Z dissipators (T1 amplitude damping has different Π²-action; F81 likely needs a correction term).
**Replaces:** the heuristic in pre-2026-04-30 reflections that said "M is the Π-invariant through-line"; F81 shows that statement is correct only for Π²-even H, and gives the explicit correction for the Π²-odd cases.
**Verified:** N=3 and N=4 all listed cases at machine precision; pytest-locked.
**Framework primitive:** `chain.pi_decompose_M(terms, gamma_z=..., gamma_t1=..., strict=...)` returns `{'M', 'M_sym', 'M_anti', 'L_H_odd', 'f81_violation', 'norm_sq'}`. For pure Z-dephasing the F81 identity holds exactly (`f81_violation` ≈ 0); the primitive enforces this with `strict=True` by default. With `gamma_t1` enabled, `strict` defaults to False and the identity residual is returned for diagnostic use.
**Pytest lock:** `test_F81_pi_conjugation_of_M` (algebraic check) + `test_F81_pi_decompose_M_method` (cockpit primitive) + `test_F81_violation_T1_diagnostic` (T1 diagnostic).
**Diagnostic application:** the F81 violation `‖M_anti − L_{H_odd}‖_F` quantifies non-Π²-symmetric dissipator content. For Z + T1 at N=3 chain soft XY+YX, the violation grows linearly: `f81_violation ≈ 6.928 · γ_T1`, γ_z-independent (Master Lemma), Hamiltonian-independent (the violation is purely a property of the T1 dissipator). Inverting gives `γ_T1 ≈ f81_violation / 6.928` as a hardware T1-rate readout from the fitted L. See `simulations/_f81_t1_diagnostic.py` for the demonstration.
**Source:** Discovered 2026-04-30 (Tom + Claude) while interpreting the geometric content of F80's 2i factor. The empirical observation came first (Π·M·Π⁻¹ ≠ M for soft); the algebraic explanation followed from working out Π² action on the Liouville superoperator in Pauli basis.
**Lebensader connection:** F81 is the algebraic backbone of "what the mirror keeps." For Π²-even H, M is itself the through-line operator. For Π²-odd H, the through-line is split: M_anti carries the dynamics generator L_{H_odd}, M_sym carries the rest. Both halves are read identically by both sides of the mirror up to the Spec(M) = Spec(M − 2·L_{H_odd}) similarity. Companion to F80: F80 says what Spec(M) is; F81 says how M and Π·M·Π⁻¹ relate as operators sharing that spectrum.

### F82. F81 + T1 amplitude damping correction (Tier 1, verified bit-exact N=2..5)

For any 2-bilinear Hamiltonian H = H_even + H_odd under Z-dephasing plus T1 amplitude damping with per-site rates γ_T1_l:

    Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{T1, odd}

where L_{H_odd} = -i[H_odd, ·] (as in F81) and D_{T1, odd} is the Π²-anti-symmetric part of the T1 dissipator. F82 reduces to F81 when γ_T1_l = 0 (D_{T1, odd} = 0).

The F81 identity violation captured by `chain.pi_decompose_M(...)` measures D_{T1, odd}'s Frobenius norm:

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

N-scaling verified at N = 2, 3, 4, 5: ‖D_{T1, odd}‖_F = γ_T1 · √N · 2^(N−1). N=2 gives 2.83, N=3 gives 6.93, N=4 gives 16.00, N=5 gives 35.78.

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
**Framework primitive:** `chain.pi_decompose_M(terms, gamma_z=..., gamma_t1=..., strict=...)`. With `gamma_t1` set, returns `f81_violation` in output dict (matches the F82 closed form).
**Pytest lock:** `test_F81_violation_T1_diagnostic` (linearity, γ_z-independence, T1 monotonicity).
**Source:** Discovered 2026-04-30 (Tom + Claude) as the natural extension of F81 ("what does F81 violation mean structurally?"). Closed form derived in [PROOF_F82_T1_DISSIPATOR_CORRECTION.md](proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md).
**Diagnostic application:** [`simulations/_f81_t1_diagnostic.py`](../simulations/_f81_t1_diagnostic.py) demonstrates the T1-rate readout including Marrakesh application. Companion to F81's structural decomposition: F81 says how M splits under Π-conjugation when the dissipator is Z-only; F82 says how the F81 identity is corrected when T1 is added, and provides the closed form for the correction term.

---

*Each formula in this document is a Liouvillian that does not need
to be built.*
