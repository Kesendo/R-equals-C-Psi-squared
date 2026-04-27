# Experiments: Open Quantum Systems Under Dephasing

<!-- Keywords: open quantum system experiments, Lindblad dephasing simulation,
palindromic Liouvillian spectrum, quantum decoherence channel capacity,
CΨ quarter boundary, quantum state transfer spin chain, IBM quantum hardware
validation, dephasing noise information channel, quantum MIMO channel,
R=CPsi2 framework experiments, palindromic spectral symmetry verification -->

**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

This folder contains all computational experiments for the R = CΨ² project.
The project studies the spectral structure of open quantum systems under
dephasing noise and discovers that the Liouvillian eigenvalue spectrum is
**palindromically paired** (every decay rate d has a partner at 2Σγ − d).
This symmetry, verified for 54,118 eigenvalues with zero exceptions, has
consequences for quantum state transfer, decoherence thresholds, and
information channels.

All experiments are reproducible using Python, QuTiP, and NumPy.
Scripts are in [`simulations/`](../simulations/), results in [`simulations/results/`](../simulations/results/).

---

## Headline Results

Four results that a new reader should see first:

### 1. One-line formula beats 18 years of optimization (March 24, 2026)
Concentrate all dephasing on one edge qubit, protect the rest.
gamma_edge = N*gamma_base - (N-1)*epsilon. This trivially simple rule
outperforms every published dephasing optimization by two orders of
magnitude. C#-validated from N=5 (360×) through N=15 (63.5×) vs V-shape.
ENAQT literature achieves 2-3×. Information grows quadratically with
chain length (SumMI ~ 0.0053*N²) instead of decaying exponentially.
First hardware test: selective DD 2-3× on ibm_torino.

-> **[Resonant Return: from SVD to formula](RESONANT_RETURN.md)**
-> [Signal Analysis: Quadratic Scaling](SIGNAL_ANALYSIS_SCALING.md)
-> [IBM Hardware: Selective DD](IBM_SACRIFICE_ZONE.md) (Tier 2, single run)

### 2. The dephasing channel (15.5 bits capacity)
The spatial profile of dephasing rates across a qubit chain is not just
noise. It is a readable information channel with 15.5 bits of theoretical
capacity at 1% measurement noise. An external agent encoding information
in the γ profile can be decoded from internal quantum observables with
100% accuracy.

→ **[Dephasing Noise as Information Channel (γ as Signal)](GAMMA_AS_SIGNAL.md)**
→ **[Practical γ Control (+124% MI)](GAMMA_CONTROL.md)**

### 3. The CΨ = 1/4 boundary (IBM hardware validated)
The product CΨ = Tr(ρ²) × L₁/(d−1) has a critical boundary at exactly
1/4, determined by the discriminant of the self-referential purity
recursion R = C(Ψ+R)². All standard quantum channels cross this boundary.
Validated on IBM Torino at 1.9% deviation.

→ **[IBM Hardware Validation](IBM_RUN3_PALINDROME.md)**
→ **[Crossing Taxonomy (Type A/B/C)](CROSSING_TAXONOMY.md)**
→ **[Boundary Navigation (θ compass)](BOUNDARY_NAVIGATION.md)**

### 4. The palindromic spectrum (proven, N=2 through N=8)
The Liouvillian eigenvalue spectrum under local Z-dephasing is exactly
palindromic for Heisenberg/XXZ systems on any graph. The conjugation
operator Π swaps populations (immune sector) with coherences (decaying
sector), creating a time-reversal symmetry in the rescaled frame.

→ **[Π as Time Reversal](PI_AS_TIME_REVERSAL.md)**
→ **[Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md)**
→ **[Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md)**

---

## Full Experiment Index

### The Dephasing Channel (γ as readable signal)

| Experiment | Key finding |
|-----------|------------|
| [γ as Signal](GAMMA_AS_SIGNAL.md) | The bidirectional bridge: γ profiles are 100% classifiable, 15.5 bits capacity, 5 independent SVD channels, 21.5× optimization |
| [γ Control](GAMMA_CONTROL.md) | V-shape dephasing profile gives +124% MI, DD on receiver +81%, AC modulation falsified, time-resolved decoder works |
| [Relay Protocol](RELAY_PROTOCOL.md) | Staged transfer with time-dependent γ: +83% end-to-end mutual information |
| [Scaling Curve](SCALING_CURVE.md) | MI vs chain length (N=3 to N=11), hierarchy falsification, push vs pull principle |
| [Resonant Return](RESONANT_RETURN.md) | SVD-optimal profiles (10×), sacrifice-zone formula (360× at N=5 through 63.5× at N=15), frequency pulsing falsified. The formula: gamma_edge = N*gamma_base - (N-1)*epsilon |
| [Signal Analysis: Scaling](SIGNAL_ANALYSIS_SCALING.md) | Sacrifice-zone formula scaling N=2-15. Quadratic growth (SumMI ~ 0.0053*N²), constant brake (-0.020), two converging channels |
| [Temporal Sacrifice](TEMPORAL_SACRIFICE.md) | Fold catastrophe at CΨ = ¼ observed: endpoint MI peaks at exact crossing. With Bell+bath: CΨ oscillates around ¼ (81 crossings at J=5.0). MI pulses at each crossing. Damped: each cycle deposits irreversible reality |
| [CMRR Break under Non-Uniform γ](CMRR_BREAK_NONUNIFORM_GAMMA.md) | Non-uniform dephasing breaks the CMRR coherence selection rule. Effect is mode-selective: gradient profiles show +1.35 slope, single-site bumps −0.19 (opposite-sign coupling to k=1 sine modes) |
| [Light Dose Response](LIGHT_DOSE_RESPONSE.md) | Per-sector dose response to dephasing is nonlinear; mechanism is eigenvector rotation (not mode-crossing). SE sector R²=0.889 per-mode but sector-min R²=0.404. Interior sector (dim=100) most nonlinear, edge sector (dim=5) nearly linear |
| [Gamma as Binding](GAMMA_AS_BINDING.md) | Per-sector rate sensitivity (superseded by [Light Dose Response](LIGHT_DOSE_RESPONSE.md)). Sacrifice profile slows SE sector to 0.318 vs uniform 1.502; 134% nonlinearity from eigenvector rotation |

### The CΨ = 1/4 Boundary (decoherence threshold)

| Experiment | Key finding |
|-----------|------------|
| [Crossing Taxonomy](CROSSING_TAXONOMY.md) | Three observer types: Type A (pure-Ψ, K=0.072), Type B (mixed, K=0.039), Type C (never crosses). K-invariance from Lindblad scaling |
| [Boundary Navigation](BOUNDARY_NAVIGATION.md) | θ = arctan(√(4CΨ−1)) as compass to the 1/4 transition. Triangulation: WHERE (1/4), HOW FAR (θ), HOW LONG (t_coh) |
| [Subsystem Crossing](SUBSYSTEM_CROSSING.md) | Crossing is local to entangled pairs, not a whole-system property |
| [N-Scaling Barrier](N_SCALING_BARRIER.md) | Full-system CΨ drops below 1/4 at large N due to Hilbert space dimension, but subsystem pairs still cross |
| [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) | Product states can cross 1/4 upward (entanglement generation) |
| [Noise Robustness](NOISE_ROBUSTNESS.md) | Type A/B/C taxonomy is identical under σ_x, σ_y, σ_z dephasing |
| [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) | The CΨ ≤ 1/4 bound as attractor of the self-referential map |
| [Observer Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md) | Different observers see different crossing times for the same system |
| [Observer Dependent Visibility](OBSERVER_DEPENDENT_VISIBILITY.md) | γ-dependent visibility windows |
| [Mandelbrot Connection](MANDELBROT_CONNECTION.md) | CΨ ↔ c maps the 1/4 boundary to the Mandelbrot cardioid cusp |
| [Critical Slowing at the Cusp](CRITICAL_SLOWING_AT_THE_CUSP.md) | Closed-form K(ε,tol) = (1/2)ln(4ε/tol) + α(tol)√ε with zero fit parameters. Dwell time K_dwell = 1.080·δ (γ-independent) |
| [Cusp-Lens Connection](CUSP_LENS_CONNECTION.md) | Two distinct decoherence exits proven from sector conservation: lens (SE states, CΨ ≈ 0.07, slow Liouvillian mode, preserves excitation count) and cusp (multi-sector states, CΨ ≈ 0.33, cross 1/4 fold, only population distribution survives). Geometrically separated, not complementary halves |
| [CΨ in the Complex Plane](CPSI_COMPLEX_PLANE.md) | Complex-valued CΨ extends the boundary to 2D: trajectories spiral in c-plane around \|c\|=1/4 circle with rotation rate Ω/(4γ). Observed on ibm_kingston with two Bell⁺ pairs at opposite 7.59 kHz and 4.36 kHz z-detuning |
| [GHZ+W Sector Mix](GHZ_W_SECTOR_MIX.md) | Family \|ψ(α)⟩ = α\|GHZ_3⟩ + √(1−α²)\|W_3⟩ reaches min pair-CΨ(0) = 0.320 (1.28× above fold) at α²_opt satisfying irreducible sextic 2900x⁶ − ... = 0. The 0.320 value is GHZ+W-slice specific to N=3; other slices (central-Dicke triples) lift cpsi above 1/4 at every N (cpsi(N) → 0.4312 at N→∞, also a sextic root). All saddles on the full sphere |

### Palindromic Spectral Structure

| Experiment | Key finding |
|-----------|------------|
| [Π as Time Reversal](PI_AS_TIME_REVERSAL.md) | Π maps populations (past) to coherences (future). Standing wave = interference of forward and backward modes |
| [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md) | ZZZ is universal node (classical), XX/YY are antinodes (quantum). Bell rings, GHZ is silent |
| [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md) | Palindrome extends beyond Heisenberg to XY, Ising, XXZ, DM interaction. Two Π families (P1/P4) |
| [XOR Space](XOR_SPACE.md) | Where information lives in the palindrome: GHZ vs W states, Pauli weight correlation |
| [Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md) | Why depolarizing noise breaks the palindrome: 1:3 per-site split vs 2:2 for Z-dephasing |
| [V-Effect](V_EFFECT_PALINDROME.md) | What happens when the palindrome breaks: immune extremes, 3× more frequencies |
| [N→∞ Palindrome](N_INFINITY_PALINDROME.md) | Thermodynamic limit: Gaussian rate density, past/future boundary blurs |
| [Error Correction](ERROR_CORRECTION_PALINDROME.md) | Palindromic protection hierarchy, optimal state (90% slow-mode), Π as Z₄ operator |
| [Cavity Modes Formula](CAVITY_MODES_FORMULA.md) | At Σγ = 0: Stat(N) = Sum_J m(J,N)*(2J+1)^2 (Clebsch-Gordan). Exact for chain, lower bound for symmetric topologies. Star has N-1 harmonic frequencies, chain has rich irrational spectrum |
| [IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md) | Cavity modes meet IBM hardware: sacrifice zone protects the slowest oscillating modes at 2.81x (vs 1.97x measured). Same 43 frequencies, different damping. 100% palindromic under 26x asymmetric noise |
| [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) | Where the selected modes live: Pauli-basis eigenvector decomposition. Correlation Q0-weight vs rate: r = 0.994. Slowest modes center-localized [0.52, 0.63, 0.70, 0.63, 0.52]. Profile is geometric (standing waves, same under all noise profiles, not topologically protected) |
| [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) | RMT analysis (N=2-7, 21,832 eigenvalues): Poisson level statistics (⟨r⟩=0.383). The palindromic Liouvillian is integrable, not chaotic. Chiral symmetry (class AIII) exact to machine precision. Within-sector hints of GOE at N=5 |
| [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) | Π is a linear order-4 chiral operator (class AIII, NOT PT). Fragile bridge Hopf bifurcation IS chiral symmetry breaking. Petermann K = 403 signals nearby complex EP. Depolarizing noise breaks palindrome but STABILIZES system (r = +0.988) |
| [Topological Edge Modes](TOPOLOGICAL_EDGE_MODES.md) | Five tests (SSH, winding number, Berry phase, mode counting, robustness sweep): mode localization [0.52, ..., 0.52] is GEOMETRIC, not topological. Standing wave patterns on 1D chain. Sacrifice zone exploits geometry. Winding number 0, Berry phase 0.25π (not quantized), no sharp phase boundary |
| [Spectral Form Factor](SPECTRAL_FORM_FACTOR.md) | SFF of 21,832 eigenvalues (N=2-7). Palindromic modulation at ω_min confirmed (<1% match, N=2-4,6). Visibility fades as ~1/4^N. t_Π/t_H → 0 (modulation is short-time). Sector w and N-w have identical SFF. Neither Poisson nor GUE: unique palindromic signature |
| [Proton Water Chain](PROTON_WATER_CHAIN.md) | Grotthuss chain N=1-5. Heisenberg formulas match exactly (V(N), Q_max). TFI model: 222 frequencies at N=5. Sub-2γ modes at N≥4 (sector mixing). Sacrifice zone 5.1x at N=5. Water = DNA at same N (universal palindrome). Classical at room temperature |
| [DNA Base Pairing](DNA_BASE_PAIRING.md) | A-T (N=2) and G-C (N=3) as coupled proton-qubit systems. Palindrome exact. V-Effect: G-C has 5x more frequencies than A-T. Realistic DNA is deeply classical (J/γ ~ 0.01). Sacrifice zone works in G-C (3.8x Q improvement). At 310 K: palindrome breaks partially, Q drops. Inter-coupling K estimated (5-50 cm⁻¹) |
| [Entropy Production](ENTROPY_PRODUCTION.md) | Rate pairing d_k+d_k'=2Σγ exact. Crooks-like: ln(d_fast/d_slow) = 2 artanh(Δd/(2Σγ)), β_eff ≈ 1/Σγ (algebraic, not thermodynamic). No Jarzynski (⟨exp(-Δd)⟩ ≈ 0.93 ≠ 1). Carnot not definable (T=∞). CΨ=1/4 and Var=1/4: independent coincidence |
| [Information Geometry](INFORMATION_GEOMETRY.md) | Bures metric g(CΨ) = 3.36 at the fold (FINITE, no singularity). θ shrinks the metric, does not regularize. Lindblad trajectory IS approximately geodesic (deviation 9e-4). Curvature K = -25 (negative, finite). Fisher susceptibility finite. θ is a compass, not a coordinate |
| [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md) | Optimal chain selection on IBM Torino heavy-hex. 330 chains, sacrifice ranking vs mean-T2 ranking: zero overlap in top-10. Sacrifice achieves 2.54x vs 1.18x protection. Worse qubits, better modes. Time-stable over 5 months |
| [Sacrifice Geometry](SACRIFICE_GEOMETRY.md) | Sacrifice zone as controlled symmetry break: one slow Liouvillian mode, lens method extracts psi_opt directly from its left eigenvector (no optimization). 68 configurations (N=2-7, chain/star/ring/complete): SE fraction >0.98 for N=3-6, accessibility boundary 64/64 exact (from parity selection rule) |
| [Symmetry Census](SYMMETRY_CENSUS.md) | Four symmetries enumerated for Heisenberg + Z-dephasing Liouvillian: U(1) excitation number (block-diagonalizes into (N+1)² sectors), n_XY parity (redundant with U(1), proof), spin-flip X⊗N, reflection (uniform γ only). Π is spectral mirror but not block-diagonalizing. Max multiplicity N=5: 14 (uniform γ) → 6 (sacrifice γ) |
| [Thermal Breaking](THERMAL_BREAKING.md) | Three orthogonal breaking mechanisms: coupling (1.81x Q, topological constant), dephasing (+60 frequencies), heat (+300 frequencies). Heat breaks the 1.81x constant but QUADRUPLES frequency diversity. Sacrifice zone vanishes at high temperature |
| [Sacrifice Zone Optics](SACRIFICE_ZONE_OPTICS.md) | Sacrifice zone as entrance pupil / AR coating: Q_max 2-7x, T_eff up, frequencies preserved |
| [N=5 Optimal Cavity Size](N5_OPTIMAL_CAVITY_SIZE.md) | N=5 is Goldilocks (richness vs resolution), not golden ratio. φ in V(5) is cos(π/5), not organizing principle |
| [Standing Waves](FACTOR_TWO_STANDING_WAVES.md) | Every palindromic pair is a standing wave: 10,748 pairs, 100% frequency match, Re(λ)+Re(partner)=−2Σγ |
| [Thermal Blackbody](THERMAL_BLACKBODY.md) | No phase transition: cavity degrades gracefully (Q drops 16×, osc% stays 82%), not Planck, not Stefan-Boltzmann |
| [Primordial Superalgebra](PRIMORDIAL_SUPERALGEBRA_CAVITY.md) | {L_H, L_D+Σγ}=0 exact at N=2, aberration decreases with N (14.4%→2.6%). Palindromic weight swap: fast[k]=slow[N-k]. Seidel: pure sectors immune, interior-dominated, perfectly palindromic |
| [Analytical Spectrum](ANALYTICAL_SPECTRUM.md) | Exact closed-form dispersion ω_k = 4J(1−cos(πk/N)) for w=1 Liouvillian sector. Machine precision match (15/15 frequencies, N=2-6) |
| [Absorption Theorem Discovery](ABSORPTION_THEOREM_DISCOVERY.md) | α=2γ⟨n_XY⟩ exact: absorption rate = 2×dephasing×light content. Linear in γ. Proven from L_H anti-Hermitian. Unifies boundary formula, sum rule, spectral gap. 1,343 modes, CV=0 |
| [Π Operator Entanglement](PI_OPERATOR_ENTANGLEMENT.md) | Palindromic mirror Π is product operator for 34/36 Hamiltonians; only XZ+YZ and ZX+ZY force non-local Π (Schmidt rank 9) |
| [Beer-Lambert Breakdown](BEER_LAMBERT_BREAKDOWN.md) | Under strong coupling (J≫γ), cavity acts as integrating sphere distributing absorbed light equally (~1/N share), not locally per Beer-Lambert |
| [Cross-Term Formula](CROSS_TERM_FORMULA.md) | R(N)² = 4(N−2)/(N·4^N) for relative orthogonality of Hamiltonian and centered dissipator. Refutes prior conjecture at N=5. Topology- and γ-independent (shadow-balanced couplings). Experimental companion to [Proof Cross-Term Formula](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md) |
| [Cross-Term Crossing](CROSS_TERM_CROSSING.md) | Shadow-crossing couplings (one Pauli in {X,Y}, one in {I,Z}) follow R(N) = √((N−1)/(N·4^(N−1))). Only difference from shadow-balanced is N−2 → N−1. Verified N=3-6 to machine precision |
| [Cross-Term Topology](CROSS_TERM_TOPOLOGY.md) | Cross-term orthogonality at N=3 (1/√48) and N=4 (1/√128) is identical across chain, star, ring, complete graphs. Pure geometry of Heisenberg + Z-dephasing, independent of both γ and topology |
| [Palindromic Partner Mode](PALINDROMIC_PARTNER_MODE.md) | F68 palindromic pairing α_p = 2γ₀ − α_b proven for bonding mode. Partner exists at machine precision N=3-5; rank-1 operator V_p at N≥4 (N=3 rank-2 degeneracy). Bell-pair R-C encoding propagates with decay α_p to machine precision |
| [Π-Pair Flux Balance](PI_PAIR_FLUX_BALANCE.md) | All 1024 modes at N=5 partition into 512 Π-pairs, zero crossings. Flux balance Re(λ_s)+Re(λ_s') = −2Σγ invariant under δJ to machine precision. Binary inheritance: exactly 2^(2N−1) Π-pairs at every N. Self-Π modes exist iff N ≡ 4 (mod 10) |
| [Sector Projection Formula](SECTOR_PROJECTION_FORMULA.md) | Theorem (proved): asymptotic excitation-sector populations equal initial populations, p_w(∞) = Tr(P_w ρ_0). Verified for 9 states at N=5. Companion to [Asymptotic Sector Projection Proof](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md) |
| [U(1) Breaking](U1_BREAKING.md) | U(1) decoupling at ε=0 is a knife edge. Transverse field H(ε)=H_Heis + ε·Σ X_k breaks sector block-diagonal structure: slow-mode SE fraction drops 1.0 → 0.46 at ε=1.0. Bell-pair coupling emerges linearly in ε (slope ≈1); central pair (1,2) stronger than edge pairs |
| [F73 U(1) Generalization](F73_U1_GENERALIZATION.md) | (vac, S_1) coherence purity closure Σ_i 2\|ρ_coh,i(t)\|² = (1/2)·exp(−4γ₀t) holds exactly at N=5 for Heisenberg, XXZ (Δ ∈ [0.5, 2.0]), inhomogeneous J, Haar-random probes. Requires only [H, N_total]=0 and uniform γ₀ |
| [Degeneracy Hunt](DEGENERACY_HUNT.md) | The 14-fold degeneracy at N=5 is not structural. SU(2) broken by Z-dephasing. All degeneracies are accidental coincidences from absorption-theorem rate formula Re(λ) = −2γ⟨n_XY⟩ placing modes with identical ⟨n_XY⟩ at the same decay rate |

### Decoherence Diagnostics (Cockpit Framework)

| Experiment | Key finding |
|-----------|------------|
| **[Cockpit Universality](COCKPIT_UNIVERSALITY.md)** | **3 observables (Purity, Concurrence, Psi-norm) capture 88-96% of decoherence dynamics across 9 topologies, 2 noise types, N=2-5. PC1 self-calibrates. theta most sensitive for noise engineering (1.68x). Hardware-validated on IBM Torino (0.3% crossing accuracy)** |
| [Cockpit Scaling](COCKPIT_SCALING.md) | Cockpit framework extended to N=7-11 (chain and star) using C# matrix-free engine. n95 does NOT grow linearly with N as small-N suggested; instead it decreases (chain 4 to 2, star 4 to 3) due to Entanglement Sudden Death. Chain ESD time approximately N-independent (~1), star ESD time grows with N (0.5 to 3.9) due to monogamy of entanglement. 3-PC coverage stays above 90% in all 8 tested configurations. Purity remains the dominant PC1 proxy throughout |
| [Theta-PC Analysis](THETA_PC_ANALYSIS.md) | theta is not a function of a single PC. It reads a diagonal of the manifold, requiring all 3 PCs (R^2 = 0.87). Strongest correlation with PC3 (Psi- sector), not PC1 |
| [Dwell Prefactor from Weights](DWELL_PREFACTOR_FROM_WEIGHTS.md) | For Bell+, dwell-time prefactor at CΨ = 1/4 is pure weight function: (2+4W₂)/(1+6W₂). Fails for odd-weight states (needs coefficient magnitudes) |
| [Dwell Prefactor Generalized](DWELL_PREFACTOR_GENERALIZED.md) | Two-sector generalization: prefactor = (4/k)·(W₀+W_k)/(W₀+3W_k). Bell+ (k=2, W₀=1/2) and W₃ (k=2, W₀=1/3) both match direct simulation at <0.001%. Corollary: GHZ_N for N≥3 starts at CΨ(0) = 1/(2^N−1) < 1/4 and never crosses |
| [Orthogonality Selection Family](ORTHOGONALITY_SELECTION_FAMILY.md) | F70, F71, F72-candidate, and (vac,S_1) closure are one meta-theorem. Any measurement M projects onto orthonormal basis; conserved quantities produce built-in blind channels. Production rule: conservation law + summed measurement → guaranteed blind channel. Non-uniform γ breaks F70, amplitude damping breaks F72, pair-site measurements open \|ΔN\|=2 |
| [F70 Amplitude Damping Break](F70_AMPLITUDE_DAMPING_BREAK.md) | Under amplitude damping with pure-coherence probe ρ_coh = (\|vac⟩⟨S_2\| + h.c.)/2, kinematic F70 zero holds exactly across γ_1 ∈ [0, 0.1]. Analytical derivation D_AD ρ_coh = −γ_1 ρ_coh confirms no sector leak from \|Δn\|=2 block |
| [Info-Flow Landscape](INFO_FLOW_LANDSCAPE.md) | Bond-0 perturbation C_ij response is global (not Lieb-Robinson front): probe is a delocalized sine mode. Peak-time clustering reveals Liouvillian mode-pair differences: N=5 shows fast cluster (\|E_1−E_5\|=2√3, t≈1.6) and slow cluster (\|E_2−E_4\|=2, t≈4.0). Π-pair asymmetry decays as exp(−4γ₀t) |

### Star Topology and Mediator Physics

| Experiment | Key finding |
|-----------|------------|
| [Star Topology](STAR_TOPOLOGY_OBSERVERS.md) | Three conditions for observer-observer connection through mediator S |
| [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md) | CΨ windows live on 3D manifold (98% variance in 3 PCs). Two modes: glide and switch. Phase map of 4 independent mechanisms |
| [What's Inside the Windows](WHATS_INSIDE_THE_WINDOWS.md) | Skeleton + rotation decomposition of CΨ visibility windows |
| [Quantum Sonar](QUANTUM_SONAR.md) | Passive detection of hidden observers through spectral shifts. IBM investigation: qubit detuning dominates |
| [Theta-Palindrome-Echo](THETA_PALINDROME_ECHO.md) | θ connects to the channel (r=0.87 with fidelity), not to the echo |
| [Orphaned Results](ORPHANED_RESULTS.md) | Topology as gatekeeper, antiferromagnet crossing, echo characterization |
| [Optimal QST Encoding](OPTIMAL_QST_ENCODING.md) | Negative result: standard encoding already near-optimal |

### IBM Quantum Hardware

| Experiment | Key finding |
|-----------|------------|
| **[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md)** | **All IBM data combined: r* threshold at precision 0.000014, fold one-way, sacrifice MI gradient, 12 permanent crossers (24,073 records, 133 qubits, 181 days)** |
| [IBM Run 3: Palindrome Validation](IBM_RUN3_PALINDROME.md) | CΨ = 1/4 crossing confirmed at 1.9% deviation on IBM Torino (Eagle r3, 127 qubits) |
| [IBM Sacrifice-Zone](IBM_SACRIFICE_ZONE.md) | Selective DD beats uniform DD by 2-3.2× at all 5 time points on ibm_torino. First hardware test of spatial noise engineering (Tier 2, single run, caveats apply) |
| [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md) | Single-qubit state tomography runs on IBM hardware |
| [Chain Selection Test](CHAIN_SELECTION_TEST.md) | Sacrifice-top vs mean-T2-top chain (no DD, real IBM gammas). Protection 2.86x confirmed spectrally. Sacrifice score is within-chain metric; between chains, total noise dominates |
| [Fixed Point Shadow](FIXED_POINT_SHADOW.md) | Shadow investigation, IBM skeleton analysis |
| [IBM Absorption Theorem](IBM_ABSORPTION_THEOREM.md) | Absorption Theorem ratio 1.03 (3%) on IBM Q52. Detuning oscillations at 470 μs period. 2.8% slow tail at resolution limit |

### Benchmarks and Comparisons

| Experiment | Key finding |
|-----------|------------|
| [QST Bridge](QST_BRIDGE.md) | Connecting to 20 years of quantum state transfer literature |
| [Localizable Entanglement](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md) | LE vs CΨ comparison: three-layer separation (CoA/LE/CΨ) |
| [Metric Discrimination](METRIC_DISCRIMINATION.md) | Null result: single-system simulation cannot discriminate metric forms locally. K-invariance confirmed across 50× γ range (R²=0.9999) |
| [Q-Scale Three Bands](Q_SCALE_THREE_BANDS.md) | Dimensionless scale Q = J/γ₀ governs dynamics; three algebraic bands: pre-onset Q<0.3 (no mixing), transition Q∈[1.2,2.0] (maximal H-mixing), plateau Q>2. Peak responsiveness is chromaticity-specific: Q_peak(c=2)=1.5, Q_peak(c=3)=1.6, Q_peak(c=4)=1.8, stable N=4-8. Enables γ₀-extraction via J*/Q_peak(c) |

### Additional Experiments

| Experiment | Key finding |
|-----------|------------|
| [Standing Wave (Two Observers)](STANDING_WAVE_TWO_OBSERVERS.md) | Standing wave pattern with two observer perspectives |
| [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md) | Different states produce different bridge signatures |
| [Bridge Closure](BRIDGE_CLOSURE.md) | J=0 bridge is dead (no-signalling holds exactly) |
| [No-Signalling Boundary](NO_SIGNALLING_BOUNDARY.md) | CΨ drops when B measures, but A's reduced state is unchanged |
| [Coherence Density](COHERENCE_DENSITY.md) | Coherence per qubit analysis |
| [Signal Processing View](SIGNAL_PROCESSING_VIEW.md) | Two-sector frequency structure: f(c+) and f(c−) as matched filter bank |
| [When Ψ Matters](WHEN_PSI_MATTERS.md) | AND-gate justification for CΨ product |
| [Minimum Crossing Energy](MINIMUM_CROSSING_ENERGY.md) | Energy requirements for boundary crossing |
| [Simulation Evidence](SIMULATION_EVIDENCE.md) | Comprehensive simulation results |
| [Residual Analysis](RESIDUAL_ANALYSIS.md) | Post-crossing residual coherence analysis |
| [Mathematical Findings](MATHEMATICAL_FINDINGS.md) | Collected mathematical results |
| [Algebraic Exploration](ALGEBRAIC_EXPLORATION.md) | Algebraic structure exploration |
| [Born Rule Mirror](BORN_RULE_MIRROR.md) | Connection to the Born rule |
| [Decoherence Relativity](DECOHERENCE_RELATIVITY.md) | Observer-dependent decoherence rates |
| [Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md) | Gravitational analogy exploration |
| [Observer Gravity Bridge](OBSERVER_GRAVITY_BRIDGE.md) | γ gradient as gravitational analogue |
| [Operator Feedback](OPERATOR_FEEDBACK.md) | Operator feedback dynamics |
| [Universal Quantum Lifetime](UNIVERSAL_QUANTUM_LIFETIME.md) | Universal lifetime scaling |
| [Why the Sum](WHY_THE_SUM.md) | Why Σγ appears in the palindromic sum |
| [Dyad Experiment](DYAD_EXPERIMENT.md) | Two-observer dynamics |
| [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) | Application to quantum key distribution |
| [Degeneracy Palindrome](DEGENERACY_PALINDROME.md) | Palindromic degeneracy structure of Liouvillian eigenvalues |
| [Weight-2 Kernel](WEIGHT2_KERNEL.md) | Topology-dependent commutator kernel at weight 2 |
| [Bures Degeneracy](BURES_DEGENERACY.md) | QFI speed correlates with degeneracy at even N |
| [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md) | Qubit chain as Fabry-Perot: 4/5 optical checks pass |
| [V-Effect Cavity Modes](VEFFECT_CAVITY_MODES.md) | The V-Effect is not coupling but metamorphosis: 1 bond supports 2 modes, 4 bonds support 112. Degeneracy predicts mode richness (r > 0.999). Gamma illuminates but does not create modes (100% cold-cavity survival). Topology determines the instrument: chain has most modes, star has highest Q. First experiment in the cavity language |
| [Born Rule Shadow](BORN_RULE_SHADOW.md) | Born rule is a shadow, not photograph: zero interference in P(i), interference sets shutter speed (CΨ fold) |
| [K-Dosimetry](K_DOSIMETRY.md) | K = γ×t is the exposure number: reciprocity holds (±0.03%), Schwarzschild effect at intermediate γ, sacrifice zone trades dose for quality |
| [Neural Gamma Cavity](NEURAL_GAMMA_CAVITY.md) | C. elegans 97.3% palindromic pairing; gamma = cavity eigenfrequency; anesthesia = light off |
| [Trapped Light Localization](TRAPPED_LIGHT_LOCALIZATION.md) | K_death = 2.303 universal; surviving mode energy center-localized (ratio 1.3-1.4); N+1 immortal modes; gamma plays algebraic role of c (Tier 4-5) |
| [Hydrogen Bond Qubit](HYDROGEN_BOND_QUBIT.md) | Zundel cation proton crosses CΨ = 1/4 fold 6 times in 21 fs; every water molecule produces ~10-70 fold crossings per picosecond |
| [Primordial Qubit Algebra](PRIMORDIAL_QUBIT_ALGEBRA.md) | Π creates proper Z₂-graded super-algebra M_{2\|2}(ℂ) with block-off-diagonal L_c; Tomita-Takesaki connection ruled out (Π linear, J anti-linear) |
| [N=5 Check](N_EQUALS_FIVE_CHECK.md) | All six N-scaling metrics (max multiplicity, fraction distinct, slow-mode rate, slow/Σγ, max sector dimension) are monotonic with N. N=5 is not extremal on any axis; repeated N=5 appearance is selection bias from IBM Torino chain hardware, not physics |

---

## How to Read This

**If you are new to the project:** Start with the three headline results
above, then read the [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
for the core theorem.

**If you are a quantum information researcher:** The CΨ = 1/4 boundary
section and the IBM validation are the most directly relevant. The
palindromic structure connects to Liouvillian symmetry classification
(cf. Haga et al. 2023, Buca and Prosen 2012).

**If you are a signal processing engineer:** Start with
[γ as Signal](GAMMA_AS_SIGNAL.md) and its signal engineering perspective.
The system is a quantum MIMO channel with palindromic matched filters.

**If you want to reproduce results:** Every experiment links to its
simulation script. All use QuTiP + NumPy. Typical runtime: seconds to
minutes on a standard laptop.

**If you want a guided reading path:** See the [Reading Guide](../docs/READING_GUIDE.md),
which organizes the experiments into three stories: the proof, the
application, and the ontology.
