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
This symmetry, verified for 87,376 eigenvalues with zero exceptions, has
consequences for quantum state transfer, decoherence thresholds, and
information channels.

All experiments are reproducible using Python, NumPy, and SciPy.
Scripts are in [`simulations/`](../simulations/), results in [`simulations/results/`](../simulations/results/).

**Index coverage:** 212 experiment writeups are tracked in this folder; all 212 are
linked below (index swept and completed 2026-08-04; three added 2026-08-15).

---

## Headline Results

Four results that a new reader should see first:

### 1. One-line formula beats 18 years of optimization (March 24, 2026)
Concentrate all dephasing on one edge qubit, protect the rest.
gamma_edge = N*gamma_base - (N-1)*epsilon. This trivially simple rule
outperforms every published dephasing optimization by two orders of
magnitude. C#-validated from N=5 (360×) through N=15 (68×) vs V-shape
(peak created Sum-MI, a transport metric, ε→0 sim ideal; ~2-3× on
hardware). ENAQT literature achieves 2-3×. Information grows quadratically with
chain length (SumMI ~ 0.0053*N²) instead of decaying exponentially.
First hardware test: selective DD 2-3× on ibm_torino.

-> **[Resonant Return: from SVD to formula](RESONANT_RETURN.md)**
-> [Signal Analysis: Quadratic Scaling](SIGNAL_ANALYSIS_SCALING.md)
-> [IBM Hardware: Selective DD](IBM_CONCENTRATOR.md) (Tier 2, single run)

### 2. The dephasing channel (15.5 bits capacity)
The spatial profile of dephasing rates across a qubit chain is not just
noise. It is a readable information channel with 15.5 bits of theoretical
capacity at 1% measurement noise. An external agent encoding information
in the γ profile can be decoded from internal quantum observables with
100% accuracy.

→ **[Dephasing Noise as Information Channel (γ as Signal)](GAMMA_AS_SIGNAL.md)**
→ **[Practical γ Control (the two-lever noise law)](GAMMA_CONTROL.md)**

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
| [γ Control](GAMMA_CONTROL.md) | Two-lever noise law: within a shape less Σγ wins; at fixed Σγ concentration at the centre wins (+46%, the V-shape only +6%; the March +124% was a Σγ confound); AC modulation falsified; the March decoder's resolution figure did not survive its control |
| [Relay Protocol](RELAY_PROTOCOL.md) | Staged transfer with time-dependent γ: +83% end-to-end mutual information |
| [Scaling Curve](SCALING_CURVE.md) | MI vs chain length (N=3 to N=11), hierarchy falsification, push vs pull principle |
| [Resonant Return](RESONANT_RETURN.md) | SVD-optimal profiles (10×), sacrifice-zone formula (360× at N=5 through 68× at N=15, peak Sum-MI transport, ε→0 sim), frequency pulsing falsified. The formula: gamma_edge = N*gamma_base - (N-1)*epsilon |
| [Signal Analysis: Scaling](SIGNAL_ANALYSIS_SCALING.md) | Sacrifice-zone formula scaling N=2-15. Quadratic growth (SumMI ~ 0.0053*N²), constant brake (-0.020), two converging channels |
| [Temporal Sacrifice](TEMPORAL_SACRIFICE.md) | Fold catastrophe at CΨ = ¼ observed: endpoint MI peaks at exact crossing. With Bell+bath: CΨ oscillates around ¼ (81 crossings at J=5.0). MI pulses at each crossing. Damped: each cycle deposits irreversible reality |
| [CMRR Break under Non-Uniform γ](CMRR_BREAK_NONUNIFORM_GAMMA.md) | Non-uniform dephasing breaks the CMRR coherence selection rule. Effect is mode-selective: gradient profiles show +1.35 slope, single-site bumps −0.19 (opposite-sign coupling to k=1 sine modes) |
| [Light Dose Response](LIGHT_DOSE_RESPONSE.md) | Per-sector dose response to dephasing is nonlinear; mechanism is eigenvector rotation (not mode-crossing). SE sector R²=0.889 per-mode but sector-min R²=0.404. Interior sector (dim=100) most nonlinear, edge sector (dim=5) nearly linear |
| [Gamma as Binding](GAMMA_AS_BINDING.md) | Per-sector rate sensitivity (superseded by [Light Dose Response](LIGHT_DOSE_RESPONSE.md)). Sacrifice profile slows SE sector to 0.318 vs uniform 1.502; 134% nonlinearity from eigenvector rotation |

### The CΨ = 1/4 Boundary (decoherence threshold)

| Experiment | Key finding |
|-----------|------------|
| [Crossing Taxonomy](CROSSING_TAXONOMY.md) | Three observer types: Type A (pure-Ψ, K=0.072), Type B (mixed, K=0.039 in the February tool's feedback model, ln(4/3)/8 = 0.03596 in standard Lindblad), Type C (never crosses). K-invariance from Lindblad scaling |
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
| [Fold and Cusp: Two Seams](FOLD_AND_CUSP_TWO_SEAMS.md) | The palindrome fold at Re = −Nγ and the coherence-horizon merge cusp coincide only at N=2 and separate as N grows; the cusp drifts −2γ → −4γ. Tier 1 for the computed geometry, Tier 2 for the reading |
| [Envelope Rise Boundary](ENVELOPE_RISE_BOUNDARY.md) | Where the full-state CΨ envelope starts to rise: certified a pure (N, Q=J/γ) observable, with an N ≥ 4 floor (N=3 never rises even at Q=2000) and a Q_c(N) contour, Q_c(4) ≈ 27, Q_c(5) ≈ 45 |
| [Operator Rigidity Across the Cusp](OPERATOR_RIGIDITY_ACROSS_CUSP.md) | EQ-031: the 120-Pauli-pair F87 trichotomy interrogated across the CΨ = 1/4 cusp and the N=3 → N=4 regime boundary. Tier 1 simulation, 0/120 category shifts; 15 truly / 46 soft / 59 hard fully stable |
| [Coherence Horizon: the EP-Sensor Debate](COHERENCE_HORIZON_EP_SENSOR_DEBATE.md) | The coherence-horizon exceptional point Q*(N) taken to the EP-sensor debate: genuinely defective by four routes at N=2-5, Petermann factor diverges in kind only, peak height a grid artifact. Geometry Tier 1, debate contribution Tier 2, no metrological verdict |
| [R=CΨ² as Decoherence Readout](RCPSI_DECOHERENCE_READOUT.md) | The Born-rule deviation C = R/Ψ² read as a γ-meter grounded in F94. Tier 2/3: leading order Tier-1 proven, invertibility gate-verified numerically rather than analytically, and convention-dependent |

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
| [IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md) | Cavity modes meet IBM hardware: sacrifice zone protects the slowest oscillating modes at 2.80x (vs 1.97x measured). Same 43 frequencies, different damping. Palindromic to the eigensolver's floor under 26x asymmetric noise |
| [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) | Where the selected modes live: Pauli-basis eigenvector decomposition. Correlation Q0-weight vs rate: r = 0.994. Slowest modes center-localized [0.52, 0.63, 0.70, 0.63, 0.52]. Profile is geometric (standing waves, same under all noise profiles, not topologically protected) |
| [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) | RMT analysis (N=2-7, 21,840 eigenvalues): Poisson level statistics (⟨r⟩=0.383). The palindromic Liouvillian is integrable, not chaotic. Chiral symmetry (class AIII) exact to machine precision. An early within-band GOE hint at N=5 was driven to a verdict and is small-sample noise |
| [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) | Π is a linear order-4 chiral operator (class AIII, NOT PT). Fragile bridge Hopf bifurcation IS chiral symmetry breaking. Petermann K = 403 signals nearby complex EP. Depolarizing noise breaks palindrome but STABILIZES system (r = +0.988) |
| [Topological Edge Modes](TOPOLOGICAL_EDGE_MODES.md) | Five tests (SSH, winding number, Berry phase, mode counting, robustness sweep): mode localization [0.52, ..., 0.52] is GEOMETRIC, not topological. Standing wave patterns on 1D chain. Sacrifice zone exploits geometry. Winding number 0, Berry phase 0.25π (not quantized), no sharp phase boundary |
| [Spectral Form Factor](SPECTRAL_FORM_FACTOR.md) | SFF of 21,840 eigenvalues (N=2-7). Palindromic modulation at ω_min confirmed (<1% match, N=2-4,6). Visibility fades as ~1/4^N. t_Π/t_H → 0 (modulation is short-time). The light-content bands at w and N-w have identical SFF. Neither Poisson nor GUE: unique palindromic signature |
| [Proton Water Chain](../docs/water/PROTON_WATER_CHAIN.md) | Grotthuss chain N=1-5. Heisenberg formulas match exactly (V(N), Q_max). TFI model: 222 frequencies at N=5. Sub-2γ modes at N≥4 (sector mixing). Sacrifice zone 5.1x at N=5. Water = DNA at same N (universal palindrome). The runs sit at Q < 1 by an input whose stated source does not carry it; where water actually sits is a band, 0.04 ≲ Q ≲ 4.6 |
| [DNA Base Pairing](DNA_BASE_PAIRING.md) | A-T (N=2) and G-C (N=3) as coupled proton-qubit systems. Palindrome exact. V-Effect: G-C has 5x more frequencies than A-T. At the parameters used, DNA is deeply classical (J/γ ~ 0.01, an unsourced denominator). Sacrifice zone works in G-C (3.8x Q improvement). At 310 K: Q drops and the palindrome breaks, the latter because the warm channel set shares an axis with the dephasing, not because of temperature. Inter-coupling K estimated (5-50 cm⁻¹) |
| [Entropy Production](ENTROPY_PRODUCTION.md) | Rate pairing d_k+d_k'=2Σγ exact. Crooks-like: ln(d_fast/d_slow) = 2 artanh(Δd/(2Σγ)), β_eff ≈ 1/Σγ (algebraic, not thermodynamic). No Jarzynski (⟨exp(-Δd)⟩ ≈ 0.93 ≠ 1). Carnot not definable (T=∞). CΨ=1/4 and Var=1/4: independent coincidence |
| [Information Geometry](INFORMATION_GEOMETRY.md) | Bures metric g(CΨ) = 3.36 at the fold (FINITE, no singularity). θ shrinks the metric, does not regularize. Lindblad trajectory IS approximately geodesic (deviation 9e-4). Curvature K = -25 (negative, finite). Fisher susceptibility finite. θ is a compass, not a coordinate |
| [Concentrator Qubit Mapping](CONCENTRATOR_MAPPING.md) | Optimal chain selection on IBM Torino heavy-hex. 330 chains, sacrifice ranking vs mean-T2 ranking: zero overlap in top-10. Sacrifice achieves 2.53x vs 1.18x protection. Worse qubits, better modes. Time-stable over 5 months |
| [Concentrator Geometry](CONCENTRATOR_GEOMETRY.md) | The concentrator as a controlled symmetry break: one slow Liouvillian mode, lens method extracts psi_opt directly from its left eigenvector (no optimization). 69 configurations (N=2-7, chain/star/ring/complete): SE fraction >0.98 for N=3-6, accessibility boundary exact with zero violations (from parity selection rule) |
| [Symmetry Census](SYMMETRY_CENSUS.md) | Four symmetries enumerated for Heisenberg + Z-dephasing Liouvillian: U(1) excitation number (block-diagonalizes into (N+1)² sectors), n_XY parity (redundant with U(1), proof), spin-flip X⊗N, reflection (uniform γ only). Π is spectral mirror but not block-diagonalizing. Max multiplicity N=5: 14 (uniform γ) → 6 (sacrifice γ) |
| [Thermal Breaking](THERMAL_BREAKING.md) | Three orthogonal breaking mechanisms: coupling (1.81x Q, a geometric constant and explicitly NOT a topological one), dephasing (+60 frequencies), heat (+300 frequencies). Heat breaks the 1.81x constant but QUADRUPLES frequency diversity. Sacrifice zone vanishes at high temperature |
| [Concentrator Optics](CONCENTRATOR_OPTICS.md) | The concentrator as entrance pupil / AR coating: Q_max 2-7x, T_eff up. The analogy holds on the absorption side only; the concentrator also moves the resonance frequencies, the ω=0 level most |
| [N=5 Optimal Cavity Size](N5_OPTIMAL_CAVITY_SIZE.md) | N=5 is Goldilocks (richness vs resolution), not golden ratio. φ in V(5) is cos(π/5), not organizing principle |
| [Standing Waves](FACTOR_TWO_STANDING_WAVES.md) | The spectrum is completely paired at every N: 21,840 eigenvalues across N=2..7, 9,921 pairs of distinct partners with Re(λ)+Re(partner)=−2Σγ plus 1,998 self-paired on the fixed locus |
| [Thermal Blackbody](THERMAL_BLACKBODY.md) | No phase transition: cavity degrades gracefully (Q drops 16×, osc% stays 82%), not Planck, not Stefan-Boltzmann |
| [Primordial Superalgebra](PRIMORDIAL_SUPERALGEBRA_CAVITY.md) | {L_H, L_D+Σγ}=0 exact at N=2, aberration decreases with N (14.4%→2.6%). Palindromic weight swap: fast[k]=slow[N-k]. Seidel: pure sectors immune, interior-dominated, perfectly palindromic |
| [Analytical Spectrum](ANALYTICAL_SPECTRUM.md) | Exact closed-form dispersion ω_k = 4J(1−cos(πk/N)) for the (0,1) coherence block. Machine precision match (15/15 frequencies, N=2-6) |
| [Absorption Theorem Discovery](ABSORPTION_THEOREM_DISCOVERY.md) | α=2γ⟨n_XY⟩ exact: absorption rate = 2×dephasing×light content. Linear in γ. Proven from L_H anti-Hermitian. Gives the boundary formula and sum rule a common reading; relocates the spectral gap. 1,342 modes, CV=0 |
| [Π Operator Entanglement](PI_OPERATOR_ENTANGLEMENT.md) | Palindromic mirror Π is a product (local) operator for every palindromic Hamiltonian; XZ+YZ and ZX+ZY just need a continuous per-site rotation, not a discrete one (corrected 2026-06-02; once read as non-local rank 9) |
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
| [The Atmosphere-Cluster](ATMOSPHERE_CLUSTER.md) | Tier 2 structural characterization of a 32-mode (N=6) / 24-mode (N=5) intrinsic F1-mirror-pair cluster: intrinsic to uniform γ, ε only shifts \|Im\|. Four findings verified; closed form \|Im\|(N, J, γ) still open |
| [Gamma Fold: a Pair of Mirrors](GAMMA_FOLD_PAIR_OF_MIRRORS.md) | The two involutions on the γ axis, the gain turn and the anti-watch turn, chained by the exact identity L_anti(γ) = L(−γ) − 2σ·Id; composing them translates by 2σ. Adopted into MirrorWorld as `GammaFold` |
| [Majorana Axis Modes](MAJORANA_AXIS_MODES.md) | Empirically verified at N=4, bit-exact: site-reflection R sorts the 94 axis modes into 58 even and 36 odd, all 18 silent modes are R-even, and Im(λ) decomposes into integer combinations of the golden-ratio Majorana dispersion {±φ, ±1/φ} |
| [Survivor Flip and Reflection-Odd](SURVIVOR_FLIP_AND_REFLECTION_ODD.md) | Gate-first verified at N=4, 6 for chain and ring: the half-filling survivor is X^⊗N-odd, reflection-odd, not fixed by staggered particle-hole, and dark with real λ; one antisymmetric Neumann density wave read three ways |
| [Slow-Mode R-Parity](SLOW_MODE_R_PARITY.md) | Slow-mode spectrum decomposed by site-reflection parity at N=4, 5, 6: the stationary subspace is exclusively R-even (dim N+1), the first slow band R-balanced, F86 lives entirely in R-even, leaving a parallel R-odd channel awaiting a probe |
| [V-Effect Boundary Localization](V_EFFECT_BOUNDARY_LOCALIZATION.md) | Computational (Tier 1-2), residual 10⁻¹⁵ at N=3, 4: the palindrome holds exactly in the extreme XY-weight sectors w = 0 and w = N for every 2-body H; all breaking is confined to 0 < w < N and needs bit_b-parity-violating terms |
| [V-Effect Fine Structure](V_EFFECT_FINE_STRUCTURE.md) | Computational (Tier 1-2): the 36 two-term combos re-run under the strict operator equation give 33 broken vs 3, against the V-Effect's 14 vs 22; the gap is a new 19-case soft-break class the spectrum-pairing test missed |
| [Two-Term Palindrome: Klein Routing](TWO_TERM_PALINDROME_KLEIN_ROUTING.md) | Computed (Tier 2), bit-exact at N=3, 4, 5: the 12 two-term Hamiltonians that break X-parity yet keep the mirror are routed onto their hidden-symmetry Q families by the Klein-Vierergruppe index of their bilinears |
| [The Palindrome Classifier](THE_PALINDROME_CLASSIFIER.md) | First reading of the settled C# palindrome classifier as a tool: it reads Hamiltonian terms, not the 4^N spectrum, so it never meets the N=8 wall. Charts the protected interior and the two coasts (field, frustration) where protection ends |
| [PTF Palindrome-Breaking Perturbations](PTF_PALINDROME_BREAKING_PERTURBATIONS.md) | Tier 2 computed: a single-site transverse field breaks the PTF time-rescaling closure Σ ln α ≈ 0 entirely (up to +5.6, RMSE blows up), but the surprise is that the field does NOT break the palindrome |
| [Z⊗N Partnership](Z_N_PARTNERSHIP.md) | Simulation only (Aer/numpy), hardware sketch open: Z⊗N is a strong symmetry whenever every H term and jump operator carries an even number of X/Y factors, so the X-basis Néel mirror pair diagnoses single transverse fields, which break it |
| [Chain Gap Sector Diagnostic](CHAIN_GAP_SECTOR_DIAGNOSTIC.md) | Tier 1 candidate: the chain dissipation-gap slow mode lives in the central diagonal popcount block, the Absorption Theorem matches it bit-exact at N=4, 5, 6, prefactor ⟨n_XY⟩ ≈ 0.55·Q²/N² to ~1% |
| [XXZ Axis: Band-Edge to Lebensader](XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md) | Tier 2 computational, bit-exact at N=4, 5 (the two-clocks charge/spin reading Tier 3): walking Δ along the XXZ axis, the slowest mode relays from the bright band-edge to the Lebensader at a handover Δ* read by bisection |
| [Filling Threshold Chaos](FILLING_THRESHOLD_CHAOS.md) | Tier 1, live-witnessed (`inspect --root fillcsr`): dissipative quantum chaos in the XXZ dephasing spectrum turns on with filling, the excitation density, not with breaking integrability. The Galois door-C question returned a clean null |
| [τ_max and the Spectral Gap](TAU_MAX_SPECTRAL_GAP.md) | Rejected: the auto-extracted τ_max = ħ/√(λ₂·J²) is wrong twice over (square root, spurious 1/J). The relaxation clock is τ = 1/λ₂ = 1/(2γ), set by γ alone; J sets frequencies, not rates |
| [The Shared Skeleton](THE_SHARED_SKELETON.md) | Six small systems (chain, water, benzene, butadiene, cyclobutadiene, hexatriene) at Q = 1.5 with their Liouvillian spectra overlaid: what is shared is law (F1 palindrome, absorption rungs, the real-axis diagonal); the fine structure is each Hamiltonian's own |
| [The One-Four Thesis](ONE_FOUR_THESIS.md) | Asks whether every period-4 structure in the repo is the same quarter-turn of i. Closed: the strong sameness form is refuted at every N in closed form (deviation ±2^(N−1)); the factorization form V_g = Ad(U_g)∘Π²∘K is derived and exact |

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
| [Star Confocal Limit](STAR_CONFOCAL_LIMIT.md) | Tier 1 derived: the star's Casimir gap gives Im_max = J·N/2 for every N ≥ 3 and every (J, γ). Saturation max\|Im λ\| = ΔE_max is universal, but the star uniquely minimises ΔE_max among connected graphs (searched N ≤ 6) |
| [The Hub Kills the Horizon](THE_HUB_KILLS_THE_HORIZON.md) | Refuted (Tier 2, numerical) at N=5, 6, 7: the wheel graph builds no coherence-horizon bridge from star to ring at any ε ∈ [0.1, 50]. The hub, not the bandwidth, is decisive; removing it restores a finite Q* |

### IBM Quantum Hardware

| Experiment | Key finding |
|-----------|------------|
| **[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md)** | **All IBM data combined: r* threshold at precision 0.000014, fold one-way, sacrifice MI gradient, 12 permanent crossers (24,073 records, 133 qubits, 181 days)** |
| [IBM Run 3: Palindrome Validation](IBM_RUN3_PALINDROME.md) | CΨ = 1/4 crossing confirmed at 1.9% deviation on IBM Torino (Eagle r3, 127 qubits) |
| [IBM Concentrator](IBM_CONCENTRATOR.md) | Selective DD beats uniform DD by 2-3.2× at all 5 time points on ibm_torino. First hardware test of spatial noise engineering (Tier 2, single run, caveats apply) |
| [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md) | Single-qubit state tomography runs on IBM hardware |
| [Chain Selection Test](CHAIN_SELECTION_TEST.md) | Sacrifice-top vs mean-T2-top chain (no DD, real IBM gammas). Protection 2.86x confirmed spectrally. Sacrifice score is within-chain metric; between chains, total noise dominates |
| [Fixed Point Shadow](FIXED_POINT_SHADOW.md) | Shadow investigation, IBM skeleton analysis |
| [IBM Absorption Theorem](IBM_ABSORPTION_THEOREM.md) | Absorption Theorem ratio 1.03 (3%) on IBM Q52. Detuning oscillations at 470 μs period. 2.8% slow tail at resolution limit |
| [IBM Concentrator Reloaded](IBM_CONCENTRATOR_RELOADED.md) | Site-resolved Absorption-Theorem pricing flown on ibm_kingston 2026-07-11 after a committed pre-registration: A-sign CONFIRMED at ≈ 5.8σ, A-magnitude off-prediction and marginal, L null-consistent |
| [IBM F129 Ramsey Fringe](IBM_F129_RAMSEY_FRINGE.md) | The F129 level collision (1,5,7) ~ (2,4,8) at n = 9 as a standing Ramsey fringe. Pre-registered before the shot, flown and confirmed on ibm_kingston 2026-07-15 |
| [IBM Block-CΨ Saturation](IBM_BLOCK_CPSI_SATURATION.md) | Tier 2, hardware-verified: block-CΨ on ibm_kingston q13-q14 reaches 88.2% of the 1/4 ceiling at t = 0, and the closed-form (1/4)·exp(−4γt) trajectory fits 5 t-points at R² = 0.9977 |
| [IBM K-Partnership Sketch](IBM_K_PARTNERSHIP_SKETCH.md) | K-partnership folding of the F67 receiver menu cross-validated on ibm_marrakesh: partner-pair deviations 15% and 46% Δ/mean on hardware versus 0.02 to 0.25% on Aer, making site and bond asymmetry measurable |
| [IBM Receiver Engineering Sketch](IBM_RECEIVER_ENGINEERING_SKETCH.md) | First receiver-engineering flight on ibm_kingston at N=5: bonding:2 beats the alt-z-bits baseline on MI(0, 4) by 2.80×, so the receiver advantage survives real hardware decoherence |
| [Marrakesh: Three Layers](MARRAKESH_THREE_LAYERS.md) | Tier 1 hardware-grounded synthesis: the ibm_marrakesh soft-break dataset re-read through three nested layers, F87 trichotomy → F80 M-spectrum → 7 observable categories; also shows the dataset's original T1-amplification hardening hypothesis is quantitatively wrong (Trotter n=3 discretization) |
| [Polarity Fingerprint (Tier-B Marrakesh)](POLARITY_FINGERPRINT_TIERB_MARRAKESH.md) | `fw.polarity_fingerprint` applied to 11 hardware-tested Hamiltonian instances (5 distinct term sets) across 3 Tier-B Marrakesh/Kingston datasets: 6/11 F112 BALANCED and 5/11 SILENT (Π²-even H against a homogeneous bath has nothing to balance; read as 11/11 until 2026-08-07), 8/11 in typed scope, so F112's typed scope is sufficient but not necessary |
| [F112 Hardware Lens (Kingston)](F112_HARDWARE_LENS_KINGSTON.md) | F112 polarity diagnostic run on Liouvillians fitted to Kingston hardware trajectories. Status: the model family is wrong for this data, so it identifies no noise channel; the channel question is open and needs a redo |
| [F113 T1 Extraction (Kingston)](F113_T1_EXTRACTION_KINGSTON.md) | First hardware use of the F113 closed form as a diagnostic: inverting it on Kingston f95 data recovers γ_T1 self-consistently to machine precision; the fit-versus-calibration gap (1.14 to 1.47) is read as absorbed non-T1 noise, not separable from epoch drift |
| [F120 Moment Tower (Kingston)](F120_MOMENT_TOWER_KINGSTON.md) | F120's moment tower flown on ibm_kingston with no two-qubit gates: the structural law is confirmed (rung 1 silent, rung 2 fires, m* ≤ 5). The first rate-layer reading of a pump-bound violation was corrected the same day to T1 telegraphing between epochs |
| [F81 Violation: Hardware Bridge](F81_VIOLATION_HARDWARE_BRIDGE.md) | Grounds the f81_violation operator diagnostic into a measurable number: it IS 2^(N−1) times the RMS identity-escape velocity. The missing heating leg was flown twice, closing the attribution and giving V_σ± = 0.02491 ± 0.00015 /µs on Marrakesh |
| [γ₀ Is Always There](GAMMA0_IS_ALWAYS_THERE.md) | Reading the carrier Q = J/γ₀ off its only lever J on ibm_kingston q13-q14: the coherent swing is born exactly at Q = 1, and the first-peak step halves as J doubles |
| [Clock Field: Site-Owned](CLOCK_FIELD_SITE_OWNED.md) | Three pre-registered questions put to four IBM calibration histories: 1/T2 is site-owned on all four (ICC 0.74-0.87), no compensation anywhere, Kingston and Fez common-mode, Torino independent. The PTF label is flagged as a reading |
| [Concentrator A/B Mechanism Test](CONCENTRATOR_AB_MECHANISM_TEST.md) | Negative result plus methods lesson, arc PARKED: three ibm_kingston runs failed to settle concentrator mechanism A vs B; created-MI was 56-96% classical-mixing artifact and measures transport, not protection |
| [Fold Shadow in Existing Hardware](FOLD_SHADOW_IN_EXISTING_HARDWARE.md) | Tier 2 honest record, deliberately NOT a Confirmations entry: a sweep of ~96 existing hardware files finds the F89d fold/price signature present but degraded, broken by the SPAM noise floor, so the pre-registered price-pair run is genuinely needed |

### Hardware Pre-Registrations (predictions frozen before the shot)

| Experiment | Key finding |
|-----------|------------|
| [Record Parity](RECORD_PARITY_HARDWARE_PREDICTION.md) | Pre-registration of the signed record-parity angle law on IBM hardware. v33, SIM-GATE STAGE, NOT FLOWN: design rounds 1 through 31 folded, sim gate built and rehearsed (joint power 0.9834 at the worst corner); binding freeze and Tom's explicit go still ahead |
| [Price Pair](PRICE_PAIR_HARDWARE_PREDICTION.md) | Pre-registered test of the F89d price law Γ(D) + Γ(D̄) = Σγ_j as a device invariant; since flown, four runs on ibm_marrakesh 2026-07-04. The law held wherever noise was local; deviations decoded into drift, quasi-static dephasing and coherent nearest-neighbour ZZ ≈ 4 kHz |
| [Staircase Null-Test](STAIRCASE_NULLTEST_HARDWARE_PREDICTION.md) | Pre-registration (v3, no hardware data): the conditional-Ramsey rate side on ibm_kingston, testing the exact dephasing null Δ(m) = Γ_tot,m(1 − 2p∞,m). Protocol, estimator and σ bands frozen |
| [Chiral Mirror](CHIRAL_MIRROR_HARDWARE_PREDICTION.md) | Tier 1 kinematic prediction from K_1 symmetry: K_1-paired sine modes give ⟨X_i⟩ → +(−1)^i, ⟨Y_i⟩ → −(−1)^i, ⟨Z_i⟩ identical. Verified in simulation to 10⁻¹⁶; hardware protocol proposed on Heron r2, not yet flown |
| [F130 Hardware Infeasibility](F130_HW_INFEASIBILITY.md) | Design-stage null result, no QPU spent: the F130 beat protocol cannot fly on Heron-class hardware. The collision survives Trotterization but the beat amplitude dies, redirecting the hardware question to a first-order observable |

### Benchmarks and Comparisons

| Experiment | Key finding |
|-----------|------------|
| [QST Bridge](QST_BRIDGE.md) | Connecting to 20 years of quantum state transfer literature |
| [Localizable Entanglement](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md) | LE vs CΨ comparison: three-layer separation (CoA/LE/CΨ) |
| [Metric Discrimination](METRIC_DISCRIMINATION.md) | Null result: single-system simulation cannot discriminate metric forms locally. K-invariance confirmed across 50× γ range (R²=0.9999) |
| [Q-Scale Three Bands](Q_SCALE_THREE_BANDS.md) | Dimensionless scale Q = J/γ₀ governs dynamics; three algebraic bands: pre-onset Q<0.3 (no mixing), transition Q∈[1.2,2.0] (maximal H-mixing), plateau Q>2. Peak responsiveness is chromaticity-specific: Q_peak(c=2)=1.5, Q_peak(c=3)=1.6, Q_peak(c=4)=1.8, stable N=4-8. Enables γ₀-extraction via J*/Q_peak(c) |
| [Naked vs Chain Benchmark](NAKED_VS_CHAIN_BENCHMARK.md) | A naked Bell pair under Z-dephasing against F67 bonding-mode chain encoding with dephasing only at the far end: protection 4.0× at N=3, 7.2× at N=4, 12.0× at N=5, at the cost of reduced initial concurrence |
| [Receiver vs γ-Sacrifice](RECEIVER_VS_GAMMA_SACRIFICE.md) | γ-profile sacrifice-zone engineering against receiver choice at uniform γ₀, via C# brecher scans at N=5, 7, 9. Tier 2: receiver choice wins by 11.5× in absolute Peak Sum-MI, 15.4× with moderate J-modulation. Carries a correction note on superseded coarse-grid Python numbers |
| [Dicke vs Endpoint Probe (JW, N=11)](DICKE_VS_ENDPOINT_PROBE_JW_N11.md) | JW algebraic comparison of Dicke against endpoint-localized probes at N=11. The critical finding is negative: the \|C_b\|² formula is F71-invariant but inverted at Center, so it cannot predict the F86 g_eff bond dependence; the signed sum shows interference at Orbit 2, not Center |
| [J-Blind Receiver Classes](J_BLIND_RECEIVER_CLASSES.md) | Tier 2: states whose observables do not respond to J at all, sorted into three mechanisms (DFS plus bond-eigenstate, H-degenerate subspace closed under L_D, M_α-polynomial). The single-mechanism reading was wrong; chromaticity and N-scaling stay open |

### The F89 Family (the octic, its monodromy, its Galois group, its seeds)

| Experiment | Key finding |
|-----------|------------|
| [Path-K Galois](F89_PATH_K_GALOIS.md) | Tier 1 derived: the algebraic certificate Gal(F₈/ℚ(i)(q)) = S₈, with path-3..6 giving S₈/S₁₈/S₃₂/S₅₃, all non-solvable, so the H_B-mixed relaxation rates are unwritable in radicals. Tier 2 for the amplitude q-dependence |
| [Path-K Diabolic](F89_PATH_K_DIABOLIC.md) | Tier 1 for the mechanism, Tier 2 for completeness: free-fermion integrability keeps the diabolic crossings at every N, Δ-verified; their placement on the physical real q axis is parity-gated, odd N ≥ 7 |
| [Monodromy Mirror](F89_MONODROMY_MIRROR.md) | Tier 1 derived: octic monodromy generates S₈ from eigenvalue braids, and the mirror's base-space face q ↦ −q̄ intertwines it; but the Re λ = −4 fibre fold does not commute with the braiding, so the mirror splits at the Galois boundary |
| [Multi-Sector Monodromy](F89_MULTI_SECTOR_MONODROMY.md) | Tier 1 derived multi-sector census: the S₈ braid is confined to the D₄ orbit at N=4 but spreads to a 12-sector joint-popcount diamond at N=5. Exclusion half derived at N=5, census evidence through N=11 |
| [Branch Locus Palindrome](F89_BRANCH_LOCUS_PALINDROME.md) | Tier 1 derived: the path-3 octic's EP and diabolic branch locus mirrors about Re λ = −4, forced by the F1 palindrome carried antiunitarily; verified on the committed octic to 4·10⁻¹³, no orphan |
| [Beta-Exotic Genericity](F89_BETA_EXOTIC_GENERICITY.md) | Asks whether each count-dropping F89 seed is a generic √-type EP2. Reduces to s₆ ≠ 0, settled unconditionally over ℚ at N = 5, 7, 9 for both R-parities; open for all N |
| [Seed Existence Reduction](F89_SEED_EXISTENCE_REDUCTION.md) | Turns the codim-1 corollary's one numerical input into an exact counting identity: a real defective seed is forced at every odd N because a path graph on odd vertices carries a zero mode, modulo a genericity check and the literal-count reading premise |
| [Topology Controls Galois Writability](F89_TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md) | Tier 1 derived: bond-graph topology controls radical-writability. The complete graph caps factor degrees at 4 (writable for all N), the star at a fixed S₉, while ring and chain grow into full S_n. Tier 2 for the exact growth forms |
| [Topology Orbit Closure](F89_TOPOLOGY_ORBIT_CLOSURE.md) | Tier 1 derived: the spatial-sum coherence S(t) depends only on the S_N orbit (topology class) of the bond set, not on bond placement; verified bit-identical at N=7 across all 14 classes and at N=4 across all 6 site pairs |
| [Seed Holonomy: the Third Clock](SEED_HOLONOMY_THIRD_CLOCK.md) | Closed: the i⁴ = 1 eigenvector-frame holonomy around the defective seed is generic complex-symmetric EP2 geometry, a third clock of its own, not Π's Z₄; gated by 24 checks plus C# pins |

### The Lattice Thread (bridged worlds, dead sets, the opening law)

| Experiment | Key finding |
|-----------|------------|
| [Lattice Opening Law](LATTICE_OPENING_LAW.md) | Verified from below, with a C# pin: on the cat pair the lattice opening is exactly max(cos²θ, sin²θ) − cosθ·sinθ·e^(−2Γt), the heavier sock minus the living spook, J-blind because the cat sector is H-dead |
| [Dead-Set Rule](LATTICE_DEAD_SET_RULE.md) | Verified from below and minted as F132 (necessity Tier 1 derived, sufficiency Tier 2 gated): three conserved structures decide which Pauli readouts die identically. The inherited hand-spotted candidate rule was falsified |
| [Dead Set at h = 0](LATTICE_DEAD_SET_H_ZERO.md) | Verified from below: the F132 dead-set law is generic along the h axis (allowed readouts ≥ 5·10⁻³, forbidden at machine zero), and at h = 0 it refines because the Majorana hopping graph disconnects into two components |
| [Dead Set: the ZZ Face](LATTICE_DEAD_SET_ZZ_FACE.md) | Verified from below: turning on ZZ makes the F132 dead set jump to its full revived size already at zz = 10⁻⁴, while magnitudes cross over as zz^m with exact integer Majorana-degree orders, worst slope deviation 0.0006 |
| [Lattice H-Thread](LATTICE_H_THREAD.md) | Verified from below: X^N is a third mirror of the F131 order-sorting law. With a field the one-sided reading leaves the world-family entirely, satisfying a mixed two-field pencil, and a second antiunitary mirror forces 48 of 63 readouts to zero |

### The Frozen Band and the Ceiling

| Experiment | Key finding |
|-----------|------------|
| [What Reaches the Ceiling](WHAT_REACHES_THE_CEILING.md) | The measurement side of F145 and F146 (proof in [PROOF_SCALAR_COUNT](../docs/proofs/PROOF_SCALAR_COUNT.md)), kept in the order the afternoon went, including two readings killed within the hour: the floor-attaining states are one spin-1 per chiral pair, counted C(⌊N/2⌋, ℓ)·R_ℓ with R the Riordan number |
| [The Exceptional Couplings](THE_EXCEPTIONAL_COUPLINGS.md) | The measuring, deliberately unnamed and given no F number: the frozen-band ceiling fails at isolated real couplings. Existence proved, the (2,2)-singlet count exact through N = 8, the complete all-rung count exact only at N = 5 |
| [η-Ceiling Reduction](ETA_CEILING_REDUCTION.md) | The frozen-band ceiling certificate moved off the block onto V₀ = ker(ad_h), certified per rung by exact GF(q) rank; the law min 𝒦 on LW_ℓ ∩ V₀ = ℓ(N−ℓ)/(N+1) closes the ceiling for all N ≥ 6, its inequality half now proved, attainment still measured |
| [XY Frozen Band](XY_FROZEN_BAND.md) | On the R₉₀ locus the XY chain carries λ = −4γ̄ not just in F140's corner blocks but across the whole band \|p − q\| ∈ {0, 2}, 3(N − 1) blocks at depth ⌊N/2⌋. Census N = 4 to 7, band-edge probes N = 8, 9, 10, ceiling by exact GF(p) rank |
| [What the R₉₀ Locus Buys](WHAT_THE_R90_LOCUS_BUYS.md) | Why one locus carries both F153 and F140: the failure of F153's size-class condition on the corner and F140's defect are the same matrix 8γ̄·P_D, reached linearly by τQ and antilinearly by the bare reflection; each map satisfies exactly one of the two identities |
| [The Spread is a Resonance](THE_SPREAD_IS_A_RESONANCE.md) | The strong-coupling Δ selection dissolves: on index-N/2 blocks flat at −σ is generic and exact (an X^N density cancellation, any profile) and the spread needs a difference-spectrum collision; the Δ = 0.5 doublet is exact over ℤ, F153's N clause is the resonance count, and blocks off half filling carry a fourth, collision-free behaviour |
| [The Endpoints are a Density Law](THE_ENDPOINTS_ARE_A_DENSITY_LAW.md) | Why the resonances saturate exactly on the size-class centres: on scalar-parity eigenspaces the reflection-symmetric compressed density is a theorem of [H, R] = 0, the parity census shows that covers everything up to N = 7, and on the R₉₀ locus Π D Π = −2γ̄·Π N_XY Π gives containment and attainment in one identity; the parity-mixed eigenspaces of N = 8 at Δ = 0 obey the same law unforced, and Δ = 1 adds a uniform density s/N per size class |
| [Ceiling: Four Non-Local Cases](CEILING_FOUR_NONLOCAL_CASES.md) | The palindrome ceiling of non-local k-body cases, narrowed 6 → 4 → 2 → 0: the last Z-middle pair is palindromized by the period-4 golden router, so the arc is closed at zero non-local cases |

### F-Registry Readings (F64 through F130)

| Experiment | Key finding |
|-----------|------------|
| [F64 Topology Generalization](F64_TOPOLOGY_GENERALIZATION.md) | Tier 1: F64's cavity-mode-exposure rate α_k = 2γ_B·\|a_B(ψ_k)\|² generalizes to chain, star, ring, complete and Y-tree at N=5 and N=7, max relative error < 0.001 in the first-order regime |
| [F86: the EP Through the Clock](F86_EP_THROUGH_THE_CLOCK.md) | Reads F86a's exceptional point through the clock's Takt and Rotation hands. Status says a seeing, not a proof and not a closed form, with γ_crit(N) and K(N) deliberately left open |
| [F87 Windowed Converse per Block](F87_WINDOWED_CONVERSE_PER_BLOCK.md) | Localizes where the F87 palindrome pairing breaks, per frequency block, via the transpose relation M(−ω) = M(ω)^T, bit-exact at N=4. A sharpening, not a proof; superseded 2026-06-09 and closed by the two-reflection monomial theorem |
| [Bipartite Chirality, Diagonal Cell](BIPARTITE_CHIRALITY_DIAGONAL_CELL.md) | F87 soft ⟺ H's hopping graph is bipartite in the dephasing basis, bit-exact at k=3 and k=4 with zero mismatches. Bipartite ⟹ soft derived; the windowed converse closed as a theorem 2026-06-10 |
| [Softness Is N-Dependent](SOFTNESS_IS_N_DEPENDENT.md) | Tier 2 computed, bit-exact: the F87 soft/hard verdict is not stamped on the Hamiltonian. The 4-body witness XXXX+XYYY+YYYX is genuinely soft at N=5 (pairErr 8.7×10⁻¹⁴) and genuinely hard at N=6 (0.2), a finite-size crossing, not a tolerance artifact |
| [F112 Non-Hermitian Basis Enumeration](F112_NONHERMITIAN_BASIS_ENUMERATION.md) | F112's non-Hermitian extension: Tier-1 derived for all N via a two-lemma structural proof, with the N=2..6 basis enumeration (8,950,568 pair F-values, all bit-exact 0 or < 1e-10) preserved as the empirical anchor |
| [F113 Break-Magnitude Formula](F113_BREAK_MAGNITUDE_FORMULA.md) | Closed form for the F112 polarity-asymmetry break when the typed scope is violated: asymmetry = (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l). Tier 1 derived at N=2, 3, 4 by constructive parameter sweep |
| [F115 Obstruction Distribution](F115_OBSTRUCTION_DISTRIBUTION.md) | Tier 2 (computed, in progress): the F115 obstruction-size distribution collapses to the d=0 layer, the Δ-bucket count and the size-3 floor are closed, and the middle sizes are located as a number-theoretic hard core, not closed |
| [F129 Family Inventory](F129_FAMILY_INVENTORY.md) | Counts every F129 level-collision family: thirteen families with exact closed forms, the family list forced at every n, counts verified on every firing n ≤ 140 plus n = 150 and 210, and since derived |
| [F130 Time-Domain Decoupling](F130_TIME_DOMAIN_DECOUPLING.md) | F130 read as dynamics for the first time: a generic mode pair couples as q², the equal-level pair as q⁴ with the q² and q³ terms exactly absent, two orders of perturbation theory removed |

### Additional Experiments

| Experiment | Key finding |
|-----------|------------|
| [The Formation Window](THE_FORMATION_WINDOW.md) | Where a k-body bound complex forms on an XXZ chain: the marginal ridge is linear (chain Δ_ridge = 1.14 + 1.28·j2, ring 1.43 + 1.75·j2), the transition broadens linearly with integrability breaking (the ring ~4× faster), the edge is body-count independent at j2=0 but k=2 is the soft outlier, all N-converged. The marginal window read as the near-threshold resonance band (Hoyle) and the aromatic-ring geometry shift |
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
| [Neural Gamma Cavity](NEURAL_GAMMA_CAVITY.md) | WITHDRAWN 2026-08-25: 97.3% pairing was the matching tolerance, Dale's law made no difference, the 18 unpaired modes were an ordering artifact. Standing: a limit cycle at the same parameters, shortest sampled period 5.74 time constants, growing towards both folds; and a zero-multiplicity excess in the wiring. The band label does not stand: the integrated model has no time constant, so its Hz are a stipulation |
| [Trapped Light Localization](TRAPPED_LIGHT_LOCALIZATION.md) | K_death = 2.303 universal; surviving mode energy center-localized (ratio 1.3-1.4); N+1 immortal modes; gamma plays algebraic role of c (Tier 4-5) |
| [Hydrogen Bond Qubit](../docs/water/HYDROGEN_BOND_QUBIT.md) | Zundel cation proton crosses CΨ = 1/4 fold 6 times in 21 fs; every water molecule produces ~10-70 fold crossings per picosecond |
| [Primordial Qubit Algebra](PRIMORDIAL_QUBIT_ALGEBRA.md) | Π creates proper Z₂-graded super-algebra M_{2\|2}(ℂ) with block-off-diagonal L_c; Tomita-Takesaki connection ruled out (Π linear, J anti-linear) |
| [N=5 Check](N_EQUALS_FIVE_CHECK.md) | All six N-scaling metrics (max multiplicity, fraction distinct, slow-mode rate, slow/Σγ, max sector dimension) are monotonic with N. N=5 is not extremal on any axis; repeated N=5 appearance is selection bias from IBM Torino chain hardware, not physics |
| [Exchange from the V-Effect](EXCHANGE_FROM_V_EFFECT.md) | Tier 1-2: the V-Effect bridge between two Heisenberg pairs generates an effective exchange J_eff = (3/8)α²/J, the 3/8 prefactor derived from Pauli algebra alone, numerics at N=4 matching second-order PT to 1.4% at α=0.05 |
| [Asymmetric Exchange from the V-Effect](ASYMMETRIC_EXCHANGE_FROM_V_EFFECT.md) | Second-order PT for two asymmetric Heisenberg pairs bridged by α: δE_GS = −3α²/(4(J_A+J_B)), matching exact N=4 diagonalization to 0.2-0.8% at α=0.025. Tier 1-2, computational plus analytical |
| [Coupling Defect: the Walk-Time Step](COUPLING_DEFECT_WALK_TIME_STEP.md) | Pre-registered and computationally verified at N=7, 20, 60, 120: a single-bond defect δJ writes a pure downstream arrival-time step −δ/(2J), exact at first order. It is NOT the PTF's α_i object, which responds smoothly and nonlocally |
| [Front Pedigree](FRONT_PEDIGREE.md) | Computationally verified by three methods: the F126 renewal ladder resolved by catch count. Front survivors are mostly rebirths, lightly caught (⟨j⟩ = 1.38 against Γt* = 5), with last rebirths spread across the whole trip |
| [Quantum Darwinism: the Pointer Door](QUANTUM_DARWINISM_POINTER_DOOR.md) | Zurek's redundancy refcount computed in our N=8 chain. The arc minted F135 (record parity trichotomy) and F136 (record letter law): transport never broadcasts, fully-witnessed worlds are exactly the stars and complete graphs, the heavy-hex bulk is dark |
| [Neural Clock: Two Hands](NEURAL_CLOCK_TWO_HANDS.md) | Splits the neural Jacobian eigenvalue into a Takt hand and a Rotation hand. Tier 2 computational: the Takt identity mean Re λ = trace(J)/d = −S is exact and wiring-independent, confirmed to six decimals on three graphs |
| [Flavor-Resolved T2 Inheritance](FLAVOR_RESOLVED_T2_INHERITANCE.md) | Tier 2 numerical inheritance, N = 1..6: the two-flavor T2 lifetime split of the carbon ring reappears on a water proton wire and a neural network. The hoped-for clean fractions 4/3, 8/7 do NOT hold; the water and neural readings are Tier 3 |
| [The Flow Between Two Singularities](THE_FLOW_BETWEEN_TWO_SINGULARITIES.md) | Seen, not yet understood; parked for future us. Numbers verified to machine precision, interpretation deliberately left open: the settled state is never reached yet already present at t = 0, and only the closed forms sit at the two edges |

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
simulation script. All use NumPy + SciPy. Typical runtime: seconds to
minutes on a standard laptop.

**If you want a guided reading path:** See the [Reading Guide](../docs/READING_GUIDE.md),
which organizes the experiments into three stories: the proof, the
application, and the ontology.
