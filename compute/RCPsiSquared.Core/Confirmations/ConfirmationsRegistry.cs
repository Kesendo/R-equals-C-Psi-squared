namespace RCPsiSquared.Core.Confirmations;

/// <summary>Hardware-confirmed framework predictions, indexed by name.
///
/// Mirrors <c>simulations/framework/confirmations.py</c>. Each entry binds a framework
/// prediction to a specific IBM Quantum run with measured values, raw-data pointer, and
/// experiment writeup. Used to look up "is this prediction hardware-confirmed?" without
/// re-deriving.
/// </summary>
public static class ConfirmationsRegistry
{
    private static readonly IReadOnlyList<Confirmation> _all = new[]
    {
        new Confirmation(
            Name: "palindrome_trichotomy",
            Date: "2026-04-26",
            Machine: "ibm_marrakesh",
            JobId: "d7mjnjjaq2pc73a1pk4g",
            Observable: "<X_0 Z_2>",
            PredictedValue:
                "continuous_lindblad γ_Z=0.1: truly=0.000 soft=-0.623 hard=+0.195 (Δ_soft-truly=-0.623); " +
                "trotter_n3 γ_Z=0.1: truly=0.000 soft=-0.723 hard=+0.327 (Δ_soft-truly=-0.723)",
            MeasuredValue:
                "truly=+0.011 soft=-0.711 hard=+0.205 (Δ_soft-truly=-0.722)",
            HardwareData: "data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json",
            ExperimentDoc: "experiments/MARRAKESH_THREE_LAYERS.md",
            FrameworkPrimitive: "classify_pauli_pair + predict_M_spectrum_pi2_odd + propagate_with_hardware_noise",
            Description:
                "Super-operator palindrome trichotomy (truly/soft/hard) tomographically distinguishable on Heron r2 hardware at N=3. " +
                "Hardware Δ(soft − truly) = -0.722 matches the Trotter n=3 prediction to 0.0014, NOT the continuous-Lindblad idealization (Δ = -0.623). " +
                "Original interpretation that T1 amplifies the soft-break is REFUTED: T1 monotonically attenuates |Δ| (γ_T1=0.5 gives Δ=-0.44). " +
                "The hardening is Trotter discretization at δt=0.267 with ‖H‖_op = 2.83·J, where ‖H·δt‖ ≈ 0.76 violates the small-step regime.",
            QubitPath: new[] { 48, 49, 50 }),

        new Confirmation(
            Name: "f25_cusp_trajectory",
            Date: "2026-04-26",
            Machine: "ibm_kingston",
            JobId: "d7mu36lqrg3c738lnda0",
            Observable: "CΨ(t) for Bell+",
            PredictedValue: "F25: CΨ(t) = f·(1+f²)/6 with f = exp(-4·γ·t)",
            MeasuredValue: "RMS residual 0.0097 vs in-situ γ_fit",
            HardwareData: "data/ibm_cusp_slowing_april2026/ (April 26 precision run)",
            ExperimentDoc: "experiments/CRITICAL_SLOWING_AT_THE_CUSP.md",
            FrameworkPrimitive: "F25 closed-form CΨ(t)",
            Description: "Bell+ trajectory through CΨ=1/4 cusp confirmed point-by-point on Kingston (19 delay points, qubits 14-15)."),

        new Confirmation(
            Name: "f57_kdwell_gamma_invariance",
            Date: "2026-04-16",
            Machine: "ibm_kingston",
            JobId: "cusp_slowing_kingston_20260416",
            Observable: "K_dwell / δ for Bell+",
            PredictedValue: "F57: γ-independent, F25 prefactor 1.0801 for pure Z-dephasing",
            MeasuredValue: "Pair A 0.6492, Pair B 0.6937 (6.3% spread despite 2.55× γ ratio)",
            HardwareData: "data/ibm_cusp_slowing_april2026/cusp_slowing_ibm_kingston_20260416_212042.json",
            ExperimentDoc: "experiments/CRITICAL_SLOWING_AT_THE_CUSP.md",
            FrameworkPrimitive: "K_dwell formula in F57",
            Description: "γ-invariance of K_dwell at CΨ=1/4 boundary verified on Kingston with two qubit pairs at 2.55× different γ. " +
                "Absolute prefactor 0.67 vs predicted 1.08 due to T1 amplitude damping."),

        new Confirmation(
            Name: "bonding_mode_receiver",
            Date: "2026-04-24",
            Machine: "ibm_kingston",
            JobId: "multiple (see external pipeline)",
            Observable: "MI(0, N-1) for bonding:2 vs alt-z-bits",
            PredictedValue: "4000-5500× over ENAQT in simulation N=5..13; ratio ≈ 1.4-3× hardware-realistic",
            MeasuredValue: "bonding:2 / alt-z-bits = 2.80× on Kingston N=5",
            HardwareData: "external (AIEvolution.UI/experiments/ibm_quantum_tomography)",
            ExperimentDoc: "experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md",
            FrameworkPrimitive: "F67 bonding-mode + F71 chain-mirror symmetry",
            Description: "Receiver-engineering bonding-mode advantage measured on Kingston. " +
                "Largest engineering lever in the framework portfolio. Confirms F67 + F71 structure on hardware."),

        new Confirmation(
            Name: "chiral_mirror_law",
            Date: "2026-04-25",
            Machine: "ibm_marrakesh",
            JobId: "k_partnership_marrakesh_20260425",
            Observable: "Bloch components of K-partner pair states",
            PredictedValue: "Chiral mirror identity from K=diag((-1)^l)",
            MeasuredValue: "retrospective verification, identity holds",
            HardwareData: "data/ibm_k_partnership_april2026/k_partnership_marrakesh_20260425_140913.json",
            ExperimentDoc: "experiments/CHIRAL_MIRROR_HARDWARE_PREDICTION.md",
            FrameworkPrimitive: "K_full chiral conjugation (Symmetry/ChiralK)",
            Description: "Chiral mirror law (K H K = -H spectral inversion) verified retrospectively on Marrakesh K-partnership data."),

        new Confirmation(
            Name: "pi_protected_xiz_yzzy",
            Date: "2026-04-26",
            Machine: "ibm_marrakesh",
            JobId: "d7n3013aq2pc73a2a18g",
            Observable: "<X_0 I Z_2>",
            PredictedValue: "protected (≈0) for YZ+ZY soft Hamiltonian",
            MeasuredValue: "+0.13 to +0.04 (within noise band, never above ±0.13)",
            HardwareData: "external (raw JSON not in repo; results documented inline in experiment doc)",
            ExperimentDoc: "review/EMERGING_QUESTIONS.md",
            FrameworkPrimitive: "pi_protected_observables",
            Description: "First-time hardware measurement of a Π-protected observable on YZ+ZY soft Hamiltonian (EQ-030). " +
                "Confirms framework primitive at hardware scale on a Hamiltonian not previously tested.",
            QubitPath: new[] { 0, 1, 2 }),

        new Confirmation(
            Name: "marrakesh_transverse_y_field_detection",
            Date: "2026-04-29",
            Machine: "ibm_marrakesh",
            JobId: "d7ornigror3c73c0c6ug",
            Observable: "Z⊗N-Mirror max-violation between ρ_a=|+−+⟩ and ρ_b=|−+−⟩",
            PredictedValue:
                "h_y=0.05 → max_violation ≈ 0.18 (linear scaling 3.5·h_y); " +
                "h_x=0.05 → max_violation ≈ 0.004 (linear scaling 0.085·h_x); " +
                "clean stack → max_violation < 1e-3",
            MeasuredValue:
                "max_violation = 0.182 (worst Pauli: Z,Z), RMS = 0.087. Matches h_y=0.05 prediction exactly. h_y_eff ≈ 0.05.",
            HardwareData: "data/ibm_zn_mirror_april2026/zn_mirror_ibm_marrakesh_20260429_102824.json",
            ExperimentDoc: "data/ibm_zn_mirror_april2026/README.md",
            FrameworkPrimitive: "chain.zn_mirror_diagnostic + zn_mirror_state",
            Description: "First hardware verification of the Z⊗N-Mirror diagnostic. Marrakesh shows effective transverse Y-field h_y_eff ≈ 0.05 at Hamiltonian level on path [48,49,50]. " +
                "NOT a transverse X-field (which would give 40× smaller violation). The Y-vs-X asymmetry predicted by the framework (Y is bit_b-odd like Z-dephasing axis, mixes more strongly) is confirmed empirically.",
            QubitPath: new[] { 48, 49, 50 }),

        new Confirmation(
            Name: "lebensader_skeleton_trace_decoupling",
            Date: "2026-04-26",
            Machine: "ibm_marrakesh",
            JobId: "d7n3013aq2pc73a2a18g + d7n3eqqt99kc73d34qtg",
            Observable: "Π-protected counts (skeleton) + θ-trajectory tails (trace)",
            PredictedValue:
                "Skeleton/trace decouple at N≥4, co-occur at N=3. " +
                "T1 makes drop measurable: bond-flipped Z-free pairs (XY+YX, IY+YI) preserve skeleton (drop≤1) and trace (long θ-tail); " +
                "Z-containing soft pairs (YZ+ZY, XZ+ZX) collapse both (drop≈28-29, no tail).",
            MeasuredValue: "drop=28 for YZ+ZY confirmed on Marrakesh. Pearson(drop, Δ∫θ)=+0.85. Bures velocity gives no third discriminator.",
            HardwareData: "external (raw JSON not in repo; tables of <X_0 I Z_2> per t and basis documented inline in experiment doc)",
            ExperimentDoc: "review/EMERGING_QUESTIONS.md",
            FrameworkPrimitive: "cockpit_panel: composes pi_protected_observables + θ-trajectory",
            Description:
                "Lebensader as Stromkabel (EQ-030 closure): skeleton (Π-protected algebraic count) and trace (θ-geometric tail) are NOT two discriminators but one bridge held together by Π·L·Π⁻¹ + L + 2Σγ·I = 0. " +
                "Hardware-confirmed across 3 of 4 bond-flipped Z-free corners; Bures velocity confirmed null as third axis.",
            QubitPath: new[] { 0, 1, 2 }),

        new Confirmation(
            Name: "block_cpsi_saturation_kingston_may2026",
            Date: "2026-05-08",
            Machine: "ibm_kingston",
            JobId: "d7ulfjdpa59c73b4rttg",
            Observable: "C_block(t) = Σ |ρ_{ab}|² over (popcount-0, popcount-1) on 2-qubit ρ from Dicke initial (|D_0⟩+|D_1⟩)/√2",
            PredictedValue:
                "Theorem 1: C_block(t=0) = 1/4 EXACTLY (universal Mandelbrot-cardioid ceiling). " +
                "Theorem 3 trajectory: C_block(t) = (1/4)·exp(-4γ·t) with γ = 1/(2·T2) for pure Z-dephasing. " +
                "From the 2026-05-08 calibration T2_min = 480 μs on q13–q14: γ_expected ≈ 1.04e-3 μs⁻¹.",
            MeasuredValue:
                "C_block(0) = 0.2205 = 88.2% of 1/4 (state-prep + tomography fidelity floor). " +
                "Trajectory at t ∈ {0, 120, 240, 360, 480} μs: {0.2205, 0.1023, 0.0458, 0.0166, 0.0074}. " +
                "Log-linear fit: γ_fit = 1.795e-3 μs⁻¹, R² = 0.9977. " +
                "T2_eff from C_block decay = 1/(2·γ_fit) = 278.5 μs vs calibrated T2_min = 480 μs (1.72× hardware speedup).",
            HardwareData: "data/ibm_block_cpsi_saturation_may2026/block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json",
            ExperimentDoc: "experiments/IBM_BLOCK_CPSI_SATURATION.md",
            FrameworkPrimitive:
                "BlockCoherenceContent.Compute + BlockCpsiClosedForm.At + IbmBlockCpsiHardwareTable (typed); " +
                "state-level Theorem 1 + Theorem 3 trajectory of PROOF_BLOCK_CPSI_QUARTER",
            Description:
                "Direct hardware confirmation of the universal-Mandelbrot-cardioid 1/4 ceiling on the (popcount-0, popcount-1) coherence block " +
                "of a 2-qubit ρ. Initial state (|D_0⟩+|D_1⟩)/√2 saturates Theorem 1 at 88.2% of the 1/4 ceiling on Kingston q13–q14. " +
                "The Theorem-3 closed-form trajectory C_block(t) = (1/4)·exp(-4γ·t) fits the 5 t-points with R² = 0.9977 (pristine log-linear decay). " +
                "Fitted γ_fit = 1.795e-3 μs⁻¹ corresponds to T2_eff = 278.5 μs from the C_block decay alone, vs T2_min = 480 μs from the morning calibration; " +
                "the 1.72× speedup quantifies the gate-noise + readout contribution beyond pure single-qubit T2 dephasing. " +
                "First Tier-2-Verified hardware anchor for the Mandelbrot 1/4 boundary on a state engineered to sit on it.",
            QubitPath: new[] { 13, 14 }),

        new Confirmation(
            Name: "f83_pi2_class_signature_marrakesh",
            Date: "2026-04-30",
            Machine: "ibm_marrakesh",
            JobId: "d7pol1e7g7gs73cf7j90",
            Observable: "Per-category 2-qubit Pauli expectations (q0=4, q2=6) for 4 F87 Π²-classes",
            PredictedValue:
                "Trotter+γ_Z patterns at the path-fit γ_Z_eff: pi2_odd_pure ⟨X₀Z₂⟩, pi2_even_nontruly ⟨X₀X₂⟩, mixed ⟨Z₀X₂⟩, truly ⟨Y₀Z₂⟩. " +
                "F83 anti-fractions = {1/2, 0, 1/6} structurally (closed form on H letters).",
            MeasuredValue:
                "4-category discrimination CONFIRMED via unique-fingerprint Paulis: " +
                "pi2_odd_pure ⟨X₀Z₂⟩=-0.849, pi2_even_nontruly ⟨X₀X₂⟩=+0.919, mixed ⟨Z₀X₂⟩=-0.721, truly ⟨Y₀Z₂⟩=+0.670. " +
                "All separations >>10σ at 4096 shots. Path-fit γ_Z_eff = 0.050.",
            HardwareData: "data/ibm_f83_signature_april2026/f83_signature_ibm_marrakesh_20260430_190035.json",
            ExperimentDoc: "data/ibm_f83_signature_april2026/README.md",
            FrameworkPrimitive: "classify_pauli_pair + predict_pi_decomposition (F83 anti-fraction closed form)",
            Description:
                "F83 4-Hamiltonian Π²-class discrimination test on path [4,5,6]. Each of the 4 F87 classes shows a unique-fingerprint Pauli observable separating it from the other 3 at >>10σ. " +
                "γ_Z_eff is path-dependent (0.05 vs 0.12 between [4,5,6] and [48,49,50]), reflecting that effective dephasing absorbs Trotter discretization, coherent gate errors, and crosstalk.",
            QubitPath: new[] { 4, 5, 6 }),

        new Confirmation(
            Name: "regime_uniformity_kingston_uniform_quantum",
            Date: "2026-05-05",
            Machine: "ibm_kingston",
            JobId: "d7sqjpiudops73976960",
            Observable: "F88b-Lens Π²-odd-memory across F87 trichotomy on uniform-quantum chain",
            PredictedValue:
                "regime-uniformity hypothesis: a uniform-quantum chain's truly-baseline should be near a uniform-classical chain's, both well below the regime-mixed value 0.0297 (Marrakesh [0,1,2] framework_snapshots). " +
                "soft-pumping should remain hardware-substrate-independent at ~0.74.",
            MeasuredValue:
                "truly = 0.0022, soft = 0.7409, pi2_even_nontruly = 0.0046, mixed = 0.6220. " +
                "truly is 13.5× lower than the regime-mixed Marrakesh [0,1,2] (0.0297) and only 1.69× higher than uniform-classical Marrakesh [48,49,50] (0.0013); the 22.8× gap between regime-mixed and uniform-classical from the original framework_snapshots vs soft_break comparison sits an order of magnitude wider than the within-uniform variation. " +
                "soft is 3.1% from Marrakesh's 0.7646: substrate-independent within shot noise.",
            HardwareData: "data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json",
            ExperimentDoc: "simulations/_f88b_lens_ibm_kingston_uniform_quantum.py",
            FrameworkPrimitive: "F88b-Lens (kernel projection + Π²-odd Pauli enumeration on reduced 2-qubit ρ)",
            Description:
                "First F87 trichotomy hardware test on a uniform-quantum CZ-coupled triple. Path [43, 56, 63] on Kingston, all three qubits classified PulseStable across the 91-day biography window (r mean 0.10 / 0.09 / 0.10, all crossing > 95%). " +
                "Result confirms regime-uniformity hypothesis: BOTH uniform-classical (Marrakesh [48,49,50] truly = 0.0013) AND uniform-quantum (Kingston [43,56,63] truly = 0.0022) give clean truly-baselines, while regime-mixed (Marrakesh [0,1,2] truly = 0.0297) is an order of magnitude dirtier. " +
                "Three findings per single 39-second billed run: regime-uniformity is the cause of the truly-baseline gap (not which side of the boundary); F87 trichotomy hardware-confirmed on a second backend (Kingston) in addition to Marrakesh; soft Π²-odd-pumping confirmed substrate-independent across two Heron-r2 chips.",
            QubitPath: new[] { 43, 56, 63 }),

        new Confirmation(
            Name: "f95_angle_steering_kingston_may2026",
            Date: "2026-05-16",
            Machine: "ibm_kingston",
            JobId: "bxyj5yd4j + bzklqwt7f",
            Observable: "arg(CΨ_com(t)) at first real-axis magnitude-minimum crossing under deliberate RZ phase injection at rate Ω (two Ω values: 0.13 and 0.25 rad/μs; per-Ω mapping in MeasuredValue)",
            PredictedValue:
                "F95 complex-CΨ angle is actively steerable via per-chunk RZ(Ω·Δt) injection on the bra-ket Hermitian off-diagonal. " +
                "Lindblad+RZ model with in-situ γ_per_us from T2 calibration predicts arg(CΨ_com) at the measured crossing time. " +
                "Pair A_mid [82,83] in-situ γ = 3.36e-3 μs⁻¹ (T2_min = 148.7 μs); Pair B_high [13,14] in-situ γ = 1.90e-3 μs⁻¹ (T2_min = 263.7 μs).",
            MeasuredValue:
                "Pair A_mid [82,83] Ω=0.13: t_cross = 1.395 μs, arg_measured = -84.70°, arg_predicted = -100.39°, residual = +15.69°. " +
                "Pair A_mid [82,83] Ω=0.25: t_cross = 1.242 μs, arg_measured = -98.64°, arg_predicted = -107.79°, residual = +9.15°. " +
                "Pair B_high [13,14] Ω=0.25: t_cross = 2.814 μs, arg_measured = -123.50°, arg_predicted = -130.31°, residual = +6.81°. " +
                "(Pair B_high Ω=0.13 produced no magnitude crossing in the delay window: higher T2 + lower Ω left the trajectory above the real-axis crossing threshold.) " +
                "Residuals monotonically shrink as Ω increases (15.7° → 9.2° at Ω 0.13 → 0.25 on Pair A_mid), consistent with a roughly constant natural-Kingston drift contribution becoming a smaller fraction of the driven rotation.",
            HardwareData: "data/ibm_f95_angle_steering_may2026/cusp_complex_phase_hardware_ibm_kingston_omega0.{130,250}_20260516_*.json (+ matching PNGs and partial-run JSONs; see README.md in that directory for details)",
            ExperimentDoc: "experiments/CPSI_COMPLEX_PLANE.md (conceptual predecessor; F95 closed form derived 2026-05-16 same day)",
            FrameworkPrimitive: "F95AngleAtQuadraticZeroPi2Inheritance (θ(c;b) = arctan(√(c/b² − 1))) + " +
                "F25 Lindblad CΨ trajectory + complex CΨ_com signed-sum-of-off-diagonals extension",
            Description:
                "First hardware verification that the complex-CΨ angle predicted by F95 is actively steerable on IBM Kingston. " +
                "Two Ω values (0.13 and 0.25 rad/μs) tested on two qubit pairs (A_mid mid-T2 [82,83], B_high high-T2 [13,14]) in 14 minutes total billed QPU. " +
                "Three of four measurement conditions produced detectable real-axis crossings; predicted arg matches measured to within 16° across all three (best 6.8°). " +
                "Confirms three open questions from CPSI_COMPLEX_PLANE.md simultaneously: " +
                "(Q1) Δφ = Ω·t reproducible: crossings occur at predicted angle modulo bounded natural drift; " +
                "(Q2) linear-in-Ω scaling holds: residual shrinks ~1.7× when Ω doubles, consistent with constant natural drift contribution; " +
                "(Q3) active steering works: the same Lindblad+RZ model fits both Ω values on the same pair without retuning. " +
                "The bra↔ket Hermitian off-diagonal of ρ on the prepared Bell-like initial state carries a well-defined complex phase, addressable on Heron r2 by per-chunk RZ injection on the steering qubit. " +
                "Hardware substrate-of-evidence: F95 is the universal polynomial-foundation algebra of the angle that appears off the d=0 mirror; this run shows that algebra is operationally controllable, not merely descriptive.",
            QubitPath: new[] { 82, 83 }),

        new Confirmation(
            Name: "gamma0_off_the_lever_kingston_may2026",
            Date: "2026-05-29",
            Machine: "ibm_kingston",
            JobId: "d8ce8l38ch0s738uorjg",
            Observable: "Transfer T(t) = P(qubit 1 excited) from |10⟩ under exchange H = J·(XX+YY)/2 + idle Z-dephasing, scanned over J (the only lever)",
            PredictedValue:
                "Q = J/γ₀ with the typed carrier γ₀ = UniversalCarrierClaim.DefaultGammaZero = 0.05. " +
                "Coherent swing (T overshoots ½, the H-clock wins) for Q > 1 (J > γ₀); overdamped creep (T → ½, the carrier wins) for Q ≤ 1; " +
                "critical damping exactly at J = γ₀ = 0.05 (Q = 1, max T = ½). Sim _q_basic_jscan.py puts the swing-death at J = γ₀.",
            MeasuredValue:
                "q13-q14, J ∈ {0.05, 0.1, 0.2, 0.4} (Q ≈ {1, 2, 4, 8}): max T = {0.335, 0.563, 0.703, 0.779}. " +
                "Coherent swing (overshoot past ½) present for J ≥ 0.1, absent for J = 0.05 (creeps to 0.335). " +
                "First transfer peak at step ∝ 1/J (the H-clock frequency tracks J). " +
                "Swing dies between J = 0.1 and J = 0.05 → γ₀ ≈ 0.05–0.1 read off the lever (idle carrier + gate error shift the threshold slightly above the pure 0.05).",
            HardwareData: "external (AIEvolution.UI/experiments/ibm_quantum_tomography/results/q_jscan_hardware_ibm_kingston_20260529_025757.json)",
            ExperimentDoc: "experiments/GAMMA0_IS_ALWAYS_THERE.md",
            FrameworkPrimitive: "UniversalCarrierClaim.DefaultGammaZero (γ₀ = 0.05) + Q = J/γ₀ regime axis; run_q_jscan.py (external pipeline)",
            Description:
                "First direct read-off of the carrier γ₀ from its only lever J on hardware. The basic Q = J/γ₀ exchange-bond J-scan: " +
                "the H-clock frequency tracks J (first transfer peak ∝ 1/J), and the coherent swing dies where J crosses γ₀ (critical damping, T = ½). " +
                "Confirms the typed UniversalCarrierClaim.DefaultGammaZero = 0.05 on Kingston q13-q14, made readable via the swing-death threshold even though γ₀ is invisible head-on (Inside-Observability: only Q = J/γ₀ is accessible from inside). " +
                "Within-state ratios (the popcount/absorption ladder) wash on hardware (distance-blind noise, nothing cancels); this between-regime J-scan is the clean carrier probe. " +
                "Sits on the same Q = J/γ₀ axis as the F86 Q-peak (Q ≈ 1.5, where |K_CC_pr| maxes): the residual F86 empirical deviation is grid/resolution-limited (fine-grid scans N≥9 pending), not a break in the structure now hardware-anchored.",
            QubitPath: new[] { 13, 14 }),

        new Confirmation(
            Name: "ibm_ep_onset_may2026",
            Date: "2026-05-31",
            Machine: "ibm_kingston",
            JobId: "d8dr7dfd0j8c73f4man0 (Part A, flow to 1/3) + d8drjbfd0j8c73f4mobg (Part B, EP onset)",
            Observable:
                "Per-site populations ⟨n_l⟩ of a single-excitation walk on a 3-site chain (q13-q14-q15, Z-basis, no tomography); " +
                "Part B revival = max ⟨n_0⟩ for t ≥ 2 μs under injected random-Z-twirl dephasing, Q = J/γ swept over {0.5, 1, 1.5, 2.5, 5, 20}",
            PredictedValue:
                "SE-walk overdamped→revival handover reading: revival pinned at the 1/N = 1/3 equipartition floor below Q ≈ 1.5 (overdamped, no memory return), " +
                "lifting off above it (the single-excitation walk's critical-damping transition; real, measured). " +
                "Whether that SE transition is itself a genuine defective EP is under a SEPARATE open review (the coherence-horizon √-EP cluster, inspect --root horizon); not asserted here. " +
                "Twirl simulate (K=16 exact statevector): revival 0.31 → 0.84 across the same Q scan. " +
                "Part A flow target: per-site populations converge to 1/N = 1/3 at late t (the post-handover flow fixed point).",
            MeasuredValue:
                "Q = {0.5, 1, 1.5, 2.5, 5, 20} → revival = {0.30, 0.36, 0.34, 0.49, 0.56, 0.70}: " +
                "floor ~1/3 for Q ≤ 1.5, liftoff 0.49 → 0.56 → 0.70 as Q crosses 2.5 → 5 → 20. " +
                "Part A late-t populations 0.34 / 0.43 / 0.34 at 20 μs (converging to 1/3). " +
                "High-Q side suppressed (0.70 vs 0.84 simulate at Q=20) by two-qubit gate error (~160 RZZ gates by 20 μs at ~0.5%); the floor and the onset are clean.",
            HardwareData:
                "data/ibm_ep_onset_may2026/ (Part A ep_onset_hardware_ibm_kingston_20260531_060943.json, " +
                "Part B ep_onset_hardware_ep_ibm_kingston_20260531_064022.json + same-day simulate JSONs)",
            ExperimentDoc: "experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md",
            FrameworkPrimitive:
                "ExceptionalPointClock (the toy 2×2 reduction: decay pinning at 4γ₀, F95 rotation angle, eigenvector overlap min(x,1/x)) + " +
                "EpField hardware node (inspect --axis ep); Q ≈ 1.5 handover marker",
            Description:
                "The single-excitation-walk overdamped→revival handover watched switching the memory on, on a real chip, populations only. " +
                "Part A (job d8dr7dfd0j8c73f4man0) at the chip's natural Q ≫ 1: the excitation sloshes 0 → 2 → 1 → 0 with ~3 μs period (the reborn memory), " +
                "the site-0 revival fades 0.84 → 0.43 over 15 μs (the forgetting), and the populations converge to 1/3 = 1/N at 20 μs (the flow target). " +
                "Part B (job d8drjbfd0j8c73f4mobg) injects dephasing via a random-Z twirl (K=16 instances; the RZ gates are virtual on IBM, so the injection is error-free) " +
                "to push Q = J/γ down through Q ≈ 1.5: the revival sits on the equipartition floor (~1/3) for Q ≤ 1.5 and lifts off as Q crosses 1.5 → 2.5 " +
                "(0.34 → 0.49 → 0.56 → 0.70): a clean overdamped→revival handover at Q ≈ 1.5 on real hardware (the SE walk's critical-damping transition; real, measured). " +
                "Whether that SE transition is itself a genuine defective EP is under a SEPARATE open review (the coherence-horizon √-EP cluster, inspect --root horizon); not asserted here. " +
                "The F86a coherence-block 'real-axis EP' this entry formerly cited was retracted 2026-06-21: the full (n,n+1)-coherence block has no eigenvalue coalescence on the real Q axis (genuine non-normality there, large but finite Petermann; the EP is not on the real axis). " +
                "The revival decay envelope is gate-cost-limited (Trotterization, ~9 μs), not T2-limited (~200 μs); only the rate is gate cost, the floor and the onset are physics. " +
                "This table is the hardware node of EpField (Diagnostics/Foundation/EpField.cs) and the overlay in simulations/ep_transition.py.",
            QubitPath: new[] { 13, 14, 15 }),

        new Confirmation(
            Name: "f120_moment_tower_kingston_june2026",
            Date: "2026-06-11",
            Machine: "ibm_kingston",
            JobId: "d8l6c7rqv2lc73863acg (Arm A) + d8l6c832d42s73cb16a0 (Arm B) + d8l6h03nn5bs738rmrug (standard-T1 arbiter)",
            Observable:
                "Energy-moment slopes d/dt⟨H_p^j⟩ from the maximally mixed state (8-basis-state average) under pure idle; " +
                "H_p = X₀+X₀Z₁+0.7·X₁X₂ (girth-2 witness, H_p² = 2.49·I+2·Z₁+1.4·XXX exact, tower t₁ ≡ 0, t₂ = [0,16,0]); " +
                "qubits q149/q13/q9, NO two-qubit gates, dynamical decoupling disabled, τ ∈ {0..150} μs, two arms permuting the middle qubit",
            PredictedValue:
                "Rung-1 null: slope⟨H_p⟩ = 0 at all τ (t₁ ≡ 0; evolution-blind, robust against all Z-flavored idle parasitics). " +
                "Rung 2 fires: slope⟨H_p²⟩ = 2·Δγ_mid (t₂ fires at the middle site only); row identity ⟨H_p²⟩ = 2.49 + 2⟨Z₁⟩ + 1.4⟨XXX⟩; " +
                "site tracking across arms. Rates-from-calibration Δγ_l = 1/T1_l was the textbook noise-model assumption (the layer the chip declined).",
            MeasuredValue:
                "Double null HELD: slope⟨H⟩ = +2.4e-4/μs (z = +1.47, Arm A) and −6.8e-6/μs (z = −0.04, Arm B); |⟨XXX⟩| ≤ 0.02 throughout. " +
                "Row identity exact in every measured row (τ=100 Arm A: 2.49+2·0.4147+1.4·0.0186 = 3.345 vs 3.3455). " +
                "Per-qubit pump slopes (Arm A / B, per μs): q149 2.327e-3 / 2.193e-3, q13 3.029e-3 / 3.090e-3, q9 5.794e-3 / 5.779e-3 " +
                "(cross-arm reproducibility 5.7% / 1.9% / 0.3%). Arbiter standard T1: q149 424.6, q13 430.3, q9 99.9 μs. " +
                "IN-SITU model test (prep-conditioned split of the same circuits; bound ⟺ s₀ ≤ 0): pump/Γ = 0.972-0.979 (q149), " +
                "0.965-0.966 (q13), 0.994-0.996 (q9) — the bound HOLDS everywhere in-situ; the 1-3% margins read the per-qubit thermal " +
                "population (q13 1.7%, q149 1.1-1.4%, q9 0.2-0.3%). CORRECTED SAME DAY: the first reading (q13 violates pump ≤ Γ at " +
                "4-6σ) compared the run against the 16-minutes-later arbiter and was an EPOCH ARTIFACT: T1 telegraphs on minute scale " +
                "(q13 ~315 μs in-run vs 430 at the arbiter; q9 ~172 vs ~75-100; q149 stable across the morning epochs). " +
                "Re-analysis: simulations/f120_prep_split_reanalysis.py. TELEGRAPH CHASE (12:02Z, job d8l8f7r2d42s73cb3q7g, " +
                "16 self-arbitrating blocks): within-job Γ flat at the shot-noise floor (no switch in 75 s), s₀ ≤ 0 in 47/48 " +
                "block-cells; across the day EVERY qubit moved, including the stable control (q149 430 → ~285 μs; q13 in a third " +
                "state ~200 μs; q9 ~108). Verdict: q13 was never special, the device T1 landscape breathes everywhere by 1.5-2x " +
                "on minute-to-hour scales, the two-level model holds within any ~minute window; the hidden variable was time.",
            HardwareData:
                "data/ibm_moment_tower_june2026/ (moment_tower_hardware_ibm_kingston_20260611T073908Z.json, " +
                "t1_arbiter_20260611T075540Z.json, the 06:33Z calibration snapshot)",
            ExperimentDoc: "experiments/F120_MOMENT_TOWER_KINGSTON.md",
            FrameworkPrimitive:
                "moment_tower + predict_pump_slope + f113_bridge_asymmetry_from_slope (framework diagnostics f120_moment_tower, " +
                "called by the pipeline script run_moment_tower.py at startup); MomentTowerPumpChannelClaim (Core/Symmetry); " +
                "PROOF_MOMENT_TOWER_PUMP_CHANNEL",
            Description:
                "F120's first hardware reading, the cleanest protocol ever sent to a QPU from this project: not one entangling gate " +
                "(X-gate preparation of all 8 basis states = exact I/d, pure idle as evolution, single-qubit basis rotations to read). " +
                "The structural law confirmed: the first energy-moment rung stays silent (the double null) while the second fires as " +
                "exactly twice the middle qubit's pump curve (row-exact identity), the girth read from hardware is 2, the firing follows " +
                "the middle-qubit identity across arms, and the per-qubit pump rates reproduce across arms to 0.3-5.7%. " +
                "The rate layer told a two-act story: the first reading (q13 violates pump ≤ Γ) was corrected the same day by the " +
                "prep-conditioned re-analysis — the 8-basis-state preparation contains the |0⟩- and |1⟩-branches, so pump AND Γ come " +
                "from the same circuits, epoch-matched, and the bound holds everywhere in-situ (worst 0.996) with margins that READ the " +
                "per-qubit thermal population. What the cross-epoch comparison had actually detected is minute-scale T1 telegraphing on " +
                "q13 (~315 ↔ 430 μs) and q9 (~172 ↔ ~75-100 μs), once visible inside a single arm. Two lessons banked: the protocol is " +
                "self-arbitrating (pump, Γ, γ↑ from one circuit set; the in-situ model test is s₀ ≤ 0), and the hardness rung of a " +
                "programmed Hamiltonian is now a quantity a chip measures about itself by decaying.",
            QubitPath: new[] { 149, 13, 9 }),

        new Confirmation(
            Name: "gamma_0_marrakesh_calibration",
            Date: "2026-04-29",
            Machine: "ibm_marrakesh",
            JobId: "d7mjnjjaq2pc73a1pk4g",
            Observable: "<X_0 Z_2> for 3 Pauli-pair Hamiltonians (truly XX+YY, soft XY+YX, hard XX+XY)",
            PredictedValue:
                "continuous-equivalent γ_Z ≈ 0.05 (continuous Lindblad fit); " +
                "Trotter-modeled γ_Z ≈ 0.1 (the same 0.1 the framework idealized prediction used)",
            MeasuredValue:
                "best-fit γ_Z = 0.05 (sweep over [0.01, 0.15] with 71 points), total residual² = 6.4e-4 across 3 Hamiltonians (continuous Lindblad, no Trotter modeling). " +
                "When Trotter n=3 is modeled, the same data fits γ_Z = 0.1 exactly via Δ matching to 0.0014.",
            HardwareData: "data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json",
            ExperimentDoc: "data/ibm_soft_break_april2026/README.md",
            FrameworkPrimitive: "ChainSystem.propagate_with_hardware_noise + 2D fit (continuous Lindblad)",
            Description:
                "Continuous-Lindblad fit of γ_Z against Marrakesh ⟨X₀Z₂⟩ data converges to 0.05, which absorbs the Trotter n=3 discretization into a lower effective γ_Z. " +
                "The 2026-04-30 follow-up showed a Trotter-modeled fit returns γ_Z = 0.1 with Δ-matching to 0.0014, so the original framework-idealized γ_Z = 0.1 was correct (not 2× too high). " +
                "The two values 0.05 (continuous) and 0.1 (Trotter) are the same data through two physics models; T1 contributes negligibly in either.",
            QubitPath: new[] { 48, 49, 50 }),

        new Confirmation(
            Name: "d_zero_sector_trichotomy_marrakesh",
            Date: "2026-05-01",
            Machine: "ibm_marrakesh",
            JobId: "d7mjnjjaq2pc73a1pk4g",
            Observable: "⟨n⟩ on (q0, q2) reduced state, from Z-basis ZZ counts only",
            PredictedValue:
                "mean_n {truly_XXYY, soft_XYYX, hard_XXXY} = {1.0000, 1.2786, 0.8752} at γ_Z=0.1, t=0.8; " +
                "spread relative to truly: soft = +0.2786, hard = -0.1248.",
            MeasuredValue:
                "mean_n {truly, soft, hard} = {1.0098, 1.2673, 0.9082}; " +
                "spread relative to truly: soft = +0.2576, hard = -0.1016; max per-case deviation in mean_n = 0.0330 (~3% relative). All three classes separated by ≥0.10 in ⟨n⟩.",
            HardwareData: "data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json",
            ExperimentDoc: "data/ibm_soft_break_april2026/README.md",
            FrameworkPrimitive: "sector_populations (d_zero diagnostic)",
            Description:
                "Retrospective d=0-axis reading of the F77 trichotomy from existing Marrakesh data. The same three Hamiltonian classes the April-26 run discriminated via ⟨X₀Z₂⟩ at 9-Pauli-pair tomography are also discriminated by ⟨n⟩ from just the Z-basis ZZ counts. " +
                "Operational consequence: the trichotomy can be read with 1/9 the Pauli-measurement cost when the question is classification, not full M-spectrum reconstruction. See simulations/d_zero_trichotomy_marrakesh.py.",
            QubitPath: new[] { 48, 49, 50 }),

        // --- Torino calibration-era runs (Feb-Mar 2026), registered 2026-06-18 ---
        // These predate the systematic April-June campaign and have no recorded IBM
        // job_id (only data-file timestamps). The first two share the Feb-9 q52 run.
        new Confirmation(
            Name: "cpsi_quarter_crossing_torino_feb2026",
            Date: "2026-02-09",
            Machine: "ibm_torino",
            JobId: "tomography_ibm_torino_20260209_131521 (data-file timestamp; no IBM job_id recorded for the Torino calibration-era runs)",
            Observable: "CΨ(t) = Tr(ρ²)·L₁/(d−1) for |+⟩ under free decoherence; crossing time t* of the CΨ = ¼ boundary",
            PredictedValue:
                "t*/T₂* = 0.936 (generalized, r = T₂*/T₁ = 0.456); 0.858 in the pure-dephasing limit (x³+x = ½)",
            MeasuredValue:
                "t* = 114.7 μs, t*/T₂* = 1.041 (11.3% above the generalized prediction); CΨ(0) = 0.885, C∞ = 0.740. " +
                "Qubit 52: T1 = 221.2 μs, T2_echo = 298.2 μs, T2*(FID) = 110.7 μs",
            HardwareData: "data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json",
            ExperimentDoc: "experiments/IBM_QUANTUM_TOMOGRAPHY.md",
            FrameworkPrimitive: "F25 closed-form CΨ(t) + the CΨ = ¼ fold (K_fold dose)",
            Description:
                "The first CΨ = ¼ crossing ever seen on hardware, found in IBM Torino single-qubit calibration tomography " +
                "(Heron r2, q52, 2026-02-09), predating the systematic April-June 2026 campaign. A QUALITATIVE confirmation that " +
                "the ¼ fold is real on a physical device, not a precision match: t*/T₂* = 1.041 sits 11% above the generalized " +
                "prediction 0.936 because the crossing was extracted from a calibration run, not a purpose-built sweep. No IBM " +
                "job_id was recorded for the Torino-era runs; the data-file timestamp is the locator.",
            QubitPath: new[] { 52 }),

        new Confirmation(
            Name: "absorption_theorem_ratio_torino",
            Date: "2026-04-04",
            Machine: "ibm_torino",
            JobId: "tomography_ibm_torino_20260209_131521 (data-file timestamp; analysis 2026-04-04 of the 2026-02-09 q52 run; no IBM job_id recorded)",
            Observable: "Re(λ) / (−2γ⟨n_XY⟩), the excess coherence-decay ratio; effective ⟨n_XY⟩",
            PredictedValue:
                "ratio = 1 (Absorption Theorem: Re(λ) = −2γ⟨n_XY⟩, with ⟨n_XY⟩ = 1 for single-qubit coherence)",
            MeasuredValue:
                "ratio = 1.03 (3% deviation) on the T2* baseline: excess α = 0.006960 μs⁻¹ vs 2γ* = 0.006771 μs⁻¹, with a 2.8% " +
                "slow tail at the resolution limit. (The 6.37 figure quoted elsewhere used the wrong T2_echo baseline; the " +
                "dephasing-relevant T2* baseline gives 1.03.)",
            HardwareData: "data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json",
            ExperimentDoc: "experiments/IBM_ABSORPTION_THEOREM.md",
            FrameworkPrimitive: "Absorption Theorem Re(λ) = −2γ⟨n_XY⟩ (PROOF_ABSORPTION_THEOREM); simulations/ibm_absorption_theorem.py",
            Description:
                "Retrospective Absorption-Theorem reading of the same 2026-02-09 Torino q52 tomography run (analysis 2026-04-04). " +
                "The single-qubit coherence decays at the predicted Re(λ) = −2γ floor (⟨n_XY⟩ = 1), ratio 1.03. Predates the " +
                "systematic registry campaign; shares the Feb-9 q52 data with cpsi_quarter_crossing_torino_feb2026.",
            QubitPath: new[] { 52 }),

        new Confirmation(
            Name: "cpsi_quarter_crossing_torino_q80_mar2026",
            Date: "2026-03-18",
            Machine: "ibm_torino",
            JobId: "palindrome_ibm_torino_20260318_191348 (data-file timestamp; no IBM job_id recorded)",
            Observable: "CΨ = ¼ crossing time t* (8-point single-qubit tomography, q80 the \"permanent crosser\" of Run 3)",
            PredictedValue:
                "t* = 15.01 μs (from the same-day in-situ Ramsey T2* = 17.36 μs)",
            MeasuredValue:
                "t* = 15.29 μs, deviation 1.9% (0.28 μs). Qubit 80: T1 = 143.1-159 μs, same-day T2* = 17.36 μs (drifted 58% in " +
                "6 days from 11.0 μs, which is why the prediction uses the same-day Ramsey, not the stale calibration)",
            HardwareData: "data/ibm_run3_march2026/palindrome_ibm_torino_20260318_191348.json",
            ExperimentDoc: "experiments/IBM_RUN3_PALINDROME.md",
            FrameworkPrimitive: "F25 closed-form CΨ(t) + the CΨ = ¼ fold; same-day in-situ Ramsey T2*",
            Description:
                "The tightest Torino-era CΨ = ¼ confirmation (1.9%): Run 3 on q80, 2026-03-18, with a same-day in-situ Ramsey " +
                "T2* (17.36 μs). Predates the systematic April-June campaign but is precision-grade because the in-situ T2* " +
                "removed the calibration-drift error (q80 T2* had drifted 58% in 6 days). The earlier two Torino rows (Feb-9 " +
                "q52) are looser; this one pins the fold to 0.28 μs.",
            QubitPath: new[] { 80 }),
    };

    public static IReadOnlyList<Confirmation> All => _all;

    public static Confirmation? Lookup(string name) =>
        _all.FirstOrDefault(c => c.Name == name);

    public static IEnumerable<string> ListNames() => _all.Select(c => c.Name);

    public static IEnumerable<Confirmation> ByMachine(string machine) =>
        _all.Where(c => c.Machine == machine);

    /// <summary>Confirmations whose <see cref="Confirmation.QubitPath"/> matches
    /// <paramref name="path"/> as a sequence (ordered, length-equal). Skips
    /// entries with no documented path.</summary>
    public static IEnumerable<Confirmation> ByPath(IReadOnlyList<int> path) =>
        _all.Where(c => c.QubitPath != null && c.QubitPath.SequenceEqual(path));

    /// <summary>Combined filter: by machine + by exact qubit path. Useful for
    /// "what has been confirmed on path X of machine Y?" pre-submit lookups.</summary>
    public static IEnumerable<Confirmation> ByMachineAndPath(string machine, IReadOnlyList<int> path) =>
        ByPath(path).Where(c => c.Machine == machine);

    /// <summary>Confirmations whose <see cref="Confirmation.QubitPath"/> shares at
    /// least one qubit with <paramref name="path"/>. Looser than
    /// <see cref="ByPath"/>; useful for "what's been touched on any of these qubits?"
    /// queries.</summary>
    public static IEnumerable<Confirmation> ByPathOverlap(IReadOnlyList<int> path)
    {
        var s = new HashSet<int>(path);
        return _all.Where(c => c.QubitPath != null && c.QubitPath.Any(s.Contains));
    }
}
