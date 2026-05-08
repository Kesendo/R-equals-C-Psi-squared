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
            Observable: "F88-Lens Π²-odd-memory across F87 trichotomy on uniform-quantum chain",
            PredictedValue:
                "regime-uniformity hypothesis: a uniform-quantum chain's truly-baseline should be near a uniform-classical chain's, both well below the regime-mixed value 0.0297 (Marrakesh [0,1,2] framework_snapshots). " +
                "soft-pumping should remain hardware-substrate-independent at ~0.74.",
            MeasuredValue:
                "truly = 0.0022, soft = 0.7409, pi2_even_nontruly = 0.0046, mixed = 0.6220. " +
                "truly is 13.5× lower than the regime-mixed Marrakesh [0,1,2] (0.0297) and only 1.69× higher than uniform-classical Marrakesh [48,49,50] (0.0013); the 22.8× gap between regime-mixed and uniform-classical from the original framework_snapshots vs soft_break comparison sits an order of magnitude wider than the within-uniform variation. " +
                "soft is 3.1% from Marrakesh's 0.7646: substrate-independent within shot noise.",
            HardwareData: "data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json",
            ExperimentDoc: "simulations/_f88_lens_ibm_kingston_uniform_quantum.py",
            FrameworkPrimitive: "F88-Lens (kernel projection + Π²-odd Pauli enumeration on reduced 2-qubit ρ)",
            Description:
                "First F87 trichotomy hardware test on a uniform-quantum CZ-coupled triple. Path [43, 56, 63] on Kingston, all three qubits classified PulseStable across the 91-day biography window (r mean 0.10 / 0.09 / 0.10, all crossing > 95%). " +
                "Result confirms regime-uniformity hypothesis: BOTH uniform-classical (Marrakesh [48,49,50] truly = 0.0013) AND uniform-quantum (Kingston [43,56,63] truly = 0.0022) give clean truly-baselines, while regime-mixed (Marrakesh [0,1,2] truly = 0.0297) is an order of magnitude dirtier. " +
                "Three findings per single 39-second billed run: regime-uniformity is the cause of the truly-baseline gap (not which side of the boundary); F87 trichotomy hardware-confirmed on a second backend (Kingston) in addition to Marrakesh; soft Π²-odd-pumping confirmed substrate-independent across two Heron-r2 chips.",
            QubitPath: new[] { 43, 56, 63 }),
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
