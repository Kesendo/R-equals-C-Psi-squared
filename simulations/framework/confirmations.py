"""Hardware-confirmed framework predictions registry."""
from __future__ import annotations


class Confirmations:
    """Hardware-confirmed framework predictions, accessible by name.

    Each entry maps a name to: date, machine, job_id, observable,
    predicted_value, measured_value, hardware_data, experiment_doc,
    framework_primitive, description.

    Methods:
        lookup(name=None) → entry or full dict
        list_names() → list[str]
        by_machine(machine) → filtered dict
    """

    _ENTRIES = {
        'palindrome_trichotomy': {
            'date': '2026-04-26',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7mjnjjaq2pc73a1pk4g',
            'observable': '<X_0 Z_2>',
            'predicted_value': {
                'continuous_lindblad_gamma_Z_0.1': {'truly': 0.000, 'soft': -0.623, 'hard': +0.195,
                                                    'delta_soft_minus_truly': -0.623},
                'trotter_n3_gamma_Z_0.1':         {'truly': 0.000, 'soft': -0.723, 'hard': +0.327,
                                                    'delta_soft_minus_truly': -0.723},
            },
            'measured_value': {'truly': +0.011, 'soft': -0.711, 'hard': +0.205,
                               'delta_soft_minus_truly': -0.722},
            'hardware_data': 'data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json',
            'experiment_doc': 'experiments/MARRAKESH_THREE_LAYERS.md',
            'framework_primitive': 'classify_pauli_pair + predict_M_spectrum_pi2_odd + propagate_with_hardware_noise',
            'description': 'Super-operator palindrome trichotomy (truly/soft/hard) tomographically distinguishable on Heron r2 hardware at N=3. Hardware Δ(soft − truly) = -0.722 matches the Trotter n=3 prediction to 0.0014, NOT the continuous-Lindblad idealization (Δ = -0.623). Original interpretation that T1 amplifies the soft-break is REFUTED: T1 monotonically attenuates |Δ| (γ_T1=0.5 gives Δ=-0.44). The hardening is Trotter discretization at δt=0.267 with ‖H‖_op = 2.83·J, where ‖H·δt‖ ≈ 0.76 violates the small-step regime. See _marrakesh_t1_amplification_test.py.',
        },
        'f25_cusp_trajectory': {
            'date': '2026-04-26',
            'machine': 'ibm_kingston',
            'job_id': 'd7mu36lqrg3c738lnda0',
            'observable': 'CΨ(t) for Bell+',
            'predicted_value': 'F25: CΨ(t) = f·(1+f²)/6 with f = exp(-4·γ·t)',
            'measured_value': 'RMS residual 0.0097 vs in-situ γ_fit',
            'hardware_data': 'data/ibm_cusp_slowing_april2026/ (April 26 precision run)',
            'experiment_doc': 'experiments/CRITICAL_SLOWING_AT_THE_CUSP.md',
            'framework_primitive': 'F25 closed-form CΨ(t)',
            'description': 'Bell+ trajectory through CΨ=1/4 cusp confirmed point-by-point on Kingston (19 delay points, qubits 14-15).',
        },
        'f57_kdwell_gamma_invariance': {
            'date': '2026-04-16',
            'machine': 'ibm_kingston',
            'job_id': 'cusp_slowing_kingston_20260416',
            'observable': 'K_dwell / δ for Bell+',
            'predicted_value': 'F57: γ-independent, F25 prefactor 1.0801 for pure Z-dephasing',
            'measured_value': 'Pair A 0.6492, Pair B 0.6937 (6.3% spread despite 2.55× γ ratio)',
            'hardware_data': 'data/ibm_cusp_slowing_april2026/cusp_slowing_ibm_kingston_20260416_212042.json',
            'experiment_doc': 'experiments/CRITICAL_SLOWING_AT_THE_CUSP.md',
            'framework_primitive': 'K_dwell formula in F57',
            'description': 'γ-invariance of K_dwell at CΨ=1/4 boundary verified on Kingston with two qubit pairs at 2.55× different γ. Absolute prefactor 0.67 vs predicted 1.08 due to T1 amplitude damping.',
        },
        'bonding_mode_receiver': {
            'date': '2026-04-24',
            'machine': 'ibm_kingston',
            'job_id': 'multiple (see external pipeline)',
            'observable': 'MI(0, N-1) for bonding:2 vs alt-z-bits',
            'predicted_value': '4000-5500× over ENAQT in simulation N=5..13; ratio ≈ 1.4-3× hardware-realistic',
            'measured_value': 'bonding:2 / alt-z-bits = 2.80× on Kingston N=5',
            'hardware_data': 'external (AIEvolution.UI/experiments/ibm_quantum_tomography)',
            'experiment_doc': 'experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md',
            'framework_primitive': 'F67 bonding-mode + F71 chain-mirror symmetry',
            'description': 'Receiver-engineering bonding-mode advantage measured on Kingston. Largest engineering lever in the framework portfolio. Confirms F67 + F71 structure on hardware.',
        },
        'chiral_mirror_law': {
            'date': '2026-04-25',
            'machine': 'ibm_marrakesh',
            'job_id': 'k_partnership_marrakesh_20260425',
            'observable': 'Bloch components of K-partner pair states',
            'predicted_value': 'Chiral mirror identity from K=diag((-1)^l)',
            'measured_value': 'retrospective verification, identity holds',
            'hardware_data': 'data/ibm_k_partnership_april2026/k_partnership_marrakesh_20260425_140913.json',
            'experiment_doc': 'experiments/CHIRAL_MIRROR_HARDWARE_PREDICTION.md',
            'framework_primitive': 'K_full chiral conjugation (symmetry module)',
            'description': 'Chiral mirror law (K H K = -H spectral inversion) verified retrospectively on Marrakesh K-partnership data.',
        },
        'pi_protected_xiz_yzzy': {
            'date': '2026-04-26',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7n3013aq2pc73a2a18g',
            'observable': '<X_0 I Z_2>',
            'predicted_value': 'protected (≈0) for YZ+ZY soft Hamiltonian',
            'measured_value': '+0.13 to +0.04 (within noise band, never above ±0.13)',
            'hardware_data': 'external (raw JSON not in repo; results documented inline in experiment_doc)',
            'experiment_doc': 'review/EMERGING_QUESTIONS.md',
            'framework_primitive': 'pi_protected_observables',
            'description': 'First-time hardware measurement of a Π-protected observable on YZ+ZY soft Hamiltonian (EQ-030). Confirms framework primitive at hardware scale on a Hamiltonian not previously tested.',
        },
        'marrakesh_transverse_y_field_detection': {
            'date': '2026-04-29',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7ornigror3c73c0c6ug',
            'observable': 'Z⊗N-Mirror max-violation between ρ_a=|+−+⟩ and ρ_b=|−+−⟩',
            'predicted_value': 'h_y=0.05 → max_violation ≈ 0.18 (linear scaling 3.5·h_y); '
                               'h_x=0.05 → max_violation ≈ 0.004 (linear scaling 0.085·h_x); '
                               'clean stack → max_violation < 1e-3',
            'measured_value': 'max_violation = 0.182 (worst Pauli: Z,Z), RMS = 0.087. '
                              'Matches h_y=0.05 prediction exactly. h_y_eff ≈ 0.05.',
            'hardware_data': 'data/ibm_zn_mirror_april2026/zn_mirror_ibm_marrakesh_20260429_102824.json',
            'experiment_doc': 'data/ibm_zn_mirror_april2026/README.md',
            'framework_primitive': 'chain.zn_mirror_diagnostic + zn_mirror_state',
            'description': 'First hardware verification of the Z⊗N-Mirror diagnostic. Marrakesh shows effective transverse Y-field h_y_eff ≈ 0.05 at Hamiltonian level on path [48,49,50]. NOT a transverse X-field (which would give 40× smaller violation). The Y-vs-X asymmetry predicted by the framework (Y is bit_b-odd like Z-dephasing axis, mixes more strongly) is confirmed empirically. Worst-violating Pauli string is Z,Z, indicating state-preparation rotation between |+−+⟩ and |−+−⟩ runs that is consistent with single-site h_y rotation.',
        },
        'gamma_0_marrakesh_calibration': {
            'date': '2026-04-29',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7mjnjjaq2pc73a1pk4g',
            'observable': '<X_0 Z_2> for 3 Pauli-pair Hamiltonians (truly XX+YY, soft XY+YX, hard XX+XY)',
            'predicted_value': 'continuous-equivalent γ_Z ≈ 0.05 (continuous Lindblad fit); '
                               'Trotter-modeled γ_Z ≈ 0.1 (the same 0.1 the framework idealized prediction used)',
            'measured_value': 'best-fit γ_Z = 0.05 (sweep over [0.01, 0.15] with 71 points), '
                              'total residual² = 6.4e-4 across 3 Hamiltonians (continuous Lindblad, '
                              'no Trotter modeling). When Trotter n=3 is modeled, the same data fits '
                              'γ_Z = 0.1 exactly via Δ matching to 0.0014 (see _marrakesh_t1_amplification_test).',
            'hardware_data': 'data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json',
            'experiment_doc': 'data/ibm_soft_break_april2026/README.md',
            'framework_primitive': 'ChainSystem.propagate_with_hardware_noise + 2D fit (continuous Lindblad)',
            'description': 'Continuous-Lindblad fit of γ_Z against Marrakesh ⟨X₀Z₂⟩ data converges to 0.05, which absorbs the Trotter n=3 discretization correction into a lower effective γ_Z. The 2026-04-30 follow-up showed that a Trotter-modeled fit returns γ_Z = 0.1 with Δ-matching to 0.0014, meaning the original framework-idealized γ_Z = 0.1 was correct (not 2× too high). The two values 0.05 (continuous) and 0.1 (Trotter) are the same data through two physics models. T1 contributes negligibly in either model.',
        },
        'lebensader_skeleton_trace_decoupling': {
            'date': '2026-04-26',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7n3013aq2pc73a2a18g + d7n3eqqt99kc73d34qtg',
            'observable': 'Π-protected counts (skeleton) + θ-trajectory tails (trace)',
            'predicted_value': 'Skeleton/trace decouple at N≥4, co-occur at N=3. T1 makes drop measurable: bond-flipped Z-free pairs (XY+YX, IY+YI) preserve skeleton (drop≤1) and trace (long θ-tail); Z-containing soft pairs (YZ+ZY, XZ+ZX) collapse both (drop≈28-29, no tail).',
            'measured_value': 'drop=28 for YZ+ZY confirmed on Marrakesh. Pearson(drop, Δ∫θ)=+0.85. Bures velocity gives no third discriminator.',
            'hardware_data': 'external (raw JSON not in repo; tables of <X_0 I Z_2> per t and basis documented inline in experiment_doc)',
            'experiment_doc': 'review/EMERGING_QUESTIONS.md',
            'framework_primitive': 'cockpit_panel — composes pi_protected_observables + θ-trajectory',
            'description': 'Lebensader as Stromkabel (EQ-030 closure): skeleton (Π-protected algebraic count) and trace (θ-geometric tail) are NOT two discriminators but one bridge held together by Π·L·Π⁻¹ + L + 2Σγ·I = 0. Hardware-confirmed across 3 of 4 bond-flipped Z-free corners; Bures velocity confirmed null as third axis. ChainSystem.cockpit_panel gives this in one call.',
        },
        'd_zero_sector_trichotomy_marrakesh': {
            'date': '2026-05-01',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7mjnjjaq2pc73a1pk4g',
            'observable': '⟨n⟩ on (q0, q2) reduced state, from Z-basis ZZ counts only',
            'predicted_value': {
                'truly_XXYY': {'p_n': [0.3085, 0.3830, 0.3085], 'mean_n': 1.0000},
                'soft_XYYX':  {'p_n': [0.1692, 0.3830, 0.4478], 'mean_n': 1.2786},
                'hard_XXXY':  {'p_n': [0.3124, 0.5000, 0.1876], 'mean_n': 0.8752},
                'spread_relative_to_truly': {'soft - truly': +0.2786, 'hard - truly': -0.1248},
            },
            'measured_value': {
                'truly_XXYY': {'p_n': [0.2859, 0.4185, 0.2957], 'mean_n': 1.0098},
                'soft_XYYX':  {'p_n': [0.1672, 0.3982, 0.4346], 'mean_n': 1.2673},
                'hard_XXXY':  {'p_n': [0.3040, 0.4839, 0.2122], 'mean_n': 0.9082},
                'spread_relative_to_truly': {'soft - truly': +0.2576, 'hard - truly': -0.1016},
                'max_per_case_deviation_in_mean_n': 0.0330,
            },
            'hardware_data': 'data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json',
            'experiment_doc': 'data/ibm_soft_break_april2026/README.md',
            'framework_primitive': 'sector_populations (d_zero diagnostic)',
            'description': 'Retrospective d=0-axis reading of the F77 trichotomy from existing Marrakesh data. The same three Hamiltonian classes that the original April-26 run discriminated via ⟨X₀Z₂⟩ at 9-Pauli-pair tomography are also discriminated by ⟨n⟩ from just the Z-basis ZZ counts. Framework predicts ⟨n⟩ values of {1.000, 1.279, 0.875} for {truly, soft, hard} on (q0,q2) reduced state of |+−+⟩ under L_H + γ_Z·L_dephase at γ_Z=0.1, t=0.8. Hardware reads {1.010, 1.267, 0.908}, max per-case deviation 0.033 (~3% relative). All three classes separated by ≥0.10 in ⟨n⟩, well above hardware noise. Operational consequence: the trichotomy can be read with 1/9 the Pauli-measurement cost when the question is classification, not full M-spectrum reconstruction. See simulations/d_zero_trichotomy_marrakesh.py for the analysis.',
        },
        'block_cpsi_saturation_kingston_may2026': {
            'date': '2026-05-08',
            'machine': 'ibm_kingston',
            'job_id': 'd7ulfjdpa59c73b4rttg',
            'observable': 'C_block(t) = Σ |ρ_{ab}|² over (popcount-0, popcount-1) on 2-qubit ρ from Dicke initial (|D_0⟩+|D_1⟩)/√2',
            'predicted_value': {
                'theorem_1_t0': 0.25,
                'theorem_3_trajectory_form': '(1/4) * exp(-4 * gamma * t), with gamma = 1/(2*T2) for pure Z-dephasing',
                'gamma_expected_from_T2_min_480us': 1.042e-3,
            },
            'measured_value': {
                't_grid_us': [0.0, 120.0, 240.0, 360.0, 480.0],
                'C_block_n0': [0.2205, 0.1023, 0.0458, 0.0166, 0.0074],
                'fraction_of_quarter_at_t0': 0.882,
                'gamma_fit_per_us': 1.795e-3,
                'R_squared': 0.9977,
                'T2_eff_from_fit_us': 278.5,
                'T2_calibration_min_us': 480.0,
                'hardware_to_T2_speedup': 1.72,
            },
            'hardware_data': 'data/ibm_block_cpsi_saturation_may2026/block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json',
            'experiment_doc': 'experiments/IBM_BLOCK_CPSI_SATURATION.md',
            'framework_primitive': 'BlockCoherenceContent.Compute + BlockCpsiClosedForm.At + IbmBlockCpsiHardwareTable (typed); state-level Theorem 1 + Theorem 3 trajectory of PROOF_BLOCK_CPSI_QUARTER',
            'description': 'Direct hardware confirmation of the universal-Mandelbrot-cardioid 1/4 ceiling on the (popcount-0, popcount-1) coherence block of a 2-qubit ρ. Initial state (|D_0⟩+|D_1⟩)/√2 saturates Theorem 1 at 88.2% of 1/4 on Kingston q13–q14 (state-prep + tomography fidelity floor). The Theorem-3 closed-form trajectory C_block(t) = (1/4)·exp(-4γt) fits the 5 t-points (0, 120, 240, 360, 480 μs) with R² = 0.9977 (pristine log-linear decay). Fitted γ_fit = 1.795e-3 μs⁻¹ corresponds to T2_eff = 278.5 μs from the C_block decay alone, vs T2_min = 480 μs from the morning calibration; the 1.72× speedup quantifies the gate-noise + readout contribution beyond pure single-qubit T2 dephasing. First Tier-2-Verified hardware anchor for the Mandelbrot 1/4 boundary on a state engineered to sit on it.',
        },
        'f83_pi2_class_signature_marrakesh': {
            'date': '2026-04-30',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7pol1e7g7gs73cf7j90',
            'observable': 'Per-category 2-qubit Pauli expectations (q0=4, q2=6) for 4 F77 Π²-classes',
            'predicted_value': 'Trotter+γ_Z patterns at the path-fit γ_Z_eff: pi2_odd_pure ⟨X₀Z₂⟩, pi2_even_nontruly ⟨X₀X₂⟩, mixed ⟨Z₀X₂⟩, truly ⟨Y₀Z₂⟩. F83 anti-fractions = {1/2, 0, 1/6} structurally (closed form on H letters).',
            'measured_value': '4-category discrimination CONFIRMED via unique-fingerprint Paulis: pi2_odd_pure ⟨X₀Z₂⟩=-0.849, pi2_even_nontruly ⟨X₀X₂⟩=+0.919, mixed ⟨Z₀X₂⟩=-0.721, truly ⟨Y₀Z₂⟩=+0.670. All separations >>10σ at 4096 shots. Path-fit γ_Z_eff = 0.050 (vs 0.120 on April 26 [48,49,50]). At γ_Z_eff, pi2_odd_pure RMS=0.039 and pi2_even_nontruly RMS=0.029 over 7 Paulis (excellent quantitative match). Truly RMS=0.188 and mixed RMS=0.163 retain residuals: ⟨Y,Z⟩-⟨Z,Y⟩ asymmetry (per-qubit T2 inequality Q4=184μs vs Q6=151μs); truly ⟨Z,Z⟩ damped 60% (candidate F82/F84 amplitude-damping signature, open for follow-up).',
            'hardware_data': 'data/ibm_f83_signature_april2026/f83_signature_ibm_marrakesh_20260430_190035.json',
            'experiment_doc': 'data/ibm_f83_signature_april2026/README.md',
            'framework_primitive': 'classify_pauli_pair + predict_pi_decomposition (F83 anti-fraction closed form)',
            'description': 'F83 4-Hamiltonian Π²-class discrimination test on path [4,5,6] (top-ranked by 2026-04-30T16:25Z calibration, score 0.0162 of 223). Each of the 4 F77 classes shows a unique-fingerprint Pauli observable separating it from the other 3 at >>10σ. The two pure-Π²-class categories (pi2_odd_pure XY+YX, pi2_even_nontruly YZ+ZY) match the Trotter+γ_Z model quantitatively at path-specific γ_Z_eff = 0.050 with RMS < 0.04. γ_Z_eff is path-dependent (0.05 vs 0.12 between [4,5,6] and [48,49,50]), reflecting that effective dephasing absorbs Trotter discretization, coherent gate errors, and crosstalk in path-specific proportions. The truly XX+YY and mixed XY+YZ categories show larger residuals attributable to per-qubit T2 inequality (asymmetric Y-basis observables) and candidate amplitude-damping signature (truly ⟨Z,Z⟩ 60% damped, consistent with F82/F84 σ⁻ breaking ⟨Z⟩-conservation for truly while pi2_odd is unaffected because its Hamiltonian doesn\'t conserve ⟨Z⟩ to begin with). See _f83_gamma_z_sweep.py and _f83_hy_field_check.py for the analysis (h_y=0.05 from zn_mirror was tested as alternative hypothesis and rejected).',
        },
    }

    @classmethod
    def lookup(cls, name=None):
        """Return all confirmations or one entry by name."""
        if name is None:
            return dict(cls._ENTRIES)
        if name in cls._ENTRIES:
            return dict(cls._ENTRIES[name])
        raise KeyError(
            f"No confirmation named {name!r}. Available: {list(cls._ENTRIES.keys())}"
        )

    @classmethod
    def list_names(cls):
        return list(cls._ENTRIES.keys())

    @classmethod
    def by_machine(cls, machine):
        return {k: dict(v) for k, v in cls._ENTRIES.items() if v.get('machine') == machine}
