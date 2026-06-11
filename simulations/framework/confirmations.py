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
            'framework_primitive': 'cockpit_panel: composes pi_protected_observables + θ-trajectory',
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
            'cross_machine': 'Reproduced on ibm_kingston 2026-05-05 (job d7sqjpiudops73976960, path [43,56,63], data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json): pi2_odd_pure ⟨X₀Z₂⟩=-0.774, pi2_even_nontruly ⟨X₀X₂⟩=+0.843, mixed ⟨Z₀X₂⟩=-0.724, truly ⟨Y₀Z₂⟩=+0.642. Same 4-way fingerprint separation, so the F83 discriminator is machine-independent. This same Kingston run is registered in the C# ConfirmationsRegistry as regime_uniformity_kingston_uniform_quantum, there read through the F88b-Lens Π²-odd-memory lens (truly=0.0022, soft=0.7409, pi2_even=0.0046, mixed=0.6220).',
            'description': 'F83 4-Hamiltonian Π²-class discrimination test on path [4,5,6] (top-ranked by 2026-04-30T16:25Z calibration, score 0.0162 of 223). Each of the 4 F77 classes shows a unique-fingerprint Pauli observable separating it from the other 3 at >>10σ. The two pure-Π²-class categories (pi2_odd_pure XY+YX, pi2_even_nontruly YZ+ZY) match the Trotter+γ_Z model quantitatively at path-specific γ_Z_eff = 0.050 with RMS < 0.04. γ_Z_eff is path-dependent (0.05 vs 0.12 between [4,5,6] and [48,49,50]), reflecting that effective dephasing absorbs Trotter discretization, coherent gate errors, and crosstalk in path-specific proportions. The truly XX+YY and mixed XY+YZ categories show larger residuals attributable to per-qubit T2 inequality (asymmetric Y-basis observables) and candidate amplitude-damping signature (truly ⟨Z,Z⟩ 60% damped, consistent with F82/F84 σ⁻ breaking ⟨Z⟩-conservation for truly while pi2_odd is unaffected because its Hamiltonian doesn\'t conserve ⟨Z⟩ to begin with). See _f83_gamma_z_sweep.py and _f83_hy_field_check.py for the analysis (h_y=0.05 from zn_mirror was tested as alternative hypothesis and rejected).',
        },
        'f95_angle_steering_kingston_may2026': {
            'date': '2026-05-16',
            'machine': 'ibm_kingston',
            'job_id': 'bxyj5yd4j + bzklqwt7f',
            'observable': 'arg(CΨ_com(t)) at first real-axis magnitude-minimum crossing under deliberate RZ phase injection at rate Ω (two Ω values: 0.13 and 0.25 rad/μs; per-Ω mapping in measured_value)',
            'predicted_value': {
                'framework_claim': 'F95 complex-CΨ angle is actively steerable via per-chunk RZ(Ω·Δt) injection on the bra-ket Hermitian off-diagonal',
                'model': 'Lindblad+RZ with in-situ γ_per_us from T2 calibration',
                'pair_A_mid_qubits': [82, 83],
                'pair_A_mid_gamma_per_us': 3.3632e-3,
                'pair_A_mid_T2_min_us': 148.67,
                'pair_B_high_qubits': [13, 14],
                'pair_B_high_gamma_per_us': 1.8964e-3,
                'pair_B_high_T2_min_us': 263.66,
            },
            'measured_value': {
                'pair_A_mid_omega_0.13': {'t_cross_us': 1.395, 'arg_measured_deg': -84.70, 'arg_predicted_deg': -100.39, 'residual_deg': 15.69},
                'pair_A_mid_omega_0.25': {'t_cross_us': 1.242, 'arg_measured_deg': -98.64, 'arg_predicted_deg': -107.79, 'residual_deg': 9.15},
                'pair_B_high_omega_0.13': 'no magnitude crossing in delay window (higher T2 + lower Ω left trajectory above threshold)',
                'pair_B_high_omega_0.25': {'t_cross_us': 2.814, 'arg_measured_deg': -123.50, 'arg_predicted_deg': -130.31, 'residual_deg': 6.81},
                'residual_shrinks_with_omega': 'on pair A_mid, residual 15.7° → 9.2° as Ω doubles from 0.13 → 0.25, consistent with roughly constant natural-Kingston drift becoming smaller fraction of driven rotation',
                'total_billed_qpu_minutes': 14,
            },
            'hardware_data': 'data/ibm_f95_angle_steering_may2026/cusp_complex_phase_hardware_ibm_kingston_omega0.{130,250}_20260516_*.json (+ matching PNGs and partial-run JSONs; see README.md in that directory for details)',
            'experiment_doc': 'experiments/CPSI_COMPLEX_PLANE.md (conceptual predecessor; F95 closed form derived 2026-05-16 same day)',
            'framework_primitive': 'F95AngleAtQuadraticZeroPi2Inheritance (θ(c;b) = arctan(√(c/b² − 1))) + F25 Lindblad CΨ trajectory + complex CΨ_com signed-sum-of-off-diagonals extension',
            'description': 'First hardware verification that the complex-CΨ angle predicted by F95 is actively steerable on IBM Kingston. Two Ω values (0.13 and 0.25 rad/μs) tested on two qubit pairs (A_mid mid-T2 [82,83], B_high high-T2 [13,14]) in 14 minutes total billed QPU. Three of four measurement conditions produced detectable real-axis crossings; predicted arg matches measured to within 16° across all three (best 6.8°). Confirms three open questions from CPSI_COMPLEX_PLANE.md simultaneously: (Q1) Δφ = Ω·t reproducible: crossings occur at predicted angle modulo bounded natural drift; (Q2) linear-in-Ω scaling holds: residual shrinks ~1.7× when Ω doubles, consistent with constant natural drift contribution; (Q3) active steering works: the same Lindblad+RZ model fits both Ω values on the same pair without retuning. The bra↔ket Hermitian off-diagonal of ρ on the prepared Bell-like initial state carries a well-defined complex phase, addressable on Heron r2 by per-chunk RZ injection on the steering qubit. Hardware substrate-of-evidence: F95 is the universal polynomial-foundation algebra of the angle that appears off the d=0 mirror; this run shows that algebra is operationally controllable, not merely descriptive.',
        },
        'regime_uniformity_kingston_uniform_quantum': {
            'date': '2026-05-05',
            'machine': 'ibm_kingston',
            'job_id': 'd7sqjpiudops73976960',
            'observable': 'F88b-Lens Π²-odd-memory across F87 trichotomy on uniform-quantum chain (path [43,56,63])',
            'predicted_value': 'regime-uniformity hypothesis: a uniform-quantum chain truly-baseline should sit near a uniform-classical one, both well below the regime-mixed value 0.0297 (Marrakesh [0,1,2] framework_snapshots); soft-pumping should stay hardware-substrate-independent at ~0.74.',
            'measured_value': 'truly = 0.0022, soft = 0.7409, pi2_even_nontruly = 0.0046, mixed = 0.6220. truly is 13.5× lower than the regime-mixed Marrakesh [0,1,2] (0.0297) and only 1.69× higher than uniform-classical Marrakesh [48,49,50] (0.0013). soft sits 3.1% from the Marrakesh value 0.7646: substrate-independent within shot noise.',
            'hardware_data': 'data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json',
            'experiment_doc': 'simulations/_f88b_lens_ibm_kingston_uniform_quantum.py',
            'framework_primitive': 'F88b-Lens (kernel projection + Π²-odd Pauli enumeration on reduced 2-qubit ρ)',
            'cross_machine': 'Same Kingston run is also the F83 cross-machine reproduction (see f83_pi2_class_signature_marrakesh cross_machine field), there read through the F83 anti-fraction fingerprint lens instead of this F88b Π²-odd-memory lens.',
            'description': 'First F87 trichotomy hardware test on a uniform-quantum CZ-coupled triple. Path [43, 56, 63] on Kingston, all three qubits PulseStable across the 91-day biography window. Confirms regime-uniformity: BOTH uniform-classical (Marrakesh [48,49,50] truly = 0.0013) AND uniform-quantum (Kingston [43,56,63] truly = 0.0022) give clean truly-baselines, while regime-mixed (Marrakesh [0,1,2] truly = 0.0297) is an order of magnitude dirtier. Three findings per single 39-second billed run: regime-uniformity is the cause of the truly-baseline gap; F87 trichotomy hardware-confirmed on a second backend (Kingston) in addition to Marrakesh; soft Π²-odd-pumping substrate-independent across two Heron-r2 chips.',
        },
        'gamma0_off_the_lever_kingston_may2026': {
            'date': '2026-05-29',
            'machine': 'ibm_kingston',
            'job_id': 'd8ce8l38ch0s738uorjg',
            'observable': 'Transfer T(t) = P(qubit 1 excited) from |10⟩ under exchange H = J·(XX+YY)/2 + idle Z-dephasing, scanned over J (the only lever); q13-q14',
            'predicted_value': 'Q = J/γ₀ with the typed carrier γ₀ = UniversalCarrierClaim.DefaultGammaZero = 0.05. Coherent swing (T overshoots ½) for Q > 1 (J > γ₀); overdamped creep (T → ½) for Q ≤ 1; critical damping exactly at J = γ₀ = 0.05 (Q = 1, max T = ½).',
            'measured_value': 'q13-q14, J ∈ {0.05, 0.1, 0.2, 0.4} (Q ≈ {1, 2, 4, 8}): max T = {0.335, 0.563, 0.703, 0.779}. Coherent swing (overshoot past ½) present for J ≥ 0.1, absent for J = 0.05 (creeps to 0.335). First transfer peak at step ∝ 1/J. Swing dies between J = 0.1 and J = 0.05, so γ₀ ≈ 0.05–0.1 read off the lever.',
            'hardware_data': 'external (AIEvolution.UI/experiments/ibm_quantum_tomography/results/q_jscan_hardware_ibm_kingston_20260529_025757.json)',
            'experiment_doc': 'experiments/GAMMA0_IS_ALWAYS_THERE.md',
            'framework_primitive': 'UniversalCarrierClaim.DefaultGammaZero (γ₀ = 0.05) + Q = J/γ₀ regime axis; run_q_jscan.py (external pipeline)',
            'description': 'First direct read-off of the carrier γ₀ from its only lever J on hardware. The H-clock frequency tracks J (first transfer peak ∝ 1/J), and the coherent swing dies where J crosses γ₀ (critical damping, T = ½). Confirms the typed UniversalCarrierClaim.DefaultGammaZero = 0.05 on Kingston q13-q14, made readable via the swing-death threshold even though γ₀ is invisible head-on (Inside-Observability: only Q = J/γ₀ accessible from inside). Sits on the same Q = J/γ₀ axis as the F86 Q-peak (Q ≈ 1.5).',
        },
        'ibm_ep_onset_may2026': {
            'date': '2026-05-31',
            'machine': 'ibm_kingston',
            'job_id': 'd8dr7dfd0j8c73f4man0 (Part A, flow to 1/3) + d8drjbfd0j8c73f4mobg (Part B, EP onset)',
            'observable': 'Per-site populations ⟨n_l⟩ of a single-excitation walk on a 3-site chain (q13-q14-q15, Z-basis, no tomography); Part B revival = max ⟨n_0⟩ for t ≥ 2 μs under injected random-Z-twirl dephasing, Q = J/γ swept over {0.5, 1, 1.5, 2.5, 5, 20}',
            'predicted_value': {
                'ep_onset_reading': 'revival pinned at the 1/N = 1/3 equipartition floor below Q_EP ≈ 1.5 (overdamped, two real decay channels, no memory), lifting off above it (the rotation born at the F86 exceptional point)',
                'twirl_simulate_revival': '0.31 → 0.84 across the same Q scan (K=16 exact-statevector twirl validation)',
                'part_A_flow_target': 'per-site populations converge to 1/N = 1/3 at late t (the post-EP flow fixed point)',
            },
            'measured_value': {
                'Q_grid': [0.5, 1.0, 1.5, 2.5, 5.0, 20.0],
                'revival': [0.30, 0.36, 0.34, 0.49, 0.56, 0.70],
                'floor_below_ep': '~1/3 for Q ≤ 1.5 (0.30, 0.36, 0.34)',
                'liftoff_above_ep': '0.49 → 0.56 → 0.70 as Q crosses 2.5 → 5 → 20',
                'part_A_late_t_populations': '0.34 / 0.43 / 0.34 at 20 μs (converging to 1/3)',
                'high_Q_suppression': '0.70 measured vs 0.84 simulate at Q=20 (two-qubit gate error, ~160 RZZ gates by 20 μs at ~0.5%); the floor and the onset are clean',
            },
            'hardware_data': 'data/ibm_ep_onset_may2026/ (Part A ep_onset_hardware_ibm_kingston_20260531_060943.json, Part B ep_onset_hardware_ep_ibm_kingston_20260531_064022.json + same-day simulate JSONs)',
            'experiment_doc': 'experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md',
            'framework_primitive': 'ExceptionalPointClock (decay pinning at 4γ₀, F95 rotation angle, eigenvector overlap min(x,1/x)) + EpField hardware node (inspect --axis ep); F86 Q_EP',
            'description': 'The F86 exceptional point watched switching the memory on, on a real chip, populations only. Part A (job d8dr7dfd0j8c73f4man0) at the chip\'s natural Q ≫ 1: the excitation sloshes 0 → 2 → 1 → 0 with ~3 μs period (the reborn memory), the site-0 revival fades 0.84 → 0.43 over 15 μs (the forgetting), and the populations converge to 1/3 = 1/N at 20 μs (the flow target). Part B (job d8drjbfd0j8c73f4mobg) injects dephasing via a random-Z twirl (K=16 instances; the RZ gates are virtual on IBM, so the injection is error-free) to push Q = J/γ down through Q_EP ≈ 1.5: the revival sits on the equipartition floor (~1/3) for Q ≤ 1.5 and lifts off as Q crosses 1.5 → 2.5 (0.34 → 0.49 → 0.56 → 0.70). The rotation born at the EP, hardware-anchored. The revival decay envelope is gate-cost-limited (Trotterization, ~9 μs), not T2-limited (~200 μs); only the rate is gate cost, the floor and the onset are physics. This table is the hardware node of EpField (compute/RCPsiSquared.Diagnostics/Foundation/EpField.cs) and the overlay in simulations/ep_transition.py.',
        },
        'f120_moment_tower_kingston_june2026': {
            'date': '2026-06-11',
            'machine': 'ibm_kingston',
            'job_id': 'd8l6c7rqv2lc73863acg (Arm A) + d8l6c832d42s73cb16a0 (Arm B) + d8l6h03nn5bs738rmrug (standard-T1 arbiter)',
            'observable': 'Energy-moment slopes d/dt⟨H_p^j⟩ from the maximally mixed state (8-basis-state average) under pure idle, H_p = X₀+X₀Z₁+0.7·X₁X₂ (girth-2 witness, H_p² = 2.49·I+2·Z₁+1.4·XXX exact); qubits q149/q13/q9, NO two-qubit gates, dynamical decoupling disabled, τ ∈ {0..150} μs, two arms permuting the middle qubit',
            'predicted_value': {
                'rung_1_null': 'slope⟨H_p⟩ = 0 at all τ (t₁ ≡ 0; evolution-blind, robust against all Z-flavored idle parasitics)',
                'rung_2_fires': 'slope⟨H_p²⟩ = 2·Δγ_mid (t₂ = [0,16,0] fires at the middle site only)',
                'row_identity': '⟨H_p²⟩ = 2.49 + 2⟨Z₁⟩ + 1.4⟨XXX⟩ row by row; ⟨XXX⟩ stays 0',
                'site_tracking': 'the firing follows the middle-qubit identity across arms; per-qubit pump slopes arm-independent',
                'rates_from_calibration': 'Δγ_l = 1/T1_l (the textbook noise-model assumption; this is the layer the chip declined)',
            },
            'measured_value': {
                'rung_1_null': 'slope⟨H⟩ = +2.4e-4/μs (z = +1.47, Arm A) and −6.8e-6/μs (z = −0.04, Arm B): the double null HELD',
                'row_identity': 'exact in every measured row (e.g. τ=100 Arm A: 2.49+2·0.4147+1.4·0.0186 = 3.345 vs 3.3455); |⟨XXX⟩| ≤ 0.02 throughout',
                'per_qubit_pump_slopes_per_us': {'q149': [2.327e-3, 2.193e-3], 'q13': [3.029e-3, 3.090e-3], 'q9': [5.794e-3, 5.779e-3]},
                'cross_arm_reproducibility': 'q9 0.3%, q13 1.9%, q149 5.7% (different chain roles, different jobs)',
                'arbiter_T1_us': {'q149': 424.6, 'q13': 430.3, 'q9': 99.9},
                'model_test': 'pump/Γ ratios 0.93-0.99 (q149, textbook), 1.30-1.33 (q13, VIOLATES pump ≤ Γ at 4-6σ: non-unital content beyond two-level amplitude damping, telegraph jump in Arm A), 0.58 (q9, deficit; its |1⟩-decay is itself non-exponential, 13.4e-3 → 7.7e-3/μs)',
            },
            'hardware_data': 'data/ibm_moment_tower_june2026/ (main + arbiter JSONs + the 06:33Z calibration snapshot)',
            'experiment_doc': 'experiments/F120_MOMENT_TOWER_KINGSTON.md',
            'framework_primitive': 'moment_tower + predict_pump_slope + f113_bridge_asymmetry_from_slope (diagnostics/f120_moment_tower, called by the pipeline script at startup); MomentTowerPumpChannelClaim (C#); PROOF_MOMENT_TOWER_PUMP_CHANNEL',
            'description': 'F120 first hardware reading, the cleanest protocol we ever sent to a QPU (not one entangling gate). '
                           'The structural law confirmed: the first energy-moment rung stays silent while the second fires as exactly '
                           'twice the middle qubit\'s pump curve, the girth read from hardware is 2, and the per-qubit pump rates '
                           'reproduce across arms to 0.3-5.7%. The rate layer returned a finding instead of a number: the pump slope '
                           'reads Tr(Z_l·D(I)), the device\'s TRUE non-unital pump vector, and it disagrees with calibrated amplitude '
                           'damping on two of three qubits. The unplanned bonus is a per-qubit noise-model test: pump ≤ Γ_T1 holds for '
                           'every two-level Lindblad model, q149 passes it, q9 dodges it (non-exponential decay, 1/T1 ill-defined), '
                           'and q13 violates it reproducibly, pointing at leakage/TLS physics beyond the textbook channel. The '
                           'hardness rung of a programmed Hamiltonian is now a quantity a chip measures about itself by decaying.',
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
