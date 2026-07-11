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
            'description': 'Super-operator palindrome trichotomy (truly/soft/hard) tomographically distinguishable on Heron r2 hardware at N=3. Hardware Δ(soft − truly) = -0.722 matches the Trotter n=3 prediction to 0.0014, NOT the continuous-Lindblad idealization (Δ = -0.623). Original interpretation that T1 amplifies the soft-break is REFUTED: T1 monotonically attenuates |Δ| (γ_T1=0.5 gives Δ=-0.44). The hardening is Trotter discretization at δt=0.267 with ‖H‖_op = 2.83·J, where ‖H·δt‖ ≈ 0.76 violates the small-step regime. See marrakesh_t1_amplification_test.py.',
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
                              'γ_Z = 0.1 exactly via Δ matching to 0.0014 (see marrakesh_t1_amplification_test).',
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
            'description': 'F83 4-Hamiltonian Π²-class discrimination test on path [4,5,6] (top-ranked by 2026-04-30T16:25Z calibration, score 0.0162 of 223). Each of the 4 F77 classes shows a unique-fingerprint Pauli observable separating it from the other 3 at >>10σ. The two pure-Π²-class categories (pi2_odd_pure XY+YX, pi2_even_nontruly YZ+ZY) match the Trotter+γ_Z model quantitatively at path-specific γ_Z_eff = 0.050 with RMS < 0.04. γ_Z_eff is path-dependent (0.05 vs 0.12 between [4,5,6] and [48,49,50]), reflecting that effective dephasing absorbs Trotter discretization, coherent gate errors, and crosstalk in path-specific proportions. The truly XX+YY and mixed XY+YZ categories show larger residuals attributable to per-qubit T2 inequality (asymmetric Y-basis observables) and candidate amplitude-damping signature (truly ⟨Z,Z⟩ 60% damped, consistent with F82/F84 σ⁻ breaking ⟨Z⟩-conservation for truly while pi2_odd is unaffected because its Hamiltonian doesn\'t conserve ⟨Z⟩ to begin with). See f83_gamma_z_sweep.py and f83_hy_field_check.py for the analysis (h_y=0.05 from zn_mirror was tested as alternative hypothesis and rejected).',
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
            'experiment_doc': 'simulations/f88b_lens_ibm_kingston_uniform_quantum.py',
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
                'se_walk_handover_reading': 'revival pinned at the 1/N = 1/3 equipartition floor below Q ≈ 1.5 (overdamped, no memory return), lifting off above it (the single-excitation walk\'s critical-damping transition; real, measured). Whether that SE transition is itself a genuine defective EP is under a SEPARATE open review (the coherence-horizon √-EP cluster, inspect --root horizon); not asserted here.',
                'twirl_simulate_revival': '0.31 → 0.84 across the same Q scan (K=16 exact-statevector twirl validation)',
                'part_A_flow_target': 'per-site populations converge to 1/N = 1/3 at late t (the post-handover flow fixed point)',
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
            'framework_primitive': 'ExceptionalPointClock (the toy 2×2 reduction: decay pinning at 4γ₀, F95 rotation angle, eigenvector overlap min(x,1/x)) + EpField hardware node (inspect --axis ep); Q ≈ 1.5 handover marker',
            'description': 'The single-excitation-walk overdamped→revival handover watched switching the memory on, on a real chip, populations only. Part A (job d8dr7dfd0j8c73f4man0) at the chip\'s natural Q ≫ 1: the excitation sloshes 0 → 2 → 1 → 0 with ~3 μs period (the reborn memory), the site-0 revival fades 0.84 → 0.43 over 15 μs (the forgetting), and the populations converge to 1/3 = 1/N at 20 μs (the flow target). Part B (job d8drjbfd0j8c73f4mobg) injects dephasing via a random-Z twirl (K=16 instances; the RZ gates are virtual on IBM, so the injection is error-free) to push Q = J/γ down through Q ≈ 1.5: the revival sits on the equipartition floor (~1/3) for Q ≤ 1.5 and lifts off as Q crosses 1.5 → 2.5 (0.34 → 0.49 → 0.56 → 0.70): a clean overdamped→revival handover at Q ≈ 1.5 on real hardware (the SE walk\'s critical-damping transition; real, measured). Whether that SE transition is itself a genuine defective EP is under a SEPARATE open review (the coherence-horizon √-EP cluster, inspect --root horizon); not asserted here. The F86a coherence-block \'real-axis EP\' this entry formerly cited was retracted 2026-06-21: the full (n,n+1)-coherence block has no eigenvalue coalescence on the real Q axis (genuine non-normality there, large but finite Petermann; the EP is not on the real axis). The revival decay envelope is gate-cost-limited (Trotterization, ~9 μs), not T2-limited (~200 μs); only the rate is gate cost, the floor and the onset are physics. This table is the hardware node of EpField (compute/RCPsiSquared.Diagnostics/Foundation/EpField.cs) and the overlay in simulations/ep_transition.py.',
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
                'model_test_in_situ': 'prep-conditioned split of the SAME circuits (pump = (s1+s0)/2, Γ = (s1−s0)/2, bound ⟺ s0 ≤ 0): pump/Γ = 0.972-0.979 (q149), 0.965-0.966 (q13), 0.994-0.996 (q9) — the bound HOLDS everywhere in-situ; margins read per-qubit thermal population (q13 1.7%, q149 1.1-1.4%, q9 0.2-0.3%)',
                'corrected_same_day': 'the first reading "q13 violates pump ≤ Γ at 4-6σ" compared the run against the 16-minutes-later arbiter and was an EPOCH ARTIFACT: T1 telegraphs on minute scale (q13 ~315 μs in-run vs 430 at the arbiter; q9 ~172 in-run vs ~75-100 at the arbiter; q149 stable across the two morning epochs). Two-level Lindblad holds within epochs; the protocol self-arbitrates via prep-splitting (no separate T1 experiment needed). Re-analysis: simulations/f120_prep_split_reanalysis.py',
                'telegraph_chase_1202Z': 'job d8l8f7r2d42s73cb3q7g (16 self-arbitrating blocks, ~75 s): within-job Γ flat at the shot-noise floor for all three qubits (no switch caught in 75 s), s₀ ≤ 0 in 47/48 block-cells (1 = expected false-positive rate); ACROSS the day every qubit moved, including the "stable" control (q149: 430 → ~285 μs; q13 in a third state ~200 μs after 430/315; q9 ~108, back near the arbiter state). Verdict: q13 was never special; the device T1 landscape breathes everywhere by 1.5-2× on minute-to-hour scales; the two-level model holds within any ~minute window; the hidden variable was time',
            },
            'hardware_data': 'data/ibm_moment_tower_june2026/ (main + arbiter JSONs + the 06:33Z calibration snapshot)',
            'experiment_doc': 'experiments/F120_MOMENT_TOWER_KINGSTON.md',
            'framework_primitive': 'moment_tower + predict_pump_slope + f113_bridge_asymmetry_from_slope (diagnostics/f120_moment_tower, called by the pipeline script at startup); MomentTowerPumpChannelClaim (C#); PROOF_MOMENT_TOWER_PUMP_CHANNEL',
            'description': 'F120 first hardware reading, the cleanest protocol we ever sent to a QPU (not one entangling gate). '
                           'The structural law confirmed: the first energy-moment rung stays silent while the second fires as exactly '
                           'twice the middle qubit\'s pump curve, the girth read from hardware is 2, and the per-qubit pump rates '
                           'reproduce across arms to 0.3-5.7%. The rate layer told a two-act story: the first reading ("q13 violates '
                           'pump ≤ Γ") was corrected the same day by the prep-conditioned re-analysis — the 8-basis-state preparation '
                           'contains the |0⟩- and |1⟩-branches, so pump AND Γ come from the same circuits, epoch-matched, and the bound '
                           'holds everywhere in-situ (worst 0.996) with 1-3% margins that READ the per-qubit thermal population. What '
                           'the cross-epoch comparison had actually detected is minute-scale T1 telegraphing on q13 (~315 ↔ 430 μs) and '
                           'q9 (~172 ↔ ~75-100 μs), once visible inside a single arm (q13, τ = 75→100). Two lessons banked: the protocol '
                           'is self-arbitrating (pump, Γ, γ↑ from one circuit set, the in-situ model test is s₀ ≤ 0), and the hardness '
                           'rung of a programmed Hamiltonian is now a quantity a chip measures about itself by decaying.',
        },
        'cpsi_quarter_crossing_torino_feb2026': {
            'date': '2026-02-09',
            'machine': 'ibm_torino',
            'job_id': 'tomography_ibm_torino_20260209_131521 (data-file timestamp; no IBM job_id recorded for the Torino calibration-era runs)',
            'observable': 'CΨ(t) = Tr(ρ²)·L₁/(d−1) for |+⟩ under free decoherence; crossing time t* of the CΨ = ¼ boundary',
            'predicted_value': 't*/T₂* = 0.936 (generalized, r = T₂*/T₁ = 0.456); 0.858 in the pure-dephasing limit (x³+x = ½)',
            'measured_value': 't* = 114.7 μs, t*/T₂* = 1.041 (11.3% above the generalized prediction); CΨ(0) = 0.885, C∞ = 0.740. Qubit 52: T1 = 221.2 μs, T2_echo = 298.2 μs, T2*(FID) = 110.7 μs',
            'hardware_data': 'data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json',
            'experiment_doc': 'experiments/IBM_QUANTUM_TOMOGRAPHY.md',
            'framework_primitive': 'F25 closed-form CΨ(t) + the CΨ = ¼ fold (K_fold dose)',
            'description': 'The first CΨ = ¼ crossing ever seen on hardware, found in IBM Torino single-qubit calibration tomography (Heron r2, q52, 2026-02-09), predating the systematic April-June 2026 campaign. A QUALITATIVE confirmation that the ¼ fold is real on a physical device, not a precision match: t*/T₂* = 1.041 sits 11% above the generalized prediction 0.936 because the crossing was extracted from a calibration run, not a purpose-built sweep. No IBM job_id was recorded for the Torino-era runs; the data-file timestamp is the locator.',
        },
        'absorption_theorem_ratio_torino': {
            'date': '2026-04-04',
            'machine': 'ibm_torino',
            'job_id': 'tomography_ibm_torino_20260209_131521 (data-file timestamp; analysis 2026-04-04 of the 2026-02-09 q52 run; no IBM job_id recorded)',
            'observable': 'Re(λ) / (−2γ⟨n_XY⟩), the excess coherence-decay ratio; effective ⟨n_XY⟩',
            'predicted_value': 'ratio = 1 (Absorption Theorem: Re(λ) = −2γ⟨n_XY⟩, with ⟨n_XY⟩ = 1 for single-qubit coherence)',
            'measured_value': 'ratio = 1.03 (3% deviation) on the T2* baseline: excess α = 0.006960 μs⁻¹ vs 2γ* = 0.006771 μs⁻¹, with a 2.8% slow tail at the resolution limit. (The 6.37 figure quoted elsewhere used the wrong T2_echo baseline; the dephasing-relevant T2* baseline gives 1.03.)',
            'hardware_data': 'data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json',
            'experiment_doc': 'experiments/IBM_ABSORPTION_THEOREM.md',
            'framework_primitive': 'Absorption Theorem Re(λ) = −2γ⟨n_XY⟩ (PROOF_ABSORPTION_THEOREM); simulations/ibm_absorption_theorem.py',
            'description': 'Retrospective Absorption-Theorem reading of the same 2026-02-09 Torino q52 tomography run (analysis 2026-04-04). The single-qubit coherence decays at the predicted Re(λ) = −2γ floor (⟨n_XY⟩ = 1), ratio 1.03. Predates the systematic registry campaign; shares the Feb-9 q52 data with cpsi_quarter_crossing_torino_feb2026.',
        },
        'cpsi_quarter_crossing_torino_q80_mar2026': {
            'date': '2026-03-18',
            'machine': 'ibm_torino',
            'job_id': 'palindrome_ibm_torino_20260318_191348 (data-file timestamp; no IBM job_id recorded)',
            'observable': 'CΨ = ¼ crossing time t* (8-point single-qubit tomography, q80 the "permanent crosser" of Run 3)',
            'predicted_value': 't* = 15.01 μs (from the same-day in-situ Ramsey T2* = 17.36 μs)',
            'measured_value': 't* = 15.29 μs, deviation 1.9% (0.28 μs). Qubit 80: T1 = 143.1-159 μs, same-day T2* = 17.36 μs (drifted 58% in 6 days from 11.0 μs, which is why the prediction uses the same-day Ramsey, not the stale calibration)',
            'hardware_data': 'data/ibm_run3_march2026/palindrome_ibm_torino_20260318_191348.json',
            'experiment_doc': 'experiments/IBM_RUN3_PALINDROME.md',
            'framework_primitive': 'F25 closed-form CΨ(t) + the CΨ = ¼ fold; same-day in-situ Ramsey T2*',
            'description': 'The tightest Torino-era CΨ = ¼ confirmation (1.9%): Run 3 on q80, 2026-03-18, with a same-day in-situ Ramsey T2* (17.36 μs). Predates the systematic April-June campaign but is precision-grade because the in-situ T2* removed the calibration-drift error (q80 T2* had drifted 58% in 6 days). The earlier two Torino rows (Feb-9 q52) are looser; this one pins the fold to 0.28 μs.',
        },
        'price_pair_locality_marrakesh_july2026': {
            'date': '2026-07-04',
            'machine': 'ibm_marrakesh',
            'job_id': 'd949n1tgc6cc73fer8sg + d949tgevtlqs73fu1v30 + d94abhvu62ks7395p4ig + d94ajjtgc6cc73fes7bg (the four-run pre-registered campaign)',
            'observable': 'coherence-pattern decay rates Γ(D) for all 7 subsets of a 3-qubit line (quadrature envelopes, readout-mitigated), the pair-sum price Γ(D)+Γ(D̄), the per-bond covariances c_ij (branch splitting), the pointwise shape-free law, and the conditional-Ramsey ζ_ij',
            'predicted_value': 'F89d fold at the ρ level: Γ(D)+Γ(D̄) = Σγ_j pattern-independent under LOCAL dephasing (all c_ij = 0); deviations decode by the pre-registered table (common-mode / single-bond / ZZ boundary rule / shape misfit)',
            'measured_value': 'dephasing covariances LOCAL (all c_ij = 0 within 2σ in clean sessions; one non-recurring 5.2σ outlier downgraded, unexplained); the one systematic deviation decoded as coherent NN ZZ, confirmed deterministically: ζ01 = −3.92±0.14 kHz, ζ12 = −3.64±0.10 kHz, ζ02 = +0.02±0.12 kHz (null next-nearest), consistent across three independent instruments (t² mask, 3:2:3 W fingerprint, conditional shift)',
            'hardware_data': 'data/ibm_price_pair_july2026/ (four hardware JSONs + null calibration + Aer parity + ZZ sim validation)',
            'experiment_doc': 'experiments/PRICE_PAIR_HARDWARE_PREDICTION.md',
            'framework_primitive': 'F89d cross-fold (dissipator half) + F1 palindrome center Σγ + F70 pattern-not-carrier + F82-F84 T1 accounting; run_price_pair.py (external pipeline)',
            'description': 'The locality premise under the F-formulas, measured: a four-run pre-registered campaign (runs 1-2 price protocol on two independent lines, run 3 the W discriminator refuting the intermediate anti-correlated-bath reading via the 3:2:3 boundary-rule fingerprint, run 4 conditional Ramsey confirming coherent NN ZZ ≈ −3.9/−3.6 kHz deterministically). Hardware noise model validated: local Z-dephasing (quasi-static + Markovian) + local T1 + coherent NN ZZ ≈ 4 kHz, no persistent correlated bath (one non-recurring 5.2σ covariance stays unexplained). Every wrong intermediate reading was falsified by its own pre-registered discriminator; the GHZ pattern-not-carrier immunity (XXX commutes with every Z_iZ_j) was derived before the first shot and carried the decoding.',
        },
        'f84_heating_leg_attribution_kingston_july2026': {
            'date': '2026-07-05',
            'machine': 'ibm_kingston + ibm_marrakesh',
            'job_id': 'd951mhkql68s73ca3u0g (Kingston [82,83,13]) + d953ti5gc6cc73ffomig (Marrakesh [93,94,95] repeat, same day)',
            'observable': 'two-legged idle T1 protocol (|111⟩ down-leg + |000⟩ up-leg, 10 delays 0-320 µs, readout-mitigated ⟨Z_l⟩(t)); per-qubit asymptote meeting test, γ↑ from the up leg, the F84 net-cooling attribution; run 2 adds the corrected-recipe V extraction',
            'predicted_value': 'pre-registered two-branch discrimination (committed pre-shot incl. the Kingston amendment): bath heating ⟹ both legs meet at one z∞ < 1 with γ↑ resolved ≥ 5σ; leg systematic ⟹ SPLIT with a flat up-leg (γ↑ < 0.0002 /µs, up-leg within ≈ 3% of +1). Instrument validated the same day on the IBM simulator (meeting test, planted-parameter recovery within 2σ, TLS double-flag, conservative error bars).',
            'measured_value': 'SPLIT on all three qubits (13.7σ / 5.1σ / 6.6σ): up-legs flat at z∞ = 0.983/0.994/0.995 (pointwise never below +0.978 over 320 µs) ⟹ p_th = (1−z∞_up)/2 = 0.83/0.31/0.23 % (the clean reading under a SPLIT; the joint-fit γ↑ = (6.8±2.1, 2.7±2.1, 2.1±2.0)·10⁻⁵ /µs is upper-bound-style, biased up by the shared-asymptote mis-specification), the bath is COLD; down-legs still rising at 320 µs with stretched (non-exponential) shape, so the sub-unity down-leg asymptotes (0.59-0.74) are extrapolation artifacts. Post-hoc (labeled as such): pinning z∞ at the up-leg level rejects the single exponential (χ²/dof 19.7/6.0/11.5) and accepts a two-rate mixture (1.45/0.32/1.91); q82\'s mixture mean rate reproduces the fresh calibration T1 (206 vs 208 µs; components 57% @ 136 µs + 43% @ 687 µs). RUN 2 (Marrakesh [93,94,95], same day, pre-registered 730ff94): no SPLIT, no χ² flag; γ↑ ≤ 2.2·10⁻⁵ /µs, populations ≤ 0.5% (up legs pointwise ≥ 0.994 over 320 µs; two of three separate up-leg fits asymptote-degenerate, the cold verdict rests on the joint fits + pointwise flatness); P3 resolved to the MOBILE-defect branch: q94 (χ²/dof 12.5 the day before; 12.1 matched-convention) returned clean (0.79) after its in-situ T1 rose 129 → 316 µs overnight (same-day calibration 228 µs); with nothing withheld, the first VALID measurement: V_σ± = 0.02491 ± 0.00015 /µs statistical (RMS velocity 0.00360 /µs), early-rate cross-check +6.9% (0.02663 ± 0.00074), a ≈7% estimator spread from residual mild convexity, both recorded; in-situ γ↓ vs same-day API calibration +7/−28/−22%.',
            'hardware_data': 'data/ibm_heating_leg_july2026/heating_leg_ibm_kingston_20260705_105648.json + heating_leg_ibm_marrakesh_20260705_132609.json (+ same-day simulator-validation JSONs and README in that directory)',
            'experiment_doc': 'experiments/F81_VIOLATION_HARDWARE_BRIDGE.md',
            'framework_primitive': 'F84 Pauli-channel cancellation + net-cooling closed form √(Σ(γ↓−γ↑)²)·2^(N−1); the identity-escape-velocity grounding (d⟨Z_l⟩/dt at ρ = I/2^N); run_heating_leg.py (external pipeline)',
            'description': 'The f81/F84 bridge\'s attribution question, answered by a pre-registered two-leg protocol in ≈ 1.6 QPU min on the f95/F113 qubits: the depressed T1-leg asymptotes seen on Marrakesh (price_pair Block B) and reproduced here are NOT thermal population (the up leg pins γ↑ ≈ 0) but a down-leg systematic: T1 fluctuates over the job, the ensemble average of exponentials is convex, and a free-asymptote single-exponential fit extrapolates a false z∞ < 1. Consequences recorded in the bridge doc: the one-leg asymptote recipe is biased on these devices; the corrected cold-bath net flux is a ≈ γ↓ with z∞ ≈ 1; the two-leg protocol (up leg for γ↑, fluctuation-robust rate estimate for γ↓) is the honest measurement of the F84 scalar. P1/P2 verdicts landed exactly by the pre-registered rules; the q82 χ² flag fired as designed. The 0-pending queue window and fresh 08:28Z calibration were Tom\'s catch.',
        },
        'concentrator_site_contrast_kingston_july2026': {
            'date': '2026-07-11',
            'machine': 'ibm_kingston',
            'job_id': 'd99a970tcv6s73dn2atg',
            'observable': 'paired-ratio slope difference slope(MP) − slope(E) of ln R_a(t) = ln[coh_a(t)/coh_0(t)] over the depth grid {1,2,3,4,6,8} Trotter steps at J·dt = 0.05; payload |+⟩ on position 2 of the 5-chain [109,108,107,106,105] (payload Q107), 2-basis Ramsey (X,Y), pooled coherence |mean(⟨X₂⟩+i⟨Y₂⟩)|; randomized Z-sink as M = 256 i.i.d. phase bindings per (arm, depth) at per-step retention e^{−0.05}, arm E on position 0 (distance 2), arm MP on the payload; three arms round-robin interleaved in ONE job (46 PUBs, 376,832 shots, billed 119 s ≈ 2.0 QPU min)',
            'predicted_value': 'pre-registered A-sign gate (v3.1, committed pre-shot after a five-round design loop): slope(MP) − slope(E) < 0 at one-sided 99.87% bootstrap AND outside the frozen null band [−0.01050, +0.01310]. From-below dressed central prediction −0.07325/step (7a exact density-matrix sim of the FLOWN circuits, 7a power margin 6.60; the operative margin is the 7b-reconciled 7.95 at projected SE 0.00921; MirrorWorld C# continuous cross-check −0.073618, the 0.5% gap = quantified Trotterization). Operative A-mag band [−0.09167, −0.05482] (7b-reconciled). THIS ENTRY REGISTERS A-SIGN ONLY, per the pre-registered verdict split.',
            'measured_value': 'slope(MP) − slope(E) = −0.05337/step; one-sided 99.87% bootstrap CI [−0.06505, −0.04107] (shot-only realized bootstrap, see description) entirely negative; far outside the null band; ≈ 5.8σ against the pre-registered projected SE 0.00921 ⟹ A-sign CONFIRMED. A-mag OFF-PREDICTION, marginal: 0.00145 above the operative band edge (0.16 projected SE, 0.33 shot-only SE) ⟹ per the committed rule "contrast confirmed, magnitude off-prediction"; the magnitude is NOT registered as confirmed. slope(E) = −0.00585 inside its null band [−0.01235, +0.00829] (L null-consistent, as pre-declared).',
            'hardware_data': 'external (AIEvolution.UI/experiments/ibm_quantum_tomography): results/reloaded_hardware_20260711_221840.json + reloaded_analyze_20260711_221936.json',
            'experiment_doc': 'experiments/IBM_CONCENTRATOR_RELOADED.md',
            'framework_primitive': 'site-resolved dephasing pricing (F1 palindrome center Σγ read per site, transport-dressed; the registered A-sign confirms the site-resolution, not the theorem magnitude, see description); 7a from-below gate simulations/concentrator_reloaded_7a.py; MirrorWorld Restless per-site γ + concentrator run mode (C# cross-check, N ≤ 9); run_concentrator_reloaded.py (external pipeline)',
            'description': 'First lifetime-rate contrast of the concentrator arc: the payload pays for a dephasing sink ON its site but not for the same dose two sites away, site-resolved on hardware at ≈ 5.8σ. SCOPE, honest: A-sign is the QUALITATIVE site-resolution, which any local-dephasing model also predicts; the theorem\'s quantitative discriminator was A-mag, which missed marginally (measured = dressed prediction uniformly compressed ≈ 0.73×; at this SNR "bare −2γ survives" and "gate-error compression of the dressed rate" are numerically degenerate, so NO bare-−2γ theorem claim is made). Instrument deviation disclosed: the persisted JSON pooled the M = 256 bindings (SamplerV2 get_counts()), so the realized bootstrap is shot-only and significance is reported against the pre-registered binding+shot projection; point estimates unaffected. The RECORD governs over the runner\'s stale-7a-band verdict printout. Arms E and MP share payload Q107, so its readout confusion and T1/T2 are common-mode and cancel in the verdict statistic. Registered 2026-07-12 (Tom\'s call per the RECORD\'s open Confirmations question): the contrast with the honest band story. Still owed elsewhere: the interior coherence lifetime under the full concentrator PROFILE (a different object).',
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
