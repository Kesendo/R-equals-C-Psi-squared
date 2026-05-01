"""R=CΨ² framework — lean cockpit package.

Layered as primitive modules + a Section-18 OOP cockpit on top.

Modules:
  pauli         — Pauli matrices, basis transforms, _build_bilinear
  symmetry      — Π conjugation, F71 chain-mirror, chiral panel, Y-parity panel
  lindblad      — Lindbladian, palindrome residual, ‖M‖² scaling
  observables   — Π-protected observables (Section 10)
  lebensader    — cockpit_panel (Section 11): skeleton + trace + cusp
  chain_system  — ChainSystem (Section 18 cockpit, the workhorse class)
  receiver      — Receiver (state-bearing F71-aware wrapper)
  confirmations — Confirmations (hardware-confirmed predictions registry)

Quick start:
    import framework as fw

    # Hardware-confirmed predictions (no doc reading)
    fw.Confirmations.list_names()
    fw.Confirmations.lookup('palindrome_trichotomy')

    # Trichotomy classification (diagnostics package)
    chain = fw.ChainSystem(N=5)
    fw.classify_pauli_pair(chain, [('Y','Z'), ('Z','Y')])  # → 'soft'

    # Receiver-engineering forecast
    r = fw.Receiver(psi)
    r.signature()['prediction']

    # Full Lebensader analysis
    fw.cockpit_panel(chain, r, terms=[('Y','Z'),('Z','Y')], gamma_t1=0.005)

The legacy primitives (Sections 1-15 of the original framework.py that aren't
part of the cockpit path) live in framework_archive.py — kept in repo for
historical reference, not imported here.
"""

# Primitives
from .pauli import (
    PAULI_LABELS, LABEL_TO_INDEX,
    ur_pauli, pauli_matrix,
    bit_a, bit_b, total_bit_a, total_bit_b_parity,
    pauli_string, site_op,
    _resolve, _k_to_indices, _indices_to_k, _pauli_label,
    _vec_to_pauli_basis_transform, pauli_basis_vector,
    _build_bilinear, _site_op_kron,
    bonding_mode_state, bonding_mode_pair_state,
    polarity_state,
)

from .symmetry import (
    pi_action, pi_squared_eigenvalue, build_pi_full,
    respects_bit_a_parity, respects_bit_b_parity, is_both_parity_even,
    chain_mirror_state, f71_symmetric_projector, f71_antisymmetric_projector,
    f71_eigenstate_class, receiver_engineering_signature, bond_mirror_basis,
    chiral_K_full, k_classify_pauli, k_classify_hamiltonian, chiral_panel,
    y_parity_panel, zn_mirror_state,
    klein_index, KLEIN_LABELS, k_partner,
)

from .lindblad import (
    lindbladian_general,
    lindbladian_z_dephasing,
    lindbladian_pauli_dephasing,
    lindbladian_z_plus_t1,
    bond_perturbation,
    palindrome_residual,
    palindrome_residual_norm_squared_factor,
    palindrome_residual_norm_squared_factor_graph,
    palindrome_residual_norm_ratio_squared,
    dissipator_c1_c2_from_pauli,
    dissipator_d2_from_pauli,
    HARDWARE_DISSIPATORS,
    HARDWARE_DISSIPATOR_D1,
    cpsi_bell_plus,
    CPSI_CUSP_K_PER_CHANNEL,
)

from .observables import pi_protected_observables

# Cockpit OOP layer
from .chain_system import ChainSystem
from .receiver import Receiver
from .confirmations import Confirmations
from .pauli_hamiltonian import PauliTerm, PauliHamiltonian

# Diagnostics: F-theorem readings as free functions
from .diagnostics import (
    classify_pauli_pair,
    predict_residual_norm_squared,
    predict_residual_norm_squared_from_terms,
    predict_M_spectrum_pi2_odd,
    pi_decompose_M,
    recover_H_odd_from_M_anti,
    predict_T1_dissipator_violation,
    estimate_T1_from_violation,
    predict_pi_decomposition,
    predict_pi_decomposition_anti_fraction,
    predict_amplitude_damping_violation,
    estimate_net_cooling_from_violation,
    pt_matrix_elements,
    pt_eigvec_shift,
    polarity_diagnostic,
    stationary_modes,
    d_zero_decomposition,
    sector_populations,
)

# Workflows: composing primitives and diagnostics into analysis flows
from .workflows import (
    residual_norm_squared,
    zn_mirror_diagnostic,
    gamma_probe_setup,
    estimate_gamma_from_cpsi,
    propagate_with_hardware_noise,
    predict_residual_with_hardware_noise,
    cockpit_panel,
    predict_signature_table,
    diagnose_hardware,
    slow_modes,
    lens_pipeline,
    ptf_alpha_fit,
    ptf_painter_panel,
    verify_k_partnership,
    bridge_panel,
    bloch_trajectory,
    polarity_crossings,
    bridge_reflection_signature,
)
