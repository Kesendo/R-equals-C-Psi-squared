"""R=CΨ² framework — lean cockpit package.

Layered as primitive modules + a Section-18 OOP cockpit on top.

Modules:
  pauli        — Pauli matrices, basis transforms, _build_bilinear
  symmetry     — Π conjugation, F71 chain-mirror, chiral panel, Y-parity panel
  lindblad     — Lindbladian, palindrome residual, ‖M‖² scaling
  observables  — Π-protected observables (Section 10)
  lebensader   — cockpit_panel (Section 11): skeleton + trace + cusp
  core         — ChainSystem, Receiver, Confirmations (Section 18 cockpit)

Quick start:
    import framework as fw

    # Hardware-confirmed predictions (no doc reading)
    fw.Confirmations.list_names()
    fw.Confirmations.lookup('palindrome_trichotomy')

    # Trichotomy classification
    chain = fw.ChainSystem(N=5)
    chain.classify_pauli_pair([('Y','Z'), ('Z','Y')])  # → 'soft'

    # Receiver-engineering forecast
    r = fw.Receiver(psi)
    r.signature()['prediction']

    # Full Lebensader analysis
    chain.cockpit_panel(r, terms=[('Y','Z'),('Z','Y')], gamma_t1=0.005)

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
)

from .symmetry import (
    pi_action, pi_squared_eigenvalue, build_pi_full,
    respects_bit_a_parity, respects_bit_b_parity, is_both_parity_even,
    chain_mirror_state, f71_symmetric_projector, f71_antisymmetric_projector,
    f71_eigenstate_class, receiver_engineering_signature, bond_mirror_basis,
    chiral_K_full, k_classify_pauli, k_classify_hamiltonian, chiral_panel,
    y_parity_panel,
)

from .lindblad import (
    lindbladian_general,
    lindbladian_z_dephasing,
    lindbladian_z_plus_t1,
    palindrome_residual,
    palindrome_residual_norm_squared_factor,
    palindrome_residual_norm_squared_factor_graph,
    palindrome_residual_norm_ratio_squared,
    dissipator_c1_c2_from_pauli,
    dissipator_d2_from_pauli,
    HARDWARE_DISSIPATORS,
    HARDWARE_DISSIPATOR_D1,
)

from .observables import pi_protected_observables

from .lebensader import cockpit_panel

# Cockpit OOP layer
from .core import ChainSystem, Receiver, Confirmations
