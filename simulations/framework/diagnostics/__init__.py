"""Diagnostics: F-theorem readings as standalone functions taking ChainSystem."""
from .f77_trichotomy import classify_pauli_pair
from .f49_frobenius_scaling import (
    predict_residual_norm_squared,
    predict_residual_norm_squared_from_terms,
)
from .f80_bloch_signwalk import predict_M_spectrum_pi2_odd
from .f81_pi_decomposition import pi_decompose_M, recover_H_odd_from_M_anti
from .f82_t1_dissipator import (
    predict_T1_dissipator_violation,
    estimate_T1_from_violation,
)
from .f83_anti_fraction import (
    predict_pi_decomposition,
    predict_pi_decomposition_anti_fraction,
)
from .f84_amplitude_damping import (
    predict_amplitude_damping_violation,
    estimate_net_cooling_from_violation,
)
from .ptf import pt_matrix_elements, pt_eigvec_shift
from .polarity import polarity_diagnostic
from .d_zero import stationary_modes, d_zero_decomposition
