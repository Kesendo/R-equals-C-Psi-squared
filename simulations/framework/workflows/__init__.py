"""Workflows: composing primitives and diagnostics into reusable analysis flows."""
from .residual_norm import residual_norm_squared
from .zn_mirror import zn_mirror_diagnostic
from .gamma_probe import gamma_probe_setup, estimate_gamma_from_cpsi
from .propagate_with_hardware_noise import propagate_with_hardware_noise
from .predict_residual_with_hardware_noise import predict_residual_with_hardware_noise
from .cockpit_panel import cockpit_panel
from .predict_signature_table import predict_signature_table
from .diagnose_hardware import diagnose_hardware
from .lens import slow_modes, lens_pipeline
from .ptf import ptf_alpha_fit, ptf_painter_panel
from .handshake import verify_k_partnership
from .bridge_panel import bridge_panel
from .bridge_dynamics import bloch_trajectory, polarity_crossings, bridge_reflection_signature
