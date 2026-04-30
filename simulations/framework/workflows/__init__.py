"""Workflows: composing primitives and diagnostics into reusable analysis flows."""
from .residual_norm import residual_norm_squared
from .zn_mirror import zn_mirror_diagnostic
from .gamma_probe import gamma_probe_setup, estimate_gamma_from_cpsi
from .propagate_with_hardware_noise import propagate_with_hardware_noise
from .predict_residual_with_hardware_noise import predict_residual_with_hardware_noise
from .cockpit_panel import cockpit_panel
from .predict_signature_table import predict_signature_table
from .diagnose_hardware import diagnose_hardware
