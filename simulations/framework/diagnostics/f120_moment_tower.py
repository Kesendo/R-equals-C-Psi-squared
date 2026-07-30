"""F120 moment-tower pump channel: the deg-1 girth ladder read linearly by a T1/pump slope."""
from __future__ import annotations

import math

import numpy as np

from ..pauli import site_op


def _n_from_H(H):
    d = H.shape[0]
    N = int(round(math.log2(d)))
    if 2 ** N != d:
        raise ValueError(f"H must be 2^N × 2^N; got shape {H.shape}")
    return N, d


def _resolve_per_site(values, N, name):
    if np.isscalar(values):
        return [float(values)] * N
    out = list(values)
    if len(out) != N:
        raise ValueError(f"{name} must have length {N}, got {len(out)}")
    return out


def moment_tower(H, N, j_max, tol=1e-9):
    """F120: the Z_l-weighted moment tower t_j(l) = Tr(Z_l H^j) and its first firing rung.

    Theorem F120 (the pump-slope law, verified in
    simulations/moment_tower_pump_channel.py): under
    L = −i[H,·] + Σ_l γ^deph_l D[Z_l] + Σ_l γ↓_l D[σ⁻_l] + Σ_l γ↑_l D[σ⁺_l]
    the slope of any measured polynomial A at the maximally mixed state is

        d/dt Tr(A ρ)|_{ρ = I/d} = (1/d) · Σ_l Δγ_l · Tr(A Z_l),   Δγ_l = γ↓_l − γ↑_l,

    so with A = H^j the slope reads the tower t_j(l) LINEARLY. Call the first j
    with t_j ≠ 0 somewhere the first nonzero rung k. Then p_{2k+1} carries the
    sum of squares (2k+1)·C(2k,k)·Σ_l t_k² > 0, so the pair is HARD and the
    hardness moment is BOUNDED: m* <= 2k+1. The certificate is ONE-SIDED:
    silence of the deg-1 tower through j_max is NOT softness (the 4-body pair
    IIXY+ZXZY has t_j = 0 everywhere, so it has no first nonzero rung at all, yet
    it is hard at m* = 11 with p₁₁ = 86507520·γ⁵).

    The bound is not an equality, which is why 'deg1_verdict' says <=. Equality
    needs k to be the F87 girth ℓ, and k is NOT ℓ: f87's ℓ
    (f87_windowed_monomial_converse.effective_ell) is 1 for a nonzero diagonal
    and otherwise the shortest odd cycle, so it is never even; k routinely is.
    At k > ℓ a higher-degree class fires BELOW 2k+1. H = Y₂ + Z₀Z₁Y₂ + Z₀Z₁Z₂ at
    N = 3 has t₁ = t₂ = 0 and t₃ = (0,0,16), so k = 3 and the bound reads 7, but
    the diagonal is nonzero so ℓ = 1 and p₅ = 7680·γ³ already fires: m* = 5. The
    dichotomy's `deg = 1` half, and with it 'hard at every γ > 0', need ℓ
    specifically and are not asserted here.

    Note also that this reads the TOWER. The rung at which a measured slope
    fires can be later, because the slope reads Σ_l Δγ_l t_j(l), which can
    cancel at a rung where t_j ≠ 0.

    Args:
        H: Hamiltonian as a 2^N × 2^N numpy array (Hermitian H gives real t_j).
        N: chain length (must match H's dimension).
        j_max: highest moment order to scan (j = 1..j_max).
        tol: firing threshold on |t_j(l)|.

    Returns:
        dict with keys:
            't':            {j: [t_j(l) per site l]} for j = 1..j_max (complex,
                            the f87_girth_dichotomy t_moment convention; exact
                            integers for integer H).
            'first_firing_rung': first j ≤ j_max with any |t_j(l)| > tol, or None.
                            NOT the F87 girth; see above.
            'deg1_verdict': 'hard, m* <= 2k+1' (with the number substituted)
                            if the tower fired, else the honest one-sided line
                            'deg-1 silent through j_max = ... (NOT a softness
                            certificate)'.
    """
    Nh, _ = _n_from_H(H)
    if Nh != N:
        raise ValueError(f"N={N} does not match H dimension 2^{Nh}")
    # Without this, j_max <= 0 returns the 'deg-1 silent' verdict from a scan that
    # never ran, shape-indistinguishable from a genuine silence certificate.
    if j_max < 1:
        raise ValueError(f"j_max must be at least 1 to scan a rung; got {j_max}")
    Z_ops = [site_op(N, l, 'Z') for l in range(N)]
    t = {}
    first_firing_rung = None
    Hp = np.eye(H.shape[0], dtype=complex)
    for j in range(1, j_max + 1):
        Hp = Hp @ H
        t[j] = [complex(np.trace(Z @ Hp)) for Z in Z_ops]
        if first_firing_rung is None and max(abs(v) for v in t[j]) > tol:
            first_firing_rung = j
    if first_firing_rung is not None:
        verdict = f"hard, m* <= {2 * first_firing_rung + 1}"
    else:
        verdict = f"deg-1 silent through j_max = {j_max} (NOT a softness certificate)"
    return {'t': t, 'first_firing_rung': first_firing_rung, 'deg1_verdict': verdict}


def predict_pump_slope(H, j, delta_gamma_l):
    """F120 closed form: the pump slope of A = H^j at the maximally mixed state.

        d/dt Tr(H^j ρ)|_{ρ = I/d} = (1/d) · Σ_l Δγ_l · t_j(l),   t_j(l) = Tr(Z_l H^j),

    with Δγ_l = γ↓_l − γ↑_l the per-site net cooling rate (σ⁻ minus σ⁺). The
    slope is exactly dephasing-blind, blind to the H generating the dynamics
    (only the measured polynomial enters), and zero at detailed balance
    (Δγ ≡ 0, the F84 vacuum-only statement).

    Args:
        H: Hermitian 2^N × 2^N Hamiltonian (the MEASURED polynomial's base).
        j: moment order (A = H^j).
        delta_gamma_l: per-site Δγ_l = γ↓_l − γ↑_l (list of length N or
            scalar for uniform).

    Returns:
        Predicted slope as float (real part; for Hermitian H the trace-pair
        is exactly real).
    """
    N, d = _n_from_H(H)
    dg = _resolve_per_site(delta_gamma_l, N, "delta_gamma_l")
    Hp = np.linalg.matrix_power(H, j)
    s = sum(dg[l] * np.trace(site_op(N, l, 'Z') @ Hp) for l in range(N)) / d
    return float(np.real(s))


def f113_bridge_asymmetry_from_slope(H, delta_gamma_l):
    """F113 bridge: the static Frobenius polarity asymmetry from the dynamic pump slope.

        ‖M₊‖² − ‖M₋‖² = −4^N · predict_pump_slope(H, 1, Δγ).

    Within F113's scope (H = Σ_l (ω_l/2)·Z_l with σ⁻ rates γ_T1,l and σ⁺
    rates γ_pump,l, so Δγ_l = γ_T1,l − γ_pump,l) this equals the F113 closed
    form (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l) and the actual asymmetry that
    polarity_coordinates / polarity_coordinates_from_hc returns (verified
    machine-exact at N = 2, 3 in simulations/moment_tower_pump_channel.py
    Block E). Outside that scope it is the slope reading −4^N·slope⟨H⟩,
    the dynamic side of the bridge.

    Args:
        H: Hermitian 2^N × 2^N Hamiltonian.
        delta_gamma_l: per-site Δγ_l = γ↓_l − γ↑_l (list of length N or scalar).

    Returns:
        Predicted asymmetry ‖M₊‖² − ‖M₋‖² as float.
    """
    N, _ = _n_from_H(H)
    return float(-(4 ** N) * predict_pump_slope(H, 1, delta_gamma_l))
