"""RK4 propagation with hardware noise (T1/T2/Tphi/X-noise/Y-noise + transverse fields)."""
from __future__ import annotations

import numpy as np

from ..pauli import _build_bilinear, pauli_matrix, _site_op_kron
from ..receiver import Receiver


def propagate_with_hardware_noise(chain, rho_0, t, terms=None,
                                   T1_l=None, Tphi_l=None,
                                   T1pump_l=None,
                                   Xnoise_l=None, Ynoise_l=None,
                                   J_zz=None,
                                   h_x_l=None, h_y_l=None, h_z_l=None):
    """Propagate ρ_0 → ρ(t) under the full hardware noise Lindbladian.

    State-level bridge primitive: where `predict_residual_with_hardware_noise`
    gives the operator-level scalar ‖M‖², this method gives ρ(t) so the
    user can compute Pauli expectations, drift vs idealized, partial
    traces, etc. The two methods together cover both halves of the
    ON_THE_INSTRUMENT picture (operator and state level).

    Builds L = -i[H,·]
            + sum_l γ_Z_l (Z_l ρ Z_l - ρ)        from Tphi (or chain.gamma_0)
            + sum_l γ_T1_l D[σ⁻_l]
            + sum_l γ_T1pump_l D[σ⁺_l]
            + sum_l γ_X_l D[X_l]
            + sum_l γ_Y_l D[Y_l]
    with D[c]ρ = cρc† − ½{c†c, ρ}.

    Args:
        rho_0: 2^N × 2^N initial density matrix or Receiver instance.
        t: propagation time.
        terms: optional Pauli-pair Hamiltonian terms; if None uses chain.H.
        T1_l, Tphi_l, T1pump_l, Xnoise_l, Ynoise_l: per-site rate lists or
            None. Tphi_l (if given) overrides chain.gamma_0; otherwise
            chain.gamma_0 uniform applies.
        J_zz: optional ZZ-crosstalk strength. Adds Σ_(i,j) J_zz · Z_i Z_j
            over the chain's bond graph as a Hamiltonian correction (not
            a dissipator). Models the always-on ZZ-coupling on Heron r2.
        h_x_l, h_y_l: optional per-site transverse field strengths;
            add Σ_l h_l · X_l (resp. Y_l) as Hamiltonian terms. These
            BREAK Z⊗N-symmetry (single X/Y has odd n_XY).
        h_z_l: optional per-site longitudinal Z-detuning ("Mini-Magnetfeld");
            adds Σ_l h_l · Z_l. This PRESERVES Z⊗N (Z commutes with Z⊗N).

    Returns:
        ρ(t) as 2^N × 2^N complex array.
    """
    from scipy.linalg import expm
    from ..lindblad import lindbladian_general

    # Resolve initial ρ
    if isinstance(rho_0, Receiver):
        rho_0_mat = rho_0.rho
    else:
        rho_0_mat = np.asarray(rho_0, dtype=complex)
    d = chain.d
    if rho_0_mat.shape != (d, d):
        raise ValueError(f"rho_0 shape {rho_0_mat.shape} != ({d},{d})")

    # Resolve Hamiltonian
    if terms is not None:
        bilinear = [(t_pair[0], t_pair[1], chain.J) for t_pair in terms]
        H = _build_bilinear(chain.N, chain.bonds, bilinear)
    else:
        H = chain.H

    # Optional ZZ-crosstalk (Hamiltonian correction)
    if J_zz is not None and J_zz != 0:
        H_zz = _build_bilinear(chain.N, chain.bonds, [('Z', 'Z', float(J_zz))])
        H = H + H_zz

    # Optional single-site Hamiltonian fields:
    # - h_x_l: transverse X-field per site h_x · X_l (breaks Z⊗N!)
    # - h_y_l: transverse Y-field per site h_y · Y_l (breaks Z⊗N!)
    # - h_z_l: longitudinal Z-detuning δ · Z_l (Mini-Magnetfeld; preserves Z⊗N)
    Xm = pauli_matrix('X')
    Ym = pauli_matrix('Y')
    Zm = pauli_matrix('Z')
    if h_x_l is not None:
        for l, h in enumerate(h_x_l):
            if h != 0:
                H = H + float(h) * _site_op_kron(Xm, l, chain.N)
    if h_y_l is not None:
        for l, h in enumerate(h_y_l):
            if h != 0:
                H = H + float(h) * _site_op_kron(Ym, l, chain.N)
    if h_z_l is not None:
        for l, h in enumerate(h_z_l):
            if h != 0:
                H = H + float(h) * _site_op_kron(Zm, l, chain.N)

    # Build c_ops
    c_ops = []
    # Tphi (Z-dephasing): use Tphi_l if provided, else chain.gamma_0 uniform
    if Tphi_l is None:
        tphi_eff = [chain.gamma_0] * chain.N
    else:
        tphi_eff = [float(g) for g in Tphi_l]
        if len(tphi_eff) != chain.N:
            raise ValueError(f"Tphi_l length {len(tphi_eff)} != N {chain.N}")
    Z = pauli_matrix('Z')
    for l, g in enumerate(tphi_eff):
        if g > 0:
            c_ops.append(np.sqrt(g) * _site_op_kron(Z, l, chain.N))

    # T1 (σ⁻)
    if T1_l is not None:
        sigma_minus = np.array([[0, 1], [0, 0]], dtype=complex)
        for l, g in enumerate(T1_l):
            if g > 0:
                c_ops.append(np.sqrt(g) * _site_op_kron(sigma_minus, l, chain.N))

    # T1 pump (σ⁺)
    if T1pump_l is not None:
        sigma_plus = np.array([[0, 0], [1, 0]], dtype=complex)
        for l, g in enumerate(T1pump_l):
            if g > 0:
                c_ops.append(np.sqrt(g) * _site_op_kron(sigma_plus, l, chain.N))

    # X-noise
    if Xnoise_l is not None:
        X = pauli_matrix('X')
        for l, g in enumerate(Xnoise_l):
            if g > 0:
                c_ops.append(np.sqrt(g) * _site_op_kron(X, l, chain.N))

    # Y-noise
    if Ynoise_l is not None:
        Y = pauli_matrix('Y')
        for l, g in enumerate(Ynoise_l):
            if g > 0:
                c_ops.append(np.sqrt(g) * _site_op_kron(Y, l, chain.N))

    # Build full Lindbladian and propagate
    L = lindbladian_general(H, c_ops)
    rho_vec = rho_0_mat.flatten('F')
    rho_t_vec = expm(L * float(t)) @ rho_vec
    rho_t = rho_t_vec.reshape(d, d, order='F')
    return 0.5 * (rho_t + rho_t.conj().T)  # symmetrize numerical drift
