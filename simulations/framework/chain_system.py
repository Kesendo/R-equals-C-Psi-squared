"""ChainSystem entity: encapsulates a quantum chain at fixed (N, γ, J, topology, H_type)."""
from __future__ import annotations

import warnings

import numpy as np

from .lindblad import (
    HARDWARE_DISSIPATORS,
    HARDWARE_DISSIPATOR_D1,
)
from .pauli import _build_bilinear, pauli_matrix, _site_op_kron
from .lebensader import cockpit_panel
from .receiver import Receiver


class ChainSystem:
    """Encapsulates a quantum chain at fixed (N, γ, J, topology, H_type).

    Caches Hamiltonian and Liouvillian on first access. Methods compose
    framework primitives into single-call answers.

    Args:
        N: number of qubits.
        gamma_0: uniform Z-dephasing rate per site (default 0.05).
        J: bond coupling (uniform; default 1.0).
        topology: 'chain' (default), 'ring', 'star', 'complete'.
        H_type: 'heisenberg' (XX+YY+ZZ) or 'xy' (XX+YY scaled by J/2).

    F-theorem readings live as free functions in `framework.diagnostics`
    (e.g. `fw.classify_pauli_pair(chain, terms)`,
    `fw.predict_residual_norm_squared_from_terms(chain, terms, gamma_t1=...)`,
    `fw.pi_decompose_M(chain, terms, ...)`). Workflow methods that
    propagate ρ(t) or hold cusp/probe state remain on this class.
    """

    _FROZEN = frozenset({'N', 'd', 'd2', 'gamma_0', 'J', 'topology', 'H_type',
                          'bonds', 'B', 'degrees', 'D2'})

    def __init__(self, N, gamma_0=0.05, J=1.0, topology='chain', H_type='heisenberg'):
        if N < 2:
            raise ValueError(f"N must be >= 2; got {N}")
        if H_type not in ('heisenberg', 'xy'):
            raise ValueError(f"H_type must be 'heisenberg' or 'xy'; got {H_type!r}")
        if N == 2:
            warnings.warn(
                "ChainSystem(N=2) is mathematically valid but structurally "
                "degenerate: F71 bond_mirror_basis has only 1 sym mode and 0 "
                "asym modes (every F71-eigenstate is capacity-suboptimal); "
                "only 1 bond exists; the drop=28 hardware anchor is "
                "unreproducible (4^2=16 Pauli strings total). Fundamental "
                "vocabulary (classify, Π, palindrome residual, Frobenius "
                "scaling) holds; structural vocabulary (F71-balance, multi-"
                "bond cockpit signatures) needs N>=3.",
                UserWarning, stacklevel=2,
            )
        # bypass __setattr__ guard during init
        object.__setattr__(self, '_initialized', False)
        self.N = N
        self.d = 2 ** N
        self.d2 = self.d * self.d
        self.gamma_0 = float(gamma_0)
        self.J = float(J)
        self.topology = topology
        self.H_type = H_type
        self._build_topology()
        self._H_cache = None
        self._L_cache = None
        object.__setattr__(self, '_initialized', True)

    def __setattr__(self, name, value):
        # Freeze structural attributes after construction. Caches (H, L) are
        # immutable in spirit because the Hamiltonian / Liouvillian depend on
        # all the frozen attrs; allowing J or gamma_0 mutation would silently
        # decouple cache from declared chain state. Make a new ChainSystem
        # instead.
        if getattr(self, '_initialized', False) and name in self._FROZEN:
            raise AttributeError(
                f"ChainSystem.{name} is immutable after construction. "
                f"Make a new ChainSystem(...) with the desired {name}."
            )
        object.__setattr__(self, name, value)

    def _build_topology(self):
        N = self.N
        if self.topology == 'chain':
            self.bonds = [(i, i + 1) for i in range(N - 1)]
        elif self.topology == 'ring':
            self.bonds = [(i, (i + 1) % N) for i in range(N)]
        elif self.topology == 'star':
            self.bonds = [(0, i) for i in range(1, N)]
        elif self.topology == 'complete':
            self.bonds = [(i, j) for i in range(N) for j in range(i + 1, N)]
        else:
            raise ValueError(f"unknown topology: {self.topology!r}")
        self.B = len(self.bonds)
        deg = [0] * N
        for i, j in self.bonds:
            deg[i] += 1
            deg[j] += 1
        self.degrees = deg
        self.D2 = sum(d * d for d in deg)

    @property
    def H(self):
        if self._H_cache is None:
            self._H_cache = self._build_H()
        return self._H_cache

    @property
    def L(self):
        if self._L_cache is None:
            self._L_cache = self._build_L()
        return self._L_cache

    def _build_H(self):
        d = self.d
        H = np.zeros((d, d), dtype=complex)
        Xm, Ym, Zm = pauli_matrix('X'), pauli_matrix('Y'), pauli_matrix('Z')
        for (i, j) in self.bonds:
            xi, xj = _site_op_kron(Xm, i, self.N), _site_op_kron(Xm, j, self.N)
            yi, yj = _site_op_kron(Ym, i, self.N), _site_op_kron(Ym, j, self.N)
            if self.H_type == 'heisenberg':
                zi, zj = _site_op_kron(Zm, i, self.N), _site_op_kron(Zm, j, self.N)
                H = H + self.J * (xi @ xj + yi @ yj + zi @ zj)
            elif self.H_type == 'xy':
                H = H + (self.J / 2.0) * (xi @ xj + yi @ yj)
        return H

    def _build_L(self):
        H = self.H
        d, d2 = self.d, self.d2
        Id = np.eye(d, dtype=complex)
        L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
        Zm = pauli_matrix('Z')
        for k in range(self.N):
            Zk = _site_op_kron(Zm, k, self.N)
            L = L + self.gamma_0 * (np.kron(Zk, Zk.conj()) - np.eye(d2))
        return L

    def residual_norm_squared(self, terms, J_scale=None, gamma_t1=None):
        """Numerical ‖M‖² for a Pauli-pair Hamiltonian on this chain.

        Builds H, the Lindbladian (Z-dephasing or Z+T1 if `gamma_t1` is set),
        and the palindrome residual M; returns ‖M‖²_F. For 'truly' classes
        with gamma_t1=0, ≈ 0 to floating-point precision.

        Replaces the recurring boilerplate
        (_build_bilinear → lindbladian_* → palindrome_residual → norm).

        Args:
            terms: list of (a, b) letter tuples.
            J_scale: optional override for chain.J (default uses self.J).
            gamma_t1: optional T1 amplitude-damping rates. Scalar (uniform),
                      list of length N, or None / 0 (pure Z-dephasing).
        """
        from .lindblad import (
            lindbladian_z_dephasing,
            lindbladian_z_plus_t1,
            palindrome_residual,
        )
        J = self.J if J_scale is None else float(J_scale)
        bilinear = [(a, b, J) for (a, b) in terms]
        H = _build_bilinear(self.N, self.bonds, bilinear)
        gamma_z_l = [self.gamma_0] * self.N
        if gamma_t1 is None:
            L = lindbladian_z_dephasing(H, gamma_z_l)
        else:
            if np.isscalar(gamma_t1):
                gamma_t1_l = [float(gamma_t1)] * self.N
            else:
                gamma_t1_l = [float(g) for g in gamma_t1]
                if len(gamma_t1_l) != self.N:
                    raise ValueError(
                        f"gamma_t1 list length {len(gamma_t1_l)} != N {self.N}"
                    )
            L = lindbladian_z_plus_t1(H, gamma_z_l, gamma_t1_l)
        Sigma_gamma = sum(gamma_z_l)
        M = palindrome_residual(L, Sigma_gamma, self.N)
        return float(np.linalg.norm(M) ** 2)

    def zn_mirror_diagnostic(self, rho_a, rho_b, tol=1e-6):
        """Z⊗N-Mirror Symmetrie-Test zwischen zwei Dichtematrizen.

        Wenn ρ_a und ρ_b Z⊗N-Partner sind (ρ_b = (Z⊗N)·ρ_a·(Z⊗N)), gilt für
        jede Pauli-String-Erwartung:

            ⟨P⟩_b = (−1)^n_XY(P) · ⟨P⟩_a

        wobei n_XY(P) = Anzahl X- oder Y-Buchstaben in P. Diese Identität
        scheitert in Anwesenheit von transverse fields h_l X_l oder h_l Y_l;
        sie gilt für XXZ, Z-Dephasing, T1 (σ⁻σ⁺ pairs), und non-uniform
        Z-Detuning (Mini-Magnetfeld δ_l Z_l).

        Liefert max-Abweichung über alle 4^N Pauli-Strings → 'preserved' wenn
        unter tol, sonst 'broken'.

        Args:
            rho_a: 2^N × 2^N Dichtematrix
            rho_b: 2^N × 2^N Dichtematrix (vermuteter Z⊗N-Partner von ρ_a)
            tol: Toleranz für 'preserved' Verdict.

        Returns:
            dict mit 'max_violation' (max Pauli-String-Abweichung),
            'verdict' ('preserved' | 'broken'), 'worst_string' (Pauli-Label
            wo die Abweichung am größten), 'worst_a', 'worst_b' (die zwei
            Erwartungen die nicht stimmen).
        """
        from .pauli import _PAULI_MATRICES, _k_to_indices, _pauli_label, bit_a
        N = self.N
        d = 2 ** N
        rho_a = np.asarray(rho_a, dtype=complex)
        rho_b = np.asarray(rho_b, dtype=complex)
        if rho_a.shape != (d, d) or rho_b.shape != (d, d):
            raise ValueError(
                f"rho_a, rho_b must be {d}×{d} for N={N}; got "
                f"{rho_a.shape}, {rho_b.shape}"
            )
        max_violation = 0.0
        worst_k = 0
        worst_a = 0.0
        worst_b = 0.0
        for k in range(4 ** N):
            indices = _k_to_indices(k, N)
            # Build Pauli string σ_indices
            P = _PAULI_MATRICES[indices[0]]
            for idx in indices[1:]:
                P = np.kron(P, _PAULI_MATRICES[idx])
            exp_a = float(np.real(np.trace(rho_a @ P)))
            exp_b = float(np.real(np.trace(rho_b @ P)))
            n_xy = sum(bit_a(idx) for idx in indices)  # X, Y have bit_a=1
            sign = (-1) ** n_xy
            violation = abs(exp_b - sign * exp_a)
            if violation > max_violation:
                max_violation = violation
                worst_k = k
                worst_a = exp_a
                worst_b = exp_b
        worst_label = _pauli_label(worst_k, N)
        return {
            'max_violation': max_violation,
            'verdict': 'preserved' if max_violation < tol else 'broken',
            'worst_string': worst_label,
            'worst_a': worst_a,
            'worst_b': worst_b,
            'tol': tol,
        }

    def gamma_probe_setup(self, gamma_assumed=None, target_precision=0.01,
                           channel='Z'):
        """Optimal γ-Sensing-Parameter via Cusp-nahe CΨ-Probe (Bell+).

        Generalizes F25 → F26 (multi-axis Pauli channels). The cusp at CΨ=1/4
        is a good probe region but not the Fisher-Information optimum; the
        optimum lies post-cusp.

        Channels:
            'Z'             — pure Z-dephasing only (γ_z = γ; F25 special case)
            'X'             — pure X-noise (γ_x = γ)
            'Y'             — pure Y-noise (γ_y = γ; symmetric with X for Bell+)
            'depolarizing'  — γ/3 on each axis

        K = γ·t is γ-invariant; channel determines K_optimal numerically.
        Reference K_cusp values per channel (F27): Z=0.0374, X=0.0867,
        Y=0.0867, depolarizing=0.0440.

        Args:
            gamma_assumed: prior γ-estimate. If None, uses self.gamma_0.
            target_precision: relative precision Δγ/γ desired.
            channel: noise channel — one of 'Z', 'X', 'Y', 'depolarizing'.

        Returns:
            dict with t_optimal, cpsi_target, K_optimal, fisher_per_shot,
            shots_needed, plus cusp parameters and the channel itself.
        """
        from scipy.optimize import minimize_scalar, brentq
        from .lindblad import cpsi_bell_plus

        if gamma_assumed is None:
            gamma_assumed = self.gamma_0
        gamma = float(gamma_assumed)

        # Map channel → axis assignment of γ
        def gammas_for_channel(g):
            if channel == 'Z':
                return (0.0, 0.0, g)
            if channel == 'X':
                return (g, 0.0, 0.0)
            if channel == 'Y':
                return (0.0, g, 0.0)
            if channel == 'depolarizing':
                return (g/3, g/3, g/3)
            raise ValueError(f"channel must be 'Z'/'X'/'Y'/'depolarizing'; got {channel!r}")

        gx, gy, gz = gammas_for_channel(gamma)

        def cpsi(t):
            return cpsi_bell_plus(gx, gy, gz, t)

        def dcpsi_dgamma(t):
            # Numerical derivative — analytic form is messy for general channel.
            eps = max(1e-7, 1e-6 * gamma)
            gx_p, gy_p, gz_p = gammas_for_channel(gamma + eps)
            gx_m, gy_m, gz_m = gammas_for_channel(gamma - eps)
            return (cpsi_bell_plus(gx_p, gy_p, gz_p, t)
                    - cpsi_bell_plus(gx_m, gy_m, gz_m, t)) / (2 * eps)

        def neg_fisher(t):
            c = cpsi(t)
            dc = dcpsi_dgamma(t)
            return -(dc * dc / max(1 - c * c, 1e-12))

        res = minimize_scalar(neg_fisher,
                              bounds=(0.001 / gamma, 5.0 / gamma),
                              method='bounded')
        t_opt = float(res.x)
        cpsi_opt = float(cpsi(t_opt))
        fisher_opt = -float(res.fun)

        # Cusp t (where CΨ = 1/4)
        try:
            t_cusp = brentq(lambda t: cpsi(t) - 0.25, 1e-6, 1000.0 / gamma)
            cpsi_cusp = 0.25
        except ValueError:
            t_cusp = float('nan')
            cpsi_cusp = float('nan')

        delta_gamma = float(target_precision) * gamma
        shots = int(np.ceil(1.0 / (fisher_opt * delta_gamma * delta_gamma)))

        return {
            'channel':          channel,
            't_optimal':        t_opt,
            'cpsi_target':      cpsi_opt,
            'K_optimal':        gamma * t_opt,
            'fisher_per_shot':  fisher_opt,
            'shots_needed':     shots,
            'target_precision': float(target_precision),
            't_cusp':           float(t_cusp),
            'K_cusp':           gamma * float(t_cusp) if not np.isnan(t_cusp) else float('nan'),
            'cpsi_cusp':        cpsi_cusp,
            'gamma_assumed':    gamma,
        }

    def estimate_gamma_from_cpsi(self, cpsi_measured, t, channel='Z'):
        """Invert F26 to extract γ from CΨ at time t (Bell+, multi-axis channels).

        For 'Z' channel uses F25 cubic root (closed form). For 'X', 'Y',
        'depolarizing' uses numerical inversion of F26.

        Bell+ initial CΨ(t=0) = 1/3, so cpsi_measured must satisfy
        0 < cpsi_measured < 1/3.

        Args:
            cpsi_measured: measured CΨ value (e.g. from hardware tomography).
            t: probe time at which CΨ was measured.
            channel: 'Z' (default), 'X', 'Y', or 'depolarizing'.

        Returns:
            Estimated γ as float.
        """
        from scipy.optimize import brentq
        from .lindblad import cpsi_bell_plus
        cpsi = float(cpsi_measured)
        if cpsi >= 1.0 / 3.0:
            raise ValueError(
                f"cpsi_measured = {cpsi} ≥ 1/3 = CΨ(t=0); cannot invert "
                f"(Bell+ never starts above this)."
            )
        if cpsi <= 0:
            raise ValueError(
                f"cpsi_measured = {cpsi} ≤ 0; physically invalid for Bell+."
            )
        if t <= 0:
            raise ValueError(f"t must be > 0; got {t}")

        if channel == 'Z':
            f = brentq(lambda fv: fv * (1 + fv * fv) / 6.0 - cpsi, 1e-12, 1.0 - 1e-12)
            return -float(np.log(f)) / (4.0 * float(t))

        # General channel: numerical inversion of F26
        def gammas_for_channel(g):
            if channel == 'X':            return (g, 0.0, 0.0)
            if channel == 'Y':            return (0.0, g, 0.0)
            if channel == 'depolarizing': return (g/3, g/3, g/3)
            raise ValueError(f"channel must be 'Z'/'X'/'Y'/'depolarizing'; got {channel!r}")

        def f(gamma):
            gx, gy, gz = gammas_for_channel(gamma)
            return cpsi_bell_plus(gx, gy, gz, t) - cpsi

        # CΨ monotonically decreases in γ ⇒ use bisection
        gamma_est = brentq(f, 1e-12, 1000.0 / float(t))
        return float(gamma_est)

    def propagate_with_hardware_noise(self, rho_0, t, terms=None,
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
        from .lindblad import lindbladian_general

        # Resolve initial ρ
        if isinstance(rho_0, Receiver):
            rho_0_mat = rho_0.rho
        else:
            rho_0_mat = np.asarray(rho_0, dtype=complex)
        d = self.d
        if rho_0_mat.shape != (d, d):
            raise ValueError(f"rho_0 shape {rho_0_mat.shape} != ({d},{d})")

        # Resolve Hamiltonian
        if terms is not None:
            bilinear = [(t_pair[0], t_pair[1], self.J) for t_pair in terms]
            H = _build_bilinear(self.N, self.bonds, bilinear)
        else:
            H = self.H

        # Optional ZZ-crosstalk (Hamiltonian correction)
        if J_zz is not None and J_zz != 0:
            H_zz = _build_bilinear(self.N, self.bonds, [('Z', 'Z', float(J_zz))])
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
                    H = H + float(h) * _site_op_kron(Xm, l, self.N)
        if h_y_l is not None:
            for l, h in enumerate(h_y_l):
                if h != 0:
                    H = H + float(h) * _site_op_kron(Ym, l, self.N)
        if h_z_l is not None:
            for l, h in enumerate(h_z_l):
                if h != 0:
                    H = H + float(h) * _site_op_kron(Zm, l, self.N)

        # Build c_ops
        c_ops = []
        # Tphi (Z-dephasing): use Tphi_l if provided, else chain.gamma_0 uniform
        if Tphi_l is None:
            tphi_eff = [self.gamma_0] * self.N
        else:
            tphi_eff = [float(g) for g in Tphi_l]
            if len(tphi_eff) != self.N:
                raise ValueError(f"Tphi_l length {len(tphi_eff)} != N {self.N}")
        Z = pauli_matrix('Z')
        for l, g in enumerate(tphi_eff):
            if g > 0:
                c_ops.append(np.sqrt(g) * _site_op_kron(Z, l, self.N))

        # T1 (σ⁻)
        if T1_l is not None:
            sigma_minus = np.array([[0, 1], [0, 0]], dtype=complex)
            for l, g in enumerate(T1_l):
                if g > 0:
                    c_ops.append(np.sqrt(g) * _site_op_kron(sigma_minus, l, self.N))

        # T1 pump (σ⁺)
        if T1pump_l is not None:
            sigma_plus = np.array([[0, 0], [1, 0]], dtype=complex)
            for l, g in enumerate(T1pump_l):
                if g > 0:
                    c_ops.append(np.sqrt(g) * _site_op_kron(sigma_plus, l, self.N))

        # X-noise
        if Xnoise_l is not None:
            X = pauli_matrix('X')
            for l, g in enumerate(Xnoise_l):
                if g > 0:
                    c_ops.append(np.sqrt(g) * _site_op_kron(X, l, self.N))

        # Y-noise
        if Ynoise_l is not None:
            Y = pauli_matrix('Y')
            for l, g in enumerate(Ynoise_l):
                if g > 0:
                    c_ops.append(np.sqrt(g) * _site_op_kron(Y, l, self.N))

        # Build full Lindbladian and propagate
        L = lindbladian_general(H, c_ops)
        rho_vec = rho_0_mat.flatten('F')
        rho_t_vec = expm(L * float(t)) @ rho_vec
        rho_t = rho_t_vec.reshape(d, d, order='F')
        return 0.5 * (rho_t + rho_t.conj().T)  # symmetrize numerical drift

    def predict_residual_with_hardware_noise(self, terms=None,
                                              T1_l=None, Tphi_l=None,
                                              T1pump_l=None,
                                              Xnoise_l=None, Ynoise_l=None):
        """Predict ‖M‖² for an IBM-like hardware noise model in closed form.

        Aggregates the per-class Frobenius forms and the cross-class d1, d2
        terms from HARDWARE_DISSIPATORS / HARDWARE_DISSIPATOR_D1. Returns a
        dict with the per-class contributions plus the total.

        All rate arguments are per-site lists (length N) of floats. Pass None
        to omit a class. Pass terms=None for purely-dissipative ‖M‖².

        σ_offset = 0 (raw Lindbladian L; do NOT subtract Π-trivial Z-dephasing
        offset, since that's already absorbed into c1=0 for Tphi/Ynoise).

        Args:
            terms: optional Pauli-pair Hamiltonian terms; if None uses no H.
            T1_l, Tphi_l, T1pump_l, Xnoise_l, Ynoise_l: per-site rate lists or
                None (= zero). Lists must be length N.

        Returns:
            dict with keys 'hamiltonian', 'per_class' (sub-dict), 'cross'
            (sub-dict), and 'total' (the sum).
        """
        from .diagnostics.f49_frobenius_scaling import (
            predict_residual_norm_squared_from_terms,
        )
        rate_args = {
            'T1':     T1_l,
            'T1pump': T1pump_l,
            'Tphi':   Tphi_l,
            'Xnoise': Xnoise_l,
            'Ynoise': Ynoise_l,
        }
        rates = {}
        for name, gl in rate_args.items():
            if gl is None or all(g == 0 for g in gl):
                continue
            if len(gl) != self.N:
                raise ValueError(
                    f"{name}_l length {len(gl)} != N {self.N}"
                )
            rates[name] = list(gl)

        factor = 4 ** (self.N - 1)

        # Hamiltonian palindrome part — reuse predict_residual_norm_squared_from_terms
        if terms:
            h_part = predict_residual_norm_squared_from_terms(self, terms)
        else:
            h_part = 0.0

        # Per-class contributions
        per_class = {}
        for name, gl in rates.items():
            spec = HARDWARE_DISSIPATORS[name]
            sum_g_sq = sum(g * g for g in gl)
            sum_g = sum(gl)
            contrib = factor * (spec['c1'] * sum_g_sq + spec['c2'] * sum_g * sum_g)
            per_class[name] = contrib

        # Cross-class contributions
        cross = {}
        active = list(rates.keys())
        for i, k1 in enumerate(active):
            for j, k2 in enumerate(active):
                if i >= j:
                    continue
                key = (k1, k2) if (k1, k2) in HARDWARE_DISSIPATOR_D1 else (k2, k1)
                d1 = HARDWARE_DISSIPATOR_D1[key]
                p1 = HARDWARE_DISSIPATORS[k1]['pauli']
                p2 = HARDWARE_DISSIPATORS[k2]['pauli']
                d2 = 32.0 * (sum(abs(x)**2 for x in p1)) * (sum(abs(x)**2 for x in p2))
                gl1, gl2 = rates[k1], rates[k2]
                sum_g1_g2 = sum(g1 * g2 for g1, g2 in zip(gl1, gl2))
                prod_sums = sum(gl1) * sum(gl2)
                contrib = factor * (d1 * sum_g1_g2 + d2 * prod_sums)
                cross[(k1, k2)] = contrib

        total = h_part + sum(per_class.values()) + sum(cross.values())
        return {
            'hamiltonian': h_part,
            'per_class':   per_class,
            'cross':       cross,
            'total':       total,
        }

    def cockpit_panel(self, receiver, terms=None, gamma_t1=None,
                      t_max=10.0, dt=0.005,
                      threshold=1e-9, cluster_tol=1e-8):
        """Full Lebensader cockpit panel: skeleton + trace + cusp + chiral + Y-parity.

        Args:
            receiver: Receiver instance providing ρ_0.
            terms: optional Pauli-pair Hamiltonian terms; if None uses chain's default H.
            gamma_t1: scalar or list of length N for T1 amplitude-damping rates.
            t_max, dt: time grid.
        """
        if receiver.N != self.N:
            raise ValueError(
                f"receiver.N ({receiver.N}) does not match chain.N ({self.N})"
            )

        if terms is not None:
            bilinear = [(t[0], t[1], self.J) for t in terms]
            H = _build_bilinear(self.N, self.bonds, bilinear)
        else:
            H = self.H

        gamma_l = [self.gamma_0] * self.N
        if gamma_t1 is None:
            gamma_t1_l = [0.0] * self.N
        elif np.isscalar(gamma_t1):
            gamma_t1_l = [float(gamma_t1)] * self.N
        else:
            gamma_t1_l = list(gamma_t1)
            if len(gamma_t1_l) != self.N:
                raise ValueError(
                    f"gamma_t1 list length {len(gamma_t1_l)} != N {self.N}"
                )

        return cockpit_panel(
            H, gamma_l, receiver.rho, self.N,
            gamma_t1_l=gamma_t1_l,
            t_max=t_max, dt=dt,
            threshold=threshold, cluster_tol=cluster_tol,
        )
