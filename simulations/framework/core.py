"""Cockpit core: ChainSystem + Receiver + Confirmations.

The OOP cockpit layer that composes framework primitives into callable
shortcuts. Each method replaces document-reading with a single call.

  ChainSystem(N, gamma_0, J, topology, H_type)   — chain setup with cached H/L
  Receiver(psi, chain=None)                       — F71-aware state wrapper
  Confirmations                                   — hardware-confirmed predictions lookup
"""
from __future__ import annotations

import warnings

import numpy as np

from .lindblad import (
    lindbladian_z_dephasing,
    palindrome_residual,
    palindrome_residual_norm_squared_factor_graph,
)
from .pauli import _build_bilinear, pauli_matrix
from .symmetry import f71_eigenstate_class, receiver_engineering_signature
from .lebensader import cockpit_panel


def _site_op_kron(op, site, N):
    """Helper: 2×2 op on `site`, identity elsewhere."""
    I2 = np.eye(2, dtype=complex)
    factors = [I2] * N
    factors[site] = op
    out = factors[0]
    for f in factors[1:]:
        out = np.kron(out, f)
    return out


# ----------------------------------------------------------------------
# ChainSystem
# ----------------------------------------------------------------------

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

    Methods:
        classify_pauli_pair(terms, J_scale=1.0) → 'truly' | 'soft' | 'hard'
        predict_residual_norm_squared(c_H, class) → ‖M‖² closed form
        cockpit_panel(receiver, terms=None, gamma_t1=None, ...)
            → full Lebensader analysis (skeleton + trace + cusp + chiral + Y-parity)
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

    def classify_pauli_pair(self, terms, J_scale=1.0, op_tol=1e-10, spec_tol=1e-6):
        """Trichotomy classification (truly / soft / hard) for a Pauli-pair H.

        Builds H from the given terms on this chain's bonds, computes the
        palindrome residual M and the spectrum-pairing error.

        Args:
            terms: list of (a, b) letter tuples, e.g. [('X','X'), ('Y','Y')].
            J_scale: bond coupling for the test Hamiltonian (default 1.0).

        Returns:
            'truly' | 'soft' | 'hard'
        """
        bilinear = [(t[0], t[1], J_scale) for t in terms]
        H_test = _build_bilinear(self.N, self.bonds, bilinear)
        L_test = lindbladian_z_dephasing(H_test, [self.gamma_0] * self.N)
        Sigma_gamma = self.N * self.gamma_0
        M = palindrome_residual(L_test, Sigma_gamma, self.N)
        op_norm = float(np.linalg.norm(M))
        evals = np.linalg.eigvals(L_test)
        used = np.zeros(len(evals), dtype=bool)
        max_err = 0.0
        for i in range(len(evals)):
            if used[i]:
                continue
            target = -evals[i] - 2 * Sigma_gamma
            dists = np.abs(evals - target)
            for j in range(len(evals)):
                if used[j]:
                    dists[j] = np.inf
            best_j = int(np.argmin(dists))
            if best_j != i:
                used[i] = True
                used[best_j] = True
            else:
                used[i] = True
            max_err = max(max_err, float(dists[best_j]))
        if op_norm < op_tol:
            return 'truly'
        if max_err < spec_tol:
            return 'soft'
        return 'hard'

    def predict_residual_norm_squared(self, c_H, hamiltonian_class='main'):
        """Closed-form ‖M(N, G)‖² = c_H · F(N, G) without computing M.

        Uses palindrome_residual_norm_squared_factor_graph with this chain's
        topology invariants (B, D2). For chains, B = N-1 and D2 = 4N - 6.
        """
        factor = palindrome_residual_norm_squared_factor_graph(
            self.N, self.B, self.D2, hamiltonian_class
        )
        return c_H * factor

    def residual_norm_squared(self, terms, J_scale=None):
        """Numerical ‖M‖² for a Pauli-pair Hamiltonian on this chain.

        Builds H, the Z-dephasing Lindbladian, and the palindrome residual M;
        returns ‖M‖²_F. For 'truly' classes, ≈ 0 to floating-point precision.

        Replaces the recurring 4-line boilerplate
        (_build_bilinear → lindbladian_z_dephasing → palindrome_residual → norm).

        Args:
            terms: list of (a, b) letter tuples.
            J_scale: optional override for chain.J (default uses self.J).
        """
        J = self.J if J_scale is None else float(J_scale)
        bilinear = [(a, b, J) for (a, b) in terms]
        H = _build_bilinear(self.N, self.bonds, bilinear)
        L = lindbladian_z_dephasing(H, [self.gamma_0] * self.N)
        Sigma_gamma = self.N * self.gamma_0
        M = palindrome_residual(L, Sigma_gamma, self.N)
        return float(np.linalg.norm(M) ** 2)

    def predict_residual_norm_squared_from_terms(self, terms, is_truly=None):
        """Closed-form ‖M‖² from terms via Frobenius identity (no L computed).

        Verified empirically across chain/ring/star/K_N at N=3..6:

            ‖M‖²_F = 2^(N+2) · n_YZ · ‖H‖²_F   (homogeneous non-truly H)
                   = 0                            (truly H)

        n_YZ = number of Y/Z letters per Pauli-pair term (the bit_b-odd letters,
        which are the Π-symmetry-breaking ones). Term list must be homogeneous
        in n_YZ; mixed-class lists must be split and added per-class.

        Topology enters only via ‖H‖²_F (cheap to compute). The graph-dependent
        c_H scaling of `predict_residual_norm_squared` is the chain/ring/star
        special case where Pauli strings are uniquely-bonded; this method is
        universal.

        Args:
            terms: Pauli-pair list (homogeneous in Y/Z count per term).
            is_truly: optional override. If None, calls classify_pauli_pair.
                      Pass True/False to skip the numerical classification when
                      the truly status is known (e.g., from prior analysis).

        Raises:
            ValueError: if the term list is not homogeneous in n_YZ_per_term.
        """
        if not terms:
            return 0.0
        # truly check first — truly Hamiltonians (e.g. Heisenberg XX+YY+ZZ) can
        # mix n_YZ per term (0+2+2) and still be palindrome-respecting; the
        # homogeneity rule only applies to the non-truly Frobenius branch.
        if is_truly is None:
            is_truly = (self.classify_pauli_pair(terms) == 'truly')
        if is_truly:
            return 0.0
        n_yz_per_term = [sum(1 for L in (a, b) if L in 'YZ') for (a, b) in terms]
        if len(set(n_yz_per_term)) > 1:
            raise ValueError(
                f"Term list is not homogeneous in Y/Z count per term "
                f"(got {n_yz_per_term}). Split into homogeneous parts and add "
                f"the predictions: ||M_total||^2 = sum_k predict(...)_k."
            )
        n_yz = n_yz_per_term[0]
        bilinear = [(a, b, self.J) for (a, b) in terms]
        H = _build_bilinear(self.N, self.bonds, bilinear)
        H_frob_sq = float(np.real(np.trace(H.conj().T @ H)))
        return (2 ** (self.N + 2)) * n_yz * H_frob_sq

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


# ----------------------------------------------------------------------
# Receiver
# ----------------------------------------------------------------------

class Receiver:
    """A quantum state vector with framework-aware F71 classification.

    Args:
        psi: complex array of length 2^N
        chain: optional ChainSystem (must match N if given)

    Properties (cached):
        N, f71_class, rho

    Methods:
        signature() → receiver_engineering_signature dict
    """

    def __init__(self, psi, chain=None, atol=1e-6):
        psi_arr = np.asarray(psi, dtype=complex)
        if psi_arr.ndim != 1:
            raise ValueError(
                f"psi must be a 1D state vector; got shape {psi_arr.shape}. "
                f"To wrap a density matrix, use Receiver.from_rho(rho, chain) instead."
            )
        self.psi = psi_arr
        self.N = int(round(np.log2(len(self.psi))))
        if 2 ** self.N != len(self.psi):
            raise ValueError(f"psi length {len(self.psi)} is not a power of 2")
        norm = float(np.linalg.norm(self.psi))
        if not np.isclose(norm, 1.0, atol=atol):
            raise ValueError(
                f"psi must be normalized (||psi||=1); got ||psi||={norm:.6g}. "
                f"Either normalize before passing or use Receiver.from_psi_unnormalized(...)."
            )
        if chain is not None and chain.N != self.N:
            raise ValueError(f"chain.N ({chain.N}) does not match psi N ({self.N})")
        self.chain = chain
        self._f71_class_cache = "unset"
        self._rho_cache = None

    @classmethod
    def from_psi_unnormalized(cls, psi, chain=None):
        """Wrap a non-normalized state vector by normalizing first."""
        psi = np.asarray(psi, dtype=complex).ravel()
        n = float(np.linalg.norm(psi))
        if n == 0.0:
            raise ValueError("psi is the zero vector; cannot normalize.")
        return cls(psi / n, chain=chain)

    @property
    def f71_class(self):
        if self._f71_class_cache == "unset":
            self._f71_class_cache = f71_eigenstate_class(self.psi)
        return self._f71_class_cache

    @property
    def rho(self):
        if self._rho_cache is None:
            rho = np.outer(self.psi, self.psi.conj())
            self._rho_cache = (rho + rho.conj().T) / 2.0
        return self._rho_cache

    def signature(self):
        return receiver_engineering_signature(self.psi)


# ----------------------------------------------------------------------
# Confirmations — hardware-confirmed framework predictions
# ----------------------------------------------------------------------

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
            'predicted_value': {'truly': 0.000, 'soft': -0.623, 'hard': +0.195,
                                'delta_soft_minus_truly': -0.623},
            'measured_value': {'truly': +0.011, 'soft': -0.711, 'hard': +0.205,
                               'delta_soft_minus_truly': -0.722},
            'hardware_data': 'data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json',
            'experiment_doc': 'experiments/V_EFFECT_FINE_STRUCTURE.md',
            'framework_primitive': 'palindrome_residual + classify_pauli_pair',
            'description': 'Super-operator palindrome trichotomy (truly/soft/hard) tomographically distinguishable on Heron r2 hardware at N=3. T1/T2 noise actually amplifies the soft-break signal.',
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
        'lebensader_skeleton_trace_decoupling': {
            'date': '2026-04-26',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7n3013aq2pc73a2a18g + d7n3eqqt99kc73d34qtg',
            'observable': 'Π-protected counts (skeleton) + θ-trajectory tails (trace)',
            'predicted_value': 'Skeleton/trace decouple at N≥4, co-occur at N=3. T1 makes drop measurable: bond-flipped Z-free pairs (XY+YX, IY+YI) preserve skeleton (drop≤1) and trace (long θ-tail); Z-containing soft pairs (YZ+ZY, XZ+ZX) collapse both (drop≈28-29, no tail).',
            'measured_value': 'drop=28 for YZ+ZY confirmed on Marrakesh. Pearson(drop, Δ∫θ)=+0.85. Bures velocity gives no third discriminator.',
            'hardware_data': 'external (raw JSON not in repo; tables of <X_0 I Z_2> per t and basis documented inline in experiment_doc)',
            'experiment_doc': 'review/EMERGING_QUESTIONS.md',
            'framework_primitive': 'cockpit_panel — composes pi_protected_observables + θ-trajectory',
            'description': 'Lebensader as Stromkabel (EQ-030 closure): skeleton (Π-protected algebraic count) and trace (θ-geometric tail) are NOT two discriminators but one bridge held together by Π·L·Π⁻¹ + L + 2Σγ·I = 0. Hardware-confirmed across 3 of 4 bond-flipped Z-free corners; Bures velocity confirmed null as third axis. ChainSystem.cockpit_panel gives this in one call.',
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
