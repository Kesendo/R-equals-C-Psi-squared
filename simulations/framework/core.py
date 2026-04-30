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
    HARDWARE_DISSIPATORS,
    HARDWARE_DISSIPATOR_D1,
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


def _pauli_pair_is_truly(a, b):
    """Syntactic per-term truly check for a single Pauli-pair bilinear.

    A bond bilinear (a, b) preserves the Π-palindrome (M = 0) iff
    a == b (any of II, XX, YY, ZZ) or {a, b} ⊆ {I, X}. Verified
    against classify_pauli_pair for all 16 single-term pairs.
    """
    if a == b:
        return True
    if a in {'I', 'X'} and b in {'I', 'X'}:
        return True
    return False


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
        from .lindblad import lindbladian_z_plus_t1
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

    def predict_residual_norm_squared_from_terms(self, terms, gamma_t1=None):
        """Closed-form ‖M‖² from terms via per-term Frobenius identity (no L computed).

        Verified empirically across chain/ring/star/K_N at N=3..6 plus the
        full V-Effect 36-combos enumeration at N=3:

            ‖M(L_Z)‖²    = sum_k [ 2^(N+2) · n_YZ_k · ‖H_k‖²_F · (1 if non-truly else 0) ]
            ‖M(L_Z+T1)‖² = ‖M(L_Z)‖² + 4^(N−1) · [ 3·Σ γT1² + 4·(Σ γT1)² ]

        Each Pauli-pair term k contributes independently because distinct
        Pauli strings are orthogonal in Frobenius, AND truly-class single
        terms (a==b, or {a,b}⊆{I,X}) contribute nothing because M_k=0
        algebraically. The V-Effect cases YY+YZ, YY+ZY, YZ+ZZ, ZY+ZZ
        enforce per-term-truly handling: gross-list classify is not enough.

        n_YZ_k = number of Y/Z letters in term k (bit_b-odd letters, which
        are the Π-symmetry-breaking ones).

        T1 contribution is Hamiltonian-independent, γ_Z-independent, and
        orthogonal to the Hamiltonian palindrome residual (no cross-term).
        Verified at N=3..6 for arbitrary {γ_T1_l} distributions.

        Topology enters only via ‖H_k‖²_F (cheap to compute). The graph-
        dependent c_H scaling of `predict_residual_norm_squared` is the
        chain/ring/star special case where Pauli strings are uniquely-bonded;
        this method is universal.

        Args:
            terms: list of (a, b) Pauli letter tuples.
            gamma_t1: optional T1 amplitude-damping rates. Scalar (uniform
                      across sites), list of length N, or None / 0 (no T1).
        """
        # Hamiltonian palindrome part. Decompose into:
        #   1. Drop per-term truly-class terms (they contribute M=0).
        #   2. Group remaining terms by n_YZ-class (Π-symmetry-breaking count).
        #   3. Within each n_YZ-class, Pauli strings can overlap (e.g.
        #      (I,Y) and (Y,I) both place Y on the middle site of a chain),
        #      so the per-class Frobenius norm is the norm of the combined
        #      sub-Hamiltonian, not a sum of per-term norms.
        #   4. Across n_YZ-classes M-residuals are orthogonal — sum the
        #      per-class predictions.
        from collections import defaultdict
        n_yz_groups = defaultdict(list)
        for (a, b) in terms:
            if _pauli_pair_is_truly(a, b):
                continue
            n_yz_k = sum(1 for L in (a, b) if L in 'YZ')
            n_yz_groups[n_yz_k].append((a, b))
        z_part = 0.0
        for n_yz, group_terms in n_yz_groups.items():
            H_group = _build_bilinear(
                self.N, self.bonds, [(a, b, self.J) for (a, b) in group_terms]
            )
            H_group_frob_sq = float(np.real(np.trace(H_group.conj().T @ H_group)))
            z_part += (2 ** (self.N + 2)) * n_yz * H_group_frob_sq

        # T1 dissipator contribution (Hamiltonian-independent, additive)
        if gamma_t1 is None:
            t1_part = 0.0
        else:
            if np.isscalar(gamma_t1):
                gamma_t1_l = [float(gamma_t1)] * self.N
            else:
                gamma_t1_l = [float(g) for g in gamma_t1]
                if len(gamma_t1_l) != self.N:
                    raise ValueError(
                        f"gamma_t1 list length {len(gamma_t1_l)} != N {self.N}"
                    )
            sum_g2 = sum(g * g for g in gamma_t1_l)
            sum_g = sum(gamma_t1_l)
            t1_part = (4 ** (self.N - 1)) * (3 * sum_g2 + 4 * sum_g * sum_g)

        return z_part + t1_part

    def predict_M_spectrum_pi2_odd(self, terms, c=1.0):
        """F80 structural-identity prediction of M's spectrum (chain Π²-odd 2-body).

        Theorem F80 (verified bit-exact at N=3..7, all four pure Π²-odd Pauli
        pairs and their sums, see PROOF_F80_BLOCH_SIGNWALK):

            Spec(M)_{nontrivial} = { ±2i · λ : λ ∈ Spec(H_non-truly) }
            mult_M(±2i·λ) = mult_H_non-truly(λ) · 2^N

        H_non-truly is the input Hamiltonian with truly-class bilinears
        (XX, YY, ZZ, or any (a, b) ⊆ {I, X}) dropped, since they contribute
        M_term = 0 and are invisible to the residual.

        The remaining bilinears must all be Π²-odd 2-body: P, Q ∈ {X, Y, Z}
        with bit_b(P) + bit_b(Q) ≡ 1 (mod 2), i.e., one of (X,Y), (X,Z),
        (Y,X), (Z,X). Π²-even non-truly bilinears (Y,Z) and (Z,Y) are not in
        F80's verified scope and are rejected.

        Args:
            terms: list of (a, b) Pauli letter tuples. Truly terms are dropped;
                Π²-odd 2-body terms are summed into H_non-truly.
            c: bond coupling magnitude (default 1.0); each term is scaled by c.

        Returns:
            dict {eigenvalue (complex, purely imaginary): multiplicity (int)}
            for M's spectrum. If H_non-truly = 0 (all terms truly), returns
            {0+0j: 4^N}; the trivially-paired truly case.

        Raises:
            ValueError: if any non-truly term is not Π²-odd 2-body, or
                contains an identity letter (single-body falls under F78).
        """
        from collections import Counter
        # Filter truly, validate remaining as Π²-odd 2-body
        non_truly_terms = []
        for (a, b) in terms:
            if _pauli_pair_is_truly(a, b):
                continue
            if 'I' in (a, b):
                raise ValueError(
                    f"term ({a},{b}) contains identity (single-body); "
                    "F80 covers chain Π²-odd 2-body only; see F78 for single-body"
                )
            from .pauli import bit_b, _resolve
            ab_idx = _resolve(a)
            bb_idx = _resolve(b)
            parity = (bit_b(ab_idx) + bit_b(bb_idx)) % 2
            if parity != 1:
                raise ValueError(
                    f"term ({a},{b}) is Π²-even non-truly (parity 0); "
                    "F80 covers Π²-odd only; Π²-even non-truly clusters "
                    "are richer and not in F80's verified scope"
                )
            non_truly_terms.append((a, b, c))

        if not non_truly_terms:
            return {0 + 0j: 4 ** self.N}

        # Build H_non-truly and compute many-body spectrum
        H_nt = _build_bilinear(self.N, self.bonds, non_truly_terms)
        H_evs = np.linalg.eigvalsh(H_nt)
        h_counts = Counter(round(float(e), 10) for e in H_evs)

        # Translate to M-spectrum: Spec(M) = ±2i · Spec(H), mult ×2^N
        # H is Hermitian so eigenvalues are real; ±2i·λ for each H eigenvalue.
        # Note: H already has ±λ pairs (particle-hole symmetry of Majorana
        # bilinear), so 2i·λ already carries the ± structure.
        bra_factor = 2 ** self.N
        m_spectrum = {}
        for h_ev, h_mult in h_counts.items():
            m_ev = 2j * h_ev
            m_spectrum[m_ev] = h_mult * bra_factor
        return m_spectrum

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
            h_part = self.predict_residual_norm_squared_from_terms(terms)
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
            'description': 'Super-operator palindrome trichotomy (truly/soft/hard) tomographically distinguishable on Heron r2 hardware at N=3. Hardware Δ(soft − truly) = -0.722 matches the Trotter n=3 prediction to 0.0014, NOT the continuous-Lindblad idealization (Δ = -0.623). Original interpretation that T1 amplifies the soft-break is REFUTED: T1 monotonically attenuates |Δ| (γ_T1=0.5 gives Δ=-0.44). The hardening is Trotter discretization at δt=0.267 with ‖H‖_op = 2.83·J, where ‖H·δt‖ ≈ 0.76 violates the small-step regime. See _marrakesh_t1_amplification_test.py.',
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
                              'γ_Z = 0.1 exactly via Δ matching to 0.0014 (see _marrakesh_t1_amplification_test).',
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
