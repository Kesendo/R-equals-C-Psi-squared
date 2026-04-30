#!/usr/bin/env python3
"""EQ-030 attack: how does T1 amplitude damping shift the Π-protected set?

Hardware finale (Marrakesh, Kingston, Fez) measured Δ⟨X₀Z₂⟩ for the
soft-broken Hamiltonian J(XY+YX) at 1.14-1.49× the continuous-Lindblad
idealization. The amplification was originally attributed to T1 amplitude
damping (plus ZZ-crosstalk). The 2026-04-30 follow-up
(_marrakesh_t1_amplification_test.py) refuted that attribution for the
Marrakesh soft-break run: T1 monotonically attenuates the soft signal;
Trotter n=3 discretization at δt=0.267 fully accounts for the hardening.
The script below remains valid as a separate question: even though T1 does
not amplify ⟨X₀Z₂⟩, it does break the strict Π-palindrome and shifts the
Π-protected count. framework.py originally only handles pure Z-dephasing
where the palindrome holds. With the `lindbladian_z_plus_t1` primitive,
we can ask: for the soft case, how does the protected-set evolve as
T1-rate is added?

Specifically:
  - At T1-rate = 0 (pure Z-dephasing): pi_protected_observables on the soft
    case finds 252 protected (e.g. for J(XY+YX) on |+−+⟩, N=3).
  - As T1-rate increases, the spectrum of L changes; the eigenmodes that
    cancel within their degenerate clusters may drift apart, breaking
    cluster-cancellation conditions for additional Pauli observables.
  - The number of protected observables should DECREASE as T1 grows.

This script sweeps T1-rate as a fraction of γ_dephasing and reports
the protected count for the soft case.
"""
import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def main():
    N = 3
    GAMMA_DEPH = 0.1
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    # |+−+⟩ initial state (Snapshot D's init)
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly  J(XX+YY)', [('X', 'X', J), ('Y', 'Y', J)]),
        ('soft   J(XY+YX)', [('X', 'Y', J), ('Y', 'X', J)]),
        ('hard   J(XX+XY)', [('X', 'X', J), ('X', 'Y', J)]),
    ]

    # T1 rates as multiples of γ_deph
    t1_ratios = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]

    print(f"Π-protected observable count vs T1 amplitude damping  (N={N}, γ_deph={GAMMA_DEPH}, |+−+⟩)")
    print()
    print(f"  γ_T1 / γ_deph   ", "    ".join(f"{r:>8.2f}" for r in t1_ratios))
    print(f"  {'Hamiltonian':<20s}  ", "    ".join(f"{'protected':>8s}" for _ in t1_ratios))
    print(f"  {'-' * 20}  ", "    ".join(f"{'-' * 8}" for _ in t1_ratios))

    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        protected_counts = []
        op_norms = []
        for ratio in t1_ratios:
            gamma_t1 = ratio * GAMMA_DEPH
            L = fw.lindbladian_z_plus_t1(
                H,
                gamma_l=[GAMMA_DEPH] * N,
                gamma_t1_l=[gamma_t1] * N,
            )
            # Build the Pauli-basis transform once
            M_basis = fw._vec_to_pauli_basis_transform(N)
            L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)

            # Eigendecomposition + cluster-sum identification (as in pi_protected_observables)
            evals, V = np.linalg.eig(L_pauli)
            Vinv = np.linalg.inv(V)
            rho_pauli = fw.pauli_basis_vector(rho_0, N)
            c = Vinv @ rho_pauli

            # Cluster eigenvalues
            cluster_tol = 1e-8
            n_eig = len(evals)
            used = np.zeros(n_eig, dtype=bool)
            clusters = []
            for i in range(n_eig):
                if used[i]:
                    continue
                cl = [i]
                used[i] = True
                for j in range(i + 1, n_eig):
                    if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                        cl.append(j)
                        used[j] = True
                clusters.append(cl)

            # Count protected
            threshold = 1e-9
            n_protected = 0
            for alpha in range(1, 4 ** N):
                max_S = 0.0
                for cl in clusters:
                    S = sum(V[alpha, k] * c[k] for k in cl)
                    max_S = max(max_S, abs(S))
                if max_S < threshold:
                    n_protected += 1
            protected_counts.append(n_protected)

            # Also note the operator residual against the pure-Z-dephasing
            # palindrome center (2*Σγ_deph). The T1 term breaks this.
            M_residual = fw.palindrome_residual(L, N * GAMMA_DEPH, N)
            op_norms.append(float(np.linalg.norm(M_residual)))

        counts_str = "    ".join(f"{c:>8d}" for c in protected_counts)
        print(f"  {label:<20s}  {counts_str}")

    print()
    print(f"Pure Z-dephasing baseline ('γ_T1/γ_deph = 0'):")
    print(f"  truly: 32 protected (XIZ, ZIX in protected set; matches Snapshot C)")
    print(f"  soft:  30 protected (XIZ, ZIX OUT of protected set; soft-break leakage)")
    print()
    print(f"What to read from the sweep:")
    print(f"  As T1 grows, do additional observables leak out of 'protected'?")
    print(f"  truly's protected count should drop (T1 itself is not Π-symmetric)")
    print(f"  soft's protected count should also drop, but the relative gap")
    print(f"  to truly tells us how T1 changes the soft-vs-truly discrimination.")
    print()
    print(f"Operator residual norms ‖Π·L·Π⁻¹ + L + 2Σγ_deph·I‖ (for reference):")
    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        op_norms = []
        for ratio in t1_ratios:
            gamma_t1 = ratio * GAMMA_DEPH
            L = fw.lindbladian_z_plus_t1(H, [GAMMA_DEPH]*N, [gamma_t1]*N)
            M_residual = fw.palindrome_residual(L, N*GAMMA_DEPH, N)
            op_norms.append(float(np.linalg.norm(M_residual)))
        norms_str = "    ".join(f"{o:>8.2e}" for o in op_norms)
        print(f"  {label:<20s}  {norms_str}")


if __name__ == "__main__":
    main()
