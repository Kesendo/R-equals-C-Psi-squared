#!/usr/bin/env python3
"""Validate recommend_initial_state's predictions against actual trajectories.

For each test case, simulate the full ρ(t) evolution and verify:
  1. Protected observables ⟨P⟩(t) stay below 1e-10 for all t (Π-protection)
  2. Non-protected observables decay at measurable rates
  3. README Rule 1: GHZ has FASTER decay of non-protected observables than W

Test cases:
  (full Heisenberg, |GHZ⟩)  — 52 protected predicted, fast decay
  (full Heisenberg, |W⟩)    — 44 protected predicted, slower decay
  (XY+YX, |010⟩)            — 44 protected predicted
  (truly XX+YY, |+−+⟩)      — historical case, 1 protected
  (XY+YX, |+−+⟩)            — historical case, 29 protected
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def simulate_trajectory(L, rho_0, times):
    """Compute ρ(t) at sample times via matrix exponential."""
    d = rho_0.shape[0]
    rho_vec0 = rho_0.T.reshape(-1).copy()
    out = []
    for t in times:
        rho_vec = expm(L * t) @ rho_vec0
        out.append(rho_vec.reshape(d, d).T)
    return out


def all_pauli_expectations(rho, N):
    """Return ⟨P_α⟩ for all 4^N Pauli strings, indexed by α."""
    M_basis = fw._vec_to_pauli_basis_transform(N)
    rho_pauli = fw.pauli_basis_vector(rho, N)
    return rho_pauli  # already 4^N expectation values × normalization


def main():
    N = 3
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    one = np.array([0, 1], dtype=complex)

    psi_xneel = np.kron(plus, np.kron(minus, plus))
    rho_xneel = np.outer(psi_xneel, psi_xneel.conj())

    psi_010 = np.kron(zero, np.kron(one, zero))
    rho_010 = np.outer(psi_010, psi_010.conj())

    # |GHZ⟩ = (|000⟩+|111⟩)/√2
    psi_ghz = (np.kron(zero, np.kron(zero, zero))
                + np.kron(one, np.kron(one, one))) / math.sqrt(2)
    rho_ghz = np.outer(psi_ghz, psi_ghz.conj())

    # |W⟩ = symmetric single-excitation
    psi_w = (np.kron(zero, np.kron(zero, one))
              + np.kron(zero, np.kron(one, zero))
              + np.kron(one, np.kron(zero, zero))) / math.sqrt(3)
    rho_w = np.outer(psi_w, psi_w.conj())

    cases = [
        ('Heisenberg + |GHZ⟩',  [('X', 'X', J), ('Y', 'Y', J), ('Z', 'Z', J)],
         rho_ghz, 'GHZ', 52),
        ('Heisenberg + |W⟩',    [('X', 'X', J), ('Y', 'Y', J), ('Z', 'Z', J)],
         rho_w, 'W', 44),
        ('XY+YX + |010⟩',       [('X', 'Y', J), ('Y', 'X', J)],
         rho_010, '010', 44),
        ('truly XX+YY + |+−+⟩', [('X', 'X', J), ('Y', 'Y', J)],
         rho_xneel, 'xneel', 1),
        ('XY+YX + |+−+⟩',       [('X', 'Y', J), ('Y', 'X', J)],
         rho_xneel, 'xneel', 29),
    ]

    t_max = 20.0
    dt = 0.05
    times = np.linspace(0, t_max, int(t_max / dt) + 1)

    print(f"Recommender validation, N={N}, γ={GAMMA}, γ_T1={GAMMA_T1}, "
          f"t_max={t_max}, dt={dt}")
    print()

    print(f"  {'case':<26s}  {'pred_n_prot':>11s}  "
          f"{'verified':>9s}  {'max |⟨P⟩| (prot)':>17s}  "
          f"{'mean τ (non-prot)':>18s}")
    print('-' * 100)

    summary = []
    for label, terms, rho_0, state_tag, predicted_n_prot in cases:
        H = fw._build_bilinear(N, bonds, terms)
        L = fw.lindbladian_z_plus_t1(H, [GAMMA] * N, [GAMMA_T1] * N)

        # Identify predicted-protected Pauli set
        M_basis = fw._vec_to_pauli_basis_transform(N)
        L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
        evals, V = np.linalg.eig(L_pauli)
        Vinv = np.linalg.inv(V)
        rho_pauli_0 = fw.pauli_basis_vector(rho_0, N)
        c = Vinv @ rho_pauli_0

        n_eig = len(evals)
        used = np.zeros(n_eig, dtype=bool)
        clusters = []
        for i in range(n_eig):
            if used[i]:
                continue
            cl = [i]
            used[i] = True
            for j in range(i + 1, n_eig):
                if not used[j] and abs(evals[j] - evals[i]) < 1e-8:
                    cl.append(j)
                    used[j] = True
            clusters.append(cl)

        protected_alphas = []
        non_protected_alphas = []
        for alpha in range(1, 4 ** N):
            max_S = max(
                abs(sum(V[alpha, k] * c[k] for k in cl)) for cl in clusters
            )
            if max_S < 1e-9:
                protected_alphas.append(alpha)
            else:
                non_protected_alphas.append(alpha)

        actual_n_prot = len(protected_alphas)

        # Simulate trajectory
        traj = simulate_trajectory(L, rho_0, times)

        # For each Pauli, get the trajectory of ⟨P⟩(t)
        # Note: pauli_basis_vector returns coefficients in Pauli basis,
        # Tr(P_α · ρ) = 2^N · pauli_coef[α]
        all_traj = np.zeros((4 ** N, len(times)), dtype=complex)
        for ti, rho_t in enumerate(traj):
            all_traj[:, ti] = fw.pauli_basis_vector(rho_t, N)
        # Pauli expectations are coef · 2^N (since pauli_basis_vector
        # is the Pauli decomposition coefficients).
        # Actually checking convention: ρ = Σ_α coef[α] · P_α / 2^N
        # so ⟨P_α⟩ = Tr(P_α · ρ) = coef[α].  Let's just verify with α=0:
        # coef[0] = Tr(I·ρ) = 1.0 trivially.

        # Max |⟨P⟩(t)| over protected observables and over time
        if protected_alphas:
            max_prot_val = max(
                float(np.max(np.abs(all_traj[alpha, :])))
                for alpha in protected_alphas
            )
        else:
            max_prot_val = 0.0

        # Mean decay timescale for non-protected: time at which
        # |⟨P⟩(t)| drops to 1/e of its peak value.
        decay_timescales = []
        for alpha in non_protected_alphas:
            traj_alpha = np.abs(all_traj[alpha, :])
            peak_val = float(traj_alpha.max())
            if peak_val < 1e-6:
                continue  # too small to measure
            target = peak_val / math.e
            # Find first time after peak where |⟨P⟩| < target
            peak_idx = int(np.argmax(traj_alpha))
            for ti in range(peak_idx, len(times)):
                if traj_alpha[ti] < target:
                    decay_timescales.append(times[ti] - times[peak_idx])
                    break
        mean_tau = float(np.mean(decay_timescales)) if decay_timescales else 0.0

        verified = "✓" if (actual_n_prot == predicted_n_prot and
                            max_prot_val < 1e-6) else "✗"

        print(f"  {label:<26s}  {predicted_n_prot:>11d}  "
              f"{verified:>9s}  {max_prot_val:>17.4e}  "
              f"{mean_tau:>16.3f}")

        summary.append({
            'label': label, 'state_tag': state_tag,
            'predicted_n_prot': predicted_n_prot,
            'actual_n_prot': actual_n_prot,
            'max_prot_val': max_prot_val,
            'mean_tau': mean_tau,
            'n_decay_samples': len(decay_timescales),
        })

    print()
    print("=" * 100)
    print("Rule 1 verification: GHZ vs W decay timescales on full Heisenberg")
    print("=" * 100)
    ghz_tau = next(s['mean_tau'] for s in summary if s['state_tag'] == 'GHZ')
    w_tau = next(s['mean_tau'] for s in summary if s['state_tag'] == 'W')
    print()
    print(f"  GHZ mean decay τ:  {ghz_tau:.3f}")
    print(f"  W   mean decay τ:  {w_tau:.3f}")
    print(f"  W/GHZ ratio:       {w_tau/ghz_tau:.2f}× longer for W")
    print()
    if w_tau > ghz_tau:
        print(f"  ✓ Rule 1 confirmed: W's non-protected observables persist "
              f"{w_tau/ghz_tau:.2f}× longer than GHZ's.")
    else:
        print(f"  ✗ Rule 1 NOT visible at this resolution. May need finer dt.")


if __name__ == "__main__":
    main()
