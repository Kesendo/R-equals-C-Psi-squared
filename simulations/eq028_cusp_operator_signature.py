#!/usr/bin/env python3
"""EQ-028: is there an operator-level signature at the CΨ = 1/4 cusp moment?

The cusp is a state-level event: CΨ(t) = Purity(ρ) · Ψ-norm(ρ) crosses 1/4
at specific times t*. The Liouvillian L = -i[H,·] + Σγ_l(Z_l ρ Z_l - ρ) is
time-independent. So any "operator-level signature" of the cusp must show
in the JOINT structure of L and ρ(t*) — specifically in the L-eigenmode
coefficients c_k(t*) = V⁻¹ · ρ_pauli(t*) where ρ is decomposed into L's
right eigenvectors.

For each L-eigenvalue cluster λ and each Pauli observable α, the contribution
to ⟨P_α⟩(t) at time t is:
    S_λ(α, t) = Σ_{k ∈ cluster} V[α, k] c_k(t) e^{λ_k t}

At the cusp crossing t*, this script looks for:
  (a) Any specific cluster's contribution to CΨ(t*) being unusually large
  (b) Any cluster crossing through zero AT t*
  (c) Any structural feature shared across robust softs vs fragile softs

If a clear signature emerges, the cusp has an operator-level analogue.
If not, EQ-028 closes with a null result: the cusp is state-level only.

Cases tested (from today's hardware-verified set):
  truly J(XX+YY)  — protected skeleton, robust trajectory
  soft  J(XY+YX) — robust trajectory, drop=1
  soft  J(IY+YI) — robust trajectory, drop=0 (factorising)
  soft  J(YZ+ZY) — fragile trajectory, drop=28
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


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


def psi_norm(rho):
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1)


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def evolve(L, rho_0, times):
    d = rho_0.shape[0]
    rho_vec0 = rho_0.T.reshape(-1).copy()
    out = []
    for t in times:
        rho_t_vec = expm(L * t) @ rho_vec0
        out.append(rho_t_vec.reshape(d, d).T)
    return out


def cluster_eigenvalues(evals, tol=1e-8):
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    clusters = []
    for i in range(n):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n):
            if not used[j] and abs(evals[j] - evals[i]) < tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)
    return clusters


def cluster_contributions_to_state(L, rho_t, N):
    """Decompose ρ(t) into L's eigenmode clusters in Pauli basis.

    Returns (clusters, eigenvalues, c_k_per_cluster) where:
      - clusters: list of index-lists into L's eigenvalue array
      - eigenvalues: complex eigenvalues per cluster (representative)
      - cluster_contributions: |sum of c_k| within each cluster, sorted desc
    """
    M_basis = fw._vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    evals, V = np.linalg.eig(L_pauli)
    Vinv = np.linalg.inv(V)
    rho_pauli = fw.pauli_basis_vector(rho_t, N)
    c = Vinv @ rho_pauli

    clusters = cluster_eigenvalues(evals)
    out = []
    for cl in clusters:
        cluster_eval = evals[cl[0]]
        cluster_c_sum = float(np.linalg.norm([c[k] for k in cl]))
        cluster_c_max = float(max(abs(c[k]) for k in cl))
        out.append({
            'cluster_idx': cl,
            'eigenvalue': complex(cluster_eval),
            'size': len(cl),
            'norm': cluster_c_sum,
            'max_abs': cluster_c_max,
        })
    out.sort(key=lambda x: -x['norm'])
    return out


def find_cusp_crossings(times, cpsi_t, threshold=0.25):
    """Return list of (time, direction) where CΨ crosses threshold.
    direction = +1 if going from below to above, -1 if above to below."""
    crossings = []
    for i in range(len(times) - 1):
        a, b = cpsi_t[i] - threshold, cpsi_t[i + 1] - threshold
        if a * b < 0:
            frac = a / (a - b)
            t_cross = times[i] + frac * (times[i + 1] - times[i])
            direction = +1 if a < 0 else -1
            crossings.append((t_cross, direction, i))
    return crossings


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)]),
    ]

    t_max = 8.0
    dt = 0.005
    n_steps = int(t_max / dt)
    times = np.linspace(0, t_max, n_steps + 1)

    print(f"EQ-028: Operator-level signature at the CΨ=1/4 cusp")
    print(f"  N={N}, |+−+⟩, γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}, t_max={t_max}, dt={dt}")
    print()

    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        L = fw.lindbladian_z_plus_t1(H, [GAMMA_DEPH] * N, [GAMMA_T1] * N)

        traj = evolve(L, rho_0, times)
        cpsi_t = np.array([cpsi(r) for r in traj])
        crossings = find_cusp_crossings(times, cpsi_t)

        print(f"{label}:")
        print(f"  Total CΨ crossings of 1/4: {len(crossings)}")
        if not crossings:
            print(f"  CΨ stays {'above' if cpsi_t[-1] > 0.25 else 'below'} 1/4 throughout.")
            print()
            continue

        # For each crossing, decompose ρ(t*) into L-eigenmode clusters
        # and report top contributions.
        for ci, (t_cross, direction, idx) in enumerate(crossings):
            rho_at_cross = traj[idx]  # close enough; could interpolate
            cluster_info = cluster_contributions_to_state(L, rho_at_cross, N)

            dir_str = "↑" if direction > 0 else "↓"
            print(f"  Crossing #{ci + 1}: t* = {t_cross:.4f} ({dir_str})")
            print(f"    Top 5 cluster contributions to ρ(t*) in Pauli/L-eigenbasis:")
            print(f"    {'cluster λ':>22s}  {'size':>4s}  {'‖c|cluster‖':>11s}  {'max |c_k|':>10s}")
            for c in cluster_info[:5]:
                ev = c['eigenvalue']
                ev_str = f"{ev.real:>+8.3f}{ev.imag:+8.3f}i"
                print(f"    {ev_str:>22s}  {c['size']:>4d}  "
                      f"{c['norm']:>11.4f}  {c['max_abs']:>10.4f}")
        print()

    # Cross-comparison: at each Hamiltonian's *first* cusp crossing,
    # which cluster eigenvalue dominates? Is there a shared structure?
    print()
    print("Cross-comparison: dominant L-eigenmode at first cusp crossing")
    print()
    print(f"  {'H':<14s}  {'t* first':>9s}  {'dom λ':>22s}  {'cluster size':>12s}  {'norm':>8s}")
    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        L = fw.lindbladian_z_plus_t1(H, [GAMMA_DEPH] * N, [GAMMA_T1] * N)
        traj = evolve(L, rho_0, times)
        cpsi_t = np.array([cpsi(r) for r in traj])
        crossings = find_cusp_crossings(times, cpsi_t)
        if not crossings:
            print(f"  {label:<14s}  {'—':>9s}  {'never crosses':>22s}  {'—':>12s}  {'—':>8s}")
            continue
        t_cross, direction, idx = crossings[0]
        rho_at = traj[idx]
        cluster_info = cluster_contributions_to_state(L, rho_at, N)
        top = cluster_info[0]
        ev = top['eigenvalue']
        ev_str = f"{ev.real:>+8.3f}{ev.imag:+8.3f}i"
        print(f"  {label:<14s}  {t_cross:>9.4f}  {ev_str:>22s}  "
              f"{top['size']:>12d}  {top['norm']:>8.4f}")

    print()
    print("Reading guide:")
    print("  If the dominant cluster at t* differs case-by-case in a structured")
    print("  way (e.g., always the same eigenvalue band, always same size),")
    print("  that's the operator signature.")
    print("  If the dominant cluster is unrelated to whether the case is robust")
    print("  or fragile, EQ-028 closes with a null result.")


if __name__ == "__main__":
    main()
