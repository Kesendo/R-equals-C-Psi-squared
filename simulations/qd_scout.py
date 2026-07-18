"""qd_scout.py — Quantum Darwinism in our own chain (the pointer door).

The pointer-reading door (project_no_wave_dies): the refcount as a computable.
One system qubit S, environment E = the other N-1 chain qubits, fragments F ⊆ E.
Metric discipline (outbound_label_adapters scar): genuine von Neumann mutual
information I(S:F) = S(rho_S) + S(rho_F) - S(rho_SF) on reduced states,
never pairwise-sum proxies.

Scenarios:
  A  star,  J=0          — Zurek's exactly solvable pure-decoherence model
                           (pipeline anchor: numeric rho_{S,E_k} vs closed form)
  B  chain, J=0          — no-broadcast: I(S:E_j) = 0 exactly for non-neighbors,
                           redundancy capped at 2 (the paper theorem, checked)
  C  chain, J>0          — the record current: redundancy only via transport
                           (C1 pointer regime J=0.25, C2 scrambling J=1.0)
  D  chain, J>0, gamma>0 — the garbage collector: local Z-dephasing erases
                           coherence records (canonical gamma_0 = 0.05)

Conventions: site 0 = MSB (Propagate/Core big-endian).
H = J * sum_bonds (X_i X_j + Y_i Y_j) + Delta * sum_bonds Z_i Z_j   (Pauli J).
Dephasing: D[rho] = sum_l gamma_l (Z_l rho Z_l - rho); entrywise mask
-2 * sum_l gamma_l * [bit_l(i) != bit_l(j)].

Redundancy convention: H_ref = 1 bit (the pointer entropy of Z_S; pinned exactly
by the X^N bridge: [H, X^N] = 0 and |+>^N is X^N-invariant, so <Z_j> = 0 for all t).
m_delta(t) = min m with mean_{|F|=m} I(S:F) >= (1-delta) * H_ref; R_delta = nE / m_delta.
"""

import numpy as np
from itertools import combinations
import time as walltime

np.set_printoptions(linewidth=200, precision=6, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

DELTA_QD = 0.10   # the QD deficit delta
N = 8             # 1 system + 7 environment qubits
DT = 0.005


def op_at(n, site, P):
    M = np.array([[1.0]], dtype=complex)
    for s in range(n):
        M = np.kron(M, P if s == site else I2)
    return M


def two_site_zz(n, i, j):
    M = np.array([[1.0]], dtype=complex)
    for s in range(n):
        M = np.kron(M, Z if s in (i, j) else I2)
    return M


def build_H(n, bonds, J, Delta):
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        if J != 0.0:
            H += J * (op_at(n, i, X) @ op_at(n, j, X) + op_at(n, i, Y) @ op_at(n, j, Y))
        if Delta != 0.0:
            H += Delta * two_site_zz(n, i, j)
    return H


def deph_mask(n, gammas):
    d = 2 ** n
    idx = np.arange(d)
    mask = np.zeros((d, d))
    for l in range(n):
        bit = (idx >> (n - 1 - l)) & 1          # site l is MSB-first
        diff = (bit[:, None] != bit[None, :])
        mask -= 2.0 * gammas[l] * diff
    return mask


def rk4_evolve(rho, H, mask, dt, steps):
    def rhs(r):
        return -1j * (H @ r - r @ H) + mask * r
    for _ in range(steps):
        k1 = rhs(rho)
        k2 = rhs(rho + 0.5 * dt * k1)
        k3 = rhs(rho + 0.5 * dt * k2)
        k4 = rhs(rho + dt * k3)
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return rho


def ptrace(rho, n, keep):
    keep = sorted(keep)
    t = rho.reshape((2,) * (2 * n))
    row = list(range(n))
    col = [(n + i if i in keep else i) for i in range(n)]
    out = keep + [n + i for i in keep]
    red = np.einsum(t, row + col, out)
    k = len(keep)
    return red.reshape((2 ** k, 2 ** k))


def vn_entropy(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-14]
    return float(-(w * np.log2(w)).sum())


def mutual_info(rho, n, A, B):
    SA = vn_entropy(ptrace(rho, n, A))
    SB = vn_entropy(ptrace(rho, n, B))
    SAB = vn_entropy(ptrace(rho, n, sorted(A + B)))
    return SA + SB - SAB


def sanity(rho, tag):
    tr = np.trace(rho).real
    herm = np.abs(rho - rho.conj().T).max()
    emin = np.linalg.eigvalsh(rho).min()
    ok = abs(tr - 1) < 1e-8 and herm < 1e-10 and emin > -1e-6   # eig_min tolerance = RK4 drift
    print(f"    [sanity {tag}] trace-1={tr-1:+.2e}  herm={herm:.2e}  eig_min={emin:+.2e}  {'OK' if ok else 'FAIL'}")
    return ok


def star_analytic_pair(t, Delta, nE):
    """Closed-form rho_{S,E_k} for the star at J=0, gamma=0, |+>^N start."""
    c = np.cos(2 * Delta * t)
    a = np.array([np.exp(-1j * Delta * t), np.exp(1j * Delta * t)]) / np.sqrt(2)
    b = np.array([np.exp(1j * Delta * t), np.exp(-1j * Delta * t)]) / np.sqrt(2)
    Paa = np.outer(a, a.conj())
    Pbb = np.outer(b, b.conj())
    Pab = np.outer(a, b.conj())
    r = c ** (nE - 1)
    top = np.hstack([0.5 * Paa, 0.5 * r * Pab])
    bot = np.hstack([0.5 * r * Pab.conj().T, 0.5 * Pbb])
    return np.vstack([top, bot])


def fragment_sweep(rho, n, S, env, t, label):
    """Full partial-information curve + redundancy at one time."""
    nE = len(env)
    HS = vn_entropy(ptrace(rho, n, [S]))
    pS = np.diag(ptrace(rho, n, [S])).real
    Hpointer = float(-(pS[pS > 1e-14] * np.log2(pS[pS > 1e-14])).sum())
    by_size = {}
    for m in range(1, nE + 1):
        vals = [mutual_info(rho, n, [S], list(F)) for F in combinations(env, m)]
        by_size[m] = (float(np.mean(vals)), float(np.min(vals)), float(np.max(vals)))
    thr = (1 - DELTA_QD) * 1.0   # H_ref = 1 bit pointer entropy (X^N-pinned)
    m_delta = next((m for m in range(1, nE + 1) if by_size[m][0] >= thr), None)
    R = (nE / m_delta) if m_delta else 0.0
    print(f"  t={t:6.3f} [{label}]  S(rho_S)={HS:.4f}  H_pointer={Hpointer:.6f}")
    curve = "    Ibar(m): " + "  ".join(f"m={m}:{by_size[m][0]:.4f}" for m in range(1, nE + 1))
    print(curve)
    print(f"    m_delta={m_delta}  R_delta={R:.3f}   (delta={DELTA_QD}, threshold={thr:.2f} bit)")
    return by_size, m_delta, R, HS


def singleton_profile(rho, n, S, env):
    return [mutual_info(rho, n, [S], [j]) for j in env]


def run_scenario(name, bonds, J, Delta, gammas, S, t_singles, t_sweeps, star_gate=False):
    print()
    print("=" * 100)
    print(f"SCENARIO {name}   J={J}  Delta={Delta}  gammas={'0' if not any(gammas) else gammas[0]}  S=site {S}")
    print("=" * 100)
    env = [j for j in range(N) if j != S]
    H = build_H(N, bonds, J, Delta)
    mask = deph_mask(N, gammas)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = np.array([1.0], dtype=complex)
    for _ in range(N):
        psi = np.kron(psi, plus)
    rho = np.outer(psi, psi.conj())

    checkpoints = sorted(set(list(t_singles) + list(t_sweeps)))
    t_now = 0.0
    results = {}
    for t_target in checkpoints:
        steps = int(round((t_target - t_now) / DT))
        rho = rk4_evolve(rho, H, mask, DT, steps)
        t_now += steps * DT
        # re-hermitize against RK4 drift
        rho = 0.5 * (rho + rho.conj().T)
        if t_target in t_singles:
            prof = singleton_profile(rho, N, S, env)
            print(f"  t={t_now:6.3f}  I(S:E_j) singletons: "
                  + "  ".join(f"j={j}:{v:.2e}" for j, v in zip(env, prof)))
        if star_gate:
            nE = len(env)
            worst = 0.0
            for k in env:
                num = ptrace(rho, N, [S, k])
                ana = star_analytic_pair(t_now, Delta, nE)
                worst = max(worst, np.abs(num - ana).max())
            print(f"    [GATE star-analytic] max|numeric - closed-form| over all E_k = {worst:.2e}"
                  f"  {'OK' if worst < 1e-7 else 'FAIL'}")
        if t_target in t_sweeps:
            sanity(rho, f"t={t_now:.2f}")
            results[t_target] = fragment_sweep(rho, N, S, env, t_now, name)
    return results, rho


def main():
    t0 = walltime.time()
    print(f"QD scout  N={N}  (1 system + {N-1} environment)  dt={DT}  delta={DELTA_QD}")
    print("Metric: genuine von Neumann I(S:F) on reduced states (no proxies).")

    # --- A: star, J=0 — Zurek anchor -------------------------------------
    S_star = 0
    star_bonds = [(0, k) for k in range(1, N)]
    run_scenario("A star J=0 (Zurek)", star_bonds, J=0.0, Delta=1.0,
                 gammas=[0.0] * N, S=S_star,
                 t_singles=[0.25, np.pi / 4, 1.5],
                 t_sweeps=[0.25, np.pi / 4, 1.5],
                 star_gate=True)
    print("  note: J=0 star is periodic (period pi/2 in the record overlap cos 2*Delta*t);")
    print("        perfect records at t = pi/4; recurrence (records erased) at t = pi/2.")

    # --- B: chain, J=0 — no-broadcast horizon ----------------------------
    S_chain = 3
    chain_bonds = [(i, i + 1) for i in range(N - 1)]
    _, rhoB = run_scenario("B chain J=0 (no-broadcast)", chain_bonds, J=0.0, Delta=1.0,
                           gammas=[0.0] * N, S=S_chain,
                           t_singles=[np.pi / 4, 2.0, 5.0],
                           t_sweeps=[np.pi / 4, 5.0])
    # Shared-dresser law at J=0: I(S:E_j) = 0 exactly iff the ZZ-neighborhoods
    # N(S)={2,4} and N(j) share no qubit, i.e. graph distance(S,j) >= 3.
    env = [j for j in range(N) if j != S_chain]
    prof = {j: mutual_info(rhoB, N, [S_chain], [j]) for j in env}
    far = {0, 6, 7}          # distance >= 3 from S=3
    near = {1, 2, 4, 5}      # distance <= 2 (neighbor or shared dresser)
    far_max = max(abs(prof[j]) for j in far)
    near_min = min(prof[j] for j in near)
    print(f"  [GATE shared-dresser] max|I| at distance>=3: {far_max:.2e} (must be ~0);"
          f"  min I at distance<=2: {near_min:.4f} (must be >0)"
          f"  {'OK' if far_max < 1e-9 and near_min > 1e-3 else 'FAIL'}")

    # --- C1: chain, pointer regime J=0.25 --------------------------------
    run_scenario("C1 chain J=0.25 (record current, pointer regime)", chain_bonds,
                 J=0.25, Delta=1.0, gammas=[0.0] * N, S=S_chain,
                 t_singles=[0.5, 1.0, 2.0, 4.0, 8.0],
                 t_sweeps=[1.0, 2.0, 4.0, 8.0])

    # --- C2: chain, hopping-dominant J=1.0, Delta=0.25 --------------------
    # NOT J=Delta: at the Heisenberg point |+>^N is an exact fixed point
    # (the outbound_label_adapters scar) and nothing evolves.
    run_scenario("C2 chain J=1.0 Delta=0.25 (hopping-dominant)", chain_bonds,
                 J=1.0, Delta=0.25, gammas=[0.0] * N, S=S_chain,
                 t_singles=[0.5, 1.0, 2.0, 4.0],
                 t_sweeps=[1.0, 2.0, 4.0])

    # --- D: the two gamma knobs, separated --------------------------------
    gS = [0.0] * N; gS[S_chain] = 0.05
    run_scenario("D1 chain J=0.25, gamma on S only (S watched from outside)", chain_bonds,
                 J=0.25, Delta=1.0, gammas=gS, S=S_chain,
                 t_singles=[0.5, 1.0, 2.0, 4.0, 8.0],
                 t_sweeps=[1.0, 4.0, 8.0])

    gE = [0.05] * N; gE[S_chain] = 0.0
    run_scenario("D2 chain J=0.25, gamma on E only (records watched away = GC)", chain_bonds,
                 J=0.25, Delta=1.0, gammas=gE, S=S_chain,
                 t_singles=[0.5, 1.0, 2.0, 4.0, 8.0],
                 t_sweeps=[1.0, 4.0, 8.0])

    run_scenario("D3 chain J=0.25, uniform gamma=0.05 (the canonical world)", chain_bonds,
                 J=0.25, Delta=1.0, gammas=[0.05] * N, S=S_chain,
                 t_singles=[1.0, 4.0, 8.0],
                 t_sweeps=[1.0, 4.0, 8.0])

    # --- E: binary tree — the pointer-pointer geometry --------------------
    # S = root (site 0), children 2i+1, 2i+2. Fan-out between chain and star.
    tree_bonds = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7)]
    _, rhoE = run_scenario("E1 tree J=0 (record horizon in tree geometry)", tree_bonds,
                           J=0.0, Delta=1.0, gammas=[0.0] * N, S=0,
                           t_singles=[np.pi / 4, 2.0, 5.0],
                           t_sweeps=[5.0])
    # Shared-dresser law in the tree: N(S)={1,2}; nonzero iff j in N(S) or
    # N(j) ∩ N(S) != {}: {1,2} direct, {3,4} share dresser 1, {5,6} share
    # dresser 2; site 7 (N(7)={3}) is exactly blind.
    envE = list(range(1, N))
    profE = {j: mutual_info(rhoE, N, [0], [j]) for j in envE}
    farE_max = abs(profE[7])
    nearE_min = min(profE[j] for j in [1, 2, 3, 4, 5, 6])
    print(f"  [GATE shared-dresser tree] |I(S:E_7)|={farE_max:.2e} (must be ~0);"
          f"  min I over sites 1-6: {nearE_min:.4f} (must be >0)"
          f"  {'OK' if farE_max < 1e-9 and nearE_min > 1e-3 else 'FAIL'}")

    run_scenario("E2 tree J=0.25 (record current on the tree)", tree_bonds,
                 J=0.25, Delta=1.0, gammas=[0.0] * N, S=0,
                 t_singles=[0.5, 1.0, 2.0, 4.0, 8.0],
                 t_sweeps=[1.0, 4.0, 8.0])

    print()
    print(f"total wall time: {walltime.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
