"""
Subsystem Crossing reproduction: the four pair-level tables of
experiments/SUBSYSTEM_CROSSING.md, rebuilt from scratch.

Setup (as in the experiment): N=4 Heisenberg RING (J=1, Pauli operators,
h=0), local Z-dephasing gamma=0.05 per site, Lindblad
L(rho) = -i[H,rho] + gamma * sum_l (Z_l rho Z_l - rho).
For each qubit pair (i,j): partial trace, l1-coherence, Psi = l1/3,
Wootters concurrence, the connected-correlator bridge
C_corr = (|<XX>_c| + |<YY>_c| + |<ZZ>_c|)/3, and the crossing of
CPsi = concurrence * Psi through 1/4.

Also prints the isolated Bell+ (N=2) baseline and the three C-readings
of its one dephasing trajectory f = e^(-4*gamma*t):
  concurrence book  C = f          -> CPsi = f^2/3,        crosses at t = 0.719
  F25 purity book   C = (1+f^2)/2  -> CPsi = f(1+f^2)/6,   crosses at t = 0.747
  constant bridge   C = 1          -> CPsi = f/3,          crosses at t = 1.438
(one trajectory, three books; the experiment's tables are the
concurrence book).

Original tables were produced by the retired delta_calc MCP tool
(February 2026); this script is the committed reproduction.
"""

import numpy as np
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, site, N):
    ops = [I2] * N
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def liouvillian(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site_op(Z, l, N)
        L += gamma * (np.kron(Zl, Zl.conj()) - np.eye(d * d, dtype=complex))
    return L


def heisenberg_ring(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    bonds = [(i, (i + 1) % N) for i in range(N)]
    for (i, j) in bonds:
        for P in (X, Y, Z):
            H += J * site_op(P, i, N) @ site_op(P, j, N)
    return H


def ptrace_pair(rho, keep, N):
    dims = [2] * N
    rho_t = rho.reshape(dims + dims)
    out = sorted(set(range(N)) - set(keep))
    for q in sorted(out, reverse=True):
        rho_t = np.trace(rho_t, axis1=q, axis2=q + rho_t.ndim // 2)
    d = 2 ** len(keep)
    return rho_t.reshape(d, d)


def l1_coherence(rho):
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))


def concurrence(rho):
    yy = np.kron(Y, Y)
    R = rho @ yy @ rho.conj() @ yy
    ev = np.sqrt(np.abs(np.sort(np.linalg.eigvals(R).real)[::-1]))
    return max(0.0, ev[0] - ev[1] - ev[2] - ev[3])


def bridge(rho):
    c = 0.0
    for P in (X, Y, Z):
        PP = np.kron(P, P)
        exp_pp = np.trace(rho @ PP).real
        exp_a = np.trace(rho @ np.kron(P, I2)).real
        exp_b = np.trace(rho @ np.kron(I2, P)).real
        c += abs(exp_pp - exp_a * exp_b)
    return c / 3.0


def crossing_time(ts, vals, level=0.25):
    for k in range(1, len(vals)):
        if vals[k - 1] > level >= vals[k]:
            frac = (vals[k - 1] - level) / (vals[k - 1] - vals[k])
            return ts[k - 1] + frac * (ts[k] - ts[k - 1])
    return None


def run_state(name, psi, N=4, gamma=0.05, t_max=5.0, dt=0.01):
    print(f"\n=== {name} (N={N}, Heisenberg ring J=1, gamma={gamma}) ===")
    rho0 = np.outer(psi, psi.conj())
    H = heisenberg_ring(N)
    L = liouvillian(H, gamma, N)
    step = expm(L * dt)
    v = rho0.reshape(-1)
    ts = np.arange(0, t_max + dt / 2, dt)
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    traj = {p: [] for p in pairs}
    for _ in ts:
        rho = v.reshape(2 ** N, 2 ** N)
        for p in pairs:
            rp = ptrace_pair(rho, list(p), N)
            traj[p].append(concurrence(rp) * l1_coherence(rp) / 3.0)
        v = step @ v
    rho0_full = rho0
    print(f"{'pair':>6} {'l1(0)':>7} {'Psi(0)':>7} {'conc(0)':>8} "
          f"{'C_corr(0)':>10} {'CPsi(0)':>8} {'t_cross':>8}")
    for p in pairs:
        rp0 = ptrace_pair(rho0_full, list(p), N)
        l1 = l1_coherence(rp0)
        tc = crossing_time(ts, traj[p])
        print(f"{str(p):>6} {l1:7.3f} {l1/3:7.3f} {concurrence(rp0):8.3f} "
              f"{bridge(rp0):10.3f} {concurrence(rp0)*l1/3:8.3f} "
              f"{tc if tc is not None else float('nan'):8.3f}")


def main():
    up = np.array([1, 0], complex)
    dn = np.array([0, 1], complex)

    def kron_all(*vs):
        r = vs[0]
        for w in vs[1:]:
            r = np.kron(r, w)
        return r

    bell = (kron_all(up, up) + kron_all(dn, dn)) / np.sqrt(2)
    ghz = (kron_all(up, up, up, up) + kron_all(dn, dn, dn, dn)) / np.sqrt(2)
    w = (kron_all(dn, up, up, up) + kron_all(up, dn, up, up)
         + kron_all(up, up, dn, up) + kron_all(up, up, up, dn)) / 2
    plus = (up + dn) / np.sqrt(2)

    run_state("GHZ", ghz)
    run_state("W", w)
    run_state("Bell+ x Bell+", np.kron(bell, bell))
    run_state("|+>^4", kron_all(plus, plus, plus, plus))

    # Isolated Bell+ baseline: one dephasing trajectory, three C-books.
    gamma = 0.05
    print("\n=== Isolated Bell+ (N=2): one trajectory f = e^(-4*gamma*t), "
          "three C-books ===")
    for label, cpsi_of_f in [
            ("concurrence book  CPsi = f^2/3", lambda f: f * f / 3),
            ("F25 purity book   CPsi = f(1+f^2)/6", lambda f: f * (1 + f * f) / 6),
            ("constant bridge   CPsi = f/3", lambda f: f / 3)]:
        ts = np.arange(0, 3.0, 0.0001)
        f = np.exp(-4 * gamma * ts)
        tc = crossing_time(ts, cpsi_of_f(f))
        print(f"  {label:38s} crosses 1/4 at t = {tc:.4f} "
              f"(K = gamma*t = {gamma*tc:.5f})")
    print("\nThe experiment's tables are the concurrence book. The 'nine "
          "times faster' ratio: 0.7192 / 0.0797 = 9.0.")


if __name__ == "__main__":
    main()
