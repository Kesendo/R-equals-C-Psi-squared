"""
Finite N=4 scalar-readout catalogue: the four pair-level tables of
experiments/SUBSYSTEM_CROSSING.md, rebuilt from scratch.  Every threshold
statement below is scoped to the selected preparation, topology, and time
grid.  It is not a local quantum/classical transition, measurement law, or
entanglement boundary.

Setup (as in the experiment): N=4 Heisenberg RING (J=1, Pauli operators,
h=0), local Z-dephasing gamma=0.05 per site, Lindblad
L(rho) = -i[H,rho] + gamma * sum_l (Z_l rho Z_l - rho).
For each qubit pair (i,j): partial trace, l1-coherence, Psi = l1/3,
Wootters concurrence, the connected-correlator bridge
C_corr = (|<XX>_c| + |<YY>_c| + |<ZZ>_c|)/3, and the crossing of
the selected pair readout R_pair = concurrence * Psi through the numerical
reference 1/4.

Also prints the isolated Bell+ (N=2) baseline and the three C-readings
of its one dephasing trajectory f = e^(-4*gamma*t):
  concurrence book  C = f          -> CPsi = f^2/3,        crosses at t = 0.719
  F25 purity book   C = (1+f^2)/2  -> CPsi = f(1+f^2)/6,   crosses at t = 0.747
  constant bridge   C = 1          -> CPsi = f/3,          crosses at t = 1.438
(one trajectory, three books; the experiment's tables are the
concurrence book).

Original tables were produced by the retired delta_calc MCP tool
(February 2026); this script is the committed reproduction.

Also reproduces (2026-07-20) the upward-crossing product-state tables of
experiments/ORPHANED_RESULTS.md section 2b/2c: |0+0+> on the N=4 ring
peaks at CPsi = 0.200 (no crossing) but crosses on the CHAIN (pair (1,2),
0.310); on the ring |+-+-> (0.284, ring neighbours) and |0+0-> (0.256,
diagonal) cross. The conflicting ring-(0,2) table in
experiments/DYNAMIC_ENTANGLEMENT.md came from the retired MCP tool
and does not reproduce under this convention (see the reproduction note
there).
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


def heisenberg_chain(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in (X, Y, Z):
            H += J * site_op(P, i, N) @ site_op(P, i + 1, N)
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


def bell_plus_dephased_pair(f):
    """Bell+ density matrix with its sole coherence scaled by 0 <= f <= 1."""
    return np.array([
        [0.5, 0.0, 0.0, f / 2],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [f / 2, 0.0, 0.0, 0.5],
    ], dtype=complex)


def selected_pair_readout(rho):
    """Named concurrence-coherence product used by this finite catalogue."""
    return concurrence(rho) * l1_coherence(rho) / 3.0


def isolated_bell_readout_books(f):
    """Three inequivalent scalar books evaluated on one Bell+ trajectory."""
    return {
        "concurrence": f * f / 3,
        "purity": f * (1 + f * f) / 6,
        "constant": f / 3,
    }


def bridge(rho):
    c = 0.0
    for P in (X, Y, Z):
        PP = np.kron(P, P)
        exp_pp = np.trace(rho @ PP).real
        exp_a = np.trace(rho @ np.kron(P, I2)).real
        exp_b = np.trace(rho @ np.kron(I2, P)).real
        c += abs(exp_pp - exp_a * exp_b)
    return c / 3.0


def crossing_time(ts, vals, level=0.25, direction="down"):
    if direction not in ("down", "up"):
        raise ValueError("direction must be 'down' or 'up'")
    for k in range(1, len(vals)):
        if direction == "down":
            crossed = vals[k - 1] > level >= vals[k]
        else:
            crossed = vals[k - 1] < level <= vals[k]
        if crossed:
            frac = (level - vals[k - 1]) / (vals[k] - vals[k - 1])
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
            value = selected_pair_readout(rp)
            traj[p].append(value)
        v = step @ v
    rho0_full = rho0
    print(f"{'pair':>6} {'l1(0)':>7} {'Psi(0)':>7} {'conc(0)':>8} "
          f"{'C_corr(0)':>10} {'R_pair(0)':>10} {'t_down_.25':>11}")
    rows = []
    for p in pairs:
        rp0 = ptrace_pair(rho0_full, list(p), N)
        l1 = l1_coherence(rp0)
        tc = crossing_time(ts, traj[p], level=0.25, direction="down")
        rows.append({
            "pair": p,
            "initial_l1": l1,
            "initial_psi": l1 / 3,
            "initial_concurrence": concurrence(rp0),
            "initial_bridge": bridge(rp0),
            "initial_readout": selected_pair_readout(rp0),
            "downward_reference_time": tc,
        })
        print(f"{str(p):>6} {l1:7.3f} {l1/3:7.3f} {concurrence(rp0):8.3f} "
              f"{bridge(rp0):10.3f} {selected_pair_readout(rp0):10.3f} "
              f"{tc if tc is not None else float('nan'):8.3f}")
    return rows


def run_upward(name, psi, topology="ring", N=4, gamma=0.05,
               t_max=3.0, dt=0.005):
    """Finite product-state table: per-pair maximum of the selected readout."""
    H = heisenberg_ring(N) if topology == "ring" else heisenberg_chain(N)
    print(f"\n=== {name} (N={N}, Heisenberg {topology} J=1, gamma={gamma}) ===")
    L = liouvillian(H, gamma, N)
    step = expm(L * dt)
    v = np.outer(psi, psi.conj()).reshape(-1)
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    best = {p: (0.0, 0.0) for p in pairs}
    t = 0.0
    while t <= t_max + dt / 2:
        rho = v.reshape(2 ** N, 2 ** N)
        for p in pairs:
            rp = ptrace_pair(rho, list(p), N)
            value = selected_pair_readout(rp)
            if value > best[p][0]:
                best[p] = (value, t)
        v = step @ v
        t += dt
    rows = []
    for p in pairs:
        c, tm = best[p]
        rows.append({"pair": p, "maximum_readout": c, "time": tm})
        print(f"  pair {p}: max selected readout = {c:.3f} at t = {tm:.3f}"
              f"{'   above selected 1/4 readout reference' if c > 0.25 else ''}")
    return rows


def main():
    print("FINITE N=4 SCALAR-READOUT CATALOGUE")
    print("Selected preparation, topology, and time grid; the selected 1/4 "
          "readout reference is not a transition.")
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
    bell_ring_rows = run_state("Bell+ x Bell+", np.kron(bell, bell))
    run_state("|+>^4", kron_all(plus, plus, plus, plus))

    # Upward crossings from product states (ORPHANED_RESULTS section 2b/2c;
    # settles the DYNAMIC_ENTANGLEMENT ring-(0,2) dispute).
    minus = (up - dn) / np.sqrt(2)
    run_upward("|0+0+>", kron_all(up, plus, up, plus), "ring")
    run_upward("|0+0+>", kron_all(up, plus, up, plus), "chain")
    run_upward("|+-+->", kron_all(plus, minus, plus, minus), "ring")
    run_upward("|0+0->", kron_all(up, plus, up, minus), "ring")

    # Isolated Bell+ baseline: one dephasing trajectory, three C-books.
    gamma = 0.05
    print("\n=== Isolated Bell+ (N=2): one trajectory f = e^(-4*gamma*t), "
          "three C-books ===")
    book_labels = {
        "concurrence": "concurrence book  R=f^2/3",
        "purity": "purity book       R=f(1+f^2)/6",
        "constant": "constant book     R=f/3",
    }
    ts = np.arange(0, 3.0, 0.0001)
    f = np.exp(-4 * gamma * ts)
    isolated_books = isolated_bell_readout_books(f)
    isolated_times = {}
    for book_name, label in book_labels.items():
        tc = crossing_time(
            ts, isolated_books[book_name], level=0.25, direction="down",
        )
        isolated_times[book_name] = tc
        print(f"  {label:38s} reaches the selected 1/4 reference at t = {tc:.4f} "
               f"(K = gamma*t = {gamma*tc:.5f})")
    counterexample = bell_plus_dephased_pair(0.5)
    print("\nEntanglement counterexample below the 1/4 readout reference: "
          f"concurrence={concurrence(counterexample):.3f}, "
          f"selected readout={selected_pair_readout(counterexample):.6f}.")
    ring_reference_time = next(
        row["downward_reference_time"] for row in bell_ring_rows
        if row["pair"] == (0, 1)
    )
    ratio = isolated_times["concurrence"] / ring_reference_time
    print("The experiment's tables use the concurrence book.")
    print("Ratio is computed from the returned ring row and isolated concurrence book: "
          f"{isolated_times['concurrence']:.4f} / {ring_reference_time:.4f} "
          f"= {ratio:.3f}. It is a finite comparison, not a universal rate ratio.")


if __name__ == "__main__":
    main()
