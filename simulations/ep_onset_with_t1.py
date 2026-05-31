#!/usr/bin/env python3
"""Would the EP onset (the sloshing memory) survive real hardware T1? A simulate-first check.

The post-EP signature we want to SEE on IBM is dynamical and cheap (no tomography): scan
Q = J/gamma_phi across Q_EP (~1.5) and watch a single excitation. Below: monotone decay
(forgetting). Above: it sloshes site to site (the reborn oscillating memory). The onset of
sloshing = the EP.

But our clean model has only Z-dephasing (number-conserving). Real hardware adds T1
(amplitude damping, sigma^-), which leaks the excitation to the ground state and drags the
true steady state to vacuum, not to 1/N. So before spending any QPU minute we ask, entirely
in simulation: does the sloshing survive a Kingston-like T1?

Model (the same channels Aer's thermal_relaxation_error puts on hardware):
    L = -i[H, .] + gamma_phi * sum_l D[Z_l] + gamma_1 * sum_l D[sigma^-_l]
H = J * sum_bonds (XX + YY) (single-excitation hopping). Units: gamma_phi = 1, so J = Q, and
gamma_1 = rho_T1 * gamma_phi. Kingston has T1 ~ T2, which puts rho_T1 ~ 1..2.

Read out per-site occupation <n_l>(t). Signatures:
    transfer  = max_t <n_{N-1}>(t)         how much reaches the far end (coherent sloshing)
    recurrence= does <n_0> rise again after its first dip (the memory comes back)
    survival  = total excitation sum_l <n_l> at the moment of best transfer (1 - what T1 ate)
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
SM = np.array([[0, 1], [0, 0]], complex)        # sigma^- : |1> -> |0>  (de-excite)
NUM = np.array([[0, 0], [0, 1]], complex)        # |1><1| = sigma^+ sigma^-


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def bond_op(N, b, P, Q):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == b else (Q if i == b + 1 else I2))
    return o


def liouvillian(N, J, g_phi, g_1):
    """vec column-stack: vec(A rho B) = (B^T (x) A) vec(rho)."""
    d = 2 ** N
    Id = np.eye(d)
    H = J * sum(bond_op(N, b, X, X) + bond_op(N, b, Y, Y) for b in range(N - 1))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))          # -i[H, .]
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += g_phi * (np.kron(Zl, Zl) - np.kron(Id, Id))   # D[Z_l] dephasing
        c = op_at(N, l, SM)
        cdc = op_at(N, l, NUM)                              # c^dagger c
        L += g_1 * (np.kron(c.conj(), c)
                    - 0.5 * np.kron(Id, cdc)
                    - 0.5 * np.kron(cdc.T, Id))             # D[sigma^-_l] amplitude damping
    return L


def occupations(N, J, g_phi, g_1, ts):
    """<n_l>(t) for a single excitation prepared at site 0."""
    d = 2 ** N
    L = liouvillian(N, J, g_phi, g_1)
    psi = np.array([1], complex)
    for i in range(N):
        psi = np.kron(psi, np.array([0, 1], complex) if i == 0 else np.array([1, 0], complex))
    rho0 = np.outer(psi, psi.conj()).flatten("F")
    w, V = np.linalg.eig(L)
    Vinv = np.linalg.inv(V)
    n_ops = [op_at(N, l, NUM) for l in range(N)]
    out = np.zeros((len(ts), N))
    for k, t in enumerate(ts):
        rho = (V @ (np.exp(w * t) * (Vinv @ rho0))).reshape(d, d, order="F")
        for l in range(N):
            out[k, l] = float(np.real(np.trace(rho @ n_ops[l])))
    return out


def signatures(N, occ, ts):
    n0, nfar = occ[:, 0], occ[:, N - 1]
    total = occ.sum(axis=1)
    transfer = float(nfar.max())
    k_transfer = int(np.argmax(nfar))
    survival = float(total[k_transfer])
    # recurrence: does n0 climb back up after its first local minimum?
    kmin = int(np.argmin(n0[: max(2, len(ts) // 2)]))
    recurrence = float(n0[kmin:].max() - n0[kmin]) if kmin < len(ts) - 1 else 0.0
    return transfer, survival, recurrence, ts[k_transfer]


def main():
    N = 3
    g_phi = 1.0
    Q_EP = 1.5
    ts = np.linspace(0.0, 16.0, 400)          # units of 1/gamma_phi
    print("=" * 86)
    print(f"  EP ONSET vs HARDWARE T1   (N={N} chain, single excitation, gamma_phi=1, Q_EP~{Q_EP})")
    print("=" * 86)
    print("  transfer = max <n_far>;  recurrence = <n_0> climb-back (sloshing);  survival = total excitation there")
    for rho_T1 in [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]:
        tag = ("ideal (our model, no T1)" if rho_T1 == 0
               else f"gamma_1/gamma_phi = {rho_T1:.2f}  "
                    + ("(injected dephasing, fast gates: J*T1>>1)" if rho_T1 <= 0.5
                       else "(idle-T2 regime: T1 ~ dephasing time)"))
        print(f"\n  --- T1 setting: {tag} ---")
        print(f"  {'Q=J/g':>7} {'regime':>14} {'transfer':>9} {'recurrence':>11} {'survival':>9}")
        for Q in [0.5, 1.0, 1.5, 2.5, 4.0]:
            occ = occupations(N, J=Q * g_phi, g_phi=g_phi, g_1=rho_T1 * g_phi, ts=ts)
            transfer, survival, recurrence, t_tr = signatures(N, occ, ts)
            regime = "overdamped" if Q < Q_EP else ("AT EP" if abs(Q - Q_EP) < 1e-9 else "sloshing")
            print(f"  {Q:7.1f} {regime:>14} {transfer:9.3f} {recurrence:11.3f} {survival:9.3f}")
    print("\n  Read: 'transfer' and 'recurrence' large = the sloshing memory is visible.")
    print("  With T1 (rho>0), 'survival' < 1 shows how much excitation T1 ate before transfer;")
    print("  if transfer/recurrence collapse from the ideal column, the EP onset is masked on hardware.")


if __name__ == "__main__":
    main()
