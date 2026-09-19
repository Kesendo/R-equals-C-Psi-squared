"""N=2 perfect-mirror null for a population-ratio diagnostic.

For |++> on the N=2 Heisenberg chain, both the unitary and Z-dephased
computational-basis populations remain 1/4.  Where the unitary population is
nonzero, writing C_i = R_i / P_i^u is merely a ratio: the ratio is definitional.
Here C_i = 1 up to floating roundoff throughout the finite gamma grid, so there
is no gamma inversion and no sensitivity claim.  The file is retained as a
negative control, not as evidence that this ratio is a physical observable.
"""
import sys
import numpy as np
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def op(P, l, N):
    m = np.array([[1.0 + 0j]])
    for k in range(N):
        m = np.kron(m, P if k == l else I2)
    return m


def build_setup():
    N = 2
    J = 1.0
    bonds = [(0, 1)]
    H = sum(
        J * (
            op(X, i, N) @ op(X, j, N)
            + op(Y, i, N) @ op(Y, j, N)
            + op(Z, i, N) @ op(Z, j, N)
        )
        for (i, j) in bonds
    )
    d = 2 ** N
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = plus
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    rho0 = np.outer(psi, psi.conj())
    identity = np.eye(d, dtype=complex)
    labels = [format(i, f"0{N}b") for i in range(d)]
    return N, J, H, d, rho0, identity, labels


def liouvillian(g, H, N, d, identity):
    Lh = -1j * (np.kron(identity, H) - np.kron(H.T, identity))
    Ld = np.zeros((d * d, d * d), dtype=complex)
    for l in range(N):
        Zl = op(Z, l, N)
        Ld += g * (np.kron(Zl.T, Zl) - np.kron(identity, identity))
    return Lh + Ld


def R_open(g, t, H, N, d, rho0, identity):
    rho = (
        expm(liouvillian(g, H, N, d, identity) * t)
        @ rho0.reshape(-1, order="F")
    ).reshape(d, d, order="F")
    return np.real(np.diag(rho))


def Psi2_closed(t, H, rho0):
    U = expm(-1j * H * t)
    rho = U @ rho0 @ U.conj().T
    return np.real(np.diag(rho))


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    N, J, H, d, rho0, identity, labels = build_setup()
    t = 0.7
    Psi2 = Psi2_closed(t, H, rho0)
    dom = int(np.argmax(Psi2))
    print(f"=== N=2 perfect-mirror null  (Heisenberg, |++>, t={t}, J={J}) ===\n", flush=True)

    # (1) A finite population-ratio table.  The ratio is definitional.
    g0 = 0.3
    R0 = R_open(g0, t, H, N, d, rho0, identity)
    print(f"[1] At gamma={g0} (Q=J/gamma={J/g0:.2f}, K=gamma*t={g0*t:.3f})", flush=True)
    print("    C_i = R_i / P_i^u is a definitional ratio, not a fitted coupling.", flush=True)
    print(f"    {'outcome':>8} {'P^u':>12} {'R (open)':>12} {'R/P^u':>10}", flush=True)
    for i in range(d):
        if Psi2[i] > 1e-9:
            print(f"    {labels[i]:>8} {Psi2[i]:>12.5f} {R0[i]:>12.5f} {R0[i]/Psi2[i]:>10.5f}", flush=True)

    # (2) Negative control: the finite gamma grid remains at unity.
    print(f"\n[2] selected outcome |{labels[dom]}>: finite null control", flush=True)
    print(f"    {'gamma':>8} {'Q=J/g':>8} {'K=g*t':>8} {'R/P^u':>10} {'ratio-1':>12}", flush=True)
    gs = [0.02, 0.05, 0.1, 0.2, 0.4, 0.8]
    devs = []
    for g in gs:
        C = R_open(g, t, H, N, d, rho0, identity)[dom] / Psi2[dom]
        devs.append(C - 1.0)
        print(f"    {g:>8.3g} {J/g:>8.2f} {g*t:>8.3f} {C:>10.5f} {C-1.0:>12.3e}", flush=True)
    max_abs_deviation = max(abs(value) for value in devs)
    print(f"    max |ratio-1| on this grid: {max_abs_deviation:.3e}", flush=True)
    print("\n[3] Verdict: no gamma inversion in this N=2 null.", flush=True)
    print("    The stationary populations make the ratio insensitive; roundoff-scale", flush=True)
    print("    deviations do not define a monotone calibration curve.", flush=True)


if __name__ == "__main__":
    main()
