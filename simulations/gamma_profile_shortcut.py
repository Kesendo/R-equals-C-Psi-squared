#!/usr/bin/env python3
"""Why can some profiles take a shortcut, reaching the closed-form rate without travelling the
loop (the Q-dressing)?

Hypothesis: a profile's rate is Q-independent (the shortcut) exactly when its slowest mode is a
SHARED eigenvector of the dephasing-Liouvillian and the Hamiltonian-Liouvillian. Then the decay
(Takt) is pure dephasing (Q-independent) and the rotation (Rotation) is pure Hamiltonian (Q*omega),
the two clock hands DECOUPLE, and the mode does not change with Q: it sits at the closed form
directly. A loop profile's slowest mode is H-dressed, so it CHANGES with Q (it must travel).

Test, per profile:
  - slowest-mode overlap |<v(Q=1.5) | v(Q=1000)>|  (~1 = frozen mode = shortcut; <1 = travels)
  - commutator residual ||[H, rho_slow]|| relative to ||H rho_slow||  (~0 = commutes = shortcut)
  - the rate at low vs high Q (equal = shortcut; converging = loop)
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def H_xy_unit(N):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        for P in (X, Y):
            t = np.array([[1]], complex)
            for i in range(N):
                t = np.kron(t, P if i in (b, b + 1) else I2)
            H += t
    return H


def L_dimless(N, Q, profile, H1):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def slowest(N, Q, profile, H1):
    L = L_dimless(N, Q, profile, H1)
    w, V = np.linalg.eig(L)
    nonkernel = np.abs(w) > 1e-7
    re = np.where(nonkernel, w.real, -np.inf)
    k = int(np.argmax(re))
    return -float(w[k].real), V[:, k]


def main():
    N = 5
    d = 2 ** N
    H1 = H_xy_unit(N)
    profiles = {
        "uniform":           [1, 1, 1, 1, 1],
        "V peaked (edges.25)": [0.25, 0.75, 3.0, 0.75, 0.25],
        "V peaked strong":   [0.1, 0.5, 3.8, 0.5, 0.1],
        "edge flat-bulk":    [0.25, 1.5, 1.5, 1.5, 0.25],
        "center-quiet":      list(np.array([1.333, 1.333, 0.5, 1.333, 0.5]) * 5 / 4.999),
        "inverse-V":         list(np.array([2.4, 0.6, 0.2, 0.6, 2.4]) * 5 / 6.2),
    }
    print(f"N={N}. Why the shortcut? frozen mode (overlap~1) + commuting ([H,rho]~0) = shortcut\n")
    print(f"  {'profile':>20}  {'rate Q=1.5':>10}  {'rate Q=1000':>11}  {'mode overlap':>12}  {'||[H,rho]||':>11}")
    for name, p in profiles.items():
        p = list(np.array(p, float) * N / np.sum(p))
        r_lo, v_lo = slowest(N, 1.5, p, H1)
        r_hi, v_hi = slowest(N, 1000.0, p, H1)
        ov = max(abs(np.vdot(v_lo, v_hi)), abs(np.vdot(v_lo, v_hi.conj()))) / (np.linalg.norm(v_lo) * np.linalg.norm(v_hi))
        # commutator residual of the slow mode with H (unvec the high-Q mode)
        rho = v_hi.reshape(d, d, order="F")
        comm = H1 @ rho - rho @ H1
        res = np.linalg.norm(comm) / (np.linalg.norm(H1 @ rho) + 1e-15)
        print(f"  {name:>20}  {r_lo:>10.5f}  {r_hi:>11.5f}  {ov:>12.6f}  {res:>11.5f}")
    print("\n  reading: shortcut profiles (rate Q-independent) should show mode overlap ~1 (the slowest")
    print("  mode does not change with Q) and a small commutator (the mode commutes with H, so the")
    print("  Hamiltonian term -iQ[H,.] does not shift its rate). Loop profiles: overlap < 1, larger comm.")


if __name__ == "__main__":
    main()
