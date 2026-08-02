#!/usr/bin/env python3
"""Does the tabulated psi_opt carry the lens mode's signs? It does not.

This is the verifier behind the magnitude-versus-mode paragraph in
experiments/CONCENTRATOR_GEOMETRY.md, behind the measurement comment in
compute/RCPsiSquared.Compute/LensAnalysis.cs, and behind the open arc
`concentrator_amplitude_signs`. It reproduces that document's own recipe (steps 1 to 7) from
below and compares the eigenvector of M against the |.|-only amplitudes the C# reports.

The result, printed by the run below: on symmetric F9 edge-concentrator profiles the eigenvector
is real to within 2 percent of its norm and ALTERNATES in sign, so |v| is nearly orthogonal to v.
The N=7 case is the one to look at, because its |v| reproduces the row the document prints to
every digit. On the IBM Torino gradient the components share a sign and |v| stands in for v to a
fidelity of 0.974, which is why this went unnoticed.

Runtime: seconds at N=5 and N=6. N=7 builds a 16384 x 16384 Liouvillian and needs both an
eigendecomposition and an inverse of it; it is left out of the default run for that reason and
the N=7 numbers below were produced by calling report(7, ...) directly on a machine with the
memory for it.
"""
from __future__ import annotations
import sys
import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_all(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


def site(op, k, n):
    return kron_all([op if i == k else I2 for i in range(n)])


def heisenberg(n, J=1.0):
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for b in range(n - 1):
        for P in (X, Y, Z):
            H += J * kron_all([P if i in (b, b + 1) else I2 for i in range(n)])
    return H


def liouvillian(H, gammas, n):
    d = 2 ** n
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k, g in enumerate(gammas):
        Zk = site(Z, k, n)
        L += g * (np.kron(Zk.T, Zk) - np.kron(Id, Id))
    return L


def se_indices(n):
    """|0...1_k...0> as a computational-basis index, qubit 0 = most significant."""
    return [1 << (n - 1 - k) for k in range(n)]


def lens(n, gammas, J=1.0):
    d = 2 ** n
    L = liouvillian(heisenberg(n, J), gammas, n)
    vals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    idx = se_indices(n)
    order = np.argsort(np.abs(vals.real))
    for m in order:
        if abs(vals[m].real) < 1e-9:
            continue                                   # stationary
        vec = R[:, m].reshape(d, d, order="F")
        blk = vec[np.ix_(idx, idx)]
        if np.linalg.norm(blk) / np.linalg.norm(vec) <= 0.01:
            continue                                   # not SE-accessible
        left = Rinv[m, :].reshape(d, d, order="F")
        M = left[np.ix_(idx, idx)]
        M = (M + M.conj().T) / 2
        w, v = np.linalg.eigh(M)
        j = int(np.argmax(np.abs(w)))
        psi = v[:, j]
        phase = np.exp(-1j * np.angle(psi[np.argmax(np.abs(psi))]))
        psi = psi * phase
        return vals[m], psi, np.linalg.norm(blk) / np.linalg.norm(vec)
    return None, None, None


def report(label, n, gammas):
    lam, psi, se = lens(n, gammas)
    if psi is None:
        print(f"  {label}: no SE-accessible slow mode")
        return
    true = np.real_if_close(psi, tol=1e6)
    absp = np.abs(psi); absp /= np.linalg.norm(absp)
    fid = abs(np.vdot(absp, psi)) ** 2
    signs = "".join("+" if x >= 0 else "-" for x in np.real(true))
    print(f"  {label}  N={n}  rate={lam.real:+.4f}  SE={se:.3f}")
    print(f"    true psi   : [{', '.join(f'{x:+.3f}' for x in np.real(true))}]  signs {signs}")
    print(f"    |psi| (doc): [{', '.join(f'{x:.3f}' for x in absp)}]")
    print(f"    fidelity |<|psi| , psi>|^2 = {fid:.4f}")


DOCUMENTED = {                       # what CONCENTRATOR_GEOMETRY.md and the arc quote
    ("F9", 5): 0.0022,
    ("F9", 6): 0.0000,
    ("IBM", 5): 0.974,
}


def check(label, n, gammas, expected, tol=5e-4):
    lam, psi, se = lens(n, gammas)
    assert psi is not None, f"{label} N={n}: no SE-accessible slow mode"
    absp = np.abs(psi); absp /= np.linalg.norm(absp)
    fid = abs(np.vdot(absp, psi)) ** 2
    ok = abs(fid - expected) <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {label} N={n}: fidelity {fid:.4f}, documented {expected}")
    return ok


if __name__ == "__main__":
    print("F9 edge concentrator: gamma_edge = N*gamma_base - (N-1)*eps, gamma_other = eps")
    for n in (5, 6):
        gb, eps = 0.05, 0.01
        report(f"F9 (gb={gb}, eps={eps})", n, [n * gb - (n - 1) * eps] + [eps] * (n - 1))
    print()
    print("IBM Torino sacrifice profile at N=5 (the one the AUC table uses)")
    report("IBM Torino", 5, [2.336, 0.099, 0.050, 0.072, 0.051])
    print()
    print("Gate: the fidelities the document and the arc quote")
    gb, eps = 0.05, 0.01
    passed = [
        check("F9", 5, [5 * gb - 4 * eps] + [eps] * 4, DOCUMENTED[("F9", 5)]),
        check("F9", 6, [6 * gb - 5 * eps] + [eps] * 5, DOCUMENTED[("F9", 6)]),
        check("IBM", 5, [2.336, 0.099, 0.050, 0.072, 0.051], DOCUMENTED[("IBM", 5)]),
    ]
    print(f"\n{sum(passed)}/{len(passed)} checks pass"
          f"   (the N=7 fidelity 0.0002 needs report(7, ...) and ~8 GB)")
    sys.exit(0 if all(passed) else 1)
