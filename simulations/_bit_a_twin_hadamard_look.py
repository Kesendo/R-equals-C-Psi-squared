"""_bit_a_twin_hadamard_look.py - looking at the general connection (Tom, 2026-05-29).

The bit_a twins of the cubic Z₂³ polarity architecture are not 52 independent
derivations. PROOF_BIT_A_TWIN_VIA_HADAMARD.md states the general connection: a BitB
claim whose content is operator-space (Π spectrum / identities) or Lindblad-spectral
(Re(λ), the dephasing-dissipator popcount) owns its bit_a twin via the X↔Z duality.
We typed it (ad9767b: 9 Absorption-descendants classified CoveredByHadamardDuality)
and proved the mechanisms elsewhere, but never LOOKED at the duality carrying a concrete
claim end to end. This is that look, bit-exact.

The realizer for the Lindblad case is the genuine Hilbert-space global Hadamard
U_H = H^⊗N (H X H = Z per site). The claims to transport: the Absorption Theorem
Re(λ) = −2γ·popcount and the full Liouvillian spectrum.

  A. U_H realizes Z↔X on the Hilbert space (per site, bit-exact).
  B. The full Liouvillian transports: Spec(L_Z(H)) = Spec(L_X(U_H H U_H†)) bit-exact,
     for any Hermitian H (so any spectral BitB claim owns its bit_a twin).
  C. The Absorption content itself: Re(λ) = −2γ·popcount for Z-dephasing equals that
     for X-dephasing (the bit_a twin of the 9 Absorption-descendants), bit-exact.
  D. The lift caveat: the spectrum is even equal for Y, but only U_H (Z↔X) is a
     Lindblad-preserving Hilbert unitary; the Klein Z↔Y element D = diag((−1)^{n_Y})
     lives on the 4^N operator space, not the 2^N Hilbert space. The duality is Z↔X.

Run: python simulations/_bit_a_twin_hadamard_look.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import numpy as np

GAMMA = 0.3
TOL = 1e-10
_ok = []

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}
HAD = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0)


def report(name, cond, extra=""):
    _ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{extra}")


def popcount(x):
    return bin(x).count("1")


def site_op(N, site, letter):
    op = PAULI[letter] if site == 0 else PAULI["I"]
    for s in range(1, N):
        op = np.kron(op, PAULI[letter] if s == site else PAULI["I"])
    return op


def two_site(N, a, b, la, lb):
    ops = ["I"] * N
    ops[a], ops[b] = la, lb
    out = PAULI[ops[0]]
    for s in range(1, N):
        out = np.kron(out, PAULI[ops[s]])
    return out


def hadamard_N(N):
    U = HAD
    for _ in range(N - 1):
        U = np.kron(U, HAD)
    return U


def heisenberg(N):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N - 1):
        for L in ("X", "Y", "Z"):
            H = H + two_site(N, a, a + 1, L, L)
    return H


def lindbladian(H, c_list, gamma):
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, I) - np.kron(I, H.T))
    for c in c_list:
        cdc = c.conj().T @ c
        L = L + gamma * (np.kron(c, c.conj()) - 0.5 * (np.kron(cdc, I) + np.kron(I, cdc.T)))
    return L


def deph(letter, N):
    return [site_op(N, l, letter) for l in range(N)]


def spec(L):
    return np.sort_complex(np.round(np.linalg.eigvals(L), 9))


def main():
    print("=" * 78)
    print("THE bit_a TWIN VIA HADAMARD - looking at the general connection, bit-exact")
    print("=" * 78)

    for N in (2, 3):
        d = 2 ** N
        U = hadamard_N(N)
        print(f"\nN={N}:")

        # A. U_H realizes Z<->X on the Hilbert space (and Y -> -Y), per site
        worst_zx = worst_xz = worst_yy = 0.0
        for l in range(N):
            Z, X, Y = site_op(N, l, "Z"), site_op(N, l, "X"), site_op(N, l, "Y")
            worst_zx = max(worst_zx, np.max(np.abs(U @ Z @ U.conj().T - X)))
            worst_xz = max(worst_xz, np.max(np.abs(U @ X @ U.conj().T - Z)))
            worst_yy = max(worst_yy, np.max(np.abs(U @ Y @ U.conj().T + Y)))
        report("A: U_H Z_l U_H† = X_l  and  U_H X_l U_H† = Z_l  (Hadamard swaps Z↔X)",
               max(worst_zx, worst_xz) < TOL, f"   max|d| = {max(worst_zx, worst_xz):.1e}")
        report("A: U_H Y_l U_H† = −Y_l  (Y flips: the third Klein face)", worst_yy < TOL)

        # B. The full Liouvillian transports: Spec(L_Z(H)) = Spec(L_X(U_H H U_H†))
        #    for a non-Hadamard-symmetric H (Heisenberg + a Z-field), so the bit_a twin
        #    holds for any spectral claim, any Hermitian H (the H-image is its X-field).
        H = heisenberg(N) + 0.7 * sum(site_op(N, l, "Z") for l in range(N))
        H_img = U @ H @ U.conj().T
        Lz = lindbladian(H, deph("Z", N), GAMMA)
        Lx_img = lindbladian(H_img, deph("X", N), GAMMA)
        dev_b = float(np.max(np.abs(spec(Lz) - spec(Lx_img))))
        report("B: Spec(L_Z(H)) = Spec(L_X(U_H H U_H†))  (the bit_a twin of any spectral claim)",
               dev_b < TOL, f"   H = Heisenberg + 0.7·ΣZ → X-field image; max|Δspec| = {dev_b:.1e}")

        # C. The Absorption content itself: Re(λ) = −2γ·popcount, equal for Z and X.
        Lz0 = lindbladian(np.zeros((d, d), complex), deph("Z", N), GAMMA)
        Lx0 = lindbladian(np.zeros((d, d), complex), deph("X", N), GAMMA)
        re_pred = np.sort(np.array([-2 * GAMMA * popcount(i ^ j) for i in range(d) for j in range(d)]))
        re_z = np.sort(np.round(np.linalg.eigvals(Lz0).real, 9))
        re_x = np.sort(np.round(np.linalg.eigvals(Lx0).real, 9))
        report("C: Re(λ) of Z-dephasing dissipator = −2γ·popcount(i⊕j)  (Absorption, bit_b)",
               float(np.max(np.abs(re_z - re_pred))) < TOL)
        report("C: Re(λ) of X-dephasing dissipator = the same multiset  (the bit_a twin)",
               float(np.max(np.abs(re_x - re_pred))) < TOL)

        # D. Lift caveat: the spectrum is even equal for Y (popcount is basis-blind),
        #    but only U_H (Z↔X) is a Lindblad-preserving Hilbert unitary. The Klein
        #    Z↔Y element D = diag((−1)^{n_Y}) lives on the 4^N operator space, not the
        #    2^N Hilbert space, so the duality that carries a Lindblad claim is Z↔X.
        Ly0 = lindbladian(np.zeros((d, d), complex), deph("Y", N), GAMMA)
        re_y = np.sort(np.round(np.linalg.eigvals(Ly0).real, 9))
        report("D: Re(λ) for Y-dephasing also = −2γ·popcount  (spectrum is basis-blind)",
               float(np.max(np.abs(re_y - re_pred))) < TOL)
        print(f"      lift: U_H is a {d}×{d} Hilbert unitary (Z↔X); the Klein Z↔Y element D "
              f"is {d*d}×{d*d} on operator space → no Lindblad-preserving lift to Y.")

    n_ok, n_tot = sum(_ok), len(_ok)
    print("\n" + "=" * 78)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'CHECK'})")
    print("=" * 78)
    print("""
The general connection, looked at:
  The global Hadamard U_H = H^⊗N swaps Z↔X on the Hilbert space, carries the Z-dephasing
  Liouvillian to the X-dephasing one for ANY Hermitian H, and so transports every
  operator-space and Lindblad-spectral BitB claim to its bit_a twin verbatim. The
  Absorption rate −2γ·popcount is the same in both bases. The nine CoveredByHadamardDuality
  corollaries are not nine bookkeeping entries; they are one unitary acting, and here it
  acts, bit-exact. The Y leg shares the spectrum but not the Hilbert lift, so the duality
  is Z↔X, exactly as the theorem bounds it.
""")


if __name__ == "__main__":
    main()
