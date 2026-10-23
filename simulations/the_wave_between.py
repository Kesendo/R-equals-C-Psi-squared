"""What a state can see of itself, and what only a second can see.

The numbers behind reflections/ON_THE_PLACE_WHERE_THE_WAVE_RECOGNIZES_ITSELF.md.
Nothing here is new physics; the script lays known facts side by side so the
reflection can point at numbers instead of words.

1. The self-encounter S_A(t) = Tr(A† e^{tL} A) of a qubit cannot tell two
   states apart that differ only in their phase about Z, while an X probe can.
   The reason is exact: rotation about Z commutes with L (field along Z, Z
   dephasing), and the self-encounter is unchanged by any unitary symmetry of L.
   A field along X breaks that symmetry and the two self-encounters part.
2. The Bell pair: both one-site Bloch vectors are zero (each site alone is
   I/2), all content sits in the nine two-site coefficients.
3. The count of Pauli strings: 3N on one site, 4^N - 1 - 3N on two or more.
4. Reciprocity with content (F119 §7): for a real H and Z dephasing the
   generator equals its own transpose, L = L^T, exactly, so for symmetric A
   and B, Tr(B e^{tL} A) = Tr(A e^{tL} B) in the SAME system: preparation and
   readout trade places. A Dzyaloshinskii-Moriya term or amplitude damping
   breaks L = L^T, and then the swapped pairing differs. (The bare identity
   Tr(B e^{tL} A) = Tr(A e^{tL^T} B) for symmetric A, B holds for every L and
   is not tested here, since it cannot fail.)
5. The mirror: in F1's class (here an XXZ pair under Z dephasing) Π L Π⁻¹ =
   -L - 2σ with Π unitary, so for any operator A (ΠA is an operator, not a
   state) S_{ΠA}(t) = e^{-2σt} S_A(-t), backwards in the formal sense only.
6. The set S_A = S_{ΠA} for one qubit under dephasing alone (H = 0): exactly
   the pure states on the equator, with off-equator states as controls.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding="utf-8")
from framework.lindblad import lindbladian_general, lindbladian_z_dephasing  # noqa: E402
from framework.symmetry import build_pi_full  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X, Y, Z]


def vec(M):
    return M.reshape(-1)


def unvec(v, d):
    return v.reshape(d, d)


def self_encounter_blind_to_phase():
    print("1. The self-encounter cannot see the phase; a second one can")
    omega, gamma = 1.0, 0.3
    L = lindbladian_z_dephasing(0.5 * omega * Z, [gamma])

    def bloch_state(theta, phi):
        a = (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta))
        return 0.5 * (I2 + a[0] * X + a[1] * Y + a[2] * Z)

    A1, A2 = bloch_state(1.0, 0.0), bloch_state(1.0, 2.0)
    print("   t    S(A1)      S(A2)      X probe on A1  X probe on A2")
    for t in (0.0, 0.5, 1.0, 2.0, 4.0):
        U = expm(L * t)
        r1, r2 = unvec(U @ vec(A1), 2), unvec(U @ vec(A2), 2)
        print(f"   {t:<4} {np.trace(A1 @ r1).real:.6f}   {np.trace(A2 @ r2).real:.6f}"
              f"   {np.trace(X @ r1).real:+.6f}      {np.trace(X @ r2).real:+.6f}")
    Lx = lindbladian_z_dephasing(0.5 * omega * X, [gamma])
    print("   control, field along X (no longer a symmetry of the phase):")
    for t in (0.5, 2.0):
        U = expm(Lx * t)
        r1, r2 = unvec(U @ vec(A1), 2), unvec(U @ vec(A2), 2)
        print(f"   {t:<4} {np.trace(A1 @ r1).real:.6f}   {np.trace(A2 @ r2).real:.6f}")


def bell_pair_is_all_between():
    print("2. The Bell pair: nothing on either site, everything between")
    b = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    R = np.outer(b, b.conj())

    def coeff(i, j):
        return np.trace(R @ np.kron(PAULI[i], PAULI[j])).real

    print("   A alone:", [float(round(coeff(i, 0), 6)) for i in (1, 2, 3)])
    print("   B alone:", [float(round(coeff(0, j), 6)) for j in (1, 2, 3)])
    print("   between:", [[float(round(coeff(i, j), 6)) for j in (1, 2, 3)] for i in (1, 2, 3)])


def count_selves_and_between():
    print("3. Pauli strings on one site against strings on two or more")
    for N in range(1, 6):
        print(f"   N = {N}: one site {3 * N}, between {4 ** N - 1 - 3 * N}")


def two_site_generator(gammas):
    H = 0.5 * (np.kron(X, X) + np.kron(Y, Y)) + 0.35 * np.kron(Z, Z)
    return lindbladian_z_dephasing(H, gammas)


def reciprocity():
    print("4. Reciprocity: preparation and readout trade places in the same system")
    rng = np.random.default_rng(1)
    Ma = rng.normal(size=(4, 4))
    Ma = Ma + Ma.T
    Mb = rng.normal(size=(4, 4))
    Mb = Mb + Mb.T
    H = 0.5 * (np.kron(X, X) + np.kron(Y, Y)) + 0.35 * np.kron(Z, Z)
    dm = 0.3 * (np.kron(X, Y) - np.kron(Y, X))
    sm = np.array([[0, 1], [0, 0]], dtype=complex)
    loss = 0.2 * (lindbladian_general(np.zeros((4, 4)), [np.kron(sm, I2), np.kron(I2, sm)]))
    cases = (("real H, Z dephasing", two_site_generator([0.2, 0.2])),
             ("with a DM term", lindbladian_z_dephasing(H + dm, [0.2, 0.2])),
             ("with amplitude damping", two_site_generator([0.2, 0.2]) + loss))
    for name, L in cases:
        U = expm(L * 0.7)
        forward = np.trace(Mb @ unvec(U @ vec(Ma), 4)).real
        swapped = np.trace(Ma @ unvec(U @ vec(Mb), 4)).real
        print(f"   {name}: max|L - L^T| = {np.abs(L - L.T).max()}, "
              f"Tr(B e^(tL) A) = {forward:.12f}, Tr(A e^(tL) B) = {swapped:.12f}")


def mirror_self_encounter():
    print("5. The mirror image meets itself as the original running backwards")
    gammas = [0.3, 0.5]
    sigma = sum(gammas)
    L = two_site_generator(gammas)
    strings = [np.kron(PAULI[a], PAULI[b]) for a, b in itertools.product(range(4), repeat=2)]
    T = np.array([vec(S) for S in strings]).T / 2
    LP = T.conj().T @ L @ T
    Pi = build_pi_full(2)
    residual = np.abs(Pi @ LP @ Pi.conj().T + LP + 2 * sigma * np.eye(16)).max()
    print(f"   palindrome residual in this basis: {residual}")
    rng = np.random.default_rng(3)
    a = rng.normal(size=16) + 1j * rng.normal(size=16)
    b = Pi @ a
    for t in (0.4, 1.3):
        s_mirror = b.conj() @ expm(LP * t) @ b
        s_back = np.exp(-2 * sigma * t) * (a.conj() @ expm(-LP * t) @ a)
        print(f"   t = {t}: S_PiA(t) = {s_mirror:.12f}, e^(-2σt) S_A(-t) = {s_back:.12f}")


def mirror_twin_exists():
    print("6. States whose self-encounter equals their mirror image's (H = 0, one qubit)")
    print("   closed form: S_A = S_PiA at every t iff a^2 + b^2 = 1 + c^2 for A = (I + aX + bY + cZ)/2,")
    print("   which among states means exactly the pure states on the equator")
    gamma = 0.3
    L = lindbladian_z_dephasing(np.zeros((2, 2), dtype=complex), [gamma])
    T = np.array([vec(S) for S in PAULI]).T / np.sqrt(2)
    LP = T.conj().T @ L @ T
    Pi = build_pi_full(1)
    cases = (("(I + Y)/2, equator", 0.5 * (I2 + Y)),
             ("(I + (X + Y)/sqrt2)/2, equator", 0.5 * (I2 + (X + Y) / np.sqrt(2))),
             ("(I + 0.5 X)/2, inside", 0.5 * (I2 + 0.5 * X)),
             ("(I + Z)/2, pole", 0.5 * (I2 + Z)))
    for name, A in cases:
        a = T.conj().T @ vec(A)
        b = Pi @ a
        diffs = [abs((a.conj() @ expm(LP * t) @ a) - (b.conj() @ expm(LP * t) @ b)) for t in (0.5, 1.7)]
        print(f"   {name}: |S_A - S_PiA| at t = 0.5, 1.7: "
              f"{diffs[0]:.3e}, {diffs[1]:.3e}")


if __name__ == "__main__":
    self_encounter_blind_to_phase()
    bell_pair_is_all_between()
    count_selves_and_between()
    reciprocity()
    mirror_self_encounter()
    mirror_twin_exists()
