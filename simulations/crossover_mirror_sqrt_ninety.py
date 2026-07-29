#!/usr/bin/env python3
"""The crossover mirror is the canonical Π turned by √(NinetyDegreeMirror).

This locks down the angle-branch reading of the XZ+YZ / ZX+ZY locality result
(experiments/PI_OPERATOR_ENTANGLEMENT.md). The local mirror found there, the continuous
per-site map M with Π = M^⊗N, is not new structure: it sits on the framework's Pi2-Z₄
rotational axis (NinetyDegreeMirrorMemoryClaim / F91), at the 45° point the discrete Z₄
skips.

Two bit-exact facts:

  PART A (the √ identity, pure 4×4 algebra).
    S := M · Π⁻¹ (the turn from the canonical mirror to the crossover mirror) is
    block-diagonal: a pure rotation, no dark{I,Z}↔light{X,Y} swap. In the light plane it
    is a 135° turn (−45° once M's free sign is taken out; (−S)² = S² either way), and
    S_light² = [[0,1],[−1,0]], a σ_x↔σ_y 90° rotation. So S_light is a square root of a
    90° turn: the crossover mirror is Π turned by HALF the 90° angle-anchor, the bisector
    between P1 (0°) and P4 (90°). M itself is order 4 (M² = −I), carrying the i⁴ = 1
    algebra of the angle-anchor.

    That reading needs the matched pair, and there are four. The Pauli-basis Liouvillian
    is exactly real, so Π and conj(Π) both mirror it, and so do M and conj(M). Pairing
    M with conj(P1) (what `crossover_map()` fixes, and what PART A uses) or conj(M) with
    P1 gives the rotation above; pairing M with P1, or conj(M) with conj(P1), gives a
    REFLECTION (det −1, a flip) whose square is the identity, which a phase Π → iΠ
    turns into −I. The 0°/90° marks for P1 and P4 are a |angle| reading and survive
    only in that branch-blind register: no single pairing makes all three rotations. The same
    phase on the matched pair keeps the 90° and reverses its sense. So the 90° is a
    property of the matched pair, not of the geometry, and its direction is free. The
    convention-free statement lives in the proof, as Ad_V² = NinetyDegreeMirror.

  PART B (the continuous dial).
    For the weighted bond H = a·XZ + b·YZ a product mirror exists at every (a,b)
    (residual ~1e-12), and its light-plane turn tracks the bond direction: M sends I to
    the direction PERPENDICULAR to (a, b). The discrete Z₄ samples this dial at the
    cardinal points (P1 at the pure-YZ end, P4 at the pure-XZ end); the crossover a=b
    sits at the 45° middle.

The clean √-of-90° identity holds exactly in the light plane {X,Y}, where the X/Y
diplexer conflict lives; the full 4×4 carries extra dark-plane phase bookkeeping (Π's
Z↦iY factors), so M is not a literal discrete-Z₄ element.

Run:  PYTHONIOENCODING=utf-8 python simulations/crossover_mirror_sqrt_ninety.py
"""
from __future__ import annotations

import sys
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, 'simulations')
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

np.set_printoptions(precision=4, suppress=True)

# Framework Pauli order [I, X, Z, Y]; dark {I,Z} = {0,2}, light {X,Y} = {1,3}.
DARK, LIGHT = [0, 2], [1, 3]
s = 1 / np.sqrt(2)


def mat(entries):
    M = np.zeros((4, 4), dtype=complex)
    for (r, c), v in entries.items():
        M[r, c] = v
    return M


# canonical Π (P1): I↦X, X↦I, Z↦iY, Y↦iZ. Equals framework build_pi_full(1).
P1 = mat({(1, 0): 1, (0, 1): 1, (3, 2): 1j, (2, 3): 1j})
# The conjugate mirror, c = diag(1,1,1,−1) on [I,X,Z,Y]. Both P1 and this are mirrors
# of every Pauli-basis Liouvillian (which is real), so they differ by convention only.
# This is the row-stack representative, the one that matches `crossover_map()`.
C_TWIST = np.diag([1, 1, 1, -1]).astype(complex)
P1_ROW = C_TWIST @ P1 @ C_TWIST


def part_a():
    print("PART A  the √(NinetyDegreeMirror) identity")
    M = fw.crossover_map()                       # the closed-form crossover mirror (row-stack)
    S = M @ np.linalg.inv(P1_ROW)                 # turn from canonical Π to crossover Π

    off = sum(abs(S[r, c]) + abs(S[c, r]) for r in DARK for c in LIGHT)
    print(f"  S = M·Π⁻¹ block-diagonal (pure rotation, no dark↔light swap)? "
          f"{off < 1e-9}  (off-block norm {off:.1e})")

    S_light = S[np.ix_(LIGHT, LIGHT)]
    # Magnitude readout: confined to [0°, 90°], so it reads 45° for any of 45/135/225/315.
    # The matrix actually built is a 135° rotation; since M and −M are both mirrors, the
    # branch is a sign choice. What is falsifiable is S_light² below, not this number.
    ang = np.degrees(np.arctan2(abs(S_light[1, 0]), abs(S_light[0, 0])))
    print(f"  S light-plane |angle| = {ang:.2f}°  (branch-blind; the signed turn is 135°)")

    # The σ_x↔σ_y 90° turn this pairing produces: X↦−Y, Y↦X.
    ninety = np.array([[0, 1], [-1, 0]], dtype=complex)
    resid = np.linalg.norm(S_light @ S_light - ninety)
    print(f"  S_light² = a σ_x↔σ_y 90° turn?  residual {resid:.2e}")
    print(f"      ⟹ S_light is a square root of a 90° turn" if resid < 1e-9
          else "      (identity FAILS)")
    # The other pairing, M against the unconjugated Π, is not the opposite turn: it is a
    # reflection, and its square is the identity. So the 90° belongs to the matched pair.
    S_mixed = M @ np.linalg.inv(P1)
    SM = S_mixed[np.ix_(LIGHT, LIGHT)]
    mixed = np.linalg.norm(SM @ SM - np.eye(2))
    print(f"      (pair M with the other Π and S² = I instead, a reflection, "
          f"residual {mixed:.2e})")
    print(f"  M² = −I (M is an 'i', order 4)? {np.allclose(M @ M, -np.eye(4), atol=1e-12)}")
    print()


def _L_pauli_weighted(a, b, N=2):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        H = H + a * fw.site_op(N, i, 'X') @ fw.site_op(N, i + 1, 'Z')
        H = H + b * fw.site_op(N, i, 'Y') @ fw.site_op(N, i + 1, 'Z')
    L = fw.lindbladian_z_dephasing(H, [0.5] * N)
    # order='C' to match PART A: this script's whole subject is that the stacking is
    # load-bearing, and the default 'F' would search for the column-stack representative
    # (in which the I-image comes out PARALLEL to the bond at a=b, not perpendicular).
    V = fw._vec_to_pauli_basis_transform(N, order='C')
    return (V.conj().T @ L @ V) / d, N * 0.5


def _block_M(p):
    B = np.array([[p[0] + 1j * p[1], p[2] + 1j * p[3]], [p[4] + 1j * p[5], p[6] + 1j * p[7]]])
    C = np.array([[p[8] + 1j * p[9], p[10] + 1j * p[11]], [p[12] + 1j * p[13], p[14] + 1j * p[15]]])
    M = np.zeros((4, 4), dtype=complex)
    for ri, r in enumerate(LIGHT):
        for ci, c in enumerate(DARK):
            M[r, c] = B[ri, ci]
    for ri, r in enumerate(DARK):
        for ci, c in enumerate(LIGHT):
            M[r, c] = C[ri, ci]
    return M


def _find_weighted_M(a, b, restarts=150, seed=0):
    Lp, sig = _L_pauli_weighted(a, b)
    rng = np.random.default_rng(seed)
    eye = np.eye(16)

    def obj(p):
        M = _block_M(p)
        Pi = np.kron(M, M)
        try:
            r = Pi @ Lp @ np.linalg.inv(Pi) + Lp + 2 * sig * eye
        except np.linalg.LinAlgError:
            return 1e6
        return float(np.linalg.norm(r))

    best = (np.inf, None)
    for _ in range(restarts):
        res = minimize(obj, rng.standard_normal(16), method='Nelder-Mead',
                       options={'maxiter': 20000, 'xatol': 1e-9, 'fatol': 1e-12})
        if res.fun < best[0]:
            best = (res.fun, _block_M(res.x))
        if best[0] < 1e-9:
            break
    return best


def part_b():
    print("PART B  the continuous dial (mirror turns with the bond a·XZ + b·YZ)")
    print(f"  {'a':>5}{'b':>6}{'residual':>12}{'bond':>8}{'I→light':>10}{'perpendicular?':>16}")
    for a, b in [(1, 1), (1, 0.5), (0.5, 1), (1, 0)]:
        val, M = _find_weighted_M(a, b, seed=int(10 * a + 100 * b))
        if M is None or val > 1e-6:
            print(f"  {a:5.2f}{b:6.2f}{val:12.2e}   (no product found)")
            continue
        img = np.degrees(np.arctan2(abs(M[3, 0]), abs(M[1, 0])))   # I → (X,Y) angle
        bond = np.degrees(np.arctan2(abs(b), abs(a)))
        perp = "yes" if abs(img - (90 - bond)) < 1.0 else "no"
        print(f"  {a:5.2f}{b:6.2f}{val:12.2e}{bond:8.1f}°{img:9.1f}°{perp:>16}")
    print("  (I↦light is perpendicular to the bond; the dial rotates continuously with (a,b))")
    print()


def main():
    print("Crossover mirror = canonical Π · √(NinetyDegreeMirror)  (Pauli order [I,X,Z,Y], γ=0.5)\n")
    part_a()
    part_b()
    print("Reading: the local crossover mirror is the √ of the framework's 90° angle-anchor,")
    print("the 45° middle of the Pi2-Z₄ that P1/P4 mark at 0°/90°, a reading that holds in")
    print("the branch-blind |angle| register only. Same axis, half the angle.")


if __name__ == "__main__":
    main()
