"""f108_bita_d_search.py - closing F108 open question (ii): the "bit_a-D" for Π_5b.

PROOF_F108_KLEIN_V4_EQUIVALENCE.md left open (the answer now sits in its §(d.4)):
  "Part 2's BitA-twin status: ... Could the bit_a-axis version of D (a Z↔X swap analog)
   intertwine Π_5b(Z) ↔ Π_5b(X) directly on operator space? [The canonical Q_zx lands on
   the X mirror's other orientation.] A search for a different operator-space involution that intertwines
   Π_5b(Z) ↔ Π_5b(X) directly (analogous to D for Z↔Y) is open. The natural ansatz
   would be an operator-space involution sitting on the bit_a axis with the Y phase fixed
   but the X/Z phases flipped; we have not formally enumerated."

This enumerates, and the answer is already in the Welle-12 Klein-V₄ group: the
involution is H, the pure X↔Z basis swap with I and Y fixed (Pi2KleinV4DephaseSwapGroup).
The F108 verifier tested Q_zx on the Z→X pairing and H on the Y→X pairing (both land on
the X mirror's other orientation), but never tested H on the Z→X pairing - exactly where it
lands on Π_5b(X) itself. H is the
"Y phase fixed, X/Z swapped" ansatz the proof guessed.

  H · Π_5b(Z) · H = Π_5b(X)   bit-exact, universal N.

Why this is the bit_a-D and not a Hilbert move: a Hilbert unitary swapping X↔Z must
send Y → −Y (since XZ = −iY forces ZX = +iY), i.e. it would be Q_zx, the Hadamard one,
which lands Π_5b(Z) on (−1)^N · Π_5b(X)⁻¹, the X mirror in its other orientation. H fixes Y while swapping X,Z, which
is NOT a Pauli automorphism - a pure operator-space involution, exactly like D.

Per-site maps from Pi5BilinearOperator.ActOnLetter, basis [I, X, Y, Z]:
  M_Z: I→+X, X→−I, Y→+iZ, Z→−iY
  M_X: I→+Z, Z→−I, X→−iY, Y→+iX
  M_Y: I→+X, X→−I, Y→−iZ, Z→+iY
Klein-V₄ from Pi2KleinV4DephaseSwapGroup (basis [I, X, Y, Z]):
  D = diag(1,1,−1,1) (negate Y) ; H = swap(X,Z) fix(I,Y) ; Q_zx = H·D

Run: python simulations/f108_bita_d_search.py
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

TOL = 1e-12
_ok = []

# basis order [I, X, Y, Z]
I_, X_, Y_, Z_ = 0, 1, 2, 3
LETTER = {I_: "I", X_: "X", Y_: "Y", Z_: "Z"}

# Per-site Π_5bilinear letter maps: dephase letter -> {in_letter: (out_letter, phase)}
PI5 = {
    "Z": {I_: (X_, 1), X_: (I_, -1), Y_: (Z_, 1j), Z_: (Y_, -1j)},
    "X": {I_: (Z_, 1), Z_: (I_, -1), X_: (Y_, -1j), Y_: (X_, 1j)},
    "Y": {I_: (X_, 1), X_: (I_, -1), Y_: (Z_, -1j), Z_: (Y_, 1j)},
}


def report(name, cond, extra=""):
    _ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{extra}")


def site_from_map(m):
    """4x4 superoperator M[out, in] = phase from a letter map."""
    M = np.zeros((4, 4), dtype=complex)
    for src, (dst, ph) in m.items():
        M[dst, src] = ph
    return M


def kron_power(M, N):
    out = M
    for _ in range(N - 1):
        out = np.kron(out, M)
    return out


def fro(A):
    return float(np.linalg.norm(A, "fro"))


def describe_letter_action(M):
    """Human-readable per-site action of a 4x4 signed-permutation superoperator."""
    parts = []
    for src in (I_, X_, Y_, Z_):
        col = M[:, src]
        nz = np.nonzero(np.abs(col) > TOL)[0]
        if len(nz) == 1:
            dst = int(nz[0])
            ph = col[dst]
            ph_s = {1: "+", -1: "−", 1j: "+i", -1j: "−i"}.get(
                complex(np.round(ph)), f"{ph:.2g}")
            parts.append(f"{LETTER[src]}→{ph_s}{LETTER[dst]}")
        else:
            parts.append(f"{LETTER[src]}→(mixed)")
    return ", ".join(parts)


def main():
    print("=" * 80)
    print("F108 open question (ii): the bit_a-D for Π_5b - is H the X↔Z involution?")
    print("=" * 80)

    M = {d: site_from_map(PI5[d]) for d in ("Z", "X", "Y")}

    # Klein-V₄ elements, basis [I, X, Y, Z]
    D = np.diag([1, 1, -1, 1]).astype(complex)          # negate Y
    H = np.zeros((4, 4), dtype=complex)                  # swap X,Z ; fix I,Y
    H[I_, I_] = 1; H[Y_, Y_] = 1; H[Z_, X_] = 1; H[X_, Z_] = 1
    Qzx = H @ D

    print("\nPer-site operators (basis I, X, Y, Z):")
    for d in ("Z", "X", "Y"):
        print(f"  Π_5b({d}): {describe_letter_action(M[d])}")
    print(f"  D     : {describe_letter_action(D)}")
    print(f"  H     : {describe_letter_action(H)}")
    print(f"  Q_zx  : {describe_letter_action(Qzx)}")

    # ---- A. Sanity: the per-site operators are what the C# source says ----
    print("\nA. Sanity on the per-site Π_5b operators")
    for d in ("Z", "X", "Y"):
        u = np.allclose(M[d].conj().T @ M[d], np.eye(4), atol=TOL)
        o4 = np.allclose(np.linalg.matrix_power(M[d], 4), np.eye(4), atol=TOL)
        report(f"M_{d} unitary and order-4 (M⁴=I)", u and o4)
    report("M_Z² = diag(−1,−1,+1,+1) on (I,X,Y,Z)",
           np.allclose(M["Z"] @ M["Z"], np.diag([-1, -1, 1, 1]), atol=TOL))
    for nm, G in (("D", D), ("H", H), ("Q_zx", Qzx)):
        inv = np.allclose(G @ G, np.eye(4), atol=TOL)
        uni = np.allclose(G.conj().T @ G, np.eye(4), atol=TOL)
        report(f"{nm} is a unitary involution ({nm}²=I)", inv and uni)

    # ---- B. Known landscape: D does Z↔Y; the two pairings that land on the other orientation ----
    print("\nB. Reproduce the known F108 Klein-V₄ landscape")
    report("D · M_Z · D = M_Y  (D does Z↔Y for Π_5b, known positive)",
           np.allclose(D @ M["Z"] @ D, M["Y"], atol=TOL),
           f"   gap = {fro(D @ M['Z'] @ D - M['Y']):.2e}")
    gap_qzx = fro(Qzx @ M["Z"] @ Qzx - M["X"])
    report("Q_zx · M_Z · Q_zx ≠ M_X  (it lands on the other orientation, below)",
           gap_qzx > 0.5, f"   gap = {gap_qzx:.3f} (proof reports 2.0 at max-norm)")
    gap_hyx = fro(H @ M["Y"] @ H - M["X"])
    report("H · M_Y · H ≠ M_X, and = -M_X^-1 (the other orientation)",
           gap_hyx > 0.5 and np.allclose(H @ M["Y"] @ H, -M["X"].conj().T, atol=TOL), f"   gap = {gap_hyx:.3f}")

    # ---- C. THE FINDING: H does Z→X for Π_5b (the pairing the proof never tested) ----
    print("\nC. The finding: H on the Z→X pairing")
    gap_hzx = fro(H @ M["Z"] @ H - M["X"])
    report("H · M_Z · H = M_X  (H IS the bit_a-D: X↔Z swap, Y fixed)",
           gap_hzx < TOL, f"   gap = {gap_hzx:.2e}")
    # and its involution partner direction
    report("H · M_X · H = M_Z  (involution, so the swap is symmetric)",
           np.allclose(H @ M["X"] @ H, M["Z"], atol=TOL))

    # ---- D. H is operator-space-only: NOT a Hilbert/Pauli automorphism ----
    print("\nD. H is operator-space-only (no Hilbert-space lift)")
    # A Hilbert unitary V with V X V† = Z, V Z V† = X must satisfy V Y V† = −Y,
    # because Y = i·X·Z  ⇒  V Y V† = i·(VXV†)(VZV†) = i·Z·X = i·(+iY) = −Y.
    # H fixes Y (H: Y→+Y), so no such V exists: H is a pure operator-space involution.
    Xp = np.array([[0, 1], [1, 0]], complex)
    Yp = np.array([[0, -1j], [1j, 0]], complex)
    Zp = np.array([[1, 0], [0, -1]], complex)
    forced = 1j * (Zp @ Xp)   # what Y must map to if X↔Z (= i·Z·X)
    report("a Hilbert X↔Z swap forces Y → −Y (so it would be Q_zx, not H)",
           np.allclose(forced, -Yp, atol=TOL))
    report("H instead fixes Y (Y→+Y) ⇒ H has no Hilbert lift, operator-space-only",
           abs(H[Y_, Y_] - 1) < TOL)

    # ---- E. Universal N: the per-site identity tensors up (Kronecker) ----
    print("\nE. Universal N: (H^⊗N)·Π_5b(Z)·(H^⊗N) = Π_5b(X)")
    for N in (1, 2, 3):
        HN = kron_power(H, N)
        PZ = kron_power(M["Z"], N)
        PX = kron_power(M["X"], N)
        gap = fro(HN @ PZ @ HN - PX)
        report(f"N={N}: H^⊗{N} · Π_5b(Z) · H^⊗{N} = Π_5b(X)  (4^{N}={4**N} dim)",
               gap < TOL, f"   gap = {gap:.2e}")

    # ---- F. The corrected Klein-V₄ picture for Π_5b ----
    print("\nF. The Π_5b Klein-V₄ action, corrected")
    print("   D : Z↔Y  (operator-space, proven Welle 14)")
    print("   H : Z↔X  (operator-space, THIS script) - the open (ii) answer")
    print("   Q_zx=H·D : X↔Y, and Z → -M_X^-1, the X mirror in its other orientation")
    report("Q_zx · M_X · Q_zx = M_Y  (Q_zx does X↔Y on Π_5b)", np.allclose(Qzx @ M["X"] @ Qzx, M["Y"], atol=TOL))
    # show Q_zx·M_Z·Qzx is -M_X^-1, the X mirror in its other orientation
    qzx_img = Qzx @ M["Z"] @ Qzx
    other_orientation = np.allclose(qzx_img, -M["X"].conj().T, atol=TOL)
    report("Q_zx·M_Z·Q_zx = -M_X^-1 (the X mirror in its other orientation)", other_orientation)

    n_ok, n_tot = sum(_ok), len(_ok)
    print("\n" + "=" * 80)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'CHECK'})")
    print("=" * 80)
    print("""
F108 open question (ii) closes POSITIVELY:
  The operator-space involution that intertwines Π_5b(Z) ↔ Π_5b(X) is H, the third
  Klein-V₄ element (Welle 12) - the pure X↔Z basis swap with I and Y fixed. It is
  exactly the "Y phase fixed, X/Z swapped" ansatz the proof guessed, and it was already
  sitting in the group. Testing each Klein element on its canonical-Π pairing (Q_zx on Z→X, H on Y→X)
  lands on the X mirror's other orientation;
  on Π_5b the roles of H and Q_zx are swapped: H does Z↔X, Q_zx (= H·D) does X↔Y, and
  each element sends the remaining variant to -M_X^-1, the X mirror in its other
  orientation. The operator-space Klein-V₄ thus acts simply transitively on the four
  oriented mirrors, as it does on the canonical {Π_Z, Π_Y, Π_X, Π_X^-1}.
""")


if __name__ == "__main__":
    main()
