"""Heisenberg-companion to the XY structural ceiling: the star survivor darkness (2026-06-20, gate-first).

RESOLVED (no bug). The typed structural-ceiling / star-frozen-seam framework is, BY DESIGN, the XY
network (hopping only): StructuralCeilingClaim says "g2 = strict_gap/2gamma of an XY network",
StructuralCeilingWitness.SectorH and StarFrozenSeamWitness.FullSlowest both build XY hopping with NO ZZ
diagonal. For the XY star that gives g2 = 4/(N-1), correctly = the XY star survivor darkness.

The rest of the survivor framework (SurvivalIncompletenessMirrorClaim, SURVIVOR_FLIP_AND_REFLECTION_ODD)
is the canonical HEISENBERG model (XX+YY+ZZ). The Heisenberg star survivor is a DIFFERENT number:
4/N, not 4/(N-1). The ZZ potential (hub V0 = -(N-1), leaves Vl = N-3 in the single-excitation sector)
shifts which ad_H-kernel commutant is darkest. This file is the Heisenberg companion: it pins 4/N and
shows the ZZ term is the entire difference. (The earlier star_commutant_gate.py / star_frozen_gate.py
used Heisenberg and read 4/N -- correctly; they were testing Heisenberg against the XY docs, the model
mismatch that opened this thread.)

What holds in BOTH models (the StarFrozenSeam qualitative content): the survivor is frozen (|Im|=0) at
every Q for N>=5; it is the [H,rho]=0 commutant only in the high-Q LIMIT (||[H,rho]||/||rho|| ~ 1/Q,
NOT zero at finite Q); it sits at the (1,1)/(N-1,N-1) popcount boundary. Only the darkness VALUE is
model-specific: XY 4/(N-1), Heisenberg 4/N.

GATES (a firing gate is the find; do not loosen):
  G1  XY star (1,1) survivor darkness -> 4/(N-1)   (reproduces the docs' g2; the docs are right for XY)
  G2  Heisenberg star (1,1) survivor darkness -> 4/N   (the canonical-model companion)
  G3  the ad_H-kernel darkest commutant flips 4/(N-1) (XY) -> 4/N (Heisenberg) when the ZZ potential is
      added -- the ZZ term is the whole difference
  G4  both frozen (|Im|~0) and approach the exact commutant as Q grows (||[H,rho]|| ~ 1/Q)
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def star_H(N, J, with_zz):
    """Single-excitation (1,1) star Hamiltonian: 2J*A (XX+YY hopping) + J*diag(V) (ZZ) if Heisenberg."""
    A = np.zeros((N, N))
    for l in range(1, N):
        A[0, l] = A[l, 0] = 1.0
    V = np.zeros(N)
    if with_zz:
        V[0] = -(N - 1)                       # hub ZZ potential (Heisenberg)
        for l in range(1, N):
            V[l] = N - 3                       # each leaf
    return 2 * J * A + J * np.diag(V)


def survivor_11(N, J, g, with_zz):
    """Slowest (1,1)-sector Liouvillian mode: returns (<n_XY>=-Re/2g, |Im|, ||[H,rho]||/||rho||)."""
    H = star_H(N, J, with_zz)
    d = N
    Id = np.eye(d)
    Lh = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    Ld = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            Ld[i + d * j, i + d * j] = 0.0 if i == j else -4 * g   # (1,1) Z-dephasing: hamming=2 -> -4g
    w, V_ = np.linalg.eig(Lh + Ld)
    re = w.real
    k = max((kk for kk in range(len(w)) if re[kk] < -1e-9), key=lambda kk: re[kk])
    raw = V_[:, k].reshape(d, d, order="F")
    rho = (raw + raw.conj().T) / 2
    if np.linalg.norm(rho) < 1e-9:
        rho = (raw - raw.conj().T) / 2j
    nr = np.linalg.norm(rho)
    comm = np.linalg.norm(H @ rho - rho @ H) / nr
    return -w[k].real / (2 * g), abs(w[k].imag), comm


def kernel_darkest(N, with_zz):
    """The witness's object: min nonzero eigenvalue of N_XY (hamming, =2 off-diag) on the ad_H kernel."""
    H = star_H(N, 1.0, with_zz)
    d = N
    Id = np.eye(d)
    adH = np.kron(H, Id) - np.kron(Id, H)
    wO, VO = np.linalg.eigh(adH)
    ker = VO[:, np.abs(wO) < 1e-7]
    nxy = np.array([0.0 if (idx // d) == (idx % d) else 2.0 for idx in range(d * d)])
    Ntil = ker.T @ (nxy[:, None] * ker)
    ew = np.linalg.eigvalsh(Ntil)
    return min(e for e in ew if e > 1e-7)


def main():
    print("=== star survivor darkness: XY (docs) vs Heisenberg (canonical) -- the ZZ term is the difference ===\n",
          flush=True)
    g1 = g2 = g3 = g4 = True
    print(f"{'N':>3} {'XY <nXY>':>9} {'4/(N-1)':>8} {'Heis <nXY>':>11} {'4/N':>7} {'XY ker':>7} {'Heis ker':>9}",
          flush=True)
    for N in (5, 6, 7, 8):
        g = 0.002
        xy, _, _ = survivor_11(N, 1.0, g, False)
        hb, _, _ = survivor_11(N, 1.0, g, True)
        xk = kernel_darkest(N, False)
        hk = kernel_darkest(N, True)
        print(f"{N:>3} {xy:>9.4f} {4/(N-1):>8.4f} {hb:>11.4f} {4/N:>7.4f} {xk:>7.4f} {hk:>9.4f}", flush=True)
        g1 &= abs(xy - 4 / (N - 1)) < 1e-3
        g2 &= abs(hb - 4 / N) < 1e-3
        g3 &= abs(xk - 4 / (N - 1)) < 1e-3 and abs(hk - 4 / N) < 1e-3

    # G4: frozen + commutant-in-the-limit (Heisenberg, N=6, Q-sweep)
    print(f"\n[G4] Heisenberg N=6 Q-sweep: frozen (|Im|~0) + commutant only in the high-Q limit (~1/Q)",
          flush=True)
    print(f"  {'Q':>8} {'<nXY>':>8} {'|Im|':>10} {'||[H,rho]||':>12}", flush=True)
    prev = None
    for Q in (1.0, 5.0, 20.0, 100.0, 500.0):
        nxy, im, comm = survivor_11(6, 1.0, 1.0 / Q, True)
        print(f"  {Q:>8.1f} {nxy:>8.4f} {im:>10.1e} {comm:>12.4f}", flush=True)
        g4 &= im < 1e-6
        prev = comm
    g4 &= prev < 0.05    # ||[H,rho]|| -> 0 at high Q (commutant only in the limit)

    print("\nGATES:", flush=True)
    print(f"  G1 XY -> 4/(N-1) (docs correct for XY)        [{'ok' if g1 else 'FIRED'}]", flush=True)
    print(f"  G2 Heisenberg -> 4/N (canonical companion)    [{'ok' if g2 else 'FIRED'}]", flush=True)
    print(f"  G3 ad_H-kernel darkest: XY 4/(N-1), Heis 4/N  [{'ok' if g3 else 'FIRED'}]", flush=True)
    print(f"  G4 frozen + commutant only in the high-Q limit[{'ok' if g4 else 'FIRED'}]", flush=True)
    ok = g1 and g2 and g3 and g4
    print("\nVERDICT:", flush=True)
    if ok:
        print("  No bug. The structural-ceiling / star-frozen-seam framework is the XY network by design", flush=True)
        print("  (g2 = 4/(N-1), correct). The canonical HEISENBERG star survivor is 4/N; the ZZ potential", flush=True)
        print("  is the entire difference. Both are frozen, both commutant only in the high-Q limit.", flush=True)
    else:
        print("  a gate FIRED -- diagnose, do not loosen.", flush=True)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
