"""Clean gate, in the (1,1) single-excitation sector: is the star survivor a [H,rho]=0 commutant?

RESOLVED 2026-06-20 (see star_survivor_heisenberg.py + THE_STAR_FROZEN_SEAM.md 'Model scope'): NO bug.
This gate is HEISENBERG (H_SE includes the ZZ potential V), so it reads the Heisenberg star survivor
darkness 4/N -- NOT the docs' 4/(N-1), which is the XY (hopping-only) value the typed StructuralCeiling/
StarFrozenSeam framework uses by design. The mismatch this gate found is a MODEL difference (the ZZ term),
not an error. Yes: the survivor is the [H,rho]=0 commutant, but only in the high-Q LIMIT (||[H,rho]|| ~ 1/Q),
frozen at every Q -- in both models.

N=5 is the MARGINAL case (g2 = 4/(N-1) = 1 = the floor): the survivor is degenerate with the band edge,
so a full-L eig grabs a non-Hermitian mix and ||[H,rho]|| is meaningless there (that is where the outer
agent got lost). N>=6 (g2 = 4/(N-1) < 1) separates cleanly. Here we build the single-excitation
Hamiltonian H_SE = 2J*A (hopping) + J*diag(V) (the ZZ potential), the (1,1) Liouvillian, and test the
slowest mode for [H_SE, rho] = 0.  The docs say: survivor = darkest commutant among DEGENERATE H-eigenstates
(the a=0 leaf manifold = the flat band), rate -2g*4/(N-1).
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def star_HSE(N, J):
    A = np.zeros((N, N))
    for l in range(1, N):
        A[0, l] = A[l, 0] = 1.0
    V = np.zeros(N)
    V[0] = -(N - 1)                      # hub ZZ potential
    for l in range(1, N):
        V[l] = N - 3                     # each leaf
    return 2 * J * A + J * np.diag(V), A


def gate(N, J, g):
    H, A = star_HSE(N, J)
    d = N
    Id = np.eye(d)
    Lh = -1j * (np.kron(Id, H) - np.kron(H.T, Id))     # -i[H_SE, .], col-stack vec
    Ld = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            Ld[i + d * j, i + d * j] = 0.0 if i == j else -4 * g   # |i><j| dephases at -4g (i!=j)
    L = Lh + Ld
    w, V_ = np.linalg.eig(L)
    re = w.real
    cand = [k for k in range(len(w)) if re[k] < -1e-9]
    slow = max(re[k] for k in cand)
    space = [k for k in cand if abs(re[k] - slow) < 1e-5 * max(1.0, abs(slow))]
    pred = -2 * g * 4 / (N - 1)

    best = None
    loc = None
    aw, aV = np.linalg.eigh(A)
    a0 = aV[:, np.abs(aw) < 1e-9]                       # the flat a=0 leaf manifold
    P0 = a0 @ a0.conj().T
    for k in space:
        raw = V_[:, k].reshape(d, d, order="F")
        for rho in ((raw + raw.conj().T) / 2, (raw - raw.conj().T) / 2j):
            nr = np.linalg.norm(rho)
            if nr < 1e-9:
                continue
            comm = np.linalg.norm(H @ rho - rho @ H) / nr
            if best is None or comm < best:
                best = comm
                loc = np.linalg.norm(P0 @ rho @ P0) / nr

    print(f"N={N}  g={g}  Q={J/g:.4g}   g2=4/(N-1)={4/(N-1):.3f}", flush=True)
    print(f"  slowest (1,1):  Re={slow:+.6f}   |Im|max={max(abs(w[k].imag) for k in space):.1e}   "
          f"degeneracy={len(space)}", flush=True)
    print(f"  predicted commutant rate -2g*4/(N-1) = {pred:+.6f}   "
          f"[match: {abs(slow-pred)<0.05*abs(pred) if pred else False}]", flush=True)
    print(f"  min ||[H_SE,rho]||/||rho|| over the slow space = {best:.3e}   "
          f"(COMMUTANT if ~0)", flush=True)
    print(f"  a=0 flat-band localization of that mode = {loc:.3f}   (1.0 = lives in the degenerate band)\n", flush=True)


def sweep(N, J):
    """Is the frozen survivor the [H,rho]=0 commutant at EVERY Q, or only in the high-Q limit?"""
    print(f"\n=== Q-sweep, N={N} (g2=4/(N-1)={4/(N-1):.3f}): does it CONVERGE to the commutant as Q grows? ===", flush=True)
    print(f"{'Q':>8} {'Re/(-2g*g2)':>12} {'|Im|':>10} {'||[H,rho]||/||rho|| (~J=1)':>26} {'a0-loc':>8}", flush=True)
    H, A = star_HSE(N, J)
    d = N
    Id = np.eye(d)
    aw, aV = np.linalg.eigh(A)
    P0 = aV[:, np.abs(aw) < 1e-9] @ aV[:, np.abs(aw) < 1e-9].conj().T
    for g in (1.0, 0.2, 0.05, 0.01, 0.002, 0.0005):
        Lh = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
        Ld = np.zeros((d * d, d * d), dtype=complex)
        for i in range(d):
            for j in range(d):
                Ld[i + d * j, i + d * j] = 0.0 if i == j else -4 * g
        w, V_ = np.linalg.eig(Lh + Ld)
        re = w.real
        cand = [k for k in range(len(w)) if re[k] < -1e-9]
        slow = max(re[k] for k in cand)
        k = max(cand, key=lambda k: re[k])
        raw = V_[:, k].reshape(d, d, order="F")
        rho = (raw + raw.conj().T) / 2
        if np.linalg.norm(rho) < 1e-9:
            rho = (raw - raw.conj().T) / 2j
        nr = np.linalg.norm(rho)
        comm = np.linalg.norm(H @ rho - rho @ H) / nr
        loc = np.linalg.norm(P0 @ rho @ P0) / nr
        pred = -2 * g * 4 / (N - 1)
        print(f"{J/g:>8.4g} {slow/pred:>12.4f} {abs(w[k].imag):>10.1e} {comm:>26.4f} {loc:>8.3f}", flush=True)


if __name__ == "__main__":
    print("=== star survivor: commutant ([H,rho]=0) or not?  (1,1) sector, high Q ===\n", flush=True)
    for N in (5, 6, 7, 8):
        gate(N, 1.0, 0.05)
    sweep(6, 1.0)
    sweep(7, 1.0)
