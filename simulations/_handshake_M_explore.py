"""Just look. Build the handshake location matrix M[b,k] = <psi_k| V_b |psi_1>, self-check it is the
real object (the K-partner null: the k=N column must vanish), then throw every sum / checksum at it
and read off what -- if anything -- is clean, conserved, integer, or otherwise checksum-like.
No target. Tom: "ich ziele auf nichts, lass uns schauen was sich ergibt."
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=4, suppress=True, linewidth=120)


def psi(k, N):                       # bonding mode k=1..N on sites i=0..N-1 (open chain, Dirichlet)
    i = np.arange(N)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (i + 1) / (N + 1))


def build_M(N):
    """M[b,k] = <psi_k| V_b |psi_1>, V_b = |a><b|+|b><a| (bond b=(a,a+1)), carrier psi_1."""
    p1 = psi(1, N)
    bonds = [(a, a + 1) for a in range(N - 1)]
    P = np.array([psi(k, N) for k in range(1, N + 1)])      # P[k-1] = psi_k
    M = np.zeros((len(bonds), N))
    for bi, (a, b) in enumerate(bonds):
        for ki in range(N):
            M[bi, ki] = P[ki][a] * p1[b] + P[ki][b] * p1[a]
    return M, bonds


def look(N):
    M, bonds = build_M(N)
    nb = M.shape[0]
    print(f"\n========================  N={N}   (M is {nb}×{N}, bonds×modes)  ========================", flush=True)
    if N <= 6:
        print("M[b,k]  (rows=bonds b, cols=modes k=1..N):", flush=True)
        print(M, flush=True)

    # Stage 0 -- self-check: is this the REAL object?  K-partner null = the k=N column must be ~0.
    kN = np.linalg.norm(M[:, -1])
    print(f"\n[self-check] ‖k=N column‖ (K-partner null, must be ~0): {kN:.2e}   "
          f"{'OK -- real object' if kN < 1e-10 else 'FAIL -- not the handshake M'}", flush=True)
    rank = np.linalg.matrix_rank(M, tol=1e-9)
    print(f"[rank] = {rank}   (expected N−2 = {N-2})   {'OK' if rank == N-2 else '??'}", flush=True)

    # now just throw sums / checksums at it and read off
    sv = np.linalg.svd(M, compute_uv=False)
    print(f"\nsingular values: {sv}", flush=True)
    print(f"  Σσ  = {sv.sum():.5f}   Σσ² = {(sv**2).sum():.5f} = ‖M‖_F²", flush=True)
    print(f"row sums  Σ_k M[b,k]  (per bond): {M.sum(axis=1)}", flush=True)
    print(f"col sums  Σ_b M[b,k]  (per mode): {M.sum(axis=0)}", flush=True)
    print(f"total Σ_bk M[b,k]   = {M.sum():.6f}", flush=True)
    print(f"Σ_bk |M[b,k]|²      = {(M**2).sum():.6f}   (= Σσ²)", flush=True)
    print(f"Σ_bk |M[b,k]|       = {np.abs(M).sum():.6f}", flush=True)

    # Gram (the 'stacked lenses' inner products), both ways
    Gb = M @ M.T                       # bond×bond  (how bonds overlap as readings)
    Gk = M.T @ M                       # mode×mode
    print(f"trace(M Mᵀ) = trace(Mᵀ M) = {np.trace(Gb):.6f}", flush=True)
    print(f"eig(M Mᵀ) (bond Gram): {np.linalg.eigvalsh(Gb)}", flush=True)

    # candidate clean numbers vs N -- look for a pattern
    return dict(N=N, kN=kN, rank=rank, fro2=(M**2).sum(), total=M.sum(),
                abssum=np.abs(M).sum(), sigsum=sv.sum(), trace=np.trace(Gb),
                maxsv=sv[0])


def main():
    print("=== handshake M[b,k] = <psi_k|V_b|psi_1> : throw every sum/checksum at it, read off ===", flush=True)
    rows = [look(N) for N in (4, 5, 6, 7, 8)]
    print("\n\n================  scan vs N -- is anything clean / conserved / integer? ================", flush=True)
    print(f"{'N':>3} {'rank':>5} {'‖M‖_F²':>12} {'Σσ':>10} {'Σ|M|':>10} {'ΣM':>10} {'tr(MMᵀ)':>12} {'σ_max':>10}", flush=True)
    for r in rows:
        print(f"{r['N']:>3} {r['rank']:>5} {r['fro2']:>12.5f} {r['sigsum']:>10.5f} {r['abssum']:>10.5f} "
              f"{r['total']:>10.5f} {r['trace']:>12.5f} {r['maxsv']:>10.5f}", flush=True)
    # a couple of guesses for a checksum law
    print("\nguesses (do any land on a clean N-law?):", flush=True)
    for r in rows:
        N = r['N']
        print(f"  N={N}: ‖M‖_F²={r['fro2']:.5f}   ‖M‖_F²·(N+1)={r['fro2']*(N+1):.5f}   "
              f"‖M‖_F²/N={r['fro2']/N:.5f}   2·‖M‖_F²={2*r['fro2']:.5f}", flush=True)


if __name__ == "__main__":
    main()
