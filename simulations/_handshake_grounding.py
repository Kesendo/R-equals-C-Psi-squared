"""Grounding F124 from below: name the physical object M_loc IS, confirm by computing.

F124:  ‖M_loc‖_F² + λ_min(M_loc M_locᵀ) = 2  on the open chain / even ring,
       with ‖M_loc‖_F² = 2 − E,  λ_min = E,  E = carrier weight on the two free ends.

Candidate object (to confirm, not assert):
  (a) M_loc M_locᵀ = the GRAM MATRIX of {V_b ψ_1} (carrier displaced across each single bond).
      The frame operator S = M_locᵀ M_loc = Σ_b |V_b ψ_1⟩⟨V_b ψ_1|.  Is it a (tight) frame?
  (b) E = λ_min = literal Born probability of the single excitation sitting at a free end?
  (c) "2" = Σ_site deg·ψ_1²  (carrier seen once per incident bond)?  -- or is that 2−E?
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=5, suppress=True, linewidth=140)


def carrier(N):
    i = np.arange(N)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * (i + 1) / (N + 1))


def modes(N):
    """All single-excitation eigenmodes of the path adjacency, carrier (k=1) first."""
    return np.array([np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (np.arange(N) + 1) / (N + 1))
                     for k in range(1, N + 1)])


def Vb_psi(N, edges, psi):
    """The displaced carriers w_b = V_b ψ  (hopping across bond b), as columns in site space."""
    W = np.zeros((N, len(edges)))
    for bi, (i, j) in enumerate(edges):
        W[i, bi] += psi[j]
        W[j, bi] += psi[i]
    return W                                   # site x bond ; column bi = V_b ψ


def build_M(N, edges, psi):
    P = modes(N)                               # mode k (row) over sites
    M = np.zeros((len(edges), N))
    for bi, (i, j) in enumerate(edges):
        for ki in range(N):
            M[bi, ki] = P[ki][i] * psi[j] + P[ki][j] * psi[i]   # <psi_k|V_b|psi_1>
    return M


def edges_of(topo, N):
    if topo == "chain":
        return [(i, i + 1) for i in range(N - 1)]
    if topo == "ring":
        return [(i, (i + 1) % N) for i in range(N)]
    if topo == "star":
        return [(0, i) for i in range(1, N)]
    raise ValueError(topo)


def deg_vector(N, edges):
    d = np.zeros(N)
    for (i, j) in edges:
        d[i] += 1
        d[j] += 1
    return d


def main():
    print("=" * 100)
    print("(a)  Is M_loc M_locᵀ literally the GRAM matrix of the displaced carriers {V_b ψ_1}?")
    print("     Is the frame operator S = M_locᵀ M_loc = Σ_b |V_b ψ_1⟩⟨V_b ψ_1| a TIGHT frame?")
    print("=" * 100)
    for N in (5, 7, 9):
        edges = edges_of("chain", N)
        psi = carrier(N)
        M = build_M(N, edges, psi)
        W = Vb_psi(N, edges, psi)              # site x bond
        Gram_direct = W.T @ W                  # <V_b ψ | V_b' ψ>, built in SITE space
        Gram_fromM = M @ M.T                    # via the mode-coordinate matrix
        S = W @ W.T                            # frame operator, N x N, in site space
        S_fromM = M.T @ M
        gram_ok = np.allclose(Gram_direct, Gram_fromM, atol=1e-12)
        frame_ok = np.allclose(S, S_fromM, atol=1e-12)
        evS = np.linalg.eigvalsh(S)
        rankS = int(np.sum(evS > 1e-10))
        tight = np.ptp(evS[evS > 1e-10]) < 1e-9           # all nonzero eigenvalues equal?
        print(f"  N={N}:  M_locM_locᵀ == Gram⟨V_bψ|V_b'ψ⟩ (site-space)?  {gram_ok}"
              f"   |   M_locᵀM_loc == Σ_b|V_bψ⟩⟨V_bψ| ?  {frame_ok}")
        print(f"        frame-operator S spectrum = [{' '.join(f'{x:.4f}' for x in evS)}]")
        print(f"        rank(S) = {rankS}  (space dim N={N})  ->  deficient by {N - rankS}; "
              f"TIGHT frame (all nonzero eig equal)?  {tight}")
        # the kernel of S = the single-excitation direction NO single-bond displacement of the carrier reaches
        w0, V0 = np.linalg.eigh(S)
        ker = V0[:, 0]
        ker = ker / ker[np.argmax(np.abs(ker))]
        print(f"        kernel(S) direction (unreachable single-exc state)/max = "
              f"[{' '.join(f'{x:+.3f}' for x in ker)}]   (eig {w0[0]:.2e})")
    print()

    print("=" * 100)
    print("(c)  ‖M_loc‖_F² = Σ_site deg·ψ_1² ?   and is that '2' or '2 − E' on the chain?")
    print("=" * 100)
    print(f"{'N':>3} {'‖M‖_F²':>11} {'Σ deg·ψ²':>11} {'2·‖ψ‖²(=2)':>11} {'E=ends':>11}  reading")
    for N in range(4, 13):
        edges = edges_of("chain", N)
        psi = carrier(N)
        M = build_M(N, edges, psi)
        fro2 = (M ** 2).sum()
        deg = deg_vector(N, edges)
        degw = float((deg * psi ** 2).sum())
        E = psi[0] ** 2 + psi[-1] ** 2
        print(f"{N:>3} {fro2:>11.6f} {degw:>11.6f} {2*(psi**2).sum():>11.6f} {E:>11.6f}  "
              f"‖M‖²=Σdeg·ψ²={fro2:.4f}=2−E ;  '2'=2·‖ψ‖² (every site twice), chain SHORT by E")
    print()

    print("=" * 100)
    print("(b)  λ_min = E = the literal probability the single excitation sits on a FREE END?")
    print("     + the staggered eigenvector = the carrier's CURRENT (discrete derivative)?")
    print("=" * 100)
    print(f"{'N':>3} {'λ_min':>12} {'P(end)=c0²+cN²':>15} {'(4/(N+1))sin²θ':>15} {'sum=2?':>9}")
    for N in range(4, 13):
        edges = edges_of("chain", N)
        psi = carrier(N)
        M = build_M(N, edges, psi)
        G = M @ M.T
        lam = np.linalg.eigvalsh(G)[0]
        E = psi[0] ** 2 + psi[-1] ** 2
        Ecf = (4.0 / (N + 1)) * np.sin(np.pi / (N + 1)) ** 2
        s = (M ** 2).sum() + lam
        print(f"{N:>3} {lam:>12.8f} {E:>15.8f} {Ecf:>15.8f} {s:>9.6f}")

    print("\n  -- staggered combination  Σ_b (−1)^b V_b ψ_1  vs the carrier's discrete derivative --")
    for N in (7, 9):
        edges = edges_of("chain", N)
        psi = carrier(N)
        W = Vb_psi(N, edges, psi)
        stag = np.array([(-1.0) ** b for b in range(len(edges))])
        combo = W @ stag                                   # Σ_b (−1)^b V_b ψ_1   in site space
        # interior prediction (−1)^m (c_{m+1} − c_{m-1})  = staggered central derivative of the carrier
        pred = np.zeros(N)
        for m in range(N):
            cp = psi[m + 1] if m + 1 < N else 0.0
            cm = psi[m - 1] if m - 1 >= 0 else 0.0
            pred[m] = (-1.0) ** m * (cp - cm)
        G = M @ M.T if False else (build_M(N, edges, psi) @ build_M(N, edges, psi).T)
        lam = np.linalg.eigvalsh(G)[0]
        E = psi[0] ** 2 + psi[-1] ** 2
        print(f"  N={N}: Σ(−1)^b V_bψ = [{' '.join(f'{x:+.4f}' for x in combo)}]")
        print(f"        (−1)^m(c_(m+1)−c_(m-1)) = [{' '.join(f'{x:+.4f}' for x in pred)}]  match={np.allclose(combo, pred, atol=1e-12)}")
        print(f"        ‖Σ(−1)^b V_bψ‖² = {combo@combo:.6f}   (N−1)·E = {(N-1)*E:.6f}   "
              f"= vᵀGv with v=(−1)^b, λ_min·‖v‖² = {lam*(N-1):.6f}")
    print()

    print("=" * 100)
    print("TOPOLOGY: the '2' is the dispersive (degree-2) bulk value; boundary is the only free part")
    print("=" * 100)
    print(f"{'topo':>6} {'N':>3} {'‖M‖_F²':>11} {'λ_min':>11} {'sum':>11} {'E=Σ(2−deg)ψ²':>14}  staggered λ_min eigvec?")
    for topo in ("chain", "ring", "star"):
        for N in (6, 8):
            edges = edges_of(topo, N)
            A = np.zeros((N, N))
            for (i, j) in edges:
                A[i, j] = A[j, i] = 1.0
            w, V = np.linalg.eigh(A)
            order = np.argsort(-w)
            psi = V[:, order[0]]
            M = build_M(N, edges, psi) if topo != "ring" else None
            # build_M uses path 'modes'; for ring/star build M generally from adjacency modes:
            P = V[:, order].T                                  # mode rows, carrier first
            Mg = np.zeros((len(edges), N))
            for bi, (i, j) in enumerate(edges):
                for ki in range(N):
                    Mg[bi, ki] = P[ki][i] * psi[j] + P[ki][j] * psi[i]
            G = Mg @ Mg.T
            ev, EV = np.linalg.eigh(G)
            lam = ev[0]
            fro2 = (Mg ** 2).sum()
            deg = deg_vector(N, edges)
            bnd = float(((2 - deg) * psi ** 2).sum())
            vmin = EV[:, 0]
            vmin = vmin / vmin[np.argmax(np.abs(vmin))]
            stag = np.array([(-1.0) ** b for b in range(len(edges))])
            is_stag = abs(abs(np.dot(vmin, stag)) / (np.linalg.norm(vmin) * np.linalg.norm(stag)) - 1.0) < 1e-6
            print(f"{topo:>6} {N:>3} {fro2:>11.6f} {lam:>11.6f} {fro2+lam:>11.6f} {bnd:>14.6f}  "
                  f"{'YES' if is_stag else 'no '}")
        print()
    print("LAND:")
    print(" M_loc M_locᵀ = Gram of {V_b ψ_1} (carrier nudged across each one bond).")
    print(" ‖M‖_F² = Tr Gram = Σ deg·ψ_1² = 2 − E   (bulk site counted twice, once per incident bond).")
    print(" λ_min  = E = P(excitation on a free end), eigvec = staggered = the carrier's CURRENT.")
    print(" Invariant = 2·‖ψ‖² = 2: the dispersive-bulk coordination, boundary deficit restored by λ_min.")
    print(" Ring: E=0, λ_min=0 EXACT (staggered is a true zero mode). Star: breaks (hub degree>2).")
    print(" Family {V_bψ_1} is rank-(N−1): NOT a frame for the N-dim space, NOT tight. Strip 'tight frame'.")


if __name__ == "__main__":
    main()
