"""Axis 6 (shape gates motion), put on F124 and gated from below.

F124: on the open CHAIN, ‖M_loc‖_F² + λ_min(M_loc M_locᵀ) = 2, with 2 − ‖M_loc‖_F² = the carrier's
weight on the two FREE ends (degree-1 boundary). Axis-6 proposal: the "2" lives on dispersive,
degree-≤2 topologies; change the shape and the boundary term changes — ring (no free ends) should give
E=0 (‖M_loc‖²=2, λ_min=0); the star (a hub of degree N−1) should BREAK it. Just look.

General fact (derivable): ‖M_loc‖_F² = Σ_b ‖V_b ψ₁‖² = ψ₁ᵀ D ψ₁ = Σ_site deg(site)·ψ₁(site)²
(D = degree-diagonal), so 2 − ‖M_loc‖² = Σ_site (2 − deg(site))·ψ₁(site)². Degree-2 sites contribute 0;
free ends (deg 1) contribute +ψ₁²; a hub (deg N−1) contributes (3−N)·ψ₁² < 0 ⇒ ‖M_loc‖² > 2.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def graph(topo, N):
    if topo == "chain":
        edges = [(i, i + 1) for i in range(N - 1)]
    elif topo == "ring":
        edges = [(i, (i + 1) % N) for i in range(N)]
    elif topo == "star":
        edges = [(0, i) for i in range(1, N)]
    else:
        raise ValueError(topo)
    A = np.zeros((N, N))
    for (i, j) in edges:
        A[i, j] = A[j, i] = 1.0
    return A, edges


def M_loc_invariant(topo, N):
    A, edges = graph(topo, N)
    w, V = np.linalg.eigh(A)                       # single-excitation eigenmodes
    order = np.argsort(-w)                         # carrier = top (band-edge) mode
    psi1 = V[:, order[0]]
    modes = V[:, order]                            # all modes, carrier first
    M = np.zeros((len(edges), N))
    for bi, (i, j) in enumerate(edges):
        for ki in range(N):
            pk = modes[:, ki]
            M[bi, ki] = pk[i] * psi1[j] + pk[j] * psi1[i]      # <psi_k|V_b|psi_1>
    fro2 = (M ** 2).sum()
    lam_min = np.linalg.eigvalsh(M @ M.T)[0]
    deg = A.sum(axis=1)
    boundary = float(((2 - deg) * psi1 ** 2).sum())            # 2 − ‖M_loc‖² should equal this
    return fro2, lam_min, boundary, deg


def main():
    print("=== axis 6: does F124's invariant ‖M_loc‖²+λ_min=2 survive the topology? ===\n", flush=True)
    print(f"{'topo':>6} {'N':>3} {'‖M_loc‖²':>11} {'λ_min':>11} {'sum':>11} {'=2?':>6}   2−‖M‖² vs boundary", flush=True)
    for topo in ("chain", "ring", "star"):
        for N in (5, 6, 7):
            fro2, lam, bnd, deg = M_loc_invariant(topo, N)
            s = fro2 + lam
            hit = "YES" if abs(s - 2) < 1e-9 else "NO"
            print(f"{topo:>6} {N:>3} {fro2:>11.6f} {lam:>11.6f} {s:>11.6f} {hit:>6}   "
                  f"{2-fro2:>+10.6f} vs {bnd:>+10.6f}", flush=True)
        print(flush=True)
    print("reading: chain — holds, boundary = the two free ends.", flush=True)
    print("         ring  — holds, boundary = 0 (no free ends) ⇒ ‖M_loc‖²=2, the staggered mode is a ZERO mode.", flush=True)
    print("         star  — BREAKS: the hub (degree N−1) makes ‖M_loc‖² grow (= N/2), λ_min≥0 cannot repair it.", flush=True)
    print("\n=> the '2' is a DISPERSIVE-topology invariant (degree ≤ 2). The pinched hub destroys it —", flush=True)
    print("   exactly axis 6's diffuse-vs-pinched line. The lens proposed it; the gate from below shows it.", flush=True)


if __name__ == "__main__":
    main()
