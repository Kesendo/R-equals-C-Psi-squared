"""ADVERSARIAL re-derivation of F124, written from scratch (not a port of the proof's verifier).

Claims under attack:
  (A) chain, band-edge carrier:  ||M||_F^2 + lam_min(MM^T) = 2, all N.
  (B) ||M||_F^2 = 2 - E,  lam_min = E,  E = (4/(N+1)) sin^2(pi/(N+1)).
  (C) the lam_min eigenvector is the staggered bond wave (-1)^b.

Probes the proof never ran:
  (P1) NON-band-edge carrier (rank 1,2 from top). Is the staggered vector STILL an
       eigenvalue (E_k)?  Is it STILL the minimum?  Does the sum stay = 2 or drop BELOW?
       -> tests whether the band-edge choice is load-bearing for the "= 2" (proof attributes
          the result to a carrier-INDEPENDENT recurrence/conserved-Q argument).
  (P2) sign of the bond-Gram off-diagonals c_a c_{a+2}: positive only for the nodeless
       band-edge mode. This (unstated) positivity is what the 'oscillation theorem' step
       silently needs to call the staggered vector the MINIMUM.
  (P3) topology gate: even ring (=2), odd ring (lam_min>0, frustrated), star (N/2).
  (P4) the "two 2's": is the headline 2 the bulk coordination number z, or ||V_b||_F^2?
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=5, suppress=True, linewidth=140)


def graph(topo, N):
    if topo == "chain":
        edges = [(i, i + 1) for i in range(N - 1)]
    elif topo == "ring":
        edges = [(i, (i + 1) % N) for i in range(N)]
    elif topo == "star":
        edges = [(0, i) for i in range(1, N)]
    elif topo == "ladder":  # 2xL ladder, bulk coordination 3
        L = N // 2
        edges = []
        for r in range(L - 1):
            edges.append((r, r + 1))
            edges.append((L + r, L + r + 1))
        for r in range(L):
            edges.append((r, L + r))
    else:
        raise ValueError(topo)
    A = np.zeros((N, N))
    for (i, j) in edges:
        A[i, j] = A[j, i] = 1.0
    return A, edges


def build_M(A, edges, carrier_rank=0):
    N = A.shape[0]
    w, V = np.linalg.eigh(A)
    order = np.argsort(-w)            # descending: rank 0 = top = band edge
    psi1 = V[:, order[carrier_rank]]
    modes = V[:, order]
    M = np.zeros((len(edges), N))
    for bi, (i, j) in enumerate(edges):
        for ki in range(N):
            pk = modes[:, ki]
            M[bi, ki] = pk[i] * psi1[j] + pk[j] * psi1[i]
    return M, psi1, w[order]


def analyse(A, edges, carrier_rank=0):
    M, psi1, evals_adj = build_M(A, edges, carrier_rank)
    fro2 = (M ** 2).sum()
    G = M @ M.T
    ev = np.linalg.eigvalsh(G)
    lam_min = ev[0]
    nb = len(edges)
    stag = (-1.0) ** np.arange(nb)
    stag_eval = float(stag @ G @ stag / (stag @ stag))            # Rayleigh quotient
    resid = float(np.linalg.norm(G @ stag - stag_eval * stag))    # 0 iff staggered is an eigenvector
    # off-diagonal signs of G (c_a c_{a+2} for chain): are they all positive?
    offdiag = np.array([G[a, a + 1] for a in range(nb - 1)])
    return dict(fro2=fro2, lam_min=lam_min, stag_eval=stag_eval, resid=resid,
                sum=fro2 + lam_min, evals_G=ev, offdiag=offdiag, psi1=psi1)


def main():
    GATE = {"pass": True}

    def gate(name, cond):
        flag = "ok " if cond else "GATE-FIRE"
        if not cond:
            GATE["pass"] = False
        print(f"   [{flag}] {name}", flush=True)

    print("=" * 92, flush=True)
    print("(A,B,C) chain, BAND-EDGE carrier: sum=2, halves=2-E and E, staggered is the min eigvec",
          flush=True)
    print("=" * 92, flush=True)
    print(f"{'N':>3} {'||M||_F^2':>11} {'lam_min':>11} {'sum':>10} {'E (closed)':>12} "
          f"{'stag_resid':>11} {'off>0?':>7}", flush=True)
    for N in range(3, 15):
        A, edges = graph("chain", N)
        r = analyse(A, edges, 0)
        E = (4.0 / (N + 1)) * np.sin(np.pi / (N + 1)) ** 2
        allpos = bool((r["offdiag"] > 0).all())
        print(f"{N:>3} {r['fro2']:>11.7f} {r['lam_min']:>11.7f} {r['sum']:>10.6f} {E:>12.7f} "
              f"{r['resid']:>11.1e} {str(allpos):>7}", flush=True)
        gate(f"N={N}: sum=2", abs(r["sum"] - 2) < 1e-10)
        gate(f"N={N}: ||M||^2 = 2-E", abs(r["fro2"] - (2 - E)) < 1e-10)
        gate(f"N={N}: lam_min = E", abs(r["lam_min"] - E) < 1e-10)
        gate(f"N={N}: staggered IS the min (resid~0 AND stag==lam_min)",
             r["resid"] < 1e-9 and abs(r["stag_eval"] - r["lam_min"]) < 1e-9)
        gate(f"N={N}: band-edge off-diagonals c_a c_a+2 all > 0", allpos)

    print("\n" + "=" * 92, flush=True)
    print("(P1,P2) NON-band-edge carrier on the chain. Proof's recurrence/Q-argument is",
          flush=True)
    print("        carrier-INDEPENDENT, so it predicts staggered eigenvalue = E_k for ANY mode k.",
          flush=True)
    print("        Question: is staggered still the MINIMUM?  Does sum stay = 2?", flush=True)
    print("=" * 92, flush=True)
    print(f"{'N':>3} {'rank':>4} {'||M||^2':>10} {'lam_min':>10} {'stag_eval':>10} {'E_k':>10} "
          f"{'sum':>9} {'resid':>9} {'off all>0':>9}", flush=True)
    for N in (7, 9, 11):
        A, edges = graph("chain", N)
        for rank in (0, 1, 2):
            r = analyse(A, edges, rank)
            k = rank + 1  # band edge is k=1 (top); rank counts down from top
            kk = N - rank  # actual sine index for descending order: top mode is sin index N
            Ek = (4.0 / (N + 1)) * np.sin(kk * np.pi / (N + 1)) ** 2
            allpos = bool((r["offdiag"] > 0).all())
            print(f"{N:>3} {rank:>4} {r['fro2']:>10.6f} {r['lam_min']:>10.6f} "
                  f"{r['stag_eval']:>10.6f} {Ek:>10.6f} {r['sum']:>9.5f} {r['resid']:>9.1e} "
                  f"{str(allpos):>9}", flush=True)
        print(flush=True)
    print("   READING: if for rank>0  stag_eval==E_k (resid~0) but lam_min<stag_eval and sum<2,",
          flush=True)
    print("   then the band-edge choice IS load-bearing for the '=2' (proof under-states it):",
          flush=True)
    print("   the recurrence makes staggered an EIGENVALUE for any carrier; only the band-edge's",
          flush=True)
    print("   nodeless (Perron) positivity makes it the MINIMUM, hence the equality not <=.",
          flush=True)

    print("\n" + "=" * 92, flush=True)
    print("(P3,P4) topology gate + the 'two 2's': headline 2 = bulk coordination z, or ||V_b||^2?",
          flush=True)
    print("=" * 92, flush=True)
    print(f"{'topo':>7} {'N':>3} {'||M||^2':>10} {'lam_min':>10} {'sum':>9} {'bulk deg z':>11}",
          flush=True)
    for topo in ("chain", "ring", "ring", "star", "ladder"):
        Ns = (6, 7) if topo == "ring" else (6,)
        for N in Ns:
            A, edges = graph(topo, N)
            r = analyse(A, edges, 0)
            deg = A.sum(axis=1)
            # crude bulk coordination = modal degree (most common interior degree)
            vals, cnts = np.unique(deg, return_counts=True)
            zbulk = vals[np.argmax(cnts)]
            par = "even" if (topo == "ring" and N % 2 == 0) else ("odd" if topo == "ring" else "")
            print(f"{topo+par:>7} {N:>3} {r['fro2']:>10.6f} {r['lam_min']:>10.6f} "
                  f"{r['sum']:>9.5f} {zbulk:>11.1f}", flush=True)
    print("\n   If 'sum' tracks bulk z (chain/ring->2, ladder->3) while ||V_b||^2=2 is fixed,",
          flush=True)
    print("   then the headline '2' is the coordination number z, not a universal constant.",
          flush=True)

    print("\n" + "=" * 92, flush=True)
    print(f"OVERALL GATE: {'ALL PASS (claims reproduce)' if GATE['pass'] else 'A GATE FIRED'}",
          flush=True)
    print("=" * 92, flush=True)


if __name__ == "__main__":
    main()
