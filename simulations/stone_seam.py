"""Drive the closure-structure SEAM (felt_time arc D), gate-first -- the seam can fail.
Seam hypothesis: the closure is ORTHOGONAL to <n_XY>; it reads STRUCTURE (defect-invariance + the
trajectory-warp geometry), not the dynamical rate magnitude. Two facets, one investigation:
  (a) GEOMETRIC value: closure ∝ dJ (linear response, so closure/dJ = Sum f is the geometric quantity),
      and it DEPENDS on the defect-bond position (a property of the mode's structure at that bond).
  (b) RIGIDITY ⊥ <n_XY>: on the RING the (p,p) are NOT filling-degenerate, so <n_XY> truly varies by
      sector -- the seam predicts ALL soft (p,p) break (OUT+coh) regardless of <n_XY>, only the
      structurally-pinned (0,1) (Re=-2gamma) holds; likewise the non-dispersive star.
Gates that can fire: closure NOT linear in dJ; closure position-INDEPENDENT; a soft ring (p,p) HOLDS,
or the break tracks <n_XY> instead of rigidity. N<=5 (full-L eig); flush per line."""
import sys
import importlib.util
import numpy as np
from scipy.linalg import eig

sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # Windows cp1252 console: keep the unicode
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds, survivor
import value_vector_felt_time as vv
_spec = importlib.util.spec_from_file_location("stone", "simulations/stone_survivor_alpha_closure.py")
st = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(st)


def readout(N, Q, pc, pr, topo="chain", dJ=0.02, defect_idx=None):
    J, g = 1.0, 1.0 / Q
    bnds = bonds(N, topo); d = 1 << N
    defect = bnds[len(bnds) // 3 if defect_idx is None else defect_idx]
    M, reMode = st.mode_embed(N, pc, pr, J, g, bnds)
    rho0 = st.herm_state(M, d)
    if rho0 is None:
        return None
    Vb = st.full_L(N, 1.0, 0.0, [defect]) - st.full_L(N, 0.0, 0.0, [defect])
    LA = st.full_L(N, J, g, bnds)
    paulis = ([st.site(st.X, i, N) for i in range(N)], [st.site(st.Y, i, N) for i in range(N)],
              [st.site(st.Z, i, N) for i in range(N)])
    tg = np.linspace(0.0, 5.0 / max(abs(reMode), 1e-3), 60)

    def P(L):
        w, R = eig(L); c = np.linalg.solve(R, rho0.flatten()); return st.purities(w, R, c, tg, paulis, N)

    pA, pB, pB2 = P(LA), P(LA + dJ * Vb), P(LA + (dJ / 2) * Vb)
    a1 = np.array([st.fit_alpha(pA[i], pB[i], tg) for i in range(N)])
    a2 = np.array([st.fit_alpha(pA[i], pB2[i], tg) for i in range(N)])
    f1, f2 = (a1 - 1) / dJ, (a2 - 1) / (dJ / 2)
    rel = (np.abs(f1) <= 10.0) & (np.abs(f1 - f2) <= 0.5 * (np.abs(f2) + 0.5))
    clo = float(np.sum(np.log(np.maximum(a1[rel], 1e-12)))) if rel.any() else float("nan")
    rf = f1[rel]; sa = float(np.sum(np.abs(rf)))
    coh = float(abs(rf.sum()) / sa) if sa > 0 else float("nan")
    return dict(clo=clo, coh=coh, sumf=float(rf.sum()), f1=f1, rel=rel, nxy=-reMode / (2 * g),
                defect=defect, win=("IN" if abs(clo) <= 0.05 else "OUT"))


def main():
    print("=== (a1) dJ-LINEARITY (N=5 chain (2,2), Q=1.5): is closure ∝ dJ? (=> closure/dJ = Sum f geometric) ===", flush=True)
    print(f"{'dJ':>6} {'Sum ln a':>9} {'clo/dJ':>8} {'coh':>5}", flush=True)
    base = None
    for dJ in (0.005, 0.01, 0.02, 0.04):
        r = readout(5, 1.5, 2, 2, "chain", dJ)
        print(f"{dJ:>6.3f} {r['clo']:>+9.4f} {r['clo']/dJ:>8.3f} {r['coh']:>5.2f}", flush=True)
        base = base or r['clo'] / dJ
    print(flush=True)

    print("=== (a2) DEFECT-POSITION (N=5 chain (2,2), Q=1.5, dJ=0.02): does closure depend on which bond? ===", flush=True)
    print(f"{'bond':>6} {'Sum ln a':>9} {'Sum f':>7} {'peak-f site':>12}", flush=True)
    for idx in range(len(bonds(5, "chain"))):
        r = readout(5, 1.5, 2, 2, "chain", 0.02, idx)
        peak = int(np.argmax(np.abs(np.where(r['rel'], r['f1'], 0.0))))
        print(f"{str(r['defect']):>6} {r['clo']:>+9.4f} {r['sumf']:>7.2f} {peak:>12}", flush=True)
    print(flush=True)

    print("=== (b) RIGIDITY ⊥ <n_XY> on RING + STAR (N=5, Q=1.5): do ALL soft (p,p) break, only (0,1) hold? ===", flush=True)
    print(f"{'topo':>6} {'sector':>7} {'<n_XY>':>7} {'Sum ln a':>9} {'win':>4} {'coh':>5}  reads", flush=True)
    for topo in ("ring", "star"):
        cands = [(p, p) for p in range(1, 5)] + [(0, 1)]
        for (pc, pr) in cands:
            r = readout(5, 1.5, pc, pr, topo, 0.02)
            if r is None:
                continue
            rigid = abs(r['nxy'] - 1.0) < 1e-6
            reads = ("structural-1 -> FROZEN" if r['win'] == "IN" else
                     f"soft (dressed <n_XY>={r['nxy']:.2f}) -> BREAKS")
            print(f"{topo:>6} ({pc},{pr}) {r['nxy']:>7.3f} {r['clo']:>+9.4f} {r['win']:>4} {r['coh']:>5.2f}  {reads}", flush=True)
        print(flush=True)


if __name__ == "__main__":
    main()
