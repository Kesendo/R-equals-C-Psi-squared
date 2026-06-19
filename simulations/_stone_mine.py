"""Mining the stone tool (felt_time arc B): feed N x Q x sector x topology, read the output numbers,
look for DERIVABLE STRUCTURE. Not adversarial re-validation (that was the two-lens review) -- broad
domain driving per use-the-tool-mine-the-data. Surfaces fields the witness computes but does not show
(per-site f, the rate <n_XY>, biorth). N<=5 (full-L eig); flushes per line."""
import sys
import importlib.util
import numpy as np
from scipy.linalg import eig

sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds
import value_vector_felt_time as vv
_spec = importlib.util.spec_from_file_location("stone", "simulations/_stone_survivor_alpha_closure.py")
st = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(st)


def readout(N, Q, pc, pr, topo="chain", dJ=0.02):
    """Full closure readout for one (N, Q, sector, topology): closure, sign-coherence, per-site alpha/f,
    the reliable mask, the mode rate <n_XY>, the first-order dRe and the prediction N*dJ*dRe/Re."""
    J, g = 1.0, 1.0 / Q
    bnds = bonds(N, topo); defect = bnds[len(bnds) // 3]; d = 1 << N
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
    rel = (np.abs(f1) <= 10.0) & (np.abs(f1 - f2) <= 0.5 * (np.abs(f2) + 0.5))   # canonical (relative) guard
    clo = float(np.sum(np.log(np.maximum(a1[rel], 1e-12)))) if rel.any() else float("nan")
    rf = f1[rel]; sa = float(np.sum(np.abs(rf)))
    coh = float(abs(rf.sum()) / sa) if sa > 0 else float("nan")
    dRe, bo = vv.re_shift(N, pc, pr, J, g, bnds, defect)
    return dict(clo=clo, coh=coh, a1=a1, f1=f1, rel=rel, nxy=-reMode / (2 * g), reMode=reMode,
                dRe=dRe, bo=bo, pred=N * dJ * dRe / reMode)


def main():
    dJ = 0.02
    print("=== 1. Q-SCALING of the survivor (2,2) closure (chain) — does it track N*dJ*dRe/Re? ===", flush=True)
    print(f"{'N':>2} {'Q':>5} {'<n_XY>':>7} {'Sum ln a':>9} {'coh':>5} {'pred':>8} {'meas/pred':>9} {'biorth':>6}", flush=True)
    for N in (4, 5):
        for Q in (1.0, 1.25, 1.5, 1.75, 2.0):
            r = readout(N, Q, 2, 2, "chain", dJ)
            ratio = r["clo"] / r["pred"] if abs(r["pred"]) > 1e-9 else float("nan")
            print(f"{N:>2} {Q:>5.2f} {r['nxy']:>7.3f} {r['clo']:>+9.3f} {r['coh']:>5.2f} {r['pred']:>+8.3f} "
                  f"{ratio:>9.2f} {r['bo']:>6.3f}", flush=True)
        print(flush=True)

    print("=== 2. SECTOR SCAN (N=5, Q=1.5, chain): which sectors break the closure? ===", flush=True)
    print(f"{'sector':>7} {'<n_XY>':>7} {'Sum ln a':>9} {'win':>4} {'coh':>5} {'rel':>5}  kind", flush=True)
    for (pc, pr) in [(1, 1), (2, 2), (3, 3), (4, 4), (0, 1)]:
        r = readout(5, 1.5, pc, pr, "chain", dJ)
        if r is None:
            continue
        win = "IN" if abs(r["clo"]) <= 0.05 else "OUT"
        kind = "frozen" if win == "IN" else ("rate-shift" if r["coh"] > 0.8 else "redistrib.")
        print(f"  ({pc},{pr}) {r['nxy']:>7.3f} {r['clo']:>+9.3f} {win:>4} {r['coh']:>5.2f} "
              f"{int(r['rel'].sum())}/5  {kind}", flush=True)
    print(flush=True)

    print("=== 3. PER-SITE alpha PROFILE (N=5, Q=1.5, survivor (2,2)): the standing-wave structure ===", flush=True)
    r = readout(5, 1.5, 2, 2, "chain", dJ)
    print(f"  per-site alpha = [{', '.join(f'{a:.4f}' for a in r['a1'])}]", flush=True)
    print(f"  per-site f     = [{', '.join(f'{f:+.3f}' for f in r['f1'])}]  (f_pred=dRe/Re={r['dRe']/r['reMode']:+.3f})", flush=True)
    print(f"  reliable       = {list(r['rel'].astype(int))}", flush=True)
    print(flush=True)

    print("=== 4. TOPOLOGY (N=5, Q=1.5): chain vs ring survivor closure ===", flush=True)
    print(f"{'topo':>6} {'sector':>7} {'<n_XY>':>7} {'Sum ln a':>9} {'win':>4} {'coh':>5}", flush=True)
    for topo, (pc, pr) in [("chain", (2, 2)), ("ring", (2, 2))]:
        r = readout(5, 1.5, pc, pr, topo, dJ)
        if r is None:
            continue
        win = "IN" if abs(r["clo"]) <= 0.05 else "OUT"
        print(f"{topo:>6} ({pc},{pr}) {r['nxy']:>7.3f} {r['clo']:>+9.3f} {win:>4} {r['coh']:>5.2f}", flush=True)


if __name__ == "__main__":
    main()
