"""Gate-first verifier: the VALUE/VECTOR felt-time split, and its CROSSOVER at the handover.
felt_time_dimensions arc, step (A). 2026-06-19, Tom + Claude.

K_decay = gamma/|Re lambda| = 1/(2<n_XY>) is the survivor's felt (dose) lifetime. The painter's
alpha reads the eigenVECTOR rotation under a single delta-J bond defect. The clean "value frozen /
vector rotates" split needs Re lambda (-> K_decay) DEFECT-INVARIANT.

ROUND 1+2 (table mode below) FIRED for the interior survivor: Re lambda moves at first order
(Re(dlam) = -0.86/-0.29/-0.12/-0.054 for N=4/5/6/8, linear), because its darkness <n_XY> is SOFT
(fractional, hopping-dependent). The (0,1) band edge is RIGID-dark (<n_XY>=1 -> Re=-2g structural,
J-independent) -> there Re(dlam) ~ 0, K_decay frozen, the clean split holds.

CROSSOVER (main mode): the survivor switches interior (p,p) -> band edge (0,1) at the handover Q*.
CLAIM: K_decay's defect-sensitivity |Re(dlam)| switches OFF exactly at that sector switch -- i.e.
the value/vector split turns ON at the handover.
GATE (can fire): if |Re(dlam)| does NOT drop to ~0 when the survivor becomes (0,1), or the drop is
not co-located with the sector switch, the crossover claim is wrong.
"""
import sys
import numpy as np
from scipy.linalg import eig

sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds, H_p, survivor  # reuse the exact machinery


def block_LV(N, prow, pcol, J, g, bnds, defect_bond):
    """The (prow,pcol) block Liouvillian L and the single-bond dJ channel V = dL/d(dJ)."""
    Hr, sr = H_p(N, prow, J, bnds)
    Hc, sc = H_p(N, pcol, J, bnds)
    Hbr, _ = H_p(N, prow, 1.0, [defect_bond])
    Hbc, _ = H_p(N, pcol, 1.0, [defect_bond])
    dr, dc = len(sr), len(sc)
    Ir, Ic = np.eye(dr), np.eye(dc)

    def liou(Mr, Mc):
        return -1j * (np.kron(Mr, Ic) - np.kron(Ir, Mc.T))

    deph = np.array([-2.0 * g * bin(sr[a] ^ sc[b]).count("1")
                     for a in range(dr) for b in range(dc)])
    return liou(Hr, Hc) + np.diag(deph), liou(Hbr, Hbc)


def slowest_idx(w):
    re = w.real
    c = np.where(re < -1e-7)[0]
    return int(c[np.argmax(re[c])])


def re_shift(N, prow, pcol, J, g, bnds, defect):
    """First-order Re(d lambda_s/d dJ) of the slowest (prow,pcol) mode = the K_decay defect-shift,
    with biorthogonality (1 = clean, ->0 = near EP)."""
    L, V = block_LV(N, prow, pcol, J, g, bnds, defect)
    w, vl, vr = eig(L, left=True, right=True)
    s = slowest_idx(w)
    l, r = vl[:, s], vr[:, s]
    dl = (l.conj() @ V @ r) / (l.conj() @ r)
    biorth = abs(l.conj() @ r) / (np.linalg.norm(l) * np.linalg.norm(r))
    return dl.real, biorth


def crossover(N, topo="chain", Qs=None):
    # Classify by <n_XY> (the physics: rigid darkness =1 vs soft <1), NOT by the (p,p) sector label,
    # which is filling-DEGENERATE on the open chain and hops arbitrarily among the tied (p,p) blocks.
    if Qs is None:
        Qs = [round(1.0 + 0.2 * k, 2) for k in range(16)]   # 1.0 .. 4.0
    J = 1.0
    bnds = bonds(N, topo)
    defect = bnds[len(bnds) // 3] if len(bnds) >= 3 else bnds[0]
    print(f"--- {topo} N={N}  (defect bond {defect}; the (p,p) label is filling-DEGENERATE and hops "
          f"-- read <n_XY>, not the label) ---", flush=True)
    print(f"{'Q':>5} {'<n_XY>':>7} {'|Re dlam|':>10} {'biorth':>6}  {'sector':>7}  regime", flush=True)
    soft, rigid, switch_Q, prev_rigid = [], [], None, None
    for Q in Qs:
        g = 1.0 / Q
        re, im, sec, nxy = survivor(N, J, g, bnds)
        dRe, bo = re_shift(N, sec[0], sec[1], J, g, bnds, defect)
        is_rigid = bool(nxy > 0.999)   # cast: numpy.bool_ breaks the `is False` switch check below
        reg = "RIGID band edge (n=1): K_decay FROZEN" if is_rigid else "SOFT survivor (n<1): K_decay moves"
        mark = ""
        if prev_rigid is False and is_rigid:
            switch_Q = Q; mark = "   <== HANDOVER (n_XY -> 1)"
        prev_rigid = is_rigid
        (rigid if is_rigid else soft).append(abs(dRe))
        print(f"{Q:>5.2f} {nxy:>7.4f} {abs(dRe):>10.4f} {bo:>6.3f}  {str(sec):>7}  {reg}{mark}", flush=True)
    ok = bool(soft and rigid and min(soft) > 0.02 and max(rigid) < 1e-6 and switch_Q)
    if soft and rigid:
        print(f"  => SOFT |Re dlam| in [{min(soft):.3f}, {max(soft):.3f}] (K_decay MOVES) | "
              f"RIGID max {max(rigid):.1e} (K_decay FROZEN) | handover Q~{switch_Q}  "
              f"[{'CONFIRMED' if ok else 'GATE FIRED'}]", flush=True)
    else:
        print("  => sweep did not bracket both regimes (widen Q)", flush=True)
    print(flush=True)
    return ok


def main():
    print("CROSSOVER at the handover: does K_decay's defect-sensitivity switch OFF when the survivor")
    print("crosses SOFT interior (p,p) -> RIGID band edge (0,1)?  |Re dlam| = K_decay's first-order")
    print("shift under a dJ defect (~0 frozen, >0 moves).\n", flush=True)
    allgood = True
    for N in (4, 5, 6):
        allgood &= crossover(N, "chain")
    print("VERDICT:", "CROSSOVER CONFIRMED -- the value/vector split turns ON at the handover Q*:"
          if allgood else "gate fired / inconclusive -- diagnose, do not loosen.")
    if allgood:
        print("  below Q* the survivor is soft-dark interior -> K_decay defect-sensitive (entangled with alpha);")
        print("  above Q* the survivor IS the rigid (0,1) band edge -> K_decay defect-invariant (clean split).")
    return allgood


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
