"""felt_time_dimensions arc, D FOLLOW-UP -- the TRAJECTORY-level ground truth (2026-06-19). Gate-first.

This is the companion to _felt_time_amplitude_law.py. That (cheap, scalable) script establishes the
EXACT functional at the eigenvalue level: dRe(b) ~ (density-mode gradient at bond b)^2 (the diffusion
Rayleigh quotient -- "amplitude^2"). HERE we confirm that the actual full 4^N PTF painter closure
Sum f(b) = closure/dJ READS that rate shift at the bonds where it reads cleanly.

Step B's link (established + two-lens reviewed): Sum f(b) ~ N * (-dRe(b)) / |reS| at a bond where the
per-site f are sign-coherent (coh ~ 1, a genuine rate shift). The closure's reliability+sign-coherence
guard suppresses the LOW-gradient bonds (there the weak rate shift is swamped by redistribution, coh<0.8),
so the closure reads a clean rate shift ONLY at the HIGH-gradient bonds. Combined with dRe ~ grad^2, the
closure is therefore ~ grad^2 (amplitude^2) wherever it reads cleanly.

GATES that can fire:
  G0  seam shape: Sum f(b) ~0 at the chain-end bonds, max interior, mirror-symmetric (the seam).
  G1  B-link at the CLEAN bonds: at every coh>0.8 bond, the measured Sum f(b) must agree in SIGN and
      in O(1) magnitude with the eigenvalue prediction N*|dRe(b)|/|reS| (the review-pinned ~2.5x near-EP
      drift is allowed; a sign flip or order-of-magnitude miss is NOT).
  G2  the clean bonds ARE the high-gradient bonds: the coh>0.8 bonds must be the ones with the largest
      grad^2 (else the "closure reads the rate shift where the gradient is strong" story is wrong).
"""
import sys
import importlib.util
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds, survivor
from value_vector_felt_time import re_shift


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


st = _load("stone", "simulations/_stone_survivor_alpha_closure.py")
sm = _load("seam", "simulations/_stone_seam.py")


def density_gradient(N, prow, pcol, J, g, bnds):
    """grad^2(b) = (n(j)-n(j+1))^2 of the slow survivor density mode (see _felt_time_amplitude_law)."""
    M, _ = st.mode_embed(N, prow, pcol, J, g, bnds)
    diag = np.real(np.diag(M))
    n = np.zeros(N)
    for s in range(1 << N):
        if abs(diag[s]) > 0:
            for j in range(N):
                if (s >> j) & 1:
                    n[j] += diag[s]
    return np.array([(n[a] - n[c]) ** 2 for a, c in bnds])


def main():
    Q, dJ = 1.5, 0.02
    print("=== D: TRAJECTORY ground truth -- the full painter closure Sum f(b) reads dRe(b) ~ grad^2 ===")
    print(f"Q={Q}, dJ={dJ}, chain. Companion to _felt_time_amplitude_law.py (the block-level law).\n",
          flush=True)
    g0_ok = g1_ok = g2_ok = True
    for N in (4, 5):   # full 4^N closure; N=6 confirms (center bond Sum f=+11.05, coh=1.0) but re-eigs 4096^2/bond
        g = 1.0 / Q
        bnds = bonds(N, "chain")
        nb = len(bnds)
        reS, im_s, sec, nxy = survivor(N, 1.0, g, bnds)
        grad2 = density_gradient(N, sec[0], sec[1], 1.0, g, bnds)
        print(f"--- N={N} survivor {sec} <n_XY>={nxy:.3f}, reS={reS:+.3f} ---", flush=True)
        print(f"{'bond':>7} {'Sum f':>9} {'win':>4} {'coh':>5} | {'dRe':>8} {'grad^2':>8} "
              f"{'N|dRe|/|reS|':>12} {'ratio':>6}", flush=True)
        sumf, coh, dpred, clean, grads = [], [], [], [], []
        for idx in range(nb):
            b = bnds[idx]
            r = sm.readout(N, Q, sec[0], sec[1], "chain", dJ, idx)
            dRe = re_shift(N, sec[0], sec[1], 1.0, g, bnds, b)[0]
            pred = N * abs(dRe) / abs(reS)
            ratio = r["sumf"] / pred if pred > 1e-9 else float("nan")
            sumf.append(r["sumf"]); coh.append(r["coh"]); dpred.append(pred)
            clean.append(r["coh"] > 0.8); grads.append(grad2[idx])
            print(f"{str(b):>7} {r['sumf']:>+9.3f} {r['win']:>4} {r['coh']:>5.2f} | {dRe:>+8.4f} "
                  f"{grad2[idx]:>8.5f} {pred:>12.3f} {ratio:>6.2f}", flush=True)
        sumf = np.array(sumf); grads = np.array(grads); clean = np.array(clean)
        # G0 seam shape
        end_mag = max(abs(sumf[0]), abs(sumf[-1]))
        interior_max = np.max(np.abs(sumf[1:-1])) if nb > 2 else abs(sumf[0])
        mirror = float(np.max(np.abs(sumf - sumf[::-1]))) if nb > 1 else 0.0
        seam = (end_mag < 0.25 * interior_max) and (mirror < 0.1 * max(interior_max, 1e-9))
        # G1 B-link: at clean bonds, sign of Sum f matches +pred (faster decay -> alpha>1 -> +) and O(1)
        link = True
        for i in range(nb):
            if clean[i]:
                ratio = sumf[i] / dpred[i] if dpred[i] > 1e-9 else float("nan")
                if not (0.3 < ratio < 4.0):     # sign-positive and within the allowed near-EP drift band
                    link = False
        # G2 clean bonds are the high-grad bonds
        if clean.any() and (~clean).any():
            g2 = grads[clean].min() > grads[~clean].max()
        else:
            g2 = True
        g0_ok &= seam; g1_ok &= link; g2_ok &= g2
        print(f"   G0 seam end/interior={end_mag/max(interior_max,1e-9):.2f} mirror={mirror:.4f} "
              f"[{'ok' if seam else 'FIRED'}]  | G1 B-link at clean bonds [{'ok' if link else 'FIRED'}]  "
              f"| G2 clean=high-grad [{'ok' if g2 else 'FIRED'}]\n", flush=True)

    ok = g0_ok and g1_ok and g2_ok
    print("VERDICT:", "the trajectory closure CONFIRMS the block law: Sum f(b) reads a clean rate shift "
          "(coh~1)\n  only at the high-gradient bonds, there matching N|dRe|/|reS| in sign and O(1) "
          "magnitude.\n  With dRe ~ grad^2 (companion script), the closure functional is amplitude^2 "
          "in the\n  survivor density-mode gradient, wherever it reads cleanly."
          if ok else "a gate fired -- diagnose (G0 seam? G1 sign/magnitude? G2 clean-vs-grad?). Do NOT loosen.")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
