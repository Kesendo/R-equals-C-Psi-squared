"""The drive-frame flight, predictions from below (experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md section 4).

Arm A of the October flight with rz(phi) as the first and rz(-phi) as the last operation on logical qubit 1 (physical
q5): a logical identity, a global turn of q5's drive frame. Through the patched replay of reality_window_crosstalk_fit
(the drive-crosstalk letter with TH tracking, the Stark letter) this prints, for phi in {0, pi/2, pi, 3pi/2}:
  (i)   no letter and the five Z letters: the odd-Y readings are invariant in phi up to the rounding of the two added
        unitaries (printed in units of machine epsilon; the law is algebraic, rz commutes with every Z letter);
  (ii)  the fitted crosstalk alone (stage 0, section 2 of the page) and all three families (section 6 of the fit):
        D(phi) = odd-Y(A(phi)) - odd-Y(A(0)) over the 36 odd-Y readings, its largest entry, and the noncentrality of the
        frame statistic T5 of section 5: the three D(phi) stacked into one 108-vector with the FULL covariance (the
        shared A(0) puts the covariance of A(0)'s parities into every cross block), multinomial covariances from the
        replay's own parities at 8192 shots per circuit, T5 = D' Sigma^+ D, with the 99 % point of chi^2 at the rank,
        and the same divided by the October drift factor 1.27;
  (iii) the share: S = <Y0 Z2> of the YZ basis for each phi and its mean over the four turns (the frame-blind part of
        the excess; a turning letter's first and second harmonics cancel in the mean of four turns spaced by pi/2).
The same for the frame turned on logical 0 (q4) and logical 2 (q6) at phi = pi (36 readings, one difference).
Output: simulations/results/reality_window/drive_frame_power.txt
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as CHI2, ncx2 as NCX2

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "simulations"))
sys.path.insert(0, str(REPO / "simulations" / "flight"))
import run_reality_window as R                       # noqa: E402
import reality_window_letters as rl                  # noqa: E402
import reality_window_crosstalk_fit as CF            # noqa: E402

XT_FIT = {("xtc", (5, 4)): 0.0079 * math.cos(-2.24), ("xts", (5, 4)): 0.0079 * math.sin(-2.24),
          ("xtc", (5, 6)): 0.0090 * math.cos(0.83), ("xts", (5, 6)): 0.0090 * math.sin(0.83),
          ("xtc", (4, 5)): 0.0106 * math.cos(-2.56), ("xts", (4, 5)): 0.0106 * math.sin(-2.56),
          ("xtc", (6, 5)): 0.0124 * math.cos(2.69), ("xts", (6, 5)): 0.0124 * math.sin(2.69)}
ALL_FIT = {("xtc", (5, 4)): 0.0064 * math.cos(-2.29), ("xts", (5, 4)): 0.0064 * math.sin(-2.29),
           ("xtc", (5, 6)): 0.0045 * math.cos(1.82), ("xts", (5, 6)): 0.0045 * math.sin(1.82),
           ("xtc", (4, 5)): 0.0018 * math.cos(2.84), ("xts", (4, 5)): 0.0018 * math.sin(2.84),
           ("xtc", (6, 5)): 0.0077 * math.cos(-2.57), ("xts", (6, 5)): 0.0077 * math.sin(-2.57),
           ("stx", (5, 4)): 1.0505, ("stx", (5, 6)): 1.1398, ("stx", (4, 5)): 0.3416, ("stx", (6, 5)): 0.4503,
           ("det", 4): -0.488, ("det", 5): 0.018, ("det", 6): -0.225, ("zz", (4, 5)): 0.011, ("zz", (5, 6)): -0.088}
FIVE = {("det", 4): 0.001, ("det", 5): 0.201, ("det", 6): 0.180, ("zz", (4, 5)): 0.002, ("zz", (5, 6)): -0.111}
SHOTS = 8192
KAPPA_OCT = 1.27
OUT = REPO / "simulations" / "results" / "reality_window" / "drive_frame_power.txt"
LINES = []


def say(s):
    print(s, flush=True)
    LINES.append(s)


def main():
    data = json.loads(CF.DATA.read_text(encoding="utf-8"))
    tot = {}
    for nm, (h1, h2) in data["counts"].items():
        t = dict(h1)
        for k, v in h2.items():
            t[k] = t.get(k, 0) + v
        tot[nm] = t
    A = R.readout_matrix_measured(tot)
    calib = R.calib_from(data["snapshot"])
    circs = R.load_circuits()
    names = sorted(c.name for c in circs if R.arm_of(c.name) == "A")
    by = {c.name: c for c in circs}
    ops = {nm: R.extract(by[nm]) for nm in names}
    idx = rl.reading_index(names)
    replay = CF.patched_replay()
    phys = R.phys_of("A")
    yz = [n for n in names if n.split("_")[1] == "YZ"][0]

    def turned(op_list, lq, phi):
        """rz(phi) first and rz(-phi) last on logical lq. The trailing rz is placed at the circuit's end, after an
        idle delay, so that it adds no time boundary to the replay's slicing (a virtual rz has no duration and its
        placement is immaterial on the device; a new boundary would re-slice the decay steps and expose the replay's
        own step error, 7e-6 on q4, instead of the law)."""
        if phi == 0.0:
            return op_list
        clock = [0.0, 0.0, 0.0]
        for name, qs, prm, dur in op_list:
            start = max(clock[q] for q in qs)
            for q in qs:
                clock[q] = start + dur
        t_end = max(clock)
        tail = [("delay", [lq], [], t_end - clock[lq])] if t_end > clock[lq] else []
        return [("rz", [lq], [phi], 0.0)] + list(op_list) + tail + [("rz", [lq], [-phi], 0.0)]

    def parities(L, lq, phi):
        out = {}
        for nm in names:
            rho = replay(turned(ops[nm], lq, phi), phys, [calib["t1_us"][q] for q in phys], [calib["t2_us"][q] for q in phys], L)
            out[nm] = rl.parities_probs(R.arm_readout(A, "A") @ rl.probs_rho(rho))
        return out

    def vec(par):
        return np.array([par[n][m] for n, m in idx])

    def block(par):
        covs = {nm: R.cov_from_parities(par[nm], SHOTS) for nm in names}
        n = len(idx)
        S = np.zeros((n, n))
        pos = {}
        for i, (nm, m) in enumerate(idx):
            pos.setdefault(nm, []).append((i, m))
        for nm, lst in pos.items():
            for i, m1 in lst:
                for j, m2 in lst:
                    S[i, j] = covs[nm][(m1, m2)]
        return S

    def stacked(par0, pars):
        """D stacked over the turns, and T = D' Sigma^+ D with the shared-A(0) cross blocks."""
        k = len(pars)
        n = len(idx)
        D = np.concatenate([vec(p) - vec(par0) for p in pars])
        S0 = block(par0)
        Sig = np.kron(np.ones((k, k)), S0)
        for i, p in enumerate(pars):
            Sig[i * n:(i + 1) * n, i * n:(i + 1) * n] += block(p)
        w, V = np.linalg.eigh(Sig)
        keep = w > 1e-9 * w.max()
        T = float(D @ (V[:, keep] / w[keep]) @ V[:, keep].T @ D)
        return D, T, int(keep.sum())

    phis = [0.0, math.pi / 2, math.pi, 3 * math.pi / 2]
    for label, L in (("no letter", None), ("five Z letters", FIVE), ("fitted crosstalk (8)", XT_FIT), ("all three families (17)", ALL_FIT)):
        say("== %s, frame turned on q5 (logical 1)" % label)
        pars = [parities(L, 1, phi) for phi in phis]
        v0 = vec(pars[0])
        for phi, p in zip(phis[1:], pars[1:]):
            d = vec(p) - v0
            if L is None or L is FIVE:
                say("  phi %.3f: rounding residual of the two added unitaries, max |D| = %.3e = %.1f eps" % (phi, float(np.abs(d).max()), float(np.abs(d).max()) / np.finfo(float).eps))
            else:
                say("  phi %.3f: largest |D| %.4f, rms %.4f" % (phi, float(np.abs(d).max()), float(np.sqrt(np.mean(d ** 2)))))
        D, T, rank = stacked(pars[0], pars[1:])
        pt = CHI2.ppf(0.99, rank)
        say("  T5 stacked over the three turns (full covariance, %d shots per circuit): noncentrality %.1f, rank %d, 99 %% point %.1f, expected T5 %.0f, power %.3f; with the October drift factor %.2f: noncentrality %.1f, expected %.0f, power %.3f"
            % (SHOTS, T, rank, pt, rank + T, float(NCX2.sf(pt, rank, T)) if T > 0 else 0.01, KAPPA_OCT, T / KAPPA_OCT, rank + T / KAPPA_OCT, float(NCX2.sf(pt, rank, T / KAPPA_OCT)) if T > 0 else 0.01))
        s = [float(p[yz][(1, 1, 0)]) for p in pars]
        say("  share: S = <Y0 Z2> at the four turns " + ", ".join("%+.4f" % x for x in s) + "; mean %+.4f; S(0) - mean %+.4f" % (float(np.mean(s)), s[0] - float(np.mean(s))))
    for lq, q in ((0, 4), (2, 6)):
        for label, L in (("five Z letters", FIVE), ("fitted crosstalk (8)", XT_FIT), ("all three families (17)", ALL_FIT)):
            p0 = parities(L, lq, 0.0)
            p1 = parities(L, lq, math.pi)
            d = vec(p1) - vec(p0)
            if L is FIVE:
                say("== %s, frame turned on q%d at pi: rounding residual max |D| = %.3e = %.1f eps" % (label, q, float(np.abs(d).max()), float(np.abs(d).max()) / np.finfo(float).eps))
                continue
            D, T, rank = stacked(p0, [p1])
            pt = CHI2.ppf(0.99, rank)
            say("== %s, frame turned on q%d at pi: largest |D| %.4f, rms %.4f; T noncentrality %.1f, rank %d, 99 %% point %.1f, power %.3f; with drift noncentrality %.1f, power %.3f"
                % (label, q, float(np.abs(d).max()), float(np.sqrt(np.mean(d ** 2))), T, rank, pt, float(NCX2.sf(pt, rank, T)), T / KAPPA_OCT, float(NCX2.sf(pt, rank, T / KAPPA_OCT))))
    # the FRAME-BLIND bound: the smallest eps, equal on (5 -> 4) and (5 -> 6) at the crosstalk-alone fit's phases, the
    # other two pairs at zero, at which the power of T5 at its 99 % point is 0.5 (with the October drift factor)
    say("== the bound: power of T5 against a uniform eps on (5 -> 4) and (5 -> 6) at the fitted phases, other pairs zero, drift 1.27")
    base_par = parities(None, 1, 0.0)
    for eps in (0.001, 0.002, 0.003, 0.005):
        L = {("xtc", (5, 4)): eps * math.cos(-2.24), ("xts", (5, 4)): eps * math.sin(-2.24),
             ("xtc", (5, 6)): eps * math.cos(0.83), ("xts", (5, 6)): eps * math.sin(0.83)}
        pars = [parities(L, 1, phi) for phi in phis]
        D, T, rank = stacked(pars[0], pars[1:])
        pt = CHI2.ppf(0.99, rank)
        say("   eps %.3f: noncentrality %.1f (with drift %.1f), power %.3f" % (eps, T, T / KAPPA_OCT, float(NCX2.sf(pt, rank, T / KAPPA_OCT))))
    OUT.write_text(chr(10).join(LINES) + chr(10), encoding="utf-8", newline=chr(10))


if __name__ == "__main__":
    main()
