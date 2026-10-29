"""The reality-window flight read beside its verdicts (exploratory, not pre-registered).

Reads the committed flight record (ibm_marrakesh job dbhcohht9uhs73c9vpj0) through the frozen analysis
pieces of simulations/flight/run_reality_window.py and prints:
  - each fit arm's share of the null chi2 (its own block of the whitened residual) and its largest readings,
  - a linear joint fit of all nineteen letters,
  - nonlinear joint fits of the three detunings, and of the detunings with static ZZ,
  - the same five letters fitted arm by arm (a static letter is one number in every arm),
  - all nineteen fitted jointly and nonlinearly from the five-letter optimum (one start, a local fit),
  - the two halves of the job against each other, the even readings' contrast per arm, D0's single-site
    readings odd beside even, and A against D3 (one logical circuit, two compiles).
The verdicts themselves are run_reality_window.py --analyze; nothing here changes them.
Output: simulations/results/reality_window/post_flight.txt
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as CHI2

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "simulations" / "flight"))
import run_reality_window as R          # noqa: E402
import reality_window_letters as rl     # noqa: E402

DATA = REPO / "data" / "ibm_reality_window_october2026" / "reality_window_ibm_marrakesh_20261029_053203.json"
SUB = [("det", 4), ("det", 5), ("det", 6), ("zz", (4, 5)), ("zz", (5, 6))]


def main():
    data = json.loads(DATA.read_text(encoding="utf-8"))
    circs = R.load_circuits()
    by_name = {c.name: c for c in circs}
    calib = R.calib_from(data["snapshot"])
    tot = {}
    for nm, (h1, h2) in data["counts"].items():
        t = dict(h1)
        for k, v in h2.items():
            t[k] = t.get(k, 0) + v
        tot[nm] = t
    inflate = R.drift_inflation(data["counts"])[0]
    A = R.readout_matrix_measured(tot)
    par, nsh = {}, {}
    for nm, t in tot.items():
        if not nm.startswith("RO_"):
            par[nm] = rl.parities_counts(t)[0]
            nsh[nm] = sum(t.values())
    fit_names = sorted(n for n in par if R.arm_of(n) in R.FIT_ARMS)
    idx = rl.reading_index(fit_names)
    basepar = {nm: R.model_parities(by_name[nm], calib, A) for nm in fit_names}
    covs = {nm: R.cov_from_parities(basepar[nm], nsh[nm]) for nm in fit_names}
    # without the bootstrap term of the readout matrix: the verdict's chi2 carries it, these readings do not
    W, rank, Wh = R.weight_matrix(idx, covs, inflate, None)
    base = R.model_vector(circs, idx, calib, A)
    r = np.array([par[n][m] for n, m in idx]) - base
    print(f"chi2_null {float(r @ W @ r):.1f} over rank {rank} (99 % point {CHI2.ppf(0.99, rank):.1f}); "
          f"drift inflation {inflate:.3f}; no bootstrap term")

    S = np.linalg.pinv(W)
    blocks = {}
    for i, (nm, m) in enumerate(idx):
        blocks.setdefault(R.arm_of(nm), []).append(i)

    def whiten(ii):
        ev, U = np.linalg.eigh(S[np.ix_(ii, ii)])
        keep = ev > ev.max() * 1e-10
        return (U[:, keep] / np.sqrt(ev[keep])).T, int(keep.sum())

    print("\neach arm alone: chi2 of its 36 odd-Y readings against the baseline, and its largest readings")
    for arm in R.FIT_ARMS:
        ii = np.array(blocks[arm])
        Wa, k = whiten(ii)
        print(f"  {arm:3s} chi2 {float(np.sum((Wa @ r[ii]) ** 2)):7.1f} over {k}")
        sd = np.sqrt(np.diag(S)[ii])
        for j in np.argsort(-np.abs(r[ii] / sd))[:3]:
            nm, m = idx[ii[j]]
            print(f"        {nm:7s} mask (q0,q2,q1) {m}  residual {r[ii[j]]:+.4f}  z {r[ii[j]] / sd[j]:+6.1f}")

    F = R.fingerprints(circs, idx, calib, A)
    keys = list(rl.LETTERS)
    M = np.array([F[k] for k in keys]).T
    a = np.linalg.lstsq(M.T @ W @ M, M.T @ W @ r, rcond=None)[0]
    res = r - M @ a
    print(f"\nlinear joint fit of all nineteen letters: chi2 {float(res @ W @ res):.1f} "
          f"(99 % point at rank - 19: {CHI2.ppf(0.99, rank - 19):.1f}); the largest in units of their bounds:")
    for k, v in sorted(zip(keys, a), key=lambda kv: -abs(kv[1] / R.PLAUSIBLE[kv[0][0]]))[:6]:
        print(f"  {str(k):22s} {v:+.4f}  ({v / R.PLAUSIBLE[k[0]]:+.1f} x bound)")

    def fit(ii, sub, starts):
        ii = np.array(ii)
        Wa, k = whiten(ii)
        idx_a = [idx[i] for i in ii]
        fun = lambda x: Wa @ (r[ii] - (R.model_vector(circs, idx_a, calib, A, dict(zip(sub, x))) - base[ii]))
        best = None
        for st in starts:
            sol = least_squares(fun, st, x_scale=[R.PLAUSIBLE[s[0]] for s in sub], max_nfev=200)
            c = float(sol.fun @ sol.fun)
            if best is None or c < best[1]:
                best = (sol.x, c)
        return best[0], best[1], k

    every = list(range(len(idx)))
    print("\nnonlinear joint fits over all fit arms")
    for sub in (SUB[:3], SUB):
        x, c, k = fit(every, sub, ([0.0] * len(sub), [0.0, 0.2, 0.18, 0.0, -0.1][:len(sub)]))
        print(f"  {len(sub)} letters: chi2 {c:.1f} (99 % point {CHI2.ppf(0.99, k - len(sub)):.1f})  "
              + "  ".join(f"{s[0]}{s[1]} {v:+.3f}" for s, v in zip(sub, x)))

        if len(sub) == 5:
            x5 = x
    # all nineteen, nonlinear, one start: the five-letter optimum and zero elsewhere (a local fit)
    st = [dict(zip(SUB, x5)).get(k, 0.0) for k in keys]
    x, c, k = fit(every, keys, (st,))
    print(f"  19 letters (one start, local): chi2 {c:.1f} (99 % point {CHI2.ppf(0.99, k - 19):.1f})  largest: "
          + "  ".join(f"{s[0]}{s[1:]} {v:+.3f}" for s, v in sorted(zip(keys, x), key=lambda kv: -abs(kv[1] / R.PLAUSIBLE[kv[0][0]]))[:4]))

    print("\nthe five letters fitted arm by arm (rad/us; D0 has no leverage on them)")
    for arm in R.FIT_ARMS:
        x, c, k = fit(blocks[arm], SUB, ([0.0] * 5, [0.0, 0.2, 0.18, 0.0, -0.1], [0.1, -0.2, -0.2, 0.0, 0.0]))
        print(f"  {arm:3s} chi2 {c:6.1f} (99 % point {CHI2.ppf(0.99, k - 5):.0f})  "
              + "  ".join(f"{s[0]}{s[1]} {v:+.3f}" for s, v in zip(SUB, x)))

    # sx amplitude, the one letter family D0 sees: on D0 alone, and beside the five letters over all arms
    SXA = [("sxa", 4), ("sxa", 5), ("sxa", 6)]
    x, c, k = fit(blocks["D0"], SXA, ([0.0] * 3,))
    print(f"\nsx amplitude on D0 alone: chi2 {c:.1f} (99 % point {CHI2.ppf(0.99, k - 3):.1f})  "
          + "  ".join(f"sxa{s[1]} {v:+.4f}" for s, v in zip(SXA, x)))
    x, c, k = fit(every, SUB + SXA, (list(x5) + [0.0] * 3,))
    print(f"the five letters and sx amplitude over all arms: chi2 {c:.1f} (99 % point {CHI2.ppf(0.99, k - 8):.1f})  "
          + "  ".join(f"{s[0]}{s[1]} {v:+.4f}" for s, v in zip(SUB + SXA, x)))

    # the readout as the job's own RO arm measured it, against the snapshot
    print("\nreadout per qubit, marginal flip rates (RO arm, averaged over the other prepared bits) against the snapshot")
    bit = {0: 4, 2: 2, 1: 1}                       # probability index 4 q0 + 2 q2 + q1
    phys = {0: 4, 2: 6, 1: 5}
    for q in (0, 2, 1):
        f = {0: [], 1: []}
        for j in range(8):
            prep = 1 if j & bit[q] else 0
            f[prep].append(sum(A[i, j] for i in range(8) if (1 if i & bit[q] else 0) != prep))
        print(f"  q{phys[q]}: P(1|0) {np.mean(f[0]):.4f} (snapshot {calib['p01'][phys[q]]:.4f}), "
              f"P(0|1) {np.mean(f[1]):.4f} (snapshot {calib['p10'][phys[q]]:.4f})")
    zz_ro = {q: 0.0 for q in (0, 2, 1)}
    for j in range(8):
        nm = f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"
        hs = [rl.counts_to_probs(data["counts"][nm][h]) for h in (0, 1)]
        ns = [sum(data["counts"][nm][h].values()) for h in (0, 1)]
        for q in (0, 2, 1):
            prep = 1 if j & bit[q] else 0
            p = [sum(hs[h][i] for i in range(8) if (1 if i & bit[q] else 0) != prep) for h in (0, 1)]
            pm = (p[0] * ns[0] + p[1] * ns[1]) / sum(ns)
            zz_ro[q] += (p[0] - p[1]) ** 2 / (pm * (1 - pm) * (1 / ns[0] + 1 / ns[1]))
    print("  flip rate, half against half over the eight RO circuits: "
          + "  ".join(f"q{phys[q]} sum z^2 {zz_ro[q]:.1f}" for q in (0, 2, 1)) + " (8 dof each)")
    # the readout matrix's own counting noise on the baseline of S in arm A (V1 leaves it out, section 5)
    rng = np.random.default_rng(R.SHUFFLE_SEED)
    sb = []
    for _ in range(R.N_BOOT):
        A_b = np.zeros((8, 8))
        for j in range(8):
            t = tot[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"]
            n_ro = sum(t.values())
            A_b[:, j] = rng.multinomial(n_ro, rl.counts_to_probs(t)) / n_ro
        sb.append(R.model_parities(by_name["A_YZ"], calib, A_b)[(1, 1, 0)])
    sd_b = float(np.std(sb, ddof=1))
    sig = R.SIGMA * np.sqrt(inflate)
    s_a = par["A_YZ"][(1, 1, 0)] - basepar["A_YZ"][(1, 1, 0)]
    print(f"  S baseline in A {basepar['A_YZ'][(1, 1, 0)]:+.4f}, its readout-bootstrap sd {sd_b:.4f}; "
          f"s_A {s_a:+.4f} = {s_a / np.sqrt(sig ** 2 + sd_b ** 2):.1f} sigma with it, {s_a / sig:.1f} without")
    print(f"  the parse gate's readings on the day's baseline (measured readout): <X0Z2> {basepar['A_XZ'][(1, 1, 0)]:+.3f}, "
          f"<Z0X2> {basepar['A_ZX'][(1, 1, 0)]:+.3f}; measured {par['A_XZ'][(1, 1, 0)]:+.3f}, {par['A_ZX'][(1, 1, 0)]:+.3f}")

    # the two halves, interleaved by the shuffle: is the excess the same in both?
    z = {0: [], 1: []}
    for h in (0, 1):
        for nm, m in idx:
            cnt = data["counts"][nm][h]
            n = sum(cnt.values())
            b = basepar[nm][m]
            z[h].append((rl.parities_counts(cnt)[0][m] - b) / np.sqrt(max(1 - b * b, 1e-9) / n))
    z0, z1 = np.array(z[0]), np.array(z[1])
    print(f"\nthe two halves: correlation of their odd-Y z-scores {np.corrcoef(z0, z1)[0, 1]:.3f}, "
          f"diagonal chi2 {float(z0 @ z0):.0f} and {float(z1 @ z1):.0f}")

    # the real sector: measured even readings against the baseline's, slope through the origin per arm
    print("\nthe even (real-sector) readings: measured / baseline, least-squares slope through the origin")
    for arm in R.FIT_ARMS:
        mb = [(par[nm][m], basepar[nm][m]) for nm in fit_names if R.arm_of(nm) == arm
              for m in rl.MASKS if not rl.odd_y(tuple(nm.split("_")[1]), m)]
        mv, bv = np.array(mb).T
        print(f"  {arm:3s} {float(mv @ bv / (bv @ bv)):.3f}")

    # D0: the single-site reading of q0 and q2 after one pulse, odd (Y basis) beside even (Z basis)
    print("\nD0, single-site residuals after one pulse: odd <Y> beside even <Z> of the same qubit")
    for q, mask in ((0, (1, 0, 0)), (2, (0, 1, 0))):
        for letter in ("Y", "Z"):
            vals = [par[nm][mask] - basepar[nm][mask] for nm in fit_names
                    if R.arm_of(nm) == "D0" and nm.split("_")[1][0 if q == 0 else 1] == letter]
            print(f"  q{q} <{letter}>: " + "  ".join(f"{v:+.3f}" for v in vals))

    # the same logical circuit, compiled twice: A (April's gates) against D3 (compiled 2026-10-28)
    ra = np.array([r[i] for i in blocks["A"]])
    rd = np.array([r[i] for i in blocks["D3"]])
    assert [idx[i][1] for i in blocks["A"]] == [idx[i][1] for i in blocks["D3"]]
    print(f"\nA against D3, reading by reading: correlation of the 36 odd-Y residuals {np.corrcoef(ra, rd)[0, 1]:.2f}")


if __name__ == "__main__":
    main()
