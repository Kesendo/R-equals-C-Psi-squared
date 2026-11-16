"""A post-flight reading of the colouring X^⊗3 behind the X^⊗3 window page's section 5.

This is a reading without a verdict, written after the flight. None of its numbers is a pre-registered quantity of
experiments/X3_WINDOW_HARDWARE_PREDICTION.md or experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md. The pre-registered
reading is run_drive_frame.colouring (section 5: X^⊗3 read in situ, a scale against the replay). This script asks what
the joint distribution of the three X outcomes behind that reading looks like, and which replay reproduces it.

What it reads from the flight's counts file. It only reads, prints to stdout and writes nothing:

1. The pooled XX-basis outcomes of B's standing variants, in normalised signs U0 = X0, U1 = -X1, U2 = X2, so that the
   ideal |+,-,+> reads +1 on all three. It prints the means, the three pair covariances, the third cumulant k3, the
   rate at which all three fail together, L3 = -ln|<X0X1X2>| and L1 = -sum ln|<X_l>|. Each comes with two bootstrap
   sigmas. The first resamples the counts and the RO arm's columns at their shots, as the runner's colouring does.
   The second resamples the RO columns at nro / kappa_RO shots and so carries the RO arm's drift factor, which the
   runner's sigma does not.
2. Robustness: the same moments under five readout treatments (the joint 8×8 inverse as in the runner, the product of
   its marginals, each half of the RO arm, no correction), for each half of the B circuits and for each variant alone.
3. The covariance of the outer pair, given that the middle qubit read as failed or as not failed.
4. Replays of B50_XX through the flight's own replay (rl.replay via run_drive_frame.replay_fn):
   - the freeze snapshot;
   - local T1 and T2 at the day's coherence: T2 from arm Q's fit, T1 and T2 divided by 4, 7 and 10, and T1 / 15 with
     T2 held at or below 2 T1;
   - at the coherence divided by 7, coherent letters (static ZZ, cz phase errors, arm Q's detunings);
   - at the same coherence, quasi-static Z noise (common, outer only, anticorrelated, independent), as Gauss-Hermite
     averages over a static detuning.
5. Every row's excess against the replay with T1 and T2 divided by 7, a local decay of about the measured size. The
   divisor at which its L3 meets the measured one exactly is found by bisection and printed beside it.

Two cautions for whoever reads the output. First, the rate of joint failure against the product of the single failure
rates separates nothing, because errors between the cz gates correlate the three qubits in the replay too; its absolute
size does separate. Second, a mixture of two product distributions on three bits has seven parameters, as many as the
outcome distribution. A mixture reading of the excess (a rare, strong failure shared by the three) is therefore a
reading and not a test.

    python drive_frame_colouring_reading.py [RESULT.json]

The default RESULT is data/ibm_drive_frame_november2026/drive_frame_ibm_marrakesh_20261113_203627.json.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
import run_drive_frame as R  # noqa: E402  read-only use: Chip, compact, readout_matrix, ro_kappa, vq, replay_fn

DEFAULT = REPO / "data" / "ibm_drive_frame_november2026" / "drive_frame_ibm_marrakesh_20261113_203627.json"
N_BOOT = 400
SEED = 20261114
REF_DIV = 7
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(9)
GH_W = GH_W / GH_W.sum()

# outcome index k = 4 q0 + 2 q2 + q1 (the runner's); the ideal outcome q0 = 0, q2 = 0, q1 = 1 reads +1 on all three
U = np.array([(1 - 2 * (k >> 2 & 1), 2 * (k & 1) - 1, 1 - 2 * (k >> 1 & 1)) for k in range(8)], dtype=float)
ALL_FAIL = [k for k in range(8) if tuple(U[k]) == (-1.0, -1.0, -1.0)]
COLS = ["U0", "U1", "U2", "cov01", "cov12", "cov02", "k3", "allfail", "L3", "L1", "gap"]
CCOLS = ["w|U1+", "cov02|U1+", "w|U1-", "U0|U1-", "U2|U1-", "cov02|U1-"]


def readings(p):
    """Moments of the normalised outcomes and the outer pair's conditional covariances, from outcome probabilities p.
    Entries may be slightly negative after an inverse RO matrix; they are kept, as in the runner's colouring."""
    p = np.asarray(p, dtype=float)
    p = p / p.sum()
    mu = U.T @ p

    def E(*i):
        return float(p @ np.prod(U[:, list(i)], axis=1))

    cov = {f"cov{i}{j}": E(i, j) - mu[i] * mu[j] for i, j in ((0, 1), (1, 2), (0, 2))}
    k3 = E(0, 1, 2) - mu[0] * cov["cov12"] - mu[1] * cov["cov02"] - mu[2] * cov["cov01"] - mu[0] * mu[1] * mu[2]
    L3 = -math.log(abs(E(0, 1, 2)))
    L1 = -sum(math.log(abs(m)) for m in mu)
    r = {"U0": mu[0], "U1": mu[1], "U2": mu[2], **cov, "k3": k3, "allfail": float(p[ALL_FAIL].sum()),
         "L3": L3, "L1": L1, "gap": L1 - L3}
    for c, tag in ((1.0, "+"), (-1.0, "-")):
        sel = U[:, 1] == c
        w = float(p[sel].sum())
        q = p[sel] / w
        m0, m2 = float(q @ U[sel, 0]), float(q @ U[sel, 2])
        r[f"w|U1{tag}"] = w
        r[f"U0|U1{tag}"], r[f"U2|U1{tag}"] = m0, m2
        r[f"cov02|U1{tag}"] = float(q @ (U[sel, 0] * U[sel, 2])) - m0 * m2
    return r


def boot_sigma(pc, A, nro, kappa=1.0):
    """Sigma of every reading: the pooled counts resampled multinomially, the RO arm's columns at nro / kappa shots."""
    rng = np.random.default_rng(SEED)
    N = int(round(pc.sum()))
    pm = pc / pc.sum()
    nr = [max(1, int(round(nro[j] / kappa))) for j in range(8)]
    rows = []
    for _ in range(N_BOOT):
        pb = rng.multinomial(N, pm) / N
        Ab = np.stack([rng.multinomial(nr[j], A[:, j] / A[:, j].sum()) / nr[j] for j in range(8)], axis=1)
        try:
            rows.append(readings(np.linalg.inv(Ab) @ pb))
        except (np.linalg.LinAlgError, ValueError):
            continue
    return {k: float(np.std([r[k] for r in rows])) for k in rows[0]}, len(rows)


def product_of_marginals(A):
    """The RO matrix each qubit would have alone (its marginal, averaged over the other two's preparations), tensored."""
    mats = {}
    for bit in (4, 2, 1):
        M = np.zeros((2, 2))
        for bp in (0, 1):
            cols = [k for k in range(8) if bool(k & bit) == bool(bp)]
            for br in (0, 1):
                M[br, bp] = np.mean([sum(A[j, k] for j in range(8) if bool(j & bit) == bool(br)) for k in cols])
        mats[bit] = M
    P = np.zeros((8, 8))
    for j in range(8):
        for k in range(8):
            P[j, k] = np.prod([mats[b][int(bool(j & b)), int(bool(k & b))] for b in (4, 2, 1)])
    return P


def ro_half(comp, h):
    Ah = np.zeros((8, 8))
    for j in range(8):
        pcj, _, n = comp[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"]
        Ah[:, j] = pcj[h] / n[h]
    return Ah


def replay(chip, T1, T2, L=None):
    rho = R.replay_fn(L)(chip.ops["B50_XX"], R.ROLE, list(T1), list(T2), L)
    p = np.clip(R.rl.probs_rho(rho), 0, None)
    return p / p.sum()


def quasi_static(chip, T1, T2, pattern, s):
    """Average over a static Z detuning (rad/us) drawn from a normal of standard deviation s: one shared draw scaled by
    the pattern's coefficients per role, or an independent draw per role."""
    if pattern == "independent":
        P = np.zeros(8)
        for a, wa in zip(GH_X, GH_W):
            for b, wb in zip(GH_X, GH_W):
                for c, wc in zip(GH_X, GH_W):
                    L = {("det", R.ROLE[0]): s * a, ("det", R.ROLE[1]): s * b, ("det", R.ROLE[2]): s * c}
                    P += wa * wb * wc * replay(chip, T1, T2, L)
        return P
    return sum(w * replay(chip, T1, T2, {("det", R.ROLE[l]): c * s * x for l, c in enumerate(pattern)})
               for x, w in zip(GH_X, GH_W))


def show(label, r, cols=COLS, signed=True):
    s = "+" if signed else ""
    cells = [f"{r[k]:{s}.5f}" if k == "allfail" else f"{r[k]:{s}.4f}" for k in cols]
    print(f"  {label:46s} " + " ".join(f"{c:>9s}" for c in cells))


def header(cols=COLS):
    print(f"  {'':46s} " + " ".join(f"{k:>9s}" for k in cols))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    data = json.loads(path.read_text(encoding="utf-8"))
    chip = R.Chip(data["backend"])
    if data.get("qpy_sha256") != R.sha256(chip.qpy):
        raise SystemExit(f"result was taken with circuits {data.get('qpy_sha256')}, this file is {R.sha256(chip.qpy)}")
    if not chip.rider:
        raise SystemExit(f"{chip.name} carries no B arms")
    comp = R.compact(data["counts"], chip.nf)
    A, nro = R.readout_matrix(comp)
    k_ro = R.ro_kappa(comp)
    names = [f"{v}_XX" for v in chip.variants("B")]
    pc = sum(comp[n][0].sum(0) for n in names)
    pm = pc / pc.sum()
    Ai = np.linalg.inv(A)

    print(f"== {chip.name} path {chip.path}: the colouring X^⊗3 behind X3_WINDOW section 5, a post-flight reading")
    print(f"   {path.name}; B variants pooled: {', '.join(names)}; shots {int(pc.sum())}; kappa_RO {k_ro:.2f}")
    print("   normalised signs U0 = X0, U1 = -X1, U2 = X2 (ideal |+,-,+> reads +1 on all three)")

    meas = readings(Ai @ pm)
    sig, nb = boot_sigma(pc, A, nro)
    sigk, _ = boot_sigma(pc, A, nro, kappa=k_ro)
    print("\n-- 1. measured, joint 8×8 inverse as in the runner, with two bootstrap sigmas "
          f"({nb} draws, seed {SEED})")
    header()
    show("measured", meas)
    show("sigma: shots and RO columns, as the runner", sig, signed=False)
    show("sigma_k: RO columns at nro / kappa_RO", sigk, signed=False)

    print("\n-- 2. robustness of the measured moments")
    header()
    show("joint 8×8 inverse (the runner)", meas)
    show("product of the joint matrix's marginals", readings(np.linalg.inv(product_of_marginals(A)) @ pm))
    show("RO arm, first half only", readings(np.linalg.inv(ro_half(comp, 0)) @ pm))
    show("RO arm, second half only", readings(np.linalg.inv(ro_half(comp, 1)) @ pm))
    show("no RO correction (raw counts)", readings(pm))
    for h in (0, 1):
        pch = sum(comp[n][0][h] for n in names)
        show(f"B circuits, half {h + 1} (joint inverse)", readings(Ai @ (pch / pch.sum())))
    for n in names:
        pcn = comp[n][0].sum(0)
        show(f"variant {n} alone (joint inverse)", readings(Ai @ (pcn / pcn.sum())))
    p0q0 = sum(A[j, 4] for j in range(8) if not j & 4)
    p0q2 = sum(A[j, 2] for j in range(8) if not j & 2)
    both = sum(A[j, 6] for j in range(8) if not j & 4 and not j & 2)
    print(f"  RO arm: P(q{chip.path[0]} reads 0 | 1) {p0q0:.4f}, P(q{chip.path[2]} reads 0 | 1) {p0q2:.4f}; both read 0"
          f" when both prepared in 1: {both:.4f} against their product {p0q0 * p0q2:.4f}")

    # the day's coherence: arm Q's fit (as the runner's colouring forms T2 from it), and local T1, T2 scaled
    vq_out = R.vq(chip, comp, A, nro, {}, None, k_ro)
    T1s, T2s = chip.T1, chip.T2
    T2q = [1.0 / max(vq_out[str(l)]["g"], 1.0 / (2 * T1s[l])) for l in range(3)]
    dq = {("det", R.ROLE[l]): vq_out[str(l)]["delta"] for l in range(3)}
    T1r, T2r = [t / REF_DIV for t in T1s], [t / REF_DIV for t in T2s]
    ref = readings(replay(chip, T1r, T2r))

    lo, hi = 1.0, 40.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if readings(replay(chip, [t / mid for t in T1s], [t / mid for t in T2s]))["L3"] < meas["L3"]:
            lo = mid
        else:
            hi = mid
    f_star = 0.5 * (lo + hi)

    print("\n-- 3. the outer pair given the middle qubit (U1 = +1: q"
          f"{chip.path[1]} read as not failed; U1 = -1: read as failed)")
    header(CCOLS)
    show("measured (joint inverse)", meas, CCOLS)
    show("sigma, as the runner", sig, CCOLS, signed=False)
    show("sigma_k", sigk, CCOLS, signed=False)
    show("no RO correction (raw counts)", readings(pm), CCOLS)
    show("replay, freeze snapshot", readings(replay(chip, T1s, T2s)), CCOLS)
    show(f"replay, T1 and T2 / {REF_DIV}", ref, CCOLS)

    rows = [("replay, freeze snapshot", readings(replay(chip, T1s, T2s))),
            ("replay, T2 from arm Q (T2 = 1/g)", readings(replay(chip, T1s, T2q)))]
    for f in (4, 7, 10):
        rows.append((f"replay, T1 and T2 / {f}", readings(replay(chip, [t / f for t in T1s], [t / f for t in T2s]))))
    rows.append(("replay, T1 / 15, T2 at most 2 T1",
                 readings(replay(chip, [t / 15 for t in T1s], [min(t2, 2 * t1 / 15) for t1, t2 in zip(T1s, T2s)]))))
    letters = []
    for z in (0.03, 0.1, 0.3):
        letters.append((f"+ static ZZ, both bonds, {z} rad/us",
                        {("zz", (R.ROLE[0], R.ROLE[1])): z, ("zz", (R.ROLE[1], R.ROLE[2])): z}))
    for kz in (0.1, -0.1, 0.2):
        letters.append((f"+ cz phase error, both pairs, {kz:+} rad",
                        {("czk", (R.ROLE[0], R.ROLE[1])): kz, ("czk", (R.ROLE[1], R.ROLE[2])): kz}))
    letters.append(("+ arm Q's three detunings", dq))
    letters.append((f"+ q{chip.path[1]} detuned by 0.5 rad/us", {("det", R.ROLE[1]): 0.5}))
    for lab, L in letters:
        rows.append((lab, readings(replay(chip, T1r, T2r, L))))
    for lab, pattern, s in (("common", (1, 1, 1), 0.1), ("common", (1, 1, 1), 0.2), ("outer only", (1, 0, 1), 0.2),
                            ("a, c against b", (1, -1, 1), 0.2), ("independent", "independent", 0.2)):
        rows.append((f"+ quasi-static Z, {lab}, s {s} rad/us", readings(quasi_static(chip, T1r, T2r, pattern, s))))

    print(f"\n-- 4. replays of B50_XX; rows marked + are at T1 and T2 / {REF_DIV}")
    header()
    show("measured (joint inverse)", meas)
    for lab, r in rows:
        show(lab, r)

    print(f"\n-- 5. excess against the replay with T1 and T2 / {REF_DIV} (its L3 {ref['L3']:.4f} against the measured "
          f"{meas['L3']:.4f}; the divisor at which they meet exactly: {f_star:.2f}); each row minus that reference")
    header()
    dmeas = {k: meas[k] - ref[k] for k in COLS}
    show("measured", dmeas)
    show("  in sigma, as the runner", {k: dmeas[k] / sig[k] if sig[k] > 0 else float("nan") for k in COLS})
    show("  in sigma_k, RO drift carried", {k: dmeas[k] / sigk[k] if sigk[k] > 0 else float("nan") for k in COLS})
    for lab, r in rows:
        show(lab, {k: r[k] - ref[k] for k in COLS})
    print(f"   arm Q (the runner's fit, readings): delta {[round(vq_out[str(l)]['delta'], 4) for l in range(3)]} rad/us, "
          f"T2 = 1/g {[round(t, 1) for t in T2q]} us; snapshot T1 {[round(t, 1) for t in T1s]} us, "
          f"T2 {[round(t, 1) for t in T2s]} us")


if __name__ == "__main__":
    main()
