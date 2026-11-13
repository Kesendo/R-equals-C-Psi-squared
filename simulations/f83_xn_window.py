"""The X^3 window of the April F83 runs: what the readings say and what a fit can and cannot close.

data/ibm_f83_signature_april2026/README.md, section "The X^⊗3 window". truly (XX + YY) and YZ+ZY commute with
X^⊗3 term by term, Z-dephasing respects it and |+,−,+⟩ is an eigenstate, so every single-site reading
anticommuting with X^⊗3 (⟨Y⟩, ⟨Z⟩ on either end) vanishes in the model, for the continuous evolution and for
the Trotter circuit alike (every Trotter factor commutes with X^⊗3; the native gates one by one do not).

  (a) the end Bloch vectors of truly and YZ+ZY on Marrakesh 2026-04-30 (path [4, 5, 6]) and Kingston 2026-05-05
      (path [43, 56, 63]), each letter pooled over the three bases that read it, σ from the binomial count;
  (b) on the flown Marrakesh circuits (the replay of simulations/reality_window_april_letters.py, both
      schedules), the 18 readings these two classes forbid (16 by X^⊗3, two of YZ+ZY's by the parity of n_Y)
      against: ideal gates; a detuning of the
      two end qubits, inside its bound; all nineteen letters of experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md
      §5, each held inside its bound (more letters than readings, so a small chi2 there decides nothing).

    python simulations/f83_xn_window.py  >  simulations/results/f83_xn_window.txt   (a few minutes)
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "flight"))
import reality_window_april_letters as ap  # noqa: E402
import reality_window_letters as rl  # noqa: E402

RUNS = {"Marrakesh 2026-04-30": HERE.parent / "data" / "ibm_f83_signature_april2026" / "f83_signature_ibm_marrakesh_20260430_190035.json",
        "Kingston 2026-05-05": HERE.parent / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_kingston_20260505_102806.json"}
CLS = ("truly_unbroken", "pi2_even_nontruly")
BOUND = {"det": 0.1, "sxa": 0.01, "sxp": 0.02, "czp": 0.05, "czk": 0.05, "zz": 0.1, "stk": 0.05}


def end_bloch(counts, end):
    vec, err = [], []
    for P in "XYZ":
        s = n = 0
        for basis, cnt in counts.items():
            if basis[0 if end == 0 else 1] != P:
                continue
            for key, c in cnt.items():
                bit = int(key[1]) if end == 0 else int(key[0])   # q0 is the right character
                s += c * (1 - 2 * bit)
                n += c
        m = s / n
        vec.append(m)
        err.append(math.sqrt(max(1 - m * m, 1e-9) / n))
    return vec, err


print("(a) end Bloch vectors (<X>, <Y>, <Z>) ± σ; <Y> and <Z> are forbidden by X^3")
for run, f in RUNS.items():
    d = json.loads(f.read_text(encoding="utf-8"))
    print(f"  {run}, path {d['path']}")
    for cls in CLS:
        for end in (0, 2):
            v, e = end_bloch(d["counts"][cls], end)
            print(f"    {cls:18s} q{end}: X {v[0]:+.3f}  Y {v[1]:+.3f} ({v[1] / e[1]:+.1f}σ)  Z {v[2]:+.3f} ({v[2] / e[2]:+.1f}σ)")


def readings2(OPS, L):
    out = {}
    for cls in CLS:
        acc = {}
        for b0 in "XYZ":
            for b2 in "XYZ":
                p8 = rl.assignment_product([ap.RO[0], ap.RO[1], (0.0, 0.0)]) @ rl.probs_rho(rl.replay(OPS[(cls, b0, b2)], ap.PATH, ap.T1, ap.T2, L))
                par = rl.parities_probs(p8)
                acc.setdefault((b0, b2), []).append(par[(1, 1, 0)])
                acc.setdefault((b0, "I"), []).append(par[(1, 0, 0)])
                acc.setdefault(("I", b2), []).append(par[(0, 1, 0)])
        for k, v in acc.items():
            out[(cls,) + k] = float(np.mean(v))
    return out


print()
print("(b) the window alone on the flown Marrakesh circuits")
M = ap.measured()
FORB = {c: ap.forbidden(c) for c in CLS}
print(f"  forbidden readings: {sum(len(v) for v in FORB.values())} ({', '.join(f'{c} {len(FORB[c])}' for c in CLS)})")
for schedule in ("alap", "asap"):
    OPS = ap.flown_ops(schedule)

    def resid(L):
        R = readings2(OPS, L)
        return np.array([(M[(c,) + k][0] - R[(c,) + k]) / M[(c,) + k][1] for c in CLS for k in FORB[c]])

    def chi2(L):
        return float(np.sum(resid(L) ** 2))

    print(f"  schedule {schedule}: ideal gates chi2 {chi2({}):.1f}")
    for name, letters in (("end-qubit detuning", [("det", 4), ("det", 6)]), ("all nineteen letters", list(rl.LETTERS))):
        lo = np.array([-BOUND[l[0]] for l in letters])
        fit = least_squares(lambda x: resid(dict(zip(letters, x))), np.zeros(len(letters)), bounds=(lo, -lo), x_scale=0.01)
        L = dict(zip(letters, fit.x))
        at = sum(abs(abs(v) - BOUND[k[0]]) < 1e-3 * BOUND[k[0]] for k, v in L.items())
        R = readings2(OPS, L)
        zq2 = R[("pi2_even_nontruly", "I", "Z")]
        print(f"    {name} ({len(letters)} letters, each inside its bound): chi2 {chi2(L):.1f}, {at} at their bound; "
              f"YZ+ZY <Z2> model {zq2:+.3f} against {M[('pi2_even_nontruly', 'I', 'Z')][0]:+.3f}")
        if len(letters) == 2:
            print(f"      detuning q4 {L[('det', 4)]:+.4f}, q6 {L[('det', 6)]:+.4f} rad/us")
