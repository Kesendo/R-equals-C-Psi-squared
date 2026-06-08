#!/usr/bin/env python3
"""Plot the palindrome classifier's protection landscape: the two coasts of the protected island, charted by
the C# trichotomy authority. See experiments/THE_PALINDROME_CLASSIFIER.md.

Reads simulations/results/tilt_sweep_csharp.tsv (columns: edge, angle_deg, verdict, pairing_error, rate_error).
That dump is produced by the C# engine: two parameter tilts driven through PauliPairTrichotomy.Classify at N=4
under Z-dephasing, recording the verdict and two breakings:

  * pairing_error: greedy max distance of each eigenvalue λ from its mirror partner -λ-2σ (real AND imaginary;
    this is the quantity SpectrumPairs thresholds at 1e-6 for the soft/hard verdict).
  * rate_error: the same mirror test on the REAL parts only (the decay rates), sorted. Comparing the two says
    whether a break lives in the rates (rate_error ≈ pairing_error) or in the coherent frequencies
    (rate_error ≈ 0 while pairing_error is large).

The two tilts:
  * field edge       H(θ) = cos(θ)·X + sin(θ)·Z   on every site   (k=1 templates, transverse → longitudinal)
  * frustration edge H(φ) = cos(φ)·XIX + sin(φ)·XXX + XXY + YXX    (k=3 templates, soft → frustrated hard)

To regenerate the dump, sweep those two template families through PauliPairTrichotomy.Classify (the k-body
overload accepts weighted PauliTerms), build L the same way (ChainKBody → PauliDephasingDissipator → Evd),
and dump the verdict plus the full and rate-only pairing errors per angle. This script only plots.
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

TSV = "simulations/results/tilt_sweep_csharp.tsv"
OUT = "simulations/results/two_edges.png"


def load():
    data = {"field": {"a": [], "e": [], "r": [], "v": []}, "frustration": {"a": [], "e": [], "r": [], "v": []}}
    with open(TSV, encoding="utf-8") as f:
        next(f)
        for line in f:
            edge, ang, verdict, err, rate = line.rstrip("\n").split("\t")
            data[edge]["a"].append(float(ang)); data[edge]["e"].append(float(err))
            data[edge]["r"].append(float(rate)); data[edge]["v"].append(verdict)
    for d in data.values():
        d["a"], d["e"], d["r"] = np.array(d["a"]), np.array(d["e"]), np.array(d["r"])
    return data


def first_hard(d):
    for ang, v in zip(d["a"], d["v"]):
        if v == "Hard":
            return ang
    return None


def main():
    d = load()
    fld, frs = d["field"], d["frustration"]
    print(f"field:       Truly→Soft→Hard, first Hard at θ={first_hard(fld):.3f}°, break(90°)={fld['e'][-1]:.3f}, rate(90°)={fld['r'][-1]:.1e}", flush=True)
    print(f"frustration: Soft→Hard,        first Hard at φ={first_hard(frs):.3f}°, break(90°)={frs['e'][-1]:.3f}, rate(90°)={frs['r'][-1]:.3f}", flush=True)
    print(f"at the first step off the protected point (~0.006°):  field {fld['e'][1]:.1e}  vs  frustration {frs['e'][1]:.1e}", flush=True)
    print(f"frustration rate vs full (median ratio): {np.median(frs['r'][1:]/frs['e'][1:]):.2f}  (≈1 means the break lives in the rates)", flush=True)

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    fa, fe = fld["a"][fld["a"] > 0], fld["e"][fld["a"] > 0]
    qa, qe = frs["a"][frs["a"] > 0], frs["e"][frs["a"] > 0]
    ax[0].loglog(fa, fe, "-o", ms=3, color="steelblue", label="field X→Z  (quadratic moat)")
    ax[0].loglog(qa, qe, "-o", ms=3, color="darkorange", label="frustration →XXX  (instant cliff)")
    ax[0].axhline(1e-6, ls="--", c="gray", lw=1, label="soft/hard threshold (1e-6)")
    ax[0].set_xlabel("tilt off the protected point (deg, log)")
    ax[0].set_ylabel("palindrome pairing error (log)")
    ax[0].set_title("onset: the field ramps up as θ²; frustration jumps at once")
    ax[0].grid(alpha=0.3, which="both"); ax[0].legend(fontsize=8)

    ax[1].plot(fld["a"], fld["e"], lw=2, color="steelblue", label="field: full break → 0.40 (deep)")
    ax[1].plot(fld["a"], fld["r"], lw=1.3, ls="--", color="steelblue", label="field: rate part → 0 (break is in the frequencies)")
    ax[1].plot(frs["a"], frs["e"], lw=2, color="darkorange", label="frustration: full break ~0.02 (shallow)")
    ax[1].plot(frs["a"], frs["r"], lw=1.3, ls="--", color="darkorange", label="frustration: rate part ≈ full (break is in the rates)")
    ax[1].set_xlabel("tilt off the protected point (deg)")
    ax[1].set_ylabel("palindrome breaking")
    ax[1].set_title("depth and where it lives: full break (solid) vs rate-only part (dashed)")
    ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8)

    fig.suptitle("Two coasts of the protection island (C# trichotomy): a gentle deep shore vs a sheer shallow cliff", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT, dpi=120, bbox_inches="tight")
    print(f"\nsaved {OUT}", flush=True)


if __name__ == "__main__":
    main()
