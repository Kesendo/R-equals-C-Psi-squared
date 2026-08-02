#!/usr/bin/env python3
"""The two coasts of the protected island, recomputed from below with an honest pairing metric.
See experiments/THE_PALINDROME_CLASSIFIER.md.

Why this script exists. The original sweep behind that experiment lives in
simulations/results/tilt_sweep_csharp.tsv, dumped by a C# runner that is no longer in the repo,
and its `pairing_error` column is a GREEDY match: each eigenvalue consumes the nearest still-free
partner, in spectrum order, with the choice forced to be an involution
(PauliPairTrichotomy.SpectrumPairs, whose own docstring says so). Greedy is an upper bound on the
real question ("is there ANY pairing of the spectrum with its mirror image, and how far is the
worst partner in the best one"), and it overshoots by more than a percent at 26 of the 190 field
angles and 56 of the 190 frustration angles. The overshoot is not cosmetic: it invents an apparent
peel-away on the field coast around 22.5 deg and moves the frustration coast's deepest break from
0.0353 at 38.0 deg to 0.1196 at 57.5 deg.

So this script recomputes both coasts from below and pairs OPTIMALLY (bottleneck assignment: the
smallest t for which a perfect matching exists using only partners within t). The sweep parameters
were reconstructed by reproducing the C# columns where greedy and optimal agree:

    N = 4, gamma = 0.05 (sigma = N*gamma = 0.2), Z-dephasing on every site, NO bond terms
    field coast        H(theta) = sum_l (cos theta X_l + sin theta Z_l)
    frustration coast  H(phi)   = sum_w (cos phi XIX + sin phi XXX + XXY + YXX) over 3-site windows

Three breakings are recorded per angle, all against the mirror lambda -> -lambda - 2 sigma:

    pairing_true  optimal (bottleneck) worst-partner distance over the full complex spectrum
    rate_error    the same test on the REAL parts alone (sorted; 1-D, so already optimal)
    freq_error    the same test on the IMAGINARY parts alone

freq_error is a control, not a measurement: L maps Hermitian operators to Hermitian ones, so its
spectrum is closed under conjugation, so the imaginary-part multiset is symmetric about 0
identically. It comes out at 1.3e-13 or below at every angle on both coasts, protected and broken
alike, which is why "the break moved into the frequencies" is not a statement this data can make.
rate_error is a relaxation of pairing_true, not an independent test: projecting a matching onto the
real axis cannot lengthen any edge, so rate_error <= pairing_true identically (measured: equality
of the maximum at 0.0 across the sweep).

The verdicts (Truly / Soft / Hard) are NOT recomputed here; they are the C# classifier's
(PauliPairTrichotomy.Classify), thresholded at 1e-6 on ITS greedy metric, and are carried over from
tilt_sweep_csharp.tsv by angle. The angle grid is carried over from the same file, which is why
`true_angle` below exists: that column is rounded to four decimals, and its small end is a radian
ladder (1e-4 to 5e-2 rad) converted to degrees, so nine of the angles are up to 0.5 percent off
their intended value if the rounded figure is taken at face value.

Outputs: simulations/results/two_coast_sweep.tsv and simulations/results/two_edges.png.
Runtime about a minute (380 eigendecompositions of a 256x256 Liouvillian plus 380 matchings).
"""
from __future__ import annotations
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from framework.diagnostics.f77_trichotomy import spectrum_pairing_error

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

N = 4
GAMMA = 0.05
SIGMA = N * GAMMA

SRC = "simulations/results/tilt_sweep_csharp.tsv"
TSV = "simulations/results/two_coast_sweep.tsv"
PNG = "simulations/results/two_edges.png"

I2 = np.eye(2, dtype=complex)
P = {"I": I2,
     "X": np.array([[0, 1], [1, 0]], dtype=complex),
     "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
     "Z": np.array([[1, 0], [0, -1]], dtype=complex)}


def kron_all(mats):
    out = np.array([[1.0 + 0j]])
    for m in mats:
        out = np.kron(out, m)
    return out


def placed(pattern, start, n=N):
    """Pauli string `pattern` placed at sites start..start+len(pattern)-1, identity elsewhere."""
    mats = [I2] * n
    for j, ch in enumerate(pattern):
        mats[start + j] = P[ch]
    return kron_all(mats)


def liouvillian(H, gamma=GAMMA, n=N):
    """L = -i[H, .] + gamma * sum_l (Z_l . Z_l - .), column-stacked."""
    d = 2 ** n
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(n):
        Zl = placed("Z", l, n)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


def H_field(theta_deg, n=N):
    t = np.deg2rad(theta_deg)
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for l in range(n):
        H += np.cos(t) * placed("X", l, n) + np.sin(t) * placed("Z", l, n)
    return H


def H_frustration(phi_deg, n=N):
    p = np.deg2rad(phi_deg)
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for s in range(n - 2):
        H += np.cos(p) * placed("XIX", s, n) + np.sin(p) * placed("XXX", s, n)
        H += placed("XXY", s, n) + placed("YXX", s, n)
    return H


def pairing_true(vals, sigma=SIGMA):
    """Bottleneck assignment distance between the spectrum and its mirror image.

    The cockpit already owns this: `framework.diagnostics.f77_trichotomy.spectrum_pairing_error`
    is the exact bottleneck (binary search on the distinct distances with a bipartite-matching
    feasibility test), and it builds its distance matrix as a SUM so the matrix is symmetric to
    the bit, which a subtracted form is not. This is a one-line adapter onto it."""
    return spectrum_pairing_error(vals, sigma)


def part_error(parts, mirrored):
    return float(np.max(np.abs(np.sort(parts) - np.sort(mirrored))))


def read_verdicts(path=SRC):
    v = {}
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            edge, ang, verdict, _pe, _re = line.rstrip("\n").split("\t")
            v[(edge, float(ang))] = verdict
    return v


def true_angle(rounded_deg):
    """Undo the source column's rounding where the angle came from a radian ladder.

    Nine of the grid's small angles are 1e-4, 2e-4, 5e-4 ... 5e-2 radians converted to degrees and
    then rounded to four decimals; taking the rounded figure at face value misses the intended
    angle by up to 0.5 percent, which shows up as a percent-level error in the break. No degree-grid
    point of this sweep collides with one of the nine, so the recovery is unambiguous."""
    for mantissa in (1, 2, 5):
        for exponent in range(-5, 0):
            radians = mantissa * 10.0 ** exponent
            if round(np.rad2deg(radians), 4) == rounded_deg:
                return np.rad2deg(radians)
    return rounded_deg


def sweep():
    verdicts = read_verdicts()
    angles = sorted({a for (e, a) in verdicts if e == "field"})
    rows = []
    for edge, Hf in (("field", H_field), ("frustration", H_frustration)):
        for ang in angles:
            vals = np.linalg.eigvals(liouvillian(Hf(true_angle(ang))))
            rows.append((edge, ang, verdicts.get((edge, ang), "?"),
                         pairing_true(vals),
                         part_error(vals.real, -vals.real - 2 * SIGMA),
                         part_error(vals.imag, -vals.imag)))
        print(f"  {edge}: {len(angles)} angles done", flush=True)
    return rows


def write_tsv(rows, path=TSV):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# two-coast protection landscape, recomputed from below by simulations/two_coast_sweep.py\n")
        f.write(f"# N={N} gamma={GAMMA} sigma={SIGMA} Z-dephasing on every site, no bond terms\n")
        f.write(f"# field:       H(theta) = sum_l (cos theta X_l + sin theta Z_l)\n")
        f.write(f"# frustration: H(phi) = sum_w (cos phi XIX + sin phi XXX + XXY + YXX), 3-site windows\n")
        f.write(f"# pairing_true = optimal (bottleneck) worst-partner distance for lambda -> -lambda-2sigma\n")
        f.write(f"# rate_error / freq_error = the same test on the real / imaginary parts alone\n")
        f.write(f"# verdict = the C# classifier's (greedy metric, threshold 1e-6), carried over from {SRC}\n")
        f.write(f"# angle_deg is that file's grid, rounded to 4 decimals; nine small angles are a radian\n")
        f.write(f"#   ladder (1e-4 to 5e-2 rad) and were evaluated at their exact value, not the rounded one\n")
        f.write(f"# the two angle-0 rows are the protected point: every column there is eigensolver noise\n")
        f.write("edge\tangle_deg\tverdict\tpairing_true\trate_error\tfreq_error\n")
        for edge, ang, verdict, pt, re_, fe in rows:
            f.write(f"{edge}\t{ang:.4f}\t{verdict}\t{pt:.6E}\t{re_:.6E}\t{fe:.6E}\n")
    print(f"wrote {path}")


def plot(rows, path=PNG):
    data = {}
    for edge, ang, verdict, pt, re_, fe in rows:
        d = data.setdefault(edge, {"a": [], "p": [], "r": [], "v": []})
        d["a"].append(ang); d["p"].append(pt); d["r"].append(re_); d["v"].append(verdict)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    titles = {"field": "the field coast: transverse X to longitudinal Z",
              "frustration": "the frustration coast: soft 3-body to frustrated XXX"}
    xlabels = {"field": "tilt angle θ (degrees)", "frustration": "frustration angle φ (degrees)"}
    for ax, edge in zip(axes, ("field", "frustration")):
        d = data[edge]
        a = np.array(d["a"]); p = np.array(d["p"]); r = np.array(d["r"])
        ax.semilogy(a, np.maximum(p, 1e-16), lw=1.8, label="mirror break (optimal pairing)")
        ax.semilogy(a, np.maximum(r, 1e-16), lw=1.2, ls="--", label="rate-only break")
        hard = [x for x, v in zip(a, d["v"]) if v == "Hard"]
        if hard:
            ax.axvline(min(hard), color="0.5", lw=0.9, ls=":")
            ax.text(min(hard), 1e-13, f" first Hard {min(hard):g}°", fontsize=8, color="0.35")
        ax.set_title(titles[edge], fontsize=11)
        ax.set_xlabel(xlabels[edge]); ax.set_ylabel("distance from the mirror")
        ax.set_ylim(1e-16, 1.0); ax.grid(alpha=0.25); ax.legend(fontsize=8, loc="lower right")
    fig.suptitle(f"the two coasts of the protected island (N={N}, γ={GAMMA}, Z-dephasing, no bonds)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f"wrote {path}")


def report(rows):
    """Counts below EXCLUDE the protected angle 0, where every column is noise."""
    for edge in ("field", "frustration"):
        d = [r for r in rows if r[0] == edge and r[1] > 0]
        deep = max(d, key=lambda x: x[3])
        pure = [x for x in d if abs(x[3] - x[4]) <= 1e-9 * max(1.0, x[4])]
        fe = max(x[5] for x in d)
        print(f"{edge:12s} deepest {deep[3]:.6f} at {deep[1]:.1f}deg (rate share {deep[4]/deep[3]:.4f})"
              f" | rates are the whole break at {len(pure)}/{len(d)} angles off the protected point"
              f" | worst freq-only {fe:.2e}")


if __name__ == "__main__":
    rows = sweep()
    write_tsv(rows)
    plot(rows)
    report(rows)
