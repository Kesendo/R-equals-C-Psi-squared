"""Read two figures of S. Reichert's thesis off their vector paths into data files.

S. Reichert, Advancing sodium triple quantum (TQ) nuclear magnetic resonance (NMR) spectroscopy and imaging,
dissertation, Heidelberg University (2024), doi:10.11588/heidok.00034464. Two of its figures are the first paper's
(Reichert et al., NMR in Biomedicine 37 (2024) e5106) drawn as vector paths, where the paper prints them as raster
images:

  Fig. 3.6a  the inversion-recovery FID of the IRTQTPPI sequence, 2 % agarose, phase cycle PC0 (the paper's Fig. 2A);
  Fig. 4.11  the fitted T1s, T1f, T1m and T1-TQ/SQ ratio against the number of coarse points added, 2 %, PC2 (the
             paper's Fig. 5).

Usage: python simulations/read_reichert_thesis_figures.py PATH/TO/Dissertation_SimonReichert_PDFA.pdf \
           > data/reichert2024_sodium_agarose/reader_output.txt
The PDF is not in the repository; the data files this writes are, under data/reichert2024_sodium_agarose/, read by
experiments/SODIUM_ONE_POOL_OR_MORE.md through simulations/sodium_one_pool_or_more.py (Part 9). Each axis is
calibrated by a straight line through its tick marks, whose positions the path's own tick segments give; the
residuals of that line and the coordinate quantum of the figure's content stream are printed.
"""
import os
import re
import sys
import numpy as np
import fitz   # PyMuPDF

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "reichert2024_sodium_agarose")
CREDIT = ("# Data: Simon Reichert's measurement (with V. Schepkin, D. Kleimaier, F. G. Zöllner and L. R. Schad), from his "
          "dissertation, Heidelberg University 2024, doi:10.11588/heidok.00034464; please cite it when using these values\n")


def segments(page):
    for d in page.get_drawings():
        for it in d["items"]:
            if it[0] == "l":
                yield it[1], it[2]


def ticks(page, axis, at, lo, hi, length=3.0):
    """Positions of the short tick segments on one axis line: axis "x" gives vertical segments starting on the
    horizontal line y = at, axis "y" horizontal segments starting on the vertical line x = at, within (lo, hi)."""
    found = set()
    for a, b in segments(page):
        if axis == "x" and abs(a.x - b.x) < 1e-6 and abs(a.y - at) < 1e-2 and 0 < abs(a.y - b.y) < length and lo < a.x < hi:
            found.add(round(a.x, 4))
        if axis == "y" and abs(a.y - b.y) < 1e-6 and abs(a.x - at) < 1e-2 and 0 < abs(a.x - b.x) < length and lo < a.y < hi:
            found.add(round(a.y, 4))
    return np.array(sorted(found))


def calibrate(pos, values, name):
    assert len(pos) == len(values), (name, len(pos), len(values))
    c = np.polyfit(pos, values, 1)
    res = np.abs(np.polyval(c, pos) - values).max()
    print(f"  {name}: {len(pos)} tick marks, {c[0]:.6f} units per pt, largest residual {res:.2e}")
    return c, res


def vertices(path):
    pts = [(path["items"][0][1].x, path["items"][0][1].y)]
    for it in path["items"]:
        assert it[0] == "l", "a path of straight segments"
        pts.append((it[2].x, it[2].y))
    return np.array(pts)


def content_quantum(doc, page, xy):
    """The coordinate step of the figure's form XObject: the most decimals its longest run of path operators
    carries, times the scale from the stream's units to the page's, read off that run against the path's vertices."""
    xrefs = [x[0] for x in page.get_xobjects()]
    assert len(xrefs) == 1, xrefs
    runs, cur = [], []
    for line in doc.xref_stream(xrefs[0]).decode("latin-1").split("\n"):
        if re.fullmatch(r"-?[\d.]+ -?[\d.]+ [ml]", line):
            cur.append(line)
        else:
            runs.append(cur)
            cur = []
    run = max(runs + [cur], key=len)
    assert len(run) == len(xy), (len(run), len(xy))
    raw = np.array([[float(v) for v in line.split()[:2]] for line in run])
    scale = np.polyfit(raw[:, 0], xy[:, 0], 1)[0]
    decimals = max(len(v.split(".")[1]) if "." in v else 0 for line in run for v in line.split()[:2])
    return 10.0 ** -decimals * scale


def fig36(doc):
    page = doc[76]
    print("Fig. 3.6a (PDF page 77): the IRTQTPPI FID, 2 %, PC0")
    fid = max(page.get_drawings(), key=lambda d: len(d["items"]))
    xy = vertices(fid)
    cx, _ = calibrate(ticks(page, "x", 273.26, 130, 315), np.arange(0, 601, 100), "time axis (ms)")
    cy, _ = calibrate(ticks(page, "y", 139.19, 100, 275)[::-1], np.linspace(-1, 1, 11), "signal axis")
    t, y = np.polyval(cx, xy[:, 0]), np.polyval(cy, xy[:, 1])
    q = content_quantum(doc, page, xy)
    print(f"  {len(t)} vertices from {t[0]:.4f} to {t[-1]:.4f} ms; the content stream's coordinate step is "
          f"{q:.5f} pt, {q * cx[0]:.4f} ms on the time axis")
    path = os.path.join(OUT, "thesis_fig3_6a_irtqtppi_fid_2pct_pc0.csv")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(CREDIT)
        f.write("# Its Fig. 3.6a (the first paper's Fig. 2A, Reichert et al., NMR in Biomedicine 37 (2024) e5106): IRTQTPPI "
                "FID, 2 % agarose, PC0, 9.4 T; read off the vector path by simulations/read_reichert_thesis_figures.py\n")
        f.write("t_ms,signal\n")
        for a, b in zip(t, y):
            f.write(f"{a:.4f},{b:.5f}\n")
    print(f"  written: data/reichert2024_sodium_agarose/{os.path.basename(path)}")


def fig411(doc):
    page = doc[101]
    print("Fig. 4.11 (PDF page 102): the fit against the coarse points added, 2 %, PC2")
    curves = [d for d in page.get_drawings() if len(d["items"]) == 12]
    assert len(curves) == 4, len(curves)
    ca, ra = calibrate(ticks(page, "x", 604.5, 105, 285), np.arange(0, 301, 50), "panel a, points added")
    cya, _ = calibrate(ticks(page, "y", 109.376, 430, 606)[::-1], np.arange(30, 51, 2), "panel a, time (ms)")
    cb, rb = calibrate(ticks(page, "x", 604.43, 300, 480), np.arange(0, 301, 50), "panel b, points added")
    cyb, _ = calibrate(ticks(page, "y", 306.707, 430, 606)[::-1], np.arange(4.5, 9.51, 0.5), "panel b, TQ/SQ (%)")
    a = [d for d in curves if d["rect"].x0 < 290]
    b = [vertices(d) for d in curves if d["rect"].x0 > 290]
    assert len(a) == 3 and len(b) == 1
    # panel a's legend: a short stroke in each curve's colour beside the label T1s, T1f or T1m
    words = [w for w in page.get_text("words") if w[4] in ("T1s", "T1f", "T1m") and w[0] < 290 and 430 < w[1] < 606]
    labels = {w[4]: (w[1] + w[3]) / 2 for w in words}
    left = min(w[0] for w in words)              # the legend's markers sit just left of its labels
    colours = {tuple(d["color"]) for d in a}
    marks = {}
    for d in page.get_drawings():
        mid = (d["rect"].y0 + d["rect"].y1) / 2
        if d["color"] and tuple(d["color"]) in colours and len(d["items"]) <= 2 and left - 15 < d["rect"].x0 < left \
                and min(labels.values()) - 8 < mid < max(labels.values()) + 8:
            marks.setdefault(tuple(d["color"]), []).append(mid)
    assert len(labels) == 3 and len(marks) == 3, (labels, marks)
    label_of = {c: min(labels, key=lambda k: abs(labels[k] - np.mean(ys))) for c, ys in marks.items()}
    assert sorted(label_of.values()) == ["T1f", "T1m", "T1s"], label_of
    times = {label_of[tuple(d["color"])]: np.polyval(cya, vertices(d)[:, 1]) for d in a}
    n = np.polyval(ca, vertices(a[0])[:, 0])
    assert all(np.abs(np.polyval(ca, vertices(d)[:, 0]) - n).max() < 1e-9 for d in a), "one x per point in panel a"
    assert np.abs(np.polyval(cb, b[0][:, 0]) - n).max() <= ra + rb, "the two panels share their points within the axes"
    ratio = np.polyval(cyb, b[0][:, 1])
    path = os.path.join(OUT, "thesis_fig4_11_fit_vs_points_added_2pct_pc2.csv")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(CREDIT)
        f.write("# Its Fig. 4.11 (the first paper's Fig. 5, Reichert et al., NMR in Biomedicine 37 (2024) e5106): IRTQTPPI "
                "fit, 2 % agarose, PC2, 9.4 T, against the coarse points added; read off the vector paths, each curve "
                "named by the legend's colour, by simulations/read_reichert_thesis_figures.py\n")
        f.write("points_added,T1s_ms,T1m_ms,T1f_ms,tq_sq_percent\n")
        for row in zip(n, times["T1s"], times["T1m"], times["T1f"], ratio):
            f.write(f"{np.round(row[0], 2) + 0.0:.2f},{row[1]:.3f},{row[2]:.3f},{row[3]:.3f},{row[4]:.3f}\n")
    print(f"  {len(n)} points added from {np.round(n[0], 2) + 0.0:.2f} to {n[-1]:.2f}; written: data/reichert2024_sodium_agarose/{os.path.basename(path)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    os.makedirs(OUT, exist_ok=True)
    doc = fitz.open(sys.argv[1])
    fig36(doc)
    fig411(doc)
