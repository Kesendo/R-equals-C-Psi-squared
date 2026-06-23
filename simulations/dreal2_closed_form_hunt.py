"""Fourth symphony: hunt per-topology closed forms for the inner degeneracy d_real(2).

DEGENERACY_PALINDROME.md (Result 7) tabulated d_real(2) but found no universal formula; the chain
is even a null result ([3,6,14,14,19,22] matches nothing). The question: do the SYMMETRIC topologies
(star S_{N-1}, complete S_N, ring C_N) have cleaner per-family closed forms, the way the structural
ceiling did (chain none, complete 4/N, star 4/(N-1))?

d_real(2) = # purely-real Liouvillian eigenvalues at grid k=2 (Re = -2*gamma_q): the weight-2
[H,A]=0 commutant modes. Read straight from the rmt CSVs. Picks up ring/star/complete N=7 the moment
the background export lands (4 points instead of 3 -> a fit can be VERIFIED, not just interpolated).

Run:  python simulations/dreal2_closed_form_hunt.py
"""
from pathlib import Path
from math import comb
import numpy as np

RESULTS = Path(__file__).parent / "results"
GAMMA_Q = 0.1
TOL = 1e-8
TOPOS = ("chain", "ring", "star", "complete")


def load(topo, N):
    name = f"rmt_eigenvalues_N{N}.csv" if topo == "chain" else f"rmt_eigenvalues_{topo}_N{N}.csv"
    p = RESULTS / name
    if not p.exists():
        return None
    data = np.loadtxt(p, delimiter="\t", skiprows=1)
    return data[:, 0] + 1j * data[:, 1]


def d_real_2(topo, N):
    ev = load(topo, N)
    if ev is None:
        return None
    real = ev[np.abs(ev.imag) < TOL].real
    return int(np.sum(np.abs(real - (-2 * GAMMA_Q)) < TOL))


# ---- collect the sequences ----
seqs = {}
for topo in TOPOS:
    row = {}
    for N in range(3, 8):
        v = d_real_2(topo, N)
        if v is not None:
            row[N] = v
    seqs[topo] = row

print("=" * 96)
print("d_real(2) sequences (from the rmt CSVs; N=7 appears when the background export finishes)")
print("=" * 96)
for topo in TOPOS:
    pts = seqs[topo]
    print(f"  {topo:9}: " + "  ".join(f"N={N}:{v}" for N, v in sorted(pts.items())))


def poly_fit_report(name, xs, ys, deg):
    """Fit a degree-deg polynomial; report integer-ness and verification on a held-out point."""
    xs, ys = np.array(xs, float), np.array(ys, float)
    if len(xs) < deg + 1:
        return f"    deg{deg}: need {deg+1} pts, have {len(xs)}"
    # fit on the first deg+1 points, verify on any remaining
    fitx, fity = xs[:deg + 1], ys[:deg + 1]
    coef = np.polyfit(fitx, fity, deg)
    pred_all = np.polyval(coef, xs)
    integerish = np.allclose(coef, np.round(coef), atol=1e-6)
    held = xs[deg + 1:]
    if len(held):
        ok = np.allclose(pred_all[deg + 1:], ys[deg + 1:], atol=1e-6)
        verdict = f"VERIFIES on N={[int(h) for h in held]}" if ok else f"FAILS on N={[int(h) for h in held]}"
    else:
        verdict = "(no held-out point yet)"
    cs = ", ".join(f"{c:.4g}" for c in coef)
    return f"    deg{deg}: coef=[{cs}] int={integerish} {verdict}"


def line_through(p0, p1):
    """Integer line a*N+b through two (N,value) points, or None if non-integer."""
    (x0, y0), (x1, y1) = p0, p1
    if (y1 - y0) % (x1 - x0) != 0:
        return None
    a = (y1 - y0) // (x1 - x0)
    b = y0 - a * x0
    return a, b


# ---- even/odd-separated hunt (the key: k=2 is center at N=4, off-center at N>=6) ----
print("\n" + "=" * 96)
print("EVEN/ODD-SEPARATED HUNT. Note: d_real(2) is a BOUNDARY value for N<=3 (k=2 >= N-1), so")
print("the 'inner' regime is N>=4. A line through two points is VERIFIED by a third (gate).")
print("=" * 96)
for topo in TOPOS:
    pts = sorted(seqs[topo].items())
    print(f"\n  {topo}:")
    for parity, lab in ((1, "odd "), (0, "even")):
        sub = [(n, v) for n, v in pts if n % 2 == parity]
        bnd = [(n, v) for n, v in sub if n <= 3]
        inner = [(n, v) for n, v in sub if n >= 4]
        tag = f"(boundary N<=3: {bnd})" if bnd else ""
        if len(inner) >= 2:
            ln = line_through(inner[0], inner[1])
            if ln:
                a, b = ln
                # verify on any further inner point AND note if boundary N=3 also fits
                preds = {n: a * n + b for n, _ in sub}
                hits = all(preds[n] == v for n, v in sub)
                gate = "VERIFIES all incl boundary" if hits else \
                       ("verifies inner" if all(preds[n] == v for n, v in inner) else "FAILS")
                third = [n for n, _ in inner[2:]]
                gtxt = f"+ N={third} held-out" if third else "(needs a 3rd inner point: N+2)"
                print(f"    {lab} inner {inner}: d_real(2) = {a}N{b:+d}  [{gate}] {gtxt} {tag}")
            else:
                print(f"    {lab} inner {inner}: non-integer slope {tag}")
        else:
            print(f"    {lab} inner {inner}: need a 2nd inner point (N=7 export pending) {tag}")

print("\n(chain odd is the proof of concept: 6,14,22 at N=3,5,7 = 4N-6, which the document missed by")
print(" mixing parities. The N=7 export gives the 2nd odd-inner point for ring/star/complete.)")
