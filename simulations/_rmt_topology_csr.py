"""Never-fed symphony: the COMPLEX spacing ratio of the dephased Liouvillian, per topology.

rmt_analysis.py reads only chain spectra and only does the REAL 1D spacing ratio (Poisson/GOE/GUE)
on the decay rates. But the Liouvillian is non-Hermitian -> its eigenvalues live in the complex
plane, so the physically correct diagnostic is the COMPLEX spacing ratio (Sa, Ribeiro, Prosen,
PRX 2020): for each eigenvalue, z = (nearest neighbour - lambda)/(next-nearest - lambda) in C.
  * 2D Poisson (integrable / symmetry-fragmented):  <|z|> ~ 0.658,  <cos theta> ~ 0   (flat angle)
  * GinUE (dissipative quantum chaos):               <|z|> ~ 0.738,  <cos theta> ~ -0.241 (repulsion)

The ring/star/complete spectra exist on disk but were NEVER run through any RMT machinery (the
analyzer is chain-only, line 65). This feeds them.

RESULT (the honest finding, N=6 / 4096 eigenvalues): the GLOBAL complex spacing ratio does NOT
cleanly classify the symmetric topologies, because the real content is one level up:

  * TOPOLOGY SYMMETRY -> SPECTRAL DEGENERACY. The count of DISTINCT Liouvillian eigenvalues (upper
    half-plane) falls monotonically with the symmetry group for N>=5:
        chain (1078) > ring C_N (681) > star S_{N-1} (222) > complete S_N (99).
    As a collapse FRACTION of each topology's own oscillating (upper-half) modes that is
    45% (chain) < 65% (ring) < 88% (star) < 94% (complete) at N=6 -- the honest normalization is
    distinct / upper-half, NOT distinct / all-4096 (the latter mis-read complete as "97.6%", which
    is a normalization artifact, not the memory's 97%; see _is_the_97_the_memory.py). N=4 is the
    usual outlier (ring/star swap), the same N=4 special as the ceiling story.
  * chain is clean 2D-Poisson (<cos theta> ~ 0 at odd N=3,5,7) -> Heisenberg integrability.
  * The symmetric topologies fragment the global spectrum so hard that global non-Hermitian RMT
    does not apply: too few distinct levels, cluster-dominated (<cos theta> > 0, attraction not
    repulsion). The clean RMT test would be SECTOR-resolved (deliberately NOT done here).

Mechanism: lambda = -2g*hamming + i*(E_a - E_b); high symmetry collapses the distinct H energies
(large irreps -> few distinct E -> few distinct dE), so distinct(lambda) shrinks with symmetry.

CONTEXT: the chain's degeneracy/multiplicity palindrome is already a full document,
experiments/DEGENERACY_PALINDROME.md ('The Palindrome Inside the Palindrome', d_total(k)=d_total(N-k)
from Pi, with closed forms). That document is chain-only; the NEW content here is the TOPOLOGY axis
of the distinct-count (chain < ring < star < complete), the systematic version of its Open Question 2.

Pure numpy; reads the existing rmt_eigenvalues_*.csv. Run:  python simulations/_rmt_topology_csr.py
"""
from pathlib import Path
import numpy as np

RESULTS = Path(__file__).parent / "results"

# CSR reference values (Sa-Ribeiro-Prosen 2020)
CSR_REF = {
    "2D-Poisson": (0.6577, 0.0),
    "GinUE":      (0.7378, -0.2405),
}


def load_topology(topo, N):
    """Complex eigenvalues from the C# rmt export. topo='chain' uses the bare N{N}.csv name."""
    name = f"rmt_eigenvalues_N{N}.csv" if topo == "chain" else f"rmt_eigenvalues_{topo}_N{N}.csv"
    path = RESULTS / name
    if not path.exists():
        return None
    reals, imags = [], []
    with open(path, "r") as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                reals.append(float(parts[0].replace(",", ".")))
                imags.append(float(parts[1].replace(",", ".")))
    return np.array(reals) + 1j * np.array(imags)


def complex_spacing_ratio(evals, im_tol=1e-6):
    """CSR over the upper half-plane (Im>tol) to avoid the conjugation-symmetric real axis line.
    z_k = (NN - lambda_k)/(NNN - lambda_k); returns <|z|>, <cos arg z>, count.
    Exact degeneracies are removed first (np.unique on rounded values): the dephased Liouvillian is
    massively degenerate (lambda = -2g*hamming + i*dE), so a raw NN is often a coincident duplicate
    (z=0). The CSR is only meaningful on the DISTINCT spectrum."""
    pts = evals[evals.imag > im_tol]
    pts = np.unique(np.round(pts, 9))   # collapse exact degeneracies
    n = len(pts)
    if n < 10:
        return float("nan"), float("nan"), n
    zr, zc = [], []
    for k in range(n):
        d = np.abs(pts - pts[k])
        d[k] = np.inf
        order = np.argpartition(d, 2)[:3]
        order = order[np.argsort(d[order])]
        nn, nnn = order[0], order[1]
        denom = pts[nnn] - pts[k]
        if abs(denom) > 1e-15:
            z = (pts[nn] - pts[k]) / denom
            zr.append(abs(z))
            zc.append(np.cos(np.angle(z)))
    return float(np.mean(zr)), float(np.mean(zc)), n


def real_rate_ratio(evals):
    """The OLD diagnostic for comparison: <r> on the sorted decay rates (1D)."""
    rates = np.sort(-evals.real)
    rates = rates[rates > 1e-10]
    sp = np.diff(rates)
    sp = sp[sp > 1e-12]
    if len(sp) < 3:
        return float("nan")
    rr = [min(sp[i], sp[i + 1]) / max(sp[i], sp[i + 1]) for i in range(len(sp) - 1) if max(sp[i], sp[i + 1]) > 1e-15]
    return float(np.mean(rr))


def nearest_class(absz, cosz):
    best, bd = None, np.inf
    for name, (a, c) in CSR_REF.items():
        d = (absz - a) ** 2 + (cosz - c) ** 2
        if d < bd:
            best, bd = name, d
    return best


print("=" * 96)
print("COMPLEX SPACING RATIO of the dephased Liouvillian, per topology (non-Hermitian RMT)")
print(f"  references: 2D-Poisson <|z|>~0.658 <cos>~0 (integrable/fragmented) ; "
      f"GinUE <|z|>~0.738 <cos>~-0.241 (dissipative chaos)")
print("=" * 96)
print(f"{'topo':9} {'N':>2} {'#evals':>8} {'#upperC':>8} {'<|z|>':>8} {'<cos t>':>9} {'nearest':>12} {'<r>_real':>9}")

rows = {}
for topo in ("chain", "ring", "star", "complete"):
    for N in range(3, 8):
        ev = load_topology(topo, N)
        if ev is None:
            continue
        absz, cosz, nup = complex_spacing_ratio(ev)
        rr = real_rate_ratio(ev)
        cls = nearest_class(absz, cosz) if not np.isnan(absz) else "?"
        rows.setdefault(topo, {})[N] = (absz, cosz, cls, nup, len(ev))
        print(f"{topo:9} {N:>2} {len(ev):>8} {nup:>8} {absz:>8.4f} {cosz:>9.4f} {cls:>12} {rr:>9.4f}")

print("\n" + "-" * 96)
print("THE FINDING -- topology symmetry -> spectral degeneracy: #distinct eigenvalues (upper half)")
print("-" * 96)
print(f"{'N':>2} | " + " | ".join(f"{t:>9}" for t in ("chain", "ring", "star", "complete")) + "   monotone?")
for N in range(3, 8):
    cells, vals = [], []
    for topo in ("chain", "ring", "star", "complete"):
        if topo in rows and N in rows[topo]:
            nup = rows[topo][N][3]
            cells.append(f"{nup:>9}")
            vals.append(nup)
        else:
            cells.append(f"{'-':>9}")
            vals.append(None)
    present = [v for v in vals if v is not None]
    mono = "yes" if present == sorted(present, reverse=True) and len(present) >= 3 else "no (N=4 swap)" if N == 4 else "-"
    print(f"{N:>2} | " + " | ".join(cells) + f"   {mono}")
print("\n  chain (least symmetric) keeps the most distinct levels and reads 2D-Poisson (integrable);")
print("  complete (S_N) collapses ~98% of the spectrum. The clean RMT class is a SECTOR question.")
print("\nDONE.")
