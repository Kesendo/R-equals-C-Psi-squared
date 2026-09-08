"""Never-fed symphony: the COMPLEX spacing ratio of the dephased Liouvillian, per topology.

rmt_analysis.py reads chain spectra and reports raw decay-rate multiset gap ratios,
without a standard-ensemble calibration. The Liouvillian is non-Hermitian: its eigenvalues live in the complex
plane, so the complex spacing ratio (Sa, Ribeiro, Prosen, PRX 2020) is the
appropriate planar diagnostic geometry: for each eigenvalue,
z = (nearest neighbour - lambda)/(next-nearest - lambda) in C. The global,
unresolved population used here is still not a universality-class diagnostic.
  * 2D Poisson (uncorrelated reference; not an integrability verdict): <|z|> ~ 0.658, <cos theta> ~ 0
  * GinUE (dissipative quantum chaos):               <|z|> ~ 0.738,  <cos theta> ~ -0.241 (repulsion)

The ring/star/complete spectra exist on disk but were NEVER run through any RMT machinery (the
analyzer is chain-only, line 65). This feeds them.

RESULT (the honest finding, N=6 / 4096 eigenvalues): the GLOBAL complex spacing ratio does NOT
cleanly classify the symmetric topologies, because the real content is one level up:

  * TOPOLOGY SYMMETRY -> SPECTRAL CLUSTERING. The count of 1e-9-rounded,
    finite-precision clusters in the upper half-plane falls monotonically with
    the symmetry group for N>=5:
        chain (1078) > ring C_N (681) > star S_{N-1} (222) > complete S_N (99).
    As a clustered FRACTION of each topology's own oscillating (upper-half) modes that is
    45% (chain) < 65% (ring) < 88% (star) < 94% (complete) at N=6 -- the honest normalization is
    clustered / upper-half, NOT clustered / all-4096 (the latter mis-read complete as "97.6%", which
    is a normalization artifact, not the memory's 97%; see is_the_97_the_memory.py). N=4 is the
    usual outlier (ring/star swap), the same N=4 special as the ceiling story.
  * chain is clean 2D-Poisson (<cos theta> ~ 0 at odd N=3,5,7), compatible with integrability or fragmentation;
    this does not prove integrability.
  * The symmetric topologies fragment the global spectrum so hard that global non-Hermitian RMT
    does not apply: too few tolerance-clustered representatives, cluster-dominated (<cos theta> > 0, attraction not
    repulsion). The clean RMT test would be SECTOR-resolved (deliberately NOT done here).

The additive rate/frequency expression applies to simultaneous Hamiltonian/dissipator eigenoperators;
it is not a formula for mixed Liouvillian modes. The cluster-count comparison below is a measured,
tolerance-dependent finite-precision pattern, not an exact degeneracy census and not a general
energy-difference mechanism for every mode.

CONTEXT: the chain's degeneracy/multiplicity palindrome is already a full document,
experiments/DEGENERACY_PALINDROME.md ('The Palindrome Inside the Palindrome', d_total(k)=d_total(N-k)
from Pi, with closed forms). That document is chain-only; the NEW content here is the TOPOLOGY axis
of the tolerance-cluster count (chain < ring < star < complete), the systematic
version of its Open Question 2.

Pure numpy; reads the existing rmt_eigenvalues_*.csv. Run:  python simulations/rmt_topology_csr.py
"""
from pathlib import Path
import numpy as np

RESULTS = Path(__file__).parent / "results"
CLUSTER_DECIMALS = 9

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
    Before CSR, finite-precision clusters are represented once by rounding to
    CLUSTER_DECIMALS decimal places. This suppresses numerical copies but is not
    an exact degeneracy test; the count and CSR are tolerance-dependent. The raw
    multiset diagnostic remains separately owned by rmt_analysis.py."""
    pts = evals[evals.imag > im_tol]
    pts = np.unique(np.round(pts, CLUSTER_DECIMALS))
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
print(f"  references: 2D-Poisson <|z|>~0.658 <cos>~0 (uncorrelated reference; not an integrability verdict) ; "
      f"GinUE <|z|>~0.738 <cos>~-0.241 (dissipative chaos)")
print(f"  cluster_decimals={CLUSTER_DECIMALS}; cluster counts and CSR are tolerance-dependent finite-precision diagnostics")
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
print("THE FINDING -- topology symmetry -> fewer finite-precision clusters (upper half)")
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
print("\n  chain (least symmetric) keeps the most 1e-9-clustered representatives and reads 2D-Poisson;")
print("  this does not prove integrability. The clean RMT class is a SECTOR question.")
complete_n6 = load_topology("complete", 6)
if complete_n6 is not None:
    clustered_upper = rows["complete"][6][3]
    total_upper = int(np.count_nonzero(complete_n6.imag > 1e-6))
    clustered_fraction = 1 - clustered_upper / total_upper
    print(f"  complete N=6: clustered_upper={clustered_upper} total_upper={total_upper} clustered_fraction={clustered_fraction:.17g}")
    upper_values = complete_n6[complete_n6.imag > 1e-6]
    sensitivity = [len(np.unique(np.round(upper_values, digits))) for digits in range(7, 11)]
    print("  cluster_sensitivity_N6_complete=" + ",".join(
        f"d{digits}:{count}" for digits, count in zip(range(7, 11), sensitivity)))
print("\nDONE.")
