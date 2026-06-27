"""
rmt_goe_hint_verdict.py

Drive the lone open RMT question to a verdict:
  RANDOM_MATRIX_THEORY.md flags a within-sector ⟨r⟩=0.513 ("GOE-like") at N=5,
  w=2/3 on ~15 frequencies as "the most important open question" -- the single
  reading that is NOT clean Poisson. Is it a real GOE signature or small-sample
  noise on top of the integrable (Poisson) frequency lattice?

Test design (gate-first):
  Within a fixed decay-rate (XY-weight) sector all eigenvalues sit on one
  vertical line -Rate + i*freq, so the spectrum is effectively 1D -> the REAL
  Wigner-Dyson spacing ratio is the right statistic (not the complex CSR).
  Three drivers:
   (1) reproduce the 0.513 (Python-rebuilt L at N=5, the historical code path).
   (2) bootstrap the Poisson sampling band for <r> at the observed sample size
       (uniform points = homogeneous Poisson process, same spacing_ratios fn).
   (3) extend to N=6,7 from the C# CSV exports (classify by rate) -> larger
       samples. Poisson => converges to 0.386; real GOE => stays ~0.536 / grows.

GATE: if the observed <r> sits inside the Poisson-n band AND the larger-N
samples converge toward 0.386, the hint is small-sample noise (verdict: artifact,
close the open question). If it stays >=0.5 and the Poisson band excludes it
across N, the hint is real (verdict: escalate).
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import rmt_analysis as rmt  # reuse build_H_chain, build_L_zdeph, load_eigenvalues, spacing_ratios

GAMMA = 0.05
RNG = np.random.default_rng(20260627)
POISSON_R = 2 * np.log(2) - 1  # 0.3863
GOE_R = 0.536


def sector_freqs_from_evals(evals, gamma, w, tol_frac=0.3):
    """Frequencies (|Im|, unique, rounded) in the XY-weight-w sector (rate=2*w*gamma)."""
    rates = -evals.real
    freqs = np.abs(evals.imag)
    tol = tol_frac * gamma
    expected = 2 * w * gamma
    mask = (rates < 1e-8) if w == 0 else (np.abs(rates - expected) < tol)
    sf = freqs[mask]
    sf = sf[sf > 1e-6]
    return np.sort(np.unique(np.round(sf, 8)))


def poisson_band(n, n_mc=20000):
    """Sampling distribution of <r> for n points drawn from a homogeneous
    Poisson process (uniform on [0,1]), via the same spacing_ratios statistic."""
    vals = np.empty(n_mc)
    for i in range(n_mc):
        pts = np.sort(RNG.random(n))
        r, _ = rmt.spacing_ratios(pts)
        vals[i] = r
    vals = vals[~np.isnan(vals)]
    return {
        "mean": float(np.mean(vals)),
        "std": float(np.std(vals)),
        "p05": float(np.percentile(vals, 5)),
        "p95": float(np.percentile(vals, 95)),
        "frac_ge_obs": vals,  # filled per-call below
    }


def main():
    print("=" * 78)
    print("RMT within-sector GOE-hint verdict  (gamma=0.05, J=1.0, Heisenberg chain)")
    print("Poisson <r>=0.386, GOE <r>=0.536, GUE <r>=0.603")
    print("=" * 78)

    # ---- (1) reproduce at N=5 via the historical Python-rebuild path ----
    print("\n[1] Reproduce N=5 sectors (Python-rebuilt L, the historical path):")
    H = rmt.build_H_chain(5, 1.0)
    L = rmt.build_L_zdeph(H, GAMMA)
    ev5 = np.linalg.eigvals(L)
    n5_obs = {}
    for w in range(6):
        sf = sector_freqs_from_evals(ev5, GAMMA, w)
        if len(sf) >= 5:
            r, ratios = rmt.spacing_ratios(sf)
            cls, _ = rmt.classify_r(r)
            n5_obs[w] = (r, len(sf))
            print(f"    w={w}: {len(sf):3d} unique freq  <r>={r:.3f}  -> {cls}")
        elif len(sf) > 0:
            print(f"    w={w}: {len(sf):3d} unique freq  (too few)")

    # ---- (2) bootstrap Poisson band at the observed sample sizes ----
    print("\n[2] Poisson sampling band for <r> at small n (Monte Carlo, 20k draws):")
    sizes = sorted({n for (_, n) in n5_obs.values()} | {12, 15, 20, 30, 60, 120})
    bands = {}
    for n in sizes:
        b = poisson_band(n)
        bands[n] = b
        print(f"    n={n:4d}:  <r>_Poisson = {b['mean']:.3f} +/- {b['std']:.3f}"
              f"   [5%,95%] = [{b['p05']:.3f}, {b['p95']:.3f}]")

    print("\n    -> Is the N=5 observed <r> inside the Poisson band at its own n?")
    for w, (r, n) in sorted(n5_obs.items()):
        b = bands[n]
        inside = b["p05"] <= r <= b["p95"]
        # one-sided p: fraction of Poisson draws >= observed
        draws = b["frac_ge_obs"]
        p_ge = float(np.mean(draws >= r))
        verdict = "WITHIN Poisson band" if inside else "OUTSIDE Poisson band"
        print(f"    w={w}: obs <r>={r:.3f} (n={n})  Poisson[5,95]=[{b['p05']:.3f},{b['p95']:.3f}]"
              f"  p(Poisson>=obs)={p_ge:.3f}  -> {verdict}")

    # ---- (3) extend to N=6,7 from the C# CSV exports (classify by rate) ----
    print("\n[3] Larger samples from C# CSV exports (full Liouvillian, classify by rate):")
    print("    Does the w=2,3 reading converge to Poisson (0.386) with more frequencies?")
    for N in (5, 6, 7):
        evals = rmt.load_eigenvalues(N)
        if evals is None:
            print(f"    N={N}: CSV not found")
            continue
        print(f"    N={N}  ({len(evals)} eigenvalues):")
        for w in range(N + 1):
            sf = sector_freqs_from_evals(evals, GAMMA, w)
            if len(sf) >= 5:
                r, _ = rmt.spacing_ratios(sf)
                cls, _ = rmt.classify_r(r)
                b = poisson_band(len(sf), n_mc=8000)
                inside = b["p05"] <= r <= b["p95"]
                tag = "Poisson-band OK" if inside else "outside Poisson band"
                print(f"        w={w}: {len(sf):4d} freq  <r>={r:.3f} -> {cls:8s}"
                      f"  (Poisson n={len(sf)}: {b['mean']:.3f}+/-{b['std']:.3f}, {tag})")

    print("\n" + "=" * 78)
    print("Read the verdict from [2] (is 0.513 a Poisson fluctuation at n~15?)")
    print("and [3] (does it converge to 0.386 with larger samples?).")
    print("=" * 78)


if __name__ == "__main__":
    main()
