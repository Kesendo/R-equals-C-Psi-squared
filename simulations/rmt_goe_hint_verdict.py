"""Raw frequency-multiset diagnostics inside average-light bands.

The legacy filename is retained. These are rate windows, not invariant
fixed-XY-weight sectors. All absolute frequencies are retained without rounding
or deduplication. The mixed degenerate population has no standard-ensemble calibration.
"""
from pathlib import Path
import numpy as np

if __package__:
    from . import rmt_analysis as rmt
else:
    import rmt_analysis as rmt

GAMMA = 0.05


def band_frequency_multiset(eigenvalues, gamma, w, tol_frac=0.3):
    """Select |d-2*w*gamma| < tol_frac*gamma, then sort the full |Im| multiset.

    Under uniform dephasing, w labels average-light bands of half-width
    tol_frac/2, not fixed-weight eigenvalue sectors. Conjugate frequencies both
    remain, including zero: this is not an independent nondegenerate population.
    """
    mask = np.abs(-eigenvalues.real - 2 * w * gamma) < tol_frac * gamma
    return np.sort(np.abs(eigenvalues.imag[mask]))


def main():
    out = ["RAW FREQUENCY-MULTISET DIAGNOSTICS IN AVERAGE-LIGHT BANDS",
           "Heisenberg chain, J=1, gamma=0.05; existing C# CSVs N=2..7.",
           "Rate window |d-2*w*gamma| < 0.3*gamma; average-light half-width 0.15, not invariant XY-weight sectors.",
           "Absolute frequencies: zero frequencies and multiplicities retained; no rounding or deduplication.",
           "Each original adjacent gap pair retained; 0/positive=0; 0/0 undefined; mean over defined ratios only.",
           "Exact parsed-double ties are counted; near-degenerate splittings are not certified physical resolution.",
           "These unresolved average-light bands have no standard-ensemble calibration and no integrability verdict."]
    for n in range(2, 8):
        eigenvalues = rmt.load_eigenvalues(n)
        if eigenvalues is None:
            raise FileNotFoundError(f"Missing eigenvalue CSV N={n}")
        out.append(f"\nN={n}: {len(eigenvalues)} eigenvalues")
        for w in range(n + 1):
            reading = rmt.spacing_population(band_frequency_multiset(eigenvalues, GAMMA, w))
            out.append(f"  w={w}:")
            out.extend(rmt.population_lines(reading))
    text = "\n".join(out)
    print(text)
    output = Path(__file__).parent / "results" / "rmt_band_multiset.txt"
    output.write_text(text, encoding="utf-8")
    print(f"\nResults: {output}")


if __name__ == "__main__":
    main()
