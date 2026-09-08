"""
rmt_analysis.py
Finite-N spectral-spacing diagnostics for the palindromic Liouvillian.
Reads eigenvalue CSVs exported by C# engine (dotnet run -- rmt).

Reports raw rate-spacing summaries, consecutive gap ratios, the centered
decay-rate reflection residual, and complex-plane bounds. The raw full multiset retains
zero gaps; no standard-ensemble calibration or irreducible-class inference.
"""

import numpy as np
from pathlib import Path

# ---- paths ----
RESULTS = Path(__file__).parent / "results"

# ---- Pauli + Liouvillian (for sector analysis N=2-5) ----
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(op, site, N):
    r = np.eye(1, dtype=complex)
    for k in range(N):
        r = np.kron(r, op if k == site else I2)
    return r


def build_H_chain(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * kron_at(P, i, N) @ kron_at(P, i + 1, N)
    return H


def build_L_zdeph(H, gamma):
    d = H.shape[0]
    N = int(np.log2(d))
    Id = np.eye(d, dtype=complex)
    d2 = d * d
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(N):
        Lk = np.sqrt(gamma) * kron_at(sz, k, N)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)
    return L


# ---- load eigenvalues ----
def load_eigenvalues(N):
    """Load complex eigenvalues from C# CSV export.
    Handles both '.' and ',' as decimal separator (German locale)."""
    path = RESULTS / f"rmt_eigenvalues_N{N}.csv"
    if not path.exists():
        return None
    reals, imags = [], []
    with open(path, 'r') as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                reals.append(float(parts[0].replace(',', '.')))
                imags.append(float(parts[1].replace(',', '.')))
    return np.array(reals) + 1j * np.array(imags)


# ---- consecutive gap ratio ----
def spacing_ratios(sorted_vals):
    """Consecutive spacing ratios r_n = min(s_n,s_{n+1})/max(s_n,s_{n+1}).
    Each original adjacent gap pair is retained: 0/positive is zero, 0/0 is
    NaN in its original position. The mean uses defined ratios only. There is
    no cutoff or deduplication. This statistic is affine-invariant in exact
    arithmetic, not invariant under nonlinear reparameterizations."""
    spacings = np.diff(sorted_vals)
    if len(spacings) < 2:
        return float('nan'), np.array([])
    maxima = np.maximum(spacings[:-1], spacings[1:])
    ratios = np.full(len(maxima), np.nan)
    defined = maxima > 0
    ratios[defined] = np.minimum(spacings[:-1], spacings[1:])[defined] / maxima[defined]
    return (float(np.mean(ratios[defined])) if np.any(defined) else float('nan')), ratios


def spacing_population(sorted_vals):
    """Population accounting for the parsed-double multiset, including ties."""
    mean, ratios = spacing_ratios(sorted_vals)
    return dict(levels=len(sorted_vals), gaps=max(0, len(sorted_vals) - 1),
                zero_gaps=int(np.count_nonzero(np.diff(sorted_vals) == 0)),
                adjacent_pairs=len(ratios), defined_ratios=int(np.count_nonzero(~np.isnan(ratios))),
                undefined_zero_zero=int(np.count_nonzero(np.isnan(ratios))), mean=mean)


def population_lines(reading):
    keys = ("levels", "gaps", "zero_gaps", "adjacent_pairs", "defined_ratios", "undefined_zero_zero")
    return ["      " + ", ".join(f"{key}={reading[key]}" for key in keys),
            f"      <r> defined = {reading['mean']:.17g}"]


def centered_rate_projection_reading(eigenvalues, sigma_gamma, tolerance=1e-10):
    """Multiplicity-sensitive reflection check on the centered real projections.

    This deliberately does not claim to match the full complex spectrum.  It does,
    however, require a perfect matching of *all* projected rates, including central
    multiplicity, so an unmatched value cannot be hidden by truncating sign lists.
    """
    centered_rates = -(np.asarray(eigenvalues, dtype=complex) + sigma_gamma).real
    reflected = -centered_rates
    # On the real line, monotone matching minimizes every L-infinity matching
    # threshold.  Equal-length sorted arrays therefore give the exact bottleneck
    # perfect-matching distance without materializing an O(n^2) cost matrix.
    bottleneck_error = (0.0 if len(centered_rates) == 0 else
                        float(np.max(np.abs(np.sort(centered_rates) -
                                            np.sort(reflected)))))
    return {
        "positive": int(np.count_nonzero(centered_rates > tolerance)),
        "negative": int(np.count_nonzero(centered_rates < -tolerance)),
        "central": int(np.count_nonzero(np.abs(centered_rates) <= tolerance)),
        "bottleneck_error": bottleneck_error,
    }
# ==================================================================
# MAIN ANALYSIS
# ==================================================================
def main():
    out = []
    out.append("=" * 70)
    out.append("FINITE-N SPECTRAL-SPACING DIAGNOSTICS")
    out.append("Palindromic Liouvillian, Heisenberg Chain, Z-dephasing")
    out.append("gamma = 0.05, J = 1.0")
    out.append("Raw full multiset: all parsed rates -Re(lambda), including zero and negative numerical rates; no cutoff or deduplication.")
    out.append("Each original adjacent gap pair is retained: 0/positive = 0; 0/0 = undefined, counted separately; mean over defined ratios only.")
    out.append("Zero gaps mean exact parsed-double ties; near-degenerate numerical splittings are not certified physical resolution.")
    out.append("No standard Poisson/GOE/GUE calibration applies to this unresolved, degenerate pooled population.")
    out.append("No nearest-ensemble verdict and no irreducible-class inference are made.")
    out.append("=" * 70)

    gamma = 0.05

    # ---- Phase 2+3: Spacing ratios for each N ----
    r_all_table = {}
    r_half_table = {}

    for N in range(2, 8):
        evals = load_eigenvalues(N)
        if evals is None:
            out.append(f"\nN={N}: CSV not found, skipping.")
            continue

        Sg = N * gamma
        rates = -evals.real
        n_total = len(evals)

        out.append(f"\n{'='*70}")
        out.append(f"N={N}: {n_total} eigenvalues")
        out.append(f"{'='*70}")

        # (a) All rates: spacing ratio
        rates_all = np.sort(rates)
        reading = spacing_population(rates_all)
        out.append("\n  (a) Raw full multiset:")
        out.extend(population_lines(reading))

        # Spacing stats
        sp = np.diff(rates_all)
        if len(sp) > 5:
            out.append(f"      Spacing stats: mean={np.mean(sp):.6f}, "
                       f"std={np.std(sp):.6f}, "
                       f"min={np.min(sp):.2e}, max={np.max(sp):.6f}")

        r_all_table[N] = reading['mean']

        # (b) Lower half only (Re < Sigma_gamma)
        rates_lower = np.sort(rates[rates < Sg])
        lower = spacing_population(rates_lower)
        out.append("\n  (b) Raw rates below N*gamma (not an irreducible sector or a decorrelation certificate):")
        out.extend(population_lines(lower))
        r_half_table[N] = lower['mean']

        # ---- Centered decay-rate projection check ----
        projection = centered_rate_projection_reading(evals, Sg)
        out.append(f"\n  Centered decay-rate projection reflection (Sg={Sg}):")
        out.append("      Scope: multiplicity-preserving perfect matching of the full real-part projection; does not match the full complex multiset.")
        out.append(f"      positive={projection['positive']}, negative={projection['negative']}, "
                   f"central={projection['central']}, multiplicity bottleneck error: "
                   f"{projection['bottleneck_error']:.2e}")

        # ---- Phase 7: Complex plane ----
        nonzero = evals[np.abs(evals) > 1e-10]
        out.append(f"\n  Complex plane:")
        out.append(f"      Fraction Re<0: "
                   f"{np.sum(nonzero.real < 0)/len(nonzero):.4f}")
        out.append(f"      |lambda| max: {np.max(np.abs(evals)):.4f}")

    # ---- Summary table ----
    out.append(f"\n{'='*70}")
    out.append(f"SUMMARY")
    out.append(f"{'='*70}")

    out.append("\n  Direct multiset spacing ratio <r> over defined original adjacent pairs:")
    out.append(f"  {'N':>3} | {'<r> all':>8} | {'<r> half':>9}")
    out.append(f"  " + "-" * 28)
    for N in sorted(r_all_table.keys()):
        ra = r_all_table.get(N, float('nan'))
        rh = r_half_table.get(N, float('nan'))
        out.append(f"  {N:3d} | {ra:8.4f} | {rh:9.4f}")

    out.append("\n  Scope: gap ratios are affine-invariant in exact arithmetic; nonlinear maps can change them.")
    out.append("  These numeric multiset summaries are sensitive to near-degenerate eigensolver splittings; no irreducible-class inference.")

    text = "\n".join(out)
    print(text)

    out_path = RESULTS / "rmt_analysis.txt"
    out_path.write_text(text, encoding="utf-8")
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
