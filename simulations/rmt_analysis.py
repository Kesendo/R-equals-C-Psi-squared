"""
rmt_analysis.py
Random Matrix Theory analysis of the palindromic Liouvillian.
Reads eigenvalue CSVs exported by C# engine (dotnet run -- rmt).

Phases:
  2: Spectral unfolding
  3: NNSD (Nearest-Neighbor Spacing Distribution)
  4: Sector decomposition (N=2-5, from Python eigenvalues)
  5: Chiral GUE test
  6: N-scaling (Brody parameter vs N)
  7: Complex plane analysis
"""

import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit, minimize_scalar
from math import comb

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


# ---- RMT distributions ----
def p_poisson(s):
    return np.exp(-s)


def p_goe(s):
    return (np.pi / 2) * s * np.exp(-np.pi * s**2 / 4)


def p_gue(s):
    return (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)


def p_brody(s, beta):
    from scipy.special import gamma as gammafn
    # Brody distribution: P(s) = c * s^beta * exp(-a * s^{beta+1})
    # with a = [Gamma((beta+2)/(beta+1))]^{beta+1}, c = (beta+1)*a
    # This ensures <s> = 1 and integral = 1.
    a = gammafn((beta + 2) / (beta + 1)) ** (beta + 1)
    c = (beta + 1) * a
    return c * s**beta * np.exp(-a * s**(beta + 1))


# ---- unfolding ----
def unfold(rates):
    """Spectral unfolding via cumulative distribution."""
    sorted_rates = np.sort(rates)
    n = len(sorted_rates)
    # Map to uniform density: x_i -> rank(x_i) / n
    unfolded = np.arange(1, n + 1) / n
    return unfolded


def nnsd(unfolded):
    """Nearest-neighbor spacing distribution from unfolded spectrum."""
    spacings = np.diff(unfolded)
    # Normalize to mean 1
    mean_s = np.mean(spacings)
    if mean_s > 0:
        spacings = spacings / mean_s
    return spacings


# ---- spacing ratio (robust, no unfolding needed) ----
def spacing_ratios(sorted_vals):
    """Consecutive spacing ratios r_n = min(s_n,s_{n+1})/max(s_n,s_{n+1}).
    Robust diagnostic: no unfolding needed.
    Reference: Poisson <r>=0.386, GOE <r>=0.536, GUE <r>=0.603."""
    spacings = np.diff(sorted_vals)
    spacings = spacings[spacings > 1e-12]  # remove degeneracies
    if len(spacings) < 3:
        return float('nan'), np.array([])
    ratios = []
    for i in range(len(spacings) - 1):
        s1, s2 = spacings[i], spacings[i + 1]
        mx = max(s1, s2)
        if mx > 1e-15:
            ratios.append(min(s1, s2) / mx)
    ratios = np.array(ratios)
    return np.mean(ratios), ratios


def classify_r(mean_r):
    """Classify spacing ratio into nearest ensemble."""
    refs = {"Poisson": 0.386, "GOE": 0.536, "GUE": 0.603}
    best = min(refs, key=lambda k: abs(refs[k] - mean_r))
    dist = abs(refs[best] - mean_r)
    return best, dist


# ---- chi-squared comparison (on unfolded spacings) ----
def chi2_compare(spacings, nbins=40):
    """Chi-squared comparison with Poisson, GOE, GUE."""
    hist, edges = np.histogram(spacings, bins=nbins, density=True,
                                range=(0, 4))
    centers = (edges[:-1] + edges[1:]) / 2
    mask = hist > 0

    results = {}
    for name, func in [("Poisson", p_poisson), ("GOE", p_goe),
                        ("GUE", p_gue)]:
        pred = np.array([func(s) for s in centers])
        chi2 = np.sum((hist[mask] - pred[mask])**2 / (pred[mask] + 1e-10))
        results[name] = chi2 / mask.sum()
    return results


# ==================================================================
# MAIN ANALYSIS
# ==================================================================
def main():
    out = []
    out.append("=" * 70)
    out.append("RANDOM MATRIX THEORY ANALYSIS")
    out.append("Palindromic Liouvillian, Heisenberg Chain, Z-dephasing")
    out.append("gamma = 0.05, J = 1.0")
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
        rates_all = np.sort(rates[rates > 1e-10])
        mean_r, ratios = spacing_ratios(rates_all)
        cls, dist = classify_r(mean_r) if not np.isnan(mean_r) else ("?", 0)

        out.append(f"\n  (a) All rates ({len(rates_all)} nonzero):")
        out.append(f"      <r> = {mean_r:.4f}  (Poisson=0.386, "
                   f"GOE=0.536, GUE=0.603)")
        out.append(f"      Nearest: {cls} (dist={dist:.3f})")

        # Spacing stats
        sp = np.diff(rates_all)
        sp = sp[sp > 1e-12]
        if len(sp) > 5:
            out.append(f"      Spacing stats: mean={np.mean(sp):.6f}, "
                       f"std={np.std(sp):.6f}, "
                       f"min={np.min(sp):.2e}, max={np.max(sp):.6f}")

        r_all_table[N] = mean_r

        # Also do NNSD chi2 on unfolded spacings
        if len(rates_all) > 30:
            sp_unf = nnsd(unfold(rates_all))
            sp_unf = sp_unf[sp_unf > 0.001]
            if len(sp_unf) > 20:
                chi2 = chi2_compare(sp_unf)
                best_chi = min(chi2, key=chi2.get)
                out.append(f"      NNSD chi2: " + ", ".join(
                    f"{k}={v:.3f}" for k, v in chi2.items())
                    + f" -> {best_chi}")

        # (b) Lower half only (Re < Sigma_gamma)
        rates_lower = np.sort(rates[(rates > 1e-10) & (rates < Sg)])
        if len(rates_lower) > 10:
            mean_r_h, _ = spacing_ratios(rates_lower)
            cls_h, dist_h = classify_r(mean_r_h) if not np.isnan(mean_r_h) else ("?", 0)
            out.append(f"\n  (b) Lower half ({len(rates_lower)} rates):")
            out.append(f"      <r> = {mean_r_h:.4f}  -> {cls_h} "
                       f"(dist={dist_h:.3f})")
            r_half_table[N] = mean_r_h
        else:
            r_half_table[N] = float('nan')

        # ---- Phase 5: Chiral GUE test ----
        evals_c = evals + Sg
        rates_c = -evals_c.real
        pos = rates_c[rates_c > 1e-10]
        neg = -rates_c[rates_c < -1e-10]
        if len(pos) > 0 and len(neg) > 0:
            pos_s = np.sort(pos)
            neg_s = np.sort(neg)
            n_pairs = min(len(pos_s), len(neg_s))
            pair_err = np.mean(np.abs(pos_s[:n_pairs] - neg_s[:n_pairs]))
            out.append(f"\n  Chiral symmetry (centered at Sg={Sg}):")
            out.append(f"      +/- pairs: {n_pairs}, "
                       f"mean error: {pair_err:.2e}")

        # ---- Phase 7: Complex plane ----
        nonzero = evals[np.abs(evals) > 1e-10]
        out.append(f"\n  Complex plane:")
        out.append(f"      Fraction Re<0: "
                   f"{np.sum(nonzero.real < 0)/len(nonzero):.4f}")
        out.append(f"      |lambda| max: {np.max(np.abs(evals)):.4f}")

    # ---- Phase 4: Sector decomposition (Python, N=2-5) ----
    out.append(f"\n{'='*70}")
    out.append(f"PHASE 4: Sector decomposition by XY-weight (Python, N=2-5)")
    out.append(f"  Within each sector, all modes have the SAME decay rate.")
    out.append(f"  NNSD is computed from FREQUENCIES (imaginary parts),")
    out.append(f"  not rates (real parts), because rates are degenerate.")
    out.append(f"{'='*70}")

    for N in range(2, 6):
        J_val = 1.0
        H = build_H_chain(N, J_val)
        L = build_L_zdeph(H, gamma)
        evals_py = np.linalg.eigvals(L)
        rates_py = -evals_py.real
        freqs_py = np.abs(evals_py.imag)

        out.append(f"\n  N={N}:")

        # Classify modes by rate into weight sectors
        tol = 0.3 * gamma
        for w_target in range(N + 1):
            expected_rate = 2 * w_target * gamma

            if w_target == 0:
                mask = rates_py < 1e-8
            else:
                mask = np.abs(rates_py - expected_rate) < tol

            sector_freqs = freqs_py[mask]
            # Only oscillating modes (nonzero frequency)
            sector_freqs = sector_freqs[sector_freqs > 1e-6]
            # Unique frequencies (remove conjugate duplicates)
            sector_freqs = np.sort(np.unique(np.round(sector_freqs, 8)))

            if len(sector_freqs) < 5:
                n_total = int(np.sum(mask))
                if n_total > 0:
                    out.append(f"    w={w_target}: {n_total} modes, "
                               f"{len(sector_freqs)} unique freq "
                               f"(too few for NNSD)")
                continue

            sp_sector = nnsd(unfold(sector_freqs))
            sp_sector = sp_sector[sp_sector > 0.001]

            if len(sector_freqs) >= 5:
                mean_r_s, _ = spacing_ratios(sector_freqs)
                cls_s, _ = classify_r(mean_r_s) if not np.isnan(mean_r_s) else ("?", 0)
                out.append(f"    w={w_target}: {len(sector_freqs)} unique freq, "
                           f"<r>={mean_r_s:.3f} -> {cls_s}")
            else:
                out.append(f"    w={w_target}: {len(sector_freqs)} unique freq "
                           f"(too few for spacing ratio)")

    # ---- Phase 6: Summary table ----
    out.append(f"\n{'='*70}")
    out.append(f"SUMMARY")
    out.append(f"{'='*70}")

    out.append(f"\n  Spacing ratio <r> (0.386=Poisson, 0.536=GOE, 0.603=GUE):")
    out.append(f"  {'N':>3} | {'<r> all':>8} | {'<r> half':>9} | "
               f"{'Class (all)':>12} | {'Class (half)':>13}")
    out.append(f"  " + "-" * 55)
    for N in sorted(r_all_table.keys()):
        ra = r_all_table.get(N, float('nan'))
        rh = r_half_table.get(N, float('nan'))
        ca, _ = classify_r(ra) if not np.isnan(ra) else ("?", 0)
        ch, _ = classify_r(rh) if not np.isnan(rh) else ("?", 0)
        out.append(f"  {N:3d} | {ra:8.4f} | {rh:9.4f} | {ca:>12} | {ch:>13}")

    text = "\n".join(out)
    print(text)

    out_path = RESULTS / "rmt_analysis.txt"
    out_path.write_text(text, encoding="utf-8")
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
