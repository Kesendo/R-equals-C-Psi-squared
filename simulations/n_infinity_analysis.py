"""
N -> infinity Limit of the Palindromic Spectrum
================================================
Seven-section analysis:
  1. Rate density histograms (full eigenvalues for small N, L_D for all N)
  2. XOR mode fraction (N=2-20)
  3. Weight sector sizes and palindromic counting
  4. Bandwidth scaling under Hamiltonian perturbation
  5. Sampled oscillation-frequency counts
  6. Past/future boundary width
  7. Z-deph vs depol Gaussian comparison

Script: simulations/n_infinity_analysis.py
Output: simulations/results/n_infinity_analysis.txt (or the path after --out)

Every number in the report is either exact (integer or rational arithmetic,
or a value the palindrome fixes: the zero rates, the endpoint 2N*gamma, the
rate mean = -Tr(L)/4^N, the vanishing rate skewness) or an eigensolver value
read through one stated error law, so the report is byte-identical whatever
the BLAS thread count or machine load:

  ERROR LAW (a measured constant, not a proven bound). np.linalg.eigvals is
  backward stable: it returns the exact spectrum of L + E with ||E|| of order
  eps * ||L||_F. An eigenvalue then moves by at most kappa_i * ||E||, kappa_i
  its condition number (1/|y_i^H x_i| for unit left/right eigenvectors), and
  L is not normal, so kappa_i > 1 is possible and no a-priori constant
  follows. What is measured, at N = 3, 4, 5 under 1, 3 and 24 OpenBLAS
  threads, in units err = eps * ||L||_F: every cluster spread (rates up to
  2.5, |Im| up to 7.7), every deviation of the zero cluster from 0 and of the
  endpoint cluster from 2N*gamma (up to 0.2), every palindrome pairing distance
  (up to 5.5); the smallest gap between distinct values is 7.2e-7, about
  2e7 err. So the effective kappa * c here is below 8. NOISE = 64 err groups
  values into clusters and is ASSERTED per cluster (spread, anchor deviation)
  and per pairing distance, not merely assumed; SIGNAL = 1e6 err is the floor
  a genuine gap must clear. A gap between the two, a cluster wider than NOISE,
  a printed value within NOISE of a rounding boundary, or a would-be "-0"
  digit string raises instead of printing. The per-N worst ratio to err is
  written to stderr (not the report, which must stay byte-stable).

  WHY 0 AND 2N*gamma ARE SNAPPED EXACTLY. In the Hilbert-Schmidt inner
  product, L_H = -i[H, .] is anti-Hermitian and the dephasing dissipator
  D = gamma * sum_k (Z_k . Z_k - 1) is Hermitian and negative semidefinite.
  If L X = lambda X with Re lambda = 0, then 0 = Re<X, L X> = <X, D X>, so
  D X = 0 and X is diagonal in the computational basis. For diagonal X the
  diagonal of [H, X] vanishes, while L X = -i[H, X] = lambda X is diagonal;
  hence lambda X = 0 and lambda = 0. So every rate is >= 0 and a rate equal
  to 0 is the eigenvalue 0 exactly. The F1 palindrome Pi L Pi^-1 = -L - 2N*gamma
  maps lambda to -lambda - 2N*gamma, carrying ker L onto the eigenvalue
  -2N*gamma exactly and making 2N*gamma the largest rate; its cluster is that
  eigenvalue exactly.
"""
import sys
import numpy as np
from fractions import Fraction
from pathlib import Path
from math import comb, factorial

OUT = Path(__file__).resolve().parent / "results" / "n_infinity_analysis.txt"
if "--out" in sys.argv:
    OUT = Path(sys.argv[sys.argv.index("--out") + 1])
f = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()


# ============================================================
# OPERATORS
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, s, N):
    d = 2 ** N
    ops = [I2] * N
    ops[s] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H_chain(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, N) @ site_op(P, i + 1, N)
    return H


def build_L(H, gamma, N):
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2, dtype=complex))
    return L


# ============================================================
# L_D RATE DISTRIBUTION (analytical, works for any N)
# ============================================================
def ld_rate_distribution(N, gamma):
    """Return (rates, counts) for L_D eigenvalues.
    Rate 2*gamma*w has multiplicity C(N,w)*2^N for w=0..N."""
    rates = []
    counts = []
    for w in range(N + 1):
        rates.append(2 * gamma * w)
        counts.append(comb(N, w) * (2 ** N))
    return np.array(rates), np.array(counts)


# ============================================================
# ERROR LAW: clusters and printing (see the module docstring)
# ============================================================
EPS = np.finfo(float).eps
NOISE_K = 64.0
SIGNAL_K = 1.0e6


def noise_of(L):
    return NOISE_K * EPS * np.linalg.norm(L)


def cluster_representatives(values, noise, exact_anchors=()):
    """Group values whose consecutive sorted gaps are <= noise; return, per
    input value, its cluster representative. A cluster within noise of an
    exact anchor (0, 2N*gamma) is represented by that anchor exactly;
    any other cluster by its mean. Raises if a gap falls between NOISE and
    SIGNAL (the law would not decide it)."""
    order = np.argsort(values, kind="stable")
    v = values[order]
    gaps = np.diff(v)
    signal = noise * (SIGNAL_K / NOISE_K)
    ambiguous = gaps[(gaps > noise) & (gaps < signal)]
    if ambiguous.size:
        raise AssertionError(f"gap {ambiguous.min():.3e} lies between NOISE "
                             f"{noise:.3e} and SIGNAL {signal:.3e}")
    reps_sorted = np.empty_like(v)
    starts = np.concatenate(([0], np.where(gaps > noise)[0] + 1, [v.size]))
    worst = 0.0
    for a, b in zip(starts[:-1], starts[1:]):
        block = v[a:b]
        # Consecutive gaps <= NOISE could chain into a wide cluster; the law
        # is on the whole cluster, so its full spread is asserted.
        spread = block[-1] - block[0]
        if spread > noise:
            raise AssertionError(f"cluster spread {spread:.3e} exceeds NOISE {noise:.3e}")
        worst = max(worst, spread)
        rep = float(np.mean(block))
        for anchor in exact_anchors:
            dev = max(abs(block[0] - anchor), abs(block[-1] - anchor))
            if dev <= noise:
                rep = float(anchor)
                worst = max(worst, dev)
        reps_sorted[a:b] = rep
    reps = np.empty_like(values)
    reps[order] = reps_sorted
    return reps, worst


def fmt(x, digits, noise=0.0):
    """Fixed-point string of x that cannot depend on rounding noise: raises
    if x lies within noise of a rounding boundary or would print as -0."""
    x = float(x)
    scaled = abs(x) * 10 ** digits
    if abs((scaled % 1.0) - 0.5) * 10.0 ** (-digits) <= noise:
        raise AssertionError(f"{x!r} within {noise:.3e} of a rounding boundary")
    s = f"{x:.{digits}f}"
    if s.startswith("-") and set(s[1:]) <= set("0."):
        raise AssertionError(f"{x!r} would print as a signed zero")
    return s


def band_mask(rates, level, noise):
    """Cluster representatives within 0.9*gamma of an L_D level (the window
    the report has always used); raises if one sits within noise of the edge."""
    edge = gamma * 0.9
    dist = np.abs(rates - level)
    if np.any(np.abs(dist - edge) <= noise):
        raise AssertionError(f"a rate sits within noise of the band edge at {level}")
    return dist < edge


# ============================================================
# MAIN
# ============================================================
gamma = 0.05

log("=" * 90)
log("N -> INFINITY: Counting asymptotics and finite-N spectral samples")
log(f"gamma = {gamma}, Heisenberg chain, Z-dephasing")
log("=" * 90)


# ############################################################
# SECTION 1: Rate density
# ############################################################
log()
log("=" * 90)
log("SECTION 1: Rate density histograms")
log("  L_D rates (exact binomial) vs full eigenvalue rates")
log("  Gaussian prediction: center = N*gamma, width = gamma*sqrt(N)")
log("=" * 90)

# Analytical L_D distribution for N up to 20
log(f"\n  L_D rate distribution (analytical):")
log(f"  {'N':>4}  {'center':>8}  {'std':>8}  {'n_levels':>10}  {'skewness':>10}  "
    f"{'kurtosis':>10}  {'total':>10}")
log(f"  {'-' * 62}")

for N in range(3, 21):
    Sg = N * gamma
    # Exact moments in weight space (rate = 2*gamma*w, a positive scale, so
    # skewness and kurtosis are the weight-space ones; rational arithmetic).
    counts = [comb(N, w) * 2 ** N for w in range(N + 1)]
    total = sum(counts)
    m1 = Fraction(sum(c * w for w, c in enumerate(counts)), total)
    m2 = sum(Fraction(c, total) * (w - m1) ** 2 for w, c in enumerate(counts))
    m3 = sum(Fraction(c, total) * (w - m1) ** 3 for w, c in enumerate(counts))
    m4 = sum(Fraction(c, total) * (w - m1) ** 4 for w, c in enumerate(counts))
    if m3 != 0:
        raise AssertionError(f"N={N}: binomial third central moment must be 0")
    kurt = m4 / m2 ** 2 - 3
    if kurt != Fraction(-2, N):
        raise AssertionError(f"N={N}: binomial excess kurtosis must be -2/N")
    mean = 2 * gamma * float(m1)
    std = 2 * gamma * float(m2) ** 0.5
    skew = 0.0
    log(f"  {N:>4}  {mean:>8.4f}  {std:>8.4f}  {N + 1:>10}  {skew:>10.6f}  "
        f"{float(kurt):>10.6f}  {total:>10}")

log(f"\n  Gaussian prediction: mean = N*gamma, std = gamma*sqrt(N)")
log(f"  Skewness = 0 (exact, binomial p=1/2)")
log(f"  Kurtosis = -2/N (exact, approaches 0 = Gaussian)")

# Full eigenvalues for N = 3-5
log(f"\n  Full Liouvillian eigenvalues:")
full_eig_data = {}
for N in [3, 4, 5]:
    H = build_H_chain(N)
    L = build_L(H, gamma, N)
    evals = np.linalg.eigvals(L)
    noise = noise_of(L)
    Sg = N * gamma
    # Rates grouped by the error law; the zero cluster (steady states) is 0
    # exactly and its palindrome partner 2N*gamma exactly.
    raw_rates = -np.real(evals)
    rates, worst_rate = cluster_representatives(raw_rates, noise,
                                                exact_anchors=(0.0, 2 * Sg))
    if not (np.min(rates) == 0.0 and np.max(rates) == 2 * Sg):
        raise AssertionError(f"N={N}: rate range must be exactly [0, 2N*gamma]")
    # |Im| grouped the same way; the zero cluster is exactly 0.
    freq_reps, worst_freq = cluster_representatives(np.abs(np.imag(evals)), noise,
                                         exact_anchors=(0.0,))

    # Moments of rate distribution. The mean is exact: sum of eigenvalues =
    # Tr(L); H is traceless in it and the dissipator diagonal is
    # -2*gamma*popcount(i XOR j), so -Tr(L)/4^N = 2*gamma * <popcount> = N*gamma.
    d = 2 ** N
    pop_sum = sum(bin(i ^ j).count("1") for i in range(d) for j in range(d))
    mean = float(2 * Fraction(gamma) * Fraction(pop_sum, d * d))
    if Fraction(pop_sum, d * d) != Fraction(N, 2):
        raise AssertionError(f"N={N}: mean popcount must be N/2")
    if abs(np.mean(raw_rates) - mean) > noise:
        raise AssertionError(f"N={N}: eigenvalue mean leaves the trace identity")
    std_r = np.std(raw_rates)
    # The palindrome makes the rate multiset symmetric about N*gamma, so the
    # skewness is exactly 0; the eigensolver value must sit inside the law.
    skew_meas = np.mean(((raw_rates - Sg) / std_r) ** 3)
    if abs(skew_meas) > 3 * noise / std_r:
        raise AssertionError(f"N={N}: rate skewness {skew_meas:.3e} exceeds the law")
    kurt_r = np.mean(((raw_rates - Sg) / std_r) ** 4) - 3

    # Theoretical
    mean_th = Sg
    std_th = gamma * np.sqrt(N)
    kurt_th = -2.0 / N

    # Palindrome check: each eigenvalue's distance to the reflected spectrum
    # is either noise (paired) or a genuine gap; the law decides, not a cut.
    n_paired = 0
    worst_pair = 0.0
    for k in range(len(evals)):
        target = -(evals[k] + 2 * Sg)
        dist = np.min(np.abs(evals - target))
        if noise < dist < noise * (SIGNAL_K / NOISE_K):
            raise AssertionError(f"N={N}: pairing distance {dist:.3e} undecided")
        if dist <= noise:
            n_paired += 1
            worst_pair = max(worst_pair, dist)
    palin_pct = 100 * n_paired / len(evals)
    err = noise / NOISE_K
    print(f"[error law] N={N}: worst / (eps*||L||_F): rate clusters "
          f"{worst_rate / err:.2f}, |Im| clusters {worst_freq / err:.2f}, "
          f"pairing {worst_pair / err:.2f} (asserted <= {NOISE_K:.0f})",
          file=sys.stderr, flush=True)

    full_eig_data[N] = (evals, rates, freq_reps, noise)

    log(f"\n    N={N} ({4 ** N} eigenvalues):")
    log(f"      Rate range: [{fmt(np.min(rates), 6)}, {fmt(np.max(rates), 6)}]")
    log(f"      Mean: {fmt(mean, 6)} (theory: {mean_th:.6f})")
    log(f"      Std:  {fmt(std_r, 6, noise)} (theory: {std_th:.6f})")
    log(f"      Skew: {fmt(0.0, 6)} (theory: 0)")
    log(f"      Kurt: {fmt(kurt_r, 6, 4 * noise / std_r)} (theory: {kurt_th:.6f})")
    log(f"      Palindromic: {palin_pct:.1f}%")

    # Distribution comparison: count eigenvalues in bins matching L_D levels
    ld_rates, ld_counts = ld_rate_distribution(N, gamma)
    log(f"\n      Rate distribution by L_D weight sector:")
    log(f"      {'w':>4}  {'L_D rate':>10}  {'L_D count':>10}  "
        f"{'eig min':>10}  {'eig max':>10}  {'bandwidth':>10}")
    log(f"      {'-' * 56}")
    for wi, (r, c) in enumerate(zip(ld_rates, ld_counts)):
        mask = band_mask(rates, r, noise)
        if np.sum(mask) > 0:
            eig_min = np.min(rates[mask])
            eig_max = np.max(rates[mask])
            bw = eig_max - eig_min
        else:
            eig_min = eig_max = r
            bw = 0
        log(f"      {wi:>4}  {r:>10.6f}  {c:>10}  "
            f"{fmt(eig_min, 6, noise):>10}  {fmt(eig_max, 6, noise):>10}  "
            f"{fmt(bw, 6, 2 * noise):>10}")


# ############################################################
# SECTION 2: XOR mode fraction
# ############################################################
log()
log("=" * 90)
log("SECTION 2: XOR mode fraction")
log("  XOR modes = N+1, total modes = 4^N")
log("=" * 90)

log(f"\n  {'N':>4}  {'XOR':>6}  {'Total':>12}  {'Fraction':>12}  {'Status':>12}")
log(f"  {'-' * 50}")

n_below_1pct = None
n_below_01pct = None

for N in range(2, 21):
    xor = N + 1
    total = 4 ** N
    frac = xor / total
    pct = 100 * frac
    log(f"  {N:>4}  {xor:>6}  {total:>12}  {pct:>11.6f}%  "
        f"{'< 1%' if pct < 1 else '':>12}")
    if pct < 1 and n_below_1pct is None:
        n_below_1pct = N
    if pct < 0.01 and n_below_01pct is None:
        n_below_01pct = N

log(f"\n  XOR fraction drops below 1% at N = {n_below_1pct}")
log(f"  XOR fraction drops below 0.01% at N = {n_below_01pct}")
log(f"  Scaling: (N+1)/4^N -> 0 exponentially")
log(f"  This operator-space fraction is not a prepared-state XOR probability")
if (n_below_1pct, n_below_01pct) != (5, 9):
    raise AssertionError("F23 threshold regression: expected first N values 5 and 9")
if 100 * (9 / 65536) <= 0.01:
    raise AssertionError("N=8 F23 row must remain above the 0.01% threshold")
log("  PASS exact thresholds: first <1% at N=5; first <0.01% at N=9")


# ############################################################
# SECTION 3: Weight sector sizes
# ############################################################
log()
log("=" * 90)
log("SECTION 3: Weight sector sizes (palindromic counting proof)")
log("  Z-deph: count(w) = C(N,w) * 2^N, partner count(N-w) = same")
log("=" * 90)

for N in [3, 4, 6, 8, 10, 12]:
    log(f"\n  N = {N}:")
    log(f"  {'w':>4}  {'count':>12}  {'N-w':>4}  {'partner':>12}  {'equal':>6}")
    log(f"  {'-' * 42}")
    for w in range(N + 1):
        c = comb(N, w) * (2 ** N)
        pw = N - w
        pc = comb(N, pw) * (2 ** N)
        log(f"  {w:>4}  {c:>12}  {pw:>4}  {pc:>12}  {'YES' if c == pc else 'NO':>6}")

for N in [3, 8, 12]:
    endpoint_fraction = (2 ** N) / (4 ** N)
    if endpoint_fraction != 2 ** (-N):
        raise AssertionError(f"N={N}: bare endpoint fraction must equal 2^-N")
log("\n  PASS bare endpoint fraction: 2^N / 4^N = 2^-N (N=3,8,12)")

# "Past is tiny" ratio
log(f"\n  'Past is tiny' ratio: count(0) / count(N/2)")
log(f"  count(0) = 2^N, count(N/2) = C(N,N/2) * 2^N")
log(f"  Ratio = 1 / C(N,N/2)")
log(f"\n  {'N':>6}  {'C(N,N/2)':>14}  {'ratio':>14}  {'log10(ratio)':>14}")
log(f"  {'-' * 52}")

for N in [4, 6, 8, 10, 20, 50, 100]:
    half = N // 2
    c = comb(N, half)
    ratio = 1.0 / c
    log(f"  {N:>6}  {c:>14}  {ratio:>14.2e}  {np.log10(ratio):>14.2f}")

log(f"\n  At N=100: the w=0 sector (all-classical, pure past) is 10^-29 times")
log(f"  smaller than the w=50 sector (half-classical, half-quantum).")
log(f"  The endpoint's fraction of the bare-dissipator Pauli-string count vanishes exponentially.")


# ############################################################
# SECTION 4: Bandwidth scaling
# ############################################################
log()
log("=" * 90)
log("SECTION 4: Bandwidth scaling under Hamiltonian perturbation")
log("  How much does L_H spread the L_D rate levels?")
log("=" * 90)

for N in [3, 4, 5]:
    if N not in full_eig_data:
        continue
    evals, rates, _, noise = full_eig_data[N]
    Sg = N * gamma
    ld_rates, ld_counts = ld_rate_distribution(N, gamma)
    num_evals = len(rates)

    log(f"\n  N = {N} ({num_evals} eigenvalues):")

    # Count distinct rates: one per cluster of the error law
    unique_rates = len(np.unique(rates))
    log(f"    Distinct rates: {unique_rates}")
    log(f"    L_D levels: {N + 1}")

    # For each L_D level, compute the band width (exactly 0 for a band that
    # is a single cluster, since clusters carry one representative)
    total_bw = 0
    n_bands = 0
    for wi, (r, c) in enumerate(zip(ld_rates, ld_counts)):
        mask = band_mask(rates, r, noise)
        n_in = np.sum(mask)
        if n_in > 0:
            bw = np.max(rates[mask]) - np.min(rates[mask])
            total_bw += bw
            if bw > 0:
                n_bands += 1

    avg_bw = total_bw / (N + 1)
    span = np.max(rates) - np.min(rates)
    log(f"    Bands with nonzero width: {n_bands}/{N + 1}")
    log(f"    Average band width: {fmt(avg_bw, 6, 2 * noise)} = "
        f"{fmt(avg_bw / gamma, 4, 2 * noise / gamma)}*gamma")
    log(f"    Total rate range: [{fmt(np.min(rates), 6)}, {fmt(np.max(rates), 6)}]")
    log(f"    Span: {fmt(span, 6)} = {fmt(span / gamma, 2)}*gamma")

log(f"\n  Finite-N rate-range summary:")
log(f"    Boundary rates are topology-independent: min=0, max=2N*gamma")
log(f"    Dynamic range: [2*gamma, 2*(N-1)*gamma]")
log(f"    Bandwidth = 2*(N-2)*gamma, linear in N")
log(f"    These N=3..5 samples do not establish a continuum limit")


# ############################################################
# SECTION 5: Standing wave frequency density
# ############################################################
log()
log("=" * 90)
log("SECTION 5: Sampled oscillation-frequency counts")
log("  Number of distinct oscillation frequencies vs N")
log("=" * 90)

for N in [3, 4, 5]:
    if N not in full_eig_data:
        continue
    evals, _, freqs, noise = full_eig_data[N]
    # Distinct nonzero frequencies: clusters of the error law other than the
    # exact-zero one, listed in numeric order
    nonzero = freqs[freqs > 0]
    unique_freqs = [fmt(fr, 6, noise) for fr in np.unique(nonzero)]
    if len(set(unique_freqs)) != len(unique_freqs):
        raise AssertionError(f"N={N}: two distinct frequencies print alike")

    log(f"\n  N = {N}:")
    log(f"    Total eigenvalues: {len(evals)}")
    log(f"    With Im != 0: {len(nonzero)}")
    log(f"    Distinct |omega|: {len(unique_freqs)}")

    if len(unique_freqs) <= 20:
        log(f"    Frequencies: {', '.join(unique_freqs)}")
    else:
        log(f"    First 10: {', '.join(unique_freqs[:10])}")
        log(f"    Last 5:  {', '.join(unique_freqs[-5:])}")

    # Frequency range
    if nonzero.size > 0:
        log(f"    Range: [{fmt(np.min(nonzero), 6, noise)}, {fmt(np.max(nonzero), 6, noise)}]")

log(f"\n  Pattern: frequency count grows rapidly with N")
log(f"  No continuous limiting spectrum follows from the N=3..5 counts.")


# ############################################################
# SECTION 6: Past/future boundary width
# ############################################################
log()
log("=" * 90)
log("SECTION 6: Past/future boundary width")
log("  What fraction of Pauli strings lie near w = N/2?")
log("  Boundary layer: |w - N/2| <= sqrt(N)")
log("=" * 90)

log(f"\n  {'N':>6}  {'fraction':>10}  {'2sigma_gauss':>12}  {'note':>20}")
log(f"  {'-' * 50}")

for N in [3, 4, 6, 10, 20, 50, 100, 500, 1000, 10000]:
    half = N / 2.0
    width = np.sqrt(N)  # boundary layer: |w - N/2| <= sqrt(N) = 2*std
    if N <= 200:
        # Exact computation for small N
        in_boundary = 0
        total = 0
        for w in range(N + 1):
            c = comb(N, w)  # factor 2^N cancels in ratio
            total += c
            if abs(w - half) <= width:
                in_boundary += c
        frac = in_boundary / total
    else:
        # Gaussian approximation: boundary = 2*sigma of Bin(N, 0.5)
        # P(|w - N/2| <= sqrt(N)) = P(|Z| <= 2) where Z is standard normal
        from scipy.stats import norm as norm_dist
        frac = norm_dist.cdf(2) - norm_dist.cdf(-2)  # ~0.9545
    gauss_2sig = 0.9545
    note = "exact" if N <= 200 else "Gaussian approx"
    log(f"  {N:>6}  {frac:>10.4f}  {gauss_2sig:>12.4f}  {note:>20}")

log(f"\n  The fraction of strings within sqrt(N) of the midpoint approaches ~0.954")
log(f"  (the 2-sigma fraction of a Gaussian), independent of N.")
log(f"  This is the bare-dissipator Pauli-string counting measure, not a state ensemble.")
log(f"  Its endpoint fractions become exponentially small while the central")
log(f"  weight window carries the stated asymptotic fraction.")


# ############################################################
# SECTION 7: Z-deph vs depol Gaussian comparison
# ############################################################
log()
log("=" * 90)
log("SECTION 7: Z-deph vs depol Gaussian comparison")
log("  Z-deph: Bin(N, 1/2) in weight -> symmetric around N/2")
log("  Depol: Bin(N, 3/4) in weight -> peaked at 3N/4, asymmetric")
log("=" * 90)

log(f"\n  Distributions (in weight space w = 0..N):")
log(f"  {'':>10}  {'Z-deph':>20}  {'Depol':>20}")
log(f"  {'':>10}  {'p = 1/2':>20}  {'p = 3/4':>20}")
log(f"  {'-' * 54}")
log(f"  {'Mean':>10}  {'N/2':>20}  {'3N/4':>20}")
log(f"  {'Variance':>10}  {'N/4':>20}  {'3N/16':>20}")
log(f"  {'Std':>10}  {'sqrt(N)/2':>20}  {'sqrt(3N)/4':>20}")
log(f"  {'Symmetric':>10}  {'YES (p=1/2)':>20}  {'NO (p=3/4)':>20}")
log(f"  {'Palindrome':>10}  {'YES':>20}  {'NO':>20}")

log(f"\n  Palindrome center (in rate space):")
log(f"  Z-deph: center = N*gamma (= mean of distribution). Symmetric. PALINDROMIC.")
log(f"  Depol: center = N*gamma. But distribution peaks at 3N/4 (in weight) =")
log(f"         (4gamma/3)*(3N/4) = N*gamma (in rate). Wait: the RATE center matches")
log(f"         the distribution peak! But the counting is still asymmetric.")

# Show the counting mismatch for depol
log(f"\n  Depol counting mismatch: count(w) / count(N-w) = 3^(N-2w)")
log(f"\n  {'N':>4}  {'w':>4}  {'N-w':>4}  {'count(w)':>12}  {'count(N-w)':>12}  {'ratio':>10}")
log(f"  {'-' * 50}")

for N in [4, 10, 20]:
    for w in [0, 1, N // 4, N // 2]:
        cw = comb(N, w) * (3 ** w)
        cnw = comb(N, N - w) * (3 ** (N - w))
        ratio = cw / cnw if cnw > 0 else float('inf')
        log(f"  {N:>4}  {w:>4}  {N - w:>4}  {cw:>12.2e}  {cnw:>12.2e}  {ratio:>10.2e}")

log(f"\n  At N=20, w=0 vs w=20: ratio = 1 / 3^20 = {1 / 3 ** 20:.2e}")
log(f"  The 'past' (w=0) has 1 string. The 'future' (w=N) has 3^N strings.")
log(f"  Under depol, the future is exponentially larger than the past.")

# Key comparison
log(f"\n  {'=' * 70}")
log(f"  SCOPE SUMMARY")
log(f"  {'=' * 70}")
log(f"""
  BARE Z-DEPHASING PAULI-STRING COUNT at large N:
    Weight-count distribution: Gaussian asymptotics around N/2
    Count symmetry: exact under w -> N-w
    One endpoint weight sector: 2^N strings, fraction 2^-N
    No state ensemble, standing wave, or continuous spectral limit follows.

  DISTINCT INTERACTING F23 OBJECT (in its stated chain scope):
    Endpoint eigenspace dimension fraction: (N+1)/4^N

  BARE DEPOLARIZING PAULI-STRING COUNT at large N:
    Weight-count distribution: Gaussian asymptotics around 3N/4
    Count symmetry: exponentially imbalanced
    Counting ratio 3^(N-2w) means sectors differ by exp(N)
    Past (w=0): 1 string. Future (w=N): 3^N strings.
    The mirror deficit grows exponentially with system size.

  SPECTRAL SCOPE:
    The Pi relation proves the palindrome for each finite N in its stated
    Hamiltonian/dephasing scope. These computations do not construct or prove
    convergence to an infinite-volume Liouvillian or spectral measure.
""")


# ############################################################
# DONE
# ############################################################
log()
log("=" * 90)
log("ANALYSIS COMPLETE")
log("=" * 90)
f.close()
print(f"\n>>> Results written to {OUT}")
