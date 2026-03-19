"""
N -> infinity Limit of the Palindromic Spectrum
================================================
Seven-section analysis:
  1. Rate density histograms (full eigenvalues for small N, L_D for all N)
  2. XOR mode fraction (N=2-20)
  3. Weight sector sizes and palindromic counting
  4. Bandwidth scaling under Hamiltonian perturbation
  5. Standing wave frequency density
  6. Past/future boundary width
  7. Z-deph vs depol Gaussian comparison

Script: simulations/n_infinity_analysis.py
Output: simulations/results/n_infinity_analysis.txt
"""
import numpy as np
from math import comb, factorial
from datetime import datetime
import time

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\n_infinity_analysis.txt"
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
# MAIN
# ============================================================
gamma = 0.05

log("=" * 90)
log("N -> INFINITY: Thermodynamic Limit of the Palindromic Spectrum")
log(f"Date: {datetime.now()}")
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
    rates, counts = ld_rate_distribution(N, gamma)
    total = np.sum(counts)
    # Moments
    probs = counts / total
    mean = np.sum(probs * rates)
    var = np.sum(probs * (rates - mean) ** 2)
    std = np.sqrt(var)
    skew = np.sum(probs * ((rates - mean) / std) ** 3) if std > 0 else 0
    kurt = np.sum(probs * ((rates - mean) / std) ** 4) - 3 if std > 0 else 0
    log(f"  {N:>4}  {mean:>8.4f}  {std:>8.4f}  {N + 1:>10}  {skew:>10.6f}  "
        f"{kurt:>10.6f}  {total:>10}")

log(f"\n  Gaussian prediction: mean = N*gamma, std = gamma*sqrt(N)")
log(f"  Skewness = 0 (exact, binomial p=1/2)")
log(f"  Kurtosis = -2/N (exact, approaches 0 = Gaussian)")

# Full eigenvalues for N = 3-5
log(f"\n  Full Liouvillian eigenvalues:")
full_eig_data = {}
for N in [3, 4, 5]:
    t0 = time.time()
    H = build_H_chain(N)
    L = build_L(H, gamma, N)
    evals = np.linalg.eigvals(L)
    dt = time.time() - t0
    rates = -np.real(evals)
    Sg = N * gamma

    # Moments of rate distribution
    mean = np.mean(rates)
    std_r = np.std(rates)
    skew_r = np.mean(((rates - mean) / std_r) ** 3) if std_r > 0 else 0
    kurt_r = np.mean(((rates - mean) / std_r) ** 4) - 3 if std_r > 0 else 0

    # Theoretical
    mean_th = Sg
    std_th = gamma * np.sqrt(N)
    kurt_th = -2.0 / N

    # Palindrome check
    n_paired = 0
    for k in range(len(evals)):
        target = -(evals[k] + 2 * Sg)
        if np.min(np.abs(evals - target)) < 1e-8:
            n_paired += 1
    palin_pct = 100 * n_paired / len(evals)

    full_eig_data[N] = (evals, rates)

    log(f"\n    N={N} ({4 ** N} eigenvalues, {dt:.2f}s):")
    log(f"      Rate range: [{np.min(rates):.6f}, {np.max(rates):.6f}]")
    log(f"      Mean: {mean:.6f} (theory: {mean_th:.6f})")
    log(f"      Std:  {std_r:.6f} (theory: {std_th:.6f})")
    log(f"      Skew: {skew_r:.6f} (theory: 0)")
    log(f"      Kurt: {kurt_r:.6f} (theory: {kurt_th:.6f})")
    log(f"      Palindromic: {palin_pct:.1f}%")

    # Distribution comparison: count eigenvalues in bins matching L_D levels
    ld_rates, ld_counts = ld_rate_distribution(N, gamma)
    log(f"\n      Rate distribution by L_D weight sector:")
    log(f"      {'w':>4}  {'L_D rate':>10}  {'L_D count':>10}  "
        f"{'eig min':>10}  {'eig max':>10}  {'bandwidth':>10}")
    log(f"      {'-' * 56}")
    for wi, (r, c) in enumerate(zip(ld_rates, ld_counts)):
        # Find eigenvalues near this L_D rate
        tol = gamma * 0.9  # half the gap between adjacent levels
        mask = np.abs(rates - r) < tol
        if np.sum(mask) > 0:
            eig_min = np.min(rates[mask])
            eig_max = np.max(rates[mask])
            bw = eig_max - eig_min
        else:
            eig_min = eig_max = r
            bw = 0
        log(f"      {wi:>4}  {r:>10.6f}  {c:>10}  "
            f"{eig_min:>10.6f}  {eig_max:>10.6f}  {bw:>10.6f}")


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
log(f"  GHZ fragility (100% XOR projection) becomes irrelevant at large N")


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
log(f"  'Pure past' vanishes exponentially. Most states are mixed.")


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
    evals, rates = full_eig_data[N]
    Sg = N * gamma
    ld_rates, ld_counts = ld_rate_distribution(N, gamma)
    num_evals = len(rates)

    log(f"\n  N = {N} ({num_evals} eigenvalues):")

    # Count distinct rates
    unique_rates = len(set(f"{r:.8f}" for r in rates))
    log(f"    Distinct rates: {unique_rates}")
    log(f"    L_D levels: {N + 1}")

    # For each L_D level, compute the band width
    total_bw = 0
    n_bands = 0
    for wi, (r, c) in enumerate(zip(ld_rates, ld_counts)):
        tol = gamma * 0.9
        mask = np.abs(rates - r) < tol
        n_in = np.sum(mask)
        if n_in > 0:
            bw = np.max(rates[mask]) - np.min(rates[mask])
            total_bw += bw
            if bw > 1e-10:
                n_bands += 1

    avg_bw = total_bw / (N + 1)
    log(f"    Bands with nonzero width: {n_bands}/{N + 1}")
    log(f"    Average band width: {avg_bw:.6f} = {avg_bw / gamma:.4f}*gamma")
    log(f"    Total rate range: [{np.min(rates):.6f}, {np.max(rates):.6f}]")
    log(f"    Span: {np.max(rates) - np.min(rates):.6f} = {(np.max(rates) - np.min(rates)) / gamma:.2f}*gamma")

log(f"\n  Bandwidth scaling prediction:")
log(f"    Boundary rates are topology-independent: min=0, max=2N*gamma")
log(f"    Dynamic range: [2*gamma, 2*(N-1)*gamma]")
log(f"    Bandwidth = 2*(N-2)*gamma, linear in N")
log(f"    As N grows, bands broaden and merge into a continuum")


# ############################################################
# SECTION 5: Standing wave frequency density
# ############################################################
log()
log("=" * 90)
log("SECTION 5: Standing wave frequency density")
log("  Number of distinct oscillation frequencies vs N")
log("=" * 90)

for N in [3, 4, 5]:
    if N not in full_eig_data:
        continue
    evals, _ = full_eig_data[N]
    freqs = np.abs(np.imag(evals))
    # Distinct nonzero frequencies
    nonzero = freqs[freqs > 1e-6]
    unique_freqs = sorted(set(f"{fr:.6f}" for fr in nonzero))

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
        log(f"    Range: [{np.min(nonzero):.6f}, {np.max(nonzero):.6f}]")

log(f"\n  Pattern: frequency count grows rapidly with N")
log(f"  N=3: discrete harmonics. N=5+: approaching continuous spectrum.")


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
log(f"  The 'sharp' classical/quantum boundary at small N becomes a smooth")
log(f"  Gaussian spread at large N. Most states are near w = N/2: half-classical,")
log(f"  half-quantum. Pure past (w=0) and pure future (w=N) become exponentially rare.")


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
log(f"  THERMODYNAMIC SUMMARY")
log(f"  {'=' * 70}")
log(f"""
  Z-DEPHASING at large N:
    Rate distribution: Gaussian, center = N*gamma, width = gamma*sqrt(N)
    Palindrome: TRIVIALLY satisfied (symmetric Gaussian)
    Past/future: boundary blurs, most states are half-classical
    XOR fraction: vanishes as (N+1)/4^N
    Standing wave: discrete -> continuous spectrum
    The palindrome is automatic but non-trivial (L_H preserves it)

  DEPOLARIZING at large N:
    Rate distribution: Gaussian, center = N*gamma, width = gamma*sqrt(N/3)
    Palindrome: EXPONENTIALLY broken
    Counting ratio 3^(N-2w) means sectors differ by exp(N)
    Past (w=0): 1 string. Future (w=N): 3^N strings.
    The mirror deficit grows exponentially with system size.

  THE ANSWER to 'does palindrome become trivially true?':
    The L_D part: YES, Gaussian symmetry makes it automatic.
    The L_H part: NO, this is the non-trivial content.
    Pi anti-commuting with [H,.] ensures L_H does not break the
    L_D palindrome. Without Pi, L_H could shift eigenvalues
    asymmetrically. The proof guarantees it never does.
    This constraint holds at every N, including N -> infinity.
""")


# ############################################################
# DONE
# ############################################################
log()
log("=" * 90)
log("ANALYSIS COMPLETE")
log(f"Date: {datetime.now()}")
log("=" * 90)
f.close()
print(f"\n>>> Results written to {OUT}")
