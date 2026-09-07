#!/usr/bin/env python3
"""
Spectral Form Factor of the Palindromic Liouvillian
=====================================================
Phase 1: Load eigenvalues from RMT CSVs (N=2-7)
Phase 2: Compute dissipative + frequency SFF
Phase 3: Identify palindromic modulation
Phase 4: Summarize reached windows with the independent-phase reference
Phase 5: Band-resolved SFF (N=3-5, via Python eigendecomposition)
Phase 6: Report raw density scale, palindromic period and sampled exceedance
Phase 7: Connection to previous results

Script: simulations/spectral_form_factor.py
Output: simulations/results/spectral_form_factor.txt
"""

import numpy as np
from scipy.linalg import eigvals
from itertools import product as iproduct
import os, sys, time as clock
from sff_window_summary import (sff_frequency, raw_multiset_density_scale,
                                independent_phase_reference, summarize_windows, format_sample)

OUT_PATH = os.environ.get("RCPSI_SFF_OUTPUT_PATH") or os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "results", "spectral_form_factor.txt")
CSV_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ========================================================================
# Eigenvalue loading
# ========================================================================
def load_eigenvalues(N):
    """Load complex eigenvalues from RMT CSV (tab-separated Re, Im)."""
    path = os.path.join(CSV_DIR, f"rmt_eigenvalues_N{N}.csv")
    data = np.loadtxt(path, delimiter='\t', skiprows=1)
    return data[:, 0] + 1j * data[:, 1]


# ========================================================================
# SFF computation
# ========================================================================
def sff_dissipative(eigenvalues, t_arr):
    """K_diss(t) = (1/N²) Σ_{j,k} exp(i(λ_j - λ_k*)t)."""
    N = len(eigenvalues)
    K = np.zeros(len(t_arr))
    # K_diss = |Σ_k exp(i λ_k t)|² when λ are eigenvalues of non-Hermitian
    # Actually: Σ_{j,k} exp(i(λ_j - conj(λ_k))t) = |Σ_k exp(i λ_k t)|²
    # if we define the sum carefully. Let's compute directly.
    for ti, t in enumerate(t_arr):
        phases = np.exp(1j * eigenvalues * t)
        # Σ_{j,k} exp(i(λ_j - λ_k*)t) = (Σ_j exp(iλ_j t))(Σ_k exp(-iλ_k* t))
        #                                = (Σ exp(iλt)) * conj(Σ exp(iλ*t))... hmm
        # More precisely: Σ_{j,k} exp(i(λ_j - conj(λ_k))t)
        # = (Σ_j exp(i λ_j t)) × (Σ_k exp(-i conj(λ_k) t))
        # = (Σ_j exp(i λ_j t)) × conj(Σ_k exp(i conj(λ_k) t))... no.
        # exp(-i conj(λ_k) t) = conj(exp(i λ_k t)) when t is real.
        # Wait: conj(exp(i λ t)) = exp(-i conj(λ) t) = exp(i (-conj(λ)) t)
        # So Σ_k exp(-i conj(λ_k) t) = conj(Σ_k exp(i λ_k t))
        # Therefore K_diss = |Σ_k exp(i λ_k t)|² / N²
        K[ti] = np.abs(np.sum(phases))**2 / N**2
    return K


def sff_connected(eigenvalues, t_arr):
    """Connected SFF: K_c(t) = K(t) - |<exp(iλt)>|² (remove disconnected)."""
    N = len(eigenvalues)
    K_full = np.zeros(len(t_arr))
    K_disc = np.zeros(len(t_arr))
    freqs = eigenvalues.imag
    for ti, t in enumerate(t_arr):
        phases = np.exp(1j * freqs * t)
        s = np.sum(phases)
        K_full[ti] = np.abs(s)**2 / N**2
        K_disc[ti] = np.abs(np.mean(phases))**2
    return K_full - K_disc


# ========================================================================
# Pauli infrastructure for sector analysis
# ========================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]


def build_liouvillian_pauli_with_sectors(N, gamma=0.05, J=1.0):
    """Build Liouvillian in Pauli basis, return L and XY-weights."""
    dim = 2**N
    num = 4**N
    all_idx = list(iproduct(range(4), repeat=N))

    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in idx[1:]:
            m = np.kron(m, PAULIS[k])
        pmats.append(m)
    pstack = np.array(pmats)

    H = np.zeros((dim, dim), dtype=complex)
    ops = [sx, sy, sz]
    for i in range(N - 1):
        for P in ops:
            opi = np.eye(1, dtype=complex)
            for k in range(N):
                opi = np.kron(opi, P if k == i else I2)
            opj = np.eye(1, dtype=complex)
            for k in range(N):
                opj = np.kron(opj, P if k == i + 1 else I2)
            H += J * (opi @ opj)

    L = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pstack[b] - pstack[b] @ H)
        L[:, b] = np.einsum('aij,ji->a', pstack, comm) / dim

    for a, idx in enumerate(all_idx):
        rate = 0.0
        for site in range(N):
            if idx[site] in (1, 2):
                rate += 2 * gamma
        L[a, a] -= rate

    # XY-weight per Pauli string
    xy_weights = np.array([sum(1 for x in idx if x in (1, 2))
                           for idx in all_idx])

    return L, xy_weights


# ========================================================================
log("=" * 72)
log("SPECTRAL FORM FACTOR OF THE PALINDROMIC LIOUVILLIAN")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)

gamma = 0.05
J = 1.0

# ========================================================================
# PHASE 1: LOAD EIGENVALUES
# ========================================================================
log()
log("=" * 72)
log("PHASE 1: LOAD EIGENVALUES FROM RMT CSVs")
log("=" * 72)
log()

all_evals = {}
for N in range(2, 8):
    try:
        ev = load_eigenvalues(N)
        all_evals[N] = ev
        n_osc = np.sum(np.abs(ev.imag) > 1e-10)
        log(f"  N={N}: {len(ev)} eigenvalues loaded, {n_osc} oscillating")
    except FileNotFoundError:
        log(f"  N={N}: CSV not found, skipping")


# ========================================================================
# PHASE 2: COMPUTE SFF
# ========================================================================
log()
log("=" * 72)
log("PHASE 2: SPECTRAL FORM FACTOR")
log("=" * 72)
log()

sff_results = {}

for N in sorted(all_evals.keys()):
    ev = all_evals[N]
    n_ev = len(ev)

    # Center: μ = λ + Σγ so palindromic pairs are at ±μ
    sigma_gamma = N * gamma
    mu = ev + sigma_gamma

    raw_scale = raw_multiset_density_scale(ev)
    raw_mean_gap = 2 * np.pi / raw_scale

    # Palindromic time: t_Pi = 2π / ω_min where ω_min = 4J(1-cos(π/N))
    omega_min = 4 * J * (1 - np.cos(np.pi / N))
    t_Pi = 2 * np.pi / omega_min

    # Bounded grid; the raw density scale is descriptive, not a physical time.
    t_max = min(3 * raw_scale, 50 * t_Pi, 200)
    n_t = min(2000, max(500, int(t_max * 50)))
    t_arr = np.linspace(0, t_max, n_t)

    t0 = clock.time()
    K_freq = sff_frequency(ev, t_arr)
    K_diss = sff_dissipative(ev, t_arr)
    elapsed = clock.time() - t0

    sff_results[N] = {
        't': t_arr, 'K_freq': K_freq, 'K_diss': K_diss,
        'raw_mean_gap': raw_mean_gap, 'raw_scale': raw_scale, 't_Pi': t_Pi,
        'omega_min': omega_min, 'n_ev': n_ev
    }

    log(f"  N={N} ({n_ev} eigenvalues, {elapsed:.1f}s):")
    log(f"    Raw multiset mean adjacent gap = {raw_mean_gap:.4f}")
    log(f"    Raw multiset density scale = {raw_scale:.2f}")
    log(f"    ω_min (slowest)  = {omega_min:.4f}")
    log(f"    t_Π (palindromic) = {t_Pi:.2f}")
    log(f"    K_freq range: [{np.min(K_freq):.4e}, {np.max(K_freq):.4f}]")
    log(f"    K_diss range: [{np.min(K_diss):.4e}, {np.max(K_diss):.4f}]")
    log()


# ========================================================================
# PHASE 3: PALINDROMIC MODULATION
# ========================================================================
log()
log("=" * 72)
log("PHASE 3: PALINDROMIC MODULATION")
log("=" * 72)
log()
log("  Each palindromic pair (μ, -μ) contributes cos(Im(μ)·t).")
log("  Expected: periodic modulation with period 2π/ω_min.")
log()

for N in sorted(sff_results.keys()):
    r = sff_results[N]
    t_arr = r['t']
    K = r['K_freq']
    t_Pi = r['t_Pi']
    omega_min = r['omega_min']

    # Remove DC and compute FFT to find modulation frequencies
    K_detrend = K - np.mean(K)
    dt = t_arr[1] - t_arr[0]
    fft_vals = np.abs(np.fft.rfft(K_detrend))
    fft_freqs = np.fft.rfftfreq(len(K_detrend), d=dt) * 2 * np.pi  # angular freq

    # Find peaks in FFT (excluding DC)
    peak_idx = np.argsort(fft_vals[1:])[::-1][:5] + 1
    peak_freqs = fft_freqs[peak_idx]
    peak_amps = fft_vals[peak_idx]

    # Find closest match among the five largest non-DC FFT amplitudes.
    if len(peak_freqs) > 0:
        best_match = np.argmin(np.abs(peak_freqs - omega_min))
        match_freq = peak_freqs[best_match]
        match_err = abs(match_freq - omega_min) / omega_min * 100

        # Also check 2*omega_min (second harmonic from |cos|²)
        best_2 = np.argmin(np.abs(peak_freqs - 2 * omega_min))
        match_2freq = peak_freqs[best_2]
        match_2err = abs(match_2freq - 2 * omega_min) / (2 * omega_min) * 100

        log(f"  N={N}:")
        log(f"    ω_min (predicted) = {omega_min:.4f}")
        log(f"    Top 3 FFT peaks:   {['%.4f' % f for f in peak_freqs[:3]]}")
        log(f"    Best match to ω_min:  {match_freq:.4f} (err {match_err:.1f}%)")
        log(f"    Best match to 2ω_min: {match_2freq:.4f} (err {match_2err:.1f}%)")

        # Modulation visibility: peak amplitude / mean
        if np.mean(K) > 0:
            visibility = peak_amps[0] / np.sum(fft_vals[1:]) if np.sum(fft_vals[1:]) > 0 else 0
            log(f"    Modulation visibility: {visibility:.4f}")
        log()


# ========================================================================
# PHASE 4: COMPARISON WITH STANDARD CLASSES
# ========================================================================
log()
log("=" * 72)
log("PHASE 4: RAW FREQUENCY SFF — SAMPLED WINDOWS")
log("=" * 72)
log()
log("  This is the raw non-unfolded frequency SFF, normalized by M^2.")
log("  The independent-phase reference 1/M is not a measured plateau or class verdict.")
log("  Degenerate frequencies can raise the actual long-time average above this reference.")
log("  The raw multiset density scale is 2*pi / mean adjacent gap(sorted abs nonzero frequencies).")
log("  It retains multiplicities (abs frequency > 1e-10), so it is multiplicity-dependent.")
log("  Below/intermediate/beyond are descriptive bins, not physical time regimes or ramp/plateau evidence.")
log("  Only the reached window is assessed; absent bins are not sampled, not zero.")
log("  Dissipative overflow/NaN remains unresolved and is not interpreted here.")
log()

for N in sorted(sff_results.keys()):
    if N < 3:
        continue
    r = sff_results[N]
    t_arr = r['t']
    K = r['K_freq']
    raw_scale = r['raw_scale']

    windows = summarize_windows(t_arr, K, raw_scale)
    K_below, K_intermediate, K_beyond = (windows[name] for name in ("below", "intermediate", "beyond"))
    slope, K_min = windows["slope"], windows["minimum"]

    log(f"  N={N}:")
    log(f"    Reached time: {t_arr[-1]:.2f}; raw scale = {raw_scale:.2f}; independent-phase reference 1/M = {independent_phase_reference(r['n_ev']):.8g}")
    log(f"    Below (t < 0.1 raw scale):  <K> = {format_sample(K_below)}")
    log(f"    Intermediate (0.1-1.0 raw scale):    <K> = {format_sample(K_intermediate)}   slope = {format_sample(slope)}")
    log(f"    Beyond (t > raw scale):        <K> = {format_sample(K_beyond)}")
    minimum_text = "not sampled" if K_min is None else f"{K_min:.4e}"
    log(f"    Min K (below bin):  {minimum_text}")

    if K_below is None or K_min is None or slope is None:
        log("    Heuristic flags: not classifiable (below bin or intermediate-bin slope not sampled)")
    elif K_below > 0.5 and abs(slope) < 0.5:
        log(f"    Heuristic flags: below-bin mean > 0.5 and abs(intermediate-bin slope) < 0.5; not a class verdict")
    elif K_min < 0.1 and slope > 0.3:
        log(f"    Heuristic flags: below-bin minimum < 0.1 and intermediate-bin slope > 0.3; not a class verdict")
    else:
        log(f"    Heuristic flags: neither below-bin-mean/slope nor minimum/slope condition; not a class verdict")
    log()


# ========================================================================
# PHASE 5: BAND-RESOLVED SFF
# ========================================================================
log()
log("=" * 72)
log("PHASE 5: BAND-RESOLVED SFF (N=3, 4, 5)")
log("=" * 72)
log()

for N in [3, 4, 5]:
    t0 = clock.time()
    L, xy_w = build_liouvillian_pauli_with_sectors(N, gamma, J)
    ev_all = eigvals(L)

    # Get unique weights
    weights = sorted(set(xy_w))
    sigma_gamma = N * gamma

    # Map eigenvalues to sectors by diagonalizing L and checking
    # which Pauli basis vectors contribute to each eigenvector
    # Simpler approach: compute eigenvalues of L restricted to each sector
    # But sectors are mixed by H (w -> w±2). So we use the full spectrum
    # and classify by the DIAGONAL rate structure.

    # Group eigenvalues by approximate decay rate band
    rates = -ev_all.real
    freqs = ev_all.imag

    log(f"  N={N} ({4**N} eigenvalues, {clock.time()-t0:.1f}s):")

    # For each light-content band w, the expected decay rate is ~2wγ.  The bin is
    # a BAND, not a weight sector: the Absorption Theorem gives
    # rate = 2*gamma*<n_XY>, so |rate - 2wγ| < band_width admits every mode whose
    # AVERAGE light content is near w, pure or mixed.
    for w in weights:
        if w == 0:
            continue
        rate_center = 2 * w * gamma
        band_width = gamma  # approximate
        in_sector = np.abs(rates - rate_center) < band_width
        sector_ev = ev_all[in_sector]

        if len(sector_ev) < 4:
            continue

        # SFF for this sector
        omega_min_w = 4 * J * (1 - np.cos(np.pi / N))
        t_max_w = min(20 * 2 * np.pi / omega_min_w, 100)
        t_w = np.linspace(0.01, t_max_w, 500)
        K_w = sff_frequency(sector_ev, t_w)

        # Characterize
        K_mean = np.mean(K_w)
        K_std = np.std(K_w)
        K_min = np.min(K_w)

        log(f"    w={w}: {len(sector_ev)} eigenvalues,"
            f" <K>={K_mean:.4f}, std={K_std:.4f}, min={K_min:.4e}")

    log()


# ========================================================================
# PHASE 6: TIMESCALES
# ========================================================================
log()
log("=" * 72)
log("PHASE 6: TIMESCALES")
log("=" * 72)
log()

log(f"  {'N':>3}  {'t_Π':>8}  {'raw scale':>10}  {'t_Π/scale':>10}  {'ω_min':>8}  {'raw gap':>8}")
log(f"  {'─'*55}")

for N in sorted(sff_results.keys()):
    r = sff_results[N]
    ratio = r['t_Pi'] / r['raw_scale'] if r['raw_scale'] > 0 else 0
    log(f"  {N:>3}  {r['t_Pi']:>8.2f}  {r['raw_scale']:>10.2f}"
        f"  {ratio:>10.4f}  {r['omega_min']:>8.4f}  {r['raw_mean_gap']:>8.4f}")

log()
log("  t_Π = palindromic time = 2π/ω_min (period of slowest mode)")
log("  Raw multiset density scale = 2*pi / raw mean adjacent gap; multiplicity-dependent.")
log("  Its ratio to t_Π is descriptive and does not define a physical time-regime boundary.")
log()

# First sampled exceedance of 1.5x the second-half-of-current-grid mean.
# This is not a physical time-scale estimate or a Poisson/physical asymptotic baseline.
log("  Heuristic: first sampled exceedance of 1.5x the second-half-of-current-grid mean:")
log("  This is not a physical time-scale estimate and not a Poisson/physical asymptotic baseline.")
for N in sorted(sff_results.keys()):
    if N < 3:
        continue
    r = sff_results[N]
    K = r['K_freq']
    t_arr = r['t']
    second_half_mean = np.mean(K[len(K)//2:]) if len(K) > 10 else 1
    threshold = 1.5 * second_half_mean
    above = np.where(K > threshold)[0]
    if len(above) > 0:
        first_sampled_exceedance_time = t_arr[above[0]]
        log(f"    N={N}: first sampled exceedance t ≈ {first_sampled_exceedance_time:.2f}  (t/raw scale = {first_sampled_exceedance_time/r['raw_scale']:.4f})")
    else:
        log(f"    N={N}: no sampled exceedance detected")


# ========================================================================
# PHASE 7: CONNECTION TO PREVIOUS RESULTS
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 7: CONNECTION TO PREVIOUS RESULTS")
log("=" * 72)
log()
log("  Poisson/no-ramp behavior is compatible with integrability or block fragmentation;")
log("  it does not prove integrability. Read the sampled SFF separately from a class claim.")
log()
log("  PT analysis said: Pi is chiral (class AIII). The palindromic")
log("  modulation in the SFF is the TIME-DOMAIN signature of the")
log("  same spectral pairing that RMT sees in LEVEL STATISTICS.")
log()
log("  Topo analysis said: geometric, not topological. The SFF is")
log("  independent of localization (it measures spectral correlations,")
log("  not spatial profiles).")
log()
log("  Analytical formulas: ω_min = 4J(1-cos(π/N)) (formula 2, k=1).")
log("  If the SFF modulation peak matches ω_min, that confirms the")
log("  (0,1) coherence block's dispersion relation in the time domain.")


# ========================================================================
# SUMMARY
# ========================================================================
log()
log()
log("=" * 72)
log("SUMMARY")
log("=" * 72)
log()

# Determine dominant behavior
n_above_threshold = 0
n_below_threshold = 0
modulation_confirmed = False

for N in sorted(sff_results.keys()):
    if N < 3:
        continue
    r = sff_results[N]
    K = r['K_freq']
    t_arr = r['t']
    below = t_arr < 0.1 * r['raw_scale']
    if np.any(below) and np.mean(K[below]) > 0.3:
        n_above_threshold += 1
    else:
        n_below_threshold += 1

log(f"  Below-bin mean K_freq > 0.3 at {n_above_threshold}/{n_above_threshold+n_below_threshold} system sizes.")
log(f"  Below-bin mean K_freq <= 0.3 (or empty below bin) at {n_below_threshold}/{n_above_threshold+n_below_threshold} system sizes.")
log()
log("  The below bin is t < 0.1*raw scale; this threshold is not a spectral-class or no-ramp test.")
log("  The sampled curves and their modulation are the SFF evidence; no integrability theorem follows.")
log("  Palindromic modulation provides time-domain fingerprint of")
log("  the spectral pairing lambda <-> -(lambda + 2*Sigma_gamma).")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
