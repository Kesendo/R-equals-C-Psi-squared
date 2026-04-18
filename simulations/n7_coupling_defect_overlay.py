#!/usr/bin/env python3
"""
N = 7 coupling-defect overlay.

Run experiment A (J_ij = 1 on all 6 bonds) and five experiments B (J_{0,1}
= J_mod, others = 1) on an open XY chain with uniform Z-dephasing gamma_0
= 0.05 per site. The palindrome Pi L Pi^{-1} = -L - 2 Sigma_gamma I is
algebraically intact in both: J multiplies XX + YY, which is Pi-invariant.
What differs is the EXPRESSION of the palindromic structure in the
spectrum and dynamics, because the eigenvectors of the palindromic pairs
depend on J. We probe that difference operationally via time evolution.

Initial state: rho_0 = |phi><phi| with phi = (|vac> + |psi_1>)/sqrt(2),
where psi_1 is the F65 bonding mode for the uniform chain. In A, phi is
a natural "lens" (an H-eigenvector plus vacuum); in B it is not, so the
difference dynamics is especially clean.

For each experiment we RK4-propagate the 128 x 128 density matrix over T
= 80 (units of J^{-1}) with 400 steps, recording site-resolved single-
qubit purity, pair-CPsi and pair mutual information for pairs (0, k),
and the global L1(rho). For each J_mod we then form Delta(t) = Obs_B -
Obs_A and analyse:

    * per-site arrival time at the 1 %% threshold,
    * max |Delta| per site vs J_mod,
    * integrated |Delta| per site,
    * FFT of Delta(t) at the middle site (3) vs J_mod,
    * symmetry of Delta(J_mod = 0.5) vs Delta(J_mod = 2.0).

Conventions (repo-standard): big-endian qubit labelling, site 0 is the
most-significant bit of the 2^N state index; site_op(op, 0, N) acts as
kron(op, I_{N-1}).

Date: 2026-04-17
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N = 7
GAMMA_0 = 0.05
T_FINAL = 80.0
N_STEPS = 400
DT = T_FINAL / N_STEPS   # 0.2
J_UNIFORM = 1.0
J_MOD_VALUES = [0.5, 0.8, 1.2, 1.5, 2.0]

ARRIVAL_THRESH = 0.01   # 1 %% absolute threshold for arrival-time detection

RESULTS_DIR = Path(__file__).parent / "results" / "n7_coupling_defect_overlay"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = RESULTS_DIR / "run_log.txt"


# ---------------------------------------------------------------------------
# Paulis and operator helpers
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    """Operator op acting on chain qubit `site`, identity elsewhere.
    Big-endian: site 0 is the outermost factor."""
    factors = [I2] * N
    factors[site] = op
    return kron_chain(*factors)


def build_H_XY(J_list, N):
    """H = sum_i (J_i / 2) (X_i X_{i+1} + Y_i Y_{i+1})."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        J_i = J_list[i]
        H += (J_i / 2.0) * (site_op(X, i, N) @ site_op(X, i + 1, N) +
                            site_op(Y, i, N) @ site_op(Y, i + 1, N))
    return H


def build_hamming_matrix(N):
    """H[a, b] = popcount(a XOR b). Used to apply full-site Z-dephasing."""
    d = 2**N
    idx = np.arange(d, dtype=np.uint32)
    xor = idx[:, None] ^ idx[None, :]
    h = np.zeros((d, d), dtype=np.int32)
    for i in range(N):
        h += ((xor >> i) & 1).astype(np.int32)
    return h


# ---------------------------------------------------------------------------
# Initial states
# ---------------------------------------------------------------------------
def single_excitation_mode(N, k=1):
    """F65 bonding mode: psi_k = sqrt(2/(N+1)) sum_i sin(pi k (i+1)/(N+1)) |1_i>.
    Big-endian: qubit i excited => state index 2^(N-1-i)."""
    psi = np.zeros(2**N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
        psi[2**(N - 1 - i)] = amp
    return psi


def vacuum(N):
    v = np.zeros(2**N, dtype=complex)
    v[0] = 1.0
    return v


def initial_rho_psi1_bonding(N):
    """rho_0 = |phi><phi| with phi = (|vac> + |psi_1>)/sqrt(2)."""
    phi = vacuum(N) + single_excitation_mode(N, 1)
    phi /= np.linalg.norm(phi)
    return np.outer(phi, phi.conj())


def initial_rho_plus_N(N):
    """rho_0 = |+>^N = (H|0>)^N tensor."""
    phi = np.ones(2**N, dtype=complex) / np.sqrt(2**N)
    return np.outer(phi, phi.conj())


def initial_rho_bell_endpoints(N):
    """rho_0 = |bell><bell|, bell = (|0...0> + |1 at q_0 and q_{N-1}>)/sqrt(2).
    Big-endian: q_0 is bit N-1, q_{N-1} is bit 0."""
    phi = np.zeros(2**N, dtype=complex)
    phi[0] = 1.0 / np.sqrt(2)
    phi[2**(N - 1) + 1] = 1.0 / np.sqrt(2)
    return np.outer(phi, phi.conj())


# ---------------------------------------------------------------------------
# Lindblad RHS + RK4
# ---------------------------------------------------------------------------
def lindblad_rhs(rho, H, hamming, gamma_0):
    """d rho / dt = -i [H, rho] - 2 gamma_0 * hamming .* rho.

    Derivation: sum_i gamma_0 (Z_i rho Z_i - rho). In the computational
    basis, Z_i rho Z_i element-(a,b) is z_i[a] z_i[b] rho[a,b] with
    z_i[x] = (-1)^{bit_i(x)}. Summed over i, the Hadamard prefactor is
    N - 2 popcount(a XOR b); the -rho part subtracts N, leaving -2 *
    popcount(a XOR b).
    """
    commutator = H @ rho - rho @ H
    dephasing = -2.0 * gamma_0 * hamming * rho
    return -1j * commutator + dephasing


def rk4_step(rho, H, hamming, gamma_0, dt):
    k1 = lindblad_rhs(rho,                H, hamming, gamma_0)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, hamming, gamma_0)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, hamming, gamma_0)
    k4 = lindblad_rhs(rho +       dt * k3, H, hamming, gamma_0)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


# ---------------------------------------------------------------------------
# Partial trace / observables (big-endian: axis q = site q)
# ---------------------------------------------------------------------------
_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def reduced_single(rho, i, N):
    row = [_LETTERS[q] for q in range(N)]
    col = [_LETTERS[q + N] for q in range(N)]
    for q in range(N):
        if q != i:
            col[q] = row[q]   # trace qubit q
    spec = ''.join(row) + ''.join(col) + '->' + row[i] + col[i]
    return np.einsum(spec, rho.reshape([2] * (2 * N)))


def reduced_pair(rho, i, j, N):
    assert i < j
    row = [_LETTERS[q] for q in range(N)]
    col = [_LETTERS[q + N] for q in range(N)]
    for q in range(N):
        if q != i and q != j:
            col[q] = row[q]
    spec = (''.join(row) + ''.join(col) + '->'
            + row[i] + row[j] + col[i] + col[j])
    r = np.einsum(spec, rho.reshape([2] * (2 * N)))
    return r.reshape(4, 4)


def site_purity_all(rho, N):
    out = np.empty(N)
    for i in range(N):
        r = reduced_single(rho, i, N)
        out[i] = float(np.real(np.trace(r @ r)))
    return out


def pair_cpsi_from_rho_AB(rho_AB):
    C = float(np.real(np.trace(rho_AB @ rho_AB)))
    diag = np.diag(np.diag(rho_AB))
    L1 = float(np.sum(np.abs(rho_AB - diag)))
    return C * L1 / 3.0


def von_neumann_entropy(r):
    """S(r) = -sum lambda_i log2 lambda_i, lambda_i eigenvalues of r."""
    rh = (r + r.conj().T) / 2.0
    eig = np.linalg.eigvalsh(rh)
    eig = np.clip(eig, 1e-15, None)
    return float(-np.sum(eig * np.log2(eig)))


def pair_cpsi_and_mi_for_refsite(rho, ref, N):
    cpsi = np.empty(N - 1)
    mi = np.empty(N - 1)
    rho_ref = reduced_single(rho, ref, N)
    s_ref = von_neumann_entropy(rho_ref)
    slot = 0
    for k in range(N):
        if k == ref:
            continue
        rho_k = reduced_single(rho, k, N)
        i, j = (ref, k) if ref < k else (k, ref)
        rho_AB = reduced_pair(rho, i, j, N)
        cpsi[slot] = pair_cpsi_from_rho_AB(rho_AB)
        s_k = von_neumann_entropy(rho_k)
        s_AB = von_neumann_entropy(rho_AB)
        mi[slot] = s_ref + s_k - s_AB
        slot += 1
    return cpsi, mi


def global_l1(rho):
    diag = np.diag(np.diag(rho))
    return float(np.sum(np.abs(rho - diag)))


# ---------------------------------------------------------------------------
# Experiment runner
# ---------------------------------------------------------------------------
def run_experiment(label, J_list, initial_rho_fn, hamming, times, ref_site=0):
    H = build_H_XY(J_list, N)
    rho = initial_rho_fn(N)

    T = len(times)
    purity = np.empty((T, N))
    cpsi_0k = np.empty((T, N - 1))
    mi_0k = np.empty((T, N - 1))
    l1 = np.empty(T)
    trace_err = np.empty(T)
    herm = np.empty(T)

    def record(step, rho):
        purity[step] = site_purity_all(rho, N)
        cpsi_0k[step], mi_0k[step] = pair_cpsi_and_mi_for_refsite(
            rho, ref_site, N)
        l1[step] = global_l1(rho)
        trace_err[step] = abs(float(np.real(np.trace(rho))) - 1.0)
        herm[step] = float(np.linalg.norm(rho - rho.conj().T))

    record(0, rho)
    for step in range(1, T):
        rho = rk4_step(rho, H, hamming, GAMMA_0, DT)
        record(step, rho)
    return dict(label=label, purity=purity, cpsi_0k=cpsi_0k, mi_0k=mi_0k,
                l1=l1, trace_err=trace_err, herm=herm)


# ---------------------------------------------------------------------------
# Overlay analysis helpers
# ---------------------------------------------------------------------------
def arrival_time(times, delta_1d, threshold):
    """First time at which |delta| exceeds threshold. NaN if never."""
    idx = np.where(np.abs(delta_1d) > threshold)[0]
    return float(times[idx[0]]) if len(idx) > 0 else float('nan')


def arrival_time_per_site(times, delta_2d, threshold):
    return np.array([arrival_time(times, delta_2d[:, i], threshold)
                     for i in range(delta_2d.shape[1])])


def max_abs_per_site(delta_2d):
    return np.max(np.abs(delta_2d), axis=0)


def integrated_abs_per_site(times, delta_2d):
    return np.trapezoid(np.abs(delta_2d), x=times, axis=0)


def fft_peak_freqs(signal, dt, k_top=5):
    n = len(signal)
    sig = signal - float(np.mean(signal))
    fft = np.fft.rfft(sig)
    freqs = np.fft.rfftfreq(n, d=dt)
    amps = np.abs(fft)
    order = np.argsort(amps)[::-1][:k_top]
    return [(float(freqs[i]), float(amps[i])) for i in order]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log_lines = []
    def log(msg=""):
        print(msg, flush=True)
        log_lines.append(msg)

    log("=" * 72)
    log("N = 7 COUPLING DEFECT OVERLAY")
    log("=" * 72)
    log(f"Model: open XY chain, N = {N}, J_uniform = {J_UNIFORM}")
    log(f"       Z-dephasing on all {N} sites, gamma_0 = {GAMMA_0}")
    log(f"       Sigma_gamma = {N * GAMMA_0:.3f}, above F18 fold threshold")
    log(f"Evolution: Lindblad RK4, dt = {DT}, T = {T_FINAL}, steps = {N_STEPS}")
    log(f"Defect:    J_{{0,1}} = J_mod for J_mod in {J_MOD_VALUES}")
    log(f"Initial state: rho_0 = |phi><phi|, "
        f"phi = (|vac> + |psi_1>)/sqrt(2)  (psi_1 = F65 bonding mode)")
    log(f"Arrival threshold: |Delta Obs| > {ARRIVAL_THRESH} absolute")
    log()

    hamming = build_hamming_matrix(N)
    times = np.arange(N_STEPS + 1) * DT

    # --- Experiment A: uniform J --------------------------------------
    t0 = time.time()
    J_A = [J_UNIFORM] * (N - 1)
    result_A = run_experiment("A (uniform)", J_A,
                              initial_rho_psi1_bonding, hamming, times)
    log(f"Experiment A (uniform J): {time.time() - t0:.2f} s   "
        f"|tr(rho)-1|_max = {result_A['trace_err'].max():.2e},   "
        f"||rho - rho*||_max = {result_A['herm'].max():.2e}")

    # --- Experiments B ------------------------------------------------
    results_B = {}
    for J_mod in J_MOD_VALUES:
        t0 = time.time()
        J_B = [J_mod] + [J_UNIFORM] * (N - 2)
        r = run_experiment(f"B(J_mod={J_mod})", J_B,
                           initial_rho_psi1_bonding, hamming, times)
        results_B[J_mod] = r
        log(f"Experiment B (J_mod = {J_mod}): {time.time() - t0:.2f} s   "
            f"|tr(rho)-1|_max = {r['trace_err'].max():.2e},   "
            f"||rho - rho*||_max = {r['herm'].max():.2e}")

    # --- Save raw time series -----------------------------------------
    np.save(RESULTS_DIR / "times.npy", times)
    np.savez(RESULTS_DIR / "experiment_A.npz",
             purity=result_A['purity'], cpsi_0k=result_A['cpsi_0k'],
             mi_0k=result_A['mi_0k'], l1=result_A['l1'])
    for J_mod, r in results_B.items():
        np.savez(RESULTS_DIR / f"experiment_B_{J_mod}.npz",
                 purity=r['purity'], cpsi_0k=r['cpsi_0k'],
                 mi_0k=r['mi_0k'], l1=r['l1'])
    log()
    log(f"Raw time series saved to {RESULTS_DIR}/")

    # --- Baseline shapes ----------------------------------------------
    log()
    log("Baseline (experiment A) snapshots:")
    for snap_t in [0.0, 5.0, 20.0, 40.0, 80.0]:
        idx = int(round(snap_t / DT))
        idx = min(idx, len(times) - 1)
        log(f"  t = {times[idx]:5.2f}:  "
            f"purity per site = [{', '.join(f'{p:.4f}' for p in result_A['purity'][idx])}]")
        log(f"               cpsi(0,k) k=1..6 = "
            f"[{', '.join(f'{p:.4f}' for p in result_A['cpsi_0k'][idx])}]")
        log(f"               L1(rho)          = {result_A['l1'][idx]:.4f}")

    # --- Overlay analysis: site purity --------------------------------
    log()
    log("=" * 72)
    log("OVERLAY ANALYSIS: site-resolved purity")
    log("=" * 72)
    log()
    log("Arrival time per site (first t at which |Delta P_i(t)| > "
        f"{ARRIVAL_THRESH})")
    log(f"  {'J_mod':>7}  " + "  ".join(f"site {i}".rjust(6) for i in range(N)))
    site_arrivals = {}
    for J_mod in J_MOD_VALUES:
        d = results_B[J_mod]['purity'] - result_A['purity']
        arr = arrival_time_per_site(times, d, ARRIVAL_THRESH)
        site_arrivals[J_mod] = arr
        cells = [f"{a:6.2f}" if not np.isnan(a) else "   nan" for a in arr]
        log(f"  {J_mod:>7}  " + "  ".join(cells))
    log()
    log("Max |Delta P_i| per site over T = [0, 80]")
    log(f"  {'J_mod':>7}  " + "  ".join(f"site {i}".rjust(8) for i in range(N)))
    site_amps = {}
    for J_mod in J_MOD_VALUES:
        d = results_B[J_mod]['purity'] - result_A['purity']
        amp = max_abs_per_site(d)
        site_amps[J_mod] = amp
        log(f"  {J_mod:>7}  " + "  ".join(f"{a:8.5f}" for a in amp))
    log()
    log("Integrated |Delta P_i| per site  (trapz over t in [0, 80])")
    log(f"  {'J_mod':>7}  " + "  ".join(f"site {i}".rjust(8) for i in range(N)))
    site_integ = {}
    for J_mod in J_MOD_VALUES:
        d = results_B[J_mod]['purity'] - result_A['purity']
        integ = integrated_abs_per_site(times, d)
        site_integ[J_mod] = integ
        log(f"  {J_mod:>7}  " + "  ".join(f"{a:8.4f}" for a in integ))

    # --- Overlay analysis: pair-CPsi(0, k) ----------------------------
    log()
    log("=" * 72)
    log("OVERLAY ANALYSIS: pair-CPsi(0, k)")
    log("=" * 72)
    log()
    log("Max |Delta CPsi_{0,k}| per pair k = 1..6")
    log(f"  {'J_mod':>7}  " +
        "  ".join(f"k = {k}".rjust(8) for k in range(1, N)))
    for J_mod in J_MOD_VALUES:
        d = results_B[J_mod]['cpsi_0k'] - result_A['cpsi_0k']
        amp = max_abs_per_site(d)
        log(f"  {J_mod:>7}  " + "  ".join(f"{a:8.5f}" for a in amp))
    log()
    log("Integrated |Delta CPsi_{0,k}|  (trapz)")
    log(f"  {'J_mod':>7}  " +
        "  ".join(f"k = {k}".rjust(8) for k in range(1, N)))
    for J_mod in J_MOD_VALUES:
        d = results_B[J_mod]['cpsi_0k'] - result_A['cpsi_0k']
        integ = integrated_abs_per_site(times, d)
        log(f"  {J_mod:>7}  " + "  ".join(f"{a:8.4f}" for a in integ))

    # --- Overlay: pair mutual info -----------------------------------
    log()
    log("=" * 72)
    log("OVERLAY ANALYSIS: pair MI(Q_0 : Q_k)")
    log("=" * 72)
    log()
    log("Max |Delta MI_{0,k}| per pair k = 1..6")
    log(f"  {'J_mod':>7}  " +
        "  ".join(f"k = {k}".rjust(8) for k in range(1, N)))
    for J_mod in J_MOD_VALUES:
        d = results_B[J_mod]['mi_0k'] - result_A['mi_0k']
        amp = max_abs_per_site(d)
        log(f"  {J_mod:>7}  " + "  ".join(f"{a:8.5f}" for a in amp))

    # --- Overlay: global L1 ------------------------------------------
    log()
    log("=" * 72)
    log("OVERLAY ANALYSIS: global L1(rho)")
    log("=" * 72)
    log()
    log(f"  {'J_mod':>7}  {'max|DL1|':>10}  {'t_argmax':>10}  {'integ':>10}")
    for J_mod in J_MOD_VALUES:
        d = results_B[J_mod]['l1'] - result_A['l1']
        m = float(np.max(np.abs(d)))
        t_m = float(times[int(np.argmax(np.abs(d)))])
        integ = float(np.trapezoid(np.abs(d), x=times))
        log(f"  {J_mod:>7}  {m:>10.5f}  {t_m:>10.2f}  {integ:>10.4f}")

    # --- FFT at middle site for Delta(site purity 3) -----------------
    log()
    log("=" * 72)
    log("FFT of Delta P_3(t) (middle site)")
    log("=" * 72)
    log()
    i_mid = 3
    for J_mod in J_MOD_VALUES:
        dp = results_B[J_mod]['purity'][:, i_mid] - result_A['purity'][:, i_mid]
        peaks = fft_peak_freqs(dp, DT, k_top=5)
        log(f"  J_mod = {J_mod}: top-5 FFT peaks of Delta P_3(t)")
        for f, a in peaks:
            period = 1.0 / f if f > 0 else float('inf')
            log(f"    f = {f:8.4f}  period = {period:8.2f}  amp = {a:10.3f}")

    # --- Symmetry check: Delta(J_mod=0.5) vs Delta(J_mod=2.0) ---------
    log()
    log("=" * 72)
    log("SYMMETRY CHECK: Delta at J_mod = 0.5 vs J_mod = 2.0")
    log("=" * 72)
    log()
    log("Per-site Pearson correlation of (Delta P_i at 0.5) vs (Delta P_i at 2.0)")
    log("  and norm ratio ||Delta(2.0)|| / ||Delta(0.5)||  per site")
    d_half = results_B[0.5]['purity'] - result_A['purity']
    d_two = results_B[2.0]['purity'] - result_A['purity']
    log(f"  {'site':>4}  {'Pearson':>10}  {'||2.0||/||0.5||':>16}  "
        f"{'max|0.5|':>10}  {'max|2.0|':>10}")
    for site in range(N):
        a = d_half[:, site]
        b = d_two[:, site]
        if np.std(a) < 1e-15 or np.std(b) < 1e-15:
            corr = float('nan')
        else:
            corr = float(np.corrcoef(a, b)[0, 1])
        ratio = float(np.linalg.norm(b) / (np.linalg.norm(a) + 1e-15))
        log(f"  {site:>4}  {corr:>+10.4f}  {ratio:>16.4f}  "
            f"{np.max(np.abs(a)):>10.5f}  {np.max(np.abs(b)):>10.5f}")

    log()
    log("=" * 72)
    log("Runs complete.")
    log("=" * 72)

    (LOG_PATH).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"\nRun log saved to: {LOG_PATH}")


if __name__ == "__main__":
    main()
