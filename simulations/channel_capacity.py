#!/usr/bin/env python3
"""
Formal Channel Capacity of the γ→Observables Bridge
======================================================
Step 1: Jacobian (response matrix) A = ∂y/∂γ
Step 2: SVD → singular values = channel gains
Step 3: Shannon capacity via waterfilling
Step 4: Capacity vs noise level curve
Step 5: Codebook verification (random codes, ML decoding)

Script:  simulations/channel_capacity.py
Output:  simulations/results/channel_capacity.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "channel_capacity.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 5; d = 32; d2 = 1024


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H(J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i) @ site_op(P, i + 1)
    return H


def build_L(H, gammas):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def ptrace_keep(rho, keep):
    keep = list(keep)
    trace_out = [q for q in range(N) if q not in keep]
    dims = [2] * N
    reshaped = rho.reshape(dims + dims)
    current_n = N
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape((d_k, d_k))


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d_r = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d_r - 1) if d_r > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


H_chain = build_H(J=1.0)

# Initial state
psi_init = plus
for _ in range(N - 1):
    psi_init = np.kron(psi_init, plus)
rho0 = np.outer(psi_init, psi_init.conj())


def extract_features(rho):
    """Extended feature vector (25 features)."""
    f = []
    for q in range(N):
        f.append(purity(ptrace_keep(rho, [q])))
    for i in range(N):
        for j in range(i + 1, N):
            f.append(cpsi(ptrace_keep(rho, [i, j])))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sz @ rho_q))))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sx @ rho_q))))
    return np.array(f)


def extract_timeseries(gammas, t_points):
    """Feature vector from multiple time points."""
    L = build_L(H_chain, gammas)
    f = []
    for t in t_points:
        rho = evolve(L, rho0, t)
        f.extend(extract_features(rho).tolist())
    return np.array(f)


# ====================================================================
# Step 1: Jacobian A = dy/dg
# ====================================================================

def compute_jacobian(gamma_ref, t_points, dg=1e-4):
    """Compute Jacobian of feature vector w.r.t. gamma perturbations."""
    y_ref = extract_timeseries(gamma_ref, t_points)
    n_feat = len(y_ref)
    A = np.zeros((n_feat, N))

    for site in range(N):
        gammas_pert = list(gamma_ref)
        gammas_pert[site] += dg
        y_pert = extract_timeseries(gammas_pert, t_points)
        A[:, site] = (y_pert - y_ref) / dg

    return A, y_ref


# ====================================================================
# Step 3: Waterfilling capacity
# ====================================================================

def waterfilling_capacity(singular_values, P_total, sigma_noise):
    """Shannon capacity via waterfilling on SVD channels.

    C = sum_i 1/2 log2(1 + P_i * lambda_i^2 / sigma^2)

    where P_i is allocated by waterfilling:
    P_i = max(0, mu - sigma^2/lambda_i^2)
    with mu chosen so sum(P_i) = P_total.
    """
    n = len(singular_values)
    gains = singular_values ** 2 / sigma_noise ** 2  # lambda_i^2 / sigma^2

    # Sort gains descending
    idx = np.argsort(gains)[::-1]
    gains_sorted = gains[idx]

    # Waterfilling: find water level mu
    # P_i = max(0, mu - 1/gain_i) with sum(P_i) = P_total
    inv_gains = 1.0 / gains_sorted

    # Try using k channels (largest k gains)
    for k in range(n, 0, -1):
        mu = (P_total + np.sum(inv_gains[:k])) / k
        powers = mu - inv_gains[:k]
        if np.all(powers >= 0):
            break

    # Capacity
    C = 0.0
    channel_info = []
    for i in range(k):
        snr_i = powers[i] * gains_sorted[i]
        c_i = 0.5 * np.log2(1 + snr_i)
        C += c_i
        channel_info.append((singular_values[idx[i]], powers[i], snr_i, c_i))

    return C, channel_info


# ====================================================================
# Step 5: Codebook verification
# ====================================================================

def codebook_test(A, y_ref, gamma_ref, sigma_noise, n_symbols, n_trials=500,
                  spread=0.02):
    """Generate random codebook and test ML decoding error rate."""
    np.random.seed(42)

    # Generate codebook: n_symbols random gamma profiles
    codebook_gammas = []
    for _ in range(n_symbols):
        delta_g = np.random.randn(N) * spread
        gammas = np.array(gamma_ref) + delta_g
        gammas = np.clip(gammas, 0.005, 0.2)
        codebook_gammas.append(gammas)

    # Compute templates (noiseless outputs)
    templates = []
    for gammas in codebook_gammas:
        delta_g = np.array(gammas) - np.array(gamma_ref)
        y_approx = y_ref + A @ delta_g  # Linear approximation
        templates.append(y_approx)
    templates = np.array(templates)

    # ML decoding with noise
    errors = 0
    for trial in range(n_trials):
        idx = trial % n_symbols
        y_true = templates[idx]
        y_noisy = y_true + np.random.randn(len(y_true)) * sigma_noise
        guess = np.argmin([np.linalg.norm(y_noisy - t) for t in templates])
        if guess != idx:
            errors += 1

    return errors / n_trials


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Formal Channel Capacity of the Bidirectional Bridge")
    log("=" * 70)
    log()

    gamma_ref = [0.05] * N
    t_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

    # Step 1: Jacobian
    log("STEP 1: JACOBIAN A = dy/dg")
    log("-" * 40)
    A, y_ref = compute_jacobian(gamma_ref, t_points)
    log(f"  Feature dimension: {A.shape[0]} ({len(t_points)} times x 25 features)")
    log(f"  Input dimension: {A.shape[1]} (N={N} gamma values)")
    log(f"  Jacobian shape: {A.shape[0]} x {A.shape[1]}")
    log()

    # Step 2: SVD
    log("STEP 2: SVD OF JACOBIAN")
    log("-" * 40)
    U, sv, Vt = np.linalg.svd(A, full_matrices=False)
    log(f"  Singular values: {', '.join(f'{s:.4f}' for s in sv)}")
    log(f"  Condition number: {sv[0]/sv[-1]:.1f}")
    log(f"  Effective rank (sv > 1% of max): {np.sum(sv > 0.01*sv[0])}")
    log()

    # Explain singular values
    log("  Each singular value = gain of one independent channel:")
    for i, s in enumerate(sv):
        log(f"    Channel {i+1}: gain = {s:.4f}, "
            f"direction = [{', '.join(f'{v:.2f}' for v in Vt[i])}]")
    log()

    # Step 3: Capacity vs noise
    log("STEP 3: CHANNEL CAPACITY (Waterfilling)")
    log("-" * 40)
    log()

    # Power constraint: P = N * spread^2 (variance of gamma perturbation)
    spreads = [0.01, 0.02, 0.03, 0.04, 0.05]
    noise_levels = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]

    log(f"  {'P_gamma':>8}" + "".join(f"  {'s=' + f'{s:.3f}':>9}" for s in noise_levels))
    log("  " + "-" * (8 + 11 * len(noise_levels)))

    for spread in spreads:
        P_total = N * spread ** 2
        row = f"  {spread:8.3f}"
        for sigma in noise_levels:
            C, _ = waterfilling_capacity(sv, P_total, sigma)
            row += f"  {C:9.2f}"
        log(row)

    log()
    log("  (Values in bits per channel use)")
    log()

    # Detailed breakdown for reference case
    spread_ref = 0.02
    sigma_ref = 0.01
    P_ref = N * spread_ref ** 2
    C_ref, channels = waterfilling_capacity(sv, P_ref, sigma_ref)

    log(f"  Reference case: spread={spread_ref}, sigma={sigma_ref}")
    log(f"  Total capacity: {C_ref:.2f} bits")
    log()
    log(f"  {'Channel':>8}  {'Gain':>8}  {'Power':>8}  {'SNR':>8}  {'Bits':>8}")
    log("  " + "-" * 44)
    for i, (gain, power, snr, bits) in enumerate(channels):
        log(f"  {i+1:>8}  {gain:8.4f}  {power:8.6f}  {snr:8.1f}  {bits:8.2f}")

    log()

    # Step 4: Capacity curve
    log("STEP 4: CAPACITY vs NOISE (spread=0.02)")
    log("-" * 40)
    log()
    log(f"  {'sigma':>8}  {'C (bits)':>8}  {'2^C symbols':>11}  {'vs empirical':>12}")
    log("  " + "-" * 44)

    P_fixed = N * 0.02 ** 2
    for sigma in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
        C, _ = waterfilling_capacity(sv, P_fixed, sigma)
        n_sym = 2 ** C
        empirical = ">> 2 bits" if C > 2 else "~ 2 bits" if C > 1.5 else "< 2 bits"
        log(f"  {sigma:8.3f}  {C:8.2f}  {n_sym:11.1f}  {empirical:>12}")

    log()

    # Step 5: Codebook verification
    log("STEP 5: CODEBOOK VERIFICATION")
    log("-" * 40)
    log()

    log(f"  {'Rate (bits)':>11}  {'n_symbols':>9}  {'sigma=0.01':>10}  {'sigma=0.05':>10}  {'sigma=0.10':>10}")
    log("  " + "-" * 56)

    for rate in [1, 2, 3, 4, 5, 6, 8]:
        n_sym = 2 ** rate
        if n_sym > 256:
            continue
        err_01 = codebook_test(A, y_ref, gamma_ref, 0.01, n_sym, n_trials=500, spread=0.02)
        err_05 = codebook_test(A, y_ref, gamma_ref, 0.05, n_sym, n_trials=500, spread=0.02)
        err_10 = codebook_test(A, y_ref, gamma_ref, 0.10, n_sym, n_trials=500, spread=0.02)
        log(f"  {rate:>11}  {n_sym:>9}  {err_01:10.1%}  {err_05:10.1%}  {err_10:10.1%}")

    log()

    # Summary
    log("=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()
    log(f"  Channel: gamma_profile (R^{N}) -> observables (R^{A.shape[0]})")
    log(f"  Singular values: {', '.join(f'{s:.3f}' for s in sv)}")
    log(f"  Effective rank: {np.sum(sv > 0.01*sv[0])}/{N}")
    log(f"  Capacity at sigma=0.01, spread=0.02: {C_ref:.2f} bits ({2**C_ref:.0f} symbols)")
    log(f"  Our empirical result: 2 bits (4 symbols, 100%)")
    log(f"  Theoretical headroom: {C_ref - 2:.1f} bits above our test")
    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
