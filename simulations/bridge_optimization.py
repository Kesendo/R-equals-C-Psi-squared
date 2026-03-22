#!/usr/bin/env python3
"""
Bridge Optimization: Make the γ-channel wider
================================================
Opt 4: Better features (all pairs + Pauli Z + decoder)
Opt 7: γ contrast sweep ([0.04,0.06] to [0.01,0.10])
Opt 5: Time series (60 features from 6 time points)
Opt 1: Initial state (Bell pairs, GHZ, Néel)

Script:  simulations/bridge_optimization.py
Output:  simulations/results/bridge_optimization.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "bridge_optimization.txt")
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

# Alice's alphabet (parameterized by contrast)
def make_alphabet(center=0.05, spread=0.02):
    lo, hi = center - spread, center + spread
    m = center
    return {
        'Grad->':   np.linspace(lo, hi, N).tolist(),
        'Grad<-':   np.linspace(hi, lo, N).tolist(),
        'Mountain': [lo, m, hi, m, lo],
        'Valley':   [hi, m, lo, m, hi],
    }


# Feature extractors
def features_baseline(rho):
    """10 features: 5 purities + 4 nn-CΨ + MI(0,4)"""
    f = []
    for q in range(N):
        f.append(purity(ptrace_keep(rho, [q])))
    for i in range(N - 1):
        f.append(cpsi(ptrace_keep(rho, [i, i + 1])))
    # MI(0,4) approximation via purities
    rho04 = ptrace_keep(rho, [0, 4])
    f.append(purity(rho04))
    return np.array(f)


def features_extended(rho):
    """25 features: 5 purities + 10 all-pair CΨ + 5 ⟨Z_i⟩ + 5 ⟨X_i⟩"""
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


def features_timeseries(rho0, gammas, t_points):
    """Features from multiple time points."""
    L = build_L(H_chain, gammas)
    f = []
    for t in t_points:
        rho = evolve(L, rho0, t)
        f.extend(features_baseline(rho).tolist())
    return np.array(f)


def classify_accuracy(templates, profiles, rho0, t_meas, noise_std, n_trials,
                      feat_fn=None):
    """Classify profiles using nearest-template in feature space."""
    M = len(profiles)
    correct = 0
    for trial in range(n_trials):
        idx = trial % M
        gammas = profiles[idx]
        if feat_fn is None:
            L = build_L(H_chain, gammas)
            rho = evolve(L, rho0, t_meas)
            feat = features_baseline(rho)
        else:
            feat = feat_fn(rho0, gammas)
        if noise_std > 0:
            feat = feat + np.random.randn(len(feat)) * noise_std
        dists = [np.linalg.norm(feat - t) for t in templates]
        if np.argmin(dists) == idx:
            correct += 1
    return correct / n_trials


def compute_distances(templates):
    """Min and max pairwise template distance."""
    M = len(templates)
    dists = []
    for i in range(M):
        for j in range(i + 1, M):
            dists.append(np.linalg.norm(templates[i] - templates[j]))
    return min(dists), max(dists)


def run_test(name, rho0, alphabet, feat_fn, t_meas=2.0):
    """Run classification test and return results."""
    profiles = list(alphabet.values())
    M = len(profiles)

    # Build templates
    templates = []
    for gammas in profiles:
        L = build_L(H_chain, gammas)
        rho = evolve(L, rho0, t_meas)
        templates.append(feat_fn(rho))
    templates = np.array(templates)

    d_min, d_max = compute_distances(templates)
    sigma_thresh = d_min / 3

    # Accuracy at key noise levels
    n_trials = 200
    acc_0 = classify_accuracy(templates, profiles, rho0, t_meas, 0.0, n_trials,
                               lambda r0, g: feat_fn(evolve(build_L(H_chain, g), r0, t_meas)))
    acc_01 = classify_accuracy(templates, profiles, rho0, t_meas, 0.01, n_trials,
                                lambda r0, g: feat_fn(evolve(build_L(H_chain, g), r0, t_meas)))
    acc_05 = classify_accuracy(templates, profiles, rho0, t_meas, 0.05, n_trials,
                                lambda r0, g: feat_fn(evolve(build_L(H_chain, g), r0, t_meas)))

    return d_min, d_max, sigma_thresh, acc_0, acc_01, acc_05


# Initial states
def make_plus_N():
    psi = plus
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())


def make_neel():
    psi = up
    for i in range(1, N):
        psi = np.kron(psi, dn if i % 2 == 1 else up)
    return np.outer(psi, psi.conj())


def make_ghz():
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[-1] = 1 / np.sqrt(2)
    return np.outer(psi, psi.conj())


# ====================================================================
# Main optimizations
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Bridge Optimization: Make the γ-channel wider")
    log("=" * 70)
    log()

    rho0_plus = make_plus_N()
    alphabet_base = make_alphabet(0.05, 0.02)
    profiles_base = list(alphabet_base.values())

    header = f"  {'Config':>30}  {'d_min':>8}  {'d_max':>8}  {'σ_thr':>8}  {'σ=0':>6}  {'σ=.01':>6}  {'σ=.05':>6}"
    sep = "  " + "-" * 80

    # ---- Opt 4: Feature comparison ----
    log("=" * 70)
    log("OPT 4: BETTER FEATURES")
    log("=" * 70)
    log()
    log(header)
    log(sep)

    for feat_name, feat_fn in [("Baseline (10 feat)", features_baseline),
                                ("Extended (25 feat)", features_extended)]:
        d_min, d_max, s_thr, a0, a1, a5 = run_test(feat_name, rho0_plus, alphabet_base, feat_fn)
        log(f"  {feat_name:>30}  {d_min:8.4f}  {d_max:8.4f}  {s_thr:8.4f}  "
            f"{a0:5.0%}  {a1:5.0%}  {a5:5.0%}")

    log()

    # ---- Opt 7: γ contrast ----
    log("=" * 70)
    log("OPT 7: γ CONTRAST")
    log("=" * 70)
    log()
    log(header)
    log(sep)

    for spread in [0.01, 0.02, 0.03, 0.04, 0.045]:
        label = f"γ∈[{0.05-spread:.3f},{0.05+spread:.3f}]"
        alph = make_alphabet(0.05, spread)
        d_min, d_max, s_thr, a0, a1, a5 = run_test(label, rho0_plus, alph, features_baseline)
        log(f"  {label:>30}  {d_min:8.4f}  {d_max:8.4f}  {s_thr:8.4f}  "
            f"{a0:5.0%}  {a1:5.0%}  {a5:5.0%}")

    log()

    # ---- Opt 5: Time series ----
    log("=" * 70)
    log("OPT 5: TIME SERIES")
    log("=" * 70)
    log()

    t_series = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    # Build templates using time series
    ts_templates = []
    for gammas in profiles_base:
        feat = features_timeseries(rho0_plus, gammas, t_series)
        ts_templates.append(feat)
    ts_templates = np.array(ts_templates)
    d_min_ts, d_max_ts = compute_distances(ts_templates)

    # Classification with time series
    n_trials = 200
    M = len(profiles_base)
    for noise_std, label in [(0.0, "σ=0"), (0.01, "σ=0.01"), (0.05, "σ=0.05")]:
        correct = 0
        for trial in range(n_trials):
            idx = trial % M
            feat = features_timeseries(rho0_plus, profiles_base[idx], t_series)
            if noise_std > 0:
                feat += np.random.randn(len(feat)) * noise_std
            dists = [np.linalg.norm(feat - t) for t in ts_templates]
            if np.argmin(dists) == idx:
                correct += 1
        acc = correct / n_trials
        if label == "σ=0":
            log(f"  Time series ({len(t_series)}×10 = {len(t_series)*10} features)")
            log(f"  d_min={d_min_ts:.4f}, d_max={d_max_ts:.4f}, σ_thresh={d_min_ts/3:.4f}")
        log(f"  {label}: {acc:.0%}")

    log()

    # ---- Opt 1: Initial states ----
    log("=" * 70)
    log("OPT 1: INITIAL STATES")
    log("=" * 70)
    log()
    log(header)
    log(sep)

    states = [
        ("|+>^5", make_plus_N()),
        ("|01010>", make_neel()),
        ("GHZ", make_ghz()),
    ]

    for sname, rho0 in states:
        d_min, d_max, s_thr, a0, a1, a5 = run_test(sname, rho0, alphabet_base, features_baseline)
        log(f"  {sname:>30}  {d_min:8.4f}  {d_max:8.4f}  {s_thr:8.4f}  "
            f"{a0:5.0%}  {a1:5.0%}  {a5:5.0%}")

    log()

    # ---- BEST COMBINATION ----
    log("=" * 70)
    log("BEST COMBINATION: Extended features + high contrast + time series")
    log("=" * 70)
    log()

    best_spread = 0.04
    best_alph = make_alphabet(0.05, best_spread)
    best_profiles = list(best_alph.values())

    # Extended features + time series
    best_templates = []
    for gammas in best_profiles:
        f_all = []
        for t in t_series:
            L = build_L(H_chain, gammas)
            rho = evolve(L, rho0_plus, t)
            f_all.extend(features_extended(rho).tolist())
        best_templates.append(np.array(f_all))
    best_templates = np.array(best_templates)

    d_min_best, d_max_best = compute_distances(best_templates)
    log(f"  Features: {len(best_templates[0])} (extended × {len(t_series)} times)")
    log(f"  γ contrast: [{0.05-best_spread:.2f}, {0.05+best_spread:.2f}]")
    log(f"  d_min = {d_min_best:.4f} (baseline: 0.024)")
    log(f"  d_max = {d_max_best:.4f} (baseline: 0.112)")
    log(f"  σ_thresh = {d_min_best/3:.4f} (baseline: 0.008)")
    log(f"  Improvement: {d_min_best/0.024:.1f}× min distance, {d_min_best/3/0.008:.1f}× noise tolerance")
    log()

    for noise_std, label in [(0.0, "σ=0"), (0.01, "σ=0.01"), (0.05, "σ=0.05"), (0.1, "σ=0.10")]:
        correct = 0
        for trial in range(n_trials):
            idx = trial % M
            f_all = []
            for t in t_series:
                L = build_L(H_chain, best_profiles[idx])
                rho = evolve(L, rho0_plus, t)
                f_all.extend(features_extended(rho).tolist())
            feat = np.array(f_all)
            if noise_std > 0:
                feat += np.random.randn(len(feat)) * noise_std
            dists = [np.linalg.norm(feat - t) for t in best_templates]
            if np.argmin(dists) == idx:
                correct += 1
        log(f"  {label}: {correct/n_trials:.0%}")

    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
