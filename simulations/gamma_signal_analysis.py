#!/usr/bin/env python3
"""
Gamma as Signal: Can Bob read what Alice encodes in the γ profile?
====================================================================
Test 1: Random vs correlated γ profiles — can the decoder distinguish?
Test 2: Time-varying γ — detection from internal observables
Test 3: Alice-Bob channel — classification accuracy and channel capacity

Script:  simulations/gamma_signal_analysis.py
Output:  simulations/results/gamma_signal_analysis.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "gamma_signal_analysis.txt")
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
            opi = site_op(P, i)
            opj = site_op(P, i + 1)
            H += J * opi @ opj
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


def von_neumann(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-15]
    return -float(np.sum(ev * np.log2(ev)))


def mutual_info(rho, kA, kB):
    rhoA = ptrace_keep(rho, [kA])
    rhoB = ptrace_keep(rho, [kB])
    rhoAB = ptrace_keep(rho, [kA, kB])
    return von_neumann(rhoA) + von_neumann(rhoB) - von_neumann(rhoAB)


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d_r = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d_r - 1) if d_r > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


H_chain = build_H(J=1.0)

# Initial state: |+>^5 (maximum single-qubit coherence)
psi_init = plus
for _ in range(N - 1):
    psi_init = np.kron(psi_init, plus)
rho0 = np.outer(psi_init, psi_init.conj())
rho0_vec = rho0.flatten()


def extract_features(gammas, t_meas=5.0):
    """Extract Bob's feature vector from a given gamma profile at time t."""
    L = build_L(H_chain, gammas)
    rho = evolve(L, rho0, t_meas)

    features = []
    # Single-qubit purities
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        features.append(purity(rho_q))

    # Nearest-neighbor pair CΨ
    for i in range(N - 1):
        rho_pair = ptrace_keep(rho, [i, i + 1])
        features.append(cpsi(rho_pair))

    # End-to-end MI
    features.append(mutual_info(rho, 0, N - 1))

    return np.array(features)


# ====================================================================
# Test 1: Random vs Correlated γ profiles
# ====================================================================

def test_1():
    log("=" * 70)
    log("TEST 1: RANDOM vs CORRELATED γ PROFILES")
    log("=" * 70)
    log()

    np.random.seed(42)
    t_meas = 5.0

    # Scenario A: uncorrelated random
    features_A = []
    for _ in range(50):
        gammas = 0.05 + np.random.uniform(-0.02, 0.02, size=N)
        features_A.append(extract_features(gammas, t_meas))
    features_A = np.array(features_A)

    # Scenario B: structured/correlated profiles
    correlated_profiles = [
        [0.03, 0.04, 0.05, 0.06, 0.07],  # Gradient
        [0.07, 0.06, 0.05, 0.04, 0.03],  # Reverse gradient
        [0.07, 0.06, 0.05, 0.06, 0.07],  # V-shape
        [0.03, 0.04, 0.07, 0.04, 0.03],  # Peak
        [0.03, 0.03, 0.03, 0.07, 0.07],  # Step
        [0.07, 0.07, 0.03, 0.03, 0.03],  # Reverse step
        [0.03, 0.07, 0.03, 0.07, 0.03],  # Alternating
        [0.07, 0.03, 0.07, 0.03, 0.07],  # Reverse alternating
    ]

    features_B = []
    for gammas in correlated_profiles:
        features_B.append(extract_features(gammas, t_meas))
    features_B = np.array(features_B)

    # Compute spread (variance) of features within each scenario
    var_A = np.mean(np.var(features_A, axis=0))
    var_B = np.mean(np.var(features_B, axis=0))

    # Compute distance between A centroid and each B profile
    centroid_A = np.mean(features_A, axis=0)
    dists_B = [np.linalg.norm(f - centroid_A) for f in features_B]

    # How many B profiles fall outside the 95th percentile of A distances?
    dists_A = [np.linalg.norm(f - centroid_A) for f in features_A]
    threshold_95 = np.percentile(dists_A, 95)

    n_distinguishable = sum(1 for d in dists_B if d > threshold_95)

    log(f"  Scenario A: 50 random profiles, var(features) = {var_A:.6f}")
    log(f"  Scenario B: {len(correlated_profiles)} structured profiles, var(features) = {var_B:.6f}")
    log(f"  95th percentile distance in A: {threshold_95:.6f}")
    log()
    log(f"  {'Profile':>20}  {'Distance':>10}  {'Distinguishable':>15}")
    log("  " + "-" * 50)
    for i, (gammas, dist) in enumerate(zip(correlated_profiles, dists_B)):
        label = f"[{','.join(f'{g:.2f}' for g in gammas)}]"
        tag = "YES" if dist > threshold_95 else "no"
        log(f"  {label:>40}  {dist:10.6f}  {tag:>15}")

    log()
    log(f"  Distinguishable from random: {n_distinguishable}/{len(correlated_profiles)}")
    log()


# ====================================================================
# Test 2: Time-varying γ detection
# ====================================================================

def test_2():
    log("=" * 70)
    log("TEST 2: TIME-VARYING γ — Detection from internal observables")
    log("=" * 70)
    log()

    dt = 0.1
    t_max = 20.0
    n_steps = int(round(t_max / dt))

    # Scenario A: constant γ
    gammas_const = [0.05] * N

    # Scenario B: γ_2 jumps at t=10
    gammas_before = [0.05] * N
    gammas_after = [0.05, 0.05, 0.10, 0.05, 0.05]

    # Scenario C: γ_2 drifts slowly
    def gammas_drift(t):
        g2 = 0.05 + 0.03 * np.sin(2 * np.pi * t / 15)
        return [0.05, 0.05, max(0.01, g2), 0.05, 0.05]

    scenarios = [
        ("Constant", lambda t: gammas_const),
        ("Jump at t=10", lambda t: gammas_after if t >= 10 else gammas_before),
        ("Slow drift", gammas_drift),
    ]

    log(f"  Tracking MI(0,4) and single-qubit purities over time")
    log()

    for name, gamma_fn in scenarios:
        log(f"  --- {name} ---")

        # Use RK4 for time-dependent gamma
        rho = rho0.copy()
        times = []
        mi_04 = []
        pur_2 = []

        for step in range(n_steps + 1):
            t = step * dt
            times.append(t)

            mi = mutual_info(rho, 0, 4)
            mi_04.append(mi)
            rho_q2 = ptrace_keep(rho, [2])
            pur_2.append(purity(rho_q2))

            if step < n_steps:
                gammas = gamma_fn(t)
                L_ops = []
                LdL_ops = []
                for k in range(N):
                    Lk = np.sqrt(gammas[k]) * site_op(sz, k)
                    L_ops.append(Lk)
                    LdL_ops.append(Lk.conj().T @ Lk)

                # RK4 step
                def rhs(r):
                    dr = -1j * (H_chain @ r - r @ H_chain)
                    for Lk, LdLk in zip(L_ops, LdL_ops):
                        dr += Lk @ r @ Lk.conj().T - 0.5 * (LdLk @ r + r @ LdLk)
                    return dr

                k1 = rhs(rho)
                k2 = rhs(rho + 0.5 * dt * k1)
                k3 = rhs(rho + 0.5 * dt * k2)
                k4 = rhs(rho + dt * k3)
                rho = rho + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
                rho = (rho + rho.conj().T) / 2
                rho /= np.trace(rho).real

        # Print snapshots
        log(f"  {'t':>6}  {'MI(0,4)':>8}  {'Pur(Q2)':>8}")
        log("  " + "-" * 26)
        for i in range(0, len(times), n_steps // 10):
            log(f"  {times[i]:6.1f}  {mi_04[i]:8.4f}  {pur_2[i]:8.4f}")

        # Detect jump: max derivative of MI
        dmi = np.diff(mi_04)
        max_jump_idx = np.argmax(np.abs(dmi))
        max_jump_t = times[max_jump_idx]
        max_jump_val = dmi[max_jump_idx]
        log(f"  Max |dMI/dt| at t={max_jump_t:.1f}, value={max_jump_val:.6f}")
        log()


# ====================================================================
# Test 3: Alice-Bob channel
# ====================================================================

def test_3():
    log("=" * 70)
    log("TEST 3: ALICE-BOB CHANNEL — Classification + Capacity")
    log("=" * 70)
    log()

    # Alice's alphabet: 4 γ profiles
    alphabet = {
        'Gradient →': [0.03, 0.04, 0.05, 0.06, 0.07],
        'Gradient ←': [0.07, 0.06, 0.05, 0.04, 0.03],
        'Mountain':    [0.03, 0.05, 0.07, 0.05, 0.03],
        'Valley':      [0.07, 0.05, 0.03, 0.05, 0.07],
    }
    M = len(alphabet)
    names = list(alphabet.keys())
    profiles = list(alphabet.values())

    log(f"  Alice's alphabet ({M} symbols):")
    for name, gammas in alphabet.items():
        log(f"    {name}: {gammas}")
    log()

    # Compute templates at multiple times
    t_values = [1.0, 2.0, 3.0, 5.0, 8.0, 12.0]

    log("  --- Classification accuracy vs measurement time ---")
    log()

    for t_meas in t_values:
        # Build templates
        templates = []
        for gammas in profiles:
            feat = extract_features(gammas, t_meas)
            templates.append(feat)
        templates = np.array(templates)

        # Test with noise: Bob's measurement has Gaussian noise
        n_trials = 200
        noise_levels = [0.0, 0.001, 0.005, 0.01, 0.05]

        header = f"  t={t_meas:.0f}" + "".join(f"  σ={s:.3f}" for s in noise_levels)
        if t_meas == t_values[0]:
            log(f"  {'t_meas':>6}" + "".join(f"  {'σ=' + f'{s:.3f}':>8}" for s in noise_levels))
            log("  " + "-" * (6 + 10 * len(noise_levels)))

        row = f"  {t_meas:6.1f}"
        for noise_std in noise_levels:
            correct = 0
            for trial in range(n_trials):
                # Alice picks uniformly
                idx = trial % M
                gammas = profiles[idx]

                # Bob measures with noise
                feat = extract_features(gammas, t_meas)
                if noise_std > 0:
                    feat += np.random.randn(len(feat)) * noise_std

                # Bob classifies: nearest template
                dists = [np.linalg.norm(feat - t) for t in templates]
                guess = np.argmin(dists)
                if guess == idx:
                    correct += 1

            accuracy = correct / n_trials
            row += f"  {accuracy:8.1%}"

        log(row)

    log()

    # Mutual information: I(X; Y)
    log("  --- Channel capacity: I(Alice; Bob) ---")
    log()

    t_meas = 5.0
    templates = np.array([extract_features(g, t_meas) for g in profiles])

    # Pairwise distances between templates
    log(f"  Template distances at t={t_meas}:")
    log(f"  {'':>12}" + "".join(f"  {n:>12}" for n in names))
    for i, ni in enumerate(names):
        row = f"  {ni:>12}"
        for j, nj in enumerate(names):
            dist = np.linalg.norm(templates[i] - templates[j])
            row += f"  {dist:12.6f}"
        log(row)

    log()

    # Feature-by-feature comparison
    feat_names = [f"Pur(Q{q})" for q in range(N)] + \
                 [f"CΨ({i},{i+1})" for i in range(N-1)] + ["MI(0,4)"]
    log(f"  Feature breakdown at t={t_meas}:")
    log(f"  {'Feature':>12}" + "".join(f"  {n:>12}" for n in names))
    log("  " + "-" * (12 + 14 * M))
    for f_idx, f_name in enumerate(feat_names):
        row = f"  {f_name:>12}"
        for p_idx in range(M):
            row += f"  {templates[p_idx][f_idx]:12.6f}"
        log(row)

    log()

    # Effective channel capacity estimate
    # If all 4 templates are perfectly distinguishable: C = log2(4) = 2 bits
    min_dist = float('inf')
    for i in range(M):
        for j in range(i + 1, M):
            d = np.linalg.norm(templates[i] - templates[j])
            if d < min_dist:
                min_dist = d
                closest = (names[i], names[j])

    log(f"  Minimum template distance: {min_dist:.6f} ({closest[0]} vs {closest[1]})")
    log(f"  At σ=0 (perfect measurements): capacity = log2({M}) = {np.log2(M):.1f} bits")
    log(f"  At σ>{min_dist/3:.4f}: classification degrades (noise > template separation/3)")
    log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Gamma as Signal: Can Bob read Alice's γ-profile?")
    log("=" * 70)
    log(f"N={N} chain, J=1.0, initial state: |+>^{N}")
    log()

    test_1()
    test_2()
    test_3()

    log("=" * 70)
    log("VERDICT")
    log("=" * 70)
    log()
    log("If Bob can classify Alice's profile with >90% accuracy:")
    log("  → The bridge is bidirectional. γ carries readable information.")
    log()
    log("If classification fails even with perfect measurements:")
    log("  → γ acts on the system but carries no readable structure.")
    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
