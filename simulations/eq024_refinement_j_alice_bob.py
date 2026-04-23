#!/usr/bin/env python3
"""
EQ-024 refinement: J-inversion of F30 Alice-Bob classification (Test 3)
=========================================================================
F30 (GAMMA_AS_SIGNAL Test 3) had Alice pick one of 4 gamma profiles
(Gradient ->, Gradient <-, Mountain, Valley) and Bob classify at 100%
accuracy for sigma <= 0.01 using 10-dim feature vector (5 purities +
4 nearest-neighbor CPsi + 1 end-to-end MI). That is the iconic
operational claim that made gamma a "readable channel".

Under gamma_0 = const, Alice cannot set per-site gamma. This script
asks: can Alice set 4 J-profiles that give analogous classification
accuracy? Direct operational test of the gamma_0 = const reframing.

Setup:
- N = 5, gamma_0 = 0.05 uniform (gamma_0 = const)
- Uniform J_baseline = 1.0, per-bond perturbations of magnitude delta
- 4 J-profiles: Gradient->, Gradient<-, Mountain (peak at interior bonds),
  Valley (dip at interior bonds)
- 10-feature vector at single time point
- Template matching with Gaussian feature noise, 200 trials per sigma

Receivers tested:
- |+>^5:     F30's OPTIMAL gamma-receiver, predicted J-BLIND (Class 3
             M_x-polynomial), accuracy expected ~25% (random guess)
- |+-+-+>:   SU(2)-breaking F71-symmetric, the J-capacity 11.92 bit
             champion from the morning RESULT
- |01010>:   another SU(2)-breaking F71-symmetric, J-capacity 11.53 bits

Output: simulations/results/eq024_refinement/j_alice_bob.{txt,json}
"""

import json
import os
import sys
import time as _time

import numpy as np
from scipy.linalg import expm

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "eq024_refinement")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_TXT = os.path.join(OUT_DIR, "j_alice_bob.txt")
OUT_JSON = os.path.join(OUT_DIR, "j_alice_bob.json")

_outf = open(OUT_TXT, "w", encoding="utf-8", buffering=1)
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
minus = (up - dn) / np.sqrt(2)

N = 5
d = 32
d2 = 1024
N_BONDS = N - 1


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


_BOND_BLOCKS = []
for b in range(N_BONDS):
    block = (site_op(sx, b) @ site_op(sx, b + 1)
             + site_op(sy, b) @ site_op(sy, b + 1)
             + site_op(sz, b) @ site_op(sz, b + 1))
    _BOND_BLOCKS.append(block)


def build_H(J_vec):
    H = np.zeros((d, d), dtype=complex)
    for b in range(N_BONDS):
        H += J_vec[b] * _BOND_BLOCKS[b]
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


def von_neumann(rho):
    """Von Neumann entropy in bits."""
    evals = np.linalg.eigvalsh(rho)
    evals = np.real(evals)
    evals = evals[evals > 1e-12]
    if len(evals) == 0:
        return 0.0
    return float(-np.sum(evals * np.log2(evals)))


def mutual_information(rho_full, sites_A, sites_B):
    rho_A = ptrace_keep(rho_full, sites_A)
    rho_B = ptrace_keep(rho_full, sites_B)
    rho_AB = ptrace_keep(rho_full, list(sites_A) + list(sites_B))
    return von_neumann(rho_A) + von_neumann(rho_B) - von_neumann(rho_AB)


def extract_features_F30_test3(rho):
    """10-dim feature vector analog to F30 Test 3:
    5 single-qubit purities + 4 nearest-neighbor CPsi + 1 end-to-end MI.
    """
    f = []
    for q in range(N):
        f.append(purity(ptrace_keep(rho, [q])))
    for b in range(N - 1):
        f.append(cpsi(ptrace_keep(rho, [b, b + 1])))
    f.append(mutual_information(rho, [0], [N - 1]))
    return np.array(f)


def kron_all(vs):
    out = vs[0]
    for v in vs[1:]:
        out = np.kron(out, v)
    return out


def _state_from_bits(bits):
    psi = np.zeros(d, dtype=complex)
    idx = int("".join(str(b) for b in bits), 2)
    psi[idx] = 1.0
    return psi


# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------

gamma_0 = 0.05
J_baseline = 1.0
delta = 0.2  # per-bond perturbation magnitude
t_measure = 2.0  # single time point (F30 Test 3 optimum is t=1-3)

# 4 J-profiles, per-bond values (4 bonds at N=5)
#   Gradient-> : linearly increasing left-to-right
#   Gradient<- : linearly decreasing (mirror)
#   Mountain   : peak at interior bonds
#   Valley     : dip at interior bonds
# All F71-symmetric except the two gradients (which are F71-mirror partners).
J_PROFILES = {
    "Gradient->":  [J_baseline + d_b for d_b in
                    [-1.5 * delta, -0.5 * delta, +0.5 * delta, +1.5 * delta]],
    "Gradient<-":  [J_baseline + d_b for d_b in
                    [+1.5 * delta, +0.5 * delta, -0.5 * delta, -1.5 * delta]],
    "Mountain":    [J_baseline + d_b for d_b in
                    [-delta, +delta, +delta, -delta]],
    "Valley":      [J_baseline + d_b for d_b in
                    [+delta, -delta, -delta, +delta]],
}

# Receivers
RECEIVERS = {
    "|+>^5 (F30-optimal-gamma, predicted J-BLIND)":
        kron_all([plus] * N),
    "|+-+-+> (morning 11.92-bit champion)":
        kron_all([plus, minus, plus, minus, plus]),
    "|01010> (morning alt)":
        _state_from_bits([0, 1, 0, 1, 0]),
}


# ------------------------------------------------------------------
# Classification pipeline
# ------------------------------------------------------------------

def compute_templates(rho0, t):
    """For a given receiver and time, return dict {profile_name: feature_vec}."""
    templates = {}
    for name, J_vec in J_PROFILES.items():
        H = build_H(J_vec)
        L = build_L(H, [gamma_0] * N)
        rho_t = evolve(L, rho0, t)
        templates[name] = extract_features_F30_test3(rho_t)
    return templates


def classify_trial(feature_vec_noisy, templates):
    """Return name of nearest template."""
    best_name = None
    best_dist = np.inf
    for name, tpl in templates.items():
        dist = np.linalg.norm(feature_vec_noisy - tpl)
        if dist < best_dist:
            best_dist = dist
            best_name = name
    return best_name


def classification_accuracy(templates, sigma, n_trials=200, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)
    names = list(templates.keys())
    correct = 0
    for trial in range(n_trials):
        true_name = names[trial % len(names)]
        true_vec = templates[true_name]
        noisy = true_vec + rng.normal(0, sigma, size=true_vec.shape)
        guessed = classify_trial(noisy, templates)
        if guessed == true_name:
            correct += 1
    return correct / n_trials


def template_distance_matrix(templates):
    names = list(templates.keys())
    n = len(names)
    D = np.zeros((n, n))
    for i, ni in enumerate(names):
        for j, nj in enumerate(names):
            D[i, j] = np.linalg.norm(templates[ni] - templates[nj])
    return D, names


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

if __name__ == "__main__":
    t_start = _time.time()

    log("EQ-024 refinement: J-inversion of F30 Alice-Bob classification (Test 3)")
    log("=" * 72)
    log(f"N = {N}, gamma_0 = {gamma_0} uniform (gamma_0 = const)")
    log(f"J_baseline = {J_baseline}, per-bond delta = {delta}")
    log(f"Measurement time t = {t_measure}")
    log(f"Feature vector: 10-dim (5 purities + 4 NN CPsi + 1 end-to-end MI)")
    log(f"Template matching, 200 trials per sigma, seed 42")
    log()
    log("J-profiles (per-bond values):")
    for name, J_vec in J_PROFILES.items():
        log(f"  {name:>12}  J = [{', '.join(f'{j:.2f}' for j in J_vec)}]")
    log()

    sigmas = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    all_results = {}

    for rx_name, psi in RECEIVERS.items():
        log("=" * 72)
        log(f"Receiver: {rx_name}")
        log("=" * 72)

        rho0 = np.outer(psi, psi.conj())
        rho0 = (rho0 + rho0.conj().T) / 2
        rho0 /= np.trace(rho0).real

        t0 = _time.time()
        templates = compute_templates(rho0, t_measure)
        dt = _time.time() - t0
        log(f"  Template compute time: {dt:.1f} s")
        log()

        # Template distance matrix
        D, names = template_distance_matrix(templates)
        log("  Template distance matrix (L2 in 10-feature space):")
        header = "    " + " ".join(f"{n:>11}" for n in names)
        log(header)
        for i, ni in enumerate(names):
            row = f"    {ni:>11}"
            for j in range(len(names)):
                row += f"  {D[i,j]:9.4f}"
            log(row)
        d_min = float(np.min(D[D > 1e-15])) if np.any(D > 1e-15) else 0.0
        log(f"  d_min (closest symbol pair): {d_min:.4e}")
        log()

        # Accuracy per sigma
        log(f"  {'sigma':>8}  {'accuracy':>10}  {'vs F30 Test 3':>18}")
        log("  " + "-" * 42)

        # F30 reference at |+>^5 Heisenberg, Test 3, same 10-feat @ t=1
        F30_ACC = {0.0: 1.00, 0.001: 1.00, 0.005: 1.00, 0.01: 1.00,
                   0.02: 0.9, 0.05: 0.595, 0.1: 0.31}

        rx_results = {"templates_d_min": d_min,
                      "template_distance_matrix": D.tolist(),
                      "template_names": names,
                      "accuracy_per_sigma": {}}

        for sigma in sigmas:
            rng = np.random.default_rng(42)
            acc = classification_accuracy(templates, sigma, n_trials=200, rng=rng)
            f30 = F30_ACC.get(sigma, None)
            f30_str = f"{f30*100:.1f}%" if f30 is not None else "n/a"
            log(f"  {sigma:8.3f}  {acc*100:9.1f}%  {f30_str:>18}")
            rx_results["accuracy_per_sigma"][str(sigma)] = float(acc)
        log()

        all_results[rx_name] = rx_results

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log("=" * 72)
    log("SUMMARY")
    log("=" * 72)
    log()
    log(f"  {'receiver':<55}  {'d_min':>10}  {'acc @ sig=0':>12}  {'acc @ sig=0.01':>14}")
    log("  " + "-" * 96)
    for name, res in all_results.items():
        d = res["templates_d_min"]
        a0 = res["accuracy_per_sigma"]["0.0"]
        a01 = res["accuracy_per_sigma"]["0.01"]
        log(f"  {name:<55}  {d:10.4e}  {a0*100:11.1f}%  {a01*100:13.1f}%")
    log()
    log(f"F30 reference (gamma-modulation at |+>^5, Heisenberg): 100% at sig<=0.01")
    log(f"Random-guess baseline (4 symbols): 25.0%")
    log()
    log(f"Total runtime: {_time.time() - t_start:.1f} s")

    out = {
        "meta": {
            "N": N,
            "gamma_0": gamma_0,
            "J_baseline": J_baseline,
            "delta": delta,
            "t_measure": t_measure,
            "features": "5 purities + 4 NN CPsi + 1 end-to-end MI (F30 Test 3)",
            "n_trials": 200,
            "seed": 42,
            "sigmas": sigmas,
            "J_profiles": {k: v for k, v in J_PROFILES.items()},
        },
        "results": all_results,
        "F30_reference": F30_ACC,
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fj:
        json.dump(out, fj, indent=2)
    log(f"JSON: {OUT_JSON}")
    log(f"Text: {OUT_TXT}")
    _outf.close()
