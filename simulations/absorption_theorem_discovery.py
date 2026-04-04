#!/usr/bin/env python3
"""
Does E = m*gamma^2 Hold in the Lindblad Cavity?
================================================
Tests candidate ratios from eigendecomposition for a
mass-energy-gamma relationship.

R=CΨ² Project, Homework #11
Source: TASK_MASS_ENERGY_GAMMA_SQUARED.md
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace')

import numpy as np
from pathlib import Path
import time
import itertools

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

TOL = 1e-8

# ── Pauli matrices ────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, Xm, Ym, Zm]
PAULI_LABELS = ['I', 'X', 'Y', 'Z']

out = []
def log(msg=""):
    print(msg)
    out.append(msg)


# ── Utilities ─────────────────────────────────────────────────────
def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]

def build_hamiltonian(N, bonds, J=1.0):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N; ops[a] = P; ops[b] = P
            H += J * kron_chain(ops)
    return H

def build_L_H(H, d):
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))

def build_L_D(N, gammas):
    d = 2**N
    Id = np.eye(d, dtype=complex)
    L_D = np.zeros((d**2, d**2), dtype=complex)
    for k in range(N):
        ops = [I2] * N; ops[k] = Zm
        Lk = np.sqrt(gammas[k]) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L_D += (np.kron(Lk, Lk.conj())
                - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T)))
    return L_D

def classify_paired(eigvals, sigma_gamma, tol=TOL):
    n = len(eigvals)
    paired = np.zeros(n, dtype=bool)
    partner_idx = -np.ones(n, dtype=int)
    used = set()
    for i in range(n):
        if i in used:
            continue
        target = -2 * sigma_gamma - eigvals[i].conjugate()
        dists = np.abs(eigvals - target)
        dists[i] = 999
        for j in used:
            dists[j] = 999
        j = np.argmin(dists)
        if dists[j] < tol:
            paired[i] = paired[j] = True
            partner_idx[i] = j
            partner_idx[j] = i
            used.add(i)
            used.add(j)
    return paired, partner_idx

def generate_pauli_basis(N):
    vecs, labels, n_xy = [], [], []
    for indices in itertools.product(range(4), repeat=N):
        labels.append(''.join(PAULI_LABELS[i] for i in indices))
        vecs.append(kron_chain([PAULIS[i] for i in indices]).flatten())
        n_xy.append(sum(1 for i in indices if i in (1, 2)))
    return labels, np.array(vecs).T, np.array(n_xy)


# ── Core computation ──────────────────────────────────────────────
def compute_all(N, J, gamma):
    """Eigendecomposition + sector weights for all modes."""
    d = 2**N
    dim = d**2
    bonds = chain_bonds(N)
    gammas_list = [gamma] * N
    sigma_gamma = sum(gammas_list)

    H = build_hamiltonian(N, bonds, J)
    L = build_L_H(H, d) + build_L_D(N, gammas_list)
    eigvals, eigvecs = np.linalg.eig(L)

    _, pauli_mat, n_xy = generate_pauli_basis(N)
    coeffs = pauli_mat.conj().T @ eigvecs / np.sqrt(d)
    wt = np.abs(coeffs)**2

    sector_w = np.zeros((N + 1, dim))
    for k in range(N + 1):
        sector_w[k] = np.sum(wt[n_xy == k], axis=0)
    totals = np.sum(sector_w, axis=0)
    sector_w /= np.maximum(totals[np.newaxis, :], 1e-30)

    # Mean X/Y count per mode: <n_XY> = sum_k k * sector_w[k]
    mean_nxy = np.zeros(dim)
    for k in range(N + 1):
        mean_nxy += k * sector_w[k]

    paired, pidx = classify_paired(eigvals, sigma_gamma)

    return dict(
        eigvals=eigvals,
        omega=np.abs(eigvals.imag),
        alpha=np.abs(eigvals.real),
        w_IZ=sector_w[0],
        w_XY=1.0 - sector_w[0],
        sector_w=sector_w,
        mean_nxy=mean_nxy,
        paired=paired, pidx=pidx,
        sigma_gamma=sigma_gamma, N=N, J=J, gamma=gamma
    )


def extract_pairs(p):
    """Extract oscillatory palindromic pairs."""
    ev = p['eigvals']
    dim = len(ev)
    pairs = []
    seen = set()
    for i in range(dim):
        j = int(p['pidx'][i])
        if not p['paired'][i] or j < 0 or i in seen:
            continue
        seen.add(i); seen.add(j)
        if abs(ev[i].real) < abs(ev[j].real):
            fi, si = i, j
        else:
            fi, si = j, i
        o = p['omega'][fi]
        if o < TOL:
            continue
        pairs.append(dict(
            omega=o,
            af=p['alpha'][fi], as_=p['alpha'][si],
            wiz_f=p['w_IZ'][fi], wiz_s=p['w_IZ'][si],
            wxy_f=p['w_XY'][fi], wxy_s=p['w_XY'][si],
            mnxy_f=p['mean_nxy'][fi], mnxy_s=p['mean_nxy'][si],
            sg=p['sigma_gamma'], gamma=p['gamma'], N=p['N']
        ))
    return pairs


# ── Candidate ratio definitions ───────────────────────────────────
RATIO_NAMES = {
    'A': 'omega / w_IZ(fast)',
    'B': 'omega / (w_IZ(fast) * alpha_fast)',
    'C': 'w_XY(fast) * omega / w_IZ(fast)',
    'D': 'omega / (2*Sg * w_IZ(fast))',
    'E': 'omega / |Dw_IZ|',
    'F': 'alpha_f * alpha_s / omega',
    'G': 'omega^2 / (alpha_f * alpha_s)',
    'H': 'omega / Sg',
    'I': '(alpha_s - alpha_f) * omega / Sg^2',
    'J': '|lam_f|^2 / (Sg^2 * w_IZ_f)',
    'K': 'alpha_f / (2*gamma * <n_XY>_f)',
    'L': 'alpha_f / (w_XY_f * 2*Sg)',
}


def candidate_ratios(pair):
    """Compute all candidate ratios for one pair."""
    o = pair['omega']
    af, as_ = pair['af'], pair['as_']
    wf, ws = pair['wiz_f'], pair['wiz_s']
    wxf, wxs = pair['wxy_f'], pair['wxy_s']
    sg = pair['sg']
    gamma = pair['gamma']
    nf = pair['mnxy_f']

    R = {}
    if wf > 1e-6:
        R['A'] = o / wf
    if wf > 1e-6 and af > 1e-6:
        R['B'] = o / (wf * af)
    if wf > 1e-6:
        R['C'] = wxf * o / wf
    if wf > 1e-6:
        R['D'] = o / (2 * sg * wf)
    dw = abs(wf - ws)
    if dw > 1e-6:
        R['E'] = o / dw
    R['F'] = af * as_ / o
    prod = af * as_
    if prod > 1e-12:
        R['G'] = o**2 / prod
    R['H'] = o / sg
    R['I'] = (as_ - af) * o / sg**2
    if wf > 1e-6:
        R['J'] = (af**2 + o**2) / (sg**2 * wf)
    if nf > 1e-6:
        R['K'] = af / (2 * gamma * nf)
    if wxf > 1e-6:
        R['L'] = af / (wxf * 2 * sg)
    return R


def ratio_stats(pairs):
    """Statistics for each ratio across pairs."""
    if not pairs:
        return {}
    all_R = {}
    for pair in pairs:
        for name, val in candidate_ratios(pair).items():
            all_R.setdefault(name, []).append(val)
    stats = {}
    for name, vals in all_R.items():
        arr = np.array(vals)
        arr = arr[np.isfinite(arr)]
        if len(arr) < 2:
            continue
        m = np.mean(arr)
        md = np.median(arr)
        s = np.std(arr)
        cv = s / abs(m) if abs(m) > 1e-15 else float('inf')
        stats[name] = dict(n=len(arr), mean=m, median=md, std=s, cv=cv)
    return stats


# ══════════════════════════════════════════════════════════════════
# STEPS 1-2: Mode Properties and Candidate Ratios
# ══════════════════════════════════════════════════════════════════
def step1_2():
    log("=" * 75)
    log("STEPS 1-2: Candidate Ratios Across All Oscillatory Pairs")
    log("=" * 75)
    log()

    # Print legend
    log("Ratio definitions:")
    for name in sorted(RATIO_NAMES):
        log(f"  {name:>2}: {RATIO_NAMES[name]}")
    log()

    best_per_N = {}

    for N in range(2, 6):
        p = compute_all(N, J=1.0, gamma=0.05)
        pairs = extract_pairs(p)
        stats = ratio_stats(pairs)

        log(f"-- N = {N}  ({len(pairs)} oscillatory pairs, "
            f"Sg = {p['sigma_gamma']:.4f}) " + "-" * 30)
        log(f"  {'R':>2}  {'n':>4}  {'mean':>14}  {'median':>14}  "
            f"{'std':>14}  {'CV':>10}")
        log(f"  " + "-" * 68)

        sorted_s = sorted(stats.items(), key=lambda x: x[1]['cv'])
        for name, s in sorted_s:
            mark = " ***" if s['cv'] < 0.05 else \
                   " **" if s['cv'] < 0.15 else \
                   " *" if s['cv'] < 0.30 else ""
            log(f"  {name:>2}  {s['n']:4d}  {s['mean']:14.6e}  "
                f"{s['median']:14.6e}  {s['std']:14.6e}  "
                f"{s['cv']:10.4f}{mark}")

        best_per_N[N] = sorted_s[0] if sorted_s else None
        log()

    log("Best (lowest CV) per N:")
    for N, best in best_per_N.items():
        if best:
            name, s = best
            log(f"  N={N}: Ratio {name}  CV = {s['cv']:.4f}  "
                f"median = {s['median']:.6e}")
    log()
    return best_per_N


# ══════════════════════════════════════════════════════════════════
# STEP 3: Gamma Sweep
# ══════════════════════════════════════════════════════════════════
def step3():
    log("=" * 75)
    log("STEP 3: Gamma Sweep (N=3, J=1.0)")
    log("=" * 75)
    log("Fit: R_median = a * gamma^b")
    log()

    gammas = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
    N, J = 3, 1.0

    all_medians = {name: [] for name in RATIO_NAMES}
    for gamma in gammas:
        p = compute_all(N, J, gamma)
        pairs = extract_pairs(p)
        stats = ratio_stats(pairs)
        for name in RATIO_NAMES:
            if name in stats:
                all_medians[name].append(stats[name]['median'])
            else:
                all_medians[name].append(np.nan)

    # Fit and report
    fits = {}
    log(f"  {'R':>2}  {'exponent b':>12}  {'prefactor a':>14}  "
        f"{'R^2':>10}  {'verdict':>20}")
    log(f"  " + "-" * 66)

    for name in sorted(RATIO_NAMES):
        vals = np.array(all_medians[name])
        g = np.array(gammas)
        valid = np.isfinite(vals) & (vals > 0)
        if np.sum(valid) < 3:
            continue

        log_g = np.log(g[valid])
        log_v = np.log(vals[valid])
        b, log_a = np.polyfit(log_g, log_v, 1)
        a = np.exp(log_a)

        pred = a * g[valid]**b
        ss_res = np.sum((vals[valid] - pred)**2)
        ss_tot = np.sum((vals[valid] - np.mean(vals[valid]))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 0

        if abs(b - 2) < 0.15:
            verdict = "gamma^2 ***"
        elif abs(b - 1) < 0.15:
            verdict = "gamma^1 (linear)"
        elif abs(b) < 0.15:
            verdict = "gamma-independent"
        elif abs(b + 1) < 0.15:
            verdict = "1/gamma"
        elif abs(b + 2) < 0.15:
            verdict = "1/gamma^2"
        else:
            verdict = f"gamma^{b:.2f}"

        log(f"  {name:>2}  {b:12.4f}  {a:14.6e}  {r2:10.6f}  {verdict:>20}")
        fits[name] = dict(b=b, a=a, r2=r2)

    # Detail for most interesting fits
    log()
    log("Detailed gamma dependence for key ratios:")
    for name in sorted(RATIO_NAMES):
        if name not in fits:
            continue
        b = fits[name]['b']
        if fits[name]['r2'] < 0.95:
            continue
        vals = all_medians[name]
        log(f"  Ratio {name} (b = {b:.4f}):")
        log(f"    gamma:  " + "  ".join(f"{g:8.4f}" for g in gammas))
        log(f"    median: " + "  ".join(
            f"{v:8.4e}" if np.isfinite(v) else "     NaN" for v in vals))
        log()

    return fits


# ══════════════════════════════════════════════════════════════════
# STEP 4: J Sweep
# ══════════════════════════════════════════════════════════════════
def step4():
    log("=" * 75)
    log("STEP 4: J Sweep (N=3, gamma=0.05)")
    log("=" * 75)
    log("Testing if ratio depends on coupling strength J")
    log()

    J_values = [0.1, 0.5, 1.0, 2.0, 5.0]
    N, gamma = 3, 0.05

    all_medians = {name: [] for name in RATIO_NAMES}
    for J in J_values:
        p = compute_all(N, J, gamma)
        pairs = extract_pairs(p)
        stats = ratio_stats(pairs)
        for name in RATIO_NAMES:
            if name in stats:
                all_medians[name].append(stats[name]['median'])
            else:
                all_medians[name].append(np.nan)

    log(f"  {'R':>2}  {'CV across J':>12}  {'mean':>14}  {'verdict':>18}")
    log(f"  " + "-" * 50)

    for name in sorted(RATIO_NAMES):
        vals = np.array(all_medians[name])
        valid = np.isfinite(vals)
        if np.sum(valid) < 2:
            continue
        v = vals[valid]
        cv = np.std(v) / abs(np.mean(v)) if abs(np.mean(v)) > 1e-15 \
            else float('inf')
        verdict = "J-INDEPENDENT" if cv < 0.1 else "J-DEPENDENT"
        log(f"  {name:>2}  {cv:12.4f}  {np.mean(v):14.6e}  {verdict:>18}")

    # Detail for each ratio
    log()
    log("Values across J:")
    for name in sorted(RATIO_NAMES):
        vals = all_medians[name]
        if all(np.isnan(v) for v in vals):
            continue
        log(f"  {name}: " + "  ".join(
            f"{v:.4e}" if np.isfinite(v) else "NaN" for v in vals))
    log()


# ══════════════════════════════════════════════════════════════════
# STEP 5: N Dependence
# ══════════════════════════════════════════════════════════════════
def step5():
    log("=" * 75)
    log("STEP 5: N Dependence (J=1.0, gamma=0.05)")
    log("=" * 75)
    log()

    N_values = [2, 3, 4, 5]
    J, gamma = 1.0, 0.05

    all_medians = {name: [] for name in RATIO_NAMES}
    for N in N_values:
        p = compute_all(N, J, gamma)
        pairs = extract_pairs(p)
        stats = ratio_stats(pairs)
        for name in RATIO_NAMES:
            if name in stats:
                all_medians[name].append(stats[name]['median'])
            else:
                all_medians[name].append(np.nan)

    log(f"  {'R':>2}  {'CV across N':>12}  {'mean':>14}  {'verdict':>18}")
    log(f"  " + "-" * 50)

    for name in sorted(RATIO_NAMES):
        vals = np.array(all_medians[name])
        valid = np.isfinite(vals)
        if np.sum(valid) < 2:
            continue
        v = vals[valid]
        cv = np.std(v) / abs(np.mean(v)) if abs(np.mean(v)) > 1e-15 \
            else float('inf')
        verdict = "N-INDEPENDENT" if cv < 0.1 else "N-DEPENDENT"
        log(f"  {name:>2}  {cv:12.4f}  {np.mean(v):14.6e}  {verdict:>18}")

    log()
    log("Values across N:")
    for name in sorted(RATIO_NAMES):
        vals = all_medians[name]
        if all(np.isnan(v) for v in vals):
            continue
        log(f"  {name}: " + "  ".join(
            f"{v:.4e}" if np.isfinite(v) else "NaN" for v in vals))
    log()


# ══════════════════════════════════════════════════════════════════
# STEP 6: Per-mode analysis (not pair-level)
# ══════════════════════════════════════════════════════════════════
def step6():
    log("=" * 75)
    log("STEP 6: Per-Mode Ratio  alpha / (2*gamma * <n_XY>)")
    log("=" * 75)
    log("Pure dissipator predicts this ratio = 1 for all modes.")
    log("Does it survive when the Hamiltonian is added?")
    log()

    for N in range(2, 6):
        p = compute_all(N, J=1.0, gamma=0.05)
        alpha = p['alpha']
        nxy = p['mean_nxy']
        gamma = p['gamma']
        dim = len(alpha)

        # Only modes with alpha > 0 and n_XY > 0
        active = (alpha > TOL) & (nxy > TOL)
        n_active = int(np.sum(active))
        if n_active < 2:
            log(f"  N={N}: insufficient active modes")
            continue

        ratio = alpha[active] / (2 * gamma * nxy[active])
        log(f"  N={N}: {n_active} active modes")
        log(f"    mean = {np.mean(ratio):.6f},  median = {np.median(ratio):.6f},  "
            f"std = {np.std(ratio):.6f},  CV = {np.std(ratio)/np.mean(ratio):.4f}")

        # Histogram summary
        bins = [0, 0.5, 0.9, 1.1, 1.5, 2.0, 5.0, 100.0]
        for lo, hi in zip(bins[:-1], bins[1:]):
            ct = int(np.sum((ratio >= lo) & (ratio < hi)))
            if ct > 0:
                log(f"    [{lo:.1f}, {hi:.1f}): {ct} modes")

    # Gamma sweep for this specific ratio
    log()
    log("  Gamma sweep of alpha/(2*gamma*<n_XY>) at N=3:")
    gammas = [0.01, 0.05, 0.1, 0.5, 1.0]
    for gamma in gammas:
        p = compute_all(3, J=1.0, gamma=gamma)
        alpha = p['alpha']
        nxy = p['mean_nxy']
        active = (alpha > TOL) & (nxy > TOL)
        if np.sum(active) < 2:
            continue
        ratio = alpha[active] / (2 * gamma * nxy[active])
        log(f"    gamma={gamma:.2f}: mean={np.mean(ratio):.6f}, "
            f"CV={np.std(ratio)/np.mean(ratio):.4f}")
    log()


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    log("Does E = m*gamma^2 Hold in the Lindblad Cavity?")
    log("=" * 75)
    log("R=CΨ² Project, Homework #11")
    log()

    t_total = time.time()

    best = step1_2()
    fits = step3()
    step4()
    step5()
    step6()

    # ── Summary ───────────────────────────────────────────────────
    log()
    log("=" * 75)
    log("VERDICT")
    log("=" * 75)
    log()

    # Find best gamma^2 candidate
    gamma2_candidates = {n: f for n, f in fits.items()
                         if abs(f['b'] - 2) < 0.2 and f['r2'] > 0.95}
    gamma1_candidates = {n: f for n, f in fits.items()
                         if abs(f['b'] - 1) < 0.2 and f['r2'] > 0.95}
    gamma0_candidates = {n: f for n, f in fits.items()
                         if abs(f['b']) < 0.2 and f['r2'] > 0.5}

    if gamma2_candidates:
        log("gamma^2 candidates found:")
        for n, f in gamma2_candidates.items():
            log(f"  Ratio {n}: b = {f['b']:.4f}, R^2 = {f['r2']:.6f}")
            log(f"    {RATIO_NAMES[n]}")
    else:
        log("No clean gamma^2 relationship found.")

    log()
    if gamma1_candidates:
        log("gamma^1 (linear) candidates found:")
        for n, f in gamma1_candidates.items():
            log(f"  Ratio {n}: b = {f['b']:.4f}, R^2 = {f['r2']:.6f}")
            log(f"    {RATIO_NAMES[n]}")
    else:
        log("No clean gamma^1 relationship found.")

    log()
    if gamma0_candidates:
        log("gamma-independent candidates:")
        for n, f in gamma0_candidates.items():
            log(f"  Ratio {n}: b = {f['b']:.4f}, R^2 = {f['r2']:.6f}")
            log(f"    {RATIO_NAMES[n]}")

    dt = time.time() - t_total
    log()
    log(f"Total time: {dt:.1f}s")

    out_path = RESULTS_DIR / "absorption_theorem_discovery.txt"
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(out))
    log(f"\n>>> Results saved to: {out_path}")
