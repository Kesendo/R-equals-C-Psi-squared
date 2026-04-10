"""
Optimal state search for the N=5 sacrifice-zone chain.
=======================================================
Question: does the "optimal protection state" story from
experiments/ERROR_CORRECTION_PALINDROME.md survive when we move
from N=3 uniform-gamma to N=5 sacrifice-gamma? And does any state
in the low-excitation subspace beat W5 full on slow-mode weight?

Reuses infrastructure from ibm_april_predictions.py.

Slow-mode definition (from error_correction_palindrome.py):
  slow_mask = (rates > 1e-6) & (rates < Sum(gamma) - 1e-6)

Search space: low-excitation subspace for N=5 --
  {|00000>, |10000>, |01000>, |00100>, |00010>, |00001>}
  6 computational basis states, Hamming weight <= 1.
  This is where the N=3 optimum lives (|000>, |100>, |010>, |001>).

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 9, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import numpy as np
from scipy import linalg
from scipy.optimize import minimize, differential_evolution

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import (
    heisenberg_H,
    build_liouvillian,
    bell_plus_center,
    w_full,
    w_center3_plus,
    partial_trace_to_pair,
    wootters_concurrence,
)


# =============================================================================
# Setup: same sacrifice-chain config as ibm_april_predictions.py
# =============================================================================
N = 5
J = 1.0
T2_us = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
gamma_phys = 1.0 / (2.0 * T2_us)
gamma_min = 0.05
gamma_sacrifice = gamma_phys / gamma_phys.min() * gamma_min
gamma_uniform = np.ones(N) * np.mean(gamma_sacrifice)
Sg_sac = float(np.sum(gamma_sacrifice))
Sg_uni = float(np.sum(gamma_uniform))

# Build the Liouvillian once per profile and cache the eigendecomposition.
# We need R, R_inv, and eigvals for the projection analysis.
PROFILES = {}
H = heisenberg_H(N, J)
for label, gammas in [("sacrifice", gamma_sacrifice), ("uniform", gamma_uniform)]:
    L = build_liouvillian(H, gammas)
    eigvals, R = linalg.eig(L)
    order = np.argsort(-eigvals.real)
    eigvals = eigvals[order]
    R = R[:, order]
    try:
        R_inv = linalg.inv(R)
    except linalg.LinAlgError:
        R_inv = linalg.pinv(R)
    rates = -eigvals.real  # positive dephasing rates
    Sg = float(np.sum(gammas))
    slow_mask = (rates > 1e-6) & (rates < Sg - 1e-6)
    osc_mask = np.abs(eigvals.imag) > 1e-6
    xor_mask = rates > 2 * Sg - 1e-6
    PROFILES[label] = dict(
        gammas=gammas, Sg=Sg, eigvals=eigvals, R=R, R_inv=R_inv,
        rates=rates, slow_mask=slow_mask, osc_mask=osc_mask, xor_mask=xor_mask,
    )


# =============================================================================
# Weight analysis: project a pure state into a profile's eigenbasis
# =============================================================================
def analyze_state(psi, profile):
    """Return dict with slow_wt, osc_wt, xor_wt, total_overlap for pure state psi."""
    rho0 = np.outer(psi, psi.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    probs = np.abs(c0) ** 2
    total = probs.sum()
    if total < 1e-15:
        return dict(slow_wt=0.0, osc_wt=0.0, xor_wt=0.0, total=0.0)
    return dict(
        slow_wt=float(probs[profile['slow_mask']].sum() / total * 100),
        osc_wt=float(probs[profile['osc_mask']].sum() / total * 100),
        xor_wt=float(probs[profile['xor_mask']].sum() / total * 100),
        total=float(total),
    )


def max_adjacent_concurrence(rho, N):
    """Maximum Wootters concurrence over all adjacent qubit pairs."""
    best = 0.0
    for i in range(N - 1):
        rho_ij = partial_trace_to_pair(rho, N, i, i + 1)
        c = wootters_concurrence(rho_ij)
        if c > best:
            best = c
    return float(best)


def mean_adjacent_concurrence(rho, N):
    """Mean Wootters concurrence over all adjacent pairs."""
    vals = []
    for i in range(N - 1):
        rho_ij = partial_trace_to_pair(rho, N, i, i + 1)
        vals.append(wootters_concurrence(rho_ij))
    return float(np.mean(vals))


# =============================================================================
# Low-excitation subspace and parametrisation
# =============================================================================
def low_excitation_basis(N):
    """Computational basis with Hamming weight <= 1.

    Returns a list of d-dim complex unit vectors:
      [|0...0>, |1000..>, |0100..>, ..., |..0001>]
    Length = N + 1.
    """
    d = 2 ** N
    basis = []
    v = np.zeros(d, dtype=complex)
    v[0] = 1.0
    basis.append(v)
    for k in range(N):
        v = np.zeros(d, dtype=complex)
        v[1 << (N - 1 - k)] = 1.0
        basis.append(v)
    return basis


BASIS = low_excitation_basis(N)
M_SUB = len(BASIS)  # 6 for N=5
# Parameters: M real parts + (M-1) imaginary parts (first imag fixed to 0)
N_PARAMS = 2 * M_SUB - 1


def state_from_params(params):
    """Map a real parameter vector to a normalized complex state in the subspace."""
    re = params[:M_SUB]
    im = np.zeros(M_SUB)
    im[1:] = params[M_SUB:]
    coeffs = re + 1j * im
    psi = np.zeros(len(BASIS[0]), dtype=complex)
    for k, v in enumerate(BASIS):
        psi = psi + coeffs[k] * v
    nrm = np.linalg.norm(psi)
    if nrm < 1e-15:
        return np.zeros_like(psi)
    return psi / nrm


# =============================================================================
# Objective: maximize slow-mode weight, penalise low concurrence / zero osc
# =============================================================================
# We optimise against the SACRIFICE profile, since that is the whole point.
SACRIFICE_PROFILE = PROFILES["sacrifice"]


def objective(params):
    psi = state_from_params(params)
    if np.linalg.norm(psi) < 1e-12:
        return 1e6
    weights = analyze_state(psi, SACRIFICE_PROFILE)
    slow = weights['slow_wt'] / 100.0
    osc = weights['osc_wt'] / 100.0
    rho = np.outer(psi, psi.conj())
    conc = max_adjacent_concurrence(rho, N)

    penalty = 0.0
    if conc < 0.01:
        penalty += 10.0 * (0.01 - conc)
    if osc < 0.01:
        penalty += 10.0 * (0.01 - osc)
    return -slow + penalty


def run_optimizer(n_trials=30, seed=42, verbose=True):
    rng = np.random.default_rng(seed)
    best_obj = np.inf
    best_psi = None
    best_params = None
    for trial in range(n_trials):
        x0 = rng.standard_normal(N_PARAMS) * 0.5
        res = minimize(
            objective, x0, method='Nelder-Mead',
            options={'maxiter': 5000, 'xatol': 1e-9, 'fatol': 1e-11},
        )
        if res.fun < best_obj:
            best_obj = res.fun
            best_psi = state_from_params(res.x)
            best_params = res.x
            if verbose:
                print(f"  trial {trial:3d}: new best obj = {res.fun:.6f}"
                      f"  (slow_wt = {-res.fun * 100:.2f}%)")
    return best_psi, best_obj, best_params


# =============================================================================
# Weight <= 2 subspace (16 basis states for N=5)
# =============================================================================
def weight_le2_basis(N):
    """Computational basis with Hamming weight <= 2."""
    d = 2 ** N
    basis = []
    for i in range(d):
        if bin(i).count('1') <= 2:
            v = np.zeros(d, dtype=complex)
            v[i] = 1.0
            basis.append(v)
    return basis


# =============================================================================
# Generic optimizer: works with any basis
# =============================================================================
def state_from_params_generic(params, basis):
    """Map a real parameter vector to a normalized complex state in a given basis."""
    m = len(basis)
    re = params[:m]
    im = np.zeros(m)
    im[1:] = params[m:]
    coeffs = re + 1j * im
    psi = np.zeros(len(basis[0]), dtype=complex)
    for k, v in enumerate(basis):
        psi = psi + coeffs[k] * v
    nrm = np.linalg.norm(psi)
    if nrm < 1e-15:
        return np.zeros_like(psi)
    return psi / nrm


def run_de_optimizer(basis, conc_floor=0.01, profile=None, seed=42,
                     maxiter=1000, popsize=20, verbose=True):
    """Differential evolution optimizer on any subspace basis."""
    if profile is None:
        profile = SACRIFICE_PROFILE
    m = len(basis)
    n_params = 2 * m - 1

    def obj(params):
        psi = state_from_params_generic(params, basis)
        if np.linalg.norm(psi) < 1e-12:
            return 1e6
        weights = analyze_state(psi, profile)
        slow = weights['slow_wt'] / 100.0
        rho = np.outer(psi, psi.conj())
        conc = max_adjacent_concurrence(rho, N)
        penalty = 0.0
        if conc < conc_floor:
            penalty += 10.0 * (conc_floor - conc)
        return -slow + penalty

    bounds = [(-2, 2)] * n_params
    result = differential_evolution(
        obj, bounds, seed=seed, maxiter=maxiter,
        tol=1e-10, polish=True, disp=verbose, popsize=popsize,
    )
    best_psi = state_from_params_generic(result.x, basis)
    return best_psi, result.fun, result.x, result


def pareto_frontier(basis, floors=(0.01, 0.05, 0.1, 0.2, 0.3), seed=42,
                    maxiter=500, popsize=15, verbose=False):
    """Sweep concurrence floors and return slow_wt vs concurrence tradeoff."""
    results = []
    for floor in floors:
        psi, obj, params, de_result = run_de_optimizer(
            basis, conc_floor=floor, seed=seed, verbose=False,
            maxiter=maxiter, popsize=popsize)
        w = analyze_state(psi, SACRIFICE_PROFILE)
        rho = np.outer(psi, psi.conj())
        conc = max_adjacent_concurrence(rho, N)
        results.append({
            'floor': floor, 'slow_wt': w['slow_wt'], 'osc_wt': w['osc_wt'],
            'xor_wt': w['xor_wt'], 'conc': conc, 'psi': psi, 'obj': obj,
            'nfev': de_result.nfev,
        })
        if verbose:
            print(f"  floor={floor:.2f}: slow_wt={w['slow_wt']:.2f}%  C_max={conc:.4f}"
                  f"  nfev={de_result.nfev}")
    return results


# =============================================================================
# Report helpers
# =============================================================================
def dominant_rate_for_state(psi, profile):
    """Rate of the eigenmode with highest overlap (non-stationary modes only)."""
    rho0 = np.outer(psi, psi.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    probs = np.abs(c0) ** 2
    nonstat = np.abs(profile['rates']) > 1e-10
    if not nonstat.any():
        return float('nan'), float('nan')
    idxs = np.where(nonstat)[0]
    dom_i = idxs[np.argmax(probs[idxs])]
    return float(profile['eigvals'][dom_i].real), float(probs[dom_i])


def print_state_row(name, psi):
    rho = np.outer(psi, psi.conj())
    conc = max_adjacent_concurrence(rho, N)
    conc_mean = mean_adjacent_concurrence(rho, N)
    row = {'name': name, 'conc_max': conc, 'conc_mean': conc_mean}
    for label in ('sacrifice', 'uniform'):
        w = analyze_state(psi, PROFILES[label])
        dom_rate, dom_ov = dominant_rate_for_state(psi, PROFILES[label])
        row[label] = dict(slow=w['slow_wt'], osc=w['osc_wt'], xor=w['xor_wt'],
                          dom_rate=dom_rate, dom_ov=dom_ov)
    print(
        f"  {name:<18}  "
        f"sac[slow={row['sacrifice']['slow']:5.1f}% "
        f"osc={row['sacrifice']['osc']:5.1f}% "
        f"dom={row['sacrifice']['dom_rate']:+6.3f}/{row['sacrifice']['dom_ov']:.2f}]  "
        f"uni[slow={row['uniform']['slow']:5.1f}% "
        f"dom={row['uniform']['dom_rate']:+6.3f}/{row['uniform']['dom_ov']:.2f}]  "
        f"C_max={conc:.3f}"
    )
    return row


def show_top_basis_components(psi, top_k=6):
    """Return a list of (basis_index, amplitude) for the largest components of psi
    in the computational basis."""
    abs_amps = np.abs(psi)
    order = np.argsort(-abs_amps)[:top_k]
    return [(int(i), complex(psi[i])) for i in order if abs_amps[i] > 1e-4]


# =============================================================================
# Main
# =============================================================================
def pure_from_rho(rho):
    """Extract pure state from a rank-1 density matrix."""
    w, v = linalg.eigh(rho)
    k = int(np.argmax(w.real))
    return v[:, k] * np.exp(-1j * np.angle(v[0, k] if abs(v[0, k]) > 1e-10 else 1.0))


if __name__ == "__main__":
    import time
    t_start = time.time()

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "results", "optimal_state_n5_sacrifice")
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 90)
    print("OPTIMAL STATE SEARCH: N=5 sacrifice chain [80, 8, 79, 53, 85]")
    print("=" * 90)
    print(f"  gamma_sacrifice = {np.round(gamma_sacrifice, 5)}")
    print(f"  Sg (sacrifice)  = {Sg_sac:.5f}")
    print(f"  Sg (uniform)    = {Sg_uni:.5f}")
    print(f"  slow-mode cut   = (rates > 0 and rates < {Sg_sac:.3f}) in sacrifice profile")
    print(f"  subspace dim    = {M_SUB} (|0..0> + {N} single-excitations)")
    print(f"  free params     = {N_PARAMS}")

    # --- Baseline states --------------------------------------------------
    print("\n--- Baseline states ---")
    psi_bell = pure_from_rho(bell_plus_center(N))
    psi_w3c = pure_from_rho(w_center3_plus(N))
    psi_w5 = pure_from_rho(w_full(N))

    rows = []
    rows.append(print_state_row("Bell+|+>", psi_bell))
    rows.append(print_state_row("W3 center|+>", psi_w3c))
    rows.append(print_state_row("W5 full", psi_w5))

    # --- DE optimizer: low-excitation subspace ----------------------------
    basis_le1 = low_excitation_basis(N)
    print(f"\n{'=' * 90}")
    print(f"PHASE 1: Differential Evolution on low-excitation subspace "
          f"(dim={len(basis_le1)}, params={2*len(basis_le1)-1})")
    print("=" * 90)

    best_psi_le1, best_obj_le1, best_x_le1, de_res_le1 = run_de_optimizer(
        basis_le1, conc_floor=0.01, seed=42, maxiter=1000, popsize=20, verbose=True)

    w_le1 = analyze_state(best_psi_le1, SACRIFICE_PROFILE)
    rho_le1 = np.outer(best_psi_le1, best_psi_le1.conj())
    conc_le1 = max_adjacent_concurrence(rho_le1, N)
    print(f"\n  Best: slow_wt={w_le1['slow_wt']:.2f}%  osc_wt={w_le1['osc_wt']:.2f}%"
          f"  xor_wt={w_le1['xor_wt']:.2f}%  C_max={conc_le1:.4f}"
          f"  nfev={de_res_le1.nfev}")

    rows.append(print_state_row("OPT(le1)", best_psi_le1))

    print("\n  Top components:")
    for idx, amp in show_top_basis_components(best_psi_le1, top_k=8):
        bits = format(idx, f'0{N}b')
        print(f"    |{bits}>  ({idx:3d})  amp = {amp.real:+.6f} {amp.imag:+.6f}j"
              f"  |amp| = {abs(amp):.6f}")

    # --- DE optimizer: weight <= 2 subspace -------------------------------
    basis_le2 = weight_le2_basis(N)
    n_params_le2 = 2 * len(basis_le2) - 1
    print(f"\n{'=' * 90}")
    print(f"PHASE 2: Differential Evolution on weight-2 subspace "
          f"(dim={len(basis_le2)}, params={n_params_le2})")
    print("=" * 90)

    best_psi_le2, best_obj_le2, best_x_le2, de_res_le2 = run_de_optimizer(
        basis_le2, conc_floor=0.01, seed=42, maxiter=1000, popsize=20, verbose=True)

    w_le2 = analyze_state(best_psi_le2, SACRIFICE_PROFILE)
    rho_le2 = np.outer(best_psi_le2, best_psi_le2.conj())
    conc_le2 = max_adjacent_concurrence(rho_le2, N)
    print(f"\n  Best: slow_wt={w_le2['slow_wt']:.2f}%  osc_wt={w_le2['osc_wt']:.2f}%"
          f"  xor_wt={w_le2['xor_wt']:.2f}%  C_max={conc_le2:.4f}"
          f"  nfev={de_res_le2.nfev}")

    rows.append(print_state_row("OPT(le2)", best_psi_le2))

    print("\n  Top components:")
    for idx, amp in show_top_basis_components(best_psi_le2, top_k=10):
        bits = format(idx, f'0{N}b')
        print(f"    |{bits}>  ({idx:3d})  amp = {amp.real:+.6f} {amp.imag:+.6f}j"
              f"  |amp| = {abs(amp):.6f}")

    # --- Pareto frontier: low-excitation subspace -------------------------
    print(f"\n{'=' * 90}")
    print("PHASE 3: Pareto frontier (low-excitation subspace)")
    print("=" * 90)
    pareto_le1 = pareto_frontier(basis_le1, floors=(0.01, 0.05, 0.1, 0.2, 0.3),
                                  seed=42, verbose=True)

    print(f"\n  {'floor':>6}  {'slow_wt':>8}  {'osc_wt':>8}  {'C_max':>8}  {'nfev':>6}")
    print(f"  {'-' * 45}")
    for p in pareto_le1:
        print(f"  {p['floor']:>6.2f}  {p['slow_wt']:>7.2f}%  {p['osc_wt']:>7.2f}%"
              f"  {p['conc']:>8.4f}  {p['nfev']:>6d}")

    # --- Pareto frontier: weight-2 subspace -------------------------
    print(f"\n{'=' * 90}")
    print("PHASE 4: Pareto frontier (weight-2 subspace)")
    print("=" * 90)
    pareto_le2 = pareto_frontier(basis_le2, floors=(0.01, 0.05, 0.1, 0.2, 0.3),
                                  seed=42, verbose=True)

    print(f"\n  {'floor':>6}  {'slow_wt':>8}  {'osc_wt':>8}  {'C_max':>8}  {'nfev':>6}")
    print(f"  {'-' * 45}")
    for p in pareto_le2:
        print(f"  {p['floor']:>6.2f}  {p['slow_wt']:>7.2f}%  {p['osc_wt']:>7.2f}%"
              f"  {p['conc']:>8.4f}  {p['nfev']:>6d}")

    # --- Summary table ----------------------------------------------------
    print(f"\n{'=' * 90}")
    print("FINAL COMPARISON (sacrifice profile)")
    print("=" * 90)
    print(f"  {'State':<18}  {'slow_wt':>8}  {'osc_wt':>8}  {'xor_wt':>8}"
          f"  {'dom_rate':>10}  {'C_max':>8}")
    print(f"  {'-' * 72}")
    for r in rows:
        sac = r['sacrifice']
        print(f"  {r['name']:<18}  {sac['slow']:>7.2f}%  {sac['osc']:>7.2f}%"
              f"  {sac['xor']:>7.2f}%  {sac['dom_rate']:>+10.4f}"
              f"  {r['conc_max']:>8.4f}")

    t_end = time.time()
    print(f"\nTotal runtime: {t_end - t_start:.1f}s")

    # --- Persist results --------------------------------------------------
    out_path = os.path.join(out_dir, "optimal_state_n5_sacrifice.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("Optimal state search for N=5 sacrifice chain\n")
        f.write(f"Computed: {np.datetime64('today')}\n")
        f.write(f"Runtime: {t_end - t_start:.1f}s\n\n")
        f.write(f"Chain:          [80, 8, 79, 53, 85]\n")
        f.write(f"Sg (sacrifice): {Sg_sac:.5f}\n")
        f.write(f"Sg (uniform):   {Sg_uni:.5f}\n")
        f.write(f"Slow-mode def:  (rates > 1e-6) & (rates < Sg - 1e-6)\n\n")

        f.write("=" * 70 + "\n")
        f.write("COMPARISON TABLE (sacrifice profile)\n")
        f.write("=" * 70 + "\n")
        f.write(f"  {'State':<18}  {'slow_wt':>8}  {'osc_wt':>8}  {'xor_wt':>8}"
                f"  {'dom_rate':>10}  {'dom_ov':>8}  {'C_max':>8}\n")
        for r in rows:
            sac = r['sacrifice']
            f.write(f"  {r['name']:<18}  {sac['slow']:>7.2f}%  {sac['osc']:>7.2f}%"
                    f"  {sac['xor']:>7.2f}%  {sac['dom_rate']:>+10.4f}"
                    f"  {sac['dom_ov']:>8.4f}  {r['conc_max']:>8.4f}\n")

        f.write(f"\n{'=' * 70}\n")
        f.write("UNIFORM PROFILE (same states)\n")
        f.write("=" * 70 + "\n")
        f.write(f"  {'State':<18}  {'slow_wt':>8}  {'dom_rate':>10}  {'dom_ov':>8}\n")
        for r in rows:
            uni = r['uniform']
            f.write(f"  {r['name']:<18}  {uni['slow']:>7.2f}%  {uni['dom_rate']:>+10.4f}"
                    f"  {uni['dom_ov']:>8.4f}\n")

        # Optimal state details
        for tag, psi, de_res in [("OPT(le1) low-excitation", best_psi_le1, de_res_le1),
                                  ("OPT(le2) weight-2", best_psi_le2, de_res_le2)]:
            f.write(f"\n{'=' * 70}\n")
            f.write(f"OPTIMAL STATE: {tag}\n")
            f.write("=" * 70 + "\n")
            w = analyze_state(psi, SACRIFICE_PROFILE)
            rho = np.outer(psi, psi.conj())
            conc = max_adjacent_concurrence(rho, N)
            f.write(f"  Method:       differential_evolution\n")
            f.write(f"  slow_wt:      {w['slow_wt']:.4f}%\n")
            f.write(f"  osc_wt:       {w['osc_wt']:.4f}%\n")
            f.write(f"  xor_wt:       {w['xor_wt']:.4f}%\n")
            f.write(f"  C_max:        {conc:.6f}\n")
            f.write(f"  nfev:         {de_res.nfev}\n")
            f.write(f"  Top components:\n")
            for idx, amp in show_top_basis_components(psi, top_k=10):
                bits = format(idx, f'0{N}b')
                f.write(f"    |{bits}>  amp = {amp.real:+.8f} {amp.imag:+.8f}j"
                        f"  |amp| = {abs(amp):.8f}\n")

        # Pareto frontiers
        for tag, pareto in [("Low-excitation (Hamming <= 1)", pareto_le1),
                             ("Weight-2 (Hamming <= 2)", pareto_le2)]:
            f.write(f"\n{'=' * 70}\n")
            f.write(f"PARETO FRONTIER: {tag}\n")
            f.write("=" * 70 + "\n")
            f.write(f"  {'floor':>6}  {'slow_wt':>8}  {'osc_wt':>8}  {'xor_wt':>8}"
                    f"  {'C_max':>8}  {'nfev':>6}\n")
            for p in pareto:
                f.write(f"  {p['floor']:>6.2f}  {p['slow_wt']:>7.2f}%  {p['osc_wt']:>7.2f}%"
                        f"  {p['xor_wt']:>7.2f}%  {p['conc']:>8.4f}  {p['nfev']:>6d}\n")
            # Top components for each Pareto point
            for p in pareto:
                f.write(f"\n  floor={p['floor']:.2f}: slow_wt={p['slow_wt']:.2f}%"
                        f"  C_max={p['conc']:.4f}\n")
                for idx, amp in show_top_basis_components(p['psi'], top_k=6):
                    bits = format(idx, f'0{N}b')
                    f.write(f"    |{bits}>  |amp| = {abs(amp):.6f}\n")

    print(f"\nSaved to {out_path}")
