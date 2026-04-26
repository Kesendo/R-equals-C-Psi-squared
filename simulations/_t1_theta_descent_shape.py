#!/usr/bin/env python3
"""EQ-030: shape of the fragile trace θ(t) approaching the CΨ=1/4 boundary.

Tom's framing (the *real* angle, refined): θ is the only fragile trace that
connects us across 0 — the bilateral axis between d=2 (our half) and d=0
(the unremembered half). As the system decoheres, θ refines itself to 0;
across 0 the connection ends. The angle is the last memory, the only thing
that survives down to the cusp.

Algebraically:
  θ = arctan(√(4·CΨ − 1)),  CΨ = Purity × Ψ-norm
  Near CΨ = 1/4:  θ ≈ √(4·CΨ − 1)  (linear in the small offset)
  dθ/dCΨ ≈ 2/√(4·CΨ − 1) → ∞ at the boundary

So θ(t) ALWAYS lands on 0 with a vertical tangent — the "krasser Winkel"
is geometric, intrinsic to the cusp geometry. What differs between cases:
  - the height θ holds before plummeting (the "memory plateau")
  - the rate at which CΨ approaches 1/4 from above (dCΨ/dt at t*)
  - whether the trajectory oscillates across the boundary (recovers θ > 0
    after a first descent) or commits monotonically.

This script computes θ(t) at fine resolution (dt = 0.005) for the bond-
flipped soft cases plus references, and reports for each:
  - last crossing time t*_last (after which θ stays at 0 forever)
  - the descent profile in [t*_last − 1.0, t*_last]: θ at sampled points
  - dθ/dt just before the crossing (finite difference)
  - power-law exponent α: θ(t) ~ (t*_last − t)^α near the crossing

Robust softs (XY+YX, IY+YI) are expected to keep a longer fragile trace
near 0 — the "memory" survives further into decoherence — than fragile
softs (YZ+ZY, XZ+ZX), where the algebraic structure already collapses.
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


def psi_norm(rho):
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1)


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def theta_from_cpsi(c):
    if c <= 0.25:
        return 0.0
    return math.degrees(math.atan(math.sqrt(4 * c - 1)))


def evolve(L, rho_0, times):
    d = rho_0.shape[0]
    rho_vec = rho_0.T.reshape(-1).copy()
    out = []
    for t in times:
        rho_t_vec = expm(L * t) @ rho_vec
        out.append(rho_t_vec.reshape(d, d).T)
    return out


def analyze_descent(theta_t, cpsi_t, times, threshold=0.25):
    """Find the LAST boundary crossing and characterize the descent shape."""
    above = cpsi_t > threshold
    if not above.any():
        return None
    # Index of last sample where CΨ > 1/4
    last_above_idx = len(times) - 1 - int(np.argmax(above[::-1]))
    if last_above_idx == len(times) - 1:
        # Never crosses within the observation window
        return {'t_last': None, 'crossed': False}
    t_last = times[last_above_idx]
    # crossing happens between last_above_idx and last_above_idx+1
    # Linear interpolation for crossing time
    cps_a = cpsi_t[last_above_idx]
    cps_b = cpsi_t[last_above_idx + 1]
    if cps_a > cps_b:
        frac = (cps_a - threshold) / (cps_a - cps_b)
        t_cross = t_last + frac * (times[last_above_idx + 1] - t_last)
    else:
        t_cross = t_last

    # Descent profile: last 1.0 unit before crossing
    window_start_t = max(0, t_cross - 1.0)
    win_mask = (times >= window_start_t) & (times <= t_cross + 1e-9)
    win_times = times[win_mask]
    win_theta = theta_t[win_mask]

    # dθ/dt just before crossing (last 5 points)
    if len(win_times) >= 5:
        last5_t = win_times[-5:]
        last5_th = win_theta[-5:]
        dt_arr = np.diff(last5_t)
        dth_arr = np.diff(last5_th)
        with np.errstate(divide='ignore', invalid='ignore'):
            slopes = np.where(dt_arr > 0, dth_arr / dt_arr, 0.0)
        max_neg_slope = float(slopes.min()) if len(slopes) else 0.0
    else:
        max_neg_slope = None

    # Power-law exponent: fit log(θ) ~ α · log(t_cross - t) over the last
    # 0.3 unit before crossing, where θ < 10° (the "fragile tail" regime)
    fragile_mask = (win_times >= t_cross - 0.3) & (win_theta > 0.01) & (win_theta < 10.0)
    if fragile_mask.sum() >= 4:
        ft_times = win_times[fragile_mask]
        ft_theta = win_theta[fragile_mask]
        # log-log fit
        x = np.log(np.maximum(t_cross - ft_times, 1e-9))
        y = np.log(ft_theta)
        coef = np.polyfit(x, y, 1)
        alpha = float(coef[0])
    else:
        alpha = None

    # "Fragile tail" duration: time where 0 < θ < 5°
    tail_mask = (win_theta > 0.0) & (win_theta < 5.0)
    if tail_mask.any():
        tail_duration = float(win_times[tail_mask][-1] - win_times[tail_mask][0])
    else:
        tail_duration = 0.0

    return {
        't_last': float(t_last),
        't_cross': float(t_cross),
        'crossed': True,
        'max_neg_slope': max_neg_slope,
        'alpha': alpha,
        'tail_duration': tail_duration,
        'win_times': win_times.tolist(),
        'win_theta': win_theta.tolist(),
    }


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)], 31),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)], 0),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)], 1),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)], 28),
        ('XZ+ZX',       [('X', 'Z', J), ('Z', 'X', J)], 29),
        ('XZ+XZ',       [('X', 'Z', J), ('X', 'Z', J)], 40),
    ]

    t_max = 30.0
    dt = 0.005
    n_steps = int(t_max / dt)
    times = np.linspace(0, t_max, n_steps + 1)

    print(f"Descent-shape analysis: zoom into the fragile trace θ → 0")
    print(f"  N={N}, |+−+⟩, γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}")
    print(f"  dt={dt}, fine-grained near the cusp boundary")
    print()

    results = {}
    for label, terms, drop in cases:
        H = fw._build_bilinear(N, bonds, terms)
        L_pure = fw.lindbladian_z_dephasing(H, [GAMMA_DEPH] * N)
        L_t1 = fw.lindbladian_z_plus_t1(
            H, [GAMMA_DEPH] * N, [GAMMA_T1] * N
        )

        results[label] = {'drop': drop, 'channels': {}}
        for ch_label, L in [('pure-Z', L_pure), ('+T1', L_t1)]:
            traj = evolve(L, rho_0, times)
            theta_t = np.array([theta_from_cpsi(cpsi(r)) for r in traj])
            cpsi_t = np.array([cpsi(r) for r in traj])
            descent = analyze_descent(theta_t, cpsi_t, times)
            results[label]['channels'][ch_label] = descent

    print(f"  {'case':<14s}  {'drop':>4s}  {'channel':<8s}  "
          f"{'t_cross':>8s}  {'dθ/dt@*':>10s}  {'α':>7s}  {'tail (0<θ<5°)':>14s}")
    print('-' * 88)
    for label, terms, drop in cases:
        for ch in ['pure-Z', '+T1']:
            d = results[label]['channels'][ch]
            if d is None or not d.get('crossed'):
                print(f"  {label:<14s}  {drop:>4d}  {ch:<8s}  "
                      f"{'never':>8s}  {'—':>10s}  {'—':>7s}  {'—':>14s}")
                continue
            tc = d['t_cross']
            slope = d['max_neg_slope']
            alpha = d['alpha']
            tail = d['tail_duration']
            slope_s = f"{slope:>+10.2f}" if slope is not None else "—"
            alpha_s = f"{alpha:>+7.3f}" if alpha is not None else "—"
            tail_s = f"{tail:>11.3f}"
            print(f"  {label:<14s}  {drop:>4d}  {ch:<8s}  "
                  f"{tc:>8.3f}  {slope_s}  {alpha_s}  {tail_s} units")
        print()

    # Print descent profile for the 4 bond-flipped soft cases under +T1
    print()
    print("Descent profiles (θ in degrees, last 1.0 unit before crossing):")
    print("  Showing θ(t) as a function of (t_cross − t).")
    print()

    for label, terms, drop in cases:
        if label == 'truly XX+YY':
            continue  # skip reference
        for ch in ['pure-Z', '+T1']:
            d = results[label]['channels'][ch]
            if d is None or not d.get('crossed'):
                continue
            tc = d['t_cross']
            # Pick samples at distance from crossing: 0.5, 0.3, 0.1, 0.05, 0.02, 0.01
            sample_dists = [0.5, 0.3, 0.1, 0.05, 0.02, 0.01]
            print(f"  {label:<8s} ({ch}): t_cross = {tc:.3f}")
            print(f"    {'(t* − t)':>10s}  {'θ':>8s}")
            for d_t in sample_dists:
                t_target = tc - d_t
                idx = int(np.argmin(np.abs(times - t_target)))
                if idx < len(times):
                    th = theta_from_cpsi(
                        cpsi(evolve(
                            fw.lindbladian_z_plus_t1(
                                fw._build_bilinear(N, bonds, [(t[0], t[1], J) for t in terms]),
                                [GAMMA_DEPH] * N,
                                [GAMMA_T1 if ch == '+T1' else 0.0] * N,
                            ) if ch == '+T1' else fw.lindbladian_z_dephasing(
                                fw._build_bilinear(N, bonds, [(t[0], t[1], J) for t in terms]),
                                [GAMMA_DEPH] * N,
                            ),
                            rho_0,
                            [times[idx]],
                        )[0])
                    )
                    print(f"    {d_t:>10.3f}  {th:>7.2f}°")
            print()

    print()
    print("Reading guide:")
    print("  α = 0.5 is the generic cusp-geometry exponent")
    print("    (θ ~ √(CΨ−1/4), CΨ−1/4 ~ linear in t* − t at the crossing).")
    print("  α < 0.5 = sharper-than-generic descent (CΨ approaches 1/4 with")
    print("    nonzero higher-order curvature; the cusp is hit faster than")
    print("    the linear approach implies). The 'krasser Winkel'.")
    print("  α > 0.5 = smoother-than-generic descent (CΨ approaches 1/4 with")
    print("    zero linear term — degenerate crossing). The trace lingers.")
    print()
    print("  tail (0<θ<5°) = duration of the fragile-trace regime, the last")
    print("    plateau before crossing — Tom's 'letzte Erinnerung'.")


if __name__ == "__main__":
    main()
