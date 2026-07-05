"""Net-cooling (R2 proxy) readout on EXISTING hardware data (the bridge, part 3).

Zero QPU, no Lindblad fit; companion of experiments/F81_VIOLATION_HARDWARE_BRIDGE.md.

The grounded object (confirmed in f81_identity_velocity_grounding.py): the
f81 discriminator is 2^(N-1) x the RMS identity-escape velocity; per site
that velocity is the net cooling flux a_l = gamma_down,l - gamma_up,l, and
on a T1 leg a_l = z_inf,l * b_l with b_l the leg rate (1/T1_total,l) and
z_inf,l the free asymptote.

This scout computes that number on the price_pair Marrakesh campaign's
Block B (|1> free decay, readout-mitigated <Z_l>(t), 10 delays), runs 1+2
(lines [2,3,4] and [93,94,95], 2026-07-04). NO Lindblad channel fit anywhere:
the inputs are two directly measured decay parameters per qubit. This is the
non-tautological version R2 of the bridge, on data already paid for.

Structural notes (checked before trusting the model):
  - idle Hamiltonian on these lines = detunings + always-on ZZ, both DIAGONAL:
    populations decouple from every coherent diagonal term, so <Z_l>(t) obeys
    the pure local T1 Bloch equation dz/dt = a_l - b_l z. The R2 leg is
    ZZ-immune (and ZZ is Pi^2-even anyway). Condition to state: this needs a
    diagonal idle H; an exchange term would mix site populations.
  - the asymptote is extrapolated from t_max << T1: error bars are propagated
    from the fit covariance and stated. This is a method demo with honest
    uncertainty, not a precision measurement.

Run: python -X utf8 simulations/f84_net_cooling_readout_marrakesh.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

REPO_ROOT = Path(__file__).resolve().parents[1]
os.chdir(REPO_ROOT)

RUNS = [
    ("run 1, line [2,3,4]", "data/ibm_price_pair_july2026/price_pair_ibm_marrakesh_20260704_073922.json"),
    ("run 2, line [93,94,95]", "data/ibm_price_pair_july2026/price_pair_ibm_marrakesh_20260704_075304.json"),
]


def bloch_decay(t, z_inf, dz, b):
    return z_inf + dz * np.exp(-b * t)


def main():
    print("=" * 92)
    print("S4 R2 proxy: f81_violation (local sigma+- family) from price_pair Block B, zero QPU")
    print("=" * 92)
    print()
    print("  Recipe per qubit: fit z(t) = z_inf + dz*exp(-b t) (free asymptote);")
    print("  identity-escape velocity a = z_inf * b; violation_R2 = 2^(N-1) * sqrt(sum_l a_l^2).")
    print()

    for label, path in RUNS:
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
        B = d["dataset"]["B"]
        eps = d["dataset"]["readout_eps"]
        shots = d["shots"]
        t1_pipeline = d["analysis"]["t1_us"]
        qubits = d["path"]
        t = np.array([s["t_us"] for s in B])
        Zs = np.array([s["z"] for s in B])  # shape (n_t, 3)
        n_q = Zs.shape[1]

        print(f"  {label}  (job {d['job_id']}, shots {shots})")
        print(f"    delays: {t.min():.1f} .. {t.max():.1f} us, {len(t)} points")

        a_list, a_err = [], []
        for q in range(n_q):
            z = Zs[:, q]
            # shot noise, amplified by readout mitigation (1/(1-e0-e1) lever arm)
            e0, e1 = eps[q]
            sig = (1.0 / np.sqrt(shots)) / max(1.0 - e0 - e1, 0.5)
            sigma = np.full_like(z, sig)
            p0 = [1.0, -2.0, 1.0 / max(t1_pipeline[q], 1.0)]
            popt, pcov = curve_fit(bloch_decay, t, z, p0=p0, sigma=sigma,
                                   absolute_sigma=True, maxfev=20000)
            z_inf, dz, b = popt
            perr = np.sqrt(np.diag(pcov))
            a = z_inf * b
            # error propagation incl. covariance between z_inf and b
            var_a = (b ** 2) * pcov[0, 0] + (z_inf ** 2) * pcov[2, 2] + 2 * z_inf * b * pcov[0, 2]
            sa = np.sqrt(max(var_a, 0.0))
            a_list.append(a)
            a_err.append(sa)
            T1_leg = 1.0 / b if b > 0 else float("inf")
            print(f"    q{qubits[q]}: z_inf = {z_inf:+.3f} ± {perr[0]:.3f},  b = {b:.5f} ± {perr[2]:.5f} /us"
                  f"  (T1_leg = {T1_leg:.0f} us; pipeline T1 = {t1_pipeline[q]:.0f} us)")
            print(f"          a = z_inf*b = {a:+.5f} ± {sa:.5f} /us   (naive 1/T1 = {1.0/t1_pipeline[q]:.5f})")

        a_arr = np.array(a_list)
        sa_arr = np.array(a_err)
        N = n_q
        viol = (2 ** (N - 1)) * float(np.sqrt(np.sum(a_arr ** 2)))
        # d viol / d a_l = 2^(N-1) * a_l / sqrt(sum a^2)
        dv = (2 ** (N - 1)) * np.sqrt(np.sum((a_arr * sa_arr) ** 2)) / np.sqrt(np.sum(a_arr ** 2))
        g_rms = float(np.sqrt(np.mean(a_arr ** 2)))
        naive = (2 ** (N - 1)) * float(np.sqrt(np.sum((1.0 / np.array(t1_pipeline)) ** 2)))
        print(f"    -> violation_R2 (N=3 line) = {viol:.5f} ± {dv:.5f}   [naive all-z_inf=1: {naive:.5f}]")
        print(f"    -> identity-escape RMS velocity = {g_rms:.5f} /us  (net-cooling reading, "
              f"T1-equivalent {1.0/g_rms:.0f} us)")
        print()


if __name__ == "__main__":
    main()
