#!/usr/bin/env python3
"""_f71_nonuniform_gamma_verification.py

F71 open question #1 -- "non-uniform gamma_i: site-dependent dephasing".

Witnesses F101: the c_1 bond-mirror deviation

    D(b) := c_1(b) - c_1(N-2-b)

is an exactly ODD function of the F71-anti-palindromic component of the per-site
Z-dephasing profile gamma, at fixed F71-palindromic component gamma_sym. F101 is
the observable-side twin of F91 (gamma spectral invariance); the J-side twin is
F100. Decompose gamma (N sites) at the site-mirror l <-> N-1-l:

    gamma_sym  = (gamma + F71(gamma)) / 2     (F71-palindromic)
    gamma_anti = (gamma - F71(gamma)) / 2     (F71-anti-palindromic)

Write the base profile gamma_base(s) = gamma_sym + s * gamma_anti_dir. Prediction:

  (a) D(b; s=0) = 0                  -- F71 survives ALL palindromic non-
                                        uniformity, incl. non-uniform gamma_sym.
  (b) D(b; -s) = -D(b; +s)           -- exactly odd in s, all orders.
  (c) D(b; s) = kappa_gamma * s + O(s^3) -- leading order linear; even powers
                                        vanish.
  (d) kappa_gamma generically DEPENDS on gamma_sym (Tier 2 empirical, no closed
      form). This run maps the gamma_sym-dependence.

c_1 is the EQ-018 closure-breaking coefficient (PROOF_C1_MIRROR_SYMMETRY): the
first-order response of Sum_i ln(alpha_i) to a single-bond dJ probe, extracted
by the alpha-rescaling fit on per-site purity. The base swept here is gamma; the
c_1-defining probe is still a dJ probe, with J held uniform.

EXACT SECTOR RESTRICTION. The chain XY + Z-dephasing Liouvillian L is exactly
block-diagonal in the (bra-excitation, ket-excitation) bigrading. The probe
states psi_k+vac live in the (popcount<=1) x (popcount<=1) operator block of
dimension (N+1)^2. build_L_sub_gamma propagates inside that block; it is
verified bit-identical to the full per-site-gamma Liouvillian sliced to the
block by the Gate-1 self-test.

PROBE STATES: psi_1+vac and psi_2+vac (the PROOF_C1-validated reflection-
symmetric states).

Usage:
  python -u _f71_nonuniform_gamma_verification.py --N 3,4,5
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (  # noqa: E402
    T_FINAL, N_STEPS, T_FIT_MAX, Z,
    build_H_XY, site_op,
    vacuum_ket, single_excitation_mode,
    density_matrix, per_site_purity, fit_alpha,
)

RESULTS_DIR = Path(__file__).parent / "results" / "f71_nonuniform_gamma_verification"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

J_UNIFORM_VAL = 1.0                     # fixed uniform coupling; c_1 probe perturbs it
DJ_EXTRACT = 0.01                       # single-bond probe that defines c_1
S_VALUES = (-0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03)  # gamma_anti sweep amplitudes
FLAT_TOL = 1e-9                         # per-site purity span below -> no signal
C1_INFORMATIVE_FLOOR = 1e-3             # |c_1| above this -> record is informative


# ---------------------------------------------------------------------------
# Exact sector restriction: the (N+1)^2-dim invariant subspace of psi_k+vac
# ---------------------------------------------------------------------------
def kept_basis(N):
    """Computational basis states with popcount <= 1: the vacuum plus the N
    single-excitation states. span(S) carries all of psi_k+vac's dynamics."""
    d = 1 << N
    return [a for a in range(d) if bin(a).count("1") <= 1]


def build_L_sub_gamma(J_list, gamma_list, N, S):
    """Liouvillian restricted to the (popcount<=1) x (popcount<=1) operator block,
    with a PER-SITE dephasing profile gamma_list (length N) instead of a scalar.

    D[rho] = sum_l gamma_list[l] (Z_l rho Z_l - rho). Verified bit-identical to
    build_L_full_gamma sliced to S x S by the Gate-1 self-test."""
    H_S = build_H_XY(J_list, N)[np.ix_(S, S)]
    m = len(S)
    I_m = np.eye(m, dtype=complex)
    L = -1j * (np.kron(I_m, H_S) - np.kron(H_S.T, I_m))
    for l in range(N):
        Zl = site_op(Z, l, N)[np.ix_(S, S)]          # Z_l restricted to span(S)
        L += gamma_list[l] * (np.kron(Zl.T, Zl) - np.kron(I_m, I_m))
    return L


def build_L_full_gamma(J_list, gamma_list, N):
    """Full 4^N Liouvillian with a per-site dephasing profile gamma_list. The
    Gate-1 reference for build_L_sub_gamma."""
    d = 1 << N
    H = build_H_XY(J_list, N)
    I_d = np.eye(d, dtype=complex)
    L = -1j * (np.kron(I_d, H) - np.kron(H.T, I_d))
    for l in range(N):
        Zl = site_op(Z, l, N)
        L += gamma_list[l] * (np.kron(Zl.T, Zl) - np.kron(I_d, I_d))
    return L


def selftest_sector_restriction():
    """Gate 1 -- exact sector restriction self-test: build_L_sub_gamma MUST equal
    the full 4^N per-site-gamma Liouvillian sliced to the (popcount<=1) operator
    block, bit-level. Hard-asserts at N=3,4."""
    print("Gate 1 -- exact sector restriction self-test:", flush=True)
    for N in (3, 4):
        d = 1 << N
        S = kept_basis(N)
        kept = [S[m] + d * S[n] for n in range(len(S)) for m in range(len(S))]
        J = [J_UNIFORM_VAL] * (N - 1)
        for label, gamma in (("uniform", [0.05] * N),
                              ("nonuniform", [0.03 + 0.02 * l for l in range(N)])):
            L_full = build_L_full_gamma(J, gamma, N)
            L_sliced = L_full[np.ix_(kept, kept)]
            L_direct = build_L_sub_gamma(J, gamma, N, S)
            diff = float(np.max(np.abs(L_sliced - L_direct)))
            assert diff < 1e-12, (f"Gate 1 FAIL N={N} {label}: "
                                  f"max|L_sub - L_full[kept]| = {diff:.2e}")
            print(f"  OK  N={N} {label:11s} dim {len(S)**2:3d} vs {d*d:5d}  "
                  f"max|L_sub - sliced full L| = {diff:.2e}", flush=True)
    print(flush=True)


# ---------------------------------------------------------------------------
# c_1 pipeline (propagation inside the sector, alpha-rescaling fit + guard)
# ---------------------------------------------------------------------------
def eig_and_inv(L):
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    return eigvals, V_R, V_Linv


def propagate(dec, rho_0, times, S, N):
    """Propagate inside the (N+1)^2-dim sector, return the full d x d
    trajectory (the non-S entries of rho(t) are exactly zero)."""
    eigvals, V_R, V_Linv = dec
    d = 1 << N
    m = len(S)
    Sij = np.ix_(S, S)
    c0 = V_Linv @ rho_0[Sij].flatten(order='F')
    out = np.zeros((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        out[k][Sij] = (V_R @ (np.exp(eigvals * t) * c0)).reshape(m, m, order='F')
    return out


def measure_c1(rho_0, decA, decBp, decBm, N, times, S):
    """EQ-018 closure-breaking coefficient c_1 via the alpha-rescaling fit.

    Flat-site guard: a site whose per-site purity P_A(i,.) is constant over the
    fit window (span < FLAT_TOL) carries no closure signal -> alpha := 1.
    Returns (c_1, max_fit_rmse, n_guarded_sites)."""
    P_A = per_site_purity(propagate(decA, rho_0, times, S, N), N)
    P_Bp = per_site_purity(propagate(decBp, rho_0, times, S, N), N)
    P_Bm = per_site_purity(propagate(decBm, rho_0, times, S, N), N)
    mask = times <= T_FIT_MAX
    ln_p = ln_m = 0.0
    rmse = 0.0
    n_guarded = 0
    for i in range(N):
        if float(np.ptp(P_A[mask, i])) < FLAT_TOL:
            n_guarded += 1            # alpha = 1 -> ln(alpha) = 0 contribution
            continue
        ap, rp = fit_alpha(times, P_A[:, i], P_Bp[:, i])
        am, rm = fit_alpha(times, P_A[:, i], P_Bm[:, i])
        ln_p += np.log(ap); ln_m += np.log(am)
        rmse = max(rmse, rp, rm)
    c_1 = (ln_p - ln_m) / (2 * DJ_EXTRACT)
    return float(c_1), float(rmse), int(n_guarded)


# ---------------------------------------------------------------------------
# Probe states / gamma_sym / gamma_anti construction
# ---------------------------------------------------------------------------
def probe_states(N):
    """PROOF_C1-validated reflection-symmetric probe states psi_k + vacuum.
    psi_1 is nodeless, psi_2 the first excited OBC sine mode. Both R-eigenmodes,
    so psi_k+vac is reflection-symmetric in the purity sense."""
    states = {}
    for k in (1, 2):
        ket = (vacuum_ket(N) + single_excitation_mode(N, k)) / np.sqrt(2.0)
        states[f"psi{k}+vac"] = density_matrix(ket)
    return states


def f71_mirror(profile):
    """F71 chain-mirror: reverse a profile array. Maps a per-site profile under
    site l <-> N-1-l (used here for gamma) or a per-bond profile under bond
    b <-> N-2-b. The reversal [::-1] is the same operation for either."""
    return np.array(profile)[::-1]


def gamma_anti_direction(N):
    """The canonical anti-palindromic 'linear ramp' direction on N sites:
    gamma_anti_dir[l] = 2l/(N-1) - 1. Anti-palindromic under the site-mirror
    l <-> N-1-l; for odd N the central site l=(N-1)/2 is F71-fixed (value 0)."""
    return np.array([2.0 * l / (N - 1) - 1.0 for l in range(N)])


def gamma_sym_profiles(N):
    """Four F71-palindromic per-site base profiles: 3 uniform magnitudes (probe
    the gamma_sym-magnitude dependence of kappa_gamma) + 1 non-uniform palindromic
    'valley' (test D=0 survival for non-uniform palindromic gamma; probe shape
    dependence). With |s| <= 0.03 and |gamma_anti_dir| <= 1, gamma_base stays
    strictly positive: min gamma = 0.05 - 0.03 = 0.02 > 0."""
    profiles = {
        "uniform_0.05": np.full(N, 0.05),
        "uniform_0.08": np.full(N, 0.08),
        "uniform_0.11": np.full(N, 0.11),
    }
    valley = np.array([0.11 - 0.05 * (1.0 - abs(2.0 * l / (N - 1) - 1.0))
                       for l in range(N)])
    profiles["valley_nonuniform"] = valley
    return profiles


def assert_decomposition(N):
    """Sanity: gamma_sym profiles palindromic, gamma_anti direction anti-
    palindromic (site-mirror l <-> N-1-l)."""
    anti = gamma_anti_direction(N)
    assert np.allclose(anti, -f71_mirror(anti)), "gamma_anti not anti-palindromic"
    for name, sym in gamma_sym_profiles(N).items():
        assert np.allclose(sym, f71_mirror(sym)), f"gamma_sym '{name}' not palindromic"


# ---------------------------------------------------------------------------
# Per-N run
# ---------------------------------------------------------------------------
def run_one_N(N, s_values):
    print(f"\n{'=' * 74}", flush=True)
    print(f"N = {N}   (full d^2 = {4**N}, sector dim = {(N+1)**2})   "
          f"sites 0..{N-1}, bonds 0..{N-2}", flush=True)
    print(f"{'=' * 74}", flush=True)
    assert_decomposition(N)
    nb = N - 1
    S = kept_basis(N)
    anti = gamma_anti_direction(N)
    sym_profiles = gamma_sym_profiles(N)
    probe = probe_states(N)
    state_names = list(probe.keys())
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    pairs = [(b, N - 2 - b) for b in range(nb) if b < N - 2 - b]
    J_uniform = [J_UNIFORM_VAL] * nb
    for nm, r in probe.items():                         # support check
        leak = float(abs(np.sum(np.abs(r)) - np.sum(np.abs(r[np.ix_(S, S)]))))
        assert leak < 1e-12, f"{nm} leaks outside the (N+1)^2 subspace: {leak}"
    print(f"  probe states: {state_names}", flush=True)
    print(f"  gamma_anti direction (linear ramp): "
          f"[{', '.join(f'{x:+.3f}' for x in anti)}]", flush=True)
    print(f"  F71 bond pairs for D: {pairs}"
          f"{'  central self-paired bond ' + str(nb // 2) if nb % 2 else ''}",
          flush=True)

    c1 = {name: {s: {b: {} for b in range(nb)} for s in s_values}
          for name in sym_profiles}
    rmse_max = 0.0
    guarded_total = 0
    n_profiles = len(sym_profiles) * len(s_values)
    done = 0
    t_run = time.time()

    for sym_name, sym in sym_profiles.items():
        for s in s_values:
            gamma_base = list(sym + s * anti)
            done += 1
            t0 = time.time()
            assert min(gamma_base) > 0.0, (
                f"gamma_base went non-positive: {sym_name} s={s}: {gamma_base}")
            decA = eig_and_inv(build_L_sub_gamma(J_uniform, gamma_base, N, S))
            for b in range(nb):
                J_bp = list(J_uniform); J_bp[b] += DJ_EXTRACT
                J_bm = list(J_uniform); J_bm[b] -= DJ_EXTRACT
                decBp = eig_and_inv(build_L_sub_gamma(J_bp, gamma_base, N, S))
                decBm = eig_and_inv(build_L_sub_gamma(J_bm, gamma_base, N, S))
                for st_name, st_rho in probe.items():
                    val, rm, ng = measure_c1(st_rho, decA, decBp, decBm,
                                             N, times, S)
                    c1[sym_name][s][b][st_name] = val
                    rmse_max = max(rmse_max, rm)
                    guarded_total += ng
            print(f"  [{done:2d}/{n_profiles}] gamma_sym={sym_name:18s} s={s:+.2f}  "
                  f"{time.time()-t0:5.1f} s", flush=True)

    print(f"  N={N} sweep done in {time.time()-t_run:.1f} s  "
          f"(max alpha-fit RMSE {rmse_max:.1e}; {guarded_total} flat-site guards "
          f"-- expect 0 for psi_k+vac)", flush=True)

    analysis = analyse_N(N, c1, s_values, pairs, state_names)
    return {
        "N": N, "j_uniform": J_UNIFORM_VAL, "dJ_extract": DJ_EXTRACT,
        "sector_dim": (N + 1) ** 2, "full_dim": 4 ** N,
        "s_values": list(s_values), "probe_states": state_names,
        "gamma_anti_direction": anti.tolist(),
        "gamma_sym_profiles": {k: v.tolist() for k, v in sym_profiles.items()},
        "c1": {name: {f"{s:+.2f}": {f"bond{b}": dict(c1[name][s][b])
                                    for b in range(nb)}
                      for s in s_values}
               for name in sym_profiles},
        "max_alpha_fit_rmse": rmse_max,
        "n_flat_site_guards": guarded_total,
        "analysis": analysis,
    }


def analyse_N(N, c1, s_values, pairs, state_names):
    """Per-N verdict: palindromic survival, oddness, leading coefficient,
    even-power suppression, gamma_sym-dependence of kappa."""
    s_arr = np.array(s_values)
    pos_s = sorted(s for s in s_values if s > 0)
    sym_names = list(c1.keys())

    records = []   # one per (sym, state, bond pair)
    for sym in sym_names:
        for st in state_names:
            c1_abs_max = max(abs(c1[sym][s][b][st])
                             for s in s_values for b in range(N - 1))
            informative = bool(c1_abs_max > C1_INFORMATIVE_FLOOR)
            for (b_lo, b_hi) in pairs:
                D = np.array([c1[sym][s][b_lo][st] - c1[sym][s][b_hi][st]
                              for s in s_values])
                D_by_s = dict(zip(s_values, D))
                oddness = max(abs(D_by_s[+s] + D_by_s[-s]) for s in pos_s)
                d3, d2, d1, d0 = np.polyfit(s_arr, D, 3)
                D_typ = max(abs(D_by_s[+s]) for s in pos_s)
                records.append({
                    "sym": sym, "state": st, "pair": (b_lo, b_hi),
                    "informative": informative,
                    "D_at_s0": float(D_by_s[0.0]),
                    "oddness_residual": float(oddness),
                    "kappa": float(d1), "const": float(d0),
                    "quad": float(d2), "cubic": float(d3),
                    "D_typ_max": float(D_typ),
                })

    inf = [r for r in records if r["informative"]]

    def _mx(rs, key):
        return max((abs(r[key]) for r in rs), default=0.0)

    survival = _mx(records, "D_at_s0")
    oddness = _mx(records, "oddness_residual")
    const_res = _mx(records, "const")
    quad_res = _mx(records, "quad")
    kappa_inf = [abs(r["kappa"]) for r in inf]
    kappa_min_inf = min(kappa_inf, default=0.0)
    kappa_max_inf = max(kappa_inf, default=0.0)
    D_typ_inf = max((r["D_typ_max"] for r in inf), default=0.0)

    gsym_spread = []
    for st in state_names:
        for pr in pairs:
            ks = [r["kappa"] for r in records
                  if r["state"] == st and r["pair"] == pr]
            ks_inf = [r["kappa"] for r in records if r["state"] == st
                      and r["pair"] == pr and r["informative"]]
            if len(ks_inf) >= 2 and max(abs(k) for k in ks_inf) > 1e-4:
                spread = max(ks) - min(ks)
                gsym_spread.append({
                    "state": st, "pair": pr, "kappa_range": float(spread),
                    "kappa_rel_spread": float(spread / max(abs(k) for k in ks_inf)),
                })
    max_rel_spread = max((j["kappa_rel_spread"] for j in gsym_spread), default=0.0)

    print(f"\n  --- N={N} verdict ---", flush=True)
    print(f"  (a) palindromic survival   max|D(s=0)|        = {survival:.2e}"
          f"   (expect ~0: F71 holds for non-uniform palindromic gamma)", flush=True)
    print(f"  (b) oddness                max|D(+s)+D(-s)|   = {oddness:.2e}"
          f"   (typ |D| at max s = {D_typ_inf:.2e})", flush=True)
    print(f"      even powers            max|const|={const_res:.1e}  "
          f"max|quad|={quad_res:.1e}   (expect ~oddness scale)", flush=True)
    print(f"  (c) leading coefficient    informative |kappa| in "
          f"[{kappa_min_inf:.4f}, {kappa_max_inf:.4f}]"
          f"   ({len(inf)} informative records)", flush=True)
    print(f"  (d) gamma_sym-dependence   max relative kappa spread across "
          f"the {len(sym_names)} gamma_sym profiles = {max_rel_spread*100:.1f}%",
          flush=True)

    print(f"\n  kappa_b table (leading coeff of D(b;s) = c_1(b) - c_1(N-2-b)):",
          flush=True)
    for st in state_names:
        print(f"    probe state {st}:", flush=True)
        for pr in pairs:
            row = []
            for sym in sym_names:
                rec = next((r for r in records if r["sym"] == sym
                            and r["state"] == st and r["pair"] == pr), None)
                row.append(f"{sym}={rec['kappa']:+.4f}" if rec else "")
            print(f"      bond pair {pr}:  " + "   ".join(row), flush=True)

    return {
        "criterion_a_palindromic_survival_maxabs": survival,
        "criterion_b_oddness_residual_maxabs": oddness,
        "even_power_const_maxabs": const_res,
        "even_power_quad_maxabs": quad_res,
        "criterion_c_kappa_informative_range": [kappa_min_inf, kappa_max_inf],
        "criterion_c_n_informative_records": len(inf),
        "criterion_d_max_rel_kappa_spread": max_rel_spread,
        "typ_D_informative": D_typ_inf,
        "records": records,
        "gsym_spread": gsym_spread,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=str, default="3,4,5",
                        help="Comma-separated N values (default 3,4,5)")
    args = parser.parse_args()
    N_list = [int(x) for x in args.N.split(",")]

    print("=" * 74, flush=True)
    print("F101 -- numerical verification: D(b) = c_1(b) - c_1(N-2-b) is exactly",
          flush=True)
    print("odd in the F71-anti-palindromic component of per-site gamma.",
          flush=True)
    print("=" * 74, flush=True)
    print(f"  J uniform = {J_UNIFORM_VAL}   dJ probe = {DJ_EXTRACT}", flush=True)
    print(f"  s sweep = {list(S_VALUES)}", flush=True)
    print(f"  N values = {N_list}   (exact (N+1)^2-dim sector restriction)\n",
          flush=True)

    selftest_sector_restriction()

    t_start = time.time()
    all_N = {}
    for N in N_list:
        all_N[str(N)] = run_one_N(N, S_VALUES)
        out_path = RESULTS_DIR / f"f71_nonuniform_gamma_N{N}.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(all_N[str(N)], f, indent=2)
        print(f"  saved {out_path}", flush=True)

    out_path = RESULTS_DIR / f"summary_N{'_'.join(str(n) for n in N_list)}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(all_N, f, indent=2)
    print(f"\nSaved summary: {out_path}", flush=True)
    print(f"Total walltime: {time.time() - t_start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
