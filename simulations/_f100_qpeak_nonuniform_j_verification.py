#!/usr/bin/env python3
"""_f100_qpeak_nonuniform_j_verification.py

F100 Q_peak half -- numerical witness for non-uniform J.

F100 (typed claim C1QPeakMirrorJParity, proof PROOF_F100_C1_QPEAK_MIRROR_J_PARITY)
states that the per-bond Q_peak deviation

    D(b) := Q_peak(b) - Q_peak(N-2-b)

is an exactly ODD function of the F71-anti-palindromic component of the bond-
coupling profile J, and zero for every palindromic J. The c_1 half of F100 is
witnessed by _f71_nonuniform_j_verification.py; PROOF_F100 carried the Q_peak
half only "by the identical R-conjugation argument". This script supplies the
missing direct numerical witness for the Q_peak half.

THE OBSERVABLE. K_b(Q, t) is F86c's per-bond Q-resonance observable: the
Hellmann-Feynman response of the (n, n+1)-block spatial-sum coherence S(t) to
the bond-b coupling,

    K_b(Q, t) = 2 Re <rho(t)| S_kernel |d rho / d J_b>,
    rho(t)    = exp(L(Q) t) rho_0,    rho_0 = Dicke block probe,
    S(t)      = <rho(t)| S_kernel |rho(t)> = Sum_i 2 |(rho_i(t))_{0,1}|^2.

It is read at the fixed F86a EP time t_peak = 1/(4 gamma_0): there the secular
t*exp(-4 gamma_0 t) Liouvillian response peaks, and K_b(Q, t_peak) is a cleanly
single-peaked function of Q. (A max-over-t reading additionally samples the
off-resonance oscillatory tail and is not single-peaked. F100's parity holds
for K_b at any FIXED t -- R-equivariance is per-t -- so the EP-time reading is
both the clean and a valid choice.)

    Q_peak(b) = argmax_Q |K_b(Q, t_peak)|,  parabola-refined on the Q-grid.

PARAMETERISATION (Q vs the J-profile). For a non-uniform profile Q must stay a
global scalar: fix a profile SHAPE J_hat (uniform = ones), let Q be the global
multiplier, bond couplings J_b(Q) = Q*gamma_0*J_hat[b], and

    L(Q) = D + Sum_b (Q gamma_0 J_hat[b]) M_H_per_bond[b].

This is forced both by the Q_EP apparatus (which needs a scalar tuning knob)
and by F100's fixed-shape J_sym/J_anti argument. The J_sym MAGNITUDE is
degenerate with Q (scaling J_hat by lambda equals scaling Q by 1/lambda), so
unlike the c_1 witness only the profile SHAPE matters here -- uniform vs the
non-uniform palindromic "valley", not three uniform magnitudes.

THE SWEEP mirrors the c_1 witness: J_hat(s) = J_hat_sym + s*J_anti_dir with
J_hat_sym palindromic and J_anti_dir the anti-palindromic linear ramp. F100
predicts:
  (a) D(b; s=0) = 0       -- the F71 Q_peak mirror survives ALL palindromic J,
                             including the non-uniform "valley" shape.
  (b) D(b; -s) = -D(b; +s) -- exactly odd in s, all orders.
  (c) D(b; s) = kappa_b s + O(s^3) -- leading-order linear; even powers vanish.
  (d) kappa_b depends on the J_sym shape.

D(b) is defined on F71 bond PAIRS (b, N-2-b). The central self-paired bond of
an odd bond count maps to itself under F71, so D(central) = 0 identically; it
carries no F100 content and its Q_peak is not extracted.

Both (a) and (b) are algebraically EXACT (R-equivariance gives
K_b(.; J_hat(-s)) = K_{N-2-b}(.; J_hat(s)) at every fixed t), so the measured
survival and oddness residuals are pure floating point -- the witness CONFIRMS
F100, it does not fit it. The discriminating content is (c): D is a real non-
zero signal (kappa_b s) whose even-power coefficients vanish. If F100 were
false the cubic fit's constant and quadratic coefficients would be non-zero at
the signal scale.

GATE. Before the non-uniform sweep, the uniform-shape run must reproduce the
F86c uniform-J mirror D(b) ~ 0 (F100's J_anti=0 corner, which the canonical C#
C2BondLQPeakScan / PerF71OrbitLQPeakTable verify bit-exact). The gate validates
this Python Q_peak engine independently of the C# ResonanceScan.

The block engine (block_L_split_xy + the eigenbasis Duhamel d rho/d J_b) is the
same algorithm as the canonical C# ResonanceScan.ScanAtQ; here it is driven
with a non-uniform L assembly and read at the fixed EP time.

Usage:
  python -u _f100_qpeak_nonuniform_j_verification.py --cases 4:1,5:1,6:1,5:2
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from math import comb
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402

RESULTS_DIR = (Path(__file__).parent / "results"
               / "f100_qpeak_nonuniform_j_verification")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

GAMMA_0 = 0.05                              # F86 standard Z-dephasing rate
T_PEAK = 1.0 / (4.0 * GAMMA_0)              # F86a universal EP time = 1/(4 gamma_0)
S_VALUES = (-0.12, -0.08, -0.04, 0.0, 0.04, 0.08, 0.12)  # J_anti sweep amplitudes
Q_GRID = np.linspace(0.6, 3.4, 141)         # global-multiplier scan, dQ = 0.02
GATE_TOL = 1e-9        # uniform-J mirror is algebraically exact -> FP-equality bound
VERIFY_TOL = 1e-6      # survival/oddness/even-power: F100 makes them exact, residual is FP
SIGNAL_FLOOR = 1e-4    # |D(b)| at max s must exceed this for the test to be non-vacuous
DEGEN_TOL = 1e-10      # eigenvalue near-degeneracy guard in the Duhamel kernel


# ---------------------------------------------------------------------------
# J_sym / J_anti construction (shapes only -- magnitude is degenerate with Q)
# ---------------------------------------------------------------------------
def f71_mirror(profile):
    """F71 chain-mirror on a bond profile: bond b <-> bond N-2-b."""
    return np.asarray(profile)[::-1]


def j_anti_direction(N):
    """Canonical anti-palindromic 'linear ramp' on N-1 bonds:
    J_anti_dir[b] = 2b/(N-2) - 1. Anti-palindromic; the c_1 witness's direction."""
    return np.array([2.0 * b / (N - 2) - 1.0 for b in range(N - 1)])


def j_sym_shapes(N):
    """Palindromic J_sym base SHAPES. 'uniform' = ones; for N >= 4 also the
    non-uniform palindromic 'valley'. No uniform-magnitude sweep: the magnitude
    is degenerate with the global Q for the Q_peak observable, so only the
    shape distinguishes profiles. (At N=3 the valley formula degenerates to
    uniform, hence N >= 4 for the second shape.)"""
    nb = N - 1
    shapes = {"uniform": np.ones(nb)}
    if N >= 4:
        valley = np.array([1.2 - 0.4 * (1.0 - abs(2.0 * b / (nb - 1) - 1.0))
                           for b in range(nb)])
        shapes["valley_nonuniform"] = valley
    return shapes


def central_bond(N):
    """Index of the F71 self-paired bond (b = N-2-b) for an odd bond count,
    else -1. Its Q_peak is not extracted: D(central) = 0 identically."""
    nb = N - 1
    return nb // 2 if nb % 2 == 1 else -1


def assert_decomposition(N):
    """Sanity: anti direction anti-palindromic, J_sym shapes palindromic."""
    anti = j_anti_direction(N)
    assert np.allclose(anti, -f71_mirror(anti)), "J_anti not anti-palindromic"
    for name, sym in j_sym_shapes(N).items():
        assert np.allclose(sym, f71_mirror(sym)), f"shape '{name}' not palindromic"


# ---------------------------------------------------------------------------
# Per-bond Q_peak engine (block-L Duhamel at the fixed EP time; non-uniform L)
# ---------------------------------------------------------------------------
def per_bond_K_curves(N, n, gamma_0, j_shape, Q_grid):
    """Per-bond resonance curve g_b(Q) = |K_b(Q, t_peak)| for one profile shape
    j_shape (length N-1), evaluated at the F86a EP time t_peak = 1/(4 gamma_0).
    Returns an (N-1, len(Q_grid)) array.

    K_b(Q, t) = 2 Re <rho(t)| S_kernel |d rho / d J_b>;
    L(Q) = D + Sum_b (Q gamma_0 j_shape[b]) M_H_per_bond[b];
    the bond derivative kernel is d L / d J_b = M_H_per_bond[b], and
    d rho / d J_b is the eigenbasis Duhamel integral. Same algorithm as the
    canonical C# ResonanceScan.ScanAtQ, driven here with a non-uniform L
    assembly and read at the fixed EP time."""
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    rho0 = fw.dicke_block_probe(N, n)
    S_kernel = fw.spatial_sum_coherence_kernel(N, n)
    nb = N - 1
    j_shape = np.asarray(j_shape, dtype=float)
    t = T_PEAK

    curves = np.zeros((nb, len(Q_grid)))
    for iQ, Q in enumerate(Q_grid):
        J_phys = Q * gamma_0 * j_shape                     # per-bond physical coupling
        L = D + sum(J_phys[b] * M_H_per_bond[b] for b in range(nb))
        evals, R = np.linalg.eig(L)
        R_inv = np.linalg.inv(R)
        c0 = R_inv @ rho0
        e = np.exp(evals * t)
        lam_j = evals[:, None]
        lam_k = evals[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            I_mat = np.where(np.abs(lam_k - lam_j) > DEGEN_TOL,
                             (e[None, :] - e[:, None]) / (lam_k - lam_j),
                             t * e[:, None])
        rho_t = R @ (e * c0)
        for b in range(nb):
            X_b = R_inv @ M_H_per_bond[b] @ R              # bond kernel in eigenbasis
            drho = R @ ((X_b * I_mat) @ c0)                # d rho / d J_b
            K = 2.0 * float(np.real(np.vdot(rho_t, S_kernel @ drho)))
            curves[b, iQ] = abs(K)
    return curves


def parabolic_peak(Q_grid, curve):
    """3-point parabolic refinement of the peak of curve(Q) on a uniform grid.
    Returns the sub-grid Q at the parabola vertex. Raises if the discrete
    maximum sits at a grid edge (peak not bracketed -> widen Q_GRID)."""
    i = int(np.argmax(curve))
    if i == 0 or i == len(curve) - 1:
        raise ValueError(f"Q_peak at grid edge (index {i} of {len(curve)}); "
                         f"widen Q_GRID")
    y0, y1, y2 = curve[i - 1], curve[i], curve[i + 1]
    denom = y0 - 2.0 * y1 + y2
    if abs(denom) < 1e-300:
        return float(Q_grid[i])
    dQ = Q_grid[1] - Q_grid[0]
    return float(Q_grid[i] + 0.5 * dQ * (y0 - y2) / denom)


def qpeak_per_bond(N, n, gamma_0, j_shape, Q_grid):
    """Per-bond Q_peak(b) = argmax_Q |K_b(Q, t_peak)|, parabola-refined.
    The central self-paired bond (odd bond count) is left NaN: it maps to
    itself under F71, so D(central) = 0 identically and it carries no F100
    content. Returns a length-(N-1) array."""
    curves = per_bond_K_curves(N, n, gamma_0, j_shape, Q_grid)
    nb = N - 1
    cb = central_bond(N)
    out = np.full(nb, np.nan)
    for b in range(nb):
        if b != cb:
            out[b] = parabolic_peak(Q_grid, curves[b])
    return out


# ---------------------------------------------------------------------------
# Gate: the uniform-J F71 Q_peak mirror (F100's J_anti = 0 corner)
# ---------------------------------------------------------------------------
def gate_uniform_mirror(c2_N_values, gamma_0):
    """Gate: at uniform J the per-bond Q_peak is F71-mirror-symmetric,
    D(b) = Q_peak(b) - Q_peak(N-2-b) ~ 0 -- the F86c uniform-J mirror that the
    canonical C# C2BondLQPeakScan / PerF71OrbitLQPeakTable verify bit-exact.
    Reproducing it here, in an independent Python engine, validates this
    witness's machinery before the non-uniform sweep. It is also F100's own
    J_anti = 0 corner."""
    print("Gate -- uniform-J F71 Q_peak mirror (F86c, F100's J_anti=0 corner):",
          flush=True)
    worst = 0.0
    for N in c2_N_values:
        qp = qpeak_per_bond(N, 1, gamma_0, np.ones(N - 1), Q_GRID)
        pairs = [(b, N - 2 - b) for b in range(N - 1) if b < N - 2 - b]
        dev = max(abs(qp[lo] - qp[hi]) for (lo, hi) in pairs)
        worst = max(worst, dev)
        print(f"  OK  N={N} c=2  per-bond Q_peak in "
              f"[{np.nanmin(qp):.4f}, {np.nanmax(qp):.4f}]   "
              f"max|D(b)| = {dev:.2e}", flush=True)
    assert worst < GATE_TOL, (f"Gate FAIL: uniform-J Q_peak mirror deviation "
                              f"{worst:.2e} exceeds {GATE_TOL:.0e}")
    print(f"  uniform-J mirror holds (max {worst:.2e} < {GATE_TOL:.0e}); "
          f"endpoint and interior bonds peak in distinct Q-bands "
          f"(F86c per-bond structure).\n", flush=True)


# ---------------------------------------------------------------------------
# Per-case run + analysis
# ---------------------------------------------------------------------------
def run_one_case(N, n, gamma_0, s_values):
    c = fw.chromaticity(N, n)
    nb = N - 1
    assert_decomposition(N)
    anti = j_anti_direction(N)
    shapes = j_sym_shapes(N)
    pairs = [(b, N - 2 - b) for b in range(nb) if b < N - 2 - b]
    block_dim = comb(N, n) * comb(N, n + 1)
    cb = central_bond(N)

    print(f"\n{'=' * 74}", flush=True)
    print(f"N = {N}   n = {n}   chromaticity c = {c}   block dim = {block_dim}"
          f"   bonds 0..{nb - 1}", flush=True)
    print(f"{'=' * 74}", flush=True)
    print(f"  J_sym shapes: {list(shapes)}", flush=True)
    print(f"  J_anti direction: [{', '.join(f'{x:+.3f}' for x in anti)}]",
          flush=True)
    central = f"  (central self-paired bond {cb}, excluded)" if cb >= 0 else ""
    print(f"  F71 bond pairs for D: {pairs}{central}", flush=True)

    qpeak = {name: {} for name in shapes}
    t_run = time.time()
    for shape_name, sym in shapes.items():
        for s in s_values:
            j_hat = sym + s * anti
            t0 = time.time()
            qp = qpeak_per_bond(N, n, gamma_0, j_hat, Q_GRID)
            qpeak[shape_name][s] = qp
            shown = ", ".join("central" if np.isnan(x) else f"{x:.4f}"
                              for x in qp)
            print(f"  shape={shape_name:18s} s={s:+.2f}  Q_peak=[{shown}]  "
                  f"({time.time() - t0:.1f} s)", flush=True)
    print(f"  case done in {time.time() - t_run:.1f} s", flush=True)

    analysis = analyse(N, n, c, qpeak, s_values, pairs, list(shapes))
    return {
        "N": N, "n": n, "chromaticity": c, "gamma_0": gamma_0,
        "block_dim": block_dim, "t_peak": T_PEAK, "s_values": list(s_values),
        "j_anti_direction": anti.tolist(),
        "j_sym_shapes": {k: v.tolist() for k, v in shapes.items()},
        "qpeak": {name: {f"{s:+.2f}": qpeak[name][s].tolist()
                         for s in s_values}
                  for name in shapes},
        "analysis": analysis,
    }


def analyse(N, n, c, qpeak, s_values, pairs, shape_names):
    """Per-case verdict: palindromic survival, oddness, leading coefficient,
    even-power suppression, shape-dependence of kappa."""
    s_arr = np.array(s_values)
    pos_s = sorted(s for s in s_values if s > 0)

    records = []
    for shape in shape_names:
        for (b_lo, b_hi) in pairs:
            D = np.array([qpeak[shape][s][b_lo] - qpeak[shape][s][b_hi]
                          for s in s_values])
            D_by_s = dict(zip(s_values, D))
            oddness = max(abs(D_by_s[+s] + D_by_s[-s]) for s in pos_s)
            d3, d2, d1, d0 = np.polyfit(s_arr, D, 3)
            D_typ = max(abs(D_by_s[+s]) for s in pos_s)
            records.append({
                "shape": shape, "pair": (b_lo, b_hi),
                "D_at_s0": float(D_by_s[0.0]),
                "oddness_residual": float(oddness),
                "kappa": float(d1), "const": float(d0),
                "quad": float(d2), "cubic": float(d3),
                "D_typ_max": float(D_typ),
                "informative": bool(D_typ > SIGNAL_FLOOR),
            })

    inf = [r for r in records if r["informative"]]

    def _mx(rs, key):
        return max((abs(r[key]) for r in rs), default=0.0)

    survival = _mx(records, "D_at_s0")
    oddness = _mx(records, "oddness_residual")
    const_res = _mx(records, "const")
    quad_res = _mx(records, "quad")
    D_typ_inf = max((r["D_typ_max"] for r in inf), default=0.0)
    kappa_inf = [abs(r["kappa"]) for r in inf]

    # shape-dependence of kappa: spread of kappa across J_sym shapes per pair
    shape_spread = []
    for pr in pairs:
        ks = [r["kappa"] for r in records if r["pair"] == pr]
        ks_inf = [r["kappa"] for r in records
                  if r["pair"] == pr and r["informative"]]
        if len(ks_inf) >= 2 and max(abs(k) for k in ks_inf) > SIGNAL_FLOOR:
            spread = max(ks) - min(ks)
            shape_spread.append({
                "pair": pr, "kappa_range": float(spread),
                "kappa_rel_spread": float(spread / max(abs(k) for k in ks_inf)),
            })
    max_rel_spread = max((j["kappa_rel_spread"] for j in shape_spread),
                         default=0.0)

    passed = (survival < VERIFY_TOL and oddness < VERIFY_TOL
              and const_res < VERIFY_TOL and quad_res < VERIFY_TOL
              and D_typ_inf > SIGNAL_FLOOR)

    print(f"\n  --- N={N} n={n} c={c} verdict ---", flush=True)
    print(f"  (a) palindromic survival   max|D(s=0)|       = {survival:.2e}"
          f"   (expect ~0: F71 Q_peak mirror holds for palindromic J)",
          flush=True)
    print(f"  (b) oddness                max|D(+s)+D(-s)|  = {oddness:.2e}",
          flush=True)
    print(f"      even powers            max|const|={const_res:.1e}  "
          f"max|quad|={quad_res:.1e}", flush=True)
    print(f"  (c) signal                 typ |D| at max s  = {D_typ_inf:.2e}"
          f"   ({len(inf)}/{len(records)} informative records)", flush=True)
    if kappa_inf:
        print(f"      leading coeff          |kappa| in "
              f"[{min(kappa_inf):.4f}, {max(kappa_inf):.4f}]", flush=True)
    print(f"  (d) shape-dependence       max relative kappa spread = "
          f"{max_rel_spread * 100:.1f}%", flush=True)
    print(f"  VERDICT: {'PASS' if passed else 'FAIL'}"
          f"   (residuals < {VERIFY_TOL:.0e}, signal > {SIGNAL_FLOOR:.0e})",
          flush=True)

    print(f"\n  kappa_b table (leading coeff of "
          f"D(b;s) = Q_peak(b) - Q_peak(N-2-b)):", flush=True)
    for pr in pairs:
        row = []
        for shape in shape_names:
            rec = next((r for r in records if r["shape"] == shape
                        and r["pair"] == pr), None)
            if rec:
                row.append(f"{shape}={rec['kappa']:+.4f}")
        print(f"    bond pair {pr}:  " + "   ".join(row), flush=True)

    return {
        "passed": bool(passed),
        "criterion_a_palindromic_survival_maxabs": float(survival),
        "criterion_b_oddness_residual_maxabs": float(oddness),
        "even_power_const_maxabs": float(const_res),
        "even_power_quad_maxabs": float(quad_res),
        "criterion_c_typ_D_informative": float(D_typ_inf),
        "n_informative_records": len(inf),
        "criterion_d_max_rel_kappa_spread": float(max_rel_spread),
        "records": records,
        "shape_spread": shape_spread,
    }


# ---------------------------------------------------------------------------
def parse_cases(spec):
    """Parse 'N:n,N:n,...' into a list of (N, n) tuples."""
    cases = []
    for tok in spec.split(","):
        N_str, n_str = tok.split(":")
        cases.append((int(N_str), int(n_str)))
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=str, default="4:1,5:1,6:1,5:2",
                        help="Comma-separated N:n cases (default 4:1,5:1,6:1,5:2; "
                             "n=1 is c=2 across the F86 Q_peak N range; 5:2 is a "
                             "block-agnostic c=3 spot check)")
    args = parser.parse_args()
    cases = parse_cases(args.cases)

    print("=" * 74, flush=True)
    print("F100 Q_peak half -- numerical witness for non-uniform J.", flush=True)
    print("D(b) = Q_peak(b) - Q_peak(N-2-b) is exactly odd in J_anti,", flush=True)
    print("zero for every palindromic J.  Observable: F86c per-bond K_b(Q,t).",
          flush=True)
    print("=" * 74, flush=True)
    print(f"  gamma_0 = {GAMMA_0}   t_peak = 1/(4 gamma_0) = {T_PEAK:.2f}",
          flush=True)
    print(f"  Q grid: [{Q_GRID[0]:.2f}, {Q_GRID[-1]:.2f}] x {len(Q_GRID)}  "
          f"(dQ = {Q_GRID[1] - Q_GRID[0]:.3f})", flush=True)
    print(f"  s sweep: {list(S_VALUES)}", flush=True)
    print(f"  cases (N, n): {cases}\n", flush=True)

    c2_N = sorted({N for (N, n) in cases if n == 1})
    if c2_N:
        gate_uniform_mirror(c2_N, GAMMA_0)

    t_start = time.time()
    all_cases = {}
    for (N, n) in cases:
        res = run_one_case(N, n, GAMMA_0, S_VALUES)
        key = f"N{N}_n{n}"
        all_cases[key] = res
        out = RESULTS_DIR / f"f100_qpeak_{key}.json"
        with open(out, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2)
        print(f"  saved {out}", flush=True)

    summary = RESULTS_DIR / "summary.json"
    with open(summary, "w", encoding="utf-8") as f:
        json.dump(all_cases, f, indent=2)

    all_pass = all(r["analysis"]["passed"] for r in all_cases.values())
    print(f"\n{'=' * 74}", flush=True)
    print(f"OVERALL: {'ALL CASES PASS' if all_pass else 'FAILURE'}  -- "
          f"F100 Q_peak half {'verified' if all_pass else 'NOT verified'} "
          f"for non-uniform J", flush=True)
    for key, r in all_cases.items():
        a = r["analysis"]
        print(f"  {key}: survival {a['criterion_a_palindromic_survival_maxabs']:.1e}"
              f"   oddness {a['criterion_b_oddness_residual_maxabs']:.1e}"
              f"   signal {a['criterion_c_typ_D_informative']:.1e}"
              f"   -> {'PASS' if a['passed'] else 'FAIL'}", flush=True)
    print(f"  saved summary: {summary}", flush=True)
    print(f"  total walltime: {time.time() - t_start:.1f} s", flush=True)
    print("=" * 74, flush=True)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
