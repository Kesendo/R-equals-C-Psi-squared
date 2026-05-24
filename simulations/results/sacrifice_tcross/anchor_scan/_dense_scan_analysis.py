"""Dense gamma_0 scan analysis: 24-point a(theta_0) curve.

Combines original 12 anchor points with 12 new dense weak-gamma points.
Fits t_cross(r) = a - b/r for r=1..N (matching existing _combined_analysis.json
convention; spec wording 'r=1..4' interpreted as r=1..N=5 to match prior fits).
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent

# All 24 points: gamma_0, theta_0 (deg), role.
# theta_0 = atan(1 / gamma_0) in degrees (since J=1, polarity-mirror parameterization
# from prior anchor-scan setup; preserved from original JSON).
POINTS = [
    # Original 12 (from _combined_analysis.json, sorted by theta)
    (3.0,    18.435,  "off-between(0,30)",       "g0_3.0.csv"),
    (2.246,  24.000,  "MIDPOINT(18.4,30)",       "g0_24deg.csv"),
    (1.7321, 30.000,  "ANCHOR-30",                "g0_30deg.csv"),
    (1.3032, 37.500,  "MIDPOINT(30,45)",         "g0_37.5deg.csv"),
    (1.0,    45.000,  "ANCHOR-45",                "g0_45deg.csv"),
    (0.7673, 52.500,  "MIDPOINT(45,60)",         "g0_52.5deg.csv"),
    (0.5774, 60.000,  "ANCHOR-60",                "g0_60deg.csv"),
    (0.5,    63.435,  "midpoint(60,76)*",         "g0_0.50.csv"),
    (0.45,   65.770,  "DENSE-65.77",              "g0_0.45.csv"),  # NEW
    (0.404,  68.000,  "MIDPOINT(60,76)/DIP",     "g0_68deg.csv"),
    (0.38,   69.180,  "DENSE-69.18",              "g0_0.38.csv"),  # NEW
    (0.37,   69.738,  "PROBE-grazing(below)",     "g0_0.37.csv"),  # PROBE
    (0.36,   70.198,  "PROBE-grazing(above)",     "g0_0.36.csv"),  # PROBE
    (0.35,   70.710,  "DENSE-70.71",              "g0_0.35.csv"),  # NEW
    (0.32,   72.260,  "DENSE-72.26",              "g0_0.32.csv"),  # NEW
    (0.30,   73.300,  "DENSE-73.30",              "g0_0.30.csv"),  # NEW
    (0.25,   75.964,  "off-between(76,87)/PEAK", "g0_0.25.csv"),
    (0.22,   77.590,  "DENSE-77.59",              "g0_0.22.csv"),  # NEW
    (0.20,   78.690,  "DENSE-78.69",              "g0_0.20.csv"),  # NEW
    (0.18,   79.790,  "DENSE-79.79",              "g0_0.18.csv"),  # NEW
    (0.1584, 81.000,  "MIDPOINT(76,87)/DIP",     "g0_81deg.csv"),
    (0.14,   82.010,  "DENSE-82.01",              "g0_0.14.csv"),  # NEW
    (0.12,   83.150,  "DENSE-83.15",              "g0_0.12.csv"),  # NEW
    (0.10,   84.290,  "DENSE-84.29",              "g0_0.10.csv"),  # NEW
    (0.08,   85.430,  "DENSE-85.43",              "g0_0.08.csv"),  # NEW
    (0.05,   87.138,  "baseline-uniform",         "g0_0.05.csv"),
]


def parse_tcross(csv_path: Path) -> list[float]:
    """Read the '# t_cross_i' comment line from a sacrifice-tcross CSV."""
    with open(csv_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("# t_cross_i"):
                # "# t_cross_i (threshold=0.2500): 0.97 1.38 1.81 3.10 3.11"
                payload = line.split(":", 1)[1].strip()
                tokens = payload.split()
                return [float(t) if t != "null" else math.nan for t in tokens]
    raise ValueError(f"No t_cross_i line found in {csv_path}")


def fit_ab(tc: list[float]) -> tuple[float, float, float]:
    """Fit t_cross(r) = a - b/r for r=1..N. Returns (a, b, R^2)."""
    n = len(tc)
    r = np.arange(1, n + 1, dtype=float)
    y = np.array(tc, dtype=float)
    mask = ~np.isnan(y)
    if mask.sum() < 2:
        return float("nan"), float("nan"), float("nan")
    x = 1.0 / r[mask]
    yv = y[mask]
    # Linear regression: y = a - b*x  ->  y = a + (-b)*x
    A = np.vstack([np.ones_like(x), -x]).T
    sol, *_ = np.linalg.lstsq(A, yv, rcond=None)
    a, b = float(sol[0]), float(sol[1])
    ypred = a - b * x
    ss_res = float(np.sum((yv - ypred) ** 2))
    ss_tot = float(np.sum((yv - yv.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return a, b, r2


def linear_envelope_residual(theta: float, a: float, points: list[tuple]) -> tuple[float, float]:
    """Compute distance below linear envelope through neighbors >3deg from theta."""
    far_left = [(t, v) for (t, v) in points if t < theta - 3.0]
    far_right = [(t, v) for (t, v) in points if t > theta + 3.0]
    if not far_left or not far_right:
        return float("nan"), float("nan")
    # Take nearest neighbors beyond 3deg radius
    tL, vL = max(far_left, key=lambda x: x[0])
    tR, vR = min(far_right, key=lambda x: x[0])
    # Linear interp at theta
    env = vL + (vR - vL) * (theta - tL) / (tR - tL)
    depth = env - a
    return env, depth


def main():
    results = []
    for g0, theta, role, fname in POINTS:
        path = ROOT / fname
        if not path.exists():
            print(f"MISSING: {fname}")
            continue
        tc = parse_tcross(path)
        a, b, r2 = fit_ab(tc)
        results.append({
            "fname": fname,
            "gamma0": g0,
            "theta": theta,
            "role": role,
            "tc": tc,
            "a": a,
            "b": b,
            "r2": r2,
        })

    # Sort by theta
    results.sort(key=lambda r: r["theta"])

    # Write JSON
    with open(ROOT / "_dense_scan_analysis.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    # Write tabular output
    print("=" * 130)
    print("24-POINT a(theta_0) TABLE (sorted by theta_0)")
    print("=" * 130)
    hdr = f"{'gamma_0':>8} {'theta':>7} {'role':<28} {'tc0':>7} {'tc1':>7} {'tc2':>7} {'tc3':>7} {'tc4':>7} {'a':>8} {'b':>8} {'R^2':>6}"
    print(hdr)
    print("-" * 130)
    for r in results:
        tc = r["tc"]
        tc_str = " ".join(f"{v:>7.4f}" if not math.isnan(v) else f"{'null':>7}" for v in tc)
        r2_str = f"{r['r2']:.4f}" if not math.isnan(r['r2']) else "  nan "
        print(f"{r['gamma0']:>8.4f} {r['theta']:>7.3f} {r['role']:<28} {tc_str} {r['a']:>8.4f} {r['b']:>8.4f} {r2_str:>6}")
    print("=" * 130)

    # Test 1: Dip characterization
    # Focus on the two known dips: theta ~ 68 deg and theta ~ 81 deg
    print("\n" + "=" * 100)
    print("TEST 1: DIP CHARACTERIZATION")
    print("=" * 100)
    ta_pairs = [(r["theta"], r["a"]) for r in results]

    for dip_label, target_theta, lo, hi in [("DIP-1 (~68deg)", 68.0, 60.0, 75.0),
                                              ("DIP-2 (~81deg)", 81.0, 76.0, 85.0)]:
        print(f"\n--- {dip_label} ---")
        local = [(t, v) for t, v in ta_pairs if lo <= t <= hi]
        local.sort()
        for t, v in local:
            mark = " <-- check" if abs(t - target_theta) < 1.5 else ""
            print(f"  theta={t:.3f}  a={v:.4f}{mark}")
        # Find minimum in this band
        min_t, min_a = min(local, key=lambda x: x[1])
        # Local envelope
        env, depth = linear_envelope_residual(min_t, min_a, ta_pairs)
        print(f"  Local minimum: theta={min_t:.3f}, a={min_a:.4f}")
        print(f"  Linear envelope at theta={min_t:.3f}: {env:.4f}  ->  depth = {depth:.4f}")
        # Parabolic interpolation if we have 3 surrounding points
        idx = next((i for i, (t, _) in enumerate(local) if t == min_t), None)
        if idx is not None and 0 < idx < len(local) - 1:
            t_m, a_m = local[idx - 1]
            t_c, a_c = local[idx]
            t_p, a_p = local[idx + 1]
            # Quadratic vertex
            denom = (t_m - t_c) * (t_m - t_p) * (t_c - t_p)
            if abs(denom) > 1e-12:
                A = (t_p * (a_c - a_m) + t_c * (a_m - a_p) + t_m * (a_p - a_c)) / denom
                B = (t_p * t_p * (a_m - a_c) + t_c * t_c * (a_p - a_m) + t_m * t_m * (a_c - a_p)) / denom
                if abs(A) > 1e-12:
                    theta_vertex = -B / (2 * A)
                    C = a_m - A * t_m * t_m - B * t_m
                    a_vertex = A * theta_vertex * theta_vertex + B * theta_vertex + C
                    print(f"  Parabolic minimum: theta_min = {theta_vertex:.3f}, a_min = {a_vertex:.4f}")
                    # Width at half-depth: a_half = a_vertex + depth/2
                    if not math.isnan(env):
                        a_half = a_vertex + depth / 2.0
                        # Solve A*t^2 + B*t + C = a_half: A*t^2 + B*t + (C - a_half) = 0
                        disc = B * B - 4 * A * (C - a_half)
                        if disc > 0 and A != 0:
                            tL = (-B - math.sqrt(disc)) / (2 * A)
                            tR = (-B + math.sqrt(disc)) / (2 * A)
                            width = abs(tR - tL)
                            print(f"  Width at half-depth: {width:.3f} deg (parabolic)")
                # Asymmetry: slope difference at vertex
                # Left slope at midpoint between left and center:
                slope_left = (a_c - a_m) / (t_c - t_m)
                slope_right = (a_p - a_c) / (t_p - t_c)
                print(f"  Slope left:  {slope_left:+.4f} a/deg  (theta {t_m:.2f}->{t_c:.2f})")
                print(f"  Slope right: {slope_right:+.4f} a/deg  (theta {t_c:.2f}->{t_p:.2f})")
                print(f"  Asymmetry (|right| - |left|): {abs(slope_right) - abs(slope_left):+.4f}")

    # Test 2: Discovery of additional features
    print("\n" + "=" * 100)
    print("TEST 2: ADDITIONAL FEATURES")
    print("=" * 100)
    # Scan for additional local minima in [60, 87]
    print("\nLooking for local minima in theta in [60, 87] (a smaller than both neighbors):")
    relevant = [(r["theta"], r["a"]) for r in results if 60 <= r["theta"] <= 87]
    relevant.sort()
    for i in range(1, len(relevant) - 1):
        t_prev, a_prev = relevant[i - 1]
        t_cur, a_cur = relevant[i]
        t_next, a_next = relevant[i + 1]
        if a_cur < a_prev and a_cur < a_next:
            print(f"  LOCAL MIN: theta={t_cur:.3f}, a={a_cur:.4f}  (neighbors: a={a_prev:.4f}, a={a_next:.4f})")
        elif a_cur > a_prev and a_cur > a_next:
            print(f"  LOCAL MAX: theta={t_cur:.3f}, a={a_cur:.4f}  (neighbors: a={a_prev:.4f}, a={a_next:.4f})")

    # Check periodic spacing: dips at 68 and 81 -> delta_theta=13. Predict 55 and 94.
    print("\nPeriodic-spacing test:")
    print("  Known dips at ~68deg and ~81deg, delta=13deg.")
    print("  Predicted additional dip at ~55deg (=68-13): check 52.5deg point.")
    print(f"     theta=52.500: a={[r['a'] for r in results if abs(r['theta']-52.5)<0.01][0]:.4f}")
    print(f"     theta=60.000: a={[r['a'] for r in results if abs(r['theta']-60.0)<0.01][0]:.4f}")
    print(f"     theta=63.435: a={[r['a'] for r in results if abs(r['theta']-63.435)<0.01][0]:.4f}")
    print(f"     theta=65.770: a={[r['a'] for r in results if abs(r['theta']-65.770)<0.01][0]:.4f}")
    print("  Predicted additional dip at ~94deg: beyond domain (max theta=87.138).")

    # Test 3: Near-baseline behavior (theta > 82deg)
    print("\n" + "=" * 100)
    print("TEST 3: NEAR-BASELINE BEHAVIOR (theta > 82deg)")
    print("=" * 100)
    near_baseline = [(r["theta"], r["a"], r["gamma0"]) for r in results if r["theta"] > 81.5]
    near_baseline.sort()
    baseline_a = [r["a"] for r in results if r["theta"] > 87][0]
    print(f"  Baseline (gamma_0=0.05, theta=87.138deg): a = {baseline_a:.4f}")
    print(f"\n  Trajectory of a as theta -> 87.138deg:")
    for t, a, g in near_baseline:
        delta = a - baseline_a
        print(f"    gamma_0={g:.4f}  theta={t:.3f}  a={a:.4f}  delta_from_baseline={delta:+.4f}")

    # Test 4: Anchor sanity
    print("\n" + "=" * 100)
    print("TEST 4: ANCHOR SANITY (smoothness of 60deg in denser scan)")
    print("=" * 100)
    # 60 deg: gamma_0=0.5774
    # Linear interp from neighbors (now denser): closest below and above
    sixty = [r for r in results if abs(r["theta"] - 60) < 0.01][0]
    below = max((r for r in results if r["theta"] < 60), key=lambda r: r["theta"])
    above = min((r for r in results if r["theta"] > 60), key=lambda r: r["theta"])
    a_interp_60 = below["a"] + (above["a"] - below["a"]) * (60 - below["theta"]) / (above["theta"] - below["theta"])
    delta_60 = sixty["a"] - a_interp_60
    print(f"\n  60-deg anchor (gamma_0=0.5774, a={sixty['a']:.4f})")
    print(f"  Linear interp from neighbors theta={below['theta']:.2f} (a={below['a']:.4f}) and theta={above['theta']:.2f} (a={above['a']:.4f}):")
    print(f"    interp(60) = {a_interp_60:.4f}, residual = {delta_60:+.4f}")

    # Compute the same residual for every interior point (excluding endpoints)
    print(f"\n  Smoothness residual for each interior point (a - linear_interp(neighbors_in_table)):")
    print(f"  {'theta':>8} {'role':<30} {'a':>8} {'interp':>8} {'residual':>10}")
    sorted_pts = sorted(results, key=lambda r: r["theta"])
    for i in range(1, len(sorted_pts) - 1):
        cur = sorted_pts[i]
        prev = sorted_pts[i - 1]
        nxt = sorted_pts[i + 1]
        a_interp = prev["a"] + (nxt["a"] - prev["a"]) * (cur["theta"] - prev["theta"]) / (nxt["theta"] - prev["theta"])
        delta = cur["a"] - a_interp
        flag = "  *<-- 60deg anchor" if abs(cur["theta"] - 60) < 0.01 else ""
        flag += "  *F99-45" if abs(cur["theta"] - 45) < 0.01 else ""
        flag += "  *F99-30" if abs(cur["theta"] - 30) < 0.01 else ""
        print(f"  {cur['theta']:>8.3f} {cur['role']:<30} {cur['a']:>8.4f} {a_interp:>8.4f} {delta:>+10.4f}{flag}")


if __name__ == "__main__":
    main()
