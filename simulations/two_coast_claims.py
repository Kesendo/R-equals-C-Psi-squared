#!/usr/bin/env python3
"""Every number experiments/THE_PALINDROME_CLASSIFIER.md quotes that the sweep dump does not carry.

simulations/results/two_coast_sweep.tsv holds the two coasts on their 190-angle grid, so the
numbers read off that grid are traceable by opening it. The document also quotes numbers that live
BELOW the grid, BESIDE it (other Hamiltonians), or ACROSS other gamma values, and those had no
committed verifier. This script is that verifier: each check prints PASS or FAIL against the value
as the document states it, and the exit code is 0 only if all of them pass.

Setting throughout, unless a check says otherwise: N=4, gamma=0.05 (sigma=0.2), Z-dephasing on
every site, the two Hamiltonian families of two_coast_sweep.py, optimal (bottleneck) pairing.

Runtime a few minutes; the gamma-scaling section dominates it.
"""
from __future__ import annotations
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from two_coast_sweep import (H_field, H_frustration, liouvillian, pairing_true, part_error,
                             placed, SIGMA, GAMMA, N)

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

RESULTS: list[tuple[bool, str]] = []


def check(ok: bool, label: str, got, quoted) -> None:
    RESULTS.append((ok, label))
    print(f"  {'PASS' if ok else 'FAIL'}  {label}: got {got}, document says {quoted}")


def close(a, b, rel=0.02) -> bool:
    return abs(a - b) <= rel * max(abs(b), 1e-300)


def pair(H, gamma=GAMMA, n=N):
    return pairing_true(np.linalg.eigvals(liouvillian(H, gamma=gamma, n=n)), n * gamma)


def rate(H, gamma=GAMMA, n=N):
    v = np.linalg.eigvals(liouvillian(H, gamma=gamma, n=n))
    return part_error(v.real, -v.real - 2 * n * gamma)


# ---------------------------------------------------------------- 1. the two quadratic laws
def section_laws():
    print("\n[1] the quadratic laws and how far each one reaches")
    cf_frust = pair(H_frustration(1e-5)) / 1e-10
    cf_field = 2 * SIGMA * (np.pi / 180) ** 2
    check(close(cf_frust, 2.79, 0.01), "frustration prefactor per deg^2", f"{cf_frust:.3f}", "2.79")
    check(close(cf_field, 1.22e-4, 0.01), "field prefactor per deg^2 (= 2 sigma in rad)",
          f"{cf_field:.4e}", "1.22e-4")

    r26 = pair(H_frustration(2.6e-4)) / (cf_frust * (2.6e-4) ** 2)
    r1e3 = pair(H_frustration(1e-3)) / (cf_frust * (1e-3) ** 2)
    check(abs(r26 - 1) <= 0.015, "frustration law within a percent at 2.6e-4 deg",
          f"{100*(r26-1):+.1f}%", "good to a percent")
    check(1.15 <= r1e3 <= 1.25, "frustration law at 1e-3 deg", f"{100*(r1e3-1):+.0f}%",
          "twenty percent low")
    f10 = pair(H_field(10.0)) / (cf_field * 100)
    f32 = pair(H_field(32.0)) / (cf_field * 32 ** 2)
    check(abs(f10 - 1) <= 0.011, "field law at 10 deg", f"{100*(f10-1):+.1f}%", "a percent")
    check(abs(f32 - 1) <= 0.11, "field law at 32 deg", f"{100*(f32-1):+.0f}%", "ten percent")


# ---------------------------------------------------------------- 2. handover and shoulder
def section_shoulder():
    print("\n[2] the handover is not monotone, and the shoulder beyond it")
    hi = pair(H_frustration(2.1e-3))
    lo = pair(H_frustration(2.7e-3))
    check(close(hi, 3.0e-5, 0.03), "local maximum at 2.1e-3 deg", f"{hi:.4e}", "3.0e-5")
    check(close(lo, 1.9e-5, 0.03), "local minimum at 2.7e-3 deg", f"{lo:.4e}", "1.9e-5")
    check(lo < hi * 0.7, "the fall between them", f"{100*(1-lo/hi):.0f}%", "by a third")
    s1 = pair(H_frustration(0.006))
    s2 = pair(H_frustration(0.06))
    check(close(s1, 3e-5, 0.15), "shoulder at 0.006 deg", f"{s1:.3e}", "3e-5")
    check(close(s2, 8e-5, 0.15), "shoulder at 0.06 deg", f"{s2:.3e}", "8e-5")


# ---------------------------------------------------------------- 3. the mechanism
def bonds(kind, n=N):
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    if kind == "none":
        return H
    for b in range(n - 1):
        for letter in {"ising": "Z", "xy": "XY", "heis": "XYZ"}[kind]:
            H += placed(letter + letter, b, n)
    return H


def zfield(profile, n=N):
    return sum(h * placed("Z", l, n) for l, h in enumerate(profile))


def section_mechanism():
    print("\n[3] a commuting longitudinal field is inert on the rates, which is not protection")
    check(rate(bonds("none") + zfield([1, 2, 3, 4])) < 1e-15,
          "no bonds, non-uniform profile", f"{rate(bonds('none') + zfield([1,2,3,4])):.1e}", "3e-17")
    check(rate(bonds("ising") + zfield([1, 2, 3, 4])) < 1e-15,
          "Ising ZZ, non-uniform profile", f"{rate(bonds('ising') + zfield([1,2,3,4])):.1e}",
          "every profile keeps the rates")
    for kind, broken in (("xy", 8.9e-2), ("heis", 1.0e-1)):
        u = rate(bonds(kind) + zfield([1, 1, 1, 1]))
        b = rate(bonds(kind) + zfield([1, 2, 3, 4]))
        check(u < 1e-13, f"{kind} uniform field", f"{u:.1e}", "1e-14")
        check(close(b, broken, 0.02), f"{kind} profile (1,2,3,4)", f"{b:.4e}", f"{broken}")

    # the counterexample: commuting, inert, and still off the mirror
    t = np.deg2rad(38.0)
    Hr = (np.cos(t) * placed("XIX", 0) + np.sin(t) * placed("XXX", 0)
          + placed("XXY", 0) + placed("YXX", 0))
    breaks = []
    for h in (0.0, 0.5, 1.0, 5.0):
        F = h * placed("Z", 3)
        assert np.linalg.norm(Hr @ F - F @ Hr) == 0.0, "the counterexample must commute exactly"
        breaks.append(rate(Hr + F))
    check(all(close(x, 1.170227e-2, 1e-4) for x in breaks),
          "frustration(38 deg) on sites 0-2 plus h*Z on site 3, h = 0, 0.5, 1, 5",
          f"{breaks[0]:.6e} at every h", "1.17e-2 at h = 0, 0.5, 1 and 5 alike")

    for eps, quoted in ((0.01, 1.8e-4), (1e-4, 1.3e-6)):
        r = rate(bonds("heis") + zfield([1, 1, 1, 1 + eps]))
        check(close(r, quoted, 0.03), f"Heisenberg with one site at 1+{eps}", f"{r:.3e}", f"{quoted}")


# ---------------------------------------------------------------- 4. the aggregate census
def section_census():
    print("\n[4] what the worst-case meter hides")
    for label, H, exact, quoted in (("field 90 deg", H_field(90.0), 230, "230 of 256 exact"),
                                    ("frustration 38 deg", H_frustration(38.0), 0, "none exact")):
        v = np.linalg.eigvals(liouvillian(H))
        d = np.abs(v[None, :] - (-v - 2 * SIGMA)[:, None]).min(axis=1)
        n_exact = int((d < 1e-9).sum())
        check(n_exact == exact, f"{label}: eigenvalues with an exact partner", n_exact, quoted)
        if exact == 230:
            at20 = int((np.abs(d - 0.20) < 1e-6).sum())
            at40 = int((np.abs(d - 0.40) < 1e-6).sum())
            check((at20, at40) == (24, 2), "field 90 deg: how the other 26 miss",
                  f"{at20} by 0.20, {at40} by 0.40", "twenty-four by 0.20 and two by 0.40")
        else:
            check(close(float(np.median(d)), 4.06e-4, 0.02),
                  "frustration 38 deg: median miss", f"{np.median(d):.3e}", "4e-4")


# ---------------------------------------------------------------- 5. ceiling and crossing
def section_ceiling():
    print("\n[5] the 2 sigma ceiling and the tolerance crossing")
    rng = np.random.default_rng(11)
    worst = 0.0
    for _ in range(8):
        A = rng.normal(size=(2 ** N, 2 ** N)) + 1j * rng.normal(size=(2 ** N, 2 ** N))
        worst = max(worst, pair(20 * (A + A.conj().T) / 2))
    check(worst <= 2 * SIGMA + 1e-12, "8 random Hermitian H stay under the ceiling",
          f"max {worst:.4f}", f"never above 2 sigma = {2*SIGMA}")
    check(close(pair(H_field(90.0)), 0.4, 1e-6), "field 90 deg saturates it",
          f"{pair(H_field(90.0)):.6f}", "0.40 = 2 sigma")

    lo, hi = 0.05, 0.15
    for _ in range(30):
        mid = (lo + hi) / 2
        if pair(H_field(mid)) < 1e-6: lo = mid
        else: hi = mid
    check(close(lo, 0.0906, 0.01), "field crosses the 1e-6 tolerance", f"{lo:.4f} deg", "about 0.09 deg")


# ---------------------------------------------------------------- 6. the budget across gamma
def section_budget():
    print("\n[6] the budget ratios, and that they are gamma-specific")
    quoted = {0.2: (11, 9.7), 0.05: (12, 151), 0.025: (25, 605)}
    for g in (0.2, 0.05, 0.025):
        f1 = pair(H_field(1.0), gamma=g)
        cf = pair(H_frustration(1e-5), gamma=g) / 1e-10
        asym = np.sqrt(cf / (2 * 4 * g * (np.pi / 180) ** 2))
        lo, hi = 0.001, 0.5
        for _ in range(30):
            mid = (lo + hi) / 2
            if pair(H_frustration(mid), gamma=g) < f1: lo = mid
            else: hi = mid
        dev, qa = 1.0 / lo, quoted[g]
        check(close(dev, qa[0], 0.06), f"gamma={g}: device-scale ratio", f"{dev:.1f}", f"{qa[0]}")
        check(close(asym, qa[1], 0.03), f"gamma={g}: asymptotic ratio", f"{asym:.4g}", f"{qa[1]}")
        if g == GAMMA:
            check(close(lo, 0.08, 0.05), "gamma=0.05: equal-cost frustration angle",
                  f"{lo:.4f} deg", "about 0.08 deg")


if __name__ == "__main__":
    print(f"N={N} gamma={GAMMA} sigma={SIGMA}; checks against experiments/THE_PALINDROME_CLASSIFIER.md")
    section_laws()
    section_shoulder()
    section_mechanism()
    section_census()
    section_ceiling()
    section_budget()
    ok = sum(1 for r, _ in RESULTS if r)
    print(f"\n{ok}/{len(RESULTS)} checks pass")
    for r, label in RESULTS:
        if not r:
            print(f"  FAILED: {label}")
    sys.exit(0 if ok == len(RESULTS) else 1)
