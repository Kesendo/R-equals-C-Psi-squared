"""docs/carbon Question 6: do periodic-palindrome deviations correspond to the
V-Effect's boundary-sector breaking?

periodic_palindrome.py establishes that per-element properties (IE_1, Pauling EN,
Allen EN) across a period are palindromic: value_k + value_{N+1-k} ~ const, with
shuffle p-values down to 1e-4. V_EFFECT_BOUNDARY_LOCALIZATION.md establishes that
the framework's palindrome break is confined to boundary sectors 0 < w < N (the
extreme sectors w=0, w=N immune, exact to 1e-15).

Question 6 asks: does the periodic palindrome's deviation pattern match the
V-Effect's boundary localization? The proposed correspondence maps closed shells
(period ends) onto the immune extreme sectors and partly-filled shells (period
middle) onto the breaking boundary sectors.

This script computes the PER-PAIR deviation (each pair sum's distance from the
period mean) that periodic_palindrome.py's aggregate CoV hides, and asks whether
the deviation localizes to the inner (partly-filled) pairs the way the V-Effect
localizes its break to boundary sectors. Investigation only.
"""
import sys

import numpy as np

sys.path.insert(0, "simulations")
import periodic_palindrome as pp

sys.stdout.reconfigure(encoding="utf-8")


def analyse(name, symbols, values):
    """Per-pair deviation of the periodic palindrome. Each pair k joins element k
    with element N-1-k; k=0 is the outermost pair (the period's two ends, the
    closed-shell side), increasing k moves toward the partly-filled middle."""
    n = len(values)
    pairs = pp.pair_sums(values)
    m = len(pairs)
    mean = float(np.mean(pairs))
    rows = []
    for k, ps in enumerate(pairs):
        dev = ps - mean
        rows.append((k, symbols[k], symbols[n - 1 - k], ps, dev, abs(dev) / mean * 100))

    print(f"{name}   (mean pair sum {mean:.3f})")
    for k, sl, sr, ps, dev, devpct in sorted(rows, key=lambda r: -abs(r[4])):
        loc = "outer" if k < m / 2.0 else "inner"
        print(f"  k={k} ({sl:>3s},{sr:>3s})  pair_sum {ps:8.3f}  "
              f"dev {dev:+8.3f} ({devpct:5.1f}%)  [{loc}]")
    inner = [abs(r[4]) for r in rows if r[0] >= m / 2.0]
    outer = [abs(r[4]) for r in rows if r[0] < m / 2.0]
    print(f"  mean |dev|:  outer pairs {np.mean(outer):7.3f}   "
          f"inner pairs {np.mean(inner):7.3f}   "
          f"(V-Effect predicts inner >> outer)")
    print()


def main():
    print("=" * 78)
    print("Q6: per-pair periodic-palindrome deviation vs the V-Effect boundary picture")
    print("=" * 78)
    print("V-Effect: the palindrome break is confined to boundary sectors 0<w<N;")
    print("the extreme sectors w=0, w=N are immune. Mapped onto a period, the")
    print("deviation should sit in the INNER (partly-filled) pairs and leave the")
    print("OUTER (closed-shell-end) pairs near-palindromic.\n")

    print("--- Layer 1: first ionization energies ---\n")
    for name, p in pp.periods.items():
        analyse(name, p["symbols"], p["IE"])

    print("--- Layer 2: Pauling electronegativities ---\n")
    for name, p in pp.en_periods.items():
        analyse(name, p["symbols"], p["EN"])

    print("--- Layer 3: Allen electronegativities ---\n")
    for name, p in pp.allen_en_periods.items():
        analyse(name, p["symbols"], p["EN"])


if __name__ == "__main__":
    main()
