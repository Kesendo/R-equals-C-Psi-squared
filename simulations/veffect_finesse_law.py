#!/usr/bin/env python3
"""The cavity finesse law of experiments/VEFFECT_CAVITY_MODES.md, checked against its own dump.

Result 4 of that document used to hedge the highest-Q mode as "Q_max ≈ 4J/γ_grid converging toward
75 as N grows". It is exact, and the limit is 80, not 75:

    Q_max(N) = 2J (1 + cos(π/N)) / γ            J = 1, γ = 0.05

This script checks it against every N the committed sweep carries,
simulations/results/veffect_cavity_modes.txt, which is where the measured Q_max values live. No
new physics is computed here: the dump is the measurement, the law is the claim, and this is the
comparison. What the law does NOT explain is why cos(π/N), the RING dispersion, fits an OPEN chain
(N-1 bonds) rather than the open chain's own cos(π/(N+1)); the document says so too.
"""
from __future__ import annotations
import os
import re
import sys
import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

DUMP = os.path.join("simulations", "results", "veffect_cavity_modes.txt")
J, GAMMA = 1.0, 0.05


def law(n, J=J, gamma=GAMMA):
    return 2 * J * (1 + np.cos(np.pi / n)) / gamma


def measured(path=DUMP):
    """Q_max of the weight-1 shell per N, read off the 'Highest-Q mode' lines of the dump."""
    out = {}
    current = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            m = re.match(r"\s*N=(\d+):", line)
            if m:
                current = int(m.group(1))
            m = re.search(r"Highest-Q mode: .*?Q=([\d.]+), weight shell k=1", line)
            if m and current is not None:
                out[current] = float(m.group(1))
    return out


if __name__ == "__main__":
    if not os.path.exists(DUMP):
        sys.exit(f"missing {DUMP}; run the veffect cavity sweep first")
    meas = measured()
    if not meas:
        sys.exit(f"no 'Highest-Q mode ... k=1' lines found in {DUMP}")
    print(f"Q_max(N) = 2J(1 + cos(pi/N))/gamma   with J={J}, gamma={GAMMA}")
    print(f"{'N':>3} {'measured':>10} {'law':>12} {'rel dev':>10}")
    ok = True
    for n in sorted(meas):
        q = law(n)
        dev = abs(meas[n] - q) / q
        ok &= dev <= 1e-3
        print(f"{n:>3} {meas[n]:>10.1f} {q:>12.4f} {dev:>10.1e}")
    print(f"\nfinesse limit 4J/gamma = {4*J/GAMMA:.0f}, approached from below and reached at no N")
    print("all N within 0.1 percent" if ok else "MISMATCH")
    sys.exit(0 if ok else 1)
