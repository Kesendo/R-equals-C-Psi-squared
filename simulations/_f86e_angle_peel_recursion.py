"""F86e: Tom's integer-peel angle recursion applied to σ_0(∞).

Tom 2026-05-20: "subtract the integer each time to get to the next
angle; the principle inherits, first −2, now minus the integer."

Recursion:  x  →  arccos(x − floor(x)) · 180/π   [degrees]
Each step extracts an integer "anchor" (the floor) and carries the
fractional part into the next arccos.

RETRACTED (2026-05-20): this test cannot work. arccos∘frac is a mixing
map; iterating it gives an aperiodic anchor sequence for ANY input of
generic precision, structured or not. "Aperiodic" is therefore the
guaranteed outcome and says nothing about whether σ_0 has hidden
structure; the controls (2√2, π, √3, e, a plain decimal) all peel
identically for exactly this reason. The script is kept as a record of
a dead end. σ_0 is properly characterised as a Schur-multiplier norm
(SigmaZeroCommutatorNormClaim, Tier1Derived), not an angle peel.

NOTE: arccos∘frac also amplifies error ~×3-5 per step, so beyond ~12
steps the digits are precision-limited.
"""
from __future__ import annotations

import math
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def peel(x0: float, steps: int = 16):
    """Iterate x → arccos(frac(x))·180/π, return (anchors, remainders, angles)."""
    anchors, rems, angles = [], [], []
    x = x0
    for _ in range(steps):
        a = math.floor(x)
        r = x - a
        anchors.append(a)
        rems.append(r)
        if not (-1.0 <= r <= 1.0):
            angles.append(float("nan"))
            break
        theta = math.degrees(math.acos(r))
        angles.append(theta)
        x = theta
    return anchors, rems, angles


def show(name: str, x0: float, steps: int = 14):
    anchors, rems, angles = peel(x0, steps)
    print(f"\n{name}  (start {x0:.10f})")
    print(f"  anchors: {anchors}")
    line = "  steps:  "
    for i, (a, r, th) in enumerate(zip(anchors, rems, angles)):
        line = (f"    step {i}: floor={a:3d}  frac={r:.8f}  "
                f"→ arccos = {th:.6f}°")
        print(line)


def main() -> None:
    print("=" * 78)
    print("Integer-peel angle recursion: x → arccos(frac(x))° , iterated")
    print("=" * 78)

    sigma_inf = 2.8628923          # best σ_0(∞) estimate (N=44 + geometric tail)
    sigma_7 = 2.0 * math.sqrt(2)   # σ_0(N=7), exact

    show("σ_0(∞) ≈ 2.8628923", sigma_inf)
    show("σ_0(7) = 2√2 (exact)", sigma_7)

    print()
    print("-" * 78)
    print("CONTROLS (generic numbers, for comparison):")
    show("π", math.pi)
    show("√3 + 1", math.sqrt(3) + 1)
    show("e", math.e)
    show("2.8642 (plain decimal near σ_0)", 2.8642)

    print()
    print("=" * 78)
    print("Reading (RETRACTED 2026-05-20): arccos∘frac is a mixing map, so an")
    print("aperiodic anchor sequence is the guaranteed output for every input")
    print("above, σ_0 and the controls alike. The test cannot distinguish")
    print("'σ_0 special' from 'σ_0 generic'; it is inconclusive by construction.")
    print("Compounding this, arccos∘frac amplifies error ~×3-5/step, so past")
    print("~step 10-12 the digits are limited by input precision (~1e-7).")
    print("σ_0 is characterised instead as a Schur-multiplier norm")
    print("(SigmaZeroCommutatorNormClaim, Tier1Derived).")


if __name__ == "__main__":
    main()
