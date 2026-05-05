"""F88-Lens applied to today's Kingston run on uniform-quantum chain [43, 56, 63].

The hardware run (job d7sqjpiudops73976960, 2026-05-05 10:28 UTC) is the
first F87 trichotomy test on a uniform-quantum CZ-coupled triple. Multi-day
biographies confirm all three qubits PulseStable across 91 days (mean r
in [0.09, 0.10], all crossing the boundary on > 95% of days, deeply
quantum-side).

Question: how does the F88-Lens Π²-odd-memory reading on this chain compare
to the regime-mixed Marrakesh [0, 1, 2] (April 26, framework_snapshots) and
the uniform-classical Marrakesh [48, 49, 50] (April 26, soft_break)? The
BOTH_SIDES_VISIBLE update flagged the 22.8× truly-baseline gap (Marrakesh
[48, 49, 50] = 0.0013 vs Marrakesh [0, 1, 2] = 0.0297) as a possible
regime-uniformity effect; this run is the first uniform-quantum data point.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _f88_lens_ibm_framework_snapshots import reconstruct_2qubit_rho, f88_lens_2qubit


REPO_ROOT = Path(__file__).resolve().parents[1]
KINGSTON_JSON = REPO_ROOT / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_kingston_20260505_102806.json"
MARRAKESH_JSON = REPO_ROOT / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_marrakesh_20260426_001101.json"


def lens_for(json_path):
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    out = {}
    for cat, exp in data["expectations"].items():
        rho = reconstruct_2qubit_rho(exp)
        out[cat] = f88_lens_2qubit(rho)
    return data, out


def main():
    print("F88-Lens cross-regime comparison: Kingston uniform-quantum vs Marrakesh uniform-classical")
    print("=" * 78)

    sb_kingston, lens_kingston = lens_for(KINGSTON_JSON)
    sb_marrakesh, lens_marrakesh = lens_for(MARRAKESH_JSON)

    print()
    print(f"  Kingston   [{', '.join(map(str, sb_kingston['path']))}] uniform-quantum")
    print(f"             job {sb_kingston['job_id']}, {sb_kingston['parameters']['shots']} shots, t={sb_kingston['parameters']['t']}")
    print(f"  Marrakesh  [48, 49, 50] uniform-classical")
    print(f"             job {sb_marrakesh['job_id']}, {sb_marrakesh['parameters']['shots']} shots, t={sb_marrakesh['parameters']['t']}")
    print()
    print(f"  {'category':<22} {'Kingston':>12} {'Marrakesh':>12} {'ratio':>8}")
    print("  " + "-" * 60)

    # Map category names: Kingston has truly_unbroken, pi2_odd_pure, pi2_even_nontruly, mixed_anti_one_sixth
    # Marrakesh has truly_unbroken, soft_broken, hard_broken
    pairs = [
        ("truly_unbroken", "truly_unbroken", "truly"),
        ("pi2_odd_pure", "soft_broken", "Π²-odd / soft"),
        ("mixed_anti_one_sixth", "hard_broken", "mixed / hard"),
    ]
    for k_cat, m_cat, label in pairs:
        if k_cat in lens_kingston and m_cat in lens_marrakesh:
            k = lens_kingston[k_cat]["pi2_odd_in_memory"]
            m = lens_marrakesh[m_cat]["pi2_odd_in_memory"]
            ratio = k / m if m > 1e-9 else float("inf")
            print(f"  {label:<22} {k:>12.4f} {m:>12.4f} {ratio:>8.2f}×")

    print()
    print("Full per-category F88-Lens output (Kingston uniform-quantum):")
    print(f"  {'category':<22} {'trace':>8} {'purity':>8} {'static':>9} {'memory':>9} {'Π²-odd/mem':>12}")
    print("  " + "-" * 70)
    for cat in ["truly_unbroken", "pi2_odd_pure", "pi2_even_nontruly", "mixed_anti_one_sixth"]:
        l = lens_kingston[cat]
        short = cat.replace("_unbroken", "").replace("_pure", "").replace("_anti_one_sixth", "")
        print(f"  {short:<22} {l['trace']:>8.4f} {l['purity']:>8.4f} "
              f"{l['static_frac']:>9.4f} {l['memory_frac']:>9.4f} {l['pi2_odd_in_memory']:>12.4f}")

    print()
    print("=" * 78)
    print("READING")
    print()
    truly_k = lens_kingston["truly_unbroken"]["pi2_odd_in_memory"]
    truly_m = lens_marrakesh["truly_unbroken"]["pi2_odd_in_memory"]
    soft_k = lens_kingston["pi2_odd_pure"]["pi2_odd_in_memory"]
    soft_m = lens_marrakesh["soft_broken"]["pi2_odd_in_memory"]
    print(f"  truly-baseline:  Kingston uniform-quantum  = {truly_k:.4f}")
    print(f"                   Marrakesh uniform-classical = {truly_m:.4f}")
    print(f"                   ratio (Kingston / Marrakesh) = {truly_k / max(truly_m, 1e-9):.2f}×")
    print()
    print(f"  soft-pumping:    Kingston = {soft_k:.4f}")
    print(f"                   Marrakesh = {soft_m:.4f}")
    print(f"                   substrate-independence: {abs(soft_k - soft_m) / max(soft_k, soft_m) * 100:.1f}% relative spread")
    print()
    print("  Interpretation guide:")
    print("    - if regime-uniformity is what cleans the truly-baseline:")
    print("        Kingston uniform-quantum truly should be similar to Marrakesh uniform-classical")
    print("        (the 22.8× framework_snapshots [0,1,2] vs [48,49,50] gap was 0.0297 vs 0.0013)")
    print("    - if quantum-side dephasing dominates:")
    print("        Kingston uniform-quantum truly should be LARGER than Marrakesh classical")
    print("        (more T2-driven dephasing accumulating in the truly-baseline noise floor)")
    print("    - if the chip's overall noise floor matters most:")
    print("        the truly-baseline reflects Kingston's calibration quality, not regime")


if __name__ == "__main__":
    main()
