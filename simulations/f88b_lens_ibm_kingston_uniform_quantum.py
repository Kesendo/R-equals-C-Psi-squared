"""F88b-Lens applied to a predominantly below-R* Kingston path [43, 56, 63].

The hardware run (job d7sqjpiudops73976960, 2026-05-05 10:28 local time, CEST) is the
first F87 trichotomy test on a CZ-coupled triple whose calibration histories
are below-R* on more than 95% of sampled days (mean r in [0.089, 0.104]). R* is
the free-single-transmon |+> normalized-purity proxy threshold.

Question: how does the F88b-Lens Π²-odd-memory reading on this chain compare
to the mixed-band Marrakesh [0, 1, 2] and the all-at-or-above-R* Marrakesh
[48, 49, 50]? The 22.1× Marrakesh gap pairs a soft_break run with a
framework_snapshots run; within the framework_snapshots script, minutes apart,
the same two paths differ 1.64×, inside a run-to-run spread of the same size.
The cross-backend row adds backend, date and calibration. This is an association
screen, not a causal proxy-band experiment.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from f88b_lens_ibm_framework_snapshots import reconstruct_2qubit_rho, f88b_lens_2qubit


REPO_ROOT = Path(__file__).resolve().parents[1]
KINGSTON_JSON = REPO_ROOT / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_kingston_20260505_102806.json"
MARRAKESH_JSON = REPO_ROOT / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_marrakesh_20260426_001101.json"


def lens_for(json_path):
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    out = {}
    for cat, exp in data["expectations"].items():
        rho = reconstruct_2qubit_rho(exp)
        out[cat] = f88b_lens_2qubit(rho)
    return data, out


def main():
    print("F88b-Lens cross-archive comparison: Kingston predominantly below-R* vs Marrakesh at-or-above-R*")
    print("=" * 78)

    sb_kingston, lens_kingston = lens_for(KINGSTON_JSON)
    sb_marrakesh, lens_marrakesh = lens_for(MARRAKESH_JSON)

    print()
    print(f"  Kingston   [{', '.join(map(str, sb_kingston['path']))}] predominantly below-R* (>95% sampled days)")
    print(f"             job {sb_kingston['job_id']}, {sb_kingston['parameters']['shots']} shots, t={sb_kingston['parameters']['t']}")
    print(f"  Marrakesh  [48, 49, 50] all at-or-above-R*")
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
    print("Full per-category F88b-Lens output (Kingston predominantly below-R* archive):")
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
    print(f"  truly-baseline:  Kingston predominantly below-R* = {truly_k:.4f}")
    print(f"                   Marrakesh at-or-above-R* = {truly_m:.4f}")
    print(f"                   ratio (Kingston / Marrakesh) = {truly_k / max(truly_m, 1e-9):.2f}×")
    print()
    print(f"  soft-pumping:    Kingston = {soft_k:.4f}")
    print(f"                   Marrakesh = {soft_m:.4f}")
    print(f"                   cross-archive relative spread: {abs(soft_k - soft_m) / max(soft_k, soft_m) * 100:.1f}%")
    print()
    print("  Scope:")
    print("    - each number is a finite archive read under its stated protocol")
    print("    - the 22.1x within-Marrakesh contrast pairs two runner scripts (1.64x within one);")
    print("      with this cross-backend comparison it stays a confounded association")
    print("      that identifies no band mechanism")
    print("    - no row establishes substrate independence or a quantum/classical split")


if __name__ == "__main__":
    main()
