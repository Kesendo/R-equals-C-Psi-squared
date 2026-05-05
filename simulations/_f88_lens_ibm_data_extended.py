"""F88-Lens applied to repository-internal IBM Marrakesh hardware data.

Two datasets in `data/`:

1. `data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json`
   F87 trichotomy replication on path [48, 49, 50] (different qubits than the
   earlier framework_snapshots run analyzed by _f88_lens_ibm_framework_snapshots.py
   which used path [0, 1, 2]). Same-day independent run; tests whether the
   F88-Lens trichotomy separation (truly ~0 / soft ~0.7 / hard ~0.3) reproduces
   on different physical qubits.

2. `data/ibm_zn_mirror_april2026/zn_mirror_ibm_marrakesh_20260429_102824.json`
   Heisenberg evolution of |+−+⟩ vs |−+−⟩ = Z⊗N |+−+⟩. Heisenberg is a truly
   Hamiltonian (XX+YY+ZZ, all Π²-even bilinears) and Z⊗N-symmetric, so both
   states should evolve to ρ's with the same Π²-odd/memory content. Hardware
   drift breaks this; we measure how cleanly F88-Lens reads it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _f88_lens_ibm_framework_snapshots import reconstruct_2qubit_rho, f88_lens_2qubit


REPO_ROOT = Path(__file__).resolve().parents[1]
SOFT_BREAK_JSON = REPO_ROOT / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_marrakesh_20260426_001101.json"
ZN_MIRROR_JSON = REPO_ROOT / "data" / "ibm_zn_mirror_april2026" / "zn_mirror_ibm_marrakesh_20260429_102824.json"


def main():
    print("F88-Lens on repository-internal IBM Marrakesh data")
    print("=" * 78)

    # ── Soft-break trichotomy replication (path [48, 49, 50]) ──
    print()
    print(f"=== Dataset 1: soft_break_marrakesh_20260426_001101 (path [48, 49, 50]) ===")
    with open(SOFT_BREAK_JSON, "r", encoding="utf-8") as fh:
        sb = json.load(fh)
    print(f"  job {sb['job_id']}, {sb['parameters']['shots']} shots/basis, t={sb['parameters']['t']}")
    print()
    print(f"  {'category':<20} {'trace':>8} {'purity':>8} {'static':>9} {'memory':>9} {'Π²-odd/mem':>12}")
    print("  " + "-" * 70)
    sb_results = {}
    for cat in ["truly_unbroken", "soft_broken", "hard_broken"]:
        rho = reconstruct_2qubit_rho(sb["expectations"][cat])
        lens = f88_lens_2qubit(rho)
        sb_results[cat] = lens
        short = cat.replace("_unbroken", "").replace("_broken", "")
        print(f"  {short:<20} {lens['trace']:>8.4f} {lens['purity']:>8.4f} "
              f"{lens['static_frac']:>9.4f} {lens['memory_frac']:>9.4f} {lens['pi2_odd_in_memory']:>12.4f}")

    print()
    print("  Replication of F87 trichotomy state-level differentiation on a")
    print("  different physical path. Earlier run (path [0, 1, 2], framework_snapshots")
    print("  20260426_105948): truly=0.0297 / soft=0.7444 / hard=0.2763")
    truly = sb_results["truly_unbroken"]["pi2_odd_in_memory"]
    soft = sb_results["soft_broken"]["pi2_odd_in_memory"]
    hard = sb_results["hard_broken"]["pi2_odd_in_memory"]
    print(f"  This run  (path [48, 49, 50]):                             "
          f"truly={truly:.4f} / soft={soft:.4f} / hard={hard:.4f}")
    if soft > 1e-6:
        print(f"  Separation truly/soft this run: {soft / max(truly, 1e-9):.1f}× (earlier run: 25×)")

    # ── Z⊗N-mirror two-state Heisenberg test ──
    print()
    print(f"=== Dataset 2: zn_mirror_marrakesh_20260429_102824 (path [48, 49, 50]) ===")
    with open(ZN_MIRROR_JSON, "r", encoding="utf-8") as fh:
        zm = json.load(fh)
    print(f"  job {zm['job_id']}, {zm['parameters']['shots']} shots/basis, t={zm['parameters']['t']}")
    print(f"  Heisenberg J={zm['parameters']['J']} on |+−+⟩ (ρ_a) vs |−+−⟩ = Z⊗N |+−+⟩ (ρ_b)")
    print(f"  truly Hamiltonian → F88 predicts Π²-odd/memory ≈ 0 for both,")
    print(f"  Z⊗N-symmetry → both states should give the same Π²-odd/memory exactly")
    print()
    print(f"  {'state':<10} {'trace':>8} {'purity':>8} {'static':>9} {'memory':>9} {'Π²-odd/mem':>12}")
    print("  " + "-" * 60)
    state_keys = list(zm["expectations"].keys())
    zm_results = {}
    for state_key in state_keys:
        rho = reconstruct_2qubit_rho(zm["expectations"][state_key])
        lens = f88_lens_2qubit(rho)
        zm_results[state_key] = lens
        print(f"  {state_key:<10} {lens['trace']:>8.4f} {lens['purity']:>8.4f} "
              f"{lens['static_frac']:>9.4f} {lens['memory_frac']:>9.4f} {lens['pi2_odd_in_memory']:>12.4f}")

    if len(state_keys) == 2:
        a_odd = zm_results[state_keys[0]]["pi2_odd_in_memory"]
        b_odd = zm_results[state_keys[1]]["pi2_odd_in_memory"]
        print()
        print(f"  Z⊗N-symmetry violation in Π²-odd/mem (|ρ_a − ρ_b| / max): "
              f"{abs(a_odd - b_odd):.4f} / {max(a_odd, b_odd):.4f}")
        print(f"  Both values should be ≈ 0 (truly Heisenberg). Hardware noise floor")
        print(f"  for path [48, 49, 50] from soft_break replication: truly = {truly:.4f}.")

    print()
    print("=== Reading ===")
    print()
    print(f"Soft-break replication (path [48, 49, 50]):")
    print(f"  truly = {truly:.4f}, soft = {soft:.4f}, hard = {hard:.4f}.")
    print(f"  Compare to earlier run on path [0, 1, 2]: 0.0297 / 0.7444 / 0.2763.")
    print(f"  Same-day, different-qubits replication of the F88-Lens F87 trichotomy")
    print(f"  state-level differentiation; the ordering and order-of-magnitude")
    print(f"  separation should hold across paths if the framework prediction is")
    print(f"  hardware-substrate-independent.")
    print()
    print(f"Z⊗N-mirror two-state symmetry test:")
    print(f"  Heisenberg is a truly (Π²-even) Hamiltonian; F88 predicts Π²-odd/mem")
    print(f"  ≈ 0 for both ρ_a and ρ_b (and equal to each other under Z⊗N-symmetry).")
    print(f"  Deviation captures hardware Π²-odd noise floor under different state")
    print(f"  preparations on the same path.")


if __name__ == "__main__":
    main()
