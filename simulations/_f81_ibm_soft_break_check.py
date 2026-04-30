"""F81 lens on the IBM Marrakesh soft-break dataset.

Companion to `_f80_ibm_soft_break_check.py`. Where F80 reads the M-spectrum
(structural strength), F81 reads the Π-decomposition: how M splits into the
Π-symmetric "remembered" part (mirror image + dissipator + Π²-even
commutator + dissipation shift) and the Π-antisymmetric "driving" part
(L_{H_odd}, the unitary commutator induced by the Π²-odd Pauli bilinears).

The three Marrakesh Hamiltonians give three qualitatively distinct
Π-decompositions:

  truly XX+YY  : M = 0,             no Π-antisymmetric drive at all
  soft  XY+YX  : 50/50, full drive, ‖M_anti‖² = full L_H Frobenius
  hard  XX+XY  : 50/50, half drive, ‖M_anti‖² = half of soft's
                                    (only XY contributes; XX is truly)

Both soft and hard satisfy the F81 50/50 condition (PROOF_F81 Step 8:
all non-truly bilinears Π²-odd). They differ in absolute scale: hard's
single Π²-odd bilinear gives half the drive of soft's two.

This is the structural reason behind the hardware trichotomy. ⟨X₀Z₂⟩
becomes non-zero exactly when the Π²-odd part of H provides drive; the
larger that drive (in Frobenius norm) the more open the "door" that γ-noise
+ Trotter discretization paint observable signals through.

This script does NOT predict ⟨X₀Z₂⟩ from M_anti (that requires full Lindblad
dynamics, which exp(L·t) handles separately). What it shows is the
structural decomposition that the hardware trichotomy directly reflects.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw
from framework.pauli import _build_bilinear


HW_PATH = Path(__file__).parent.parent / "data" / "ibm_soft_break_april2026" / \
    "soft_break_ibm_marrakesh_20260426_001101.json"


def main():
    with open(HW_PATH) as f:
        hw = json.load(f)
    N = hw["parameters"]["N"]
    J = hw["parameters"]["J"]
    t = hw["parameters"]["t"]
    print("=" * 72)
    print("F81 lens on IBM Marrakesh soft_break (2026-04-26)")
    print("=" * 72)
    print(f"  job_id: {hw['job_id']}")
    print(f"  backend: {hw['backend']}, path: {hw['path']}")
    print(f"  N={N}, J={J}, t={t}, n_trotter={hw['parameters']['n_trotter']}")
    print()

    chain = fw.ChainSystem(N=N, J=J)
    cases = [
        ("truly_unbroken", [("X", "X"), ("Y", "Y")], "truly"),
        ("soft_broken",   [("X", "Y"), ("Y", "X")], "soft (Π²-odd, pure)"),
        ("hard_broken",   [("X", "X"), ("X", "Y")], "hard (mixed even+odd)"),
    ]

    print("Hardware ⟨X₀Z₂⟩ (q0, q2) and F81 Π-decomposition of M:")
    print()
    print(f"  {'Hamiltonian':<18} | {'⟨X₀Z₂⟩':>9} | {'‖M‖²':>9} | {'‖M_sym‖²':>10} | "
          f"{'‖M_anti‖²':>11} | {'anti/total':>11}")
    print("  " + "-" * 90)
    decompositions = {}
    for name, terms, label in cases:
        d = chain.pi_decompose_M(terms, gamma_z=0.0)
        decompositions[name] = d
        xz_hw = hw["expectations"][name]["X,Z"]
        ns = d["norm_sq"]
        anti_frac = ns["M_anti"] / ns["M"] if ns["M"] > 1e-12 else 0.0
        print(f"  {name:<18} | {xz_hw:+9.4f} | {ns['M']:>9.2f} | {ns['M_sym']:>10.2f} | "
              f"{ns['M_anti']:>11.2f} | {anti_frac:>11.4f}")
    print()

    # The F81 identity is enforced internally by pi_decompose_M (it raises
    # RuntimeError if violated). If we got here, all three identities hold.
    print("F81 identity Π·M·Π⁻¹ = M − 2·L_{H_odd} verified internally for all three.")
    print()

    # Read the structural picture
    print("=" * 72)
    print("Structural reading")
    print("=" * 72)

    soft_M_anti = decompositions["soft_broken"]["norm_sq"]["M_anti"]
    hard_M_anti = decompositions["hard_broken"]["norm_sq"]["M_anti"]
    truly_M = decompositions["truly_unbroken"]["norm_sq"]["M"]

    print(f"""
  truly_unbroken (XX+YY): M = 0 entirely. No Π²-odd Pauli bilinears in H,
    so L_{{H_odd}} = 0, hence M_anti = 0. Π closes the palindrome perfectly;
    no drive, no remembered residue. Hardware ⟨X₀Z₂⟩ ≈ 0 (measured: +0.011)
    is a direct consequence: nothing in M to project onto observables.

  soft_broken (XY+YX): Pure Π²-odd. M_anti = L_H = -i[H, ·] entirely.
    The 50/50 split ‖M_sym‖² = ‖M_anti‖² = {soft_M_anti:.0f} is structural
    (analytical, see PROOF_F81 Step 8). The whole dynamics generator drives
    M's Π-antisymmetric component. Hardware ⟨X₀Z₂⟩ = -0.711 is the strongest
    signature in the trichotomy: maximally open door.

  hard_broken (XX+XY): Mixed truly + Π²-odd. The truly XX bilinear contributes
    zero to M (Master Lemma); the XY bilinear drives the entire residual.
    Internally, the same 50/50 ‖M_sym‖² = ‖M_anti‖² split holds (PROOF_F81
    Step 8 covers truly + Π²-odd combinations, not just pure Π²-odd), but
    hard's M is structurally smaller than soft's because hard contributes
    only one Π²-odd bilinear (XY) versus soft's two (XY+YX): ‖M‖²_hard =
    {decompositions['hard_broken']['norm_sq']['M']:.0f} = half of soft's
    {decompositions['soft_broken']['norm_sq']['M']:.0f}. Hardware ⟨X₀Z₂⟩ =
    +0.205, weaker in magnitude than soft (-0.711) and with opposite sign.

  Hardware reading: the magnitude ranking |⟨X₀Z₂⟩|: 0.711 (soft) > 0.205
    (hard) > 0.011 (truly) tracks the ranking of ‖M_anti‖²: {soft_M_anti:.0f}
    (soft) > {hard_M_anti:.0f} (hard) > {truly_M:.0f} (truly = 0). The
    observable magnitude is not linearly proportional to ‖M_anti‖ (the
    Lindblad evolution exp(L·t) is non-linear in L), but the QUALITATIVE
    ordering is exact: more Π-antisymmetric drive → stronger soft-break
    signature.

  What F81 adds beyond F80: F80 says Spec(M) (the strength dial). F81 says
    *which half of M does the driving* (Π-antisymmetric component) and
    *which half is remembered* (Π-symmetric component). Together they give
    a complete structural picture of M.
""")


if __name__ == "__main__":
    main()
