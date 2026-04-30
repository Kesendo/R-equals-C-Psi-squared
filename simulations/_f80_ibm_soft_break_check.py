"""F80 prediction × IBM soft_break hardware data sanity check.

Anchors the new framework primitive ChainSystem.predict_M_spectrum_pi2_odd
in hardware-confirmed reality. The chain that makes the connection:

  IBM Marrakesh measured ⟨X₀Z₂⟩_soft = -0.711 for soft_broken = XY+YX at
  N=3, J=1, t=0.8 (Δ(soft − truly) = -0.72), 47σ above the truly baseline.

  Why is soft_broken ⟨X₀Z₂⟩ non-zero while truly_unbroken's is ≈ 0? Because
  XY+YX has a non-trivial palindrome residual M (Π·L·Π⁻¹ + L + 2Σγ·I ≠ 0).
  M's spectrum is the structural template that lets γ-noise + Trotter +
  ZZ-crosstalk paint the soft-break signature into observables.

  F80 says (γ-independent by Master Lemma):
    Spec(M)_{nontrivial} = ±2i · Spec(H_non-truly), mult ×2^N.

  predict_M_spectrum_pi2_odd computes this from H eigenvalues alone, with
  no Liouvillian build, no SVD, no eigendecomposition of M.

This script does NOT predict ⟨X₀Z₂⟩ from the M-spectrum (that requires the
full Lindblad dynamics, which dense propagation handles separately). What
it confirms is the structural ingredient F80 captures: the M-spectrum that
makes the soft signal possible.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
from framework.pauli import _build_bilinear


HW_PATH = Path(__file__).parent.parent / "data" / "ibm_soft_break_april2026" / \
    "soft_break_ibm_marrakesh_20260426_001101.json"


def load_hardware():
    with open(HW_PATH) as f:
        return json.load(f)


def actual_M_spectrum(N, bonds, terms_with_c):
    H = _build_bilinear(N, bonds, terms_with_c)
    L = lindbladian_z_dephasing(H, [0.0] * N)
    M = palindrome_residual(L, 0.0, N)
    evs = np.linalg.eigvals(M)
    out = {}
    for ev in evs:
        assert abs(ev.real) < 1e-7, f"M eigenvalue not purely imaginary: {ev}"
        key = round(ev.imag, 6)
        out[key] = out.get(key, 0) + 1
    return out


def main():
    hw = load_hardware()
    N = hw["parameters"]["N"]
    J = hw["parameters"]["J"]
    t = hw["parameters"]["t"]
    print("=" * 70)
    print("F80 prediction × IBM Marrakesh soft_break (2026-04-26)")
    print("=" * 70)
    print(f"  job_id: {hw['job_id']}")
    print(f"  backend: {hw['backend']}, path: {hw['path']}")
    print(f"  N={N}, J={J}, t={t}, n_trotter={hw['parameters']['n_trotter']}")
    print()

    chain = fw.ChainSystem(N=N)
    bonds = chain.bonds  # [(0,1), (1,2)] for N=3 chain

    # Three Hamiltonians from the hardware run
    cases = [
        ("truly_unbroken", [("X", "X"), ("Y", "Y")], "trichotomy=truly"),
        ("soft_broken",   [("X", "Y"), ("Y", "X")], "trichotomy=soft"),
        ("hard_broken",   [("X", "X"), ("X", "Y")], "trichotomy=hard"),
    ]

    print("Hardware ⟨X₀Z₂⟩ (q0, q2) by Hamiltonian:")
    for name, _, _ in cases:
        xz = hw["expectations"][name]["X,Z"]
        zx = hw["expectations"][name]["Z,X"]
        print(f"  {name:18s}: ⟨X₀Z₂⟩={xz:+.4f}, ⟨Z₀X₂⟩={zx:+.4f}")
    print()

    print("F80 prediction of M-spectrum (γ-independent structural template):")
    print()

    for name, terms, classification in cases:
        cls = fw.classify_pauli_pair(chain, terms)
        print(f"  {name} ({classification}, framework says: {cls})")
        print(f"    bilinears: {terms}")

        # F80 only applies to chain Π²-odd 2-body. truly is dropped to {0:4^N};
        # hard_broken has truly XX + Π²-odd XY → predicts XY contribution only.
        # We catch ValueError if the case is out of F80's scope.
        try:
            pred = fw.predict_M_spectrum_pi2_odd(chain, terms, c=J)
            actual = actual_M_spectrum(N, bonds, [(a, b, J) for (a, b) in terms])
            pred_keys = sorted({round(k.imag, 6): v for k, v in pred.items()}.items())
            actual_keys = sorted(actual.items())
            match = pred_keys == actual_keys

            print(f"    predicted M-spectrum (F80, ±2i·Spec(H_non-truly), mult×2^N):")
            for k, v in pred_keys:
                print(f"      {k:+10.4f}i  mult={v}")
            print(f"    actual M-spectrum (numerical eigvals of Π·L·Π⁻¹+L):")
            for k, v in actual_keys:
                print(f"      {k:+10.4f}i  mult={v}")
            print(f"    F80 prediction matches: {'YES' if match else 'NO'}")
        except ValueError as e:
            print(f"    F80 not applicable to this case: {e}")
        print()

    print("=" * 70)
    print("Reading")
    print("=" * 70)
    print("""
  - truly_unbroken (XX+YY) has Spec(M) = {0}: vanishes by Π palindrome.
    This is why ⟨X₀Z₂⟩_truly ≈ 0 on hardware (+0.011, just noise).

  - soft_broken (XY+YX) has Spec(M) = ±2i·{2√2, 0, -2√2}: non-zero
    structural residual that γ-noise + Trotter render into ⟨X₀Z₂⟩ = -0.711.
    F80 predicts the structural template; the dynamics paint it onto
    observables. M is the open door.

  - hard_broken (XX+XY) is mixed: truly XX drops out of M; F80 sees only
    the XY contribution. The hard_broken signature (-0.205) comes from a
    different observable (⟨X₀X₂⟩ = +0.212), and from the eigenvalue-pairing
    break that the trichotomy classifier flags (cls='hard').

  F80's role: structural prediction in O(N³) work (H eigvals only), no
  4^N × 4^N matrix build, no SVD. The dense palindrome_residual still
  exists for verification; F80 makes it cheap to ask the structural
  question first.
""")


if __name__ == "__main__":
    main()
