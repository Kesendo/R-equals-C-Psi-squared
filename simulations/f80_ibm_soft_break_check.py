"""F80 prediction on the three Hamiltonians of the IBM soft_break run.

Checks the framework primitive ChainSystem.predict_M_spectrum_pi2_odd against
the M-spectrum built from the dense Liouvillian, for the three Hamiltonians
IBM Marrakesh flew (truly XX+YY, soft XY+YX, hard XX+XY at N=3, J=1).

  F80 says (γ-independent by Master Lemma):
    Spec(M)_{nontrivial} = ±2i · Spec(H_non-truly), mult ×2^N.

  predict_M_spectrum_pi2_odd computes this from H eigenvalues alone, with
  no Liouvillian build, no SVD, no eigendecomposition of M.

The M-spectrum does not explain the run's ⟨X₀Z₂⟩ values. XY+YX is XX+YY
with the middle qubit turned by (X+Y)/√2, which keeps the Z-dephasing and
leaves q0 and q2 alone, so soft from |+−+⟩ is truly from |+,−i,+⟩ on every
two-qubit Pauli: soft's predicted value (-0.623 continuous, -0.723 Trotter,
against -0.711 measured) is truly's evolution from a turned middle state, and M is nonzero for XY+YX only because Π is written in one frame
(experiments/MARRAKESH_THREE_LAYERS.md, Layer 3).
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

  - soft_broken (XY+YX) has Spec(M) = ±2i·{2√2, 0, -2√2}. This does not
    make its ⟨X₀Z₂⟩ = -0.711: XY+YX is XX+YY with the middle qubit turned,
    and truly from |+,−i,+⟩ has the same predicted value (-0.623 continuous,
    -0.723 Trotter); the zeros of the table come from symmetries of
    state and generator, not from M.

  - hard_broken (XX+XY) is mixed: truly XX drops out of M; F80 sees only
    the XY contribution.

  F80's role: structural prediction in O(N³) work (H eigvals only), no
  4^N × 4^N matrix build, no SVD. The dense palindrome_residual still
  exists for verification; F80 makes it cheap to ask the structural
  question first.
""")


if __name__ == "__main__":
    main()
