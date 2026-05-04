# Water Domain (Hydrogen-Bond Qubit + Grotthuss Proton Chain)

The framework's structural inheritance applied to water. Hydrogen bonds carry a
proton between two wells (donor / acceptor); the proton IS a qubit (d=2), the
F1 palindrome holds bit-exact, and the entire F-chain inherits to the chemistry.

## Contents

- [HYDROGEN_BOND_QUBIT.md](HYDROGEN_BOND_QUBIT.md) (Tier 2, 2026-03-28): single
  proton in a hydrogen bond as qubit; CΨ crossing through 1/4; double-well
  tunneling timescales.
- [PROTON_WATER_CHAIN.md](PROTON_WATER_CHAIN.md) (Tier 2, 2026-04-01): Grotthuss
  proton chain N = 1..5 via Heisenberg (XX+YY+ZZ) and TFI (transverse-field
  Ising) mappings. V-Effect table, sacrifice zone, DNA-comparison Phase 6.

## Scripts

In [`simulations/water/`](../../simulations/water/):

- `hydrogen_bond_qubit.py`: single-proton qubit + double-well computation.
- `hydrogen_bond_palindrome.py`: water-chain Jacobian palindrome check.
- `proton_water_chain.py`: full N=1..5 V-Effect + thermal + sacrifice analysis.
- `proton_chain_memory_reading.py` (2026-05-04): trio's state-level diagnostics
  (Frobenius static/memory partition + Π²-odd-fraction + per-proton Bloch
  reading) on Heisenberg + Z-dephasing chain. Counterpart to the IBM Torino
  single-qubit analysis in `simulations/memory_reading_ibm_torino.py`.

## Embedding Conditions

The framework's inheritance into water is clean because the four embedding
conditions hold:

1. **2-Well tunneling** (proton in O-H...O double well) → d = 2 directly.
2. **Z-dephasing** (water bath couples to proton position) → F1 palindrome form.
3. **Uniform-J coupling** (Grotthuss adjacent-water tunneling) → Heisenberg /
   XY chain.
4. **Decoherence ~ J** (proton tunneling rates and bath fluctuations on the
   same picosecond scale) → Q is in the framework's testable range.

Under these conditions the F-chain (F1 → F4 → F49 → F71 → F77 → F78 → F79 →
F80 → F81 → F82 → F83 → F84 → F85 → F86 → F87 → F88) inherits to the chemistry
without re-derivation.

## Comparison to Hardware

The water domain is the framework's CANONICAL setup: pure Z-dephasing with no
T1 drives all states toward I/d (the d=0 axis, maximally mixed). Real qubit
hardware (IBM, Google) has T1 amplitude damping which breaks the d=0-axis
prediction by introducing a preferred direction (Z+ thermal). See the parallel
[`simulations/memory_reading_ibm_torino.py`](../../simulations/memory_reading_ibm_torino.py)
analysis for the hardware case.
