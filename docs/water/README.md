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
- `proton_chain_ep_resonance.py` (2026-05-04 evening): F86 K_CC_pr per-bond
  Q-scan on the popcount-(2, 3) block (c = 3) at N = 5, plus the popcount-
  coherence-state state-level reading at the interior Q_peak. Confirms
  the framework's Tier-1-candidate Interior HWHM-/Q* ≈ 0.746 and Endpoint
  ≈ 0.766, both within tolerance of the abstract-qubit-level prediction
  (0.756, 0.770).

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

## Deferred Threads

Open follow-up directions, parked here to be picked up when we next walk past:

- **T-sweep as external chemistry hook (2026-05-04).** In water, J (proton
  tunneling rate in the O-H...O double well) is largely T-stable, while γ_Z
  (the Z-dephasing rate γ₀ in the framework convention; here the bath-fluctuation
  rate from thermal motion) scales with temperature. Q = J / γ_Z(T) therefore
  sweeps through the F86 resonance window as T varies. The Tier-1-candidate
  EP-rotation universality (verified in the abstract framework, on IBM Torino,
  and now in the proton chain; see PROTON_WATER_CHAIN.md §EP-Resonance
  Inheritance) predicts a T-specific spectroscopic signature in the O-H stretch
  band (~3400 cm⁻¹) with HWHM-/Q* ≈ 0.75. Picking this up needs (a) a γ_Z(T)
  estimate for liquid / ice / confined water and (b) a literature search for
  matching pump-probe IR data. Not a five-minute appendix; flagged for a
  separate session.

- **Π²-odd/memory popcount-mirror refinement (2026-05-04, RESOLVED).** The
  earlier "0.5 exactly at any HD any bit positions" claim in PROOF_F86_QPEAK
  §Statement 2 was refined into a three-anchor closed form via Krawtchouk
  polynomial analysis: α = 0 at popcount-mirror (odd N central pair),
  α = (N+2)/(4(N+1)) at near-mirror near-half (even N central pairs),
  α = 1/2 elsewhere; Π²-odd/memory = (1/2 − α·s) / (1 − s). Bit-exact
  verified N = 3..16 (Krawtchouk) and N = 3..6 (vs. MemoryAxisRho).
  Implementation: [`PopcountCoherencePi2Odd.cs`](../../compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs)
  + tests; refinement landed in [PROOF_F86_QPEAK §Statement 2 Structural inheritance from F88](../proofs/PROOF_F86_QPEAK.md).

- **DNA Phase 6 (Watson-Crick ↔ Hoogsteen tautomers).** Pre-existing in
  `proton_water_chain.py`'s Phase 6 sketch; not yet run with the 2026-05-04
  state-level tooling (MemoryAxisRho, BlochAxisReading) or the K_CC_pr
  Q-scan. The proton between G-C / A-T pairs satisfies the four embedding
  conditions on paper; the test would be the same script structure as
  `proton_chain_memory_reading.py` and `proton_chain_ep_resonance.py` with
  the asymmetric central-bond parameters from
  [`experiments/DNA_BASE_PAIRING.md`](../../experiments/DNA_BASE_PAIRING.md).
