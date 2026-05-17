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

## Findings since May 4

### F86b 3/8 K-intermediate anchor → 1/4 asymptote bridge (2026-05-17 evening)

The F86b 3/8 Dicke-K-intermediate anchor (Tier 1 derived this morning via X⊗N-
eigenbasis decomposition, commit `b9ba5f6`, `compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs`)
states: the Dicke superposition `ψ = (|D_{N/2−1}⟩ + |D_{N/2}⟩) / √2` on even N has
Π²-odd Frobenius² total `α_total = 3/8` at t = 0 exactly.

Testing the inheritance to the proton water chain (Heisenberg + Z-dephasing,
`simulations/water/proton_chain_dicke_anchor.py`) finds the closed-form bit-
exactly at t = 0 — F86b inheritance to the chemistry-grounded substrate
confirmed, as expected from the four embedding conditions.

The new finding is the *long-time* behaviour under L evolution:

```
α(∞)_KIntermediate(N even) = (N + 2) / [4·(N + 1)]
```

Bit-exact verified N = 4..16 via the script. Two ingredient closed forms enter the
derivation:

1. `‖P_{N/2−1}_odd‖² = C(N, N/2−1) / 2` — Π²-odd Frobenius² of the sub-mid sector
   projector is **exactly half its rank**. Verified N = 4..16 via direct Krawtchouk
   enumeration over Pauli-string supports.
2. `‖P_{N/2}_odd‖² = 0` — mid-popcount Krawtchouk parity vanishing (odd-k
   `K_k(N/2; N) = 0`). Standard.

Combined via the kernel projection `ρ_∞ = (1/2)/C(N, m)·P_m + (1/2)/C(N, m+1)·P_{m+1}`:

```
α(∞) = C(N, N/2) / [2·C(N+1, N/2)] = (N+2) / (4(N+1))
```

The asymptote at N → ∞ is exactly **1/4** — the `HalfAsStructuralFixedPoint²` =
Mandelbrot cardioid maxval = Theorem 2 ceiling = CΨ fold boundary documented in
`compute/RCPsiSquared.Core/Symmetry/QuarterAsBilinearMaxval.cs` and tied to
`compute/RCPsiSquared.Core/Symmetry/HalfAsStructuralFixedPoint.cs`. The morning's
*static* 3/8 anchor and the framework's *universal* 1/4 are connected by an
explicit N-dependent decay curve traversed by KIntermediate states under Heisenberg-
XY + Z-dephasing dynamics:

```
N=4:  α(0) = 3/8 = 0.375  → α(∞) = 3/10 ≈ 0.300   (Δ to 1/4 = 1/20)
N=10: α(0) = 3/8 = 0.375  → α(∞) = 3/11 ≈ 0.273   (Δ to 1/4 = 1/44)
N=20: α(0) = 3/8 = 0.375  → α(∞) = 11/42 ≈ 0.262  (Δ to 1/4 = 1/84)
N→∞:  α(0) = 3/8          → α(∞) → 1/4            [universal boundary]
```

The water-chain experiment uncovered this because it asked a question that the
abstract framework hadn't asked: not "what is α_total at t = 0" (the F86b
question, closed in static form) but "what is α_total at t = ∞ under truly-class
evolution". The four embedding conditions guarantee the inheritance the other way
too: the closed form is universal across any truly-class Hamiltonian + Z-dephasing
on N qubits, not specific to chain XY — the bond topology drops out because the
long-time limit projects onto `ker L = span(P_0, …, P_N)` for any connected graph
(per `F4`).

**Two new closed forms** (Tier 1 derived, bit-exact verified N=4..16):

- `‖P_{N/2−1}_odd‖² = C(N, N/2−1) / 2`
- `α(∞)_KIntermediate(N even) = (N + 2) / [4·(N + 1)]`

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
