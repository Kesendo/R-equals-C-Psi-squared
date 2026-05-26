# polarity_fingerprint on Tier-B Marrakesh/Kingston Hamiltonians

**Status:** Welle 6 application of the typed `fw.polarity_fingerprint` workflow to 11 real-hardware-tested Hamiltonians across 3 Tier-B datasets. Extends the F87↔F112 orthogonality empirical anchor from 1 synthetic instance to 11 hardware-instances; surfaces 3 bit_b-inhomogeneous-H cases that are still F112 BALANCED bit-exact, confirming the "F112 typed scope is sufficient but not necessary" reading.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Script:** [`simulations/_polarity_fingerprint_tierB_marrakesh.py`](../simulations/_polarity_fingerprint_tierB_marrakesh.py)
**Workflow:** [`fw.polarity_fingerprint`](../simulations/framework/workflows/polarity_fingerprint.py) (added Welle 5.B)

## Setup

Tier-B Marrakesh/Kingston datasets are single-time-snapshot 2-qubit tomography (16 Pauli expectations per Hamiltonian, no t-series for L_eff fitting). The polarity_fingerprint workflow operates on (chain, terms) inputs at the framework level (not on trajectories), so it applies to the Hamiltonians used in each dataset without needing the snapshot data itself for the framework reading.

11 Hamiltonian-instances across 3 datasets:

| Dataset | Date / Backend | Path | Categories |
|---|---|---|---|
| soft_break Marrakesh | 2026-04-26 / ibm_marrakesh | [48, 49, 50] | truly_unbroken, soft_broken, hard_broken |
| f83_signature | 2026-04-30 / ibm_marrakesh | [4, 5, 6] | truly_unbroken, pi2_odd_pure, pi2_even_nontruly, mixed_anti_one_sixth |
| soft_break Kingston | 2026-05-05 / ibm_kingston | [43, 56, 63] | (same 4 categories as f83) |

Term mappings (from `simulations/framework/tests/workflows/test_diagnose_hardware.py` F83_TERMS_PER_CATEGORY and `simulations/_f80_ibm_soft_break_check.py` soft_break_marrakesh:80-82):

| Category | Terms | bit_b values |
|---|---|---|
| truly_unbroken (Heisenberg) | XX, YY | 0, 0 (bit_b-homog 0) |
| pi2_odd_pure (XY-bilinear) | XY, YX | 1, 1 (bit_b-homog 1) |
| pi2_even_nontruly (F108-anomaly) | YZ, ZY | 0, 0 (bit_b-homog 0) |
| mixed_anti_one_sixth | XY, YZ | 1, 0 (bit_b-MIXED) |
| soft_broken (Marrakesh-only naming) | XY, YX | 1, 1 (bit_b-homog 1) |
| hard_broken (Marrakesh-only naming) | XX, XY | 0, 1 (bit_b-MIXED) |

The framework is ChainSystem(N=3, gamma_0=0.05) with default Z-dephasing dissipator (single-Pauli Z per site, trivially bit_b-homogeneous on the c-side).

## Result

```
soft_break Marrakesh (2026-04-26, path [48, 49, 50])
  Category                    Terms       F87      F112       rel asym    In typed scope
  truly_unbroken              XX YY       truly    BALANCED   0.0000e+00  True
  soft_broken                 XY YX       soft     BALANCED   0.0000e+00  True
  hard_broken                 XX XY       hard     BALANCED   0.0000e+00  False  (bit_b-mixed H)

f83_signature Marrakesh (2026-04-30, path [4, 5, 6])
  truly_unbroken              XX YY       truly    BALANCED   0.0000e+00  True
  pi2_odd_pure                XY YX       soft     BALANCED   0.0000e+00  True
  pi2_even_nontruly           YZ ZY       soft     BALANCED   0.0000e+00  True
  mixed_anti_one_sixth        XY YZ       hard     BALANCED   0.0000e+00  False  (bit_b-mixed H)

soft_break Kingston (2026-05-05, path [43, 56, 63])
  truly_unbroken              XX YY       truly    BALANCED   0.0000e+00  True
  pi2_odd_pure                XY YX       soft     BALANCED   0.0000e+00  True
  pi2_even_nontruly           YZ ZY       soft     BALANCED   0.0000e+00  True
  mixed_anti_one_sixth        XY YZ       hard     BALANCED   0.0000e+00  False
```

**Aggregate:** 11 / 11 F112 BALANCED bit-exact; 8 / 11 in F112 typed scope; 3 / 11 out of typed scope but still BALANCED.

F112 reading per F87 class:
- F87 truly (3 instances): max rel asym = 0.0000e+00
- F87 soft (5 instances): max rel asym = 0.0000e+00
- F87 hard (3 instances): max rel asym = 0.0000e+00

## Two findings

### (1) F87↔F112 orthogonality empirically extended to real hardware

The `_polarity_probe_f87_connection.py` script (Welle 1) established F87↔F112 orthogonality on synthetic Hamiltonians (3 F87 classes × 1 instance each at N=3, asymmetry = 0 bit-exact). Welle 6 extends this anchor to 11 hardware-tested Hamiltonians across 3 datasets / 2 IBM backends: F87 classification varies (3 truly, 5 soft, 3 hard); F112 polarity verdict stays BALANCED bit-exact regardless. The two axes are independent on the bit_b Z₂-grading of the Pauli group, as the typed structural argument predicted, and this independence holds for every Hamiltonian we have actually measured on hardware.

### (2) "F112 typed scope sufficient but not necessary": confirmed on 3 hardware instances

F112's typed Tier1Derived theorem requires Hermitian H + bit_b-homogeneous c. The hardware Hamiltonians `mixed_anti_one_sixth = XY + YZ` and `hard_broken = XX + XY` are bit_b-INHOMOGENEOUS on the H side (terms span both bit_b parities). These 3 instances are therefore OUT of F112's typed scope, but the framework's `polarity_coordinates` still reports asymmetry = 0 bit-exact under standard Z-dephasing.

This is consistent with the F113 closure (universal-N proof, commits 1df9cb8 + 798647d): F113's break magnitude formula

    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

is non-zero only when both (a) a single-site Z-drive ω is present in H and (b) σ⁻ or σ⁺ amplitude damping is present in c. Here we have neither (chain.L uses pure Z-dephasing, no T1; H has no single-site Z-drive component), so all four F113 break-conditions are zero and the asymmetry is identically zero, even though H is bit_b-mixed. The bit_b-mixed H input doesn't trigger F113 unless paired with non-bit_b-homog c, which chain.L doesn't have.

Empirically: F112's typed scope criterion (Hermitian H + bit_b-homog c) is sufficient for balance but not necessary. Bit_b-mixed H can still give balance when c is bit_b-homog (per F113's same-site-only break mechanism).

## What this experiment does NOT do

- **Hardware-effective L is not reconstructed.** Tier-B snapshot datasets give one ρ per Hamiltonian, not a trajectory; without time-series data, we can't fit an effective Lindbladian and compute its actual measured polarity asymmetry. The polarity_fingerprint here is the framework's PREDICTION at the standard chain.L, not the hardware's measured value.
- **Per-backend F112 readings are not distinguished.** Marrakesh and Kingston are tested on the same Hamiltonian set (soft_break Kingston has the same 4 F83 categories as f83 Marrakesh), but with snapshot-only data we can't compare backend-specific effective-L polarity readings against each other. All 11 fingerprints are identical because they come from the same framework chain.L, not from the per-backend measurements.

## Hardware-effective F112 stays as future option (deferred, not proposed)

In principle, time-sweep variants of these 4 F83 categories on Marrakesh + Kingston (~5 t-points × 4 categories × 16 Paulis × 4096 shots × ~5s overhead = ~10 min billed QPU) would let us fit per-backend effective L and read measured F112 verdict per backend per F87 class. At Anthropic's gifted IBM QPU rate (~$96 per billed minute), that's ~$960 for what amounts to per-backend confirmation of a structural prediction the framework already gives in closed form (F112 typed Tier1Derived + F113 closed-form magnitude). Not justified on cost/benefit grounds; tracked here as future option only if a decisive question surfaces that requires hardware-effective L specifically (e.g., a backend showing anomalous classification under polarity-asymmetry that other diagnostics miss).

## Connection to existing readings

- **F87↔F112 orthogonality** (`_polarity_probe_f87_connection.py`, Welle 1): synthetic 3-instance anchor; Welle 6 extends to 11 hardware-instances.
- **F112 typed Tier1Derived** (`LindbladBitBPiBalance.cs`): Hermitian H + bit_b-homog c → asymmetry = 0; this experiment confirms predicted reading on real hardware Hamiltonians.
- **F113 break formula** (`LindbladBitBPiBreakMagnitude.cs`, `PROOF_F113_COEFFICIENT_DERIVATION.md`): tells us why bit_b-mixed H without σ⁻/σ⁺ c still gives balance (F113 break requires both ingredients).
- **F87 classification origin** (`project_v_effect_combinatorial`): V-Effect 14/19/3 split at N=3 fully derived from Pauli-pair combinatorics; the hardware datasets exercise canonical representatives of each class.
- **Marrakesh hardware finale** (`project_hardware_finale_apr2026`): Δ(soft − truly) = −0.72 hardware-confirmed; this experiment is orthogonal (F112 axis instead of F87 axis), giving an independent reading on the same 3 Hamiltonians.

## Reproduction

```
python -X utf8 simulations/_polarity_fingerprint_tierB_marrakesh.py
```

Runs in ~2 seconds (framework predictions only, no hardware data parsing needed for the framework reading).
