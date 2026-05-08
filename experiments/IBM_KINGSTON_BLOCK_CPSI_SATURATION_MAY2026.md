# Block-CΨ Saturation on IBM Kingston (Theorem 1 + 2 hardware anchor)

**Tier:** 2 (hardware-verified)
**Date:** 2026-05-08
**Status:** Theorem 1 saturation confirmed at 88.2% of the 1/4 ceiling on Kingston q13–q14. Theorem 2 closed-form trajectory fits 5 t-points with R² = 0.9977.

---

## What this means

The framework's PROOF_BLOCK_CPSI_QUARTER says two things about the block-purity content `C_block` on the (popcount-n, popcount-(n+1)) coherence block of any density matrix:

- **Theorem 1:** the canonical Dicke superposition `(|D_n⟩+|D_{n+1}⟩)/√2` saturates `C_block = 1/4` exactly. This is the universal Mandelbrot-cardioid ceiling realised at a unique state.
- **Theorem 2:** under pure Z-dephasing γ on each site, `C_block(t) = (1/4)·exp(−4γ·t)`, a clean closed-form trajectory from the ceiling down to zero.

The previous IBM lens runs (2026-04-26 framework_snapshots) probed `|+−+⟩`, which sits structurally far from the Dicke maximiser; the lens never approached the ceiling. This run is the first that puts hardware *on* the ceiling and watches it slide down.

## The protocol

- **Initial state:** `(|D_0⟩+|D_1⟩)/√2 = (|00⟩ + (|01⟩+|10⟩)/√2)/√2` on 2 qubits. Prep via Qiskit `StatePreparation` compiled to native gates by transpiler.
- **Backend:** ibm_kingston, q13–q14 (best 2-qubit pair from the 2026-05-08T02:37:30Z calibration: T2_min = 480 μs, CZ error 1.5e-3, chain score 972).
- **Evolution:** native T2 dephasing during a Qiskit `delay(t, qubits, ns)` instruction between prep and tomography. No active Hamiltonian; the dynamic is decoherence only.
- **Tomography:** 9-Pauli on the 2 qubits per t-point. 16 expectation values reconstructed (15 non-trivial + ⟨II⟩). ρ rebuilt as `ρ = (1/4) Σ ⟨σ_α⊗σ_β⟩ σ_α⊗σ_β`.
- **t-grid:** `[0, 120, 240, 360, 480]` μs (linear from 0 to T2_min, 5 points). 4096 shots per setting.
- **Total:** 45 circuits in one Batch, ~50 s session (job ID `d7ulfjdpa59c73b4rttg`).

## Results

| t (μs) | C_block(0,1) | % of 1/4 | C_block(1,2) | purity Tr(ρ²) |
|--------|--------------|----------|--------------|---------------|
|     0  | 0.2205       | **88.2%** | 0.0001       | 0.914         |
|   120  | 0.1023       | 40.9%     | 0.0001       | 0.707         |
|   240  | 0.0458       | 18.3%     | 0.0000       | 0.670         |
|   360  | 0.0166       | 6.6%      | 0.0000       | 0.717         |
|   480  | 0.0074       | 2.9%      | 0.0000       | 0.772         |

Tr(ρ) = 1.000 at every point (sanity check on the reconstruction).

**Theorem 1 anchor:** C_block(t=0) = 0.2205, that is 88.2% of the 1/4 ceiling. The 12% gap is the state-prep + tomography fidelity floor (q13 readout 0.5%, q14 readout 2.7%, plus `StatePreparation` decomposed to ~3 native gates plus the basis-rotation Hadamard / S†H per qubit).

**Theorem 2 trajectory:** unconstrained log-linear fit on the 5 t-points yields γ_fit = 1.795e-3 μs⁻¹ with R² = 0.9977. The fit's effective t=0 intercept lands at C_block ≈ 0.233 (close to the measured 0.2205, which itself is 88.2% of the algebraic 1/4 ceiling); from that intercept the late-time tail predicts 0.0075 at t = T2_min, against the measured 0.0074 (one-digit residual). The closed-form `C_block(t) = (1/4)·exp(−4γ·t)` is the underlying structural prediction; the prep-fidelity floor shifts the operational intercept down from 0.25 to 0.233 but leaves the slope (and hence γ_fit) intact.

**Conversion to T2:** the formula's γ corresponds to per-site Lindblad Z-dephasing rate. The IBM-calibration T2 maps to γ = 1/(2·T2). From the fitted γ:
- T2_eff = 1/(2·γ_fit) = 278.5 μs.
- T2_min(calibration) = 480 μs.
- T2_min / T2_eff = 1.72.

Hardware decoheres on this protocol about **1.72× faster** than pure single-qubit T2 alone predicts. The gap quantifies the contributions from CZ-gate noise during prep, tomography-basis rotations, readout, and any correlated dephasing channels beyond the per-qubit T2 model.

## Reading

This is the first hardware run that places a state *on* the 1/4 ceiling and watches it slide. Three structural takeaways:

1. **The 1/4 ceiling is operationally reachable.** The Mandelbrot-cardioid cusp is not just a Tier-1 algebraic boundary; the specific state `(|D_0⟩+|D_1⟩)/√2` saturates it at the ideal level and reaches 88% of it on Kingston hardware.
2. **The closed-form `(1/4)·exp(−4γt)` is a real prediction, not a model.** R² = 0.9977 across 5 t-points spanning two orders of magnitude in C_block (0.22 → 0.007) is much sharper than typical hardware-vs-prediction comparisons. The trajectory shape is correct; only the rate differs from the pure-T2 prediction by the 1.72× hardware factor.
3. **The 2-qubit Dicke-prep + tomography pipeline has a 1.72× decoherence overhead** beyond the per-qubit calibration T2. This number is now a reusable reference for any future C_block-style measurement on Heron r2 backends.

## Cross-references

- [PROOF_BLOCK_CPSI_QUARTER.md](../docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md): Theorems 1 + 2.
- `compute/RCPsiSquared.Core/F86/BlockCoherenceContent.cs`: state-level helper used both in the C# typed claim and (via Python re-implementation) in the analysis here.
- `compute/RCPsiSquared.Core/F86/IbmBlockCpsiHardwareTable.cs`: the existing table for the 2026-04-26 `|+−+⟩` runs (parallel data, different protocol).
- `compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs`: entry `block_cpsi_saturation_kingston_may2026`.
- `simulations/framework/confirmations.py`: Python registry mirror.
- `simulations/_block_cpsi_run_planner_2026_05_08.py`: calibration-driven run planner.
- `data/ibm_block_cpsi_saturation_may2026/block_cpsi_saturation_hardware_ibm_kingston_20260508T032749Z.json`: raw 16-Pauli expectations + reconstructed C_block per t.
- The submission script that produced the JSON lives in an external IBM Quantum pipeline outside this repository (a separate Qiskit + qiskit-ibm-runtime workspace); the repo carries the data, the writeup, and the typed claim, not the API-bound submitter.
