# IBM Marrakesh F83 4-Hamiltonian Signature Test, April 2026

Live hardware verification of F83's structural prediction that the four
F77 Π²-classes (truly / pure Π²-odd / pure Π²-even-non-truly / mixed)
produce **operationally distinguishable Pauli-expectation patterns** on
N=3 chains, and a first hardware look at whether any deviation from
Trotter+γ_Z=0.1 carries the F82/F84 amplitude-damping signature.

Run on `ibm_marrakesh` (Heron r2) on path **[4, 5, 6]** (top-ranked
by 2026-04-30T16:25Z calibration, score 0.0162, rank 1 of 223
3-qubit chains) using the same pipeline as the April 26 soft_break
run with `CATEGORIES` extended to four entries.

## Files

| File | Description |
|------|-------------|
| `f83_signature_ibm_marrakesh_20260430_190035.json` | Raw counts and reconstructed 2-qubit Pauli expectations for **4 Hamiltonians × 9 tomography bases** on (q0=4, q2=6); 36 circuits, 4096 shots/circuit. |

## Experiment summary

- **Backend:** ibm_marrakesh (Heron r2). Same backend as April 26 soft_break and other April 2026 runs.
- **Date:** 2026-04-30, ~19:00 UTC.
- **Job ID:** `d7pol1e7g7gs73cf7j90`.
- **Path:** [4, 5, 6]. Calibration: T1 = 206 / 144 / 351 μs, T2 = 184 / 121 / 151 μs, RO err = 0.42 / 0.38 / 0.37 %, CZ err = 0.633 % (4-5) / 0.287 % (5-6).
- **Initial state:** \|+−+⟩ X-Néel (Hadamard on q0, q2; Hadamard+Z on q1).
- **Four Hamiltonians:**
  - `truly_unbroken` = XX+YY (truly; M = 0)
  - `pi2_odd_pure` = XY+YX (pure Π²-odd; F83 anti-fraction = 1/2)
  - `pi2_even_nontruly` = YZ+ZY (pure Π²-even non-truly; F83 anti-fraction = 0)
  - `mixed_anti_one_sixth` = XY+YZ (equal-Frobenius mixed; F83 anti-fraction = 1/6)
- **Evolution:** uniform J = 1, t = 0.8, n_trotter = 3 (first-order Trotter via Qiskit `PauliEvolutionGate`).
- **Tomography:** 9 Pauli bases on (q0, q2), 4096 shots/basis. Total 36 circuits.
- **QPU cost:** ~3 minutes wall-clock.

## Results: hardware vs Trotter+γ_Z=0.1 prediction

The Trotter+γ_Z=0.1 model fit the April 26 soft_break ⟨X₀Z₂⟩ to within
0.0014 on path [48, 49, 50]. Here, on path [4, 5, 6], the same model is
the prediction baseline. Aer with [4, 5, 6]-calibrated T1/T2/CZ noise
matched the Trotter prediction to within 0.001 for the anchor observable.

| Pauli (P_0, P_2) | Trotter pred | Aer (calibrated) | **Hardware [4,5,6]** | Δ(HW − Trotter) |
|------------------|--------------|------------------|-----------------------|------------------|
| ⟨X₀Z₂⟩ truly | 0.000 | -0.013 | **-0.001** | 0.001 |
| ⟨X₀Z₂⟩ pi2_odd | -0.723 | -0.722 | **-0.849** | **-0.126** ⚠ |
| ⟨X₀Z₂⟩ pi2_even_nt | 0.000 | +0.023 | **+0.030** | 0.030 |
| ⟨X₀Z₂⟩ mixed | -0.065 | +0.153 | **+0.154** | 0.219 |
| ⟨Z₀X₂⟩ pi2_odd | -0.506 | (not measured all) | **-0.566** | -0.061 |
| ⟨Z₀X₂⟩ mixed | -0.621 | — | **-0.721** | -0.100 |
| ⟨X₀X₂⟩ pi2_even_nt | +0.726 | +0.825 | **+0.919** | **+0.193** ⚠ |
| ⟨Y₀Z₂⟩ truly | +0.381 | — | **+0.670** | **+0.289** ⚠ |
| ⟨Y₀Y₂⟩ truly | +0.323 | — | **+0.535** | **+0.212** |
| ⟨Z₀Z₂⟩ truly | +0.489 | +0.225 | +0.215 | -0.275 |
| ⟨Z₀Z₂⟩ pi2_odd | +0.199 | +0.222 | +0.210 | +0.011 |
| ⟨Z₀Z₂⟩ pi2_even_nt | 0.000 | -0.015 | -0.004 | -0.004 |
| ⟨Z₀Z₂⟩ mixed | -0.395 | -0.252 | -0.229 | +0.166 |

## What this run confirms

**The four F77 Π²-classes are operationally distinguishable on hardware
via Pauli-expectation pattern.** Each category has at least one
unique-fingerprint observable that separates it from all three other
categories at >>10σ (statistical error 0.0156 at 4096 shots):

- `truly_unbroken`: ⟨Y₀Z₂⟩ = +0.67 (others all near 0); ⟨X₀Z₂⟩ ≈ 0
- `pi2_odd_pure`: ⟨X₀Z₂⟩ = −0.85 (uniquely large negative)
- `pi2_even_nontruly`: ⟨X₀X₂⟩ = +0.92 (uniquely large positive)
- `mixed_anti_one_sixth`: ⟨Z₀X₂⟩ = −0.72, ⟨X₀Z₂⟩ = +0.15 (sign-flip on Z,X axis vs pi2_odd)

This validates F83's framework-level claim that the Π²-class structure
of the Hamiltonian (which is closed-form computable from the bilinear
letters alone) maps onto distinct hardware-observable patterns.

## What this run leaves open

**The pi2_odd_pure ⟨X₀Z₂⟩ magnitude is 0.13 stronger than Trotter
predicts** (−0.849 vs −0.723) and 0.13 stronger than Aer with
calibration-derived noise (−0.722). The amplification is **not in the
Aer noise model**. Two candidate explanations:

1. **F82/F84 vacuum amplitude-damping signature.** F82 says
   ‖D_T1_odd‖_F = γ_T1 · √N · 2^(N-1), i.e. T1 amplitude damping leaks
   into M_anti (the dynamics-generator side of M). Q5 has the shortest
   T1 in the path (144 μs), and the F84 reading is that the
   σ⁻-asymmetric (vacuum) component of amplitude damping is what reaches
   M_anti. This would predict a stronger soft signal than pure Z-dephasing.
   But Aer's `thermal_relaxation_error` channel also models amplitude
   damping and gives −0.72, not −0.85. So either (a) Aer's channel
   composition with depolarization is missing structure, or (b) the
   F82/F84 interpretation needs an additional element not in the simple
   noise sum.

2. **Coherent calibration drift.** Calibration CSV was downloaded at
   16:25 UTC; job ran ~2.5 h later. Heron r2 calibration is refreshed
   ~hourly, so drift is plausible. A coherent ZZ-coupling shift or
   single-qubit phase miscalibration would amplify (or attenuate) the
   soft signal without breaking the Π²-class symmetry pattern.

The April 26 run on [48, 49, 50] gave ⟨X₀Z₂⟩ = −0.711 for the same
Hamiltonian XY+YX. The 2026-04-30 run on [4, 5, 6] gives −0.849. The
difference 0.138 is path-dependent. We can't resolve (1) vs (2) from
this run alone; a controlled re-run on [48, 49, 50] with the same
calibration window would isolate the path effect.

**The Y-basis observables drift systematically high.** ⟨Y₀Z₂⟩ truly:
predicted 0.38, measured 0.67. ⟨Y₀Y₂⟩ truly: 0.32 vs 0.54. The
Y-basis preparation is `sdg + H` followed by computational measurement;
this insertion may carry coherent gate-error not modeled by depolarization.

## What this run does NOT establish

- That F83's anti-fraction r = ‖H_even_nontruly‖² / ‖H_odd‖² is invertible
  from data. The prediction is structural (closed-form on H), not fitted
  from observables. We measured signatures, not r.
- That F82/F84 explains the X₀Z₂ amplification. This is a candidate, not
  a confirmation. Needs a controlled T1-versus-path experiment.
- That the discrimination scales to N > 3. Open.
- That arbitrary mixed Hamiltonians (other than 1/2, 0, 1/6) follow the
  same pattern. Three F83-discrete cases tested; the continuous tunability
  of r remains a structural prediction, not a hardware-mapped one.

## Reproducing the run

```bash
cd "D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography"
python run_soft_break.py --hardware --backend ibm_marrakesh --path 4,5,6 --shots 4096
```

The run-time on Marrakesh queue + execution: ~3 minutes total.

## Reading

- [PROOF_F83_PI_DECOMPOSITION_RATIO](../../docs/proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md): the closed-form anti-fraction.
- [PROOF_F82_T1_DISSIPATOR_CORRECTION](../../docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md): the T1-into-M_anti closed form, candidate for the X₀Z₂ amplification.
- [PROOF_F84_AMPLITUDE_DAMPING](../../docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md): only σ⁻/σ⁺ break the Π palindrome among single-qubit dissipators.
- [ON_THE_RESIDUAL](../../reflections/ON_THE_RESIDUAL.md): consolidating reflection on F80–F85.
- [ibm_soft_break_april2026/](../ibm_soft_break_april2026/): the April 26 anchor run on path [48, 49, 50] that this test extends.
- [`_f83_signature_predictions.py`](../../simulations/_f83_signature_predictions.py): closed-form Trotter+γ_Z=0.1 predictions per category.
- [`_f83_aer_preflight.py`](../../simulations/_f83_aer_preflight.py): Aer noise simulation with [4, 5, 6] calibration values; matched Trotter to 0.001 on the anchor observable but missed the +0.13 hardware drift.
- [`_2qubit_dissipator_exploration.py`](../../simulations/_2qubit_dissipator_exploration.py) (commit ba8e861): partial F86, single-bond closed form for σ-channels; overlap structure open.
