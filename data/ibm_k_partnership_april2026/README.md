# IBM Marrakesh K-Partnership Cross-Validation, April 2026

Raw data from the K-partnership cross-validation experiment on IBM Marrakesh
(Heron r2), testing the receiver-menu folding under K-partnership at N=5.

Companion to [`experiments/IBM_K_PARTNERSHIP_SKETCH.md`](../../experiments/IBM_K_PARTNERSHIP_SKETCH.md).

## What is in this directory

| File | Description |
|------|-------------|
| `k_partnership_marrakesh_aer_20260425_140311.json` | Aer pre-flight with Marrakesh noise model on path [48, 49, 50, 51, 58] |
| `k_partnership_marrakesh_20260425_140913.json` | Live Marrakesh QPU run, same 5 receivers, same path |

## Experiment summary

- **Backend:** ibm_marrakesh (Heron r2)
- **Date:** 2026-04-25, 14:09 UTC
- **Protocol:** Trotterised XX+YY-only Heisenberg evolution (no ZZ; the K-symmetric subset of Heisenberg) at t = 0.8, 3 Trotter steps, on a 5-qubit linear path; 2-qubit state tomography on (qubit 48, qubit 58) with 9 Pauli bases at 8192 shots each
- **Path:** [48, 49, 50, 51, 58], top-rated by `rank_paths` on the Marrakesh calibration
- **Wall-clock:** ~5 minutes QPU + queue (Kingston was down at run time, hence Marrakesh)
- **Calibration reference:** `ClaudeTasks/IBM_R2_calibration_ibm_marrakesh/ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv`

## Headline results

K-partner pair deviations on hardware (mirror-pair MI(0, 4)):

| Pair | Hardware Δ | Relativ | Aer Δ |
|------|------------|---------|-------|
| (bonding:1, bonding:5) | 0.0277 | **60 %** | 0.0001 (0.25 %) |
| (bonding:2, bonding:4) | 0.0108 | **14 %** | 0.0001 (0.02 %) |
| bonding:3 (self-partner anchor) | 0.1700 absolute | — | — |

Per-receiver MI(0, 4) on live Marrakesh:

| Receiver | MI(0, 4) | S(ρ_AB) | Purity |
|----------|----------|---------|--------|
| bonding:1 | 0.0460 | 1.594 | 0.390 |
| bonding:2 | 0.0768 | 1.670 | 0.349 |
| bonding:3 | 0.1700 | 1.665 | 0.356 |
| bonding:4 | 0.0660 | 1.536 | 0.410 |
| bonding:5 | 0.0737 | 1.385 | 0.475 |

## Interpretation

K-partnership **breaks on Marrakesh hardware** by 14 to 60 % across the two non-self K-partner pairs, far above the Aer prediction of 0.02 to 0.25 %. The break is the hardware Site- and Bond-Asymmetrie diagnostic that the [PROOF_K_PARTNERSHIP](../../docs/proofs/PROOF_K_PARTNERSHIP.md) theorem predicts: K-equivariance of the Lindblad generator requires the Site-Dephasing γ_ℓ-profile and the bond hopping to be K-symmetric (γ_0 = γ_4, γ_1 = γ_3, hopping symmetric). Where they are not, K-partner pairs deviate.

**Source of the asymmetry on the chosen path:**

| Asymmetry source | Site 0 / Bond (0,1) | Site 4 / Bond (3,4) | Asymmetry factor |
|------------------|---------------------|---------------------|------------------|
| T1 (μs) | 238.1 | 259.5 | 1.09× |
| T2 (μs) | 231.7 | 241.4 | 1.04× |
| Readout error | 5.85 % | 3.06 % | 1.91× |
| CZ error | 0.876 % | 0.149 % | **5.88×** |

The 6× CZ-error asymmetry between the end-bonds is the dominant K-symmetry breaker on this Marrakesh path. T1/T2 alone (4-9 %) would not produce a 60 % K-partner spread; the 6× bond-error asymmetry can. This is a direct hardware diagnostic that the K-partner-pair comparison surfaces but the conventional Site-T2 averages do not.

## Mode-specific spread

The (bonding:1, bonding:5) pair spreads more than (bonding:2, bonding:4) by ~4×. Plausible explanation: bonding:5 has an alternating sign pattern (k=N highest mode), which makes it more sensitive to bond-by-bond accumulated phase errors than bonding:2 (which has node-at-center structure and less accumulated phase). The K-partner spread is therefore not just a hardware property; it is a hardware property weighted by the Mode's spatial coherence pattern.

## Caveats

- bonding:1 and bonding:5 hardware MI values (0.046, 0.074) are slightly *higher* than the Aer prediction (0.040, 0.039). At very low MI absolute values, two-qubit tomography reconstruction noise can systematically overestimate the signal (positivity-violating ρ_2q gets clipped to physical, biasing eigenvalues). Treat MI < 0.1 as approximate; the K-partner ratio is more robust than absolute values in this regime.
- Single run, one path, one calibration snapshot. Repeating on a different Marrakesh path or at a different calibration time would test whether the 60 % vs 14 % pattern is reproducible or path-specific.
- Heisenberg ZZ excluded by `--xx-only`. A separate Heisenberg-with-ZZ run on the same path would calibrate the boundary breakdown documented in PROOF_K_PARTNERSHIP "XXZ boundary scaling".

## Reproducing the run

```bash
python run_receiver_engineering.py --hardware \
  --backend ibm_marrakesh \
  --no-alt-z-bits \
  --k-list "1,2,3,4,5" \
  --xx-only \
  --path "48,49,50,51,58" \
  --calibration <path-to-marrakesh-csv>
```

Pipeline location: `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\` (external; not in the R=CΨ² repo). Cost ~5 QPU minutes on Marrakesh.

## Integrity

- **Drift check:** none available because Run 1 (2026-04-24) was on Kingston with full Heisenberg; the XX-only K-partnership run on Marrakesh has no prior anchor on the same backend. Self-consistency: bonding:2 and bonding:4 are K-partners and gave 0.0768 vs 0.0660 (14 % spread) — consistent with each other to within hardware noise scale.
- **Path uniqueness:** [48, 49, 50, 51, 58] was the top-ranked 5-qubit linear path by the calibration scoring (score 185.3). Second-best [1, 2, 3, 4, 5] scored 159.8.
- **Queue conditions:** Kingston was reported down at run time; Marrakesh was online. Backend status was checked before run.
