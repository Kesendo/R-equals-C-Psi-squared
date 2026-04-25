# IBM Marrakesh K-Partnership Cross-Validation, April 2026

Raw data from a five-receiver K-partnership cross-validation on IBM Marrakesh
(Heron r2). Companion to [`experiments/IBM_K_PARTNERSHIP_SKETCH.md`](../../experiments/IBM_K_PARTNERSHIP_SKETCH.md).

## Files

| File | Description |
|------|-------------|
| `k_partnership_marrakesh_aer_20260425_140311.json` | Aer pre-flight with Marrakesh noise model on path [48, 49, 50, 51, 58] |
| `k_partnership_marrakesh_20260425_140913.json` | Live Marrakesh QPU run, same 5 receivers, same path |

## Experiment summary

- **Backend:** ibm_marrakesh (Heron r2). Kingston was down at run time.
- **Date:** 2026-04-25, 14:09 UTC
- **Protocol:** uniform J=1, t=0.8, 3 Trotter steps, **XX+YY hopping only** (no ZZ; the K-symmetric subset of XXZ for the single-excitation sector). Five F67 bonding-mode receivers (k = 1..5). 9 Pauli bases × 8192 shots, 2-qubit state tomography on (qubit 48, qubit 58).
- **Path:** [48, 49, 50, 51, 58], top-rated 5-qubit linear path on Marrakesh by `rank_paths` (score 185.3).
- **Wall-clock:** ~5 minutes QPU + queue.
- **Calibration:** `ClaudeTasks/IBM_R2_calibration_ibm_marrakesh/ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv`

## Results

K-partner pair deviations (Δ / mean(MI_k, MI_{N+1-k})):

| Pair | Hardware Δ/mean | Aer Δ/mean |
|------|-----------------|------------|
| (bonding:1, bonding:5) | 46 % | 0.25 % |
| (bonding:2, bonding:4) | 15 % | 0.02 % |
| bonding:3 (self) | — | — |

Per-receiver MI(0, 4) on live Marrakesh:

| Receiver | MI(0, 4) | S(ρ_AB) | Purity |
|----------|----------|---------|--------|
| bonding:1 | 0.0460 | 1.594 | 0.390 |
| bonding:2 | 0.0768 | 1.670 | 0.349 |
| bonding:3 | 0.1700 | 1.665 | 0.356 |
| bonding:4 | 0.0660 | 1.536 | 0.410 |
| bonding:5 | 0.0737 | 1.385 | 0.475 |

## Reading

K-partnership ([PROOF_K_PARTNERSHIP](../../docs/proofs/PROOF_K_PARTNERSHIP.md)) predicts pair-identity to within shot-noise under uniform γ_l-profile and K-symmetric Hamiltonian. Aer with the Marrakesh noise model gives 0.02 to 0.25 %. Live Marrakesh shows 15 to 46 %; the difference is the Hardware-γ-profile inhomogeneity along the chosen path that the noise model does not capture.

This is the same physics that [GAMMA_AS_SIGNAL](../../experiments/GAMMA_AS_SIGNAL.md) reads via the spatial-sum channel and [CMRR_BREAK_NONUNIFORM_GAMMA](../../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md) reads via the (vac, S_1)-coherence kernel. The K-difference channel is a third route to the same observation: non-uniform γ_l breaks a closure that the symmetric setup preserves; the size of the break measures the inhomogeneity.

The run does not establish a new fundamental phenomenon beyond what is already in PROOF_K_PARTNERSHIP. It is a hardware existence proof: the receiver-menu folding ⌈N/2⌉ is operationally distinguishable on Marrakesh (different K-classes give different MI; K-partner pairs sit much closer to each other than to other classes), and the size of the K-partner deviation is consistent with the Marrakesh path-asymmetry.

## Caveats

- bonding:1 and bonding:5 hardware MI values (0.046, 0.074) sit slightly above the Aer prediction (0.040, 0.039). At MI < 0.1, two-qubit tomography reconstruction noise can systematically bias eigenvalues (positivity-projected ρ). The K-partner Δ/mean within the pair is the cleaner observable than absolute HW/Aer ratios in this regime.
- Single run, one path, one calibration snapshot. The 15 % vs 46 % pair-asymmetry is not a falsifiable structural prediction; it reflects this specific Marrakesh path's Site-/Bond-Asymmetrie. A different path (or different calibration day) gives different numbers.
- ZZ excluded by `--xx-only`. With full Heisenberg (ZZ included), Aer pre-flight already shows 7-16× K-partner spread on the open chain, dominated by the F65-boundary-breakdown of K (PROOF_K_PARTNERSHIP "XXZ boundary scaling"); that is a separate effect from the hardware-γ-profile reading.

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

Pipeline: `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\` (external; not in the R=CΨ² repo). Cost ~5 QPU minutes.
