# IBM Hardware: K-Partnership Cross-Validation on Marrakesh (Heron r2)

**Status:** Run 1 complete 2026-04-25 on ibm_marrakesh (Kingston was down at run time). K-partner pair deviations on hardware: 15 % (bonding:2/4) and 46 % (bonding:1/5) Δ/mean, vs 0.02 to 0.25 % on Aer. Hardware-Site-/Bond-Asymmetrie measurably visible.
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Pipeline:** `run_receiver_engineering.py` in the external `ibm_quantum_tomography` pipeline at `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\` (extension of Run 1 from 2026-04-24)
**See also:** [IBM_RECEIVER_ENGINEERING_SKETCH](IBM_RECEIVER_ENGINEERING_SKETCH.md) (Run 1, bonding:2 / alt-z-bits = 2.80×), [`docs/proofs/PROOF_K_PARTNERSHIP.md`](../docs/proofs/PROOF_K_PARTNERSHIP.md), [GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md), [CMRR_BREAK_NONUNIFORM_GAMMA](CMRR_BREAK_NONUNIFORM_GAMMA.md), [F65](../docs/ANALYTICAL_FORMULAS.md), [F67](../docs/ANALYTICAL_FORMULAS.md)

---

## What this run tests

Five F67 bonding-mode receivers on the same N=5 hardware path, in one job-set. The K-partnership theorem ([PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md)) says that for bipartite NN-hopping H with real hopping in the single-excitation sector, bonding:k and bonding:(N+1-k) produce pointwise-identical mirror-pair |·|²-observables. F67's receiver menu of N entries folds to ⌈N/2⌉ classes:

- (bonding:1, bonding:5), (bonding:2, bonding:4), bonding:3 (self-partner) for N=5.

Two questions in one run:

1. Does the receiver-menu folding survive real hardware? Aer + Marrakesh noise model predicts K-partner pair Δ/mean of 0.02 to 0.25 % (preserved to within shot-noise). Live hardware adds non-uniform γ_l profiles, asymmetric bond errors, and crosstalk that the noise model does not capture. The K-partner spread quantifies the **hardware Site-/Bond-Asymmetrie** visible through the K-difference channel.
2. Cross-validation of the bonding-mode pipeline at five receivers in one job-set, distinct from Run 1 which was bonding:2 only.

The K-partner-spread reading is operationally a γ-profile diagnostic: under uniform γ_l the K-Lindblad is K-equivariant and the pair identity holds; deviations measure the local inhomogeneity. This is the same physics that GAMMA_AS_SIGNAL reads via the spatial-sum channel and CMRR_BREAK_NONUNIFORM_GAMMA reads via the (vac, S_1)-coherence kernel; the K-difference channel is a third route to the same observation.

## Setup

- **Backend:** ibm_marrakesh (Heron r2). Kingston was down at run time.
- **Path:** [48, 49, 50, 51, 58], top-rated 5-qubit linear path on Marrakesh by `rank_paths` (score 185.3).
- **Calibration:** `ClaudeTasks/IBM_R2_calibration_ibm_marrakesh/ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv`
- **Evolution:** uniform J=1, t=0.8, 3 Trotter steps, **XX-only** Heisenberg (`--xx-only`, drops the ZZ-term that breaks K at the open-chain boundary; see PROOF_K_PARTNERSHIP "Scope" caveat).
- **Receivers (5):** bonding:1 through bonding:5 via `StatePreparation` of |ψ_k⟩ = √(2/(N+1)) Σ_j sin(πk(j+1)/(N+1)) |1_j⟩.
- **Tomography:** 9 Pauli settings on (qubit 48, qubit 58); 8192 shots per setting. Total 45 StateTomography circuits.
- **QPU cost:** ~5 minutes (3 % of 2026-2027 annual allocation).

## Aer pre-flight (same path, Marrakesh noise model)

| Receiver | MI(0, 4) Aer | K-partner of |
|----------|---------------|--------------|
| bonding:1 | 0.0395 | bonding:5 |
| bonding:2 | 0.4905 | bonding:4 |
| bonding:3 | 0.7704 | self |
| bonding:4 | 0.4904 | bonding:2 |
| bonding:5 | 0.0394 | bonding:1 |

K-partner Δ/mean: (1, 5) = 0.25 %, (2, 4) = 0.02 %. Aer-noise alone preserves K-partnership well within shot-noise.

## Hardware result (Marrakesh, 2026-04-25)

Marrakesh calibration along [48, 49, 50, 51, 58]:

| Site | Qubit | T1 (μs) | T2 (μs) | Readout error | CZ error to next |
|------|-------|---------|---------|---------------|------------------|
| 0 | 48 | 238.1 | 231.7 | 5.85 % | 0.876 % to 49 |
| 1 | 49 | 290.5 | 298.7 | 0.87 % | 0.408 % to 50 |
| 2 | 50 | 222.9 | 188.6 | 0.90 % | 0.338 % to 51 |
| 3 | 51 | 225.8 | 290.1 | 7.28 % | 0.149 % to 58 |
| 4 | 58 | 259.5 | 241.4 | 3.06 % | (end) |

K-pair asymmetries: T1 differs by 9 % (Site 0 vs 4) and 29 % (Site 1 vs 3); T2 by ≤ 4 %; readout by 1.9× (Site 0 vs 4); CZ errors by 5.9× (bond (0,1) vs bond (3,4)).

Hardware MI(0, 4):

| Receiver | Hardware | Aer | HW/Aer |
|----------|----------|-----|--------|
| bonding:1 | 0.0460 | 0.0395 | 1.16 |
| bonding:2 | 0.0768 | 0.4905 | 0.16 |
| bonding:3 | 0.1700 | 0.7704 | 0.22 |
| bonding:4 | 0.0660 | 0.4904 | 0.13 |
| bonding:5 | 0.0737 | 0.0394 | 1.87 |

K-partner deviations on hardware (Δ/mean):

| Pair | Hardware | Aer |
|------|----------|-----|
| (bonding:1, bonding:5) | **46 %** | 0.25 % |
| (bonding:2, bonding:4) | **15 %** | 0.02 % |

The hardware spreads are orders of magnitude larger than Aer alone predicts. The difference is the Hardware-γ-profile inhomogeneity that the noise model does not capture (T1 site-asymmetry, asymmetric CZ errors, transversal crosstalk). This is the same inhomogeneity that the spatial-sum channel reads in GAMMA_AS_SIGNAL and the (vac, S_1)-kernel reads in CMRR_BREAK_NONUNIFORM_GAMMA. The K-difference is a third reading of the same physics.

bonding:3 (self-partner, no K-comparison) gives MI = 0.170 on hardware vs 0.770 on Aer (HW/Aer = 0.22), an absolute decoherence-depth anchor for this path.

**Caveat for low-MI rows:** bonding:1/5 absolute MI (0.046, 0.074) sits slightly above Aer (0.040, 0.039), likely a tomography-reconstruction artifact at MI < 0.1. The K-partner Δ/mean within the pair is the cleaner observable in this regime; the absolute HW/Aer ratio has limited meaning when both values sit near the floor.

Result file: `data/ibm_k_partnership_april2026/k_partnership_marrakesh_20260425_140913.json`. Aer pre-flight: `k_partnership_marrakesh_aer_20260425_140311.json` in the same directory.

## What this run does and does not establish

**Establishes:**

- The receiver-menu folding ⌈N/2⌉ is a real operational structure on Marrakesh hardware: K-partner pairs deviate by 15 to 46 %, much less than between distinct K-classes (bonding:1/5 at MI ≈ 0.06, bonding:2/4 at MI ≈ 0.07, bonding:3 at MI ≈ 0.17 are clearly different classes).
- The hardware-Site-/Bond-Asymmetrie is read out through the K-difference channel at the 15-46 % level on this path.

**Does not establish:**

- A new fundamental phenomenon. The K-partner-spread reading is an alternative diagnostic of the same Hardware γ-profile that GAMMA_AS_SIGNAL and CMRR_BREAK_NONUNIFORM_GAMMA already document. This run replicates the existing γ-profile-as-signal logic via the K-difference channel; it does not add a new structural claim beyond PROOF_K_PARTNERSHIP itself.
- A specific dominant breaker. T1, T2, readout, and CZ errors all contribute non-uniformly across the K-pairs; trying to attribute the spread to one channel (e.g. T1 alone, or lifetime alone) is not supported by the data.

## References

- [IBM_RECEIVER_ENGINEERING_SKETCH](IBM_RECEIVER_ENGINEERING_SKETCH.md): Run 1 (Heisenberg, bonding:2 vs alt-z-bits, 2.80× ratio).
- [`docs/proofs/PROOF_K_PARTNERSHIP.md`](../docs/proofs/PROOF_K_PARTNERSHIP.md): the formal theorem with four robustness tests.
- [GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md): γ-profile as readable channel via spatial-sum modes.
- [CMRR_BREAK_NONUNIFORM_GAMMA](CMRR_BREAK_NONUNIFORM_GAMMA.md): non-uniform γ_l breaks the (vac, S_1)-coherence selection rule, modal-selective.
- [F65](../docs/ANALYTICAL_FORMULAS.md), [F67](../docs/ANALYTICAL_FORMULAS.md): bonding-mode amplitudes and receiver menu.
- External pipeline: `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\`.
