# IBM Kingston Cusp-Precision Run, April 26 2026

The dense point-by-point confirmation of the CΨ = 1/4 crossing (F25), curated copy of an
external full run. Single Bell⁺ pair, sampled densely right across the cusp.

## What is in this directory

| File | Description |
|------|-------------|
| `cusp_precision_ibm_kingston_20260426_115939.json` | The 19-delay trajectory: per-delay real CΨ, the F25 fit, and the pointwise residuals |

## Why this run is here

This is the experiment with extremely many data points exactly at 1/4. Where the April-16
[cusp-slowing run](../ibm_cusp_slowing_april2026/) used only six delays per pair to bracket the
crossing (but saved full density matrices, so it carries phase), this run sampled the crossing
**densely**: 19 delays at factors 0.05 to 3.0 of the predicted t_cross, with about eight of them
packed across the fold (factors 0.85 to 1.15, CΨ ≈ 0.22 to 0.26). It is the point-by-point F25
confirmation.

**It carries no phase.** The JSON stores only a real scalar `cpsi` per delay (the F25 magnitude),
not the 4×4 density matrices. So this trajectory lives on the **real axis** of the c-plane: it is
the dense 1D crossing (the Ω = 0, head-on case), not a 2D spiral. The spiral structure (a complex
CΨ leaving the real axis) is the sparse, phase-carrying April-16 run; this run is the dense,
real-only magnitude crossing. The two are complementary, not the same experiment.

## Experiment summary

- **Backend:** ibm_kingston (Heron r2)
- **Date:** 2026-04-26
- **Job ID:** `d7mu36lqrg3c738lnda0`
- **Pair:** qubits 14, 15 (T2_min ≈ 303 μs, γ ≈ 0.00165 1/μs); γ_calib 0.015, γ_fit 0.014977
- **Protocol:** Bell⁺ (H + CX), native delay, 9-Pauli 2-qubit tomography, 2048 shots per basis, 171 circuits total
- **Delays:** 19, at factors [0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.96, 0.99, 1.0, 1.01, 1.04, 1.08, 1.15, 1.3, 1.5, 1.8, 2.3, 3.0] of t_cross
- **Result:** F25 = f·(1+f²)/6, f = exp(−4γt), fits with γ the only free parameter; `fit_rms` = 0.0097; `F_CROSS` = 0.8612, `K_CROSS` = 0.037350 (γ·t at the crossing)

## Structure of the JSON

```
{
  "mode": "hardware", "backend": "ibm_kingston", "job_id": "d7mu36lqrg3c738lnda0",
  "pair": { "qubits": [14, 15], "gamma_per_us": ... },
  "gamma_calib": 0.015, "gamma_fit": 0.014977,
  "delay_factors": [ ... 19 ... ],
  "cpsi_data": [ { "factor": ..., "t_us": ..., "cpsi": ... }, ... 19 ... ],
  "pointwise_residual_calib": [ ... ], "pointwise_residual_fit": [ ... ],
  "fit_rms": 0.0097, "F_CROSS": 0.8612, "K_CROSS": 0.037350
}
```

## Provenance

Curated copy of the full run, which lives external to this repo in the IBM tomography directory
(`ibm_quantum_tomography/run_cusp_precision.py`, results folder). See
[CRITICAL_SLOWING_AT_THE_CUSP.md](../../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) (the April-26
precision update) for the full analysis. The JSON here is the authoritative record of the values
the repo's plots read.
